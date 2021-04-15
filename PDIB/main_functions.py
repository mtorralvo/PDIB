"""Main functions used to create the models. Includes functions to superimpose chains, check clashes and the main recursive algorithm."""

from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import NeighborSearch
from Bio.PDB import Structure as pdb_struct
from Bio.PDB import Model as pdb_model
from Bio.PDB import PDBIO
import random
from process_input_files import id_generator
import logging

ppb = PPBuilder()
parser = PDBParser(PERMISSIVE = True, QUIET = True)

def superimpose(chainA, chainB, moving = None):

    """
    Superimposes two homolog chains applying the rotation matrix to the newly added chain, returning the rotated copy. It also returns the RMSD of the superimposition.
    chainA: chain of the previous structure (Chain object from PDB.Chain).
    chainB: homolog chain (Chain object from PDB.Chain) of the pairwise interaction to be added.
    moving: non-homolog chain (Chain object from PDB.Chain) that will be transformed using the rotation matrix.
    """

    residuesA = list(chainA.get_residues())
    residuesB = list(chainA.get_residues())

    # Only keep alpha carbons (or P in the case of nucleic acids).
    filteredA = [res for res in residuesA if ('CA' in res or 'P' in res)]
    filteredB = [res for res in residuesB if ('CA' in res or 'P' in res)]

    # Both chains must have the same length in order to be superimposed.
    if len(filteredA) > len(filteredB):
        filteredA = filteredA[:len(filteredB)]
    elif len(filteredA) < len(filteredB):
        filteredB = filteredB[:len(filteredA)]

    # Get alpha carbons in chainA
    atomsA = []
    for res in filteredA:
        for atom in res.get_atoms():
            if atom.id == 'CA' or atom.id == 'P':
                atomsA.append(atom)
                break

    # Get alpha carbons in chainB
    atomsB = []
    for res in filteredB:
        for atom in res.get_atoms():
            if atom.id == 'CA' or atom.id == 'P':
                atomsB.append(atom)
                break

    superimposer = Superimposer()           # Superimposer object

    superimposer.set_atoms(atomsA, atomsB)  # Superimposes both chains

    rmsd = superimposer.rms                 # RMSD

    # If a moving chain is specified, apply the rotation matrix to it.
    if moving:
        superimposer.apply(list(moving.get_atoms()))
        return rmsd, moving

    return rmsd  # Else, return the RMSD  


def check_clash(model, chain_to_add, clash_distance=2.5):

    """
    Checks wether a newly added chain has clashes with the rest of the structure. Returns a set containing the clashing chains and a boolean indicating whether the newly added chain is clashing with itself or not.
    model: previous structure.
    chain_to_add: new chain (Chain object from PDB.Chain).
    clash_distance: indicates a distance threshold (in Angstroms) to consider two atoms as clashing.
    """

    # Initialization of NeighbourSearch, that allows to find clashes
    neighbor_object = NeighborSearch(list(model.get_atoms()))

    structure_clashing_chains = set()
    total_clashes = 0

    # For each atom, clashes are compute
    for atom in chain_to_add.get_atoms():
        clashes = neighbor_object.search(atom.get_coord(), clash_distance)
        
        if len(clashes) > 0:    # If clashes are found...

            # Increase the total number of clashes and add the conflicting chains to the set
            for clash in clashes:
                structure_clashing_chains.add(clash.get_parent().get_parent().id)
                total_clashes += 1
    
    # In case the new chain is conflicting with several chains
    if len(structure_clashing_chains) > 1 and total_clashes > 20:
        return structure_clashing_chains, False
    

    elif len(structure_clashing_chains) == 1 and total_clashes > 20 and chain_to_add.id[1] == list(structure_clashing_chains)[0][1]:

        clash_chain = model[0][list(structure_clashing_chains)[0]]  # Conflictive chain
        RMSD = superimpose(clash_chain, chain_to_add)               # Compute RMSD

        # If the RMSD is lower than 3.0, the chain is clashing with itself.
        if RMSD <= 3.0:
            return structure_clashing_chains, True

        # Else the chain is clashing with another chain of the structure
        else:
            return structure_clashing_chains, False

    elif total_clashes > 20:
        return structure_clashing_chains, False

    # No clashes found
    else:
        return None, False

def check_structure_exists(structure, all_structures):

    """Returns a boolean indicating wether the input structure is in the set of all structures. Returns True if all the chains of the structure matches with all the chains of any structure present in all_structures."""

    structure = structure.copy()    # Copy the structure in order to preserve the original
    all_structures = all_structures.copy()

    # Generate a list containing all the chain IDs of the input structure
    chain_ids_structure = tuple(sorted([x.id[1] for x in structure.get_chains()]))

    for each_structure in all_structures:

        # Chain IDs of another structure
        chain_ids_each_structure = tuple(sorted([x.id[1] for x in each_structure.get_chains()]))

        # Checks wether both structures have the same chain IDs
        if chain_ids_structure == chain_ids_each_structure:

            # Selects a chain of the input structure
            chain_str = list(structure.get_chains())[0]
            id_str = chain_str.id[1]

            # Searches for chains with the same ID
            for chain_each_str in each_structure.get_chains():

                id_each_str = chain_each_str.id[1]

                # If a chain with the same name is found, superimpose both chains
                if id_str == id_each_str:

                    RMSD, structure = superimpose(chain_str, chain_each_str, moving=structure)

                    # If the superimposition yields a low RMSD (<3A) they are the same chain
                    if  RMSD > 3.0:
                        continue

                    partners = set()

                    # Creates a list of partners
                    for searching_partner in each_structure.get_chains():
                        partner_found = False

                        for possible_partner in structure.get_chains():

                            if partner_found is True:
                                break

                            if possible_partner.id[1] == searching_partner.id[1] and possible_partner not in partners:

                                    RMSD = superimpose(searching_partner, possible_partner)

                                    if RMSD <= 3.0:
                                        partners.add(possible_partner)
                                        partner_found = True

                        # If all chains have a partner returns True
                        if len(partners) == len(list(each_structure.get_chains())):
                            return True  

    # If not all chains match returns False
    return False    

def create_model(current_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = 0, complex = False, in_a_branch = False, non_brancheable_clashes=set(), tried_branch_structures=list(), stoich=None, verbose=False):

    """Recursively adds a new chain to the current structure. If looks for clashes and repeated chains in order to obtain only feasible models.
    current_structure: chain obtained in the previous recursion level.
    stored_structures: structures generated by the program.
    info_files: contains the paths of the pairwise interactions and the chains contained in each one.
    num_chains: number of chains of the final model. Specified by the user.
    num_models: maximum number of models that will be generated.
    all_models: boolean. If True, all the possible structures will be generated."""

    # If not all the models are required and the number of models required has been achieved, stop the recursion
    if all_models is False and num_models == len(stored_structures):
        return stored_structures

    if complex:
        recursion_level += 1

    # List of pairwise interaction files. 
    # A random behaviour is implemented as not all the conformations will be explored if all_models=False
    all_files= list(info_files.keys())
    random.shuffle(all_files)

    current_chains = list(current_structure.get_chains())

    something_added = False

    for structure_chain in current_chains:

        # If complex=True, the structure does not have enough chains. If in_a_branch=True, the program will have replaced a chain in the previous level.
        # The program will only test those chains that have been added in the previous level of recursion.
        if (complex and structure_chain.id[-1] != str(recursion_level-1)) or (in_a_branch and structure_chain.id[-1] != str(recursion_level)):
            continue

        # Iterates through all the pairwise interaction files
        for file in all_files:

            moving_chain = None     # Chain different to the receptor chain of the structure
            common_chain = None     # Chain homolog to the receptor chain of the structure

            # If the file has a chain with the same ID as the structure chain:
            if structure_chain.id[1] in [x[1] for x in info_files[file]]:   

                # Converts the pairwise interaction to a structure
                structure_to_add = parser.get_structure("structure", file) 

                # Translates the chain ID of the file to the same format as the ID of the current chain
                for chain in structure_to_add.get_chains():
                    original_id = chain.id
                    chain.id = [x for x in info_files[file] if x[0] == original_id][0]

                # If both chains in the file have the same ID, it contains a homodimer
                if info_files[file][0][1] == info_files[file][1][1]:

                    # Each chain can either be a moving or a common chain
                    moving_chains = [structure_to_add[0][info_files[file][0]], structure_to_add[0][info_files[file][1]]]
                    common_chains = [structure_to_add[0][info_files[file][1]], structure_to_add[0][info_files[file][0]]]

                else:   # If the file contains an heterodimer one will be rotated and the other superimposed

                    for chain in structure_to_add.get_chains():

                        if structure_chain.id[1] == chain.id[1]:
                            common_chains = [chain]
                        else:
                            moving_chains = [chain]

                for i in range(len(common_chains)):

                    # Copy the chains so the originals are not changed
                    if len(common_chains) > 1:
                        moving_chain = moving_chains[i].copy()
                        common_chain = common_chains[i].copy()
                    else:
                        moving_chain = moving_chains[i]
                        common_chain = common_chains[i]

                    RMSD, moved_chain = superimpose(structure_chain, common_chain, moving_chain)    # Superimpose the common chains and rotate the moving one
                    clashing_chains, same_chain= check_clash(current_structure, moved_chain)        # Look for clashes

                    # Added indicates wether a new chain has been added or not
                    added = False
                    if not clashing_chains and RMSD <= 3.0: # No clashes and good superimposition
                        
                        partial_id = moved_chain.id
                        existing_names = [x.id for x in current_structure[0].get_chains()]
                        
                        while not added:
                        
                            random_id = (id_generator(6),str(recursion_level)) # Assigns a random ID and stores the recursion level
                        
                            if partial_id + random_id not in existing_names:
                                moved_chain.id = tuple(list(partial_id) + list(random_id))
                                current_structure[0].add(moved_chain)
                                added = True
                                something_added = True


                    # Clashing chains, but more than one model is required or stoichiometry must be fulfilled 
                    if clashing_chains and not same_chain and (all_models or num_models > 1 or stoich):

                        open_branch = False

                        if in_a_branch:

                            clash_ids = set([(chain, moved_chain.id[1]) for chain in clashing_chains])
                            if len(non_brancheable_clashes.intersection(clash_ids)) == 0:
                                open_branch = True

                        else:
                            open_branch = True

                        if open_branch:

                            added_chain = moved_chain.copy()

                            # New ID containing random characters and the recursion level
                            added_chain.id = tuple(list(added_chain.id) + [id_generator(6), str(recursion_level)])

                            # New branch is generated
                            branch_structure = current_structure.copy()

                            # The conflicting chains are removed
                            for chain in clashing_chains:
                                branch_structure[0].detach_child(chain)
                                non_brancheable_clashes.add((added_chain.id, chain[1]))

                            # The chain that was clashing is added
                            branch_structure[0].add(added_chain)

                            # If the current structure has not already been tried:
                            if check_structure_exists(branch_structure, tried_branch_structures) is False:

                                logging.info('Chain %s in file %s conflicts with the rest of the structure. Opening a new branch.' % (added_chain.id[1], file))

                                # Branch added to the list of explored branch
                                tried_branch_structures.append(branch_structure)

                                # Go to the next recursion level
                                stored_structures = create_model(branch_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = recursion_level, in_a_branch = True, non_brancheable_clashes=non_brancheable_clashes,
                                    tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

                                # If the number of models is fulfilled
                                if all_models is False and num_models == len(stored_structures):
                                    return stored_structures


    if something_added and len(list(current_structure.get_chains())) < num_chains:

        if in_a_branch:
            stored_structures = create_model(current_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = recursion_level, complex = True, in_a_branch = True, non_brancheable_clashes=non_brancheable_clashes,
                tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)
        else:
            stored_structures = create_model(current_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = recursion_level, complex = True, in_a_branch = False, non_brancheable_clashes=non_brancheable_clashes,
                tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

        if all_models is False and num_models == len(stored_structures):
            return stored_structures

    else:

        # Check if stoichiometry is fulfilled
        if stoich:
            final_stoich = {}
            for chain in current_structure.get_chains():
                chain_id = chain.id[0]

                if chain_id not in final_stoich:
                    final_stoich[chain_id] = 1
                else:
                    final_stoich[chain_id] += 1                 

            divisors = set()
            for key, value in final_stoich.items():
                try:
                    divisor = value / stoich[key]
                    divisors.add(divisor)
                
                except KeyError:
                    logging.error('Stoichiometry file not valid, could not found chain %s' % key)
                    quit()

            if len(final_stoich) == len(stoich) and len(divisors) == 1 and check_structure_exists(current_structure, stored_structures) is False:
                stored_structures.append(current_structure)

                logging.info('Saving model.')

        elif check_structure_exists(current_structure, stored_structures) is False:

            stored_structures.append(current_structure)
            logging.info('Saving model.')

    return stored_structures