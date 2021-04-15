#!/usr/bin/env python
# coding=utf-8

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
                        if len(partners) == len(list(created_structure.get_chains())):
                            return True  

    # If not all chains match returns False
    return False    

def create_model(current_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = 0, complex = False, in_a_branch = False, non_brancheable_clashes=set(), tried_branch_structures=list(), stoich=None, verbose=False):

    """Build model
    AÑADIR MÁS INFO"""

    # 
    if all_models is False and num_models == len(stored_structures):
        return stored_structures

    if complex:
        recursion_level += 1

    # files that have to be used. The shuffle is for generating alternative paths, so that running the program many times may generate slightly different resukts
    all_files= list(info_files.keys())
    random.shuffle(all_files)


    # #    POSIBLES MEJORAS print the conformation of the current complex:
    # if verbose:
    #     print_topology_of_complex(current_str)

    current_chains = list(current_structure.get_chains())

    # iterate through the chains of the current structure
    something_added = False

    for structure_chain in current_chains:

        # If complex = true means that there are not enough chains and I need to add more. But I don't want to add chains to those which I have already tested. Therefore; I am only going to use those chains added at previous recursion
        # If I am in a branch, I have not added but replaced a chain. Therefore, I still want to add a chains to those on the same recursion level
        # This way we are not doing not necessary comparison

        if (complex and structure_chain.id[-1] != str(recursion_level-1)) or (in_a_branch and structure_chain.id[-1] != str(recursion_level)):
            continue

        # if verbose:
        #     print('Trying to add chains through the common chain: ', chain1.id)

        # iterate through the PDB files of the directory to compare them to the current struct
        for file in all_files:

            # initialize the common and rotating chain
            moving_chain = None
            common_chain = None

            ## AÑADIR AQUÍ QUE 

            # if the chain id of the current structure is found in the second PDB file:
            if structure_chain.id[1] in [x[1] for x in info_files[file]]:

                # get the structure
                structure_to_add = parser.get_structure("structure", file) 

                # change the chain id to our format
                for chain in structure_to_add.get_chains():
                    original_id = chain.id
                    chain.id = [x for x in info_files[file] if x[0] == original_id][0]

                # if it is a homodimer set both pssibilities of rotating / common

                if info_files[file][0][1] == info_files[file][1][1]:

                    # generate an array with the two possibilities
                    moving_chains = [structure_to_add[0][info_files[file][0]], structure_to_add[0][info_files[file][1]]]
                    common_chains = [structure_to_add[0][info_files[file][1]], structure_to_add[0][info_files[file][0]]]

                else:
                    # decide which is (rotating / common)
                    for chain in structure_to_add.get_chains():

                        if structure_chain.id[1] == chain.id[1]:
                            common_chains = [chain]
                        else:
                            moving_chains = [chain]

                # go through each possible rotating and common chains

                for i in range(len(common_chains)):

                    # define each of the rotating/common, without perturbing the original ones
                    if len(common_chains) > 1:
                        moving_chain = moving_chains[i].copy()
                        common_chain = common_chains[i].copy()
                    else:
                        moving_chain = moving_chains[i]
                        common_chain = common_chains[i]

                    RMSD, moved_chain = superimpose(structure_chain, common_chain, moving_chain)
                    clashing_chains, same_chain= check_clash(current_structure, moved_chain)

                    # something is added if there's no clashes and the RMSD is very low, indicating that the two chains compared are actually the same
                    added = False
                    if not clashing_chains and RMSD <= 3.0:
                        partial_id = moved_chain.id
                        existing_names = [x.id for x in current_structure[0].get_chains()]
                        while not added:
                            random_id = (id_generator(6),str(recursion_level)) # random ID + a number that indicates the recursion level at which this chain has been added
                            if partial_id + random_id not in existing_names:
                                moved_chain.id = tuple(list(partial_id) + list(random_id))
                                current_structure[0].add(moved_chain)
                                added = True
                                something_added = True


                    # when there's a aberrant clash and you fulfill one of the branch-opening  conditions
                    if clashing_chains and not same_chain and (all_models or num_models > 1 or stoich):

                        # a branch complex will be created if the rotating chain is not one of the previously branch-opening clashing chains
                        # or if the clashes happen against the chain that opened this branch

                        open_branch = False

                        if in_a_branch:

                            clash_ids = set([(chain, moved_chain.id[1]) for chain in clashing_chains])
                            # check if the clashes are not an exception
                            if len(non_brancheable_clashes.intersection(clash_ids)) == 0:
                                open_branch = True

                        else:
                            open_branch = True

                        if open_branch:

                            # set the id of the chain you were trying to add:
                            added_chain = moved_chain.copy()
                            added_chain.id = tuple(list(added_chain.id) + [id_generator(6), str(recursion_level)])

                            # make a new branch:
                            branch_structure = current_structure.copy()

                            # remove the chains that chains that are clashing:
                            for chain in clashing_chains:
                                branch_structure[0].detach_child(chain)
                                non_brancheable_clashes.add((added_chain.id, chain[1]))

                            # add the chain that was clashing:
                            branch_structure[0].add(added_chain)

                            # check if the branch_new_str is not in tried_branch_structures:
                            if check_structure_exists(branch_structure, tried_branch_structures) is False:

                                # indicate that a branch is opening
                                if verbose:
                                    print('%s of %s  clashes against the complex. Opening a new branch...' % (added_chain.id[1], file))

                                # add this branch to the ones already tested:
                                tried_branch_structures.append(branch_structure)

                                # create a new structure based on this branch:
                                stored_structures = create_model(branch_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = recursion_level, in_a_branch = True, non_brancheable_clashes=non_brancheable_clashes,
                                    tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

                                # return as soon as possible:
                                if all_models is False and num_models == len(stored_structures):
                                    return stored_structures


    if something_added and len(list(current_structure.get_chains())) < num_chains:

        if in_a_branch:
            stored_structures = create_model(current_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = recursion_level, complex = True, in_a_branch = True, non_brancheable_clashes=non_brancheable_clashes,
                tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)
        else:
            stored_structures = create_model(current_structure, stored_structures, info_files, num_chains, num_models, all_models, recursion_level = recursion_level, complex = True, in_a_branch = False, non_brancheable_clashes=non_brancheable_clashes,
                tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)


        # return as soon as possible:
        if all_models is False and num_models == len(stored_structures):
            return stored_structures

    else:
        
        print("trying to save model")
        print("this structure has already been saved")
        print(check_structure_exists(current_structure, stored_structures))

        # check stoichiometry
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
                divisor = value / stoich[key]
                divisors.add(divisor)

            if len(final_stoich) == len(stoich) and len(divisors) == 1 and check_structure_exists(current_structure, stored_structures) is False:
                stored_structures.append(current_structure)

                # if verbose:
                #     print("saving model")
                #     print_topology_of_complex(current_str)
                #     print('\n\n')

        elif check_structure_exists(current_structure, stored_structures) is False:

            stored_structures.append(current_structure)
            # if verbose:
            #     print("saving model")
            #     print_topology_of_complex(current_str)
            #     print('\n')

    return stored_structures