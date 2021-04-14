#!/usr/bin/env python
# coding=utf-8


from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import NeighborSearch
from Bio.PDB import Structure as pdb_struct
from Bio.PDB import Model as pdb_model
from Bio.PDB import PDBIO
import random
from process_input_files import id_generator

ppb = PPBuilder()
parser = PDBParser(PERMISSIVE = True, QUIET = True)

def superimpose(chainA, chainB, moving = None):

    """Superimpose one chain of a structure with an equivalent chain  chains and add another with the rotation parameters obtained. Return structure object with added chain, information about clashes and a flag for having added something.
    Keyword arguments:
    eq_chain1 -- common chain in the current structure (curr_struct)
    eq_chain2 -- common chain in the structure from which a chain wants to be added (moving_chain)
    moving_chain -- chain that may be added to the current complex
    rec_level_complex -- recursion level of building the complex
    filename2 -- name of the file that contains the moving_chain """

    residuesA = list(chainA.get_residues())
    residuesB = list(chainA.get_residues())

    # filter non-standard residues or errors
    filteredA = [res for res in residuesA if ('CA' in res or 'P' in res)]
    filteredB = [res for res in residuesB if ('CA' in res or 'P' in res)]

    # superimposer will try to pair atoms to move them from the position in one chain to the other, so it is important to have the same number of atoms
    if len(filteredA) > len(filteredB):
        filteredA = filteredA[:len(filteredB)]
    elif len(filteredA) < len(filteredB):
        filteredB = filteredB[:len(filteredA)]

    atomsA = []
    for res in filteredA:
        for atom in res.get_atoms():
            if atom.id == 'CA' or atom.id == 'P':
                atomsA.append(atom)
                break

    atomsB = []
    for res in filteredB:
        for atom in res.get_atoms():
            if atom.id == 'CA' or atom.id == 'P':
                atomsB.append(atom)
                break

    
    superimposer = Superimposer()

    superimposer.set_atoms(atomsA, atomsB)

    rmsd = superimposer.rms

    if moving:
        superimposer.apply(list(moving.get_atoms()))
        return rmsd, moving

    return rmsd    


def check_clash(model, chain_to_add, clash_distance=2.5):

    """Check if there is a steric clash between a rotating chain and current structure. Return False (no clash), 1 (clash between two different chains) or 2 (same chain). Also returns the ids of the chains in structure that are clashing
    Keyword arguments:
    structure -- whole structure
    rotating_chain -- chain to be analyzed with respect to structure
    distance_for_clash -- threshold to consider clash between atoms. Default = 2.5 (Armstrongs)
    Clash criteria: at least 20 of the atoms are at a lower distance than distance_for_clash
    Same chain criteria: RMSD between them <= 3.0 """

    # initialize the neighbor search
    neighbor_object = NeighborSearch(list(model.get_atoms()))

    structure_clashing_chains = set() # It is important that it is a set because no duplicates
    total_clashes = 0

    for atom in chain_to_add.get_atoms():
        clashes = neighbor_object.search(atom.get_coord(), clash_distance)
        if len(clashes) > 0:
            for clash in clashes:
                structure_clashing_chains.add(clash.get_parent().get_parent().id)
                total_clashes += 1 # Number of atoms that clash
    
    if len(structure_clashing_chains) > 1 and total_clashes > 20:
        # a clash against different chains:
        return structure_clashing_chains, False
    
    elif len(structure_clashing_chains) == 1 and total_clashes > 20 and chain_to_add.id[1] == list(structure_clashing_chains)[0][1]:

        # a clash MAYBE because you are trying to superimpose something in the place it was already

        # define the clashing chain:
        clash_chain = model[0][list(structure_clashing_chains)[0]]
        RMSD = superimpose(clash_chain, chain_to_add)

        if RMSD <= 3.0:
            # it is the same chain
            return structure_clashing_chains, True

        else:
            # it is another chain or the same with different structure
            return structure_clashing_chains, False

    elif total_clashes > 20:
        # it is ine chain and the previous conditions are not fullfilled
        return structure_clashing_chains, False

    else:
        # no clash
        return None, False

def check_structure_exists(structure, all_structures):

    """Ask if structure is already in created_structures (a list of structures). Returns a boolean.
    Considerations:
    Return True if all of the chains in structure are in one of the structures in created_structures and RMSD <= 3.0, meaning they are the same structure """

    # make a deepcopy of these objects
    structure = structure.copy()
    all_structures = all_structures.copy()

    # Get the ids of the chains in structure
    chain_ids_structure = tuple(sorted([x.id[1] for x in structure.get_chains()]))

    # loop through each of the contents of created_structures:
    for each_structure in all_structures:

        # get the chains of created_structure
        chain_ids_each_structure = tuple(sorted([x.id[1] for x in each_structure.get_chains()]))

        # ask if the number of each and ids of the chains are the same:
        if chain_ids_structure == chain_ids_each_structure:

            # pick one chain in structure to compare with the chains in created
            chain_str = list(structure.get_chains())[0]
            id_str = chain_str.id[1]

            # try to find a partner in created_structure:
            for chain_each_str in each_structure.get_chains():

                id_each_str = chain_each_str.id[1]

                # if they have the same id they are potential partners. Superimpose these next. The id_created_str has also to be avaliable in possible_partners
                if id_str == id_each_str:

                    RMSD, structure = superimpose(chain_str, chain_each_str, moving=structure)

                    # if I have superimposed same ID but different structure, try another chain
                    if  RMSD > 3.0:
                        continue

                    # if the previous chain_str and chain_created_str are real partners they should also result in haveing all they cross-superimposed chains with partners
                    partners = set()

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

                        if len(partners) == len(list(created_structure.get_chains())):
                            return True  # all chains have a partner, which means that the structure is in the created_structures

    # if you didn't find any match return false:
    return False    

def create_model(current_structure, stored_structures, info_files, num_chains, num_models, exhaustive, tree_level = 0, climb = False, in_a_branch = False, non_brancheable_clashes=set(), tried_branch_structures=list(), stoich=None, verbose=False):

    """Build model
    AÑADIR MÁS INFO"""

    # 
    if exhaustive is False and num_models == len(stored_structures):
        return stored_structures

    if climb:
        tree_level += 1

    # files that have to be used. The shuffle is for generating alternative paths, so that running the program many times may generate slightly different resukts
    all_files= list(info_files.keys())
    random.shuffle(all_files)


    # # a boolean that indicates if there's something added at this level
    something_added = False

    # #    POSIBLES MEJORAS print the conformation of the current complex:
    # if verbose:
    #     print_topology_of_complex(current_str)

    current_chains = list(current_structure.get_chains())

    # iterate through the chains of the current structure
    for structure_chain in current_chains:

        # when I climb the tree, I only want to add chains to those added in the previous node
        # If I am in a branch, I have not added but replaced a chain. Therefore, I still want to add a chain in the same node
        # This way we are not doing not necessary comparison

        if (climb and structure_chain.id[-1] != str(tree_level-1)) or (in_a_branch and structure_chain.id[-1] != str(tree_level)):
            continue

        # if verbose:
        #     print('Trying to add chains through the common chain: ', chain1.id)

        # iterate through the PDB files of the directory to compare them to the current struct
        for file in all_files:

            # initialize the common and rotating chain
            moving_chain = None
            equal_chain = None

            ## AÑADIR AQUÍ QUE 

            # if the chain id of the current structure is found in the second PDB file:
            if structure_chain.id[1] in [x[1] for x in info_files[file]]:

                # get the structure
                structure_to_add = parser.get_structure("structure", file) 


                # change the chain id to our format
                for chain in structure_to_add.get_chains():
                    original_id = chain.id
                    chain.id = [x for x in info_files[file] if x[0] == original_id][0]
                    #chain.transform(np.array([[0.36, 0.48, -0.8], [-0.8, 0.6, 0], [0.48, 0.64, 0.6]]), np.array([2, 2, 2]))


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
                    if len(moving_chains) > 1:
                        moving_chain = moving_chains[i].copy()
                        common_chain = common_chains[i].copy()
                    else:
                        moving_chain = moving_chains[i]
                        common_chain = common_chains[i]

                    RMSD, moved_chain = superimpose(structure_chain, common_chain, moving_chain)
                    clashing_chains, same_chain= check_clash(current_structure, moved_chain)

                    # something is added if there's no clashes and the RMSD is very low, indicating that the two chains are actually the same
                    added = False
                    if not clashing_chains and RMSD <= 3.0:
                        partial_id = moved_chain.id
                        existing_names = [x.id for x in current_structure[0].get_chains()]
                        while not added:
                            random_id = (id_generator(6),str(tree_level)) # random ID + a number that indicates the recursion level at which this chain has been added
                            if partial_id + random_id not in existing_names:
                                moved_chain.id = tuple(list(partial_id) + list(random_id))
                                current_structure[0].add(moved_chain)
                                added = True


                    # when there's a aberrant clash and you fulfill one of the branch-opening  conditions
                    if clashing_chains and not same_chain and (exhaustive or num_models > 1 or stoich):

                        # a branch complex will be created if the rotating chain is not one of the previously branch-opening clashing chains
                        # or if the clashes happen against the chain that opened this branch

                        open_branch = False

                        if in_a_branch:

                            clash_ids = set([(x, moved.id[1]) for x in clashing_chains])
                            # check if the clashes are not an exception
                            if len(non_brancheable_clashes.intersection(clash_ids)) == 0:
                                open_branch = True

                        else:
                            open_branch = True

                        if open_branch:

                            # set the id of the chain you were trying to add:
                            added_chain = moved_chain.copy()
                            added_chain.id = tuple(list(added_chain.id) + [id_generator(6), str(tree_level)])

                            # make a new branch:
                            branch_structure = current_structure.copy()

                            # remove the chains that chains that are clashing:
                            for clashing_chain in clashing_chains:
                                branch_structure[0].detach_child(clashing_chain)

                            # create a set of tuples (added_chain, clashing_chains_ids) of branch opening exceptions
                            for x in clashing_chains:
                                non_brancheable_clashes.add((added_chain.id, x[1]))

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
                                stored_structures = create_model(branch_structure, stored_structures, info_files, num_chains, num_models, exhaustive, tree_level = tree_level, in_a_branch = True, non_brancheable_clashes=non_brancheable_clashes,
                                    tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

                                # return as soon as possible:
                                if exhaustive is False and num_models == len(stored_structures):
                                    return stored_structures

    if added and len(list(current_structure.get_chains())) <= num_chains:

        if in_a_branch:
            stored_structures = create_model(current_structure, stored_structures, info_files, num_chains, num_models, exhaustive, tree_level = tree_level, climb = True, in_a_branch = True, non_brancheable_clashes=non_brancheable_clashes,
                tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)
        else:
            stored_structures = create_model(current_structure, stored_structures, info_files, num_chains, num_models, exhaustive, tree_level = tree_level, climb = True, in_a_branch = False, non_brancheable_clashes=non_brancheable_clashes,
                tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)


        # return as soon as possible:
        if exhaustive is False and num_models == len(stored_structures):
            return stored_structures

    else:
        if verbose:
            print("trying to save model")

        # check stoichiometry
        if stoich:
            final_stoich = {}
            for chain in current_structure.get_chains():
                chain_id = chain.id[1]

                if chain_id not in final_stoich:
                    final_stoich[chain_id] = 1
                else:
                    final_stoich[chain_id] += 1

            divisors = set()
            for key, value in final_stoich.items():
                divisor = value / stoich[key]
                divisors.add(divisor)

            if len(divisors) == 1 and check_structure_exists(current_structure, stored_structures) is False:
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

