#!/usr/bin/env python
# coding=utf-8


from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import NeighborSearch
from Bio import pairwise2
ppb = PPBuilder()
parser = PDBParser(PERMISSIVE = True, QUIET = True)
import string
import random

def id_generator (n, chars=string.ascii_uppercase + string.digits):
    """Generates a random ID of n number of elements."""
    return ''.join(random.choice(chars) for _ in range(n))

def process_input (list_files):
    """Process input generating two dictionaries from a list of PDB files.
    The info_seqs dictionary contains unique sequences as values and the key is the new ID assigned to this unique sequence.
    The info_files dictionary has, as key, the name of the file and as value a list of tuples (one for each of the chains)
    with the original chain ID found in the file and the new ID assigned to this unique sequence."""

    info_files = {}
    info_seqs = {}

    # loop through files found in the list of files
    for filename in list_files:

        structure = parser.get_structure("structure", filename) # get the structure of the PDB file

        n_chains = len(structure[0].child_list)

        if n_chains != 2:
            # ESTE ERROR INCLUYE EL ARCHIVO EN INFO_FILES? HABRÍA QUE QUITARLO
            logging.error("There should be 2 chains in the pdb file. %d chains found in %s instead." % (n_chains, filename))

        chains = list(structure.get_chains())

        match_IDs = [] # it will contain a tuple for each of the chains of the file with the original ID and the new ID assigned to it

        # loop through chains of the current file
        for chain in chains:
            print(chain.id)
            sequence = "".join([str(pp.get_sequence()) for pp in ppb.build_peptides(chain)])

            # if sequence is still empty it may be because it is a DNA sequence instead of a nucleotide one
            if sequence  == "":
                sequence = "".join([nt.resname[2] for nt in chain.get_residues() if nt.resname[2] in 'ACGTU' and len(nt.resname) == 3])

            ## HABRÍA QUE COMPROBAR QUE LA SECUENCIA YA NO ESTÁ VACÍA

            # assign a unique ID to the sequence
            seq_id = None
            identity_threshold = 95 # identity threshold for which above it sequences will be considered the same
            highest_identity = identity_threshold

            # if info_seqs is not empty, loop through each of the keys to align the current sequence with the ones stored in this dict.
            if info_seqs:
                for chain_id in info_seqs:

                    # align my sequence with all the ones stored in the info_seqs
                    alignment = pairwise2.align.globalxx(sequence, info_seqs[chain_id])

                    # calculate the percentage identity in the alignment
                    match_percent = (alignment[0][2] / min(len(sequence), len(info_seqs[chain_id]))) * 100

                    if match_percent > identity_threshold:
                        seq_id = chain_id # assign the same ID
                        highest_identity = match_percent # update the highest identity

                # if sequences are not considered the same, assign a new ID to it and make sure this ID as a key is unique in the info_seqs dict
                if seq_id is None:
                    n = 0
                    while n < 1:
                        seq_id = id_generator(6)
                        if seq_id not in info_seqs:
                            n = 1
                            info_seqs[seq_id] = sequence
                        else:
                            n = 0

            # if the info_seqs is still empty, add the sequence and generate a random ID to it
            else:
                seq_id = id_generator(6)
                info_seqs[seq_id] = sequence

            # create the tuple of original ID and new ID for each of the chains of the file
            ids = (chain.id, seq_id)
            match_IDs.append(ids)

        # save info of each file in the files_info dictionary
        info_files[filename] = match_IDs

    return info_files, info_seqs


def superimpose(chainA, chainB, moving_chain = None):

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

    if moving_chain:
        superimposer.apply(list(moving_chain.get_atoms()))
        return rmsd, moving_chain

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
    neighbor_object = NeighborSearch(list(structure.get_atoms()))

    structure_clashing_chains = set() # It is important that it is a set because no duplicates
    total_clashes = 0

    # CAMBIAR A: (PARA QUE SE PAREZCA MENOS A MIKI)   
    # clashes = 0
    # length = 0
    # for atom in new_chain_atoms:
    #     length +=1 
    #     if neighbor_object.search(atom.coord, 2):
    #         clashes += 1
    
    # proportion_of_clashes = float(clashes/length)
    
    # if proportion_of_clashes < clash_perc:

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
        clash_chain = structure[0][list(structure_clashing_chains)[0]]
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
        chain_ids_each_structure = tuple(sorted([x.id[1] for x in created_structure.get_chains()]))

        # ask if the number of each and ids of the chains are the same:
        if chain_ids_structure == chain_ids_each_structure:

            # pick one chain in structure to compare with the chains in created
            chain_str = list(structure.get_chains())[0]
            id_str = chain_str.id[1]

            # try to find a partner in created_structure:
            for chain_created_str in created_structure.get_chains():

                id_created_str = chain_created_str.id[1]

                # if they have the same id they are potential partners. Superimpose these next. The id_created_str has also to be avaliable in possible_partners
                if id_str == id_created_str:

                    RMSD, structure = superimpose(chain_str, chain_created_str, structure)

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
                                    # get list of residues:
                                    res_partner1 = list(searching_partner.get_residues())
                                    res_partner2 = list(possible_partner.get_residues())

                                    # get the atoms of the previous list, ONLY belonging to common RESIDUES! to be then able to superimpose
                                    # so first we obtain a list of the common residues
                                    common_res_p1 = get_list_of_common_res(res_partner1, res_partner2)
                                    common_res_p2 = get_list_of_common_res(res_partner2, res_partner1)

                                    # then we obtain a list of coordinates
                                    common_coords_p1 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_p1)])
                                    common_coords_p2 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_p2)])

                                    rms = rmsd.kabsch_rmsd(common_coords_p2, common_coords_p1)
                                    if rms <= 3.0:
                                        partners.add(possible_partner)
                                        partner_found = True

                        if len(partners) == len(list(created_structure.get_chains())):
                            return True  # all chains have a partner, which means that the structure is in the created_structures

    # if you didn't find any match return false:
    return False    

def create_model(current_structure, stored_structures, info_files, num_chains, num_models, tree_level = 0, climb = False, exhaustive, in_a_branch = False, non_brancheable_clashes=set(), tried_branch_structures=[], stoich=None, verbose=False):

    """Build model
    AÑADIR MÁS INFO"""

    # 
    if exhaustive is False and num_models == len(stored_structures):
        return stored_structures

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

        # when you are in deeper recursion levels check that the chains are the ones of the previous level, for avoiding many comparisons
        interesting_complex_rec_level = None

        if this_is_a_complex_recursion:
            # we are interested in analyzing chains that have been added JUST in the previous recursion levels
            interesting_complex_rec_level = rec_level_complex - 1

        elif this_is_a_branch:
            # we are interested in analyzing chains of the current rec_level_complex
            interesting_complex_rec_level = rec_level_complex

        # remember that chain1.id will be atuple of ONLY two elements if rec_level_complex=0
        if rec_level_complex > 0 and chain1.id[-1] != str(interesting_complex_rec_level):
            continue

        if verbose:
            print('Trying to add chains through the common chain: ', chain1.id)

        id_chain1 = chain1.id[1]

        # iterate through the PDB files of the directory to compare them to the current struct
        for file in all_files:

            # initialize the common and rotating chain
            rotating_chain = None
            common_chain2 = None

            ## AÑADIR AQUÍ QUE 

            # if the chain id of the current structure is found in the second PDB file:
            if id_chain1 in [x[1] for x in PDB_dict[filename2]]:

                # get the structure
                structure2 = p.get_structure("pr2", mydir + filename2)

                if len(list(structure2.get_chains())) != 2:
                    raise TwoChainException(filename2)

                # change the chain id names
                for chain in structure2.get_chains():
                    curr_id = chain.id  # the so-called chain accession, usually a letter (A,B,C,D ...)
                    chain.id = [x for x in PDB_dict[filename2] if x[0] == curr_id][0]

                # if it is a homodimer set both pssibilities of rotating / common
                if PDB_dict[filename2][0][1] == PDB_dict[filename2][1][1]:

                    # generate an array with the two possibilities

                    rotating_chains = [structure2[0][PDB_dict[filename2][0]], structure2[0][PDB_dict[filename2][1]]]
                    common_chains2 = [structure2[0][PDB_dict[filename2][1]], structure2[0][PDB_dict[filename2][0]]]

                else:
                    # decide which is (rotating / common)
                    for chain2 in structure2.get_chains():
                        id_chain2 = chain2.id[1]

                        if id_chain1 == id_chain2:
                            common_chain2 = chain2
                        else:
                            rotating_chain = chain2

                    # define the array
                    rotating_chains = [rotating_chain]
                    common_chains2 = [common_chain2]

                # go through each possible rotating and common chains

                for I in range(0, len(rotating_chains)):

                    # define each of the rotating/common, without perturbing the original ones
                    if len(rotating_chains) > 1:
                        rotating_chain = rotating_chains[I].copy()
                        common_chain2 = common_chains2[I].copy()
                    else:
                        rotating_chain = rotating_chains[I]
                        common_chain2 = common_chains2[I]

                    # try to add the new chain to the complex


                    RMSD, moving_chain = superimpose(chain1, common_chain2, rotating_chain)
                    clashing_chains, same_chain= is_steric_clash(curr_struct, moving_chain)

                    # something is added if there's no clashes and the RMSD is very low, indicating that the two chains are actually the same
                    added = False
                    if not clashing_chains and rms <= 3.0:
                        my_id = moving_chain.id
                        chain_names = [x.id for x in curr_struct[0].get_chains()]
                        while not added:
                            rand = (create_random_chars_id(6),str(rec_level_complex)) # random ID + a number that indicates the recursion level at which this chain has been added
                            if my_id + rand not in chain_names:
                                moving_chain.id = tuple(list(my_id) + list(rand))
                                curr_struct[0].add(moving_chain)
                                added = True


                    # when there's a aberrant clash and you fulfill one of the branch-opening  conditions
                    if clashing_chains and not same_chain and (exhaustive or num_models > 1 or stoich):

                        # a branch complex will be created if the rotating chain is not one of the previously branch-opening clashing chains
                        # or if the clashes happen against the chain that opened this branch

                        open_branch = False

                        if this_is_a_branch:

                            clash_ids = set([(x, added_chain.id[1]) for x in clashing_chains])

                            # check if the clashes are not an exception
                            if len(non_brancheable_clashes.intersection(clash_ids)) == 0:
                                open_branch = True

                        else:
                            open_branch = True

                        if open_branch:

                            # set the id of the chain you were trying to add:
                            added_chain = added_chain.copy()
                            added_chain.id = tuple(list(added_chain.id) + [create_random_chars_id(6), str(rec_level_complex)])

                            # make a new branch:
                            branch_new_str = current_str.copy()

                            # remove the chains that chains that are clashing:
                            for clashing_chain in clashing_chains:
                                branch_new_str[0].detach_child(clashing_chain)

                            # create a set of tuples (added_chain, clashing_chains_ids) of branch opening exceptions
                            for x in clashing_chains:
                                non_brancheable_clashes.add((added_chain.id, x[1]))

                            # add the chain that was clashing:
                            branch_new_str[0].add(added_chain)

                            # check if the branch_new_str is not in tried_branch_structures:
                            if structure_in_created_structures(branch_new_str, tried_branch_structures) is False:

                                # indicate that a branch is opening
                                if verbose:
                                    print('%s of %s  clashes against the complex. Opening a new branch...' % (added_chain.id[1], filename2))

                                # add this branch to the ones already tested:
                                tried_branch_structures.append(branch_new_str)

                                # create a new structure based on this branch:
                                saved_models = build_complex(saved_models, branch_new_str, mydir, PDB_dict, num_models, exhaustive, n_chains, this_is_a_branch=True, non_brancheable_clashes=non_brancheable_clashes,
                                              rec_level_complex=rec_level_complex, rec_level_branch=rec_level_branch, tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

                                # return as soon as possible:
                                if exhaustive is False and num_models == len(saved_models):
                                    return saved_models

    if something_added and len(list(current_str.get_chains())) <= n_chains:

        if this_is_a_branch:
            set_this_is_a_branch = True
        else:
            set_this_is_a_branch = False

        saved_models = build_complex(saved_models, current_str, mydir, PDB_dict, num_models, exhaustive, n_chains, this_is_a_complex_recursion=True, this_is_a_branch=set_this_is_a_branch, non_brancheable_clashes=non_brancheable_clashes,
            rec_level_complex=rec_level_complex, rec_level_branch=rec_level_branch, tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

        # return as soon as possible:
        if exhaustive is False and num_models == len(saved_models):
            return saved_models

    else:
        if verbose:
            print("trying to save model")

        # check stoichiometry
        if stoich:
            final_stoich = {}
            for chain in current_str.get_chains():
                chain_id = chain.id[1]

                if chain_id not in final_stoich:
                    final_stoich[chain_id] = 1
                else:
                    final_stoich[chain_id] += 1

            divisors = set()
            for key, value in final_stoich.items():
                divisor = value / stoich[key]
                divisors.add(divisor)

            if len(divisors) == 1 and structure_in_created_structures(current_str, saved_models) is False:
                saved_models.append(current_str)

                if verbose:
                    print("saving model")
                    print_topology_of_complex(current_str)
                    print('\n\n')

        elif structure_in_created_structures(current_str, saved_models) is False:

            saved_models.append(current_str)
            if verbose:
                print("saving model")
                print_topology_of_complex(current_str)
                print('\n')

    return saved_models
