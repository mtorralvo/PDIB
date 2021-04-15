"""
Contains functions for converting the input to a set of dictionaries, containing information about the files,
the sequences in each file and the required stoichiometry. It also generates random IDs so there is no conflict between
chain names.
"""

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2
import random
from tempfile import NamedTemporaryFile
from tempfile import TemporaryDirectory
import string
import sys
import logging
import os
import re


def id_generator(n, chars=string.ascii_uppercase + string.digits):
    """Generates a random ID of n number of elements."""
    return ''.join(random.choice(chars) for _ in range(n))

def process_input (list_files):
    """
    Processes the input generating two dictionaries from a list of PDB files.
    The info_seqs dictionary contains unique sequences as values and the key is the new ID assigned to this unique sequence.
    The info_files dictionary has, as key, the name of the file and as value a list of tuples (one for each of the chains)
    with the original chain ID found in the file and the new ID assigned to this unique sequence.
    """

    info_files = {}
    info_seqs = {}
    parser = PDBParser(PERMISSIVE=1)
    ppb = PPBuilder()
    check_user = 0
    dir_created = 0

    # Loop through files found in the list of files
    for filename in list_files:

        structure = parser.get_structure('structure', filename) # get the structure of the PDB file

        n_chains = len(structure[0].child_list)

        if n_chains != 2:
            if check_user == 0:
                user_allows = input('It seems that some files contain more than two chains. These may be caused by DNA or RNA sequences being represented as two chains.'
                    'PDIB will try and treat both chains as one, generating files in the process that will be later removed. '
                    'Do you want to continue? (y/n): ')
                check_user = 1
            
            if user_allows == 'y':

                if dir_created == 0:
                    os.makedirs('PDIB_tmp', exist_ok=True)
                    dir_created = 1
                
                dna_res = [' DA', ' DT', ' DC', ' DG']
                temp = NamedTemporaryFile(dir='PDIB_tmp', delete=False)
                flag = 0
                
                for line in open(filename):
                    
                    if line[17:20] in dna_res:
                        
                        if flag == 0:
                            dna_name = line[21]
                            flag = 1
                        
                        if line[21] == dna_name:
                            temp.write(line.encode())
                        
                        else:
                            temp.write(f'{line[0:21] + dna_name + line[22:]}'.encode())
                    
                    else:
                        temp.write(line.encode())
                
                filename = temp.name
                structure = parser.get_structure('structure', filename) # get the structure of the PDB file

            else:
                logging.error('There should be 2 chains in the pdb file. %d chains found in %s instead. Exiting the program.' % (n_chains, filename))
                quit()
        
        chains = list(structure.get_chains())

        match_IDs = [] # It will contain a tuple for each of the chains of the file with the original ID and the new ID assigned to it

        # Loop through chains of the current file
        for chain in chains:
            
            sequence = ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(chain)])

            # If sequence is still empty it may be because it is a DNA sequence instead of a nucleotide one
            if sequence  == '':
                sequence = ''.join([nt.resname[2] for nt in chain.get_residues() if nt.resname[2] in 'ACGTU' and len(nt.resname) == 3])

            # Assign a unique ID to the sequence
            seq_id = None
            identity_threshold = 95 # identity threshold for which above it sequences will be considered the same
            highest_identity = identity_threshold

            # If info_seqs is not empty, loop through each of the keys to align the current sequence with the ones stored in this dict.
            if info_seqs:
                for chain_id in info_seqs:

                    # Align my sequence with all the ones stored in the info_seqs
                    alignment = pairwise2.align.globalxx(sequence, info_seqs[chain_id])

                    # Calculate the percentage identity in the alignment
                    match_percent = (alignment[0][2] / min(len(sequence), len(info_seqs[chain_id]))) * 100

                    if match_percent > identity_threshold:
                        seq_id = chain_id # assign the same ID
                        highest_identity = match_percent # update the highest identity

                # If sequences are not considered the same, assign a new ID to it and make sure this ID as a key is unique in the info_seqs dict
                if seq_id is None:
                    n = 0
                    while n < 1:
                        seq_id = id_generator(6)
                        if seq_id not in info_seqs:
                            n = 1
                            info_seqs[seq_id] = sequence
                        else:
                            n = 0

            # If the info_seqs is still empty, add the sequence and generate a random ID to it
            else:
                seq_id = id_generator(6)
                info_seqs[seq_id] = sequence

            # Create the tuple of original ID and new ID for each of the chains of the file
            ids = (chain.id, seq_id)
            match_IDs.append(ids)

        # Save info of each file in the files_info dictionary
        info_files[filename] = match_IDs

    return info_files, info_seqs

def get_stoichiometry(path):
    """Converts the stoichiometry file to a dictionary"""
    
    st_dict = {}
    with open(path) as fh:
    
        for line in fh:
            if re.match('\s*\w+\s*(:|=)\s*\d+', line):
                splitted = re.split('(:|=)', line)
                st_dict[splitted[0].strip()] = int(splitted[2].strip())
            else:
                logging.warning("The stoichiometry file must follow the format 'Chain : number of repeats'. Each chain must be separated by a line break. Error found in line '%s'" % line.strip())
    
    return st_dict
  

