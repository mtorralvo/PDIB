import argparse
import re
import os
import logging
import sys
import glob
import random
import string

parser = argparse.ArgumentParser(description='Make protein origami and stuff.')


parser.add_argument( '-i', '--input-directory',
                    dest = "input_directory",
                    action = "store",
                    default = "./",
                    help = "Input directory containing files")

parser.add_argument('-s', '--stoichiometry',
                    type = str,
                    help = "Path containing the stechiometry of the complex")

parser.add_argument('-o', '--output-directory',
                    dest = 'output_directory',
                    action = "store",
                    default = "./",
                    required = True,
                    help = "Output directory to store the models")

parser.add_argument('-f', '--force',
                    dest = 'force',
                    action = 'store_true',
                    default = False,
                    help = "If it is True, the program can overwrite all the content of the output directory")

parser.add_argument('-v', '--verbose',
                    dest = 'verbose',
                    action = 'store_true',
                    default = False,
                    help = "Gives information during the execution of the script")

args = parser.parse_args() # list of arguments

if args.verbose:
    logging.basicConfig(stream=sys.stderr,
                        level=logging.DEBUG,
                        format = '%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')

regex = re.compile(".*\w+_\w+_\w+.pdb(.gz)*")


input_dir = args.input_directory #read input directory

# Check if input is a directory and if True, save files matching the regex in a list
if os.path.isdir(input_dir):
    list_files = [ os.path.join(input_dir, f) for f in os.listdir(input_dir) if regex.match(f)]

    if len(os.listdir(input_dir)) == 0:
        raise ValueError("The input folder is empty")

    elif len(list_files) == 0:
        raise ValueError("Input must be of the format <name>_<chain1>_<chain2>.pdb(.gz)")

elif os.path.isfile(input_dir):
    raise ValueError("Input must be a directory")

if args.verbose:
    sys.stderr.write("%d input files found.\n" %len(list_files))

# Generate the output folders in case they do not exist and the force option is activated
if not os.path.isdir(args.output_directory) or args.force:
    os.makedirs("./%s/structures" % args.output_directory, exist_ok = True)
    os.makedirs("./%s/analysis" % args.output_directory, exist_ok = True)

else:
    raise OSError("%s already exists, specify a different output directory or enable option -f to override the already existing one." % args.output_directory)

logging.info('doing something')

print(list_files)


from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2
ppb = PPBuilder()

def id_generator (n, chars=string.ascii_uppercase + string.digits):
    """Generates a random ID of n number of elements"""
    return ''.join(random.choice(chars) for _ in range(n))

def process_input (list_files):
    """Process input generating two dictionaries from a list of PDB files.
    The info_seqs dictionary contains unique sequences as values and the key is the new ID assigned to this unique sequence.
    The info_files dictionary has, as key, the name of the file and as value a list of tuples (one for each of the chains)
    with the original chain ID found in the file and the new ID assigned to this unique sequence."""

    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    info_files = {}
    info_seqs = {}

    # loop through files found in the list of files
    for filename in list_files:
        structure = parser.get_structure("Name", filename) # get the structure of the PDB file

        n_chains = len(structure[0].child_list)

        if n_chains != 2:
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

    print(info_seqs)
    #print(info_files)
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
                logging.warning('The stoichiometry file must follow the format "Chain : number of repeats". Each chain must be separated by a line break. Error found in line "%s"' % line.strip())

    return st_dict


def get_interacting_chains_alt(pathslist):  # Using file names
    """Constructs a dictionary with each chain as an entry and the chains it is interacting with as values"""
    int_dict = {}

    for path in pathslist:
        path = os.path.basename(path)       # Get filename
        path = re.sub(".pdb*", "", path)    # Remove extension
        splitted = path.split("_")
        int_dict.setdefault(splitted[1], []).append(splitted[2])    # First chain key, second chain entry
        int_dict.setdefault(splitted[2], []).append(splitted[1])    # Second chain key, first chain entry

    return int_dict



if __name__ == '__main__':
    #structures = [process_input(filename) for filename in list_files]
    structures = [process_input(list_files)]
    #print(structures[0].get_atoms().__eq__)
#    print(structure[0].get_chains)
    #superimpose(structures[0], structures[1])





"""
structure
    model
        chain
            residue
                atoms

        A  B  C

A_B
B_C
"""


    #Comprobar formato pdb
    #Qué tipos de moléculas hay
    #Comprobar todos los pares de interacciones (si sobra tiempo hacerlas de novo)

#Superimposición

    #Coger un dímero
    #Superimponer otro dímero
    #...

#Modeller (alternativa puede que más fácil)

#Evaluación del macrocomplejo

#GUI
