import argparse
import re
import os
import logging
import sys
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(description='Make protein origami and stuff.')

parser.add_argument('-i', '--input',
                    dest = 'input',
                    type = str,
                    required = False,
                    nargs='+')

parser.add_argument('-s', '--stoichiometry',
                    type = str)

parser.add_argument('-o', '--output-directory',
                    dest = 'output_directory',
                    type = str,
                    required = True)

parser.add_argument('-f', '--force',
                    dest = 'force',
                    action = 'store_true',
                    default = False)

parser.add_argument('-v', '--verbose',
                    dest = 'verbose',
                    action = 'store_true',
                    default = False)

parser.add_argument('-d', '--directory',
                    dest = 'directory',
                    default = False)

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(stream=sys.stderr,
                        level=logging.DEBUG,
                        format = '%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')

regex = re.compile(".*\w+_\w+_\w+.pdb(.gz)*")

if not args.directory:
    for filein in args.input:
        if not regex.match(filein):
            raise ValueError("Input file must be of the format <name>_<chain1>_<chain2>.pdb(.gz)")

if not os.path.isdir(args.output_directory) or args.force:
    os.makedirs("./%s/structures" % args.output_directory, exist_ok = True)
    os.makedirs("./%s/analysis" % args.output_directory, exist_ok = True)

else:
    raise OSError("%s already exists, specify a different output directory or enable option -f to override the already existing." % args.output_directory)

logging.info('doing something')


def process_input(filename):
    """Processes the input converting a pdb file to a model object of Bio.PDB"""

    parser = PDBParser(PERMISSIVE = True, QUIET = True)

    with open(filename) as fh:
        structure = parser.get_structure("Name", fh)

    n_chains = len(structure[0].child_list)

    if n_chains != 2:
        logging.error("There should be 2 chains in your pdb file. %d of chains found in %s." % (n_chains, filename))

    return structure[0]


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