import argparse
import re
import os
import logging
import sys

parser = argparse.ArgumentParser(description='Make protein origami and stuff.')

parser.add_argument('-i', '--input',
                    dest = 'input',
                    type = str,
                    required = True,
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

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(stream=sys.stderr,
                        level=logging.DEBUG,
                        format = '%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')

regex = re.compile(".*\w+_\w+_\w+.pdb(.gz)*")

print(args.input)
for filein in args.input:
    if not regex.match(filein):
        raise ValueError("Input file must be of the format <name>_<chain1>_<chain2>.pdb(.gz)")

if not os.path.isdir(args.output_directory) or args.force:
    os.makedirs("./%s/structures" % args.output_directory, exist_ok = True)
    os.makedirs("./%s/analysis" % args.output_directory, exist_ok = True)

else:
    raise OSError("%s already exists, specify a different output directory or enable option -f to override the already existing." % args.output_directory)

logging.info('doing something')


#Procesado del Input
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
def process_input(filename):

    parser = PDBParser(PERMISSIVE = True, QUIET = True)

    with open(filename) as fh:
        structure = parser.get_structure("Name", fh)


    #print([atom.coord for atom in structure[0].get_atoms()])
    n_chains = len(structure[0].child_list)

    if n_chains != 2:
        logging.error("There should be 2 chains in your pdb file. %d of chains found in %s instead." % (n_chains, filename))

    return structure

def superimpose(structure1, structure2):

    superimposer = Superimposer()
    atoms1 = list(structure1.get_atoms())
    atoms2 = list(structure2.get_atoms())
    if len(atoms1) > len(atoms2):
        atoms1 = atoms1[:len(atoms2)]

    elif len(atoms2) > len(atoms1):
        atoms2 = atoms2[:len(atoms1)]

    superimposer.set_atoms(atoms1, atoms2)
    print(dir(superimposer))
    #print(superimpose(atoms1, atoms2))




if __name__ == '__main__':
    structures = [process_input(filename) for filename in args.input]
    print(structures[0].get_atoms().__eq__)
    superimpose(structures[0], structures[1])





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
