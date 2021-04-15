#!/usr/bin/env python
# coding=utf-8

'''Main program. Handles the user input and calls the main functions'''

import argparse
import re
import os
import logging
import sys
import glob
import random
import string
import numpy as np
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

parser = argparse.ArgumentParser(description='Make protein origami and stuff.')

parser.add_argument( '-i', '--input-directory',
    dest = 'input_directory',
    action = 'store',
    required = True,                    
    default = './',
    help = 'Input directory containing PDB files')

parser.add_argument('-o', '--output-directory',
    dest = 'output_directory',
    action = 'store',
    default = './',
    required = True,
    help = 'Output directory to store the models')

parser.add_argument('-f', '--force',
    dest = 'force',
    action = 'store_true',
    default = False,
    help = 'Overwrites the output directory if already exists')

parser.add_argument('-v', '--verbose',
    dest = 'verbose',
    action = 'store_true',
    default = False,
    help = 'Gives information during the execution of the script')

parser.add_argument('-all', '--all_models',
    dest='all_models',
    action='store_true',
    default=False,
    help='Try all the models. It might take a long time! It is not recommended (specially when the complex is large) unless you suspect that more than one model is possible')

parser.add_argument('-sto', '--stoichiometry',
   dest='stoich',
   action='store',
   default=None,
   help='File with sequence ID and stoichiometry')

parser.add_argument('-n_models', '--number-of-models',
   dest='num_models',
   type=int,
   action='store',
   default=1,
   help='Number of models to be generated (default 1)')

# Default 60 because 62 is the maximum number of chains that can be stored in a single pdb file (non multicharacter chain ids)
parser.add_argument('-n_chains', '--number-of-chains',
   dest='num_chains',
   type=int,
   action='store',
   default=60,
   help='Expected number of chains')

parser.add_argument('-e', '--energies',
   dest='energies',
   action='store_true',
   help='Compute the energies of the resulting models and their z-scores')

parser.add_argument('-r', '--refine',
   dest='refine',
   action='store_true',
   help='Use PyRosetta in order to improve the obtained models')

args = parser.parse_args() # List of arguments

from process_input_files import *
from main_functions import *
from compute_energy import *

if args.verbose:
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.DEBUG,
        format = '%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')

if not args.stoich or args.num_chains:
    warning = input('No stoichiometry file nor number of chains specified. The program will keep running until 60 chains are added. Do you want to continue (y/n): ')
    
    if warning in ['y', 'Y', 'yes', 'Yes']:
        pass
    
    else:
        logging.warning('Exiting the program.')
        quit()

input_dir = args.input_directory # Read input directory

regex = re.compile(".*\.pdb(.gz)*")

list_files = []
# Check if input is a directory and if True, save files in a list
if os.path.isdir(input_dir):
    list_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if regex.match(f)]

    if len(os.listdir(input_dir)) == 0:
        raise ValueError('Input folder is empty')

    elif len(list_files) < 2:
        raise ValueError('No files found in the specified directory. Make sure they are stored as [filename].pdb')

elif os.path.isfile(input_dir):
    raise ValueError('Input must be a directory')

logging.info('%d input files found.\n' %len(list_files))

# Generate the output folders in case they do not exist and the force option is activated
if not os.path.isdir(args.output_directory) or args.force:
    os.makedirs('./%s/structures' % args.output_directory, exist_ok = True)
    os.makedirs('./%s/analysis' % args.output_directory, exist_ok = True)

else:
    raise OSError('%s already exists, specify a different output directory or enable option -f to override the already existing one.' % args.output_directory)

info_files, info_seqs = process_input(list_files)   # Generate dictionaries with info about the sequencues and the files they are in

# Generate stoichiometry dictionary
stoichiometry = None
if args.stoich:
    stoichiometry = get_stoichiometry(args.stoich)  

parser = PDBParser(PERMISSIVE=1)

first_file = random.choice(list(info_files.keys()))
first_structure = parser.get_structure('structure', first_file)

# Change the chain id names
for chain in first_structure.get_chains():
    original_id = chain.id  # Original ID
    chain.id = [x for x in info_files[first_file] if x[0] == original_id][0]    # New ID

final_models = []

final_models = create_model(first_structure, final_models, info_files, args.num_chains, args.num_models, args.all_models, stoich=stoichiometry, verbose=False)


if len(final_models) == 0:
        sys.stderr.write('No complex could be obtained with the provided files. \n')

chain_alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
            'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'
            , 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4'
            , '5', '6', '7', '8', '9']

model_paths = []

logging.info('Storing the final models.')

for final_model in final_models:

    # PDB format limits the number of possible chains that can be added to a file, as the chain ID is assigned
    # only one position. In order to solve that, the program assigns a character from chain_alphabet to each chain
    # in the model. If the model has more than 62 chains (the length of the alphabet), t is splitted in several files.
    structure = pdb_struct.Structure('id')  # empty structure:

    chain_counter = 0
    model_counter = 0
    legend = ''

    chains = list(final_model.get_chains()) # Get all the chains in the final model

    for chain in chains:

        # Assign a new name from the alphabet
        new_id = chain_alphabet[chain_counter]

        # Create a legend so the user can map the old chain names to the new ones
        legend += 'CHAIN    %s   %s   %s\n' % (new_id, chain.id[0], info_seqs[chain.id[1]])
      
        chain.id = new_id   # Change the ID
        chain_counter += 1  # One more chain

        # Create new model in the first chain
        if chain_counter == 1:
            model_counter += 1
            model = pdb_model.Model(model_counter)

        model.add(chain)

        # If the chain number exceeds the number of characters in the alphabet, reset the counter and
        # add the structure
        if chain_counter == len(chain_alphabet):
            
            chain_counter = 0

            model = model.copy()
            structure.add(model)

            # As IDs can't be duplicated in a model, delete the chain that have already been processed
            for chain_m in model.get_chains():
                final_model[0].detach_child(chain_m.id)

    # Repeat for the final chain
    current_models_strcuture = [x.id for x in structure.get_models()]
    if model.id not in current_models_strcuture:
        structure.add(model)

    # Model name: model_n.pdb where n is increased for each model
    written_id = 0
    PDB_name = 'model_' + str(written_id) + '.pdb'
    while PDB_name in os.listdir(args.output_directory + '/structures/') :
        written_id += 1
        PDB_name = 'model_'+str(written_id)+'.pdb'

    logging.info(f'Model {PDB_name} stored at {args.output_directory}')

    # Specify the output path and append it to a list of paths
    model_path = args.output_directory + '/structures/' + PDB_name
    model_paths.append(model_path)

    # Save the coordinates
    io = PDBIO()
    io.set_structure(structure)
    io.save(model_path)

    # Save the legend
    fh = open(model_path, 'a')
    fh.write('CHAIN HEADER    current id   molecule id    sequence\n')
    fh.write(legend)
    fh.write('\n')
    fh.close()

# Remove the temporary directory if exists
import shutil
shutil.rmtree('PDIB_tmp', ignore_errors=True)

# Run relaxation protocol with PyRosetta
if args.refine:
    try:
        from refine import *
    
    except ModuleNotFoundError:
        logging.warning('Module PyRossetta not installed. No refined models will be generated.')
    
    else:
        
        os.makedirs('./%s/refined' % args.output_directory, exist_ok = True)
        
        for path in model_paths:
            refine(path, args.output_directory)
            model_paths.append(path)

# Compute z-scores with ProSa2003
if args.energies:
    compute_energy(model_paths, args.output_directory)

