from Bio.PDB import NeighborSearch
import subprocess

def check_clashes(model, new_chain, clash_perc = 0.3):
    
    """Checks if there is any clash between the backbone of the existing model and the chain to be added."""
    
    model_atoms = [atom for atom in model.get_atoms() if atom.id in ("CA", "C1")]
    new_chain_atoms = (atom for atom in new_chain.get_atoms() if atom.id in ("CA", "C1"))
    
    neighbor_object = NeighborSearch(model_atoms)
    
    clashes = 0
    length = 0
    for atom in new_chain_atoms:
        length +=1 
        if neighbor_object.search(atom.coord, 2):
            clashes += 1
    
    proportion_of_clashes = float(clashes/length)
    
    if proportion_of_clashes < clash_perc:
        return False
    
    else:
        return True
    

import tempfile
import os
def compute_energy(model_paths, logging):

    print(model_paths)
    
    # Create the commands to analyse energies and produce z-scores in a temporary file
    commands = tempfile.NamedTemporaryFile()

    for path in model_paths:
        commands.write(f"read pdb '{path}' {path}\n".encode())
    
    commands.write(b"init zscore\n")
    commands.write(b"zscore *")

    process = subprocess.Popen(["prosa2003", commands.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    commands.close()

    stdout, stderr = process.communicate()

    logging.info(stdout)
    logging.error(stderr)

    return







