"""This module makes use of PyRosseta to refine the obtained models."""
from pyrosetta.toolbox.cleaning import cleanATOM
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.scoring import get_score_function
import pyrosetta
import os
import logging

def refine(model_path, output_path, relax_type = "f"):
    """Refines the model using PyRosetta."""

    logging.info(f"Refining the model at {model_path}, this may take several time.")
    
    pyrosetta.init()
    
    basepath = os.path.splitext(model_path)[0]
    basename = os.path.basename(model_path)

    cleanATOM(model_path)                 # Keep only ATOM and TER
    new_path = basepath + ".clean.pdb"    # Path of the clean pdb

    rosetta_object = pose_from_pdb(new_path)    # Creates a pose object used by pyRosetta

    score_function = get_score_function()       # Creates a scoring function

    if relax_type == "c":
        relax_model = ClassicRelax()            # Classic relax (time expensive)
    
    if relax_type == "f":
        relax_model = FastRelax()               # Fast relax
    
    relax_model.set_scorefxn(score_function)
    relax_model.apply(rosetta_object)

    rosetta_object.dump_pdb(f"{output_path}/refined/{basename}_refined.pdb")

    os.remove(new_path)

    logging.info(f"Refined model saved at {output_path}/refined/{basename}_refined.pdb")

if __name__ == '__main__':
    
    import sys
    from sys import argv
    
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.DEBUG,
        format = '%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')

    refine(argv[1], logging)