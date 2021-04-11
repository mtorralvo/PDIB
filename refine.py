"""This module makes use of PyRosseta to refine the obtained models."""

from pyrosetta.toolbox.cleaning import cleanATOM
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.scoring import get_score_function
import pyrosetta
import os

def refine(model_path, logging, relax_type = "f"):
    """Refines the model using PyRosetta."""

    logging.info(f"Refining the model at {model_path}, this may take several time.")
    
    pyrosetta.init()
    
    basepath = os.path.splitext(model_path)[0]

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

    rosetta_object.dump_pdb(f"{basepath}_refined.pdb")

    logging.info(f"Refined model saved at {basepath}_refined.pdb")
