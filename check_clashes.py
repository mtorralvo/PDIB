from Bio.PDB import NeighborSearch

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
    