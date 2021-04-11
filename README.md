# PYT-project

Proposed steps: 
1. Obtain the structure and chains for the PDB files found in the directory
2. Store information about unique chains found in the directory: check if it is stored (sequence alignment) or assign a unique identifier
3. Pairwise interactions with a common chain matched in step 2 with the same ID: superimpose them by fixing one and moving the other, obtaining
a rotation matrix that is used to copy and move the other chain (the one with unmatched ID) from the moving structure to the fixed and add it to the model.
4. Recursive approach: 
  - First level of recursion: a pairwise interaction is chosen and all the pairs containing one of these two initial molecules are tried to be added in the complex. 
  - Second level of recursion: addition of more subunits 
5. RMSD calculation: computed between the core atoms of each of the residues. Search bibliography to set a cutoff and classify a structural superimposition of two molecules as the same or not. 
6. Detecting clashes: each time a new chain is added to the complex, the model NeighbourSearch is used to detect atoms located at closer distance than a specific threshold (define threshold based on bibliography).
7. Producing several complexes from same input structures
8. Refine the model
9. GUI


## Model scoring
By making use of ProSa2003, the program allows the user to compute the z-scores of the generated models. These scores are calculated by summing the energies (both surface, paired and combined) for all the residues of a model and comparing it to the energies of several models with different folds.

This generates a file with all the models sorted by their z-scores. A low z-score usually indicates that the model is good and stable, while higher or even positive z-scores account for incorrect model with high pseudoenergies.

## Clash detections
Of all the possible chains that can be added to the developing model, not all of them result in reliable models when added. It may happen that, when adding a new chain to the model, the superimposition is succesful but the chain presents clashes with other parts of the structure.

We implemented a function for detecting clashes between chains and discarding the chains that conflict with the structure. By default, we considered a clash whenever more than a 3% of the CA and CB atoms of the added chain where closer than 2Å from any CA or CB of the rest of the structure, although these thresholds can be modified.

## Model refinement
The program makes use of PyRosetta to refine the models, allowing for some flexibility. PyRosetta is an interactive Python-based interface to the Rosetta molecular modeling suite. Our program allows to choose between different relaxing methods, such as the classic relax protocol or the fast relax protocol. This relaxing methods allow for all-atom refinement searching the local conformational space around the starting structure. This refinement allows the program to yield better models with lower energies.


