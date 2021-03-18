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
