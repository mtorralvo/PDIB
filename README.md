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

## Installation

The program can be installed from the PDIB/ directory using the following command:

```
python setup.py install
```

## Command line arguments

```
usage: project2.py [-h] [-i INPUT [INPUT ...]] [-s STOICHIOMETRY] -o
                   OUTPUT_DIRECTORY [-f] [-v] [-d DIRECTORY]

Make protein complexes from protein-protein and protein-DNA interactions.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input files in PDB format. Use this option to use
                        individual files as input. To work with an entire
                        directory use the option -d.
  -s STOICHIOMETRY, --stoichiometry STOICHIOMETRY
                        Specifies the stoichiometry of the complex. Each chain
                        must be separated by newlines. Each line must contain
                        the number of the chainfollowed by a colon (:) or an
                        equal sign (=) and the number of copies of that chain
                        in the final complex.
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Specifies the output directory.
  -f, --force           Force the overwrite of the output directory if it
                        already exists.
  -v, --verbose         Enables detailed information about the program run.
  -d DIRECTORY, --directory DIRECTORY
                        Specifies the directory where all the input files are
                        to be taken. To work with individual files use the
                        option -i.
```

## Usage example

An example is provided with the program and can be run by using the following command:
```
pdib -d example -o . -f -v > log.txt
```
This will create two folders inside the current directory: "Structures" and "Analysis". The first one contains all the structures that PDIB could find meeting the requirements of the user. The latter contains ...........

As the verbose (-v) option is enabled, the program will print information about the process to the standard output, that is redirected to a new file named "log.txt" in this example.


## Program overview

Protein-DNA Interaction Builder (PDIB) is a python program that aims to generate protein macrocomplexes from known interactions between pairs of chains. The program is able to model both protein-protein, protein-DNA and protein-RNA interactions.

As an input, it takes PDB files containing pairs of interacting chains that will form the final structure, in order to know the relative position of each of the subunits. By making use of a recursive approach, PDIB explores different ways of arranging the specified chains, yielding several models that preserve the interactions between pairs of chains as specified in the interaction files.

The program also """""allows""""" the user to use an stechiometry file, specifying the number of times that each chain should appear in the final complex, and the program will try to find solutions that meet the user requirements.

In order to yield better and more reliable results, the program can use PyRosetta to perform a relax pipeline that optimizes the energy of the models. In addition, ProSa2003 can be runned from PDIB to retrieve the z-scores of the generated model, so the user can identify those with a better energy profile.

The program can be run from the command line or using a graphical user interface (GUI) that allows to run PDIB with all its functionalities in a more user-friendly way.


## Input files

PDIB will look for all the files having a .pdb extension in the directory specified by -d. The filenames must match the following pattern:
<Name>_<Chain1>_<Chain2>.pdb(.gz)

Where name is the name of the structure and Chain1 and Chain2 are the interacting subunits present in the file.
"""""""""""Optionally""""""""", the user may provide an stoichiometry file, indicating the number of times each chain appears in the final complex. The stoichiometry file must have the following format:

```
Chain1: 1
Chain2: 1
```
Where the number at the right indicates the number of appearances of that chain. The colons (":") can be replaced with an equal sign ("=") if the user pleases. The names appearing in the stoichiometry file must match those indicated in the .pdb filenames.


### Input preprocessing

First, the program looks for all the PDB files contained in the specified directory, keeping only those that match the pattern previously mentioned (<Name>_<Chain1>_<Chain2>.pdb(.gz)). Chains are randomly assigned a new unique name in order to avoid conflicts during the runtime. Then, a dictionary is generated containing each chain as a key and a list of its possible interactions as a value. """In case""" the user provides a stoichiometry file, it is converted to a dictionary as well.

## Main algorithm

PDIB uses a recursive approach to explore the different possible structures that satisfy the conditions provided by the user. The model starts with an interacting pair of chains, and recursively explores all the possible chains that could be added to the structure, attending to the available interactions specified by the .pdb filenames.

For each possible interaction, the pairwise interaction file containing both itself and the new chain that is to be added. The first is superimposed with its homolog in the structure, and the resulting transformation matrix is applied to the added chain so it is aligned with the structure in the manner specified by the pairwise interaction. The program then checks that the new chain has no clashes with the backbone of the rest of the structure, discarding the new chain if it is the case.

## Clash detections

Of all the possible chains that can be added to the developing model, not all of them result in reliable models when added. It may happen that, when adding a new chain to the model, the superimposition is succesful but the chain presents clashes with other parts of the structure.

We implemented a function for detecting clashes between chains and discarding the chains that conflict with the structure. By default, we considered a clash whenever more than a 3% of the CA and CB atoms of the added chain where closer than 2Å from any CA or CB of the rest of the structure, although these thresholds can be modified.

## Model scoring
By making use of ProSa2003, the program allows the user to compute the z-scores of the generated models. These scores are calculated by summing the energies (both surface, paired and combined) for all the residues of a model and comparing it to the energies of several models with different folds.

This generates a file with all the models sorted by their z-scores. A low z-score usually indicates that the model is good and stable, while higher or even positive z-scores account for incorrect model with high pseudoenergies.


## Model refinement
The program makes use of PyRosetta to refine the models, allowing for some flexibility. PyRosetta is an interactive Python-based interface to the Rosetta molecular modeling suite. Our program allows to choose between different relaxing methods, such as the classic relax protocol or the fast relax protocol. This relaxing methods allow for all-atom refinement searching the local conformational space around the starting structure. This refinement allows the program to yield better models with lower energies.


