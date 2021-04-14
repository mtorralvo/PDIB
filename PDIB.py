from main_functions import *
from refine import *
from process_input_files import *
from Bio.PDB import PDBIO


def generate_output(structure, name):
    io = PDBIO()
    io.set_structure(structure)
    io.save("structures/" + name, write_end=True)

if __name__ == '__main__':

  
    # Generate dictionaries containing alternative IDs and the sequences of each chain
    info_files, info_seqs = process_input(list_files)   

    print(info_files)
    quit()

    model_seed = random.choice(list(info_files.keys()))  # First interaction is chosen randomly

    parser = PDBParser()
    seed_structure = parser.get_structure('seed', model_seed)

    stoich = get_stoichiometry(args.stoichiometry)

    new_id_stoich = {}
    for chain in seed_structure.get_chains():
        if chain.id in stoich.keys():
            curr_id = chain.id  # the so-called chain accession, usually a letter (A,B,C,D ...)
            new_id = [x for x in info_files[model_seed] if x[0] == curr_id][0]
            chain.id = new_id
            new_id_stoich[chain.id[1]] = stoich[curr_id]

    structures = create_model(seed_structure, [], info_files, 5, 2, args, exhaustive=False, verbose=True, stoich=new_id_stoich)

    with open("output.txt", "w+") as output:
        for atom in structures[0].get_atoms():
            output.write(str(atom.coord))

    for structure in structures:
        generate_output(structure)

