### GUI ###
from tkinter import *
from tkinter import Checkbutton
from tkinter import Entry
import os
import re
import sys
from Bio.PDB import *
import tkinter.filedialog
import tkinter.messagebox
import collections
import numpy as np
from process_input_files import *
from main_functions import *
from compute_energy import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
#from refine import *


class Application(tkinter.Frame):

    # Get the input files from the input directory and insert them in the list box of the PDB files
    def input_dir(self):

        self.folder_path = tkinter.filedialog.askdirectory(title="Select directory with PDB files")
        if self.folder_path:
            for file in os.listdir(self.folder_path):
                name = re.search('.*\.pdb', file)
                if name:
                    file_name = self.folder_path + "/" + name.group()
                    self.pdb_files_paths.append(file_name)
                    self.list_files.append(name.group())

        self.seq_listbox.delete(0, END)
        index = 1
        for record in self.list_files:
            self.seq_listbox.insert(index, record)
            index += 1

    # Select the output directory where the models and the z-scores are going to be saved by the program
    def output_dir(self):
        self.output_dir = tkinter.filedialog.askdirectory(title="Select output directory")

    # The required inputs are the input directory with PDB files, the output directory where we store the models and the number of models to generate
    def create_input_options(self):

        self.options_frame = tkinter.LabelFrame(self, text="Required Inputs", width=20)

        input_options_frame = tkinter.Frame(self.options_frame, width=390, borderwidth=2)

        L1 = Label(input_options_frame, text="Select input folder")
        L1.pack(side=LEFT)
        input_dir_entry = Button(input_options_frame, text="Browse", command=self.input_dir)
        input_dir_entry.pack(side=RIGHT)
        input_options_frame.pack(fill=BOTH)

        output_options_frame = tkinter.Frame(self.options_frame, width=390, borderwidth=2)

        L2 = Label(output_options_frame, text="Select output folder")
        L2.pack(side=LEFT)
        output_dir_entry = Button(output_options_frame, text="Browse", command=self.output_dir)
        output_dir_entry.pack(side=RIGHT)
        output_options_frame.pack(fill=BOTH)

        choose_frame = Frame(self.options_frame, width=390, borderwidth=2)

        L3 = Label(choose_frame, text="Choose between these options:")
        L3.pack(side=LEFT)
        choose_frame.pack(fill=BOTH)

        number_models_frame = Frame(self.options_frame, width=390, borderwidth=2)

        self.number_models = IntVar(value=2)
        L5 = Label(number_models_frame, text="Enter number of models")
        L5.pack(side=LEFT)
        number_models_entry = Entry(number_models_frame, width=5, textvariable=self.number_models)
        number_models_entry.pack(side=RIGHT)

        number_models_frame.pack(fill=BOTH)

        all_models_frame = Frame(self.options_frame, width=390, borderwidth=2)

        self.all_models = BooleanVar(value=False)

        exhaustive = Checkbutton(all_models_frame, text="Create all possible models", variable=self.all_models, onvalue=True, offvalue=False)
        exhaustive.pack(side=LEFT)

        all_models_frame.pack(fill=BOTH)

        self.options_frame.grid(row=0, column=0)

    def create_left_frame(self):
        self.left_frame = tkinter.LabelFrame(self, text="PDB Files List", padx=5, pady=5, background="green")
        self.create_seq_listbox()
        self.left_frame.grid(row=1, column=0)


    def create_seq_listbox(self):

        frame = tkinter.Frame(self.left_frame)

        scrollbar = tkinter.Scrollbar(frame, orient=tkinter.VERTICAL)

        self.seq_listbox = tkinter.Listbox(frame, selectmode=tkinter.SINGLE,
						   height=20,
						   yscrollcommand = scrollbar.set)

        scrollbar.config(command=self.seq_listbox.yview)

        scrollbar.pack( side=tkinter.RIGHT, fill=tkinter.Y)

        self.seq_listbox.pack( side=tkinter.LEFT, expand=True, fill=tkinter.BOTH)

        frame.pack( fill=tkinter.BOTH )

    # As optional inputs we have the refinement of the model, obtaining the z-scores or/and specifying the stoichiometry
    def create_optional_inputs (self):
        self.optional_inputs_frame = tkinter.LabelFrame(self, text="Optional Inputs", width=20)

        optional_inputs_frame = tkinter.Frame(self.optional_inputs_frame, width=390, borderwidth=2)

        self.refine_model = BooleanVar(value=False)
        self.get_z_scores = BooleanVar(value=False)

        refine = Checkbutton(optional_inputs_frame, text="Refine models", variable=self.refine_model, onvalue=True, offvalue=False)
        refine.pack(side=LEFT)

        get_z_scores = Checkbutton(optional_inputs_frame, text="Z-scores", variable=self.get_z_scores, onvalue=True, offvalue=False)
        get_z_scores.pack(side=LEFT)

        optional_inputs_frame.pack(fill=BOTH)

        stoich_frame = tkinter.Frame(self.optional_inputs_frame, width=390, borderwidth=2)

        L4 = Label(stoich_frame, text="Select stoichiometry file")
        L4.pack(side=LEFT)
        stoich_file_entry = Button(stoich_frame, text="Browse", command=self.stoichiometry_file)
        stoich_file_entry.pack(side=RIGHT)
        stoich_frame.pack(fill=BOTH)

        number_chains_frame = Frame(self.optional_inputs_frame, width=390, borderwidth=2)

        self.number_chains = IntVar(value=60)
        L6 = Label(number_chains_frame, text="Enter number of maximum chains")
        L6.pack(side=LEFT)
        number_chains_entry = Entry(number_chains_frame, width=5, textvariable=self.number_chains)
        number_chains_entry.pack(side=RIGHT)

        number_chains_frame.pack(fill=BOTH)

        self.optional_inputs_frame.grid(row=0, column=1)


    # Obtain stoichiometry from txt file
    def stoichiometry_file(self):
        self.file_path = tkinter.filedialog.askopenfilename(title="Select a txt file with stoichiometry", filetypes = (("text files","txt"),))


    def create_run(self):
        self.create_run = tkinter.LabelFrame(self)

        run_frame = tkinter.Frame(self.create_run)

        run = Button(run_frame, text="Run PDIB", command=self.build_complex)
        run.pack(side=RIGHT)

        run_frame.pack(fill=BOTH)

        self.create_run.grid(row=1, column=1)

    def build_complex(self):

        n_models = self.number_models.get()
        all_models = self.all_models.get()
        n_chains = self.number_chains.get()
        refine = self.refine_model.get()
        z_scores = self.get_z_scores.get()
        stoich = self.stoichiometry_file
        out_dir = self.output_dir

        #if os.listdir(out_dir):
        #    os.makedirs('./%s/structures' % out_dir, exist_ok = True)
        #    os.makedirs('./%s/analysis' % out_dir, exist_ok = True)

        # The files from the input directory are stored in a list
        list_files = self.pdb_files_paths

        info_files, info_seqs = process_input(list_files)

        # In case there is a stoichiometry file
        stoichiometry = None
        if stoich:
            stoichiometry = get_stoichiometry(stoich)

        parser = PDBParser(PERMISSIVE=1)

        first_file = random.choice(list(info_files.keys()))
        first_structure = parser.get_structure("structure", first_file)

        # change the chain id names
        for chain in first_structure.get_chains():
            original_id = chain.id  # Original ID
            chain.id = [x for x in info_files[first_file] if x[0] == original_id][0]

        final_models = []

        #final_models = create_model(first_structure, final_models, info_files, num_chains, n_models, stoich=stoichiometry, verbose=False)
        exhaustive = False
        final_models = create_model(first_structure, final_models, info_files, n_chains, n_models, exhaustive, stoich=stoichiometry, verbose=False)

        if len(final_models) == 0:
                sys.stderr.write("No complex could be obtained with the provided files. \n")

        chain_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                    "T", "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"
                    , "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "0", "1", "2", "3", "4"
                    , "5", "6", "7", "8", "9"]

        model_paths = []

        logging.info('Storing the final models.')

        for final_model in final_models:

            # PDB format limits the number of possible chains that can be added to a file, as the chain ID is assigned
            # only one position. In order to solve that, the program assigns a character from chain_alphabet to each chain
            # in the model. If the model has more than 62 chains (the length of the alphabet), t is splitted in several files.
            structure = pdb_struct.Structure('id')  # empty structure:


            chain_counter = 0
            model_counter = 0
            legend = ""

            chains = list(final_model.get_chains()) # Get all the chains in the final model

            for chain in chains:

                # Assign a new name from the alphabet
                new_id = chain_alphabet[chain_counter]

                # Create a legend so the user can map the old chain names to the new ones
                legend += 'CHAIN    %s   %s   %s\n' % (new_id, chain.id[1], info_seqs[chain.id[1]])

                chain.id = new_id # Change the ID
                chain_counter += 1 # One more chain

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

        # # Model name: model_n.pdb where n is increased for each model
            written_id = 0
            PDB_name = 'model_' + str(written_id) + '.pdb'
            while PDB_name in os.listdir(out_dir + '/structures/') :
                written_id += 1
                PDB_name = 'model_'+str(written_id)+'.pdb'

            logging.info(f'Model {PDB_name} stored at {out_dir}')

            # Specify the output path and append it to a list of paths
            model_path = out_dir + '/structures/' + PDB_name
            model_paths.append(model_path)

            # Save the coordinates
            io = PDBIO()
            io.set_structure(structure)
            io.save(model_path)

            # Save the legend
            fd = open(model_path, 'a')
            fd.write('CHAIN HEADER    current id   molecule id    sequence\n')
            fd.write(legend)
            fd.write('\n')
            fd.close()

        # Remove the temporary directory if exists
        import shutil
        shutil.rmtree("PDIB_tmp", ignore_errors=True)

        #if refine:
        #    try:
                #from refine import *

            #except ModuleNotFoundError:
            #    logging.warning("Module PyRossetta not installed. No refined models will be generated.")

            #else:
            #    os.makedirs("./%s/refined" % out_dir, exist_ok = True)
            #    for path in model_paths:

            #        refine(path, out_dir)
            #        model_paths.append(path)

        # Compute the z-scores with ProSa2003
        if z_scores:
            compute_energy(model_paths, out_dir)


    def createWidgets(self):

        self.create_left_frame()
        self.create_input_options()
        self.create_optional_inputs()
        self.create_run()

        self.grid(row=0)

    def __init__(self, master=None, **kwargs):

        tkinter.Frame.__init__(self, master, **kwargs)

        self.master.wm_title("Protein-DNA Interaction Builder")
        self.master.resizable(width=False, height=False)

        self.config( width=600 )
        self.config( height=600 )

        self.left_frame = None
        self.seq_listbox = None
        self.input_options = None
        self.optional_inputs = None
        self.run = None
        self.pdb_files_paths = []
        self.list_files = []
        self.number_models = None
        self.all_models = None
        self.number_chains = None
        self.refine_model = None
        self.get_z_scores = None
        self.stoichiometry_file = None

        self.createWidgets()

root = tkinter.Tk()
app = Application(master=root,padx=10, pady=10)

app.mainloop()
