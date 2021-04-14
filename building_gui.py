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

        number_models_frame = Frame(self.options_frame, width=390, borderwidth=2)

        self.number_models = IntVar(value=3)

        L3 = Label(number_models_frame, text="Enter number of models")
        L3.pack(side=LEFT)
        number_models_entry = Entry(number_models_frame, width=5, textvariable=self.number_models)
        number_models_entry.pack(side=RIGHT)

        number_models_frame.pack(fill=BOTH)

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
        stoich_file_entry = Button(stoich_frame, text="Browse", command=self.get_stoichiometry)
        stoich_file_entry.pack(side=RIGHT)
        stoich_frame.pack(fill=BOTH)

        self.optional_inputs_frame.grid(row=0, column=1)


    # Obtain stoichiometry from txt file
    def get_stoichiometry(self):
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
        refine = self.refine_model.get()
        z_scores = self.get_z_scores.get()
        stoich = self.get_stoichiometry.get()
        out_dir = self.output_dir

        list_files = self.pdb_files_paths


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

        self.createWidgets()

root = tkinter.Tk()
app = Application(master=root,padx=10, pady=10)

app.mainloop()
