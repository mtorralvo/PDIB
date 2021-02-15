import argparse
import re
import os

parser = argparse.ArgumentParser(description='Make protein origami and stuff.')

parser.add_argument('-i', '--input',
                    dest = 'input',
                    type = str,
                    required = True)

parser.add_argument('-s', '--stechiometry',
                    type = str,
                    required = False)

parser.add_argument('-o', '--output-directory',
                    dest = 'output_directory',
                    type = str,
                    required = True)

parser.add_argument('-f', '--force',
                    dest = 'force')

args = parser.parse_args()

regex = re.compile("\w+_\w+_\w+.pdb(.gz)*")
if not regex.match(args.input):
    raise IOError("Input file must have the format <..........>")

os.makedirs("./%s/structures" % args.output_directory, exist_ok = True)
os.makedirs("./%s/analysis" % args.output_directory, exist_ok = True)
