import argparse
import re
import os
import logging
import sys

parser = argparse.ArgumentParser(description='Make protein origami and stuff.')

parser.add_argument('-i', '--input',
                    dest = 'input',
                    type = str,
                    required = True)

parser.add_argument('-s', '--stechiometry',
                    type = str)

parser.add_argument('-o', '--output-directory',
                    dest = 'output_directory',
                    type = str,
                    required = True)

parser.add_argument('-f', '--force',
                    dest = 'force',
                    action = 'store_true',
                    default = False)

parser.add_argument('-v', '--verbose',
                    dest = 'verbose',
                    action = 'store_true',
                    default = False)

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(stream=sys.stderr,
                        level=logging.DEBUG,
                        format = '%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')

regex = re.compile("\w+_\w+_\w+.pdb(.gz)*")
if not regex.match(args.input):
    raise ValueError("Input file must be of the format <name>_<chain1>_<chain2>.pdb(.gz)")

if not os.path.isdir(args.output_directory) or args.force:
    os.makedirs("./%s/structures" % args.output_directory, exist_ok = True)
    os.makedirs("./%s/analysis" % args.output_directory, exist_ok = True)

else:
    raise OSError("%s already exists, specify a different output directory or enable option -f to override the already existing." % args.output_directory)

logging.info('doing something')
