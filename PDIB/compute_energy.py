import os
import tempfile
import subprocess
import logging

def compute_energy(model_paths, output_path):
    """Compute the combined energies and z-scores of the generated models in order to compare them"""

    # Create the commands to analyse energies and produce z-scores in a temporary file
    commands = tempfile.NamedTemporaryFile(dir='.', delete=False)

    for path in model_paths:
        commands.write(f"read pdb {path} {path}\n".encode())
    
    commands.write(b"init zscore\n")
    commands.write(f"zscore * {output_path}/analysis/zscores\n".encode())
    commands.write(b"exit")

    commands.seek(0)
  
    try:
        process = subprocess.Popen(["prosa2003", commands.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    except FileNotFoundError:
        logging.error("ProSa command not found. Check that you have ProSa installed and the alias for calling the program is 'prosa2003'")
        
    else:
        process.communicate()
    
    finally:
        commands.close()
        os.remove(commands.name)

    return

    if __name__ == '__main__':
    
        import sys
        from sys import argv
        
        logging.basicConfig(
            stream=sys.stderr,
            level=logging.DEBUG,
            format = '%(asctime)s %(message)s',
            datefmt='%m/%d/%Y %H:%M:%S')


        compute_energy(argv[1], logging)