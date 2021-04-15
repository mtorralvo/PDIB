from distutils.core import setup

setup(name='PYT-SBI',
    version='1.0',
    description='PYT-SBI project',
    author='Mireia Codina Tobías, Ángel Monsalve Fernández, María Torralvo Márquez',
    url='https://github.com/angelmf97/PYT-project',
    packages=['PDIB'],
    scripts=['PDIB/PDIB.py', 'PDIB/refine.py', 'PDIB/compute_energy.py', 'PDIB/process_input_files.py', 'PDIB/main_functions.py', 'PDIB/gui.py']
    )
