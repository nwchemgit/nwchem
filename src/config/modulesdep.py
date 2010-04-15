#!/usr/bin/env python
'''Look for the MODULES := $(pathsubst ... section of the make_nwchem_config
file and convert it into a DOT file so we can visualize the dependencies.'''

import os

# Create the dot input file from the make_nwchem_config file.
dotfilename = 'modulesdep.dot'
dotfile = open(dotfilename, 'w')
dotfile.write('digraph G {\n')
start_parse = False
for line in open('make_nwchem_config'):
    line = line.strip()
    if not line:
        continue
    elif 'will be hard to generate automatically' in line:
        start_parse = True
        continue
    elif line.startswith('MODULES :=') and start_parse:
        parts = line.split(',')
        if len(parts) == 3:
            modules = parts[1].strip().split()
            for module in modules[1:]:
                dotfile.write('%s -> %s\n' % (modules[0], module))
        else:
            start_parse = False
            continue
dotfile.write('}')
dotfile.close()

# Run dot.
os.system('dot -Tpng %s -omodulesdep.png' % dotfilename)
