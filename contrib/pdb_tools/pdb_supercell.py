#!/bin/env python
#
# Build a super cell of a given PDB file. 
# Because of peculiarities of NWChem we need to produce a
# super cell by creating more atoms with translated coordinates.
# Then we need to write all the non-water atoms out first 
# followed by all water atoms.
#
def read_pdb_file(filename):
    '''
    Read the contents of the PDB file and return a list of lines.
    '''
    fp = open(filename,'r')
    pdb = fp.readlines()
    fp.close()
    return pdb

def get_cell_dimensions(pdb):
    '''
    Take the PDB file as a list of lines. 
    Read the cell dimensions from the first line (starting with "CRYST1")
    and return this as a tuple.
    '''
    line = pdb[0]
    elements = line.split()
    x = float(elements[1])
    y = float(elements[2])
    z = float(elements[3])
    return (x,y,z)

def gen_translation(direction,pdb):
    '''
    Generate the translation vector based on the choice of direction (given as a string),
    and the cell dimension of the PDB file.
    Return this a tuple.
    '''
    (x,y,z) = get_cell_dimensions(pdb)
    if direction == "x":
        trans = (x,0.0,0.0)
    elif direction == "y":
        trans = (0.0,y,0.0)
    elif direction == "z":
        trans = (0.0,0.0,z)
    else:
        print("Interesting coordinate system you are using. Your choice of direction is %s\n"%direction)
    return trans

def new_cell_dimensions(direction,pdb):
    '''
    Given the PDB data and the direction compute the new cell dimensions
    after extending the system in the given dimension. 
    Return the new dimensions as a tuple.
    '''
    (x,y,z) = get_cell_dimensions(pdb)
    if direction == "x":
        cell = (2*x,y,z)
    elif direction == "y":
        cell = (x,2*y,z)
    elif direction == "z":
        cell = (x,y,2*z)
    else:
        print("Interesting coordinate system you are using. Your choice of direction is %s\n"%direction)
    return cell

def split_pdb(pdb):
    '''
    Given the PDB file, separate the contents into three parts:
    - the cell definition
    - the solute atoms
    - the solvent atoms
    Return these results as a tuple of lists of lines.
    '''
    cryst   = []
    solute  = []
    solvent = []
    for line in pdb:
        elements = line.split()
        if elements[0] == "CRYST1":
            cryst.append(line)
        elif elements[0] == "ATOM":
            if elements[3] == "HOH":
                solvent.append(line)
            else:
                solute.append(line)
        elif elements[0] == "TER":
            solute.append(line)
    return (cryst,solute,solvent)

def translate_atoms(atomsin,translation):
    '''
    Take a given list of atoms and create a new list of atoms
    where all atoms are translated by the given translations.
    Return the new list of atoms. The input list remains unchanged.
    '''
    (dx,dy,dz) = translation
    atomsout = []
    for atom in atomsin:
        elements = atom.split()
        if elements[0] == "ATOM":
            sx = atom[30:38]
            sy = atom[38:46]
            sz = atom[46:54]
            x  = float(sx)+dx
            y  = float(sy)+dy
            z  = float(sz)+dz
            new_atom = "%s%8.3f%8.3f%8.3f%s"%(atom[:30],x,y,z,atom[54:])
        else:
            new_atom = atom
        atomsout.append(new_atom)
    return atomsout

def new_pdb(cryst,oldsolute,oldsolvent,newsolute,newsolvent):
    '''
    Given the different sections of the old PDB and the translated parts of the 
    solute and solvent lists construct a new PDB file. 
    Return the list of lines of the new PDB file.
    '''
    end = ["END"]
    newpdb = cryst+oldsolute+newsolute+oldsolvent+newsolvent+end
    return newpdb

def new_cell(oldcryst,newcell):
    '''
    Given the old crystal line and the new cell dimensions construct 
    a new crystal line. 
    Return the list of new crystal lines.
    '''
    inline = oldcryst[0]
    (x,y,z) = newcell
    outline = "%s%9.3f%9.3f%9.3f%s"%(inline[:6],x,y,z,inline[33:])
    outlist = []
    outlist.append(outline)
    return outlist

def write_pdb(filename,pdb):
    '''
    Given the name of the output file and the list of PDB lines 
    write the data to the PDB file.
    '''
    fp = open(filename,'w')
    for line in pdb:
        fp.write(line)
    fp.close()

def parse_arguments():
    '''
    Parse command line arguments.
    '''
    from argparse import ArgumentParser
    prs = ArgumentParser()
    prs.add_argument("input",help="the input PDF file")
    prs.add_argument("output",help="the output PDF file")
    prs.add_argument("direction",help="the direction to expand along")
    args = prs.parse_args()
    return args

def execute_with_arguments(args):
    inputfile  = args.input
    outputfile = args.output
    direction  = args.direction
    inputdata = read_pdb_file(inputfile)
    trans = gen_translation(direction,inputdata)
    newcell = new_cell_dimensions(direction,inputdata)
    (cryst,solute,solvent) = split_pdb(inputdata)
    newcryst = new_cell(cryst,newcell)
    newsolute = translate_atoms(solute,trans)
    newsolvent = translate_atoms(solvent,trans)
    outputdata = new_pdb(newcryst,solute,solvent,newsolute,newsolvent)
    write_pdb(outputfile,outputdata)

def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
