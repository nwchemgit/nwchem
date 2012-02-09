'''
Created on Feb 5, 2012

@author: marat
'''
from atom_params import *
from generic_atom import *
from generic_residue import *

class MySystem(object):
    '''
    classdocs
    '''


    def __init__(self, name=None,atoms=None):
        '''
        Constructor
        '''
        self.name  = "system" if name==None else name
        self.atoms = [] if atoms==None else atoms
        self.residues = {}
        self.natoms=len(self.atoms)

    def info1(self):
        output = self.name +"\n"
        for x in self.atoms:
            output = output + str(x)+"\n"
        return output
    
    @classmethod        
    def fromPDBfile(cls,filename):
        '''
        alternative constructor from PDB file
        '''
        cls = MySystem()
#        print cls
#        print cls.atoms
        fp = open(str(filename),'r')
        
        for line in fp.readlines():
            if line.startswith('ATOM'):
                a=GenericAtom.fromPDBrecord(line)
                cls.AddAtom(a)
        return cls
    
    def toPDBfile(self):
        for residue in self.residues.itervalues():
            print residue
            
    def AddAtom(self,a1):
        self.atoms.append(a1)
        rmap = self.residues
        tag = a1.groupTag()
#        print tag
        if tag not in rmap:
#            print "tag not found"
            rmap[tag]=GenericResidue()
        rmap[tag].AddAtom(a1)

        
    def info(self):
        for residue in self.residues.itervalues():
            print residue           
    
if __name__ == '__main__':
#    sim0 = MySystem("test")
#    sim1 = MySystem.fromPDBfile("test.pdb")  
#    sys.exit()

    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

#    try:
#        b = ResAtom.fromPDBrecord(aline1)
#    except SystemExit:
#        print "error reading PDB line"
#        sys.exit(1)
#    sim0.AddAtom(b)
#    a = ResAtom.fromPDBrecord(aline2)
#    sim0.AddAtom(a)
#    
    sim1 = MySystem.fromPDBfile("test.pdb")
    sim1.info()
#    sim1 = MySystem.fromPDBfile("test.pdb")
# 
#    sim1.toPDBfile()
