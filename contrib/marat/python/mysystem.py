'''
Created on Feb 5, 2012

@author: marat
'''
from ResAtom import *
from Residue import *
from atom_dictionary import *
class MySystem:
    '''
    classdocs
    '''


    def __init__(self, name="system",atoms=[]):
        '''
        Constructor
        '''
        self.name  = name
        self.atoms = atoms
        self.residue = {}
        self.residues = {}
    def __str__(self):
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
        fp = open(str(filename),'r')
        
        for line in fp.readlines():
            if line.startswith('ATOM'):
                a=ResAtom.fromPDBrecord(line)
                cls.AddAtom(a)
        return cls
       
    def AddAtom(self,a1):
        self.atoms.append(a1)
        rmap = self.residues
        if not rmap.has_key(a1.resid):
            rmap[a1.resid]=Residue()
        rmap[a1.resid].AddAtom(a1)
                        
    
if __name__ == '__main__':
    sim0 = MySystem("test")
    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

    try:
        b = ResAtom.fromPDBrecord(aline1)
    except SystemExit:
        print "error reading PDB line"
        sys.exit(1)
    sim0.AddAtom(b)
    a = ResAtom.fromPDBrecord(aline2)
    sim0.AddAtom(a)
    sim1 = MySystem.fromPDBfile("test.pdb")
    print sim1.residues

    
#    def AddAtom1(self,a1):
#        self.atoms.append(a1)
#        rmap = self.residue
#        if not rmap.has_key(a1.resid):
#            rmap[a1.resid]={}
#            rmap[a1.resid]["atoms"]=[]
#            rmap[a1.resid]["name"]=a1.resname
#        if rmap[a1.resid]["name"]!=a1.resname:
#            print "different names for the same residue index"
#            sys.exit(1)
#        rmap[a1.resid]["atoms"].append(a1)
#        print "added atom", a1.resid    