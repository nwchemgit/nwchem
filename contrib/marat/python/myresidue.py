'''
Created on Feb 5, 2012

@author: marat
'''
from myatom import *
class MyResidue:
    '''
    classdocs
    '''


    def __init__(self,name="UNK",id=0,atoms=[]):
        '''
        Constructor
        '''
        self.name = name
        self.id = id
        self.atoms = atoms
        
    def __str__(self):
        output = self.name +"\n"
        for x in self.atoms:
            output = output + str(x)+"\n"
        return output
            
    def AddAtom(self,a1):
        self.atoms.append(a1)
        print self.atoms
            
if __name__ == '__main__':
    res0 = MyResidue("test")
    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

    try:
        b = Atom.fromPDBrecord(aline1)
    except SystemExit:
        print "error reading PDB line"
        sys.exit(1)
    res0.AddAtom(b)
#    a = Atom.fromPDBrecord(aline2)
#    sim0.AddAtom(a)
    print res0