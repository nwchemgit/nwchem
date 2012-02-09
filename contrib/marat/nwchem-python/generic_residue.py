'''
Created on Feb 5, 2012

@author: marat
'''
import sys
from GenericAtom import *

class GenericResidue(object):
    '''
    classdocs
    '''


    def __init__(self,atoms=None):
        '''
        Default constructor for Residue class
        atoms list of atoms in the residue
        name residue name
        '''
        if atoms:
            self.atoms = atoms
        else:
            self.atoms=[]
            

    def AddAtom(self,a):
        self.atoms.append(a)
        
if __name__ == '__main__':
    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

    res0 = GenericResidue()
    print res0
    a = GenericAtom.fromPDBrecord(aline2)
    print a
    res0.AddAtom(a)
    print res0.atoms[0]
    
#    b = ResAtom.fromPDBrecord(aline1)
#    res0.AddAtom(a)
#    res0.AddAtom(b)
#    print res0.toPDBrecord(atom_start=1)