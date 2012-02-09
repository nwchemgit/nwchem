'''
Created on Feb 7, 2012

@author: marat
'''
import sys
from pdbrecord import PDBAtomRecord 
from numpy import array,linalg
from atom_params import *

class GenericAtom(object):
    '''
    classdocs
    '''

    def __init__(self,d):
        '''
        Constructor
        '''
        self.dct = d

        if d:
            if type(d) is not type({}):
                print "wrong type ", type(d)
                print "expecting", type({})
                sys.exit(1)
            else:
                self.dct  = d
        try:
            self.coord=array(self.dct.pop("coord"))
        except:
            self.coord=None
                
    @classmethod                
    def fromPDBrecord(cls,buf):
        '''
        alternative constructor from PDB record
        '''
        d=PDBAtomRecord.dct(buf)
        if d:
            return cls(d)
        return None
    
    @staticmethod
    def bondlength(a1,a2):
        dr=a1.coord-a2.coord
        return linalg.norm(dr)
    
    def covRadius(self):
        name = self.dct["name"]
        return AtomParams.covRadius(name)   
   
    def groupTag(self):
        '''
        returns that identifies group association of an atom
        based on residue name and residue id (e.g. "ASP_1")
        Default values of "UNK" and "0" are used should 
        residue name and residue id be absent
        '''
        resname = self.dct.get("resname","UNK").strip()
        resid = str(self.dct.get("resid",0)).strip()
        return "_".join((resname,resid))
        
    def __str__(self):
        return str(self.dct) + " " + str(self.coord)

        
    @staticmethod    
    def bonded(a1,a2):
        dr = linalg.norm(a1.coord-a2.coord)
        return dr <= (a1.covRadius()+a2.covRadius())
            
if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    aline2="ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39           I"
    aline3="ATOM      2  O1  IO3     1      -0.985  -1.156   1.840       -0.80           O"
    aline4="ATOM      2  O1                 -0.985  -1.156   1.840       -0.80           O"

    a=GenericAtom.fromPDBrecord(aline2)   
    print a.groupTag()
    print a.coord
    print a.dct     
    b=GenericAtom.fromPDBrecord(aline3)
    print b.coord
    print b.dct         
    print GenericAtom.bondlength(a,b),GenericAtom.bonded(a,b)
    print a.covRadius()+b.covRadius()
    c=GenericAtom.fromPDBrecord(aline1)
    print GenericAtom.bondlength(a,c),GenericAtom.bonded(a,c)
    c=GenericAtom.fromPDBrecord(aline4)    
    print c.groupTag()