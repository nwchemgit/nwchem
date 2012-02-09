'''
Created on Feb 7, 2012

@author: marat
'''
import sys
from pdbrecord import PDBAtomRecord 
from numpy import array,linalg

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

            
if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    aline2="ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39           I"
    aline3="ATOM      2  O1  IO3     1      -0.985  -1.156   1.840       -0.80           O"

    a=GenericAtom.fromPDBrecord(aline2)   
    print a.coord
    print a.dct     
    b=GenericAtom.fromPDBrecord(aline3)
    print b.coord
    print b.dct         
    print GenericAtom.bondlength(a,b)