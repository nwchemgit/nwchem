'''
Created on Feb 4, 2012

@author: marat
'''
import sys
from numpy import *

import string

class Atom(object):
    '''
    simple class to handle atom specifications
    '''

    def __init__(self,name="X",coords=[0.0,0.0,0.0]):
        '''
        Constructor
        '''
        self.name = name
        self.coords = array(coords)

    def __str__(self):
#        output = "%12.6f %12.6f %12.6f"%(self.coords[0],self.coords[1],self.coords[2])
        output = self.name 
        for x in self.coords:
            output = output + "  %12.6f"%x
        return output
    
    def toXYZrecord(self):
        output = self.name 
        for x in self.coords:
            output = output + "  %12.6f"%x 
        return output+"\n"
    
    @classmethod
    def fromXYZrecord(cls,aline):
        AtomStr  = string.split(aline)[0:4]
        name = AtomStr[0]
        x = string.atof(AtomStr[1])
        y = string.atof(AtomStr[2])
        z = string.atof(AtomStr[3])
        return cls(name,[x,y,z])

    @classmethod        
    def fromPDBrecord(cls,string):
        if string.startswith('ATOM'):
            name = string[12:16].strip()
            x = float(string[30:38].strip())
            y = float(string[38:46].strip())
            z = float(string[46:54].strip())
            return cls(name,[x,y,z])
        else:
            sys.exit(1)
    @staticmethod
    def bondlength(a1,a2):
        dr=a1.coords-a2.coords
        return linalg.norm(dr)
    
                
if __name__ == '__main__':
    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

    try:
        b = Atom.fromPDBrecord(aline1)
    except SystemExit:
        print "error reading PDB line"
        sys.exit(1)
    try:
        a = Atom.fromPDBrecord(aline2)
    except SystemExit:
        print "error reading PDB line"
        sys.exit(1)        
    print "bondlength is", Atom.bondlength(a,b)
    print "finished I am"