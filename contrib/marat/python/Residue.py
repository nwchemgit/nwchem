'''
Created on Feb 5, 2012

@author: marat
'''
import sys
from ResAtom import *

class Residue:
    '''
    classdocs
    '''


    def __init__(self,params={}):
        '''
        Default constructor for Residue class
        atoms list of atoms in the residue
        name residue name
        '''
        if params.has_key("atoms"):
            self.atoms = params["coords"]
        else:
            self.atoms=[]
            
        if params.has_key("name"):
            self.name = params["name"]
        else:
            self.name = ""
        
    def __str__(self):
        output = ""
        for a in self.atoms:
            output=output + str(a)+"\n"
        return output
    
    def toPDBrecord(self,id_atom=1,id_res=1):
        output = ""
        i=id_atom-1
        for a in self.atoms:
            i = i + 1
            output=output + a.toPDBrecord(i) +"\n"
        return output    


    def AddAtom(self,a):
        if self.name =="" :
            self.name = a.resname
        else:
            if a.resname != self.name:
                print "different names for the same residue index"
                sys.exit(1)                             
        self.atoms.append(a)
        
if __name__ == '__main__':
    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

    res0 = Residue()
    a = ResAtom.fromPDBrecord(aline2)
    b = ResAtom.fromPDBrecord(aline1)
    res0.AddAtom(a)
    res0.AddAtom(b)
    print res0.toPDBrecord(id_atom=4)