'''
Created on Feb 7, 2012

@author: marat
'''
import sys
from pdbrecord import PDBAtomRecord 
import numpy 
import math
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
            self.coord=numpy.array(self.dct.pop("coord"))
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

    @classmethod
    def fromXYZrecord(cls,aline):
        atomstr  = string.split(aline)[0:4]
        if len(atomstr) < 4:
            return None
            
        name = atomstr[0]
        coord=[float(x) for x  in atomstr[1:]]
        return cls({"grouptag":"XYZ", "coord":coord,"name":name})
        
    @staticmethod
    def bondlength(a1,a2):
        dr=a1.coord-a2.coord
        return numpy.linalg.norm(dr)

    @staticmethod
    def angle(a1,a2,a3):
        u=a1.coord-a2.coord
        v=a3.coord-a2.coord
        c = numpy.dot(u,v)/numpy.linalg.norm(u)/numpy.linalg.norm(v) 
        return  math.degrees(math.acos(c))

    def translate(self,v):
        self.coord=self.coord-v
            
    def covRadius(self):
        name = self.dct["name"]
        return AtomParams.covRadius(name)  
     
    def elemName(self):
        name = self.dct["name"]
        return AtomParams.elementName(name)   

    def name(self):
        return self.dct["name"]
    
    def res_name(self):
        '''
        '''
        return  self.dct.get("resname","UNK").strip()
                   
    def groupTag(self):
        '''
        returns that identifies group association of an atom
        based on residue name and residue id (e.g. "ASP_1")
        Default values of "UNK" and "0" are used should 
        residue name and residue id be absent
        '''
        tag = self.dct.get("grouptag")
        if tag is None:
            resname = self.dct.get("resname","UNK").strip()
            resid = str(self.dct.get("resid",0)).strip()
            tag ="_".join((resname,resid))

        return tag
 
    def setGroupTag(self,tag):
        '''
        returns that identifies group association of an atom
        based on residue name and residue id (e.g. "ASP_1")
        Default values of "UNK" and "0" are used should 
        residue name and residue id be absent
        '''
        self.dct["grouptag"]=tag
        
    def resTag(self):
        '''
        returns that identifies res association of an atom
        based on residue name and residue id (e.g. "ASP_1")
        '''

        resname = self.dct.get("resname","UNK").strip()
        resid = str(self.dct.get("resid",0)).strip()
        tag ="_".join((resname,resid))

        return tag    
    def setBond(self,i):
        '''
        returns that identifies group association of an atom
        based on residue name and residue id (e.g. "ASP_1")
        Default values of "UNK" and "0" are used should 
        residue name and residue id be absent
        '''
        bond = self.dct.setdefault("bond",set())
        bond.add(i)
    
        
    def __str__(self):
        return str(self.dct) + " " + str(self.coord)

        
    @staticmethod    
    def bonded(a1,a2):
        dr = numpy.linalg.norm(a1.coord-a2.coord)
        return dr <= (a1.covRadius()+a2.covRadius())
            
if __name__ == '__main__':
#    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
#    aline2="ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39           I"
#    aline3="ATOM      2  O1  IO3     1      -0.985  -1.156   1.840       -0.80           O"
#    aline4="ATOM      2  O1                 -0.985  -1.156   1.840       -0.80           O"
#
#    a=GenericAtom.fromPDBrecord(aline2)   
#    print a.groupTag()
#    print a.coord
#    print a.dct     
#    b=GenericAtom.fromPDBrecord(aline3)
#    print b.coord
#    print b.dct         
#    print GenericAtom.bondlength(a,b),GenericAtom.bonded(a,b)
#    print a.covRadius()+b.covRadius()
#    c=GenericAtom.fromPDBrecord(aline1)
#    print GenericAtom.bondlength(a,c),GenericAtom.bonded(a,c)
#    c=GenericAtom.fromPDBrecord(aline4)    
#    print c.groupTag()
    aline6="O1                 -0.985  -1.156   1.140  0.0 0.0  "
    aline5="O1                 -0.985  -1.156   "
    d=GenericAtom.fromXYZrecord(aline6)
    d.translate([1.0,1.0,1.0])
    print d
