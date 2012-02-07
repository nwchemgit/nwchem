'''
Created on Feb 5, 2012

@author: marat
'''
from myatom import Atom

from atom_dictionary import *

class ResAtom(Atom):
    '''
    Extends Atom class to include extra information about parent residue (name and id).
    Responsibilities
    1. Create Atom object from PDB string 
    '''


    def __init__(self,params):
        '''
        Constructor uses dictionary to accept arguments
        '''
        if params.has_key("coords"):
            coords = params["coords"]
        else:
            coords=[0.0,0.0,0.0]
            
        if params.has_key("name"):
            name = params["name"]
        else:
            name = "atom"      
              
        if params.has_key("resname"):
            self.resname = params["resname"]
        else:
            self.resname = "UNK"                    
            
        Atom.__init__(self,name,coords)               
        
        if params.has_key("resid"):
            self.resid = params["resid"]
        else:
            self.resid = 1              

        if params.has_key("element"):
            self.element = params["element"]
        else:
            if name != "atom":
                self.element = AtomDictionary.elementName(name)
            else:
                self.element = "atom"    
                                       
        
    def __str__(self):
        return self.toPDBrecord()

    @classmethod        
    def fromPDBrecord(cls,string):
        '''
        alternative constructor from PDB record
        '''
        if string.startswith('ATOM'):
            d = {}
            aname = string[12:16] #this actually means characters from 12 - 15 inclusive
            d['name']=aname
            rname = string[17:20]
            d['resname']=rname
            resid = int(string[22:26].strip())
            d['resid']=resid
            x = float(string[30:38].strip())
            y = float(string[38:46].strip())
            z = float(string[46:54].strip())
            d['coords']=[x,y,z]
            return cls(d)
        else:
            sys.exit(1)

    def toPDBrecord(self,id_atom=1,id_res=0):
 
        a1=" "
        a2=2*" "
        a3=3*" "
        a4=4*" "
        a5=5*" "
        a6=6*" "
        a8=8*" "
        a22=22*" "
 
        pdbformat="%-6s%5s%1s%4s%1s%3s%1s%4s%4s%8.3f%8.3f%8.3f%22s%2s"

        if id_res==0:
            resi = str(self.resid)
        else:
            resi = id_res
        return pdbformat%('ATOM',str(id_atom),a1,self.name,a1,
                                           self.resname,a2,str(resi),a4,
                                           self.coords[0],self.coords[1],self.coords[2],a22,
                                           self.element)  
 
                                
if __name__ == '__main__':
    
    print "creating first atom"
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    print aline1
    a = ResAtom.fromPDBrecord(aline1)   
    print "it should come out as this"
    print a
    aline2="ATOM    589 2HG  GLU    18     -12.634  -3.023  -3.475  1.00  0.00           H"

    print "creating second atom"
    aline2="ATOM    589 2HG  GLU    18     -12.634  -3.023  -3.475  1.00  0.00           H"
    print aline2
    b = ResAtom.fromPDBrecord(aline2)   
    print "it should come out as this"
    print b

    print "The distance between these two atoms is", Atom.bondlength(a, b)
    
    print "PDB record for second atom with starting index 5 and resid 23" 
    print b.toPDBrecord(id_atom=5,id_res=23)

    
#    PDB ATOM RECORD FORMAT
#    COLUMNS      DATA TYPE        FIELD      DEFINITION
#    ------------------------------------------------------
#     1 -  6      Record name      "ATOM    "
#     7 - 11      Integer          serial     Atom serial number.
#    13 - 16      Atom             name       Atom name.
#    17           Character        altLoc     Alternate location indicator.
#    18 - 20      Residue name     resName    Residue name.
#    22           Character        chainID    Chain identifier.
#    23 - 26      Integer          resSeq     Residue sequence number.
#    27           AChar            iCode      Code for insertion of residues.
#    31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
#                                             Angstroms
#    39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
#                                             Angstroms
#    47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
#                                             Angstroms
#    55 - 60      Real(6.2)        occupancy  Occupancy.
#    61 - 66      Real(6.2)        tempFactor Temperature factor.
#    77 - 78      LString(2)       element    Element symbol, right-justified.
#    79 - 80      LString(2)       charge     Charge on the atom.

#    Example
#             1         2         3         4         5         6         7         8
#    12345678901234567890123456789012345678901234567890123456789012345678901234567890
#    MODEL        1
#    ATOM      1  N   ALA     1      11.104   6.134  -6.504  1.00  0.00           N
#    ATOM      2  CA  ALA     1      11.639   6.071  -5.147  1.00  0.00           C
#    ...
#    ...
#    ATOM    293 1HG  GLU    18     -14.861  -4.847   0.361  1.00  0.00           H
#    ATOM    294 2HG  GLU    18     -13.518  -3.769   0.084  1.00  0.00           H
                                          

