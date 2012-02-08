'''
Created on Feb 7, 2012

@author: marat
'''
import string 
from array import *
class PDBRecord:
    '''
    classdocs
    '''
    irecord=[0,6]
    inum=[6,11]
    iname=[12,16]
    iresname=[17,20]
    iresid=[22,26]
    ix=[30,38]
    iy=[38,46]
    iz=[46,54]
    iel=[76,78]
    

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

    def __init__(self,line):
        '''
        Constructor
        '''
        self.buf=line


    def recordName(self):
        return PDBRecord.genrecord(self.buf,self.irecord)

    def atomName(self):
        return PDBRecord.genrecord(self.buf,self.iname)

    def resName(self):
        return PDBRecord.genrecord(self.buf,self.iresname)

    def elemName(self):
        return PDBRecord.genrecord(self.buf,self.iel)
    
    def resId(self):
        return int(PDBRecord.genrecord(self.buf,self.iresid))
    
    def atomCoord(self):
        return [float(PDBRecord.genrecord(self.buf,self.ix)),
                float(PDBRecord.genrecord(self.buf,self.iy)),
                float(PDBRecord.genrecord(self.buf,self.iz))]

              
                            
    @staticmethod     
    def genrecord(buf,ir):
        return buf[ir[0]:ir[1]] 

        
    def __str__(self):
        return self.buf
    

if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    pdb1 = PDBRecord(aline1)
    print aline1
    print pdb1.recordName()
    print pdb1.atomName()    
    print pdb1.resName() 
    print pdb1.resId()
    print pdb1.elemName()     
    print pdb1.atomCoord()          

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
                                          

        