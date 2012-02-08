'''
Created on Feb 7, 2012

@author: marat
'''
import string 
from array import *
import inspect

def rrecordName(buf):
#        print inspect.getsource(PDBAtomParser.recordName)
    return buf
    
class PDBAtomParser:
    '''
    class PDBAtomParser parses ATOM or HETATM line in PDB file.
    It can only properly handle ATOM or HETATM records.
    Before getting individual fields one should test
    whether the record is atom based using PDBAtomParser.isAtomType() method.
    Perhaps the best way to use this class at this point is to generate 
    dictionary using PDBAtomParser.getDict(aline1)
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

    def __init__(self):
        '''
        Constructor
        '''
        pass

    @staticmethod     
    def genrecord(buf,ir):
        '''
        returns value of a record as a string
        if string is empty exception is raised
        '''
        value = buf[ir[0]:ir[1]]
        if value.strip()=='':
            value=None
            raise ValueError("empty field while parsing PDB file")
        return value 
    
    @staticmethod   
    def recordName(buf):
#        print inspect.getsource(PDBAtomParser.recordName)
        try:
            return PDBAtomParser.genrecord(buf,PDBAtomParser.irecord)
        except:
            raise ValueError

    @staticmethod     
    def atomId(buf):
        value= PDBAtomParser.genrecord(buf,PDBAtomParser.inum)
        try:
            return int(value)
        except:
            raise ValueError("bad field"+value)
                
    @staticmethod   
    def atomName(buf):
        try:
            return PDBAtomParser.genrecord(buf,PDBAtomParser.iname)
        except:
            raise ValueError
            
    @staticmethod   
    def resName(buf):
        try:
            return PDBAtomParser.genrecord(buf,PDBAtomParser.iresname)
        except:
            raise ValueError  
          
    @staticmethod   
    def elemName(buf):
        try:
            return PDBAtomParser.genrecord(buf,PDBAtomParser.iel)
        except:
            raise ValueError
            
    @staticmethod     
    def resId(buf):
        value = PDBAtomParser.genrecord(buf,PDBAtomParser.iresid)
        try:
            return int(value)
        except:
            raise ValueError("bad field"+str(value))
        
    @staticmethod        
    def toInt(f):
        value = f
        try:
            return int(value)
        except:
            raise ValueError("bad field"+str(value))
                
    @staticmethod
    def resId0(buf):
        return int(PDBAtomParser.genrecord(buf,PDBAtomParser.iresid))
 

            
    @staticmethod     
    def atomCoord(buf):
        value = (PDBAtomParser.genrecord(buf,PDBAtomParser.ix),
                 PDBAtomParser.genrecord(buf,PDBAtomParser.iy),
                 PDBAtomParser.genrecord(buf,PDBAtomParser.iz))
        try:
            return [float(x) for x in value]
        except:
            raise ValueError("bad field"+str(value))
        
#        try:
#            return [float(PDBAtomParser.genrecord(buf,PDBAtomParser.ix)),
#                    float(PDBAtomParser.genrecord(buf,PDBAtomParser.iy)),
#                    float(PDBAtomParser.genrecord(buf,PDBAtomParser.iz))]
#        except:
#            raise ValueError
        
    @staticmethod              
    def getDict(buf):
        d = {}
        d["record"]=PDBAtomParser.recordName(buf)
#        if d["record"] not in ("ATOM  ", "HETATM"):
#            raise RuntimeError('this is the error message')
#            print d["record"]
#            return d  
        try:
            value=PDBAtomParser.elemName(buf)  
        except:
            pass
        else:
            d["element"]=value

        try:
            value=PDBAtomParser.atomCoord(buf)  
        except:
            pass
        else:
            d["coords"]=value            
            
        try:
            value=PDBAtomParser.atomName(buf) 
        except:
            pass
        else:
            d["name"]=value                    

        try:
            value=PDBAtomParser.resId(buf) 
        except:
            pass
        else:
            d["resid"]=value                    

        try:
            value=PDBAtomParser.atomId(buf)
        except:
            pass
        else:
            d["atomid"]=value          
            
        try:
            value=PDBAtomParser.resName(buf)
        except:
            pass
        else:
            d["resname"]=value      
        

        return d   
                          

    
    @staticmethod     
    def isAtomType(buf):
        return buf[0:4]!="ATOM" or buf[0:6]!="HETATM"
              
        
    def __str__(self):
        return self.buf
    

if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    A8     -13.363  -4.163  -2.372  1.00  0.00           H"
    aline2="ATOM    5A8      GLU    18     -13.363  -4.163  -2.372  1.00  0.00"

#    print PDBAtomParser.atomId(aline1) 
    
#    print aline1
#    print PDBAtomParser.recordName(aline1)
#    print PDBAtomParser.atomName(aline2)    
#    print PDBAtomParser.resName(aline1) 
#    print PDBAtomParser.resId0(aline1)
#    print PDBAtomParser.elemName(aline2)     
#    print PDBAtomParser.atomCoord(aline1)        
#    print PDBAtomParser.isAtomType(aline1)   
     
#    print PDBAtomParser.getDict("567 rty")
#    print PDBAtomParser.getDict(aline2)
    print PDBAtomParser.getDict(aline1)


#
#    try:
#        print PDBAtomParser.getDict("567 rty")
#    except:
#        print "Cannpot do that"


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
                                          

        