'''
Created on Feb 7, 2012

@author: marat
'''
    
class PDBAtomParser(object):
    '''
    class PDBAtomParser parses ATOM or HETATM line in PDB file.
    It can only properly handle ATOM or HETATM records.
    Before getting individual fields one should test
    whether the record is atom based using PDBAtomParser.isAtomType() method.
    Perhaps the best way to use this class at this point is to generate 
    dictionary using PDBAtomParser.getDict(aline1)
    '''

    irec=dict(record=[0,6],atomid=[6,11],name=[12,16],resname=[17,20],
              resid=[22,26],coord=[30,54],element=[76,78])
    atype=dict(record="string",atomid="int",name="string",resname="string",
               resid="int",coord="float array",element="string")

    def __init__(self):
        '''
        Constructor
        '''
        pass

    
    @staticmethod
    def record(name,buf):
#        print inspect.getsource(PDBAtomParser.recordName)        
        ir=PDBAtomParser.irec[name]
        atype=PDBAtomParser.atype[name]
        value = buf[ir[0]:ir[1]]
        if value.strip()=='':
            return None
        if atype=='int':
            try:
                value = int(value)
            except:
                value = None
        elif atype=='float':
            try:
                value = float(value)
            except:
                value = None
        elif atype=='string':
            pass
        elif atype=='float array':
            try:
                value = [float(x) for x  in value.split()]
            except:
                value = None
        else:
            raise ValueError("unknown type definition")
                                
        return value 
    
    @staticmethod
    def getDict(buf):
#        d = dict((name,PDBAtomParser.record(name,buf)) 
#                 for name in PDBAtomParser.irec.iterkeys() )
        d={}
        for name in PDBAtomParser.irec.iterkeys():
            value=PDBAtomParser.record(name,buf)
            if value!=None:
                d[name]=value

        return d   
   
    @staticmethod     
    def isAtomType(buf):
        return buf[0:4]!="ATOM" or buf[0:6]!="HETATM"
              
        
    def __str__(self):
        return self.buf
    

if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    aline2="ATOM    588      GLU           -13.363  -4.163  -2.372  1.00  0.00"
    aline3="ATTM    588      GLU           -13.363  -4.163  -2.372  1.00  0.00"
    
    print PDBAtomParser.record("name",aline2) 
    print PDBAtomParser.record("name",aline1) 
    print PDBAtomParser.getDict(aline1)
    print PDBAtomParser.getDict(aline2)
    print PDBAtomParser.getDict(aline3)


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
                                          

        