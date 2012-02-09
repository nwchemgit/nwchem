'''
Created on Feb 7, 2012

@author: marat
'''
    
class PDBAtomRecord(object):
    '''
    Class PDBAtomRecord processes ATOM or HETATM line in PDB file.
    For all other records None value will be returned.
    PDBAtomRecord.test is used to verify PDB line is of ATOM or HETATM type.
    PDBAtomRecord.dct will return all available records in a dictionary with the keys 
    indicated below
    '''
#   Dictionary for column locations (see at the end for more info)
    irec=dict(record=[0,6],atomid=[6,11],name=[12,16],resname=[17,20],
              resid=[22,26],coord=[30,54],element=[76,78])
#   Dictionary for types of fields
    atype=dict(record="string",atomid="int",name="string",resname="string",
               resid="int",coord="float array",element="string")

#    def __init__(self):
#        '''
#        Constructor
#        '''
#        pass

    
    @staticmethod
    def field(name,buf):
        '''
        returns value of the "name" field in the provided "buf" buffer
        always returns None value if buffer is not of ATOM or HETATM type
        '''
#        print inspect.getsource(PDBAtomRecord.fieldName)   
        if not PDBAtomRecord.test(buf):
            return None     
        ir=PDBAtomRecord.irec[name]
        atype=PDBAtomRecord.atype[name]
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
    def dct(buf):
        '''
        returns all the available fields in a dictionary
        '''
        if not PDBAtomRecord.test(buf):
            return None
        
        d={}
        for name in PDBAtomRecord.irec.iterkeys():
            value=PDBAtomRecord.field(name,buf)
            if value!=None:
                d[name]=value

        return d   
   
    @staticmethod     
    def test(buf):
        '''
        tests whether buffer is of ATOM or HETATM type
        '''
        return buf[0:6] in ("ATOM  ","HETATM")
    

if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    aline2="ATOM    588      GLU           -13.363  -4.163  -2.372  1.00  0.00"
    aline3="ATTM    588      GLU           -13.363  -4.163  -2.372  1.00  0.00"
    
    print PDBAtomRecord.field("name",'') 
    print PDBAtomRecord.field("name",aline3) 
    print PDBAtomRecord.dct(aline1)
    print PDBAtomRecord.dct(aline2)
    print PDBAtomRecord.dct(aline3)
    print PDBAtomRecord.dct('')


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
                                          

        