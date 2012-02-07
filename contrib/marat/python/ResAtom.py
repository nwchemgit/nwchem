'''
Created on Feb 5, 2012

@author: marat
'''
from myatom import *
import string 
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
            element = params["element"]
        else:
            if name != "atom":
                self.element = AtomDictionary.elementName(name)
            else:
                self.element = "atom"    
                                       
        
    def __str__(self):
        output = self.name + "  " 
        output = output + self.resname + "  " 
        output = output + str(self.resid) + "  " 
        for x in self.coords:
            output = output + "  %12.6f"%x
        return output + "  " + self.element

    @classmethod        
    def fromPDBrecord(cls,string):
        '''
        alternative constructor from PDB record
        '''
        if string.startswith('ATOM'):
            d = {}
            aname = string[12:16]
            print "here",string[12:16]
            print "here",aname
            d['name']=aname
            rname = string[17:20].strip()
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

    def toPDBrecord(self):
 
        a1=" "
        a2=2*" "
        a3=3*" "
        a4=4*" "
        a5=5*" "
        a6=6*" "
        a8=8*" "
        a22=22*" "
        
        rec = 78*" "
        recname=a6
        serial=a5
        atom=a4
        resname=a3
        chain=a1 
        seq=a5
        x=a8
        y=a8
        z=a8

        recname = "ATOM"
        pdbformat="%-6s%5s%1s%4s%1s%3s%1s%4s%4s%8.3f%8.3f%8.3f%22s%2s"

        print pdbformat%('ATOM',str(1),a1,self.name,a1,
                                           self.resname,a2,str(self.resid),a4,
                                           self.coords[0],self.coords[1],self.coords[2],a22,
                                           self.element)  
                
#        print '"{:<30}"'.format('left aligned')
#        print '{:<4}'.format('ATOM')
#        i=1
#        serial="%4d"%i
        
        print recname+serial
#        if string.startswith('ATOM'):
#            d = {}
#            aname = string[12:16].strip()
#            d['name']=aname
#            rname = string[17:20].strip()
#            d['resname']=rname
#            resid = int(string[22:26].strip())
#            d['resid']=resid
#            x = float(string[30:38].strip())
#            y = float(string[38:46].strip())
#            z = float(string[46:54].strip())
#            d['coords']=[x,y,z]
#            return cls(d)
#        else:
#            sys.exit(1)
                                
if __name__ == '__main__':
    aline1="ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H"
    aline2="ATOM    589 2HG  GLU    18     -12.634  -3.023  -3.475  1.00  0.00           H"
#    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
#    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

    a = ResAtom.fromPDBrecord(aline2)
    b = ResAtom.fromPDBrecord(aline1)
    print a
    print b
    print Atom.bondlength(a, b)
    print aline1[0:78]
    b.toPDBrecord()
    a.toPDBrecord()
    