'''
Created on Feb 5, 2012

@author: marat
'''
from atom_params import *
from generic_atom import *
from generic_residue import *

class MySystem(object):
    '''
    classdocs
    '''


    def __init__(self, name=None,atoms=None):
        '''
        Constructor
        '''
        self.name  = "system" if name==None else name
        self.atoms = [] if atoms==None else atoms
        self.residues = {}
        self.natoms=len(self.atoms)

    def info1(self):
        output = self.name +"\n"
        for x in self.atoms:
            output = output + str(x)+"\n"
        return output
    
    @classmethod        
    def fromPDBfile(cls,filename):
        '''
        alternative constructor from PDB file
        '''
        cls = MySystem()
#        print cls
#        print cls.atoms
        fp = open(str(filename),'r')
        
        for line in fp.readlines():
            if line.startswith('ATOM'):
                a=GenericAtom.fromPDBrecord(line)
                cls.AddAtom(a)
        return cls

    @classmethod        
    def fromXYZfile(cls,filename):
        '''
        alternative constructor from PDB file
        '''
        cls = MySystem()

        fp = open(str(filename),'r')
        
        for line in fp.readlines():
            a=GenericAtom.fromXYZrecord(line)
            if a is not None:
                cls.AddAtom(a)
        return cls
        
    def toPDBfile(self):
        for residue in self.residues.itervalues():
            print residue
            
    def AddAtom(self,a1):
        self.atoms.append(a1)
        rmap = self.residues
        tag = a1.groupTag()
        if tag not in rmap:
            rmap[tag]=GenericResidue()
        rmap[tag].AddAtom(a1)

    def connectAtoms(self):
        for i in range(len(self.atoms)):
            for j in range(i+1,len(self.atoms)):
                a1=self.atoms[i]
                a2=self.atoms[j]
                if GenericAtom.bonded(a1, a2):
                    a1.setBond(j)
                    a2.setBond(i)
                    
        for a in self.atoms:
            print a                    

    def groupAtoms(self):
        rmap = self.residues
        ir = 0
        for i in range(len(self.atoms)):
            a1=self.atoms[i]
            if a1.groupTag()!="XYZ":
                continue
            ir = ir +1
#            tag = '{:0>3}'.format(ir)
            tag="00"+str(ir)
            tag=tag[-3:]
            a1.setGroupTag(tag)
            rmap[tag]=GenericResidue()
            rmap[tag].AddAtom(a1)
            rmap["XYZ"].delAtom(a1)
            for j in range(i+1,len(self.atoms)):                
                a2=self.atoms[j]
                if GenericAtom.bonded(a1, a2):
                    a2.setGroupTag(tag)
                    rmap[tag].AddAtom(a2)
                    rmap["XYZ"].delAtom(a2)
        del rmap["XYZ"]
        for a in self.atoms:
            print a                       
                    
    def info(self):
        for tag,res in self.residues.iteritems():
            print tag 
            print res       
    
if __name__ == '__main__':
#    sim0 = MySystem("test")
#    sim1 = MySystem.fromPDBfile("test.pdb")  
#    sys.exit()

    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"

#    try:
#        b = ResAtom.fromPDBrecord(aline1)
#    except SystemExit:
#        print "error reading PDB line"
#        sys.exit(1)
#    sim0.AddAtom(b)
#    a = ResAtom.fromPDBrecord(aline2)
#    sim0.AddAtom(a)
#    
    sim1 = MySystem.fromXYZfile("test.xyz")
    print sim1.residues["XYZ"]
#    sim1.connectAtoms()
    sim1.groupAtoms()
    sim1.info()

#    sim1 = MySystem.fromPDBfile("test.pdb")
# 
#    sim1.toPDBfile()
