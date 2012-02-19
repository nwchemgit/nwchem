'''
Created on Feb 5, 2012

@author: marat
'''
import sys
from generic_atom import *

class GenericResidue(object):
    '''
    classdocs
    '''

    dict_name=dict(HHO="WAT",IOOO="IO3")
    
    def __init__(self,atoms=None,name="UNK",offset=0):
        '''
        Default constructor for Residue class
        atoms list of atoms in the residue
        name residue name
        '''
        self.name = name
        self.offset = offset
        if atoms:
            self.atoms = atoms
        else:
            self.atoms=[]
            
    def signature(self):
        s=sorted(x.elemName() for x in self.atoms)
        return ''.join(s)
            
    def guess_name(self):
        s=sorted(x.elemName() for x in self.atoms)
        name = GenericResidue.dict_name.get(self.signature(),None)
        if name:
            self.name = name
            
    def toPDBrecord(self,resid,offset=None):
 
        a1=" "
        a2=2*" "
        a3=3*" "
        a4=4*" "
        a5=5*" "
        a6=6*" "
        a8=8*" "
        a22=22*" "
 
        pdbformat="%-6s%5s%1s%4s%1s%3s%1s%4s%4s%8.3f%8.3f%8.3f%22s%2s"
        pdbformat="%-6s%5d%1s%4s%1s%3s%1s%4s%4s%8.3f%8.3f%8.3f%22s%2s"
        
        output=''
        if offset is not None:
            atomid=offset
        else:
            atomid=self.offset
        for a in self.atoms:
            atomid = atomid + 1
            output = output + pdbformat%('ATOM',atomid,a1,a.name(),a1,
                                           self.name,a2,str(resid),a4,
                                           a.coord[0],a.coord[1],a.coord[2],a22,
                                           a.elemName())  + "\n"

        return output
    
    @classmethod        
    def fromPDBfile(cls,filename):
        '''
        alternative constructor from PDB file
        '''
        cls = GenericResidue()
        fp = open(str(filename),'r')
        
        for line in fp.readlines():
            if line.startswith('ATOM'):
                a=GenericAtom.fromPDBrecord(line)
                cls.AddAtom(a)
        fp.close
        return cls            

    def AddAtom(self,a):
        self.atoms.append(a)
        
    def delAtom(self,a):
        self.atoms.remove(a)
        
#    def __str__(self):
#        output = ""
#        for a in self.atoms:
#            output = output + str(a) + "\n"
#        return output

    def connectAtoms(self):
        for i in range(len(self.atoms)):
            for j in range(i+1,len(self.atoms)):
                a1=self.atoms[i]
                a2=self.atoms[j]
                if GenericAtom.bonded(a1, a2):
                    a1.setBond(j)
                    a2.setBond(i)
                    
    def get_bonded(self,a0,elem=None):
        al =[]
        for a in self.byFilter():
            if a!=a0 and GenericAtom.bonded(a, a0):
                al.append(a)
        return al
                                    
    @staticmethod
    def distance(res1,res2):
        rmin=100
        for a1 in res1.atoms:
            for a2 in res2.atoms:
                r = GenericAtom.bondlength(a1, a2)
                if r < rmin:
                    a1_min = a1
                    a2_min = a2
                    rmin = r
        return rmin,a1_min,a2_min
    
    @staticmethod
    def hbonded(res1,res2):
#        rOH=2.0
#        OHO=143
        elems = set(['O', 'H'])
        rOH=2.27
        OHO=138
#       find minimum distance pair
        rmin=100
        for a1 in res1.atoms:
            el1=a1.elemName() 
            if el1 in elems:
                for a2 in res2.atoms:
                    el2=a2.elemName()
                    if el2 in elems and el1!=el2:
                        r = GenericAtom.bondlength(a1, a2)
                        if r < rmin:
                            a1_min = a1
                            a2_min = a2
                            rmin = r
        a1=a1_min
        a2=a2_min
        r=rmin                

        if r > rOH:
            return False
        if a1.elemName()=='H':
            res1,res2=res2,res1
            a1,a2=a2,a1
        elif a2.elemName()=='H':
            pass
        else:
            return False
        a3 = res2.get_bonded(a2, 'O')[0]
        angle = GenericAtom.angle(a1, a2, a3)
        return angle>OHO

    
    @staticmethod
    def hbonded1(res1,res2,**kwargs):
#        rOH=2.0
#        OHO=143
        rOH = kwargs.get('rOH',2.27)
        OHO = kwargs.get('OHO',138)
#        rOH=2.27
#        OHO=138        
        elems = set(['O', 'H'])
        ohbond = False
        nhb=0
#       find minimum distance pair
        plist=[]
        for a1 in res1.atoms:
            el1=a1.elemName() 
            if el1 in elems:
                for a2 in res2.atoms:
                    el2=a2.elemName()
                    if el2 in elems and el1!=el2:
                        r = GenericAtom.bondlength(a1, a2)
                        if r < rOH:
                            plist.append([a1,a2])
#        print "plist",plist
        for a1,a2 in plist:
            if a1.elemName()=='H':
                res1,res2=res2,res1
                a1,a2=a2,a1
            elif a2.elemName()=='H':
                pass
            else:
                continue
            a3 = res2.get_bonded(a2, 'O')[0]
            angle = GenericAtom.angle(a1, a2, a3)
            if angle>OHO:
                hbond = True
                nhb = nhb +1
                
        return nhb

    @staticmethod
    def hbonded1_directed(res1,res2):
#        rOH=2.0
#        OHO=143
        elems = set(['O', 'H'])
        rOH=2.27
        OHO=138
        ohbond = False
        nhb=0
#       find minimum distance pair
        plist=[]
        hbond = []
        for a1 in res1.atoms:
            el1=a1.elemName() 
            if el1 in elems:
                for a2 in res2.atoms:
                    el2=a2.elemName()
                    if el2 in elems and el1!=el2:
                        r = GenericAtom.bondlength(a1, a2)
                        if r < rOH:
                            plist.append([a1,a2])
#        print "plist",plist
        for a1,a2 in plist:
            if a1.elemName()=='H':
                res1,res2=res2,res1
                a1,a2=a2,a1
                h=1
            elif a2.elemName()=='H':
                h=-1
            else:
                continue
            a3 = res2.get_bonded(a2, 'O')[0]
            angle = GenericAtom.angle(a1, a2, a3)
            if angle>OHO:
                hbond.append(h)
                
        return hbond
    

    @staticmethod
    def spec_bonded(res1,res2):
        if set([res1.name,res2.name])!=set(["IO3","WAT"]):
            return False
        elems = set(['O', 'I'])
#       find minimum distance pair
        rmin=100
        for a1 in res1.atoms:
            el1=a1.elemName() 
            if el1 in elems:
                for a2 in res2.atoms:
                    el2=a2.elemName()
                    if el2 in elems and el1!=el2:
                        r = GenericAtom.bondlength(a1, a2)
                        if r < rmin:
                            rmin = r
 
        return rmin < 2.9                    

    @staticmethod
    def touching(res1,res2,rcut):
        return rcut > GenericResidue.distance(res1, res2)[0]
                

#    @staticmethod
#    def hbonded1(res1,res2):
##        rOH=2.0
##        OHO=143
#        rOH=2.27
#        OHO=138
#        (r,a1,a2)=GenericResidue.distance(res1, res2)
#        if r > rOH:
#            return False
#        if a1.elemName()=='H':
#            res1,res2=res2,res1
#            a1,a2=a2,a1
#        elif a2.elemName()=='H':
#            pass
#        else:
#            return False
#        a3 = res2.get_bonded(a2, 'O')[0]
#        angle = GenericAtom.angle(a1, a2, a3)
#        print r,angle
#        return angle>OHO
        
    @staticmethod
    def distance1(res1,res2):
        dr=100
        for a1 in res1.byFilter():
            for a2 in res2.byFilter():
                dr=min(dr,GenericAtom.bondlength(a1, a2))
                print dr,a1.name(),a2.name()
        return dr
        
    def byElement(self,name):
        for a in self.atoms:
            if a.elemName() is name:
                yield a

    def byFilter(self,elem=None):
        for a in self.atoms:
            if elem is None:
                yield a
            elif a.elemName() is elem:
                yield a
                            
    def size(self):
        return len(self.atoms)
    
if __name__ == '__main__':
#    aline1 = "ATOM      3  O2  IO3     1      -1.182   1.410   0.573       -0.80     O"
#    aline2 = "ATOM      1  I1  IO3     1      -1.555  -0.350   0.333        1.39     I"
#
#    res0 = GenericResidue()
#    print res0
#    a = GenericAtom.fromPDBrecord(aline2)
#    print a
#    res0.AddAtom(a)
#    print res0.size()
    
    res0 = GenericResidue.fromPDBfile("io3.pdb")
    print res0
 
    res1 = GenericResidue.fromPDBfile("h2o-1.pdb")
    print res1.signature()
    print res0.signature()
    res1.guess_name()
    print res1.name
    res0.guess_name()
    print res0.name
    
    print "distance test"
    (r,a1,a2)=GenericResidue.distance(res0, res1)
    print r, a1.name(), a2.name()
    print res1.get_bonded(a2, "O")
    name = None
    print (filter(lambda a: name is None or a.elemName()==name,res1.atoms ))
    print GenericResidue.hbonded(res0,res1)
    print "HERE COMES PDB RECORD"
    print res0.toPDBrecord(1)
#    b = ResAtom.fromPDBrecord(aline1)
#    res0.AddAtom(a)
#    res0.AddAtom(b)
#    print res0.toPDBrecord(atom_start=1)