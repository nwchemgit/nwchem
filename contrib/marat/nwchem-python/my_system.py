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
        self.reslist=[]
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

        fp = open(str(filename),'r')
        
        for line in fp.readlines():
            if line.startswith('ATOM'):
                a=GenericAtom.fromPDBrecord(line)
                cls.AddAtom(a)
#        cls.reslist=list(cls.residues.itervalues())
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
        cls.reslist=[]
        cls.residues={}
        return cls
        
    def toPDBfile(self):
        for residue in self.residues.itervalues():
            print residue
            
    def AddAtom(self,a1):
        self.atoms.append(a1)
        rmap = self.residues
        tag = a1.groupTag()
        if tag not in rmap:
            rmap[tag]=GenericResidue(name=tag)
            self.reslist.append(rmap[tag])
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
        reslist=self.reslist
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
            rmap[tag]=GenericResidue(name=tag)
            rmap[tag].AddAtom(a1)
            reslist.append(rmap[tag])
            for j in range(i+1,len(self.atoms)):                
                a2=self.atoms[j]
                if GenericAtom.bonded(a1, a2):
                    a2.setGroupTag(tag)
                    rmap[tag].AddAtom(a2)
               
                    
    def info(self):
        for tag,res in self.residues.iteritems():
            print tag 
            print res       
    
if __name__ == '__main__':
#    sim0 = MySystem("test")
#    sim1 = MySystem.fromPDBfile("test.pdb")  
#    sys.exit()
    import numpy
    import pygraphviz as pgv
    import networkx as nx
    


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
    sim1 = MySystem.fromXYZfile("w10-1.xyz")
    sim1.groupAtoms()
    sim2 = MySystem.fromPDBfile("test.pdb")
    sim3 = MySystem.fromXYZfile("w2-test.xyz")
    sim3.groupAtoms()
#    print sim1.residues["XYZ"]
#    sim1.connectAtoms()

    
    sim=sim1
    nr = len(sim.reslist)
    print nr
    hbond = numpy.zeros(shape=(nr,nr),dtype=int)


    for i in range(len(sim.reslist)):
        ri=sim.reslist[i]
        print ri.name     
        for j in range(i+1,len(sim.reslist)):
            rj=sim.reslist[j]
            hbond[i][j]=GenericResidue.hbonded(ri,rj)
            hbond[j][i]=hbond[i][j]

    
    G=nx.Graph(hbond)

    for i in range(nr):
        print i+1, hbond[i], sum(hbond[i])
        
    print nx.connected_components(G)
    print nx.clustering(G)
    nx.write_dot(G,'file.dot')
    print "looking for cliques"
    print list(nx.find_cliques(G))
    print nx.number_connected_components(G)
    print sorted(nx.degree(G).values())
        
#    print len(hbond)
    
#    sim1.info()
#    r=list(sim1.residues.itervalues())
#    print r

#    sim1 = MySystem.fromPDBfile("test.pdb")
# 
#    sim1.toPDBfile()
