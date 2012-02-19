'''
Created on Feb 5, 2012

@author: marat
'''
from atom_params import *
from generic_atom import *
from generic_residue import *
import numpy
import networkx as nx
import itertools as it
    
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

    def natoms(self):
        return len(self.atoms)

    def nres(self):
        return len(self.reslist)
            
    @classmethod        
    def from_file(cls,filename):
        '''
        alternative constructor from file
        '''
        # figure out format
        fmt=string.strip(filename)[-3:]
        if fmt=="xyz":
            cls=MySystem.fromXYZfile(filename)
            cls.groupAtoms()
        elif fmt=="pdb":
            cls=MySystem.fromPDBfile(filename)
        else:
            raise ValueError("unknown file format")
        
        return cls    
  
    
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
        
    def toPDBfile(self,filename):
        fp = open(str(filename),'w')
        i=0
        for residue in self.reslist:
            i = i+1
            fp.write(residue.toPDBrecord(resid=i))

    def toPDBfile1(self,filename,rlist,comment=None):
        fp = open(str(filename),'w')
        if comment:
            fp.write("#comment: %s\n"%comment)
        i=0
        offset = 0
        for ir in rlist:
            residue = self.reslist[ir]
            i = i+1
            fp.write(residue.toPDBrecord(resid=i,offset=offset))
            offset = offset + residue.size()
        fp.close()  
                      
    def AddAtom(self,a1):
        self.atoms.append(a1)
        rmap = self.residues
        tag = a1.groupTag()
        offset = len(self.atoms)
        resname = a1.res_name()
        if tag not in rmap:
            rmap[tag]=GenericResidue(name=resname,offset=offset)
            self.reslist.append(rmap[tag])
        rmap[tag].AddAtom(a1)

    def create_graph_directed(self,name):
        import networkx as nx
        import matplotlib.pyplot as plt    
        G=nx.MultiDiGraph()
        for i,r in enumerate(self.reslist):
            G.add_node(i+1,name=r.name)
        print G.nodes()
        solvent = [n for n,d in G.nodes_iter(data=True) if d['name'] in ['WAT','HOH','WTR' ]]
        solute = [n for n,d in G.nodes_iter(data=True) if d['name'] not in ['WAT','HOH','WTR' ]]
        print "solute",solute
        
        nr = len(self.reslist)
        for i in range(nr):
            ri = self.reslist[i]
            for j in range(i+1,nr):
                rj = self.reslist[j]
                hbond = GenericResidue.hbonded1_directed(ri, rj)
                if len(hbond) > 1:
                    hname = "hbond2"
                else:
                    hname = "hbond"
                for h in hbond:
                    if h <0:
                        G.add_edge(j+1,i+1,name=hname) 
                    else:
                        G.add_edge(i+1,j+1,name=hname)                         
                if GenericResidue.spec_bonded(ri, rj):
                    G.add_edge(j+1,i+1,name="special")                     
         
        print "all",[(u,v) for u,v,d in G.edges_iter(data=True)]                   
        esolute = [(u,v,d) for u,v,d in G.edges_iter(data=True) if (u in solute or v in solute) and  d['name']=='hbond'  ]
        print "solute-solute",esolute

        esolute2 = [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solute or v in solute) and  d['name']=='hbond2'  ]
        print "double solute",esolute2
                
        esolvent= [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solvent and v in solvent) and  d['name']=='hbond']
        print "solvent-solvent",esolvent       

        especial= [(u,v) for u,v,d in G.edges_iter(data=True) if d['name']=='special' ]
        print especial      
        
        G=nx.relabel_nodes(G,lambda x: str(x))        
        pos0=nx.spectral_layout(G)
        pos=nx.spring_layout(G,iterations=500,pos=pos0)
#        pos=nx.shell_layout(G)
#        pos=nx.graphviz_layout(G,root=1,prog='dot')
        nx.draw(G)
#        nx.draw_networkx_nodes(G,pos,node_size=500,nodelist=solute)
#        nx.draw_networkx_nodes(G,pos,node_size=300,nodelist=solvent,node_color='b')
#
#        nx.draw_networkx_edges(G,pos,edgelist=esolute,
#                            width=3,edge_color='red',style='dashed')
#        nx.draw_networkx_edges(G,pos,edgelist=esolvent,
#                            width=3,edge_color='blue')   
#        nx.draw_networkx_edges(G,pos,edgelist=especial,
#                            width=3,edge_color='c',style='dashed')     
#
#        nx.draw_networkx_edges(G,pos,edgelist=esolute2,
#                            width=8,edge_color='red',style='dashed')             
#              
#        nx.draw_networkx_labels(G,pos)
#             
#        plt.axis('off')
        plt.savefig(name)
#        T=nx.dfs_tree(G)
#        print(sorted(T.edges(data=True)))
        nx.write_dot(G,"shell.dot")
#        print "density=",nx.density(G)
        nx.write_gml(G, "test.gml")
#        plt.show() # display
#        return G


    def create_graph(self,name):
        import networkx as nx
        import matplotlib.pyplot as plt    
        G=nx.MultiGraph()
        for i,r in enumerate(self.reslist):
            G.add_node(i+1,name=r.name)
        print G.nodes()
        solvent = [n for n,d in G.nodes_iter(data=True) if d['name'] in ['WAT','HOH','WTR' ]]
        solute = [n for n,d in G.nodes_iter(data=True) if d['name'] not in ['WAT','HOH','WTR' ]]
        print "solute",solute
        
        h = self.hbond_matrix1()
        nr = numpy.size(h,0)
        for i in range(nr):
            for j in range(i+1,nr):
                if h[i][j]==1:
                    G.add_edge(i+1,j+1,name="hbond") 
                if h[i][j]==2:
                    G.add_edge(i+1,j+1,name="hbond2")     
                    G.add_edge(i+1,j+1,name="hbond2")          
                elif GenericResidue.spec_bonded(self.reslist[i], self.reslist[j]):
                    G.add_edge(i+1,j+1,name="special")                     
         
        print "all",[(u,v) for u,v,d in G.edges_iter(data=True)]                   
        esolute = [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solute or v in solute) and  d['name']=='hbond'  ]
        print esolute

        esolute2 = [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solute or v in solute) and  d['name']=='hbond2'  ]
        print esolute2
                
        esolvent= [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solvent and v in solvent) and  d['name']=='hbond']
        print esolvent       

        especial= [(u,v) for u,v,d in G.edges_iter(data=True) if d['name']=='special' ]
        print especial      
                
        pos0=nx.spectral_layout(G)
        pos=nx.spring_layout(G,iterations=500,pos=pos0)
#        pos=nx.graphviz_layout(G,root=1,prog='dot')
#        nx.draw(G)
        nx.draw_networkx_nodes(G,pos,node_size=500,nodelist=solute)
        nx.draw_networkx_nodes(G,pos,node_size=300,nodelist=solvent,node_color='b')

        nx.draw_networkx_edges(G,pos,edgelist=esolute,
                            width=3,edge_color='red',style='dashed')
        nx.draw_networkx_edges(G,pos,edgelist=esolvent,
                            width=3,edge_color='blue')   
        nx.draw_networkx_edges(G,pos,edgelist=especial,
                            width=3,edge_color='c',style='dashed')     

        nx.draw_networkx_edges(G,pos,edgelist=esolute2,
                            width=8,edge_color='red',style='dashed')             
              
        nx.draw_networkx_labels(G,pos)
             
        plt.axis('off')
        plt.savefig(name)
#        T=nx.dfs_tree(G)
#        print(sorted(T.edges(data=True)))
#        nx.write_dot(G,"shell.dot")
#        print "density=",nx.density(G)
  
#        plt.show() # display
#        return G

    def num_hbonds(self,**kwargs):
        '''
        '''
        nu=0
        nv=0
        h = self.hbond_matrix1(**kwargs)
        nr = numpy.size(h,0)
        for i in range(nr):
            for j in range(i+1,nr):
                if (self.reslist[i].name in  ['WAT','HOH','WTR' ] 
                    and  self.reslist[j].name in  ['WAT','HOH','WTR' ]):
                    nv = nv + h[i][j]         
                else:
                    nu = nu + h[i][j]  
        return nu+nv,nu,nv
                        
    def graph(self):
        '''
        '''
        import networkx as nx
        G=nx.MultiGraph()
        for i,r in enumerate(self.reslist):
            G.add_node(i+1,name=r.name)

        solvent = [n for n,d in G.nodes_iter(data=True) if d['name'] in ['WAT','HOH','WTR' ]]
        solute = [n for n,d in G.nodes_iter(data=True) if d['name'] not in ['WAT','HOH','WTR' ]]
        
        h = self.hbond_matrix1()
        nr = numpy.size(h,0)
        for i in range(nr):
            for j in range(i+1,nr):
                if h[i][j]==1:
                    G.add_edge(i+1,j+1,name="hbond") 
                if h[i][j]==2:
                    G.add_edge(i+1,j+1,name="hbond2")     
                    G.add_edge(i+1,j+1,name="hbond2")          
                elif GenericResidue.spec_bonded(self.reslist[i], self.reslist[j]):
                    G.add_edge(i+1,j+1,name="special")                     
                        
        esolute = [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solute or v in solute) and  d['name']=='hbond'  ]
        esolute2 = [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solute or v in solute) and  d['name']=='hbond2'  ]                
        esolvent= [(u,v) for u,v,d in G.edges_iter(data=True) if (u in solvent and v in solvent) and  d['name']=='hbond']       
        especial= [(u,v) for u,v,d in G.edges_iter(data=True) if d['name']=='special' ]  
                
        return G
  
            
    def create_simple_graph(self,name):
        import networkx as nx
        G=nx.Graph()
        for i,r in enumerate(self.reslist):
            G.add_node(i+1,name=r.name)
        print G.nodes()

        
        h = self.hbond_matrix()
        nr = numpy.size(h,0)
        for i in range(nr):
            for j in range(i+1,nr):
                if GenericResidue.touching(self.reslist[i], self.reslist[j],2.8):
                    print 'adding edge',i+1,j+1
                    G.add_edge(i+1,j+1,name="special")                     
 
#        pos0=nx.spectral_layout(G)
#        pos=nx.spring_layout(G,iterations=500,pos=pos0)

#        pos=nx.spring_layout(G,iterations=500)
        
#        nx.draw_shell(G)

        pos0=nx.spectral_layout(G)
        pos=nx.spring_layout(G,iterations=500,pos=pos0)
#        pos=nx.graphviz_layout(G,root=1,prog='dot')
#        nx.draw(G)
        nx.draw_networkx(G,pos)
                     
        plt.axis('off')
        plt.savefig(name)       
        

    def connectAtoms(self):
        for i in range(len(self.atoms)):
            for j in range(i+1,len(self.atoms)):
                a1=self.atoms[i]
                a2=self.atoms[j]
                if GenericAtom.bonded(a1, a2):
                    a1.setBond(j)
                    a2.setBond(i)
                                   

    def hbond_matrix(self):
        nr = len(self.reslist)
        hbond = numpy.zeros(shape=(nr,nr),dtype=int)    
    
        for i in range(nr):
            ri=self.reslist[i]   
            for j in range(i+1,nr):
                rj=self.reslist[j]
                hbond[i][j]=GenericResidue.hbonded(ri,rj)
                hbond[j][i]=hbond[i][j]
        return hbond

    def hbond_matrix1(self,**kwargs):
        nr = len(self.reslist)
        hbond = numpy.zeros(shape=(nr,nr),dtype=int)    
        for i in range(nr):
            ri=self.reslist[i]   
            for j in range(i+1,nr):
                rj=self.reslist[j]
                hbond[i][j]=GenericResidue.hbonded1(ri,rj,**kwargs)
#                hbond[i][j]=hbonded(ri,rj)
                hbond[j][i]=hbond[i][j]
        return hbond
    
    def dist_matrix(self,rcut):
        nr = len(self.reslist)
        dm = numpy.zeros(shape=(nr,nr),dtype=int)    
    
        for i in range(nr):
            ri=self.reslist[i]   
            for j in range(i+1,nr):
                rj=self.reslist[j]
                dm[i][j]=GenericResidue.touching(self.reslist[i], self.reslist[j],2.8)
                dm[j][i]=dm[i][j]
        return dm         

                                        
    def groupAtoms(self):
        rmap = self.residues
        reslist=self.reslist
        ir = 0
        offset=0
        for i in range(len(self.atoms)):
            a1=self.atoms[i]
            if a1.groupTag()!="XYZ":
                continue
            ir = ir +1
            tag="00"+str(ir)
            tag=tag[-3:]
            a1.setGroupTag(tag)
            rmap[tag]=GenericResidue(name=tag,offset=offset)
            rmap[tag].AddAtom(a1)
            offset = offset + 1
            reslist.append(rmap[tag])
            for j in range(i+1,len(self.atoms)):                
                a2=self.atoms[j]
                if GenericAtom.bonded(a1, a2):
                    a2.setGroupTag(tag)
                    rmap[tag].AddAtom(a2)
                    offset = offset + 1
                    
            for r in self.reslist:
                r.guess_name()
                    
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
    import matplotlib.pyplot as plt    
    import random
    
#    
    sim1 = MySystem.from_file("w10-11.xyz")
    print sim1.num_hbonds(rOH=2.18, HOH=138)
#    sim1.create_graph_directed("test.png")
#    sim1.toPDBfile("shell-1.pdb")       
#    sim1.create_graph("shell.png")
#    G=sim1.graph()
#    print nx.density(G)
#    print sim1.num_hbonds()[0]
#    dm=sim1.dist_matrix(2.8)
#    chain0=[0]
#    level=2
#    chain1=[]
#    print dm
#    for i in chain0:
#        for j in range(i+1,25):
#            if dm[i][j]==1:
#                chain1.append([i,j])
#    print chain1
#    
#    chain0=chain1
#    chain1=[]
#
#    for i in chain0[0]:
#        for j in range(i+1,25):
#            if dm[i][j]==1:
#                chain1.append([i,j])
#    print chain1
# 
#    slist=set() 
#    while len(slist)<15:
#        alist = set()
#        while len(alist)<7:
#            i=random.randint(3, 20)
#            alist.add(i)     
#        
#        slist.add((0,1,2)+tuple(sorted(alist)))
#        
#    for i,s in enumerate(slist):
#        filename = "cw9-%d.pdb"%(i)
#        comment="-".join(["%s" % el for el in s])
#        print s,filename,comment
#        sim1.toPDBfile1(filename,s,comment)

#    it = iter(slist)
#    sim1.toPDBfile1("test-123.pdb",next(it),comment="my best trs")
    
#
#    nx.draw_spring(H)
#    plt.show()

    print "finished"
    sys.exit(0)
    

#    sim1 = MySystem.fromXYZfile("w10-1.xyz")
#    sim1.groupAtoms()
#    sim2 = MySystem.fromPDBfile("test.pdb")
#    sim3 = MySystem.fromXYZfile("w2-test.xyz")
#    sim3.groupAtoms()
#    print sim1.residues["XYZ"]
#    sim1.connectAtoms()

    
#    G=nx.Graph(hbond)
#
#    for i in range(nr):
#        print i+1, hbond[i], sum(hbond[i])
#        
#    print nx.connected_components(G)
#    print nx.clustering(G)
#    nx.write_dot(G,'file.dot')
#    pos=nx.spring_layout(G)
#    colors=range(20)
 #   nx.draw(G,pos,node_color='#A0CBE2',edge_cmap=plt.cm.Blues,with_labels=True)
#    nx.draw_spring(G)
#    plt.show()
#    plt.savefig("path.png")
#    print "looking for cliques"
#    print list(nx.find_cliques(G))
#    print nx.number_connected_components(G)
#    print sorted(nx.degree(G).values())
#    
#    sim1.toPDBfile("mytest.pdb")    
#    try:
#        import numpy.linalg
#        eigenvalues=numpy.linalg.eigvals
#    except ImportError:
#        raise ImportError("numpy can not be imported.")
#
#    L=nx.generalized_laplacian(G)
#    e=eigenvalues(L)
#    print e
#    print "pagerank"
#    d= nx.degree(G)
#    for w in sorted(d, key=d.get, reverse=True):
#      print w, d[w]

#    T=nx.minimum_spanning_tree(G)
#    print(sorted(T.edges(data=True)))
#    print(nx.cycle_basis(G))
 #   nx.draw_spring(G)
 #   plt.show()
#    print len(hbond)
    
#    sim1.info()
#    r=list(sim1.residues.itervalues())
#    print r

#    sim1 = MySystem.fromPDBfile("test.pdb")
# 
#    sim1.toPDBfile()
