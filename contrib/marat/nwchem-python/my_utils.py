'''
Created on Feb 9, 2012

@author: marat
'''

def extract_random_clusters(filename,nres,prefix=None,nconf=15,with_residues=None):
    '''
    extracts random clusters from file
    specified by filename and saves them
    into PDB files specified by prefix
    Input parameters:
    filename - name of coordinate file (xyz or PDB)
    prefix   - prefix for created clusters
    nres     - number of residues 
    nconf=15 - number of configurations to generate
    with_residues = [] - list of residues that should always be included
    '''
    import random
    import string
    from  my_system import MySystem
    
    if with_residues is None:
        with_residues=[]
    if prefix is None:
        prefix = string.strip(filename)[:-4]+"-"+str(nres)  
    
    sim1 = MySystem.fromPDBfile(filename)
    nr = sim1.nres()
    
    reslist=set() 
    while len(reslist)<nconf:
        alist = set(with_residues)
        while len(alist)<nres:
            i=random.randint(0,nr-1)
            alist.add(i)     
        
        reslist.add(tuple(sorted(alist)))
        
    print reslist
    
    for i,s in enumerate(reslist):
        filename = "%s-%d.pdb"%(prefix,i)
        comment="-".join(["%s" % el for el in s])
        sim1.toPDBfile1(filename,s,comment)
        print "generated cluster",s,"as",filename


    
def make_unique(seq, idfun=None):  
    # order preserving 
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

if __name__ == '__main__':
    extract_random_clusters("shell-1.pdb",11,with_residues=[0,1,2])
    
    