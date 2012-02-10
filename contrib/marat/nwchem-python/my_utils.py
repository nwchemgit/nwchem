'''
Created on Feb 9, 2012

@author: marat
'''
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
    print make_unique([1,2,2,2,3,4,5,6,6,6,6,5,1])