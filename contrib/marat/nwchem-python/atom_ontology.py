'''
Created on Feb 7, 2012

@author: marat
'''

class AtomOntology:
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    @staticmethod    
    def atomName(mydict):
        if mydict.has_key("name"):
            name = mydict["name"]
        else:
            name = None     
        return name
    
    @staticmethod  
    def resName(mydict):
        if mydict.has_key("resname"):
            name = mydict["resname"]
        else:
            name = None     
        return name
            
    @staticmethod  
    def resId(mydict):
        if mydict.has_key("resid"):
            i = mydict["resid"]
        else:
            i = None     
        return i

    @staticmethod  
    def atomId(mydict):
        if mydict.has_key("atomid"):
            i = mydict["atomid"]
        else:
            i = None     
        return i
    
    @staticmethod 
    def atomCoords(mydict):
        if mydict.has_key("coords"):
            coords = mydict["coords"]
        else:
            coords = None     
        return coords    
    
    
        