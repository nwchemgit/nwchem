'''
Created on Feb 7, 2012

@author: marat
'''
import sys
from types import *

class GenericAtom:
    '''
    classdocs
    '''


    def __init__(self,d=None):
        '''
        Constructor
        '''
        if d is None:
            self.dict  = {}
        else:
            if type(d) is not type({}):
                print "wrong type ", type(d)
                print "expecting", type({})
                sys.exit(1)
            else:
                self.dict  = d
            
if __name__ == '__main__':
    a=GenericAtom({})
    print a
        