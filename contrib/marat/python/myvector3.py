'''
Created on Feb 4, 2012

@author: marat
'''
import math

class Vector(object):
    '''
    classdocs
    '''


    def __init__(self,x=0.0,y=0.0,z=0.0):
        '''
        Constructor
        '''
        self.coords=[x,y,z]
        self.x=x
        self.y=y
        self.z=z
    
    def length(self):
        dr = 0.0
        for x in self.coords:
            dr = dr + x**2
        dr = math.sqrt(dr)
        return dr

    def __str__(self):
        output = "%12.6f %12.6f %12.6f"%(self.coords[0],self.coords[1],self.coords[2])
        return output
    
    def __repr__(self):
        output = "Vector(%g,%g,%g)"%(self.coords[0],self.coords[1],self.coords[2])
        return output
        
    def __add__(self,y):
        z = Vector()
        for i in range(len(self.coords)):
            z.coords[i] = self.coords[i]+y.coords[i]
        return z        
    
    def __sub__(self,y):
        z = Vector()
        for i in range(len(self.coords)):
            z.coords[i] = self.coords[i]-y.coords[i]
        return z
     

if __name__ == '__main__':
    a = Vector(1.0,1.0,0.0)
    b = Vector(1.0,1.0,0.0)
    print a
    c=a-b
    print c
    print "Length ", c.length()