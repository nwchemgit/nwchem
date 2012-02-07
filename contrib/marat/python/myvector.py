'''
Created on Feb 4, 2012

@author: marat
'''
import math
from numpy import array

class Vector:
    def __init__(self, coord=[],name="Vector"):
        self.size = len(coord)
        self.coord = array(coord)

    def __add__(self, value):
        res = []
        for x in self.coord: res.append(x + value)
        return Vector(res)
    __radd__ = __add__

    def __mul__(self, value):
        res = []
        for x in self.coord: res.append(x * value)
        return Vector(res)
    __rmul__ = __mul__
    
    def __repr__(self): 
        return `self.coord`

    def norm1(self, start=0):
        return reduce(lambda x,y: x + abs(y), self.coord, start)
    
    def norm2(self):
        start = 0.0
        rd = reduce(lambda x,y: x + y**2, self.coord, start)
        return math.sqrt(rd)
    
    def length(self):
        return self.norm2()
    
    def prod(self, start=1):
        return reduce(lambda x,y: x * y, self.coord, start)



    
def test():
    x = Vector([2, 4, 6])
    print x + 3, 3 + x
    print x * 4, 4 * x
    print x.norm1(), x.prod()

    y = Vector([1, 2, 3])
    print x + y
    print x * y
    print x * y * 2
    print x * y * x
    
    z = Vector([1,1,1])
    print z.length(), math.sqrt(3.0), z.size

if __name__ == '__main__': test()    # run my self-test code
