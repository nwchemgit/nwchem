'''
Created on Feb 3, 2012

@author: marat
'''

import string
from myvector import Vector

class XYZ_file:
    '''
    classdocs
    '''

    def __init__(self,name=None):
        if name==None:
            self.FileName = " "
            self.NumAtom = 0
            self.AtomName = []
            self.AtomPos = []
            self.AtomVel = []
            self.AtomId = []
            print "am here"
        else:
            print "loading file", name
            self.LoadFile(name)
            
    def LoadFile(self,FileName):
        fp = open(str(FileName),'r')
        self.FileName = FileName
        
        lines = fp.readlines()
        NumAtom = string.atoi(lines[0])
        self.NumAtom = NumAtom
        AtomStr = []*4
        self.AtomName = [None]*NumAtom
        self.AtomId = [None]*NumAtom
        self.AtomPos = [None]*NumAtom
        self.AtomVel = [None]*NumAtom

        x = 0.0
        y = 0.0
        z = 0.0
        
        for i in range(NumAtom):
            self.AtomId[i] = i+1
            AtomStr  = string.split(lines[i+2])[0:4]
            self.AtomName[i] = AtomStr[0]
            x = string.atof(AtomStr[1])
            y = string.atof(AtomStr[2])
            z = string.atof(AtomStr[3])
            self.AtomPos[i] = Vector(x,y,z)
            self.AtomVel[i] = Vector()
            
    def AddAtom(self,Name,x,y,z):
        self.AtomId.append(self.NumAtom+1)
        self.AtomName.append(Name)
        self.AtomPos.append(Vector(x,y,z))
        self.NumAtom = len(self.AtomId)
                
    def MoveAtom(self,i,dr):
        self.AtomPos[i-1] = self.AtomPos[i-1] + dr

    def SetAtomVel(self,i,v):
        self.AtomVel[i-1] = v
        
    def BondVector(self,i1,i2):
        dr = self.AtomPos[i2-1] - self.AtomPos[i1-1]
        return dr

    def BondLength(self,i1,i2):
        dr = self.AtomPos[i2-1] - self.AtomPos[i1-1]
        return dr.length()

    def WriteFile(self,FileName):
        fp = open(str(FileName),'w')
        fp.write(str(self.NumAtom))
        fp.write("\n")
        fp.write("molecule")
        fp.write("\n")       
        for i in range(self.NumAtom):
            fp.write(self.AtomName[i])
            fp.write("   ")
            fp.write(str(self.AtomPos[i]))
#            fp.write("  ")
#            fp.write(str(self.AtomVel[i]))         
            fp.write("\n")

    def AppendFile(self,FileName):
        fp = open(str(FileName),'a')
        fp.write(str(self.NumAtom))
        fp.write("\n")
        fp.write("molecule")
        fp.write("\n")       
        for i in range(self.NumAtom):
            fp.write(self.AtomName[i])
            fp.write("   ")
            fp.write(str(self.AtomPos[i]))
#            fp.write("  ")
#            fp.write(str(self.AtomVel[i]))         
            fp.write("\n")


if __name__ == '__main__':
    a = XYZ_file("test.xyz")
    print a.AtomPos[0]
    print a.BondLength(1, 2)
    a.WriteFile("test1.xyz")
