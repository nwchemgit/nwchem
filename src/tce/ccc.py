# Ground-state coupled-cluster ansatz generator for OCE
# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)
# $Id: ccc.py,v 1.1 2003-07-25 17:55:26 sohirata Exp $

import string
import copy
import sys

def factorial(n):
   if (n == 0):
      return 1
   else:
      return n*factorial(n-1)
 
def newindex(type,number):
   index = string.join([type,repr(number)],"")
   return index
 
def newtensor(type,indexes):
   tensor = string.join([type,"("],"")
   for index in indexes:
      tensor = string.join([tensor,index])
   tensor = string.join([tensor,")"])
   return tensor

print "\nInput ranks of cluster amplitudes (list integers separated by spaces; default: [0])?"
stringorders = sys.stdin.readline()
stringorders = stringorders[0:len(stringorders)-1]
stringorders = string.split(stringorders)
torders = [0]
for stringorder in stringorders:
   torders.append(int(stringorder))
print " ... Ranks of cluster amplitudes:",torders

print "\nInput ranks of excitation amplitudes (list integers separated by spaces; default: [0])?"
stringorders = sys.stdin.readline()
stringorders = stringorders[0:len(stringorders)-1]
stringorders = string.split(stringorders)
xorders = [0]
for stringorder in stringorders:
   xorders.append(int(stringorder))
print " ... Ranks of excitation amplitudes:",xorders

print "\nR0 = 1 (y or n; default y)?"
yesno = sys.stdin.readline()
if ((yesno[0] == "n") or (yesno[0] == "N")):
   r0is1 = 0
   print " ... no"
else:
   r0is1 = 1
   print " ... yes"

print "\nInput ranks of deexcitation amplitudes (list integers separated by spaces; default: [0])?"
stringorders = sys.stdin.readline()
stringorders = stringorders[0:len(stringorders)-1]
stringorders = string.split(stringorders)
yorders = [0]
for stringorder in stringorders:
   yorders.append(int(stringorder))
print " ... Ranks of deexcitation amplitudes:",yorders

print "\nL0 = 1 (y or n; default y)?"
yesno = sys.stdin.readline()
if ((yesno[0] == "n") or (yesno[0] == "N")):
   l0is1 = 0
   print " ... no"
else:
   l0is1 = 1
   print " ... yes"

print "\nInput a left projection order (give one integer; default: 0)?"
projection = sys.stdin.readline()
if (projection == "\n"):
   leftprojection = 0
else:
   leftprojection = int(projection[0:len(projection)-1])
print " ... Left projection:",leftprojection

print "\nInput a right projection order (give one integer; default: 0)?"
projection = sys.stdin.readline()
if (projection == "\n"):
   rightprojection = 0
else:
   rightprojection = int(projection[0:len(projection)-1])
print " ... Right projection:",rightprojection

print "\nInput ranks of operators for expectation value evaluation (list integers separated by spaces; default: [1,2])?"
stringorders = sys.stdin.readline()
stringorders = stringorders[0:len(stringorders)-1]
stringorders = string.split(stringorders)
horders = []
for stringorder in stringorders:
   horders.append(int(stringorder))
if (not horders):
   horders = [1,2]
print " ... Ranks of operators:",horders

ansatz = "<"+repr(leftprojection)+"|"
yans = ""
for order in yorders:
   if ((order == 0) and l0is1):
      if (yans):
         yans = yans + "+1"
      else:
         yans = yans + " (1"
   else:
      if (yans):
         yans = yans + "+L"+repr(order)
      else:
         yans = yans + " (L"+repr(order)
if (yans):
   ansatz = ansatz + yans + ")"
if (horders == [1,2]):
   ansatz = ansatz + " H"
elif (horders == [1]):
   ansatz = ansatz + " F"
elif (horders == [2]):
   ansatz = ansatz + " V"
tans = ""
for order in torders:
   if (order > 0):
      if (tans):
         tans = tans + "+T"+repr(order)
      else:
         tans = tans + "(T"+repr(order)
if (tans):
   ansatz = ansatz + " exp"+tans+")"
xans = ""
for order in xorders:
   if ((order == 0) and r0is1):
      if (xans):
         xans = xans + "+1"
      else:
         xans = xans + " (1"
   else:
      if (xans):
         xans = xans + "+R"+repr(order)
      else:
         xans = xans + " (R"+repr(order)
if (xans):
   ansatz = ansatz + xans + ")"
ansatz = ansatz + " |"+repr(rightprojection)+">"
print ""
print ansatz
print ""

print "Input a filename?"
filename = sys.stdin.readline()
filename = filename[0:len(filename)-1]
file = open(filename,"w")

for t1 in torders:
   for t2 in torders:
      if (t2 < t1):
         continue
      for t3 in torders:
         if (t3 < t2):
            continue
         for t4 in torders:
            if (t4 < t3):
               continue
            for x in xorders:
               for y in yorders:
                  if ((t1+t2+t3+t4+x+y == 0) and (leftprojection == 0) and (rightprojection == 0)):
                     continue
                  if ((rightprojection+t1+t2+t3+t4+x > leftprojection+y+2) or (rightprojection+t1+t2+t3+t4+x < leftprojection+y-2)):
                     continue
                  for hamiltonian in horders:
                     counter = 0
                     summation = []
                     tensors = []
                     operator = []
                     # left projection
                     if (leftprojection > 0):
                        curly = []
                        for i in range(leftprojection):
                           counter = counter + 1
                           j = newindex('h',counter)
                           curly.append(j+"+")
                        pointer = len(curly)
                        for i in range(leftprojection):
                           counter = counter + 1
                           a = newindex('p',counter)
                           curly.insert(pointer,a)
                        operator.append(curly)   
                     # right projection
#
# 6/18/03 we promoted right projection logic here, just to reserve
#         the same consequtive sets of indexes for the externals.
#         If externals have different indexes, the factorization will
#         break.  The actual insertion of right projection occurs later.
#
                     if (rightprojection > 0):
                        rightcurly = []
                        for i in range(rightprojection):
                           counter = counter + 1
                           j = newindex('p',counter)
                           rightcurly.append(j+"+")
                        pointer = len(rightcurly)
                        for i in range(rightprojection):
                           counter = counter + 1
                           a = newindex('h',counter)
                           rightcurly.insert(pointer,a)
#                    tlinesused = 0
#                    tlinesleft = 0
#                    if (rightprojection > 0):
#                       tlinesused = tlinesused + 1
#                       tlinesleft = tlinesleft + rightprojection*2 - 1
#                    if (t1 > 0):
#                       tlinesused = tlinesused + 1
#                       tlinesleft = tlinesleft + t1*2 - 1
#                    if (t2 > 0):
#                       tlinesused = tlinesused + 1
#                       tlinesleft = tlinesleft + t2*2 - 1
#                    if (t3 > 0):
#                       tlinesused = tlinesused + 1
#                       tlinesleft = tlinesleft + t3*2 - 1
#                    if (t4 > 0):
#                       tlinesused = tlinesused + 1
#                       tlinesleft = tlinesleft + t4*2 - 1
#                    if (x > 0):
#                       tlinesused = tlinesused + 1
#                       tlinesleft = tlinesleft + x*2 - 1
                     factor = 1
                     if ((not l0is1) or (y > 0)):
                        factor = factor * factorial(y) * factorial(y)
                        curly = []
                        indexes = []
                        for i in range(y):
                           counter = counter + 1
                           a = newindex('h',counter)
                           curly.append(a+"+")
                           summation.append(a)
                           indexes.append(a)
                        pointer = len(curly)
                        for i in range(y):
                           counter = counter + 1
                           j = newindex('p',counter)
                           curly.insert(pointer,j)
                           summation.append(j)
                           indexes.append(j)
                        operator.append(curly)   
                        tensors.append(newtensor('y',indexes))
                     if (hamiltonian == 1):
                        # f operator
#                       if (tlinesused > 2):
#                          continue
#                       if ((tlinesleft + (2 - tlinesused) < 2*(leftprojection+y)) or \
#                           (tlinesleft - (2 - tlinesused) > 2*(leftprojection+y))):
#                          continue
                        factor = factor * 1
                        counter = counter + 1
                        g1 = newindex('g',counter)
                        counter = counter + 1
                        g2 = newindex('g',counter)
                        summation.append(g1)
                        summation.append(g2)
                        tensors.append(newtensor('f',[g1,g2]))
                        curly = [g1+"+",g2]
                        operator.append(curly)
                     elif (hamiltonian == 2):
                        # v operator
#                       if (tlinesused > 4):
#                          continue
#                       if ((tlinesleft + (4 - tlinesused) < 2*(leftprojection+y)) or \
#                           (tlinesleft - (4 - tlinesused) > 2*(leftprojection+y))):
#                          continue
                        factor = factor * 4
                        counter = counter + 1
                        g1 = newindex('g',counter)
                        counter = counter + 1
                        g2 = newindex('g',counter)
                        counter = counter + 1
                        g3 = newindex('g',counter)
                        counter = counter + 1
                        g4 = newindex('g',counter)
                        summation.append(g1)
                        summation.append(g2)
                        summation.append(g3)
                        summation.append(g4)
                        tensors.append(newtensor('v',[g1,g2,g3,g4]))
                        curly = [g1+"+",g2+"+",g4,g3]
                        operator.append(curly)
                     else:
                        raise RuntimeError, "unsupported operator rank"
                     list = []
                     if (t1 > 0):
                        list.append(t1)
                        factor = factor * factorial(t1) * factorial(t1)
                        curly = []
                        indexes = []
                        for i in range(t1):
                           counter = counter + 1
                           a = newindex('p',counter)
                           curly.append(a+"+")
                           summation.append(a)
                           indexes.append(a)
                        pointer = len(curly)
                        for i in range(t1):
                           counter = counter + 1
                           j = newindex('h',counter)
                           curly.insert(pointer,j)
                           summation.append(j)
                           indexes.append(j)
                        operator.append(curly)   
                        tensors.append(newtensor('t',indexes))
                     if (t2 > 0):
                        list.append(t2)
                        factor = factor * factorial(t2) * factorial(t2)
                        curly = []
                        indexes = []
                        for i in range(t2):
                           counter = counter + 1
                           a = newindex('p',counter)
                           curly.append(a+"+")
                           summation.append(a)
                           indexes.append(a)
                        pointer = len(curly)
                        for i in range(t2):
                           counter = counter + 1
                           j = newindex('h',counter)
                           curly.insert(pointer,j)
                           summation.append(j)
                           indexes.append(j)
                        operator.append(curly)   
                        tensors.append(newtensor('t',indexes))
                     if (t3 > 0):
                        list.append(t3)
                        factor = factor * factorial(t3) * factorial(t3)
                        curly = []
                        indexes = []
                        for i in range(t3):
                           counter = counter + 1
                           a = newindex('p',counter)
                           curly.append(a+"+")
                           summation.append(a)
                           indexes.append(a)
                        pointer = len(curly)
                        for i in range(t3):
                           counter = counter + 1
                           j = newindex('h',counter)
                           curly.insert(pointer,j)
                           summation.append(j)
                           indexes.append(j)
                        operator.append(curly)   
                        tensors.append(newtensor('t',indexes))
                     if (t4 > 0):
                        list.append(t4)
                        factor = factor * factorial(t4) * factorial(t4)
                        curly = []
                        indexes = []
                        for i in range(t4):
                           counter = counter + 1
                           a = newindex('p',counter)
                           curly.append(a+"+")
                           summation.append(a)
                           indexes.append(a)
                        pointer = len(curly)
                        for i in range(t4):
                           counter = counter + 1
                           j = newindex('h',counter)
                           curly.insert(pointer,j)
                           summation.append(j)
                           indexes.append(j)
                        operator.append(curly)   
                        tensors.append(newtensor('t',indexes))
                     if ((not r0is1) or (x > 0)):
                        factor = factor * factorial(x) * factorial(x)
                        curly = []
                        indexes = []
                        for i in range(x):
                           counter = counter + 1
                           a = newindex('p',counter)
                           curly.append(a+"+")
                           summation.append(a)
                           indexes.append(a)
                        pointer = len(curly)
                        for i in range(x):
                           counter = counter + 1
                           j = newindex('h',counter)
                           curly.insert(pointer,j)
                           summation.append(j)
                           indexes.append(j)
                        operator.append(curly)   
                        tensors.append(newtensor('x',indexes))
                     # right projection
                     if (rightprojection > 0):
                        operator.append(rightcurly)   
                     while (list):
                        i = list[0]
                        n = list.count(i)
                        factor = factor * factorial(n) 
                        numberofi = list.count(i)
                        for dummy in range(numberofi):
                           list.remove(i)
                     ansatz = "1.0/"+repr(factor)+".0"
                     ansatz = string.join([ansatz,"Sum("])
                     for index in summation:
                        ansatz = string.join([ansatz,index])
                     ansatz = string.join([ansatz,")"])
                     for tensor in tensors:
                        ansatz = string.join([ansatz,tensor])
                     for curly in operator:
                        ansatz = string.join([ansatz,"{"])
                        for index in curly:
                           ansatz = string.join([ansatz,index])
                        ansatz = string.join([ansatz,"}"])
                     print ansatz
                     file.write(ansatz)
                     file.write("\n")
