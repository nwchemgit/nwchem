# Tensor Contraction Engine v.1.0
# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)
# $Id$

import string
import types
import copy

def readfromfile(filename):
   """Converts the content of a file to a list of tensor contractions"""
   newlist = ListTensorContractions()
   file = open(filename,"r")
   alwaystrue = 1
   while (alwaystrue):
      selfexpr = file.readline()
      if (selfexpr == ""):
         file.close()
         return newlist
      else:
         selfexpr = string.split(selfexpr)
         donewithfactors = 0
         pointer = 1
         coefficients = []
         permutations = []
         while (not donewithfactors):
            if (selfexpr[pointer] == "+"):
               parity = 1
               pointer = pointer + 1
            elif (selfexpr[pointer] == "-"):
               parity = -1
               pointer = pointer + 1
            else:
               # neither "+" or "-"; assume "+" and do not increment pointer
               parity = 1
            coefficients.append(string.atof(selfexpr[pointer]) * parity)
            pointer = pointer + 1
            if (selfexpr[pointer] == "]"):
               permutations.append([])
               donewithfactors = 1
            elif ((selfexpr[pointer] == "*") or (selfexpr[pointer] == "P") or (selfexpr[pointer] == "p")):
               if (selfexpr[pointer] == "*"): 
                  pointer = pointer + 1
               indexes = []
               while (selfexpr[pointer + 1] != ")"):
                  pointer = pointer + 1
                  if (selfexpr[pointer] != "=>"):
                     if (selfexpr[pointer][0] == "h"):
                        type = "hole"
                     elif (selfexpr[pointer][0] == "p"):
                        type = "particle"
                     else:
                        type = "general"
                     label = string.atoi(selfexpr[pointer][1:])
                     newindex = Index(type,label)
                     indexes.append(newindex)
               permutations.append(indexes)
               pointer = pointer + 2
               if (selfexpr[pointer] == "]"):
                  donewithfactors = 1
            else:
               permutations.append([])
         factor = Factor(coefficients,permutations)
         summation = []
         if ("Sum" in selfexpr):
            indexes = []
            pointer = selfexpr.index("Sum")
            while (selfexpr[pointer + 1] != ")"):
               pointer = pointer + 1
               if (selfexpr[pointer] != "("):
                  if (selfexpr[pointer][0] == "h"):
                     type = "hole"
                  elif (selfexpr[pointer][0] == "p"):
                     type = "particle"
                  else:
                     type = "general"
                  label = string.atoi(selfexpr[pointer][1:])
                  newindex = Index(type,label)
                  indexes.append(newindex)
            summation = Summation(indexes)
         tensors=[]
         ntensors = selfexpr[selfexpr.index("]"):].count("*")
         if (ntensors != selfexpr[selfexpr.index("]"):].count("(")):
            return "Wrong input format"
         if ("Sum" in selfexpr):
            ntensors = ntensors - 1
            offset = 2
         else:
            offset = 1
         if (ntensors > 0):
            for itensor in range(0,ntensors):
               counter = 0
               for pointer in range(selfexpr.index("]"),len(selfexpr)):
                  if (selfexpr[pointer] == "*"):
                     counter = counter + 1
                     if (counter == itensor + offset):
                        if (selfexpr[pointer+1][len(selfexpr[pointer+1])-1] == "+"):
                           tensortype = selfexpr[pointer + 1][0:len(selfexpr[pointer+1])-1]
                           tensorconj = 1
                        else:
                           tensortype = selfexpr[pointer + 1]
                           tensorconj = 0
                        tensorlabel = itensor + 1
                        anotherpointer = pointer + 2
                        indexes = []
                        while (selfexpr[anotherpointer + 1] != ")"):
                           anotherpointer = anotherpointer + 1
                           if (selfexpr[anotherpointer] != "("):
                              if (selfexpr[anotherpointer][0] == "h"):
                                  type = "hole"
                              elif (selfexpr[anotherpointer][0] == "p"):
                                  type = "particle"
                              else:
                                  type = "general"
                              label = string.atoi(selfexpr[anotherpointer][1:])
                              indexes.append(Index(type,label))
               tensors.append(Tensor(tensortype,indexes,tensorlabel,tensorconj))
      newlist.list.append(TensorContraction(factor,summation,tensors))

def rationaltofractional(x):
   """Returns a fractional num/den that equals to a rational number x"""
   fracx = x - int(x)
   if (fracx < 0.0):
      sign = -1
   else:
      sign = 1
   xx = abs(fracx)
   alwaystrue = 1
   den = 1
   num1 = 0
   num2 = 1
   while (alwaystrue):
      newnum2 = -1
      for num in range(num1,num2+1): 
         y = float(num)/float(den)
         if (abs(y-xx) < 0.0000001):
            return [num*sign+int(x)*den,den]
         elif (y < xx):
            newnum1 = num
         elif (y > xx and (newnum2 == -1)):
            newnum2 = num
      den = den + 1
      num1 = max(0,newnum1)
      num2 = min(den,newnum2+1)
      if (den > 1000):
         print x
         raise RuntimeError, "irrational number"

def factorial(n):
   """Returns a factorial of an integer n"""
   if (n == 0):
      return 1
   else:
      return n*factorial(n-1)

def permutation(n):
   """Returns a list of all permutation of n integers"""
   if (n == 1):
      result = [[1]]
      return result
   else:
      result = permutation(n-1)
      newresult = []
      for shorterpermutation in result:
         for position in range(0,n):
            newpermutation = copy.deepcopy(shorterpermutation)
            newpermutation.insert(position,n)
            newresult.append(newpermutation)
      return newresult

def permutationwithparity(n):
   """Returns a list of all permutation of n integers, with its first element being the parity"""
   if (n == 1):
      result = [[1,1]]
      return result
   else:
      result = permutationwithparity(n-1)
      newresult = []
      for shorterpermutation in result:
         for position in range(1,n+1):
            parity = shorterpermutation[0]
            for swaps in range(n-position):
               parity = - parity
            newpermutation = copy.deepcopy(shorterpermutation)
            newpermutation.insert(position,n)
            newpermutation[0] = parity
            newresult.append(newpermutation)
      return newresult

def restrictedpermutationwithparity(lista,listb,listc):
   """Returns all permutations of a combined list of three input lists of indexes except those change the orders among the input lists"""
   result = []
   combined = lista + listb + listc
   if (len(combined) < 2):
      return [[1,"empty"]]
   permutations = permutationwithparity(len(combined))
   for permutation in permutations:
      newpermutation = []
      newpermutation.append(permutation[0])
      for position in range(1,len(permutation)):
         newpermutation.append(combined[permutation[position]-1])
      rejected = 0
      for nindexa in range(1,len(newpermutation)):
         for nindexb in range(1,len(newpermutation)):
            if (nindexa < nindexb):
               indexa = newpermutation[nindexa]
               indexb = newpermutation[nindexb]
# DEEXCITATION EXTENSION FROM HERE ...
#              if (indexa.type != indexb.type):
               if (indexa.isparticle() and indexb.ishole()):
                  rejected = 1
# ... TO HERE
               if (indexa.isin(lista) and indexb.isin(lista) and indexa.isgreaterthan(indexb)):
                  rejected = 1
               if (indexa.isin(listb) and indexb.isin(listb) and indexa.isgreaterthan(indexb)):
                  rejected = 1
               if (indexa.isin(listc) and indexb.isin(listc) and indexa.isgreaterthan(indexb)):
                  rejected = 1
      if (not rejected):
         result.append(newpermutation)

   if (not result):
      return [[1,"empty"]]
   else:
      return result

def arethesamepermutation(permutationa,permutationb):
   """Returns true if two permutations are identical"""
   # input parameters are lists of indexes objects
   # [a,b,c,c,a,b] means the permutation a -> c; b -> a ; c -> b
   # so, [a,b,c,c,a,b] == [b,c,a,a,b,c]
   if ((permutationa == []) and (permutationb == [])):
      return 1
   if (permutationa == []):
      for nindex in range(len(permutationb)/2):
         if (not permutationb[nindex].isidenticalto(permutationb[nindex+len(permutationb)/2])):
            return 0
      return 1
   if (permutationb == []):
      for nindex in range(len(permutationa)/2):
         if (not permutationa[nindex].isidenticalto(permutationa[nindex+len(permutationa)/2])):
            print permutationa[nindex].show(),permutationa[nindex+len(permutationa)/2].show()
            return 0
      return 1
   if (len(permutationa) != len(permutationb)):
      return 0
   for nindexa in range(len(permutationa)/2):
      indexa = permutationa[nindexa]
      indexx = permutationa[nindexa+len(permutationa)/2]
      found = 0
      for nindexb in range(len(permutationb)/2):
         indexb = permutationb[nindexb]
         indexy = permutationb[nindexb+len(permutationb)/2]
         if ((indexa.isidenticalto(indexb)) and (indexx.isidenticalto(indexy))):
            found = 1
      if (not found):
         return 0
   return 1
 
def combinepermutations(one,two):
   """Connects two permutations of indexes"""
   if (len(one) == 0):
      return two
   elif (len(two) == 0):
      return one
   elif (len(one) != len(two)):
      return "Internal error"
   three = []
   for n in range(len(one)/2):
      three.append(one[n])
   for n in range(len(one)/2,len(one)):
      for m in range(len(two)/2):
         if (one[n].isidenticalto(two[m])):
            three.append(two[m+len(two)/2])
   return three

def performpermutation(list,permutation,reverse=0):
   """Performs a permutation operation to a list of indexes"""
   result = []
   for indexa in list:
      exist = 0
      if (reverse):
         for nindexb in range(len(permutation)/2,len(permutation)):
            indexb = permutation[nindexb]
            if (indexa.isidenticalto(indexb)):
               result.append(permutation[nindexb-len(permutation)/2])
               exist = 1
      else:
         for nindexb in range(len(permutation)/2):
            indexb = permutation[nindexb]
            if (indexa.isidenticalto(indexb)):
               result.append(permutation[nindexb+len(permutation)/2])
               exist = 1
      if (not exist):
         result.append(indexa)

   return result

def parityofpermutation(permutation):
   """Returns the parity of a given permutation"""
   origin = permutation[0:len(permutation)/2]
   destination = permutation[len(permutation)/2:len(permutation)]
   allpermutations = permutationwithparity(len(permutation)/2)
   for onepermutation in allpermutations:
      found = 1
      for i in range(len(permutation)/2):
         if (not origin[i].isidenticalto(destination[onepermutation[i+1]-1])):
            found = 0
      if (found):
         return onepermutation[0]
   raise RuntimeError, "permutation not found"

def arethesamelists(list1,list2):
   """Returns true if two lists of indexes are the same (permutation allowed)"""
   if (len(list1) != len(list2)):
      return 0
   copylist2 = copy.deepcopy(list2)
   for index in list1:
      exist = 0
      for nindex in range(len(copylist2)):
         if (index.isidenticalto(copylist2[nindex])):
            exist = 1
            del copylist2[nindex]
      if (not exist):
         return 0
   return 1

def picknfromlist(n,list):
   """Returns a list of all possible n choices from a given list"""
   integerlist = []
   for i in range(len(list)):
      integerlist.append(i)
   integerchoices = [[]]
   for i in range(n):
      integerchoices = pick1fromlist(integerchoices,integerlist)
   result = []
   for integerchoice in integerchoices:
      choice = []
      for i in integerchoice:
         choice.append(copy.deepcopy(list[i]))
      result.append(choice)
   return result

def pick1fromlist(choices,list):
   """Appends one additional non-overlapping choice of an integer from the list to the exisiting choices"""
   newchoices = []
   for choice in choices:
      for i in list:
         overlap = 0
         for j in choice:
            if (i <= j):
               overlap = 1
         if (not overlap):
            newchoice = copy.deepcopy(choice)
            newchoice.append(i)
            newchoices.append(newchoice)
          
   return newchoices

def sortindexes(self):
   """Sorts a list of indexes in an ascending order with no regard to parity"""
   selfcopy = copy.deepcopy(self)
   alwaystrue = 1
   while (alwaystrue):
      done = 1
      for nindexa in range(len(selfcopy)):
         if (not done):
            break
         indexa = selfcopy[nindexa]
         for nindexb in range(len(selfcopy)):
            if (not done):
               break
            if (nindexb <= nindexa):
               continue
            indexb = selfcopy[nindexb]
            if (indexa.isgreaterthan(indexb)):
               swap = copy.deepcopy(indexa)
               selfcopy[nindexa] = copy.deepcopy(indexb)
               selfcopy[nindexb] = swap
               done = 0
      if (done):
         return selfcopy

def printindexes(self):
   """Prints a list of indexes"""
   show = "("
   for index in self:
      show = string.join([show,index.show()])
   show = string.join([show,")"])
   print show

def printpermutation(self):
   """Prints a permutation of indexes"""
   show = "("
   for nindex in range(len(self)/2):
      index = self[nindex]
      show = string.join([show,index.show()])
   show = string.join([show,"=>"])
   for nindex in range(len(self)/2,len(self)):
      index = self[nindex]
      show = string.join([show,index.show()])
   show = string.join([show,")"])
   print show

def expand(nestedlist):
   """Expands a nested list into a non-nested list"""
   result = []
   if (nestedlist):
      for member in nestedlist:
         if (isinstance(member,types.ListType)):
            result = result + expand(member)
         else:
            result.append(member)
      return result
   else:
      return result

def createasmallfactor(before,after,all):
   """Creates a factor of 1 + P(before -> after)"""
   beforeandafter = copy.deepcopy(all)
   parity = 0.5
   for i in range(len(all)):
      foundinbefore = -1
      foundinafter = -1
      for j in range(len(before)):
         if (before[j].isidenticalto(all[i])):
            foundinbefore = j 
         if (after[j].isidenticalto(all[i])):
            foundinafter = j 
      if (foundinbefore != -1):
         beforeandafter.append(after[foundinbefore])
         parity = - parity
      elif (foundinafter != -1):
         beforeandafter.append(before[foundinafter])
      else:
         beforeandafter.append(all[i])
   factor = Factor([0.5],[[]])
   factor.add(Factor([parity],[beforeandafter]))
   return factor

def createfactor(permutables,all):
   """Creates a factor object with all possible permutations of indexes in the given permutable list"""
   factor = Factor([0.0],[[]])
   permutations = permutationwithparity(len(permutables))
   for permutation in permutations:
      before = []
      after = []
      parity = permutation[0]
      permutation = permutation[1:len(permutation)]
      for index in all:
         before.append(index)
         replaced = 0
         for ianother in range(len(permutables)):
            another = permutables[ianother]
            if (another.isidenticalto(index)):
               replacedindex = copy.deepcopy(permutables[permutation[ianother]-1])
               replaced = 1
         if (replaced):
            after.append(replacedindex)
         else:
            after.append(index)
      beforeandafter = before + after
      factor.add(Factor([float(parity)/float(factorial(len(permutation)))],[beforeandafter]))
   return factor

def permutationsoffoursets(seedfactor,onesuper,twosuper,onesub,twosub,targetsuper,targetsub):
   """Returns a set of permutations as a factor object that permutes indices all possible meaningful ways"""
   # "meaningful" means excluding permutations of indices within given four sets.

   newfactor = seedfactor.duplicate()
   nswapsuper = min(len(onesuper),len(twosuper)) + 1
   nswapsub = min(len(onesub),len(twosub)) + 1
   for iswapsuper in range(nswapsuper):
      swaponesuperlists = picknfromlist(iswapsuper,onesuper)
      swaptwosuperlists = picknfromlist(iswapsuper,twosuper)
      for swaponesuperlist in swaponesuperlists:
         for swaptwosuperlist in swaptwosuperlists:
            for iswapsub in range(nswapsub):
               swaponesublists = picknfromlist(iswapsub,onesub)
               swaptwosublists = picknfromlist(iswapsub,twosub)
               for swaponesublist in swaponesublists:
                  for swaptwosublist in swaptwosublists:
                     permutation = []
                     permutation = permutation + copy.deepcopy(targetsuper) + copy.deepcopy(targetsub)
                     parity = 1.0
                     for targetindex in targetsuper:
                        swapped = 0
                        for nindexone in range(len(swaponesuperlist)):
                           indexone = swaponesuperlist[nindexone]
                           if (targetindex.isidenticalto(indexone)):
                              permutation.append(copy.deepcopy(swaptwosuperlist[nindexone]))
                              swapped = 1
                              parity = parity * (-1.0)
                        for nindextwo in range(len(swaponesuperlist)):
                           indextwo = swaptwosuperlist[nindextwo]
                           if (targetindex.isidenticalto(indextwo)):
                              permutation.append(copy.deepcopy(swaponesuperlist[nindextwo]))
                              swapped = 1
                        if (not swapped):
                           permutation.append(targetindex)
                     for targetindex in targetsub:
                        swapped = 0
                        for nindexone in range(len(swaponesublist)):
                           indexone = swaponesublist[nindexone]
                           if (targetindex.isidenticalto(indexone)):
                              permutation.append(copy.deepcopy(swaptwosublist[nindexone]))
                              swapped = 1
                              parity = parity * (-1.0)
                        for nindextwo in range(len(swaponesublist)):
                           indextwo = swaptwosublist[nindextwo]
                           if (targetindex.isidenticalto(indextwo)):
                              permutation.append(copy.deepcopy(swaponesublist[nindextwo]))
                              swapped = 1
                        if (not swapped):
                           permutation.append(targetindex)
                     identity = 1
                     for nindex in range(len(permutation)/2):
                        if (not permutation[nindex].isidenticalto(permutation[nindex+len(permutation)/2])):
                           identity = 0
                     if (not identity):
                        newfactor.add(Factor([parity*seedfactor.coefficients[0]],[permutation]))
   newfactor = newfactor.canonicalize(onesuper)
   newfactor = newfactor.canonicalize(onesub)
   newfactor = newfactor.canonicalize(twosuper)
   newfactor = newfactor.canonicalize(twosub)
   return newfactor
 
def writetofile(list,filename):
   """Writes a list to a given file"""
   file = open(filename,"w")
   for n in list:
      file.write(n)
      file.write("\n")

class Index:

   def __init__(self,type="unknown",label=0):
      """Creates a hole/particle/general index of tensors"""
      self.type = type
      self.label = label

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = string.join([self.type[0], repr(self.label)], "")
      return show

   def tex(self):
      """Returns a LaTeX string of the content"""
      show = string.join([self.type[0],"_{",repr(self.label),"}"], "")
      return show

   def isidenticalto(self,another):
      """Returns true if self and another indexes are identical"""
      if ((self.type == another.type) and (self.label == another.label)):
         return 1
      else:
         return 0
 
   def isin(self,list):
      """Returns true if an index is in the list"""
      for index in list:
         if (self.isidenticalto(index)):
            return 1
      return 0

   def isgeneral(self):
      """Returns true if self is a general (as opposed to particle/hole) index"""
      if (self.type == "general"):
         return 1
      else:
         return 0

   def ishole(self):
      """Returns true if self is a hole index"""
      if (self.type == "hole"):
         return 1
      else:
         return 0

   def isparticle(self):
      """Returns true if self is a particle index"""
      if (self.type == "particle"):
         return 1
      else:
         return 0
 
   def isgreaterthan(self,another):
      """Returns true if self should be to the right of another in the canonical order"""
 
      if ((self.type == 'hole') and (another.type == 'particle')):
         return 0
      elif ((self.type == 'hole') and (another.type == 'general')):
         return 0
      elif ((self.type == 'particle') and (another.type == 'hole')):
         return 1
      elif ((self.type == 'particle') and (another.type == 'general')):
         return 0
      elif ((self.type == 'general') and (another.type == 'hole')):
         return 1
      elif ((self.type == 'general') and (another.type == 'particle')):
         return 1
 
      # at this point self.type = another.type
      if (self.label > another.label):
         return 1
      else:
         return 0
 
      return 0

   def duplicate(self):
      """Returns a deepcopy of self"""
      duplicate = Index(self.type,self.label)
      return duplicate

class Factor:

   def __init__(self,coefficients=[],permutations=[]):
      """Creates a numerical and permutation factor of an operator sequence"""
      self.coefficients = coefficients
      self.permutations = copy.deepcopy(permutations)

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self,verbose=1):
      """Returns a human-friendly string of contests"""
      if (verbose):
         show = "["
         for n in range(len(self.coefficients)):
            coefficient = self.coefficients[n]
            num = rationaltofractional(coefficient)[0]
            den = rationaltofractional(coefficient)[1]
            if (num >= 0):
               show = string.join([show,"+",repr(num)])
            elif (num < 0):
               show = string.join([show,"-",repr(-num)])
            if (den != 1):
               show = string.join([show,"/",repr(den)],"")
            if (self.permutations[n]):
               show = string.join([show,"* P("])
               for nindex in range(len(self.permutations[n])/2):
                  index = self.permutations[n][nindex]
                  show = string.join([show,index.show()])
               show = string.join([show,"=>"])
               for nindex in range(len(self.permutations[n])/2,len(self.permutations[n])):
                  index = self.permutations[n][nindex]
                  show = string.join([show,index.show()])
               show = string.join([show,")"])
         show = string.join([show,"]"])
      else:
         num = rationaltofractional(self.coefficients[0])[0]
         den = rationaltofractional(self.coefficients[0])[1]
         show = repr(num)
         if (den != 1):
            show = string.join([show,"/",repr(den)],"")
         if (len(self.coefficients) > 1):
            show = string.join([show,"* P(",repr(len(self.coefficients)),")"])
      return show
 
   def tex(self,verbose=1):
      """Returns a LaTeX string of contests"""
      if (verbose):
         coefficient = self.coefficients[0]
         for n in range(len(self.coefficients)):
            if (abs(self.coefficients[n]) != abs(coefficient)):
               raise RuntimeError, "unrealistic factor"
         if (abs(coefficient) < 1.0):
            fraction = abs(int(1.0/coefficient))
            if (1.0/float(fraction) != abs(coefficient)):
               print " !!! WARNING !!! inaccurate arithmatic"
            if (fraction == 1):
               frac = ""
            else:
               frac = string.join(["\\frac{1}{",str(fraction),"}"],"")
         else:
            fraction = abs(int(coefficient))
            if (fraction == 1):
               frac = ""
            else:
               frac = str(fraction)
         if (coefficient >= 0.0):
            show = string.join(["+",frac])
         elif (coefficient < 0.0):
            show = string.join(["-",frac])
         if (len(self.coefficients) > 1):
            show = string.join([show,"\\left("],"")
            for n in range(len(self.coefficients)):
               if (self.permutations[n]):
                  if (self.coefficients[n]/coefficient > 0.0):
                     show = string.join([show,"+"],"")
                  else:
                     show = string.join([show,"-"],"")
                  show = string.join([show,"P^{"])
                  super = ""
                  sub = ""
                  for nindex in range(len(self.permutations[n])/2):
                     indexa = self.permutations[n][nindex]
                     indexb = self.permutations[n][nindex+len(self.permutations[n])/2]
                     if (indexa.type == 'hole'):
                        super = string.join([super,indexa.tex()])
                        sub = string.join([sub,indexb.tex()])
                     elif (indexa.type == 'particle'):
                        super = string.join([super,indexb.tex()])
                        sub = string.join([sub,indexa.tex()])
                  show = string.join([show,super,"}_{",sub,"}"],"")
               else:
                  if (self.coefficients[n]/coefficient < 0.0):
                     show = string.join([show,"-"],"")
                  show = string.join([show,"1"],"")
            show = string.join([show,"\\right)"])
      else:
         coefficient = self.coefficients[0]
         if (abs(coefficient) < 1.0):
            fraction = abs(int(1.0/coefficient))
            if (1.0/float(fraction) != abs(coefficient)):
               print " !!! WARNING !!! inaccurate arithmatic"
            if (fraction == 1):
               frac = ""
            else:
               frac = string.join(["\\frac{1}{",str(fraction),"}"],"")
         else:
            fraction = abs(int(coefficient))
            if (fraction == 1):
               frac = ""
            else:
               frac = str(fraction)
         if (self.coefficients[0] >= 0.0):
            show = string.join(["+",frac])
         elif (self.coefficients[0] < 0.0):
            show = string.join(["-",frac])
         if (len(self.coefficients) > 1):
            show = string.join([show,"P_{",repr(len(self.coefficients)),"}"])
      return show

   def duplicate(self):
      """Returns a deepcopy"""
      duplicate = Factor(copy.deepcopy(self.coefficients),copy.deepcopy(self.permutations))
      return duplicate

   def multiply(self,factor):
      """Multiply a factor to all coefficients"""
      for n in range(len(self.coefficients)):
         self.coefficients[n] = self.coefficients[n] * factor

   def add(self,another,factor=1.0):
      """Add two Factors together"""
      for m in range(len(another.coefficients)):
         done = 0
         for n in range(len(self.coefficients)):
            if (arethesamepermutation(self.permutations[n],another.permutations[m])):
               self.coefficients[n] = self.coefficients[n] + another.coefficients[m] * factor
               done = 1
         if (not done):
            self.coefficients.append(another.coefficients[m] * factor)
            self.permutations.append(another.permutations[m])

   def isthesameas(self,another):
      """Returns true if two factors are the same to a common scalar multiplier"""
      if (len(self.coefficients) != len(another.coefficients)):
         return 0
      alreadyused = []
      for nself in range(len(self.coefficients)):
         selfpermutation = self.permutations[nself]
         selfcoefficient = self.coefficients[nself]
         found = 0
         for nanother in range(len(another.coefficients)):
            if (nanother in alreadyused):
               continue
            anotherpermutation = another.permutations[nanother]
            anothercoefficient = another.coefficients[nanother]
            if (arethesamepermutation(selfpermutation,anotherpermutation)):
               ratio = anothercoefficient/selfcoefficient
               if ((nself != 0) and (ratio != previousratio)):
                  return 0
               previousratio = ratio
               found = 1
               alreadyused.append(nanother)
         if (not found):
            return 0
      return ratio

   def checkparity(self):
      """Check if the parity of permutations seems correct"""
      for n in range(len(self.coefficients)):
         if (self.permutations[n] == []):
            denominator = self.coefficients[n]
      for n in range(len(self.coefficients)):
         if (self.permutations[n]):
            parity = self.coefficients[n]/denominator
            permutations = permutationwithparity(len(self.permutations[n])/2)
            for permutation in permutations:
               same = 1
               for m in range(len(self.permutations[n])/2):
                  if (not self.permutations[n][m].isidenticalto( \
                          self.permutations[n][permutation[m+1]-1+len(self.permutations[n])/2])):
                     same = 0
               if (same):
                  if (parity != permutation[0]):
                     return 0
      return 1

   def canonicalize(self,permutables):
      """Reorder the permutable indexes in the ascending order"""
      # "Permutable indexes" are the target indexes appearing in one tensor
      # By canonicalizing permutation indexes, one can break down
      # permutation operator into two smaller ones plus take advantage
      # of restricted target index ranges
 
      another = self.duplicate()
      # sort the left half in the ascending order regardless of permutables
      for nfactor in range(len(another.coefficients)):
         coefficient = another.coefficients[nfactor]
         permutation = another.permutations[nfactor]
         done = 0
         while (not done):
            done = 1
            for nindexa in range(len(permutation)/2):
               for nindexb in range(len(permutation)/2):
                  if (nindexa >= nindexb):
                     continue
                  indexa = permutation[nindexa]
                  indexb = permutation[nindexb]
                  if (indexa.isgreaterthan(indexb)):
                     another.permutations[nfactor][nindexb] = copy.deepcopy(indexa)
                     another.permutations[nfactor][nindexa] = copy.deepcopy(indexb)
                     indexc = another.permutations[nfactor][nindexa+len(permutation)/2]
                     indexd = another.permutations[nfactor][nindexb+len(permutation)/2]
                     another.permutations[nfactor][nindexb+len(permutation)/2] = copy.deepcopy(indexc)
                     another.permutations[nfactor][nindexa+len(permutation)/2] = copy.deepcopy(indexd)
                     done = 0
      # sort the right half in the ascending order within the permutables
      for nfactor in range(len(another.coefficients)):
         coefficient = another.coefficients[nfactor]
         permutation = another.permutations[nfactor]
         done = 0
         while (not done):
            done = 1
            for nindexa in range(len(permutation)/2,len(permutation)):
               for nindexb in range(len(permutation)/2,len(permutation)):
                  if (nindexa >= nindexb):
                     continue
                  indexa = permutation[nindexa]
                  indexb = permutation[nindexb]
                  if (permutation[nindexa-len(permutation)/2].isin(permutables) and permutation[nindexb-len(permutation)/2].isin(permutables)):
                     if (indexa.isgreaterthan(indexb)):
                        another.permutations[nfactor][nindexb] = copy.deepcopy(indexa)
                        another.permutations[nfactor][nindexa] = copy.deepcopy(indexb)
                        another.coefficients[nfactor] = (-1.0) * another.coefficients[nfactor]
                        done = 0
      # simplify and bring to the head the identity
      for npermutation in range(len(another.permutations)):
         permutation = another.permutations[npermutation]
         identity = 1
         for nindex in range(len(permutation)/2):
            if (not permutation[nindex].isidenticalto(permutation[nindex+len(permutation)/2])):
               identity = 0
         if (identity):
            another.permutations[npermutation] = copy.deepcopy(another.permutations[0])
            another.permutations[0] = []
            swap = another.coefficients[npermutation]
            another.coefficients[npermutation] = copy.deepcopy(another.coefficients[0])
            another.coefficients[0] = swap

      return another

   def product(self,another):
      """Returns a product of two permutation operators"""
      factorobjectiscreated = 0
      for iself in range(len(self.permutations)):
         selfpermutation = self.permutations[iself]
         selfcoefficient = self.coefficients[iself]
         for ianother in range(len(another.permutations)):
            anotherpermutation = another.permutations[ianother]
            anothercoefficient = another.coefficients[ianother]
            productpermutation = combinepermutations(selfpermutation,anotherpermutation)
            productcoefficient = selfcoefficient * anothercoefficient
            if (factorobjectiscreated):
               product.add(Factor([productcoefficient],[productpermutation]))
            else:
               product = Factor([productcoefficient],[productpermutation])
               factorobjectiscreated = 1
      return product

   def normalize(self):
      """Make the absolute values of coefficients sum to unity"""
      another = self.duplicate()
      nterms = len(another.coefficients)
      if (nterms == 0):
         raise RuntimeError, "unable to normalize a factor"
      another.multiply(1.0/nterms)
      return another 

class Summation:

   def __init__(self,indexes=[]):
      """Creates a summation"""
      self.indexes = indexes

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = "Sum ("
      for index in self.indexes:
         show = string.join([show, index.show()])
      show = string.join([show,")"])
      return show

   def tex(self):
      """Returns a LaTeX string of the content"""
      show = ""
      for index in self.indexes:
         if (show):
            show = string.join([show,","],"")
         else:
            show = "\\sum_{"
         show = string.join([show,index.tex()])
      show = string.join([show,"}"])
      return show

   def hastheindex(self,another):
      """Returns true if the summation has the input index"""
      has = 0
      for index in self.indexes:
         if (index.isidenticalto(another)):
            has = 1
      return has

   def duplicate(self):
      """Returns a deepcopy of itself"""
      duplicate = Summation([])
      for index in self.indexes:
         duplicate.indexes.append(index.duplicate())
      return duplicate

   def canonicalize(self):
      """Sort the indexes in the canonical order"""
      another = self.duplicate()
      for nindexa in range(len(another.indexes)):
         indexa = another.indexes[nindexa]
         for nindexb in range(len(another.indexes)):
            indexb = another.indexes[nindexb]
            if (nindexa <= nindexb):
               continue
            if (not indexa.isgreaterthan(indexb,another)):
               swap = another.indexes[nindexa]
               another.indexes[nindexa] = copy.deepcopy(another.indexes[nindexb])
               another.indexes[nindexb] = copy.deepcopy(swap)
      return another

class Tensor:

   def __init__(self,type="unknown",indexes=[],label=0,conjugate=0,irrep=[]):
      """Creates an integral/amplitude"""
      self.type = type
      self.indexes = indexes
      self.label = label
      self.conjugate = conjugate
      if ((type == "i") or (type == "j")):
         self.irrep = irrep
      else:
         self.irrep = [self.type]

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = self.type
      if ((self.type == "i") or (self.type == "j")):
         show = string.join([show, repr(self.label)],"")
      if (self.conjugate):
         show = string.join([show, "+"],"")
      show = string.join([show, "("])
      for index in self.indexes:
         show = string.join([show, index.show()])
      show = string.join([show,")_"])
      for each in self.irrep:
         show = string.join([show, each],"")
      return show
 
   def tex(self):
      """Returns a LaTeX string of the content"""
      intermediatelist = ["chi","xi","kappa","eta","zeta","gamma","theta","lambda","pi","sigma"]
      persistentintermediatelist = ["Theta","Gamma","Lambda","Xi","Pi","Sigma","Upsilon","Omega"]
      if (self.type == "i"):
         if (self.label < 10):
            show = string.join(["\\",intermediatelist[self.label]],"")
         else:
            show = string.join(["\\",intermediatelist[0],"_{",repr(self.label),"}"],"")
      elif (self.type == "j"):
         if (self.label < 8):
            show = string.join(["\\",persistentintermediatelist[self.label]],"")
         else:
            show = string.join(["\\",persistentintermediatelist[0],"_{",repr(self.label),"}"],"")
      else:
         show = self.type
      show = string.join([show, "^{"])
      for index in self.indexes[0:len(self.indexes)/2]:
         show = string.join([show, index.tex()])
      show = string.join([show,"}"])
      show = string.join([show, "_{"])
      for index in self.indexes[len(self.indexes)/2:len(self.indexes)]:
         show = string.join([show, index.tex()])
      show = string.join([show,"}"])
      if (self.conjugate):
         show = string.join(["\\left(",show,"\\right)^{\\dagger}"],"")
      return show
 
   def textable(self,name=""):
      """Returns a LaTeX string of the content"""
      if (self.type == "i"):
         if (not name):
            show = "r"
         else:
            show = string.join(["\\left(\\xi_{",name,"}\\right)"],"")
      elif (self.type == "j"):
         show = string.join(["\\left(\\Xi_{",str(self.label),"}\\right)"],"")
      else:
         show = self.type
      show = string.join([show, "^{"])
      for index in self.indexes[0:len(self.indexes)/2]:
         show = string.join([show, index.tex()])
      show = string.join([show,"}"])
      show = string.join([show, "_{"])
      for index in self.indexes[len(self.indexes)/2:len(self.indexes)]:
         show = string.join([show, index.tex()])
      show = string.join([show,"}"])
      if (self.conjugate):
         show = string.join(["\\left(",show,"\\right)^{\\dagger}"],"")
      return show
 
   def duplicate(self):
      """Makes a copy of itself"""
      duplicate = Tensor()
      duplicate.type = self.type
      duplicate.indexes = copy.deepcopy(self.indexes)
      duplicate.label = self.label
      duplicate.conjugate = self.conjugate
      duplicate.irrep = copy.deepcopy(self.irrep)
      return duplicate
 
   def usesindexlabel(self,label):
      """Returns true if the input index label is already in use"""
      for index in self.indexes:
         if (index.label == label):
            return 1
      return 0

   def relabels(self,oldlabel,newlabel):
      """Renames an index label (followed by index sort if the tensor is intermediate)"""
      for nindex in range(len(self.indexes)):
         index = self.indexes[nindex]
         if (index.label == oldlabel):
            self.indexes[nindex].label = newlabel
      if ((self.type == "i") or (self.type == "j")):
         parity = self.sortindexes()

   def swapindexes(self,indexone,indextwo):
      """Swaps indexes"""
      for nindex in range(len(self.indexes)):
         index = self.indexes[nindex]
         if (index.isidenticalto(indexone)):
            self.indexes[nindex] = indextwo
         elif (index.isidenticalto(indextwo)):
            self.indexes[nindex] = indexone
      if ((self.type == "i") or (self.type == "j")):
         parity = self.sortindexes()

   def sortindexes(self):
      """Sort super and sub indexes of tensors in ascending order"""
      parity = 1
      alwaystrue = 1
      while (alwaystrue):
         done = 1
         for nindex in range(len(self.indexes)/2):
            indextype = self.indexes[nindex].type
            indexlabel = self.indexes[nindex].label
            for nanother in range(len(self.indexes)/2):
               if (nindex < nanother):
                  anothertype = self.indexes[nanother].type
                  anotherlabel = self.indexes[nanother].label
                  if ((indextype == 'particle') and (anothertype == 'hole')):
                     done = 0
                     swapa = nindex
                     swapb = nanother
                  elif ((indextype == anothertype) and (indexlabel >= anotherlabel)):
                     done = 0
                     swapa = nindex
                     swapb = nanother
         if (done):
            alwaystrue = 1
            while (alwaystrue):
               done = 1
               for nindex in range(len(self.indexes)/2,len(self.indexes)):
                  indextype = self.indexes[nindex].type
                  indexlabel = self.indexes[nindex].label
                  for nanother in range(len(self.indexes)/2,len(self.indexes)):
                     if (nindex < nanother):
                        anothertype = self.indexes[nanother].type
                        anotherlabel = self.indexes[nanother].label
                        if ((indextype == 'particle') and (anothertype == 'hole')):
                           done = 0
                           swapa = nindex
                           swapb = nanother
                        elif ((indextype == anothertype) and (indexlabel >= anotherlabel)):
                           done = 0
                           swapa = nindex
                           swapb = nanother
               if (done):
                  return parity
               else:
                  swap = self.indexes[swapa]
                  self.indexes[swapa] = self.indexes[swapb]
                  self.indexes[swapb] = swap
                  parity = parity * (-1)
         else:
            swap = self.indexes[swapa]
            self.indexes[swapa] = self.indexes[swapb]
            self.indexes[swapb] = swap
            parity = parity * (-1)

   def isidenticalto(self,another):
      """Returns true if two tensors are identical"""
      if (self.type != another.type):
         return 0
      if (self.label != another.label):
         return 0
      if (self.conjugate != another.conjugate):
         return 0
      if (len(self.indexes) != len(another.indexes)):
         return 0
      for nindex in range(len(self.indexes)):
         selfindex = self.indexes[nindex]
         anotherindex = another.indexes[nindex]
         if (not selfindex.isidenticalto(anotherindex)):
            return 0
      return 1

   def hasindexeswith(self,another):
      """Returns a non-redundant set of indexes that are in self and another tensors"""
      result = []
      for index in self.indexes:
         alreadyincluded = 0
         for resultindex in result:
            if (index.isidenticalto(resultindex)):
               alreadyincluded = 1
         if (not alreadyincluded):
            result.append(index)
      for index in another.indexes:
         alreadyincluded = 0
         for resultindex in result:
            if (index.isidenticalto(resultindex)):
               alreadyincluded = 1
         if (not alreadyincluded):
            result.append(index)
      return result

   def hascommonindexeswith(self,another):
      """Returns a set of indexes that are commonly in self and another tensors"""
      result = []
      for index in self.indexes:
         for anotherindex in another.indexes:
            if (index.isidenticalto(anotherindex)):
               result.append(index)
      return result

   def contracts(self,another,label):
      """Return the tensor obtained by a contraction of self and another tensors"""

      if ((not self.conjugate) and (not another.conjugate)):
         indexes = self.indexes[0:len(self.indexes)/2] + another.indexes[0:len(another.indexes)/2] \
                 + self.indexes[len(self.indexes)/2:]  + another.indexes[len(another.indexes)/2:]
      elif ((not self.conjugate) and (another.conjugate)):
         indexes = self.indexes[0:len(self.indexes)/2] + another.indexes[len(another.indexes)/2:] \
                 + self.indexes[len(self.indexes)/2:]  + another.indexes[0:len(another.indexes)/2]
      elif ((self.conjugate) and (not another.conjugate)):
         indexes = self.indexes[len(self.indexes)/2:]  + another.indexes[0:len(another.indexes)/2] \
                 + self.indexes[0:len(self.indexes)/2] + another.indexes[len(another.indexes)/2:]
      elif ((self.conjugate) and (another.conjugate)):
         indexes = self.indexes[len(self.indexes)/2:]  + another.indexes[len(another.indexes)/2:] \
                 + self.indexes[0:len(self.indexes)/2] + another.indexes[0:len(another.indexes)/2]

      product = self.irrep + another.irrep

      # Eliminate any common indexes between super/subindexes
      alwaystrue = 1
      while (alwaystrue):
         common = 0
         halflength = len(indexes)/2
         for isuper in range(0,halflength):
            superindex = indexes[isuper]
            for isub in range(halflength,2*halflength):
               subindex = indexes[isub]
               if (superindex.isidenticalto(subindex)):
                  common = 1
                  super = isuper
                  sub = isub
         if (not common):
            # reorder the super and sub indexes individually
            # (I think that the order of indexes of intermediates is immaterial, so we can decide for our convenience.
            #  However, do not mix up super and sub indexes, since otherwise contraction will be screwed up.)
            intermediate = Tensor("i",indexes,label,0,product)
            parity = intermediate.sortindexes()
            return intermediate
         else:
            if (super == sub - halflength):
               # preserve the order of deletion
               del indexes[sub]
               del indexes[super]
            else:
               indexes[super] = indexes[sub-halflength]
               # preserve the order of deletion
               del indexes[sub]
               del indexes[sub-halflength]
 
   def isequivalentto(self,another,tensorcontraction):
      """Returns true if two amplitudes are of the same type and are connected in the same way"""
 
      if (self.type != another.type):
         return 0
      if (len(self.indexes) != len(another.indexes)):
         return 0
      # super hole, super particle, super general, sub hole, sub particle, sub general
      selftargets = [0,0,0,0,0,0]
      for iselfindex in range(len(self.indexes)/2):
         selfindex = self.indexes[iselfindex]
         if (not tensorcontraction.summation.hastheindex(selfindex)):
            if (selfindex.type == 'hole'):
               selftargets[0] = selftargets[0] + 1
            elif (selfindex.type == 'particle'):
               selftargets[1] = selftargets[1] + 1
            elif (selfindex.type == 'general'):
               selftargets[2] = selftargets[2] + 1
      for iselfindex in range(len(self.indexes)/2,len(self.indexes)):
         selfindex = self.indexes[iselfindex]
         if (not tensorcontraction.summation.hastheindex(selfindex)):
            if (selfindex.type == 'hole'):
               selftargets[3] = selftargets[3] + 1
            elif (selfindex.type == 'particle'):
               selftargets[4] = selftargets[4] + 1
            elif (selfindex.type == 'general'):
               selftargets[5] = selftargets[5] + 1
      anothertargets = [0,0,0,0,0,0]
      for ianotherindex in range(len(another.indexes)/2):
         anotherindex = another.indexes[ianotherindex]
         if (not tensorcontraction.summation.hastheindex(anotherindex)):
            if (anotherindex.type == 'hole'):
               anothertargets[0] = anothertargets[0] + 1
            elif (anotherindex.type == 'particle'):
               anothertargets[1] = anothertargets[1] + 1
            elif (anotherindex.type == 'general'):
               anothertargets[2] = anothertargets[2] + 1
      for ianotherindex in range(len(another.indexes)/2,len(another.indexes)):
         anotherindex = another.indexes[ianotherindex]
         if (not tensorcontraction.summation.hastheindex(anotherindex)):
            if (anotherindex.type == 'hole'):
               anothertargets[3] = anothertargets[3] + 1
            elif (anotherindex.type == 'particle'):
               anothertargets[4] = anothertargets[4] + 1
            elif (anotherindex.type == 'general'):
               anothertargets[5] = anothertargets[5] + 1
      if (selftargets != anothertargets):
         return 0
      else:
         return 1

   def pythongen(self):
      """Returns a character expression in Python syntax"""
      pythoncode = self.type
      if (self.type == "i"):
         pythoncode = string.join([pythoncode,repr(self.label)],"")
      else:
         pythoncode = string.join([pythoncode,repr(len(self.indexes)/2)],"")
      pythoncode = string.join([pythoncode,"[","("*len(self.indexes)],"")
      if (self.indexes):
         for nindex in range(len(self.indexes)):
            index = self.indexes[nindex]
            if (nindex == 0):
               pythoncode = string.join([pythoncode,index.show(),")"],"")
            else:
               pythoncode = string.join([pythoncode,"*N+",index.show(),")"],"")
      else:
         pythoncode = string.join([pythoncode,"0"],"")
      pythoncode = string.join([pythoncode,"]"],"")
      return pythoncode

   def fortran90(self,types,permutation=[],reverse=0,suffix=""):
      """Returns a character expression in Fortran90 syntax"""
      f90code = self.type
      if (self.type in types[3]):
         f90code = string.join([f90code,repr(self.label),suffix],"")
      else:
         f90code = string.join([f90code,repr(len(self.indexes)/2),suffix],"")
      f90code = string.join([f90code,"(","("*len(self.indexes)],"")
      if (self.indexes):
         for nindex in range(len(self.indexes)):
            if (permutation == []):
               index = self.indexes[nindex]
            else:
               if (reverse):
                  beforeindex = self.indexes[nindex]
                  index = beforeindex
                  for nbeforeindex in range(len(permutation)/2):
                     if (beforeindex.isidenticalto(permutation[nbeforeindex + len(permutation)/2])):
                        index = permutation[nbeforeindex]
               else:
                  beforeindex = self.indexes[nindex]
                  index = beforeindex
                  for nbeforeindex in range(len(permutation)/2):
                     if (beforeindex.isidenticalto(permutation[nbeforeindex])):
                        index = permutation[nbeforeindex+len(permutation)/2]
            if (nindex == 0):
               f90code = string.join([f90code,index.show(),"-1)"],"")
            elif (nindex == len(self.indexes) - 1):
               f90code = string.join([f90code,"*N+",index.show(),")"],"")
            else:
               f90code = string.join([f90code,"*N+",index.show(),"-1)"],"")
      else:
         f90code = string.join([f90code,"1"],"")
      f90code = string.join([f90code,")"],"")
      return f90code

   def fortran90x(self,types):
      """Anti-symmetrizer"""

      newcode = Code("Fortran90","")
      
      if (not self.indexes):
         return newcode

      # generate loops over indexes
      for index in self.indexes:
         newcode.insertdoloop(index)
      
      # generate anti-symmetrizer statement
      
      nindexes = len(self.indexes)/2
      npermutations = factorial(nindexes)
      permutations = permutationwithparity(nindexes)
      newline = string.join(["TMP = 0.0d0"])
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      for terma in range(npermutations):
         for termb in range(npermutations):
            newindexes = []
            for nindexa in range(1,nindexes+1):
               newindexes.append(self.indexes[permutations[terma][nindexa]-1])
            for nindexb in range(1,nindexes+1):
               newindexes.append(self.indexes[nindexes+permutations[termb][nindexb]-1])
            newtensor = Tensor(self.type,newindexes,self.label,self.irrep)
            if (permutations[terma][0]*permutations[termb][0] == 1):
               parity = "+"
            else:
               parity = "-"
            newline = string.join(["TMP = TMP",parity,newtensor.fortran90(types)])
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
      newline = string.join([self.fortran90(types)," = TMP / ",repr(npermutations**2),".0d0"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = len(newcode.statements)

      return newcode

   def fortran77x(self,subroutinename="NONAME"):
      """Anti-symmetrizer"""

      errquit = 0

      newcode = Code("Fortran77","ANTISYM_"+subroutinename)
 
      # header
      newline = "!" + self.show()
      newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
 
      # insert include statements
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)
      newline = '#include "sym.fh"'
      newcode.add("headers",newline)
      newline = '#include "errquit.fh"'
      newcode.add("headers",newline)
      newline = '#include "tce.fh"'
      newcode.add("headers",newline)
 
      # declaration
      newint = "d_a"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "k_a_offset"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      
      # parallel declaration
      newint = "NXTASK"
      newcode.add("integers",newint)
      newcode.add("externals",newint)
      newint = "next"
      newcode.add("integers",newint)
      newint = "nprocs"
      newcode.add("integers",newint)
      newint = "count"
      newcode.add("integers",newint)

      super = self.indexes[0:len(self.indexes)/2]
      sub = self.indexes[len(self.indexes)/2:len(self.indexes)]
      all = super + sub
 
      # parallel related
      newline = "nprocs = GA_NNODES()"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "count = 0"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "next = NXTASK(nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # loop over output tensor indexes
      newcode.inserttileddoloops(sub[0:1])
      newcode.inserttileddoloops(super[0:1])
      for index in super:
         if (not index.isidenticalto(super[0])):
            newint = index.show()+"b"
            newcode.add("integers",newint)
            newline = string.join([newint," = ",super[0].show(),"b"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
      for index in sub:
         if (not index.isidenticalto(sub[0])):
            newint = index.show()+"b"
            newcode.add("integers",newint)
            newline = string.join([newint," = ",sub[0].show(),"b"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
      
      # parallel related
      newline = "IF (next.eq.count) THEN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.setamark(1)
      newcode.pointer = newcode.pointer + 1
      newline = "next = NXTASK(nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END IF"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "count = count + 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.getamark(1) + 1
 
      # symmetry of output tensor
      newcode.inserttiledifsymmetry(super,sub,self.irrep)

      # allocate for original and antisymmetrized tensors
      newcode.add("integers","dim0")
      newline = ""
      for index in all:
         if (newline == ""):
            newline = string.join(["dim0 = int_mb(k_range+",index.show(),"b-1)"],"")
         else:
            newline = string.join([newline," * int_mb(k_range+",index.show(),"b-1)"],"")
      if (newline == ""):
         newline = "dim0 = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.add("integers","l_b")
      newcode.add("integers","k_b")
      newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dim0,'noname',l_b,k_b)) CALL ERRQUIT('",\
                             subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1
      newcode.add("integers","l_a")
      newcode.add("integers","k_a")
      newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dim0,'noname',l_a,k_a)) CALL ERRQUIT('",\
                             subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1

      # get a diagonal block
      newline = "IF (dim0 .gt. 0) THEN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END IF"
      newcode.statements.insert(newcode.pointer,newline)
      arguments = ""
      argumentsend = ""
      for nindex in range(len(all)-1,-1,-1):
         if (all[nindex].type == "hole"):
            boffset = "b - 1"
         else:
            boffset = "b - noab - 1"
         if (arguments == ""):
            arguments = string.join(["d_a,dbl_mb(k_a),dim0,int_mb(k_a_offset + ",all[nindex].show(),boffset],"")
         else:
            if (all[nindex+1].type == "hole"):
               arguments = string.join([arguments," + noab * (",all[nindex].show(),boffset],"")
            else:
               arguments = string.join([arguments," + nvab * (",all[nindex].show(),boffset],"")
            argumentsend = string.join([argumentsend,")"],"")
      arguments = string.join([arguments,argumentsend,")"],"")
      newline = string.join(["CALL GET_BLOCK(",arguments,")"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # do loops over indexes in the diagonal block
      for nindex in range(len(all)):
         index = all[nindex]
         newint = index.show()
         newcode.add("integers",newint)
         newline = string.join(["DO ",newint," = 1,int_mb(k_range+",index.show(),"b-1)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newline = "END DO"
         newcode.statements.insert(newcode.pointer,newline)
         if (nindex == 0):
            newcode.setamark(2)

      # permutations
      newline = ""
      newlineend = ""
      newcode.add("integers","idim0")
      for nindex in range(len(all)-1,-1,-1):
         if (newline == ""):
            newline = string.join(["idim0 = ",all[nindex].show()],"")
         else:
            newline = string.join([newline," + int_mb(k_range+",all[nindex+1].show(),"b-1)"\
                                   " * ((",all[nindex].show()," - 1)"],"")
            newlineend = string.join([newlineend,")"],"")
      newline = string.join([newline,newlineend],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      oldline = ""
      integerlabel = 0
      # the following does not mean this is for super alone, len(super)
      # is used just to get half the number of indexes
      permutations = permutationwithparity(len(super))
      for permutationa in permutations:
         for permutationb in permutations:
            integerlabel = integerlabel + 1
            newint = "idim"+repr(integerlabel)
            newcode.add("integers",newint)
            newindexes = []
            for nindexa in range(1,len(super)+1):
               newindexes.append(self.indexes[permutationa[nindexa]-1])
            for nindexb in range(1,len(super)+1):
               newindexes.append(self.indexes[len(super)+permutationb[nindexb]-1])
            newline = ""
            newlineend = ""
            for nindex in range(len(all)-1,-1,-1):
               if (newline == ""):
                  newline = string.join([newint," = ",newindexes[nindex].show()],"")
               else:
                  newline = string.join([newline," + int_mb(k_range+",newindexes[nindex+1].show(),"b-1)"\
                                         " * ((",newindexes[nindex].show()," - 1)"],"")
                  newlineend = string.join([newlineend,")"],"")
            newline = string.join([newline,newlineend],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            if (permutationa[0] * permutationb[0] == 1):
               if (oldline):
                  oldline = string.join([oldline," + dbl_mb(k_a + ",newint," - 1)"],"")
               else:
                  oldline = string.join(["dbl_mb(k_b + idim0 - 1) = (dbl_mb(k_a + ",newint," - 1)"],"")
            else:
               if (oldline):
                  oldline = string.join([oldline," - dbl_mb(k_a + ",newint," - 1)"],"")
               else:
                  oldline = string.join(["dbl_mb(k_b + idim0 - 1) = (-dbl_mb(k_a + ",newint," - 1)"],"")
      factor = float(factorial(len(super))*factorial(len(sub)))
      newline = string.join([oldline,") * ",repr(1.0/factor),"d0"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      newcode.pointer = newcode.getamark(2) + 1
      arguments = ""
      argumentsend = ""
      for nindex in range(len(all)-1,-1,-1):
         if (all[nindex].type == "hole"):
            boffset = "b - 1"
         else:
            boffset = "b - noab - 1"
         if (arguments == ""):
            arguments = string.join(["d_a,dbl_mb(k_b),dim0,int_mb(k_a_offset + ",all[nindex].show(),boffset],"")
         else:
            if (all[nindex+1].type == "hole"):
               arguments = string.join([arguments," + noab * (",all[nindex].show(),boffset],"")
            else:
               arguments = string.join([arguments," + nvab * (",all[nindex].show(),boffset],"")
            argumentsend = string.join([argumentsend,")"],"")
      arguments = string.join([arguments,argumentsend,")"],"")
      newline = string.join(["CALL PUT_BLOCK(",arguments,")"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = string.join(["IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1
      newline = string.join(["IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1

      newcode.pointer = len(newcode.statements)
      newline = "next = NXTASK(-nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "call GA_SYNC()"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      return newcode

   def fortran77y(self,globaltargetindexes,types,subroutinename="NONAME"):
      """Precompute offsets and size"""

      errquit = 0

      newcode = Code("Fortran77","OFFSET_"+subroutinename)
 
      # header
      newline = "!" + self.show()
      newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
 
      # insert include statements
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)
      newline = '#include "sym.fh"'
      newcode.add("headers",newline)
      newline = '#include "errquit.fh"'
      newcode.add("headers",newline)
      newline = '#include "tce.fh"'
      newcode.add("headers",newline)
 
      # declaration
      newint = "d_a"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "l_a_offset"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "k_a_offset"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "size"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)

      # allocate offsets
      arguments = ""
      for index in self.indexes:
         if (index.type == "hole"):
            factor = "noab"
         else:
            factor = "nvab"
         if (arguments):
            arguments = string.join([arguments,"*",factor],"")
         else:
            arguments = factor
      if (arguments == ""):
         arguments = "1"
      newline = string.join(["IF (.not.MA_PUSH_GET(mt_int,",arguments,",'noname',l_a_offset,k_a_offset)) CALL ERRQUIT('",\
                             subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1

      # classify indexes
      globalsuper = []
      localsuper = []
      globalsub = []
      localsub = []
      super = []
      sub = []
      for nindex in range(len(self.indexes)/2):
         index = self.indexes[nindex]
         super.append(index)
         if (index.isin(globaltargetindexes)):
            globalsuper.append(index)
         else:
            localsuper.append(index)
      for nindex in range(len(self.indexes)/2,len(self.indexes)):
         index = self.indexes[nindex]
         sub.append(index)
         if (index.isin(globaltargetindexes)):
            globalsub.append(index)
         else:
            localsub.append(index)

      # do loops and if's
      newline = "size = 0"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.inserttileddoloops(globalsuper)
      newcode.inserttileddoloops(localsuper)
      newcode.inserttileddoloops(globalsub)
      newcode.inserttileddoloops(localsub)
      newcode.inserttiledifsymmetry(super,sub,self.irrep)
      newcode.inserttiledifrestricted(super+sub)

      # offsets and size (offset first!)
      newline = ""
      newlineend = ""
      all = globalsuper + localsuper + globalsub + localsub
      for nindex in range(len(all)-1,-1,-1):
         if (all[nindex].type == "hole"):
            boffset = "b - 1"
         else:
            boffset = "b - noab - 1"
         if (newline == ""):
            newline = string.join(["int_mb(k_a_offset + ",all[nindex].show(),boffset],"")
         else:
            if (all[nindex+1].type == "hole"):
               newline = string.join([newline," + noab * (",all[nindex].show(),boffset],"")
            else:
               newline = string.join([newline," + nvab * (",all[nindex].show(),boffset],"")
            newlineend = string.join([newlineend,")"],"")
      if (newline == ""):
         newline = "int_mb(k_a_offset"
      newline = string.join([newline,newlineend,") = size"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = ""
      for index in all:
         if (newline == ""):
            newline = string.join(["size = size + int_mb(k_range+",index.show(),"b-1)"],"")
         else:
            newline = string.join([newline," * int_mb(k_range+",index.show(),"b-1)"],"")
      if (newline == ""):
         newline = "size = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      
      newcode.pointer = len(newcode.statements)
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      return newcode

   def utchemy(self,globaltargetindexes,subroutinename="NONAME",relativistic=0):
      """Precompute offsets and size"""

      errquit = 0

      if (relativistic):
         newcode = Code("Fortran90","R4D_OFFSET_"+subroutinename)
      else:
         newcode = Code("Fortran90","OFFSET_"+subroutinename)
 
      # header
      newline = "!" + self.show()
      newcode.add("headers",newline)
      newline = "USE UT_SYS_MODULE"
      newcode.add("headers",newline)
      newline = "USE UT_MOLINP_MODULE"
      newcode.add("headers",newline)
      if (relativistic):
         newline = "USE UT_R4DTCE_MODULE"
      else:
         newline = "USE UT_TCE_MODULE"
      newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
 
      # insert include statements
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)
 
      # declaration
      newint = "d_a"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "a_offset"
      newcode.add("integerarrays",newint)
      newcode.add("arguments",newint)
      newint = "size"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)

      # classify indexes
      globalsuper = []
      localsuper = []
      globalsub = []
      localsub = []
      super = []
      sub = []
      for nindex in range(len(self.indexes)/2):
         index = self.indexes[nindex]
         super.append(index)
         if (index.isin(globaltargetindexes)):
            globalsuper.append(index)
         else:
            localsuper.append(index)
      for nindex in range(len(self.indexes)/2,len(self.indexes)):
         index = self.indexes[nindex]
         sub.append(index)
         if (index.isin(globaltargetindexes)):
            globalsub.append(index)
         else:
            localsub.append(index)

      # do loops and if's
      newline = "size = 0"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.inserttileddoloops(globalsuper)
      newcode.inserttileddoloops(localsuper)
      newcode.inserttileddoloops(globalsub)
      newcode.inserttileddoloops(localsub)
      newcode.inserttiledifsymmetry(super,sub,self.irrep,relativistic)
      if (not relativistic):
         newcode.inserttiledifrestricted(super+sub)

      # offsets and size (offset first!)
      newline = ""
      newlineend = ""
      all = globalsuper + localsuper + globalsub + localsub
      for nindex in range(len(all)-1,-1,-1):
         if (all[nindex].type == "hole"):
            boffset = "b - 1"
            bfirstoffset = "b"
         else:
            boffset = "b - noab - 1"
            bfirstoffset = "b - noab"
         if (newline == ""):
            newline = string.join(["a_offset(",all[nindex].show(),bfirstoffset],"")
         else:
            if (all[nindex+1].type == "hole"):
               newline = string.join([newline," + noab * (",all[nindex].show(),boffset],"")
            else:
               newline = string.join([newline," + nvab * (",all[nindex].show(),boffset],"")
            newlineend = string.join([newlineend,")"],"")
      if (newline == ""):
         newline = "int_mb(k_a_offset"
      newline = string.join([newline,newlineend,") = size"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = ""
      for index in all:
         if (newline == ""):
            newline = string.join(["size = size + range(",index.show(),"b)"],"")
         else:
            newline = string.join([newline," * range(",index.show(),"b)"],"")
      if (newline == ""):
         newline = "size = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      
      newcode.pointer = len(newcode.statements)
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      return newcode

class TensorContraction:

   def __init__(self,factor=[],summation=[],tensors=[]):
      """Creates a tensor contraction"""
      self.factor = factor
      self.summation = summation
      self.tensors = tensors

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self,verbose=1):
      """Returns a human-friendly string of the content"""
      show = self.factor.show(verbose)
      if (self.summation):
         if (len(self.summation.indexes) > 0):
            show = string.join([show, "*", self.summation.show()])
      for selftensor in self.tensors:
         show = string.join([show, "*", selftensor.show()])
      return show 

   def tex(self,verbose=1):
      """Returns a LaTeX string of the content"""
      show = self.factor.tex(verbose)
#     if (self.summation):
#        if (len(self.summation.indexes) > 0):
#           show = string.join([show, self.summation.tex()])
      for selftensor in self.tensors:
         show = string.join([show, selftensor.tex()])
      return show 

   def duplicate(self):
      """Returns a deepcopy"""
      duplicate = TensorContraction(copy.deepcopy(self.factor),copy.deepcopy(self.summation),copy.deepcopy(self.tensors))
      return duplicate

   def findthebestbreakdown(self,excitation=[],verbose=0):
      """Returns the best breakdown in terms of flop costs and then memory costs"""
      ntensors = len(self.tensors)
      nbreakdowns = factorial(ntensors)
      breakdowns = permutation(ntensors)
      print " ... there are %d ways of breaking down the sequence into elementary tensor contractions" % (nbreakdowns)
      if (ntensors == 1):
         ibest = 0
         minoperationcost = [0,0,0]
         minmemorycost = [0,0,0]
      else:
         minoperationcost = [99999,99999,99999]
         minmemorycost = [99999,99999,99999]
         minaggoperationcost = [99999,99999,99999]
         minaggmemorycost = [99999,99999,99999]
         ibest = 0
         discouraged = 0
         for itensor in range(ntensors):
            if ((itensor != 0) and (itensor != ntensors-1) and \
                (self.tensors[breakdowns[ibest][itensor]-1].type not in excitation)):
               discouraged = 1
#        if ((len(breakdowns[ibest]) >= 3) and (self.tensors[breakdowns[ibest][0]-1].type in excitation) \
#                                          and (self.tensors[breakdowns[ibest][1]-1].type in excitation)):
#           discouraged = 1
#        if ((len(breakdowns[ibest]) >= 3) and (self.tensors[breakdowns[ibest][0]-1].type not in excitation) \
#                                          and (self.tensors[breakdowns[ibest][1]-1].type not in excitation)):
#           discouraged = 1
         for ibreakdown in range(0,nbreakdowns):
            newdiscouraged = 0
#           if ((len(breakdowns[ibreakdown]) >= 3) and (self.tensors[breakdowns[ibreakdown][0]-1].type in excitation) \
#                                             and (self.tensors[breakdowns[ibreakdown][1]-1].type in excitation)):
#              newdiscouraged = 1
#           if ((len(breakdowns[ibreakdown]) >= 3) and (self.tensors[breakdowns[ibreakdown][0]-1].type not in excitation) \
#                                             and (self.tensors[breakdowns[ibreakdown][1]-1].type not in excitation)):
#              newdiscouraged = 1
            for itensor in range(ntensors):
               if ((itensor != 0) and (itensor != ntensors-1) and \
                   (self.tensors[breakdowns[ibreakdown][itensor]-1].type not in excitation)):
                  newdiscouraged = 1
            tensorone = self.tensors[breakdowns[ibreakdown][0]-1]
#           if ((tensorone.type != 'f') and (tensorone.type != 'v')):
#              continue
#
# 6/18/03 The above two lines are commented out for EOM-CC left hand side,
#         wherein (L*T)*V can be the cheapest by many orders.  I hope that
#         the above is still a natural outcome for CC.
#
#           if (tensorone.type in excitation):
#              continue
#
# 6/24/03 In the above, we restrict the first tensor to be
#         either an integral tensor (V or F) or a deexcitation tensor (L).
#
            if (verbose):
               stringbreakdowna = ""
               stringbreakdownb = ""
               for itensor in range(ntensors): 
                  if (itensor == 0):
                     stringbreakdownb = stringbreakdownb + self.tensors[breakdowns[ibreakdown][itensor]-1].type \
                                      + repr(len(self.tensors[breakdowns[ibreakdown][itensor]-1].indexes)/2)
                  elif (itensor == 1):
                     stringbreakdownb = stringbreakdownb + "*" + self.tensors[breakdowns[ibreakdown][itensor]-1].type \
                                      + repr(len(self.tensors[breakdowns[ibreakdown][itensor]-1].indexes)/2)
                  else:
                     stringbreakdowna = stringbreakdowna + "("
                     stringbreakdownb = stringbreakdownb + ")*" + self.tensors[breakdowns[ibreakdown][itensor]-1].type \
                                      + repr(len(self.tensors[breakdowns[ibreakdown][itensor]-1].indexes)/2)
               print " ... ", ibreakdown + 1, breakdowns[ibreakdown], stringbreakdowna + stringbreakdownb
            # The triplets are the exponents of N, O, V
            maxoperationcost = [0,0,0]
            maxmemorycost = [0,0,0]
            aggoperationcost = [0,0,0]
            aggmemorycost = [0,0,0]
            for ielementary in range(0,ntensors-1):
               tensortwo = self.tensors[breakdowns[ibreakdown][ielementary+1]-1]
               allindexes = tensorone.hasindexeswith(tensortwo)
               commonindexes = tensorone.hascommonindexeswith(tensortwo)
               operationcost = [0,0,0]
               memorycost = [0,0,0]
               for index in allindexes:
                  if (index.isgeneral()):
                     operationcost[0] = operationcost[0] + 1
                     memorycost[0] = memorycost[0] + 1
                  elif (index.ishole()):
                     operationcost[1] = operationcost[1] + 1
                     memorycost[1] = memorycost[1] + 1
                  elif (index.isparticle()):
                     operationcost[2] = operationcost[2] + 1
                     memorycost[2] = memorycost[2] + 1
               for index in commonindexes:
                  if (index.isgeneral()):
                     memorycost[0] = memorycost[0] - 1
                  elif (index.ishole()):
                     memorycost[1] = memorycost[1] - 1
                  elif (index.isparticle()):
                     memorycost[2] = memorycost[2] - 1
               aggoperationcost[0] = aggoperationcost[0] + operationcost[0]
               aggoperationcost[1] = aggoperationcost[1] + operationcost[1]
               aggoperationcost[2] = aggoperationcost[2] + operationcost[2]
               aggmemorycost[0] = aggmemorycost[0] + memorycost[0]
               aggmemorycost[1] = aggmemorycost[1] + memorycost[1]
               aggmemorycost[2] = aggmemorycost[2] + memorycost[2]
               # if (verbose):
               #    print " ...... ",tensorone.show(), tensortwo.show(), operationcost, memorycost
               tensorone = tensorone.contracts(tensortwo,0)
               if (operationcost[0]+operationcost[1]+operationcost[2] > maxoperationcost[0]+maxoperationcost[1]+maxoperationcost[2]):
                  maxoperationcost = operationcost
               elif (operationcost[0]+operationcost[1]+operationcost[2] == maxoperationcost[0]+maxoperationcost[1]+maxoperationcost[2]):
                  if (operationcost[0] > maxoperationcost[0]):
                     maxoperationcost = operationcost
                  elif (operationcost[0] == maxoperationcost[0]):
                     if (operationcost[2] > maxoperationcost[2]):
                        maxoperationcost = operationcost
               if (memorycost[0]+memorycost[1]+memorycost[2] > maxmemorycost[0]+maxmemorycost[1]+maxmemorycost[2]):
                  maxmemorycost = memorycost
               elif (memorycost[0]+memorycost[1]+memorycost[2] == maxmemorycost[0]+maxmemorycost[1]+maxmemorycost[2]):
                  if (memorycost[0] > maxmemorycost[0]):
                     maxmemorycost = memorycost
                  elif (memorycost[0] == maxmemorycost[0]):
                     if (memorycost[2] > maxmemorycost[2]):
                        maxmemorycost = memorycost
            if (verbose):
               if (newdiscouraged):
                  print " ...... this breakdown: ", maxoperationcost, maxmemorycost, aggoperationcost, aggmemorycost," (discouraged)"
               else:
                  print " ...... this breakdown: ", maxoperationcost, maxmemorycost, aggoperationcost, aggmemorycost
            if (maxoperationcost[0]+maxoperationcost[1]+maxoperationcost[2] < minoperationcost[0]+minoperationcost[1]+minoperationcost[2]):
               minoperationcost = maxoperationcost
               minmemorycost = maxmemorycost
               minaggoperationcost = aggoperationcost
               minaggmemorycost = aggmemorycost
               ibest = ibreakdown
               discouraged = newdiscouraged
            elif (maxoperationcost[0]+maxoperationcost[1]+maxoperationcost[2] == minoperationcost[0]+minoperationcost[1]+minoperationcost[2]):
               if (newdiscouraged < discouraged):
                  minoperationcost = maxoperationcost
                  minmemorycost = maxmemorycost
                  minaggoperationcost = aggoperationcost
                  minaggmemorycost = aggmemorycost
                  ibest = ibreakdown
                  discouraged = newdiscouraged
               elif (newdiscouraged == discouraged):
                  if (maxoperationcost[0] < minoperationcost[0]):
                     minoperationcost = maxoperationcost
                     minmemorycost = maxmemorycost
                     minaggoperationcost = aggoperationcost
                     minaggmemorycost = aggmemorycost
                     ibest = ibreakdown
                     discouraged = newdiscouraged
                  elif (maxoperationcost[0] == minoperationcost[0]):
                     if (maxoperationcost[2] < minoperationcost[2]):
                        minoperationcost = maxoperationcost
                        minmemorycost = maxmemorycost
                        minaggoperationcost = aggoperationcost
                        minaggmemorycost = aggmemorycost
                        ibest = ibreakdown
                        discouraged = newdiscouraged
                     elif (maxoperationcost[2] == minoperationcost[2]):
                        if (maxmemorycost[0]+maxmemorycost[1]+maxmemorycost[2] < minmemorycost[0]+minmemorycost[1]+minmemorycost[2]):
                           minoperationcost = maxoperationcost
                           minmemorycost = maxmemorycost
                           minaggoperationcost = aggoperationcost
                           minaggmemorycost = aggmemorycost
                           ibest = ibreakdown
                           discouraged = newdiscouraged
                        elif (maxmemorycost[0]+maxmemorycost[1]+maxmemorycost[2] == minmemorycost[0]+minmemorycost[1]+minmemorycost[2]):
                           if (maxmemorycost[0] < minmemorycost[0]):
                              minoperationcost = maxoperationcost
                              minmemorycost = maxmemorycost
                              minaggoperationcost = aggoperationcost
                              minaggmemorycost = aggmemorycost
                              ibest = ibreakdown
                              discouraged = newdiscouraged
                           elif (maxmemorycost[0] == minmemorycost[0]):
                              if (maxmemorycost[2] < minmemorycost[2]):
                                 minoperationcost = maxoperationcost
                                 minmemorycost = maxmemorycost
                                 minaggoperationcost = aggoperationcost
                                 minaggmemorycost = aggmemorycost
                                 ibest = ibreakdown
                                 discouraged = newdiscouraged
                              elif (aggoperationcost[0]+aggoperationcost[1]+aggoperationcost[2] < minaggoperationcost[0]+minaggoperationcost[1]+minaggoperationcost[2]):
                                 minaggoperationcost = aggoperationcost
                                 minaggmemorycost = aggmemorycost
                                 ibest = ibreakdown
                                 discouraged = newdiscouraged
                              elif (aggoperationcost[0]+aggoperationcost[1]+aggoperationcost[2] == minaggoperationcost[0]+minaggoperationcost[1]+minaggoperationcost[2]):
                                 if (aggoperationcost[0] < minaggoperationcost[0]):
                                    minaggoperationcost = aggoperationcost
                                    minaggmemorycost = aggmemorycost
                                    ibest = ibreakdown
                                    discouraged = newdiscouraged
                                 elif (aggoperationcost[0] == minaggoperationcost[0]):
                                    if (aggoperationcost[2] < minaggoperationcost[2]):
                                       minaggoperationcost = aggoperationcost
                                       minaggmemorycost = aggmemorycost
                                       ibest = ibreakdown
                                       discouraged = newdiscouraged
                                    elif (aggoperationcost[2] == minaggoperationcost[2]):
                                       if (aggmemorycost[0]+aggmemorycost[1]+aggmemorycost[2] < minaggmemorycost[0]+minaggmemorycost[1]+minaggmemorycost[2]):
                                          minaggoperationcost = aggoperationcost
                                          minaggmemorycost = aggmemorycost
                                          ibest = ibreakdown
                                          discouraged = newdiscouraged
                                       elif (aggmemorycost[0]+aggmemorycost[1]+aggmemorycost[2] == minaggmemorycost[0]+minaggmemorycost[1]+minaggmemorycost[2]):
                                          if (aggmemorycost[0] < minaggmemorycost[0]):
                                             minaggoperationcost = aggoperationcost
                                             minaggmemorycost = aggmemorycost
                                             ibest = ibreakdown
                                             discouraged = newdiscouraged
                                          elif (aggmemorycost[0] == minaggmemorycost[0]):
                                             if (aggmemorycost[2] < minaggmemorycost[2]):
                                                minaggoperationcost = aggoperationcost
                                                minaggmemorycost = aggmemorycost
                                                ibest = ibreakdown
                                                discouraged = newdiscouraged
            if (verbose):
               print " ...... current best:   ", minoperationcost, minmemorycost, minaggoperationcost, minaggmemorycost

      print " ... the best breakdown is %s with operationcost=N%d O%d V%d, memorycost=N%d O%d V%d " % (breakdowns[ibest], \
      minoperationcost[0], minoperationcost[1], minoperationcost[2], minmemorycost[0], minmemorycost[1], minmemorycost[2])

      if ((len(breakdowns[ibest]) >= 3) and (self.tensors[breakdowns[ibest][0]-1].type in excitation) \
                                        and (self.tensors[breakdowns[ibest][1]-1].type in excitation)):
         raise RuntimeError, "null contraction appeared"

      if ((len(breakdowns[ibest]) >= 3) and (self.tensors[breakdowns[ibest][0]-1].type not in excitation) \
                                        and (self.tensors[breakdowns[ibest][1]-1].type not in excitation)):
         raise RuntimeError, "illegal contraction order has been selected"

      return breakdowns[ibest]

   def breakdown(self,label=-1,excitation=[],verbose=0):
      """Breaks down the tensor contraction into elementary tensor contractions according to a given order"""
      another = self.canonicalize(1)
      another = another.symmetrize(0)
      order = another.findthebestbreakdown(excitation,verbose)
      if (label < 0):
         if (len(order) == 1):
            label = 0
         else:
            label = len(order)-2
      result = OperationTree()
      tensorone = another.tensors[order[0]-1]
      targettensor = another.tensors[0]
      if (len(another.tensors) > 1):
         for i in range(len(another.tensors)-1):
            targettensor = targettensor.contracts(another.tensors[i+1],0)
      targetsuper = copy.deepcopy(targettensor.indexes[0:len(targettensor.indexes)/2])
      targetsub   = copy.deepcopy(targettensor.indexes[len(targettensor.indexes)/2:len(targettensor.indexes)])
      if (len(order) == 1):
         tensorthree = tensorone.duplicate()
         tensorthree.label = label
         tensorthree.type = "i"
         # reorder the super and sub indexes individually of tensorthree
         # (I think that the order of indexes of intermediates is immaterial, so we can decide for our convenience.
         # However, do not mix up super and sub indexes, since otherwise contraction will be screwed up.)
         alwaystrue = 1
         while (alwaystrue):
            done = 1
            for nindex in range(len(tensorthree.indexes)/2):
               indexlabel = tensorthree.indexes[nindex].label
               for nanother in range(len(tensorthree.indexes)/2):
                  if (nindex < nanother):
                     anotherlabel = tensorthree.indexes[nanother].label
                     if (indexlabel >= anotherlabel):
                        done = 0
                        swapa = nindex
                        swapb = nanother
            if (done):
               alwaystrue = 1
               while (alwaystrue):
                  done = 1
                  for nindex in range(len(tensorthree.indexes)/2,len(tensorthree.indexes)):
                     indexlabel = tensorthree.indexes[nindex].label
                     for nanother in range(len(tensorthree.indexes)/2,len(tensorthree.indexes)):
                        if (nindex < nanother):
                           anotherlabel = tensorthree.indexes[nanother].label
                           if (indexlabel >= anotherlabel):
                              done = 0
                              swapa = nindex
                              swapb = nanother
                  if (done):
                     elementary = ElementaryTensorContraction(another.factor,another.summation,[tensorthree,tensorone])
                     result = OperationTree(elementary,[],[result])
                     return result
                  else:
                     swap = tensorthree.indexes[swapa]
                     tensorthree.indexes[swapa] = tensorthree.indexes[swapb]
                     tensorthree.indexes[swapb] = swap
            else:
               swap = tensorthree.indexes[swapa]
               tensorthree.indexes[swapa] = tensorthree.indexes[swapb]
               tensorthree.indexes[swapb] = swap
      else:

         # breakdown of a permutation operator
         suggestedfactors = []
         factorproduct = Factor([1.0],[[]])
         for contraction in range(0,len(order)-1):
            tensortwo = another.tensors[order[contraction+1]-1]
            onesuper   = []
            onesub     = []
            twosuper   = []
            twosub     = []
            for index in tensorone.indexes:
               if (index.isin(targetsuper)):
                  onesuper.append(index)
            for index in tensorone.indexes:
               if (index.isin(targetsub)):
                  onesub.append(index)
            for index in tensortwo.indexes:
               if (index.isin(targetsuper)):
                  twosuper.append(index)
            for index in tensortwo.indexes:
               if (index.isin(targetsub)):
                  twosub.append(index)
            nswapsuper = min(len(onesuper),len(twosuper)) + 1
            nswapsub = min(len(onesub),len(twosub)) + 1
            if (contraction == len(order)-2):
               newfactor = Factor([another.factor.coefficients[0]],[[]])
            else:
               newfactor = Factor([1.0],[[]])
            newfactor = permutationsoffoursets(newfactor,onesuper,twosuper,onesub,twosub,targetsuper,targetsub)
#           for iswapsuper in range(nswapsuper):
#              swaponesuperlists = picknfromlist(iswapsuper,onesuper)
#              swaptwosuperlists = picknfromlist(iswapsuper,twosuper)
#              for swaponesuperlist in swaponesuperlists:
#                 for swaptwosuperlist in swaptwosuperlists:
#                    for iswapsub in range(nswapsub):
#                       swaponesublists = picknfromlist(iswapsub,onesub)
#                       swaptwosublists = picknfromlist(iswapsub,twosub)
#                       for swaponesublist in swaponesublists:
#                          for swaptwosublist in swaptwosublists:
#                             permutation = []
#                             permutation = permutation + copy.deepcopy(targetsuper) + copy.deepcopy(targetsub)
#                             parity = 1.0
#                             for targetindex in targetsuper:
#                                swapped = 0
#                                for nindexone in range(len(swaponesuperlist)):
#                                   indexone = swaponesuperlist[nindexone]
#                                   if (targetindex.isidenticalto(indexone)):
#                                      permutation.append(copy.deepcopy(swaptwosuperlist[nindexone]))
#                                      swapped = 1
#                                      parity = parity * (-1.0)
#                                for nindextwo in range(len(swaponesuperlist)):
#                                   indextwo = swaptwosuperlist[nindextwo]
#                                   if (targetindex.isidenticalto(indextwo)):
#                                      permutation.append(copy.deepcopy(swaponesuperlist[nindextwo]))
#                                      swapped = 1
#                                if (not swapped):
#                                   permutation.append(targetindex)
#                             for targetindex in targetsub:
#                                swapped = 0
#                                for nindexone in range(len(swaponesublist)):
#                                   indexone = swaponesublist[nindexone]
#                                   if (targetindex.isidenticalto(indexone)):
#                                      permutation.append(copy.deepcopy(swaptwosublist[nindexone]))
#                                      swapped = 1
#                                      parity = parity * (-1.0)
#                                for nindextwo in range(len(swaponesublist)):
#                                   indextwo = swaptwosublist[nindextwo]
#                                   if (targetindex.isidenticalto(indextwo)):
#                                      permutation.append(copy.deepcopy(swaponesublist[nindextwo]))
#                                      swapped = 1
#                                if (not swapped):
#                                   permutation.append(targetindex)
#                             identity = 1
#                             for nindex in range(len(permutation)/2):
#                                if (not permutation[nindex].isidenticalto(permutation[nindex+len(permutation)/2])):
#                                   identity = 0
#                             if (not identity):
#                                if (contraction == len(order)-2):
#                                   newfactor.add(Factor([parity*another.factor.coefficients[0]],[permutation]))
#                                else:
#                                   newfactor.add(Factor([parity],[permutation]))
#           newfactor = newfactor.canonicalize(onesuper)
#           newfactor = newfactor.canonicalize(onesub)
#           newfactor = newfactor.canonicalize(twosuper)
#           newfactor = newfactor.canonicalize(twosub)
            suggestedfactors.append(newfactor)
            factorproduct = factorproduct.product(newfactor)
            tensorthree = tensorone.contracts(tensortwo,label)
            tensorone = copy.deepcopy(tensorthree)
         for tensor in self.tensors:
            super = []
            for nindex in range(len(tensor.indexes)/2):
               index = tensor.indexes[nindex]
               common = 0
               if (self.summation):
                  for anotherindex in self.summation.indexes:
                     if (index.isidenticalto(anotherindex)):
                        common = 1
               if (not common):
                  super.append(tensor.indexes[nindex])
            factorproduct = factorproduct.canonicalize(super)
            sub = []
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               index = tensor.indexes[nindex]
               common = 0
               if (self.summation):
                  for anotherindex in self.summation.indexes:
                     if (index.isidenticalto(anotherindex)):
                        common = 1
               if (not common):
                  sub.append(tensor.indexes[nindex])
            factorproduct = factorproduct.canonicalize(sub)
         if (factorproduct.isthesameas(another.factor)):
            factorisbrokendown = 1
            print " ... the suggested decomposition of the permutation operator is valid"
         else:
            factorisbrokendown = 0
            raise RuntimeError, "Cannot breakdown permutation operator"

         # breakdown of a multiple tensor contraction
         tensorone = another.tensors[order[0]-1]
         for contraction in range(0,len(order)-1):
            tensortwo = another.tensors[order[contraction+1]-1]
            tensorthree = tensorone.contracts(tensortwo,label)
            label = label - 1
            if (factorisbrokendown):
               factor = suggestedfactors[contraction]
            else:
               if (contraction == len(order)-2):
                  factor = copy.deepcopy(another.factor)
               else:
                  factor = Factor([1.0],[[]])
            summation = Summation(tensorone.hascommonindexeswith(tensortwo))
            elementary = ElementaryTensorContraction(factor,summation,[tensorthree,tensorone,tensortwo])
            result = OperationTree(elementary,[],[result])
            tensorone = copy.deepcopy(tensorthree)
         result.sortindexes()
         return result

   def canonicalize(self,verbose=0):
      """ Canonicalizes the permutation operators"""
 
      if (verbose):
         print " ... canonicalizing permutation operator expressions"
      another = self.duplicate()
      for tensor in self.tensors:
         super = []
         for nindex in range(len(tensor.indexes)/2):
            index = tensor.indexes[nindex]
            common = 0
            if (self.summation):
               for anotherindex in self.summation.indexes:
                  if (index.isidenticalto(anotherindex)):
                     common = 1
            if (not common):
               super.append(tensor.indexes[nindex])
         another.factor = another.factor.canonicalize(super)
         sub = []
         for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
            index = tensor.indexes[nindex]
            common = 0
            if (self.summation):
               for anotherindex in self.summation.indexes:
                  if (index.isidenticalto(anotherindex)):
                     common = 1
            if (not common):
               sub.append(tensor.indexes[nindex])
         another.factor = another.factor.canonicalize(sub)

      return another
 
   def symmetrize(self,verbose=0):
      """Introduces a permutation operator that permutes equivalent target indexes of permutable amplitudes"""

      symmetrized = 0

      # target indexes
      targetsuper = []
      targetsub = []
      for tensor in self.tensors:
         for nindex in range(len(tensor.indexes)/2):
            index = tensor.indexes[nindex]
            common = 0
            if (self.summation):
               for anotherindex in self.summation.indexes:
                  if (index.isidenticalto(anotherindex)):
                     common = 1
            if (not common):
               targetsuper.append(tensor.indexes[nindex])
         for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
            index = tensor.indexes[nindex]
            common = 0
            if (self.summation):
               for anotherindex in self.summation.indexes:
                  if (index.isidenticalto(anotherindex)):
                     common = 1
            if (not common):
               targetsub.append(tensor.indexes[nindex])

      # identify equivalent tensors
#     another = self.duplicate()
#     for itensora in range(len(another.tensors)):
#        tensora = another.tensors[itensora]
#        super = []
#        sub = []
#        for itensorb in range(len(another.tensors)):
#           tensorb = another.tensors[itensorb]
#           if (itensorb <= itensora):
#              continue
#           if (tensora.isequivalentto(tensorb,another)):
#              if ((not super) and (not sub)):
#                 for nindex in range(len(tensora.indexes)/2):
#                    index = tensora.indexes[nindex]
#                    if (not another.summation.hastheindex(index)):
#                       super.append(index)
#                 for nindex in range(len(tensora.indexes)/2,len(tensora.indexes)):
#                    index = tensora.indexes[nindex]
#                    if (not another.summation.hastheindex(index)):
#                       sub.append(index)
#              for nindex in range(len(tensorb.indexes)/2):
#                 index = tensorb.indexes[nindex]
#                 if (not another.summation.hastheindex(index)):
#                    super.append(index)
#              for nindex in range(len(tensorb.indexes)/2,len(tensorb.indexes)):
#                 index = tensorb.indexes[nindex]
#                 if (not another.summation.hastheindex(index)):
#                    sub.append(index)
#        # multiply permutation operator of norm unity
#        if (super):
#           symmetrized = 1
#           factor = createfactor(super,targetsuper+targetsub)
#           if (verbose):
#              print factor.show()," will be multiplied"
#           another.factor = factor.product(another.factor)
#        if (sub):
#           symmetrized = 1
#           factor = createfactor(sub,targetsuper+targetsub)
#           if (verbose):
#              print factor.show()," will be multiplied"
#           another.factor = factor.product(another.factor)

      # identify equivalent tensors
      another = self.duplicate()
      for itensora in range(len(another.tensors)):
         tensora = another.tensors[itensora]
         before = []
         for nindex in range(len(tensora.indexes)/2):
            index = tensora.indexes[nindex]
            if (another.summation):
               if (not another.summation.hastheindex(index)):
                  before.append(index)
         for nindex in range(len(tensora.indexes)/2,len(tensora.indexes)):
            index = tensora.indexes[nindex]
            if (another.summation):
               if (not another.summation.hastheindex(index)):
                  before.append(index)
         for itensorb in range(len(another.tensors)):
            tensorb = another.tensors[itensorb]
            if (itensorb <= itensora):
               continue
            if (tensora.isequivalentto(tensorb,another)):
               after = []
               for nindex in range(len(tensorb.indexes)/2):
                  index = tensorb.indexes[nindex]
                  if (another.summation):
                     if (not another.summation.hastheindex(index)):
                        after.append(index)
               for nindex in range(len(tensorb.indexes)/2,len(tensorb.indexes)):
                  index = tensorb.indexes[nindex]
                  if (another.summation):
                     if (not another.summation.hastheindex(index)):
                        after.append(index)

               # multiply permutation operator of norm unity
               if (before):
                  symmetrized = 1
                  factor = createasmallfactor(before,after,targetsuper+targetsub)
                  if (verbose):
                     print factor.show()," will be multiplied"
                  another.factor = factor.product(another.factor)

      if (symmetrized):
         print " ... expression has been symmetrized"
      if (verbose):
         if (symmetrized):
            print self.show()
            print another.show()
            print ""

      return another

class ListTensorContractions:
 
   def __init__(self):
      """Creates a list of tensor contractions"""
      self.list = []
 
   def __str__(self):
      """Prints the content"""
      print ""
      for line in self.show():
         print line
      return ""
 
   def show(self,verbose=1):
      """Prints the tensor contractions"""
      show = []
      for tensorcontraction in self.list:
         show.append(tensorcontraction.show(verbose))
      return show
 
   def tex(self,verbose=1):
      """Prints the tensor contractions"""
      show = []
      for ntensorcontraction in range(len(self.list)):
         tensorcontraction = self.list[ntensorcontraction]
         if (ntensorcontraction == 0):
            show.append("\\begin{eqnarray}")
            show.append(string.join(["&&",tensorcontraction.tex(verbose),"\\nonumber\\\\"],""))
         elif (ntensorcontraction == len(self.list)-1):
            show.append(string.join(["&&",tensorcontraction.tex(verbose),"\\nonumber"],""))
            show.append("\\end{eqnarray}")
         else:
            show.append(string.join(["&&",tensorcontraction.tex(verbose),"\\nonumber\\\\"],""))
      return show

   def excitationtensortypes(self,de=0):
      """Returns a list of tensor types that are de/excitation"""
      alltypes = []
      for tensorcontraction in self.list:
         for tensor in tensorcontraction.tensors:
            if (tensor.type not in alltypes):
               alltypes.append(tensor.type)
      excitationtypes = []
      for type in alltypes:
         excitation = 1
         deexcitation = 1
         for tensorcontraction in self.list:
            for tensor in tensorcontraction.tensors:
               if (tensor.type == type):
                  for nindex in range(len(tensor.indexes)/2):
                     if (tensor.indexes[nindex].ishole()):
                        excitation = 0
                  for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
                     if (tensor.indexes[nindex].isparticle()):
                        excitation = 0
                  for nindex in range(len(tensor.indexes)/2):
                     if (tensor.indexes[nindex].isparticle()):
                        deexcitation = 0
                  for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
                     if (tensor.indexes[nindex].ishole()):
                        deexcitation = 0
         if ((not de) and (excitation)):
            excitationtypes.append(type)
         if ((de) and (deexcitation)):
            excitationtypes.append(type)
         
      return excitationtypes

   def findthebestbreakdown(self,verbose=0):
      """Iteratively calls findthebestbreakdown() for all tensor contractions"""
      for tensorcontraction in self.list:
         print tensorcontraction.show()
         tensorcontraction.findthebestbreakdown(self.excitationtensortypes(),verbose)

   def breakdown(self,verbose=0):
      """Iteratively calls breakdown() for all tensor contractions"""
      children = []
      for tensorcontraction in self.list:
         if (len(tensorcontraction.tensors) == 1):
            label = 0
         else:
            label = len(tensorcontraction.tensors) - 2
         newchild = tensorcontraction.breakdown(label,self.excitationtensortypes(),verbose)
         children.append(newchild)
      result = OperationTree(NoOperation(),[],children)
      return result

   def canonicalize(self):
      """Calls canonicalize() for all tensor contractions"""
      for ntensorcontraction in range(len(self.list)):
         tensorcontraction = self.list[ntensorcontraction]
         self.list[ntensorcontraction] = tensorcontraction.canonicalize(0)
      return self

   def symmetrize(self,verbose=1):
      """Calls symmetrize() for all tensor contractions"""
      for ntensorcontraction in range(len(self.list)):
         tensorcontraction = self.list[ntensorcontraction]
         self.list[ntensorcontraction] = tensorcontraction.symmetrize(verbose)
      return self

   def pythongen(self,filename="NONAME"):
      """Genrates a python code for debugging purposes"""

      pythoncode = []
      newline = "# This is a Python program generated by Tensor Contraction Engine v.1.0"
      pythoncode.append(newline)
      newline = "# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)"
      pythoncode.append(newline)

      tensors = ["i0"]
      for tensorcontraction in self.list:
         for tensor in tensorcontraction.tensors:
            newtensor = tensor.type+repr(len(tensor.indexes)/2)
            if (not (newtensor in tensors)):
               tensors.append(newtensor)
      newline = string.join(["def ",filename,"(N,nall,nocc"],"")
      for tensor in tensors:
         newline = string.join([newline,",",tensor],"")
      newline = string.join([newline,"):"],"")
      pythoncode.append(newline)

      for tensorcontraction in self.list:

         # generate a target tensor
         super = []
         sub = []
         for tensor in tensorcontraction.tensors:
            for nindex in range(len(tensor.indexes)/2):
               index = tensor.indexes[nindex]
               common = 0
               if (tensorcontraction.summation):
                  for another in tensorcontraction.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  super.append(tensor.indexes[nindex]) 
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               index = tensor.indexes[nindex]
               common = 0
               if (tensorcontraction.summation):
                  for another in tensorcontraction.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  sub.append(tensor.indexes[nindex]) 
         intermediate = Tensor("i",super+sub,0)
         parity = intermediate.sortindexes()

         # generate loops over target indexes
         indent = 1
         for index in intermediate.indexes:
            spin = string.join(["spin",repr(index.label)],"")
            newline = string.join(["for ",spin," in range(2):"],"")
            newline = string.join([" "*indent,newline],"")
            pythoncode.append(newline)
            indent = indent + 1
            if (index.type == 'hole'):
               newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nocc[",spin,"]):"],"")
               newline = string.join([" "*indent,newline],"")
               pythoncode.append(newline)
               indent = indent + 1
            elif (index.type == 'particle'):
               newline = string.join(["for ",index.show()," in range(",spin,"*nall[0]+nocc[",spin,"],",spin,"*nall[0]+nall[",spin,"]):"],"")
               newline = string.join([" "*indent,newline],"")
               pythoncode.append(newline)
               indent = indent + 1
            elif (index.type == 'general'):
               newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nall[",spin,"]):"],"")
               newline = string.join([" "*indent,newline],"")
               pythoncode.append(newline)
               indent = indent + 1

         # generate loops over common indexes
         if (tensorcontraction.summation):
            for index in tensorcontraction.summation.indexes:
               spin = string.join(["spin",repr(index.label)],"")
               newline = string.join(["for ",spin," in range(2):"],"")
               newline = string.join([" "*indent,newline],"")
               pythoncode.append(newline)
               indent = indent + 1
               if (index.type == 'hole'):
                  newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nocc[",spin,"]):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.append(newline)
                  indent = indent + 1
               elif (index.type == 'particle'):
                  newline = string.join(["for ",index.show()," in range(",spin,"*nall[0]+nocc[",spin,"],",spin,"*nall[0]+nall[",spin,"]):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.append(newline)
                  indent = indent + 1
               elif (index.type == 'general'):
                  newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nall[",spin,"]):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.append(newline)
                  indent = indent + 1
         newline = string.join([intermediate.pythongen(),"=",intermediate.pythongen(),"+",\
                                "(",repr(tensorcontraction.factor),")"],"")
         for tensor in tensorcontraction.tensors:
            newline = string.join([newline,"*",tensor.pythongen()],"")
         newline = string.join([" "*indent,newline],"")
         pythoncode.append(newline)
         
      # dump the code to a file
      writetofile(pythoncode,string.join([filename,".py.out"],""))

class NoOperation:

   def __init__(self):
      """Creates a dummy operation"""

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self,verbose=1):
      """Returns a dummy operation expression"""
      if (verbose):
         show = "No operation"
         return show
      else:
         return

   def tex(self):
      """Returns nothing"""
      return

   def isoperation(self):
      """Returns false"""
      return 0

   def usesindexlabel(self,label):
      """Returns false"""
      return 0

   def relabels(self,oldlabel,newlabel):
      """Just returns"""

   def swapindexes(self,indexone,indextwo):
      """Just returns"""

class ElementaryTensorContraction:

   def __init__(self,factor=[],summation=[],tensors=[]):
      """Creates an elementary tensor contraction: tensor A = factor * tensor B * tensor C"""
      self.factor = factor
      self.summation = summation
      self.tensors = tensors
 
   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self,verbose=1):
      """Returns a human-friendly string of the content"""
      show = string.join([self.tensors[0].show(),"+ ="])
      show = string.join([show,self.factor.show(verbose)])
      if (self.summation):
         if (len(self.summation.indexes) > 0):
            show = string.join([show, "*", self.summation.show()])
      show = string.join([show, "*", self.tensors[1].show()])
      if (len(self.tensors) == 3):
         show = string.join([show, "*", self.tensors[2].show()])
      return show 

   def tex(self,verbose=1):
      """Returns a LaTeX string of the content"""
      show = string.join([self.tensors[0].tex(),"=",self.tensors[0].tex()])
      show = string.join([show,self.factor.tex(verbose)])
#     if (self.summation):
#        if (len(self.summation.indexes) > 0):
#           show = string.join([show, self.summation.tex()])
      show = string.join([show, self.tensors[1].tex()])
      if (len(self.tensors) == 3):
         show = string.join([show, self.tensors[2].tex()])
      return show 

   def textable(self,name=""):
      """Returns a LaTeX string of the content"""
      show = self.factor.tex(0)
      show = string.join([show, self.tensors[1].textable(name)])
      if (len(self.tensors) == 3):
         show = string.join([show, self.tensors[2].textable(name)])
      return show 

   def isoperation(self):
      """Returns true"""
      return 1

   def usesindexlabel(self,label):
      """Returns true if the index label is already in use"""
      for tensor in self.tensors:
         if (tensor.usesindexlabel(label)):
            return 1
      return 0

   def relabels(self,oldlabel,newlabel):
      """Renames an index label"""
      for tensor in self.tensors:
         tensor.relabels(oldlabel,newlabel)
      if (self.summation):
         for nindex in range(len(self.summation.indexes)):
            index = self.summation.indexes[nindex]
            if (index.label == oldlabel):
               self.summation.indexes[nindex].label = newlabel

   def swapindexes(self,indexone,indextwo):
      """Swaps indexes"""
      for tensor in self.tensors:
         tensor.swapindexes(indexone,indextwo)
         if (self.summation):
            for nindex in range(len(self.summation.indexes)):
               index = self.summation.indexes[nindex]
               if (index.isidenticalto(indexone)):
                  self.summation.indexes[nindex] = copy.deepcopy(indextwo)
               elif (index.isidenticalto(indextwo)):
                  self.summation.indexes[nindex] = copy.deepcopy(indexone)

   def sortindexes(self):
      """Sorts the indexes of tensors taking account of parities"""
      for tensor in self.tensors:
         self.factor.multiply(tensor.sortindexes())

   def canonicalize(self,globaltargetindexes):
      """Canonicalizes the expression"""

      # reorder tensors
      superone = []
      subone = []
      for nindex in range(len(self.tensors[1].indexes)/2):
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            superone.append(index)
      for nindex in range(len(self.tensors[1].indexes)/2,len(self.tensors[1].indexes)):
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            subone.append(index)
      nsuper1 = len(superone)
      nsub1 = len(subone)
      if (len(self.tensors) == 3):
         supertwo = []
         subtwo = []
         for nindex in range(len(self.tensors[2].indexes)/2):
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               supertwo.append(index)
         for nindex in range(len(self.tensors[2].indexes)/2,len(self.tensors[2].indexes)):
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               subtwo.append(index)
         nsuper2 = len(supertwo)
         nsub2 = len(subtwo)
         if (self.tensors[1].conjugate < self.tensors[2].conjugate):
            swap = 1
         elif (self.tensors[1].conjugate > self.tensors[2].conjugate):
            swap = 0
         elif ((self.tensors[1].type == 'f') and (self.tensors[2].type != 'f')):
            swap = 1
         elif ((self.tensors[1].type != 'f') and (self.tensors[2].type == 'f')):
            swap = 0
         elif ((self.tensors[1].type == 'v') and (self.tensors[2].type != 'v')):
            swap = 1
         elif ((self.tensors[1].type != 'v') and (self.tensors[2].type == 'v')):
            swap = 0
         elif ((self.tensors[1].type == 'i') and (self.tensors[2].type != 'i')):
            swap = 1
         elif ((self.tensors[1].type != 'i') and (self.tensors[2].type == 'i')):
            swap = 0
         elif (self.tensors[1].type > self.tensors[2].type):
            swap = 1
         elif (self.tensors[1].type < self.tensors[2].type):
            swap = 0
         elif (len(self.tensors[1].indexes) < len(self.tensors[2].indexes)):
            swap = 1
         elif (len(self.tensors[1].indexes) > len(self.tensors[2].indexes)):
            swap = 0
         elif (nsuper1 + nsub1 < nsuper2 + nsub2):
            swap = 1
         elif (nsuper1 + nsub1 > nsuper2 + nsub2):
            swap = 0
         elif (nsuper1 < nsuper2):
            swap = 1
         elif (nsuper1 > nsuper2):
            swap = 0
         else:
            swap = 0
         if (swap):
            swaptensor = copy.deepcopy(self.tensors[1])
            self.tensors[1] = copy.deepcopy(self.tensors[2])
            self.tensors[2] = copy.deepcopy(swaptensor)

      indexesintheoriginalorder = []
      for index in self.tensors[1].indexes:
         if (index.isin(globaltargetindexes)):
            indexesintheoriginalorder.append(index)
      if (len(self.tensors) == 3):
         for index in self.tensors[2].indexes:
            if (index.isin(globaltargetindexes)):
               indexesintheoriginalorder.append(index)

      indexesintheneworder = []
      for index in indexesintheoriginalorder:
         minimum = "empty"
         for targetindex in globaltargetindexes:
            if ((targetindex.type == index.type) and (not targetindex.isin(indexesintheneworder))):
               if (minimum == "empty"):
                  minimum = copy.deepcopy(targetindex)
               else:
                  if (minimum.label > targetindex.label):
                     minimum = copy.deepcopy(targetindex)
         if (minimum == "empty"):
            print self
            printindexes(globaltargetindexes) 
            raise RuntimeError, "unable to canonicalize"
         indexesintheneworder.append(minimum)

      for nindex in range(len(self.tensors[0].indexes)):
         index = self.tensors[0].indexes[nindex]
         for n in range(len(indexesintheoriginalorder)):
            if (index.isidenticalto(indexesintheoriginalorder[n])):
               self.tensors[0].indexes[nindex] = copy.deepcopy(indexesintheneworder[n])
      self.factor.multiply(self.tensors[0].sortindexes())
      for nindex in range(len(self.tensors[1].indexes)):
         index = self.tensors[1].indexes[nindex]
         for n in range(len(indexesintheoriginalorder)):
            if (index.isidenticalto(indexesintheoriginalorder[n])):
               self.tensors[1].indexes[nindex] = copy.deepcopy(indexesintheneworder[n])
      self.factor.multiply(self.tensors[1].sortindexes())
      if (len(self.tensors) == 3):
         for nindex in range(len(self.tensors[2].indexes)):
            index = self.tensors[2].indexes[nindex]
            for n in range(len(indexesintheoriginalorder)):
               if (index.isidenticalto(indexesintheoriginalorder[n])):
                  self.tensors[2].indexes[nindex] = copy.deepcopy(indexesintheneworder[n])
         self.factor.multiply(self.tensors[2].sortindexes())
      for npermutation in range(len(self.factor.permutations)):
         permutation = self.factor.permutations[npermutation]
         newpermutationorigin = []
         newpermutationdestination = []
         for nindex in range(len(permutation)/2):
            index = permutation[nindex]
            for n in range(len(indexesintheoriginalorder)):
               if (index.isidenticalto(indexesintheoriginalorder[n])):
                  newpermutationorigin.append(indexesintheneworder[n])
         for nindex in range(len(permutation)/2,len(permutation)):
            index = permutation[nindex]
            for n in range(len(indexesintheoriginalorder)):
               if (index.isidenticalto(indexesintheoriginalorder[n])):
                  newpermutationdestination.append(indexesintheneworder[n])
         for index in globaltargetindexes:
            if (not index.isin(newpermutationorigin)):
               newpermutationorigin.append(index)
               newpermutationdestination.append(index)
         self.factor.permutations[npermutation] = copy.deepcopy(newpermutationorigin + newpermutationdestination)
      newsuperone = []
      newsubone = []
      for nindex in range(len(self.tensors[1].indexes)/2):
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            newsuperone.append(index)
      for nindex in range(len(self.tensors[1].indexes)/2,len(self.tensors[1].indexes)):
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            newsubone.append(index)
      if (len(self.tensors) == 3):
         newsupertwo = []
         newsubtwo = []
         for nindex in range(len(self.tensors[2].indexes)/2):
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               newsupertwo.append(index)
         for nindex in range(len(self.tensors[2].indexes)/2,len(self.tensors[2].indexes)):
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               newsubtwo.append(index)
      self.factor = self.factor.canonicalize(newsuperone)
      self.factor = self.factor.canonicalize(newsubone)
      if (len(self.tensors) == 3):
         self.factor = self.factor.canonicalize(newsupertwo)
         self.factor = self.factor.canonicalize(newsubtwo)

   def isthesameas(self,another,globaltargetindexes,verbose=0):
      """Returns a ratio when the two contractions are essentially the same (returns 0 otherwise)"""
      if (verbose):
         print " ... comparing"
         print self
         print another
      self.canonicalize(globaltargetindexes)
      another.canonicalize(globaltargetindexes)
      if ((len(self.tensors) == 2) or (len(another.tensors) == 2)):
         return 0
      if ((self.tensors[1].type == "i") or (self.tensors[2].type == "i")):
         return 0
      if ((another.tensors[1].type == "i") or (another.tensors[2].type == "i")):
         return 0
      if ((self.tensors[1].type == "j") or (self.tensors[2].type == "j")):
         return 0
      if ((another.tensors[1].type == "j") or (another.tensors[2].type == "j")):
         return 0
      if (self.tensors[1].type != another.tensors[1].type):
         return 0
      if (self.tensors[2].type != another.tensors[2].type):
         return 0
      if (len(self.tensors[1].indexes) != len(another.tensors[1].indexes)):
         return 0
      if (len(self.tensors[2].indexes) != len(another.tensors[2].indexes)):
         return 0
      if (len(self.summation.indexes) != len(another.summation.indexes)):
         return 0
      ratio = self.factor.isthesameas(another.factor)
      if (ratio == 0):
         return 0
      reserved = []
      for tensor in self.tensors:
         for index in tensor.indexes:
            label = index.label
            if (label not in reserved):
               reserved.append(label)
      for tensor in another.tensors:
         for index in tensor.indexes:
            label = index.label
            if (label not in reserved):
               reserved.append(label)
      newlabel = 0
      for label in reserved:
         if (label > newlabel):
            newlabel = label
      selfcopy = copy.deepcopy(self)
      anothercopy = copy.deepcopy(another)
      for ntensor in range(len(self.tensors)):
         tensor = self.tensors[ntensor]
         for nindex in range(len(tensor.indexes)):
            selfindex = self.tensors[ntensor].indexes[nindex]
            anotherindex = another.tensors[ntensor].indexes[nindex]
            if (selfindex.type != anotherindex.type):
               return 0
            if (selfindex.label != anotherindex.label):
               newlabel = newlabel + 1
               selfcopy.relabels(selfindex.label,newlabel)
               anothercopy.relabels(anotherindex.label,newlabel)
      for nindex in range(len(selfcopy.tensors[1].indexes)/2):
         for mindex in range(len(selfcopy.tensors[1].indexes)/2):
            if (mindex > nindex):
               if (selfcopy.tensors[1].indexes[nindex].isgreaterthan(selfcopy.tensors[1].indexes[mindex])):
                  swap = copy.deepcopy(selfcopy.tensors[1].indexes[nindex])
                  selfcopy.tensors[1].indexes[nindex] = copy.deepcopy(selfcopy.tensors[1].indexes[mindex])
                  selfcopy.tensors[1].indexes[mindex] = copy.deepcopy(swap)
                  selfcopy.factor.multiply(-1.0)
      for nindex in range(len(selfcopy.tensors[1].indexes)/2,len(selfcopy.tensors[1].indexes)):
         for mindex in range(len(selfcopy.tensors[1].indexes)/2,len(selfcopy.tensors[1].indexes)):
            if (mindex > nindex):
               if (selfcopy.tensors[1].indexes[nindex].isgreaterthan(selfcopy.tensors[1].indexes[mindex])):
                  swap = copy.deepcopy(selfcopy.tensors[1].indexes[nindex])
                  selfcopy.tensors[1].indexes[nindex] = copy.deepcopy(selfcopy.tensors[1].indexes[mindex])
                  selfcopy.tensors[1].indexes[mindex] = copy.deepcopy(swap)
                  selfcopy.factor.multiply(-1.0)
      for nindex in range(len(selfcopy.tensors[2].indexes)/2):
         for mindex in range(len(selfcopy.tensors[2].indexes)/2):
            if (mindex > nindex):
               if (selfcopy.tensors[2].indexes[nindex].isgreaterthan(selfcopy.tensors[2].indexes[mindex])):
                  swap = copy.deepcopy(selfcopy.tensors[2].indexes[nindex])
                  selfcopy.tensors[2].indexes[nindex] = copy.deepcopy(selfcopy.tensors[2].indexes[mindex])
                  selfcopy.tensors[2].indexes[mindex] = copy.deepcopy(swap)
                  selfcopy.factor.multiply(-1.0)
      for nindex in range(len(selfcopy.tensors[2].indexes)/2,len(selfcopy.tensors[2].indexes)):
         for mindex in range(len(selfcopy.tensors[2].indexes)/2,len(selfcopy.tensors[2].indexes)):
            if (mindex > nindex):
               if (selfcopy.tensors[2].indexes[nindex].isgreaterthan(selfcopy.tensors[2].indexes[mindex])):
                  swap = copy.deepcopy(selfcopy.tensors[2].indexes[nindex])
                  selfcopy.tensors[2].indexes[nindex] = copy.deepcopy(selfcopy.tensors[2].indexes[mindex])
                  selfcopy.tensors[2].indexes[mindex] = copy.deepcopy(swap)
                  selfcopy.factor.multiply(-1.0)
      for nindex in range(len(anothercopy.tensors[1].indexes)/2):
         for mindex in range(len(anothercopy.tensors[1].indexes)/2):
            if (mindex > nindex):
               if (anothercopy.tensors[1].indexes[nindex].isgreaterthan(anothercopy.tensors[1].indexes[mindex])):
                  swap = copy.deepcopy(anothercopy.tensors[1].indexes[nindex])
                  anothercopy.tensors[1].indexes[nindex] = copy.deepcopy(anothercopy.tensors[1].indexes[mindex])
                  anothercopy.tensors[1].indexes[mindex] = copy.deepcopy(swap)
                  anothercopy.factor.multiply(-1.0)
      for nindex in range(len(anothercopy.tensors[1].indexes)/2,len(anothercopy.tensors[1].indexes)):
         for mindex in range(len(anothercopy.tensors[1].indexes)/2,len(anothercopy.tensors[1].indexes)):
            if (mindex > nindex):
               if (anothercopy.tensors[1].indexes[nindex].isgreaterthan(anothercopy.tensors[1].indexes[mindex])):
                  swap = copy.deepcopy(anothercopy.tensors[1].indexes[nindex])
                  anothercopy.tensors[1].indexes[nindex] = copy.deepcopy(anothercopy.tensors[1].indexes[mindex])
                  anothercopy.tensors[1].indexes[mindex] = copy.deepcopy(swap)
                  anothercopy.factor.multiply(-1.0)
      for nindex in range(len(anothercopy.tensors[2].indexes)/2):
         for mindex in range(len(anothercopy.tensors[2].indexes)/2):
            if (mindex > nindex):
               if (anothercopy.tensors[2].indexes[nindex].isgreaterthan(anothercopy.tensors[2].indexes[mindex])):
                  swap = copy.deepcopy(anothercopy.tensors[2].indexes[nindex])
                  anothercopy.tensors[2].indexes[nindex] = copy.deepcopy(anothercopy.tensors[2].indexes[mindex])
                  anothercopy.tensors[2].indexes[mindex] = copy.deepcopy(swap)
                  anothercopy.factor.multiply(-1.0)
      for nindex in range(len(anothercopy.tensors[2].indexes)/2,len(anothercopy.tensors[2].indexes)):
         for mindex in range(len(anothercopy.tensors[2].indexes)/2,len(anothercopy.tensors[2].indexes)):
            if (mindex > nindex):
               if (anothercopy.tensors[2].indexes[nindex].isgreaterthan(anothercopy.tensors[2].indexes[mindex])):
                  swap = copy.deepcopy(anothercopy.tensors[2].indexes[nindex])
                  anothercopy.tensors[2].indexes[nindex] = copy.deepcopy(anothercopy.tensors[2].indexes[mindex])
                  anothercopy.tensors[2].indexes[mindex] = copy.deepcopy(swap)
                  anothercopy.factor.multiply(-1.0)
      reserved = []
      for tensor in selfcopy.tensors:
         for index in tensor.indexes:
            label = index.label
            if (label not in reserved):
               reserved.append(label)
      for tensor in anothercopy.tensors:
         for index in tensor.indexes:
            label = index.label
            if (label not in reserved):
               reserved.append(label)
      newlabel = 0
      for label in reserved:
         if (label > newlabel):
            newlabel = label
      for ntensor in range(len(selfcopy.tensors)):
         tensor = selfcopy.tensors[ntensor]
         for nindex in range(len(tensor.indexes)):
            selfindex = selfcopy.tensors[ntensor].indexes[nindex]
            anotherindex = anothercopy.tensors[ntensor].indexes[nindex]
            if (selfindex.type != anotherindex.type):
               return 0
            if (selfindex.label != anotherindex.label):
               newlabel = newlabel + 1
               selfcopy.relabels(selfindex.label,newlabel)
               anothercopy.relabels(anotherindex.label,newlabel)
      for nindex in range(len(selfcopy.summation.indexes)):
         for mindex in range(len(selfcopy.summation.indexes)):
            if (mindex > nindex):
               if (selfcopy.summation.indexes[nindex].isgreaterthan(selfcopy.summation.indexes[mindex])):
                  swap = copy.deepcopy(selfcopy.summation.indexes[nindex])
                  selfcopy.summation.indexes[nindex] = copy.deepcopy(selfcopy.summation.indexes[mindex])
                  selfcopy.summation.indexes[mindex] = copy.deepcopy(swap)
      for nindex in range(len(anothercopy.summation.indexes)):
         for mindex in range(len(anothercopy.summation.indexes)):
            if (mindex > nindex):
               if (anothercopy.summation.indexes[nindex].isgreaterthan(anothercopy.summation.indexes[mindex])):
                  swap = copy.deepcopy(anothercopy.summation.indexes[nindex])
                  anothercopy.summation.indexes[nindex] = copy.deepcopy(anothercopy.summation.indexes[mindex])
                  anothercopy.summation.indexes[mindex] = copy.deepcopy(swap)
      if (verbose):
         print " ... after canonicalization and relabeling"
         print selfcopy
         print anothercopy
      for ntensor in range(len(selfcopy.tensors)):
         tensor = selfcopy.tensors[ntensor]
         for nindex in range(len(tensor.indexes)):
            selfindex = selfcopy.tensors[ntensor].indexes[nindex]
            anotherindex = anothercopy.tensors[ntensor].indexes[nindex]
            if (selfindex.label != anotherindex.label):
               return 0
      for nindex in range(len(selfcopy.summation.indexes)):
         selfindex = selfcopy.summation.indexes[nindex]
         anotherindex = anothercopy.summation.indexes[nindex]
         if (selfindex.type != anotherindex.type):
            return 0
         if (selfindex.label != anotherindex.label):
            return 0
      if (verbose):
         print " ... essentially the same"
      return selfcopy.factor.isthesameas(anothercopy.factor)

   def fortran77(self,globaltargetindexes,types,subroutinename="NONAME"):
      """Suggests an implementation in Fortran77 for an elementary tensor contraction C = A * B"""

      verbose = 0

      if (verbose):
         print self
         print "globaltargetindexes"
         printindexes(globaltargetindexes)

      if (len(self.tensors) == 3):
         three = 1
      else:
         three = 0

      errquit = 0

      newcode = Code("Fortran77",subroutinename)

      # header
      newline = "!" + self.show(0)
      newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
      
      # insert include statements
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)
      newline = '#include "sym.fh"'
      newcode.add("headers",newline)
      newline = '#include "errquit.fh"'
      newcode.add("headers",newline)
      newline = '#include "tce.fh"'
      newcode.add("headers",newline)

      # declaration
      newint = "d_a"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "k_a_offset"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      if (three):
         newint = "d_b"
         newcode.add("integers",newint)
         newcode.add("arguments",newint)
         newint = "k_b_offset"
         newcode.add("integers",newint)
         newcode.add("arguments",newint)
      newint = "d_c"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "k_c_offset"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "NXTASK"
      newcode.add("integers",newint)
      newcode.add("externals",newint)
      newint = "next"
      newcode.add("integers",newint)
      newint = "nprocs"
      newcode.add("integers",newint)
      newint = "count"
      newcode.add("integers",newint)
            
      # Tensor 0
      superglobalzero = []
      subglobalzero = []
      superlocalzero = []
      sublocalzero = []
      for nindex in range(len(self.tensors[0].indexes)/2):
         index = self.tensors[0].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            superglobalzero.append(index)
         else:
            superlocalzero.append(index)
      for nindex in range(len(self.tensors[0].indexes)/2,len(self.tensors[0].indexes)):
         index = self.tensors[0].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            subglobalzero.append(index)
         else:
            sublocalzero.append(index)
      if (verbose):
         print "Tensor 0 superglobal"
         printindexes(superglobalzero)
         print "Tensor 0 subglobal"
         printindexes(subglobalzero)
         print "Tensor 0 superlocal"
         printindexes(superlocalzero)
         print "Tensor 0 sublocal"
         printindexes(sublocalzero)

      # Tensor 1
      superglobalone = []
      subglobalone = []
      superlocalone = []
      sublocalone = []
      supercommonone = []
      subcommonone = []
      if (not self.tensors[1].conjugate):
         superrangeone = range(len(self.tensors[1].indexes)/2)
         subrangeone = range(len(self.tensors[1].indexes)/2,len(self.tensors[1].indexes))
      else:
         superrangeone = range(len(self.tensors[1].indexes)/2,len(self.tensors[1].indexes))
         subrangeone = range(len(self.tensors[1].indexes)/2)
      for nindex in superrangeone:
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            superglobalone.append(index)
         elif (self.summation):
            if (index.isin(self.summation.indexes)):
               supercommonone.append(index)
            else:
               superlocalone.append(index)
         else:
            superlocalone.append(index)
      for nindex in subrangeone:
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            subglobalone.append(index)
         elif (self.summation):
            if (index.isin(self.summation.indexes)):
               subcommonone.append(index)
            else:
               sublocalone.append(index)
         else:
            sublocalone.append(index)
      if (verbose):
         print "Tensor 1 superglobal"
         printindexes(superglobalone)
         print "Tensor 1 subglobal"
         printindexes(subglobalone)
         print "Tensor 1 superlocal"
         printindexes(superlocalone)
         print "Tensor 1 sublocal"
         printindexes(sublocalone)
         print "Tensor 1 supercommon"
         printindexes(supercommonone)
         print "Tensor 1 subcommon"
         printindexes(subcommonone)

      # Tensor 2
      superglobaltwo = []
      subglobaltwo = []
      superlocaltwo = []
      sublocaltwo = []
      supercommontwo = []
      subcommontwo = []
      if (three):
         if (not self.tensors[2].conjugate):
            superrangetwo = range(len(self.tensors[2].indexes)/2)
            subrangetwo = range(len(self.tensors[2].indexes)/2,len(self.tensors[2].indexes))
         else:
            superrangetwo = range(len(self.tensors[2].indexes)/2,len(self.tensors[2].indexes))
            subrangetwo = range(len(self.tensors[2].indexes)/2)
         for nindex in superrangetwo:
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               superglobaltwo.append(index)
            elif (self.summation):
               if (index.isin(self.summation.indexes)):
                  supercommontwo.append(index)
               else:
                  superlocaltwo.append(index)
            else:
               superlocaltwo.append(index)
         for nindex in subrangetwo:
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               subglobaltwo.append(index)
            elif (self.summation):
               if (index.isin(self.summation.indexes)):
                  subcommontwo.append(index)
               else:
                  sublocaltwo.append(index)
            else:
               sublocaltwo.append(index)
         if (verbose):
            print "Tensor 2 superglobal"
            printindexes(superglobaltwo)
            print "Tensor 2 subglobal"
            printindexes(subglobaltwo)
            print "Tensor 2 superlocal"
            printindexes(superlocaltwo)
            print "Tensor 2 sublocal"
            printindexes(sublocaltwo)
            print "Tensor 2 supercommon"
            printindexes(supercommontwo)
            print "Tensor 2 subcommon"
            printindexes(subcommontwo)
      if (len(supercommonone) > len(subcommontwo)):
         supercommon = supercommonone
      else:
         supercommon = subcommontwo
      if (len(subcommonone) > len(supercommontwo)):
         subcommon = subcommonone
      else:
         subcommon = supercommontwo

      # parallel related
      newline = "nprocs = GA_NNODES()"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "count = 0"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "next = NXTASK(nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # loop over output tensor indexes
      newcode.inserttileddoloops(superglobalone)
      if (three):
         newcode.inserttileddoloops(superglobaltwo)
      newcode.inserttileddoloops(superlocalone)
      if (three):
         newcode.inserttileddoloops(superlocaltwo)
      newcode.inserttileddoloops(subglobalone)
      if (three):
         newcode.inserttileddoloops(subglobaltwo)
      newcode.inserttileddoloops(sublocalone)
      if (three):
         newcode.inserttileddoloops(sublocaltwo)

      # parallel related
      newline = "IF (next.eq.count) THEN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.setamark(4)
      newcode.pointer = newcode.pointer + 1
      newline = "next = NXTASK(nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END IF"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "count = count + 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.getamark(4) + 1

      # spin restriction on output tensor
      newcode.inserttiledifrestricted(superglobalzero + superlocalzero + subglobalzero + sublocalzero)

      # symmetry of output tensor
      super = superglobalzero + superlocalzero
      sub = subglobalzero + sublocalzero
      newcode.inserttiledifsymmetry(super,sub,self.tensors[0].irrep)

      # loop over summation indexes
      newcode.inserttileddoloops(supercommon)
      newcode.inserttileddoloops(subcommon)

      # symmetry of input tensor 1
      if (three):
         super = superglobalone + superlocalone + supercommonone
         sub = subglobalone + sublocalone + subcommonone
         newcode.inserttiledifsymmetry(super,sub,self.tensors[1].irrep)

      # spin restriction on input tensor one
      indexesone = superglobalone + superlocalone + supercommonone + subglobalone + sublocalone + subcommonone
      if (indexesone):
         newline = "IF ((restricted).and.("
         conjugation = ""
         for index in indexesone:
            newint = string.join([index.show(),"b"],"")
            newline = string.join([newline,conjugation,"int_mb(k_spin+",newint,"-1)"],"")
            conjugation = "+"
         newline = string.join([newline,".eq.",repr(2*len(indexesone)),")) THEN"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         for index in indexesone:
            newint = string.join([index.show(),"b_1"],"")
            newcode.add("integers",newint)
            newline = string.join([newint," = int_mb(k_alpha+",index.show(),"b-1)"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
         newline = "ELSE"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         for index in indexesone:
            newint = string.join([index.show(),"b_1"],"")
            newcode.add("integers",newint)
            newline = string.join([newint," = ",index.show(),"b"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
         newline = "END IF"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      # spin restriction on input tensor one
      if (three):
         indexestwo = superglobaltwo + superlocaltwo + supercommontwo + subglobaltwo + sublocaltwo + subcommontwo
         if (indexestwo):
            newline = "IF ((restricted).and.("
            conjugation = ""
            for index in indexestwo:
               newint = string.join([index.show(),"b"],"")
               newline = string.join([newline,conjugation,"int_mb(k_spin+",newint,"-1)"],"")
               conjugation = "+"
            newline = string.join([newline,".eq.",repr(2*len(indexestwo)),")) THEN"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for index in indexestwo:
               newint = string.join([index.show(),"b_2"],"")
               newcode.add("integers",newint)
               newline = string.join([newint," = int_mb(k_alpha+",index.show(),"b-1)"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "ELSE"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for index in indexestwo:
               newint = string.join([index.show(),"b_2"],"")
               newcode.add("integers",newint)
               newline = string.join([newint," = ",index.show(),"b"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "END IF"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1

      # create MA's for tensor 1
      newcode.add("integers","dim_common")
      newline = ""
      for index in supercommonone + subcommonone:
         if (newline == ""):
            newline = string.join(["dim_common = int_mb(k_range+",index.show(),"b-1)"],"")
         else:
            newline = string.join([newline," * int_mb(k_range+",index.show(),"b-1)"],"")
      if (newline == ""):
         newline = "dim_common = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.add("integers","dima_sort")
      newline = ""
      for index in superglobalone + superlocalone + subglobalone + sublocalone:
         if (newline == ""):
            newline = string.join(["dima_sort = int_mb(k_range+",index.show(),"b-1)"],"")
         else:
            newline = string.join([newline," * int_mb(k_range+",index.show(),"b-1)"],"")
      if (newline == ""):
         newline = "dima_sort = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.add("integers","dima")
      newline = "dima = dim_common * dima_sort"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # create MA's for tensor 2
      if (three):
         newcode.add("integers","dimb_sort")
         newline = ""
         for index in superglobaltwo + superlocaltwo + subglobaltwo + sublocaltwo:
            if (newline == ""):
               newline = string.join(["dimb_sort = int_mb(k_range+",index.show(),"b-1)"],"")
            else:
               newline = string.join([newline," * int_mb(k_range+",index.show(),"b-1)"],"")
         if (newline == ""):
            newline = "dimb_sort = 1"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newcode.add("integers","dimb")
         newline = "dimb = dim_common * dimb_sort"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      if (three):
         newline = "IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN"
      else:
         newline = "IF (dima .gt. 0) THEN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END IF"
      newcode.statements.insert(newcode.pointer,newline)

      # allocate sorted and unsorted tensor 1
      newcode.add("integers","l_a_sort")
      newcode.add("integers","k_a_sort")
      newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL ERRQUIT('",\
                             subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1
      newcode.add("integers","l_a")
      newcode.add("integers","k_a")
      newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('",\
                             subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1

      # mapping to a permutation symmetry unique block
      if (self.tensors[1].type in types[3]):
         superpermutations = restrictedpermutationwithparity(superlocalone,supercommonone,[])
      else:
         superpermutations = restrictedpermutationwithparity(superglobalone,superlocalone,supercommonone)
      if (self.tensors[1].type in types[3]):
         subpermutations = restrictedpermutationwithparity(sublocalone,subcommonone,[])
      else:
         subpermutations = restrictedpermutationwithparity(subglobalone,sublocalone,subcommonone)
      newcode.pointer = newcode.pointer - 1
      newcode.setamark(1)
      newcode.pointer = newcode.pointer + 1
      ifblock = 0
      for superpermutation in superpermutations:
         superline = ""
         if (self.tensors[1].type in types[3]):
            if (superpermutation[1] == "empty"):
               superpermutedindexes = superglobalone + sortindexes(superlocalone + supercommonone)
               superfactor = 1
            else:
               superpermutedindexes = superglobalone + superpermutation[1:]
               superfactor = parityofpermutation(superglobalone + sortindexes(superlocalone + supercommonone) + \
                                                 superglobalone + superpermutation[1:])
         else:
            if (superpermutation[1] == "empty"):
               superpermutedindexes = self.tensors[1].indexes[0:len(self.tensors[1].indexes)/2]
               superfactor = 1
            else:
               superpermutedindexes = superpermutation[1:]
               superfactor = parityofpermutation(self.tensors[1].indexes[0:len(self.tensors[1].indexes)/2] + \
                                                 superpermutation[1:])
         if (superpermutation[1] != "empty"):
            for nindex in range(len(superpermutedindexes)-1):
               indexa = superpermutedindexes[nindex]
               indexb = superpermutedindexes[nindex+1]
               if ((self.tensors[1].type in types[3]) and (indexa.isin(superglobalone) or indexb.isin(superglobalone))):
                  continue
               if (indexa.isin(superglobalone) and indexb.isin(superglobalone)):
                  continue
               if (indexa.isin(superlocalone) and indexb.isin(superlocalone)):
                  continue
               if (indexa.isin(supercommonone) and indexb.isin(supercommonone)):
                  continue
               if (indexa.isgreaterthan(indexb)):
                  inequality = " .lt. "
               else:
                  inequality = " .le. "
               if (superline):
                  superline = string.join([superline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
               else:
                  superline = string.join([superline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
         for subpermutation in subpermutations:
            subline = superline
            if (self.tensors[1].type in types[3]):
               if (subpermutation[1] == "empty"):
                  subpermutedindexes = subglobalone + sortindexes(sublocalone + subcommonone)
                  subfactor = 1
               else:
                  subpermutedindexes = subglobalone + subpermutation[1:]
                  subfactor = parityofpermutation(subglobalone + sortindexes(sublocalone + subcommonone) + \
                                                  subglobalone + subpermutation[1:])
            else:
               if (subpermutation[1] == "empty"):
                  subpermutedindexes = self.tensors[1].indexes[len(self.tensors[1].indexes)/2:len(self.tensors[1].indexes)]
                  subfactor = 1
               else:
                  subpermutedindexes = subpermutation[1:]
                  subfactor = parityofpermutation(self.tensors[1].indexes[len(self.tensors[1].indexes)/2:len(self.tensors[1].indexes)] + \
                                                  subpermutation[1:])
            if (subpermutation[1] != "empty"):
               for nindex in range(len(subpermutedindexes)-1):
                  indexa = subpermutedindexes[nindex]
                  indexb = subpermutedindexes[nindex+1]
                  if ((self.tensors[1].type in types[3]) and (indexa.isin(subglobalone) or indexb.isin(subglobalone))):
                     continue
                  if (indexa.isin(subglobalone) and indexb.isin(subglobalone)):
                     continue
                  if (indexa.isin(sublocalone) and indexb.isin(sublocalone)):
                     continue
                  if (indexa.isin(subcommonone) and indexb.isin(subcommonone)):
                     continue
                  if (indexa.isgreaterthan(indexb)):
                     inequality = " .lt. "
                  else:
                     inequality = " .le. "
                  if (subline):
                     subline = string.join([subline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
                  else:
                     subline = string.join([subline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
            if (subline):
               subline = string.join([subline,") THEN"],"")
               newcode.pointer = newcode.getamark(1) + 1
               newcode.statements.insert(newcode.pointer,subline)
               newcode.pointer = newcode.pointer + 1
               newline = "END IF"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.setamark(1)
               ifblock = 1

            permutedindexes = superpermutedindexes + subpermutedindexes

            # get a block
            arguments = ""
            argumentsend = ""
            for nindex in range(len(permutedindexes)-1,-1,-1):
               if (permutedindexes[nindex].type == "hole"):
                  boffset = "b_1 - 1"
               else:
                  if (self.tensors[1].type in types[4]):
                     boffset = "b_1 - 1"
                  else:
                     boffset = "b_1 - noab - 1"
               if (arguments == ""):
                  arguments = string.join(["d_a,dbl_mb(k_a),dima,int_mb(k_a_offset + ",permutedindexes[nindex].show(),boffset],"")
               else:
                  if (self.tensors[1].type in types[4]):
                     arguments = string.join([arguments," + (noab+nvab) * (",permutedindexes[nindex].show(),boffset],"")
                  else:
                     if (permutedindexes[nindex+1].type == "hole"):
                        arguments = string.join([arguments," + noab * (",permutedindexes[nindex].show(),boffset],"")
                     else:
                        arguments = string.join([arguments," + nvab * (",permutedindexes[nindex].show(),boffset],"")
                  argumentsend = string.join([argumentsend,")"],"")
            if (not arguments):
               arguments = "d_a,dbl_mb(k_a),dima,int_mb(k_a_offset)"
            else:
               arguments = string.join([arguments,argumentsend,")"],"")
            newline = string.join(["CALL GET_BLOCK(",arguments,")"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1

            # sort indexes of tensor 1
            doloopofa = 0
            for nindex in range(len(permutedindexes)):
               index = permutedindexes[nindex]
               newint = index.show()
               newcode.add("integers",newint)
               newline = string.join(["DO ",newint," = 1,int_mb(k_range+",newint,"b-1)"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "END DO"
               newcode.statements.insert(newcode.pointer,newline)
               doloopofa = 1
               if ((not ifblock) and (nindex == 0)):
                  newcode.setamark(1)
            newline = ""
            newlineend = ""
            newcode.add("integers","idima")
            for nindex in range(len(permutedindexes)-1,-1,-1):
               if (newline == ""):
                  newline = string.join(["idima = ",permutedindexes[nindex].show()],"")
               else:
                  newline = string.join([newline," + int_mb(k_range+",permutedindexes[nindex+1].show(),"b-1)"\
                                         " * ((",permutedindexes[nindex].show()," - 1)"],"")
                  newlineend = string.join([newlineend,")"],"")
            newline = string.join([newline,newlineend],"")
            if (not newline):
               newline = "idima = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            newline = ""
            newlineend = ""
            newcode.add("integers","idima_sort")
            sorted = supercommonone + subcommonone + superglobalone + superlocalone + subglobalone + sublocalone
            for nindex in range(len(sorted)):
               if (newline == ""):
                  newline = string.join(["idima_sort = ",sorted[nindex].show()],"")
               else:
                  newline = string.join([newline," + int_mb(k_range+",sorted[nindex-1].show(),"b-1)"\
                                         " * ((",sorted[nindex].show()," - 1)"],"")
                  newlineend = string.join([newlineend,")"],"")
            newline = string.join([newline,newlineend],"")
            if (not newline):
               newline = "idima_sort = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            if (superfactor * subfactor == 1):
               newline = "dbl_mb(k_a_sort + idima_sort - 1) = dbl_mb(k_a + idima - 1)"
            else:
               newline = "dbl_mb(k_a_sort + idima_sort - 1) = - dbl_mb(k_a + idima - 1)"
            newcode.statements.insert(newcode.pointer,newline)
            if (not doloopofa):
               newcode.setamark(1)
            newcode.pointer = newcode.pointer + 1

      newcode.pointer = newcode.getamark(1) + 1
      newline = string.join(["IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1

      # allocate sorted and unsorted tensor 2
      if (three):
         newcode.add("integers","l_b_sort")
         newcode.add("integers","k_b_sort")
         newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL ERRQUIT('",\
                                subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1
         newcode.add("integers","l_b")
         newcode.add("integers","k_b")
         newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('",\
                                subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1
   
         # mapping to a permutation symmetry unique block
         if (self.tensors[2].type in types[3]):
            superpermutations = restrictedpermutationwithparity(superlocaltwo,supercommontwo,[])
         else:
            superpermutations = restrictedpermutationwithparity(superglobaltwo,superlocaltwo,supercommontwo)
         if (self.tensors[2].type in types[3]):
            subpermutations = restrictedpermutationwithparity(sublocaltwo,subcommontwo,[])
         else:
            subpermutations = restrictedpermutationwithparity(subglobaltwo,sublocaltwo,subcommontwo)
         newcode.pointer = newcode.pointer - 1
         newcode.setamark(2)
         newcode.pointer = newcode.pointer + 1
         ifblock = 0
         for superpermutation in superpermutations:
            superline = ""
            if (self.tensors[2].type in types[3]):
               if (superpermutation[1] == "empty"):
                  superpermutedindexes = superglobaltwo + sortindexes(superlocaltwo + supercommontwo)
                  superfactor = 1
               else:
                  superpermutedindexes = superglobaltwo + superpermutation[1:]
                  superfactor = parityofpermutation(superglobaltwo + sortindexes(superlocaltwo + supercommontwo) + \
                                                    superglobaltwo + superpermutation[1:])
            else:
               if (superpermutation[1] == "empty"):
                  superpermutedindexes = self.tensors[2].indexes[0:len(self.tensors[2].indexes)/2]
                  superfactor = 1
               else:
                  superpermutedindexes = superpermutation[1:]
                  superfactor = parityofpermutation(self.tensors[2].indexes[0:len(self.tensors[2].indexes)/2] + \
                                                    superpermutation[1:])
            if (superpermutation[1] != "empty"):
               for nindex in range(len(superpermutedindexes)-1):
                  indexa = superpermutedindexes[nindex]
                  indexb = superpermutedindexes[nindex+1]
                  if ((self.tensors[2].type in types[3]) and (indexa.isin(superglobaltwo) or indexb.isin(superglobaltwo))):
                     continue
                  if (indexa.isin(superglobaltwo) and indexb.isin(superglobaltwo)):
                     continue
                  if (indexa.isin(superlocaltwo) and indexb.isin(superlocaltwo)):
                     continue
                  if (indexa.isin(supercommontwo) and indexb.isin(supercommontwo)):
                     continue
                  if (indexa.isgreaterthan(indexb)):
                     inequality = " .lt. "
                  else:
                     inequality = " .le. "
                  if (superline):
                     superline = string.join([superline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
                  else:
                     superline = string.join([superline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
            for subpermutation in subpermutations:
               subline = superline
               if (self.tensors[2].type in types[3]):
                  if (subpermutation[1] == "empty"):
                     subpermutedindexes = subglobaltwo + sortindexes(sublocaltwo + subcommontwo)
                     subfactor = 1
                  else:
                     subpermutedindexes = subglobaltwo + subpermutation[1:]
                     subfactor = parityofpermutation(subglobaltwo + sortindexes(sublocaltwo + subcommontwo) + \
                                                     subglobaltwo + subpermutation[1:])
               else:
                  if (subpermutation[1] == "empty"):
                     subpermutedindexes = self.tensors[2].indexes[len(self.tensors[2].indexes)/2:len(self.tensors[2].indexes)]
                     subfactor = 1
                  else:
                     subpermutedindexes = subpermutation[1:]
                     subfactor = parityofpermutation(self.tensors[2].indexes[len(self.tensors[2].indexes)/2:len(self.tensors[2].indexes)] + \
                                                     subpermutation[1:])
               if (subpermutation[1] != "empty"):
                  for nindex in range(len(subpermutedindexes)-1):
                     indexa = subpermutedindexes[nindex]
                     indexb = subpermutedindexes[nindex+1]
                     if ((self.tensors[2].type in types[3]) and (indexa.isin(subglobaltwo) or indexb.isin(subglobaltwo))):
                        continue
                     if (indexa.isin(subglobaltwo) and indexb.isin(subglobaltwo)):
                        continue
                     if (indexa.isin(sublocaltwo) and indexb.isin(sublocaltwo)):
                        continue
                     if (indexa.isin(subcommontwo) and indexb.isin(subcommontwo)):
                        continue
                     if (indexa.isgreaterthan(indexb)):
                        inequality = " .lt. "
                     else:
                        inequality = " .le. "
                     if (subline):
                        subline = string.join([subline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
                     else:
                        subline = string.join([subline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
               if (subline):
                  subline = string.join([subline,") THEN"],"")
                  newcode.pointer = newcode.getamark(2) + 1
                  newcode.statements.insert(newcode.pointer,subline)
                  newcode.pointer = newcode.pointer + 1
                  newline = "END IF"
                  newcode.statements.insert(newcode.pointer,newline)
                  newcode.setamark(2)
                  ifblock = 1
   
               permutedindexes = superpermutedindexes + subpermutedindexes

               # get a block
               arguments = ""
               argumentsend = ""
               for nindex in range(len(permutedindexes)-1,-1,-1):
                  if (permutedindexes[nindex].type == "hole"):
                     boffset = "b_2 - 1"
                  else:
                     if (self.tensors[2].type in types[4]):
                        boffset = "b_2 - 1"
                     else:
                        boffset = "b_2 - noab - 1"
                  if (arguments == ""):
                     arguments = string.join(["d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset + ",permutedindexes[nindex].show(),boffset],"")
                  else:
                     if (self.tensors[2].type in types[4]):
                        arguments = string.join([arguments," + (noab+nvab) * (",permutedindexes[nindex].show(),boffset],"")
                     else:
                        if (permutedindexes[nindex+1].type == "hole"):
                           arguments = string.join([arguments," + noab * (",permutedindexes[nindex].show(),boffset],"")
                        else:
                           arguments = string.join([arguments," + nvab * (",permutedindexes[nindex].show(),boffset],"")
                     argumentsend = string.join([argumentsend,")"],"")
               if (not arguments):
                  arguments = "d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset)"
               else:
                  arguments = string.join([arguments,argumentsend,")"],"")
               newline = string.join(["CALL GET_BLOCK(",arguments,")"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
   
               # sort indexes of tensor 2
               doloopofb = 0
               for nindex in range(len(permutedindexes)):
                  index = permutedindexes[nindex]
                  newint = index.show()
                  newcode.add("integers",newint)
                  newline = string.join(["DO ",newint," = 1,int_mb(k_range+",newint,"b-1)"],"")
                  newcode.statements.insert(newcode.pointer,newline)
                  newcode.pointer = newcode.pointer + 1
                  newline = "END DO"
                  newcode.statements.insert(newcode.pointer,newline)
                  doloopofb = 1
                  if ((not ifblock) and (nindex == 0)):
                     newcode.setamark(2)
               newline = ""
               newlineend = ""
               newcode.add("integers","idimb")
               for nindex in range(len(permutedindexes)-1,-1,-1):
                  if (newline == ""):
                     newline = string.join(["idimb = ",permutedindexes[nindex].show()],"")
                  else:
                     newline = string.join([newline," + int_mb(k_range+",permutedindexes[nindex+1].show(),"b-1)"\
                                            " * ((",permutedindexes[nindex].show()," - 1)"],"")
                     newlineend = string.join([newlineend,")"],"")
               newline = string.join([newline,newlineend],"")
               if (not newline):
                  newline = "idimb = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = ""
               newlineend = ""
               newcode.add("integers","idimb_sort")
               # note the sub - super order !
               sorted = subcommontwo + supercommontwo + superglobaltwo + superlocaltwo + subglobaltwo + sublocaltwo
               for nindex in range(len(sorted)):
                  if (newline == ""):
                     newline = string.join(["idimb_sort = ",sorted[nindex].show()],"")
                  else:
                     newline = string.join([newline," + int_mb(k_range+",sorted[nindex-1].show(),"b-1)"\
                                            " * ((",sorted[nindex].show()," - 1)"],"")
                     newlineend = string.join([newlineend,")"],"")
               newline = string.join([newline,newlineend],"")
               if (not newline):
                  newline = "idimb_sort = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               if (superfactor * subfactor == 1):
                  newline = "dbl_mb(k_b_sort + idimb_sort - 1) = dbl_mb(k_b + idimb - 1)"
               else:
                  newline = "dbl_mb(k_b_sort + idimb_sort - 1) = - dbl_mb(k_b + idimb - 1)"
               newcode.statements.insert(newcode.pointer,newline)
               if (not doloopofb):
                  newcode.setamark(2)
               newcode.pointer = newcode.pointer + 1

         newcode.pointer = newcode.getamark(2) + 1
         newline = string.join(["IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1

         # factor
# DEEXCITATION EXTENSION FROM HERE ...
         supercommonhole = []
         supercommonparticle = []
         for index in supercommon:
            if (index.ishole()):
               supercommonhole.append(index)
            elif (index.isparticle()):
               supercommonparticle.append(index)
            else:
               raise RuntimeError, "a general summation index appeared"
         subcommonhole = []
         subcommonparticle = []
         for index in subcommon:
            if (index.ishole()):
               subcommonhole.append(index)
            elif (index.isparticle()):
               subcommonparticle.append(index)
            else:
               raise RuntimeError, "a general summation index appeared"
         factorialforsuperhole = 0
         if (len(supercommonhole) > 1):
            factorialforsuperhole = 1
            newint = "nsuperh("+repr(len(supercommonhole))+")"
            newcode.add("integers",newint)
            newint = "isuperh"
            newcode.add("integers",newint)
            for isuperh in range(len(supercommonhole)):
               newline = "nsuperh("+repr(isuperh+1)+") = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "isuperh = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for nindex in range(len(supercommonhole)-1):
               newline = string.join(["IF (",supercommonhole[nindex].show(),"b .eq. ",supercommonhole[nindex+1].show(),"b) THEN"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "nsuperh(isuperh) = nsuperh(isuperh) + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "ELSE"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "isuperh = isuperh + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "END IF"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
         factorialforsuperparticle = 0
         if (len(supercommonparticle) > 1):
            factorialforsuperparticle = 1
            newint = "nsuperp("+repr(len(supercommonparticle))+")"
            newcode.add("integers",newint)
            newint = "isuperp"
            newcode.add("integers",newint)
            for isuperp in range(len(supercommonparticle)):
               newline = "nsuperp("+repr(isuperp+1)+") = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "isuperp = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for nindex in range(len(supercommonparticle)-1):
               newline = string.join(["IF (",supercommonparticle[nindex].show(),"b .eq. ",supercommonparticle[nindex+1].show(),"b) THEN"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "nsuperp(isuperp) = nsuperp(isuperp) + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "ELSE"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "isuperp = isuperp + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "END IF"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
         factorialforsubhole = 0
         if (len(subcommonhole) > 1):
            factorialforsubhole = 1
            newint = "nsubh("+repr(len(subcommonhole))+")"
            newcode.add("integers",newint)
            newint = "isubh"
            newcode.add("integers",newint)
            for isubh in range(len(subcommonhole)):
               newline = "nsubh("+repr(isubh+1)+") = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "isubh = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for nindex in range(len(subcommonhole)-1):
               newline = string.join(["IF (",subcommonhole[nindex].show(),"b .eq. ",subcommonhole[nindex+1].show(),"b) THEN"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "nsubh(isubh) = nsubh(isubh) + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "ELSE"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "isubh = isubh + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "END IF"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
         factorialforsubparticle = 0
         if (len(subcommonparticle) > 1):
            factorialforsubparticle = 1
            newint = "nsubp("+repr(len(subcommonparticle))+")"
            newcode.add("integers",newint)
            newint = "isubp"
            newcode.add("integers",newint)
            for isubp in range(len(subcommonparticle)):
               newline = "nsubp("+repr(isubp+1)+") = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "isubp = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for nindex in range(len(subcommonparticle)-1):
               newline = string.join(["IF (",subcommonparticle[nindex].show(),"b .eq. ",subcommonparticle[nindex+1].show(),"b) THEN"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "nsubp(isubp) = nsubp(isubp) + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "ELSE"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "isubp = isubp + 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "END IF"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
         if (factorialforsuperhole or factorialforsuperparticle or factorialforsubhole or factorialforsubparticle):
            newdbl = "FACTORIAL"
            newcode.add("doubles",newdbl)
            newcode.add("externals",newdbl)
            factor = repr(factorial(len(supercommonparticle))*factorial(len(supercommonhole))* \
                          factorial(len(subcommonparticle))*factorial(len(subcommonhole)))+".0d0"
            if (factorialforsuperhole):
               for isuperh in range(len(supercommonhole)):
                  factor = factor + "/FACTORIAL(nsuperh("+repr(isuperh+1)+"))"
            if (factorialforsuperparticle):
               for isuperp in range(len(supercommonparticle)):
                  factor = factor + "/FACTORIAL(nsuperp("+repr(isuperp+1)+"))"
            if (factorialforsubhole):
               for isubh in range(len(subcommonhole)):
                  factor = factor + "/FACTORIAL(nsubh("+repr(isubh+1)+"))"
            if (factorialforsubparticle):
               for isubp in range(len(subcommonparticle)):
                  factor = factor + "/FACTORIAL(nsubp("+repr(isubp+1)+"))"
         else:
            factor = "1.0d0"
# ... TO HERE
#        factorialforsuper = 0
#        if (len(supercommon) > 1):
#           factorialforsuper = 1
#           newint = "nsuper"
#           newcode.add("integers",newint)
#           newline = "nsuper = 1"
#           newcode.statements.insert(newcode.pointer,newline)
#           newcode.pointer = newcode.pointer + 1
#           for nindex in range(len(supercommon)-1):
#              newline = string.join(["IF (",supercommon[nindex].show(),"b .ne. ",supercommon[nindex+1].show(),\
#                                     "b) nsuper = nsuper + 1"],"")
#              newcode.statements.insert(newcode.pointer,newline)
#              newcode.pointer = newcode.pointer + 1
#        factorialforsub = 0
#        if (len(subcommon) > 1):
#           factorialforsub = 1
#           newint = "nsub"
#           newcode.add("integers",newint)
#           newline = "nsub = 1"
#           newcode.statements.insert(newcode.pointer,newline)
#           newcode.pointer = newcode.pointer + 1
#           for nindex in range(len(subcommon)-1):
#              newline = string.join(["IF (",subcommon[nindex].show(),"b .ne. ",subcommon[nindex+1].show(),\
#                                     "b) nsub = nsub + 1"],"")
#              newcode.statements.insert(newcode.pointer,newline)
#              newcode.pointer = newcode.pointer + 1
#        if (factorialforsuper and factorialforsub):
#           newdbl = "FACTORIAL"
#           newcode.add("doubles",newdbl)
#           newcode.add("externals",newdbl)
#           factor = "FACTORIAL(nsuper)*FACTORIAL(nsub)"
#        elif (factorialforsuper):
#           newdbl = "FACTORIAL"
#           newcode.add("doubles",newdbl)
#           newcode.add("externals",newdbl)
#           factor = "FACTORIAL(nsuper)"
#        elif (factorialforsub):
#           newdbl = "FACTORIAL"
#           newcode.add("doubles",newdbl)
#           newcode.add("externals",newdbl)
#           factor = "FACTORIAL(nsub)"
#        else:
#           factor = "1.0d0"
   
         # perform contraction and store the result
         newcode.add("integers","l_c_sort")
         newcode.add("integers","k_c_sort")
         if (three):
            newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dima_sort*dimb_sort,'noname',l_c_sort,k_c_sort)) CALL ERRQUIT('",\
                                   subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         else:
            newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dima_sort,'noname',l_c_sort,k_c_sort)) CALL ERRQUIT('",\
                                   subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1
         newline = string.join(["CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,",factor,\
                                ",dbl_mb(k_a_sort),dim_common,dbl_mb(k_b_sort),dim_common,0.0d0,dbl_mb(k_c_sort),dima_sort)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      # create an MA for unsorted tensor 0
      newcode.add("integers","l_c")
      newcode.add("integers","k_c")
      if (three):
         newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dima_sort*dimb_sort,'noname',l_c,k_c)) CALL ERRQUIT('",\
                                subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      else:
         newline = string.join(["IF (.not.MA_PUSH_GET(mt_dbl,dima_sort,'noname',l_c,k_c)) CALL ERRQUIT('",\
                                subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1
         
      # mapping to a permutation symmetry unique block
# DEEXCITATION EXTENSION FROM HERE ...
      superlocalparticleone  = []
      superlocalparticletwo  = []
      for index in superlocalone:
         if (index.isparticle()):
            superlocalparticleone.append(index)
      for index in superlocaltwo:
         if (index.isparticle()):
            superlocalparticletwo.append(index)
      superlocalparticlezero = sortindexes( superlocalparticleone + superlocalparticletwo )
      sublocalholeone  = []
      sublocalholetwo  = []
      for index in sublocalone:
         if (index.ishole()):
            sublocalholeone.append(index)
      for index in sublocaltwo:
         if (index.ishole()):
            sublocalholetwo.append(index)
      sublocalholezero = sortindexes( sublocalholeone + sublocalholetwo )
      deexcitationfactor = Factor([1.0],[[]])
      if ((superlocalparticleone and superlocalparticletwo) or \
          (sublocalholeone and sublocalholetwo)):
         deexcitationfactor = permutationsoffoursets(deexcitationfactor,superlocalparticleone, superlocalparticletwo, \
                                                                        sublocalholeone,       sublocalholetwo, \
                                                                        superlocalparticlezero,sublocalholezero)
         deexcitationfactor = deexcitationfactor.normalize()
# ... TO HERE
      newcode.pointer = newcode.pointer - 1
      newcode.setamark(3)
      newcode.pointer = newcode.pointer + 1
      ifblock = 0
# DEEXCITATION EXTENSION FROM HERE ...
      if ((len(deexcitationfactor.permutations) > 1) and (len(self.factor.permutations) > 1)):
         raise RuntimeError, "A logical error in code generator regarding deexcitation operator"
      if (len(deexcitationfactor.permutations) > 1):
         deexcitationfactor.multiply(self.factor.coefficients[0])
         currentfactor = deexcitationfactor.duplicate()
      else:
         currentfactor = self.factor.duplicate()
      for npermutation in range(len(currentfactor.permutations)):
         permutation = currentfactor.permutations[npermutation]
# ... TO HERE
         indexesintheoriginalorder = copy.deepcopy(superglobalzero + superlocalzero + subglobalzero + sublocalzero)
         permutedindexes = performpermutation(indexesintheoriginalorder,permutation,1)
         newline = ""
         for nindex in range(len(permutedindexes)-1):
            # no IF across super and sub indexes
            if (nindex == len(permutedindexes)/2-1):
               continue
            indexa = permutedindexes[nindex]
            indexb = permutedindexes[nindex+1]
            if (indexa.isin(superglobalone) and indexb.isin(superglobalone)):
               continue
            if (indexa.isin(superglobaltwo) and indexb.isin(superglobaltwo)):
               continue
            if (indexa.isin(superlocalone) and indexb.isin(superlocalone)):
               continue
            if (indexa.isin(superlocaltwo) and indexb.isin(superlocaltwo)):
               continue
            if (indexa.isin(subglobalone) and indexb.isin(subglobalone)):
               continue
            if (indexa.isin(subglobaltwo) and indexb.isin(subglobaltwo)):
               continue
            if (indexa.isin(sublocalone) and indexb.isin(sublocalone)):
               continue
            if (indexa.isin(sublocaltwo) and indexb.isin(sublocaltwo)):
               continue
            if (indexa.isin(superglobalone) and indexb.isin(superlocalone)):
               continue
            if (indexa.isin(superglobalone) and indexb.isin(superlocaltwo)):
               continue
            if (indexa.isin(superglobaltwo) and indexb.isin(superlocalone)):
               continue
            if (indexa.isin(superglobaltwo) and indexb.isin(superlocaltwo)):
               continue
            if (indexa.isin(superlocalone) and indexb.isin(superglobalone)):
               continue
            if (indexa.isin(superlocalone) and indexb.isin(superglobaltwo)):
               continue
            if (indexa.isin(superlocaltwo) and indexb.isin(superglobalone)):
               continue
            if (indexa.isin(superlocaltwo) and indexb.isin(superglobaltwo)):
               continue
            if (indexa.isin(subglobalone) and indexb.isin(sublocalone)):
               continue
            if (indexa.isin(subglobalone) and indexb.isin(sublocaltwo)):
               continue
            if (indexa.isin(subglobaltwo) and indexb.isin(sublocalone)):
               continue
            if (indexa.isin(subglobaltwo) and indexb.isin(sublocaltwo)):
               continue
            if (indexa.isin(sublocalone) and indexb.isin(subglobalone)):
               continue
            if (indexa.isin(sublocalone) and indexb.isin(subglobaltwo)):
               continue
            if (indexa.isin(sublocaltwo) and indexb.isin(subglobalone)):
               continue
            if (indexa.isin(sublocaltwo) and indexb.isin(subglobaltwo)):
               continue
            if (indexa.type != indexb.type):
               continue
            if (verbose):
               print indexa, " .le. ",indexb
            inequality = " .le. "
            if (newline):
               newline = string.join([newline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
            else:
               newline = string.join(["IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
         if (newline):
            newline = string.join([newline,") THEN"],"")
            newcode.pointer = newcode.getamark(3) + 1
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            newline = "END IF"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.setamark(3)
            ifblock = 1

         # sort indexes of tensor 0
         doloop = 0
         for nindex in range(len(permutedindexes)):
            index = permutedindexes[nindex]
            newint = index.show()
            newcode.add("integers",newint)
            newline = string.join(["DO ",newint," = 1,int_mb(k_range+",newint,"b-1)"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            newline = "END DO"
            newcode.statements.insert(newcode.pointer,newline)
            if (nindex == 0):
               if (ifblock):
                  newcode.setamark(4)
                  doloop = 1
               else:
                  newcode.setamark(3)
                  doloop = 1
         newline = ""
         newlineend = ""
         if (three):
            newcode.add("integers","idimc_sort")
         sorted = superglobalone + superlocalone + subglobalone + sublocalone \
                + superglobaltwo + superlocaltwo + subglobaltwo + sublocaltwo
         for nindex in range(len(sorted)):
            if (newline == ""):
               if (three):
                  newline = string.join(["idimc_sort = ",sorted[nindex].show()],"")
               else:
                  newline = string.join(["idima_sort = ",sorted[nindex].show()],"")
            else:
               newline = string.join([newline," + int_mb(k_range+",sorted[nindex-1].show(),"b-1)"\
                                      " * ((",sorted[nindex].show()," - 1)"],"")
               newlineend = string.join([newlineend,")"],"")
         newline = string.join([newline,newlineend],"")
         if (not newline):
            if (three):
               newline = "idimc_sort = 1"
            else:
               newline = "idima_sort = 1"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newline = ""
         newlineend = ""
         newcode.add("integers","idimc")
         for nindex in range(len(permutedindexes)-1,-1,-1):
            if (newline == ""):
               newline = string.join(["idimc = ",permutedindexes[nindex].show()],"")
            else:
               newline = string.join([newline," + int_mb(k_range+",permutedindexes[nindex+1].show(),"b-1)"\
                                      " * ((",permutedindexes[nindex].show()," - 1)"],"")
               newlineend = string.join([newlineend,")"],"")
         if (not newline):
            newline = "idimc = 1"
         newline = string.join([newline,newlineend],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         if (three):
# DEEXCITATION EXTENSION FROM HERE ...
            if (currentfactor.coefficients[npermutation] == 1.0):
               newline = "dbl_mb(k_c + idimc - 1) = dbl_mb(k_c_sort + idimc_sort - 1)"
            elif (currentfactor.coefficients[npermutation] == - 1.0):
               newline = "dbl_mb(k_c + idimc - 1) = - dbl_mb(k_c_sort + idimc_sort - 1)"
            else:
               newline = string.join(["dbl_mb(k_c + idimc - 1) = ",\
                                      repr(rationaltofractional(currentfactor.coefficients[npermutation])[0]),".0d0/",\
                                      repr(rationaltofractional(currentfactor.coefficients[npermutation])[1]),\
                                      ".0d0 * dbl_mb(k_c_sort + idimc_sort - 1)"],"")
         else:
            if (currentfactor.coefficients[npermutation] == 1.0):
               newline = "dbl_mb(k_c + idimc - 1) = dbl_mb(k_a_sort + idima_sort - 1)"
            elif (currentfactor.coefficients[npermutation] == - 1.0):
               newline = "dbl_mb(k_c + idimc - 1) = - dbl_mb(k_a_sort + idima_sort - 1)"
            else:
               newline = string.join(["dbl_mb(k_c + idimc - 1) = ",\
                                      repr(rationaltofractional(currentfactor.coefficients[npermutation])[0]),".0d0/",\
                                      repr(rationaltofractional(currentfactor.coefficients[npermutation])[1]),\
                                      ".0d0 * dbl_mb(k_a_sort + idima_sort - 1)"],"")
# ... TO HERE
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

         # accumulate a block
         if (ifblock):
            newcode.pointer = newcode.getamark(4) + 1
         elif (doloop):
            newcode.pointer = newcode.getamark(3) + 1
         arguments = ""
         argumentsend = ""
         for nindex in range(len(permutedindexes)-1,-1,-1):
            if (permutedindexes[nindex].type == "hole"):
               boffset = "b - 1"
            else:
               boffset = "b - noab - 1"
            if (arguments == ""):
               if (three):
                  arguments = string.join(["d_c,dbl_mb(k_c),dima_sort*dimb_sort,int_mb(k_c_offset + ",permutedindexes[nindex].show(),boffset],"")
               else:
                  arguments = string.join(["d_c,dbl_mb(k_c),dima_sort,int_mb(k_c_offset + ",permutedindexes[nindex].show(),boffset],"")
            else:
               if (permutedindexes[nindex+1].type == "hole"):
                  arguments = string.join([arguments," + noab * (",permutedindexes[nindex].show(),boffset],"")
               else:
                  arguments = string.join([arguments," + nvab * (",permutedindexes[nindex].show(),boffset],"")
               argumentsend = string.join([argumentsend,")"],"")
         if (arguments == ""):
            if (three):
               arguments = "d_c,dbl_mb(k_c),dima_sort*dimb_sort,int_mb(k_c_offset"
            else:
               arguments = "d_c,dbl_mb(k_c),dima_sort,int_mb(k_c_offset"
            argumentsend = ""
         arguments = string.join([arguments,argumentsend,")"],"")
         newline = string.join(["CALL ADD_BLOCK(",arguments,")"],"")
         newcode.statements.insert(newcode.pointer,newline)
         if (not ifblock):
            newcode.setamark(3)
         newcode.pointer = newcode.pointer + 1

      newcode.pointer = newcode.getamark(3) + 1
      newline = string.join(["IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1
      if (three):
         newline = string.join(["IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1
         newline = string.join(["IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1
         newline = string.join(["IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1
      else:
         newline = string.join(["IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('",subroutinename,"',",repr(errquit),",MA_ERR)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         errquit = errquit + 1

      # close the subroutine
      newcode.pointer = len(newcode.statements)
      newline = "next = NXTASK(-nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "call GA_SYNC()"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      return newcode

   def utchem(self,globaltargetindexes,subroutinename="NONAME",relativistic=0):
      """Suggests an implementation in Fortran77 for an elementary tensor contraction C = A * B"""

      if (len(self.tensors) == 3):
         three = 1
      else:
         three = 0

      errquit = 0

      if (relativistic):
         newcode = Code("Fortran90","R4D_"+subroutinename)
      else:
         newcode = Code("Fortran90",subroutinename)

      # header
      newline = "!" + self.show(0)
      newcode.add("headers",newline)
      newline = "USE UT_SYS_MODULE"
      newcode.add("headers",newline)
      newline = "USE UT_MOLINP_MODULE"
      newcode.add("headers",newline)
      if (relativistic):
         newline = "USE UT_R4DTCE_MODULE"
      else:
         newline = "USE UT_TCE_MODULE"
      newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
      
      # insert include statements
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)

      # declaration
      newint = "d_a"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "a_offset"
      newcode.add("integerarrays",newint)
      newcode.add("arguments",newint)
      if (three):
         newint = "d_b"
         newcode.add("integers",newint)
         newcode.add("arguments",newint)
         newint = "b_offset"
         newcode.add("integerarrays",newint)
         newcode.add("arguments",newint)
      newint = "d_c"
      newcode.add("integers",newint)
      newcode.add("arguments",newint)
      newint = "c_offset"
      newcode.add("integerarrays",newint)
      newcode.add("arguments",newint)
      newint = "NXTASK"
      newcode.add("integers",newint)
      newcode.add("externals",newint)
      newint = "next"
      newcode.add("integers",newint)
      newint = "nprocs"
      newcode.add("integers",newint)
      newint = "count"
      newcode.add("integers",newint)
            
      # Tensor 0
      superglobalzero = []
      subglobalzero = []
      superlocalzero = []
      sublocalzero = []
      for nindex in range(len(self.tensors[0].indexes)/2):
         index = self.tensors[0].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            superglobalzero.append(index)
         else:
            superlocalzero.append(index)
      for nindex in range(len(self.tensors[0].indexes)/2,len(self.tensors[0].indexes)):
         index = self.tensors[0].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            subglobalzero.append(index)
         else:
            sublocalzero.append(index)

      # Tensor 1
      superglobalone = []
      subglobalone = []
      superlocalone = []
      sublocalone = []
      supercommonone = []
      subcommonone = []
      if (not self.tensors[1].conjugate):
         superrangeone = range(len(self.tensors[1].indexes)/2)
         subrangeone = range(len(self.tensors[1].indexes)/2,len(self.tensors[1].indexes))
      else:
         superrangeone = range(len(self.tensors[1].indexes)/2,len(self.tensors[1].indexes))
         subrangeone = range(len(self.tensors[1].indexes)/2)
      for nindex in superrangeone:
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            superglobalone.append(index)
         elif (self.summation):
            if (index.isin(self.summation.indexes)):
               supercommonone.append(index)
            else:
               superlocalone.append(index)
         else:
            superlocalone.append(index)
      for nindex in subrangeone:
         index = self.tensors[1].indexes[nindex]
         if (index.isin(globaltargetindexes)):
            subglobalone.append(index)
         elif (self.summation):
            if (index.isin(self.summation.indexes)):
               subcommonone.append(index)
            else:
               sublocalone.append(index)
         else:
            sublocalone.append(index)

      # Tensor 2
      superglobaltwo = []
      subglobaltwo = []
      superlocaltwo = []
      sublocaltwo = []
      supercommontwo = []
      subcommontwo = []
      if (three):
         if (not self.tensors[2].conjugate):
            superrangetwo = range(len(self.tensors[2].indexes)/2)
            subrangetwo = range(len(self.tensors[2].indexes)/2,len(self.tensors[2].indexes))
         else:
            superrangetwo = range(len(self.tensors[2].indexes)/2,len(self.tensors[2].indexes))
            subrangetwo = range(len(self.tensors[2].indexes)/2)
         for nindex in superrangetwo:
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               superglobaltwo.append(index)
            elif (self.summation):
               if (index.isin(self.summation.indexes)):
                  supercommontwo.append(index)
               else:
                  superlocaltwo.append(index)
            else:
               superlocaltwo.append(index)
         for nindex in subrangetwo:
            index = self.tensors[2].indexes[nindex]
            if (index.isin(globaltargetindexes)):
               subglobaltwo.append(index)
            elif (self.summation):
               if (index.isin(self.summation.indexes)):
                  subcommontwo.append(index)
               else:
                  sublocaltwo.append(index)
            else:
               sublocaltwo.append(index)
      if (len(supercommonone) > len(subcommontwo)):
         supercommon = supercommonone
      else:
         supercommon = subcommontwo
      if (len(subcommonone) > len(supercommontwo)):
         subcommon = subcommonone
      else:
         subcommon = supercommontwo

      # parallel related
      newline = "nprocs = GA_NNODES()"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "count = 0"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "next = NXTASK(nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # loop over output tensor indexes
      newcode.inserttileddoloops(superglobalone)
      if (three):
         newcode.inserttileddoloops(superglobaltwo)
      newcode.inserttileddoloops(superlocalone)
      if (three):
         newcode.inserttileddoloops(superlocaltwo)
      newcode.inserttileddoloops(subglobalone)
      if (three):
         newcode.inserttileddoloops(subglobaltwo)
      newcode.inserttileddoloops(sublocalone)
      if (three):
         newcode.inserttileddoloops(sublocaltwo)

      # parallel related
      newline = "IF (next == count) THEN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.setamark(4)
      newcode.pointer = newcode.pointer + 1
      newline = "next = NXTASK(nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END IF"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "count = count + 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.getamark(4) + 1

      # spin restriction on output tensor
      if (not relativistic):
         newcode.inserttiledifrestricted(superglobalzero + superlocalzero + subglobalzero + sublocalzero)

      # symmetry of output tensor
      super = superglobalzero + superlocalzero
      sub = subglobalzero + sublocalzero
      newcode.inserttiledifsymmetry(super,sub,self.tensors[0].irrep,relativistic)

      # loop over summation indexes
      newcode.inserttileddoloops(supercommon)
      newcode.inserttileddoloops(subcommon)

      # symmetry of input tensor 1
      if (three):
         super = superglobalone + superlocalone + supercommonone
         sub = subglobalone + sublocalone + subcommonone
         newcode.inserttiledifsymmetry(super,sub,self.tensors[1].irrep,relativistic)

      # spin restriction on input tensor one
      indexesone = superglobalone + superlocalone + supercommonone + subglobalone + sublocalone + subcommonone
      if (indexesone and (not relativistic)):
         newline = "IF ((restricted).and.("
         conjugation = ""
         for index in indexesone:
            newint = string.join([index.show(),"b"],"")
            newline = string.join([newline,conjugation,"spin(",newint,")"],"")
            conjugation = "+"
         newline = string.join([newline," == ",repr(2*len(indexesone)),")) THEN"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         for index in indexesone:
            if (relativistic):
               newint = string.join([index.show(),"b"],"")
            else:
               newint = string.join([index.show(),"b_1"],"")
            newcode.add("integers",newint)
            newline = string.join([newint," = alpha(",index.show(),"b)"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
         newline = "ELSE"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         for index in indexesone:
            if (relativistic):
               newint = string.join([index.show(),"b"],"")
            else:
               newint = string.join([index.show(),"b_1"],"")
            newcode.add("integers",newint)
            newline = string.join([newint," = ",index.show(),"b"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
         newline = "END IF"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      # spin restriction on input tensor one
      if (three):
         indexestwo = superglobaltwo + superlocaltwo + supercommontwo + subglobaltwo + sublocaltwo + subcommontwo
         if (indexestwo and (not relativistic)):
            newline = "IF ((restricted).and.("
            conjugation = ""
            for index in indexestwo:
               newint = string.join([index.show(),"b"],"")
               newline = string.join([newline,conjugation,"spin(",newint,")"],"")
               conjugation = "+"
            newline = string.join([newline," == ",repr(2*len(indexestwo)),")) THEN"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for index in indexestwo:
               if (relativistic):
                  newint = string.join([index.show(),"b"],"")
               else:
                  newint = string.join([index.show(),"b_2"],"")
               newcode.add("integers",newint)
               newline = string.join([newint," = alpha(",index.show(),"b)"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "ELSE"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for index in indexestwo:
               if (relativistic):
                  newint = string.join([index.show(),"b"],"")
               else:
                  newint = string.join([index.show(),"b_2"],"")
               newcode.add("integers",newint)
               newline = string.join([newint," = ",index.show(),"b"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
            newline = "END IF"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1

      # create MA's for tensor 1
      newcode.add("integers","dim_common")
      newline = ""
      for index in supercommonone + subcommonone:
         if (newline == ""):
            newline = string.join(["dim_common = range(",index.show(),"b)"],"")
         else:
            newline = string.join([newline," * range(",index.show(),"b)"],"")
      if (newline == ""):
         newline = "dim_common = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.add("integers","dima_sort")
      newline = ""
      for index in superglobalone + superlocalone + subglobalone + sublocalone:
         if (newline == ""):
            newline = string.join(["dima_sort = range(",index.show(),"b)"],"")
         else:
            newline = string.join([newline," * range(",index.show(),"b)"],"")
      if (newline == ""):
         newline = "dima_sort = 1"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newcode.add("integers","dima")
      newline = "dima = dim_common * dima_sort"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # create MA's for tensor 2
      if (three):
         newcode.add("integers","dimb_sort")
         newline = ""
         for index in superglobaltwo + superlocaltwo + subglobaltwo + sublocaltwo:
            if (newline == ""):
               newline = string.join(["dimb_sort = range(",index.show(),"b)"],"")
            else:
               newline = string.join([newline," * range(",index.show(),"b)"],"")
         if (newline == ""):
            newline = "dimb_sort = 1"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newcode.add("integers","dimb")
         newline = "dimb = dim_common * dimb_sort"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      if (three):
         newline = "IF ((dima > 0) .and. (dimb > 0)) THEN"
      else:
         newline = "IF (dima > 0) THEN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END IF"
      newcode.statements.insert(newcode.pointer,newline)

      # allocate sorted and unsorted tensor 1
      if (relativistic):
         newcode.add("doublecomplexallocatables","a_sort")
      else:
         newcode.add("doubleallocatables","a_sort")
      newline = "ALLOCATE(a_sort(dima))"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      if (relativistic):
         newcode.add("doublecomplexallocatables","a")
      else:
         newcode.add("doubleallocatables","a")
      newline = "ALLOCATE(a(dima))"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # mapping to a permutation symmetry unique block
      if ((self.tensors[1].type == "i") or (self.tensors[1].type == "j")):
         superpermutations = restrictedpermutationwithparity(superlocalone,supercommonone,[])
      else:
         superpermutations = restrictedpermutationwithparity(superglobalone,superlocalone,supercommonone)
      if ((self.tensors[1].type == "i") or (self.tensors[1].type == "j")):
         subpermutations = restrictedpermutationwithparity(sublocalone,subcommonone,[])
      else:
         subpermutations = restrictedpermutationwithparity(subglobalone,sublocalone,subcommonone)
      newcode.pointer = newcode.pointer - 1
      newcode.setamark(1)
      newcode.pointer = newcode.pointer + 1
      ifblock = 0
      for superpermutation in superpermutations:
         superline = ""
         if ((self.tensors[1].type == "i") or (self.tensors[1].type == "j")):
            if (superpermutation[1] == "empty"):
               superpermutedindexes = superglobalone + sortindexes(superlocalone + supercommonone)
               superfactor = 1
            else:
               superpermutedindexes = superglobalone + superpermutation[1:]
               superfactor = parityofpermutation(superglobalone + sortindexes(superlocalone + supercommonone) + \
                                                 superglobalone + superpermutation[1:])
         else:
            if (superpermutation[1] == "empty"):
               superpermutedindexes = self.tensors[1].indexes[0:len(self.tensors[1].indexes)/2]
               superfactor = 1
            else:
               superpermutedindexes = superpermutation[1:]
               superfactor = parityofpermutation(self.tensors[1].indexes[0:len(self.tensors[1].indexes)/2] + \
                                                 superpermutation[1:])
         if (superpermutation[1] != "empty"):
            for nindex in range(len(superpermutedindexes)-1):
               indexa = superpermutedindexes[nindex]
               indexb = superpermutedindexes[nindex+1]
               if (((self.tensors[1].type == "i") or (self.tensors[1].type == "j")) and (indexa.isin(superglobalone) or indexb.isin(superglobalone))):
                  continue
               if (indexa.isin(superglobalone) and indexb.isin(superglobalone)):
                  continue
               if (indexa.isin(superlocalone) and indexb.isin(superlocalone)):
                  continue
               if (indexa.isin(supercommonone) and indexb.isin(supercommonone)):
                  continue
               if (indexa.isgreaterthan(indexb)):
                  inequality = " < "
               else:
                  inequality = " <= "
               if (superline):
                  superline = string.join([superline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
               else:
                  superline = string.join([superline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
         for subpermutation in subpermutations:
            subline = superline
            if ((self.tensors[1].type == "i") or (self.tensors[1].type == "j")):
               if (subpermutation[1] == "empty"):
                  subpermutedindexes = subglobalone + sortindexes(sublocalone + subcommonone)
                  subfactor = 1
               else:
                  subpermutedindexes = subglobalone + subpermutation[1:]
                  subfactor = parityofpermutation(subglobalone + sortindexes(sublocalone + subcommonone) + \
                                                  subglobalone + subpermutation[1:])
            else:
               if (subpermutation[1] == "empty"):
                  subpermutedindexes = self.tensors[1].indexes[len(self.tensors[1].indexes)/2:len(self.tensors[1].indexes)]
                  subfactor = 1
               else:
                  subpermutedindexes = subpermutation[1:]
                  subfactor = parityofpermutation(self.tensors[1].indexes[len(self.tensors[1].indexes)/2:len(self.tensors[1].indexes)] + \
                                                  subpermutation[1:])
            if (subpermutation[1] != "empty"):
               for nindex in range(len(subpermutedindexes)-1):
                  indexa = subpermutedindexes[nindex]
                  indexb = subpermutedindexes[nindex+1]
                  if (((self.tensors[1].type == "i") or (self.tensors[1].type == "j")) and (indexa.isin(subglobalone) or indexb.isin(subglobalone))):
                     continue
                  if (indexa.isin(subglobalone) and indexb.isin(subglobalone)):
                     continue
                  if (indexa.isin(sublocalone) and indexb.isin(sublocalone)):
                     continue
                  if (indexa.isin(subcommonone) and indexb.isin(subcommonone)):
                     continue
                  if (indexa.isgreaterthan(indexb)):
                     inequality = " < "
                  else:
                     inequality = " <= "
                  if (subline):
                     subline = string.join([subline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
                  else:
                     subline = string.join([subline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
            if (subline):
               subline = string.join([subline,") THEN"],"")
               newcode.pointer = newcode.getamark(1) + 1
               newcode.statements.insert(newcode.pointer,subline)
               newcode.pointer = newcode.pointer + 1
               newline = "END IF"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.setamark(1)
               ifblock = 1

            permutedindexes = superpermutedindexes + subpermutedindexes

            # get a block
            arguments = ""
            argumentsend = ""
            for nindex in range(len(permutedindexes)-1,-1,-1):
               if (permutedindexes[nindex].type == "hole"):
                  if (relativistic):
                     boffset = "b - 1"
                  else:
                     boffset = "b_1 - 1"
               else:
                  if ((self.tensors[1].type == "f") or (self.tensors[1].type == "v")):
                     if (relativistic):
                        boffset = "b - 1"
                     else:
                        boffset = "b_1 - 1"
                  else:
                     if (relativistic):
                        boffset = "b - noab - 1"
                     else:
                        boffset = "b_1 - noab - 1"
               if (arguments == ""):
                  arguments = string.join(["d_a,a,dima,a_offset(1 + ",permutedindexes[nindex].show(),boffset],"")
               else:
                  if ((self.tensors[1].type == "f") or (self.tensors[1].type == "v")):
                     arguments = string.join([arguments," + (noab+nvab) * (",permutedindexes[nindex].show(),boffset],"")
                  else:
                     if (permutedindexes[nindex+1].type == "hole"):
                        arguments = string.join([arguments," + noab * (",permutedindexes[nindex].show(),boffset],"")
                     else:
                        arguments = string.join([arguments," + nvab * (",permutedindexes[nindex].show(),boffset],"")
                  argumentsend = string.join([argumentsend,")"],"")
            if (not arguments):
               arguments = "d_a,a,dima,a_offset"
            else:
               arguments = string.join([arguments,argumentsend,")"],"")
            if (relativistic):
               newline = string.join(["CALL R4D_GET_BLOCK(",arguments,")"],"")
            else:
               newline = string.join(["CALL GET_BLOCK(",arguments,")"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1

            # sort indexes of tensor 1
            doloopofa = 0
            for nindex in range(len(permutedindexes)):
               index = permutedindexes[nindex]
               newint = index.show()
               newcode.add("integers",newint)
               newline = string.join(["DO ",newint," = 1,range(",newint,"b)"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = "END DO"
               newcode.statements.insert(newcode.pointer,newline)
               doloopofa = 1
               if ((not ifblock) and (nindex == 0)):
                  newcode.setamark(1)
            newline = ""
            newlineend = ""
            newcode.add("integers","idima")
            for nindex in range(len(permutedindexes)-1,-1,-1):
               if (newline == ""):
                  newline = string.join(["idima = ",permutedindexes[nindex].show()],"")
               else:
                  newline = string.join([newline," + range(",permutedindexes[nindex+1].show(),"b)"\
                                         " * ((",permutedindexes[nindex].show()," - 1)"],"")
                  newlineend = string.join([newlineend,")"],"")
            newline = string.join([newline,newlineend],"")
            if (not newline):
               newline = "idima = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            newline = ""
            newlineend = ""
            newcode.add("integers","idima_sort")
            sorted = supercommonone + subcommonone + superglobalone + superlocalone + subglobalone + sublocalone
            for nindex in range(len(sorted)):
               if (newline == ""):
                  newline = string.join(["idima_sort = ",sorted[nindex].show()],"")
               else:
                  newline = string.join([newline," + range(",sorted[nindex-1].show(),"b)"\
                                         " * ((",sorted[nindex].show()," - 1)"],"")
                  newlineend = string.join([newlineend,")"],"")
            newline = string.join([newline,newlineend],"")
            if (not newline):
               newline = "idima_sort = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            if (superfactor * subfactor == 1):
               newline = "a_sort(idima_sort) = a(idima)"
            else:
               newline = "a_sort(idima_sort) = - a(idima)"
            newcode.statements.insert(newcode.pointer,newline)
            if (not doloopofa):
               newcode.setamark(1)
            newcode.pointer = newcode.pointer + 1

      newcode.pointer = newcode.getamark(1) + 1
      newline = "DEALLOCATE(a)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      errquit = errquit + 1

      # allocate sorted and unsorted tensor 2
      if (three):
         if (relativistic):
            newcode.add("doublecomplexallocatables","b_sort")
         else:
            newcode.add("doubleallocatables","b_sort")
         newline = "ALLOCATE(b_sort(dimb))"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         if (relativistic):
            newcode.add("doublecomplexallocatables","b")
         else:
            newcode.add("doubleallocatables","b")
         newline = "ALLOCATE(b(dimb))"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
   
         # mapping to a permutation symmetry unique block
         if ((self.tensors[2].type == "i") or (self.tensors[2].type == "j")):
            superpermutations = restrictedpermutationwithparity(superlocaltwo,supercommontwo,[])
         else:
            superpermutations = restrictedpermutationwithparity(superglobaltwo,superlocaltwo,supercommontwo)
         if ((self.tensors[2].type == "i") or (self.tensors[2].type == "j")):
            subpermutations = restrictedpermutationwithparity(sublocaltwo,subcommontwo,[])
         else:
            subpermutations = restrictedpermutationwithparity(subglobaltwo,sublocaltwo,subcommontwo)
         newcode.pointer = newcode.pointer - 1
         newcode.setamark(2)
         newcode.pointer = newcode.pointer + 1
         ifblock = 0
         for superpermutation in superpermutations:
            superline = ""
            if ((self.tensors[2].type == "i") or (self.tensors[2].type == "j")):
               if (superpermutation[1] == "empty"):
                  superpermutedindexes = superglobaltwo + sortindexes(superlocaltwo + supercommontwo)
                  superfactor = 1
               else:
                  superpermutedindexes = superglobaltwo + superpermutation[1:]
                  superfactor = parityofpermutation(superglobaltwo + sortindexes(superlocaltwo + supercommontwo) + \
                                                    superglobaltwo + superpermutation[1:])
            else:
               if (superpermutation[1] == "empty"):
                  superpermutedindexes = self.tensors[2].indexes[0:len(self.tensors[2].indexes)/2]
                  superfactor = 1
               else:
                  superpermutedindexes = superpermutation[1:]
                  superfactor = parityofpermutation(self.tensors[2].indexes[0:len(self.tensors[2].indexes)/2] + \
                                                    superpermutation[1:])
            if (superpermutation[1] != "empty"):
               for nindex in range(len(superpermutedindexes)-1):
                  indexa = superpermutedindexes[nindex]
                  indexb = superpermutedindexes[nindex+1]
                  if (((self.tensors[2].type == "i") or (self.tensors[2].type == "j")) and (indexa.isin(superglobaltwo) or indexb.isin(superglobaltwo))):
                     continue
                  if (indexa.isin(superglobaltwo) and indexb.isin(superglobaltwo)):
                     continue
                  if (indexa.isin(superlocaltwo) and indexb.isin(superlocaltwo)):
                     continue
                  if (indexa.isin(supercommontwo) and indexb.isin(supercommontwo)):
                     continue
                  if (indexa.isgreaterthan(indexb)):
                     inequality = " < "
                  else:
                     inequality = " <= "
                  if (superline):
                     superline = string.join([superline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
                  else:
                     superline = string.join([superline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
            for subpermutation in subpermutations:
               subline = superline
               if ((self.tensors[2].type == "i") or (self.tensors[2].type == "j")):
                  if (subpermutation[1] == "empty"):
                     subpermutedindexes = subglobaltwo + sortindexes(sublocaltwo + subcommontwo)
                     subfactor = 1
                  else:
                     subpermutedindexes = subglobaltwo + subpermutation[1:]
                     subfactor = parityofpermutation(subglobaltwo + sortindexes(sublocaltwo + subcommontwo) + \
                                                     subglobaltwo + subpermutation[1:])
               else:
                  if (subpermutation[1] == "empty"):
                     subpermutedindexes = self.tensors[2].indexes[len(self.tensors[2].indexes)/2:len(self.tensors[2].indexes)]
                     subfactor = 1
                  else:
                     subpermutedindexes = subpermutation[1:]
                     subfactor = parityofpermutation(self.tensors[2].indexes[len(self.tensors[2].indexes)/2:len(self.tensors[2].indexes)] + \
                                                     subpermutation[1:])
               if (subpermutation[1] != "empty"):
                  for nindex in range(len(subpermutedindexes)-1):
                     indexa = subpermutedindexes[nindex]
                     indexb = subpermutedindexes[nindex+1]
                     if (((self.tensors[2].type == "i") or (self.tensors[2].type == "j")) and (indexa.isin(subglobaltwo) or indexb.isin(subglobaltwo))):
                        continue
                     if (indexa.isin(subglobaltwo) and indexb.isin(subglobaltwo)):
                        continue
                     if (indexa.isin(sublocaltwo) and indexb.isin(sublocaltwo)):
                        continue
                     if (indexa.isin(subcommontwo) and indexb.isin(subcommontwo)):
                        continue
                     if (indexa.isgreaterthan(indexb)):
                        inequality = " < "
                     else:
                        inequality = " <= "
                     if (subline):
                        subline = string.join([subline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
                     else:
                        subline = string.join([subline,"IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
               if (subline):
                  subline = string.join([subline,") THEN"],"")
                  newcode.pointer = newcode.getamark(2) + 1
                  newcode.statements.insert(newcode.pointer,subline)
                  newcode.pointer = newcode.pointer + 1
                  newline = "END IF"
                  newcode.statements.insert(newcode.pointer,newline)
                  newcode.setamark(2)
                  ifblock = 1
   
               permutedindexes = superpermutedindexes + subpermutedindexes

               # get a block
               arguments = ""
               argumentsend = ""
               for nindex in range(len(permutedindexes)-1,-1,-1):
                  if (permutedindexes[nindex].type == "hole"):
                     if (relativistic):
                        boffset = "b - 1"
                     else:
                        boffset = "b_2 - 1"
                  else:
                     if ((self.tensors[2].type == "f") or (self.tensors[2].type == "v")):
                        if (relativistic):
                           boffset = "b - 1"
                        else:
                           boffset = "b_2 - 1"
                     else:
                        if (relativistic):
                           boffset = "b - noab - 1"
                        else:
                           boffset = "b_2 - noab - 1"
                  if (arguments == ""):
                     arguments = string.join(["d_b,b,dimb,b_offset(1 + ",permutedindexes[nindex].show(),boffset],"")
                  else:
                     if ((self.tensors[2].type == "f") or (self.tensors[2].type == "v")):
                        arguments = string.join([arguments," + (noab+nvab) * (",permutedindexes[nindex].show(),boffset],"")
                     else:
                        if (permutedindexes[nindex+1].type == "hole"):
                           arguments = string.join([arguments," + noab * (",permutedindexes[nindex].show(),boffset],"")
                        else:
                           arguments = string.join([arguments," + nvab * (",permutedindexes[nindex].show(),boffset],"")
                     argumentsend = string.join([argumentsend,")"],"")
               if (not arguments):
                  arguments = "d_b,b,dimb,b_offset"
               else:
                  arguments = string.join([arguments,argumentsend,")"],"")
               if (relativistic):
                  newline = string.join(["CALL R4D_GET_BLOCK(",arguments,")"],"")
               else:
                  newline = string.join(["CALL GET_BLOCK(",arguments,")"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
   
               # sort indexes of tensor 2
               doloopofb = 0
               for nindex in range(len(permutedindexes)):
                  index = permutedindexes[nindex]
                  newint = index.show()
                  newcode.add("integers",newint)
                  newline = string.join(["DO ",newint," = 1,range(",newint,"b)"],"")
                  newcode.statements.insert(newcode.pointer,newline)
                  newcode.pointer = newcode.pointer + 1
                  newline = "END DO"
                  newcode.statements.insert(newcode.pointer,newline)
                  doloopofb = 1
                  if ((not ifblock) and (nindex == 0)):
                     newcode.setamark(2)
               newline = ""
               newlineend = ""
               newcode.add("integers","idimb")
               for nindex in range(len(permutedindexes)-1,-1,-1):
                  if (newline == ""):
                     newline = string.join(["idimb = ",permutedindexes[nindex].show()],"")
                  else:
                     newline = string.join([newline," + range(",permutedindexes[nindex+1].show(),"b)"\
                                            " * ((",permutedindexes[nindex].show()," - 1)"],"")
                     newlineend = string.join([newlineend,")"],"")
               newline = string.join([newline,newlineend],"")
               if (not newline):
                  newline = "idimb = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newline = ""
               newlineend = ""
               newcode.add("integers","idimb_sort")
               # note the sub - super order !
               sorted = subcommontwo + supercommontwo + superglobaltwo + superlocaltwo + subglobaltwo + sublocaltwo
               for nindex in range(len(sorted)):
                  if (newline == ""):
                     newline = string.join(["idimb_sort = ",sorted[nindex].show()],"")
                  else:
                     newline = string.join([newline," + range(",sorted[nindex-1].show(),"b)"\
                                            " * ((",sorted[nindex].show()," - 1)"],"")
                     newlineend = string.join([newlineend,")"],"")
               newline = string.join([newline,newlineend],"")
               if (not newline):
                  newline = "idimb_sort = 1"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               if (superfactor * subfactor == 1):
                  newline = "b_sort(idimb_sort) = b(idimb)"
               else:
                  newline = "b_sort(idimb_sort) = - b(idimb)"
               newcode.statements.insert(newcode.pointer,newline)
               if (not doloopofb):
                  newcode.setamark(2)
               newcode.pointer = newcode.pointer + 1

         newcode.pointer = newcode.getamark(2) + 1
         newline = "DEALLOCATE(b)"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

         # factor
         factorialforsuper = 0
         if (len(supercommon) > 1):
            factorialforsuper = 1
            newint = "nsuper"
            newcode.add("integers",newint)
            newline = "nsuper = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for nindex in range(len(supercommon)-1):
               newline = string.join(["IF (",supercommon[nindex].show(),"b /= ",supercommon[nindex+1].show(),\
                                      "b) nsuper = nsuper + 1"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
         factorialforsub = 0
         if (len(subcommon) > 1):
            factorialforsub = 1
            newint = "nsub"
            newcode.add("integers",newint)
            newline = "nsub = 1"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            for nindex in range(len(subcommon)-1):
               newline = string.join(["IF (",subcommon[nindex].show(),"b /= ",subcommon[nindex+1].show(),\
                                      "b) nsub = nsub + 1"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
         if (factorialforsuper and factorialforsub):
            newdbl = "FACTORIAL"
            newcode.add("doubles",newdbl)
            newcode.add("externals",newdbl)
            factor = "FACTORIAL(nsuper)*FACTORIAL(nsub)"
         elif (factorialforsuper):
            newdbl = "FACTORIAL"
            newcode.add("doubles",newdbl)
            newcode.add("externals",newdbl)
            factor = "FACTORIAL(nsuper)"
         elif (factorialforsub):
            newdbl = "FACTORIAL"
            newcode.add("doubles",newdbl)
            newcode.add("externals",newdbl)
            factor = "FACTORIAL(nsub)"
         else:
            factor = "1.0d0"
   
         # perform contraction and store the result
         if (relativistic):
            newcode.add("doublecomplexallocatables","c_sort")
         else:
            newcode.add("doubleallocatables","c_sort")
         if (three):
            newline = "ALLOCATE(c_sort(dima_sort*dimb_sort))"
         else:
            newline = "ALLOCATE(c_sort(dima_sort))"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         if (relativistic):
            newline = string.join(["CALL ZGEMM('T','N',dima_sort,dimb_sort,dim_common,DCMPLX(",factor,\
                                   ",0.0d0),a_sort,dim_common,b_sort,dim_common,DCMPLX(0.0d0,0.0d0),c_sort,dima_sort)"],"")
         else:
            newline = string.join(["CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,",factor,\
                                   ",a_sort,dim_common,b_sort,dim_common,0.0d0,c_sort,dima_sort)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      # create an MA for unsorted tensor 0
      if (relativistic):
         newcode.add("doublecomplexallocatables","c")
      else:
         newcode.add("doubleallocatables","c")
      if (three):
         newline = "ALLOCATE(c(dima_sort*dimb_sort))"
      else:
         newline = "ALLOCATE(c(dima_sort))"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
         
      # mapping to a permutation symmetry unique block
      newcode.pointer = newcode.pointer - 1
      newcode.setamark(3)
      newcode.pointer = newcode.pointer + 1
      ifblock = 0
      for npermutation in range(len(self.factor.permutations)):
         permutation = self.factor.permutations[npermutation]
         indexesintheoriginalorder = copy.deepcopy(superglobalzero + superlocalzero + subglobalzero + sublocalzero)
         permutedindexes = performpermutation(indexesintheoriginalorder,permutation,1)
         newline = ""
         for nindex in range(len(permutedindexes)-1):
            # no IF across super and sub indexes
            if (nindex == len(permutedindexes)/2-1):
               continue
            indexa = permutedindexes[nindex]
            indexb = permutedindexes[nindex+1]
            if (indexa.isin(superglobalone) and indexb.isin(superglobalone)):
               continue
            if (indexa.isin(superglobaltwo) and indexb.isin(superglobaltwo)):
               continue
            if (indexa.isin(superlocalone) and indexb.isin(superlocalone)):
               continue
            if (indexa.isin(superlocaltwo) and indexb.isin(superlocaltwo)):
               continue
            if (indexa.isin(subglobalone) and indexb.isin(subglobalone)):
               continue
            if (indexa.isin(subglobaltwo) and indexb.isin(subglobaltwo)):
               continue
            if (indexa.isin(sublocalone) and indexb.isin(sublocalone)):
               continue
            if (indexa.isin(sublocaltwo) and indexb.isin(sublocaltwo)):
               continue
            if (indexa.isin(superglobalone) and indexb.isin(superlocalone)):
               continue
            if (indexa.isin(superglobalone) and indexb.isin(superlocaltwo)):
               continue
            if (indexa.isin(superglobaltwo) and indexb.isin(superlocalone)):
               continue
            if (indexa.isin(superglobaltwo) and indexb.isin(superlocaltwo)):
               continue
            if (indexa.isin(superlocalone) and indexb.isin(superglobalone)):
               continue
            if (indexa.isin(superlocalone) and indexb.isin(superglobaltwo)):
               continue
            if (indexa.isin(superlocaltwo) and indexb.isin(superglobalone)):
               continue
            if (indexa.isin(superlocaltwo) and indexb.isin(superglobaltwo)):
               continue
            if (indexa.isin(subglobalone) and indexb.isin(sublocalone)):
               continue
            if (indexa.isin(subglobalone) and indexb.isin(sublocaltwo)):
               continue
            if (indexa.isin(subglobaltwo) and indexb.isin(sublocalone)):
               continue
            if (indexa.isin(subglobaltwo) and indexb.isin(sublocaltwo)):
               continue
            if (indexa.isin(sublocalone) and indexb.isin(subglobalone)):
               continue
            if (indexa.isin(sublocalone) and indexb.isin(subglobaltwo)):
               continue
            if (indexa.isin(sublocaltwo) and indexb.isin(subglobalone)):
               continue
            if (indexa.isin(sublocaltwo) and indexb.isin(subglobaltwo)):
               continue
            if (indexa.type != indexb.type):
               continue
            inequality = " <= "
            if (newline):
               newline = string.join([newline," .and. (",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
            else:
               newline = string.join(["IF ((",indexa.show(),"b",inequality,indexb.show(),"b)"],"")
         if (newline):
            newline = string.join([newline,") THEN"],"")
            newcode.pointer = newcode.getamark(3) + 1
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            newline = "END IF"
            newcode.statements.insert(newcode.pointer,newline)
            newcode.setamark(3)
            ifblock = 1

         # sort indexes of tensor 0
         doloop = 0
         for nindex in range(len(permutedindexes)):
            index = permutedindexes[nindex]
            newint = index.show()
            newcode.add("integers",newint)
            newline = string.join(["DO ",newint," = 1,range(",newint,"b)"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1
            newline = "END DO"
            newcode.statements.insert(newcode.pointer,newline)
            if (nindex == 0):
               if (ifblock):
                  newcode.setamark(4)
                  doloop = 1
               else:
                  newcode.setamark(3)
                  doloop = 1
         newline = ""
         newlineend = ""
         if (three):
            newcode.add("integers","idimc_sort")
         sorted = superglobalone + superlocalone + subglobalone + sublocalone \
                + superglobaltwo + superlocaltwo + subglobaltwo + sublocaltwo
         for nindex in range(len(sorted)):
            if (newline == ""):
               if (three):
                  newline = string.join(["idimc_sort = ",sorted[nindex].show()],"")
               else:
                  newline = string.join(["idima_sort = ",sorted[nindex].show()],"")
            else:
               newline = string.join([newline," + range(",sorted[nindex-1].show(),"b)"\
                                      " * ((",sorted[nindex].show()," - 1)"],"")
               newlineend = string.join([newlineend,")"],"")
         newline = string.join([newline,newlineend],"")
         if (not newline):
            if (three):
               newline = "idimc_sort = 1"
            else:
               newline = "idima_sort = 1"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newline = ""
         newlineend = ""
         newcode.add("integers","idimc")
         for nindex in range(len(permutedindexes)-1,-1,-1):
            if (newline == ""):
               newline = string.join(["idimc = ",permutedindexes[nindex].show()],"")
            else:
               newline = string.join([newline," + range(",permutedindexes[nindex+1].show(),"b)"\
                                      " * ((",permutedindexes[nindex].show()," - 1)"],"")
               newlineend = string.join([newlineend,")"],"")
         if (not newline):
            newline = "idimc = 1"
         newline = string.join([newline,newlineend],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         if (three):
            if (self.factor.coefficients[npermutation] == 1.0):
               newline = "c(idimc) = c_sort(idimc_sort)"
            elif (self.factor.coefficients[npermutation] == - 1.0):
               newline = "c(idimc) = - c_sort(idimc_sort)"
            else:
               newline = string.join(["c(idimc) = ",repr(self.factor.coefficients[npermutation]),\
                                      "d0 * c_sort(idimc_sort)"],"")
         else:
            if (self.factor.coefficients[npermutation] == 1.0):
               newline = "c(idimc) = a_sort(idima_sort)"
            elif (self.factor.coefficients[npermutation] == - 1.0):
               newline = "c(idimc) = - a_sort(idima_sort)"
            else:
               newline = string.join(["c(idimc) = ",repr(self.factor.coefficients[npermutation]),\
                                      "d0 * a_sort(idima_sort)"],"")
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

         # accumulate a block
         if (ifblock):
            newcode.pointer = newcode.getamark(4) + 1
         elif (doloop):
            newcode.pointer = newcode.getamark(3) + 1
         arguments = ""
         argumentsend = ""
         for nindex in range(len(permutedindexes)-1,-1,-1):
            if (permutedindexes[nindex].type == "hole"):
               boffset = "b - 1"
            else:
               boffset = "b - noab - 1"
            if (arguments == ""):
               if (three):
                  arguments = string.join(["d_c,c,dima_sort*dimb_sort,c_offset(1 + ",permutedindexes[nindex].show(),boffset],"")
               else:
                  arguments = string.join(["d_c,c,dima_sort,c_offset(1 + ",permutedindexes[nindex].show(),boffset],"")
            else:
               if (permutedindexes[nindex+1].type == "hole"):
                  arguments = string.join([arguments," + noab * (",permutedindexes[nindex].show(),boffset],"")
               else:
                  arguments = string.join([arguments," + nvab * (",permutedindexes[nindex].show(),boffset],"")
               argumentsend = string.join([argumentsend,")"],"")
         if (arguments == ""):
            if (three):
               arguments = "d_c,c,dima_sort*dimb_sort,c_offset(1"
            else:
               arguments = "d_c,c,dima_sort,c_offset(1"
            argumentsend = ""
         arguments = string.join([arguments,argumentsend,")"],"")
         if (relativistic):
            newline = string.join(["CALL R4D_ADD_BLOCK(",arguments,")"],"")
         else:
            newline = string.join(["CALL ADD_BLOCK(",arguments,")"],"")
         newcode.statements.insert(newcode.pointer,newline)
         if (not ifblock):
            newcode.setamark(3)
         newcode.pointer = newcode.pointer + 1

      newcode.pointer = newcode.getamark(3) + 1
      newline = "DEALLOCATE(c)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      if (three):
         newline = "DEALLOCATE(c_sort)"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newline = "DEALLOCATE(b_sort)"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
         newline = "DEALLOCATE(a_sort)"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1
      else:
         newline = "DEALLOCATE(a_sort)"
         newcode.statements.insert(newcode.pointer,newline)
         newcode.pointer = newcode.pointer + 1

      # close the subroutine
      newcode.pointer = len(newcode.statements)
      newline = "next = NXTASK(-nprocs,1)"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "call GA_SYNC()"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      return newcode

class OperationTree:
 
   def __init__(self,contraction=NoOperation(),common=[],children=[],sisters=[]):
      """Creates a tree of contraction operations"""
      self.contraction = contraction
      self.common = common
      self.children = children
      self.sisters = sisters

   def isoperation(self):
      """Returns true if the operation has a valid contraction"""
      return self.contraction.isoperation()
 
   def __str__(self):
      """Prints the content"""
      print ""
      for line in self.show():
         print line
      return ""

   def show(self,ntab=0,verbose=1):
      """Returns the contents as a list of strings"""
      show = []
      if (self.sisters):
         for sister in self.sisters:
            show.append(string.join(["    "*ntab,sister.show(verbose)],""))
      if (self.contraction.isoperation()):
         show.append(string.join(["    "*(ntab-1),self.contraction.show(verbose)],""))
      for child in self.children:
         if (child.isoperation()):
            show = show + child.show(ntab+1,verbose)
      return show

   def tex(self,ntab=0,verbose=1):
      """Returns the contents as a list of strings in LaTeX format"""
      show = []
      show.append("\\begin{eqnarray}")
      if (self.contraction.isoperation()):
         show.append(string.join(["&&\\hspace{",repr(ntab-1),"cm}",self.contraction.tex(verbose),"\\nonumber\\\\"],""))
      if (self.sisters):
         for sister in self.sisters:
            show.append(string.join(["&&\\hspace{",repr(ntab),"cm}",sister.tex(verbose),"\\nonumber\\\\"],""))
      for child in self.children:
         if (child.isoperation()):
            show = show + child.texa(ntab+1,verbose)
      originallength = len(show[len(show)-1])
      show[len(show)-1] = show[len(show)-1][0:originallength-2]
      show.append("\\end{eqnarray}")
      return show

   def texa(self,ntab=0,verbose=1):
      """Returns the contents as a list of strings in LaTeX format"""
      show = []
      if (self.contraction.isoperation()):
         show.append(string.join(["&&\\hspace{",repr(ntab-1),"cm}",self.contraction.tex(verbose),"\\nonumber\\\\"],""))
      if (self.sisters):
         for sister in self.sisters:
            show.append(string.join(["&&\\hspace{",repr(ntab),"cm}",sister.tex(verbose),"\\nonumber\\\\"],""))
      for child in self.children:
         if (child.isoperation()):
            show = show + child.texa(ntab+1,verbose)
      return show

   def textable(self,name=""):
      """Returns the contents as a list of strings in LaTeX table format"""
      show = []
      show.append("\\begin{table}")
      show.append("\\begin{tabular}{l}")
      if (self.sisters):
         for sister in self.sisters:
            show.append(string.join(["$\\displaystyle ",sister.tensors[0].textable()," = ",sister.textable(),"$\\\\"],""))
      newline = ""
      counter = 0
      for child in self.children:
         if (child.isoperation()):
            counter = counter + 1
            newname = string.join([name,str(counter)],"")
            show = show + child.textablea(newname)
            if (child.contraction.isoperation()):
               if (not newline):
                  newline = string.join(["$\\displaystyle ",child.contraction.tensors[0].textable(name)," = "],"")
               newline = string.join([newline,child.contraction.textable(newname)],"")
      if (newline):
         newline = string.join([newline,"$\\\\"],"")
         show.append(newline)
      show.append("\\end{tabular}")
      show.append("\\end{table}")
      return show

   def textablea(self,name):
      """Returns the contents as a list of strings in LaTeX table format"""
      show = []
      if (self.sisters):
         for sister in self.sisters:
            show.append(string.join(["$\\displaystyle ",sister.tensors[0].textable()," = ",sister.textable(),"$\\\\"],""))
      newline = ""
      counter = 0
      for child in self.children:
         if (child.isoperation()):
            counter = counter + 1
            newname = string.join([name,str(counter)],"")
            show = show + child.textablea(newname)
            if (child.contraction.isoperation()):
               if (not newline):
                  newline = string.join(["$\\displaystyle ",child.contraction.tensors[0].textable(name)," = "],"")
               newline = string.join([newline,child.contraction.textable(newname)],"")
      if (newline):
         newline = string.join([newline,"$\\\\"],"")
         show.append(newline)
      return show
 
   def tensorslist(self,list=[]):
      """Returns the non-redundant list of the names of tensors appearing in the tree"""
      if (self.contraction.isoperation()):
         for tensor in self.contraction.tensors:
            if ((tensor.type == 'i') or (tensor.type == 'j')):
               tensorname = tensor.type + repr(tensor.label)
            else:
               tensorname = tensor.type + repr(len(tensor.indexes)/2)
            if (tensorname not in list):
               list.append(tensorname)
      if (self.children):
         for child in self.children:
            list = child.tensorslist(list)
      return list

   def usesindexlabel(self,label):
      """Recursively examines whether the index label is used"""
      if (self.contraction.usesindexlabel(label)):
         return 1
      if (self.common):
         if (self.common.usesindexlabel(label)):
            return 1
      for child in self.children:
         if (child.usesindexlabel(label)):
            return 1
      return 0 

   def relabelsone(self,oldlabel,newlabel):
      """Renames an index label in a whole operation tree"""
      self.contraction.relabels(oldlabel,newlabel)
      if (self.common):
         self.common.relabels(oldlabel,newlabel)
      for child in self.children:
         child.relabelsone(oldlabel,newlabel)

   def relabelstwo(self,another,selflabel,anotherlabel,reserved=[]):
      """Renames index labels that are among summation indexes and are arbitrary so that the two input operation trees look more alike"""
      newlabel = 0
      exist = 1
      while (exist):
         newlabel = newlabel + 1
         exist = 0
         if (self.usesindexlabel(newlabel)):
            exist = 1
         if (another.usesindexlabel(newlabel)):
            exist = 1
         if (newlabel in reserved):
            exist = 1
      self.relabelsone(selflabel,newlabel)
      another.relabelsone(anotherlabel,newlabel)

   def swapindexes(self,indexone,indextwo):
      """Swaps two indexes in a whole operation tree"""
      self.contraction.swapindexes(indexone,indextwo)
      if (self.common):
         self.common.swapindexes(indexone,indextwo)
      for child in self.children:
         child.swapindexes(indexone,indextwo)

   def isfactorizablewith(self,another,reserved=[],verbose=0):
      """Returns true if two elementary tensor contractions can be factorized"""

      selfcopy = OperationTree()
      anothercopy = OperationTree()
      selfcopy.contraction = copy.deepcopy(self.contraction)
      selfcopy.common = copy.deepcopy(self.common)
      selfcopy.children = copy.deepcopy(self.children)
      anothercopy.contraction = copy.deepcopy(another.contraction)
      anothercopy.common = copy.deepcopy(another.common)
      anothercopy.children = copy.deepcopy(another.children)

      # check if the two operator trees have valid contraction operations
      if (not (selfcopy.contraction.isoperation() and anothercopy.contraction.isoperation())):
         if (verbose):
            print "not a valid contraction operation"
         return 0

      # check if the output tensors have the identical form
      # (do not check the output tensor types since they are always "i")
      nselfindex = len(selfcopy.contraction.tensors[0].indexes)
      if (nselfindex != len(anothercopy.contraction.tensors[0].indexes)):
         if (verbose):
            print "output tensors imcompatible"
         return 0
      for nindex in range(nselfindex):
         if (not selfcopy.contraction.tensors[0].indexes[nindex].isidenticalto(anothercopy.contraction.tensors[0].indexes[nindex])):
            if (verbose):
               print "output tensors imcompatible"
            return 0

      # check if the factors (permutation operators) are compatible
      ratio = selfcopy.contraction.factor.isthesameas(anothercopy.contraction.factor)
      if (ratio == 0):
         if (verbose):
            print "factors incompatible"
         return 0

      # check if the summation indexes have the same number of holes and particles
      nselfhole = 0
      nselfparticle = 0
      nselfgeneral = 0
      if (not selfcopy.contraction.summation):
         if (verbose):
            print "summations incompatible"
         return 0
      for index in selfcopy.contraction.summation.indexes:
         if (index.type == 'hole'):
            nselfhole = nselfhole + 1
         elif (index.type == 'particle'):
            nselfparticle = nselfparticle + 1
         elif (index.type == 'general'):
            nselfgeneral = nselfgeneral + 1
      nanotherhole = 0
      nanotherparticle = 0
      nanothergeneral = 0
      if (not anothercopy.contraction.summation):
         if (verbose):
            print "summations incompatible"
         return 0
      for index in anothercopy.contraction.summation.indexes:
         if (index.type == 'hole'):
            nanotherhole = nanotherhole + 1
         elif (index.type == 'particle'):
            nanotherparticle = nanotherparticle + 1
         elif (index.type == 'general'):
            nanothergeneral = nanothergeneral + 1
      if ((nselfhole != nanotherhole) or (nselfparticle != nanotherparticle) or (nselfgeneral != nanothergeneral)):
         if (verbose):
            print "summations incompatible"
         return 0

      # find a common tensor
      # (note that the common tensor cannot be an intermediate, and must be identical to the common tensor
      #  used in prior factorizations)
      for nselftensor in range(1,3):
         selftensor = selfcopy.contraction.tensors[nselftensor]
         for nanothertensor in range(1,3):
            anothertensor = anothercopy.contraction.tensors[nanothertensor]
            found = 1
            if (selftensor.type != anothertensor.type):
               found = 0
            elif (selftensor.type == 'i'):
               found = 0
            elif (len(selftensor.indexes) != len(anothertensor.indexes)):
               found = 0
            else:
               for nindex in range(len(selftensor.indexes)):
                  selfindex = selftensor.indexes[nindex]
                  anotherindex = anothertensor.indexes[nindex]
                  if (selfindex.type != anotherindex.type):
                     found = 0
                  elif (selfindex.label != anotherindex.label):
                     if ((not selfcopy.contraction.summation.hastheindex(selfindex)) or \
                         (not anothercopy.contraction.summation.hastheindex(anotherindex))):
                        if (len(anothercopy.contraction.factor.permutations) > 1):
                           before = []
                           for permutation in anothercopy.contraction.factor.permutations:
                              if (len(permutation)/2 > len(before)):
                                 before = copy.deepcopy(permutation[0:len(permutation)/2])
                           beforeandafter = before
                           for nindex in range(len(before)):
                              if (before[nindex].isidenticalto(selfindex)):
                                 beforeandafter.append(anotherindex)
                              elif (before[nindex].isidenticalto(anotherindex)):
                                 beforeandafter.append(selfindex)
                              else:
                                 beforeandafter.append(before[nindex])
                           newfactor = Factor([1.0],[beforeandafter])
                           newfactor = anothercopy.contraction.factor.product(newfactor)
                           ratio = ratio * anothercopy.contraction.factor.isthesameas(newfactor)
                           if (ratio != 0):
                              anothercopy.swapindexes(selfindex,anotherindex)
                           else:
                              found = 0
                        else:
                           found = 0
                     else:
                        selfcopy.relabelstwo(anothercopy,selfindex.label,anotherindex.label,reserved)
               if (selfcopy.common):
                  if (not selfcopy.common.isidenticalto(selftensor)):
                     found = 0

            if (found):
               commonself = selftensor
               if (nselftensor == 1):
                  nuncommonself = 2
               else:
                  nuncommonself = 1
               uncommonself = selfcopy.contraction.tensors[nuncommonself]
               commonanother = anothertensor
               if (nanothertensor == 1):
                  nuncommonanother = 2
               else:
                  nuncommonanother = 1
               uncommonanother = anothercopy.contraction.tensors[nuncommonanother]

               newlabel = selfcopy.contraction.tensors[0].label + 1
               if (uncommonanother.type != "i"):
                  newintermediate = uncommonanother.duplicate()
                  newintermediate.type = "i"
                  newintermediate.label = newlabel
                  parity = newintermediate.sortindexes()
                  newcontraction = ElementaryTensorContraction(Factor([1.0],[[]]),[],[newintermediate,uncommonanother])
                  newchild = OperationTree(newcontraction,[],[])
                  anothercopy.contraction.tensors[nuncommonself] = newintermediate
                  anothercopy.children.insert(0,newchild)
               if (uncommonself.type != "i"):
                  newintermediate = uncommonself.duplicate()
                  newintermediate.type = "i"
                  newintermediate.label = newlabel
                  parity = newintermediate.sortindexes()
                  newcontraction = ElementaryTensorContraction(Factor([1.0],[[]]),[],[newintermediate,uncommonself])
                  newchild = OperationTree(newcontraction,[],[])
                  selfcopy.contraction.tensors[nuncommonself] = newintermediate
                  selfcopy.children.insert(0,newchild)
               selfcopy.common = commonself
               for child in anothercopy.children:
                  if (child.contraction.isoperation()):
                     child.contraction.factor.multiply(ratio)
                  selfcopy.children.append(child)

               # now we overwrite the self operation tree
               # (nothing is done to another after all, since it will be deleted anyway)
               self.contraction = copy.deepcopy(selfcopy.contraction)
               self.common = copy.deepcopy(selfcopy.common)
               self.children = copy.deepcopy(selfcopy.children)
               return 1

      # no common tensor found
      if (verbose):
         print "no common tensor"
      return 0

   def factorize(self,reserved=[],verbose=0):
      """Factors a common operation"""
      deletelist = []
      for nchilda in range(len(self.children)):
         if (nchilda in deletelist):
            continue
         childa = self.children[nchilda]
         if (childa.contraction.isoperation()):
            for nchildb in range(len(self.children)):
               if (nchildb in deletelist):
                  continue
               if (nchilda >= nchildb):
                  continue
               childb = self.children[nchildb]
               if (childb.contraction.isoperation()):
                  childc = copy.deepcopy(childa)
                  if (childa.isfactorizablewith(childb,reserved)):
                     if (verbose):
                        print childc
                        print " ... and"
                        print childb
                        print " ... have been consolidated into"
                        print childa
                        print ""
                     self.children[nchilda] = childa
                     deletelist.append(nchildb)
      if (deletelist):
         print " ... %d terms have been consolidated" %(len(deletelist))
         for nchildb in range(len(self.children)-1,-1,-1):
            if (nchildb in deletelist):
               del self.children[nchildb]
         return 1
      else:
         return 0

   def fullyfactorize(self,verbose=0,iteration=0,generation=1,reserved=[],globaltargetindexes=[]):
      """Performs factorize() recursively until fully factorized"""
      if ((iteration == 0) and (generation == 1)):
         if (self.sisters):
            raise RuntimeError, "factorize before reuseintermediate"
         reserved = []
         globaltargetindexes = []
         for targetindex in self.children[0].contraction.tensors[0].indexes:
            reserved.append(targetindex.label)
            globaltargetindexes.append(targetindex)
         print " ... commencing full factorization"
         print " ... initial contraction cost %d" %(self.operationcost())
      print " ... tensor contraction tier %d" %(generation)
      self.canonicalize(globaltargetindexes)
      while (self.factorize(reserved,verbose)):
         self.canonicalize(globaltargetindexes)
         iteration = iteration + 1
         print " ... iteration %d cost %d" %(iteration, self.operationcost())
      if (self.children):
         for child in self.children:
            child.fullyfactorize(verbose,0,generation+1,reserved,globaltargetindexes)
      if (generation == 1):
         print " ... exiting full factorization"
         print " ... final contraction cost %d" %(self.operationcost())
      self.sortindexes()
      return self

   def sortindexes(self):
      """Perform sortindex() for the whole operation tree"""
      if (self.contraction.isoperation()):
         self.contraction.sortindexes()
      if (self.children):
         for child in self.children:
            child.sortindexes()
      return self

   def operationcost(self,cost=0):
      """Returns a contraction cost"""
      if (self.contraction.isoperation()):
         if (len(self.contraction.tensors) > 2):
            cost = cost + 1
      cost = cost + len(self.sisters)
      for child in self.children:
         cost = child.operationcost(cost)
      return cost

   def canonicalize(self,globaltargetindexes):
      """Canonicalizes the whole operation tree expression"""

      self.canonicalizea(globaltargetindexes)

   def canonicalizea(self,globaltargetindexes):
      """Subsidiary to canonicalize()"""

      if (self.contraction.isoperation()):
         self.contraction.canonicalize(globaltargetindexes)

      for nchild in range(len(self.children)):
         self.children[nchild].canonicalizea(globaltargetindexes)

   def reuseintermediates(self):
      """Eliminates redundant intermediates from the tree"""

      print " ... commencing common subexpression elimination"

      if (self.children[0].contraction.isoperation()):
         globaltargetindexes = copy.deepcopy(self.children[0].contraction.tensors[0].indexes)
      else:
         return "The tree top must be an addition"

      intermediatesreused = self.collectuniqueintermediates()
      print " ... %d new intermediates are created" % (len(intermediatesreused))

      numberofreuse = self.renameintermediates(intermediatesreused,globaltargetindexes)
      print " ... %d reusable intermediates are found" % (numberofreuse)

      self.sisters = self.sisters + intermediatesreused

      print " ... final contraction cost %d" % (self.operationcost())

      # sort sisters
      for isister in range(len(self.sisters)):
         for jsister in range(len(self.sisters)):
            if ((isister > jsister) and (self.sisters[isister].tensors[0].label < self.sisters[jsister].tensors[0].label)):
               swap = copy.deepcopy(self.sisters[isister])
               self.sisters[isister] = copy.deepcopy(self.sisters[jsister])
               self.sisters[jsister] = copy.deepcopy(swap)

      return self

   def renameintermediates(self,intermediatesreused,globaltargetindexes):
      """Recursively relabels intermediates that are to be reused"""
      numberofreuse = 0
      for nchild in range(len(self.children)):
         child = self.children[nchild]
         if (child.contraction.isoperation()):
            if (len(child.contraction.tensors) == 3):
               if ((child.contraction.tensors[1].type != "i") and (child.contraction.tensors[2].type != "i")):
                  for nintermediate in range(len(intermediatesreused)):
                     intermediate = intermediatesreused[nintermediate]
                     ratio = intermediate.isthesameas(child.contraction,globaltargetindexes)
                     if (ratio != 0):
                        numberofreuse = numberofreuse + 1
                        if (intermediate.tensors[0].type != "j"):
                           found = 0
                           newlabel = 0
                           while (not found):
                              newlabel = newlabel + 1
                              alreadyused = 0
                              for anotherintermediate in intermediatesreused:
                                 if ((anotherintermediate.tensors[0].type == "j") and (anotherintermediate.tensors[0].label == newlabel)):
                                    alreadyused = 1
                              if (not alreadyused):
                                 found = 1
                           intermediatesreused[nintermediate].tensors[0].type = "j"
                           intermediatesreused[nintermediate].tensors[0].label = newlabel
                        del self.children[nchild].contraction.tensors[2]
                        self.children[nchild].contraction.summation = []
                        self.children[nchild].contraction.tensors[1] = copy.deepcopy(child.contraction.tensors[0])
                        self.children[nchild].contraction.tensors[1].type = "j"
                        self.children[nchild].contraction.tensors[1].label = intermediatesreused[nintermediate].tensors[0].label
                        self.children[nchild].contraction.factor = Factor([ratio],[[]])
                        break
         numberofreuse = numberofreuse + self.children[nchild].renameintermediates(intermediatesreused,globaltargetindexes)
      return numberofreuse

   def collectintermediates(self):
      """Collects intermediates from the tree edges"""
      listofcontractions = []
      for child in self.children:
         if (child.contraction.isoperation()):
            if (len(child.contraction.tensors) == 3):
               if ((child.contraction.tensors[1].type != "i") and (child.contraction.tensors[2].type != "i")):
                  listofcontractions.append(copy.deepcopy(child.contraction))
         listofcontractions = listofcontractions + child.collectintermediates()
      return listofcontractions 

   def collectuniqueintermediates(self):
      """Collects unique and reusable intermediates from the tree edges"""

      if (self.children[0].contraction.isoperation()):
         globaltargetindexes = copy.deepcopy(self.children[0].contraction.tensors[0].indexes)
      else:
         return "The tree top must be an addition"

      listofcontractions = self.collectintermediates()

      intermediatesreused = []
      for icontraction in range(len(listofcontractions)):
         reused = 0
         for jcontraction in range(len(listofcontractions)):
            if ((jcontraction > icontraction) and (not reused)):
               ratio = listofcontractions[icontraction].isthesameas(listofcontractions[jcontraction],globaltargetindexes)
               if (ratio != 0):
                  intermediatesreused.append(listofcontractions[icontraction])
                  reused = 1

      intermediatesunique = []
      for icontraction in range(len(intermediatesreused)):
         unique = 1
         for jcontraction in range(len(intermediatesreused)):
            if (jcontraction > icontraction):
               ratio = intermediatesreused[icontraction].isthesameas(intermediatesreused[jcontraction],globaltargetindexes)
               if (ratio != 0):
                  unique = 0
         if (unique):
            intermediatesunique.append(intermediatesreused[icontraction])

      return intermediatesunique
 
   def tensortypes(self,types=[[],[],[],[],[]]):
      """Returns a list of tensor types classified according to excitation, deexcitation, and general"""
      alltypes = []
      excitationtypes = []
      deexcitationtypes = []
      intermediatetypes = []
      generaltypes = []
      for child in self.children:
         result = child.tensortypes()
         for type in result[0]:
            if (type not in alltypes):
               alltypes.append(type)
         for type in result[1]:
            if (type not in excitationtypes):
               excitationtypes.append(type)
         for type in result[2]:
            if (type not in deexcitationtypes):
               deexcitationtypes.append(type)
         for type in result[3]:
            if (type not in intermediatetypes):
               intermediatetypes.append(type)
         for type in result[4]:
            if (type not in generaltypes):
               generaltypes.append(type)
      if (self.contraction.isoperation()):
         for tensor in self.contraction.tensors:
            if (tensor.type not in alltypes):
               alltypes.append(tensor.type)
            excitation = 1
            for nindex in range(len(tensor.indexes)/2):
               if (tensor.indexes[nindex].ishole()):
                  excitation = 0
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               if (tensor.indexes[nindex].isparticle()):
                  excitation = 0
            deexcitation = 1
            for nindex in range(len(tensor.indexes)/2):
               if (tensor.indexes[nindex].isparticle()):
                  deexcitation = 0
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               if (tensor.indexes[nindex].ishole()):
                  deexcitation = 0
            if ((tensor.type == "i") or (tensor.type == "j")):
               if (tensor.type not in intermediatetypes):
                  intermediatetypes.append(tensor.type)
            elif (excitation):
               if (tensor.type not in excitationtypes):
                  excitationtypes.append(tensor.type)
            elif (deexcitation):
               if (tensor.type not in deexcitationtypes):
                  deexcitationtypes.append(tensor.type)
            else:
               if (tensor.type not in generaltypes):
                  generaltypes.append(tensor.type)
      for sister in self.sisters:
         if (sister.isoperation()):
            for tensor in sister.tensors:
               if (tensor.type not in alltypes):
                  alltypes.append(tensor.type)
               excitation = 1
               for nindex in range(len(tensor.indexes)/2):
                  if (tensor.indexes[nindex].ishole()):
                     excitation = 0
               for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
                  if (tensor.indexes[nindex].isparticle()):
                     excitation = 0
               deexcitation = 1
               for nindex in range(len(tensor.indexes)/2):
                  if (tensor.indexes[nindex].isparticle()):
                     deexcitation = 0
               for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
                  if (tensor.indexes[nindex].ishole()):
                     deexcitation = 0
               if ((tensor.type == "i") or (tensor.type == "j")):
                  if (tensor.type not in intermediatetypes):
                     intermediatetypes.append(tensor.type)
               elif (excitation):
                  if (tensor.type not in excitationtypes):
                     excitationtypes.append(tensor.type)
               elif (deexcitation):
                  if (tensor.type not in deexcitationtypes):
                     deexcitationtypes.append(tensor.type)
               else:
                  if (tensor.type not in generaltypes):
                     generaltypes.append(tensor.type)
      newexcitationtypes = []
      for type in excitationtypes:
         if (type not in generaltypes):
            newexcitationtypes.append(type)
      newdeexcitationtypes = []
      for type in deexcitationtypes:
         if (type not in generaltypes):
            newdeexcitationtypes.append(type)
    
      return [alltypes,newexcitationtypes,newdeexcitationtypes,intermediatetypes,generaltypes]

   def fortran77(self,filename="NONAME",excitation=[],deexcitation=[],intermediate=[],general=[]):
      """Suggests an implementation in Fortran77 for the whole operation tree"""

      print " ... generating a Fortran77 code"
      print " "
      all = excitation+deexcitation+intermediate+general
      if (not all):
         types = self.tensortypes()
      else:
         types = [all,excitation,deexcitation,intermediate,general]
      for type in types[0]:
         if (type in types[1]):
            print " '"+type+"' is an excitaion tensor"
            if ((type == 'v') and (not all)):
               print " !!! WARNING !!! unusual naming convention"
            if ((type == 'f') and (not all)):
               print " !!! WARNING !!! unusual naming convention"
            if ((type == 'd') and (not all)):
               print " !!! WARNING !!! unusual naming convention"
         elif (type in types[2]):
            print " '"+type+"' is a de-excitaion tensor"
            if ((type == 'v') and (not all)):
               print " !!! WARNING !!! unusual naming convention"
            if ((type == 'f') and (not all)):
               print " !!! WARNING !!! unusual naming convention"
            if ((type == 'd') and (not all)):
               print " !!! WARNING !!! unusual naming convention"
         elif (type in types[3]):
            print " '"+type+"' is an intermediate tensor"
         elif (type in types[4]):
            print " '"+type+"' is a general tensor"
         else:
            raise RuntimeError, "unknown tensor type"
      print " "

      newlistofcodes = ListofCodes()

      # callees (tensor contraction subroutines called from the main)
      callees = ListofCodes()

      # copy of self will be reduced as we write the program
      selfcopy = OperationTree()
      selfcopy.contraction = copy.deepcopy(self.contraction)
      selfcopy.common = copy.deepcopy(self.common)
      selfcopy.children = copy.deepcopy(self.children)
      selfcopy.sisters = copy.deepcopy(self.sisters)

      # target indexes
      if (selfcopy.children[0].contraction.isoperation()):
         globaltargetindexes = copy.deepcopy(selfcopy.children[0].contraction.tensors[0].indexes)
         print "Target indices"
         printindexes(globaltargetindexes)
         print " "
      else:
         return "The tree top must be an addition"

      newcode = Code("Fortran77",filename)

      # header
      for newline in self.show(0,0):
         newline = string.join(["!",newline],"")
         newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)
      newline = '#include "util.fh"'
      newcode.add("headers",newline)
      newline = '#include "errquit.fh"'
      newcode.add("headers",newline)
      newline = '#include "tce.fh"'
      newcode.add("headers",newline)
      
      # loop over the tree
      newcode.join(selfcopy.fortran77a(filename,globaltargetindexes,types,callees).expand())

      # antisymmetrizer
      newcode.pointer = len(newcode.statements)
#     if (globaltargetindexes):
#        newline = string.join(["CALL ANTISYM_",filename,"(d_i0,k_i0_offset)"],"")
#        newcode.statements.insert(newcode.pointer,newline)
#        newcode.pointer = newcode.pointer + 1
#        antisymmetrizer = selfcopy.children[0].contraction.tensors[0].fortran77x(filename)

      # close the subroutine
#     newline = "CALL RECONCILEFILE(d_i0,size_i0)"
#     newcode.statements.insert(newcode.pointer,newline)
#     newcode.pointer = newcode.pointer + 1
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # append the callees
      newlistofcodes.add(newcode)
      newlistofcodes.join(callees)
#     if (globaltargetindexes):
#        newlistofcodes.add(antisymmetrizer)

      newlistofcodes.list[0].sortarguments()
      newlistofcodes.list[0].removeredundantio()
      return newlistofcodes

   def fortran77a(self,subroutinename,globaltargetindexes,types,callees):
      """Returns a part of program that is generated by recursively interpreting the tree"""

      newcode = Code("Fortran77",subroutinename)

      # check if we need to proceed
      if ((not self.children) and (not self.sisters)):
         return newcode
      else:
         empty = 1
         for child in self.children:
            if (child.contraction.isoperation()):
               empty = 0
         if (empty):
            return newcode

      # loop over sisters
      if (self.sisters):
         counter = 0
         for sister in self.sisters:
            counter = counter + 1
            name = string.join([subroutinename,"__",repr(counter)],"")
 
            # Tensor 1
            superglobalzero = []
            subglobalzero = []
            superlocalzero = []
            sublocalzero = []
            for nindex in range(len(sister.tensors[0].indexes)/2):
               index = sister.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  superglobalzero.append(index)
               else:
                  superlocalzero.append(index)
            for nindex in range(len(sister.tensors[0].indexes)/2, \
                                len(sister.tensors[0].indexes)):
               index = sister.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  subglobalzero.append(index)
               else:
                  sublocalzero.append(index)

            d_c = string.join(["d_j",repr(sister.tensors[0].label)],"")
            newcode.add("integers",d_c)
            l_c_offset = string.join(["l_j",repr(sister.tensors[0].label),"_offset"],"")
            k_c_offset = string.join(["k_j",repr(sister.tensors[0].label),"_offset"],"")
            size_c = string.join(["size_j",repr(sister.tensors[0].label)],"")
            newcode.add("integers",k_c_offset)
            newcode.add("integers",l_c_offset)
            newcode.add("integers",size_c)
 
            # generate the contraction callee
            callee = sister.fortran77(globaltargetindexes,types,name)
            callees.add(callee)

            # Tensor 2
            d_a = string.join(["d_",sister.tensors[1].type,repr(len(sister.tensors[1].indexes)/2)],"") 
            newcode.add("integers",d_a)
            newcode.add("arguments",d_a)
            k_a_offset = string.join(["k_",sister.tensors[1].type,\
                         repr(len(sister.tensors[1].indexes)/2),"_offset"],"")
            newcode.add("integers",k_a_offset)
            newcode.add("arguments",k_a_offset)

            # Tensor 3
            d_b = string.join(["d_",sister.tensors[2].type,repr(len(sister.tensors[2].indexes)/2)],"") 
            newcode.add("integers",d_b)
            newcode.add("arguments",d_b)
            k_b_offset = string.join(["k_",sister.tensors[2].type,\
                         repr(len(sister.tensors[2].indexes)/2),"_offset"],"")
            newcode.add("integers",k_b_offset)
            newcode.add("arguments",k_b_offset)

            # dump the caller
            newline = string.join(["CALL OFFSET_",name,"(",d_c,",",l_c_offset,",",k_c_offset,",",size_c,")"],"")
            newcode.statements.insert(0,newline)
            newchar = "filename"
            newcode.add("characters",newchar)
            filename = string.join([name,"_j",repr(sister.tensors[0].label)],"")
            newline = string.join(["CALL TCE_FILENAME('",filename,"',filename)"],"")
            newcode.statements.insert(0,newline)
            newline = string.join(["CALL CREATEFILE(filename,",d_c,",",size_c,")"],"")
            newcode.statements.insert(0,newline)
#           newline = string.join(["CALL DRATOGA(",d_b,")"],"")
#           newcode.statements.insert(0,newline)
#           newline = string.join(["CALL DRATOGA(",d_a,")"],"")
#           newcode.statements.insert(0,newline)
            callee = sister.tensors[0].fortran77y(globaltargetindexes,types,name)
            callees.add(callee)
            argument = string.join([d_a,",",k_a_offset],"")
            argument = string.join([argument,",",d_b,",",k_b_offset],"")
            argument = string.join([argument,",",d_c,",",k_c_offset],"")
            newline = string.join(["CALL ",name,"(",argument,")"],"")
            newcode.statements.insert(0,newline)
#           newline = string.join(["CALL GATODRA(",d_a,")"],"")
#           newcode.statements.insert(0,newline)
#           newline = string.join(["CALL GATODRA(",d_b,")"],"")
#           newcode.statements.insert(0,newline)
            newline = string.join(["CALL RECONCILEFILE(",d_c,",",size_c,")"],"")
            newcode.statements.insert(0,newline)
 
      # get a filename for intermediate storage
      d_c = string.join(["d_i",repr(self.children[0].contraction.tensors[0].label)],"")
      newcode.add("integers",d_c)
      if (not self.contraction.isoperation()):
         newcode.add("arguments",d_c)
      l_c_offset = string.join(["l_i",repr(self.children[0].contraction.tensors[0].label),"_offset"],"")
      k_c_offset = string.join(["k_i",repr(self.children[0].contraction.tensors[0].label),"_offset"],"")
      size_c = string.join(["size_i",repr(self.children[0].contraction.tensors[0].label)],"")
      newcode.add("integers",k_c_offset)
#     newcode.add("integers",size_c)
      if (self.contraction.isoperation()):
         newcode.add("integers",l_c_offset)
      else:
         newcode.add("arguments",k_c_offset)
#        newcode.add("arguments",size_c)

      # loop over children
      if (self.contraction.isoperation()):
         createfile = 1
      else:
         createfile = 0
      counter = 0
      for nchild in range(len(self.children)):
         child = self.children[nchild]
         if (child.contraction.isoperation()):
            counter = counter + 1
            name = string.join([subroutinename,"_",repr(counter)],"")
 
            # Tensor 1
            superglobalzero = []
            subglobalzero = []
            superlocalzero = []
            sublocalzero = []
            for nindex in range(len(child.contraction.tensors[0].indexes)/2):
               index = child.contraction.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  superglobalzero.append(index)
               else:
                  superlocalzero.append(index)
            for nindex in range(len(child.contraction.tensors[0].indexes)/2, \
                                len(child.contraction.tensors[0].indexes)):
               index = child.contraction.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  subglobalzero.append(index)
               else:
                  sublocalzero.append(index)

            # generate the contraction callee
            callee = child.contraction.fortran77(globaltargetindexes,types,name)
            callees.add(callee)

            if (child.contraction.tensors[1].type == "i"):
               d_a = string.join(["d_i",repr(child.contraction.tensors[1].label)],"") 
               size_a = string.join(["size_i",repr(child.contraction.tensors[1].label)],"") 
               newcode.add("integers",d_a)
               k_a_offset = string.join(["k_i",repr(child.contraction.tensors[1].label),"_offset"],"")
               l_a_offset = string.join(["l_i",repr(child.contraction.tensors[1].label),"_offset"],"")
               newcode.add("integers",k_a_offset)
            elif (child.contraction.tensors[1].type == "j"):
               d_a = string.join(["d_j",repr(child.contraction.tensors[1].label)],"") 
               size_a = string.join(["size_j",repr(child.contraction.tensors[1].label)],"") 
               newcode.add("integers",d_a)
               k_a_offset = string.join(["k_j",repr(child.contraction.tensors[1].label),"_offset"],"")
               l_a_offset = string.join(["l_j",repr(child.contraction.tensors[1].label),"_offset"],"")
               newcode.add("integers",k_a_offset)
            else:
               d_a = string.join(["d_",child.contraction.tensors[1].type,repr(len(child.contraction.tensors[1].indexes)/2)],"") 
               newcode.add("integers",d_a)
               newcode.add("arguments",d_a)
               k_a_offset = string.join(["k_",child.contraction.tensors[1].type,\
                            repr(len(child.contraction.tensors[1].indexes)/2),"_offset"],"")
               newcode.add("integers",k_a_offset)
               newcode.add("arguments",k_a_offset)
            if (len(child.contraction.tensors) == 3):
               if (child.contraction.tensors[2].type == "i"):
                  d_b = string.join(["d_i",repr(child.contraction.tensors[2].label)],"") 
                  size_b = string.join(["size_i",repr(child.contraction.tensors[2].label)],"") 
                  newcode.add("integers",d_b)
                  k_b_offset = string.join(["k_i",repr(child.contraction.tensors[2].label),"_offset"],"")
                  l_b_offset = string.join(["l_i",repr(child.contraction.tensors[2].label),"_offset"],"")
                  newcode.add("integers",k_b_offset)
               elif (child.contraction.tensors[2].type == "j"):
                  d_b = string.join(["d_j",repr(child.contraction.tensors[2].label)],"") 
                  size_b = string.join(["size_j",repr(child.contraction.tensors[2].label)],"") 
                  newcode.add("integers",d_b)
                  k_b_offset = string.join(["k_j",repr(child.contraction.tensors[2].label),"_offset"],"")
                  l_b_offset = string.join(["l_j",repr(child.contraction.tensors[2].label),"_offset"],"")
                  newcode.add("integers",k_b_offset)
               else:
                  d_b = string.join(["d_",child.contraction.tensors[2].type,repr(len(child.contraction.tensors[2].indexes)/2)],"") 
                  newcode.add("integers",d_b)
                  newcode.add("arguments",d_b)
                  k_b_offset = string.join(["k_",child.contraction.tensors[2].type,\
                               repr(len(child.contraction.tensors[2].indexes)/2),"_offset"],"")
                  newcode.add("integers",k_b_offset)
                  newcode.add("arguments",k_b_offset)
            else:
               d_b = ""
            
            # dump the code
            if (createfile):
               newint = string.join(["size_i",repr(self.children[0].contraction.tensors[0].label)],"")
               newcode.add("integers",newint)
               newline = string.join(["CALL OFFSET_",name,"(",d_c,",",l_c_offset,",",k_c_offset,",",size_c,")"],"")
               newcode.statements.insert(0,newline)
               newchar = "filename"
               newcode.add("characters",newchar)
               filename = string.join([name,"_i",repr(child.contraction.tensors[0].label)],"")
               newline = string.join(["CALL TCE_FILENAME('",filename,"',filename)"],"")
               newcode.statements.insert(0,newline)
               newline = string.join(["CALL CREATEFILE(filename,",d_c,",",size_c,")"],"")
               newcode.statements.insert(0,newline)
               callee = child.contraction.tensors[0].fortran77y(globaltargetindexes,types,name)
               callees.add(callee)
               createfile = 0
            newcode.statements.insert(0,child.fortran77a(name,globaltargetindexes,types,callees))
#           newline = string.join(["CALL DRATOGA(",d_c,")"],"")
#           newcode.statements.insert(0,newline)
#           if (d_b):
#              newline = string.join(["CALL DRATOGA(",d_b,")"],"")
#              newcode.statements.insert(0,newline)
#           newline = string.join(["CALL DRATOGA(",d_a,")"],"")
#           newcode.statements.insert(0,newline)
            if (child.contraction.tensors[1].type == "i"):
               newline = string.join(["CALL RECONCILEFILE(",d_a,",",size_a,")"],"")
               newcode.statements.insert(0,newline)
            if (d_b):
               if (child.contraction.tensors[2].type == "i"):
                  newline = string.join(["CALL RECONCILEFILE(",d_b,",",size_b,")"],"")
                  newcode.statements.insert(0,newline)
            argument = string.join([d_a,",",k_a_offset],"")
            if (d_b):
               argument = string.join([argument,",",d_b,",",k_b_offset],"")
            argument = string.join([argument,",",d_c,",",k_c_offset],"")
            newline = string.join(["CALL ",name,"(",argument,")"],"")
            newcode.statements.insert(0,newline)
#           newline = string.join(["CALL GATODRA(",d_a,")"],"")
#           newcode.statements.insert(0,newline)
#           if (d_b):
#              newline = string.join(["CALL GATODRA(",d_b,")"],"")
#              newcode.statements.insert(0,newline)
#           newline = string.join(["CALL GATODRA(",d_c,")"],"")
#           newcode.statements.insert(0,newline)
            if (child.contraction.tensors[1].type == "i"):
               newline = string.join(["CALL DELETEFILE(",d_a,")"],"")
               newcode.statements.insert(0,newline)
               newline = string.join(["IF (.not.MA_POP_STACK(",l_a_offset,")) CALL ERRQUIT('",subroutinename,"',-1,MA_ERR)"],"")
               newcode.statements.insert(0,newline)
            if (d_b):
               if (child.contraction.tensors[2].type == "i"):
                  newline = string.join(["CALL DELETEFILE(",d_b,")"],"")
                  newcode.statements.insert(0,newline)
                  newline = string.join(["IF (.not.MA_POP_STACK(",l_b_offset,")) CALL ERRQUIT('",subroutinename,"',-1,MA_ERR)"],"")
                  newcode.statements.insert(0,newline)

      if (self.sisters):
         counter = 0
         for isister in range(len(self.sisters)-1,-1,-1):
            sister = self.sisters[isister]
            counter = counter + 1
            d_c = string.join(["d_j",repr(sister.tensors[0].label)],"")
            l_c_offset = string.join(["l_j",repr(sister.tensors[0].label),"_offset"],"")
            k_c_offset = string.join(["k_j",repr(sister.tensors[0].label),"_offset"],"")
            name = string.join([subroutinename,"_r_",repr(counter)],"")
            newline = string.join(["CALL DELETEFILE(",d_c,")"],"")
            newcode.statements.insert(0,newline)
            newline = string.join(["IF (.not.MA_POP_STACK(",l_c_offset,")) CALL ERRQUIT('",subroutinename,"',-1,MA_ERR)"],"")
            newcode.statements.insert(0,newline)
 
      newcode.reverse()
      return newcode

   def utchem(self,filename="NONAME",relativistic=0):
      """Suggests an implementation in Fortran90 for the whole operation tree"""

      newlistofcodes = ListofCodes()

      # callees (tensor contraction subroutines called from the main)
      callees = ListofCodes()

      # copy of self will be reduced as we write the program
      selfcopy = OperationTree()
      selfcopy.contraction = copy.deepcopy(self.contraction)
      selfcopy.common = copy.deepcopy(self.common)
      selfcopy.children = copy.deepcopy(self.children)
      selfcopy.sisters = copy.deepcopy(self.sisters)

      # target indexes
      if (selfcopy.children[0].contraction.isoperation()):
         globaltargetindexes = copy.deepcopy(selfcopy.children[0].contraction.tensors[0].indexes)
      else:
         return "The tree top must be an addition"

      if (relativistic):
         newcode = Code("Fortran90","R4D_"+filename)
      else:
         newcode = Code("Fortran90",filename)

      # header
      for newline in self.show(0,0):
         newline = string.join(["!",newline],"")
         newcode.add("headers",newline)
      newline = "USE UT_SYS_MODULE"
      newcode.add("headers",newline)
      newline = "USE UT_MOLINP_MODULE"
      newcode.add("headers",newline)
      if (relativistic):
         newline = "USE UT_R4DTCE_MODULE"
      else:
         newline = "USE UT_TCE_MODULE"
      newcode.add("headers",newline)
      newline = "IMPLICIT NONE"
      newcode.add("headers",newline)
      newline = '#include "global.fh"'
      newcode.add("headers",newline)
      newline = '#include "mafdecls.fh"'
      newcode.add("headers",newline)
      
      # loop over the tree
      newcode.join(selfcopy.utchema(filename,globaltargetindexes,callees,relativistic).expand())

      # antisymmetrizer
      newcode.pointer = len(newcode.statements)
#     if (globaltargetindexes):
#        newline = string.join(["CALL ANTISYM_",filename,"(d_i0,i0_offset)"],"")
#        newcode.statements.insert(newcode.pointer,newline)
#        newcode.pointer = newcode.pointer + 1
#        antisymmetrizer = selfcopy.children[0].contraction.tensors[0].utchemx(filename)

      # close the subroutine
#     newline = "CALL RECONCILEFILE(d_i0,size_i0)"
#     newcode.statements.insert(newcode.pointer,newline)
#     newcode.pointer = newcode.pointer + 1
      newline = "RETURN"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1
      newline = "END"
      newcode.statements.insert(newcode.pointer,newline)
      newcode.pointer = newcode.pointer + 1

      # append the callees
      newlistofcodes.add(newcode)
      newlistofcodes.join(callees)
#     if (globaltargetindexes):
#        newlistofcodes.add(antisymmetrizer)

      newlistofcodes.list[0].sortarguments()
      return newlistofcodes

   def utchema(self,subroutinename,globaltargetindexes,callees,relativistic=0):
      """Returns a part of program that is generated by recursively interpreting the tree"""

      newcode = Code("Fortran90",subroutinename)

      # check if we need to proceed
      if ((not self.children) and (not self.sisters)):
         return newcode
      else:
         empty = 1
         for child in self.children:
            if (child.contraction.isoperation()):
               empty = 0
         if (empty):
            return newcode

      # loop over sisters
      if (self.sisters):
         counter = 0
         for sister in self.sisters:
            counter = counter + 1
            name = string.join([subroutinename,"__",repr(counter)],"")
 
            # Tensor 1
            superglobalzero = []
            subglobalzero = []
            superlocalzero = []
            sublocalzero = []
            for nindex in range(len(sister.tensors[0].indexes)/2):
               index = sister.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  superglobalzero.append(index)
               else:
                  superlocalzero.append(index)
            for nindex in range(len(sister.tensors[0].indexes)/2, \
                                len(sister.tensors[0].indexes)):
               index = sister.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  subglobalzero.append(index)
               else:
                  sublocalzero.append(index)

            d_c = string.join(["d_j",repr(sister.tensors[0].label)],"")
            newcode.add("integers",d_c)
            c_offset = string.join(["j",repr(sister.tensors[0].label),"_offset"],"")
            size_c = string.join(["size_j",repr(sister.tensors[0].label)],"")
            newcode.add("integerarrays",c_offset)
            newcode.add("integers",size_c)
 
            # generate the contraction callee
            callee = sister.utchem(globaltargetindexes,name,relativistic)
            callees.add(callee)

            # Tensor 2
            d_a = string.join(["d_",sister.tensors[1].type,repr(len(sister.tensors[1].indexes)/2)],"") 
            newcode.add("integers",d_a)
            newcode.add("arguments",d_a)
            a_offset = string.join([sister.tensors[1].type,\
                         repr(len(sister.tensors[1].indexes)/2),"_offset"],"")
            newcode.add("integerarrays",a_offset)
            newcode.add("arguments",a_offset)

            # Tensor 3
            d_b = string.join(["d_",sister.tensors[2].type,repr(len(sister.tensors[2].indexes)/2)],"") 
            newcode.add("integers",d_b)
            newcode.add("arguments",d_b)
            b_offset = string.join([sister.tensors[2].type,\
                         repr(len(sister.tensors[2].indexes)/2),"_offset"],"")
            newcode.add("integerarrays",b_offset)
            newcode.add("arguments",b_offset)
 
            # allocate offsets
            arguments = ""
            for index in sister.tensors[0].indexes:
               if (index.type == "hole"):
                  factor = "noab"
               else:
                  factor = "nvab"
               if (arguments):
                  arguments = string.join([arguments,"*",factor],"")
               else:
                  arguments = factor
            if (arguments == ""):
               arguments = "1"
            if (relativistic):
               newline = string.join(["R4D_ALLOCATE(",c_offset,"(",arguments,"))"],"")
            else:
               newline = string.join(["ALLOCATE(",c_offset,"(",arguments,"))"],"")
            newcode.statements.insert(newcode.pointer,newline)
            newcode.pointer = newcode.pointer + 1

            # dump the caller
            if (relativistic):
               newline = string.join(["CALL R4D_OFFSET_",name,"(",d_c,",",c_offset,",",size_c,")"],"")
            else:
               newline = string.join(["CALL OFFSET_",name,"(",d_c,",",c_offset,",",size_c,")"],"")
            newcode.statements.insert(0,newline)
            newchar = "filename"
            newcode.add("characters",newchar)
            filename = string.join([name,"_j",repr(sister.tensors[0].label)],"")
            if (relativistic):
               newline = string.join(["CALL UT_R4DTCE_FILENAME('",filename,"',filename)"],"")
            else:
               newline = string.join(["CALL UT_TCE_FILENAME('",filename,"',filename)"],"")
            newcode.statements.insert(0,newline)
            if (relativistic):
               newline = string.join(["CALL R4D_CREATEFILE(filename,",d_c,",",size_c,")"],"")
            else:
               newline = string.join(["CALL CREATEFILE(filename,",d_c,",",size_c,")"],"")
            newcode.statements.insert(0,newline)
            callee = sister.tensors[0].utchemy(globaltargetindexes,name,relativistic)
            callees.add(callee)
            argument = string.join([d_a,",",a_offset],"")
            argument = string.join([argument,",",d_b,",",b_offset],"")
            argument = string.join([argument,",",d_c,",",c_offset],"")
            if (relativistic):
               newline = string.join(["CALL R4D_",name,"(",argument,")"],"")
            else:
               newline = string.join(["CALL ",name,"(",argument,")"],"")
            newcode.statements.insert(0,newline)
            if (relativistic):
               newline = string.join(["CALL R4D_RECONCILEFILE(",d_c,",",size_c,")"],"")
            else:
               newline = string.join(["CALL RECONCILEFILE(",d_c,",",size_c,")"],"")
            newcode.statements.insert(0,newline)
 
      # get a filename for intermediate storage
      d_c = string.join(["d_i",repr(self.children[0].contraction.tensors[0].label)],"")
      newcode.add("integers",d_c)
      if (not self.contraction.isoperation()):
         newcode.add("arguments",d_c)
      c_offset = string.join(["i",repr(self.children[0].contraction.tensors[0].label),"_offset"],"")
      size_c = string.join(["size_i",repr(self.children[0].contraction.tensors[0].label)],"")
#     newcode.add("integers",size_c)
      if (self.contraction.isoperation()):
         newcode.add("integerallocatables",c_offset)
      else:
         newcode.add("integerarrays",c_offset)
         newcode.add("arguments",c_offset)
#        newcode.add("arguments",size_c)

      # loop over children
      if (self.contraction.isoperation()):
         createfile = 1
      else:
         createfile = 0
      counter = 0
      for nchild in range(len(self.children)):
         child = self.children[nchild]
         if (child.contraction.isoperation()):
            counter = counter + 1
            name = string.join([subroutinename,"_",repr(counter)],"")
 
            # Tensor 1
            superglobalzero = []
            subglobalzero = []
            superlocalzero = []
            sublocalzero = []
            for nindex in range(len(child.contraction.tensors[0].indexes)/2):
               index = child.contraction.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  superglobalzero.append(index)
               else:
                  superlocalzero.append(index)
            for nindex in range(len(child.contraction.tensors[0].indexes)/2, \
                                len(child.contraction.tensors[0].indexes)):
               index = child.contraction.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  subglobalzero.append(index)
               else:
                  sublocalzero.append(index)

            # generate the contraction callee
            callee = child.contraction.utchem(globaltargetindexes,name,relativistic)
            callees.add(callee)

            if (child.contraction.tensors[1].type == "i"):
               d_a = string.join(["d_i",repr(child.contraction.tensors[1].label)],"") 
               size_a = string.join(["size_i",repr(child.contraction.tensors[1].label)],"") 
               newcode.add("integers",d_a)
               a_offset = string.join(["i",repr(child.contraction.tensors[1].label),"_offset"],"")
               newcode.add("integerarrays",a_offset)
            elif (child.contraction.tensors[1].type == "j"):
               d_a = string.join(["d_j",repr(child.contraction.tensors[1].label)],"") 
               size_a = string.join(["size_j",repr(child.contraction.tensors[1].label)],"") 
               newcode.add("integers",d_a)
               a_offset = string.join(["j",repr(child.contraction.tensors[1].label),"_offset"],"")
               newcode.add("integerarrays",a_offset)
            else:
               d_a = string.join(["d_",child.contraction.tensors[1].type,repr(len(child.contraction.tensors[1].indexes)/2)],"") 
               newcode.add("integers",d_a)
               newcode.add("arguments",d_a)
               a_offset = string.join([child.contraction.tensors[1].type,\
                            repr(len(child.contraction.tensors[1].indexes)/2),"_offset"],"")
               newcode.add("integerarrays",a_offset)
               newcode.add("arguments",a_offset)
            if (len(child.contraction.tensors) == 3):
               if (child.contraction.tensors[2].type == "i"):
                  d_b = string.join(["d_i",repr(child.contraction.tensors[2].label)],"") 
                  size_b = string.join(["size_i",repr(child.contraction.tensors[2].label)],"") 
                  newcode.add("integers",d_b)
                  b_offset = string.join(["i",repr(child.contraction.tensors[2].label),"_offset"],"")
                  newcode.add("integerallocatables",b_offset)
               elif (child.contraction.tensors[2].type == "j"):
                  d_b = string.join(["d_j",repr(child.contraction.tensors[2].label)],"") 
                  size_b = string.join(["size_j",repr(child.contraction.tensors[2].label)],"") 
                  newcode.add("integers",d_b)
                  b_offset = string.join(["j",repr(child.contraction.tensors[2].label),"_offset"],"")
                  newcode.add("integerallocatables",b_offset)
               else:
                  d_b = string.join(["d_",child.contraction.tensors[2].type,repr(len(child.contraction.tensors[2].indexes)/2)],"") 
                  newcode.add("integers",d_b)
                  newcode.add("arguments",d_b)
                  b_offset = string.join([child.contraction.tensors[2].type,\
                               repr(len(child.contraction.tensors[2].indexes)/2),"_offset"],"")
                  newcode.add("integerarrays",b_offset)
                  newcode.add("arguments",b_offset)
            else:
               d_b = ""
            
            # dump the code
            if (createfile):
 
               # allocate offsets
               arguments = ""
               for index in child.contraction.tensors[0].indexes:
                  if (index.type == "hole"):
                     factor = "noab"
                  else:
                     factor = "nvab"
                  if (arguments):
                     arguments = string.join([arguments,"*",factor],"")
                  else:
                     arguments = factor
               if (arguments == ""):
                  arguments = "1"
               newline = string.join(["ALLOCATE(",c_offset,"(",arguments,"))"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1

               newint = string.join(["size_i",repr(self.children[0].contraction.tensors[0].label)],"")
               newcode.add("integers",newint)
               if (relativistic):
                  newline = string.join(["CALL R4D_OFFSET_",name,"(",d_c,",",c_offset,",",size_c,")"],"")
               else:
                  newline = string.join(["CALL OFFSET_",name,"(",d_c,",",c_offset,",",size_c,")"],"")
               newcode.statements.insert(0,newline)
               newchar = "filename"
               newcode.add("characters",newchar)
               filename = string.join([name,"_i",repr(child.contraction.tensors[0].label)],"")
               if (relativistic):
                  newline = string.join(["CALL UT_R4DTCE_FILENAME('",filename,"',filename)"],"")
               else:
                  newline = string.join(["CALL UT_TCE_FILENAME('",filename,"',filename)"],"")
               newcode.statements.insert(0,newline)
               if (relativistic):
                  newline = string.join(["CALL R4D_CREATEFILE(filename,",d_c,",",size_c,")"],"")
               else:
                  newline = string.join(["CALL CREATEFILE(filename,",d_c,",",size_c,")"],"")
               newcode.statements.insert(0,newline)
               callee = child.contraction.tensors[0].utchemy(globaltargetindexes,name,relativistic)
               callees.add(callee)
               createfile = 0
            newcode.statements.insert(0,child.utchema(name,globaltargetindexes,callees,relativistic))
            if (child.contraction.tensors[1].type == "i"):
               if (relativistic):
                  newline = string.join(["CALL R4D_RECONCILEFILE(",d_a,",",size_a,")"],"")
               else:
                  newline = string.join(["CALL RECONCILEFILE(",d_a,",",size_a,")"],"")
               newcode.statements.insert(0,newline)
            if (d_b):
               if (child.contraction.tensors[2].type == "i"):
                  if (relativistic):
                     newline = string.join(["CALL R4D_RECONCILEFILE(",d_b,",",size_b,")"],"")
                  else:
                     newline = string.join(["CALL RECONCILEFILE(",d_b,",",size_b,")"],"")
                  newcode.statements.insert(0,newline)
            argument = string.join([d_a,",",a_offset],"")
            if (d_b):
               argument = string.join([argument,",",d_b,",",b_offset],"")
            argument = string.join([argument,",",d_c,",",c_offset],"")
            if (relativistic):
               newline = string.join(["CALL R4D_",name,"(",argument,")"],"")
            else:
               newline = string.join(["CALL ",name,"(",argument,")"],"")
            newcode.statements.insert(0,newline)
            if (child.contraction.tensors[1].type == "i"):
               if (relativistic):
                  newline = string.join(["CALL R4D_DELETEFILE(",d_a,")"],"")
               else:
                  newline = string.join(["CALL DELETEFILE(",d_a,")"],"")
               newcode.statements.insert(0,newline)
               newline = string.join(["DEALLOCATE(",a_offset,")"],"")
               newcode.statements.insert(0,newline)
            if (d_b):
               if (child.contraction.tensors[2].type == "i"):
                  if (relativistic):
                     newline = string.join(["CALL R4D_DELETEFILE(",d_b,")"],"")
                  else:
                     newline = string.join(["CALL DELETEFILE(",d_b,")"],"")
                  newcode.statements.insert(0,newline)
                  newline = string.join(["DEALLOCATE(",b_offset,")"],"")
                  newcode.statements.insert(0,newline)

      if (self.sisters):
         counter = 0
         for isister in range(len(self.sisters)-1,-1,-1):
            sister = self.sisters[isister]
            counter = counter + 1
            d_c = string.join(["d_j",repr(sister.tensors[0].label)],"")
            c_offset = string.join(["j",repr(sister.tensors[0].label),"_offset"],"")
            name = string.join([subroutinename,"_r_",repr(counter)],"")
            if (relativistic):
               newline = string.join(["CALL R4D_DELETEFILE(",d_c,")"],"")
            else:
               newline = string.join(["CALL DELETEFILE(",d_c,")"],"")
            newcode.statements.insert(0,newline)
            newline = string.join(["DEALLOCATE(",c_offset,")"],"")
            newcode.statements.insert(0,newline)
 
      newcode.reverse()
      return newcode

   def pythongen(self,filename="NONAME"):
      """Genrates a python code for debugging purposes"""

      pythoncode = []

      tensornames = self.tensorslist([])
      tensorargument = ""
      for tensorname in tensornames:
         tensorargument = string.join([tensorargument,",",tensorname],"")
      newline = string.join(["def ",filename,"(N,nall,nocc",tensorargument,"):"],"")
      pythoncode.append(newline)
      newline = " # This is a Python program generated by Tensor Contraction Engine v.1.0"
      pythoncode.append(newline)
      newline = " # (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)"
      pythoncode.append(newline)
 
      # copy of self will be reduced as we write the program
      selfcopy = OperationTree()
      selfcopy.contraction = copy.deepcopy(self.contraction)
      selfcopy.common = copy.deepcopy(self.common)
      selfcopy.children = copy.deepcopy(self.children)
 
      # loop over the tree
      newcode = selfcopy.pythongena(0)
      newcodeexpanded = expand(newcode)
      for newline in newcodeexpanded:
         pythoncode.append(newline)
         
      # dump the code to a file
      writetofile(pythoncode,string.join([filename,".py.out"],""))

   def pythongena(self,pointer=0):
      """Recursive subprogram called by pythongen"""

      pythoncode = []

      # check if we need to proceed
      if (not self.children):
         return pythoncode
      else:
         empty = 1
         for child in self.children:
            if (child.contraction.isoperation()):
               empty = 0
         if (empty):
            return pythoncode

      # loop over children
      zeroscratch = 1
      for child in self.children:
         if (child.contraction.isoperation()):

            # recursive pythongena() call
            pythoncode.insert(pointer,child.pythongena(pointer))
            pointer = len(pythoncode)
            
            # generate loops over target indexes
            indent = 1
            for index in child.contraction.tensors[0].indexes:
               spin = string.join(["spin",repr(index.label)],"")
               newline = string.join(["for ",spin," in range(2):"],"")
               newline = string.join([" "*indent,newline],"")
               pythoncode.insert(pointer,newline)
               pointer = pointer + 1
               indent = indent + 1
               if (index.type == 'hole'):
                  newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nocc[",spin,"]):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.insert(pointer,newline)
                  pointer = pointer + 1
                  indent = indent + 1
               elif (index.type == 'particle'):
                  newline = string.join(["for ",index.show()," in range(",spin,"*nall[0]+nocc[",spin,"],",spin,"*nall[0]+nall[",spin,"]):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.insert(pointer,newline)
                  pointer = pointer + 1
                  indent = indent + 1
               elif (index.type == 'general'):
                  newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nall[",spin,"]):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.insert(pointer,newline)
                  pointer = pointer + 1
                  indent = indent + 1
            if (zeroscratch):
               newline = string.join([child.contraction.tensors[0].pythongen(),"=0.0"],"")
               newline = string.join([" "*indent,newline],"")
               pythoncode.insert(pointer,newline)
               pointer = pointer + 1
               zeroscratch = 0

            # generate loops over common indexes
            if (child.contraction.summation):
               for index in child.contraction.summation.indexes:
                  spin = string.join(["spin",repr(index.label)],"")
                  newline = string.join(["for ",spin," in range(2):"],"")
                  newline = string.join([" "*indent,newline],"")
                  pythoncode.insert(pointer,newline)
                  pointer = pointer + 1
                  indent = indent + 1
                  if (index.type == 'hole'):
                     newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nocc[",spin,"]):"],"")
                     newline = string.join([" "*indent,newline],"")
                     pythoncode.insert(pointer,newline)
                     pointer = pointer + 1
                     indent = indent + 1
                  elif (index.type == 'particle'):
                     newline = string.join(["for ",index.show()," in range(",spin,"*nall[0]+nocc[",spin,"],",spin,"*nall[0]+nall[",spin,"]):"],"")
                     newline = string.join([" "*indent,newline],"")
                     pythoncode.insert(pointer,newline)
                     pointer = pointer + 1
                     indent = indent + 1
                  elif (index.type == 'general'):
                     newline = string.join(["for ",index.show()," in range(",spin,"*nall[0],",spin,"*nall[0]+nall[",spin,"]):"],"")
                     newline = string.join([" "*indent,newline],"")
                     pythoncode.insert(pointer,newline)
                     pointer = pointer + 1
                     indent = indent + 1
            newline = string.join([child.contraction.tensors[0].pythongen(),"=",child.contraction.tensors[0].pythongen(),"+",\
                                   "(",repr(child.contraction.factor),")"],"")
            for ntensor in range(len(child.contraction.tensors)):
               if (ntensor != 0):
                  tensor = child.contraction.tensors[ntensor]
                  newline = string.join([newline,"*",tensor.pythongen()],"")
            newline = string.join([" "*indent,newline],"")
            pythoncode.insert(pointer,newline)
            pointer = pointer + 1

      return pythoncode

   def fortran90(self,filename="NONAME",mode="nopermutation",excitation=[],deexcitation=[],intermediate=[],general=[]):
      """Genrates a partial Fortran90 code for debugging purposes"""
      # Mode = "permutation"   : writes a code which takes index permutation into account
      # Mode = "nopermutation" : writes a code without index permutation considered
      # Mode = "analysis"      : stdouts a plan of implementation with index permutation

      print " ... generating a Fortran90 code"
      print " "
      all = excitation+deexcitation+intermediate+general
      if (not all):
         types = self.tensortypes()
      else:
         types = [all,excitation,deexcitation,intermediate,general]
      for type in types[0]:
         if (type in types[1]):
            print " '"+type+"' is an excitaion tensor"
            if ((type == 'v') and (not all)):
               raise RuntimeError, "unusual naming convention"
            if ((type == 'f') and (not all)):
               raise RuntimeError, "unusual naming convention"
            if ((type == 'd') and (not all)):
               raise RuntimeError, "unusual naming convention"
         elif (type in types[2]):
            print " '"+type+"' is a de-excitaion tensor"
            if ((type == 'v') and (not all)):
               raise RuntimeError, "unusual naming convention"
            if ((type == 'f') and (not all)):
               raise RuntimeError, "unusual naming convention"
            if ((type == 'd') and (not all)):
               raise RuntimeError, "unusual naming convention"
         elif (type in types[3]):
            print " '"+type+"' is an intermediate tensor"
         elif (type in types[4]):
            print " '"+type+"' is a general tensor"
         else:
            raise RuntimeError, "unknown tensor type"
      print " "

      f90code = Code("Fortran90",filename)
      newline = "IMPLICIT NONE"
      f90code.add("headers",newline)

      # copy of self will be reduced as we write the program
      selfcopy = OperationTree()
      selfcopy.contraction = copy.deepcopy(self.contraction)
      selfcopy.common = copy.deepcopy(self.common)
      selfcopy.children = copy.deepcopy(self.children)

      # target indexes
      if (selfcopy.children[0].contraction.isoperation()):
         globaltargetindexes = copy.deepcopy(selfcopy.children[0].contraction.tensors[0].indexes)
      else:
         return "The tree top must be an addition"
      if (mode == "analysis"):
         show = "Global target indexes: "
         for index in globaltargetindexes:
            show = string.join([show,index.show()])
         print show
         print ""
 
      # loop over the tree
      f90code.statements.insert(f90code.pointer,selfcopy.fortran90a(globaltargetindexes,types,mode))
      
      # add an antisymmetrizer (only for the target intermediate)
      if (mode == "nopermutation"):
         f90code.statements.append(selfcopy.children[0].contraction.tensors[0].fortran90x(types))
      
      # close the subroutine
      newline = "RETURN"
      f90code.statements.append(newline)
      newline = "END SUBROUTINE"
      f90code.statements.append(newline)

      # headers
      f90code.add("arguments","N")
      f90code.add("integers","N")
      f90code.add("arguments","nocc")
      f90code.add("integers","nocc")
      if (mode == "nopermutation"):
         f90code.add("doubles","TMP")
      for n in self.tensorslist([]):
         f90code.add("arguments",n)
         f90code.add("doublearrays",n)
         if (mode == "permutation"):
            f90code.add("arguments",n+"e")
            f90code.add("doublearrays",n+"e")

      # dump the code to a file
      f90code = f90code.expand()
      f90code.sortarguments()
      f90code.writetofile(filename)
      if (mode == "analysis"):
         return "No Fortran code is dumped"
      else:
         return f90code

   def fortran90a(self,globaltargetindexes,types,mode="nopermutation"):
      """Recursive subprogram called by fortran90"""

      newcode = Code("Fortran90","")

      # check if we need to proceed
      if (not self.children):
         return newcode
      else:
         empty = 1
         for child in self.children:
            if (child.contraction.isoperation()):
               empty = 0
         if (empty):
            return newcode

      # loop over children
      zeroscratch = 1
      for child in self.children:
         if (child.contraction.isoperation()):

            # recursive fortran90a() call
            newcode.statements.insert(newcode.pointer,child.fortran90a(globaltargetindexes,types,mode))
            newcode.pointer = len(newcode.statements)
            
            # Tensor 1
            superglobalzero = []
            subglobalzero = []
            superlocalzero = []
            sublocalzero = []
            for nindex in range(len(child.contraction.tensors[0].indexes)/2):
               index = child.contraction.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  superglobalzero.append(index)
               else:
                  superlocalzero.append(index)
            for nindex in range(len(child.contraction.tensors[0].indexes)/2,len(child.contraction.tensors[0].indexes)):
               index = child.contraction.tensors[0].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  subglobalzero.append(index)
               else:
                  sublocalzero.append(index)

            # Tensor 2
            superglobalone = []
            subglobalone = []
            superlocalone = []
            sublocalone = []
            supercommonone = []
            subcommonone = []
            for nindex in range(len(child.contraction.tensors[1].indexes)/2):
               index = child.contraction.tensors[1].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  superglobalone.append(index)
               elif (child.contraction.summation):
                  if (index.isin(child.contraction.summation.indexes)):
                     supercommonone.append(index)
                  else:
                     superlocalone.append(index)
               else:
                  superlocalone.append(index)
            for nindex in range(len(child.contraction.tensors[1].indexes)/2,len(child.contraction.tensors[1].indexes)):
               index = child.contraction.tensors[1].indexes[nindex]
               if (index.isin(globaltargetindexes)):
                  subglobalone.append(index)
               elif (child.contraction.summation):
                  if (index.isin(child.contraction.summation.indexes)):
                     subcommonone.append(index)
                  else:
                     sublocalone.append(index)
               else:
                  sublocalone.append(index)

            # Tensor 3
            superglobaltwo = []
            subglobaltwo = []
            superlocaltwo = []
            sublocaltwo = []
            supercommontwo = []
            subcommontwo = []
            if (len(child.contraction.tensors) > 2):
               for nindex in range(len(child.contraction.tensors[2].indexes)/2):
                  index = child.contraction.tensors[2].indexes[nindex]
                  if (index.isin(globaltargetindexes)):
                     superglobaltwo.append(index)
                  elif (child.contraction.summation):
                     if (index.isin(child.contraction.summation.indexes)):
                        supercommontwo.append(index)
                     else:
                        superlocaltwo.append(index)
                  else:
                     superlocaltwo.append(index)
               for nindex in range(len(child.contraction.tensors[2].indexes)/2,len(child.contraction.tensors[2].indexes)):
                  index = child.contraction.tensors[2].indexes[nindex]
                  if (index.isin(globaltargetindexes)):
                     subglobaltwo.append(index)
                  elif (child.contraction.summation):
                     if (index.isin(child.contraction.summation.indexes)):
                        subcommontwo.append(index)
                     else:
                        sublocaltwo.append(index)
                  else:
                     sublocaltwo.append(index)
            if (len(supercommonone) > len(subcommontwo)):
               supercommon = supercommonone
            else:
               supercommon = subcommontwo
            if (len(subcommonone) > len(supercommontwo)):
               subcommon = subcommonone
            else:
               subcommon = supercommontwo
# DEEXCITATION EXTENTION FROM HERE ...
            superlocalparticleone  = []
            superlocalparticletwo  = []
            superlocalholeone  = []
            superlocalholetwo  = []
            supercommonparticleone  = []
            supercommonparticletwo  = []
            supercommonholeone  = []
            supercommonholetwo  = []
            for index in superlocalone:
               if (index.isparticle()):
                  superlocalparticleone.append(index)
               elif (index.ishole()):
                  superlocalholeone.append(index)
            for index in superlocaltwo:
               if (index.isparticle()):
                  superlocalparticletwo.append(index)
               elif (index.ishole()):
                  superlocalholetwo.append(index)
            for index in supercommonone:
               if (index.isparticle()):
                  supercommonparticleone.append(index)
               elif (index.ishole()):
                  supercommonholeone.append(index)
            for index in supercommontwo:
               if (index.isparticle()):
                  supercommonparticletwo.append(index)
               elif (index.ishole()):
                  supercommonholetwo.append(index)
            superlocalparticlezero = sortindexes( superlocalparticleone + superlocalparticletwo )
            superlocalholezero = sortindexes( superlocalholeone + superlocalholetwo )
            sublocalparticleone  = []
            sublocalparticletwo  = []
            sublocalholeone  = []
            sublocalholetwo  = []
            subcommonparticleone  = []
            subcommonparticletwo  = []
            subcommonholeone  = []
            subcommonholetwo  = []
            for index in sublocalone:
               if (index.isparticle()):
                  sublocalparticleone.append(index)
               elif (index.ishole()):
                  sublocalholeone.append(index)
            for index in sublocaltwo:
               if (index.isparticle()):
                  sublocalparticletwo.append(index)
               elif (index.ishole()):
                  sublocalholetwo.append(index)
            for index in subcommonone:
               if (index.isparticle()):
                  subcommonparticleone.append(index)
               elif (index.ishole()):
                  subcommonholeone.append(index)
            for index in subcommontwo:
               if (index.isparticle()):
                  subcommonparticletwo.append(index)
               elif (index.ishole()):
                  subcommonholetwo.append(index)
            sublocalparticlezero = sortindexes( sublocalparticleone + sublocalparticletwo )
            sublocalholezero = sortindexes( sublocalholeone + sublocalholetwo )
            supercommonhole = []
            supercommonparticle = []
            subcommonhole = []
            subcommonparticle = []
            for index in supercommon:
               if (index.ishole()):
                  supercommonhole.append(index)
               elif (index.isparticle()):
                  supercommonparticle.append(index)
               else:
                  raise RuntimeError, "a general common index appeared"
            for index in subcommon:
               if (index.ishole()):
                  subcommonhole.append(index)
               elif (index.isparticle()):
                  subcommonparticle.append(index)
               else:
                  raise RuntimeError, "a general common index appeared"
# ... TO HERE

            if (mode == "analysis"):

               # Zero scratch?
               if (zeroscratch):
                  print child.contraction.tensors[0].show(),"will be zeroscratched"
                  zeroscratch = 0

               # Structure of tensor 1
               tensorzerocompressed = 0
               show = "Storage of tensor 1:"
               if (superglobalzero):
                  show = string.join([show,"["])
                  for nindex in range(len(superglobalzero)):
                     index = superglobalzero[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (subglobalzero):
                  show = string.join([show,"["])
                  for nindex in range(len(subglobalzero)):
                     index = subglobalzero[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (superlocalzero):
                  show = string.join([show,"["])
                  for nindex in range(len(superlocalzero)):
                     index = superlocalzero[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (sublocalzero):
                  show = string.join([show,"["])
                  for nindex in range(len(sublocalzero)):
                     index = sublocalzero[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if ((superglobalone) and (superglobaltwo)):
                  tensorzerocompressed = 1
                  show = string.join([show,"Compress ["])
                  for nindex in range(len(superglobalone)):
                     index = superglobalone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"] ["])
                  for nindex in range(len(superglobaltwo)):
                     index = superglobaltwo[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if ((subglobalone) and (subglobaltwo)):
                  tensorzerocompressed = 1
                  show = string.join([show,"Compress ["])
                  for nindex in range(len(subglobalone)):
                     index = subglobalone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"] ["])
                  for nindex in range(len(subglobaltwo)):
                     index = subglobaltwo[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               print show

               # Structure of tensor 2
               tensoroneexpanded = 0
               show = "Storage of tensor 2:"
               if (child.contraction.tensors[1].type in types[3]):
                  if (superglobalone):
                     show = string.join([show,"["])
                     for nindex in range(len(superglobalone)):
                        index = superglobalone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (subglobalone):
                     show = string.join([show,"["])
                     for nindex in range(len(subglobalone)):
                        index = subglobalone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
               if (child.contraction.tensors[1].type in types[3]):
                  superremainderone = sortindexes(superlocalone + supercommonone)
               else:
                  superremainderone = sortindexes(superglobalone + superlocalone + supercommonone)
               if (superremainderone):
                  show = string.join([show,"["])
                  for nindex in range(len(superremainderone)):
                     index = superremainderone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (child.contraction.tensors[1].type in types[3]):
                  subremainderone = sublocalone + subcommonone
               else:
                  subremainderone = subglobalone + sublocalone + subcommonone
               if (subremainderone):
                  show = string.join([show,"["])
                  for nindex in range(len(subremainderone)):
                     index = subremainderone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (child.contraction.tensors[1].type in types[3]):
                  if (superlocalparticleone and supercommonparticleone):
                     tensoroneexpanded = 1
                     show = string.join([show,"Expand ["])
                     for nindex in range(len(superlocalparticleone)):
                        index = superlocalparticleone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"] ["])
                     for nindex in range(len(supercommonparticleone)):
                        index = supercommonparticleone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (superlocalholeone and supercommonholeone):
                     tensoroneexpanded = 1
                     show = string.join([show,"Expand ["])
                     for nindex in range(len(superlocalholeone)):
                        index = superlocalholeone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"] ["])
                     for nindex in range(len(supercommonholeone)):
                        index = supercommonholeone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
               elif ((superglobalone and superlocalone) or \
                     (superlocalone  and supercommonone) or \
                     (supercommonone and superglobalone)):
                  tensoroneexpanded = 1
                  show = string.join([show,"Expand"])
                  if (superglobalone):
                     show = string.join([show,"["])
                     for nindex in range(len(superglobalone)):
                        index = superglobalone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (superlocalone):
                     show = string.join([show,"["])
                     for nindex in range(len(superlocalone)):
                        index = superlocalone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (supercommonone):
                     show = string.join([show,"["])
                     for nindex in range(len(supercommonone)):
                        index = supercommonone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
               if (child.contraction.tensors[1].type in types[3]):
                  if ((sublocalparticleone) and (subcommonparticleone)):
                     tensoroneexpanded = 1
                     show = string.join([show,"Expand ["])
                     for nindex in range(len(sublocalparticleone)):
                        index = sublocalparticleone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"] ["])
                     for nindex in range(len(subcommonparticleone)):
                        index = subcommonparticleone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if ((sublocalholeone) and (subcommonholeone)):
                     tensoroneexpanded = 1
                     show = string.join([show,"Expand ["])
                     for nindex in range(len(sublocalholeone)):
                        index = sublocalholeone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"] ["])
                     for nindex in range(len(subcommonholeone)):
                        index = subcommonholeone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
               elif ((subglobalone and sublocalone) or \
                     (sublocalone  and subcommonone) or \
                     (subcommonone and subglobalone)):
                  tensoroneexpanded = 1
                  show = string.join([show,"Expand"])
                  if (subglobalone):
                     show = string.join([show,"["])
                     for nindex in range(len(subglobalone)):
                        index = subglobalone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (sublocalone):
                     show = string.join([show,"["])
                     for nindex in range(len(sublocalone)):
                        index = sublocalone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (subcommonone):
                     show = string.join([show,"["])
                     for nindex in range(len(subcommonone)):
                        index = subcommonone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
               print show

               # Structure of tensor 3
               if (len(child.contraction.tensors) > 2):
                  tensortwoexpanded = 0
                  show = "Storage of tensor 3:"
                  if (child.contraction.tensors[2].type in types[3]):
                     if (superglobaltwo):
                        show = string.join([show,"["])
                        for nindex in range(len(superglobaltwo)):
                           index = superglobaltwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (subglobaltwo):
                        show = string.join([show,"["])
                        for nindex in range(len(subglobaltwo)):
                           index = subglobaltwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                  if (child.contraction.tensors[2].type in types[3]):
                     superremaindertwo = sortindexes(superlocaltwo + supercommontwo)
                  else:
                     superremaindertwo = sortindexes(superglobaltwo + superlocaltwo + supercommontwo)
                  if (superremaindertwo):
                     show = string.join([show,"["])
                     for nindex in range(len(superremaindertwo)):
                        index = superremaindertwo[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (child.contraction.tensors[2].type in types[3]):
                     subremaindertwo = sublocaltwo + subcommontwo
                  else:
                     subremaindertwo = subglobaltwo + sublocaltwo + subcommontwo
                  if (subremaindertwo):
                     show = string.join([show,"["])
                     for nindex in range(len(subremaindertwo)):
                        index = subremaindertwo[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                  if (child.contraction.tensors[2].type in types[3]):
                     if (superlocalparticletwo and supercommonparticletwo):
                        tensortwoexpanded = 1
                        show = string.join([show,"Expand ["])
                        for nindex in range(len(superlocalparticletwo)):
                           index = superlocalparticletwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"] ["])
                        for nindex in range(len(supercommonparticletwo)):
                           index = supercommonparticletwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (superlocalholetwo and supercommonholetwo):
                        tensortwoexpanded = 1
                        show = string.join([show,"Expand ["])
                        for nindex in range(len(superlocalholetwo)):
                           index = superlocalholetwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"] ["])
                        for nindex in range(len(supercommonholetwo)):
                           index = supercommonholetwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                  elif ((superglobaltwo and superlocaltwo) or \
                        (superlocaltwo  and supercommontwo) or \
                        (supercommontwo and superglobaltwo)):
                     tensortwoexpanded = 1
                     show = string.join([show,"Expand"])
                     if (superglobaltwo):
                        show = string.join([show,"["])
                        for nindex in range(len(superglobaltwo)):
                           index = superglobaltwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (superlocaltwo):
                        show = string.join([show,"["])
                        for nindex in range(len(superlocaltwo)):
                           index = superlocaltwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (supercommontwo):
                        show = string.join([show,"["])
                        for nindex in range(len(supercommontwo)):
                           index = supercommontwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                  if (child.contraction.tensors[2].type in types[3]):
                     if (sublocalparticletwo and subcommonparticletwo):
                        tensortwoexpanded = 1
                        show = string.join([show,"Expand ["])
                        for nindex in range(len(sublocalparticletwo)):
                           index = sublocalparticletwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"] ["])
                        for nindex in range(len(subcommonparticletwo)):
                           index = subcommonparticletwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (sublocalholetwo and subcommonholetwo):
                        tensortwoexpanded = 1
                        show = string.join([show,"Expand ["])
                        for nindex in range(len(sublocalholetwo)):
                           index = sublocalholetwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"] ["])
                        for nindex in range(len(subcommonholetwo)):
                           index = subcommonholetwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                  elif ((subglobaltwo and sublocaltwo) or \
                        (sublocaltwo  and subcommontwo) or \
                        (subcommontwo and subglobaltwo)):
                     tensortwoexpanded = 1
                     show = string.join([show,"Expand"])
                     if (subglobaltwo):
                        show = string.join([show,"["])
                        for nindex in range(len(subglobaltwo)):
                           index = subglobaltwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (sublocaltwo):
                        show = string.join([show,"["])
                        for nindex in range(len(sublocaltwo)):
                           index = sublocaltwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                     if (subcommontwo):
                        show = string.join([show,"["])
                        for nindex in range(len(subcommontwo)):
                           index = subcommontwo[nindex]
                           if (nindex > 0):
                              show = string.join([show,"<"])
                           show = string.join([show,index.show()])
                        show = string.join([show,"]"])
                  print show
# DEEXCITATION EXTENSION FROM HERE ...
                  if (superlocalparticleone and superlocalparticletwo):
                     show = "["
                     for nindex in range(len(superlocalparticleone)):
                        index = superlocalparticleone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]["])
                     for nindex in range(len(superlocalparticletwo)):
                        index = superlocalparticletwo[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"] Compress ["])
                     for nindex in range(len(superlocalparticlezero)):
                        index = superlocalparticlezero[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                     print show
                  if (sublocalholeone and sublocalholetwo):
                     show = "["
                     for nindex in range(len(sublocalholeone)):
                        index = sublocalholeone[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]["])
                     for nindex in range(len(sublocalholetwo)):
                        index = sublocalholetwo[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"] Compress ["])
                     for nindex in range(len(sublocalholezero)):
                        index = sublocalholezero[nindex]
                        if (nindex > 0):
                           show = string.join([show,"<"])
                        show = string.join([show,index.show()])
                     show = string.join([show,"]"])
                     print show
                  if ((superlocalparticleone and superlocalparticletwo) or \
                      (sublocalholeone and sublocalholetwo)):
                     newfactor = Factor([1.0],[[]])
                     newfactor = permutationsoffoursets(newfactor,superlocalparticleone, superlocalparticletwo, \
                                                                  sublocalholeone,       sublocalholetwo, \
                                                                  superlocalparticlezero,sublocalholezero)
                     newfactor = newfactor.normalize()
                     print newfactor
# ... TO HERE

               # Summation indexes
               show = "Summation composite indexes:"
               if (supercommonhole):
                  show = string.join([show,"["])
                  for nindex in range(len(supercommonhole)):
                     index = supercommonhole[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (supercommonparticle):
                  show = string.join([show,"["])
                  for nindex in range(len(supercommonparticle)):
                     index = supercommonparticle[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (subcommonhole):
                  show = string.join([show,"["])
                  for nindex in range(len(subcommonhole)):
                     index = subcommonhole[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (subcommonparticle):
                  show = string.join([show,"["])
                  for nindex in range(len(subcommonparticle)):
                     index = subcommonparticle[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if ((not supercommon) and (not subcommon)):
                  show = string.join([show,"none"])
               else:
                  factor = factorial(len(supercommonhole)) * factorial(len(supercommonparticle)) \
                         * factorial(len(subcommonhole)) * factorial(len(subcommonparticle))
                  show = string.join([show,"with a factor of",repr(factor)])
               print show

               # Target indexes
               show = "Target composite indexes:"
               if (superglobalone):
                  show = string.join([show,"["])
                  for nindex in range(len(superglobalone)):
                     index = superglobalone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (subglobalone):
                  show = string.join([show,"["])
                  for nindex in range(len(subglobalone)):
                     index = subglobalone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (superglobaltwo):
                  show = string.join([show,"["])
                  for nindex in range(len(superglobaltwo)):
                     index = superglobaltwo[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (subglobaltwo):
                  show = string.join([show,"["])
                  for nindex in range(len(subglobaltwo)):
                     index = subglobaltwo[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (superlocalone):
                  show = string.join([show,"["])
                  for nindex in range(len(superlocalone)):
                     index = superlocalone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (sublocalone):
                  show = string.join([show,"["])
                  for nindex in range(len(sublocalone)):
                     index = sublocalone[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (superlocaltwo):
                  show = string.join([show,"["])
                  for nindex in range(len(superlocaltwo)):
                     index = superlocaltwo[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               if (sublocaltwo):
                  show = string.join([show,"["])
                  for nindex in range(len(sublocaltwo)):
                     index = sublocaltwo[nindex]
                     if (nindex > 0):
                        show = string.join([show,"<"])
                     show = string.join([show,index.show()])
                  show = string.join([show,"]"])
               print show
               print child.contraction
               print ""

            # zero scratch
            if ((zeroscratch) and (mode == "nopermutation")):
               for index in child.contraction.tensors[0].indexes:
                  newcode.insertdoloop(index)
               newdbl = child.contraction.tensors[0].fortran90(types)
               newline = string.join([newdbl,"=0.0d0"],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = len(newcode.statements)
               zeroscratch = 0

            expanded = [0,0,0]
            # expand index ranges of tensor 1
            if (mode == "permutation"):
               if (child.contraction.tensors[1].type in types[3]):
                  if ((superlocalone and supercommonone) or (sublocalone and subcommonone)):
                     expanded[1] = 1
                     super = sortindexes(superlocalone + supercommonone)
                     sub = sortindexes(sublocalone + subcommonone)
                     for index in superglobalone:
                        newcode.insertdoloop(index)
                     if (len(superglobalone) > 1):
                        newcode.insertif(superglobalone,1)
                     for index in super:
                        newcode.insertdoloop(index)
                     if (len(super) > 1):
                        newcode.insertif(super,1)
                     for index in subglobalone:
                        newcode.insertdoloop(index)
                     if (len(subglobalone) > 1):
                        newcode.insertif(subglobalone,1)
                     for index in sub:
                        newcode.insertdoloop(index)
                     if (len(sub) > 1):
                        newcode.insertif(sub,1)
                     if (superlocalone and supercommonone):
                        nsuperpermutations = factorial(len(super))
                        superpermutations = permutationwithparity(len(super))
                     else:
                        nsuperpermutations = 1
                     if (sublocalone and subcommonone):
                        nsubpermutations = factorial(len(sub))
                        subpermutations = permutationwithparity(len(sub))
                     else:
                        nsubpermutations = 1
                     for nsuperpermutation in range(nsuperpermutations):
                        for nsubpermutation in range(nsubpermutations):
                           permutation = super + sub
                           parity = 1
                           rejected = 0
                           if (superlocalone and supercommonone):
                              parity = parity * superpermutations[nsuperpermutation][0]
                              for nindex in range(1,len(super)+1):
                                 permutation.append(super[superpermutations[nsuperpermutation][nindex]-1])
                              if (len(superlocalone) > 1):
                                 for nindexa in range(len(superlocalone)):
                                    for nindexb in range(len(superlocalone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = superlocalone[nindexa]
                                       indexb = superlocalone[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                              if (len(supercommonone) > 1):
                                 for nindexa in range(len(supercommonone)):
                                    for nindexb in range(len(supercommonone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = supercommonone[nindexa]
                                       indexb = supercommonone[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + super
                           if (sublocalone and subcommonone):
                              parity = parity * subpermutations[nsubpermutation][0]
                              for nindex in range(1,len(sub)+1):
                                 permutation.append(sub[subpermutations[nsubpermutation][nindex]-1])
                              if (len(sublocalone) > 1):
                                 for nindexa in range(len(sublocalone)):
                                    for nindexb in range(len(sublocalone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = sublocalone[nindexa]
                                       indexb = sublocalone[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                              if (len(subcommonone) > 1):
                                 for nindexa in range(len(subcommonone)):
                                    for nindexb in range(len(subcommonone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = subcommonone[nindexa]
                                       indexb = subcommonone[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + sub
                           for nindexa in range(len(permutation)/2,len(permutation)):
                              for nindexb in range(len(permutation)/2,len(permutation)):
                                 if (nindexa >= nindexb):
                                    continue
                                 indexa = permutation[nindexa]
                                 indexb = permutation[nindexb]
                                 if ((indexa.isin(super)) and (indexb.isin(super)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                                 if ((indexa.isin(sub)) and (indexb.isin(sub)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                           if (not rejected):
                              if (parity == 1):
                                 sign = " + "
                              else:
                                 sign = " - "
                              newdbl = child.contraction.tensors[1].fortran90(types,permutation,0,"e")
                              newline = string.join([newdbl,"=",sign,child.contraction.tensors[1].fortran90(types)],"")
                              newcode.statements.insert(newcode.pointer,newline)
                              newcode.pointer = newcode.pointer + 1
                     newcode.pointer = len(newcode.statements)
               else:
                  if ((superglobalone and superlocalone) or (superglobalone and supercommonone) or (superlocalone and supercommonone) \
                   or (subglobalone and sublocalone) or (subglobalone and subcommonone) or (sublocalone and subcommonone)):
                     expanded[1] = 1
                     super = sortindexes(superglobalone + superlocalone + supercommonone)
                     sub = sortindexes(subglobalone + sublocalone + subcommonone)
                     for index in superglobalone:
                        newcode.insertdoloop(index)
                     for index in superlocalone:
                        newcode.insertdoloop(index)
                     for index in supercommonone:
                        newcode.insertdoloop(index)
                     if (len(super) > 1):
                        newcode.insertif(super,1)
                     for index in subglobalone:
                        newcode.insertdoloop(index)
                     for index in sublocalone:
                        newcode.insertdoloop(index)
                     for index in subcommonone:
                        newcode.insertdoloop(index)
                     if (len(sub) > 1):
                        newcode.insertif(sub,1)
                     if ((superglobalone and superlocalone) or (superglobalone and supercommonone) or (superlocalone and supercommonone)):
                        nsuperpermutations = factorial(len(super))
                        superpermutations = permutationwithparity(len(super))
                     else:
                        nsuperpermutations = 1
                     if ((subglobalone and sublocalone) or (subglobalone and subcommonone) or (sublocalone and subcommonone)):
                        nsubpermutations = factorial(len(sub))
                        subpermutations = permutationwithparity(len(sub))
                     else:
                        nsubpermutations = 1
                     for nsuperpermutation in range(nsuperpermutations):
                        for nsubpermutation in range(nsubpermutations):
                           permutation = super + sub
                           parity = 1
                           rejected = 0
                           if ((superglobalone and superlocalone) or \
                               (superglobalone and supercommonone) or \
                               (superlocalone and supercommonone)):
                              parity = parity * superpermutations[nsuperpermutation][0]
                              for nindex in range(1,len(super)+1):
                                 permutation.append(super[superpermutations[nsuperpermutation][nindex]-1])
                              if (len(superglobalone) > 1):
                                 for nindexa in range(len(superglobalone)):
                                    for nindexb in range(len(superglobalone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = superglobalone[nindexa]
                                       indexb = superglobalone[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                              if (len(superlocalone) > 1):
                                 for nindexa in range(len(superlocalone)):
                                    for nindexb in range(len(superlocalone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = superlocalone[nindexa]
                                       indexb = superlocalone[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                              if (len(supercommonone) > 1):
                                 for nindexa in range(len(supercommonone)):
                                    for nindexb in range(len(supercommonone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = supercommonone[nindexa]
                                       indexb = supercommonone[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + super
                           if ((subglobalone and sublocalone) or \
                               (subglobalone and subcommonone) or \
                               (sublocalone and subcommonone)):
                              parity = parity * subpermutations[nsubpermutation][0]
                              for nindex in range(1,len(sub)+1):
                                 permutation.append(sub[subpermutations[nsubpermutation][nindex]-1])
                              if (len(subglobalone) > 1):
                                 for nindexa in range(len(subglobalone)):
                                    for nindexb in range(len(subglobalone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = subglobalone[nindexa]
                                       indexb = subglobalone[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                              if (len(sublocalone) > 1):
                                 for nindexa in range(len(sublocalone)):
                                    for nindexb in range(len(sublocalone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = sublocalone[nindexa]
                                       indexb = sublocalone[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                              if (len(subcommonone) > 1):
                                 for nindexa in range(len(subcommonone)):
                                    for nindexb in range(len(subcommonone)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = subcommonone[nindexa]
                                       indexb = subcommonone[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + sub
                           for nindexa in range(len(permutation)/2,len(permutation)):
                              for nindexb in range(len(permutation)/2,len(permutation)):
                                 if (nindexa >= nindexb):
                                    continue
                                 indexa = permutation[nindexa]
                                 indexb = permutation[nindexb]
                                 if ((indexa.isin(super)) and (indexb.isin(super)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                                 if ((indexa.isin(sub)) and (indexb.isin(sub)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                           if (not rejected):
                              if (parity == 1):
                                 sign = " + "
                              else:
                                 sign = " - "
                              newdbl = child.contraction.tensors[1].fortran90(types,permutation,0,"e")
                              newline = string.join([newdbl,"=",sign,child.contraction.tensors[1].fortran90(types)],"")
                              newcode.statements.insert(newcode.pointer,newline)
                              newcode.pointer = newcode.pointer + 1
                     newcode.pointer = len(newcode.statements)

            # expand index ranges of tensor 2
            if ((mode == "permutation") and (len(child.contraction.tensors) > 2)):
               if (child.contraction.tensors[2].type in types[3]):
                  if ((superlocaltwo and supercommontwo) or (sublocaltwo and subcommontwo)):
                     expanded[2] = 1
                     super = sortindexes(superlocaltwo + supercommontwo)
                     sub = sortindexes(sublocaltwo + subcommontwo)
                     for index in superglobaltwo:
                        newcode.insertdoloop(index)
                     if (len(superglobaltwo) > 1):
                        newcode.insertif(superglobaltwo,1)
                     for index in super:
                        newcode.insertdoloop(index)
                     if (len(super) > 1):
                        newcode.insertif(super,1)
                     for index in subglobaltwo:
                        newcode.insertdoloop(index)
                     if (len(subglobaltwo) > 1):
                        newcode.insertif(subglobaltwo,1)
                     for index in sub:
                        newcode.insertdoloop(index)
                     if (len(sub) > 1):
                        newcode.insertif(sub,1)
                     if (superlocaltwo and supercommontwo):
                        nsuperpermutations = factorial(len(super))
                        superpermutations = permutationwithparity(len(super))
                     else:
                        nsuperpermutations = 1
                     if (sublocaltwo and subcommontwo):
                        nsubpermutations = factorial(len(sub))
                        subpermutations = permutationwithparity(len(sub))
                     else:
                        nsubpermutations = 1
                     for nsuperpermutation in range(nsuperpermutations):
                        for nsubpermutation in range(nsubpermutations):
                           permutation = super + sub
                           parity = 1
                           rejected = 0
                           if (superlocaltwo and supercommontwo):
                              parity = parity * superpermutations[nsuperpermutation][0]
                              for nindex in range(1,len(super)+1):
                                 permutation.append(super[superpermutations[nsuperpermutation][nindex]-1])
                              if (len(superlocaltwo) > 1):
                                 for nindexa in range(len(superlocaltwo)):
                                    for nindexb in range(len(superlocaltwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = superlocaltwo[nindexa]
                                       indexb = superlocaltwo[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                              if (len(supercommontwo) > 1):
                                 for nindexa in range(len(supercommontwo)):
                                    for nindexb in range(len(supercommontwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = supercommontwo[nindexa]
                                       indexb = supercommontwo[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + super
                           if (sublocaltwo and subcommontwo):
                              parity = parity * subpermutations[nsubpermutation][0]
                              for nindex in range(1,len(sub)+1):
                                 permutation.append(sub[subpermutations[nsubpermutation][nindex]-1])
                              if (len(sublocaltwo) > 1):
                                 for nindexa in range(len(sublocaltwo)):
                                    for nindexb in range(len(sublocaltwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = sublocaltwo[nindexa]
                                       indexb = sublocaltwo[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                              if (len(subcommontwo) > 1):
                                 for nindexa in range(len(subcommontwo)):
                                    for nindexb in range(len(subcommontwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = subcommontwo[nindexa]
                                       indexb = subcommontwo[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + sub
                           for nindexa in range(len(permutation)/2,len(permutation)):
                              for nindexb in range(len(permutation)/2,len(permutation)):
                                 if (nindexa >= nindexb):
                                    continue
                                 indexa = permutation[nindexa]
                                 indexb = permutation[nindexb]
                                 if ((indexa.isin(super)) and (indexb.isin(super)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                                 if ((indexa.isin(sub)) and (indexb.isin(sub)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                           if (not rejected):
                              if (parity == 1):
                                 sign = " + "
                              else:
                                 sign = " - "
                              newdbl = child.contraction.tensors[2].fortran90(types,permutation,0,"e")
                              newline = string.join([newdbl,"=",sign,child.contraction.tensors[2].fortran90(types)],"")
                              newcode.statements.insert(newcode.pointer,newline)
                              newcode.pointer = newcode.pointer + 1
                     newcode.pointer = len(newcode.statements)
               else:
                  if ((superglobaltwo and superlocaltwo) or (superglobaltwo and supercommontwo) or (superlocaltwo and supercommontwo) \
                   or (subglobaltwo and sublocaltwo) or (subglobaltwo and subcommontwo) or (sublocaltwo and subcommontwo)):
                     expanded[2] = 1
                     super = sortindexes(superglobaltwo + superlocaltwo + supercommontwo)
                     sub = sortindexes(subglobaltwo + sublocaltwo + subcommontwo)
                     for index in superglobaltwo:
                        newcode.insertdoloop(index)
                     for index in superlocaltwo:
                        newcode.insertdoloop(index)
                     for index in supercommontwo:
                        newcode.insertdoloop(index)
                     if (len(super) > 1):
                        newcode.insertif(super,1)
                     for index in subglobaltwo:
                        newcode.insertdoloop(index)
                     for index in sublocaltwo:
                        newcode.insertdoloop(index)
                     for index in subcommontwo:
                        newcode.insertdoloop(index)
                     if (len(sub) > 1):
                        newcode.insertif(sub,1)
                     if ((superglobaltwo and superlocaltwo) or (superglobaltwo and supercommontwo) or (superlocaltwo and supercommontwo)):
                        nsuperpermutations = factorial(len(super))
                        superpermutations = permutationwithparity(len(super))
                     else:
                        nsuperpermutations = 1
                     if ((subglobaltwo and sublocaltwo) or (subglobaltwo and subcommontwo) or (sublocaltwo and subcommontwo)):
                        nsubpermutations = factorial(len(sub))
                        subpermutations = permutationwithparity(len(sub))
                     else:
                        nsubpermutations = 1
                     for nsuperpermutation in range(nsuperpermutations):
                        for nsubpermutation in range(nsubpermutations):
                           permutation = super + sub
                           parity = 1
                           rejected = 0
                           if ((superglobaltwo and superlocaltwo) or \
                               (superglobaltwo and supercommontwo) or \
                               (superlocaltwo and supercommontwo)):
                              parity = parity * superpermutations[nsuperpermutation][0]
                              for nindex in range(1,len(super)+1):
                                 permutation.append(super[superpermutations[nsuperpermutation][nindex]-1])
                              if (len(superglobaltwo) > 1):
                                 for nindexa in range(len(superglobaltwo)):
                                    for nindexb in range(len(superglobaltwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = superglobaltwo[nindexa]
                                       indexb = superglobaltwo[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                              if (len(superlocaltwo) > 1):
                                 for nindexa in range(len(superlocaltwo)):
                                    for nindexb in range(len(superlocaltwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = superlocaltwo[nindexa]
                                       indexb = superlocaltwo[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                              if (len(supercommontwo) > 1):
                                 for nindexa in range(len(supercommontwo)):
                                    for nindexb in range(len(supercommontwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = supercommontwo[nindexa]
                                       indexb = supercommontwo[nindexb]
                                       for nindex in range(len(super)):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(super+sub)+nindexc].isgreaterthan(permutation[len(super+sub)+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + super
                           if ((subglobaltwo and sublocaltwo) or \
                               (subglobaltwo and subcommontwo) or \
                               (sublocaltwo and subcommontwo)):
                              parity = parity * subpermutations[nsubpermutation][0]
                              for nindex in range(1,len(sub)+1):
                                 permutation.append(sub[subpermutations[nsubpermutation][nindex]-1])
                              if (len(subglobaltwo) > 1):
                                 for nindexa in range(len(subglobaltwo)):
                                    for nindexb in range(len(subglobaltwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = subglobaltwo[nindexa]
                                       indexb = subglobaltwo[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                              if (len(sublocaltwo) > 1):
                                 for nindexa in range(len(sublocaltwo)):
                                    for nindexb in range(len(sublocaltwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = sublocaltwo[nindexa]
                                       indexb = sublocaltwo[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                              if (len(subcommontwo) > 1):
                                 for nindexa in range(len(subcommontwo)):
                                    for nindexb in range(len(subcommontwo)):
                                       if (nindexa >= nindexb):
                                          continue
                                       indexa = subcommontwo[nindexa]
                                       indexb = subcommontwo[nindexb]
                                       for nindex in range(len(permutation)/2):
                                          if (permutation[nindex].isidenticalto(indexa)):
                                             nindexc = nindex
                                          elif (permutation[nindex].isidenticalto(indexb)):
                                             nindexd = nindex
                                       if (permutation[len(permutation)/2+nindexc].isgreaterthan(permutation[len(permutation)/2+nindexd])):
                                          rejected = 1
                           else:
                              permutation = permutation + sub
                           for nindexa in range(len(permutation)/2,len(permutation)):
                              for nindexb in range(len(permutation)/2,len(permutation)):
                                 if (nindexa >= nindexb):
                                    continue
                                 indexa = permutation[nindexa]
                                 indexb = permutation[nindexb]
                                 if ((indexa.isin(super)) and (indexb.isin(super)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                                 if ((indexa.isin(sub)) and (indexb.isin(sub)) and (indexa.type == "particle") and (indexb.type == "hole")):
                                    rejected = 1
                           if (not rejected):
                              if (parity == 1):
                                 sign = " + "
                              else:
                                 sign = " - "
                              newdbl = child.contraction.tensors[2].fortran90(types,permutation,0,"e")
                              newline = string.join([newdbl,"=",sign,child.contraction.tensors[2].fortran90(types)],"")
                              newcode.statements.insert(newcode.pointer,newline)
                              newcode.pointer = newcode.pointer + 1
                     newcode.pointer = len(newcode.statements)
                  
            # generate loops over target indexes
            if (mode == "nopermutation"):
               for index in child.contraction.tensors[0].indexes:
                  newcode.insertdoloop(index)
               newline = "TMP=0.0d0"
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               pointersave = newcode.pointer
               for npermutation in range(len(child.contraction.factor.permutations)):
                  permutation = child.contraction.factor.permutations[npermutation]
                  newdbl = child.contraction.tensors[0].fortran90(types,permutation,1)
                  newline = string.join([newdbl,"=",newdbl,"+",\
                     "(",repr(child.contraction.factor.coefficients[npermutation]),"d0)*TMP"],"")
                  newcode.statements.insert(newcode.pointer,newline)
                  newcode.pointer = newcode.pointer + 1
               newcode.pointer = pointersave
            elif (mode == "permutation"):
               for index in superglobalzero:
                  newcode.insertdoloop(index)
               if (len(superglobalzero) > 1):
                  newcode.insertif(superglobalzero,1)
               for index in superlocalzero:
                  newcode.insertdoloop(index)
               if (len(superlocalzero) > 1):
                  newcode.insertif(superlocalzero,1)
               for index in subglobalzero:
                  newcode.insertdoloop(index)
               if (len(subglobalzero) > 1):
                  newcode.insertif(subglobalzero,1)
               for index in sublocalzero:
                  newcode.insertdoloop(index)
               if (len(sublocalzero) > 1):
                  newcode.insertif(sublocalzero,1)
               if ((zeroscratch) and (len(child.contraction.tensors) > 2)):
                  newdbl = child.contraction.tensors[0].fortran90(types)
                  newline = string.join([newdbl,"=0.0d0"],"")
                  newcode.statements.insert(newcode.pointer,newline)
                  newcode.pointer = newcode.pointer + 1
                  zeroscratch = 0

            # generate loops over common indexes
            if (mode == "nopermutation"):
               if (child.contraction.summation):
                  for index in child.contraction.summation.indexes:
                     newcode.insertdoloop(index)
               newline = "TMP=TMP+"
               for ntensor in range(len(child.contraction.tensors)):
                  if (ntensor == 1):
                     newdbl = child.contraction.tensors[ntensor].fortran90(types)
                     newline = string.join([newline,newdbl],"")
                  elif (ntensor > 1):
                     newdbl = child.contraction.tensors[ntensor].fortran90(types)
                     newline = string.join([newline,"*",newdbl],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newcode.pointer = len(newcode.statements)
            elif (mode == "permutation"):
               factor = 1
               for index in supercommon:
                  newcode.insertdoloop(index)
# DEEXCITATION EXTENSION FROM HERE ...
#              if (len(supercommon) > 1):
#                 newcode.insertif(supercommon,1)
#                 factor = factor * factorial(len(supercommon))
#              if (len(subcommon) > 1):
#                 newcode.insertif(subcommon,1)
#                 factor = factor * factorial(len(subcommon))
# In the following, a multiplicative factor that compensates
# the triangular summation index is determined.  In so doing,
# we must distinguish hole and particle types when deexcitation
# operator is present.
               supercommonhole = []
               supercommonparticle = []
               for index in supercommon:
                  if (index.ishole()):
                     supercommonhole.append(index)
                  elif (index.isparticle()):
                     supercommonparticle.append(index)
                  else:
                     raise RuntimeError, "a general summation index appeared"
               if (len(supercommonhole) > 1):
                  newcode.insertif(supercommonhole,1)
                  factor = factor * factorial(len(supercommonhole))
               if (len(supercommonparticle) > 1):
                  newcode.insertif(supercommonparticle,1)
                  factor = factor * factorial(len(supercommonparticle))
               for index in subcommon:
                  newcode.insertdoloop(index)
               subcommonhole = []
               subcommonparticle = []
               for index in subcommon:
                  if (index.ishole()):
                     subcommonhole.append(index)
                  elif (index.isparticle()):
                     subcommonparticle.append(index)
                  else:
                     raise RuntimeError, "a general summation index appeared"
               if (len(subcommonhole) > 1):
                  newcode.insertif(subcommonhole,1)
                  factor = factor * factorial(len(subcommonhole))
               if (len(subcommonparticle) > 1):
                  newcode.insertif(subcommonparticle,1)
                  factor = factor * factorial(len(subcommonparticle))
# In the following, permutation of local (internal) indices will be 
# performed.  This occurs when there is a deexcitation and when and 
# only when there is no permutation operation of target indices.
# When these two exclusive permutation operations coexist, an error results.
               superlocalparticleone  = []
               superlocalparticletwo  = []
               for index in superlocalone:
                  if (index.isparticle()):
                     superlocalparticleone.append(index)
               for index in superlocaltwo:
                  if (index.isparticle()):
                     superlocalparticletwo.append(index)
               superlocalparticlezero = sortindexes( superlocalparticleone + superlocalparticletwo )
               sublocalholeone  = []
               sublocalholetwo  = []
               for index in sublocalone:
                  if (index.ishole()):
                     sublocalholeone.append(index)
               for index in sublocaltwo:
                  if (index.ishole()):
                     sublocalholetwo.append(index)
               sublocalholezero = sortindexes( sublocalholeone + sublocalholetwo )
               deexcitationfactor = Factor([1.0],[[]])
               if ((superlocalparticleone and superlocalparticletwo) or \
                   (sublocalholeone and sublocalholetwo)):
                  deexcitationfactor = permutationsoffoursets(deexcitationfactor,superlocalparticleone, superlocalparticletwo, \
                                                                                 sublocalholeone,       sublocalholetwo, \
                                                                                 superlocalparticlezero,sublocalholezero)
                  deexcitationfactor = deexcitationfactor.normalize()
# ... TO HERE
               newdbl = child.contraction.tensors[0].fortran90(types)
               if ((zeroscratch) and (len(child.contraction.tensors) == 2)):
                  newline = string.join([newdbl,"="],"")
                  zeroscratch = 0
               else:
                  newline = string.join([newdbl,"=",newdbl,"+"],"")
# DEEXCITATION EXTENSION FROM HERE ...
               if ((len(deexcitationfactor.permutations) > 1) and (len(child.contraction.factor.permutations) > 1)):
                  raise RuntimeError, "A logical error in code generator regarding deexcitation operator"
               if (len(deexcitationfactor.permutations) > 1):
                  deexcitationfactor.multiply(child.contraction.factor.coefficients[0])
                  currentfactor = deexcitationfactor.duplicate()
               else:
                  currentfactor = child.contraction.factor.duplicate()
               for npermutation in range(len(currentfactor.permutations)):
                  permutation = currentfactor.permutations[npermutation]
                  if (npermutation > 0):
                     newline = string.join([newline,"+"],"")
                  for ntensor in range(len(child.contraction.tensors)):
                     newfactor = float(factor) * currentfactor.coefficients[npermutation]
# ... TO HERE
                     if (expanded[ntensor]):
                        suffix = "e"
                     else:
                        suffix = ""
                     if (ntensor == 1):
                        newdbl = child.contraction.tensors[ntensor].fortran90(types,permutation,0,suffix)
                        newline = string.join([newline,"(",repr(newfactor),"d0)*",newdbl],"")
                     elif (ntensor > 1):
                        newdbl = child.contraction.tensors[ntensor].fortran90(types,permutation,0,suffix)
                        newline = string.join([newline,"*",newdbl],"")
               newcode.statements.insert(newcode.pointer,newline)
               newcode.pointer = newcode.pointer + 1
               newcode.pointer = len(newcode.statements)

      return newcode

class Code:
 
   def __init__(self,language,name):
      """Creates an empty code of program"""
      self.language = language
      self.name = name
      self.statements = []
      self.pointer = 0
      self.headers = []
      self.arguments = []
      self.integers = []
      self.integerarrays = []
      self.integerallocatables = []
      self.doubles = []
      self.doublearrays = []
      self.doubleallocatables = []
      self.doublecomplexs = []
      self.doublecomplexarrays = []
      self.doublecomplexallocatables = []
      self.logicals = []
      self.logicalarrays = []
      self.characters = []
      self.externals = []

      # Comment lines
      if (self.language == "Fortran77"):
         self.comment = "C     "
         self.indent  = "      "
         self.nlang = 0
      elif (self.language == "Fortran90"):
         self.comment = "! "
         self.indent  = ""
         self.nlang = 1
      elif (self.language == "Python"):
         self.comment = "# "
         self.indent  = ""
         self.nlang = 2
      else:
         return "Unknown language"

      # Standard headers
      newline = "!$Id$"
      self.headers.append(newline)
      newline = "!This is a " + self.language + " program generated by Tensor Contraction Engine v.1.0"
      self.headers.append(newline)
      newline = "!Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)"
      self.headers.append(newline)

   def __str__(self):
      """Prints code"""
      print ""
      for line in self.wrap():
         print line
      return ""

   def isnested(self):
      """Returns true if self.statements contains a nested code object"""
      for member in self.statements:
         if (isinstance(member,Code)):
            return 1
      return 0

   def expand(self):
      """Expands a code object with a nested statement list into a non-nested code"""
      result = Code(self.language,self.name)
      result.pointer = self.pointer
      for n in self.headers:
         result.add("headers",n)
      for n in self.arguments:
         result.add("arguments",n)
      for n in self.integers:
         result.add("integers",n)
      for n in self.integerarrays:
         result.add("integerarrays",n)
      for n in self.integerallocatables:
         result.add("integerallocatables",n)
      for n in self.doubles:
         result.add("doubles",n)
      for n in self.doublearrays:
         result.add("doublearrays",n)
      for n in self.doubleallocatables:
         result.add("doubleallocatables",n)
      for n in self.doublecomplexs:
         result.add("doublecomplexs",n)
      for n in self.doublecomplexarrays:
         result.add("doublecomplexarrays",n)
      for n in self.doublecomplexallocatables:
         result.add("doublecomplexallocatables",n)
      for n in self.logicals:
         result.add("logicals",n)
      for n in self.logicalarrays:
         result.add("logicalarrays",n)
      for n in self.characters:
         result.add("characters",n)
      for n in self.externals:
         result.add("externals",n)
      for member in self.statements:
         if (isinstance(member,Code)):
            result.join(member.expand())
         else:
            result.statements.append(member)
      return result

   def setamark(self,number):
      """Inserts a special symbol with an identifier number"""
      self.deleteamark(number)
      statement = self.statements[self.pointer]
      statement = "#marker"+repr(number)+"#"+statement
      self.statements[self.pointer]=statement

   def getamark(self,number):
      """Returns the position of an input marker"""
      for nstatement in range(len(self.statements)):
         statement = self.statements[nstatement]
         if (statement[0:7] == "#marker"):
            if (number == int(string.split(statement[7:],"#")[0])):
               position = string.index(statement[7:],"#")
               self.statements[nstatement] = statement[position+8:]
               return nstatement

      raise ValueError, "Maker not found"

   def deleteamark(self,number):
      """Delete a marker"""
      for nstatement in range(len(self.statements)):
         statement = self.statements[nstatement]
         if (statement[0:7] == "#marker"):
            if (number == int(string.split(statement[7:],"#")[0])):
               position = string.index(statement[7:],"#")
               self.statements[nstatement] = statement[position+8:]

   def deleteallmarks(self):
      """Deletes all existing marks"""
      for nstatement in range(len(self.statements)):
         statement = self.statements[nstatement]
         if (statement[0:7] == "#marker"):
            position = string.index(statement[7:],"#")
            self.statements[nstatement] = statement[position+8:]

   def show(self): 
      """Returns an output of the contents"""

      self.deleteallmarks()

      if (self.isnested()):
         return "This code object is nested; first use expand()"

      show = []

      if (self.language == "Fortran77"):

         # add the headers and declarations
         pointer = 0
         subroutine = string.join([self.indent,"SUBROUTINE ",self.name],"")
         if (self.arguments):
            subroutine = string.join([subroutine,"("],"")
            for n in range(len(self.arguments)):
               argument = self.arguments[n]
               if (n != 0):
                  argument = string.join([",",argument],"")
               subroutine = string.join([subroutine,argument],"")
            subroutine = string.join([subroutine,")"],"")
         show.insert(pointer,subroutine)
         pointer = pointer + 1
         for n in self.headers:
            if (n[0] == "#"):
               show.insert(pointer,n)
            elif (n[0] == "!"):
               show.insert(pointer,string.join([self.comment,n[1:]],""))
            else:
               show.insert(pointer,string.join([self.indent,n],""))
            pointer = pointer + 1
         for n in self.integers:
            show.insert(pointer,string.join([self.indent,"INTEGER ",n],""))
            pointer = pointer + 1
         for n in self.integerarrays:
            show.insert(pointer,string.join([self.indent,"INTEGER ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.doubles:
            show.insert(pointer,string.join([self.indent,"DOUBLE PRECISION ",n],""))
            pointer = pointer + 1
         for n in self.doublearrays:
            show.insert(pointer,string.join([self.indent,"DOUBLE PRECISION ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.doublecomplexs:
            show.insert(pointer,string.join([self.indent,"DOUBLE COMPLEX ",n],""))
            pointer = pointer + 1
         for n in self.doublecomplexarrays:
            show.insert(pointer,string.join([self.indent,"DOUBLE COMPLEX ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.logicals:
            show.insert(pointer,string.join([self.indent,"LOGICAL ",n],""))
            pointer = pointer + 1
         for n in self.logicalarrays:
            show.insert(pointer,string.join([self.indent,"LOGICAL ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.characters:
            if (n in self.arguments):
               show.insert(pointer,string.join([self.indent,"CHARACTER*(*) ",n],""))
               pointer = pointer + 1
            else:
               show.insert(pointer,string.join([self.indent,"CHARACTER*255 ",n],""))
               pointer = pointer + 1
         for n in self.externals:
            show.insert(pointer,string.join([self.indent,"EXTERNAL ",n],""))
            pointer = pointer + 1

         # add the statements
         for n in self.statements:
            show.insert(pointer,string.join([self.indent,n],""))
            pointer = pointer + 1

      elif (self.language == "Fortran90"):

         # add the headers and declarations
         pointer = 0
         subroutine = string.join([self.indent,"SUBROUTINE ",self.name],"")
         if (self.arguments):
            subroutine = string.join([subroutine,"("],"")
            for n in range(len(self.arguments)):
               argument = self.arguments[n]
               if (n != 0):
                  argument = string.join([",",argument],"")
               subroutine = string.join([subroutine,argument],"")
            subroutine = string.join([subroutine,")"],"")
         show.insert(pointer,subroutine)
         pointer = pointer + 1
         for n in self.headers:
            if ((n[0] == "#") or (n[0] == "!")):
               show.insert(pointer,n)
            else:
               show.insert(pointer,string.join([self.indent,n],""))
            pointer = pointer + 1
         for n in self.integers:
            show.insert(pointer,string.join([self.indent,"INTEGER :: ",n],""))
            pointer = pointer + 1
         for n in self.integerarrays:
            show.insert(pointer,string.join([self.indent,"INTEGER :: ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.integerallocatables:
            show.insert(pointer,string.join([self.indent,"INTEGER, ALLOCATABLE :: ",n,"(:)"],""))
            pointer = pointer + 1
         for n in self.doubles:
#           show.insert(pointer,string.join([self.indent,"DOUBLE PRECISION :: ",n],""))
            show.insert(pointer,string.join([self.indent,"REAL*8 :: ",n],""))
            pointer = pointer + 1
         for n in self.doublearrays:
#           show.insert(pointer,string.join([self.indent,"DOUBLE PRECISION :: ",n,"(*)"],""))
            show.insert(pointer,string.join([self.indent,"REAL*8 :: ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.doubleallocatables:
#           show.insert(pointer,string.join([self.indent,"DOUBLE PRECISION, ALLOCATABLE :: ",n,"(:)"],""))
            show.insert(pointer,string.join([self.indent,"REAL*8, ALLOCATABLE :: ",n,"(:)"],""))
            pointer = pointer + 1
         for n in self.doublecomplexs:
#           show.insert(pointer,string.join([self.indent,"DOUBLE COMPLEX :: ",n],""))
            show.insert(pointer,string.join([self.indent,"COMPLEX*16 :: ",n],""))
            pointer = pointer + 1
         for n in self.doublecomplexarrays:
#           show.insert(pointer,string.join([self.indent,"DOUBLE COMPLEX :: ",n,"(*)"],""))
            show.insert(pointer,string.join([self.indent,"COMPLEX*16 :: ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.doublecomplexallocatables:
#           show.insert(pointer,string.join([self.indent,"DOUBLE COMPLEX, ALLOCATABLE :: ",n,"(:)"],""))
            show.insert(pointer,string.join([self.indent,"COMPLEX*16, ALLOCATABLE :: ",n,"(:)"],""))
            pointer = pointer + 1
         for n in self.logicals:
            show.insert(pointer,string.join([self.indent,"LOGICAL :: ",n],""))
            pointer = pointer + 1
         for n in self.logicalarrays:
            show.insert(pointer,string.join([self.indent,"LOGICAL :: ",n,"(*)"],""))
            pointer = pointer + 1
         for n in self.characters:
            show.insert(pointer,string.join([self.indent,"CHARACTER(LEN=255) :: ",n],""))
            pointer = pointer + 1
         for n in self.externals:
            show.insert(pointer,string.join([self.indent,"EXTERNAL :: ",n],""))
            pointer = pointer + 1

         # add the statements
         for n in self.statements:
            show.insert(pointer,string.join([self.indent,n],""))
            pointer = pointer + 1

      return show

   def wrap(self):
      """Wraps around long statements; calls show()"""

      show = self.show()
      
      if (self.language == "Fortran77"):
         show72 = []
         for n in show:
            if ((n[0] == "C") or (n[0] == "c")):
               done = 1
            else:
               done = 0
            while (not done):
               if (len(n) > 72):
                  show72.append(n[0:72])
                  n = string.join(["     &",n[72:]],"")
               else:
                  done = 1
            show72.append(n)
         return show72
      elif (self.language == "Fortran90"):
         show132 = []
         for n in show:
            if (n[0] == "!"):
               done = 1
            else:
               done = 0
            while (not done):
               if (len(n) > 132):
                  show132.append(string.join([n[0:131],"&"],""))
                  n = "&"+n[131:]
               else:
                  done = 1
            show132.append(n)
         return show132
      else:
         return show

   def join(self,another):
      """Join two code objects together"""
      if (self.language != another.language):
         return "Cannot join two codes"
      for n in another.headers:
         self.add("headers",n)
      for n in another.arguments:
         self.add("arguments",n)
      for n in another.integers:
         self.add("integers",n)
      for n in another.integerarrays:
         self.add("integerarrays",n)
      for n in another.integerallocatables:
         self.add("integerallocatables",n)
      for n in another.doubles:
         self.add("doubles",n)
      for n in another.doublearrays:
         self.add("doublearrays",n)
      for n in another.doubleallocatables:
         self.add("doubleallocatables",n)
      for n in another.doublecomplexs:
         self.add("doublecomplexs",n)
      for n in another.doublecomplexarrays:
         self.add("doublecomplexarrays",n)
      for n in another.doublecomplexallocatables:
         self.add("doublecomplexallocatables",n)
      for n in another.logicals:
         self.add("logicals",n)
      for n in another.logicalarrays:
         self.add("logicalarrays",n)
      for n in another.characters:
         self.add("characters",n)
      for n in another.statements:
         self.statements.append(n)
      for n in another.externals:
         self.add("externals",n)
      return self

   def add(self,towhat,what):
      """Add a new integer/double/logical etc to an existing list; checks redundancy"""
      if (towhat == "integers"):
         redundant = 0
         for n in self.integers:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.integers.append(what)
      elif (towhat == "integerarrays"):
         redundant = 0
         for n in self.integerarrays:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.integerarrays.append(what)
      elif (towhat == "integerallocatables"):
         redundant = 0
         for n in self.integerallocatables:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.integerallocatables.append(what)
      elif (towhat == "doubles"):
         redundant = 0
         for n in self.doubles:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.doubles.append(what)
      elif (towhat == "doublearrays"):
         redundant = 0
         for n in self.doublearrays:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.doublearrays.append(what)
      elif (towhat == "doubleallocatables"):
         redundant = 0
         for n in self.doubleallocatables:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.doubleallocatables.append(what)
      elif (towhat == "doublecomplexs"):
         redundant = 0
         for n in self.doublecomplexs:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.doublecomplexs.append(what)
      elif (towhat == "doublecomplexarrays"):
         redundant = 0
         for n in self.doublecomplexarrays:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.doublecomplexarrays.append(what)
      elif (towhat == "doublecomplexallocatables"):
         redundant = 0
         for n in self.doublecomplexallocatables:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.doublecomplexallocatables.append(what)
      elif (towhat == "logicals"):
         redundant = 0
         for n in self.logicals:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.logicals.append(what)
      elif (towhat == "logicalarrays"):
         redundant = 0
         for n in self.logicalarrays:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.logicalarrays.append(what)
      elif (towhat == "characters"):
         redundant = 0
         for n in self.characters:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.characters.append(what)
      elif (towhat == "externals"):
         redundant = 0
         for n in self.externals:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.externals.append(what)
      elif (towhat == "arguments"):
         redundant = 0
         for n in self.arguments:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.arguments.append(what)
      elif (towhat == "headers"):
         redundant = 0
         for n in self.headers:
            if (n == what):
               redundant = 1
         if (not redundant):
            self.headers.append(what)
      else:
         raise ValueError, "Unknown variable type"

   def insertdoloop(self,index):
      """Inserts a DO-ENDDO pair into a code"""
      if (index.type == 'hole'):
         newint = index.show()
         self.add("integers",newint)
         newline = string.join([self.indent,"DO ",newint,"=1,nocc"],"")
         self.statements.insert(self.pointer,newline)
         self.pointer = self.pointer + 1
         newline = "END DO"
         self.statements.insert(self.pointer,newline)
      elif (index.type == 'particle'):
         newint = index.show()
         self.add("integers",newint)
         newline = string.join([self.indent,"DO ",newint,"=nocc+1,N"],"")
         self.statements.insert(self.pointer,newline)
         self.pointer = self.pointer + 1
         newline = "END DO"
         self.statements.insert(self.pointer,newline)
      elif (index.type == 'general'):
         newint = index.show()
         self.add("integers",newint)
         newline = string.join([self.indent,"DO ",newint,"=1,N"],"")
         self.statements.insert(self.pointer,newline)
         self.pointer = self.pointer + 1
         newline = "END DO"
         self.statements.insert(self.pointer,newline)

   def insertif(self,list,holeisalwayslessthanparticle=0):
      """Inserts an IF sentence for skipping permutation redundant block"""
      for nindex in range(len(list)-1):
         indexa = list[nindex]
         indexb = list[nindex+1]
         if (holeisalwayslessthanparticle):
            if ((indexa.type == "hole") and (indexb.type == "particle")):
# DEEXCITATION EXTENSION FROM HERE ...
               continue
#              return
# ... TO HERE
            elif ((indexa.type == "particle") and (indexb.type == "hole")):
               raise ValueError, "A particle, hole sequence in a tensor"
         newline = string.join(["IF (",indexa.show(),">=",indexb.show(),") CYCLE"],"")
         self.statements.insert(self.pointer,newline)
         self.pointer = self.pointer + 1

   def inserttileddoloops(self,list):
      """Inserts a nested restricted DO-ENDDO pair"""
      for nindex in range(len(list)):
         index = list[nindex]
         newint = string.join([index.show(),"b"],"")
         self.add("integers",newint)
         if (nindex == 0):
            if (index.type == 'hole'):
               newline = string.join(["DO ",newint," = 1,noab"],"")
            if (index.type == 'particle'):
               newline = string.join(["DO ",newint," = noab+1,noab+nvab"],"")
            self.statements.insert(self.pointer,newline)
            self.pointer = self.pointer + 1
            newline = "END DO"
            self.statements.insert(self.pointer,newline)
            previousint = newint
            previoustype = index.type
         else:
            if (index.type == 'hole'):
               if (previoustype == 'hole'):
                  newline = string.join(["DO ",newint," = ",previousint,",noab"],"")
               else:
                  raise ValueError, "non-canonical expression found"
            if (index.type == 'particle'):
               if (previoustype == 'hole'):
                  newline = string.join(["DO ",newint," = noab+1,noab+nvab"],"")
               else:
                  newline = string.join(["DO ",newint," = ",previousint,",noab+nvab"],"")
            self.statements.insert(self.pointer,newline)
            self.pointer = self.pointer + 1
            newline = "END DO"
            self.statements.insert(self.pointer,newline)
            previousint = newint
            previoustype = index.type
 
   def inserttiledifrestricted(self,indexes):
      """Inserts an IF-ENDIF pair for screening all beta intermediates"""

      if (not indexes):
         return

      # spin symmetry
      newline = "IF ((.not.restricted).or.("
      conjugation = ""
      for index in indexes:
         newint = string.join([index.show(),"b"],"")
         if (self.language == "Fortran77"):
            newline = string.join([newline,conjugation,"int_mb(k_spin+",newint,"-1)"],"")
         elif (self.language == "Fortran90"):
            newline = string.join([newline,conjugation,"spin(",newint,")"],"")
         conjugation = "+"
      if (self.language == "Fortran77"):
         newline = string.join([newline,".ne.",repr(2*len(indexes)),")) THEN"],"")
      elif (self.language == "Fortran90"):
         newline = string.join([newline," /= ",repr(2*len(indexes)),")) THEN"],"")
      self.statements.insert(self.pointer,newline)
      self.pointer = self.pointer + 1
      newline = "END IF"
      self.statements.insert(self.pointer,newline)
 
   def inserttiledifsymmetry(self,super,sub,irrep,relativistic=0):
      """Inserts an IF-ENDIF pair for screening spin/spatial symmetry"""

# SYMMETRY EXTENSION FROM HERE ...
#     if ((not super) and (not sub)):
#        return
# ... TO HERE
      if ((super and (not sub)) or ((not super) and sub)):
         raise ValueError, "asymmetric IF encountered"

      # spin symmetry
# SYMMETRY EXTENSION FROM HERE ...
      if ((super or sub) and (not relativistic)):
#     if (not relativistic):
# ... TO HERE
         newline = "IF ("
         conjugation = ""
         for index in super:
            newint = string.join([index.show(),"b"],"")
            if (self.language == "Fortran77"):
               newline = string.join([newline,conjugation,"int_mb(k_spin+",newint,"-1)"],"")
            elif (self.language == "Fortran90"):
               newline = string.join([newline,conjugation,"spin(",newint,")"],"")
            conjugation = "+"
         if (self.language == "Fortran77"):
            conjugation = " .eq. "
         elif (self.language == "Fortran90"):
            conjugation = " == "
         for index in sub:
            newint = string.join([index.show(),"b"],"")
            if (self.language == "Fortran77"):
               newline = string.join([newline,conjugation,"int_mb(k_spin+",newint,"-1)"],"")
            elif (self.language == "Fortran90"):
               newline = string.join([newline,conjugation,"spin(",newint,")"],"")
            conjugation = "+"
         newline = string.join([newline,") THEN"],"")
         self.statements.insert(self.pointer,newline)
         self.pointer = self.pointer + 1
         newline = "END IF"
         self.statements.insert(self.pointer,newline)

      # spatial symmetry
# SYMMETRY EXTENSION FROM HERE ...
      if ((not super) and (not sub)):
         newline = "IF (0"
      else:
# ... TO HERE
         all = super + sub
         newline = "IF ("
         conjugation = ""
         for nindex in range(len(all)-1):
            index = all[nindex]
            newint = string.join([index.show(),"b"],"")
            if (self.language == "Fortran77"):
               newline = string.join([newline,conjugation,"ieor(int_mb(k_sym+",newint,"-1)"],"")
            elif (self.language == "Fortran90"):
               newline = string.join([newline,conjugation,"ieor(sym(",newint,")"],"")
            conjugation = ","
         index = all[len(all)-1]
         newint = string.join([index.show(),"b"],"")
         if (self.language == "Fortran77"):
            newline = string.join([newline,",int_mb(k_sym+",newint,"-1)"],"")
         elif (self.language == "Fortran90"):
            newline = string.join([newline,",sym(",newint,")"],"")
         for nindex in range(len(all)-1):
            newline = string.join([newline,")"],"")
      if (self.language == "Fortran77"):
         newline = string.join([newline," .eq. "],"")
      elif (self.language == "Fortran90"):
         newline = string.join([newline," == "],"")
      if (len(irrep) == 1):
         newline = string.join([newline,"irrep_",irrep[0],") THEN"],"")
      else:
         conjugation = ""
         for neach in range(len(irrep)-1):
            each = irrep[neach]
            newline = string.join([newline,conjugation,"ieor(irrep_",each],"")
            conjugation = ","
         newline = string.join([newline,",irrep_",irrep[len(irrep)-1]],"")
         for neach in range(len(irrep)-1):
            newline = string.join([newline,")"],"")
         newline = string.join([newline,") THEN"],"")
      self.statements.insert(self.pointer,newline)
      self.pointer = self.pointer + 1
      newline = "END IF"
      self.statements.insert(self.pointer,newline)

   def inserttiledifpermutation(self,list,holeisalwayslessthanparticle=0):
      """Inserts an IF sentence for skipping permutation redundant block"""
      for nindex in range(len(list)-1):
         indexa = list[nindex]
         indexb = list[nindex+1]
         if (holeisalwayslessthanparticle):
            if ((indexa.type == "hole") and (indexb.type == "particle")):
# DEEXCITATION EXTENSION FROM HERE ...
               continue
#              return
# ... TO HERE
            elif ((indexa.type == "particle") and (indexb.type == "hole")):
               raise ValueError, "A particle, hole sequence in a tensor"
         if (self.language == "Fortran77"):
            newline = string.join(["IF (",indexa.show(),"b .le. ",indexb.show(),"b) THEN"],"")
         elif (self.language == "Fortran90"):
            newline = string.join(["IF (",indexa.show(),"b <= ",indexb.show(),"b) THEN"],"")
         self.statements.insert(self.pointer,newline)
         self.pointer = self.pointer + 1
         newline = "END IF"
         self.statements.insert(self.pointer,newline)

   def writetofile(self,filename):
      """Writes a list to a given file"""

      if (self.isnested()):
         return "This code object is nested; first use expand()"

      if (self.language == "Fortran77"):
         file = open(filename+".F","w")
      elif (self.language == "Fortran90"):
         file = open(filename+".f90","w")
      for n in self.wrap():
         file.write(n)
         file.write("\n")

   def sortarguments(self):
      """Sorts the list of arguments in an ascending order"""

      done = 0
      while (not done):
         done = 1
         for n in range(len(self.arguments)):
            for m in range(len(self.arguments)):
               if (n >= m):
                  continue
               if (self.arguments[n] > self.arguments[m]):
                  swap = self.arguments[n]
                  self.arguments[n] = self.arguments[m]
                  self.arguments[m] = swap
                  done = 0

   def reverse(self):
      """Reverses the execution seqeuence of the whole code"""

      newstatements = []
      for statement in self.statements:
         newstatements.insert(0,statement)
      self.statements = copy.deepcopy(newstatements)

   def removeredundantio(self):
      """Removes unneccesary GATODRA-DRATOGA pairs"""

      newstatements = []
      temporary = []
      arguments = []
      type = []
      erased = []
      for statement in self.statements:
         if (string.find(statement,"DRATOGA") != -1):
            temporary.append(statement)
            endline = len(statement)-1
            beginline = string.find(statement,"(")+1
            arguments.append(statement[beginline:endline])
            type.append(1)
            erased.append(0)
         elif (string.find(statement,"GATODRA") != -1):
            temporary.append(statement)
            endline = len(statement)-1
            beginline = string.find(statement,"(")+1
            arguments.append(statement[beginline:endline])
            type.append(-1)
            erased.append(0)
         else:
            if (temporary):
               for i in range(len(temporary)):
                  if (type[i] == 1):
                     for j in range(len(temporary)):
                        if (type[j] == -1) and (arguments[i] == arguments[j]) and (not erased[j]):
                           erased[i] = 1
                           erased[j] = 1
               for i in range(len(temporary)):
                  if (not erased[i]):
                     newstatements.append(temporary[i])
               temporary = []
               arguments = []
               type = []
               erased = []
            newstatements.append(statement)
      self.statements = copy.deepcopy(newstatements)

class ListofCodes:
 
   def __init__(self):
      """Creates an empty code of program"""
      self.list = []

   def add(self,code):
      """Adds a code"""
      self.list.append(code)

   def join(self,another):
      """Joins two lists of codes"""
      for code in another.list:
         self.list.append(code)

   def show(self):
      """Calls wrap() of each Code object"""
      show = []
      for code in self.list:
         show = show + code.wrap()
      return show

   def __str__(self):
      """Prints code"""
      print ""
      for line in self.show():
         print line
      return ""

   def writetofile(self,filename):
      """Writes a list to a given file"""

      for code in self.list:
         if (code.isnested()):
            return "This code object is nested; first use expand()"

      if (self.list[0].language == "Fortran77"):
         file = open(filename+".F","w")
      elif (self.list[0].language == "Fortran90"):
         file = open(filename+".F90","w")
      for code in self.list:
         for n in code.wrap():
            file.write(n)
            file.write("\n")
