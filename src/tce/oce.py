# Operator Contraction Engine v.1.0
# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)
# $Id$

import string

import copy

import sys

def readfromfile(filename):
   """Converts the content of a file to a ListOperatorSequences object"""

   result = ListOperatorSequences()
   file = open(filename,"r")
   alwaystrue = 1
   while (alwaystrue):
      line = file.readline()
      if (line == ""):
         file.close()
         return result
      else:
         line = line[0:len(line)-1]
         result.add(stringtooperatorsequence(line))

def stringtooperatorsequence(expression):
   """Converts a string to an operatorsequence object"""
   # Syntax of the string is rather loosely defined as:
   # (1) Numerical factor (with no permutation allowed) (optional), summation (optional), amplitudes (optional), normal ordered sequence,
   # (2) Numerical factor can be an arithmatic expression such as (1.0/4.0),
   # (3) Summation starts with either "SUM" or "sum" followed by a parenthesis of indexes,
   # (4) Indexes can be either in one-letter notation (a-h, A-H for virtuals, i-o, I-O for occupieds, p-z, P-Z for either, case matters) 
   #     or in OCE notation (p1,p2 for virtuals, h3,h4 for occupieds, g5,g6 for either, no overlap in numbering)
   # (5) Amplitudes start with "t" or any name followed by a dagger ("+") indicating complex conjugate (optional) and a parenthesis of indexes,
   # (6) Normal ordered operator sequence must exist even when it is empty "{}".
   # (7) An example is: (1.0/16.0) Sum (p q r s c d k l) v(p q r s) t(c d k l) {i+ j+ b a}{p+ q+ s r}{c+ d+ l k}
   
   sequences = expression[expression.index("{"):]
   expression = expression[0:expression.index("{")]
   operatorlist = []
   # first we decipher normal ordered operator sequence and define operators with/without daggers
   newsequences = []
   while (string.find(sequences,"{") != -1):
      sequences = sequences[0:string.find(sequences,"{")] + sequences[string.find(sequences,"{")+1:]
   if (sequences[len(sequences)-1] != "}"):
      raise RuntimeError, "Syntax error: the string must end with a normal ordered operator sequence"
   sequences = sequences[0:len(sequences)-1]
   sequences = string.split(sequences,"}")
   for sequence in sequences:
      newsequence = []
      sequence = string.split(sequence)
      for index in sequence:
         if (index[len(index)-1] == "+"):
            dagger = "creation"
            index = index[0:len(index)-1]
         else:
            dagger = "annihilation"
         if (len(index) == 1):
            # one letter notation (vir: a-h & A-H; occ: i-o & I-O; gen: p-z & P-Z)
            if (((index >= "a") and (index <= "h")) or ((index >= "A") and (index <= "H"))):
               newsequence.append(Operator("particle",dagger,string.ascii_letters.index(index)+1))
            elif (((index >= "i") and (index <= "o")) or ((index >= "I") and (index <= "O"))):
               newsequence.append(Operator("hole",dagger,string.ascii_letters.index(index)+1))
            elif (((index >= "p") and (index <= "z")) or ((index >= "P") and (index <= "Z"))):
               newsequence.append(Operator("general",dagger,string.ascii_letters.index(index)+1))
         else:
            if (index[0] == "p"):
               newsequence.append(Operator("particle",dagger,int(index[1:])))
            elif (index[0] == "h"):
               newsequence.append(Operator("hole",dagger,int(index[1:])))
            elif (index[0] == "g"):
               newsequence.append(Operator("general",dagger,int(index[1:])))
            else:
               print "Syntax error: an operator not recognized"
               stop
      operatorlist = operatorlist + newsequence
      newsequences.append(newsequence)

   breakdown = string.split(expression)

   # get a numerical factor if any
   numericalfactor = ""
   for element in breakdown:
      if (((element[0] >= 'a') and (element[0] <= 'z')) or \
          ((element[0] >= 'A') and (element[0] <= 'Z'))):
         break
      numericalfactor = string.join([numericalfactor, element])
   if (numericalfactor == ""):
      numericalfactor = Factor([1.0],[[]])
   else:
      numericalfactor = eval(numericalfactor)
      numericalfactor = Factor([numericalfactor],[[]])

   # get a summation if any
   summationindexes = ""
   remainder = ""
   join = 0
   for element in breakdown:
      if ((element[0:3] == "SUM") or (element[0:3] == "sum") or (element[0:3] == "Sum")):
         join = 1
      if (join == 1):
         summationindexes = string.join([summationindexes,element])
      elif (join == 2):
         remainder = string.join([remainder,element])
      if ((join == 1) and (")" in element)):
         join = 2
   index = string.find(summationindexes,"sum")
   if (index != -1):
      summationindexes = summationindexes[0:index] + summationindexes[index+3:]
   index = string.find(summationindexes,"SUM")
   if (index != -1):
      summationindexes = summationindexes[0:index] + summationindexes[index+3:]
   index = string.find(summationindexes,"Sum")
   if (index != -1):
      summationindexes = summationindexes[0:index] + summationindexes[index+3:]
   index = string.find(summationindexes,"(")
   if (index != -1):
      summationindexes = summationindexes[0:index] + summationindexes[index+1:]
   index = string.find(summationindexes,")")
   if (index != -1):
      summationindexes = summationindexes[0:index] + summationindexes[index+1:]
   summationindexes = string.split(summationindexes)
   summation = Summation([])
   for index in summationindexes:
      if (len(index) == 1):
         # one letter notation (vir: a-h & A-H; occ: i-o & I-O; gen: p-z & P-Z)
         if (((index >= "a") and (index <= "h")) or ((index >= "A") and (index <= "H"))):
            for indexinthelist in operatorlist:
               if ((indexinthelist.type == "particle") and (indexinthelist.index == string.ascii_letters.index(index)+1)):
                  summation.indexes.append(indexinthelist)
                  break
         elif (((index >= "i") and (index <= "o")) or ((index >= "I") and (index <= "O"))):
            for indexinthelist in operatorlist:
               if ((indexinthelist.type == "hole") and (indexinthelist.index == string.ascii_letters.index(index)+1)):
                  summation.indexes.append(indexinthelist)
                  break
         elif (((index >= "p") and (index <= "z")) or ((index >= "P") and (index <= "Z"))):
            for indexinthelist in operatorlist:
               if ((indexinthelist.type == "general") and (indexinthelist.index == string.ascii_letters.index(index)+1)):
                  summation.indexes.append(indexinthelist)
                  break
      else:
         if (index[0] == "p"):
            for indexinthelist in operatorlist:
               if ((indexinthelist.type == "particle") and (indexinthelist.index == int(index[1:]))):
                  summation.indexes.append(indexinthelist)
                  break
         elif (index[0] == "h"):
            for indexinthelist in operatorlist:
               if ((indexinthelist.type == "hole") and (indexinthelist.index == int(index[1:]))):
                  summation.indexes.append(indexinthelist)
                  break
         elif (index[0] == "g"):
            for indexinthelist in operatorlist:
               if ((indexinthelist.type == "general") and (indexinthelist.index == int(index[1:]))):
                  summation.indexes.append(indexinthelist)
                  break
         else:
            print "syntax error"
            stop

   # get amplitudes
   remainder = string.split(remainder,")")
   tobeamplitudes = remainder[0:len(remainder)-1]
   remainder = remainder[len(remainder)-1]
   amplitudes = []
   for tobeamplitude in tobeamplitudes:
      newamplitude = Amplitude()
      tobeamplitude = string.split(tobeamplitude,"(")
      type = tobeamplitude[0]
      conjugate = 1
      lastdaggerposition = len(type)
      for i in range(len(type)-1,-1,-1):
         if (type[i] == "+"):
            conjugate = - conjugate
            lastdaggerposition = i
      if (conjugate == -1):
         newamplitude.conjugate = 1
      else:
         newamplitude.conjugate = 0
      newamplitude.type = string.strip(type[0:lastdaggerposition])
      index = 0
      for amplitude in amplitudes:
         if ((amplitude.type == newamplitude.type) and (amplitude.index > index)):
            index = amplitude.index
      newamplitude.index = index + 1
      tobeamplitudeindexes = string.split(tobeamplitude[1])
      newamplitude.indexes = []
      for index in tobeamplitudeindexes:
         if (len(index) == 1):
            # one letter notation (vir: a-h & A-H; occ: i-o & I-O; gen: p-z & P-Z)
            if (((index >= "a") and (index <= "h")) or ((index >= "A") and (index <= "H"))):
               for indexinthelist in operatorlist:
                  if ((indexinthelist.type == "particle") and (indexinthelist.index == string.ascii_letters.index(index)+1)):
                     newamplitude.indexes.append(indexinthelist)
                     break
            elif (((index >= "i") and (index <= "o")) or ((index >= "I") and (index <= "O"))):
               for indexinthelist in operatorlist:
                  if ((indexinthelist.type == "hole") and (indexinthelist.index == string.ascii_letters.index(index)+1)):
                     newamplitude.indexes.append(indexinthelist)
                     break
            elif (((index >= "p") and (index <= "z")) or ((index >= "P") and (index <= "Z"))):
               for indexinthelist in operatorlist:
                  if ((indexinthelist.type == "general") and (indexinthelist.index == string.ascii_letters.index(index)+1)):
                     newamplitude.indexes.append(indexinthelist)
                     break
         else:
            if (index[0] == "p"):
               for indexinthelist in operatorlist:
                  if ((indexinthelist.type == "particle") and (indexinthelist.index == int(index[1:]))):
                     newamplitude.indexes.append(indexinthelist)
                     break
            elif (index[0] == "h"):
               for indexinthelist in operatorlist:
                  if ((indexinthelist.type == "hole") and (indexinthelist.index == int(index[1:]))):
                     newamplitude.indexes.append(indexinthelist)
                     break
            elif (index[0] == "g"):
               for indexinthelist in operatorlist:
                  if ((indexinthelist.type == "general") and (indexinthelist.index == int(index[1:]))):
                     newamplitude.indexes.append(indexinthelist)
                     break
            else:
               print "syntax error"
               stop
      amplitudes.append(newamplitude)

   newoperatorsequence = OperatorSequence(numericalfactor,summation,amplitudes,newsequences)
   return newoperatorsequence

def combinepermutations(one,two):
   """Connects two permutations of indexes"""
   if (len(one) != len(two)):
      print "Internal error"
      stop
   three = []
   for n in range(len(one)/2):
      three.append(one[n])
   for n in range(len(one)/2,len(one)):
      for m in range(len(two)/2):
         if (one[n].isidenticalto(two[m])):
            three.append(two[m+len(two)/2])
   return three

def isidenticalto(one,two):
   """Returns true if two permutations of indexes are identical"""
   if ((one == []) and (two == [])):
      return 1
   if (one == []):
      for m in range(len(two)/2):
         if (not two[m].isidenticalto(two[m+len(two)/2])):
            return 0
      return 1
   if (two == []):
      for m in range(len(one)/2):
         if (not one[m].isidenticalto(one[m+len(one)/2])):
            return 0
      return 1
   if (len(one) != len(two)):
      return 0
   for n in range(len(one)/2):
      found = 0
      for m in range(len(two)/2):
         if ((one[n].isidenticalto(two[m])) and (one[n+len(one)/2].isidenticalto(two[m+len(two)/2]))):
            found = 1
      if (not found):
         return 0
   return 1

class Operator:

   def __init__(self,type="unknown",dagger="unknown",index=0):
      """Creates a second-quantized hole/particle/general creation/annihilation operator"""
      self.type = type
      self.dagger = dagger
      self.index = index

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = string.join([self.type[0], repr(self.index)], "")
      if (self.dagger == "creation"):
         show = string.join([show, "+"], "")
      return show

   def tex(self):
      """Returns a LaTex form of output"""
      show = string.join([self.type[0],"_{",repr(self.index),"}"], "")
      if (self.dagger == "creation"):
         show = string.join([show, "^{\dagger}"], "")
      return show

   def duplicate(self):
      """Returns a deepcopy of self"""
      duplicate = Operator(self.type,self.dagger,self.index)
      return duplicate
      
   def isidenticalto(self,another):
      """Checks if two second-quantized operators are identical"""
      if ((self.type == another.type) and (self.dagger == another.dagger) and (self.index == another.index)):
         return 1
      else:
         return 0

   def issimilarto(self,another):
      """Checks if two second-quantized operators are similar"""
      if ((self.type == another.type) and (self.dagger == another.dagger)):
         return 1
      else:
         return 0
 
   def isin(self,list):
      """Returns true if an operator is in the list"""
      for index in list:
         if (self.isidenticalto(index)):
            return 1
      return 0

   def showwithoutdagger(self):
      """Returns a human-friendly string of the content"""
      show = string.join([self.type[0], repr(self.index)], "")
      return show

   def texwithoutdagger(self):
      """Returns a human-friendly string of the content"""
      show = string.join([self.type[0],"_{",repr(self.index),"}"], "")
      return show

   def isgreaterthan(self,another,operatorsequence):
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

      # at this point, self.type = another.type
      if ((not operatorsequence.summation.hastheindex(self)) and (not operatorsequence.summation.hastheindex(another))):
         if (self.index > another.index):
            return 1
         else:
            return 0
      elif (operatorsequence.summation.hastheindex(self) and (not operatorsequence.summation.hastheindex(another))):
         return 0
      elif ((not operatorsequence.summation.hastheindex(self)) and operatorsequence.summation.hastheindex(another)):
         return 1
      else:
         # at this point, self.type = another.type and both are summed over
         selfconnectivity = []
         anotherconnectivity = []
         for namplitude in range(len(operatorsequence.amplitudes)):
            amplitude = operatorsequence.amplitudes[namplitude]
            if (amplitude.hastheindex(self)):
               selfconnectivity.append(namplitude)
         for namplitude in range(len(operatorsequence.amplitudes)):
            amplitude = operatorsequence.amplitudes[namplitude]
            if (amplitude.hastheindex(another)):
               anotherconnectivity.append(namplitude)
         selfconnectivity.sort()
         anotherconnectivity.sort()
         if (selfconnectivity < anotherconnectivity):
            return 1
         elif (anotherconnectivity < selfconnectivity):
            return 0

      return 0

class Summation:

   def __init__(self,indexes=[]):
      """Creates a summation"""
      self.indexes = indexes

   def __str__(self):
      """Print the amplitude"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = "Sum ("
      for index in self.indexes:
         show = string.join([show, index.showwithoutdagger()])
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
         show = string.join([show,index.texwithoutdagger()])
      show = string.join([show,"}"])
      return show

   def duplicate(self):
      """Returns a deepcopy of itself"""
      duplicate = Summation([])
      for index in self.indexes:
         duplicate.indexes.append(index.duplicate())
      return duplicate

   def hasthesameform(self,another):
      """Checks if two summations have the same numbers of holes, particles, and generals"""
      nself = 0
      for operator in self.indexes:
         if (operator.type == "hole"):
            nself = nself + 1
      nanother = 0
      for operator in another.indexes:
         if (operator.type == "hole"):
            nanother = nanother + 1
      if (nself != nanother):
         return 0
      nself = 0
      for operator in self.indexes:
         if (operator.type == "particle"):
            nself = nself + 1
      nanother = 0
      for operator in another.indexes:
         if (operator.type == "particle"):
            nanother = nanother + 1
      if (nself != nanother):
         return 0
      nself = 0
      for operator in self.indexes:
         if (operator.type == "general"):
            nself = nself + 1
      nanother = 0
      for operator in another.indexes:
         if (operator.type == "general"):
            nanother = nanother + 1
      if (nself != nanother):
         return 0
      return 1

   def isidenticalto(self,another):
      """Returns true if two summations are identical"""
      if (len(self.indexes) != len(another.indexes)):
         return 0
      else:
         for nindex in range(len(self.indexes)):
            selfindex = self.indexes[nindex]
            anotherindex = another.indexes[nindex]
            if (not selfindex.isidenticalto(anotherindex)):
               return 0
      return 1
 
   def hastheindex(self,another):
      """Returns true if the summation has the input index"""
      has = 0
      for index in self.indexes:
         if (index.isidenticalto(another)):
            has = 1
      return has

class Amplitude:

   def __init__(self,type="unknown",indexes=[],index=0,conjugate=0):
      """Creates an integral/amplitude"""
      self.type = type
      self.indexes = indexes
      self.index = index
      self.conjugate = conjugate

   def __str__(self):
      """Print the amplitude"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = self.type
      if (self.conjugate):
         show = string.join([show, "+"],"")
      show = string.join([show, "("])
      for index in self.indexes:
         show = string.join([show, index.showwithoutdagger()])
      show = string.join([show,")"])
      return show
 
   def hastheindex(self,another):
      """Returns true if the summation has the input index"""
      has = 0
      for index in self.indexes:
         if (index.isidenticalto(another)):
            has = 1
      return has

   def tex(self):
      """Returns a LaTeX string of the content"""
      show = self.type
      show = string.join([show, "^{"])
      for index in self.indexes[0:len(self.indexes)/2]:
         show = string.join([show, index.texwithoutdagger()])
      show = string.join([show,"}"])
      show = string.join([show, "_{"])
      for index in self.indexes[len(self.indexes)/2:len(self.indexes)]:
         show = string.join([show, index.texwithoutdagger()])
      show = string.join([show,"}"])
      if (self.conjugate):
         show = string.join(["\\left(",show,"\\right)^{\\dagger}"],"")
      return show
      
   def duplicate(self):
      """Returns a deepcopy of itself"""
      duplicate = Amplitude(self.type,[],self.index,self.conjugate)
      for index in self.indexes:
         duplicate.indexes.append(index.duplicate())
      return duplicate

   def hasthesameform(self,another):
      """Checks if two amplitude sets have the same numbers of holes, particles, and generals"""
      if (self.type != another.type):
         return 0
      if (len(self.indexes) != len(another.indexes)):
         return 0
      if (self.conjugate != another.conjugate):
         return 0
      nself = 0
      for operator in self.indexes:
         if (operator.type == "hole"):
            nself = nself + 1
      nanother = 0
      for operator in another.indexes:
         if (operator.type == "hole"):
            nanother = nanother + 1
      if (nself != nanother):
         return 0
      nself = 0
      for operator in self.indexes:
         if (operator.type == "particle"):
            nself = nself + 1
      nanother = 0
      for operator in another.indexes:
         if (operator.type == "particle"):
            nanother = nanother + 1
      if (nself != nanother):
         return 0
      nself = 0
      for operator in self.indexes:
         if (operator.type == "general"):
            nself = nself + 1
      nanother = 0
      for operator in another.indexes:
         if (operator.type == "general"):
            nanother = nanother + 1
      if (nself != nanother):
         return 0
      return 1

   def isidenticalto(self,another):
      """Returns true if two amplitudes are identical"""
      if (self.type != another.type):
         return 0
      elif (len(self.indexes) != len(another.indexes)):
         return 0
      elif (self.conjugate != another.conjugate):
         return 0
      else:
         for nindex in range(len(self.indexes)):
            selfindex = self.indexes[nindex]
            anotherindex = another.indexes[nindex]
            if (not selfindex.isidenticalto(anotherindex)):
               return 0
      return 1

   def isgreaterthan(self,another,operatorsequence):
      """Returns true if self should be to the right of another in the canonical order"""
      
      # count the number of like amplitudes in operatorsequence
      nself = 0
      for amplitude in operatorsequence.amplitudes:
         if ((amplitude.type == self.type) and (len(amplitude.indexes) == len(self.indexes)) and (amplitude.conjugate == self.conjugate)):
            nself = nself + 1
      nanother = 0
      for amplitude in operatorsequence.amplitudes:
         if ((amplitude.type == another.type) and (len(amplitude.indexes) == len(another.indexes)) and (amplitude.conjugate == another.conjugate)):
            nanother = nanother + 1
      if (nself > nanother):
         return 0
      elif (nself < nanother):
         return 1
      # conjugate
      if (self.conjugate > another.conjugate):
         return 0
      elif (self.conjugate < another.conjugate):
         return 1
      # type
      if (self.type > another.type):
         return 1
      elif (self.type < another.type):
         return 0
      # number of indexes
      if (len(self.indexes) < len(another.indexes)):
         return 1
      elif (len(self.indexes) > len(another.indexes)):
         return 0
      # number of external indexes
      nself = 0
      iself = 9999999999
      for operator in self.indexes:
         if (not operatorsequence.summation.hastheindex(operator)):
            nself = nself + 1
            if (iself > operator.index):
               iself = operator.index
      nanother = 0
      ianother = 9999999999
      for operator in another.indexes:
         if (not operatorsequence.summation.hastheindex(operator)):
            nanother = nanother + 1
            if (ianother > operator.index):
               ianother = operator.index
      if (nself < nanother):
         return 1
      elif (nself > nanother):
         return 0
      # earliest external indexes
#     if (nself > 0):
#        if (ianother > iself):
#           return 0
#        elif (ianother < iself):
#           return 1
      # connectivity
      selfconnectivity = []
      anotherconnectivity = []
      for operator in self.indexes:
         if (operatorsequence.summation.hastheindex(operator)):
            for amplitude in operatorsequence.amplitudes:
               if (amplitude.hastheindex(operator)):
                  amplitudesymbol = amplitude.type + repr(len(amplitude.indexes))
                  if (amplitude.conjugate):
                     amplitudesymbol = amplitudesymbol + "+"
                  selfconnectivity.append(amplitudesymbol)
      for operator in another.indexes:
         if (operatorsequence.summation.hastheindex(operator)):
            for amplitude in operatorsequence.amplitudes:
               if (amplitude.hastheindex(operator)):
                  amplitudesymbol = amplitude.type + repr(len(amplitude.indexes))
                  if (amplitude.conjugate):
                     amplitudesymbol = amplitudesymbol + "+"
                  anotherconnectivity.append(amplitudesymbol)
      selfconnectivity.sort()
      anotherconnectivity.sort()
      if (selfconnectivity < anotherconnectivity):
         return 1
      elif (anotherconnectivity < selfconnectivity):
         return 0
      return 0

   def canonicalize(self,operatorsequence):
      """Reorder the indexes in the canonical order"""

      another = self.duplicate()
      parity = 1
      done = 0
      while (not done):
         done = 1
         # reorder super indexes
         for noperatora in range(len(another.indexes)/2):
            for noperatorb in range(len(another.indexes)/2):
               if (noperatora >= noperatorb):
                  continue
               operatora = another.indexes[noperatora]
               operatorb = another.indexes[noperatorb]
               if (operatora.isgreaterthan(operatorb,operatorsequence)):
                  another.indexes[noperatorb] = copy.deepcopy(operatora)
                  another.indexes[noperatora] = copy.deepcopy(operatorb)
                  parity = parity * (-1)
                  done = 0
      done = 0
      while (not done):
         done = 1
         # reorder sub indexes
         for noperatora in range(len(another.indexes)/2,len(another.indexes)):
            for noperatorb in range(len(another.indexes)/2,len(another.indexes)):
               if (noperatora >= noperatorb):
                  continue
               operatora = another.indexes[noperatora]
               operatorb = another.indexes[noperatorb]
               if (operatora.isgreaterthan(operatorb,operatorsequence)):
                  another.indexes[noperatorb] = copy.deepcopy(operatora)
                  another.indexes[noperatora] = copy.deepcopy(operatorb)
                  parity = parity * (-1)
                  done = 0

      return [another,parity]

class Factor:

   def __init__(self,coefficients=[],permutations=[]):
      """Creates a numerical and permutation factor of an operator sequence"""
      self.coefficients = coefficients
      self.permutations = copy.deepcopy(permutations)

   def __str__(self):
      """Prints the content"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of contests"""
      show = "["
      for n in range(len(self.coefficients)):
         coefficient = self.coefficients[n]
         # str() rounds a float after 12 digits, while repr() after 17,
         # so the former tends to give a more pleasant expression.
#        num = rationaltofractional(coefficient)[0]
#        den = rationaltofractional(coefficient)[1]
#        if (num >= 0):
#           show = string.join([show,"+",repr(num)])
#        elif (num < 0):
#           show = string.join([show,"-",repr(-num)])
#        if (den != 1):
#           show = string.join([show,"/",repr(den)],"")
         if (coefficient >= 0.0):
            show = string.join([show,"+",str(coefficient)])
         elif (coefficient < 0.0):
            show = string.join([show,"-",str(-coefficient)])
         if (self.permutations[n]):
            show = string.join([show,"* P("])
            for noperator in range(len(self.permutations[n])/2):
               operator = self.permutations[n][noperator]
               show = string.join([show,operator.showwithoutdagger()])
            show = string.join([show,"=>"])
            for noperator in range(len(self.permutations[n])/2,len(self.permutations[n])):
               operator = self.permutations[n][noperator]
               show = string.join([show,operator.showwithoutdagger()])
            show = string.join([show,")"])
      show = string.join([show,"]"])
      return show

   def tex(self):
      """Returns a LaTeX string of contests"""
      coefficient = self.coefficients[0]
      for n in range(len(self.coefficients)):
         if (abs(self.coefficients[n]) != abs(coefficient)):
            raise RuntimeError, "unrealistic factor"
      fraction = abs(int(1.0/coefficient))
      if (1.0/float(fraction) != abs(coefficient)):
         print " !!! WARNING !!! inaccurate arithmatic"
      if (fraction == 1):
         frac = ""
      else:
         frac = string.join(["\\frac{1}{",str(fraction),"}"],"")
      if (coefficient >= 0.0):
         show = string.join(["+",frac])
      elif (coefficient < 0.0):
         show = string.join(["-",frac])
      if (len(self.coefficients) > 1):
         show = string.join([show,"\\left("],"")
         for n in range(len(self.coefficients)):
            if (self.coefficients[n]/coefficient > 0.0):
               show = string.join([show,"+"],"")
            else:
               show = string.join([show,"-"],"")
            if (self.permutations[n]):
               show = string.join([show,"P^{"])
               for nindex in range(len(self.permutations[n])/2,3*len(self.permutations[n])/4):
                  index = self.permutations[n][nindex]
                  show = string.join([show,index.texwithoutdagger()])
               for nindex in range(len(self.permutations[n])/4,len(self.permutations[n])/2):
                  index = self.permutations[n][nindex]
                  show = string.join([show,index.texwithoutdagger()])
               show = string.join([show,"}_{"])
               for nindex in range(len(self.permutations[n])/4):
                  index = self.permutations[n][nindex]
                  show = string.join([show,index.texwithoutdagger()])
               for nindex in range(3*len(self.permutations[n])/4,len(self.permutations[n])):
                  index = self.permutations[n][nindex]
                  show = string.join([show,index.texwithoutdagger()])
               show = string.join([show,"}"])
            else:
               show = string.join([show,"1"],"")
         show = string.join([show,"\\right)"])
      return show

   def multiply(self,factor):
      """Multiply a factor to all coefficients"""
      for n in range(len(self.coefficients)):
         self.coefficients[n] = self.coefficients[n] * factor

   def add(self,another,factor=1.0):
      """Add two Factors together"""
      for m in range(len(another.coefficients)):
         done = 0
         for n in range(len(self.coefficients)):
            if (isidenticalto(self.permutations[n],another.permutations[m])):
               if ((self.coefficients[n] < 0.0) and (another.coefficients[m] * factor > 0.0)):
                  print " ! Warning ! cancellation of terms occurred "
               if ((self.coefficients[n] > 0.0) and (another.coefficients[m] * factor < 0.0)):
                  print " ! Warning ! cancellation of terms occurred "
               self.coefficients[n] = self.coefficients[n] + another.coefficients[m] * factor
               done = 1
         if (not done):
            self.coefficients.append(another.coefficients[m] * factor)
            self.permutations.append(another.permutations[m])

class OperatorSequence:

   def __init__(self,factor=[],summation=[],amplitudes=[],sequence=[]):
      """Creates a sequence of normal ordered second-quantized operators with some numerical factor, amplitudes, and summation"""
      self.factor = factor
      self.summation = summation
      self.amplitudes = amplitudes
      self.sequence = sequence

   def __str__(self):
      """Prints the sequence of operator contractions"""
      return self.show()

   def show(self):
      """Returns a human-friendly string of the content"""
      show = self.factor.show()
      if (self.summation):
         if (len(self.summation.indexes) > 0):
            show = string.join([show, "*", self.summation.show()])
      for index in self.amplitudes:
         show = string.join([show, "*", index.show()])
      if (self.sequence):
         show = string.join([show, "* <0|"])
         for sequence in self.sequence:
            show = string.join([show, "{"])
            for operator in sequence:
               show = string.join([show, operator.show()])
            show = string.join([show, "}"])
         show = string.join([show, "|0>"])
      return show

   def tex(self):
      """Returns a LaTeX string of the content"""
      show = self.factor.tex()
#     if (self.summation):
#        if (len(self.summation.indexes) > 0):
#           show = string.join([show, self.summation.tex()])
      for index in self.amplitudes:
         show = string.join([show, index.tex()])
      if (self.sequence):
         show = string.join([show, "\\langle 0 |"])
         for sequence in self.sequence:
            show = string.join([show, "\{"])
            for operator in sequence:
               show = string.join([show, operator.tex()])
            show = string.join([show, "\}"])
         show = string.join([show, "|0\\rangle"])
      return show

   def duplicate(self):
      """Makes a copy of itself"""
      duplicate = OperatorSequence()
      duplicate.factor = copy.deepcopy(self.factor)
      duplicate.summation = copy.deepcopy(self.summation)
      duplicate.amplitudes = copy.deepcopy(self.amplitudes)
      duplicate.sequence = copy.deepcopy(self.sequence)
      return duplicate

   def writetofile(self,filename):
      """Writes the output to a given file"""
      file = open(filename,"w")
      file.write(self.show())
      file.write("\n")

   def removeemptycurly(self):
      """Eliminates all empty curly brackets (curly means a sequence of normal ordered operator in {})"""

      hasempty = 0
      for ncurly in range(len(self.sequence)):
         curly = self.sequence[ncurly]
         if (not curly):
            del self.sequence[ncurly]
            hasempty = 1
            break

      if (hasempty):
         self.removeemptycurly()
      else:
         return self

   def alreadycontracted(self):
      """Checks if an operator sequence object is fully contracted"""

      # first, we delete all empty {} just in case
      self.removeemptycurly()

      # already fully contracted?
      if (not self.sequence):
         return 1
      else:
         return 0

   def isunabletocontract(self):
      """Counts the number of operators and determine if it is possible to give nonzero contraction at the end"""

      # count the number of hole/particle/general creation/annihilation operators
      nholecreation = 0
      nholeannihilation = 0
      nparticlecreation = 0
      nparticleannihilation = 0
      ngeneralcreation = 0
      ngeneralannihilation = 0
      for sequence in self.sequence:
         for operator in sequence:
            if ((operator.type == "hole") and (operator.dagger == "creation")):
               nholecreation = nholecreation + 1
            elif ((operator.type == "hole") and (operator.dagger == "annihilation")):
               nholeannihilation = nholeannihilation + 1
            if ((operator.type == "particle") and (operator.dagger == "creation")):
               nparticlecreation = nparticlecreation + 1
            elif ((operator.type == "particle") and (operator.dagger == "annihilation")):
               nparticleannihilation = nparticleannihilation + 1
            if ((operator.type == "general") and (operator.dagger == "creation")):
               ngeneralcreation = ngeneralcreation + 1
            elif ((operator.type == "general") and (operator.dagger == "annihilation")):
               ngeneralannihilation = ngeneralannihilation + 1

      # see if enough operators remain for contractions to survive
      uncontractable = 0
      if (nholecreation + ngeneralcreation < nholeannihilation):
         uncontractable = 1
      if (nholeannihilation + ngeneralannihilation < nholecreation):
         uncontractable = 1
      if (nparticlecreation + ngeneralcreation < nparticleannihilation):
         uncontractable = 1
      if (nparticleannihilation + ngeneralannihilation < nparticlecreation):
         uncontractable = 1
      return uncontractable

   def performcontraction(self):
      """Perform a contraction of the left-most operator"""

      # result will be a list of new operator sequence objects
      result = ListOperatorSequences()

      # already fully contracted?
      if (self.alreadycontracted()):
         newsequence = self.duplicate()
         result.add(newsequence)
         return result

      # no way to contract?
      elif ((len(self.sequence) == 1) or (self.isunabletocontract())):
         return result

      # get the left-most operator
      leftmost = self.sequence[0][0] 
      
      # loop over other {}
      for ncurly in range(len(self.sequence)):
         curly = self.sequence[ncurly]
         if (ncurly == 0):
            continue
         for noperator in range(len(self.sequence[ncurly])):
            operator = curly[noperator]
            
            # only allowed contractions are {h+}{h} and {p}{p+}
            if (leftmost.dagger == operator.dagger):
               continue
            elif ((leftmost.type == "hole") and (operator.type == "particle")):
               continue
            elif ((leftmost.type == "particle") and (operator.type == "hole")):
               continue
            elif ((leftmost.type == "hole") and (leftmost.dagger == "annihilation")):
               continue
            elif ((leftmost.type == "particle") and (leftmost.dagger == "creation")):
               continue
            elif ((operator.type == "hole") and (operator.dagger == "creation")):
               continue
            elif ((operator.type == "particle") and (operator.dagger == "annihilation")):
               continue

            # check if the indexes can be made to match by virtue of summation
            exist = "neither"
            for index in self.summation.indexes:
               if (leftmost.isidenticalto(index)):
                  exist = "leftmost"
            if (exist == "neither"):
               for index in self.summation.indexes:
                  if (operator.isidenticalto(index)):
                     exist = "operator"
            if (exist == "leftmost"):

               # now contraction is possible --- add a new operator sequence object to result
               newsequence = self.duplicate()

               # delete leftmost from the summation indexes
               for index in newsequence.summation.indexes:
                  if (leftmost.isidenticalto(index)):
                     del newsequence.summation.indexes[newsequence.summation.indexes.index(index)]

               # count the number of operators between leftmost and the current operator and determine the parity
               length = len(self.sequence[0][1:]) + curly.index(operator)
               for anothercurly in self.sequence[1:]:
                  if (anothercurly == curly):
                     break
                  else:
                     length = length + len(anothercurly)
               parity = (-1)**length
               newsequence.factor.multiply(parity)

               # delete the contracted pair from the sequence
               del newsequence.sequence[0][0]
               del newsequence.sequence[ncurly][noperator]

               # replace any appearance of leftmost by operator
               if (operator.type == 'general'):
                  for namplitude in range(len(newsequence.amplitudes)):
                     amplitude = newsequence.amplitudes[namplitude]
                     for nindex in range(len(amplitude.indexes)):
                        index = amplitude.indexes[nindex]
                        if (operator.isidenticalto(index)):
                           newsequence.amplitudes[namplitude].indexes[nindex] = leftmost
                  for nindex in range(len(newsequence.summation.indexes)):
                     index = newsequence.summation.indexes[nindex]
                     if (operator.isidenticalto(index)):
                        newsequence.summation.indexes[nindex] = leftmost
               else:
                  for namplitude in range(len(newsequence.amplitudes)):
                     amplitude = newsequence.amplitudes[namplitude]
                     for nindex in range(len(amplitude.indexes)):
                        index = amplitude.indexes[nindex]
                        if (leftmost.isidenticalto(index)):
                           newsequence.amplitudes[namplitude].indexes[nindex] = operator
                  for nindex in range(len(newsequence.summation.indexes)):
                     index = newsequence.summation.indexes[nindex]
                     if (leftmost.isidenticalto(index)):
                        newsequence.summation.indexes[nindex] = operator

               # cleanup the empty brackets
               newsequence.removeemptycurly()

               # add to the result
               result.add(newsequence)
               
            elif (exist == "operator"): 

               # contraction is again possible --- add a new operator sequence object to result
               newsequence = self.duplicate()

               # delete operator from the summation indexes
               for index in newsequence.summation.indexes:
                  if (operator.isidenticalto(index)):
                     del newsequence.summation.indexes[newsequence.summation.indexes.index(index)]

               # count the number of operators between leftmost and the current operator and determine the parity
               length = len(self.sequence[0][1:]) + curly.index(operator)
               for anothercurly in self.sequence[1:]:
                  if (anothercurly == curly):
                     break
                  else:
                     length = length + len(anothercurly)
               parity = (-1)**length
               newsequence.factor.multiply(parity)

               # delete the contracted pair from the sequence
               del newsequence.sequence[0][0]
               del newsequence.sequence[ncurly][noperator]

               # replace any appearance of operator by leftmost
               if (leftmost.type == 'general'):
                  for namplitude in range(len(newsequence.amplitudes)):
                     amplitude = newsequence.amplitudes[namplitude]
                     for nindex in range(len(amplitude.indexes)):
                        index = amplitude.indexes[nindex]
                        if (leftmost.isidenticalto(index)):
                           newsequence.amplitudes[namplitude].indexes[nindex] = operator
                  for nindex in range(len(newsequence.summation.indexes)):
                     index = newsequence.summation.indexes[nindex]
                     if (leftmost.isidenticalto(index)):
                        newsequence.summation.indexes[nindex] = operator
               else:
                  for namplitude in range(len(newsequence.amplitudes)):
                     amplitude = newsequence.amplitudes[namplitude]
                     for nindex in range(len(amplitude.indexes)):
                        index = amplitude.indexes[nindex]
                        if (operator.isidenticalto(index)):
                           newsequence.amplitudes[namplitude].indexes[nindex] = leftmost
                  for nindex in range(len(newsequence.summation.indexes)):
                     index = newsequence.summation.indexes[nindex]
                     if (operator.isidenticalto(index)):
                        newsequence.summation.indexes[nindex] = leftmost

               # cleanup the empty brackets
               newsequence.removeemptycurly()

               # add to the result
               result.add(newsequence)

            else:
               break

      return result

   def performfullcontraction(self):
      """Performs full contraction of a given operator sequence and returns a list of tensor contractions"""

      print self.show()
      print " ... commencing full operator contraction"

      # result will be a list of tensor contractions (operator sequence objects with empty operator sequence)
      result = ListOperatorSequences()
      result.add(self)

      # see if already fully contracted
      done = self.alreadycontracted()

      # recursive execution of performcontraction()
      iteration = 0
      while (not done):
         iteration = iteration + 1
         newresult = ListOperatorSequences()
         for halfwaycontracted in result.list:
            newaddition = halfwaycontracted.performcontraction()
            if (newaddition):
               newresult.join(newaddition)
         newresult.simplifyone()
         numberofterms = len(newresult.list)
         print " ... iteration = %d, number of terms = %d" %(iteration, numberofterms)
         done = 1
         for halfwaycontracted in newresult.list:
            if (not halfwaycontracted.alreadycontracted()):
               done = 0
         result = newresult.duplicate()

      return result

   def hasthesameform(self,another):
      """Checks if two operator sequences have the same form for possible consolidation"""
      if (not self.summation.hasthesameform(another.summation)):
         return 0
      if (len(self.amplitudes) != len(another.amplitudes)):
         return 0
      else:
         for namplitude in range(len(self.amplitudes)):
            if (not self.amplitudes[namplitude].hasthesameform(another.amplitudes[namplitude])):
               return 0
      if (len(self.sequence) != len(another.sequence)):
         return 0
      else:
         for nsequence in range(len(self.sequence)):
            nself = 0
            for operator in self.sequence[nsequence]:
               if (operator.type == "hole"):
                   nself = nself + 1
            nanother = 0
            for operator in another.sequence[nsequence]:
               if (operator.type == "hole"):
                  nanother = nanother + 1
            if (nself != nanother):
               return 0
            nself = 0
            for operator in self.sequence[nsequence]:
               if (operator.type == "particle"):
                  nself = nself + 1
            nanother = 0
            for operator in another.sequence[nsequence]:
               if (operator.type == "particle"):
                  nanother = nanother + 1
            if (nself != nanother):
               return 0
            nself = 0
            for operator in self.sequence[nsequence]:
               if (operator.type == "general"):
                  nself = nself + 1
            nanother = 0
            for operator in another.sequence[nsequence]:
               if (operator.type == "general"):
                  nanother = nanother + 1
            if (nself != nanother):
               return 0
      return 1

   def isidenticalto(self,another):
      """Returns true if two operator sequences are identical except for the factor"""

      if (not self.summation.isidenticalto(another.summation)):
         return 0
      if (len(self.amplitudes) != len(another.amplitudes)):
         return 0
      if (len(self.sequence) != len(another.sequence)):
         return 0
      for namplitude in range(len(self.amplitudes)):
         if (not self.amplitudes[namplitude].isidenticalto(another.amplitudes[namplitude])):
            return 0
      for nsequence in range(len(self.sequence)):
         if (len(self.sequence[nsequence]) != len(another.sequence[nsequence])):
            return 0
         else:
            for noperator in range(len(self.sequence[nsequence])):
               if (not self.sequence[nsequence][noperator].isidenticalto(another.sequence[nsequence][noperator])):
                  return 0
      return 1

   def has(self,index):
      """Checks if a certain index is included in an operator sequence"""

      # see if the index is in summation indexes
      for another in self.summation.indexes:
         if (another.isidenticalto(index)):
            return 1

      # see if the index is in amplitude indexes
      for amplitude in self.amplitudes:
         for another in amplitude.indexes:
            if (another.isidenticalto(index)):
               return 1

      # see if the index is in the operator sequences
      for sequence in self.sequence:
         for another in sequence:
            if (another.isidenticalto(index)):
               return 1

      # not included
      return 0

   def relabels(self,another):
      """Relabels the operator indexes to help consolidate terms"""

      if (not self.hasthesameform(another)):
         return another

      else:

         # find a lone index in summation indexes
         for index in self.summation.indexes:
            if (not another.has(index)):
               for anotherindex in another.summation.indexes:
                  if ((not self.has(anotherindex)) and (index.issimilarto(anotherindex))):
                     
                     # at this point, we know that we should relabel anotherindex by index everywhere in another
                     for nyetanother in range(len(another.summation.indexes)):
                        yetanother = another.summation.indexes[nyetanother]
                        if (yetanother.isidenticalto(anotherindex)):
                           another.summation.indexes[nyetanother] = copy.deepcopy(index)
  
                     for namplitude in range(len(another.amplitudes)):
                        amplitude = another.amplitudes[namplitude]
                        for nyetanother in range(len(amplitude.indexes)):
                           yetanother = amplitude.indexes[nyetanother]
                           if (yetanother.isidenticalto(anotherindex)):
                              another.amplitudes[namplitude].indexes[nyetanother] = copy.deepcopy(index)
 
                     for nsequence in range(len(another.sequence)):
                        sequence = another.sequence[nsequence]
                        for nyetanother in range(len(sequence)):
                           yetanother = sequence[nyetanother]
                           if (yetanother.isidenticalto(anotherindex)):
                              another.sequence[nsequence][nyetanother] = copy.deepcopy(index)

                     return another

         # find a lone index in amplitude indexes
         for selfamplitude in self.amplitudes:
            for index in selfamplitude.indexes:
               if (not another.has(index)):
                  for anotheramplitude in another.amplitudes:
                     for anotherindex in anotheramplitude.indexes:
                        if ((not self.has(anotherindex)) and (index.issimilarto(anotherindex))):
                     
                           # at this point, we know that we should relabel anotherindex by index everywhere in another
                           for nyetanother in range(len(another.summation.indexes)):
                              yetanother = another.summation.indexes[nyetanother]
                              if (yetanother.isidenticalto(anotherindex)):
                                 another.summation.indexes[nyetanother] = copy.deepcopy(index)

                           for namplitude in range(len(another.amplitudes)):
                              amplitude = another.amplitudes[namplitude]
                              for nyetanother in range(len(amplitude.indexes)):
                                 yetanother = amplitude.indexes[nyetanother]
                                 if (yetanother.isidenticalto(anotherindex)):
                                    another.amplitudes[namplitude].indexes[nyetanother] = copy.deepcopy(index)

                           for nsequence in range(len(another.sequence)):
                              sequence = another.sequence[nsequence]
                              for nyetanother in range(len(sequence)):
                                 yetanother = sequence[nyetanother]
                                 if (yetanother.isidenticalto(anotherindex)):
                                    another.sequence[nsequence][nyetanother] = copy.deepcopy(index)

                           return another

         # find a lone index in operator sequences
         for selfsequence in self.sequence:
            for index in selfsequence:
               if (not another.has(index)):
                  for anothersequence in another.sequence:
                     for anotherindex in anothersequence:
                        if ((not self.has(anotherindex)) and (index.issimilarto(anotherindex))):
                        
                           # at this point, we know that we should relabel anotherindex by index everywhere in another
                           for nyetanother in range(len(another.summation.indexes)):
                              yetanother = another.summation.indexes[nyetanother]
                              if (yetanother.isidenticalto(anotherindex)):
                                 another.summation.indexes[nyetanother] = copy.deepcopy(index)
  
                           for namplitude in range(len(another.amplitudes)):
                              amplitude = another.amplitudes[namplitude]
                              for nyetanother in range(len(amplitude.indexes)):
                                 yetanother = amplitude.indexes[nyetanother]
                                 if (yetanother.isidenticalto(anotherindex)):
                                    another.amplitudes[namplitude].indexes[nyetanother] = copy.deepcopy(index)

                           for nsequence in range(len(another.sequence)):
                              sequence = another.sequence[nsequence]
                              for nyetanother in range(len(sequence)):
                                 yetanother = sequence[nyetanother]
                                 if (yetanother.isidenticalto(anotherindex)):
                                    another.sequence[nsequence][nyetanother] = copy.deepcopy(index)

                           return another

         return another

   def fullyrelabels(self,another):
      """Relabels the operator indexes to help consolidate terms"""

      if (not self.hasthesameform(another)):
         return another
      else:
         done = 0
         while (not done):
            another = self.relabels(another)
            done = self.hasnomismatch(another)
    
         return another

   def hasnomismatch(self,another):
      """Returns 1 if there is no index in self that does not exist in another"""

      if (not self.hasthesameform(another)):
         return 0

      nomismatch = 1
      for index in self.summation.indexes:
         if (not another.has(index)):
             nomismatch = 0
      for selfamplitude in self.amplitudes:
         for index in selfamplitude.indexes:
            if (not another.has(index)):
               nomismatch = 0
      for selfsequence in self.sequence:
         for index in selfsequence:
            if (not another.has(index)):
               nomismatch = 0

      return nomismatch

   def canmerge(self,another):
      """Returns 1 if another operator sequence can be merged to itself"""

      # do they have any index mismatch?
      if (not self.hasnomismatch(another)):
         return 0

      # do they have the identical operator sequences?
      for nsequence in range(len(self.sequence)):
         selfsequence = self.sequence[nsequence]
         anothersequence = another.sequence[nsequence]
         for noperator in range(len(selfsequence)):
            selfoperator = selfsequence[noperator]
            anotheroperator = anothersequence[noperator]
            if (not selfoperator.isidenticalto(anotheroperator)):
               return 0

      # do they have the summation indexes that do not differ by more than just permutation?
      for selfindex in self.summation.indexes:
         exist = 0
         for anotherindex in another.summation.indexes:
            if (anotherindex.isidenticalto(selfindex)):
               exist = 1
         if (not exist):
            return 0

      # do they have the amplitude indexes that do not differ by more than just permutation?
      for namplitude in range(len(self.amplitudes)):
         selfamplitude = self.amplitudes[namplitude].indexes
         anotheramplitude = another.amplitudes[namplitude].indexes
         for selfindex in selfamplitude:
            exist = 0
            for anotherindex in anotheramplitude:
               if (anotherindex.isidenticalto(selfindex)):
                  exist = 1
            if (not exist):
               return 0 

      return 1

   def merges(self,another):
      """Merges another operator sequence to itself when possible"""

      # parity of a permutation can be computed as the product of parities of all pairwise permutations
      parity = 1.0

      # determine the parity for amplitudes
      for namplitude in range(len(self.amplitudes)): 
         selfamplitude = self.amplitudes[namplitude]
         anotheramplitude = another.amplitudes[namplitude] 
         for nselfindexa in range(len(selfamplitude.indexes)):
            selfindexa = selfamplitude.indexes[nselfindexa]
            for nselfindexb in range(len(selfamplitude.indexes)):
               if (nselfindexb <= nselfindexa):
                  continue
               selfindexb = selfamplitude.indexes[nselfindexb]
               for nanotherindexa in range(len(anotheramplitude.indexes)):
                  anotherindexa = anotheramplitude.indexes[nanotherindexa]
                  if (anotherindexa.isidenticalto(selfindexa)):
                     for nanotherindexb in range(len(anotheramplitude.indexes)):
                        anotherindexb = anotheramplitude.indexes[nanotherindexb]
                        if (anotherindexb.isidenticalto(selfindexb)):
                           if (nanotherindexb < nanotherindexa):
                              parity = parity * (-1.0)

      self.factor.add(another.factor, parity)

      return self

   def swapoperators(self,indexa,indexb):
      """Swap indexa and indexb everywhere they appear in self"""

      for nindex in range(len(self.summation.indexes)):
         index = self.summation.indexes[nindex]
         if (index.isidenticalto(indexa)):
            self.summation.indexes[nindex] = indexb
         elif (index.isidenticalto(indexb)):
            self.summation.indexes[nindex] = indexa

      for namplitude in range(len(self.amplitudes)):
         amplitude = self.amplitudes[namplitude]
         for nindex in range(len(amplitude.indexes)):
            index = amplitude.indexes[nindex]
            if (index.isidenticalto(indexa)):
               self.amplitudes[namplitude].indexes[nindex] = indexb
            elif (index.isidenticalto(indexb)):
               self.amplitudes[namplitude].indexes[nindex] = indexa

      for nsequence in range(len(self.sequence)):
         sequence = self.sequence[nsequence]
         for nindex in range(len(sequence)):
            index = sequence[nindex]
            if (index.isidenticalto(indexa)):
               self.sequence[nsequence][nindex] = indexb
            elif (index.isidenticalto(indexb)):
               self.sequence[nsequence][nindex] = indexa

      return self

   def swapamplitudes(self,namplitudea,namplitudeb):
      """Swap two amplitudes in self"""

      swap = copy.deepcopy(self.amplitudes[namplitudea])
      self.amplitudes[namplitudea] = copy.deepcopy(self.amplitudes[namplitudeb])
      self.amplitudes[namplitudeb] = copy.deepcopy(swap)
      
      return self

   def targetindexpermutation(self):
      """Returns a list of all possible permutations and redundancy of target indexes of self"""

      # generate a target tensor
      super = []
      sub = []
      for tensor in self.amplitudes:
         for nindex in range(len(tensor.indexes)/2):
            index = tensor.indexes[nindex]
            common = 0
            if (self.summation):
               for another in self.summation.indexes:
                  if (index.isidenticalto(another)):
                     common = 1
            if (not common):
               super.append(tensor.indexes[nindex])
         for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
            index = tensor.indexes[nindex]
            common = 0
            if (self.summation):
               for another in self.summation.indexes:
                  if (index.isidenticalto(another)):
                     common = 1
            if (not common):
               sub.append(tensor.indexes[nindex])

      # permutation
      result = ListOperatorSequences()
      result.add(self.duplicate())
      for nsupera in range(len(super)-1):
         result = result.targetsuperpermutation(nsupera)
      for nsuba in range(len(sub)-1):
         result = result.targetsubpermutation(nsuba)
      for operatorsequence in result.list:
         newsuper = []
         newsub = []
         for tensor in operatorsequence.amplitudes:
            for nindex in range(len(tensor.indexes)/2):
               index = tensor.indexes[nindex]
               common = 0
               if (operatorsequence.summation):
                  for another in operatorsequence.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  newsuper.append(tensor.indexes[nindex])
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               index = tensor.indexes[nindex]
               common = 0
               if (operatorsequence.summation):
                  for another in operatorsequence.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  newsub.append(tensor.indexes[nindex])
         for ncoeff in range(len(operatorsequence.factor.coefficients)):
            operatorsequence.factor.coefficients[ncoeff] = self.factor.coefficients[ncoeff]
            newpermutation = newsuper + newsub + super + sub
            if (operatorsequence.factor.permutations[ncoeff] == []):
               operatorsequence.factor.permutations[ncoeff] = newpermutation
            else:
               operatorsequence.factor.permutations[ncoeff] = combinepermutations(newpermutation, operatorsequence.factor.permutations[ncoeff])
      return result

   def canonicalize(self):
      """Reorder amplitudes and common indexes in the canonical order"""
      # In canonical order, amplitudes are ordered in alphabetical then size-ascending order.
      # Then same amplitudes are ordered in ascending order in target index labels.
      # Then all common indexes (that are summation indexes) are renamed in the ascending order.

      another = self.duplicate()
      
      # reorder amplitudes
      done = 0
      while (not done):
         done = 1
         for namplitudea in range(len(another.amplitudes)):
            for namplitudeb in range(len(another.amplitudes)):
               if (namplitudea >= namplitudeb):
                  continue
               amplitudea = another.amplitudes[namplitudea]
               amplitudeb = another.amplitudes[namplitudeb]
               if (amplitudea.isgreaterthan(amplitudeb,another)):
                  another.swapamplitudes(namplitudea,namplitudeb)
                  done = 0

      # reorder indexes
      for namplitude in range(len(another.amplitudes)):
         amplitude = another.amplitudes[namplitude]
         result = amplitude.canonicalize(another)
         another.amplitudes[namplitude] = copy.deepcopy(result[0])
         parity = result[1]
         another.factor.multiply(parity)

      # relabel summation indexes in the order of appearance
      labelsinuse = []
      for amplitude in another.amplitudes:
         for operator in amplitude.indexes:
            if (not another.summation.hastheindex(operator)):
               labelsinuse.append(operator.index)
      for sequence in another.sequence:
         for operator in sequence:
            if (not another.summation.hastheindex(operator)):
               labelsinuse.append(operator.index)
      oldlabels = []
      newlabels = []
      newlabel = 0
      for amplitude in another.amplitudes:
         for operator in amplitude.indexes:
            if (another.summation.hastheindex(operator)):
               if (operator.index not in oldlabels):
                  oldlabels.append(operator.index)
                  newlabel = newlabel + 1
                  while (newlabel in labelsinuse):
                     newlabel = newlabel + 1
                  newlabels.append(newlabel)
      for operator in another.summation.indexes:
         if (operator.index in oldlabels):
            operator.index = newlabels[oldlabels.index(operator.index)]
      for amplitude in another.amplitudes:
         for operator in amplitude.indexes:
            if (operator.index in oldlabels):
               operator.index = newlabels[oldlabels.index(operator.index)]
      for sequence in another.sequence:
         for operator in sequence:
            if (operator.index in oldlabels):
               operator.index = newlabels[oldlabels.index(operator.index)]

      # reorder summation indexes
      for nindexa in range(len(another.summation.indexes)):
         indexa = another.summation.indexes[nindexa]
         for nindexb in range(len(another.summation.indexes)):
            indexb = another.summation.indexes[nindexb]
            if (nindexa <= nindexb):
               continue
            if (indexa.index < indexb.index):
               swap = another.summation.indexes[nindexa]
               another.summation.indexes[nindexa] = copy.deepcopy(another.summation.indexes[nindexb])
               another.summation.indexes[nindexb] = copy.deepcopy(swap)

      return another

   def isacycliccontraction(self):
      """Returns 1 if self is a cyclic contraction"""
      ncontractions = 0
      for namplitudea in range(len(self.amplitudes)):
         for namplitudeb in range(len(self.amplitudes)):
            if (namplitudea > namplitudeb):
               amplitudea = self.amplitudes[namplitudea]
               amplitudeb = self.amplitudes[namplitudeb]
               for operator in self.summation.indexes:
                  if (operator.isin(amplitudea.indexes) and operator.isin(amplitudeb.indexes)):
                     ncontractions = ncontractions + 1
                     break
      if (ncontractions > len(self.amplitudes) - 1):
         return 1
      else:
         return 0

   def isdisconnected(self,withrespectto=[]):
      """Returns 1 if disconnected; if (withrespectto) connectivity among the given amplitude types is tested"""

      if (self.alreadycontracted()):
 
         # make a connectedness table
         connectedness = [0]*len(self.amplitudes)
         connectedness[0] = 1
         for iteration in range(len(self.amplitudes)):
            for namplitudea in range(len(self.amplitudes)):
               if (connectedness[namplitudea] == 1):
                  amplitudea = self.amplitudes[namplitudea]
                  if ((withrespectto) and (amplitudea.type not in withrespectto)):
                     continue
                  for namplitudeb in range(len(self.amplitudes)):
                     if (connectedness[namplitudeb] == 0):
                        amplitudeb = self.amplitudes[namplitudeb]
                        if ((withrespectto) and (amplitudeb.type not in withrespectto)):
                           continue

                        # see if they have at least one common index
                        exist = 0
                        for indexa in amplitudea.indexes:
                           for indexb in amplitudeb.indexes:
                              if (indexa.isidenticalto(indexb)):
                                 exist = 1
                        if (exist):
                           connectedness[namplitudeb] = 1

         if (0 not in connectedness):
            # connected!
            return 0
         else:
            if (withrespectto):
               for namplitudea in range(len(self.amplitudes)):
                  amplitudea = self.amplitudes[namplitudea]
                  if (amplitudea.type in withrespectto):
                     if (connectedness[namplitudea] == 0):
                        return 1
               return 0
               
            else:
               return 1
 
      else:
         return 0

   def isunlinked(self):
      """Returns 1 if unlinked"""

      if (not self.isdisconnected()):
         return 0

      for seed in range(len(self.amplitudes)): 
         # make a connectedness table
         connectedness = [0]*len(self.amplitudes)
         connectedness[seed] = 1
         for iteration in range(len(self.amplitudes)):
            for namplitudea in range(len(self.amplitudes)):
               if (connectedness[namplitudea] == 1):
                  amplitudea = self.amplitudes[namplitudea]
                  for namplitudeb in range(len(self.amplitudes)):
                     if (connectedness[namplitudeb] == 0):
                        amplitudeb = self.amplitudes[namplitudeb]
 
                        # see if they have at least one common index
                        exist = 0
                        for indexa in amplitudea.indexes:
                           for indexb in amplitudeb.indexes:
                              if (indexa.isidenticalto(indexb)):
                                 exist = 1
                        if (exist):
                           connectedness[namplitudeb] = 1
 
         if (0 in connectedness):
            # disconnected
            closed = 1
            for namplitudeb in range(len(self.amplitudes)):
               amplitudeb = self.amplitudes[namplitudeb]
               if (connectedness[namplitudeb] == 1):
                  for indexa in amplitudeb.indexes:
                     if (not self.summation.hastheindex(indexa)):
                        closed = 0
            if (closed):
               # disconnected & closed = unlinked
               return 1

      return 0
      
   def iszero(self):
      """True if the numerical factor is computationally zero"""
      threshold = 1.0e-12
      zero = 1
      for coefficient in self.factor.coefficients:
         if (abs(coefficient) > threshold):
            zero = 0
      return zero
 
class ListOperatorSequences:

   def __init__(self):
      """Creates a list of operator sequence objects"""
      self.list = []

   def __str__(self):
      """Prints the sequences of operator contractions"""
      print ""
      for line in self.show():
         print line
      return ""

   def show(self):
      """Returns a human-friendly string of the content"""
      show = []
      for operatorsequence in self.list:
         if (operatorsequence == "deleted"):
            show.append("deleted")
         else:
            show.append(operatorsequence.show())
      return show

   def tex(self):
      """Returns a LaTeX string of the content"""
      show = []
      for noperatorsequence in range(len(self.list)):
         operatorsequence = self.list[noperatorsequence]
         if (operatorsequence != "deleted"):
            if (noperatorsequence == 0):
               show.append("\\begin{eqnarray}")
               show.append(string.join(["&&",operatorsequence.tex(),"\\nonumber\\\\"],""))
            elif (noperatorsequence == len(self.list)-1):
               show.append(string.join(["&&",operatorsequence.tex(),"\\nonumber"],""))
               show.append("\\end{eqnarray}")
            else:
               show.append(string.join(["&&",operatorsequence.tex(),"\\nonumber\\\\"],""))
      return show

   def duplicate(self):
      """Makes a copy of itself"""
      duplicate = ListOperatorSequences()
      for operatorsequence in self.list:
         duplicate.list.append(operatorsequence.duplicate())
      return duplicate

   def writetofile(self,filename):
      """Writes the output to a given file"""
      file = open(filename,"w")
      for operatorsequence in self.list:
         file.write(operatorsequence.show())
         file.write("\n")

   def appendtofile(self,filename):
      """Writes the output to a given file"""
      file = open(filename,"a")
      for operatorsequence in self.list:
         file.write(operatorsequence.show())
         file.write("\n")

   def add(self,newoperatorsequence):
      """Adds a new operator sequence member to the list"""
      self.list.append(newoperatorsequence)

   def join(self,another):
      """Joins two list operator sequences"""
      for operatorsequence in another.list:
         self.list.append(operatorsequence)

   def performcontraction(self):
      """Perform a contraction of the left-most operator"""

      # result will be a list of new operator sequence objects
      result = ListOperatorSequences()

      # loop over operator sequences
      for operatorsequence in self.list:

         # call performcontraction()
         result.join(operatorsequence.performcontraction())

      return result

   def simplifythreesub(self,verbose=0):
      """Simplify the list by consolidating operator sequences using permutation of operators"""

      if (len(self.list) == 1):
         return self

      # pick up a pair of operator sequences
      for nsequencea in range(len(self.list)):
#        if (verbose):
#           print 'processing ',nsequencea,' / ',range(len(self.list))
         sequencea = self.list[nsequencea]
         for nsequenceb in range(len(self.list)):
            sequenceb = self.list[nsequenceb]
            if (nsequenceb <= nsequencea):
               continue
            if (sequencea.hasthesameform(sequenceb)):
               sequencec = sequencea.fullyrelabels(sequenceb)
               if (sequencea.canmerge(sequencec)):
                  self.add(sequencea.merges(sequencec))
                  # It is extremely important that the following two statements are executed in this order
                  del self.list[nsequenceb]
                  del self.list[nsequencea]
                  return self
               elif (sequencea.hasnomismatch(sequencec)):
                  permutation = ListOperatorSequences()
                  permutation.add(sequencec)
                  # permutation of tensors
                  for namplitude in range(len(sequencec.amplitudes)):
                     permutation = permutation.amplitudepermutation(namplitude)
                  for sequenced in permutation.list:
                     if (sequencea.canmerge(sequenced)):
                        self.add(sequencea.merges(sequenced))
                        # It is extremely important that the following two statements are executed in this order
                        del self.list[nsequenceb]
                        del self.list[nsequencea]
                        return self
                  # permutation of summation indexes
                  for noperator in range(len(sequencec.summation.indexes)):
                     permutation = permutation.operatorpermutation(noperator)
                  for sequenced in permutation.list:
                     if (sequencea.canmerge(sequenced)):
                        self.add(sequencea.merges(sequenced))
                        # It is extremely important that the following two statements are executed in this order
                        del self.list[nsequenceb]
                        del self.list[nsequencea]
                        return self

      return self

   def simplifyfoursub(self,quick=0):
      """Identify the permutation symmetry among the target indexes"""

      if (len(self.list) == 1):
         return self

      # pick up a pair of operator sequences
      for nsequencea in range(len(self.list)):
         sequencea = self.list[nsequencea]
         for nsequenceb in range(len(self.list)):
            sequenceb = self.list[nsequenceb]
            if (nsequenceb <= nsequencea):
               continue
            if (sequencea.hasthesameform(sequenceb)):
               sequencec = sequencea.fullyrelabels(sequenceb)
               if (sequencea.hasnomismatch(sequencec)):
                  permutation = sequencec.targetindexpermutation()
                  # permutation of target indexes
                  for sequenced in permutation.list:
                     if (sequencea.canmerge(sequenced)):
                        self.add(sequencea.merges(sequenced))
                        # Important that the following two statements are executed in this order
                        del self.list[nsequenceb]
                        del self.list[nsequencea]
                        return self
                  if (quick):
                     continue
                  seed = ListOperatorSequences()
                  seed.add(sequencec)
                  # permutation of tensors
                  for namplitude in range(len(sequencec.amplitudes)):
                     seed = seed.amplitudepermutation(namplitude)
                  for noperator in range(len(sequencec.summation.indexes)):
                     seed = seed.operatorpermutation(noperator)
                  permutation = ListOperatorSequences()
                  for sequenced in seed.list:
                     permutation.join(sequenced.targetindexpermutation())
                  for sequenced in permutation.list:
                     if (sequencea.canmerge(sequenced)):
                        self.add(sequencea.merges(sequenced))
                        # It is extremely important that the following two statements are executed in this order
                        del self.list[nsequenceb]
                        del self.list[nsequencea]
                        return self

      return self

   def simplify(self,verbose=0):
      """Call simplyone through four"""
      #self.simplifyone(1)
      if (self.containscycliccontractions()):
         print " ! Warning! a cyclic contraction is found"
#        self.simplifythree(verbose)
      self.simplifytwo(verbose)
      # the followings do not seem to affect the result, yet it costs enormous memory & time
      # self.simplifyfour(1)
      self = copy.deepcopy(self.deletezero())
      return self

   def simplifyone(self,verbose=0):
      """Consolidate the identical operator sequences"""
      if (len(self.list) == 0):
         return self
      if (verbose):
         print " ... canonicalizing the expressions"
      self = self.canonicalize()
      if (len(self.list) == 1):
         return self
      if (verbose):
         print " ... consolidating terms"
      originallength = len(self.list)
      # pick up a pair of operator sequences
      for nsequencea in range(len(self.list)):
         if (verbose):
            if ((nsequencea/100)*100 == nsequencea):
               print "simplifying:",nsequencea,"/",len(self.list)
         sequencea = self.list[nsequencea]
         if (sequencea == "deleted"):
            continue
         for nsequenceb in range(len(self.list)):
            sequenceb = self.list[nsequenceb]
            if (nsequenceb <= nsequencea):
               continue
            if (sequenceb == "deleted"):
               continue
            if (sequencea.isidenticalto(sequenceb)):
               sequencea.factor.add(sequenceb.factor,1)
               self.list[nsequencea] = copy.deepcopy(sequencea)
               self.list[nsequenceb] = "deleted"
      numberofdeleted = self.list.count("deleted")
      for dummy in range(numberofdeleted):
         self.list.remove("deleted")
      if (verbose):
         print " ... %d terms have been consolidated" %(originallength - len(self.list))
         print " ... number of terms = %d" %(len(self.list))
      return self

   def simplifytwo(self,verbose=0):
      """Identify the permutation symmetries of target indexes"""
      if (len(self.list) == 0):
         return self
      self = self.canonicalize()
      if (len(self.list) == 1):
         return self
      print " ... identifying permutation symmetry among target indexes"
      originallength = len(self.list)
      # pick up a pair of operator sequences
      for nsequencea in range(len(self.list)):
         if (verbose):
            if ((nsequencea/10)*10 == nsequencea):
               print "permutation-simplifying:",nsequencea,"/",len(self.list)
         sequencea = self.list[nsequencea]
         if (sequencea == "deleted"):
            continue
         for nsequenceb in range(len(self.list)):
            if (nsequenceb <= nsequencea):
               continue
            sequenceb = self.list[nsequenceb]
            if (sequenceb == "deleted"):
               continue
            if (sequencea.hasthesameform(sequenceb)):
               permutation = sequenceb.targetindexpermutation()
               permutation = permutation.canonicalize()
               for sequencec in permutation.list:
                  if (sequencea.isidenticalto(sequencec)):
                     sequencea.factor.add(sequencec.factor,1)
                     self.list[nsequencea] = copy.deepcopy(sequencea)
                     self.list[nsequenceb] = "deleted"
                     break
      numberofdeleted = self.list.count("deleted")
      for dummy in range(numberofdeleted):
         self.list.remove("deleted")
      print " ... %d terms have been consolidated" %(originallength - len(self.list))
      print " ... number of terms = %d" %(len(self.list))
      return self

   def simplifythree(self,verbose=0):
      """Aggressively consolidate identical terms"""
      if (len(self.list) == 0):
         return self
      elif (len(self.list) == 1):
         return self
      if (verbose):
         print " ... aggressively consolidating terms"
      done = 0
      iteration = 0
      originallength = len(self.list)
      while (not done):
         iteration = iteration + 1
         if (verbose):
            print 'iteration ',iteration,' number of terms ',len(self.list)
         beforesimplify = len(self.list)
         self = self.simplifythreesub(verbose)
         if (len(self.list) < beforesimplify):
            done = 0
         else:
            done = 1
      print " ... %d terms have been consolidated" %(originallength - len(self.list))
      print " ... number of terms = %d" %(len(self.list))
      return self

   def simplifyfour(self,verbose=0):
      """Aggressively identify the permutation symmetries of target indexes"""
      if (len(self.list) == 0):
         return self
      elif (len(self.list) == 1):
         return self
      if (verbose):
         print " ... aggressively identifying permutation symmetry among target indexes"
      originallength = len(self.list)
      # quick merge to reduce the number of terms
      done = 0
      iteration = 0
      while (not done):
         iteration = iteration + 1
         beforesimplify = len(self.list)
         self = self.simplifyfoursub(1)
         if (len(self.list) < beforesimplify):
            done = 0
         else:
            done = 1
      # more exhaustive merge
      done = 0
      iteration = 0
      while (not done):
         iteration = iteration + 1
         beforesimplify = len(self.list)
         self = self.simplifyfoursub(0)
         if (len(self.list) < beforesimplify):
            done = 0
         else:
            done = 1
      if (originallength - len(self.list) > 0):
         print " ... ***** warning *****"
         print " ... %d terms have been consolidated" %(originallength - len(self.list))
         print " ... number of terms = %d" %(len(self.list))
      return self
      
   def performfullcontraction(self):
      """Performs full contraction of a list of operator sequences and returns a list of tensor contractions"""

      # result will be a list of tensor contractions (operator sequence objects with empty operator sequence)
      self = self.simplifyone()
      result = ListOperatorSequences()

      # loop over operator sequences
      for operatorsequence in self.list:

         # call performcontraction()
         result.join(operatorsequence.performfullcontraction())

      return result

   def operatorpermutation(self,noperatora=0):
      """Return all possible permutation of operators in self"""

      result = ListOperatorSequences()

      for operatorsequence in self.list:
         result.add(operatorsequence)
         if (operatorsequence.summation.indexes):
            operatora = operatorsequence.summation.indexes[noperatora]
            for noperatorb in range(len(operatorsequence.summation.indexes)):
               if (noperatorb <= noperatora):
                  continue
               operatorb = operatorsequence.summation.indexes[noperatorb]
               if (not operatora.issimilarto(operatorb)):
                  continue
               permutation = operatorsequence.duplicate()
               permutation.swapoperators(operatora,operatorb)
               result.add(permutation)

      return result

   def amplitudepermutation(self,namplitudea=0):
      """Return all possible permutation of amplitude in self"""

      result = ListOperatorSequences()

      for operatorsequence in self.list:
         result.add(operatorsequence)
         amplitudea = operatorsequence.amplitudes[namplitudea]
         for namplitudeb in range(len(operatorsequence.amplitudes)):
            if (namplitudeb <= namplitudea):
               continue
            amplitudeb = operatorsequence.amplitudes[namplitudeb]
            if (amplitudea.type != amplitudeb.type):
               continue
            permutation = operatorsequence.duplicate()
            permutation.swapamplitudes(namplitudea,namplitudeb)
            result.add(permutation)

      return result

   def targetsuperpermutation(self,nsupera=0):
      """Return all possible permutation of target super indexes in self"""

      result = ListOperatorSequences()

      for operatorsequence in self.list:

         super = []
         sub = []
         for tensor in operatorsequence.amplitudes:
            for nindex in range(len(tensor.indexes)/2):
               index = tensor.indexes[nindex]
               common = 0
               if (operatorsequence.summation):
                  for another in operatorsequence.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  super.append(tensor.indexes[nindex])
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               index = tensor.indexes[nindex]
               common = 0
               if (operatorsequence.summation):
                  for another in operatorsequence.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  sub.append(tensor.indexes[nindex])

         result.add(operatorsequence)
         supera = super[nsupera]
         for nsuperb in range(len(super)):
            if (nsuperb <= nsupera):
               continue
            superb = super[nsuperb]
            permutation = operatorsequence.duplicate()
            permutation.swapoperators(supera,superb)
            result.add(permutation)

      return result

   def targetsubpermutation(self,nsuba=0):
      """Return all possible permutation of target sub indexes in self"""

      result = ListOperatorSequences()

      for operatorsequence in self.list:

         super = []
         sub = []
         for tensor in operatorsequence.amplitudes:
            for nindex in range(len(tensor.indexes)/2):
               index = tensor.indexes[nindex]
               common = 0
               if (operatorsequence.summation):
                  for another in operatorsequence.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  super.append(tensor.indexes[nindex])
            for nindex in range(len(tensor.indexes)/2,len(tensor.indexes)):
               index = tensor.indexes[nindex]
               common = 0
               if (operatorsequence.summation):
                  for another in operatorsequence.summation.indexes:
                     if (index.isidenticalto(another)):
                        common = 1
               if (not common):
                  sub.append(tensor.indexes[nindex])

         result.add(operatorsequence)
         suba = sub[nsuba]
         for nsubb in range(len(sub)):
            if (nsubb <= nsuba):
               continue
            subb = sub[nsubb]
            permutation = operatorsequence.duplicate()
            permutation.swapoperators(suba,subb)
            result.add(permutation)

      return result

   def canonicalize(self):
      """Reorder amplitudes and common indexes in the canonical order"""

      for noperatorsequence in range(len(self.list)):
         operatorsequence = self.list[noperatorsequence]
         self.list[noperatorsequence] = operatorsequence.canonicalize()
      return self

   def deletedisconnected(self,withrespectto=[]):
      """Deletes disconnected terms"""

      result = ListOperatorSequences()

      originallength = len(self.list)

      # for a fully contracted sequence ...
      for noperatorsequence in range(len(self.list)):
         operatorsequence = self.list[noperatorsequence]
         if (not operatorsequence.isdisconnected(withrespectto)):
            result.add(operatorsequence)

      newlength = len(result.list)

      print " ... %d disconnected terms have been deleted" %(originallength - newlength)

      return result

   def deleteunlinked(self):
      """Deletes unlinked terms"""

      result = ListOperatorSequences()

      originallength = len(self.list)

      # for a fully contracted sequence ...
      for noperatorsequence in range(len(self.list)):
         operatorsequence = self.list[noperatorsequence]
         if (not operatorsequence.isunlinked()):
            result.add(operatorsequence)

      newlength = len(result.list)

      print " ... %d unlinked terms have been deleted" %(originallength - newlength)

      return result

   def containscycliccontractions(self):
      """Returns 1 if self contains a cyclic contraction"""
      for operatorsequence in self.list:
         if (operatorsequence.isacycliccontraction()):
            return 1
      return 0

   def relabelamplitudes(self,old,new):
      """Relabels amplitude"""
      for noperatorsequence in range(len(self.list)):
         operatorsequence = self.list[noperatorsequence]
         for namplitude in range(len(operatorsequence.amplitudes)):
            amplitude = operatorsequence.amplitudes[namplitude]
            if (amplitude.type == old):
               self.list[noperatorsequence].amplitudes[namplitude].type = copy.deepcopy(new)
      return self

   def deletezero(self):
      """Deletes computationally zero terms"""

      result = ListOperatorSequences()

      originallength = len(self.list)

      # for a fully contracted sequence ...
      for noperatorsequence in range(len(self.list)):
         operatorsequence = self.list[noperatorsequence]
         if (not operatorsequence.iszero()):
            result.add(operatorsequence)

      newlength = len(result.list)

      if (originallength != newlength):
         print " !!! WARNING !!! %d computationally zero terms have been deleted" %(originallength - newlength)

      return result
