# Ground-state coupled-cluster ansatz generator for OCE
# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)

from Tkinter import *
import tkMessageBox
import tkSimpleDialog
import tkFileDialog

import string
import copy
import sys
import oce
import tce

class Window:

   def __init__(self,master):
   
      self.frame = Frame(master)
      self.frame.pack(fill=X)
      self.progress = 0
      
      Label(self.frame, text = "Tensor Contraction Engine, Version 1.0").grid(row=0,columnspan=6)
      Label(self.frame, text = "Copyright (c) 2003, Battelle & Pacific Northwest National Laboratory").grid(row=1,columnspan=6)
      self.l = IntVar()
      self.cl0 = Radiobutton(self.frame,text = "<0|", variable=self.l,value=0)
      self.cl0.grid(row=2,column=0,sticky=W)
      self.cl1 = Radiobutton(self.frame,text = "<S|", variable=self.l,value=1)
      self.cl1.grid(row=3,column=0,sticky=W)
      self.cl2 = Radiobutton(self.frame,text = "<D|", variable=self.l,value=2)
      self.cl2.grid(row=4,column=0,sticky=W)
      self.cl3 = Radiobutton(self.frame,text = "<T|", variable=self.l,value=3)
      self.cl3.grid(row=5,column=0,sticky=W)
      self.cl4 = Radiobutton(self.frame,text = "<Q|", variable=self.l,value=4)
      self.cl4.grid(row=6,column=0,sticky=W)
      self.y0 = IntVar()
      self.cy0 = Checkbutton(self.frame,text = "L0 = 1", variable=self.y0,onvalue=1,offvalue=0)
      self.cy0.grid(row=2,column=1,sticky=W)
      self.cy0.select()
      self.y1 = IntVar()
      self.cy1 = Checkbutton(self.frame,text = "L1 Operator", variable=self.y1,onvalue=1,offvalue=0)
      self.cy1.grid(row=3,column=1,sticky=W)
      self.y2 = IntVar()
      self.cy2 = Checkbutton(self.frame,text = "L2 Operator", variable=self.y2,onvalue=2,offvalue=0)
      self.cy2.grid(row=4,column=1,sticky=W)
      self.y3 = IntVar()
      self.cy3 = Checkbutton(self.frame,text = "L3 Operator", variable=self.y3,onvalue=3,offvalue=0)
      self.cy3.grid(row=5,column=1,sticky=W)
      self.y4 = IntVar()
      self.cy4 = Checkbutton(self.frame,text = "L4 Operator", variable=self.y4,onvalue=4,offvalue=0)
      self.cy4.grid(row=6,column=1,sticky=W)
      self.h1 = IntVar()
      self.ch1 = Checkbutton(self.frame,text = "<p|f|q>{p+q}", variable=self.h1,onvalue=1,offvalue=0)
      self.ch1.grid(row=3,column=2,sticky=W)
      self.ch1.select()
      self.h2 = IntVar()
      self.ch2 = Checkbutton(self.frame,text = "1/4<pq||rs>{p+q+sr}", variable=self.h2,onvalue=2,offvalue=0)
      self.ch2.grid(row=4,column=2,sticky=W)
      self.ch2.select()
      self.t1 = IntVar()
      self.c1 = Checkbutton(self.frame,text = "T1 Operator", variable=self.t1,onvalue=1,offvalue=0)
      self.c1.grid(row=3,column=3,sticky=W)
      self.t2 = IntVar()
      self.c2 = Checkbutton(self.frame,text = "T2 Operator", variable=self.t2,onvalue=2,offvalue=0)
      self.c2.grid(row=4,column=3,sticky=W)
      self.t3 = IntVar()
      self.c3 = Checkbutton(self.frame,text = "T3 Operator", variable=self.t3,onvalue=3,offvalue=0)
      self.c3.grid(row=5,column=3,sticky=W)
      self.t4 = IntVar()
      self.c4 = Checkbutton(self.frame,text = "T4 Operator", variable=self.t4,onvalue=4,offvalue=0)
      self.c4.grid(row=6,column=3,sticky=W)
      self.x0 = IntVar()
      self.cx0 = Checkbutton(self.frame,text = "R0 = 1", variable=self.x0,onvalue=1,offvalue=0)
      self.cx0.grid(row=2,column=4,sticky=W)
      self.cx0.select()
      self.x1 = IntVar()
      self.cx1 = Checkbutton(self.frame,text = "R1 Operator", variable=self.x1,onvalue=1,offvalue=0)
      self.cx1.grid(row=3,column=4,sticky=W)
      self.x2 = IntVar()
      self.cx2 = Checkbutton(self.frame,text = "R2 Operator", variable=self.x2,onvalue=2,offvalue=0)
      self.cx2.grid(row=4,column=4,sticky=W)
      self.x3 = IntVar()
      self.cx3 = Checkbutton(self.frame,text = "R3 Operator", variable=self.x3,onvalue=3,offvalue=0)
      self.cx3.grid(row=5,column=4,sticky=W)
      self.x4 = IntVar()
      self.cx4 = Checkbutton(self.frame,text = "R4 Operator", variable=self.x4,onvalue=4,offvalue=0)
      self.cx4.grid(row=6,column=4,sticky=W)
      self.r = IntVar()
      self.cr0 = Radiobutton(self.frame,text = "|0>", variable=self.r,value=0)
      self.cr0.grid(row=2,column=5,sticky=W)
      self.cr1 = Radiobutton(self.frame,text = "|S>", variable=self.r,value=1)
      self.cr1.grid(row=3,column=5,sticky=W)
      self.cr2 = Radiobutton(self.frame,text = "|D>", variable=self.r,value=2)
      self.cr2.grid(row=4,column=5,sticky=W)
      self.cr3 = Radiobutton(self.frame,text = "|T>", variable=self.r,value=3)
      self.cr3.grid(row=5,column=5,sticky=W)
      self.cr4 = Radiobutton(self.frame,text = "|Q>", variable=self.r,value=4)
      self.cr4.grid(row=6,column=5,sticky=W)
      self.ccv1 = IntVar()
      self.cc1 = Checkbutton(self.frame,text = "L Is Connected", variable=self.ccv1,onvalue=1,offvalue=0)
      self.cc1.grid(row=7,column=1,sticky=W)
      self.cc1.select()
      self.ccv2 = IntVar()
      self.cc2 = Checkbutton(self.frame,text = "H Is Connected", variable=self.ccv2,onvalue=1,offvalue=0)
      self.cc2.grid(row=7,column=2,sticky=W)
      self.cc2.select()
      self.ccv3 = IntVar()
      self.cc3 = Checkbutton(self.frame,text = "T Is Connected", variable=self.ccv3,onvalue=1,offvalue=0)
      self.cc3.grid(row=7,column=3,sticky=W)
      self.cc3.select()
      self.ccv4 = IntVar()
      self.cc4 = Checkbutton(self.frame,text = "R Is Connected", variable=self.ccv4,onvalue=1,offvalue=0)
      self.cc4.grid(row=7,column=4,sticky=W)
      self.cc4.select()
      self.ccv5 = IntVar()
      self.cc5 = Checkbutton(self.frame,text = "All Are Linked", variable=self.ccv5,onvalue=1,offvalue=0)
      self.cc5.grid(row=7,column=5,sticky=W)
      self.cc5.select()
      self.buttonOK = Button(self.frame, text = "Generate Ansatz", command = self.go, width=100).grid(row=8,columnspan=6)
      self.buttonSKIP = Button(self.frame, text = "Skip", command = self.skip, width=20).grid(row=9,column=0,columnspan=3)
      self.buttonCA = Button(self.frame, text = "Clear All", command = self.clearall, width=20).grid(row=9,column=3,columnspan=3)

      self.frame2 = Frame(master)
      self.frame2.pack(fill=X)
      self.scrollbar = Scrollbar(self.frame2)
      self.scrollbar.pack(side=RIGHT, fill=Y)
      self.listbox = Listbox(self.frame2, yscrollcommand=self.scrollbar.set)
      self.scrollbar.config(command=self.listbox.yview)
      self.oce_result = ""
      self.tce_result = ""

      self.frame3 = Frame(master)
      self.frame3.pack(fill=X)
      self.status = Label(self.frame3, text = "", bd=1, relief=SUNKEN, anchor=W)
      self.status.pack(side=BOTTOM, fill=X)
      
      menu = Menu(master)
      master.config(menu=menu)
      filemenu = Menu(master)
      menu.add_cascade(label="File",menu=filemenu)
      filemenu.add_command(label="Import second quantized expressions",command=self.read_oceinput)
      filemenu.add_command(label="Import tensor contraction expressions",command=self.read_tceinput)
      filemenu.add_command(label="Exit",command=self.exit)
      helpmenu = Menu(menu)
      menu.add_cascade(label="Help",menu=helpmenu)
      helpmenu.add_command(label="About",command=self.about)
      helpmenu.add_command(label="References",command=self.references)

   def read_oceinput(self):
      
      filename = tkFileDialog.askopenfilename()
      if (filename != ""):
         self.oce_result = oce.readfromfile(filename)
         self.listbox.delete(0,END)
         for line in self.oce_result.show():
             self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         self.status.config(text = repr(len(self.oce_result.show()))+" second quantized expressions imported")
         self.status.update_idletasks()
      self.progress = 1
      self.buttonOK = Button(self.frame, text = "Perform Operator Contractions", command = self.go, width = 100).grid(row=8,columnspan=6)
      return

   def read_tceinput(self):
      
      filename = tkFileDialog.askopenfilename()
      if (filename != ""):
         self.tce_result = tce.readfromfile(filename)
         self.listbox.delete(0,END)
         for line in self.tce_result.show():
             self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         self.status.config(text = repr(len(self.tce_result.show()))+" multiple tensor contraction expressions imported")
         self.status.update_idletasks()
      self.progress = 61
      self.buttonOK = Button(self.frame, text = "Perform Strength Reduction", command = self.go, width = 100).grid(row=8,columnspan=6)
      return

   def exit(self):
      sys.exit()

   def about(self):
      tkMessageBox.showinfo("About", "Tensor Contraction Engine, Version 1.0\nCopyright (c) 2003, Battelle & Pacific Northwest National Laboratory\n")

   def references(self):
      tkMessageBox.showinfo("References", "S. Hirata, J. Phys. Chem. A 107, 9887-9897 (2003).\n")

   def go(self):

      if (self.progress == 0):
         horders = []
         if (self.h1.get() != 0):
            horders.append(self.h1.get())
         if (self.h2.get() != 0):
            horders.append(self.h2.get())
         if (not horders):
	    horders.append(0)
         torders = [0]
         if (self.t1.get() != 0):
            torders.append(self.t1.get())
         if (self.t2.get() != 0):
            torders.append(self.t2.get())
         if (self.t3.get() != 0):
            torders.append(self.t3.get())
         if (self.t4.get() != 0):
            torders.append(self.t4.get())
         xorders = [0]
         if (self.x1.get() != 0):
            xorders.append(self.x1.get())
         if (self.x2.get() != 0):
            xorders.append(self.x2.get())
         if (self.x3.get() != 0):
            xorders.append(self.x3.get())
         if (self.x4.get() != 0):
            xorders.append(self.x4.get())
         yorders = [0]
         if (self.y1.get() != 0):
            yorders.append(self.y1.get())
         if (self.y2.get() != 0):
            yorders.append(self.y2.get())
         if (self.y3.get() != 0):
            yorders.append(self.y3.get())
         if (self.y4.get() != 0):
            yorders.append(self.y4.get())
         r0is1 = self.x0.get()
         l0is1 = self.y0.get()
         leftprojection = self.l.get()
         rightprojection = self.r.get()
         
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
         
         self.status.config(text = ansatz)
         self.status.update_idletasks()
   
         listofansatz = []
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
                              elif (hamiltonian != 0):
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
                              listofansatz.append(ansatz)
         
         self.oce_result = oce.ListOperatorSequences()
         for line in listofansatz:
            self.oce_result.add(oce.stringtooperatorsequence(line))
         for line in self.oce_result.show():
            self.listbox.insert(END,line)
         self.listbox.pack(fill=X)

         self.progress = 1
         self.buttonOK = Button(self.frame, text = "Perform Operator Contractions", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 1):

         self.oce_result = self.oce_result.performfullcontraction()
         self.listbox.delete(0,END)
         for line in self.oce_result.show():
             self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         self.status.config(text = repr(len(self.oce_result.show()))+" multiple tensor contractions")
         self.status.update_idletasks()
         self.progress = 2
         self.buttonOK = Button(self.frame, text = "Delete Disconnected Diagrams", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 2):

         originallen = len(self.oce_result.show())
         connectedtensors = []
         if (self.ccv1.get() == 1):
            connectedtensors.append('y')
         if (self.ccv2.get() == 1):
            connectedtensors.append('f')
            connectedtensors.append('v')
         if (self.ccv3.get() == 1):
            connectedtensors.append('t')
         if (self.ccv4.get() == 1):
            connectedtensors.append('x')
         if (connectedtensors):
            self.oce_result = self.oce_result.deletedisconnected(connectedtensors)
            self.listbox.delete(0,END)
            for line in self.oce_result.show():
               self.listbox.insert(END, line)
            self.listbox.pack(fill=X)
            newlen = len(self.oce_result.show())
            self.status.config(text = repr(newlen)+" multiple tensor contractions ("+repr(originallen-newlen)+" disconnected diagrams deleted)")
            self.status.update_idletasks()
         self.progress = 3
         self.buttonOK = Button(self.frame, text = "Delete Unlinked Diagrams", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 3):

         originallen = len(self.oce_result.show())
         if (self.ccv5.get() == 1):
            linked = 1
         else:
            linked = 0
         if (linked):
            originallen = len(self.oce_result.show())
            self.oce_result = self.oce_result.deleteunlinked()
            self.listbox.delete(0,END)
            for line in self.oce_result.show():
               self.listbox.insert(END, line)
            self.listbox.pack(fill=X)
            newlen = len(self.oce_result.show())
            self.status.config(text = repr(newlen)+" multiple tensor contractions ("+repr(originallen-newlen)+" unlinked diagrams deleted)")
            self.status.update_idletasks()
         self.progress = 4
         self.buttonOK = Button(self.frame, text = "Extract Permutation Symmetry", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 4):
         originallen = len(self.oce_result.show())
         self.oce_result = self.oce_result.simplify()
         self.listbox.delete(0,END)
         for line in self.oce_result.show():
            self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         newlen = len(self.oce_result.show())
         self.status.config(text = repr(newlen)+" multiple tensor contractions ("+repr(originallen-newlen)+" diagrams merged)")
         if tkMessageBox.askyesno("TCE", "Relabel operators?\n"):
            oldname = tkSimpleDialog.askstring("TCE", "Old label for operator?\n")
            newname = tkSimpleDialog.askstring("TCE", "New label for operator?\n")
            while (oldname and newname):
               self.oce_result.relabelamplitudes(oldname,newname)
               self.listbox.delete(0,END)
               for line in self.oce_result.show():
                  self.listbox.insert(END,line)
               self.listbox.pack(fill=X)
               oldname = tkSimpleDialog.askstring("TCE", "Old label for operator?\n")
               newname = tkSimpleDialog.askstring("TCE", "New label for operator?\n")
         self.status.update_idletasks()
         self.progress = 5
         self.buttonOK = Button(self.frame, text = "Save Working Equations...", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 5):
         if tkMessageBox.askyesno("TCE", "Save working equations?\n"):
            filename = tkFileDialog.asksaveasfilename()
            if (filename != ""):
               self.oce_result.writetofile(filename)
               self.status.config(text = "Working equation saved")
               self.status.update_idletasks()
         self.progress = 6
         self.buttonOK = Button(self.frame, text = "Perform Strength Reduction", command = self.go, width = 100).grid(row=8,columnspan=6)
         return
            
      if (self.progress == 6):
         self.tce_result = tce.ListTensorContractions()
         for selfexpr in self.oce_result.list:
            selfexpr = string.split(selfexpr.show())
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
                        newindex = tce.Index(type,label)
                        indexes.append(newindex)
                  permutations.append(indexes)
                  pointer = pointer + 2
                  if (selfexpr[pointer] == "]"):
                     donewithfactors = 1
               else:
                  permutations.append([])
            factor = tce.Factor(coefficients,permutations)
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
                     newindex = tce.Index(type,label)
                     indexes.append(newindex)
               summation = tce.Summation(indexes)
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
                                 indexes.append(tce.Index(type,label))
                  tensors.append(tce.Tensor(tensortype,indexes,tensorlabel,tensorconj))
            self.tce_result.list.append(tce.TensorContraction(factor,summation,tensors))
            
         self.tce_result = self.tce_result.breakdown()
         self.listbox.delete(0,END)
         for line in self.tce_result.show():
            self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         self.status.config(text = repr(len(self.tce_result.children))+" first order binary contractions")
         self.status.update_idletasks()
         self.progress = 7
         self.buttonOK = Button(self.frame, text = "Perform Factorization", command = self.go, width = 100).grid(row=8,columnspan=6)
         return
            
      if (self.progress == 61):
         self.tce_result = self.tce_result.breakdown()
         self.listbox.delete(0,END)
         for line in self.tce_result.show():
            self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         self.status.config(text = repr(len(self.tce_result.children))+" first order binary contractions")
         self.status.update_idletasks()
         self.progress = 7
         self.buttonOK = Button(self.frame, text = "Perform Factorization", command = self.go, width = 100).grid(row=8,columnspan=6)
         return
         
      if (self.progress == 7):
         originallen = len(self.tce_result.children)
         self.tce_result = self.tce_result.fullyfactorize()
         self.listbox.delete(0,END)
         for line in self.tce_result.show():
            self.listbox.insert(END, line)
         self.listbox.pack(fill=X)
         newlen = len(self.tce_result.children)
         self.status.config(text = repr(newlen)+" first order binary contractions ("+repr(originallen-newlen)+" first order contractions merged)")
         self.status.update_idletasks()
         self.progress = 8
         self.buttonOK = Button(self.frame, text = "Generate Fortran77 Code", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 8):
         types = self.tce_result.tensortypes()
         message = ""
         excitation = ""
         excitationold = ""
         for type in types[1]:
            excitationold = excitationold + type + " "
            if (excitation):
               excitation = excitation + "and '" + type + "' "
            else:
               excitation = "'" + type + "' "
         if (excitation):
            message = message + excitation + "are excitation tensors\n"
         deexcitation = ""
         deexcitationold = ""
         for type in types[2]:
            deexcitationold = deexcitationold + type + " "
            if (deexcitation):
               deexcitation = deexcitation + "and '" + type + "' "
            else:
               deexcitation = "'" + type + "' "
         if (deexcitation):
            message = message + deexcitation + "are deexcitation tensors\n"
         intermediate = ""
         intermediateold = ""
         for type in types[3]:
            intermediateold = intermediateold + type + " "
            if (intermediate):
               intermediate = intermediate + "and '" + type + "' "
            else:
               intermediate = "'" + type + "' "
         if (intermediate):
            message = message + intermediate + "are intermediate tensors\n"
         general = ""
         generalold = ""
         for type in types[4]:
            generalold = generalold + type + " "
            if (general):
               general = general + "and '" + type + "' "
            else:
               general = "'" + type + "' "
         if (general):
            message = message + general + "are general tensors\n"
         if tkMessageBox.askyesno("TCE", message):
            subroutinename = tkSimpleDialog.askstring("TCE", "Stub Name of Fortran77 code?\n")
            if (subroutinename):
               self.tce_result = self.tce_result.fortran77(subroutinename)
               self.listbox.delete(0,END)
               for line in self.tce_result.show():
                  self.listbox.insert(END, line)
               self.listbox.pack(fill=X)
               self.status.config(text = repr(len(self.tce_result.show()))+" lines of Fortran code generated")
               self.status.update_idletasks()
               self.progress = 9
               self.buttonOK = Button(self.frame, text = "Save Fortran77 Code...", command = self.go, width = 100).grid(row=8,columnspan=6)
         else:
            excitationlist = []
            deexcitationlist = []
            intermediatelist = []
            generallist = []
            dialogresult = tkSimpleDialog.askstring("TCE","Excitation Tensors",initialvalue=excitationold)
            if (dialogresult):
               excitationlist = string.split(dialogresult)
            dialogresult = tkSimpleDialog.askstring("TCE","Deexcitation Tensors",initialvalue=deexcitationold)
            if (dialogresult):
               deexcitationlist = string.split(dialogresult)
            dialogresult = tkSimpleDialog.askstring("TCE","Intermediate Tensors",initialvalue=intermediateold)
            if (dialogresult):
               intermediatelist = string.split(dialogresult)
            dialogresult = tkSimpleDialog.askstring("TCE","General Tensors",initialvalue=generalold)
            if (dialogresult):
               generallist = string.split(dialogresult)
            subroutinename = tkSimpleDialog.askstring("TCE", "Stub Name of Fortran77 code?\n")
            if (subroutinename):
               self.tce_result = self.tce_result.fortran77(subroutinename,excitationlist,deexcitationlist,intermediatelist,generallist)
               self.listbox.delete(0,END)
               for line in self.tce_result.show():
                  self.listbox.insert(END, line)
               self.listbox.pack(fill=X)
               self.status.config(text = repr(len(self.tce_result.show()))+" lines of Fortran code generated")
               self.status.update_idletasks()
               self.progress = 9
               self.buttonOK = Button(self.frame, text = "Save Fortran77 Code...", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 9):
         if tkMessageBox.askyesno("TCE", "Save Fortran77 code?\n"):
            filename = tkFileDialog.asksaveasfilename()
            if (filename != ""):
               self.tce_result.writetofile(filename)
         self.status.config(text = "")
         self.status.update_idletasks()
         self.progress = 0
         self.buttonOK = Button(self.frame, text = "Generate Ansatz", command = self.go, width = 100).grid(row=8,columnspan=6)
         self.listbox.delete(0,END)
         return

   def skip(self):
   
      if (self.progress == 0):
         tkMessageBox.showwarning("TCE", "Unable to skip")
         return

      if (self.progress == 1):
         tkMessageBox.showwarning("TCE", "Unable to skip")
         return

      if (self.progress == 2):
         self.progress = 3
         self.buttonOK = Button(self.frame, text = "Delete Unlinked Diagrams", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 3):
         self.progress = 4
         self.buttonOK = Button(self.frame, text = "Extract Permutation Symmetry", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 4):
         tkMessageBox.showwarning("TCE", "Unable to skip")
         return

      if (self.progress == 5):
         self.progress = 6
         self.buttonOK = Button(self.frame, text = "Perform Strength Reduction", command = self.go, width = 100).grid(row=8,columnspan=6)
         return
            
      if (self.progress == 6):
         tkMessageBox.showwarning("TCE", "Unable to skip")
         return
         
      if (self.progress == 7):
         self.progress = 8
         self.buttonOK = Button(self.frame, text = "Generate Fortran77 Code", command = self.go, width = 100).grid(row=8,columnspan=6)
         return

      if (self.progress == 8):
         tkMessageBox.showwarning("TCE", "Unable to skip")
         return

      if (self.progress == 9):
         self.progress = 0
         self.buttonOK = Button(self.frame, text = "Generate Ansatz", command = self.go, width = 100).grid(row=8,columnspan=6)
         self.listbox.delete(0,END)
         return

   def clearall(self):
   
      self.cl0.select()
      self.cr0.select()
      self.cy0.select()
      self.cy1.deselect()
      self.cy2.deselect()
      self.cy3.deselect()
      self.cy4.deselect()
      self.c1.deselect()
      self.c2.deselect()
      self.c3.deselect()
      self.c4.deselect()
      self.ch1.select()
      self.ch2.select()
      self.cx0.select()
      self.cx1.deselect()
      self.cx2.deselect()
      self.cx3.deselect()
      self.cx4.deselect()
      self.cc1.select()
      self.cc2.select()
      self.cc3.select()
      self.cc4.select()
      self.cc5.select()
      self.listbox.delete(0,END)
      self.progress = 0
      self.buttonOK = Button(self.frame, text = "Generate Ansatz", command = self.go, width = 100).grid(row=8,columnspan=6)
      self.status.config(text = "")
      self.status.update_idletasks()

root = Tk()

window = Window(root)

root.mainloop()
