
from nwchem import *
INT = 1010
CHAR = 1000
LOGICAL = 1011
DBL = 1013
from Tkinter import *
import signal

class Xrtdb:
  def __init__(self):

    # The class is designed to have only one instance,
    # so all of the X root stuff is in here.

    # Ignore SIGCHLD so that ARMCI is not confused by the editor
    signal.signal(signal.SIGCHLD, signal.SIG_DFL)    

    self.root = Tk(baseName = 'Xrtdb', className='xrtdb')

    frame = Frame(self.root)
    frame.pack(expand=YES,fill=BOTH)

    self.scroll=Scrollbar(frame)
    self.listbox = Listbox(frame,bd=2,relief="sunken", 
                                yscrollcommand=(self.scroll.set))
    self.scroll.config(command=self.listbox.yview)

    list = []
    try:
      name = rtdb_first()
      while (name):
        list.append(name)
        name = rtdb_next()
    except NWChemError:
      list.sort()

    i = 0
    for name in list:
      self.listbox.insert(i,name)
      i = i + 1
    
    self.button=Label(frame);
    self.button.quit=Button(self.button, text="Quit", command=frame.quit)
    self.button.new=Button(self.button, text="New", command=self.new)
    self.button.edit=Button(self.button, text="Edit", command=self.edit)
    self.button.delete=Button(self.button, text="Delete", command=self.delete)
    self.button.help=Button(self.button, text="Help",  command=self.help)

    self.text=Text(frame,bd=2,relief="sunken",height=1,width=35)

    self.listbox.bind("<ButtonRelease-1>", func=self.select) 
    self.listbox.bind("<Double-ButtonRelease-1>",func=self.selectdouble)

    self.button.pack(fill=X)
    self.text.pack(fill=X)
    self.scroll.pack(side=LEFT,fill=Y)
    self.listbox.pack(fill=BOTH,expand=1)
    self.button.quit.pack(side=LEFT,padx=5)
    self.button.delete.pack(side=LEFT,padx=5)
    self.button.edit.pack(side=LEFT,padx=5)
    self.button.new.pack(side=LEFT,padx=5)
    self.button.help.pack(side=LEFT,padx=5)

    #wm title . "$filename"
    #wm iconname . "Xrtdb"

    self.root.mainloop()

  def new(self):
    from SimpleDialog import *
    import tkSimpleDialog

    self.selection_clear()

    dialog = SimpleDialog(self.root,text="Select data type",
                          buttons=["int","dbl","logical","char","cancel"])
    status = dialog.go()
    if (status == 0):
      ma_type = INT
    elif (status == 1):
      ma_type = DBL
    elif (status == 2):
      ma_type = LOGICAL
    elif (status == 3):
      ma_type = 1000
    else:
      self.set_text("")
      return

    name = tkSimpleDialog.askstring("xrtdb","Enter name")
    if (not name):
      self.set_text("")
      return

    self.set_text("Edit %s of type %s" % (name,self.ma_type_name(ma_type)))
    self.delete_tmpfile() 
    self.write_tmpfile("Replace this file with your data of type %s" % 
                        self.ma_type_name(ma_type), 1000)
    self.process_tmpfile(name, ma_type, 1)

  def edit(self):
    from SimpleDialog import *

    try:
      name = self.selection_value()
    except:
      self.set_text("Select an entry then edit")
      return

    self.set_text("Editing %s" % name)
    try:
      values = rtdb_get(name)
      (ma_type, nelem, date) = rtdb_get_info(name)
    except NWChemError:
      self.set_text("Failed to get %s" % name)
      self.selection_clear()
      return
    
    try:
      self.write_tmpfile(values, ma_type)
    except:
      self.set_text("Failed writing file for edit")
      self.selection_clear()
      return
      
    self.process_tmpfile(name, ma_type, 0)

  def delete(self):
    try:
      name = self.selection_value()
    except:
      self.set_text("Select an entry then delete")
      return

    try:
      rtdb_delete(name)
      self.set_text("Deleted %s" % name)
      self.listbox.delete(self.selection_index())
    except NWChemError:
      self.set_text("Failed to delete %s" % name)
    self.selection_clear()

  def help(self):
    from SimpleDialog import *
    dialog = SimpleDialog(self.root,
                          text="Quit \t- immediately exits.\n"
                               "Delete \t- deletes selected entry.\n"
                               "Edit \t- edits selected entry, optional save.\n"
                               "New \t- makes new entry, prompts for type+name.\n"
                               "Help \t- this is all that there is.\n"
                               "\n"
                               "Single click on item to select it.\n"
                               "Double click on item to edit it.\n"
                               "\n"
                               "When editing\n"
                               "\t- enter one value per line.\n"
                               "\t- use true/false for logicals.\n"
                               "\t- use 'e' in floating point exponents.\n"
                               "\t- character strings are terminated at EOL.",
                          buttons=["Done"])
    dialog.go()    
  
  def select(self,event):
    name = self.selection_value()
    try:
      (ma_type, nelem, date) = rtdb_get_info(name)
      type = self.ma_type_name(ma_type)
      self.set_text("%s = %s(%d) %s" % (name, type, nelem, date))
    except NWChemError:
      self.set_text("Failed to get info for %s" % name)
      self.selection_clear()

  def selectdouble(self,event):
    self.select(event)
    self.edit()

  def selection_clear(self):
    self.listbox.select_clear(0,"end")

  def selection_index(self):
    item = self.listbox.curselection()
    return int(item[0])

  def selection_value(self):
    return self.listbox.get(self.selection_index())
    

  def set_text(self,string):
    self.text.delete("0.0","end")
    self.text.insert("0.0",string)

  def edit_tmpfile(self):
    import os
    if (os.system('emacs xrtdbtmp.txt')):
       raise "Edit failed"

  def delete_tmpfile(self):
    import os
    try:
      os.remove('xrtdbtmp.txt')
    except:
      pass

  def write_tmpfile(self, values, ma_type):
    import types
    tmpfile = open('xrtdbtmp.txt','w+')
    if (type(values) == type([])):
      for value in values:
        tmpfile.write(self.value_to_string(value,ma_type)+'\n')
    else:
       tmpfile.write(self.value_to_string(values,ma_type)+'\n')
    tmpfile.close()

  def read_tmpfile(self, ma_type):
    import string
    result = []
    tmpfile = open('xrtdbtmp.txt','r')

    line = tmpfile.readline()
    while (line):
      if (string.split(line)):
        for value in string.split(line):
           if (value):
             try:
               result.append(self.string_to_value(string.strip(value),ma_type))
             except:
               typename = self.ma_type_name(ma_type)
               self.set_text("'%s' is not a valid %s" % (value,typename))
               self.selection_clear()
               tmpfile.close()
               raise "Invalid type"
      else:
        result.append(" ")
      line = tmpfile.readline()
    tmpfile.close()
    return result

  def process_tmpfile(self, name, ma_type, isnew):
    from SimpleDialog import *
    # Edit tmpfile
    # Prompt for saving changes
    # Try to read the results
    # Store in the database
    # If it is a new entry update the listbox

    try:
      self.edit_tmpfile()
    except "Edit failed":
      self.set_text("Edit failed")
      self.selection_clear()

    dialog = SimpleDialog(self.root,text="Save edits?",buttons=["Yes","No"])
    status = dialog.go()    

    if (status == 0):
      try:
        result = self.read_tmpfile(ma_type)
      except "Invalid type":
        return  # msg set in read_tmpfile

      try:
        rtdb_put(name,result,ma_type)
        self.set_text("%s edits saved" % name)
      except NWChemError:
        self.set_text("%s save failed!" % name)
      if (isnew):
	n = 0
	try:
	  while (name > self.listbox.get(n)):
	    n = n + 1
	except:
	  pass
	self.listbox.insert(n,name)
	self.listbox.see(n)
    else:
      self.set_text("%s nothing saved" % name)
      
    self.selection_clear()

  def value_to_string(self,value,ma_type):
    if (ma_type == INT):
      return "%d" % value
    elif (ma_type == DBL):
      return "%21.15e" % value
    elif (ma_type == LOGICAL):
      if (value):
        return "true"
      else:
        return "false"
    elif (ma_type == 1000):  # since Tk overwrites defn of CHAR
      return value
    else:
      return value

  def string_to_value(self,value,ma_type):
    import string
    if (ma_type == INT):
      return string.atoi(value)
    elif (ma_type == DBL):
      return string.atof(value)
    elif (ma_type == LOGICAL):
      if (value == "true"):
        return 1
      elif (value == "false"):
        return 0
      else:
        raise NWChemError,'invalid value'
    elif (ma_type == 1000):  # since Tk overwrites defn of CHAR
      return value
    else:
      raise NWChemError,'invalid type'

  def ma_type_name(self, ma_type):
    if (ma_type == INT):
      return "int"
    elif (ma_type == DBL):
      return "double"
    elif (ma_type == LOGICAL):
      return "logical"
    elif (ma_type == 1000):  # since Tk overwrites defn of CHAR
      return "char"
    else:
      raise NWChemError,'invalid type'
