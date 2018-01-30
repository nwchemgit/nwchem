# Since each sub-task assumes it can read parameters from stdin,
# we have to call them like so:
#
# Load up the nwchem library.
# db = nwchem_init("test", "1500 mb")
#
# Set parameters and do a DFT run.
# dft(db, "xc b3lyp")
# nwtask(db, "energy")

import atexit
from nwproto import *
from rtdb import *

# NWChem Initialization.
# Send it the output prefix name, e.g. "h2o"
# and receive an rtdb class object.
def nwchem_init(name, mem="total 400 mb", scratch="./scratch", perm="./perm"):
	def to_strbuf(x):
		s = strbuf()
		s[:len(x)] = x
		s[len(x)] = chr(0)
		return s
	inp = map(lambda x: byref(to_strbuf(x)), [name, mem, scratch, perm])
	db = nwlib.nwchem_init_(*inp) # calls pbegin, etc.
	z = byref(c_int(0))
	atexit.register(nwlib.nwchem_dtor_, byref(c_int(db)), z, z) # "swept and put in order."
	return rtdb(db)

# cast the function to run only on the head node.
def single_op(f):
	def run_single(*args):
	    if nwlib.ga_nodeid_() == 0:
		nwlib.rtdb_parallel(0)
		f(*args)
		nwlib.rtdb_parallel(1)
	return run_single

# Abstract wrapper for an interface function.
def nw_interface(fname):
	try:
		nwf = getattr(nwlib, fname+"_") # fortran name mangling
	except AttributeError:
		print "Symbol not found: %s_"%fname
		return lambda x,y: None
	int_fn(nwf, c_int_p) # declare the function to ctypes.
	def fn(db, inp):
		nwlib.push_inp_cstring(inp)
		nwf(byref(db._n))
	return single_op(fn)

# Wrap all top-level functions.
toplev = ["geom", "bsse", "bas",
#	"python",
	"cosmo", "intgrl", "scf",
	"mp2", "drdy", "stepper",
	"mepgs", "tropt", "driver",
#	"string",
	"dft", "occup",
	"pre", "md", "argos",
	"esp", "et", "transp",
	"ana", "dia",
	"gradients", "ccsd", "oniom",
	"mcscf", "plnwv", "bq",
	"cons", "dplot", "prop",
	"speech", "nwpw", "smd",
	"rism", "qmmm", "ccca",
	"rel", "nbo", "vscf",
	"raman", "dntmc", "freq_vib",
	"hess", "tddft", "mymd", "mymc" ]

for fname in toplev:
	globals()[fname] = nw_interface(fname+"_input")

# Some need re-naming...
# First is the name we shall call it by.
top_cust = [ ("string_input", "string_input"),
		("time_input", "input_time"),
		("set_input", "input_set"),
		("unset_input", "input_unset"),
		("title_input", "input_title"),
		("qcharge", "input_qcharge"),
		("charge", "input_charge") ]
for fname, alt_name in top_cust:
	globals()[fname] = nw_interface(alt_name)

# The awkward squad:
#            else if (inp_compare(.false.,test,'print')) then
#               call util_print_input(db,' ')
#               call util_print_rtdb_load(db,' ') ! High level print
#            else if (inp_compare(.false.,test,'noprint')) then
#               call util_print_input(db,' ')

# tce -> tce(db, inp, 'tce')
# uccsd -> tce(db, inp, 'uccsd')
# etc.
def do_tce(db, inp, *tp):
	if 'tce' in tp: # sic
		db['tce:model'] = 'ccsd'
	for t in tp:
		db['tce:module'] = tp
	nwlib.push_inp_cstring(inp)
	nwlib.tce_input_(byref(db._n))
tce = single_op(do_tce) # wrap with single_op incantations

def nwtask(db, task):
	def read_task():
		nwlib.push_inp_cstring("task " + task)
		nwlib.task_input_(byref(db._n))
		nwlib.util_print_rtdb_load_(byref(db._n), "", 0)
	single_op(read_task)()
	nwlib.task_(byref(db._n))

