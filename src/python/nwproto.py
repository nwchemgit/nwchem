# This file declares the top-level functions that can be called inside nwchem.
#
# More are available inside nwlib, but these represent all those accessible
# from the usual config file.
#
import os
from ctypes import *

# Shorthand for setting function prototypes
def decl_fn(a, *args):
	a.argtypes = args[:-1]
	a.restype = args[-1]

def void_fn(a, *args):
	decl_fn(a, *(args+(None,)))

def int_fn(a, *args):
	decl_fn(a, *(args+(c_int,)))

def dbl_fn(a, *args):
	decl_fn(a, *(args+(c_double,)))

strbuf = c_char*120
strbuf_p = POINTER(strbuf)
c_int_p = POINTER(c_int)
cwd = os.path.dirname(os.path.abspath(__file__))

# The library.
nwlib = cdll.LoadLibrary(os.path.join(cwd, "../../lib/MACX64/libnwchem.so"))

# Function type declaration
int_fn(nwlib.push_inp_cstring, c_char_p)
int_fn(nwlib.input_parse_, c_int_p) # just in case.
int_fn(nwlib.nwchem_init_, strbuf_p, strbuf_p, strbuf_p, strbuf_p)
void_fn(nwlib.nwchem_dtor_, c_int_p, c_int_p, c_int_p)
void_fn(nwlib.task_input_, c_int_p)
void_fn(nwlib.task_, c_int_p)
int_fn(nwlib.tce_input_, c_int) # handled specially by tce()
void_fn(nwlib.util_print_rtdb_load_, c_int_p, c_char_p, c_int)

