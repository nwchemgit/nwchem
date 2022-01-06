from nwproto import *
import numpy as np
#import numpy.ctypeslib as nc
from ctypes import *

#int rtdb_parallel(const int mode)
int_fn(nwlib.rtdb_parallel, c_int)
int_fn(nwlib.rtdb_parallel, c_int)
int_fn(nwlib.rtdb_open, c_char_p,  c_char_p, c_int_p)
int_fn(nwlib.rtdb_clone, c_int, c_char_p)
int_fn(nwlib.rtdb_getfname, c_int, c_char_p)
int_fn(nwlib.rtdb_close, c_int, c_char_p)
int_fn(nwlib.rtdb_put, c_int, c_char_p, c_int, c_int, c_void_p)
int_fn(nwlib.rtdb_get, c_int, c_char_p, c_int, c_int, c_void_p)
int_fn(nwlib.rtdb_get_info, c_int, c_char_p, c_int_p, c_int_p, c_char*26)
int_fn(nwlib.rtdb_first, c_int, c_int, c_char_p)
int_fn(nwlib.rtdb_next, c_int, c_int, c_char_p)
int_fn(nwlib.rtdb_print, c_int, c_int)
int_fn(nwlib.rtdb_delete, c_int, c_char_p)

# from tools/ma/macommon.h
#ma_tp = {1000:c_char,
#	1001:c_int,
#	1002:c_longlong,
#	1003:c_float,
#	1004:c_double,
#	1005:c_longdouble,
#	1006:c_float*2, # complex
#	1007:c_double*2, # complex
#	1008:c_longdouble*2 } # complex
ma_tp = {1000:np.char,
	 1001:np.int32,
	 1002:np.int64,
	 1003:np.float32,
	 1004:np.float64,
	 1005:np.float128,
	 1006:np.complex64,
	 1007:np.complex128,
	 1008:np.complex256,
	 1009:np.char,
	 1010:np.int,
	 1011:np.bool,
	 1012:np.float32,
	 1013:np.float64,
	 1014:np.complex64,
	 1015:np.complex128,
	 1016:np.int64 }
#define MT_F_BYTE     (MT_BASE + 9)
#define MT_F_INT      (MT_BASE + 10)
#define MT_F_LOG      (MT_BASE + 11)
#define MT_F_REAL     (MT_BASE + 12)
#define MT_F_DBL      (MT_BASE + 13)
#define MT_F_SCPL     (MT_BASE + 14)
#define MT_F_DCPL     (MT_BASE + 15)
#define MT_C_LONGLONG (MT_BASE + 16)

def lookup_ma_dtype(t):
	#if t == np.char:
	#	return 1000
	if t == np.int32:
		return 1001
	elif t == np.int64:
		return 1002
	elif t == np.float32:
		return 1003
	elif t == np.float64:
		return 1004
	elif t == np.float128:
		return 1005
	elif t == np.complex64:
		return 1006
	elif t == np.complex128:
		return 1007
	elif t == np.complex256:
		return 1008
	elif t == np.bool:
		return 1011
	raise KeyError

class rtdb:
#    mode     = 'new'     Open only if it does not exist already
#               'old',    Open only if it does exist already
#               'unknown' Create new or open existing (preserving contents)
#               'empty'   Create new or open existing (deleting contents)
#               'scratch' Create new or open existing (deleting contents)
#                         and automatically delete upon closing.  Also, items
#                         cached in memory are not written to disk.
	def __init__(self, name, mode="unknown"):
		self._n = None
		self._fname = ""
		self._opened = False
		if type(name) == type(0):
		    self._n = c_int(name)
		elif type(name) == type("string"):
		    p = c_int()
		    if nwlib.rtdb_open(name, mode, byref(p)) == 0:
			raise RuntimeError, "Error opening RTDB."
		    self._n = p
		    self._opened = True
		else:
		    raise RuntimeError, "usage: rtdb(name or id, mode)"
		self._fname = name
	def __del__(self):
		if self._opened:
		    nwlib.rtdb_close(self._n, "keep")
	# return the filename of the rtdb
	def getfname(self):
		# unsafe method
		#fname = (c_char*1024)()
		#nwlib.rtdb_getfname(self._n, fname)
		#return fname.value
		return self._fname

	# Copy the database to a new file, to be named
	# <self.fname()>.suffix
	def clone(self, suf):
		if nwlib.rtdb_clone(self._n, suf) == 0:
			raise RuntimeError, "Error cloning RTDB."
	# store the string or array x in the rtdb at the given key
	def put(self, key, x):
		# If the key exists, coerce the type to the key's type.
		try:
			tp, _, _ = self.stat(key)
		except KeyError:
			tp = None
		
		if type(x) == type("string"):
			if tp != None or tp != 1000:
				raise ValueError, "Can't overwrite a non-string with a string."
			tp = 1000
			n = len(x)
			x = cast(create_string_buffer(x), c_void_p)
		else:
			if tp != None:
				x = x.astype(ma_tp[tp])
			else:
				tp = lookup_ma_dtype(x.dtype)
			n = np.prod(x.shape)
			x = x.ctypes.data_as(c_void_p)
		if nwlib.rtdb_put(self._n, key, tp, n, x) == 0:
			raise KeyError
	# returns ma_type, number, and the array
	def get(self, key):
		tp, n, date = self.stat(key)
		if tp == 1000:
		    x = (c_char*(n+1))()
		    if nwlib.rtdb_get(self._n, key, tp, n, cast(x, c_void_p)) == 0:
			raise KeyError
		    x = x.value
		else:
		    #x = (ma_tp[tp]*n)()
		    x = np.zeros(n, dtype=ma_tp[tp])
		    if nwlib.rtdb_get(self._n, key, tp, n, x.ctypes.data_as(c_void_p)) == 0:
			raise KeyError
		return x
	# stat of the object - returns ma_type, number, and date
	def stat(self, key):
		tp = c_int()
		n = c_int()
		date = (c_char*26)()
		if nwlib.rtdb_get_info(self._n, key, byref(tp), byref(n), date) == 0:
			raise KeyError
		return tp.value, n.value, date.value
	# val = True => print values and keys.
	def ls(self, val=False):
		nwlib.rtdb_print(self._n, int(val))
	# make the rtdb act like a dictionary.
	def __getitem__(self, key):
		return self.get(key)
	def __setitem__(self, key, val):
		return self.put(key, val)
	def __delitem__(self, key):
		nwlib.rtdb_delete(self._n, key)
	#def __getattribute__(self, name):
        #  if name in ("A", "B", "C"):
        #    return "Immutable value of %s" % name
        #  else:
        #    # This should trigger the default behavior for any other
        #    # attribute name.
        #    raise AttributeError()
	#def __setattr__(self, name, value):
        #  if name in ("A", "B", "C"):
        #    raise AttributeError("%s is an immutable attribute.")
        #  else:
        #    dict.__getitem__(self, 
	#def __delattr__(self, name):
	# return an iterator over my keys.
	def iterkeys(self):
		name = (c_char*1024)()
		x = nwlib.rtdb_first(self._n, 1024, name)
		while x:
			yield name.value
			x = nwlib.rtdb_next(self._n, 1024, name)
		raise StopIteration
	def __iter__(self):
		return self.iterkeys()
	def keys(self):
		return [k for k in self.iterkeys()]
