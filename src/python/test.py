#!/usr/bin/env python

import sys
#from mpi4py import MPI
from nwchem import *

def main(argv):
	#mpi_init(argv)
	db = nwchem_init(1500)

	inp = """geometry autosym
  O    0.0    0.0    -0.02
  H   -0.74   0.0    -0.76
  H    0.74   0.0    -0.76
end
"""
	#nwlib.push_inp_cstring(inp)
	# This will not work next, since push_inp_cstring calls read_inp()
        #nwlib.input_parse_(byref(c_int(db)))

	geom(db, inp)
	bas(db,  """basis
  H library cc-pvdz
  O library cc-pvdz
end
""")
	driver(db, "driver; clear; end;")

	scf(db, "scf; print low; end;")
	nwtask(db, "scf optimize\n")
	
	x = db['geometry:geometry:coords'].reshape((3,3))
	y = db["task:energy"]
	if nwlib.ga_nodeid_() == 0:
		print x
		print y
	
def test_newgeom(db):
	x = db["geometry:geometry:coords"]
	x[3:] *= 1.01
	db["geometry:geometry:coords"] = x
	nwtask(db, "scf energy")
	return db["scf:energy"]

if __name__=="__main__":
	main(sys.argv)

