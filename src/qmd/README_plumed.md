#Notes on NWChem-Plumed interface

Used Plumed version 2.6

## Plumed Installation

https://www.plumed.org/doc-v2.6/user-doc/html/_installation.html

* 32bit integer BLAS/LAPACK (e.g. MKL,OpenBLAS) are required both in NWChem and Plumed.
In order to get Plumed to link the correct Blas/Lapack libraries, the variales LDFLAGS and LIBS must be supplied to
autconf.
* No need of MPI for now. It is better to disable the MPI option with --disable-mpi
Example on EMSL Tahoma 
```
FC=ifort  LDFLAGS=-L$MKLROOT/lib/intel64 LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -lpthread -lm" ./configure --prefix=/home/edo/tahoma/apps/plumed262.intel20u2 --disable-mpi
make
make install
```

## Compilation of the NWChem Interface

* modify PATH to point to the location of the plumed command (`/home/edo/tahoma/apps/plumed262.intel20u2/bin` in the 
example above)
* set the env. variable `USE_PLUMED=y`
* include the qmd module in `NWCHEM_MODULES`

## Files of the NWChem-Plumed interface

* input and output files reside in `permanent_dir`
* the Plumed input file can be either suffix.`plumed.dat` our just `plumed.dat`
* the Plumed output file will be name suffix.`plumed.out` (and its content will be copied to stdout, too)