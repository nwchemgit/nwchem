#
# $Id$
#

LIBRARY = libdriver.a

OBJ_OPTIMIZE = driver_input.o 

OBJ =oniom.o  opt_drv.o lbfgs.o  lbfgs-b.o  lbfgs_driver.o lbfgs_utils.o lbfgs-b_driver.o

USES_BLAS = opt_drv.F\
	oniom.F\
	lbfgs-b.F\
	lbfgs-b_driver.F \
	lbfgs.F

include ../config/makefile.h
include ../config/makelib.h

#FFLAGS  = -qsource -qxref=full -bloadmap:load.map 
