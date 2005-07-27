#
# $Id: GNUmakefile,v 1.23 2005-07-27 23:05:13 marat Exp $
#

LIBRARY = libdriver.a

OBJ_OPTIMIZE = driver_input.o 

OBJ = oniom.o opt_drv.o lbfgs.o lbfgs-b.o lbfgs_driver.o lbfgs-b_driver.o

USES_BLAS = opt_drv.F oniom.F lbfgs-b.F

include ../config/makefile.h
include ../config/makelib.h

#FFLAGS  = -qsource -qxref=full -bloadmap:load.map 
