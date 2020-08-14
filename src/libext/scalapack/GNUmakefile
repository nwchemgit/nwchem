#$Id$
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../config/makefile.h

ifndef USE_MPI
nompierror:
	@echo "  "
	@echo "  " You must set USE_MPI and
	@echo "  " have a working MPI installation
	@echo "  " with mpif90 in your PATH
	@echo "  "  to compile NWChem
	@echo
	@exit 1
endif # USE_MPI

../lib/libnwc_scalapack.a:
	./build_scalapa.sh 
LIB_TARGETS += scalapack

include ../../config/makelib.h

clean:
	@rm -rf  scalapack* ../lib/libnwc_sca*



