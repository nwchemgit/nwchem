# This makefile has to be processed by GNU make
#
     TOPDIR = tools 
MAKESUBDIRS = for dir in $(SUBDIRS); do $(MAKE)  -C $$dir $@ || exit 1 ; done
    SUBDIRS = ma global tcgmsg-mpi LinAlg/lapack+blas server tcgmsg armci/src pario

include ../config/makefile.h

ifdef NWCHEM_TARGET_CPU
    TARGET_CPU = $(NWCHEM_TARGET_CPU)
endif

ifndef TARGET
info:
	@echo ''
	@echo '*** This is the Global Arrays package ****'
	@echo "*** you need to define TARGET machine, for example"
	@echo "***     setenv TARGET SUN"
	@echo "*** or "
	@echo "***     make TARGET=SUN"
	@echo ''
endif

export TARGET TARGET_CPU BINDIR


tools:   lib/$(TARGET) include
	@echo ''
	@echo '*** Building Parallel Tools ****'
	@echo ''
ifdef USE_ARMCI
		@echo '>>> making ARMCI remote memory copy library <<< '
		(cd ./armci/src; $(MAKE) LIB_DISTRIB=../../lib )
endif
	@echo '>>> making Memory Allocator (MA) library <<< '
	(cd ./ma; $(MAKE))
        ifdef USE_MPI
	      @echo '>>> making TCGMSG library on top of MPI <<< '
	      (cd ./tcgmsg-mpi; $(MAKE) MACHINE=$(TARGET) )
        else
	      @echo '>>> making TCGMSG library <<<'
	      (cd ./tcgmsg; $(MAKE) MACHINE=$(TARGET))
        endif
	@echo '>>> making Linear Algebra utility libraries <<< '
	(cd ./LinAlg/lapack+blas; $(MAKE))
	(cd ./LinAlg/summa; $(MAKE))
	@echo '>>> making Global Arrays library <<< '
	(cd ./global; $(MAKE) DIAG=PAR)
	(cd ./pario; $(MAKE))

lib/$(TARGET):
	$(MKDIR) -p $@

include:
	$(MKDIR) -p $@

include_stamp: tools

clean:
	$(MAKESUBDIRS) $^
	$(RM) -r include
ifdef USE_ARMCI
		(cd ./armci/src; $(MAKE) LIB_DISTRIB=../../lib clean )
endif

realclean: clean

.PHONY: cleanF
cleanF:
	@echo Not done here
source:
	
depend:
	@echo No dependency analysis here

