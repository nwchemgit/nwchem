#
# $Id: GNUmakefile,v 1.9 1999-07-27 21:04:56 d3e129 Exp $
#
# This makefile has to be processed by GNU make
#
     TOPDIR = ${NWCHEM_TOP}/src/tools 
MAKESUBDIRS = for dir in $(SUBDIRS); do $(MAKE)  -C $$dir $@ || exit 1 ; done
    SUBDIRS = ma global tcgmsg-mpi LinAlg/lapack+blas tcgmsg armci/src pario
  USE_ARMCI = yes
ifndef USE_MPI
  EXTRA_INC = "TCG_INCLUDE=../../include"
      EXTRA = "MSG_COMMS=TCGMSG"
     MP_VER = TCG
else
     MP_VER = MPI
endif

include ../config/makefile.h

ifdef NWCHEM_TARGET_CPU
    TARGET_CPU = $(NWCHEM_TARGET_CPU)
endif

ifdef OLD_GA
    unexport USE_ARMCI
    GA_VER = O
else
    export USE_ARMCI
    GA_VER = N
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

STAMP = $(TARGET)-$(MP_VER)-$(GA_VER).stamp


tools:  $(STAMP) lib/$(TARGET) include
	@echo ''
	@echo '*** Building Parallel Tools ****'
	@echo ''
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
ifndef  OLD_GA
ifdef   USE_ARMCI
		@echo '>>> making ARMCI remote memory copy library <<< '
		(cd ./armci/src; \
		$(MAKE) LIB_DISTRIB=../../lib $(EXTRA) $(EXTRA_INC) )
endif
endif
	@echo '>>> making Global Arrays library <<< '
	(cd ./global; $(MAKE) DIAG=PAR )
	(cd ./pario; $(MAKE))

lib/$(TARGET):
	$(MKDIR) -p $@

include: 
	$(MKDIR) -p $@

include_stamp: tools

clean:
	$(MAKESUBDIRS) $^
	$(RM) -r include/* *.stamp
ifdef USE_ARMCI
		(cd ./armci/src; $(MAKE) LIB_DISTRIB=../../lib clean )
endif

realclean: clean


$(STAMP): 
	$(MAKE) clean
	echo "" > $(STAMP)

.PHONY: cleanF
cleanF:
	@echo Not done here
source:
	
depend:
	@echo No dependency analysis here
sngl_to_dbl dbl_to_sngl:
	@echo No blas conversion done here
#
# This is a convenience target that will make the TAGS file for current 
# checked out source tree.  This is only useful if you know something 
# about emacs.  Note: find, grep and etags must be in your path.
#
tags_file:
	find . \( -name "*.[cfFh]" -o -name "*.fh" \) -print | grep -v "\./include" | grep -v "\./tools/include" | etags -
tags_clean:
	find . -name TAGS -print -exec rm -f "{}" ";"
