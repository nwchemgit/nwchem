#$Id: GNUmakefile,v 1.17 1995-10-24 19:56:27 rg240 Exp $
 LIBRARY = libutil.a
#LIBRARY = libinput.a

     OBJ = input_parse.o input_mem_size.o memory_input.o input_set.o \
           input_start_opt.o input_title.o task_input.o scf_input.o \
           get_input_filename.o vectors_input.o stepper_input.o input_echo.o \
	   dft_input.o md_input.o gradients_input.o memory_defaults.o \
	   ccsd_input.o mcscf_input.o

include ../config/makefile.h
include ../config/makelib.h

ifeq ($(TARGET),SGITFP)
get_input_filename.o:	get_input_filename.F
	f77 -g -d8 $(INCLUDES) $(DEFINES)  -c get_input_filename.F -o get_input_filename.o	
endif
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
ccsd_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
dft_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas_staticP.fh $(INCDIR)/bas.fh 
get_input_filename.o :	 $(INCDIR)/inp.fh $(INCDIR)/global.fh 
gradients_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
input_echo.o :	 $(INCDIR)/inp.fh 
input_mem_size.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
input_parse.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
input_set.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
input_start_opt.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
input_title.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh 
md_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
memory_defaults.o :	 $(INCDIR)/inp.fh 
memory_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
scf_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
stepper_input.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh 
task_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh 
vectors_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
