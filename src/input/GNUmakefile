#$Id: GNUmakefile,v 1.22 1995-11-19 19:41:34 d3g681 Exp $

 LIBRARY = libutil.a

     OBJ = input_parse.o input_mem_size.o memory_input.o input_set.o \
           input_start_opt.o input_title.o task_input.o scf_input.o \
           get_input_filename.o vectors_input.o stepper_input.o input_echo.o \
	   dft_input.o md_input.o gradients_input.o memory_defaults.o \
	   ccsd_input.o mcscf_input.o plnwv_input.o

     HEADERS = ../global/src/types.f2c.h

include ../config/makefile.h
include ../config/makelib.h

ifeq ($(TARGET),SGITFP)
$(LIBRARY_PATH)(get_input_filename.o):	get_input_filename.F
	f77 -g -d8 $(INCLUDES) $(DEFINES)  -c get_input_filename.F
endif

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(ccsd_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(dft_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas_staticP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(get_input_filename.o):	 $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(gradients_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(input_echo.o):	 $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(input_mem_size.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(input_parse.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(input_set.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(input_start_opt.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(input_title.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(mcscf_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(md_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(memory_defaults.o):	 $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(memory_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(plnwv_input.o):	 $(INCDIR)/types.f2c.h $(INCDIR)/rtdb.h $(INCDIR)/macommon.h 
$(LIBRARY_PATH)(scf_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(stepper_input.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(task_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(vectors_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
