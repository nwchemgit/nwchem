include ./Makefile.config

ARCHDIR=$(ARCH)
#MOBJECTS=$(addprefix $(ARCHDIR)/,workspace.o)
#OBJECTS=$(addprefix $(ARCHDIR)/,$(filter-out workspace.o printversion.o,$(patsubst %.f,%.o,$(wildcard *.f))))
OBJECTS=$(addprefix $(ARCHDIR)/,$(filter-out printversion.o,$(patsubst %.f,%.o,$(wildcard *.f))))
COBJECTS=$(addprefix $(ARCHDIR)/,$(patsubst %.c,%.o,$(wildcard *.c)))
INCLUDES=$(filter-out date.inc,$(wildcard *.inc))

.PHONY: all
all: dir prog

.PHONY: dir 
dir:
	@set -e ; if [ ! -d $(ARCHDIR) ] ; then mkdir $(ARCHDIR) ; fi

.PHONY: prog
prog:  $(ARCHDIR)/lucia.x

.PHONY: clean
clean: 
	@rm  $(ARCHDIR)/*.o ; rm $(ARCHDIR)/lucia.x 

$(ARCHDIR)/lucia.x  : $(MOBJECTS) $(OBJECTS) $(COBJECTS)
		make printversion.o
		$(F90) $(OPTFLAGS) $(LLIBS) -o $(ARCHDIR)/lucia.x $(ARCHDIR)/*.o\
                       $(MATHLIBS)
#$(ARCHDIR)/lucia.x  : printversion.o $(OBJECTS) $(COBJECTS)
#		$(F90) $(OPTFLAGS) $(LLIBS) -o $(ARCHDIR)/lucia.x $(ARCHDIR)/*.o\
#                       -llapack -lcblas -latlas -lsvml -L/home/andreas/prog/TURBOMOLE/blaslib/i786-pc-linux-gnu/ -lblaslib -Vaxlib
#$(ARCHDIR)/lucia.x  : printversion.o $(OBJECTS) $(COBJECTS)
#		$(F90) $(OPTFLAGS) $(LLIBS) -o $(ARCHDIR)/lucia.x $(ARCHDIR)/*.o\
#                       -llapack -lcblas -latlas -lsvml -Vaxlib
# Object rules

date.inc: 
	./date.sh

printversion.o: date.inc
	$(F90) $(OPTFLAGS) -c printversion.f -o $(ARCHDIR)/printversion.o
	rm date.inc

$(MOBJECTS): $(ARCHDIR)/%.o: %.f 
	$(F90) $(INCLUDE_DIR) $(OPTFLAGS) $(MODFLAG) -c $< -o $@
$(OBJECTS): $(ARCHDIR)/%.o: %.f $(MOBJECTS) $(INCLUDES) 
	$(F90) $(INCLUDE_DIR) $(OPTFLAGS) -c $< -o $@
$(COBJECTS): $(ARCHDIR)/%.o: %.c
	$(CC) $(INCLUDE_DIR) $(CPP_FLAGS) -c $< -o $@
#$(OBJECTS):
#		src=$(subst $(ARCHDIR)/,,$*.f)
#		$(F90) $(INCLUDE_DIR) $(OPTFLAGS) -c $(src) -o $*.o
#$(COBJECTS):
#		src=$(subst $(ARCHDIR)/,,$*.c)
#		$(CC) $(INCLUDE_DIR) $(CPP_FLAGS) $(OPTFLAGS) -c $(src) -o $*.o
#
#
# Suffix rules
#
#.SUFFIXES : .f .F .c .o
#
#.f.o:
#		src=$(subst $(ARCHDIR)/,,$*.f)
#		$(F90) $(INCLUDE_DIR) $(OPTFLAGS) -c $(src) -o $*.o
#.F.o:
#		src=$(subst $(ARCHDIR)/,,$*.F)
#		$(F90) $(INCLUDE_DIR) $(CPP_FLAGS) $(OPTFLAGS) -c $(src) -o $*.o 
#.c.o:
#		src=$(subst $(ARCHDIR)/,,$*.c)
#		$(CC) $(INCLUDE_DIR) $(CPP_FLAGS) $(OPTFLAGS) -c $(src) -o $*.o
# 
