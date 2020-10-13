OBJ_OPTIMIZE=inp.o \
	     pdb_large.o \
	     pdb_readcommandlinearguments.o 


all: $(OBJ_OPTIMIZE)
	$(FC) -o pdb_large $(OBJ_OPTIMIZE)

.F.o:
	$(FC) -c $<
