OBJ_OPTIMIZE=inp.o \
	     pdb_amber2nwchem.o \
	     pdb_readcommandlinearguments.o \
	     pre_atnum.o \
	     pre_namiup.o \
	     pre_protein.o \
	     util_print.o


all: $(OBJ_OPTIMIZE)
	$(FC) -o pdb_amber2nwchem $(OBJ_OPTIMIZE)

.F.o:
	$(FC) -c $<
