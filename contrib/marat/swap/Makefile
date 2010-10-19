#
# Makefile For Response/Local Potential Calculation
#
EXE     = swap_pdb
SRCS = \
main.f smd_pdb.f
OBJS =$(SRCS:.f=.o)
${EXE}: $(OBJS) ;
	${FC} -o ${EXE} ${OBJS}

clean:
	rm -rf *.o swap_pdb
