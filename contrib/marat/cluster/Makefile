#
# Makefile For Response/Local Potential Calculation
#
EXE     = cluster
SRCS = \
main.f smd_pdb.f
OBJS =$(SRCS:.f=.o)
${EXE}: $(OBJS) ;
	${FC} -o ${EXE} ${OBJS}

clean:
	rm -rf *.o cluster
