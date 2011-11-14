EXE	= rdf
SRCS = \
main.f
OBJS =$(SRCS:.f=.o)

.f.o :
	gfortran -g -c $<
${EXE}: $(OBJS) ;
	gfortran -g -o ${EXE} ${OBJS}
	dsymutil ${EXE}

clean:
	rm -f *.o rdf *.dSYM
