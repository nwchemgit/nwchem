EXE	= rdf
SRCS = \
main.f
OBJS =$(SRCS:.f=.o)

.f.o :
	${FC} -g -c $<
${EXE}: $(OBJS) ;
	${FC} -g -o ${EXE} ${OBJS}
	dsymutil ${EXE}

clean:
	rm -f *.o rdf *.dSYM
