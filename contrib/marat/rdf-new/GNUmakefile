EXE	= rdf
SRCS = \
main.f
OBJS =$(SRCS:.f=.o)
${EXE}: $(OBJS) ;
	${FC} -g -o ${EXE} ${OBJS}

clean:
	rm -f *.o rdf 
