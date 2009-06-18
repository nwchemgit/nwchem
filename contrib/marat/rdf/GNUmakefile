EXE	= rdf
SRCS = \
main.f
OBJS =$(SRCS:.f=.o)
${EXE}: $(OBJS) ;
	${FC} -o ${EXE} ${OBJS}

clean:
	rm -f *.o rdf 
