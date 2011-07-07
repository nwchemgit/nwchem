EXE	= amber-helper
SRCS = \
main.F test.F
OBJS =$(SRCS:.f=.o)
${EXE}: $(OBJS) ;
	${FC} -o ${EXE} ${OBJS}

clean:
	rm -f *.o amber-helper
