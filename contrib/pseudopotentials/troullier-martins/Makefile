#
# a makefile
#
CC=cc
F2C=f2c
F77=xlf
CDEBUGFLAGS=-g -O
CFLAGS= $(CDEBUGFLAGS)
#FFLAGS = -Zv -Wf"-dp" -Wf"-a stack"
FDEBUGFLAGS=
#FFLAGS = -Nn802 $(FDEBUGFLAGS)
FFLAGS = -O
#LIBS = -lF77 -lI77 -lm
LIBS =
LDEBUGFLAGS=-s
LDFLAGS = $(LDEBUGFLAGS) $(LIBS)

PROGS = kbconv atm

.SUFFIXES : .o .c .f

.f.o:
	${F77} ${FFLAGS} -c $*.f
.c.o:
	${CC} ${CFLAGS} -c $*.c

#run: all
#	/bin/csh -f run.cgscf >RESULTS

all : $(PROGS)

clean :
	rm -f *.o *.s kbconv atm core

OBJECTS1 = kbconv.o craylib.o mylib.o tridib.o
OBJECTS2 = atm.o craylib.o mylib.o

kbconv: $(OBJECTS1)
	${F77} -o kbconv ${CFLAGS} ${OBJECTS1} ${LDFLAGS}

atm: ${OBJECTS2}
	${F77} -o atm ${CFLAGS} ${OBJECTS2} ${LDFLAGS}
