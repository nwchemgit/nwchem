#!/bin/bash
sed -e 's/NXTVAL_/NXTVAL_OFF_/' $1/tcgmsg/tcgmsg-mpi/sndrcv.h > $1/tcgmsg/tcgmsg-mpi/sndrcv.h.mod
mv $1/tcgmsg/tcgmsg-mpi/sndrcv.h.mod  $1/tcgmsg/tcgmsg-mpi/sndrcv.h
sed -e 's/NXTVAL_/NXTVAL_OFF_/' $1/tcgmsg/tcgmsg-mpi/srftoc.h > $1/tcgmsg/tcgmsg-mpi/srftoc.h.mod
mv $1/tcgmsg/tcgmsg-mpi/srftoc.h.mod  $1/tcgmsg/tcgmsg-mpi/srftoc.h
sed -e 's/NXTVAL_(/NXTVAL_OFF_(/' $1/tcgmsg/tcgmsg-mpi/misc.c > $1/tcgmsg/tcgmsg-mpi/misc.c.mod
mv $1/tcgmsg/tcgmsg-mpi/misc.c.mod $1/tcgmsg/tcgmsg-mpi/misc.c
sed -e 's/NXTVAL_(/NXTVAL_OFF_(/' $1/tcgmsg/tcgmsg-mpi/nxtval.c > $1/tcgmsg/tcgmsg-mpi/nxtval.c.mod
mv $1/tcgmsg/tcgmsg-mpi/nxtval.c.mod $1/tcgmsg/tcgmsg-mpi/nxtval.c
sed -e 's/NXTVAL_(/NXTVAL_OFF_(/' $1/tcgmsg/tcgmsg-mpi/nxtval-armci.c > $1/tcgmsg/tcgmsg-mpi/nxtval-armci.c.mod
mv $1/tcgmsg/tcgmsg-mpi/nxtval-armci.c.mod $1/tcgmsg/tcgmsg-mpi/nxtval-armci.c
sed -e 's/ NXTVAL_(/ NXTVAL_OFF_(/' $1/tcgmsg/capi.c > $1/tcgmsg/capi.c.mod
mv $1/tcgmsg/capi.c.mod $1/tcgmsg/capi.c
sed -e 's/ NXTVAL_(/ NXTVAL_OFF_(/' $1/tcgmsg/fapi.c > $1/tcgmsg/fapi.c.mod
mv $1/tcgmsg/fapi.c.mod $1/tcgmsg/fapi.c
