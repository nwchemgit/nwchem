#!/bin/sh
# eine kleine TURBOMOLE-Leihgabe

host=`hostname`
set `date`
cat > date.inc <<EOF
      character, parameter :: luvers*80 = 
     ,  "Version compiled $3 $2 $6 at $4 on "//
     ,  "$host"
EOF
