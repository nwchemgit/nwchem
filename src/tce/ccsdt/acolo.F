c
c
c auxiliary active-type logical functions
c acolo and acolo_o 
c
      LOGICAL FUNCTION ACOLO(p4b,p5b,p6b,h1b,h2b,h3b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER p4b
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER iocc
      INTEGER ivirt
      acolo = .false.
      iocc=int_mb(k_active+h1b-1)+int_mb(k_active+h2b-1)+ 
     &     int_mb(k_active+h3b-1)
      ivirt=int_mb(k_active+p4b-1)+int_mb(k_active+p5b-1)+
     &     int_mb(k_active+p6b-1) 
      if(iocc.ge.numact.and.ivirt.ge.numact) acolo=.true. 
      RETURN
      END
c $Id$
