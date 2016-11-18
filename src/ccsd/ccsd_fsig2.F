      subroutine ccsd_fsig2(nocc,nvir,bv,u,t1,ht1,hiu)
C     $Id$
      implicit none
      integer nocc,nvir,bv,u
      double precision t1(nocc,nvir),
     &                 ht1(nocc,nvir),hiu(nocc,nvir)
c
      External ddot
      double precision ddot
C
#include "ccsd_debug.fh"
#include "ccsd_len.fh"
c
c term s2 n**3
      if (dos(2).ne.0)then
c       print *,' s2 '
       ht1(u,bv)=ht1(u,bv)-ddot(nocc,hiu(1,u),1,t1(1,bv),1)
      endif
c
      return
      end
