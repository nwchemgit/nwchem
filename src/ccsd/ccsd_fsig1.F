      subroutine ccsd_fsig1(av,bv,nocc,nvir,nbf,iprt,
     &                      t1,ht1,habe)
      implicit none
c
      integer av,bv,nocc,nvir,nbf,iprt
c
       double precision t1(nocc,nvir),ht1(nocc,nvir),habe(nvir,nvir)
c
      integer nodes,iam
#include "ccsd_len.fh"
#include "ccsd_debug.fh"
#include "global.fh"
c
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
C
c term s1 n**3
      if (dos(1).ne.0)then
c       print *,' s1 '
       call daxpy(nocc,habe(bv,av),t1(1,av),1,ht1(1,bv),1)
      endif
c
      return
      end
