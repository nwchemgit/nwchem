      subroutine ccsd_mkt2(ncor,nocc,nvir,eorb,g_t2,
     $                     g_moint,t2,iprt)
      implicit none
      integer ncor,nocc,nvir,g_t2,g_moint,iprt
      double precision eorb(*),t2(*)
c
      integer nodes,iam
c
      integer jlo,jhi,alo,ahi
      integer a,b,i,j,av,icnt
C
#include "ccsd_len.fh"
#include "global.fh"
C
      Nodes = GA_NNodes()
      IAm = GA_NodeID()

      if (iprt.gt.5)print *,' making starting t2 '
c
      call ga_distribution(g_moint,iam,jlo,jhi,alo,ahi)
      av=0
      do a=ncor+nocc+1,ncor+nocc+nvir
       av=av+1
       if (av.ge.alo.and.av.le.ahi)then
        call ga_get(g_moint,1,lnoov,av,av,t2,lnoov)
        icnt=0
        do b=ncor+nocc+1,ncor+nocc+nvir
         do i=ncor+1,ncor+nocc
          do j=ncor+1,ncor+nocc
           icnt=icnt+1
           t2(icnt)=t2(icnt)/(eorb(i)+eorb(j)-eorb(a)-eorb(b))
          enddo
         enddo
        enddo
        if (iprt.gt.50)then
         print *,'t2 start ',av
         write(*,435)(t2(i),i=1,lnoov)
 435     format(4e15.5)
        endif
        call ga_put(g_t2,1,lnoov,av,av,t2,lnoov)
       endif
      enddo
c
      return
      end

