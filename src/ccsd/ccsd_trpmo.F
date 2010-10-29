      subroutine ccsd_trpmo(ncor,nocc,nvir,nact,nbf,scr,scra,t1,g_nt2,
     $                       g_objv,g_ncoul,g_nexch,
     $                       iprt)
      implicit none
C     $Id$
      integer ncor,nocc,nvir,nact,nbf,g_objv,g_nt2,
     $        g_ncoul,g_nexch,iprt
      double precision scr(*),scra(*),t1(*)
c
#include "ccsd_len.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer iam
      integer ihi,ilo,jhi,jlo,i,j,k,b,av,ad1,ad2
      integer nvp1,iin,ad1g,jin
c
      if (iprt.gt.5)print *,' sort_ccsd '
c
      iam = ga_nodeid()
      call ga_sync
#ifdef CRAYXT
      call util_mpinap(100)
#endif
c
      nvp1=nvir+1
      lnoo=nocc*nocc
      lnov=nocc*nvir
c
      call ga_distribution(g_objv,iam,jlo,jhi,ilo,ihi)
      do av=1,nvir
       if (av.ge.ilo.and.av.le.ihi)then
cedo         do i=1,nocc
cstaggering to avoid congestion in ga_get
         do iin=ga_nodeid()+1,ga_nodeid()+nocc
            i=mod(iin,nocc)+1
           ad1g=(i-1)*nact+nocc+av
cold           do j=1,nocc
           do jin=ga_nodeid()+1,ga_nodeid()+nocc
            j=mod(jin,nocc)+1
             ad1=(j-1)*nact
             call ga_get(g_nexch,ad1+1,ad1+nact,
     A            ad1g,ad1g,scra,nact)
           do b=1,nvir
c             ad1=(j-1)*nact+nocc+b
             ad2=(b-1)*lnoo+(i-1)*nocc+j
             scr(ad2)=scra(nocc+b)
           enddo
           do k=1,nocc 
c             ad1=(j-1)*nact+k
             ad2=(k-1)*lnoo+(i-1)*nocc+j
             scr(lnoov+ad2)=scra(k)
           enddo
        enddo
        enddo
        if (iprt.gt.50)then
         print *,'mo object ',av
         write(*,435)(scr(i),i=1,lnoov+lnooo)
 435     format(4e15.5)
        endif
        call ga_put(g_objv,1,lnoov+lnooo,av,av,scr,lnoov+lnooo)
c
       endif
      enddo
      call ga_sync
c
c - put vectors in g_objv for use in triples, overwrite [vvoo] ints
      lnobj=lnoov+lnooo
      if (iam.eq.0)then
       call ga_put(g_objv,lnobj+1,lnobj+lnov,nvp1,nvp1,t1,lnov)
      endif
c
      call ga_distribution(g_objv,iam,jlo,jhi,ilo,ihi)
      do av=1,nvir
       if (av.ge.ilo.and.av.le.ihi)then
cedo         do i=1,nocc
cstaggering to avoid congestion in ga_get
         do iin=ga_nodeid()+1,ga_nodeid()+nocc
            i=mod(iin,nocc)+1
             ad1g=(i-1)*nvir+av
cold             call ga_get(g_nt2,1,lnov,ad1g,ad1g,scra,lnov)
             do jin=ga_nodeid()+1,ga_nodeid()+nocc
                j=mod(jin,nocc)+1
                ad1=(j-1)*nvir
             call ga_get(g_nt2,ad1+1,ad1+nvir,ad1g,ad1g,scra,lnov)
                do b=1,nvir
                   ad2=(b-1)*lnoo+(i-1)*nocc+j
                   scr(ad2)=scra(b)
                enddo
             enddo
          enddo
          call ga_put(g_objv,lnobj+1,lnobj+lnoov,av,av,scr,lnoov)
       endif
      enddo
      call ga_sync
c
      return
      end

