      subroutine ccsd_pampt2(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     $                      iprt,cmo,eorb,t1,t2,tau,ht1,ht2,
     $                      scra,scrb,hiu,hia,habe,giu,gabe,
     $                      g_t2,g_ht2,g_moint,
     $                      eccsd,max2e,mem2,g_ncoul,g_nexch,tklst)
C     $Id: ccsd_pampt2.F,v 2.8 1997-04-14 06:34:59 gg502 Exp $
      implicit none
      integer basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,iprt,
     &        g_t2,g_ht2,g_moint,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
c
#include "ccsd_len.fh"
      double precision cmo(nbf,nbf),eorb(nbf),t1(*),t2(lnoov),
     &                 tau(lnoov),ht1(*),ht2(lnoov),scra(*),scrb(*),
     &                 hiu(*),hia(*),habe(*),giu(*),gabe(*)
c
      integer nodes,iam
#include "ccsd_debug.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer i,j,a,b,av,bv,ad1,ad2,offt2,icnt,nvp1,ilo,ihi,
     &        jlo,jhi,inode
      double precision eccsd,eccsd1,eccsd2
c
      integer max2e,mem2
      double precision zip,one,two
c
      data zip/0.0d00/,one/1.0d00/,two/2.0d00/
c
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
      call ga_sync()
      call qenter('pampt',0)
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
c
      offt2=(idiis-1)*lnoov
      nvp1=nvir+1
      lnobj=lnoov+lnooo
c
c ----------------------------------------------------------
c - alternative vvvv term
c ----------------------------------------------------------
c
      call ccsd_itf(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &              hiu,giu,habe,gabe,hia,ht1,idiis,
     &               cmo,eorb,iprt,g_t2,g_ht2,g_ncoul,g_nexch,tklst)
c
c ----------------------------------------------------------
c symmetrize the sigma vector and evaluated error vectors
c ----------------------------------------------------------
c
      if (iprt.gt.5)print *,'start symmetrization'
      icnt=0
      do a=ncor+nocc+1,ncor+nocc+nvir
       do i=ncor+1,ncor+nocc
        icnt=icnt+1
        ht1(icnt)=ht1(icnt)/(eorb(i)-eorb(a))
       enddo
      enddo
      call daxpy(lnov,-one,ht1,1,t1,1)
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,ht1,1)
       call ga_put(g_ht2,offt2+1,offt2+lnov,nvp1,nvp1,t1,1)
      endif
c
      inode=-1
      do bv=1,nvir
       do av=1,bv
c
c -- parallel task --
        inode=inode+1
        if (mod(inode,nodes).eq.iam)then
         ad1=offt2+(av-1)*lnoo
         ad2=offt2+(bv-1)*lnoo
         call ga_get(g_ht2,ad1+1,ad1+lnoo,bv,bv,scra,lnoo)
         call ga_get(g_ht2,ad2+1,ad2+lnoo,av,av,scrb,lnoo)
c         print *,' scra ',(scra(i),i=1,lnoo)
c         print *,' scrb ',(scrb(i),i=1,lnoo)
         do i=1,nocc
          do j=1,nocc
           scra((i-1)*nocc+j)=scra((i-1)*nocc+j)+scrb((j-1)*nocc+i)
          enddo
         enddo
         do i=1,nocc
          do j=1,nocc
           scrb((i-1)*nocc+j)=scra((j-1)*nocc+i)
          enddo
         enddo
         call ga_put(g_ht2,ad1+1,ad1+lnoo,bv,bv,scra,lnoo)
         call ga_put(g_ht2,ad2+1,ad2+lnoo,av,av,scrb,lnoo)
        endif
c -- end parallel task --
c
       enddo
      enddo
      call ga_sync()
c
      if (iprt.gt.50)then
       print *,' sym ht1 '
       write(6,189)(ht1(i),i=1,lnov)
       do bv=1,nvir
        call ga_get(g_ht2,offt2+1,offt2+lnoov,bv,bv,scra,1)
        print *,'sym sig2 for bv ',bv
        write(6,189)(scra(i),i=1,lnoov)
  189   format(4e14.4)
       enddo
      endif
c
      do b=ncor+nocc+1,ncor+nocc+nvir
       bv=b-nocc-ncor
c
c -- parallel task --
       if (bv.ge.ilo.and.bv.le.ihi)then
        call ga_get(g_t2,offt2+1,offt2+lnoov,bv,bv,scra,lnoov)
        call ga_get(g_ht2,offt2+1,offt2+lnoov,bv,bv,scrb,lnoov)
        icnt=0
        do a=ncor+nocc+1,ncor+nocc+nvir
         do i=ncor+1,ncor+nocc
          do j=ncor+1,ncor+nocc
           icnt=icnt+1
           scrb(icnt)=scrb(icnt)/(eorb(i)+eorb(j)-eorb(a)-eorb(b))
          enddo
         enddo
        enddo
c -- compute error vector in t2 location
        call daxpy(lnoov,-one,scrb,1,scra,1)
        call ga_put(g_t2,offt2+1,offt2+lnoov,bv,bv,scrb,lnoov)
        call ga_put(g_ht2,offt2+1,offt2+lnoov,bv,bv,scra,lnoov)
       endif
c -- end parallel task --
c
      enddo
c
c ----------------------------------------------------------
c compute energy
c ----------------------------------------------------------
c
      if (iprt.gt.5)print *,'start energy'
      eccsd1=zip
      eccsd2=zip
      do bv=1,nvir
c
c -- parallel task --
       if (bv.ge.ilo.and.bv.le.ihi)then
        call ga_get(g_t2,offt2+1,offt2+lnoov,bv,bv,scra,lnoov)
        call ccsd_mktau(nocc,nvir,bv,ht1,scra,scrb)
        call ga_get(g_moint,1,lnoov,bv,bv,scra,lnoov)
        do av=1,nvir
         do i=1,nocc
          do j=1,nocc
           ad1=(av-1)*lnoo+(i-1)*nocc+j
           ad2=(av-1)*lnoo+(j-1)*nocc+i
           eccsd2=eccsd2+(two*scrb(ad1)-scrb(ad2))*scra(ad1)
c$$$       write(6,129)nocc+bv,nocc+av,i,j,scrb(ad1),scrb(ad2),
c$$$     & scra(ad1),eccsd2
c$$$  129  format('<<<<<',4i3,4e15.5)
          enddo
         enddo
        enddo
       endif
c -- end parallel task --
c
      enddo
      call ga_dgop(msg_cc_ecc,eccsd2,1, '+')
      eccsd=eccsd1+eccsd1+eccsd2
      if (iprt.gt.5.and.iam.eq.0)then
       write(6,1238)eccsd1,eccsd2,eccsd
 1238  format(/,' eccsd1 ',e20.12,/
     &          ' eccsd2 ',e20.12,/
     &          ' eccsd  ',e20.12,/)
      endif
c
      call ga_sync()
      call qexit('pampt',0)
      return
      end
