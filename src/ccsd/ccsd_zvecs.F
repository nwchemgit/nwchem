      subroutine ccsd_zvecs(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     $                      iprt,cmo,eorb,t1,hz1,hz2,
     $                      scra,scrb,z1,hiu,hia,habe,giu,gabe,
     $                      g_t2,g_z2,g_hz2,g_moint,
     $                      eccsd,max2e,mem2,g_ncoul,g_nexch,tklst)
C     $Id: ccsd_zvecs.F,v 2.3 1997-04-14 06:35:06 gg502 Exp $
      implicit none
      integer basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,iprt,
     &        g_t2,g_z2,g_hz2,g_moint,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
c
#include "ccsd_len.fh"
#include "ccsd_debug.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                 hz1(*),hz2(lnoov),scra(*),scrb(*),z1(*),
     &                 hiu(*),hia(*),habe(*),giu(*),gabe(*)
c
c
      integer i,j,a,b,av,bv,ad1,ad2,offt2,icnt,nvp1,ilo,ihi,
     &        jlo,jhi, inode
      double precision eccsd
c
      integer max2e,mem2
      Integer Nodes, IAm
      double precision one
      data one/1.0d00/
c
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      call ga_sync()
      call qenter('zvecs',0)
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
      call ccsd_zitf(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &               t1,z1,hiu,giu,habe,gabe,hia,hz1,idiis,
     &               cmo,eorb,iprt,g_t2,g_z2,g_hz2,g_ncoul,g_nexch,
     &               tklst)
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
        hz1(icnt)=hz1(icnt)/(eorb(i)-eorb(a))
       enddo
      enddo
      call daxpy(lnov,-one,hz1,1,z1,1)
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_z2,offt2+1,offt2+lnov,nvp1,nvp1,hz1,1)
       call ga_put(g_hz2,offt2+1,offt2+lnov,nvp1,nvp1,z1,1)
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
         call ga_get(g_hz2,ad1+1,ad1+lnoo,bv,bv,scra,lnoo)
         call ga_get(g_hz2,ad2+1,ad2+lnoo,av,av,scrb,lnoo)
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
         call ga_put(g_hz2,ad1+1,ad1+lnoo,bv,bv,scra,lnoo)
         call ga_put(g_hz2,ad2+1,ad2+lnoo,av,av,scrb,lnoo)
        endif
c -- end parallel task --
c
       enddo
      enddo
      call ga_sync()
c
      if (iprt.gt.50)then
       print *,' sym hz1 '
       write(6,189)(hz1(i),i=1,lnov)
       do bv=1,nvir
        call ga_get(g_hz2,offt2+1,offt2+lnoov,bv,bv,scra,1)
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
        call ga_get(g_z2,offt2+1,offt2+lnoov,bv,bv,scra,lnoov)
        call ga_get(g_hz2,offt2+1,offt2+lnoov,bv,bv,scrb,lnoov)
        icnt=0
        do a=ncor+nocc+1,ncor+nocc+nvir
         do i=ncor+1,ncor+nocc
          do j=ncor+1,ncor+nocc
           icnt=icnt+1
           scrb(icnt)=scrb(icnt)/(eorb(i)+eorb(j)-eorb(a)-eorb(b))
          enddo
         enddo
        enddo
c -- compute error vector in z2 location
        call daxpy(lnoov,-one,scrb,1,scra,1)
        call ga_put(g_z2,offt2+1,offt2+lnoov,bv,bv,scrb,lnoov)
        call ga_put(g_hz2,offt2+1,offt2+lnoov,bv,bv,scra,lnoov)
       endif
c -- end parallel task --
c
      enddo
      call ga_sync()
      call qexit('zvecs',0)
      return
      end
