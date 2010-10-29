      subroutine ccsd_zvecs(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     $                      iprt,cmo,eorb,t1,hz1,
     $                      scra,scrb,z1,hiu,hia,habe,giu,gabe,
     $                      g_nt2,g_z2,g_hz2,g_nz2,g_nhz2,
     $                      g_ncoul,g_nexch,tklst, Tol2E)
C     $Id$
      implicit none
      integer basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,iprt,
     &        g_nt2,g_z2,g_hz2,g_nz2,g_nhz2,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
c
      double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                 hz1(*),scra(*),scrb(*),z1(*),
     &                 hiu(*),hia(*),habe(*),giu(*),gabe(*)
      Double Precision Tol2e
c
#include "ccsd_len.fh"
#include "ccsd_debug.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "ccsdps.fh"
#include "util.fh"
c
      integer i,j,a,b,bv,ad1,ad2,ad3,ad4,offt2,icnt,nvp1,
     &        ilo,ihi,jlo,jhi,lhoo,lhoov
c
      Integer Nodes, IAm
      Logical oprint
      double precision one
      data one/1.0d00/
c
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      call ga_sync()
      if (occsdps) then
         call pstat_on(ps_zvecs)
      else
         call qenter('zvecs',0)
      endif
c
      oprint = util_print('information', print_low)
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
      lnov=nocc*nvir
      lnvv=nvir*nvir
      offt2=(idiis-1)*lhoov
c
      call ga_get(g_z2,offt2+1,offt2+lnov,nvir+1,nvir+1,
     &            z1,lnov)
c
c ------------------------------------------------------------
c fill square GAs
c ------------------------------------------------------------
c
      call ga_distribution(g_z2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          call ga_get(g_z2,offt2+1,offt2+lhoov,a,a,
     &                 scra,lhoov)
          do i=1,nocc
            call dfill(lnov,0.0d00,scrb,1)
            do j=1,i
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                ad2=(j-1)*nvir+b
                scrb(ad2)=scra(ad1)
              enddo
            enddo
            ad1=(i-1)*nvir+a
            call ga_put(g_nz2,1,i*nvir,ad1,ad1,scrb,i*nvir)
          enddo
        endif
      enddo
c
      call ga_sync()
c
      call ga_distribution(g_nz2,iam,jlo,jhi,ilo,ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.ilo.and.ad1+1.le.ihi)then
          do j=1,i-1
            ad2=(j-1)*nvir
            if (ad2+1.ge.jlo.and.ad2+1.le.jhi)then
              call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              do a=1,nvir
              do b=1,nvir
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                scrb(ad4)=scra(ad3)
              enddo
              enddo
              call ga_put(g_nz2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_zero(g_nhz2)
      call ccsd_pzamp(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                 tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,hia,
     &                 t1,z1,hz1,idiis,g_nt2,g_nz2,g_nhz2,
     &                 g_ncoul,g_nexch,tklst)
c
      call ga_distribution(g_nhz2,iam,jlo,jhi,ilo,ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.ilo.and.ad1+1.le.ihi)then
          do j=1,i
            ad2=(j-1)*nvir
            if (ad2+1.ge.jlo.and.ad2+1.le.jhi)then
              call ga_get(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              call ga_get(g_nhz2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
              do a=1,nvir
              do b=1,nvir
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                scra(ad3)=scra(ad3)+scrb(ad4)
              enddo
              enddo
              call ga_put(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_sync()
c
      call ga_distribution(g_hz2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          do i=1,nocc
            ad1=(i-1)*nvir+a
            call ga_get(g_nhz2,1,lnov,ad1,ad1,scrb,lnov)
            do j=1,i
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                ad2=(j-1)*nvir+b
                scra(ad1)=scrb(ad2)
              enddo
            enddo
          enddo
          call ga_acc(g_hz2,offt2+1,offt2+lhoov,a,a,
     &                 scra,lhoov,1.0d00)
        endif
      enddo
c
      call ga_sync()
c
c
c ----------------------------------------------------------
c symmetrize the sigma vector and evaluated error vectors
c ----------------------------------------------------------
c
      offt2=(idiis-1)*lhoov
      nvp1=nvir+1
c
      if (iprt.gt.5)print *,'start symmetrization'
      icnt=0
      do a=ncor+nocc+1,ncor+nocc+nvir
       do i=ncor+1,ncor+nocc
        icnt=icnt+1
        hz1(icnt)=-hz1(icnt)/(eorb(i)-eorb(a))
       enddo
      enddo
      call daxpy(lnov,-one,hz1,1,z1,1)
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_z2,offt2+1,offt2+lnov,nvp1,nvp1,hz1,1)
       call ga_put(g_hz2,offt2+1,offt2+lnov,nvp1,nvp1,z1,1)
      endif
c
      call ga_sync()
c
      do b=ncor+nocc+1,ncor+nocc+nvir
       bv=b-nocc-ncor
c
c -- parallel task --
       if (bv.ge.ilo.and.bv.le.ihi)then
        call ga_get(g_z2,offt2+1,offt2+lhoov,bv,bv,scra,lhoov)
        call ga_get(g_hz2,offt2+1,offt2+lhoov,bv,bv,scrb,lhoov)
        icnt=0
        do a=ncor+nocc+1,ncor+nocc+nvir
         do i=ncor+1,ncor+nocc
          do j=ncor+1,i
           icnt=icnt+1
           scrb(icnt)=-scrb(icnt)/(eorb(i)+eorb(j)-eorb(a)-eorb(b))
          enddo
         enddo
        enddo
c -- compute error vector in z2 location
        call daxpy(lhoov,-one,scrb,1,scra,1)
        call ga_put(g_z2,offt2+1,offt2+lhoov,bv,bv,scrb,lhoov)
        call ga_put(g_hz2,offt2+1,offt2+lhoov,bv,bv,scra,lhoov)
       endif
c -- end parallel task --
c
      enddo
      call ga_sync()
c
c     if (iprt.gt.50.and.iam.eq.0)then
      if (iam.eq.0.and.oprint)then
       print *,' sym hz1 '
       write(6,189)(hz1(i),i=1,lnov)
       do bv=1,nvir
        call ga_get(g_z2,offt2+1,offt2+lhoov,bv,bv,scra,1)
        print *,'sym sig2 for bv ',bv
        write(6,189)(scra(i),i=1,lhoov)
  189   format(4e14.4)
       enddo
      endif
c
      call ga_sync()
c
      if (occsdps) then
         call pstat_off(ps_zvecs)
      else
         call qexit('zvecs',0)
      endif
      return
      end
