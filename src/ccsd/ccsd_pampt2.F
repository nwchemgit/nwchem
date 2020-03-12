      subroutine ccsd_pampt2(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     $                      iprt,cmo,eorb,t1,ht1,
     $                      scra,scrb,hiu,hia,habe,giu,gabe,
     $                      g_t2,g_ht2,g_nt2,g_nht2,
     $                      eccsd,g_ncoul,g_nexch,
     &                      tklst, Tol2E, occd, eccsdos,use_ccsd_omp)
C     $Id$
      implicit none
      integer basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,iprt,
     &        g_t2,g_ht2,g_ncoul,g_nexch,g_nt2,g_nht2,
     &        tklst(nsh*(nsh+1)/2,2)
c
#include "ccsd_len.fh"
      double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                 ht1(*),scra(*),scrb(*),
     &                 hiu(*),hia(*),habe(*),giu(*),gabe(*)
      double precision tol2e
      logical occd
      logical, optional, intent(in) :: use_ccsd_omp
c
      integer nodes,iam
#include "mafdecls.fh"
#include "ccsd_debug.fh"
#include "bas.fh"
#include "msgids.fh"
#include "global.fh"
#include "ccsdps.fh"
c
      integer i,j,a,b,av,bv,ij,ad1,ad2,ad3,ad4,offt2,icnt,nvp1,
     &        ilo,ihi,jlo,jhi,lhoo,lhoov
      integer iin,jin
      double precision eccsd,fac, eccsdos
#ifdef USE_F90INTERFACE
      interface
      subroutine ccsd_pampt3(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                       tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,
     &                       hia,t1,ht1,idiis,g_nt2,g_nht2,
     &                       g_ncoul,g_nexch,tklst,use_ccsd_omp)
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,
     &        iprt,iam,nodes,idiis,g_nt2,g_nht2,
     &        g_ncoul,g_nexch,tklst(*)

      double precision tol2e,cmo(nbf,*),eorb(nbf),
     &                 hiu(*),giu(*),habe(*),gabe(*),hia(*),
     &                 t1(*),ht1(*)
      logical, optional, intent(in) :: use_ccsd_omp
      end subroutine
      end interface
#endif         
c
      double precision zip,one,half,two
c
      data zip/0.0d00/,one/1.0d00/,half/0.5d00/,two/2.0d00/
c
      if (.not.present(use_ccsd_omp)) then
          call errquit('ccsd_pampt2: use_ccsd_omp not present!',0,0)
      endif
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
      call ga_sync()
      if (occsdps) then
         call pstat_on(ps_pampt)
      else
         call qenter('pampt',0)
      endif
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
      lnov=nocc*nvir
      lnvv=nvir*nvir
      offt2=(idiis-1)*lhoov
c
#if 0
      call ga_get(g_t2,offt2+1,offt2+lnov,nvir+1,nvir+1,
     &            t1,lnov)
#else
      if(ga_nodeid().eq.0)
     G      call ga_get(g_t2,offt2+1,offt2+lnov,nvir+1,nvir+1,
     &            t1,lnov)
      call ga_brdcst(65537, t1, 
     G     MA_sizeof(MT_DBL,1,MT_BYTE)*lnov, 0)
#endif
c
c ------------------------------------------------------------
c fill square GAs
c ------------------------------------------------------------
c
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          call ga_get(g_t2,offt2+1,offt2+lhoov,a,a,
     &                 scra,lhoov)
cstaggering to avoid congestion in ga_put
         do iin=ga_nodeid()+1,ga_nodeid()+nocc
            i=mod(iin,nocc)+1
            call dfill(nvir,0.0d00,scrb,1)
            do jin=ga_nodeid()+1,ga_nodeid()+i
               j=mod(jin,i)+1
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                scrb(b)=scra(ad1)
              enddo
              ad1=(i-1)*nvir+a
              call ga_put(g_nt2,(j-1)*nvir+1,j*nvir,ad1,ad1,scrb,nvir)
            enddo
          enddo
        endif
      enddo
c
      call ga_sync()
c
      call ga_distribution(g_nt2,iam,jlo,jhi,ilo,ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.ilo.and.ad1+1.le.ihi)then
cstagger          do j=1,i-1
          do jin=ga_nodeid()+1,ga_nodeid()+i-1
             j=mod(jin,i-1)+1
            ad2=(j-1)*nvir
            if (ad2+1.ge.jlo.and.ad2+1.le.jhi)then
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              do a=1,nvir
              do b=1,nvir
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                scrb(ad4)=scra(ad3)
              enddo
              enddo
              call ga_put(g_nt2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_zero(g_nht2)
      call ccsd_pampt3(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                 tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,hia,
     &                 t1,ht1,idiis,g_nt2,g_nht2,
     &                 g_ncoul,g_nexch,tklst,use_ccsd_omp)
      if (iam.eq.0.and.iprt.gt.50)then
       print *,'after pampt3 sym ht1 '
       write(6,189)(ht1(i),i=1,lnov)
       do bv=1,lnov
        call ga_get(g_nht2,1,lnov,bv,bv,scra,1)
        print *,'after pampt3 sym sig2 for bv ',bv
        write(6,189)(scra(i),i=1,lnov)
       enddo
      endif
      call ga_sync()
c
      call ga_distribution(g_nht2,iam,jlo,jhi,ilo,ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.ilo.and.ad1+1.le.ihi)then
          do jin=ga_nodeid()+1,ga_nodeid()+i
             j=mod(jin,i)+1
            ad2=(j-1)*nvir
            if (ad2+1.ge.jlo.and.ad2+1.le.jhi)then
              call ga_get(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              call ga_get(g_nht2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
              do a=1,nvir
              do b=1,nvir
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                scra(ad3)=scra(ad3)+scrb(ad4)
              enddo
              enddo
              call ga_put(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_sync()
      call ga_distribution(g_ht2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          do iin=1,nocc
             i=mod(iin,nocc)+1
            ad1=(i-1)*nvir+a
            call ga_get(g_nht2,1,lnov,ad1,ad1,scrb,lnov)
            do j=1,i
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                ad2=(j-1)*nvir+b
                scra(ad1)=scrb(ad2)
              enddo
            enddo
          enddo
          call ga_acc(g_ht2,offt2+1,offt2+lhoov,a,a,
     &                 scra,lhoov,1.0d00)
        endif
      enddo
      call ga_sync()
c
c
c ----------------------------------------------------------
c symmetrize the sigma vector and evaluated error vectors
c ----------------------------------------------------------
c
      offt2=(idiis-1)*lhoov
      nvp1=nvir+1

*     RJH
      if (occd) call dfill(lnov, 0.0d0, ht1, 1)
*
c
      if (iam.eq.0.and.iprt.gt.5)print *,'start symmetrization'
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
      icnt=0
      do a=ncor+nocc+1,ncor+nocc+nvir
       do i=ncor+1,ncor+nocc
        icnt=icnt+1
        ht1(icnt)=ht1(icnt)/(eorb(i)-eorb(a))
       enddo
      enddo
      call daxpy(lnov,-one,ht1,1,t1,1)
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,ht1,lnov)
       call ga_put(g_ht2,offt2+1,offt2+lnov,nvp1,nvp1,t1,lnov)
      endif
      call ga_sync()
c
      if (iam.eq.0.and.iprt.gt.50)then
       print *,' sym ht1 '
       write(6,189)(ht1(i),i=1,lnov)
       do bv=1,nvir
        call ga_get(g_ht2,offt2+1,offt2+lhoov,bv,bv,scra,lhoov)
        print *,'sym sig2 for bv ',bv
        write(6,189)(scra(i),i=1,lhoov)
  189   format(4e14.4)
       enddo
      endif
      call ga_sync()
c
      do b=ncor+nocc+1,ncor+nocc+nvir
       bv=b-nocc-ncor
c
c -- parallel task --
       if (bv.ge.ilo.and.bv.le.ihi)then
        call ga_get(g_t2,offt2+1,offt2+lhoov,bv,bv,scra,lhoov)
        call ga_get(g_ht2,offt2+1,offt2+lhoov,bv,bv,scrb,lhoov)
        icnt=0
        do a=ncor+nocc+1,ncor+nocc+nvir
         do i=ncor+1,ncor+nocc
          do j=ncor+1,i
           icnt=icnt+1
           scrb(icnt)=scrb(icnt)/(eorb(i)+eorb(j)-eorb(a)-eorb(b))
          enddo
         enddo
        enddo
c -- compute error vector in t2 location
        call daxpy(lhoov,-one,scrb,1,scra,1)
        call ga_put(g_t2,offt2+1,offt2+lhoov,bv,bv,scrb,lhoov)
        call ga_put(g_ht2,offt2+1,offt2+lhoov,bv,bv,scra,lhoov)
       endif
c -- end parallel task --
c
      enddo
      call ga_sync()
      if (iam.eq.0.and.iprt.gt.50)then
       do bv=1,nvir
        call ga_get(g_t2,offt2+1,offt2+lhoov,bv,bv,scra,lhoov)
        print *,'sym sig2 for bv ',bv
        write(6,189)(scra(i),i=1,lhoov)
       enddo
      endif
      call ga_sync()
c
c
c ----------------------------------------------------------
c compute energy
c ----------------------------------------------------------
c
      if (iam.eq.0.and.iprt.gt.5)print *,'start energy'
      eccsd=zip
      eccsdos=zip
      call ga_zero(g_nt2)
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          call ga_get(g_t2,offt2+1,offt2+lhoov,a,a,
     &                 scra,lhoov)
          do i=1,nocc
            do j=1,i
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                ad2=(j-1)*nvir+b
                scrb(ad2)=scra(ad1)
              enddo
            enddo
            ad1=(i-1)*nvir+a
            call ga_put(g_nt2,1,i*nvir,ad1,ad1,scrb,i*nvir)
          enddo
        endif
      enddo
c
      call ga_sync()
      call ga_distribution(g_nexch,iam,jlo,jhi,ilo,ihi)
      do i=1,nocc
        ad1=(i-1)*nact
        if (ad1+1.ge.ilo.and.ad1+1.le.ihi)then
         do j=1,i
          ad2=(j-1)*nact
          if (ad2+1.ge.jlo.and.ad2+1.le.jhi)then
            ij=(i-1)*nocc+j
            fac=one
            if(i.eq.j)fac=half
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+nocc+1,ad1+nact,
     &                  scra,nvir)
            call ga_get(g_nt2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,scrb,nvir)
            do av=1,nvir
            do bv=1,nvir
              ad3=(av-1)*nvir+bv
              scrb(ad3)=scrb(ad3)+ht1((av-1)*nocc+i)*ht1((bv-1)*nocc+j)
            enddo
            enddo
            do av=1,nvir
            do bv=1,nvir
              ad3=(av-1)*nvir+bv
              ad4=(bv-1)*nvir+av
              eccsd=eccsd+fac*(two*scrb(ad3)-scrb(ad4))*scra(ad3)
              eccsd=eccsd+fac*(two*scrb(ad4)-scrb(ad3))*scra(ad4)

c             opposite spin component of ccsd energy for spin component scaled CCSD

              eccsdos=eccsdos+fac*scra(ad3)*scrb(ad3)
              eccsdos=eccsdos+fac*scra(ad4)*scrb(ad4)
            enddo
            enddo
          endif
         enddo
        endif
      enddo
      call ga_dgop(msg_cc_ecc,eccsd,1, '+')
      call ga_dgop(msg_cc_ecc,eccsdos,1, '+')
      if (iprt.gt.5.and.iam.eq.0)then
       write(6,1238)eccsd
 1238  format(/, ' eccsd  ',e20.12,/)
      endif
c
      if (occsdps) then
         call pstat_off(ps_pampt)
      else
         call qexit('pampt',0)
      endif
      return
      end
