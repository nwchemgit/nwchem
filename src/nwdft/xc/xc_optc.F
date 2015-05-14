c
C$Id$
c
      subroutine xc_optc(rho, delrho, 
     &                      Amat, Cmat, nq, Ec, qwght,ldew,func)
      implicit none
#include "errquit.fh"
c      
#include "dft2drv.fh"
#include "cdft.fh"
#include "mafdecls.fh"
c
      double precision Ec
      integer nq
      logical ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c
c References:
c
c    Handy NC, Cohen AJ, Mol Phys 99 (7); 607-615 2001
c
      integer l_rho,k_rho,l_delrho,k_delrho
      integer l_amat, k_amat, l_cmat, k_cmat
      double precision c1,c2
!      parameter (c1=1.1015d0,c2=0.6625d0)
      parameter (c1=1.d0,c2=0d0)
c
c***************************************************************************
c
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
c
         call errquit(' optc closed shell not coded ',0, INPUT_ERR)
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         
      lcfac(5)=.true.
      nlcfac(5)=.true.
      cfac(6)=1d0
c
c     assuming cfac(5) set to 1 ...
c     
         cfac(5)=c1
c
c     c1*Ec[a,b]
c
            call xc_perdew91(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Cmat, nq, 2, 
     &           Ec, qwght, ldew, func)
            if (.not. ma_push_get(mt_dbl, nq*3, 
     .           ' rho aux', l_rho, k_rho))
     .           call errquit(' xcoptc: push_get failed ',0, MA_ERR)
            if (.not. ma_push_get(mt_dbl, nq*6, 
     .           ' rho aux', l_delrho, k_delrho))
     .           call errquit(' xcoptc: pus_get failed ',1, MA_ERR)
            if (.not. ma_push_get(mt_dbl, nq*2, 
     .           ' amat aux', l_amat, k_amat))
     .           call errquit(' xcoptc: push_get failed ',2, MA_ERR)
            if (.not. ma_push_get(mt_dbl, nq*3, 
     .           ' cmat aux', l_cmat, k_cmat))
     .           call errquit(' xcoptc: push_get failed ',3, MA_ERR)
            call dcopy(nq*3,0.d0,0,dbl_mb(k_rho),1)
            call dcopy(nq*6,0.d0,0,dbl_mb(k_delrho),1)
            call dcopy(nq*2,0.d0,0,dbl_mb(k_amat),1)
            call dcopy(nq*3,0.d0,0,dbl_mb(k_cmat),1)
            call  dcopy(nq,rho(1,2),1,dbl_mb(k_rho),1)
            call  dcopy(nq,rho(1,2),1,dbl_mb(k_rho+nq),1)
            call  dcopy(nq*3,delrho(1,1,1),1,dbl_mb(k_delrho),1)
         cfac(5)=c2-c1
c
c     (c2-c1)*Ec[a,0]
c
            call xc_perdew91(tol_rho, cfac, lcfac, nlcfac, 
     .        dbl_mb(k_rho),  dbl_mb(k_delrho), 
     &        dbl_mb(k_amat), dbl_mb(k_cmat), nq, 2, 
     &           Ec, qwght, ldew, func)
c
            call daxpy(nq,1.0d0,dbl_mb(k_amat+(D1_RA-1)*nq),1,
     &                 amat(1,D1_RA),1)
            call daxpy(nq,1.0d0,dbl_mb(k_cmat+(D1_GAA-1)*nq),1,
     &                 cmat(1,D1_GAA),1)
c
            call dcopy(nq*3,0.d0,0,dbl_mb(k_rho),1)
            call dcopy(nq*6,0.d0,0,dbl_mb(k_delrho),1)
            call dcopy(nq*2,0.d0,0,dbl_mb(k_amat),1)
            call dcopy(nq*3,0.d0,0,dbl_mb(k_cmat),1)
            call  dcopy(nq,rho(1,3),1,dbl_mb(k_rho),1)
            call  dcopy(nq,rho(1,3),1,dbl_mb(k_rho+2*nq),1)
            call  dcopy(nq*3,delrho(1,1,2),1,dbl_mb(k_delrho+3*nq),1)
         cfac(5)=c2-c1
c
c     (c2-c1)*Ec[0,b]
c
            call xc_perdew91(tol_rho, cfac, lcfac, nlcfac, 
     .        dbl_mb(k_rho),  dbl_mb(k_delrho), 
     &        dbl_mb(k_amat), dbl_mb(k_cmat), nq, 2, 
     &           Ec, qwght, ldew, func)
c
            call daxpy(nq,1.0d0,dbl_mb(k_amat+(D1_RB-1)*nq),1,
     &                 amat(1,D1_RB),1)
            call daxpy(nq,1.0d0,dbl_mb(k_cmat+(D1_GBB-1)*nq),1,
     &                 cmat(1,D1_GBB),1)
c
            if(.not.ma_chop_stack(l_rho)) call
     .           errquit(' xcoptc: chop_stack failed',0, MA_ERR)
c
      endif
c
      lcfac(5)=.false.
      nlcfac(5)=.false.
      cfac(5)=0d0
      cfac(6)=0d0
      return
      end
