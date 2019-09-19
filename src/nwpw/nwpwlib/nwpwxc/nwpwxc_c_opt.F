C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_opt.F
C> The OPT correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The OPT correlation functional
C>
C> The OPTimized correlation functional [1].
C>
C> [1] N.C. Handy, A.J. Cohen, "Dynamic correlation",
C>     Mol. Phys. <b>99</b>, 607-615 (2001), DOI:
C>     <A HREF="https://doi.org/10.1080/00268970010023435">
C>     10.1080/00268970010023435</A>.
c
C$Id$
c
      subroutine nwpwxc_c_opt(tol_rho,ipol,nq,wght,rho,rgamma,func,
     &                      Amat, Cmat)
      implicit none
c      
#include "nwpwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision func(nq)   !< [Output] The value of the functional
      double precision Amat(nq,*) !< [Output] The derivative wrt rho
      double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
c
c References:
c
c    Handy NC, Cohen AJ, Mol Phys 99 (7); 607-615 2001
c
      integer l_rho,k_rho,l_delrho,k_delrho
      double precision c1,c2
!      parameter (c1=1.1015d0,c2=0.6625d0)
      parameter (c1=1.d0,c2=0d0)
c
      integer iq   ! counter over grid points
      integer ii   ! offset
      integer iqt  ! upper limit
      integer num  ! the number of grid points in current batch
      integer maxp ! the maximum number of points in a batch
      parameter (maxp = 16)
      double precision funcl(maxp)
      double precision rhol(maxp,2)
      double precision rgammal(maxp,3)
      double precision Amatl(maxp,2)
      double precision Cmatl(maxp,3)
      double precision fac
c
c***************************************************************************
c
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do iq = 1, nq, maxp
           iqt = min(nq,iq+maxp-1)
           num = iqt-iq+1
c
c          c1*Ec[a,b]
c
           call dfill(maxp,0.0d0,funcl,1)
           call dfill(maxp*2,0.0d0,rhol,1)
           call dfill(maxp*3,0.0d0,rgammal,1)
           call dfill(maxp*2,0.0d0,Amatl,1)
           call dfill(maxp*3,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              rhol(ii+1,R_T) = rho(iq+ii,R_T)
              rgammal(ii+1,G_TT) = rgamma(iq+ii,G_TT)
           enddo
           fac = c1*wght
           call nwpwxc_c_pw91lda(tol_rho,1,maxp,fac,rhol,funcl,Amatl)
           call nwpwxc_c_p91(tol_rho,1,maxp,fac,rhol,rgammal,funcl,
     +                     Amatl,Cmatl)
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
              Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
              Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
     +                           + Cmatl(ii+1,D1_GAA)
              Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
     +                           + Cmatl(ii+1,D1_GAB)
           enddo
c
c          {(c2-c1)*Ec[a,0]}(a) + {(c2-c1)*Ec[0,b]}(a)
c
           fac = (c2-c1)*wght
           call dfill(maxp,0.0d0,funcl,1)
           call dfill(maxp*2,0.0d0,Amatl,1)
           call dfill(maxp*3,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              rhol(ii+1,R_A) = 0.5d0*rhol(ii+1,R_T)
              rgammal(ii+1,G_AA) = 0.25d0*rgammal(ii+1,G_TT)
           enddo
           call nwpwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl,Amatl)
           call nwpwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl,
     +                     Amatl,Cmatl)
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + 2*funcl(ii+1)
              Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA)
     +                          + Amatl(ii+1,D1_RA) + Amatl(ii+1,D1_RB)
              Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
     +                           + (2.0d0)*Cmatl(ii+1,D1_GAA) 
              Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
     +                           + (2.0d0)*Cmatl(ii+1,D1_GAB)
           enddo
         enddo
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do iq = 1, nq, maxp
           iqt = min(nq,iq+maxp-1)
           num = iqt-iq+1
c
c          c1*Ec[a,b]
c
           call dfill(maxp,0.0d0,funcl,1)
           call dfill(maxp*2,0.0d0,rhol,1)
           call dfill(maxp*3,0.0d0,rgammal,1)
           call dfill(maxp*2,0.0d0,Amatl,1)
           call dfill(maxp*3,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              rhol(ii+1,R_A) = rho(iq+ii,R_A)
              rhol(ii+1,R_B) = rho(iq+ii,R_B)
              rgammal(ii+1,G_AA) = rgamma(iq+ii,G_AA)
              rgammal(ii+1,G_AB) = rgamma(iq+ii,G_AB)
              rgammal(ii+1,G_BB) = rgamma(iq+ii,G_BB)
           enddo
           fac = c1*wght
           call nwpwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl,Amatl)
           call nwpwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl,
     +                     Amatl,Cmatl)
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
              Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
              Amat(iq+ii,D1_RB) = Amat(iq+ii,D1_RB) + Amatl(ii+1,D1_RB)
              Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
     +                           + Cmatl(ii+1,D1_GAA)
              Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
     +                           + Cmatl(ii+1,D1_GAB)
              Cmat(iq+ii,D1_GBB) = Cmat(iq+ii,D1_GBB)
     +                           + Cmatl(ii+1,D1_GBB)
           enddo
c
c          (c2-c1)*Ec[a,0]
c
           fac = (c2-c1)*wght
           call dfill(maxp  ,0.0d0,funcl,1)
           call dfill(maxp*2,0.0d0,Amatl,1)
           call dfill(maxp*3,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              rhol(ii+1,R_A) = rho(iq+ii,R_A)
              rhol(ii+1,R_B) = 0.0d0
              rgammal(ii+1,G_AA) = rgamma(iq+ii,G_AA)
              rgammal(ii+1,G_AB) = 0.0d0
              rgammal(ii+1,G_BB) = 0.0d0
           enddo
           call nwpwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl,Amatl)
           call nwpwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl,
     +                     Amatl,Cmatl)
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
              Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
              Amat(iq+ii,D1_RB) = Amat(iq+ii,D1_RB) + Amatl(ii+1,D1_RB)
              Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
     +                           + Cmatl(ii+1,D1_GAA)
              Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
     +                           + Cmatl(ii+1,D1_GAB)
              Cmat(iq+ii,D1_GBB) = Cmat(iq+ii,D1_GBB)
     +                           + Cmatl(ii+1,D1_GBB)
           enddo
c
c          (c2-c1)*Ec[0,b]
c
           fac = (c2-c1)*wght
           call dfill(maxp  ,0.0d0,funcl,1)
           call dfill(maxp*2,0.0d0,Amatl,1)
           call dfill(maxp*3,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              rhol(ii+1,R_A) = 0.0d0
              rhol(ii+1,R_B) = rho(iq+ii,R_B)
              rgammal(ii+1,G_AA) = 0.0d0
              rgammal(ii+1,G_AB) = 0.0d0
              rgammal(ii+1,G_BB) = rgamma(iq+ii,G_BB)
           enddo
           call nwpwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl,Amatl)
           call nwpwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl,
     +                     Amatl,Cmatl)
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
              Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
              Amat(iq+ii,D1_RB) = Amat(iq+ii,D1_RB) + Amatl(ii+1,D1_RB)
              Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
     +                           + Cmatl(ii+1,D1_GAA)
              Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
     +                           + Cmatl(ii+1,D1_GAB)
              Cmat(iq+ii,D1_GBB) = Cmat(iq+ii,D1_GBB)
     +                           + Cmatl(ii+1,D1_GBB)
           enddo
         enddo
c
      endif
c
      return
      end
C>
C> @}
