#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_opt.F
C> The OPT correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The OPT correlation functional
C>
C> The OPTimized correlation functional [1]. This functional uses
C> the Stoll partitioning of the correlation energy [2]. This is 
C> formulated as (see Eq.(3) in [2]):
C> \f{eqnarray*}
C>   E_c \approx \int (\rho_\alpha+\rho_\beta)
C>       \epsilon_c(\rho_\alpha,\rho_\beta)dr
C>      -  \int\rho_\alpha\epsilon_c(\rho_\alpha,0)dr
C>      -  \int\rho_\beta\epsilon_c(0,\rho_\beta)dr
C> \f}
C> However, this equation can be interpreted in two ways. The second
C> and third terms can be interpreted as the full functional evaluated
C> at zero density for one of the spin components. Alternatively,
C> they can be interpreted as functionals of one component of the 
C> spin density only (and not a functional of the other spin component
C> at all). The text in [2] clearly states that the latter 
C> interpretation is the correct one. In practice this means that for
C> the second term all derivatives with respect to \f$\rho_\beta,
C> \gamma_{\alpha\beta}\f$ and \f$\gamma_{\beta\beta}\f$ are zero.
C> Similar consequences follow for the third term.
C>
C> [1] N.C. Handy, A.J. Cohen, "Dynamic correlation",
C>     Mol. Phys. <b>99</b>, 607-615 (2001), DOI:
C>     <A HREF="https://doi.org/10.1080/00268970010023435">
C>     10.1080/00268970010023435</A>.
C>
C> [2] H. Stoll, C.M.E. Pavlidou, H. Preu&szlig;,
C>     "On the calculation of correlation energies in the spin-density
C>     functional formalism", Theor. Chim. Acta <b>49</b>, 143-149
C>     (1978), DOI:
C>     <A HREF="https://doi.org/10.1007/PL00020511">
C>     10.1007/PL00020511</A>.
c
C$Id$
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      subroutine nwxc_c_opt_p(tol_rho,ipol,nq,wght,rho,rgamma,func)
#else
      subroutine nwxc_c_opt(tol_rho, ipol, nq, wght, rho, rgamma, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_c_opt_d2(tol_rho,ipol,nq,wght,rho,rgamma,func)
#else
      subroutine nwxc_c_opt_d3(tol_rho,ipol,nq,wght,rho,rgamma,func)
#endif
c      
#include "nwad.fh"
c
      implicit none
c      
#include "nwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      type(nwad_dble)::func(nq)   !< [Output] The value of the functional
c     double precision Amat(nq,*) !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
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
      type(nwad_dble)::funcl(maxp)
      type(nwad_dble)::rhol(maxp,2)
      type(nwad_dble)::rgammal(maxp,3)
c     double precision Amatl(maxp,2)
c     double precision Cmatl(maxp,3)
      double precision fac
      double precision tmp(2)
      integer lennw
c
c***************************************************************************
c
c
      lennw = loc(tmp(2))-loc(tmp(1))
      lennw = (loc(funcl(2))-loc(funcl(1))+lennw-1)/lennw
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
c          call dfill(maxp*lennw,0.0d0,funcl,1)
c          call dfill(maxp*2*lennw,0.0d0,rhol,1)
c          call dfill(maxp*3*lennw,0.0d0,rgammal,1)
c          call dfill(maxp*2*lennw,0.0d0,Amatl,1)
c          call dfill(maxp*3*lennw,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              funcl(ii+1)        = 0.0d0
              rhol(ii+1,R_T)     = rho(iq+ii,R_T)
              rgammal(ii+1,G_TT) = rgamma(iq+ii,G_TT)
           enddo
           fac = c1*wght
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_c_pw91lda_p(tol_rho,1,maxp,fac,rhol,funcl)
           call nwxc_c_p91_p(tol_rho,1,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda(tol_rho,1,maxp,fac,rhol,funcl)
           call nwxc_c_p91(tol_rho,1,maxp,fac,rhol,rgammal,funcl)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_c_pw91lda_d2(tol_rho,1,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d2(tol_rho,1,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda_d3(tol_rho,1,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d3(tol_rho,1,maxp,fac,rhol,rgammal,funcl)
#endif
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
c             Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
c             Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
c    +                           + Cmatl(ii+1,D1_GAA)
c             Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
c    +                           + Cmatl(ii+1,D1_GAB)
           enddo
c
c          {(c2-c1)*Ec[a,0]}(a) + {(c2-c1)*Ec[0,b]}(a)
c
           fac = (c2-c1)*wght
c          call dfill(maxp*lennw,0.0d0,funcl,1)
c          call dfill(maxp*2*lennw,0.0d0,Amatl,1)
c          call dfill(maxp*3*lennw,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              funcl(ii+1)        = 0.0d0
              rhol(ii+1,R_A)     = 0.5d0*rhol(ii+1,R_T)
              rhol(ii+1,R_B)     = 0.0d0
              rgammal(ii+1,G_AA) = 0.25d0*rgammal(ii+1,G_TT)
              rgammal(ii+1,G_AB) = 0.0d0
              rgammal(ii+1,G_BB) = 0.0d0
           enddo
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_c_pw91lda_p(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_p(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_c_pw91lda_d2(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d2(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda_d3(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d3(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + 2.0d0*funcl(ii+1)
c             Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA)
c    +                          + Amatl(ii+1,D1_RA) + Amatl(ii+1,D1_RB)
c             Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
c    +                           + (2.0d0)*Cmatl(ii+1,D1_GAA) 
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
c          call dfill(maxp*lennw,0.0d0,funcl,1)
c          call dfill(maxp*2*lennw,0.0d0,rhol,1)
c          call dfill(maxp*3*lennw,0.0d0,rgammal,1)
c          call dfill(maxp*2,0.0d0,Amatl,1)
c          call dfill(maxp*3,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              funcl(ii+1)        = 0.0d0
              rhol(ii+1,R_A)     = 0.0d0
              rhol(ii+1,R_B)     = 0.0d0
              rgammal(ii+1,G_AA) = 0.0d0
              rgammal(ii+1,G_AB) = 0.0d0
              rgammal(ii+1,G_BB) = 0.0d0

              if (rho(iq+ii,R_A).gt.0.5d0*tol_rho) then
                rhol(ii+1,R_A)     = rho(iq+ii,R_A)
                rgammal(ii+1,G_AA) = rgamma(iq+ii,G_AA)
              endif
              if (rho(iq+ii,R_B).gt.0.5d0*tol_rho) then
                rhol(ii+1,R_B)     = rho(iq+ii,R_B)
                rgammal(ii+1,G_BB) = rgamma(iq+ii,G_BB)
              endif
              if ((rho(iq+ii,R_A).gt.0.5d0*tol_rho).and.
     &            (rho(iq+ii,R_B).gt.0.5d0*tol_rho)) then
                rgammal(ii+1,G_AB) = rgamma(iq+ii,G_AB)
              endif
           enddo
           fac = c1*wght
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_c_pw91lda_p(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_p(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_c_pw91lda_d2(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d2(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda_d3(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d3(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
c             Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
c             Amat(iq+ii,D1_RB) = Amat(iq+ii,D1_RB) + Amatl(ii+1,D1_RB)
c             Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
c    +                           + Cmatl(ii+1,D1_GAA)
c             Cmat(iq+ii,D1_GAB) = Cmat(iq+ii,D1_GAB)
c    +                           + Cmatl(ii+1,D1_GAB)
c             Cmat(iq+ii,D1_GBB) = Cmat(iq+ii,D1_GBB)
c    +                           + Cmatl(ii+1,D1_GBB)
           enddo
c
c          (c2-c1)*Ec[a,0]
c
           fac = (c2-c1)*wght
c          call dfill(maxp*lennw  ,0.0d0,funcl,1)
c          call dfill(maxp*2*lennw,0.0d0,Amatl,1)
c          call dfill(maxp*3*lennw,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              funcl(ii+1)        = 0.0d0
              rhol(ii+1,R_A)     = 0.0d0
              rhol(ii+1,R_B)     = 0.0d0
              rgammal(ii+1,G_AA) = 0.0d0
              rgammal(ii+1,G_AB) = 0.0d0
              rgammal(ii+1,G_BB) = 0.0d0
              if (rho(iq+ii,R_A).gt.0.5d0*tol_rho) then
                rhol(ii+1,R_A)     = rho(iq+ii,R_A)
                rgammal(ii+1,G_AA) = rgamma(iq+ii,G_AA)
              endif
           enddo
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_c_pw91lda_p(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_p(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_c_pw91lda_d2(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d2(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda_d3(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d3(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
c             Amat(iq+ii,D1_RA) = Amat(iq+ii,D1_RA) + Amatl(ii+1,D1_RA)
c             Cmat(iq+ii,D1_GAA) = Cmat(iq+ii,D1_GAA)
c    +                           + Cmatl(ii+1,D1_GAA)
           enddo
c
c          (c2-c1)*Ec[0,b]
c
           fac = (c2-c1)*wght
c          call dfill(maxp*lennw  ,0.0d0,funcl,1)
c          call dfill(maxp*2*lennw,0.0d0,Amatl,1)
c          call dfill(maxp*3*lennw,0.0d0,Cmatl,1)
           do ii = 0, num - 1
              funcl(ii+1)        = 0.0d0
              rhol(ii+1,R_A)     = 0.0d0
              rhol(ii+1,R_B)     = 0.0d0
              rgammal(ii+1,G_AA) = 0.0d0
              rgammal(ii+1,G_AB) = 0.0d0
              rgammal(ii+1,G_BB) = 0.0d0
              if (rho(iq+ii,R_B).gt.0.5d0*tol_rho) then
                rhol(ii+1,R_B)     = rho(iq+ii,R_B)
                rgammal(ii+1,G_BB) = rgamma(iq+ii,G_BB)
              endif
           enddo
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_c_pw91lda_p(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_p(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_c_pw91lda_d2(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d2(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#else
           call nwxc_c_pw91lda_d3(tol_rho,2,maxp,fac,rhol,funcl)
           call nwxc_c_p91_d3(tol_rho,2,maxp,fac,rhol,rgammal,funcl)
#endif
           do ii = 0, num - 1
              func(iq+ii) = func(iq+ii) + funcl(ii+1)
c             Amat(iq+ii,D1_RB) = Amat(iq+ii,D1_RB) + Amatl(ii+1,D1_RB)
c             Cmat(iq+ii,D1_GBB) = Cmat(iq+ii,D1_GBB)
c    +                           + Cmatl(ii+1,D1_GBB)
           enddo
         enddo
c
      endif
c
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_opt.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_opt.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_opt.F"
#endif
#undef NWAD_PRINT
C>
C> @}
