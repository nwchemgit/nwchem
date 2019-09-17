#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_op.F
C> The OP correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the OP correlation functional
C>
C> The OP correlation functional [1,2] is a functional designed to
C> have few optimized parameters (only one in this case) and has
C> a "progressive" form. 
C>
C> ### References ###
C>
C> [1] T. Tsuneda, T. Suzumura, K. Hirao,
C>     "A new one-parameter progressive Colle-Salvetti-type correlation
C>     functional", J. Chem. Phys. <b>110</b>, 10664-10678 (1999), DOI:
C>     <a href="https://doi.org/10.1063/1.479012">
C>     10.1063/1.479012</a>.
C>
C> [2] T. Tsuneda, T. Suzumura, K. Hirao,
C>     "A reexamination of exchange energy functionals",
C>     J. Chem. Phys. <b>111</b>, 5656-5667 (1999), DOI:
C>     <a href="https://doi.org/10.1063/1.479954">
C>     10.1063/1.479954</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_op_p(kop,param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                     func)
#else
      Subroutine nwxc_c_op(kop,param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                     func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_op_d2(kop,param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                     func)
#else
      Subroutine nwxc_c_op_d3(kop,param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                     func)
#endif
c
C$Id$
c
#include "nwad.fh"
c
      implicit none
c      
#include "nwxc_param.fh"
c
c     Input and other parameters
c
      external kop             !< [Input] Subroutine to evaluate the
      !< GGA exchange enhancement factor (see Eq.(14) in [1]).
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)!< [Input] Parameters of functional
      type(nwad_dble)::QABOP
#else
      double precision param(*)!< [Input] Parameters of functional
      double precision QABOP
#endif
#else
      double precision param(*)!< [Input] Parameters of functional
      !< - param(1): \f$ q^{\alpha\beta}_{OP} \f$ see Eq.(27) in [1] and
      !<   Table III in [2].
      double precision QABOP
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      type(nwad_dble)::rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      type(nwad_dble)::func(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c
      double precision QAB88OP,QABPBOP
c     Parameter (QAB88OP=2.3670D0,QABPBOP=2.3789D0)
c
c References:
c    Tsuneda, Suzumura, Hirao, JCP 110, 10664 (1999)
c    Tsuneda, Suzumura, Hirao, JCP 111, 5656 (1999)
c
c***************************************************************************
c
      integer n
      type(nwad_dble)::rho13, rho43, gamma, x
      type(nwad_dble)::kalpha,kbeta, rho13a, rho13b,rhoa,rhob
c     type(nwad_dble)::banb, hbab, hbabx
      type(nwad_dble)::banb, hbab
      double precision dhdab,dhdabx,dkadra,dkbdrb,dkadxa,dkbdxb,
     A     dbabdra,dbabdrb,dbabdga,dbabdgb,dkadga,dkbdgb,
     A     dbabdka,dbabdkb
c
c     hbabx(x) = (1.5214d0*x + 0.5764d0)/
c    /           (x**2.0d0*(x**2.0d0+1.1284d0*x+0.3183d0))
c     dhdabx(x) = -(4.5642d0*x**4+5.7391d0*x**3+
c    +     2.4355*x**2+0.3669d0*x)/
c    /           ((x**4+1.1284d0*x**3+0.3183d0*x**2)**2)
c
c     if(whichf.eq.'be88') then
c        QABOP=QAB88OP
c     endif
c     if(whichf.eq.'pb96') then
c        QABOP=QABPBOP
c     endif
      QABOP = param(1)
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rhoa=rho(n,R_T)*0.5d0
            rho13a = (rhoa)**(1.d0/3.d0)
            rho43 = rho13a**4.0d0
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = 0.25d0 * gamma
            if (sqrt(gamma).gt.tol_rho)then
               x = sqrt(gamma) / rho43
               call kop(tol_rho,x,kalpha)
c              dkadra = -(4d0/3d0)*x*dkadxa/rhoa
c              dkadga = (dkadxa/rho43)*0.5d0/sqrt(gamma)
            else
               x=0d0
               call kop(tol_rho,x,kalpha)
c              dkadra = 0d0
c              dkadga = 0d0
            endif
c
c
          
            banb = qabop * rho13a * kalpha *0.5d0

            if(banb.ne.0.0d0) then
c              dbabdra = banb*0.5d0*
c    /              (1d0/(3d0*rhoa)+dkadra/kalpha)

c              dbabdga = banb/kalpha*dkadga*0.5d0
               
               hbab = hbabx(banb)
c              dhdab = dhdabx(banb)
            else
c              dbabdra =0d0
c              dbabdga =0d0
               hbab = 0d0
c              dhdab = 0d0
            endif

c           Ec = Ec - rhoa**2*hbab*qwght(n)*fac
            func(n) = func(n) - rhoa**2.0d0*hbab*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) - 
c    -           (rhoa*hbab + rhoa**2*dhdab*dbabdra)*wght

c
c           if (x.gt.tol_rho) then
c               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - 
c    -              rhoa**2*dhdab*dbabdga*wght
c            endif
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
            if (rho(n,R_A).ge.tol_rho*0.5d0)  then
c
c           Spin alpha:
c
               rhoa=rho(n,R_A)
               rho13a = rhoa**(1.d0/3.d0)
               rho43 = rho13a*rhoa
               gamma = rgamma(n,G_AA)
c              gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
               if (sqrt(gamma).gt.tol_rho)then
                  x = sqrt(gamma) / rho43
                  call kop(tol_rho,x,kalpha)
c                 dkadra = -(4d0/3d0)*x*dkadxa/rhoa
c                 dkadga = dkadxa*0.5d0/(rho43*dsqrt(gamma))
               else
                  x = 0d0
               endif
            else
               rhoa=0d0
               rho13a=0d0
               x = 0d0
            endif
            if(x.eq.0d0) then
               call kop(tol_rho,x,kalpha)
c              dkadra = 0d0
c              dkadga = 0d0
            endif
c     
c           Spin beta:
c
            if (rho(n,R_B).ge.tol_rho*0.5d0) then
c
               rhob=rho(n,R_B)
               rho13b = rhob**(1.d0/3.d0)
               rho43 = rho13b*rhob
               gamma = rgamma(n,G_BB)
c              gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
               if (sqrt(gamma).gt.tol_rho)then
                  x = sqrt(gamma) / rho43
                  call kop(tol_rho,x, kbeta)
c                 dkbdrb = -(4d0/3d0)*x*dkbdxb/rhob
c                 dkbdgb = dkbdxb*0.5d0/(rho43*dsqrt(gamma))
               else
                  x = 0d0
               endif
            else
               if(rho13a.eq.0.0d0) goto 20
               rhob=0d0
               rho13b=0d0
               x=0d0
            endif
            if(x.eq.0d0) then
               call kop(tol_rho,x, kbeta)
c              dkbdrb = 0d0
c              dkbdgb=  0d0
            endif

            banb = qabop*(rho13a*kalpha*rho13b*kbeta)/
     /           (rho13a*kalpha+rho13b*kbeta)

            if(banb.ne.0.0d0) then
c              dbabdra = banb*kbeta*rho13b/
c    /              (rho13a*kalpha+rho13b*kbeta)*
c    /              (1d0/(3d0*rhoa)+dkadra/kalpha)
c              dbabdrb = banb*kalpha*rho13a/
c    /              (rho13a*kalpha+rho13b*kbeta)*
c    /              (1d0/(3d0*rhob)+dkbdrb/kbeta)

c              dbabdga = banb*rho13b*kbeta/
c    /              ((rho13a*kalpha+rho13b*kbeta)*kalpha)*
c    *              dkadga
c              dbabdgb = banb*rho13a*kalpha/
c    /              ((rho13a*kalpha+rho13b*kbeta)*kbeta)*
c    *              dkbdgb

               hbab = hbabx(banb)
c              dhdab = dhdabx(banb)
            else
c              dbabdra =0d0
c              dbabdrb =0d0
c              dbabdga =0d0
c              dbabdgb =0d0
               hbab = 0d0
c              dhdab = 0d0
            endif

c           Ec = Ec - rhoa*rhob*hbab*qwght(n)*fac
            func(n) = func(n) - rhoa*rhob*hbab*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) - 
c    -           (rhob*hbab + rhoa*rhob*dhdab*dbabdra)*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB) - 
c    -           (rhoa*hbab + rhoa*rhob*dhdab*dbabdrb)*wght
c
c
c           if (x.gt.tol_rho) then
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - 
c    -              rhoa*rhob*dhdab*dbabdga*wght
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) - 
c    -              rhoa*rhob*dhdab*dbabdgb*wght
c           endif
               
c
c
 20      continue
c
      endif
c
      return
c
      contains
c
c     The combination of statement functions and derived types with
c     overloaded operators is not properly supported in the Fortran 
c     standard (apparently). Therefore the statement functions from
c     the original subroutine had to be transformed into contained 
c     functions.
c
c     WARNING: This code is EXTREMELY EVIL! Although there appears to be
c     only one function there are actually three with the same name,
c     each one returning results of a different data type. The trick is
c     that depending on the data type the the subroutine that contains
c     these functions changes its name and hence these different
c     functions of the same name do not lead to conflicts. The 
c     alternative would have been to add a forest of conditional 
c     compilation constructs (#ifdef's) to change the function names
c     in the declarations and all places where they are used. That 
c     would have been extremely ugly, so we are between a rock and a
c     hard place on this one.
c
      function hbabx(x) result(s)
#include "nwad.fh"
        implicit none
        type(nwad_dble), intent(in) :: x
        type(nwad_dble)             :: s
        s = (1.5214d0*x + 0.5764d0)/
     .      (x**2.0d0*(x**2.0d0+1.1284d0*x+0.3183d0))
      end function
c
      end
C>
C> \brief The Becke88 exchange GGA enhancement factor
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      subroutine nwxc_k_becke88_p(tol_rho,x, kalpha)
#else
      subroutine nwxc_k_becke88(tol_rho,x, kalpha)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_k_becke88_d2(tol_rho,x, kalpha)
#else
      subroutine nwxc_k_becke88_d3(tol_rho,x, kalpha)
#endif
#include "nwad.fh"
      implicit none
c
      double precision tol_rho
      type(nwad_dble)::x,kalpha
      double precision dkadxa
c
      double precision BETA, C
      Parameter (BETA = 0.0042D0)
      type(nwad_dble)::g,gdenom
      double precision dgdenom,dg
c     double precision arcsinh, darcsinh
c     arcsinh(x)=log(x+dsqrt(1d0+x*x))
c     darcsinh(x)=1d0/dsqrt(1d0+x*x)
c
c
c     Uniform electron gas constant
c
      C =  3d0*(0.75d0/acos(-1d0))**(1d0/3d0)

      if (x.gt.0d0)then
         gdenom = 1d0 + 6d0*BETA*x*asinh(x)
         g = 2d0*BETA*x*x / gdenom
c        dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
c        dg = g*(2d0/x-dgdenom/gdenom)
         
         kalpha= C + g
c        dkadxa = dg
         
      else
         kalpha= C
c        dkadxa = 0d0
      endif
      return
      end
C>
C> \brief The PBE96 exchange GGA enhancement factor
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      subroutine nwxc_k_pbe96_p(tol_rho,x, kalpha)
#else
      subroutine nwxc_k_pbe96(tol_rho,x, kalpha)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_k_pbe96_d2(tol_rho,x, kalpha)
#else
      subroutine nwxc_k_pbe96_d3(tol_rho,x, kalpha)
#endif
#include "nwad.fh"
      implicit none
c
      double precision tol_rho
      type(nwad_dble)::x,kalpha,deno
      double precision dkadxa
c
      double precision pi,um, uk, umk
      parameter(um=0.21951d0, uk=0.804d0, umk=um/uk)
      double precision C
      double precision forty8
c
c
c     Uniform electron gas constant
c
      pi = acos(-1.d0)
      C =  3d0*(0.75d0/pi)**(1d0/3d0)

      if (x.gt.0d0)then
         forty8=1d0/((48d0*pi*pi)**(2d0/3d0))
         deno=1d0/(1d0+um*x*x*forty8/uk)
         kalpha= C * (1d0 + uk - uk *deno)
c        dkadxa = C * (2d0*um*x*deno*deno*
c    *        forty8)
         
      else
         kalpha= C
c        dkadxa = 0d0
      endif
      return
      end
C>
C> \brief The Dirac exchange GGA enhancement factor
C>
C> Of course the Dirac functional is the exchange part for LDA
C> so this subroutine just returns a constant.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      subroutine nwxc_k_dirac_p(tol_rho,x, kalpha)
#else
      subroutine nwxc_k_dirac(tol_rho,x, kalpha)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_k_dirac_d2(tol_rho,x, kalpha)
#else
      subroutine nwxc_k_dirac_d3(tol_rho,x, kalpha)
#endif
#include "nwad.fh"
      implicit none
c
      double precision tol_rho
      type(nwad_dble)::x,kalpha
      double precision dkadxa
c
      double precision pi
      double precision C
c
c
c     Uniform electron gas constant
c
      pi = acos(-1.d0)
      C =  3d0*(0.75d0/pi)**(1d0/3d0)

      kalpha= C
c     dkadxa = 0d0

      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_op.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_op.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_op.F"
#endif
#undef NWAD_PRINT
C>
C> @}
