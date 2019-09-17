C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_op.F
C> The OP correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
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
      Subroutine nwpwxc_c_op(kop,param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                     func,Amat,Cmat)
c
C$Id: xc_op.F 20247 2011-04-28 18:58:49Z d3y133 $
c
      implicit none
c      
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      external kop             !< [Input] Subroutine to evaluate the
      !< GGA exchange enhancement factor (see Eq.(14) in [1]).
      double precision param(*)!< [Input] Parameters of functional
      !< - param(1): \f$ q^{\alpha\beta}_{OP} \f$ see Eq.(27) in [1] and
      !<   Table III in [2].
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)    !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c
      double precision QABOP,QAB88OP,QABPBOP
c     Parameter (QAB88OP=2.3670D0,QABPBOP=2.3789D0)
c
c References:
c    Tsuneda, Suzumura, Hirao, JCP 110, 10664 (1999)
c    Tsuneda, Suzumura, Hirao, JCP 111, 5656 (1999)
c
c***************************************************************************
c
      integer n
      double precision rho13, rho43, gamma, x
      double precision kalpha,kbeta, rho13a, rho13b,rhoa,rhob
      double precision banb, hbab, hbabx
      double precision dhdab,dhdabx,dkadra,dkbdrb,dkadxa,dkbdxb,
     A     dbabdra,dbabdrb,dbabdga,dbabdgb,dkadga,dkbdgb,
     A     dbabdka,dbabdkb
c
      hbabx(x) = (1.5214d0*x + 0.5764d0)/
     /           (x**2*(x**2+1.1284d0*x+0.3183d0))
      dhdabx(x) = -(4.5642d0*x**4+5.7391d0*x**3+
     +     2.4355*x**2+0.3669d0*x)/
     /           ((x**4+1.1284d0*x**3+0.3183d0*x**2)**2)
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
            rho43 = rho13a**4 
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = 0.25d0 * gamma
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
               call kop(tol_rho,x,kalpha, dkadxa)
               dkadra = -(4d0/3d0)*x*dkadxa/rhoa
               dkadga = (dkadxa/rho43)*0.5d0/dsqrt(gamma)
            else
               x=0d0
               call kop(tol_rho,x,kalpha, dkadxa)
               dkadra = 0d0
               dkadga = 0d0
            endif
c
c
          
            banb = qabop * rho13a * kalpha *0.5d0

            if(banb.ne.0) then
               dbabdra = banb*0.5d0*
     /              (1d0/(3d0*rhoa)+dkadra/kalpha)

               dbabdga = banb/kalpha*dkadga*0.5d0
               
               hbab = hbabx(banb)
               dhdab = dhdabx(banb)
            else
               dbabdra =0d0
               dbabdga =0d0
               hbab = 0d0
               dhdab = 0d0
            endif

c           Ec = Ec - rhoa**2*hbab*qwght(n)*fac
            func(n) = func(n) - rhoa**2*hbab*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) - 
     -           (rhoa*hbab + rhoa**2*dhdab*dbabdra)*wght

c
            if (x.gt.tol_rho) then
                Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - 
     -              rhoa**2*dhdab*dbabdga*wght
             endif
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
               if (dsqrt(gamma).gt.tol_rho)then
                  x = dsqrt(gamma) / rho43
                  call kop(tol_rho,x,kalpha, dkadxa)
                  
                  dkadra = -(4d0/3d0)*x*dkadxa/rhoa
                  dkadga = dkadxa*0.5d0/(rho43*dsqrt(gamma))
               else
                  x = 0d0
               endif
            else
               rhoa=0d0
               rho13a=0d0
               x = 0d0
            endif
            if(x.eq.0d0) then
               call kop(tol_rho,x, kalpha, dkadxa)
               dkadra = 0d0
               dkadga = 0d0
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
               if (dsqrt(gamma).gt.tol_rho)then
                  x = dsqrt(gamma) / rho43
                  call kop(tol_rho,x, kbeta, dkbdxb)
            
                  dkbdrb = -(4d0/3d0)*x*dkbdxb/rhob
                  dkbdgb = dkbdxb*0.5d0/(rho43*dsqrt(gamma))
               else
                  x = 0d0
               endif
            else
               if(rho13a.eq.0) goto 20
               rhob=0d0
               rho13b=0d0
               x=0d0
            endif
            if(x.eq.0d0) then
               call kop(tol_rho,x, kbeta, dkbdxb)
               dkbdrb = 0d0
               dkbdgb=  0d0
            endif

            banb = qabop*(rho13a*kalpha*rho13b*kbeta)/
     /           (rho13a*kalpha+rho13b*kbeta)

            if(banb.ne.0) then
               dbabdra = banb*kbeta*rho13b/
     /              (rho13a*kalpha+rho13b*kbeta)*
     /              (1d0/(3d0*rhoa)+dkadra/kalpha)
               dbabdrb = banb*kalpha*rho13a/
     /              (rho13a*kalpha+rho13b*kbeta)*
     /              (1d0/(3d0*rhob)+dkbdrb/kbeta)

               dbabdga = banb*rho13b*kbeta/
     /              ((rho13a*kalpha+rho13b*kbeta)*kalpha)*
     *              dkadga
               dbabdgb = banb*rho13a*kalpha/
     /              ((rho13a*kalpha+rho13b*kbeta)*kbeta)*
     *              dkbdgb

               hbab = hbabx(banb)
               dhdab = dhdabx(banb)
            else
               dbabdra =0d0
               dbabdrb =0d0
               dbabdga =0d0
               dbabdgb =0d0
               hbab = 0d0
               dhdab = 0d0
            endif

c           Ec = Ec - rhoa*rhob*hbab*qwght(n)*fac
            func(n) = func(n) - rhoa*rhob*hbab*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) - 
     -           (rhob*hbab + rhoa*rhob*dhdab*dbabdra)*wght
            Amat(n,D1_RB) = Amat(n,D1_RB) - 
     -           (rhoa*hbab + rhoa*rhob*dhdab*dbabdrb)*wght
c
c
            if (x.gt.tol_rho) then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - 
     -              rhoa*rhob*dhdab*dbabdga*wght
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) - 
     -              rhoa*rhob*dhdab*dbabdgb*wght
            endif
               
c
c
 20      continue
c
      endif
c
      return
      end
C>
C> \brief The Becke88 exchange GGA enhancement factor
C>
      subroutine nwpwxc_k_becke88(tol_rho,x, kalpha, dkadxa)
      implicit none
c
      double precision tol_rho,x
      double precision kalpha,dkadxa
c
      double precision BETA, C
      Parameter (BETA = 0.0042D0)
      double precision g,gdenom,dgdenom,dg
      double precision arcsinh, darcsinh
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c
c
c     Uniform electron gas constant
c
      C =  3d0*(0.75d0/acos(-1d0))**(1d0/3d0)

      if (x.gt.0d0)then
         gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
         g = 2d0*BETA*x*x / gdenom
         dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
         dg = g*(2d0/x-dgdenom/gdenom)
         
         kalpha= C + g
         dkadxa = dg
         
      else
         kalpha= C
         dkadxa = 0d0
      endif
      return
      end
C>
C> \brief The PBE96 exchange GGA enhancement factor
C>
      subroutine nwpwxc_k_pbe96(tol_rho,x, kalpha, dkadxa)
      implicit none
c
      double precision tol_rho,x
      double precision kalpha,dkadxa
c
      double precision pi,um, uk, umk
      parameter(um=0.21951d0, uk=0.804d0, umk=um/uk)
      double precision C
      double precision forty8,deno
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
         dkadxa = C * (2d0*um*x*deno*deno*
     *        forty8)
         
      else
         kalpha= C
         dkadxa = 0d0
      endif
      return
      end
C>
C> \brief The Dirac exchange GGA enhancement factor
C>
C> Of course the Dirac functional is the exchange part for LDA
C> so this subroutine just returns a constant.
C>
      subroutine nwpwxc_k_dirac(tol_rho,x, kalpha, dkadxa)
      implicit none
c
      double precision tol_rho,x
      double precision kalpha,dkadxa
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
      dkadxa = 0d0

      return
      end
C>
C> @}
