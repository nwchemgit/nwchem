      Subroutine xc_xpbe96(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,
     &                     ldew,func)     
c
c    $$
c
      implicit none
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
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
c     Sampling Matrices for the XC Potential & Energy
c
      double precision amat(nq,ipol), cmat(nq,3,ipol)

c
c     Compute the partial derivatives of the PBE exchange functional 
c
      double precision TOLL, tol_rho, p2, p2t,
     &     ONE, TWO, SIX, ONEHALF, FOURTHIRDS
      double precision pi,tp3,um,uk,ul
      double precision fx,fxp,x,fach

      Parameter (TOLL=1.E-20)
c
c     p2 = -(3/4)*(6/pi)**(1/3)
c
      parameter (p2 = -0.9305257363490993d+00)

      parameter (one=1.d0, two=2.d0, six=6.d0)
      parameter (onehalf=0.5d0, fourthirds=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c
c***************************************************************************

      integer n
      double precision cd13a, cd43a, gna,  sa,  ca,
     &                 pa, gnb, sb, cb,  pb, cd13b, 
     &                 cd43b, gaa, g, gbb, gamma,
     &     cd13,cd43, exlda, vxlda
c
      fx(x)=uk-uk/(1.d0+ul*x*x)
      fxp(x)=2.d0*x*um/(1.d0+ul*x*x)**2
      tp3=0.5/((3.d0*pi*pi)**(1.d0/3.d0))
      p2t=-0.75d0*(3.d0/pi)**(1.d0/3.d0)
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            cd13 = rho(n,1)**(1.d0/3.d0)
            cd43 = cd13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gna = dsqrt( gamma )
               sa=gna*tp3/cd43
            else
               gna = 0.d0
               sa=0.d0
            endif
            ca=fx(sa)
            pa=fxp(sa)
c
            exlda=cd43*p2t
            vxlda=cd13*p2t*fourthirds
            if (lfac)then
               amat(n,1) = amat(n,1) + vxlda*fac
            endif
c
            if (nlfac)then
               amat(n,1) = amat(n,1) + vxlda*fac*(ca - sa*pa)
            endif
c
            if (lfac)then
               ex = ex + exlda*qwght(n)*fac
               if(ldew)func(n)=func(n)+exlda*fac
            endif
c
            if (nlfac)then
               ex = ex + exlda*ca*qwght(n)*fac
               if(ldew)func(n)=func(n)+exlda*ca*fac
            endif
c
            if (gna.gt.tol_rho)  then
               g = tp3*p2t*pa/gna*fac
               cmat(n,1,1) = cmat(n,1,1) + g*delrho(n,1,1)
               cmat(n,2,1) = cmat(n,2,1) + g*delrho(n,2,1)
               cmat(n,3,1) = cmat(n,3,1) + g*delrho(n,3,1)
            endif
c
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
        fach=fac*0.5d0
         do 20 n = 1, nq
C
C          alpha bit
C
            if (rho(n,2).gt.tol_rho) then
              cd13a = (two*rho(n,2))**(1.d0/3.d0)
              cd43a = cd13a**4 
            else
              cd13a=0.d0
              cd43a=0.d0
            endif
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            if (gamma.gt.tol_rho)then
              gna = two*dsqrt( gamma )
              sa=gna*tp3/cd43a
            else
              gna = 0.d0
              sa=0.d0
            endif
            ca=fx(sa)
            pa=fxp(sa)
c
            exlda=cd43a*p2t
            vxlda=cd13a*p2t*fourthirds
            if (lfac)then
              amat(n,1) = amat(n,1) + vxlda*fac
            endif
c           
            if (nlfac)then
              amat(n,1) = amat(n,1) + vxlda*fac*(ca - sa*pa)
            endif
c           
            if (lfac)then
              ex = ex + exlda*qwght(n)*fach
              if(ldew)func(n)=func(n)+exlda*fach
            endif
c           
            if (nlfac)then
              ex = ex + exlda*ca*qwght(n)*fach
              if(ldew)func(n)=func(n)+exlda*ca*fach
            endif
c           
            if (gna.gt.tol_rho)  then
              g = 2.d0*tp3*p2t*pa/gna*fac
              cmat(n,1,1) = cmat(n,1,1) + g*delrho(n,1,1)
              cmat(n,2,1) = cmat(n,2,1) + g*delrho(n,2,1)
              cmat(n,3,1) = cmat(n,3,1) + g*delrho(n,3,1)
            endif
C
C          beta bit
C
            if (rho(n,3).gt.tol_rho) then
              cd13b = (two*rho(n,3))**(1.d0/3.d0)
              cd43b = cd13b**4 
            else
              cd13b=0.d0
              cd43b=0.d0
            endif
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2)
            if (gamma.gt.tol_rho)then
              gnb = two*dsqrt( gamma )
              sb=gnb*tp3/cd43b
            else
              gnb = 0.d0
              sb=0.d0
            endif
            cb=fx(sb)
            pb=fxp(sb)
c
            exlda=cd43b*p2t
            vxlda=cd13b*p2t*fourthirds
            if (lfac)then
              amat(n,2) = amat(n,2) + vxlda*fac
            endif
c           
            if (nlfac)then
              amat(n,2) = amat(n,2) + vxlda*fac*(cb- sb*pb)
            endif
c           
            if (lfac)then
              ex = ex + exlda*qwght(n)*fach
              if(ldew)func(n)=func(n)+exlda*fach
            endif
c           
            if (nlfac)then
              ex = ex + exlda*cb*qwght(n)*fach
              if(ldew)func(n)=func(n)+exlda*cb*fach
            endif
c           
            if (gnb.gt.0d0)  then
              g = 2.d0*tp3*p2t*pb/gnb*fac
              cmat(n,1,2) = cmat(n,1,2) + g*delrho(n,1,2)
              cmat(n,2,2) = cmat(n,2,2) + g*delrho(n,2,2)
              cmat(n,3,2) = cmat(n,3,2) + g*delrho(n,3,2)
            endif
c           
   20     continue
c
      endif
c
      return
      end
