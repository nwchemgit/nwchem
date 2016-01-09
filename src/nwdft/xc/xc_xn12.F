c $Id: xc_eval_fnl.F 27518 2015-09-17 09:34:57Z alogsdail $
corg      Subroutine N12X(NGrid,RA,RB,D1RA,D1RB,F,D1F,FG,ijzy)
      subroutine xc_xn12(tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     amat, cmat, nq, ipol, ex, 
     &                     qwght, ldew, func, ijzy)
      Implicit Real*8(A-H,O-Z)
************************************************************************
*                                                                      *
*  N12 evaluates the exchange part of the N12 and N12-SX               *
*        functionals on the grid.                                      *
*                                                                      *
*    OUTPUT:                                                           *
*      F      - Functional values                                      *
*      D1F    - First derivatives with respect to RA, RB, GA, GB       *
*                                                                      *
*    INPUT:                                                            *
*       ijzy - 1 N12                                                   *
*       ijzy - 2 N12-SX                                                *
*                                                                      *
*      RA,B     - Spin densities                                       *
*      D1RA,B   - Spin density gradients                               *
*      NGrid    - number of grid points                                *
*                                                                      *
*  Analytic second derivatives are not available yet!!!                *
*                                                                      *
*  RP (09/12)                                                          *
*                                                                      *
************************************************************************
#include "errquit.fh"
#include "dft2drv.fh"
c      
      integer ijzy
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew,   uselc
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
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision tol_rho, pi
c      Real*8 RA(*),RB(*),D1RA(NGrid,3),D1RB(NGrid,3),F(*),D1F(NGrid,*)
      Save One,Two,Three,Four,Five,Six,Seven,Eight,Nine
      Data One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Five/5.0d0/,
     $  Six/6.0d0/,Seven/7.0d0/,Eight/8.0d0/,Nine/9.0d0/
C
C
      TinyDFT = tol_rho
      G    = 0.004d+0
      ome  = 2.5d+0        
      If (ijzy.eq.1) then
c
c      N12
c      
       CC00 =  1.00000D+00
       CC01 =  5.07880D-01
       CC02 =  1.68233D-01
       CC03 =  1.28887D-01
       CC10 =  8.60211D-02
       CC11 = -1.71008D+01
       CC12 =  6.50814D+01
       CC13 = -7.01726D+01
       CC20 = -3.90755D-01
       CC21 =  5.13392D+01
       CC22 = -1.66220D+02
       CC23 =  1.42738D+02
       CC30 =  4.03611D-01
       CC31 = -3.44631D+01
       CC32 =  7.61661D+01
       CC33 = -2.41834D+00
      Else If (ijzy.eq.2) then
c
c      N12-SX 
c      
       CC00 =  6.81116D-01
       CC01 =  1.88858D+00
       CC02 =  1.78590D+00
       CC03 =  8.79456D-01
       CC10 = -8.12270D-02
       CC11 = -1.08723D+00
       CC12 = -4.18682D+00
       CC13 = -3.00000D+01
       CC20 =  5.36236D-01
       CC21 = -5.45678D+00
       CC22 =  3.00000D+01
       CC23 =  5.51105D+01
       CC30 = -7.09913D-01
       CC31 =  1.30001D+01
       CC32 = -7.24877D+01
       CC33 =  2.98363D+01
      End If
      F12   = Two * Six
      F24   = Four * Six
      F28   = Four * Seven
      F2o3  = Two / Three
      F3o2  = Three / Two
      F1o3  = One / Three
      F4o3  = Four / Three
      F7o3  = Seven / Three
      F8o3  = Eight / Three
      F10o3 = F2o3 * Five
      F28o9 = F28 / Nine
      Pi    = ACos(-One)
C
      Ax = -F3o2*(F4o3*PI)**(-F1o3)
C
C     Alpha contributions ...
C
      Do 10 n = 1, nq
        If(rho(n,2).gt.TinyDFT) then
          pX    = rho(n,2)
cold          GamX2 = D1RA(iG,1)**2 + D1RA(iG,2)**2 + D1RA(iG,3)**2
          GamX2 =(delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E     = Ax*pX**F4o3
          FV    = ome*pX**F1o3/(One+ome*pX**F1o3)
c
       FN12 =CC00  +CC01      *U+CC02      *U**2+CC03      *U**3+
     $    CC10*FV   +CC11*FV   *U+CC12*FV   *U**2+CC13*FV   *U**3+
     $    CC20*FV**2+CC21*FV**2*U+CC22*FV**2*U**2+CC23*FV**2*U**3+
     $    CC30*FV**3+CC31*FV**3*U+CC32*FV**3*U**2+CC33*FV**3*U**3
c
       ex = ex + E*FN12*qwght(n)
       if(ldew) func(n)=func(n) + E*FN12
c
c First Derivative
c
          ER = F4o3*E/pX
          S    = Sqrt(S2)
          GamX = Sqrt(GamX2)
          SR   = -F4o3*S/pX
          SG   = S/GamX
          US   = Two*G*S/((One+G*S*S)**2)
c        
          dFVdR =  (ome/(Three*pX**F2o3))
     $                *(One+ome*pX**F1o3)**(-Two)
c
        dFN12dU=CC01+2.0d0*CC02      *U+3.0d0*CC03      *U**2+
     $     CC11*FV   +2.0d0*CC12*FV   *U+3.0d0*CC13*FV   *U**2+
     $     CC21*FV**2+2.0d0*CC22*FV**2*U+3.0d0*CC23*FV**2*U**2+
     $     CC31*FV**3+2.0d0*CC32*FV**3*U+3.0d0*CC33*FV**3*U**2
        dFN12dV=
     $     CC10+CC11   *U+CC12   *U**2+CC13   *U**3+
     $     2.0d0*CC20*FV   +2.0d0*CC21*FV*U+
     $     2.0d0*CC22*FV*U**2+2.0d0*CC23*FV*U**3+
     $     3.0d0*CC30*FV**2+3.0d0*CC31*FV**2*U+
     $     3.0d0*CC32*FV**2*U**2+3.0d0*CC33*FV**2*U**3
c
        dFN12dR  = dFN12dU*US*SR+dFN12dV*dFVdR
        dFN12dG  = dFN12dU*US*SG
c
            amat(n,D1_RA) = amat(n,D1_RA) + ER*FN12
     $                      + E*dFN12dR
            cmat(n,D1_GAA) = cmat(n,D1_GAA) + E*dFN12dG/(Two*GamX)
          endIf
   10   Continue
C
C     Beta contributions ...
C
      Do 20 n = 1, nq
        If(rho(n,3).gt.TinyDFT) then
          pX    = rho(n,3)
          GamX2 =(delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2))
cold          GamX2 = D1RB(iG,1)**2 + D1RB(iG,2)**2 + D1RB(iG,3)**2
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E     = Ax*pX**F4o3
          FV    = ome*pX**F1o3/(One+ome*pX**F1o3)
c
       FN12 =CC00  +CC01      *U+CC02      *U**2+CC03      *U**3+
     $    CC10*FV   +CC11*FV   *U+CC12*FV   *U**2+CC13*FV   *U**3+
     $    CC20*FV**2+CC21*FV**2*U+CC22*FV**2*U**2+CC23*FV**2*U**3+
     $    CC30*FV**3+CC31*FV**3*U+CC32*FV**3*U**2+CC33*FV**3*U**3
c
       ex = ex + E*FN12*qwght(n)
       if(ldew) func(n)=func(n) + E*FN12
c
c First Derivative
c
          ER = F4o3*E/pX
          S    = Sqrt(S2)
          GamX = Sqrt(GamX2)
          SR   = -F4o3*S/pX
          SG   = S/GamX
          US   = Two*G*S/((One+G*S*S)**2)
c        
          dFVdR =  (ome/(Three*pX**F2o3))
     $                *(One+ome*pX**F1o3)**(-Two)
c
        dFN12dU=CC01+2.0d0*CC02      *U+3.0d0*CC03      *U**2+
     $     CC11*FV   +2.0d0*CC12*FV   *U+3.0d0*CC13*FV   *U**2+
     $     CC21*FV**2+2.0d0*CC22*FV**2*U+3.0d0*CC23*FV**2*U**2+
     $     CC31*FV**3+2.0d0*CC32*FV**3*U+3.0d0*CC33*FV**3*U**2
c
        dFN12dV=
     $     CC10+CC11   *U+CC12   *U**2+CC13   *U**3+
     $     2.0d0*CC20*FV   +2.0d0*CC21*FV*U+
     $     2.0d0*CC22*FV*U**2+2.0d0*CC23*FV*U**3+
     $     3.0d0*CC30*FV**2+3.0d0*CC31*FV**2*U+
     $     3.0d0*CC32*FV**2*U**2+3.0d0*CC33*FV**2*U**3
c
        dFN12dR  = dFN12dU*US*SR+dFN12dV*dFVdR
        dFN12dG  = dFN12dU*US*SG
c
            amat(n,D1_RB) = amat(n,D1_RB) + ER*FN12
     $                      + E*dFN12dR
            cmat(n,D1_GBB) = cmat(n,D1_GBB) + E*dFN12dG/(Two*GamX)
          endIf
   20   Continue
      Return
      End
      Subroutine xc_xn12_d2()
      call errquit(' not coded ',0,0)
      return
      end
