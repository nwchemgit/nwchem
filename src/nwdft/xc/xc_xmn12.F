      subroutine xc_xmn12l(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat, nq, ipol, ex, 
     &     qwght, ldew, func)
      implicit none
#include "errquit.fh"
#include "dft2drv.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew,   uselc
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision delrho(nq,3,ipol)
      double precision tau(nq,ipol)
      double precision qwght(nq)
      double precision Amat(nq,ipol), Cmat(nq,*),mmat(nq,*)
      double precision tol_rho
      external  xc_xmn12l_os
      call xc_os2cs_mgga(xc_xmn12l_os,
     &     tol_rho, fac,lfac,nlfac, rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq, ipol, ex, 
     &     qwght, ldew, func)
      return
      end
      subroutine xc_xmn12l_os(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq,  ex, 
     &                     qwght, ldew, func, fact_cs)
      implicit none
#include "errquit.fh"
#include "dft2drv.fh"
c
c wrapper for n12 open-shell
      double precision fac, Ex
      double precision fact_cs
      integer nq
      logical lfac, nlfac,ldew,   uselc
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq)
      double precision delrho(nq,3), tau(nq)
      double precision qwght(nq)
      double precision Amat(nq), Cmat(nq), mmat(nq)
      double precision tol_rho

      call xc_xmn12(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     N     nq,  ex,   qwght, ldew, func,
     I      1    ,fact_cs)
      return
      end
      subroutine xc_xmn12sx(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat, nq, ipol, ex, 
     &     qwght, ldew, func)
      implicit none
#include "errquit.fh"
#include "dft2drv.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew,   uselc
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision delrho(nq,3,ipol)
      double precision tau(nq,ipol)
      double precision qwght(nq)
      double precision Amat(nq,ipol), Cmat(nq,*),mmat(nq,*)
      double precision tol_rho
      external  xc_xmn12sx_os
      call xc_os2cs_mgga(xc_xmn12sx_os,
     &     tol_rho, fac,lfac,nlfac, rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq, ipol, ex, 
     &     qwght, ldew, func)
      return
      end
      subroutine xc_xmn12sx_os(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq,  ex, 
     &                     qwght, ldew, func, fact_cs)
      implicit none
#include "errquit.fh"
#include "dft2drv.fh"
c
c wrapper for n12 open-shell
      double precision fac, Ex
      double precision fact_cs
      integer nq
      logical lfac, nlfac,ldew,   uselc
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq)
      double precision delrho(nq,3), tau(nq)
      double precision qwght(nq)
      double precision Amat(nq), Cmat(nq), mmat(nq)
      double precision tol_rho

      call xc_xmn12(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     N     nq,  ex,   qwght, ldew, func,
     I      2    ,fact_cs)
      return
      end
      Subroutine xc_xmn12(tol_rho, fac,lfac,nlfac,
     R     rho, delrho, tau,
     &     Amat, Cmat, mmat,
     n     nq,  Ex, 
     &     qwght, ldew, func,
     I     ijzy, fact_cs)
cSubroutine MN12x(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
**************************************************************
*                                                            *
*  xc_xmn12 evaluates the exchange part of MN12-L and MN12-SX   *
*        functionals on the grid.                            *
*                                                            *
*  OUTPUT:                                                   *
*     F   - Functional values                                *
*     D1F - First derivatives with respect to RA, RB, GA, GB *
*              TA, TB                                        *
*                                                            *
*  INPUT:                                                    *
*       ijzy - 1 MN12-L                                      *
*       ijzy - 2 MN12-SX                                     *
*                                                            *
*     RA,B   - Spin densities                                *
*     D1RA,B - Spin density gradients                        *
*     TA,B   - Spin kinetic energy densities                 *
*     NGrid  - number of grids                               *
*                                                            *
*  RP (09/12)                                                *
*                                                            *
**************************************************************
      
c
      implicit none
c
#include "errquit.fh"
c      
      double precision fac, Ex
      double precision fact_cs
      integer nq
      logical lfac, nlfac,ldew,   uselc
      double precision func(*)  ! value of the functional [output]
      double precision qwght(nq)
#if 0      
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision delrho(nq,3,ipol)
      double precision tau(nq,ipol)
      double precision Amat(nq,ipol), Cmat(nq,*), Mmat(nq,*)
#else      
      double precision rho(nq)
      double precision delrho(nq,3)
      double precision tau(nq)
      double precision Amat(nq), Cmat(nq), Mmat(nq)
#endif      
      double precision tol_rho, pi
c
c     kinetic energy density   or  tau
c
      double precision tau_in
c
c      functional derivatives
c
      double precision dWdT, dTdR, dTdTau
c
c     Intermediate derivative results, etc.
c
      integer n, ijzy
c
      double precision Ax, s, s2
c
      double precision F1o3,F2o3,F3o5,F4o3,F5o3,F48,F81
      double precision f10o3,f12,f24,f28,f28o9,f3o2,f6,
     f     f7o3,f8o3
      double precision Fx1, Fx2
c
      double precision cc000,cc001,cc002,cc003,cc004,cc005
      double precision cc010,cc011,cc012,cc013,cc014
      double precision cc020,cc021,cc022,cc023
      double precision cc030,cc031,cc032
      double precision cc100,cc101,cc102,cc103,cc104
      double precision cc110,cc111,cc112,cc113
      double precision cc120,cc121,cc122
      double precision cc200,cc201,cc202,cc203
      double precision cc210,cc211,cc212
      double precision cc300,cc301,cc302
      double precision rho43, rho13, rhoo, rho53
      double precision gamx2, gamx
      double precision TauUEG, Tsig, Wsig
      double precision px,u,g,ome
      double precision sg,s3,us
      double precision dFMN12dFT, dFMN12dFU, dFMN12dFV, dFMN12dG,
     d     dFMN12dR, dFMN12dT
      double precision dFTdR, dFTdTau, dFVdR
      double precision e,er,sr,ft,fmn12,fu,fv
#if 0
      INTEGER NGrid
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC      
      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
#endif      
      double precision One, Two, Three, Four, Five, Six, Seven, Eight,
     c     nine
      Save One,Two,Three,Four,Five,Six,Seven,Eight,Nine
      Data One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Five/5.0d0/,
     $     Six/6.0d0/,Seven/7.0d0/,Eight/8.0d0/,Nine/9.0d0/
c
c      if(ipol.eq.1) call errquit(' not reaaaady ',0,0)
c
        G   = 0.004d+0
        ome = 2.5d+0
        CC000 = 0d0
        CC001 = 0d0
        CC002 = 0d0
        CC003 = 0d0
        CC004 = 0d0
        CC005 = 0d0
        CC010 = 0d0
        CC011 = 0d0
        CC012 = 0d0
        CC013 = 0d0
        CC014 = 0d0
        CC020 = 0d0
        CC021 = 0d0
        CC022 = 0d0
        CC023 = 0d0
        CC030 = 0d0
        CC031 = 0d0
        CC032 = 0d0
        CC100 = 0d0
        CC101 = 0d0
        CC102 = 0d0
        CC103 = 0d0
        CC104 = 0d0
        CC110 = 0d0
        CC111 = 0d0
        CC112 = 0d0
        CC113 = 0d0
        CC120 = 0d0
        CC121 = 0d0
        CC122 = 0d0
        CC200 = 0d0
        CC201 = 0d0
        CC202 = 0d0
        CC203 = 0d0
        CC210 = 0d0
        CC211 = 0d0
        CC212 = 0d0
        CC300 = 0d0
        CC301 = 0d0
        CC302 = 0d0
      If (ijzy.eq.1) then
c
c       MN12-L
c      
        CC000 =  6.735981D-01
        CC001 = -2.270598D+00
        CC002 = -2.613712D+00
        CC003 =  3.993609D+00
        CC004 =  4.635575D+00
        CC005 =  1.250676D+00
        CC010 =  8.444920D-01
        CC011 = -1.301173D+01
        CC012 = -1.777730D+01
        CC013 = -4.627211D+00
        CC014 =  5.976605D+00
        CC020 =  1.142897D+00
        CC021 = -2.040226D+01
        CC022 = -2.382843D+01
        CC023 =  7.119109D+00
        CC030 = -2.335726D+01
        CC031 = -1.622633D+01
        CC032 =  1.482732D+01
        CC100 =  1.449285D+00
        CC101 =  1.020598D+01
        CC102 =  4.407450D+00
        CC103 = -2.008193D+01
        CC104 = -1.253561D+01
        CC110 = -5.435031D+00
        CC111 =  1.656736D+01
        CC112 =  2.000229D+01
        CC113 = -2.513105D+00
        CC120 =  9.658436D+00
        CC121 = -3.825281D+00
        CC122 = -2.500000D+01
        CC200 = -2.070080D+00
        CC201 = -9.951913D+00
        CC202 =  8.731211D-01
        CC203 =  2.210891D+01
        CC210 =  8.822633D+00
        CC211 =  2.499949D+01
        CC212 =  2.500000D+01
        CC300 =  6.851693D-01
        CC301 = -7.406948D-02
        CC302 = -6.788000D-01
      Else If (ijzy.eq.2) then
c
c       MN12-SX
c      
        CC000 =  5.226556D-01
        CC001 = -2.681208D-01
        CC002 = -4.670705D+00
        CC003 =  3.067320D+00
        CC004 =  4.095370D+00
        CC005 =  2.653023D+00
        CC010 =  5.165969D-01
        CC011 = -2.035442D+01
        CC012 = -9.946472D+00
        CC013 =  2.938637D+00
        CC014 =  1.131100D+01
        CC020 =  4.752452D+00
        CC021 = -3.061331D+00
        CC022 = -2.523173D+01
        CC023 =  1.710903D+01
        CC030 = -2.357480D+01
        CC031 = -2.727754D+01
        CC032 =  1.603291D+01
        CC100 =  1.842503D+00
        CC101 =  1.927120D+00
        CC102 =  1.107987D+01
        CC103 = -1.182087D+01
        CC104 = -1.117768D+01
        CC110 = -5.821000D+00
        CC111 =  2.266545D+01
        CC112 =  8.246708D+00
        CC113 = -4.778364D+00
        CC120 =  5.329122D-01
        CC121 = -6.666755D+00
        CC122 =  1.671429D+00
        CC200 = -3.311409D+00
        CC201 =  3.415913D-01
        CC202 = -6.413076D+00
        CC203 =  1.038584D+01
        CC210 =  9.026277D+00
        CC211 =  1.929689D+01
        CC212 =  2.669232D+01
        CC300 =  1.517278D+00
        CC301 = -3.442503D+00
        CC302 =  1.100161D+00
      else
         call errquit(' xc_xmn12: ijzy not valid ',ijzy,0)
      End If
c      Tol_Rho = 1.0d-8
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
      PI    = Four * ATan(One)
C
C     Local Spin Density factor.
C
      Ax = -F3o2*(F4o3*PI)**(-F1o3)
      Do 10 n = 1, nq
c
c alpha component
c
        If(rho(n).gt.tol_rho.and.tau(n).gt.tol_rho) then
          F3o5 = Three/Five
          F6   = Six
          F5o3 = Five/Three
          pX    = rho(n)
          Gamx2 = delrho(n,1)*delrho(n,1) +
     y         delrho(n,2)*delrho(n,2)  +
     z         delrho(n,3)*delrho(n,3)  
c          GamX2 = D1RA(i)
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E = Ax*pX**F4o3
c          tau_in    = tau(n)*2D0
          tau_in    = tau(n)
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3) 
          Tsig =TauUEG/tau_in
          Wsig =(Tsig - One)/(Tsig + One)
c          
          FU = U
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FT = Wsig          
c
          FMN12 = CC000 + FT*CC001 + FT**2*CC002 + FT**3*CC003 + 
     $     FT**4*CC004 + FT**5*CC005 + 
     $    FU*CC010 + FT*FU*CC011 + FT**2*FU*CC012 + FT**3*FU*CC013 + 
     $     FT**4*FU*CC014 + FU**2*CC020 + 
     $    FT*FU**2*CC021 + FT**2*FU**2*CC022 + FT**3*FU**2*CC023 + 
     $     FU**3*CC030 + FT*FU**3*CC031 + 
     $    FT**2*FU**3*CC032 + FV*CC100 + FV*FT*CC101 + 
     $     FV*FT**2*CC102 + FV*FT**3*CC103 + 
     $    FV*FT**4*CC104 + FV*FU*CC110 + FV*FT*FU*CC111 + 
     $     FV*FT**2*FU*CC112 + FV*FT**3*FU*CC113 + 
     $    FV*FU**2*CC120 + FV*FT*FU**2*CC121 + FV*FT**2*FU**2*CC122 + 
     $     FV**2*CC200 + FV**2*FT*CC201 + 
     $    FV**2*FT**2*CC202 + FV**2*FT**3*CC203 + FV**2*FU*CC210 + 
     $     FV**2*FT*FU*CC211 + 
     $    FV**2*FT**2*FU*CC212 + FV**3*CC300 + FV**3*FT*CC301 + 
     $     FV**3*FT**2*CC302
c
c          F(i) = F(i) + E*FMN12
          Ex = Ex + E*FMN12*qwght(n)*fact_cs
          if(ldew) func(n)=func(n) + E*FMN12*fact_cs

c
c          If(MaxDrF.ge.1) then
          If(.true.) then
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
            dWdT = Two/((One + Tsig)**2)
            dTdR = ((F6*PI*PI)**F2o3)*(pX**F2o3)/tau_in
            dTdTau = -TauUEG/tau_in**2
c
            dFTdR = dWdT*dTdR
            dFTdTau=dWdT*dTdTau
c     
          dFMN12dFV = CC100 + FT*CC101 + FT**2*CC102 + FT**3*CC103 + 
     $     FT**4*CC104 + FU*CC110 + 
     $    FT*FU*CC111 + FT**2*FU*CC112 + FT**3*FU*CC113 + 
     $     FU**2*CC120 + FT*FU**2*CC121 + 
     $    FT**2*FU**2*CC122 + Two*FV*CC200 + Two*FV*FT*CC201 + 
     $     Two*FV*FT**2*CC202 + Two*FV*FT**3*CC203 + 
     $    Two*FV*FU*CC210 + Two*FV*FT*FU*CC211+Two*FV*FT**2*FU*CC212+ 
     $     Three*FV**2*CC300 + Three*FV**2*FT*CC301 + 
     $    Three*FV**2*FT**2*CC302
c
          dFMN12dFU = CC010 + FT*CC011 + FT**2*CC012 + FT**3*CC013 + 
     $     FT**4*CC014 + Two*FU*CC020 + 
     $    Two*FT*FU*CC021 + Two*FT**2*FU*CC022 + Two*FT**3*FU*CC023 + 
     $     Three*FU**2*CC030 + Three*FT*FU**2*CC031 + 
     $    Three*FT**2*FU**2*CC032 + FV*CC110 + FV*FT*CC111 + 
     $     FV*FT**2*CC112 + FV*FT**3*CC113 + 
     $    Two*FV*FU*CC120 + Two*FV*FT*FU*CC121+Two*FV*FT**2*FU*CC122+ 
     $     FV**2*CC210 + FV**2*FT*CC211 + 
     $    FV**2*FT**2*CC212
c    
          dFMN12dFT = CC001 + Two*FT*CC002 + Three*FT**2*CC003 + 
     $     Four*FT**3*CC004 + Five*FT**4*CC005 + FU*CC011 + 
     $    Two*FT*FU*CC012 + Three*FT**2*FU*CC013+Four*FT**3*FU*CC014+ 
     $     FU**2*CC021 + Two*FT*FU**2*CC022 + 
     $    Three*FT**2*FU**2*CC023 + FU**3*CC031 + Two*FT*FU**3*CC032 + 
     $     FV*CC101 + Two*FV*FT*CC102 + 
     $    Three*FV*FT**2*CC103 + Four*FV*FT**3*CC104 + FV*FU*CC111 + 
     $     Two*FV*FT*FU*CC112 + Three*FV*FT**2*FU*CC113 + 
     $    FV*FU**2*CC121 + Two*FV*FT*FU**2*CC122 + FV**2*CC201 + 
     $     Two*FV**2*FT*CC202 + 
     $    Three*FV**2*FT**2*CC203 + FV**2*FU*CC211 + 
     $     Two*FV**2*FT*FU*CC212 + FV**3*CC301 + Two*FV**3*FT*CC302
c
        dFMN12dR = dFMN12dFV*dFVdR+dFMN12dFU*US*SR+dFMN12dFT*dFTdR
        dFMN12dG = dFMN12dFU*US*SG
        dFMN12dT = dFMN12dFT*dFTdTau
c
#if 0        
            Amat(n,1) = Amat(n,1) + ER*FMN12
     $                      + E*dFMN12dR
            Cmat(n,1)=  Cmat(n,1) + E*dFMN12dG/(Two*GamX)
            Mmat(n,1)=  Mmat(n,1) + E*dFMN12dT
#else            
            Amat(n) = Amat(n) + ER*FMN12
     $                      + E*dFMN12dR
            Cmat(n)=  Cmat(n) + E*dFMN12dG/(Two*GamX)
            Mmat(n)=  Mmat(n) + E*dFMN12dT
#endif            
            endIf
         endIf
#if 0         
c
c beta component
c
        If(rho(n,3).gt.tol_rho.and.tau(n,2).gt.tol_rho) then
          F3o5 = Three/Five
          F6   = Six
          F5o3 = Five/Three
c          pX    = RB(i)
c          GamX2 = D1RB(i) 
          pX    = rho(n,3)
          Gamx2 = delrho(n,1,2)*delrho(n,1,2) +
     y         delrho(n,2,2)*delrho(n,2,2)  +
     z         delrho(n,3,2)*delrho(n,3,2)  
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E = Ax*pX**F4o3
          tau_in    = tau(n,2)*2d0
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3) 
          Tsig =TauUEG/tau_in
          Wsig =(Tsig - One)/(Tsig + One)
c          
          FU = U
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FT = Wsig          
c
          FMN12 = CC000 + FT*CC001 + FT**2*CC002 + FT**3*CC003 + 
     $     FT**4*CC004 + FT**5*CC005 + 
     $    FU*CC010 + FT*FU*CC011 + FT**2*FU*CC012 + FT**3*FU*CC013 + 
     $     FT**4*FU*CC014 + FU**2*CC020 + 
     $    FT*FU**2*CC021 + FT**2*FU**2*CC022 + FT**3*FU**2*CC023 + 
     $     FU**3*CC030 + FT*FU**3*CC031 + 
     $    FT**2*FU**3*CC032 + FV*CC100 + FV*FT*CC101 + 
     $     FV*FT**2*CC102 + FV*FT**3*CC103 + 
     $    FV*FT**4*CC104 + FV*FU*CC110 + FV*FT*FU*CC111 + 
     $     FV*FT**2*FU*CC112 + FV*FT**3*FU*CC113 + 
     $    FV*FU**2*CC120 + FV*FT*FU**2*CC121 + FV*FT**2*FU**2*CC122 + 
     $     FV**2*CC200 + FV**2*FT*CC201 + 
     $    FV**2*FT**2*CC202 + FV**2*FT**3*CC203 + FV**2*FU*CC210 + 
     $     FV**2*FT*FU*CC211 + 
     $    FV**2*FT**2*FU*CC212 + FV**3*CC300 + FV**3*FT*CC301 + 
     $     FV**3*FT**2*CC302
c
c          F(i) = F(i) + E*FMN12
          Ex = Ex + E*FMN12*qwght(n)
          if(ldew) func(n)=func(n) + E*FMN12
c
c          If(MaxDrF.ge.1) then
          If(.true.) then
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
            dWdT = Two/((One + Tsig)**2)
            dTdR = ((F6*PI*PI)**F2o3)*(pX**F2o3)/tau_in
            dTdTau = -TauUEG/tau_in**2
c
            dFTdR = dWdT*dTdR
            dFTdTau=dWdT*dTdTau
c     
          dFMN12dFV = CC100 + FT*CC101 + FT**2*CC102 + FT**3*CC103 + 
     $     FT**4*CC104 + FU*CC110 + 
     $    FT*FU*CC111 + FT**2*FU*CC112 + FT**3*FU*CC113 + 
     $     FU**2*CC120 + FT*FU**2*CC121 + 
     $    FT**2*FU**2*CC122 + Two*FV*CC200 + Two*FV*FT*CC201 + 
     $     Two*FV*FT**2*CC202 + Two*FV*FT**3*CC203 + 
     $    Two*FV*FU*CC210 + Two*FV*FT*FU*CC211+Two*FV*FT**2*FU*CC212+ 
     $     Three*FV**2*CC300 + Three*FV**2*FT*CC301 + 
     $    Three*FV**2*FT**2*CC302
c
          dFMN12dFU = CC010 + FT*CC011 + FT**2*CC012 + FT**3*CC013 + 
     $     FT**4*CC014 + Two*FU*CC020 + 
     $    Two*FT*FU*CC021 + Two*FT**2*FU*CC022 + Two*FT**3*FU*CC023 + 
     $     Three*FU**2*CC030 + Three*FT*FU**2*CC031 + 
     $    Three*FT**2*FU**2*CC032 + FV*CC110 + FV*FT*CC111 + 
     $     FV*FT**2*CC112 + FV*FT**3*CC113 + 
     $    Two*FV*FU*CC120 + Two*FV*FT*FU*CC121+Two*FV*FT**2*FU*CC122+ 
     $     FV**2*CC210 + FV**2*FT*CC211 + 
     $    FV**2*FT**2*CC212
c    
          dFMN12dFT = CC001 + Two*FT*CC002 + Three*FT**2*CC003 + 
     $     Four*FT**3*CC004 + Five*FT**4*CC005 + FU*CC011 + 
     $    Two*FT*FU*CC012 + Three*FT**2*FU*CC013+Four*FT**3*FU*CC014+ 
     $     FU**2*CC021 + Two*FT*FU**2*CC022 + 
     $    Three*FT**2*FU**2*CC023 + FU**3*CC031 + Two*FT*FU**3*CC032 + 
     $     FV*CC101 + Two*FV*FT*CC102 + 
     $    Three*FV*FT**2*CC103 + Four*FV*FT**3*CC104 + FV*FU*CC111 + 
     $     Two*FV*FT*FU*CC112 + Three*FV*FT**2*FU*CC113 + 
     $    FV*FU**2*CC121 + Two*FV*FT*FU**2*CC122 + FV**2*CC201 + 
     $     Two*FV**2*FT*CC202 + 
     $    Three*FV**2*FT**2*CC203 + FV**2*FU*CC211 + 
     $     Two*FV**2*FT*FU*CC212 + FV**3*CC301 + Two*FV**3*FT*CC302
c
        dFMN12dR = dFMN12dFV*dFVdR+dFMN12dFU*US*SR+dFMN12dFT*dFTdR
        dFMN12dG = dFMN12dFU*US*SG
        dFMN12dT = dFMN12dFT*dFTdTau
c
            Amat(n,2) = Amat(n,2)   + ER*FMN12 +
     $                       E*dFMN12dR
            Cmat(n,3)=  Cmat(n,3)   + E*dFMN12dG/(Two*GamX)
            Mmat(n,2)=  Mmat(n,2)   + E*dFMN12dT
            endIf
         endIf
#endif         
 10    Continue
       
      Return
      End
      Subroutine xc_xmn12_d2()
      call errquit(' not coded ',0,0)
      return
      end
