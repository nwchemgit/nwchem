      subroutine xc_xmn15l(tol_rho, fac,lfac,nlfac,
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
      external  xc_xmn15l_os
      call xc_os2cs_mgga(xc_xmn15l_os,
     &     tol_rho, fac,lfac,nlfac, rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq, ipol, ex, 
     &     qwght, ldew, func)
      return
      end
      subroutine xc_xmn15l_os(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq,  ex, 
     &                     qwght, ldew, func, fact_cs)
      implicit none
#include "errquit.fh"
#include "dft2drv.fh"
c
c wrapper for n15 open-shell
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

      call xc_xmn15_0(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     N     nq,  ex,   qwght, ldew, func,
     I      1    ,fact_cs)
      return
      end
      subroutine xc_xmn15(tol_rho, fac,lfac,nlfac,
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
      external  xc_xmn15_os
      call xc_os2cs_mgga(xc_xmn15_os,
     &     tol_rho, fac,lfac,nlfac, rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq, ipol, ex, 
     &     qwght, ldew, func)
      return
      end
      subroutine xc_xmn15_os(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     n     nq,  ex, 
     &                     qwght, ldew, func, fact_cs)
      implicit none
#include "errquit.fh"
#include "dft2drv.fh"
c
c wrapper for n15 open-shell
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

      call xc_xmn15_0(tol_rho, fac,lfac,nlfac,
     r     rho, delrho, tau,
     &     amat, cmat, mmat,
     N     nq,  ex,   qwght, ldew, func,
     I      2    ,fact_cs)
      return
      end
      Subroutine xc_xmn15_0(tol_rho, fac,lfac,nlfac,
     R     rho, delrho, tau,
     &     Amat, Cmat, mmat,
     n     nq,  Ex, 
     &     qwght, ldew, func,
     I     ijzy, fact_cs)
c      Subroutine MN15x(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
**************************************************************
*                                                            *
*  MN15x evaluates the exchange part of MN15-L and MN15      *
*        functionals on the grid.                            *
*                                                            *
*  OUTPUT:                                                   *
*     F   - Functional values                                *
*     D1F - First derivatives with respect to RA, RB, GA, GB *
*              TA, TB                                        *
*                                                            *
*  INPUT:                                                    *
*       ijzy - 1 MN15-L                                      *
*       ijzy - 2 MN15                                        *
*                                                            *
*     RA,B   - Spin densities                                *
*     D1RA,B - Spin density gradients                        *
*     TA,B   - Spin kinetic energy densities                 *
*     NGrid  - number of grids                               *
*                                                            *
*  RP (09/12)                                                *
*                                                            *
**************************************************************
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
      double precision dFMN15dFT, dFMN15dFU, dFMN15dFV, dFMN15dG,
     d     dFMN15dR, dFMN15dT
      double precision dFTdR, dFTdTau, dFVdR
      double precision e,er,sr,ft,fmn15,fu,fv
c
c
c      INTEGER NGrid
c      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC      
c      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
c     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
      double precision One, Two, Three, Four, Five, Six, Seven, Eight,
     c     nine
      Save One,Two,Three,Four,Five,Six,Seven,Eight,Nine
      Data One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Five/5.0d0/,
     $  Six/6.0d0/,Seven/7.0d0/,Eight/8.0d0/,Nine/9.0d0/
c
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
c       MN15-L
c      
        CC000 =  0.670864162d0
        CC001 = -0.822003903d0
        CC002 = -1.022407046d0
        CC003 =  1.689460986d0
        CC004 = -0.00562032d0
        CC005 = -0.110293849d0
        CC010 =  0.972245178d0
        CC011 = -6.697641991d0
        CC012 = -4.322814495d0
        CC013 = -6.786641376d0
        CC014 = -5.687461462d0
        CC020 =  9.419643818d0
        CC021 =  11.83939406d0
        CC022 =  5.086951311d0
        CC023 =  4.302369948d0
        CC030 = -8.07344065d0
        CC031 =  2.429988978d0
        CC032 =  11.09485698d0
        CC100 =  1.247333909d0
        CC101 =  3.700485291d0
        CC102 =  0.867791614d0
        CC103 = -0.591190518d0
        CC104 = -0.295305435d0
        CC110 = -5.825759145d0
        CC111 =  2.537532196d0
        CC112 =  3.143390933d0
        CC113 =  2.939126332d0
        CC120 =  0.599342114d0
        CC121 =  2.241702738d0
        CC122 =  2.035713838d0
        CC200 = -1.525344043d0
        CC201 = -2.325875691d0
        CC202 =  1.141940663d0
        CC203 = -1.563165026d0
        CC210 =  7.882032871d0
        CC211 =  11.93400684d0
        CC212 =  9.852928303d0
        CC300 =  0.584030245d0
        CC301 = -0.720941131d0
        CC302 = -2.836037078d0

      Else If (ijzy.eq.2) then
c
c       MN15
c      
        CC000  =  0.073852235d0
        CC001  = -0.839976156d0
        CC002  = -3.082660125d0
        CC003  = -1.02881285d0
        CC004  = -0.811697255d0
        CC005  = -0.063404387d0
        CC010  =  2.54805518d0
        CC011  = -5.031578906d0
        CC012  =  0.31702159d0
        CC013  =  2.981868205d0
        CC014  = -0.749503735d0
        CC020  =  0.231825661d0
        CC021  =  1.261961411d0
        CC022  =  1.665920815d0
        CC023  =  7.483304941d0
        CC030  = -2.544245723d0
        CC031  =  1.384720031d0
        CC032  =  6.902569885d0
        CC100  =  1.657399451d0
        CC101  =  2.98526709d0
        CC102  =  6.89391326d0
        CC103  =  2.489813993d0
        CC104  =  1.454724691d0
        CC110  = -5.054324071d0
        CC111  =  2.35273334d0
        CC112  =  1.299104132d0
        CC113  =  1.203168217d0
        CC120  =  0.121595877d0
        CC121  =  8.048348238d0
        CC122  =  21.91203659d0
        CC200  = -1.852335832d0
        CC201  = -3.4722735d0
        CC202  = -1.564591493d0
        CC203  = -2.29578769d0
        CC210  =  3.666482991d0
        CC211  =  10.87074639d0
        CC212  =  9.696691388d0
        CC300  =  0.630701064d0
        CC301  = -0.505825216d0
        CC302  = -3.562354535d0
      else
         call errquit(' xc_xmn15: ijzy not valid ',ijzy,0)

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
c        If(RA(i).gt.Tiny.and.TA(i).gt.Tiny) then
          F3o5 = Three/Five
          F6   = Six
          F5o3 = Five/Three
c          pX    = RA(i)
c          GamX2 = D1RA(i)
          pX    = rho(n)
          Gamx2 = delrho(n,1)*delrho(n,1) +
     y         delrho(n,2)*delrho(n,2)  +
     z         delrho(n,3)*delrho(n,3)  
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E = Ax*pX**F4o3
          tau_in    = tau(n)
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3) 
          Tsig =TauUEG/Tau_In
          Wsig =(Tsig - One)/(Tsig + One)
c          
          FU = U
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FT = Wsig          
c
          FMN15 = CC000 + FT*CC001 + FT**2*CC002 + FT**3*CC003 + 
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
c          F(i) = F(i) + E*FMN15
          Ex = Ex + E*FMN15*qwght(n)*fact_cs
          if(ldew) func(n)=func(n) + E*FMN15*fact_cs
c
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
            dTdR = ((F6*PI*PI)**F2o3)*(pX**F2o3)/Tau_In
            dTdTau = -TauUEG/Tau_In**2
c
            dFTdR = dWdT*dTdR
            dFTdTau=dWdT*dTdTau
c     
          dFMN15dFV = CC100 + FT*CC101 + FT**2*CC102 + FT**3*CC103 + 
     $     FT**4*CC104 + FU*CC110 + 
     $    FT*FU*CC111 + FT**2*FU*CC112 + FT**3*FU*CC113 + 
     $     FU**2*CC120 + FT*FU**2*CC121 + 
     $    FT**2*FU**2*CC122 + Two*FV*CC200 + Two*FV*FT*CC201 + 
     $     Two*FV*FT**2*CC202 + Two*FV*FT**3*CC203 + 
     $    Two*FV*FU*CC210 + Two*FV*FT*FU*CC211+Two*FV*FT**2*FU*CC212+ 
     $     Three*FV**2*CC300 + Three*FV**2*FT*CC301 + 
     $    Three*FV**2*FT**2*CC302
c
          dFMN15dFU = CC010 + FT*CC011 + FT**2*CC012 + FT**3*CC013 + 
     $     FT**4*CC014 + Two*FU*CC020 + 
     $    Two*FT*FU*CC021 + Two*FT**2*FU*CC022 + Two*FT**3*FU*CC023 + 
     $     Three*FU**2*CC030 + Three*FT*FU**2*CC031 + 
     $    Three*FT**2*FU**2*CC032 + FV*CC110 + FV*FT*CC111 + 
     $     FV*FT**2*CC112 + FV*FT**3*CC113 + 
     $    Two*FV*FU*CC120 + Two*FV*FT*FU*CC121+Two*FV*FT**2*FU*CC122+ 
     $     FV**2*CC210 + FV**2*FT*CC211 + 
     $    FV**2*FT**2*CC212
c    
          dFMN15dFT = CC001 + Two*FT*CC002 + Three*FT**2*CC003 + 
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
        dFMN15dR = dFMN15dFV*dFVdR+dFMN15dFU*US*SR+dFMN15dFT*dFTdR
        dFMN15dG = dFMN15dFU*US*SG
        dFMN15dT = dFMN15dFT*dFTdTau
c
            Amat(n) = Amat(n) + ER*FMN15
     $                      + E*dFMN15dR
            Cmat(n)=  Cmat(n) + E*dFMN15dG/(Two*GamX)
            Mmat(n)=  Mmat(n) + E*dFMN15dT
            endIf
          endIf
   10   Continue
      Return
      End
      Subroutine xc_xmn15_d2()
      call errquit(' not coded ',0,0)
      return
      end
