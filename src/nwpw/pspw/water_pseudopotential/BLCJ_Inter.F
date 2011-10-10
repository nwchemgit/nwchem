**********************************************************
*                                                        *
*              BLCJ_Inter module                         *
*                                                        *
*          BLCJ Electron-Water pseudopotential           *
*                                                        *
*          used by BLCJ_Water module                     *
*                                                        *
*    -- developed by Eric J. Bylaska on February 6, 2001 *
*                                                        *
**********************************************************



*     **********************************
*     *	                               *
*     *        BLCJ_Intermolecular     *
*     *                                *
*     **********************************

      real*8 function BLCJ_Intermolecular(ro, rh1, rh2,
     >                                    rpo,rph1,rph2)
      implicit none
      real*8 ro(3),  rh1(3),  rh2(3)
      real*8 rpo(3), rph1(3), rph2(3)

*     **** BLCJ Intermolecular parameters ****
      real*8 q1,q2,q3,DELTA
      real*8 AH, alphaH, AO,alphaO
      parameter (q1=0.6d0,q2=0.6d0,q3=-1.2d0)
      parameter (DELTA=0.22183756d0)
      parameter (AH=1.0077d0, alphaH=1.7360d0)
      parameter (AO=5110.7d0, alphaO=2.6301d0)


*     **** local variables ****
      real*8 r3(3),rp3(3)
      real*8 r11,r12,r13
      real*8 r21,r22,r23
      real*8 r31,r32,r33
      real*8 roo,ro1,ro2,r1o,r2o
      real*8 ewater_water

*     **** external functions ****
      real*8   Inter_Morse, Inter_Dispersion
      external Inter_Morse, Inter_Dispersion
  
*     **** distance values ****
      r3(1) = ro(1) + DELTA*(rh1(1) + rh2(1) - 2.0d0*ro(1))
      r3(2) = ro(2) + DELTA*(rh1(2) + rh2(2) - 2.0d0*ro(2))
      r3(3) = ro(3) + DELTA*(rh1(3) + rh2(3) - 2.0d0*ro(3))

      rp3(1) = rpo(1) + DELTA*(rph1(1) + rph2(1) - 2.0d0*rpo(1))
      rp3(2) = rpo(2) + DELTA*(rph1(2) + rph2(2) - 2.0d0*rpo(2))
      rp3(3) = rpo(3) + DELTA*(rph1(3) + rph2(3) - 2.0d0*rpo(3))

      r11 = (rh1(1)-rph1(1))**2+(rh1(2)-rph1(2))**2+(rh1(3)-rph1(3))**2
      r12 = (rh1(1)-rph2(1))**2+(rh1(2)-rph2(2))**2+(rh1(3)-rph2(3))**2
      r13 = (rh1(1)-rp3(1))**2 +(rh1(2)-rp3(2))**2 +(rh1(3)-rp3(3))**2

      r21 = (rh2(1)-rph1(1))**2+(rh2(2)-rph1(2))**2+(rh2(3)-rph1(3))**2
      r22 = (rh2(1)-rph2(1))**2+(rh2(2)-rph2(2))**2+(rh2(3)-rph2(3))**2
      r23 = (rh2(1)-rp3(1))**2 +(rh2(2)-rp3(2))**2 +(rh2(3)-rp3(3))**2

      r31 = (r3(1)-rph1(1))**2+(r3(2)-rph1(2))**2+(r3(3)-rph1(3))**2
      r32 = (r3(1)-rph2(1))**2+(r3(2)-rph2(2))**2+(r3(3)-rph2(3))**2
      r33 = (r3(1)-rp3(1))**2 +(r3(2)-rp3(2))**2 +(r3(3)-rp3(3))**2

      roo = (ro(1)- rpo(1))**2+(ro(2)- rpo(2))**2+(ro(3)- rpo(3))**2

      ro1 = (ro(1)-rph1(1))**2+(ro(2)-rph1(2))**2+(ro(3)-rph1(3))**2
      ro2 = (ro(1)-rph2(1))**2+(ro(2)-rph2(2))**2+(ro(3)-rph2(3))**2

      r1o = (rh1(1)-rpo(1))**2+(rh1(2)-rpo(2))**2+(rh1(3)-rpo(3))**2
      r2o = (rh2(1)-rpo(1))**2+(rh2(2)-rpo(2))**2+(rh2(3)-rpo(3))**2

      r11 = dsqrt(r11)
      r12 = dsqrt(r12)
      r13 = dsqrt(r13)

      r21 = dsqrt(r21)
      r22 = dsqrt(r22)
      r23 = dsqrt(r23)

      r31 = dsqrt(r31)
      r32 = dsqrt(r32)
      r33 = dsqrt(r33)

      roo = dsqrt(roo)
      ro1 = dsqrt(ro1)
      ro2 = dsqrt(ro2)
      r1o = dsqrt(r1o)
      r2o = dsqrt(r2o)

*     **** Coulomb interaction ****
      ewater_water = q1*q1/r11 + q1*q2/r12 + q1*q3/r13
     >             + q2*q1/r21 + q2*q2/r22 + q2*q3/r23
     >             + q3*q1/r31 + q3*q2/r32 + q3*q3/r33

*     **** Exponential repulsion ****
      ewater_water = ewater_water
     >             + AH*dexp(-alphaH*r11)
     >             + AH*dexp(-alphaH*r12)
     >             + AH*dexp(-alphaH*r21)
     >             + AH*dexp(-alphaH*r22)
     >             + AO*dexp(-alphaO*roo)

*     **** Morse ****
      ewater_water = ewater_water
     >             + Inter_Morse(ro1) + Inter_Morse(ro2)
     >             + Inter_Morse(r1o) + Inter_Morse(r2o)

*     **** Dispersion ****
      ewater_water = ewater_water
     >             + Inter_Dispersion(roo)

      BLCJ_Intermolecular = ewater_water
      return
      end




*     **********************************
*     *	                               *
*     * BLCJ_Intermolecular_Fwater     *
*     *                                *
*     **********************************

      subroutine  BLCJ_Intermolecular_Fwater(ro, rh1, rh2, 
     >                                       fo, fh1, fh2,
     >                                       rpo,rph1,rph2,
     >                                       fpo,fph1,fph2)
      implicit none
      real*8 ro(3),  rh1(3),  rh2(3)
      real*8 fo(3),  fh1(3),  fh2(3)
      real*8 rpo(3), rph1(3), rph2(3)
      real*8 fpo(3), fph1(3), fph2(3)

*     **** BLCJ Intermolecular parameters ****
      real*8 q1,q2,q3,DELTA
      real*8 AH, alphaH, AO,alphaO
      parameter (q1=0.6d0,q2=0.6d0,q3=-1.2d0)
      parameter (DELTA=0.22183756d0)
      parameter (AH=1.0077d0, alphaH=1.7360d0)
      parameter (AO=5110.7d0, alphaO=2.6301d0)


*     **** local variables ****
      real*8 r3(3),rp3(3)
      real*8 r11,r12,r13
      real*8 r21,r22,r23
      real*8 r31,r32,r33
      real*8 roo,ro1,ro2,r1o,r2o
      real*8 tmp,d3do

*     **** external functions ****
      real*8   DInter_Morse, DInter_Dispersion
      external DInter_Morse, DInter_Dispersion
  
      d3do = (1.0d0-2.0d0*DELTA)

*     **** distance values ****
      r3(1) = ro(1) + DELTA*(rh1(1) + rh2(1) - 2.0d0*ro(1))
      r3(2) = ro(2) + DELTA*(rh1(2) + rh2(2) - 2.0d0*ro(2))
      r3(3) = ro(3) + DELTA*(rh1(3) + rh2(3) - 2.0d0*ro(3))

      rp3(1) = rpo(1) + DELTA*(rph1(1) + rph2(1) - 2.0d0*rpo(1))
      rp3(2) = rpo(2) + DELTA*(rph1(2) + rph2(2) - 2.0d0*rpo(2))
      rp3(3) = rpo(3) + DELTA*(rph1(3) + rph2(3) - 2.0d0*rpo(3))

      r11 = (rh1(1)-rph1(1))**2+(rh1(2)-rph1(2))**2+(rh1(3)-rph1(3))**2
      r12 = (rh1(1)-rph2(1))**2+(rh1(2)-rph2(2))**2+(rh1(3)-rph2(3))**2
      r13 = (rh1(1)-rp3(1))**2 +(rh1(2)-rp3(2))**2 +(rh1(3)-rp3(3))**2

      r21 = (rh2(1)-rph1(1))**2+(rh2(2)-rph1(2))**2+(rh2(3)-rph1(3))**2
      r22 = (rh2(1)-rph2(1))**2+(rh2(2)-rph2(2))**2+(rh2(3)-rph2(3))**2
      r23 = (rh2(1)-rp3(1))**2 +(rh2(2)-rp3(2))**2 +(rh2(3)-rp3(3))**2

      r31 = (r3(1)-rph1(1))**2+(r3(2)-rph1(2))**2+(r3(3)-rph1(3))**2
      r32 = (r3(1)-rph2(1))**2+(r3(2)-rph2(2))**2+(r3(3)-rph2(3))**2
      r33 = (r3(1)-rp3(1))**2 +(r3(2)-rp3(2))**2 +(r3(3)-rp3(3))**2

      roo = (ro(1)- rpo(1))**2+(ro(2)- rpo(2))**2+(ro(3)- rpo(3))**2

      ro1 = (ro(1)-rph1(1))**2+(ro(2)-rph1(2))**2+(ro(3)-rph1(3))**2
      ro2 = (ro(1)-rph2(1))**2+(ro(2)-rph2(2))**2+(ro(3)-rph2(3))**2

      r1o = (rh1(1)-rpo(1))**2+(rh1(2)-rpo(2))**2+(rh1(3)-rpo(3))**2
      r2o = (rh2(1)-rpo(1))**2+(rh2(2)-rpo(2))**2+(rh2(3)-rpo(3))**2

      r11 = dsqrt(r11)
      r12 = dsqrt(r12)
      r13 = dsqrt(r13)

      r21 = dsqrt(r21)
      r22 = dsqrt(r22)
      r23 = dsqrt(r23)

      r31 = dsqrt(r31)
      r32 = dsqrt(r32)
      r33 = dsqrt(r33)

      roo = dsqrt(roo)
      ro1 = dsqrt(ro1)
      ro2 = dsqrt(ro2)
      r1o = dsqrt(r1o)
      r2o = dsqrt(r2o)

*     **** Coulomb interaction ****
      tmp = -q1*q1/r11**2
      fh1(1)  = fh1(1) - ((rh1(1)-rph1(1))/r11) * tmp
      fh1(2)  = fh1(2) - ((rh1(2)-rph1(2))/r11) * tmp
      fh1(3)  = fh1(3) - ((rh1(3)-rph1(3))/r11) * tmp

      fph1(1)  = fph1(1) - ((rph1(1)-rh1(1))/r11) * tmp
      fph1(2)  = fph1(2) - ((rph1(2)-rh1(2))/r11) * tmp
      fph1(3)  = fph1(3) - ((rph1(3)-rh1(3))/r11) * tmp

      tmp = -q1*q2/r12**2
      fh1(1)  = fh1(1) - ((rh1(1)-rph2(1))/r12) * tmp
      fh1(2)  = fh1(2) - ((rh1(2)-rph2(2))/r12) * tmp
      fh1(3)  = fh1(3) - ((rh1(3)-rph2(3))/r12) * tmp

      fph2(1)  = fph2(1) - ((rph2(1)-rh1(1))/r12) * tmp
      fph2(2)  = fph2(2) - ((rph2(2)-rh1(2))/r12) * tmp
      fph2(3)  = fph2(3) - ((rph2(3)-rh1(3))/r12) * tmp

      tmp = -q1*q3/r13**2
      fh1(1)  = fh1(1) - ((rh1(1)-rp3(1))/r13) * tmp * DELTA
      fh1(2)  = fh1(2) - ((rh1(2)-rp3(2))/r13) * tmp * DELTA
      fh1(3)  = fh1(3) - ((rh1(3)-rp3(3))/r13) * tmp * DELTA

      fpo(1)  = fpo(1) - ((rp3(1)-rh1(1))/r13) * tmp * d3do
      fpo(2)  = fpo(2) - ((rp3(2)-rh1(2))/r13) * tmp * d3do
      fpo(3)  = fpo(3) - ((rp3(3)-rh1(3))/r13) * tmp * d3do


      tmp = -q2*q1/r21**2
      fh2(1)  = fh2(1) - ((rh2(1)-rph1(1))/r21) * tmp
      fh2(2)  = fh2(2) - ((rh2(2)-rph1(2))/r21) * tmp
      fh2(3)  = fh2(3) - ((rh2(3)-rph1(3))/r21) * tmp

      fph1(1)  = fph1(1) - ((rph1(1)-rh2(1))/r21) * tmp
      fph1(2)  = fph1(2) - ((rph1(2)-rh2(2))/r21) * tmp
      fph1(3)  = fph1(3) - ((rph1(3)-rh2(3))/r21) * tmp

      tmp = -q2*q2/r21**2
      fh2(1)  = fh2(1) - ((rh2(1)-rph2(1))/r22) * tmp
      fh2(2)  = fh2(2) - ((rh2(2)-rph2(2))/r22) * tmp
      fh2(3)  = fh2(3) - ((rh2(3)-rph2(3))/r22) * tmp

      fph2(1)  = fph2(1) - ((rph2(1)-rh2(1))/r22) * tmp
      fph2(2)  = fph2(2) - ((rph2(2)-rh2(2))/r22) * tmp
      fph2(3)  = fph2(3) - ((rph2(3)-rh2(3))/r22) * tmp

      tmp = -q2*q3/r23**2
      fh2(1)  = fh2(1) - ((rh2(1)-rp3(1))/r23) * tmp * DELTA
      fh2(2)  = fh2(2) - ((rh2(2)-rp3(2))/r23) * tmp * DELTA
      fh2(3)  = fh2(3) - ((rh2(3)-rp3(3))/r23) * tmp * DELTA

      fpo(1)  = fpo(1) - ((rp3(1)-rh2(1))/r23) * tmp * d3do
      fpo(2)  = fpo(2) - ((rp3(2)-rh2(2))/r23) * tmp * d3do
      fpo(3)  = fpo(3) - ((rp3(3)-rh2(3))/r23) * tmp * d3do


      tmp = -q3*q1/r31**2
      fo(1)  = fo(1) - ((r3(1)-rph1(1))/r31) * tmp * d3do
      fo(2)  = fo(2) - ((r3(2)-rph1(2))/r31) * tmp * d3do
      fo(3)  = fo(3) - ((r3(3)-rph1(3))/r31) * tmp * d3do

      fph1(1)  = fph1(1) - ((rph1(1)-r3(1))/r31) * tmp * DELTA
      fph1(2)  = fph1(2) - ((rph1(2)-r3(2))/r31) * tmp * DELTA
      fph1(3)  = fph1(3) - ((rph1(3)-r3(3))/r31) * tmp * DELTA

      tmp = -q3*q2/r32**2
      fo(1)  = fo(1) - ((r3(1)-rph2(1))/r32) * tmp * d3do
      fo(2)  = fo(2) - ((r3(2)-rph2(2))/r32) * tmp * d3do
      fo(3)  = fo(3) - ((r3(3)-rph2(3))/r32) * tmp * d3do

      fph2(1)  = fph2(1) - ((rph2(1)-r3(1))/r32) * tmp * DELTA
      fph2(2)  = fph2(2) - ((rph2(2)-r3(2))/r32) * tmp * DELTA
      fph2(3)  = fph2(3) - ((rph2(3)-r3(3))/r32) * tmp * DELTA

      tmp = -q3*q3/r33**2
      fo(1)  = fo(1) - ((r3(1)-rp3(1))/r33) * tmp * d3do
      fo(2)  = fo(2) - ((r3(2)-rp3(2))/r33) * tmp * d3do
      fo(3)  = fo(3) - ((r3(3)-rp3(3))/r33) * tmp * d3do

      fpo(1)  = fpo(1) - ((rp3(1)-r3(1))/r33) * tmp * d3do
      fpo(2)  = fpo(2) - ((rp3(2)-r3(2))/r33) * tmp * d3do
      fpo(3)  = fpo(3) - ((rp3(3)-r3(3))/r33) * tmp * d3do



*     **** Exponential repulsion ****

      tmp = -alphaH*AH*dexp(-alphaH*r11)
      fh1(1)  = fh1(1) - ((rh1(1)-rph1(1))/r11) * tmp
      fh1(2)  = fh1(2) - ((rh1(2)-rph1(2))/r11) * tmp
      fh1(3)  = fh1(3) - ((rh1(3)-rph1(3))/r11) * tmp

      fph1(1)  = fph1(1) - ((rph1(1)-rh1(1))/r11) * tmp
      fph1(2)  = fph1(2) - ((rph1(2)-rh1(2))/r11) * tmp
      fph1(3)  = fph1(3) - ((rph1(3)-rh1(3))/r11) * tmp

      tmp = -alphaH*AH*dexp(-alphaH*r12)
      fh1(1)  = fh1(1) - ((rh1(1)-rph2(1))/r12) * tmp
      fh1(2)  = fh1(2) - ((rh1(2)-rph2(2))/r12) * tmp
      fh1(3)  = fh1(3) - ((rh1(3)-rph2(3))/r12) * tmp

      fph2(1)  = fph2(1) - ((rph2(1)-rh1(1))/r12) * tmp
      fph2(2)  = fph2(2) - ((rph2(2)-rh1(2))/r12) * tmp
      fph2(3)  = fph2(3) - ((rph2(3)-rh1(3))/r12) * tmp

      tmp = -alphaH*AH*dexp(-alphaH*r21)
      fh2(1)  = fh2(1) - ((rh2(1)-rph1(1))/r21) * tmp
      fh2(2)  = fh2(2) - ((rh2(2)-rph1(2))/r21) * tmp
      fh2(3)  = fh2(3) - ((rh2(3)-rph1(3))/r21) * tmp

      fph1(1)  = fph1(1) - ((rph1(1)-rh2(1))/r21) * tmp
      fph1(2)  = fph1(2) - ((rph1(2)-rh2(2))/r21) * tmp
      fph1(3)  = fph1(3) - ((rph1(3)-rh2(3))/r21) * tmp

      tmp = -alphaH*AH*dexp(-alphaH*r22)
      fh2(1)  = fh2(1) - ((rh2(1)-rph2(1))/r22) * tmp
      fh2(2)  = fh2(2) - ((rh2(2)-rph2(2))/r22) * tmp
      fh2(3)  = fh2(3) - ((rh2(3)-rph2(3))/r22) * tmp

      fph2(1)  = fph2(1) - ((rph2(1)-rh2(1))/r22) * tmp
      fph2(2)  = fph2(2) - ((rph2(2)-rh2(2))/r22) * tmp
      fph2(3)  = fph2(3) - ((rph2(3)-rh2(3))/r22) * tmp

      tmp = -alphaO*AO*dexp(-alphaO*roo)
      fo(1)  = fo(1) - ((ro(1)-rpo(1))/roo) * tmp
      fo(2)  = fo(2) - ((ro(2)-rpo(2))/roo) * tmp
      fo(3)  = fo(3) - ((ro(3)-rpo(3))/roo) * tmp

      fpo(1)  = fpo(1) - ((rpo(1)-ro(1))/roo) * tmp
      fpo(2)  = fpo(2) - ((rpo(2)-ro(2))/roo) * tmp
      fpo(3)  = fpo(3) - ((rpo(3)-ro(3))/roo) * tmp


*     **** Morse ****
      tmp = DInter_Morse(ro1)
      fo(1)  = fo(1) - ((ro(1)-rph1(1))/ro1) * tmp
      fo(2)  = fo(2) - ((ro(2)-rph1(2))/ro1) * tmp
      fo(3)  = fo(3) - ((ro(3)-rph1(3))/ro1) * tmp

      fph1(1) = fph1(1) - ((rph1(1)-ro(1))/ro1) * tmp
      fph1(2) = fph1(2) - ((rph1(2)-ro(2))/ro1) * tmp
      fph1(3) = fph1(3) - ((rph1(3)-ro(3))/ro1) * tmp

      tmp = DInter_Morse(ro2)
      fo(1)  = fo(1) - ((ro(1)-rph2(1))/ro1) * tmp
      fo(2)  = fo(2) - ((ro(2)-rph2(2))/ro1) * tmp
      fo(3)  = fo(3) - ((ro(3)-rph2(3))/ro1) * tmp

      fph2(1) = fph2(1) - ((rph2(1)-ro(1))/ro1) * tmp
      fph2(2) = fph2(2) - ((rph2(2)-ro(2))/ro1) * tmp
      fph2(3) = fph2(3) - ((rph2(3)-ro(3))/ro1) * tmp

      tmp = DInter_Morse(r1o)
      fpo(1)  = fpo(1) - ((rpo(1)-rh1(1))/r1o) * tmp
      fpo(2)  = fpo(2) - ((rpo(2)-rh1(2))/r1o) * tmp
      fpo(3)  = fpo(3) - ((rpo(3)-rh1(3))/r1o) * tmp

      fh1(1) = fh1(1) - ((rh1(1)-rpo(1))/r1o) * tmp
      fh1(2) = fh1(2) - ((rh1(2)-rpo(2))/r1o) * tmp
      fh1(3) = fh1(3) - ((rh1(3)-rpo(3))/r1o) * tmp

      tmp = DInter_Morse(r2o)
      fpo(1)  = fpo(1) - ((rpo(1)-rh2(1))/r2o) * tmp
      fpo(2)  = fpo(2) - ((rpo(2)-rh2(2))/r2o) * tmp
      fpo(3)  = fpo(3) - ((rpo(3)-rh2(3))/r2o) * tmp

      fh2(1) = fh2(1) - ((rh2(1)-rpo(1))/r2o) * tmp
      fh2(2) = fh2(2) - ((rh2(2)-rpo(2))/r2o) * tmp
      fh2(3) = fh2(3) - ((rh2(3)-rpo(3))/r2o) * tmp


*     **** Dispersion ****
      tmp = DInter_Dispersion(roo)
      fo(1)  = fo(1) - ((ro(1)-rpo(1))/roo) * tmp
      fo(2)  = fo(2) - ((ro(2)-rpo(2))/roo) * tmp
      fo(3)  = fo(3) - ((ro(3)-rpo(3))/roo) * tmp

      fpo(1) = fpo(1) - ((rpo(1)-ro(1))/roo) * tmp
      fpo(2) = fpo(2) - ((rpo(2)-ro(2))/roo) * tmp
      fpo(3) = fpo(3) - ((rpo(3)-ro(3))/roo) * tmp

      return
      end



      
      real*8 FUNCTION Inter_Morse(r)
      implicit none
      real*8 r

*     **** BLCJ Intermolecular parameters ****
      real*8 A, alpha, Rm
      parameter (A=0.003066d0, alpha=3.89556d0, Rm=3.0951d0)

      Inter_Morse = A*dexp(-alpha*(r-Rm))*(dexp(-alpha*(r-Rm))-2.0d0)
      return
      end

      real*8 FUNCTION DInter_Morse(r)
      implicit none
      real*8 r

*     **** BLCJ Intermolecular parameters ****
      real*8 A, alpha, Rm
      parameter (A=0.003066d0, alpha=3.89556d0, Rm=3.0951d0)

      DInter_Morse = -alpha*A*dexp(-alpha*(r-Rm))
     >             *   (2.0d0*dexp(-alpha*(r-Rm))-2.0d0)
      return
      end


      real*8 FUNCTION Inter_Dispersion(r)
      implicit none
      real*8 r

*     **** BLCJ Intermolecular parameters ****
      real*8 C6,C8,C10
      parameter (C6=62.437d0, C8=1343.7d0, C10=50050d0)    

*     **** external functions ****
      real*8 Fd,Gd
      external Fd,Gd      

      Inter_Dispersion = Fd(r)*( C6 *(Gd( 6,r)/r)**6
     >                         + C8 *(Gd( 8,r)/r)**8
     >                         + C10*(Gd(10,r)/r)**10 )

      return 
      end

      real*8 FUNCTION DInter_Dispersion(r)
      implicit none
      real*8 r

*     **** BLCJ Intermolecular parameters ****
      real*8 C6,C8,C10
      parameter (C6=62.437d0, C8=1343.7d0, C10=50050d0)    

*     **** external functions ****
      real*8   Fd,Gd,DFd,DGd
      external Fd,Gd,DFd,DGd      

      DInter_Dispersion = DFd(r)*( C6 *(Gd( 6,r)/r)**6
     >                         + C8 *(Gd( 8,r)/r)**8
     >                         + C10*(Gd(10,r)/r)**10 )
     >  + Fd(r)*( 6*C6 *(Gd(6,r)/r)**5*(DGd(6,r)/r -Gd(6,r)/r**2)
     >         +  8*C8 *(Gd(8,r)/r)**7*(DGd(8,r)/r -Gd(8,r)/r**2)
     >        + 10*C10*(Gd(10,r)/r)**9*(DGd(10,r)/r-Gd(10,r)/r**2 ))

      return 
      end



      real*8 FUNCTION Fd(r)
      implicit none
      real*8 r

*     **** BLCJ Intermolecular parameters ****  
      real*8 c,beta
      parameter (c=0.94835d0, beta= 2.326d0)
  
      Fd = 1.0d0 - (c*r)**beta * dexp(-c*r)
      return
      end

      real*8 FUNCTION DFd(r)
      implicit none
      real*8 r

*     **** BLCJ Intermolecular parameters ****  
      real*8 c,beta
      parameter (c=0.94835d0, beta= 2.326d0)
  
      DFd = (-beta/r + c) * (c*r)**beta * dexp(-c*r)
      return
      end




      real*8 FUNCTION Gd(n,r)
      implicit none
      integer n
      real*8  r

*     **** BLCJ Intermolecular parameters ****
      real*8 a,b
      parameter (a=1.9915d0, b=0.98030d0)

      real*8 x
      x = dble(n)

      Gd = 1.0d0 - dexp(-(a*r/x + b*r*r/dsqrt(x)))
      return
      end

      real*8 FUNCTION DGd(n,r)
      implicit none
      integer n
      real*8  r

*     **** BLCJ Intermolecular parameters ****
      real*8 a,b
      parameter (a=1.9915d0, b=0.98030d0)

      real*8 x
      x = dble(n)

      DGd = (a/x +2*r*b/dsqrt(x))*dexp(-(a*r/x + b*r*r/dsqrt(x)))
      return
      end
c $Id$
