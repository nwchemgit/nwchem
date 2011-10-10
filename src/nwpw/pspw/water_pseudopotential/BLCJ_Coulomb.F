*     **********************************
*     *	                               *
*     *          BLCJ_ion_Coulomb      *
*     *                                *
*     **********************************

*   This routine calculate the Coulomb contribution
*  to the BLCJ water psp.
*

      subroutine BLCJ_ion_Coulomb(RO,R1,R2,
     >                        rion,zv,
     >                        ec)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      real*8 rion(3),zv
      real*8 ec

*     **** Parameters ****
#include "waterpsp_param.fh"

*     **** local variables ****
      real*8 R3(3)
      real*8 x1,y1,z1,d1
      real*8 x2,y2,z2,d2
      real*8 x3,y3,z3,d3

*     **** define R3 = center of negative charge ***
      R3(1) = RO(1) + DELTA*(R1(1)+R2(1)-2.0d0*RO(1))
      R3(2) = RO(2) + DELTA*(R1(2)+R2(2)-2.0d0*RO(2))
      R3(3) = RO(3) + DELTA*(R1(3)+R2(3)-2.0d0*RO(3))


      x1 = rion(1)-R1(1)
      x2 = rion(1)-R2(1)
      x3 = rion(1)-R3(1)

      y1 = rion(2)-R1(2)
      y2 = rion(2)-R2(2)
      y3 = rion(2)-R3(2)

      z1 = rion(3)-R1(3)
      z2 = rion(3)-R2(3)
      z3 = rion(3)-R3(3)

      d1 = dsqrt(x1*x1 + y1*y1 + z1*z1)
      d2 = dsqrt(x2*x2 + y2*y2 + z2*z2)
      d3 = dsqrt(x3*x3 + y3*y3 + z3*z3)

      ec = zv*( q1/d1 
     >        + q2/d2 
     >        + q3/d3)    

      return 
      end



*     **********************************
*     *	                               *
*     *     BLCJ_ion_Coulomb_Fion      *
*     *                                *
*     **********************************

*   This routine calculate the Coulomb contribution
*  to the BLCJ water psp.
*

      subroutine BLCJ_ion_Coulomb_Fion(RO,R1,R2,
     >                        rion,zv,fion)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      real*8 rion(3),zv
      real*8 fion(3)

*     **** Parameters ****
#include "waterpsp_param.fh"


*     **** local variables ****
      real*8 R3(3)
      real*8 x1,y1,z1,d1
      real*8 x2,y2,z2,d2
      real*8 x3,y3,z3,d3
      real*8 dec1,dec2,dec3

     
*     **** define R3 = center of negative charge ***
      R3(1) = RO(1) + DELTA*(R1(1)+R2(1)-2*RO(1))
      R3(2) = RO(2) + DELTA*(R1(2)+R2(2)-2*RO(2))
      R3(3) = RO(3) + DELTA*(R1(3)+R2(3)-2*RO(3))


      x1 = rion(1)-R1(1)
      x2 = rion(1)-R2(1)
      x3 = rion(1)-R3(1)

      y1 = rion(2)-R1(2)
      y2 = rion(2)-R2(2)
      y3 = rion(2)-R3(2)

      z1 = rion(3)-R1(3)
      z2 = rion(3)-R2(3)
      z3 = rion(3)-R3(3)

      d1 = dsqrt(x1*x1 + y1*y1 + z1*z1)
      d2 = dsqrt(x2*x2 + y2*y2 + z2*z2)
      d3 = dsqrt(x3*x3 + y3*y3 + z3*z3)


      dec1 = -zv*q1/(d1*d1) 
      dec2 = -zv*q2/(d2*d2)
      dec3 = -zv*q3/(d3*d3)

      fion(1) = -(x1/d1)*dec1 - (x2/d2)*dec2 - (x3/d3)*dec3
      fion(2) = -(y1/d1)*dec1 - (y2/d2)*dec2 - (y3/d3)*dec3
      fion(3) = -(z1/d1)*dec1 - (z2/d2)*dec2 - (z3/d3)*dec3

      return 
      end


*     **********************************
*     *	                               *
*     *    BLCJ_ion_Coulomb_Fwater     *
*     *                                *
*     **********************************

*   This routine calculate the Coulomb contribution
*  to the BLCJ water psp.
*

      subroutine BLCJ_ion_Coulomb_Fwater(RO,R1,R2,
     >                        rion,zv,fo,f1,f2)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      real*8 rion(3),zv
      real*8 fo(3),f1(3),f2(3)

*     **** Parameters ****
#include "waterpsp_param.fh"

*     **** local variables ****
      real*8 R3(3)
      real*8 x1,y1,z1,d1
      real*8 x2,y2,z2,d2
      real*8 x3,y3,z3,d3
      real*8 dec1,dec2,dec3
      real*8 d3do
  
     
*     **** define R3 = center of negative charge ***
      R3(1) = RO(1) + DELTA*(R1(1)+R2(1)-2.0d0*RO(1))
      R3(2) = RO(2) + DELTA*(R1(2)+R2(2)-2.0d0*RO(2))
      R3(3) = RO(3) + DELTA*(R1(3)+R2(3)-2.0d0*RO(3))
      d3do  = (1.0d0-2.0d0*DELTA)

      x1 = rion(1)-R1(1)
      x2 = rion(1)-R2(1)
      x3 = rion(1)-R3(1)

      y1 = rion(2)-R1(2)
      y2 = rion(2)-R2(2)
      y3 = rion(2)-R3(2)

      z1 = rion(3)-R1(3)
      z2 = rion(3)-R2(3)
      z3 = rion(3)-R3(3)

      d1 = dsqrt(x1*x1 + y1*y1 + z1*z1)
      d2 = dsqrt(x2*x2 + y2*y2 + z2*z2)
      d3 = dsqrt(x3*x3 + y3*y3 + z3*z3)
      
      dec1 = -zv*q1/(d1*d1) 
      dec2 = -zv*q2/(d2*d2)
      dec3 = -zv*q3/(d3*d3)

      fo(1) =   d3do*(x3/d3)*dec3
      fo(2) =   d3do*(y3/d3)*dec3
      fo(3) =   d3do*(z3/d3)*dec3

      f1(1) =   (x1/d1)*dec1 + DELTA*(x3/d3)*dec3
      f1(2) =   (y1/d1)*dec1 + DELTA*(y3/d3)*dec3
      f1(3) =   (z1/d1)*dec1 + DELTA*(z3/d3)*dec3

      f2(1) =   (x2/d2)*dec2 + DELTA*(x3/d3)*dec3
      f2(2) =   (y2/d2)*dec2 + DELTA*(y3/d3)*dec3
      f2(3) =   (z2/d2)*dec2 + DELTA*(z3/d3)*dec3

      return 
      end


c $Id$
