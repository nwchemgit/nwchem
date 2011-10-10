*     **********************************
*     *	                               *
*     *       BLCJ_ion_Polarization    *
*     *                                *
*     **********************************
*
*   This routine calculate the Polarization contribution
*  to the BLCJ water psp.
*
      subroutine BLCJ_ion_Polarization(RO,R1,R2,
     >                        rion,zv,
     >                        ep)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      real*8 rion(3)
      real*8 zv
      real*8 ep

*     **** BLCJ Polarization Parameters ****
#include "waterpsp_param.fh"


*     **** local variables ****
      real*8 xo,yo,zo,d4
     
     
      xo = rion(1)-RO(1)
      yo = rion(2)-RO(2)
      zo = rion(3)-RO(3)

      d4 = (xo*xo + yo*yo + zo*zo + Rp*Rp)**2
        
      ep = 0.5d0*zv*alpha/d4    

      return 
      end


 
 
*     **********************************
*     *	                               *
*     *   BLCJ_ion_Polarization_Fion   *
*     *                                *
*     **********************************

*   This routine calculate the Polarization contribution
*  to the BLCJ water psp.
*

      subroutine BLCJ_ion_Polarization_Fion(RO,R1,R2,
     >                        rion,zv,
     >                        fion)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      real*8 rion(3)
      real*8 zv
      real*8 fion(3)

*     **** Parameters ****
#include "waterpsp_param.fh"


*     **** local variables ****
      real*8 xo,yo,zo,d2,d4
      real*8 dep     
     
      xo = rion(1)-RO(1)
      yo = rion(2)-RO(2)
      zo = rion(3)-RO(3)
c      do = dsqrt(xo*xo + yo*yo + zo*zo)

      d2 = (xo*xo + yo*yo + zo*zo + Rp*Rp)
      d4 = (d2)**2
        
      dep = -2.0d0*alpha*zv/(d4*d2)    

      fion(1) = -(xo)*dep
      fion(2) = -(yo)*dep
      fion(3) = -(zo)*dep

      return 
      end


*     **********************************
*     *	                               *
*     *  BLCJ_ion_Polarization_Fwater  *
*     *                                *
*     **********************************

*   This routine calculate the Polarization contribution
*  to the BLCJ water psp.
*

      subroutine BLCJ_ion_Polarization_Fwater(RO,R1,R2,
     >                        rion,zv,
     >                        fo)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      real*8 rion(3)
      real*8 zv
      real*8 fo(3)

*     **** Parameters ****
#include "waterpsp_param.fh"


*     **** local variables ****
      real*8 xo,yo,zo,d2,d4
      real*8 dep     
     
      xo = rion(1)-RO(1)
      yo = rion(2)-RO(2)
      zo = rion(3)-RO(3)
c      do = dsqrt(xo*xo + yo*yo + zo*zo)

      d2 = (xo*xo + yo*yo + zo*zo + Rp*Rp)
      d4 = (d2)**2
        
      dep = -2.0d0*alpha*zv/(d4*d2)    

      fo(1) = (xo)*dep
      fo(2) = (yo)*dep
      fo(3) = (zo)*dep

      return 
      end
c $Id$
