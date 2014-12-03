*     **********************************
*     *                                *
*     *          BLCJ_LR_init          *
*     *                                *
*     **********************************

*   This reads in Gaussian Basis and initializes
*  variables for the LR contribution
*
      subroutine BLCJ_LR_init(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

* *** parameters ****
#include "waterpsp_param.fh"


*     **** local variables ***
      logical value
      integer taskid,MASTER
      parameter (MASTER=0)

      value = btdb_get(rtdb,'water_pseudopotential:DELTA',
     >                 mt_dbl,1,DELTA)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:Rcoul',
     >                 mt_dbl,1,Rcoul)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:q1',
     >                 mt_dbl,1,q1)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:q2',
     >                 mt_dbl,1,q2)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:q3',
     >                 mt_dbl,1,q3)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:Rp',
     >                 mt_dbl,1,Rp)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:Rpe',
     >                 mt_dbl,1,Rpe)
      value = value.and.
     >        btdb_get(rtdb,'water_pseudopotential:alpha',
     >                 mt_dbl,1,alpha)
        if (.not.value)
     >    call errquit('init_BLCJ_LR:rtdb error',0, RTDB_ERR)

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
        write(*,*)
        write(*,*) 
     >  'Charge and Polarization Parameters Water Pseudopotential'
        write(*,*) 
     >  '--------------------------------------------------------'
        write(*,*)
        write(*,*) 'qO:',q3,'  qH1:',q1,'  qH2:',q2,'   Rcoul:',Rcoul
        write(*,*) 'alpha:',alpha,'  Rp:',Rp, '  Rpe:',Rpe
      end if



      return
      end

*     **********************************
*     *	                               *
*     *          BLCJ_LR               *
*     *                                *
*     **********************************

*   This routine calculate the LR contribution 
*  to the BLCJ water psp.
*

      subroutine BLCJ_LR(RO,R1,R2,
     >                        n2ft3d,rgrid,
     >                        Vwpsp)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      integer n2ft3d
      real*8 rgrid(3,n2ft3d)
      real*8 Vwpsp(n2ft3d)
      real*8 util_erf
      external util_erf

*     *** parameters ****
#include "waterpsp_param.fh"


*     **** local variables ****
      integer ii,k
      real*8 Rwd(3,3),q(3)
      real*8 x,y,z,d,xx,rx
      real*8 sqrtpi

      sqrtpi = dsqrt(4.0d0*datan(1.0d0))

*     **** define center of charges ***
      q(1) = q1
      Rwd(1,1) = R1(1)
      Rwd(2,1) = R1(2)
      Rwd(3,1) = R1(3)

      q(2) = q2
      Rwd(1,2) = R2(1)
      Rwd(2,2) = R2(2)
      Rwd(3,2) = R2(3)
     
      q(3) = q3
      Rwd(1,3) = RO(1) + DELTA*(R1(1)+R2(1)-2*RO(1))
      Rwd(2,3) = RO(2) + DELTA*(R1(2)+R2(2)-2*RO(2))
      Rwd(3,3) = RO(3) + DELTA*(R1(3)+R2(3)-2*RO(3))


      do k=1,n2ft3d

*       **** LR Coulomb contribution ***
        do ii=1,3
           x = rgrid(1,k)-Rwd(1,ii)
           y = rgrid(2,k)-Rwd(2,ii)
           z = rgrid(3,k)-Rwd(3,ii)
           d = dsqrt(x*x + y*y + z*z)

           if (d.gt.1.0d-8) then
             Vwpsp(k) = Vwpsp(k) - q(ii)*util_erf(d/Rcoul)/d 
           else
             Vwpsp(k) = Vwpsp(k) - 2.0d0*q(ii)/(sqrtpi*Rcoul)
           end if

        end do

*       **** LR polarization contribution ***
        x = rgrid(1,k)-RO(1)
        y = rgrid(2,k)-RO(2)
        z = rgrid(3,k)-RO(3)
        d  = dsqrt(x*x + y*y + z*z)
        xx = d**4
        rx = (Rpe)**4
        if (xx .gt. 1.0d-8) then
           Vwpsp(k) = Vwpsp(k) - 0.5d0*alpha*util_erf(xx/rx)/xx
        else
           Vwpsp(k) = Vwpsp(k) - alpha/(sqrtpi*rx)
        end if

      end do

      return 
      end


*     **********************************
*     *	                               *
*     *          BLCJ_LR_Fwater        *
*     *                                *
*     **********************************

*   This routine calculates the water force of LR contribution 
*  to the BLCJ water psp.
*    

      subroutine BLCJ_LR_Fwater(RO,R1,R2,
     >                        n2ft3d,rgrid,rho,dv,
     >                        fo,f1,f2)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      integer n2ft3d
      real*8 rgrid(3,*),rho(*)
      real*8 dv
      real*8 fo(3),f1(3),f2(3)
      real*8 util_erf
      external util_erf

*     *** parameters ****
#include "waterpsp_param.fh"


*     **** local variables ****
      integer k
      real*8  R3(3)
      real*8  xo,yo,zo,d,depo,xx,rx
      real*8  x1,y1,z1,d1,dec1
      real*8  x2,y2,z2,d2,dec2
      real*8  x3,y3,z3,d3,dec3
      real*8  sqrtpi
      real*8  d3do

      sqrtpi = dsqrt(4.0d0*datan(1.0d0))

*     **** define R3 = center of negative charge ***
      R3(1) = RO(1) + DELTA*(R1(1)+R2(1)-2*RO(1))
      R3(2) = RO(2) + DELTA*(R1(2)+R2(2)-2*RO(2))
      R3(3) = RO(3) + DELTA*(R1(3)+R2(3)-2*RO(3))

      d3do=(1.0d0-2.0d0*DELTA)

      fo(1) = 0.0d0
      fo(2) = 0.0d0
      fo(3) = 0.0d0

      f1(1) = 0.0d0
      f1(2) = 0.0d0
      f1(3) = 0.0d0

      f2(1) = 0.0d0
      f2(2) = 0.0d0
      f2(3) = 0.0d0


      do k=1,n2ft3d

        xo = rgrid(1,k)-RO(1)
        x1 = rgrid(1,k)-R1(1)
        x2 = rgrid(1,k)-R2(1)
        x3 = rgrid(1,k)-R3(1)

        yo = rgrid(2,k)-RO(2)
        y1 = rgrid(2,k)-R1(2)
        y2 = rgrid(2,k)-R2(2)
        y3 = rgrid(2,k)-R3(2)

        zo = rgrid(3,k)-RO(3)
        z1 = rgrid(3,k)-R1(3)
        z2 = rgrid(3,k)-R2(3)
        z3 = rgrid(3,k)-R3(3)

        d  = dsqrt(xo*xo + yo*yo + zo*zo)
        d1 = dsqrt(x1*x1 + y1*y1 + z1*z1)
        d2 = dsqrt(x2*x2 + y2*y2 + z2*z2)
        d3 = dsqrt(x3*x3 + y3*y3 + z3*z3)

        xx = d**4
        rx = (Rpe)**4
        if (xx .gt. 1.0d-8) then
           depo = dv*rho(k)*0.5d0*alpha*
     >         ( -(2.0d0/sqrtpi)*(xx/rx)*dexp(-(xx/rx)**2)
     >          + util_erf(xx/rx))/(xx**2)    
           depo = depo*(4.0d0*d*d)
        else
           depo = 0.0d0
        end if
        

        if (d1.gt.1.0d-8) then
           dec1 = dv*rho(k)*q1*
     >          ( -(2.0d0/sqrtpi)*(d1/Rcoul)*dexp(-(d1/Rcoul)**2)
     >          + util_erf(d1/Rcoul))/(d1**3)
        else
           dec1 = 0.0d0
        end if
        if (d2.gt.1.0d-8) then
           dec2 = dv*rho(k)*q2*
     >          ( -(2.0d0/sqrtpi)*(d2/Rcoul)*dexp(-(d2/Rcoul)**2)
     >          + util_erf(d2/Rcoul))/(d2**3)
        else
           dec2 = 0.0d0
        end if
        if (d3.gt.1.0d-8) then
           dec3 = dv*rho(k)*q3*
     >          ( -(2.0d0/sqrtpi)*(d3/Rcoul)*dexp(-(d3/Rcoul)**2)
     >          + util_erf(d3/Rcoul))/(d3**3)
        else
          dec3 = 0.0d0
        end if

        fo(1) = fo(1) + d3do*(x3)*dec3 + (xo)*depo
        fo(2) = fo(2) + d3do*(y3)*dec3 + (yo)*depo
        fo(3) = fo(3) + d3do*(z3)*dec3 + (zo)*depo

        f1(1) = f1(1) + (x1)*dec1 + DELTA*(x3)*dec3
        f1(2) = f1(2) + (y1)*dec1 + DELTA*(y3)*dec3
        f1(3) = f1(3) + (z1)*dec1 + DELTA*(z3)*dec3

        f2(1) = f2(1) + (x2)*dec2 + DELTA*(x3)*dec3
        f2(2) = f2(2) + (y2)*dec2 + DELTA*(y3)*dec3
        f2(3) = f2(3) + (z2)*dec2 + DELTA*(z3)*dec3

      end do
      call D3dB_Vector_SumAll(3,fo)
      call D3dB_Vector_SumAll(3,f1)
      call D3dB_Vector_SumAll(3,f2)

      return 
      end






c $Id$
