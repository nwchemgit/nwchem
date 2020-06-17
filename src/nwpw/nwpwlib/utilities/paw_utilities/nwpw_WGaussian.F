c
c $Id$
c

*     ******************************************************
*     *                                                    *
*     *             nwpw_WGaussian                         *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral
*
*                                            //
*    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                            || ------------------------------------  dr dr'
*                                            //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
*            /
*
      real*8 function nwpw_WGaussian(la,ma,sa,lb,mb,sb,Rab)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,Rab(3)
 
*     *** local variables ***
      integer l,m,fac,fac2
      real*8 c,x,y,pi,tmp,mtmp,alpha
      real*8 cos_theta,phi,R

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm
      external nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm

      call nwpw_timing_start(49)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = dsqrt(0.25d0*(sa*sa + sb*sb))
      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R

      if ((dabs(Rab(2)).lt.1.0d-9).and.(dabs(Rab(1)).lt.1.0d-9)) then
         phi = 0.0d0
      else
         phi = datan2(Rab(2),Rab(1))
      end if

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      tmp = 0.0d0
      do l = abs(la-lb), (la+lb), 2
         mtmp = 0.0d0
         do m=-l,l
            mtmp = mtmp + nwpw_gaunt(.false.,l,m,la,ma,lb,mb)
     >                   *Tesseral_lm(l,m,cos_theta,phi)
         end do
         tmp = tmp + fac * mtmp * nwpw_GaussBessel(la+lb,l,alpha,R)
         fac = -fac
      end do
      call nwpw_timing_end(49)

      nwpw_WGaussian = c * tmp
      return
      end 


*     ******************************************************
*     *                                                    *
*     *             nwpw_dWGaussian                        *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral and it's derivative wrt to Rab
*
*                                            //
*   dWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                            || ------------------------------------  dr dr'
*                                            //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
*            /
*
      subroutine nwpw_dWGaussian(la,ma,sa,lb,mb,sb,Rab,W,dW)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,Rab(3)
      real*8  W,dW(3)
 
*     *** local variables ***
      integer l,m,fac
      real*8 c,x,y,pi,tmp,mtmp,mtmpx,mtmpy,mtmpz,alpha
      real*8 cos_theta,phi,R,gg1,gg2,gg3,Tx,Ty,Tz

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel,Tesseral_lm
      external nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel,Tesseral_lm

      call nwpw_timing_start(55)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = dsqrt(0.25d0*(sa*sa + sb*sb))

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R
      phi       = datan2(Rab(2),Rab(1))

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      W = 0.0d0
      dW(1) = 0.0d0
      dW(2) = 0.0d0
      dW(3) = 0.0d0
      do l = abs(la-lb), (la+lb), 2
         mtmp  = 0.0d0
         mtmpx = 0.0d0
         mtmpy = 0.0d0
         mtmpz = 0.0d0
         do m=-l,l
            gg1 = nwpw_gaunt(.false.,l,m,la,ma,lb,mb)
            call dTesseral_lm(l,m,cos_theta,phi,Tx,Ty,Tz)
            mtmp = mtmp   + gg1*Tesseral_lm(l,m,cos_theta,phi)
            mtmpx = mtmpx + gg1*Tx
            mtmpy = mtmpy + gg1*Ty
            mtmpz = mtmpz + gg1*Tz
         end do
         gg2 = nwpw_GaussBessel(la+lb,l,alpha,R)
         gg3 = nwpw_dGaussBessel(la+lb,l,alpha,R)
         W     = W     + fac*(mtmp*gg2)
         dW(1) = dW(1) + fac*(mtmpx*gg2 + mtmp*(Rab(1)/R)*gg3)
         dW(2) = dW(2) + fac*(mtmpy*gg2 + mtmp*(Rab(2)/R)*gg3)
         dW(3) = dW(3) + fac*(mtmpz*gg2 + mtmp*(Rab(3)/R)*gg3)
         fac = -fac
      end do
      W = W*c
      dW(1) = dW(1)*c
      dW(2) = dW(2)*c
      dW(3) = dW(3)*c
      call nwpw_timing_end(55)

      return
      end 

*     ******************************************************
*     *                                                    *
*     *             nwpw_UGaussian                         *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron one center Gaussian integral
*
*                                            //
*    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r) * g(lb,mb,sb;r')  
*                                            || ------------------------------------  dr dr'
*                                            //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
*            /
*
*
*   Note - this routine is equivalent to find_self_energy_coeff in the paw code
*
      real*8 function nwpw_UGaussian(la,ma,sa,lb,mb,sb)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb
 
*     *** local variables ***
      real*8 x,y,twopi,s,U

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial

      call nwpw_timing_start(48)
      U = 0.0d0
      if ((la.eq.lb).and.(ma.eq.mb)) then
         twopi = 8.0d0*datan(1.0d0)
         x = dble((2*la+1)*nwpw_doublefactorial(2*la+1))
         y = (dsqrt(0.5d0*(sa*sa+sb*sb)))**(2*la+1)
         U = 4.0d0*dsqrt(twopi)/(x*y) 
      end if
      call nwpw_timing_end(48)

      nwpw_UGaussian = U
      return
      end 


*************************** real combo versions ***********************************

*     ******************************************************
*     *                                                    *
*     *             nwpw_WGaussian3                        *
*     *                                                    *
*     ******************************************************
*
*     Calculates the 4 terms sum
*
*        WGaussian(la,ma,sa,lb,mb,sb,Rab)
*      - WGaussian(la,ma,sa,lb,mb,sm,Rab)
*      - WGaussian(la,ma,sm,lb,mb,sb,Rab)
*      + WGaussian(la,ma,sm,lb,mb,sm,Rab)
*
*    of the two electron two center Gaussian integral
*
*                                            //
*    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                            || ------------------------------------  dr dr'
*                                            //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
*            /
*
      real*8 function nwpw_WGaussian3(la,ma,sa,lb,mb,sb,sm,Rab)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,sm,Rab(3)
 
*     *** local variables ***
      integer l,m,fac,fac2
      real*8 c,x,y,pi,mtmp
      real*8 alpha_ab,alpha_am,alpha_mb,alpha_mm
      real*8 tmp
      real*8 cos_theta,phi,R,gg2

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm
      external nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm

      call nwpw_timing_start(49)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha_ab = dsqrt(0.25d0*(sa*sa + sb*sb))
      alpha_am = dsqrt(0.25d0*(sa*sa + sm*sm))
      alpha_mb = dsqrt(0.25d0*(sm*sm + sb*sb))
      alpha_mm = dsqrt(0.25d0*(sm*sm + sm*sm))

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R

      if ((dabs(Rab(2)).lt.1.0d-9).and.(dabs(Rab(1)).lt.1.0d-9)) then
         phi = 0.0d0
      else
         phi = datan2(Rab(2),Rab(1))
      end if

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      tmp = 0.0d0
      do l = abs(la-lb), (la+lb), 2
         mtmp = 0.0d0
         do m=-l,l
            mtmp = mtmp + nwpw_gaunt(.false.,l,m,la,ma,lb,mb)
     >                   *Tesseral_lm(l,m,cos_theta,phi)
         end do
         gg2 = nwpw_GaussBessel(la+lb,l,alpha_ab,R)
     >       - nwpw_GaussBessel(la+lb,l,alpha_am,R)
     >       - nwpw_GaussBessel(la+lb,l,alpha_mb,R)
     >       + nwpw_GaussBessel(la+lb,l,alpha_mm,R)

         tmp = tmp + fac*mtmp*gg2
         fac = -fac
      end do
      call nwpw_timing_end(49)

      nwpw_WGaussian3 = c*tmp
      return
      end 

*     ******************************************************
*     *                                                    *
*     *             nwpw_WGaussian2                        *
*     *                                                    *
*     ******************************************************
*
*     Calculates the 4 term sum
*
*        WGaussian(la,ma,sa,lb,mb,sa,Rab)
*      - WGaussian(la,ma,sa,lb,mb,sm,Rab)
*      - WGaussian(la,ma,sm,lb,mb,sa,Rab)
*      + WGaussian(la,ma,sm,lb,mb,sm,Rab)
*
*    of the two electron two center Gaussian integral
*
*                                            //
*    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                            || ------------------------------------  dr dr'
*                                            //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
*            /
*
      real*8 function nwpw_WGaussian2(la,ma,sa,lb,mb,sm,Rab)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,sm,Rab(3)
 
*     *** local variables ***
      integer l,m,fac,fac2
      real*8 c,x,y,pi,mtmp
      real*8 alpha_aa,alpha_am,alpha_mm
      real*8 tmp
      real*8 cos_theta,phi,R,gg2

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm
      external nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm

      call nwpw_timing_start(49)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha_aa = dsqrt(0.25d0*(sa*sa + sa*sa))
      alpha_am = dsqrt(0.25d0*(sa*sa + sm*sm))
      alpha_mm = dsqrt(0.25d0*(sm*sm + sm*sm))

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R

      if ((dabs(Rab(2)).lt.1.0d-9).and.(dabs(Rab(1)).lt.1.0d-9)) then
         phi = 0.0d0
      else
         phi = datan2(Rab(2),Rab(1))
      end if

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      tmp = 0.0d0
      do l = abs(la-lb), (la+lb), 2
         mtmp = 0.0d0
         do m=-l,l
            mtmp = mtmp + nwpw_gaunt(.false.,l,m,la,ma,lb,mb)
     >                   *Tesseral_lm(l,m,cos_theta,phi)
         end do

         gg2 = nwpw_GaussBessel(la+lb,l,alpha_aa,R)
     >       + nwpw_GaussBessel(la+lb,l,alpha_mm,R)
     >       - 2.0d0*nwpw_GaussBessel(la+lb,l,alpha_am,R)

         tmp = tmp + fac*mtmp*gg2
         
         fac = -fac
      end do
      call nwpw_timing_end(49)

      nwpw_WGaussian2 = c*tmp
      return
      end 

*     ******************************************************
*     *                                                    *
*     *             nwpw_dWGaussian3                       *
*     *                                                    *
*     ******************************************************
*
*     Calculates the 4 terms sum
*
*        WGaussian(la,ma,sa,lb,mb,sb,Rab)
*      - WGaussian(la,ma,sa,lb,mb,sm,Rab)
*      - WGaussian(la,ma,sm,lb,mb,sb,Rab)
*      + WGaussian(la,ma,sm,lb,mb,sm,Rab)
*
*     of the two electron two center Gaussian integral and it's derivative wrt to Rab
*
*                                            //
*   dWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                            || ------------------------------------  dr dr'
*                                            //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
*            /
*
      subroutine nwpw_dWGaussian3(la,ma,sa,lb,mb,sb,sm,Rab,W,dW)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,sm,Rab(3)
      real*8  W,dW(3)
 
*     *** local variables ***
      integer l,m,fac
      real*8 c,x,y,pi,tmp,mtmp,mtmpx,mtmpy,mtmpz,alpha
      real*8 cos_theta,phi,R,gg1,gg2,gg3,Tx,Ty,Tz
      real*8 W_ab,W_am,W_mb,W_mm
      real*8 dW_ab(3),dW_am(3),dW_mb(3),dW_mm(3)
      real*8 alpha_ab,alpha_am,alpha_mb,alpha_mm

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel,Tesseral_lm
      external nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel,Tesseral_lm

      call nwpw_timing_start(55)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha_ab = dsqrt(0.25d0*(sa*sa + sb*sb))
      alpha_am = dsqrt(0.25d0*(sa*sa + sm*sm))
      alpha_mb = dsqrt(0.25d0*(sm*sm + sb*sb))
      alpha_mm = dsqrt(0.25d0*(sm*sm + sm*sm))

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R
      phi       = datan2(Rab(2),Rab(1))

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      W = 0.0d0
      dW(1) = 0.0d0
      dW(2) = 0.0d0
      dW(3) = 0.0d0

      do l = abs(la-lb), (la+lb), 2
         mtmp  = 0.0d0
         mtmpx = 0.0d0
         mtmpy = 0.0d0
         mtmpz = 0.0d0
         do m=-l,l
            gg1 = nwpw_gaunt(.false.,l,m,la,ma,lb,mb)
            call dTesseral_lm(l,m,cos_theta,phi,Tx,Ty,Tz)
            mtmp = mtmp   + gg1*Tesseral_lm(l,m,cos_theta,phi)
            mtmpx = mtmpx + gg1*Tx
            mtmpy = mtmpy + gg1*Ty
            mtmpz = mtmpz + gg1*Tz
         end do

         gg2 = nwpw_GaussBessel(la+lb,l,alpha_ab,R)
     >       - nwpw_GaussBessel(la+lb,l,alpha_am,R)
     >       - nwpw_GaussBessel(la+lb,l,alpha_mb,R)
     >       + nwpw_GaussBessel(la+lb,l,alpha_mm,R)

         gg3 = nwpw_dGaussBessel(la+lb,l,alpha_ab,R)
     >       - nwpw_dGaussBessel(la+lb,l,alpha_am,R)
     >       - nwpw_dGaussBessel(la+lb,l,alpha_mb,R)
     >       + nwpw_dGaussBessel(la+lb,l,alpha_mm,R)

         W     = W     + fac*(mtmp*gg2)
         dW(1) = dW(1) + fac*(mtmpx*gg2 + mtmp*(Rab(1)/R)*gg3)
         dW(2) = dW(2) + fac*(mtmpy*gg2 + mtmp*(Rab(2)/R)*gg3)
         dW(3) = dW(3) + fac*(mtmpz*gg2 + mtmp*(Rab(3)/R)*gg3)

         fac = -fac
      end do
      W = W*c
      dW(1) = dW(1)*c
      dW(2) = dW(2)*c
      dW(3) = dW(3)*c
      call nwpw_timing_end(55)

      return
      end 


*************************** real combo versions ***********************************



*************************** complex versions ***********************************


*     ******************************************************
*     *                                                    *
*     *             nwpw_CWGaussian                        *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral
*
*                                             //
*    CWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                             || ------------------------------------  dr dr'
*                                             //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *conjg(Ylm(rhat)) dr = 1
*            /
*
      complex*16 function nwpw_CWGaussian(la,ma,sa,lb,mb,sb,Rab)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,Rab(3)
 
*     *** local variables ***
      integer l,m,fac,fac2
      real*8 c,x,y,pi,alpha
      real*8 cos_theta,phi,R
      complex*16 mtmp,tmp

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel
      external nwpw_gaunt,nwpw_GaussBessel
      complex*16 Yspherical_lm
      external   Yspherical_lm
      real*8   nwpw_gaunt_sub,gen_gaunt_coeff_sub
      external nwpw_gaunt_sub,gen_gaunt_coeff_sub


      call nwpw_timing_start(49)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = dsqrt(0.25d0*(sa*sa + sb*sb))
      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R

      if ((dabs(Rab(2)).lt.1.0d-9).and.(dabs(Rab(1)).lt.1.0d-9)) then
         phi = 0.0d0
      else
         phi = datan2(Rab(2),Rab(1))
      end if

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      !phase_factor = (-1)**(m1+l1)/sigma**(l1+l2+1)
      !if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
      !if (mod(lb+mb,2).eq.1) then 

      
      !if (mod(lb+mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
      if (mod(mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      m = ma + mb
      tmp = 0.0d0
      do l = abs(la-lb), (la+lb), 2
         if (abs(m).le.l) then
c            write(*,*) "l,m=",l,m,la,ma,lb,mb,
c     >                 nwpw_gaunt(.true.,l,m,la,ma,lb,-mb),
c     >                 nwpw_gaunt_sub(.true.,l,m,la,ma,lb,-mb),
c     >                 gen_gaunt_coeff_sub(l,m,la,ma,lb,-mb)
            mtmp = nwpw_gaunt(.true.,l,m,la,ma,lb,-mb)
     >            *YSpherical_lm(l,m,cos_theta,phi)
            tmp = tmp + fac * mtmp * nwpw_GaussBessel(la+lb,l,alpha,R)
         end if
         fac = -fac
      end do
      call nwpw_timing_end(49)

      nwpw_CWGaussian = c * tmp
      return
      end 


*     ******************************************************
*     *                                                    *
*     *             nwpw_dCWGaussian                       *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral and it's derivative wrt to Rab
*
*                                             //
*   dCWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                             || ------------------------------------  dr dr'
*                                             //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l * congj(Ylm(rhat)) dr = 1
*            /
*
      subroutine nwpw_dCWGaussian(la,ma,sa,lb,mb,sb,Rab,W,dW)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,Rab(3)
      complex*16  W,dW(3)
 
*     *** local variables ***
      integer l,m,fac
      real*8 c,x,y,pi,alpha
      real*8 cos_theta,phi,R,gg1,gg2,gg3
      complex*16 mtmp,mtmpx,mtmpy,mtmpz,Tx,Ty,Tz

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel
      external nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel
      complex*16 Yspherical_lm
      external   Yspherical_lm

      call nwpw_timing_start(55)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = dsqrt(0.25d0*(sa*sa + sb*sb))

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R
      phi       = datan2(Rab(2),Rab(1))

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      !if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
      !if (mod(lb+mb,2).eq.1) then 

      if (mod(mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      W = dcmplx(0.0d0,0.0d0)
      dW(1) = dcmplx(0.0d0,0.0d0)
      dW(2) = dcmplx(0.0d0,0.0d0)
      dW(3) = dcmplx(0.0d0,0.0d0)
      m = ma + mb
      do l = abs(la-lb), (la+lb), 2
         if (abs(m).le.l) then
            gg1 = nwpw_gaunt(.true.,l,m,la,ma,lb,-mb)
            call dYspherical_lm(l,m,cos_theta,phi,Tx,Ty,Tz)
            mtmp  = gg1*Yspherical_lm(l,m,cos_theta,phi)
            mtmpx = gg1*Tx
            mtmpy = gg1*Ty
            mtmpz = gg1*Tz
            gg2 = nwpw_GaussBessel(la+lb,l,alpha,R)
            gg3 = nwpw_dGaussBessel(la+lb,l,alpha,R)
            W     = W     + fac*(mtmp*gg2)
            dW(1) = dW(1) + fac*(mtmpx*gg2 + mtmp*(Rab(1)/R)*gg3)
            dW(2) = dW(2) + fac*(mtmpy*gg2 + mtmp*(Rab(2)/R)*gg3)
            dW(3) = dW(3) + fac*(mtmpz*gg2 + mtmp*(Rab(3)/R)*gg3)
         end if
         fac = -fac
      end do
      W = W*c
      dW(1) = dW(1)*c
      dW(2) = dW(2)*c
      dW(3) = dW(3)*c
      call nwpw_timing_end(55)

      return
      end 


*************************** complex versions ***********************************

*************************** complex combo versions ***********************************

*     ******************************************************
*     *                                                    *
*     *             nwpw_CWGaussian3                       *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral
*
*                                             //
*    CWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                             || ------------------------------------  dr dr'
*                                             //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *conjg(Ylm(rhat)) dr = 1
*            /
*
      complex*16 function nwpw_CWGaussian3(la,ma,sa,lb,mb,sb,sm,Rab)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,sm,Rab(3)
 
*     *** local variables ***
      integer l,m,fac,fac2
      real*8 c,x,y,pi
      real*8 alpha_ab,alpha_am,alpha_mb,alpha_mm
      real*8 cos_theta,phi,R,gg2
      complex*16 mtmp,tmp

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel
      external nwpw_gaunt,nwpw_GaussBessel
      complex*16 Yspherical_lm
      external   Yspherical_lm
      real*8   nwpw_gaunt_sub,gen_gaunt_coeff_sub
      external nwpw_gaunt_sub,gen_gaunt_coeff_sub


      call nwpw_timing_start(49)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha_ab = dsqrt(0.25d0*(sa*sa + sb*sb))
      alpha_am = dsqrt(0.25d0*(sa*sa + sm*sm))
      alpha_mb = dsqrt(0.25d0*(sm*sm + sb*sb))
      alpha_mm = dsqrt(0.25d0*(sm*sm + sm*sm))
      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R

      if ((dabs(Rab(2)).lt.1.0d-9).and.(dabs(Rab(1)).lt.1.0d-9)) then
         phi = 0.0d0
      else
         phi = datan2(Rab(2),Rab(1))
      end if

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      !phase_factor = (-1)**(m1+l1)/sigma**(l1+l2+1)
      !if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
      !if (mod(lb+mb,2).eq.1) then 

      
      !if (mod(lb+mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
      if (mod(mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      m = ma + mb
      tmp = dcmplx(0.0d0,0.0d0)
      do l = abs(la-lb), (la+lb), 2
         if (abs(m).le.l) then
c            write(*,*) "l,m=",l,m,la,ma,lb,mb,
c     >                 nwpw_gaunt(.true.,l,m,la,ma,lb,-mb),
c     >                 nwpw_gaunt_sub(.true.,l,m,la,ma,lb,-mb),
c     >                 gen_gaunt_coeff_sub(l,m,la,ma,lb,-mb)
            mtmp = nwpw_gaunt(.true.,l,m,la,ma,lb,-mb)
     >            *YSpherical_lm(l,m,cos_theta,phi)

            gg2 = nwpw_GaussBessel(la+lb,l,alpha_ab,R)
     >          - nwpw_GaussBessel(la+lb,l,alpha_am,R)
     >          - nwpw_GaussBessel(la+lb,l,alpha_mb,R)
     >          + nwpw_GaussBessel(la+lb,l,alpha_mm,R)

            tmp = tmp + fac*mtmp*gg2

         end if
         fac = -fac
      end do
      call nwpw_timing_end(49)

      nwpw_CWGaussian3 = c*tmp
      return
      end 


*     ******************************************************
*     *                                                    *
*     *             nwpw_CWGaussian2                       *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral
*
*                                             //
*    CWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                             || ------------------------------------  dr dr'
*                                             //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l *conjg(Ylm(rhat)) dr = 1
*            /
*
      complex*16 function nwpw_CWGaussian2(la,ma,sa,lb,mb,sm,Rab)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sm,Rab(3)
 
*     *** local variables ***
      integer l,m,fac,fac2
      real*8 c,x,y,pi
      real*8 alpha_aa,alpha_am,alpha_mm
      real*8 cos_theta,phi,R,gg2
      complex*16 mtmp,tmp

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel
      external nwpw_gaunt,nwpw_GaussBessel
      complex*16 Yspherical_lm
      external   Yspherical_lm
      real*8   nwpw_gaunt_sub,gen_gaunt_coeff_sub
      external nwpw_gaunt_sub,gen_gaunt_coeff_sub


      call nwpw_timing_start(49)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha_aa = dsqrt(0.25d0*(sa*sa + sa*sa))
      alpha_am = dsqrt(0.25d0*(sa*sa + sm*sm))
      alpha_mm = dsqrt(0.25d0*(sm*sm + sm*sm))
      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R

      if ((dabs(Rab(2)).lt.1.0d-9).and.(dabs(Rab(1)).lt.1.0d-9)) then
         phi = 0.0d0
      else
         phi = datan2(Rab(2),Rab(1))
      end if

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      !phase_factor = (-1)**(m1+l1)/sigma**(l1+l2+1)
      !if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
      !if (mod(lb+mb,2).eq.1) then 

      
      !if (mod(lb+mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
      if (mod(mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      m = ma + mb
      tmp = dcmplx(0.0d0,0.0d0)
      do l = abs(la-lb), (la+lb), 2
         if (abs(m).le.l) then
c            write(*,*) "l,m=",l,m,la,ma,lb,mb,
c     >                 nwpw_gaunt(.true.,l,m,la,ma,lb,-mb),
c     >                 nwpw_gaunt_sub(.true.,l,m,la,ma,lb,-mb),
c     >                 gen_gaunt_coeff_sub(l,m,la,ma,lb,-mb)
            mtmp = nwpw_gaunt(.true.,l,m,la,ma,lb,-mb)
     >            *YSpherical_lm(l,m,cos_theta,phi)

            gg2 = nwpw_GaussBessel(la+lb,l,alpha_aa,R)
     >          + nwpw_GaussBessel(la+lb,l,alpha_mm,R)
     >          - 2.0d0*nwpw_GaussBessel(la+lb,l,alpha_am,R)

            tmp = tmp + fac*mtmp*gg2

         end if
         fac = -fac
      end do
      call nwpw_timing_end(49)

      nwpw_CWGaussian2 = c*tmp
      return
      end 



*     ******************************************************
*     *                                                    *
*     *             nwpw_dCWGaussian3                      *
*     *                                                    *
*     ******************************************************
*
*     Calculates the two electron two center Gaussian integral and it's derivative wrt to Rab
*
*                                             //
*   dCWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
*                                             || ------------------------------------  dr dr'
*                                             //                |r-r'|
*
*     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 
*
*          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
*           
*     The normalization constant C_l is defined such at
*            /
*            | g(l,m,s;r) * |r|**l * congj(Ylm(rhat)) dr = 1
*            /
*
      subroutine nwpw_dCWGaussian3(la,ma,sa,lb,mb,sb,sm,Rab,W,dW)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,sm,Rab(3)
      complex*16  W,dW(3)
 
*     *** local variables ***
      integer l,m,fac
      real*8 c,x,y,pi
      real*8 alpha_ab,alpha_am,alpha_mb,alpha_mm
      real*8 cos_theta,phi,R,gg1,gg2,gg3
      complex*16 mtmp,mtmpx,mtmpy,mtmpz,Tx,Ty,Tz

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel
      external nwpw_gaunt,nwpw_GaussBessel,nwpw_dGaussBessel
      complex*16 Yspherical_lm
      external   Yspherical_lm

      call nwpw_timing_start(55)
      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha_ab = dsqrt(0.25d0*(sa*sa + sb*sb))
      alpha_am = dsqrt(0.25d0*(sa*sa + sm*sm))
      alpha_mb = dsqrt(0.25d0*(sm*sm + sb*sb))
      alpha_mm = dsqrt(0.25d0*(sm*sm + sm*sm))

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R
      phi       = datan2(Rab(2),Rab(1))

      if (mod(2*la+lb,2).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      !if (mod((abs(la-lb)+la+lb)/2,2).eq.1) then 
      !if (mod(lb+mb,2).eq.1) then 

      if (mod(mb+(abs(la-lb)+la+lb)/2,2).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      W = dcmplx(0.0d0,0.0d0)
      dW(1) = dcmplx(0.0d0,0.0d0)
      dW(2) = dcmplx(0.0d0,0.0d0)
      dW(3) = dcmplx(0.0d0,0.0d0)

      m = ma + mb
      do l = abs(la-lb), (la+lb), 2
         if (abs(m).le.l) then
            gg1 = nwpw_gaunt(.true.,l,m,la,ma,lb,-mb)
            call dYspherical_lm(l,m,cos_theta,phi,Tx,Ty,Tz)
            mtmp  = gg1*Yspherical_lm(l,m,cos_theta,phi)
            mtmpx = gg1*Tx
            mtmpy = gg1*Ty
            mtmpz = gg1*Tz

            gg2 = nwpw_GaussBessel(la+lb,l,alpha_ab,R)
     >          - nwpw_GaussBessel(la+lb,l,alpha_am,R)
     >          - nwpw_GaussBessel(la+lb,l,alpha_mb,R)
     >          + nwpw_GaussBessel(la+lb,l,alpha_mm,R)

            gg3 = nwpw_dGaussBessel(la+lb,l,alpha_ab,R)
     >          - nwpw_dGaussBessel(la+lb,l,alpha_am,R)
     >          - nwpw_dGaussBessel(la+lb,l,alpha_mb,R)
     >          + nwpw_dGaussBessel(la+lb,l,alpha_mm,R)

            W     = W     + fac*(mtmp*gg2)
            dW(1) = dW(1) + fac*(mtmpx*gg2 + mtmp*(Rab(1)/R)*gg3)
            dW(2) = dW(2) + fac*(mtmpy*gg2 + mtmp*(Rab(2)/R)*gg3)
            dW(3) = dW(3) + fac*(mtmpz*gg2 + mtmp*(Rab(3)/R)*gg3)
         end if
         fac = -fac
      end do
      W = W*c
      dW(1) = dW(1)*c
      dW(2) = dW(2)*c
      dW(3) = dW(3)*c
      call nwpw_timing_end(55)

      return
      end 


*************************** complex combo versions ***********************************
