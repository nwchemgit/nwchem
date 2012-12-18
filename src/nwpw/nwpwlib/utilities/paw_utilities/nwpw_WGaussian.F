c
c $Id: nwpw_WGaussian.F 21176 2011-10-10 06:35:49Z d3y133 $
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
      integer l,m,fac
      real*8 c,x,y,pi,tmp,mtmp,alpha
      real*8 cos_theta,phi,R

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm
      external nwpw_gaunt,nwpw_GaussBessel,Tesseral_lm

      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = 0.25d0*(sa*sa + sb*sb)
      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R
      phi       = datan2(Rab(2),Rab(1))

      if (mod(2*la+lb,1).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,1).eq.1) then 
         fac = -1
      else
         fac = 1
      end if
      
      tmp = 0.0d0
      do l = abs(la-lb), (la+lb), 2
         mtmp = 0.0d0
         do m=-l,l
            mtmp = mtmp + nwpw_gaunt(.false.,l,m,la,ma,lb,mb)
     <                   *Tesseral_lm(l,m,cos_theta,phi)
         end do
         tmp = tmp + fac * mtmp * nwpw_GaussBessel(la+lb,l,alpha,R)
         fac = -fac
      end do

      nwpw_WGaussian = c * tmp
      return
      end 


*     ******************************************************
*     *                                                    *
*     *             nwpw_dWGaussian                        *
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

      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = 0.25d0*(sa*sa + sb*sb)

      R = dsqrt(Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3))
      cos_theta = Rab(3)/R
      phi       = datan2(Rab(2),Rab(1))

      if (mod(2*la+lb,1).eq.1) then 
         c = -32.0d0*pi/(x*y)
      else
         c = 32.0d0*pi/(x*y)
      end if

      if (mod((abs(la-lb)+la+lb)/2,1).eq.1) then 
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
      real*8 function nwpw_UGaussian(la,ma,sa,lb,mb,sb)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb
 
*     *** local variables ***
      real*8 x,y,twopi,s,U

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial

      U = 0.0d0
      if ((la.eq.lb).and.(ma.eq.mb)) then
         twopi = 8.0d0*datan(1.0d0)
         x = dble((2*la+1)*nwpw_doublefactorial(2*la+1))
         y = (dsqrt(0.5d0*(sa*sa+sb*sb)))**(2*la+1)
         U = 4.0d0*dsqrt(twopi)/(x*y) 
      end if

      nwpw_UGaussian = U
      return
      end 

