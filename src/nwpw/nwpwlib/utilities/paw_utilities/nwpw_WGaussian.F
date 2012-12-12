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
      real*8 function nwpw_WGaussian(la,ma,sa,Ra, lb,mb,sb,Rb)
      implicit none
      integer la,ma,lb,mb
      real*8  sa,sb,Ra(3),Rb(3)
 
*     *** local variables ***
      integer l,m,fac
      real*8 c,x,y,pi,tmp,mtmp,alpha
      real*8 Rx,Ry,Rz,R

*     **** external functions ****
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   nwpw_gaunt,nwpw_GaussBessel,Tesseral3_lm
      external nwpw_gaunt,nwpw_GaussBessel,Tesseral3_lm

      pi = 4.0d0*datan(1.0d0)
      x = dble(nwpw_doublefactorial(2*la+1))
      y = dble(nwpw_doublefactorial(2*lb+1))
      alpha = 0.25d0*(sa*sa + sb*sb)
      Rx = Ra(1) - Rb(1)
      Ry = Ra(2) - Rb(2)
      Rz = Ra(3) - Rb(3)
      R = dsqrt(Rx*Rx + Ry*Ry + Rz*Rz)

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
            mtmp = mtmp + nwpw_gaunt(l,m,la,ma,lb,mb)
     <                   *Tesseral3_lm(l,m,Rx,Ry,Rz)
         end do
         tmp = tmp + fac * mtmp * nwpw_GaussBessel(la+lb,l,alpha,R)
         fac = -fac
      end do

      nwpw_WGaussian = c * tmp
      return
      end 
