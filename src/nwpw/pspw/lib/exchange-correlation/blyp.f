*
* $Id: blyp.f,v 1.3 2005-09-16 01:15:49 bylaska Exp $
*


*     blyp restricted  calc.    
*
*
*

*      subroutine gen_BLYP_BW_restricted(n2ft3d,rho_in,agr_in,xce,xcp,fn,fdn)
*      input:  n2ft3d                  grid
*              rho_in                  density
*              agr_in                  absolute gradient of density
*              x_parameter:            scale parameter for exchange
*              c_parameter:            scale parameter for correlation
*      output: xce                     exchange correlation energy density
*              fn                      d(n*exc)/dn
*              fdn                     d(n*exc)/d(|grad n|)

      subroutine gen_BLYP_BW_restricted(n2ft3d,rho_in,agr_in,
     &                               x_parameter,
     &                               c_parameter,xce,fn,fdn)
      implicit none

      integer   n2ft3d
      real*8    rho_in(n2ft3d)
      real*8    agr_in(n2ft3d)
      real*8    xce(n2ft3d)
      real*8    fn(n2ft3d)
      real*8    fdn(n2ft3d)
      real*8    x_parameter, c_parameter


*******local declarations***************************************
      integer   i
      real*8    pi
      real*8    thrd
      real*8    two_thrd
      real*8    beta, lda_c,Cf
      real*8    Fc,Gc,C1,C2
      real*8    n, n_thrd,n_fv,n_fr,n_tw
      real*8    agr,agr2, chi, chi2,chiSQ,sd
      real*8    K
      real*8    F1, F2      
      real*8    xe,fdnx,fdagrx
      real*8    ce,fdnc, fdagrc,dfdnc_lda
      

*******density cutoff parameters********************
      real*8 DNS_CUT, ETA
      parameter (DNS_CUT        =      1.0d-18)
      parameter (ETA            =      1.0d-20)
*******LYP correlation parameters a, b, c, d********
      real*8 a,b,c,d
      parameter (a              =      0.04918d0)
      parameter (b              =      0.132d0)
      parameter (c              =      0.2533d0)
      parameter (d              =      0.349d0)
****************************************************
      pi 	= 4.0d0*datan(1.0d0)
      thrd 	= dble(1.0d0/3.0d0)
      two_thrd 	= dble(2.0d0**thrd)
      beta 	= dble(0.0042d0)
      lda_c 	= dble(3.0d0/(4.0d0*pi))
      lda_c 	= dble(lda_c**thrd)
*      lda_c 	= dble(lda_c*3.0d0/4.0d0)
      lda_c 	= dble(lda_c*3.0d0/2.0d0)


      do i=1,n2ft3d
       n        = rho_in(i) + ETA
       agr      = agr_in(i)
       n_thrd 	= n**thrd
       if (rho_in(i).lt.DNS_CUT) then
         xe     = 0.0d0
         fdnx   = 0.0d0
         fdagrx = 0.0d0
       else
******************************************************************
*     *******calc. becke exchange energy density, fnx, fdnx*******
*****************************************************************
       sd       = 1.0d0/(n_thrd*n)
       chi 	= two_thrd*agr*sd
       chi2	= chi*chi
       chiSQ    = dsqrt(1.0d0+chi2)   



       K 	= 6.0d0*beta*dlog(chi+chiSQ)
       F1 	= chi2/(1.0d0+chi*K)
       xe 	= -n_thrd*(lda_c+beta*F1)/two_thrd
       F2 	= (2.0d0 + chi*K-(chi2)*6.0d0*beta
     &		/chiSQ)
     &      	/((1.0d0+chi*K)*(1.0d0+chi*K))
       fdnx 	= -(n_thrd/two_thrd)*dble(4.0d0/3.0d0)
     &  	*(lda_c+beta*(F1-chi2*F2))
       fdagrx 	= -beta*chi*F2 
       end if


*******final result for restricted LYP*****************************
       agr2 = agr*agr
       n_fv = n_thrd*n_thrd*n_thrd*n_thrd*n_thrd
       n_fr = n_thrd*n_thrd*n_thrd*n_thrd
       n_tw = n_thrd*n_thrd
       Cf = dble(3.0d0*pi*pi)
       Cf = dble(3.0d0*(Cf**(2.0d0/3.0d0))/10.0d0)
       Fc = dble(1.0d0/(1.0d0+d/n_thrd))
       Gc = Fc*dexp(-c/n_thrd)/n_fv
       C1 = dble(1.0d0+7.0d0*(c+d*Fc)/(n_thrd*3.0d0))
       C2 = (d*Fc-5.0d0*n_thrd/3.0d0+c/3.0d0)/n_fr
       dfdnc_lda =  -a*(Fc*(d*Fc/(n_thrd*3.0d0) + 1.0d0)
     &      +b*Cf*8.0d0*n_fv/3.0d0)


       ce = -a*Fc - a*b*Cf*n_fv + a*b*agr2*Gc*C1/(n*24.0d0)
       fdnc = dfdnc_lda+dble(a*b/24.0d0)*agr2*Gc
     &        *(C2*C1-(Fc*d/(3.0d0*n_fr))*(7.0d0*c/(3.0d0*Fc*d)
     &        +1.0d0-(d*Fc/n_thrd)))
       fdagrc = a*b*agr*Gc*C1/12.0d0



       xce(i) 	= x_parameter*xe   + c_parameter*ce
       fn(i)	= x_parameter*fdnx  + c_parameter*fdnc
       fdn(i) 	= x_parameter*fdagrx + c_parameter*fdagrc

      end do
      return
      end
