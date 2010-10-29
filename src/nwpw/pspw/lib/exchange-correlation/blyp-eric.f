*
* $Id$
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


      
*     *****Becke exchange parameter*******
      real*8 beta
      parameter (beta = 0.0042d0)

*     *****local variables****************
      integer i
      real*8 n,fnx,fnx_lda,fdnx 
      real*8 xe
      real*8 agr 
      real*8 n_ponethird, n_pfourthirds,n_ptwothirds,n_pfivethirds 
      real*8 n3,n_peightthirds
      real*8 chi, chidn,chiddn 
      real*8 pi
      real*8 lda_c,c2
      real*8 H,F,Fdchi




      real*8 agr2
      real*8 fnc, fnc_1, fnc_2, fnc_3, fnc_4
      real*8 fdnc, fdnc_1, fdnc_2
      real*8 ce,ce_1,ce_2,ce_3
      real*8 a, b, c, d, Cf
      real*8 onethird
      real*8 c_onethird, d_onethird, c_twothirds,
     &       cc_onethird
      real*8 d_fourninths, dd_twoninths,
     &       cc_oneninth, c_fourteenninths
      real*8 ab_fourths, ab_seventyseconds
      real*8 n_onethird, n_fourthirds, n_fivethirds,
     &       n_fivethirds_2, n_seventhirds, n_eightthirds,
     &       n_ninethirds, n_eleventhirds,n_twelvethirds,
     &       n_thirteenthirds
      real*8 c_n_onethird, d_n_onethird, d_n_fourthirds
      real*8 exp_c_n_onethird
      real*8 F1, dF1dn, d2F1dn2,d2F1dn2_1, d2F1dn2_2, F12, F14
      real*8 G1, dG1dn,dG1dn_1, dG1dn_2, d2G1dn2,d2G1dn2_1, d2G1dn2_2,
     &       d2G1dn2_3, d2G1dn2_4,d2G1dn2_5,d2G1dn2_6

*******density cutoff parameters********************
      real*8 DNS_CUT, ETA
      parameter (DNS_CUT        =      1.0d-20)
      parameter (ETA            =      1.0d-20)
*******LYP correlation parameters a, b, c, d********
      parameter (a              =      0.04918d0)
      parameter (b              =      0.132d0)
      parameter (c              =      0.2533d0)
      parameter (d              =      0.349d0)
****************************************************

*******other local assignments************************
*******LDA parameter*******************************
      Cf             = (3.0d0/10.0d0)*(3.0d0*pi*pi)**2.0d0/3.0d0


      
      
  
     
      pi = 4.0d0*datan(1.0d0)
*     define lda constant***************************** 
      lda_c = (3.0d0/2.0d0)*(3.0d0/(4.0d0 * pi))**(1.0d0/3.0d0)
*     ***********************************************

      do i=1,n2ft3d
       if (rho_in(i).gt.ETA) then
       n        = rho_in(i) 
       agr      = agr_in(i)
******************************************************************
*     *******calc. becke exchange energy density, fnx, fdnx*******
*****************************************************************

       n_ponethird     = n**(1.0d0/3.0d0)
       n_ptwothirds    = n_ponethird*n_ponethird
       n_pfourthirds   = n_ponethird*n
       n_pfivethirds   = n_ponethird*n_pfourthirds
       n_peightthirds  = n_pfourthirds*n_pfourthirds
       n3 = n**3
       c2   = 2**(1.0d0/3.0d0)

*      **calculate chi, chidn, chiddn     
       chi            = agr/n_pfourthirds
       chidn          = -1.0*(4.0d0/3.0d0)*chi/n
       chiddn         = 1/n_pfourthirds

*      **calculate H
       H              = beta*c2*n_pfourthirds

*      **calculate F and dF/dchi
       F              = chi*chi/(1.0d0+6.0d0*beta*c2*chi*log(c2*chi 
     &                  + dsqrt(1.0d0+c2*c2*chi*chi)))
       Fdchi          = 2.0d0*F/chi - (F*F/(chi*chi))
     &                  *(6.0d0*beta*c2*log(c2*chi+sqrt(1.0d0 
     &                  +c2*c2*chi*chi))+6.0d0*beta*c2*c2*chi
     &                  /dsqrt(1.0d0 + c2*c2*chi*chi))

*      *calculate fnx and fdnx 
       fnx_lda        = -1.0d0*lda_c/c2 * (4.0d0/3.0d0) * n_ponethird
       fnx            = fnx_lda - (4.0d0/3.0d0)*H*F/n - H*Fdchi*chidn
       fdnx           = -1.0d0*H*Fdchi*chiddn

*      *calculate exchange energy density
       xe             = -1.0d0*(lda_c/c2)*n_ponethird 
     &                - (c2*beta*n_ponethird*chi*chi)
     &                /(1.0d0+6.0d0*beta*c2*chi*log(c2*chi 
     &                + dsqrt(1.0d0+c2*c2*chi*chi))) 
       
******calculate LYP correlation energy***************
       agr2 = agr*agr

       ce = (a*(agr2*b*(7*c*(d +n_ponethird) + 10*d*n_ponethird 
     >                    + 3*n_ptwothirds) - 
     >      72.0d0*(b*Cf + dexp(c/n_ponethird))
     >            *(d +n_ponethird)*n3))/
     >  (72.0d0*dexp(c/n_ponethird)*(d +n_ponethird)**2*n_peightthirds)

       fnc =  -(a*(72*(d + n_ponethird)*(b*Cf*(c*(d + n_ponethird) 
     >            + 4*d*n_ponethird + 
     >             3*n_ptwothirds) + dexp(c/n_ponethird)*(4*d 
     >             + 3*n_ponethird)*
     >           n_ponethird)*n3 + agr2*b*
     >        (-7*c**2*(d + n_ponethird)**2 + 40*d**2*n_ptwothirds 
     >            + 69*d*n + 15*n_pfourthirds + 
     >          c*(25*d**2*n_ponethird + 64*d*n_ptwothirds + 39*n))))/
     >  (216.0d0*Exp(c/n_ponethird)*(d + n_ponethird)**3*n3)

       fdnc =  (a*agr*b*(7*c*(d + n_ponethird) 
     >   + 10*d*n_ponethird + 3*n_ptwothirds))/
     >  (36.0d0*dexp(c/n_ponethird)*(d + n_ponethird)**2*n_pfivethirds)


*******final result for restricted LYP*****************************

c       write(*,*) "blyp:",i,n,agr,ce,fnc,fdnc



       xce(i)          = x_parameter*xe   + c_parameter*ce
       fn(i)           = x_parameter*fnx  + c_parameter*fnc
       fdn(i)          = x_parameter*fdnx + c_parameter*fdnc

       else
       xce(i) = 0.0d0
       fn(i)  = 0.0d0
       fdn(i) = 0.0d0
       end if

      end do
      return
      end
