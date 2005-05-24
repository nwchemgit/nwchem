*
* $Id: blyp.f,v 1.2 2005-05-24 17:36:27 bylaska Exp $
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
      real*8 n_ponethird, n_pfourthirds
      real*8 chi, chidn,chiddn
      real*8 pi
      real*8 lda_c,c2
      real*8 H,F,Fdchi




      real*8 agr2
      real*8 fnc, fnc_1, fnc_2, fnc_3, fnc_4
      real*8 fdnc, fdnc_1, fdnc_2
      real*8 ce,ce_1,ce_2,ce_3
      real*8 a, b, c, d, Cf
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
      parameter (DNS_CUT        =      1.0d-40)
      parameter (ETA            =      1.0d-70)
*******LYP correlation parameters a, b, c, d********
      parameter (a              =      0.04918d0)
      parameter (b              =      0.132d0)
      parameter (c              =      0.2533d0)
      parameter (d              =      0.349d0)
****************************************************

*******other local assignments************************
      c_onethird     =  c*1.0d0/3.0d0
      c_twothirds    =  c*2.0d0/3.0d0
      cc_onethird    =  c*c*1.0d0/3.0d0
      c_fourteenninths = c*14.0d0/9.0d0
      cc_oneninth    =  c*c*1.0d0/9.0d0
      ab_fourths     =  (a*b)/4.0d0
      ab_seventyseconds = (a*b)/72.0d0

      d_onethird     =  d*1.0d0/3.0d0
      d_fourninths   =  d*4.0d0/9.0d0
      dd_twoninths   =  d*d*2.0d0/9.0d0
*******LDA parameter*******************************
      Cf             = (3.0d0/10.0d0)*(3.0d0*pi*pi)**2.0d0/3.0d0






      pi = 4.0d0*datan(1.0d0)
*     define lda constant*****************************
      lda_c = (3.0d0/2.0d0)*(3.0d0/(4.0d0 * pi))**(1.0d0/3.0d0)
*     ***********************************************

      do i=1,n2ft3d
       n        = rho_in(i) + ETA
       agr      = agr_in(i) + ETA
******************************************************************
*     *******calc. becke exchange energy density, fnx, fdnx*******
*****************************************************************

       n_ponethird     = n**(1.0d0/3.0d0)
       n_pfourthirds   = n**(4.0d0/3.0d0)
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

       n_onethird     =  n**(-1.0d0/3.0d0)
       n_fourthirds   =  n**(-4.0d0/3.0d0)
       n_fivethirds   =  n**(-5.0d0/3.0d0)
       n_fivethirds_2 =  n**(+5.0d0/3.0d0)
       n_seventhirds  =  n**(-7.0d0/3.0d0)
       n_eightthirds  =  n**(-8.0d0/3.0d0)
       n_ninethirds   =  n**(-3.0d0)
       n_eleventhirds =  n**(-11.0d0/3.0d0)
       n_twelvethirds =  n**(-4.0d0)
       n_thirteenthirds = n**(-13.0d0/3.0d0)
       c_n_onethird   =  c*n_onethird
       d_n_onethird   = d*n_onethird
       d_n_fourthirds = d*n_fourthirds
       exp_c_n_onethird = exp(-c_n_onethird)


       F1             = 1.0d0/(1.0d0 + d_n_onethird)
       F12            = F1*F1
       F14            = F12*F12

       G1             = F1 * n_fivethirds * exp_c_n_onethird

       dF1dn          = 1.0d0/3.0d0*d_n_fourthirds*F12

       d2F1dn2_1      = (dd_twoninths*n_eightthirds/F1)
       d2F1dn2_2      = (d_fourninths*n_seventhirds/F12)
       d2F1dn2        = (d2F1dn2_1 - d2F1dn2_2)*F14

       dG1dn_1        = dF1dn*n_fivethirds
     &                  -F1*(5.0d0/3.0d0)*n_eightthirds
       dG1dn_2        = F1*c_onethird*n_ninethirds
       dG1dn          = (dG1dn_1 + dG1dn_2)*exp_c_n_onethird

       d2G1dn2_1      = d2F1dn2*n_fivethirds
       d2G1dn2_2      = c_twothirds*n_ninethirds
     &                  -(10.0d0/3.0d0)*n_eightthirds
       d2G1dn2_3      = (40.0d0/9.0d0)*n_eleventhirds
       d2G1dn2_4      = c_fourteenninths*n_twelvethirds
       d2G1dn2_5      = cc_oneninth*n_thirteenthirds
       d2G1dn2_6      = d2G1dn2_3-d2G1dn2_4 + d2G1dn2_5
       d2G1dn2        = (d2G1dn2_1+dF1dn*d2G1dn2_2+F1*d2G1dn2_6)*
     &                  exp_c_n_onethird


       fnc_1           = -a*(dF1dn*n+F1)
       fnc_2           =  a*b*Cf*n_fivethirds*(dG1dn*n+(8.0d0/3.0d0)*G1)
       fnc_3           =  ab_fourths*(d2G1dn2*n*agr2+3.0d0*dG1dn*agr2)
       fnc_4           =  ab_seventyseconds*(3.0d0*d2G1dn2*n*agr2+
     &                   5.0d0*dG1dn*agr2)

       fdnc_1          =  ab_fourths*(2.0d0*dG1dn*n*agr+4.0d0*G1*agr)
       fdnc_2          = ab_seventyseconds*(6.0d0*dG1dn*n*agr
     &                   +4.0d0*G1*agr)

       ce_1            =  -a*F1-a*b*Cf*n_fivethirds_2*G1
       ce_2            =   ab_fourths*(dG1dn*agr2+(2.0d0*G1*agr2)/n)
       ce_3            =   ab_seventyseconds*(3.0d0*dG1dn*agr2+
     &                     (2.0d0*G1*agr2)/n)

*******final result for restricted LYP*****************************
       fnc             =  fnc_1-fnc_2+fnc_3+fnc_4
       fdnc            =  fdnc_1+fdnc_2
       ce              = ce_1+ce_2+ce_3


       xce(i)          = x_parameter*xe   + c_parameter*ce
       fn(i)           = x_parameter*fnx  + c_parameter*fnc
       fdn(i)          = x_parameter*fdnx + c_parameter*fdnc

      end do
      return
      end
