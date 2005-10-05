*
* $Id: blyp.f,v 1.4 2005-10-05 21:19:08 bylaska Exp $
*

*     **************************************************
*     *                                                *
*     *           gen_BLYP_BW_unrestricted             *
*     *                                                *
*     **************************************************

      subroutine gen_BLYP_BW_unrestricted(n2ft3d,
     &                          dn_in,agr_in,
     &                          x_parameter,c_parameter,
     &                          xce,fn,fdn)

      implicit none

      integer   n2ft3d
      real*8    dn_in(n2ft3d,2)
      real*8    agr_in(n2ft3d,3)
      real*8    xce(n2ft3d)
      real*8    fn(n2ft3d,2)
      real*8    fdn(n2ft3d,3)
      real*8    x_parameter, c_parameter


      real*8 ETA,ETA2,DNS_CUT
      parameter (ETA = 0.0d0)
      parameter (ETA2 = 1.0d-20)
      parameter (DNS_CUT        =      1.0d-18)

*     ***LYP parameters****************
      real*8 a,b,c,d
      parameter (a = 0.04918d0)
      parameter (b = 0.132d0)
      parameter (c = 0.2533d0)
      parameter (d = 0.349d0)
*     ***Becke88 parameters*************
      real*8 beta
      parameter (beta = 0.0042d0)
      
      real*8 thrd
      parameter (thrd = 1.0d0/3.0d0)

      integer j
      real*8 twthrd,frthrd,fvthrd,snthrd,etthrd
      real*8 pi, Cf,lda_c

      real*8 nup,ndn,n,agrup,agrdn,agr
      real*8 agr2, agrup2,agrdn2
      real*8 n2,nup2,ndn2
      real*8 gamma,gamma_u,gamma_uu
      real*8 gamma_d,gamma_dd,gamma_du
      real*8 F,F_u,F_uu
      real*8 F_d,F_dd,F_du
      real*8 n13d
      real*8 enthrd
      real*8 G,G_u,G_uu
      real*8 G_d,G_dd,G_du
      real*8 fc,fclda,H0,Hu,Hd
      real*8 ce
      real*8 Hu_u,Hu_d,Hd_d,Hd_u
      real*8 fc_u,fc_d
      real*8 H0_u,H0_d
      real*8 fclda_u,fclda_d
      real*8 fc_agr, fc_agrup,fc_agrdn

      real*8 chiup,chidn
      real*8 chiup2,chidn2,chiupSQ,chidnSQ
      real*8 Kup,Kdn,F1up,F1dn,F2up,F2dn
      real*8 xeup,xedn,xe
      real*8 fdnxup,fdnxdn,fdagrxup,fdagrxdn
      real*8  Q,Q1,P,P1,n_m,n2_m,n3_m,n4_m
      real*8  n13d_m,n13d2_m,n13d3_m,d3
      real*8  n_thrd,nupndn
      real*8  n_mthrd,n_mfrthrd,n_mfvthrd
      real*8   nup_etthrd,nup_fvthrd
      real*8   ndn_etthrd,ndn_fvthrd
 
      
      twthrd = thrd*2.0d0
      frthrd = thrd*4.0d0
      fvthrd = thrd*5.0d0
      snthrd = thrd*7.0d0
      etthrd = thrd*8.0d0

      pi = 4.0d0*datan(1.0d0)
      Cf = dble(3.0d0*pi*pi)
      Cf = dble(3.0d0*(Cf**(2.0d0/3.0d0))/10.0d0)

      lda_c =(3.0d0/2.0d0)*(3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)


      do j=1,n2ft3d
       nup     = dn_in(j,1) + 0.5d0*ETA2 
       ndn     = dn_in(j,2) + 0.5d0*ETA2 

       agrup   = agr_in(j,1)
       agrdn   = agr_in(j,2)
       agrup2  = agrup*agrup
       agrdn2  = agrdn*agrdn

       n       = nup + ndn
       agr     = agr_in(j,3)
       agr2    = agr*agr

       n2      = n*n
       nup2    = nup*nup
       ndn2    = ndn*ndn

       
       if ((dn_in(j,1)+dn_in(j,2)).lt.DNS_CUT) then
          xe       = 0.0d0
          fdnxup   = 0.0d0
          fdnxdn   = 0.0d0
          fdagrxup = 0.0d0
          fdagrxdn = 0.0d0
       else
      
*      *******exchange part***************


*      **************UP*******************
       chiup = agrup/nup**(4.0d0/3.0d0)
       chiup2 = chiup*chiup
       chiupSQ = dsqrt(1.0d0+chiup2)

       Kup = 6.0d0*beta*dlog(chiup+chiupSQ)
       F1up = chiup2/(1.0d0 + chiup*Kup)
       xeup = -nup**thrd*(lda_c + beta*F1up)
       F2up = (2.0d0 + chiup*Kup 
     &        - 6.0d0*beta*chiup2/chiupSQ)
     &        /(1.0d0+chiup*Kup)**2.0d0
       fdnxup = -nup**(thrd)*(4.0d0/3.0d0)
     &          *(lda_c+beta*(F1up-chiup2*F2up))
       fdagrxup = -beta*chiup*F2up 
*      ************END UP*****************


*      *************DOWN******************
       chidn = agrdn/ndn**(4.0d0/3.0d0)
       chidn2 = chidn*chidn
       chidnSQ = dsqrt(1.0d0+chidn2)

       Kdn = 6.0d0*beta*dlog(chidn+chidnSQ)
       F1dn = chidn2/(1.0d0 + chidn*Kdn)
       xedn = -ndn**thrd*(lda_c + beta*F1dn)
       F2dn = (2.0d0 + chidn*Kdn
     &        - 6.0d0*beta*chidn2/chidnSQ)
     &        /(1.0d0+chidn*Kdn)**2.0d0
       fdnxdn = -ndn**(thrd)*(4.0d0/3.0d0)
     &          *(lda_c+beta*(F1dn-chidn2*F2dn))
       fdagrxdn = -beta*chidn*F2dn


*      ***********END DOWN****************

       xe = (xeup*nup + xedn*ndn)/n

*      *******end excange part************
       end if


*      *******correlation part************
       n_thrd   = n**thrd
       n_m     = 1.0d0/n
       n2_m    = n_m*n_m
       n3_m    = n2_m*n_m
       n4_m    = n3_m*n_m
       nupndn  = nup*ndn
       n_mthrd = 1.0d0/n_thrd
       n_mfrthrd = n_mthrd*n_m
       n_mfvthrd = n_mfrthrd*n_mthrd

       nup_etthrd = nup**etthrd
       nup_fvthrd = nup**fvthrd
       ndn_etthrd = ndn**etthrd
       ndn_fvthrd = ndn**fvthrd




       gamma   = (4.0d0*nupndn)*n2_m
       gamma_u = 4.0d0*ndn*n2_m - 8.0d0*nupndn*n3_m
       gamma_d = 4.0d0*nup*n2_m - 8.0d0*nupndn*n3_m

       gamma_uu = -16.0d0*ndn*n3_m + 24.0d0*nupndn*n4_m
       gamma_dd = -16.0d0*nup*n3_m + 24.0d0*nupndn*n4_m

       gamma_du =  (6.0d0*gamma-4.0d0)*n2_m



       d3         = d/3.0d0
       n13d_m     = 1.0d0/(1.0d0 + d*n_mthrd)
       n13d2_m    = n13d_m*n13d_m
       n13d3_m    = n13d2_m*n13d_m
       F        = gamma*n13d_m
       F_u      = gamma_u*n13d_m + d3*gamma*n_mfrthrd*n13d2_m
       F_d      = gamma_d*n13d_m + d3*gamma*n_mfrthrd*n13d2_m

       F_uu   = gamma_uu*n13d_m + d3*(gamma_u+gamma_u)*n_mfrthrd*n13d2_m
     &        - (4.0d0/9.0d0)*d*gamma*n_mfrthrd*n_m*n13d2_m
     &        + (2.0d0/9.0d0)*d*d*gamma*n_mfrthrd*n_mfrthrd*n13d3_m

       F_dd   = gamma_dd*n13d_m + d3*(gamma_d+gamma_d)*n_mfrthrd*n13d2_m
     &        - (4.0d0/9.0d0)*d*gamma*n_mfrthrd*n_m*n13d2_m
     &        + (2.0d0/9.0d0)*d*d*gamma*n_mfrthrd*n_mfrthrd*n13d3_m

       F_du   = gamma_du*n13d_m + d3*(gamma_u+gamma_d)*n_mfrthrd*n13d2_m
     &        - (4.0d0/9.0d0)*d*gamma*n_mfrthrd*n_m*n13d2_m
     &        + (2.0d0/9.0d0)*d*d*gamma*n_mfrthrd*n_mfrthrd*n13d3_m


       enthrd   = dexp(-c*n_mthrd)
       Q        = enthrd*n_mfvthrd
       Q1       = (1.0d0/3.0d0)*c*n_mfrthrd  *enthrd
     &          - (5.0d0/3.0d0)*n_mfvthrd*n_m*enthrd
       G        = F*Q

       P  = (1.0d0/3.0d0)*c*n_mfrthrd - (5.0d0/3.0d0)*n_m
       P1 = ((-4.0d0/9.0d0)*c*n_mfrthrd + (5.d0/3.0d0)*n_m)*n_m

       G_u      = F_u*Q + G*P
       G_d      = F_d*Q + G*P

       G_uu     = F_uu*Q + F_u*Q1 + G_u*P + G*P1
       G_dd     = F_dd*Q + F_d*Q1 + G_d*P + G*P1
       G_du     = F_du*Q + F_d*Q1 + G_u*P + G*P1



       fclda = -a*F*n - 2.0d0*a*b*G*Cf*(2.0d0**twthrd)
     &         *(nup_etthrd + ndn_etthrd)

       fclda_u = -a*F_u*n - a*F 
     &           - 2.0d0*a*b*G_u*Cf*(2.0d0**twthrd)
     &           *(nup_etthrd + ndn_etthrd) 
     &           - 2.0d0*a*b*G*Cf*(2.0d0**twthrd)
     &           *(8.0d0/3.0d0)*nup_fvthrd

       fclda_d = -a*F_d*n - a*F 
     &           - 2.0d0*a*b*G_d*Cf*(2.0d0**twthrd)
     &           *(nup_etthrd + ndn_etthrd) 
     &           - 2.0d0*a*b*G*Cf*(2.0d0**twthrd)
     &           *(8.0d0/3.0d0)*ndn_fvthrd



       H0 = (a*b/2.0d0)*(G 
     &    + (1.0d0/3.0d0)*(nup*G_d + ndn*G_u) 
     &    + (1.0d0/4.0d0)*(nup*G_u + ndn*G_d))
       H0_u=(a*b/2.0d0)*(G_u
     &     + (1.0d0/3.0d0)*(G_d + nup*G_du + ndn*G_uu)
     &     + (1.0d0/4.0d0)*(G_u + nup*G_uu + ndn*G_du))
       H0_d=(a*b/2.0d0)*(G_d
     &     + (1.0d0/3.0d0)*(G_u + ndn*G_du + nup*G_dd)
     &     + (1.0d0/4.0d0)*(G_d + ndn*G_dd + nup*G_du))


       Hu = (a*b/18.0d0)*(G + (15.0d0/4.0d0)*nup*G_u 
     &                      -  (9.0d0/4.0d0)*ndn*G_d
     &                      -        (3.0d0)*nup*G_d 
     &                      +  (3.0d0/2.0d0)*ndn*G_u)

       Hu_u = (a*b/18.0d0)*(G_u + (15.0d0/4.0d0)*(G_u + nup*G_uu)
     &                          -  (9.0d0/4.0d0)*(ndn*G_du)
     &                          -        (3.0d0)*(G_d + nup*G_du) 
     &                          +  (3.0d0/2.0d0)*(ndn*G_uu))

       Hu_d = (a*b/18.0d0)*(G_d + (15.0d0/4.0d0)*(nup*G_du)
     &                          -  (9.0d0/4.0d0)*(G_d+ndn*G_dd)
     &                          -        (3.0d0)*(nup*G_dd)
     &                          +  (3.0d0/2.0d0)*(G_u + ndn*G_du))



       Hd = (a*b/18.0d0)*(G + (15.0d0/4.0d0)*ndn*G_d 
     &                      -  (9.0d0/4.0d0)*nup*G_u
     &                      -        (3.0d0)*ndn*G_u 
     &                      +  (3.0d0/2.0d0)*nup*G_d)

       Hd_d = (a*b/18.0d0)*(G_d + (15.0d0/4.0d0)*(G_d + ndn*G_dd)
     &                          -  (9.0d0/4.0d0)*(nup*G_du)
     &                          -        (3.0d0)*(G_u + ndn*G_du) 
     &                          +  (3.0d0/2.0d0)*(nup*G_dd))

       Hd_u = (a*b/18.0d0)*(G_u + (15.0d0/4.0d0)*(ndn*G_du)
     &                          -  (9.0d0/4.0d0)*(G_u+nup*G_uu)
     &                          -        (3.0d0)*(ndn*G_uu)
     &                          +  (3.0d0/2.0d0)*(G_d + nup*G_du))


       fc = fclda + H0*agr2 + Hu*agrup2 + Hd*agrdn2 

*    ***calculate derivatives w.r.t up and down density
       fc_u = fclda_u + H0_u*agr2 + Hu_u*agrup2 + Hd_u*agrdn2
       fc_d = fclda_d + H0_d*agr2 + Hu_d*agrup2 + Hd_d*agrdn2

*    ***calculate derivatives w.r.t. up,down and total density gradients 
       fc_agr   = 2.0d0*H0*agr
       fc_agrup = 2.0d0*Hu*agrup
       fc_agrdn = 2.0d0*Hd*agrdn


*      ******DEBUG*********
       fc = fclda + 0.001d0*agr2 + 0.002d0*(agrup2 + agrdn2)
       fc_u = fclda_u
       fc_d = fclda_d
       fc_agr   = 2.0d0*agr  *0.001d0
       fc_agrup = 2.0d0*agrup*0.002d0
       fc_agrdn = 2.0d0*agrdn*0.002d0
*      ******DEBUG*********
   
*    ***correlation energy dentsity
       ce = fc/n




*      *******end correlation part********

c       write(*,*) "n,nup,ndn,agrup,agrdn,:",j,n,nup,ndn,agrup,agrdn,agr
c       write(*,*) "xe,ce         :",j,xe,ce
c       write(*,*) "fdnxup,fdncup     :",j,fdnxup,fc_u
c       write(*,*) "fdnxdn,fdncdn     :",j,fdnxdn,fc_d
c       write(*,*) "fdagrxup,fdargrcup:",j,fdagrxup,
c     >                                fc_agrup,fc_agr
c       write(*,*) "fdagrxdn,fdargrcdn:",j,fdagrxdn,
c     >                                fc_agrdn,fc_agr
c
c       write(*,*) "restricted fdagrc",fc_agr+0.5d0*(fc_agrdn+fc_agrup)
c       write(*,*) "fc_lda,fdnc_lda:",fclda,fclda_u,fclda_d
c       write(*,*) "Hu,Hd,H0:",Hu,Hd,H0
c       write(*,*) "Ho:",H0+0.25d0*Hu+0.25d0*Hd
c       write(*,*) "Ho_n:",0.5d0*(H0_u+H0_d+0.25d0*(Hu_u+Hu_d+Hd_d+Hd_u))
c       write(*,*) "F,G:",F,G
c       write(*,*) "F_u,G_u:",F_u,G_u
c       write(*,*) "F_d,G_d:",F_d,G_d
c       write(*,*) "G_uu,G_dd,G_du:",G_uu,G_dd,G_du
c       write(*,*) "F_uu,F_dd,F_du:",F_uu,F_dd,F_du
c       write(*,*) "Gamma's:",gamma,gamma_u,gamma_d,
c     >             gamma_uu,gamma_dd,gamma_du
c       write(*,*) "n13d",n13d
c       write(*,*) "0.5/nup2,0.5/ndn2:",0.5d0/nup2,0.5d0/ndn2
c       write(*,*) "fc:",fc
c       write(*,*) "0.5d0*(F_uu+F_du)",0.5d0*(F_uu + F_du)
c       write(*,*) "0.5d0*(G_uu+G_du)",0.5d0*(G_uu + G_du)
c       write(*,*) "nup*G_u...",nup*G_u,ndn*G_d,nup*G_d,ndn*G_u
c       write(*,*)



*      ***return blyp exchange correlation values*** 
       xce(j) = x_parameter*xe + c_parameter*ce
       fn(j,1) = x_parameter*fdnxup + c_parameter*fc_u
       fn(j,2) = x_parameter*fdnxdn + c_parameter*fc_d

       fdn(j,1) = x_parameter*fdagrxup + c_parameter*fc_agrup
       fdn(j,2) = x_parameter*fdagrxdn + c_parameter*fc_agrdn
       fdn(j,3) =                        c_parameter*fc_agr

   


      end do
      return
      end

*     **************************************************
*     *                                                *
*     *            gen_BLYP_BW_restricted              *
*     *                                                *
*     **************************************************

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
      real*8    Fc,Gc,C1,C2
      real*8    n, n_thrd,n_fv,n_fr,n_tw
      real*8    n_m,n_mthrd,n_mfrthrd,n_mfvthrd
      real*8    agr,agr2, chi, chi2,chiSQ,sd
      real*8    K
      real*8    F1, F2      
      real*8    xe,fdnx,fdagrx
      real*8    ce,fdnc, fdagrc,fdnc_lda

      real*8    fc_lda,Ho,Ho_n,Gc_n,Gc_nn,Fc_n,Fc_nn
      real*8    P,P_n
      
******* constants **********************************
      real*8 pi,thrd,two_thrd
      parameter (pi=3.14159265358979311599d0)
      parameter (thrd=1.0d0/3.0d0)
      parameter (two_thrd=1.25992104989487319066d0)  ! two_thrd = 2**thrd
    

*******density cutoff parameters********************
      real*8 DNS_CUT, ETA
      parameter (DNS_CUT        =      1.0d-18)
      parameter (ETA            =      1.0d-20)
****** Becke constants *****************************
      real*8 lda_c,beta
      parameter (beta = 0.0042d0)
      parameter (lda_c = 0.93052573634910018540d0)   ! lda_c = (3/2)*(3/(4*pi))**(1/3)
*******LYP correlation parameters a, b, c, d********
      real*8 a,b,c,d,Cf
      parameter (a              =      0.04918d0)
      parameter (b              =      0.132d0)
      parameter (c              =      0.2533d0)
      parameter (d              =      0.349d0)
      parameter (Cf =2.87123400018819108225d0)     ! Cf = (3/10)*(3*pi*pi)**(2/3)

****** collated LYP parameters *********************
      real*8 ho1,ho2,ho3,thrd_d,abCf,p2,p3,p4,p5,p6
      parameter (ho1=(19.0d0/36.0d0)*a*b)
      parameter (ho2=( 7.0d0/24.0d0)*a*b)
      parameter (ho3=(59.0d0/72.0d0)*a*b)
      parameter (thrd_d=thrd*d)
      parameter (abCf=a*b*Cf)
      parameter (p2=-5.0d0*thrd,p3=c*thrd)
      parameter (p4=-4.0d0*thrd*thrd_d)
      parameter (p5=5.0d0*thrd)
      parameter (p6=-4.0d0*thrd*thrd*c)
      
****************************************************



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

       n_m     = 1.0d0/n
       n_mthrd = 1.0d0/n_thrd
       n_mfrthrd = n_mthrd*n_m
       n_mfvthrd = n_mfrthrd*n_mthrd

       n_fv = n_thrd*n_thrd*n_thrd*n_thrd*n_thrd
       n_fr = n_thrd*n_thrd*n_thrd*n_thrd
       n_tw = n_thrd*n_thrd


       Fc = (1.0d0/(1.0d0+d*n_mthrd))
       Fc_n = thrd_d*n_mfrthrd*Fc*Fc
       Fc_nn = thrd_d*(-4.0d0*thrd)*n_mfrthrd*n_m*Fc*Fc
     >       + thrd_d*n_mfrthrd*2.0d0*Fc*Fc_n

       Gc = Fc*dexp(-c*n_mthrd)*n_mfvthrd

       P  = (thrd_d*Fc*n_mfrthrd - 5.0d0*thrd*n_m + c*thrd*n_mfrthrd)

       P_n = thrd_d*Fc_n*n_mfrthrd
     >     - thrd_d*Fc  *n_mfrthrd*(4.0d0*thrd*n_m)
     >     + 5.0d0*thrd*n_m*n_m
     >     - 4.0d0*thrd*thrd*c*n_mfrthrd*n_m

c       P  = (thrd_d*Fc*n_mfrthrd + p2*n_m + p3*n_mfrthrd)
c       P_n = thrd_d*Fc_n*n_mfrthrd
c     >     + p4*Fc*n_mfrthrd*n_m
c     >     + p5*n_m*n_m
c     >     + p6*n_mfrthrd*n_m

       Gc_n  = Gc*P
       Gc_nn = Gc*P*P + Gc*P_n

       
       fc_lda   = -a*Fc*n 
     >          - abCf*Gc*(n_fr*n_fr)
       fdnc_lda = -a*Fc_n*n 
     >          - a*Fc
     >          - abCf*Gc_n*n_fr*n_fr
     >          - 8.0d0*thrd*abCf*Gc*n_fv

c      Ho   = (19.0d0/36.0d0)*(a*b*Gc)   + (7.0d0/24.0d0)*a*b*Gc_n*n
c      Ho_n = (59.0d0/72.0d0)*(a*b*Gc_n) + (7.0d0/24.0d0)*a*b*Gc_nn*n
       Ho   = ho1*Gc   + ho2*Gc_n*n
       Ho_n = ho3*Gc_n + ho2*Gc_nn*n

       ce = (fc_lda + Ho*agr2)*n_m
       fdnc = fdnc_lda + Ho_n*agr2

       fdagrc = 2.0d0*Ho*agr

c       write(*,*) "n,agr         :",i,n,agr
c       write(*,*) "xe,ce         :",i,xe,ce
c       write(*,*) "fdnx,fdnc     :",i,fdnx,fdnc
c       write(*,*) "fdagrx,fdargrc:",i,fdagrx,fdagrc
c       write(*,*) "fc_lda,fdnc_lda:",fc_lda,fdnc_lda
c       write(*,*) "Ho :",Ho
c       write(*,*) "Ho_n:",Ho_n
c       write(*,*) "F,G:",Fc,Gc
c       write(*,*) "F_n,G_n:",Fc_n,Gc_n
c       write(*,*) "G_nn:",Gc_nn
c       write(*,*) "F_nn:",Fc_nn
c       write(*,*) "n13d:",1.0d0+d*n_mthrd
c       write(*,*) "fc:",fc_lda + Ho*agr2
c       write(*,*)

 


       xce(i) 	= x_parameter*xe   + c_parameter*ce
       fn(i)	= x_parameter*fdnx  + c_parameter*fdnc
       fdn(i) 	= x_parameter*fdagrx + c_parameter*fdagrc

      end do
      return
      end
