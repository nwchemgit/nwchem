
*      **********************************************
*      *                                            *
*      *               FsrxHSE                      *
*      *                                            *
*      **********************************************
c
c HSE evaluates the Heyd et al. Screened Coulomb
c Exchange Functional
c
c Calculates the enhancement factor
c
      subroutine FsrxHSE(s,nu,Fxhse,d10Fxhse,d01Fxhse)
      implicit none
      double precision s,nu,Fxhse,d10Fxhse,d01Fxhse
c

      double precision  A,B,C,D,E
      double precision  ha2,ha3,ha4,ha5,ha6,ha7
      double precision  hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9
      double precision  smax,strans,sconst
c
      double precision  zero,one,two,three,four,five,six,seven,eight
      double precision  nine,ten
      double precision  fifteen,sixteen

      double precision  H,hnum,hden 
      double precision  d1H,d1hnum,d1hden 
      double precision  s2,s3,s4,s5,s6,s7,s8,s9
      double precision  Fs, d1Fs
      double precision  zeta, lambda, eta, kf, chi, lambda2
      double precision  d1zeta,d1lambda,d1eta,d1nu,d1chi,d1lambda2
      double precision  EGs,d1EGs
      double precision  nu2,L2,L3,nu3,nu4,nu5,nu6
      double precision  Js,Ks,Ms,Ns
      double precision  d1Js,d1Ks,d1Ms,d1Ns

      double precision  tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
      double precision  tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15
      double precision  Fxhse1,Fxhse2,Fxhse3,Fxhse4,Fxhse5,Fxhse6
      double precision  d1Fxhse1,d1Fxhse2,d1Fxhse3,d1Fxhse4,d1Fxhse5
      double precision  d1Fxhse6,d1Fxhse7

      double precision  r42,r27,r12,r15,r14,r18,r20,r30,r56,r72
      double precision  r16,r32,r24,r48,r11,r64,r35
      double precision  pi,pi2,srpi,s02
      double precision  f12,f13,f32,f52,f72,f92

c
c     Constants for HJS hole
c
      Data A,B,C,D,E
     &     / 7.57211D-1,-1.06364D-1,-1.18649D-1,
     &       6.09650D-1,-4.77963D-2 /
c
c     Constants for fit of H(s) (PBE hole)
c     Taken from JCTC_5_754 (2009)
c
      Data ha2,ha3,ha4,ha5,ha6,ha7
     &     / 1.59941D-2,8.52995D-2,-1.60368D-1,1.52645D-1,
     &      -9.71263D-2,4.22061D-2 /
c
      Data hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9
     &     / 5.33319D0,-12.4780D0,11.0988D0,-5.11013D0,
     &      1.71468D0,-6.10380D-1,3.07555D-1,-7.70547D-2,
     &      3.34840D-2 /

c
c     Whole numbers used during evaluation
c
      Data zero,one,two,three,four,five,six,seven,eight,nine,ten
     &     / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
       
      Data r11,r12,r14,r15,r16,r18,r20,r24,r27,r30,r32
     &     / 11D0,12D0,14D0,15D0,16D0,18D0,20D0,24D0,27d0,30D0,32D0 /

      Data r35,r42,r48,r56,r64,r72
     &     / 35D0,42D0,48D0,56D0,64D0,72D0 /
c
c     Fractions used during evaluation
c
      Data f12
     &     / 0.5D0 /
c
c     General constants
c
      f13   = one/three
      f32   = three/two
      f52   = five/two
      f72   = seven/two
      f92   = nine/two
      pi    = ACos(-one)
      pi2   = pi*pi
      srpi = dsqrt(pi)
c
c
c     Calculate prelim variables
c
      s2 = s*s
      s02 = s2/four
      s3 = s2*s
      s4 = s3*s
      s5 = s4*s
      s6 = s5*s
      s7 = s6*s
      s8 = s7*s
      s9 = s8*s

c
c     Calculate H(s) the model exhange hole
c
      hnum = ha2*s2 + ha3*s3 + ha4*s4 + ha5*s5 + ha6*s6 + ha7*s7 
      hden = one + hb1*s + hb2*s2 + hb3*s3 + hb4*s4 + hb5*s5 +
     &       hb6*s6 + hb7*s7 + hb8*s8 + hb9*s9
      H = hnum/hden

c
c     Calculate helper variables
c
      zeta = s2*H
      eta = A + zeta
      lambda = D + zeta

c      kf = (three*pi2*rho)**f13 
c      if (ipol.eq.1) then
c         kf = (three*pi2*rho)**f13 
c      else
c         kf = (six*pi2*rho)**f13 
c      endif
c      nu = omega/kf

      chi = nu/dsqrt(lambda+nu**two)
      lambda2 = (one+chi)*(lambda+nu**two)

c
c     Calculate F(H(s)) for the model exhange hole
c
      Fs = one-s2/(r27*C*(one+s02))-zeta/(two*C)

c
c     Calculate EG(s) 
c
      EGs = -(two/five)*C*Fs*lambda - (four/r15)*B*lambda**two -
     &      (six/five)*A*lambda**three - 
     &      (four/five)*srpi*lambda**(seven/two) -
     &      (r12/five)*(lambda**(seven/two))*(dsqrt(zeta)-dsqrt(eta))
 
c
c     Calculate the denominators needed
c

      nu2 = nu*nu
      Js = (dsqrt(zeta+nu2)+dsqrt(eta+nu2))*(dsqrt(zeta+nu2)+nu) 
      Ks = (dsqrt(zeta+nu2)+dsqrt(eta+nu2))*(dsqrt(eta+nu2)+nu) 
      Ms = (dsqrt(zeta+nu2)+dsqrt(lambda+nu2))*(dsqrt(lambda+nu2)+nu) 
      Ns = (dsqrt(eta+nu2)+dsqrt(lambda+nu2))*(dsqrt(lambda+nu2)+nu) 

c
c       The final value for the enhancement factor is
c
      tmp1 = one + f12*chi
      tmp2 = one + (nine/eight)*chi + (three/eight)*chi**two 
      Fxhse1  = A*(zeta/Js + eta/Ks) 
      Fxhse2  = -(four/nine)*B/lambda2
      Fxhse3  = -(four/nine)*C*Fs*tmp1/lambda2**two
      Fxhse4  = -(eight/nine)*EGs*tmp2/lambda2**three
      Fxhse5  = two*zeta*dlog(one -D/Ms)
      Fxhse6  = -two*eta*dlog(one -(D-A)/Ns)

      Fxhse = Fxhse1+Fxhse2+Fxhse3+Fxhse4+Fxhse5+Fxhse6
c
c     Calculate the first derivative of H with respect to the
c     reduced density gradient, s.
c
      d1hnum = two*ha2*s + three*ha3*s2 + four*ha4*s3 +
     &          five*ha5*s4 + six*ha6*s5 + seven*ha7*s6

      d1hden  = hb1 + two*hb2*s +three*hb3*s2 + four*hb4*s3 +
     &          five*hb5*s4 + six*hb6*s5 + seven*hb7*s6 +
     &          eight*hb8*s7 + nine*hb9*s8 
      d1H =   (hden*d1hnum -hnum*d1hden)/hden**two

c
c     calculate first derivative of variables needed with respect to s
c 
      d1zeta = two*s*H + s2*d1H
      d1eta  = d1zeta
      d1lambda = d1zeta
      d1chi = -f12*nu*d1zeta/(lambda + nu2)**f32
      d1lambda2 = d1chi*(lambda + nu**two) + (one+chi)*d1lambda
      !d1lambda2 = (d1lambda*(one-chi)+lambda*d1chi)/(one-chi)**two

c   
c     calculate the first derivative of Fs with respect to s
c   
      d1Fs = -two*s/(r27*C*(one+s02)**two) - d1zeta/(two*C)

c
c     Calculate the first derivate of EGs with respect to s
c
      d1EGs = -(two/five)*C*(d1Fs*lambda + Fs*d1lambda) -
     &        (eight/r15)*B*lambda*d1lambda -
     &        (r18/five)*A*lambda*lambda*d1lambda -
     &        (r14/five)*srpi*d1lambda*lambda**f52 -
     &        (r42/five)*(lambda**f52)*
     &        d1lambda*(dsqrt(zeta)-dsqrt(eta))-
     &        (six/five)*(lambda**(seven/two))*
     &        (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))

c
c     Calculate the first derivate of denominators needed with respect
c     to s
c
      tmp1 = (dsqrt(zeta+nu2)+nu)/(dsqrt(eta+nu2)) 
      tmp2 = (dsqrt(eta+nu2)+nu)/(dsqrt(zeta+nu2))

      d1Js = f12*d1zeta*(two+tmp1+tmp2)
      d1Ks = d1Js

      tmp3 = (dsqrt(zeta+nu2)+nu)/(dsqrt(lambda+nu2))
      tmp4 = (dsqrt(lambda+nu2)+nu)/(dsqrt(zeta+nu2)) 
      d1Ms = f12*d1zeta*(two +tmp3+tmp4)

      tmp5 = (dsqrt(lambda+nu2)+nu)/(dsqrt(eta+nu2))
      tmp6 = (dsqrt(eta+nu2)+nu)/(dsqrt(lambda+nu2))
      d1Ns = f12*d1zeta*(two + tmp5+tmp6)
c
c
c     Calculate the derivative of the 08-Fxhse with respect to s
c
      L2 = lambda2*lambda2
      L3 = lambda2*lambda2*lambda2
      d1Fxhse1  = A*( (Js*d1zeta - zeta*d1Js)/(Js*Js) +
     &                (Ks*d1zeta - eta*d1Ks)/(Ks*Ks) ) 

      d1Fxhse2  = (four/nine)*B*d1lambda2/L2 

      tmp9 = d1lambda2/lambda2
      tmp7 = d1Fs - two*Fs*tmp9
      tmp8 = one + f12*chi
      tmp10 =  f12*Fs*d1chi

      d1Fxhse3 = -(four*C/(nine*L2))*(tmp7*tmp8+tmp10)

c      d1Fxhse3  = -(four/nine)*(C/(L2*L2))* 
c     &           (L2*d1Fs - two*lambda2*Fs*d1lambda2 +
c     &           f12*(L2*(d1Fs*chi+Fs*d1chi) -
c     &           two*lambda2*chi*Fs*d1chi))

c       tmp4 = d1chi/(one+chi) + d1lambda/(lambda+nu2)
c       tmp4 = (one-chi)*d1lambda2/lambda
c       tmp1 = (eight/three+three*chi+chi*chi)*tmp4
c       tmp2 = d1chi + (two/three)*chi*d1chi
c       tmp3 = (eight/nine + chi + f13*chi*chi)
c       d1FXhse4 = ((tmp1-tmp2)*EGs-tmp3*d1EGs)/L3

        tmp7 = one + (nine/eight)*chi+(three/eight)*chi*chi
        tmp8 = (nine/eight)*d1chi + (six/eight)*chi*d1chi

       d1Fxhse4 = -(eight/(nine*L3))*((d1EGs-three*EGs*tmp9)*tmp7
     &           + EGs*tmp8)
c      d1Fxhse4  = -(eight/nine)*(L3*d1EGs - 
c     &            three*EGs*L2*d1lambda2)/(L3*L3) -
c     &             (L3*(d1EGs*chi + EGs*d1chi) -
c     &            three*EGs*chi*L2*d1lambda2)/(L3*L3)-
c     &           (L3*(d1EGs*chi*chi+two*EGs*chi*d1chi)-
c     &            three*EGs*chi*chi*L2*d1lambda2)/(three*L3*L3) 
c
      d1Fxhse5  = two*d1zeta*dlog(one-D/Ms) +
     &           two*zeta*D*d1Ms/(Ms*Ms*(one-D/Ms)) 

      d1Fxhse6  = -two*d1eta*dlog(one- (D-A)/Ns) -
     &           two*eta*(D-A)*d1Ns/(Ns*Ns*(one-(D-A)/Ns)) 
c
      d10Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6
c
c     Calculate the derivative of 08-Fxhse with respect to nu
c
      nu3 = nu2*nu
c
      d1Fxhse1 = -((A*(nu + dsqrt(eta + nu2))*zeta)/
     &            (dsqrt(eta + nu2)*dsqrt(nu2 + zeta)*
     &            (nu + dsqrt(nu2 + zeta))*
     &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))))
c
      d1Fxhse2 = -((A*eta*(nu/dsqrt(eta + nu2) + nu/
     &            dsqrt(nu2 + zeta)))/
     &            ((nu + dsqrt(eta + nu2))*
     &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -
     &            (A*eta*(one + nu/dsqrt(eta + nu2)))/
     &            ((nu + dsqrt(eta + nu2))**two*
     &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta)))
c
      d1Fxhse3 = (four*B)/(nine*(lambda + nu2)**(f32))
c
      d1Fxhse4 = (two*C*Fs)/(three*(lambda + nu2)**(f52))
c
      d1Fxhse5 = (five*EGs*(lambda**two + four*nu3*
     &            (nu + dsqrt(lambda + nu2)) +
     &            lambda*nu*(five*nu + three*dsqrt(lambda + nu2))))/
     &   (three*(lambda + nu2)**four*(nu + dsqrt(lambda + nu2))**three)
c
      d1Fxhse6 = (two*D*zeta*(nu + dsqrt(nu2 + zeta)))/
     &            (dsqrt(lambda + nu2)*dsqrt(nu2 + zeta)*
     &            (-D + lambda + (nu + dsqrt(lambda + nu2))*
     &            (nu + dsqrt(nu2 + zeta))))
c
      d1Fxhse7 = (two*(A - D)*eta*(nu + dsqrt(eta + nu2)))/
     &            (dsqrt(eta + nu2)*dsqrt(lambda + nu2)*
     &            (A - D + lambda + nu2 + nu*dsqrt(eta + nu2) +
     &            nu*dsqrt(lambda + nu2) +
     &            dsqrt(eta + nu2)*dsqrt(lambda + nu2)))
c
      d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+
     &           d1Fxhse7
c   
      return
      end


*    ************************************
*    *					*
*    *	    gen_HSE_BW_restricted	*
*    *					*
*    ************************************
*
*   This routine calculates the HSE exchange-correlation 
*   potential(xcp) and energy density(xce).
*
*
*   Entry - n2ft3d     : number of grid points
*           rho_in(*) :  density (nup+ndn)
*           agr_in(*): |grad rho_in|
*           x_parameter: scale parameter for exchange
*           c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : HSE exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce/d|grad n|
*
      subroutine gen_HSE_BW_restricted(n2ft3d,rho_in,agr_in,
     >                                x_parameter,c_parameter,
     >                                xce,fn,fdn)
      implicit none

      integer    n2ft3d
      real*8     rho_in(n2ft3d)
      real*8     agr_in(n2ft3d)
      real*8     x_parameter,c_parameter
      real*8     xce(n2ft3d)
      real*8     fn(n2ft3d)
      real*8     fdn(n2ft3d)

      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT,ETA
      parameter (DNS_CUT = 1.0d-20)
      parameter (ETA     = 1.0d-20)

c     ***** HSEPBE96 GGA exchange constants ******
      real*8 MU,KAPPA,OMEGA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)
      parameter (OMEGA = 0.2070000000000000d0)
 
c     ****** PBE96 GGA correlation constants ******
      real*8 GAMMA,BETA,BOG
      parameter (GAMMA	= 0.031090690869655d0)
      parameter (BETA	= 0.066724550603149d0)
      parameter (BOG    = BETA/GAMMA)


c     ****** Perdew-Wang92 LDA correlation coefficients *******
      real*8 A_1,A1_1,B1_1,B2_1,B3_1,B4_1     
      parameter (A_1  = 0.0310907d0)
      parameter (A1_1 =	0.2137000d0)
      parameter (B1_1 =	7.5957000d0)
      parameter (B2_1 =	3.5876000d0)
      parameter (B3_1 =	1.6382000d0)
      parameter (B4_1 =	0.4929400d0)

      real*8 A_2,A1_2,B1_2,B2_2,B3_2,B4_2     
      parameter (A_2  =  0.01554535d0)
      parameter (A1_2 =	 0.20548000d0)
      parameter (B1_2 =	14.11890000d0)
      parameter (B2_2 =	 6.19770000d0)
      parameter (B3_2 =	 3.36620000d0)
      parameter (B4_2 =	 0.62517000d0)
      
      real*8 A_3,A1_3,B1_3,B2_3,B3_3,B4_3     
      parameter (A_3  =  0.0168869d0)
      parameter (A1_3 =	 0.1112500d0)
      parameter (B1_3 =	10.3570000d0)
      parameter (B2_3 =	 3.6231000d0)
      parameter (B3_3 =	 0.8802600d0)
      parameter (B4_3 =	 0.4967100d0)

c     **** other constants ****
      real*8 onethird,fourthird,sevensixths
      parameter (onethird=1.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)
      parameter (sevensixths=7.0d0/6.0d0)

c     **** local variables ****
      integer i
      real*8 n,agr
      real*8 kf,ks,s,P0,n_onethird,pi,rs_scale,nu
      real*8 fdnx_const
      real*8 rs,rss,t,t2,t4,t6
      real*8 Q0,Q1,Q2,Q3,Q4,Q5,Q8,Q9,B
      real*8 Ht
      real*8 B_ec,Hrs,H_B
      real*8 F,Fs
      real*8 Fxhse,Fxhse_s,Fxhse_nu

      real*8 ex_lda,ec_lda
      real*8 ec_lda_rs
      real*8 ex,ec,H
      real*8 fnx,fdnx,fnc,fdnc


      pi         = 4.0d0*datan(1.0d0)
      rs_scale   = (0.75d0/pi)**onethird
      fdnx_const = -3.0d0/(8.0d0*pi)
      
      do i=1,n2ft3d
         n     = rho_in(i)+ETA
         agr   = agr_in(i)
        
c        ***** calculate unpolarized Exchange energies and potentials *****
         n_onethird = (3.0d0*n/pi)**onethird
         ex_lda     = -0.75d0*n_onethird

         kf = (3.0d0*pi*pi*n)**onethird
         s  = agr/(2.0d0*kf*n)
         P0 = 1.0d0 + (MU/KAPPA)*s*s

c        if (n.gt.DNS_CUT) then
c           F   = (1.0d0 + KAPPA - KAPPA/P0)
c           Fs  = 2.0d0*MU/(P0*P0)*s
c        else
c           F   = 1.0d0
c           Fs  = 0.0d0
c        end if
         F   = (1.0d0 + KAPPA - KAPPA/P0)
         Fs  = 2.0d0*MU/(P0*P0)*s

*        **** shortrange-HSEPBE96 ****
         nu = OMEGA/kf
         call FsrxHSE(s,nu,Fxhse,Fxhse_s,Fxhse_nu)

         ex   = ex_lda*(F - 0.25d0*Fxhse)
         fnx  = fourthird*(ex - ex_lda*(Fs - 0.25d0*Fxhse_s)*s)
         fnx  = fnx + 0.25d0*onethird*ex_lda*Fxhse_nu*nu
         fdnx = fdnx_const*(Fs - 0.25d0*Fxhse_s)
			

*        *********************************************************************
c        ***** calculate unpolarized correlation energies and potentials *****
*        *********************************************************************

c        **** calculate rs and t ****
         rs    = rs_scale/(n**onethird)
         rss   = dsqrt(rs)

         kf = (3.0d0*pi*pi*n)**onethird
         ks = dsqrt(4.0d0*kf/pi)
         t  = agr/(2.0*ks*n)


c        **** unpolarized LDA correlation energy ****
c        **** ec_p = correlation energy          ****
c        ****   ec_p_rs = dec_p/drs              ****
c        ****   uc_p    = dec_p/dn               ****
         call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ec_lda,ec_lda_rs)
c        **** PBE96 correlation energy  corrections ****
         t2 = t*t
         t4 = t2*t2
         B = -ec_lda/GAMMA
         B = BOG/(exp(B)-1.0d0+ETA)
         Q4 = 1.0d0 + B*t2
         Q5 = 1.0d0 + B*t2 + B*B*t4
         H = GAMMA*dlog(1.0d0 + BOG*Q4*t2/Q5)


c        **** PBE96 correlation fdn and fdnc derivatives ****
         t6   = t4*t2

         B_ec = (B/BETA)*(BOG+B)

         Q8  = Q5*Q5+BOG*Q4*Q5*t2
         Q9  = 1.0d0+2*B*t2
         H_B  = -BETA*B*t6*(2.0d0+B*t2)/Q8 
         Hrs  = H_B*B_ec*ec_lda_rs

         Ht  = 2.0d0*BETA*Q9/Q8*t

         ec   = ec_lda + H
         fnc = ec  - (onethird*rs*ec_lda_rs)
     >             - (onethird*rs*Hrs)
     >             - (sevensixths*t*Ht)
         fdnc = 0.5d0* Ht/ks

         xce(i) = x_parameter*ex   + c_parameter*ec
         fn(i)  = x_parameter*fnx  + c_parameter*fnc
         fdn(i) = x_parameter*fdnx + c_parameter*fdnc
         

c       write(*,*) "pbe96:",i,ec,fnc,fdnc


      end do

      return
      end
      
