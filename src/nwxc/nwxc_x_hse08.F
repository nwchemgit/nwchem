#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_hse08.F
C> The range separated HSE exchange enhancement factor
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the range separated HSE exchange enhancement factor
C>
C> This routine calculates the HSE exchange enhancement factor [1,2].
C>
C> ### References ###
C>
C> [1] J. Heyd, G.E. Scuceria, M. Ernzerhof, "Hybrid functionals based
C> on a screened Coulomb potential", J. Chem. Phys. <b>118</b>, 
C> 8207 (2003), DOI: <a href="https://doi.org/10.1063/1.1564060">
C> 10.1063/1.1564060</a>.
C>
C> [2] J. Heyd, G.E. Scuceria, M. Ernzerhof, "Erratum: Hybrid 
C> functionals based on a screened Coulomb potential", J. Chem. Phys.
C> <b>124</b>, 219906 (2006), DOI:
C> <a href="https://doi.org/10.1063/1.2204597">
C> 10.1063/1.2204597</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_HSE08_p(cam_omega,ipol,rho,s,Fxhse)
#else
      Subroutine nwxc_x_HSE08(cam_omega,ipol,rho,s,Fxhse)
c    &           d10Fxhse,d01Fxhse)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_HSE08_d2(cam_omega,ipol,rho,s,Fxhse)
c    &           d10Fxhse,d01Fxhse,d20Fxhse,d02Fxhse,d11Fxhse)
#else
      Subroutine nwxc_x_HSE08_d3(cam_omega,ipol,rho,s,Fxhse)
c    &           d10Fxhse,d01Fxhse,d20Fxhse,d02Fxhse,d11Fxhse,
c    &           d30Fxhse,d21Fxhse,d12Fxhse,d03Fxhse)
#endif
c
#include "nwad.fh"
c
      implicit none
c
c HSE evaluates the Heyd et al. Screened Coulomb
c Exchange Functional
c
c Calculates the enhancement factor
c
      integer ipol
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::cam_omega
#else
      double precision cam_omega
#endif
#else
      double precision cam_omega
#endif
      type(nwad_dble),intent(in) ::rho,s
      type(nwad_dble),intent(out)::Fxhse
c     double precision d10Fxhse,d01Fxhse
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision d20Fxhse,d02Fxhse,d11Fxhse
#endif
#if defined(THIRD_DERIV)
c     double precision d30Fxhse,d21Fxhse,d12Fxhse,d03Fxhse
#endif

      double precision  A,B,C,D,E
      double precision  ha2,ha3,ha4,ha5,ha6,ha7
      double precision  hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9
      double precision  smax,strans,sconst
c
      double precision  zero,one,two,three,four,five,six,seven,eight
      double precision  nine,ten
      double precision  fifteen,sixteen

      type(nwad_dble):: H,hnum,hden 
      double precision  d1H,d1hnum,d1hden 
      type(nwad_dble):: s2,s3,s4,s5,s6,s7,s8,s9,s02
      type(nwad_dble):: zeta,eta,kf,Fs,nu,nu2,lambda,lambda2,chi
      double precision  d1Fs
      double precision  d1zeta,d1lambda,d1eta,d1nu
      double precision  d10chi,d10lambda2,d01chi,d01lambda2
      type(nwad_dble):: EGs
      double precision  d1EGs
      double precision  L2,L3,nu3,nu4,nu5,nu6,nu7,nu8
      type(nwad_dble):: Js,Ks,Ms,Ns
      double precision  Js2,Ks2,Ms2,Ns2
      double precision  d10Js,d10Ks,d10Ms,d10Ns
      double precision  d01Js,d01Ks,d01Ms,d01Ns

      type(nwad_dble):: tmp1,tmp2
      double precision  tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
      double precision  tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15
      type(nwad_dble):: Fxhse1,Fxhse2,Fxhse3,Fxhse4,Fxhse5,Fxhse6
      double precision  d1Fxhse1,d1Fxhse2,d1Fxhse3,d1Fxhse4,d1Fxhse5
      double precision  d1Fxhse6,d1Fxhse7

      double precision  r42,r27,r12,r15,r14,r18,r20,r30,r56,r72
      double precision  r16,r32,r24,r48,r11,r64,r35,r60,r36,r63
      double precision  r189,r256,r120,r210,r336,r504,r21,r105,r126
      double precision  pi,pi2,srpi
      double precision  f12,f13,f32,f52,f72,f92
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision  d2H,d2hnum,d2hden
c     double precision  d2Fs
c     double precision  d2zeta,d2lambda,d2eta,d2nu
c     double precision  d20chi,d20lambda2,d02chi,d02lambda2
c     double precision  d11chi,d11lambda2
c     double precision  d2EGs
c     double precision  L5
c     double precision  Js3,Ks3,Ms4,Ns4
c     double precision  d20Js,d20Ks,d20Ms,d20Ns
c     double precision  d02Js,d02Ks,d02Ms,d02Ns
c     double precision  d11Js,d11Ks,d11Ms,d11Ns
c     double precision  d2Fxhse1,d2Fxhse2,d2Fxhse3,d2Fxhse4,d2Fxhse5
c     double precision  d2Fxhse6,d2Fxhse7
c     double precision  d11Fxhse1,d11Fxhse2,d11Fxhse3,d11Fxhse4
c     double precision  d11Fxhse5,d11Fxhse6,d11Fxhse7
#endif
#if defined(THIRD_DERIV)
c     double precision  d3H,d3hnum,d3hden
c     double precision  d3Fs
c     double precision  d3zeta,d3lambda,d3eta,d3nu,d30chi,d30lambda2
c     double precision  d03chi,d03lambda2,d21chi,d12chi
c     double precision  d21lambda2,d12lambda2
c     double precision  L4,L6
c     double precision  Js4,Ks4,Ms3,Ns3
c     double precision  d3EGs
c     double precision  d30Js,d30Ks,d30Ms,d30Ns
c     double precision  d03Js,d03Ks,d03Ms,d03Ns
c     double precision  d21Js,d21Ks,d21Ms,d21Ns
c     double precision  d12Js,d12Ks,d12Ms,d12Ns
c     double precision  d3Fxhse1,d3Fxhse2,d3Fxhse3,d3Fxhse4,d3Fxhse5
c     double precision  d3Fxhse6,d3Fxhse7
c     double precision  d21Fxhse1,d21Fxhse2,d21Fxhse3,d21Fxhse4
c     double precision  d21Fxhse5,d21Fxhse6,d21Fxhse7
c     double precision  d12Fxhse1,d12Fxhse2,d12Fxhse3,d12Fxhse4
c     double precision  d12Fxhse5,d12Fxhse6,d12Fxhse7
#endif
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

      Data r21,r36,r60,r63,r120,r189,r210,r256,r336,r504
     1     / 21.D0,36.D0,60.D0,63.D0,120.D0,189.D0,210.D0,256.D0,
     2       336.D0,504.D0 /

      Data r105,r126
     1     / 105.D0,126.D0 /
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
!*********************************
! Calculate the Enhancement Factor
!*********************************
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
      if (ipol.eq.1) then
         kf = (three*pi2*rho)**f13 
      else
         kf = (six*pi2*rho)**f13 
      endif
      nu = cam_omega/kf
      chi = nu/sqrt(lambda+nu**2)
      lambda2 = (one+chi)*(lambda+nu**2)

c
c     Calculate F(H(s)) for the model exhange hole
c
      Fs = one-s2/(r27*C*(one+s02))-zeta/(two*C)

c
c     Calculate EG(s) 
c
      EGs = -(two/five)*C*Fs*lambda - (four/r15)*B*lambda**2 -
     &      (six/five)*A*lambda**3 - 
     &      (four/five)*srpi*lambda**(seven/two) -
     &      (r12/five)*(lambda**(seven/two))*(sqrt(zeta)-sqrt(eta))
 
c
c     Calculate the denominators needed
c

      nu2 = nu*nu
      Js = (sqrt(zeta+nu2)+sqrt(eta+nu2))*(sqrt(zeta+nu2)+nu) 
      Ks = (sqrt(zeta+nu2)+sqrt(eta+nu2))*(sqrt(eta+nu2)+nu) 
      Ms = (sqrt(zeta+nu2)+sqrt(lambda+nu2))*(sqrt(lambda+nu2)+nu) 
      Ns = (sqrt(eta+nu2)+sqrt(lambda+nu2))*(sqrt(lambda+nu2)+nu) 

c
c       The final value for the enhancement factor is
c
      tmp1 = one + f12*chi
      tmp2 = one + (nine/eight)*chi + (three/eight)*chi**2 
      Fxhse1  = A*(zeta/Js + eta/Ks) 
      Fxhse2  = -(four/nine)*B/lambda2
      Fxhse3  = -(four/nine)*C*Fs*tmp1/lambda2**2
      Fxhse4  = -(eight/nine)*EGs*tmp2/lambda2**3
      Fxhse5  = two*zeta*log(one -D/Ms)
      Fxhse6  = -two*eta*log(one -(D-A)/Ns)
c
      Fxhse = Fxhse1+Fxhse2+Fxhse3+Fxhse4+Fxhse5+Fxhse6
c
!*********************************************************
! Calculate the First Derivative of the Enhancement Factor
!*********************************************************
c
c     Calculate the first derivative of H with respect to the
c     reduced density gradient, s.
c
c     d1hnum = two*ha2*s + three*ha3*s2 + four*ha4*s3 +
c    &          five*ha5*s4 + six*ha6*s5 + seven*ha7*s6

c     d1hden  = hb1 + two*hb2*s +three*hb3*s2 + four*hb4*s3 +
c    &          five*hb5*s4 + six*hb6*s5 + seven*hb7*s6 +
c    &          eight*hb8*s7 + nine*hb9*s8 
c     d1H =   (hden*d1hnum -hnum*d1hden)/hden**two

c
c     calculate first derivative of variables needed with respect to s
c 
c     d1zeta = two*s*H + s2*d1H
c     d1eta  = d1zeta
c     d1lambda = d1zeta
c     d10chi = -f12*nu*d1zeta/(lambda + nu2)**f32
c     d10lambda2 = d10chi*(lambda + nu**two) + (one+chi)*d1lambda
c     !d1lambda2 = (d1lambda*(one-chi)+lambda*d1chi)/(one-chi)**two

c   
c     calculate the first derivative of Fs with respect to s
c   
c     d1Fs = -two*s/(r27*C*(one+s02)**two) - d1zeta/(two*C)

c
c     Calculate the first derivate of EGs with respect to s
c
c     d1EGs = -(two/five)*C*(d1Fs*lambda + Fs*d1lambda) -
c    &        (eight/r15)*B*lambda*d1lambda -
c    &        (r18/five)*A*lambda*lambda*d1lambda -
c    &        (r14/five)*srpi*d1lambda*lambda**f52 -
c    &        (r42/five)*(lambda**f52)*
c    &        d1lambda*(dsqrt(zeta)-dsqrt(eta))-
c    &        (six/five)*(lambda**(seven/two))*
c    &        (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))

c
c     Calculate the first derivate of denominators needed with respect
c     to s
c
c     tmp1 = (dsqrt(zeta+nu2)+nu)/(dsqrt(eta+nu2)) 
c     tmp2 = (dsqrt(eta+nu2)+nu)/(dsqrt(zeta+nu2))

c     d10Js = f12*d1zeta*(two+tmp1+tmp2)
c     d10Ks = d10Js

c     tmp3 = (dsqrt(zeta+nu2)+nu)/(dsqrt(lambda+nu2))
c     tmp4 = (dsqrt(lambda+nu2)+nu)/(dsqrt(zeta+nu2)) 
c     d10Ms = f12*d1zeta*(two +tmp3+tmp4)

c     tmp5 = (dsqrt(lambda+nu2)+nu)/(dsqrt(eta+nu2))
c     tmp6 = (dsqrt(eta+nu2)+nu)/(dsqrt(lambda+nu2))
c     d10Ns = f12*d1zeta*(two + tmp5+tmp6)
c
c
c     Calculate the derivative of the 08-Fxhse with respect to s
c
c     L2 = lambda2*lambda2
c     L3 = lambda2*lambda2*lambda2
c     Js2 = Js*Js
c     Ks2 = Ks*Ks
c     Ms2 = Ms*Ms
c     Ns2 = Ns*Ns
c
c     d1Fxhse1  = A*( (Js*d1zeta - zeta*d10Js)/(Js2) +
c    &                (Ks*d1zeta - eta*d10Ks)/(Ks2) ) 

c     d1Fxhse2  = (four/nine)*B*d10lambda2/L2 

c     tmp9 = d10lambda2/lambda2
c     tmp7 = d1Fs - two*Fs*tmp9
c     tmp8 = one + f12*chi
c     tmp10 =  f12*Fs*d10chi

c     d1Fxhse3 = -(four*C/(nine*L2))*(tmp7*tmp8+tmp10)

c     tmp7 = one + (nine/eight)*chi+(three/eight)*chi*chi
c     tmp8 = (nine/eight)*d10chi + (six/eight)*chi*d10chi

c     d1Fxhse4 = -(eight/(nine*L3))*((d1EGs-three*EGs*tmp9)*tmp7
c    &           + EGs*tmp8)

c     d1Fxhse5  = two*d1zeta*dlog(one-D/Ms) +
c    &           two*zeta*D*d10Ms/(Ms2*(one-D/Ms)) 

c     d1Fxhse6  = -two*d1eta*dlog(one- (D-A)/Ns) -
c    &           two*eta*(D-A)*d10Ns/(Ns2*(one-(D-A)/Ns)) 
c
c     d10Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6
c
c     Calculate the derivative of 08-Fxhse with respect to nu
c
c     nu3 = nu2*nu
c
c     Define all derivatives with respect to nu necessary to 
c     evaluate enhancement factor derivatives.
c
c     d01chi = lambda/(nu2 + lambda)**f32
c     d01lambda2 = (lambda*d01chi)/(one - chi)**two
c     d01Js = (eta*(nu + dsqrt(nu2 + zeta)) + 
c    1        (nu + dsqrt(nu2 + eta))*
c    2        (zeta + two*nu*(nu + dsqrt(nu2 + zeta))))/
c    3        (dsqrt(nu2 + eta)*dsqrt(nu2 + zeta))
c     d01Ks = d01Js 
c     d01Ms = (lambda*(nu + dsqrt(nu2 + zeta)) + 
c    1        (nu + dsqrt(nu2 + lambda))*
c    2        (zeta + two*nu*(nu + dsqrt(nu2 + zeta))))/
c    3        (dsqrt(nu2 + lambda)*dsqrt(nu2 + zeta))
c     d01Ns = (eta*(nu + dsqrt(nu2 + lambda)) + 
c    1        (nu + dsqrt(nu2 + eta))*
c    2        (lambda + two*nu*(nu + dsqrt(nu2 + lambda))))/
c    3        (dsqrt(nu2 + eta)*dsqrt(nu2 + lambda))
c
c     d1Fxhse1 = -((A*(nu + dsqrt(eta + nu2))*zeta)/
c    &            (dsqrt(eta + nu2)*dsqrt(nu2 + zeta)*
c    &            (nu + dsqrt(nu2 + zeta))*
c    &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))))
c
c     d1Fxhse2 = -((A*eta*(nu/dsqrt(eta + nu2) + nu/
c    &            dsqrt(nu2 + zeta)))/
c    &            ((nu + dsqrt(eta + nu2))*
c    &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -
c    &            (A*eta*(one + nu/dsqrt(eta + nu2)))/
c    &            ((nu + dsqrt(eta + nu2))**two*
c    &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta)))
c
c     d1Fxhse3 = (four*B)/(nine*(lambda + nu2)**(f32))
c
c     d1Fxhse4 = (two*C*Fs)/(three*(lambda + nu2)**(f52))
c
c     d1Fxhse5 = (five*EGs*(lambda**two + four*nu3*
c    &            (nu + dsqrt(lambda + nu2)) +
c    &            lambda*nu*(five*nu + three*dsqrt(lambda + nu2))))/
c    &   (three*(lambda + nu2)**four*(nu + dsqrt(lambda + nu2))**three)
c
c     d1Fxhse6 = (two*D*zeta*(nu + dsqrt(nu2 + zeta)))/
c    &            (dsqrt(lambda + nu2)*dsqrt(nu2 + zeta)*
c    &            (-D + lambda + (nu + dsqrt(lambda + nu2))*
c    &            (nu + dsqrt(nu2 + zeta))))
c
c     d1Fxhse7 = (two*(A - D)*eta*(nu + dsqrt(eta + nu2)))/
c    &            (dsqrt(eta + nu2)*dsqrt(lambda + nu2)*
c    &            (A - D + lambda + nu2 + nu*dsqrt(eta + nu2) +
c    &            nu*dsqrt(lambda + nu2) +
c    &            dsqrt(eta + nu2)*dsqrt(lambda + nu2)))
c
c     d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+
c    &           d1Fxhse7
c   
!**********************************************************
! Calculate the Second Derivative of the Enhancement Factor
!**********************************************************
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Calculate the second derivative of H
c
c     d2hnum =  two*ha2+six*ha3*s+r12*ha4*s2+r20*ha5*s3+
c    &           r30*ha6*s4 + r42*ha7*s5

c     d2hden  =  two*hb2+six*hb3*s+r12*hb4*s2+r20*hb5*s3+
c    &           r30*hb6*s4+r42*hb7*s5+r56*hb8*s6+r72*hb9*s7

c     d2H     =  (hden*(hden*d2hnum-hnum*d2hden) -
c    &           two*(hden*d1hnum-hnum*d1hden)*d1hden)/hden**three   

c
c     calculate second derivative of variables needed
c 
c     d2zeta    = two*H +four*s*d1H + s2*d2H
c     d2eta     = d2zeta
c     d2lambda  = d2zeta

c     d20chi    = -f12*(nu/(lambda+nu2)**(f32))*
c    &           (-f32*d1zeta*d1zeta/(lambda+nu2) +d2zeta)

c     d20lambda2 =(one-chi)*(d2lambda-d2lambda*chi+lambda*d20chi)+
c    &           two*d10chi*(d1lambda-d1lambda*chi+lambda*d10chi)
c     d20lambda2 = d20lambda2/(one-chi)**three
c   
c     calculate the second derivative of Fs
c   
c     d2Fs = -(two/(r27*C))*(one-three*s02)/
c    &        ((one+s02)**three)-d2zeta/(two*C)  
c
c     Calculate the second derivative of EGs
c
c     tmp7 = -(two/five)*C*(d2Fs*lambda+two*d1Fs*d1lambda+Fs*d2lambda)
c     tmp8 = -(eight/r15)*B*(d1lambda*d1lambda+lambda*d2lambda)
c     tmp9 =-(r18/five)*A*lambda*(two*d1lambda*d1lambda+lambda*d2lambda)
c     tmp10 = -(r14/five)*srpi*(f52*(lambda**f32)*d1lambda*d1lambda
c    &          +(lambda**f52)*d2lambda)
c     tmp11 = -(r42/five)*(f52*(lambda**f32)*d1lambda*d1lambda
c    &          +(lambda**f52)*d2lambda)*(dsqrt(zeta)-dsqrt(eta))
c     tmp12 = -(r42/five)*(lambda**f52)*d1lambda*
c    &          (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))
c     tmp13 = -(six/five)*(lambda**(seven/two))*
c    &          (-d1zeta*d1zeta/(two*zeta**f32)+d2zeta/dsqrt(zeta)
c    &          +d1eta*d1eta/(two*eta**f32)-d2eta/dsqrt(eta))

c     d2EGs = tmp7+tmp8+tmp9+tmp10+tmp11+tmp12+tmp13
c
c     Calculate the second derivative of denominators needed
c
c     tmp7 = two/(dsqrt(zeta+nu2)*dsqrt(eta+nu2))
c     tmp8 = (dsqrt(eta+nu2)+nu)/(zeta+nu2)**f32
c     tmp9 = (dsqrt(zeta+nu2)+nu)/(eta+nu2)**f32
c     
c     d20Js = f12*(d2zeta*(two+tmp1+tmp2) +
c    &        f12*d1zeta*d1zeta*(tmp7-tmp8-tmp9)) 

c     d20Ks = d20Js

c     tmp10 = two/(dsqrt(zeta+nu2)*dsqrt(lambda+nu2))
c     tmp11 = (dsqrt(lambda+nu2)+nu)/(zeta+nu2)**f32
c     tmp12 = (dsqrt(zeta+nu2)+nu)/(lambda+nu2)**f32
c     d20Ms = f12*(d2zeta*(two+tmp3+tmp4) +
c    &        f12*d1zeta*d1zeta*(tmp10-tmp11-tmp12)) 

c     tmp13 = two/(dsqrt(eta+nu2)*dsqrt(lambda+nu2))
c     tmp14 = (dsqrt(lambda+nu2)+nu)/(eta+nu2)**f32
c     tmp15 = (dsqrt(eta+nu2)+nu)/(lambda+nu2)**f32
c     d20Ns = f12*(d2zeta*(two+tmp5+tmp6) +
c    &        f12*d1zeta*d1zeta*(tmp13-tmp14-tmp15)) 
c
c     Calculate the second derivative of the Fxhse with respect to the
c     reduced density gradient, s.
c
c     Ms4 = Ms2*Ms2
c     Ns4 = Ns2*Ns2
c
c     tmp1      = A*(Js2*d2zeta -two*Js*d10Js*d1zeta 
c    &            + two*d10Js*d10Js*zeta-Js*d20Js*zeta)/
c    &            (Js2*Js) 
c     tmp2       = A*(Ks2*d2zeta -two*Ks*d10Ks*d1zeta 
c    &            + two*d10Ks*d10Ks*eta-Ks*d20Ks*eta)/
c    &            (Ks2*Ks) 
c     d2Fxhse1  = tmp1 +tmp2 

c     d2Fxhse2  = (four/nine)*B*(L2*d20lambda2-
c    &            two*lambda2*d10lambda2*d10lambda2)/(L2*L2)

c     tmp4 = d10lambda2/lambda2
c     tmp5 = d20lambda2/lambda2
c     d2Fxhse3  = -((four*C)/(nine*L2))*(d2Fs +
c    &           six*Fs*tmp4**two - two*Fs*tmp5 -
c    &           four*d1Fs*tmp4 +
c    &          f12*((d2Fs*chi+two*d1Fs*d10chi+Fs*d20chi)-
c    &           four*tmp4*(d1Fs*chi+Fs*d10chi) +
c    &           six*chi*Fs*tmp4*tmp4)-chi*Fs*tmp5)

c     tmp1 = -EGs*(eight/nine)/L3
c     tmp2 = (-d1EGs + three*EGs*d10lambda2/lambda2)*(eight/nine)/L3
c     tmp3 = (-d2EGs - r12*EGs*d10lambda2**two/L2 + 
c    &   three*(two*d1EGs*d10lambda2 + EGs*d20lambda2)/lambda2)
c    &   *(eight/nine)/L3
c     tmp4 = (one + (nine/eight)*chi + (three/eight)*chi**two)
c     tmp5 = (three*(three + two*chi)*d10chi)/eight
c     tmp6 = (three*(two*d10chi**two + (three + two*chi)*d20chi))/eight

c     d2Fxhse4 = (tmp1*tmp6+two*tmp2*tmp5+tmp3*tmp4)

c     tmp1      = (one-D/Ms) 

c     d2Fxhse5  = two*d2zeta*dlog(tmp1)+four*d1zeta*D*d10Ms/(Ms2*tmp1)
c    &          - two*zeta*(D*d10Ms/(Ms2*tmp1))**two +
c    &          two*zeta*D*(Ms2*d20Ms-two*Ms*d10Ms*d10Ms)/
c    &          (tmp1*Ms**four)

c     tmp1      = (one-(D-A)/Ns) 
c     d2Fxhse6  = -two*d2eta*dlog(tmp1)-
c    &             four*d1eta*(D-A)*d10Ns/(Ns2*tmp1)+
c    &             two*eta*((D-A)*d10Ns/(Ns2*tmp1))**two -
c    &          two*eta*(D-A)*(Ns2*d20Ns-two*Ns*d10Ns*d10Ns)/
c    &          (tmp1*Ns**four)

c     d20Fxhse = d2Fxhse1+d2Fxhse2+d2Fxhse3+d2Fxhse4+d2Fxhse5+d2Fxhse6
c
c     Calculate the second derivative of Fxhse with respect to the
c     parameter nu (cam_omega/kf).
c
c     nu5 = nu3*nu2
c     nu6 = nu5*nu
c
c     Calculate second derivatives with respect to nu necessary to
c     simplify the appearance of the code.
c
c     d02chi = (-three*nu*lambda)/(nu2 + lambda)**f52
c     d02lambda2 = -((lambda*(two*d01chi**two - 
c    &             (-one + chi)*d02chi))/
c    &             (-one + chi)**three)
c     d02Js = two*nu*(one/dsqrt(nu2 + eta) + one/dsqrt(nu2 + zeta))*
c    &        (one + nu/dsqrt(nu2 + zeta)) + 
c    &        (zeta*(dsqrt(nu2 + eta) + 
c    &        dsqrt(nu2 + zeta)))/(nu2 + zeta)**f32 + 
c    &        (nu + dsqrt(nu2 + zeta))*
c    &        (one/dsqrt(nu2 + eta) + one/dsqrt(nu2 + zeta) + 
c    &        nu2*(-(nu2 + eta)**(-f32) - 
c    &        (nu2 + zeta)**(-f32)))
c     d02Ks = d02Js
c     d02Ms = two*nu*(one + nu/dsqrt(nu2 + lambda))*
c    &        (one/dsqrt(nu2 + lambda) + one/dsqrt(nu2 + zeta))
c    &        + (lambda*(dsqrt(nu2 + lambda) +
c    &        dsqrt(nu2 + zeta)))/(nu2 + lambda)**f32 +
c    &        (nu + dsqrt(nu2 + lambda))*
c    &        (one/dsqrt(nu2 + lambda) + one/dsqrt(nu2 + zeta) +
c    &        nu2*(-(nu2 + lambda)**(-f32) -
c    &        (nu2 + zeta)**(-f32)))
c     d02Ns = two*nu*(one/dsqrt(nu2 + eta) + 
c    &        one/dsqrt(nu2 + lambda))*
c    &        (one + nu/dsqrt(nu2 + lambda)) + 
c    &        (lambda*(dsqrt(nu2 + eta) + 
c    &        dsqrt(nu2 + lambda)))/(nu2 + lambda)**f32
c    &        + (nu + dsqrt(nu2 + lambda))*
c    &        (one/dsqrt(nu2 + eta) + one/dsqrt(nu2 + lambda) + 
c    &        nu2*(-(nu2 + eta)**(-f32) - 
c    &        (nu2 + lambda)**(-f32))) 
c
c     L5 = L2*L3
c     Js3 = Js2*Js
c     Ks3 = Ks2*Ks
c
c     d2Fxhse1 = A*((zeta*(two*d01Js**two - Js*d02Js))/Js3 + 
c    &    (eta*(two*d01Ks**two - Ks*d02Ks))/Ks3)
c
c     d2Fxhse2 = (four*B*(-two*d01lambda2**two + lambda2*d02lambda2))/
c    &    (nine*L3)
c
c     d2Fxhse3 = -(ten*C*Fs*nu)/(three*(lambda+nu2)**f72)
c
c     d2Fxhse4 = (EGs*(-four*
c    &    (eight + nine*chi + three*chi**two)*d01lambda2**two - 
c    &    lambda2**two*(two*d01chi**two + (three + two*chi)*d02chi) + 
c    &    lambda2*(six*(three + two*chi)*d01chi*d01lambda2 + 
c    &    (eight + nine*chi + three*chi**two)*d02lambda2)))/
c    &    (three*L5)
c
c     d2Fxhse5 = (two*D*zeta*((D - two*Ms)*d01Ms**two + 
c    &    Ms*(-D + Ms)*d02Ms))/((D - Ms)**two*Ms**two)
c
c     d2Fxhse6 = (two*(A - D)*eta*((-A + D - two*Ns)*d01Ns**two +
c    &  Ns*(A - D + Ns)*d02Ns))/
c    &  (Ns**two*(A - D + Ns)**two)
c
c     d02Fxhse = d2Fxhse1 + d2Fxhse2 + d2Fxhse3 + d2Fxhse4 + 
c    &           d2Fxhse5 + d2Fxhse6
c
c     Calculate the mixed partial derivative of the enhancement factor 
c     with respect to both the reduced density gradient, s, and the 
c     parameter nu (cam_omega/kF).  Note that the enhancement factor
c     satisfies continuity requirements and therefore it does not matter
c     what order the derivatives are taken in.
c     
c     nu4 = nu3*nu
c
c     Calculate mixed partial derivatives for the variables depending
c     on both s and nu.
c
c     d11chi = (three*nu2*d1lambda)/
c    &         (two*(nu2 + lambda)**f52) - 
c    &         d1lambda/(two*(nu2 + lambda)**f32) 
c     d11lambda2 = (d1lambda*d01chi)/
c    &             (one - chi)**two + 
c    &             (two*lambda*d01chi*
c    &             d10chi)/(one - chi)**three + 
c    &             (lambda*d11chi)/
c    &             (one - chi)**two
c     d11Js = (-(nu*eta**two*d1zeta) + 
c    &        nu*zeta*((-nu2 - zeta)*
c    &        d1eta + 
c    &        nu*(nu + dsqrt(nu2 + eta))*
c    &        d1zeta) + 
c    &        eta*((nu2 + zeta)*
c    &        (nu + dsqrt(nu2 + zeta))*d1eta
c    &        + (-nu3 + (nu + dsqrt(nu2 + eta))*zeta)*
c    &        d1zeta))/
c    &        (two*(nu2 + eta)**f32*(nu2 + zeta)**f32)
c     d11Ks = d11Js
c     d11Ms = (-(nu*lambda**two*d1zeta) + 
c    &        nu*zeta*((-nu2 - zeta)*
c    &        d1lambda + 
c    &        nu*(nu + dsqrt(nu2 + lambda))*
c    &        d1zeta) + 
c    &        lambda*((nu2 + zeta)*
c    &        (nu + dsqrt(nu2 + zeta))*
c    &        d1lambda + 
c    &        (-nu3 + (nu + dsqrt(nu2 + lambda))*zeta)*
c    &        d1zeta))/
c    &        (two*(nu2 + lambda)**f32*(nu2 + zeta)**f32)
c     d11Ns = (-(nu*eta**two*d1lambda) + 
c    &        nu*lambda*((-nu2 - lambda)*
c    &        d1eta + 
c    &        nu*(nu + dsqrt(nu2 + eta))*
c    &        d1lambda) + 
c    &        eta*((nu2 + lambda)*
c    &        (nu + dsqrt(nu2 + lambda))*
c    &        d1eta + 
c    &        (-nu3 + (nu + dsqrt(nu2 + eta))*lambda)*
c    &        d1lambda))/
c    &        (two*(nu2 + eta)**f32*(nu2 + lambda)**f32)
c
c     d11Fxhse1 = -((A*(-two*Ks3*zeta*d01Js*d10Js + 
c    &   Js*Ks3*(d1zeta*d01Js + zeta*d11Js) + 
c    &   Js3*(-two*eta*d01Ks*d10Ks + 
c    &   Ks*(d1eta*d01Ks + eta*d11Ks))))/(Js3*Ks3))
c
c     d11Fxhse2 = (-two*B*(two*nu2*(nu + dsqrt(nu2 + lambda)) +
c    &   lambda*(two*nu + dsqrt(nu2 + lambda)))*d1lambda)/
c    &   (three*(nu2 + lambda)**three*
c    &   (nu + dsqrt(nu2 + lambda))**two)
c
c     d11Fxhse3 = (C*(two*(nu2 + lambda)*d1Fs - five*Fs*d1lambda))/
c    &   (three*(nu2 + lambda)**(f72))
c
c     d11Fxhse4 = (five*(eight*nu4*(nu + dsqrt(nu2 + lambda)) +
c    &   lambda**two*(four*nu + dsqrt(nu2 + lambda)) +
c    &   four*nu2*lambda*(three*nu + two*dsqrt(nu2 + lambda)))*
c    &   (two*(nu2 + lambda)*d1EGs - seven*EGs*d1lambda))/
c    &   (six*(nu2 + lambda)**five*(nu + dsqrt(nu2 + lambda))**four)
c
c     d11Fxhse5 = (two*D*(D*zeta*d01Ms*d10Ms + 
c    &   Ms**two*(d1zeta*d01Ms + zeta*d11Ms) - 
c    &   Ms*(D*d1zeta*d01Ms + 
c    &   zeta*(two*d01Ms*d10Ms + D*d11Ms))))/((D - Ms)**two*Ms**two)
c
c     d11Fxhse6 = (two*(A - D)*((-A + D)*eta*d01Ns*d10Ns + 
c    &   Ns**two*(d1eta*d01Ns + eta*d11Ns) + 
c    &   Ns*((A - D)*d1eta*d01Ns + 
c    &   eta*(-two*d01Ns*d10Ns + (A - D)*d11Ns))))/
c    &   (Ns**two*(A - D + Ns)**two)
c
c     d11Fxhse = d11Fxhse1 + d11Fxhse2 + d11Fxhse3 + d11Fxhse4 + 
c    &           d11Fxhse5 + d11Fxhse6
#endif
!*********************************************************
! Calculate the Third Derivative of the Enhancement Factor
!*********************************************************
#ifdef THIRD_DERIV
c
c    Calculate the third order derivative of H with respect to s.
c
c     d3hnum = six*ha3+r24*ha4*s+r60*ha5*s2+
c    1         r120*ha6*s3 + r210*ha7*s4

c     d3hden  =  six*hb3+r24*hb4*s+r60*hb5*s2+
c    1           r120*hb6*s3+r210*hb7*s4+r336*hb8*s5+r504*hb9*s6

c     d3H = (-six*hnum*d1hden**three +
c    1  six*hden*d1hden*(d1hnum*d1hden + hnum*d2hden) +
c    2  hden**three*d3hnum - hden**two*(three*d1hden*d2hnum +
c    3  three*d1hnum*d2hden + hnum*d3hden))/hden**four
c
c    Calculate the third order derivatives of the "helper variables"
c
c     d3zeta = six*d1H + six*s*d2H + s2*d3H
c     d3eta = d3zeta
c     d3lambda = d3zeta
c     d30chi = -(nu*(r15*d1lambda*d1lambda*d1lambda
c    1    - r18*(nu2 + lambda)*d1lambda*
c    2   d2lambda + four*(nu2 + lambda)**two*d3lambda))/
c    3   (eight*(nu2 + lambda)**f72) 
c     d30lambda2 = d3lambda/(one - chi) + 
c    1             (three*d2lambda*d10chi)/
c    2             (-one + chi)**two - 
c    3             (three*d1lambda*
c    4             (two*d10chi**two - 
c    5             (-one + chi)*d20chi))/
c    6             (-one + chi)**three + 
c    7             (lambda*(six*d10chi**three - 
c    8             six*(-one + chi)*d10chi*
c    9             d20chi + 
c    A             (-one + chi)**two*d30chi))/
c    B             (-one + chi)**four
c
c    Calculate the third order derivative of Fs
c
c     d3Fs = -(r256*s*(s2 - four) +
c    1       nine*(s2 + four)**four*d3zeta)/
c    2       (r18*C*(s2 + four)**four)

c
c    Calculate the third order derivative of EGs
c 
c     tmp1 = (-two*C*(three*d1lambda*d2Fs + 
c    1  three*d1Fs*d2lambda + lambda*d3Fs + 
c    2  Fs*d3lambda))/five
c     tmp2 = (-four*B*(six*d1lambda*d2lambda + 
c    1  two*lambda*d3lambda))/r15
c     tmp3 = (-six*A*(six*d1lambda**three + 
c    1  r18*lambda*d1lambda*d2lambda + 
c    2  three*lambda**two*d3lambda))/five
c     tmp4 = (-four*srpi*((r105*dsqrt(lambda)*
c    1  d1lambda**three)/eight + (r105*lambda**f32*d1lambda*
c    2  d2lambda)/four + (seven*lambda**f52*d3lambda)/two))/five
c     tmp5 = (-r36*(-d1eta/(two*dsqrt(eta)) + 
c    1  d1zeta/(two*dsqrt(zeta)))*
c    2  ((r35*lambda**f32*d1lambda**two)/four +
c    3  (seven*lambda**f52*d2lambda)/two))/five 
c    4  - (r126*lambda**f52*d1lambda*(d1eta**two/(four*eta**f32) - 
c    5  d1zeta**two/(four*zeta**f32) - d2eta/(two*dsqrt(eta)) + 
c    6  d2zeta/(two*dsqrt(zeta))))/five - 
c    7  (r12*(-dsqrt(eta) + dsqrt(zeta))*
c    8  ((r105*dsqrt(lambda)*d1lambda**three)/
c    9  eight + (r105*lambda**f32*d1lambda*d2lambda)/four + 
c    A  (seven*lambda**f52*d3lambda)/two))/five
c    B  - (r12*lambda**f72*((-three*d1eta**three)/(eight*eta**f52) + 
c    C  (three*d1zeta**three)/(eight*zeta**f52) + 
c    D  (three*d1eta*d2eta)/(four*eta**f32) - 
c    E  (three*d1zeta*d2zeta)/(four*zeta**f32) - 
c    F  d3eta/(two*dsqrt(eta)) + d3zeta/(two*dsqrt(zeta))))/five

c     d3EGs = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 

c
c    Calculate derivatives of denominators needed for third derivatives
c
c     d30Js = (three*(nu + dsqrt(nu2 + zeta))*
c    1        d1eta**three - 
c    2        (three*(nu2 + eta)*d1eta**two*
c    3        d1zeta)/dsqrt(nu2 + zeta) - 
c    4        (three*(nu2 + eta)*d1eta*
c    5        ((nu2 + eta)*d1zeta**two + 
c    6        two*(nu2 + zeta)*
c    7        ((zeta + nu*(nu + dsqrt(nu2 + zeta)))*
c    8        d2eta - 
c    9        (nu2 + eta)*d2zeta)))/
c    A        (nu2 + zeta)**f32 + 
c    B        ((nu2 + eta)**two*
c    C        (three*(eta + nu*(nu + dsqrt(nu2 + eta)))*
c    D        d1zeta**three + 
c    E        six*(nu2 + zeta)*d1zeta*
c    F        ((nu2 + zeta)*d2eta - 
c    G        (eta + nu*(nu + dsqrt(nu2 + eta)))*
c    H        d2zeta) + 
c    I        four*(nu2 + zeta)**two*
c    J        ((zeta + nu*(nu + dsqrt(nu2 + zeta)))*
c    K        d3eta + 
c    L        (nu2 + eta + nu*dsqrt(nu2 + eta) + 
c    M        two*dsqrt(nu2 + eta)*
c    N        dsqrt(nu2 + zeta))*
c    O        d3zeta)))/
c    P        (nu2 + zeta)**f52)/(eight*(nu2 + eta)**f52) 
c     d30Ks = d30Js
c     d30Ms = (three*(nu + dsqrt(nu2 + zeta))*
c    1        d1lambda**three - 
c    2        (three*(nu2 + lambda)*d1lambda**two*
c    3        d1zeta)/dsqrt(nu2 + zeta) - 
c    4        (three*(nu2 + lambda)*d1lambda*
c    5        ((nu2 + lambda)*d1zeta**two + 
c    6        two*(nu2 + zeta)*
c    7        ((zeta + nu*(nu + dsqrt(nu2 + zeta)))*
c    8        d2lambda - 
c    9        (nu2 + lambda)*d2zeta)))/
c    A        (nu2 + zeta)**f32 + 
c    B        ((nu2 + lambda)**two*
c    C        (three*(lambda + nu*(nu + dsqrt(nu2 + lambda)))*
c    D        d1zeta**three + 
c    E        six*(nu2 + zeta)*d1zeta*
c    F        ((nu2 + zeta)*d2lambda - 
c    G        (lambda + 
c    H        nu*(nu + dsqrt(nu2 + lambda)))*
c    I        d2zeta) + 
c    J        four*(nu2 + zeta)**two*
c    K        ((nu2 + zeta + nu*dsqrt(nu2 + zeta) + 
c    L        two*dsqrt(nu2 + lambda)*
c    M        dsqrt(nu2 + zeta))*
c    N        d3lambda + 
c    O        (lambda + 
c    P        nu*(nu + dsqrt(nu2 + lambda)))*
c    Q        d3zeta)))/
c    R        (nu2 + zeta)**f52)/(eight*(nu2 + lambda)**f52)
c     d30Ns = (three*(nu + dsqrt(nu2 + lambda))*
c    1        d1eta**three - 
c    2        (three*(nu2 + eta)*d1eta**two*
c    3        d1lambda)/dsqrt(nu2 + lambda)
c    4        - (three*(nu2 + eta)*d1eta*
c    5        ((nu2 + eta)*d1lambda**two + 
c    6        two*(nu2 + lambda)*
c    7        ((lambda + 
c    8        nu*(nu + dsqrt(nu2 + lambda)))*
c    9        d2eta - 
c    A        (nu2 + eta)*d2lambda)))/
c    B        (nu2 + lambda)**f32 + 
c    C        ((nu2 + eta)**two*
c    D        (three*(eta + nu*(nu + dsqrt(nu2 + eta)))*
c    E        d1lambda**three + 
c    F        six*(nu2 + lambda)*d1lambda*
c    G        ((nu2 + lambda)*d2eta - 
c    H        (eta + nu*(nu + dsqrt(nu2 + eta)))*
c    I        d2lambda) + 
c    J        four*(nu2 + lambda)**two*
c    K        ((lambda + 
c    L        nu*(nu + dsqrt(nu2 + lambda)))*
c    M        d3eta + 
c    N        (nu2 + eta + nu*dsqrt(nu2 + eta) + 
c    O        two*dsqrt(nu2 + eta)*
c    P        dsqrt(nu2 + lambda))*
c    Q        d3lambda)))/
c    R        (nu2 + lambda)**f52)/(eight*(nu2 + eta)**f52)
c
c     Calculate derivatives of the enhancement factor with respect
c     to s.
c
c     L4 = L2*L2
c     L6 = L4*L2
c     Ms3 = Ms2*Ms
c     Ns3 = Ns2*Ns
c     Js4 = Js2*Js2
c     Ks4 = Ks2*Ks2
c
c     d3Fxhse1 = (A*(Js3*Ks4*d3zeta - 
c    1    six*Ks4*zeta*d10Js**three + 
c    2    six*Js*Ks4*d10Js*
c    3    (d1zeta*d10Js + zeta*d20Js) - 
c    4    Js**two*Ks4*(three*d2zeta*d10Js + 
c    5    three*d1zeta*d20Js + 
c    6    zeta*d30Js) + Js4*
c    7    (Ks3*d3eta - six*eta*d10Ks**three + 
c    8    six*Ks*d10Ks*(d1eta*d10Ks + eta*d20Ks) - 
c    9    Ks**two*(three*d2eta*d10Ks + three*d1eta*d20Ks + 
c    A    eta*d30Ks))))/(Js4*Ks4)
c
c     d3Fxhse2 = (four*B*(six*d10lambda2**three - 
c    1    six*lambda2*d10lambda2*d20lambda2 + 
c    2    lambda2**two*d30lambda2))/(nine*L4)
c
c     d3Fxhse3 = (two*C*(r24*(two + chi)*Fs*d10lambda2**three - 
c    1   r18*lambda2*d10lambda2*
c    2   ((two + chi)*d1Fs*d10lambda2 + 
c    3   Fs*(d10chi*d10lambda2 + (two +
c    4   chi)*d20lambda2)) - 
c    5   L3*((two + chi)*d3Fs + three*d2Fs*d10chi + 
c    6   three*d1Fs*d20chi + Fs*d30chi) + 
c    7   two*L2*(three*(two + chi)*d2Fs*d10lambda2 + 
c    8   three*d1Fs*(two*d10chi*d10lambda2 + 
c    9   (two + chi)*d20lambda2) + 
c    A   Fs*(three*d10lambda2*d20chi + three*d10chi*d20lambda2 + 
c    B   (two + chi)*d30lambda2))))/(nine*L5)
c
c     tmp1 = -EGs*(eight/nine)/L3
c     tmp2 = (-d1EGs + (three*EGs*d10lambda2)/(lambda2))*(eight/nine)/L3
c     tmp3 = (-d2EGs + (six*d1EGs*d10lambda2)/(lambda2)
c    1  - (EGs*((r12*d10lambda2**two)/L2 - 
c    2  (three*d20lambda2)/lambda2)))*(eight/nine)/L3
c     tmp4 = (-d3EGs + (nine*d2EGs*d10lambda2)/lambda2 - 
c    1  three*(d1EGs*((r12*d10lambda2**two)/L2 - 
c    2  (three*d20lambda2)/lambda2)) - (EGs*
c    3  ((-r60*d10lambda2**three)/L3 + 
c    4  (r36*d10lambda2*d20lambda2)/L2
c    5  - (three*d30lambda2)/lambda2)))*(eight/nine)/L3
c     tmp5 = (one + (nine/eight)*chi + (three/eight)*chi*chi)
c     tmp6 = (nine*d10chi)/eight + (three*chi*d10chi)/four
c     tmp7 = (nine*d20chi)/eight + (three*(two*d10chi**two + 
c    1  two*chi*d20chi))/eight
c     tmp8 = (nine*d30chi)/eight + (three*(six*d10chi*
c    1  d20chi + two*chi*d30chi))/eight


c     d3Fxhse4 = (tmp1*tmp8+3d0*tmp2*tmp7+3d0*tmp3*tmp6+
c    1           tmp4*tmp5)
c
c     d3Fxhse5 = two*dlog(one - D/Ms)*d3zeta - 
c    1    (six*D*d2zeta*d10Ms)/(D*Ms - Ms**two) + 
c    2  (six*D*d1zeta*((D - two*Ms)*d10Ms**two + Ms*(-D + Ms)*d20Ms))/
c    3  ((D - Ms)**two*Ms**two) - 
c    4  (two*D*zeta*
c    5  (two*(D**two - three*D*Ms + three*Ms**two)*d10Ms**three - 
c    6  three*Ms*(D**two - three*D*Ms + two*Ms**two)*d10Ms*d20Ms + 
c    7  (D - Ms)**two*Ms**two*d30Ms))/((D - Ms)**three*Ms3) 
c
c     d3Fxhse6 = -two*dlog((A - D + Ns)/Ns)*d3eta + 
c    1  (six*(A - D)*d2eta*d10Ns)/(Ns*(A - D + Ns)) + 
c    2  (six*(A - D)*d1eta*((-A + D - two*Ns)*d10Ns**two + 
c    3  Ns*(A - D + Ns)*d20Ns))/(Ns**two*(A - D + Ns)**two) + 
c    4  (two*(A - D)*eta*(two*((A - D)**two + three*(A - D)*Ns
c    5  + three*Ns**two)*d10Ns**three - 
c    6  three*Ns*((A - D)**two + three*(A - D)*Ns +
c    7  two*Ns**two)*d10Ns*d20Ns + 
c    8  Ns**two*(A - D + Ns)**two*d30Ns))/
c    9  (Ns3*(A - D + Ns)**three)
c
c     d30Fxhse = d3Fxhse1 + d3Fxhse2 + d3Fxhse3 + d3Fxhse4 +
c    1           d3Fxhse5 + d3Fxhse6      
c
c
c     Calculate the third derivative of Fxhse with respect to the
c     parameter nu (cam_omega/kf).
c
c     nu7 = nu6*nu
c     nu8 = nu7*nu
c
c     Calculate third derivatives with respect to nu necessary to
c     simplify the appearance of the code.
c
c     d03chi = (three*(four*nu2 - lambda)*lambda)/(nu2 + lambda)**f72
c     d03lambda2 = (lambda*(six*d01chi**three - 
c    1             six*(-one + chi)*d01chi*
c    2             d02chi + 
c    3             (-one + chi)**two*d03chi))/
c    4             (-one + chi)**four
c     d03Js = (three*(-(nu*eta**three*zeta) + 
c    1        nu4*(nu + dsqrt(nu2 + eta))*zeta**two - 
c    2        nu*eta*zeta*
c    3        (two*nu4 - two*nu*dsqrt(nu2 + eta)*zeta + 
c    4        zeta**two) + 
c    5        eta**two*(two*nu2*zeta*dsqrt(nu2 + zeta) + 
c    6        nu4*(nu + dsqrt(nu2 + zeta)) + 
c    7        zeta**two*
c    8        (two*nu + dsqrt(nu2 + eta) + 
c    9        dsqrt(nu2 + zeta)))))/
c    A        ((nu2 + eta)**f52*(nu2 + zeta)**f52)
c     d03Ks = d03Js
c     d03Ms = (three*(-(nu*lambda**three*zeta) + 
c    1        nu4*(nu + dsqrt(nu2 + lambda))*zeta**two - 
c    2        nu*lambda*zeta*
c    3        (two*nu4 - two*nu*dsqrt(nu2 + lambda)*zeta + 
c    4        zeta**two) + 
c    5        lambda**two*(two*nu2*zeta*
c    6        dsqrt(nu2 + zeta) + 
c    7        nu4*(nu + dsqrt(nu2 + zeta)) + 
c    8        zeta**two*
c    9        (two*nu + dsqrt(nu2 + lambda) + 
c    A        dsqrt(nu2 + zeta)))))/
c    B        ((nu2 + lambda)**f52*(nu2 + zeta)**f52)
c     d03Ns = (three*(-(nu*eta**three*lambda) + 
c    1        nu4*(nu + dsqrt(nu2 + eta))*lambda**two - 
c    2        nu*eta*lambda*
c    3        (two*nu4 - two*nu*dsqrt(nu2 + eta)*lambda + 
c    4        lambda**two) + 
c    5        eta**two*(two*nu2*lambda*
c    6        dsqrt(nu2 + lambda) + 
c    7        nu4*(nu + dsqrt(nu2 + lambda)) + 
c    8        lambda**two*
c    9        (two*nu + dsqrt(nu2 + eta) + 
c    A        dsqrt(nu2 + lambda)))))/
c    B        ((nu2 + eta)**f52*(nu2 + lambda)**f52)
c
c     d3Fxhse1 = A*(-((zeta*(six*d01Js**three - six*Js*d01Js*d02Js + 
c    1           Js**two*d03Js))/Js4) - 
c    2    (eta*(six*d01Ks**three - six*Ks*d01Ks*d02Ks + 
c    3         Ks**two*d03Ks))/Ks4)
c
c     d3Fxhse2 = (four*B*(six*d01lambda2**three - 
c    1      six*lambda2*d01lambda2*d02lambda2 + 
c    2      lambda2**two*d03lambda2))/(nine*L4)
c
c     d3Fxhse3 = (two*C*Fs*(r24*(two + chi)*
c    1  d01lambda2**three - r18*lambda2*d01lambda2*
c    2  (d01chi*d01lambda2 + (two + chi)*d02lambda2)
c    3  - L3*d03chi + two*L2*(three*d01lambda2*d02chi + 
c    4  three*d01chi*d02lambda2 + 
c    5  (two + chi)*d03lambda2)))/(nine*L5)
c
c     d3Fxhse4 = (EGs*(r20*(eight + nine*chi + three*chi**two)*
c    1  d01lambda2**three - r12*lambda2*d01lambda2*
c    2  (three*(three + two*chi)*d01chi*d01lambda2 + 
c    3  (eight + nine*chi + three*chi**two)*d02lambda2) - 
c    4  L3*(six*d01chi*d02chi + (three + two*chi)*d03chi) + 
c    5  L2*(r18*d01chi**two*d01lambda2 + nine*(three + two*chi)*
c    6  d01lambda2*d02chi + nine*(three + two*chi)*d01chi*
c    7  d02lambda2 + (eight + nine*chi + three*chi**two)*
c    8  d03lambda2)))/(three*L6)
c
c     d3Fxhse5 = (-two*D*zeta*(two*(D**two - three*D*Ms + 
c    1  three*Ms**two)*d01Ms**three - 
c    2  three*Ms*(D**two - three*D*Ms + two*Ms**two)*d01Ms*d02Ms + 
c    3  (D - Ms)**two*Ms**two*d03Ms))/((D - Ms)**three*Ms3)
c
c     d3Fxhse6 = (two*(A - D)*eta*(two*((A - D)**two + 
c    1  three*(A - D)*Ns + three*Ns**two)*d01Ns**three - 
c    2  three*Ns*((A - D)**two + 
c    3  three*(A - D)*Ns + two*Ns**two)*d01Ns*d02Ns + 
c    4  Ns**two*(A - D + Ns)**two*d03Ns))/
c    5  (Ns3*(A - D + Ns)**three)
c
c     d03Fxhse = d3Fxhse1 + d3Fxhse2 + d3Fxhse3 + d3Fxhse4 +
c    1           d3Fxhse5 + d3Fxhse6
c
c     Calculate the mixed third derivative of Fxhse (dnu ds^2)
c
c
c     Calculate mixed partial derivatives for the variables depending
c     on both s and nu.
c
c     d21chi = (three*(-four*nu2 + lambda)*d1lambda**two + 
c    1         two*(two*nu4 + nu2*lambda - lambda**two)*d2lambda)/
c    2         (four*(nu2 + lambda)**f72)
c     d21lambda2 = ((-one + chi)**two*d2lambda*d01chi + 
c    1         two*(-one + chi)*d1lambda*(-two*d01chi*d10chi + 
c    2         (-one + chi)*d11chi) + 
c    3         lambda*(two*d01chi*(three*d10chi**two - 
c    4         (-one + chi)*d20chi) + 
c    5    (-one + chi)*(-four*d10chi*d11chi + (-one + chi)*d21chi)))/
c    6         (-one + chi)**four
c     d21Js = (-(nu*(nu2 + eta)*(nu2 + zeta)**two*d1eta**two) + 
c    1        three*nu2*(nu2 + zeta)**f52*d1eta**two - 
c    2        (nu2 + eta)*(nu2 + zeta)**f52*d1eta**two +
c    3        three*nu*(nu2 + zeta)**three*d1eta**two - 
c    4        two*nu*(nu2 + eta)**two*(nu2 + zeta)*d1eta*d1zeta - 
c    5        two*nu*(nu2 + eta)*(nu2 + zeta)**two*d1eta*d1zeta + 
c    6        three*nu2*(nu2 + eta)**f52*d1zeta**two + 
c    7        three*nu*(nu2 + eta)**three*d1zeta**two - 
c    8        nu*(nu2 + eta)**two*(nu2 + zeta)*d1zeta**two - 
c    9        (nu2 + eta)**f52*(nu2 + zeta)*d1zeta**two + 
c    A        two*nu*(nu2 + eta)**two*(nu2 + zeta)**two*d2eta - 
c    B        two*nu2*(nu2 + eta)*(nu2 + zeta)**f52*d2eta + 
c    C        two*(nu2 + eta)**two*(nu2 + zeta)**f52*d2eta - 
c    D        two*nu*(nu2 + eta)*(nu2 + zeta)**three*d2eta - 
c    E        two*nu2*(nu2 + eta)**f52*(nu2 + zeta)*d2zeta - 
c    F        two*nu*(nu2 + eta)**three*(nu2 + zeta)*d2zeta + 
c    G        two*nu*(nu2 + eta)**two*(nu2 + zeta)**two*d2zeta + 
c    H        two*(nu2 + eta)**f52*(nu2 + zeta)**two*d2zeta)/
c    I        (four*(nu2 + eta)**f52*(nu2 + zeta)**f52)
c     d21Ks = d21Js
c     d21Ms = -(nu*(nu2 + lambda)*(nu2 + zeta)**two*d1lambda**two - 
c    1        three*nu2*(nu2 + zeta)**f52*d1lambda**two + 
c    2        (nu2 + lambda)*(nu2 + zeta)**f52*d1lambda**two -
c    3        three*nu*(nu2 + zeta)**three*d1lambda**two + 
c    4      two*nu*(nu2 + lambda)**two*(nu2 + zeta)*d1lambda*d1zeta + 
c    5      two*nu*(nu2 + lambda)*(nu2 + zeta)**two*d1lambda*d1zeta - 
c    6        three*nu2*(nu2 + lambda)**f52*d1zeta**two -
c    7        three*nu*(nu2 + lambda)**three*d1zeta**two + 
c    8        nu*(nu2 + lambda)**two*(nu2 + zeta)*d1zeta**two + 
c    9        (nu2 + lambda)**f52*(nu2 + zeta)*d1zeta**two - 
c    A        two*nu*(nu2 + lambda)**two*(nu2 + zeta)**two*d2lambda + 
c    B        two*nu2*(nu2 + lambda)*(nu2 + zeta)**f52*d2lambda - 
c    C        two*(nu2 + lambda)**two*(nu2 + zeta)**f52*d2lambda + 
c    D        two*nu*(nu2 + lambda)*(nu2 + zeta)**three*d2lambda + 
c    E        two*nu2*(nu2 + lambda)**f52*(nu2 + zeta)*d2zeta + 
c    F        two*nu*(nu2 + lambda)**three*(nu2 + zeta)*d2zeta - 
c    G        two*nu*(nu2 + lambda)**two*(nu2 + zeta)**two*d2zeta - 
c    H        two*(nu2 + lambda)**f52*(nu2 + zeta)**two*d2zeta)/
c    I        (four*(nu2 + lambda)**f52*(nu2 + zeta)**f52) 
c     d21Ns = -(nu*(nu2 + eta)*(nu2 + lambda)**two*d1eta**two -
c    1        three*nu2*(nu2 + lambda)**f52*d1eta**two + 
c    2        (nu2 + eta)*(nu2 + lambda)**f52*d1eta**two -
c    3        three*nu*(nu2 + lambda)**three*d1eta**two + 
c    4        two*nu*(nu2 + eta)**two*(nu2 + lambda)*d1eta*d1lambda + 
c    5        two*nu*(nu2 + eta)*(nu2 + lambda)**two*d1eta*d1lambda - 
c    6        three*nu2*(nu2 + eta)**f52*d1lambda**two -
c    7        three*nu*(nu2 + eta)**three*d1lambda**two + 
c    8        nu*(nu2 + eta)**two*(nu2 + lambda)*d1lambda**two + 
c    9        (nu2 + eta)**f52*(nu2 + lambda)*d1lambda**two - 
c    A        two*nu*(nu2 + eta)**two*(nu2 + lambda)**two*d2eta + 
c    B        two*nu2*(nu2 + eta)*(nu2 + lambda)**f52*d2eta - 
c    C        two*(nu2 + eta)**two*(nu2 + lambda)**f52*d2eta + 
c    D        two*nu*(nu2 + eta)*(nu2 + lambda)**three*d2eta + 
c    E        two*nu2*(nu2 + eta)**f52*(nu2 + lambda)*d2lambda + 
c    F        two*nu*(nu2 + eta)**three*(nu2 + lambda)*d2lambda - 
c    G        two*nu*(nu2 + eta)**two*(nu2 + lambda)**two*d2lambda - 
c    H        two*(nu2 + eta)**f52*(nu2 + lambda)**two*d2lambda)/
c    I        (four*(nu2 + eta)**f52*(nu2 + lambda)**f52)
c
c     d21Fxhse1 = -((A*(six*Ks4*zeta*d01Js*d10Js**two - 
c    1   two*Js*Ks4*(two*d1zeta*d01Js*d10Js + 
c    2   zeta*(two*d10Js*d11Js + d01Js*d20Js)) + 
c    3   Js**two*Ks4*(d2zeta*d01Js + two*d1zeta*d11Js + 
c    4   zeta*d21Js) + Js4*
c    5   (six*eta*d01Ks*d10Ks**two - 
c    6   two*Ks*(two*d1eta*d01Ks*d10Ks + 
c    7   two*eta*d10Ks*d11Ks + eta*d01Ks*d20Ks) + 
c    8   Ks**two*(d2eta*d01Ks + two*d1eta*d11Ks + 
c    9   eta*d21Ks))))/(Js4*Ks4))
c
c     d21Fxhse2 = (four*B*(d01lambda2*(six*d10lambda2**two - 
c    1   two*lambda2*d20lambda2) + 
c    2   lambda2*(-four*d10lambda2*d11lambda2 + lambda2*d21lambda2)))/
c    3   (nine*L4)
c
c     d21Fxhse3 = (two*C*(r24*(two + chi)*Fs*d01lambda2*
c    1   d10lambda2**two - six*lambda2*
c    2   (two*(two + chi)*d1Fs*d01lambda2*d10lambda2 + 
c    3   Fs*(d10lambda2*(d01chi*d10lambda2 + two*(two + chi)*
c    4   d11lambda2) + d01lambda2*(two*d10chi*d10lambda2 + 
c    5   (two + chi)*d20lambda2))) - L3*(d2Fs*d01chi + two*d1Fs*
c    6   d11chi + Fs*d21chi) + two*L2*((two + chi)*d2Fs*d01lambda2 + 
c    7   two*d1Fs*(d01lambda2*d10chi + d01chi*d10lambda2 + 
c    8   (two + chi)*d11lambda2) + Fs*(two*d10lambda2*d11chi + 
c    9   two*d10chi*d11lambda2 + d01lambda2*d20chi + d01chi*
c    A   d20lambda2 + two*d21lambda2 + chi*d21lambda2))))/
c    B  (nine*L5)
c
c     d21Fxhse4 = (-two*lambda2**two*(lambda2*d1EGs -
c    1   three*EGs*d10lambda2)*
c    2   (two*d01chi*d10chi + (three + two*chi)*d11chi) + 
c    3   six*(three + two*chi)*lambda2*d10chi*
c    4   (-four*EGs*d01lambda2*d10lambda2 + 
c    5   lambda2*(d1EGs*d01lambda2 + EGs*d11lambda2)) + 
c    6   three*EGs*lambda2**two*d01lambda2*
c    7   (two*d10chi**two + (three + two*chi)*d20chi) - 
c    8   (three + two*chi)*lambda2*d01chi*
c    9   (lambda2**two*d2EGs + r12*EGs*d10lambda2**two - 
c    A   three*lambda2*(two*d1EGs*d10lambda2 + EGs*d20lambda2)) - 
c    B   EGs*L3*(four*d10chi*d11chi + 
c    C   two*d01chi*d20chi + (three + two*chi)*d21chi) + 
c    D   (eight + nine*chi + three*chi**two)*
c    E   (r20*EGs*d01lambda2*d10lambda2**two - 
c    F   four*lambda2*(two*d1EGs*d01lambda2*d10lambda2 + 
c    G   EGs*(two*d10lambda2*d11lambda2 + 
c    H   d01lambda2*d20lambda2)) + 
c    I   lambda2**two*(d2EGs*d01lambda2 + two*d1EGs*d11lambda2 + 
c    J   EGs*d21lambda2)))/(three*L6)
c
c     d21Fxhse5 = (-two*D*(two*D**two*zeta*d01Ms*d10Ms**two - 
c    1   D*Ms*(two*D*d1zeta*d01Ms*d10Ms + 
c    2   zeta*(two*D*d10Ms*d11Ms + 
c    3   d01Ms*(six*d10Ms**two + D*d20Ms))) + 
c    4   Ms**four*(d2zeta*d01Ms + two*d1zeta*d11Ms + 
c    5   zeta*d21Ms) - two*Ms3*
c    6   (D*d2zeta*d01Ms + 
c    7   two*d1zeta*(d01Ms*d10Ms + D*d11Ms) + 
c    8   zeta*(two*d10Ms*d11Ms + d01Ms*d20Ms + 
c    9   D*d21Ms)) + Ms**two*
c    A   (D**two*d2zeta*d01Ms + 
c    B   two*D*d1zeta*(three*d01Ms*d10Ms + D*d11Ms) + 
c    C   zeta*(three*d01Ms*(two*d10Ms**two + D*d20Ms) + 
c    D   D*(six*d10Ms*d11Ms + D*d21Ms)))))/((D - Ms)**three*Ms3)
c
c     d21Fxhse6 = (two*(A - D)*(Ns**two*(A - D + Ns)**two*d2eta*d01Ns + 
c    1   two*(A - D)*Ns*(A - D + Ns)*d1eta*d01Ns*d10Ns - 
c    2   four*Ns*(A - D + Ns)**two*d1eta*d01Ns*d10Ns + 
c    3   two*Ns**two*(A - D + Ns)**two*d1eta*d11Ns + 
c    4   eta*(d01Ns*(two*((A - D)**two + three*(A - D)*Ns +
c    5   three*Ns**two)*d10Ns**two - 
c    6   Ns*((A - D)**two + three*(A - D)*Ns + two*Ns**two)*d20Ns) + 
c    7   Ns*(A - D + Ns)*(-two*(A - D + two*Ns)*d10Ns*d11Ns + 
c    8   Ns*(A - D + Ns)*d21Ns))))/(Ns3*(A - D + Ns)**three)
c
c     d21Fxhse = d21Fxhse1 + d21Fxhse2 + d21Fxhse3 + d21Fxhse4 +
c    1           d21Fxhse5 + d21Fxhse6
c
c     Calculate the mixed third derivative of Fxhse (dnu^2 ds)
c
c
c     Calculate mixed partial derivatives for the variables depending
c     on both s and nu.
c
c     d12chi = (three*nu*(-two*nu2 + three*lambda)*d1lambda)/
c    1         (two*(nu2 + lambda)**f72)
c     d12lambda2 = ((-one + chi)*d1lambda*(-two*d01chi**two +
c    1   (-one + chi)*d02chi) + 
c    2   lambda*(six*d01chi**two*d10chi - 
c    3   four*(-one + chi)*d01chi*d11chi + 
c    4   (-one + chi)*(-two*d02chi*d10chi + (-one + chi)*d12chi)))/
c    5   (-one + chi)**four
c     d12Js = (eta**three*(two*nu2 - zeta)*d1zeta + 
c    1   nu2*zeta*((three*nu4 + five*nu2*zeta + two*zeta**two)*d1eta - 
c    2   three*nu3*(nu + dsqrt(nu2 + eta))*d1zeta) - 
c    3   eta*((nu2 + zeta)*(zeta**two + three*nu3*(nu + 
c    4   dsqrt(nu2 + zeta)) + 
c    5   nu*zeta*(two*nu + three*dsqrt(nu2 + zeta)))*
c    6   d1eta - nu2*(three*nu4 - nu*(five*nu + 
c    7   six*dsqrt(nu2 + eta))*zeta + zeta**two)*d1zeta) + 
c    8   eta**two*(five*nu4*d1zeta + zeta**two*(d1eta + d1zeta) + 
c    9   nu*zeta*(nu*d1eta - three*(nu + dsqrt(nu2 + eta))*d1zeta)))/
c    A   (two*(nu2 + eta)**f52*(nu2 + zeta)**f52)
c     d12Ks = d12Js
c     d12Ms = (lambda**three*(two*nu2 - zeta)*d1zeta + 
c    1   nu2*zeta*((three*nu4 + five*nu2*zeta + 
c    2   two*zeta**two)*d1lambda - three*nu3*(nu + dsqrt(nu2 + lambda))*
c    3   d1zeta) - lambda*((nu2 + zeta)*(zeta**two + 
c    4   three*nu3*(nu + dsqrt(nu2 + zeta)) + 
c    5   nu*zeta*(two*nu + three*dsqrt(nu2 + zeta)))*d1lambda - 
c    6   nu2*(three*nu4 - nu*(five*nu + six*dsqrt(nu2 + lambda))*zeta + 
c    7   zeta**two)*d1zeta) + lambda**two*(five*nu4*d1zeta + 
c    8   zeta**two*(d1lambda + d1zeta) + 
c    9   nu*zeta*(nu*d1lambda - three*(nu + dsqrt(nu2 + lambda))*
c    A   d1zeta)))/(two*(nu2 + lambda)**f52*(nu2 + zeta)**f52)
c     d12Ns = (eta**three*(two*nu2 - lambda)*d1lambda + nu2*lambda*
c    1   ((three*nu4 + five*nu2*lambda + two*lambda**two)*d1eta - 
c    2   three*nu3*(nu + dsqrt(nu2 + eta))*d1lambda) - 
c    3   eta*((nu2 + lambda)*(lambda**two + 
c    4   three*nu3*(nu + dsqrt(nu2 + lambda)) + 
c    5   nu*lambda*(two*nu + three*dsqrt(nu2 + lambda)))*d1eta - 
c    6   nu2*(three*nu4 - nu*(five*nu + six*dsqrt(nu2 + eta))*lambda + 
c    7   lambda**two)*d1lambda) + eta**two*(five*nu4*d1lambda + 
c    8   lambda**two*(d1eta + d1lambda)+ nu*lambda*(nu*d1eta - 
c    9   three*(nu + dsqrt(nu2 + eta))*d1lambda)))/
c    A   (two*(nu2 + eta)**f52*(nu2 + lambda)**f52)
c
c     d12Fxhse1 = -((A*(six*Ks4*zeta*d01Js**two*d10Js - 
c    1   two*Js*Ks4*(d1zeta*d01Js**two + 
c    2   zeta*(d02Js*d10Js + two*d01Js*d11Js)) + 
c    3   Js**two*Ks4*(d1zeta*d02Js + zeta*d12Js) + 
c    4   Js4*(six*eta*d01Ks**two*d10Ks - 
c    5   two*Ks*(d1eta*d01Ks**two + eta*d02Ks*d10Ks + 
c    6   two*eta*d01Ks*d11Ks) + 
c    7   Ks**two*(d1eta*d02Ks + eta*d12Ks))))/(Js4*Ks4))
c
c     d12Fxhse2 = (four*B*(six*d01lambda2**two*d10lambda2 - 
c    1   four*lambda2*d01lambda2*d11lambda2 + 
c    2   lambda2*(-two*d02lambda2*d10lambda2 + lambda2*d12lambda2)))/
c    3   (nine*L4)
c
c     d12Fxhse3 = (-five*C*nu*(two*(nu2 + lambda)*d1Fs - 
c    1   seven*Fs*d1lambda))/(three*(nu2 + lambda)**f92)
c
c     d12Fxhse4 = (r20*(eight + nine*chi + three*chi**two)*EGs*
c    1   d01lambda2**two*d10lambda2 - 
c    2   four*lambda2*((eight + nine*chi + three*chi**two)*
c    3   d1EGs*d01lambda2**two + 
c    4   EGs*(three*(three + two*chi)*d01lambda2**two*d10chi + 
c    5   (eight + nine*chi + three*chi**two)*d02lambda2*d10lambda2 + 
c    6   two*d01lambda2*(three*(three + two*chi)*d01chi*d10lambda2 + 
c    7   (eight + nine*chi + three*chi**two)*d11lambda2))) - 
c    8   L3*(d1EGs*(two*d01chi**two + 
c    9   (three + two*chi)*d02chi) + 
c    A   EGs*(two*d02chi*d10chi + four*d01chi*d11chi + 
c    B   (three + two*chi)*d12chi)) + 
c    C   lambda2**two*(d1EGs*(six*(three + two*chi)*d01chi*d01lambda2 + 
c    D   (eight + nine*chi + three*chi**two)*d02lambda2) + 
c    E   EGs*(three*(three + two*chi)*d02lambda2*d10chi + 
c    F   six*d01chi**two*d10lambda2 + nine*d02chi*d10lambda2 + 
c    G   six*chi*d02chi*d10lambda2 + 
c    H   r18*d01lambda2*d11chi + 
c    I   r12*chi*d01lambda2*d11chi + 
c    J   six*d01chi*(two*d01lambda2*d10chi + 
c    K   (three + two*chi)*d11lambda2) + eight*d12lambda2 + 
c    L   nine*chi*d12lambda2 + three*chi**two*d12lambda2)))/
c    M   (three*L6)
c
c     d12Fxhse5 = (-two*D*(two*D**two*zeta*d01Ms**two*d10Ms - 
c    1   D*Ms*(D*d1zeta*d01Ms**two + zeta*(six*d01Ms**two*
c    2   d10Ms + D*d02Ms*d10Ms + two*D*d01Ms*d11Ms)) + 
c    3   Ms4*(d1zeta*d02Ms + zeta*d12Ms) - 
c    4   two*Ms3*(d1zeta*(d01Ms**two + D*d02Ms) + 
c    5   zeta*(d02Ms*d10Ms + two*d01Ms*d11Ms + D*d12Ms)) + 
c    6   Ms2*(D*d1zeta*(three*d01Ms**two + D*d02Ms) + 
c    7   zeta*(six*d01Ms**two*d10Ms + three*D*d02Ms*d10Ms + 
c    8   six*D*d01Ms*d11Ms + D**two*d12Ms))))/
c    9  ((D - Ms)**three*Ms3)
c
c     d12Fxhse6 = (two*(A - D)*(Ns*(A - D + Ns)*d1eta*
c    1   ((-A + D - two*Ns)*d01Ns**two + Ns*(A - D + Ns)*d02Ns) + 
c    2   (A - D)*Ns*eta*(-two*d01Ns**two + (A - D + Ns)*d02Ns)*
c    3   d10Ns + two*(A - D)*(A - D + Ns)*eta*d01Ns*
c    4   (-two*d01Ns*d10Ns + Ns*d11Ns) + 
c    5   (A - D + Ns)**two*eta*(six*d01Ns**two*d10Ns - 
c    6   four*Ns*d01Ns*d11Ns + 
c    7   Ns*(-two*d02Ns*d10Ns + Ns*d12Ns))))/
c    8   (Ns3*(A - D + Ns)**three)
c
c     d12Fxhse = d12Fxhse1 + d12Fxhse2 + d12Fxhse3 + d12Fxhse4 +
c    1           d12Fxhse5 + d12Fxhse6
#endif
      
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_hse08.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_hse08.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_hse08.F"
#endif
#undef NWAD_PRINT
c $Id$
C> @}
