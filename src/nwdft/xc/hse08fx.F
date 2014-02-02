#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine HSE08Fx(ipol,rho,s,Fxhse,d10Fxhse,d01Fxhse)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine HSE08Fx_d2(ipol,rho,s,Fxhse,d10Fxhse,d01Fxhse,
     &                      d20Fxhse,d02Fxhse,d11Fxhse)
#else
      Subroutine HSE08Fx_d3(ipol,rho,s,Fxhse,d10Fxhse,d01Fxhse,
     &                   d20Fxhse,d02Fxhse,d11Fxhse,d30Fxhse,
     &                   d21Fxhse,d12Fxhse,d03Fxhse)
#endif
c
c$Id$
c
      implicit none
c
ccase...start
#include "case.fh"
ccase...end
c
c HSE evaluates the Heyd et al. Screened Coulomb
c Exchange Functional
c
c Calculates the enhancement factor
c
      integer ipol
      double precision rho,s,Fxhse,d10Fxhse,d01Fxhse
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d20Fxhse,d02Fxhse,d11Fxhse
#endif
c
#ifdef THIRD_DERIV
      double precision d30Fxhse,d03Fxhse,d21Fxhse,d12Fxhse
#endif

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
      double precision  zeta, lambda, eta, kf, nu, chi, lambda2
      double precision  d1zeta,d1lambda,d1eta,d1nu,d10chi,d10lambda2
      double precision  d01chi,d01lambda2
      double precision  EGs,d1EGs
      double precision  nu2,L2,L3,nu3,nu4,nu5,nu6,nu7,nu8
      double precision  Js,Ks,Ms,Ns
      double precision  Js2,Ks2,Ms2,Ns2
      double precision  d10Js,d10Ks,d10Ms,d10Ns
      double precision  d01Js,d01Ks,d01Ms,d01Ns

      double precision  tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
      double precision  tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15
      double precision  Fxhse1,Fxhse2,Fxhse3,Fxhse4,Fxhse5,Fxhse6
      double precision  d1Fxhse1,d1Fxhse2,d1Fxhse3,d1Fxhse4,d1Fxhse5
      double precision  d1Fxhse6,d1Fxhse7

      double precision  r42,r27,r12,r15,r14,r18,r20,r30,r56,r72
      double precision  r16,r32,r24,r48,r11,r64,r35,r60,r36,r63
      double precision  r189,r256,r120,r210,r336,r504,r21,r105,r126
      double precision  pi,pi2,srpi,s02
      double precision  f12,f13,f32,f52,f72,f92
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision  d2H,d2hnum,d2hden
      double precision  d2Fs
      double precision  d2zeta,d2lambda,d2eta,d2nu,d20chi,d20lambda2
      double precision  d02chi,d02lambda2,d11chi,d11lambda2
      double precision  d2EGs
      double precision  L5
      double precision  Js3,Ks3,Ms4,Ns4
      double precision  d20Js,d20Ks,d20Ms,d20Ns
      double precision  d02Js,d02Ks,d02Ms,d02Ns
      double precision  d11Js,d11Ks,d11Ms,d11Ns
      double precision  d2Fxhse1,d2Fxhse2,d2Fxhse3,d2Fxhse4,d2Fxhse5
      double precision  d2Fxhse6,d2Fxhse7
      double precision  d11Fxhse1,d11Fxhse2,d11Fxhse3,d11Fxhse4
      double precision  d11Fxhse5,d11Fxhse6,d11Fxhse7
#endif
c
#ifdef THIRD_DERIV
      double precision  d3H,d3hnum,d3hden
      double precision  d3Fs
      double precision  d3zeta,d3lambda,d3eta,d3nu,d30chi,d30lambda2
      double precision  d03chi,d03lambda2,d21chi,d12chi
      double precision  d21lambda2,d12lambda2
      double precision  L4,L6
      double precision  Js4,Ks4,Ms3,Ns3
      double precision  d3EGs
      double precision  d30Js,d30Ks,d30Ms,d30Ns
      double precision  d03Js,d03Ks,d03Ms,d03Ns
      double precision  d21Js,d21Ks,d21Ms,d21Ns
      double precision  d12Js,d12Ks,d12Ms,d12Ns
      double precision  d3Fxhse1,d3Fxhse2,d3Fxhse3,d3Fxhse4,d3Fxhse5
      double precision  d3Fxhse6,d3Fxhse7
      double precision  d21Fxhse1,d21Fxhse2,d21Fxhse3,d21Fxhse4
      double precision  d21Fxhse5,d21Fxhse6,d21Fxhse7
      double precision  d12Fxhse1,d12Fxhse2,d12Fxhse3,d12Fxhse4
      double precision  d12Fxhse5,d12Fxhse6,d12Fxhse7
#endif
c
c     Constants for HJS hole
c
      Data A,B,C,D,E
     1     / 7.57211D-1,-1.06364D-1,-1.18649D-1,
     2       6.09650D-1,-4.77963D-2 /
c
c     Constants for fit of H(s) (PBE hole)
c     Taken from JCTC_5_754 (2009)
c
      Data ha2,ha3,ha4,ha5,ha6,ha7
     1     / 1.59941D-2,8.52995D-2,-1.60368D-1,1.52645D-1,
     2      -9.71263D-2,4.22061D-2 /
c
      Data hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9
     1     / 5.33319D0,-12.4780D0,11.0988D0,-5.11013D0,
     2      1.71468D0,-6.10380D-1,3.07555D-1,-7.70547D-2,
     3      3.34840D-2 /

c
c     Whole numbers used during evaluation
c
      Data zero,one,two,three,four,five,six,seven,eight,nine,ten
     1     / 0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0 /
       
      Data r11,r12,r14,r15,r16,r18,r20,r24,r27,r30,r32
     1     / 11.D0,12.D0,14.D0,15.D0,16.D0,18.D0,20.D0,24.D0,27.d0,
     2       30.D0,32.D0 /

      Data r35,r42,r48,r56,r64,r72
     1     / 35.D0,42.D0,48.D0,56.D0,64.D0,72.D0 /

      Data r21,r36,r60,r63,r120,r189,r210,r256,r336,r504
     1     / 21.D0,36.D0,60.D0,63.D0,120.D0,189.D0,210.D0,256.D0,
     2       336.D0,504.D0 /

      Data r105,r126
     1     / 105.D0,126.D0 /
c
c     Fractions used during evaluation
c
      Data f12
     1     / 0.5D0 /
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
     1       hb6*s6 + hb7*s7 + hb8*s8 + hb9*s9
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
      chi = nu/dsqrt(lambda+nu**two)
      lambda2 = (one+chi)*(lambda+nu**two)

c
c     Calculate F(H(s)) for the model exhange hole
c
      Fs = one-s2/(r27*C*(one+s02))-zeta/(two*C)

c
c     Calculate EGs 
c
      EGs = -(two/five)*C*Fs*lambda - (four/r15)*B*lambda*lambda -
     1      (six/five)*A*lambda*lambda*lambda - 
     2      (four/five)*srpi*lambda**(seven/two) -
     3      (r12/five)*(lambda**(seven/two))*(dsqrt(zeta)-dsqrt(eta))
 
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

!*********************************************************
! Calculate the First Derivative of the Enhancement Factor
!*********************************************************

c
c     Calculate the first derivative of H with respect to the
c     reduced density gradient, s.
c
      d1hnum = two*ha2*s + three*ha3*s2 + four*ha4*s3 +
     1          five*ha5*s4 + six*ha6*s5 + seven*ha7*s6

      d1hden  = hb1 + two*hb2*s +three*hb3*s2 + four*hb4*s3 +
     1          five*hb5*s4 + six*hb6*s5 + seven*hb7*s6 +
     2          eight*hb8*s7 + nine*hb9*s8 
      d1H =   (hden*d1hnum -hnum*d1hden)/hden**two

c
c     calculate first derivative of variables needed with respect to s
c 
      d1zeta = two*s*H + s2*d1H
      d1eta  = d1zeta
      d1lambda = d1zeta
      d10chi = -f12*nu*d1zeta/(lambda + nu2)**f32
      d10lambda2 = d10chi*(lambda + nu**two) + (one+chi)*d1lambda

c   
c     calculate the first derivative of Fs with respect to s
c   
      d1Fs = -two*s/(r27*C*(one+s02)**two) - d1zeta/(two*C)

c
c     Calculate the first derivate of EGs with respect to s
c
      d1EGs = -(two/five)*C*(d1Fs*lambda + Fs*d1lambda) -
     1        (eight/r15)*B*lambda*d1lambda -
     2        (r18/five)*A*lambda*lambda*d1lambda -
     3        (r14/five)*srpi*d1lambda*lambda**f52 -
     4        (r42/five)*(lambda**f52)*
     5        d1lambda*(dsqrt(zeta)-dsqrt(eta))-
     6        (six/five)*(lambda**(seven/two))*
     7        (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))

c
c     Calculate the first derivate of denominators needed with respect
c     to s
c
      tmp1 = (dsqrt(zeta+nu2)+nu)/(dsqrt(eta+nu2)) 
      tmp2 = (dsqrt(eta+nu2)+nu)/(dsqrt(zeta+nu2))

      d10Js = f12*d1zeta*(two+tmp1+tmp2)
      d10Ks = d10Js

      tmp3 = (dsqrt(zeta+nu2)+nu)/(dsqrt(lambda+nu2))
      tmp4 = (dsqrt(lambda+nu2)+nu)/(dsqrt(zeta+nu2)) 
      d10Ms = f12*d1zeta*(two +tmp3+tmp4)

      tmp5 = (dsqrt(lambda+nu2)+nu)/(dsqrt(eta+nu2))
      tmp6 = (dsqrt(eta+nu2)+nu)/(dsqrt(lambda+nu2))
      d10Ns = f12*d1zeta*(two + tmp5+tmp6)
c
c
c     Calculate the derivative of the 08-Fxhse with respect to s
c
      L2 = lambda2*lambda2
      L3 = lambda2*lambda2*lambda2
      Js2 = Js*Js
      Ks2 = Ks*Ks
      Ms2 = Ms*Ms
      Ns2 = Ns*Ns
c
      d1Fxhse1  = A*( (Js*d1zeta - zeta*d10Js)/(Js2) +
     1                (Ks*d1zeta - eta*d10Ks)/(Ks2) ) 
c
      d1Fxhse2  = (four/nine)*B*d10lambda2/L2 
c
      tmp9 = d10lambda2/lambda2
      tmp7 = d1Fs - two*Fs*tmp9
      tmp8 = one + f12*chi
      tmp10 =  f12*Fs*d10chi

      d1Fxhse3 = -(four*C/(nine*L2))*(tmp7*tmp8+tmp10)
c
       tmp7 = one + (nine/eight)*chi+(three/eight)*chi*chi
       tmp8 = (nine/eight)*d10chi + (six/eight)*chi*d10chi

      d1Fxhse4 = -(eight/(nine*L3))*((d1EGs-three*EGs*tmp9)*tmp7
     1           + EGs*tmp8)
c
      d1Fxhse5  = two*d1zeta*dlog(one-D/Ms) +
     1           two*zeta*D*d10Ms/(Ms2*(one-D/Ms)) 
c
      d1Fxhse6  = -two*d1eta*dlog(one- (D-A)/Ns) -
     1           two*eta*(D-A)*d10Ns/(Ns2*(one-(D-A)/Ns)) 
c
      d10Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6
c
c
c     Calculate the derivative of 08-Fxhse with respect to nu
c
      nu3 = nu2*nu
c
c     Define all derivatives with respect to nu necessary to 
c     evaluate enhancement factor derivatives.
c
      d01chi = lambda/(nu2 + lambda)**f32
      d01lambda2 = (lambda*d01chi)/(one - chi)**two
      d01Js = (eta*(nu + dsqrt(nu2 + zeta)) + 
     1        (nu + dsqrt(nu2 + eta))*
     2        (zeta + two*nu*(nu + dsqrt(nu2 + zeta))))/
     3        (dsqrt(nu2 + eta)*dsqrt(nu2 + zeta))
      d01Ks = d01Js 
      d01Ms = (lambda*(nu + dsqrt(nu2 + zeta)) + 
     1        (nu + dsqrt(nu2 + lambda))*
     2        (zeta + two*nu*(nu + dsqrt(nu2 + zeta))))/
     3        (dsqrt(nu2 + lambda)*dsqrt(nu2 + zeta))
      d01Ns = (eta*(nu + dsqrt(nu2 + lambda)) + 
     1        (nu + dsqrt(nu2 + eta))*
     2        (lambda + two*nu*(nu + dsqrt(nu2 + lambda))))/
     3        (dsqrt(nu2 + eta)*dsqrt(nu2 + lambda))
c
      d1Fxhse1 = -((A*(nu + dsqrt(eta + nu2))*zeta)/
     1            (dsqrt(eta + nu2)*dsqrt(nu2 + zeta)*
     2            (nu + dsqrt(nu2 + zeta))*
     3            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))))
c
      d1Fxhse2 = -((A*eta*(nu/dsqrt(eta + nu2) + nu/
     1            dsqrt(nu2 + zeta)))/
     2            ((nu + dsqrt(eta + nu2))*
     3            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -
     4            (A*eta*(one + nu/dsqrt(eta + nu2)))/
     5            ((nu + dsqrt(eta + nu2))**two*
     6            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta)))
c
      d1Fxhse3 = (four*B)/(nine*(lambda + nu2)**(f32))
c
      d1Fxhse4 = (two*C*Fs)/(three*(lambda + nu2)**(f52))
c
      d1Fxhse5 = (five*EGs*(lambda**two + four*nu3*
     1            (nu + dsqrt(lambda + nu2)) +
     2            lambda*nu*(five*nu + three*dsqrt(lambda + nu2))))/
     3   (three*(lambda + nu2)**four*(nu + dsqrt(lambda + nu2))**three)
c
      d1Fxhse6 = (two*D*zeta*(nu + dsqrt(nu2 + zeta)))/
     1            (dsqrt(lambda + nu2)*dsqrt(nu2 + zeta)*
     2            (-D + lambda + (nu + dsqrt(lambda + nu2))*
     3            (nu + dsqrt(nu2 + zeta))))
c
      d1Fxhse7 = (two*(A - D)*eta*(nu + dsqrt(eta + nu2)))/
     1            (dsqrt(eta + nu2)*dsqrt(lambda + nu2)*
     2            (A - D + lambda + nu2 + nu*dsqrt(eta + nu2) +
     3            nu*dsqrt(lambda + nu2) +
     4            dsqrt(eta + nu2)*dsqrt(lambda + nu2)))
c
      d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+
     1           d1Fxhse7
c   

!**********************************************************
! Calculate the Second Derivative of the Enhancement Factor
!**********************************************************
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Calculate the second derivative of H
c
      d2hnum =  two*ha2+six*ha3*s+r12*ha4*s2+r20*ha5*s3+
     1           r30*ha6*s4 + r42*ha7*s5

      d2hden  =  two*hb2+six*hb3*s+r12*hb4*s2+r20*hb5*s3+
     1           r30*hb6*s4+r42*hb7*s5+r56*hb8*s6+r72*hb9*s7

      d2H     =  (hden*(hden*d2hnum-hnum*d2hden) -
     1           two*(hden*d1hnum-hnum*d1hden)*d1hden)/hden**three   

c
c     calculate second derivative of variables needed
c 
      d2zeta    = two*H +four*s*d1H + s2*d2H
      d2eta     = d2zeta
      d2lambda  = d2zeta

      d20chi     = -f12*(nu/(lambda+nu2)**(f32))*
     1           (-f32*d1zeta*d1zeta/(lambda+nu2) +d2zeta)

      d20lambda2 =(one-chi)*(d2lambda-d2lambda*chi+lambda*d20chi)+
     1           two*d10chi*(d1lambda-d1lambda*chi+lambda*d10chi)
      d20lambda2 = d20lambda2/(one-chi)**three
c   
c     calculate the second derivative of Fs
c   
      d2Fs = -(two/(r27*C))*(one-three*s02)/
     1        ((one+s02)**three)-d2zeta/(two*C)  

c
c     Calculate the second derivative of EGs
c
      tmp7 = -(two/five)*C*(d2Fs*lambda+two*d1Fs*d1lambda+Fs*d2lambda)
      tmp8 = -(eight/r15)*B*(d1lambda*d1lambda+lambda*d2lambda)
      tmp9 =-(r18/five)*A*lambda*(two*d1lambda*d1lambda+lambda*d2lambda)
      tmp10 = -(r14/five)*srpi*(f52*(lambda**f32)*d1lambda*d1lambda
     1          +(lambda**f52)*d2lambda)
      tmp11 = -(r42/five)*(f52*(lambda**f32)*d1lambda*d1lambda
     1          +(lambda**f52)*d2lambda)*(dsqrt(zeta)-dsqrt(eta))
      tmp12 = -(r42/five)*(lambda**f52)*d1lambda*
     1          (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))
      tmp13 = -(six/five)*(lambda**(seven/two))*
     1          (-d1zeta*d1zeta/(two*zeta**f32)+d2zeta/dsqrt(zeta)
     2          +d1eta*d1eta/(two*eta**f32)-d2eta/dsqrt(eta))

      d2EGs = tmp7+tmp8+tmp9+tmp10+tmp11+tmp12+tmp13
c
c     Calculate the second derivative of denominators needed
c
      tmp7 = two/(dsqrt(zeta+nu2)*dsqrt(eta+nu2))
      tmp8 = (dsqrt(eta+nu2)+nu)/(zeta+nu2)**f32
      tmp9 = (dsqrt(zeta+nu2)+nu)/(eta+nu2)**f32
      
      d20Js = f12*(d2zeta*(two+tmp1+tmp2) +
     1       f12*d1zeta*d1zeta*(tmp7-tmp8-tmp9)) 

      d20Ks = d20Js

      tmp10 = two/(dsqrt(zeta+nu2)*dsqrt(lambda+nu2))
      tmp11 = (dsqrt(lambda+nu2)+nu)/(zeta+nu2)**f32
      tmp12 = (dsqrt(zeta+nu2)+nu)/(lambda+nu2)**f32
      d20Ms = f12*(d2zeta*(two+tmp3+tmp4) +
     1       f12*d1zeta*d1zeta*(tmp10-tmp11-tmp12)) 

      tmp13 = two/(dsqrt(eta+nu2)*dsqrt(lambda+nu2))
      tmp14 = (dsqrt(lambda+nu2)+nu)/(eta+nu2)**f32
      tmp15 = (dsqrt(eta+nu2)+nu)/(lambda+nu2)**f32
      d20Ns = f12*(d2zeta*(two+tmp5+tmp6) +
     1       f12*d1zeta*d1zeta*(tmp13-tmp14-tmp15)) 


c
c     Calculate the second derivative of the Fxhse with respect to the
c     reduced density gradient, s.
c
      Ms4 = Ms2*Ms2
      Ns4 = Ns2*Ns2
c
      tmp1      = A*(Js2*d2zeta -two*Js*d10Js*d1zeta 
     1            + two*d10Js*d10Js*zeta-Js*d20Js*zeta)/
     2            (Js2*Js) 
      tmp2       = A*(Ks2*d2zeta -two*Ks*d10Ks*d1zeta 
     1            + two*d10Ks*d10Ks*eta-Ks*d20Ks*eta)/
     2            (Ks2*Ks) 
      d2Fxhse1  = tmp1 +tmp2 
c
      d2Fxhse2  = (four/nine)*B*(L2*d20lambda2-
     1            two*lambda2*d10lambda2*d10lambda2)/(L2*L2)
c
      tmp4 = d10lambda2/lambda2
      tmp5 = d20lambda2/lambda2
      d2Fxhse3  = -((four*C)/(nine*L2))*(d2Fs +
     1           six*Fs*tmp4**two - two*Fs*tmp5 -
     2           four*d1Fs*tmp4 +
     3          f12*((d2Fs*chi+two*d1Fs*d10chi+Fs*d20chi)-
     4           four*tmp4*(d1Fs*chi+Fs*d10chi) +
     5           six*chi*Fs*tmp4*tmp4)-chi*Fs*tmp5)
c
      tmp1 = -EGs*(eight/nine)/L3
      tmp2 = (-d1EGs + three*EGs*d10lambda2/lambda2)*(eight/nine)/L3
      tmp3 = (-d2EGs - r12*EGs*d10lambda2**two/L2 + 
     1   three*(two*d1EGs*d10lambda2 + EGs*d20lambda2)/lambda2)
     2   *(eight/nine)/L3
      tmp4 = (one + (nine/eight)*chi + (three/eight)*chi**two)
      tmp5 = (three*(three + two*chi)*d10chi)/eight
      tmp6 = (three*(two*d10chi**two + (three + two*chi)*d20chi))/eight

      d2Fxhse4 = (tmp1*tmp6+two*tmp2*tmp5+tmp3*tmp4)
c
      tmp1      = (one-D/Ms) 

      d2Fxhse5  = two*d2zeta*dlog(tmp1)+four*d1zeta*D*d10Ms/(Ms2*tmp1)
     1          - two*zeta*(D*d10Ms/(Ms2*tmp1))**two +
     2          two*zeta*D*(Ms2*d20Ms-two*Ms*d10Ms*d10Ms)/
     3          (tmp1*Ms**four)
c
      tmp1      = (one-(D-A)/Ns) 
      d2Fxhse6  = -two*d2eta*dlog(tmp1)-
     1             four*d1eta*(D-A)*d10Ns/(Ns2*tmp1)+
     2             two*eta*((D-A)*d10Ns/(Ns2*tmp1))**two -
     3       two*eta*(D-A)*(Ns2*d20Ns-two*Ns*d10Ns*d10Ns)/
     4             (tmp1*Ns**four)
c
      d20Fxhse = d2Fxhse1+d2Fxhse2+d2Fxhse3+d2Fxhse4+d2Fxhse5+d2Fxhse6
c
c
c     Calculate the second derivative of Fxhse with respect to the
c     parameter nu (cam_omega/kf).
c
      nu5 = nu3*nu2
      nu6 = nu5*nu
c
c     Calculate second derivatives with respect to nu necessary to
c     simplify the appearance of the code.
c
      d02chi = (-three*nu*lambda)/(nu2 + lambda)**f52
      d02lambda2 = -((lambda*(two*d01chi**two - 
     1             (-one + chi)*d02chi))/
     2             (-one + chi)**three)
      d02Js = two*nu*(one/dsqrt(nu2 + eta) + one/dsqrt(nu2 + zeta))*
     1        (one + nu/dsqrt(nu2 + zeta)) + 
     2        (zeta*(dsqrt(nu2 + eta) + 
     3        dsqrt(nu2 + zeta)))/(nu2 + zeta)**f32 + 
     4        (nu + dsqrt(nu2 + zeta))*
     5        (one/dsqrt(nu2 + eta) + one/dsqrt(nu2 + zeta) + 
     6        nu2*(-(nu2 + eta)**(-f32) - 
     7        (nu2 + zeta)**(-f32)))
      d02Ks = d02Js
      d02Ms = two*nu*(one + nu/dsqrt(nu2 + lambda))*
     1        (one/dsqrt(nu2 + lambda) + one/dsqrt(nu2 + zeta))
     2        + (lambda*(dsqrt(nu2 + lambda) +
     3        dsqrt(nu2 + zeta)))/(nu2 + lambda)**f32 +
     4        (nu + dsqrt(nu2 + lambda))*
     5        (one/dsqrt(nu2 + lambda) + one/dsqrt(nu2 + zeta) +
     6        nu2*(-(nu2 + lambda)**(-f32) -
     7        (nu2 + zeta)**(-f32)))
      d02Ns = two*nu*(one/dsqrt(nu2 + eta) + 
     1        one/dsqrt(nu2 + lambda))*
     2        (one + nu/dsqrt(nu2 + lambda)) + 
     3        (lambda*(dsqrt(nu2 + eta) + 
     4        dsqrt(nu2 + lambda)))/(nu2 + lambda)**f32
     5        + (nu + dsqrt(nu2 + lambda))*
     6        (one/dsqrt(nu2 + eta) + one/dsqrt(nu2 + lambda) + 
     7        nu2*(-(nu2 + eta)**(-f32) - 
     8        (nu2 + lambda)**(-f32))) 
c
      L5 = L2*L3
      Js3 = Js2*Js
      Ks3 = Ks2*Ks
c
      d2Fxhse1 = A*((zeta*(two*d01Js**two - Js*d02Js))/Js3 + 
     1    (eta*(two*d01Ks**two - Ks*d02Ks))/Ks3)
c
      d2Fxhse2 = (four*B*(-two*d01lambda2**two + lambda2*d02lambda2))/
     1    (nine*L3)
c
      d2Fxhse3 = -(ten*C*Fs*nu)/(three*(lambda+nu2)**f72)
c
      d2Fxhse4 = (EGs*(-four*
     1    (eight + nine*chi + three*chi**two)*d01lambda2**two - 
     2    lambda2**two*(two*d01chi**two + (three + two*chi)*d02chi) + 
     3    lambda2*(six*(three + two*chi)*d01chi*d01lambda2 + 
     4    (eight + nine*chi + three*chi**two)*d02lambda2)))/
     5    (three*L5)
c
      d2Fxhse5 = (two*D*zeta*((D - two*Ms)*d01Ms**two + 
     1    Ms*(-D + Ms)*d02Ms))/((D - Ms)**two*Ms**two)
c
      d2Fxhse6 = (two*(A - D)*eta*((-A + D - two*Ns)*d01Ns**two +
     1  Ns*(A - D + Ns)*d02Ns))/
     2  (Ns**two*(A - D + Ns)**two)
c
      d02Fxhse = d2Fxhse1 + d2Fxhse2 + d2Fxhse3 + d2Fxhse4 + 
     1           d2Fxhse5 + d2Fxhse6
c
c
c     Calculate the mixed partial derivative of the enhancement factor 
c     with respect to both the reduced density gradient, s, and the 
c     parameter nu (cam_omega/kF).  Note that the enhancement factor
c     satisfies continuity requirements and therefore it does not matter
c     what order the derivatives are taken in.
c     
      nu4 = nu3*nu
c
c     Calculate mixed partial derivatives for the variables depending
c     on both s and nu.
c
      d11chi = (three*nu2*d1lambda)/
     1         (two*(nu2 + lambda)**f52) - 
     2         d1lambda/(two*(nu2 + lambda)**f32) 
      d11lambda2 = (d1lambda*d01chi)/
     1             (one - chi)**two + 
     2             (two*lambda*d01chi*
     3             d10chi)/(one - chi)**three + 
     4             (lambda*d11chi)/
     5             (one - chi)**two
      d11Js = (-(nu*eta**two*d1zeta) + 
     1        nu*zeta*((-nu2 - zeta)*
     2        d1eta + 
     3        nu*(nu + dsqrt(nu2 + eta))*
     4        d1zeta) + 
     5        eta*((nu2 + zeta)*
     6        (nu + dsqrt(nu2 + zeta))*d1eta
     7        + (-nu3 + (nu + dsqrt(nu2 + eta))*zeta)*
     8        d1zeta))/
     9        (two*(nu2 + eta)**f32*(nu2 + zeta)**f32)
      d11Ks = d11Js
      d11Ms = (-(nu*lambda**two*d1zeta) + 
     1        nu*zeta*((-nu2 - zeta)*
     2        d1lambda + 
     3        nu*(nu + dsqrt(nu2 + lambda))*
     4        d1zeta) + 
     5        lambda*((nu2 + zeta)*
     6        (nu + dsqrt(nu2 + zeta))*
     7        d1lambda + 
     8        (-nu3 + (nu + dsqrt(nu2 + lambda))*zeta)*
     9        d1zeta))/
     A        (two*(nu2 + lambda)**f32*(nu2 + zeta)**f32)
      d11Ns = (-(nu*eta**two*d1lambda) + 
     1        nu*lambda*((-nu2 - lambda)*
     2        d1eta + 
     3        nu*(nu + dsqrt(nu2 + eta))*
     4        d1lambda) + 
     5        eta*((nu2 + lambda)*
     6        (nu + dsqrt(nu2 + lambda))*
     7        d1eta + 
     8        (-nu3 + (nu + dsqrt(nu2 + eta))*lambda)*
     9        d1lambda))/
     A        (two*(nu2 + eta)**f32*(nu2 + lambda)**f32)
c
      d11Fxhse1 = -((A*(-two*Ks3*zeta*d01Js*d10Js + 
     1   Js*Ks3*(d1zeta*d01Js + zeta*d11Js) + 
     2   Js3*(-two*eta*d01Ks*d10Ks + 
     3   Ks*(d1eta*d01Ks + eta*d11Ks))))/(Js3*Ks3))
c
      d11Fxhse2 = (-two*B*(two*nu2*(nu + dsqrt(nu2 + lambda)) +
     1   lambda*(two*nu + dsqrt(nu2 + lambda)))*d1lambda)/
     2   (three*(nu2 + lambda)**three*
     3   (nu + dsqrt(nu2 + lambda))**two)
c
      d11Fxhse3 = (C*(two*(nu2 + lambda)*d1Fs - five*Fs*d1lambda))/
     1   (three*(nu2 + lambda)**(f72))
c
      d11Fxhse4 = (five*(eight*nu4*(nu + dsqrt(nu2 + lambda)) +
     1   lambda**two*(four*nu + dsqrt(nu2 + lambda)) +
     2   four*nu2*lambda*(three*nu + two*dsqrt(nu2 + lambda)))*
     3   (two*(nu2 + lambda)*d1EGs - seven*EGs*d1lambda))/
     4   (six*(nu2 + lambda)**five*(nu + dsqrt(nu2 + lambda))**four)
c
      d11Fxhse5 = (two*D*(D*zeta*d01Ms*d10Ms + 
     1   Ms**two*(d1zeta*d01Ms + zeta*d11Ms) - 
     2   Ms*(D*d1zeta*d01Ms + 
     3   zeta*(two*d01Ms*d10Ms + D*d11Ms))))/((D - Ms)**two*Ms**two)
c
      d11Fxhse6 = (two*(A - D)*((-A + D)*eta*d01Ns*d10Ns + 
     1   Ns**two*(d1eta*d01Ns + eta*d11Ns) + 
     2   Ns*((A - D)*d1eta*d01Ns + 
     3   eta*(-two*d01Ns*d10Ns + (A - D)*d11Ns))))/
     4   (Ns**two*(A - D + Ns)**two)
c
      d11Fxhse = d11Fxhse1 + d11Fxhse2 + d11Fxhse3 + d11Fxhse4 + 
     1           d11Fxhse5 + d11Fxhse6
c
#endif
!*********************************************************
! Calculate the Third Derivative of the Enhancement Factor
!*********************************************************
#ifdef THIRD_DERIV
c
c    Calculate the third order derivative of H with respect to s.
c
      d3hnum = six*ha3+r24*ha4*s+r60*ha5*s2+
     1         r120*ha6*s3 + r210*ha7*s4

      d3hden  =  six*hb3+r24*hb4*s+r60*hb5*s2+
     1           r120*hb6*s3+r210*hb7*s4+r336*hb8*s5+r504*hb9*s6

      d3H = (-six*hnum*d1hden**three +
     1  six*hden*d1hden*(d1hnum*d1hden + hnum*d2hden) +
     2  hden**three*d3hnum - hden**two*(three*d1hden*d2hnum +
     3  three*d1hnum*d2hden + hnum*d3hden))/hden**four
c
c    Calculate the third order derivatives of the "helper variables"
c
      d3zeta = six*d1H + six*s*d2H + s2*d3H
      d3eta = d3zeta
      d3lambda = d3zeta
      d30chi = -(nu*(r15*d1lambda*d1lambda*d1lambda
     1    - r18*(nu2 + lambda)*d1lambda*
     2   d2lambda + four*(nu2 + lambda)**two*d3lambda))/
     3   (eight*(nu2 + lambda)**f72) 
      d30lambda2 = d3lambda/(one - chi) + 
     1             (three*d2lambda*d10chi)/
     2             (-one + chi)**two - 
     3             (three*d1lambda*
     4             (two*d10chi**two - 
     5             (-one + chi)*d20chi))/
     6             (-one + chi)**three + 
     7             (lambda*(six*d10chi**three - 
     8             six*(-one + chi)*d10chi*
     9             d20chi + 
     A             (-one + chi)**two*d30chi))/
     B             (-one + chi)**four
c
c    Calculate the third order derivative of Fs
c
      d3Fs = -(r256*s*(s2 - four) +
     1       nine*(s2 + four)**four*d3zeta)/
     2       (r18*C*(s2 + four)**four)

c
c    Calculate the third order derivative of EGs
c 
      tmp1 = (-two*C*(three*d1lambda*d2Fs + 
     1  three*d1Fs*d2lambda + lambda*d3Fs + 
     2  Fs*d3lambda))/five
      tmp2 = (-four*B*(six*d1lambda*d2lambda + 
     1  two*lambda*d3lambda))/r15
      tmp3 = (-six*A*(six*d1lambda**three + 
     1  r18*lambda*d1lambda*d2lambda + 
     2  three*lambda**two*d3lambda))/five
      tmp4 = (-four*srpi*((r105*dsqrt(lambda)*
     1  d1lambda**three)/eight + (r105*lambda**f32*d1lambda*
     2  d2lambda)/four + (seven*lambda**f52*d3lambda)/two))/five
      tmp5 = (-r36*(-d1eta/(two*dsqrt(eta)) + 
     1  d1zeta/(two*dsqrt(zeta)))*
     2  ((r35*lambda**f32*d1lambda**two)/four +
     3  (seven*lambda**f52*d2lambda)/two))/five 
     4  - (r126*lambda**f52*d1lambda*(d1eta**two/(four*eta**f32) - 
     5  d1zeta**two/(four*zeta**f32) - d2eta/(two*dsqrt(eta)) + 
     6  d2zeta/(two*dsqrt(zeta))))/five - 
     7  (r12*(-dsqrt(eta) + dsqrt(zeta))*
     8  ((r105*dsqrt(lambda)*d1lambda**three)/
     9  eight + (r105*lambda**f32*d1lambda*d2lambda)/four + 
     A  (seven*lambda**f52*d3lambda)/two))/five
     B  - (r12*lambda**f72*((-three*d1eta**three)/(eight*eta**f52) + 
     C  (three*d1zeta**three)/(eight*zeta**f52) + 
     D  (three*d1eta*d2eta)/(four*eta**f32) - 
     E  (three*d1zeta*d2zeta)/(four*zeta**f32) - 
     F  d3eta/(two*dsqrt(eta)) + d3zeta/(two*dsqrt(zeta))))/five

      d3EGs = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 

c
c    Calculate derivatives of denominators needed for third derivatives
c
      d30Js = (three*(nu + dsqrt(nu2 + zeta))*
     1        d1eta**three - 
     2        (three*(nu2 + eta)*d1eta**two*
     3        d1zeta)/dsqrt(nu2 + zeta) - 
     4        (three*(nu2 + eta)*d1eta*
     5        ((nu2 + eta)*d1zeta**two + 
     6        two*(nu2 + zeta)*
     7        ((zeta + nu*(nu + dsqrt(nu2 + zeta)))*
     8        d2eta - 
     9        (nu2 + eta)*d2zeta)))/
     A        (nu2 + zeta)**f32 + 
     B        ((nu2 + eta)**two*
     C        (three*(eta + nu*(nu + dsqrt(nu2 + eta)))*
     D        d1zeta**three + 
     E        six*(nu2 + zeta)*d1zeta*
     F        ((nu2 + zeta)*d2eta - 
     G        (eta + nu*(nu + dsqrt(nu2 + eta)))*
     H        d2zeta) + 
     I        four*(nu2 + zeta)**two*
     J        ((zeta + nu*(nu + dsqrt(nu2 + zeta)))*
     K        d3eta + 
     L        (nu2 + eta + nu*dsqrt(nu2 + eta) + 
     M        two*dsqrt(nu2 + eta)*
     N        dsqrt(nu2 + zeta))*
     O        d3zeta)))/
     P        (nu2 + zeta)**f52)/(eight*(nu2 + eta)**f52) 
      d30Ks = d30Js
      d30Ms = (three*(nu + dsqrt(nu2 + zeta))*
     1        d1lambda**three - 
     2        (three*(nu2 + lambda)*d1lambda**two*
     3        d1zeta)/dsqrt(nu2 + zeta) - 
     4        (three*(nu2 + lambda)*d1lambda*
     5        ((nu2 + lambda)*d1zeta**two + 
     6        two*(nu2 + zeta)*
     7        ((zeta + nu*(nu + dsqrt(nu2 + zeta)))*
     8        d2lambda - 
     9        (nu2 + lambda)*d2zeta)))/
     A        (nu2 + zeta)**f32 + 
     B        ((nu2 + lambda)**two*
     C        (three*(lambda + nu*(nu + dsqrt(nu2 + lambda)))*
     D        d1zeta**three + 
     E        six*(nu2 + zeta)*d1zeta*
     F        ((nu2 + zeta)*d2lambda - 
     G        (lambda + 
     H        nu*(nu + dsqrt(nu2 + lambda)))*
     I        d2zeta) + 
     J        four*(nu2 + zeta)**two*
     K        ((nu2 + zeta + nu*dsqrt(nu2 + zeta) + 
     L        two*dsqrt(nu2 + lambda)*
     M        dsqrt(nu2 + zeta))*
     N        d3lambda + 
     O        (lambda + 
     P        nu*(nu + dsqrt(nu2 + lambda)))*
     Q        d3zeta)))/
     R        (nu2 + zeta)**f52)/(eight*(nu2 + lambda)**f52)
      d30Ns = (three*(nu + dsqrt(nu2 + lambda))*
     1        d1eta**three - 
     2        (three*(nu2 + eta)*d1eta**two*
     3        d1lambda)/dsqrt(nu2 + lambda)
     4        - (three*(nu2 + eta)*d1eta*
     5        ((nu2 + eta)*d1lambda**two + 
     6        two*(nu2 + lambda)*
     7        ((lambda + 
     8        nu*(nu + dsqrt(nu2 + lambda)))*
     9        d2eta - 
     A        (nu2 + eta)*d2lambda)))/
     B        (nu2 + lambda)**f32 + 
     C        ((nu2 + eta)**two*
     D        (three*(eta + nu*(nu + dsqrt(nu2 + eta)))*
     E        d1lambda**three + 
     F        six*(nu2 + lambda)*d1lambda*
     G        ((nu2 + lambda)*d2eta - 
     H        (eta + nu*(nu + dsqrt(nu2 + eta)))*
     I        d2lambda) + 
     J        four*(nu2 + lambda)**two*
     K        ((lambda + 
     L        nu*(nu + dsqrt(nu2 + lambda)))*
     M        d3eta + 
     N        (nu2 + eta + nu*dsqrt(nu2 + eta) + 
     O        two*dsqrt(nu2 + eta)*
     P        dsqrt(nu2 + lambda))*
     Q        d3lambda)))/
     R        (nu2 + lambda)**f52)/(eight*(nu2 + eta)**f52)
c
c     Calculate derivatives of the enhancement factor with respect
c     to s.
c
      L4 = L2*L2
      L6 = L4*L2
      Ms3 = Ms2*Ms
      Ns3 = Ns2*Ns
      Js4 = Js2*Js2
      Ks4 = Ks2*Ks2
c
      d3Fxhse1 = (A*(Js3*Ks4*d3zeta - 
     1    six*Ks4*zeta*d10Js**three + 
     2    six*Js*Ks4*d10Js*
     3    (d1zeta*d10Js + zeta*d20Js) - 
     4    Js**two*Ks4*(three*d2zeta*d10Js + 
     5    three*d1zeta*d20Js + 
     6    zeta*d30Js) + Js4*
     7    (Ks3*d3eta - six*eta*d10Ks**three + 
     8    six*Ks*d10Ks*(d1eta*d10Ks + eta*d20Ks) - 
     9    Ks**two*(three*d2eta*d10Ks + three*d1eta*d20Ks + 
     A    eta*d30Ks))))/(Js4*Ks4)
c
      d3Fxhse2 = (four*B*(six*d10lambda2**three - 
     1    six*lambda2*d10lambda2*d20lambda2 + 
     2    lambda2**two*d30lambda2))/(nine*L4)
c
      d3Fxhse3 = (two*C*(r24*(two + chi)*Fs*d10lambda2**three - 
     1   r18*lambda2*d10lambda2*
     2   ((two + chi)*d1Fs*d10lambda2 + 
     3   Fs*(d10chi*d10lambda2 + (two +
     4   chi)*d20lambda2)) - 
     5   L3*((two + chi)*d3Fs + three*d2Fs*d10chi + 
     6   three*d1Fs*d20chi + Fs*d30chi) + 
     7   two*L2*(three*(two + chi)*d2Fs*d10lambda2 + 
     8   three*d1Fs*(two*d10chi*d10lambda2 + 
     9   (two + chi)*d20lambda2) + 
     A   Fs*(three*d10lambda2*d20chi + three*d10chi*d20lambda2 + 
     B   (two + chi)*d30lambda2))))/(nine*L5)
c
      tmp1 = -EGs*(eight/nine)/L3
      tmp2 = (-d1EGs + (three*EGs*d10lambda2)/(lambda2))*(eight/nine)/L3
      tmp3 = (-d2EGs + (six*d1EGs*d10lambda2)/(lambda2)
     1  - (EGs*((r12*d10lambda2**two)/L2 - 
     2  (three*d20lambda2)/lambda2)))*(eight/nine)/L3
      tmp4 = (-d3EGs + (nine*d2EGs*d10lambda2)/lambda2 - 
     1  three*(d1EGs*((r12*d10lambda2**two)/L2 - 
     2  (three*d20lambda2)/lambda2)) - (EGs*
     3  ((-r60*d10lambda2**three)/L3 + 
     4  (r36*d10lambda2*d20lambda2)/L2
     5  - (three*d30lambda2)/lambda2)))*(eight/nine)/L3
      tmp5 = (one + (nine/eight)*chi + (three/eight)*chi*chi)
      tmp6 = (nine*d10chi)/eight + (three*chi*d10chi)/four
      tmp7 = (nine*d20chi)/eight + (three*(two*d10chi**two + 
     1  two*chi*d20chi))/eight
      tmp8 = (nine*d30chi)/eight + (three*(six*d10chi*
     1  d20chi + two*chi*d30chi))/eight


      d3Fxhse4 = (tmp1*tmp8+3d0*tmp2*tmp7+3d0*tmp3*tmp6+
     1           tmp4*tmp5)
c
      d3Fxhse5 = two*dlog(one - D/Ms)*d3zeta - 
     1    (six*D*d2zeta*d10Ms)/(D*Ms - Ms**two) + 
     2  (six*D*d1zeta*((D - two*Ms)*d10Ms**two + Ms*(-D + Ms)*d20Ms))/
     3  ((D - Ms)**two*Ms**two) - 
     4  (two*D*zeta*
     5  (two*(D**two - three*D*Ms + three*Ms**two)*d10Ms**three - 
     6  three*Ms*(D**two - three*D*Ms + two*Ms**two)*d10Ms*d20Ms + 
     7  (D - Ms)**two*Ms**two*d30Ms))/((D - Ms)**three*Ms3) 
c
      d3Fxhse6 = -two*dlog((A - D + Ns)/Ns)*d3eta + 
     1  (six*(A - D)*d2eta*d10Ns)/(Ns*(A - D + Ns)) + 
     2  (six*(A - D)*d1eta*((-A + D - two*Ns)*d10Ns**two + 
     3  Ns*(A - D + Ns)*d20Ns))/(Ns**two*(A - D + Ns)**two) + 
     4  (two*(A - D)*eta*(two*((A - D)**two + three*(A - D)*Ns
     5  + three*Ns**two)*d10Ns**three - 
     6  three*Ns*((A - D)**two + three*(A - D)*Ns +
     7  two*Ns**two)*d10Ns*d20Ns + 
     8  Ns**two*(A - D + Ns)**two*d30Ns))/
     9  (Ns3*(A - D + Ns)**three)
c
      d30Fxhse = d3Fxhse1 + d3Fxhse2 + d3Fxhse3 + d3Fxhse4 +
     1           d3Fxhse5 + d3Fxhse6      
c
c
c     Calculate the third derivative of Fxhse with respect to the
c     parameter nu (cam_omega/kf).
c
      nu7 = nu6*nu
      nu8 = nu7*nu
c
c     Calculate third derivatives with respect to nu necessary to
c     simplify the appearance of the code.
c
      d03chi = (three*(four*nu2 - lambda)*lambda)/(nu2 + lambda)**f72
      d03lambda2 = (lambda*(six*d01chi**three - 
     1             six*(-one + chi)*d01chi*
     2             d02chi + 
     3             (-one + chi)**two*d03chi))/
     4             (-one + chi)**four
      d03Js = (three*(-(nu*eta**three*zeta) + 
     1        nu4*(nu + dsqrt(nu2 + eta))*zeta**two - 
     2        nu*eta*zeta*
     3        (two*nu4 - two*nu*dsqrt(nu2 + eta)*zeta + 
     4        zeta**two) + 
     5        eta**two*(two*nu2*zeta*dsqrt(nu2 + zeta) + 
     6        nu4*(nu + dsqrt(nu2 + zeta)) + 
     7        zeta**two*
     8        (two*nu + dsqrt(nu2 + eta) + 
     9        dsqrt(nu2 + zeta)))))/
     A        ((nu2 + eta)**f52*(nu2 + zeta)**f52)
      d03Ks = d03Js
      d03Ms = (three*(-(nu*lambda**three*zeta) + 
     1        nu4*(nu + dsqrt(nu2 + lambda))*zeta**two - 
     2        nu*lambda*zeta*
     3        (two*nu4 - two*nu*dsqrt(nu2 + lambda)*zeta + 
     4        zeta**two) + 
     5        lambda**two*(two*nu2*zeta*
     6        dsqrt(nu2 + zeta) + 
     7        nu4*(nu + dsqrt(nu2 + zeta)) + 
     8        zeta**two*
     9        (two*nu + dsqrt(nu2 + lambda) + 
     A        dsqrt(nu2 + zeta)))))/
     B        ((nu2 + lambda)**f52*(nu2 + zeta)**f52)
      d03Ns = (three*(-(nu*eta**three*lambda) + 
     1        nu4*(nu + dsqrt(nu2 + eta))*lambda**two - 
     2        nu*eta*lambda*
     3        (two*nu4 - two*nu*dsqrt(nu2 + eta)*lambda + 
     4        lambda**two) + 
     5        eta**two*(two*nu2*lambda*
     6        dsqrt(nu2 + lambda) + 
     7        nu4*(nu + dsqrt(nu2 + lambda)) + 
     8        lambda**two*
     9        (two*nu + dsqrt(nu2 + eta) + 
     A        dsqrt(nu2 + lambda)))))/
     B        ((nu2 + eta)**f52*(nu2 + lambda)**f52)
c
      d3Fxhse1 = A*(-((zeta*(six*d01Js**three - six*Js*d01Js*d02Js + 
     1           Js**two*d03Js))/Js4) - 
     2    (eta*(six*d01Ks**three - six*Ks*d01Ks*d02Ks + 
     3         Ks**two*d03Ks))/Ks4)
c
      d3Fxhse2 = (four*B*(six*d01lambda2**three - 
     1      six*lambda2*d01lambda2*d02lambda2 + 
     2      lambda2**two*d03lambda2))/(nine*L4)
c
      d3Fxhse3 = (two*C*Fs*(r24*(two + chi)*
     1  d01lambda2**three - r18*lambda2*d01lambda2*
     2  (d01chi*d01lambda2 + (two + chi)*d02lambda2)
     3  - L3*d03chi + two*L2*(three*d01lambda2*d02chi + 
     4  three*d01chi*d02lambda2 + 
     5  (two + chi)*d03lambda2)))/(nine*L5)
c
      d3Fxhse4 = (EGs*(r20*(eight + nine*chi + three*chi**two)*
     1  d01lambda2**three - r12*lambda2*d01lambda2*
     2  (three*(three + two*chi)*d01chi*d01lambda2 + 
     3  (eight + nine*chi + three*chi**two)*d02lambda2) - 
     4  L3*(six*d01chi*d02chi + (three + two*chi)*d03chi) + 
     5  L2*(r18*d01chi**two*d01lambda2 + nine*(three + two*chi)*
     6  d01lambda2*d02chi + nine*(three + two*chi)*d01chi*
     7  d02lambda2 + (eight + nine*chi + three*chi**two)*
     8  d03lambda2)))/(three*L6)
c
      d3Fxhse5 = (-two*D*zeta*(two*(D**two - three*D*Ms + 
     1  three*Ms**two)*d01Ms**three - 
     2  three*Ms*(D**two - three*D*Ms + two*Ms**two)*d01Ms*d02Ms + 
     3  (D - Ms)**two*Ms**two*d03Ms))/((D - Ms)**three*Ms3)
c
      d3Fxhse6 = (two*(A - D)*eta*(two*((A - D)**two + 
     1  three*(A - D)*Ns + three*Ns**two)*d01Ns**three - 
     2  three*Ns*((A - D)**two + 
     3  three*(A - D)*Ns + two*Ns**two)*d01Ns*d02Ns + 
     4  Ns**two*(A - D + Ns)**two*d03Ns))/
     5  (Ns3*(A - D + Ns)**three)
c
      d03Fxhse = d3Fxhse1 + d3Fxhse2 + d3Fxhse3 + d3Fxhse4 +
     1           d3Fxhse5 + d3Fxhse6
c
c     Calculate the mixed third derivative of Fxhse (dnu ds^2)
c
c
c     Calculate mixed partial derivatives for the variables depending
c     on both s and nu.
c
      d21chi = (three*(-four*nu2 + lambda)*d1lambda**two + 
     1         two*(two*nu4 + nu2*lambda - lambda**two)*d2lambda)/
     2         (four*(nu2 + lambda)**f72)
      d21lambda2 = ((-one + chi)**two*d2lambda*d01chi + 
     1         two*(-one + chi)*d1lambda*(-two*d01chi*d10chi + 
     2         (-one + chi)*d11chi) + 
     3         lambda*(two*d01chi*(three*d10chi**two - 
     4         (-one + chi)*d20chi) + 
     5    (-one + chi)*(-four*d10chi*d11chi + (-one + chi)*d21chi)))/
     6         (-one + chi)**four
      d21Js = (-(nu*(nu2 + eta)*(nu2 + zeta)**two*d1eta**two) + 
     1        three*nu2*(nu2 + zeta)**f52*d1eta**two - 
     2        (nu2 + eta)*(nu2 + zeta)**f52*d1eta**two +
     3        three*nu*(nu2 + zeta)**three*d1eta**two - 
     4        two*nu*(nu2 + eta)**two*(nu2 + zeta)*d1eta*d1zeta - 
     5        two*nu*(nu2 + eta)*(nu2 + zeta)**two*d1eta*d1zeta + 
     6        three*nu2*(nu2 + eta)**f52*d1zeta**two + 
     7        three*nu*(nu2 + eta)**three*d1zeta**two - 
     8        nu*(nu2 + eta)**two*(nu2 + zeta)*d1zeta**two - 
     9        (nu2 + eta)**f52*(nu2 + zeta)*d1zeta**two + 
     A        two*nu*(nu2 + eta)**two*(nu2 + zeta)**two*d2eta - 
     B        two*nu2*(nu2 + eta)*(nu2 + zeta)**f52*d2eta + 
     C        two*(nu2 + eta)**two*(nu2 + zeta)**f52*d2eta - 
     D        two*nu*(nu2 + eta)*(nu2 + zeta)**three*d2eta - 
     E        two*nu2*(nu2 + eta)**f52*(nu2 + zeta)*d2zeta - 
     F        two*nu*(nu2 + eta)**three*(nu2 + zeta)*d2zeta + 
     G        two*nu*(nu2 + eta)**two*(nu2 + zeta)**two*d2zeta + 
     H        two*(nu2 + eta)**f52*(nu2 + zeta)**two*d2zeta)/
     I        (four*(nu2 + eta)**f52*(nu2 + zeta)**f52)
      d21Ks = d21Js
      d21Ms = -(nu*(nu2 + lambda)*(nu2 + zeta)**two*d1lambda**two - 
     1        three*nu2*(nu2 + zeta)**f52*d1lambda**two + 
     2        (nu2 + lambda)*(nu2 + zeta)**f52*d1lambda**two -
     3        three*nu*(nu2 + zeta)**three*d1lambda**two + 
     4      two*nu*(nu2 + lambda)**two*(nu2 + zeta)*d1lambda*d1zeta + 
     5      two*nu*(nu2 + lambda)*(nu2 + zeta)**two*d1lambda*d1zeta - 
     6        three*nu2*(nu2 + lambda)**f52*d1zeta**two -
     7        three*nu*(nu2 + lambda)**three*d1zeta**two + 
     8        nu*(nu2 + lambda)**two*(nu2 + zeta)*d1zeta**two + 
     9        (nu2 + lambda)**f52*(nu2 + zeta)*d1zeta**two - 
     A        two*nu*(nu2 + lambda)**two*(nu2 + zeta)**two*d2lambda + 
     B        two*nu2*(nu2 + lambda)*(nu2 + zeta)**f52*d2lambda - 
     C        two*(nu2 + lambda)**two*(nu2 + zeta)**f52*d2lambda + 
     D        two*nu*(nu2 + lambda)*(nu2 + zeta)**three*d2lambda + 
     E        two*nu2*(nu2 + lambda)**f52*(nu2 + zeta)*d2zeta + 
     F        two*nu*(nu2 + lambda)**three*(nu2 + zeta)*d2zeta - 
     G        two*nu*(nu2 + lambda)**two*(nu2 + zeta)**two*d2zeta - 
     H        two*(nu2 + lambda)**f52*(nu2 + zeta)**two*d2zeta)/
     I        (four*(nu2 + lambda)**f52*(nu2 + zeta)**f52) 
      d21Ns = -(nu*(nu2 + eta)*(nu2 + lambda)**two*d1eta**two -
     1        three*nu2*(nu2 + lambda)**f52*d1eta**two + 
     2        (nu2 + eta)*(nu2 + lambda)**f52*d1eta**two -
     3        three*nu*(nu2 + lambda)**three*d1eta**two + 
     4        two*nu*(nu2 + eta)**two*(nu2 + lambda)*d1eta*d1lambda + 
     5        two*nu*(nu2 + eta)*(nu2 + lambda)**two*d1eta*d1lambda - 
     6        three*nu2*(nu2 + eta)**f52*d1lambda**two -
     7        three*nu*(nu2 + eta)**three*d1lambda**two + 
     8        nu*(nu2 + eta)**two*(nu2 + lambda)*d1lambda**two + 
     9        (nu2 + eta)**f52*(nu2 + lambda)*d1lambda**two - 
     A        two*nu*(nu2 + eta)**two*(nu2 + lambda)**two*d2eta + 
     B        two*nu2*(nu2 + eta)*(nu2 + lambda)**f52*d2eta - 
     C        two*(nu2 + eta)**two*(nu2 + lambda)**f52*d2eta + 
     D        two*nu*(nu2 + eta)*(nu2 + lambda)**three*d2eta + 
     E        two*nu2*(nu2 + eta)**f52*(nu2 + lambda)*d2lambda + 
     F        two*nu*(nu2 + eta)**three*(nu2 + lambda)*d2lambda - 
     G        two*nu*(nu2 + eta)**two*(nu2 + lambda)**two*d2lambda - 
     H        two*(nu2 + eta)**f52*(nu2 + lambda)**two*d2lambda)/
     I        (four*(nu2 + eta)**f52*(nu2 + lambda)**f52)
c
      d21Fxhse1 = -((A*(six*Ks4*zeta*d01Js*d10Js**two - 
     1   two*Js*Ks4*(two*d1zeta*d01Js*d10Js + 
     2   zeta*(two*d10Js*d11Js + d01Js*d20Js)) + 
     3   Js**two*Ks4*(d2zeta*d01Js + two*d1zeta*d11Js + 
     4   zeta*d21Js) + Js4*
     5   (six*eta*d01Ks*d10Ks**two - 
     6   two*Ks*(two*d1eta*d01Ks*d10Ks + 
     7   two*eta*d10Ks*d11Ks + eta*d01Ks*d20Ks) + 
     8   Ks**two*(d2eta*d01Ks + two*d1eta*d11Ks + 
     9   eta*d21Ks))))/(Js4*Ks4))
c
      d21Fxhse2 = (four*B*(d01lambda2*(six*d10lambda2**two - 
     1   two*lambda2*d20lambda2) + 
     2   lambda2*(-four*d10lambda2*d11lambda2 + lambda2*d21lambda2)))/
     3   (nine*L4)
c
      d21Fxhse3 = (two*C*(r24*(two + chi)*Fs*d01lambda2*
     1   d10lambda2**two - six*lambda2*
     2   (two*(two + chi)*d1Fs*d01lambda2*d10lambda2 + 
     3   Fs*(d10lambda2*(d01chi*d10lambda2 + two*(two + chi)*
     4   d11lambda2) + d01lambda2*(two*d10chi*d10lambda2 + 
     5   (two + chi)*d20lambda2))) - L3*(d2Fs*d01chi + two*d1Fs*
     6   d11chi + Fs*d21chi) + two*L2*((two + chi)*d2Fs*d01lambda2 + 
     7   two*d1Fs*(d01lambda2*d10chi + d01chi*d10lambda2 + 
     8   (two + chi)*d11lambda2) + Fs*(two*d10lambda2*d11chi + 
     9   two*d10chi*d11lambda2 + d01lambda2*d20chi + d01chi*
     A   d20lambda2 + two*d21lambda2 + chi*d21lambda2))))/
     B  (nine*L5)
c
      d21Fxhse4 = (-two*lambda2**two*(lambda2*d1EGs -
     1   three*EGs*d10lambda2)*
     2   (two*d01chi*d10chi + (three + two*chi)*d11chi) + 
     3   six*(three + two*chi)*lambda2*d10chi*
     4   (-four*EGs*d01lambda2*d10lambda2 + 
     5   lambda2*(d1EGs*d01lambda2 + EGs*d11lambda2)) + 
     6   three*EGs*lambda2**two*d01lambda2*
     7   (two*d10chi**two + (three + two*chi)*d20chi) - 
     8   (three + two*chi)*lambda2*d01chi*
     9   (lambda2**two*d2EGs + r12*EGs*d10lambda2**two - 
     A   three*lambda2*(two*d1EGs*d10lambda2 + EGs*d20lambda2)) - 
     B   EGs*L3*(four*d10chi*d11chi + 
     C   two*d01chi*d20chi + (three + two*chi)*d21chi) + 
     D   (eight + nine*chi + three*chi**two)*
     E   (r20*EGs*d01lambda2*d10lambda2**two - 
     F   four*lambda2*(two*d1EGs*d01lambda2*d10lambda2 + 
     G   EGs*(two*d10lambda2*d11lambda2 + 
     H   d01lambda2*d20lambda2)) + 
     I   lambda2**two*(d2EGs*d01lambda2 + two*d1EGs*d11lambda2 + 
     J   EGs*d21lambda2)))/(three*L6)
c
      d21Fxhse5 = (-two*D*(two*D**two*zeta*d01Ms*d10Ms**two - 
     1   D*Ms*(two*D*d1zeta*d01Ms*d10Ms + 
     2   zeta*(two*D*d10Ms*d11Ms + 
     3   d01Ms*(six*d10Ms**two + D*d20Ms))) + 
     4   Ms**four*(d2zeta*d01Ms + two*d1zeta*d11Ms + 
     5   zeta*d21Ms) - two*Ms3*
     6   (D*d2zeta*d01Ms + 
     7   two*d1zeta*(d01Ms*d10Ms + D*d11Ms) + 
     8   zeta*(two*d10Ms*d11Ms + d01Ms*d20Ms + 
     9   D*d21Ms)) + Ms**two*
     A   (D**two*d2zeta*d01Ms + 
     B   two*D*d1zeta*(three*d01Ms*d10Ms + D*d11Ms) + 
     C   zeta*(three*d01Ms*(two*d10Ms**two + D*d20Ms) + 
     D   D*(six*d10Ms*d11Ms + D*d21Ms)))))/((D - Ms)**three*Ms3)
c
      d21Fxhse6 = (two*(A - D)*(Ns**two*(A - D + Ns)**two*d2eta*d01Ns + 
     1   two*(A - D)*Ns*(A - D + Ns)*d1eta*d01Ns*d10Ns - 
     2   four*Ns*(A - D + Ns)**two*d1eta*d01Ns*d10Ns + 
     3   two*Ns**two*(A - D + Ns)**two*d1eta*d11Ns + 
     4   eta*(d01Ns*(two*((A - D)**two + three*(A - D)*Ns +
     5   three*Ns**two)*d10Ns**two - 
     6   Ns*((A - D)**two + three*(A - D)*Ns + two*Ns**two)*d20Ns) + 
     7   Ns*(A - D + Ns)*(-two*(A - D + two*Ns)*d10Ns*d11Ns + 
     8   Ns*(A - D + Ns)*d21Ns))))/(Ns3*(A - D + Ns)**three)
c
      d21Fxhse = d21Fxhse1 + d21Fxhse2 + d21Fxhse3 + d21Fxhse4 +
     1           d21Fxhse5 + d21Fxhse6
c
c     Calculate the mixed third derivative of Fxhse (dnu^2 ds)
c
c
c     Calculate mixed partial derivatives for the variables depending
c     on both s and nu.
c
      d12chi = (three*nu*(-two*nu2 + three*lambda)*d1lambda)/
     1         (two*(nu2 + lambda)**f72)
      d12lambda2 = ((-one + chi)*d1lambda*(-two*d01chi**two +
     1   (-one + chi)*d02chi) + 
     2   lambda*(six*d01chi**two*d10chi - 
     3   four*(-one + chi)*d01chi*d11chi + 
     4   (-one + chi)*(-two*d02chi*d10chi + (-one + chi)*d12chi)))/
     5   (-one + chi)**four
      d12Js = (eta**three*(two*nu2 - zeta)*d1zeta + 
     1   nu2*zeta*((three*nu4 + five*nu2*zeta + two*zeta**two)*d1eta - 
     2   three*nu3*(nu + dsqrt(nu2 + eta))*d1zeta) - 
     3   eta*((nu2 + zeta)*(zeta**two + three*nu3*(nu + 
     4   dsqrt(nu2 + zeta)) + 
     5   nu*zeta*(two*nu + three*dsqrt(nu2 + zeta)))*
     6   d1eta - nu2*(three*nu4 - nu*(five*nu + 
     7   six*dsqrt(nu2 + eta))*zeta + zeta**two)*d1zeta) + 
     8   eta**two*(five*nu4*d1zeta + zeta**two*(d1eta + d1zeta) + 
     9   nu*zeta*(nu*d1eta - three*(nu + dsqrt(nu2 + eta))*d1zeta)))/
     A   (two*(nu2 + eta)**f52*(nu2 + zeta)**f52)
      d12Ks = d12Js
      d12Ms = (lambda**three*(two*nu2 - zeta)*d1zeta + 
     1   nu2*zeta*((three*nu4 + five*nu2*zeta + 
     2   two*zeta**two)*d1lambda - three*nu3*(nu + dsqrt(nu2 + lambda))*
     3   d1zeta) - lambda*((nu2 + zeta)*(zeta**two + 
     4   three*nu3*(nu + dsqrt(nu2 + zeta)) + 
     5   nu*zeta*(two*nu + three*dsqrt(nu2 + zeta)))*d1lambda - 
     6   nu2*(three*nu4 - nu*(five*nu + six*dsqrt(nu2 + lambda))*zeta + 
     7   zeta**two)*d1zeta) + lambda**two*(five*nu4*d1zeta + 
     8   zeta**two*(d1lambda + d1zeta) + 
     9   nu*zeta*(nu*d1lambda - three*(nu + dsqrt(nu2 + lambda))*
     A   d1zeta)))/(two*(nu2 + lambda)**f52*(nu2 + zeta)**f52)
      d12Ns = (eta**three*(two*nu2 - lambda)*d1lambda + nu2*lambda*
     1   ((three*nu4 + five*nu2*lambda + two*lambda**two)*d1eta - 
     2   three*nu3*(nu + dsqrt(nu2 + eta))*d1lambda) - 
     3   eta*((nu2 + lambda)*(lambda**two + 
     4   three*nu3*(nu + dsqrt(nu2 + lambda)) + 
     5   nu*lambda*(two*nu + three*dsqrt(nu2 + lambda)))*d1eta - 
     6   nu2*(three*nu4 - nu*(five*nu + six*dsqrt(nu2 + eta))*lambda + 
     7   lambda**two)*d1lambda) + eta**two*(five*nu4*d1lambda + 
     8   lambda**two*(d1eta + d1lambda)+ nu*lambda*(nu*d1eta - 
     9   three*(nu + dsqrt(nu2 + eta))*d1lambda)))/
     A   (two*(nu2 + eta)**f52*(nu2 + lambda)**f52)
c
      d12Fxhse1 = -((A*(six*Ks4*zeta*d01Js**two*d10Js - 
     1   two*Js*Ks4*(d1zeta*d01Js**two + 
     2   zeta*(d02Js*d10Js + two*d01Js*d11Js)) + 
     3   Js**two*Ks4*(d1zeta*d02Js + zeta*d12Js) + 
     4   Js4*(six*eta*d01Ks**two*d10Ks - 
     5   two*Ks*(d1eta*d01Ks**two + eta*d02Ks*d10Ks + 
     6   two*eta*d01Ks*d11Ks) + 
     7   Ks**two*(d1eta*d02Ks + eta*d12Ks))))/(Js4*Ks4))
c
      d12Fxhse2 = (four*B*(six*d01lambda2**two*d10lambda2 - 
     1   four*lambda2*d01lambda2*d11lambda2 + 
     2   lambda2*(-two*d02lambda2*d10lambda2 + lambda2*d12lambda2)))/
     3   (nine*L4)
c
      d12Fxhse3 = (-five*C*nu*(two*(nu2 + lambda)*d1Fs - 
     1   seven*Fs*d1lambda))/(three*(nu2 + lambda)**f92)
c
      d12Fxhse4 = (r20*(eight + nine*chi + three*chi**two)*EGs*
     1   d01lambda2**two*d10lambda2 - 
     2   four*lambda2*((eight + nine*chi + three*chi**two)*
     3   d1EGs*d01lambda2**two + 
     4   EGs*(three*(three + two*chi)*d01lambda2**two*d10chi + 
     5   (eight + nine*chi + three*chi**two)*d02lambda2*d10lambda2 + 
     6   two*d01lambda2*(three*(three + two*chi)*d01chi*d10lambda2 + 
     7   (eight + nine*chi + three*chi**two)*d11lambda2))) - 
     8   L3*(d1EGs*(two*d01chi**two + 
     9   (three + two*chi)*d02chi) + 
     A   EGs*(two*d02chi*d10chi + four*d01chi*d11chi + 
     B   (three + two*chi)*d12chi)) + 
     C   lambda2**two*(d1EGs*(six*(three + two*chi)*d01chi*d01lambda2 + 
     D   (eight + nine*chi + three*chi**two)*d02lambda2) + 
     E   EGs*(three*(three + two*chi)*d02lambda2*d10chi + 
     F   six*d01chi**two*d10lambda2 + nine*d02chi*d10lambda2 + 
     G   six*chi*d02chi*d10lambda2 + 
     H   r18*d01lambda2*d11chi + 
     I   r12*chi*d01lambda2*d11chi + 
     J   six*d01chi*(two*d01lambda2*d10chi + 
     K   (three + two*chi)*d11lambda2) + eight*d12lambda2 + 
     L   nine*chi*d12lambda2 + three*chi**two*d12lambda2)))/
     M   (three*L6)
c
      d12Fxhse5 = (-two*D*(two*D**two*zeta*d01Ms**two*d10Ms - 
     1   D*Ms*(D*d1zeta*d01Ms**two + zeta*(six*d01Ms**two*
     2   d10Ms + D*d02Ms*d10Ms + two*D*d01Ms*d11Ms)) + 
     3   Ms4*(d1zeta*d02Ms + zeta*d12Ms) - 
     4   two*Ms3*(d1zeta*(d01Ms**two + D*d02Ms) + 
     5   zeta*(d02Ms*d10Ms + two*d01Ms*d11Ms + D*d12Ms)) + 
     6   Ms2*(D*d1zeta*(three*d01Ms**two + D*d02Ms) + 
     7   zeta*(six*d01Ms**two*d10Ms + three*D*d02Ms*d10Ms + 
     8   six*D*d01Ms*d11Ms + D**two*d12Ms))))/
     9  ((D - Ms)**three*Ms3)
c
      d12Fxhse6 = (two*(A - D)*(Ns*(A - D + Ns)*d1eta*
     1   ((-A + D - two*Ns)*d01Ns**two + Ns*(A - D + Ns)*d02Ns) + 
     2   (A - D)*Ns*eta*(-two*d01Ns**two + (A - D + Ns)*d02Ns)*
     3   d10Ns + two*(A - D)*(A - D + Ns)*eta*d01Ns*
     4   (-two*d01Ns*d10Ns + Ns*d11Ns) + 
     5   (A - D + Ns)**two*eta*(six*d01Ns**two*d10Ns - 
     6   four*Ns*d01Ns*d11Ns + 
     7   Ns*(-two*d02Ns*d10Ns + Ns*d12Ns))))/
     8   (Ns3*(A - D + Ns)**three)
c
      d12Fxhse = d12Fxhse1 + d12Fxhse2 + d12Fxhse3 + d12Fxhse4 +
     1           d12Fxhse5 + d12Fxhse6
#endif

      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "hse08fx.F"
c
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "hse08fx.F"
#endif
