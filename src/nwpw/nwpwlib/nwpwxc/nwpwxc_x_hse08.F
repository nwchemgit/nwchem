#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_hse08.F
C> The range separated HSE exchange enhancement factor
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
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
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_HSE08(cam_omega,ipol,rho,s,Fxhse,
     &           d10Fxhse,d01Fxhse)
#else
      Subroutine nwpwxc_x_HSE08_d2(cam_omega,ipol,rho,s,Fxhse,
     &           d10Fxhse,d01Fxhse,d20Fxhse,d02Fxhse,d11Fxhse)
#endif

c
      implicit none
c
c HSE evaluates the Heyd et al. Screened Coulomb
c Exchange Functional
c
c Calculates the enhancement factor
c
      integer ipol
      double precision cam_omega
      double precision rho,s,Fxhse,d10Fxhse,d01Fxhse
c
#ifdef SECOND_DERIV
      double precision  d20Fxhse,d02Fxhse,d11Fxhse
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
#ifdef SECOND_DERIV
      double precision  d2H,d2hnum,d2hden
      double precision  d2Fs
      double precision  d2zeta,d2lambda,d2eta,d2nu,d2chi,d2lambda2
      double precision  d2EGs
      double precision  d2Js,d2Ks,d2Ms,d2Ns
      double precision  d2Fxhse1,d2Fxhse2,d2Fxhse3,d2Fxhse4,d2Fxhse5
      double precision  d2Fxhse6,d2Fxhse7
      double precision  d11Fxhse1,d11Fxhse2,d11Fxhse3,d11Fxhse4
      double precision  d11Fxhse5,d11Fxhse6,d11Fxhse7
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
c
c      d1Fxhse1 = -((A*zeta*(nu/dsqrt(eta + nu2) +
c     &           nu/dsqrt(nu2 + zeta)))/((nu + dsqrt(nu2 + zeta))*
c     &           (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -
c     &           (A*zeta*(one + nu/dsqrt(nu2 + zeta)))/
c     &           ((nu + dsqrt(nu2 + zeta))**two*(dsqrt(eta + nu2) +
c     &           dsqrt(nu2 + zeta)))
c
c      d1Fxhse2 = -((A*eta*(nu/dsqrt(eta + nu2) +
c     &           nu/dsqrt(nu2 + zeta)))/((nu + dsqrt(eta + nu2))*
c     &           (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -
c     &           (A*eta*(one + nu/dsqrt(eta + nu2)))/
c     &           ((nu + dsqrt(eta + nu2))**two*(dsqrt(eta + nu2) +
c     &           dsqrt(nu2 + zeta)))
c
c      d1Fxhse3 = (four*B*(-(nu2/(lambda + nu2)**(three/two)) +
c     &           one/dsqrt(lambda + nu2)))/(nine*(lambda + nu2)*
c     &           (one + nu/dsqrt(lambda + nu2))**two) +
c     &           (eight*B*nu)/(nine*(lambda + nu2)**two*
c     &           (one + nu/dsqrt(lambda + nu2)))
c
c      d1Fxhse4 = (eight*C*Fs*(-(nu2/(lambda + nu2)**(three/two)) +
c     &           one/dsqrt(lambda + nu2))*
c     & (one + nu/(two*dsqrt(lambda + nu2))))/(nine*(lambda + nu2)**two*
c     &           (one + nu/dsqrt(lambda + nu2))**three) -
c     &           (four*C*Fs*(-nu2/(two*(lambda + nu2)**(three/two)) +
c     &        one/(two*dsqrt(lambda + nu2))))/(nine*(lambda + nu2)**two*
c     &           (one + nu/dsqrt(lambda + nu2))**two) +
c     &           (r16*C*Fs*nu*(one + nu/(two*dsqrt(lambda + nu2))))/
c     &           (nine*(lambda + nu2)**three*
c     &           (one + nu/dsqrt(lambda + nu2))**two)
c
c      d1Fxhse5 = (-eight*EGs*((-three*nu3)/(four*(lambda + nu2)**two) -
c     &           (nine*nu2)/(eight*(lambda + nu2)**(three/two)) +
c     &           (three*nu)/(four*(lambda + nu2)) +
c     & nine/(eight*dsqrt(lambda + nu2))))/(nine*(lambda + nu2)**three*
c     &           (one + nu/dsqrt(lambda + nu2))**three) +
c     &           (eight*EGs*(-(nu2/(lambda + nu2)**(three/two)) +
c     &           one/dsqrt(lambda + nu2))*
c     &           (one + (three*nu2)/(eight*(lambda + nu2)) +
c     &           (nine*nu)/(eight*dsqrt(lambda + nu2))))/
c     &           (three*(lambda + nu2)**three*
c     &           (one + nu/dsqrt(lambda + nu2))**four)
c     &   + (r16*EGs*nu*(one + (three*nu2)/(eight*(lambda + nu2)) +
c     &           (nine*nu)/(eight*dsqrt(lambda + nu2))))/
c     &           (three*(lambda + nu2)**four*
c     &           (one + nu/dsqrt(lambda + nu2))**three)
c
c      d1Fxhse6 = (two*zeta*((D*(nu/dsqrt(lambda + nu2) +
c     &           nu/dsqrt(nu2 + zeta)))/
c     &           ((nu + dsqrt(lambda + nu2))*(dsqrt(lambda + nu2) +
c     &           dsqrt(nu2 + zeta))**two) +
c     &           (D*(one + nu/dsqrt(lambda + nu2)))/
c     &           ((nu + dsqrt(lambda + nu2))**two*(dsqrt(lambda + nu2) +
c     &           dsqrt(nu2 + zeta)))))/
c     &           (one - D/((nu + dsqrt(lambda + nu2))*
c     &           (dsqrt(lambda + nu2) + dsqrt(nu2 + zeta))))
c
c      d1Fxhse7 = (-two*eta*(((-A + D)*(nu/dsqrt(eta + nu2) +
c     &           nu/dsqrt(lambda + nu2)))/
c     &           ((nu + dsqrt(lambda + nu2))*(dsqrt(eta + nu2) +
c     &           dsqrt(lambda + nu2))**two) +
c     &           ((-A + D)*(one + nu/dsqrt(lambda + nu2)))/
c     &           ((nu + dsqrt(lambda + nu2))**two*(dsqrt(eta + nu2) +
c     &           dsqrt(lambda + nu2)))))/
c     &           (one - (-A + D)/((nu + dsqrt(lambda + nu2))*
c     &           (dsqrt(eta + nu2) + dsqrt(lambda + nu2))))
c
c
      d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+
     &           d1Fxhse7
c   
#ifdef SECOND_DERIV
c
c     Calculate the second derivative of H
c
      d2hnum =  two*ha2+six*ha3*s+r12*ha4*s2+r20*ha5*s3+
     &           r30*ha6*s4 + r42*ha7*s5

      d2hden  =  two*hb2+six*hb3*s+r12*hb4*s2+r20*hb5*s3+
     &           r30*hb6*s4+r42*hb7*s5+r56*hb8*s6+r72*hb9*s7

      d2H     =  (hden*(hden*d2hnum-hnum*d2hden) -
     &           two*(hden*d1hnum-hnum*d1hden)*d1hden)/hden**three   

c
c     calculate second derivative of variables needed
c 
      d2zeta    = two*H +four*s*d1H + s2*d2H
      d2eta     = d2zeta
      d2lambda  = d2zeta

      d2chi     = -f12*(nu/(lambda+nu2)**(f32))*
     &           (-f32*d1zeta*d1zeta/(lambda+nu2) +d2zeta)

      d2lambda2 =(one-chi)*(d2lambda-d2lambda*chi+lambda*d2chi)+
     &           two*d1chi*(d1lambda-d1lambda*chi+lambda*d1chi)
      d2lambda2 = d2lambda2/(one-chi)**three

c   
c     calculate the second derivative of Fs
c   
      d2Fs = -(two/(r27*C))*(one-three*s02)/
     &        ((one+s02)**three)-d2zeta/(two*C)  

c
c     Calculate the second derivate of EGs
c
      tmp7 = -(two/five)*C*(d2Fs*lambda+two*d1Fs*d1lambda+Fs*d2lambda)
      tmp8 = -(eight/r15)*B*(d1lambda*d1lambda+lambda*d2lambda)
      tmp9 =-(r18/five)*A*lambda*(two*d1lambda*d1lambda+lambda*d2lambda)
      tmp10 = -(r14/five)*srpi*(f52*(lambda**f32)*d1lambda*d1lambda
     &          +(lambda**f52)*d2lambda)
      tmp11 = -(r42/five)*(f52*(lambda**f32)*d1lambda*d1lambda
     &          +(lambda**f52)*d2lambda)*(dsqrt(zeta)-dsqrt(eta))
      tmp12 = -(r42/five)*(lambda**f52)*d1lambda*
     &          (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))
      tmp13 = -(six/five)*(lambda**(seven/two))*
     &          (-d1zeta*d1zeta/(two*zeta**f32)+d2zeta/dsqrt(zeta)
     &          +d1eta*d1eta/(two*eta**f32)-d2eta/dsqrt(eta))

      d2EGs = tmp7+tmp8+tmp9+tmp10+tmp11+tmp12+tmp13
c
c     Calculate the second derivate of denominators needed
c
      tmp7 = two/(dsqrt(zeta+nu2)*dsqrt(eta+nu2))
      tmp8 = (dsqrt(eta+nu2)+nu)/(zeta+nu2)**f32
      tmp9 = (dsqrt(zeta+nu2)+nu)/(eta+nu2)**f32
      
      d2Js = f12*(d2zeta*(two+tmp1+tmp2) +
     &       f12*d1zeta*d1zeta*(tmp7-tmp8-tmp9)) 

      d2Ks = d2Js

      tmp10 = two/(dsqrt(zeta+nu2)*dsqrt(lambda+nu2))
      tmp11 = (dsqrt(lambda+nu2)+nu)/(zeta+nu2)**f32
      tmp12 = (dsqrt(zeta+nu2)+nu)/(lambda+nu2)**f32
      d2Ms = f12*(d2zeta*(two+tmp3+tmp4) +
     &       f12*d1zeta*d1zeta*(tmp10-tmp11-tmp12)) 

      tmp13 = two/(dsqrt(eta+nu2)*dsqrt(lambda+nu2))
      tmp14 = (dsqrt(lambda+nu2)+nu)/(eta+nu2)**f32
      tmp15 = (dsqrt(eta+nu2)+nu)/(lambda+nu2)**f32
      d2Ns = f12*(d2zeta*(two+tmp5+tmp6) +
     &       f12*d1zeta*d1zeta*(tmp13-tmp14-tmp15)) 


c
c     Calculate the second derivative of the Fxhse with respect to the
c     reduced density gradient, s.
c
      tmp1      = A*(Js*Js*d2zeta -two*Js*d1Js*d1zeta 
     &            + two*d1Js*d1Js*zeta-Js*d2Js*zeta)/
     &            (Js*Js*Js) 
      tmp2       = A*(Ks*Ks*d2zeta -two*Ks*d1Ks*d1zeta 
     &            + two*d1Ks*d1Ks*eta-Ks*d2Ks*eta)/
     &            (Ks*Ks*Ks) 
      d2Fxhse1  = tmp1 +tmp2 

      d2Fxhse2  = (four/nine)*B*(L2*d2lambda2-
     &            two*lambda2*d1lambda2*d1lambda2)/(L2*L2)

      tmp4 = d1lambda2/lambda2
      tmp5 = d2lambda2/lambda2
      d2Fxhse3  = -((four*C)/(nine*L2))*(d2Fs +
     &           six*Fs*tmp4**two - two*Fs*tmp5 -
     &           four*d1Fs*tmp4 +
     &          f12*((d2Fs*chi+two*d1Fs*d1chi+Fs*d2chi)-
     &           four*tmp4*(d1Fs*chi+Fs*d1chi) +
     &           six*chi*Fs*tmp4*tmp4)-chi*Fs*tmp5)

      tmp1  =  -d2chi -(two/three)*(d1chi*d1chi+chi*d2chi)
     &      + (six*d1chi+four*chi*d1chi
     &      -tmp4*(r32/three+r12*chi+four*chi*chi))*tmp4
     &      +((eight/three) +three*chi +chi*chi)*tmp5

      tmp2  =  -(two +(four/three)*chi)*d1chi
     &        +(two*chi*chi+(r16/three)+six*chi)*tmp4

      tmp3  =  -((eight/nine)+chi+(one/three)*chi*chi)

      d2Fxhse4  = (EGs*tmp1 + d1EGs*tmp2+ d2EGs*tmp3)/L3

      tmp1      = (one-D/Ms) 

      d2Fxhse5  = two*d2zeta*dlog(tmp1)+four*d1zeta*D*d1Ms/(Ms*Ms*tmp1)-
     &            two*zeta*(D*d1Ms/(Ms*Ms*tmp1))**two +
     &          two*zeta*D*(Ms*Ms*d2Ms-two*Ms*d1Ms*d1Ms)/(tmp1*Ms**four)

      tmp1      = (one-(D-A)/Ns) 
      d2Fxhse6  = -two*d2eta*dlog(tmp1)-
     &             four*d1eta*(D-A)*d1Ns/(Ns*Ns*tmp1)+
     &             two*eta*((D-A)*d1Ns/(Ns*Ns*tmp1))**two -
     &       two*eta*(D-A)*(Ns*Ns*d2Ns-two*Ns*d1Ns*d1Ns)/(tmp1*Ns**four)

      d20Fxhse = d2Fxhse1+d2Fxhse2+d2Fxhse3+d2Fxhse4+d2Fxhse5+d2Fxhse6
c
c     Calculate the second derivative of Fxhse with respect to the
c     parameter nu (cam_omega/kf).
c
      nu5 = nu3*nu2
      nu6 = nu5*nu
c
      d2Fxhse1 = A*zeta*((two*(one + nu/dsqrt(nu2 + zeta))*
     &            (nu/dsqrt(eta + nu2) + nu/dsqrt(nu2 + zeta)))/
     &            ((nu + dsqrt(nu2 + zeta))**two*(dsqrt(eta + nu2) +
     &            dsqrt(nu2 + zeta))**two) +
     &            ((two*(one + nu/dsqrt(nu2 + zeta))**two)/
     &            (nu + dsqrt(nu2 + zeta))**three -
     &      (-(nu2/(nu2 + zeta)**(three/two)) + one/dsqrt(nu2 + zeta))/
     &            (nu + dsqrt(nu2 + zeta))**two)/
     &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta)) +
     &    ((two*(nu/dsqrt(eta + nu2) + nu/dsqrt(nu2 + zeta))**two)/
     &    (dsqrt(eta + nu**two) + dsqrt(nu2 + zeta))**three -
     &       (-(nu2/(eta + nu2)**(three/two)) + one/dsqrt(eta + nu2) -
     &        nu2/(nu2 + zeta)**(three/two) + one/dsqrt(nu2 + zeta))/
     &        (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)/
     &             (nu + dsqrt(nu2 + zeta)))
c
      d2Fxhse2 = A*eta*((two*(one + nu/dsqrt(eta + nu2))*
     &            (nu/dsqrt(eta + nu2) + nu/dsqrt(nu2 + zeta)))/
     &            ((nu + dsqrt(eta + nu2))**two*(dsqrt(eta + nu2) +
     &            dsqrt(nu2 + zeta))**two) +
     &            ((two*(one + nu/dsqrt(eta + nu2))**two)/
     &            (nu + dsqrt(eta + nu2))**three -
     &       (-(nu2/(eta + nu2)**(three/two)) + one/dsqrt(eta + nu2))/
     &            (nu + dsqrt(eta + nu2))**two)/(dsqrt(eta + nu2) +
     &            dsqrt(nu2 + zeta)) +
     &       ((two*(nu/dsqrt(eta + nu2) + nu/dsqrt(nu2 + zeta))**two)/
     &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**three -
     &       (-(nu2/(eta + nu2)**(three/two)) + one/dsqrt(eta + nu2) -
     &         nu2/(nu2 + zeta)**(three/two) + one/dsqrt(nu2 + zeta))/
     &            (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)/
     &            (nu + dsqrt(eta + nu2)))
c
      d2Fxhse3 = (-four*B*((four*nu*(-(nu2/(lambda + nu2)**(three/two))
     &            + one/dsqrt(lambda + nu2)))/((lambda + nu2)**two*
     &            (one + nu/dsqrt(lambda + nu2))**two) +
     &            ((eight*nu2)/(lambda + nu2)**three -
     &         two/(lambda + nu2)**two)/(one + nu/dsqrt(lambda + nu2)) +
     &            ((two*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2))**two)/
     &            (one + nu/dsqrt(lambda + nu2))**three -
     &            ((three*nu3)/(lambda + nu2)**(five/two) -
     &            (three*nu)/(lambda + nu2)**(three/two))/
     &        (one + nu/dsqrt(lambda + nu2))**two)/(lambda + nu2)))/nine
c
      d2Fxhse4 =(-four*C*Fs*((-four*(-(nu2/(lambda + nu2)**(three/two))
     &            + one/dsqrt(lambda + nu2))*
     &            ((-nu2/(two*(lambda + nu2)**(three/two)) +
     &            one/(two*dsqrt(lambda + nu2)))/(lambda + nu2)**two -
     &            (four*nu*(one + nu/(two*dsqrt(lambda + nu2))))/
     &   (lambda + nu2)**three))/(one + nu/dsqrt(lambda + nu2))**three +
     &            ((one + nu/(two*dsqrt(lambda + nu2)))*
     &            ((six*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2))**two)/
     &            (one + nu/dsqrt(lambda + nu2))**four -
     &            (two*((three*nu3)/(lambda + nu2)**(five/two) -
     &            (three*nu)/(lambda + nu2)**(three/two)))/
     &            (one + nu/dsqrt(lambda + nu2))**three))/
     &            (lambda + nu2)**two +
     &            ((-eight*nu*(-nu2/(two*(lambda + nu2)**(three/two)) +
     &          one/(two*dsqrt(lambda + nu2))))/(lambda + nu2)**three +
     &   ((r24*nu2)/(lambda + nu2)**four - four/(lambda + nu2)**three)*
     &            (one + nu/(two*dsqrt(lambda + nu2))) +
     &            (-(nu/(lambda + nu2)**(three/two)) +
     &            (nu*((three*nu2)/(lambda + nu2)**(five/two) -
     &        (lambda + nu2)**(-three/two)))/two)/(lambda + nu2)**two)/
     &            (one + nu/dsqrt(lambda + nu2))**two))/nine
c
      d2Fxhse5 = (-eight*EGs*((r48*nu2)/(lambda + nu2)**five -
     &            six/(lambda + nu2)**four)*
     &            (one + (three*nu2)/(eight*(lambda + nu2)) +
     &            (nine*nu)/(eight*dsqrt(lambda + nu2))))/
     &            (nine*(one + nu/dsqrt(lambda + nu2))**three) +
     &          (r32*EGs*nu*(((-three*nu3)/(four*(lambda + nu2)**two) -
     &            (nine*nu2)/(eight*(lambda + nu2)**(three/two)) +
     &            (three*nu)/(four*(lambda + nu2)) +
     &            nine/(eight*dsqrt(lambda + nu2)))/
     &            (one + nu/dsqrt(lambda + nu2))**three -
     &            (three*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2))*
     &            (one + (three*nu2)/(eight*(lambda + nu2)) +
     &            (nine*nu)/(eight*dsqrt(lambda + nu2))))/
     &            (one + nu/dsqrt(lambda + nu2))**four))/
     &            (three*(lambda + nu2)**four) -
     &          (eight*EGs*((-six*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2))*
     &            ((-three*nu3)/(four*(lambda + nu2)**two) -
     &            (nine*nu2)/(eight*(lambda + nu2)**(three/two)) +
     &            (three*nu)/(four*(lambda + nu2)) +
     &            nine/(eight*dsqrt(lambda + nu2))))/
     &            (one + nu/dsqrt(lambda + nu2))**four +
     &            ((-three*nu2)/(lambda + nu2)**two -
     &            (nine*nu)/(four*(lambda + nu2)**(three/two)) +
     &            three/(four*(lambda + nu2)) +
     &            (three*nu2*((eight*nu2)/(lambda + nu2)**three -
     &            two/(lambda + nu2)**two))/eight +
     &            (nine*nu*((three*nu2)/(lambda + nu2)**(five/two) -
     &            (lambda + nu2)**(-three/two)))/eight)/
     &            (one + nu/dsqrt(lambda + nu2))**three +
     &            (one + (three*nu2)/(eight*(lambda + nu2)) +
     &            (nine*nu)/(eight*dsqrt(lambda + nu2)))*
     &            ((r12*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2))**two)/
     &            (one + nu/dsqrt(lambda + nu2))**five -
     &            (three*((three*nu3)/(lambda + nu2)**(five/two) -
     &            (three*nu)/(lambda + nu2)**(three/two)))/
     &            (one + nu/dsqrt(lambda + nu2))**four)))/
     &            (nine*(lambda + nu2)**three)
c
      d2Fxhse6 = two*zeta*(-(((D*(nu/dsqrt(lambda + nu2) +
     &            nu/dsqrt(nu2 + zeta)))/
     &            ((nu + dsqrt(lambda + nu2))*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))**two) +
     &            (D*(one + nu/dsqrt(lambda + nu2)))/
     &         ((nu + dsqrt(lambda + nu2))**two*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))))**two/
     &            (one - D/((nu + dsqrt(lambda + nu2))*
     &            (dsqrt(lambda + nu2) + dsqrt(nu2 + zeta))))**two) +
     &  ((-two*D*(nu/dsqrt(lambda + nu2) + nu/dsqrt(nu2 + zeta))**two)/
     &            ((nu + dsqrt(lambda + nu2))*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))**three) +
     &            (D*(-(nu2/(lambda + nu2)**(three/two)) +
     &        one/dsqrt(lambda + nu2) - nu2/(nu2 + zeta)**(three/two) +
     &            one/dsqrt(nu2 + zeta)))/
     &            ((nu + dsqrt(lambda + nu2))*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))**two) -
     &            (two*D*(one + nu/dsqrt(lambda + nu2))*
     &            (nu/dsqrt(lambda + nu2) + nu/dsqrt(nu2 + zeta)))/
     &        ((nu + dsqrt(lambda + nu2))**two*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))**two) -
     &            (two*D*(one + nu/dsqrt(lambda + nu2))**two)/
     &        ((nu + dsqrt(lambda + nu2))**three*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))) +
     &            (D*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2)))/
     &          ((nu + dsqrt(lambda + nu2))**two*(dsqrt(lambda + nu2) +
     &            dsqrt(nu2 + zeta))))/
     &            (one - D/((nu + dsqrt(lambda + nu2))*
     &            (dsqrt(lambda + nu2) + dsqrt(nu2 + zeta)))))
c
      d2Fxhse7 = -two*eta*(-((((-A + D)*(nu/dsqrt(eta + nu2) +
     &            nu/dsqrt(lambda + nu2)))/
     &            ((nu + dsqrt(lambda + nu2))*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))**two) +
     &            ((-A + D)*(one + nu/dsqrt(lambda + nu2)))/
     &            ((nu + dsqrt(lambda + nu2))**two*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))))**two/
     & (one - (-A + D)/((nu + dsqrt(lambda + nu2))*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))))**two) +
     &            ((-two*(-A + D)*(nu/dsqrt(eta + nu2) +
     &            nu/dsqrt(lambda + nu2))**two)/
     &            ((nu + dsqrt(lambda + nu2))*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))**three) -
     &            (two*(-A + D)*(one + nu/dsqrt(lambda + nu2))*
     &            (nu/dsqrt(eta + nu2) + nu/dsqrt(lambda + nu2)))/
     &            ((nu + dsqrt(lambda + nu2))**two*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))**two) +
     &            ((-A + D)*(-(nu2/(eta + nu2)**(three/two)) +
     &        one/dsqrt(eta + nu2) - nu2/(lambda + nu2)**(three/two) +
     &            one/dsqrt(lambda + nu2)))/
     &            ((nu + dsqrt(lambda + nu2))*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))**two) -
     &            (two*(-A + D)*(one + nu/dsqrt(lambda + nu2))**two)/
     &          ((nu + dsqrt(lambda + nu2))**three*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))) +
     &            ((-A + D)*(-(nu2/(lambda + nu2)**(three/two)) +
     &            one/dsqrt(lambda + nu2)))/
     &            ((nu + dsqrt(lambda + nu2))**two*(dsqrt(eta + nu2) +
     &            dsqrt(lambda + nu2))))/
     &            (one - (-A + D)/((nu + dsqrt(lambda + nu2))*
     &            (dsqrt(eta + nu2) + dsqrt(lambda + nu2)))))
c
      d02Fxhse = d2Fxhse1+d2Fxhse2+d2Fxhse3+d2Fxhse4+d2Fxhse5+d2Fxhse6+
     &           d2Fxhse7
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
      d11Fxhse1 = (A*(nu*zeta**three*d1eta - two*nu2*(nu2 + eta)*
     &           (eta + two*nu*(nu + dsqrt(nu2 + eta)))*
     &           (nu + dsqrt(nu2 + zeta))*d1zeta +
     &           nu*(eta + two*nu*(nu + dsqrt(nu2 + eta)))*zeta*
     &        (nu*(nu + dsqrt(nu2 + zeta))*d1eta - (nu2 + eta)*d1zeta) +
     &           zeta**two*((eta*(nu + dsqrt(nu2 + zeta)) +
     &          nu*(three*nu2 + two*dsqrt(nu2 + eta)*dsqrt(nu2 + zeta) +
     &           two*nu*(dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))))*
     &           d1eta + (nu2 + eta)*(nu + dsqrt(nu2 + eta))*d1zeta)))/
     &           (two*(nu2 + eta)**(f32)*(nu2 + zeta)**(f32)*
     &            (nu + dsqrt(nu2 + zeta))**two*
     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))**two)
c
      d11Fxhse2 = (A*((-nu2 - zeta)*
     &            (-(eta**two*(nu + dsqrt(nu2 + zeta))) +
     &            nu*eta*(zeta + two*nu*(nu + dsqrt(nu2 + zeta))) +
     &            two*nu2*(nu + dsqrt(nu2 + eta))*(zeta +
     &            two*nu*(nu + dsqrt(nu2 + zeta))))*d1eta +
     &            eta*(nu2 + eta)*(nu*eta + (nu + dsqrt(nu2 + eta))*
     &            (zeta + two*nu*(nu + dsqrt(nu2 + zeta))))*d1zeta))/
     &            (two*(nu2 + eta)**(f32)*(nu + dsqrt(nu2 + eta))**two*
     &            (nu2 + zeta)**(f32)*
     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))**two)
c
      d11Fxhse3 = (-two*B*(two*nu2*(nu + dsqrt(nu2 + lambda)) +
     &            lambda*(two*nu + dsqrt(nu2 + lambda)))*d1lambda)/
     &            (three*(nu2 + lambda)**three*
     &            (nu + dsqrt(nu2 + lambda))**two)
c
      d11Fxhse4 = (C*(two*(nu2 + lambda)*d1Fs - five*Fs*d1lambda))/
     &            (three*(nu2 + lambda)**(f72))
c
      d11Fxhse5 = (five*(eight*nu4*(nu + dsqrt(nu2 + lambda)) +
     &         lambda**two*(four*nu + dsqrt(nu2 + lambda)) +
     &            four*nu2*lambda*(three*nu + two*dsqrt(nu2 + lambda)))*
     &            (two*(nu2 + lambda)*d1EGs - seven*EGs*d1lambda))/
     &       (six*(nu2 + lambda)**five*(nu + dsqrt(nu2 + lambda))**four)
c
      d11Fxhse6 = (D*((-nu - two*dsqrt(nu2 + lambda))*zeta**three*
     &            d1lambda +
     &            two*nu2*(nu2 + lambda)*(-D + lambda + two*nu*
     &            (nu + dsqrt(nu2 + lambda)))*(nu + dsqrt(nu2 + zeta))*
     &            d1zeta + zeta**two*((D*(nu + dsqrt(nu2 + zeta)) -
     &            three*lambda*(nu + dsqrt(nu2 + zeta)) -
     &            nu*(five*nu2 + six*nu*dsqrt(nu2 + lambda) +
     &            four*nu*dsqrt(nu2 + zeta) + four*dsqrt(nu2 + lambda)*
     &            dsqrt(nu2 + zeta)))*
     &   d1lambda + (nu2 + lambda)*(nu + dsqrt(nu2 + lambda))*d1zeta) +
     &            zeta*(-(nu2*(-D + three*lambda + four*nu*
     & (nu + dsqrt(nu2 + lambda)))*(nu + dsqrt(nu2 + zeta))*d1lambda) +
     &            (nu2 + lambda)*(two*nu*(nu + dsqrt(nu2 + lambda))*
     &            (two*nu + dsqrt(nu2 + zeta)) -
     &            D*(nu + two*dsqrt(nu2 + zeta)) +
     &            lambda*(nu + two*dsqrt(nu2 + zeta)))*d1zeta)))/
     &            ((nu2 + lambda)**(f32)*(nu2 + zeta)**(f32)*
     &            (-D + lambda + (nu + dsqrt(nu2 + lambda))*
     &            (nu + dsqrt(nu2 + zeta)))**two)
c
      d11Fxhse7 = ((A - D)*((nu2 + lambda)*(eta**two*
     &            (nu + dsqrt(nu2 + lambda)) +
     &            two*nu2*(nu + dsqrt(nu2 + eta))*
     &        (A - D + two*nu2 + lambda + two*nu*dsqrt(nu2 + lambda)) +
     &           eta*(A*nu - D*nu + four*nu3 + two*A*dsqrt(nu2 + eta) -
     &            two*D*dsqrt(nu2 + eta) + two*nu2*dsqrt(nu2 + eta) +
     &            (nu + two*dsqrt(nu2 + eta))*lambda +
     &            four*nu2*dsqrt(nu2 + lambda) +
     &            two*nu*dsqrt(nu2 + eta)*dsqrt(nu2 + lambda)))*d1eta -
     &            eta*(nu2 + eta)*
     &            (eta*(nu + two*dsqrt(nu2 + lambda)) +
     &            (nu + dsqrt(nu2 + eta))*
     &            (A - D + four*nu2 + three*lambda +
     &            four*nu*dsqrt(nu2 + lambda)))*d1lambda))/
     &            ((nu2 + eta)**(f32)*(nu2 + lambda)**(f32)*
     &            (A - D + nu2 + nu*dsqrt(nu2 + eta) +
     &            lambda + nu*dsqrt(nu2 + lambda) +
     &            dsqrt(nu2 + eta)*dsqrt(nu2 + lambda))**two)
c
c
c      d11Fxhse1 = (A*zeta*(nu/dsqrt(nu2 + eta) +
c     &            nu/dsqrt(nu2 + zeta))*d1zeta)/
c     &      (two*dsqrt(nu2 + zeta)*(nu + dsqrt(nu2 + zeta))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))**two) -
c     &   (A*(nu/dsqrt(nu2 + eta) + nu/dsqrt(nu2 + zeta))*d1zeta)/
c     &            ((nu + dsqrt(nu2 + zeta))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**two) +
c     &            (A*zeta*(one + nu/dsqrt(nu2 + zeta))*d1zeta)/
c     &            (dsqrt(nu2 + zeta)*(nu + dsqrt(nu2 + zeta))**three*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))) +
c     &            (A*nu*zeta*d1zeta)/
c     & (two*(nu2 + zeta)**(three/two)*(nu + dsqrt(nu2 + zeta))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))) -
c     &            (A*(one + nu/dsqrt(nu2 + zeta))*d1zeta)/
c     &            ((nu + dsqrt(nu2 + zeta))**two*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))) -
c     &            (A*zeta*(-(nu*d1eta)/(two*(nu2 + eta)**(three/two)) -
c     &            (nu*d1zeta)/(two*(nu2 + zeta)**(three/two))))/
c     &            ((nu + dsqrt(nu2 + zeta))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**two) +
c     &        (two*A*zeta*(nu/dsqrt(nu2 + eta) + nu/dsqrt(nu2 + zeta))*
c     &            (d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &            ((nu + dsqrt(nu2 + zeta))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**three) +
c     &            (A*zeta*(one + nu/dsqrt(nu2 + zeta))*
c     &            (d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &            ((nu + dsqrt(nu2 + zeta))**two*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**two)
c
c      d11Fxhse2 = (A*eta*(nu/dsqrt(nu2 + eta) +
c     &            nu/dsqrt(nu2 + zeta))*d1eta)/
c     &            (two*dsqrt(nu2 + eta)*(nu + dsqrt(nu2 + eta))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))**two) -
c     &        (A*(nu/dsqrt(nu2 + eta) + nu/dsqrt(nu2 + zeta))*d1eta)/
c     &            ((nu + dsqrt(nu2 + eta))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**two) +
c     &            (A*eta*(one + nu/dsqrt(nu2 + eta))*d1eta)/
c     &            (dsqrt(nu2 + eta)*(nu + dsqrt(nu2 + eta))**three*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))) +
c     &            (A*nu*eta*d1eta)/
c     &    (two*(nu2 + eta)**(three/two)*(nu + dsqrt(nu2 + eta))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + zeta))) -
c     &            (A*(one + nu/dsqrt(nu2 + eta))*d1eta)/
c     &            ((nu + dsqrt(nu2 + eta))**two*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))) -
c     &            (A*eta*(-(nu*d1eta)/(two*(nu2 + eta)**(three/two)) -
c     &            (nu*d1zeta)/(two*(nu2 + zeta)**(three/two))))/
c     &            ((nu + dsqrt(nu2 + eta))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**two) +
c     &        (two*A*eta*(nu/dsqrt(nu2 + eta) + nu/dsqrt(nu2 + zeta))*
c     &            (d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &            ((nu + dsqrt(nu2 + eta))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**three) +
c     &            (A*eta*(one + nu/dsqrt(nu2 + eta))*
c     &            (d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &            ((nu + dsqrt(nu2 + eta))**two*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + zeta))**two)
c
c      d11Fxhse3 = (four*B*nu*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*d1lambda)/
c     &            (nine*(nu2 + lambda)**(five/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) +
c     &        (four*B*nu2*d1lambda)/(nine*(nu2 + lambda)**(seven/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) -
c     &            (four*B*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*d1lambda)/
c     &            (nine*(nu2 + lambda)**two*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) -
c     &            (r16*B*nu*d1lambda)/(nine*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))) +
c     &            (four*B*((three*nu2*d1lambda)/
c     &            (two*(nu2 + lambda)**(five/two)) -
c     &            d1lambda/(two*(nu2 + lambda)**(three/two))))/
c     &       (nine*(nu2 + lambda)*(one + nu/dsqrt(nu2 + lambda))**two)
c
c      d11Fxhse4 = (eight*C*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*
c     &            (one + nu/(two*dsqrt(nu2 + lambda)))*d1Fs)/
c     &            (nine*(nu2 + lambda)**two*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &            (four*C*(-nu2/(two*(nu2 + lambda)**(three/two)) +
c     &            one/(two*dsqrt(nu2 + lambda)))*d1Fs)/
c     &            (nine*(nu2 + lambda)**two*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) +
c     &            (r16*C*nu*(one + nu/(two*dsqrt(nu2 + lambda)))*d1Fs)/
c     &            (nine*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) +
c     &            (four*C*nu*Fs*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*(one + nu/
c     &            (two*dsqrt(nu2 + lambda)))*
c     &            d1lambda)/(three*(nu2 + lambda)**(seven/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**four) -
c     &          (four*C*nu*Fs*(-nu2/(two*(nu2 + lambda)**(three/two)) +
c     &            one/(two*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (nine*(nu2 + lambda)**(seven/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &            (two*C*nu*Fs*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*d1lambda)/
c     &            (nine*(nu2 + lambda)**(seven/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) +
c     &    (r16*C*nu2*Fs*(one + nu/(two*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (nine*(nu2 + lambda)**(nine/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &            (r16*C*Fs*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &   one/dsqrt(nu2 + lambda))*(one + nu/(two*dsqrt(nu2 + lambda)))*
c     &            d1lambda)/(nine*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &       (four*C*nu2*Fs*d1lambda)/(nine*(nu2 + lambda)**(nine/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) +
c     &            (eight*C*Fs*(-nu2/(two*(nu2 + lambda)**(three/two)) +
c     &            one/(two*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (nine*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) -
c     &     (r16*C*nu*Fs*(one + nu/(two*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (three*(nu2 + lambda)**four*
c     &            (one + nu/dsqrt(nu2 + lambda))**two) +
c     &            (eight*C*Fs*(one + nu/(two*dsqrt(nu2 + lambda)))*
c     &         ((three*nu2*d1lambda)/(two*(nu2 + lambda)**(five/two)) -
c     &            d1lambda/(two*(nu2 + lambda)**(three/two))))/
c     &            (nine*(nu2 + lambda)**two*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &            (four*C*Fs*((three*nu2*d1lambda)/
c     &            (four*(nu2 + lambda)**(five/two)) -
c     &            d1lambda/(four*(nu2 + lambda)**(three/two))))/
c     &            (nine*(nu2 + lambda)**two*
c     &            (one + nu/dsqrt(nu2 + lambda))**two)
c
c      d11Fxhse5 = (-eight*((-three*nu3)/(four*(nu2 + lambda)**two) -
c     &            (nine*nu2)/(eight*(nu2 + lambda)**(three/two)) +
c     &            (three*nu)/(four*(nu2 + lambda)) +
c     &            nine/(eight*dsqrt(nu2 + lambda)))*d1EGs)/
c     &            (nine*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) +
c     &            (eight*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*
c     &            (one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*d1EGs)/
c     &            (three*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**four) +
c     &            (r16*nu*(one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*d1EGs)/
c     &            (three*(nu2 + lambda)**four*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &        (four*nu*EGs*((-three*nu3)/(four*(nu2 + lambda)**two) -
c     &            (nine*nu2)/(eight*(nu2 + lambda)**(three/two)) +
c     &            (three*nu)/(four*(nu2 + lambda)) +
c     &            nine/(eight*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (three*(nu2 + lambda)**(nine/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**four) +
c     &            (eight*EGs*((-three*nu3)/(four*(nu2 + lambda)**two) -
c     &            (nine*nu2)/(eight*(nu2 + lambda)**(three/two)) +
c     &            (three*nu)/(four*(nu2 + lambda)) +
c     &            nine/(eight*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (three*(nu2 + lambda)**four*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) +
c     &            (r16*nu*EGs*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*
c     &            (one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (three*(nu2 + lambda)**(nine/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**five) +
c     &       (eight*nu2*EGs*(one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            ((nu2 + lambda)**(r11/two)*
c     &            (one + nu/dsqrt(nu2 + lambda))**four) -
c     &            (eight*EGs*(-(nu2/(nu2 + lambda)**(three/two)) +
c     &            one/dsqrt(nu2 + lambda))*
c     &            (one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*d1lambda)/
c     &    ((nu2 + lambda)**four*(one + nu/dsqrt(nu2 + lambda))**four) -
c     &          (r64*nu*EGs*(one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*d1lambda)/
c     &            (three*(nu2 + lambda)**five*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) -
c     &            (eight*EGs*((three*nu3*d1lambda)/
c     &            (two*(nu2 + lambda)**three) +
c     &            (r27*nu2*d1lambda)/(r16*(nu2 + lambda)**(five/two)) -
c     &            (three*nu*d1lambda)/(four*(nu2 + lambda)**two) -
c     &            (nine*d1lambda)/(r16*(nu2 + lambda)**(three/two))))/
c     &            (nine*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**three) +
c     &           (eight*EGs*(one + (three*nu2)/(eight*(nu2 + lambda)) +
c     &            (nine*nu)/(eight*dsqrt(nu2 + lambda)))*
c     &         ((three*nu2*d1lambda)/(two*(nu2 + lambda)**(five/two)) -
c     &            d1lambda/(two*(nu2 + lambda)**(three/two))))/
c     &            (three*(nu2 + lambda)**three*
c     &            (one + nu/dsqrt(nu2 + lambda))**four) +
c     &            (r16*nu*EGs*((-three*nu2*d1lambda)/
c     &            (eight*(nu2 + lambda)**two) -
c     &          (nine*nu*d1lambda)/(r16*(nu2 + lambda)**(three/two))))/
c     &            (three*(nu2 + lambda)**four*
c     &            (one + nu/dsqrt(nu2 + lambda))**three)
c
c      d11Fxhse6 = (two*((D*(nu/dsqrt(nu2 + lambda) +
c     &            nu/dsqrt(nu2 + zeta)))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))**two) +
c     &            (D*(one + nu/dsqrt(nu2 + lambda)))/
c     &       ((nu + dsqrt(nu2 + lambda))**two*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))))*
c     &            d1zeta)/(one - D/((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + lambda) + dsqrt(nu2 + zeta)))) +
c     &            (two*zeta*(-(D*(nu/dsqrt(nu2 + lambda) +
c     &            nu/dsqrt(nu2 + zeta))*d1lambda)/
c     &      (two*dsqrt(nu2 + lambda)*(nu + dsqrt(nu2 + lambda))**two*
c     &            (dsqrt(nu2 + lambda) + dsqrt(nu2 + zeta))**two) -
c     &            (D*(one + nu/dsqrt(nu2 + lambda))*d1lambda)/
c     &        (dsqrt(nu2 + lambda)*(nu + dsqrt(nu2 + lambda))**three*
c     &            (dsqrt(nu2 + lambda) + dsqrt(nu2 + zeta))) -
c     &            (D*nu*d1lambda)/
c     &            (two*(nu2 + lambda)**(three/two)*
c     &           (nu + dsqrt(nu2 + lambda))**two*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))) +
c     &           (D*(-(nu*d1lambda)/(two*(nu2 + lambda)**(three/two)) -
c     &            (nu*d1zeta)/(two*(nu2 + zeta)**(three/two))))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))**two) -
c     &          (two*D*(nu/dsqrt(nu2 + lambda) + nu/dsqrt(nu2 + zeta))*
c     &            (d1lambda/(two*dsqrt(nu2 + lambda)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))**three) -
c     &            (D*(one + nu/dsqrt(nu2 + lambda))*
c     &            (d1lambda/(two*dsqrt(nu2 + lambda)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &          ((nu + dsqrt(nu2 + lambda))**two*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))**two)))/
c     &            (one - D/((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + lambda) + dsqrt(nu2 + zeta)))) -
c     &            (two*zeta*((D*(nu/dsqrt(nu2 + lambda) +
c     &            nu/dsqrt(nu2 + zeta)))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))**two) +
c     &            (D*(one + nu/dsqrt(nu2 + lambda)))/
c     &          ((nu + dsqrt(nu2 + lambda))**two*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))))*
c     &            ((D*d1lambda)/
c     &        (two*dsqrt(nu2 + lambda)*(nu + dsqrt(nu2 + lambda))**two*
c     &            (dsqrt(nu2 + lambda) + dsqrt(nu2 + zeta))) +
c     &            (D*(d1lambda/(two*dsqrt(nu2 + lambda)) +
c     &            d1zeta/(two*dsqrt(nu2 + zeta))))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + lambda) +
c     &            dsqrt(nu2 + zeta))**two)))/
c     &            (one - D/((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + lambda) + dsqrt(nu2 + zeta))))**two
c
c      d11Fxhse7 = (-two*(((-A + D)*(nu/dsqrt(nu2 + eta) +
c     &            nu/dsqrt(nu2 + lambda)))/
c     &            ((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))**two) +
c     &            ((-A + D)*(one + nu/dsqrt(nu2 + lambda)))/
c     &            ((nu + dsqrt(nu2 + lambda))**two*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + lambda))))*d1eta)/(one - (-A + D)/
c     &            ((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda)))) -
c     &            (two*eta*(-((-A + D)*(nu/dsqrt(nu2 + eta) +
c     &            nu/dsqrt(nu2 + lambda))*d1lambda)/
c     &       (two*dsqrt(nu2 + lambda)*(nu + dsqrt(nu2 + lambda))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))**two) -
c     &            ((-A + D)*(one + nu/dsqrt(nu2 + lambda))*d1lambda)/
c     &         (dsqrt(nu2 + lambda)*(nu + dsqrt(nu2 + lambda))**three*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))) -
c     &            ((-A + D)*nu*d1lambda)/
c     &            (two*(nu2 + lambda)**(three/two)*
c     &            (nu + dsqrt(nu2 + lambda))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))) +
c     &         ((-A + D)*(-(nu*d1eta)/(two*(nu2 + eta)**(three/two)) -
c     &            (nu*d1lambda)/(two*(nu2 + lambda)**(three/two))))/
c     &            ((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))**two) -
c     &            (two*(-A + D)*(nu/dsqrt(nu2 + eta) +
c     &            nu/dsqrt(nu2 + lambda))*
c     &            (d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1lambda/(two*dsqrt(nu2 + lambda))))/
c     &            ((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))**three) -
c     &            ((-A + D)*(one + nu/dsqrt(nu2 + lambda))*
c     &            (d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1lambda/(two*dsqrt(nu2 + lambda))))/
c     &            ((nu + dsqrt(nu2 + lambda))**two*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + lambda))**two)))/
c     &            (one - (-A + D)/((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda)))) +
c     &            (two*eta*(((-A + D)*(nu/dsqrt(nu2 + eta) +
c     &            nu/dsqrt(nu2 + lambda)))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + lambda))**two) +
c     &            ((-A + D)*(one + nu/dsqrt(nu2 + lambda)))/
c     &            ((nu + dsqrt(nu2 + lambda))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))))*
c     &            (((-A + D)*d1lambda)/
c     &       (two*dsqrt(nu2 + lambda)*(nu + dsqrt(nu2 + lambda))**two*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))) +
c     &            ((-A + D)*(d1eta/(two*dsqrt(nu2 + eta)) +
c     &            d1lambda/(two*dsqrt(nu2 + lambda))))/
c     &            ((nu + dsqrt(nu2 + lambda))*(dsqrt(nu2 + eta) +
c     &            dsqrt(nu2 + lambda))**two)))/
c     &            (one - (-A + D)/((nu + dsqrt(nu2 + lambda))*
c     &            (dsqrt(nu2 + eta) + dsqrt(nu2 + lambda))))**two
c
c
      d11Fxhse = d11Fxhse1+d11Fxhse2+d11Fxhse3+d11Fxhse4+d11Fxhse5+
     &           d11Fxhse6+d11Fxhse7

#endif

      
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_x_hse08.F"
#endif
c $Id$
C> @}
