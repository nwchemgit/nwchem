      Subroutine xc_perdew91(tol_rho,fac,
     *     rho,delrho,Amat,Cmat,nq,ipol,Ec,qwght)

C     $Id: xc_perdew91.F,v 1.13 1996-09-19 19:16:22 d3h449 Exp $
      Implicit none

      double precision TOLL,EXPTOL
      double precision alpha,beta,pgamma,delta,beta10
      double precision zzz,fff,CINF
      double precision CX37,aa

      double precision ONE,ONE3,THREE
      double precision FOURTHIRDS,SEV6
      double precision FIVE3,TWO3,FIVE6
      double precision pi

      Parameter (TOLL=1.D-20,EXPTOL=40.d0)
      Parameter (alpha=0.023266D0,beta = 7.389D-6,
     ,   pgamma=8.723d0,delta=0.472d0, beta10=10000.d0*beta)
      parameter (zzz=0.001667d0,fff=0.002568d0)
      parameter(CINF=zzz+fff)
      parameter(CX37=-zzz*3.d0/7.d0)
      parameter(aa=0.09d0)

      Parameter (ONE=1.D0,ONE3=1.d0/3.d0,THREE=3.d0)
      Parameter (FOURTHIRDS=4.D0/3.D0,SEV6=7.d0/6.d0)
      parameter (FIVE3=5.d0/3.d0,TWO3=2.d0/3.d0,FIVE6=5.d0/6.d0)
      parameter (pi=3.1415926535897932385d0)

      integer ipol,nq

      logical grad_corr
c
c     logical grad_corr true if only GCA part needed
c
c
      double precision fac,tol_rho

c     Correlation energy
      double precision Ec

c--   > Charge Density 

      double precision rho(ipol*(ipol+1)/2,nq)

c--   > Charge Density Gradient

      double precision delrho(ipol,3,nq)

c--   > Quadrature Weights

      double precision qwght(nq)

c--   > Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq),Cmat(ipol,3,nq)

      double precision dp81drho(2)
      
c     
c     J.P.Perdew, J.A.Chevary, S.H.Vosko, K.A.Jackson, M.R.Pederson, 
c     D.J.~Singh and C.~Fiolhais,  Phys. Rev. B  46, 6671 (1992).     
c     
C
      integer n
      double precision expo1,expo2,dexpo1,dexpo2
      double precision TPP,nu,bb,twoab,twoab2,atwob2
      double precision rsfact,rhoval,rho13,rs,rs2,rs3,
     *     drsdrho,arho,drho
      double precision gammaval,dsqgamma,dtdgamma
      double precision anum,aden,CC,CCPRIME
      double precision egas,degasda,degasdb
      double precision g,g2,g3
      double precision gzeta,gzprime,x
      double precision dgzetadag,dgzetadbg,gggr,ggg
      double precision akf,aks,dksdrho,dkfdrho
      double precision t,t2,t3,t4,dtdrhoa,dtdrhob
      double precision A,A1p2At2,aden1,alognd,dAdrhoa,dAdrhob
      double precision func,H0,H1,dH1drhoa,dH0drhoa,
     *     dH1drhob,dH0drhob
      double precision dH1dgamma,dH0dgamma
      double precision gaa,gbb
      double precision zeta
c***************************************************************************
      gzeta(x)=0.5d0*((ONE+x)**TWO3+(ONE-x)**TWO3)
      gzprime(x)=0.5d0*TWO3*((ONE+x)**(-ONE3)-(ONE-x)**(-ONE3))
      TPP=(3.d0*pi*pi)**ONE3
      nu=(16.d0/pi)*TPP
      bb=nu*CINF
      twoab=(aa+aa)/bb
      twoab2=twoab/bb
      atwob2=1.d0/twoab2
      grad_corr=.false.
      if(abs(fac-1.d0).gt.1.e-8) grad_corr=.true.
        rsfact=(0.75d0/pi)**ONE3
        if( ipol.eq.1 )then

c       ======> SPIN-RESTRICTED <======

        do 10 n = 1,nq
          rhoval = rho(1,n)
          if(rhoval.le.tol_rho) goto 10
          arho=1.d0/rhoval
          rho13  = rhoval**ONE3
          rs=rsfact/rho13
          rs2=rs*rs
          rs3=rs2*rs
          drsdrho = -ONE3*rs*arho
          gammaval=delrho(1,1,n)*delrho(1,1,n)+
     +         delrho(1,2,n)*delrho(1,2,n)+
     +         delrho(1,3,n)*delrho(1,3,n)
          dsqgamma=dsqrt(gammaval)
C         
C         C(n)
C         
          anum=fff+alpha*rs+beta*rs2
          aden=1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
          CC=zzz+(anum*aden)
          CCPRIME=aden*((alpha+beta*(rs+rs)) -
     -         (anum*aden)*(pgamma+(rs+rs)*delta+beta10*rs2*THREE))
          CCPRIME=CCPRIME*drsdrho

          egas=0.d0
          dp81drho(1)=0.d0
C
C         call for LDA bit
C
C          call xc_pw91lda(1.d0,rhoval,dp81drho,1,1,egas,1.d0)
          call xc_p81(1.d0,rhoval,dp81drho,1,1,egas,1.d0)
C          call xc_vwn(1.d0,rhoval,dp81drho,1,1,egas,1.d0,.false.)

          degasda=arho*(dp81drho(1) - egas*arho)
          
          akf=TPP*rho13
          aks=sqrt(4.d0*akf/pi)
          dksdrho=aks*arho/6.d0
          dkfdrho=akf*arho*ONE3


          t=0.5*dsqgamma*arho/aks
          t2=t*t
          t3=t2*t
          t4=t3*t
          dtdrhoa=-t*arho*SEV6
          dtdgamma=0.5d0*t/gammaval

          dexpo1=twoab2*(egas*arho)
          dexpo2=(10.d0*(aks/akf)*t)**2
          if(dexpo1.gt.EXPTOL) then
            expo1=0.
          else
            expo1=exp(-dexpo1)
          endif
          if(dexpo2.gt.EXPTOL) then
            expo2=0.
          else
            expo2=exp(-dexpo2)
          endif

          A=twoab/(expo1-ONE)
          A1p2At2=ONE+2.d0*A*t2
          aden1=ONE/(ONE+A*t2+A*A*t4)
          alognd=(t2+A*t4)*aden1
          if(twoab*alognd.gt.TOLL) then
            H0=atwob2*log(ONE+twoab*alognd)
          else
            H0=0.d0
          endif
          H1=nu*(CC-CINF-CX37)*t2*expo2


          func=H0+H1

          if(.not.grad_corr) then
            Ec = Ec + (egas+func*rhoval)*qwght(n)*fac
          else
            Ec = Ec + (func*rhoval)*qwght(n)*fac
          endif
c
c         df/drho = func + rhoval*dfunc/drho
c
          dH1drhoa = H1*(
     *         CCPRIME/(CC-CINF-CX37) +
     *         2.d0*dtdrhoa/t -
     *         2.d0*dexpo2*(
     *         dtdrhoa/t + dksdrho/aks - dkfdrho/akf))

          dAdrhoa = A/(expo1-ONE)*expo1*
     *         twoab2*degasda
          dH0drhoa = 
     =         bb*aden1*A1p2At2*(
     (         2.d0*t*dtdrhoa + dAdrhoa*t4/A1p2At2 -
     -         alognd*(2.d0*A*t*dtdrhoa +
     +         dAdrhoa*t2))/(ONE+twoab*alognd)
          if(.not.grad_corr) then
            Amat(1,n) = Amat(1,n) +  (dp81drho(1) + func +
     +           rhoval*(dH0drhoa+dH1drhoa))*fac
          else
            Amat(1,n) = Amat(1,n) +  (func +
     +           rhoval*(dH0drhoa+dH1drhoa))*fac
          endif
c
C         df/dgamma
c

          if(dsqgamma.gt.TOLL) then

            dH1dgamma = H1/gammaval*(ONE - dexpo2 )

            dH0dgamma = bb/(ONE+twoab*alognd)*
     *           aden1*2.d0*t*dtdgamma*A1p2At2*
     *           (ONE-alognd*A)

            gaa = (dH0dgamma +dH1dgamma)*rhoval*fac*2.d0
            
            Cmat(1,1,n) = Cmat(1,1,n) + gaa*delrho(1,1,n)
            Cmat(1,2,n) = Cmat(1,2,n) + gaa*delrho(1,2,n)
            Cmat(1,3,n) = Cmat(1,3,n) + gaa*delrho(1,3,n)
          endif


   10   continue

      else

c       ======> SPIN-UNRESTRICTED <======
        do 20 n = 1,nq
          rhoval = rho(1,n)
          if(rhoval.le.tol_rho) goto 20
          arho=1.d0/rhoval
          rho13  = rhoval**ONE3
          rs=rsfact/rho13
          rs2=rs*rs
          rs3=rs2*rs
          drsdrho = -ONE3*rs*arho
          gammaval=delrho(1,1,n)*delrho(1,1,n)+
     +         delrho(2,1,n)*delrho(2,1,n)+
     +         delrho(1,2,n)*delrho(1,2,n)+
     +         delrho(2,2,n)*delrho(2,2,n)+
     +         delrho(1,3,n)*delrho(1,3,n)+
     +         delrho(2,3,n)*delrho(2,3,n)+(
     x         delrho(1,1,n)*delrho(2,1,n)+
     +         delrho(1,2,n)*delrho(2,2,n)+
     +         delrho(1,3,n)*delrho(2,3,n))*2.d0
          dsqgamma=max(dsqrt(gammaval),tol_rho)
          drho=rho(2,n)-rho(3,n)
          zeta=drho*arho

C         
C         C(n)
C         
          anum=fff+alpha*rs+beta*rs2
          aden=1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
          CC=zzz+(anum*aden)
          CCPRIME=aden*((alpha+beta*(rs+rs)) -
     -         (anum*aden)*(pgamma+(rs+rs)*delta+beta10*rs2*THREE))
          CCPRIME=CCPRIME*drsdrho

          egas=0.d0
          dp81drho(1)=0.d0
          dp81drho(2)=0.d0
C
C         call for LDA bit
C
C          call xc_pw91lda(1.d0,rho(1,n),dp81drho,1,2,egas,1.d0)
          call xc_p81(1.d0,rho(1,n),dp81drho,1,2,egas,1.d0)
C          call xc_vwn(1.d0,rhoval,dp81drho,1,1,egas,1.d0,.false.)

          degasda=arho*(dp81drho(1) - egas*arho)
          degasdb=arho*(dp81drho(2) - egas*arho)

          g=max(tol_rho,gzeta(zeta))
          g2=g*g
          g3=g2*g
          akf=TPP*rho13
          aks=sqrt(4.d0*akf/pi)
          dksdrho=aks*arho/6.d0
          dkfdrho=akf*arho*ONE3

          if(abs(1.d0-zeta).lt.tol_rho) then
            dgzetadag=0.d0
            dgzetadbg=-2.d0*(gzprime(0.d0)/g)*arho
          else
            gggr=(gzprime(zeta)/g)*arho
            dgzetadag=(1.d0-zeta)*gggr
            dgzetadbg=-(1.d0+zeta)*gggr
          endif

          t=0.5*dsqgamma*arho/(g*aks)
          t2=t*t
          t3=t2*t
          t4=t3*t
          dtdrhoa=-t*(dgzetadag+dksdrho/aks+arho)
          dtdrhob=-t*(dgzetadbg+dksdrho/aks+arho)
          dtdgamma=0.5d0*t/gammaval

          dexpo1=twoab2*(egas*arho)/g3
          dexpo2=(10.d0*g2*(aks/akf)*t)**2
          expo1=1.e-20
          if(dexpo1.lt.EXPTOL) expo1=exp(-dexpo1)
          expo2=1.e-20
          if(dexpo2.lt.EXPTOL) expo2=exp(-dexpo2)

          A=twoab/(expo1-ONE)
          A1p2At2=ONE+2.d0*A*t2
          aden1=ONE/(ONE+A*t2+A*A*t4)
          alognd=(t2+A*t4)*aden1
          H0=atwob2*g3*log(ONE+twoab*alognd)
          H1=nu*(CC-CINF-CX37)*g3*t2*expo2

          func=H0+H1

          if(.not.grad_corr) then
            Ec = Ec + (egas+func*rhoval)*qwght(n)*fac
          else
            Ec = Ec + (func*rhoval)*qwght(n)*fac
          endif
c
c         df/drho = func + rhoval*dfunc/drho
c
          dH1drhoa = H1*(
     *         CCPRIME/(CC-CINF-CX37) +
     *         THREE*dgzetadag + 2.d0*dtdrhoa/t -
     *         2.d0*dexpo2*(2.d0*dgzetadag +
     *         dtdrhoa/t + dksdrho/aks - dkfdrho/akf))

          dH1drhob = H1*(
     *         CCPRIME/(CC-CINF-CX37) +
     *         THREE*dgzetadbg + 2.d0*dtdrhob/t -
     *         2.d0*dexpo2*(2.d0*dgzetadbg +
     *         dtdrhob/t + dksdrho/aks - dkfdrho/akf))

          dAdrhoa = A/(expo1-ONE)*expo1*
     *         twoab2/g3*(-THREE*dgzetadag + degasda)

          dAdrhob = A/(expo1-ONE)*expo1*
     *         twoab2/g3*(-THREE*dgzetadbg + degasdb)

          dH0drhoa = Three*H0*dgzetadag +
     +         g3*bb*aden1*(
     (         2.d0*t*dtdrhoa*A1p2At2 + dAdrhoa*t4 -
     -         alognd*A1p2At2*(2.d0*A*t*dtdrhoa +
     +         dAdrhoa*t2))/(ONE+twoab*alognd)

          dH0drhob = Three*H0*dgzetadbg +
     +         g3*bb*aden1*(
     (         2.d0*t*dtdrhob*A1p2At2 + dAdrhob*t4 -
     -         alognd*A1p2At2*(2.d0*A*t*dtdrhob +
     +         dAdrhob*t2))/(ONE+twoab*alognd)

          if(.not.grad_corr) then
            Amat(1,n) = Amat(1,n) +  (dp81drho(1) + func +
     +           rhoval*(dH0drhoa+dH1drhoa))*fac

            Amat(2,n) = Amat(2,n) +  (dp81drho(2) + func +
     +           rhoval*(dH0drhob+dH1drhob))*fac
          else
            Amat(1,n) = Amat(1,n) +  (func +
     +           rhoval*(dH0drhoa+dH1drhoa))*fac

            Amat(2,n) = Amat(2,n) +  (func +
     +           rhoval*(dH0drhob+dH1drhob))*fac
          endif
          
c
C- (func*add)*dzetadb
c

          if(dsqgamma.gt.TOLL) then

            dH1dgamma = H1/gammaval*(
     (           ONE - dexpo2 )

            dH0dgamma = g3*bb/(ONE+twoab*alognd)*
     *           aden1*2.d0*t*dtdgamma*A1p2At2*
     *           (ONE-alognd*A)

            gaa = 2.d0*(dH0dgamma +dH1dgamma)*rhoval*fac
            gbb = gaa*fac
            
            Cmat(1,1,n) = Cmat(1,1,n) + 
     +           gaa*(delrho(1,1,n)+delrho(2,1,n))
            Cmat(1,2,n) = Cmat(1,2,n) + 
     +           gaa*(delrho(1,2,n)+delrho(2,2,n))
            Cmat(1,3,n) = Cmat(1,3,n) + 
     +           gaa*(delrho(1,3,n)+delrho(2,3,n))
            Cmat(2,1,n) = Cmat(2,1,n) + 
     +           gbb*(delrho(2,1,n)+delrho(1,1,n))
            Cmat(2,2,n) = Cmat(2,2,n) + 
     +           gbb*(delrho(2,2,n)+delrho(1,2,n))
            Cmat(2,3,n) = Cmat(2,3,n) + 
     +           gbb*(delrho(2,3,n)+delrho(1,3,n))
          endif

   20   continue
      end if
      return
      end

