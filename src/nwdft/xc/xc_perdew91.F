      Subroutine xc_perdew91(rho,delrho,gamma,Amat,Cmat,nq,ipol,Ec,
     *     qwght)

C     $Id: xc_perdew91.F,v 1.1 1995-05-01 08:08:07 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Parameter (TOLL=1.E-40,EXPTOL=80.d0)
      Parameter (alpha=0.023266D0,beta = 7.389D-6,
     ,   pgamma=8.723d0,delta=0.472d0, beta10=10000.d0*beta)
      parameter (zzz=0.001667d0,fff=0.002568d0)
      parameter(CINF=zzz+fff)
      parameter(CX37=-zzz*3.d0/7.d0)
      parameter(akf_par=1.91915829d0)
      parameter(aa=0.09d0)

      Parameter (ONE=1.D0,ONE3=1.d0/3.d0,THREE=3.d0)
      Parameter (FOURTHIRDS=4.D0/3.D0,SEV6=7.d0/6.d0)
      parameter (FIVE3=5.d0/3.d0,TWO3=2.d0/3.d0,FIVE6=5.d0/6.d0)
      parameter (pi=3.1415926535897932385d0)

c--   > Charge Density 

      Dimension rho(ipol*(ipol+1)/2,nq)

c--   > Charge Density Gradient

      Dimension delrho(ipol,3,nq)

c--   > Charge Density Gradient Dot Products

      Dimension gamma(ipol*(ipol+1)/2,nq)

c--   > Quadrature Weights

      Dimension qwght(nq)

c--   > Sampling Matrices for the XC Potential & Energy

      Dimension Amat(ipol,nq),Cmat(ipol,3,nq)
c     
c     
c     Perdew, Fiolhais, Vosko, ..
c     
c***************************************************************************
      gzeta(x)=0.5d0*((ONE+x)**TWO3+(ONE-x)**TWO3)
      gzprime(x)=0.5d0*TWO3*((ONE+x)**(-ONE3)-(ONE-x)**(-ONE3))
      nu=(16.d0/pi)*(3.d0*pi*pi)**ONE3
      bb=nu*CINF
      twoab=(aa+aa)/bb
      twoab2=twoab/bb
      atwob2=1.d0/twoab2
c       
c       call to Perdew81 for LDA bit 
c       
        rsfact=(0.75d0/pi)**ONE3
        if( ipol.eq.1 )then

c       ======> SPIN-RESTRICTED <======

        do 10 n = 1,nq

          rhoval = rho(1,n)
          rho13  = rhoval**ONE3
          rho43  = rhoval*rho13
          rho76  = rhoval**SEV6 
          arho=1.d0/rhoval
          rs=rsfact/rho13
          rs2=rs*rs
          rs3=rs2*rs
          drsdrho = -ONE3*rs*arho
          gammaval=gamma(1,n)
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

            
c          AAadd = ?

          Amat(1,n) = Amat(1,n) + AAadd

          if(dsqgamma.gt.0.d0) then

C            g =?
            
            Cmat(1,1,n) = Cmat(1,1,n) + g*delrho(1,1,n)
            Cmat(1,2,n) = Cmat(1,2,n) + g*delrho(1,2,n)
            Cmat(1,3,n) = Cmat(1,3,n) + g*delrho(1,3,n)
          endif
          Ec      = Ec + func*qwght(n)
          
   10   continue

      else

c       ======> SPIN-UNRESTRICTED <======

        do 20 n = 1,nq
          rhoval = rho(1,n)
          arho=1.d0/rhoval
          rho13  = rhoval**ONE3
          rho43  = rhoval*rho13
          rho76  = rhoval**SEV6 
          rs=rsfact/rho13
          rs2=rs*rs
          rs3=rs2*rs
          drsdrho = -ONE3*rs*arho
          gammaval=gamma(1,n)+gamma(2,n)+2.d0*gamma(3,n)
          dsqgamma=dsqrt(gammaval)
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
C
C         call for LDA bit
C
          call xc_p81(rhoval,Amat(1,n),1,ipol,egas,1.d0)

          g=gzeta(zeta)
          g2=g*g
          g3=g2*g
          akf=akf_par/rs
          aks=sqrt(4.d0*akf/pi)

          adg=gzprime(zeta)
          dzetada=(1.d0-zeta)*arho*adg
          dzetadb=-(1.d0+zeta)*arho*adg

          t=0.5*dsqgamma*arho/(g*aks)
          t2=t*t
          t3=t2*t
          t4=t3*t

          A=twoab/(exp(-twoab2*egas*arho/g3)-ONE)
          H0=atwob2*g3*log(ONE+twoab*(t2+A*t4)/
     *         (1.d0+A*t2+A*A*t4))
          H1=nu*(CC-CINF-CX37)*g3*t2*
     *         exp(-(10.d0*g2*(aks/akf)*t)**2)

          func=H0+H1

          Ec = Ec + (egas+func*rhoval)*qwght(n)

          dphidrho =  -phi * ( SEV6/rhoval  +   (CCPRIME/CC) ) 
c
c         df/drho = func + rhoval*dfunc/drho
c
          AAadd = (  (CCPRIME/CC) -
     -      (dphidrho +   FOURTHIRDS/rhoval) )*func
C
C         wrong! put g and gzprime
C
          Amat(1,n) = Amat(1,n) + AAadd + func 
c
C derivatives of H0+H1 times rhoval
C- (func*add)*dzetada
c
          Amat(2,n) = Amat(2,n) + AAadd + func 
          
c
C- (func*add)*dzetadb
c
          if(dsqgamma.gt.TOLL) then

            gaa = 2.d0*func/gammaval*(1.d0-phi*0.5d0)
            gbb = gaa
            gab = gaa
            
            Cmat(1,1,n) = Cmat(1,1,n) + gaa*delrho(1,1,n)+
     +           gab*delrho(2,1,n)
            Cmat(1,2,n) = Cmat(1,2,n) + gaa*delrho(1,2,n)+
     +           gab*delrho(2,2,n)
            Cmat(1,3,n) = Cmat(1,3,n) + gaa*delrho(1,3,n)+
     +           gab*delrho(2,3,n)
            Cmat(2,1,n) = Cmat(2,1,n) + gbb*delrho(2,1,n)+
     +           gab*delrho(1,1,n)
            Cmat(2,2,n) = Cmat(2,2,n) + gbb*delrho(2,2,n)+
     +           gab*delrho(1,2,n)
            Cmat(2,3,n) = Cmat(2,3,n) + gbb*delrho(2,3,n)+
     +           gab*delrho(1,3,n)
          endif

   20   continue

      end if

      end

