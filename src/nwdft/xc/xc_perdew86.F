      Subroutine xc_perdew86(fac,rho,delrho,Amat,Cmat,
     *     nq,ipol,Ec,qwght)

C     $Id: xc_perdew86.F,v 1.7 1996-09-19 19:16:21 d3h449 Exp $
      implicit none

      double precision fac ! [input]
      integer nq,ipol      ! [input]
      double precision Ec  ! [input/output]

c
c--   > Charge Density 
c

      double precision rho(ipol*(ipol+1)/2,nq)

c--   > Charge Density Gradient

      double precision delrho(ipol,3,nq)

c--   > Quadrature Weights

      double precision qwght(nq)

c--   > Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq),Cmat(ipol,3,nq)
      double precision TOLL,EXPTOL,
     ,     alpha,beta,pgamma,delta,beta10,
     ,     ftilde,zzz,fff,pfff,CINF,
     ,     ONE,ONETHIRD,THREE,FOURTHIRDS,SEV6,
     ,     FIVE3,TWO3,FIVE6,pi

      Parameter (TOLL=1.E-40,EXPTOL=80.d0)
      Parameter (alpha=0.023266D0,beta = 7.389D-6,
     ,   pgamma=8.723d0,delta=0.472d0, beta10=10000.d0*beta)
      parameter (ftilde=0.11d0,zzz=0.001667d0,fff=0.002568d0)
      parameter(pfff=1.745d0,CINF=zzz+fff)

      Parameter (ONE=1.D0,ONETHIRD=1.d0/3.d0,THREE=3.d0)
      Parameter (FOURTHIRDS=4.D0/3.D0,SEV6=7.d0/6.d0)
      parameter (FIVE3=5.d0/3.d0,TWO3=2.d0/3.d0,FIVE6=5.d0/6.d0)
      parameter (pi=3.1415926535897932385d0)
c     
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact,rs,rs2,rs3
      double precision rhoval,rho13,rho43,rho76,arho
      double precision drsdrho,gammaval,dsqgamma
      double precision anum,aden,CC,CCPRIME,expfac,phi,func,dphidrho
      double precision AAadd,g,zeta,add,adp,dzetada,dzetadb
      double precision gaa,gab,gbb
      double precision dzeta,dprime
      double precision x
      dzeta(x)=2.d0**ONETHIRD*
     *     dsqrt( ((ONE+x)*.5d0)**FIVE3 + ((ONE-x)*.5d0)**FIVE3 )
      dprime(x)=(2.d0**ONETHIRD*FIVE6)*
     *     ( ((ONE+x)*.5d0)*TWO3-((ONE-x)*.5d0)*TWO3 )/
     /     dsqrt( ((ONE+x)*.5d0)**FIVE3+((ONE-x)*.5d0)**FIVE3)

        rsfact=(0.75d0/pi)**ONETHIRD
        if( ipol.eq.1 )then

c       ======> SPIN-RESTRICTED <======

        do 10 n = 1,nq

          rhoval = rho(1,n)
          rho13  = rhoval**ONETHIRD
          rho43  = rhoval*rho13
          rho76  = rhoval**SEV6 
          rs=rsfact/rho13
          rs2=rs*rs
          rs3=rs2*rs
          drsdrho = -ONETHIRD*rs/rhoval
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
            
C         
C         phi(n,gradn)
C         
          expfac=0.d0
          phi=(pfff*ftilde)*(CINF/CC)*dsqgamma/rho76
          if(phi.lt.EXPTOL) expfac=exp(-phi)
C         
C         functional
C         
          func=expfac*CC*gammaval/rho43*fac

          dphidrho =  -phi * ( SEV6/rhoval  +   (CCPRIME/CC) ) 

          AAadd = (  (CCPRIME/CC) -
     -      (dphidrho +   FOURTHIRDS/rhoval) )*func

          Amat(1,n) = Amat(1,n) + AAadd*fac

          if(dsqgamma.gt.0.d0) then

            g =func/gammaval*(2.d0-phi)*fac
            
            Cmat(1,1,n) = Cmat(1,1,n) + g*delrho(1,1,n)
            Cmat(1,2,n) = Cmat(1,2,n) + g*delrho(1,2,n)
            Cmat(1,3,n) = Cmat(1,3,n) + g*delrho(1,3,n)
          endif
          Ec      = Ec + func*qwght(n)*fac
          
   10   continue

      else

c       ======> SPIN-UNRESTRICTED <======

        do 20 n = 1,nq
          rhoval = rho(1,n)
          arho=1.d0/rhoval
          rho13  = rhoval**ONETHIRD
          rho43  = rhoval*rho13
          rho76  = rhoval**SEV6 
          rs=rsfact/rho13
          rs2=rs*rs
          rs3=rs2*rs
          drsdrho = -ONETHIRD*rs*arho
          gammaval=delrho(1,1,n)*delrho(1,1,n)+
     +         delrho(1,2,n)*delrho(1,2,n)+
     +         delrho(1,3,n)*delrho(1,3,n)+
     2         delrho(2,1,n)*delrho(2,1,n)+
     +         delrho(2,2,n)*delrho(2,2,n)+
     +         delrho(2,3,n)*delrho(2,3,n)+(
     x         delrho(1,1,n)*delrho(2,1,n)+
     +         delrho(1,2,n)*delrho(2,2,n)+
     +         delrho(1,3,n)*delrho(2,3,n))*2.d0
          dsqgamma=sqrt(gammaval)
          zeta=(rho(2,n)-rho(3,n))*arho
          add=1.d0/dzeta(zeta)
          adp=dprime(zeta)
          dzetada=(1.d0-zeta)*arho*adp
          dzetadb=-(1.d0+zeta)*arho*adp
C         
C         C(n)
C         
          anum=fff+alpha*rs+beta*rs2
          aden=1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
          CC=zzz+(anum*aden)
          CCPRIME=aden*((alpha+beta*(rs+rs)) -
     -         (anum*aden)*(pgamma+(rs+rs)*delta+beta10*rs2*THREE))
          CCPRIME=CCPRIME*drsdrho
            
C         
C         phi(n,gradn)
C         
          expfac=0.d0
          phi=(pfff*ftilde)*(CINF/CC)*dsqgamma/rho76
          if(phi.lt.EXPTOL) expfac=exp(-phi)
C         
C         functional
C         
           func=expfac*CC*gammaval/rho43*add

          dphidrho =  -phi * ( SEV6/rhoval  +   (CCPRIME/CC) ) 

          AAadd = (  (CCPRIME/CC) -
     -      (dphidrho +   FOURTHIRDS/rhoval) )*func
          Amat(1,n) = Amat(1,n) + (AAadd - (func*add)*dzetada)*fac
          Amat(2,n) = Amat(2,n) + (AAadd - (func*add)*dzetadb)*fac

          if(dsqgamma.gt.TOLL) then

            gaa = func/gammaval*(2.d0-phi)*fac
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

          Ec      = Ec + func*qwght(n)*fac
          
   20   continue

      end if

      end

      Subroutine xc_p81(fac,rho,Amat,nq,ipol,Ec,qwght)
c
c     Ceperley Alder LDA from
c     Perdew Zunger PRB 23, 5048 (1981)
c
      implicit none
      integer nq,ipol
      double precision Ec,fac
c--   > Charge Density 

      double precision rho(ipol*(ipol+1)/2,nq)

c--   > Quadrature Weights

      double precision qwght(nq)

c--   > Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq)
c
      double precision AU,BU,CU,DU,
     *     GU,B1U,B2U,AP,BP,CP,DP,
     *     GP,B1P,B2P,pi,TOLL,ONE3,SEV6,
     *     FOUR3,TWO3
      parameter(
     *     AU=0.0311d0,BU=-0.048d0,CU=0.0020d0,DU=-0.0116d0,
     *     GU=-.1423d0,B1U=1.0529d0,B2U=0.3334d0,
     *     AP=0.01555d0,BP=-0.0269d0,CP=0.0007d0,DP=-0.0048d0,
     *     GP=-.0843d0,B1P=1.3981d0,B2P=0.2611d0)
      parameter (pi=3.1415926535897932385d0)
      Parameter (TOLL=1.D-40)
      Parameter (ONE3=1.d0/3.d0,SEV6=7.d0/6.d0)
      Parameter (FOUR3=4.D0/3.D0,TWO3=2.d0/3.d0)

      integer n
      double precision rhoval,rs,alnrs,eU,vU,drs,adU,ePol,
     ,     fz,zeta,eP,vP,adP,AAadd
      double precision x,fzeta,fzzprime,rsfact
      fzeta(x)=( (1.+x)**FOUR3 + (1.-x)**FOUR3 - 2.d0 ) /
     /     (2.d0**FOUR3-2.d0)
      fzzprime(x)=FOUR3*( (1.d0+x)**ONE3 - (1.d0-x)**ONE3 ) /
     /     (2.d0**FOUR3-2.d0)
C

      rsfact=(0.75d0/pi)**ONE3

      if( ipol.eq.1 )then

c       ======> SPIN-RESTRICTED <======
        do n=1,nq
          
          if(rho(1,n).gt.TOLL) then
            
            rhoval=rho(1,n)
            rs=rsfact/rhoval**ONE3

            if (rs.lt.1.d0) then

              alnrs=log(rs)
              
              eU=AU*alnrs+BU+CU*rs*alnrs+DU*rs
              vU=AU*alnrs + (BU-ONE3*AU) + TWO3*CU*rs*alnrs + 
     +             ONE3*(DU+DU-CU)*rs
            else

              drs=sqrt(rs)

              adU=1.d0/(1.d0+B1U*drs+B2U*rs)

              eU=GU*adU

              vU=eU*(1.d0+SEV6*B1U*drs+FOUR3*B2U*rs)*adU

            endif

              Ec=Ec+eU*qwght(n)*rhoval*fac

            Amat(1,n) = Amat(1,n) + vU*fac

          endif

        enddo

      else

c       ======> SPIN-UNRESTRICTED <======
        do n=1,nq
          
          if(rho(1,n).gt.TOLL) then
            rhoval=rho(1,n)
            zeta=(rho(2,n)-rho(3,n))/rhoval
            fz=fzeta(zeta)
            rs=rsfact/rhoval**ONE3

            if (rs.lt.1.d0) then

              alnrs=log(rs)
              
              eU=AU*alnrs+BU+CU*rs*alnrs+DU*rs
              eP=AP*alnrs+BP+CP*rs*alnrs+DP*rs
              
              vU=AU*alnrs + (BU-ONE3*AU) + TWO3*CU*rs*alnrs + 
     +             ONE3*(DU+DU-CU)*rs
              vP=AP*alnrs + (BP-ONE3*AP) + TWO3*CP*rs*alnrs + 
     +             ONE3*(DP+DP-CP)*rs
            else

              drs=sqrt(rs)

              adU=1.d0/(1.d0+B1U*drs+B2U*rs)
              adP=1.d0/(1.d0+B1P*drs+B2P*rs)

              eU=GU*adU
              eP=GP*adP

              vU=eU*(1.d0+SEV6*B1U*drs+FOUR3*B2U*rs)*adU
              vP=eP*(1.d0+SEV6*B1P*drs+FOUR3*B2P*rs)*adP

            endif

            AAadd = vU + fz*(vP-vU) 

            ePOL=(eP-eU)*fzzprime(zeta)

            Ec=Ec+(eU+fz*(eP-eU))*qwght(n)*rhoval*fac
            
            Amat(1,n) = Amat(1,n) + (AAadd + (1.d0-zeta)*ePOL)*fac
            Amat(2,n) = Amat(2,n) + (AAadd - (1.d0+zeta)*ePOL)*fac

          endif

        enddo

      endif

      return
      end

