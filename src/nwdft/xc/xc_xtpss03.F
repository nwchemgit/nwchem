c    Tao,Perdew, Staroverov, Scuseria exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density
c     References:
c     [a] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         PRL 91, 146401 (2003).


      Subroutine xc_xtpss03(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau,Mmat)


c   
c$Id: xc_xtpss03.F,v 1.2 2004-10-05 18:23:46 pollack Exp $
c
      implicit none
c
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
 

      double precision tol_rho, pi

c
      integer n
      double precision rrho, rho43, rho13, gamma
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,3,ipol), Mmat(nq,*)

      double precision tauN, tauW, tauU

       double precision  p, qtil, x,  al, mt, z
       double precision   rho53, rho83
       double precision   F83, F23, F53, F43, F13
       double precision   G920 
       double precision  b,c,e,es
       double precision    C1, C2, C3
       double precision    kap, mu
       double precision xa, xb,xc,xd
       double precision x1,x2,x3,x4,x5,x6,x7

c       double precision    CX1, CX2, CX3, CX4
c       double precision   bigD, uk, ruk


       double precision   P32, Ax

c      functional derivatives below FFFFFFFFFFFF

       double precision dzdn, dpdn, daldn, dqtdn
       double precision til1, til2
       double precision phaA, phaB
       double precision dtil3dn, dtil2dn, dtil1dn
       double precision ax1, bx1, dax1dn, dbx1dn, dx1dn
       double precision dx2dn
       double precision axd, bxd, dxddn, dxcdn, dx3dn
       double precision dx4dn, dx5dn, dx6dn, dx7dn
       double precision  dxdn
       double precision xmany, dxmanydn
       double precision dmtdn, derivn




       double precision dzdg, dpdg, daldg, dqtdg
       double precision dtil2dg, dax1dg, dbx1dg
       double precision dx1dg, dx2dg
       double precision dag, dbg, dxcdg, dxddg,dx3dg
       double precision dx4dg, dx5dg, dx6dg, dx7dg
       double precision dxmanydg, dxdg, dmtdg, derivg

       double precision dzdt, dpdt, daldt, dqtdt
       double precision dtildt, dax1dt, dax3dt
       double precision dx1dt, dx2dt, dx3dt, dx4dt
       double precision dx5dt, dx6dt, dx7dt
       double precision dxmanydt, dxdt, dmtdt, derivt

c      functional derivatives above FFFFFFFFFFFF


       parameter( pi = 3.1415926535897932384626433832795d0 )

        parameter(kap=0.8040d0, mu=0.21951d0)
        parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
        parameter (F83=8.d0/3.d0, F23=2.d0/3.d0, F53=5.d0/3.d0)
        parameter (G920  =9.d0/20.d0 )

       parameter(b=0.40d0, c=1.59096d0, e=1.537d0)
        parameter (C1  =  10.d0/81.d0, 
     &             C2  = 146.d0/2025.d0,
     &             C3  = -73.d0/405.d0 )

 

         
c
      Ax = (-0.75d0)*(3d0/pi)**F13
      P32 = (3.d0*pi**2)**F23
      es=dsqrt(e)


c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c

         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10

c  rho43= n*e_x^unif=exchange energy per electron for uniform electron gas
c       = n* Ax*n^(1/3)   or n*C*n^(1/3) 

            rho43 = Ax*rho(n,1)**F43  ! Ax*n^4/3
            rrho = 1d0/rho(n,1)       ! reciprocal of rho
            rho13 = F43*rho43*rrho   !functional deriv of rho43

            rho53 = rho(n,1)**F53
            rho83 = rho(n,1)**F83


C  Below we just sum up the LDA contribution to the functional
c            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
               Amat(n,1) = Amat(n,1) + rho13*fac
c            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
c            gam12 = dsqrt(gamma)


            tauN = tau(n,1,1) + 
     &             tau(n,2,1) + 
     &             tau(n,3,1)

            tauW=0.125d0*gamma*rrho

            tauU=0.3d0*P32*rho53 
c
c     Evaluate the Fx, i.e. mt(x) = Fx - 1
c           
            p=gamma/(rho83*P32*4.d0)
            al=(tauN-tauW)/tauU
            al=dabs(al)
            qtil=(G920*(al-1.d0)/((1.d0+b*al*(al-1.d0))**.5d0))
     &      + F23*p
            z=tauW/tauN
            
c          xa= (C1)
           xb=(c*z**2)/( (1+z**2)**2 )
           x1=(C1 + xb)*p
           x2=C2*qtil*qtil
           xc=C3*qtil
           xd=(0.5d0*(.6d0*z)**2  + .5d0*p*p)**.5d0
           x3=xc*xd
           x4=C1*C1*p*p/kap
           x5=2.d0*es*C1*(.6d0*z)**2
           x6= e*mu*p*p*p
c           x7 = (1.d0+es*p)**2
            x7 = (1.d0+es*p)**(-2.d0)

c           x=(x1+x2+x3 +x4 +x5+x6)/x7
            x=(x1+x2+x3 +x4 +x5+x6)*x7

c            if (x.lt.tol_rho) goto 10

c      functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF

C Derivatives wrt n, density   below
c         dzdn=-tauW*rrho/tauN
         dzdn=-z*rrho
         dpdn = -p*rrho*F83

c             phaA= F53*dpdn
c             phaB= -1.d0*dzdn/z**2
c         daldn=F53*p*phaB+phaA*(-1.d0+1.d0/z)
      daldn=F53*( -p*dzdn/z**2 +dpdn*(-1.d0+1.d0/z) )

           til1=al-1.d0
           til2=(1.d0+b*al*(al-1.d0))**-0.5d0
c          til3 not explicitly needed
c            dtil3dn=F23*dpdn
            dtil1dn=daldn
            dtil2dn=(2.d0*b*al*daldn-b*daldn)*
     &            (-.5d0)*til2**3.d0

c         dqtdn =G920*(til2*dtil1dn+til1*dtil2dn)+dtil3dn
         dqtdn =G920*(til2*dtil1dn+til1*dtil2dn)+F23*dpdn
             
               ax1=c*p*z*z
               bx1=(1+z*z)**(-2.d0)
              dax1dn=2.d0*c*p*z*dzdn+c*z*z*dpdn
              dbx1dn=-4.d0*z*dzdn/(1.d0+z*z)**3
           dx1dn=C1*dpdn+ax1*dbx1dn+bx1*dax1dn
           dx2dn=2.d0*C2*qtil*dqtdn

             axd=0.18d0*2.d0*dzdn*z
             bxd=p*dpdn
             dxddn=.5d0*((0.18d0*z*z+.5d0*p*p)**(-.5d0))*(axd+bxd)           
             dxcdn=C3*dqtdn
            dx3dn=xc*dxddn+xd*dxcdn

            dx4dn=(2.d0*p*dpdn*C1**2)/kap
            dx5dn=4.d0*es*C1*0.36d0*z*dzdn
            dx6dn=3.d0*e*mu*p*p*dpdn
            dx7dn=-2.d0*es*dpdn/(1.d0+es*p)**3

            xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydn= dx1dn+dx2dn+dx3dn+dx4dn+dx5dn+dx6dn
            dxdn=x7*dxmanydn+xmany*dx7dn
C Derivatives wrt n, density   above

C Derivatives wrt gamma,    below

            dpdg=1.d0/(4.d0*P32*rho83)
            dzdg=rrho*0.125d0/tauN
            daldg=F53*( -p*dzdg/(z*z)+dpdg*(-1.d0+1.d0/z) )
c   recall     til1=al-1.d0
c     "        til2=(1.d0+b*al*(al-1.d0))**-0.5d0
c     "        til3 not explicitly needed
c     "        dtil1dg=daldg
               dtil2dg=0.5d0*(2.d0*b*al*daldg-b*daldg)*til2**3.d0               
            dqtdg=G920*(til1*dtil2dg+til2*daldg)+F23*dpdg
          

c   recall    ax1=c*p*z*z
c             bx1=(1.d0+z*z)**(-2.d0)
              dax1dg=c*z*z*dpdg+2.d0*c*p*z*dzdg
              dbx1dg=-4.d0*z*dzdg*(1.d0+z*z)**-3.d0
            dx1dg=C1*dpdg+ax1*dbx1dg+bx1*dax1dg
            dx2dg=C2*2.d0*qtil*dqtdg

c   recall    xc=C3*qtil
c             xd=(0.5d0*(.6d0*z)**2  + .5d0*p*p)**.5d0
              dxcdg=C3*dqtdg
               dag= 0.36d0*z*dzdg
               dbg=p*dpdg
              dxddg=.5d0*(dag+dbg)*(0.18d0*z*z+.5d0*p*p)**(-.5d0)
             dx3dg=xc*dxddg+xd*dxcdg
             dx4dg=2.d0*C1*C1*p*dpdg/kap
             dx5dg=1.44d0*es*C1*z*dzdg
             dx6dg=3.d0*e*mu*p*p*dpdg
             dx7dg=-2.d0*es*dpdg*(1.d0+p*es)**-3.d0

c    recall    xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydg= dx1dg+dx2dg+dx3dg+dx4dg+dx5dg+dx6dg
            dxdg=x7*dxmanydg+xmany*dx7dg


            
C Derivatives wrt tau,    below
c ttttttttttttttttttttttttttttttttttttttttttttttttt
          dpdt=0.d0
          dzdt= -tauW/(tauN*tauN)
          daldt=1.d0/tauU
            dtildt=-.5d0*( (2.d0*al-1.d0)*b*daldt )*
     &             (1.d0+b*al*al-b*al)**(-1.5d0)
          dqtdt=G920*( (al-1.d0)*dtildt+ 
     &               daldt*(1.d0+b*al*al-b*al)**(-.5d0) )


            dax1dt=-4.d0*z*dzdt*(1.d0+z**2)**-3
        dx1dt=c*p*dax1dt*z*z+2*c*p*z*dzdt/(1.d0+z*z)**2
        dx2dt=2*c2*qtil*dqtdt
            dax3dt=0.18d0*z*dzdt/(0.18d0*z*z+p*p*.5d0)**(.5d0)
        dx3dt=C3*qtil*dax3dt +
     &   C3*dqtdt*(.18d0*z*z+.5d0*p*p)**(.5d0)
        dx4dt=0.d0
        dx5dt=1.44d0*es*C1*z*dzdt
        dx6dt=0.d0
        dx7dt=0.d0
        
c    recall    xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydt= dx1dt+dx2dt+dx3dt+dx5dt
            dxdt=x7*dxmanydt
c ttttttttttttttttttttttttttttttttttttttttttttttttttt

            mt = kap - kap/(1.d0 + x/kap)
            
            Ex = Ex + mt*rho43*qwght(n)*fac

            dmtdn=dxdn/(1.d0+x/kap)**2
            derivn=mt*rho13+rho43*dmtdn

            dmtdg=dxdg/(1.d0+x/kap)**2
            derivg = rho43*dmtdg

            dmtdt=dxdt/(1.d0+x/kap)**2
            derivt = rho43*dmtdt

            Amat(n,1) = Amat(n,1) + 1.0d0*derivn*fac
            Cmat(n,1)=  Cmat(n,1) + 1.0d0*derivg*fac
            Mmat(n,1)=  Mmat(n,1) + 0.5d0*derivt*fac

10      continue







c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]


         do 20 n = 1, nq
           if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = Ax*(2d0*rho(n,2))**F43 ! spin scaled
            rrho = 0.5d0/rho(n,2)          ! spin scaled
            rho13 = F43*rho43*rrho  !spin scaled & (1/2)factor
            rho53 = (2.d0*rho(n,2))**F53   ! spin scaled
            rho83 = (2.d0*rho(n,2))**F83   ! spin scaled
c  note all the "rho"  quantities ARE spin scaled, for later use

c   in metaGGAs there is no distinct local and nonlocal piece,
c   so that is why we comment out the if statement, and
c   always calculate the local part, the LDA
c            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
                Amat(n,1) = Amat(n,1) + rho13*fac
c            endif

            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
c NOTE  gamma above  is not spin scaled.  that is why
c                           -there is 4.d0*gamma in p
c                           -there is 2.0 in the gam12 term

c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 25


c   below note factor of two for spin scaling
            tauN = 2.d0*( tau(n,1,1) + 
     &                    tau(n,2,1) + 
     &                    tau(n,3,1)   )
            
            tauW=0.125*(4.d0*gamma)*rrho
c            tauW=(4.d0*gamma)/(8.d0*2.d0*rho(n,2))
c                      this 2.0 above for spin scaling

            tauU=0.3d0*P32*rho53 ! rho53 spin scaled already
c
c     Evaluate the GC part of Fx, i.e. mt(x) = Fx - 1
c           
            p=(4.d0*gamma)/(rho83*P32*4.d0)
            al=(tauN-tauW)/tauU
            al=dabs(al)
            qtil=(G920*(al-1.d0)/((1.d0+b*al*(al-1.d0))**.5d0))
     &      + F23*p
            z=tauW/tauN
            
c           xa= (C1)
           xb=(c*z**2)/( (1+z**2)**2 )
           x1=(C1+xb)*p
           x2=C2*qtil*qtil
           xc=C3*qtil
           xd=dsqrt(0.5d0*(.6d0*z)**2  + .5d0*p*p)
           x3=xc*xd
           x4=C1*C1*p**2/kap
           x5=2.d0*es*C1*(.6d0*z)**2
           x6= e*mu*p*p*p
           x7 = (1.d0+es*p)**-2

           x=(x1+x2+x3 +x4 +x5+x6)*x7

c            if (x.lt.tol_rho) goto 25

c      functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF

C Derivatives wrt n, density   below
         dzdn=-tauW*rrho*2.d0/tauN
         dpdn = -p*rrho*2.d0*F83

             phaA= F53*dpdn
             phaB= -1.d0*dzdn/z**2
         daldn=F53*p*phaB+phaA*(-1.d0+1.d0/z)

           til1=al-1.d0
           til2=(1.d0+b*al*(al-1.d0))**-0.5d0
c          til3 not explicitly needed
            dtil3dn=F23*dpdn
            dtil1dn=daldn
            dtil2dn=(2.d0*b*al*daldn-b*daldn)*
     &            (-.5d0)*til2**3

         dqtdn =G920*(til2*dtil1dn+til1*dtil2dn)+dtil3dn
             
               ax1=c*p*z*z
               bx1=(1+z*z)**(-2.d0)
              dax1dn=2.d0*c*p*z*dzdn+c*z*z*dpdn
              dbx1dn=-4.d0*z*dzdn/(1.d0+z**2)**3
           dx1dn=C1*dpdn+ax1*dbx1dn+bx1*dax1dn
           dx2dn=2.d0*C2*qtil*dqtdn

             axd=0.18d0*2.d0*dzdn*z
             bxd=p*dpdn
             dxddn=.5d0*((0.18d0*z*z+.5d0*p*p)**(-.5d0))*(axd+bxd)           
             dxcdn=C3*dqtdn
            dx3dn=xc*dxddn+xd*dxcdn

            dx4dn=(2.d0*p*dpdn*C1**2)/kap
            dx5dn=4.d0*es*C1*0.36d0*z*dzdn
            dx6dn=3.d0*e*mu*p*p*dpdn
            dx7dn=-2.d0*es*dpdn/(1.d0+es*p)**3

            xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydn= dx1dn+dx2dn+dx3dn+dx4dn+dx5dn+dx6dn
            dxdn=x7*dxmanydn+xmany*dx7dn
C Derivatives wrt n, density   above

C Derivatives wrt gamma,    below

            dpdg=4.d0/(4.d0*P32*rho83)
            dzdg=4.d0*rrho*0.125d0/tauN
            daldg=F53*( -p*dzdg/(z*z)+dpdg*(-1.d0+1.d0/z) )
c   recall     til1=al-1.d0
c     "        til2=(1.d0+b*al*(al-1.d0))**-0.5d0
c     "        til3 not explicitly needed
c     "        dtil1dg=daldg
               dtil2dg=0.5d0*(2.d0*b*al*daldg-b*daldg)*til2**3.d0               
            dqtdg=G920*(til1*dtil2dg+til2*daldg)+F23*dpdg
          

c   recall    ax1=c*p*z*z
c             bx1=(1.d0+z*z)**(-2.d0)
              dax1dg=c*z*z*dpdg+2.d0*c*p*z*dzdg
              dbx1dg=-4.d0*z*dzdg*(1.d0+z*z)**-3.d0
            dx1dg=C1*dpdg+ax1*dbx1dg+bx1*dax1dg
            dx2dg=C2*2.d0*qtil*dqtdg

c   recall    xc=C3*qtil
c             xd=(0.5d0*(.6d0*z)**2  + .5d0*p*p)**.5d0
              dxcdg=C3*dqtdg
               dag= 0.36d0*z*dzdg
               dbg=p*dpdg
              dxddg=.5d0*(dag+dbg)*(0.18d0*z*z+.5d0*p*p)**(-.5d0)
             dx3dg=xc*dxddg+xd*dxcdg
             dx4dg=2.d0*C1*C1*p*dpdg/kap
             dx5dg=1.44d0*es*C1*z*dzdg
             dx6dg=3.d0*e*mu*p*p*dpdg
             dx7dg=-2.d0*es*dpdg*(1.d0+p*es)**-3.d0

c    recall    xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydg= dx1dg+dx2dg+dx3dg+dx4dg+dx5dg+dx6dg
            dxdg=x7*dxmanydg+xmany*dx7dg


C Derivatives wrt tau,    below
c ttttttttttttttttttttttttttttttttttttttttttttttttt
          dpdt=0.d0
          dzdt= -tauW*2.d0/(tauN*tauN)
          daldt=2.d0/tauU
            dtildt=-.5d0*( (2.d0*al-1.d0)*b*daldt )*
     &             (1.d0+b*al*al-b*al)**(-1.5d0)
          dqtdt=G920*( (al-1.d0)*dtildt+ 
     &               daldt*(1.d0+b*al*al-b*al)**(-.5d0) )


            dax1dt=-4.d0*z*dzdt*(1.d0+z**2)**-3
        dx1dt=c*p*dax1dt*z*z+2*c*p*z*dzdt/(1.d0+z*z)**2
        dx2dt=2*c2*qtil*dqtdt
            dax3dt=0.18d0*z*dzdt/(0.18d0*z*z+p*p*.5d0)**(.5d0)
        dx3dt=C3*qtil*dax3dt +
     &   C3*dqtdt*(.18d0*z*z+.5d0*p*p)**(.5d0)
        dx4dt=0.d0
        dx5dt=1.44d0*es*C1*z*dzdt
        dx6dt=0.d0
        dx7dt=0.d0
        
c    recall    xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydt= dx1dt+dx2dt+dx3dt+dx5dt
            dxdt=x7*dxmanydt
c ttttttttttttttttttttttttttttttttttttttttttttttttttt




            mt = kap - kap/(1.d0 + x/kap)
            
            Ex = Ex + mt*rho43*qwght(n)*0.5d0*fac

            dmtdn=dxdn/(1+x/kap)**2
c            derivn=mt*rho13+ 0.5d0*rho43*dmtdn

            dmtdg=dxdg/(1.d0 + x/kap)**2
c            derivg=rho43*dmtdg

            dmtdt=dxdt/(1.d0+x/kap)**2
            derivt = rho43*dmtdt


c   rho13 has (1/2) built in.  see .5 in second term for sp.sc
            Amat(n,1) = Amat(n,1)+fac*(mt*rho13+ 0.5d0*rho43*dmtdn)
            Cmat(n,1) = Cmat(n,1) + 0.5d0*rho43*dmtdg*fac

            Mmat(n,1) = Mmat(n,1) + 0.5d0*.5d0*derivt*fac




          
c
c     Beta               BETA           BETA
c

25         continue

c
c     Beta
c
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = Ax*(2d0*rho(n,3))**F43 ! spin scaled
            rrho = 0.5d0/rho(n,3)          ! spin scaled
            rho13 = F43*rho43*rrho  !spin scaled & (1/2)factor
            rho53 = (2.d0*rho(n,3))**F53   ! spin scaled
            rho83 = (2.d0*rho(n,3))**F83   ! spin scaled
c  note all the "rho"  quantities ARE spin scaled, for later use

c            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
                Amat(n,2) = Amat(n,2) + rho13*fac
c            endif

            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
c NOTE  gamma above  is not spin scaled.  that is why
c                           -there is 4.d0*gamma in p
c                           -there is 2.0 in the gam12 term

c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 20


c   below note factor of two for spin scaling
            tauN = 2.d0*( tau(n,1,2) + 
     &                    tau(n,2,2) + 
     &                    tau(n,3,2)   )
            

            tauW=(4.d0*gamma)*rrho/8.d0
c            tauW=(4.d0*gamma)/(8.d0*2.d0*rho(n,3))
c                      this 2.0 above for spin scaling
c             but rrho is alread spin scaled

            tauU=0.3d0*P32*rho53 ! rho53 spin scaled already
c
c           
c           4gamma for spin scaling below
            p=(4.d0*gamma)/(rho83*P32*4.d0)
            al=(tauN-tauW)/tauU
            al=dabs(al)
            qtil=(G920*(al-1.d0)/((1.d0+b*al*(al-1.d0))**.5d0))
     &      + F23*p
            z=tauW/tauN
            
c           xa= (C1)
           xb=(c*z**2)/( (1+z**2)**2 )
           x1=(C1+xb)*p
           x2=C2*qtil*qtil
           xc=C3*qtil
c          note below we left off the p**2 term
           xd=dsqrt(0.5d0*(.6d0*z)**2  + .5d0*p*p)
           x3=xc*xd
           x4=C1*C1*p**2/kap
           x5=2.d0*es*C1*(.6d0*z)**2
           x6= e*mu*p*p*p
           x7 = (1.d0+es*p)**-2

           x=(x1+x2+x3 +x4 +x5+x6)*x7


c            if (x.lt.tol_rho) goto 20



c      functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF

C Derivatives wrt n, density   below
         dzdn=-tauW*rrho*2.d0/tauN
         dpdn = -p*rrho*2.d0*F83

             phaA= F53*dpdn
             phaB= -1.d0*dzdn/z**2
         daldn=F53*p*phaB+phaA*(-1.d0+1.d0/z)

           til1=al-1.d0
           til2=(1.d0+b*al*(al-1.d0))**-0.5d0
c          til3 not explicitly needed
            dtil3dn=F23*dpdn
            dtil1dn=daldn
            dtil2dn=(2.d0*b*al*daldn-b*daldn)*
     &            (-.5d0)*til2**3

         dqtdn =G920*(til2*dtil1dn+til1*dtil2dn)+dtil3dn
             
               ax1=c*p*z*z
               bx1=(1+z*z)**(-2.d0)
              dax1dn=2.d0*c*p*z*dzdn+c*z*z*dpdn
              dbx1dn=-4.d0*z*dzdn/(1.d0+z**2)**3
           dx1dn=C1*dpdn+ax1*dbx1dn+bx1*dax1dn
           dx2dn=2.d0*C2*qtil*dqtdn

             axd=0.18d0*2.d0*dzdn*z
             bxd=p*dpdn
             dxddn=.5d0*((0.18d0*z*z+.5d0*p*p)**(-.5d0))*(axd+bxd)           
             dxcdn=C3*dqtdn
            dx3dn=xc*dxddn+xd*dxcdn

            dx4dn=(2.d0*p*dpdn*C1**2)/kap
            dx5dn=4.d0*es*C1*0.36d0*z*dzdn
            dx6dn=3.d0*e*mu*p*p*dpdn
            dx7dn=-2.d0*es*dpdn/(1.d0+es*p)**3

            xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydn= dx1dn+dx2dn+dx3dn+dx4dn+dx5dn+dx6dn
            dxdn=x7*dxmanydn+xmany*dx7dn
C Derivatives wrt n, density   above


C Derivatives wrt gamma,    below

            dpdg=4.d0/(4.d0*P32*rho83)
            dzdg=4.d0*rrho*0.125d0/tauN
            daldg=F53*( -p*dzdg/(z*z)+dpdg*(-1.d0+1.d0/z) )
c   recall     til1=al-1.d0
c     "        til2=(1.d0+b*al*(al-1.d0))**-0.5d0
c     "        til3 not explicitly needed
c     "        dtil1dg=daldg
               dtil2dg=0.5d0*(2.d0*b*al*daldg-b*daldg)*til2**3.d0               
            dqtdg=G920*(til1*dtil2dg+til2*daldg)+F23*dpdg
          

c   recall    ax1=c*p*z*z
c             bx1=(1.d0+z*z)**(-2.d0)
              dax1dg=c*z*z*dpdg+2.d0*c*p*z*dzdg
              dbx1dg=-4.d0*z*dzdg*(1.d0+z*z)**-3.d0
            dx1dg=C1*dpdg+ax1*dbx1dg+bx1*dax1dg
            dx2dg=C2*2.d0*qtil*dqtdg

c   recall    xc=C3*qtil
c             xd=(0.5d0*(.6d0*z)**2  + .5d0*p*p)**.5d0
              dxcdg=C3*dqtdg
               dag= 0.36d0*z*dzdg
               dbg=p*dpdg
              dxddg=.5d0*(dag+dbg)*(0.18d0*z*z+.5d0*p*p)**(-.5d0)
             dx3dg=xc*dxddg+xd*dxcdg
             dx4dg=2.d0*C1*C1*p*dpdg/kap
             dx5dg=1.44d0*es*C1*z*dzdg
             dx6dg=3.d0*e*mu*p*p*dpdg
             dx7dg=-2.d0*es*dpdg*(1.d0+p*es)**-3.d0

c    recall    xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydg= dx1dg+dx2dg+dx3dg+dx4dg+dx5dg+dx6dg
            dxdg=x7*dxmanydg+xmany*dx7dg

C Derivatives wrt tau,    below
c ttttttttttttttttttttttttttttttttttttttttttttttttt
          dpdt=0.d0
          dzdt= -tauW*2.d0/(tauN*tauN)
          daldt=2.d0/tauU
            dtildt=-.5d0*( (2.d0*al-1.d0)*b*daldt )*
     &             (1.d0+b*al*al-b*al)**(-1.5d0)
          dqtdt=G920*( (al-1.d0)*dtildt+ 
     &               daldt*(1.d0+b*al*al-b*al)**(-.5d0) )


            dax1dt=-4.d0*z*dzdt*(1.d0+z**2)**-3
        dx1dt=c*p*dax1dt*z*z+2*c*p*z*dzdt/(1.d0+z*z)**2
        dx2dt=2*c2*qtil*dqtdt
            dax3dt=0.18d0*z*dzdt/(0.18d0*z*z+p*p*.5d0)**(.5d0)
        dx3dt=C3*qtil*dax3dt +
     &   C3*dqtdt*(.18d0*z*z+.5d0*p*p)**(.5d0)
        dx4dt=0.d0
        dx5dt=1.44d0*es*C1*z*dzdt
        dx6dt=0.d0
        dx7dt=0.d0
        
c    recall    xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydt= dx1dt+dx2dt+dx3dt+dx5dt
            dxdt=x7*dxmanydt
c ttttttttttttttttttttttttttttttttttttttttttttttttttt



            mt = kap - kap/(1.d0 + x/kap)
            
            Ex = Ex + mt*rho43*qwght(n)*0.5d0*fac

            dmtdn=dxdn/(1+x/kap)**2
            derivn=mt*rho13+ 0.5d0*rho43*dmtdn
c above  rho13 has (1/2) built in.  see .5 in second term for sp.sc

            dmtdg=dxdg/(1.d0 + x/kap)**2
c            derivg=rho43*dmtdg

            dmtdt=dxdt/(1.d0+x/kap)**2
            derivt = rho43*dmtdt


            Amat(n,2) = Amat(n,2) + derivn*fac
           Cmat(n,3) = Cmat(n,3) + 0.5d0*rho43*dmtdg*fac

            Mmat(n,2) = Mmat(n,2) + 0.5d0*.5d0*derivt*fac




c
20      continue
      endif
c
      return
      end




      Subroutine xc_xtpss03_d2()
      call errquit(' not coded ',0,0)
      return
      end


