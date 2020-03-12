c   VS98  exchange functional   
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauU - uniform-gas KE density
c                              ijzy - 1  VS98
c                              ijzy - 2  M06-L  
c                              ijzy - 3  M06-HF
c                              ijzy - 4  M06
c                              ijzy - 6  revM06-L
c                              ijzy - 7  revM06
c                              
c     References: 
c 
c     [a] T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998). 
c     [b] Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101 (2006).



      Subroutine xc_xvs98(tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau, Mmat, ijzy)


c   
c$Id$
c
      implicit none
c
#include "errquit.fh"
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
      integer n, ijzy
      double precision rrho, rho43, rho13, rhoo, rho53, rho83
      double precision  Gamma
c
c     kinetic energy density or tau
c
      double precision tau(nq,ipol), Mmat(nq,*)

      double precision tauN,tauu,DTol
      double precision Tiny, f13, f43, f53, f83, f113 
      double precision gx, gg, x, z, kx,xk,zk
      double precision One, Two, Three, Four, Five, Six, Seven, Eight
      double precision Nine, F10, F11
      double precision cf, Axlsda, r1, r2, r3, r4, r5, r6

c      functional derivatives below FFFFFFFFFFFF

       double precision dxdr, dxdg, dzdr, dzdt, dgdx, dgdz

c      functional derivatives above FFFFFFFFFFFF


cedo       parameter( pi = 3.1415926535897932384626433832795d0 )
         
       parameter (cf = 9.115599720d0, Axlsda = -0.9305257363491d0 )
       parameter (gg  = 0.00186726d0) 
       parameter (f13=1.d0/3.d0,f43=4.0d0/3.0d0,f53=5.0d0/3.0d0)
       parameter (f83=8.d0/3.0d0, F113=11.0d0/3.d0)
       parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &             Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &             Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)
      pi=acos(-1d0)


      if (ijzy.eq.1) then
c
c     Parameters for VS98
c
        r1=  -9.800683d-01
        r2=  -3.556788d-03
        r3=   6.250326d-03
        r4=  -2.354518d-05
        r5=  -1.282732d-04
        r6=   3.574822d-04
      elseif (ijzy.eq.2) then
c
c     Parameters for M06-L
c
        r1 =   6.012244D-01*Axlsda
        r2 =   4.748822D-03*Axlsda
        r3 =  -8.635108D-03*Axlsda
        r4 =  -9.308062D-06*Axlsda
        r5 =   4.482811D-05*Axlsda
        r6 =   0.000000D+00
      elseif (ijzy.eq.3) then
c
c     Parameters for M06-HF
c
        r1 =   -1.179732D-01*Axlsda
        r2 =   -2.500000D-03*Axlsda
        r3 =   -1.180065D-02*Axlsda
        r4 =   0.000000D+00
        r5 =   0.000000D+00
        r6 =   0.000000D+00
      elseif (ijzy.eq.4) then
c
c     Parameters for M06
c
        r1 =   1.422057D-01*Axlsda
        r2 =   7.370319D-04*Axlsda
        r3 =   -1.601373D-02*Axlsda
        r4 =   0.000000D+00
        r5 =   0.000000D+00
        r6 =   0.000000D+00
      elseif (ijzy.eq.6) then
c
c     Parameters for revM06-L
c
        r1 =   -4.23227252D-01*Axlsda
        r2 =   0.000000D+00*Axlsda
        r3 =   3.724234D-03*Axlsda
        r4 =   0.000000D+00*Axlsda
        r5 =   0.000000D+00*Axlsda
        r6 =   0.000000D+00
      elseif (ijzy.eq.7) then
c
c     Parameters for revM06
c
        r1 =   -5.523940140D-02*Axlsda
        r2 =   0.000000D+00*Axlsda
        r3 =  -3.782631233D-03*Axlsda
        r4 =   0.000000D+00*Axlsda
        r5 =   0.000000D+00*Axlsda
        r6 =   0.000000D+00
      else
        call errquit("xc_xvs98: illegal value of ijzy",ijzy,UERR)
      endif
      
      DTol = tol_rho
      
         
c     


c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c
         do 10 n = 1, nq
            rhoo = rho(n,1)/Two   ! rho_sigma
            if (rhoo.lt.DTol) goto 10
            rho43 = rhoo**f43  
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**f53
            rho83 = rho53*rhoo 
c
             
            tauN = tau(n,1)
            if(taun.lt.dtol) goto 10
            tauu=tauN
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
         if(sqrt(gamma).lt.dtol) goto 10
            gamma = gamma/Four
            x = gamma/rho83
            dxdr = -f83*x*rrho
            dxdg = One/rho83
            z = tauu/rho53 - cf
            dzdr = -f53 * tauu/rho83
            dzdt = One/rho53
            kx = One + gg*x + gg*z
            xk = x/kx
            zk = z/kx
            call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
           
            Ex = Ex + Two*rho43*gx*qwght(n)
            if(ldew) func(n)=func(n)+ Two*rho43*gx
c
c     functional derivatives 
c
            Amat(n,1) = Amat(n,1) + f43*rho13*gx + 
     &                  rho43*(dgdx*dxdr + dgdz*dzdr) 
            Cmat(n,1)=  Cmat(n,1) + rho43*(dgdx*dxdg) 
            Mmat(n,1)=  Mmat(n,1) + rho43*(dgdz*dzdt)

10      continue


c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======


         do 20 n = 1, nq
           if (rho(n,1).lt.DTol) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,2).lt.DTol) goto 25
             rhoo = rho(n,2)
             rho43 = rhoo**f43
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**f53
             rho83 = rho53*rhoo

c

             tauN = tau(n,1)
             tauu = tauN*Two

            Gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)

            x = gamma/rho83
            dxdr = -f83*x*rrho
            dxdg = One/rho83
            z = tauu/rho53 - cf
            dzdr = -f53 * tauu/rho83
            dzdt = One/rho53
            kx = One + gg*x + gg*z
            xk = x/kx
            zk = z/kx
            call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
            Ex = Ex + rho43*gx*qwght(n)
            if(ldew) func(n)=func(n)+ rho43*gx
c
c     functional derivatives
c
            Amat(n,1) = Amat(n,1) + f43*rho13*gx + 
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)    
            Cmat(n,1)=  Cmat(n,1) + rho43*(dgdx*dxdg)      
            Mmat(n,1)=  Mmat(n,1) + rho43*(dgdz*dzdt)


c
c     Beta               BETA           BETA
c

25         continue

c
c     Beta
c
            if (rho(n,3).lt.DTol) goto 20
             rhoo = rho(n,3)
             rho43 = rhoo**f43
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**f53
             rho83 = rho53*rhoo

c

             tauN = tau(n,2)
             tauu = Two*tauN


            Gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)

            x = gamma/rho83
            dxdr = -f83*x*rrho
            dxdg = One/rho83
            z = tauu/rho53 - cf
            dzdr = -f53 * tauu/rho83
            dzdt = One/rho53
            kx = One + gg*x + gg*z
            xk = x/kx
            zk = z/kx
            call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)

            Ex = Ex + rho43*gx*qwght(n)
            if(ldew) func(n)=func(n)+ rho43*gx
c
c     functional derivatives
c
            Amat(n,2) = Amat(n,2) + f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
            Cmat(n,3)=  Cmat(n,3) + rho43*(dgdx*dxdg)
            Mmat(n,2)=  Mmat(n,2) + rho43*(dgdz*dzdt)

c
20      continue
      endif
c
      return
      end


      Subroutine xc_xvs98_d2()
      implicit none
      call errquit(' xvs98: d2 not coded ',0,0)
      return
      end

      Subroutine gvt4(g,dgdx,dgdz,xk,zk,k,c,ct,r1,r2,r3,r4,r5,r6)
      Implicit none
c
c     Evaluate the GVT4 form in the VS98 functional
c
c
c    some working variables
      double precision g,dgdx,dgdz,xk,zk,k,c,ct,r1,r2,r3,r4,r5,r6
      double precision sxk,szk,sk,sc,sct,sr1,sr2,sr3,sr4,sr5,sr6,sk2
      double precision One, Two, Three, Four, Six
      Data One/1.0d0/, Two/2.0d0/, Three/3.0d0/, Four/4.0d0/, Six/6.0d0/
C
      sxk = xk
      szk = zk
      sk = k
      sc = c
      sct = ct
      sr1 = r1
      sr2 = r2
      sr3 = r3
      sr4 = r4
      sr5 = r5
      sr6 = r6
      sk2 = sk*sk
      g =  (sr1 + sr2*sxk + sr3*szk
     $  +sr4*sxk*sxk + sr5*szk*sxk + sr6*szk*szk)/sk
      dgdx =   (-sr1*sc
     $  +sr2*(One-Two*sc*sxk)
     $  -Two*sr3*szk*sc
     $  +sr4*(Two*sxk-Three*sxk*sxk*sc)
     $  +sr5*(szk -Three*szk*sxk*sc)
     $  -Three*sr6*szk*szk*sc )/sk2
      dgdz =   (-sr1*sct
     $  -Two*sr2*sxk*sct
     $  +sr3*(One-Two*szk*sct)
     $  -Three*sr4*sxk*sxk*sct
     $  +sr5*(sxk-Three*sxk*szk*sct)
     $  +sr6*(Two*szk-Three*szk*szk*sct))/sk2

      return
      end


