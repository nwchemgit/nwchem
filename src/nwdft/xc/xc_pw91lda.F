      subroutine xc_pw91lda(fac,rho,Amat,nq,ipol,Ec,qwght)
c  uniform-gas correlation of Perdew and Wang 1991
c  input: seitz radius (rs), relative spin polarization (zet)
c  output: correlation energy per electron (ec), up- and down-spin
c     potentials (vcup,vcdn), derivatives of ec wrt rs (ecrs) & zet (eczet)
c  output: correlation contribution (alfc) to the spin stiffness
      implicit none
      double precision fac ! [input]
      integer nq,ipol      ! [input]
      double precision Ec  ! [input/output]

      double precision rho(ipol*(ipol+1)/2,nq) ! Charge Density 
      double precision qwght(nq)!Quadrature Weights
      double precision Amat(ipol,nq)!Sampling Matrices for XC V & E
      
      double precision eu,ep,eurs,eprs,alfm,alfrsm,alfc,epw91
      double precision zet
      double precision gam,fzz,thrd,thrd4,f
      double precision rhoval,arho,pi,rsfact
      integer n
      double precision vcup,vcdn,comm,eczet,fz,ecrs,z4,rs,
     ,     rho13

      data gam,fzz/0.5198421d0,1.709921d0/
      data thrd,thrd4/0.333333333333d0,1.333333333333d0/
      parameter (pi=3.1415926535897932385d0)
C
C
C
      rsfact=(0.75d0/pi)**thrd
      do n = 1,nq
        rhoval = rho(1,n)
        rho13  = rhoval**thrd
        arho=1.d0/rhoval
        rs=rsfact/rho13
        if (ipol.eq.1) then
          zet=0.d0
          z4 =0.d0
        else
          zet=(rho(2,n)-rho(3,n))*arho
          z4 = zet**4
        endif

        call xc_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,
     *       1.6382d0,0.49294d0,1.00d0,rs,eu,eurs)
        call xc_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,
     *       3.3662d0,0.62517d0,1.00d0,rs,ep,eprs)
        call xc_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0,
     *       0.88026d0,0.49671d0,1.00d0,rs,alfm,alfrsm)
c       alfm is minus the spin stiffness alfc
        alfc = -alfm
        if(ipol.eq.2) then
          f = ((1.d0+zet)**thrd4+(1.d0-zet)**thrd4-2.d0)/gam
          epw91=(eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz)
        else
          epw91=eu
        endif
        ec = ec+epw91*qwght(n)*rhoval*fac
c       energy done. now the potential:
        ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
        comm = epw91 -rs*ecrs/3.d0-zet*eczet
        if(ipol.eq.2) then
          fz = thrd4*((1.d0+zet)**thrd-(1.d0-zet)**thrd)/gam
          eczet = 4.d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu
     *         -(1.d0-z4)*alfm/fzz)
          vcup = comm + eczet
          vcdn = comm - eczet
          Amat(1,n) = Amat(1,n) + vcup*fac
          Amat(2,n) = Amat(2,n) + vcdn*fac
        else
          Amat(1,n) = Amat(1,n) + comm*fac
        endif

      enddo
      return
      end
      subroutine xc_gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)
c  called by subroutine xc_pw91lsd
      implicit none
      double precision a,a1,b1,b2,b3,b4,p,rs ! [input]
      double precision gg,ggrs ! [output]
c
c     local
c
      double precision p1,q0,rs12,rs32,rsp,q1,q2,q3
      p1 = p + 1.d0
      q0 = -2.d0*a*(1.d0+a1*rs)
      rs12 = dsqrt(rs)
      rs32 = rs12**3
      rsp = rs**p
      q1 = 2.d0*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
      q2 = dlog(1.d0+1.d0/q1)
      gg = q0*q2
      q3 = a*(b1/rs12+2.d0*b2+3.d0*b3*rs12+2.d0*b4*p1*rsp)
      ggrs = -2.d0*a*a1*q2-q0*q3/(q1**2+q1)
      return
      end
