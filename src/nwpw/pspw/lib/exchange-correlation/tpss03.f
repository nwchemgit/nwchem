      subroutine gen_TPSS03_restricted(n2ft3d,rho_in,agr_in,tau_in,
     >                               x_parameter,c_parameter,
     >                               xce,fn,fdn,fdtau)
      implicit none
      integer n2ft3d
      real*8 rho_in(*),agr_in(*),tau_in(*)     !*** input  ****
      real*8 x_parameter,c_parameter
      real*8 xce(*),fn(*),fdn(*),fdtau(*)      !*** output ****

********* local declarations******
      integer i
      real*8 n,agr,tau
      real*8 n_13,n_53,n_83,inv_n,agr2,tauW,tauU
      real*8 Cx,P23,es
      real*8 p,z,alpha,qb,p2,p3,z2,qb2,fa,thresA,dalpha_dfa
      real*8 x1a,x1,x2,x3a,x3b,x3,x4,x5,x6,x7,xnum,x
      real*8 dp_dn,dz_dn,dalpha_dn,qba,qbb,dqba_dn,dqbb_dn,dqb_dn
      real*8 dx1_dn,dx2_dn,dx3_dn,dx4_dn,dx5_dn,dx6_dn,dx7_dn 
      real*8 dx3a_dn,dx3b_dn,dxnum_dn,dx_dn
      real*8 dp_dagr,dz_dagr,dalpha_dagr,dqba_dagr,dqbb_dagr,dqb_dagr
      real*8 dx1_dagr,dx2_dagr,dx3_dagr,dx4_dagr,dx5_dagr,dx6_dagr
      real*8 dx7_dagr,dx3a_dagr,dx3b_dagr,dxnum_dagr,dx_dagr
      real*8 dz_dtau,dalpha_dtau,dqb_dtau
      real*8 dx1_dtau,dx2_dtau,dx3_dtau,dx5_dtau,dxnum_dtau,dx_dtau
      real*8 Fx,ex0,nex0,dFx_dx,dFx_dn,dFx_dagr,dFx_dtau
      real*8 ex,fnx,fdnx,fdtaux
      real*8 ec,fnc,fdnc,fdtauc

********* constants ***************
      real*8 pi,thrd,twthrd,frthrd,fvthrd,etthrd,C920,C1,C2,C3
      parameter (pi     = 3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (twthrd = 2.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
      parameter (C920   = 9.0d0/20.0d0)
      parameter (C1     = 10.0d0/81.0d0)
      parameter (C2     = 146.0d0/2025.0d0)
      parameter (C3     = -73.0d0/405.0d0)

********* density cutoff parameters******
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)

********* VS98 constants ******
      real*8 kappa,mu,b,c,e
      parameter (kappa = 0.8040d0)
      parameter (mu = 0.21951d0)
      parameter (b=0.40d0)
      parameter (c=1.59096d0)
      parameter (e=1.537d0)

      Cx = (-0.75d0)*(3.0d0/pi)**thrd
      P23 = (3.0d0*pi**2.0d0)**twthrd
      es=dsqrt(e)

      do i=1,n2ft3d
        n       = rho_in(i) + ETA
        agr     = agr_in(i) + ETA
        tau     = tau_in(i) + ETA

        n_13 = n**thrd
        n_53 = n_13*n_13*n
        n_83 = n_53*n
        inv_n = 1.0d0/n
        agr2 = agr*agr
        tauW = 0.125d0*agr2*inv_n
        tauU = 0.3d0*P23*n_53
        ex0 = Cx*n_13
        p = agr2/(4.0d0*P23*n_83)
        z = tauW/tau

        if (z .gt. 1.0d0) then
           z = 1.0d0
           dz_dn = 0.0d0
           dz_dagr = 0.0d0
           dz_dtau = 0.0d0
        else 
           z = tauW/tau
           dz_dn     = -z*inv_n
           dz_dagr = 2.0d0*z/agr
           dz_dtau = -z/tau
        end if

        p2 = p*p
        p3 = p2*p
        z2 = z*z
        dp_dn     = -etthrd*p*inv_n
        dp_dagr = 2.0d0*p/agr
        !dp_dtau = 0.0d0

        !fa = (tau - tauW)/tauU
        fa = fvthrd*p*(1.0d0/z - 1.0d0)
        thresA = 0.001d0
        alpha     = 0.0d0
        dalpha_dn = 0.0d0
        dalpha_dagr = 0.0d0
        dalpha_dtau = 0.0d0

        if (fa .ge. thresA) then

          alpha     = fa
          dalpha_dn = fvthrd*(-p*dz_dn/z2 + dp_dn*(1.d0/z - 1.0d0))
          dalpha_dagr = (alpha/p)*dp_dagr - fvthrd*(p/z2)*dz_dagr
          dalpha_dtau = 1.0d0/tauU
          !dalpha_dtau = fvthrd*p*(-1.0d0/z2)*dz_dtau

        else if (fa .le. 0.0d0) then

          alpha     = 0.0d0
          dalpha_dn = 0.0d0
          dalpha_dagr = 0.0d0
          dalpha_dtau = 0.0d0

        else if (fa .gt. 0.0d0 .and. fa .lt. thresA) then
          alpha       = 2.0d0*fa*fa/thresA 
     &                  - fa*fa*fa/(thresA*thresA)
          dalpha_dfa  = (4.0d0*fa/thresA 
     &                   - 3.0d0*fa*fa/(thresA*thresA))
          dalpha_dn   = dalpha_dfa*fvthrd*(-p*dz_dn/z2 
     &                + dp_dn*(-1.d0 + 1.d0/z))
          dalpha_dagr   = dalpha_dfa*((alpha/p)*dp_dagr 
     &                - fvthrd*(p/z2)*dz_dagr)
          dalpha_dtau  = dalpha_dfa/tauU
          !dalpha_dtau  = dalpha_dfa*fvthrd*p*(-1.0d0/z2)*dz_dtau

        end if

        qb = C920*(alpha - 1.0d0)/dsqrt(1.0d0+b*alpha*(alpha - 1.0d0))
     &     + twthrd*p

        qb2 = qb*qb

        x1a = (c*z2)/((1.0d0 + z2)**2.0d0)
        x1  = (C1 + x1a)*p
        x2  = C2*qb2
        x3a = C3*qb
        x3b = dsqrt(0.5d0*(0.36d0*z2 + p2))
        x3  = x3a*x3b
        x4  = C1*C1*p2/kappa
        x5  = 2.0d0*es*C1*0.36d0*z2
        x6  = e*mu*p3
        x7  = 1.0d0/((1.0d0 + es*p)**2.0d0)

        xnum = x1 + x2 + x3 + x4 + x5 + x6
        x=xnum*x7

********* dxdn *********

        qba = alpha - 1.0d0
        qbb = 1.0d0/dsqrt(1.0d0 + b*alpha*(alpha - 1.0d0))
        dqba_dn = dalpha_dn
        dqbb_dn = b*dalpha_dn*(2.0d0*alpha - 1.0d0)*
     &            (-0.5d0)*(qbb**3.0d0)
        dqb_dn = C920*(qbb*dqba_dn + qba*dqbb_dn) + twthrd*dp_dn

        dx1_dn = (x1/p)*dp_dn 
     &        + 2.0d0*c*p*z*dz_dn/((1.0d0 + z2)**3.0d0)*(1.0d0 - z2)
        dx2_dn = 2.0d0*C2*qb*dqb_dn

        dx3a_dn = C3*dqb_dn
        dx3b_dn = 0.5d0/x3b*(0.36d0*z*dz_dn + p*dp_dn)
        dx3_dn  = x3b*dx3a_dn+x3a*dx3b_dn

        dx4_dn  = (2.0d0*x4/p)*dp_dn
        dx5_dn  = (2.0d0*x5/z)*dz_dn
        dx6_dn  = (3.0d0*x6/p)*dp_dn
        dx7_dn  = -2.0d0*es*dp_dn/(1.0d0 + es*p)**3.0d0

        dxnum_dn = dx1_dn + dx2_dn + dx3_dn + dx4_dn + dx5_dn + dx6_dn
        dx_dn    = x7*dxnum_dn + xnum*dx7_dn

********* dxdagr *********
        dqba_dagr = dalpha_dagr
        dqbb_dagr = -0.5d0*dalpha_dagr*b*(2.0d0*alpha - 1.0d0)
     &              *qbb**3.0d0
        dqb_dagr = C920*(qba*dqbb_dagr + qbb*dqba_dagr) 
     &             + twthrd*dp_dagr
        dx1_dagr = (x1/p)*dp_dagr 
     &         + 2.0d0*c*p*z*dz_dagr/((1.0d0 + z2)**3.0d0)*(1.0d0 - z2)

        dx2_dagr = C2*2.0d0*qb*dqb_dagr

        dx3a_dagr = C3*dqb_dagr
        dx3b_dagr = 0.5d0/x3b*( 0.36d0*z*dz_dagr + p*dp_dagr)
        dx3_dagr = x3b*dx3a_dagr + x3a*dx3b_dagr

        dx4_dagr = (2.0d0*x4/p)*dp_dagr
        dx5_dagr = (2.0d0*x5/z)*dz_dagr
        dx6_dagr = (3.0d0*x6/p)*dp_dagr

        dx7_dagr = -2.0d0*es*dp_dagr/(1.0d0 + es*p)**3.0d0

        dxnum_dagr = dx1_dagr + dx2_dagr + dx3_dagr 
     &             + dx4_dagr + dx5_dagr + dx6_dagr
        dx_dagr = x7*dxnum_dagr + xnum*dx7_dagr

********* dxdtau *********

        dqb_dtau = C920*dalpha_dtau*qbb*(1.0d0
     &           - 0.5d0*b*qba*qbb*qbb*(2.0d0*alpha - 1.0d0))

        dx1_dtau = c*p*dz_dtau*2.0d0*z*(1.0d0 - z2)/
     &             ((1.0d0 + z2)**3.0d0)
        dx2_dtau = 2.0d0*c2*qb*dqb_dtau
        dx3_dtau = x3*(dqb_dtau/qb 
     &           + 0.5d0*0.36d0*z*dz_dtau/(x3b*x3b))
        dx5_dtau = 2.0d0*(x5/z)*dz_dtau

        dxnum_dtau= dx1_dtau + dx2_dtau + dx3_dtau + dx5_dtau
        dx_dtau = x7*dxnum_dtau

        Fx = 1.0d0 + kappa - kappa/(1.0d0 + x/kappa)

        ex = ex0*Fx
        nex0 = n*ex0

        dFx_dx = 1.0d0/(1.0d0 + x/kappa)**2.0d0
        dFx_dn = dFx_dx*dx_dn
        fnx = frthrd*ex + nex0*dFx_dn

        dFx_dagr = dFx_dx*dx_dagr
        fdnx = nex0*dFx_dagr

        dFx_dtau = dFx_dx*dx_dtau
        fdtaux = nex0*dFx_dtau

        ec = 0.0d0
        fnc = 0.0d0
        fdnc = 0.0d0
        fdtauc = 0.0d0
 
        xce(i)   = x_parameter*ex     + c_parameter*ec
        fn(i)    = x_parameter*fnx    + c_parameter*fnc
        fdn(i)   = x_parameter*fdnx   + c_parameter*fdnc
        fdtau(i) = x_parameter*fdtaux + c_parameter*fdtauc

      end do
      return
      end
