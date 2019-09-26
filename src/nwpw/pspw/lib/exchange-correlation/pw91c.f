*     ************************************************
*     *                                              *
*     *                gen_PW91_c_rz                 *
*     *                                              *
*     ************************************************
      subroutine gen_PW91_c_rz(tol,rs,zeta,pwc,dpwc_drs,dpwc_dzeta)
      implicit none
*     ***** input *****
      real*8 tol,rs,zeta
*     ***** output *****
      real*8 pwc,dpwc_drs,dpwc_dzeta
*     ***** local declarations *****
      real*8 gammaI,fzzI,zeta2,zeta3,zeta4
      real*8 eu,deu_drs
      real*8 ep,dep_drs
      real*8 alpham,dam_drs
      real*8 gz,hz,dgz,dhz
      real*8 fzeta,df_dzeta
*     ***** constants *****
      real*8 pi,thrd,twthrd,frthrd,fvthrd,etthrd
      parameter (pi     =  3.14159265358979311599d0)
      parameter (thrd   =  1.0d0/3.0d0)
*     ***** PW91 constants *****
      real*8 eps0c1,eps0c2,eps0c3,eps0c4,eps0c5,eps0c6
      real*8 eps1c1,eps1c2,eps1c3,eps1c4,eps1c5,eps1c6
      real*8 epsc1,epsc2,epsc3,epsc4,epsc5,epsc6 
      parameter (eps0c1 = 0.03109070d0)
      parameter (eps0c2 = 0.21370d0)
      parameter (eps0c3 = 7.5957d0)
      parameter (eps0c4 = 3.5876d0)
      parameter (eps0c5 = 1.6382d0)
      parameter (eps0c6 = 0.49294d0)
      parameter (eps1c1 = 0.01554535d0)
      parameter (eps1c2 = 0.20548d0)
      parameter (eps1c3 = 14.1189d0)
      parameter (eps1c4 = 6.1977d0)
      parameter (eps1c5 = 3.3662d0)
      parameter (eps1c6 = 0.62517d0)
      parameter (epsc1  = 0.01688686394d0)
      parameter (epsc2  = 0.11125d0)
      parameter (epsc3  = 10.3570d0)
      parameter (epsc4  = 3.6231d0)
      parameter (epsc5  = 0.88026d0)
      parameter (epsc6  = 0.49671d0)

      fzzI   = 9.0d0*(2.0d0**thrd - 1.0d0)/4.0d0
      gammaI = 1.0d0/(2.0d0*2.0d0**thrd - 2.0d0)
      
      call nwpw_pw91c_zeta(tol,gammaI,zeta,fzeta,df_dzeta)

      call nwpw_pw91c_rs(eps0c1,eps0c2,eps0c3,eps0c4,eps0c5,eps0c6,
     >                  rs,eu,deu_drs)

      call nwpw_pw91c_rs(eps1c1,eps1c2,eps1c3,eps1c4,eps1c5,eps1c6,
     >                  rs,ep,dep_drs)

      call nwpw_pw91c_rs(epsc1,epsc2,epsc3,epsc4,epsc5,epsc6,
     >                  rs,alpham,dam_drs)

      zeta2 = zeta*zeta
      zeta3 = zeta2*zeta
      zeta4 = zeta3*zeta

      gz  = fzeta*zeta4
      hz  = fzzI*(fzeta - gz)
      dgz = df_dzeta*zeta4 + 4.0d0*fzeta*zeta3
      dhz = fzzI*(df_dzeta - dgz)

      pwc        = eu*(1.0d0 - gz) + ep*gz - alpham*hz
      dpwc_drs   = deu_drs*(1.0d0 - gz) + dep_drs*gz - dam_drs*hz
      dpwc_dzeta = (ep - eu)*dgz - alpham*dhz

      return
      end
*     ************************************************
*     *                                              *
*     *                nwpw_pw91c_zeta               *
*     *                                              *
*     ************************************************
      subroutine nwpw_pw91c_zeta(tol,s,zeta,fz,dfdz)
      implicit none
*     ***** input *****
      real*8 tol,s,zeta
*     ***** output *****
      real*8 fz,dfdz
*     ***** local declarations *****
      real*8 small,omz,opz,omz2,opz2,omz_13,opz_13
*     ***** constants *****
      real*8 thrd,frthrd
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)

      small =  tol
      fz    = -2.0d0
      dfdz  =  0.0d0
      omz   =  1.0d0 - zeta
      opz   =  1.0d0 + zeta
      omz2  =  omz**2.0d0
      opz2  =  opz**2.0d0

      if (omz .gt. small) then
        omz_13 = omz**thrd
        fz     = fz   + omz*omz_13
        dfdz   = dfdz - omz_13
      end if

      if (opz .gt. small) then
        opz_13 = opz**thrd
        fz     = fz   + opz*opz_13
        dfdz   = dfdz + opz_13
      end if

      fz    = fz*s
      dfdz  = frthrd*dfdz*s

      return
      end
*     ************************************************
*     *                                              *
*     *                nwpw_pw91c_rs                 *
*     *                                              *
*     ************************************************
      subroutine nwpw_pw91c_rs(a,a1,b1,b2,b3,b4,rs,v,dvdr)
      implicit none
*     ***** input *****
      real*8 a,a1,b1,b2,b3,b4,rs
*     ***** output *****
      real*8 v,dvdr
*     ***** local declarations *****
      real*8 q0,rs_12,rs_32,q1,q2
      real*8 dq0_drs,dq1_drs,dq2_drs
 
      q0    = -2.0d0*a*(1.0d0 + a1*rs)
      rs_12 =  dsqrt(rs)
      rs_32 =  rs*rs_12

      q1 = 2.0d0*a*(b1*rs_12+b2*rs+b3*rs_32+b4*rs*rs)
      q2 = dlog(1.0d0 + 1.0d0/q1)
      v  = q0*q2

      dq0_drs = -2.0d0*a*a1
      dq1_drs =  a*(b1/rs_12 + 2.0d0*b2 + 3.0d0*b3*rs_12 + 4.0d0*b4*rs)
      dq2_drs = -dq1_drs/(q1 + q1**2.0d0)
      dvdr    =  dq0_drs*q2 + q0*dq2_drs

      return
      end
