C> \ingroup nwxc
C> @{
C>
C> \file nwxcm_c_pkzb99.F
C> The nwxcm_c_pkzb99 functional
C>
C> @}
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the nwxcm_c_pkzb99 functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\
C>   {\it t_3} &=& {{1}\over{\rho_\alpha}}\\\\
C>   {\it t_4} &=& {{1}\over{{\it t_1}}}\\\\
C>   {\it t_5} &=& {{1}\over{\rho_\beta}}\\\\
C>   f &=& 1.0\,{\it t_1}\,\left({\it nwxc\_c\_mpbe}\left({\it t_2} 
C>    , \rho_\alpha , \rho_\beta , \sigma_{\alpha\alpha} , 
C>    \sigma_{\alpha\beta} , \sigma_{\beta\beta}\right)\,
C>    \left({{0.53\,\left(0.125\,{\it t_5}\,\sigma_{\beta\beta}
C>    +0.125\,{\it t_3}\,\sigma_{\alpha\alpha}\right)^2}\over{
C>    \left(\tau_\beta+\tau_\alpha\right)^2}}+1.0\right)
C>    -{{0.02390625\,{\it nwxc\_c\_mpbe}\left({\it t_2} , 0.0 
C>    , \rho_\beta , 0.0 , 0.0 , \sigma_{\beta\beta}\right)\,{
C>    \it t_5}\,{\it t_4}\,\sigma_{\beta\beta}^2}\over{\tau_\beta^2}}
C>    -{{0.02390625\,{\it nwxc\_c\_mpbe}\left({\it t_2} 
C>    , \rho_\alpha , 0.0 , \sigma_{\alpha\alpha} , 0.0 , 0.0\right)
C>    \,{\it t_3}\,{\it t_4}\,\sigma_{\alpha\alpha}^2}
C>    \over{\tau_\alpha^2}}\right)\\\\
C>   g &=& 0\\\\
C>   G &=& 0.0\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc-Ds [4,5,6].
C>
C> ### References ###
C>
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82
C>    , 2544 (1999)  , DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.82.2544 ">
C> 10.1103/PhysRevLett.82.2544 </a>
C>
C> [2] Maxima, a computer algebra system,
C> <a href="http://maxima.sourceforge.net/">
C> http://maxima.sourceforge.net/</a>
C>
C> [3] CLISP 2.49 (2010-07-07)
C>
C> [4] autoxc-Ds, revision 27701 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_pkzb99(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,tau,fnc,Amat,Cmat,Mmat) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision sr(NCOL_RHO)
      double precision sg(NCOL_GAMMA)
      double precision st(NCOL_TAU)
      double precision s1f
      double precision s1a(NCOL_AMAT)
      double precision s1c(NCOL_CMAT)
      double precision s1m(NCOL_MMAT)
      double precision s2f
      double precision s2a(NCOL_AMAT)
      double precision s2c(NCOL_CMAT)
      double precision s2m(NCOL_MMAT)
      double precision s3f
      double precision s3a(NCOL_AMAT)
      double precision s3c(NCOL_CMAT)
      double precision s3m(NCOL_MMAT)
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the 
                                    !< density 
      integer ipol                  !< [Input] The number of spin 
                                    !< channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the 
                                    !< functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the 
                                             !< density gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic 
                                             !< energy density 
      double precision fnc(nq)      !< [Output] The value of the 
                                    !< functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative 
                                            !< wrt tau 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            if (taua.gt.tol_rho) then 
              t1 = param(1)
              sr(R_A) = 0.0d+0
              sr(R_B) = rhoa
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammaaa
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call nwxcm_c_mpbe(t1,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,s1c)
              t2 = s1f
              t3 = gammaaa**2
              t4 = 1/rhoa**2
              t5 = 1/taua**2
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call nwxcm_c_mpbe(t1,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,s2c)
              t6 = s2f
              sr(R_T) = 2.0d0*rhoa
              sg(G_TT) = 4.0d0*gammaaa
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call nwxcm_c_mpbe(t1,tol_rho,ipol,1,1.0d0,sr,sg,s3f,s3a,
     +          s3c)
              t7 = s3f
              t8 = 8.28125d-3*t3*t4*t5+1.0d+0
              t9 = s3f*t8-1.1953125d-2*t3*t4*t5*s2f-1.1953125d-2*s1f
     +          *t3*t4*t5
              t10 = 1/rhoa**3
              t11 = 1/taua**3
              fnc(iq) = 2.0d+0*rhoa*t9*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t9*wght+2.0d+0*rhoa*(s3a(D1_RA)
     +          *t8-8.28125d-3*t10*t3*t5*s3f-1.1953125d-2*t3*t4*t5*s2a
     +          (D1_RA)+1.79296875d-2*t10*t3*t5*s2f+5.9765625d-3*t10
     +          *s1f*t3*t5)*wght+Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = 2.0d+0*rhoa*(s3c(D1_GAA)*t8+8.28125d-3
     +          *gammaaa*t4*t5*s3f+5.0d-1*t4*t5*(-2.390625d-2*t3*s2c
     +          (D1_GAA)-4.78125d-2*gammaaa*s2f))*wght+Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = 2.0d+0*rhoa*s3c(D1_GAB)*t8*wght+Cmat
     +          (iq,D1_GAB)
              Mmat(iq,D1_TA) = 2.0d+0*rhoa*(2.390625d-2*t11*t3*t4*s2f
     +          -8.28125d-3*t11*t3*t4*s3f)*wght+Mmat(iq,D1_TA)
            else 
              t1 = param(1)
              t2 = 4.0d-50*rhoa
              sr(R_A) = 0.0d+0
              sr(R_B) = rhoa
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t2
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call nwxcm_c_mpbe(t1,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,s1c)
              t3 = s1f
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t2
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call nwxcm_c_mpbe(t1,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,s2c)
              t4 = s2f
              sr(R_T) = 2.0d0*rhoa
              sg(G_TT) = 4.0d0*t2
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call nwxcm_c_mpbe(t1,tol_rho,ipol,1,1.0d0,sr,sg,s3f,s3a,
     +          s3c)
              t5 = s3f
              t6 = 1.1325d+0*s3f-1.9125d-1*s2f-1.9125d-1*s1f
              t7 = 1/rhoa
              fnc(iq) = 2.0d+0*rhoa*t6*wght+fnc(iq)
              Amat(iq,D1_RA) = 2.0d+0*rhoa*(-9.5625d-2*s2f*t7
     +          +9.5625d-2*s1f*t7+1.1325d+0*s3a(D1_RA)-1.9125d-1*s2a
     +          (D1_RA))*wght+1.0d+0*t6*wght+Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            endif 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,s1c)
              t3 = s1f
              t4 = gammaaa**2
              t5 = 1/rhoa
              t6 = 1/t1
              t7 = 1/taua**2
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammabb
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,s2c)
              t8 = s2f
              t9 = gammabb**2
              t10 = 1/rhob
              t11 = 1/taub**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = gammaaa
              sg(G_AB) = gammaab
              sg(G_BB) = gammabb
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call nwxcm_c_mpbe(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,s3a,
     +          s3c)
              t12 = s3f
              t13 = 1.25d-1*gammaaa*t5+1.25d-1*gammabb*t10
              t14 = t13**2
              t15 = taub+taua
              t16 = 1/t15**2
              t17 = 5.3d-1*t14*t16+1.0d+0
              t18 = -2.390625d-2*t10*t11*t6*s2f*t9-2.390625d-2*s1f*t4
     +          *t5*t6*t7+s3f*t17
              t19 = 1.0d+0*t18*wght
              t20 = 1/t1**2
              t21 = 2.390625d-2*t20*s1f*t4*t5*t7
              t22 = 1/rhoa**2
              t23 = 2.390625d-2*t10*t11*t20*s2f*t9
              t24 = 1/rhob**2
              t25 = -1.0600000000000001d+0*s3f*t14/t15**3
              fnc(iq) = 1.0d+0*t1*t18*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*(-2.390625d-2*s1a(D1_RA)*t4
     +          *t5*t6*t7+2.390625d-2*t22*s1f*t4*t6*t7+t23-1.325d-1
     +          *gammaaa*s3f*t13*t16*t22+t21+s3a(D1_RA)*t17)*wght+t19
     +          +Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*(-2.390625d-2*t10*t11*t6*s2a
     +          (D1_RB)*t9+2.390625d-2*t11*t24*t6*s2f*t9-1.325d-1
     +          *gammabb*s3f*t13*t16*t24+t23+t21+s3a(D1_RB)*t17)*wght
     +          +t19+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = 1.0d+0*t1*((-2.390625d-2*s1c(D1_GAA)
     +          *t4-4.78125d-2*gammaaa*s1f)*t5*t6*t7+1.325d-1*s3f*t13
     +          *t16*t5+s3c(D1_GAA)*t17)*wght+Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = 1.0d+0*t1*s3c(D1_GAB)*t17*wght+Cmat
     +          (iq,D1_GAB)
              Cmat(iq,D1_GBB) = 1.0d+0*t1*(t10*t11*t6*(-2.390625d-2
     +          *s2c(D1_GBB)*t9-4.78125d-2*gammabb*s2f)+s3c(D1_GBB)
     +          *t17+1.325d-1*t10*s3f*t13*t16)*wght+Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = 1.0d+0*t1*(4.78125d-2*s1f*t4*t5*t6/taua
     +          **3+t25)*wght+Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = 1.0d+0*t1*(4.78125d-2*t10*t6*s2f*t9
     +          /taub**3+t25)*wght+Mmat(iq,D1_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhob
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t3
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,s1c)
              t4 = s1f
              t5 = 1/t1
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,s2c)
              t6 = s2f
              t7 = gammaaa**2
              t8 = 1/rhoa
              t9 = 1/taua**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = gammaaa
              sg(G_AB) = 0
              sg(G_BB) = t3
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call nwxcm_c_mpbe(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,s3a,
     +          s3c)
              t10 = s3f
              t11 = 1.25d-1*gammaaa*t8+5.0d-51
              t12 = t11**2
              t13 = taua+1.0d-50
              t14 = 1/t13**2
              t15 = 5.3d-1*t12*t14+1.0d+0
              t16 = -2.390625d-2*t5*s2f*t7*t8*t9-3.825d-1*rhob*s1f*t5
     +          +s3f*t15
              t17 = 1.0d+0*t16*wght
              t18 = 1/t1**2
              t19 = 3.825d-1*rhob*t18*s1f
              t20 = 2.390625d-2*t18*s2f*t7*t8*t9
              t21 = 1/rhoa**2
              fnc(iq) = 1.0d+0*t1*t16*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*(-2.390625d-2*t5*s2a(D1_RA)
     +          *t7*t8*t9+2.390625d-2*t21*t5*s2f*t7*t9-1.325d-1
     +          *gammaaa*s3f*t11*t14*t21+t20+t19+s3a(D1_RA)*t15)*wght
     +          +t17+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*(-3.825d-1*rhob*s1a(D1_RB)*t5
     +          -3.825d-1*s1f*t5+t20+t19+s3a(D1_RB)*t15)*wght+t17+Amat
     +          (iq,D1_RB)
              Cmat(iq,D1_GAA) = 1.0d+0*t1*(t5*(-2.390625d-2*s2c(D1_GAA
     +          )*t7-4.78125d-2*gammaaa*s2f)*t8*t9+1.325d-1*s3f*t11
     +          *t14*t8+s3c(D1_GAA)*t15)*wght+Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = 1.0d+0*t1*(4.78125d-2*t5*s2f*t7*t8/taua
     +          **3-1.0600000000000001d+0*s3f*t12/t13**3)*wght+Mmat
     +          (iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhoa
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t3
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,s1c)
              t4 = s1f
              t5 = 1/t1
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammabb
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,s2c)
              t6 = s2f
              t7 = gammabb**2
              t8 = 1/rhob
              t9 = 1/taub**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = t3
              sg(G_AB) = 0
              sg(G_BB) = gammabb
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call nwxcm_c_mpbe(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,s3a,
     +          s3c)
              t10 = s3f
              t11 = 1.25d-1*gammabb*t8+5.0d-51
              t12 = t11**2
              t13 = taub+1.0d-50
              t14 = 1/t13**2
              t15 = 5.3d-1*t12*t14+1.0d+0
              t16 = -2.390625d-2*t5*s2f*t7*t8*t9-3.825d-1*rhoa*s1f*t5
     +          +s3f*t15
              t17 = 1.0d+0*t16*wght
              t18 = 1/t1**2
              t19 = 3.825d-1*rhoa*t18*s1f
              t20 = 2.390625d-2*t18*s2f*t7*t8*t9
              t21 = 1/rhob**2
              fnc(iq) = 1.0d+0*t1*t16*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*(-3.825d-1*rhoa*s1a(D1_RA)*t5
     +          -3.825d-1*s1f*t5+t20+t19+s3a(D1_RA)*t15)*wght+t17+Amat
     +          (iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*(-2.390625d-2*t5*s2a(D1_RB)
     +          *t7*t8*t9+2.390625d-2*t21*t5*s2f*t7*t9-1.325d-1
     +          *gammabb*s3f*t11*t14*t21+t20+t19+s3a(D1_RB)*t15)*wght
     +          +t17+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = 1.0d+0*t1*(t5*(-2.390625d-2*s2c(D1_GBB
     +          )*t7-4.78125d-2*gammabb*s2f)*t8*t9+1.325d-1*s3f*t11
     +          *t14*t8+s3c(D1_GBB)*t15)*wght+Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = 1.0d+0*t1*(4.78125d-2*t5*s2f*t7*t8/taub
     +          **3-1.0600000000000001d+0*s3f*t12/t13**3)*wght+Mmat
     +          (iq,D1_TB)
            else 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhoa
              t4 = 4.0d-50*rhob
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = t3
              sg(G_AB) = 0
              sg(G_BB) = t4
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call nwxcm_c_mpbe(t2,tol_rho,ipol,1,1.0d0,sr,sg,s1f,s1a,
     +          s1c)
              t5 = s1f
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t3
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,s2c)
              t6 = s2f
              t7 = 1/t1
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t4
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call nwxcm_c_mpbe(t2,tol_rho,2,1,1.0d0,sr,sg,s3f,s3a,s3c)
              t8 = s3f
              t9 = -3.825d-1*rhob*t7*s3f-3.825d-1*rhoa*s2f*t7
     +          +1.1325d+0*s1f
              t10 = 1/t1**2
              t11 = 3.825d-1*rhoa*t10*s2f
              t12 = 3.825d-1*rhob*t10*s3f
              t13 = 1.0d+0*t9*wght
              fnc(iq) = 1.0d+0*t1*t9*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*(-3.825d-1*rhoa*s2a(D1_RA)*t7
     +          -3.825d-1*s2f*t7+1.1325d+0*s1a(D1_RA)+t12+t11)*wght
     +          +t13+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*(-3.825d-1*rhob*t7*s3a(D1_RB)
     +          -3.825d-1*t7*s3f+1.1325d+0*s1a(D1_RB)+t12+t11)*wght
     +          +t13+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            endif 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            else 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            endif 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            else 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            endif 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_pkzb99 
C>
C> \brief Evaluate the nwxcm_c_pkzb99 functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\
C>   {\it t_3} &=& {{1}\over{\rho_\alpha}}\\\\
C>   {\it t_4} &=& {{1}\over{{\it t_1}}}\\\\
C>   {\it t_5} &=& {{1}\over{\rho_\beta}}\\\\
C>   f &=& 1.0\,{\it t_1}\,\left({\it nwxc\_c\_mpbe}\left({\it t_2} 
C>    , \rho_\alpha , \rho_\beta , \sigma_{\alpha\alpha} , 
C>    \sigma_{\alpha\beta} , \sigma_{\beta\beta}\right)\,
C>    \left({{0.53\,\left(0.125\,{\it t_5}\,\sigma_{\beta\beta}
C>    +0.125\,{\it t_3}\,\sigma_{\alpha\alpha}\right)^2}\over{
C>    \left(\tau_\beta+\tau_\alpha\right)^2}}+1.0\right)
C>    -{{0.02390625\,{\it nwxc\_c\_mpbe}\left({\it t_2} , 0.0 
C>    , \rho_\beta , 0.0 , 0.0 , \sigma_{\beta\beta}\right)\,{
C>    \it t_5}\,{\it t_4}\,\sigma_{\beta\beta}^2}\over{\tau_\beta^2}}
C>    -{{0.02390625\,{\it nwxc\_c\_mpbe}\left({\it t_2} 
C>    , \rho_\alpha , 0.0 , \sigma_{\alpha\alpha} , 0.0 , 0.0\right)
C>    \,{\it t_3}\,{\it t_4}\,\sigma_{\alpha\alpha}^2}
C>    \over{\tau_\alpha^2}}\right)\\\\
C>   g &=& 0\\\\
C>   G &=& 0.0\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc-Ds [4,5,6].
C>
C> ### References ###
C>
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82
C>    , 2544 (1999)  , DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.82.2544 ">
C> 10.1103/PhysRevLett.82.2544 </a>
C>
C> [2] Maxima, a computer algebra system,
C> <a href="http://maxima.sourceforge.net/">
C> http://maxima.sourceforge.net/</a>
C>
C> [3] CLISP 2.49 (2010-07-07)
C>
C> [4] autoxc-Ds, revision 27701 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_pkzb99_d2(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision sr(NCOL_RHO)
      double precision sg(NCOL_GAMMA)
      double precision st(NCOL_TAU)
      double precision s1f
      double precision s1a(NCOL_AMAT)
      double precision s1c(NCOL_CMAT)
      double precision s1m(NCOL_MMAT)
      double precision s1a2(NCOL_AMAT2)
      double precision s1c2(NCOL_CMAT2)
      double precision s1m2(NCOL_MMAT2)
      double precision s2f
      double precision s2a(NCOL_AMAT)
      double precision s2c(NCOL_CMAT)
      double precision s2m(NCOL_MMAT)
      double precision s2a2(NCOL_AMAT2)
      double precision s2c2(NCOL_CMAT2)
      double precision s2m2(NCOL_MMAT2)
      double precision s3f
      double precision s3a(NCOL_AMAT)
      double precision s3c(NCOL_CMAT)
      double precision s3m(NCOL_MMAT)
      double precision s3a2(NCOL_AMAT2)
      double precision s3c2(NCOL_CMAT2)
      double precision s3m2(NCOL_MMAT2)
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the 
                                    !< density 
      integer ipol                  !< [Input] The number of spin 
                                    !< channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the 
                                    !< functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the 
                                             !< density gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic 
                                             !< energy density 
      double precision fnc(nq)      !< [Output] The value of the 
                                    !< functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative 
                                            !< wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt tau and 
                                             !< possibly rho or rgamma 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            if (taua.gt.tol_rho) then 
              t1 = param(1)
              sr(R_A) = 0.0d+0
              sr(R_B) = rhoa
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammaaa
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call nwxcm_c_mpbe_d2(t1,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1c,s1c2)
              t2 = s1f
              t3 = gammaaa**2
              t4 = 1/rhoa**2
              t5 = 1/taua**2
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call nwxcm_c_mpbe_d2(t1,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2c,s2c2)
              t6 = s2f
              sr(R_T) = 2.0d0*rhoa
              sg(G_TT) = 4.0d0*gammaaa
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call nwxcm_c_mpbe_d2(t1,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3c,s3c2)
              t7 = s3f
              t8 = 8.28125d-3*t3*t4*t5+1.0d+0
              t9 = s3f*t8-1.1953125d-2*t3*t4*t5*s2f-1.1953125d-2*s1f
     +          *t3*t4*t5
              t10 = 1/rhoa**3
              t11 = -8.28125d-3*t10*t3*t5*s3f
              t12 = s2a(D1_RA)
              t13 = s3a(D1_RA)
              t14 = t13*t8
              t15 = 1.79296875d-2*t10*t3*t5*s2f-1.1953125d-2*t12*t3*t4
     +          *t5+5.9765625d-3*t10*s1f*t3*t5+t14+t11
              t16 = 8.28125d-3*gammaaa*t4*t5*s3f
              t17 = s2c(D1_GAA)
              t18 = -4.78125d-2*gammaaa*s2f-2.390625d-2*t17*t3
              t19 = s3c(D1_GAA)
              t20 = t19*t8
              t21 = s3c(D1_GAB)
              t22 = 1/taua**3
              t23 = -8.28125d-3*t22*t3*t4*s3f
              t24 = 2.390625d-2*t22*t3*t4*s2f+t23
              t25 = 1/rhoa**4
              t26 = -1.6562500000000002d-2*t10*t13*t3*t5
              t27 = s1a(D1_RA)
              t28 = -8.28125d-3*t10*t19*t3*t5
              t29 = 8.28125d-3*gammaaa*t13*t4*t5
              t30 = s1c(D1_GAA)
              t31 = 4.140625d-3*t4*t5*s3f
              t32 = 1.6562500000000002d-2*gammaaa*t19*t4*t5
              t33 = 8.28125d-3*t10*t22*t3*s3f
              t34 = -8.28125d-3*t13*t22*t3*t4
              t35 = -8.28125d-3*gammaaa*t22*t4*s3f
              t36 = -8.28125d-3*t19*t22*t3*t4
              t37 = 1/rhoa
              t38 = 1/taua**4
              fnc(iq) = 2.0d+0*rhoa*t9*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t9*wght+2.0d+0*rhoa*t15*wght
     +          +Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = 2.0d+0*rhoa*(5.0d-1*t18*t4*t5+t20+t16)
     +          *wght+Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = 2.0d+0*rhoa*t21*t8*wght+Cmat(iq,D1_GAB)
              Mmat(iq,D1_TA) = 2.0d+0*rhoa*t24*wght+Mmat(iq,D1_TA)
              Amat2(iq,D2_RA_RA) = 2.0d+0*rhoa*(s3a2(D2_RA_RA)*t8
     +          +2.0703125000000003d-2*t25*t3*t5*s3f-1.1953125d-2*t3
     +          *t4*t5*s2a2(D2_RA_RA)-4.1835937500000003d-2*t25*t3*t5
     +          *s2f-5.9765625d-3*s1f*t25*t3*t5+3.5859375d-2*t10*t12
     +          *t3*t5+t26)*wght+2.0d+0*t15*wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 2.0d+0*rhoa*(s3a2(D2_RA_RB)*t8
     +          +4.140625d-3*t25*t3*t5*s3f-1.1953125d-2*t25*t3*t5*s2f
     +          +5.9765625d-3*t10*t27*t3*t5-1.1953125d-2*s1f*t25*t3*t5
     +          +5.9765625d-3*t10*t12*t3*t5+t26)*wght+1.0d+0*
     +          (5.9765625d-3*t10*t3*t5*s2f-1.1953125d-2*t27*t3*t4*t5
     +          +1.79296875d-2*t10*s1f*t3*t5+t14+t11)*wght+1.0d+0*t15
     +          *wght+Amat2(iq,D2_RA_RB)
              Cmat2(iq,D2_RA_GAA) = 2.0d+0*rhoa*(s3c2(D2_RA_GAA)*t8
     +          -1.2421875d-2*gammaaa*t10*t5*s3f-1.1953125d-2*t3*t4*t5
     +          *s2c2(D2_RA_GAA)+3.5859375d-2*gammaaa*t10*t5*s2f
     +          -2.390625d-2*gammaaa*t12*t4*t5+1.79296875d-2*t10*t17
     +          *t3*t5+t29+t28)*wght+1.0d+0*(-2.390625d-2*gammaaa*t4
     +          *t5*s2f-1.1953125d-2*t17*t3*t4*t5+t20+t16)*wght+Cmat2
     +          (iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = 2.0d+0*rhoa*(s3c2(D2_RA_GAB)*t8
     +          -8.28125d-3*t10*t21*t3*t5)*wght+1.0d+0*t21*t8*wght
     +          +Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = 2.0d+0*rhoa*(s3c2(D2_RA_GBB)*t8
     +          -4.140625d-3*gammaaa*t10*t5*s3f+5.9765625d-3*t10*t3
     +          *t30*t5+1.1953125d-2*gammaaa*t10*s1f*t5+t29+t28)*wght
     +          +1.0d+0*(-1.1953125d-2*t3*t30*t4*t5-2.390625d-2
     +          *gammaaa*s1f*t4*t5+t20+t16)*wght+Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_GAA_GAA) = 2.0d+0*rhoa*(s3c2(D2_GAA_GAA)*t8
     +          +5.0d-1*t4*t5*(-2.390625d-2*t3*s2c2(D2_GAA_GAA)
     +          -4.78125d-2*s2f-9.5625d-2*gammaaa*t17)+t32+t31)*wght
     +          +Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = 2.0d+0*rhoa*(s3c2(D2_GAA_GAB)*t8
     +          +8.28125d-3*gammaaa*t21*t4*t5)*wght+Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = 2.0d+0*rhoa*(s3c2(D2_GAA_GBB)*t8
     +          +t32+t31)*wght+Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = 2.0d+0*rhoa*s3c2(D2_GAB_GAB)*t8
     +          *wght+Cmat2(iq,D2_GAB_GAB)
              Mmat2(iq,D2_RA_TA) = 2.0d+0*rhoa*(-3.5859375d-2*t10*t22
     +          *t3*s2f+2.390625d-2*t12*t22*t3*t4+t34+t33)*wght+1.0d+0
     +          *t24*wght+Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = 1.0d+0*(2.390625d-2*s1f*t22*t3*t4
     +          +t23)*wght+2.0d+0*rhoa*(t34+t33-1.1953125d-2*t10*s1f
     +          *t22*t3)*wght+Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = 2.0d+0*rhoa*(-t18*t22*t4+t36+t35)
     +          *wght+Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = 2.0d+0*rhoa*(t36+t35)*wght+Mmat2
     +          (iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
     +          -1.6562500000000002d-2*t21*t22*t3*t37*wght
              Mmat2(iq,D2_TA_TA) = 2.0d+0*rhoa*(1.2421875d-2*t3*t38*t4
     +          *s3f-7.171875d-2*t3*t38*t4*s2f)*wght+Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = 2.484375d-2*t3*t37*t38*s3f*wght
     +          +Mmat2(iq,D2_TA_TB)
            else 
              t1 = param(1)
              t2 = 4.0d-50*rhoa
              sr(R_A) = 0.0d+0
              sr(R_B) = rhoa
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t2
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call nwxcm_c_mpbe_d2(t1,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1c,s1c2)
              t3 = s1f
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t2
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call nwxcm_c_mpbe_d2(t1,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2c,s2c2)
              t4 = s2f
              sr(R_T) = 2.0d0*rhoa
              sg(G_TT) = 4.0d0*t2
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call nwxcm_c_mpbe_d2(t1,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3c,s3c2)
              t5 = s3f
              t6 = 1.1325d+0*s3f-1.9125d-1*s2f-1.9125d-1*s1f
              t7 = s2a(D1_RA)
              t8 = 1.1325d+0*s3a(D1_RA)
              t9 = 1/rhoa
              t10 = -9.5625d-2*s2f*t9+9.5625d-2*s1f*t9+t8-1.9125d-1*t7
              t11 = 1/rhoa**2
              t12 = s1a(D1_RA)
              fnc(iq) = 2.0d+0*rhoa*t6*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t6*wght+2.0d+0*rhoa*t10*wght
     +          +Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RA_RA) = 2.0d+0*rhoa*(-1.9125d-1*t7*t9
     +          +1.1325d+0*s3a2(D2_RA_RA)-1.9125d-1*s2a2(D2_RA_RA)
     +          +9.5625d-2*t11*s2f-9.5625d-2*t11*s1f)*wght+2.0d+0*t10
     +          *wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 2.0d+0*rhoa*(9.5625d-2*t7*t9
     +          +9.5625d-2*t12*t9+1.1325d+0*s3a2(D2_RA_RB))*wght
     +          +1.0d+0*(9.5625d-2*s2f*t9-9.5625d-2*s1f*t9+t8
     +          -1.9125d-1*t12)*wght+1.0d+0*t10*wght+Amat2(iq,D2_RA_RB)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            endif 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1c,s1c2)
              t3 = s1f
              t4 = gammaaa**2
              t5 = 1/rhoa
              t6 = 1/t1
              t7 = 1/taua**2
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammabb
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2c,s2c2)
              t8 = s2f
              t9 = gammabb**2
              t10 = 1/rhob
              t11 = 1/taub**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = gammaaa
              sg(G_AB) = gammaab
              sg(G_BB) = gammabb
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3c,s3c2)
              t12 = s3f
              t13 = 1.25d-1*gammaaa*t5+1.25d-1*gammabb*t10
              t14 = t13**2
              t15 = taub+taua
              t16 = 1/t15**2
              t17 = 5.3d-1*t14*t16+1.0d+0
              t18 = -2.390625d-2*t10*t11*t6*s2f*t9-2.390625d-2*s1f*t4
     +          *t5*t6*t7+s3f*t17
              t19 = 1.0d+0*t18*wght
              t20 = 1/t1**2
              t21 = 2.390625d-2*t20*s1f*t4*t5*t7
              t22 = 1/rhoa**2
              t23 = s1a(D1_RA)
              t24 = 2.390625d-2*t10*t11*t20*s2f*t9
              t25 = s3a(D1_RA)
              t26 = -2.390625d-2*t23*t4*t5*t6*t7+2.390625d-2*t22*s1f
     +          *t4*t6*t7+t17*t25+t24-1.325d-1*gammaaa*s3f*t13*t16*t22
     +          +t21
              t27 = 1/rhob**2
              t28 = s2a(D1_RB)
              t29 = s3a(D1_RB)
              t30 = 2.390625d-2*t11*t27*t6*s2f*t9-2.390625d-2*t10*t11
     +          *t28*t6*t9+t17*t29-1.325d-1*gammabb*s3f*t13*t16*t27
     +          +t24+t21
              t31 = s1c(D1_GAA)
              t32 = -2.390625d-2*t31*t4-4.78125d-2*gammaaa*s1f
              t33 = 1.325d-1*s3f*t13*t16*t5
              t34 = s3c(D1_GAA)
              t35 = t17*t34
              t36 = s3c(D1_GAB)
              t37 = s2c(D1_GBB)
              t38 = -2.390625d-2*t37*t9-4.78125d-2*gammabb*s2f
              t39 = 1.325d-1*t10*s3f*t13*t16
              t40 = s3c(D1_GBB)
              t41 = t17*t40
              t42 = 1/taua**3
              t43 = 1/t15**3
              t44 = -1.0600000000000001d+0*s3f*t14*t43
              t45 = 4.78125d-2*s1f*t4*t42*t5*t6+t44
              t46 = 1/taub**3
              t47 = 4.78125d-2*t10*t46*t6*s2f*t9+t44
              t48 = 1/t1**3
              t49 = -4.78125d-2*s1f*t4*t48*t5*t7
              t50 = 1/rhoa**3
              t51 = -4.78125d-2*t10*t11*t48*s2f*t9
              t52 = 1/rhob**3
              t53 = 1.0d+0*(-2.390625d-2*t31*t4*t5*t6*t7-4.78125d-2
     +          *gammaaa*s1f*t5*t6*t7+t35+t33)*wght
              t54 = 4.78125d-2*gammaaa*t20*s1f*t5*t7
              t55 = 2.390625d-2*t20*t31*t4*t5*t7
              t56 = 1.0d+0*t17*t36*wght
              t57 = 1.0d+0*(-2.390625d-2*t10*t11*t37*t6*t9-4.78125d-2
     +          *gammabb*t10*t11*t6*s2f+t41+t39)*wght
              t58 = 4.78125d-2*gammabb*t10*t11*t20*s2f
              t59 = 2.390625d-2*t10*t11*t20*t37*t9
              t60 = 1.0d+0*t45*wght
              t61 = -4.78125d-2*t20*s1f*t4*t42*t5
              t62 = -1.0600000000000001d+0*t14*t25*t43
              t63 = 2.65d-1*gammaaa*s3f*t13*t22*t43
              t64 = 1.0d+0*t47*wght
              t65 = -4.78125d-2*t10*t20*t46*s2f*t9
              t66 = -1.0600000000000001d+0*t14*t29*t43
              t67 = 2.65d-1*gammabb*s3f*t13*t27*t43
              t68 = -1.0600000000000001d+0*t14*t34*t43
              t69 = -2.65d-1*s3f*t13*t43*t5
              t70 = -1.0600000000000001d+0*t1*t14*t36*t43*wght
              t71 = -1.0600000000000001d+0*t14*t40*t43
              t72 = -2.65d-1*t10*s3f*t13*t43
              t73 = 1/t15**4
              t74 = 3.18d+0*s3f*t14*t73
              fnc(iq) = 1.0d+0*t1*t18*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t26*wght+t19+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t30*wght+t19+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = 1.0d+0*t1*(t32*t5*t6*t7+t35+t33)*wght
     +          +Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = 1.0d+0*t1*t17*t36*wght+Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = 1.0d+0*t1*(t10*t11*t38*t6+t41+t39)
     +          *wght+Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = 1.0d+0*t1*t45*wght+Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = 1.0d+0*t1*t47*wght+Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*(-4.78125d-2*s1f*t4*t50
     +          *t6*t7-2.390625d-2*s1a2(D2_RA_RA)*t4*t5*t6*t7
     +          +4.78125d-2*t22*t23*t4*t6*t7+4.78125d-2*t20*t23*t4*t5
     +          *t7-4.78125d-2*t20*t22*s1f*t4*t7+t51+2.65d-1*gammaaa
     +          *s3f*t13*t16*t50+t49+1.6562500000000002d-2*s3f*t16*t4
     +          /rhoa**4-2.65d-1*gammaaa*t13*t16*t22*t25+s3a2(D2_RA_RA
     +          )*t17)*wght+2.0d+0*t26*wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*(-2.390625d-2*t11*t20*t27
     +          *s2f*t9+2.390625d-2*t10*t11*t20*t28*t9+2.390625d-2*t20
     +          *t23*t4*t5*t7-2.390625d-2*t20*t22*s1f*t4*t7+t51+t49
     +          -1.325d-1*gammaaa*t13*t16*t22*t29-1.325d-1*gammabb*t13
     +          *t16*t25*t27+1.6562500000000002d-2*gammaaa*gammabb*s3f
     +          *t16*t22*t27+s3a2(D2_RA_RB)*t17)*wght+1.0d+0*t30*wght
     +          +1.0d+0*t26*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*(-2.390625d-2*t10*t11*t6
     +          *s2a2(D2_RB_RB)*t9-4.78125d-2*t11*t52*t6*s2f*t9
     +          -4.78125d-2*t11*t20*t27*s2f*t9+4.78125d-2*t11*t27*t28
     +          *t6*t9+4.78125d-2*t10*t11*t20*t28*t9
     +          +1.6562500000000002d-2*s3f*t16*t9/rhob**4+2.65d-1
     +          *gammabb*s3f*t13*t16*t52+t51+t49-2.65d-1*gammabb*t13
     +          *t16*t27*t29+s3a2(D2_RB_RB)*t17)*wght+2.0d+0*t30*wght
     +          +Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = 1.0d+0*t1*(-2.390625d-2*s1c2
     +          (D2_RA_GAA)*t4*t5*t6*t7-4.78125d-2*gammaaa*t23*t5*t6
     +          *t7+2.390625d-2*t22*t31*t4*t6*t7+4.78125d-2*gammaaa
     +          *t22*s1f*t6*t7+t55+t54-1.6562500000000002d-2*gammaaa
     +          *s3f*t16*t50+1.325d-1*t13*t16*t25*t5-1.325d-1*gammaaa
     +          *t13*t16*t22*t34-1.325d-1*s3f*t13*t16*t22+s3c2
     +          (D2_RA_GAA)*t17)*wght+t53+Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = 1.0d+0*t1*(s3c2(D2_RA_GAB)*t17
     +          -1.325d-1*gammaaa*t13*t16*t22*t36)*wght+t56+Cmat2(iq,
     +          D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = 1.0d+0*t1*(t59+t58-1.325d-1
     +          *gammaaa*t13*t16*t22*t40+1.325d-1*t10*t13*t16*t25
     +          -1.6562500000000002d-2*gammaaa*t10*s3f*t16*t22+s3c2
     +          (D2_RA_GBB)*t17)*wght+t57+Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = 1.0d+0*t1*(t55+t54+1.325d-1*t13
     +          *t16*t29*t5-1.6562500000000002d-2*gammabb*s3f*t16*t27
     +          *t5-1.325d-1*gammabb*t13*t16*t27*t34+s3c2(D2_RB_GAA)
     +          *t17)*wght+t53+Cmat2(iq,D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = 1.0d+0*t1*(s3c2(D2_RB_GAB)*t17
     +          -1.325d-1*gammabb*t13*t16*t27*t36)*wght+t56+Cmat2(iq,
     +          D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = 1.0d+0*t1*(-2.390625d-2*t10*t11*t6
     +          *s2c2(D2_RB_GBB)*t9+2.390625d-2*t11*t27*t37*t6*t9
     +          +4.78125d-2*gammabb*t11*t27*t6*s2f-4.78125d-2*gammabb
     +          *t10*t11*t28*t6+t59+t58-1.6562500000000002d-2*gammabb
     +          *s3f*t16*t52-1.325d-1*gammabb*t13*t16*t27*t40+1.325d-1
     +          *t10*t13*t16*t29-1.325d-1*s3f*t13*t16*t27+s3c2
     +          (D2_RB_GBB)*t17)*wght+t57+Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t1*((-2.390625d-2*s1c2
     +          (D2_GAA_GAA)*t4-9.5625d-2*gammaaa*t31-4.78125d-2*s1f)
     +          *t5*t6*t7+2.65d-1*t13*t16*t34*t5+1.6562500000000002d-2
     +          *s3f*t16*t22+s3c2(D2_GAA_GAA)*t17)*wght+Cmat2(iq,
     +          D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = 1.0d+0*t1*(1.325d-1*t13*t16*t36
     +          *t5+s3c2(D2_GAA_GAB)*t17)*wght+Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = 1.0d+0*t1*(1.325d-1*t13*t16*t40
     +          *t5+1.6562500000000002d-2*t10*s3f*t16*t5+1.325d-1*t10
     +          *t13*t16*t34+s3c2(D2_GAA_GBB)*t17)*wght+Cmat2(iq,
     +          D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = 1.0d+0*t1*s3c2(D2_GAB_GAB)*t17
     +          *wght+Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = 1.325d-1*t1*t10*t13*t16*t36*wght
     +          +1.0d+0*t1*s3c2(D2_GAB_GBB)*t17*wght+Cmat2(iq,
     +          D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t1*(t10*t11*t6*(
     +          -2.390625d-2*s2c2(D2_GBB_GBB)*t9-4.78125d-2*s2f
     +          -9.5625d-2*gammabb*t37)+2.65d-1*t10*t13*t16*t40
     +          +1.6562500000000002d-2*s3f*t16*t27+s3c2(D2_GBB_GBB)
     +          *t17)*wght+Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = 1.0d+0*t1*(t63+t62+t61+4.78125d-2
     +          *t23*t4*t42*t5*t6-4.78125d-2*t22*s1f*t4*t42*t6)*wght
     +          +t60+Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = 1.0d+0*t1*(t65+t63+t62)*wght+t64
     +          +Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = 1.0d+0*t1*(t67+t66+t61)*wght+t60
     +          +Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = 1.0d+0*t1*(-4.78125d-2*t27*t46*t6
     +          *s2f*t9+4.78125d-2*t10*t28*t46*t6*t9+t67+t66+t65)*wght
     +          +t64+Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = 1.0d+0*t1*(t69+t68-2*t32*t42*t5*t6
     +          )*wght+Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = 1.0d+0*t1*(t69+t68)*wght+Mmat2(iq,
     +          D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = t70+Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = t70+Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = 1.0d+0*t1*(t72+t71)*wght+Mmat2(iq,
     +          D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = 1.0d+0*t1*(t72+t71-2*t10*t38*t46
     +          *t6)*wght+Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = 1.0d+0*t1*(t74-1.434375d-1*s1f*t4
     +          *t5*t6/taua**4)*wght+Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = 3.18d+0*t1*s3f*t14*t73*wght+Mmat2
     +          (iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = 1.0d+0*t1*(t74-1.434375d-1*t10*t6
     +          *s2f*t9/taub**4)*wght+Mmat2(iq,D2_TB_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhob
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t3
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1c,s1c2)
              t4 = s1f
              t5 = 1/t1
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2c,s2c2)
              t6 = s2f
              t7 = gammaaa**2
              t8 = 1/rhoa
              t9 = 1/taua**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = gammaaa
              sg(G_AB) = 0
              sg(G_BB) = t3
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3c,s3c2)
              t10 = s3f
              t11 = 1.25d-1*gammaaa*t8+5.0d-51
              t12 = t11**2
              t13 = taua+1.0d-50
              t14 = 1/t13**2
              t15 = 5.3d-1*t12*t14+1.0d+0
              t16 = -2.390625d-2*t5*s2f*t7*t8*t9-3.825d-1*rhob*s1f*t5
     +          +s3f*t15
              t17 = 1.0d+0*t16*wght
              t18 = 1/t1**2
              t19 = 3.825d-1*rhob*t18*s1f
              t20 = 2.390625d-2*t18*s2f*t7*t8*t9
              t21 = 1/rhoa**2
              t22 = s2a(D1_RA)
              t23 = s3a(D1_RA)
              t24 = -2.390625d-2*t22*t5*t7*t8*t9+2.390625d-2*t21*t5
     +          *s2f*t7*t9+t15*t23-1.325d-1*gammaaa*s3f*t11*t14*t21
     +          +t20+t19
              t25 = s1a(D1_RB)
              t26 = s3a(D1_RB)
              t27 = -3.825d-1*s1f*t5-3.825d-1*rhob*t25*t5+t15*t26+t20
     +          +t19
              t28 = s2c(D1_GAA)
              t29 = -2.390625d-2*t28*t7-4.78125d-2*gammaaa*s2f
              t30 = 1.325d-1*s3f*t11*t14*t8
              t31 = s3c(D1_GAA)
              t32 = t15*t31
              t33 = 1/taua**3
              t34 = 1/t13**3
              t35 = 4.78125d-2*t33*t5*s2f*t7*t8-1.0600000000000001d+0
     +          *s3f*t12*t34
              t36 = 1/t1**3
              t37 = -7.65d-1*rhob*t36*s1f
              t38 = -4.78125d-2*t36*s2f*t7*t8*t9
              t39 = 1/rhoa**3
              t40 = 1.0d+0*(-2.390625d-2*t28*t5*t7*t8*t9-4.78125d-2
     +          *gammaaa*t5*s2f*t8*t9+t32+t30)*wght
              t41 = 4.78125d-2*gammaaa*t18*s2f*t8*t9
              t42 = 2.390625d-2*t18*t28*t7*t8*t9
              t43 = 1.0d+0*t35*wght
              t44 = -4.78125d-2*t18*t33*s2f*t7*t8
              fnc(iq) = 1.0d+0*t1*t16*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t24*wght+t17+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t27*wght+t17+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = 1.0d+0*t1*(t29*t5*t8*t9+t32+t30)*wght
     +          +Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = 1.0d+0*t1*t35*wght+Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*(-2.390625d-2*t5*s2a2
     +          (D2_RA_RA)*t7*t8*t9+4.78125d-2*t18*t22*t7*t8*t9
     +          -4.78125d-2*t39*t5*s2f*t7*t9-4.78125d-2*t18*t21*s2f*t7
     +          *t9+4.78125d-2*t21*t22*t5*t7*t9+1.6562500000000002d-2
     +          *s3f*t14*t7/rhoa**4+2.65d-1*gammaaa*s3f*t11*t14*t39
     +          +t38+t37-2.65d-1*gammaaa*t11*t14*t21*t23+s3a2(D2_RA_RA
     +          )*t15)*wght+2.0d+0*t24*wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*(2.390625d-2*t18*t22*t7
     +          *t8*t9-2.390625d-2*t18*t21*s2f*t7*t9+3.825d-1*t18*s1f
     +          +t38+t37-1.325d-1*gammaaa*t11*t14*t21*t26+3.825d-1
     +          *rhob*t18*t25+s3a2(D2_RA_RB)*t15)*wght+1.0d+0*t27*wght
     +          +1.0d+0*t24*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*(-3.825d-1*rhob*s1a2
     +          (D2_RB_RB)*t5-7.65d-1*t25*t5+7.65d-1*t18*s1f+t38+t37
     +          +7.65d-1*rhob*t18*t25+s3a2(D2_RB_RB)*t15)*wght+2.0d+0
     +          *t27*wght+Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = 1.0d+0*t1*(-2.390625d-2*t5*s2c2
     +          (D2_RA_GAA)*t7*t8*t9-4.78125d-2*gammaaa*t22*t5*t8*t9
     +          +2.390625d-2*t21*t28*t5*t7*t9+4.78125d-2*gammaaa*t21
     +          *t5*s2f*t9+1.325d-1*t11*t14*t23*t8+t42+t41
     +          -1.6562500000000002d-2*gammaaa*s3f*t14*t39-1.325d-1
     +          *gammaaa*t11*t14*t21*t31-1.325d-1*s3f*t11*t14*t21+s3c2
     +          (D2_RA_GAA)*t15)*wght+t40+Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = 1.0d+0*t1*(1.325d-1*t11*t14*t26*t8
     +          +t42+t41+s3c2(D2_RB_GAA)*t15)*wght+t40+Cmat2(iq,
     +          D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t1*(t5*(-2.390625d-2*s2c2
     +          (D2_GAA_GAA)*t7-4.78125d-2*s2f-9.5625d-2*gammaaa*t28)
     +          *t8*t9+2.65d-1*t11*t14*t31*t8+1.6562500000000002d-2
     +          *s3f*t14*t21+s3c2(D2_GAA_GAA)*t15)*wght+Cmat2(iq,
     +          D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = 1.0d+0*t1*(4.78125d-2*t22*t33*t5*t7
     +          *t8-4.78125d-2*t21*t33*t5*s2f*t7+t44
     +          -1.0600000000000001d+0*t12*t23*t34+2.65d-1*gammaaa*s3f
     +          *t11*t21*t34)*wght+t43+Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = 1.0d+0*t1*(t44
     +          -1.0600000000000001d+0*t12*t26*t34)*wght+t43+Mmat2(iq,
     +          D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = 1.0d+0*t1*(-2*t29*t33*t5*t8
     +          -2.65d-1*s3f*t11*t34*t8-1.0600000000000001d+0*t12*t31
     +          *t34)*wght+Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = 1.0d+0*t1*(3.18d+0*s3f*t12/t13**4
     +          -1.434375d-1*t5*s2f*t7*t8/taua**4)*wght+Mmat2(iq,
     +          D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhoa
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t3
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1c,s1c2)
              t4 = s1f
              t5 = 1/t1
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammabb
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2c,s2c2)
              t6 = s2f
              t7 = gammabb**2
              t8 = 1/rhob
              t9 = 1/taub**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = t3
              sg(G_AB) = 0
              sg(G_BB) = gammabb
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3c,s3c2)
              t10 = s3f
              t11 = 1.25d-1*gammabb*t8+5.0d-51
              t12 = t11**2
              t13 = taub+1.0d-50
              t14 = 1/t13**2
              t15 = 5.3d-1*t12*t14+1.0d+0
              t16 = -2.390625d-2*t5*s2f*t7*t8*t9-3.825d-1*rhoa*s1f*t5
     +          +s3f*t15
              t17 = 1.0d+0*t16*wght
              t18 = 1/t1**2
              t19 = 3.825d-1*rhoa*t18*s1f
              t20 = s1a(D1_RA)
              t21 = 2.390625d-2*t18*s2f*t7*t8*t9
              t22 = s3a(D1_RA)
              t23 = -3.825d-1*s1f*t5-3.825d-1*rhoa*t20*t5+t15*t22+t21
     +          +t19
              t24 = 1/rhob**2
              t25 = s2a(D1_RB)
              t26 = s3a(D1_RB)
              t27 = -2.390625d-2*t25*t5*t7*t8*t9+2.390625d-2*t24*t5
     +          *s2f*t7*t9+t15*t26-1.325d-1*gammabb*s3f*t11*t14*t24
     +          +t21+t19
              t28 = s2c(D1_GBB)
              t29 = -2.390625d-2*t28*t7-4.78125d-2*gammabb*s2f
              t30 = 1.325d-1*s3f*t11*t14*t8
              t31 = s3c(D1_GBB)
              t32 = t15*t31
              t33 = 1/taub**3
              t34 = 1/t13**3
              t35 = 4.78125d-2*t33*t5*s2f*t7*t8-1.0600000000000001d+0
     +          *s3f*t12*t34
              t36 = 1/t1**3
              t37 = -7.65d-1*rhoa*t36*s1f
              t38 = -4.78125d-2*t36*s2f*t7*t8*t9
              t39 = 1/rhob**3
              t40 = 1.0d+0*(-2.390625d-2*t28*t5*t7*t8*t9-4.78125d-2
     +          *gammabb*t5*s2f*t8*t9+t32+t30)*wght
              t41 = 4.78125d-2*gammabb*t18*s2f*t8*t9
              t42 = 2.390625d-2*t18*t28*t7*t8*t9
              t43 = 1.0d+0*t35*wght
              t44 = -4.78125d-2*t18*t33*s2f*t7*t8
              fnc(iq) = 1.0d+0*t1*t16*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t23*wght+t17+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t27*wght+t17+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = 1.0d+0*t1*(t29*t5*t8*t9+t32+t30)*wght
     +          +Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = 1.0d+0*t1*t35*wght+Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*(-3.825d-1*rhoa*s1a2
     +          (D2_RA_RA)*t5-7.65d-1*t20*t5+7.65d-1*t18*s1f+t38+t37
     +          +7.65d-1*rhoa*t18*t20+s3a2(D2_RA_RA)*t15)*wght+2.0d+0
     +          *t23*wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*(2.390625d-2*t18*t25*t7
     +          *t8*t9-2.390625d-2*t18*t24*s2f*t7*t9+3.825d-1*t18*s1f
     +          +t38+t37-1.325d-1*gammabb*t11*t14*t22*t24+3.825d-1
     +          *rhoa*t18*t20+s3a2(D2_RA_RB)*t15)*wght+1.0d+0*t27*wght
     +          +1.0d+0*t23*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*(-2.390625d-2*t5*s2a2
     +          (D2_RB_RB)*t7*t8*t9+4.78125d-2*t18*t25*t7*t8*t9
     +          -4.78125d-2*t39*t5*s2f*t7*t9-4.78125d-2*t18*t24*s2f*t7
     +          *t9+4.78125d-2*t24*t25*t5*t7*t9+1.6562500000000002d-2
     +          *s3f*t14*t7/rhob**4+2.65d-1*gammabb*s3f*t11*t14*t39
     +          +t38+t37-2.65d-1*gammabb*t11*t14*t24*t26+s3a2(D2_RB_RB
     +          )*t15)*wght+2.0d+0*t27*wght+Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = 1.0d+0*t1*(1.325d-1*t11*t14*t22*t8
     +          +t42+t41+s3c2(D2_RA_GBB)*t15)*wght+t40+Cmat2(iq,
     +          D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = 1.0d+0*t1*(-2.390625d-2*t5*s2c2
     +          (D2_RB_GBB)*t7*t8*t9-4.78125d-2*gammabb*t25*t5*t8*t9
     +          +2.390625d-2*t24*t28*t5*t7*t9+4.78125d-2*gammabb*t24
     +          *t5*s2f*t9+1.325d-1*t11*t14*t26*t8+t42+t41
     +          -1.6562500000000002d-2*gammabb*s3f*t14*t39-1.325d-1
     +          *gammabb*t11*t14*t24*t31-1.325d-1*s3f*t11*t14*t24+s3c2
     +          (D2_RB_GBB)*t15)*wght+t40+Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t1*(t5*(-2.390625d-2*s2c2
     +          (D2_GBB_GBB)*t7-4.78125d-2*s2f-9.5625d-2*gammabb*t28)
     +          *t8*t9+2.65d-1*t11*t14*t31*t8+1.6562500000000002d-2
     +          *s3f*t14*t24+s3c2(D2_GBB_GBB)*t15)*wght+Cmat2(iq,
     +          D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = 1.0d+0*t1*(t44
     +          -1.0600000000000001d+0*t12*t22*t34)*wght+t43+Mmat2(iq,
     +          D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = 1.0d+0*t1*(4.78125d-2*t25*t33*t5*t7
     +          *t8-4.78125d-2*t24*t33*t5*s2f*t7+t44
     +          -1.0600000000000001d+0*t12*t26*t34+2.65d-1*gammabb*s3f
     +          *t11*t24*t34)*wght+t43+Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = 1.0d+0*t1*(-2*t29*t33*t5*t8
     +          -2.65d-1*s3f*t11*t34*t8-1.0600000000000001d+0*t12*t31
     +          *t34)*wght+Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = 1.0d+0*t1*(3.18d+0*s3f*t12/t13**4
     +          -1.434375d-1*t5*s2f*t7*t8/taub**4)*wght+Mmat2(iq,
     +          D2_TB_TB)
            else 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhoa
              t4 = 4.0d-50*rhob
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = t3
              sg(G_AB) = 0
              sg(G_BB) = t4
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,ipol,1,1.0d0,sr,sg,s1f,
     +          s1a,s1a2,s1c,s1c2)
              t5 = s1f
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t3
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2c,s2c2)
              t6 = s2f
              t7 = 1/t1
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t4
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call nwxcm_c_mpbe_d2(t2,tol_rho,2,1,1.0d0,sr,sg,s3f,s3a,
     +          s3a2,s3c,s3c2)
              t8 = s3f
              t9 = -3.825d-1*rhob*t7*s3f-3.825d-1*rhoa*s2f*t7
     +          +1.1325d+0*s1f
              t10 = 1/t1**2
              t11 = 3.825d-1*rhoa*t10*s2f
              t12 = 3.825d-1*rhob*t10*s3f
              t13 = s2a(D1_RA)
              t14 = -3.825d-1*s2f*t7-3.825d-1*rhoa*t13*t7+1.1325d+0
     +          *s1a(D1_RA)+t12+t11
              t15 = 1.0d+0*t9*wght
              t16 = s3a(D1_RB)
              t17 = -3.825d-1*t7*s3f-3.825d-1*rhob*t16*t7+1.1325d+0
     +          *s1a(D1_RB)+t12+t11
              t18 = 1/t1**3
              t19 = -7.65d-1*rhoa*t18*s2f
              t20 = -7.65d-1*rhob*t18*s3f
              fnc(iq) = 1.0d+0*t1*t9*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t14*wght+t15+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t17*wght+t15+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*(-3.825d-1*rhoa*s2a2
     +          (D2_RA_RA)*t7-7.65d-1*t13*t7+7.65d-1*t10*s2f+1.1325d+0
     +          *s1a2(D2_RA_RA)+t20+t19+7.65d-1*rhoa*t10*t13)*wght
     +          +2.0d+0*t14*wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*(3.825d-1*t10*s3f
     +          +3.825d-1*t10*s2f+1.1325d+0*s1a2(D2_RA_RB)+t20+t19
     +          +3.825d-1*rhob*t10*t16+3.825d-1*rhoa*t10*t13)*wght
     +          +1.0d+0*t17*wght+1.0d+0*t14*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*(-3.825d-1*rhob*t7*s3a2
     +          (D2_RB_RB)+7.65d-1*t10*s3f-7.65d-1*t16*t7+1.1325d+0
     +          *s1a2(D2_RB_RB)+t20+t19+7.65d-1*rhob*t10*t16)*wght
     +          +2.0d+0*t17*wght+Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            endif 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            else 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            endif 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            else 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            endif 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_pkzb99_d2 
C>
C> \brief Evaluate the nwxcm_c_pkzb99 functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\
C>   {\it t_3} &=& {{1}\over{\rho_\alpha}}\\\\
C>   {\it t_4} &=& {{1}\over{{\it t_1}}}\\\\
C>   {\it t_5} &=& {{1}\over{\rho_\beta}}\\\\
C>   f &=& 1.0\,{\it t_1}\,\left({\it nwxc\_c\_mpbe}\left({\it t_2} 
C>    , \rho_\alpha , \rho_\beta , \sigma_{\alpha\alpha} , 
C>    \sigma_{\alpha\beta} , \sigma_{\beta\beta}\right)\,
C>    \left({{0.53\,\left(0.125\,{\it t_5}\,\sigma_{\beta\beta}
C>    +0.125\,{\it t_3}\,\sigma_{\alpha\alpha}\right)^2}\over{
C>    \left(\tau_\beta+\tau_\alpha\right)^2}}+1.0\right)
C>    -{{0.02390625\,{\it nwxc\_c\_mpbe}\left({\it t_2} , 0.0 
C>    , \rho_\beta , 0.0 , 0.0 , \sigma_{\beta\beta}\right)\,{
C>    \it t_5}\,{\it t_4}\,\sigma_{\beta\beta}^2}\over{\tau_\beta^2}}
C>    -{{0.02390625\,{\it nwxc\_c\_mpbe}\left({\it t_2} 
C>    , \rho_\alpha , 0.0 , \sigma_{\alpha\alpha} , 0.0 , 0.0\right)
C>    \,{\it t_3}\,{\it t_4}\,\sigma_{\alpha\alpha}^2}
C>    \over{\tau_\alpha^2}}\right)\\\\
C>   g &=& 0\\\\
C>   G &=& 0.0\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc-Ds [4,5,6].
C>
C> ### References ###
C>
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82
C>    , 2544 (1999)  , DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.82.2544 ">
C> 10.1103/PhysRevLett.82.2544 </a>
C>
C> [2] Maxima, a computer algebra system,
C> <a href="http://maxima.sourceforge.net/">
C> http://maxima.sourceforge.net/</a>
C>
C> [3] CLISP 2.49 (2010-07-07)
C>
C> [4] autoxc-Ds, revision 27701 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_pkzb99_d3(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,tau,fnc,Amat,Amat2,Amat3, Cmat,Cmat2,Cmat3,Mmat,Mmat2,
     +   Mmat3) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision sr(NCOL_RHO)
      double precision sg(NCOL_GAMMA)
      double precision st(NCOL_TAU)
      double precision s1f
      double precision s1a(NCOL_AMAT)
      double precision s1c(NCOL_CMAT)
      double precision s1m(NCOL_MMAT)
      double precision s1a2(NCOL_AMAT2)
      double precision s1c2(NCOL_CMAT2)
      double precision s1m2(NCOL_MMAT2)
      double precision s1a3(NCOL_AMAT3)
      double precision s1c3(NCOL_CMAT3)
      double precision s1m3(NCOL_MMAT3)
      double precision s2f
      double precision s2a(NCOL_AMAT)
      double precision s2c(NCOL_CMAT)
      double precision s2m(NCOL_MMAT)
      double precision s2a2(NCOL_AMAT2)
      double precision s2c2(NCOL_CMAT2)
      double precision s2m2(NCOL_MMAT2)
      double precision s2a3(NCOL_AMAT3)
      double precision s2c3(NCOL_CMAT3)
      double precision s2m3(NCOL_MMAT3)
      double precision s3f
      double precision s3a(NCOL_AMAT)
      double precision s3c(NCOL_CMAT)
      double precision s3m(NCOL_MMAT)
      double precision s3a2(NCOL_AMAT2)
      double precision s3c2(NCOL_CMAT2)
      double precision s3m2(NCOL_MMAT2)
      double precision s3a3(NCOL_AMAT3)
      double precision s3c3(NCOL_CMAT3)
      double precision s3m3(NCOL_MMAT3)
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the 
                                    !< density 
      integer ipol                  !< [Input] The number of spin 
                                    !< channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the 
                                    !< functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the 
                                             !< density gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic 
                                             !< energy density 
      double precision fnc(nq)      !< [Output] The value of the 
                                    !< functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative 
                                            !< wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt tau and 
                                             !< possibly rho or rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,NCOL_AMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt rho 
      double precision Cmat3(nq,NCOL_CMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat3(nq,NCOL_MMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt tau and 
                                             !< possibly rho or rgamma 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            if (taua.gt.tol_rho) then 
              t1 = param(1)
              sr(R_A) = 0.0d+0
              sr(R_B) = rhoa
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammaaa
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s1c3,1)
              call nwxcm_c_mpbe_d3(t1,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1a3,s1c,s1c2,s1c3)
              t2 = s1f
              t3 = gammaaa**2
              t4 = 1/rhoa**2
              t5 = 1/taua**2
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s2a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s2c3,1)
              call nwxcm_c_mpbe_d3(t1,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2a3,s2c,s2c2,s2c3)
              t6 = s2f
              sr(R_T) = 2.0d0*rhoa
              sg(G_TT) = 4.0d0*gammaaa
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s3a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s3c3,1)
              call nwxcm_c_mpbe_d3(t1,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3a3,s3c,s3c2,s3c3)
              t7 = s3f
              t8 = 8.28125d-3*t3*t4*t5+1.0d+0
              t9 = s3f*t8-1.1953125d-2*t3*t4*t5*s2f-1.1953125d-2*s1f
     +          *t3*t4*t5
              t10 = 1/rhoa**3
              t11 = -8.28125d-3*t10*t3*t5*s3f
              t12 = s2a(D1_RA)
              t13 = s3a(D1_RA)
              t14 = t13*t8
              t15 = 1.79296875d-2*t10*t3*t5*s2f-1.1953125d-2*t12*t3*t4
     +          *t5+5.9765625d-3*t10*s1f*t3*t5+t14+t11
              t16 = 8.28125d-3*gammaaa*t4*t5*s3f
              t17 = s2c(D1_GAA)
              t18 = -4.78125d-2*gammaaa*s2f-2.390625d-2*t17*t3
              t19 = s3c(D1_GAA)
              t20 = t19*t8
              t21 = s3c(D1_GAB)
              t22 = 1/taua**3
              t23 = -8.28125d-3*t22*t3*t4*s3f
              t24 = 2.390625d-2*t22*t3*t4*s2f+t23
              t25 = 1/rhoa**4
              t26 = -1.6562500000000002d-2*t10*t13*t3*t5
              t27 = s2a2(D2_RA_RA)
              t28 = s3a2(D2_RA_RA)
              t29 = t28*t8+2.0703125000000003d-2*t25*t3*t5*s3f
     +          -4.1835937500000003d-2*t25*t3*t5*s2f-1.1953125d-2*t27
     +          *t3*t4*t5-5.9765625d-3*s1f*t25*t3*t5+3.5859375d-2*t10
     +          *t12*t3*t5+t26
              t30 = s1a(D1_RA)
              t31 = s3a2(D2_RA_RB)
              t32 = t31*t8+4.140625d-3*t25*t3*t5*s3f-1.1953125d-2*t25
     +          *t3*t5*s2f+5.9765625d-3*t10*t3*t30*t5-1.1953125d-2*s1f
     +          *t25*t3*t5+5.9765625d-3*t10*t12*t3*t5+t26
              t33 = -8.28125d-3*t10*t19*t3*t5
              t34 = 8.28125d-3*gammaaa*t13*t4*t5
              t35 = s2c2(D2_RA_GAA)
              t36 = s3c2(D2_RA_GAA)
              t37 = t36*t8-1.2421875d-2*gammaaa*t10*t5*s3f
     +          +3.5859375d-2*gammaaa*t10*t5*s2f-1.1953125d-2*t3*t35
     +          *t4*t5-2.390625d-2*gammaaa*t12*t4*t5+1.79296875d-2*t10
     +          *t17*t3*t5+t34+t33
              t38 = -8.28125d-3*t10*t21*t3*t5
              t39 = s3c2(D2_RA_GAB)
              t40 = t39*t8+t38
              t41 = s1c(D1_GAA)
              t42 = -4.140625d-3*gammaaa*t10*t5*s3f
              t43 = s3c2(D2_RA_GBB)
              t44 = t43*t8+5.9765625d-3*t10*t3*t41*t5+1.1953125d-2
     +          *gammaaa*t10*s1f*t5+t42+t34+t33
              t45 = 4.140625d-3*t4*t5*s3f
              t46 = 1.6562500000000002d-2*gammaaa*t19*t4*t5
              t47 = s2c2(D2_GAA_GAA)
              t48 = -4.78125d-2*s2f-2.390625d-2*t3*t47-9.5625d-2
     +          *gammaaa*t17
              t49 = s3c2(D2_GAA_GAA)
              t50 = t49*t8
              t51 = s3c2(D2_GAA_GAB)
              t52 = t51*t8+8.28125d-3*gammaaa*t21*t4*t5
              t53 = s3c2(D2_GAA_GBB)
              t54 = t53*t8+t46+t45
              t55 = s3c2(D2_GAB_GAB)
              t56 = 8.28125d-3*t10*t22*t3*s3f
              t57 = -8.28125d-3*t13*t22*t3*t4
              t58 = -3.5859375d-2*t10*t22*t3*s2f+t57+t56+2.390625d-2
     +          *t12*t22*t3*t4
              t59 = t57+t56-1.1953125d-2*t10*s1f*t22*t3
              t60 = -8.28125d-3*gammaaa*t22*t4*s3f
              t61 = -8.28125d-3*t19*t22*t3*t4
              t62 = t61+t60
              t63 = 1/rhoa
              t64 = 1/taua**4
              t65 = 1/rhoa**5
              t66 = 2.0703125000000003d-2*t19*t25*t3*t5
              t67 = 8.28125d-3*gammaaa*t28*t4*t5
              t68 = s3c2(D2_RB_GAA)
              t69 = s3c2(D2_RB_GAB)
              t70 = -8.28125d-3*t10*t3*t49*t5
              t71 = 4.140625d-3*t13*t4*t5
              t72 = -8.28125d-3*t10*t3*t5*t51
              t73 = 8.28125d-3*gammaaa*t39*t4*t5
              t74 = s1c2(D2_GAA_GAA)
              t75 = 1.2421875d-2*t19*t4*t5
              t76 = 4.140625d-3*t21*t4*t5
              t77 = 1.6562500000000002d-2*gammaaa*t4*t5*t51
              t78 = -2.0703125000000003d-2*t22*t25*t3*s3f
              t79 = 1.6562500000000002d-2*t10*t13*t22*t3
              t80 = -8.28125d-3*t22*t28*t3*t4
              t81 = 1.2421875d-2*gammaaa*t10*t22*s3f
              t82 = 8.28125d-3*t10*t19*t22*t3
              t83 = -8.28125d-3*gammaaa*t13*t22*t4
              t84 = -8.28125d-3*t22*t3*t36*t4
              t85 = 1.0d+0*t62*wght
              t86 = 2.0d+0*rhoa*(8.28125d-3*t10*t21*t22*t3-8.28125d-3
     +          *t22*t3*t39*t4)*wght-8.28125d-3*t21*t22*t3*t4*wght
              t87 = 4.140625d-3*gammaaa*t10*t22*s3f
              t88 = -8.28125d-3*t22*t3*t4*t43
              t89 = -4.140625d-3*t22*t4*s3f
              t90 = -1.6562500000000002d-2*gammaaa*t19*t22*t4
              t91 = -8.28125d-3*t22*t3*t4*t49
              t92 = 2.0d+0*rhoa*(-8.28125d-3*t22*t3*t4*t51-8.28125d-3
     +          *gammaaa*t21*t22*t4)*wght
              t93 = 1/taua**5
              fnc(iq) = 2.0d+0*rhoa*t9*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t9*wght+2.0d+0*rhoa*t15*wght
     +          +Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = 2.0d+0*rhoa*(5.0d-1*t18*t4*t5+t20+t16)
     +          *wght+Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = 2.0d+0*rhoa*t21*t8*wght+Cmat(iq,D1_GAB)
              Mmat(iq,D1_TA) = 2.0d+0*rhoa*t24*wght+Mmat(iq,D1_TA)
              Amat2(iq,D2_RA_RA) = 2.0d+0*rhoa*t29*wght+2.0d+0*t15
     +          *wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*(5.9765625d-3*t10*t3*t5*s2f
     +          -1.1953125d-2*t3*t30*t4*t5+1.79296875d-2*t10*s1f*t3*t5
     +          +t14+t11)*wght+2.0d+0*rhoa*t32*wght+1.0d+0*t15*wght
     +          +Amat2(iq,D2_RA_RB)
              Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-2.390625d-2*gammaaa*t4*t5
     +          *s2f-1.1953125d-2*t17*t3*t4*t5+t20+t16)*wght+2.0d+0
     +          *rhoa*t37*wght+Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = 1.0d+0*t21*t8*wght+2.0d+0*rhoa*t40
     +          *wght+Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = 1.0d+0*(-1.1953125d-2*t3*t4*t41*t5
     +          -2.390625d-2*gammaaa*s1f*t4*t5+t20+t16)*wght+2.0d+0
     +          *rhoa*t44*wght+Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_GAA_GAA) = 2.0d+0*rhoa*(t50+5.0d-1*t4*t48*t5
     +          +t46+t45)*wght+Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = 2.0d+0*rhoa*t52*wght+Cmat2(iq,
     +          D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = 2.0d+0*rhoa*t54*wght+Cmat2(iq,
     +          D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = 2.0d+0*rhoa*t55*t8*wght+Cmat2(iq,
     +          D2_GAB_GAB)
              Mmat2(iq,D2_RA_TA) = 2.0d+0*rhoa*t58*wght+1.0d+0*t24
     +          *wght+Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = 2.0d+0*rhoa*t59*wght+1.0d+0*
     +          (2.390625d-2*s1f*t22*t3*t4+t23)*wght+Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = 2.0d+0*rhoa*(t61+t60-t18*t22*t4)
     +          *wght+Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = 2.0d+0*rhoa*t62*wght+Mmat2(iq,
     +          D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
     +          -1.6562500000000002d-2*t21*t22*t3*t63*wght
              Mmat2(iq,D2_TA_TA) = 2.0d+0*rhoa*(1.2421875d-2*t3*t4*t64
     +          *s3f-7.171875d-2*t3*t4*s2f*t64)*wght+Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = 2.484375d-2*t3*t63*t64*s3f*wght
     +          +Mmat2(iq,D2_TA_TB)
              Amat3(iq,D3_RA_RA_RA) = 2.0d+0*rhoa*(s3a3(D3_RA_RA_RA)
     +          *t8-7.453125d-2*t3*t5*t65*s3f+1.3447265625d-1*t3*t5
     +          *s2f*t65+8.96484375d-3*s1f*t3*t5*t65-1.1953125d-2*t3
     +          *t4*t5*s2a3(D3_RA_RA_RA)-2.484375d-2*t10*t28*t3*t5
     +          +5.37890625d-2*t10*t27*t3*t5+6.2109375d-2*t13*t25*t3
     +          *t5-1.2550781249999998d-1*t12*t25*t3*t5)*wght+3.0d+0
     +          *t29*wght+Amat3(iq,D3_RA_RA_RA)
              Amat3(iq,D3_RA_RA_RB) = 2.0d+0*rhoa*(s3a3(D3_RA_RA_RB)
     +          *t8-8.28125d-3*t3*t5*t65*s3f+3.287109375d-2*t3*t5*s2f
     +          *t65+1.494140625d-2*s1f*t3*t5*t65
     +          -1.6562500000000002d-2*t10*t3*t31*t5-5.9765625d-3*t25
     +          *t3*t30*t5-8.28125d-3*t10*t28*t3*t5+5.9765625d-3*t10
     +          *t27*t3*t5+2.8984375d-2*t13*t25*t3*t5-2.390625d-2*t12
     +          *t25*t3*t5)*wght+2.0d+0*t32*wght+1.0d+0*t29*wght+Amat3
     +          (iq,D3_RA_RA_RB)
              Cmat3(iq,D3_RA_RA_GAA) = 2.0d+0*rhoa*(s3c3(D3_RA_RA_GAA)
     +          *t8+3.3125000000000004d-2*gammaaa*t25*t5*s3f+t67+t66
     +          -1.1953125d-2*t3*t4*t5*s2c3(D3_RA_RA_GAA)
     +          -8.367187500000001d-2*gammaaa*t25*t5*s2f-2.390625d-2
     +          *gammaaa*t27*t4*t5-1.6562500000000002d-2*t10*t3*t36*t5
     +          +3.5859375d-2*t10*t3*t35*t5-4.1835937500000003d-2*t17
     +          *t25*t3*t5-2.484375d-2*gammaaa*t10*t13*t5+7.171875d-2
     +          *gammaaa*t10*t12*t5)*wght+2.0d+0*t37*wght+Cmat3(iq,
     +          D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_RA_GAB) = 2.0d+0*rhoa*(s3c3(D3_RA_RA_GAB)
     +          *t8-1.6562500000000002d-2*t10*t3*t39*t5
     +          +2.0703125000000003d-2*t21*t25*t3*t5)*wght+2.0d+0*t40
     +          *wght+Cmat3(iq,D3_RA_RA_GAB)
              Cmat3(iq,D3_RA_RA_GBB) = 2.0d+0*rhoa*(s3c3(D3_RA_RA_GBB)
     +          *t8+8.28125d-3*gammaaa*t25*t5*s3f+t67+t66
     +          -1.6562500000000002d-2*t10*t3*t43*t5-5.9765625d-3*t25
     +          *t3*t41*t5-1.1953125d-2*gammaaa*s1f*t25*t5-8.28125d-3
     +          *gammaaa*t10*t13*t5)*wght+2.0d+0*t44*wght+Cmat3(iq,
     +          D3_RA_RA_GBB)
              Cmat3(iq,D3_RA_RB_GAA) = 2.0d+0*rhoa*(s3c3(D3_RA_RB_GAA)
     +          *t8+4.140625d-3*gammaaa*t25*t5*s3f-8.28125d-3*t10*t3
     +          *t5*t68-2.390625d-2*gammaaa*t25*t5*s2f+8.28125d-3
     +          *gammaaa*t31*t4*t5-8.28125d-3*t10*t3*t36*t5
     +          +5.9765625d-3*t10*t3*t35*t5+4.140625d-3*t19*t25*t3*t5
     +          -1.1953125d-2*t17*t25*t3*t5-1.6562500000000002d-2
     +          *gammaaa*t10*t13*t5+1.1953125d-2*gammaaa*t10*t12*t5)
     +          *wght+1.0d+0*(t68*t8+1.1953125d-2*gammaaa*t10*t5*s2f
     +          +5.9765625d-3*t10*t17*t3*t5+t42+t34+t33)*wght+1.0d+0
     +          *t37*wght+Cmat3(iq,D3_RA_RB_GAA)
              Cmat3(iq,D3_RA_RB_GAB) = 2.0d+0*rhoa*(s3c3(D3_RA_RB_GAB)
     +          *t8-8.28125d-3*t10*t3*t5*t69-8.28125d-3*t10*t3*t39*t5
     +          +4.140625d-3*t21*t25*t3*t5)*wght+1.0d+0*(t69*t8+t38)
     +          *wght+1.0d+0*t40*wght+Cmat3(iq,D3_RA_RB_GAB)
              Cmat3(iq,D3_RA_GAA_GAA) = 2.0d+0*rhoa*(s3c3
     +          (D3_RA_GAA_GAA)*t8+t71+t70-8.28125d-3*t10*t5*s3f
     +          -1.1953125d-2*t3*t4*t5*s2c3(D3_RA_GAA_GAA)
     +          +3.5859375d-2*t10*t5*s2f+1.79296875d-2*t10*t3*t47*t5
     +          +1.6562500000000002d-2*gammaaa*t36*t4*t5-4.78125d-2
     +          *gammaaa*t35*t4*t5-2.390625d-2*t12*t4*t5-2.484375d-2
     +          *gammaaa*t10*t19*t5+7.171875d-2*gammaaa*t10*t17*t5)
     +          *wght+1.0d+0*(-2.390625d-2*t4*t5*s2f+t50-1.1953125d-2
     +          *t3*t4*t47*t5-4.78125d-2*gammaaa*t17*t4*t5+t46+t45)
     +          *wght+Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_RA_GAA_GAB) = 2.0d+0*rhoa*(s3c3
     +          (D3_RA_GAA_GAB)*t8+t73+t72-1.2421875d-2*gammaaa*t10
     +          *t21*t5)*wght+1.0d+0*t52*wght+Cmat3(iq,D3_RA_GAA_GAB)
              Cmat3(iq,D3_RA_GAA_GBB) = 2.0d+0*rhoa*(s3c3
     +          (D3_RA_GAA_GBB)*t8+t71-4.140625d-3*t10*t5*s3f
     +          -8.28125d-3*t10*t3*t5*t53+8.28125d-3*gammaaa*t4*t43*t5
     +          +8.28125d-3*gammaaa*t36*t4*t5-1.6562500000000002d-2
     +          *gammaaa*t10*t19*t5)*wght+1.0d+0*t54*wght+Cmat3(iq,
     +          D3_RA_GAA_GBB)
              Cmat3(iq,D3_RA_GAB_GAB) = 2.0d+0*rhoa*(s3c3
     +          (D3_RA_GAB_GAB)*t8-8.28125d-3*t10*t3*t5*t55)*wght
     +          +1.0d+0*t55*t8*wght+Cmat3(iq,D3_RA_GAB_GAB)
              Cmat3(iq,D3_RA_GAB_GBB) = 2.0d+0*rhoa*(s3c3
     +          (D3_RA_GAB_GBB)*t8+t73+t72-4.140625d-3*gammaaa*t10*t21
     +          *t5)*wght+1.0d+0*t51*t8*wght+8.28125d-3*gammaaa*t21*t4
     +          *t5*wght+Cmat3(iq,D3_RA_GAB_GBB)
              Cmat3(iq,D3_RA_GBB_GBB) = 2.0d+0*rhoa*(s3c3
     +          (D3_RA_GBB_GBB)*t8+5.9765625d-3*t10*t3*t5*t74+t71+t70
     +          +1.6562500000000002d-2*gammaaa*t4*t43*t5+2.390625d-2
     +          *gammaaa*t10*t41*t5+1.1953125d-2*t10*s1f*t5-8.28125d-3
     +          *gammaaa*t10*t19*t5)*wght+1.0d+0*(-1.1953125d-2*t3*t4
     +          *t5*t74+t50-4.78125d-2*gammaaa*t4*t41*t5-2.390625d-2
     +          *s1f*t4*t5+t46+t45)*wght+Cmat3(iq,D3_RA_GBB_GBB)
              Cmat3(iq,D3_GAA_GAA_GAA) = 2.0d+0*rhoa*(s3c3
     +          (D3_GAA_GAA_GAA)*t8+t75+5.0d-1*t4*t5*(-2.390625d-2*t3
     +          *s2c3(D3_GAA_GAA_GAA)-1.434375d-1*gammaaa*t47
     +          -1.434375d-1*t17)+2.484375d-2*gammaaa*t4*t49*t5)*wght
     +          +Cmat3(iq,D3_GAA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAB) = 2.0d+0*rhoa*(s3c3
     +          (D3_GAA_GAA_GAB)*t8+t77+t76)*wght+Cmat3(iq,
     +          D3_GAA_GAA_GAB)
              Cmat3(iq,D3_GAA_GAA_GBB) = 2.0d+0*rhoa*(s3c3
     +          (D3_GAA_GAA_GBB)*t8+t75+1.6562500000000002d-2*gammaaa
     +          *t4*t5*t53+8.28125d-3*gammaaa*t4*t49*t5)*wght+Cmat3
     +          (iq,D3_GAA_GAA_GBB)
              Cmat3(iq,D3_GAA_GAB_GAB) = 2.0d+0*rhoa*(s3c3
     +          (D3_GAA_GAB_GAB)*t8+8.28125d-3*gammaaa*t4*t5*t55)*wght
     +          +Cmat3(iq,D3_GAA_GAB_GAB)
              Cmat3(iq,D3_GAA_GAB_GBB) = 2.0d+0*rhoa*(s3c3
     +          (D3_GAA_GAB_GBB)*t8+t77+t76)*wght+Cmat3(iq,
     +          D3_GAA_GAB_GBB)
              Cmat3(iq,D3_GAB_GAB_GAB) = 2.0d+0*rhoa*s3c3
     +          (D3_GAB_GAB_GAB)*t8*wght+Cmat3(iq,D3_GAB_GAB_GAB)
              Mmat3(iq,D3_RA_RA_TA) = 2.0d+0*rhoa*(t80+t79+t78
     +          +8.367187500000001d-2*t22*t25*t3*s2f+2.390625d-2*t22
     +          *t27*t3*t4-7.171875d-2*t10*t12*t22*t3)*wght+2.0d+0*t58
     +          *wght+Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = 2.0d+0*rhoa*(t80+t79+t78
     +          +1.1953125d-2*s1f*t22*t25*t3)*wght+2.0d+0*t59*wght
     +          +Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_RB_TA) = 2.0d+0*rhoa*(t79-4.140625d-3*t22
     +          *t25*t3*s3f+2.390625d-2*t22*t25*t3*s2f-8.28125d-3*t22
     +          *t3*t31*t4-1.1953125d-2*t10*t12*t22*t3)*wght+1.0d+0*(
     +          -1.1953125d-2*t10*t22*t3*s2f+t57+t56)*wght+1.0d+0*t58
     +          *wght+Mmat3(iq,D3_RA_RB_TA)
              Mmat3(iq,D3_RA_GAA_TA) = 2.0d+0*rhoa*(t84+t83+t82+t81
     +          -7.171875d-2*gammaaa*t10*t22*s2f+2.390625d-2*t22*t3
     +          *t35*t4+4.78125d-2*gammaaa*t12*t22*t4-3.5859375d-2*t10
     +          *t17*t22*t3)*wght+1.0d+0*(t61+t60+4.78125d-2*gammaaa
     +          *t22*t4*s2f+2.390625d-2*t17*t22*t3*t4)*wght+Mmat3(iq,
     +          D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = 2.0d+0*rhoa*(t84+t83+t82+t81)
     +          *wght+t85+Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_RA_GAB_TA) = t86+Mmat3(iq,D3_RA_GAB_TA)
              Mmat3(iq,D3_RA_GAB_TB) = t86+Mmat3(iq,D3_RA_GAB_TB)
              Mmat3(iq,D3_RA_GBB_TA) = 2.0d+0*rhoa*(t88+t87+t83+t82)
     +          *wght+t85+Mmat3(iq,D3_RA_GBB_TA)
              Mmat3(iq,D3_RA_GBB_TB) = 2.0d+0*rhoa*(t88+t87+t83+t82
     +          -1.1953125d-2*t10*t22*t3*t41-2.390625d-2*gammaaa*t10
     +          *s1f*t22)*wght+1.0d+0*(t61+t60+2.390625d-2*t22*t3*t4
     +          *t41+4.78125d-2*gammaaa*s1f*t22*t4)*wght+Mmat3(iq,
     +          D3_RA_GBB_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = 2.0d+0*rhoa*(t91+t90+t89-t22
     +          *t4*t48)*wght+Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = 2.0d+0*rhoa*(t91+t90+t89)*wght
     +          +Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_GAA_GAB_TA) = t92+Mmat3(iq,D3_GAA_GAB_TA)
              Mmat3(iq,D3_GAA_GAB_TB) = t92+Mmat3(iq,D3_GAA_GAB_TB)
              Mmat3(iq,D3_GAA_GBB_TA) = 2.0d+0*rhoa*(t90+t89
     +          -8.28125d-3*t22*t3*t4*t53)*wght+Mmat3(iq,D3_GAA_GBB_TA)
              Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
     +          -1.6562500000000002d-2*t22*t3*t55*t63*wght
              Mmat3(iq,D3_TA_TA_TA) = 2.0d+0*rhoa*(2.86875d-1*t3*t4
     +          *s2f*t93-2.484375d-2*t3*t4*s3f*t93)*wght+Mmat3(iq,
     +          D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)-4.96875d-2
     +          *t3*t63*s3f*t93*wght
            else 
              t1 = param(1)
              t2 = 4.0d-50*rhoa
              sr(R_A) = 0.0d+0
              sr(R_B) = rhoa
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t2
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s1c3,1)
              call nwxcm_c_mpbe_d3(t1,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1a3,s1c,s1c2,s1c3)
              t3 = s1f
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t2
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s2a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s2c3,1)
              call nwxcm_c_mpbe_d3(t1,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2a3,s2c,s2c2,s2c3)
              t4 = s2f
              sr(R_T) = 2.0d0*rhoa
              sg(G_TT) = 4.0d0*t2
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s3a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s3c3,1)
              call nwxcm_c_mpbe_d3(t1,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3a3,s3c,s3c2,s3c3)
              t5 = s3f
              t6 = 1.1325d+0*s3f-1.9125d-1*s2f-1.9125d-1*s1f
              t7 = s2a(D1_RA)
              t8 = 1.1325d+0*s3a(D1_RA)
              t9 = 1/rhoa
              t10 = -9.5625d-2*s2f*t9+9.5625d-2*s1f*t9+t8-1.9125d-1*t7
              t11 = s2a2(D2_RA_RA)
              t12 = 1/rhoa**2
              t13 = -1.9125d-1*t7*t9+1.1325d+0*s3a2(D2_RA_RA)
     +          +9.5625d-2*t12*s2f-9.5625d-2*t12*s1f-1.9125d-1*t11
              t14 = s1a(D1_RA)
              t15 = 9.5625d-2*t7*t9+9.5625d-2*t14*t9+1.1325d+0*s3a2
     +          (D2_RA_RB)
              t16 = 1/rhoa**3
              fnc(iq) = 2.0d+0*rhoa*t6*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t6*wght+2.0d+0*rhoa*t10*wght
     +          +Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RA_RA) = 2.0d+0*rhoa*t13*wght+2.0d+0*t10
     +          *wght+Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*(9.5625d-2*s2f*t9-9.5625d-2
     +          *s1f*t9+t8-1.9125d-1*t14)*wght+2.0d+0*rhoa*t15*wght
     +          +1.0d+0*t10*wght+Amat2(iq,D2_RA_RB)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Amat3(iq,D3_RA_RA_RA) = 2.0d+0*rhoa*(-2.86875d-1*t11*t9
     +          +2.86875d-1*t12*t7+1.1325d+0*s3a3(D3_RA_RA_RA)
     +          -1.9125d-1*s2a3(D3_RA_RA_RA)-1.434375d-1*t16*s2f
     +          +1.434375d-1*t16*s1f)*wght+3.0d+0*t13*wght+Amat3(iq,
     +          D3_RA_RA_RA)
              Amat3(iq,D3_RA_RA_RB) = 2.0d+0*rhoa*(9.5625d-2*t11*t9
     +          +1.1325d+0*s3a3(D3_RA_RA_RB)-4.781250000000001d-2*t16
     +          *s2f+4.781249999999999d-2*t16*s1f-9.5625d-2*t12*t14)
     +          *wght+2.0d+0*t15*wght+1.0d+0*t13*wght+Amat3(iq,
     +          D3_RA_RA_RB)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
              Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
              Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
              Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
              Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
              Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
              Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
              Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
              Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
              Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
              Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
              Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
              Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
              Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
              Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
              Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
              Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
              Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            endif 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s1c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1a3,s1c,s1c2,s1c3)
              t3 = s1f
              t4 = gammaaa**2
              t5 = 1/rhoa
              t6 = 1/t1
              t7 = 1/taua**2
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammabb
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s2a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s2c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2a3,s2c,s2c2,s2c3)
              t8 = s2f
              t9 = gammabb**2
              t10 = 1/rhob
              t11 = 1/taub**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = gammaaa
              sg(G_AB) = gammaab
              sg(G_BB) = gammabb
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s3a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s3c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3a3,s3c,s3c2,s3c3)
              t12 = s3f
              t13 = 1.25d-1*gammaaa*t5+1.25d-1*gammabb*t10
              t14 = t13**2
              t15 = taub+taua
              t16 = 1/t15**2
              t17 = 5.3d-1*t14*t16+1.0d+0
              t18 = -2.390625d-2*t10*t11*t6*s2f*t9-2.390625d-2*s1f*t4
     +          *t5*t6*t7+s3f*t17
              t19 = 1.0d+0*t18*wght
              t20 = 1/t1**2
              t21 = 2.390625d-2*t20*s1f*t4*t5*t7
              t22 = 1/rhoa**2
              t23 = s1a(D1_RA)
              t24 = 2.390625d-2*t10*t11*t20*s2f*t9
              t25 = s3a(D1_RA)
              t26 = -2.390625d-2*t23*t4*t5*t6*t7+2.390625d-2*t22*s1f
     +          *t4*t6*t7+t17*t25+t24-1.325d-1*gammaaa*s3f*t13*t16*t22
     +          +t21
              t27 = 1/rhob**2
              t28 = s2a(D1_RB)
              t29 = s3a(D1_RB)
              t30 = 2.390625d-2*t11*t27*t6*s2f*t9-2.390625d-2*t10*t11
     +          *t28*t6*t9+t17*t29-1.325d-1*gammabb*s3f*t13*t16*t27
     +          +t24+t21
              t31 = s1c(D1_GAA)
              t32 = -2.390625d-2*t31*t4-4.78125d-2*gammaaa*s1f
              t33 = 1.325d-1*s3f*t13*t16*t5
              t34 = s3c(D1_GAA)
              t35 = t17*t34
              t36 = s3c(D1_GAB)
              t37 = s2c(D1_GBB)
              t38 = -2.390625d-2*t37*t9-4.78125d-2*gammabb*s2f
              t39 = 1.325d-1*t10*s3f*t13*t16
              t40 = s3c(D1_GBB)
              t41 = t17*t40
              t42 = 1/taua**3
              t43 = 1/t15**3
              t44 = -1.0600000000000001d+0*s3f*t14*t43
              t45 = 4.78125d-2*s1f*t4*t42*t5*t6+t44
              t46 = 1/taub**3
              t47 = 4.78125d-2*t10*t46*t6*s2f*t9+t44
              t48 = 1/t1**3
              t49 = -4.78125d-2*s1f*t4*t48*t5*t7
              t50 = 1/rhoa**3
              t51 = s1a2(D2_RA_RA)
              t52 = -4.78125d-2*t10*t11*t48*s2f*t9
              t53 = 1/rhoa**4
              t54 = s3a2(D2_RA_RA)
              t55 = -2.390625d-2*t4*t5*t51*t6*t7-4.78125d-2*s1f*t4*t50
     +          *t6*t7+4.78125d-2*t22*t23*t4*t6*t7+4.78125d-2*t20*t23
     +          *t4*t5*t7-4.78125d-2*t20*t22*s1f*t4*t7+t17*t54
     +          +1.6562500000000002d-2*s3f*t16*t4*t53+t52+2.65d-1
     +          *gammaaa*s3f*t13*t16*t50+t49-2.65d-1*gammaaa*t13*t16
     +          *t22*t25
              t56 = s3a2(D2_RA_RB)
              t57 = -2.390625d-2*t11*t20*t27*s2f*t9+2.390625d-2*t10
     +          *t11*t20*t28*t9+2.390625d-2*t20*t23*t4*t5*t7
     +          -2.390625d-2*t20*t22*s1f*t4*t7+t17*t56+t52+t49
     +          -1.325d-1*gammaaa*t13*t16*t22*t29-1.325d-1*gammabb*t13
     +          *t16*t25*t27+1.6562500000000002d-2*gammaaa*gammabb*s3f
     +          *t16*t22*t27
              t58 = 1/rhob**3
              t59 = s2a2(D2_RB_RB)
              t60 = 1/rhob**4
              t61 = s3a2(D2_RB_RB)
              t62 = -4.78125d-2*t11*t58*t6*s2f*t9-4.78125d-2*t11*t20
     +          *t27*s2f*t9+1.6562500000000002d-2*s3f*t16*t60*t9
     +          -2.390625d-2*t10*t11*t59*t6*t9+4.78125d-2*t11*t27*t28
     +          *t6*t9+4.78125d-2*t10*t11*t20*t28*t9+t17*t61+2.65d-1
     +          *gammabb*s3f*t13*t16*t58+t52+t49-2.65d-1*gammabb*t13
     +          *t16*t27*t29
              t63 = 1.0d+0*(-2.390625d-2*t31*t4*t5*t6*t7-4.78125d-2
     +          *gammaaa*s1f*t5*t6*t7+t35+t33)*wght
              t64 = 4.78125d-2*gammaaa*t20*s1f*t5*t7
              t65 = 2.390625d-2*t20*t31*t4*t5*t7
              t66 = s1c2(D2_RA_GAA)
              t67 = s3c2(D2_RA_GAA)
              t68 = -2.390625d-2*t4*t5*t6*t66*t7-4.78125d-2*gammaaa
     +          *t23*t5*t6*t7+2.390625d-2*t22*t31*t4*t6*t7+4.78125d-2
     +          *gammaaa*t22*s1f*t6*t7+t17*t67+t65+t64
     +          -1.6562500000000002d-2*gammaaa*s3f*t16*t50+1.325d-1
     +          *t13*t16*t25*t5-1.325d-1*gammaaa*t13*t16*t22*t34
     +          -1.325d-1*s3f*t13*t16*t22
              t69 = 1.0d+0*t17*t36*wght
              t70 = s3c2(D2_RA_GAB)
              t71 = t17*t70-1.325d-1*gammaaa*t13*t16*t22*t36
              t72 = 1.0d+0*(-2.390625d-2*t10*t11*t37*t6*t9-4.78125d-2
     +          *gammabb*t10*t11*t6*s2f+t41+t39)*wght
              t73 = 4.78125d-2*gammabb*t10*t11*t20*s2f
              t74 = 2.390625d-2*t10*t11*t20*t37*t9
              t75 = s3c2(D2_RA_GBB)
              t76 = t17*t75+t74+t73-1.325d-1*gammaaa*t13*t16*t22*t40
     +          +1.325d-1*t10*t13*t16*t25-1.6562500000000002d-2
     +          *gammaaa*t10*s3f*t16*t22
              t77 = s3c2(D2_RB_GAA)
              t78 = t17*t77+t65+t64+1.325d-1*t13*t16*t29*t5
     +          -1.6562500000000002d-2*gammabb*s3f*t16*t27*t5-1.325d-1
     +          *gammabb*t13*t16*t27*t34
              t79 = s3c2(D2_RB_GAB)
              t80 = t17*t79-1.325d-1*gammabb*t13*t16*t27*t36
              t81 = s2c2(D2_RB_GBB)
              t82 = s3c2(D2_RB_GBB)
              t83 = -2.390625d-2*t10*t11*t6*t81*t9+2.390625d-2*t11*t27
     +          *t37*t6*t9+t17*t82+4.78125d-2*gammabb*t11*t27*t6*s2f
     +          +t74+t73-4.78125d-2*gammabb*t10*t11*t28*t6
     +          -1.6562500000000002d-2*gammabb*s3f*t16*t58-1.325d-1
     +          *gammabb*t13*t16*t27*t40+1.325d-1*t10*t13*t16*t29
     +          -1.325d-1*s3f*t13*t16*t27
              t84 = s1c2(D2_GAA_GAA)
              t85 = -2.390625d-2*t4*t84-9.5625d-2*gammaaa*t31
     +          -4.78125d-2*s1f
              t86 = 1.6562500000000002d-2*s3f*t16*t22
              t87 = 2.65d-1*t13*t16*t34*t5
              t88 = s3c2(D2_GAA_GAA)
              t89 = t17*t88
              t90 = s3c2(D2_GAA_GAB)
              t91 = t17*t90+1.325d-1*t13*t16*t36*t5
              t92 = s3c2(D2_GAA_GBB)
              t93 = t17*t92+1.325d-1*t13*t16*t40*t5
     +          +1.6562500000000002d-2*t10*s3f*t16*t5+1.325d-1*t10*t13
     +          *t16*t34
              t94 = s3c2(D2_GAB_GAB)
              t95 = s3c2(D2_GAB_GBB)
              t96 = s2c2(D2_GBB_GBB)
              t97 = -2.390625d-2*t9*t96-4.78125d-2*s2f-9.5625d-2
     +          *gammabb*t37
              t98 = 1.6562500000000002d-2*s3f*t16*t27
              t99 = 2.65d-1*t10*t13*t16*t40
              t100 = s3c2(D2_GBB_GBB)
              t101 = t100*t17
              t102 = 1.0d+0*t45*wght
              t103 = -4.78125d-2*t20*s1f*t4*t42*t5
              t104 = -1.0600000000000001d+0*t14*t25*t43
              t105 = 2.65d-1*gammaaa*s3f*t13*t22*t43
              t106 = 4.78125d-2*t23*t4*t42*t5*t6-4.78125d-2*t22*s1f*t4
     +          *t42*t6+t105+t104+t103
              t107 = 1.0d+0*t47*wght
              t108 = -4.78125d-2*t10*t20*t46*s2f*t9
              t109 = t108+t105+t104
              t110 = -1.0600000000000001d+0*t14*t29*t43
              t111 = 2.65d-1*gammabb*s3f*t13*t27*t43
              t112 = t111+t110+t103
              t113 = -4.78125d-2*t27*t46*t6*s2f*t9+4.78125d-2*t10*t28
     +          *t46*t6*t9+t111+t110+t108
              t114 = -1.0600000000000001d+0*t14*t34*t43
              t115 = -2.65d-1*s3f*t13*t43*t5
              t116 = t115+t114
              t117 = -1.0600000000000001d+0*t1*t14*t36*t43*wght
              t118 = -1.0600000000000001d+0*t14*t40*t43
              t119 = -2.65d-1*t10*s3f*t13*t43
              t120 = t119+t118
              t121 = 1/taua**4
              t122 = 1/t15**4
              t123 = 3.18d+0*s3f*t122*t14
              t124 = 1/taub**4
              t125 = 1/t1**4
              t126 = 1.434375d-1*t125*s1f*t4*t5*t7
              t127 = 1.434375d-1*t10*t11*t125*s2f*t9
              t128 = 2.0d+0*t57*wght
              t129 = -9.5625d-2*gammaaa*s1f*t48*t5*t7
              t130 = -4.78125d-2*t31*t4*t48*t5*t7
              t131 = -9.5625d-2*gammabb*t10*t11*t48*s2f
              t132 = -4.78125d-2*t10*t11*t37*t48*t9
              t133 = 1.0d+0*(t89+t87+t86-2.390625d-2*t4*t5*t6*t7*t84
     +          -9.5625d-2*gammaaa*t31*t5*t6*t7-4.78125d-2*s1f*t5*t6
     +          *t7)*wght
              t134 = 4.78125d-2*t20*s1f*t5*t7
              t135 = 9.5625d-2*gammaaa*t20*t31*t5*t7
              t136 = 2.390625d-2*t20*t4*t5*t7*t84
              t137 = 1.0d+0*t91*wght
              t138 = 1.0d+0*t93*wght
              t139 = 1.0d+0*t17*t94*wght
              t140 = 1.325d-1*t10*t13*t16*t36*wght
              t141 = 1.0d+0*t17*t95*wght
              t142 = 1.0d+0*(t99+t98-2.390625d-2*t10*t11*t6*t9*t96
     +          -4.78125d-2*t10*t11*t6*s2f-9.5625d-2*gammabb*t10*t11
     +          *t37*t6+t101)*wght
              t143 = 4.78125d-2*t10*t11*t20*s2f
              t144 = 9.5625d-2*gammabb*t10*t11*t20*t37
              t145 = 2.390625d-2*t10*t11*t20*t9*t96
              t146 = 9.5625d-2*s1f*t4*t42*t48*t5
              t147 = -3.3125000000000004d-2*s3f*t4*t43*t53
              t148 = -1.0600000000000001d+0*t14*t43*t54
              t149 = -5.3d-1*gammaaa*s3f*t13*t43*t50
              t150 = 5.3d-1*gammaaa*t13*t22*t25*t43
              t151 = 9.5625d-2*t10*t46*t48*s2f*t9
              t152 = -1.0600000000000001d+0*t14*t43*t56
              t153 = -3.3125000000000004d-2*gammaaa*gammabb*s3f*t22
     +          *t27*t43
              t154 = 2.65d-1*gammabb*t13*t25*t27*t43
              t155 = 2.65d-1*gammaaa*t13*t22*t29*t43
              t156 = -1.0600000000000001d+0*t14*t43*t61
              t157 = -3.3125000000000004d-2*s3f*t43*t60*t9
              t158 = -5.3d-1*gammabb*s3f*t13*t43*t58
              t159 = 5.3d-1*gammabb*t13*t27*t29*t43
              t160 = 1.0d+0*(4.78125d-2*t31*t4*t42*t5*t6+9.5625d-2
     +          *gammaaa*s1f*t42*t5*t6+t115+t114)*wght
              t161 = -9.5625d-2*gammaaa*t20*s1f*t42*t5
              t162 = -4.78125d-2*t20*t31*t4*t42*t5
              t163 = 3.3125000000000004d-2*gammaaa*s3f*t43*t50
              t164 = -1.0600000000000001d+0*t14*t43*t67
              t165 = 2.65d-1*s3f*t13*t22*t43
              t166 = 2.65d-1*gammaaa*t13*t22*t34*t43
              t167 = -2.65d-1*t13*t25*t43*t5
              t168 = 1.0d+0*t116*wght
              t169 = -1.0600000000000001d+0*t14*t36*t43*wght
              t170 = 1.0d+0*t1*(2.65d-1*gammaaa*t13*t22*t36*t43
     +          -1.0600000000000001d+0*t14*t43*t70)*wght+t169
              t171 = 1.0d+0*t120*wght
              t172 = -1.0600000000000001d+0*t14*t43*t75
              t173 = 3.3125000000000004d-2*gammaaa*t10*s3f*t22*t43
              t174 = -2.65d-1*t10*t13*t25*t43
              t175 = 2.65d-1*gammaaa*t13*t22*t40*t43
              t176 = 1.0d+0*(4.78125d-2*t10*t37*t46*t6*t9+9.5625d-2
     +          *gammabb*t10*t46*t6*s2f+t119+t118)*wght
              t177 = -9.5625d-2*gammabb*t10*t20*t46*s2f
              t178 = -4.78125d-2*t10*t20*t37*t46*t9
              t179 = -1.0600000000000001d+0*t14*t43*t77
              t180 = 3.3125000000000004d-2*gammabb*s3f*t27*t43*t5
              t181 = 2.65d-1*gammabb*t13*t27*t34*t43
              t182 = -2.65d-1*t13*t29*t43*t5
              t183 = 1.0d+0*t1*(2.65d-1*gammabb*t13*t27*t36*t43
     +          -1.0600000000000001d+0*t14*t43*t79)*wght+t169
              t184 = -1.0600000000000001d+0*t14*t43*t82
              t185 = 3.3125000000000004d-2*gammabb*s3f*t43*t58
              t186 = 2.65d-1*s3f*t13*t27*t43
              t187 = 2.65d-1*gammabb*t13*t27*t40*t43
              t188 = -2.65d-1*t10*t13*t29*t43
              t189 = -3.3125000000000004d-2*s3f*t22*t43
              t190 = -1.0600000000000001d+0*t14*t43*t88
              t191 = -5.3d-1*t13*t34*t43*t5
              t192 = 1.0d+0*t1*(-1.0600000000000001d+0*t14*t43*t90
     +          -2.65d-1*t13*t36*t43*t5)*wght
              t193 = 1.0d+0*t1*(-1.0600000000000001d+0*t14*t43*t92
     +          -2.65d-1*t13*t40*t43*t5-3.3125000000000004d-2*t10*s3f
     +          *t43*t5-2.65d-1*t10*t13*t34*t43)*wght
              t194 = -1.0600000000000001d+0*t1*t14*t43*t94*wght
              t195 = -1.0600000000000001d+0*t1*t14*t43*t95*wght
     +          -2.65d-1*t1*t10*t13*t36*t43*wght
              t196 = -1.0600000000000001d+0*t100*t14*t43
              t197 = -3.3125000000000004d-2*s3f*t27*t43
              t198 = -5.3d-1*t10*t13*t40*t43
              t199 = 3.18d+0*t122*t14*t34
              t200 = 7.95d-1*s3f*t122*t13*t5
              t201 = 1.0d+0*t1*(t200+t199)*wght
              t202 = 3.18d+0*t1*t122*t14*t36*wght
              t203 = 3.18d+0*t122*t14*t40
              t204 = 7.95d-1*t10*s3f*t122*t13
              t205 = 1.0d+0*t1*(t204+t203)*wght
              t206 = 1/t15**5
              t207 = -1.272d+1*s3f*t14*t206
              t208 = -1.272d+1*t1*s3f*t14*t206*wght
              fnc(iq) = 1.0d+0*t1*t18*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t26*wght+t19+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t30*wght+t19+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = 1.0d+0*t1*(t32*t5*t6*t7+t35+t33)*wght
     +          +Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = 1.0d+0*t1*t17*t36*wght+Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = 1.0d+0*t1*(t10*t11*t38*t6+t41+t39)
     +          *wght+Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = 1.0d+0*t1*t45*wght+Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = 1.0d+0*t1*t47*wght+Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*t55*wght+2.0d+0*t26*wght
     +          +Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*t57*wght+1.0d+0*t30*wght
     +          +1.0d+0*t26*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*t62*wght+2.0d+0*t30*wght
     +          +Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = 1.0d+0*t1*t68*wght+t63+Cmat2(iq,
     +          D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = 1.0d+0*t1*t71*wght+t69+Cmat2(iq,
     +          D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = 1.0d+0*t1*t76*wght+t72+Cmat2(iq,
     +          D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = 1.0d+0*t1*t78*wght+t63+Cmat2(iq,
     +          D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = 1.0d+0*t1*t80*wght+t69+Cmat2(iq,
     +          D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = 1.0d+0*t1*t83*wght+t72+Cmat2(iq,
     +          D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t1*(t89+t87+t86+t5*t6*t7
     +          *t85)*wght+Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = 1.0d+0*t1*t91*wght+Cmat2(iq,
     +          D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = 1.0d+0*t1*t93*wght+Cmat2(iq,
     +          D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = 1.0d+0*t1*t17*t94*wght+Cmat2(iq,
     +          D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = 1.0d+0*t1*t17*t95*wght+1.325d-1
     +          *t1*t10*t13*t16*t36*wght+Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t1*(t99+t98+t10*t11*t6*t97
     +          +t101)*wght+Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = 1.0d+0*t1*t106*wght+t102+Mmat2(iq,
     +          D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = 1.0d+0*t1*t109*wght+t107+Mmat2(iq,
     +          D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = 1.0d+0*t1*t112*wght+t102+Mmat2(iq,
     +          D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = 1.0d+0*t1*t113*wght+t107+Mmat2(iq,
     +          D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = 1.0d+0*t1*(-2*t32*t42*t5*t6+t115
     +          +t114)*wght+Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = 1.0d+0*t1*t116*wght+Mmat2(iq,
     +          D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = t117+Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = t117+Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = 1.0d+0*t1*t120*wght+Mmat2(iq,
     +          D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = 1.0d+0*t1*(-2*t10*t38*t46*t6+t119
     +          +t118)*wght+Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = 1.0d+0*t1*(t123-1.434375d-1*t121
     +          *s1f*t4*t5*t6)*wght+Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = 3.18d+0*t1*s3f*t122*t14*wght+Mmat2
     +          (iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = 1.0d+0*t1*(t123-1.434375d-1*t10
     +          *t124*t6*s2f*t9)*wght+Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = 1.0d+0*t1*(1.434375d-1*s1f*t4
     +          *t53*t6*t7+7.171875d-2*t22*t4*t51*t6*t7-1.434375d-1
     +          *t23*t4*t50*t6*t7-2.390625d-2*s1a3(D3_RA_RA_RA)*t4*t5
     +          *t6*t7+7.171875d-2*t20*t4*t5*t51*t7+1.434375d-1*t20
     +          *s1f*t4*t50*t7-1.434375d-1*t23*t4*t48*t5*t7
     +          +1.434375d-1*t22*s1f*t4*t48*t7-1.434375d-1*t20*t22*t23
     +          *t4*t7-3.975d-1*gammaaa*t13*t16*t22*t54+4.96875d-2*t16
     +          *t25*t4*t53-7.95d-1*gammaaa*s3f*t13*t16*t53+7.95d-1
     +          *gammaaa*t13*t16*t25*t50-9.9375d-2*s3f*t16*t4/rhoa**5
     +          +s3a3(D3_RA_RA_RA)*t17+t127+t126)*wght+3.0d+0*t55*wght
     +          +Amat3(iq,D3_RA_RA_RA)
              Amat3(iq,D3_RA_RA_RB) = 1.0d+0*t1*(4.78125d-2*t11*t27
     +          *t48*s2f*t9-4.78125d-2*t10*t11*t28*t48*t9+2.390625d-2
     +          *t20*t4*t5*t51*t7+4.78125d-2*t20*s1f*t4*t50*t7
     +          -9.5625d-2*t23*t4*t48*t5*t7+9.5625d-2*t22*s1f*t4*t48
     +          *t7-4.78125d-2*t20*t22*t23*t4*t7-2.65d-1*gammaaa*t13
     +          *t16*t22*t56-1.325d-1*gammabb*t13*t16*t27*t54
     +          +1.6562500000000002d-2*t16*t29*t4*t53+2.65d-1*gammaaa
     +          *t13*t16*t29*t50-3.3125000000000004d-2*gammaaa*gammabb
     +          *s3f*t16*t27*t50+3.3125000000000004d-2*gammaaa*gammabb
     +          *t16*t22*t25*t27+s3a3(D3_RA_RA_RB)*t17+t127+t126)*wght
     +          +1.0d+0*t55*wght+t128+Amat3(iq,D3_RA_RA_RB)
              Amat3(iq,D3_RA_RB_RB) = 1.0d+0*t1*(4.78125d-2*t11*t20
     +          *t58*s2f*t9+9.5625d-2*t11*t27*t48*s2f*t9
     +          +1.6562500000000002d-2*t16*t25*t60*t9+2.390625d-2*t10
     +          *t11*t20*t59*t9-9.5625d-2*t10*t11*t28*t48*t9
     +          -4.78125d-2*t11*t20*t27*t28*t9-4.78125d-2*t23*t4*t48
     +          *t5*t7+4.78125d-2*t22*s1f*t4*t48*t7-1.325d-1*gammaaa
     +          *t13*t16*t22*t61+2.65d-1*gammabb*t13*t16*t25*t58
     +          -3.3125000000000004d-2*gammaaa*gammabb*s3f*t16*t22*t58
     +          -2.65d-1*gammabb*t13*t16*t27*t56+3.3125000000000004d-2
     +          *gammaaa*gammabb*t16*t22*t27*t29+s3a3(D3_RA_RB_RB)*t17
     +          +t127+t126)*wght+1.0d+0*t62*wght+t128+Amat3(iq,
     +          D3_RA_RB_RB)
              Amat3(iq,D3_RB_RB_RB) = 1.0d+0*t1*(-2.390625d-2*t10*t11
     +          *t6*s2a3(D3_RB_RB_RB)*t9+1.434375d-1*t11*t6*t60*s2f*t9
     +          +1.434375d-1*t11*t20*t58*s2f*t9+1.434375d-1*t11*t27
     +          *t48*s2f*t9+4.96875d-2*t16*t29*t60*t9+7.171875d-2*t11
     +          *t27*t59*t6*t9-1.434375d-1*t11*t28*t58*t6*t9
     +          +7.171875d-2*t10*t11*t20*t59*t9-1.434375d-1*t10*t11
     +          *t28*t48*t9-1.434375d-1*t11*t20*t27*t28*t9-9.9375d-2
     +          *s3f*t16*t9/rhob**5-3.975d-1*gammabb*t13*t16*t27*t61
     +          -7.95d-1*gammabb*s3f*t13*t16*t60+7.95d-1*gammabb*t13
     +          *t16*t29*t58+s3a3(D3_RB_RB_RB)*t17+t127+t126)*wght
     +          +3.0d+0*t62*wght+Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*t1*(4.78125d-2*t22*t4*t6
     +          *t66*t7+4.78125d-2*t20*t4*t5*t66*t7-4.78125d-2*gammaaa
     +          *t5*t51*t6*t7-4.78125d-2*t31*t4*t50*t6*t7-9.5625d-2
     +          *gammaaa*s1f*t50*t6*t7-2.390625d-2*s1c3(D3_RA_RA_GAA)
     +          *t4*t5*t6*t7+9.5625d-2*gammaaa*t22*t23*t6*t7+9.5625d-2
     +          *gammaaa*t20*t23*t5*t7-4.78125d-2*t20*t22*t31*t4*t7
     +          -9.5625d-2*gammaaa*t20*t22*s1f*t7-2.65d-1*gammaaa*t13
     +          *t16*t22*t67+1.325d-1*t13*t16*t5*t54
     +          +1.6562500000000002d-2*t16*t34*t4*t53
     +          +6.625000000000001d-2*gammaaa*s3f*t16*t53+2.65d-1
     +          *gammaaa*t13*t16*t34*t50-3.3125000000000004d-2*gammaaa
     +          *t16*t25*t50+2.65d-1*s3f*t13*t16*t50-2.65d-1*t13*t16
     +          *t22*t25+s3c3(D3_RA_RA_GAA)*t17+t130+t129)*wght+2.0d+0
     +          *t68*wght+Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_RA_GAB) = 2.0d+0*t71*wght+1.0d+0*t1*(
     +          -2.65d-1*gammaaa*t13*t16*t22*t70+1.6562500000000002d-2
     +          *t16*t36*t4*t53+2.65d-1*gammaaa*t13*t16*t36*t50+s3c3
     +          (D3_RA_RA_GAB)*t17)*wght+Cmat3(iq,D3_RA_RA_GAB)
              Cmat3(iq,D3_RA_RA_GBB) = 2.0d+0*t76*wght+1.0d+0*t1*(
     +          -2.65d-1*gammaaa*t13*t16*t22*t75+1.325d-1*t10*t13*t16
     +          *t54+1.6562500000000002d-2*t16*t4*t40*t53+2.65d-1
     +          *gammaaa*t13*t16*t40*t50+3.3125000000000004d-2*gammaaa
     +          *t10*s3f*t16*t50-3.3125000000000004d-2*gammaaa*t10*t16
     +          *t22*t25+s3c3(D3_RA_RA_GBB)*t17+t132+t131)*wght+Cmat3
     +          (iq,D3_RA_RA_GBB)
              Cmat3(iq,D3_RA_RB_GAA) = 1.0d+0*t78*wght+1.0d+0*t1*(
     +          -1.325d-1*gammaaa*t13*t16*t22*t77+2.390625d-2*t20*t4
     +          *t5*t66*t7+4.78125d-2*gammaaa*t20*t23*t5*t7
     +          -2.390625d-2*t20*t22*t31*t4*t7-4.78125d-2*gammaaa*t20
     +          *t22*s1f*t7-1.325d-1*gammabb*t13*t16*t27*t67+1.325d-1
     +          *t13*t16*t5*t56-1.6562500000000002d-2*gammaaa*t16*t29
     +          *t50-1.6562500000000002d-2*gammabb*t16*t25*t27*t5
     +          +1.6562500000000002d-2*gammaaa*gammabb*t16*t22*t27*t34
     +          -1.325d-1*t13*t16*t22*t29+1.6562500000000002d-2
     +          *gammabb*s3f*t16*t22*t27+s3c3(D3_RA_RB_GAA)*t17+t130
     +          +t129)*wght+1.0d+0*t68*wght+Cmat3(iq,D3_RA_RB_GAA)
              Cmat3(iq,D3_RA_RB_GAB) = 1.0d+0*t80*wght+1.0d+0*t1*(
     +          -1.325d-1*gammaaa*t13*t16*t22*t79-1.325d-1*gammabb*t13
     +          *t16*t27*t70+1.6562500000000002d-2*gammaaa*gammabb*t16
     +          *t22*t27*t36+s3c3(D3_RA_RB_GAB)*t17)*wght+1.0d+0*t71
     +          *wght+Cmat3(iq,D3_RA_RB_GAB)
              Cmat3(iq,D3_RA_RB_GBB) = 1.0d+0*t1*(2.390625d-2*t10*t11
     +          *t20*t81*t9-2.390625d-2*t11*t20*t27*t37*t9-1.325d-1
     +          *gammaaa*t13*t16*t22*t82-4.78125d-2*gammabb*t11*t20
     +          *t27*s2f-1.325d-1*gammabb*t13*t16*t27*t75
     +          -1.6562500000000002d-2*gammabb*t16*t25*t58+1.325d-1
     +          *t10*t13*t16*t56+1.6562500000000002d-2*gammaaa*gammabb
     +          *t16*t22*t27*t40-1.6562500000000002d-2*gammaaa*t10*t16
     +          *t22*t29+4.78125d-2*gammabb*t10*t11*t20*t28-1.325d-1
     +          *t13*t16*t25*t27+1.6562500000000002d-2*gammaaa*s3f*t16
     +          *t22*t27+s3c3(D3_RA_RB_GBB)*t17+t132+t131)*wght+1.0d+0
     +          *t83*wght+1.0d+0*t76*wght+Cmat3(iq,D3_RA_RB_GBB)
              Cmat3(iq,D3_RB_RB_GAA) = 1.0d+0*t1*
     +          (1.6562500000000002d-2*t16*t34*t60*t9-2.65d-1*gammabb
     +          *t13*t16*t27*t77+1.325d-1*t13*t16*t5*t61
     +          +3.3125000000000004d-2*gammabb*s3f*t16*t5*t58+2.65d-1
     +          *gammabb*t13*t16*t34*t58-3.3125000000000004d-2*gammabb
     +          *t16*t27*t29*t5+s3c3(D3_RB_RB_GAA)*t17+t130+t129)*wght
     +          +2.0d+0*t78*wght+Cmat3(iq,D3_RB_RB_GAA)
              Cmat3(iq,D3_RB_RB_GAB) = 1.0d+0*t1*
     +          (1.6562500000000002d-2*t16*t36*t60*t9-2.65d-1*gammabb
     +          *t13*t16*t27*t79+2.65d-1*gammabb*t13*t16*t36*t58+s3c3
     +          (D3_RB_RB_GAB)*t17)*wght+2.0d+0*t80*wght+Cmat3(iq,
     +          D3_RB_RB_GAB)
              Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*t1*(4.78125d-2*t11*t27
     +          *t6*t81*t9+4.78125d-2*t10*t11*t20*t81*t9-2.390625d-2
     +          *t10*t11*t6*s2c3(D3_RB_RB_GBB)*t9
     +          +1.6562500000000002d-2*t16*t40*t60*t9-4.78125d-2*t11
     +          *t37*t58*t6*t9-4.78125d-2*t11*t20*t27*t37*t9-2.65d-1
     +          *gammabb*t13*t16*t27*t82-9.5625d-2*gammabb*t11*t58*t6
     +          *s2f-9.5625d-2*gammabb*t11*t20*t27*s2f+1.325d-1*t10
     +          *t13*t16*t61+6.625000000000001d-2*gammabb*s3f*t16*t60
     +          -4.78125d-2*gammabb*t10*t11*t59*t6+9.5625d-2*gammabb
     +          *t11*t27*t28*t6+2.65d-1*gammabb*t13*t16*t40*t58
     +          -3.3125000000000004d-2*gammabb*t16*t29*t58+2.65d-1*s3f
     +          *t13*t16*t58-2.65d-1*t13*t16*t27*t29+9.5625d-2*gammabb
     +          *t10*t11*t20*t28+s3c3(D3_RB_RB_GBB)*t17+t132+t131)
     +          *wght+2.0d+0*t83*wght+Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*t1*(-1.325d-1*gammaaa
     +          *t13*t16*t22*t88+2.390625d-2*t22*t4*t6*t7*t84
     +          -9.5625d-2*gammaaa*t5*t6*t66*t7-2.390625d-2*s1c3
     +          (D3_RA_GAA_GAA)*t4*t5*t6*t7-4.78125d-2*t23*t5*t6*t7
     +          +9.5625d-2*gammaaa*t22*t31*t6*t7+4.78125d-2*t22*s1f*t6
     +          *t7+2.65d-1*t13*t16*t5*t67-3.3125000000000004d-2
     +          *gammaaa*t16*t34*t50-3.3125000000000004d-2*s3f*t16*t50
     +          -2.65d-1*t13*t16*t22*t34+1.6562500000000002d-2*t16*t22
     +          *t25+s3c3(D3_RA_GAA_GAA)*t17+t136+t135+t134)*wght+t133
     +          +Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_RA_GAA_GAB) = 1.0d+0*t1*(-1.325d-1*gammaaa
     +          *t13*t16*t22*t90+1.325d-1*t13*t16*t5*t70
     +          -1.6562500000000002d-2*gammaaa*t16*t36*t50-1.325d-1
     +          *t13*t16*t22*t36+s3c3(D3_RA_GAA_GAB)*t17)*wght+t137
     +          +Cmat3(iq,D3_RA_GAA_GAB)
              Cmat3(iq,D3_RA_GAA_GBB) = 1.0d+0*t1*(-1.325d-1*gammaaa
     +          *t13*t16*t22*t92+1.325d-1*t13*t16*t5*t75+1.325d-1*t10
     +          *t13*t16*t67-1.6562500000000002d-2*gammaaa*t16*t40*t50
     +          +1.6562500000000002d-2*t10*t16*t25*t5-1.325d-1*t13*t16
     +          *t22*t40-1.6562500000000002d-2*gammaaa*t10*t16*t22*t34
     +          -1.6562500000000002d-2*t10*s3f*t16*t22+s3c3
     +          (D3_RA_GAA_GBB)*t17)*wght+t138+Cmat3(iq,D3_RA_GAA_GBB)
              Cmat3(iq,D3_RA_GAB_GAB) = 1.0d+0*t1*(s3c3(D3_RA_GAB_GAB)
     +          *t17-1.325d-1*gammaaa*t13*t16*t22*t94)*wght+t139+Cmat3
     +          (iq,D3_RA_GAB_GAB)
              Cmat3(iq,D3_RA_GAB_GBB) = 1.0d+0*t1*(-1.325d-1*gammaaa
     +          *t13*t16*t22*t95+1.325d-1*t10*t13*t16*t70
     +          -1.6562500000000002d-2*gammaaa*t10*t16*t22*t36+s3c3
     +          (D3_RA_GAB_GBB)*t17)*wght+t141+t140+Cmat3(iq,
     +          D3_RA_GAB_GBB)
              Cmat3(iq,D3_RA_GBB_GBB) = 1.0d+0*t1*(2.65d-1*t10*t13*t16
     +          *t75-3.3125000000000004d-2*gammaaa*t10*t16*t22*t40
     +          +1.6562500000000002d-2*t16*t25*t27-1.325d-1*gammaaa
     +          *t100*t13*t16*t22+s3c3(D3_RA_GBB_GBB)*t17+t145+t144
     +          +t143)*wght+t142+Cmat3(iq,D3_RA_GBB_GBB)
              Cmat3(iq,D3_RB_GAA_GAA) = 1.0d+0*t1*(-1.325d-1*gammabb
     +          *t13*t16*t27*t88+2.65d-1*t13*t16*t5*t77
     +          -3.3125000000000004d-2*gammabb*t16*t27*t34*t5
     +          +1.6562500000000002d-2*t16*t22*t29+s3c3(D3_RB_GAA_GAA)
     +          *t17+t136+t135+t134)*wght+t133+Cmat3(iq,D3_RB_GAA_GAA)
              Cmat3(iq,D3_RB_GAA_GAB) = 1.0d+0*t1*(-1.325d-1*gammabb
     +          *t13*t16*t27*t90+1.325d-1*t13*t16*t5*t79
     +          -1.6562500000000002d-2*gammabb*t16*t27*t36*t5+s3c3
     +          (D3_RB_GAA_GAB)*t17)*wght+t137+Cmat3(iq,D3_RB_GAA_GAB)
              Cmat3(iq,D3_RB_GAA_GBB) = 1.0d+0*t1*(-1.325d-1*gammabb
     +          *t13*t16*t27*t92+1.325d-1*t13*t16*t5*t82+1.325d-1*t10
     +          *t13*t16*t77-1.6562500000000002d-2*gammabb*t16*t34*t58
     +          -1.6562500000000002d-2*gammabb*t16*t27*t40*t5
     +          +1.6562500000000002d-2*t10*t16*t29*t5
     +          -1.6562500000000002d-2*s3f*t16*t27*t5-1.325d-1*t13*t16
     +          *t27*t34+s3c3(D3_RB_GAA_GBB)*t17)*wght+t138+Cmat3(iq,
     +          D3_RB_GAA_GBB)
              Cmat3(iq,D3_RB_GAB_GAB) = 1.0d+0*t1*(s3c3(D3_RB_GAB_GAB)
     +          *t17-1.325d-1*gammabb*t13*t16*t27*t94)*wght+t139+Cmat3
     +          (iq,D3_RB_GAB_GAB)
              Cmat3(iq,D3_RB_GAB_GBB) = 1.0d+0*t1*(-1.325d-1*gammabb
     +          *t13*t16*t27*t95+1.325d-1*t10*t13*t16*t79
     +          -1.6562500000000002d-2*gammabb*t16*t36*t58-1.325d-1
     +          *t13*t16*t27*t36+s3c3(D3_RB_GAB_GBB)*t17)*wght+t141
     +          +t140+Cmat3(iq,D3_RB_GAB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*t1*(2.390625d-2*t11*t27
     +          *t6*t9*t96-2.390625d-2*t10*t11*t6*s2c3(D3_RB_GBB_GBB)
     +          *t9+2.65d-1*t10*t13*t16*t82-9.5625d-2*gammabb*t10*t11
     +          *t6*t81+4.78125d-2*t11*t27*t6*s2f+9.5625d-2*gammabb
     +          *t11*t27*t37*t6-4.78125d-2*t10*t11*t28*t6
     +          -3.3125000000000004d-2*gammabb*t16*t40*t58
     +          -3.3125000000000004d-2*s3f*t16*t58-2.65d-1*t13*t16*t27
     +          *t40+1.6562500000000002d-2*t16*t27*t29-1.325d-1
     +          *gammabb*t100*t13*t16*t27+s3c3(D3_RB_GBB_GBB)*t17+t145
     +          +t144+t143)*wght+t142+Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GAA_GAA_GAA) = 1.0d+0*t1*(3.975d-1*t13*t16
     +          *t5*t88+t5*t6*t7*(-1.434375d-1*gammaaa*t84-2.390625d-2
     +          *s1c3(D3_GAA_GAA_GAA)*t4-1.434375d-1*t31)+4.96875d-2
     +          *t16*t22*t34+s3c3(D3_GAA_GAA_GAA)*t17)*wght+Cmat3(iq,
     +          D3_GAA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAB) = 1.0d+0*t1*(2.65d-1*t13*t16*t5
     +          *t90+1.6562500000000002d-2*t16*t22*t36+s3c3
     +          (D3_GAA_GAA_GAB)*t17)*wght+Cmat3(iq,D3_GAA_GAA_GAB)
              Cmat3(iq,D3_GAA_GAA_GBB) = 1.0d+0*t1*(2.65d-1*t13*t16*t5
     +          *t92+1.325d-1*t10*t13*t16*t88+3.3125000000000004d-2
     +          *t10*t16*t34*t5+1.6562500000000002d-2*t16*t22*t40+s3c3
     +          (D3_GAA_GAA_GBB)*t17)*wght+Cmat3(iq,D3_GAA_GAA_GBB)
              Cmat3(iq,D3_GAA_GAB_GAB) = 1.0d+0*t1*(1.325d-1*t13*t16
     +          *t5*t94+s3c3(D3_GAA_GAB_GAB)*t17)*wght+Cmat3(iq,
     +          D3_GAA_GAB_GAB)
              Cmat3(iq,D3_GAA_GAB_GBB) = 1.0d+0*t1*(1.325d-1*t13*t16
     +          *t5*t95+1.325d-1*t10*t13*t16*t90+1.6562500000000002d-2
     +          *t10*t16*t36*t5+s3c3(D3_GAA_GAB_GBB)*t17)*wght+Cmat3
     +          (iq,D3_GAA_GAB_GBB)
              Cmat3(iq,D3_GAA_GBB_GBB) = 1.0d+0*t1*(2.65d-1*t10*t13
     +          *t16*t92+3.3125000000000004d-2*t10*t16*t40*t5+1.325d-1
     +          *t100*t13*t16*t5+1.6562500000000002d-2*t16*t27*t34
     +          +s3c3(D3_GAA_GBB_GBB)*t17)*wght+Cmat3(iq,D3_GAA_GBB_GBB)
              Cmat3(iq,D3_GAB_GAB_GAB) = 1.0d+0*t1*s3c3(D3_GAB_GAB_GAB
     +          )*t17*wght+Cmat3(iq,D3_GAB_GAB_GAB)
              Cmat3(iq,D3_GAB_GAB_GBB) = 1.325d-1*t1*t10*t13*t16*t94
     +          *wght+1.0d+0*t1*s3c3(D3_GAB_GAB_GBB)*t17*wght+Cmat3
     +          (iq,D3_GAB_GAB_GBB)
              Cmat3(iq,D3_GAB_GBB_GBB) = 2.65d-1*t1*t10*t13*t16*t95
     +          *wght+1.6562500000000002d-2*t1*t16*t27*t36*wght+1.0d+0
     +          *t1*s3c3(D3_GAB_GBB_GBB)*t17*wght+Cmat3(iq,
     +          D3_GAB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = 1.0d+0*t1*(t10*t11*t6*(
     +          -1.434375d-1*gammabb*t96-2.390625d-2*s2c3
     +          (D3_GBB_GBB_GBB)*t9-1.434375d-1*t37)+4.96875d-2*t16
     +          *t27*t40+s3c3(D3_GBB_GBB_GBB)*t17+3.975d-1*t10*t100
     +          *t13*t16)*wght+Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RA_RA_TA) = 1.0d+0*t1*(4.78125d-2*t4*t42*t5
     +          *t51*t6+9.5625d-2*s1f*t4*t42*t50*t6-9.5625d-2*t22*t23
     +          *t4*t42*t6-9.5625d-2*t20*t23*t4*t42*t5+9.5625d-2*t20
     +          *t22*s1f*t4*t42+t150+t149+t148+t147+t146)*wght+2.0d+0
     +          *t106*wght+Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = 1.0d+0*t1*(t151+t150+t149+t148
     +          +t147)*wght+2.0d+0*t109*wght+Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_RB_TA) = 1.0d+0*t1*(-4.78125d-2*t20*t23
     +          *t4*t42*t5+4.78125d-2*t20*t22*s1f*t4*t42+t155+t154
     +          +t153+t152+t146)*wght+1.0d+0*t112*wght+1.0d+0*t106
     +          *wght+Mmat3(iq,D3_RA_RB_TA)
              Mmat3(iq,D3_RA_RB_TB) = 1.0d+0*t1*(4.78125d-2*t20*t27
     +          *t46*s2f*t9-4.78125d-2*t10*t20*t28*t46*t9+t155+t154
     +          +t153+t152+t151)*wght+1.0d+0*t113*wght+1.0d+0*t109
     +          *wght+Mmat3(iq,D3_RA_RB_TB)
              Mmat3(iq,D3_RB_RB_TA) = 1.0d+0*t1*(t159+t158+t157+t156
     +          +t146)*wght+2.0d+0*t112*wght+Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = 1.0d+0*t1*(9.5625d-2*t46*t58*t6
     +          *s2f*t9+9.5625d-2*t20*t27*t46*s2f*t9+4.78125d-2*t10
     +          *t46*t59*t6*t9-9.5625d-2*t27*t28*t46*t6*t9-9.5625d-2
     +          *t10*t20*t28*t46*t9+t159+t158+t157+t156+t151)*wght
     +          +2.0d+0*t113*wght+Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RA_GAA_TA) = 1.0d+0*t1*(4.78125d-2*t4*t42*t5
     +          *t6*t66+9.5625d-2*gammaaa*t23*t42*t5*t6-4.78125d-2*t22
     +          *t31*t4*t42*t6-9.5625d-2*gammaaa*t22*s1f*t42*t6+t167
     +          +t166+t165+t164+t163+t162+t161)*wght+t160+Mmat3(iq,
     +          D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = 1.0d+0*t1*(t167+t166+t165+t164
     +          +t163)*wght+t168+Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_RA_GAB_TA) = t170+Mmat3(iq,D3_RA_GAB_TA)
              Mmat3(iq,D3_RA_GAB_TB) = t170+Mmat3(iq,D3_RA_GAB_TB)
              Mmat3(iq,D3_RA_GBB_TA) = 1.0d+0*t1*(t175+t174+t173+t172)
     +          *wght+t171+Mmat3(iq,D3_RA_GBB_TA)
              Mmat3(iq,D3_RA_GBB_TB) = 1.0d+0*t1*(t178+t177+t175+t174
     +          +t173+t172)*wght+t176+Mmat3(iq,D3_RA_GBB_TB)
              Mmat3(iq,D3_RB_GAA_TA) = 1.0d+0*t1*(t182+t181+t180+t179
     +          +t162+t161)*wght+t160+Mmat3(iq,D3_RB_GAA_TA)
              Mmat3(iq,D3_RB_GAA_TB) = 1.0d+0*t1*(t182+t181+t180+t179)
     +          *wght+t168+Mmat3(iq,D3_RB_GAA_TB)
              Mmat3(iq,D3_RB_GAB_TA) = t183+Mmat3(iq,D3_RB_GAB_TA)
              Mmat3(iq,D3_RB_GAB_TB) = t183+Mmat3(iq,D3_RB_GAB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = 1.0d+0*t1*(t188+t187+t186+t185
     +          +t184)*wght+t171+Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = 1.0d+0*t1*(4.78125d-2*t10*t46
     +          *t6*t81*t9-4.78125d-2*t27*t37*t46*t6*t9-9.5625d-2
     +          *gammabb*t27*t46*t6*s2f+9.5625d-2*gammabb*t10*t28*t46
     +          *t6+t188+t187+t186+t185+t184+t178+t177)*wght+t176
     +          +Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_RA_TA_TA) = 1.0d+0*t1*(-2*t42*t5*t6*t85+t191
     +          +t190+t189)*wght+Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = 1.0d+0*t1*(t191+t190+t189)*wght
     +          +Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = t192+Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_RB_TA_TA) = t192+Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = t193+Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = t193+Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = t194+Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = t194+Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_GAA_GAB_TA) = t195+Mmat3(iq,D3_GAA_GAB_TA)
              Mmat3(iq,D3_GAA_GAB_TB) = t195+Mmat3(iq,D3_GAA_GAB_TB)
              Mmat3(iq,D3_GAA_GBB_TA) = 1.0d+0*t1*(t198+t197+t196)
     +          *wght+Mmat3(iq,D3_GAA_GBB_TA)
              Mmat3(iq,D3_GAA_GBB_TB) = 1.0d+0*t1*(-2*t10*t46*t6*t97
     +          +t198+t197+t196)*wght+Mmat3(iq,D3_GAA_GBB_TB)
              Mmat3(iq,D3_GAB_GAB_TA) = 1.0d+0*t1*(6*t121*t32*t5*t6
     +          +t200+t199)*wght+Mmat3(iq,D3_GAB_GAB_TA)
              Mmat3(iq,D3_GAB_GAB_TB) = t201+Mmat3(iq,D3_GAB_GAB_TB)
              Mmat3(iq,D3_GAB_GBB_TA) = t201+Mmat3(iq,D3_GAB_GBB_TA)
              Mmat3(iq,D3_GAB_GBB_TB) = t202+Mmat3(iq,D3_GAB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = t202+Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = t202+Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = t205+Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = t205+Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = 1.0d+0*t1*(6*t10*t124*t38*t6
     +          +t204+t203)*wght+Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_GAB_TA_TA) = 1.0d+0*t1*(5.7375d-1*s1f*t4*t5
     +          *t6/taua**5+t207)*wght+Mmat3(iq,D3_GAB_TA_TA)
              Mmat3(iq,D3_GAB_TA_TB) = t208+Mmat3(iq,D3_GAB_TA_TB)
              Mmat3(iq,D3_GAB_TB_TB) = t208+Mmat3(iq,D3_GAB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = 1.0d+0*t1*(5.7375d-1*t10*t6*s2f
     +          *t9/taub**5+t207)*wght+Mmat3(iq,D3_GBB_TA_TA)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhob
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t3
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s1c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1a3,s1c,s1c2,s1c3)
              t4 = s1f
              t5 = 1/t1
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = gammaaa
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s2a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s2c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2a3,s2c,s2c2,s2c3)
              t6 = s2f
              t7 = gammaaa**2
              t8 = 1/rhoa
              t9 = 1/taua**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = gammaaa
              sg(G_AB) = 0
              sg(G_BB) = t3
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s3a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s3c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3a3,s3c,s3c2,s3c3)
              t10 = s3f
              t11 = 1.25d-1*gammaaa*t8+5.0d-51
              t12 = t11**2
              t13 = taua+1.0d-50
              t14 = 1/t13**2
              t15 = 5.3d-1*t12*t14+1.0d+0
              t16 = -2.390625d-2*t5*s2f*t7*t8*t9-3.825d-1*rhob*s1f*t5
     +          +s3f*t15
              t17 = 1.0d+0*t16*wght
              t18 = 1/t1**2
              t19 = 3.825d-1*rhob*t18*s1f
              t20 = 2.390625d-2*t18*s2f*t7*t8*t9
              t21 = 1/rhoa**2
              t22 = s2a(D1_RA)
              t23 = s3a(D1_RA)
              t24 = -2.390625d-2*t22*t5*t7*t8*t9+2.390625d-2*t21*t5
     +          *s2f*t7*t9+t15*t23-1.325d-1*gammaaa*s3f*t11*t14*t21
     +          +t20+t19
              t25 = s1a(D1_RB)
              t26 = s3a(D1_RB)
              t27 = -3.825d-1*s1f*t5-3.825d-1*rhob*t25*t5+t15*t26+t20
     +          +t19
              t28 = s2c(D1_GAA)
              t29 = -2.390625d-2*t28*t7-4.78125d-2*gammaaa*s2f
              t30 = 1.325d-1*s3f*t11*t14*t8
              t31 = s3c(D1_GAA)
              t32 = t15*t31
              t33 = 1/taua**3
              t34 = 1/t13**3
              t35 = 4.78125d-2*t33*t5*s2f*t7*t8-1.0600000000000001d+0
     +          *s3f*t12*t34
              t36 = 1/t1**3
              t37 = -7.65d-1*rhob*t36*s1f
              t38 = -4.78125d-2*t36*s2f*t7*t8*t9
              t39 = 1/rhoa**3
              t40 = s2a2(D2_RA_RA)
              t41 = 1/rhoa**4
              t42 = s3a2(D2_RA_RA)
              t43 = -2.390625d-2*t40*t5*t7*t8*t9+4.78125d-2*t18*t22*t7
     +          *t8*t9-4.78125d-2*t39*t5*s2f*t7*t9-4.78125d-2*t18*t21
     +          *s2f*t7*t9+4.78125d-2*t21*t22*t5*t7*t9
     +          +1.6562500000000002d-2*s3f*t14*t41*t7+t15*t42+2.65d-1
     +          *gammaaa*s3f*t11*t14*t39+t38+t37-2.65d-1*gammaaa*t11
     +          *t14*t21*t23
              t44 = s3a2(D2_RA_RB)
              t45 = 2.390625d-2*t18*t22*t7*t8*t9-2.390625d-2*t18*t21
     +          *s2f*t7*t9+t15*t44+3.825d-1*t18*s1f+t38+t37-1.325d-1
     +          *gammaaa*t11*t14*t21*t26+3.825d-1*rhob*t18*t25
              t46 = s1a2(D2_RB_RB)
              t47 = s3a2(D2_RB_RB)
              t48 = -3.825d-1*rhob*t46*t5-7.65d-1*t25*t5+t15*t47
     +          +7.65d-1*t18*s1f+t38+t37+7.65d-1*rhob*t18*t25
              t49 = 1.0d+0*(-2.390625d-2*t28*t5*t7*t8*t9-4.78125d-2
     +          *gammaaa*t5*s2f*t8*t9+t32+t30)*wght
              t50 = 4.78125d-2*gammaaa*t18*s2f*t8*t9
              t51 = 2.390625d-2*t18*t28*t7*t8*t9
              t52 = s2c2(D2_RA_GAA)
              t53 = s3c2(D2_RA_GAA)
              t54 = -2.390625d-2*t5*t52*t7*t8*t9-4.78125d-2*gammaaa
     +          *t22*t5*t8*t9+2.390625d-2*t21*t28*t5*t7*t9+4.78125d-2
     +          *gammaaa*t21*t5*s2f*t9+1.325d-1*t11*t14*t23*t8+t15*t53
     +          +t51+t50-1.6562500000000002d-2*gammaaa*s3f*t14*t39
     +          -1.325d-1*gammaaa*t11*t14*t21*t31-1.325d-1*s3f*t11*t14
     +          *t21
              t55 = s3c2(D2_RB_GAA)
              t56 = 1.325d-1*t11*t14*t26*t8+t15*t55+t51+t50
              t57 = s2c2(D2_GAA_GAA)
              t58 = -2.390625d-2*t57*t7-4.78125d-2*s2f-9.5625d-2
     +          *gammaaa*t28
              t59 = 1.6562500000000002d-2*s3f*t14*t21
              t60 = 2.65d-1*t11*t14*t31*t8
              t61 = s3c2(D2_GAA_GAA)
              t62 = t15*t61
              t63 = 1.0d+0*t35*wght
              t64 = -4.78125d-2*t18*t33*s2f*t7*t8
              t65 = 4.78125d-2*t22*t33*t5*t7*t8-4.78125d-2*t21*t33*t5
     +          *s2f*t7+t64-1.0600000000000001d+0*t12*t23*t34+2.65d-1
     +          *gammaaa*s3f*t11*t21*t34
              t66 = t64-1.0600000000000001d+0*t12*t26*t34
              t67 = -1.0600000000000001d+0*t12*t31*t34
              t68 = -2.65d-1*s3f*t11*t34*t8
              t69 = 1/taua**4
              t70 = 1/t13**4
              t71 = 1/t1**4
              t72 = 2.2949999999999998d+0*rhob*s1f*t71
              t73 = 1.434375d-1*s2f*t7*t71*t8*t9
              t74 = 2.0d+0*t45*wght
              t75 = -9.5625d-2*gammaaa*t36*s2f*t8*t9
              t76 = -4.78125d-2*t28*t36*t7*t8*t9
              t77 = 1.0d+0*(-2.390625d-2*t5*t57*t7*t8*t9-4.78125d-2*t5
     +          *s2f*t8*t9-9.5625d-2*gammaaa*t28*t5*t8*t9+t62+t60+t59)
     +          *wght
              t78 = 4.78125d-2*t18*s2f*t8*t9
              t79 = 9.5625d-2*gammaaa*t18*t28*t8*t9
              t80 = 2.390625d-2*t18*t57*t7*t8*t9
              t81 = 9.5625d-2*t33*t36*s2f*t7*t8
              t82 = 1.0d+0*(4.78125d-2*t28*t33*t5*t7*t8+9.5625d-2
     +          *gammaaa*t33*t5*s2f*t8+t68+t67)*wght
              t83 = -9.5625d-2*gammaaa*t18*t33*s2f*t8
              t84 = -4.78125d-2*t18*t28*t33*t7*t8
              fnc(iq) = 1.0d+0*t1*t16*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t24*wght+t17+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t27*wght+t17+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = 1.0d+0*t1*(t29*t5*t8*t9+t32+t30)*wght
     +          +Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = 1.0d+0*t1*t35*wght+Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*t43*wght+2.0d+0*t24*wght
     +          +Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*t45*wght+1.0d+0*t27*wght
     +          +1.0d+0*t24*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*t48*wght+2.0d+0*t27*wght
     +          +Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = 1.0d+0*t1*t54*wght+t49+Cmat2(iq,
     +          D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = 1.0d+0*t1*t56*wght+t49+Cmat2(iq,
     +          D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t1*(t5*t58*t8*t9+t62+t60
     +          +t59)*wght+Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = 1.0d+0*t1*t65*wght+t63+Mmat2(iq,
     +          D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = 1.0d+0*t1*t66*wght+t63+Mmat2(iq,
     +          D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = 1.0d+0*t1*(-2*t29*t33*t5*t8+t68
     +          +t67)*wght+Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = 1.0d+0*t1*(3.18d+0*s3f*t12*t70
     +          -1.434375d-1*t5*s2f*t69*t7*t8)*wght+Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = 1.0d+0*t1*(-2.390625d-2*t5*s2a3
     +          (D3_RA_RA_RA)*t7*t8*t9+7.171875d-2*t18*t40*t7*t8*t9
     +          -1.434375d-1*t22*t36*t7*t8*t9+1.434375d-1*t41*t5*s2f
     +          *t7*t9+1.434375d-1*t18*t39*s2f*t7*t9+1.434375d-1*t21
     +          *t36*s2f*t7*t9+7.171875d-2*t21*t40*t5*t7*t9
     +          -1.434375d-1*t22*t39*t5*t7*t9-1.434375d-1*t18*t21*t22
     +          *t7*t9+t73+t72+4.96875d-2*t14*t23*t41*t7-9.9375d-2*s3f
     +          *t14*t7/rhoa**5-3.975d-1*gammaaa*t11*t14*t21*t42
     +          -7.95d-1*gammaaa*s3f*t11*t14*t41+7.95d-1*gammaaa*t11
     +          *t14*t23*t39+s3a3(D3_RA_RA_RA)*t15)*wght+3.0d+0*t43
     +          *wght+Amat3(iq,D3_RA_RA_RA)
              Amat3(iq,D3_RA_RA_RB) = 1.0d+0*t1*(2.390625d-2*t18*t40
     +          *t7*t8*t9-9.5625d-2*t22*t36*t7*t8*t9+4.78125d-2*t18
     +          *t39*s2f*t7*t9+9.5625d-2*t21*t36*s2f*t7*t9-4.78125d-2
     +          *t18*t21*t22*t7*t9+t73+t72+1.6562500000000002d-2*t14
     +          *t26*t41*t7-2.65d-1*gammaaa*t11*t14*t21*t44-7.65d-1
     +          *t36*s1f+2.65d-1*gammaaa*t11*t14*t26*t39-7.65d-1*rhob
     +          *t25*t36+s3a3(D3_RA_RA_RB)*t15)*wght+1.0d+0*t43*wght
     +          +t74+Amat3(iq,D3_RA_RA_RB)
              Amat3(iq,D3_RA_RB_RB) = 1.0d+0*t1*(-4.78125d-2*t22*t36
     +          *t7*t8*t9+4.78125d-2*t21*t36*s2f*t7*t9+t73+t72
     +          -1.325d-1*gammaaa*t11*t14*t21*t47+3.825d-1*rhob*t18
     +          *t46-1.53d+0*t36*s1f-1.53d+0*rhob*t25*t36+7.65d-1*t18
     +          *t25+s3a3(D3_RA_RB_RB)*t15)*wght+1.0d+0*t48*wght+t74
     +          +Amat3(iq,D3_RA_RB_RB)
              Amat3(iq,D3_RB_RB_RB) = 1.0d+0*t1*(t73+t72
     +          -1.1474999999999999d+0*t46*t5-3.825d-1*rhob*s1a3
     +          (D3_RB_RB_RB)*t5+1.1474999999999999d+0*rhob*t18*t46
     +          -2.2949999999999998d+0*t36*s1f-2.2949999999999998d+0
     +          *rhob*t25*t36+2.2949999999999998d+0*t18*t25+s3a3
     +          (D3_RB_RB_RB)*t15)*wght+3.0d+0*t48*wght+Amat3(iq,
     +          D3_RB_RB_RB)
              Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*t1*(-2.390625d-2*t5*s2c3
     +          (D3_RA_RA_GAA)*t7*t8*t9+4.78125d-2*t18*t52*t7*t8*t9
     +          -4.78125d-2*gammaaa*t40*t5*t8*t9+9.5625d-2*gammaaa*t18
     +          *t22*t8*t9+4.78125d-2*t21*t5*t52*t7*t9-4.78125d-2*t28
     +          *t39*t5*t7*t9-4.78125d-2*t18*t21*t28*t7*t9-9.5625d-2
     +          *gammaaa*t39*t5*s2f*t9-9.5625d-2*gammaaa*t18*t21*s2f
     +          *t9+9.5625d-2*gammaaa*t21*t22*t5*t9+1.325d-1*t11*t14
     +          *t42*t8+t76+t75+1.6562500000000002d-2*t14*t31*t41*t7
     +          -2.65d-1*gammaaa*t11*t14*t21*t53+6.625000000000001d-2
     +          *gammaaa*s3f*t14*t41+2.65d-1*gammaaa*t11*t14*t31*t39
     +          -3.3125000000000004d-2*gammaaa*t14*t23*t39+2.65d-1*s3f
     +          *t11*t14*t39-2.65d-1*t11*t14*t21*t23+s3c3(D3_RA_RA_GAA
     +          )*t15)*wght+2.0d+0*t54*wght+Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
              Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
              Cmat3(iq,D3_RA_RB_GAA) = 1.0d+0*t1*(2.390625d-2*t18*t52
     +          *t7*t8*t9+4.78125d-2*gammaaa*t18*t22*t8*t9-2.390625d-2
     +          *t18*t21*t28*t7*t9-4.78125d-2*gammaaa*t18*t21*s2f*t9
     +          +1.325d-1*t11*t14*t44*t8+t76+t75-1.325d-1*gammaaa*t11
     +          *t14*t21*t55-1.6562500000000002d-2*gammaaa*t14*t26*t39
     +          -1.325d-1*t11*t14*t21*t26+s3c3(D3_RA_RB_GAA)*t15)*wght
     +          +1.0d+0*t56*wght+1.0d+0*t54*wght+Cmat3(iq,D3_RA_RB_GAA)
              Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
              Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
              Cmat3(iq,D3_RB_RB_GAA) = 1.0d+0*t1*(1.325d-1*t11*t14*t47
     +          *t8+t76+t75+s3c3(D3_RB_RB_GAA)*t15)*wght+2.0d+0*t56
     +          *wght+Cmat3(iq,D3_RB_RB_GAA)
              Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
              Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*t1*(-2.390625d-2*t5
     +          *s2c3(D3_RA_GAA_GAA)*t7*t8*t9-9.5625d-2*gammaaa*t5*t52
     +          *t8*t9-4.78125d-2*t22*t5*t8*t9+2.390625d-2*t21*t5*t57
     +          *t7*t9+4.78125d-2*t21*t5*s2f*t9+9.5625d-2*gammaaa*t21
     +          *t28*t5*t9+t80+2.65d-1*t11*t14*t53*t8+t79+t78-1.325d-1
     +          *gammaaa*t11*t14*t21*t61-3.3125000000000004d-2*gammaaa
     +          *t14*t31*t39-3.3125000000000004d-2*s3f*t14*t39-2.65d-1
     +          *t11*t14*t21*t31+1.6562500000000002d-2*t14*t21*t23
     +          +s3c3(D3_RA_GAA_GAA)*t15)*wght+t77+Cmat3(iq,
     +          D3_RA_GAA_GAA)
              Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
              Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
              Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
              Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
              Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
              Cmat3(iq,D3_RB_GAA_GAA) = 1.0d+0*t1*(t80+2.65d-1*t11*t14
     +          *t55*t8+t79+t78+1.6562500000000002d-2*t14*t21*t26+s3c3
     +          (D3_RB_GAA_GAA)*t15)*wght+t77+Cmat3(iq,D3_RB_GAA_GAA)
              Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
              Cmat3(iq,D3_RB_GAA_GBB) = Cmat3(iq,D3_RB_GAA_GBB)
              Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
              Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GAA_GAA_GAA) = 1.0d+0*t1*(t5*(-2.390625d-2
     +          *s2c3(D3_GAA_GAA_GAA)*t7-1.434375d-1*gammaaa*t57
     +          -1.434375d-1*t28)*t8*t9+3.975d-1*t11*t14*t61*t8
     +          +4.96875d-2*t14*t21*t31+s3c3(D3_GAA_GAA_GAA)*t15)*wght
     +          +Cmat3(iq,D3_GAA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
              Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
              Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
              Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
              Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
              Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
              Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
              Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RA_RA_TA) = 1.0d+0*t1*(t81+4.78125d-2*t33
     +          *t40*t5*t7*t8-9.5625d-2*t18*t22*t33*t7*t8+9.5625d-2
     +          *t33*t39*t5*s2f*t7+9.5625d-2*t18*t21*t33*s2f*t7
     +          -9.5625d-2*t21*t22*t33*t5*t7-3.3125000000000004d-2*s3f
     +          *t34*t41*t7-1.0600000000000001d+0*t12*t34*t42-5.3d-1
     +          *gammaaa*s3f*t11*t34*t39+5.3d-1*gammaaa*t11*t21*t23
     +          *t34)*wght+2.0d+0*t65*wght+Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_RB_TA) = 1.0d+0*t1*(t81-4.78125d-2*t18
     +          *t22*t33*t7*t8+4.78125d-2*t18*t21*t33*s2f*t7
     +          -1.0600000000000001d+0*t12*t34*t44+2.65d-1*gammaaa*t11
     +          *t21*t26*t34)*wght+1.0d+0*t66*wght+1.0d+0*t65*wght
     +          +Mmat3(iq,D3_RA_RB_TA)
              Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
              Mmat3(iq,D3_RB_RB_TA) = 1.0d+0*t1*(t81
     +          -1.0600000000000001d+0*t12*t34*t47)*wght+2.0d+0*t66
     +          *wght+Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RA_GAA_TA) = 1.0d+0*t1*(t84+t83+4.78125d-2
     +          *t33*t5*t52*t7*t8+9.5625d-2*gammaaa*t22*t33*t5*t8
     +          -2.65d-1*t11*t23*t34*t8-4.78125d-2*t21*t28*t33*t5*t7
     +          -9.5625d-2*gammaaa*t21*t33*t5*s2f
     +          -1.0600000000000001d+0*t12*t34*t53
     +          +3.3125000000000004d-2*gammaaa*s3f*t34*t39+2.65d-1
     +          *gammaaa*t11*t21*t31*t34+2.65d-1*s3f*t11*t21*t34)*wght
     +          +t82+Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
              Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
              Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
              Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
              Mmat3(iq,D3_RB_GAA_TA) = 1.0d+0*t1*(t84+t83-2.65d-1*t11
     +          *t26*t34*t8-1.0600000000000001d+0*t12*t34*t55)*wght
     +          +t82+Mmat3(iq,D3_RB_GAA_TA)
              Mmat3(iq,D3_RB_GAA_TB) = Mmat3(iq,D3_RB_GAA_TB)
              Mmat3(iq,D3_RB_GAB_TA) = Mmat3(iq,D3_RB_GAB_TA)
              Mmat3(iq,D3_RB_GAB_TB) = Mmat3(iq,D3_RB_GAB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_RA_TA_TA) = 1.0d+0*t1*(-2*t33*t5*t58*t8
     +          -5.3d-1*t11*t31*t34*t8-1.0600000000000001d+0*t12*t34
     +          *t61-3.3125000000000004d-2*s3f*t21*t34)*wght+Mmat3(iq,
     +          D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
              Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
              Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
              Mmat3(iq,D3_GAA_GBB_TB) = Mmat3(iq,D3_GAA_GBB_TB)
              Mmat3(iq,D3_GAB_GAB_TA) = 1.0d+0*t1*(7.95d-1*s3f*t11*t70
     +          *t8+6*t29*t5*t69*t8+3.18d+0*t12*t31*t70)*wght+Mmat3
     +          (iq,D3_GAB_GAB_TA)
              Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
              Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
              Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_GAB_TA_TA) = 1.0d+0*t1*(5.7375d-1*t5*s2f*t7
     +          *t8/taua**5-1.272d+1*s3f*t12/t13**5)*wght+Mmat3(iq,
     +          D3_GAB_TA_TA)
              Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
              Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhoa
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t3
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s1c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s1f,s1a,
     +          s1a2,s1a3,s1c,s1c2,s1c3)
              t4 = s1f
              t5 = 1/t1
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = gammabb
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s2a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s2c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2a3,s2c,s2c2,s2c3)
              t6 = s2f
              t7 = gammabb**2
              t8 = 1/rhob
              t9 = 1/taub**2
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = t3
              sg(G_AB) = 0
              sg(G_BB) = gammabb
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s3a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s3c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,ipol,1,1.0d0,sr,sg,s3f,
     +          s3a,s3a2,s3a3,s3c,s3c2,s3c3)
              t10 = s3f
              t11 = 1.25d-1*gammabb*t8+5.0d-51
              t12 = t11**2
              t13 = taub+1.0d-50
              t14 = 1/t13**2
              t15 = 5.3d-1*t12*t14+1.0d+0
              t16 = -2.390625d-2*t5*s2f*t7*t8*t9-3.825d-1*rhoa*s1f*t5
     +          +s3f*t15
              t17 = 1.0d+0*t16*wght
              t18 = 1/t1**2
              t19 = 3.825d-1*rhoa*t18*s1f
              t20 = s1a(D1_RA)
              t21 = 2.390625d-2*t18*s2f*t7*t8*t9
              t22 = s3a(D1_RA)
              t23 = -3.825d-1*s1f*t5-3.825d-1*rhoa*t20*t5+t15*t22+t21
     +          +t19
              t24 = 1/rhob**2
              t25 = s2a(D1_RB)
              t26 = s3a(D1_RB)
              t27 = -2.390625d-2*t25*t5*t7*t8*t9+2.390625d-2*t24*t5
     +          *s2f*t7*t9+t15*t26-1.325d-1*gammabb*s3f*t11*t14*t24
     +          +t21+t19
              t28 = s2c(D1_GBB)
              t29 = -2.390625d-2*t28*t7-4.78125d-2*gammabb*s2f
              t30 = 1.325d-1*s3f*t11*t14*t8
              t31 = s3c(D1_GBB)
              t32 = t15*t31
              t33 = 1/taub**3
              t34 = 1/t13**3
              t35 = 4.78125d-2*t33*t5*s2f*t7*t8-1.0600000000000001d+0
     +          *s3f*t12*t34
              t36 = 1/t1**3
              t37 = -7.65d-1*rhoa*t36*s1f
              t38 = s1a2(D2_RA_RA)
              t39 = -4.78125d-2*t36*s2f*t7*t8*t9
              t40 = s3a2(D2_RA_RA)
              t41 = -3.825d-1*rhoa*t38*t5-7.65d-1*t20*t5+t15*t40
     +          +7.65d-1*t18*s1f+t39+t37+7.65d-1*rhoa*t18*t20
              t42 = s3a2(D2_RA_RB)
              t43 = 2.390625d-2*t18*t25*t7*t8*t9-2.390625d-2*t18*t24
     +          *s2f*t7*t9+t15*t42+3.825d-1*t18*s1f+t39+t37-1.325d-1
     +          *gammabb*t11*t14*t22*t24+3.825d-1*rhoa*t18*t20
              t44 = 1/rhob**3
              t45 = s2a2(D2_RB_RB)
              t46 = 1/rhob**4
              t47 = s3a2(D2_RB_RB)
              t48 = -2.390625d-2*t45*t5*t7*t8*t9+4.78125d-2*t18*t25*t7
     +          *t8*t9-4.78125d-2*t44*t5*s2f*t7*t9-4.78125d-2*t18*t24
     +          *s2f*t7*t9+4.78125d-2*t24*t25*t5*t7*t9
     +          +1.6562500000000002d-2*s3f*t14*t46*t7+t15*t47+2.65d-1
     +          *gammabb*s3f*t11*t14*t44+t39+t37-2.65d-1*gammabb*t11
     +          *t14*t24*t26
              t49 = 1.0d+0*(-2.390625d-2*t28*t5*t7*t8*t9-4.78125d-2
     +          *gammabb*t5*s2f*t8*t9+t32+t30)*wght
              t50 = 4.78125d-2*gammabb*t18*s2f*t8*t9
              t51 = 2.390625d-2*t18*t28*t7*t8*t9
              t52 = s3c2(D2_RA_GBB)
              t53 = 1.325d-1*t11*t14*t22*t8+t15*t52+t51+t50
              t54 = s2c2(D2_RB_GBB)
              t55 = s3c2(D2_RB_GBB)
              t56 = -2.390625d-2*t5*t54*t7*t8*t9-4.78125d-2*gammabb
     +          *t25*t5*t8*t9+2.390625d-2*t24*t28*t5*t7*t9+4.78125d-2
     +          *gammabb*t24*t5*s2f*t9+1.325d-1*t11*t14*t26*t8+t15*t55
     +          +t51+t50-1.6562500000000002d-2*gammabb*s3f*t14*t44
     +          -1.325d-1*gammabb*t11*t14*t24*t31-1.325d-1*s3f*t11*t14
     +          *t24
              t57 = s2c2(D2_GBB_GBB)
              t58 = -2.390625d-2*t57*t7-4.78125d-2*s2f-9.5625d-2
     +          *gammabb*t28
              t59 = 1.6562500000000002d-2*s3f*t14*t24
              t60 = 2.65d-1*t11*t14*t31*t8
              t61 = s3c2(D2_GBB_GBB)
              t62 = t15*t61
              t63 = 1.0d+0*t35*wght
              t64 = -4.78125d-2*t18*t33*s2f*t7*t8
              t65 = t64-1.0600000000000001d+0*t12*t22*t34
              t66 = 4.78125d-2*t25*t33*t5*t7*t8-4.78125d-2*t24*t33*t5
     +          *s2f*t7+t64-1.0600000000000001d+0*t12*t26*t34+2.65d-1
     +          *gammabb*s3f*t11*t24*t34
              t67 = -1.0600000000000001d+0*t12*t31*t34
              t68 = -2.65d-1*s3f*t11*t34*t8
              t69 = 1/taub**4
              t70 = 1/t13**4
              t71 = 1/t1**4
              t72 = 2.2949999999999998d+0*rhoa*s1f*t71
              t73 = 1.434375d-1*s2f*t7*t71*t8*t9
              t74 = 2.0d+0*t43*wght
              t75 = -9.5625d-2*gammabb*t36*s2f*t8*t9
              t76 = -4.78125d-2*t28*t36*t7*t8*t9
              t77 = 1.0d+0*(-2.390625d-2*t5*t57*t7*t8*t9-4.78125d-2*t5
     +          *s2f*t8*t9-9.5625d-2*gammabb*t28*t5*t8*t9+t62+t60+t59)
     +          *wght
              t78 = 4.78125d-2*t18*s2f*t8*t9
              t79 = 9.5625d-2*gammabb*t18*t28*t8*t9
              t80 = 2.390625d-2*t18*t57*t7*t8*t9
              t81 = 9.5625d-2*t33*t36*s2f*t7*t8
              t82 = 1.0d+0*(4.78125d-2*t28*t33*t5*t7*t8+9.5625d-2
     +          *gammabb*t33*t5*s2f*t8+t68+t67)*wght
              t83 = -9.5625d-2*gammabb*t18*t33*s2f*t8
              t84 = -4.78125d-2*t18*t28*t33*t7*t8
              fnc(iq) = 1.0d+0*t1*t16*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t23*wght+t17+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t27*wght+t17+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = 1.0d+0*t1*(t29*t5*t8*t9+t32+t30)*wght
     +          +Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = 1.0d+0*t1*t35*wght+Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*t41*wght+2.0d+0*t23*wght
     +          +Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*t43*wght+1.0d+0*t27*wght
     +          +1.0d+0*t23*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*t48*wght+2.0d+0*t27*wght
     +          +Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = 1.0d+0*t1*t53*wght+t49+Cmat2(iq,
     +          D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = 1.0d+0*t1*t56*wght+t49+Cmat2(iq,
     +          D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t1*(t5*t58*t8*t9+t62+t60
     +          +t59)*wght+Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = 1.0d+0*t1*t65*wght+t63+Mmat2(iq,
     +          D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = 1.0d+0*t1*t66*wght+t63+Mmat2(iq,
     +          D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = 1.0d+0*t1*(-2*t29*t33*t5*t8+t68
     +          +t67)*wght+Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = 1.0d+0*t1*(3.18d+0*s3f*t12*t70
     +          -1.434375d-1*t5*s2f*t69*t7*t8)*wght+Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = 1.0d+0*t1*(t73+t72-3.825d-1*rhoa
     +          *s1a3(D3_RA_RA_RA)*t5-1.1474999999999999d+0*t38*t5
     +          -2.2949999999999998d+0*t36*s1f+1.1474999999999999d+0
     +          *rhoa*t18*t38-2.2949999999999998d+0*rhoa*t20*t36
     +          +2.2949999999999998d+0*t18*t20+s3a3(D3_RA_RA_RA)*t15)
     +          *wght+3.0d+0*t41*wght+Amat3(iq,D3_RA_RA_RA)
              Amat3(iq,D3_RA_RA_RB) = 1.0d+0*t1*(-4.78125d-2*t25*t36
     +          *t7*t8*t9+4.78125d-2*t24*t36*s2f*t7*t9+t73+t72
     +          -1.325d-1*gammabb*t11*t14*t24*t40-1.53d+0*t36*s1f
     +          +3.825d-1*rhoa*t18*t38-1.53d+0*rhoa*t20*t36+7.65d-1
     +          *t18*t20+s3a3(D3_RA_RA_RB)*t15)*wght+1.0d+0*t41*wght
     +          +t74+Amat3(iq,D3_RA_RA_RB)
              Amat3(iq,D3_RA_RB_RB) = 1.0d+0*t1*(2.390625d-2*t18*t45
     +          *t7*t8*t9-9.5625d-2*t25*t36*t7*t8*t9+4.78125d-2*t18
     +          *t44*s2f*t7*t9+9.5625d-2*t24*t36*s2f*t7*t9-4.78125d-2
     +          *t18*t24*t25*t7*t9+t73+t72+1.6562500000000002d-2*t14
     +          *t22*t46*t7+2.65d-1*gammabb*t11*t14*t22*t44-2.65d-1
     +          *gammabb*t11*t14*t24*t42-7.65d-1*t36*s1f-7.65d-1*rhoa
     +          *t20*t36+s3a3(D3_RA_RB_RB)*t15)*wght+1.0d+0*t48*wght
     +          +t74+Amat3(iq,D3_RA_RB_RB)
              Amat3(iq,D3_RB_RB_RB) = 1.0d+0*t1*(-2.390625d-2*t5*s2a3
     +          (D3_RB_RB_RB)*t7*t8*t9+7.171875d-2*t18*t45*t7*t8*t9
     +          -1.434375d-1*t25*t36*t7*t8*t9+1.434375d-1*t46*t5*s2f
     +          *t7*t9+1.434375d-1*t18*t44*s2f*t7*t9+1.434375d-1*t24
     +          *t36*s2f*t7*t9+7.171875d-2*t24*t45*t5*t7*t9
     +          -1.434375d-1*t25*t44*t5*t7*t9-1.434375d-1*t18*t24*t25
     +          *t7*t9+t73+t72+4.96875d-2*t14*t26*t46*t7-9.9375d-2*s3f
     +          *t14*t7/rhob**5-3.975d-1*gammabb*t11*t14*t24*t47
     +          -7.95d-1*gammabb*s3f*t11*t14*t46+7.95d-1*gammabb*t11
     +          *t14*t26*t44+s3a3(D3_RB_RB_RB)*t15)*wght+3.0d+0*t48
     +          *wght+Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
              Cmat3(iq,D3_RA_RA_GBB) = 1.0d+0*t1*(1.325d-1*t11*t14*t40
     +          *t8+t76+t75+s3c3(D3_RA_RA_GBB)*t15)*wght+2.0d+0*t53
     +          *wght+Cmat3(iq,D3_RA_RA_GBB)
              Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
              Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
              Cmat3(iq,D3_RA_RB_GBB) = 1.0d+0*t1*(2.390625d-2*t18*t54
     +          *t7*t8*t9+4.78125d-2*gammabb*t18*t25*t8*t9-2.390625d-2
     +          *t18*t24*t28*t7*t9-4.78125d-2*gammabb*t18*t24*s2f*t9
     +          +1.325d-1*t11*t14*t42*t8+t76+t75-1.325d-1*gammabb*t11
     +          *t14*t24*t52-1.6562500000000002d-2*gammabb*t14*t22*t44
     +          -1.325d-1*t11*t14*t22*t24+s3c3(D3_RA_RB_GBB)*t15)*wght
     +          +1.0d+0*t56*wght+1.0d+0*t53*wght+Cmat3(iq,D3_RA_RB_GBB)
              Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
              Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
              Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*t1*(-2.390625d-2*t5*s2c3
     +          (D3_RB_RB_GBB)*t7*t8*t9+4.78125d-2*t18*t54*t7*t8*t9
     +          -4.78125d-2*gammabb*t45*t5*t8*t9+9.5625d-2*gammabb*t18
     +          *t25*t8*t9+4.78125d-2*t24*t5*t54*t7*t9-4.78125d-2*t28
     +          *t44*t5*t7*t9-4.78125d-2*t18*t24*t28*t7*t9-9.5625d-2
     +          *gammabb*t44*t5*s2f*t9-9.5625d-2*gammabb*t18*t24*s2f
     +          *t9+9.5625d-2*gammabb*t24*t25*t5*t9+1.325d-1*t11*t14
     +          *t47*t8+t76+t75+1.6562500000000002d-2*t14*t31*t46*t7
     +          -2.65d-1*gammabb*t11*t14*t24*t55+6.625000000000001d-2
     +          *gammabb*s3f*t14*t46+2.65d-1*gammabb*t11*t14*t31*t44
     +          -3.3125000000000004d-2*gammabb*t14*t26*t44+2.65d-1*s3f
     +          *t11*t14*t44-2.65d-1*t11*t14*t24*t26+s3c3(D3_RB_RB_GBB
     +          )*t15)*wght+2.0d+0*t56*wght+Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
              Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
              Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
              Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
              Cmat3(iq,D3_RA_GBB_GBB) = 1.0d+0*t1*(t80+2.65d-1*t11*t14
     +          *t52*t8+t79+t78+1.6562500000000002d-2*t14*t22*t24+s3c3
     +          (D3_RA_GBB_GBB)*t15)*wght+t77+Cmat3(iq,D3_RA_GBB_GBB)
              Cmat3(iq,D3_RB_GAA_GAA) = Cmat3(iq,D3_RB_GAA_GAA)
              Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
              Cmat3(iq,D3_RB_GAA_GBB) = Cmat3(iq,D3_RB_GAA_GBB)
              Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
              Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*t1*(-2.390625d-2*t5
     +          *s2c3(D3_RB_GBB_GBB)*t7*t8*t9-9.5625d-2*gammabb*t5*t54
     +          *t8*t9-4.78125d-2*t25*t5*t8*t9+2.390625d-2*t24*t5*t57
     +          *t7*t9+4.78125d-2*t24*t5*s2f*t9+9.5625d-2*gammabb*t24
     +          *t28*t5*t9+t80+2.65d-1*t11*t14*t55*t8+t79+t78-1.325d-1
     +          *gammabb*t11*t14*t24*t61-3.3125000000000004d-2*gammabb
     +          *t14*t31*t44-3.3125000000000004d-2*s3f*t14*t44-2.65d-1
     +          *t11*t14*t24*t31+1.6562500000000002d-2*t14*t24*t26
     +          +s3c3(D3_RB_GBB_GBB)*t15)*wght+t77+Cmat3(iq,
     +          D3_RB_GBB_GBB)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
              Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
              Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
              Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
              Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
              Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
              Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
              Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = 1.0d+0*t1*(t5*(-2.390625d-2
     +          *s2c3(D3_GBB_GBB_GBB)*t7-1.434375d-1*gammabb*t57
     +          -1.434375d-1*t28)*t8*t9+3.975d-1*t11*t14*t61*t8
     +          +4.96875d-2*t14*t24*t31+s3c3(D3_GBB_GBB_GBB)*t15)*wght
     +          +Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = 1.0d+0*t1*(t81
     +          -1.0600000000000001d+0*t12*t34*t40)*wght+2.0d+0*t65
     +          *wght+Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
              Mmat3(iq,D3_RA_RB_TB) = 1.0d+0*t1*(t81-4.78125d-2*t18
     +          *t25*t33*t7*t8+4.78125d-2*t18*t24*t33*s2f*t7
     +          -1.0600000000000001d+0*t12*t34*t42+2.65d-1*gammabb*t11
     +          *t22*t24*t34)*wght+1.0d+0*t66*wght+1.0d+0*t65*wght
     +          +Mmat3(iq,D3_RA_RB_TB)
              Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = 1.0d+0*t1*(t81+4.78125d-2*t33
     +          *t45*t5*t7*t8-9.5625d-2*t18*t25*t33*t7*t8+9.5625d-2
     +          *t33*t44*t5*s2f*t7+9.5625d-2*t18*t24*t33*s2f*t7
     +          -9.5625d-2*t24*t25*t33*t5*t7-3.3125000000000004d-2*s3f
     +          *t34*t46*t7-1.0600000000000001d+0*t12*t34*t47-5.3d-1
     +          *gammabb*s3f*t11*t34*t44+5.3d-1*gammabb*t11*t24*t26
     +          *t34)*wght+2.0d+0*t66*wght+Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
              Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
              Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
              Mmat3(iq,D3_RA_GBB_TB) = 1.0d+0*t1*(t84+t83-2.65d-1*t11
     +          *t22*t34*t8-1.0600000000000001d+0*t12*t34*t52)*wght
     +          +t82+Mmat3(iq,D3_RA_GBB_TB)
              Mmat3(iq,D3_RB_GAA_TA) = Mmat3(iq,D3_RB_GAA_TA)
              Mmat3(iq,D3_RB_GAA_TB) = Mmat3(iq,D3_RB_GAA_TB)
              Mmat3(iq,D3_RB_GAB_TA) = Mmat3(iq,D3_RB_GAB_TA)
              Mmat3(iq,D3_RB_GAB_TB) = Mmat3(iq,D3_RB_GAB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = 1.0d+0*t1*(t84+t83+4.78125d-2
     +          *t33*t5*t54*t7*t8+9.5625d-2*gammabb*t25*t33*t5*t8
     +          -2.65d-1*t11*t26*t34*t8-4.78125d-2*t24*t28*t33*t5*t7
     +          -9.5625d-2*gammabb*t24*t33*t5*s2f
     +          -1.0600000000000001d+0*t12*t34*t55
     +          +3.3125000000000004d-2*gammabb*s3f*t34*t44+2.65d-1
     +          *gammabb*t11*t24*t31*t34+2.65d-1*s3f*t11*t24*t34)*wght
     +          +t82+Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
              Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
              Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
              Mmat3(iq,D3_GAA_GBB_TB) = 1.0d+0*t1*(-2*t33*t5*t58*t8
     +          -5.3d-1*t11*t31*t34*t8-1.0600000000000001d+0*t12*t34
     +          *t61-3.3125000000000004d-2*s3f*t24*t34)*wght+Mmat3(iq,
     +          D3_GAA_GBB_TB)
              Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
              Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
              Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
              Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = 1.0d+0*t1*(7.95d-1*s3f*t11*t70
     +          *t8+6*t29*t5*t69*t8+3.18d+0*t12*t31*t70)*wght+Mmat3
     +          (iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_GAB_TA_TA) = Mmat3(iq,D3_GAB_TA_TA)
              Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
              Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = 1.0d+0*t1*(5.7375d-1*t5*s2f*t7
     +          *t8/taub**5-1.272d+1*s3f*t12/t13**5)*wght+Mmat3(iq,
     +          D3_GBB_TA_TA)
            else 
              t1 = rhob+rhoa
              t2 = param(1)
              t3 = 4.0d-50*rhoa
              t4 = 4.0d-50*rhob
              sr(R_A) = rhoa
              sr(R_B) = rhob
              sg(G_AA) = t3
              sg(G_AB) = 0
              sg(G_BB) = t4
              s1f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s1c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s1c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s1c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,ipol,1,1.0d0,sr,sg,s1f,
     +          s1a,s1a2,s1a3,s1c,s1c2,s1c3)
              t5 = s1f
              sr(R_A) = rhoa
              sr(R_B) = 0.0d+0
              sg(G_AA) = t3
              sg(G_AB) = 0.0d+0
              sg(G_BB) = 0.0d+0
              s2f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s2a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s2a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s2a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s2c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s2c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s2c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s2f,s2a,
     +          s2a2,s2a3,s2c,s2c2,s2c3)
              t6 = s2f
              t7 = 1/t1
              sr(R_A) = 0.0d+0
              sr(R_B) = rhob
              sg(G_AA) = 0.0d+0
              sg(G_AB) = 0.0d+0
              sg(G_BB) = t4
              s3f = 0.0d0
              call dcopy(NCOL_AMAT,0.0d0,0,s3a,1)
              call dcopy(NCOL_AMAT2,0.0d0,0,s3a2,1)
              call dcopy(NCOL_AMAT3,0.0d0,0,s3a3,1)
              call dcopy(NCOL_CMAT,0.0d0,0,s3c,1)
              call dcopy(NCOL_CMAT2,0.0d0,0,s3c2,1)
              call dcopy(NCOL_CMAT3,0.0d0,0,s3c3,1)
              call nwxcm_c_mpbe_d3(t2,tol_rho,2,1,1.0d0,sr,sg,s3f,s3a,
     +          s3a2,s3a3,s3c,s3c2,s3c3)
              t8 = s3f
              t9 = -3.825d-1*rhob*t7*s3f-3.825d-1*rhoa*s2f*t7
     +          +1.1325d+0*s1f
              t10 = 1/t1**2
              t11 = 3.825d-1*rhoa*t10*s2f
              t12 = 3.825d-1*rhob*t10*s3f
              t13 = s2a(D1_RA)
              t14 = -3.825d-1*s2f*t7-3.825d-1*rhoa*t13*t7+1.1325d+0
     +          *s1a(D1_RA)+t12+t11
              t15 = 1.0d+0*t9*wght
              t16 = s3a(D1_RB)
              t17 = -3.825d-1*t7*s3f-3.825d-1*rhob*t16*t7+1.1325d+0
     +          *s1a(D1_RB)+t12+t11
              t18 = 1/t1**3
              t19 = -7.65d-1*rhoa*t18*s2f
              t20 = -7.65d-1*rhob*t18*s3f
              t21 = s2a2(D2_RA_RA)
              t22 = -3.825d-1*rhoa*t21*t7-7.65d-1*t13*t7+7.65d-1*t10
     +          *s2f+1.1325d+0*s1a2(D2_RA_RA)+t20+t19+7.65d-1*rhoa*t10
     +          *t13
              t23 = 3.825d-1*t10*s3f+3.825d-1*t10*s2f+1.1325d+0*s1a2
     +          (D2_RA_RB)+t20+t19+3.825d-1*rhob*t10*t16+3.825d-1*rhoa
     +          *t10*t13
              t24 = s3a2(D2_RB_RB)
              t25 = 7.65d-1*t10*s3f-3.825d-1*rhob*t24*t7-7.65d-1*t16
     +          *t7+1.1325d+0*s1a2(D2_RB_RB)+t20+t19+7.65d-1*rhob*t10
     +          *t16
              t26 = 1/t1**4
              t27 = 2.2949999999999998d+0*rhoa*t26*s2f
              t28 = 2.2949999999999998d+0*rhob*t26*s3f
              t29 = 2.0d+0*t23*wght
              fnc(iq) = 1.0d+0*t1*t9*wght+fnc(iq)
              Amat(iq,D1_RA) = 1.0d+0*t1*t14*wght+t15+Amat(iq,D1_RA)
              Amat(iq,D1_RB) = 1.0d+0*t1*t17*wght+t15+Amat(iq,D1_RB)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = 1.0d+0*t1*t22*wght+2.0d+0*t14*wght
     +          +Amat2(iq,D2_RA_RA)
              Amat2(iq,D2_RA_RB) = 1.0d+0*t1*t23*wght+1.0d+0*t17*wght
     +          +1.0d+0*t14*wght+Amat2(iq,D2_RA_RB)
              Amat2(iq,D2_RB_RB) = 1.0d+0*t1*t25*wght+2.0d+0*t17*wght
     +          +Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
              Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
              Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
              Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
              Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
              Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
              Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
              Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = 1.0d+0*t1*(-3.825d-1*rhoa*s2a3
     +          (D3_RA_RA_RA)*t7-1.1474999999999999d+0*t21*t7
     +          -2.2949999999999998d+0*t18*s2f+1.1325d+0*s1a3
     +          (D3_RA_RA_RA)+t28+t27+1.1474999999999999d+0*rhoa*t10
     +          *t21-2.2949999999999998d+0*rhoa*t13*t18
     +          +2.2949999999999998d+0*t10*t13)*wght+3.0d+0*t22*wght
     +          +Amat3(iq,D3_RA_RA_RA)
              Amat3(iq,D3_RA_RA_RB) = 1.0d+0*t1*(-7.65d-1*t18*s3f
     +          -1.53d+0*t18*s2f+1.1325d+0*s1a3(D3_RA_RA_RB)+t28+t27
     +          +3.825d-1*rhoa*t10*t21-7.65d-1*rhob*t16*t18-1.53d+0
     +          *rhoa*t13*t18+7.65d-1*t10*t13)*wght+1.0d+0*t22*wght
     +          +t29+Amat3(iq,D3_RA_RA_RB)
              Amat3(iq,D3_RA_RB_RB) = 1.0d+0*t1*(-1.53d+0*t18*s3f
     +          -7.65d-1*t18*s2f+1.1325d+0*s1a3(D3_RA_RB_RB)+t28+t27
     +          +3.825d-1*rhob*t10*t24-1.53d+0*rhob*t16*t18-7.65d-1
     +          *rhoa*t13*t18+7.65d-1*t10*t16)*wght+1.0d+0*t25*wght
     +          +t29+Amat3(iq,D3_RA_RB_RB)
              Amat3(iq,D3_RB_RB_RB) = 1.0d+0*t1*(-3.825d-1*rhob*t7
     +          *s3a3(D3_RB_RB_RB)-2.2949999999999998d+0*t18*s3f
     +          -1.1474999999999999d+0*t24*t7+1.1325d+0*s1a3
     +          (D3_RB_RB_RB)+t28+t27+1.1474999999999999d+0*rhob*t10
     +          *t24-2.2949999999999998d+0*rhob*t16*t18
     +          +2.2949999999999998d+0*t10*t16)*wght+3.0d+0*t25*wght
     +          +Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
              Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
              Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
              Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
              Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
              Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
              Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
              Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
              Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
              Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
              Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
              Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
              Cmat3(iq,D3_RB_GAA_GAA) = Cmat3(iq,D3_RB_GAA_GAA)
              Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
              Cmat3(iq,D3_RB_GAA_GBB) = Cmat3(iq,D3_RB_GAA_GBB)
              Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
              Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
              Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
              Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
              Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
              Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
              Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
              Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
              Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
              Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
              Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
              Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
              Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
              Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
              Mmat3(iq,D3_RB_GAA_TA) = Mmat3(iq,D3_RB_GAA_TA)
              Mmat3(iq,D3_RB_GAA_TB) = Mmat3(iq,D3_RB_GAA_TB)
              Mmat3(iq,D3_RB_GAB_TA) = Mmat3(iq,D3_RB_GAB_TA)
              Mmat3(iq,D3_RB_GAB_TB) = Mmat3(iq,D3_RB_GAB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
              Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
              Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
              Mmat3(iq,D3_GAA_GBB_TB) = Mmat3(iq,D3_GAA_GBB_TB)
              Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
              Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
              Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
              Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_GAB_TA_TA) = Mmat3(iq,D3_GAB_TA_TA)
              Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
              Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            endif 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            else 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RA) = Amat(iq,D1_RA)
              Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
              Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
              Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
              Mmat2(iq,D2_RA_TA) = Mmat2(iq,D2_RA_TA)
              Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
              Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)
              Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
              Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
              Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
              Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
              Mmat3(iq,D3_RA_RA_TA) = Mmat3(iq,D3_RA_RA_TA)
              Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
              Mmat3(iq,D3_RA_GAA_TA) = Mmat3(iq,D3_RA_GAA_TA)
              Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
              Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
              Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
              Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)
              Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
              Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
              Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
              Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
              Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            endif 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            if (taua.gt.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
              Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
              Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
              Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
              Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
              Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
              Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            else 
              fnc(iq) = fnc(iq)
              Amat(iq,D1_RB) = Amat(iq,D1_RB)
              Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
              Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
              Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
              Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
              Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
              Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
              Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
              Mmat2(iq,D2_RB_TB) = Mmat2(iq,D2_RB_TB)
              Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
              Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)
              Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
              Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
              Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
              Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
              Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
              Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
              Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
              Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
              Mmat3(iq,D3_RB_RB_TB) = Mmat3(iq,D3_RB_RB_TB)
              Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
              Mmat3(iq,D3_RB_GBB_TB) = Mmat3(iq,D3_RB_GBB_TB)
              Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
              Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
              Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
              Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
              Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
              Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
              Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
              Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)
              Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
              Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
              Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
              Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
            endif 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_pkzb99_d3 
C> @}
