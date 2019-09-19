C> \ingroup nwxc
C> @{
C>
C> \file nwxcm_c_pbe.F
C> The nwxcm_c_pbe functional
C>
C> @}
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the nwxcm_c_pbe functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\
C>   {\it t_2} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_\alpha , \rho_\beta\right)\\\\
C>   {\it t_3} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_4} &=& \rho_\alpha-\rho_\beta\\\\
C>   {\it t_5} &=& {{1}\over{{\it t_3}}}\\\\
C>   {\it t_6} &=& 0.5\,\left({\it t_4}\,{\it t_5}
C>    +1.0\right)^{{{2}\over{3}}}+0.5\,\left(1.0-{\it t_4}\,{
C>    \it t_5}\right)^{{{2}\over{3}}}\\\\
C>   {\it t_7} &=& {\it t_6}^3\\\\
C>   {\it t_8} &=& {\it t_3}^{1.166666666666667}\\\\
C>   {\it t_9} &=& {{1}\over{{\it t_8}^{2.0}}}\\\\
C>   {\it t_{10}} &=& {{1}\over{{\it t_6}^{2.0}}}\\\\
C>   {\it t_{11}} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\
C>   {\it t_{12}} &=& {{1}\over{e^{{{32.16396844291483\,{\it t_2}
C>    \,{\it t_5}}\over{{\it t_7}}}}}}-1.0\\\\
C>   {\it t_{13}} &=& {{2.041389378054956\,{\it t_1}\,{\it t_9}\,{
C>    \it t_{10}}\,{\it t_{11}}}\over{{\it t_{12}}}}\\\\
C>   {\it t_{14}} &=& {\it t_1}^2\\\\
C>   {\it t_{15}} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_s , 0.0\right)\\\\
C>   {\it t_{16}} &=& \rho_s^{1.166666666666667}\\\\
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}^{2.0}}}\\\\
C>   {\it t_{18}} &=& {{1}\over{e^{{{64.32793688582966\,{
C>    \it t_{15}}}\over{\rho_s}}}}}-1.0\\\\
C>   {\it t_{19}} &=& {{3.240503646201146\,{\it t_1}\,{\it t_{17}}
C>    \,\sigma_{ss}}\over{{\it t_{18}}}}\\\\
C>   f &=& 1.0\,\left(0.03109069086965489\,{\it t_3}\,{\it t_7}
C>    \,\log \left({{2.041389378054956\,{\it t_1}\,{\it t_9}\,{
C>    \it t_{10}}\,{\it t_{11}}\,\left({\it t_{13}}
C>    +1.0\right)}\over{{{4.1672705928356\,{\it t_{14}}\,{
C>    \it t_{11}}^2}\over{{\it t_8}^{4.0}\,{\it t_6}^{4.0}\,{
C>    \it t_{12}}^2}}+{\it t_{13}}+1.0}}+1.0\right)+{
C>    \it t_2}\right)\\\\
C>   g &=& 0\\\\
C>   G &=& 1.0\,\left(0.01554534543482744\,\rho_s\,\log 
C>    \left({{3.240503646201146\,{\it t_1}\,{\it t_{17}}\,
C>    \sigma_{ss}\,\left({\it t_{19}}+1.0\right)}
C>    \over{{{10.50086388104292\,{\it t_{14}}\,\sigma_{ss}^2}
C>    \over{{\it t_{16}}^{4.0}\,{\it t_{18}}^2}}+{\it t_{19}}+1.0}}
C>    +1.0\right)+{\it t_{15}}\right)\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc [4,5,6].
C>
C> ### References ###
C>
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77
C>    , 3865 (1996)  , DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865 ">
C> 10.1103/PhysRevLett.77.3865 </a>
C>
C> [2] Maxima, a computer algebra system,
C> <a href="http://maxima.sourceforge.net/">
C> http://maxima.sourceforge.net/</a>
C>
C> [3] CLISP 2.49 (2010-07-07)
C>
C> [4] autoxc, revision 27694 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_pbe(param,tol_rho,ipol,nq,wght, rho,rgamma,
     +   fnc,Amat,Cmat) 
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
      double precision fnc(nq
     +   )      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
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
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            sr(R_T) = 2.0d0*rhoa
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call nwxcm_c_pw91lda(t1,tol_rho,ipol,1,1.0d0,sr,s1f,s1a)
            t2 = s1f
            t3 = rhoa**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhoa
            t6 = exp(-1.6081984221457415d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 1.6202518231005725d+0*gammaaa*t1*t4*t8
            t10 = t1**2
            t11 = gammaaa**2
            t12 = 1/t3**4.0d+0
            t13 = 1/t7**2
            t14 = s1a(D1_RA)
            t15 = 1.6202518231005728d+0*gammaaa*t1*t4*t8
            t16 = t15+1.0d+0
            t17 = t15+2.625215970260729d+0*t10*t11*t12*t13+1.0d+0
            t18 = 1/t17
            t19 = 1.6202518231005728d+0*gammaaa*t1*t16*t18*t4+1.0d+0
            t20 = 1/t19
            t21 = rhoa**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -1.890293793617335d+0*gammaaa*t1*t21*t22*t8
            t24 = 8.040992110728707d+0*s1f/rhoa**2
     +        -1.6081984221457415d+1*t14*t5
            t25 = -1.6202518231005728d+0*gammaaa*t1*t13*t24*t4*t6
            t26 = 1/t17**2
            fnc(iq) = 1.0d+0*(6.218138173930978d-2*rhoa*log
     +        (1.6202518231005725d+0*gammaaa*t1*t4*(t9+1.0d+0)/(t9
     +        +2.6252159702607286d+0*t10*t11*t12*t13+1.0d+0)+1.0d+0)
     +        +s1f)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(6.218138173930978d-2*rhoa*t20*(
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t26*t4*(
     +        -5.250431940521458d+0*t10*t11*t12*t24*t6/t7**3
     +        -6.125503930608368d+0*t10*t11*t13*t21/t3**5.0d+0+t25+t23
     +        )+1.6202518231005728d+0*gammaaa*t1*t18*(t25+t23)*t4
     +        -1.890293793617335d+0*gammaaa*t1*t16*t18*t21*t22)
     +        +3.109069086965489d-2*log(t19)+t14)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 6.218138173930978d-2*rhoa*t20*(
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t26*t4*
     +        (4.050629557751432d-1*t1*t4*t8+1.3126079851303646d+0
     +        *gammaaa*t10*t12*t13)+6.563039925651822d-1*gammaaa*t10
     +        *t12*t18*t8+4.050629557751432d-1*t1*t16*t18*t4)*wght
     +        +Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 6.218138173930978d-2*rhoa*t20*(
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t26*t4*
     +        (8.101259115502863d-1*t1*t4*t8+2.625215970260729d+0
     +        *gammaaa*t10*t12*t13)+1.3126079851303646d+0*gammaaa*t10
     +        *t12*t18*t8+8.101259115502863d-1*t1*t16*t18*t4)*wght
     +        +Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = rhob
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call nwxcm_c_pw91lda(t1,tol_rho,ipol,1,1.0d0,sr,s1f,s1a)
            t2 = s1f
            t3 = rhob+rhoa
            t4 = rhoa-rhob
            t5 = 1/t3
            t6 = 1.0d+0-t4*t5
            t7 = t4*t5+1.0d+0
            t8 = 5.0d-1*t7**6.666666666666666d-1+5.0d-1*t6*
     +        *6.666666666666666d-1
            t9 = t8**3
            t10 = gammabb+2.0d+0*gammaab+gammaaa
            t11 = t3**1.1666666666666667d+0
            t12 = 1/t11**2.0d+0
            t13 = 1/t8**2.0d+0
            t14 = 1/t9
            t15 = exp(-3.216396844291483d+1*t14*s1f*t5)
            t16 = t15-1.0d+0
            t17 = 1/t16
            t18 = 2.0413893780549558d+0*t1*t10*t12*t13*t17
            t19 = t18+1.0d+0
            t20 = t1**2
            t21 = t10**2
            t22 = 1/t11**4.0d+0
            t23 = 1/t8**4.0d+0
            t24 = 1/t16**2
            t25 = 4.1672705928355996d+0*t20*t21*t22*t23*t24+t18+1.0d+0
            t26 = 1/t25
            t27 = 2.0413893780549558d+0*t1*t10*t12*t13*t19*t26+1.0d+0
            t28 = log(t27)
            t29 = s1a(D1_RA)
            t30 = 1/t27
            t31 = 1/t3**2
            t32 = t31*t4
            t33 = -t5
            t34 = 1/t6**3.333333333333333d-1
            t35 = -t31*t4
            t36 = 1/t7**3.333333333333333d-1
            t37 = 3.333333333333333d-1*t36*(t5+t35)
     +        +3.333333333333333d-1*(t33+t32)*t34
            t38 = 1/t8**3.0d+0
            t39 = t3**1.6666666666666674d-1
            t40 = 1/t11**3.0d+0
            t41 = -4.763241882128231d+0*t1*t10*t13*t19*t26*t39*t40
            t42 = -4.0827787561099116d+0*t1*t10*t12*t17*t37*t38
            t43 = -4.763241882128231d+0*t1*t10*t13*t17*t39*t40
            t44 = 1/t8**4
            t45 = 3.216396844291483d+1*t14*s1f*t31
            t46 = 9.649190532874448d+1*s1f*t37*t44*t5
     +        -3.216396844291483d+1*t14*t29*t5+t45
            t47 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t46
            t48 = 1/t25**2
            t49 = 1/t8**5.0d+0
            t50 = -1.9447262766566134d+1*t20*t21*t23*t24*t39/t11**5.0d+0
            t51 = 1/t16**3
            t52 = t8**2
            t53 = 3.109069086965489d-2*t28*t9
            t54 = s1a(D1_RB)
            t55 = 3.333333333333333d-1*t34*(t5+t32)
     +        +3.333333333333333d-1*(t35+t33)*t36
            t56 = -4.0827787561099116d+0*t1*t10*t12*t17*t38*t55
            t57 = 9.649190532874448d+1*s1f*t44*t5*t55
     +        -3.216396844291483d+1*t14*t5*t54+t45
            t58 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t57
            t59 = 3.109069086965489d-2*t3*t30*(-2.0413893780549558d+0
     +        *t1*t10*t12*t13*t19*(8.334541185671199d+0*t10*t20*t22
     +        *t23*t24+2.0413893780549558d+0*t1*t12*t13*t17)*t48
     +        +4.1672705928355996d+0*t10*t17*t20*t22*t23*t26
     +        +2.0413893780549558d+0*t1*t12*t13*t19*t26)*t9*wght
            fnc(iq) = 1.0d+0*(3.109069086965489d-2*t28*t3*t9+s1f)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(3.109069086965489d-2*t3*t30*(
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t48*(
     +        -8.334541185671199d+0*t15*t20*t21*t22*t23*t46*t51+t50
     +        -1.6669082371342397d+1*t20*t21*t22*t24*t37*t49+t47+t43
     +        +t42)+2.0413893780549558d+0*t1*t10*t12*t13*t26*(t47+t43
     +        +t42)+t41-4.0827787561099116d+0*t1*t10*t12*t19*t26*t37
     +        *t38)*t9+t53+9.327207260896466d-2*t28*t3*t37*t52+t29)
     +        *wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(3.109069086965489d-2*t3*t30*(
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t48*(t58
     +        -8.334541185671199d+0*t15*t20*t21*t22*t23*t51*t57+t56
     +        -1.6669082371342397d+1*t20*t21*t22*t24*t49*t55+t50+t43)
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t58+t56+t43)
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t38*t55+t41)
     +        *t9+9.327207260896466d-2*t28*t3*t52*t55+t54+t53)*wght
     +        +Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t59+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t3*t30*(
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*
     +        (1.6669082371342397d+1*t10*t20*t22*t23*t24
     +        +4.0827787561099116d+0*t1*t12*t13*t17)*t48
     +        +8.334541185671199d+0*t10*t17*t20*t22*t23*t26
     +        +4.0827787561099116d+0*t1*t12*t13*t19*t26)*t9*wght+Cmat
     +        (iq,D1_GAB)
            Cmat(iq,D1_GBB) = t59+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call nwxcm_c_pw91lda(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a)
            t2 = s1f
            t3 = rhoa**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhoa
            t6 = exp(-6.432793688582966d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 3.240503646201146d+0*gammaaa*t1*t4*t8
            t10 = t9+1.0d+0
            t11 = t1**2
            t12 = gammaaa**2
            t13 = 1/t3**4.0d+0
            t14 = 1/t7**2
            t15 = t9+1.0500863881042921d+1*t11*t12*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 3.240503646201146d+0*gammaaa*t1*t10*t16*t4+1.0d+0
            t18 = log(t17)
            t19 = s1a(D1_RA)
            t20 = 1/t17
            t21 = rhoa**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -7.561175174469341d+0*gammaaa*t1*t21*t22*t8
            t24 = 6.432793688582966d+1*s1f/rhoa**2
     +        -6.432793688582966d+1*t19*t5
            t25 = -3.240503646201146d+0*gammaaa*t1*t14*t24*t4*t6
            t26 = 1/t15**2
            fnc(iq) = 1.0d+0*(s1f+1.5545345434827443d-2*rhoa*t18)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(1.5545345434827443d-2*rhoa*t20*(
     +        -3.240503646201146d+0*gammaaa*t1*t10*t26*t4*(
     +        -2.1001727762085842d+1*t11*t12*t13*t24*t6/t7**3
     +        -4.900403144486697d+1*t11*t12*t14*t21/t3**5.0d+0+t25+t23
     +        )+3.240503646201146d+0*gammaaa*t1*t16*(t25+t23)*t4
     +        -7.561175174469341d+0*gammaaa*t1*t10*t16*t21*t22)+t19
     +        +1.5545345434827443d-2*t18)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.5545345434827443d-2*rhoa*t20*(
     +        -3.240503646201146d+0*gammaaa*t1*t10*t26*t4*
     +        (3.240503646201146d+0*t1*t4*t8+2.1001727762085842d+1
     +        *gammaaa*t11*t13*t14)+1.0500863881042921d+1*gammaaa*t11
     +        *t13*t16*t8+3.240503646201146d+0*t1*t10*t16*t4)*wght
     +        +Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhob
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call nwxcm_c_pw91lda(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a)
            t2 = s1f
            t3 = rhob**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhob
            t6 = exp(-6.432793688582966d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 3.240503646201146d+0*gammabb*t1*t4*t8
            t10 = t9+1.0d+0
            t11 = t1**2
            t12 = gammabb**2
            t13 = 1/t3**4.0d+0
            t14 = 1/t7**2
            t15 = t9+1.0500863881042921d+1*t11*t12*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 3.240503646201146d+0*gammabb*t1*t10*t16*t4+1.0d+0
            t18 = log(t17)
            t19 = s1a(D1_RB)
            t20 = 1/t17
            t21 = rhob**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -7.561175174469341d+0*gammabb*t1*t21*t22*t8
            t24 = 6.432793688582966d+1*s1f/rhob**2
     +        -6.432793688582966d+1*t19*t5
            t25 = -3.240503646201146d+0*gammabb*t1*t14*t24*t4*t6
            t26 = 1/t15**2
            fnc(iq) = 1.0d+0*(s1f+1.5545345434827443d-2*rhob*t18)*wght
     +        +fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*(1.5545345434827443d-2*rhob*t20*(
     +        -3.240503646201146d+0*gammabb*t1*t10*t26*t4*(
     +        -2.1001727762085842d+1*t11*t12*t13*t24*t6/t7**3
     +        -4.900403144486697d+1*t11*t12*t14*t21/t3**5.0d+0+t25+t23
     +        )+3.240503646201146d+0*gammabb*t1*t16*(t25+t23)*t4
     +        -7.561175174469341d+0*gammabb*t1*t10*t16*t21*t22)+t19
     +        +1.5545345434827443d-2*t18)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.5545345434827443d-2*rhob*t20*(
     +        -3.240503646201146d+0*gammabb*t1*t10*t26*t4*
     +        (3.240503646201146d+0*t1*t4*t8+2.1001727762085842d+1
     +        *gammabb*t11*t13*t14)+1.0500863881042921d+1*gammabb*t11
     +        *t13*t16*t8+3.240503646201146d+0*t1*t10*t16*t4)*wght
     +        +Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_pbe 
C>
C> \brief Evaluate the nwxcm_c_pbe functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\
C>   {\it t_2} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_\alpha , \rho_\beta\right)\\\\
C>   {\it t_3} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_4} &=& \rho_\alpha-\rho_\beta\\\\
C>   {\it t_5} &=& {{1}\over{{\it t_3}}}\\\\
C>   {\it t_6} &=& 0.5\,\left({\it t_4}\,{\it t_5}
C>    +1.0\right)^{{{2}\over{3}}}+0.5\,\left(1.0-{\it t_4}\,{
C>    \it t_5}\right)^{{{2}\over{3}}}\\\\
C>   {\it t_7} &=& {\it t_6}^3\\\\
C>   {\it t_8} &=& {\it t_3}^{1.166666666666667}\\\\
C>   {\it t_9} &=& {{1}\over{{\it t_8}^{2.0}}}\\\\
C>   {\it t_{10}} &=& {{1}\over{{\it t_6}^{2.0}}}\\\\
C>   {\it t_{11}} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\
C>   {\it t_{12}} &=& {{1}\over{e^{{{32.16396844291483\,{\it t_2}
C>    \,{\it t_5}}\over{{\it t_7}}}}}}-1.0\\\\
C>   {\it t_{13}} &=& {{2.041389378054956\,{\it t_1}\,{\it t_9}\,{
C>    \it t_{10}}\,{\it t_{11}}}\over{{\it t_{12}}}}\\\\
C>   {\it t_{14}} &=& {\it t_1}^2\\\\
C>   {\it t_{15}} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_s , 0.0\right)\\\\
C>   {\it t_{16}} &=& \rho_s^{1.166666666666667}\\\\
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}^{2.0}}}\\\\
C>   {\it t_{18}} &=& {{1}\over{e^{{{64.32793688582966\,{
C>    \it t_{15}}}\over{\rho_s}}}}}-1.0\\\\
C>   {\it t_{19}} &=& {{3.240503646201146\,{\it t_1}\,{\it t_{17}}
C>    \,\sigma_{ss}}\over{{\it t_{18}}}}\\\\
C>   f &=& 1.0\,\left(0.03109069086965489\,{\it t_3}\,{\it t_7}
C>    \,\log \left({{2.041389378054956\,{\it t_1}\,{\it t_9}\,{
C>    \it t_{10}}\,{\it t_{11}}\,\left({\it t_{13}}
C>    +1.0\right)}\over{{{4.1672705928356\,{\it t_{14}}\,{
C>    \it t_{11}}^2}\over{{\it t_8}^{4.0}\,{\it t_6}^{4.0}\,{
C>    \it t_{12}}^2}}+{\it t_{13}}+1.0}}+1.0\right)+{
C>    \it t_2}\right)\\\\
C>   g &=& 0\\\\
C>   G &=& 1.0\,\left(0.01554534543482744\,\rho_s\,\log 
C>    \left({{3.240503646201146\,{\it t_1}\,{\it t_{17}}\,
C>    \sigma_{ss}\,\left({\it t_{19}}+1.0\right)}
C>    \over{{{10.50086388104292\,{\it t_{14}}\,\sigma_{ss}^2}
C>    \over{{\it t_{16}}^{4.0}\,{\it t_{18}}^2}}+{\it t_{19}}+1.0}}
C>    +1.0\right)+{\it t_{15}}\right)\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc [4,5,6].
C>
C> ### References ###
C>
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77
C>    , 3865 (1996)  , DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865 ">
C> 10.1103/PhysRevLett.77.3865 </a>
C>
C> [2] Maxima, a computer algebra system,
C> <a href="http://maxima.sourceforge.net/">
C> http://maxima.sourceforge.net/</a>
C>
C> [3] CLISP 2.49 (2010-07-07)
C>
C> [4] autoxc, revision 27694 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_pbe_d2(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,fnc,Amat,Amat2,Cmat,Cmat2) 
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
      double precision fnc(nq
     +   )      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
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
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            sr(R_T) = 2.0d0*rhoa
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call nwxcm_c_pw91lda_d2(t1,tol_rho,ipol,1,1.0d0,sr,s1f,
     +        s1a,s1a2)
            t2 = s1f
            t3 = rhoa**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhoa
            t6 = exp(-1.6081984221457415d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 1.6202518231005725d+0*gammaaa*t1*t4*t8
            t10 = t1**2
            t11 = gammaaa**2
            t12 = 1/t3**4.0d+0
            t13 = 1/t7**2
            t14 = s1a(D1_RA)
            t15 = 1.6202518231005728d+0*gammaaa*t1*t4*t8
            t16 = t15+1.0d+0
            t17 = t15+2.625215970260729d+0*t10*t11*t12*t13+1.0d+0
            t18 = 1/t17
            t19 = 1.6202518231005728d+0*gammaaa*t1*t16*t18*t4+1.0d+0
            t20 = 1/t19
            t21 = rhoa**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -1.890293793617335d+0*gammaaa*t1*t21*t22*t8
            t24 = 1/rhoa**2
            t25 = 8.040992110728707d+0*s1f*t24-1.6081984221457415d+1
     +        *t14*t5
            t26 = -1.6202518231005728d+0*gammaaa*t1*t13*t25*t4*t6
            t27 = t26+t23
            t28 = 1/t17**2
            t29 = 1/t3**5.0d+0
            t30 = 1/t7**3
            t31 = -5.250431940521458d+0*t10*t11*t12*t25*t30*t6
     +        -6.125503930608368d+0*t10*t11*t13*t21*t29+t26+t23
            t32 = -1.6202518231005728d+0*gammaaa*t1*t16*t28*t31*t4
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t27*t4
     +        -1.890293793617335d+0*gammaaa*t1*t16*t18*t21*t22
            t33 = log(t19)
            t34 = 4.050629557751432d-1*t1*t4*t8+1.3126079851303646d+0
     +        *gammaaa*t10*t12*t13
            t35 = 6.563039925651822d-1*gammaaa*t10*t12*t18*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t34*t4
     +        +4.050629557751432d-1*t1*t16*t18*t4
            t36 = 8.101259115502863d-1*t1*t4*t8+2.625215970260729d+0
     +        *gammaaa*t10*t12*t13
            t37 = 1.3126079851303646d+0*gammaaa*t10*t12*t18*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t36*t4
     +        +8.101259115502863d-1*t1*t16*t18*t4
            t38 = s1a2(D2_RA_RA)
            t39 = 6.218138173930978d-2*t32
            t40 = rhoa**3.333333333333335d-1
            t41 = 3.3080141388303363d+0*gammaaa*t1*t12*t16*t18*t40
            t42 = 1/rhoa**8.333333333333333d-1
            t43 = -1.5752448280144465d-1*gammaaa*t1*t16*t18*t22*t42
            t44 = -3.78058758723467d+0*gammaaa*t1*t18*t21*t22*t27
            t45 = 3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t28*t31
            t46 = -3.2405036462011455d+0*gammaaa*t1*t27*t28*t31*t4
            t47 = t25**2
            t48 = exp(-3.216396844291483d+1*s1f*t5)
            t49 = 3.2405036462011455d+0*gammaaa*t1*t30*t4*t47*t48
            t50 = 3.3080141388303363d+0*gammaaa*t1*t12*t40*t8
            t51 = -1.5752448280144465d-1*gammaaa*t1*t22*t42*t8
            t52 = 1.80027980344508d-1*gammaaa*t1*t24*t4*t8
            t53 = 3.78058758723467d+0*gammaaa*t1*t13*t21*t22*t25*t6
            t54 = -1.6202518231005728d+0*gammaaa*t1*t13*t4*t47*t6
            t55 = 1/rhoa**3
            t56 = 1.6081984221457415d+1*t14*t24
            t57 = t56-1.0721322814304943d+1*s1f*t55
     +        -1.6081984221457415d+1*t38*t5
            t58 = -1.6202518231005728d+0*gammaaa*t1*t13*t4*t57*t6
            t59 = 1.5751295821564373d+1*t10*t11*t12*t47*t48/t7**4
            t60 = 1/t3**6.0d+0
            t61 = 1.7866053130941079d+1*t10*t11*t13*t40*t60
            t62 = -5.10458660884031d-1*t10*t11*t13*t29*t42
            t63 = 2.4502015722433473d+1*t10*t11*t21*t25*t29*t30*t6
            t64 = -5.250431940521458d+0*t10*t11*t12*t30*t47*t6
            t65 = 1/t17**3
            t66 = 3.2405036462011455d+0*gammaaa*t1*t16*t31**2*t4*t65
            t67 = 1/t19**2
            t68 = -6.218138173930978d-2*rhoa*t32**2*t67
            t69 = s1a2(D2_RA_RB)
            t70 = -1.8002798034450807d-1*gammaaa*t1*t24*t4*t8
            t71 = -1.6081984221457415d+1*t5*t69+t56
     +        -5.360661407152472d+0*s1f*t55
            t72 = -1.6202518231005728d+0*gammaaa*t1*t13*t4*t6*t71
            t73 = -4.7257344840433374d-1*t1*t21*t22*t8
            t74 = -4.050629557751432d-1*t1*t13*t25*t4*t6
            t75 = 1.0d+0*(t20*(6.218138173930978d-2*rhoa*(
     +        -6.563039925651822d-1*gammaaa*t10*t12*t28*t31*t8
     +        -7.65687991326046d-1*gammaaa*t10*t18*t21*t29*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(t74+t73
     +        -2.625215970260729d+0*gammaaa*t10*t12*t25*t30*t6
     +        -3.062751965304184d+0*gammaaa*t10*t13*t21*t29)
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*(t74+t73)
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t31*t34*t4*t65
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t28*t34*t4
     +        -4.050629557751432d-1*t1*t16*t28*t31*t4
     +        +4.050629557751432d-1*t1*t18*t27*t4+1.890293793617335d+0
     +        *gammaaa*t1*t16*t21*t22*t28*t34-4.7257344840433374d-1*t1
     +        *t16*t18*t21*t22)+3.109069086965489d-2*t35)
     +        -6.218138173930978d-2*rhoa*t32*t35*t67)*wght
            t76 = -9.451468968086675d-1*t1*t21*t22*t8
            t77 = -8.101259115502863d-1*t1*t13*t25*t4*t6
            t78 = -8.101259115502863d-1*t1*t16*t28*t34*t4
            t79 = t1**3
            t80 = -1.3126079851303646d+0*gammaaa*t10*t12*t28*t34*t8
            t81 = 6.218138173930978d-2*rhoa*t20*(t80
     +        +3.281519962825911d-1*t10*t12*t18*t8
     +        -5.316888702309608d-1*gammaaa*t13*t16*t28*t60*t79+t78
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t34**2*t4*t65)
     +        *wght-6.218138173930978d-2*rhoa*t35**2*t67*wght
            fnc(iq) = 1.0d+0*(6.218138173930978d-2*rhoa*log
     +        (1.6202518231005725d+0*gammaaa*t1*t4*(t9+1.0d+0)/(t9
     +        +2.6252159702607286d+0*t10*t11*t12*t13+1.0d+0)+1.0d+0)
     +        +s1f)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(3.109069086965489d-2*t33
     +        +6.218138173930978d-2*rhoa*t20*t32+t14)*wght+Amat(iq,
     +        D1_RA)
            Cmat(iq,D1_GAA) = 6.218138173930978d-2*rhoa*t20*t35*wght
     +        +Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 6.218138173930978d-2*rhoa*t20*t37*wght
     +        +Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t68+t20*(6.218138173930978d-2
     +        *rhoa*(t66-1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*
     +        (t64+t63+t62+t61-5.250431940521458d+0*t10*t11*t12*t30
     +        *t57*t6+t59+t58+t54+t53+t52+t51+t50+t49
     +        +5.833813267246063d-1*t10*t11*t12*t13*t24)
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*(t58+t54+t53
     +        +t52+t51+t50+t49)+t46+t45+t44+t43+t41
     +        +1.80027980344508d-1*gammaaa*t1*t16*t18*t24*t4)+t39)
     +        -1.0363563623218294d-2*t33*t5+t38)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t20*(6.218138173930978d-2
     +        *rhoa*(-1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(t72
     +        -5.250431940521458d+0*t10*t11*t12*t30*t6*t71+t70+t64+t63
     +        +t62+t61+t59+t54+t53+t51+t50+t49-5.833813267246064d-1
     +        *t10*t11*t12*t13*t24)+1.6202518231005728d+0*gammaaa*t1
     +        *t18*t4*(t72+t70+t54+t53+t51+t50+t49)+t66+t46+t45+t44
     +        +t43+t41-1.8002798034450807d-1*gammaaa*t1*t16*t18*t24*t4
     +        )+t39)+t69+t68+1.0363563623218296d-2*t33*t5)*wght+Amat2
     +        (iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t75+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(t20*(6.218138173930978d-2
     +        *rhoa*(-1.3126079851303646d+0*gammaaa*t10*t12*t28*t31*t8
     +        -1.531375982652092d+0*gammaaa*t10*t18*t21*t29*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(t77+t76
     +        -5.250431940521458d+0*gammaaa*t10*t12*t25*t30*t6
     +        -6.125503930608368d+0*gammaaa*t10*t13*t21*t29)
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*(t77+t76)
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t31*t36*t4*t65
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t28*t36*t4
     +        -8.101259115502863d-1*t1*t16*t28*t31*t4
     +        +8.101259115502863d-1*t1*t18*t27*t4+1.890293793617335d+0
     +        *gammaaa*t1*t16*t21*t22*t28*t36-9.451468968086675d-1*t1
     +        *t16*t18*t21*t22)+3.109069086965489d-2*t37)
     +        -6.218138173930978d-2*rhoa*t32*t37*t67)*wght+Cmat2(iq,
     +        D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t75+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t81+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = 6.218138173930978d-2*rhoa*t20*(t80
     +        -6.563039925651822d-1*gammaaa*t10*t12*t28*t36*t8
     +        +6.563039925651822d-1*t10*t12*t18*t8
     +        -1.0633777404619216d+0*gammaaa*t13*t16*t28*t60*t79+t78
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t34*t36*t4*t65
     +        -4.050629557751432d-1*t1*t16*t28*t36*t4)*wght
     +        -6.218138173930978d-2*rhoa*t35*t37*t67*wght+Cmat2(iq,
     +        D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t81+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 6.218138173930978d-2*rhoa*t20*(
     +        -2.625215970260729d+0*gammaaa*t10*t12*t28*t36*t8
     +        +1.3126079851303646d+0*t10*t12*t18*t8
     +        -2.1267554809238431d+0*gammaaa*t13*t16*t28*t60*t79
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t36**2*t4*t65
     +        -1.6202518231005728d+0*t1*t16*t28*t36*t4)*wght
     +        -6.218138173930978d-2*rhoa*t37**2*t67*wght+Cmat2(iq,
     +        D2_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = rhob
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call nwxcm_c_pw91lda_d2(t1,tol_rho,ipol,1,1.0d0,sr,s1f,
     +        s1a,s1a2)
            t2 = s1f
            t3 = rhob+rhoa
            t4 = rhoa-rhob
            t5 = 1/t3
            t6 = 1.0d+0-t4*t5
            t7 = t4*t5+1.0d+0
            t8 = 5.0d-1*t7**6.666666666666666d-1+5.0d-1*t6*
     +        *6.666666666666666d-1
            t9 = t8**3
            t10 = gammabb+2.0d+0*gammaab+gammaaa
            t11 = t3**1.1666666666666667d+0
            t12 = 1/t11**2.0d+0
            t13 = 1/t8**2.0d+0
            t14 = 1/t9
            t15 = exp(-3.216396844291483d+1*t14*s1f*t5)
            t16 = t15-1.0d+0
            t17 = 1/t16
            t18 = 2.0413893780549558d+0*t1*t10*t12*t13*t17
            t19 = t18+1.0d+0
            t20 = t1**2
            t21 = t10**2
            t22 = 1/t11**4.0d+0
            t23 = 1/t8**4.0d+0
            t24 = 1/t16**2
            t25 = 4.1672705928355996d+0*t20*t21*t22*t23*t24+t18+1.0d+0
            t26 = 1/t25
            t27 = 2.0413893780549558d+0*t1*t10*t12*t13*t19*t26+1.0d+0
            t28 = log(t27)
            t29 = s1a(D1_RA)
            t30 = 1/t27
            t31 = 1/t3**2
            t32 = t31*t4
            t33 = -t5
            t34 = t33+t32
            t35 = 1/t6**3.333333333333333d-1
            t36 = -t31*t4
            t37 = t5+t36
            t38 = 1/t7**3.333333333333333d-1
            t39 = 3.333333333333333d-1*t37*t38+3.333333333333333d-1
     +        *t34*t35
            t40 = 1/t8**3.0d+0
            t41 = t3**1.6666666666666674d-1
            t42 = 1/t11**3.0d+0
            t43 = -4.763241882128231d+0*t1*t10*t13*t19*t26*t41*t42
            t44 = -4.0827787561099116d+0*t1*t10*t12*t17*t39*t40
            t45 = -4.763241882128231d+0*t1*t10*t13*t17*t41*t42
            t46 = 1/t8**4
            t47 = 3.216396844291483d+1*t14*s1f*t31
            t48 = 9.649190532874448d+1*s1f*t39*t46*t5
     +        -3.216396844291483d+1*t14*t29*t5+t47
            t49 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t48
            t50 = t49+t45+t44
            t51 = 1/t25**2
            t52 = 1/t8**5.0d+0
            t53 = 1/t11**5.0d+0
            t54 = -1.9447262766566134d+1*t20*t21*t23*t24*t41*t53
            t55 = 1/t16**3
            t56 = -8.334541185671199d+0*t15*t20*t21*t22*t23*t48*t55
     +        +t54-1.6669082371342397d+1*t20*t21*t22*t24*t39*t52+t49
     +        +t45+t44
            t57 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t56
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*t50+t43
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t39*t40
            t58 = t8**2
            t59 = 3.109069086965489d-2*t28*t9
            t60 = s1a(D1_RB)
            t61 = t5+t32
            t62 = t36+t33
            t63 = 3.333333333333333d-1*t38*t62+3.333333333333333d-1
     +        *t35*t61
            t64 = -4.0827787561099116d+0*t1*t10*t12*t17*t40*t63
            t65 = 9.649190532874448d+1*s1f*t46*t5*t63
     +        -3.216396844291483d+1*t14*t5*t60+t47
            t66 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t65
            t67 = t66+t64+t45
            t68 = t66-8.334541185671199d+0*t15*t20*t21*t22*t23*t55*t65
     +        +t64-1.6669082371342397d+1*t20*t21*t22*t24*t52*t63+t54+t45
            t69 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t68
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t40*t63+t43
            t70 = 8.334541185671199d+0*t10*t20*t22*t23*t24
     +        +2.0413893780549558d+0*t1*t12*t13*t17
            t71 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t70
     +        +4.1672705928355996d+0*t10*t17*t20*t22*t23*t26
     +        +2.0413893780549558d+0*t1*t12*t13*t19*t26
            t72 = 3.109069086965489d-2*t3*t30*t71*t9*wght
            t73 = 1.6669082371342397d+1*t10*t20*t22*t23*t24
     +        +4.0827787561099116d+0*t1*t12*t13*t17
            t74 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t73
     +        +8.334541185671199d+0*t10*t17*t20*t22*t23*t26
     +        +4.0827787561099116d+0*t1*t12*t13*t19*t26
            t75 = s1a2(D2_RA_RA)
            t76 = t39**2
            t77 = 1/t6**1.3333333333333333d+0
            t78 = 1/t3**3
            t79 = -2*t4*t78
            t80 = 2*t31
            t81 = 1/t7**1.3333333333333333d+0
            t82 = 2*t4*t78
            t83 = -2*t31
            t84 = 3.333333333333333d-1*t38*(t83+t82)
     +        -1.111111111111111d-1*t37**2*t81+3.333333333333333d-1
     +        *t35*(t80+t79)-1.111111111111111d-1*t34**2*t77
            t85 = t3**3.333333333333335d-1
            t86 = 1.667134658744881d+1*t1*t10*t13*t19*t22*t26*t85
            t87 = 1/t3**8.333333333333333d-1
            t88 = -7.938736470213722d-1*t1*t10*t13*t19*t26*t42*t87
            t89 = t48**2
            t90 = exp(-6.432793688582966d+1*t14*s1f*t5)
            t91 = 4.0827787561099116d+0*t1*t10*t12*t13*t55*t89*t90
            t92 = 1.2248336268329735d+1*t1*t10*t12*t17*t23*t76
            t93 = -4.0827787561099116d+0*t1*t10*t12*t17*t40*t84
            t94 = 1.9052967528512924d+1*t1*t10*t17*t39*t40*t41*t42
            t95 = 1.667134658744881d+1*t1*t10*t13*t17*t22*t85
            t96 = -7.938736470213722d-1*t1*t10*t13*t17*t42*t87
            t97 = 8.165557512219823d+0*t1*t10*t12*t15*t24*t39*t40*t48
            t98 = 9.526483764256461d+0*t1*t10*t13*t15*t24*t41*t42*t48
            t99 = 1/t8**5
            t100 = -6.432793688582966d+1*t14*s1f*t78
            t101 = -3.859676213149779d+2*s1f*t5*t76*t99
     +        +9.649190532874448d+1*s1f*t46*t5*t84
     +        -3.216396844291483d+1*t14*t5*t75+1.9298381065748896d+2
     +        *t29*t39*t46*t5-1.9298381065748896d+2*s1f*t31*t39*t46
     +        +6.432793688582966d+1*t14*t29*t31+t100
            t102 = -2.0413893780549558d+0*t1*t10*t101*t12*t13*t15*t24
            t103 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t89
            t104 = 1/t16**4
            t105 = 1/t8**6.0d+0
            t106 = 1/t11**6.0d+0
            t107 = 1.1344236613830246d+2*t106*t20*t21*t23*t24*t85
            t108 = -3.241210461094357d+0*t20*t21*t23*t24*t53*t87
            t109 = 1/t25**3
            t110 = 1/t27**2
            t111 = s1a2(D2_RA_RB)
            t112 = -1.111111111111111d-1*t37*t62*t81
     +        +6.666666666666666d-1*t38*t4*t78-6.666666666666666d-1
     +        *t35*t4*t78-1.111111111111111d-1*t34*t61*t77
            t113 = 4.0827787561099116d+0*t1*t10*t12*t13*t48*t55*t65*t90
            t114 = 1.2248336268329735d+1*t1*t10*t12*t17*t23*t39*t63
            t115 = -4.0827787561099116d+0*t1*t10*t112*t12*t17*t40
            t116 = 9.526483764256461d+0*t1*t10*t17*t40*t41*t42*t63
            t117 = 9.526483764256461d+0*t1*t10*t17*t39*t40*t41*t42
            t118 = 4.0827787561099116d+0*t1*t10*t12*t15*t24*t40*t48*t63
            t119 = 4.763241882128231d+0*t1*t10*t13*t15*t24*t41*t42*t48
            t120 = -3.859676213149779d+2*s1f*t39*t5*t63*t99
     +        +9.649190532874448d+1*t29*t46*t5*t63
     +        -9.649190532874448d+1*s1f*t31*t46*t63
     +        +9.649190532874448d+1*t39*t46*t5*t60
     +        +3.216396844291483d+1*t14*t31*t60+9.649190532874448d+1
     +        *t112*s1f*t46*t5-3.216396844291483d+1*t111*t14*t5
     +        -9.649190532874448d+1*s1f*t31*t39*t46
     +        +3.216396844291483d+1*t14*t29*t31+t100
            t121 = -2.0413893780549558d+0*t1*t10*t12*t120*t13*t15*t24
            t122 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t48*t65
            t123 = 4.0827787561099116d+0*t1*t10*t12*t15*t24*t39*t40*t65
            t124 = 4.763241882128231d+0*t1*t10*t13*t15*t24*t41*t42*t65
            t125 = s1a2(D2_RB_RB)
            t126 = t63**2
            t127 = 3.333333333333333d-1*t35*(t83+t79)
     +        +3.333333333333333d-1*t38*(t82+t80)-1.111111111111111d-1
     +        *t62**2*t81-1.111111111111111d-1*t61**2*t77
            t128 = t65**2
            t129 = 4.0827787561099116d+0*t1*t10*t12*t128*t13*t55*t90
            t130 = 1.2248336268329735d+1*t1*t10*t12*t126*t17*t23
            t131 = -4.0827787561099116d+0*t1*t10*t12*t127*t17*t40
            t132 = 1.9052967528512924d+1*t1*t10*t17*t40*t41*t42*t63
            t133 = 8.165557512219823d+0*t1*t10*t12*t15*t24*t40*t63*t65
            t134 = 9.526483764256461d+0*t1*t10*t13*t15*t24*t41*t42*t65
            t135 = -3.859676213149779d+2*t126*s1f*t5*t99
     +        +1.9298381065748896d+2*t46*t5*t60*t63
     +        -1.9298381065748896d+2*s1f*t31*t46*t63
     +        +6.432793688582966d+1*t14*t31*t60+9.649190532874448d+1
     +        *t127*s1f*t46*t5-3.216396844291483d+1*t125*t14*t5+t100
            t136 = -2.0413893780549558d+0*t1*t10*t12*t13*t135*t15*t24
            t137 = -2.0413893780549558d+0*t1*t10*t12*t128*t13*t15*t24
            t138 = 3.109069086965489d-2*t71*t9
            t139 = 4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t70
            t140 = -4.763241882128231d+0*t1*t13*t19*t26*t41*t42
            t141 = -9.723631383283067d+0*t10*t17*t20*t23*t26*t41*t53
            t142 = -4.0827787561099116d+0*t1*t12*t17*t39*t40
            t143 = -4.763241882128231d+0*t1*t13*t17*t41*t42
            t144 = -2.0413893780549558d+0*t1*t12*t13*t15*t24*t48
            t145 = -3.889452553313227d+1*t10*t20*t23*t24*t41*t53
            t146 = 1.0d+0*(t30*(3.109069086965489d-2*t3*
     +        (4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t56*t70
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t50*t51*t70
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t39*t40*t51*t70
     +        -4.1672705928355996d+0*t10*t17*t20*t22*t23*t51*t56
     +        -2.0413893780549558d+0*t1*t12*t13*t19*t51*t56
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -1.6669082371342397d+1*t10*t15*t20*t22*t23*t48*t55
     +        -3.3338164742684795d+1*t10*t20*t22*t24*t39*t52+t145+t144
     +        +t143+t142)-8.334541185671199d+0*t10*t17*t20*t22*t26*t39
     +        *t52+2.0413893780549558d+0*t1*t12*t13*t26*t50
     +        -4.0827787561099116d+0*t1*t12*t19*t26*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t12*t13*(t144+t143+t142)
     +        *t26+t141+t140+t139)*t9+9.327207260896466d-2*t3*t39*t58
     +        *t71+t138)-3.109069086965489d-2*t110*t3*t57*t71*t9)*wght
            t147 = 3.109069086965489d-2*t74*t9
            t148 = 4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t73
            t149 = -9.526483764256461d+0*t1*t13*t19*t26*t41*t42
            t150 = -1.9447262766566134d+1*t10*t17*t20*t23*t26*t41*t53
            t151 = -8.165557512219823d+0*t1*t12*t17*t39*t40
            t152 = -9.526483764256461d+0*t1*t13*t17*t41*t42
            t153 = -4.0827787561099116d+0*t1*t12*t13*t15*t24*t48
            t154 = -7.778905106626454d+1*t10*t20*t23*t24*t41*t53
            t155 = -4.0827787561099116d+0*t1*t12*t17*t40*t63
            t156 = -2.0413893780549558d+0*t1*t12*t13*t15*t24*t65
            t157 = 1.0d+0*(t30*(3.109069086965489d-2*t3*
     +        (4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t68*t70
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t51*t67*t70
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t63*t70
     +        -4.1672705928355996d+0*t10*t17*t20*t22*t23*t51*t68
     +        -2.0413893780549558d+0*t1*t12*t13*t19*t51*t68
     +        +2.0413893780549558d+0*t1*t12*t13*t26*t67
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -1.6669082371342397d+1*t10*t15*t20*t22*t23*t55*t65
     +        -3.3338164742684795d+1*t10*t20*t22*t24*t52*t63+t156+t155
     +        +t145+t143)-8.334541185671199d+0*t10*t17*t20*t22*t26*t52
     +        *t63-4.0827787561099116d+0*t1*t12*t19*t26*t40*t63
     +        +2.0413893780549558d+0*t1*t10*t12*t13*(t156+t155+t143)
     +        *t26+t141+t140+t139)*t9+9.327207260896466d-2*t3*t58*t63
     +        *t71+t138)-3.109069086965489d-2*t110*t3*t69*t71*t9)*wght
            t158 = -8.165557512219823d+0*t1*t12*t17*t40*t63
            t159 = -4.0827787561099116d+0*t1*t12*t13*t15*t24*t65
            t160 = -4.0827787561099116d+0*t1*t12*t13*t19*t51*t70
            t161 = t1**3
            t162 = -8.334541185671199d+0*t10*t17*t20*t22*t23*t51*t70
            t163 = 3.109069086965489d-2*t3*t30*(4.0827787561099116d+0
     +        *t1*t10*t109*t12*t13*t19*t70**2-1.7014043847390745d+1
     +        *t10*t105*t106*t161*t19*t24*t51+8.334541185671199d+0*t17
     +        *t20*t22*t23*t26+t162+t160)*t9*wght-3.109069086965489d-2
     +        *t110*t3*t71**2*t9*wght
            t164 = 3.109069086965489d-2*t3*t30*(4.0827787561099116d+0
     +        *t1*t10*t109*t12*t13*t19*t70*t73-4.1672705928355996d+0
     +        *t10*t17*t20*t22*t23*t51*t73-2.0413893780549558d+0*t1
     +        *t12*t13*t19*t51*t73-3.402808769478149d+1*t10*t105*t106
     +        *t161*t19*t24*t51+1.6669082371342397d+1*t17*t20*t22*t23
     +        *t26+t162+t160)*t9*wght-3.109069086965489d-2*t110*t3*t71
     +        *t74*t9*wght
            fnc(iq) = 1.0d+0*(3.109069086965489d-2*t28*t3*t9+s1f)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(3.109069086965489d-2*t3*t30*t57
     +        *t9+t59+9.327207260896466d-2*t28*t3*t39*t58+t29)*wght
     +        +Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(3.109069086965489d-2*t3*t30*t69
     +        *t9+9.327207260896466d-2*t28*t3*t58*t63+t60+t59)*wght
     +        +Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t72+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t3*t30*t74*t9*wght
     +        +Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = t72+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t30*(3.109069086965489d-2*t3
     +        *t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(t98
     +        +t97+t96+t95+t94+t93+t92+t91+2.50036235570136d+1*t104
     +        *t20*t21*t22*t23*t89*t90-8.334541185671199d+0*t15*t20
     +        *t21*t22*t23*t55*t89-1.6669082371342397d+1*t20*t21*t22
     +        *t24*t52*t84+8.334541185671199d+1*t105*t20*t21*t22*t24
     +        *t76+7.778905106626454d+1*t15*t20*t21*t23*t41*t48*t53
     +        *t55+6.667632948536959d+1*t15*t20*t21*t22*t39*t48*t52
     +        *t55-8.334541185671199d+0*t101*t15*t20*t21*t22*t23*t55
     +        +1.5557810213252907d+2*t20*t21*t24*t39*t41*t52*t53+t108
     +        +t107+t103+t102)+2.0413893780549558d+0*t1*t10*t12*t13
     +        *t26*(t98+t97+t96+t95+t94+t93+t92+t91+t103+t102)+t88+t86
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t40*t84
     +        +1.2248336268329735d+1*t1*t10*t12*t19*t23*t26*t76
     +        +4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t56**2
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t50*t51*t56
     +        +9.526483764256461d+0*t1*t10*t13*t19*t41*t42*t51*t56
     +        +8.165557512219823d+0*t1*t10*t12*t19*t39*t40*t51*t56
     +        -9.526483764256461d+0*t1*t10*t13*t26*t41*t42*t50
     +        -8.165557512219823d+0*t1*t10*t12*t26*t39*t40*t50
     +        +1.9052967528512924d+1*t1*t10*t19*t26*t39*t40*t41*t42)
     +        +6.218138173930978d-2*t57*t9+1.8654414521792934d-1*t3
     +        *t39*t57*t58)-3.109069086965489d-2*t110*t3*t57**2*t9
     +        +9.327207260896466d-2*t28*t3*t58*t84
     +        +1.8654414521792934d-1*t28*t3*t76*t8+t75
     +        +1.8654414521792934d-1*t28*t39*t58)*wght+Amat2(iq,
     +        D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t30*(3.109069086965489d-2*t3
     +        *t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(t96
     +        +t95+2.50036235570136d+1*t104*t20*t21*t22*t23*t48*t65
     +        *t90+3.889452553313227d+1*t15*t20*t21*t23*t41*t53*t55
     +        *t65+3.3338164742684795d+1*t15*t20*t21*t22*t39*t52*t55
     +        *t65-8.334541185671199d+0*t15*t20*t21*t22*t23*t48*t55
     +        *t65+3.3338164742684795d+1*t15*t20*t21*t22*t48*t52*t55
     +        *t63+7.778905106626454d+1*t20*t21*t24*t41*t52*t53*t63
     +        +8.334541185671199d+1*t105*t20*t21*t22*t24*t39*t63
     +        +3.889452553313227d+1*t15*t20*t21*t23*t41*t48*t53*t55
     +        -8.334541185671199d+0*t120*t15*t20*t21*t22*t23*t55
     +        +7.778905106626454d+1*t20*t21*t24*t39*t41*t52*t53
     +        -1.6669082371342397d+1*t112*t20*t21*t22*t24*t52+t124
     +        +t123+t122+t121+t119+t118+t117+t116+t115+t114+t113+t108
     +        +t107)+2.0413893780549558d+0*t1*t10*t12*t13*t26*(t96+t95
     +        +t124+t123+t122+t121+t119+t118+t117+t116+t115+t114+t113)
     +        +t88+t86+4.0827787561099116d+0*t1*t10*t109*t12*t13*t19
     +        *t56*t68-2.0413893780549558d+0*t1*t10*t12*t13*t50*t51
     +        *t68+4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t68
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t39*t40*t51*t68
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t51*t56*t67
     +        -4.763241882128231d+0*t1*t10*t13*t26*t41*t42*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t26*t39*t40*t67
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t56*t63
     +        -4.0827787561099116d+0*t1*t10*t12*t26*t40*t50*t63
     +        +9.526483764256461d+0*t1*t10*t19*t26*t40*t41*t42*t63
     +        +1.2248336268329735d+1*t1*t10*t12*t19*t23*t26*t39*t63
     +        +4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t56
     +        -4.763241882128231d+0*t1*t10*t13*t26*t41*t42*t50
     +        +9.526483764256461d+0*t1*t10*t19*t26*t39*t40*t41*t42
     +        -4.0827787561099116d+0*t1*t10*t112*t12*t19*t26*t40)
     +        +3.109069086965489d-2*t69*t9+3.109069086965489d-2*t57*t9
     +        +9.327207260896466d-2*t3*t39*t58*t69
     +        +9.327207260896466d-2*t3*t57*t58*t63)
     +        -3.109069086965489d-2*t110*t3*t57*t69*t9
     +        +1.8654414521792934d-1*t28*t3*t39*t63*t8
     +        +9.327207260896466d-2*t28*t58*t63+9.327207260896466d-2
     +        *t28*t39*t58+9.327207260896466d-2*t112*t28*t3*t58+t111)
     +        *wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(t30*(3.109069086965489d-2*t3
     +        *t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(t96
     +        +t95+2.50036235570136d+1*t104*t128*t20*t21*t22*t23*t90
     +        +6.667632948536959d+1*t15*t20*t21*t22*t52*t55*t63*t65
     +        +7.778905106626454d+1*t15*t20*t21*t23*t41*t53*t55*t65
     +        +1.5557810213252907d+2*t20*t21*t24*t41*t52*t53*t63
     +        -8.334541185671199d+0*t135*t15*t20*t21*t22*t23*t55
     +        -8.334541185671199d+0*t128*t15*t20*t21*t22*t23*t55
     +        -1.6669082371342397d+1*t127*t20*t21*t22*t24*t52
     +        +8.334541185671199d+1*t105*t126*t20*t21*t22*t24+t137
     +        +t136+t134+t133+t132+t131+t130+t129+t108+t107)
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t96+t95+t137
     +        +t136+t134+t133+t132+t131+t130+t129)+t88+t86
     +        +4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t68**2
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t51*t67*t68
     +        +8.165557512219823d+0*t1*t10*t12*t19*t40*t51*t63*t68
     +        +9.526483764256461d+0*t1*t10*t13*t19*t41*t42*t51*t68
     +        -8.165557512219823d+0*t1*t10*t12*t26*t40*t63*t67
     +        -9.526483764256461d+0*t1*t10*t13*t26*t41*t42*t67
     +        +1.9052967528512924d+1*t1*t10*t19*t26*t40*t41*t42*t63
     +        -4.0827787561099116d+0*t1*t10*t12*t127*t19*t26*t40
     +        +1.2248336268329735d+1*t1*t10*t12*t126*t19*t23*t26)
     +        +6.218138173930978d-2*t69*t9+1.8654414521792934d-1*t3
     +        *t58*t63*t69)-3.109069086965489d-2*t110*t3*t69**2*t9
     +        +1.8654414521792934d-1*t126*t28*t3*t8
     +        +1.8654414521792934d-1*t28*t58*t63+9.327207260896466d-2
     +        *t127*t28*t3*t58+t125)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t146+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(t30*(3.109069086965489d-2*t3
     +        *(4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t56*t73
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t50*t51*t73
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t39*t40*t51*t73
     +        -8.334541185671199d+0*t10*t17*t20*t22*t23*t51*t56
     +        -4.0827787561099116d+0*t1*t12*t13*t19*t51*t56
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -3.3338164742684795d+1*t10*t15*t20*t22*t23*t48*t55
     +        -6.667632948536959d+1*t10*t20*t22*t24*t39*t52+t154+t153
     +        +t152+t151)-1.6669082371342397d+1*t10*t17*t20*t22*t26
     +        *t39*t52+4.0827787561099116d+0*t1*t12*t13*t26*t50
     +        -8.165557512219823d+0*t1*t12*t19*t26*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t12*t13*(t153+t152+t151)
     +        *t26+t150+t149+t148)*t9+9.327207260896466d-2*t3*t39*t58
     +        *t74+t147)-3.109069086965489d-2*t110*t3*t57*t74*t9)*wght
     +        +Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t146+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t157+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = 1.0d+0*(t30*(3.109069086965489d-2*t3
     +        *(4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t68*t73
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t51*t67*t73
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t63*t73
     +        -8.334541185671199d+0*t10*t17*t20*t22*t23*t51*t68
     +        -4.0827787561099116d+0*t1*t12*t13*t19*t51*t68
     +        +4.0827787561099116d+0*t1*t12*t13*t26*t67
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -3.3338164742684795d+1*t10*t15*t20*t22*t23*t55*t65
     +        -6.667632948536959d+1*t10*t20*t22*t24*t52*t63+t159+t158
     +        +t154+t152)-1.6669082371342397d+1*t10*t17*t20*t22*t26
     +        *t52*t63-8.165557512219823d+0*t1*t12*t19*t26*t40*t63
     +        +2.0413893780549558d+0*t1*t10*t12*t13*(t159+t158+t152)
     +        *t26+t150+t149+t148)*t9+9.327207260896466d-2*t3*t58*t63
     +        *t74+t147)-3.109069086965489d-2*t110*t3*t69*t74*t9)*wght
     +        +Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t157+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t163+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t164+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t163+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = -3.109069086965489d-2*t110*t3*t74*
     +        *2*t9*wght+3.109069086965489d-2*t3*t30*
     +        (4.0827787561099116d+0*t1*t10*t109*t12*t13*t19*t73**2
     +        -1.6669082371342397d+1*t10*t17*t20*t22*t23*t51*t73
     +        -8.165557512219823d+0*t1*t12*t13*t19*t51*t73
     +        -6.805617538956298d+1*t10*t105*t106*t161*t19*t24*t51
     +        +3.3338164742684795d+1*t17*t20*t22*t23*t26)*t9*wght
     +        +Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = t164+Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t163+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call nwxcm_c_pw91lda_d2(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a,
     +        s1a2)
            t2 = s1f
            t3 = rhoa**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhoa
            t6 = exp(-6.432793688582966d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 3.240503646201146d+0*gammaaa*t1*t4*t8
            t10 = t9+1.0d+0
            t11 = t1**2
            t12 = gammaaa**2
            t13 = 1/t3**4.0d+0
            t14 = 1/t7**2
            t15 = t9+1.0500863881042921d+1*t11*t12*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 3.240503646201146d+0*gammaaa*t1*t10*t16*t4+1.0d+0
            t18 = log(t17)
            t19 = s1a(D1_RA)
            t20 = 1/t17
            t21 = rhoa**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -7.561175174469341d+0*gammaaa*t1*t21*t22*t8
            t24 = 1/rhoa**2
            t25 = 6.432793688582966d+1*s1f*t24-6.432793688582966d+1
     +        *t19*t5
            t26 = -3.240503646201146d+0*gammaaa*t1*t14*t25*t4*t6
            t27 = t26+t23
            t28 = 1/t15**2
            t29 = 1/t3**5.0d+0
            t30 = 1/t7**3
            t31 = -2.1001727762085842d+1*t11*t12*t13*t25*t30*t6
     +        -4.900403144486697d+1*t11*t12*t14*t21*t29+t26+t23
            t32 = -3.240503646201146d+0*gammaaa*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*gammaaa*t1*t16*t27*t4
     +        -7.561175174469341d+0*gammaaa*t1*t10*t16*t21*t22
            t33 = 3.240503646201146d+0*t1*t4*t8+2.1001727762085842d+1
     +        *gammaaa*t11*t13*t14
            t34 = 1.0500863881042921d+1*gammaaa*t11*t13*t16*t8
     +        -3.240503646201146d+0*gammaaa*t1*t10*t28*t33*t4
     +        +3.240503646201146d+0*t1*t10*t16*t4
            t35 = s1a2(D2_RA_RA)
            t36 = rhoa**3.333333333333335d-1
            t37 = 1/rhoa**8.333333333333333d-1
            t38 = t25**2
            t39 = exp(-1.2865587377165932d+2*s1f*t5)
            t40 = 6.481007292402292d+0*gammaaa*t1*t30*t38*t39*t4
            t41 = 2.6464113110642695d+1*gammaaa*t1*t13*t36*t8
            t42 = -1.2601958624115575d+0*gammaaa*t1*t22*t37*t8
            t43 = 1.5122350348938682d+1*gammaaa*t1*t14*t21*t22*t25*t6
            t44 = -3.240503646201146d+0*gammaaa*t1*t14*t38*t4*t6
            t45 = -6.432793688582966d+1*t35*t5+1.2865587377165932d+2
     +        *t19*t24-1.2865587377165932d+2*s1f/rhoa**3
            t46 = -3.240503646201146d+0*gammaaa*t1*t14*t4*t45*t6
            t47 = 1/t3**6.0d+0
            t48 = 1/t15**3
            t49 = 1/t17**2
            t50 = -7.561175174469341d+0*t1*t21*t22*t8
            t51 = -3.240503646201146d+0*t1*t14*t25*t4*t6
            fnc(iq) = 1.0d+0*(s1f+1.5545345434827443d-2*rhoa*t18)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(1.5545345434827443d-2*rhoa*t20
     +        *t32+t19+1.5545345434827443d-2*t18)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.5545345434827443d-2*rhoa*t20*t34*wght
     +        +Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t20*(1.5545345434827443d-2
     +        *rhoa*(-3.240503646201146d+0*gammaaa*t1*t10*t28*t4*
     +        (6.300518328625753d+1*t11*t12*t13*t38*t39/t7**4
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t45*t6
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t38*t6
     +        +1.9601612577946786d+2*t11*t12*t21*t25*t29*t30*t6
     +        +2.8585685009505735d+2*t11*t12*t14*t36*t47+t46+t44+t43
     +        +t42+t41+t40-8.167338574144498d+0*t11*t12*t14*t29*t37)
     +        +6.481007292402292d+0*gammaaa*t1*t10*t31**2*t4*t48
     +        +3.240503646201146d+0*gammaaa*t1*t16*t4*(t46+t44+t43+t42
     +        +t41+t40)-6.481007292402292d+0*gammaaa*t1*t27*t28*t31*t4
     +        -1.2601958624115575d+0*gammaaa*t1*t10*t16*t22*t37
     +        +2.6464113110642695d+1*gammaaa*t1*t10*t13*t16*t36
     +        +1.5122350348938682d+1*gammaaa*t1*t10*t21*t22*t28*t31
     +        -1.5122350348938682d+1*gammaaa*t1*t16*t21*t22*t27)
     +        +3.1090690869654886d-2*t32)-1.5545345434827443d-2*rhoa
     +        *t32**2*t49+t35)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(t20*(1.5545345434827443d-2
     +        *rhoa*(-1.0500863881042921d+1*gammaaa*t11*t13*t28*t31*t8
     +        -2.4502015722433484d+1*gammaaa*t11*t16*t21*t29*t8
     +        -3.240503646201146d+0*gammaaa*t1*t10*t28*t4*(
     +        -4.2003455524171684d+1*gammaaa*t11*t13*t25*t30*t6+t51
     +        +t50-9.800806288973394d+1*gammaaa*t11*t14*t21*t29)
     +        +3.240503646201146d+0*gammaaa*t1*t16*t4*(t51+t50)
     +        +6.481007292402292d+0*gammaaa*t1*t10*t31*t33*t4*t48
     +        -3.240503646201146d+0*gammaaa*t1*t27*t28*t33*t4
     +        -3.240503646201146d+0*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*t1*t16*t27*t4+7.561175174469341d+0
     +        *gammaaa*t1*t10*t21*t22*t28*t33-7.561175174469341d+0*t1
     +        *t10*t16*t21*t22)+1.5545345434827443d-2*t34)
     +        -1.5545345434827443d-2*rhoa*t32*t34*t49)*wght+Cmat2(iq,
     +        D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.5545345434827443d-2*rhoa*t20*(
     +        -2.1001727762085842d+1*gammaaa*t11*t13*t28*t33*t8
     +        +2.1001727762085842d+1*t11*t13*t16*t8
     +        +6.481007292402292d+0*gammaaa*t1*t10*t33**2*t4*t48
     +        -6.805617538956301d+1*gammaaa*t1**3*t10*t14*t28*t47
     +        -6.481007292402292d+0*t1*t10*t28*t33*t4)*wght
     +        -1.5545345434827443d-2*rhoa*t34**2*t49*wght+Cmat2(iq,
     +        D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhob
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call nwxcm_c_pw91lda_d2(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a,
     +        s1a2)
            t2 = s1f
            t3 = rhob**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhob
            t6 = exp(-6.432793688582966d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 3.240503646201146d+0*gammabb*t1*t4*t8
            t10 = t9+1.0d+0
            t11 = t1**2
            t12 = gammabb**2
            t13 = 1/t3**4.0d+0
            t14 = 1/t7**2
            t15 = t9+1.0500863881042921d+1*t11*t12*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 3.240503646201146d+0*gammabb*t1*t10*t16*t4+1.0d+0
            t18 = log(t17)
            t19 = s1a(D1_RB)
            t20 = 1/t17
            t21 = rhob**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -7.561175174469341d+0*gammabb*t1*t21*t22*t8
            t24 = 1/rhob**2
            t25 = 6.432793688582966d+1*s1f*t24-6.432793688582966d+1
     +        *t19*t5
            t26 = -3.240503646201146d+0*gammabb*t1*t14*t25*t4*t6
            t27 = t26+t23
            t28 = 1/t15**2
            t29 = 1/t3**5.0d+0
            t30 = 1/t7**3
            t31 = -2.1001727762085842d+1*t11*t12*t13*t25*t30*t6
     +        -4.900403144486697d+1*t11*t12*t14*t21*t29+t26+t23
            t32 = -3.240503646201146d+0*gammabb*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*gammabb*t1*t16*t27*t4
     +        -7.561175174469341d+0*gammabb*t1*t10*t16*t21*t22
            t33 = 3.240503646201146d+0*t1*t4*t8+2.1001727762085842d+1
     +        *gammabb*t11*t13*t14
            t34 = 1.0500863881042921d+1*gammabb*t11*t13*t16*t8
     +        -3.240503646201146d+0*gammabb*t1*t10*t28*t33*t4
     +        +3.240503646201146d+0*t1*t10*t16*t4
            t35 = s1a2(D2_RB_RB)
            t36 = rhob**3.333333333333335d-1
            t37 = 1/rhob**8.333333333333333d-1
            t38 = t25**2
            t39 = exp(-1.2865587377165932d+2*s1f*t5)
            t40 = 6.481007292402292d+0*gammabb*t1*t30*t38*t39*t4
            t41 = 2.6464113110642695d+1*gammabb*t1*t13*t36*t8
            t42 = -1.2601958624115575d+0*gammabb*t1*t22*t37*t8
            t43 = 1.5122350348938682d+1*gammabb*t1*t14*t21*t22*t25*t6
            t44 = -3.240503646201146d+0*gammabb*t1*t14*t38*t4*t6
            t45 = -6.432793688582966d+1*t35*t5+1.2865587377165932d+2
     +        *t19*t24-1.2865587377165932d+2*s1f/rhob**3
            t46 = -3.240503646201146d+0*gammabb*t1*t14*t4*t45*t6
            t47 = 1/t3**6.0d+0
            t48 = 1/t15**3
            t49 = 1/t17**2
            t50 = -7.561175174469341d+0*t1*t21*t22*t8
            t51 = -3.240503646201146d+0*t1*t14*t25*t4*t6
            fnc(iq) = 1.0d+0*(s1f+1.5545345434827443d-2*rhob*t18)*wght
     +        +fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*(1.5545345434827443d-2*rhob*t20
     +        *t32+t19+1.5545345434827443d-2*t18)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.5545345434827443d-2*rhob*t20*t34*wght
     +        +Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(t20*(1.5545345434827443d-2
     +        *rhob*(-3.240503646201146d+0*gammabb*t1*t10*t28*t4*
     +        (6.300518328625753d+1*t11*t12*t13*t38*t39/t7**4
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t45*t6
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t38*t6
     +        +1.9601612577946786d+2*t11*t12*t21*t25*t29*t30*t6
     +        +2.8585685009505735d+2*t11*t12*t14*t36*t47+t46+t44+t43
     +        +t42+t41+t40-8.167338574144498d+0*t11*t12*t14*t29*t37)
     +        +6.481007292402292d+0*gammabb*t1*t10*t31**2*t4*t48
     +        +3.240503646201146d+0*gammabb*t1*t16*t4*(t46+t44+t43+t42
     +        +t41+t40)-6.481007292402292d+0*gammabb*t1*t27*t28*t31*t4
     +        -1.2601958624115575d+0*gammabb*t1*t10*t16*t22*t37
     +        +2.6464113110642695d+1*gammabb*t1*t10*t13*t16*t36
     +        +1.5122350348938682d+1*gammabb*t1*t10*t21*t22*t28*t31
     +        -1.5122350348938682d+1*gammabb*t1*t16*t21*t22*t27)
     +        +3.1090690869654886d-2*t32)-1.5545345434827443d-2*rhob
     +        *t32**2*t49+t35)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(t20*(1.5545345434827443d-2
     +        *rhob*(-1.0500863881042921d+1*gammabb*t11*t13*t28*t31*t8
     +        -2.4502015722433484d+1*gammabb*t11*t16*t21*t29*t8
     +        -3.240503646201146d+0*gammabb*t1*t10*t28*t4*(
     +        -4.2003455524171684d+1*gammabb*t11*t13*t25*t30*t6+t51
     +        +t50-9.800806288973394d+1*gammabb*t11*t14*t21*t29)
     +        +3.240503646201146d+0*gammabb*t1*t16*t4*(t51+t50)
     +        +6.481007292402292d+0*gammabb*t1*t10*t31*t33*t4*t48
     +        -3.240503646201146d+0*gammabb*t1*t27*t28*t33*t4
     +        -3.240503646201146d+0*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*t1*t16*t27*t4+7.561175174469341d+0
     +        *gammabb*t1*t10*t21*t22*t28*t33-7.561175174469341d+0*t1
     +        *t10*t16*t21*t22)+1.5545345434827443d-2*t34)
     +        -1.5545345434827443d-2*rhob*t32*t34*t49)*wght+Cmat2(iq,
     +        D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5545345434827443d-2*rhob*t20*(
     +        -2.1001727762085842d+1*gammabb*t11*t13*t28*t33*t8
     +        +2.1001727762085842d+1*t11*t13*t16*t8
     +        +6.481007292402292d+0*gammabb*t1*t10*t33**2*t4*t48
     +        -6.805617538956301d+1*gammabb*t1**3*t10*t14*t28*t47
     +        -6.481007292402292d+0*t1*t10*t28*t33*t4)*wght
     +        -1.5545345434827443d-2*rhob*t34**2*t49*wght+Cmat2(iq,
     +        D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_pbe_d2 
C>
C> \brief Evaluate the nwxcm_c_pbe functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\
C>   {\it t_2} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_\alpha , \rho_\beta\right)\\\\
C>   {\it t_3} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_4} &=& \rho_\alpha-\rho_\beta\\\\
C>   {\it t_5} &=& {{1}\over{{\it t_3}}}\\\\
C>   {\it t_6} &=& 0.5\,\left({\it t_4}\,{\it t_5}
C>    +1.0\right)^{{{2}\over{3}}}+0.5\,\left(1.0-{\it t_4}\,{
C>    \it t_5}\right)^{{{2}\over{3}}}\\\\
C>   {\it t_7} &=& {\it t_6}^3\\\\
C>   {\it t_8} &=& {\it t_3}^{1.166666666666667}\\\\
C>   {\it t_9} &=& {{1}\over{{\it t_8}^{2.0}}}\\\\
C>   {\it t_{10}} &=& {{1}\over{{\it t_6}^{2.0}}}\\\\
C>   {\it t_{11}} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\
C>   {\it t_{12}} &=& {{1}\over{e^{{{32.16396844291483\,{\it t_2}
C>    \,{\it t_5}}\over{{\it t_7}}}}}}-1.0\\\\
C>   {\it t_{13}} &=& {{2.041389378054956\,{\it t_1}\,{\it t_9}\,{
C>    \it t_{10}}\,{\it t_{11}}}\over{{\it t_{12}}}}\\\\
C>   {\it t_{14}} &=& {\it t_1}^2\\\\
C>   {\it t_{15}} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_s , 0.0\right)\\\\
C>   {\it t_{16}} &=& \rho_s^{1.166666666666667}\\\\
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}^{2.0}}}\\\\
C>   {\it t_{18}} &=& {{1}\over{e^{{{64.32793688582966\,{
C>    \it t_{15}}}\over{\rho_s}}}}}-1.0\\\\
C>   {\it t_{19}} &=& {{3.240503646201146\,{\it t_1}\,{\it t_{17}}
C>    \,\sigma_{ss}}\over{{\it t_{18}}}}\\\\
C>   f &=& 1.0\,\left(0.03109069086965489\,{\it t_3}\,{\it t_7}
C>    \,\log \left({{2.041389378054956\,{\it t_1}\,{\it t_9}\,{
C>    \it t_{10}}\,{\it t_{11}}\,\left({\it t_{13}}
C>    +1.0\right)}\over{{{4.1672705928356\,{\it t_{14}}\,{
C>    \it t_{11}}^2}\over{{\it t_8}^{4.0}\,{\it t_6}^{4.0}\,{
C>    \it t_{12}}^2}}+{\it t_{13}}+1.0}}+1.0\right)+{
C>    \it t_2}\right)\\\\
C>   g &=& 0\\\\
C>   G &=& 1.0\,\left(0.01554534543482744\,\rho_s\,\log 
C>    \left({{3.240503646201146\,{\it t_1}\,{\it t_{17}}\,
C>    \sigma_{ss}\,\left({\it t_{19}}+1.0\right)}
C>    \over{{{10.50086388104292\,{\it t_{14}}\,\sigma_{ss}^2}
C>    \over{{\it t_{16}}^{4.0}\,{\it t_{18}}^2}}+{\it t_{19}}+1.0}}
C>    +1.0\right)+{\it t_{15}}\right)\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc [4,5,6].
C>
C> ### References ###
C>
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77
C>    , 3865 (1996)  , DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865 ">
C> 10.1103/PhysRevLett.77.3865 </a>
C>
C> [2] Maxima, a computer algebra system,
C> <a href="http://maxima.sourceforge.net/">
C> http://maxima.sourceforge.net/</a>
C>
C> [3] CLISP 2.49 (2010-07-07)
C>
C> [4] autoxc, revision 27694 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_pbe_d3(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,fnc,Amat,Amat2,Amat3, Cmat,Cmat2,Cmat3) 
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
      double precision fnc(nq
     +   )      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,NCOL_AMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt rho 
      double precision Cmat3(nq,NCOL_CMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
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
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            sr(R_T) = 2.0d0*rhoa
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
            call nwxcm_c_pw91lda_d3(t1,tol_rho,ipol,1,1.0d0,sr,s1f,
     +        s1a,s1a2,s1a3)
            t2 = s1f
            t3 = rhoa**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhoa
            t6 = exp(-1.6081984221457415d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 1.6202518231005725d+0*gammaaa*t1*t4*t8
            t10 = t1**2
            t11 = gammaaa**2
            t12 = 1/t3**4.0d+0
            t13 = 1/t7**2
            t14 = s1a(D1_RA)
            t15 = 1.6202518231005728d+0*gammaaa*t1*t4*t8
            t16 = t15+1.0d+0
            t17 = t15+2.625215970260729d+0*t10*t11*t12*t13+1.0d+0
            t18 = 1/t17
            t19 = 1.6202518231005728d+0*gammaaa*t1*t16*t18*t4+1.0d+0
            t20 = 1/t19
            t21 = rhoa**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -1.890293793617335d+0*gammaaa*t1*t21*t22*t8
            t24 = 1/rhoa**2
            t25 = 8.040992110728707d+0*s1f*t24-1.6081984221457415d+1
     +        *t14*t5
            t26 = -1.6202518231005728d+0*gammaaa*t1*t13*t25*t4*t6
            t27 = t26+t23
            t28 = 1/t17**2
            t29 = 1/t3**5.0d+0
            t30 = 1/t7**3
            t31 = -5.250431940521458d+0*t10*t11*t12*t25*t30*t6
     +        -6.125503930608368d+0*t10*t11*t13*t21*t29+t26+t23
            t32 = -1.6202518231005728d+0*gammaaa*t1*t16*t28*t31*t4
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t27*t4
     +        -1.890293793617335d+0*gammaaa*t1*t16*t18*t21*t22
            t33 = log(t19)
            t34 = 4.050629557751432d-1*t1*t4*t8+1.3126079851303646d+0
     +        *gammaaa*t10*t12*t13
            t35 = 6.563039925651822d-1*gammaaa*t10*t12*t18*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t34*t4
     +        +4.050629557751432d-1*t1*t16*t18*t4
            t36 = 8.101259115502863d-1*t1*t4*t8+2.625215970260729d+0
     +        *gammaaa*t10*t12*t13
            t37 = 1.3126079851303646d+0*gammaaa*t10*t12*t18*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t36*t4
     +        +8.101259115502863d-1*t1*t16*t18*t4
            t38 = s1a2(D2_RA_RA)
            t39 = 6.218138173930978d-2*t32
            t40 = rhoa**3.333333333333335d-1
            t41 = 3.3080141388303363d+0*gammaaa*t1*t12*t16*t18*t40
            t42 = 1/rhoa**8.333333333333333d-1
            t43 = -1.5752448280144465d-1*gammaaa*t1*t16*t18*t22*t42
            t44 = -3.78058758723467d+0*gammaaa*t1*t18*t21*t22*t27
            t45 = 3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t28*t31
            t46 = -3.2405036462011455d+0*gammaaa*t1*t27*t28*t31*t4
            t47 = t25**2
            t48 = exp(-3.216396844291483d+1*s1f*t5)
            t49 = 3.2405036462011455d+0*gammaaa*t1*t30*t4*t47*t48
            t50 = 3.3080141388303363d+0*gammaaa*t1*t12*t40*t8
            t51 = -1.5752448280144465d-1*gammaaa*t1*t22*t42*t8
            t52 = 1.80027980344508d-1*gammaaa*t1*t24*t4*t8
            t53 = 3.78058758723467d+0*gammaaa*t1*t13*t21*t22*t25*t6
            t54 = -1.6202518231005728d+0*gammaaa*t1*t13*t4*t47*t6
            t55 = 1/rhoa**3
            t56 = 1.6081984221457415d+1*t14*t24
            t57 = t56-1.0721322814304943d+1*s1f*t55
     +        -1.6081984221457415d+1*t38*t5
            t58 = -1.6202518231005728d+0*gammaaa*t1*t13*t4*t57*t6
            t59 = t58+t54+t53+t52+t51+t50+t49
            t60 = 1/t7**4
            t61 = 1.5751295821564373d+1*t10*t11*t12*t47*t48*t60
            t62 = 1/t3**6.0d+0
            t63 = 1.7866053130941079d+1*t10*t11*t13*t40*t62
            t64 = -5.10458660884031d-1*t10*t11*t13*t29*t42
            t65 = 2.4502015722433473d+1*t10*t11*t21*t25*t29*t30*t6
            t66 = -5.250431940521458d+0*t10*t11*t12*t30*t47*t6
            t67 = t66+t65+t64+t63+t61-5.250431940521458d+0*t10*t11*t12
     +        *t30*t57*t6+t58+t54+t53+t52+t51+t50+t49
     +        +5.833813267246063d-1*t10*t11*t12*t13*t24
            t68 = 1/t17**3
            t69 = t31**2
            t70 = 3.2405036462011455d+0*gammaaa*t1*t16*t4*t68*t69
            t71 = t70-1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*t67
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*t59+t46+t45+t44
     +        +t43+t41+1.80027980344508d-1*gammaaa*t1*t16*t18*t24*t4
            t72 = 6.218138173930978d-2*rhoa*t71+t39
            t73 = 1/t19**2
            t74 = t32**2
            t75 = -6.218138173930978d-2*rhoa*t73*t74
            t76 = s1a2(D2_RA_RB)
            t77 = -1.8002798034450807d-1*gammaaa*t1*t24*t4*t8
            t78 = -1.6081984221457415d+1*t5*t76+t56
     +        -5.360661407152472d+0*s1f*t55
            t79 = -1.6202518231005728d+0*gammaaa*t1*t13*t4*t6*t78
            t80 = t79+t77+t54+t53+t51+t50+t49
            t81 = t79-5.250431940521458d+0*t10*t11*t12*t30*t6*t78+t77
     +        +t66+t65+t64+t63+t61+t54+t53+t51+t50+t49
     +        -5.833813267246064d-1*t10*t11*t12*t13*t24
            t82 = -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*t81
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*t80+t70+t46+t45
     +        +t44+t43+t41-1.8002798034450807d-1*gammaaa*t1*t16*t18
     +        *t24*t4
            t83 = 6.218138173930978d-2*rhoa*t82+t39
            t84 = -4.7257344840433374d-1*t1*t21*t22*t8
            t85 = -4.050629557751432d-1*t1*t13*t25*t4*t6
            t86 = t85+t84
            t87 = t85+t84-2.625215970260729d+0*gammaaa*t10*t12*t25*t30
     +        *t6-3.062751965304184d+0*gammaaa*t10*t13*t21*t29
            t88 = -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*t87
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*t86
     +        -6.563039925651822d-1*gammaaa*t10*t12*t28*t31*t8
     +        -7.65687991326046d-1*gammaaa*t10*t18*t21*t29*t8
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t31*t34*t4*t68
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t28*t34*t4
     +        -4.050629557751432d-1*t1*t16*t28*t31*t4
     +        +4.050629557751432d-1*t1*t18*t27*t4+1.890293793617335d+0
     +        *gammaaa*t1*t16*t21*t22*t28*t34-4.7257344840433374d-1*t1
     +        *t16*t18*t21*t22
            t89 = 6.218138173930978d-2*rhoa*t88+3.109069086965489d-2*t35
            t90 = 1.0d+0*(t20*t89-6.218138173930978d-2*rhoa*t32*t35
     +        *t73)*wght
            t91 = -9.451468968086675d-1*t1*t21*t22*t8
            t92 = -8.101259115502863d-1*t1*t13*t25*t4*t6
            t93 = t92+t91
            t94 = t92+t91-5.250431940521458d+0*gammaaa*t10*t12*t25*t30
     +        *t6-6.125503930608368d+0*gammaaa*t10*t13*t21*t29
            t95 = -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*t94
     +        +1.6202518231005728d+0*gammaaa*t1*t18*t4*t93
     +        -1.3126079851303646d+0*gammaaa*t10*t12*t28*t31*t8
     +        -1.531375982652092d+0*gammaaa*t10*t18*t21*t29*t8
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t31*t36*t4*t68
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t28*t36*t4
     +        -8.101259115502863d-1*t1*t16*t28*t31*t4
     +        +8.101259115502863d-1*t1*t18*t27*t4+1.890293793617335d+0
     +        *gammaaa*t1*t16*t21*t22*t28*t36-9.451468968086675d-1*t1
     +        *t16*t18*t21*t22
            t96 = 6.218138173930978d-2*rhoa*t95+3.109069086965489d-2*t37
            t97 = t35**2
            t98 = t34**2
            t99 = -8.101259115502863d-1*t1*t16*t28*t34*t4
            t100 = t1**3
            t101 = -1.3126079851303646d+0*gammaaa*t10*t12*t28*t34*t8
            t102 = t99+3.2405036462011455d+0*gammaaa*t1*t16*t4*t68*t98
     +        +3.281519962825911d-1*t10*t12*t18*t8
     +        -5.316888702309608d-1*gammaaa*t100*t13*t16*t28*t62+t101
            t103 = 6.218138173930978d-2*rhoa*t102*t20*wght
     +        -6.218138173930978d-2*rhoa*t73*t97*wght
            t104 = t99-6.563039925651822d-1*gammaaa*t10*t12*t28*t36*t8
     +        +6.563039925651822d-1*t10*t12*t18*t8
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t34*t36*t4*t68
     +        -1.0633777404619216d+0*gammaaa*t100*t13*t16*t28*t62
     +        -4.050629557751432d-1*t1*t16*t28*t36*t4+t101
            t105 = t37**2
            t106 = t36**2
            t107 = -2.625215970260729d+0*gammaaa*t10*t12*t28*t36*t8
     +        +1.3126079851303646d+0*t10*t12*t18*t8
     +        +3.2405036462011455d+0*gammaaa*t1*t106*t16*t4*t68
     +        -2.1267554809238431d+0*gammaaa*t100*t13*t16*t28*t62
     +        -1.6202518231005728d+0*t1*t16*t28*t36*t4
            t108 = s1a3(D3_RA_RA_RA)
            t109 = rhoa**5.000000000000002d-1
            t110 = -7.718699657270787d+0*gammaaa*t1*t109*t16*t18*t29
            t111 = 1/rhoa**6.666666666666665d-1
            t112 = 8.270035347075845d-1*gammaaa*t1*t111*t12*t16*t18
            t113 = 1/rhoa**1.8333333333333332d+0
            t114 = 9.924042416491009d+0*gammaaa*t1*t12*t18*t27*t40
            t115 = -4.7257344840433396d-1*gammaaa*t1*t18*t22*t27*t42
            t116 = -9.924042416491009d+0*gammaaa*t1*t12*t16*t28*t31*t40
            t117 = 4.7257344840433396d-1*gammaaa*t1*t16*t22*t28*t31*t42
            t118 = 1.1341762761704009d+1*gammaaa*t1*t21*t22*t27*t28*t31
            t119 = t25**3
            t120 = exp(-4.824595266437224d+1*s1f*t5)
            t121 = -9.721510938603435d+0*gammaaa*t1*t119*t120*t4*t60
            t122 = -1.1341762761704009d+1*gammaaa*t1*t21*t22*t30*t47*t48
            t123 = 1.6081984221457415d+1*s1f*t24-3.216396844291483d+1
     +        *t14*t5
            t124 = 3.2405036462011455d+0*gammaaa*t1*t123*t30*t4*t47*t48
            t125 = 3.2405036462011455d+0*gammaaa*t1*t119*t30*t4*t48
            t126 = 9.721510938603435d+0*gammaaa*t1*t25*t30*t4*t48*t57
            t127 = -7.718699657270787d+0*gammaaa*t1*t109*t29*t8
            t128 = 8.270035347075845d-1*gammaaa*t1*t111*t12*t8
            t129 = -5.644627300385094d-1*gammaaa*t1*t113*t22*t8
            t130 = -5.400839410335243d-1*gammaaa*t1*t4*t55*t8
            t131 = -9.924042416491009d+0*gammaaa*t1*t12*t13*t25*t40*t6
            t132 = 4.7257344840433396d-1*gammaaa*t1*t13*t22*t25*t42*t6
            t133 = 5.670881380852004d+0*gammaaa*t1*t13*t21*t22*t47*t6
            t134 = 5.670881380852004d+0*gammaaa*t1*t13*t21*t22*t57*t6
            t135 = -1.6202518231005728d+0*gammaaa*t1*t119*t13*t4*t6
            t136 = -5.40083941033524d-1*gammaaa*t1*t13*t24*t25*t4*t6
            t137 = -4.8607554693017174d+0*gammaaa*t1*t13*t25*t4*t57*t6
            t138 = 1/rhoa**4
            t139 = -3.216396844291482d+1*t14*t55-1.6081984221457415d+1
     +        *t108*t5+2.412297633218612d+1*t24*t38
     +        +2.412297633218612d+1*t138*s1f
            t140 = -1.6202518231005728d+0*gammaaa*t1*t13*t139*t4*t6
            t141 = -6.300518328625749d+1*t10*t11*t119*t12*t120/t7**5
            t142 = -1.1025907075095062d+2*t10*t11*t21*t29*t47*t48*t60
            t143 = 1.5751295821564373d+1*t10*t11*t12*t123*t47*t48*t60
            t144 = 1.5751295821564373d+1*t10*t11*t119*t12*t48*t60
            t145 = 1/t3**7.0d+0
            t146 = -6.253118595829379d+1*t10*t109*t11*t13*t145
            t147 = 4.4665132827352727d+0*t10*t11*t111*t13*t62
            t148 = 3.6753023583650213d+1*t10*t11*t21*t29*t30*t47*t6
            t149 = -5.250431940521458d+0*t10*t11*t119*t12*t30*t6
            t150 = -1.1341762761704009d+1*gammaaa*t1*t16*t21*t22*t68*t69
            t151 = 1/t17**4
            t152 = -9.721510938603435d+0*gammaaa*t1*t151*t16*t31**3*t4
            t153 = -t32*t72
            t154 = -3.109069086965489d-2*t74
            t155 = 1/t19**3
            t156 = 1.2436276347861956d-1*rhoa*t155*t32**3
            t157 = s1a3(D3_RA_RA_RB)
            t158 = 6.481007292402291d+0*gammaaa*t1*t25*t30*t4*t48*t78
            t159 = 3.2405036462011455d+0*gammaaa*t1*t25*t30*t4*t48*t57
            t160 = 2.756678449025281d-1*gammaaa*t1*t113*t22*t8
            t161 = 1.8002798034450812d-1*gammaaa*t1*t4*t55*t8
            t162 = 3.78058758723467d+0*gammaaa*t1*t13*t21*t22*t6*t78
            t163 = 1.890293793617335d+0*gammaaa*t1*t13*t21*t22*t57*t6
            t164 = 1.8002798034450812d-1*gammaaa*t1*t13*t24*t25*t4*t6
            t165 = -3.2405036462011455d+0*gammaaa*t1*t13*t25*t4*t6*t78
            t166 = -1.6202518231005728d+0*gammaaa*t1*t13*t25*t4*t57*t6
            t167 = 1.6081984221457415d+1*t24*t76-2.1442645628609888d+1
     +        *t14*t55-1.6081984221457415d+1*t157*t5
     +        +8.040992110728707d+0*t24*t38+8.040992110728705d+0*t138
     +        *s1f
            t168 = -1.6202518231005728d+0*gammaaa*t1*t13*t167*t4*t6
            t169 = 6.218138173930978d-2*t88
            t170 = -3.3080141388303363d+0*gammaaa*t1*t12*t16*t28*t34*t40
            t171 = 1.5752448280144465d-1*gammaaa*t1*t16*t22*t28*t34*t42
            t172 = 8.27003534707584d-1*t1*t12*t16*t18*t40
            t173 = -3.938112070036116d-2*t1*t16*t18*t22*t42
            t174 = 1.339953984820581d+0*gammaaa*t10*t18*t40*t62*t8
            t175 = -6.380733261050388d-2*gammaaa*t10*t18*t29*t42*t8
            t176 = -3.78058758723467d+0*gammaaa*t1*t18*t21*t22*t86
            t177 = 3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t28*t87
            t178 = 3.78058758723467d+0*gammaaa*t1*t21*t22*t27*t28*t34
            t179 = -9.451468968086675d-1*t1*t18*t21*t22*t27
            t180 = -3.2405036462011455d+0*gammaaa*t1*t27*t28*t4*t87
            t181 = -7.56117517446934d+0*gammaaa*t1*t16*t21*t22*t31*t34
     +        *t68
            t182 = 9.451468968086675d-1*t1*t16*t21*t22*t28*t31
            t183 = 1.531375982652092d+0*gammaaa*t10*t21*t28*t29*t31*t8
            t184 = -3.2405036462011455d+0*gammaaa*t1*t28*t31*t4*t86
            t185 = 6.481007292402291d+0*gammaaa*t1*t16*t31*t4*t68*t87
            t186 = 6.481007292402291d+0*gammaaa*t1*t27*t31*t34*t4*t68
            t187 = -8.101259115502863d-1*t1*t27*t28*t31*t4
            t188 = 8.101259115502863d-1*t1*t30*t4*t47*t48
            t189 = 8.27003534707584d-1*t1*t12*t40*t8
            t190 = -3.938112070036116d-2*t1*t22*t42*t8
            t191 = 4.5006995086127d-2*t1*t24*t4*t8
            t192 = 9.451468968086675d-1*t1*t13*t21*t22*t25*t6
            t193 = -4.050629557751432d-1*t1*t13*t4*t47*t6
            t194 = -4.050629557751432d-1*t1*t13*t4*t57*t6
            t195 = 7.875647910782188d+0*gammaaa*t10*t12*t47*t48*t60
            t196 = 8.93302656547054d+0*gammaaa*t10*t13*t40*t62
            t197 = -2.552293304420155d-1*gammaaa*t10*t13*t29*t42
            t198 = 1.2251007861216737d+1*gammaaa*t10*t21*t25*t29*t30*t6
            t199 = -2.625215970260729d+0*gammaaa*t10*t12*t30*t47*t6
            t200 = -9.721510938603435d+0*gammaaa*t1*t151*t16*t34*t4*t69
            t201 = 8.101259115502863d-1*t1*t16*t4*t68*t69
            t202 = 1.3126079851303646d+0*gammaaa*t10*t12*t68*t69*t8
            t203 = -1.2436276347861956d-1*rhoa*t32*t88
            t204 = 1.2436276347861956d-1*rhoa*t155*t35*t74
            t205 = 1.0d+0*(t20*(6.218138173930978d-2*rhoa*(
     +        -6.563039925651822d-1*gammaaa*t10*t12*t28*t67*t8
     +        +7.292266584057578d-2*gammaaa*t10*t12*t18*t24*t8
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t34*t4*t67*t68
     +        -4.050629557751432d-1*t1*t16*t28*t4*t67
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -2.625215970260729d+0*gammaaa*t10*t12*t30*t57*t6
     +        +2.9169066336230315d-1*gammaaa*t10*t12*t13*t24+t199+t198
     +        +t197+t196+t195+t194+t193+t192+t191+t190+t189+t188)
     +        -1.6202518231005728d+0*gammaaa*t1*t28*t34*t4*t59
     +        +4.050629557751432d-1*t1*t18*t4*t59-1.80027980344508d-1
     +        *gammaaa*t1*t16*t24*t28*t34*t4+4.5006995086127d-2*t1*t16
     +        *t18*t24*t4+1.6202518231005728d+0*gammaaa*t1*t18*(t194
     +        +t193+t192+t191+t190+t189+t188)*t4+t202+t201+t200+t187
     +        +t186+t185+t184+t183+t182+t181+t180+t179+t178+t177+t176
     +        +t175+t174+t173+t172+t171+t170)-1.0363563623218294d-2
     +        *t35*t5+t169)+(t203-t35*t72)*t73+t204)*wght
            t206 = 6.218138173930978d-2*t95
            t207 = -3.3080141388303363d+0*gammaaa*t1*t12*t16*t28*t36*t40
            t208 = 1.5752448280144465d-1*gammaaa*t1*t16*t22*t28*t36*t42
            t209 = 1.6540070694151682d+0*t1*t12*t16*t18*t40
            t210 = -7.876224140072232d-2*t1*t16*t18*t22*t42
            t211 = 2.679907969641162d+0*gammaaa*t10*t18*t40*t62*t8
            t212 = -1.2761466522100776d-1*gammaaa*t10*t18*t29*t42*t8
            t213 = -3.78058758723467d+0*gammaaa*t1*t18*t21*t22*t93
            t214 = 3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t28*t94
            t215 = 3.78058758723467d+0*gammaaa*t1*t21*t22*t27*t28*t36
            t216 = -1.890293793617335d+0*t1*t18*t21*t22*t27
            t217 = -3.2405036462011455d+0*gammaaa*t1*t27*t28*t4*t94
            t218 = -7.56117517446934d+0*gammaaa*t1*t16*t21*t22*t31*t36
     +        *t68
            t219 = 1.890293793617335d+0*t1*t16*t21*t22*t28*t31
            t220 = 3.062751965304184d+0*gammaaa*t10*t21*t28*t29*t31*t8
            t221 = -3.2405036462011455d+0*gammaaa*t1*t28*t31*t4*t93
            t222 = 6.481007292402291d+0*gammaaa*t1*t16*t31*t4*t68*t94
            t223 = 6.481007292402291d+0*gammaaa*t1*t27*t31*t36*t4*t68
            t224 = -1.6202518231005728d+0*t1*t27*t28*t31*t4
            t225 = 1.6202518231005728d+0*t1*t30*t4*t47*t48
            t226 = 1.6540070694151682d+0*t1*t12*t40*t8
            t227 = -7.876224140072232d-2*t1*t22*t42*t8
            t228 = 9.0013990172254d-2*t1*t24*t4*t8
            t229 = 1.890293793617335d+0*t1*t13*t21*t22*t25*t6
            t230 = -8.101259115502863d-1*t1*t13*t4*t47*t6
            t231 = -8.101259115502863d-1*t1*t13*t4*t57*t6
            t232 = 1.5751295821564373d+1*gammaaa*t10*t12*t47*t48*t60
            t233 = 1.7866053130941079d+1*gammaaa*t10*t13*t40*t62
            t234 = -5.10458660884031d-1*gammaaa*t10*t13*t29*t42
            t235 = 2.4502015722433473d+1*gammaaa*t10*t21*t25*t29*t30*t6
            t236 = -5.250431940521458d+0*gammaaa*t10*t12*t30*t47*t6
            t237 = -9.721510938603435d+0*gammaaa*t1*t151*t16*t36*t4*t69
            t238 = 1.6202518231005728d+0*t1*t16*t4*t68*t69
            t239 = 2.625215970260729d+0*gammaaa*t10*t12*t68*t69*t8
            t240 = -1.2436276347861956d-1*rhoa*t32*t95
            t241 = 1.2436276347861956d-1*rhoa*t155*t37*t74
            t242 = -4.5006995086127016d-2*t1*t24*t4*t8
            t243 = -4.050629557751432d-1*t1*t13*t4*t6*t78
            t244 = -9.001399017225403d-2*t1*t24*t4*t8
            t245 = -8.101259115502863d-1*t1*t13*t4*t6*t78
            t246 = 9.451468968086675d-1*t1*t16*t21*t22*t28*t34
            t247 = 1.531375982652092d+0*gammaaa*t10*t21*t28*t29*t34*t8
            t248 = 8.101259115502863d-1*t1*t18*t4*t86
            t249 = -8.101259115502863d-1*t1*t16*t28*t4*t87
            t250 = -1.3126079851303646d+0*gammaaa*t10*t12*t28*t8*t87
            t251 = -8.101259115502863d-1*t1*t27*t28*t34*t4
            t252 = 1.6202518231005728d+0*t1*t16*t31*t34*t4*t68
            t253 = 2.625215970260729d+0*gammaaa*t10*t12*t31*t34*t68*t8
            t254 = 1.0d+0*(t20*(6.218138173930978d-2*rhoa*
     +        (3.2405036462011455d+0*gammaaa*t1*t27*t4*t68*t98
     +        -3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t68*t98
     +        -9.721510938603435d+0*gammaaa*t1*t151*t16*t31*t4*t98
     +        +6.481007292402291d+0*gammaaa*t1*t16*t34*t4*t68*t87
     +        -3.2405036462011455d+0*gammaaa*t1*t28*t34*t4*t86
     +        -3.281519962825911d-1*t10*t12*t28*t31*t8
     +        -3.82843995663023d-1*t10*t18*t21*t29*t8
     +        +1.0633777404619216d+0*gammaaa*t100*t13*t16*t31*t62*t68
     +        -5.316888702309608d-1*gammaaa*t100*t13*t27*t28*t62
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -6.563039925651822d-1*t10*t12*t25*t30*t6
     +        -7.65687991326046d-1*t10*t13*t21*t29)
     +        +6.203036819361207d-1*gammaaa*t100*t13*t145*t16*t21*t28
     +        +t253+t252+t251+t250+t249+t248+t247+t246)
     +        +3.109069086965489d-2*t102)+1.2436276347861956d-1*rhoa
     +        *t155*t32*t97+t73*(-t35*t89-6.218138173930978d-2*rhoa
     +        *t35*t88-6.218138173930978d-2*rhoa*t102*t32))*wght
            t255 = 1.2436276347861956d-1*rhoa*t155*t32*t35*t37
            t256 = t20*(6.218138173930978d-2*rhoa*(
     +        -6.563039925651822d-1*gammaaa*t10*t12*t28*t8*t94
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t34*t4*t68*t94
     +        -4.050629557751432d-1*t1*t16*t28*t4*t94
     +        -1.6202518231005728d+0*gammaaa*t1*t28*t34*t4*t93
     +        +4.050629557751432d-1*t1*t18*t4*t93
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t36*t4*t68*t87
     +        -1.6202518231005728d+0*gammaaa*t1*t28*t36*t4*t86
     +        +1.3126079851303646d+0*gammaaa*t10*t12*t31*t36*t68*t8
     +        +7.65687991326046d-1*gammaaa*t10*t21*t28*t29*t36*t8
     +        -6.563039925651822d-1*t10*t12*t28*t31*t8
     +        -7.65687991326046d-1*t10*t18*t21*t29*t8
     +        +2.1267554809238431d+0*gammaaa*t100*t13*t16*t31*t62*t68
     +        +3.2405036462011455d+0*gammaaa*t1*t27*t34*t36*t4*t68
     +        +8.101259115502863d-1*t1*t16*t31*t36*t4*t68
     +        -3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t34*t36*t68
     +        -1.0633777404619216d+0*gammaaa*t100*t13*t27*t28*t62
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -1.3126079851303646d+0*t10*t12*t25*t30*t6
     +        -1.531375982652092d+0*t10*t13*t21*t29)
     +        -9.721510938603435d+0*gammaaa*t1*t151*t16*t31*t34*t36*t4
     +        -4.050629557751432d-1*t1*t27*t28*t36*t4
     +        +4.7257344840433374d-1*t1*t16*t21*t22*t28*t36
     +        +1.2406073638722415d+0*gammaaa*t100*t13*t145*t16*t21*t28
     +        +t253+t252+t251+t250+t249+t248+t247+t246)
     +        +3.109069086965489d-2*t104)
            t257 = -6.218138173930978d-2*rhoa*t104*t32
            t258 = t1**4
            t259 = 1/t3**8.0d+0
            t260 = 6.218138173930978d-2*rhoa*t20*(3.937823955391094d+0
     +        *gammaaa*t10*t12*t68*t8*t98+2.4303777346508587d+0*t1*t16
     +        *t4*t68*t98-9.844559888477734d-1*t10*t12*t28*t34*t8
     +        +3.1901332213857647d+0*gammaaa*t100*t13*t16*t34*t62*t68
     +        -3.987666526732206d-1*t100*t13*t16*t28*t62
     +        -9.721510938603435d+0*gammaaa*t1*t151*t16*t34**3*t4
     +        -6.461023959854985d-1*gammaaa*t258*t259*t28*t30)*wght
     +        -1.8654414521792934d-1*rhoa*t102*t35*t73*wght
     +        +1.2436276347861956d-1*rhoa*t155*t35**3*wght
            t261 = -9.721510938603435d+0*gammaaa*t1*t151*t16*t36*t4*t98
            t262 = 1.6202518231005728d+0*t1*t16*t4*t68*t98
            t263 = 1.6202518231005728d+0*t1*t16*t34*t36*t4*t68
            t264 = 4.2535109618476863d+0*gammaaa*t100*t13*t16*t34*t62
     +        *t68
            t265 = 1.0633777404619216d+0*gammaaa*t100*t13*t16*t36*t62
     +        *t68
            t266 = -7.975333053464412d-1*t100*t13*t16*t28*t62
            t267 = 2.625215970260729d+0*gammaaa*t10*t12*t68*t8*t98
            t268 = 2.625215970260729d+0*gammaaa*t10*t12*t34*t36*t68*t8
            t269 = -1.3126079851303646d+0*t10*t12*t28*t34*t8
            t270 = -3.281519962825911d-1*t10*t12*t28*t36*t8
            t271 = 1.2436276347861956d-1*rhoa*t155*t37*t97*wght
            t272 = t73*(-6.218138173930978d-2*rhoa*t102*t37*wght
     +        -1.2436276347861956d-1*rhoa*t104*t35*wght)
            fnc(iq) = 1.0d+0*(6.218138173930978d-2*rhoa*log
     +        (1.6202518231005725d+0*gammaaa*t1*t4*(t9+1.0d+0)/(t9
     +        +2.6252159702607286d+0*t10*t11*t12*t13+1.0d+0)+1.0d+0)
     +        +s1f)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(3.109069086965489d-2*t33
     +        +6.218138173930978d-2*rhoa*t20*t32+t14)*wght+Amat(iq,
     +        D1_RA)
            Cmat(iq,D1_GAA) = 6.218138173930978d-2*rhoa*t20*t35*wght
     +        +Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 6.218138173930978d-2*rhoa*t20*t37*wght
     +        +Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t75+t20*t72
     +        -1.0363563623218294d-2*t33*t5+t38)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t20*t83+t76+t75
     +        +1.0363563623218296d-2*t33*t5)*wght+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t90+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(t20*t96-6.218138173930978d-2
     +        *rhoa*t32*t37*t73)*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t90+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t103+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = -6.218138173930978d-2*rhoa*t35*t37
     +        *t73*wght+6.218138173930978d-2*rhoa*t104*t20*wght+Cmat2
     +        (iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t103+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = -6.218138173930978d-2*rhoa*t105*t73
     +        *wght+6.218138173930978d-2*rhoa*t107*t20*wght+Cmat2(iq,
     +        D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*((-1.2436276347861956d-1
     +        *rhoa*t32*t71+t154+t153)*t73+t20*(9.327207260896466d-2
     +        *t71+6.218138173930978d-2*rhoa*(9.721510938603435d+0
     +        *gammaaa*t1*t27*t4*t68*t69+9.721510938603435d+0*gammaaa
     +        *t1*t16*t31*t4*t67*t68-4.8607554693017174d+0*gammaaa*t1
     +        *t27*t28*t4*t67+5.670881380852004d+0*gammaaa*t1*t16*t21
     +        *t22*t28*t67-1.6202518231005728d+0*gammaaa*t1*t16*t28*t4
     +        *(-1.0719631878564648d+2*t10*t11*t25*t30*t40*t6*t62
     +        +4.7253887464693123d+1*t10*t11*t12*t25*t48*t57*t60
     +        +3.6753023583650213d+1*t10*t11*t21*t29*t30*t57*t6
     +        -1.5751295821564373d+1*t10*t11*t12*t25*t30*t57*t6
     +        +3.0627519653041857d+0*t10*t11*t25*t29*t30*t42*t6
     +        -3.500287960347638d+0*t10*t11*t12*t24*t25*t30*t6
     +        -5.250431940521458d+0*t10*t11*t12*t139*t30*t6
     +        -1.7501439801738194d+0*t10*t11*t12*t13*t55
     +        -3.870978178370565d+0*t10*t11*t113*t13*t29+t149+t148
     +        +t147+t146+t144+t143+t142+t141+t140+t137+t136+t135+t134
     +        +t133+t132+t131+t130+t129+t128+t127+t126+t125+t124+t122
     +        +t121)-4.8607554693017174d+0*gammaaa*t1*t28*t31*t4*t59
     +        -5.670881380852004d+0*gammaaa*t1*t18*t21*t22*t59
     +        -5.400839410335243d-1*gammaaa*t1*t16*t18*t4*t55
     +        -5.40083941033524d-1*gammaaa*t1*t16*t24*t28*t31*t4
     +        +5.40083941033524d-1*gammaaa*t1*t18*t24*t27*t4
     +        +1.6202518231005728d+0*gammaaa*t1*(t140+t137+t136+t135
     +        +t134+t133+t132+t131+t130+t129+t128+t127+t126+t125+t124
     +        +t122+t121)*t18*t4-5.644627300385094d-1*gammaaa*t1*t113
     +        *t16*t18*t22+t152+t150+t118+t117+t116+t115+t114+t112
     +        +t110)-3.1090690869654886d-2*t32*t5)
     +        +1.554534543482745d-2*t24*t33+t156+t108)*wght+Amat3(iq,
     +        D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 1.0d+0*(t73*(
     +        -1.2436276347861956d-1*rhoa*t32*t82+t154+t153)+t20*
     +        (6.218138173930978d-2*t82+6.218138173930978d-2*rhoa*
     +        (6.481007292402291d+0*gammaaa*t1*t16*t31*t4*t68*t81
     +        -3.2405036462011455d+0*gammaaa*t1*t27*t28*t4*t81
     +        +3.78058758723467d+0*gammaaa*t1*t16*t21*t22*t28*t81
     +        -3.2405036462011455d+0*gammaaa*t1*t28*t31*t4*t80
     +        -3.78058758723467d+0*gammaaa*t1*t18*t21*t22*t80
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*
     +        (3.1502591643128747d+1*t10*t11*t12*t25*t48*t60*t78
     +        +2.4502015722433473d+1*t10*t11*t21*t29*t30*t6*t78
     +        -1.0500863881042916d+1*t10*t11*t12*t25*t30*t6*t78
     +        -1.0719631878564646d+2*t10*t11*t25*t30*t40*t6*t62
     +        +1.5751295821564373d+1*t10*t11*t12*t25*t48*t57*t60
     +        +1.2251007861216737d+1*t10*t11*t21*t29*t30*t57*t6
     +        -5.250431940521458d+0*t10*t11*t12*t25*t30*t57*t6
     +        +3.062751965304186d+0*t10*t11*t25*t29*t30*t42*t6
     +        +1.166762653449213d+0*t10*t11*t12*t24*t25*t30*t6
     +        -5.250431940521458d+0*t10*t11*t12*t167*t30*t6
     +        +5.833813267246065d-1*t10*t11*t12*t13*t55
     +        +1.5739142043924287d+0*t10*t11*t113*t13*t29+t168+t166
     +        +t165+t164+t163+t162+t161+t160+t159+t158+t149+t148+t147
     +        +t146+t144+t143+t142+t141+t135+t133+t132+t131+t128+t127
     +        +t125+t124+t122+t121)+9.721510938603437d+0*gammaaa*t1
     +        *t27*t4*t68*t69+3.2405036462011455d+0*gammaaa*t1*t16*t31
     +        *t4*t67*t68-1.6202518231005728d+0*gammaaa*t1*t27*t28*t4
     +        *t67+1.890293793617335d+0*gammaaa*t1*t16*t21*t22*t28*t67
     +        -1.6202518231005728d+0*gammaaa*t1*t28*t31*t4*t59
     +        -1.890293793617335d+0*gammaaa*t1*t18*t21*t22*t59
     +        +1.8002798034450812d-1*gammaaa*t1*t16*t18*t4*t55
     +        +1.8002798034450812d-1*gammaaa*t1*t16*t24*t28*t31*t4
     +        -1.8002798034450812d-1*gammaaa*t1*t18*t24*t27*t4
     +        +1.6202518231005728d+0*gammaaa*t1*(t168+t166+t165+t164
     +        +t163+t162+t161+t160+t159+t158+t135+t133+t132+t131+t128
     +        +t127+t125+t124+t122+t121)*t18*t4+2.756678449025281d-1
     +        *gammaaa*t1*t113*t16*t18*t22+t152+t150+t118+t117+t116
     +        +t115+t114+t112+t110)+3.109069086965489d-2*t71
     +        +1.0363563623218297d-2*t32*t5)-5.181781811609151d-3*t24
     +        *t33+t157+t156)*wght+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t205+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = 1.0d+0*(t20*(6.218138173930978d-2
     +        *rhoa*(-1.3126079851303646d+0*gammaaa*t10*t12*t28*t67*t8
     +        +1.4584533168115157d-1*gammaaa*t10*t12*t18*t24*t8
     +        +3.2405036462011455d+0*gammaaa*t1*t16*t36*t4*t67*t68
     +        -8.101259115502863d-1*t1*t16*t28*t4*t67
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -5.250431940521458d+0*gammaaa*t10*t12*t30*t57*t6
     +        +5.833813267246063d-1*gammaaa*t10*t12*t13*t24+t236+t235
     +        +t234+t233+t232+t231+t230+t229+t228+t227+t226+t225)
     +        -1.6202518231005728d+0*gammaaa*t1*t28*t36*t4*t59
     +        +8.101259115502863d-1*t1*t18*t4*t59-1.80027980344508d-1
     +        *gammaaa*t1*t16*t24*t28*t36*t4+9.0013990172254d-2*t1*t16
     +        *t18*t24*t4+1.6202518231005728d+0*gammaaa*t1*t18*(t231
     +        +t230+t229+t228+t227+t226+t225)*t4+t239+t238+t237+t224
     +        +t223+t222+t221+t220+t219+t218+t217+t216+t215+t214+t213
     +        +t212+t211+t210+t209+t208+t207)-1.0363563623218294d-2
     +        *t37*t5+t206)+(t240-t37*t72)*t73+t241)*wght+Cmat3(iq,
     +        D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t205+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = 1.0d+0*(t73*(t203-t35*t83)+t20*
     +        (6.218138173930978d-2*rhoa*(-6.563039925651822d-1
     +        *gammaaa*t10*t12*t28*t8*t81+3.2405036462011455d+0
     +        *gammaaa*t1*t16*t34*t4*t68*t81-4.050629557751432d-1*t1
     +        *t16*t28*t4*t81-1.6202518231005728d+0*gammaaa*t1*t28*t34
     +        *t4*t80+4.050629557751432d-1*t1*t18*t4*t80
     +        -7.292266584057581d-2*gammaaa*t10*t12*t18*t24*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -2.625215970260729d+0*gammaaa*t10*t12*t30*t6*t78+t243
     +        +t242-2.916906633623032d-1*gammaaa*t10*t12*t13*t24+t199
     +        +t198+t197+t196+t195+t193+t192+t190+t189+t188)
     +        +1.8002798034450807d-1*gammaaa*t1*t16*t24*t28*t34*t4
     +        +1.6202518231005728d+0*gammaaa*t1*t18*(t243+t242+t193
     +        +t192+t190+t189+t188)*t4-4.5006995086127016d-2*t1*t16
     +        *t18*t24*t4+t202+t201+t200+t187+t186+t185+t184+t183+t182
     +        +t181+t180+t179+t178+t177+t176+t175+t174+t173+t172+t171
     +        +t170)+1.0363563623218296d-2*t35*t5+t169)+t204)*wght
     +        +Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = 1.0d+0*(t73*(t240-t37*t83)+t20*
     +        (6.218138173930978d-2*rhoa*(-1.3126079851303646d+0
     +        *gammaaa*t10*t12*t28*t8*t81+3.2405036462011455d+0
     +        *gammaaa*t1*t16*t36*t4*t68*t81-8.101259115502863d-1*t1
     +        *t16*t28*t4*t81-1.6202518231005728d+0*gammaaa*t1*t28*t36
     +        *t4*t80+8.101259115502863d-1*t1*t18*t4*t80
     +        -1.458453316811516d-1*gammaaa*t10*t12*t18*t24*t8
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -5.250431940521458d+0*gammaaa*t10*t12*t30*t6*t78+t245
     +        +t244-5.833813267246064d-1*gammaaa*t10*t12*t13*t24+t236
     +        +t235+t234+t233+t232+t230+t229+t227+t226+t225)
     +        +1.8002798034450807d-1*gammaaa*t1*t16*t24*t28*t36*t4
     +        +1.6202518231005728d+0*gammaaa*t1*t18*(t245+t244+t230
     +        +t229+t227+t226+t225)*t4-9.001399017225403d-2*t1*t16*t18
     +        *t24*t4+t239+t238+t237+t224+t223+t222+t221+t220+t219
     +        +t218+t217+t216+t215+t214+t213+t212+t211+t210+t209+t208
     +        +t207)+1.0363563623218296d-2*t37*t5+t206)+t241)*wght
     +        +Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = t254+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = 1.0d+0*(t73*(
     +        -6.218138173930978d-2*rhoa*t35*t95-t37*t89+t257)+t256
     +        +t255)*wght+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t254+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = 1.0d+0*(t73*(-t37*t96
     +        -6.218138173930978d-2*rhoa*t37*t95-6.218138173930978d-2
     +        *rhoa*t107*t32)+t20*(6.218138173930978d-2*rhoa*(
     +        -2.625215970260729d+0*gammaaa*t10*t12*t28*t8*t94
     +        +6.481007292402291d+0*gammaaa*t1*t16*t36*t4*t68*t94
     +        -1.6202518231005728d+0*t1*t16*t28*t4*t94
     +        -3.2405036462011455d+0*gammaaa*t1*t28*t36*t4*t93
     +        +1.6202518231005728d+0*t1*t18*t4*t93
     +        +5.250431940521458d+0*gammaaa*t10*t12*t31*t36*t68*t8
     +        +3.062751965304184d+0*gammaaa*t10*t21*t28*t29*t36*t8
     +        -1.3126079851303646d+0*t10*t12*t28*t31*t8
     +        -1.531375982652092d+0*t10*t18*t21*t29*t8
     +        +4.2535109618476863d+0*gammaaa*t100*t13*t16*t31*t62*t68
     +        +3.2405036462011455d+0*t1*t16*t31*t36*t4*t68
     +        +3.2405036462011455d+0*gammaaa*t1*t106*t27*t4*t68
     +        -3.78058758723467d+0*gammaaa*t1*t106*t16*t21*t22*t68
     +        -2.1267554809238431d+0*gammaaa*t100*t13*t27*t28*t62
     +        -1.6202518231005728d+0*gammaaa*t1*t16*t28*t4*(
     +        -2.625215970260729d+0*t10*t12*t25*t30*t6
     +        -3.062751965304184d+0*t10*t13*t21*t29)
     +        -1.6202518231005728d+0*t1*t27*t28*t36*t4
     +        -9.721510938603435d+0*gammaaa*t1*t106*t151*t16*t31*t4
     +        +1.890293793617335d+0*t1*t16*t21*t22*t28*t36
     +        +2.481214727744483d+0*gammaaa*t100*t13*t145*t16*t21*t28)
     +        +3.109069086965489d-2*t107)+1.2436276347861956d-1*rhoa
     +        *t105*t155*t32)*wght+Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = 1.0d+0*(t73*(-t35*t96
     +        -6.218138173930978d-2*rhoa*t37*t88+t257)+t256+t255)*wght
     +        +Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t254+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t260+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = 6.218138173930978d-2*rhoa*t20*(
     +        -1.292204791970997d+0*gammaaa*t258*t259*t28*t30+t270
     +        +t269+t268+t267+t266+t265+t264+t263+t262+t261)*wght+t272
     +        +t271+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t260+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = t73*(-1.2436276347861956d-1
     +        *rhoa*t104*t37*wght-6.218138173930978d-2*rhoa*t107*t35
     +        *wght)+6.218138173930978d-2*rhoa*t20*
     +        (5.250431940521458d+0*gammaaa*t10*t12*t34*t36*t68*t8
     +        +1.3126079851303646d+0*gammaaa*t10*t106*t12*t68*t8
     +        -1.3126079851303646d+0*t10*t12*t28*t36*t8
     +        +4.2535109618476863d+0*gammaaa*t100*t13*t16*t36*t62*t68
     +        +3.2405036462011455d+0*t1*t16*t34*t36*t4*t68
     +        +8.101259115502863d-1*t1*t106*t16*t4*t68
     +        -1.5950666106928824d+0*t100*t13*t16*t28*t62
     +        -9.721510938603435d+0*gammaaa*t1*t106*t151*t16*t34*t4
     +        -2.5844095839419934d+0*gammaaa*t258*t259*t28*t30+t269
     +        +t264)*wght+1.2436276347861956d-1*rhoa*t105*t155*t35
     +        *wght+Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = 6.218138173930978d-2*rhoa*t20*(
     +        -1.2922047919709967d+0*gammaaa*t258*t259*t28*t30+t270
     +        +t269+t268+t267+t266+t265+t264+t263+t262+t261)*wght+t272
     +        +t271+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = 6.218138173930978d-2*rhoa*t20*
     +        (7.875647910782188d+0*gammaaa*t10*t106*t12*t68*t8
     +        -3.937823955391094d+0*t10*t12*t28*t36*t8
     +        +1.2760532885543058d+1*gammaaa*t100*t13*t16*t36*t62*t68
     +        +4.8607554693017174d+0*t1*t106*t16*t4*t68
     +        -3.1901332213857647d+0*t100*t13*t16*t28*t62
     +        -9.721510938603435d+0*gammaaa*t1*t151*t16*t36**3*t4
     +        -5.168819167883988d+0*gammaaa*t258*t259*t28*t30)*wght
     +        -1.8654414521792934d-1*rhoa*t107*t37*t73*wght
     +        +1.2436276347861956d-1*rhoa*t155*t37**3*wght+Cmat3(iq,
     +        D3_GAB_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = rhob
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
            call nwxcm_c_pw91lda_d3(t1,tol_rho,ipol,1,1.0d0,sr,s1f,
     +        s1a,s1a2,s1a3)
            t2 = s1f
            t3 = rhob+rhoa
            t4 = rhoa-rhob
            t5 = 1/t3
            t6 = 1.0d+0-t4*t5
            t7 = t4*t5+1.0d+0
            t8 = 5.0d-1*t7**6.666666666666666d-1+5.0d-1*t6*
     +        *6.666666666666666d-1
            t9 = t8**3
            t10 = gammabb+2.0d+0*gammaab+gammaaa
            t11 = t3**1.1666666666666667d+0
            t12 = 1/t11**2.0d+0
            t13 = 1/t8**2.0d+0
            t14 = 1/t9
            t15 = exp(-3.216396844291483d+1*t14*s1f*t5)
            t16 = t15-1.0d+0
            t17 = 1/t16
            t18 = 2.0413893780549558d+0*t1*t10*t12*t13*t17
            t19 = t18+1.0d+0
            t20 = t1**2
            t21 = t10**2
            t22 = 1/t11**4.0d+0
            t23 = 1/t8**4.0d+0
            t24 = 1/t16**2
            t25 = 4.1672705928355996d+0*t20*t21*t22*t23*t24+t18+1.0d+0
            t26 = 1/t25
            t27 = 2.0413893780549558d+0*t1*t10*t12*t13*t19*t26+1.0d+0
            t28 = log(t27)
            t29 = s1a(D1_RA)
            t30 = 1/t27
            t31 = 1/t3**2
            t32 = t31*t4
            t33 = -t5
            t34 = t33+t32
            t35 = 1/t6**3.333333333333333d-1
            t36 = -t31*t4
            t37 = t5+t36
            t38 = 1/t7**3.333333333333333d-1
            t39 = 3.333333333333333d-1*t37*t38+3.333333333333333d-1
     +        *t34*t35
            t40 = 1/t8**3.0d+0
            t41 = t3**1.6666666666666674d-1
            t42 = 1/t11**3.0d+0
            t43 = -4.763241882128231d+0*t1*t10*t13*t19*t26*t41*t42
            t44 = -4.0827787561099116d+0*t1*t10*t12*t17*t39*t40
            t45 = -4.763241882128231d+0*t1*t10*t13*t17*t41*t42
            t46 = 1/t8**4
            t47 = 3.216396844291483d+1*t14*s1f*t31
            t48 = 9.649190532874448d+1*s1f*t39*t46*t5
     +        -3.216396844291483d+1*t14*t29*t5+t47
            t49 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t48
            t50 = t49+t45+t44
            t51 = 1/t25**2
            t52 = 1/t8**5.0d+0
            t53 = 1/t11**5.0d+0
            t54 = -1.9447262766566134d+1*t20*t21*t23*t24*t41*t53
            t55 = 1/t16**3
            t56 = -8.334541185671199d+0*t15*t20*t21*t22*t23*t48*t55
     +        +t54-1.6669082371342397d+1*t20*t21*t22*t24*t39*t52+t49
     +        +t45+t44
            t57 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t56
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*t50+t43
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t39*t40
            t58 = t8**2
            t59 = 3.109069086965489d-2*t28*t9
            t60 = s1a(D1_RB)
            t61 = t5+t32
            t62 = t36+t33
            t63 = 3.333333333333333d-1*t38*t62+3.333333333333333d-1
     +        *t35*t61
            t64 = -4.0827787561099116d+0*t1*t10*t12*t17*t40*t63
            t65 = 9.649190532874448d+1*s1f*t46*t5*t63
     +        -3.216396844291483d+1*t14*t5*t60+t47
            t66 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t65
            t67 = t66+t64+t45
            t68 = t66-8.334541185671199d+0*t15*t20*t21*t22*t23*t55*t65
     +        +t64-1.6669082371342397d+1*t20*t21*t22*t24*t52*t63+t54+t45
            t69 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t68
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t40*t63+t43
            t70 = 8.334541185671199d+0*t10*t20*t22*t23*t24
     +        +2.0413893780549558d+0*t1*t12*t13*t17
            t71 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t70
     +        +4.1672705928355996d+0*t10*t17*t20*t22*t23*t26
     +        +2.0413893780549558d+0*t1*t12*t13*t19*t26
            t72 = 3.109069086965489d-2*t3*t30*t71*t9*wght
            t73 = 1.6669082371342397d+1*t10*t20*t22*t23*t24
     +        +4.0827787561099116d+0*t1*t12*t13*t17
            t74 = -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*t73
     +        +8.334541185671199d+0*t10*t17*t20*t22*t23*t26
     +        +4.0827787561099116d+0*t1*t12*t13*t19*t26
            t75 = s1a2(D2_RA_RA)
            t76 = t39**2
            t77 = t34**2
            t78 = 1/t6**1.3333333333333333d+0
            t79 = 1/t3**3
            t80 = -2*t4*t79
            t81 = 2*t31
            t82 = t81+t80
            t83 = t37**2
            t84 = 1/t7**1.3333333333333333d+0
            t85 = 2*t4*t79
            t86 = -2*t31
            t87 = t86+t85
            t88 = 3.333333333333333d-1*t38*t87-1.111111111111111d-1
     +        *t83*t84+3.333333333333333d-1*t35*t82
     +        -1.111111111111111d-1*t77*t78
            t89 = t3**3.333333333333335d-1
            t90 = 1.667134658744881d+1*t1*t10*t13*t19*t22*t26*t89
            t91 = 1/t3**8.333333333333333d-1
            t92 = -7.938736470213722d-1*t1*t10*t13*t19*t26*t42*t91
            t93 = t48**2
            t94 = exp(-6.432793688582966d+1*t14*s1f*t5)
            t95 = 4.0827787561099116d+0*t1*t10*t12*t13*t55*t93*t94
            t96 = 1.2248336268329735d+1*t1*t10*t12*t17*t23*t76
            t97 = -4.0827787561099116d+0*t1*t10*t12*t17*t40*t88
            t98 = 1.9052967528512924d+1*t1*t10*t17*t39*t40*t41*t42
            t99 = 1.667134658744881d+1*t1*t10*t13*t17*t22*t89
            t100 = -7.938736470213722d-1*t1*t10*t13*t17*t42*t91
            t101 = 8.165557512219823d+0*t1*t10*t12*t15*t24*t39*t40*t48
            t102 = 9.526483764256461d+0*t1*t10*t13*t15*t24*t41*t42*t48
            t103 = 1/t8**5
            t104 = -6.432793688582966d+1*t14*s1f*t79
            t105 = 9.649190532874448d+1*s1f*t46*t5*t88
     +        -3.859676213149779d+2*t103*s1f*t5*t76
     +        -3.216396844291483d+1*t14*t5*t75+1.9298381065748896d+2
     +        *t29*t39*t46*t5-1.9298381065748896d+2*s1f*t31*t39*t46
     +        +6.432793688582966d+1*t14*t29*t31+t104
            t106 = -2.0413893780549558d+0*t1*t10*t105*t12*t13*t15*t24
            t107 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t93
            t108 = t99+t98+t97+t96+t95+t107+t106+t102+t101+t100
            t109 = 1/t16**4
            t110 = 1/t8**6.0d+0
            t111 = 1/t11**6.0d+0
            t112 = 1.1344236613830246d+2*t111*t20*t21*t23*t24*t89
            t113 = -3.241210461094357d+0*t20*t21*t23*t24*t53*t91
            t114 = t99+t98+t97+t96+t95+2.50036235570136d+1*t109*t20
     +        *t21*t22*t23*t93*t94-8.334541185671199d+0*t15*t20*t21
     +        *t22*t23*t55*t93-1.6669082371342397d+1*t20*t21*t22*t24
     +        *t52*t88+8.334541185671199d+1*t110*t20*t21*t22*t24*t76
     +        +7.778905106626454d+1*t15*t20*t21*t23*t41*t48*t53*t55
     +        +6.667632948536959d+1*t15*t20*t21*t22*t39*t48*t52*t55
     +        -8.334541185671199d+0*t105*t15*t20*t21*t22*t23*t55
     +        +1.5557810213252907d+2*t20*t21*t24*t39*t41*t52*t53+t113
     +        +t112+t107+t106+t102+t101+t100
            t115 = 1/t25**3
            t116 = t56**2
            t117 = t92+t90-4.0827787561099116d+0*t1*t10*t12*t19*t26
     +        *t40*t88+1.2248336268329735d+1*t1*t10*t12*t19*t23*t26
     +        *t76-4.0827787561099116d+0*t1*t10*t12*t13*t50*t51*t56
     +        +9.526483764256461d+0*t1*t10*t13*t19*t41*t42*t51*t56
     +        +8.165557512219823d+0*t1*t10*t12*t19*t39*t40*t51*t56
     +        -2.0413893780549558d+0*t1*t10*t114*t12*t13*t19*t51
     +        -9.526483764256461d+0*t1*t10*t13*t26*t41*t42*t50
     +        -8.165557512219823d+0*t1*t10*t12*t26*t39*t40*t50
     +        +1.9052967528512924d+1*t1*t10*t19*t26*t39*t40*t41*t42
     +        +2.0413893780549558d+0*t1*t10*t108*t12*t13*t26
     +        +4.0827787561099116d+0*t1*t10*t115*t116*t12*t13*t19
            t118 = 6.218138173930978d-2*t57*t9+3.109069086965489d-2
     +        *t117*t3*t9+1.8654414521792934d-1*t3*t39*t57*t58
            t119 = 1/t27**2
            t120 = t57**2
            t121 = s1a2(D2_RA_RB)
            t122 = -1.111111111111111d-1*t37*t62*t84
     +        +6.666666666666666d-1*t38*t4*t79-6.666666666666666d-1
     +        *t35*t4*t79-1.111111111111111d-1*t34*t61*t78
            t123 = 4.0827787561099116d+0*t1*t10*t12*t13*t48*t55*t65*t94
            t124 = 1.2248336268329735d+1*t1*t10*t12*t17*t23*t39*t63
            t125 = -4.0827787561099116d+0*t1*t10*t12*t122*t17*t40
            t126 = 9.526483764256461d+0*t1*t10*t17*t40*t41*t42*t63
            t127 = 9.526483764256461d+0*t1*t10*t17*t39*t40*t41*t42
            t128 = 4.0827787561099116d+0*t1*t10*t12*t15*t24*t40*t48*t63
            t129 = 4.763241882128231d+0*t1*t10*t13*t15*t24*t41*t42*t48
            t130 = 9.649190532874448d+1*t29*t46*t5*t63
     +        -3.859676213149779d+2*t103*s1f*t39*t5*t63
     +        -9.649190532874448d+1*s1f*t31*t46*t63
     +        +9.649190532874448d+1*t39*t46*t5*t60
     +        +3.216396844291483d+1*t14*t31*t60+9.649190532874448d+1
     +        *t122*s1f*t46*t5-3.216396844291483d+1*t121*t14*t5
     +        -9.649190532874448d+1*s1f*t31*t39*t46
     +        +3.216396844291483d+1*t14*t29*t31+t104
            t131 = -2.0413893780549558d+0*t1*t10*t12*t13*t130*t15*t24
            t132 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t48*t65
            t133 = 4.0827787561099116d+0*t1*t10*t12*t15*t24*t39*t40*t65
            t134 = 4.763241882128231d+0*t1*t10*t13*t15*t24*t41*t42*t65
            t135 = t99+t134+t133+t132+t131+t129+t128+t127+t126+t125
     +        +t124+t123+t100
            t136 = t99+2.50036235570136d+1*t109*t20*t21*t22*t23*t48
     +        *t65*t94+3.889452553313227d+1*t15*t20*t21*t23*t41*t53
     +        *t55*t65+3.3338164742684795d+1*t15*t20*t21*t22*t39*t52
     +        *t55*t65-8.334541185671199d+0*t15*t20*t21*t22*t23*t48
     +        *t55*t65+3.3338164742684795d+1*t15*t20*t21*t22*t48*t52
     +        *t55*t63+7.778905106626454d+1*t20*t21*t24*t41*t52*t53
     +        *t63+8.334541185671199d+1*t110*t20*t21*t22*t24*t39*t63
     +        +3.889452553313227d+1*t15*t20*t21*t23*t41*t48*t53*t55
     +        -8.334541185671199d+0*t130*t15*t20*t21*t22*t23*t55
     +        +7.778905106626454d+1*t20*t21*t24*t39*t41*t52*t53
     +        -1.6669082371342397d+1*t122*t20*t21*t22*t24*t52+t134
     +        +t133+t132+t131+t129+t128+t127+t126+t125+t124+t123+t113
     +        +t112+t100
            t137 = t92+t90+4.0827787561099116d+0*t1*t10*t115*t12*t13
     +        *t19*t56*t68-2.0413893780549558d+0*t1*t10*t12*t13*t50
     +        *t51*t68+4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51
     +        *t68+4.0827787561099116d+0*t1*t10*t12*t19*t39*t40*t51
     +        *t68-2.0413893780549558d+0*t1*t10*t12*t13*t51*t56*t67
     +        -4.763241882128231d+0*t1*t10*t13*t26*t41*t42*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t26*t39*t40*t67
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t56*t63
     +        -4.0827787561099116d+0*t1*t10*t12*t26*t40*t50*t63
     +        +9.526483764256461d+0*t1*t10*t19*t26*t40*t41*t42*t63
     +        +1.2248336268329735d+1*t1*t10*t12*t19*t23*t26*t39*t63
     +        +4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t56
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t136*t19*t51
     +        -4.763241882128231d+0*t1*t10*t13*t26*t41*t42*t50
     +        +9.526483764256461d+0*t1*t10*t19*t26*t39*t40*t41*t42
     +        -4.0827787561099116d+0*t1*t10*t12*t122*t19*t26*t40
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t135*t26
            t138 = 3.109069086965489d-2*t69*t9+3.109069086965489d-2
     +        *t57*t9+3.109069086965489d-2*t137*t3*t9
     +        +9.327207260896466d-2*t3*t39*t58*t69
     +        +9.327207260896466d-2*t3*t57*t58*t63
            t139 = s1a2(D2_RB_RB)
            t140 = t63**2
            t141 = t61**2
            t142 = t86+t80
            t143 = t62**2
            t144 = t85+t81
            t145 = -1.111111111111111d-1*t143*t84-1.111111111111111d-1
     +        *t141*t78+3.333333333333333d-1*t144*t38
     +        +3.333333333333333d-1*t142*t35
            t146 = t65**2
            t147 = 4.0827787561099116d+0*t1*t10*t12*t13*t146*t55*t94
            t148 = 1.2248336268329735d+1*t1*t10*t12*t140*t17*t23
            t149 = -4.0827787561099116d+0*t1*t10*t12*t145*t17*t40
            t150 = 1.9052967528512924d+1*t1*t10*t17*t40*t41*t42*t63
            t151 = 8.165557512219823d+0*t1*t10*t12*t15*t24*t40*t63*t65
            t152 = 9.526483764256461d+0*t1*t10*t13*t15*t24*t41*t42*t65
            t153 = 1.9298381065748896d+2*t46*t5*t60*t63
     +        -1.9298381065748896d+2*s1f*t31*t46*t63
     +        +6.432793688582966d+1*t14*t31*t60+9.649190532874448d+1
     +        *t145*s1f*t46*t5-3.859676213149779d+2*t103*t140*s1f*t5
     +        -3.216396844291483d+1*t139*t14*t5+t104
            t154 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t153*t24
            t155 = -2.0413893780549558d+0*t1*t10*t12*t13*t146*t15*t24
            t156 = t99+t155+t154+t152+t151+t150+t149+t148+t147+t100
            t157 = t99+2.50036235570136d+1*t109*t146*t20*t21*t22*t23
     +        *t94+6.667632948536959d+1*t15*t20*t21*t22*t52*t55*t63
     +        *t65+7.778905106626454d+1*t15*t20*t21*t23*t41*t53*t55
     +        *t65+1.5557810213252907d+2*t20*t21*t24*t41*t52*t53*t63
     +        -8.334541185671199d+0*t15*t153*t20*t21*t22*t23*t55
     +        -8.334541185671199d+0*t146*t15*t20*t21*t22*t23*t55
     +        -1.6669082371342397d+1*t145*t20*t21*t22*t24*t52
     +        +8.334541185671199d+1*t110*t140*t20*t21*t22*t24+t155
     +        +t154+t152+t151+t150+t149+t148+t147+t113+t112+t100
            t158 = t68**2
            t159 = t92+t90-4.0827787561099116d+0*t1*t10*t12*t13*t51
     +        *t67*t68+8.165557512219823d+0*t1*t10*t12*t19*t40*t51*t63
     +        *t68+9.526483764256461d+0*t1*t10*t13*t19*t41*t42*t51*t68
     +        -8.165557512219823d+0*t1*t10*t12*t26*t40*t63*t67
     +        -9.526483764256461d+0*t1*t10*t13*t26*t41*t42*t67
     +        +1.9052967528512924d+1*t1*t10*t19*t26*t40*t41*t42*t63
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t157*t19*t51
     +        -4.0827787561099116d+0*t1*t10*t12*t145*t19*t26*t40
     +        +1.2248336268329735d+1*t1*t10*t12*t140*t19*t23*t26
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t156*t26
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t158*t19
            t160 = 6.218138173930978d-2*t69*t9+3.109069086965489d-2
     +        *t159*t3*t9+1.8654414521792934d-1*t3*t58*t63*t69
            t161 = t69**2
            t162 = 3.109069086965489d-2*t71*t9
            t163 = 4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t70
            t164 = -4.763241882128231d+0*t1*t13*t19*t26*t41*t42
            t165 = -9.723631383283067d+0*t10*t17*t20*t23*t26*t41*t53
            t166 = -4.0827787561099116d+0*t1*t12*t17*t39*t40
            t167 = -4.763241882128231d+0*t1*t13*t17*t41*t42
            t168 = -2.0413893780549558d+0*t1*t12*t13*t15*t24*t48
            t169 = t168+t167+t166
            t170 = -3.889452553313227d+1*t10*t20*t23*t24*t41*t53
            t171 = -1.6669082371342397d+1*t10*t15*t20*t22*t23*t48*t55
     +        -3.3338164742684795d+1*t10*t20*t22*t24*t39*t52+t170+t168
     +        +t167+t166
            t172 = 4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t56
     +        *t70-2.0413893780549558d+0*t1*t10*t12*t13*t50*t51*t70
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t39*t40*t51*t70
     +        -4.1672705928355996d+0*t10*t17*t20*t22*t23*t51*t56
     +        -2.0413893780549558d+0*t1*t12*t13*t19*t51*t56
     +        -8.334541185671199d+0*t10*t17*t20*t22*t26*t39*t52
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t171*t19*t51
     +        +2.0413893780549558d+0*t1*t12*t13*t26*t50
     +        -4.0827787561099116d+0*t1*t12*t19*t26*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t169*t26+t165+t164
     +        +t163
            t173 = 3.109069086965489d-2*t172*t3*t9
     +        +9.327207260896466d-2*t3*t39*t58*t71+t162
            t174 = 1.0d+0*(t173*t30-3.109069086965489d-2*t119*t3*t57
     +        *t71*t9)*wght
            t175 = 3.109069086965489d-2*t74*t9
            t176 = 4.763241882128231d+0*t1*t10*t13*t19*t41*t42*t51*t73
            t177 = -9.526483764256461d+0*t1*t13*t19*t26*t41*t42
            t178 = -1.9447262766566134d+1*t10*t17*t20*t23*t26*t41*t53
            t179 = -8.165557512219823d+0*t1*t12*t17*t39*t40
            t180 = -9.526483764256461d+0*t1*t13*t17*t41*t42
            t181 = -4.0827787561099116d+0*t1*t12*t13*t15*t24*t48
            t182 = t181+t180+t179
            t183 = -7.778905106626454d+1*t10*t20*t23*t24*t41*t53
            t184 = -3.3338164742684795d+1*t10*t15*t20*t22*t23*t48*t55
     +        -6.667632948536959d+1*t10*t20*t22*t24*t39*t52+t183+t181
     +        +t180+t179
            t185 = 4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t56
     +        *t73-2.0413893780549558d+0*t1*t10*t12*t13*t50*t51*t73
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t39*t40*t51*t73
     +        -8.334541185671199d+0*t10*t17*t20*t22*t23*t51*t56
     +        -4.0827787561099116d+0*t1*t12*t13*t19*t51*t56
     +        -1.6669082371342397d+1*t10*t17*t20*t22*t26*t39*t52
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t184*t19*t51
     +        +4.0827787561099116d+0*t1*t12*t13*t26*t50
     +        -8.165557512219823d+0*t1*t12*t19*t26*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t182*t26+t178+t177
     +        +t176
            t186 = 3.109069086965489d-2*t185*t3*t9
     +        +9.327207260896466d-2*t3*t39*t58*t74+t175
            t187 = -4.0827787561099116d+0*t1*t12*t17*t40*t63
            t188 = -2.0413893780549558d+0*t1*t12*t13*t15*t24*t65
            t189 = t188+t187+t167
            t190 = -1.6669082371342397d+1*t10*t15*t20*t22*t23*t55*t65
     +        -3.3338164742684795d+1*t10*t20*t22*t24*t52*t63+t188+t187
     +        +t170+t167
            t191 = 4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t68
     +        *t70-2.0413893780549558d+0*t1*t10*t12*t13*t51*t67*t70
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t63*t70
     +        -4.1672705928355996d+0*t10*t17*t20*t22*t23*t51*t68
     +        -2.0413893780549558d+0*t1*t12*t13*t19*t51*t68
     +        +2.0413893780549558d+0*t1*t12*t13*t26*t67
     +        -8.334541185671199d+0*t10*t17*t20*t22*t26*t52*t63
     +        -4.0827787561099116d+0*t1*t12*t19*t26*t40*t63
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t190*t51
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t189*t26+t165+t164
     +        +t163
            t192 = 3.109069086965489d-2*t191*t3*t9
     +        +9.327207260896466d-2*t3*t58*t63*t71+t162
            t193 = 1.0d+0*(t192*t30-3.109069086965489d-2*t119*t3*t69
     +        *t71*t9)*wght
            t194 = -8.165557512219823d+0*t1*t12*t17*t40*t63
            t195 = -4.0827787561099116d+0*t1*t12*t13*t15*t24*t65
            t196 = t195+t194+t180
            t197 = -3.3338164742684795d+1*t10*t15*t20*t22*t23*t55*t65
     +        -6.667632948536959d+1*t10*t20*t22*t24*t52*t63+t195+t194
     +        +t183+t180
            t198 = 4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t68
     +        *t73-2.0413893780549558d+0*t1*t10*t12*t13*t51*t67*t73
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t63*t73
     +        -8.334541185671199d+0*t10*t17*t20*t22*t23*t51*t68
     +        -4.0827787561099116d+0*t1*t12*t13*t19*t51*t68
     +        +4.0827787561099116d+0*t1*t12*t13*t26*t67
     +        -1.6669082371342397d+1*t10*t17*t20*t22*t26*t52*t63
     +        -8.165557512219823d+0*t1*t12*t19*t26*t40*t63
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t197*t51
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t196*t26+t178+t177
     +        +t176
            t199 = 3.109069086965489d-2*t198*t3*t9
     +        +9.327207260896466d-2*t3*t58*t63*t74+t175
            t200 = t71**2
            t201 = t70**2
            t202 = -4.0827787561099116d+0*t1*t12*t13*t19*t51*t70
            t203 = t1**3
            t204 = -8.334541185671199d+0*t10*t17*t20*t22*t23*t51*t70
            t205 = -1.7014043847390745d+1*t10*t110*t111*t19*t203*t24
     +        *t51+8.334541185671199d+0*t17*t20*t22*t23*t26+t204+t202
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t201
            t206 = 3.109069086965489d-2*t205*t3*t30*t9*wght
     +        -3.109069086965489d-2*t119*t200*t3*t9*wght
            t207 = 4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t70
     +        *t73-4.1672705928355996d+0*t10*t17*t20*t22*t23*t51*t73
     +        -2.0413893780549558d+0*t1*t12*t13*t19*t51*t73
     +        -3.402808769478149d+1*t10*t110*t111*t19*t203*t24*t51
     +        +1.6669082371342397d+1*t17*t20*t22*t23*t26+t204+t202
            t208 = 3.109069086965489d-2*t207*t3*t30*t9*wght
     +        -3.109069086965489d-2*t119*t3*t71*t74*t9*wght
            t209 = t74**2
            t210 = t73**2
            t211 = -1.6669082371342397d+1*t10*t17*t20*t22*t23*t51*t73
     +        -8.165557512219823d+0*t1*t12*t13*t19*t51*t73
     +        -6.805617538956298d+1*t10*t110*t111*t19*t203*t24*t51
     +        +3.3338164742684795d+1*t17*t20*t22*t23*t26
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t210
            t212 = s1a3(D3_RA_RA_RA)
            t213 = t39**3
            t214 = 1/t6**2.3333333333333334d+0
            t215 = 1/t3**4
            t216 = 6*t215*t4
            t217 = -6*t79
            t218 = 1/t7**2.3333333333333334d+0
            t219 = -6*t215*t4
            t220 = 6*t79
            t221 = -3.333333333333333d-1*t37*t84*t87
     +        -3.333333333333333d-1*t34*t78*t82+3.333333333333333d-1*
     +        (t220+t219)*t38+1.4814814814814814d-1*t218*t37**3
     +        +3.333333333333333d-1*(t217+t216)*t35
     +        +1.4814814814814814d-1*t214*t34**3
            t222 = t3**5.000000000000002d-1
            t223 = -7.779961740809445d+1*t1*t10*t13*t19*t222*t26*t53
            t224 = 1/t3**6.666666666666665d-1
            t225 = 8.335673293724408d+0*t1*t10*t13*t19*t22*t224*t26
            t226 = 1/t3**1.8333333333333332d+0
            t227 = 6.615613725178101d-1*t1*t10*t13*t19*t226*t26*t42
            t228 = t48**3
            t229 = exp(-9.649190532874448d+1*t14*s1f*t5)
            t230 = -1.2248336268329735d+1*t1*t10*t109*t12*t13*t228*t229
            t231 = 1.2248336268329735d+1*t1*t10*t105*t12*t13*t48*t55*t94
            t232 = 6.432793688582966d+1*t14*s1f*t31
            t233 = 1.9298381065748896d+2*s1f*t39*t46*t5
     +        -6.432793688582966d+1*t14*t29*t5+t232
            t234 = 4.0827787561099116d+0*t1*t10*t12*t13*t233*t55*t93*t94
            t235 = -2.449667253665947d+1*t1*t10*t12*t39*t40*t55*t93*t94
            t236 = -2.8579451292769387d+1*t1*t10*t13*t41*t42*t55*t93*t94
            t237 = 4.0827787561099116d+0*t1*t10*t12*t13*t228*t55*t94
            t238 = -4.899334507331894d+1*t1*t10*t12*t17*t213*t52
            t239 = 3.6745008804989204d+1*t1*t10*t12*t17*t23*t39*t88
            t240 = -8.573835387830815d+1*t1*t10*t17*t23*t41*t42*t76
            t241 = -4.0827787561099116d+0*t1*t10*t12*t17*t221*t40
            t242 = 2.8579451292769387d+1*t1*t10*t17*t40*t41*t42*t88
            t243 = -1.0002807952469285d+2*t1*t10*t17*t22*t39*t40*t89
            t244 = 4.7632418821282335d+0*t1*t10*t17*t39*t40*t42*t91
            t245 = -7.779961740809445d+1*t1*t10*t13*t17*t222*t53
            t246 = 8.335673293724408d+0*t1*t10*t13*t17*t22*t224
            t247 = 6.615613725178101d-1*t1*t10*t13*t17*t226*t42
            t248 = -3.6745008804989204d+1*t1*t10*t12*t15*t23*t24*t48*t76
            t249 = 1.2248336268329735d+1*t1*t10*t12*t15*t24*t40*t48*t88
            t250 = -5.715890258553877d+1*t1*t10*t15*t24*t39*t40*t41
     +        *t42*t48
            t251 = -5.001403976234643d+1*t1*t10*t13*t15*t22*t24*t48*t89
            t252 = 2.3816209410641168d+0*t1*t10*t13*t15*t24*t42*t48*t91
            t253 = -6.124168134164868d+0*t1*t10*t105*t12*t13*t15*t24*t48
            t254 = 1.2248336268329735d+1*t1*t10*t105*t12*t15*t24*t39*t40
            t255 = 1.4289725646384693d+1*t1*t10*t105*t13*t15*t24*t41*t42
            t256 = 1/t8**6
            t257 = 1.9298381065748896d+2*t14*s1f*t215
            t258 = 2.8947571598623345d+2*t29*t46*t5*t88
     +        -1.1579028639449339d+3*t103*s1f*t39*t5*t88
     +        -2.8947571598623345d+2*s1f*t31*t46*t88
     +        +5.789514319724669d+2*s1f*t39*t46*t79
     +        -1.9298381065748896d+2*t14*t29*t79-1.1579028639449339d+3
     +        *t103*t29*t5*t76+1.1579028639449339d+3*t103*s1f*t31*t76
     +        +2.8947571598623345d+2*t39*t46*t5*t75
     +        +9.649190532874448d+1*t14*t31*t75+9.649190532874448d+1
     +        *s1f*t221*t46*t5+1.9298381065748896d+3*s1f*t213*t256*t5
     +        -3.216396844291483d+1*t14*t212*t5-5.789514319724669d+2
     +        *t29*t31*t39*t46+t257
            t259 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t258
            t260 = 1.2248336268329735d+1*t1*t10*t12*t15*t24*t39*t40*t93
            t261 = 1.4289725646384693d+1*t1*t10*t13*t15*t24*t41*t42*t93
            t262 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t228*t24
            t263 = 1/t16**5
            t264 = 1/t8**7.0d+0
            t265 = 1/t11**7.0d+0
            t266 = -7.940965629681173d+2*t20*t21*t222*t23*t24*t265
            t267 = 5.672118306915125d+1*t111*t20*t21*t224*t23*t24
            t268 = 2.7010087175786307d+0*t20*t21*t226*t23*t24*t53
            t269 = 1/t25**4
            t270 = -3.109069086965489d-2*t120*t9
            t271 = 1/t27**3
            t272 = s1a3(D3_RA_RA_RB)
            t273 = 1.8654414521792934d-1*t57*t58*t63
            t274 = 1.8654414521792934d-1*t39*t58*t69
            t275 = 6.218138173930978d-2*t137*t9
            t276 = -5.715890258553877d+1*t1*t10*t19*t23*t26*t39*t41
     +        *t42*t63
            t277 = 1.9052967528512924d+1*t1*t10*t122*t19*t26*t40*t41*t42
            t278 = -1.111111111111111d-1*t62*t84*t87
     +        -4.444444444444444d-1*t37*t4*t79*t84
     +        +1.4814814814814814d-1*t218*t62*t83-1.111111111111111d-1
     +        *t61*t78*t82+3.333333333333333d-1*t38*(2*t79+t219)
     +        +4.444444444444444d-1*t34*t4*t78*t79
     +        +3.333333333333333d-1*t35*(t216-2*t79)
     +        +1.4814814814814814d-1*t214*t61*t77
            t279 = 1.9052967528512924d+1*t1*t10*t26*t40*t41*t42*t50*t63
            t280 = -1.9052967528512924d+1*t1*t10*t19*t40*t41*t42*t51
     +        *t56*t63
            t281 = -9.526483764256461d+0*t1*t10*t13*t135*t26*t41*t42
            t282 = 9.526483764256461d+0*t1*t10*t13*t136*t19*t41*t42*t51
            t283 = 1.9052967528512924d+1*t1*t10*t26*t39*t40*t41*t42*t67
            t284 = 9.526483764256461d+0*t1*t10*t13*t41*t42*t51*t56*t67
            t285 = -1.9052967528512924d+1*t1*t10*t19*t39*t40*t41*t42
     +        *t51*t68
            t286 = 9.526483764256461d+0*t1*t10*t13*t41*t42*t50*t51*t68
            t287 = -1.9052967528512924d+1*t1*t10*t115*t13*t19*t41*t42
     +        *t56*t68
            t288 = -1.2248336268329735d+1*t1*t10*t109*t12*t13*t229*t65
     +        *t93
            t289 = 8.165557512219823d+0*t1*t10*t12*t13*t130*t48*t55*t94
            t290 = 1.9298381065748896d+2*s1f*t46*t5*t63
     +        -6.432793688582966d+1*t14*t5*t60+t232
            t291 = 4.0827787561099116d+0*t1*t10*t12*t13*t290*t55*t93*t94
            t292 = -1.6331115024439646d+1*t1*t10*t12*t39*t40*t48*t55
     +        *t65*t94
            t293 = -1.9052967528512924d+1*t1*t10*t13*t41*t42*t48*t55
     +        *t65*t94
            t294 = 4.0827787561099116d+0*t1*t10*t105*t12*t13*t55*t65*t94
            t295 = 4.0827787561099116d+0*t1*t10*t12*t13*t55*t65*t93*t94
            t296 = -8.165557512219823d+0*t1*t10*t12*t40*t55*t63*t93*t94
            t297 = -9.526483764256461d+0*t1*t10*t13*t41*t42*t55*t93*t94
            t298 = -4.899334507331894d+1*t1*t10*t12*t17*t52*t63*t76
            t299 = 1.2248336268329735d+1*t1*t10*t12*t17*t23*t63*t88
            t300 = 2.449667253665947d+1*t1*t10*t12*t122*t17*t23*t39
            t301 = -5.715890258553877d+1*t1*t10*t17*t23*t39*t41*t42*t63
            t302 = -2.8579451292769387d+1*t1*t10*t17*t23*t41*t42*t76
            t303 = 1.9052967528512924d+1*t1*t10*t122*t17*t40*t41*t42
            t304 = -4.0827787561099116d+0*t1*t10*t12*t17*t278*t40
            t305 = 9.526483764256461d+0*t1*t10*t17*t40*t41*t42*t88
            t306 = -3.334269317489762d+1*t1*t10*t17*t22*t40*t63*t89
            t307 = 1.5877472940427445d+0*t1*t10*t17*t40*t42*t63*t91
            t308 = -6.668538634979524d+1*t1*t10*t17*t22*t39*t40*t89
            t309 = 3.175494588085489d+0*t1*t10*t17*t39*t40*t42*t91
            t310 = -2.449667253665947d+1*t1*t10*t12*t15*t23*t24*t39
     +        *t48*t63
            t311 = 8.165557512219823d+0*t1*t10*t12*t122*t15*t24*t40*t48
            t312 = -1.9052967528512924d+1*t1*t10*t15*t24*t40*t41*t42
     +        *t48*t63
            t313 = -1.9052967528512924d+1*t1*t10*t15*t24*t39*t40*t41
     +        *t42*t48
            t314 = -3.334269317489762d+1*t1*t10*t13*t15*t22*t24*t48*t89
            t315 = 1.5877472940427445d+0*t1*t10*t13*t15*t24*t42*t48*t91
            t316 = -4.0827787561099116d+0*t1*t10*t12*t13*t130*t15*t24
     +        *t48
            t317 = 8.165557512219823d+0*t1*t10*t12*t130*t15*t24*t39*t40
            t318 = 9.526483764256461d+0*t1*t10*t13*t130*t15*t24*t41*t42
            t319 = 4.0827787561099116d+0*t1*t10*t105*t12*t15*t24*t40*t63
            t320 = 4.763241882128231d+0*t1*t10*t105*t13*t15*t24*t41*t42
            t321 = 7.719352426299558d+2*t103*s1f*t31*t39*t63
            t322 = -1.9298381065748896d+2*t122*s1f*t31*t46
            t323 = -1.9298381065748896d+2*t29*t31*t46*t63
            t324 = -1.9298381065748896d+2*t31*t39*t46*t60
            t325 = 6.432793688582966d+1*t121*t14*t31
            t326 = -3.859676213149779d+2*t103*s1f*t5*t63*t88
     +        +9.649190532874448d+1*t46*t5*t60*t88
     +        -9.649190532874448d+1*s1f*t31*t46*t88
     +        +1.9298381065748896d+2*s1f*t46*t63*t79
     +        -6.432793688582966d+1*t14*t60*t79+3.859676213149779d+2
     +        *s1f*t39*t46*t79-1.2865587377165932d+2*t14*t29*t79
     +        +1.9298381065748896d+3*s1f*t256*t5*t63*t76
     +        -3.859676213149779d+2*t103*t5*t60*t76
     +        +3.859676213149779d+2*t103*s1f*t31*t76
     +        +9.649190532874448d+1*t46*t5*t63*t75
     +        +3.216396844291483d+1*t14*t31*t75-7.719352426299558d+2
     +        *t103*t29*t39*t5*t63+1.9298381065748896d+2*t121*t39*t46
     +        *t5+1.9298381065748896d+2*t122*t29*t46*t5
     +        +9.649190532874448d+1*s1f*t278*t46*t5
     +        -7.719352426299558d+2*t103*t122*s1f*t39*t5
     +        -3.216396844291483d+1*t14*t272*t5-1.9298381065748896d+2
     +        *t29*t31*t39*t46+t325+t324+t323+t322+t321+t257
            t327 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t326
            t328 = -1.2248336268329735d+1*t1*t10*t12*t15*t23*t24*t65*t76
            t329 = 8.165557512219823d+0*t1*t10*t12*t15*t24*t39*t40*t48
     +        *t65
            t330 = 9.526483764256461d+0*t1*t10*t13*t15*t24*t41*t42*t48
     +        *t65
            t331 = -2.0413893780549558d+0*t1*t10*t105*t12*t13*t15*t24
     +        *t65
            t332 = 4.0827787561099116d+0*t1*t10*t12*t15*t24*t40*t65*t88
            t333 = -1.9052967528512924d+1*t1*t10*t15*t24*t39*t40*t41
     +        *t42*t65
            t334 = -1.667134658744881d+1*t1*t10*t13*t15*t22*t24*t65*t89
            t335 = 7.938736470213722d-1*t1*t10*t13*t15*t24*t42*t65*t91
            t336 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t65*t93
            t337 = 4.0827787561099116d+0*t1*t10*t12*t15*t24*t40*t63*t93
            t338 = 4.763241882128231d+0*t1*t10*t13*t15*t24*t41*t42*t93
            t339 = -2.3336715319879361d+2*t109*t20*t21*t23*t41*t48*t53
     +        *t65*t94
            t340 = -7.778905106626453d+2*t110*t20*t21*t24*t39*t41*t53
     +        *t63
            t341 = 1.5557810213252907d+2*t122*t20*t21*t24*t41*t52*t53
            t342 = -3.1115620426505813d+2*t15*t20*t21*t41*t48*t52*t53
     +        *t55*t63
            t343 = 7.778905106626454d+1*t130*t15*t20*t21*t23*t41*t53*t55
            t344 = -3.1115620426505813d+2*t15*t20*t21*t39*t41*t52*t53
     +        *t55*t65
            t345 = 7.778905106626454d+1*t15*t20*t21*t23*t41*t48*t53
     +        *t55*t65
            t346 = 3.7308829043585867d-1*t28*t39*t63*t8
            t347 = 1.8654414521792934d-1*t122*t28*t58
            t348 = s1a3(D3_RA_RB_RB)
            t349 = -4.444444444444444d-1*t4*t62*t79*t84
     +        -1.111111111111111d-1*t144*t37*t84+4.444444444444444d-1
     +        *t4*t61*t78*t79-6.666666666666666d-1*t38*t79
     +        +6.666666666666666d-1*t35*t79-1.111111111111111d-1*t142
     +        *t34*t78-2.0d+0*t215*t38*t4+2.0d+0*t215*t35*t4
     +        +1.4814814814814814d-1*t143*t218*t37
     +        +1.4814814814814814d-1*t141*t214*t34
            t350 = -1.2248336268329735d+1*t1*t10*t109*t12*t13*t146
     +        *t229*t48
            t351 = -1.6331115024439646d+1*t1*t10*t12*t40*t48*t55*t63
     +        *t65*t94
            t352 = 8.165557512219823d+0*t1*t10*t12*t13*t130*t55*t65*t94
            t353 = 4.0827787561099116d+0*t1*t10*t12*t13*t290*t48*t55
     +        *t65*t94
            t354 = 4.0827787561099116d+0*t1*t10*t12*t13*t153*t48*t55*t94
            t355 = 4.0827787561099116d+0*t1*t10*t12*t13*t146*t48*t55*t94
            t356 = -8.165557512219823d+0*t1*t10*t12*t146*t39*t40*t55*t94
            t357 = -9.526483764256461d+0*t1*t10*t13*t146*t41*t42*t55*t94
            t358 = -4.899334507331894d+1*t1*t10*t12*t140*t17*t39*t52
            t359 = 2.449667253665947d+1*t1*t10*t12*t122*t17*t23*t63
            t360 = 1.2248336268329735d+1*t1*t10*t12*t145*t17*t23*t39
            t361 = -2.8579451292769387d+1*t1*t10*t140*t17*t23*t41*t42
            t362 = -4.0827787561099116d+0*t1*t10*t12*t17*t349*t40
            t363 = 9.526483764256461d+0*t1*t10*t145*t17*t40*t41*t42
            t364 = -6.668538634979524d+1*t1*t10*t17*t22*t40*t63*t89
            t365 = 3.175494588085489d+0*t1*t10*t17*t40*t42*t63*t91
            t366 = -3.334269317489762d+1*t1*t10*t17*t22*t39*t40*t89
            t367 = 1.5877472940427445d+0*t1*t10*t17*t39*t40*t42*t91
            t368 = -1.2248336268329735d+1*t1*t10*t12*t140*t15*t23*t24
     +        *t48
            t369 = 4.0827787561099116d+0*t1*t10*t12*t145*t15*t24*t40*t48
            t370 = -1.667134658744881d+1*t1*t10*t13*t15*t22*t24*t48*t89
            t371 = 7.938736470213722d-1*t1*t10*t13*t15*t24*t42*t48*t91
            t372 = 8.165557512219823d+0*t1*t10*t12*t130*t15*t24*t40*t63
            t373 = 3.859676213149779d+2*s1f*t46*t63*t79
     +        -1.2865587377165932d+2*t14*t60*t79+1.9298381065748896d+2
     +        *s1f*t39*t46*t79-6.432793688582966d+1*t14*t29*t79
     +        -7.719352426299558d+2*t103*t39*t5*t60*t63
     +        -1.9298381065748896d+2*t31*t46*t60*t63
     +        +1.9298381065748896d+2*t121*t46*t5*t63
     +        -7.719352426299558d+2*t103*t122*s1f*t5*t63
     +        +1.9298381065748896d+2*t122*t46*t5*t60
     +        +9.649190532874448d+1*t139*t39*t46*t5
     +        +9.649190532874448d+1*s1f*t349*t46*t5
     +        +9.649190532874448d+1*t145*t29*t46*t5
     +        +1.9298381065748896d+3*t140*s1f*t256*t39*t5
     +        -3.859676213149779d+2*t103*t145*s1f*t39*t5
     +        -3.216396844291483d+1*t14*t348*t5-3.859676213149779d+2
     +        *t103*t140*t29*t5-9.649190532874448d+1*t145*s1f*t31*t46
     +        +t325+t324+t323+t322+t321+3.859676213149779d+2*t103*t140
     +        *s1f*t31+3.216396844291483d+1*t139*t14*t31+t257
            t374 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t373
            t375 = -2.449667253665947d+1*t1*t10*t12*t15*t23*t24*t39
     +        *t63*t65
            t376 = 8.165557512219823d+0*t1*t10*t12*t15*t24*t40*t48*t63
     +        *t65
            t377 = -4.0827787561099116d+0*t1*t10*t12*t13*t130*t15*t24
     +        *t65
            t378 = 8.165557512219823d+0*t1*t10*t12*t122*t15*t24*t40*t65
            t379 = -1.9052967528512924d+1*t1*t10*t15*t24*t40*t41*t42
     +        *t63*t65
            t380 = -3.334269317489762d+1*t1*t10*t13*t15*t22*t24*t65*t89
            t381 = 1.5877472940427445d+0*t1*t10*t13*t15*t24*t42*t65*t91
            t382 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t153*t24
     +        *t48
            t383 = 4.0827787561099116d+0*t1*t10*t12*t15*t153*t24*t39*t40
            t384 = 4.763241882128231d+0*t1*t10*t13*t15*t153*t24*t41*t42
            t385 = -2.0413893780549558d+0*t1*t10*t12*t13*t146*t15*t24
     +        *t48
            t386 = 4.0827787561099116d+0*t1*t10*t12*t146*t15*t24*t39*t40
            t387 = 4.763241882128231d+0*t1*t10*t13*t146*t15*t24*t41*t42
            t388 = s1a3(D3_RB_RB_RB)
            t389 = t63**3
            t390 = -3.333333333333333d-1*t144*t62*t84
     +        -3.333333333333333d-1*t142*t61*t78+1.4814814814814814d-1
     +        *t218*t62**3+1.4814814814814814d-1*t214*t61**3
     +        +3.333333333333333d-1*(t219+t217)*t38
     +        +3.333333333333333d-1*(t220+t216)*t35
            t391 = t65**3
            t392 = -1.2248336268329735d+1*t1*t10*t109*t12*t13*t229*t391
            t393 = 1.2248336268329735d+1*t1*t10*t12*t13*t153*t55*t65*t94
            t394 = 4.0827787561099116d+0*t1*t10*t12*t13*t146*t290*t55
     +        *t94
            t395 = -2.449667253665947d+1*t1*t10*t12*t146*t40*t55*t63*t94
            t396 = -2.8579451292769387d+1*t1*t10*t13*t146*t41*t42*t55
     +        *t94
            t397 = 4.0827787561099116d+0*t1*t10*t12*t13*t391*t55*t94
            t398 = -4.899334507331894d+1*t1*t10*t12*t17*t389*t52
            t399 = 3.6745008804989204d+1*t1*t10*t12*t145*t17*t23*t63
            t400 = -8.573835387830815d+1*t1*t10*t140*t17*t23*t41*t42
            t401 = -4.0827787561099116d+0*t1*t10*t12*t17*t390*t40
            t402 = 2.8579451292769387d+1*t1*t10*t145*t17*t40*t41*t42
            t403 = -1.0002807952469285d+2*t1*t10*t17*t22*t40*t63*t89
            t404 = 4.7632418821282335d+0*t1*t10*t17*t40*t42*t63*t91
            t405 = -3.6745008804989204d+1*t1*t10*t12*t140*t15*t23*t24
     +        *t65
            t406 = 1.2248336268329735d+1*t1*t10*t12*t145*t15*t24*t40*t65
            t407 = -5.715890258553877d+1*t1*t10*t15*t24*t40*t41*t42
     +        *t63*t65
            t408 = -5.001403976234643d+1*t1*t10*t13*t15*t22*t24*t65*t89
            t409 = 2.3816209410641168d+0*t1*t10*t13*t15*t24*t42*t65*t91
            t410 = -6.124168134164868d+0*t1*t10*t12*t13*t15*t153*t24*t65
            t411 = 1.2248336268329735d+1*t1*t10*t12*t15*t153*t24*t40*t63
            t412 = 1.4289725646384693d+1*t1*t10*t13*t15*t153*t24*t41*t42
            t413 = 5.789514319724669d+2*s1f*t46*t63*t79
     +        -1.9298381065748896d+2*t14*t60*t79-5.789514319724669d+2
     +        *t31*t46*t60*t63+2.8947571598623345d+2*t139*t46*t5*t63
     +        -1.1579028639449339d+3*t103*t145*s1f*t5*t63
     +        +2.8947571598623345d+2*t145*t46*t5*t60
     +        -1.1579028639449339d+3*t103*t140*t5*t60
     +        +9.649190532874448d+1*s1f*t390*t46*t5
     +        +1.9298381065748896d+3*s1f*t256*t389*t5
     +        -3.216396844291483d+1*t14*t388*t5-2.8947571598623345d+2
     +        *t145*s1f*t31*t46+1.1579028639449339d+3*t103*t140*s1f
     +        *t31+9.649190532874448d+1*t139*t14*t31+t257
            t414 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t413
            t415 = 1.2248336268329735d+1*t1*t10*t12*t146*t15*t24*t40*t63
            t416 = 1.4289725646384693d+1*t1*t10*t13*t146*t15*t24*t41*t42
            t417 = -2.0413893780549558d+0*t1*t10*t12*t13*t15*t24*t391
            t418 = -1.667134658744881d+1*t1*t10*t13*t19*t22*t51*t70*t89
            t419 = 7.938736470213722d-1*t1*t10*t13*t19*t42*t51*t70*t91
            t420 = 1.9052967528512924d+1*t1*t19*t26*t39*t40*t41*t42
            t421 = 1.667134658744881d+1*t1*t13*t19*t22*t26*t89
            t422 = -7.938736470213722d-1*t1*t13*t19*t26*t42*t91
            t423 = 3.889452553313227d+1*t10*t17*t20*t26*t39*t41*t52*t53
            t424 = 3.403270984149074d+1*t10*t111*t17*t20*t23*t26*t89
            t425 = -1.6206052305471785d+0*t10*t17*t20*t23*t26*t53*t91
            t426 = -9.526483764256461d+0*t1*t13*t26*t41*t42*t50
            t427 = 9.526483764256461d+0*t1*t13*t19*t41*t42*t51*t56
            t428 = 1.9447262766566134d+1*t10*t17*t20*t23*t41*t51*t53*t56
            t429 = 4.0827787561099116d+0*t1*t12*t13*t55*t93*t94
            t430 = 1.2248336268329735d+1*t1*t12*t17*t23*t76
            t431 = -4.0827787561099116d+0*t1*t12*t17*t40*t88
            t432 = 1.9052967528512924d+1*t1*t17*t39*t40*t41*t42
            t433 = 1.667134658744881d+1*t1*t13*t17*t22*t89
            t434 = -7.938736470213722d-1*t1*t13*t17*t42*t91
            t435 = 8.165557512219823d+0*t1*t12*t15*t24*t39*t40*t48
            t436 = 9.526483764256461d+0*t1*t13*t15*t24*t41*t42*t48
            t437 = -2.0413893780549558d+0*t1*t105*t12*t13*t15*t24
            t438 = -2.0413893780549558d+0*t1*t12*t13*t15*t24*t93
            t439 = 3.1115620426505813d+2*t10*t20*t24*t39*t41*t52*t53
            t440 = 2.268847322766049d+2*t10*t111*t20*t23*t24*t89
            t441 = -6.482420922188714d+0*t10*t20*t23*t24*t53*t91
            t442 = 1.5557810213252907d+2*t10*t15*t20*t23*t41*t48*t53*t55
            t443 = 1.0d+0*(t30*(3.109069086965489d-2*t3*t9*(
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (5.00072471140272d+1*t10*t109*t20*t22*t23*t93*t94
     +        -1.6669082371342397d+1*t10*t15*t20*t22*t23*t55*t93
     +        -3.3338164742684795d+1*t10*t20*t22*t24*t52*t88
     +        +1.6669082371342397d+2*t10*t110*t20*t22*t24*t76
     +        +1.333526589707392d+2*t10*t15*t20*t22*t39*t48*t52*t55
     +        -1.6669082371342397d+1*t10*t105*t15*t20*t22*t23*t55+t442
     +        +t441+t440+t439+t438+t437+t436+t435+t434+t433+t432+t431
     +        +t430+t429)+4.0827787561099116d+0*t1*t10*t12*t19*t40*t51
     +        *t70*t88-8.334541185671199d+0*t10*t17*t20*t22*t26*t52
     +        *t88-4.0827787561099116d+0*t1*t12*t19*t26*t40*t88
     +        -1.2248336268329735d+1*t1*t10*t12*t19*t23*t51*t70*t76
     +        +1.2248336268329735d+1*t1*t12*t19*t23*t26*t76
     +        +2.50036235570136d+1*t10*t110*t17*t20*t22*t26*t76
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t50*t56*t70
     +        -1.9052967528512924d+1*t1*t10*t115*t13*t19*t41*t42*t56
     +        *t70-1.6331115024439646d+1*t1*t10*t115*t12*t19*t39*t40
     +        *t56*t70+9.526483764256461d+0*t1*t10*t13*t41*t42*t50*t51
     +        *t70+8.165557512219823d+0*t1*t10*t12*t39*t40*t50*t51*t70
     +        -1.9052967528512924d+1*t1*t10*t19*t39*t40*t41*t42*t51
     +        *t70-2.0413893780549558d+0*t1*t10*t108*t12*t13*t51*t70
     +        -1.2248336268329735d+1*t1*t10*t116*t12*t13*t19*t269*t70
     +        +4.0827787561099116d+0*t1*t10*t114*t115*t12*t13*t19*t70
     +        +1.6669082371342397d+1*t10*t17*t20*t22*t39*t51*t52*t56
     +        -4.0827787561099116d+0*t1*t12*t13*t50*t51*t56
     +        +8.165557512219823d+0*t1*t12*t19*t39*t40*t51*t56
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t169*t51*t56
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t171*t19*t56
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t171*t50*t51
     +        +9.526483764256461d+0*t1*t10*t13*t171*t19*t41*t42*t51
     +        +8.165557512219823d+0*t1*t10*t12*t171*t19*t39*t40*t51
     +        -4.1672705928355996d+0*t10*t114*t17*t20*t22*t23*t51
     +        -2.0413893780549558d+0*t1*t114*t12*t13*t19*t51
     +        -8.165557512219823d+0*t1*t12*t26*t39*t40*t50
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t438+t437
     +        +t436+t435+t434+t433+t432+t431+t430+t429)+t428+t427+t426
     +        +t425+t424+t423+t422+t421+t420-9.526483764256461d+0*t1
     +        *t10*t13*t169*t26*t41*t42+t419+t418-8.165557512219823d+0
     +        *t1*t10*t12*t169*t26*t39*t40+2.0413893780549558d+0*t1
     +        *t108*t12*t13*t26+8.334541185671199d+0*t10*t115*t116*t17
     +        *t20*t22*t23+4.0827787561099116d+0*t1*t115*t116*t12*t13
     +        *t19)+6.218138173930978d-2*t172*t9+9.327207260896466d-2
     +        *t3*t58*t71*t88+1.8654414521792934d-1*t3*t71*t76*t8
     +        +1.8654414521792934d-1*t39*t58*t71+1.8654414521792934d-1
     +        *t172*t3*t39*t58)+t119*(-6.218138173930978d-2*t172*t3
     +        *t57*t9-t118*t71)+6.218138173930978d-2*t120*t271*t3*t71
     +        *t9)*wght
            t444 = -1.667134658744881d+1*t1*t10*t13*t19*t22*t51*t73*t89
            t445 = 7.938736470213722d-1*t1*t10*t13*t19*t42*t51*t73*t91
            t446 = 3.334269317489762d+1*t1*t13*t19*t22*t26*t89
            t447 = -1.5877472940427445d+0*t1*t13*t19*t26*t42*t91
            t448 = 6.806541968298148d+1*t10*t111*t17*t20*t23*t26*t89
            t449 = -3.241210461094357d+0*t10*t17*t20*t23*t26*t53*t91
            t450 = 8.165557512219823d+0*t1*t12*t13*t55*t93*t94
            t451 = 2.449667253665947d+1*t1*t12*t17*t23*t76
            t452 = -8.165557512219823d+0*t1*t12*t17*t40*t88
            t453 = 3.8105935057025847d+1*t1*t17*t39*t40*t41*t42
            t454 = 3.334269317489762d+1*t1*t13*t17*t22*t89
            t455 = -1.5877472940427445d+0*t1*t13*t17*t42*t91
            t456 = 1.6331115024439646d+1*t1*t12*t15*t24*t39*t40*t48
            t457 = 1.9052967528512924d+1*t1*t13*t15*t24*t41*t42*t48
            t458 = -4.0827787561099116d+0*t1*t105*t12*t13*t15*t24
            t459 = -4.0827787561099116d+0*t1*t12*t13*t15*t24*t93
            t460 = 4.537694645532098d+2*t10*t111*t20*t23*t24*t89
            t461 = -1.2964841844377428d+1*t10*t20*t23*t24*t53*t91
            t462 = 4.0827787561099116d+0*t1*t12*t13*t48*t55*t65*t94
            t463 = 1.2248336268329735d+1*t1*t12*t17*t23*t39*t63
            t464 = -4.0827787561099116d+0*t1*t12*t122*t17*t40
            t465 = 9.526483764256461d+0*t1*t17*t40*t41*t42*t63
            t466 = 9.526483764256461d+0*t1*t17*t39*t40*t41*t42
            t467 = 4.0827787561099116d+0*t1*t12*t15*t24*t40*t48*t63
            t468 = 4.763241882128231d+0*t1*t13*t15*t24*t41*t42*t48
            t469 = -2.0413893780549558d+0*t1*t12*t13*t130*t15*t24
            t470 = -2.0413893780549558d+0*t1*t12*t13*t15*t24*t48*t65
            t471 = 4.0827787561099116d+0*t1*t12*t15*t24*t39*t40*t65
            t472 = 4.763241882128231d+0*t1*t13*t15*t24*t41*t42*t65
            t473 = 1.0d+0*(t30*(3.109069086965489d-2*t3*t9*(
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (5.00072471140272d+1*t10*t109*t20*t22*t23*t48*t65*t94
     +        +7.778905106626454d+1*t10*t15*t20*t23*t41*t53*t55*t65
     +        +6.667632948536959d+1*t10*t15*t20*t22*t39*t52*t55*t65
     +        -1.6669082371342397d+1*t10*t15*t20*t22*t23*t48*t55*t65
     +        +6.667632948536959d+1*t10*t15*t20*t22*t48*t52*t55*t63
     +        +1.5557810213252907d+2*t10*t20*t24*t41*t52*t53*t63
     +        +1.6669082371342397d+2*t10*t110*t20*t22*t24*t39*t63
     +        +7.778905106626454d+1*t10*t15*t20*t23*t41*t48*t53*t55
     +        -1.6669082371342397d+1*t10*t130*t15*t20*t22*t23*t55
     +        +1.5557810213252907d+2*t10*t20*t24*t39*t41*t52*t53
     +        -3.3338164742684795d+1*t10*t122*t20*t22*t24*t52+t472
     +        +t471+t470+t469+t468+t467+t466+t465+t464+t463+t462+t441
     +        +t440+t434+t433)-1.2248336268329735d+1*t1*t10*t12*t13
     +        *t19*t269*t56*t68*t70+4.0827787561099116d+0*t1*t10*t115
     +        *t12*t13*t50*t68*t70-9.526483764256461d+0*t1*t10*t115
     +        *t13*t19*t41*t42*t68*t70-8.165557512219823d+0*t1*t10
     +        *t115*t12*t19*t39*t40*t68*t70+4.0827787561099116d+0*t1
     +        *t10*t115*t12*t13*t56*t67*t70+4.763241882128231d+0*t1
     +        *t10*t13*t41*t42*t51*t67*t70+4.0827787561099116d+0*t1
     +        *t10*t12*t39*t40*t51*t67*t70-8.165557512219823d+0*t1*t10
     +        *t115*t12*t19*t40*t56*t63*t70+4.0827787561099116d+0*t1
     +        *t10*t12*t40*t50*t51*t63*t70-9.526483764256461d+0*t1*t10
     +        *t19*t40*t41*t42*t51*t63*t70-1.2248336268329735d+1*t1
     +        *t10*t12*t19*t23*t39*t51*t63*t70-9.526483764256461d+0*t1
     +        *t10*t115*t13*t19*t41*t42*t56*t70+4.763241882128231d+0
     +        *t1*t10*t13*t41*t42*t50*t51*t70-9.526483764256461d+0*t1
     +        *t10*t19*t39*t40*t41*t42*t51*t70+4.0827787561099116d+0
     +        *t1*t10*t12*t122*t19*t40*t51*t70-2.0413893780549558d+0
     +        *t1*t10*t12*t13*t135*t51*t70+4.0827787561099116d+0*t1
     +        *t10*t115*t12*t13*t136*t19*t70+8.334541185671199d+0*t10
     +        *t115*t17*t20*t22*t23*t56*t68+4.0827787561099116d+0*t1
     +        *t115*t12*t13*t19*t56*t68+9.723631383283067d+0*t10*t17
     +        *t20*t23*t41*t51*t53*t68+8.334541185671199d+0*t10*t17
     +        *t20*t22*t39*t51*t52*t68-2.0413893780549558d+0*t1*t12
     +        *t13*t50*t51*t68+4.763241882128231d+0*t1*t13*t19*t41*t42
     +        *t51*t68+4.0827787561099116d+0*t1*t12*t19*t39*t40*t51
     +        *t68-2.0413893780549558d+0*t1*t10*t12*t13*t169*t51*t68
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t171*t19*t68
     +        -2.0413893780549558d+0*t1*t12*t13*t51*t56*t67
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t171*t51*t67
     +        -4.763241882128231d+0*t1*t13*t26*t41*t42*t67
     +        -4.0827787561099116d+0*t1*t12*t26*t39*t40*t67
     +        +8.334541185671199d+0*t10*t17*t20*t22*t51*t52*t56*t63
     +        +4.0827787561099116d+0*t1*t12*t19*t40*t51*t56*t63
     +        +1.9447262766566134d+1*t10*t17*t20*t26*t41*t52*t53*t63
     +        +4.0827787561099116d+0*t1*t10*t12*t171*t19*t40*t51*t63
     +        -4.0827787561099116d+0*t1*t12*t26*t40*t50*t63
     +        +9.526483764256461d+0*t1*t19*t26*t40*t41*t42*t63
     +        -4.0827787561099116d+0*t1*t10*t12*t169*t26*t40*t63
     +        +1.2248336268329735d+1*t1*t12*t19*t23*t26*t39*t63
     +        +2.50036235570136d+1*t10*t110*t17*t20*t22*t26*t39*t63
     +        +9.723631383283067d+0*t10*t17*t20*t23*t41*t51*t53*t56
     +        +4.763241882128231d+0*t1*t13*t19*t41*t42*t51*t56
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t189*t51*t56
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t190*t56
     +        +1.9447262766566134d+1*t10*t17*t20*t26*t39*t41*t52*t53
     +        -8.334541185671199d+0*t10*t122*t17*t20*t22*t26*t52
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t190*t50*t51
     +        +4.763241882128231d+0*t1*t10*t13*t19*t190*t41*t42*t51
     +        +4.763241882128231d+0*t1*t10*t13*t171*t19*t41*t42*t51
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t190*t39*t40*t51
     +        -4.1672705928355996d+0*t10*t136*t17*t20*t22*t23*t51
     +        -2.0413893780549558d+0*t1*t12*t13*t136*t19*t51
     +        -4.763241882128231d+0*t1*t13*t26*t41*t42*t50
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t472+t471
     +        +t470+t469+t468+t467+t466+t465+t464+t463+t462+t434+t433)
     +        +t425+t424+t422+t421+9.526483764256461d+0*t1*t19*t26*t39
     +        *t40*t41*t42-4.763241882128231d+0*t1*t10*t13*t189*t26
     +        *t41*t42-4.763241882128231d+0*t1*t10*t13*t169*t26*t41
     +        *t42+t419+t418-4.0827787561099116d+0*t1*t10*t12*t189*t26
     +        *t39*t40-4.0827787561099116d+0*t1*t12*t122*t19*t26*t40
     +        +2.0413893780549558d+0*t1*t12*t13*t135*t26)
     +        +3.109069086965489d-2*t191*t9+3.109069086965489d-2*t172
     +        *t9+1.8654414521792934d-1*t3*t39*t63*t71*t8
     +        +9.327207260896466d-2*t58*t63*t71+9.327207260896466d-2
     +        *t39*t58*t71+9.327207260896466d-2*t122*t3*t58*t71
     +        +9.327207260896466d-2*t172*t3*t58*t63
     +        +9.327207260896466d-2*t191*t3*t39*t58)+t119*(
     +        -3.109069086965489d-2*t172*t3*t69*t9
     +        -3.109069086965489d-2*t191*t3*t57*t9-t138*t71)
     +        +6.218138173930978d-2*t271*t3*t57*t69*t71*t9)*wght
            t474 = 1.9052967528512924d+1*t1*t19*t26*t40*t41*t42*t63
            t475 = 3.889452553313227d+1*t10*t17*t20*t26*t41*t52*t53*t63
            t476 = 8.165557512219823d+0*t1*t12*t13*t48*t55*t65*t94
            t477 = 2.449667253665947d+1*t1*t12*t17*t23*t39*t63
            t478 = -8.165557512219823d+0*t1*t12*t122*t17*t40
            t479 = 1.9052967528512924d+1*t1*t17*t40*t41*t42*t63
            t480 = 8.165557512219823d+0*t1*t12*t15*t24*t40*t48*t63
            t481 = -4.0827787561099116d+0*t1*t12*t13*t130*t15*t24
            t482 = -4.0827787561099116d+0*t1*t12*t13*t15*t24*t48*t65
            t483 = 8.165557512219823d+0*t1*t12*t15*t24*t39*t40*t65
            t484 = 9.526483764256461d+0*t1*t13*t15*t24*t41*t42*t65
            t485 = 3.1115620426505813d+2*t10*t20*t24*t41*t52*t53*t63
            t486 = 1.5557810213252907d+2*t10*t15*t20*t23*t41*t53*t55*t65
            t487 = -9.526483764256461d+0*t1*t13*t26*t41*t42*t67
            t488 = 9.526483764256461d+0*t1*t13*t19*t41*t42*t51*t68
            t489 = 1.9447262766566134d+1*t10*t17*t20*t23*t41*t51*t53*t68
            t490 = 4.0827787561099116d+0*t1*t12*t13*t146*t55*t94
            t491 = 1.2248336268329735d+1*t1*t12*t140*t17*t23
            t492 = -4.0827787561099116d+0*t1*t12*t145*t17*t40
            t493 = 8.165557512219823d+0*t1*t12*t15*t24*t40*t63*t65
            t494 = -2.0413893780549558d+0*t1*t12*t13*t15*t153*t24
            t495 = -2.0413893780549558d+0*t1*t12*t13*t146*t15*t24
            t496 = 1.0d+0*(t30*(3.109069086965489d-2*t3*t9*(
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (5.00072471140272d+1*t10*t109*t146*t20*t22*t23*t94
     +        +1.333526589707392d+2*t10*t15*t20*t22*t52*t55*t63*t65
     +        -1.6669082371342397d+1*t10*t15*t153*t20*t22*t23*t55
     +        -1.6669082371342397d+1*t10*t146*t15*t20*t22*t23*t55
     +        -3.3338164742684795d+1*t10*t145*t20*t22*t24*t52+t495
     +        +t494+t493+t492+t491+t490+t486+t485+t484+t479+t441+t440
     +        +t434+t433+1.6669082371342397d+2*t10*t110*t140*t20*t22
     +        *t24)+8.165557512219823d+0*t1*t10*t115*t12*t13*t67*t68
     +        *t70-1.6331115024439646d+1*t1*t10*t115*t12*t19*t40*t63
     +        *t68*t70-1.9052967528512924d+1*t1*t10*t115*t13*t19*t41
     +        *t42*t68*t70+8.165557512219823d+0*t1*t10*t12*t40*t51*t63
     +        *t67*t70+9.526483764256461d+0*t1*t10*t13*t41*t42*t51*t67
     +        *t70-1.9052967528512924d+1*t1*t10*t19*t40*t41*t42*t51
     +        *t63*t70+4.0827787561099116d+0*t1*t10*t12*t145*t19*t40
     +        *t51*t70-1.2248336268329735d+1*t1*t10*t12*t140*t19*t23
     +        *t51*t70-2.0413893780549558d+0*t1*t10*t12*t13*t156*t51
     +        *t70-1.2248336268329735d+1*t1*t10*t12*t13*t158*t19*t269
     +        *t70+4.0827787561099116d+0*t1*t10*t115*t12*t13*t157*t19
     +        *t70-4.0827787561099116d+0*t1*t12*t13*t51*t67*t68
     +        +1.6669082371342397d+1*t10*t17*t20*t22*t51*t52*t63*t68
     +        +8.165557512219823d+0*t1*t12*t19*t40*t51*t63*t68
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t189*t51*t68
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t19*t190*t68
     +        -8.165557512219823d+0*t1*t12*t26*t40*t63*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t190*t51*t67
     +        +8.165557512219823d+0*t1*t10*t12*t19*t190*t40*t51*t63
     +        -8.165557512219823d+0*t1*t10*t12*t189*t26*t40*t63
     +        -8.334541185671199d+0*t10*t145*t17*t20*t22*t26*t52
     +        +9.526483764256461d+0*t1*t10*t13*t19*t190*t41*t42*t51
     +        -4.1672705928355996d+0*t10*t157*t17*t20*t22*t23*t51
     +        -2.0413893780549558d+0*t1*t12*t13*t157*t19*t51
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t495+t494
     +        +t493+t492+t491+t490+t484+t479+t434+t433)+t489+t488+t487
     +        +t475+t474+t425+t424+t422+t421-9.526483764256461d+0*t1
     +        *t10*t13*t189*t26*t41*t42+t419+t418
     +        -4.0827787561099116d+0*t1*t12*t145*t19*t26*t40
     +        +1.2248336268329735d+1*t1*t12*t140*t19*t23*t26
     +        +2.50036235570136d+1*t10*t110*t140*t17*t20*t22*t26
     +        +2.0413893780549558d+0*t1*t12*t13*t156*t26
     +        +8.334541185671199d+0*t10*t115*t158*t17*t20*t22*t23
     +        +4.0827787561099116d+0*t1*t115*t12*t13*t158*t19)
     +        +6.218138173930978d-2*t191*t9+1.8654414521792934d-1*t140
     +        *t3*t71*t8+1.8654414521792934d-1*t58*t63*t71
     +        +9.327207260896466d-2*t145*t3*t58*t71
     +        +1.8654414521792934d-1*t191*t3*t58*t63)+t119*(
     +        -6.218138173930978d-2*t191*t3*t69*t9-t160*t71)
     +        +6.218138173930978d-2*t161*t271*t3*t71*t9)*wght
            t497 = 8.165557512219823d+0*t1*t12*t13*t146*t55*t94
            t498 = 2.449667253665947d+1*t1*t12*t140*t17*t23
            t499 = -8.165557512219823d+0*t1*t12*t145*t17*t40
            t500 = 3.8105935057025847d+1*t1*t17*t40*t41*t42*t63
            t501 = 1.6331115024439646d+1*t1*t12*t15*t24*t40*t63*t65
            t502 = 1.9052967528512924d+1*t1*t13*t15*t24*t41*t42*t65
            t503 = -4.0827787561099116d+0*t1*t12*t13*t15*t153*t24
            t504 = -4.0827787561099116d+0*t1*t12*t13*t146*t15*t24
            t505 = 3.109069086965489d-2*t205*t9
            t506 = -9.526483764256461d+0*t1*t10*t115*t13*t19*t201*t41
     +        *t42
            t507 = 8.165557512219823d+0*t1*t12*t19*t39*t40*t51*t70
            t508 = 9.526483764256461d+0*t1*t13*t19*t41*t42*t51*t70
            t509 = 3.969943564391174d+1*t10*t110*t19*t203*t24*t265*t41
     +        *t51
            t510 = 1.6669082371342397d+1*t10*t17*t20*t22*t39*t51*t52*t70
            t511 = 1.9447262766566134d+1*t10*t17*t20*t23*t41*t51*t53*t70
            t512 = -1.9447262766566134d+1*t17*t20*t23*t26*t41*t53
            t513 = -3.889452553313227d+1*t20*t23*t24*t41*t53
            t514 = 4.0827787561099116d+0*t1*t12*t13*t169*t26
            t515 = -4.0827787561099116d+0*t1*t12*t13*t171*t19*t51
            t516 = -8.334541185671199d+0*t10*t17*t171*t20*t22*t23*t51
            t517 = -4.0827787561099116d+0*t1*t12*t13*t50*t51*t70
            t518 = 8.165557512219823d+0*t1*t115*t12*t13*t19*t56*t70
            t519 = 1.6669082371342397d+1*t10*t115*t17*t20*t22*t23*t56
     +        *t70
            t520 = 1.0d+0*(t119*(-3.109069086965489d-2*t172*t3*t71*t9
     +        -3.109069086965489d-2*t205*t3*t57*t9-t173*t71)+t30*
     +        (3.109069086965489d-2*t3*(-4.0827787561099116d+0*t1*t10
     +        *t12*t13*t169*t51*t70+8.165557512219823d+0*t1*t10*t115
     +        *t12*t13*t171*t19*t70-8.334541185671199d+0*t17*t20*t22
     +        *t23*t51*t56-1.2248336268329735d+1*t1*t10*t12*t13*t19
     +        *t201*t269*t56+3.402808769478149d+1*t10*t110*t111*t115
     +        *t19*t203*t24*t56-2.0413893780549558d+0*t1*t10*t12*t13
     +        *t19*t51*(-1.6669082371342397d+1*t15*t20*t22*t23*t48*t55
     +        -3.3338164742684795d+1*t20*t22*t24*t39*t52+t513)
     +        -1.6669082371342397d+1*t17*t20*t22*t26*t39*t52+t519+t518
     +        +t517+t516+t515+t514+t512+t511+t510
     +        -1.7014043847390745d+1*t10*t110*t111*t203*t24*t50*t51
     +        +3.402808769478149d+1*t10*t111*t19*t203*t24*t264*t39*t51
     +        +t509+t508+t507+t506+4.0827787561099116d+0*t1*t10*t115
     +        *t12*t13*t201*t50-8.165557512219823d+0*t1*t10*t115*t12
     +        *t19*t201*t39*t40)*t9+9.327207260896466d-2*t205*t3*t39
     +        *t58+t505)+6.218138173930978d-2*t200*t271*t3*t57*t9)*wght
            t521 = 6.218138173930978d-2*t271*t3*t57*t71*t74*t9
            t522 = 3.109069086965489d-2*t207*t9
            t523 = -9.526483764256461d+0*t1*t10*t115*t13*t19*t41*t42
     +        *t70*t73
            t524 = 4.763241882128231d+0*t1*t13*t19*t41*t42*t51*t73
            t525 = 7.939887128782348d+1*t10*t110*t19*t203*t24*t265*t41
     +        *t51
            t526 = 9.723631383283067d+0*t10*t17*t20*t23*t41*t51*t53*t73
            t527 = -3.889452553313227d+1*t17*t20*t23*t26*t41*t53
            t528 = -7.778905106626454d+1*t20*t23*t24*t41*t53
            t529 = t30*(3.109069086965489d-2*t3*(
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t19*t269*t56*t70
     +        *t73+4.0827787561099116d+0*t1*t10*t115*t12*t13*t50*t70
     +        *t73-8.165557512219823d+0*t1*t10*t115*t12*t19*t39*t40
     +        *t70*t73+8.334541185671199d+0*t10*t115*t17*t20*t22*t23
     +        *t56*t73+4.0827787561099116d+0*t1*t115*t12*t13*t19*t56
     +        *t73+8.334541185671199d+0*t10*t17*t20*t22*t39*t51*t52
     +        *t73-2.0413893780549558d+0*t1*t12*t13*t50*t51*t73
     +        +4.0827787561099116d+0*t1*t12*t19*t39*t40*t51*t73
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t169*t51*t73
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t171*t19*t73
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t182*t51*t70
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t184*t19*t70
     +        -1.6669082371342397d+1*t17*t20*t22*t23*t51*t56
     +        +6.805617538956298d+1*t10*t110*t111*t115*t19*t203*t24
     +        *t56-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -3.3338164742684795d+1*t15*t20*t22*t23*t48*t55+t528
     +        -6.667632948536959d+1*t20*t22*t24*t39*t52)+t527+t526
     +        +t525+t524+t523-3.3338164742684795d+1*t17*t20*t22*t26
     +        *t39*t52+t519+t518+t517+t516+t515+t514+t511+t510
     +        -3.402808769478149d+1*t10*t110*t111*t203*t24*t50*t51
     +        +6.805617538956298d+1*t10*t111*t19*t203*t24*t264*t39*t51
     +        -4.1672705928355996d+0*t10*t17*t184*t20*t22*t23*t51
     +        -2.0413893780549558d+0*t1*t12*t13*t184*t19*t51+t508+t507
     +        +2.0413893780549558d+0*t1*t12*t13*t182*t26)*t9
     +        +9.327207260896466d-2*t207*t3*t39*t58+t522)
            t530 = -3.109069086965489d-2*t207*t3*t57*t9
            t531 = 3.109069086965489d-2*t211*t9
            t532 = -9.526483764256461d+0*t1*t10*t115*t13*t19*t210*t41
     +        *t42
            t533 = 1.9052967528512924d+1*t1*t13*t19*t41*t42*t51*t73
            t534 = 1.5879774257564697d+2*t10*t110*t19*t203*t24*t265
     +        *t41*t51
            t535 = 3.889452553313227d+1*t10*t17*t20*t23*t41*t51*t53*t73
            t536 = -7.778905106626454d+1*t17*t20*t23*t26*t41*t53
            t537 = -1.5557810213252907d+2*t20*t23*t24*t41*t53
            t538 = 8.165557512219823d+0*t1*t12*t19*t40*t51*t63*t70
            t539 = 1.6669082371342397d+1*t10*t17*t20*t22*t51*t52*t63*t70
            t540 = 4.0827787561099116d+0*t1*t12*t13*t189*t26
            t541 = -4.0827787561099116d+0*t1*t12*t13*t19*t190*t51
            t542 = -8.334541185671199d+0*t10*t17*t190*t20*t22*t23*t51
            t543 = -4.0827787561099116d+0*t1*t12*t13*t51*t67*t70
            t544 = 8.165557512219823d+0*t1*t115*t12*t13*t19*t68*t70
            t545 = 1.6669082371342397d+1*t10*t115*t17*t20*t22*t23*t68
     +        *t70
            t546 = 1.0d+0*(t119*(-3.109069086965489d-2*t191*t3*t71*t9
     +        -3.109069086965489d-2*t205*t3*t69*t9-t192*t71)+t30*
     +        (3.109069086965489d-2*t3*(-4.0827787561099116d+0*t1*t10
     +        *t12*t13*t189*t51*t70+8.165557512219823d+0*t1*t10*t115
     +        *t12*t13*t19*t190*t70-8.334541185671199d+0*t17*t20*t22
     +        *t23*t51*t68-1.2248336268329735d+1*t1*t10*t12*t13*t19
     +        *t201*t269*t68+3.402808769478149d+1*t10*t110*t111*t115
     +        *t19*t203*t24*t68-1.7014043847390745d+1*t10*t110*t111
     +        *t203*t24*t51*t67+4.0827787561099116d+0*t1*t10*t115*t12
     +        *t13*t201*t67-2.0413893780549558d+0*t1*t10*t12*t13*t19
     +        *t51*(-1.6669082371342397d+1*t15*t20*t22*t23*t55*t65
     +        -3.3338164742684795d+1*t20*t22*t24*t52*t63+t513)
     +        -1.6669082371342397d+1*t17*t20*t22*t26*t52*t63
     +        +3.402808769478149d+1*t10*t111*t19*t203*t24*t264*t51*t63
     +        -8.165557512219823d+0*t1*t10*t115*t12*t19*t201*t40*t63
     +        +t545+t544+t543+t542+t541+t540+t539+t538+t512+t511+t509
     +        +t508+t506)*t9+9.327207260896466d-2*t205*t3*t58*t63+t505
     +        )+6.218138173930978d-2*t200*t271*t3*t69*t9)*wght
            t547 = 6.218138173930978d-2*t271*t3*t69*t71*t74*t9
            t548 = t30*(3.109069086965489d-2*t3*(
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t19*t269*t68*t70
     +        *t73+4.0827787561099116d+0*t1*t10*t115*t12*t13*t67*t70
     +        *t73-8.165557512219823d+0*t1*t10*t115*t12*t19*t40*t63
     +        *t70*t73+8.334541185671199d+0*t10*t115*t17*t20*t22*t23
     +        *t68*t73+4.0827787561099116d+0*t1*t115*t12*t13*t19*t68
     +        *t73-2.0413893780549558d+0*t1*t12*t13*t51*t67*t73
     +        +8.334541185671199d+0*t10*t17*t20*t22*t51*t52*t63*t73
     +        +4.0827787561099116d+0*t1*t12*t19*t40*t51*t63*t73
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t189*t51*t73
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t190*t73
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t196*t51*t70
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t197*t70
     +        -1.6669082371342397d+1*t17*t20*t22*t23*t51*t68
     +        +6.805617538956298d+1*t10*t110*t111*t115*t19*t203*t24
     +        *t68-3.402808769478149d+1*t10*t110*t111*t203*t24*t51*t67
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -3.3338164742684795d+1*t15*t20*t22*t23*t55*t65
     +        -6.667632948536959d+1*t20*t22*t24*t52*t63+t528)
     +        -3.3338164742684795d+1*t17*t20*t22*t26*t52*t63
     +        +6.805617538956298d+1*t10*t111*t19*t203*t24*t264*t51*t63
     +        +t545+t544+t543+t542+t541+t540+t539+t538+t527+t526+t525
     +        +t524+t523+t511-4.1672705928355996d+0*t10*t17*t197*t20
     +        *t22*t23*t51-2.0413893780549558d+0*t1*t12*t13*t19*t197
     +        *t51+t508+2.0413893780549558d+0*t1*t12*t13*t196*t26)*t9
     +        +9.327207260896466d-2*t207*t3*t58*t63+t522)
            t549 = -3.109069086965489d-2*t207*t3*t69*t9
            t550 = t1**4
            t551 = 1/t11**8.0d+0
            t552 = 1/t8**8.0d+0
            t553 = 6.218138173930978d-2*t271*t3*t71**3*t9*wght
     +        -9.327207260896466d-2*t119*t205*t3*t71*t9*wght
     +        +3.109069086965489d-2*t3*t30*(-1.2248336268329735d+1*t1
     +        *t10*t12*t13*t19*t269*t70**3-2.50036235570136d+1*t17*t20
     +        *t22*t23*t51*t70+1.0208426308434447d+2*t10*t110*t111
     +        *t115*t19*t203*t24*t70-1.0419686516347423d+2*t10*t51*t55
     +        *t550*t551*t552-5.104213154217223d+1*t110*t111*t19*t203
     +        *t24*t51+2.50036235570136d+1*t10*t115*t17*t20*t201*t22
     +        *t23+1.2248336268329735d+1*t1*t115*t12*t13*t19*t201)*t9
     +        *wght
            t554 = -1.2248336268329735d+1*t1*t10*t12*t13*t19*t201*t269
     +        *t73
            t555 = 8.165557512219823d+0*t1*t115*t12*t13*t19*t201
            t556 = 8.165557512219823d+0*t1*t115*t12*t13*t19*t70*t73
            t557 = 1.3611235077912595d+2*t10*t110*t111*t115*t19*t203
     +        *t24*t70
            t558 = 3.402808769478149d+1*t10*t110*t111*t115*t19*t203
     +        *t24*t73
            t559 = -1.0208426308434447d+2*t110*t111*t19*t203*t24*t51
            t560 = 1.6669082371342397d+1*t10*t115*t17*t20*t201*t22*t23
            t561 = 1.6669082371342397d+1*t10*t115*t17*t20*t22*t23*t70
     +        *t73
            t562 = -3.3338164742684795d+1*t17*t20*t22*t23*t51*t70
            t563 = -8.334541185671199d+0*t17*t20*t22*t23*t51*t73
            t564 = 6.218138173930978d-2*t200*t271*t3*t74*t9*wght
            t565 = t119*(-3.109069086965489d-2*t205*t3*t74*t9*wght
     +        -6.218138173930978d-2*t207*t3*t71*t9*wght)
            t566 = 6.218138173930978d-2*t209*t271*t3*t71*t9*wght
            t567 = -1.2248336268329735d+1*t1*t10*t12*t13*t19*t210*t269
     +        *t70
            t568 = 4.0827787561099116d+0*t1*t115*t12*t13*t19*t210
            t569 = 1.6331115024439646d+1*t1*t115*t12*t13*t19*t70*t73
            t570 = 1.3611235077912595d+2*t10*t110*t111*t115*t19*t203
     +        *t24*t73
            t571 = -2.0416852616868894d+2*t110*t111*t19*t203*t24*t51
            t572 = 8.334541185671199d+0*t10*t115*t17*t20*t210*t22*t23
            t573 = 3.3338164742684795d+1*t10*t115*t17*t20*t22*t23*t70
     +        *t73
            t574 = -3.3338164742684795d+1*t17*t20*t22*t23*t51*t73
            t575 = t119*(-6.218138173930978d-2*t207*t3*t74*t9*wght
     +        -3.109069086965489d-2*t211*t3*t71*t9*wght)
            t576 = 3.109069086965489d-2*t3*t30*(t563+t562+t561+t560
     +        +t559+t558+t557+t556+t555+t554-2.0839373032694844d+2*t10
     +        *t51*t55*t550*t551*t552)*t9*wght+t565+t564
            fnc(iq) = 1.0d+0*(3.109069086965489d-2*t28*t3*t9+s1f)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(3.109069086965489d-2*t3*t30*t57
     +        *t9+t59+9.327207260896466d-2*t28*t3*t39*t58+t29)*wght
     +        +Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(3.109069086965489d-2*t3*t30*t69
     +        *t9+9.327207260896466d-2*t28*t3*t58*t63+t60+t59)*wght
     +        +Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t72+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t3*t30*t74*t9*wght
     +        +Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = t72+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-3.109069086965489d-2*t119
     +        *t120*t3*t9+9.327207260896466d-2*t28*t3*t58*t88
     +        +1.8654414521792934d-1*t28*t3*t76*t8+t75
     +        +1.8654414521792934d-1*t28*t39*t58+t118*t30)*wght+Amat2
     +        (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(-3.109069086965489d-2*t119*t3
     +        *t57*t69*t9+1.8654414521792934d-1*t28*t3*t39*t63*t8
     +        +9.327207260896466d-2*t28*t58*t63+9.327207260896466d-2
     +        *t28*t39*t58+9.327207260896466d-2*t122*t28*t3*t58+t138
     +        *t30+t121)*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-3.109069086965489d-2*t119
     +        *t161*t3*t9+1.8654414521792934d-1*t140*t28*t3*t8
     +        +1.8654414521792934d-1*t28*t58*t63+9.327207260896466d-2
     +        *t145*t28*t3*t58+t160*t30+t139)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t174+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(t186*t30
     +        -3.109069086965489d-2*t119*t3*t57*t74*t9)*wght+Cmat2(iq,
     +        D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t174+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t193+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = 1.0d+0*(t199*t30
     +        -3.109069086965489d-2*t119*t3*t69*t74*t9)*wght+Cmat2(iq,
     +        D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t193+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t206+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t208+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t206+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 3.109069086965489d-2*t211*t3*t30*t9
     +        *wght-3.109069086965489d-2*t119*t209*t3*t9*wght+Cmat2
     +        (iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = t208+Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t206+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -3.500507297981904d+2*t109*t20*t21*t23*t41*t53*t93*t94
     +        -3.0004348268416314d+2*t109*t20*t21*t22*t39*t52*t93*t94
     +        +2.50036235570136d+1*t109*t20*t21*t22*t23*t233*t93*t94
     +        +7.50108706710408d+1*t105*t109*t20*t21*t22*t23*t48*t94
     +        +2.50036235570136d+1*t109*t20*t21*t22*t228*t23*t94
     +        +1.1668357659939681d+2*t15*t20*t21*t23*t41*t53*t55*t93
     +        +1.0001449422805439d+2*t15*t20*t21*t22*t39*t52*t55*t93
     +        +1.944726276656614d+1*t15*t20*t21*t23*t48*t53*t55*t91
     +        +3.889452553313228d+1*t20*t21*t24*t39*t52*t53*t91
     +        -6.806541968298148d+2*t111*t15*t20*t21*t23*t48*t55*t89
     +        -1.3613083936596296d+3*t111*t20*t21*t24*t39*t52*t89
     +        +1.0001449422805439d+2*t15*t20*t21*t22*t48*t52*t55*t88
     +        +2.3336715319879361d+2*t20*t21*t24*t41*t52*t53*t88
     +        +2.50036235570136d+2*t110*t20*t21*t22*t24*t39*t88
     +        -5.00072471140272d+2*t110*t15*t20*t21*t22*t48*t55*t76
     +        -1.166835765993968d+3*t110*t20*t21*t24*t41*t53*t76
     +        -9.334686127951745d+2*t15*t20*t21*t39*t41*t48*t52*t53
     +        *t55+1.1668357659939681d+2*t105*t15*t20*t21*t23*t41*t53
     +        *t55+1.0001449422805439d+2*t105*t15*t20*t21*t22*t39*t52
     +        *t55-2.50036235570136d+1*t105*t15*t20*t21*t22*t23*t48
     +        *t55-8.334541185671199d+0*t15*t20*t21*t22*t23*t258*t55
     +        -8.334541185671199d+0*t15*t20*t21*t22*t228*t23*t55
     +        -1.6669082371342397d+1*t20*t21*t22*t221*t24*t52+t268
     +        +t267+t266-5.00072471140272d+2*t20*t21*t213*t22*t24*t264
     +        -1.0001449422805439d+2*t20*t21*t22*t228*t229*t23*t263
     +        +t262+t261+t260+t259+t255+t254+t253+t252+t251+t250+t249
     +        +t248+t247+t246+t245+t244+t243+t242+t241+t240+t239+t238
     +        +t237+t236+t235+t234+t231+t230)+2.3816209410641168d+0*t1
     +        *t10*t13*t19*t42*t51*t56*t91-2.3816209410641168d+0*t1
     +        *t10*t13*t26*t42*t50*t91+4.7632418821282335d+0*t1*t10
     +        *t19*t26*t39*t40*t42*t91-5.001403976234643d+1*t1*t10*t13
     +        *t19*t22*t51*t56*t89+5.001403976234643d+1*t1*t10*t13*t22
     +        *t26*t50*t89-1.0002807952469285d+2*t1*t10*t19*t22*t26
     +        *t39*t40*t89+1.2248336268329735d+1*t1*t10*t12*t19*t40
     +        *t51*t56*t88-1.2248336268329735d+1*t1*t10*t12*t26*t40
     +        *t50*t88+2.8579451292769387d+1*t1*t10*t19*t26*t40*t41
     +        *t42*t88+3.6745008804989204d+1*t1*t10*t12*t19*t23*t26
     +        *t39*t88-3.6745008804989204d+1*t1*t10*t12*t19*t23*t51
     +        *t56*t76+3.6745008804989204d+1*t1*t10*t12*t23*t26*t50
     +        *t76-8.573835387830815d+1*t1*t10*t19*t23*t26*t41*t42*t76
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t19*t269*t56**3
     +        +2.8579451292769387d+1*t1*t10*t13*t41*t42*t50*t51*t56
     +        +2.449667253665947d+1*t1*t10*t12*t39*t40*t50*t51*t56
     +        -5.715890258553877d+1*t1*t10*t19*t39*t40*t41*t42*t51*t56
     +        -6.124168134164868d+0*t1*t10*t108*t12*t13*t51*t56
     +        +1.2248336268329735d+1*t1*t10*t114*t115*t12*t13*t19*t56
     +        -4.899334507331894d+1*t1*t10*t12*t19*t213*t26*t52
     +        -6.124168134164868d+0*t1*t10*t114*t12*t13*t50*t51
     +        +1.4289725646384693d+1*t1*t10*t114*t13*t19*t41*t42*t51
     +        +1.2248336268329735d+1*t1*t10*t114*t12*t19*t39*t40*t51
     +        +5.715890258553877d+1*t1*t10*t26*t39*t40*t41*t42*t50
     +        +1.2248336268329735d+1*t1*t10*t115*t116*t12*t13*t50
     +        -1.4289725646384693d+1*t1*t10*t108*t13*t26*t41*t42
     +        -2.8579451292769387d+1*t1*t10*t115*t116*t13*t19*t41*t42
     +        -1.2248336268329735d+1*t1*t10*t108*t12*t26*t39*t40
     +        -2.449667253665947d+1*t1*t10*t115*t116*t12*t19*t39*t40
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t221*t26*t40
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t262+t261
     +        +t260+t259+t255+t254+t253+t252+t251+t250+t249+t248+t247
     +        +t246+t245+t244+t243+t242+t241+t240+t239+t238+t237+t236
     +        +t235+t234+t231+t230)+t227+t225+t223)
     +        +9.327207260896466d-2*t117*t9+2.7981621782689403d-1*t3
     +        *t57*t58*t88+5.596324356537881d-1*t3*t57*t76*t8
     +        +5.596324356537881d-1*t39*t57*t58+2.7981621782689403d-1
     +        *t117*t3*t39*t58)+t119*(-6.218138173930978d-2*t117*t3
     +        *t57*t9-9.327207260896466d-2*t120*t3*t39*t58-t118*t57
     +        +t270)+6.218138173930978d-2*t271*t3*t57**3*t9
     +        +5.596324356537881d-1*t28*t3*t39*t8*t88
     +        +2.7981621782689403d-1*t28*t58*t88+5.596324356537881d-1
     +        *t28*t76*t8+9.327207260896466d-2*t221*t28*t3*t58
     +        +1.8654414521792934d-1*t213*t28*t3+t212)*wght+Amat3(iq,
     +        D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (2.50036235570136d+1*t109*t20*t21*t22*t23*t65*t93*t94
     +        -1.0001449422805439d+2*t109*t20*t21*t22*t52*t63*t93*t94
     +        -1.1668357659939681d+2*t109*t20*t21*t23*t41*t53*t93*t94
     +        +2.50036235570136d+1*t109*t20*t21*t22*t23*t290*t93*t94
     +        -2.0002898845610878d+2*t109*t20*t21*t22*t39*t48*t52*t65
     +        *t94+2.50036235570136d+1*t105*t109*t20*t21*t22*t23*t65
     +        *t94+5.00072471140272d+1*t109*t130*t20*t21*t22*t23*t48
     +        *t94-8.334541185671199d+0*t15*t20*t21*t22*t23*t55*t65
     +        *t93-1.0001449422805439d+2*t20*t21*t22*t229*t23*t263*t65
     +        *t93+3.3338164742684795d+1*t15*t20*t21*t22*t52*t55*t63
     +        *t93+3.889452553313227d+1*t15*t20*t21*t23*t41*t53*t55
     +        *t93+6.482420922188714d+0*t15*t20*t21*t23*t53*t55*t65
     +        *t91+1.2964841844377428d+1*t20*t21*t24*t52*t53*t63*t91
     +        +1.2964841844377428d+1*t15*t20*t21*t23*t48*t53*t55*t91
     +        +2.5929683688754857d+1*t20*t21*t24*t39*t52*t53*t91
     +        -2.268847322766049d+2*t111*t15*t20*t21*t23*t55*t65*t89
     +        -4.537694645532098d+2*t111*t20*t21*t24*t52*t63*t89
     +        -4.537694645532098d+2*t111*t15*t20*t21*t23*t48*t55*t89
     +        -9.075389291064196d+2*t111*t20*t21*t24*t39*t52*t89
     +        +3.3338164742684795d+1*t15*t20*t21*t22*t52*t55*t65*t88
     +        +8.334541185671199d+1*t110*t20*t21*t22*t24*t63*t88
     +        +7.778905106626454d+1*t20*t21*t24*t41*t52*t53*t88
     +        -1.6669082371342397d+2*t110*t15*t20*t21*t22*t55*t65*t76
     +        -5.00072471140272d+2*t20*t21*t22*t24*t264*t63*t76
     +        -3.889452553313226d+2*t110*t20*t21*t24*t41*t53*t76
     +        +6.667632948536959d+1*t15*t20*t21*t22*t39*t48*t52*t55
     +        *t65-8.334541185671199d+0*t105*t15*t20*t21*t22*t23*t55
     +        *t65+3.3338164742684795d+1*t105*t15*t20*t21*t22*t52*t55
     +        *t63-3.3338164742684795d+2*t110*t15*t20*t21*t22*t39*t48
     +        *t55*t63-3.1115620426505813d+2*t15*t20*t21*t39*t41*t48
     +        *t52*t53*t55+3.889452553313227d+1*t105*t15*t20*t21*t23
     +        *t41*t53*t55+6.667632948536959d+1*t122*t15*t20*t21*t22
     +        *t48*t52*t55+6.667632948536959d+1*t130*t15*t20*t21*t22
     +        *t39*t52*t55-1.6669082371342397d+1*t130*t15*t20*t21*t22
     +        *t23*t48*t55-8.334541185671199d+0*t15*t20*t21*t22*t23
     +        *t326*t55-1.6669082371342397d+1*t20*t21*t22*t24*t278*t52
     +        +1.6669082371342397d+2*t110*t122*t20*t21*t22*t24*t39
     +        +t345+t344+t343+t342+t341+t340+t339+t338+t337+t336+t335
     +        +t334+t333+t332+t331+t330+t329+t328+t327+t320+t319+t318
     +        +t317+t316+t315+t314+t313+t312+t311+t310+t309+t308+t307
     +        +t306+t305+t304+t303+t302+t301+t300+t299+t298+t297+t296
     +        +t295+t294+t293+t292+t291+t289+t288+t268+t267+t266+t247
     +        +t246+t245)+7.938736470213722d-1*t1*t10*t13*t19*t42*t51
     +        *t68*t91-7.938736470213722d-1*t1*t10*t13*t26*t42*t67*t91
     +        +1.5877472940427445d+0*t1*t10*t19*t26*t40*t42*t63*t91
     +        +1.5877472940427445d+0*t1*t10*t13*t19*t42*t51*t56*t91
     +        -1.5877472940427445d+0*t1*t10*t13*t26*t42*t50*t91
     +        +3.175494588085489d+0*t1*t10*t19*t26*t39*t40*t42*t91
     +        -1.667134658744881d+1*t1*t10*t13*t19*t22*t51*t68*t89
     +        +1.667134658744881d+1*t1*t10*t13*t22*t26*t67*t89
     +        -3.334269317489762d+1*t1*t10*t19*t22*t26*t40*t63*t89
     +        -3.334269317489762d+1*t1*t10*t13*t19*t22*t51*t56*t89
     +        +3.334269317489762d+1*t1*t10*t13*t22*t26*t50*t89
     +        -6.668538634979524d+1*t1*t10*t19*t22*t26*t39*t40*t89
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t68*t88
     +        -4.0827787561099116d+0*t1*t10*t12*t26*t40*t67*t88
     +        +1.2248336268329735d+1*t1*t10*t12*t19*t23*t26*t63*t88
     +        +9.526483764256461d+0*t1*t10*t19*t26*t40*t41*t42*t88
     +        -1.2248336268329735d+1*t1*t10*t12*t19*t23*t51*t68*t76
     +        +1.2248336268329735d+1*t1*t10*t12*t23*t26*t67*t76
     +        -4.899334507331894d+1*t1*t10*t12*t19*t26*t52*t63*t76
     +        -2.8579451292769387d+1*t1*t10*t19*t23*t26*t41*t42*t76
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t50*t56*t68
     +        -1.6331115024439646d+1*t1*t10*t115*t12*t19*t39*t40*t56
     +        *t68+8.165557512219823d+0*t1*t10*t12*t39*t40*t50*t51*t68
     +        -2.0413893780549558d+0*t1*t10*t108*t12*t13*t51*t68
     +        -1.2248336268329735d+1*t1*t10*t116*t12*t13*t19*t269*t68
     +        +4.0827787561099116d+0*t1*t10*t114*t115*t12*t13*t19*t68
     +        +8.165557512219823d+0*t1*t10*t12*t39*t40*t51*t56*t67
     +        -2.0413893780549558d+0*t1*t10*t114*t12*t13*t51*t67
     +        +4.0827787561099116d+0*t1*t10*t115*t116*t12*t13*t67
     +        +8.165557512219823d+0*t1*t10*t12*t40*t50*t51*t56*t63
     +        -2.449667253665947d+1*t1*t10*t12*t19*t23*t39*t51*t56*t63
     +        +4.0827787561099116d+0*t1*t10*t114*t12*t19*t40*t51*t63
     +        +2.449667253665947d+1*t1*t10*t12*t23*t26*t39*t50*t63
     +        -4.0827787561099116d+0*t1*t10*t108*t12*t26*t40*t63
     +        -8.165557512219823d+0*t1*t10*t115*t116*t12*t19*t40*t63
     +        +9.526483764256461d+0*t1*t10*t13*t41*t42*t50*t51*t56
     +        -1.9052967528512924d+1*t1*t10*t19*t39*t40*t41*t42*t51
     +        *t56+8.165557512219823d+0*t1*t10*t12*t122*t19*t40*t51
     +        *t56-4.0827787561099116d+0*t1*t10*t12*t13*t135*t51*t56
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t136*t19*t56
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t136*t50*t51
     +        +4.763241882128231d+0*t1*t10*t114*t13*t19*t41*t42*t51
     +        +8.165557512219823d+0*t1*t10*t12*t136*t19*t39*t40*t51
     +        +1.9052967528512924d+1*t1*t10*t26*t39*t40*t41*t42*t50
     +        -8.165557512219823d+0*t1*t10*t12*t122*t26*t40*t50
     +        -4.763241882128231d+0*t1*t10*t108*t13*t26*t41*t42
     +        -9.526483764256461d+0*t1*t10*t115*t116*t13*t19*t41*t42
     +        -8.165557512219823d+0*t1*t10*t12*t135*t26*t39*t40
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t278*t40
     +        +2.449667253665947d+1*t1*t10*t12*t122*t19*t23*t26*t39
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t338+t337
     +        +t336+t335+t334+t333+t332+t331+t330+t329+t328+t327+t320
     +        +t319+t318+t317+t316+t315+t314+t313+t312+t311+t310+t309
     +        +t308+t307+t306+t305+t304+t303+t302+t301+t300+t299+t298
     +        +t297+t296+t295+t294+t293+t292+t291+t289+t288+t247+t246
     +        +t245)+t287+t286+t285+t284+t283+t282+t281+t280+t279+t277
     +        +t276+t227+t225+t223)+3.109069086965489d-2*t117*t9
     +        +9.327207260896466d-2*t3*t58*t69*t88
     +        +1.8654414521792934d-1*t3*t69*t76*t8
     +        +3.7308829043585867d-1*t3*t39*t57*t63*t8
     +        +9.327207260896466d-2*t117*t3*t58*t63
     +        +1.8654414521792934d-1*t39*t57*t58+1.8654414521792934d-1
     +        *t122*t3*t57*t58+1.8654414521792934d-1*t137*t3*t39*t58
     +        +t275+t274+t273)+t119*(-6.218138173930978d-2*t137*t3*t57
     +        *t9-t118*t69-9.327207260896466d-2*t120*t3*t58*t63+t270)
     +        +6.218138173930978d-2*t120*t271*t3*t69*t9
     +        +1.8654414521792934d-1*t28*t3*t63*t8*t88
     +        +9.327207260896466d-2*t28*t58*t88+1.8654414521792934d-1
     +        *t28*t76*t8+3.7308829043585867d-1*t122*t28*t3*t39*t8
     +        +1.8654414521792934d-1*t28*t3*t63*t76
     +        +9.327207260896466d-2*t278*t28*t3*t58+t347+t346+t272)
     +        *wght+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -2.0002898845610878d+2*t109*t20*t21*t22*t48*t52*t63*t65
     +        *t94+2.50036235570136d+1*t109*t20*t21*t22*t23*t290*t48
     +        *t65*t94+5.00072471140272d+1*t109*t130*t20*t21*t22*t23
     +        *t65*t94-1.1668357659939681d+2*t109*t146*t20*t21*t23*t41
     +        *t53*t94-1.0001449422805439d+2*t109*t146*t20*t21*t22*t39
     +        *t52*t94+2.50036235570136d+1*t109*t153*t20*t21*t22*t23
     +        *t48*t94+2.50036235570136d+1*t109*t146*t20*t21*t22*t23
     +        *t48*t94+1.2964841844377428d+1*t15*t20*t21*t23*t53*t55
     +        *t65*t91+2.5929683688754857d+1*t20*t21*t24*t52*t53*t63
     +        *t91+6.482420922188714d+0*t15*t20*t21*t23*t48*t53*t55
     +        *t91+1.2964841844377428d+1*t20*t21*t24*t39*t52*t53*t91
     +        -4.537694645532098d+2*t111*t15*t20*t21*t23*t55*t65*t89
     +        -9.075389291064196d+2*t111*t20*t21*t24*t52*t63*t89
     +        -2.268847322766049d+2*t111*t15*t20*t21*t23*t48*t55*t89
     +        -4.537694645532098d+2*t111*t20*t21*t24*t39*t52*t89
     +        -3.1115620426505813d+2*t15*t20*t21*t41*t52*t53*t55*t63
     +        *t65+6.667632948536959d+1*t15*t20*t21*t22*t48*t52*t55
     +        *t63*t65-3.3338164742684795d+2*t110*t15*t20*t21*t22*t39
     +        *t55*t63*t65+6.667632948536959d+1*t122*t15*t20*t21*t22
     +        *t52*t55*t65-1.6669082371342397d+1*t130*t15*t20*t21*t22
     +        *t23*t55*t65+6.667632948536959d+1*t130*t15*t20*t21*t22
     +        *t52*t55*t63+1.6669082371342397d+2*t110*t122*t20*t21*t22
     +        *t24*t63+3.889452553313227d+1*t15*t153*t20*t21*t23*t41
     +        *t53*t55+3.889452553313227d+1*t146*t15*t20*t21*t23*t41
     +        *t53*t55+3.3338164742684795d+1*t145*t15*t20*t21*t22*t48
     +        *t52*t55+3.3338164742684795d+1*t15*t153*t20*t21*t22*t39
     +        *t52*t55+3.3338164742684795d+1*t146*t15*t20*t21*t22*t39
     +        *t52*t55-8.334541185671199d+0*t15*t153*t20*t21*t22*t23
     +        *t48*t55-8.334541185671199d+0*t146*t15*t20*t21*t22*t23
     +        *t48*t55-1.6669082371342397d+2*t110*t140*t15*t20*t21*t22
     +        *t48*t55-8.334541185671199d+0*t15*t20*t21*t22*t23*t373
     +        *t55+7.778905106626454d+1*t145*t20*t21*t24*t41*t52*t53
     +        -3.889452553313226d+2*t110*t140*t20*t21*t24*t41*t53
     +        -1.6669082371342397d+1*t20*t21*t22*t24*t349*t52
     +        -1.0001449422805439d+2*t146*t20*t21*t22*t229*t23*t263
     +        *t48-5.00072471140272d+2*t140*t20*t21*t22*t24*t264*t39
     +        +8.334541185671199d+1*t110*t145*t20*t21*t22*t24*t39+t387
     +        +t386+t385+t384+t383+t382+t381+t380+t379+t378+t377+t376
     +        +t375+t374+t372+t371+t370+t369+t368+t367+t366+t365+t364
     +        +t363+t362+t361+t360+t359+t358+t357+t356+t355+t354+t353
     +        +t352+t351+t350+t345+t344+t343+t342+t341+t340+t339+t333
     +        +t330+t318+t312+t303+t301+t293+t268+t267+t266+t247+t246
     +        +t245)+1.5877472940427445d+0*t1*t10*t13*t19*t42*t51*t68
     +        *t91-1.5877472940427445d+0*t1*t10*t13*t26*t42*t67*t91
     +        +3.175494588085489d+0*t1*t10*t19*t26*t40*t42*t63*t91
     +        +7.938736470213722d-1*t1*t10*t13*t19*t42*t51*t56*t91
     +        -7.938736470213722d-1*t1*t10*t13*t26*t42*t50*t91
     +        +1.5877472940427445d+0*t1*t10*t19*t26*t39*t40*t42*t91
     +        -3.334269317489762d+1*t1*t10*t13*t19*t22*t51*t68*t89
     +        +3.334269317489762d+1*t1*t10*t13*t22*t26*t67*t89
     +        -6.668538634979524d+1*t1*t10*t19*t22*t26*t40*t63*t89
     +        -1.667134658744881d+1*t1*t10*t13*t19*t22*t51*t56*t89
     +        +1.667134658744881d+1*t1*t10*t13*t22*t26*t50*t89
     +        -3.334269317489762d+1*t1*t10*t19*t22*t26*t39*t40*t89
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t56*t67*t68
     +        +9.526483764256461d+0*t1*t10*t13*t41*t42*t51*t67*t68
     +        +8.165557512219823d+0*t1*t10*t12*t39*t40*t51*t67*t68
     +        -1.6331115024439646d+1*t1*t10*t115*t12*t19*t40*t56*t63
     +        *t68+8.165557512219823d+0*t1*t10*t12*t40*t50*t51*t63*t68
     +        -1.9052967528512924d+1*t1*t10*t19*t40*t41*t42*t51*t63
     +        *t68-2.449667253665947d+1*t1*t10*t12*t19*t23*t39*t51*t63
     +        *t68+8.165557512219823d+0*t1*t10*t12*t122*t19*t40*t51
     +        *t68-4.0827787561099116d+0*t1*t10*t12*t13*t135*t51*t68
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t136*t19*t68
     +        +8.165557512219823d+0*t1*t10*t12*t40*t51*t56*t63*t67
     +        +1.9052967528512924d+1*t1*t10*t26*t40*t41*t42*t63*t67
     +        +2.449667253665947d+1*t1*t10*t12*t23*t26*t39*t63*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t136*t51*t67
     +        -8.165557512219823d+0*t1*t10*t12*t122*t26*t40*t67
     +        +8.165557512219823d+0*t1*t10*t12*t136*t19*t40*t51*t63
     +        -8.165557512219823d+0*t1*t10*t12*t135*t26*t40*t63
     +        +2.449667253665947d+1*t1*t10*t12*t122*t19*t23*t26*t63
     +        +4.0827787561099116d+0*t1*t10*t12*t145*t19*t40*t51*t56
     +        -1.2248336268329735d+1*t1*t10*t12*t140*t19*t23*t51*t56
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t156*t51*t56
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t158*t19*t269*t56
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t157*t19*t56
     +        -4.899334507331894d+1*t1*t10*t12*t140*t19*t26*t39*t52
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t157*t50*t51
     +        +4.763241882128231d+0*t1*t10*t13*t157*t19*t41*t42*t51
     +        +4.0827787561099116d+0*t1*t10*t12*t157*t19*t39*t40*t51
     +        -4.0827787561099116d+0*t1*t10*t12*t145*t26*t40*t50
     +        +1.2248336268329735d+1*t1*t10*t12*t140*t23*t26*t50
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t158*t50
     +        +9.526483764256461d+0*t1*t10*t145*t19*t26*t40*t41*t42
     +        -2.8579451292769387d+1*t1*t10*t140*t19*t23*t26*t41*t42
     +        -4.763241882128231d+0*t1*t10*t13*t156*t26*t41*t42
     +        -9.526483764256461d+0*t1*t10*t115*t13*t158*t19*t41*t42
     +        -4.0827787561099116d+0*t1*t10*t12*t156*t26*t39*t40
     +        -8.165557512219823d+0*t1*t10*t115*t12*t158*t19*t39*t40
     +        -4.0827787561099116d+0*t1*t10*t12*t19*t26*t349*t40
     +        +1.2248336268329735d+1*t1*t10*t12*t145*t19*t23*t26*t39
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t387+t386
     +        +t385+t384+t383+t382+t381+t380+t379+t378+t377+t376+t375
     +        +t374+t372+t371+t370+t369+t368+t367+t366+t365+t364+t363
     +        +t362+t361+t360+t359+t358+t357+t356+t355+t354+t353+t352
     +        +t351+t350+t333+t330+t318+t312+t303+t301+t293+t247+t246
     +        +t245)+t287+t286+t285+t284+t283+t282+t281+t280+t279+t277
     +        +t276+t227+t225+t223)+3.109069086965489d-2*t159*t9
     +        +3.7308829043585867d-1*t3*t39*t63*t69*t8
     +        +1.8654414521792934d-1*t140*t3*t57*t8
     +        +1.8654414521792934d-1*t58*t63*t69+1.8654414521792934d-1
     +        *t122*t3*t58*t69+1.8654414521792934d-1*t137*t3*t58*t63
     +        +9.327207260896466d-2*t145*t3*t57*t58
     +        +9.327207260896466d-2*t159*t3*t39*t58+t275+t274+t273)
     +        +t119*(-3.109069086965489d-2*t57*t69*t9
     +        -3.109069086965489d-2*t137*t3*t69*t9
     +        -3.109069086965489d-2*t159*t3*t57*t9
     +        -9.327207260896466d-2*t3*t57*t58*t63*t69-t138*t69)
     +        +6.218138173930978d-2*t161*t271*t3*t57*t9
     +        +3.7308829043585867d-1*t122*t28*t3*t63*t8
     +        +1.8654414521792934d-1*t145*t28*t3*t39*t8
     +        +1.8654414521792934d-1*t140*t28*t8+9.327207260896466d-2
     +        *t28*t3*t349*t58+9.327207260896466d-2*t145*t28*t58
     +        +1.8654414521792934d-1*t140*t28*t3*t39+t348+t347+t346)
     +        *wght+Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (7.50108706710408d+1*t109*t153*t20*t21*t22*t23*t65*t94
     +        -3.0004348268416314d+2*t109*t146*t20*t21*t22*t52*t63*t94
     +        -3.500507297981904d+2*t109*t146*t20*t21*t23*t41*t53*t94
     +        +2.50036235570136d+1*t109*t20*t21*t22*t23*t391*t94
     +        +2.50036235570136d+1*t109*t146*t20*t21*t22*t23*t290*t94
     +        +1.944726276656614d+1*t15*t20*t21*t23*t53*t55*t65*t91
     +        +3.889452553313228d+1*t20*t21*t24*t52*t53*t63*t91
     +        -6.806541968298148d+2*t111*t15*t20*t21*t23*t55*t65*t89
     +        -1.3613083936596296d+3*t111*t20*t21*t24*t52*t63*t89
     +        -9.334686127951745d+2*t15*t20*t21*t41*t52*t53*t55*t63
     +        *t65+1.0001449422805439d+2*t145*t15*t20*t21*t22*t52*t55
     +        *t65-2.50036235570136d+1*t15*t153*t20*t21*t22*t23*t55
     +        *t65-5.00072471140272d+2*t110*t140*t15*t20*t21*t22*t55
     +        *t65+1.0001449422805439d+2*t15*t153*t20*t21*t22*t52*t55
     +        *t63+1.0001449422805439d+2*t146*t15*t20*t21*t22*t52*t55
     +        *t63+2.50036235570136d+2*t110*t145*t20*t21*t22*t24*t63
     +        +1.1668357659939681d+2*t15*t153*t20*t21*t23*t41*t53*t55
     +        +1.1668357659939681d+2*t146*t15*t20*t21*t23*t41*t53*t55
     +        -8.334541185671199d+0*t15*t20*t21*t22*t23*t413*t55
     +        -8.334541185671199d+0*t15*t20*t21*t22*t23*t391*t55
     +        +2.3336715319879361d+2*t145*t20*t21*t24*t41*t52*t53
     +        -1.166835765993968d+3*t110*t140*t20*t21*t24*t41*t53
     +        -1.6669082371342397d+1*t20*t21*t22*t24*t390*t52+t417
     +        +t416+t415+t414+t412+t411+t410+t409+t408+t407+t406+t405
     +        +t404+t403+t402+t401+t400+t399+t398+t397+t396+t395+t394
     +        +t393+t392-1.0001449422805439d+2*t20*t21*t22*t229*t23
     +        *t263*t391-5.00072471140272d+2*t20*t21*t22*t24*t264*t389
     +        +t268+t267+t266+t247+t246+t245)+2.3816209410641168d+0*t1
     +        *t10*t13*t19*t42*t51*t68*t91-2.3816209410641168d+0*t1
     +        *t10*t13*t26*t42*t67*t91+4.7632418821282335d+0*t1*t10
     +        *t19*t26*t40*t42*t63*t91-5.001403976234643d+1*t1*t10*t13
     +        *t19*t22*t51*t68*t89+5.001403976234643d+1*t1*t10*t13*t22
     +        *t26*t67*t89-1.0002807952469285d+2*t1*t10*t19*t22*t26
     +        *t40*t63*t89-1.2248336268329735d+1*t1*t10*t12*t13*t19
     +        *t269*t68**3+2.449667253665947d+1*t1*t10*t12*t40*t51*t63
     +        *t67*t68+2.8579451292769387d+1*t1*t10*t13*t41*t42*t51
     +        *t67*t68-5.715890258553877d+1*t1*t10*t19*t40*t41*t42*t51
     +        *t63*t68+1.2248336268329735d+1*t1*t10*t12*t145*t19*t40
     +        *t51*t68-3.6745008804989204d+1*t1*t10*t12*t140*t19*t23
     +        *t51*t68-6.124168134164868d+0*t1*t10*t12*t13*t156*t51
     +        *t68+1.2248336268329735d+1*t1*t10*t115*t12*t13*t157*t19
     +        *t68+5.715890258553877d+1*t1*t10*t26*t40*t41*t42*t63*t67
     +        -6.124168134164868d+0*t1*t10*t12*t13*t157*t51*t67
     +        -1.2248336268329735d+1*t1*t10*t12*t145*t26*t40*t67
     +        +3.6745008804989204d+1*t1*t10*t12*t140*t23*t26*t67
     +        +1.2248336268329735d+1*t1*t10*t115*t12*t13*t158*t67
     +        +1.2248336268329735d+1*t1*t10*t12*t157*t19*t40*t51*t63
     +        -1.2248336268329735d+1*t1*t10*t12*t156*t26*t40*t63
     +        -2.449667253665947d+1*t1*t10*t115*t12*t158*t19*t40*t63
     +        +3.6745008804989204d+1*t1*t10*t12*t145*t19*t23*t26*t63
     +        -4.899334507331894d+1*t1*t10*t12*t19*t26*t389*t52
     +        +1.4289725646384693d+1*t1*t10*t13*t157*t19*t41*t42*t51
     +        +2.8579451292769387d+1*t1*t10*t145*t19*t26*t40*t41*t42
     +        -8.573835387830815d+1*t1*t10*t140*t19*t23*t26*t41*t42
     +        -1.4289725646384693d+1*t1*t10*t13*t156*t26*t41*t42
     +        -2.8579451292769387d+1*t1*t10*t115*t13*t158*t19*t41*t42
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t417+t416
     +        +t415+t414+t412+t411+t410+t409+t408+t407+t406+t405+t404
     +        +t403+t402+t401+t400+t399+t398+t397+t396+t395+t394+t393
     +        +t392+t247+t246+t245)-4.0827787561099116d+0*t1*t10*t12
     +        *t19*t26*t390*t40+t227+t225+t223)+9.327207260896466d-2
     +        *t159*t9+5.596324356537881d-1*t140*t3*t69*t8
     +        +5.596324356537881d-1*t58*t63*t69+2.7981621782689403d-1
     +        *t145*t3*t58*t69+2.7981621782689403d-1*t159*t3*t58*t63)
     +        +t119*(-6.218138173930978d-2*t159*t3*t69*t9
     +        -3.109069086965489d-2*t161*t9-t160*t69
     +        -9.327207260896466d-2*t161*t3*t58*t63)
     +        +6.218138173930978d-2*t271*t3*t69**3*t9
     +        +5.596324356537881d-1*t145*t28*t3*t63*t8
     +        +5.596324356537881d-1*t140*t28*t8+9.327207260896466d-2
     +        *t28*t3*t390*t58+2.7981621782689403d-1*t145*t28*t58
     +        +1.8654414521792934d-1*t28*t3*t389+t388)*wght+Amat3(iq,
     +        D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t443+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (1.0001449422805439d+2*t10*t109*t20*t22*t23*t93*t94
     +        -3.3338164742684795d+1*t10*t15*t20*t22*t23*t55*t93
     +        -6.667632948536959d+1*t10*t20*t22*t24*t52*t88
     +        +3.3338164742684795d+2*t10*t110*t20*t22*t24*t76
     +        +3.1115620426505813d+2*t10*t15*t20*t23*t41*t48*t53*t55
     +        +2.667053179414784d+2*t10*t15*t20*t22*t39*t48*t52*t55
     +        -3.3338164742684795d+1*t10*t105*t15*t20*t22*t23*t55
     +        +6.223124085301163d+2*t10*t20*t24*t39*t41*t52*t53+t461
     +        +t460+t459+t458+t457+t456+t455+t454+t453+t452+t451+t450)
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t40*t51*t73*t88
     +        -1.6669082371342397d+1*t10*t17*t20*t22*t26*t52*t88
     +        -8.165557512219823d+0*t1*t12*t19*t26*t40*t88
     +        -1.2248336268329735d+1*t1*t10*t12*t19*t23*t51*t73*t76
     +        +2.449667253665947d+1*t1*t12*t19*t23*t26*t76
     +        +5.00072471140272d+1*t10*t110*t17*t20*t22*t26*t76
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t50*t56*t73
     +        -1.9052967528512924d+1*t1*t10*t115*t13*t19*t41*t42*t56
     +        *t73-1.6331115024439646d+1*t1*t10*t115*t12*t19*t39*t40
     +        *t56*t73+9.526483764256461d+0*t1*t10*t13*t41*t42*t50*t51
     +        *t73+8.165557512219823d+0*t1*t10*t12*t39*t40*t50*t51*t73
     +        -1.9052967528512924d+1*t1*t10*t19*t39*t40*t41*t42*t51
     +        *t73-2.0413893780549558d+0*t1*t10*t108*t12*t13*t51*t73
     +        -1.2248336268329735d+1*t1*t10*t116*t12*t13*t19*t269*t73
     +        +4.0827787561099116d+0*t1*t10*t114*t115*t12*t13*t19*t73
     +        +3.889452553313227d+1*t10*t17*t20*t23*t41*t51*t53*t56
     +        +3.3338164742684795d+1*t10*t17*t20*t22*t39*t51*t52*t56
     +        -8.165557512219823d+0*t1*t12*t13*t50*t51*t56
     +        +1.9052967528512924d+1*t1*t13*t19*t41*t42*t51*t56
     +        +1.6331115024439646d+1*t1*t12*t19*t39*t40*t51*t56
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t182*t51*t56
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t184*t19*t56
     +        +7.778905106626454d+1*t10*t17*t20*t26*t39*t41*t52*t53
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t184*t50*t51
     +        +9.526483764256461d+0*t1*t10*t13*t184*t19*t41*t42*t51
     +        +8.165557512219823d+0*t1*t10*t12*t184*t19*t39*t40*t51
     +        -8.334541185671199d+0*t10*t114*t17*t20*t22*t23*t51
     +        -4.0827787561099116d+0*t1*t114*t12*t13*t19*t51
     +        -1.9052967528512924d+1*t1*t13*t26*t41*t42*t50
     +        -1.6331115024439646d+1*t1*t12*t26*t39*t40*t50
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t459+t458
     +        +t457+t456+t455+t454+t453+t452+t451+t450)+t449+t448+t447
     +        +t446+t445+t444+3.8105935057025847d+1*t1*t19*t26*t39*t40
     +        *t41*t42-9.526483764256461d+0*t1*t10*t13*t182*t26*t41
     +        *t42-8.165557512219823d+0*t1*t10*t12*t182*t26*t39*t40
     +        +4.0827787561099116d+0*t1*t108*t12*t13*t26
     +        +1.6669082371342397d+1*t10*t115*t116*t17*t20*t22*t23
     +        +8.165557512219823d+0*t1*t115*t116*t12*t13*t19)
     +        +6.218138173930978d-2*t185*t9+9.327207260896466d-2*t3
     +        *t58*t74*t88+1.8654414521792934d-1*t3*t74*t76*t8
     +        +1.8654414521792934d-1*t39*t58*t74+1.8654414521792934d-1
     +        *t185*t3*t39*t58)+t119*(-6.218138173930978d-2*t185*t3
     +        *t57*t9-t118*t74)+6.218138173930978d-2*t120*t271*t3*t74
     +        *t9)*wght+Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t443+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = t473+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (1.0001449422805439d+2*t10*t109*t20*t22*t23*t48*t65*t94
     +        +1.333526589707392d+2*t10*t15*t20*t22*t39*t52*t55*t65
     +        -3.3338164742684795d+1*t10*t15*t20*t22*t23*t48*t55*t65
     +        +1.333526589707392d+2*t10*t15*t20*t22*t48*t52*t55*t63
     +        +3.3338164742684795d+2*t10*t110*t20*t22*t24*t39*t63
     +        -3.3338164742684795d+1*t10*t130*t15*t20*t22*t23*t55
     +        -6.667632948536959d+1*t10*t122*t20*t22*t24*t52+t486+t485
     +        +t484+t483+t482+t481+t480+t479+t478+t477+t476+t461+t460
     +        +t455+t454+t442+t439+t436+t432)-1.2248336268329735d+1*t1
     +        *t10*t12*t13*t19*t269*t56*t68*t73+4.0827787561099116d+0
     +        *t1*t10*t115*t12*t13*t50*t68*t73-9.526483764256461d+0*t1
     +        *t10*t115*t13*t19*t41*t42*t68*t73-8.165557512219823d+0
     +        *t1*t10*t115*t12*t19*t39*t40*t68*t73
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t56*t67*t73
     +        +4.763241882128231d+0*t1*t10*t13*t41*t42*t51*t67*t73
     +        +4.0827787561099116d+0*t1*t10*t12*t39*t40*t51*t67*t73
     +        -8.165557512219823d+0*t1*t10*t115*t12*t19*t40*t56*t63
     +        *t73+4.0827787561099116d+0*t1*t10*t12*t40*t50*t51*t63
     +        *t73-9.526483764256461d+0*t1*t10*t19*t40*t41*t42*t51*t63
     +        *t73-1.2248336268329735d+1*t1*t10*t12*t19*t23*t39*t51
     +        *t63*t73-9.526483764256461d+0*t1*t10*t115*t13*t19*t41
     +        *t42*t56*t73+4.763241882128231d+0*t1*t10*t13*t41*t42*t50
     +        *t51*t73-9.526483764256461d+0*t1*t10*t19*t39*t40*t41*t42
     +        *t51*t73+4.0827787561099116d+0*t1*t10*t12*t122*t19*t40
     +        *t51*t73-2.0413893780549558d+0*t1*t10*t12*t13*t135*t51
     +        *t73+4.0827787561099116d+0*t1*t10*t115*t12*t13*t136*t19
     +        *t73+1.6669082371342397d+1*t10*t115*t17*t20*t22*t23*t56
     +        *t68+8.165557512219823d+0*t1*t115*t12*t13*t19*t56*t68
     +        +1.6669082371342397d+1*t10*t17*t20*t22*t39*t51*t52*t68
     +        -4.0827787561099116d+0*t1*t12*t13*t50*t51*t68
     +        +8.165557512219823d+0*t1*t12*t19*t39*t40*t51*t68
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t182*t51*t68
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t184*t19*t68
     +        -4.0827787561099116d+0*t1*t12*t13*t51*t56*t67
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t184*t51*t67
     +        -8.165557512219823d+0*t1*t12*t26*t39*t40*t67
     +        +1.6669082371342397d+1*t10*t17*t20*t22*t51*t52*t56*t63
     +        +8.165557512219823d+0*t1*t12*t19*t40*t51*t56*t63
     +        +4.0827787561099116d+0*t1*t10*t12*t184*t19*t40*t51*t63
     +        -8.165557512219823d+0*t1*t12*t26*t40*t50*t63
     +        -4.0827787561099116d+0*t1*t10*t12*t182*t26*t40*t63
     +        +2.449667253665947d+1*t1*t12*t19*t23*t26*t39*t63
     +        +5.00072471140272d+1*t10*t110*t17*t20*t22*t26*t39*t63
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t196*t51*t56
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t19*t197*t56
     +        -1.6669082371342397d+1*t10*t122*t17*t20*t22*t26*t52
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t197*t50*t51
     +        +4.763241882128231d+0*t1*t10*t13*t19*t197*t41*t42*t51
     +        +4.763241882128231d+0*t1*t10*t13*t184*t19*t41*t42*t51
     +        +4.0827787561099116d+0*t1*t10*t12*t19*t197*t39*t40*t51
     +        -8.334541185671199d+0*t10*t136*t17*t20*t22*t23*t51
     +        -4.0827787561099116d+0*t1*t12*t13*t136*t19*t51+t489+t488
     +        +t487+2.0413893780549558d+0*t1*t10*t12*t13*t26*(t484
     +        +t483+t482+t481+t480+t479+t478+t477+t476+t455+t454+t436
     +        +t432)+t475+t474+t449+t448+t447+t446+t445+t444+t428+t427
     +        +t426+t423+t420-4.763241882128231d+0*t1*t10*t13*t196*t26
     +        *t41*t42-4.763241882128231d+0*t1*t10*t13*t182*t26*t41
     +        *t42-4.0827787561099116d+0*t1*t10*t12*t196*t26*t39*t40
     +        -8.165557512219823d+0*t1*t12*t122*t19*t26*t40
     +        +4.0827787561099116d+0*t1*t12*t13*t135*t26)
     +        +3.109069086965489d-2*t198*t9+3.109069086965489d-2*t185
     +        *t9+1.8654414521792934d-1*t3*t39*t63*t74*t8
     +        +9.327207260896466d-2*t58*t63*t74+9.327207260896466d-2
     +        *t39*t58*t74+9.327207260896466d-2*t122*t3*t58*t74
     +        +9.327207260896466d-2*t185*t3*t58*t63
     +        +9.327207260896466d-2*t198*t3*t39*t58)+t119*(
     +        -3.109069086965489d-2*t185*t3*t69*t9
     +        -3.109069086965489d-2*t198*t3*t57*t9-t138*t74)
     +        +6.218138173930978d-2*t271*t3*t57*t69*t74*t9)*wght+Cmat3
     +        (iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = t473+Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = t496+Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = 1.0d+0*(t30*(3.109069086965489d-2
     +        *t3*t9*(-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*
     +        (1.0001449422805439d+2*t10*t109*t146*t20*t22*t23*t94
     +        +2.667053179414784d+2*t10*t15*t20*t22*t52*t55*t63*t65
     +        +3.1115620426505813d+2*t10*t15*t20*t23*t41*t53*t55*t65
     +        +6.223124085301163d+2*t10*t20*t24*t41*t52*t53*t63
     +        -3.3338164742684795d+1*t10*t15*t153*t20*t22*t23*t55
     +        -3.3338164742684795d+1*t10*t146*t15*t20*t22*t23*t55
     +        -6.667632948536959d+1*t10*t145*t20*t22*t24*t52+t504+t503
     +        +t502+t501+t500+t499+t498+t497+t461+t460+t455+t454
     +        +3.3338164742684795d+2*t10*t110*t140*t20*t22*t24)
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t67*t68*t73
     +        -1.6331115024439646d+1*t1*t10*t115*t12*t19*t40*t63*t68
     +        *t73-1.9052967528512924d+1*t1*t10*t115*t13*t19*t41*t42
     +        *t68*t73+8.165557512219823d+0*t1*t10*t12*t40*t51*t63*t67
     +        *t73+9.526483764256461d+0*t1*t10*t13*t41*t42*t51*t67*t73
     +        -1.9052967528512924d+1*t1*t10*t19*t40*t41*t42*t51*t63
     +        *t73+4.0827787561099116d+0*t1*t10*t12*t145*t19*t40*t51
     +        *t73-1.2248336268329735d+1*t1*t10*t12*t140*t19*t23*t51
     +        *t73-2.0413893780549558d+0*t1*t10*t12*t13*t156*t51*t73
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t158*t19*t269*t73
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t157*t19*t73
     +        -8.165557512219823d+0*t1*t12*t13*t51*t67*t68
     +        +3.3338164742684795d+1*t10*t17*t20*t22*t51*t52*t63*t68
     +        +1.6331115024439646d+1*t1*t12*t19*t40*t51*t63*t68
     +        +3.889452553313227d+1*t10*t17*t20*t23*t41*t51*t53*t68
     +        +1.9052967528512924d+1*t1*t13*t19*t41*t42*t51*t68
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t196*t51*t68
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t19*t197*t68
     +        -1.6331115024439646d+1*t1*t12*t26*t40*t63*t67
     +        -4.0827787561099116d+0*t1*t10*t12*t13*t197*t51*t67
     +        -1.9052967528512924d+1*t1*t13*t26*t41*t42*t67
     +        +7.778905106626454d+1*t10*t17*t20*t26*t41*t52*t53*t63
     +        +8.165557512219823d+0*t1*t10*t12*t19*t197*t40*t51*t63
     +        +3.8105935057025847d+1*t1*t19*t26*t40*t41*t42*t63
     +        -8.165557512219823d+0*t1*t10*t12*t196*t26*t40*t63
     +        -1.6669082371342397d+1*t10*t145*t17*t20*t22*t26*t52
     +        +9.526483764256461d+0*t1*t10*t13*t19*t197*t41*t42*t51
     +        -8.334541185671199d+0*t10*t157*t17*t20*t22*t23*t51
     +        -4.0827787561099116d+0*t1*t12*t13*t157*t19*t51
     +        +2.0413893780549558d+0*t1*t10*t12*t13*t26*(t504+t503
     +        +t502+t501+t500+t499+t498+t497+t455+t454)+t449+t448+t447
     +        +t446+t445+t444-9.526483764256461d+0*t1*t10*t13*t196*t26
     +        *t41*t42-8.165557512219823d+0*t1*t12*t145*t19*t26*t40
     +        +2.449667253665947d+1*t1*t12*t140*t19*t23*t26
     +        +5.00072471140272d+1*t10*t110*t140*t17*t20*t22*t26
     +        +4.0827787561099116d+0*t1*t12*t13*t156*t26
     +        +1.6669082371342397d+1*t10*t115*t158*t17*t20*t22*t23
     +        +8.165557512219823d+0*t1*t115*t12*t13*t158*t19)
     +        +6.218138173930978d-2*t198*t9+1.8654414521792934d-1*t140
     +        *t3*t74*t8+1.8654414521792934d-1*t58*t63*t74
     +        +9.327207260896466d-2*t145*t3*t58*t74
     +        +1.8654414521792934d-1*t198*t3*t58*t63)+t119*(
     +        -6.218138173930978d-2*t198*t3*t69*t9-t160*t74)
     +        +6.218138173930978d-2*t161*t271*t3*t74*t9)*wght+Cmat3
     +        (iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = t496+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = t520+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = 1.0d+0*(t119*(
     +        -3.109069086965489d-2*t185*t3*t71*t9-t173*t74+t530)+t529
     +        +t521)*wght+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t520+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = 1.0d+0*(t119*(
     +        -3.109069086965489d-2*t185*t3*t74*t9
     +        -3.109069086965489d-2*t211*t3*t57*t9-t186*t74)+t30*
     +        (3.109069086965489d-2*t3*(3.3338164742684795d+1*t10*t115
     +        *t17*t20*t22*t23*t56*t73+1.6331115024439646d+1*t1*t115
     +        *t12*t13*t19*t56*t73+3.3338164742684795d+1*t10*t17*t20
     +        *t22*t39*t51*t52*t73-8.165557512219823d+0*t1*t12*t13*t50
     +        *t51*t73+1.6331115024439646d+1*t1*t12*t19*t39*t40*t51
     +        *t73-4.0827787561099116d+0*t1*t10*t12*t13*t182*t51*t73
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t184*t19*t73
     +        -3.3338164742684795d+1*t17*t20*t22*t23*t51*t56
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t19*t210*t269*t56
     +        +1.3611235077912595d+2*t10*t110*t111*t115*t19*t203*t24
     +        *t56-2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -6.667632948536959d+1*t15*t20*t22*t23*t48*t55+t537
     +        -1.333526589707392d+2*t20*t22*t24*t39*t52)+t536+t535
     +        +t534+t533+t532-6.667632948536959d+1*t17*t20*t22*t26*t39
     +        *t52-6.805617538956298d+1*t10*t110*t111*t203*t24*t50*t51
     +        +1.3611235077912595d+2*t10*t111*t19*t203*t24*t264*t39
     +        *t51-1.6669082371342397d+1*t10*t17*t184*t20*t22*t23*t51
     +        -8.165557512219823d+0*t1*t12*t13*t184*t19*t51
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t210*t50
     +        -8.165557512219823d+0*t1*t10*t115*t12*t19*t210*t39*t40
     +        +8.165557512219823d+0*t1*t12*t13*t182*t26)*t9
     +        +9.327207260896466d-2*t211*t3*t39*t58+t531)
     +        +6.218138173930978d-2*t209*t271*t3*t57*t9)*wght+Cmat3
     +        (iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = 1.0d+0*(t119*(
     +        -3.109069086965489d-2*t172*t3*t74*t9-t186*t71+t530)+t529
     +        +t521)*wght+Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t520+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = t546+Cmat3(iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = 1.0d+0*(t119*(
     +        -3.109069086965489d-2*t198*t3*t71*t9-t192*t74+t549)+t548
     +        +t547)*wght+Cmat3(iq,D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = t546+Cmat3(iq,D3_RB_GAA_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = 1.0d+0*(t119*(
     +        -3.109069086965489d-2*t198*t3*t74*t9
     +        -3.109069086965489d-2*t211*t3*t69*t9-t199*t74)+t30*
     +        (3.109069086965489d-2*t3*(3.3338164742684795d+1*t10*t115
     +        *t17*t20*t22*t23*t68*t73+1.6331115024439646d+1*t1*t115
     +        *t12*t13*t19*t68*t73-8.165557512219823d+0*t1*t12*t13*t51
     +        *t67*t73+3.3338164742684795d+1*t10*t17*t20*t22*t51*t52
     +        *t63*t73+1.6331115024439646d+1*t1*t12*t19*t40*t51*t63
     +        *t73-4.0827787561099116d+0*t1*t10*t12*t13*t196*t51*t73
     +        +8.165557512219823d+0*t1*t10*t115*t12*t13*t19*t197*t73
     +        -3.3338164742684795d+1*t17*t20*t22*t23*t51*t68
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t19*t210*t269*t68
     +        +1.3611235077912595d+2*t10*t110*t111*t115*t19*t203*t24
     +        *t68-6.805617538956298d+1*t10*t110*t111*t203*t24*t51*t67
     +        +4.0827787561099116d+0*t1*t10*t115*t12*t13*t210*t67
     +        -2.0413893780549558d+0*t1*t10*t12*t13*t19*t51*(
     +        -6.667632948536959d+1*t15*t20*t22*t23*t55*t65
     +        -1.333526589707392d+2*t20*t22*t24*t52*t63+t537)
     +        -6.667632948536959d+1*t17*t20*t22*t26*t52*t63
     +        +1.3611235077912595d+2*t10*t111*t19*t203*t24*t264*t51
     +        *t63-8.165557512219823d+0*t1*t10*t115*t12*t19*t210*t40
     +        *t63+t536+t535+t534+t533+t532-1.6669082371342397d+1*t10
     +        *t17*t197*t20*t22*t23*t51-8.165557512219823d+0*t1*t12
     +        *t13*t19*t197*t51+8.165557512219823d+0*t1*t12*t13*t196
     +        *t26)*t9+9.327207260896466d-2*t211*t3*t58*t63+t531)
     +        +6.218138173930978d-2*t209*t271*t3*t69*t9)*wght+Cmat3
     +        (iq,D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = 1.0d+0*(t119*(
     +        -3.109069086965489d-2*t191*t3*t74*t9-t199*t71+t549)+t548
     +        +t547)*wght+Cmat3(iq,D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t546+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t553+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = 3.109069086965489d-2*t3*t30*
     +        (t563+t562+t561+t560+t559+t558+t557+t556+t555+t554
     +        -2.0839373032694847d+2*t10*t51*t55*t550*t551*t552)*t9
     +        *wght+t565+t564+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t553+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = 3.109069086965489d-2*t3*t30*
     +        (t574+t573+t572+t571+t570+t569+t568+t567+t562+t557
     +        -4.167874606538969d+2*t10*t51*t55*t550*t551*t552)*t9
     +        *wght+t575+t566+Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = t576+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = t553+Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = 6.218138173930978d-2*t271*t3
     +        *t74**3*t9*wght-9.327207260896466d-2*t119*t211*t3*t74*t9
     +        *wght+3.109069086965489d-2*t3*t30*(
     +        -1.2248336268329735d+1*t1*t10*t12*t13*t19*t269*t73**3
     +        -1.0001449422805439d+2*t17*t20*t22*t23*t51*t73
     +        +4.083370523373779d+2*t10*t110*t111*t115*t19*t203*t24
     +        *t73-8.335749213077939d+2*t10*t51*t55*t550*t551*t552
     +        -4.083370523373779d+2*t110*t111*t19*t203*t24*t51
     +        +5.00072471140272d+1*t10*t115*t17*t20*t210*t22*t23
     +        +2.449667253665947d+1*t1*t115*t12*t13*t19*t210)*t9*wght
     +        +Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = 3.109069086965489d-2*t3*t30*
     +        (t574+t573+t572+t571+t570+t569+t568+t567+t562+t557
     +        -4.1678746065389694d+2*t10*t51*t55*t550*t551*t552)*t9
     +        *wght+t575+t566+Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = t576+Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t553+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
            call nwxcm_c_pw91lda_d3(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a,
     +        s1a2,s1a3)
            t2 = s1f
            t3 = rhoa**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhoa
            t6 = exp(-6.432793688582966d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 3.240503646201146d+0*gammaaa*t1*t4*t8
            t10 = t9+1.0d+0
            t11 = t1**2
            t12 = gammaaa**2
            t13 = 1/t3**4.0d+0
            t14 = 1/t7**2
            t15 = t9+1.0500863881042921d+1*t11*t12*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 3.240503646201146d+0*gammaaa*t1*t10*t16*t4+1.0d+0
            t18 = log(t17)
            t19 = s1a(D1_RA)
            t20 = 1/t17
            t21 = rhoa**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -7.561175174469341d+0*gammaaa*t1*t21*t22*t8
            t24 = 1/rhoa**2
            t25 = 6.432793688582966d+1*s1f*t24-6.432793688582966d+1
     +        *t19*t5
            t26 = -3.240503646201146d+0*gammaaa*t1*t14*t25*t4*t6
            t27 = t26+t23
            t28 = 1/t15**2
            t29 = 1/t3**5.0d+0
            t30 = 1/t7**3
            t31 = -2.1001727762085842d+1*t11*t12*t13*t25*t30*t6
     +        -4.900403144486697d+1*t11*t12*t14*t21*t29+t26+t23
            t32 = -3.240503646201146d+0*gammaaa*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*gammaaa*t1*t16*t27*t4
     +        -7.561175174469341d+0*gammaaa*t1*t10*t16*t21*t22
            t33 = 3.240503646201146d+0*t1*t4*t8+2.1001727762085842d+1
     +        *gammaaa*t11*t13*t14
            t34 = 1.0500863881042921d+1*gammaaa*t11*t13*t16*t8
     +        -3.240503646201146d+0*gammaaa*t1*t10*t28*t33*t4
     +        +3.240503646201146d+0*t1*t10*t16*t4
            t35 = s1a2(D2_RA_RA)
            t36 = rhoa**3.333333333333335d-1
            t37 = 1/rhoa**8.333333333333333d-1
            t38 = t25**2
            t39 = exp(-1.2865587377165932d+2*s1f*t5)
            t40 = 6.481007292402292d+0*gammaaa*t1*t30*t38*t39*t4
            t41 = 2.6464113110642695d+1*gammaaa*t1*t13*t36*t8
            t42 = -1.2601958624115575d+0*gammaaa*t1*t22*t37*t8
            t43 = 1.5122350348938682d+1*gammaaa*t1*t14*t21*t22*t25*t6
            t44 = -3.240503646201146d+0*gammaaa*t1*t14*t38*t4*t6
            t45 = 1/rhoa**3
            t46 = -6.432793688582966d+1*t35*t5-1.2865587377165932d+2
     +        *s1f*t45+1.2865587377165932d+2*t19*t24
            t47 = -3.240503646201146d+0*gammaaa*t1*t14*t4*t46*t6
            t48 = t47+t44+t43+t42+t41+t40
            t49 = 1/t7**4
            t50 = 1/t3**6.0d+0
            t51 = -2.1001727762085842d+1*t11*t12*t13*t30*t46*t6
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t38*t6
     +        +1.9601612577946786d+2*t11*t12*t21*t25*t29*t30*t6
     +        +2.8585685009505735d+2*t11*t12*t14*t36*t50
     +        +6.300518328625753d+1*t11*t12*t13*t38*t39*t49+t47+t44
     +        +t43+t42+t41+t40-8.167338574144498d+0*t11*t12*t14*t29*t37
            t52 = 1/t15**3
            t53 = t31**2
            t54 = 6.481007292402292d+0*gammaaa*t1*t10*t4*t52*t53
     +        -3.240503646201146d+0*gammaaa*t1*t10*t28*t4*t51
     +        +3.240503646201146d+0*gammaaa*t1*t16*t4*t48
     +        -6.481007292402292d+0*gammaaa*t1*t27*t28*t31*t4
     +        -1.2601958624115575d+0*gammaaa*t1*t10*t16*t22*t37
     +        +2.6464113110642695d+1*gammaaa*t1*t10*t13*t16*t36
     +        +1.5122350348938682d+1*gammaaa*t1*t10*t21*t22*t28*t31
     +        -1.5122350348938682d+1*gammaaa*t1*t16*t21*t22*t27
            t55 = 1.5545345434827443d-2*rhoa*t54+3.1090690869654886d-2
     +        *t32
            t56 = 1/t17**2
            t57 = t32**2
            t58 = -7.561175174469341d+0*t1*t21*t22*t8
            t59 = -3.240503646201146d+0*t1*t14*t25*t4*t6
            t60 = t59+t58
            t61 = -4.2003455524171684d+1*gammaaa*t11*t13*t25*t30*t6
     +        +t59+t58-9.800806288973394d+1*gammaaa*t11*t14*t21*t29
            t62 = -1.0500863881042921d+1*gammaaa*t11*t13*t28*t31*t8
     +        -2.4502015722433484d+1*gammaaa*t11*t16*t21*t29*t8
     +        -3.240503646201146d+0*gammaaa*t1*t10*t28*t4*t61
     +        +3.240503646201146d+0*gammaaa*t1*t16*t4*t60
     +        +6.481007292402292d+0*gammaaa*t1*t10*t31*t33*t4*t52
     +        -3.240503646201146d+0*gammaaa*t1*t27*t28*t33*t4
     +        -3.240503646201146d+0*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*t1*t16*t27*t4+7.561175174469341d+0
     +        *gammaaa*t1*t10*t21*t22*t28*t33-7.561175174469341d+0*t1
     +        *t10*t16*t21*t22
            t63 = 1.5545345434827443d-2*rhoa*t62+1.5545345434827443d-2
     +        *t34
            t64 = t34**2
            t65 = t33**2
            t66 = t1**3
            t67 = -2.1001727762085842d+1*gammaaa*t11*t13*t28*t33*t8
     +        +2.1001727762085842d+1*t11*t13*t16*t8
     +        -6.805617538956301d+1*gammaaa*t10*t14*t28*t50*t66
     +        +6.481007292402292d+0*gammaaa*t1*t10*t4*t52*t65
     +        -6.481007292402292d+0*t1*t10*t28*t33*t4
            t68 = s1a3(D3_RA_RA_RA)
            t69 = rhoa**5.000000000000002d-1
            t70 = 1/rhoa**6.666666666666665d-1
            t71 = 1/rhoa**1.8333333333333332d+0
            t72 = t25**3
            t73 = exp(-1.9298381065748896d+2*s1f*t5)
            t74 = -1.9443021877206873d+1*gammaaa*t1*t4*t49*t72*t73
            t75 = -4.536705104681604d+1*gammaaa*t1*t21*t22*t30*t38*t39
            t76 = 1.2865587377165932d+2*s1f*t24-1.2865587377165932d+2
     +        *t19*t5
            t77 = 6.481007292402292d+0*gammaaa*t1*t30*t38*t39*t4*t76
            t78 = 6.481007292402292d+0*gammaaa*t1*t30*t39*t4*t72
            t79 = 1.9443021877206873d+1*gammaaa*t1*t25*t30*t39*t4*t46
            t80 = -1.2349919451633258d+2*gammaaa*t1*t29*t69*t8
            t81 = 1.3232056555321353d+1*gammaaa*t1*t13*t70*t8
            t82 = 1.0501632186762977d+0*gammaaa*t1*t22*t71*t8
            t83 = -7.939233933192807d+1*gammaaa*t1*t13*t14*t25*t36*t6
            t84 = 3.780587587234672d+0*gammaaa*t1*t14*t22*t25*t37*t6
            t85 = 2.268352552340802d+1*gammaaa*t1*t14*t21*t22*t38*t6
            t86 = 2.268352552340802d+1*gammaaa*t1*t14*t21*t22*t46*t6
            t87 = -3.240503646201146d+0*gammaaa*t1*t14*t4*t6*t72
            t88 = -9.721510938603437d+0*gammaaa*t1*t14*t25*t4*t46*t6
            t89 = -6.432793688582966d+1*t5*t68-3.859676213149779d+2
     +        *t19*t45+1.9298381065748896d+2*t24*t35
     +        +3.859676213149779d+2*s1f/rhoa**4
            t90 = -3.240503646201146d+0*gammaaa*t1*t14*t4*t6*t89
            t91 = 1/t3**7.0d+0
            t92 = 1/t15**4
            t93 = 1/t17**3
            t94 = 6.481007292402292d+0*t1*t30*t38*t39*t4
            t95 = 2.6464113110642695d+1*t1*t13*t36*t8
            t96 = -1.2601958624115575d+0*t1*t22*t37*t8
            t97 = 1.5122350348938682d+1*t1*t14*t21*t22*t25*t6
            t98 = -3.240503646201146d+0*t1*t14*t38*t4*t6
            t99 = -3.240503646201146d+0*t1*t14*t4*t46*t6
            fnc(iq) = 1.0d+0*(s1f+1.5545345434827443d-2*rhoa*t18)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(1.5545345434827443d-2*rhoa*t20
     +        *t32+t19+1.5545345434827443d-2*t18)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.5545345434827443d-2*rhoa*t20*t34*wght
     +        +Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-1.5545345434827443d-2*rhoa
     +        *t56*t57+t20*t55+t35)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(t20*t63
     +        -1.5545345434827443d-2*rhoa*t32*t34*t56)*wght+Cmat2(iq,
     +        D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.5545345434827443d-2*rhoa*t20*t67
     +        *wght-1.5545345434827443d-2*rhoa*t56*t64*wght+Cmat2(iq,
     +        D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(3.1090690869654886d-2*rhoa
     +        *t32**3*t93+t20*(1.5545345434827443d-2*rhoa*(
     +        -1.9443021877206873d+1*gammaaa*t1*t10*t31**3*t4*t92
     +        -3.240503646201146d+0*gammaaa*t1*t10*t28*t4*(
     +        -2.0009979506654013d+3*t11*t12*t14*t69*t91+t90
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t6*t89+t88+t87
     +        +t86+t85+t84+t83+t82+t81+t80+t79+t78+t77
     +        +6.300518328625753d+1*t11*t12*t13*t38*t39*t49*t76+t75
     +        +t74-2.520207331450301d+2*t11*t12*t13*t72*t73/t7**5
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t6*t72
     +        +6.300518328625753d+1*t11*t12*t13*t39*t49*t72
     +        +6.806115478453748d+0*t11*t12*t14*t29*t71
     +        +1.4292842504752873d+2*t11*t12*t14*t50*t70
     +        -1.715141100570344d+3*t11*t12*t25*t30*t36*t50*t6
     +        +2.940241886692018d+2*t11*t12*t21*t29*t30*t46*t6
     +        -6.300518328625753d+1*t11*t12*t13*t25*t30*t46*t6
     +        +2.940241886692018d+2*t11*t12*t21*t29*t30*t38*t6
     +        +4.900403144486699d+1*t11*t12*t25*t29*t30*t37*t6
     +        +1.8901554985877259d+2*t11*t12*t13*t25*t39*t46*t49
     +        -8.820725660076054d+2*t11*t12*t21*t29*t38*t39*t49)
     +        +3.240503646201146d+0*gammaaa*t1*t16*t4*(t90+t88+t87+t86
     +        +t85+t84+t83+t82+t81+t80+t79+t78+t77+t75+t74)
     +        +1.0501632186762977d+0*gammaaa*t1*t10*t16*t22*t71
     +        +1.3232056555321353d+1*gammaaa*t1*t10*t13*t16*t70
     +        -1.2349919451633258d+2*gammaaa*t1*t10*t16*t29*t69
     +        +1.9443021877206873d+1*gammaaa*t1*t27*t4*t52*t53
     +        -4.536705104681604d+1*gammaaa*t1*t10*t21*t22*t52*t53
     +        +1.9443021877206873d+1*gammaaa*t1*t10*t31*t4*t51*t52
     +        -9.721510938603437d+0*gammaaa*t1*t27*t28*t4*t51
     +        +2.268352552340802d+1*gammaaa*t1*t10*t21*t22*t28*t51
     +        -9.721510938603437d+0*gammaaa*t1*t28*t31*t4*t48
     +        -2.268352552340802d+1*gammaaa*t1*t16*t21*t22*t48
     +        +3.780587587234672d+0*gammaaa*t1*t10*t22*t28*t31*t37
     +        -3.780587587234672d+0*gammaaa*t1*t16*t22*t27*t37
     +        -7.939233933192807d+1*gammaaa*t1*t10*t13*t28*t31*t36
     +        +7.939233933192807d+1*gammaaa*t1*t13*t16*t27*t36
     +        +4.536705104681604d+1*gammaaa*t1*t21*t22*t27*t28*t31)
     +        +4.6636036304482326d-2*t54)+t68+t56*(
     +        -1.5545345434827443d-2*t57-t32*t55-3.1090690869654886d-2
     +        *rhoa*t32*t54))*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(t20*
     +        (1.5545345434827443d-2*rhoa*(-3.240503646201146d+0
     +        *gammaaa*t1*t10*t28*t4*(t99+t98+t97+t96+t95+t94
     +        -4.2003455524171684d+1*gammaaa*t11*t13*t30*t46*t6
     +        -4.2003455524171684d+1*gammaaa*t11*t13*t30*t38*t6
     +        +3.9203225155893573d+2*gammaaa*t11*t21*t25*t29*t30*t6
     +        +5.717137001901147d+2*gammaaa*t11*t14*t36*t50
     +        +1.2601036657251505d+2*gammaaa*t11*t13*t38*t39*t49
     +        -1.6334677148288995d+1*gammaaa*t11*t14*t29*t37)
     +        +3.240503646201146d+0*gammaaa*t1*t16*t4*(t99+t98+t97+t96
     +        +t95+t94)-1.9443021877206873d+1*gammaaa*t1*t10*t33*t4
     +        *t53*t92+2.1001727762085842d+1*gammaaa*t11*t13*t52*t53
     +        *t8-1.0500863881042921d+1*gammaaa*t11*t13*t28*t51*t8
     +        +8.575705502851718d+1*gammaaa*t11*t16*t36*t50*t8
     +        -4.083669287072249d+0*gammaaa*t11*t16*t29*t37*t8
     +        +4.900403144486697d+1*gammaaa*t11*t21*t28*t29*t31*t8
     +        +1.2962014584804582d+1*gammaaa*t1*t10*t31*t4*t52*t61
     +        -6.481007292402292d+0*gammaaa*t1*t27*t28*t4*t61
     +        +1.5122350348938682d+1*gammaaa*t1*t10*t21*t22*t28*t61
     +        -6.481007292402292d+0*gammaaa*t1*t28*t31*t4*t60
     +        -1.5122350348938682d+1*gammaaa*t1*t16*t21*t22*t60
     +        +6.481007292402292d+0*t1*t10*t4*t52*t53
     +        +6.481007292402292d+0*gammaaa*t1*t10*t33*t4*t51*t52
     +        +1.2962014584804582d+1*gammaaa*t1*t27*t31*t33*t4*t52
     +        -3.0244700697877364d+1*gammaaa*t1*t10*t21*t22*t31*t33
     +        *t52-3.240503646201146d+0*t1*t10*t28*t4*t51
     +        -3.240503646201146d+0*gammaaa*t1*t28*t33*t4*t48
     +        +3.240503646201146d+0*t1*t16*t4*t48-6.481007292402292d+0
     +        *t1*t27*t28*t31*t4+1.2601958624115575d+0*gammaaa*t1*t10
     +        *t22*t28*t33*t37-1.2601958624115575d+0*t1*t10*t16*t22
     +        *t37-2.6464113110642695d+1*gammaaa*t1*t10*t13*t28*t33
     +        *t36+2.6464113110642695d+1*t1*t10*t13*t16*t36
     +        +1.5122350348938682d+1*gammaaa*t1*t21*t22*t27*t28*t33
     +        +1.5122350348938682d+1*t1*t10*t21*t22*t28*t31
     +        -1.5122350348938682d+1*t1*t16*t21*t22*t27)
     +        +3.1090690869654886d-2*t62)+3.1090690869654886d-2*rhoa
     +        *t34*t57*t93+t56*(-3.1090690869654886d-2*rhoa*t32*t62
     +        -t34*t55))*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(3.1090690869654886d-2
     +        *rhoa*t32*t64*t93+t20*(1.5545345434827443d-2*rhoa*(
     +        -1.9443021877206873d+1*gammaaa*t1*t10*t31*t4*t65*t92
     +        +1.5879774257564702d+2*gammaaa*t10*t14*t21*t28*t66*t91
     +        -2.1001727762085842d+1*gammaaa*t11*t13*t28*t61*t8
     +        +4.2003455524171684d+1*gammaaa*t11*t13*t31*t33*t52*t8
     +        +4.900403144486697d+1*gammaaa*t11*t21*t28*t29*t33*t8
     +        -2.1001727762085842d+1*t11*t13*t28*t31*t8
     +        -4.900403144486697d+1*t11*t16*t21*t29*t8
     +        +1.36112350779126d+2*gammaaa*t10*t14*t31*t50*t52*t66
     +        -6.805617538956301d+1*gammaaa*t14*t27*t28*t50*t66
     +        +6.481007292402292d+0*gammaaa*t1*t27*t4*t52*t65
     +        -1.5122350348938682d+1*gammaaa*t1*t10*t21*t22*t52*t65
     +        +1.2962014584804582d+1*gammaaa*t1*t10*t33*t4*t52*t61
     +        -6.481007292402292d+0*t1*t10*t28*t4*t61
     +        -6.481007292402292d+0*gammaaa*t1*t28*t33*t4*t60
     +        +6.481007292402292d+0*t1*t16*t4*t60-3.240503646201146d+0
     +        *gammaaa*t1*t10*t28*t4*(-4.2003455524171684d+1*t11*t13
     +        *t25*t30*t6-9.800806288973394d+1*t11*t14*t21*t29)
     +        +1.2962014584804582d+1*t1*t10*t31*t33*t4*t52
     +        -6.481007292402292d+0*t1*t27*t28*t33*t4
     +        +1.5122350348938682d+1*t1*t10*t21*t22*t28*t33)
     +        +1.5545345434827443d-2*t67)+t56*(-1.5545345434827443d-2
     +        *rhoa*t32*t67-t34*t63-1.5545345434827443d-2*rhoa*t34*t62
     +        ))*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 3.1090690869654886d-2*rhoa*t34*
     +        *3*t93*wght+1.5545345434827443d-2*rhoa*t20*(
     +        -1.9443021877206873d+1*gammaaa*t1*t10*t33**3*t4*t92
     +        +6.300518328625753d+1*gammaaa*t11*t13*t52*t65*t8
     +        -6.300518328625753d+1*t11*t13*t28*t33*t8
     +        +4.0833705233737805d+2*gammaaa*t10*t14*t33*t50*t52*t66
     +        -2.0416852616868902d+2*t10*t14*t28*t50*t66
     +        +1.9443021877206873d+1*t1*t10*t4*t52*t65
     +        -6.616088534891508d+2*gammaaa*t1**4*t28*t30/t3**8.0d+0)
     +        *wght-4.6636036304482326d-2*rhoa*t34*t56*t67*wght+Cmat3
     +        (iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhob
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call dcopy(NCOL_AMAT2,0.0d0,0,s1a2,1)
            call dcopy(NCOL_AMAT3,0.0d0,0,s1a3,1)
            call nwxcm_c_pw91lda_d3(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a,
     +        s1a2,s1a3)
            t2 = s1f
            t3 = rhob**1.1666666666666667d+0
            t4 = 1/t3**2.0d+0
            t5 = 1/rhob
            t6 = exp(-6.432793688582966d+1*s1f*t5)
            t7 = t6-1.0d+0
            t8 = 1/t7
            t9 = 3.240503646201146d+0*gammabb*t1*t4*t8
            t10 = t9+1.0d+0
            t11 = t1**2
            t12 = gammabb**2
            t13 = 1/t3**4.0d+0
            t14 = 1/t7**2
            t15 = t9+1.0500863881042921d+1*t11*t12*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 3.240503646201146d+0*gammabb*t1*t10*t16*t4+1.0d+0
            t18 = log(t17)
            t19 = s1a(D1_RB)
            t20 = 1/t17
            t21 = rhob**1.6666666666666674d-1
            t22 = 1/t3**3.0d+0
            t23 = -7.561175174469341d+0*gammabb*t1*t21*t22*t8
            t24 = 1/rhob**2
            t25 = 6.432793688582966d+1*s1f*t24-6.432793688582966d+1
     +        *t19*t5
            t26 = -3.240503646201146d+0*gammabb*t1*t14*t25*t4*t6
            t27 = t26+t23
            t28 = 1/t15**2
            t29 = 1/t3**5.0d+0
            t30 = 1/t7**3
            t31 = -2.1001727762085842d+1*t11*t12*t13*t25*t30*t6
     +        -4.900403144486697d+1*t11*t12*t14*t21*t29+t26+t23
            t32 = -3.240503646201146d+0*gammabb*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*gammabb*t1*t16*t27*t4
     +        -7.561175174469341d+0*gammabb*t1*t10*t16*t21*t22
            t33 = 3.240503646201146d+0*t1*t4*t8+2.1001727762085842d+1
     +        *gammabb*t11*t13*t14
            t34 = 1.0500863881042921d+1*gammabb*t11*t13*t16*t8
     +        -3.240503646201146d+0*gammabb*t1*t10*t28*t33*t4
     +        +3.240503646201146d+0*t1*t10*t16*t4
            t35 = s1a2(D2_RB_RB)
            t36 = rhob**3.333333333333335d-1
            t37 = 1/rhob**8.333333333333333d-1
            t38 = t25**2
            t39 = exp(-1.2865587377165932d+2*s1f*t5)
            t40 = 6.481007292402292d+0*gammabb*t1*t30*t38*t39*t4
            t41 = 2.6464113110642695d+1*gammabb*t1*t13*t36*t8
            t42 = -1.2601958624115575d+0*gammabb*t1*t22*t37*t8
            t43 = 1.5122350348938682d+1*gammabb*t1*t14*t21*t22*t25*t6
            t44 = -3.240503646201146d+0*gammabb*t1*t14*t38*t4*t6
            t45 = 1/rhob**3
            t46 = -6.432793688582966d+1*t35*t5-1.2865587377165932d+2
     +        *s1f*t45+1.2865587377165932d+2*t19*t24
            t47 = -3.240503646201146d+0*gammabb*t1*t14*t4*t46*t6
            t48 = t47+t44+t43+t42+t41+t40
            t49 = 1/t7**4
            t50 = 1/t3**6.0d+0
            t51 = -2.1001727762085842d+1*t11*t12*t13*t30*t46*t6
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t38*t6
     +        +1.9601612577946786d+2*t11*t12*t21*t25*t29*t30*t6
     +        +2.8585685009505735d+2*t11*t12*t14*t36*t50
     +        +6.300518328625753d+1*t11*t12*t13*t38*t39*t49+t47+t44
     +        +t43+t42+t41+t40-8.167338574144498d+0*t11*t12*t14*t29*t37
            t52 = 1/t15**3
            t53 = t31**2
            t54 = 6.481007292402292d+0*gammabb*t1*t10*t4*t52*t53
     +        -3.240503646201146d+0*gammabb*t1*t10*t28*t4*t51
     +        +3.240503646201146d+0*gammabb*t1*t16*t4*t48
     +        -6.481007292402292d+0*gammabb*t1*t27*t28*t31*t4
     +        -1.2601958624115575d+0*gammabb*t1*t10*t16*t22*t37
     +        +2.6464113110642695d+1*gammabb*t1*t10*t13*t16*t36
     +        +1.5122350348938682d+1*gammabb*t1*t10*t21*t22*t28*t31
     +        -1.5122350348938682d+1*gammabb*t1*t16*t21*t22*t27
            t55 = 1.5545345434827443d-2*rhob*t54+3.1090690869654886d-2
     +        *t32
            t56 = 1/t17**2
            t57 = t32**2
            t58 = -7.561175174469341d+0*t1*t21*t22*t8
            t59 = -3.240503646201146d+0*t1*t14*t25*t4*t6
            t60 = t59+t58
            t61 = -4.2003455524171684d+1*gammabb*t11*t13*t25*t30*t6
     +        +t59+t58-9.800806288973394d+1*gammabb*t11*t14*t21*t29
            t62 = -1.0500863881042921d+1*gammabb*t11*t13*t28*t31*t8
     +        -2.4502015722433484d+1*gammabb*t11*t16*t21*t29*t8
     +        -3.240503646201146d+0*gammabb*t1*t10*t28*t4*t61
     +        +3.240503646201146d+0*gammabb*t1*t16*t4*t60
     +        +6.481007292402292d+0*gammabb*t1*t10*t31*t33*t4*t52
     +        -3.240503646201146d+0*gammabb*t1*t27*t28*t33*t4
     +        -3.240503646201146d+0*t1*t10*t28*t31*t4
     +        +3.240503646201146d+0*t1*t16*t27*t4+7.561175174469341d+0
     +        *gammabb*t1*t10*t21*t22*t28*t33-7.561175174469341d+0*t1
     +        *t10*t16*t21*t22
            t63 = 1.5545345434827443d-2*rhob*t62+1.5545345434827443d-2
     +        *t34
            t64 = t34**2
            t65 = t33**2
            t66 = t1**3
            t67 = -2.1001727762085842d+1*gammabb*t11*t13*t28*t33*t8
     +        +2.1001727762085842d+1*t11*t13*t16*t8
     +        -6.805617538956301d+1*gammabb*t10*t14*t28*t50*t66
     +        +6.481007292402292d+0*gammabb*t1*t10*t4*t52*t65
     +        -6.481007292402292d+0*t1*t10*t28*t33*t4
            t68 = s1a3(D3_RB_RB_RB)
            t69 = rhob**5.000000000000002d-1
            t70 = 1/rhob**6.666666666666665d-1
            t71 = 1/rhob**1.8333333333333332d+0
            t72 = t25**3
            t73 = exp(-1.9298381065748896d+2*s1f*t5)
            t74 = -1.9443021877206873d+1*gammabb*t1*t4*t49*t72*t73
            t75 = -4.536705104681604d+1*gammabb*t1*t21*t22*t30*t38*t39
            t76 = 1.2865587377165932d+2*s1f*t24-1.2865587377165932d+2
     +        *t19*t5
            t77 = 6.481007292402292d+0*gammabb*t1*t30*t38*t39*t4*t76
            t78 = 6.481007292402292d+0*gammabb*t1*t30*t39*t4*t72
            t79 = 1.9443021877206873d+1*gammabb*t1*t25*t30*t39*t4*t46
            t80 = -1.2349919451633258d+2*gammabb*t1*t29*t69*t8
            t81 = 1.3232056555321353d+1*gammabb*t1*t13*t70*t8
            t82 = 1.0501632186762977d+0*gammabb*t1*t22*t71*t8
            t83 = -7.939233933192807d+1*gammabb*t1*t13*t14*t25*t36*t6
            t84 = 3.780587587234672d+0*gammabb*t1*t14*t22*t25*t37*t6
            t85 = 2.268352552340802d+1*gammabb*t1*t14*t21*t22*t38*t6
            t86 = 2.268352552340802d+1*gammabb*t1*t14*t21*t22*t46*t6
            t87 = -3.240503646201146d+0*gammabb*t1*t14*t4*t6*t72
            t88 = -9.721510938603437d+0*gammabb*t1*t14*t25*t4*t46*t6
            t89 = -6.432793688582966d+1*t5*t68-3.859676213149779d+2
     +        *t19*t45+1.9298381065748896d+2*t24*t35
     +        +3.859676213149779d+2*s1f/rhob**4
            t90 = -3.240503646201146d+0*gammabb*t1*t14*t4*t6*t89
            t91 = 1/t3**7.0d+0
            t92 = 1/t15**4
            t93 = 1/t17**3
            t94 = 6.481007292402292d+0*t1*t30*t38*t39*t4
            t95 = 2.6464113110642695d+1*t1*t13*t36*t8
            t96 = -1.2601958624115575d+0*t1*t22*t37*t8
            t97 = 1.5122350348938682d+1*t1*t14*t21*t22*t25*t6
            t98 = -3.240503646201146d+0*t1*t14*t38*t4*t6
            t99 = -3.240503646201146d+0*t1*t14*t4*t46*t6
            fnc(iq) = 1.0d+0*(s1f+1.5545345434827443d-2*rhob*t18)*wght
     +        +fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*(1.5545345434827443d-2*rhob*t20
     +        *t32+t19+1.5545345434827443d-2*t18)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.5545345434827443d-2*rhob*t20*t34*wght
     +        +Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-1.5545345434827443d-2*rhob
     +        *t56*t57+t20*t55+t35)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(t20*t63
     +        -1.5545345434827443d-2*rhob*t32*t34*t56)*wght+Cmat2(iq,
     +        D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5545345434827443d-2*rhob*t20*t67
     +        *wght-1.5545345434827443d-2*rhob*t56*t64*wght+Cmat2(iq,
     +        D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(3.1090690869654886d-2*rhob
     +        *t32**3*t93+t20*(1.5545345434827443d-2*rhob*(
     +        -1.9443021877206873d+1*gammabb*t1*t10*t31**3*t4*t92
     +        -3.240503646201146d+0*gammabb*t1*t10*t28*t4*(
     +        -2.0009979506654013d+3*t11*t12*t14*t69*t91+t90
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t6*t89+t88+t87
     +        +t86+t85+t84+t83+t82+t81+t80+t79+t78+t77
     +        +6.300518328625753d+1*t11*t12*t13*t38*t39*t49*t76+t75
     +        +t74-2.520207331450301d+2*t11*t12*t13*t72*t73/t7**5
     +        -2.1001727762085842d+1*t11*t12*t13*t30*t6*t72
     +        +6.300518328625753d+1*t11*t12*t13*t39*t49*t72
     +        +6.806115478453748d+0*t11*t12*t14*t29*t71
     +        +1.4292842504752873d+2*t11*t12*t14*t50*t70
     +        -1.715141100570344d+3*t11*t12*t25*t30*t36*t50*t6
     +        +2.940241886692018d+2*t11*t12*t21*t29*t30*t46*t6
     +        -6.300518328625753d+1*t11*t12*t13*t25*t30*t46*t6
     +        +2.940241886692018d+2*t11*t12*t21*t29*t30*t38*t6
     +        +4.900403144486699d+1*t11*t12*t25*t29*t30*t37*t6
     +        +1.8901554985877259d+2*t11*t12*t13*t25*t39*t46*t49
     +        -8.820725660076054d+2*t11*t12*t21*t29*t38*t39*t49)
     +        +3.240503646201146d+0*gammabb*t1*t16*t4*(t90+t88+t87+t86
     +        +t85+t84+t83+t82+t81+t80+t79+t78+t77+t75+t74)
     +        +1.0501632186762977d+0*gammabb*t1*t10*t16*t22*t71
     +        +1.3232056555321353d+1*gammabb*t1*t10*t13*t16*t70
     +        -1.2349919451633258d+2*gammabb*t1*t10*t16*t29*t69
     +        +1.9443021877206873d+1*gammabb*t1*t27*t4*t52*t53
     +        -4.536705104681604d+1*gammabb*t1*t10*t21*t22*t52*t53
     +        +1.9443021877206873d+1*gammabb*t1*t10*t31*t4*t51*t52
     +        -9.721510938603437d+0*gammabb*t1*t27*t28*t4*t51
     +        +2.268352552340802d+1*gammabb*t1*t10*t21*t22*t28*t51
     +        -9.721510938603437d+0*gammabb*t1*t28*t31*t4*t48
     +        -2.268352552340802d+1*gammabb*t1*t16*t21*t22*t48
     +        +3.780587587234672d+0*gammabb*t1*t10*t22*t28*t31*t37
     +        -3.780587587234672d+0*gammabb*t1*t16*t22*t27*t37
     +        -7.939233933192807d+1*gammabb*t1*t10*t13*t28*t31*t36
     +        +7.939233933192807d+1*gammabb*t1*t13*t16*t27*t36
     +        +4.536705104681604d+1*gammabb*t1*t21*t22*t27*t28*t31)
     +        +4.6636036304482326d-2*t54)+t68+t56*(
     +        -1.5545345434827443d-2*t57-t32*t55-3.1090690869654886d-2
     +        *rhob*t32*t54))*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*(t20*
     +        (1.5545345434827443d-2*rhob*(-3.240503646201146d+0
     +        *gammabb*t1*t10*t28*t4*(t99+t98+t97+t96+t95+t94
     +        -4.2003455524171684d+1*gammabb*t11*t13*t30*t46*t6
     +        -4.2003455524171684d+1*gammabb*t11*t13*t30*t38*t6
     +        +3.9203225155893573d+2*gammabb*t11*t21*t25*t29*t30*t6
     +        +5.717137001901147d+2*gammabb*t11*t14*t36*t50
     +        +1.2601036657251505d+2*gammabb*t11*t13*t38*t39*t49
     +        -1.6334677148288995d+1*gammabb*t11*t14*t29*t37)
     +        +3.240503646201146d+0*gammabb*t1*t16*t4*(t99+t98+t97+t96
     +        +t95+t94)-1.9443021877206873d+1*gammabb*t1*t10*t33*t4
     +        *t53*t92+2.1001727762085842d+1*gammabb*t11*t13*t52*t53
     +        *t8-1.0500863881042921d+1*gammabb*t11*t13*t28*t51*t8
     +        +8.575705502851718d+1*gammabb*t11*t16*t36*t50*t8
     +        -4.083669287072249d+0*gammabb*t11*t16*t29*t37*t8
     +        +4.900403144486697d+1*gammabb*t11*t21*t28*t29*t31*t8
     +        +1.2962014584804582d+1*gammabb*t1*t10*t31*t4*t52*t61
     +        -6.481007292402292d+0*gammabb*t1*t27*t28*t4*t61
     +        +1.5122350348938682d+1*gammabb*t1*t10*t21*t22*t28*t61
     +        -6.481007292402292d+0*gammabb*t1*t28*t31*t4*t60
     +        -1.5122350348938682d+1*gammabb*t1*t16*t21*t22*t60
     +        +6.481007292402292d+0*t1*t10*t4*t52*t53
     +        +6.481007292402292d+0*gammabb*t1*t10*t33*t4*t51*t52
     +        +1.2962014584804582d+1*gammabb*t1*t27*t31*t33*t4*t52
     +        -3.0244700697877364d+1*gammabb*t1*t10*t21*t22*t31*t33
     +        *t52-3.240503646201146d+0*t1*t10*t28*t4*t51
     +        -3.240503646201146d+0*gammabb*t1*t28*t33*t4*t48
     +        +3.240503646201146d+0*t1*t16*t4*t48-6.481007292402292d+0
     +        *t1*t27*t28*t31*t4+1.2601958624115575d+0*gammabb*t1*t10
     +        *t22*t28*t33*t37-1.2601958624115575d+0*t1*t10*t16*t22
     +        *t37-2.6464113110642695d+1*gammabb*t1*t10*t13*t28*t33
     +        *t36+2.6464113110642695d+1*t1*t10*t13*t16*t36
     +        +1.5122350348938682d+1*gammabb*t1*t21*t22*t27*t28*t33
     +        +1.5122350348938682d+1*t1*t10*t21*t22*t28*t31
     +        -1.5122350348938682d+1*t1*t16*t21*t22*t27)
     +        +3.1090690869654886d-2*t62)+3.1090690869654886d-2*rhob
     +        *t34*t57*t93+t56*(-3.1090690869654886d-2*rhob*t32*t62
     +        -t34*t55))*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*(3.1090690869654886d-2
     +        *rhob*t32*t64*t93+t20*(1.5545345434827443d-2*rhob*(
     +        -1.9443021877206873d+1*gammabb*t1*t10*t31*t4*t65*t92
     +        +1.5879774257564702d+2*gammabb*t10*t14*t21*t28*t66*t91
     +        -2.1001727762085842d+1*gammabb*t11*t13*t28*t61*t8
     +        +4.2003455524171684d+1*gammabb*t11*t13*t31*t33*t52*t8
     +        +4.900403144486697d+1*gammabb*t11*t21*t28*t29*t33*t8
     +        -2.1001727762085842d+1*t11*t13*t28*t31*t8
     +        -4.900403144486697d+1*t11*t16*t21*t29*t8
     +        +1.36112350779126d+2*gammabb*t10*t14*t31*t50*t52*t66
     +        -6.805617538956301d+1*gammabb*t14*t27*t28*t50*t66
     +        +6.481007292402292d+0*gammabb*t1*t27*t4*t52*t65
     +        -1.5122350348938682d+1*gammabb*t1*t10*t21*t22*t52*t65
     +        +1.2962014584804582d+1*gammabb*t1*t10*t33*t4*t52*t61
     +        -6.481007292402292d+0*t1*t10*t28*t4*t61
     +        -6.481007292402292d+0*gammabb*t1*t28*t33*t4*t60
     +        +6.481007292402292d+0*t1*t16*t4*t60-3.240503646201146d+0
     +        *gammabb*t1*t10*t28*t4*(-4.2003455524171684d+1*t11*t13
     +        *t25*t30*t6-9.800806288973394d+1*t11*t14*t21*t29)
     +        +1.2962014584804582d+1*t1*t10*t31*t33*t4*t52
     +        -6.481007292402292d+0*t1*t27*t28*t33*t4
     +        +1.5122350348938682d+1*t1*t10*t21*t22*t28*t33)
     +        +1.5545345434827443d-2*t67)+t56*(-1.5545345434827443d-2
     +        *rhob*t32*t67-t34*t63-1.5545345434827443d-2*rhob*t34*t62
     +        ))*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 3.1090690869654886d-2*rhob*t34*
     +        *3*t93*wght+1.5545345434827443d-2*rhob*t20*(
     +        -1.9443021877206873d+1*gammabb*t1*t10*t33**3*t4*t92
     +        +6.300518328625753d+1*gammabb*t11*t13*t52*t65*t8
     +        -6.300518328625753d+1*t11*t13*t28*t33*t8
     +        +4.0833705233737805d+2*gammabb*t10*t14*t33*t50*t52*t66
     +        -2.0416852616868902d+2*t10*t14*t28*t50*t66
     +        +1.9443021877206873d+1*t1*t10*t4*t52*t65
     +        -6.616088534891508d+2*gammabb*t1**4*t28*t30/t3**8.0d+0)
     +        *wght-4.6636036304482326d-2*rhob*t34*t56*t67*wght+Cmat3
     +        (iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_pbe_d3 
C> @}
