C> \ingroup nwxc
C> @{
C>
C> \file nwxcm_c_mpbe.F
C> The nwxcm_c_mpbe functional
C>
C> @}
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the nwxcm_c_mpbe functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\
C>   {\it t_2} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_\alpha , \rho_\beta\right)\\\\
C>   {\it t_3} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\
C>   {\it t_5} &=& \rho_\alpha-\rho_\beta\\\\
C>   {\it t_6} &=& 0.5\,\left({\it t_5}\,{\it t_4}
C>    +1.0\right)^{{{2}\over{3}}}+0.5\,\left(1.0-{\it t_5}\,{
C>    \it t_4}\right)^{{{2}\over{3}}}\\\\
C>   {\it t_7} &=& {\it t_6}^{3.0}\\\\
C>   {\it t_8} &=& {{1}\over{{\it t_3}^{2.333333333333333}}}\\\\
C>   {\it t_9} &=& {{1}\over{{\it t_6}^2}}\\\\
C>   {\it t_{10}} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\
C>   {\it t_{11}} &=& {{1}\over{e^{{{32.16396844291483\,{\it t_2}
C>    \,{\it t_4}}\over{{\it t_7}}}}}}-1.0\\\\
C>   {\it t_{12}} &=& {{2.041389378054956\,{\it t_1}\,{\it t_8}\,{
C>    \it t_9}\,{\it t_{10}}}\over{{\it t_{11}}}}\\\\
C>   {\it t_{13}} &=& {\it t_1}^2\\\\
C>   {\it t_{14}} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_s , 0.0\right)\\\\
C>   {\it t_{15}} &=& {{1}\over{\rho_s}}\\\\
C>   {\it t_{16}} &=& {{1}\over{\rho_s^{2.333333333333333}}}\\\\
C>   {\it t_{17}} &=& {{1}\over{e^{64.32793688582967\,{\it t_{14}}
C>    \,{\it t_{15}}}}}-1.0\\\\
C>   {\it t_{18}} &=& {{3.240503646201146\,{\it t_1}\,{\it t_{16}}
C>    \,\sigma_{ss}}\over{{\it t_{17}}}}\\\\
C>   f &=& 0.03109069086965489\,{\it t_7}\,\log 
C>    \left({{2.041389378054956\,{\it t_1}\,{\it t_8}\,{\it t_9}
C>    \,{\it t_{10}}\,\left({\it t_{12}}+1.0\right)}
C>    \over{{{4.1672705928356\,{\it t_{13}}\,{\it t_{10}}^2}
C>    \over{{\it t_3}^{4.666666666666667}\,{\it t_6}^4\,{
C>    \it t_{11}}^2}}+{\it t_{12}}+1.0}}+1.0\right)+1.0\,{\it t_2}
C>    \,{\it t_4}\\\\
C>   g &=& 0\\\\
C>   G &=& 0.01554534543482744\,\log \left({{3.240503646201146\,{
C>    \it t_1}\,{\it t_{16}}\,\sigma_{ss}\,\left({\it t_{18}}
C>    +1.0\right)}\over{{{10.50086388104292\,{\it t_{13}}\,
C>    \sigma_{ss}^2}\over{\rho_s^{4.666666666666667}\,{\it t_{17}}^2}}
C>    +{\it t_{18}}+1.0}}+1.0\right)+1.0\,{\it t_{14}}\,{
C>    \it t_{15}}\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc [4,5,6].
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
C> [4] autoxc, revision 27694 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_mpbe(param,tol_rho,ipol,nq,wght, rho,rgamma,
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
            t3 = 1/rhoa
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = exp(-1.6081984221457415d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 1.6202518231005725d+0*gammaaa*t1*t4*t7
            t9 = t1**2
            t10 = gammaaa**2
            t11 = 1/rhoa**4.666666666666667d+0
            t12 = 1/t6**2
            t13 = 2.6252159702607286d+0*t10*t11*t12*t9
            t14 = 1/rhoa**2
            t15 = s1a(D1_RA)
            t16 = 1.6202518231005728d+0*gammaaa*t1*t4*t7
            t17 = t16+1.0d+0
            t18 = t16+t13+1.0d+0
            t19 = 1/t18
            t20 = 1/(1.6202518231005728d+0*gammaaa*t1*t17*t19*t4+1.0d+0)
            t21 = 1/rhoa**3.3333333333333337d+0
            t22 = -1.8902937936173353d+0*gammaaa*t1*t21*t7
            t23 = 8.040992110728707d+0*t14*s1f-1.6081984221457415d+1
     +        *t15*t3
            t24 = -1.6202518231005728d+0*gammaaa*t1*t12*t23*t4*t5
            t25 = 1/t18**2
            fnc(iq) = (3.109069086965489d-2*log(1.6202518231005725d+0
     +        *gammaaa*t1*t4*(t8+1.0d+0)/(t8+t13+1.0d+0)+1.0d+0)
     +        +5.0d-1*s1f*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (3.109069086965489d-2*t20*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t25*t4*(
     +        -5.250431940521457d+0*t10*t11*t23*t5*t9/t6**3
     +        -6.12550393060837d+0*t10*t12*t9/rhoa*
     +        *5.666666666666667d+0+t24+t22)+1.6202518231005728d+0
     +        *gammaaa*t1*t19*(t24+t22)*t4-1.8902937936173353d+0
     +        *gammaaa*t1*t17*t19*t21)+5.0d-1*t15*t3-2.5d-1*t14*s1f)
     +        *wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 3.109069086965489d-2*t20*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t25*t4*
     +        (1.3126079851303643d+0*gammaaa*t11*t12*t9
     +        +4.050629557751432d-1*t1*t4*t7)+6.563039925651821d-1
     +        *gammaaa*t11*t19*t7*t9+4.050629557751432d-1*t1*t17*t19
     +        *t4)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t20*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t25*t4*
     +        (2.6252159702607286d+0*gammaaa*t11*t12*t9
     +        +8.101259115502863d-1*t1*t4*t7)+1.3126079851303643d+0
     +        *gammaaa*t11*t19*t7*t9+8.101259115502863d-1*t1*t17*t19
     +        *t4)*wght+Cmat(iq,D1_GAB)
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
            t4 = 1/t3
            t5 = rhoa-rhob
            t6 = 1.0d+0-t4*t5
            t7 = t4*t5+1.0d+0
            t8 = 5.0d-1*t7**6.666666666666666d-1+5.0d-1*t6*
     +        *6.666666666666666d-1
            t9 = t8**3.0d+0
            t10 = gammabb+2.0d+0*gammaab+gammaaa
            t11 = 1/t3**2.3333333333333334d+0
            t12 = 1/t8**2
            t13 = 1/t9
            t14 = exp(-3.216396844291483d+1*t13*s1f*t4)
            t15 = t14-1.0d+0
            t16 = 1/t15
            t17 = 2.0413893780549558d+0*t1*t10*t11*t12*t16
            t18 = t17+1.0d+0
            t19 = t1**2
            t20 = t10**2
            t21 = 1/t3**4.666666666666667d+0
            t22 = 1/t8**4
            t23 = 1/t15**2
            t24 = 4.1672705928355996d+0*t19*t20*t21*t22*t23+t17+1.0d+0
            t25 = 1/t24
            t26 = 2.0413893780549558d+0*t1*t10*t11*t12*t18*t25+1.0d+0
            t27 = log(t26)
            t28 = 1/t3**2
            t29 = -1.0d+0*s1f*t28
            t30 = s1a(D1_RA)
            t31 = 1/t26
            t32 = t28*t5
            t33 = -t4
            t34 = 1/t6**3.333333333333333d-1
            t35 = -t28*t5
            t36 = 1/t7**3.333333333333333d-1
            t37 = 3.333333333333333d-1*t36*(t4+t35)
     +        +3.333333333333333d-1*(t33+t32)*t34
            t38 = 1/t8**3
            t39 = 1/t3**3.3333333333333337d+0
            t40 = -4.763241882128231d+0*t1*t10*t12*t18*t25*t39
            t41 = -4.0827787561099116d+0*t1*t10*t11*t16*t37*t38
            t42 = -4.763241882128231d+0*t1*t10*t12*t16*t39
            t43 = 1/t8**4.0d+0
            t44 = 3.216396844291483d+1*t13*s1f*t28
            t45 = t44+9.649190532874448d+1*s1f*t37*t4*t43
     +        -3.216396844291483d+1*t13*t30*t4
            t46 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t45
            t47 = 1/t24**2
            t48 = 1/t8**5
            t49 = -1.9447262766566134d+1*t19*t20*t22*t23/t3*
     +        *5.666666666666667d+0
            t50 = 1/t15**3
            t51 = t8**2.0d+0
            t52 = s1a(D1_RB)
            t53 = 3.333333333333333d-1*t34*(t4+t32)
     +        +3.333333333333333d-1*(t35+t33)*t36
            t54 = -4.0827787561099116d+0*t1*t10*t11*t16*t38*t53
            t55 = 9.649190532874448d+1*s1f*t4*t43*t53
     +        -3.216396844291483d+1*t13*t4*t52+t44
            t56 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t55
            t57 = 3.109069086965489d-2*t31*(-2.0413893780549558d+0*t1
     +        *t10*t11*t12*t18*(8.334541185671199d+0*t10*t19*t21*t22
     +        *t23+2.0413893780549558d+0*t1*t11*t12*t16)*t47
     +        +4.1672705928355996d+0*t10*t16*t19*t21*t22*t25
     +        +2.0413893780549558d+0*t1*t11*t12*t18*t25)*t9*wght
            fnc(iq) = (3.109069086965489d-2*t27*t9+1.0d+0*s1f*t4)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = (3.109069086965489d-2*t31*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t47*(
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t45*t50+t49
     +        -1.6669082371342397d+1*t19*t20*t21*t23*t37*t48+t46+t42
     +        +t41)+2.0413893780549558d+0*t1*t10*t11*t12*t25*(t46+t42
     +        +t41)+t40-4.0827787561099116d+0*t1*t10*t11*t18*t25*t37
     +        *t38)*t9+9.327207260896466d-2*t27*t37*t51+1.0d+0*t30*t4
     +        +t29)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (3.109069086965489d-2*t31*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t47*(t56
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t50*t55+t54
     +        -1.6669082371342397d+1*t19*t20*t21*t23*t48*t53+t49+t42)
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t56+t54+t42)
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t38*t53+t40)
     +        *t9+9.327207260896466d-2*t27*t51*t53+1.0d+0*t4*t52+t29)
     +        *wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t57+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t31*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*
     +        (1.6669082371342397d+1*t10*t19*t21*t22*t23
     +        +4.0827787561099116d+0*t1*t11*t12*t16)*t47
     +        +8.334541185671199d+0*t10*t16*t19*t21*t22*t25
     +        +4.0827787561099116d+0*t1*t11*t12*t18*t25)*t9*wght+Cmat
     +        (iq,D1_GAB)
            Cmat(iq,D1_GBB) = t57+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhoa
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call nwxcm_c_pw91lda(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a)
            t2 = s1f
            t3 = 1/rhoa
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = exp(-6.432793688582967d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 3.240503646201146d+0*gammaaa*t1*t4*t7
            t9 = t8+1.0d+0
            t10 = t1**2
            t11 = gammaaa**2
            t12 = 1/rhoa**4.666666666666667d+0
            t13 = 1/t6**2
            t14 = t8+1.0500863881042923d+1*t10*t11*t12*t13+1.0d+0
            t15 = 1/t14
            t16 = 3.240503646201146d+0*gammaaa*t1*t15*t4*t9+1.0d+0
            t17 = 1/rhoa**2
            t18 = s1a(D1_RA)
            t19 = 1/t16
            t20 = 1/rhoa**3.3333333333333337d+0
            t21 = -7.561175174469341d+0*gammaaa*t1*t20*t7
            t22 = 6.432793688582967d+1*t17*s1f-6.432793688582967d+1
     +        *t18*t3
            t23 = -3.240503646201146d+0*gammaaa*t1*t13*t22*t4*t5
            t24 = 1/t14**2
            fnc(iq) = (1.0d+0*s1f*t3+1.5545345434827443d-2*log(t16))
     +        *wght+fnc(iq)
            Amat(iq,D1_RA) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammaaa*t1*t24*t4*(
     +        -2.1001727762085845d+1*t10*t11*t12*t22*t5/t6**3+t23+t21
     +        -4.9004031444866974d+1*t10*t11*t13/rhoa*
     +        *5.666666666666667d+0)*t9-7.561175174469341d+0*gammaaa
     +        *t1*t15*t20*t9+3.240503646201146d+0*gammaaa*t1*t15*(t23
     +        +t21)*t4)+1.0d+0*t18*t3-1.0d+0*t17*s1f)*wght+Amat(iq,
     +        D1_RA)
            Cmat(iq,D1_GAA) = 1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammaaa*t1*t24*t4*
     +        (3.240503646201146d+0*t1*t4*t7+2.1001727762085845d+1
     +        *gammaaa*t10*t12*t13)*t9+3.240503646201146d+0*t1*t15*t4
     +        *t9+1.0500863881042925d+1*gammaaa*t10*t12*t15*t7)*wght
     +        +Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            sr(R_A) = rhob
            sr(R_B) = 0.0d+0
            s1f = 0.0d0
            call dcopy(NCOL_AMAT,0.0d0,0,s1a,1)
            call nwxcm_c_pw91lda(t1,tol_rho,2,1,1.0d0,sr,s1f,s1a)
            t2 = s1f
            t3 = 1/rhob
            t4 = 1/rhob**2.3333333333333334d+0
            t5 = exp(-6.432793688582967d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 3.240503646201146d+0*gammabb*t1*t4*t7
            t9 = t8+1.0d+0
            t10 = t1**2
            t11 = gammabb**2
            t12 = 1/rhob**4.666666666666667d+0
            t13 = 1/t6**2
            t14 = t8+1.0500863881042923d+1*t10*t11*t12*t13+1.0d+0
            t15 = 1/t14
            t16 = 3.240503646201146d+0*gammabb*t1*t15*t4*t9+1.0d+0
            t17 = 1/rhob**2
            t18 = s1a(D1_RB)
            t19 = 1/t16
            t20 = 1/rhob**3.3333333333333337d+0
            t21 = -7.561175174469341d+0*gammabb*t1*t20*t7
            t22 = 6.432793688582967d+1*t17*s1f-6.432793688582967d+1
     +        *t18*t3
            t23 = -3.240503646201146d+0*gammabb*t1*t13*t22*t4*t5
            t24 = 1/t14**2
            fnc(iq) = (1.0d+0*s1f*t3+1.5545345434827443d-2*log(t16))
     +        *wght+fnc(iq)
            Amat(iq,D1_RB) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammabb*t1*t24*t4*(
     +        -2.1001727762085845d+1*t10*t11*t12*t22*t5/t6**3+t23+t21
     +        -4.9004031444866974d+1*t10*t11*t13/rhob*
     +        *5.666666666666667d+0)*t9-7.561175174469341d+0*gammabb
     +        *t1*t15*t20*t9+3.240503646201146d+0*gammabb*t1*t15*(t23
     +        +t21)*t4)+1.0d+0*t18*t3-1.0d+0*t17*s1f)*wght+Amat(iq,
     +        D1_RB)
            Cmat(iq,D1_GBB) = 1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammabb*t1*t24*t4*
     +        (3.240503646201146d+0*t1*t4*t7+2.1001727762085845d+1
     +        *gammabb*t10*t12*t13)*t9+3.240503646201146d+0*t1*t15*t4
     +        *t9+1.0500863881042925d+1*gammabb*t10*t12*t15*t7)*wght
     +        +Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_mpbe 
C>
C> \brief Evaluate the nwxcm_c_mpbe functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\
C>   {\it t_2} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_\alpha , \rho_\beta\right)\\\\
C>   {\it t_3} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\
C>   {\it t_5} &=& \rho_\alpha-\rho_\beta\\\\
C>   {\it t_6} &=& 0.5\,\left({\it t_5}\,{\it t_4}
C>    +1.0\right)^{{{2}\over{3}}}+0.5\,\left(1.0-{\it t_5}\,{
C>    \it t_4}\right)^{{{2}\over{3}}}\\\\
C>   {\it t_7} &=& {\it t_6}^{3.0}\\\\
C>   {\it t_8} &=& {{1}\over{{\it t_3}^{2.333333333333333}}}\\\\
C>   {\it t_9} &=& {{1}\over{{\it t_6}^2}}\\\\
C>   {\it t_{10}} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\
C>   {\it t_{11}} &=& {{1}\over{e^{{{32.16396844291483\,{\it t_2}
C>    \,{\it t_4}}\over{{\it t_7}}}}}}-1.0\\\\
C>   {\it t_{12}} &=& {{2.041389378054956\,{\it t_1}\,{\it t_8}\,{
C>    \it t_9}\,{\it t_{10}}}\over{{\it t_{11}}}}\\\\
C>   {\it t_{13}} &=& {\it t_1}^2\\\\
C>   {\it t_{14}} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_s , 0.0\right)\\\\
C>   {\it t_{15}} &=& {{1}\over{\rho_s}}\\\\
C>   {\it t_{16}} &=& {{1}\over{\rho_s^{2.333333333333333}}}\\\\
C>   {\it t_{17}} &=& {{1}\over{e^{64.32793688582967\,{\it t_{14}}
C>    \,{\it t_{15}}}}}-1.0\\\\
C>   {\it t_{18}} &=& {{3.240503646201146\,{\it t_1}\,{\it t_{16}}
C>    \,\sigma_{ss}}\over{{\it t_{17}}}}\\\\
C>   f &=& 0.03109069086965489\,{\it t_7}\,\log 
C>    \left({{2.041389378054956\,{\it t_1}\,{\it t_8}\,{\it t_9}
C>    \,{\it t_{10}}\,\left({\it t_{12}}+1.0\right)}
C>    \over{{{4.1672705928356\,{\it t_{13}}\,{\it t_{10}}^2}
C>    \over{{\it t_3}^{4.666666666666667}\,{\it t_6}^4\,{
C>    \it t_{11}}^2}}+{\it t_{12}}+1.0}}+1.0\right)+1.0\,{\it t_2}
C>    \,{\it t_4}\\\\
C>   g &=& 0\\\\
C>   G &=& 0.01554534543482744\,\log \left({{3.240503646201146\,{
C>    \it t_1}\,{\it t_{16}}\,\sigma_{ss}\,\left({\it t_{18}}
C>    +1.0\right)}\over{{{10.50086388104292\,{\it t_{13}}\,
C>    \sigma_{ss}^2}\over{\rho_s^{4.666666666666667}\,{\it t_{17}}^2}}
C>    +{\it t_{18}}+1.0}}+1.0\right)+1.0\,{\it t_{14}}\,{
C>    \it t_{15}}\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc [4,5,6].
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
C> [4] autoxc, revision 27694 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_mpbe_d2(param,tol_rho,ipol,nq,wght, rho,
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
            t3 = 1/rhoa
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = exp(-1.6081984221457415d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 1.6202518231005725d+0*gammaaa*t1*t4*t7
            t9 = t1**2
            t10 = gammaaa**2
            t11 = 1/rhoa**4.666666666666667d+0
            t12 = 1/t6**2
            t13 = 2.6252159702607286d+0*t10*t11*t12*t9
            t14 = 1/rhoa**2
            t15 = s1a(D1_RA)
            t16 = 1.6202518231005728d+0*gammaaa*t1*t4*t7
            t17 = t16+1.0d+0
            t18 = t16+t13+1.0d+0
            t19 = 1/t18
            t20 = 1.6202518231005728d+0*gammaaa*t1*t17*t19*t4+1.0d+0
            t21 = 1/t20
            t22 = 1/rhoa**3.3333333333333337d+0
            t23 = -1.8902937936173353d+0*gammaaa*t1*t22*t7
            t24 = 8.040992110728707d+0*t14*s1f-1.6081984221457415d+1
     +        *t15*t3
            t25 = -1.6202518231005728d+0*gammaaa*t1*t12*t24*t4*t5
            t26 = t25+t23
            t27 = 1/t18**2
            t28 = 1/rhoa**5.666666666666667d+0
            t29 = 1/t6**3
            t30 = -5.250431940521457d+0*t10*t11*t24*t29*t5*t9
     +        -6.12550393060837d+0*t10*t12*t28*t9+t25+t23
            t31 = -1.6202518231005728d+0*gammaaa*t1*t17*t27*t30*t4
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t26*t4
     +        -1.8902937936173353d+0*gammaaa*t1*t17*t19*t22
            t32 = 1.3126079851303643d+0*gammaaa*t11*t12*t9
     +        +4.050629557751432d-1*t1*t4*t7
            t33 = 6.563039925651821d-1*gammaaa*t11*t19*t7*t9
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t32*t4
     +        +4.050629557751432d-1*t1*t17*t19*t4
            t34 = 2.6252159702607286d+0*gammaaa*t11*t12*t9
     +        +8.101259115502863d-1*t1*t4*t7
            t35 = 1.3126079851303643d+0*gammaaa*t11*t19*t7*t9
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t34*t4
     +        +8.101259115502863d-1*t1*t17*t19*t4
            t36 = 1/rhoa**3
            t37 = 2.5d-1*s1f*t36
            t38 = -5.0d-1*t14*t15
            t39 = s1a2(D2_RA_RA)
            t40 = 1/rhoa**4.333333333333334d+0
            t41 = -3.7805875872346706d+0*gammaaa*t1*t19*t22*t26
            t42 = 3.7805875872346706d+0*gammaaa*t1*t17*t22*t27*t30
            t43 = -3.2405036462011455d+0*gammaaa*t1*t26*t27*t30*t4
            t44 = t24**2
            t45 = exp(-3.216396844291483d+1*s1f*t3)
            t46 = 3.2405036462011455d+0*gammaaa*t1*t29*t4*t44*t45
            t47 = 3.3305176363733985d+0*gammaaa*t1*t40*t7
            t48 = 3.7805875872346706d+0*gammaaa*t1*t12*t22*t24*t5
            t49 = -1.6202518231005728d+0*gammaaa*t1*t12*t4*t44*t5
            t50 = 1.6081984221457415d+1*t14*t15
            t51 = t50-1.6081984221457415d+1*t3*t39
     +        -1.0721322814304943d+1*s1f*t36
            t52 = -1.6202518231005728d+0*gammaaa*t1*t12*t4*t5*t51
            t53 = 1.5751295821564373d+1*t10*t11*t44*t45*t9/t6**4
            t54 = 1/rhoa**6.666666666666667d+0
            t55 = 2.450201572243348d+1*t10*t24*t28*t29*t5*t9
            t56 = -5.250431940521457d+0*t10*t11*t29*t44*t5*t9
            t57 = 1/t18**3
            t58 = 3.2405036462011455d+0*gammaaa*t1*t17*t30**2*t4*t57
            t59 = 1/t20**2
            t60 = -3.109069086965489d-2*t31**2*t59
            t61 = log(t20)
            t62 = s1a2(D2_RA_RB)
            t63 = 2.970461675684382d+0*gammaaa*t1*t40*t7
            t64 = -1.6081984221457415d+1*t3*t62+t50
     +        -5.360661407152472d+0*s1f*t36
            t65 = -1.6202518231005728d+0*gammaaa*t1*t12*t4*t5*t64
            t66 = -4.725734484043338d-1*t1*t22*t7
            t67 = -4.050629557751432d-1*t1*t12*t24*t4*t5
            t68 = (3.109069086965489d-2*t21*(-1.6202518231005728d+0
     +        *gammaaa*t1*t17*t27*t4*(-2.6252159702607286d+0*gammaaa
     +        *t11*t24*t29*t5*t9-3.062751965304185d+0*gammaaa*t12*t28
     +        *t9+t67+t66)-6.563039925651821d-1*gammaaa*t11*t27*t30*t7
     +        *t9-7.656879913260463d-1*gammaaa*t19*t28*t7*t9
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t4*(t67+t66)
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t30*t32*t4*t57
     +        -1.6202518231005728d+0*gammaaa*t1*t26*t27*t32*t4
     +        -4.050629557751432d-1*t1*t17*t27*t30*t4
     +        +4.050629557751432d-1*t1*t19*t26*t4
     +        +1.8902937936173353d+0*gammaaa*t1*t17*t22*t27*t32
     +        -4.725734484043338d-1*t1*t17*t19*t22)
     +        -3.109069086965489d-2*t31*t33*t59)*wght
            t69 = -9.451468968086676d-1*t1*t22*t7
            t70 = -8.101259115502863d-1*t1*t12*t24*t4*t5
            t71 = -8.101259115502863d-1*t1*t17*t27*t32*t4
            t72 = t1**3
            t73 = 1/rhoa**7.0d+0
            t74 = -1.3126079851303643d+0*gammaaa*t11*t27*t32*t7*t9
            t75 = 3.109069086965489d-2*t21*(3.2815199628259106d-1*t11
     +        *t19*t7*t9+t74-5.316888702309608d-1*gammaaa*t12*t17*t27
     +        *t72*t73+t71+3.2405036462011455d+0*gammaaa*t1*t17*t32**2
     +        *t4*t57)*wght-3.109069086965489d-2*t33**2*t59*wght
            fnc(iq) = (3.109069086965489d-2*log(1.6202518231005725d+0
     +        *gammaaa*t1*t4*(t8+1.0d+0)/(t8+t13+1.0d+0)+1.0d+0)
     +        +5.0d-1*s1f*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (3.109069086965489d-2*t21*t31+5.0d-1*t15
     +        *t3-2.5d-1*t14*s1f)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 3.109069086965489d-2*t21*t33*wght+Cmat
     +        (iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t21*t35*wght+Cmat
     +        (iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (3.109069086965489d-2*t21*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*
     +        (1.7938975796781656d+1*t10*t12*t54*t9
     +        -5.250431940521457d+0*t10*t11*t29*t5*t51*t9+t56+t55+t53
     +        +t52+t49+t48+t47+t46)+t58+1.6202518231005728d+0*gammaaa
     +        *t1*t19*t4*(t52+t49+t48+t47+t46)+t43+t42+t41
     +        +3.3305176363733985d+0*gammaaa*t1*t17*t19*t40)
     +        -5.181781811609147d-3*t14*t61+t60+5.0d-1*t3*t39+t38+t37)
     +        *wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = (3.109069086965489d-2*t21*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -5.250431940521457d+0*t10*t11*t29*t5*t64*t9
     +        +1.6772213143332443d+1*t10*t12*t54*t9+t65+t63+t56+t55
     +        +t53+t49+t48+t46)+1.6202518231005728d+0*gammaaa*t1*t19
     +        *t4*(t65+t63+t49+t48+t46)+t58+t43+t42+t41
     +        +2.970461675684382d+0*gammaaa*t1*t17*t19*t40)+5.0d-1*t3
     +        *t62+5.181781811609147d-3*t14*t61+t60+t38+t37)*wght
     +        +Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t68+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = (3.109069086965489d-2*t21*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -5.250431940521457d+0*gammaaa*t11*t24*t29*t5*t9
     +        -6.12550393060837d+0*gammaaa*t12*t28*t9+t70+t69)
     +        -1.3126079851303643d+0*gammaaa*t11*t27*t30*t7*t9
     +        -1.5313759826520926d+0*gammaaa*t19*t28*t7*t9
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t4*(t70+t69)
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t30*t34*t4*t57
     +        -1.6202518231005728d+0*gammaaa*t1*t26*t27*t34*t4
     +        -8.101259115502863d-1*t1*t17*t27*t30*t4
     +        +8.101259115502863d-1*t1*t19*t26*t4
     +        +1.8902937936173353d+0*gammaaa*t1*t17*t22*t27*t34
     +        -9.451468968086676d-1*t1*t17*t19*t22)
     +        -3.109069086965489d-2*t31*t35*t59)*wght+Cmat2(iq,
     +        D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t68+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t75+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = 3.109069086965489d-2*t21*(
     +        -6.563039925651821d-1*gammaaa*t11*t27*t34*t7*t9
     +        +6.563039925651821d-1*t11*t19*t7*t9+t74
     +        -1.0633777404619216d+0*gammaaa*t12*t17*t27*t72*t73+t71
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t32*t34*t4*t57
     +        -4.050629557751432d-1*t1*t17*t27*t34*t4)*wght
     +        -3.109069086965489d-2*t33*t35*t59*wght+Cmat2(iq,
     +        D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t75+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 3.109069086965489d-2*t21*(
     +        -2.6252159702607286d+0*gammaaa*t11*t27*t34*t7*t9
     +        +1.3126079851303643d+0*t11*t19*t7*t9
     +        -2.1267554809238431d+0*gammaaa*t12*t17*t27*t72*t73
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t34**2*t4*t57
     +        -1.6202518231005728d+0*t1*t17*t27*t34*t4)*wght
     +        -3.109069086965489d-2*t35**2*t59*wght+Cmat2(iq,D2_GAB_GAB)
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
            t4 = 1/t3
            t5 = rhoa-rhob
            t6 = 1.0d+0-t4*t5
            t7 = t4*t5+1.0d+0
            t8 = 5.0d-1*t7**6.666666666666666d-1+5.0d-1*t6*
     +        *6.666666666666666d-1
            t9 = t8**3.0d+0
            t10 = gammabb+2.0d+0*gammaab+gammaaa
            t11 = 1/t3**2.3333333333333334d+0
            t12 = 1/t8**2
            t13 = 1/t9
            t14 = exp(-3.216396844291483d+1*t13*s1f*t4)
            t15 = t14-1.0d+0
            t16 = 1/t15
            t17 = 2.0413893780549558d+0*t1*t10*t11*t12*t16
            t18 = t17+1.0d+0
            t19 = t1**2
            t20 = t10**2
            t21 = 1/t3**4.666666666666667d+0
            t22 = 1/t8**4
            t23 = 1/t15**2
            t24 = 4.1672705928355996d+0*t19*t20*t21*t22*t23+t17+1.0d+0
            t25 = 1/t24
            t26 = 2.0413893780549558d+0*t1*t10*t11*t12*t18*t25+1.0d+0
            t27 = log(t26)
            t28 = 1/t3**2
            t29 = -1.0d+0*s1f*t28
            t30 = s1a(D1_RA)
            t31 = 1/t26
            t32 = t28*t5
            t33 = -t4
            t34 = t33+t32
            t35 = 1/t6**3.333333333333333d-1
            t36 = -t28*t5
            t37 = t4+t36
            t38 = 1/t7**3.333333333333333d-1
            t39 = 3.333333333333333d-1*t37*t38+3.333333333333333d-1
     +        *t34*t35
            t40 = 1/t8**3
            t41 = 1/t3**3.3333333333333337d+0
            t42 = -4.763241882128231d+0*t1*t10*t12*t18*t25*t41
            t43 = -4.0827787561099116d+0*t1*t10*t11*t16*t39*t40
            t44 = -4.763241882128231d+0*t1*t10*t12*t16*t41
            t45 = 1/t8**4.0d+0
            t46 = 3.216396844291483d+1*t13*s1f*t28
            t47 = t46+9.649190532874448d+1*s1f*t39*t4*t45
     +        -3.216396844291483d+1*t13*t30*t4
            t48 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t47
            t49 = t48+t44+t43
            t50 = 1/t24**2
            t51 = 1/t8**5
            t52 = 1/t3**5.666666666666667d+0
            t53 = -1.9447262766566134d+1*t19*t20*t22*t23*t52
            t54 = 1/t15**3
            t55 = -8.334541185671199d+0*t14*t19*t20*t21*t22*t47*t54
     +        +t53-1.6669082371342397d+1*t19*t20*t21*t23*t39*t51+t48
     +        +t44+t43
            t56 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t55
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*t49+t42
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t39*t40
            t57 = t8**2.0d+0
            t58 = s1a(D1_RB)
            t59 = t4+t32
            t60 = t36+t33
            t61 = 3.333333333333333d-1*t38*t60+3.333333333333333d-1
     +        *t35*t59
            t62 = -4.0827787561099116d+0*t1*t10*t11*t16*t40*t61
            t63 = 9.649190532874448d+1*s1f*t4*t45*t61
     +        -3.216396844291483d+1*t13*t4*t58+t46
            t64 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t63
            t65 = t64+t62+t44
            t66 = t64-8.334541185671199d+0*t14*t19*t20*t21*t22*t54*t63
     +        +t62-1.6669082371342397d+1*t19*t20*t21*t23*t51*t61+t53+t44
            t67 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t66
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t40*t61+t42
            t68 = 8.334541185671199d+0*t10*t19*t21*t22*t23
     +        +2.0413893780549558d+0*t1*t11*t12*t16
            t69 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t68
     +        +4.1672705928355996d+0*t10*t16*t19*t21*t22*t25
     +        +2.0413893780549558d+0*t1*t11*t12*t18*t25
            t70 = 3.109069086965489d-2*t31*t69*t9*wght
            t71 = 1.6669082371342397d+1*t10*t19*t21*t22*t23
     +        +4.0827787561099116d+0*t1*t11*t12*t16
            t72 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t71
     +        +8.334541185671199d+0*t10*t16*t19*t21*t22*t25
     +        +4.0827787561099116d+0*t1*t11*t12*t18*t25
            t73 = 1/t3**3
            t74 = 2.0d+0*s1f*t73
            t75 = s1a2(D2_RA_RA)
            t76 = t39**2
            t77 = 1/t6**1.3333333333333333d+0
            t78 = -2*t5*t73
            t79 = 2*t28
            t80 = 1/t7**1.3333333333333333d+0
            t81 = 2*t5*t73
            t82 = -2*t28
            t83 = 3.333333333333333d-1*t38*(t82+t81)
     +        -1.111111111111111d-1*t37**2*t80+3.333333333333333d-1
     +        *t35*(t79+t78)-1.111111111111111d-1*t34**2*t77
            t84 = 1/t3**4.333333333333334d+0
            t85 = 1.5877472940427437d+1*t1*t10*t12*t18*t25*t84
            t86 = t47**2
            t87 = exp(-6.432793688582966d+1*t13*s1f*t4)
            t88 = 4.0827787561099116d+0*t1*t10*t11*t12*t54*t86*t87
            t89 = 1.2248336268329735d+1*t1*t10*t11*t16*t22*t76
            t90 = -4.0827787561099116d+0*t1*t10*t11*t16*t40*t83
            t91 = 1.9052967528512924d+1*t1*t10*t16*t39*t40*t41
            t92 = 1.5877472940427437d+1*t1*t10*t12*t16*t84
            t93 = 8.165557512219823d+0*t1*t10*t11*t14*t23*t39*t40*t47
            t94 = 9.526483764256461d+0*t1*t10*t12*t14*t23*t41*t47
            t95 = 1/t8**5.0d+0
            t96 = -6.432793688582966d+1*t13*s1f*t73
            t97 = t96-3.859676213149779d+2*s1f*t4*t76*t95
     +        +9.649190532874448d+1*s1f*t4*t45*t83
     +        -3.216396844291483d+1*t13*t4*t75+1.9298381065748896d+2
     +        *t30*t39*t4*t45-1.9298381065748896d+2*s1f*t28*t39*t45
     +        +6.432793688582966d+1*t13*t28*t30
            t98 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t97
            t99 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t86
            t100 = 1/t15**4
            t101 = 1/t8**6
            t102 = 1.1020115567720809d+2*t19*t20*t22*t23/t3*
     +        *6.666666666666667d+0
            t103 = 1/t24**3
            t104 = 1/t26**2
            t105 = s1a2(D2_RA_RB)
            t106 = -1.111111111111111d-1*t37*t60*t80
     +        -1.111111111111111d-1*t34*t59*t77+6.666666666666666d-1
     +        *t38*t5*t73-6.666666666666666d-1*t35*t5*t73
            t107 = 4.0827787561099116d+0*t1*t10*t11*t12*t47*t54*t63*t87
            t108 = 1.2248336268329735d+1*t1*t10*t11*t16*t22*t39*t61
            t109 = -4.0827787561099116d+0*t1*t10*t106*t11*t16*t40
            t110 = 9.526483764256461d+0*t1*t10*t16*t40*t41*t61
            t111 = 9.526483764256461d+0*t1*t10*t16*t39*t40*t41
            t112 = 4.0827787561099116d+0*t1*t10*t11*t14*t23*t40*t47*t61
            t113 = 4.763241882128231d+0*t1*t10*t12*t14*t23*t41*t47
            t114 = t96-3.859676213149779d+2*s1f*t39*t4*t61*t95
     +        +9.649190532874448d+1*t30*t4*t45*t61
     +        -9.649190532874448d+1*s1f*t28*t45*t61
     +        +9.649190532874448d+1*t39*t4*t45*t58
     +        +3.216396844291483d+1*t13*t28*t58+9.649190532874448d+1
     +        *t106*s1f*t4*t45-9.649190532874448d+1*s1f*t28*t39*t45
     +        -3.216396844291483d+1*t105*t13*t4+3.216396844291483d+1
     +        *t13*t28*t30
            t115 = -2.0413893780549558d+0*t1*t10*t11*t114*t12*t14*t23
            t116 = 4.0827787561099116d+0*t1*t10*t11*t14*t23*t39*t40*t63
            t117 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t47*t63
            t118 = 4.763241882128231d+0*t1*t10*t12*t14*t23*t41*t63
            t119 = s1a2(D2_RB_RB)
            t120 = t61**2
            t121 = 3.333333333333333d-1*t35*(t82+t78)
     +        +3.333333333333333d-1*t38*(t81+t79)-1.111111111111111d-1
     +        *t60**2*t80-1.111111111111111d-1*t59**2*t77
            t122 = t63**2
            t123 = 4.0827787561099116d+0*t1*t10*t11*t12*t122*t54*t87
            t124 = 1.2248336268329735d+1*t1*t10*t11*t120*t16*t22
            t125 = -4.0827787561099116d+0*t1*t10*t11*t121*t16*t40
            t126 = 1.9052967528512924d+1*t1*t10*t16*t40*t41*t61
            t127 = 8.165557512219823d+0*t1*t10*t11*t14*t23*t40*t61*t63
            t128 = 9.526483764256461d+0*t1*t10*t12*t14*t23*t41*t63
            t129 = t96-3.859676213149779d+2*t120*s1f*t4*t95
     +        +1.9298381065748896d+2*t4*t45*t58*t61
     +        -1.9298381065748896d+2*s1f*t28*t45*t61
     +        +6.432793688582966d+1*t13*t28*t58+9.649190532874448d+1
     +        *t121*s1f*t4*t45-3.216396844291483d+1*t119*t13*t4
            t130 = -2.0413893780549558d+0*t1*t10*t11*t12*t129*t14*t23
            t131 = -2.0413893780549558d+0*t1*t10*t11*t12*t122*t14*t23
            t132 = 4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t68
            t133 = -4.763241882128231d+0*t1*t12*t18*t25*t41
            t134 = -9.723631383283067d+0*t10*t16*t19*t22*t25*t52
            t135 = -4.0827787561099116d+0*t1*t11*t16*t39*t40
            t136 = -4.763241882128231d+0*t1*t12*t16*t41
            t137 = -2.0413893780549558d+0*t1*t11*t12*t14*t23*t47
            t138 = -3.889452553313227d+1*t10*t19*t22*t23*t52
            t139 = (t31*(3.109069086965489d-2*(4.0827787561099116d+0
     +        *t1*t10*t103*t11*t12*t18*t55*t68-2.0413893780549558d+0
     +        *t1*t10*t11*t12*t49*t50*t68+4.0827787561099116d+0*t1*t10
     +        *t11*t18*t39*t40*t50*t68-4.1672705928355996d+0*t10*t16
     +        *t19*t21*t22*t50*t55-2.0413893780549558d+0*t1*t11*t12
     +        *t18*t50*t55-2.0413893780549558d+0*t1*t10*t11*t12*t18
     +        *t50*(-1.6669082371342397d+1*t10*t14*t19*t21*t22*t47*t54
     +        -3.3338164742684795d+1*t10*t19*t21*t23*t39*t51+t138+t137
     +        +t136+t135)-8.334541185671199d+0*t10*t16*t19*t21*t25*t39
     +        *t51+2.0413893780549558d+0*t1*t11*t12*t25*t49
     +        -4.0827787561099116d+0*t1*t11*t18*t25*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*(t137+t136+t135)
     +        *t25+t134+t133+t132)*t9+9.327207260896466d-2*t39*t57*t69
     +        )-3.109069086965489d-2*t104*t56*t69*t9)*wght
            t140 = 4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t71
            t141 = -9.526483764256461d+0*t1*t12*t18*t25*t41
            t142 = -1.9447262766566134d+1*t10*t16*t19*t22*t25*t52
            t143 = -8.165557512219823d+0*t1*t11*t16*t39*t40
            t144 = -9.526483764256461d+0*t1*t12*t16*t41
            t145 = -4.0827787561099116d+0*t1*t11*t12*t14*t23*t47
            t146 = -7.778905106626454d+1*t10*t19*t22*t23*t52
            t147 = -4.0827787561099116d+0*t1*t11*t16*t40*t61
            t148 = -2.0413893780549558d+0*t1*t11*t12*t14*t23*t63
            t149 = (t31*(3.109069086965489d-2*(4.0827787561099116d+0
     +        *t1*t10*t103*t11*t12*t18*t66*t68-2.0413893780549558d+0
     +        *t1*t10*t11*t12*t50*t65*t68+4.0827787561099116d+0*t1*t10
     +        *t11*t18*t40*t50*t61*t68-4.1672705928355996d+0*t10*t16
     +        *t19*t21*t22*t50*t66-2.0413893780549558d+0*t1*t11*t12
     +        *t18*t50*t66+2.0413893780549558d+0*t1*t11*t12*t25*t65
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -1.6669082371342397d+1*t10*t14*t19*t21*t22*t54*t63
     +        -3.3338164742684795d+1*t10*t19*t21*t23*t51*t61+t148+t147
     +        +t138+t136)-8.334541185671199d+0*t10*t16*t19*t21*t25*t51
     +        *t61-4.0827787561099116d+0*t1*t11*t18*t25*t40*t61
     +        +2.0413893780549558d+0*t1*t10*t11*t12*(t148+t147+t136)
     +        *t25+t134+t133+t132)*t9+9.327207260896466d-2*t57*t61*t69
     +        )-3.109069086965489d-2*t104*t67*t69*t9)*wght
            t150 = -8.165557512219823d+0*t1*t11*t16*t40*t61
            t151 = -4.0827787561099116d+0*t1*t11*t12*t14*t23*t63
            t152 = -4.0827787561099116d+0*t1*t11*t12*t18*t50*t68
            t153 = t1**3
            t154 = 1/t3**7.0d+0
            t155 = -8.334541185671199d+0*t10*t16*t19*t21*t22*t50*t68
            t156 = 3.109069086965489d-2*t31*(4.0827787561099116d+0*t1
     +        *t10*t103*t11*t12*t18*t68**2-1.7014043847390745d+1*t10
     +        *t101*t153*t154*t18*t23*t50+8.334541185671199d+0*t16*t19
     +        *t21*t22*t25+t155+t152)*t9*wght-3.109069086965489d-2
     +        *t104*t69**2*t9*wght
            t157 = 3.109069086965489d-2*t31*(4.0827787561099116d+0*t1
     +        *t10*t103*t11*t12*t18*t68*t71-4.1672705928355996d+0*t10
     +        *t16*t19*t21*t22*t50*t71-2.0413893780549558d+0*t1*t11
     +        *t12*t18*t50*t71-3.402808769478149d+1*t10*t101*t153*t154
     +        *t18*t23*t50+1.6669082371342397d+1*t16*t19*t21*t22*t25
     +        +t155+t152)*t9*wght-3.109069086965489d-2*t104*t69*t72*t9
     +        *wght
            fnc(iq) = (3.109069086965489d-2*t27*t9+1.0d+0*s1f*t4)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = (3.109069086965489d-2*t31*t56*t9
     +        +9.327207260896466d-2*t27*t39*t57+1.0d+0*t30*t4+t29)
     +        *wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (3.109069086965489d-2*t31*t67*t9
     +        +9.327207260896466d-2*t27*t57*t61+1.0d+0*t4*t58+t29)
     +        *wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t70+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t31*t72*t9*wght
     +        +Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = t70+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(t99+t98
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t54*t97+t94
     +        +t93+t92+t91+t90+t89+t88+2.50036235570136d+1*t100*t19
     +        *t20*t21*t22*t86*t87-8.334541185671199d+0*t14*t19*t20
     +        *t21*t22*t54*t86-1.6669082371342397d+1*t19*t20*t21*t23
     +        *t51*t83+8.334541185671199d+1*t101*t19*t20*t21*t23*t76
     +        +7.778905106626454d+1*t14*t19*t20*t22*t47*t52*t54
     +        +6.667632948536959d+1*t14*t19*t20*t21*t39*t47*t51*t54
     +        +1.5557810213252907d+2*t19*t20*t23*t39*t51*t52+t102)
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t99+t98+t94
     +        +t93+t92+t91+t90+t89+t88)+t85-4.0827787561099116d+0*t1
     +        *t10*t11*t18*t25*t40*t83+1.2248336268329735d+1*t1*t10
     +        *t11*t18*t22*t25*t76+4.0827787561099116d+0*t1*t10*t103
     +        *t11*t12*t18*t55**2-4.0827787561099116d+0*t1*t10*t11*t12
     +        *t49*t50*t55+9.526483764256461d+0*t1*t10*t12*t18*t41*t50
     +        *t55+8.165557512219823d+0*t1*t10*t11*t18*t39*t40*t50*t55
     +        -9.526483764256461d+0*t1*t10*t12*t25*t41*t49
     +        -8.165557512219823d+0*t1*t10*t11*t25*t39*t40*t49
     +        +1.9052967528512924d+1*t1*t10*t18*t25*t39*t40*t41)
     +        +1.8654414521792934d-1*t39*t56*t57)-3.109069086965489d-2
     +        *t104*t56**2*t9+9.327207260896466d-2*t27*t57*t83
     +        +1.8654414521792934d-1*t27*t76*t8+1.0d+0*t4*t75+t74
     +        -2.0d+0*t28*t30)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(t92
     +        +2.50036235570136d+1*t100*t19*t20*t21*t22*t47*t63*t87
     +        +3.889452553313227d+1*t14*t19*t20*t22*t52*t54*t63
     +        +3.3338164742684795d+1*t14*t19*t20*t21*t39*t51*t54*t63
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t47*t54*t63
     +        +3.3338164742684795d+1*t14*t19*t20*t21*t47*t51*t54*t61
     +        +7.778905106626454d+1*t19*t20*t23*t51*t52*t61
     +        +8.334541185671199d+1*t101*t19*t20*t21*t23*t39*t61
     +        +3.889452553313227d+1*t14*t19*t20*t22*t47*t52*t54
     +        -8.334541185671199d+0*t114*t14*t19*t20*t21*t22*t54
     +        +7.778905106626454d+1*t19*t20*t23*t39*t51*t52
     +        -1.6669082371342397d+1*t106*t19*t20*t21*t23*t51+t118
     +        +t117+t116+t115+t113+t112+t111+t110+t109+t108+t107+t102)
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t92+t118+t117
     +        +t116+t115+t113+t112+t111+t110+t109+t108+t107)+t85
     +        +4.0827787561099116d+0*t1*t10*t103*t11*t12*t18*t55*t66
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t49*t50*t66
     +        +4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t66
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t39*t40*t50*t66
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t50*t55*t65
     +        -4.763241882128231d+0*t1*t10*t12*t25*t41*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t25*t39*t40*t65
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t55*t61
     +        -4.0827787561099116d+0*t1*t10*t11*t25*t40*t49*t61
     +        +9.526483764256461d+0*t1*t10*t18*t25*t40*t41*t61
     +        +1.2248336268329735d+1*t1*t10*t11*t18*t22*t25*t39*t61
     +        +4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t55
     +        -4.763241882128231d+0*t1*t10*t12*t25*t41*t49
     +        +9.526483764256461d+0*t1*t10*t18*t25*t39*t40*t41
     +        -4.0827787561099116d+0*t1*t10*t106*t11*t18*t25*t40)
     +        +9.327207260896466d-2*t39*t57*t67+9.327207260896466d-2
     +        *t56*t57*t61)-3.109069086965489d-2*t104*t56*t67*t9
     +        +1.8654414521792934d-1*t27*t39*t61*t8+t74+t28*(-1.0d+0
     +        *t58-1.0d+0*t30)+9.327207260896466d-2*t106*t27*t57
     +        +1.0d+0*t105*t4)*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(t92
     +        +2.50036235570136d+1*t100*t122*t19*t20*t21*t22*t87
     +        +6.667632948536959d+1*t14*t19*t20*t21*t51*t54*t61*t63
     +        +7.778905106626454d+1*t14*t19*t20*t22*t52*t54*t63
     +        +1.5557810213252907d+2*t19*t20*t23*t51*t52*t61
     +        -8.334541185671199d+0*t129*t14*t19*t20*t21*t22*t54
     +        -8.334541185671199d+0*t122*t14*t19*t20*t21*t22*t54
     +        -1.6669082371342397d+1*t121*t19*t20*t21*t23*t51
     +        +8.334541185671199d+1*t101*t120*t19*t20*t21*t23+t131
     +        +t130+t128+t127+t126+t125+t124+t123+t102)
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t92+t131+t130
     +        +t128+t127+t126+t125+t124+t123)+t85
     +        +4.0827787561099116d+0*t1*t10*t103*t11*t12*t18*t66**2
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t50*t65*t66
     +        +8.165557512219823d+0*t1*t10*t11*t18*t40*t50*t61*t66
     +        +9.526483764256461d+0*t1*t10*t12*t18*t41*t50*t66
     +        -8.165557512219823d+0*t1*t10*t11*t25*t40*t61*t65
     +        -9.526483764256461d+0*t1*t10*t12*t25*t41*t65
     +        +1.9052967528512924d+1*t1*t10*t18*t25*t40*t41*t61
     +        -4.0827787561099116d+0*t1*t10*t11*t121*t18*t25*t40
     +        +1.2248336268329735d+1*t1*t10*t11*t120*t18*t22*t25)
     +        +1.8654414521792934d-1*t57*t61*t67)-3.109069086965489d-2
     +        *t104*t67**2*t9+1.8654414521792934d-1*t120*t27*t8+t74
     +        -2.0d+0*t28*t58+9.327207260896466d-2*t121*t27*t57+1.0d+0
     +        *t119*t4)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t139+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = (t31*(3.109069086965489d-2*
     +        (4.0827787561099116d+0*t1*t10*t103*t11*t12*t18*t55*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t49*t50*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t39*t40*t50*t71
     +        -8.334541185671199d+0*t10*t16*t19*t21*t22*t50*t55
     +        -4.0827787561099116d+0*t1*t11*t12*t18*t50*t55
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -3.3338164742684795d+1*t10*t14*t19*t21*t22*t47*t54
     +        -6.667632948536959d+1*t10*t19*t21*t23*t39*t51+t146+t145
     +        +t144+t143)-1.6669082371342397d+1*t10*t16*t19*t21*t25
     +        *t39*t51+4.0827787561099116d+0*t1*t11*t12*t25*t49
     +        -8.165557512219823d+0*t1*t11*t18*t25*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*(t145+t144+t143)
     +        *t25+t142+t141+t140)*t9+9.327207260896466d-2*t39*t57*t72
     +        )-3.109069086965489d-2*t104*t56*t72*t9)*wght+Cmat2(iq,
     +        D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t139+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t149+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = (t31*(3.109069086965489d-2*
     +        (4.0827787561099116d+0*t1*t10*t103*t11*t12*t18*t66*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t50*t65*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t61*t71
     +        -8.334541185671199d+0*t10*t16*t19*t21*t22*t50*t66
     +        -4.0827787561099116d+0*t1*t11*t12*t18*t50*t66
     +        +4.0827787561099116d+0*t1*t11*t12*t25*t65
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -3.3338164742684795d+1*t10*t14*t19*t21*t22*t54*t63
     +        -6.667632948536959d+1*t10*t19*t21*t23*t51*t61+t151+t150
     +        +t146+t144)-1.6669082371342397d+1*t10*t16*t19*t21*t25
     +        *t51*t61-8.165557512219823d+0*t1*t11*t18*t25*t40*t61
     +        +2.0413893780549558d+0*t1*t10*t11*t12*(t151+t150+t144)
     +        *t25+t142+t141+t140)*t9+9.327207260896466d-2*t57*t61*t72
     +        )-3.109069086965489d-2*t104*t67*t72*t9)*wght+Cmat2(iq,
     +        D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t149+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t156+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t157+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t156+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = -3.109069086965489d-2*t104*t72**2
     +        *t9*wght+3.109069086965489d-2*t31*(4.0827787561099116d+0
     +        *t1*t10*t103*t11*t12*t18*t71**2-1.6669082371342397d+1
     +        *t10*t16*t19*t21*t22*t50*t71-8.165557512219823d+0*t1*t11
     +        *t12*t18*t50*t71-6.805617538956298d+1*t10*t101*t153*t154
     +        *t18*t23*t50+3.3338164742684795d+1*t16*t19*t21*t22*t25)
     +        *t9*wght+Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = t157+Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t156+Cmat2(iq,D2_GBB_GBB)
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
            t3 = 1/rhoa
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = exp(-6.432793688582967d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 3.240503646201146d+0*gammaaa*t1*t4*t7
            t9 = t8+1.0d+0
            t10 = t1**2
            t11 = gammaaa**2
            t12 = 1/rhoa**4.666666666666667d+0
            t13 = 1/t6**2
            t14 = t8+1.0500863881042923d+1*t10*t11*t12*t13+1.0d+0
            t15 = 1/t14
            t16 = 3.240503646201146d+0*gammaaa*t1*t15*t4*t9+1.0d+0
            t17 = 1/rhoa**2
            t18 = s1a(D1_RA)
            t19 = 1/t16
            t20 = 1/rhoa**3.3333333333333337d+0
            t21 = -7.561175174469341d+0*gammaaa*t1*t20*t7
            t22 = 6.432793688582967d+1*t17*s1f-6.432793688582967d+1
     +        *t18*t3
            t23 = -3.240503646201146d+0*gammaaa*t1*t13*t22*t4*t5
            t24 = t23+t21
            t25 = 1/t14**2
            t26 = 1/rhoa**5.666666666666667d+0
            t27 = 1/t6**3
            t28 = -2.1001727762085845d+1*t10*t11*t12*t22*t27*t5
     +        -4.9004031444866974d+1*t10*t11*t13*t26+t23+t21
            t29 = -3.240503646201146d+0*gammaaa*t1*t25*t28*t4*t9
     +        -7.561175174469341d+0*gammaaa*t1*t15*t20*t9
     +        +3.240503646201146d+0*gammaaa*t1*t15*t24*t4
            t30 = 3.240503646201146d+0*t1*t4*t7+2.1001727762085845d+1
     +        *gammaaa*t10*t12*t13
            t31 = -3.240503646201146d+0*gammaaa*t1*t25*t30*t4*t9
     +        +3.240503646201146d+0*t1*t15*t4*t9+1.0500863881042925d+1
     +        *gammaaa*t10*t12*t15*t7
            t32 = 1/rhoa**3
            t33 = s1a2(D2_RA_RA)
            t34 = 1/rhoa**4.333333333333334d+0
            t35 = t22**2
            t36 = exp(-1.2865587377165935d+2*s1f*t3)
            t37 = 6.481007292402292d+0*gammaaa*t1*t27*t35*t36*t4
            t38 = 2.520391724823114d+1*gammaaa*t1*t34*t7
            t39 = 1.5122350348938685d+1*gammaaa*t1*t13*t20*t22*t5
            t40 = -3.240503646201146d+0*gammaaa*t1*t13*t35*t4*t5
            t41 = -6.432793688582967d+1*t3*t33-1.2865587377165935d+2
     +        *s1f*t32+1.2865587377165935d+2*t17*t18
            t42 = -3.240503646201146d+0*gammaaa*t1*t13*t4*t41*t5
            t43 = 1/t14**3
            t44 = 1/t16**2
            t45 = -7.561175174469341d+0*t1*t20*t7
            t46 = -3.240503646201146d+0*t1*t13*t22*t4*t5
            fnc(iq) = (1.0d+0*s1f*t3+1.5545345434827443d-2*log(t16))
     +        *wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t18*t3+1.5545345434827443d-2*t19
     +        *t29-1.0d+0*t17*s1f)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.5545345434827443d-2*t19*t31*wght+Cmat
     +        (iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammaaa*t1*t25*t4*
     +        (6.300518328625754d+1*t10*t11*t12*t35*t36/t6**4
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t41*t5
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t35*t5
     +        +1.960161257794679d+2*t10*t11*t22*t26*t27*t5+t42+t40+t39
     +        +t38+t37+2.7768951152091287d+2*t10*t11*t13/rhoa*
     +        *6.666666666666667d+0)*t9+6.481007292402292d+0*gammaaa
     +        *t1*t28**2*t4*t43*t9+2.520391724823114d+1*gammaaa*t1*t15
     +        *t34*t9+1.5122350348938685d+1*gammaaa*t1*t20*t25*t28*t9
     +        +3.240503646201146d+0*gammaaa*t1*t15*t4*(t42+t40+t39+t38
     +        +t37)-6.481007292402292d+0*gammaaa*t1*t24*t25*t28*t4
     +        -1.5122350348938685d+1*gammaaa*t1*t15*t20*t24)
     +        -1.5545345434827443d-2*t29**2*t44+1.0d+0*t3*t33+2.0d+0
     +        *s1f*t32-2.0d+0*t17*t18)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammaaa*t1*t25*t4*(
     +        -4.200345552417169d+1*gammaaa*t10*t12*t22*t27*t5+t46+t45
     +        -9.800806288973395d+1*gammaaa*t10*t13*t26)*t9
     +        +6.481007292402292d+0*gammaaa*t1*t28*t30*t4*t43*t9
     +        -3.240503646201146d+0*t1*t25*t28*t4*t9
     +        +7.561175174469341d+0*gammaaa*t1*t20*t25*t30*t9
     +        -7.561175174469341d+0*t1*t15*t20*t9
     +        -1.0500863881042925d+1*gammaaa*t10*t12*t25*t28*t7
     +        -2.450201572243349d+1*gammaaa*t10*t15*t26*t7
     +        +3.240503646201146d+0*gammaaa*t1*t15*t4*(t46+t45)
     +        -3.240503646201146d+0*gammaaa*t1*t24*t25*t30*t4
     +        +3.240503646201146d+0*t1*t15*t24*t4)
     +        -1.5545345434827443d-2*t29*t31*t44)*wght+Cmat2(iq,
     +        D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.5545345434827443d-2*t19*
     +        (6.481007292402292d+0*gammaaa*t1*t30**2*t4*t43*t9
     +        -6.481007292402292d+0*t1*t25*t30*t4*t9
     +        -6.805617538956302d+1*gammaaa*t1**3*t13*t25*t9/rhoa*
     +        *7.0d+0-2.100172776208585d+1*gammaaa*t10*t12*t25*t30*t7
     +        +2.100172776208585d+1*t10*t12*t15*t7)*wght
     +        -1.5545345434827443d-2*t31**2*t44*wght+Cmat2(iq,
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
            t3 = 1/rhob
            t4 = 1/rhob**2.3333333333333334d+0
            t5 = exp(-6.432793688582967d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 3.240503646201146d+0*gammabb*t1*t4*t7
            t9 = t8+1.0d+0
            t10 = t1**2
            t11 = gammabb**2
            t12 = 1/rhob**4.666666666666667d+0
            t13 = 1/t6**2
            t14 = t8+1.0500863881042923d+1*t10*t11*t12*t13+1.0d+0
            t15 = 1/t14
            t16 = 3.240503646201146d+0*gammabb*t1*t15*t4*t9+1.0d+0
            t17 = 1/rhob**2
            t18 = s1a(D1_RB)
            t19 = 1/t16
            t20 = 1/rhob**3.3333333333333337d+0
            t21 = -7.561175174469341d+0*gammabb*t1*t20*t7
            t22 = 6.432793688582967d+1*t17*s1f-6.432793688582967d+1
     +        *t18*t3
            t23 = -3.240503646201146d+0*gammabb*t1*t13*t22*t4*t5
            t24 = t23+t21
            t25 = 1/t14**2
            t26 = 1/rhob**5.666666666666667d+0
            t27 = 1/t6**3
            t28 = -2.1001727762085845d+1*t10*t11*t12*t22*t27*t5
     +        -4.9004031444866974d+1*t10*t11*t13*t26+t23+t21
            t29 = -3.240503646201146d+0*gammabb*t1*t25*t28*t4*t9
     +        -7.561175174469341d+0*gammabb*t1*t15*t20*t9
     +        +3.240503646201146d+0*gammabb*t1*t15*t24*t4
            t30 = 3.240503646201146d+0*t1*t4*t7+2.1001727762085845d+1
     +        *gammabb*t10*t12*t13
            t31 = -3.240503646201146d+0*gammabb*t1*t25*t30*t4*t9
     +        +3.240503646201146d+0*t1*t15*t4*t9+1.0500863881042925d+1
     +        *gammabb*t10*t12*t15*t7
            t32 = 1/rhob**3
            t33 = s1a2(D2_RB_RB)
            t34 = 1/rhob**4.333333333333334d+0
            t35 = t22**2
            t36 = exp(-1.2865587377165935d+2*s1f*t3)
            t37 = 6.481007292402292d+0*gammabb*t1*t27*t35*t36*t4
            t38 = 2.520391724823114d+1*gammabb*t1*t34*t7
            t39 = 1.5122350348938685d+1*gammabb*t1*t13*t20*t22*t5
            t40 = -3.240503646201146d+0*gammabb*t1*t13*t35*t4*t5
            t41 = -6.432793688582967d+1*t3*t33-1.2865587377165935d+2
     +        *s1f*t32+1.2865587377165935d+2*t17*t18
            t42 = -3.240503646201146d+0*gammabb*t1*t13*t4*t41*t5
            t43 = 1/t14**3
            t44 = 1/t16**2
            t45 = -7.561175174469341d+0*t1*t20*t7
            t46 = -3.240503646201146d+0*t1*t13*t22*t4*t5
            fnc(iq) = (1.0d+0*s1f*t3+1.5545345434827443d-2*log(t16))
     +        *wght+fnc(iq)
            Amat(iq,D1_RB) = (1.0d+0*t18*t3+1.5545345434827443d-2*t19
     +        *t29-1.0d+0*t17*s1f)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.5545345434827443d-2*t19*t31*wght+Cmat
     +        (iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammabb*t1*t25*t4*
     +        (6.300518328625754d+1*t10*t11*t12*t35*t36/t6**4
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t41*t5
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t35*t5
     +        +1.960161257794679d+2*t10*t11*t22*t26*t27*t5+t42+t40+t39
     +        +t38+t37+2.7768951152091287d+2*t10*t11*t13/rhob*
     +        *6.666666666666667d+0)*t9+6.481007292402292d+0*gammabb
     +        *t1*t28**2*t4*t43*t9+2.520391724823114d+1*gammabb*t1*t15
     +        *t34*t9+1.5122350348938685d+1*gammabb*t1*t20*t25*t28*t9
     +        +3.240503646201146d+0*gammabb*t1*t15*t4*(t42+t40+t39+t38
     +        +t37)-6.481007292402292d+0*gammabb*t1*t24*t25*t28*t4
     +        -1.5122350348938685d+1*gammabb*t1*t15*t20*t24)
     +        -1.5545345434827443d-2*t29**2*t44+1.0d+0*t3*t33+2.0d+0
     +        *s1f*t32-2.0d+0*t17*t18)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammabb*t1*t25*t4*(
     +        -4.200345552417169d+1*gammabb*t10*t12*t22*t27*t5+t46+t45
     +        -9.800806288973395d+1*gammabb*t10*t13*t26)*t9
     +        +6.481007292402292d+0*gammabb*t1*t28*t30*t4*t43*t9
     +        -3.240503646201146d+0*t1*t25*t28*t4*t9
     +        +7.561175174469341d+0*gammabb*t1*t20*t25*t30*t9
     +        -7.561175174469341d+0*t1*t15*t20*t9
     +        -1.0500863881042925d+1*gammabb*t10*t12*t25*t28*t7
     +        -2.450201572243349d+1*gammabb*t10*t15*t26*t7
     +        +3.240503646201146d+0*gammabb*t1*t15*t4*(t46+t45)
     +        -3.240503646201146d+0*gammabb*t1*t24*t25*t30*t4
     +        +3.240503646201146d+0*t1*t15*t24*t4)
     +        -1.5545345434827443d-2*t29*t31*t44)*wght+Cmat2(iq,
     +        D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5545345434827443d-2*t19*
     +        (6.481007292402292d+0*gammabb*t1*t30**2*t4*t43*t9
     +        -6.481007292402292d+0*t1*t25*t30*t4*t9
     +        -6.805617538956302d+1*gammabb*t1**3*t13*t25*t9/rhob*
     +        *7.0d+0-2.100172776208585d+1*gammabb*t10*t12*t25*t30*t7
     +        +2.100172776208585d+1*t10*t12*t15*t7)*wght
     +        -1.5545345434827443d-2*t31**2*t44*wght+Cmat2(iq,
     +        D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_mpbe_d2 
C>
C> \brief Evaluate the nwxcm_c_mpbe functional [1]
C>
C> \f{eqnarray*}{
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\
C>   {\it t_2} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_\alpha , \rho_\beta\right)\\\\
C>   {\it t_3} &=& \rho_\beta+\rho_\alpha\\\\
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\
C>   {\it t_5} &=& \rho_\alpha-\rho_\beta\\\\
C>   {\it t_6} &=& 0.5\,\left({\it t_5}\,{\it t_4}
C>    +1.0\right)^{{{2}\over{3}}}+0.5\,\left(1.0-{\it t_5}\,{
C>    \it t_4}\right)^{{{2}\over{3}}}\\\\
C>   {\it t_7} &=& {\it t_6}^{3.0}\\\\
C>   {\it t_8} &=& {{1}\over{{\it t_3}^{2.333333333333333}}}\\\\
C>   {\it t_9} &=& {{1}\over{{\it t_6}^2}}\\\\
C>   {\it t_{10}} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\
C>   {\it t_{11}} &=& {{1}\over{e^{{{32.16396844291483\,{\it t_2}
C>    \,{\it t_4}}\over{{\it t_7}}}}}}-1.0\\\\
C>   {\it t_{12}} &=& {{2.041389378054956\,{\it t_1}\,{\it t_8}\,{
C>    \it t_9}\,{\it t_{10}}}\over{{\it t_{11}}}}\\\\
C>   {\it t_{13}} &=& {\it t_1}^2\\\\
C>   {\it t_{14}} &=& {\it nwxc\_c\_pw91lda}\left({\it t_1} 
C>    , \rho_s , 0.0\right)\\\\
C>   {\it t_{15}} &=& {{1}\over{\rho_s}}\\\\
C>   {\it t_{16}} &=& {{1}\over{\rho_s^{2.333333333333333}}}\\\\
C>   {\it t_{17}} &=& {{1}\over{e^{64.32793688582967\,{\it t_{14}}
C>    \,{\it t_{15}}}}}-1.0\\\\
C>   {\it t_{18}} &=& {{3.240503646201146\,{\it t_1}\,{\it t_{16}}
C>    \,\sigma_{ss}}\over{{\it t_{17}}}}\\\\
C>   f &=& 0.03109069086965489\,{\it t_7}\,\log 
C>    \left({{2.041389378054956\,{\it t_1}\,{\it t_8}\,{\it t_9}
C>    \,{\it t_{10}}\,\left({\it t_{12}}+1.0\right)}
C>    \over{{{4.1672705928356\,{\it t_{13}}\,{\it t_{10}}^2}
C>    \over{{\it t_3}^{4.666666666666667}\,{\it t_6}^4\,{
C>    \it t_{11}}^2}}+{\it t_{12}}+1.0}}+1.0\right)+1.0\,{\it t_2}
C>    \,{\it t_4}\\\\
C>   g &=& 0\\\\
C>   G &=& 0.01554534543482744\,\log \left({{3.240503646201146\,{
C>    \it t_1}\,{\it t_{16}}\,\sigma_{ss}\,\left({\it t_{18}}
C>    +1.0\right)}\over{{{10.50086388104292\,{\it t_{13}}\,
C>    \sigma_{ss}^2}\over{\rho_s^{4.666666666666667}\,{\it t_{17}}^2}}
C>    +{\it t_{18}}+1.0}}+1.0\right)+1.0\,{\it t_{14}}\,{
C>    \it t_{15}}\\\\
C> \f}
C>
C> Code generated with Maxima 5.34.0 [2,3]
C> driven by autoxc [4,5,6].
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
C> [4] autoxc, revision 27694 2015-10-10
C>
C> [5] rewrap.py revision 27697 2015-10-10
C>
C> [6] call_subroutine.py revision 27339 2015-08-13
C>
      subroutine nwxcm_c_mpbe_d3(param,tol_rho,ipol,nq,wght, rho,
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
            t3 = 1/rhoa
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = exp(-1.6081984221457415d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 1.6202518231005725d+0*gammaaa*t1*t4*t7
            t9 = t1**2
            t10 = gammaaa**2
            t11 = 1/rhoa**4.666666666666667d+0
            t12 = 1/t6**2
            t13 = 2.6252159702607286d+0*t10*t11*t12*t9
            t14 = 1/rhoa**2
            t15 = s1a(D1_RA)
            t16 = 1.6202518231005728d+0*gammaaa*t1*t4*t7
            t17 = t16+1.0d+0
            t18 = t16+t13+1.0d+0
            t19 = 1/t18
            t20 = 1.6202518231005728d+0*gammaaa*t1*t17*t19*t4+1.0d+0
            t21 = 1/t20
            t22 = 1/rhoa**3.3333333333333337d+0
            t23 = -1.8902937936173353d+0*gammaaa*t1*t22*t7
            t24 = 8.040992110728707d+0*t14*s1f-1.6081984221457415d+1
     +        *t15*t3
            t25 = -1.6202518231005728d+0*gammaaa*t1*t12*t24*t4*t5
            t26 = t25+t23
            t27 = 1/t18**2
            t28 = 1/rhoa**5.666666666666667d+0
            t29 = 1/t6**3
            t30 = -5.250431940521457d+0*t10*t11*t24*t29*t5*t9
     +        -6.12550393060837d+0*t10*t12*t28*t9+t25+t23
            t31 = -1.6202518231005728d+0*gammaaa*t1*t17*t27*t30*t4
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t26*t4
     +        -1.8902937936173353d+0*gammaaa*t1*t17*t19*t22
            t32 = 1.3126079851303643d+0*gammaaa*t11*t12*t9
     +        +4.050629557751432d-1*t1*t4*t7
            t33 = 6.563039925651821d-1*gammaaa*t11*t19*t7*t9
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t32*t4
     +        +4.050629557751432d-1*t1*t17*t19*t4
            t34 = 2.6252159702607286d+0*gammaaa*t11*t12*t9
     +        +8.101259115502863d-1*t1*t4*t7
            t35 = 1.3126079851303643d+0*gammaaa*t11*t19*t7*t9
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t34*t4
     +        +8.101259115502863d-1*t1*t17*t19*t4
            t36 = 1/rhoa**3
            t37 = 2.5d-1*s1f*t36
            t38 = -5.0d-1*t14*t15
            t39 = s1a2(D2_RA_RA)
            t40 = 1/rhoa**4.333333333333334d+0
            t41 = -3.7805875872346706d+0*gammaaa*t1*t19*t22*t26
            t42 = 3.7805875872346706d+0*gammaaa*t1*t17*t22*t27*t30
            t43 = -3.2405036462011455d+0*gammaaa*t1*t26*t27*t30*t4
            t44 = t24**2
            t45 = exp(-3.216396844291483d+1*s1f*t3)
            t46 = 3.2405036462011455d+0*gammaaa*t1*t29*t4*t44*t45
            t47 = 3.3305176363733985d+0*gammaaa*t1*t40*t7
            t48 = 3.7805875872346706d+0*gammaaa*t1*t12*t22*t24*t5
            t49 = -1.6202518231005728d+0*gammaaa*t1*t12*t4*t44*t5
            t50 = 1.6081984221457415d+1*t14*t15
            t51 = t50-1.6081984221457415d+1*t3*t39
     +        -1.0721322814304943d+1*s1f*t36
            t52 = -1.6202518231005728d+0*gammaaa*t1*t12*t4*t5*t51
            t53 = t52+t49+t48+t47+t46
            t54 = 1/t6**4
            t55 = 1.5751295821564373d+1*t10*t11*t44*t45*t54*t9
            t56 = 1/rhoa**6.666666666666667d+0
            t57 = 2.450201572243348d+1*t10*t24*t28*t29*t5*t9
            t58 = -5.250431940521457d+0*t10*t11*t29*t44*t5*t9
            t59 = 1.7938975796781656d+1*t10*t12*t56*t9
     +        -5.250431940521457d+0*t10*t11*t29*t5*t51*t9+t58+t57+t55
     +        +t52+t49+t48+t47+t46
            t60 = 1/t18**3
            t61 = t30**2
            t62 = 3.2405036462011455d+0*gammaaa*t1*t17*t4*t60*t61
            t63 = t62-1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*t59
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t4*t53+t43+t42+t41
     +        +3.3305176363733985d+0*gammaaa*t1*t17*t19*t40
            t64 = 1/t20**2
            t65 = t31**2
            t66 = -3.109069086965489d-2*t64*t65
            t67 = log(t20)
            t68 = s1a2(D2_RA_RB)
            t69 = 2.970461675684382d+0*gammaaa*t1*t40*t7
            t70 = -1.6081984221457415d+1*t3*t68+t50
     +        -5.360661407152472d+0*s1f*t36
            t71 = -1.6202518231005728d+0*gammaaa*t1*t12*t4*t5*t70
            t72 = t71+t69+t49+t48+t46
            t73 = -5.250431940521457d+0*t10*t11*t29*t5*t70*t9
     +        +1.6772213143332443d+1*t10*t12*t56*t9+t71+t69+t58+t57
     +        +t55+t49+t48+t46
            t74 = -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*t73
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t4*t72+t62+t43+t42
     +        +t41+2.970461675684382d+0*gammaaa*t1*t17*t19*t40
            t75 = -4.725734484043338d-1*t1*t22*t7
            t76 = -4.050629557751432d-1*t1*t12*t24*t4*t5
            t77 = t76+t75
            t78 = -2.6252159702607286d+0*gammaaa*t11*t24*t29*t5*t9
     +        -3.062751965304185d+0*gammaaa*t12*t28*t9+t76+t75
            t79 = -6.563039925651821d-1*gammaaa*t11*t27*t30*t7*t9
     +        -7.656879913260463d-1*gammaaa*t19*t28*t7*t9
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*t78
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t4*t77
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t30*t32*t4*t60
     +        -1.6202518231005728d+0*gammaaa*t1*t26*t27*t32*t4
     +        -4.050629557751432d-1*t1*t17*t27*t30*t4
     +        +4.050629557751432d-1*t1*t19*t26*t4
     +        +1.8902937936173353d+0*gammaaa*t1*t17*t22*t27*t32
     +        -4.725734484043338d-1*t1*t17*t19*t22
            t80 = (3.109069086965489d-2*t21*t79-3.109069086965489d-2
     +        *t31*t33*t64)*wght
            t81 = -9.451468968086676d-1*t1*t22*t7
            t82 = -8.101259115502863d-1*t1*t12*t24*t4*t5
            t83 = t82+t81
            t84 = -5.250431940521457d+0*gammaaa*t11*t24*t29*t5*t9
     +        -6.12550393060837d+0*gammaaa*t12*t28*t9+t82+t81
            t85 = -1.3126079851303643d+0*gammaaa*t11*t27*t30*t7*t9
     +        -1.5313759826520926d+0*gammaaa*t19*t28*t7*t9
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*t84
     +        +1.6202518231005728d+0*gammaaa*t1*t19*t4*t83
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t30*t34*t4*t60
     +        -1.6202518231005728d+0*gammaaa*t1*t26*t27*t34*t4
     +        -8.101259115502863d-1*t1*t17*t27*t30*t4
     +        +8.101259115502863d-1*t1*t19*t26*t4
     +        +1.8902937936173353d+0*gammaaa*t1*t17*t22*t27*t34
     +        -9.451468968086676d-1*t1*t17*t19*t22
            t86 = t33**2
            t87 = t32**2
            t88 = -8.101259115502863d-1*t1*t17*t27*t32*t4
            t89 = t1**3
            t90 = 1/rhoa**7.0d+0
            t91 = -1.3126079851303643d+0*gammaaa*t11*t27*t32*t7*t9
            t92 = t91-5.316888702309608d-1*gammaaa*t12*t17*t27*t89*t90
     +        +3.2815199628259106d-1*t11*t19*t7*t9+t88
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t4*t60*t87
            t93 = 3.109069086965489d-2*t21*t92*wght
     +        -3.109069086965489d-2*t64*t86*wght
            t94 = t91-1.0633777404619216d+0*gammaaa*t12*t17*t27*t89
     +        *t90-6.563039925651821d-1*gammaaa*t11*t27*t34*t7*t9
     +        +6.563039925651821d-1*t11*t19*t7*t9+t88
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t32*t34*t4*t60
     +        -4.050629557751432d-1*t1*t17*t27*t34*t4
            t95 = t35**2
            t96 = t34**2
            t97 = 3.2405036462011455d+0*gammaaa*t1*t17*t4*t60*t96
     +        -2.1267554809238431d+0*gammaaa*t12*t17*t27*t89*t90
     +        -2.6252159702607286d+0*gammaaa*t11*t27*t34*t7*t9
     +        +1.3126079851303643d+0*t11*t19*t7*t9
     +        -1.6202518231005728d+0*t1*t17*t27*t34*t4
            t98 = 1/rhoa**4
            t99 = -3.75d-1*s1f*t98
            t100 = 7.5d-1*t15*t36
            t101 = s1a3(D3_RA_RA_RA)
            t102 = 1/rhoa**5.333333333333334d+0
            t103 = 1.134176276170401d+1*gammaaa*t1*t22*t26*t27*t30
            t104 = t24**3
            t105 = exp(-4.824595266437224d+1*s1f*t3)
            t106 = -9.721510938603435d+0*gammaaa*t1*t104*t105*t4*t54
            t107 = -1.134176276170401d+1*gammaaa*t1*t22*t29*t44*t45
            t108 = 1.6081984221457415d+1*t14*s1f-3.216396844291483d+1
     +        *t15*t3
            t109 = 3.2405036462011455d+0*gammaaa*t1*t108*t29*t4*t44*t45
            t110 = 3.2405036462011455d+0*gammaaa*t1*t104*t29*t4*t45
            t111 = 9.721510938603435d+0*gammaaa*t1*t24*t29*t4*t45*t51
            t112 = -7.996242793635233d+0*gammaaa*t1*t102*t7
            t113 = -9.991552909120195d+0*gammaaa*t1*t12*t24*t40*t5
            t114 = 5.670881380852005d+0*gammaaa*t1*t12*t22*t44*t5
            t115 = 5.670881380852005d+0*gammaaa*t1*t12*t22*t5*t51
            t116 = -1.6202518231005728d+0*gammaaa*t1*t104*t12*t4*t5
            t117 = -4.8607554693017174d+0*gammaaa*t1*t12*t24*t4*t5*t51
            t118 = 2.412297633218612d+1*s1f*t98+2.412297633218612d+1
     +        *t14*t39-3.216396844291482d+1*t15*t36
     +        -1.6081984221457415d+1*t101*t3
            t119 = -1.6202518231005728d+0*gammaaa*t1*t118*t12*t4*t5
            t120 = -6.300518328625749d+1*t10*t104*t105*t11*t9/t6**5
            t121 = -1.1025907075095065d+2*t10*t28*t44*t45*t54*t9
            t122 = 1.5751295821564373d+1*t10*t108*t11*t44*t45*t54*t9
            t123 = 1.5751295821564373d+1*t10*t104*t11*t45*t54*t9
            t124 = 1/rhoa**7.666666666666667d+0
            t125 = 3.675302358365022d+1*t10*t28*t29*t44*t5*t9
            t126 = -5.250431940521457d+0*t10*t104*t11*t29*t5*t9
            t127 = -1.134176276170401d+1*gammaaa*t1*t17*t22*t60*t61
            t128 = 1/t18**4
            t129 = -9.721510938603435d+0*gammaaa*t1*t128*t17*t30**3*t4
            t130 = 1/t20**3
            t131 = 6.218138173930978d-2*t130*t31**3
            t132 = s1a3(D3_RA_RA_RB)
            t133 = 6.481007292402291d+0*gammaaa*t1*t24*t29*t4*t45*t70
            t134 = 3.2405036462011455d+0*gammaaa*t1*t24*t29*t4*t45*t51
            t135 = -6.436000297316163d+0*gammaaa*t1*t102*t7
            t136 = -9.271440987742162d+0*gammaaa*t1*t12*t24*t40*t5
            t137 = 3.7805875872346706d+0*gammaaa*t1*t12*t22*t5*t70
            t138 = 1.8902937936173353d+0*gammaaa*t1*t12*t22*t5*t51
            t139 = -3.2405036462011455d+0*gammaaa*t1*t12*t24*t4*t5*t70
            t140 = -1.6202518231005728d+0*gammaaa*t1*t12*t24*t4*t5*t51
            t141 = 8.040992110728705d+0*s1f*t98+1.6081984221457415d+1
     +        *t14*t68+8.040992110728707d+0*t14*t39
     +        -2.1442645628609888d+1*t15*t36-1.6081984221457415d+1
     +        *t132*t3
            t142 = -1.6202518231005728d+0*gammaaa*t1*t12*t141*t4*t5
            t143 = -6.218138173930978d-2*t31*t79
            t144 = 1/rhoa**6.666666666666667d+0
            t145 = 1.276146652210076d+0*gammaaa*t144*t19*t7*t9
            t146 = -3.7805875872346706d+0*gammaaa*t1*t19*t22*t77
            t147 = 3.7805875872346706d+0*gammaaa*t1*t17*t22*t27*t78
            t148 = 3.7805875872346706d+0*gammaaa*t1*t22*t26*t27*t32
            t149 = -9.451468968086676d-1*t1*t19*t22*t26
            t150 = -3.2405036462011455d+0*gammaaa*t1*t26*t27*t4*t78
            t151 = -7.561175174469341d+0*gammaaa*t1*t17*t22*t30*t32*t60
            t152 = 9.451468968086676d-1*t1*t17*t22*t27*t30
            t153 = 1.5313759826520926d+0*gammaaa*t27*t28*t30*t7*t9
            t154 = -3.2405036462011455d+0*gammaaa*t1*t27*t30*t4*t77
            t155 = 6.481007292402291d+0*gammaaa*t1*t17*t30*t4*t60*t78
            t156 = 6.481007292402291d+0*gammaaa*t1*t26*t30*t32*t4*t60
            t157 = -8.101259115502863d-1*t1*t26*t27*t30*t4
            t158 = 8.101259115502863d-1*t1*t29*t4*t44*t45
            t159 = 8.326294090933496d-1*t1*t40*t7
            t160 = 9.451468968086676d-1*t1*t12*t22*t24*t5
            t161 = -4.050629557751432d-1*t1*t12*t4*t44*t5
            t162 = -4.050629557751432d-1*t1*t12*t4*t5*t51
            t163 = 7.875647910782186d+0*gammaaa*t11*t44*t45*t54*t9
            t164 = 1.225100786121674d+1*gammaaa*t24*t28*t29*t5*t9
            t165 = -2.6252159702607286d+0*gammaaa*t11*t29*t44*t5*t9
            t166 = -9.721510938603435d+0*gammaaa*t1*t128*t17*t32*t4*t61
            t167 = 8.101259115502863d-1*t1*t17*t4*t60*t61
            t168 = 1.3126079851303643d+0*gammaaa*t11*t60*t61*t7*t9
            t169 = 6.218138173930978d-2*t130*t33*t65
            t170 = (t21*(3.109069086965489d-2*(-1.6202518231005728d+0
     +        *gammaaa*t1*t17*t27*t4*(8.969487898390828d+0*gammaaa*t12
     +        *t56*t9-2.6252159702607286d+0*gammaaa*t11*t29*t5*t51*t9
     +        +t165+t164+t163+t162+t161+t160+t159+t158)
     +        -6.563039925651821d-1*gammaaa*t11*t27*t59*t7*t9
     +        +7.292266584057577d-2*gammaaa*t19*t56*t7*t9
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t32*t4*t59*t60
     +        -4.050629557751432d-1*t1*t17*t27*t4*t59
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t32*t4*t53
     +        +4.050629557751432d-1*t1*t19*t4*t53
     +        -3.3305176363733985d+0*gammaaa*t1*t17*t27*t32*t40
     +        +8.326294090933496d-1*t1*t17*t19*t40
     +        +1.6202518231005728d+0*gammaaa*t1*(t162+t161+t160+t159
     +        +t158)*t19*t4+t168+t167+t166+t157+t156+t155+t154+t153
     +        +t152+t151+t150+t149+t148+t147+t146+t145)
     +        -5.181781811609147d-3*t14*t33)+(t143
     +        -3.109069086965489d-2*t33*t63)*t64+t169)*wght
            t171 = -6.218138173930978d-2*t31*t85
            t172 = 2.552293304420152d+0*gammaaa*t144*t19*t7*t9
            t173 = -3.7805875872346706d+0*gammaaa*t1*t19*t22*t83
            t174 = 3.7805875872346706d+0*gammaaa*t1*t17*t22*t27*t84
            t175 = 3.7805875872346706d+0*gammaaa*t1*t22*t26*t27*t34
            t176 = -1.8902937936173353d+0*t1*t19*t22*t26
            t177 = -3.2405036462011455d+0*gammaaa*t1*t26*t27*t4*t84
            t178 = -7.561175174469341d+0*gammaaa*t1*t17*t22*t30*t34*t60
            t179 = 1.8902937936173353d+0*t1*t17*t22*t27*t30
            t180 = 3.062751965304185d+0*gammaaa*t27*t28*t30*t7*t9
            t181 = -3.2405036462011455d+0*gammaaa*t1*t27*t30*t4*t83
            t182 = 6.481007292402291d+0*gammaaa*t1*t17*t30*t4*t60*t84
            t183 = 6.481007292402291d+0*gammaaa*t1*t26*t30*t34*t4*t60
            t184 = -1.6202518231005728d+0*t1*t26*t27*t30*t4
            t185 = 1.6202518231005728d+0*t1*t29*t4*t44*t45
            t186 = 1.6652588181866992d+0*t1*t40*t7
            t187 = 1.8902937936173353d+0*t1*t12*t22*t24*t5
            t188 = -8.101259115502863d-1*t1*t12*t4*t44*t5
            t189 = -8.101259115502863d-1*t1*t12*t4*t5*t51
            t190 = 1.5751295821564373d+1*gammaaa*t11*t44*t45*t54*t9
            t191 = 2.450201572243348d+1*gammaaa*t24*t28*t29*t5*t9
            t192 = -5.250431940521457d+0*gammaaa*t11*t29*t44*t5*t9
            t193 = -9.721510938603435d+0*gammaaa*t1*t128*t17*t34*t4*t61
            t194 = 1.6202518231005728d+0*t1*t17*t4*t60*t61
            t195 = 2.6252159702607286d+0*gammaaa*t11*t60*t61*t7*t9
            t196 = 6.218138173930978d-2*t130*t35*t65
            t197 = 7.426154189210955d-1*t1*t40*t7
            t198 = -4.050629557751432d-1*t1*t12*t4*t5*t70
            t199 = 1.485230837842191d+0*t1*t40*t7
            t200 = -8.101259115502863d-1*t1*t12*t4*t5*t70
            t201 = 9.451468968086676d-1*t1*t17*t22*t27*t32
            t202 = 1/rhoa**8.0d+0
            t203 = 1.5313759826520926d+0*gammaaa*t27*t28*t32*t7*t9
            t204 = 8.101259115502863d-1*t1*t19*t4*t77
            t205 = -8.101259115502863d-1*t1*t17*t27*t4*t78
            t206 = -1.3126079851303643d+0*gammaaa*t11*t27*t7*t78*t9
            t207 = -8.101259115502863d-1*t1*t26*t27*t32*t4
            t208 = 1.6202518231005728d+0*t1*t17*t30*t32*t4*t60
            t209 = 2.6252159702607286d+0*gammaaa*t11*t30*t32*t60*t7*t9
            t210 = (t64*(-3.109069086965489d-2*t31*t92
     +        -6.218138173930978d-2*t33*t79)+3.109069086965489d-2*t21*
     +        (1.0633777404619216d+0*gammaaa*t12*t17*t30*t60*t89*t90
     +        -5.316888702309608d-1*gammaaa*t12*t26*t27*t89*t90
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -6.563039925651821d-1*t11*t24*t29*t5*t9
     +        -7.656879913260463d-1*t12*t28*t9)-3.2815199628259106d-1
     +        *t11*t27*t30*t7*t9-3.8284399566302313d-1*t19*t28*t7*t9
     +        +6.20303681936121d-1*gammaaa*t12*t17*t202*t27*t89
     +        +3.2405036462011455d+0*gammaaa*t1*t26*t4*t60*t87
     +        -3.7805875872346706d+0*gammaaa*t1*t17*t22*t60*t87
     +        -9.721510938603435d+0*gammaaa*t1*t128*t17*t30*t4*t87
     +        +6.481007292402291d+0*gammaaa*t1*t17*t32*t4*t60*t78
     +        -3.2405036462011455d+0*gammaaa*t1*t27*t32*t4*t77+t209
     +        +t208+t207+t206+t205+t204+t203+t201)
     +        +6.218138173930978d-2*t130*t31*t86)*wght
            t211 = (t64*(-3.109069086965489d-2*t31*t94
     +        -3.109069086965489d-2*t33*t85-3.109069086965489d-2*t35
     +        *t79)+3.109069086965489d-2*t21*(2.1267554809238431d+0
     +        *gammaaa*t12*t17*t30*t60*t89*t90-1.0633777404619216d+0
     +        *gammaaa*t12*t26*t27*t89*t90-1.6202518231005728d+0
     +        *gammaaa*t1*t17*t27*t4*(-1.3126079851303643d+0*t11*t24
     +        *t29*t5*t9-1.5313759826520926d+0*t12*t28*t9)
     +        -6.563039925651821d-1*gammaaa*t11*t27*t7*t84*t9
     +        +1.3126079851303643d+0*gammaaa*t11*t30*t34*t60*t7*t9
     +        +7.656879913260463d-1*gammaaa*t27*t28*t34*t7*t9
     +        -6.563039925651821d-1*t11*t27*t30*t7*t9
     +        -7.656879913260463d-1*t19*t28*t7*t9
     +        +1.2406073638722419d+0*gammaaa*t12*t17*t202*t27*t89
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t32*t4*t60*t84
     +        -4.050629557751432d-1*t1*t17*t27*t4*t84
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t32*t4*t83
     +        +4.050629557751432d-1*t1*t19*t4*t83
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t34*t4*t60*t78
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t34*t4*t77
     +        +3.2405036462011455d+0*gammaaa*t1*t26*t32*t34*t4*t60
     +        +8.101259115502863d-1*t1*t17*t30*t34*t4*t60
     +        -3.7805875872346706d+0*gammaaa*t1*t17*t22*t32*t34*t60
     +        -9.721510938603435d+0*gammaaa*t1*t128*t17*t30*t32*t34*t4
     +        -4.050629557751432d-1*t1*t26*t27*t34*t4
     +        +4.725734484043338d-1*t1*t17*t22*t27*t34+t209+t208+t207
     +        +t206+t205+t204+t203+t201)+6.218138173930978d-2*t130*t31
     +        *t33*t35)*wght
            t212 = t1**4
            t213 = 1/rhoa**9.333333333333333d+0
            t214 = -9.327207260896466d-2*t33*t64*t92*wght
     +        +3.109069086965489d-2*t21*(3.1901332213857647d+0*gammaaa
     +        *t12*t17*t32*t60*t89*t90-3.987666526732206d-1*t12*t17
     +        *t27*t89*t90+3.937823955391093d+0*gammaaa*t11*t60*t7*t87
     +        *t9-9.844559888477733d-1*t11*t27*t32*t7*t9
     +        +2.4303777346508587d+0*t1*t17*t4*t60*t87
     +        -9.721510938603435d+0*gammaaa*t1*t128*t17*t32**3*t4
     +        -6.461023959854982d-1*gammaaa*t212*t213*t27*t29)*wght
     +        +6.218138173930978d-2*t130*t33**3*wght
            t215 = -9.721510938603435d+0*gammaaa*t1*t128*t17*t34*t4*t87
            t216 = 1.6202518231005728d+0*t1*t17*t4*t60*t87
            t217 = 1.6202518231005728d+0*t1*t17*t32*t34*t4*t60
            t218 = 4.2535109618476863d+0*gammaaa*t12*t17*t32*t60*t89*t90
            t219 = 1.0633777404619216d+0*gammaaa*t12*t17*t34*t60*t89*t90
            t220 = -7.975333053464412d-1*t12*t17*t27*t89*t90
            t221 = 2.6252159702607286d+0*gammaaa*t11*t60*t7*t87*t9
            t222 = 2.6252159702607286d+0*gammaaa*t11*t32*t34*t60*t7*t9
            t223 = -1.3126079851303643d+0*t11*t27*t32*t7*t9
            t224 = -3.2815199628259106d-1*t11*t27*t34*t7*t9
            t225 = 6.218138173930978d-2*t130*t35*t86*wght
            t226 = t64*(-6.218138173930978d-2*t33*t94*wght
     +        -3.109069086965489d-2*t35*t92*wght)
            fnc(iq) = (3.109069086965489d-2*log(1.6202518231005725d+0
     +        *gammaaa*t1*t4*(t8+1.0d+0)/(t8+t13+1.0d+0)+1.0d+0)
     +        +5.0d-1*s1f*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (3.109069086965489d-2*t21*t31+5.0d-1*t15
     +        *t3-2.5d-1*t14*s1f)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 3.109069086965489d-2*t21*t33*wght+Cmat
     +        (iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t21*t35*wght+Cmat
     +        (iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-5.181781811609147d-3*t14*t67+t66
     +        +3.109069086965489d-2*t21*t63+5.0d-1*t3*t39+t38+t37)
     +        *wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = (3.109069086965489d-2*t21*t74+5.0d-1
     +        *t3*t68+5.181781811609147d-3*t14*t67+t66+t38+t37)*wght
     +        +Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t80+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = (3.109069086965489d-2*t21*t85
     +        -3.109069086965489d-2*t31*t35*t64)*wght+Cmat2(iq,
     +        D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t80+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t93+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = 3.109069086965489d-2*t21*t94*wght
     +        -3.109069086965489d-2*t33*t35*t64*wght+Cmat2(iq,
     +        D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t93+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 3.109069086965489d-2*t21*t97*wght
     +        -3.109069086965489d-2*t64*t95*wght+Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (t99+t21*(3.109069086965489d-2*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -1.0763385478068993d+2*t10*t24*t29*t5*t56*t9
     +        +4.725388746469312d+1*t10*t11*t24*t45*t51*t54*t9
     +        +3.675302358365022d+1*t10*t28*t29*t5*t51*t9
     +        -1.5751295821564373d+1*t10*t11*t24*t29*t5*t51*t9
     +        -5.250431940521457d+0*t10*t11*t118*t29*t5*t9
     +        -6.368579483410288d+1*t10*t12*t124*t9+t126+t125+t123
     +        +t122+t121+t120+t119+t117+t116+t115+t114+t113+t112+t111
     +        +t110+t109+t107+t106)+9.721510938603435d+0*gammaaa*t1
     +        *t26*t4*t60*t61+9.721510938603435d+0*gammaaa*t1*t17*t30
     +        *t4*t59*t60-4.8607554693017174d+0*gammaaa*t1*t26*t27*t4
     +        *t59+5.670881380852005d+0*gammaaa*t1*t17*t22*t27*t59
     +        -4.8607554693017174d+0*gammaaa*t1*t27*t30*t4*t53
     +        -5.670881380852005d+0*gammaaa*t1*t19*t22*t53
     +        -9.991552909120195d+0*gammaaa*t1*t17*t27*t30*t40
     +        +9.991552909120195d+0*gammaaa*t1*t19*t26*t40
     +        +1.6202518231005728d+0*gammaaa*t1*(t119+t117+t116+t115
     +        +t114+t113+t112+t111+t110+t109+t107+t106)*t19*t4
     +        -7.996242793635233d+0*gammaaa*t1*t102*t17*t19+t129+t127
     +        +t103)-1.5545345434827443d-2*t14*t31)
     +        +1.5545345434827446d-2*t36*t67-9.327207260896466d-2*t31
     +        *t63*t64-7.5d-1*t14*t39+5.0d-1*t101*t3+t131+t100)*wght
     +        +Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = (t99+t21*(3.109069086965489d-2*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*
     +        (3.1502591643128747d+1*t10*t11*t24*t45*t54*t70*t9
     +        +2.450201572243348d+1*t10*t28*t29*t5*t70*t9
     +        -1.0500863881042914d+1*t10*t11*t24*t29*t5*t70*t9
     +        -1.0296680416689308d+2*t10*t24*t29*t5*t56*t9
     +        +1.5751295821564373d+1*t10*t11*t24*t45*t51*t54*t9
     +        +1.225100786121674d+1*t10*t28*t29*t5*t51*t9
     +        -5.250431940521457d+0*t10*t11*t24*t29*t5*t51*t9
     +        -5.250431940521457d+0*t10*t11*t141*t29*t5*t9
     +        -5.590737714444147d+1*t10*t12*t124*t9+t142+t140+t139
     +        +t138+t137+t136+t135+t134+t133+t126+t125+t123+t122+t121
     +        +t120+t116+t114+t110+t109+t107+t106)
     +        +6.481007292402291d+0*gammaaa*t1*t17*t30*t4*t60*t73
     +        -3.2405036462011455d+0*gammaaa*t1*t26*t27*t4*t73
     +        +3.7805875872346706d+0*gammaaa*t1*t17*t22*t27*t73
     +        -3.2405036462011455d+0*gammaaa*t1*t27*t30*t4*t72
     +        -3.7805875872346706d+0*gammaaa*t1*t19*t22*t72
     +        +9.721510938603437d+0*gammaaa*t1*t26*t4*t60*t61
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t30*t4*t59*t60
     +        -1.6202518231005728d+0*gammaaa*t1*t26*t27*t4*t59
     +        +1.8902937936173353d+0*gammaaa*t1*t17*t22*t27*t59
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t30*t4*t53
     +        -1.8902937936173353d+0*gammaaa*t1*t19*t22*t53
     +        -9.271440987742162d+0*gammaaa*t1*t17*t27*t30*t40
     +        +9.271440987742162d+0*gammaaa*t1*t19*t26*t40
     +        +1.6202518231005728d+0*gammaaa*t1*(t142+t140+t139+t138
     +        +t137+t136+t135+t134+t133+t116+t114+t110+t109+t107+t106)
     +        *t19*t4-6.436000297316163d+0*gammaaa*t1*t102*t17*t19
     +        +t129+t127+t103)+5.181781811609149d-3*t14*t31)+t64*(
     +        -6.218138173930978d-2*t31*t74-3.109069086965489d-2*t31
     +        *t63)+2.5d-1*t14*(-2.0d+0*t68-1.0d+0*t39)
     +        -5.18178181160915d-3*t36*t67+5.0d-1*t132*t3+t131+t100)
     +        *wght+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t170+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = (t21*(3.109069086965489d-2*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*
     +        (1.7938975796781656d+1*gammaaa*t12*t56*t9
     +        -5.250431940521457d+0*gammaaa*t11*t29*t5*t51*t9+t192
     +        +t191+t190+t189+t188+t187+t186+t185)
     +        -1.3126079851303643d+0*gammaaa*t11*t27*t59*t7*t9
     +        +1.4584533168115155d-1*gammaaa*t19*t56*t7*t9
     +        +3.2405036462011455d+0*gammaaa*t1*t17*t34*t4*t59*t60
     +        -8.101259115502863d-1*t1*t17*t27*t4*t59
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t34*t4*t53
     +        +8.101259115502863d-1*t1*t19*t4*t53
     +        -3.3305176363733985d+0*gammaaa*t1*t17*t27*t34*t40
     +        +1.6652588181866992d+0*t1*t17*t19*t40
     +        +1.6202518231005728d+0*gammaaa*t1*(t189+t188+t187+t186
     +        +t185)*t19*t4+t195+t194+t193+t184+t183+t182+t181+t180
     +        +t179+t178+t177+t176+t175+t174+t173+t172)
     +        -5.181781811609147d-3*t14*t35)+(t171
     +        -3.109069086965489d-2*t35*t63)*t64+t196)*wght+Cmat3(iq,
     +        D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t170+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = (t21*(3.109069086965489d-2*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -2.6252159702607286d+0*gammaaa*t11*t29*t5*t70*t9
     +        +8.386106571666222d+0*gammaaa*t12*t56*t9+t198+t197+t165
     +        +t164+t163+t161+t160+t158)-6.563039925651821d-1*gammaaa
     +        *t11*t27*t7*t73*t9-7.292266584057578d-2*gammaaa*t19*t56
     +        *t7*t9+3.2405036462011455d+0*gammaaa*t1*t17*t32*t4*t60
     +        *t73-4.050629557751432d-1*t1*t17*t27*t4*t73
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t32*t4*t72
     +        +4.050629557751432d-1*t1*t19*t4*t72-2.970461675684382d+0
     +        *gammaaa*t1*t17*t27*t32*t40+7.426154189210955d-1*t1*t17
     +        *t19*t40+1.6202518231005728d+0*gammaaa*t1*t19*(t198+t197
     +        +t161+t160+t158)*t4+t168+t167+t166+t157+t156+t155+t154
     +        +t153+t152+t151+t150+t149+t148+t147+t146+t145)
     +        +5.181781811609147d-3*t14*t33)+t64*(t143
     +        -3.109069086965489d-2*t33*t74)+t169)*wght+Cmat3(iq,
     +        D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = (t21*(3.109069086965489d-2*(
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -5.250431940521457d+0*gammaaa*t11*t29*t5*t70*t9
     +        +1.6772213143332443d+1*gammaaa*t12*t56*t9+t200+t199+t192
     +        +t191+t190+t188+t187+t185)-1.3126079851303643d+0*gammaaa
     +        *t11*t27*t7*t73*t9-1.4584533168115157d-1*gammaaa*t19*t56
     +        *t7*t9+3.2405036462011455d+0*gammaaa*t1*t17*t34*t4*t60
     +        *t73-8.101259115502863d-1*t1*t17*t27*t4*t73
     +        -1.6202518231005728d+0*gammaaa*t1*t27*t34*t4*t72
     +        +8.101259115502863d-1*t1*t19*t4*t72-2.970461675684382d+0
     +        *gammaaa*t1*t17*t27*t34*t40+1.485230837842191d+0*t1*t17
     +        *t19*t40+1.6202518231005728d+0*gammaaa*t1*t19*(t200+t199
     +        +t188+t187+t185)*t4+t195+t194+t193+t184+t183+t182+t181
     +        +t180+t179+t178+t177+t176+t175+t174+t173+t172)
     +        +5.181781811609147d-3*t14*t35)+t64*(t171
     +        -3.109069086965489d-2*t35*t74)+t196)*wght+Cmat3(iq,
     +        D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = t210+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = t211+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t210+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = (t64*(-3.109069086965489d-2*t31
     +        *t97-6.218138173930978d-2*t35*t85)+3.109069086965489d-2
     +        *t21*(3.2405036462011455d+0*gammaaa*t1*t26*t4*t60*t96
     +        -3.7805875872346706d+0*gammaaa*t1*t17*t22*t60*t96
     +        -9.721510938603435d+0*gammaaa*t1*t128*t17*t30*t4*t96
     +        +4.2535109618476863d+0*gammaaa*t12*t17*t30*t60*t89*t90
     +        -2.1267554809238431d+0*gammaaa*t12*t26*t27*t89*t90
     +        -1.6202518231005728d+0*gammaaa*t1*t17*t27*t4*(
     +        -2.6252159702607286d+0*t11*t24*t29*t5*t9
     +        -3.062751965304185d+0*t12*t28*t9)-2.6252159702607286d+0
     +        *gammaaa*t11*t27*t7*t84*t9+5.250431940521457d+0*gammaaa
     +        *t11*t30*t34*t60*t7*t9+3.062751965304185d+0*gammaaa*t27
     +        *t28*t34*t7*t9-1.3126079851303643d+0*t11*t27*t30*t7*t9
     +        -1.5313759826520926d+0*t19*t28*t7*t9
     +        +2.4812147277444838d+0*gammaaa*t12*t17*t202*t27*t89
     +        +6.481007292402291d+0*gammaaa*t1*t17*t34*t4*t60*t84
     +        -1.6202518231005728d+0*t1*t17*t27*t4*t84
     +        -3.2405036462011455d+0*gammaaa*t1*t27*t34*t4*t83
     +        +1.6202518231005728d+0*t1*t19*t4*t83
     +        +3.2405036462011455d+0*t1*t17*t30*t34*t4*t60
     +        -1.6202518231005728d+0*t1*t26*t27*t34*t4
     +        +1.8902937936173353d+0*t1*t17*t22*t27*t34)
     +        +6.218138173930978d-2*t130*t31*t95)*wght+Cmat3(iq,
     +        D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = t211+Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t210+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t214+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = 3.109069086965489d-2*t21*(
     +        -1.2922047919709964d+0*gammaaa*t212*t213*t27*t29+t224
     +        +t223+t222+t221+t220+t219+t218+t217+t216+t215)*wght+t226
     +        +t225+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t214+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = t64*(-3.109069086965489d-2*t33
     +        *t97*wght-6.218138173930978d-2*t35*t94*wght)
     +        +3.109069086965489d-2*t21*(1.3126079851303643d+0*gammaaa
     +        *t11*t60*t7*t9*t96+8.101259115502863d-1*t1*t17*t4*t60
     +        *t96-9.721510938603435d+0*gammaaa*t1*t128*t17*t32*t4*t96
     +        +4.2535109618476863d+0*gammaaa*t12*t17*t34*t60*t89*t90
     +        -1.5950666106928824d+0*t12*t17*t27*t89*t90
     +        +5.250431940521457d+0*gammaaa*t11*t32*t34*t60*t7*t9
     +        -1.3126079851303643d+0*t11*t27*t34*t7*t9
     +        +3.2405036462011455d+0*t1*t17*t32*t34*t4*t60
     +        -2.584409583941992d+0*gammaaa*t212*t213*t27*t29+t223
     +        +t218)*wght+6.218138173930978d-2*t130*t33*t95*wght+Cmat3
     +        (iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = 3.109069086965489d-2*t21*(
     +        -1.292204791970996d+0*gammaaa*t212*t213*t27*t29+t224
     +        +t223+t222+t221+t220+t219+t218+t217+t216+t215)*wght+t226
     +        +t225+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = -9.327207260896466d-2*t35*t64
     +        *t97*wght+3.109069086965489d-2*t21*(7.875647910782186d+0
     +        *gammaaa*t11*t60*t7*t9*t96+4.8607554693017174d+0*t1*t17
     +        *t4*t60*t96+1.2760532885543058d+1*gammaaa*t12*t17*t34
     +        *t60*t89*t90-3.1901332213857647d+0*t12*t17*t27*t89*t90
     +        -3.937823955391093d+0*t11*t27*t34*t7*t9
     +        -9.721510938603435d+0*gammaaa*t1*t128*t17*t34**3*t4
     +        -5.168819167883986d+0*gammaaa*t212*t213*t27*t29)*wght
     +        +6.218138173930978d-2*t130*t35**3*wght+Cmat3(iq,
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
            t4 = 1/t3
            t5 = rhoa-rhob
            t6 = 1.0d+0-t4*t5
            t7 = t4*t5+1.0d+0
            t8 = 5.0d-1*t7**6.666666666666666d-1+5.0d-1*t6*
     +        *6.666666666666666d-1
            t9 = t8**3.0d+0
            t10 = gammabb+2.0d+0*gammaab+gammaaa
            t11 = 1/t3**2.3333333333333334d+0
            t12 = 1/t8**2
            t13 = 1/t9
            t14 = exp(-3.216396844291483d+1*t13*s1f*t4)
            t15 = t14-1.0d+0
            t16 = 1/t15
            t17 = 2.0413893780549558d+0*t1*t10*t11*t12*t16
            t18 = t17+1.0d+0
            t19 = t1**2
            t20 = t10**2
            t21 = 1/t3**4.666666666666667d+0
            t22 = 1/t8**4
            t23 = 1/t15**2
            t24 = 4.1672705928355996d+0*t19*t20*t21*t22*t23+t17+1.0d+0
            t25 = 1/t24
            t26 = 2.0413893780549558d+0*t1*t10*t11*t12*t18*t25+1.0d+0
            t27 = log(t26)
            t28 = 1/t3**2
            t29 = -1.0d+0*s1f*t28
            t30 = s1a(D1_RA)
            t31 = 1/t26
            t32 = t28*t5
            t33 = -t4
            t34 = t33+t32
            t35 = 1/t6**3.333333333333333d-1
            t36 = -t28*t5
            t37 = t4+t36
            t38 = 1/t7**3.333333333333333d-1
            t39 = 3.333333333333333d-1*t37*t38+3.333333333333333d-1
     +        *t34*t35
            t40 = 1/t8**3
            t41 = 1/t3**3.3333333333333337d+0
            t42 = -4.763241882128231d+0*t1*t10*t12*t18*t25*t41
            t43 = -4.0827787561099116d+0*t1*t10*t11*t16*t39*t40
            t44 = -4.763241882128231d+0*t1*t10*t12*t16*t41
            t45 = 1/t8**4.0d+0
            t46 = 3.216396844291483d+1*t13*s1f*t28
            t47 = t46+9.649190532874448d+1*s1f*t39*t4*t45
     +        -3.216396844291483d+1*t13*t30*t4
            t48 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t47
            t49 = t48+t44+t43
            t50 = 1/t24**2
            t51 = 1/t8**5
            t52 = 1/t3**5.666666666666667d+0
            t53 = -1.9447262766566134d+1*t19*t20*t22*t23*t52
            t54 = 1/t15**3
            t55 = -8.334541185671199d+0*t14*t19*t20*t21*t22*t47*t54
     +        +t53-1.6669082371342397d+1*t19*t20*t21*t23*t39*t51+t48
     +        +t44+t43
            t56 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t55
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*t49+t42
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t39*t40
            t57 = t8**2.0d+0
            t58 = s1a(D1_RB)
            t59 = t4+t32
            t60 = t36+t33
            t61 = 3.333333333333333d-1*t38*t60+3.333333333333333d-1
     +        *t35*t59
            t62 = -4.0827787561099116d+0*t1*t10*t11*t16*t40*t61
            t63 = 9.649190532874448d+1*s1f*t4*t45*t61
     +        -3.216396844291483d+1*t13*t4*t58+t46
            t64 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t63
            t65 = t64+t62+t44
            t66 = t64-8.334541185671199d+0*t14*t19*t20*t21*t22*t54*t63
     +        +t62-1.6669082371342397d+1*t19*t20*t21*t23*t51*t61+t53+t44
            t67 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t66
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t40*t61+t42
            t68 = 8.334541185671199d+0*t10*t19*t21*t22*t23
     +        +2.0413893780549558d+0*t1*t11*t12*t16
            t69 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t68
     +        +4.1672705928355996d+0*t10*t16*t19*t21*t22*t25
     +        +2.0413893780549558d+0*t1*t11*t12*t18*t25
            t70 = 3.109069086965489d-2*t31*t69*t9*wght
            t71 = 1.6669082371342397d+1*t10*t19*t21*t22*t23
     +        +4.0827787561099116d+0*t1*t11*t12*t16
            t72 = -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*t71
     +        +8.334541185671199d+0*t10*t16*t19*t21*t22*t25
     +        +4.0827787561099116d+0*t1*t11*t12*t18*t25
            t73 = 1/t3**3
            t74 = 2.0d+0*s1f*t73
            t75 = s1a2(D2_RA_RA)
            t76 = t39**2
            t77 = t34**2
            t78 = 1/t6**1.3333333333333333d+0
            t79 = -2*t5*t73
            t80 = 2*t28
            t81 = t80+t79
            t82 = t37**2
            t83 = 1/t7**1.3333333333333333d+0
            t84 = 2*t5*t73
            t85 = -2*t28
            t86 = t85+t84
            t87 = 3.333333333333333d-1*t38*t86-1.111111111111111d-1
     +        *t82*t83+3.333333333333333d-1*t35*t81
     +        -1.111111111111111d-1*t77*t78
            t88 = 1/t3**4.333333333333334d+0
            t89 = 1.5877472940427437d+1*t1*t10*t12*t18*t25*t88
            t90 = t47**2
            t91 = exp(-6.432793688582966d+1*t13*s1f*t4)
            t92 = 4.0827787561099116d+0*t1*t10*t11*t12*t54*t90*t91
            t93 = 1.2248336268329735d+1*t1*t10*t11*t16*t22*t76
            t94 = -4.0827787561099116d+0*t1*t10*t11*t16*t40*t87
            t95 = 1.9052967528512924d+1*t1*t10*t16*t39*t40*t41
            t96 = 1.5877472940427437d+1*t1*t10*t12*t16*t88
            t97 = 8.165557512219823d+0*t1*t10*t11*t14*t23*t39*t40*t47
            t98 = 9.526483764256461d+0*t1*t10*t12*t14*t23*t41*t47
            t99 = 1/t8**5.0d+0
            t100 = -6.432793688582966d+1*t13*s1f*t73
            t101 = -3.859676213149779d+2*s1f*t4*t76*t99
     +        +9.649190532874448d+1*s1f*t4*t45*t87
     +        -3.216396844291483d+1*t13*t4*t75+1.9298381065748896d+2
     +        *t30*t39*t4*t45-1.9298381065748896d+2*s1f*t28*t39*t45
     +        +6.432793688582966d+1*t13*t28*t30+t100
            t102 = -2.0413893780549558d+0*t1*t10*t101*t11*t12*t14*t23
            t103 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t90
            t104 = t98+t97+t96+t95+t94+t93+t92+t103+t102
            t105 = 1/t15**4
            t106 = 1/t8**6
            t107 = 1/t3**6.666666666666667d+0
            t108 = 1.1020115567720809d+2*t107*t19*t20*t22*t23
            t109 = t98+t97+t96+t95+t94+t93+t92+2.50036235570136d+1
     +        *t105*t19*t20*t21*t22*t90*t91-8.334541185671199d+0*t14
     +        *t19*t20*t21*t22*t54*t90-1.6669082371342397d+1*t19*t20
     +        *t21*t23*t51*t87+8.334541185671199d+1*t106*t19*t20*t21
     +        *t23*t76+7.778905106626454d+1*t14*t19*t20*t22*t47*t52
     +        *t54+6.667632948536959d+1*t14*t19*t20*t21*t39*t47*t51
     +        *t54-8.334541185671199d+0*t101*t14*t19*t20*t21*t22*t54
     +        +1.5557810213252907d+2*t19*t20*t23*t39*t51*t52+t108+t103
     +        +t102
            t110 = 1/t24**3
            t111 = t55**2
            t112 = t89-4.0827787561099116d+0*t1*t10*t11*t18*t25*t40
     +        *t87+1.2248336268329735d+1*t1*t10*t11*t18*t22*t25*t76
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t49*t50*t55
     +        +9.526483764256461d+0*t1*t10*t12*t18*t41*t50*t55
     +        +8.165557512219823d+0*t1*t10*t11*t18*t39*t40*t50*t55
     +        -2.0413893780549558d+0*t1*t10*t109*t11*t12*t18*t50
     +        -9.526483764256461d+0*t1*t10*t12*t25*t41*t49
     +        -8.165557512219823d+0*t1*t10*t11*t25*t39*t40*t49
     +        +1.9052967528512924d+1*t1*t10*t18*t25*t39*t40*t41
     +        +2.0413893780549558d+0*t1*t10*t104*t11*t12*t25
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t111*t12*t18
            t113 = 3.109069086965489d-2*t112*t9+1.8654414521792934d-1
     +        *t39*t56*t57
            t114 = 1/t26**2
            t115 = t56**2
            t116 = -1.0d+0*t58-1.0d+0*t30
            t117 = s1a2(D2_RA_RB)
            t118 = -1.111111111111111d-1*t37*t60*t83
     +        -1.111111111111111d-1*t34*t59*t78+6.666666666666666d-1
     +        *t38*t5*t73-6.666666666666666d-1*t35*t5*t73
            t119 = 4.0827787561099116d+0*t1*t10*t11*t12*t47*t54*t63*t91
            t120 = 1.2248336268329735d+1*t1*t10*t11*t16*t22*t39*t61
            t121 = -4.0827787561099116d+0*t1*t10*t11*t118*t16*t40
            t122 = 9.526483764256461d+0*t1*t10*t16*t40*t41*t61
            t123 = 9.526483764256461d+0*t1*t10*t16*t39*t40*t41
            t124 = 4.0827787561099116d+0*t1*t10*t11*t14*t23*t40*t47*t61
            t125 = 4.763241882128231d+0*t1*t10*t12*t14*t23*t41*t47
            t126 = -3.859676213149779d+2*s1f*t39*t4*t61*t99
     +        +9.649190532874448d+1*t30*t4*t45*t61
     +        -9.649190532874448d+1*s1f*t28*t45*t61
     +        +9.649190532874448d+1*t39*t4*t45*t58
     +        +3.216396844291483d+1*t13*t28*t58+9.649190532874448d+1
     +        *t118*s1f*t4*t45-9.649190532874448d+1*s1f*t28*t39*t45
     +        -3.216396844291483d+1*t117*t13*t4+3.216396844291483d+1
     +        *t13*t28*t30+t100
            t127 = -2.0413893780549558d+0*t1*t10*t11*t12*t126*t14*t23
            t128 = 4.0827787561099116d+0*t1*t10*t11*t14*t23*t39*t40*t63
            t129 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t47*t63
            t130 = 4.763241882128231d+0*t1*t10*t12*t14*t23*t41*t63
            t131 = t96+t130+t129+t128+t127+t125+t124+t123+t122+t121
     +        +t120+t119
            t132 = t96+2.50036235570136d+1*t105*t19*t20*t21*t22*t47
     +        *t63*t91+3.889452553313227d+1*t14*t19*t20*t22*t52*t54
     +        *t63+3.3338164742684795d+1*t14*t19*t20*t21*t39*t51*t54
     +        *t63-8.334541185671199d+0*t14*t19*t20*t21*t22*t47*t54
     +        *t63+3.3338164742684795d+1*t14*t19*t20*t21*t47*t51*t54
     +        *t61+7.778905106626454d+1*t19*t20*t23*t51*t52*t61
     +        +8.334541185671199d+1*t106*t19*t20*t21*t23*t39*t61
     +        +3.889452553313227d+1*t14*t19*t20*t22*t47*t52*t54
     +        -8.334541185671199d+0*t126*t14*t19*t20*t21*t22*t54
     +        +7.778905106626454d+1*t19*t20*t23*t39*t51*t52
     +        -1.6669082371342397d+1*t118*t19*t20*t21*t23*t51+t130
     +        +t129+t128+t127+t125+t124+t123+t122+t121+t120+t119+t108
            t133 = t89+4.0827787561099116d+0*t1*t10*t11*t110*t12*t18
     +        *t55*t66-2.0413893780549558d+0*t1*t10*t11*t12*t49*t50
     +        *t66+4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t66
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t39*t40*t50*t66
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t50*t55*t65
     +        -4.763241882128231d+0*t1*t10*t12*t25*t41*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t25*t39*t40*t65
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t55*t61
     +        -4.0827787561099116d+0*t1*t10*t11*t25*t40*t49*t61
     +        +9.526483764256461d+0*t1*t10*t18*t25*t40*t41*t61
     +        +1.2248336268329735d+1*t1*t10*t11*t18*t22*t25*t39*t61
     +        +4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t55
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t132*t18*t50
     +        -4.763241882128231d+0*t1*t10*t12*t25*t41*t49
     +        +9.526483764256461d+0*t1*t10*t18*t25*t39*t40*t41
     +        -4.0827787561099116d+0*t1*t10*t11*t118*t18*t25*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t131*t25
            t134 = 3.109069086965489d-2*t133*t9+9.327207260896466d-2
     +        *t39*t57*t67+9.327207260896466d-2*t56*t57*t61
            t135 = s1a2(D2_RB_RB)
            t136 = t61**2
            t137 = t59**2
            t138 = t85+t79
            t139 = t60**2
            t140 = t84+t80
            t141 = -1.111111111111111d-1*t139*t83-1.111111111111111d-1
     +        *t137*t78+3.333333333333333d-1*t140*t38
     +        +3.333333333333333d-1*t138*t35
            t142 = t63**2
            t143 = 4.0827787561099116d+0*t1*t10*t11*t12*t142*t54*t91
            t144 = 1.2248336268329735d+1*t1*t10*t11*t136*t16*t22
            t145 = -4.0827787561099116d+0*t1*t10*t11*t141*t16*t40
            t146 = 1.9052967528512924d+1*t1*t10*t16*t40*t41*t61
            t147 = 8.165557512219823d+0*t1*t10*t11*t14*t23*t40*t61*t63
            t148 = 9.526483764256461d+0*t1*t10*t12*t14*t23*t41*t63
            t149 = -3.859676213149779d+2*t136*s1f*t4*t99
     +        +1.9298381065748896d+2*t4*t45*t58*t61
     +        -1.9298381065748896d+2*s1f*t28*t45*t61
     +        +6.432793688582966d+1*t13*t28*t58+9.649190532874448d+1
     +        *t141*s1f*t4*t45-3.216396844291483d+1*t13*t135*t4+t100
            t150 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t149*t23
            t151 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t142*t23
            t152 = t96+t151+t150+t148+t147+t146+t145+t144+t143
            t153 = t96+2.50036235570136d+1*t105*t142*t19*t20*t21*t22
     +        *t91+6.667632948536959d+1*t14*t19*t20*t21*t51*t54*t61
     +        *t63+7.778905106626454d+1*t14*t19*t20*t22*t52*t54*t63
     +        +1.5557810213252907d+2*t19*t20*t23*t51*t52*t61
     +        -8.334541185671199d+0*t14*t149*t19*t20*t21*t22*t54
     +        -8.334541185671199d+0*t14*t142*t19*t20*t21*t22*t54
     +        -1.6669082371342397d+1*t141*t19*t20*t21*t23*t51
     +        +8.334541185671199d+1*t106*t136*t19*t20*t21*t23+t151
     +        +t150+t148+t147+t146+t145+t144+t143+t108
            t154 = t66**2
            t155 = t89-4.0827787561099116d+0*t1*t10*t11*t12*t50*t65
     +        *t66+8.165557512219823d+0*t1*t10*t11*t18*t40*t50*t61*t66
     +        +9.526483764256461d+0*t1*t10*t12*t18*t41*t50*t66
     +        -8.165557512219823d+0*t1*t10*t11*t25*t40*t61*t65
     +        -9.526483764256461d+0*t1*t10*t12*t25*t41*t65
     +        +1.9052967528512924d+1*t1*t10*t18*t25*t40*t41*t61
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t153*t18*t50
     +        -4.0827787561099116d+0*t1*t10*t11*t141*t18*t25*t40
     +        +1.2248336268329735d+1*t1*t10*t11*t136*t18*t22*t25
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t152*t25
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t154*t18
            t156 = 3.109069086965489d-2*t155*t9+1.8654414521792934d-1
     +        *t57*t61*t67
            t157 = t67**2
            t158 = 4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t68
            t159 = -4.763241882128231d+0*t1*t12*t18*t25*t41
            t160 = -9.723631383283067d+0*t10*t16*t19*t22*t25*t52
            t161 = -4.0827787561099116d+0*t1*t11*t16*t39*t40
            t162 = -4.763241882128231d+0*t1*t12*t16*t41
            t163 = -2.0413893780549558d+0*t1*t11*t12*t14*t23*t47
            t164 = t163+t162+t161
            t165 = -3.889452553313227d+1*t10*t19*t22*t23*t52
            t166 = -1.6669082371342397d+1*t10*t14*t19*t21*t22*t47*t54
     +        -3.3338164742684795d+1*t10*t19*t21*t23*t39*t51+t165+t163
     +        +t162+t161
            t167 = 4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t55
     +        *t68-2.0413893780549558d+0*t1*t10*t11*t12*t49*t50*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t39*t40*t50*t68
     +        -4.1672705928355996d+0*t10*t16*t19*t21*t22*t50*t55
     +        -2.0413893780549558d+0*t1*t11*t12*t18*t50*t55
     +        -8.334541185671199d+0*t10*t16*t19*t21*t25*t39*t51
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t166*t18*t50
     +        +2.0413893780549558d+0*t1*t11*t12*t25*t49
     +        -4.0827787561099116d+0*t1*t11*t18*t25*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t164*t25+t160+t159
     +        +t158
            t168 = 3.109069086965489d-2*t167*t9+9.327207260896466d-2
     +        *t39*t57*t69
            t169 = (t168*t31-3.109069086965489d-2*t114*t56*t69*t9)*wght
            t170 = 4.763241882128231d+0*t1*t10*t12*t18*t41*t50*t71
            t171 = -9.526483764256461d+0*t1*t12*t18*t25*t41
            t172 = -1.9447262766566134d+1*t10*t16*t19*t22*t25*t52
            t173 = -8.165557512219823d+0*t1*t11*t16*t39*t40
            t174 = -9.526483764256461d+0*t1*t12*t16*t41
            t175 = -4.0827787561099116d+0*t1*t11*t12*t14*t23*t47
            t176 = t175+t174+t173
            t177 = -7.778905106626454d+1*t10*t19*t22*t23*t52
            t178 = -3.3338164742684795d+1*t10*t14*t19*t21*t22*t47*t54
     +        -6.667632948536959d+1*t10*t19*t21*t23*t39*t51+t177+t175
     +        +t174+t173
            t179 = 4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t55
     +        *t71-2.0413893780549558d+0*t1*t10*t11*t12*t49*t50*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t39*t40*t50*t71
     +        -8.334541185671199d+0*t10*t16*t19*t21*t22*t50*t55
     +        -4.0827787561099116d+0*t1*t11*t12*t18*t50*t55
     +        -1.6669082371342397d+1*t10*t16*t19*t21*t25*t39*t51
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t178*t18*t50
     +        +4.0827787561099116d+0*t1*t11*t12*t25*t49
     +        -8.165557512219823d+0*t1*t11*t18*t25*t39*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t176*t25+t172+t171
     +        +t170
            t180 = 3.109069086965489d-2*t179*t9+9.327207260896466d-2
     +        *t39*t57*t72
            t181 = -4.0827787561099116d+0*t1*t11*t16*t40*t61
            t182 = -2.0413893780549558d+0*t1*t11*t12*t14*t23*t63
            t183 = t182+t181+t162
            t184 = -1.6669082371342397d+1*t10*t14*t19*t21*t22*t54*t63
     +        -3.3338164742684795d+1*t10*t19*t21*t23*t51*t61+t182+t181
     +        +t165+t162
            t185 = 4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t66
     +        *t68-2.0413893780549558d+0*t1*t10*t11*t12*t50*t65*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t61*t68
     +        -4.1672705928355996d+0*t10*t16*t19*t21*t22*t50*t66
     +        -2.0413893780549558d+0*t1*t11*t12*t18*t50*t66
     +        +2.0413893780549558d+0*t1*t11*t12*t25*t65
     +        -8.334541185671199d+0*t10*t16*t19*t21*t25*t51*t61
     +        -4.0827787561099116d+0*t1*t11*t18*t25*t40*t61
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t184*t50
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t183*t25+t160+t159
     +        +t158
            t186 = 3.109069086965489d-2*t185*t9+9.327207260896466d-2
     +        *t57*t61*t69
            t187 = (t186*t31-3.109069086965489d-2*t114*t67*t69*t9)*wght
            t188 = -8.165557512219823d+0*t1*t11*t16*t40*t61
            t189 = -4.0827787561099116d+0*t1*t11*t12*t14*t23*t63
            t190 = t189+t188+t174
            t191 = -3.3338164742684795d+1*t10*t14*t19*t21*t22*t54*t63
     +        -6.667632948536959d+1*t10*t19*t21*t23*t51*t61+t189+t188
     +        +t177+t174
            t192 = 4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t66
     +        *t71-2.0413893780549558d+0*t1*t10*t11*t12*t50*t65*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t61*t71
     +        -8.334541185671199d+0*t10*t16*t19*t21*t22*t50*t66
     +        -4.0827787561099116d+0*t1*t11*t12*t18*t50*t66
     +        +4.0827787561099116d+0*t1*t11*t12*t25*t65
     +        -1.6669082371342397d+1*t10*t16*t19*t21*t25*t51*t61
     +        -8.165557512219823d+0*t1*t11*t18*t25*t40*t61
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t191*t50
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t190*t25+t172+t171
     +        +t170
            t193 = 3.109069086965489d-2*t192*t9+9.327207260896466d-2
     +        *t57*t61*t72
            t194 = t69**2
            t195 = t68**2
            t196 = -4.0827787561099116d+0*t1*t11*t12*t18*t50*t68
            t197 = t1**3
            t198 = 1/t3**7.0d+0
            t199 = -8.334541185671199d+0*t10*t16*t19*t21*t22*t50*t68
            t200 = -1.7014043847390745d+1*t10*t106*t18*t197*t198*t23
     +        *t50+8.334541185671199d+0*t16*t19*t21*t22*t25+t199+t196
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t195
            t201 = 3.109069086965489d-2*t200*t31*t9*wght
     +        -3.109069086965489d-2*t114*t194*t9*wght
            t202 = 4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t68
     +        *t71-4.1672705928355996d+0*t10*t16*t19*t21*t22*t50*t71
     +        -2.0413893780549558d+0*t1*t11*t12*t18*t50*t71
     +        -3.402808769478149d+1*t10*t106*t18*t197*t198*t23*t50
     +        +1.6669082371342397d+1*t16*t19*t21*t22*t25+t199+t196
            t203 = 3.109069086965489d-2*t202*t31*t9*wght
     +        -3.109069086965489d-2*t114*t69*t72*t9*wght
            t204 = t72**2
            t205 = t71**2
            t206 = -1.6669082371342397d+1*t10*t16*t19*t21*t22*t50*t71
     +        -8.165557512219823d+0*t1*t11*t12*t18*t50*t71
     +        -6.805617538956298d+1*t10*t106*t18*t197*t198*t23*t50
     +        +3.3338164742684795d+1*t16*t19*t21*t22*t25
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t205
            t207 = 1/t3**4
            t208 = -6.0d+0*s1f*t207
            t209 = s1a3(D3_RA_RA_RA)
            t210 = t39**3
            t211 = 1/t6**2.3333333333333334d+0
            t212 = 6*t207*t5
            t213 = -6*t73
            t214 = 1/t7**2.3333333333333334d+0
            t215 = -6*t207*t5
            t216 = 6*t73
            t217 = -3.333333333333333d-1*t37*t83*t86
     +        -3.333333333333333d-1*t34*t78*t81+3.333333333333333d-1*
     +        (t216+t215)*t38+1.4814814814814814d-1*t214*t37**3
     +        +3.333333333333333d-1*(t213+t212)*t35
     +        +1.4814814814814814d-1*t211*t34**3
            t218 = 1/t3**5.333333333333334d+0
            t219 = -6.880238274185223d+1*t1*t10*t12*t18*t218*t25
            t220 = t47**3
            t221 = exp(-9.649190532874448d+1*t13*s1f*t4)
            t222 = -1.2248336268329735d+1*t1*t10*t105*t11*t12*t220*t221
            t223 = 1.2248336268329735d+1*t1*t10*t101*t11*t12*t47*t54*t91
            t224 = -2.449667253665947d+1*t1*t10*t11*t39*t40*t54*t90*t91
            t225 = 6.432793688582966d+1*t13*s1f*t28
            t226 = 1.9298381065748896d+2*s1f*t39*t4*t45
     +        -6.432793688582966d+1*t13*t30*t4+t225
            t227 = 4.0827787561099116d+0*t1*t10*t11*t12*t226*t54*t90*t91
            t228 = -2.8579451292769387d+1*t1*t10*t12*t41*t54*t90*t91
            t229 = 4.0827787561099116d+0*t1*t10*t11*t12*t220*t54*t91
            t230 = -4.899334507331894d+1*t1*t10*t11*t16*t210*t51
            t231 = 3.6745008804989204d+1*t1*t10*t11*t16*t22*t39*t87
            t232 = -8.573835387830815d+1*t1*t10*t16*t22*t41*t76
            t233 = -4.0827787561099116d+0*t1*t10*t11*t16*t217*t40
            t234 = 2.8579451292769387d+1*t1*t10*t16*t40*t41*t87
            t235 = -9.526483764256461d+1*t1*t10*t16*t39*t40*t88
            t236 = -6.880238274185223d+1*t1*t10*t12*t16*t218
            t237 = -3.6745008804989204d+1*t1*t10*t11*t14*t22*t23*t47*t76
            t238 = 1.2248336268329735d+1*t1*t10*t11*t14*t23*t40*t47*t87
            t239 = -5.715890258553877d+1*t1*t10*t14*t23*t39*t40*t41*t47
            t240 = -4.763241882128231d+1*t1*t10*t12*t14*t23*t47*t88
            t241 = 1.2248336268329735d+1*t1*t10*t101*t11*t14*t23*t39*t40
            t242 = -6.124168134164868d+0*t1*t10*t101*t11*t12*t14*t23*t47
            t243 = 1.4289725646384693d+1*t1*t10*t101*t12*t14*t23*t41
            t244 = 1/t8**6.0d+0
            t245 = 1.9298381065748896d+2*t13*s1f*t207
            t246 = -1.1579028639449339d+3*s1f*t39*t4*t87*t99
     +        -1.1579028639449339d+3*t30*t4*t76*t99
     +        +1.1579028639449339d+3*s1f*t28*t76*t99
     +        +2.8947571598623345d+2*t30*t4*t45*t87
     +        -2.8947571598623345d+2*s1f*t28*t45*t87
     +        +2.8947571598623345d+2*t39*t4*t45*t75
     +        +9.649190532874448d+1*t13*t28*t75+5.789514319724669d+2
     +        *s1f*t39*t45*t73-1.9298381065748896d+2*t13*t30*t73
     +        +9.649190532874448d+1*s1f*t217*t4*t45
     +        -5.789514319724669d+2*t28*t30*t39*t45
     +        +1.9298381065748896d+3*s1f*t210*t244*t4
     +        -3.216396844291483d+1*t13*t209*t4+t245
            t247 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t246
            t248 = 1.2248336268329735d+1*t1*t10*t11*t14*t23*t39*t40*t90
            t249 = 1.4289725646384693d+1*t1*t10*t12*t14*t23*t41*t90
            t250 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t220*t23
            t251 = 1/t15**5
            t252 = 1/t8**7
            t253 = -7.346743711813873d+2*t19*t20*t22*t23/t3*
     +        *7.666666666666667d+0
            t254 = 1/t24**4
            t255 = 1/t26**3
            t256 = 2.0d+0*t58
            t257 = -2.0d+0*t117
            t258 = s1a3(D3_RA_RA_RB)
            t259 = -5.715890258553877d+1*t1*t10*t18*t22*t25*t39*t41*t61
            t260 = 1.9052967528512924d+1*t1*t10*t118*t18*t25*t40*t41
            t261 = -1.111111111111111d-1*t60*t83*t86
     +        -4.444444444444444d-1*t37*t5*t73*t83
     +        +1.4814814814814814d-1*t214*t60*t82-1.111111111111111d-1
     +        *t59*t78*t81+4.444444444444444d-1*t34*t5*t73*t78
     +        +1.4814814814814814d-1*t211*t59*t77+3.333333333333333d-1
     +        *t38*(2*t73+t215)+3.333333333333333d-1*t35*(t212-2*t73)
            t262 = 1.9052967528512924d+1*t1*t10*t25*t40*t41*t49*t61
            t263 = -1.9052967528512924d+1*t1*t10*t18*t40*t41*t50*t55*t61
            t264 = -9.526483764256461d+0*t1*t10*t12*t131*t25*t41
            t265 = 9.526483764256461d+0*t1*t10*t12*t132*t18*t41*t50
            t266 = 1.9052967528512924d+1*t1*t10*t25*t39*t40*t41*t65
            t267 = 9.526483764256461d+0*t1*t10*t12*t41*t50*t55*t65
            t268 = -1.9052967528512924d+1*t1*t10*t18*t39*t40*t41*t50*t66
            t269 = 9.526483764256461d+0*t1*t10*t12*t41*t49*t50*t66
            t270 = -1.9052967528512924d+1*t1*t10*t110*t12*t18*t41*t55
     +        *t66
            t271 = -1.2248336268329735d+1*t1*t10*t105*t11*t12*t221*t63
     +        *t90
            t272 = 8.165557512219823d+0*t1*t10*t11*t12*t126*t47*t54*t91
            t273 = 1.9298381065748896d+2*s1f*t4*t45*t61
     +        -6.432793688582966d+1*t13*t4*t58+t225
            t274 = 4.0827787561099116d+0*t1*t10*t11*t12*t273*t54*t90*t91
            t275 = -1.6331115024439646d+1*t1*t10*t11*t39*t40*t47*t54
     +        *t63*t91
            t276 = -1.9052967528512924d+1*t1*t10*t12*t41*t47*t54*t63*t91
            t277 = 4.0827787561099116d+0*t1*t10*t101*t11*t12*t54*t63*t91
            t278 = 4.0827787561099116d+0*t1*t10*t11*t12*t54*t63*t90*t91
            t279 = -8.165557512219823d+0*t1*t10*t11*t40*t54*t61*t90*t91
            t280 = -9.526483764256461d+0*t1*t10*t12*t41*t54*t90*t91
            t281 = -4.899334507331894d+1*t1*t10*t11*t16*t51*t61*t76
            t282 = 1.2248336268329735d+1*t1*t10*t11*t16*t22*t61*t87
            t283 = 2.449667253665947d+1*t1*t10*t11*t118*t16*t22*t39
            t284 = -5.715890258553877d+1*t1*t10*t16*t22*t39*t41*t61
            t285 = -2.8579451292769387d+1*t1*t10*t16*t22*t41*t76
            t286 = 1.9052967528512924d+1*t1*t10*t118*t16*t40*t41
            t287 = -4.0827787561099116d+0*t1*t10*t11*t16*t261*t40
            t288 = 9.526483764256461d+0*t1*t10*t16*t40*t41*t87
            t289 = -3.1754945880854873d+1*t1*t10*t16*t40*t61*t88
            t290 = -6.350989176170975d+1*t1*t10*t16*t39*t40*t88
            t291 = -2.449667253665947d+1*t1*t10*t11*t14*t22*t23*t39
     +        *t47*t61
            t292 = 8.165557512219823d+0*t1*t10*t11*t118*t14*t23*t40*t47
            t293 = -1.9052967528512924d+1*t1*t10*t14*t23*t40*t41*t47*t61
            t294 = -1.9052967528512924d+1*t1*t10*t14*t23*t39*t40*t41*t47
            t295 = -3.1754945880854873d+1*t1*t10*t12*t14*t23*t47*t88
            t296 = 8.165557512219823d+0*t1*t10*t11*t126*t14*t23*t39*t40
            t297 = -4.0827787561099116d+0*t1*t10*t11*t12*t126*t14*t23
     +        *t47
            t298 = 9.526483764256461d+0*t1*t10*t12*t126*t14*t23*t41
            t299 = 4.0827787561099116d+0*t1*t10*t101*t11*t14*t23*t40*t61
            t300 = 4.763241882128231d+0*t1*t10*t101*t12*t14*t23*t41
            t301 = 7.719352426299558d+2*s1f*t28*t39*t61*t99
            t302 = -1.9298381065748896d+2*t118*s1f*t28*t45
            t303 = -1.9298381065748896d+2*t28*t30*t45*t61
            t304 = -1.9298381065748896d+2*t28*t39*t45*t58
            t305 = 6.432793688582966d+1*t117*t13*t28
            t306 = -3.859676213149779d+2*s1f*t4*t61*t87*t99
     +        -3.859676213149779d+2*t4*t58*t76*t99
     +        +3.859676213149779d+2*s1f*t28*t76*t99
     +        -7.719352426299558d+2*t30*t39*t4*t61*t99
     +        -7.719352426299558d+2*t118*s1f*t39*t4*t99
     +        +9.649190532874448d+1*t4*t45*t58*t87
     +        -9.649190532874448d+1*s1f*t28*t45*t87
     +        +1.9298381065748896d+3*s1f*t244*t4*t61*t76
     +        +9.649190532874448d+1*t4*t45*t61*t75
     +        +3.216396844291483d+1*t13*t28*t75+1.9298381065748896d+2
     +        *s1f*t45*t61*t73-6.432793688582966d+1*t13*t58*t73
     +        +3.859676213149779d+2*s1f*t39*t45*t73
     +        -1.2865587377165932d+2*t13*t30*t73+1.9298381065748896d+2
     +        *t117*t39*t4*t45+1.9298381065748896d+2*t118*t30*t4*t45
     +        +9.649190532874448d+1*s1f*t261*t4*t45
     +        -1.9298381065748896d+2*t28*t30*t39*t45
     +        -3.216396844291483d+1*t13*t258*t4+t305+t304+t303+t302
     +        +t301+t245
            t307 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t306
            t308 = -1.2248336268329735d+1*t1*t10*t11*t14*t22*t23*t63*t76
            t309 = 4.0827787561099116d+0*t1*t10*t11*t14*t23*t40*t63*t87
            t310 = -1.9052967528512924d+1*t1*t10*t14*t23*t39*t40*t41*t63
            t311 = 8.165557512219823d+0*t1*t10*t11*t14*t23*t39*t40*t47
     +        *t63
            t312 = 9.526483764256461d+0*t1*t10*t12*t14*t23*t41*t47*t63
            t313 = -2.0413893780549558d+0*t1*t10*t101*t11*t12*t14*t23
     +        *t63
            t314 = -1.5877472940427437d+1*t1*t10*t12*t14*t23*t63*t88
            t315 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t63*t90
            t316 = 4.0827787561099116d+0*t1*t10*t11*t14*t23*t40*t61*t90
            t317 = 4.763241882128231d+0*t1*t10*t12*t14*t23*t41*t90
            t318 = -2.3336715319879361d+2*t105*t19*t20*t22*t47*t52*t63
     +        *t91
            t319 = -7.778905106626453d+2*t106*t19*t20*t23*t39*t52*t61
            t320 = 1.5557810213252907d+2*t118*t19*t20*t23*t51*t52
            t321 = -3.1115620426505813d+2*t14*t19*t20*t47*t51*t52*t54
     +        *t61
            t322 = 7.778905106626454d+1*t126*t14*t19*t20*t22*t52*t54
            t323 = -3.1115620426505813d+2*t14*t19*t20*t39*t51*t52*t54
     +        *t63
            t324 = 7.778905106626454d+1*t14*t19*t20*t22*t47*t52*t54*t63
            t325 = s1a3(D3_RA_RB_RB)
            t326 = -4.444444444444444d-1*t5*t60*t73*t83
     +        -1.111111111111111d-1*t140*t37*t83+4.444444444444444d-1
     +        *t5*t59*t73*t78-1.111111111111111d-1*t138*t34*t78
     +        -6.666666666666666d-1*t38*t73+6.666666666666666d-1*t35
     +        *t73-2.0d+0*t207*t38*t5+2.0d+0*t207*t35*t5
     +        +1.4814814814814814d-1*t139*t214*t37
     +        +1.4814814814814814d-1*t137*t211*t34
            t327 = -1.2248336268329735d+1*t1*t10*t105*t11*t12*t142
     +        *t221*t47
            t328 = -1.6331115024439646d+1*t1*t10*t11*t40*t47*t54*t61
     +        *t63*t91
            t329 = 8.165557512219823d+0*t1*t10*t11*t12*t126*t54*t63*t91
            t330 = 4.0827787561099116d+0*t1*t10*t11*t12*t273*t47*t54
     +        *t63*t91
            t331 = 4.0827787561099116d+0*t1*t10*t11*t12*t149*t47*t54*t91
            t332 = -8.165557512219823d+0*t1*t10*t11*t142*t39*t40*t54*t91
            t333 = 4.0827787561099116d+0*t1*t10*t11*t12*t142*t47*t54*t91
            t334 = -9.526483764256461d+0*t1*t10*t12*t142*t41*t54*t91
            t335 = -4.899334507331894d+1*t1*t10*t11*t136*t16*t39*t51
            t336 = 2.449667253665947d+1*t1*t10*t11*t118*t16*t22*t61
            t337 = 1.2248336268329735d+1*t1*t10*t11*t141*t16*t22*t39
            t338 = -2.8579451292769387d+1*t1*t10*t136*t16*t22*t41
            t339 = -4.0827787561099116d+0*t1*t10*t11*t16*t326*t40
            t340 = 9.526483764256461d+0*t1*t10*t141*t16*t40*t41
            t341 = -6.350989176170975d+1*t1*t10*t16*t40*t61*t88
            t342 = -3.1754945880854873d+1*t1*t10*t16*t39*t40*t88
            t343 = -1.2248336268329735d+1*t1*t10*t11*t136*t14*t22*t23
     +        *t47
            t344 = 4.0827787561099116d+0*t1*t10*t11*t14*t141*t23*t40*t47
            t345 = -1.5877472940427437d+1*t1*t10*t12*t14*t23*t47*t88
            t346 = 8.165557512219823d+0*t1*t10*t11*t126*t14*t23*t40*t61
            t347 = -7.719352426299558d+2*t39*t4*t58*t61*t99
     +        -7.719352426299558d+2*t118*s1f*t4*t61*t99
     +        -3.859676213149779d+2*t141*s1f*t39*t4*t99
     +        -3.859676213149779d+2*t136*t30*t4*t99
     +        +3.859676213149779d+2*t136*s1f*t28*t99
     +        +3.859676213149779d+2*s1f*t45*t61*t73
     +        -1.2865587377165932d+2*t13*t58*t73+1.9298381065748896d+2
     +        *s1f*t39*t45*t73-6.432793688582966d+1*t13*t30*t73
     +        -1.9298381065748896d+2*t28*t45*t58*t61
     +        +1.9298381065748896d+2*t117*t4*t45*t61
     +        +1.9298381065748896d+2*t118*t4*t45*t58
     +        +9.649190532874448d+1*t135*t39*t4*t45
     +        +9.649190532874448d+1*s1f*t326*t4*t45
     +        +9.649190532874448d+1*t141*t30*t4*t45
     +        -9.649190532874448d+1*t141*s1f*t28*t45
     +        +1.9298381065748896d+3*t136*s1f*t244*t39*t4
     +        -3.216396844291483d+1*t13*t325*t4+t305+t304+t303+t302
     +        +t301+3.216396844291483d+1*t13*t135*t28+t245
            t348 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t347
            t349 = -2.449667253665947d+1*t1*t10*t11*t14*t22*t23*t39
     +        *t61*t63
            t350 = 8.165557512219823d+0*t1*t10*t11*t118*t14*t23*t40*t63
            t351 = -1.9052967528512924d+1*t1*t10*t14*t23*t40*t41*t61*t63
            t352 = 8.165557512219823d+0*t1*t10*t11*t14*t23*t40*t47*t61
     +        *t63
            t353 = -4.0827787561099116d+0*t1*t10*t11*t12*t126*t14*t23
     +        *t63
            t354 = -3.1754945880854873d+1*t1*t10*t12*t14*t23*t63*t88
            t355 = 4.0827787561099116d+0*t1*t10*t11*t14*t149*t23*t39*t40
            t356 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t149*t23
     +        *t47
            t357 = 4.763241882128231d+0*t1*t10*t12*t14*t149*t23*t41
            t358 = 4.0827787561099116d+0*t1*t10*t11*t14*t142*t23*t39*t40
            t359 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t142*t23
     +        *t47
            t360 = 4.763241882128231d+0*t1*t10*t12*t14*t142*t23*t41
            t361 = s1a3(D3_RB_RB_RB)
            t362 = t61**3
            t363 = -3.333333333333333d-1*t140*t60*t83
     +        -3.333333333333333d-1*t138*t59*t78+1.4814814814814814d-1
     +        *t214*t60**3+1.4814814814814814d-1*t211*t59**3
     +        +3.333333333333333d-1*(t215+t213)*t38
     +        +3.333333333333333d-1*(t216+t212)*t35
            t364 = t63**3
            t365 = -1.2248336268329735d+1*t1*t10*t105*t11*t12*t221*t364
            t366 = 1.2248336268329735d+1*t1*t10*t11*t12*t149*t54*t63*t91
            t367 = -2.449667253665947d+1*t1*t10*t11*t142*t40*t54*t61*t91
            t368 = 4.0827787561099116d+0*t1*t10*t11*t12*t142*t273*t54
     +        *t91
            t369 = -2.8579451292769387d+1*t1*t10*t12*t142*t41*t54*t91
            t370 = 4.0827787561099116d+0*t1*t10*t11*t12*t364*t54*t91
            t371 = -4.899334507331894d+1*t1*t10*t11*t16*t362*t51
            t372 = 3.6745008804989204d+1*t1*t10*t11*t141*t16*t22*t61
            t373 = -8.573835387830815d+1*t1*t10*t136*t16*t22*t41
            t374 = -4.0827787561099116d+0*t1*t10*t11*t16*t363*t40
            t375 = 2.8579451292769387d+1*t1*t10*t141*t16*t40*t41
            t376 = -9.526483764256461d+1*t1*t10*t16*t40*t61*t88
            t377 = -3.6745008804989204d+1*t1*t10*t11*t136*t14*t22*t23
     +        *t63
            t378 = 1.2248336268329735d+1*t1*t10*t11*t14*t141*t23*t40*t63
            t379 = -5.715890258553877d+1*t1*t10*t14*t23*t40*t41*t61*t63
            t380 = -4.763241882128231d+1*t1*t10*t12*t14*t23*t63*t88
            t381 = 1.2248336268329735d+1*t1*t10*t11*t14*t149*t23*t40*t61
            t382 = -6.124168134164868d+0*t1*t10*t11*t12*t14*t149*t23*t63
            t383 = 1.4289725646384693d+1*t1*t10*t12*t14*t149*t23*t41
            t384 = -1.1579028639449339d+3*t141*s1f*t4*t61*t99
     +        -1.1579028639449339d+3*t136*t4*t58*t99
     +        +1.1579028639449339d+3*t136*s1f*t28*t99
     +        +5.789514319724669d+2*s1f*t45*t61*t73
     +        -1.9298381065748896d+2*t13*t58*t73-5.789514319724669d+2
     +        *t28*t45*t58*t61+2.8947571598623345d+2*t135*t4*t45*t61
     +        +2.8947571598623345d+2*t141*t4*t45*t58
     +        +9.649190532874448d+1*s1f*t363*t4*t45
     +        -2.8947571598623345d+2*t141*s1f*t28*t45
     +        +1.9298381065748896d+3*s1f*t244*t362*t4
     +        -3.216396844291483d+1*t13*t361*t4+9.649190532874448d+1
     +        *t13*t135*t28+t245
            t385 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t384
            t386 = 1.2248336268329735d+1*t1*t10*t11*t14*t142*t23*t40*t61
            t387 = 1.4289725646384693d+1*t1*t10*t12*t14*t142*t23*t41
            t388 = -2.0413893780549558d+0*t1*t10*t11*t12*t14*t23*t364
            t389 = -1.5877472940427437d+1*t1*t10*t12*t18*t50*t68*t88
            t390 = 1.9052967528512924d+1*t1*t18*t25*t39*t40*t41
            t391 = 1.5877472940427437d+1*t1*t12*t18*t25*t88
            t392 = 3.889452553313227d+1*t10*t16*t19*t25*t39*t51*t52
            t393 = 1/t3**6.666666666666667d+0
            t394 = 3.241210461094356d+1*t10*t16*t19*t22*t25*t393
            t395 = -9.526483764256461d+0*t1*t12*t25*t41*t49
            t396 = 9.526483764256461d+0*t1*t12*t18*t41*t50*t55
            t397 = 1.9447262766566134d+1*t10*t16*t19*t22*t50*t52*t55
            t398 = 4.0827787561099116d+0*t1*t11*t12*t54*t90*t91
            t399 = 1.2248336268329735d+1*t1*t11*t16*t22*t76
            t400 = -4.0827787561099116d+0*t1*t11*t16*t40*t87
            t401 = 1.9052967528512924d+1*t1*t16*t39*t40*t41
            t402 = 1.5877472940427437d+1*t1*t12*t16*t88
            t403 = 8.165557512219823d+0*t1*t11*t14*t23*t39*t40*t47
            t404 = 9.526483764256461d+0*t1*t12*t14*t23*t41*t47
            t405 = -2.0413893780549558d+0*t1*t101*t11*t12*t14*t23
            t406 = -2.0413893780549558d+0*t1*t11*t12*t14*t23*t90
            t407 = 3.1115620426505813d+2*t10*t19*t23*t39*t51*t52
            t408 = 2.2040231135441618d+2*t10*t107*t19*t22*t23
            t409 = 1.5557810213252907d+2*t10*t14*t19*t22*t47*t52*t54
            t410 = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (5.00072471140272d+1*t10*t105*t19*t21*t22*t90*t91
     +        -1.6669082371342397d+1*t10*t14*t19*t21*t22*t54*t90
     +        -3.3338164742684795d+1*t10*t19*t21*t23*t51*t87
     +        +1.6669082371342397d+2*t10*t106*t19*t21*t23*t76
     +        +1.333526589707392d+2*t10*t14*t19*t21*t39*t47*t51*t54
     +        -1.6669082371342397d+1*t10*t101*t14*t19*t21*t22*t54+t409
     +        +t408+t407+t406+t405+t404+t403+t402+t401+t400+t399+t398)
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t68*t87
     +        -8.334541185671199d+0*t10*t16*t19*t21*t25*t51*t87
     +        -4.0827787561099116d+0*t1*t11*t18*t25*t40*t87
     +        -1.2248336268329735d+1*t1*t10*t11*t18*t22*t50*t68*t76
     +        +1.2248336268329735d+1*t1*t11*t18*t22*t25*t76
     +        +2.50036235570136d+1*t10*t106*t16*t19*t21*t25*t76
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t49*t55*t68
     +        -1.9052967528512924d+1*t1*t10*t110*t12*t18*t41*t55*t68
     +        -1.6331115024439646d+1*t1*t10*t11*t110*t18*t39*t40*t55
     +        *t68+9.526483764256461d+0*t1*t10*t12*t41*t49*t50*t68
     +        +8.165557512219823d+0*t1*t10*t11*t39*t40*t49*t50*t68
     +        -1.9052967528512924d+1*t1*t10*t18*t39*t40*t41*t50*t68
     +        -2.0413893780549558d+0*t1*t10*t104*t11*t12*t50*t68
     +        -1.2248336268329735d+1*t1*t10*t11*t111*t12*t18*t254*t68
     +        +4.0827787561099116d+0*t1*t10*t109*t11*t110*t12*t18*t68
     +        +1.6669082371342397d+1*t10*t16*t19*t21*t39*t50*t51*t55
     +        -4.0827787561099116d+0*t1*t11*t12*t49*t50*t55
     +        +8.165557512219823d+0*t1*t11*t18*t39*t40*t50*t55
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t164*t50*t55
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t166*t18*t55
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t166*t49*t50
     +        +9.526483764256461d+0*t1*t10*t12*t166*t18*t41*t50
     +        +8.165557512219823d+0*t1*t10*t11*t166*t18*t39*t40*t50
     +        -4.1672705928355996d+0*t10*t109*t16*t19*t21*t22*t50
     +        -2.0413893780549558d+0*t1*t109*t11*t12*t18*t50
     +        -8.165557512219823d+0*t1*t11*t25*t39*t40*t49
     +        -9.526483764256461d+0*t1*t10*t12*t164*t25*t41
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t406+t405
     +        +t404+t403+t402+t401+t400+t399+t398)
     +        -8.165557512219823d+0*t1*t10*t11*t164*t25*t39*t40+t397
     +        +t396+t395+t394+t392+t391+t390+t389
     +        +2.0413893780549558d+0*t1*t104*t11*t12*t25
     +        +8.334541185671199d+0*t10*t110*t111*t16*t19*t21*t22
     +        +4.0827787561099116d+0*t1*t11*t110*t111*t12*t18)
     +        +9.327207260896466d-2*t57*t69*t87+1.8654414521792934d-1
     +        *t69*t76*t8+1.8654414521792934d-1*t167*t39*t57)+t114*(
     +        -6.218138173930978d-2*t167*t56*t9-t113*t69)
     +        +6.218138173930978d-2*t115*t255*t69*t9)*wght
            t411 = -1.5877472940427437d+1*t1*t10*t12*t18*t50*t71*t88
            t412 = 3.1754945880854873d+1*t1*t12*t18*t25*t88
            t413 = 6.482420922188712d+1*t10*t16*t19*t22*t25*t393
            t414 = 8.165557512219823d+0*t1*t11*t12*t54*t90*t91
            t415 = 2.449667253665947d+1*t1*t11*t16*t22*t76
            t416 = -8.165557512219823d+0*t1*t11*t16*t40*t87
            t417 = 3.8105935057025847d+1*t1*t16*t39*t40*t41
            t418 = 3.1754945880854873d+1*t1*t12*t16*t88
            t419 = 1.6331115024439646d+1*t1*t11*t14*t23*t39*t40*t47
            t420 = 1.9052967528512924d+1*t1*t12*t14*t23*t41*t47
            t421 = -4.0827787561099116d+0*t1*t101*t11*t12*t14*t23
            t422 = -4.0827787561099116d+0*t1*t11*t12*t14*t23*t90
            t423 = 4.4080462270883236d+2*t10*t107*t19*t22*t23
            t424 = 4.0827787561099116d+0*t1*t11*t12*t47*t54*t63*t91
            t425 = 1.2248336268329735d+1*t1*t11*t16*t22*t39*t61
            t426 = -4.0827787561099116d+0*t1*t11*t118*t16*t40
            t427 = 9.526483764256461d+0*t1*t16*t40*t41*t61
            t428 = 9.526483764256461d+0*t1*t16*t39*t40*t41
            t429 = 4.0827787561099116d+0*t1*t11*t14*t23*t40*t47*t61
            t430 = 4.763241882128231d+0*t1*t12*t14*t23*t41*t47
            t431 = -2.0413893780549558d+0*t1*t11*t12*t126*t14*t23
            t432 = 4.0827787561099116d+0*t1*t11*t14*t23*t39*t40*t63
            t433 = -2.0413893780549558d+0*t1*t11*t12*t14*t23*t47*t63
            t434 = 4.763241882128231d+0*t1*t12*t14*t23*t41*t63
            t435 = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (5.00072471140272d+1*t10*t105*t19*t21*t22*t47*t63*t91
     +        +7.778905106626454d+1*t10*t14*t19*t22*t52*t54*t63
     +        +6.667632948536959d+1*t10*t14*t19*t21*t39*t51*t54*t63
     +        -1.6669082371342397d+1*t10*t14*t19*t21*t22*t47*t54*t63
     +        +6.667632948536959d+1*t10*t14*t19*t21*t47*t51*t54*t61
     +        +1.5557810213252907d+2*t10*t19*t23*t51*t52*t61
     +        +1.6669082371342397d+2*t10*t106*t19*t21*t23*t39*t61
     +        +7.778905106626454d+1*t10*t14*t19*t22*t47*t52*t54
     +        -1.6669082371342397d+1*t10*t126*t14*t19*t21*t22*t54
     +        +1.5557810213252907d+2*t10*t19*t23*t39*t51*t52
     +        -3.3338164742684795d+1*t10*t118*t19*t21*t23*t51+t434
     +        +t433+t432+t431+t430+t429+t428+t427+t426+t425+t424+t408
     +        +t402)-1.2248336268329735d+1*t1*t10*t11*t12*t18*t254*t55
     +        *t66*t68+4.0827787561099116d+0*t1*t10*t11*t110*t12*t49
     +        *t66*t68-9.526483764256461d+0*t1*t10*t110*t12*t18*t41
     +        *t66*t68-8.165557512219823d+0*t1*t10*t11*t110*t18*t39
     +        *t40*t66*t68+4.0827787561099116d+0*t1*t10*t11*t110*t12
     +        *t55*t65*t68+4.763241882128231d+0*t1*t10*t12*t41*t50*t65
     +        *t68+4.0827787561099116d+0*t1*t10*t11*t39*t40*t50*t65
     +        *t68-8.165557512219823d+0*t1*t10*t11*t110*t18*t40*t55
     +        *t61*t68+4.0827787561099116d+0*t1*t10*t11*t40*t49*t50
     +        *t61*t68-9.526483764256461d+0*t1*t10*t18*t40*t41*t50*t61
     +        *t68-1.2248336268329735d+1*t1*t10*t11*t18*t22*t39*t50
     +        *t61*t68-9.526483764256461d+0*t1*t10*t110*t12*t18*t41
     +        *t55*t68+4.763241882128231d+0*t1*t10*t12*t41*t49*t50*t68
     +        -9.526483764256461d+0*t1*t10*t18*t39*t40*t41*t50*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t118*t18*t40*t50*t68
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t131*t50*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t132*t18*t68
     +        +8.334541185671199d+0*t10*t110*t16*t19*t21*t22*t55*t66
     +        +4.0827787561099116d+0*t1*t11*t110*t12*t18*t55*t66
     +        +9.723631383283067d+0*t10*t16*t19*t22*t50*t52*t66
     +        +8.334541185671199d+0*t10*t16*t19*t21*t39*t50*t51*t66
     +        -2.0413893780549558d+0*t1*t11*t12*t49*t50*t66
     +        +4.763241882128231d+0*t1*t12*t18*t41*t50*t66
     +        +4.0827787561099116d+0*t1*t11*t18*t39*t40*t50*t66
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t164*t50*t66
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t166*t18*t66
     +        -2.0413893780549558d+0*t1*t11*t12*t50*t55*t65
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t166*t50*t65
     +        -4.763241882128231d+0*t1*t12*t25*t41*t65
     +        -4.0827787561099116d+0*t1*t11*t25*t39*t40*t65
     +        +8.334541185671199d+0*t10*t16*t19*t21*t50*t51*t55*t61
     +        +4.0827787561099116d+0*t1*t11*t18*t40*t50*t55*t61
     +        +1.9447262766566134d+1*t10*t16*t19*t25*t51*t52*t61
     +        +4.0827787561099116d+0*t1*t10*t11*t166*t18*t40*t50*t61
     +        -4.0827787561099116d+0*t1*t11*t25*t40*t49*t61
     +        +9.526483764256461d+0*t1*t18*t25*t40*t41*t61
     +        -4.0827787561099116d+0*t1*t10*t11*t164*t25*t40*t61
     +        +1.2248336268329735d+1*t1*t11*t18*t22*t25*t39*t61
     +        +2.50036235570136d+1*t10*t106*t16*t19*t21*t25*t39*t61
     +        +9.723631383283067d+0*t10*t16*t19*t22*t50*t52*t55
     +        +4.763241882128231d+0*t1*t12*t18*t41*t50*t55
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t183*t50*t55
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t184*t55
     +        +1.9447262766566134d+1*t10*t16*t19*t25*t39*t51*t52
     +        -8.334541185671199d+0*t10*t118*t16*t19*t21*t25*t51
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t184*t49*t50
     +        +4.763241882128231d+0*t1*t10*t12*t18*t184*t41*t50
     +        +4.763241882128231d+0*t1*t10*t12*t166*t18*t41*t50
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t184*t39*t40*t50
     +        -4.1672705928355996d+0*t10*t132*t16*t19*t21*t22*t50
     +        -2.0413893780549558d+0*t1*t11*t12*t132*t18*t50
     +        -4.763241882128231d+0*t1*t12*t25*t41*t49
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t434+t433
     +        +t432+t431+t430+t429+t428+t427+t426+t425+t424+t402)
     +        +9.526483764256461d+0*t1*t18*t25*t39*t40*t41
     +        -4.763241882128231d+0*t1*t10*t12*t183*t25*t41
     +        -4.763241882128231d+0*t1*t10*t12*t164*t25*t41
     +        -4.0827787561099116d+0*t1*t10*t11*t183*t25*t39*t40
     +        -4.0827787561099116d+0*t1*t11*t118*t18*t25*t40+t394+t391
     +        +t389+2.0413893780549558d+0*t1*t11*t12*t131*t25)
     +        +1.8654414521792934d-1*t39*t61*t69*t8
     +        +9.327207260896466d-2*t118*t57*t69+9.327207260896466d-2
     +        *t167*t57*t61+9.327207260896466d-2*t185*t39*t57)+t114*(
     +        -3.109069086965489d-2*t167*t67*t9-3.109069086965489d-2
     +        *t185*t56*t9-t134*t69)+6.218138173930978d-2*t255*t56*t67
     +        *t69*t9)*wght
            t436 = 1.9052967528512924d+1*t1*t18*t25*t40*t41*t61
            t437 = 3.889452553313227d+1*t10*t16*t19*t25*t51*t52*t61
            t438 = 8.165557512219823d+0*t1*t11*t12*t47*t54*t63*t91
            t439 = 2.449667253665947d+1*t1*t11*t16*t22*t39*t61
            t440 = -8.165557512219823d+0*t1*t11*t118*t16*t40
            t441 = 1.9052967528512924d+1*t1*t16*t40*t41*t61
            t442 = 8.165557512219823d+0*t1*t11*t14*t23*t40*t47*t61
            t443 = -4.0827787561099116d+0*t1*t11*t12*t126*t14*t23
            t444 = 8.165557512219823d+0*t1*t11*t14*t23*t39*t40*t63
            t445 = -4.0827787561099116d+0*t1*t11*t12*t14*t23*t47*t63
            t446 = 9.526483764256461d+0*t1*t12*t14*t23*t41*t63
            t447 = 3.1115620426505813d+2*t10*t19*t23*t51*t52*t61
            t448 = 1.5557810213252907d+2*t10*t14*t19*t22*t52*t54*t63
            t449 = -9.526483764256461d+0*t1*t12*t25*t41*t65
            t450 = 9.526483764256461d+0*t1*t12*t18*t41*t50*t66
            t451 = 1.9447262766566134d+1*t10*t16*t19*t22*t50*t52*t66
            t452 = 4.0827787561099116d+0*t1*t11*t12*t142*t54*t91
            t453 = 1.2248336268329735d+1*t1*t11*t136*t16*t22
            t454 = -4.0827787561099116d+0*t1*t11*t141*t16*t40
            t455 = 8.165557512219823d+0*t1*t11*t14*t23*t40*t61*t63
            t456 = -2.0413893780549558d+0*t1*t11*t12*t14*t149*t23
            t457 = -2.0413893780549558d+0*t1*t11*t12*t14*t142*t23
            t458 = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (5.00072471140272d+1*t10*t105*t142*t19*t21*t22*t91
     +        +1.333526589707392d+2*t10*t14*t19*t21*t51*t54*t61*t63
     +        -1.6669082371342397d+1*t10*t14*t149*t19*t21*t22*t54
     +        -1.6669082371342397d+1*t10*t14*t142*t19*t21*t22*t54
     +        -3.3338164742684795d+1*t10*t141*t19*t21*t23*t51+t457
     +        +t456+t455+t454+t453+t452+t448+t447+t446+t441+t408+t402
     +        +1.6669082371342397d+2*t10*t106*t136*t19*t21*t23)
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t65*t66*t68
     +        -1.6331115024439646d+1*t1*t10*t11*t110*t18*t40*t61*t66
     +        *t68-1.9052967528512924d+1*t1*t10*t110*t12*t18*t41*t66
     +        *t68+8.165557512219823d+0*t1*t10*t11*t40*t50*t61*t65*t68
     +        +9.526483764256461d+0*t1*t10*t12*t41*t50*t65*t68
     +        -1.9052967528512924d+1*t1*t10*t18*t40*t41*t50*t61*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t141*t18*t40*t50*t68
     +        -1.2248336268329735d+1*t1*t10*t11*t136*t18*t22*t50*t68
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t152*t50*t68
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t154*t18*t254*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t153*t18*t68
     +        -4.0827787561099116d+0*t1*t11*t12*t50*t65*t66
     +        +1.6669082371342397d+1*t10*t16*t19*t21*t50*t51*t61*t66
     +        +8.165557512219823d+0*t1*t11*t18*t40*t50*t61*t66
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t183*t50*t66
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t18*t184*t66
     +        -8.165557512219823d+0*t1*t11*t25*t40*t61*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t184*t50*t65
     +        +8.165557512219823d+0*t1*t10*t11*t18*t184*t40*t50*t61
     +        -8.165557512219823d+0*t1*t10*t11*t183*t25*t40*t61
     +        -8.334541185671199d+0*t10*t141*t16*t19*t21*t25*t51
     +        +9.526483764256461d+0*t1*t10*t12*t18*t184*t41*t50
     +        -4.1672705928355996d+0*t10*t153*t16*t19*t21*t22*t50
     +        -2.0413893780549558d+0*t1*t11*t12*t153*t18*t50
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t457+t456
     +        +t455+t454+t453+t452+t446+t441+t402)+t451+t450+t449+t437
     +        +t436-9.526483764256461d+0*t1*t10*t12*t183*t25*t41
     +        -4.0827787561099116d+0*t1*t11*t141*t18*t25*t40+t394+t391
     +        +t389+1.2248336268329735d+1*t1*t11*t136*t18*t22*t25
     +        +2.50036235570136d+1*t10*t106*t136*t16*t19*t21*t25
     +        +2.0413893780549558d+0*t1*t11*t12*t152*t25
     +        +8.334541185671199d+0*t10*t110*t154*t16*t19*t21*t22
     +        +4.0827787561099116d+0*t1*t11*t110*t12*t154*t18)
     +        +1.8654414521792934d-1*t136*t69*t8+9.327207260896466d-2
     +        *t141*t57*t69+1.8654414521792934d-1*t185*t57*t61)+t114*(
     +        -6.218138173930978d-2*t185*t67*t9-t156*t69)
     +        +6.218138173930978d-2*t157*t255*t69*t9)*wght
            t459 = 8.165557512219823d+0*t1*t11*t12*t142*t54*t91
            t460 = 2.449667253665947d+1*t1*t11*t136*t16*t22
            t461 = -8.165557512219823d+0*t1*t11*t141*t16*t40
            t462 = 3.8105935057025847d+1*t1*t16*t40*t41*t61
            t463 = 1.6331115024439646d+1*t1*t11*t14*t23*t40*t61*t63
            t464 = 1.9052967528512924d+1*t1*t12*t14*t23*t41*t63
            t465 = -4.0827787561099116d+0*t1*t11*t12*t14*t149*t23
            t466 = -4.0827787561099116d+0*t1*t11*t12*t14*t142*t23
            t467 = -9.526483764256461d+0*t1*t10*t110*t12*t18*t195*t41
            t468 = 8.165557512219823d+0*t1*t11*t18*t39*t40*t50*t68
            t469 = 9.526483764256461d+0*t1*t12*t18*t41*t50*t68
            t470 = 1/t3**8.0d+0
            t471 = 3.969943564391174d+1*t10*t106*t18*t197*t23*t470*t50
            t472 = 1.6669082371342397d+1*t10*t16*t19*t21*t39*t50*t51*t68
            t473 = 1.9447262766566134d+1*t10*t16*t19*t22*t50*t52*t68
            t474 = -1.9447262766566134d+1*t16*t19*t22*t25*t52
            t475 = -3.889452553313227d+1*t19*t22*t23*t52
            t476 = 4.0827787561099116d+0*t1*t11*t12*t164*t25
            t477 = -4.0827787561099116d+0*t1*t11*t12*t166*t18*t50
            t478 = -8.334541185671199d+0*t10*t16*t166*t19*t21*t22*t50
            t479 = -4.0827787561099116d+0*t1*t11*t12*t49*t50*t68
            t480 = 8.165557512219823d+0*t1*t11*t110*t12*t18*t55*t68
            t481 = 1.6669082371342397d+1*t10*t110*t16*t19*t21*t22*t55
     +        *t68
            t482 = (t114*(-3.109069086965489d-2*t167*t69*t9
     +        -3.109069086965489d-2*t200*t56*t9-t168*t69)+t31*
     +        (3.109069086965489d-2*(-4.0827787561099116d+0*t1*t10*t11
     +        *t12*t164*t50*t68+8.165557512219823d+0*t1*t10*t11*t110
     +        *t12*t166*t18*t68-8.334541185671199d+0*t16*t19*t21*t22
     +        *t50*t55-1.2248336268329735d+1*t1*t10*t11*t12*t18*t195
     +        *t254*t55+3.402808769478149d+1*t10*t106*t110*t18*t197
     +        *t198*t23*t55-2.0413893780549558d+0*t1*t10*t11*t12*t18
     +        *t50*(-1.6669082371342397d+1*t14*t19*t21*t22*t47*t54
     +        -3.3338164742684795d+1*t19*t21*t23*t39*t51+t475)
     +        -1.6669082371342397d+1*t16*t19*t21*t25*t39*t51
     +        -1.7014043847390745d+1*t10*t106*t197*t198*t23*t49*t50
     +        +3.402808769478149d+1*t10*t18*t197*t198*t23*t252*t39*t50
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t195*t49+t481
     +        +t480+t479+t478+t477+t476+t474+t473+t472+t471+t469+t468
     +        +t467-8.165557512219823d+0*t1*t10*t11*t110*t18*t195*t39
     +        *t40)*t9+9.327207260896466d-2*t200*t39*t57)
     +        +6.218138173930978d-2*t194*t255*t56*t9)*wght
            t483 = 6.218138173930978d-2*t255*t56*t69*t72*t9
            t484 = -9.526483764256461d+0*t1*t10*t110*t12*t18*t41*t68*t71
            t485 = 4.763241882128231d+0*t1*t12*t18*t41*t50*t71
            t486 = 7.939887128782348d+1*t10*t106*t18*t197*t23*t470*t50
            t487 = 9.723631383283067d+0*t10*t16*t19*t22*t50*t52*t71
            t488 = -3.889452553313227d+1*t16*t19*t22*t25*t52
            t489 = -7.778905106626454d+1*t19*t22*t23*t52
            t490 = t31*(3.109069086965489d-2*(-1.2248336268329735d+1
     +        *t1*t10*t11*t12*t18*t254*t55*t68*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t49*t68*t71
     +        -8.165557512219823d+0*t1*t10*t11*t110*t18*t39*t40*t68
     +        *t71+8.334541185671199d+0*t10*t110*t16*t19*t21*t22*t55
     +        *t71+4.0827787561099116d+0*t1*t11*t110*t12*t18*t55*t71
     +        +8.334541185671199d+0*t10*t16*t19*t21*t39*t50*t51*t71
     +        -2.0413893780549558d+0*t1*t11*t12*t49*t50*t71
     +        +4.0827787561099116d+0*t1*t11*t18*t39*t40*t50*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t164*t50*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t166*t18*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t176*t50*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t178*t18*t68
     +        -1.6669082371342397d+1*t16*t19*t21*t22*t50*t55
     +        +6.805617538956298d+1*t10*t106*t110*t18*t197*t198*t23
     +        *t55-2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -3.3338164742684795d+1*t14*t19*t21*t22*t47*t54
     +        -6.667632948536959d+1*t19*t21*t23*t39*t51+t489)
     +        -3.3338164742684795d+1*t16*t19*t21*t25*t39*t51
     +        -3.402808769478149d+1*t10*t106*t197*t198*t23*t49*t50
     +        +6.805617538956298d+1*t10*t18*t197*t198*t23*t252*t39*t50
     +        -4.1672705928355996d+0*t10*t16*t178*t19*t21*t22*t50
     +        -2.0413893780549558d+0*t1*t11*t12*t178*t18*t50+t488+t487
     +        +t486+t485+t484+t481+t480+t479+t478+t477+t476+t473+t472
     +        +t469+t468+2.0413893780549558d+0*t1*t11*t12*t176*t25)*t9
     +        +9.327207260896466d-2*t202*t39*t57)
            t491 = -3.109069086965489d-2*t202*t56*t9
            t492 = -9.526483764256461d+0*t1*t10*t110*t12*t18*t205*t41
            t493 = 1.9052967528512924d+1*t1*t12*t18*t41*t50*t71
            t494 = 1.5879774257564697d+2*t10*t106*t18*t197*t23*t470*t50
            t495 = 3.889452553313227d+1*t10*t16*t19*t22*t50*t52*t71
            t496 = -7.778905106626454d+1*t16*t19*t22*t25*t52
            t497 = -1.5557810213252907d+2*t19*t22*t23*t52
            t498 = 8.165557512219823d+0*t1*t11*t18*t40*t50*t61*t68
            t499 = 1.6669082371342397d+1*t10*t16*t19*t21*t50*t51*t61*t68
            t500 = 4.0827787561099116d+0*t1*t11*t12*t183*t25
            t501 = -4.0827787561099116d+0*t1*t11*t12*t18*t184*t50
            t502 = -8.334541185671199d+0*t10*t16*t184*t19*t21*t22*t50
            t503 = -4.0827787561099116d+0*t1*t11*t12*t50*t65*t68
            t504 = 8.165557512219823d+0*t1*t11*t110*t12*t18*t66*t68
            t505 = 1.6669082371342397d+1*t10*t110*t16*t19*t21*t22*t66
     +        *t68
            t506 = (t114*(-3.109069086965489d-2*t185*t69*t9
     +        -3.109069086965489d-2*t200*t67*t9-t186*t69)+t31*
     +        (3.109069086965489d-2*(-4.0827787561099116d+0*t1*t10*t11
     +        *t12*t183*t50*t68+8.165557512219823d+0*t1*t10*t11*t110
     +        *t12*t18*t184*t68-8.334541185671199d+0*t16*t19*t21*t22
     +        *t50*t66-1.2248336268329735d+1*t1*t10*t11*t12*t18*t195
     +        *t254*t66+3.402808769478149d+1*t10*t106*t110*t18*t197
     +        *t198*t23*t66-1.7014043847390745d+1*t10*t106*t197*t198
     +        *t23*t50*t65+4.0827787561099116d+0*t1*t10*t11*t110*t12
     +        *t195*t65-2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -1.6669082371342397d+1*t14*t19*t21*t22*t54*t63
     +        -3.3338164742684795d+1*t19*t21*t23*t51*t61+t475)
     +        -1.6669082371342397d+1*t16*t19*t21*t25*t51*t61
     +        +3.402808769478149d+1*t10*t18*t197*t198*t23*t252*t50*t61
     +        -8.165557512219823d+0*t1*t10*t11*t110*t18*t195*t40*t61
     +        +t505+t504+t503+t502+t501+t500+t499+t498+t474+t473+t471
     +        +t469+t467)*t9+9.327207260896466d-2*t200*t57*t61)
     +        +6.218138173930978d-2*t194*t255*t67*t9)*wght
            t507 = 6.218138173930978d-2*t255*t67*t69*t72*t9
            t508 = t31*(3.109069086965489d-2*(-1.2248336268329735d+1
     +        *t1*t10*t11*t12*t18*t254*t66*t68*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t65*t68*t71
     +        -8.165557512219823d+0*t1*t10*t11*t110*t18*t40*t61*t68
     +        *t71+8.334541185671199d+0*t10*t110*t16*t19*t21*t22*t66
     +        *t71+4.0827787561099116d+0*t1*t11*t110*t12*t18*t66*t71
     +        -2.0413893780549558d+0*t1*t11*t12*t50*t65*t71
     +        +8.334541185671199d+0*t10*t16*t19*t21*t50*t51*t61*t71
     +        +4.0827787561099116d+0*t1*t11*t18*t40*t50*t61*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t183*t50*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t184*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t190*t50*t68
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t191*t68
     +        -1.6669082371342397d+1*t16*t19*t21*t22*t50*t66
     +        +6.805617538956298d+1*t10*t106*t110*t18*t197*t198*t23
     +        *t66-3.402808769478149d+1*t10*t106*t197*t198*t23*t50*t65
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -3.3338164742684795d+1*t14*t19*t21*t22*t54*t63
     +        -6.667632948536959d+1*t19*t21*t23*t51*t61+t489)
     +        -3.3338164742684795d+1*t16*t19*t21*t25*t51*t61
     +        +6.805617538956298d+1*t10*t18*t197*t198*t23*t252*t50*t61
     +        +t505+t504+t503+t502+t501+t500-4.1672705928355996d+0*t10
     +        *t16*t19*t191*t21*t22*t50-2.0413893780549558d+0*t1*t11
     +        *t12*t18*t191*t50+t499+t498+t488+t487+t486+t485+t484
     +        +t473+t469+2.0413893780549558d+0*t1*t11*t12*t190*t25)*t9
     +        +9.327207260896466d-2*t202*t57*t61)
            t509 = -3.109069086965489d-2*t202*t67*t9
            t510 = t1**4
            t511 = 1/t3**9.333333333333333d+0
            t512 = 1/t8**8
            t513 = 6.218138173930978d-2*t255*t69**3*t9*wght
     +        -9.327207260896466d-2*t114*t200*t69*t9*wght
     +        +3.109069086965489d-2*t31*(-1.2248336268329735d+1*t1*t10
     +        *t11*t12*t18*t254*t68**3-2.50036235570136d+1*t16*t19*t21
     +        *t22*t50*t68+1.0208426308434447d+2*t10*t106*t110*t18
     +        *t197*t198*t23*t68-1.0419686516347423d+2*t10*t50*t510
     +        *t511*t512*t54-5.104213154217223d+1*t106*t18*t197*t198
     +        *t23*t50+2.50036235570136d+1*t10*t110*t16*t19*t195*t21
     +        *t22+1.2248336268329735d+1*t1*t11*t110*t12*t18*t195)*t9
     +        *wght
            t514 = -1.2248336268329735d+1*t1*t10*t11*t12*t18*t195*t254
     +        *t71
            t515 = 8.165557512219823d+0*t1*t11*t110*t12*t18*t195
            t516 = 8.165557512219823d+0*t1*t11*t110*t12*t18*t68*t71
            t517 = 1.3611235077912595d+2*t10*t106*t110*t18*t197*t198
     +        *t23*t68
            t518 = 3.402808769478149d+1*t10*t106*t110*t18*t197*t198
     +        *t23*t71
            t519 = -1.0208426308434447d+2*t106*t18*t197*t198*t23*t50
            t520 = 1.6669082371342397d+1*t10*t110*t16*t19*t195*t21*t22
            t521 = 1.6669082371342397d+1*t10*t110*t16*t19*t21*t22*t68
     +        *t71
            t522 = -3.3338164742684795d+1*t16*t19*t21*t22*t50*t68
            t523 = -8.334541185671199d+0*t16*t19*t21*t22*t50*t71
            t524 = 6.218138173930978d-2*t194*t255*t72*t9*wght
            t525 = t114*(-3.109069086965489d-2*t200*t72*t9*wght
     +        -6.218138173930978d-2*t202*t69*t9*wght)
            t526 = 6.218138173930978d-2*t204*t255*t69*t9*wght
            t527 = -1.2248336268329735d+1*t1*t10*t11*t12*t18*t205*t254
     +        *t68
            t528 = 4.0827787561099116d+0*t1*t11*t110*t12*t18*t205
            t529 = 1.6331115024439646d+1*t1*t11*t110*t12*t18*t68*t71
            t530 = 1.3611235077912595d+2*t10*t106*t110*t18*t197*t198
     +        *t23*t71
            t531 = -2.0416852616868894d+2*t106*t18*t197*t198*t23*t50
            t532 = 8.334541185671199d+0*t10*t110*t16*t19*t205*t21*t22
            t533 = 3.3338164742684795d+1*t10*t110*t16*t19*t21*t22*t68
     +        *t71
            t534 = -3.3338164742684795d+1*t16*t19*t21*t22*t50*t71
            t535 = t114*(-6.218138173930978d-2*t202*t72*t9*wght
     +        -3.109069086965489d-2*t206*t69*t9*wght)
            t536 = 3.109069086965489d-2*t31*(-2.0839373032694844d+2
     +        *t10*t50*t510*t511*t512*t54+t523+t522+t521+t520+t519
     +        +t518+t517+t516+t515+t514)*t9*wght+t525+t524
            fnc(iq) = (3.109069086965489d-2*t27*t9+1.0d+0*s1f*t4)*wght
     +        +fnc(iq)
            Amat(iq,D1_RA) = (3.109069086965489d-2*t31*t56*t9
     +        +9.327207260896466d-2*t27*t39*t57+1.0d+0*t30*t4+t29)
     +        *wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (3.109069086965489d-2*t31*t67*t9
     +        +9.327207260896466d-2*t27*t57*t61+1.0d+0*t4*t58+t29)
     +        *wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t70+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 3.109069086965489d-2*t31*t72*t9*wght
     +        +Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = t70+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (-3.109069086965489d-2*t114*t115*t9
     +        +9.327207260896466d-2*t27*t57*t87+1.8654414521792934d-1
     +        *t27*t76*t8+1.0d+0*t4*t75+t74+t113*t31-2.0d+0*t28*t30)
     +        *wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = (-3.109069086965489d-2*t114*t56*t67
     +        *t9+1.8654414521792934d-1*t27*t39*t61*t8+t74
     +        +9.327207260896466d-2*t118*t27*t57+1.0d+0*t117*t4+t134
     +        *t31+t116*t28)*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-3.109069086965489d-2*t114*t157*t9
     +        +1.8654414521792934d-1*t136*t27*t8+t74-2.0d+0*t28*t58
     +        +9.327207260896466d-2*t141*t27*t57+1.0d+0*t135*t4+t156
     +        *t31)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t169+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = (t180*t31-3.109069086965489d-2*t114
     +        *t56*t72*t9)*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t169+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t187+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = (t193*t31-3.109069086965489d-2*t114
     +        *t67*t72*t9)*wght+Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t187+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t201+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t203+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t201+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 3.109069086965489d-2*t206*t31*t9
     +        *wght-3.109069086965489d-2*t114*t204*t9*wght+Cmat2(iq,
     +        D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = t203+Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t201+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -3.500507297981904d+2*t105*t19*t20*t22*t52*t90*t91
     +        -3.0004348268416314d+2*t105*t19*t20*t21*t39*t51*t90*t91
     +        +2.50036235570136d+1*t105*t19*t20*t21*t22*t226*t90*t91
     +        +7.50108706710408d+1*t101*t105*t19*t20*t21*t22*t47*t91
     +        +2.50036235570136d+1*t105*t19*t20*t21*t22*t220*t91
     +        +1.1668357659939681d+2*t14*t19*t20*t22*t52*t54*t90
     +        +1.0001449422805439d+2*t14*t19*t20*t21*t39*t51*t54*t90
     +        +1.0001449422805439d+2*t14*t19*t20*t21*t47*t51*t54*t87
     +        +2.3336715319879361d+2*t19*t20*t23*t51*t52*t87
     +        +2.50036235570136d+2*t106*t19*t20*t21*t23*t39*t87
     +        -5.00072471140272d+2*t106*t14*t19*t20*t21*t47*t54*t76
     +        -1.166835765993968d+3*t106*t19*t20*t23*t52*t76
     +        -9.334686127951745d+2*t14*t19*t20*t39*t47*t51*t52*t54
     +        +1.1668357659939681d+2*t101*t14*t19*t20*t22*t52*t54
     +        +1.0001449422805439d+2*t101*t14*t19*t20*t21*t39*t51*t54
     +        -2.50036235570136d+1*t101*t14*t19*t20*t21*t22*t47*t54
     +        -6.612069340632486d+2*t107*t14*t19*t20*t22*t47*t54
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t246*t54
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t220*t54
     +        -1.322413868126497d+3*t107*t19*t20*t23*t39*t51
     +        -1.6669082371342397d+1*t19*t20*t21*t217*t23*t51+t253
     +        -5.00072471140272d+2*t19*t20*t21*t210*t23*t252
     +        -1.0001449422805439d+2*t19*t20*t21*t22*t220*t221*t251
     +        +t250+t249+t248+t247+t243+t242+t241+t240+t239+t238+t237
     +        +t236+t235+t234+t233+t232+t231+t230+t229+t228+t227+t224
     +        +t223+t222)-4.763241882128231d+1*t1*t10*t12*t18*t50*t55
     +        *t88+4.763241882128231d+1*t1*t10*t12*t25*t49*t88
     +        -9.526483764256461d+1*t1*t10*t18*t25*t39*t40*t88
     +        +1.2248336268329735d+1*t1*t10*t11*t18*t40*t50*t55*t87
     +        -1.2248336268329735d+1*t1*t10*t11*t25*t40*t49*t87
     +        +2.8579451292769387d+1*t1*t10*t18*t25*t40*t41*t87
     +        +3.6745008804989204d+1*t1*t10*t11*t18*t22*t25*t39*t87
     +        -3.6745008804989204d+1*t1*t10*t11*t18*t22*t50*t55*t76
     +        +3.6745008804989204d+1*t1*t10*t11*t22*t25*t49*t76
     +        -8.573835387830815d+1*t1*t10*t18*t22*t25*t41*t76
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t18*t254*t55**3
     +        +2.8579451292769387d+1*t1*t10*t12*t41*t49*t50*t55
     +        +2.449667253665947d+1*t1*t10*t11*t39*t40*t49*t50*t55
     +        -5.715890258553877d+1*t1*t10*t18*t39*t40*t41*t50*t55
     +        -6.124168134164868d+0*t1*t10*t104*t11*t12*t50*t55
     +        +1.2248336268329735d+1*t1*t10*t109*t11*t110*t12*t18*t55
     +        -4.899334507331894d+1*t1*t10*t11*t18*t210*t25*t51
     +        -6.124168134164868d+0*t1*t10*t109*t11*t12*t49*t50
     +        +1.4289725646384693d+1*t1*t10*t109*t12*t18*t41*t50
     +        +1.2248336268329735d+1*t1*t10*t109*t11*t18*t39*t40*t50
     +        +5.715890258553877d+1*t1*t10*t25*t39*t40*t41*t49
     +        +1.2248336268329735d+1*t1*t10*t11*t110*t111*t12*t49
     +        -1.4289725646384693d+1*t1*t10*t104*t12*t25*t41
     +        -2.8579451292769387d+1*t1*t10*t110*t111*t12*t18*t41
     +        -1.2248336268329735d+1*t1*t10*t104*t11*t25*t39*t40
     +        -2.449667253665947d+1*t1*t10*t11*t110*t111*t18*t39*t40
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t217*t25*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t250+t249
     +        +t248+t247+t243+t242+t241+t240+t239+t238+t237+t236+t235
     +        +t234+t233+t232+t231+t230+t229+t228+t227+t224+t223+t222)
     +        +t219)+2.7981621782689403d-1*t56*t57*t87
     +        +5.596324356537881d-1*t56*t76*t8+2.7981621782689403d-1
     +        *t112*t39*t57)+t114*(-6.218138173930978d-2*t112*t56*t9
     +        -9.327207260896466d-2*t115*t39*t57-t113*t56)
     +        +6.218138173930978d-2*t255*t56**3*t9
     +        +5.596324356537881d-1*t27*t39*t8*t87-3.0d+0*t28*t75
     +        +6.0d+0*t30*t73+9.327207260896466d-2*t217*t27*t57+1.0d+0
     +        *t209*t4+1.8654414521792934d-1*t210*t27+t208)*wght+Amat3
     +        (iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (2.50036235570136d+1*t105*t19*t20*t21*t22*t63*t90*t91
     +        -1.0001449422805439d+2*t105*t19*t20*t21*t51*t61*t90*t91
     +        -1.1668357659939681d+2*t105*t19*t20*t22*t52*t90*t91
     +        +2.50036235570136d+1*t105*t19*t20*t21*t22*t273*t90*t91
     +        -2.0002898845610878d+2*t105*t19*t20*t21*t39*t47*t51*t63
     +        *t91+2.50036235570136d+1*t101*t105*t19*t20*t21*t22*t63
     +        *t91+5.00072471140272d+1*t105*t126*t19*t20*t21*t22*t47
     +        *t91-8.334541185671199d+0*t14*t19*t20*t21*t22*t54*t63
     +        *t90-1.0001449422805439d+2*t19*t20*t21*t22*t221*t251*t63
     +        *t90+3.3338164742684795d+1*t14*t19*t20*t21*t51*t54*t61
     +        *t90+3.889452553313227d+1*t14*t19*t20*t22*t52*t54*t90
     +        +3.3338164742684795d+1*t14*t19*t20*t21*t51*t54*t63*t87
     +        +8.334541185671199d+1*t106*t19*t20*t21*t23*t61*t87
     +        +7.778905106626454d+1*t19*t20*t23*t51*t52*t87
     +        -1.6669082371342397d+2*t106*t14*t19*t20*t21*t54*t63*t76
     +        -5.00072471140272d+2*t19*t20*t21*t23*t252*t61*t76
     +        -3.889452553313226d+2*t106*t19*t20*t23*t52*t76
     +        +6.667632948536959d+1*t14*t19*t20*t21*t39*t47*t51*t54
     +        *t63-8.334541185671199d+0*t101*t14*t19*t20*t21*t22*t54
     +        *t63-2.2040231135441618d+2*t107*t14*t19*t20*t22*t54*t63
     +        +3.3338164742684795d+1*t101*t14*t19*t20*t21*t51*t54*t61
     +        -3.3338164742684795d+2*t106*t14*t19*t20*t21*t39*t47*t54
     +        *t61-4.4080462270883236d+2*t107*t19*t20*t23*t51*t61
     +        -3.1115620426505813d+2*t14*t19*t20*t39*t47*t51*t52*t54
     +        +3.889452553313227d+1*t101*t14*t19*t20*t22*t52*t54
     +        +6.667632948536959d+1*t118*t14*t19*t20*t21*t47*t51*t54
     +        +6.667632948536959d+1*t126*t14*t19*t20*t21*t39*t51*t54
     +        -1.6669082371342397d+1*t126*t14*t19*t20*t21*t22*t47*t54
     +        -4.4080462270883236d+2*t107*t14*t19*t20*t22*t47*t54
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t306*t54
     +        -8.816092454176647d+2*t107*t19*t20*t23*t39*t51
     +        -1.6669082371342397d+1*t19*t20*t21*t23*t261*t51
     +        +1.6669082371342397d+2*t106*t118*t19*t20*t21*t23*t39
     +        +t324+t323+t322+t321+t320+t319+t318+t317+t316+t315+t314
     +        +t313+t312+t311+t310+t309+t308+t307+t300+t299+t298+t297
     +        +t296+t295+t294+t293+t292+t291+t290+t289+t288+t287+t286
     +        +t285+t284+t283+t282+t281+t280+t279+t278+t277+t276+t275
     +        +t274+t272+t271+t253+t236)-1.5877472940427437d+1*t1*t10
     +        *t12*t18*t50*t66*t88+1.5877472940427437d+1*t1*t10*t12
     +        *t25*t65*t88-3.1754945880854873d+1*t1*t10*t18*t25*t40
     +        *t61*t88-3.1754945880854873d+1*t1*t10*t12*t18*t50*t55
     +        *t88+3.1754945880854873d+1*t1*t10*t12*t25*t49*t88
     +        -6.350989176170975d+1*t1*t10*t18*t25*t39*t40*t88
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t66*t87
     +        -4.0827787561099116d+0*t1*t10*t11*t25*t40*t65*t87
     +        +1.2248336268329735d+1*t1*t10*t11*t18*t22*t25*t61*t87
     +        +9.526483764256461d+0*t1*t10*t18*t25*t40*t41*t87
     +        -1.2248336268329735d+1*t1*t10*t11*t18*t22*t50*t66*t76
     +        +1.2248336268329735d+1*t1*t10*t11*t22*t25*t65*t76
     +        -4.899334507331894d+1*t1*t10*t11*t18*t25*t51*t61*t76
     +        -2.8579451292769387d+1*t1*t10*t18*t22*t25*t41*t76
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t49*t55*t66
     +        -1.6331115024439646d+1*t1*t10*t11*t110*t18*t39*t40*t55
     +        *t66+8.165557512219823d+0*t1*t10*t11*t39*t40*t49*t50*t66
     +        -2.0413893780549558d+0*t1*t10*t104*t11*t12*t50*t66
     +        -1.2248336268329735d+1*t1*t10*t11*t111*t12*t18*t254*t66
     +        +4.0827787561099116d+0*t1*t10*t109*t11*t110*t12*t18*t66
     +        +8.165557512219823d+0*t1*t10*t11*t39*t40*t50*t55*t65
     +        -2.0413893780549558d+0*t1*t10*t109*t11*t12*t50*t65
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t111*t12*t65
     +        +8.165557512219823d+0*t1*t10*t11*t40*t49*t50*t55*t61
     +        -2.449667253665947d+1*t1*t10*t11*t18*t22*t39*t50*t55*t61
     +        +4.0827787561099116d+0*t1*t10*t109*t11*t18*t40*t50*t61
     +        +2.449667253665947d+1*t1*t10*t11*t22*t25*t39*t49*t61
     +        -4.0827787561099116d+0*t1*t10*t104*t11*t25*t40*t61
     +        -8.165557512219823d+0*t1*t10*t11*t110*t111*t18*t40*t61
     +        +9.526483764256461d+0*t1*t10*t12*t41*t49*t50*t55
     +        -1.9052967528512924d+1*t1*t10*t18*t39*t40*t41*t50*t55
     +        +8.165557512219823d+0*t1*t10*t11*t118*t18*t40*t50*t55
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t131*t50*t55
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t132*t18*t55
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t132*t49*t50
     +        +4.763241882128231d+0*t1*t10*t109*t12*t18*t41*t50
     +        +8.165557512219823d+0*t1*t10*t11*t132*t18*t39*t40*t50
     +        +1.9052967528512924d+1*t1*t10*t25*t39*t40*t41*t49
     +        -8.165557512219823d+0*t1*t10*t11*t118*t25*t40*t49
     +        -4.763241882128231d+0*t1*t10*t104*t12*t25*t41
     +        -9.526483764256461d+0*t1*t10*t110*t111*t12*t18*t41
     +        -8.165557512219823d+0*t1*t10*t11*t131*t25*t39*t40
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t261*t40
     +        +2.449667253665947d+1*t1*t10*t11*t118*t18*t22*t25*t39
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t317+t316
     +        +t315+t314+t313+t312+t311+t310+t309+t308+t307+t300+t299
     +        +t298+t297+t296+t295+t294+t293+t292+t291+t290+t289+t288
     +        +t287+t286+t285+t284+t283+t282+t281+t280+t279+t278+t277
     +        +t276+t275+t274+t272+t271+t236)+t270+t269+t268+t267+t266
     +        +t265+t264+t263+t262+t260+t259+t219)
     +        +9.327207260896466d-2*t57*t67*t87+1.8654414521792934d-1
     +        *t67*t76*t8+3.7308829043585867d-1*t39*t56*t61*t8
     +        +9.327207260896466d-2*t112*t57*t61+1.8654414521792934d-1
     +        *t118*t56*t57+1.8654414521792934d-1*t133*t39*t57)+t114*(
     +        -6.218138173930978d-2*t133*t56*t9-t113*t67
     +        -9.327207260896466d-2*t115*t57*t61)+6.218138173930978d-2
     +        *t115*t255*t67*t9+1.8654414521792934d-1*t27*t61*t8*t87
     +        +3.7308829043585867d-1*t118*t27*t39*t8
     +        +1.8654414521792934d-1*t27*t61*t76+t28*(t257-1.0d+0*t75)
     +        +(4.0d+0*t30+t256)*t73+9.327207260896466d-2*t261*t27*t57
     +        +1.0d+0*t258*t4+t208)*wght+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -2.0002898845610878d+2*t105*t19*t20*t21*t47*t51*t61*t63
     +        *t91+2.50036235570136d+1*t105*t19*t20*t21*t22*t273*t47
     +        *t63*t91+5.00072471140272d+1*t105*t126*t19*t20*t21*t22
     +        *t63*t91-1.1668357659939681d+2*t105*t142*t19*t20*t22*t52
     +        *t91-1.0001449422805439d+2*t105*t142*t19*t20*t21*t39*t51
     +        *t91+2.50036235570136d+1*t105*t149*t19*t20*t21*t22*t47
     +        *t91+2.50036235570136d+1*t105*t142*t19*t20*t21*t22*t47
     +        *t91-3.1115620426505813d+2*t14*t19*t20*t51*t52*t54*t61
     +        *t63+6.667632948536959d+1*t14*t19*t20*t21*t47*t51*t54
     +        *t61*t63-3.3338164742684795d+2*t106*t14*t19*t20*t21*t39
     +        *t54*t61*t63+6.667632948536959d+1*t118*t14*t19*t20*t21
     +        *t51*t54*t63-1.6669082371342397d+1*t126*t14*t19*t20*t21
     +        *t22*t54*t63-4.4080462270883236d+2*t107*t14*t19*t20*t22
     +        *t54*t63+6.667632948536959d+1*t126*t14*t19*t20*t21*t51
     +        *t54*t61-8.816092454176647d+2*t107*t19*t20*t23*t51*t61
     +        +1.6669082371342397d+2*t106*t118*t19*t20*t21*t23*t61
     +        +3.889452553313227d+1*t14*t149*t19*t20*t22*t52*t54
     +        +3.889452553313227d+1*t14*t142*t19*t20*t22*t52*t54
     +        +3.3338164742684795d+1*t14*t141*t19*t20*t21*t47*t51*t54
     +        +3.3338164742684795d+1*t14*t149*t19*t20*t21*t39*t51*t54
     +        +3.3338164742684795d+1*t14*t142*t19*t20*t21*t39*t51*t54
     +        -8.334541185671199d+0*t14*t149*t19*t20*t21*t22*t47*t54
     +        -8.334541185671199d+0*t14*t142*t19*t20*t21*t22*t47*t54
     +        -2.2040231135441618d+2*t107*t14*t19*t20*t22*t47*t54
     +        -1.6669082371342397d+2*t106*t136*t14*t19*t20*t21*t47*t54
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t347*t54
     +        +7.778905106626454d+1*t141*t19*t20*t23*t51*t52
     +        -3.889452553313226d+2*t106*t136*t19*t20*t23*t52
     +        -4.4080462270883236d+2*t107*t19*t20*t23*t39*t51
     +        -1.6669082371342397d+1*t19*t20*t21*t23*t326*t51
     +        -1.0001449422805439d+2*t142*t19*t20*t21*t22*t221*t251
     +        *t47-5.00072471140272d+2*t136*t19*t20*t21*t23*t252*t39
     +        +8.334541185671199d+1*t106*t141*t19*t20*t21*t23*t39+t360
     +        +t359+t358+t357+t356+t355+t354+t353+t352+t351+t350+t349
     +        +t348+t346+t345+t344+t343+t342+t341+t340+t339+t338+t337
     +        +t336+t335+t334+t333+t332+t331+t330+t329+t328+t327+t324
     +        +t323+t322+t321+t320+t319+t318+t312+t310+t298+t293+t286
     +        +t284+t276+t253+t236)-3.1754945880854873d+1*t1*t10*t12
     +        *t18*t50*t66*t88+3.1754945880854873d+1*t1*t10*t12*t25
     +        *t65*t88-6.350989176170975d+1*t1*t10*t18*t25*t40*t61*t88
     +        -1.5877472940427437d+1*t1*t10*t12*t18*t50*t55*t88
     +        +1.5877472940427437d+1*t1*t10*t12*t25*t49*t88
     +        -3.1754945880854873d+1*t1*t10*t18*t25*t39*t40*t88
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t55*t65*t66
     +        +9.526483764256461d+0*t1*t10*t12*t41*t50*t65*t66
     +        +8.165557512219823d+0*t1*t10*t11*t39*t40*t50*t65*t66
     +        -1.6331115024439646d+1*t1*t10*t11*t110*t18*t40*t55*t61
     +        *t66+8.165557512219823d+0*t1*t10*t11*t40*t49*t50*t61*t66
     +        -1.9052967528512924d+1*t1*t10*t18*t40*t41*t50*t61*t66
     +        -2.449667253665947d+1*t1*t10*t11*t18*t22*t39*t50*t61*t66
     +        +8.165557512219823d+0*t1*t10*t11*t118*t18*t40*t50*t66
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t131*t50*t66
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t132*t18*t66
     +        +8.165557512219823d+0*t1*t10*t11*t40*t50*t55*t61*t65
     +        +1.9052967528512924d+1*t1*t10*t25*t40*t41*t61*t65
     +        +2.449667253665947d+1*t1*t10*t11*t22*t25*t39*t61*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t132*t50*t65
     +        -8.165557512219823d+0*t1*t10*t11*t118*t25*t40*t65
     +        +8.165557512219823d+0*t1*t10*t11*t132*t18*t40*t50*t61
     +        -8.165557512219823d+0*t1*t10*t11*t131*t25*t40*t61
     +        +2.449667253665947d+1*t1*t10*t11*t118*t18*t22*t25*t61
     +        +4.0827787561099116d+0*t1*t10*t11*t141*t18*t40*t50*t55
     +        -1.2248336268329735d+1*t1*t10*t11*t136*t18*t22*t50*t55
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t152*t50*t55
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t154*t18*t254*t55
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t153*t18*t55
     +        -4.899334507331894d+1*t1*t10*t11*t136*t18*t25*t39*t51
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t153*t49*t50
     +        +4.763241882128231d+0*t1*t10*t12*t153*t18*t41*t50
     +        +4.0827787561099116d+0*t1*t10*t11*t153*t18*t39*t40*t50
     +        -4.0827787561099116d+0*t1*t10*t11*t141*t25*t40*t49
     +        +1.2248336268329735d+1*t1*t10*t11*t136*t22*t25*t49
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t154*t49
     +        +9.526483764256461d+0*t1*t10*t141*t18*t25*t40*t41
     +        -2.8579451292769387d+1*t1*t10*t136*t18*t22*t25*t41
     +        -4.763241882128231d+0*t1*t10*t12*t152*t25*t41
     +        -9.526483764256461d+0*t1*t10*t110*t12*t154*t18*t41
     +        -4.0827787561099116d+0*t1*t10*t11*t152*t25*t39*t40
     +        -8.165557512219823d+0*t1*t10*t11*t110*t154*t18*t39*t40
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t326*t40
     +        +1.2248336268329735d+1*t1*t10*t11*t141*t18*t22*t25*t39
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t360+t359
     +        +t358+t357+t356+t355+t354+t353+t352+t351+t350+t349+t348
     +        +t346+t345+t344+t343+t342+t341+t340+t339+t338+t337+t336
     +        +t335+t334+t333+t332+t331+t330+t329+t328+t327+t312+t310
     +        +t298+t293+t286+t284+t276+t236)+t270+t269+t268+t267+t266
     +        +t265+t264+t263+t262+t260+t259+t219)
     +        +3.7308829043585867d-1*t39*t61*t67*t8
     +        +1.8654414521792934d-1*t136*t56*t8+1.8654414521792934d-1
     +        *t118*t57*t67+1.8654414521792934d-1*t133*t57*t61
     +        +9.327207260896466d-2*t141*t56*t57+9.327207260896466d-2
     +        *t155*t39*t57)+t114*(-3.109069086965489d-2*t133*t67*t9
     +        -3.109069086965489d-2*t155*t56*t9-9.327207260896466d-2
     +        *t56*t57*t61*t67-t134*t67)+6.218138173930978d-2*t157
     +        *t255*t56*t9+3.7308829043585867d-1*t118*t27*t61*t8
     +        +1.8654414521792934d-1*t141*t27*t39*t8+(t256-2*t116)*t73
     +        +9.327207260896466d-2*t27*t326*t57+1.0d+0*t325*t4
     +        +1.8654414521792934d-1*t136*t27*t39+(t257-1.0d+0*t135)
     +        *t28+t208)*wght+Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (7.50108706710408d+1*t105*t149*t19*t20*t21*t22*t63*t91
     +        -3.0004348268416314d+2*t105*t142*t19*t20*t21*t51*t61*t91
     +        -3.500507297981904d+2*t105*t142*t19*t20*t22*t52*t91
     +        +2.50036235570136d+1*t105*t19*t20*t21*t22*t364*t91
     +        +2.50036235570136d+1*t105*t142*t19*t20*t21*t22*t273*t91
     +        -9.334686127951745d+2*t14*t19*t20*t51*t52*t54*t61*t63
     +        +1.0001449422805439d+2*t14*t141*t19*t20*t21*t51*t54*t63
     +        -2.50036235570136d+1*t14*t149*t19*t20*t21*t22*t54*t63
     +        -6.612069340632486d+2*t107*t14*t19*t20*t22*t54*t63
     +        -5.00072471140272d+2*t106*t136*t14*t19*t20*t21*t54*t63
     +        +1.0001449422805439d+2*t14*t149*t19*t20*t21*t51*t54*t61
     +        +1.0001449422805439d+2*t14*t142*t19*t20*t21*t51*t54*t61
     +        -1.322413868126497d+3*t107*t19*t20*t23*t51*t61
     +        +2.50036235570136d+2*t106*t141*t19*t20*t21*t23*t61
     +        +1.1668357659939681d+2*t14*t149*t19*t20*t22*t52*t54
     +        +1.1668357659939681d+2*t14*t142*t19*t20*t22*t52*t54
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t384*t54
     +        -8.334541185671199d+0*t14*t19*t20*t21*t22*t364*t54
     +        +2.3336715319879361d+2*t141*t19*t20*t23*t51*t52
     +        -1.166835765993968d+3*t106*t136*t19*t20*t23*t52
     +        -1.6669082371342397d+1*t19*t20*t21*t23*t363*t51+t388
     +        +t387+t386+t385+t383+t382+t381+t380+t379+t378+t377+t376
     +        +t375+t374+t373+t372+t371+t370+t369+t368+t367+t366+t365
     +        -1.0001449422805439d+2*t19*t20*t21*t22*t221*t251*t364
     +        -5.00072471140272d+2*t19*t20*t21*t23*t252*t362+t253+t236
     +        )-4.763241882128231d+1*t1*t10*t12*t18*t50*t66*t88
     +        +4.763241882128231d+1*t1*t10*t12*t25*t65*t88
     +        -9.526483764256461d+1*t1*t10*t18*t25*t40*t61*t88
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t18*t254*t66**3
     +        +2.449667253665947d+1*t1*t10*t11*t40*t50*t61*t65*t66
     +        +2.8579451292769387d+1*t1*t10*t12*t41*t50*t65*t66
     +        -5.715890258553877d+1*t1*t10*t18*t40*t41*t50*t61*t66
     +        +1.2248336268329735d+1*t1*t10*t11*t141*t18*t40*t50*t66
     +        -3.6745008804989204d+1*t1*t10*t11*t136*t18*t22*t50*t66
     +        -6.124168134164868d+0*t1*t10*t11*t12*t152*t50*t66
     +        +1.2248336268329735d+1*t1*t10*t11*t110*t12*t153*t18*t66
     +        +5.715890258553877d+1*t1*t10*t25*t40*t41*t61*t65
     +        -6.124168134164868d+0*t1*t10*t11*t12*t153*t50*t65
     +        -1.2248336268329735d+1*t1*t10*t11*t141*t25*t40*t65
     +        +3.6745008804989204d+1*t1*t10*t11*t136*t22*t25*t65
     +        +1.2248336268329735d+1*t1*t10*t11*t110*t12*t154*t65
     +        +1.2248336268329735d+1*t1*t10*t11*t153*t18*t40*t50*t61
     +        -1.2248336268329735d+1*t1*t10*t11*t152*t25*t40*t61
     +        -2.449667253665947d+1*t1*t10*t11*t110*t154*t18*t40*t61
     +        +3.6745008804989204d+1*t1*t10*t11*t141*t18*t22*t25*t61
     +        -4.899334507331894d+1*t1*t10*t11*t18*t25*t362*t51
     +        +1.4289725646384693d+1*t1*t10*t12*t153*t18*t41*t50
     +        +2.8579451292769387d+1*t1*t10*t141*t18*t25*t40*t41
     +        -8.573835387830815d+1*t1*t10*t136*t18*t22*t25*t41
     +        -1.4289725646384693d+1*t1*t10*t12*t152*t25*t41
     +        -2.8579451292769387d+1*t1*t10*t110*t12*t154*t18*t41
     +        -4.0827787561099116d+0*t1*t10*t11*t18*t25*t363*t40
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t388+t387
     +        +t386+t385+t383+t382+t381+t380+t379+t378+t377+t376+t375
     +        +t374+t373+t372+t371+t370+t369+t368+t367+t366+t365+t236)
     +        +t219)+5.596324356537881d-1*t136*t67*t8
     +        +2.7981621782689403d-1*t141*t57*t67
     +        +2.7981621782689403d-1*t155*t57*t61)+t114*(
     +        -6.218138173930978d-2*t155*t67*t9-t156*t67
     +        -9.327207260896466d-2*t157*t57*t61)+6.218138173930978d-2
     +        *t255*t67**3*t9+5.596324356537881d-1*t141*t27*t61*t8
     +        +6.0d+0*t58*t73+9.327207260896466d-2*t27*t363*t57+1.0d+0
     +        *t361*t4+1.8654414521792934d-1*t27*t362-3.0d+0*t135*t28
     +        +t208)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t410+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (1.0001449422805439d+2*t10*t105*t19*t21*t22*t90*t91
     +        -3.3338164742684795d+1*t10*t14*t19*t21*t22*t54*t90
     +        -6.667632948536959d+1*t10*t19*t21*t23*t51*t87
     +        +3.3338164742684795d+2*t10*t106*t19*t21*t23*t76
     +        +3.1115620426505813d+2*t10*t14*t19*t22*t47*t52*t54
     +        +2.667053179414784d+2*t10*t14*t19*t21*t39*t47*t51*t54
     +        -3.3338164742684795d+1*t10*t101*t14*t19*t21*t22*t54
     +        +6.223124085301163d+2*t10*t19*t23*t39*t51*t52+t423+t422
     +        +t421+t420+t419+t418+t417+t416+t415+t414)
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t40*t50*t71*t87
     +        -1.6669082371342397d+1*t10*t16*t19*t21*t25*t51*t87
     +        -8.165557512219823d+0*t1*t11*t18*t25*t40*t87
     +        -1.2248336268329735d+1*t1*t10*t11*t18*t22*t50*t71*t76
     +        +2.449667253665947d+1*t1*t11*t18*t22*t25*t76
     +        +5.00072471140272d+1*t10*t106*t16*t19*t21*t25*t76
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t49*t55*t71
     +        -1.9052967528512924d+1*t1*t10*t110*t12*t18*t41*t55*t71
     +        -1.6331115024439646d+1*t1*t10*t11*t110*t18*t39*t40*t55
     +        *t71+9.526483764256461d+0*t1*t10*t12*t41*t49*t50*t71
     +        +8.165557512219823d+0*t1*t10*t11*t39*t40*t49*t50*t71
     +        -1.9052967528512924d+1*t1*t10*t18*t39*t40*t41*t50*t71
     +        -2.0413893780549558d+0*t1*t10*t104*t11*t12*t50*t71
     +        -1.2248336268329735d+1*t1*t10*t11*t111*t12*t18*t254*t71
     +        +4.0827787561099116d+0*t1*t10*t109*t11*t110*t12*t18*t71
     +        +3.889452553313227d+1*t10*t16*t19*t22*t50*t52*t55
     +        +3.3338164742684795d+1*t10*t16*t19*t21*t39*t50*t51*t55
     +        -8.165557512219823d+0*t1*t11*t12*t49*t50*t55
     +        +1.9052967528512924d+1*t1*t12*t18*t41*t50*t55
     +        +1.6331115024439646d+1*t1*t11*t18*t39*t40*t50*t55
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t176*t50*t55
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t178*t18*t55
     +        +7.778905106626454d+1*t10*t16*t19*t25*t39*t51*t52
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t178*t49*t50
     +        +9.526483764256461d+0*t1*t10*t12*t178*t18*t41*t50
     +        +8.165557512219823d+0*t1*t10*t11*t178*t18*t39*t40*t50
     +        -8.334541185671199d+0*t10*t109*t16*t19*t21*t22*t50
     +        -4.0827787561099116d+0*t1*t109*t11*t12*t18*t50
     +        -1.9052967528512924d+1*t1*t12*t25*t41*t49
     +        -1.6331115024439646d+1*t1*t11*t25*t39*t40*t49
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t422+t421
     +        +t420+t419+t418+t417+t416+t415+t414)+t413+t412+t411
     +        +3.8105935057025847d+1*t1*t18*t25*t39*t40*t41
     +        -9.526483764256461d+0*t1*t10*t12*t176*t25*t41
     +        -8.165557512219823d+0*t1*t10*t11*t176*t25*t39*t40
     +        +4.0827787561099116d+0*t1*t104*t11*t12*t25
     +        +1.6669082371342397d+1*t10*t110*t111*t16*t19*t21*t22
     +        +8.165557512219823d+0*t1*t11*t110*t111*t12*t18)
     +        +9.327207260896466d-2*t57*t72*t87+1.8654414521792934d-1
     +        *t72*t76*t8+1.8654414521792934d-1*t179*t39*t57)+t114*(
     +        -6.218138173930978d-2*t179*t56*t9-t113*t72)
     +        +6.218138173930978d-2*t115*t255*t72*t9)*wght+Cmat3(iq,
     +        D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t410+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = t435+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (1.0001449422805439d+2*t10*t105*t19*t21*t22*t47*t63*t91
     +        +1.333526589707392d+2*t10*t14*t19*t21*t39*t51*t54*t63
     +        -3.3338164742684795d+1*t10*t14*t19*t21*t22*t47*t54*t63
     +        +1.333526589707392d+2*t10*t14*t19*t21*t47*t51*t54*t61
     +        +3.3338164742684795d+2*t10*t106*t19*t21*t23*t39*t61
     +        -3.3338164742684795d+1*t10*t126*t14*t19*t21*t22*t54
     +        -6.667632948536959d+1*t10*t118*t19*t21*t23*t51+t448+t447
     +        +t446+t445+t444+t443+t442+t441+t440+t439+t438+t423+t418
     +        +t409+t407+t404+t401)-1.2248336268329735d+1*t1*t10*t11
     +        *t12*t18*t254*t55*t66*t71+4.0827787561099116d+0*t1*t10
     +        *t11*t110*t12*t49*t66*t71-9.526483764256461d+0*t1*t10
     +        *t110*t12*t18*t41*t66*t71-8.165557512219823d+0*t1*t10
     +        *t11*t110*t18*t39*t40*t66*t71+4.0827787561099116d+0*t1
     +        *t10*t11*t110*t12*t55*t65*t71+4.763241882128231d+0*t1
     +        *t10*t12*t41*t50*t65*t71+4.0827787561099116d+0*t1*t10
     +        *t11*t39*t40*t50*t65*t71-8.165557512219823d+0*t1*t10*t11
     +        *t110*t18*t40*t55*t61*t71+4.0827787561099116d+0*t1*t10
     +        *t11*t40*t49*t50*t61*t71-9.526483764256461d+0*t1*t10*t18
     +        *t40*t41*t50*t61*t71-1.2248336268329735d+1*t1*t10*t11
     +        *t18*t22*t39*t50*t61*t71-9.526483764256461d+0*t1*t10
     +        *t110*t12*t18*t41*t55*t71+4.763241882128231d+0*t1*t10
     +        *t12*t41*t49*t50*t71-9.526483764256461d+0*t1*t10*t18*t39
     +        *t40*t41*t50*t71+4.0827787561099116d+0*t1*t10*t11*t118
     +        *t18*t40*t50*t71-2.0413893780549558d+0*t1*t10*t11*t12
     +        *t131*t50*t71+4.0827787561099116d+0*t1*t10*t11*t110*t12
     +        *t132*t18*t71+1.6669082371342397d+1*t10*t110*t16*t19*t21
     +        *t22*t55*t66+8.165557512219823d+0*t1*t11*t110*t12*t18
     +        *t55*t66+1.6669082371342397d+1*t10*t16*t19*t21*t39*t50
     +        *t51*t66-4.0827787561099116d+0*t1*t11*t12*t49*t50*t66
     +        +8.165557512219823d+0*t1*t11*t18*t39*t40*t50*t66
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t176*t50*t66
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t178*t18*t66
     +        -4.0827787561099116d+0*t1*t11*t12*t50*t55*t65
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t178*t50*t65
     +        -8.165557512219823d+0*t1*t11*t25*t39*t40*t65
     +        +1.6669082371342397d+1*t10*t16*t19*t21*t50*t51*t55*t61
     +        +8.165557512219823d+0*t1*t11*t18*t40*t50*t55*t61
     +        +4.0827787561099116d+0*t1*t10*t11*t178*t18*t40*t50*t61
     +        -8.165557512219823d+0*t1*t11*t25*t40*t49*t61
     +        -4.0827787561099116d+0*t1*t10*t11*t176*t25*t40*t61
     +        +2.449667253665947d+1*t1*t11*t18*t22*t25*t39*t61
     +        +5.00072471140272d+1*t10*t106*t16*t19*t21*t25*t39*t61
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t190*t50*t55
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t18*t191*t55
     +        -1.6669082371342397d+1*t10*t118*t16*t19*t21*t25*t51
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t191*t49*t50
     +        +4.763241882128231d+0*t1*t10*t12*t18*t191*t41*t50
     +        +4.763241882128231d+0*t1*t10*t12*t178*t18*t41*t50
     +        +4.0827787561099116d+0*t1*t10*t11*t18*t191*t39*t40*t50
     +        -8.334541185671199d+0*t10*t132*t16*t19*t21*t22*t50
     +        -4.0827787561099116d+0*t1*t11*t12*t132*t18*t50+t451+t450
     +        +t449+2.0413893780549558d+0*t1*t10*t11*t12*t25*(t446
     +        +t445+t444+t443+t442+t441+t440+t439+t438+t418+t404+t401)
     +        +t437+t436+t413+t412+t411-4.763241882128231d+0*t1*t10
     +        *t12*t190*t25*t41-4.763241882128231d+0*t1*t10*t12*t176
     +        *t25*t41-4.0827787561099116d+0*t1*t10*t11*t190*t25*t39
     +        *t40-8.165557512219823d+0*t1*t11*t118*t18*t25*t40+t397
     +        +t396+t395+t392+t390+4.0827787561099116d+0*t1*t11*t12
     +        *t131*t25)+1.8654414521792934d-1*t39*t61*t72*t8
     +        +9.327207260896466d-2*t118*t57*t72+9.327207260896466d-2
     +        *t179*t57*t61+9.327207260896466d-2*t192*t39*t57)+t114*(
     +        -3.109069086965489d-2*t179*t67*t9-3.109069086965489d-2
     +        *t192*t56*t9-t134*t72)+6.218138173930978d-2*t255*t56*t67
     +        *t72*t9)*wght+Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = t435+Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = t458+Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = (t31*(3.109069086965489d-2*t9*(
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*
     +        (1.0001449422805439d+2*t10*t105*t142*t19*t21*t22*t91
     +        +2.667053179414784d+2*t10*t14*t19*t21*t51*t54*t61*t63
     +        +3.1115620426505813d+2*t10*t14*t19*t22*t52*t54*t63
     +        +6.223124085301163d+2*t10*t19*t23*t51*t52*t61
     +        -3.3338164742684795d+1*t10*t14*t149*t19*t21*t22*t54
     +        -3.3338164742684795d+1*t10*t14*t142*t19*t21*t22*t54
     +        -6.667632948536959d+1*t10*t141*t19*t21*t23*t51+t466+t465
     +        +t464+t463+t462+t461+t460+t459+t423+t418
     +        +3.3338164742684795d+2*t10*t106*t136*t19*t21*t23)
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t65*t66*t71
     +        -1.6331115024439646d+1*t1*t10*t11*t110*t18*t40*t61*t66
     +        *t71-1.9052967528512924d+1*t1*t10*t110*t12*t18*t41*t66
     +        *t71+8.165557512219823d+0*t1*t10*t11*t40*t50*t61*t65*t71
     +        +9.526483764256461d+0*t1*t10*t12*t41*t50*t65*t71
     +        -1.9052967528512924d+1*t1*t10*t18*t40*t41*t50*t61*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t141*t18*t40*t50*t71
     +        -1.2248336268329735d+1*t1*t10*t11*t136*t18*t22*t50*t71
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t152*t50*t71
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t154*t18*t254*t71
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t153*t18*t71
     +        -8.165557512219823d+0*t1*t11*t12*t50*t65*t66
     +        +3.3338164742684795d+1*t10*t16*t19*t21*t50*t51*t61*t66
     +        +1.6331115024439646d+1*t1*t11*t18*t40*t50*t61*t66
     +        +3.889452553313227d+1*t10*t16*t19*t22*t50*t52*t66
     +        +1.9052967528512924d+1*t1*t12*t18*t41*t50*t66
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t190*t50*t66
     +        +8.165557512219823d+0*t1*t10*t11*t110*t12*t18*t191*t66
     +        -1.6331115024439646d+1*t1*t11*t25*t40*t61*t65
     +        -4.0827787561099116d+0*t1*t10*t11*t12*t191*t50*t65
     +        -1.9052967528512924d+1*t1*t12*t25*t41*t65
     +        +7.778905106626454d+1*t10*t16*t19*t25*t51*t52*t61
     +        +8.165557512219823d+0*t1*t10*t11*t18*t191*t40*t50*t61
     +        +3.8105935057025847d+1*t1*t18*t25*t40*t41*t61
     +        -8.165557512219823d+0*t1*t10*t11*t190*t25*t40*t61
     +        -1.6669082371342397d+1*t10*t141*t16*t19*t21*t25*t51
     +        +9.526483764256461d+0*t1*t10*t12*t18*t191*t41*t50
     +        -8.334541185671199d+0*t10*t153*t16*t19*t21*t22*t50
     +        -4.0827787561099116d+0*t1*t11*t12*t153*t18*t50
     +        +2.0413893780549558d+0*t1*t10*t11*t12*t25*(t466+t465
     +        +t464+t463+t462+t461+t460+t459+t418)+t413+t412+t411
     +        -9.526483764256461d+0*t1*t10*t12*t190*t25*t41
     +        -8.165557512219823d+0*t1*t11*t141*t18*t25*t40
     +        +2.449667253665947d+1*t1*t11*t136*t18*t22*t25
     +        +5.00072471140272d+1*t10*t106*t136*t16*t19*t21*t25
     +        +4.0827787561099116d+0*t1*t11*t12*t152*t25
     +        +1.6669082371342397d+1*t10*t110*t154*t16*t19*t21*t22
     +        +8.165557512219823d+0*t1*t11*t110*t12*t154*t18)
     +        +1.8654414521792934d-1*t136*t72*t8+9.327207260896466d-2
     +        *t141*t57*t72+1.8654414521792934d-1*t192*t57*t61)+t114*(
     +        -6.218138173930978d-2*t192*t67*t9-t156*t72)
     +        +6.218138173930978d-2*t157*t255*t72*t9)*wght+Cmat3(iq,
     +        D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = t458+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = t482+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = (t114*(-3.109069086965489d-2
     +        *t179*t69*t9-t168*t72+t491)+t490+t483)*wght+Cmat3(iq,
     +        D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t482+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = (t114*(-3.109069086965489d-2
     +        *t179*t72*t9-3.109069086965489d-2*t206*t56*t9-t180*t72)
     +        +t31*(3.109069086965489d-2*(3.3338164742684795d+1*t10
     +        *t110*t16*t19*t21*t22*t55*t71+1.6331115024439646d+1*t1
     +        *t11*t110*t12*t18*t55*t71+3.3338164742684795d+1*t10*t16
     +        *t19*t21*t39*t50*t51*t71-8.165557512219823d+0*t1*t11*t12
     +        *t49*t50*t71+1.6331115024439646d+1*t1*t11*t18*t39*t40
     +        *t50*t71-4.0827787561099116d+0*t1*t10*t11*t12*t176*t50
     +        *t71+8.165557512219823d+0*t1*t10*t11*t110*t12*t178*t18
     +        *t71-3.3338164742684795d+1*t16*t19*t21*t22*t50*t55
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t18*t205*t254*t55
     +        +1.3611235077912595d+2*t10*t106*t110*t18*t197*t198*t23
     +        *t55-2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -6.667632948536959d+1*t14*t19*t21*t22*t47*t54
     +        -1.333526589707392d+2*t19*t21*t23*t39*t51+t497)
     +        -6.667632948536959d+1*t16*t19*t21*t25*t39*t51
     +        -6.805617538956298d+1*t10*t106*t197*t198*t23*t49*t50
     +        +1.3611235077912595d+2*t10*t18*t197*t198*t23*t252*t39
     +        *t50-1.6669082371342397d+1*t10*t16*t178*t19*t21*t22*t50
     +        -8.165557512219823d+0*t1*t11*t12*t178*t18*t50+t496+t495
     +        +t494+t493+t492+4.0827787561099116d+0*t1*t10*t11*t110
     +        *t12*t205*t49-8.165557512219823d+0*t1*t10*t11*t110*t18
     +        *t205*t39*t40+8.165557512219823d+0*t1*t11*t12*t176*t25)
     +        *t9+9.327207260896466d-2*t206*t39*t57)
     +        +6.218138173930978d-2*t204*t255*t56*t9)*wght+Cmat3(iq,
     +        D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = (t114*(-3.109069086965489d-2
     +        *t167*t72*t9-t180*t69+t491)+t490+t483)*wght+Cmat3(iq,
     +        D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t482+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = t506+Cmat3(iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = (t114*(-3.109069086965489d-2
     +        *t192*t69*t9-t186*t72+t509)+t508+t507)*wght+Cmat3(iq,
     +        D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = t506+Cmat3(iq,D3_RB_GAA_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = (t114*(-3.109069086965489d-2
     +        *t192*t72*t9-3.109069086965489d-2*t206*t67*t9-t193*t72)
     +        +t31*(3.109069086965489d-2*(3.3338164742684795d+1*t10
     +        *t110*t16*t19*t21*t22*t66*t71+1.6331115024439646d+1*t1
     +        *t11*t110*t12*t18*t66*t71-8.165557512219823d+0*t1*t11
     +        *t12*t50*t65*t71+3.3338164742684795d+1*t10*t16*t19*t21
     +        *t50*t51*t61*t71+1.6331115024439646d+1*t1*t11*t18*t40
     +        *t50*t61*t71-4.0827787561099116d+0*t1*t10*t11*t12*t190
     +        *t50*t71+8.165557512219823d+0*t1*t10*t11*t110*t12*t18
     +        *t191*t71-3.3338164742684795d+1*t16*t19*t21*t22*t50*t66
     +        -1.2248336268329735d+1*t1*t10*t11*t12*t18*t205*t254*t66
     +        +1.3611235077912595d+2*t10*t106*t110*t18*t197*t198*t23
     +        *t66-6.805617538956298d+1*t10*t106*t197*t198*t23*t50*t65
     +        +4.0827787561099116d+0*t1*t10*t11*t110*t12*t205*t65
     +        -2.0413893780549558d+0*t1*t10*t11*t12*t18*t50*(
     +        -6.667632948536959d+1*t14*t19*t21*t22*t54*t63
     +        -1.333526589707392d+2*t19*t21*t23*t51*t61+t497)
     +        -6.667632948536959d+1*t16*t19*t21*t25*t51*t61
     +        +1.3611235077912595d+2*t10*t18*t197*t198*t23*t252*t50
     +        *t61-8.165557512219823d+0*t1*t10*t11*t110*t18*t205*t40
     +        *t61-1.6669082371342397d+1*t10*t16*t19*t191*t21*t22*t50
     +        -8.165557512219823d+0*t1*t11*t12*t18*t191*t50+t496+t495
     +        +t494+t493+t492+8.165557512219823d+0*t1*t11*t12*t190*t25
     +        )*t9+9.327207260896466d-2*t206*t57*t61)
     +        +6.218138173930978d-2*t204*t255*t67*t9)*wght+Cmat3(iq,
     +        D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = (t114*(-3.109069086965489d-2
     +        *t185*t72*t9-t193*t69+t509)+t508+t507)*wght+Cmat3(iq,
     +        D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t506+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t513+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = 3.109069086965489d-2*t31*(
     +        -2.0839373032694847d+2*t10*t50*t510*t511*t512*t54+t523
     +        +t522+t521+t520+t519+t518+t517+t516+t515+t514)*t9*wght
     +        +t525+t524+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t513+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = 3.109069086965489d-2*t31*(
     +        -4.167874606538969d+2*t10*t50*t510*t511*t512*t54+t534
     +        +t533+t532+t531+t530+t529+t528+t527+t522+t517)*t9*wght
     +        +t535+t526+Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = t536+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = t513+Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = 6.218138173930978d-2*t255*t72*
     +        *3*t9*wght-9.327207260896466d-2*t114*t206*t72*t9*wght
     +        +3.109069086965489d-2*t31*(-1.2248336268329735d+1*t1*t10
     +        *t11*t12*t18*t254*t71**3-1.0001449422805439d+2*t16*t19
     +        *t21*t22*t50*t71+4.083370523373779d+2*t10*t106*t110*t18
     +        *t197*t198*t23*t71-8.335749213077939d+2*t10*t50*t510
     +        *t511*t512*t54-4.083370523373779d+2*t106*t18*t197*t198
     +        *t23*t50+5.00072471140272d+1*t10*t110*t16*t19*t205*t21
     +        *t22+2.449667253665947d+1*t1*t11*t110*t12*t18*t205)*t9
     +        *wght+Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = 3.109069086965489d-2*t31*(
     +        -4.1678746065389694d+2*t10*t50*t510*t511*t512*t54+t534
     +        +t533+t532+t531+t530+t529+t528+t527+t522+t517)*t9*wght
     +        +t535+t526+Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = t536+Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t513+Cmat3(iq,D3_GBB_GBB_GBB)
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
            t3 = 1/rhoa
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = exp(-6.432793688582967d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 3.240503646201146d+0*gammaaa*t1*t4*t7
            t9 = t8+1.0d+0
            t10 = t1**2
            t11 = gammaaa**2
            t12 = 1/rhoa**4.666666666666667d+0
            t13 = 1/t6**2
            t14 = t8+1.0500863881042923d+1*t10*t11*t12*t13+1.0d+0
            t15 = 1/t14
            t16 = 3.240503646201146d+0*gammaaa*t1*t15*t4*t9+1.0d+0
            t17 = 1/rhoa**2
            t18 = s1a(D1_RA)
            t19 = 1/t16
            t20 = 1/rhoa**3.3333333333333337d+0
            t21 = -7.561175174469341d+0*gammaaa*t1*t20*t7
            t22 = 6.432793688582967d+1*t17*s1f-6.432793688582967d+1
     +        *t18*t3
            t23 = -3.240503646201146d+0*gammaaa*t1*t13*t22*t4*t5
            t24 = t23+t21
            t25 = 1/t14**2
            t26 = 1/rhoa**5.666666666666667d+0
            t27 = 1/t6**3
            t28 = -2.1001727762085845d+1*t10*t11*t12*t22*t27*t5
     +        -4.9004031444866974d+1*t10*t11*t13*t26+t23+t21
            t29 = -3.240503646201146d+0*gammaaa*t1*t25*t28*t4*t9
     +        -7.561175174469341d+0*gammaaa*t1*t15*t20*t9
     +        +3.240503646201146d+0*gammaaa*t1*t15*t24*t4
            t30 = 3.240503646201146d+0*t1*t4*t7+2.1001727762085845d+1
     +        *gammaaa*t10*t12*t13
            t31 = -3.240503646201146d+0*gammaaa*t1*t25*t30*t4*t9
     +        +3.240503646201146d+0*t1*t15*t4*t9+1.0500863881042925d+1
     +        *gammaaa*t10*t12*t15*t7
            t32 = 1/rhoa**3
            t33 = s1a2(D2_RA_RA)
            t34 = 1/rhoa**4.333333333333334d+0
            t35 = t22**2
            t36 = exp(-1.2865587377165935d+2*s1f*t3)
            t37 = 6.481007292402292d+0*gammaaa*t1*t27*t35*t36*t4
            t38 = 2.520391724823114d+1*gammaaa*t1*t34*t7
            t39 = 1.5122350348938685d+1*gammaaa*t1*t13*t20*t22*t5
            t40 = -3.240503646201146d+0*gammaaa*t1*t13*t35*t4*t5
            t41 = -6.432793688582967d+1*t3*t33-1.2865587377165935d+2
     +        *s1f*t32+1.2865587377165935d+2*t17*t18
            t42 = -3.240503646201146d+0*gammaaa*t1*t13*t4*t41*t5
            t43 = t42+t40+t39+t38+t37
            t44 = 1/t6**4
            t45 = 1/rhoa**6.666666666666667d+0
            t46 = -2.1001727762085845d+1*t10*t11*t12*t27*t41*t5
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t35*t5
     +        +1.960161257794679d+2*t10*t11*t22*t26*t27*t5
     +        +2.7768951152091287d+2*t10*t11*t13*t45
     +        +6.300518328625754d+1*t10*t11*t12*t35*t36*t44+t42+t40
     +        +t39+t38+t37
            t47 = 1/t14**3
            t48 = t28**2
            t49 = 6.481007292402292d+0*gammaaa*t1*t4*t47*t48*t9
     +        -3.240503646201146d+0*gammaaa*t1*t25*t4*t46*t9
     +        +2.520391724823114d+1*gammaaa*t1*t15*t34*t9
     +        +1.5122350348938685d+1*gammaaa*t1*t20*t25*t28*t9
     +        +3.240503646201146d+0*gammaaa*t1*t15*t4*t43
     +        -6.481007292402292d+0*gammaaa*t1*t24*t25*t28*t4
     +        -1.5122350348938685d+1*gammaaa*t1*t15*t20*t24
            t50 = 1/t16**2
            t51 = t29**2
            t52 = -7.561175174469341d+0*t1*t20*t7
            t53 = -3.240503646201146d+0*t1*t13*t22*t4*t5
            t54 = t53+t52
            t55 = t53+t52-4.200345552417169d+1*gammaaa*t10*t12*t22*t27
     +        *t5-9.800806288973395d+1*gammaaa*t10*t13*t26
            t56 = -3.240503646201146d+0*gammaaa*t1*t25*t4*t55*t9
     +        +6.481007292402292d+0*gammaaa*t1*t28*t30*t4*t47*t9
     +        -3.240503646201146d+0*t1*t25*t28*t4*t9
     +        +7.561175174469341d+0*gammaaa*t1*t20*t25*t30*t9
     +        -7.561175174469341d+0*t1*t15*t20*t9
     +        -1.0500863881042925d+1*gammaaa*t10*t12*t25*t28*t7
     +        -2.450201572243349d+1*gammaaa*t10*t15*t26*t7
     +        +3.240503646201146d+0*gammaaa*t1*t15*t4*t54
     +        -3.240503646201146d+0*gammaaa*t1*t24*t25*t30*t4
     +        +3.240503646201146d+0*t1*t15*t24*t4
            t57 = t31**2
            t58 = t30**2
            t59 = t1**3
            t60 = 1/rhoa**7.0d+0
            t61 = -6.805617538956302d+1*gammaaa*t13*t25*t59*t60*t9
     +        +6.481007292402292d+0*gammaaa*t1*t4*t47*t58*t9
     +        -6.481007292402292d+0*t1*t25*t30*t4*t9
     +        -2.100172776208585d+1*gammaaa*t10*t12*t25*t30*t7
     +        +2.100172776208585d+1*t10*t12*t15*t7
            t62 = 1/rhoa**4
            t63 = s1a3(D3_RA_RA_RA)
            t64 = 1/rhoa**5.333333333333334d+0
            t65 = t22**3
            t66 = exp(-1.92983810657489d+2*s1f*t3)
            t67 = -1.9443021877206879d+1*gammaaa*t1*t4*t44*t65*t66
            t68 = -4.5367051046816054d+1*gammaaa*t1*t20*t27*t35*t36
            t69 = 1.2865587377165935d+2*t17*s1f-1.2865587377165935d+2
     +        *t18*t3
            t70 = 6.481007292402292d+0*gammaaa*t1*t27*t35*t36*t4*t69
            t71 = 6.481007292402292d+0*gammaaa*t1*t27*t36*t4*t65
            t72 = 1.9443021877206879d+1*gammaaa*t1*t22*t27*t36*t4*t41
            t73 = -1.0921697474233495d+2*gammaaa*t1*t64*t7
            t74 = -7.561175174469342d+1*gammaaa*t1*t13*t22*t34*t5
            t75 = 2.2683525523408027d+1*gammaaa*t1*t13*t20*t35*t5
            t76 = 2.2683525523408027d+1*gammaaa*t1*t13*t20*t41*t5
            t77 = -3.240503646201146d+0*gammaaa*t1*t13*t4*t5*t65
            t78 = -9.721510938603439d+0*gammaaa*t1*t13*t22*t4*t41*t5
            t79 = -6.432793688582967d+1*t3*t63+3.85967621314978d+2*s1f
     +        *t62+1.92983810657489d+2*t17*t33-3.85967621314978d+2*t18
     +        *t32
            t80 = -3.240503646201146d+0*gammaaa*t1*t13*t4*t5*t79
            t81 = 1/t14**4
            t82 = 1/t16**3
            t83 = 6.481007292402292d+0*t1*t27*t35*t36*t4
            t84 = 2.520391724823114d+1*t1*t34*t7
            t85 = 1.5122350348938685d+1*t1*t13*t20*t22*t5
            t86 = -3.240503646201146d+0*t1*t13*t35*t4*t5
            t87 = -3.240503646201146d+0*t1*t13*t4*t41*t5
            fnc(iq) = (1.0d+0*s1f*t3+1.5545345434827443d-2*log(t16))
     +        *wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t18*t3+1.5545345434827443d-2*t19
     +        *t29-1.0d+0*t17*s1f)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.5545345434827443d-2*t19*t31*wght+Cmat
     +        (iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = (-1.5545345434827443d-2*t50*t51
     +        +1.5545345434827443d-2*t19*t49+1.0d+0*t3*t33+2.0d+0*s1f
     +        *t32-2.0d+0*t17*t18)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = (1.5545345434827443d-2*t19*t56
     +        -1.5545345434827443d-2*t29*t31*t50)*wght+Cmat2(iq,
     +        D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.5545345434827443d-2*t19*t61*wght
     +        -1.5545345434827443d-2*t50*t57*wght+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = (1.5545345434827443d-2*t19*(
     +        -1.9443021877206879d+1*gammaaa*t1*t28**3*t4*t81*t9
     +        -3.240503646201146d+0*gammaaa*t1*t25*t4*(t80
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t5*t79+t78+t77
     +        +t76+t75+t74+t73+t72+t71+t70+6.300518328625754d+1*t10
     +        *t11*t12*t35*t36*t44*t69+t68+t67-2.5202073314503015d+2
     +        *t10*t11*t12*t65*t66/t6**5-2.1001727762085845d+1*t10*t11
     +        *t12*t27*t5*t65+6.300518328625754d+1*t10*t11*t12*t36*t44
     +        *t65-1.666137069125477d+3*t10*t11*t22*t27*t45*t5
     +        +2.940241886692019d+2*t10*t11*t26*t27*t41*t5
     +        -6.300518328625754d+1*t10*t11*t12*t22*t27*t41*t5
     +        +2.940241886692019d+2*t10*t11*t26*t27*t35*t5
     +        +1.890155498587726d+2*t10*t11*t12*t22*t36*t41*t44
     +        -8.820725660076056d+2*t10*t11*t26*t35*t36*t44
     +        -1.851263410139419d+3*t10*t11*t13/rhoa*
     +        *7.666666666666667d+0)*t9-1.0921697474233495d+2*gammaaa
     +        *t1*t15*t64*t9-4.5367051046816054d+1*gammaaa*t1*t20*t47
     +        *t48*t9+1.9443021877206879d+1*gammaaa*t1*t28*t4*t46*t47
     +        *t9+2.2683525523408027d+1*gammaaa*t1*t20*t25*t46*t9
     +        -7.561175174469342d+1*gammaaa*t1*t25*t28*t34*t9
     +        +3.240503646201146d+0*gammaaa*t1*t15*t4*(t80+t78+t77+t76
     +        +t75+t74+t73+t72+t71+t70+t68+t67)+1.9443021877206879d+1
     +        *gammaaa*t1*t24*t4*t47*t48-9.721510938603439d+0*gammaaa
     +        *t1*t24*t25*t4*t46-9.721510938603439d+0*gammaaa*t1*t25
     +        *t28*t4*t43-2.2683525523408027d+1*gammaaa*t1*t15*t20*t43
     +        +7.561175174469342d+1*gammaaa*t1*t15*t24*t34
     +        +4.5367051046816054d+1*gammaaa*t1*t20*t24*t25*t28)
     +        +3.1090690869654886d-2*t29**3*t82+1.0d+0*t3*t63-6.0d+0
     +        *s1f*t62-4.6636036304482326d-2*t29*t49*t50-3.0d+0*t17
     +        *t33+6.0d+0*t18*t32)*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammaaa*t1*t25*t4*(t87+t86+t85+t84
     +        +t83-4.200345552417169d+1*gammaaa*t10*t12*t27*t41*t5
     +        -4.200345552417169d+1*gammaaa*t10*t12*t27*t35*t5
     +        +3.920322515589358d+2*gammaaa*t10*t22*t26*t27*t5
     +        +5.553790230418257d+2*gammaaa*t10*t13*t45
     +        +1.2601036657251508d+2*gammaaa*t10*t12*t35*t36*t44)*t9
     +        -1.9443021877206879d+1*gammaaa*t1*t30*t4*t48*t81*t9
     +        +1.2962014584804585d+1*gammaaa*t1*t28*t4*t47*t55*t9
     +        +1.5122350348938685d+1*gammaaa*t1*t20*t25*t55*t9
     +        +6.481007292402292d+0*t1*t4*t47*t48*t9
     +        +6.481007292402292d+0*gammaaa*t1*t30*t4*t46*t47*t9
     +        -3.024470069787737d+1*gammaaa*t1*t20*t28*t30*t47*t9
     +        -3.240503646201146d+0*t1*t25*t4*t46*t9
     +        -2.520391724823114d+1*gammaaa*t1*t25*t30*t34*t9
     +        +2.520391724823114d+1*t1*t15*t34*t9
     +        +1.5122350348938685d+1*t1*t20*t25*t28*t9
     +        +3.240503646201146d+0*gammaaa*t1*t15*t4*(t87+t86+t85+t84
     +        +t83)+2.100172776208585d+1*gammaaa*t10*t12*t47*t48*t7
     +        -1.0500863881042925d+1*gammaaa*t10*t12*t25*t46*t7
     +        +4.900403144486698d+1*gammaaa*t10*t25*t26*t28*t7
     +        +8.167338574144497d+1*gammaaa*t10*t15*t7/rhoa*
     +        *6.666666666666667d+0-6.481007292402292d+0*gammaaa*t1
     +        *t24*t25*t4*t55-6.481007292402292d+0*gammaaa*t1*t25*t28
     +        *t4*t54-1.5122350348938685d+1*gammaaa*t1*t15*t20*t54
     +        +1.2962014584804585d+1*gammaaa*t1*t24*t28*t30*t4*t47
     +        -3.240503646201146d+0*gammaaa*t1*t25*t30*t4*t43
     +        +3.240503646201146d+0*t1*t15*t4*t43-6.481007292402292d+0
     +        *t1*t24*t25*t28*t4+1.5122350348938685d+1*gammaaa*t1*t20
     +        *t24*t25*t30-1.5122350348938685d+1*t1*t15*t20*t24)
     +        +3.1090690869654886d-2*t31*t51*t82+t50*(
     +        -3.1090690869654886d-2*t29*t56-1.5545345434827443d-2*t31
     +        *t49))*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = (1.5545345434827443d-2*t19*(
     +        -1.9443021877206879d+1*gammaaa*t1*t28*t4*t58*t81*t9
     +        +1.3611235077912603d+2*gammaaa*t13*t28*t47*t59*t60*t9
     +        +1.5879774257564705d+2*gammaaa*t13*t25*t59*t9/rhoa*
     +        *8.0d+0-1.5122350348938685d+1*gammaaa*t1*t20*t47*t58*t9
     +        +1.2962014584804585d+1*gammaaa*t1*t30*t4*t47*t55*t9
     +        -6.481007292402292d+0*t1*t25*t4*t55*t9
     +        -3.240503646201146d+0*gammaaa*t1*t25*t4*(
     +        -4.200345552417169d+1*t10*t12*t22*t27*t5
     +        -9.800806288973395d+1*t10*t13*t26)*t9
     +        +1.2962014584804585d+1*t1*t28*t30*t4*t47*t9
     +        +1.5122350348938685d+1*t1*t20*t25*t30*t9
     +        -2.100172776208585d+1*gammaaa*t10*t12*t25*t55*t7
     +        +4.20034555241717d+1*gammaaa*t10*t12*t28*t30*t47*t7
     +        +4.900403144486698d+1*gammaaa*t10*t25*t26*t30*t7
     +        -2.100172776208585d+1*t10*t12*t25*t28*t7
     +        -4.900403144486698d+1*t10*t15*t26*t7
     +        -6.805617538956302d+1*gammaaa*t13*t24*t25*t59*t60
     +        +6.481007292402292d+0*gammaaa*t1*t24*t4*t47*t58
     +        -6.481007292402292d+0*gammaaa*t1*t25*t30*t4*t54
     +        +6.481007292402292d+0*t1*t15*t4*t54-6.481007292402292d+0
     +        *t1*t24*t25*t30*t4)+3.1090690869654886d-2*t29*t57*t82
     +        +t50*(-1.5545345434827443d-2*t29*t61
     +        -3.1090690869654886d-2*t31*t56))*wght+Cmat3(iq,
     +        D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 1.5545345434827443d-2*t19*(
     +        -1.9443021877206879d+1*gammaaa*t1*t30**3*t4*t81*t9
     +        +4.0833705233737816d+2*gammaaa*t13*t30*t47*t59*t60*t9
     +        -2.0416852616868908d+2*t13*t25*t59*t60*t9
     +        +1.9443021877206879d+1*t1*t4*t47*t58*t9
     +        +6.300518328625755d+1*gammaaa*t10*t12*t47*t58*t7
     +        -6.300518328625755d+1*t10*t12*t25*t30*t7
     +        -6.616088534891511d+2*gammaaa*t1**4*t25*t27/rhoa*
     +        *9.333333333333333d+0)*wght+3.1090690869654886d-2*t31**3
     +        *t82*wght-4.6636036304482326d-2*t31*t50*t61*wght+Cmat3
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
            t3 = 1/rhob
            t4 = 1/rhob**2.3333333333333334d+0
            t5 = exp(-6.432793688582967d+1*s1f*t3)
            t6 = t5-1.0d+0
            t7 = 1/t6
            t8 = 3.240503646201146d+0*gammabb*t1*t4*t7
            t9 = t8+1.0d+0
            t10 = t1**2
            t11 = gammabb**2
            t12 = 1/rhob**4.666666666666667d+0
            t13 = 1/t6**2
            t14 = t8+1.0500863881042923d+1*t10*t11*t12*t13+1.0d+0
            t15 = 1/t14
            t16 = 3.240503646201146d+0*gammabb*t1*t15*t4*t9+1.0d+0
            t17 = 1/rhob**2
            t18 = s1a(D1_RB)
            t19 = 1/t16
            t20 = 1/rhob**3.3333333333333337d+0
            t21 = -7.561175174469341d+0*gammabb*t1*t20*t7
            t22 = 6.432793688582967d+1*t17*s1f-6.432793688582967d+1
     +        *t18*t3
            t23 = -3.240503646201146d+0*gammabb*t1*t13*t22*t4*t5
            t24 = t23+t21
            t25 = 1/t14**2
            t26 = 1/rhob**5.666666666666667d+0
            t27 = 1/t6**3
            t28 = -2.1001727762085845d+1*t10*t11*t12*t22*t27*t5
     +        -4.9004031444866974d+1*t10*t11*t13*t26+t23+t21
            t29 = -3.240503646201146d+0*gammabb*t1*t25*t28*t4*t9
     +        -7.561175174469341d+0*gammabb*t1*t15*t20*t9
     +        +3.240503646201146d+0*gammabb*t1*t15*t24*t4
            t30 = 3.240503646201146d+0*t1*t4*t7+2.1001727762085845d+1
     +        *gammabb*t10*t12*t13
            t31 = -3.240503646201146d+0*gammabb*t1*t25*t30*t4*t9
     +        +3.240503646201146d+0*t1*t15*t4*t9+1.0500863881042925d+1
     +        *gammabb*t10*t12*t15*t7
            t32 = 1/rhob**3
            t33 = s1a2(D2_RB_RB)
            t34 = 1/rhob**4.333333333333334d+0
            t35 = t22**2
            t36 = exp(-1.2865587377165935d+2*s1f*t3)
            t37 = 6.481007292402292d+0*gammabb*t1*t27*t35*t36*t4
            t38 = 2.520391724823114d+1*gammabb*t1*t34*t7
            t39 = 1.5122350348938685d+1*gammabb*t1*t13*t20*t22*t5
            t40 = -3.240503646201146d+0*gammabb*t1*t13*t35*t4*t5
            t41 = -6.432793688582967d+1*t3*t33-1.2865587377165935d+2
     +        *s1f*t32+1.2865587377165935d+2*t17*t18
            t42 = -3.240503646201146d+0*gammabb*t1*t13*t4*t41*t5
            t43 = t42+t40+t39+t38+t37
            t44 = 1/t6**4
            t45 = 1/rhob**6.666666666666667d+0
            t46 = -2.1001727762085845d+1*t10*t11*t12*t27*t41*t5
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t35*t5
     +        +1.960161257794679d+2*t10*t11*t22*t26*t27*t5
     +        +2.7768951152091287d+2*t10*t11*t13*t45
     +        +6.300518328625754d+1*t10*t11*t12*t35*t36*t44+t42+t40
     +        +t39+t38+t37
            t47 = 1/t14**3
            t48 = t28**2
            t49 = 6.481007292402292d+0*gammabb*t1*t4*t47*t48*t9
     +        -3.240503646201146d+0*gammabb*t1*t25*t4*t46*t9
     +        +2.520391724823114d+1*gammabb*t1*t15*t34*t9
     +        +1.5122350348938685d+1*gammabb*t1*t20*t25*t28*t9
     +        +3.240503646201146d+0*gammabb*t1*t15*t4*t43
     +        -6.481007292402292d+0*gammabb*t1*t24*t25*t28*t4
     +        -1.5122350348938685d+1*gammabb*t1*t15*t20*t24
            t50 = 1/t16**2
            t51 = t29**2
            t52 = -7.561175174469341d+0*t1*t20*t7
            t53 = -3.240503646201146d+0*t1*t13*t22*t4*t5
            t54 = t53+t52
            t55 = t53+t52-4.200345552417169d+1*gammabb*t10*t12*t22*t27
     +        *t5-9.800806288973395d+1*gammabb*t10*t13*t26
            t56 = -3.240503646201146d+0*gammabb*t1*t25*t4*t55*t9
     +        +6.481007292402292d+0*gammabb*t1*t28*t30*t4*t47*t9
     +        -3.240503646201146d+0*t1*t25*t28*t4*t9
     +        +7.561175174469341d+0*gammabb*t1*t20*t25*t30*t9
     +        -7.561175174469341d+0*t1*t15*t20*t9
     +        -1.0500863881042925d+1*gammabb*t10*t12*t25*t28*t7
     +        -2.450201572243349d+1*gammabb*t10*t15*t26*t7
     +        +3.240503646201146d+0*gammabb*t1*t15*t4*t54
     +        -3.240503646201146d+0*gammabb*t1*t24*t25*t30*t4
     +        +3.240503646201146d+0*t1*t15*t24*t4
            t57 = t31**2
            t58 = t30**2
            t59 = t1**3
            t60 = 1/rhob**7.0d+0
            t61 = -6.805617538956302d+1*gammabb*t13*t25*t59*t60*t9
     +        +6.481007292402292d+0*gammabb*t1*t4*t47*t58*t9
     +        -6.481007292402292d+0*t1*t25*t30*t4*t9
     +        -2.100172776208585d+1*gammabb*t10*t12*t25*t30*t7
     +        +2.100172776208585d+1*t10*t12*t15*t7
            t62 = 1/rhob**4
            t63 = s1a3(D3_RB_RB_RB)
            t64 = 1/rhob**5.333333333333334d+0
            t65 = t22**3
            t66 = exp(-1.92983810657489d+2*s1f*t3)
            t67 = -1.9443021877206879d+1*gammabb*t1*t4*t44*t65*t66
            t68 = -4.5367051046816054d+1*gammabb*t1*t20*t27*t35*t36
            t69 = 1.2865587377165935d+2*t17*s1f-1.2865587377165935d+2
     +        *t18*t3
            t70 = 6.481007292402292d+0*gammabb*t1*t27*t35*t36*t4*t69
            t71 = 6.481007292402292d+0*gammabb*t1*t27*t36*t4*t65
            t72 = 1.9443021877206879d+1*gammabb*t1*t22*t27*t36*t4*t41
            t73 = -1.0921697474233495d+2*gammabb*t1*t64*t7
            t74 = -7.561175174469342d+1*gammabb*t1*t13*t22*t34*t5
            t75 = 2.2683525523408027d+1*gammabb*t1*t13*t20*t35*t5
            t76 = 2.2683525523408027d+1*gammabb*t1*t13*t20*t41*t5
            t77 = -3.240503646201146d+0*gammabb*t1*t13*t4*t5*t65
            t78 = -9.721510938603439d+0*gammabb*t1*t13*t22*t4*t41*t5
            t79 = -6.432793688582967d+1*t3*t63+3.85967621314978d+2*s1f
     +        *t62+1.92983810657489d+2*t17*t33-3.85967621314978d+2*t18
     +        *t32
            t80 = -3.240503646201146d+0*gammabb*t1*t13*t4*t5*t79
            t81 = 1/t14**4
            t82 = 1/t16**3
            t83 = 6.481007292402292d+0*t1*t27*t35*t36*t4
            t84 = 2.520391724823114d+1*t1*t34*t7
            t85 = 1.5122350348938685d+1*t1*t13*t20*t22*t5
            t86 = -3.240503646201146d+0*t1*t13*t35*t4*t5
            t87 = -3.240503646201146d+0*t1*t13*t4*t41*t5
            fnc(iq) = (1.0d+0*s1f*t3+1.5545345434827443d-2*log(t16))
     +        *wght+fnc(iq)
            Amat(iq,D1_RB) = (1.0d+0*t18*t3+1.5545345434827443d-2*t19
     +        *t29-1.0d+0*t17*s1f)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.5545345434827443d-2*t19*t31*wght+Cmat
     +        (iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = (-1.5545345434827443d-2*t50*t51
     +        +1.5545345434827443d-2*t19*t49+1.0d+0*t3*t33+2.0d+0*s1f
     +        *t32-2.0d+0*t17*t18)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = (1.5545345434827443d-2*t19*t56
     +        -1.5545345434827443d-2*t29*t31*t50)*wght+Cmat2(iq,
     +        D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5545345434827443d-2*t19*t61*wght
     +        -1.5545345434827443d-2*t50*t57*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = (1.5545345434827443d-2*t19*(
     +        -1.9443021877206879d+1*gammabb*t1*t28**3*t4*t81*t9
     +        -3.240503646201146d+0*gammabb*t1*t25*t4*(t80
     +        -2.1001727762085845d+1*t10*t11*t12*t27*t5*t79+t78+t77
     +        +t76+t75+t74+t73+t72+t71+t70+6.300518328625754d+1*t10
     +        *t11*t12*t35*t36*t44*t69+t68+t67-2.5202073314503015d+2
     +        *t10*t11*t12*t65*t66/t6**5-2.1001727762085845d+1*t10*t11
     +        *t12*t27*t5*t65+6.300518328625754d+1*t10*t11*t12*t36*t44
     +        *t65-1.666137069125477d+3*t10*t11*t22*t27*t45*t5
     +        +2.940241886692019d+2*t10*t11*t26*t27*t41*t5
     +        -6.300518328625754d+1*t10*t11*t12*t22*t27*t41*t5
     +        +2.940241886692019d+2*t10*t11*t26*t27*t35*t5
     +        +1.890155498587726d+2*t10*t11*t12*t22*t36*t41*t44
     +        -8.820725660076056d+2*t10*t11*t26*t35*t36*t44
     +        -1.851263410139419d+3*t10*t11*t13/rhob*
     +        *7.666666666666667d+0)*t9-1.0921697474233495d+2*gammabb
     +        *t1*t15*t64*t9-4.5367051046816054d+1*gammabb*t1*t20*t47
     +        *t48*t9+1.9443021877206879d+1*gammabb*t1*t28*t4*t46*t47
     +        *t9+2.2683525523408027d+1*gammabb*t1*t20*t25*t46*t9
     +        -7.561175174469342d+1*gammabb*t1*t25*t28*t34*t9
     +        +3.240503646201146d+0*gammabb*t1*t15*t4*(t80+t78+t77+t76
     +        +t75+t74+t73+t72+t71+t70+t68+t67)+1.9443021877206879d+1
     +        *gammabb*t1*t24*t4*t47*t48-9.721510938603439d+0*gammabb
     +        *t1*t24*t25*t4*t46-9.721510938603439d+0*gammabb*t1*t25
     +        *t28*t4*t43-2.2683525523408027d+1*gammabb*t1*t15*t20*t43
     +        +7.561175174469342d+1*gammabb*t1*t15*t24*t34
     +        +4.5367051046816054d+1*gammabb*t1*t20*t24*t25*t28)
     +        +3.1090690869654886d-2*t29**3*t82+1.0d+0*t3*t63-6.0d+0
     +        *s1f*t62-4.6636036304482326d-2*t29*t49*t50-3.0d+0*t17
     +        *t33+6.0d+0*t18*t32)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = (1.5545345434827443d-2*t19*(
     +        -3.240503646201146d+0*gammabb*t1*t25*t4*(t87+t86+t85+t84
     +        +t83-4.200345552417169d+1*gammabb*t10*t12*t27*t41*t5
     +        -4.200345552417169d+1*gammabb*t10*t12*t27*t35*t5
     +        +3.920322515589358d+2*gammabb*t10*t22*t26*t27*t5
     +        +5.553790230418257d+2*gammabb*t10*t13*t45
     +        +1.2601036657251508d+2*gammabb*t10*t12*t35*t36*t44)*t9
     +        -1.9443021877206879d+1*gammabb*t1*t30*t4*t48*t81*t9
     +        +1.2962014584804585d+1*gammabb*t1*t28*t4*t47*t55*t9
     +        +1.5122350348938685d+1*gammabb*t1*t20*t25*t55*t9
     +        +6.481007292402292d+0*t1*t4*t47*t48*t9
     +        +6.481007292402292d+0*gammabb*t1*t30*t4*t46*t47*t9
     +        -3.024470069787737d+1*gammabb*t1*t20*t28*t30*t47*t9
     +        -3.240503646201146d+0*t1*t25*t4*t46*t9
     +        -2.520391724823114d+1*gammabb*t1*t25*t30*t34*t9
     +        +2.520391724823114d+1*t1*t15*t34*t9
     +        +1.5122350348938685d+1*t1*t20*t25*t28*t9
     +        +3.240503646201146d+0*gammabb*t1*t15*t4*(t87+t86+t85+t84
     +        +t83)+2.100172776208585d+1*gammabb*t10*t12*t47*t48*t7
     +        -1.0500863881042925d+1*gammabb*t10*t12*t25*t46*t7
     +        +4.900403144486698d+1*gammabb*t10*t25*t26*t28*t7
     +        +8.167338574144497d+1*gammabb*t10*t15*t7/rhob*
     +        *6.666666666666667d+0-6.481007292402292d+0*gammabb*t1
     +        *t24*t25*t4*t55-6.481007292402292d+0*gammabb*t1*t25*t28
     +        *t4*t54-1.5122350348938685d+1*gammabb*t1*t15*t20*t54
     +        +1.2962014584804585d+1*gammabb*t1*t24*t28*t30*t4*t47
     +        -3.240503646201146d+0*gammabb*t1*t25*t30*t4*t43
     +        +3.240503646201146d+0*t1*t15*t4*t43-6.481007292402292d+0
     +        *t1*t24*t25*t28*t4+1.5122350348938685d+1*gammabb*t1*t20
     +        *t24*t25*t30-1.5122350348938685d+1*t1*t15*t20*t24)
     +        +3.1090690869654886d-2*t31*t51*t82+t50*(
     +        -3.1090690869654886d-2*t29*t56-1.5545345434827443d-2*t31
     +        *t49))*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = (1.5545345434827443d-2*t19*(
     +        -1.9443021877206879d+1*gammabb*t1*t28*t4*t58*t81*t9
     +        +1.3611235077912603d+2*gammabb*t13*t28*t47*t59*t60*t9
     +        +1.5879774257564705d+2*gammabb*t13*t25*t59*t9/rhob*
     +        *8.0d+0-1.5122350348938685d+1*gammabb*t1*t20*t47*t58*t9
     +        +1.2962014584804585d+1*gammabb*t1*t30*t4*t47*t55*t9
     +        -6.481007292402292d+0*t1*t25*t4*t55*t9
     +        -3.240503646201146d+0*gammabb*t1*t25*t4*(
     +        -4.200345552417169d+1*t10*t12*t22*t27*t5
     +        -9.800806288973395d+1*t10*t13*t26)*t9
     +        +1.2962014584804585d+1*t1*t28*t30*t4*t47*t9
     +        +1.5122350348938685d+1*t1*t20*t25*t30*t9
     +        -2.100172776208585d+1*gammabb*t10*t12*t25*t55*t7
     +        +4.20034555241717d+1*gammabb*t10*t12*t28*t30*t47*t7
     +        +4.900403144486698d+1*gammabb*t10*t25*t26*t30*t7
     +        -2.100172776208585d+1*t10*t12*t25*t28*t7
     +        -4.900403144486698d+1*t10*t15*t26*t7
     +        -6.805617538956302d+1*gammabb*t13*t24*t25*t59*t60
     +        +6.481007292402292d+0*gammabb*t1*t24*t4*t47*t58
     +        -6.481007292402292d+0*gammabb*t1*t25*t30*t4*t54
     +        +6.481007292402292d+0*t1*t15*t4*t54-6.481007292402292d+0
     +        *t1*t24*t25*t30*t4)+3.1090690869654886d-2*t29*t57*t82
     +        +t50*(-1.5545345434827443d-2*t29*t61
     +        -3.1090690869654886d-2*t31*t56))*wght+Cmat3(iq,
     +        D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 1.5545345434827443d-2*t19*(
     +        -1.9443021877206879d+1*gammabb*t1*t30**3*t4*t81*t9
     +        +4.0833705233737816d+2*gammabb*t13*t30*t47*t59*t60*t9
     +        -2.0416852616868908d+2*t13*t25*t59*t60*t9
     +        +1.9443021877206879d+1*t1*t4*t47*t58*t9
     +        +6.300518328625755d+1*gammabb*t10*t12*t47*t58*t7
     +        -6.300518328625755d+1*t10*t12*t25*t30*t7
     +        -6.616088534891511d+2*gammabb*t1**4*t25*t27/rhob*
     +        *9.333333333333333d+0)*wght+3.1090690869654886d-2*t31**3
     +        *t82*wght-4.6636036304482326d-2*t31*t50*t61*wght+Cmat3
     +        (iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_c_mpbe_d3 
C> @}
