C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_sop.F 
C> The nwxcm_c_sop functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_sop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{1.8610514726982\,{\it t_3}
C>    +1.8610514726982\,{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \left({\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,
C>    \left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}+0.5764\right)\,\rho_\beta}\over{{\it t_5}\,
C>    \left(3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}+11.99591941974898\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111
C>    , 5656 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_sop(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
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
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = t1*t2
            t4 = t3**2.0d+0
            t5 = 1/t4
            t6 = 1.0500052408963245d+0*t1*t2
            t7 = 1.415701855281521d+0*t1*t2+5.764d-1
            t8 = rhoa**2
            t9 = t6+8.658781460080348d-1*t4+3.183d-1
            t10 = 1/t9
            t11 = rhoa**1.3333333333333333d+0
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t5*t7*t8*wght/(t6+8.
     1         658781460080349d-1*t4+3.183d-1)
            Amat(iq,D1_RA) = 1.385405033612855d+1*t10*t5*(-8.33616803355
     1         3906d-2*rhoa*t7*wght-1.966921425173462d-2*t1*t11*wght)+1.
     2         1548969154728159d+0*(2.8862604866934494d-1*t1**2/t2+1.750
     3         0087348272073d-1*t1/rhoa**6.666666666666666d-1)*t5*t7*t8*
     4         wght/t9**2+3.849656384909386d-1*t1*t10*t11*t7*wght/t3**3.
     5         0d+0+Amat(iq,D1_RA)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = rhob**3.333333333333333d-1
            t4 = 1.8610514726982d+0*t3+1.8610514726982d+0*t2
            t5 = 1/t4
            t6 = t1*t2*t3*t5
            t7 = t6**2.0d+0
            t8 = 1/t7
            t9 = 1.1995919419748984d+1*t7+3.908227599821866d+0*t1*t2*t3*
     1         t5+3.183d-1
            t10 = 1/t9
            t11 = 5.269388045346497d+0*t1*t2*t3*t5+5.764d-1
            t12 = 1/t6**3.0d+0
            t13 = 1/t2
            t14 = 1/t4**2
            t15 = 1/rhoa**6.666666666666666d-1
            t16 = 3.333333333333333d-1*t1*t15*t3*t5-6.203504908994d-1*t1
     1         *t13*t14*t3
            t17 = 1/t9**2
            t18 = 1/rhob**6.666666666666666d-1
            t19 = 1/t3
            t20 = 3.333333333333333d-1*t1*t18*t2*t5-6.203504908994d-1*t1
     1         *t14*t19*t2
            fnc(iq) = fnc(iq)-8.336168033553906d-2*rhoa*rhob*t10*t11*t8*
     1         wght
            Amat(iq,D1_RA) = t10*t8*(-8.336168033553906d-2*rhoa*rhob*(1.
     1         7564626817821655d+0*t1*t15*t3*t5-3.268867460670129d+0*t1*
     2         t13*t14*t3)*wght-8.336168033553906d-2*rhob*t11*wght)+8.33
     3         6168033553906d-2*rhoa*rhob*t11*t17*(2.3991838839497967d+1
     4         *t1*t16*t2*t3*t5+1.302742533273955d+0*t1*t15*t3*t5-2.4244
     5         70910096078d+0*t1*t13*t14*t3)*t8*wght+1.6672336067107812d
     6         -1*rhoa*rhob*t10*t11*t12*t16*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = t10*t8*(-8.336168033553906d-2*rhoa*rhob*(1.
     1         7564626817821655d+0*t1*t18*t2*t5-3.268867460670129d+0*t1*
     2         t14*t19*t2)*wght-8.336168033553906d-2*rhoa*t11*wght)+8.33
     3         6168033553906d-2*rhoa*rhob*t11*t17*(2.3991838839497967d+1
     4         *t1*t2*t20*t3*t5+1.302742533273955d+0*t1*t18*t2*t5-2.4244
     5         70910096078d+0*t1*t14*t19*t2)*t8*wght+1.6672336067107812d
     6         -1*rhoa*rhob*t10*t11*t12*t20*wght+Amat(iq,D1_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_sop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{1.8610514726982\,{\it t_3}
C>    +1.8610514726982\,{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \left({\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,
C>    \left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}+0.5764\right)\,\rho_\beta}\over{{\it t_5}\,
C>    \left(3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}+11.99591941974898\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111
C>    , 5656 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_sop_d2(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
c 
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
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = t1*t2
            t4 = t3**2.0d+0
            t5 = 1/t4
            t6 = 1.0500052408963245d+0*t1*t2
            t7 = 1.415701855281521d+0*t1*t2+5.764d-1
            t8 = rhoa**2
            t9 = 1/t3**3.0d+0
            t10 = t6+8.658781460080348d-1*t4+3.183d-1
            t11 = 1/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = rhoa**6.666666666666666d-1
            t14 = 1/t13
            t15 = 1.7500087348272073d-1*t1*t14
            t16 = t1**2
            t17 = 1/t2
            t18 = 2.8862604866934494d-1*t16*t17+t15
            t19 = 1/t10**2
            t20 = -8.336168033553906d-2*rhoa*t7*wght-1.966921425173462d-
     1         2*t1*t12*wght
            t21 = 1/t3**4.0d+0
            t22 = 1/t10**3
            t23 = -8.955510855647054d-2*t1*t14*t20
            t24 = 8.336168033553906d-2*rhoa*t18*t7*wght
            t25 = 1.966921425173462d-2*t1*t12*t18*wght
            t26 = 1/rhoa**1.6666666666666669d+0
            t27 = 1/t12
            t28 = 2.886260486693449d-1*t16*t17+t15
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t5*t7*t8*wght/(t6+8.
     1         658781460080349d-1*t4+3.183d-1)
            Amat(iq,D1_RA) = 3.849656384909386d-1*t1*t11*t12*t7*t9*wght+
     1         1.1548969154728159d+0*t18*t19*t5*t7*t8*wght+1.38540503361
     2         2855d+1*t11*t20*t5+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = 1.385405033612855d+1*t19*t5*(8.33616803
     1         3553906d-2*(-2.4052170722445415d-1*t16*t27-1.750008734827
     2         2073d-1*t1*t26)*t7*t8*wght+t25+t24-t18*t20)+5.15662015617
     3         7406d+1*t11*t9*(-8.673617379884035d-19*t1*t2*t7*wght+1.76
     4         14786175345715d-3*t13*t16*wght+t23)-7.699312769818772d-1*
     5         t1*t12*t18*t19*t7*t9*wght-2.3097938309456317d+0*t18**2*t2
     6         2*t5*t7*t8*wght-1.9248281924546928d-1*t11*t13*t16*t21*t7*
     7         wght-2.724982843156285d-1*t1*t11*t2*t5*wght+Amat2(iq,D2_R
     8         A_RA)
            Amat2(iq,D2_RA_RB) = 1.385405033612855d+1*t19*t5*(8.33616803
     1         3553906d-2*(1.4431302433467244d-1*t16*t27+5.8333624494240
     2         24d-2*t1*t26)*t7*t8*wght-t20*t28+t25+t24)+5.1566201561774
     3         06d+1*t19*t9*(-7.465464331898997d-3*t1*t12*t28*t7*wght-7.
     4         465464331898997d-3*t1*t12*t18*t7*wght)+5.156620156177406d
     5         +1*t11*t9*(9.953952442531996d-3*t1*t2*t7*wght+1.761478617
     6         5345715d-3*t13*t16*wght+t23)+1.385405033612855d+1*t11*t5*
     7         (-8.336168033553906d-2*t7*wght-4.5894833254047446d-2*t1*t
     8         2*wght)-2.3097938309456317d+0*t18*t22*t28*t5*t7*t8*wght-1
     9         .924828192454693d-1*t11*t13*t16*t21*t7*wght+Amat2(iq,D2_R
     :         A_RB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = rhob**3.333333333333333d-1
            t4 = 1.8610514726982d+0*t3+1.8610514726982d+0*t2
            t5 = 1/t4
            t6 = t1*t2*t3*t5
            t7 = t6**2.0d+0
            t8 = 1/t7
            t9 = 1.1995919419748984d+1*t7+3.908227599821866d+0*t1*t2*t3*
     1         t5+3.183d-1
            t10 = 1/t9
            t11 = 5.269388045346497d+0*t1*t2*t3*t5+5.764d-1
            t12 = 1/t6**3.0d+0
            t13 = 1/t2
            t14 = 1/t4**2
            t15 = 1/rhoa**6.666666666666666d-1
            t16 = 3.333333333333333d-1*t1*t15*t3*t5-6.203504908994d-1*t1
     1         *t13*t14*t3
            t17 = 1/t9**2
            t18 = 2.3991838839497967d+1*t1*t16*t2*t3*t5+1.30274253327395
     1         5d+0*t1*t15*t3*t5-2.424470910096078d+0*t1*t13*t14*t3
            t19 = 1.7564626817821655d+0*t1*t15*t3*t5-3.268867460670129d+
     1         0*t1*t13*t14*t3
            t20 = -8.336168033553906d-2*rhoa*rhob*t19*wght-8.33616803355
     1         3906d-2*rhob*t11*wght
            t21 = 1/rhob**6.666666666666666d-1
            t22 = 1/t3
            t23 = 3.333333333333333d-1*t1*t2*t21*t5-6.203504908994d-1*t1
     1         *t14*t2*t22
            t24 = 2.3991838839497967d+1*t1*t2*t23*t3*t5+1.30274253327395
     1         5d+0*t1*t2*t21*t5-2.424470910096078d+0*t1*t14*t2*t22
            t25 = 1.7564626817821655d+0*t1*t2*t21*t5-3.268867460670129d+
     1         0*t1*t14*t2*t22
            t26 = -8.336168033553906d-2*rhoa*rhob*t25*wght-8.33616803355
     1         3906d-2*rhoa*t11*wght
            t27 = 1/t9**3
            t28 = 1/t6**4.0d+0
            t29 = 1/rhoa
            t30 = 1/t4**3
            t31 = 1/rhoa**1.6666666666666669d+0
            t32 = 7.69669463118253d-1*t1*t29*t3*t30-2.222222222222222d-1
     1         *t1*t3*t31*t5
            t33 = 1.111111111111111d-1*t1*t15*t21*t5+7.69669463118253d-1
     1         *t1*t13*t22*t30-2.0678349696646664d-1*t1*t14*t15*t22-2.06
     2         78349696646664d-1*t1*t13*t14*t21
            t34 = 1/rhob**1.6666666666666669d+0
            t35 = 1/rhob
            t36 = 7.69669463118253d-1*t1*t2*t30*t35-2.222222222222222d-1
     1         *t1*t2*t34*t5
            fnc(iq) = fnc(iq)-8.336168033553906d-2*rhoa*rhob*t10*t11*t8*
     1         wght
            Amat(iq,D1_RA) = 8.336168033553906d-2*rhoa*rhob*t11*t17*t18*
     1         t8*wght+1.6672336067107812d-1*rhoa*rhob*t10*t11*t12*t16*w
     2         ght+t10*t20*t8+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.336168033553906d-2*rhoa*rhob*t11*t17*t24*
     1         t8*wght+1.6672336067107812d-1*rhoa*rhob*t10*t11*t12*t23*w
     2         ght+t10*t26*t8+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = t17*t8*(8.336168033553906d-2*rhoa*rhob*
     1         t11*(2.3991838839497967d+1*t1*t2*t3*t32*t5-8.684950221826
     2         366d-1*t1*t3*t31*t5+7.997279613165988d+0*t1*t15*t16*t3*t5
     3         +3.008043438498834d+0*t1*t29*t3*t30-1.4883349001661852d+1
     4         *t1*t13*t14*t16*t3)*wght+8.336168033553906d-2*rhoa*rhob*t
     5         18*t19*wght+8.336168033553906d-2*rhob*t11*t18*wght-t18*t2
     6         0)+t10*t8*(-8.336168033553906d-2*rhoa*rhob*(4.05568706782
     7         35785d+0*t1*t29*t3*t30-1.17097512118811d+0*t1*t3*t31*t5)*
     8         wght-1.6672336067107812d-1*rhob*t19*wght)+t10*t12*(1.6672
     9         336067107812d-1*rhoa*rhob*t11*t32*wght+1.6672336067107812
     :         d-1*rhoa*rhob*t16*t19*wght+1.6672336067107812d-1*rhob*t11
     ;         *t16*wght-2.0d+0*t16*t20)-1.6672336067107812d-1*rhoa*rhob
     <         *t11*t18**2*t27*t8*wght-5.001700820132344d-1*rhoa*rhob*t1
     =         0*t11*t16**2*t28*wght-3.3344672134215625d-1*rhoa*rhob*t11
     >         *t12*t16*t17*t18*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t17*t8*(8.336168033553906d-2*rhoa*rhob*
     1         t11*(2.3991838839497967d+1*t1*t2*t3*t33*t5+7.997279613165
     2         988d+0*t1*t16*t2*t21*t5+4.342475110913183d-1*t1*t15*t21*t
     3         5+3.008043438498834d+0*t1*t13*t22*t30-1.4883349001661852d
     4         +1*t1*t14*t16*t2*t22-8.081569700320259d-1*t1*t14*t15*t22-
     5         8.081569700320259d-1*t1*t13*t14*t21)*wght+8.3361680335539
     6         06d-2*rhoa*rhob*t18*t25*wght+8.336168033553906d-2*rhoa*t1
     7         1*t18*wght-t20*t24)+t10*t8*(-8.336168033553906d-2*rhoa*rh
     8         ob*(5.85487560594055d-1*t1*t15*t21*t5+4.0556870678235785d
     9         +0*t1*t13*t22*t30-1.0896224868900428d+0*t1*t14*t15*t22-1.
     :         0896224868900428d+0*t1*t13*t14*t21)*wght-8.33616803355390
     ;         6d-2*rhob*t25*wght-8.336168033553906d-2*rhoa*t19*wght-8.3
     <         36168033553906d-2*t11*wght)+t10*t12*(1.6672336067107812d-
     =         1*rhoa*rhob*t11*t33*wght+1.6672336067107812d-1*rhoa*rhob*
     >         t16*t25*wght+1.6672336067107812d-1*rhoa*t11*t16*wght-2.0d
     ?         +0*t20*t23)+t12*t17*(-1.6672336067107812d-1*rhoa*rhob*t11
     @         *t16*t24*wght-1.6672336067107812d-1*rhoa*rhob*t11*t18*t23
     1         *wght)-1.6672336067107812d-1*rhoa*rhob*t11*t18*t24*t27*t8
     2         *wght-5.001700820132344d-1*rhoa*rhob*t10*t11*t16*t23*t28*
     3         wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = t17*t8*(8.336168033553906d-2*rhoa*rhob*
     1         t11*(2.3991838839497967d+1*t1*t2*t3*t36*t5-8.684950221826
     2         366d-1*t1*t2*t34*t5+7.997279613165988d+0*t1*t2*t21*t23*t5
     3         +3.008043438498834d+0*t1*t2*t30*t35-1.4883349001661852d+1
     4         *t1*t14*t2*t22*t23)*wght+8.336168033553906d-2*rhoa*rhob*t
     5         24*t25*wght+8.336168033553906d-2*rhoa*t11*t24*wght-t24*t2
     6         6)+t10*t8*(-8.336168033553906d-2*rhoa*rhob*(4.05568706782
     7         35785d+0*t1*t2*t30*t35-1.17097512118811d+0*t1*t2*t34*t5)*
     8         wght-1.6672336067107812d-1*rhoa*t25*wght)+t10*t12*(1.6672
     9         336067107812d-1*rhoa*rhob*t11*t36*wght+1.6672336067107812
     :         d-1*rhoa*rhob*t23*t25*wght+1.6672336067107812d-1*rhoa*t11
     ;         *t23*wght-2.0d+0*t23*t26)-1.6672336067107812d-1*rhoa*rhob
     <         *t11*t24**2*t27*t8*wght-5.001700820132344d-1*rhoa*rhob*t1
     =         0*t11*t23**2*t28*wght-3.3344672134215625d-1*rhoa*rhob*t11
     >         *t12*t17*t23*t24*wght+Amat2(iq,D2_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_sop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{1.8610514726982\,{\it t_3}
C>    +1.8610514726982\,{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \left({\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,
C>    \left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}+0.5764\right)\,\rho_\beta}\over{{\it t_5}\,
C>    \left(3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{
C>    \it t_3}+11.99591941974898\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111
C>    , 5656 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_sop_d3(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2,Amat3) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho 
c 
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
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = t1*t2
            t4 = t3**2.0d+0
            t5 = 1/t4
            t6 = 1.0500052408963245d+0*t1*t2
            t7 = 1.415701855281521d+0*t1*t2+5.764d-1
            t8 = rhoa**2
            t9 = 1/t3**3.0d+0
            t10 = t6+8.658781460080348d-1*t4+3.183d-1
            t11 = 1/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = rhoa**6.666666666666666d-1
            t14 = 1/t13
            t15 = 1.7500087348272073d-1*t1*t14
            t16 = t1**2
            t17 = 1/t2
            t18 = 2.8862604866934494d-1*t16*t17+t15
            t19 = 1/t10**2
            t20 = -8.336168033553906d-2*rhoa*t7*wght-1.966921425173462d-
     1         2*t1*t12*wght
            t21 = 1/t3**4.0d+0
            t22 = t18**2
            t23 = 1/t10**3
            t24 = -8.955510855647054d-2*t1*t14*t20
            t25 = -8.673617379884035d-19*t1*t2*t7*wght+1.761478617534571
     1         5d-3*t13*t16*wght+t24
            t26 = 8.336168033553906d-2*rhoa*t18*t7*wght
            t27 = 1.966921425173462d-2*t1*t12*t18*wght
            t28 = 1/rhoa**1.6666666666666669d+0
            t29 = 1/t12
            t30 = -2.4052170722445415d-1*t16*t29-1.7500087348272073d-1*t
     1         1*t28
            t31 = 8.336168033553906d-2*t30*t7*t8*wght+t27+t26-t18*t20
            t32 = 2.886260486693449d-1*t16*t17+t15
            t33 = -8.336168033553906d-2*t7*wght-4.5894833254047446d-2*t1
     1         *t2*wght
            t34 = 1.4431302433467244d-1*t16*t29+5.833362449424024d-2*t1*
     1         t28
            t35 = t1**3
            t36 = 1/t3**5.0d+0
            t37 = 1/t10**4
            t38 = 1/rhoa**2.6666666666666666d+0
            t39 = 1/rhoa**2.3333333333333334d+0
            t40 = -2.3662411321981533d-4*t35*wght
            t41 = -1.3433266283470582d-1*t1*t14*t25
            t42 = -8.955510855647054d-2*t1*t14*t31
            t43 = -1.6672336067107812d-1*rhoa*t22*t7*wght
            t44 = -3.933842850346924d-2*t1*t12*t22*wght
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t5*t7*t8*wght/(t6+8.
     1         658781460080349d-1*t4+3.183d-1)
            Amat(iq,D1_RA) = 3.849656384909386d-1*t1*t11*t12*t7*t9*wght+
     1         1.1548969154728159d+0*t18*t19*t5*t7*t8*wght+1.38540503361
     2         2855d+1*t11*t20*t5+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = -7.699312769818772d-1*t1*t12*t18*t19*t7
     1         *t9*wght-2.3097938309456317d+0*t22*t23*t5*t7*t8*wght-1.92
     2         48281924546928d-1*t11*t13*t16*t21*t7*wght-2.7249828431562
     3         85d-1*t1*t11*t2*t5*wght+5.156620156177406d+1*t11*t25*t9+1
     4         .385405033612855d+1*t19*t31*t5+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.385405033612855d+1*t19*t5*(8.33616803
     1         3553906d-2*t34*t7*t8*wght-t20*t32+t27+t26)+5.156620156177
     2         406d+1*t19*t9*(-7.465464331898997d-3*t1*t12*t32*t7*wght-7
     3         .465464331898997d-3*t1*t12*t18*t7*wght)+5.156620156177406
     4         d+1*t11*t9*(9.953952442531996d-3*t1*t2*t7*wght+1.76147861
     5         75345715d-3*t13*t16*wght+t24)-2.3097938309456317d+0*t18*t
     6         23*t32*t5*t7*t8*wght-1.924828192454693d-1*t11*t13*t16*t21
     7         *t7*wght+1.385405033612855d+1*t11*t33*t5+Amat2(iq,D2_RA_R
     8         B)
            Amat3(iq,D3_RA_RA_RA) = 1.385405033612855d+1*t19*t5*(8.33616
     1         8033553906d-2*(4.1690429252238714d-1*t16*t39+3.4027947621
     2         64014d-1*t1*t38)*t7*t8*wght+1.6672336067107812d-1*rhoa*t3
     3         0*t7*wght+3.933842850346924d-2*t1*t12*t30*wght+5.90076427
     4         5520386d-2*t1*t18*t2*wght-t20*t30)+1.385405033612855d+1*t
     5         23*t5*(-3.3344672134215625d-1*t18*t30*t7*t8*wght+t44+t43-
     6         2*t18*t31)+5.156620156177406d+1*t19*t9*(-1.49309286637979
     7         93d-2*t1*t12*t30*t7*wght+1.734723475976807d-18*t1*t18*t2*
     8         t7*wght-3.522957235069143d-3*t13*t16*t18*wght+t42-t18*t25
     9         )+1.9193471071598367d+2*t11*t21*(1.0028557030015108d-3*t1
     :         6*t17*t7*wght+t41+t40)+5.156620156177406d+1*t11*t9*(-4.14
     ;         748018438835d-4*t1*t14*t7*wght+1.7614786175345712d-3*t16*
     <         t17*wght+8.955510855647055d-2*t1*t20*t28)+2.3097938309456
     =         323d+0*t1*t12*t22*t23*t7*t9*wght+6.929381492836896d+0*t18
     >         **3*t37*t5*t7*t8*wght+1.2832187949697954d-1*t11*t35*t36*t
     ?         7*wght+5.77448457736408d-1*t13*t16*t18*t19*t21*t7*wght+2.
     @         876370778887191d-1*t1*t11*t14*t5*wght+Amat3(iq,D3_RA_RA_R
     1         A)
            Amat3(iq,D3_RA_RA_RB) = 1.385405033612855d+1*t19*t5*(8.33616
     1         8033553906d-2*(-9.620868288978159d-2*t16*t39-4.8611353745
     2         20021d-2*t1*t38)*t7*t8*wght+8.336168033553906d-2*rhoa*t34
     3         *t7*wght+8.336168033553906d-2*rhoa*t30*t7*wght+8.33616803
     4         3553906d-2*t18*t7*wght+1.966921425173462d-2*t1*t12*t34*wg
     5         ht+1.966921425173462d-2*t1*t2*t32*wght+1.966921425173462d
     6         -2*t1*t12*t30*wght+4.5894833254047446d-2*t1*t18*t2*wght-t
     7         20*t34-t18*t33)+1.385405033612855d+1*t23*t5*(-3.334467213
     8         4215625d-1*t18*t34*t7*t8*wght+t44+t43-2*t31*t32)+5.156620
     9         156177406d+1*t19*t9*(-1.4930928663797993d-2*t1*t12*t34*t7
     :         *wght-1.9907904885063993d-2*t1*t18*t2*t7*wght-3.522957235
     ;         069143d-3*t13*t16*t18*wght+t42-t25*t32)+5.156620156177406
     <         d+1*t23*t9*(2.9861857327595986d-2*t1*t12*t18*t32*t7*wght+
     =         1.4930928663797993d-2*t1*t12*t22*t7*wght)+1.9193471071598
     >         367d+2*t19*t21*(1.0028557030015103d-3*t13*t16*t32*t7*wght
     ?         +2.0057114060030207d-3*t13*t16*t18*t7*wght)+1.91934710715
     @         98367d+2*t11*t21*(-1.671426171669184d-3*t16*t17*t7*wght+t
     1         41+t40)+5.156620156177406d+1*t11*t9*(4.147480184388324d-4
     2         *t1*t14*t7*wght+4.6972763134255235d-3*t16*t17*wght-8.9555
     3         10855647054d-2*t1*t14*t33-2.985170285215685d-2*t1*t20*t28
     4         )+6.929381492836896d+0*t22*t32*t37*t5*t7*t8*wght+1.283218
     5         794969795d-1*t11*t35*t36*t7*wght-3.7846983932726175d-1*t1
     6         *t11*t14*t5*wght+Amat3(iq,D3_RA_RA_RB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = rhob**3.333333333333333d-1
            t4 = 1.8610514726982d+0*t3+1.8610514726982d+0*t2
            t5 = 1/t4
            t6 = t1*t2*t3*t5
            t7 = t6**2.0d+0
            t8 = 1/t7
            t9 = 1.1995919419748984d+1*t7+3.908227599821866d+0*t1*t2*t3*
     1         t5+3.183d-1
            t10 = 1/t9
            t11 = 5.269388045346497d+0*t1*t2*t3*t5+5.764d-1
            t12 = 1/t6**3.0d+0
            t13 = 1/t2
            t14 = 1/t4**2
            t15 = 1/rhoa**6.666666666666666d-1
            t16 = 3.333333333333333d-1*t1*t15*t3*t5-6.203504908994d-1*t1
     1         *t13*t14*t3
            t17 = 1/t9**2
            t18 = 2.3991838839497967d+1*t1*t16*t2*t3*t5+1.30274253327395
     1         5d+0*t1*t15*t3*t5-2.424470910096078d+0*t1*t13*t14*t3
            t19 = 1.7564626817821655d+0*t1*t15*t3*t5-3.268867460670129d+
     1         0*t1*t13*t14*t3
            t20 = -8.336168033553906d-2*rhoa*rhob*t19*wght-8.33616803355
     1         3906d-2*rhob*t11*wght
            t21 = 1/rhob**6.666666666666666d-1
            t22 = 1/t3
            t23 = 3.333333333333333d-1*t1*t2*t21*t5-6.203504908994d-1*t1
     1         *t14*t2*t22
            t24 = 2.3991838839497967d+1*t1*t2*t23*t3*t5+1.30274253327395
     1         5d+0*t1*t2*t21*t5-2.424470910096078d+0*t1*t14*t2*t22
            t25 = 1.7564626817821655d+0*t1*t2*t21*t5-3.268867460670129d+
     1         0*t1*t14*t2*t22
            t26 = -8.336168033553906d-2*rhoa*rhob*t25*wght-8.33616803355
     1         3906d-2*rhoa*t11*wght
            t27 = 1/t9**3
            t28 = t18**2
            t29 = 1/t6**4.0d+0
            t30 = t16**2
            t31 = 1/rhoa
            t32 = 1/t4**3
            t33 = 1/rhoa**1.6666666666666669d+0
            t34 = 4.0556870678235785d+0*t1*t3*t31*t32-1.17097512118811d+
     1         0*t1*t3*t33*t5
            t35 = -8.336168033553906d-2*rhoa*rhob*t34*wght-1.66723360671
     1         07812d-1*rhob*t19*wght
            t36 = 7.69669463118253d-1*t1*t3*t31*t32-2.222222222222222d-1
     1         *t1*t3*t33*t5
            t37 = 1.6672336067107812d-1*rhoa*rhob*t11*t36*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t16*t19*wght+1.6672336067107812d-1
     2         *rhob*t11*t16*wght-2.0d+0*t16*t20
            t38 = 2.3991838839497967d+1*t1*t2*t3*t36*t5-8.68495022182636
     1         6d-1*t1*t3*t33*t5+7.997279613165988d+0*t1*t15*t16*t3*t5+3
     2         .008043438498834d+0*t1*t3*t31*t32-1.4883349001661852d+1*t
     3         1*t13*t14*t16*t3
            t39 = 8.336168033553906d-2*rhoa*rhob*t11*t38*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t18*t19*wght+8.336168033553906d-2*rh
     2         ob*t11*t18*wght-t18*t20
            t40 = 5.85487560594055d-1*t1*t15*t21*t5+4.0556870678235785d+
     1         0*t1*t13*t22*t32-1.0896224868900428d+0*t1*t14*t15*t22-1.0
     2         896224868900428d+0*t1*t13*t14*t21
            t41 = -8.336168033553906d-2*rhoa*rhob*t40*wght-8.33616803355
     1         3906d-2*rhob*t25*wght-8.336168033553906d-2*rhoa*t19*wght-
     2         8.336168033553906d-2*t11*wght
            t42 = -1.6672336067107812d-1*rhoa*rhob*t11*t16*t24*wght-1.66
     1         72336067107812d-1*rhoa*rhob*t11*t18*t23*wght
            t43 = 1.111111111111111d-1*t1*t15*t21*t5+7.69669463118253d-1
     1         *t1*t13*t22*t32-2.0678349696646664d-1*t1*t14*t15*t22-2.06
     2         78349696646664d-1*t1*t13*t14*t21
            t44 = 1.6672336067107812d-1*rhoa*rhob*t11*t43*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t16*t25*wght+1.6672336067107812d-1
     2         *rhoa*t11*t16*wght-2.0d+0*t20*t23
            t45 = 2.3991838839497967d+1*t1*t2*t3*t43*t5+7.99727961316598
     1         8d+0*t1*t16*t2*t21*t5+4.342475110913183d-1*t1*t15*t21*t5+
     2         3.008043438498834d+0*t1*t13*t22*t32-1.4883349001661852d+1
     3         *t1*t14*t16*t2*t22-8.081569700320259d-1*t1*t14*t15*t22-8.
     4         081569700320259d-1*t1*t13*t14*t21
            t46 = 8.336168033553906d-2*rhoa*rhob*t11*t45*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t18*t25*wght+8.336168033553906d-2*rh
     2         oa*t11*t18*wght-t20*t24
            t47 = t24**2
            t48 = t23**2
            t49 = 1/rhob**1.6666666666666669d+0
            t50 = 1/rhob
            t51 = 4.0556870678235785d+0*t1*t2*t32*t50-1.17097512118811d+
     1         0*t1*t2*t49*t5
            t52 = -8.336168033553906d-2*rhoa*rhob*t51*wght-1.66723360671
     1         07812d-1*rhoa*t25*wght
            t53 = 7.69669463118253d-1*t1*t2*t32*t50-2.222222222222222d-1
     1         *t1*t2*t49*t5
            t54 = 1.6672336067107812d-1*rhoa*rhob*t11*t53*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t23*t25*wght+1.6672336067107812d-1
     2         *rhoa*t11*t23*wght-2.0d+0*t23*t26
            t55 = 2.3991838839497967d+1*t1*t2*t3*t5*t53+3.00804343849883
     1         4d+0*t1*t2*t32*t50-8.684950221826366d-1*t1*t2*t49*t5+7.99
     2         7279613165988d+0*t1*t2*t21*t23*t5-1.4883349001661852d+1*t
     3         1*t14*t2*t22*t23
            t56 = 8.336168033553906d-2*rhoa*rhob*t11*t55*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t24*t25*wght+8.336168033553906d-2*rh
     2         oa*t11*t24*wght-t24*t26
            t57 = 1/t9**4
            t58 = 1/t6**5.0d+0
            t59 = 1/t4**4
            t60 = 1/rhoa**2
            t61 = 1/rhoa**2.3333333333333334d+0
            t62 = 1/rhoa**2.6666666666666666d+0
            t63 = 3.7037037037037035d-1*t1*t3*t5*t62+1.378556646443111d-
     1         1*t1*t14*t3*t61-7.69669463118253d-1*t1*t3*t32*t60-1.43239
     2         44878270575d+0*t1*t3*t33*t59
            t64 = -1.4323944878270575d+0*t1*t22*t31*t59-7.40740740740740
     1         7d-2*t1*t21*t33*t5+1.378556646443111d-1*t1*t14*t22*t33+2.
     2         565564877060843d-1*t1*t21*t31*t32
            t65 = -1.4323944878270575d+0*t1*t13*t50*t59+2.56556487706084
     1         36d-1*t1*t15*t32*t50-7.407407407407407d-2*t1*t15*t49*t5+1
     2         .378556646443111d-1*t1*t13*t14*t49+5.551115123125783d-17*
     3         t1*t13*t32/rhob**1.3333333333333333d+0
            t66 = 1/rhob**2.6666666666666666d+0
            t67 = 1/rhob**2.3333333333333334d+0
            t68 = 1/rhob**2
            t69 = -7.69669463118253d-1*t1*t2*t32*t68+1.378556646443111d-
     1         1*t1*t14*t2*t67+3.7037037037037035d-1*t1*t2*t5*t66-1.4323
     2         944878270575d+0*t1*t2*t49*t59
            fnc(iq) = fnc(iq)-8.336168033553906d-2*rhoa*rhob*t10*t11*t8*
     1         wght
            Amat(iq,D1_RA) = 8.336168033553906d-2*rhoa*rhob*t11*t17*t18*
     1         t8*wght+1.6672336067107812d-1*rhoa*rhob*t10*t11*t12*t16*w
     2         ght+t10*t20*t8+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.336168033553906d-2*rhoa*rhob*t11*t17*t24*
     1         t8*wght+1.6672336067107812d-1*rhoa*rhob*t10*t11*t12*t23*w
     2         ght+t10*t26*t8+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = -1.6672336067107812d-1*rhoa*rhob*t11*t2
     1         7*t28*t8*wght-5.001700820132344d-1*rhoa*rhob*t10*t11*t29*
     2         t30*wght-3.3344672134215625d-1*rhoa*rhob*t11*t12*t16*t17*
     3         t18*wght+t17*t39*t8+t10*t35*t8+t10*t12*t37+Amat2(iq,D2_RA
     4         _RA)
            Amat2(iq,D2_RA_RB) = -1.6672336067107812d-1*rhoa*rhob*t11*t1
     1         8*t24*t27*t8*wght-5.001700820132344d-1*rhoa*rhob*t10*t11*
     2         t16*t23*t29*wght+t17*t46*t8+t10*t41*t8+t10*t12*t44+t12*t1
     3         7*t42+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = -1.6672336067107812d-1*rhoa*rhob*t11*t2
     1         7*t47*t8*wght-5.001700820132344d-1*rhoa*rhob*t10*t11*t29*
     2         t48*wght-3.3344672134215625d-1*rhoa*rhob*t11*t12*t17*t23*
     3         t24*wght+t17*t56*t8+t10*t52*t8+t10*t12*t54+Amat2(iq,D2_RB
     4         _RB)
            Amat3(iq,D3_RA_RA_RA) = t17*t8*(8.336168033553906d-2*rhoa*rh
     1         ob*t11*(2.3991838839497967d+1*t1*t2*t3*t5*t63+1.447491703
     2         6377277d+0*t1*t3*t5*t62+5.387713133546839d-1*t1*t14*t3*t6
     3         1-3.008043438498834d+0*t1*t3*t32*t60-5.598123671158411d+0
     4         *t1*t3*t33*t59+1.5994559226331975d+1*t1*t15*t3*t36*t5-5.3
     5         31519742110659d+0*t1*t16*t3*t33*t5-2.9766698003323705d+1*
     6         t1*t13*t14*t3*t36+1.846578571881605d+1*t1*t16*t3*t31*t32)
     7         *wght+1.6672336067107812d-1*rhoa*rhob*t19*t38*wght+1.6672
     8         336067107812d-1*rhob*t11*t38*wght+8.336168033553906d-2*rh
     9         oa*rhob*t18*t34*wght+1.6672336067107812d-1*rhob*t18*t19*w
     :         ght-t20*t38-2*t18*t35)+t10*t12*(1.6672336067107812d-1*rho
     ;         a*rhob*t11*t63*wght+3.3344672134215625d-1*rhoa*rhob*t19*t
     <         36*wght+3.3344672134215625d-1*rhob*t11*t36*wght+1.6672336
     =         067107812d-1*rhoa*rhob*t16*t34*wght+3.3344672134215625d-1
     >         *rhob*t16*t19*wght-2.0d+0*t20*t36-4.0d+0*t16*t35)+t10*t8*
     ?         (-8.336168033553906d-2*rhoa*rhob*(1.9516252019801836d+0*t
     @         1*t3*t5*t62+7.264149912600285d-1*t1*t14*t3*t61-4.05568706
     1         78235785d+0*t1*t3*t32*t60-7.547842390376116d+0*t1*t3*t33*
     2         t59)*wght-2.500850410066172d-1*rhob*t34*wght)+t27*t8*(-3.
     3         3344672134215625d-1*rhoa*rhob*t11*t18*t38*wght-1.66723360
     4         67107812d-1*rhoa*rhob*t19*t28*wght-1.6672336067107812d-1*
     5         rhob*t11*t28*wght-2*t18*t39)+t12*t17*(-3.3344672134215625
     6         d-1*rhoa*rhob*t11*t16*t38*wght-3.3344672134215625d-1*rhoa
     7         *rhob*t11*t18*t36*wght-3.3344672134215625d-1*rhoa*rhob*t1
     8         6*t18*t19*wght-3.3344672134215625d-1*rhob*t11*t16*t18*wgh
     9         t-2.0d+0*t16*t39-t18*t37)+t10*t29*(-1.0003401640264689d+0
     :         *rhoa*rhob*t11*t16*t36*wght-5.001700820132344d-1*rhoa*rho
     ;         b*t19*t30*wght-5.001700820132344d-1*rhob*t11*t30*wght-3.0
     <         d+0*t16*t37)+5.001700820132344d-1*rhoa*rhob*t11*t18**3*t5
     =         7*t8*wght+2.0006803280529378d+0*rhoa*rhob*t10*t11*t16**3*
     >         t58*wght+1.5005102460397032d+0*rhoa*rhob*t11*t17*t18*t29*
     ?         t30*wght+1.0003401640264689d+0*rhoa*rhob*t11*t12*t16*t27*
     @         t28*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = t17*t8*(8.336168033553906d-2*rhoa*rh
     1         ob*t11*(2.3991838839497967d+1*t1*t2*t3*t5*t64-5.598123671
     2         158411d+0*t1*t22*t31*t59+7.997279613165988d+0*t1*t15*t3*t
     3         43*t5+7.997279613165988d+0*t1*t2*t21*t36*t5-2.89498340727
     4         54554d-1*t1*t21*t33*t5+2.6657598710553293d+0*t1*t15*t16*t
     5         21*t5-1.4883349001661852d+1*t1*t13*t14*t3*t43-1.488334900
     6         1661852d+1*t1*t14*t2*t22*t36+5.387713133546839d-1*t1*t14*
     7         t22*t33+1.002681146166278d+0*t1*t21*t31*t32+1.84657857188
     8         1605d+1*t1*t13*t16*t22*t32-4.9611163338872843d+0*t1*t14*t
     9         15*t16*t22-4.9611163338872843d+0*t1*t13*t14*t16*t21)*wght
     :         +8.336168033553906d-2*rhoa*rhob*t19*t45*wght+8.3361680335
     ;         53906d-2*rhob*t11*t45*wght+8.336168033553906d-2*rhoa*rhob
     <         *t18*t40*wght+8.336168033553906d-2*rhoa*rhob*t25*t38*wght
     =         +8.336168033553906d-2*rhoa*t11*t38*wght+8.336168033553906
     >         d-2*rhob*t18*t25*wght+8.336168033553906d-2*rhoa*t18*t19*w
     ?         ght+8.336168033553906d-2*t11*t18*wght-t20*t45-t18*t41-t24
     @         *t35)+t10*t12*(1.6672336067107812d-1*rhoa*rhob*t11*t64*wg
     1         ht+1.6672336067107812d-1*rhoa*rhob*t19*t43*wght+1.6672336
     2         067107812d-1*rhob*t11*t43*wght+1.6672336067107812d-1*rhoa
     3         *rhob*t16*t40*wght+1.6672336067107812d-1*rhoa*rhob*t25*t3
     4         6*wght+1.6672336067107812d-1*rhoa*t11*t36*wght+1.66723360
     5         67107812d-1*rhob*t16*t25*wght+1.6672336067107812d-1*rhoa*
     6         t16*t19*wght+1.6672336067107812d-1*t11*t16*wght-2.0d+0*t2
     7         0*t43-2.0d+0*t16*t41-2.0d+0*t23*t35)+t10*t8*(-8.336168033
     8         553906d-2*rhoa*rhob*(-7.547842390376116d+0*t1*t22*t31*t59
     9         -3.9032504039603666d-1*t1*t21*t33*t5+7.264149912600285d-1
     :         *t1*t14*t22*t33+1.3518956892745262d+0*t1*t21*t31*t32)*wgh
     ;         t-1.6672336067107812d-1*rhob*t40*wght-8.336168033553906d-
     <         2*rhoa*t34*wght-1.6672336067107812d-1*t19*wght)+t27*t8*(-
     =         3.3344672134215625d-1*rhoa*rhob*t11*t18*t45*wght-1.667233
     >         6067107812d-1*rhoa*rhob*t25*t28*wght-1.6672336067107812d-
     ?         1*rhoa*t11*t28*wght-2*t24*t39)+t12*t17*(-3.33446721342156
     @         25d-1*rhoa*rhob*t11*t16*t45*wght-3.3344672134215625d-1*rh
     1         oa*rhob*t11*t18*t43*wght-3.3344672134215625d-1*rhoa*rhob*
     2         t16*t18*t25*wght-3.3344672134215625d-1*rhoa*t11*t16*t18*w
     3         ght-2.0d+0*t23*t39-t24*t37)+t10*t29*(-1.0003401640264689d
     4         +0*rhoa*rhob*t11*t16*t43*wght-5.001700820132344d-1*rhoa*r
     5         hob*t25*t30*wght-5.001700820132344d-1*rhoa*t11*t30*wght-3
     6         .0d+0*t23*t37)+t17*t29*(5.001700820132344d-1*rhoa*rhob*t1
     7         1*t24*t30*wght+1.0003401640264689d+0*rhoa*rhob*t11*t16*t1
     8         8*t23*wght)+t12*t27*(3.3344672134215625d-1*rhoa*rhob*t11*
     9         t23*t28*wght+6.668934426843125d-1*rhoa*rhob*t11*t16*t18*t
     :         24*wght)+5.001700820132344d-1*rhoa*rhob*t11*t24*t28*t57*t
     ;         8*wght+2.0006803280529378d+0*rhoa*rhob*t10*t11*t23*t30*t5
     <         8*wght+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = t17*t8*(8.336168033553906d-2*rhoa*rh
     1         ob*t11*(2.3991838839497967d+1*t1*t2*t3*t5*t65-5.598123671
     2         158411d+0*t1*t13*t50*t59+1.846578571881605d+1*t1*t16*t2*t
     3         32*t50+1.002681146166278d+0*t1*t15*t32*t50-5.331519742110
     4         659d+0*t1*t16*t2*t49*t5-2.8949834072754554d-1*t1*t15*t49*
     5         t5+1.5994559226331975d+1*t1*t2*t21*t43*t5+5.3877131335468
     6         39d-1*t1*t13*t14*t49-2.9766698003323705d+1*t1*t14*t2*t22*
     7         t43)*wght+8.336168033553906d-2*rhoa*rhob*t18*t51*wght+1.6
     8         672336067107812d-1*rhoa*rhob*t25*t45*wght+1.6672336067107
     9         812d-1*rhoa*t11*t45*wght+1.6672336067107812d-1*rhoa*t18*t
     :         25*wght-t20*t55-2*t24*t41)+t10*t12*(1.6672336067107812d-1
     ;         *rhoa*rhob*t11*t65*wght+1.6672336067107812d-1*rhoa*rhob*t
     <         16*t51*wght+3.3344672134215625d-1*rhoa*rhob*t25*t43*wght+
     =         3.3344672134215625d-1*rhoa*t11*t43*wght+3.334467213421562
     >         5d-1*rhoa*t16*t25*wght-2.0d+0*t20*t53-4.0d+0*t23*t41)+t10
     ?         *t8*(-8.336168033553906d-2*rhoa*rhob*(-7.547842390376116d
     @         +0*t1*t13*t50*t59+1.3518956892745262d+0*t1*t15*t32*t50-3.
     1         9032504039603666d-1*t1*t15*t49*t5+7.264149912600285d-1*t1
     2         *t13*t14*t49)*wght-8.336168033553906d-2*rhob*t51*wght-1.6
     3         672336067107812d-1*rhoa*t40*wght-1.6672336067107812d-1*t2
     4         5*wght)+t27*t8*(-1.6672336067107812d-1*rhoa*rhob*t11*t18*
     5         t55*wght-1.6672336067107812d-1*rhoa*rhob*t11*t24*t45*wght
     6         -1.6672336067107812d-1*rhoa*rhob*t18*t24*t25*wght-1.66723
     7         36067107812d-1*rhoa*t11*t18*t24*wght-2*t24*t46)+t12*t17*(
     8         -1.6672336067107812d-1*rhoa*rhob*t11*t16*t55*wght-1.66723
     9         36067107812d-1*rhoa*rhob*t11*t18*t53*wght-1.6672336067107
     :         812d-1*rhoa*rhob*t11*t23*t45*wght-1.6672336067107812d-1*r
     ;         hoa*rhob*t11*t24*t43*wght-1.6672336067107812d-1*rhoa*rhob
     <         *t16*t24*t25*wght-1.6672336067107812d-1*rhoa*rhob*t18*t23
     =         *t25*wght-1.6672336067107812d-1*rhoa*t11*t16*t24*wght-1.6
     >         672336067107812d-1*rhoa*t11*t18*t23*wght-2.0d+0*t23*t46-t
     ?         24*t44)+t10*t29*(-5.001700820132344d-1*rhoa*rhob*t11*t16*
     @         t53*wght-5.001700820132344d-1*rhoa*rhob*t11*t23*t43*wght-
     1         5.001700820132344d-1*rhoa*rhob*t16*t23*t25*wght-5.0017008
     2         20132344d-1*rhoa*t11*t16*t23*wght-3.0d+0*t23*t44)+t12*t27
     3         *(3.3344672134215625d-1*rhoa*rhob*t11*t18*t23*t24*wght-2*
     4         t24*t42)+t17*t29*(5.001700820132344d-1*rhoa*rhob*t11*t16*
     5         t23*t24*wght-3.0d+0*t23*t42)+5.001700820132344d-1*rhoa*rh
     6         ob*t11*t18*t47*t57*t8*wght+2.0006803280529378d+0*rhoa*rho
     7         b*t10*t11*t16*t48*t58*wght+Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = t17*t8*(8.336168033553906d-2*rhoa*rh
     1         ob*t11*(2.3991838839497967d+1*t1*t2*t3*t5*t69-3.008043438
     2         498834d+0*t1*t2*t32*t68+5.387713133546839d-1*t1*t14*t2*t6
     3         7+1.4474917036377277d+0*t1*t2*t5*t66-5.598123671158411d+0
     4         *t1*t2*t49*t59+1.5994559226331975d+1*t1*t2*t21*t5*t53-2.9
     5         766698003323705d+1*t1*t14*t2*t22*t53+1.846578571881605d+1
     6         *t1*t2*t23*t32*t50-5.331519742110659d+0*t1*t2*t23*t49*t5)
     7         *wght+1.6672336067107812d-1*rhoa*rhob*t25*t55*wght+1.6672
     8         336067107812d-1*rhoa*t11*t55*wght+8.336168033553906d-2*rh
     9         oa*rhob*t24*t51*wght+1.6672336067107812d-1*rhoa*t24*t25*w
     :         ght-t26*t55-2*t24*t52)+t10*t12*(1.6672336067107812d-1*rho
     ;         a*rhob*t11*t69*wght+3.3344672134215625d-1*rhoa*rhob*t25*t
     <         53*wght+3.3344672134215625d-1*rhoa*t11*t53*wght+1.6672336
     =         067107812d-1*rhoa*rhob*t23*t51*wght+3.3344672134215625d-1
     >         *rhoa*t23*t25*wght-2.0d+0*t26*t53-4.0d+0*t23*t52)+t10*t8*
     ?         (-8.336168033553906d-2*rhoa*rhob*(-4.0556870678235785d+0*
     @         t1*t2*t32*t68+7.264149912600285d-1*t1*t14*t2*t67+1.951625
     1         2019801836d+0*t1*t2*t5*t66-7.547842390376116d+0*t1*t2*t49
     2         *t59)*wght-2.500850410066172d-1*rhoa*t51*wght)+t27*t8*(-3
     3         .3344672134215625d-1*rhoa*rhob*t11*t24*t55*wght-1.6672336
     4         067107812d-1*rhoa*rhob*t25*t47*wght-1.6672336067107812d-1
     5         *rhoa*t11*t47*wght-2*t24*t56)+t12*t17*(-3.334467213421562
     6         5d-1*rhoa*rhob*t11*t23*t55*wght-3.3344672134215625d-1*rho
     7         a*rhob*t11*t24*t53*wght-3.3344672134215625d-1*rhoa*rhob*t
     8         23*t24*t25*wght-3.3344672134215625d-1*rhoa*t11*t23*t24*wg
     9         ht-2.0d+0*t23*t56-t24*t54)+t10*t29*(-1.0003401640264689d+
     :         0*rhoa*rhob*t11*t23*t53*wght-5.001700820132344d-1*rhoa*rh
     ;         ob*t25*t48*wght-5.001700820132344d-1*rhoa*t11*t48*wght-3.
     <         0d+0*t23*t54)+5.001700820132344d-1*rhoa*rhob*t11*t24**3*t
     =         57*t8*wght+2.0006803280529378d+0*rhoa*rhob*t10*t11*t23**3
     >         *t58*wght+1.5005102460397032d+0*rhoa*rhob*t11*t17*t24*t29
     ?         *t48*wght+1.0003401640264689d+0*rhoa*rhob*t11*t12*t23*t27
     @         *t47*wght+Amat3(iq,D3_RB_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
