C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_dldf.F 
C> The nwxcm_c_dldf functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\sqrt{\rho_\alpha}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_5} &=& -6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_4}+2.225569421150687\,{\it t_3}+0.8004286349993635\,{\it t_2}+0.1897004325747558\,{\it t_1}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_3}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_4}+3.844746237447211\,{\it t_3}+1.644733775567609\,{\it t_2}+0.2405871291288192\,{\it t_1}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_3}+1.0\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\sqrt{\rho_\beta}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\beta^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_\beta^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{10}} &=& -6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_9}+2.225569421150687\,{\it t_8}+0.8004286349993635\,{\it t_7}+0.1897004325747558\,{\it t_6}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_8}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_9}+3.844746237447211\,{\it t_8}+1.644733775567609\,{\it t_7}+0.2405871291288192\,{\it t_6}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_8}+1.0\right)\\\\ 
C>   {\it t_{11}} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}^4\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{11}}^4}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{17}} &=& \left({\it t_{13}}\,{\it t_{16}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{13}}\,{\it t_{16}}\right)^{{{4}\over{3}}}-2.0\\\\ 
C>   {\it t_{18}} &=& -1.923661050931536\,{\it t_{14}}\,{\it t_{15}}\,{\it t_{17}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{{\it t_{11}}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\sqrt{{\it t_{11}}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{11}}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{24}} &=& {\it t_{23}}\,\sigma_{\beta\beta}+{\it t_{22}}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{0.0031\,{\it t_{24}}+1.0}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{0.06\,{\it t_{22}}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{0.06\,{\it t_{23}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\sqrt{\rho_s}}}\\\\ 
C>   {\it t_{30}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{32}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{33}} &=& {{1}\over{0.06\,{\it t_{32}}\,\sigma_{ss}+1.0}}\\\\ 
C>   f &=& {{0.5\,{\it t_{10}}\,\rho_\beta\,\left(0.06\,{\it t_{23}}\,\sigma_{\beta\beta}\,{\it t_{27}}\,\left(0.133402758\,{\it t_{23}}\,\sigma_{\beta\beta}\,{\it t_{27}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_\beta-{{0.25\,\sigma_{\beta\beta}}\over{\rho_\beta}}\right)}\over{\tau_\beta}}+{{0.5\,{\it t_5}\,\rho_\alpha\,\left(0.06\,{\it t_{22}}\,\sigma_{\alpha\alpha}\,{\it t_{26}}\,\left(0.133402758\,{\it t_{22}}\,\sigma_{\alpha\alpha}\,{\it t_{26}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_\alpha-{{0.25\,\sigma_{\alpha\alpha}}\over{\rho_\alpha}}\right)}\over{\tau_\alpha}}+1.0\,\left({\it t_{11}}\,\left(-0.0598079686361971\,{\it t_{14}}\,{\it t_{15}}\,\left(0.1274696188700087\,{\it t_{12}}+1.0\right)\,{\it t_{17}}\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{21}}+3.844746237447211\,{\it t_{12}}+1.644733775567609\,{\it t_{20}}+0.2405871291288192\,{\it t_{19}}}}+1.0\right)+0.01975163132168077\,\left(0.06901399211255825\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.923661050931536\,{\it t_{17}}\right)\,\log \left({{29.60881320326807}\over{8.15741470348764\,{\it t_{21}}+2.247591863577616\,{\it t_{12}}+0.4300972471276643\,{\it t_{20}}+0.1911512595127337\,{\it t_{19}}}}+1.0\right)-0.0621814\,\left(0.1325688999052018\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.0\right)\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{21}}+2.225569421150687\,{\it t_{12}}+0.8004286349993635\,{\it t_{20}}+0.1897004325747558\,{\it t_{19}}}}+1.0\right)\right)-{\it t_{10}}\,\rho_\beta-{\it t_5}\,\rho_\alpha\right)\,\left(0.0031\,{\it t_{24}}\,{\it t_{25}}\,\left(5.9515308-0.03459689187\,{\it t_{24}}\,{\it t_{25}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.5\,\left(-6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{31}}+2.225569421150687\,{\it t_{30}}+0.8004286349993635\,{\it t_{29}}+0.1897004325747558\,{\it t_{28}}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_{30}}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{31}}+3.844746237447211\,{\it t_{30}}+1.644733775567609\,{\it t_{29}}+0.2405871291288192\,{\it t_{28}}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_{30}}+1.0\right)\right)\,\rho_s\,\left(0.06\,{\it t_{32}}\,\sigma_{ss}\,{\it t_{33}}\,\left(0.133402758\,{\it t_{32}}\,\sigma_{ss}\,{\it t_{33}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_s-{{0.25\,\sigma_{ss}}\over{\rho_s}}\right)}\over{\tau_s}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
      subroutine nwxcm_c_dldf(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.2d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.9515308d+0-6.919378374d-2*gammaaa*t1*t3
            t5 = 6.2d-3*gammaaa*t1*t3*t4+1.0d+0
            t6 = 1/rhoa**6.666666666666666d-1
            t7 = rhoa**5.0d-1
            t8 = 1/t7
            t9 = 1/rhoa**3.333333333333333d-1
            t10 = 1/rhoa**1.6666666666666666d-1
            t11 = 3.8447462374472113d+0*t9+1.644733775567609d+0*t8+2.405
     1         8712912881916d-1*t6+1.1120374863094686d+1*t10
            t12 = 3.216395899738507d+1/t11+1.0d+0
            t13 = log(t12)
            t14 = 1.274696188700087d-1*t9+1.0d+0
            t15 = 2.2255694211506874d+0*t9+8.004286349993635d-1*t8+1.897
     1         0043257475583d-1*t6+5.982550435771079d+0*t10
            t16 = 1.6081979498692536d+1/t15+1.0d+0
            t17 = log(t16)
            t18 = 1.3256889990520176d-1*t9+1.0d+0
            t19 = -6.903522198342671d-18*t17*t18-3.1090699999999993d-2*t
     1         13*t14
            t20 = 1.766435620181429d+0*t9+5.659885156639418d-1*t8+1.1950
     1         3784087549d-1*t6+5.329846514438383d+0*t10
            t21 = 1.6081979498692536d+1/t20+1.0d+0
            t22 = log(t21)
            t23 = 1.0522000558389212d-1*t9+1.0d+0
            t24 = -1.243628d-1*rhoa*t22*t23-2*rhoa*t19
            t25 = 6.0d-2*gammaaa*t1+1.0d+0
            t26 = 1/t25
            t27 = 1.3340275799999998d-1*gammaaa*t1*t26-2.5960897d+0
            t28 = 6.0d-2*gammaaa*t1*t26*t27+1.0d+0
            t29 = 1/taua
            t30 = 1/rhoa
            t31 = 2.0d+0*taua-2.5d-1*gammaaa*t30
            t32 = 1/rhoa**1.6666666666666669d+0
            t33 = 1/t7**3
            t34 = 1/rhoa**1.3333333333333333d+0
            t35 = 1/rhoa**1.1666666666666667d+0
            t36 = 1.1102230246251565d-16*t18*(-9.970917392951799d-1*t35-
     1         7.418564737168958d-1*t34-4.002143174996817d-1*t33-1.26466
     2         9550498372d-1*t32)/(t15**2*t16)+1.0d+0*t14*(-1.8533958105
     3         157808d+0*t35-1.2815820791490704d+0*t34-8.223668877838045
     4         d-1*t33-1.6039141941921276d-1*t32)/(t11**2*t12)+3.0506411
     5         443514264d-19*t17*t34+1.3210398931339265d-3*t13*t34
            t37 = gammaaa**2
            t38 = 1/t2**2
            t39 = 1/rhoa**6.333333333333333d+0
            t40 = 1/rhoa**3.6666666666666664d+0
            t41 = 1/t25**2
            t42 = 1/rhoa**5.333333333333333d+0
            t43 = (t29*(1.0d+0*rhoa*t19*t31*(6.0d-2*gammaaa*t1*t26*(1.33
     1         40275799999998d-1*t1*t26-8.004165479999998d-3*gammaaa*t41
     2         *t42)-3.6d-3*gammaaa*t27*t41*t42+6.0d-2*t1*t26*t27)-2.5d-
     3         1*t19*t28)+1.0d+0*t24*(6.2d-3*gammaaa*t1*t3*(4.2900145918
     4         8d-4*gammaaa*t38*t42-6.919378374d-2*t1*t3)-3.843999999999
     5         9996d-5*gammaaa*t38*t4*t42+6.2d-3*t1*t3*t4))*wght
            fnc(iq) = (1.0d+0*t24*t5+1.0d+0*rhoa*t19*t28*t29*t31)*wght+f
     1         nc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t5*(4.361818170142819d-3*t22*t9-2*r
     1         hoa*t36+2.0d+0*rhoa*t23*(-8.883077524063971d-1*t35-5.8881
     2         18733938097d-1*t34-2.829942578319709d-1*t33-7.96691893916
     3         9932d-2*t32)/(t20**2*t21)-1.243628d-1*t22*t23-2*t19)+t29*
     4         (1.0d+0*rhoa*t19*t31*(6.0d-2*gammaaa*t1*t26*(2.1344441279
     5         999993d-2*t37*t39*t41-3.557406879999999d-1*gammaaa*t26*t4
     6         0)+9.599999999999997d-3*t27*t37*t39*t41-1.599999999999999
     7         8d-1*gammaaa*t26*t27*t40)+1.0d+0*rhoa*t28*t31*t36+1.0d+0*
     8         t19*t28*t31)+1.0d+0*t24*(6.2d-3*gammaaa*t1*t3*(1.84516756
     9         64d-1*gammaaa*t3*t40-1.1440038911679998d-3*t37*t38*t39)-1
     :         .653333333333333d-2*gammaaa*t3*t4*t40+1.0250666666666663d
     ;         -4*t37*t38*t39*t4)+2.5d-1*gammaaa*t19*t28*t29*t30)*wght+A
     <         mat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t43+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t43+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = (2.0d+0*rhoa*t19*t28*t29-1.0d+0*rhoa*t19*t2
     1         8*t31/taua**2)*wght+Amat(iq,D1_TA)
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = gammabb*t2+gammaaa*t1
            t4 = 3.1d-3*t3+1.0d+0
            t5 = 1/t4
            t6 = 5.9515308d+0-3.459689187d-2*t3*t5
            t7 = 3.1d-3*t3*t5*t6+1.0d+0
            t8 = 1/rhoa**6.666666666666666d-1
            t9 = rhoa**5.0d-1
            t10 = 1/t9
            t11 = 1/rhoa**3.333333333333333d-1
            t12 = 1/rhoa**1.6666666666666666d-1
            t13 = 2.4058712912881916d-1*t8+1.1120374863094686d+1*t12+3.8
     1         447462374472113d+0*t11+1.644733775567609d+0*t10
            t14 = 3.216395899738507d+1/t13+1.0d+0
            t15 = log(t14)
            t16 = 1.274696188700087d-1*t11+1.0d+0
            t17 = 1.8970043257475583d-1*t8+5.982550435771079d+0*t12+2.22
     1         55694211506874d+0*t11+8.004286349993635d-1*t10
            t18 = 1.6081979498692536d+1/t17+1.0d+0
            t19 = log(t18)
            t20 = 1.3256889990520176d-1*t11+1.0d+0
            t21 = -6.903522198342671d-18*t19*t20-3.1090699999999993d-2*t
     1         15*t16
            t22 = 1/rhob**6.666666666666666d-1
            t23 = rhob**5.0d-1
            t24 = 1/t23
            t25 = 1/rhob**3.333333333333333d-1
            t26 = 1/rhob**1.6666666666666666d-1
            t27 = 1.1120374863094686d+1*t26+3.8447462374472113d+0*t25+1.
     1         644733775567609d+0*t24+2.4058712912881916d-1*t22
            t28 = 3.216395899738507d+1/t27+1.0d+0
            t29 = log(t28)
            t30 = 1.274696188700087d-1*t25+1.0d+0
            t31 = 5.982550435771079d+0*t26+2.2255694211506874d+0*t25+8.0
     1         04286349993635d-1*t24+1.8970043257475583d-1*t22
            t32 = 1.6081979498692536d+1/t31+1.0d+0
            t33 = log(t32)
            t34 = 1.3256889990520176d-1*t25+1.0d+0
            t35 = -6.903522198342671d-18*t33*t34-3.1090699999999993d-2*t
     1         29*t30
            t36 = rhob+rhoa
            t37 = 1/t36**3.333333333333333d-1
            t38 = 1.3256889990520176d-1*t37+1.0d+0
            t39 = rhoa-rhob
            t40 = t39**4
            t41 = 1/t36**4
            t42 = 1/t36
            t43 = 1.0d+0-t39*t42
            t44 = t39*t42+1.0d+0
            t45 = t44**1.3333333333333333d+0+t43**1.3333333333333333d+0-
     1         2.0d+0
            t46 = -1.923661050931536d+0*t40*t41*t45
            t47 = t46+1.0d+0
            t48 = 1/t36**6.666666666666666d-1
            t49 = t36**5.0d-1
            t50 = 1/t49
            t51 = 1/t36**1.6666666666666666d-1
            t52 = 5.982550435771079d+0*t51+8.004286349993635d-1*t50+1.89
     1         70043257475583d-1*t48+2.2255694211506874d+0*t37
            t53 = 1.6081979498692536d+1/t52+1.0d+0
            t54 = log(t53)
            t55 = -6.21814d-2*t38*t47*t54
            t56 = 6.901399211255825d-2*t37+1.0d+0
            t57 = t46+1.923661050931536d+0*t45
            t58 = 8.157414703487641d+0*t51+4.3009724712766434d-1*t50+1.9
     1         115125951273373d-1*t48+2.247591863577616d+0*t37
            t59 = 2.9608813203268075d+1/t58+1.0d+0
            t60 = log(t59)
            t61 = 1.9751631321680768d-2*t56*t57*t60
            t62 = 1.274696188700087d-1*t37+1.0d+0
            t63 = 1.1120374863094686d+1*t51+1.644733775567609d+0*t50+2.4
     1         058712912881916d-1*t48+3.8447462374472113d+0*t37
            t64 = 3.216395899738507d+1/t63+1.0d+0
            t65 = log(t64)
            t66 = -5.980796863619711d-2*t40*t41*t45*t62*t65
            t67 = t36*(t66+t61+t55)-rhob*t35-rhoa*t21
            t68 = 6.0d-2*gammaaa*t1+1.0d+0
            t69 = 1/t68
            t70 = 1.3340275799999998d-1*gammaaa*t1*t69-2.5960897d+0
            t71 = 6.0d-2*gammaaa*t1*t69*t70+1.0d+0
            t72 = 1/taua
            t73 = 1/rhoa
            t74 = 2.0d+0*taua-2.5d-1*gammaaa*t73
            t75 = 6.0d-2*gammabb*t2+1.0d+0
            t76 = 1/t75
            t77 = 1.3340275799999998d-1*gammabb*t2*t76-2.5960897d+0
            t78 = 6.0d-2*gammabb*t2*t76*t77+1.0d+0
            t79 = 1/taub
            t80 = 1/rhob
            t81 = 2.0d+0*taub-2.5d-1*gammabb*t80
            t82 = 1/rhoa**1.6666666666666669d+0
            t83 = 1/t9**3
            t84 = 1/rhoa**1.3333333333333333d+0
            t85 = 1/rhoa**1.1666666666666667d+0
            t86 = 1.1102230246251565d-16*t20*(-9.970917392951799d-1*t85-
     1         7.418564737168958d-1*t84-4.002143174996817d-1*t83-1.26466
     2         9550498372d-1*t82)/(t17**2*t18)+1.0d+0*t16*(-1.8533958105
     3         157808d+0*t85-1.2815820791490704d+0*t84-8.223668877838045
     4         d-1*t83-1.6039141941921276d-1*t82)/(t13**2*t14)+3.0506411
     5         443514264d-19*t19*t84+1.3210398931339265d-3*t15*t84
            t87 = 1/t36**1.6666666666666669d+0
            t88 = 1/t49**3
            t89 = 1/t36**1.3333333333333333d+0
            t90 = 1/t36**1.1666666666666667d+0
            t91 = 1.0d+0*t38*t47*(-9.970917392951799d-1*t90-7.4185647371
     1         68958d-1*t89-4.002143174996817d-1*t88-1.264669550498372d-
     2         1*t87)/(t52**2*t53)
            t92 = -5.848223622634647d-1*t56*t57*(-1.3595691172479402d+0*
     1         t90-7.491972878592053d-1*t89-2.1504862356383217d-1*t88-1.
     2         2743417300848914d-1*t87)/(t58**2*t59)
            t93 = 1.923661050931536d+0*t40*t41*t45*t62*(-1.8533958105157
     1         808d+0*t90-1.2815820791490704d+0*t89-8.223668877838045d-1
     2         *t88-1.6039141941921276d-1*t87)/(t63**2*t64)
            t94 = 1/t36**2
            t95 = t39*t94
            t96 = -t42
            t97 = t43**3.333333333333333d-1
            t98 = -t39*t94
            t99 = t44**3.333333333333333d-1
            t100 = 1.3333333333333333d+0*(t98+t42)*t99+1.333333333333333
     1         3d+0*(t96+t95)*t97
            t101 = -1.923661050931536d+0*t100*t40*t41
            t102 = 1/t36**5
            t103 = 7.694644203726144d+0*t102*t40*t45
            t104 = t39**3
            t105 = -7.694644203726144d+0*t104*t41*t45
            t106 = 2.747773264188437d-3*t47*t54*t89
            t107 = -4.5437964274821163d-4*t57*t60*t89
            t108 = 2.541232989148493d-3*t40*t45*t65/t36**5.3333333333333
     1         33d+0
            t109 = 2.3923187454478842d-1*t102*t40*t45*t62*t65
            t110 = 1/rhoa**3.6666666666666664d+0
            t111 = 1/t4**2
            t112 = gammaaa**2
            t113 = 1/t68**2
            t114 = 1/rhoa**6.333333333333333d+0
            t115 = 1/rhob**1.6666666666666669d+0
            t116 = 1/t23**3
            t117 = 1/rhob**1.3333333333333333d+0
            t118 = 1/rhob**1.1666666666666667d+0
            t119 = 1.1102230246251565d-16*(-9.970917392951799d-1*t118-7.
     1         418564737168958d-1*t117-4.002143174996817d-1*t116-1.26466
     2         9550498372d-1*t115)*t34/(t31**2*t32)+3.0506411443514264d-
     3         19*t117*t33+1.0d+0*(-1.8533958105157808d+0*t118-1.2815820
     4         791490704d+0*t117-8.223668877838045d-1*t116-1.60391419419
     5         21276d-1*t115)*t30/(t27**2*t28)+1.3210398931339265d-3*t11
     6         7*t29
            t120 = 1.3333333333333333d+0*(t98+t96)*t99+1.333333333333333
     1         3d+0*(t95+t42)*t97
            t121 = -1.923661050931536d+0*t120*t40*t41
            t122 = 7.694644203726144d+0*t104*t41*t45
            t123 = 1/rhob**3.6666666666666664d+0
            t124 = gammabb**2
            t125 = 1/t75**2
            t126 = 1/rhob**6.333333333333333d+0
            t127 = 1/rhoa**5.333333333333333d+0
            t128 = 1/rhob**5.333333333333333d+0
            fnc(iq) = (5.0d-1*rhob*t35*t78*t79*t81+5.0d-1*rhoa*t21*t71*t
     1         72*t74+1.0d+0*t67*t7)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t7*(t36*(t93+t92+t91-2.392318745447
     1         8842d-1*t104*t41*t45*t62*t65-5.980796863619711d-2*t100*t4
     2         0*t41*t62*t65+1.9751631321680768d-2*(t105+t103+t101+1.923
     3         661050931536d+0*t100)*t56*t60-6.21814d-2*(t105+t103+t101)
     4         *t38*t54+t109+t108+t107+t106)-rhoa*t86+t66+t61+t55+6.9035
     5         22198342671d-18*t19*t20+3.1090699999999993d-2*t15*t16)+t7
     6         2*(5.0d-1*rhoa*t71*t74*t86+5.0d-1*t21*t71*t74+5.0d-1*rhoa
     7         *t21*(-1.5999999999999998d-1*gammaaa*t110*t69*t70+9.59999
     8         9999999997d-3*t112*t113*t114*t70+6.0d-2*gammaaa*t1*t69*(2
     9         .1344441279999993d-2*t112*t113*t114-3.557406879999999d-1*
     :         gammaaa*t110*t69))*t74)+1.25d-1*gammaaa*t21*t71*t72*t73+1
     ;         .0d+0*(-8.266666666666665d-3*gammaaa*t110*t5*t6+2.5626666
     <         66666666d-5*gammaaa*t110*t111*t3*t6+3.1d-3*t3*t5*(9.22583
     =         7832d-2*gammaaa*t110*t5-2.86000972792d-4*gammaaa*t110*t11
     >         1*t3))*t67)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t7*(t36*(t93+t92+t91+2.392318745447
     1         8842d-1*t104*t41*t45*t62*t65-5.980796863619711d-2*t120*t4
     2         0*t41*t62*t65+1.9751631321680768d-2*(t122+t121+1.92366105
     3         0931536d+0*t120+t103)*t56*t60-6.21814d-2*(t122+t121+t103)
     4         *t38*t54+t109+t108+t107+t106)+t66+t61+t55+6.9035221983426
     5         71d-18*t33*t34+3.1090699999999993d-2*t29*t30-rhob*t119)+t
     6         79*(5.0d-1*t35*t78*t81+5.0d-1*rhob*t119*t78*t81+5.0d-1*rh
     7         ob*t35*(-1.5999999999999998d-1*gammabb*t123*t76*t77+9.599
     8         999999999997d-3*t124*t125*t126*t77+6.0d-2*gammabb*t2*t76*
     9         (2.1344441279999993d-2*t124*t125*t126-3.557406879999999d-
     :         1*gammabb*t123*t76))*t81)+1.25d-1*gammabb*t35*t78*t79*t80
     ;         +1.0d+0*(-8.266666666666665d-3*gammabb*t123*t5*t6+2.56266
     <         6666666666d-5*gammabb*t111*t123*t3*t6+3.1d-3*t3*t5*(9.225
     =         837832d-2*gammabb*t123*t5-2.86000972792d-4*gammabb*t111*t
     >         123*t3))*t67)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t72*(5.0d-1*rhoa*t21*(6.0d-2*t1*t69*t70-3
     1         .6d-3*gammaaa*t113*t127*t70+6.0d-2*gammaaa*t1*t69*(1.3340
     2         275799999998d-1*t1*t69-8.004165479999998d-3*gammaaa*t113*
     3         t127))*t74-1.25d-1*t21*t71)+1.0d+0*(3.1d-3*t1*t5*t6-9.61d
     4         -6*t1*t111*t3*t6+3.1d-3*t3*t5*(1.07250364797d-4*t1*t111*t
     5         3-3.459689187d-2*t1*t5))*t67)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t79*(5.0d-1*rhob*t35*(6.0d-2*t2*t76*t77-3
     1         .6d-3*gammabb*t125*t128*t77+6.0d-2*gammabb*t2*t76*(1.3340
     2         275799999998d-1*t2*t76-8.004165479999998d-3*gammabb*t125*
     3         t128))*t81-1.25d-1*t35*t78)+1.0d+0*(3.1d-3*t2*t5*t6-9.61d
     4         -6*t111*t2*t3*t6+3.1d-3*t3*t5*(1.07250364797d-4*t111*t2*t
     5         3-3.459689187d-2*t2*t5))*t67)*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TA) = (1.0d+0*rhoa*t21*t71*t72-5.0d-1*rhoa*t21*t7
     1         1*t74/taua**2)*wght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = (1.0d+0*rhob*t35*t78*t79-5.0d-1*rhob*t35*t7
     1         8*t81/taub**2)*wght+Mmat(iq,D1_TB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-2*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 1.3340275799999998d-1*gammaaa*t1*t3-2.5960897d+0
            t5 = 6.0d-2*gammaaa*t1*t3*t4+1.0d+0
            t6 = 1/rhoa**6.666666666666666d-1
            t7 = rhoa**5.0d-1
            t8 = 1/t7
            t9 = 1/rhoa**3.333333333333333d-1
            t10 = 1/rhoa**1.6666666666666666d-1
            t11 = 3.8447462374472113d+0*t9+1.644733775567609d+0*t8+2.405
     1         8712912881916d-1*t6+1.1120374863094686d+1*t10
            t12 = 3.216395899738507d+1/t11+1.0d+0
            t13 = log(t12)
            t14 = 1.274696188700087d-1*t9+1.0d+0
            t15 = 2.2255694211506874d+0*t9+8.004286349993635d-1*t8+1.897
     1         0043257475583d-1*t6+5.982550435771079d+0*t10
            t16 = 1.6081979498692536d+1/t15+1.0d+0
            t17 = log(t16)
            t18 = 1.3256889990520176d-1*t9+1.0d+0
            t19 = -6.903522198342671d-18*t17*t18-3.1090699999999993d-2*t
     1         13*t14
            t20 = 1/taua
            t21 = 1/rhoa
            t22 = 2.0d+0*taua-2.5d-1*gammaaa*t21
            t23 = 1/rhoa**1.6666666666666669d+0
            t24 = 1/t7**3
            t25 = 1/rhoa**1.3333333333333333d+0
            t26 = 1/rhoa**1.1666666666666667d+0
            t27 = gammaaa**2
            t28 = 1/t2**2
            t29 = 1/rhoa**6.333333333333333d+0
            t30 = 1/rhoa**3.6666666666666664d+0
            t31 = 1/rhoa**5.333333333333333d+0
            fnc(iq) = 5.0d-1*rhoa*t19*t20*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = t20*(5.0d-1*rhoa*t22*(1.1102230246251565d-1
     1         6*t18*(-9.970917392951799d-1*t26-7.418564737168958d-1*t25
     2         -4.002143174996817d-1*t24-1.264669550498372d-1*t23)/(t15*
     3         *2*t16)+1.0d+0*t14*(-1.8533958105157808d+0*t26-1.28158207
     4         91490704d+0*t25-8.223668877838045d-1*t24-1.60391419419212
     5         76d-1*t23)/(t11**2*t12)+3.0506411443514264d-19*t17*t25+1.
     6         3210398931339265d-3*t13*t25)*t5*wght+5.0d-1*t19*t22*t5*wg
     7         ht+5.0d-1*rhoa*t19*t22*(-1.5999999999999998d-1*gammaaa*t3
     8         *t30*t4+9.599999999999997d-3*t27*t28*t29*t4+6.0d-2*gammaa
     9         a*t1*t3*(2.1344441279999993d-2*t27*t28*t29-3.557406879999
     :         999d-1*gammaaa*t3*t30))*wght)+1.25d-1*gammaaa*t19*t20*t21
     ;         *t5*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t20*(5.0d-1*rhoa*t19*t22*(-3.6d-3*gammaaa*
     1         t28*t31*t4+6.0d-2*t1*t3*t4+6.0d-2*gammaaa*t1*t3*(1.334027
     2         5799999998d-1*t1*t3-8.004165479999998d-3*gammaaa*t28*t31)
     3         )*wght-1.25d-1*t19*t5*wght)+Cmat(iq,D1_GAA)
            Mmat(iq,D1_TA) = -5.0d-1*rhoa*t19*t22*t5*wght/taua**2+1.0d+0
     1         *rhoa*t19*t20*t5*wght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 6.0d-2*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 1.3340275799999998d-1*gammabb*t1*t3-2.5960897d+0
            t5 = 6.0d-2*gammabb*t1*t3*t4+1.0d+0
            t6 = 1/rhob**6.666666666666666d-1
            t7 = rhob**5.0d-1
            t8 = 1/t7
            t9 = 1/rhob**3.333333333333333d-1
            t10 = 1/rhob**1.6666666666666666d-1
            t11 = 3.8447462374472113d+0*t9+1.644733775567609d+0*t8+2.405
     1         8712912881916d-1*t6+1.1120374863094686d+1*t10
            t12 = 3.216395899738507d+1/t11+1.0d+0
            t13 = log(t12)
            t14 = 1.274696188700087d-1*t9+1.0d+0
            t15 = 2.2255694211506874d+0*t9+8.004286349993635d-1*t8+1.897
     1         0043257475583d-1*t6+5.982550435771079d+0*t10
            t16 = 1.6081979498692536d+1/t15+1.0d+0
            t17 = log(t16)
            t18 = 1.3256889990520176d-1*t9+1.0d+0
            t19 = -6.903522198342671d-18*t17*t18-3.1090699999999993d-2*t
     1         13*t14
            t20 = 1/taub
            t21 = 1/rhob
            t22 = 2.0d+0*taub-2.5d-1*gammabb*t21
            t23 = 1/rhob**1.6666666666666669d+0
            t24 = 1/t7**3
            t25 = 1/rhob**1.3333333333333333d+0
            t26 = 1/rhob**1.1666666666666667d+0
            t27 = gammabb**2
            t28 = 1/t2**2
            t29 = 1/rhob**6.333333333333333d+0
            t30 = 1/rhob**3.6666666666666664d+0
            t31 = 1/rhob**5.333333333333333d+0
            fnc(iq) = 5.0d-1*rhob*t19*t20*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RB) = t20*(5.0d-1*rhob*t22*(1.1102230246251565d-1
     1         6*t18*(-9.970917392951799d-1*t26-7.418564737168958d-1*t25
     2         -4.002143174996817d-1*t24-1.264669550498372d-1*t23)/(t15*
     3         *2*t16)+1.0d+0*t14*(-1.8533958105157808d+0*t26-1.28158207
     4         91490704d+0*t25-8.223668877838045d-1*t24-1.60391419419212
     5         76d-1*t23)/(t11**2*t12)+3.0506411443514264d-19*t17*t25+1.
     6         3210398931339265d-3*t13*t25)*t5*wght+5.0d-1*t19*t22*t5*wg
     7         ht+5.0d-1*rhob*t19*t22*(-1.5999999999999998d-1*gammabb*t3
     8         *t30*t4+9.599999999999997d-3*t27*t28*t29*t4+6.0d-2*gammab
     9         b*t1*t3*(2.1344441279999993d-2*t27*t28*t29-3.557406879999
     :         999d-1*gammabb*t3*t30))*wght)+1.25d-1*gammabb*t19*t20*t21
     ;         *t5*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t20*(5.0d-1*rhob*t19*t22*(-3.6d-3*gammabb*
     1         t28*t31*t4+6.0d-2*t1*t3*t4+6.0d-2*gammabb*t1*t3*(1.334027
     2         5799999998d-1*t1*t3-8.004165479999998d-3*gammabb*t28*t31)
     3         )*wght-1.25d-1*t19*t5*wght)+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TB) = -5.0d-1*rhob*t19*t22*t5*wght/taub**2+1.0d+0
     1         *rhob*t19*t20*t5*wght+Mmat(iq,D1_TB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\sqrt{\rho_\alpha}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_5} &=& -6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_4}+2.225569421150687\,{\it t_3}+0.8004286349993635\,{\it t_2}+0.1897004325747558\,{\it t_1}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_3}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_4}+3.844746237447211\,{\it t_3}+1.644733775567609\,{\it t_2}+0.2405871291288192\,{\it t_1}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_3}+1.0\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\sqrt{\rho_\beta}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\beta^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_\beta^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{10}} &=& -6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_9}+2.225569421150687\,{\it t_8}+0.8004286349993635\,{\it t_7}+0.1897004325747558\,{\it t_6}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_8}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_9}+3.844746237447211\,{\it t_8}+1.644733775567609\,{\it t_7}+0.2405871291288192\,{\it t_6}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_8}+1.0\right)\\\\ 
C>   {\it t_{11}} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}^4\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{11}}^4}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{17}} &=& \left({\it t_{13}}\,{\it t_{16}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{13}}\,{\it t_{16}}\right)^{{{4}\over{3}}}-2.0\\\\ 
C>   {\it t_{18}} &=& -1.923661050931536\,{\it t_{14}}\,{\it t_{15}}\,{\it t_{17}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{{\it t_{11}}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\sqrt{{\it t_{11}}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{11}}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{24}} &=& {\it t_{23}}\,\sigma_{\beta\beta}+{\it t_{22}}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{0.0031\,{\it t_{24}}+1.0}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{0.06\,{\it t_{22}}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{0.06\,{\it t_{23}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\sqrt{\rho_s}}}\\\\ 
C>   {\it t_{30}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{32}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{33}} &=& {{1}\over{0.06\,{\it t_{32}}\,\sigma_{ss}+1.0}}\\\\ 
C>   f &=& {{0.5\,{\it t_{10}}\,\rho_\beta\,\left(0.06\,{\it t_{23}}\,\sigma_{\beta\beta}\,{\it t_{27}}\,\left(0.133402758\,{\it t_{23}}\,\sigma_{\beta\beta}\,{\it t_{27}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_\beta-{{0.25\,\sigma_{\beta\beta}}\over{\rho_\beta}}\right)}\over{\tau_\beta}}+{{0.5\,{\it t_5}\,\rho_\alpha\,\left(0.06\,{\it t_{22}}\,\sigma_{\alpha\alpha}\,{\it t_{26}}\,\left(0.133402758\,{\it t_{22}}\,\sigma_{\alpha\alpha}\,{\it t_{26}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_\alpha-{{0.25\,\sigma_{\alpha\alpha}}\over{\rho_\alpha}}\right)}\over{\tau_\alpha}}+1.0\,\left({\it t_{11}}\,\left(-0.0598079686361971\,{\it t_{14}}\,{\it t_{15}}\,\left(0.1274696188700087\,{\it t_{12}}+1.0\right)\,{\it t_{17}}\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{21}}+3.844746237447211\,{\it t_{12}}+1.644733775567609\,{\it t_{20}}+0.2405871291288192\,{\it t_{19}}}}+1.0\right)+0.01975163132168077\,\left(0.06901399211255825\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.923661050931536\,{\it t_{17}}\right)\,\log \left({{29.60881320326807}\over{8.15741470348764\,{\it t_{21}}+2.247591863577616\,{\it t_{12}}+0.4300972471276643\,{\it t_{20}}+0.1911512595127337\,{\it t_{19}}}}+1.0\right)-0.0621814\,\left(0.1325688999052018\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.0\right)\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{21}}+2.225569421150687\,{\it t_{12}}+0.8004286349993635\,{\it t_{20}}+0.1897004325747558\,{\it t_{19}}}}+1.0\right)\right)-{\it t_{10}}\,\rho_\beta-{\it t_5}\,\rho_\alpha\right)\,\left(0.0031\,{\it t_{24}}\,{\it t_{25}}\,\left(5.9515308-0.03459689187\,{\it t_{24}}\,{\it t_{25}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.5\,\left(-6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{31}}+2.225569421150687\,{\it t_{30}}+0.8004286349993635\,{\it t_{29}}+0.1897004325747558\,{\it t_{28}}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_{30}}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{31}}+3.844746237447211\,{\it t_{30}}+1.644733775567609\,{\it t_{29}}+0.2405871291288192\,{\it t_{28}}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_{30}}+1.0\right)\right)\,\rho_s\,\left(0.06\,{\it t_{32}}\,\sigma_{ss}\,{\it t_{33}}\,\left(0.133402758\,{\it t_{32}}\,\sigma_{ss}\,{\it t_{33}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_s-{{0.25\,\sigma_{ss}}\over{\rho_s}}\right)}\over{\tau_s}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
      subroutine nwxcm_c_dldf_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd derivative wrt tau 
                                             !< and possibly rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
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
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\sqrt{\rho_\alpha}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_5} &=& -6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_4}+2.225569421150687\,{\it t_3}+0.8004286349993635\,{\it t_2}+0.1897004325747558\,{\it t_1}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_3}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_4}+3.844746237447211\,{\it t_3}+1.644733775567609\,{\it t_2}+0.2405871291288192\,{\it t_1}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_3}+1.0\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\sqrt{\rho_\beta}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\beta^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_\beta^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{10}} &=& -6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_9}+2.225569421150687\,{\it t_8}+0.8004286349993635\,{\it t_7}+0.1897004325747558\,{\it t_6}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_8}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_9}+3.844746237447211\,{\it t_8}+1.644733775567609\,{\it t_7}+0.2405871291288192\,{\it t_6}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_8}+1.0\right)\\\\ 
C>   {\it t_{11}} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}^4\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{11}}^4}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{17}} &=& \left({\it t_{13}}\,{\it t_{16}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{13}}\,{\it t_{16}}\right)^{{{4}\over{3}}}-2.0\\\\ 
C>   {\it t_{18}} &=& -1.923661050931536\,{\it t_{14}}\,{\it t_{15}}\,{\it t_{17}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{{\it t_{11}}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\sqrt{{\it t_{11}}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{11}}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{24}} &=& {\it t_{23}}\,\sigma_{\beta\beta}+{\it t_{22}}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{0.0031\,{\it t_{24}}+1.0}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{0.06\,{\it t_{22}}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{0.06\,{\it t_{23}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\sqrt{\rho_s}}}\\\\ 
C>   {\it t_{30}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{32}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{33}} &=& {{1}\over{0.06\,{\it t_{32}}\,\sigma_{ss}+1.0}}\\\\ 
C>   f &=& {{0.5\,{\it t_{10}}\,\rho_\beta\,\left(0.06\,{\it t_{23}}\,\sigma_{\beta\beta}\,{\it t_{27}}\,\left(0.133402758\,{\it t_{23}}\,\sigma_{\beta\beta}\,{\it t_{27}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_\beta-{{0.25\,\sigma_{\beta\beta}}\over{\rho_\beta}}\right)}\over{\tau_\beta}}+{{0.5\,{\it t_5}\,\rho_\alpha\,\left(0.06\,{\it t_{22}}\,\sigma_{\alpha\alpha}\,{\it t_{26}}\,\left(0.133402758\,{\it t_{22}}\,\sigma_{\alpha\alpha}\,{\it t_{26}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_\alpha-{{0.25\,\sigma_{\alpha\alpha}}\over{\rho_\alpha}}\right)}\over{\tau_\alpha}}+1.0\,\left({\it t_{11}}\,\left(-0.0598079686361971\,{\it t_{14}}\,{\it t_{15}}\,\left(0.1274696188700087\,{\it t_{12}}+1.0\right)\,{\it t_{17}}\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{21}}+3.844746237447211\,{\it t_{12}}+1.644733775567609\,{\it t_{20}}+0.2405871291288192\,{\it t_{19}}}}+1.0\right)+0.01975163132168077\,\left(0.06901399211255825\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.923661050931536\,{\it t_{17}}\right)\,\log \left({{29.60881320326807}\over{8.15741470348764\,{\it t_{21}}+2.247591863577616\,{\it t_{12}}+0.4300972471276643\,{\it t_{20}}+0.1911512595127337\,{\it t_{19}}}}+1.0\right)-0.0621814\,\left(0.1325688999052018\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.0\right)\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{21}}+2.225569421150687\,{\it t_{12}}+0.8004286349993635\,{\it t_{20}}+0.1897004325747558\,{\it t_{19}}}}+1.0\right)\right)-{\it t_{10}}\,\rho_\beta-{\it t_5}\,\rho_\alpha\right)\,\left(0.0031\,{\it t_{24}}\,{\it t_{25}}\,\left(5.9515308-0.03459689187\,{\it t_{24}}\,{\it t_{25}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.5\,\left(-6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{31}}+2.225569421150687\,{\it t_{30}}+0.8004286349993635\,{\it t_{29}}+0.1897004325747558\,{\it t_{28}}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_{30}}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{31}}+3.844746237447211\,{\it t_{30}}+1.644733775567609\,{\it t_{29}}+0.2405871291288192\,{\it t_{28}}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_{30}}+1.0\right)\right)\,\rho_s\,\left(0.06\,{\it t_{32}}\,\sigma_{ss}\,{\it t_{33}}\,\left(0.133402758\,{\it t_{32}}\,\sigma_{ss}\,{\it t_{33}}-2.5960897\right)+1.0\right)\,\left(2.0\,\tau_s-{{0.25\,\sigma_{ss}}\over{\rho_s}}\right)}\over{\tau_s}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
      subroutine nwxcm_c_dldf_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd derivative wrt tau 
                                             !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,NCOL_AMAT3)  !< [Output] The 3rd derivative wrt rho 
      double precision Cmat3(nq,NCOL_CMAT3)  !< [Output] The 3rd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat3(nq,NCOL_MMAT3)  !< [Output] The 3rd derivative wrt tau 
                                             !< and possibly rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
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
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
