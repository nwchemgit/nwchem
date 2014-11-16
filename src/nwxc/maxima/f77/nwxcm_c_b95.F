C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_b95.F 
C> The nwxcm_c_b95 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_b95 functional [1] 
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
C>   {\it t_{24}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{\sqrt{\rho_s}}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   f &=& {{0.1097020523068037\,{\it t_{10}}\,{\it t_6}\,\left(2.0\,\tau_\beta-{{0.25\,\sigma_{\beta\beta}}\over{\rho_\beta}}\right)}\over{\left({\it t_{24}}\,{\it t_{23}}\,\sigma_{\beta\beta}+1.0\right)^2}}+{{0.1097020523068037\,{\it t_5}\,{\it t_1}\,\left(2.0\,\tau_\alpha-{{0.25\,\sigma_{\alpha\alpha}}\over{\rho_\alpha}}\right)}\over{\left({\it t_{24}}\,{\it t_{22}}\,\sigma_{\alpha\alpha}+1.0\right)^2}}+{{1.0\,\left({\it t_{11}}\,\left(-0.0598079686361971\,{\it t_{14}}\,{\it t_{15}}\,\left(0.1274696188700087\,{\it t_{12}}+1.0\right)\,{\it t_{17}}\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{21}}+3.844746237447211\,{\it t_{12}}+1.644733775567609\,{\it t_{20}}+0.2405871291288192\,{\it t_{19}}}}+1.0\right)+0.01975163132168077\,\left(0.06901399211255825\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.923661050931536\,{\it t_{17}}\right)\,\log \left({{29.60881320326807}\over{8.15741470348764\,{\it t_{21}}+2.247591863577616\,{\it t_{12}}+0.4300972471276643\,{\it t_{20}}+0.1911512595127337\,{\it t_{19}}}}+1.0\right)-0.0621814\,\left(0.1325688999052018\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.0\right)\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{21}}+2.225569421150687\,{\it t_{12}}+0.8004286349993635\,{\it t_{20}}+0.1897004325747558\,{\it t_{19}}}}+1.0\right)\right)-{\it t_{10}}\,\rho_\beta-{\it t_5}\,\rho_\alpha\right)}\over{{\it param}\left(1\right)\,\left({\it t_{23}}\,\sigma_{\beta\beta}+{\it t_{22}}\,\sigma_{\alpha\alpha}\right)+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.1097020523068037\,\left(-6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{28}}+2.225569421150687\,{\it t_{27}}+0.8004286349993635\,{\it t_{26}}+0.1897004325747558\,{\it t_{25}}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_{27}}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{28}}+3.844746237447211\,{\it t_{27}}+1.644733775567609\,{\it t_{26}}+0.2405871291288192\,{\it t_{25}}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_{27}}+1.0\right)\right)\,{\it t_{25}}\,\left(2.0\,\tau_s-{{0.25\,\sigma_{ss}}\over{\rho_s}}\right)}\over{\left({{{\it t_{24}}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0\right)^2}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 104, 1040 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.4708298 "> 
C> 10.1063/1.4708298 </a> 
C> 
      subroutine nwxcm_c_b95(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 2*gammaaa*t1*t2+1.0d+0
            t4 = 1/t3
            t5 = 1/rhoa**6.666666666666666d-1
            t6 = rhoa**5.0d-1
            t7 = 1/t6
            t8 = 1/rhoa**3.333333333333333d-1
            t9 = 1/rhoa**1.6666666666666666d-1
            t10 = 1.1120374863094686d+1*t9+3.8447462374472113d+0*t8+1.64
     1         4733775567609d+0*t7+2.4058712912881916d-1*t5
            t11 = 3.216395899738507d+1/t10+1.0d+0
            t12 = log(t11)
            t13 = 1.274696188700087d-1*t8+1.0d+0
            t14 = 5.982550435771079d+0*t9+2.2255694211506874d+0*t8+8.004
     1         286349993635d-1*t7+1.8970043257475583d-1*t5
            t15 = 1.6081979498692536d+1/t14+1.0d+0
            t16 = log(t15)
            t17 = 1.3256889990520176d-1*t8+1.0d+0
            t18 = -6.903522198342671d-18*t16*t17-3.1090699999999993d-2*t
     1         12*t13
            t19 = 5.329846514438383d+0*t9+1.766435620181429d+0*t8+5.6598
     1         85156639418d-1*t7+1.19503784087549d-1*t5
            t20 = 1.6081979498692536d+1/t19+1.0d+0
            t21 = log(t20)
            t22 = 1.0522000558389212d-1*t8+1.0d+0
            t23 = -1.243628d-1*rhoa*t21*t22-2*rhoa*t18
            t24 = param(2)
            t25 = gammaaa*t2*t24+1.0d+0
            t26 = 1/t25**2
            t27 = 2.0d+0*taua-2.5d-1*gammaaa/rhoa
            t28 = 1/rhoa**1.6666666666666669d+0
            t29 = 1/t6**3
            t30 = 1/rhoa**1.3333333333333333d+0
            t31 = 1/rhoa**1.1666666666666667d+0
            t32 = 1.1102230246251565d-16*t17*(-9.970917392951799d-1*t31-
     1         7.418564737168958d-1*t30-4.002143174996817d-1*t29-1.26466
     2         9550498372d-1*t28)/(t14**2*t15)+1.0d+0*t13*(-1.8533958105
     3         157808d+0*t31-1.2815820791490704d+0*t30-8.223668877838045
     4         d-1*t29-1.6039141941921276d-1*t28)/(t10**2*t11)+3.0506411
     5         443514264d-19*t16*t30+1.3210398931339265d-3*t12*t30
            t33 = 1/t3**2
            t34 = 1/t25**3
            t35 = (-4.38808209227215d-1*t18*t24*t27*t34/rhoa**3.33333333
     1         33333337d+0-2.0d+0*t1*t2*t23*t33-5.485102615340187d-2*t18
     2         *t26*t28)*wght
            fnc(iq) = (2.194041046136075d-1*t18*t26*t27*t5+1.0d+0*t23*t4
     1         )*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t4*(4.361818170142819d-3*t21*t8-2*r
     1         hoa*t32+2.0d+0*rhoa*t22*(-8.883077524063971d-1*t31-5.8881
     2         18733938097d-1*t30-2.829942578319709d-1*t29-7.96691893916
     3         9932d-2*t28)/(t19**2*t20)-1.243628d-1*t21*t22-2*t18)+2.19
     4         4041046136075d-1*t26*t27*t32*t5+1.1701552246059066d+0*gam
     5         maaa*t18*t24*t27*t34/rhoa**4.333333333333333d+0+5.3333333
     6         33333333d+0*gammaaa*t1*t23*t33/rhoa**3.6666666666666664d+
     7         0-1.4626940307573832d-1*t18*t26*t27*t28+5.485102615340187
     8         d-2*gammaaa*t18*t2*t26)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t35+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t35+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = 4.38808209227215d-1*t18*t26*t5*wght+Amat(iq
     1         ,D1_TA)
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
            t1 = param(1)
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = t1*(gammabb*t3+gammaaa*t2)+1.0d+0
            t5 = 1/t4
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
            t20 = 1/rhob**6.666666666666666d-1
            t21 = rhob**5.0d-1
            t22 = 1/t21
            t23 = 1/rhob**3.333333333333333d-1
            t24 = 1/rhob**1.6666666666666666d-1
            t25 = 1.1120374863094686d+1*t24+3.8447462374472113d+0*t23+1.
     1         644733775567609d+0*t22+2.4058712912881916d-1*t20
            t26 = 3.216395899738507d+1/t25+1.0d+0
            t27 = log(t26)
            t28 = 1.274696188700087d-1*t23+1.0d+0
            t29 = 5.982550435771079d+0*t24+2.2255694211506874d+0*t23+8.0
     1         04286349993635d-1*t22+1.8970043257475583d-1*t20
            t30 = 1.6081979498692536d+1/t29+1.0d+0
            t31 = log(t30)
            t32 = 1.3256889990520176d-1*t23+1.0d+0
            t33 = -6.903522198342671d-18*t31*t32-3.1090699999999993d-2*t
     1         27*t28
            t34 = rhob+rhoa
            t35 = 1/t34**3.333333333333333d-1
            t36 = 1.3256889990520176d-1*t35+1.0d+0
            t37 = rhoa-rhob
            t38 = t37**4
            t39 = 1/t34**4
            t40 = 1/t34
            t41 = 1.0d+0-t37*t40
            t42 = t37*t40+1.0d+0
            t43 = t42**1.3333333333333333d+0+t41**1.3333333333333333d+0-
     1         2.0d+0
            t44 = -1.923661050931536d+0*t38*t39*t43
            t45 = t44+1.0d+0
            t46 = 1/t34**6.666666666666666d-1
            t47 = t34**5.0d-1
            t48 = 1/t47
            t49 = 1/t34**1.6666666666666666d-1
            t50 = 5.982550435771079d+0*t49+8.004286349993635d-1*t48+1.89
     1         70043257475583d-1*t46+2.2255694211506874d+0*t35
            t51 = 1.6081979498692536d+1/t50+1.0d+0
            t52 = log(t51)
            t53 = -6.21814d-2*t36*t45*t52
            t54 = 6.901399211255825d-2*t35+1.0d+0
            t55 = t44+1.923661050931536d+0*t43
            t56 = 8.157414703487641d+0*t49+4.3009724712766434d-1*t48+1.9
     1         115125951273373d-1*t46+2.247591863577616d+0*t35
            t57 = 2.9608813203268075d+1/t56+1.0d+0
            t58 = log(t57)
            t59 = 1.9751631321680768d-2*t54*t55*t58
            t60 = 1.274696188700087d-1*t35+1.0d+0
            t61 = 1.1120374863094686d+1*t49+1.644733775567609d+0*t48+2.4
     1         058712912881916d-1*t46+3.8447462374472113d+0*t35
            t62 = 3.216395899738507d+1/t61+1.0d+0
            t63 = log(t62)
            t64 = -5.980796863619711d-2*t38*t39*t43*t60*t63
            t65 = t34*(t64+t59+t53)-rhob*t33-rhoa*t19
            t66 = param(2)
            t67 = gammaaa*t2*t66+1.0d+0
            t68 = 1/t67**2
            t69 = 2.0d+0*taua-2.5d-1*gammaaa/rhoa
            t70 = gammabb*t3*t66+1.0d+0
            t71 = 1/t70**2
            t72 = 2.0d+0*taub-2.5d-1*gammabb/rhob
            t73 = 1/rhoa**1.6666666666666669d+0
            t74 = 1/t7**3
            t75 = 1/rhoa**1.3333333333333333d+0
            t76 = 1/rhoa**1.1666666666666667d+0
            t77 = 1.1102230246251565d-16*t18*(-9.970917392951799d-1*t76-
     1         7.418564737168958d-1*t75-4.002143174996817d-1*t74-1.26466
     2         9550498372d-1*t73)/(t15**2*t16)+1.0d+0*t14*(-1.8533958105
     3         157808d+0*t76-1.2815820791490704d+0*t75-8.223668877838045
     4         d-1*t74-1.6039141941921276d-1*t73)/(t11**2*t12)+3.0506411
     5         443514264d-19*t17*t75+1.3210398931339265d-3*t13*t75
            t78 = 1/t34**1.6666666666666669d+0
            t79 = 1/t47**3
            t80 = 1/t34**1.3333333333333333d+0
            t81 = 1/t34**1.1666666666666667d+0
            t82 = 1.0d+0*t36*t45*(-9.970917392951799d-1*t81-7.4185647371
     1         68958d-1*t80-4.002143174996817d-1*t79-1.264669550498372d-
     2         1*t78)/(t50**2*t51)
            t83 = -5.848223622634647d-1*t54*t55*(-1.3595691172479402d+0*
     1         t81-7.491972878592053d-1*t80-2.1504862356383217d-1*t79-1.
     2         2743417300848914d-1*t78)/(t56**2*t57)
            t84 = 1.923661050931536d+0*t38*t39*t43*t60*(-1.8533958105157
     1         808d+0*t81-1.2815820791490704d+0*t80-8.223668877838045d-1
     2         *t79-1.6039141941921276d-1*t78)/(t61**2*t62)
            t85 = 1/t34**2
            t86 = t37*t85
            t87 = -t40
            t88 = t41**3.333333333333333d-1
            t89 = -t37*t85
            t90 = t42**3.333333333333333d-1
            t91 = 1.3333333333333333d+0*(t89+t40)*t90+1.3333333333333333
     1         d+0*(t87+t86)*t88
            t92 = -1.923661050931536d+0*t38*t39*t91
            t93 = 1/t34**5
            t94 = 7.694644203726144d+0*t38*t43*t93
            t95 = t37**3
            t96 = -7.694644203726144d+0*t39*t43*t95
            t97 = 2.747773264188437d-3*t45*t52*t80
            t98 = -4.5437964274821163d-4*t55*t58*t80
            t99 = 2.541232989148493d-3*t38*t43*t63/t34**5.33333333333333
     1         3d+0
            t100 = 2.3923187454478842d-1*t38*t43*t60*t63*t93
            t101 = 1/t4**2
            t102 = 1/t67**3
            t103 = 1/rhob**1.6666666666666669d+0
            t104 = 1/t21**3
            t105 = 1/rhob**1.3333333333333333d+0
            t106 = 1/rhob**1.1666666666666667d+0
            t107 = 1.1102230246251565d-16*(-9.970917392951799d-1*t106-7.
     1         418564737168958d-1*t105-4.002143174996817d-1*t104-1.26466
     2         9550498372d-1*t103)*t32/(t29**2*t30)+3.0506411443514264d-
     3         19*t105*t31+1.0d+0*(-1.8533958105157808d+0*t106-1.2815820
     4         791490704d+0*t105-8.223668877838045d-1*t104-1.60391419419
     5         21276d-1*t103)*t28/(t25**2*t26)+1.3210398931339265d-3*t10
     6         5*t27
            t108 = 1.3333333333333333d+0*(t89+t87)*t90+1.333333333333333
     1         3d+0*(t86+t40)*t88
            t109 = -1.923661050931536d+0*t108*t38*t39
            t110 = 7.694644203726144d+0*t39*t43*t95
            t111 = 1/t70**3
            fnc(iq) = (1.0970205230680374d-1*t20*t33*t71*t72+1.097020523
     1         0680374d-1*t19*t6*t68*t69+1.0d+0*t5*t65)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t5*(t34*(t99+t98+t97+1.975163132168
     1         0768d-2*t54*t58*(t96+t94+t92+1.923661050931536d+0*t91)-6.
     2         21814d-2*t36*t52*(t96+t94+t92)-2.3923187454478842d-1*t39*
     3         t43*t60*t63*t95-5.980796863619711d-2*t38*t39*t60*t63*t91+
     4         t84+t83+t82+t100)-rhoa*t77+t64+t59+t53+6.903522198342671d
     5         -18*t17*t18+3.1090699999999993d-2*t13*t14)+1.097020523068
     6         0374d-1*t6*t68*t69*t77-7.313470153786916d-2*t19*t68*t69*t
     7         73+5.850776123029533d-1*gammaaa*t102*t19*t66*t69/rhoa**4.
     8         333333333333333d+0+2.7425513076700936d-2*gammaaa*t19*t2*t
     9         68+2.6666666666666666d+0*gammaaa*t1*t101*t65/rhoa**3.6666
     :         666666666664d+0)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t5*(t34*(t99+t98+t97+2.392318745447
     1         8842d-1*t39*t43*t60*t63*t95+1.9751631321680768d-2*t54*t58
     2         *(t94+t110+t109+1.923661050931536d+0*t108)-6.21814d-2*t36
     3         *t52*(t94+t110+t109)+t84+t83+t82-5.980796863619711d-2*t10
     4         8*t38*t39*t60*t63+t100)+t64+t59+t53+6.903522198342671d-18
     5         *t31*t32+3.1090699999999993d-2*t27*t28-rhob*t107)-7.31347
     6         0153786916d-2*t103*t33*t71*t72+1.0970205230680374d-1*t107
     7         *t20*t71*t72+5.850776123029533d-1*gammabb*t111*t33*t66*t7
     8         2/rhob**4.333333333333333d+0+2.7425513076700936d-2*gammab
     9         b*t3*t33*t71+2.6666666666666666d+0*gammabb*t1*t101*t65/rh
     :         ob**3.6666666666666664d+0)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (-2.7425513076700936d-2*t19*t68*t73-2.1940
     1         41046136075d-1*t102*t19*t66*t69/rhoa**3.3333333333333337d
     2         +0-1.0d+0*t1*t101*t2*t65)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (-2.194041046136075d-1*t111*t33*t66*t72/rh
     1         ob**3.3333333333333337d+0-2.7425513076700936d-2*t103*t33*
     2         t71-1.0d+0*t1*t101*t3*t65)*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TA) = 2.194041046136075d-1*t19*t6*t68*wght+Mmat(i
     1         q,D1_TA)
            Mmat(iq,D1_TB) = 2.194041046136075d-1*t20*t33*t71*wght+Mmat(
     1         iq,D1_TB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(2)
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = gammaaa*t1*t2+1.0d+0
            t4 = 1/t3**2
            t5 = 1/rhoa**6.666666666666666d-1
            t6 = rhoa**5.0d-1
            t7 = 1/t6
            t8 = 1/rhoa**3.333333333333333d-1
            t9 = 1/rhoa**1.6666666666666666d-1
            t10 = 1.1120374863094686d+1*t9+3.8447462374472113d+0*t8+1.64
     1         4733775567609d+0*t7+2.4058712912881916d-1*t5
            t11 = 3.216395899738507d+1/t10+1.0d+0
            t12 = log(t11)
            t13 = 1.274696188700087d-1*t8+1.0d+0
            t14 = 5.982550435771079d+0*t9+2.2255694211506874d+0*t8+8.004
     1         286349993635d-1*t7+1.8970043257475583d-1*t5
            t15 = 1.6081979498692536d+1/t14+1.0d+0
            t16 = log(t15)
            t17 = 1.3256889990520176d-1*t8+1.0d+0
            t18 = -6.903522198342671d-18*t16*t17-3.1090699999999993d-2*t
     1         12*t13
            t19 = 2.0d+0*taua-2.5d-1*gammaaa/rhoa
            t20 = 1/t3**3
            t21 = 1/rhoa**1.6666666666666669d+0
            t22 = 1/t6**3
            t23 = 1/rhoa**1.3333333333333333d+0
            t24 = 1/rhoa**1.1666666666666667d+0
            fnc(iq) = 1.0970205230680374d-1*t18*t19*t4*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0970205230680374d-1*t19*(1.11022302462515
     1         65d-16*t17*(-9.970917392951799d-1*t24-7.418564737168958d-
     2         1*t23-4.002143174996817d-1*t22-1.264669550498372d-1*t21)/
     3         (t14**2*t15)+1.0d+0*t13*(-1.8533958105157808d+0*t24-1.281
     4         5820791490704d+0*t23-8.223668877838045d-1*t22-1.603914194
     5         1921276d-1*t21)/(t10**2*t11)+3.0506411443514264d-19*t16*t
     6         23+1.3210398931339265d-3*t12*t23)*t4*t5*wght-7.3134701537
     7         86916d-2*t18*t19*t21*t4*wght+2.7425513076700936d-2*gammaa
     8         a*t18*t2*t4*wght+5.850776123029533d-1*gammaaa*t1*t18*t19*
     9         t20*wght/rhoa**4.333333333333333d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -2.7425513076700936d-2*t18*t21*t4*wght-2.1
     1         94041046136075d-1*t1*t18*t19*t20*wght/rhoa**3.33333333333
     2         33337d+0+Cmat(iq,D1_GAA)
            Mmat(iq,D1_TA) = 2.194041046136075d-1*t18*t4*t5*wght+Mmat(iq
     1         ,D1_TA)
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = gammabb*t1*t2+1.0d+0
            t4 = 1/t3**2
            t5 = 1/rhob**6.666666666666666d-1
            t6 = rhob**5.0d-1
            t7 = 1/t6
            t8 = 1/rhob**3.333333333333333d-1
            t9 = 1/rhob**1.6666666666666666d-1
            t10 = 1.1120374863094686d+1*t9+3.8447462374472113d+0*t8+1.64
     1         4733775567609d+0*t7+2.4058712912881916d-1*t5
            t11 = 3.216395899738507d+1/t10+1.0d+0
            t12 = log(t11)
            t13 = 1.274696188700087d-1*t8+1.0d+0
            t14 = 5.982550435771079d+0*t9+2.2255694211506874d+0*t8+8.004
     1         286349993635d-1*t7+1.8970043257475583d-1*t5
            t15 = 1.6081979498692536d+1/t14+1.0d+0
            t16 = log(t15)
            t17 = 1.3256889990520176d-1*t8+1.0d+0
            t18 = -6.903522198342671d-18*t16*t17-3.1090699999999993d-2*t
     1         12*t13
            t19 = 2.0d+0*taub-2.5d-1*gammabb/rhob
            t20 = 1/t3**3
            t21 = 1/rhob**1.6666666666666669d+0
            t22 = 1/t6**3
            t23 = 1/rhob**1.3333333333333333d+0
            t24 = 1/rhob**1.1666666666666667d+0
            fnc(iq) = 1.0970205230680374d-1*t18*t19*t4*t5*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0970205230680374d-1*t19*(1.11022302462515
     1         65d-16*t17*(-9.970917392951799d-1*t24-7.418564737168958d-
     2         1*t23-4.002143174996817d-1*t22-1.264669550498372d-1*t21)/
     3         (t14**2*t15)+1.0d+0*t13*(-1.8533958105157808d+0*t24-1.281
     4         5820791490704d+0*t23-8.223668877838045d-1*t22-1.603914194
     5         1921276d-1*t21)/(t10**2*t11)+3.0506411443514264d-19*t16*t
     6         23+1.3210398931339265d-3*t12*t23)*t4*t5*wght-7.3134701537
     7         86916d-2*t18*t19*t21*t4*wght+2.7425513076700936d-2*gammab
     8         b*t18*t2*t4*wght+5.850776123029533d-1*gammabb*t1*t18*t19*
     9         t20*wght/rhob**4.333333333333333d+0+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -2.7425513076700936d-2*t18*t21*t4*wght-2.1
     1         94041046136075d-1*t1*t18*t19*t20*wght/rhob**3.33333333333
     2         33337d+0+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TB) = 2.194041046136075d-1*t18*t4*t5*wght+Mmat(iq
     1         ,D1_TB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_b95 functional [1] 
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
C>   {\it t_{24}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{\sqrt{\rho_s}}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   f &=& {{0.1097020523068037\,{\it t_{10}}\,{\it t_6}\,\left(2.0\,\tau_\beta-{{0.25\,\sigma_{\beta\beta}}\over{\rho_\beta}}\right)}\over{\left({\it t_{24}}\,{\it t_{23}}\,\sigma_{\beta\beta}+1.0\right)^2}}+{{0.1097020523068037\,{\it t_5}\,{\it t_1}\,\left(2.0\,\tau_\alpha-{{0.25\,\sigma_{\alpha\alpha}}\over{\rho_\alpha}}\right)}\over{\left({\it t_{24}}\,{\it t_{22}}\,\sigma_{\alpha\alpha}+1.0\right)^2}}+{{1.0\,\left({\it t_{11}}\,\left(-0.0598079686361971\,{\it t_{14}}\,{\it t_{15}}\,\left(0.1274696188700087\,{\it t_{12}}+1.0\right)\,{\it t_{17}}\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{21}}+3.844746237447211\,{\it t_{12}}+1.644733775567609\,{\it t_{20}}+0.2405871291288192\,{\it t_{19}}}}+1.0\right)+0.01975163132168077\,\left(0.06901399211255825\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.923661050931536\,{\it t_{17}}\right)\,\log \left({{29.60881320326807}\over{8.15741470348764\,{\it t_{21}}+2.247591863577616\,{\it t_{12}}+0.4300972471276643\,{\it t_{20}}+0.1911512595127337\,{\it t_{19}}}}+1.0\right)-0.0621814\,\left(0.1325688999052018\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.0\right)\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{21}}+2.225569421150687\,{\it t_{12}}+0.8004286349993635\,{\it t_{20}}+0.1897004325747558\,{\it t_{19}}}}+1.0\right)\right)-{\it t_{10}}\,\rho_\beta-{\it t_5}\,\rho_\alpha\right)}\over{{\it param}\left(1\right)\,\left({\it t_{23}}\,\sigma_{\beta\beta}+{\it t_{22}}\,\sigma_{\alpha\alpha}\right)+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.1097020523068037\,\left(-6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{28}}+2.225569421150687\,{\it t_{27}}+0.8004286349993635\,{\it t_{26}}+0.1897004325747558\,{\it t_{25}}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_{27}}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{28}}+3.844746237447211\,{\it t_{27}}+1.644733775567609\,{\it t_{26}}+0.2405871291288192\,{\it t_{25}}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_{27}}+1.0\right)\right)\,{\it t_{25}}\,\left(2.0\,\tau_s-{{0.25\,\sigma_{ss}}\over{\rho_s}}\right)}\over{\left({{{\it t_{24}}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0\right)^2}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 104, 1040 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.4708298 "> 
C> 10.1063/1.4708298 </a> 
C> 
      subroutine nwxcm_c_b95_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_c_b95 functional [1] 
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
C>   {\it t_{24}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{\sqrt{\rho_s}}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   f &=& {{0.1097020523068037\,{\it t_{10}}\,{\it t_6}\,\left(2.0\,\tau_\beta-{{0.25\,\sigma_{\beta\beta}}\over{\rho_\beta}}\right)}\over{\left({\it t_{24}}\,{\it t_{23}}\,\sigma_{\beta\beta}+1.0\right)^2}}+{{0.1097020523068037\,{\it t_5}\,{\it t_1}\,\left(2.0\,\tau_\alpha-{{0.25\,\sigma_{\alpha\alpha}}\over{\rho_\alpha}}\right)}\over{\left({\it t_{24}}\,{\it t_{22}}\,\sigma_{\alpha\alpha}+1.0\right)^2}}+{{1.0\,\left({\it t_{11}}\,\left(-0.0598079686361971\,{\it t_{14}}\,{\it t_{15}}\,\left(0.1274696188700087\,{\it t_{12}}+1.0\right)\,{\it t_{17}}\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{21}}+3.844746237447211\,{\it t_{12}}+1.644733775567609\,{\it t_{20}}+0.2405871291288192\,{\it t_{19}}}}+1.0\right)+0.01975163132168077\,\left(0.06901399211255825\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.923661050931536\,{\it t_{17}}\right)\,\log \left({{29.60881320326807}\over{8.15741470348764\,{\it t_{21}}+2.247591863577616\,{\it t_{12}}+0.4300972471276643\,{\it t_{20}}+0.1911512595127337\,{\it t_{19}}}}+1.0\right)-0.0621814\,\left(0.1325688999052018\,{\it t_{12}}+1.0\right)\,\left({\it t_{18}}+1.0\right)\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{21}}+2.225569421150687\,{\it t_{12}}+0.8004286349993635\,{\it t_{20}}+0.1897004325747558\,{\it t_{19}}}}+1.0\right)\right)-{\it t_{10}}\,\rho_\beta-{\it t_5}\,\rho_\alpha\right)}\over{{\it param}\left(1\right)\,\left({\it t_{23}}\,\sigma_{\beta\beta}+{\it t_{22}}\,\sigma_{\alpha\alpha}\right)+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.1097020523068037\,\left(-6.903522198342671 \times 10^{-18}\,\log \left({{16.08197949869254}\over{5.98255043577108\,{\it t_{28}}+2.225569421150687\,{\it t_{27}}+0.8004286349993635\,{\it t_{26}}+0.1897004325747558\,{\it t_{25}}}}+1.0\right)\,\left(0.1325688999052018\,{\it t_{27}}+1.0\right)-0.0310907\,\log \left({{32.16395899738507}\over{11.12037486309469\,{\it t_{28}}+3.844746237447211\,{\it t_{27}}+1.644733775567609\,{\it t_{26}}+0.2405871291288192\,{\it t_{25}}}}+1.0\right)\,\left(0.1274696188700087\,{\it t_{27}}+1.0\right)\right)\,{\it t_{25}}\,\left(2.0\,\tau_s-{{0.25\,\sigma_{ss}}\over{\rho_s}}\right)}\over{\left({{{\it t_{24}}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0\right)^2}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 104, 1040 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.4708298 "> 
C> 10.1063/1.4708298 </a> 
C> 
      subroutine nwxcm_c_b95_d3(param,tol_rho,ipol,nq,wght, 
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
