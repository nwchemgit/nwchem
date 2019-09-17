C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_camlsd.F 
C> The nwxcm_x_camlsd functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_camlsd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_4} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}}}\\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_7} &=& {\it t_3}^2\\\\ 
C>   {\it t_8} &=& \rho_\alpha^{{{2}\over{3}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{10}} &=& e^ {- 15.19266624115207\,{\it t_8}\,{
C>    \it t_9} }\\\\ 
C>   {\it t_{11}} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\beta^{{{2}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- 15.19266624115207\,{\it t_{13}}\,{
C>    \it t_9} }\\\\ 
C>   {\it t_{15}} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{15}}}}\\\\ 
C>   {\it t_{17}} &=& \rho_s^{{{2}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& e^ {- 15.19266624115207\,{\it t_{17}}\,{
C>    \it t_9} }\\\\ 
C>   f &=& -0.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{12}}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_{12}}\,\left({
C>    \it t_{14}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{14}}
C>    -1.0\right)}\over{{\it t_{13}}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_{11}}\right)\right)+{\it t_1}+1.0\right)
C>    -0.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_5}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_5}\,\left({
C>    \it t_{10}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{10}}
C>    -1.0\right)}\over{{\it t_8}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_4}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{16}}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_{16}}\,\left({
C>    \it t_{18}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{18}}
C>    -1.0\right)}\over{{\it t_{17}}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_{15}}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  
C>    , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_camlsd(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhoa**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhoa**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t1*t15*wght
            Amat(iq,D1_RA) = (-9.305257363490993d-1*t1*(1.14025105647148
     1         35d-1*t14*t16*t2*t3-3.4207531694144505d-1*t2*t3*t5*(2.598
     2         5180598138485d+0*t11*t6*t9+2.5655648770608375d-1*t3*t5*(2
     3         .1940410461360632d-2*t12*t7/rhoa**1.6666666666666669d+0-1
     4         .0128444160768044d+1*t10*t11*t5+3.333333333333333d-1*t11/
     5         rhoa)-8.551882923536125d-2*t13*t16*t3))-1.240700981798799
     6         d+0*t15*t4)*wght+Amat(iq,D1_RA)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = -param(1)
            t3 = param(2)
            t4 = param(3)
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t5
            t7 = 1/t4
            t8 = t4**2
            t9 = rhoa**6.666666666666666d-1
            t10 = 1/t9
            t11 = 1/t8
            t12 = exp(-1.5192666241152067d+1*t11*t9)
            t13 = t12-1.0d+0
            t14 = -3.291061569204095d-2*t10*t13*t8+t12-1.5d+0
            t15 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t5*t7)
     1         +2.5655648770608375d-1*t14*t4*t6
            t16 = -3.4207531694144505d-1*t15*t3*t4*t6+t2+1.0d+0
            t17 = rhob**1.3333333333333333d+0
            t18 = rhob**3.333333333333333d-1
            t19 = 1/t18
            t20 = rhob**6.666666666666666d-1
            t21 = 1/t20
            t22 = exp(-1.5192666241152067d+1*t11*t20)
            t23 = t22-1.0d+0
            t24 = -3.291061569204095d-2*t21*t23*t8+t22-1.5d+0
            t25 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t18*t7
     1         )+2.5655648770608375d-1*t19*t24*t4
            t26 = -3.4207531694144505d-1*t19*t25*t3*t4+t2+1.0d+0
            t27 = 1/t1
            t28 = 5.641895835477563d-1
            t29 = 1/t17
            fnc(iq) = (-9.305257363490993d-1*t17*t26-9.305257363490993d-
     1         1*t1*t16)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490993d-1*t1*(1.14025105647148
     1         35d-1*t15*t27*t3*t4-3.4207531694144505d-1*t3*t4*t6*(2.565
     2         5648770608375d-1*t4*t6*(2.1940410461360632d-2*t13*t8/rhoa
     3         **1.6666666666666669d+0-1.0128444160768044d+1*t11*t12*t6+
     4         3.333333333333333d-1*t12/rhoa)+4.605753341764586d+0*t10*t
     5         12*t28*t7-8.551882923536125d-2*t14*t27*t4))-1.24070098179
     6         8799d+0*t16*t5)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490993d-1*t17*(1.1402510564714
     1         835d-1*t25*t29*t3*t4-3.4207531694144505d-1*t19*t3*t4*(2.5
     2         655648770608375d-1*t19*t4*(2.1940410461360632d-2*t23*t8/r
     3         hob**1.6666666666666669d+0-1.0128444160768044d+1*t11*t19*
     4         t22+3.333333333333333d-1*t22/rhob)+4.605753341764586d+0*t
     5         21*t22*t28*t7-8.551882923536125d-2*t24*t29*t4))-1.2407009
     6         81798799d+0*t18*t26)*wght+Amat(iq,D1_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhoa**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhoa**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t1*t15*wght
            Amat(iq,D1_RA) = -9.305257363490993d-1*t1*(1.140251056471483
     1         5d-1*t14*t16*t2*t3-3.4207531694144505d-1*t2*t3*t5*(2.5985
     2         180598138485d+0*t11*t6*t9+2.5655648770608375d-1*t3*t5*(2.
     3         1940410461360632d-2*t12*t7/rhoa**1.6666666666666669d+0-1.
     4         0128444160768044d+1*t10*t11*t5+3.333333333333333d-1*t11/r
     5         hoa)-8.551882923536125d-2*t13*t16*t3))*wght-1.24070098179
     6         8799d+0*t15*t4*wght+Amat(iq,D1_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhob**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhob**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t1*t15*wght
            Amat(iq,D1_RB) = -9.305257363490993d-1*t1*(1.140251056471483
     1         5d-1*t14*t16*t2*t3-3.4207531694144505d-1*t2*t3*t5*(2.5985
     2         180598138485d+0*t11*t6*t9+2.5655648770608375d-1*t3*t5*(2.
     3         1940410461360632d-2*t12*t7/rhob**1.6666666666666669d+0-1.
     4         0128444160768044d+1*t10*t11*t5+3.333333333333333d-1*t11/r
     5         hob)-8.551882923536125d-2*t13*t16*t3))*wght-1.24070098179
     6         8799d+0*t15*t4*wght+Amat(iq,D1_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_camlsd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_4} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}}}\\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_7} &=& {\it t_3}^2\\\\ 
C>   {\it t_8} &=& \rho_\alpha^{{{2}\over{3}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{10}} &=& e^ {- 15.19266624115207\,{\it t_8}\,{
C>    \it t_9} }\\\\ 
C>   {\it t_{11}} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\beta^{{{2}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- 15.19266624115207\,{\it t_{13}}\,{
C>    \it t_9} }\\\\ 
C>   {\it t_{15}} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{15}}}}\\\\ 
C>   {\it t_{17}} &=& \rho_s^{{{2}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& e^ {- 15.19266624115207\,{\it t_{17}}\,{
C>    \it t_9} }\\\\ 
C>   f &=& -0.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{12}}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_{12}}\,\left({
C>    \it t_{14}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{14}}
C>    -1.0\right)}\over{{\it t_{13}}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_{11}}\right)\right)+{\it t_1}+1.0\right)
C>    -0.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_5}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_5}\,\left({
C>    \it t_{10}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{10}}
C>    -1.0\right)}\over{{\it t_8}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_4}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{16}}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_{16}}\,\left({
C>    \it t_{18}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{18}}
C>    -1.0\right)}\over{{\it t_{17}}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_{15}}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  
C>    , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_camlsd_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhoa**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhoa**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            t17 = 5.641895835477563d-1
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 1/rhoa
            t20 = 2.1940410461360632d-2*t12*t18*t7-1.0128444160768044d+1
     1         *t10*t11*t5+3.333333333333333d-1*t11*t19
            t21 = 4.605753341764586d+0*t11*t17*t6*t9+2.5655648770608375d
     1         -1*t20*t3*t5-8.551882923536125d-2*t13*t16*t3
            t22 = 1.1402510564714835d-1*t14*t16*t2*t3-3.4207531694144505
     1         d-1*t2*t21*t3*t5
            t23 = 1/rhoa**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t1*t15*wght
            Amat(iq,D1_RA) = (-1.240700981798799d+0*t15*t4-9.30525736349
     1         0993d-1*t1*t22)*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = (-9.305257363490993d-1*t1*(-3.420753169
     1         4144505d-1*t2*t3*t5*(2.5655648770608375d-1*t3*t5*(1.02585
     2         38111779629d+2*t11*t9/t3**4-3.6567350768934387d-2*t12*t7/
     3         rhoa**2.6666666666666666d+0-5.555555555555555d-1*t11/rhoa
     4         **2)-3.0705022278430566d+0*t11*t17*t18*t6+1.1402510564714
     5         832d-1*t13*t23*t3-1.710376584707225d-1*t16*t20*t3-4.66491
     6         1554033342d+1*t11*t17*t19/t3**3)-1.5203347419619778d-1*t1
     7         4*t2*t23*t3+2.280502112942967d-1*t16*t2*t21*t3)-4.1356699
     8         3932933d-1*t15*t9-2.481401963597598d+0*t22*t4)*wght+Amat2
     9         (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = -param(1)
            t3 = param(2)
            t4 = param(3)
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t5
            t7 = 1/t4
            t8 = t4**2
            t9 = rhoa**6.666666666666666d-1
            t10 = 1/t9
            t11 = 1/t8
            t12 = exp(-1.5192666241152067d+1*t11*t9)
            t13 = t12-1.0d+0
            t14 = -3.291061569204095d-2*t10*t13*t8+t12-1.5d+0
            t15 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t5*t7)
     1         +2.5655648770608375d-1*t14*t4*t6
            t16 = -3.4207531694144505d-1*t15*t3*t4*t6+t2+1.0d+0
            t17 = rhob**1.3333333333333333d+0
            t18 = rhob**3.333333333333333d-1
            t19 = 1/t18
            t20 = rhob**6.666666666666666d-1
            t21 = 1/t20
            t22 = exp(-1.5192666241152067d+1*t11*t20)
            t23 = t22-1.0d+0
            t24 = -3.291061569204095d-2*t21*t23*t8+t22-1.5d+0
            t25 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t18*t7
     1         )+2.5655648770608375d-1*t19*t24*t4
            t26 = -3.4207531694144505d-1*t19*t25*t3*t4+t2+1.0d+0
            t27 = 1/t1
            t28 = 5.641895835477563d-1
            t29 = 1/rhoa**1.6666666666666669d+0
            t30 = 1/rhoa
            t31 = 2.1940410461360632d-2*t13*t29*t8-1.0128444160768044d+1
     1         *t11*t12*t6+3.333333333333333d-1*t12*t30
            t32 = 4.605753341764586d+0*t10*t12*t28*t7+2.5655648770608375
     1         d-1*t31*t4*t6-8.551882923536125d-2*t14*t27*t4
            t33 = 1.1402510564714835d-1*t15*t27*t3*t4-3.4207531694144505
     1         d-1*t3*t32*t4*t6
            t34 = 1/t17
            t35 = 1/rhob**1.6666666666666669d+0
            t36 = 1/rhob
            t37 = 2.1940410461360632d-2*t23*t35*t8+3.333333333333333d-1*
     1         t22*t36-1.0128444160768044d+1*t11*t19*t22
            t38 = 4.605753341764586d+0*t21*t22*t28*t7+2.5655648770608375
     1         d-1*t19*t37*t4-8.551882923536125d-2*t24*t34*t4
            t39 = 1.1402510564714835d-1*t25*t3*t34*t4-3.4207531694144505
     1         d-1*t19*t3*t38*t4
            t40 = 1/rhoa**2.3333333333333334d+0
            t41 = 1/t4**3
            t42 = 1/t4**4
            t43 = 1/rhob**2.3333333333333334d+0
            fnc(iq) = (-9.305257363490993d-1*t17*t26-9.305257363490993d-
     1         1*t1*t16)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.240700981798799d+0*t16*t5-9.30525736349
     1         0993d-1*t1*t33)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490993d-1*t17*t39-1.2407009817
     1         98799d+0*t18*t26)*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = (-9.305257363490993d-1*t1*(-3.420753169
     1         4144505d-1*t3*t4*t6*(2.5655648770608375d-1*t4*t6*(-3.6567
     2         350768934387d-2*t13*t8/rhoa**2.6666666666666666d+0+1.0258
     3         538111779629d+2*t10*t12*t42-5.555555555555555d-1*t12/rhoa
     4         **2)-3.0705022278430566d+0*t12*t28*t29*t7-4.6649115540333
     5         42d+1*t12*t28*t30*t41+1.1402510564714832d-1*t14*t4*t40-1.
     6         710376584707225d-1*t27*t31*t4)-1.5203347419619778d-1*t15*
     7         t3*t4*t40+2.280502112942967d-1*t27*t3*t32*t4)-2.481401963
     8         597598d+0*t33*t5-4.13566993932933d-1*t10*t16)*wght+Amat2(
     9         iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363490993d-1*t17*(-3.42075316
     1         94144505d-1*t19*t3*t4*(2.5655648770608375d-1*t19*t4*(-3.6
     2         567350768934387d-2*t23*t8/rhob**2.6666666666666666d+0+1.0
     3         258538111779629d+2*t21*t22*t42-5.555555555555555d-1*t22/r
     4         hob**2)-3.0705022278430566d+0*t22*t28*t35*t7+1.1402510564
     5         714832d-1*t24*t4*t43-4.664911554033342d+1*t22*t28*t36*t41
     6         -1.710376584707225d-1*t34*t37*t4)-1.5203347419619778d-1*t
     7         25*t3*t4*t43+2.280502112942967d-1*t3*t34*t38*t4)-2.481401
     8         963597598d+0*t18*t39-4.13566993932933d-1*t21*t26)*wght+Am
     9         at2(iq,D2_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhoa**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhoa**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            t17 = 5.641895835477563d-1
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 1/rhoa
            t20 = 2.1940410461360632d-2*t12*t18*t7-1.0128444160768044d+1
     1         *t10*t11*t5+3.333333333333333d-1*t11*t19
            t21 = 4.605753341764586d+0*t11*t17*t6*t9+2.5655648770608375d
     1         -1*t20*t3*t5-8.551882923536125d-2*t13*t16*t3
            t22 = 1.1402510564714835d-1*t14*t16*t2*t3-3.4207531694144505
     1         d-1*t2*t21*t3*t5
            t23 = 1/rhoa**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t1*t15*wght
            Amat(iq,D1_RA) = -1.240700981798799d+0*t15*t4*wght-9.3052573
     1         63490993d-1*t1*t22*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = -9.305257363490993d-1*t1*(-3.4207531694
     1         144505d-1*t2*t3*t5*(2.5655648770608375d-1*t3*t5*(1.025853
     2         8111779629d+2*t11*t9/t3**4-3.6567350768934387d-2*t12*t7/r
     3         hoa**2.6666666666666666d+0-5.555555555555555d-1*t11/rhoa*
     4         *2)-3.0705022278430566d+0*t11*t17*t18*t6+1.14025105647148
     5         32d-1*t13*t23*t3-1.710376584707225d-1*t16*t20*t3-4.664911
     6         554033342d+1*t11*t17*t19/t3**3)-1.5203347419619778d-1*t14
     7         *t2*t23*t3+2.280502112942967d-1*t16*t2*t21*t3)*wght-4.135
     8         66993932933d-1*t15*t9*wght-2.481401963597598d+0*t22*t4*wg
     9         ht+Amat2(iq,D2_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhob**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhob**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            t17 = 5.641895835477563d-1
            t18 = 1/rhob**1.6666666666666669d+0
            t19 = 1/rhob
            t20 = 2.1940410461360632d-2*t12*t18*t7-1.0128444160768044d+1
     1         *t10*t11*t5+3.333333333333333d-1*t11*t19
            t21 = 4.605753341764586d+0*t11*t17*t6*t9+2.5655648770608375d
     1         -1*t20*t3*t5-8.551882923536125d-2*t13*t16*t3
            t22 = 1.1402510564714835d-1*t14*t16*t2*t3-3.4207531694144505
     1         d-1*t2*t21*t3*t5
            t23 = 1/rhob**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t1*t15*wght
            Amat(iq,D1_RB) = -1.240700981798799d+0*t15*t4*wght-9.3052573
     1         63490993d-1*t1*t22*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RB_RB) = -9.305257363490993d-1*t1*(-3.4207531694
     1         144505d-1*t2*t3*t5*(2.5655648770608375d-1*t3*t5*(1.025853
     2         8111779629d+2*t11*t9/t3**4-3.6567350768934387d-2*t12*t7/r
     3         hob**2.6666666666666666d+0-5.555555555555555d-1*t11/rhob*
     4         *2)-3.0705022278430566d+0*t11*t17*t18*t6+1.14025105647148
     5         32d-1*t13*t23*t3-1.710376584707225d-1*t16*t20*t3-4.664911
     6         554033342d+1*t11*t17*t19/t3**3)-1.5203347419619778d-1*t14
     7         *t2*t23*t3+2.280502112942967d-1*t16*t2*t21*t3)*wght-4.135
     8         66993932933d-1*t15*t9*wght-2.481401963597598d+0*t22*t4*wg
     9         ht+Amat2(iq,D2_RB_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_camlsd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_4} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}}}\\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_7} &=& {\it t_3}^2\\\\ 
C>   {\it t_8} &=& \rho_\alpha^{{{2}\over{3}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{10}} &=& e^ {- 15.19266624115207\,{\it t_8}\,{
C>    \it t_9} }\\\\ 
C>   {\it t_{11}} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\beta^{{{2}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- 15.19266624115207\,{\it t_{13}}\,{
C>    \it t_9} }\\\\ 
C>   {\it t_{15}} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{15}}}}\\\\ 
C>   {\it t_{17}} &=& \rho_s^{{{2}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& e^ {- 15.19266624115207\,{\it t_{17}}\,{
C>    \it t_9} }\\\\ 
C>   f &=& -0.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{12}}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_{12}}\,\left({
C>    \it t_{14}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{14}}
C>    -1.0\right)}\over{{\it t_{13}}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_{11}}\right)\right)+{\it t_1}+1.0\right)
C>    -0.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_5}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_5}\,\left({
C>    \it t_{10}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{10}}
C>    -1.0\right)}\over{{\it t_8}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_4}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left(
C>    -0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{16}}\,
C>    \left(0.2565564877060837\,{\it t_3}\,{\it t_{16}}\,\left({
C>    \it t_{18}}-{{0.03291061569204095\,{\it t_7}\,\left({\it t_{18}}
C>    -1.0\right)}\over{{\it t_{17}}}}-1.5\right)+1.77245385090552
C>    \,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{
C>    \it t_{15}}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  
C>    , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_camlsd_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhoa**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhoa**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            t17 = 5.641895835477563d-1
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 1/rhoa
            t20 = 2.1940410461360632d-2*t12*t18*t7-1.0128444160768044d+1
     1         *t10*t11*t5+3.333333333333333d-1*t11*t19
            t21 = 4.605753341764586d+0*t11*t17*t6*t9+2.5655648770608375d
     1         -1*t20*t3*t5-8.551882923536125d-2*t13*t16*t3
            t22 = 1.1402510564714835d-1*t14*t16*t2*t3-3.4207531694144505
     1         d-1*t2*t21*t3*t5
            t23 = 1/rhoa**2.3333333333333334d+0
            t24 = 1/t3**3
            t25 = 1/rhoa**2.6666666666666666d+0
            t26 = 1/rhoa**2
            t27 = 1/t3**4
            t28 = 1.0258538111779629d+2*t11*t27*t9-3.6567350768934387d-2
     1         *t12*t25*t7-5.555555555555555d-1*t11*t26
            t29 = -3.0705022278430566d+0*t11*t17*t18*t6+2.56556487706083
     1         75d-1*t28*t3*t5+1.1402510564714832d-1*t13*t23*t3-1.710376
     2         584707225d-1*t16*t20*t3-4.664911554033342d+1*t11*t17*t19*
     3         t24
            t30 = -3.4207531694144505d-1*t2*t29*t3*t5-1.5203347419619778
     1         d-1*t14*t2*t23*t3+2.280502112942967d-1*t16*t2*t21*t3
            t31 = 1/rhoa**3.3333333333333337d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t1*t15*wght
            Amat(iq,D1_RA) = (-1.240700981798799d+0*t15*t4-9.30525736349
     1         0993d-1*t1*t22)*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = (-4.13566993932933d-1*t15*t9-2.48140196
     1         3597598d+0*t22*t4-9.305257363490993d-1*t1*t30)*wght+Amat2
     2         (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat3(iq,D3_RA_RA_RA) = (-1.240700981798799d+0*t22*t9-9.3052
     1         57363490993d-1*t1*(-3.4207531694144505d-1*t2*t3*t5*(2.565
     2         5648770608375d-1*t3*t5*(9.751293538382502d-2*t12*t7/rhoa*
     3         *3.6666666666666664d+0-1.0390303043627082d+3*t11*t19/t3**
     4         6-6.839025407853085d+1*t11*t18*t27+5.626913422648913d+0*t
     5         10*t11*t23+1.4814814814814814d+0*t11/rhoa**3)+5.117503713
     6         071762d+0*t11*t17*t25*t6-2.660585798433461d-1*t13*t3*t31-
     7         2.5655648770608375d-1*t16*t28*t3+3.42075316941445d-1*t20*
     8         t23*t3+4.724829618994838d+2*t11*t16*t17/t3**5+7.774852590
     9         055571d+1*t11*t17*t24*t26)+3.547447731244615d-1*t14*t2*t3
     :         *t31+3.4207531694144505d-1*t16*t2*t29*t3-4.56100422588593
     ;         35d-1*t2*t21*t23*t3)-3.722102945396397d+0*t30*t4+2.757113
     <         2928862196d-1*t15*t18)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = -param(1)
            t3 = param(2)
            t4 = param(3)
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t5
            t7 = 1/t4
            t8 = t4**2
            t9 = rhoa**6.666666666666666d-1
            t10 = 1/t9
            t11 = 1/t8
            t12 = exp(-1.5192666241152067d+1*t11*t9)
            t13 = t12-1.0d+0
            t14 = -3.291061569204095d-2*t10*t13*t8+t12-1.5d+0
            t15 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t5*t7)
     1         +2.5655648770608375d-1*t14*t4*t6
            t16 = -3.4207531694144505d-1*t15*t3*t4*t6+t2+1.0d+0
            t17 = rhob**1.3333333333333333d+0
            t18 = rhob**3.333333333333333d-1
            t19 = 1/t18
            t20 = rhob**6.666666666666666d-1
            t21 = 1/t20
            t22 = exp(-1.5192666241152067d+1*t11*t20)
            t23 = t22-1.0d+0
            t24 = -3.291061569204095d-2*t21*t23*t8+t22-1.5d+0
            t25 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t18*t7
     1         )+2.5655648770608375d-1*t19*t24*t4
            t26 = -3.4207531694144505d-1*t19*t25*t3*t4+t2+1.0d+0
            t27 = 1/t1
            t28 = 5.641895835477563d-1
            t29 = 1/rhoa**1.6666666666666669d+0
            t30 = 1/rhoa
            t31 = 2.1940410461360632d-2*t13*t29*t8-1.0128444160768044d+1
     1         *t11*t12*t6+3.333333333333333d-1*t12*t30
            t32 = 4.605753341764586d+0*t10*t12*t28*t7+2.5655648770608375
     1         d-1*t31*t4*t6-8.551882923536125d-2*t14*t27*t4
            t33 = 1.1402510564714835d-1*t15*t27*t3*t4-3.4207531694144505
     1         d-1*t3*t32*t4*t6
            t34 = 1/t17
            t35 = 1/rhob**1.6666666666666669d+0
            t36 = 1/rhob
            t37 = 2.1940410461360632d-2*t23*t35*t8+3.333333333333333d-1*
     1         t22*t36-1.0128444160768044d+1*t11*t19*t22
            t38 = 4.605753341764586d+0*t21*t22*t28*t7+2.5655648770608375
     1         d-1*t19*t37*t4-8.551882923536125d-2*t24*t34*t4
            t39 = 1.1402510564714835d-1*t25*t3*t34*t4-3.4207531694144505
     1         d-1*t19*t3*t38*t4
            t40 = 1/rhoa**2.3333333333333334d+0
            t41 = 1/t4**3
            t42 = 1/rhoa**2.6666666666666666d+0
            t43 = 1/rhoa**2
            t44 = 1/t4**4
            t45 = -3.6567350768934387d-2*t13*t42*t8+1.0258538111779629d+
     1         2*t10*t12*t44-5.555555555555555d-1*t12*t43
            t46 = -3.0705022278430566d+0*t12*t28*t29*t7+2.56556487706083
     1         75d-1*t4*t45*t6-4.664911554033342d+1*t12*t28*t30*t41+1.14
     2         02510564714832d-1*t14*t4*t40-1.710376584707225d-1*t27*t31
     3         *t4
            t47 = -3.4207531694144505d-1*t3*t4*t46*t6-1.5203347419619778
     1         d-1*t15*t3*t4*t40+2.280502112942967d-1*t27*t3*t32*t4
            t48 = 1/rhob**2.3333333333333334d+0
            t49 = 1/rhob**2.6666666666666666d+0
            t50 = 1/rhob**2
            t51 = -3.6567350768934387d-2*t23*t49*t8-5.555555555555555d-1
     1         *t22*t50+1.0258538111779629d+2*t21*t22*t44
            t52 = -3.0705022278430566d+0*t22*t28*t35*t7+2.56556487706083
     1         75d-1*t19*t4*t51+1.1402510564714832d-1*t24*t4*t48-4.66491
     2         1554033342d+1*t22*t28*t36*t41-1.710376584707225d-1*t34*t3
     3         7*t4
            t53 = -3.4207531694144505d-1*t19*t3*t4*t52-1.520334741961977
     1         8d-1*t25*t3*t4*t48+2.280502112942967d-1*t3*t34*t38*t4
            t54 = 1/rhoa**3.3333333333333337d+0
            t55 = 1/t4**5
            t56 = 1/t4**6
            t57 = 1/rhob**3.3333333333333337d+0
            fnc(iq) = (-9.305257363490993d-1*t17*t26-9.305257363490993d-
     1         1*t1*t16)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.240700981798799d+0*t16*t5-9.30525736349
     1         0993d-1*t1*t33)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490993d-1*t17*t39-1.2407009817
     1         98799d+0*t18*t26)*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = (-2.481401963597598d+0*t33*t5-9.3052573
     1         63490993d-1*t1*t47-4.13566993932933d-1*t10*t16)*wght+Amat
     2         2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363490993d-1*t17*t53-2.481401
     1         963597598d+0*t18*t39-4.13566993932933d-1*t21*t26)*wght+Am
     2         at2(iq,D2_RB_RB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490993d-1*t1*(-3.420753
     1         1694144505d-1*t3*t4*t6*(2.5655648770608375d-1*t4*t6*(9.75
     2         1293538382502d-2*t13*t8/rhoa**3.6666666666666664d+0-1.039
     3         0303043627082d+3*t12*t30*t56-6.839025407853085d+1*t12*t29
     4         *t44+5.626913422648913d+0*t11*t12*t40+1.4814814814814814d
     5         +0*t12/rhoa**3)+5.117503713071762d+0*t12*t28*t42*t7+4.724
     6         829618994838d+2*t12*t27*t28*t55-2.660585798433461d-1*t14*
     7         t4*t54-2.5655648770608375d-1*t27*t4*t45+7.774852590055571
     8         d+1*t12*t28*t41*t43+3.42075316941445d-1*t31*t4*t40)+3.547
     9         447731244615d-1*t15*t3*t4*t54+3.4207531694144505d-1*t27*t
     :         3*t4*t46-4.5610042258859335d-1*t3*t32*t4*t40)-3.722102945
     ;         396397d+0*t47*t5-1.240700981798799d+0*t10*t33+2.757113292
     <         8862196d-1*t16*t29)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-9.305257363490993d-1*t17*(-3.42075
     1         31694144505d-1*t19*t3*t4*(2.5655648770608375d-1*t19*t4*(9
     2         .751293538382502d-2*t23*t8/rhob**3.6666666666666664d+0-1.
     3         0390303043627082d+3*t22*t36*t56+5.626913422648913d+0*t11*
     4         t22*t48-6.839025407853085d+1*t22*t35*t44+1.48148148148148
     5         14d+0*t22/rhob**3)+5.117503713071762d+0*t22*t28*t49*t7-2.
     6         660585798433461d-1*t24*t4*t57+4.724829618994838d+2*t22*t2
     7         8*t34*t55-2.5655648770608375d-1*t34*t4*t51+7.774852590055
     8         571d+1*t22*t28*t41*t50+3.42075316941445d-1*t37*t4*t48)+3.
     9         547447731244615d-1*t25*t3*t4*t57+3.4207531694144505d-1*t3
     :         *t34*t4*t52-4.5610042258859335d-1*t3*t38*t4*t48)-3.722102
     ;         945396397d+0*t18*t53-1.240700981798799d+0*t21*t39+2.75711
     <         32928862196d-1*t26*t35)*wght+Amat3(iq,D3_RB_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhoa**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhoa**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            t17 = 5.641895835477563d-1
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 1/rhoa
            t20 = 2.1940410461360632d-2*t12*t18*t7-1.0128444160768044d+1
     1         *t10*t11*t5+3.333333333333333d-1*t11*t19
            t21 = 4.605753341764586d+0*t11*t17*t6*t9+2.5655648770608375d
     1         -1*t20*t3*t5-8.551882923536125d-2*t13*t16*t3
            t22 = 1.1402510564714835d-1*t14*t16*t2*t3-3.4207531694144505
     1         d-1*t2*t21*t3*t5
            t23 = 1/rhoa**2.3333333333333334d+0
            t24 = 1/t3**3
            t25 = 1/rhoa**2.6666666666666666d+0
            t26 = 1/rhoa**2
            t27 = 1/t3**4
            t28 = 1.0258538111779629d+2*t11*t27*t9-3.6567350768934387d-2
     1         *t12*t25*t7-5.555555555555555d-1*t11*t26
            t29 = -3.0705022278430566d+0*t11*t17*t18*t6+2.56556487706083
     1         75d-1*t28*t3*t5+1.1402510564714832d-1*t13*t23*t3-1.710376
     2         584707225d-1*t16*t20*t3-4.664911554033342d+1*t11*t17*t19*
     3         t24
            t30 = -3.4207531694144505d-1*t2*t29*t3*t5-1.5203347419619778
     1         d-1*t14*t2*t23*t3+2.280502112942967d-1*t16*t2*t21*t3
            t31 = 1/rhoa**3.3333333333333337d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t1*t15*wght
            Amat(iq,D1_RA) = -1.240700981798799d+0*t15*t4*wght-9.3052573
     1         63490993d-1*t1*t22*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = -4.13566993932933d-1*t15*t9*wght-2.4814
     1         01963597598d+0*t22*t4*wght-9.305257363490993d-1*t1*t30*wg
     2         ht+Amat2(iq,D2_RA_RA)
            Amat3(iq,D3_RA_RA_RA) = -1.240700981798799d+0*t22*t9*wght-9.
     1         305257363490993d-1*t1*(-3.4207531694144505d-1*t2*t3*t5*(2
     2         .5655648770608375d-1*t3*t5*(9.751293538382502d-2*t12*t7/r
     3         hoa**3.6666666666666664d+0-1.0390303043627082d+3*t11*t19/
     4         t3**6-6.839025407853085d+1*t11*t18*t27+5.626913422648913d
     5         +0*t10*t11*t23+1.4814814814814814d+0*t11/rhoa**3)+5.11750
     6         3713071762d+0*t11*t17*t25*t6-2.660585798433461d-1*t13*t3*
     7         t31-2.5655648770608375d-1*t16*t28*t3+3.42075316941445d-1*
     8         t20*t23*t3+4.724829618994838d+2*t11*t16*t17/t3**5+7.77485
     9         2590055571d+1*t11*t17*t24*t26)+3.547447731244615d-1*t14*t
     :         2*t3*t31+3.4207531694144505d-1*t16*t2*t29*t3-4.5610042258
     ;         859335d-1*t2*t21*t23*t3)*wght-3.722102945396397d+0*t30*t4
     <         *wght+2.7571132928862196d-1*t15*t18*wght+Amat3(iq,D3_RA_R
     =         A_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(2)
            t3 = param(3)
            t4 = rhob**3.333333333333333d-1
            t5 = 1/t4
            t6 = 1/t3
            t7 = t3**2
            t8 = rhob**6.666666666666666d-1
            t9 = 1/t8
            t10 = 1/t7
            t11 = exp(-1.5192666241152067d+1*t10*t8)
            t12 = t11-1.0d+0
            t13 = -3.291061569204095d-2*t12*t7*t9+t11-1.5d+0
            t14 = 1.7724538509055202d+0*erf(3.8977770897207636d+0*t4*t6)
     1         +2.5655648770608375d-1*t13*t3*t5
            t15 = -3.4207531694144505d-1*t14*t2*t3*t5-param(1)+1.0d+0
            t16 = 1/t1
            t17 = 5.641895835477563d-1
            t18 = 1/rhob**1.6666666666666669d+0
            t19 = 1/rhob
            t20 = 2.1940410461360632d-2*t12*t18*t7-1.0128444160768044d+1
     1         *t10*t11*t5+3.333333333333333d-1*t11*t19
            t21 = 4.605753341764586d+0*t11*t17*t6*t9+2.5655648770608375d
     1         -1*t20*t3*t5-8.551882923536125d-2*t13*t16*t3
            t22 = 1.1402510564714835d-1*t14*t16*t2*t3-3.4207531694144505
     1         d-1*t2*t21*t3*t5
            t23 = 1/rhob**2.3333333333333334d+0
            t24 = 1/t3**3
            t25 = 1/rhob**2.6666666666666666d+0
            t26 = 1/rhob**2
            t27 = 1/t3**4
            t28 = 1.0258538111779629d+2*t11*t27*t9-3.6567350768934387d-2
     1         *t12*t25*t7-5.555555555555555d-1*t11*t26
            t29 = -3.0705022278430566d+0*t11*t17*t18*t6+2.56556487706083
     1         75d-1*t28*t3*t5+1.1402510564714832d-1*t13*t23*t3-1.710376
     2         584707225d-1*t16*t20*t3-4.664911554033342d+1*t11*t17*t19*
     3         t24
            t30 = -3.4207531694144505d-1*t2*t29*t3*t5-1.5203347419619778
     1         d-1*t14*t2*t23*t3+2.280502112942967d-1*t16*t2*t21*t3
            t31 = 1/rhob**3.3333333333333337d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t1*t15*wght
            Amat(iq,D1_RB) = -1.240700981798799d+0*t15*t4*wght-9.3052573
     1         63490993d-1*t1*t22*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RB_RB) = -4.13566993932933d-1*t15*t9*wght-2.4814
     1         01963597598d+0*t22*t4*wght-9.305257363490993d-1*t1*t30*wg
     2         ht+Amat2(iq,D2_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = -1.240700981798799d+0*t22*t9*wght-9.
     1         305257363490993d-1*t1*(-3.4207531694144505d-1*t2*t3*t5*(2
     2         .5655648770608375d-1*t3*t5*(9.751293538382502d-2*t12*t7/r
     3         hob**3.6666666666666664d+0-1.0390303043627082d+3*t11*t19/
     4         t3**6-6.839025407853085d+1*t11*t18*t27+5.626913422648913d
     5         +0*t10*t11*t23+1.4814814814814814d+0*t11/rhob**3)+5.11750
     6         3713071762d+0*t11*t17*t25*t6-2.660585798433461d-1*t13*t3*
     7         t31-2.5655648770608375d-1*t16*t28*t3+3.42075316941445d-1*
     8         t20*t23*t3+4.724829618994838d+2*t11*t16*t17/t3**5+7.77485
     9         2590055571d+1*t11*t17*t24*t26)+3.547447731244615d-1*t14*t
     :         2*t3*t31+3.4207531694144505d-1*t16*t2*t29*t3-4.5610042258
     ;         859335d-1*t2*t21*t23*t3)*wght-3.722102945396397d+0*t30*t4
     <         *wght+2.7571132928862196d-1*t15*t18*wght+Amat3(iq,D3_RB_R
     =         B_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
