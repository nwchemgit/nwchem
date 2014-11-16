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
C>   {\it t_{10}} &=& e^ {- 15.19266624115207\,{\it t_8}\,{\it t_9} }\\\\ 
C>   {\it t_{11}} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\beta^{{{2}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- 15.19266624115207\,{\it t_{13}}\,{\it t_9} }\\\\ 
C>   {\it t_{15}} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{15}}}}\\\\ 
C>   {\it t_{17}} &=& \rho_s^{{{2}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& e^ {- 15.19266624115207\,{\it t_{17}}\,{\it t_9} }\\\\ 
C>   f &=& -.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{12}}\,\left(.2565564877060837\,{\it t_3}\,{\it t_{12}}\,\left({\it t_{14}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{14}}-1.0\right)}\over{{\it t_{13}}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_{11}}\right)\right)+{\it t_1}+1.0\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_5}\,\left(.2565564877060837\,{\it t_3}\,{\it t_5}\,\left({\it t_{10}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{10}}-1.0\right)}\over{{\it t_8}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_4}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{16}}\,\left(.2565564877060837\,{\it t_3}\,{\it t_{16}}\,\left({\it t_{18}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{18}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_{15}}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
      subroutine nwxcm_x_camlsd(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
            Amat(iq,D1_RA) = -1.8610514726981986d+0*t1*(1.14025105647148
     1         35d-1*t14*t16*t2*t3-3.4207531694144505d-1*t2*t3*t5*(2.598
     2         5180598138485d+0*t11*t6*t9+2.5655648770608375d-1*t3*t5*(2
     3         .1940410461360632d-2*t12*t7/rhoa**1.6666666666666669d+0-1
     4         .0128444160768044d+1*t10*t11*t5+3.333333333333333d-1*t11/
     5         rhoa)-8.551882923536125d-2*t13*t16*t3))*wght-2.4814019635
     6         97598d+0*t15*t4*wght+Amat(iq,D1_RA)
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
C>   {\it t_{10}} &=& e^ {- 15.19266624115207\,{\it t_8}\,{\it t_9} }\\\\ 
C>   {\it t_{11}} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\beta^{{{2}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- 15.19266624115207\,{\it t_{13}}\,{\it t_9} }\\\\ 
C>   {\it t_{15}} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{15}}}}\\\\ 
C>   {\it t_{17}} &=& \rho_s^{{{2}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& e^ {- 15.19266624115207\,{\it t_{17}}\,{\it t_9} }\\\\ 
C>   f &=& -.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{12}}\,\left(.2565564877060837\,{\it t_3}\,{\it t_{12}}\,\left({\it t_{14}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{14}}-1.0\right)}\over{{\it t_{13}}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_{11}}\right)\right)+{\it t_1}+1.0\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_5}\,\left(.2565564877060837\,{\it t_3}\,{\it t_5}\,\left({\it t_{10}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{10}}-1.0\right)}\over{{\it t_8}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_4}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{16}}\,\left(.2565564877060837\,{\it t_3}\,{\it t_{16}}\,\left({\it t_{18}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{18}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_{15}}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
      subroutine nwxcm_x_camlsd_d2(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
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
C>   {\it t_{10}} &=& e^ {- 15.19266624115207\,{\it t_8}\,{\it t_9} }\\\\ 
C>   {\it t_{11}} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& \rho_\beta^{{{2}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- 15.19266624115207\,{\it t_{13}}\,{\it t_9} }\\\\ 
C>   {\it t_{15}} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{15}}}}\\\\ 
C>   {\it t_{17}} &=& \rho_s^{{{2}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& e^ {- 15.19266624115207\,{\it t_{17}}\,{\it t_9} }\\\\ 
C>   f &=& -.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{12}}\,\left(.2565564877060837\,{\it t_3}\,{\it t_{12}}\,\left({\it t_{14}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{14}}-1.0\right)}\over{{\it t_{13}}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_{11}}\right)\right)+{\it t_1}+1.0\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_5}\,\left(.2565564877060837\,{\it t_3}\,{\it t_5}\,\left({\it t_{10}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{10}}-1.0\right)}\over{{\it t_8}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_4}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left(-0.342075316941445\,{\it t_2}\,{\it t_3}\,{\it t_{16}}\,\left(.2565564877060837\,{\it t_3}\,{\it t_{16}}\,\left({\it t_{18}}-{{.03291061569204095\,{\it t_7}\,\left({\it t_{18}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)+1.77245385090552\,\mathrm{erf}\left(3.897777089720764\,{\it t_6}\,{\it t_{15}}\right)\right)+{\it t_1}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
      subroutine nwxcm_x_camlsd_d3(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2,Amat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
