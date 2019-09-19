C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_bnl.F 
C> The nwxcm_x_bnl functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_bnl functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^2}}\\\\ 
C>   {\it t_6} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_6}}}\\\\ 
C>   {\it t_8} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_8}}}\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363491001\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_7}
C>    \,\left(0.2565564877060844\,{\it t_1}\,
C>    \left(0.06582123138408225\,\left({\it t_1}\,{
C>    \it t_7}\right)^{2.0}-2.0\right)\,{\it t_7}\,\left(1.0
C>    -{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_\beta^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{
C>    \it t_1}\,{\it t_7}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_6}\right)\right)\right)-0.9305257363491001
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229
C>    \,{\it t_1}\,{\it t_3}\,\left(0.2565564877060844\,{\it t_1}
C>    \,\left(0.06582123138408225\,\left({\it t_1}\,{
C>    \it t_3}\right)^{2.0}-2.0\right)\,{\it t_3}\,\left(1.0
C>    -{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_\alpha^{{{2}\over{3}}}}}}\right)-0.2565564877060844
C>    \,{\it t_1}\,{\it t_3}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_2}\right)\right)\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491001\,\rho_s^{{{4}\over{3}}}\,\left(1.0
C>    -0.1710376584707229\,{\it t_1}\,{\it t_9}\,
C>    \left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225
C>    \,\left({\it t_1}\,{\it t_9}\right)^{2.0}-2.0\right)\,{
C>    \it t_9}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_s^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{
C>    \it t_1}\,{\it t_9}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_8}\right)\right)\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] R Baer, D Neuhauser, E Livshits, Phys.Chem.Chem.Phys. 9
C>    , 2932 (2007)  , DOI: 
C> <a href="https://doi.org/10.1039/B617919C "> 
C> 10.1039/B617919C </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_bnl(param,tol_rho,ipol,nq,wght, 
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
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = 6.582123138408225d-2*(t2*t4)**2.0d+0-2.0d+0
            t7 = rhoa**6.666666666666666d-1
            t8 = exp(-1.5192666241151986d+1*t7/t2**2)
            t9 = 1.0d+0-t8
            t10 = 2.565564877060844d-1*t2*t4*t6*t9+3.544907701811032d+0*
     1         erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2*
     2         t4
            t11 = 1.0d+0-1.7103765847072294d-1*t10*t2*t4
            t12 = 1/t1
            t13 = 1/t7
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*wght
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363491001d-1*t1*(5.7012552
     1         82357431d-2*t10*t12*t2-1.7103765847072294d-1*t2*t4*(-8.55
     2         1882923536147d-2*t12*t2*t6*t9-1.125790929359309d-2*t2**3*
     3         t9/rhoa**2+2.598518059813836d+0*t13*t5*t6*t8+5.1970361196
     4         27672d+0*t13*t5*t8+8.551882923536147d-2*t12*t2))-1.240700
     5         9817987999d+0*t11*t3)*wght+Amat(iq,D1_RA)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = 1/t8
            t10 = rhoa**6.666666666666666d-1
            t11 = exp(-1.5192666241151986d+1*t10*t9)
            t12 = 1.0d+0-t11
            t13 = 2.565564877060844d-1*t12*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t14 = 1.0d+0-1.7103765847072294d-1*t13*t2*t4
            t15 = rhob**1.3333333333333333d+0
            t16 = rhob**3.333333333333333d-1
            t17 = 1/t16
            t18 = t17*t2
            t19 = 6.582123138408225d-2*t18**2.0d+0-2.0d+0
            t20 = rhob**6.666666666666666d-1
            t21 = exp(-1.5192666241151986d+1*t20*t9)
            t22 = 1.0d+0-t21
            t23 = 3.544907701811032d+0*erf(3.8977770897207537d+0*t16*t5)
     1         +2.565564877060844d-1*t17*t19*t2*t22-2.565564877060844d-1
     2         *t17*t2
            t24 = 1.0d+0-1.7103765847072294d-1*t17*t2*t23
            t25 = 1/t1
            t26 = 5.641895835477563d-1
            t27 = 1/t10
            t28 = 1/t15
            t29 = 1/t20
            fnc(iq) = 1.0d+0*(-9.305257363491001d-1*t15*t24-9.3052573634
     1         91001d-1*t1*t14)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363491001d-1*t1*(5.7012552
     1         82357431d-2*t13*t2*t25-1.7103765847072294d-1*t2*t4*(-1.12
     2         5790929359309d-2*t12*t6*t8/rhoa**1.6666666666666669d+0+2.
     3         598518059813836d+0*t11*t27*t5*t7-8.551882923536147d-2*t12
     4         *t2*t25*t7+9.211506683529127d+0*t11*t26*t27*t5+8.55188292
     5         3536147d-2*t2*t25))-1.2407009817987999d+0*t14*t3)*wght+Am
     6         at(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363491001d-1*t15*(5.701255
     1         282357431d-2*t2*t23*t28-1.7103765847072294d-1*t17*t2*(-1.
     2         125790929359309d-2*t18*t22*t8/rhob**1.6666666666666669d+0
     3         +9.211506683529127d+0*t21*t26*t29*t5+2.598518059813836d+0
     4         *t19*t21*t29*t5-8.551882923536147d-2*t19*t2*t22*t28+8.551
     5         882923536147d-2*t2*t28))-1.2407009817987999d+0*t16*t24)*w
     6         ght+Amat(iq,D1_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = rhoa**6.666666666666666d-1
            t10 = exp(-1.5192666241151986d+1*t9/t8)
            t11 = 1.0d+0-t10
            t12 = 2.565564877060844d-1*t11*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t13 = 1.0d+0-1.7103765847072294d-1*t12*t2*t4
            t14 = 1/t1
            t15 = 1/t9
            fnc(iq) = fnc(iq)-9.305257363491001d-1*t1*t13*wght
            Amat(iq,D1_RA) = -9.305257363491001d-1*t1*(5.701255282357431
     1         d-2*t12*t14*t2-1.7103765847072294d-1*t2*t4*(-1.1257909293
     2         59309d-2*t11*t6*t8/rhoa**1.6666666666666669d+0+2.59851805
     3         9813836d+0*t10*t15*t5*t7-8.551882923536147d-2*t11*t14*t2*
     4         t7+5.197036119627672d+0*t10*t15*t5+8.551882923536147d-2*t
     5         14*t2))*wght-1.2407009817987999d+0*t13*t3*wght+Amat(iq,D1
     6         _RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhob**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = rhob**6.666666666666666d-1
            t10 = exp(-1.5192666241151986d+1*t9/t8)
            t11 = 1.0d+0-t10
            t12 = 2.565564877060844d-1*t11*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t13 = 1.0d+0-1.7103765847072294d-1*t12*t2*t4
            t14 = 1/t1
            t15 = 1/t9
            fnc(iq) = fnc(iq)-9.305257363491001d-1*t1*t13*wght
            Amat(iq,D1_RB) = -9.305257363491001d-1*t1*(5.701255282357431
     1         d-2*t12*t14*t2-1.7103765847072294d-1*t2*t4*(-1.1257909293
     2         59309d-2*t11*t6*t8/rhob**1.6666666666666669d+0+2.59851805
     3         9813836d+0*t10*t15*t5*t7-8.551882923536147d-2*t11*t14*t2*
     4         t7+5.197036119627672d+0*t10*t15*t5+8.551882923536147d-2*t
     5         14*t2))*wght-1.2407009817987999d+0*t13*t3*wght+Amat(iq,D1
     6         _RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_bnl functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^2}}\\\\ 
C>   {\it t_6} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_6}}}\\\\ 
C>   {\it t_8} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_8}}}\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363491001\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_7}
C>    \,\left(0.2565564877060844\,{\it t_1}\,
C>    \left(0.06582123138408225\,\left({\it t_1}\,{
C>    \it t_7}\right)^{2.0}-2.0\right)\,{\it t_7}\,\left(1.0
C>    -{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_\beta^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{
C>    \it t_1}\,{\it t_7}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_6}\right)\right)\right)-0.9305257363491001
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229
C>    \,{\it t_1}\,{\it t_3}\,\left(0.2565564877060844\,{\it t_1}
C>    \,\left(0.06582123138408225\,\left({\it t_1}\,{
C>    \it t_3}\right)^{2.0}-2.0\right)\,{\it t_3}\,\left(1.0
C>    -{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_\alpha^{{{2}\over{3}}}}}}\right)-0.2565564877060844
C>    \,{\it t_1}\,{\it t_3}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_2}\right)\right)\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491001\,\rho_s^{{{4}\over{3}}}\,\left(1.0
C>    -0.1710376584707229\,{\it t_1}\,{\it t_9}\,
C>    \left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225
C>    \,\left({\it t_1}\,{\it t_9}\right)^{2.0}-2.0\right)\,{
C>    \it t_9}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_s^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{
C>    \it t_1}\,{\it t_9}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_8}\right)\right)\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] R Baer, D Neuhauser, E Livshits, Phys.Chem.Chem.Phys. 9
C>    , 2932 (2007)  , DOI: 
C> <a href="https://doi.org/10.1039/B617919C "> 
C> 10.1039/B617919C </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_bnl_d2(param,tol_rho,ipol,nq,wght, 
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
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = 6.582123138408225d-2*(t2*t4)**2.0d+0-2.0d+0
            t7 = rhoa**6.666666666666666d-1
            t8 = exp(-1.5192666241151986d+1*t7/t2**2)
            t9 = 1.0d+0-t8
            t10 = 2.565564877060844d-1*t2*t4*t6*t9+3.544907701811032d+0*
     1         erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2*
     2         t4
            t11 = 1.0d+0-1.7103765847072294d-1*t10*t2*t4
            t12 = 1/t1
            t13 = t2**3
            t14 = 5.641895835477563d-1
            t15 = 1/t7
            t16 = -8.551882923536147d-2*t12*t2*t6*t9-1.125790929359309d-
     1         2*t13*t9/rhoa**2+2.598518059813836d+0*t15*t5*t6*t8+9.2115
     2         06683529127d+0*t14*t15*t5*t8+8.551882923536147d-2*t12*t2
            t17 = 5.701255282357431d-2*t10*t12*t2-1.7103765847072294d-1*
     1         t16*t2*t4
            t18 = 1/rhoa**2.3333333333333334d+0
            t19 = 1/rhoa**1.6666666666666669d+0
            t20 = 1/t13
            t21 = 1/rhoa
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*wght
            Amat(iq,D1_RA) = 1.0d+0*(-1.2407009817987999d+0*t11*t3-9.305
     1         257363491001d-1*t1*t17)*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363491001d-1*t1*(-1.71
     1         03765847072294d-1*t2*t4*(1.1402510564714863d-1*t18*t2*t6*
     2         t9+2.6268455018383874d-2*t13*t9/rhoa**3-2.598518059813836
     3         d+0*t19*t5*t6*t8-2.6318945069571614d+1*t20*t21*t6*t8-6.14
     4         1004455686084d+0*t14*t19*t5*t8-9.32982310806659d+1*t14*t2
     5         0*t21*t8-2.2805021129429726d-1*t18*t2*t8-1.14025105647148
     6         63d-1*t18*t2)-7.601673709809909d-2*t10*t18*t2+1.140251056
     7         4714863d-1*t12*t16*t2)-2.4814019635975998d+0*t17*t3-4.135
     8         669939329333d-1*t11*t15)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = 1/t8
            t10 = rhoa**6.666666666666666d-1
            t11 = exp(-1.5192666241151986d+1*t10*t9)
            t12 = 1.0d+0-t11
            t13 = 2.565564877060844d-1*t12*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t14 = 1.0d+0-1.7103765847072294d-1*t13*t2*t4
            t15 = rhob**1.3333333333333333d+0
            t16 = rhob**3.333333333333333d-1
            t17 = 1/t16
            t18 = t17*t2
            t19 = 6.582123138408225d-2*t18**2.0d+0-2.0d+0
            t20 = rhob**6.666666666666666d-1
            t21 = exp(-1.5192666241151986d+1*t20*t9)
            t22 = 1.0d+0-t21
            t23 = 3.544907701811032d+0*erf(3.8977770897207537d+0*t16*t5)
     1         +2.565564877060844d-1*t17*t19*t2*t22-2.565564877060844d-1
     2         *t17*t2
            t24 = 1.0d+0-1.7103765847072294d-1*t17*t2*t23
            t25 = 1/t1
            t26 = 1/rhoa**1.6666666666666669d+0
            t27 = 5.641895835477563d-1
            t28 = 1/t10
            t29 = -1.125790929359309d-2*t12*t26*t6*t8+2.598518059813836d
     1         +0*t11*t28*t5*t7-8.551882923536147d-2*t12*t2*t25*t7+9.211
     2         506683529127d+0*t11*t27*t28*t5+8.551882923536147d-2*t2*t2
     3         5
            t30 = 5.701255282357431d-2*t13*t2*t25-1.7103765847072294d-1*
     1         t2*t29*t4
            t31 = 1/t15
            t32 = 1/rhob**1.6666666666666669d+0
            t33 = 1/t20
            t34 = -1.125790929359309d-2*t18*t22*t32*t8+9.211506683529127
     1         d+0*t21*t27*t33*t5+2.598518059813836d+0*t19*t21*t33*t5-8.
     2         551882923536147d-2*t19*t2*t22*t31+8.551882923536147d-2*t2
     3         *t31
            t35 = 5.701255282357431d-2*t2*t23*t31-1.7103765847072294d-1*
     1         t17*t2*t34
            t36 = 1/rhoa**2.3333333333333334d+0
            t37 = t2**3
            t38 = 1/t37
            t39 = 1/rhoa
            t40 = 1/rhob**2.3333333333333334d+0
            t41 = 1/rhob
            fnc(iq) = 1.0d+0*(-9.305257363491001d-1*t15*t24-9.3052573634
     1         91001d-1*t1*t14)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363491001d-1*t1*t30-1.2407
     1         009817987999d+0*t14*t3)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363491001d-1*t15*t35-1.240
     1         7009817987999d+0*t16*t24)*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363491001d-1*t1*(-1.71
     1         03765847072294d-1*t2*t4*(3.752636431197697d-3*t12*t6*t8/r
     2         hoa**2.6666666666666666d+0-2.598518059813836d+0*t11*t26*t
     3         5*t7-2.6318945069571614d+1*t11*t38*t39*t7+1.1402510564714
     4         863d-1*t12*t2*t36*t7-1.1402510564714863d-1*t11*t6/rhoa**2
     5         -6.141004455686084d+0*t11*t26*t27*t5-9.32982310806659d+1*
     6         t11*t27*t38*t39+2.251581858718618d-2*t12*t37/rhoa**3-1.14
     7         02510564714861d-1*t11*t2*t36-1.1402510564714863d-1*t2*t36
     8         )-7.601673709809909d-2*t13*t2*t36+1.1402510564714863d-1*t
     9         2*t25*t29)-2.4814019635975998d+0*t3*t30-4.135669939329333
     :         d-1*t14*t28)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-9.305257363491001d-1*t15*(-1.7
     1         103765847072294d-1*t17*t2*(3.752636431197697d-3*t18*t22*t
     2         8/rhob**2.6666666666666666d+0-6.141004455686084d+0*t21*t2
     3         7*t32*t5-2.598518059813836d+0*t19*t21*t32*t5-9.3298231080
     4         6659d+1*t21*t27*t38*t41-2.6318945069571614d+1*t19*t21*t38
     5         *t41+1.1402510564714863d-1*t19*t2*t22*t40-1.1402510564714
     6         861d-1*t2*t21*t40-1.1402510564714863d-1*t2*t40+2.25158185
     7         8718618d-2*t22*t37/rhob**3-1.1402510564714863d-1*t18*t21/
     8         rhob**2)-7.601673709809909d-2*t2*t23*t40+1.14025105647148
     9         63d-1*t2*t31*t34)-2.4814019635975998d+0*t16*t35-4.1356699
     :         39329333d-1*t24*t33)*wght+Amat2(iq,D2_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = rhoa**6.666666666666666d-1
            t10 = exp(-1.5192666241151986d+1*t9/t8)
            t11 = 1.0d+0-t10
            t12 = 2.565564877060844d-1*t11*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t13 = 1.0d+0-1.7103765847072294d-1*t12*t2*t4
            t14 = 1/t1
            t15 = 1/rhoa**1.6666666666666669d+0
            t16 = 5.641895835477563d-1
            t17 = 1/t9
            t18 = -1.125790929359309d-2*t11*t15*t6*t8+2.598518059813836d
     1         +0*t10*t17*t5*t7-8.551882923536147d-2*t11*t14*t2*t7+9.211
     2         506683529127d+0*t10*t16*t17*t5+8.551882923536147d-2*t14*t
     3         2
            t19 = 5.701255282357431d-2*t12*t14*t2-1.7103765847072294d-1*
     1         t18*t2*t4
            t20 = 1/rhoa**2.3333333333333334d+0
            t21 = t2**3
            t22 = 1/t21
            t23 = 1/rhoa
            fnc(iq) = fnc(iq)-9.305257363491001d-1*t1*t13*wght
            Amat(iq,D1_RA) = -1.2407009817987999d+0*t13*t3*wght-9.305257
     1         363491001d-1*t1*t19*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = -9.305257363491001d-1*t1*(-1.7103765847
     1         072294d-1*t2*t4*(3.752636431197697d-3*t11*t6*t8/rhoa**2.6
     2         666666666666666d+0-2.598518059813836d+0*t10*t15*t5*t7-2.6
     3         318945069571614d+1*t10*t22*t23*t7+1.1402510564714863d-1*t
     4         11*t2*t20*t7-1.1402510564714863d-1*t10*t6/rhoa**2-6.14100
     5         4455686084d+0*t10*t15*t16*t5-9.32982310806659d+1*t10*t16*
     6         t22*t23+2.251581858718618d-2*t11*t21/rhoa**3-1.1402510564
     7         714861d-1*t10*t2*t20-1.1402510564714863d-1*t2*t20)-7.6016
     8         73709809909d-2*t12*t2*t20+1.1402510564714863d-1*t14*t18*t
     9         2)*wght-2.4814019635975998d+0*t19*t3*wght-4.1356699393293
     :         33d-1*t13*t17*wght+Amat2(iq,D2_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhob**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = rhob**6.666666666666666d-1
            t10 = exp(-1.5192666241151986d+1*t9/t8)
            t11 = 1.0d+0-t10
            t12 = 2.565564877060844d-1*t11*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t13 = 1.0d+0-1.7103765847072294d-1*t12*t2*t4
            t14 = 1/t1
            t15 = 1/rhob**1.6666666666666669d+0
            t16 = 5.641895835477563d-1
            t17 = 1/t9
            t18 = -1.125790929359309d-2*t11*t15*t6*t8+2.598518059813836d
     1         +0*t10*t17*t5*t7-8.551882923536147d-2*t11*t14*t2*t7+9.211
     2         506683529127d+0*t10*t16*t17*t5+8.551882923536147d-2*t14*t
     3         2
            t19 = 5.701255282357431d-2*t12*t14*t2-1.7103765847072294d-1*
     1         t18*t2*t4
            t20 = 1/rhob**2.3333333333333334d+0
            t21 = t2**3
            t22 = 1/t21
            t23 = 1/rhob
            fnc(iq) = fnc(iq)-9.305257363491001d-1*t1*t13*wght
            Amat(iq,D1_RB) = -1.2407009817987999d+0*t13*t3*wght-9.305257
     1         363491001d-1*t1*t19*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RB_RB) = -9.305257363491001d-1*t1*(-1.7103765847
     1         072294d-1*t2*t4*(3.752636431197697d-3*t11*t6*t8/rhob**2.6
     2         666666666666666d+0-2.598518059813836d+0*t10*t15*t5*t7-2.6
     3         318945069571614d+1*t10*t22*t23*t7+1.1402510564714863d-1*t
     4         11*t2*t20*t7-1.1402510564714863d-1*t10*t6/rhob**2-6.14100
     5         4455686084d+0*t10*t15*t16*t5-9.32982310806659d+1*t10*t16*
     6         t22*t23+2.251581858718618d-2*t11*t21/rhob**3-1.1402510564
     7         714861d-1*t10*t2*t20-1.1402510564714863d-1*t2*t20)-7.6016
     8         73709809909d-2*t12*t2*t20+1.1402510564714863d-1*t14*t18*t
     9         2)*wght-2.4814019635975998d+0*t19*t3*wght-4.1356699393293
     :         33d-1*t13*t17*wght+Amat2(iq,D2_RB_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_bnl functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^2}}\\\\ 
C>   {\it t_6} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_6}}}\\\\ 
C>   {\it t_8} &=& \rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_8}}}\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363491001\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_7}
C>    \,\left(0.2565564877060844\,{\it t_1}\,
C>    \left(0.06582123138408225\,\left({\it t_1}\,{
C>    \it t_7}\right)^{2.0}-2.0\right)\,{\it t_7}\,\left(1.0
C>    -{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_\beta^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{
C>    \it t_1}\,{\it t_7}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_6}\right)\right)\right)-0.9305257363491001
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229
C>    \,{\it t_1}\,{\it t_3}\,\left(0.2565564877060844\,{\it t_1}
C>    \,\left(0.06582123138408225\,\left({\it t_1}\,{
C>    \it t_3}\right)^{2.0}-2.0\right)\,{\it t_3}\,\left(1.0
C>    -{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_\alpha^{{{2}\over{3}}}}}}\right)-0.2565564877060844
C>    \,{\it t_1}\,{\it t_3}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_2}\right)\right)\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491001\,\rho_s^{{{4}\over{3}}}\,\left(1.0
C>    -0.1710376584707229\,{\it t_1}\,{\it t_9}\,
C>    \left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225
C>    \,\left({\it t_1}\,{\it t_9}\right)^{2.0}-2.0\right)\,{
C>    \it t_9}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}
C>    \,\rho_s^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{
C>    \it t_1}\,{\it t_9}+3.544907701811032\,\mathrm{erf}
C>    \left(3.897777089720754\,{\it t_4}\,{
C>    \it t_8}\right)\right)\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] R Baer, D Neuhauser, E Livshits, Phys.Chem.Chem.Phys. 9
C>    , 2932 (2007)  , DOI: 
C> <a href="https://doi.org/10.1039/B617919C "> 
C> 10.1039/B617919C </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_bnl_d3(param,tol_rho,ipol,nq,wght, 
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
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = 6.582123138408225d-2*(t2*t4)**2.0d+0-2.0d+0
            t7 = rhoa**6.666666666666666d-1
            t8 = exp(-1.5192666241151986d+1*t7/t2**2)
            t9 = 1.0d+0-t8
            t10 = 2.565564877060844d-1*t2*t4*t6*t9+3.544907701811032d+0*
     1         erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2*
     2         t4
            t11 = 1.0d+0-1.7103765847072294d-1*t10*t2*t4
            t12 = 1/t1
            t13 = t2**3
            t14 = 1/rhoa**2
            t15 = 5.641895835477563d-1
            t16 = 1/t7
            t17 = -8.551882923536147d-2*t12*t2*t6*t9-1.125790929359309d-
     1         2*t13*t14*t9+2.598518059813836d+0*t16*t5*t6*t8+9.21150668
     2         3529127d+0*t15*t16*t5*t8+8.551882923536147d-2*t12*t2
            t18 = 5.701255282357431d-2*t10*t12*t2-1.7103765847072294d-1*
     1         t17*t2*t4
            t19 = 1/rhoa**2.3333333333333334d+0
            t20 = 1/rhoa**1.6666666666666669d+0
            t21 = 1/t13
            t22 = 1/rhoa
            t23 = 1.1402510564714863d-1*t19*t2*t6*t9+2.6268455018383874d
     1         -2*t13*t9/rhoa**3-2.598518059813836d+0*t20*t5*t6*t8-2.631
     2         8945069571614d+1*t21*t22*t6*t8-6.141004455686084d+0*t15*t
     3         20*t5*t8-9.32982310806659d+1*t15*t21*t22*t8-2.28050211294
     4         29726d-1*t19*t2*t8-1.1402510564714863d-1*t19*t2
            t24 = -1.7103765847072294d-1*t2*t23*t4-7.601673709809909d-2*
     1         t10*t19*t2+1.1402510564714863d-1*t12*t17*t2
            t25 = 1/rhoa**3.3333333333333337d+0
            t26 = 1/rhoa**2.6666666666666666d+0
            t27 = 1/t2**5
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*wght
            Amat(iq,D1_RA) = 1.0d+0*(-1.2407009817987999d+0*t11*t3-9.305
     1         257363491001d-1*t1*t18)*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-2.4814019635975998d+0*t18*t3-9
     1         .305257363491001d-1*t1*t24-4.135669939329333d-1*t11*t16)*
     2         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-9.305257363491001d-1*t1*(-1
     1         .7103765847072294d-1*t2*t4*(-2.660585798433468d-1*t2*t25*
     2         t6*t9-8.380888029674854d-2*t13*t9/rhoa**4+5.4857603484958
     3         75d+0*t26*t5*t6*t8+2.665699655074761d+2*t12*t27*t6*t8+5.2
     4         63789013914323d+1*t14*t21*t6*t8+1.0235007426143473d+1*t15
     5         *t26*t5*t8+3.4646907464184473d+0*t26*t5*t8+9.449659237989
     6         532d+2*t12*t15*t27*t8+9.12200845177189d-1*t2*t25*t8+1.554
     7         9705180110981d+2*t14*t15*t21*t8+2.660585798433468d-1*t2*t
     8         25)+1.773723865622312d-1*t10*t2*t25+1.7103765847072294d-1
     9         *t12*t2*t23-2.2805021129429726d-1*t17*t19*t2)-3.722102945
     :         3964d+0*t24*t3+2.757113292886222d-1*t11*t20-1.24070098179
     ;         87999d+0*t16*t18)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = 1/t8
            t10 = rhoa**6.666666666666666d-1
            t11 = exp(-1.5192666241151986d+1*t10*t9)
            t12 = 1.0d+0-t11
            t13 = 2.565564877060844d-1*t12*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t14 = 1.0d+0-1.7103765847072294d-1*t13*t2*t4
            t15 = rhob**1.3333333333333333d+0
            t16 = rhob**3.333333333333333d-1
            t17 = 1/t16
            t18 = t17*t2
            t19 = 6.582123138408225d-2*t18**2.0d+0-2.0d+0
            t20 = rhob**6.666666666666666d-1
            t21 = exp(-1.5192666241151986d+1*t20*t9)
            t22 = 1.0d+0-t21
            t23 = 3.544907701811032d+0*erf(3.8977770897207537d+0*t16*t5)
     1         +2.565564877060844d-1*t17*t19*t2*t22-2.565564877060844d-1
     2         *t17*t2
            t24 = 1.0d+0-1.7103765847072294d-1*t17*t2*t23
            t25 = 1/t1
            t26 = 5.701255282357431d-2*t13*t2*t25
            t27 = 8.551882923536147d-2*t2*t25
            t28 = 1/rhoa**1.6666666666666669d+0
            t29 = -8.551882923536147d-2*t12*t2*t25*t7
            t30 = 5.641895835477563d-1
            t31 = 1/t10
            t32 = 9.211506683529127d+0*t11*t30*t31*t5
            t33 = 2.598518059813836d+0*t11*t31*t5*t7
            t34 = -1.125790929359309d-2*t12*t28*t6*t8+t33+t32+t29+t27
            t35 = t26-1.7103765847072294d-1*t2*t34*t4
            t36 = 1/t15
            t37 = 5.701255282357431d-2*t2*t23*t36
            t38 = 8.551882923536147d-2*t2*t36
            t39 = 1/rhob**1.6666666666666669d+0
            t40 = -8.551882923536147d-2*t19*t2*t22*t36
            t41 = 1/t20
            t42 = 9.211506683529127d+0*t21*t30*t41*t5
            t43 = 2.598518059813836d+0*t19*t21*t41*t5
            t44 = -1.125790929359309d-2*t18*t22*t39*t8+t43+t42+t40+t38
            t45 = t37-1.7103765847072294d-1*t17*t2*t44
            t46 = 1/rhoa**2.3333333333333334d+0
            t47 = t2**3
            t48 = 1/rhoa**3
            t49 = 1/rhoa**2.6666666666666666d+0
            t50 = 1/rhoa**2
            t51 = 1/t47
            t52 = 1/rhoa
            t53 = 3.752636431197697d-3*t12*t49*t6*t8-2.6318945069571614d
     1         +1*t11*t51*t52*t7-2.598518059813836d+0*t11*t28*t5*t7+1.14
     2         02510564714863d-1*t12*t2*t46*t7-1.1402510564714863d-1*t11
     3         *t50*t6-9.32982310806659d+1*t11*t30*t51*t52-6.14100445568
     4         6084d+0*t11*t28*t30*t5+2.251581858718618d-2*t12*t47*t48-1
     5         .1402510564714861d-1*t11*t2*t46-1.1402510564714863d-1*t2*
     6         t46
            t54 = -1.7103765847072294d-1*t2*t4*t53-7.601673709809909d-2*
     1         t13*t2*t46+1.1402510564714863d-1*t2*t25*t34
            t55 = 1/rhob**2.3333333333333334d+0
            t56 = 1/rhob**3
            t57 = 1/rhob**2.6666666666666666d+0
            t58 = 1/rhob**2
            t59 = 1/rhob
            t60 = 3.752636431197697d-3*t18*t22*t57*t8-9.32982310806659d+
     1         1*t21*t30*t51*t59-2.6318945069571614d+1*t19*t21*t51*t59-1
     2         .1402510564714863d-1*t18*t21*t58+2.251581858718618d-2*t22
     3         *t47*t56+1.1402510564714863d-1*t19*t2*t22*t55-1.140251056
     4         4714861d-1*t2*t21*t55-1.1402510564714863d-1*t2*t55-6.1410
     5         04455686084d+0*t21*t30*t39*t5-2.598518059813836d+0*t19*t2
     6         1*t39*t5
            t61 = -1.7103765847072294d-1*t17*t2*t60-7.601673709809909d-2
     1         *t2*t23*t55+1.1402510564714863d-1*t2*t36*t44
            t62 = 1/rhoa**3.3333333333333337d+0
            t63 = 1/t2**5
            t64 = 1/rhob**3.3333333333333337d+0
            fnc(iq) = 1.0d+0*(-9.305257363491001d-1*t15*t24-9.3052573634
     1         91001d-1*t1*t14)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363491001d-1*t1*t35-1.2407
     1         009817987999d+0*t14*t3)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363491001d-1*t15*t45-1.240
     1         7009817987999d+0*t16*t24)*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363491001d-1*t1*t54-2.
     1         4814019635975998d+0*t3*t35-4.135669939329333d-1*t14*t31)*
     2         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-9.305257363491001d-1*t15*t61-2
     1         .4814019635975998d+0*t16*t45-4.135669939329333d-1*t24*t41
     2         )*wght+Amat2(iq,D2_RB_RB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-9.305257363491001d-1*t1*(-1
     1         .7103765847072294d-1*t2*t4*(1.1548969154728157d+0*t11*t46
     2         *t6*t9-5.003515241596929d-3*t12*t6*t8/rhoa**3.66666666666
     3         66664d+0+2.665699655074761d+2*t11*t25*t63*t7-2.6605857984
     4         33468d-1*t12*t2*t62*t7+5.263789013914323d+1*t11*t50*t51*t
     5         7+5.485760348495875d+0*t11*t49*t5*t7+9.449659237989532d+2
     6         *t11*t25*t30*t63+7.981757395300404d-1*t11*t2*t62+2.660585
     7         798433468d-1*t2*t62+1.1402510564714863d-1*t11*t48*t6+1.55
     8         49705180110981d+2*t11*t30*t50*t51+1.0235007426143473d+1*t
     9         11*t30*t49*t5+2.3097938309456315d+0*t11*t49*t5-7.88053650
     :         5515162d-2*t12*t47/rhoa**4)+1.773723865622312d-1*t13*t2*t
     ;         62+1.7103765847072294d-1*t2*t25*t53-2.2805021129429726d-1
     <         *t2*t34*t46)-3.7221029453964d+0*t3*t54-1.2407009817987999
     =         d+0*t31*(t26-1.7103765847072294d-1*t2*t4*(-1.125790929359
     >         309d-2*t12*t47*t50+t33+t32+t29+t27))+2.757113292886222d-1
     ?         *t14*t28)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(-9.305257363491001d-1*t15*(-
     1         1.7103765847072294d-1*t17*t2*(1.1548969154728157d+0*t18*t
     2         21*t55*t9-5.003515241596929d-3*t18*t22*t8/rhob**3.6666666
     3         666666664d+0-2.660585798433468d-1*t19*t2*t22*t64+7.981757
     4         395300404d-1*t2*t21*t64+2.660585798433468d-1*t2*t64+9.449
     5         659237989532d+2*t21*t30*t36*t63+2.665699655074761d+2*t19*
     6         t21*t36*t63+1.5549705180110981d+2*t21*t30*t51*t58+5.26378
     7         9013914323d+1*t19*t21*t51*t58+1.0235007426143473d+1*t21*t
     8         30*t5*t57+5.485760348495875d+0*t19*t21*t5*t57+2.309793830
     9         9456315d+0*t21*t5*t57+1.1402510564714863d-1*t18*t21*t56-7
     :         .880536505515162d-2*t22*t47/rhob**4)+1.773723865622312d-1
     ;         *t2*t23*t64+1.7103765847072294d-1*t2*t36*t60-2.2805021129
     <         429726d-1*t2*t44*t55)-3.7221029453964d+0*t16*t61-1.240700
     =         9817987999d+0*t41*(t37-1.7103765847072294d-1*t17*t2*(-1.1
     >         25790929359309d-2*t22*t47*t58+t43+t42+t40+t38))+2.7571132
     ?         92886222d-1*t24*t39)*wght+Amat3(iq,D3_RB_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = 1/t8
            t10 = rhoa**6.666666666666666d-1
            t11 = exp(-1.5192666241151986d+1*t10*t9)
            t12 = 1.0d+0-t11
            t13 = 2.565564877060844d-1*t12*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t14 = 1.0d+0-1.7103765847072294d-1*t13*t2*t4
            t15 = 1/t1
            t16 = 5.701255282357431d-2*t13*t15*t2
            t17 = 8.551882923536147d-2*t15*t2
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = -8.551882923536147d-2*t12*t15*t2*t7
            t20 = 5.641895835477563d-1
            t21 = 1/t10
            t22 = 9.211506683529127d+0*t11*t20*t21*t5
            t23 = 2.598518059813836d+0*t11*t21*t5*t7
            t24 = -1.125790929359309d-2*t12*t18*t6*t8+t23+t22+t19+t17
            t25 = t16-1.7103765847072294d-1*t2*t24*t4
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = t2**3
            t28 = 1/rhoa**3
            t29 = 1/rhoa**2.6666666666666666d+0
            t30 = 1/rhoa**2
            t31 = 1/t27
            t32 = 1/rhoa
            t33 = 3.752636431197697d-3*t12*t29*t6*t8-2.598518059813836d+
     1         0*t11*t18*t5*t7-2.6318945069571614d+1*t11*t31*t32*t7+1.14
     2         02510564714863d-1*t12*t2*t26*t7-1.1402510564714863d-1*t11
     3         *t30*t6-6.141004455686084d+0*t11*t18*t20*t5-9.32982310806
     4         659d+1*t11*t20*t31*t32+2.251581858718618d-2*t12*t27*t28-1
     5         .1402510564714861d-1*t11*t2*t26-1.1402510564714863d-1*t2*
     6         t26
            t34 = -1.7103765847072294d-1*t2*t33*t4-7.601673709809909d-2*
     1         t13*t2*t26+1.1402510564714863d-1*t15*t2*t24
            t35 = 1/rhoa**3.3333333333333337d+0
            t36 = 1/t2**5
            fnc(iq) = fnc(iq)-9.305257363491001d-1*t1*t14*wght
            Amat(iq,D1_RA) = -1.2407009817987999d+0*t14*t3*wght-9.305257
     1         363491001d-1*t1*t25*wght+Amat(iq,D1_RA)
            Amat2(iq,D2_RA_RA) = -9.305257363491001d-1*t1*t34*wght-2.481
     1         4019635975998d+0*t25*t3*wght-4.135669939329333d-1*t14*t21
     2         *wght+Amat2(iq,D2_RA_RA)
            Amat3(iq,D3_RA_RA_RA) = -9.305257363491001d-1*t1*(-1.7103765
     1         847072294d-1*t2*t4*(1.1548969154728157d+0*t11*t26*t6*t9-5
     2         .003515241596929d-3*t12*t6*t8/rhoa**3.6666666666666664d+0
     3         +5.485760348495875d+0*t11*t29*t5*t7+2.665699655074761d+2*
     4         t11*t15*t36*t7-2.660585798433468d-1*t12*t2*t35*t7+5.26378
     5         9013914323d+1*t11*t30*t31*t7+1.1402510564714863d-1*t11*t2
     6         8*t6+1.0235007426143473d+1*t11*t20*t29*t5+2.3097938309456
     7         315d+0*t11*t29*t5+9.449659237989532d+2*t11*t15*t20*t36+7.
     8         981757395300404d-1*t11*t2*t35+2.660585798433468d-1*t2*t35
     9         +1.5549705180110981d+2*t11*t20*t30*t31-7.880536505515162d
     :         -2*t12*t27/rhoa**4)+1.773723865622312d-1*t13*t2*t35+1.710
     ;         3765847072294d-1*t15*t2*t33-2.2805021129429726d-1*t2*t24*
     <         t26)*wght-1.2407009817987999d+0*t21*(t16-1.71037658470722
     =         94d-1*t2*(-1.125790929359309d-2*t12*t27*t30+t23+t22+t19+t
     >         17)*t4)*wght-3.7221029453964d+0*t3*t34*wght+2.75711329288
     ?         6222d-1*t14*t18*wght+Amat3(iq,D3_RA_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(1)
            t3 = rhob**3.333333333333333d-1
            t4 = 1/t3
            t5 = 1/t2
            t6 = t2*t4
            t7 = 6.582123138408225d-2*t6**2.0d+0-2.0d+0
            t8 = t2**2
            t9 = 1/t8
            t10 = rhob**6.666666666666666d-1
            t11 = exp(-1.5192666241151986d+1*t10*t9)
            t12 = 1.0d+0-t11
            t13 = 2.565564877060844d-1*t12*t2*t4*t7+3.544907701811032d+0
     1         *erf(3.8977770897207537d+0*t3*t5)-2.565564877060844d-1*t2
     2         *t4
            t14 = 1.0d+0-1.7103765847072294d-1*t13*t2*t4
            t15 = 1/t1
            t16 = 5.701255282357431d-2*t13*t15*t2
            t17 = 8.551882923536147d-2*t15*t2
            t18 = 1/rhob**1.6666666666666669d+0
            t19 = -8.551882923536147d-2*t12*t15*t2*t7
            t20 = 5.641895835477563d-1
            t21 = 1/t10
            t22 = 9.211506683529127d+0*t11*t20*t21*t5
            t23 = 2.598518059813836d+0*t11*t21*t5*t7
            t24 = -1.125790929359309d-2*t12*t18*t6*t8+t23+t22+t19+t17
            t25 = t16-1.7103765847072294d-1*t2*t24*t4
            t26 = 1/rhob**2.3333333333333334d+0
            t27 = t2**3
            t28 = 1/rhob**3
            t29 = 1/rhob**2.6666666666666666d+0
            t30 = 1/rhob**2
            t31 = 1/t27
            t32 = 1/rhob
            t33 = 3.752636431197697d-3*t12*t29*t6*t8-2.598518059813836d+
     1         0*t11*t18*t5*t7-2.6318945069571614d+1*t11*t31*t32*t7+1.14
     2         02510564714863d-1*t12*t2*t26*t7-1.1402510564714863d-1*t11
     3         *t30*t6-6.141004455686084d+0*t11*t18*t20*t5-9.32982310806
     4         659d+1*t11*t20*t31*t32+2.251581858718618d-2*t12*t27*t28-1
     5         .1402510564714861d-1*t11*t2*t26-1.1402510564714863d-1*t2*
     6         t26
            t34 = -1.7103765847072294d-1*t2*t33*t4-7.601673709809909d-2*
     1         t13*t2*t26+1.1402510564714863d-1*t15*t2*t24
            t35 = 1/rhob**3.3333333333333337d+0
            t36 = 1/t2**5
            fnc(iq) = fnc(iq)-9.305257363491001d-1*t1*t14*wght
            Amat(iq,D1_RB) = -1.2407009817987999d+0*t14*t3*wght-9.305257
     1         363491001d-1*t1*t25*wght+Amat(iq,D1_RB)
            Amat2(iq,D2_RB_RB) = -9.305257363491001d-1*t1*t34*wght-2.481
     1         4019635975998d+0*t25*t3*wght-4.135669939329333d-1*t14*t21
     2         *wght+Amat2(iq,D2_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = -9.305257363491001d-1*t1*(-1.7103765
     1         847072294d-1*t2*t4*(1.1548969154728157d+0*t11*t26*t6*t9-5
     2         .003515241596929d-3*t12*t6*t8/rhob**3.6666666666666664d+0
     3         +5.485760348495875d+0*t11*t29*t5*t7+2.665699655074761d+2*
     4         t11*t15*t36*t7-2.660585798433468d-1*t12*t2*t35*t7+5.26378
     5         9013914323d+1*t11*t30*t31*t7+1.1402510564714863d-1*t11*t2
     6         8*t6+1.0235007426143473d+1*t11*t20*t29*t5+2.3097938309456
     7         315d+0*t11*t29*t5+9.449659237989532d+2*t11*t15*t20*t36+7.
     8         981757395300404d-1*t11*t2*t35+2.660585798433468d-1*t2*t35
     9         +1.5549705180110981d+2*t11*t20*t30*t31-7.880536505515162d
     :         -2*t12*t27/rhob**4)+1.773723865622312d-1*t13*t2*t35+1.710
     ;         3765847072294d-1*t15*t2*t33-2.2805021129429726d-1*t2*t24*
     <         t26)*wght-1.2407009817987999d+0*t21*(t16-1.71037658470722
     =         94d-1*t2*(-1.125790929359309d-2*t12*t27*t30+t23+t22+t19+t
     >         17)*t4)*wght-3.7221029453964d+0*t3*t34*wght+2.75711329288
     ?         6222d-1*t14*t18*wght+Amat3(iq,D3_RB_RB_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
