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
C>   f &=& 1.0\,\left(-0.9305257363491001\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_7}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_7}\right)^{2.0}-2.0\right)\,{\it t_7}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_\beta^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_7}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_6}\right)\right)\right)-0.9305257363491001\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_3}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_3}\right)^{2.0}-2.0\right)\,{\it t_3}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_\alpha^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_3}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_2}\right)\right)\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491001\,\rho_s^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_9}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_9}\right)^{2.0}-2.0\right)\,{\it t_9}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_s^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_9}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_8}\right)\right)\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] R Baer, D Neuhauser, E Livshits, Phys.Chem.Chem.Phys. 9, 2932 (2007)  , DOI: 
C> <a href="http://dx.doi.org/10.1039/B617919C "> 
C> 10.1039/B617919C </a> 
C> 
      subroutine nwxcm_x_bnl(param,tol_rho,ipol,nq,wght, 
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
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t13*wght
            Amat(iq,D1_RA) = -1.8610514726982003d+0*t1*(5.70125528235743
     1         1d-2*t12*t14*t2-1.7103765847072294d-1*t2*t4*(-1.125790929
     2         359309d-2*t11*t6*t8/rhoa**1.6666666666666669d+0+2.5985180
     3         59813836d+0*t10*t15*t5*t7-8.551882923536147d-2*t11*t14*t2
     4         *t7+5.197036119627672d+0*t10*t15*t5+8.551882923536147d-2*
     5         t14*t2))*wght-2.4814019635975998d+0*t13*t3*wght+Amat(iq,D
     6         1_RA)
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
C>   f &=& 1.0\,\left(-0.9305257363491001\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_7}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_7}\right)^{2.0}-2.0\right)\,{\it t_7}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_\beta^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_7}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_6}\right)\right)\right)-0.9305257363491001\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_3}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_3}\right)^{2.0}-2.0\right)\,{\it t_3}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_\alpha^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_3}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_2}\right)\right)\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491001\,\rho_s^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_9}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_9}\right)^{2.0}-2.0\right)\,{\it t_9}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_s^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_9}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_8}\right)\right)\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] R Baer, D Neuhauser, E Livshits, Phys.Chem.Chem.Phys. 9, 2932 (2007)  , DOI: 
C> <a href="http://dx.doi.org/10.1039/B617919C "> 
C> 10.1039/B617919C </a> 
C> 
      subroutine nwxcm_x_bnl_d2(param,tol_rho,ipol,nq,wght, 
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
C>   f &=& 1.0\,\left(-0.9305257363491001\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_7}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_7}\right)^{2.0}-2.0\right)\,{\it t_7}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_\beta^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_7}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_6}\right)\right)\right)-0.9305257363491001\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_3}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_3}\right)^{2.0}-2.0\right)\,{\it t_3}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_\alpha^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_3}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_2}\right)\right)\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491001\,\rho_s^{{{4}\over{3}}}\,\left(1.0-0.1710376584707229\,{\it t_1}\,{\it t_9}\,\left(0.2565564877060844\,{\it t_1}\,\left(0.06582123138408225\,\left({\it t_1}\,{\it t_9}\right)^{2.0}-2.0\right)\,{\it t_9}\,\left(1.0-{{1}\over{e^{15.19266624115199\,{\it t_5}\,\rho_s^{{{2}\over{3}}}}}}\right)-0.2565564877060844\,{\it t_1}\,{\it t_9}+3.544907701811032\,\mathrm{erf}\left(3.897777089720754\,{\it t_4}\,{\it t_8}\right)\right)\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] R Baer, D Neuhauser, E Livshits, Phys.Chem.Chem.Phys. 9, 2932 (2007)  , DOI: 
C> <a href="http://dx.doi.org/10.1039/B617919C "> 
C> 10.1039/B617919C </a> 
C> 
      subroutine nwxcm_x_bnl_d3(param,tol_rho,ipol,nq,wght, 
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
