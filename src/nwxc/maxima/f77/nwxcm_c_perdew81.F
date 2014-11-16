C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_perdew81.F 
C> The nwxcm_c_perdew81 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_perdew81 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_1}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.8292885914166397\,{\it t_3}+0.2068248536658599\,{\it t_2}+1.0}}\\\\ 
C>   {\it t_5} &=& 0.6203504908994\,{\it t_2}\\\\ 
C>   {\it t_6} &=& {\it nwxcm\_heaviside}\left({\it t_5}-1.0\right)\\\\ 
C>   {\it t_7} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994\,{\it t_2}\right)\\\\ 
C>   {\it t_8} &=& \log {\it t_5}\\\\ 
C>   {\it t_9} &=& 0.0012407009817988\,{\it t_2}\,{\it t_8}+0.0311\,{\it t_8}-0.007196065694433039\,{\it t_2}-0.048\\\\ 
C>   {\it t_{10}} &=& -\rho_\beta\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}+\rho_\alpha\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{13}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_{13}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{13}}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{0.8292885914166397\,{\it t_{15}}+0.2068248536658599\,{\it t_{14}}+1.0}}\\\\ 
C>   {\it t_{17}} &=& 0.6203504908994\,{\it t_{14}}\\\\ 
C>   {\it t_{18}} &=& {\it nwxcm\_heaviside}\left({\it t_{17}}-1.0\right)\\\\ 
C>   {\it t_{19}} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994\,{\it t_{14}}\right)\\\\ 
C>   {\it t_{20}} &=& \log {\it t_{17}}\\\\ 
C>   {\it t_{21}} &=& 0.0012407009817988\,{\it t_{14}}\,{\it t_{20}}+0.0311\,{\it t_{20}}-0.007196065694433039\,{\it t_{14}}-0.048\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({\it t_{11}}\,{\it t_{12}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{11}}\,{\it t_{12}}\right)^{{{4}\over{3}}}-2.0\right)\,\left(-{\it t_7}\,{\it t_9}+{\it t_7}\,\left(4.3424534362958 \times 10^{-4}\,{\it t_2}\,{\it t_8}+0.01555\,{\it t_8}-0.00297768235631712\,{\it t_2}-0.0269\right)-{{0.0843\,{\it t_6}}\over{1.101176160755631\,{\it t_3}+0.1619735131738333\,{\it t_2}+1.0}}+0.1423\,{\it t_4}\,{\it t_6}\right)+{\it t_7}\,{\it t_9}-0.1423\,{\it t_4}\,{\it t_6}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{13}}\,\left(1.923661050931536\,\left(\left({{\rho_s+{\it t_{10}}}\over{{\it t_{13}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left(-{\it t_{19}}\,{\it t_{21}}+{\it t_{19}}\,\left(4.3424534362958 \times 10^{-4}\,{\it t_{14}}\,{\it t_{20}}+0.01555\,{\it t_{20}}-0.00297768235631712\,{\it t_{14}}-0.0269\right)-{{0.0843\,{\it t_{18}}}\over{1.101176160755631\,{\it t_{15}}+0.1619735131738333\,{\it t_{14}}+1.0}}+0.1423\,{\it t_{16}}\,{\it t_{18}}\right)+{\it t_{19}}\,{\it t_{21}}-0.1423\,{\it t_{16}}\,{\it t_{18}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, A Zunger, Phys.Rev. B 23, 5048 (1981)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.23.5048 "> 
C> 10.1103/PhysRevB.23.5048 </a> 
C> 
      subroutine nwxcm_c_perdew81(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 7.937005259840998d-1
            t2 = 1/rhoa**3.333333333333333d-1
            t3 = 6.203504908994d-1*t1*t2
            t4 = nwxcm_heaviside(t3-1.0d+0)
            t5 = 2.0682485366585993d-1*t1*t2+7.388121430614919d-1/rhoa**
     1         1.6666666666666666d-1+1.0d+0
            t6 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1*t2)
            t7 = log(t3)
            t8 = t6*(1.2407009817988d-3*t1*t2*t7+3.11d-2*t7-7.1960656944
     1         33039d-3*t1*t2-4.8d-2)-1.423d-1*t4/t5
            t9 = 1/rhoa**1.3333333333333333d+0
            fnc(iq) = 2.0d+0*rhoa*t8*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*rhoa*(t6*(-4.1356699393293334d-4*t1*
     1         t7*t9+1.9851215708780795d-3*t1*t9-1.0366666666666666d-2/r
     2         hoa)+1.423d-1*t4*(-6.894161788861997d-2*t1*t9-1.231353571
     3         7691529d-1/rhoa**1.1666666666666667d+0)/t5**2)*wght+2.0d+
     4         0*t8*wght+Amat(iq,D1_RA)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob+rhoa
            t2 = 1/t1**3.333333333333333d-1
            t3 = 1/t1**1.6666666666666666d-1
            t4 = 8.292885914166397d-1*t3+2.0682485366585993d-1*t2+1.0d+0
            t5 = 1/t4
            t6 = 6.203504908994d-1*t2
            t7 = nwxcm_heaviside(t6-1.0d+0)
            t8 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t2)
            t9 = log(t6)
            t10 = 1.2407009817988d-3*t2*t9+3.11d-2*t9-7.196065694433039d
     1         -3*t2-4.8d-2
            t11 = rhoa-rhob
            t12 = 1/t1
            t13 = 1.0d+0-t11*t12
            t14 = t11*t12+1.0d+0
            t15 = t14**1.3333333333333333d+0+t13**1.3333333333333333d+0-
     1         2.0d+0
            t16 = 1.101176160755631d+0*t3+1.6197351317383332d-1*t2+1.0d+
     1         0
            t17 = t8*(4.3424534362957995d-4*t2*t9+1.555d-2*t9-2.97768235
     1         63171195d-3*t2-2.69d-2)-t10*t8+1.423d-1*t5*t7-8.43d-2*t7/
     2         t16
            t18 = t10*t8-1.423d-1*t5*t7+1.923661050931536d+0*t15*t17
            t19 = 1/t1**1.3333333333333333d+0
            t20 = 1/t1**1.1666666666666667d+0
            t21 = -1.382147652361066d-1*t20-6.894161788861997d-2*t19
            t22 = 1/t4**2
            t23 = 1.423d-1*t21*t22*t7
            t24 = -4.1356699393293334d-4*t19*t9+1.9851215708780795d-3*t1
     1         9-1.0366666666666666d-2*t12
            t25 = t24*t8
            t26 = 1.923661050931536d+0*t15*(t8*(-1.4474844787652666d-4*t
     1         19*t9+8.478123375625131d-4*t19-5.183333333333333d-3*t12)-
     2         t24*t8-1.423d-1*t21*t22*t7+8.43d-2*(-1.8352936012593848d-
     3         1*t20-5.399117105794444d-2*t19)*t7/t16**2)
            t27 = 1/t1**2
            t28 = t11*t27
            t29 = -t12
            t30 = t13**3.333333333333333d-1
            t31 = -t11*t27
            t32 = t14**3.333333333333333d-1
            t33 = 1.0d+0*t18*wght
            fnc(iq) = 1.0d+0*t1*t18*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*(1.923661050931536d+0*t17*(1.3333
     1         333333333333d+0*(t31+t12)*t32+1.3333333333333333d+0*(t29+
     2         t28)*t30)+t26+t25+t23)*wght+t33+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*t1*(1.923661050931536d+0*t17*(1.3333
     1         333333333333d+0*(t31+t29)*t32+1.3333333333333333d+0*(t28+
     2         t12)*t30)+t26+t25+t23)*wght+t33+Amat(iq,D1_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhob+rhoa
            t2 = 1/t1**3.333333333333333d-1
            t3 = 1/t1**1.6666666666666666d-1
            t4 = 8.292885914166397d-1*t3+2.0682485366585993d-1*t2+1.0d+0
            t5 = 1/t4
            t6 = 6.203504908994d-1*t2
            t7 = nwxcm_heaviside(t6-1.0d+0)
            t8 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t2)
            t9 = log(t6)
            t10 = 1.2407009817988d-3*t2*t9+3.11d-2*t9-7.196065694433039d
     1         -3*t2-4.8d-2
            t11 = rhoa-rhob
            t12 = 1/t1
            t13 = t11*t12+1.0d+0
            t14 = t13**1.3333333333333333d+0-2.0d+0
            t15 = 1.101176160755631d+0*t3+1.6197351317383332d-1*t2+1.0d+
     1         0
            t16 = t8*(4.3424534362957995d-4*t2*t9+1.555d-2*t9-2.97768235
     1         63171195d-3*t2-2.69d-2)-t10*t8+1.423d-1*t5*t7-8.43d-2*t7/
     2         t15
            t17 = t10*t8-1.423d-1*t5*t7+1.923661050931536d+0*t14*t16
            t18 = 1/t1**1.3333333333333333d+0
            t19 = 1/t1**1.1666666666666667d+0
            t20 = -1.382147652361066d-1*t19-6.894161788861997d-2*t18
            t21 = 1/t4**2
            t22 = -4.1356699393293334d-4*t18*t9+1.9851215708780795d-3*t1
     1         8-1.0366666666666666d-2*t12
            fnc(iq) = 1.0d+0*t1*t17*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*(1.923661050931536d+0*t14*(t8*(-1
     1         .4474844787652666d-4*t18*t9+8.478123375625131d-4*t18-5.18
     2         3333333333333d-3*t12)-t22*t8-1.423d-1*t20*t21*t7+8.43d-2*
     3         (-1.8352936012593848d-1*t19-5.399117105794444d-2*t18)*t7/
     4         t15**2)+t22*t8+1.423d-1*t20*t21*t7+2.5648814012420484d+0*
     5         (t12-t11/t1**2)*t13**3.333333333333333d-1*t16)*wght+1.0d+
     6         0*t17*wght+Amat(iq,D1_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 7.937005259840998d-1
            t2 = 1/rhob**3.333333333333333d-1
            t3 = 6.203504908994d-1*t1*t2
            t4 = nwxcm_heaviside(t3-1.0d+0)
            t5 = 8.908987181403393d-1
            t6 = 1/rhob**1.6666666666666666d-1
            t7 = 8.292885914166397d-1*t5*t6+2.0682485366585993d-1*t1*t2+
     1         1.0d+0
            t8 = 1/t7
            t9 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1*t2)
            t10 = log(t3)
            t11 = 1.2407009817988d-3*t1*t10*t2-7.196065694433039d-3*t1*t
     1         2+3.11d-2*t10-4.8d-2
            t12 = 1.101176160755631d+0*t5*t6+1.6197351317383332d-1*t1*t2
     1         +1.0d+0
            t13 = -1.923661050931536d+0*((4.3424534362957995d-4*t1*t10*t
     1         2-2.9776823563171195d-3*t1*t2+1.555d-2*t10-2.69d-2)*t9-t1
     2         1*t9+1.423d-1*t4*t8-8.43d-2*t4/t12)+t11*t9-1.423d-1*t4*t8
            t14 = 1/rhob**1.3333333333333333d+0
            t15 = 4.4544935907016964d-1
            t16 = 1/rhob**1.1666666666666667d+0
            t17 = -2.764295304722132d-1*t15*t16-6.894161788861997d-2*t1*
     1         t14
            t18 = 1/t7**2
            t19 = 1/rhob
            t20 = -1.0366666666666666d-2*t19-4.1356699393293334d-4*t1*t1
     1         0*t14+1.9851215708780795d-3*t1*t14
            fnc(iq) = 2.0d+0*rhob*t13*wght+fnc(iq)
            Amat(iq,D1_RB) = 2.0d+0*rhob*(-1.923661050931536d+0*(-t20*t9
     1         +(-5.183333333333333d-3*t19-1.4474844787652666d-4*t1*t10*
     2         t14+8.478123375625131d-4*t1*t14)*t9-1.423d-1*t17*t18*t4+8
     3         .43d-2*(-3.6705872025187697d-1*t15*t16-5.399117105794444d
     4         -2*t1*t14)*t4/t12**2)+t20*t9+1.423d-1*t17*t18*t4)*wght+2.
     5         0d+0*t13*wght+Amat(iq,D1_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_perdew81 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_1}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.8292885914166397\,{\it t_3}+0.2068248536658599\,{\it t_2}+1.0}}\\\\ 
C>   {\it t_5} &=& 0.6203504908994\,{\it t_2}\\\\ 
C>   {\it t_6} &=& {\it nwxcm\_heaviside}\left({\it t_5}-1.0\right)\\\\ 
C>   {\it t_7} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994\,{\it t_2}\right)\\\\ 
C>   {\it t_8} &=& \log {\it t_5}\\\\ 
C>   {\it t_9} &=& 0.0012407009817988\,{\it t_2}\,{\it t_8}+0.0311\,{\it t_8}-0.007196065694433039\,{\it t_2}-0.048\\\\ 
C>   {\it t_{10}} &=& -\rho_\beta\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}+\rho_\alpha\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{13}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_{13}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{13}}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{0.8292885914166397\,{\it t_{15}}+0.2068248536658599\,{\it t_{14}}+1.0}}\\\\ 
C>   {\it t_{17}} &=& 0.6203504908994\,{\it t_{14}}\\\\ 
C>   {\it t_{18}} &=& {\it nwxcm\_heaviside}\left({\it t_{17}}-1.0\right)\\\\ 
C>   {\it t_{19}} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994\,{\it t_{14}}\right)\\\\ 
C>   {\it t_{20}} &=& \log {\it t_{17}}\\\\ 
C>   {\it t_{21}} &=& 0.0012407009817988\,{\it t_{14}}\,{\it t_{20}}+0.0311\,{\it t_{20}}-0.007196065694433039\,{\it t_{14}}-0.048\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({\it t_{11}}\,{\it t_{12}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{11}}\,{\it t_{12}}\right)^{{{4}\over{3}}}-2.0\right)\,\left(-{\it t_7}\,{\it t_9}+{\it t_7}\,\left(4.3424534362958 \times 10^{-4}\,{\it t_2}\,{\it t_8}+0.01555\,{\it t_8}-0.00297768235631712\,{\it t_2}-0.0269\right)-{{0.0843\,{\it t_6}}\over{1.101176160755631\,{\it t_3}+0.1619735131738333\,{\it t_2}+1.0}}+0.1423\,{\it t_4}\,{\it t_6}\right)+{\it t_7}\,{\it t_9}-0.1423\,{\it t_4}\,{\it t_6}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{13}}\,\left(1.923661050931536\,\left(\left({{\rho_s+{\it t_{10}}}\over{{\it t_{13}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left(-{\it t_{19}}\,{\it t_{21}}+{\it t_{19}}\,\left(4.3424534362958 \times 10^{-4}\,{\it t_{14}}\,{\it t_{20}}+0.01555\,{\it t_{20}}-0.00297768235631712\,{\it t_{14}}-0.0269\right)-{{0.0843\,{\it t_{18}}}\over{1.101176160755631\,{\it t_{15}}+0.1619735131738333\,{\it t_{14}}+1.0}}+0.1423\,{\it t_{16}}\,{\it t_{18}}\right)+{\it t_{19}}\,{\it t_{21}}-0.1423\,{\it t_{16}}\,{\it t_{18}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, A Zunger, Phys.Rev. B 23, 5048 (1981)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.23.5048 "> 
C> 10.1103/PhysRevB.23.5048 </a> 
C> 
      subroutine nwxcm_c_perdew81_d2(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
C> \brief Evaluate the nwxcm_c_perdew81 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_1}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.8292885914166397\,{\it t_3}+0.2068248536658599\,{\it t_2}+1.0}}\\\\ 
C>   {\it t_5} &=& 0.6203504908994\,{\it t_2}\\\\ 
C>   {\it t_6} &=& {\it nwxcm\_heaviside}\left({\it t_5}-1.0\right)\\\\ 
C>   {\it t_7} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994\,{\it t_2}\right)\\\\ 
C>   {\it t_8} &=& \log {\it t_5}\\\\ 
C>   {\it t_9} &=& 0.0012407009817988\,{\it t_2}\,{\it t_8}+0.0311\,{\it t_8}-0.007196065694433039\,{\it t_2}-0.048\\\\ 
C>   {\it t_{10}} &=& -\rho_\beta\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}+\rho_\alpha\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{13}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_{13}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{13}}^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{0.8292885914166397\,{\it t_{15}}+0.2068248536658599\,{\it t_{14}}+1.0}}\\\\ 
C>   {\it t_{17}} &=& 0.6203504908994\,{\it t_{14}}\\\\ 
C>   {\it t_{18}} &=& {\it nwxcm\_heaviside}\left({\it t_{17}}-1.0\right)\\\\ 
C>   {\it t_{19}} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994\,{\it t_{14}}\right)\\\\ 
C>   {\it t_{20}} &=& \log {\it t_{17}}\\\\ 
C>   {\it t_{21}} &=& 0.0012407009817988\,{\it t_{14}}\,{\it t_{20}}+0.0311\,{\it t_{20}}-0.007196065694433039\,{\it t_{14}}-0.048\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({\it t_{11}}\,{\it t_{12}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{11}}\,{\it t_{12}}\right)^{{{4}\over{3}}}-2.0\right)\,\left(-{\it t_7}\,{\it t_9}+{\it t_7}\,\left(4.3424534362958 \times 10^{-4}\,{\it t_2}\,{\it t_8}+0.01555\,{\it t_8}-0.00297768235631712\,{\it t_2}-0.0269\right)-{{0.0843\,{\it t_6}}\over{1.101176160755631\,{\it t_3}+0.1619735131738333\,{\it t_2}+1.0}}+0.1423\,{\it t_4}\,{\it t_6}\right)+{\it t_7}\,{\it t_9}-0.1423\,{\it t_4}\,{\it t_6}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{13}}\,\left(1.923661050931536\,\left(\left({{\rho_s+{\it t_{10}}}\over{{\it t_{13}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left(-{\it t_{19}}\,{\it t_{21}}+{\it t_{19}}\,\left(4.3424534362958 \times 10^{-4}\,{\it t_{14}}\,{\it t_{20}}+0.01555\,{\it t_{20}}-0.00297768235631712\,{\it t_{14}}-0.0269\right)-{{0.0843\,{\it t_{18}}}\over{1.101176160755631\,{\it t_{15}}+0.1619735131738333\,{\it t_{14}}+1.0}}+0.1423\,{\it t_{16}}\,{\it t_{18}}\right)+{\it t_{19}}\,{\it t_{21}}-0.1423\,{\it t_{16}}\,{\it t_{18}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, A Zunger, Phys.Rev. B 23, 5048 (1981)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.23.5048 "> 
C> 10.1103/PhysRevB.23.5048 </a> 
C> 
      subroutine nwxcm_c_perdew81_d3(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
