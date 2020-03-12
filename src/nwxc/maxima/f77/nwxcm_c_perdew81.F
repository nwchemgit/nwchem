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
C>   {\it t_4} &=& {{1}\over{0.8292885914166397\,{\it t_3}
C>    +0.2068248536658599\,{\it t_2}+1.0}}\\\\ 
C>   {\it t_5} &=& 0.6203504908994\,{\it t_2}\\\\ 
C>   {\it t_6} &=& {\it nwxcm\_heaviside}\left({\it t_5}-1.0\right)\\\\ 
C>   {\it t_7} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994
C>    \,{\it t_2}\right)\\\\ 
C>   {\it t_8} &=& \log {\it t_5}\\\\ 
C>   {\it t_9} &=& 0.0012407009817988\,{\it t_2}\,{\it t_8}+0.0311
C>    \,{\it t_8}-0.007196065694433039\,{\it t_2}-0.048\\\\ 
C>   {\it t_{10}} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& 0.6203504908994\,{\it t_{12}}\\\\ 
C>   {\it t_{14}} &=& {\it nwxcm\_heaviside}\left({\it t_{13}}
C>    -1.0\right)\\\\ 
C>   {\it t_{15}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{0.8292885914166397\,{\it t_{15}}
C>    +0.2068248536658599\,{\it t_{12}}+1.0}}\\\\ 
C>   {\it t_{17}} &=& {\it nwxcm\_heaviside}\left(1.0
C>    -0.6203504908994\,{\it t_{12}}\right)\\\\ 
C>   {\it t_{18}} &=& \log {\it t_{13}}\\\\ 
C>   {\it t_{19}} &=& 0.0012407009817988\,{\it t_{18}}\,{\it t_{12}}
C>    -0.007196065694433039\,{\it t_{12}}+0.0311\,{\it t_{18}}
C>    -0.048\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({
C>    \it t_{10}}\,{\it t_{11}}+1.0\right)^{{{4}\over{3}}}+\left(1.0
C>    -{\it t_{10}}\,{\it t_{11}}\right)^{{{4}\over{3}}}-2.0\right)
C>    \,\left(-{\it t_7}\,{\it t_9}+{\it t_7}\,
C>    \left(4.3424534362958 \times 10^{-4}\,{\it t_2}\,{\it t_8}
C>    +0.01555\,{\it t_8}-0.00297768235631712\,{\it t_2}-0.0269\right)
C>    -{{0.0843\,{\it t_6}}\over{1.101176160755631\,{\it t_3}
C>    +0.1619735131738333\,{\it t_2}+1.0}}+0.1423\,{\it t_4}\,{
C>    \it t_6}\right)+{\it t_7}\,{\it t_9}-0.1423\,{\it t_4}\,{
C>    \it t_6}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(0.9999999999999991\,\left(-{\it t_{17}}\,{
C>    \it t_{19}}+{\it t_{17}}\,\left(4.3424534362958 \times 10^{-4}
C>    \,{\it t_{18}}\,{\it t_{12}}-0.00297768235631712\,{\it t_{12}}
C>    +0.01555\,{\it t_{18}}-0.0269\right)-{{0.0843\,{
C>    \it t_{14}}}\over{1.101176160755631\,{\it t_{15}}
C>    +0.1619735131738333\,{\it t_{12}}+1.0}}+0.1423\,{\it t_{14}}
C>    \,{\it t_{16}}\right)+{\it t_{17}}\,{\it t_{19}}-0.1423\,{
C>    \it t_{14}}\,{\it t_{16}}\right)\,\rho_s\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, A Zunger, Phys.Rev. B 23, 5048 (1981)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.23.5048 "> 
C> 10.1103/PhysRevB.23.5048 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_perdew81(param,tol_rho,ipol,nq,wght, 
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
            t9 = 3.968502629920499d-1
            t10 = 1/rhoa**1.3333333333333333d+0
            fnc(iq) = 2.0d+0*rhoa*t8*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*rhoa*(t6*(-4.1356699393293334d-4*t10
     1         *t7*t9+1.9851215708780795d-3*t10*t9-5.183333333333333d-3/
     2         rhoa)+1.423d-1*t4*(-6.894161788861997d-2*t10*t9-6.1567678
     3         58845764d-2/rhoa**1.1666666666666667d+0)/t5**2)*wght+1.0d
     4         +0*t8*wght+Amat(iq,D1_RA)
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
            t1 = 1/rhoa**3.333333333333333d-1
            t2 = 6.203504908994d-1*t1
            t3 = nwxcm_heaviside(t2-1.0d+0)
            t4 = 1/rhoa**1.6666666666666666d-1
            t5 = 8.292885914166397d-1*t4+2.0682485366585993d-1*t1+1.0d+0
            t6 = 1/t5
            t7 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1)
            t8 = log(t2)
            t9 = 1.2407009817988d-3*t1*t8+3.11d-2*t8-7.196065694433039d-
     1         3*t1-4.8d-2
            t10 = 1.101176160755631d+0*t4+1.6197351317383332d-1*t1+1.0d+
     1         0
            t11 = 9.999999999999991d-1*(-t7*t9+t7*(4.3424534362957995d-4
     1         *t1*t8+1.555d-2*t8-2.9776823563171195d-3*t1-2.69d-2)+1.42
     2         3d-1*t3*t6-8.43d-2*t3/t10)+t7*t9-1.423d-1*t3*t6
            t12 = 1/rhoa**1.3333333333333333d+0
            t13 = 1/rhoa**1.1666666666666667d+0
            t14 = -1.382147652361066d-1*t13-6.894161788861997d-2*t12
            t15 = 1/t5**2
            t16 = 1/rhoa
            t17 = -4.1356699393293334d-4*t12*t8-1.0366666666666666d-2*t1
     1         6+1.9851215708780795d-3*t12
            fnc(iq) = 1.0d+0*rhoa*t11*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*rhoa*(9.999999999999991d-1*(t7*(-1.4
     1         474844787652666d-4*t12*t8-5.183333333333333d-3*t16+8.4781
     2         23375625131d-4*t12)-t17*t7-1.423d-1*t14*t15*t3+8.43d-2*(-
     3         1.8352936012593848d-1*t13-5.399117105794444d-2*t12)*t3/t1
     4         0**2)+t17*t7+1.423d-1*t14*t15*t3)*wght+1.0d+0*t11*wght+Am
     5         at(iq,D1_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**3.333333333333333d-1
            t2 = 6.203504908994d-1*t1
            t3 = nwxcm_heaviside(t2-1.0d+0)
            t4 = 1/rhob**1.6666666666666666d-1
            t5 = 8.292885914166397d-1*t4+2.0682485366585993d-1*t1+1.0d+0
            t6 = 1/t5
            t7 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1)
            t8 = log(t2)
            t9 = 1.2407009817988d-3*t1*t8+3.11d-2*t8-7.196065694433039d-
     1         3*t1-4.8d-2
            t10 = 1.101176160755631d+0*t4+1.6197351317383332d-1*t1+1.0d+
     1         0
            t11 = 9.999999999999991d-1*(-t7*t9+t7*(4.3424534362957995d-4
     1         *t1*t8+1.555d-2*t8-2.9776823563171195d-3*t1-2.69d-2)+1.42
     2         3d-1*t3*t6-8.43d-2*t3/t10)+t7*t9-1.423d-1*t3*t6
            t12 = 1/rhob**1.3333333333333333d+0
            t13 = 1/rhob**1.1666666666666667d+0
            t14 = -1.382147652361066d-1*t13-6.894161788861997d-2*t12
            t15 = 1/t5**2
            t16 = 1/rhob
            t17 = -4.1356699393293334d-4*t12*t8-1.0366666666666666d-2*t1
     1         6+1.9851215708780795d-3*t12
            fnc(iq) = 1.0d+0*rhob*t11*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*rhob*(9.999999999999991d-1*(t7*(-1.4
     1         474844787652666d-4*t12*t8-5.183333333333333d-3*t16+8.4781
     2         23375625131d-4*t12)-t17*t7-1.423d-1*t14*t15*t3+8.43d-2*(-
     3         1.8352936012593848d-1*t13-5.399117105794444d-2*t12)*t3/t1
     4         0**2)+t17*t7+1.423d-1*t14*t15*t3)*wght+1.0d+0*t11*wght+Am
     5         at(iq,D1_RB)
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
C>   {\it t_4} &=& {{1}\over{0.8292885914166397\,{\it t_3}
C>    +0.2068248536658599\,{\it t_2}+1.0}}\\\\ 
C>   {\it t_5} &=& 0.6203504908994\,{\it t_2}\\\\ 
C>   {\it t_6} &=& {\it nwxcm\_heaviside}\left({\it t_5}-1.0\right)\\\\ 
C>   {\it t_7} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994
C>    \,{\it t_2}\right)\\\\ 
C>   {\it t_8} &=& \log {\it t_5}\\\\ 
C>   {\it t_9} &=& 0.0012407009817988\,{\it t_2}\,{\it t_8}+0.0311
C>    \,{\it t_8}-0.007196065694433039\,{\it t_2}-0.048\\\\ 
C>   {\it t_{10}} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& 0.6203504908994\,{\it t_{12}}\\\\ 
C>   {\it t_{14}} &=& {\it nwxcm\_heaviside}\left({\it t_{13}}
C>    -1.0\right)\\\\ 
C>   {\it t_{15}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{0.8292885914166397\,{\it t_{15}}
C>    +0.2068248536658599\,{\it t_{12}}+1.0}}\\\\ 
C>   {\it t_{17}} &=& {\it nwxcm\_heaviside}\left(1.0
C>    -0.6203504908994\,{\it t_{12}}\right)\\\\ 
C>   {\it t_{18}} &=& \log {\it t_{13}}\\\\ 
C>   {\it t_{19}} &=& 0.0012407009817988\,{\it t_{18}}\,{\it t_{12}}
C>    -0.007196065694433039\,{\it t_{12}}+0.0311\,{\it t_{18}}
C>    -0.048\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({
C>    \it t_{10}}\,{\it t_{11}}+1.0\right)^{{{4}\over{3}}}+\left(1.0
C>    -{\it t_{10}}\,{\it t_{11}}\right)^{{{4}\over{3}}}-2.0\right)
C>    \,\left(-{\it t_7}\,{\it t_9}+{\it t_7}\,
C>    \left(4.3424534362958 \times 10^{-4}\,{\it t_2}\,{\it t_8}
C>    +0.01555\,{\it t_8}-0.00297768235631712\,{\it t_2}-0.0269\right)
C>    -{{0.0843\,{\it t_6}}\over{1.101176160755631\,{\it t_3}
C>    +0.1619735131738333\,{\it t_2}+1.0}}+0.1423\,{\it t_4}\,{
C>    \it t_6}\right)+{\it t_7}\,{\it t_9}-0.1423\,{\it t_4}\,{
C>    \it t_6}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(0.9999999999999991\,\left(-{\it t_{17}}\,{
C>    \it t_{19}}+{\it t_{17}}\,\left(4.3424534362958 \times 10^{-4}
C>    \,{\it t_{18}}\,{\it t_{12}}-0.00297768235631712\,{\it t_{12}}
C>    +0.01555\,{\it t_{18}}-0.0269\right)-{{0.0843\,{
C>    \it t_{14}}}\over{1.101176160755631\,{\it t_{15}}
C>    +0.1619735131738333\,{\it t_{12}}+1.0}}+0.1423\,{\it t_{14}}
C>    \,{\it t_{16}}\right)+{\it t_{17}}\,{\it t_{19}}-0.1423\,{
C>    \it t_{14}}\,{\it t_{16}}\right)\,\rho_s\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, A Zunger, Phys.Rev. B 23, 5048 (1981)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.23.5048 "> 
C> 10.1103/PhysRevB.23.5048 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_perdew81_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = 7.937005259840998d-1
            t2 = 1/rhoa**3.333333333333333d-1
            t3 = 6.203504908994d-1*t1*t2
            t4 = nwxcm_heaviside(t3-1.0d+0)
            t5 = 8.908987181403393d-1
            t6 = 1/rhoa**1.6666666666666666d-1
            t7 = 8.292885914166397d-1*t5*t6+2.0682485366585993d-1*t1*t2+
     1         1.0d+0
            t8 = 1/t7
            t9 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1*t2)
            t10 = log(t3)
            t11 = 1.2407009817988d-3*t1*t10*t2-7.196065694433039d-3*t1*t
     1         2+3.11d-2*t10-4.8d-2
            t12 = t11*t9-1.423d-1*t4*t8
            t13 = 3.968502629920499d-1
            t14 = 1/rhoa**1.3333333333333333d+0
            t15 = -6.894161788861997d-2*t13*t14-6.156767858845764d-2/rho
     1         a**1.1666666666666667d+0
            t16 = 1/t7**2
            t17 = (-4.1356699393293334d-4*t10*t13*t14+1.9851215708780795
     1         d-3*t13*t14-5.183333333333333d-3/rhoa)*t9+1.423d-1*t15*t1
     2         6*t4
            t18 = 2.0d+0*t17*wght
            t19 = -2.846d-1*t15**2*t4/t7**3
            t20 = 1.9842513149602492d-1
            t21 = 1/rhoa**2.3333333333333334d+0
            t22 = 1.423d-1*t16*(9.192215718482661d-2*t20*t21+3.591447917
     1         660031d-2/rhoa**2.1666666666666665d+0)*t4
            t23 = 1/rhoa**2
            t24 = (2.5916666666666666d-3*t23+5.514226585772444d-4*t10*t2
     1         0*t21-2.5089730965264617d-3*t20*t21)*t9
            t25 = (4.3424534362957995d-4*t1*t10*t2-2.9776823563171195d-3
     1         *t1*t2+1.555d-2*t10-2.69d-2)*t9-t11*t9+1.423d-1*t4*t8-8.4
     2         3d-2*t4/(1.101176160755631d+0*t5*t6+1.6197351317383332d-1
     3         *t1*t2+1.0d+0)
            fnc(iq) = 2.0d+0*rhoa*t12*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*rhoa*t17*wght+1.0d+0*t12*wght+Amat(i
     1         q,D1_RA)
            Amat2(iq,D2_RA_RA) = 2.0d+0*rhoa*(4.2748023354034137d-1*t23*
     1         t25+t24+t22+t19)*wght+t18+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 2.0d+0*rhoa*(-4.2748023354034137d-1*t23
     1         *t25+t24+t22+t19)*wght+t18+Amat2(iq,D2_RA_RB)
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
            t26 = -1.8352936012593848d-1*t20-5.399117105794444d-2*t19
            t27 = 1/t16**2
            t28 = t8*(-1.4474844787652666d-4*t19*t9+8.478123375625131d-4
     1         *t19-5.183333333333333d-3*t12)-t24*t8+8.43d-2*t26*t27*t7-
     2         1.423d-1*t21*t22*t7
            t29 = 1.923661050931536d+0*t15*t28
            t30 = 1/t1**2
            t31 = t11*t30
            t32 = -t12
            t33 = t32+t31
            t34 = t13**3.333333333333333d-1
            t35 = -t11*t30
            t36 = t35+t12
            t37 = t14**3.333333333333333d-1
            t38 = 1.3333333333333333d+0*t36*t37+1.3333333333333333d+0*t3
     1         3*t34
            t39 = 1.923661050931536d+0*t17*t38+t29+t25+t23
            t40 = 1.0d+0*t18*wght
            t41 = t31+t12
            t42 = t35+t32
            t43 = 1.3333333333333333d+0*t37*t42+1.3333333333333333d+0*t3
     1         4*t41
            t44 = 1.923661050931536d+0*t17*t43+t29+t25+t23
            t45 = t21**2
            t46 = 1/t4**3
            t47 = -2.846d-1*t45*t46*t7
            t48 = 1/t1**2.3333333333333334d+0
            t49 = 1/t1**2.1666666666666665d+0
            t50 = 1.6125055944212438d-1*t49+9.192215718482661d-2*t48
            t51 = 1.423d-1*t22*t50*t7
            t52 = 5.514226585772444d-4*t48*t9-2.5089730965264617d-3*t48+
     1         1.0366666666666666d-2*t30
            t53 = t52*t8
            t54 = 1.923661050931536d+0*t15*(t8*(1.9299793050203554d-4*t4
     1         8*t9-1.082166967457842d-3*t48+5.183333333333333d-3*t30)-t
     2         52*t8-1.423d-1*t22*t50*t7+8.43d-2*t27*(2.141175868135949d
     3         -1*t49+7.198822807725924d-2*t48)*t7+2.846d-1*t45*t46*t7-1
     4         .686d-1*t26**2*t7/t16**3)
            t55 = 1/t13**6.666666666666666d-1
            t56 = 1/t1**3
            t57 = -2*t11*t56
            t58 = 2*t30
            t59 = 1/t14**6.666666666666666d-1
            t60 = 2*t11*t56
            t61 = -2*t30
            fnc(iq) = 1.0d+0*t1*t18*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*t39*wght+t40+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*t1*t44*wght+t40+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*t1*(1.923661050931536d+0*t17*(1.
     1         3333333333333333d+0*t37*(t61+t60)+4.444444444444444d-1*t3
     2         6**2*t59+1.3333333333333333d+0*t34*(t58+t57)+4.4444444444
     3         44444d-1*t33**2*t55)+t54+t53+t51+t47+3.847322101863072d+0
     4         *t28*t38)*wght+2.0d+0*t39*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*t1*(1.923661050931536d+0*t17*(4.
     1         444444444444444d-1*t36*t42*t59+2.6666666666666666d+0*t11*
     2         t37*t56-2.6666666666666666d+0*t11*t34*t56+4.4444444444444
     3         44d-1*t33*t41*t55)+t54+t53+t51+t47+1.923661050931536d+0*t
     4         28*t43+1.923661050931536d+0*t28*t38)*wght+1.0d+0*t44*wght
     5         +1.0d+0*t39*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*t1*(1.923661050931536d+0*t17*(1.
     1         3333333333333333d+0*t34*(t61+t57)+1.3333333333333333d+0*t
     2         37*(t60+t58)+4.444444444444444d-1*t42**2*t59+4.4444444444
     3         44444d-1*t41**2*t55)+t54+t53+t51+t47+3.847322101863072d+0
     4         *t28*t43)*wght+2.0d+0*t44*wght+Amat2(iq,D2_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**3.333333333333333d-1
            t2 = 6.203504908994d-1*t1
            t3 = nwxcm_heaviside(t2-1.0d+0)
            t4 = 1/rhoa**1.6666666666666666d-1
            t5 = 8.292885914166397d-1*t4+2.0682485366585993d-1*t1+1.0d+0
            t6 = 1/t5
            t7 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1)
            t8 = log(t2)
            t9 = 1.2407009817988d-3*t1*t8+3.11d-2*t8-7.196065694433039d-
     1         3*t1-4.8d-2
            t10 = 1.101176160755631d+0*t4+1.6197351317383332d-1*t1+1.0d+
     1         0
            t11 = 9.999999999999991d-1*(-t7*t9+t7*(4.3424534362957995d-4
     1         *t1*t8+1.555d-2*t8-2.9776823563171195d-3*t1-2.69d-2)+1.42
     2         3d-1*t3*t6-8.43d-2*t3/t10)+t7*t9-1.423d-1*t3*t6
            t12 = 1/rhoa**1.3333333333333333d+0
            t13 = 1/rhoa**1.1666666666666667d+0
            t14 = -1.382147652361066d-1*t13-6.894161788861997d-2*t12
            t15 = 1/t5**2
            t16 = -1.8352936012593848d-1*t13-5.399117105794444d-2*t12
            t17 = 1/t10**2
            t18 = 1/rhoa
            t19 = -4.1356699393293334d-4*t12*t8-1.0366666666666666d-2*t1
     1         8+1.9851215708780795d-3*t12
            t20 = 9.999999999999991d-1*(t7*(-1.4474844787652666d-4*t12*t
     1         8-5.183333333333333d-3*t18+8.478123375625131d-4*t12)-t19*
     2         t7+8.43d-2*t16*t17*t3-1.423d-1*t14*t15*t3)+t19*t7+1.423d-
     3         1*t14*t15*t3
            t21 = t14**2
            t22 = 1/t5**3
            t23 = 1/rhoa**2.3333333333333334d+0
            t24 = 1/rhoa**2.1666666666666665d+0
            t25 = 1.6125055944212438d-1*t24+9.192215718482661d-2*t23
            t26 = 1/rhoa**2
            t27 = 5.514226585772444d-4*t23*t8+1.0366666666666666d-2*t26-
     1         2.5089730965264617d-3*t23
            fnc(iq) = 1.0d+0*rhoa*t11*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*rhoa*t20*wght+1.0d+0*t11*wght+Amat(i
     1         q,D1_RA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*rhoa*(9.999999999999991d-1*(t7*(
     1         1.9299793050203554d-4*t23*t8+5.183333333333333d-3*t26-1.0
     2         82166967457842d-3*t23)-t27*t7-1.423d-1*t15*t25*t3+8.43d-2
     3         *t17*(2.141175868135949d-1*t24+7.198822807725924d-2*t23)*
     4         t3+2.846d-1*t21*t22*t3-1.686d-1*t16**2*t3/t10**3)+t27*t7+
     5         1.423d-1*t15*t25*t3-2.846d-1*t21*t22*t3)*wght+2.0d+0*t20*
     6         wght+Amat2(iq,D2_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**3.333333333333333d-1
            t2 = 6.203504908994d-1*t1
            t3 = nwxcm_heaviside(t2-1.0d+0)
            t4 = 1/rhob**1.6666666666666666d-1
            t5 = 8.292885914166397d-1*t4+2.0682485366585993d-1*t1+1.0d+0
            t6 = 1/t5
            t7 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1)
            t8 = log(t2)
            t9 = 1.2407009817988d-3*t1*t8+3.11d-2*t8-7.196065694433039d-
     1         3*t1-4.8d-2
            t10 = 1.101176160755631d+0*t4+1.6197351317383332d-1*t1+1.0d+
     1         0
            t11 = 9.999999999999991d-1*(-t7*t9+t7*(4.3424534362957995d-4
     1         *t1*t8+1.555d-2*t8-2.9776823563171195d-3*t1-2.69d-2)+1.42
     2         3d-1*t3*t6-8.43d-2*t3/t10)+t7*t9-1.423d-1*t3*t6
            t12 = 1/rhob**1.3333333333333333d+0
            t13 = 1/rhob**1.1666666666666667d+0
            t14 = -1.382147652361066d-1*t13-6.894161788861997d-2*t12
            t15 = 1/t5**2
            t16 = -1.8352936012593848d-1*t13-5.399117105794444d-2*t12
            t17 = 1/t10**2
            t18 = 1/rhob
            t19 = -4.1356699393293334d-4*t12*t8-1.0366666666666666d-2*t1
     1         8+1.9851215708780795d-3*t12
            t20 = 9.999999999999991d-1*(t7*(-1.4474844787652666d-4*t12*t
     1         8-5.183333333333333d-3*t18+8.478123375625131d-4*t12)-t19*
     2         t7+8.43d-2*t16*t17*t3-1.423d-1*t14*t15*t3)+t19*t7+1.423d-
     3         1*t14*t15*t3
            t21 = t14**2
            t22 = 1/t5**3
            t23 = 1/rhob**2.3333333333333334d+0
            t24 = 1/rhob**2.1666666666666665d+0
            t25 = 1.6125055944212438d-1*t24+9.192215718482661d-2*t23
            t26 = 1/rhob**2
            t27 = 5.514226585772444d-4*t23*t8+1.0366666666666666d-2*t26-
     1         2.5089730965264617d-3*t23
            fnc(iq) = 1.0d+0*rhob*t11*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*rhob*t20*wght+1.0d+0*t11*wght+Amat(i
     1         q,D1_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*rhob*(9.999999999999991d-1*(t7*(
     1         1.9299793050203554d-4*t23*t8+5.183333333333333d-3*t26-1.0
     2         82166967457842d-3*t23)-t27*t7-1.423d-1*t15*t25*t3+8.43d-2
     3         *t17*(2.141175868135949d-1*t24+7.198822807725924d-2*t23)*
     4         t3+2.846d-1*t21*t22*t3-1.686d-1*t16**2*t3/t10**3)+t27*t7+
     5         1.423d-1*t15*t25*t3-2.846d-1*t21*t22*t3)*wght+2.0d+0*t20*
     6         wght+Amat2(iq,D2_RB_RB)
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
C>   {\it t_4} &=& {{1}\over{0.8292885914166397\,{\it t_3}
C>    +0.2068248536658599\,{\it t_2}+1.0}}\\\\ 
C>   {\it t_5} &=& 0.6203504908994\,{\it t_2}\\\\ 
C>   {\it t_6} &=& {\it nwxcm\_heaviside}\left({\it t_5}-1.0\right)\\\\ 
C>   {\it t_7} &=& {\it nwxcm\_heaviside}\left(1.0-0.6203504908994
C>    \,{\it t_2}\right)\\\\ 
C>   {\it t_8} &=& \log {\it t_5}\\\\ 
C>   {\it t_9} &=& 0.0012407009817988\,{\it t_2}\,{\it t_8}+0.0311
C>    \,{\it t_8}-0.007196065694433039\,{\it t_2}-0.048\\\\ 
C>   {\it t_{10}} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_s^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& 0.6203504908994\,{\it t_{12}}\\\\ 
C>   {\it t_{14}} &=& {\it nwxcm\_heaviside}\left({\it t_{13}}
C>    -1.0\right)\\\\ 
C>   {\it t_{15}} &=& {{1}\over{\rho_s^{{{1}\over{6}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{0.8292885914166397\,{\it t_{15}}
C>    +0.2068248536658599\,{\it t_{12}}+1.0}}\\\\ 
C>   {\it t_{17}} &=& {\it nwxcm\_heaviside}\left(1.0
C>    -0.6203504908994\,{\it t_{12}}\right)\\\\ 
C>   {\it t_{18}} &=& \log {\it t_{13}}\\\\ 
C>   {\it t_{19}} &=& 0.0012407009817988\,{\it t_{18}}\,{\it t_{12}}
C>    -0.007196065694433039\,{\it t_{12}}+0.0311\,{\it t_{18}}
C>    -0.048\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({
C>    \it t_{10}}\,{\it t_{11}}+1.0\right)^{{{4}\over{3}}}+\left(1.0
C>    -{\it t_{10}}\,{\it t_{11}}\right)^{{{4}\over{3}}}-2.0\right)
C>    \,\left(-{\it t_7}\,{\it t_9}+{\it t_7}\,
C>    \left(4.3424534362958 \times 10^{-4}\,{\it t_2}\,{\it t_8}
C>    +0.01555\,{\it t_8}-0.00297768235631712\,{\it t_2}-0.0269\right)
C>    -{{0.0843\,{\it t_6}}\over{1.101176160755631\,{\it t_3}
C>    +0.1619735131738333\,{\it t_2}+1.0}}+0.1423\,{\it t_4}\,{
C>    \it t_6}\right)+{\it t_7}\,{\it t_9}-0.1423\,{\it t_4}\,{
C>    \it t_6}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(0.9999999999999991\,\left(-{\it t_{17}}\,{
C>    \it t_{19}}+{\it t_{17}}\,\left(4.3424534362958 \times 10^{-4}
C>    \,{\it t_{18}}\,{\it t_{12}}-0.00297768235631712\,{\it t_{12}}
C>    +0.01555\,{\it t_{18}}-0.0269\right)-{{0.0843\,{
C>    \it t_{14}}}\over{1.101176160755631\,{\it t_{15}}
C>    +0.1619735131738333\,{\it t_{12}}+1.0}}+0.1423\,{\it t_{14}}
C>    \,{\it t_{16}}\right)+{\it t_{17}}\,{\it t_{19}}-0.1423\,{
C>    \it t_{14}}\,{\it t_{16}}\right)\,\rho_s\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, A Zunger, Phys.Rev. B 23, 5048 (1981)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.23.5048 "> 
C> 10.1103/PhysRevB.23.5048 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_perdew81_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = 7.937005259840998d-1
            t2 = 1/rhoa**3.333333333333333d-1
            t3 = 6.203504908994d-1*t1*t2
            t4 = nwxcm_heaviside(t3-1.0d+0)
            t5 = 8.908987181403393d-1
            t6 = 1/rhoa**1.6666666666666666d-1
            t7 = 8.292885914166397d-1*t5*t6+2.0682485366585993d-1*t1*t2+
     1         1.0d+0
            t8 = 1/t7
            t9 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1*t2)
            t10 = log(t3)
            t11 = 1.2407009817988d-3*t1*t10*t2-7.196065694433039d-3*t1*t
     1         2+3.11d-2*t10-4.8d-2
            t12 = t11*t9-1.423d-1*t4*t8
            t13 = 3.968502629920499d-1
            t14 = 1/rhoa**1.3333333333333333d+0
            t15 = 4.4544935907016964d-1
            t16 = 1/rhoa**1.1666666666666667d+0
            t17 = -1.382147652361066d-1*t15*t16-6.894161788861997d-2*t13
     1         *t14
            t18 = 1/t7**2
            t19 = 1/rhoa
            t20 = -5.183333333333333d-3*t19-4.1356699393293334d-4*t10*t1
     1         3*t14+1.9851215708780795d-3*t13*t14
            t21 = t20*t9+1.423d-1*t17*t18*t4
            t22 = 2.0d+0*t21*wght
            t23 = 1/t7**3
            t24 = -2.846d-1*t17**2*t23*t4
            t25 = 1.9842513149602492d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 9.192215718482661d-2*t25*t26+3.591447917660031d-2/rhoa
     1         **2.1666666666666665d+0
            t28 = 1.423d-1*t18*t27*t4
            t29 = 1/rhoa**2
            t30 = (2.5916666666666666d-3*t29+5.514226585772444d-4*t10*t2
     1         5*t26-2.5089730965264617d-3*t25*t26)*t9
            t31 = 1.101176160755631d+0*t5*t6+1.6197351317383332d-1*t1*t2
     1         +1.0d+0
            t32 = (4.3424534362957995d-4*t1*t10*t2-2.9776823563171195d-3
     1         *t1*t2+1.555d-2*t10-2.69d-2)*t9-t11*t9+1.423d-1*t4*t8-8.4
     2         3d-2*t4/t31
            t33 = 4.2748023354034137d-1*t29*t32+t30+t28+t24
            t34 = -4.2748023354034137d-1*t29*t32+t30+t28+t24
            t35 = 8.538000000000001d-1*t17**3*t4/t7**4
            t36 = -8.538000000000001d-1*t17*t23*t27*t4
            t37 = 9.921256574801247d-2
            t38 = 1/rhoa**3.3333333333333337d+0
            t39 = 1.423d-1*t18*(-2.1448503343126213d-1*t37*t38-3.8907352
     1         441317d-2/rhoa**3.1666666666666665d+0)*t4
            t40 = 1/rhoa**3
            t41 = (-2.5916666666666666d-3*t40-1.2866528700135701d-3*t10*
     1         t37*t38+5.670463005702662d-3*t37*t38)*t9
            t42 = -t20*t9+(-2.5916666666666666d-3*t19-1.4474844787652666
     1         d-4*t10*t13*t14+8.478123375625131d-4*t13*t14)*t9+8.43d-2*
     2         (-1.8352936012593848d-1*t15*t16-5.399117105794444d-2*t13*
     3         t14)*t4/t31**2-1.423d-1*t17*t18*t4
            fnc(iq) = 2.0d+0*rhoa*t12*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*rhoa*t21*wght+1.0d+0*t12*wght+Amat(i
     1         q,D1_RA)
            Amat2(iq,D2_RA_RA) = 2.0d+0*rhoa*t33*wght+t22+Amat2(iq,D2_RA
     1         _RA)
            Amat2(iq,D2_RA_RB) = 2.0d+0*rhoa*t34*wght+t22+Amat2(iq,D2_RA
     1         _RB)
            Amat3(iq,D3_RA_RA_RA) = 2.0d+0*rhoa*(1.282440700621024d+0*t2
     1         9*t42+t41-1.282440700621024d+0*t32*t40+t39+t36+t35)*wght+
     2         3.0d+0*t33*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 2.0d+0*rhoa*(-4.2748023354034137d-1*
     1         t29*t42+t41+4.2748023354034137d-1*t32*t40+t39+t36+t35)*wg
     2         ht+2.0d+0*t34*wght+1.0d+0*t33*wght+Amat3(iq,D3_RA_RA_RB)
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
            t26 = -1.8352936012593848d-1*t20-5.399117105794444d-2*t19
            t27 = 1/t16**2
            t28 = t8*(-1.4474844787652666d-4*t19*t9+8.478123375625131d-4
     1         *t19-5.183333333333333d-3*t12)-t24*t8+8.43d-2*t26*t27*t7-
     2         1.423d-1*t21*t22*t7
            t29 = 1.923661050931536d+0*t15*t28
            t30 = 1/t1**2
            t31 = t11*t30
            t32 = -t12
            t33 = t32+t31
            t34 = t13**3.333333333333333d-1
            t35 = -t11*t30
            t36 = t35+t12
            t37 = t14**3.333333333333333d-1
            t38 = 1.3333333333333333d+0*t36*t37+1.3333333333333333d+0*t3
     1         3*t34
            t39 = 1.923661050931536d+0*t17*t38+t29+t25+t23
            t40 = 1.0d+0*t18*wght
            t41 = t31+t12
            t42 = t35+t32
            t43 = 1.3333333333333333d+0*t37*t42+1.3333333333333333d+0*t3
     1         4*t41
            t44 = 1.923661050931536d+0*t17*t43+t29+t25+t23
            t45 = t21**2
            t46 = 1/t4**3
            t47 = -2.846d-1*t45*t46*t7
            t48 = 1/t1**2.3333333333333334d+0
            t49 = 1/t1**2.1666666666666665d+0
            t50 = 1.6125055944212438d-1*t49+9.192215718482661d-2*t48
            t51 = 1.423d-1*t22*t50*t7
            t52 = 5.514226585772444d-4*t48*t9-2.5089730965264617d-3*t48+
     1         1.0366666666666666d-2*t30
            t53 = t52*t8
            t54 = 1/t16**3
            t55 = 2.141175868135949d-1*t49+7.198822807725924d-2*t48
            t56 = t8*(1.9299793050203554d-4*t48*t9-1.082166967457842d-3*
     1         t48+5.183333333333333d-3*t30)-t52*t8+8.43d-2*t27*t55*t7-1
     2         .686d-1*t26**2*t54*t7-1.423d-1*t22*t50*t7+2.846d-1*t45*t4
     3         6*t7
            t57 = 1.923661050931536d+0*t15*t56
            t58 = t33**2
            t59 = 1/t13**6.666666666666666d-1
            t60 = 1/t1**3
            t61 = -2*t11*t60
            t62 = 2*t30
            t63 = t62+t61
            t64 = t36**2
            t65 = 1/t14**6.666666666666666d-1
            t66 = 2*t11*t60
            t67 = -2*t30
            t68 = t67+t66
            t69 = 1.3333333333333333d+0*t37*t68+4.444444444444444d-1*t64
     1         *t65+1.3333333333333333d+0*t34*t63+4.444444444444444d-1*t
     2         58*t59
            t70 = 1.923661050931536d+0*t17*t69+t57+t53+t51+t47+3.8473221
     1         01863072d+0*t28*t38
            t71 = 4.444444444444444d-1*t36*t42*t65+2.6666666666666666d+0
     1         *t11*t37*t60-2.6666666666666666d+0*t11*t34*t60+4.44444444
     2         4444444d-1*t33*t41*t59
            t72 = 1.923661050931536d+0*t17*t71+t57+t53+t51+t47+1.9236610
     1         50931536d+0*t28*t43+1.923661050931536d+0*t28*t38
            t73 = t41**2
            t74 = t67+t61
            t75 = t42**2
            t76 = t66+t62
            t77 = 1.3333333333333333d+0*t37*t76+4.444444444444444d-1*t65
     1         *t75+1.3333333333333333d+0*t34*t74+4.444444444444444d-1*t
     2         59*t73
            t78 = 1.923661050931536d+0*t17*t77+t57+t53+t51+t47+3.8473221
     1         01863072d+0*t28*t43
            t79 = t21**3
            t80 = 1/t4**4
            t81 = 8.538000000000001d-1*t7*t79*t80
            t82 = -8.538000000000001d-1*t21*t46*t50*t7
            t83 = 1/t1**3.3333333333333337d+0
            t84 = 1/t1**3.1666666666666665d+0
            t85 = -3.493762121246028d-1*t84-2.1448503343126213d-1*t83
            t86 = 1.423d-1*t22*t7*t85
            t87 = -1.2866528700135701d-3*t83*t9+5.670463005702662d-3*t83
     1         -2.0733333333333331d-2*t60
            t88 = t8*t87
            t89 = 1.923661050931536d+0*t15*(t8*(-4.503285045047496d-4*t8
     1         3*t9+2.4607236139009533d-3*t83-1.0366666666666666d-2*t60)
     2         -t8*t87-1.423d-1*t22*t7*t85+8.43d-2*t27*t7*(-4.6392143809
     3         61223d-1*t84-1.6797253218027158d-1*t83)-8.538000000000001
     4         d-1*t7*t79*t80-5.058d-1*t26*t54*t55*t7+8.538000000000001d
     5         -1*t21*t46*t50*t7+5.058d-1*t26**3*t7/t16**4)
            t90 = 1/t13**1.6666666666666669d+0
            t91 = 1/t1**4
            t92 = 6*t11*t91
            t93 = -6*t60
            t94 = 1/t14**1.6666666666666669d+0
            t95 = -6*t11*t91
            t96 = 6*t60
            t97 = 2.0d+0*t72*wght
            t98 = 3.847322101863072d+0*t28*t71
            fnc(iq) = 1.0d+0*t1*t18*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*t39*wght+t40+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*t1*t44*wght+t40+Amat(iq,D1_RB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*t1*t70*wght+2.0d+0*t39*wght+Amat
     1         2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*t1*t72*wght+1.0d+0*t44*wght+1.0d
     1         +0*t39*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*t1*t78*wght+2.0d+0*t44*wght+Amat
     1         2(iq,D2_RB_RB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*t1*(1.923661050931536d+0*t17*
     1         (1.3333333333333333d+0*t37*(t96+t95)-2.962962962962963d-1
     2         *t36**3*t94+1.3333333333333333d+0*t34*(t93+t92)-2.9629629
     3         62962963d-1*t33**3*t90+1.3333333333333333d+0*t36*t65*t68+
     4         1.3333333333333333d+0*t33*t59*t63)+t89+t88+t86+t82+t81+5.
     5         770983152794609d+0*t28*t69+5.770983152794609d+0*t38*t56)*
     6         wght+3.0d+0*t70*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 1.0d+0*t1*(t98+1.923661050931536d+0*
     1         t17*(1.3333333333333333d+0*t37*(t95+2*t60)-2.962962962962
     2         963d-1*t42*t64*t94+1.3333333333333333d+0*t34*(t92-2*t60)-
     3         2.962962962962963d-1*t41*t58*t90+4.444444444444444d-1*t42
     4         *t65*t68+1.7777777777777776d+0*t11*t36*t60*t65+4.44444444
     5         4444444d-1*t41*t59*t63-1.7777777777777776d+0*t11*t33*t59*
     6         t60)+t89+t88+t86+t82+t81+1.923661050931536d+0*t28*t69+1.9
     7         23661050931536d+0*t43*t56+3.847322101863072d+0*t38*t56)*w
     8         ght+1.0d+0*t70*wght+t97+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = 1.0d+0*t1*(t98+1.923661050931536d+0*
     1         t17*(-2.962962962962963d-1*t36*t75*t94-8*t11*t37*t91+8*t1
     2         1*t34*t91-2.962962962962963d-1*t33*t73*t90+4.444444444444
     3         444d-1*t36*t65*t76+4.444444444444444d-1*t33*t59*t74+1.777
     4         7777777777776d+0*t11*t42*t60*t65-1.7777777777777776d+0*t1
     5         1*t41*t59*t60-2.6666666666666666d+0*t37*t60+2.66666666666
     6         66666d+0*t34*t60)+t89+t88+t86+t82+t81+1.923661050931536d+
     7         0*t28*t77+3.847322101863072d+0*t43*t56+1.923661050931536d
     8         +0*t38*t56)*wght+1.0d+0*t78*wght+t97+Amat3(iq,D3_RA_RB_RB
     9         )
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*t1*(1.923661050931536d+0*t17*
     1         (1.3333333333333333d+0*t34*(t96+t92)+1.3333333333333333d+
     2         0*t37*(t95+t93)-2.962962962962963d-1*t42**3*t94-2.9629629
     3         62962963d-1*t41**3*t90+1.3333333333333333d+0*t42*t65*t76+
     4         1.3333333333333333d+0*t41*t59*t74)+t89+t88+t86+t82+t81+5.
     5         770983152794609d+0*t28*t77+5.770983152794609d+0*t43*t56)*
     6         wght+3.0d+0*t78*wght+Amat3(iq,D3_RB_RB_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**3.333333333333333d-1
            t2 = 6.203504908994d-1*t1
            t3 = nwxcm_heaviside(t2-1.0d+0)
            t4 = 1/rhoa**1.6666666666666666d-1
            t5 = 8.292885914166397d-1*t4+2.0682485366585993d-1*t1+1.0d+0
            t6 = 1/t5
            t7 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1)
            t8 = log(t2)
            t9 = 1.2407009817988d-3*t1*t8+3.11d-2*t8-7.196065694433039d-
     1         3*t1-4.8d-2
            t10 = 1.101176160755631d+0*t4+1.6197351317383332d-1*t1+1.0d+
     1         0
            t11 = 9.999999999999991d-1*(-t7*t9+t7*(4.3424534362957995d-4
     1         *t1*t8+1.555d-2*t8-2.9776823563171195d-3*t1-2.69d-2)+1.42
     2         3d-1*t3*t6-8.43d-2*t3/t10)+t7*t9-1.423d-1*t3*t6
            t12 = 1/rhoa**1.3333333333333333d+0
            t13 = 1/rhoa**1.1666666666666667d+0
            t14 = -1.382147652361066d-1*t13-6.894161788861997d-2*t12
            t15 = 1/t5**2
            t16 = -1.8352936012593848d-1*t13-5.399117105794444d-2*t12
            t17 = 1/t10**2
            t18 = 1/rhoa
            t19 = -4.1356699393293334d-4*t12*t8-1.0366666666666666d-2*t1
     1         8+1.9851215708780795d-3*t12
            t20 = 9.999999999999991d-1*(t7*(-1.4474844787652666d-4*t12*t
     1         8-5.183333333333333d-3*t18+8.478123375625131d-4*t12)-t19*
     2         t7+8.43d-2*t16*t17*t3-1.423d-1*t14*t15*t3)+t19*t7+1.423d-
     3         1*t14*t15*t3
            t21 = t14**2
            t22 = 1/t5**3
            t23 = 1/rhoa**2.3333333333333334d+0
            t24 = 1/rhoa**2.1666666666666665d+0
            t25 = 1.6125055944212438d-1*t24+9.192215718482661d-2*t23
            t26 = 1/t10**3
            t27 = 2.141175868135949d-1*t24+7.198822807725924d-2*t23
            t28 = 1/rhoa**2
            t29 = 5.514226585772444d-4*t23*t8+1.0366666666666666d-2*t28-
     1         2.5089730965264617d-3*t23
            t30 = 9.999999999999991d-1*(t7*(1.9299793050203554d-4*t23*t8
     1         +5.183333333333333d-3*t28-1.082166967457842d-3*t23)-t29*t
     2         7+8.43d-2*t17*t27*t3-1.686d-1*t16**2*t26*t3-1.423d-1*t15*
     3         t25*t3+2.846d-1*t21*t22*t3)+t29*t7+1.423d-1*t15*t25*t3-2.
     4         846d-1*t21*t22*t3
            t31 = t14**3
            t32 = 1/t5**4
            t33 = 1/rhoa**3.3333333333333337d+0
            t34 = 1/rhoa**3.1666666666666665d+0
            t35 = -3.493762121246028d-1*t34-2.1448503343126213d-1*t33
            t36 = 1/rhoa**3
            t37 = -1.2866528700135701d-3*t33*t8-2.0733333333333331d-2*t3
     1         6+5.670463005702662d-3*t33
            fnc(iq) = 1.0d+0*rhoa*t11*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*rhoa*t20*wght+1.0d+0*t11*wght+Amat(i
     1         q,D1_RA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*rhoa*t30*wght+2.0d+0*t20*wght+Am
     1         at2(iq,D2_RA_RA)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*rhoa*(9.999999999999991d-1*(t
     1         7*(-4.503285045047496d-4*t33*t8-1.0366666666666666d-2*t36
     2         +2.4607236139009533d-3*t33)-t37*t7-1.423d-1*t15*t3*t35+8.
     3         43d-2*t17*t3*(-4.639214380961223d-1*t34-1.679725321802715
     4         8d-1*t33)-8.538000000000001d-1*t3*t31*t32-5.058d-1*t16*t2
     5         6*t27*t3+8.538000000000001d-1*t14*t22*t25*t3+5.058d-1*t16
     6         **3*t3/t10**4)+t37*t7+1.423d-1*t15*t3*t35+8.5380000000000
     7         01d-1*t3*t31*t32-8.538000000000001d-1*t14*t22*t25*t3)*wgh
     8         t+3.0d+0*t30*wght+Amat3(iq,D3_RA_RA_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**3.333333333333333d-1
            t2 = 6.203504908994d-1*t1
            t3 = nwxcm_heaviside(t2-1.0d+0)
            t4 = 1/rhob**1.6666666666666666d-1
            t5 = 8.292885914166397d-1*t4+2.0682485366585993d-1*t1+1.0d+0
            t6 = 1/t5
            t7 = nwxcm_heaviside(1.0d+0-6.203504908994d-1*t1)
            t8 = log(t2)
            t9 = 1.2407009817988d-3*t1*t8+3.11d-2*t8-7.196065694433039d-
     1         3*t1-4.8d-2
            t10 = 1.101176160755631d+0*t4+1.6197351317383332d-1*t1+1.0d+
     1         0
            t11 = 9.999999999999991d-1*(-t7*t9+t7*(4.3424534362957995d-4
     1         *t1*t8+1.555d-2*t8-2.9776823563171195d-3*t1-2.69d-2)+1.42
     2         3d-1*t3*t6-8.43d-2*t3/t10)+t7*t9-1.423d-1*t3*t6
            t12 = 1/rhob**1.3333333333333333d+0
            t13 = 1/rhob**1.1666666666666667d+0
            t14 = -1.382147652361066d-1*t13-6.894161788861997d-2*t12
            t15 = 1/t5**2
            t16 = -1.8352936012593848d-1*t13-5.399117105794444d-2*t12
            t17 = 1/t10**2
            t18 = 1/rhob
            t19 = -4.1356699393293334d-4*t12*t8-1.0366666666666666d-2*t1
     1         8+1.9851215708780795d-3*t12
            t20 = 9.999999999999991d-1*(t7*(-1.4474844787652666d-4*t12*t
     1         8-5.183333333333333d-3*t18+8.478123375625131d-4*t12)-t19*
     2         t7+8.43d-2*t16*t17*t3-1.423d-1*t14*t15*t3)+t19*t7+1.423d-
     3         1*t14*t15*t3
            t21 = t14**2
            t22 = 1/t5**3
            t23 = 1/rhob**2.3333333333333334d+0
            t24 = 1/rhob**2.1666666666666665d+0
            t25 = 1.6125055944212438d-1*t24+9.192215718482661d-2*t23
            t26 = 1/t10**3
            t27 = 2.141175868135949d-1*t24+7.198822807725924d-2*t23
            t28 = 1/rhob**2
            t29 = 5.514226585772444d-4*t23*t8+1.0366666666666666d-2*t28-
     1         2.5089730965264617d-3*t23
            t30 = 9.999999999999991d-1*(t7*(1.9299793050203554d-4*t23*t8
     1         +5.183333333333333d-3*t28-1.082166967457842d-3*t23)-t29*t
     2         7+8.43d-2*t17*t27*t3-1.686d-1*t16**2*t26*t3-1.423d-1*t15*
     3         t25*t3+2.846d-1*t21*t22*t3)+t29*t7+1.423d-1*t15*t25*t3-2.
     4         846d-1*t21*t22*t3
            t31 = t14**3
            t32 = 1/t5**4
            t33 = 1/rhob**3.3333333333333337d+0
            t34 = 1/rhob**3.1666666666666665d+0
            t35 = -3.493762121246028d-1*t34-2.1448503343126213d-1*t33
            t36 = 1/rhob**3
            t37 = -1.2866528700135701d-3*t33*t8-2.0733333333333331d-2*t3
     1         6+5.670463005702662d-3*t33
            fnc(iq) = 1.0d+0*rhob*t11*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*rhob*t20*wght+1.0d+0*t11*wght+Amat(i
     1         q,D1_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*rhob*t30*wght+2.0d+0*t20*wght+Am
     1         at2(iq,D2_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*rhob*(9.999999999999991d-1*(t
     1         7*(-4.503285045047496d-4*t33*t8-1.0366666666666666d-2*t36
     2         +2.4607236139009533d-3*t33)-t37*t7-1.423d-1*t15*t3*t35+8.
     3         43d-2*t17*t3*(-4.639214380961223d-1*t34-1.679725321802715
     4         8d-1*t33)-8.538000000000001d-1*t3*t31*t32-5.058d-1*t16*t2
     5         6*t27*t3+8.538000000000001d-1*t14*t22*t25*t3+5.058d-1*t16
     6         **3*t3/t10**4)+t37*t7+1.423d-1*t15*t3*t35+8.5380000000000
     7         01d-1*t3*t31*t32-8.538000000000001d-1*t14*t22*t25*t3)*wgh
     8         t+3.0d+0*t30*wght+Amat3(iq,D3_RB_RB_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
