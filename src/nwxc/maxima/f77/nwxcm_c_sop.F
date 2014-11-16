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
C>   {\it t_4} &=& {{1}\over{1.8610514726982\,{\it t_3}+1.8610514726982\,{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \left({\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}+0.5764\right)\,\rho_\beta}\over{{\it t_5}\,\left(3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}+11.99591941974898\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111, 5656 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
      subroutine nwxcm_c_sop(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = rhoa**3.333333333333333d-1
            t3 = t1*t2
            t4 = t3**2.0d+0
            t5 = 1/t4
            t6 = 8.658781460080349d-1*t4+1.0500052408963245d+0*t1*t2+3.1
     1         83d-1
            t7 = 1/t6
            t8 = 1.415701855281521d+0*t1*t2+5.764d-1
            t9 = rhoa**2
            t10 = rhoa**1.3333333333333333d+0
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t5*t7*t8*t9*wght
            Amat(iq,D1_RA) = t5*t7*(-2.3097938309456317d+0*rhoa*t8*wght-
     1         5.449965686312571d-1*t1*t10*wght)+1.1548969154728159d+0*(
     2         5.772520973386899d-1*t1**2/t2+3.500017469654415d-1*t1/rho
     3         a**6.666666666666666d-1)*t5*t8*t9*wght/t6**2+7.6993127698
     4         18772d-1*t1*t10*t7*t8*wght/t3**3.0d+0+Amat(iq,D1_RA)
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
C>   {\it t_4} &=& {{1}\over{1.8610514726982\,{\it t_3}+1.8610514726982\,{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \left({\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}+0.5764\right)\,\rho_\beta}\over{{\it t_5}\,\left(3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}+11.99591941974898\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111, 5656 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
      subroutine nwxcm_c_sop_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_c_sop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{1.8610514726982\,{\it t_3}+1.8610514726982\,{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \left({\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}+0.5764\right)\,\rho_\beta}\over{{\it t_5}\,\left(3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_4}\,{\it t_3}+11.99591941974898\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111, 5656 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
      subroutine nwxcm_c_sop_d3(param,tol_rho,ipol,nq,wght, 
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
