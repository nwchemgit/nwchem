C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_m06_2x.F 
C> The nwxcm_x_m06_2x functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_m06_2x functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_4} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_5} &=& {\it t_4}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_6} &=& 2.0\,\tau_\alpha+{\it t_4}\\\\ 
C>   {\it t_7} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{16}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_{19}} &=& {\it t_{18}}-2.0\,\tau_\beta\\\\ 
C>   {\it t_{20}} &=& 2.0\,\tau_\beta+{\it t_{18}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{22}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{23}} &=& {\it t_{22}}-2.0\,\tau_s\\\\ 
C>   {\it t_{24}} &=& 2.0\,\tau_s+{\it t_{22}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{17}}\,
C>    \sigma_{\beta\beta}}\over{0.00449276922095889\,{\it t_{17}}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left({{{\it t_{16}}
C>    \,{\it t_{19}}}\over{{\it t_{20}}}}+{{{\it t_{15}}\,{
C>    \it t_{19}}^2}\over{{\it t_{20}}^2}}+{{{\it t_{14}}\,{
C>    \it t_{19}}^3}\over{{\it t_{20}}^3}}+{{{\it t_{13}}\,{
C>    \it t_{19}}^4}\over{{\it t_{20}}^4}}+{{{\it t_{12}}\,{
C>    \it t_{19}}^5}\over{{\it t_{20}}^5}}+{{{\it t_{11}}\,{
C>    \it t_{19}}^6}\over{{\it t_{20}}^6}}+{{{\it t_{10}}\,{
C>    \it t_{19}}^7}\over{{\it t_{20}}^7}}+{{{\it t_9}\,{
C>    \it t_{19}}^8}\over{{\it t_{20}}^8}}+{{{\it t_8}\,{
C>    \it t_{19}}^9}\over{{\it t_{20}}^9}}+{{{\it t_7}\,{
C>    \it t_{19}}^{10}}\over{{\it t_{20}}^{10}}}+{{{\it t_3}\,{
C>    \it t_{19}}^{11}}\over{{\it t_{20}}^{11}}}+{\it t_2}\right)
C>    -0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.00449276922095889\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left({{{\it t_{16}}
C>    \,{\it t_5}}\over{{\it t_6}}}+{{{\it t_{15}}\,{
C>    \it t_5}^2}\over{{\it t_6}^2}}+{{{\it t_{14}}\,{
C>    \it t_5}^3}\over{{\it t_6}^3}}+{{{\it t_{13}}\,{
C>    \it t_5}^4}\over{{\it t_6}^4}}+{{{\it t_{12}}\,{
C>    \it t_5}^5}\over{{\it t_6}^5}}+{{{\it t_{11}}\,{
C>    \it t_5}^6}\over{{\it t_6}^6}}+{{{\it t_{10}}\,{
C>    \it t_5}^7}\over{{\it t_6}^7}}+{{{\it t_9}\,{\it t_5}^8}
C>    \over{{\it t_6}^8}}+{{{\it t_8}\,{\it t_5}^9}\over{{\it t_6}^9}}
C>    +{{{\it t_7}\,{\it t_5}^{10}}\over{{\it t_6}^{10}}}+{{{
C>    \it t_3}\,{\it t_5}^{11}}\over{{\it t_6}^{11}}}+{
C>    \it t_2}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{21}}\,\sigma_{ss}}
C>    \over{0.00449276922095889\,{\it t_{21}}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left({{{\it t_{16}}\,{\it t_{23}}}\over{{
C>    \it t_{24}}}}+{{{\it t_{15}}\,{\it t_{23}}^2}\over{{
C>    \it t_{24}}^2}}+{{{\it t_{14}}\,{\it t_{23}}^3}\over{{
C>    \it t_{24}}^3}}+{{{\it t_{13}}\,{\it t_{23}}^4}\over{{
C>    \it t_{24}}^4}}+{{{\it t_{12}}\,{\it t_{23}}^5}\over{{
C>    \it t_{24}}^5}}+{{{\it t_{11}}\,{\it t_{23}}^6}\over{{
C>    \it t_{24}}^6}}+{{{\it t_{10}}\,{\it t_{23}}^7}\over{{
C>    \it t_{24}}^7}}+{{{\it t_9}\,{\it t_{23}}^8}\over{{
C>    \it t_{24}}^8}}+{{{\it t_8}\,{\it t_{23}}^9}\over{{
C>    \it t_{24}}^9}}+{{{\it t_7}\,{\it t_{23}}^{10}}\over{{
C>    \it t_{24}}^{10}}}+{{{\it t_3}\,{\it t_{23}}^{11}}\over{{
C>    \it t_{24}}^{11}}}+{\it t_2}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 123, 194101 (2006)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2370993 "> 
C> 10.1063/1.2370993 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_m06_2x(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(11)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(10)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(9)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(8)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(7)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(6)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(5)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(4)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(3)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(2)
            t40 = 1/t10
            t41 = t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t30*t31*t
     1         32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t20+t15*t1
     2         6*t17+t12*t13*t14+param(1)
            t42 = rhoa**6.666666666666666d-1
            t43 = 1/t10**12
            t44 = 1/t2**2
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t4*t41*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-1.6711932865
     1         267184d+2*t42*t43*t6*t9-1.5192666241151986d+1*t38*t39*t42
     2         *t8+3.038533248230397d+1*t36*t38*t42*t8+1.671193286526718
     3         4d+2*t11*t13*t42*t6+1.5192666241151986d+1*t39*t40*t42-3.0
     4         38533248230397d+1*t35*t36*t37*t42+4.5577998723455965d+1*t
     5         33*t35*t37*t42-4.5577998723455965d+1*t32*t33*t34*t42+6.07
     6         7066496460795d+1*t30*t32*t34*t42-6.077066496460795d+1*t29
     7         *t30*t31*t42+7.596333120575993d+1*t27*t29*t31*t42-7.59633
     8         3120575993d+1*t26*t27*t28*t42+9.115599744691193d+1*t24*t2
     9         6*t28*t42-9.115599744691193d+1*t23*t24*t25*t42+1.06348663
     :         68806391d+2*t21*t23*t25*t42-1.0634866368806391d+2*t20*t21
     ;         *t22*t42+1.2154132992921589d+2*t18*t20*t22*t42-1.21541329
     <         92921589d+2*t17*t18*t19*t42+1.3673399617036788d+2*t15*t17
     =         *t19*t42-1.3673399617036788d+2*t14*t15*t16*t42+1.51926662
     >         41151986d+2*t12*t14*t16*t42-1.5192666241151986d+2*t11*t12
     ?         *t13*t42)-9.305257363490997d-1*t41*(4.327658698487368d-5*
     @         gammaaa**2*t44/rhoa**6.333333333333333d+0-9.6324972097358
     1         63d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)*t5-1.240700
     2         9817987995d+0*rhoa**3.333333333333333d-1*t4*t41)*wght+Ama
     3         t(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t41*(
     1         3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gammaaa
     2         *t44/rhoa**5.333333333333333d+0)*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         -2.2d+1*t43*t6*t9+t38*(-2.0d+0*t39*t8-4.0d+0*t36*t8)+t11*
     2         (-2.2d+1*t13*t6-2.0d+1*t12*t13)-2.0d+0*t39*t40+t35*(-4.0d
     3         +0*t36*t37-6.0d+0*t33*t37)+t32*(-6.0d+0*t33*t34-8.0d+0*t3
     4         0*t34)+t29*(-8.0d+0*t30*t31-1.0d+1*t27*t31)+t26*(-1.0d+1*
     5         t27*t28-1.2d+1*t24*t28)+t23*(-1.2d+1*t24*t25-1.4d+1*t21*t
     6         25)+t20*(-1.4d+1*t21*t22-1.6d+1*t18*t22)+t17*(-1.6d+1*t18
     7         *t19-1.8d+1*t15*t19)+t14*(-1.8d+1*t15*t16-2.0d+1*t12*t16)
     8         )*wght
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
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(1)
            t7 = param(12)
            t8 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(11)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(10)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(9)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(8)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(7)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(6)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(5)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(4)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(3)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(2)
            t41 = 1/t11
            t42 = t40*t41*t9+t10*t12*t7+t6+t37*t38*t39+t34*t35*t36+t31*t
     1         32*t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t21+t1
     2         6*t17*t18+t13*t14*t15
            t43 = 1/rhob**2.6666666666666666d+0
            t44 = 4.49276922095889d-3*gammabb*t43+1.0d+0
            t45 = 1/t44
            t46 = 3.6121864536509485d-3*gammabb*t43*t45+1.0d+0
            t47 = rhob**1.3333333333333333d+0
            t48 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t49 = t48-2.0d+0*taub
            t50 = t49**11
            t51 = 2.0d+0*taub+t48
            t52 = 1/t51**11
            t53 = t49**10
            t54 = 1/t51**10
            t55 = t49**9
            t56 = 1/t51**9
            t57 = t49**8
            t58 = 1/t51**8
            t59 = t49**7
            t60 = 1/t51**7
            t61 = t49**6
            t62 = 1/t51**6
            t63 = t49**5
            t64 = 1/t51**5
            t65 = t49**4
            t66 = 1/t51**4
            t67 = t49**3
            t68 = 1/t51**3
            t69 = t49**2
            t70 = 1/t51**2
            t71 = 1/t51
            t72 = t40*t49*t71+t37*t69*t70+t50*t52*t7+t34*t67*t68+t31*t65
     1         *t66+t28*t63*t64+t25*t61*t62+t22*t59*t60+t6+t19*t57*t58+t
     2         16*t55*t56+t13*t53*t54
            t73 = rhoa**6.666666666666666d-1
            t74 = 1/t11**12
            t75 = 1/t2**2
            t76 = rhob**6.666666666666666d-1
            t77 = 1/t51**12
            t78 = 1/t44**2
            fnc(iq) = (-9.305257363490997d-1*t46*t47*t72-9.3052573634909
     1         97d-1*t4*t42*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-1.5192666241
     1         151986d+1*t39*t40*t73*t9+3.038533248230397d+1*t37*t39*t73
     2         *t9-1.6711932865267184d+2*t10*t7*t73*t74+1.67119328652671
     3         84d+2*t12*t14*t7*t73+1.5192666241151986d+1*t40*t41*t73-3.
     4         038533248230397d+1*t36*t37*t38*t73+4.5577998723455965d+1*
     5         t34*t36*t38*t73-4.5577998723455965d+1*t33*t34*t35*t73+6.0
     6         77066496460795d+1*t31*t33*t35*t73-6.077066496460795d+1*t3
     7         0*t31*t32*t73+7.596333120575993d+1*t28*t30*t32*t73-7.5963
     8         33120575993d+1*t27*t28*t29*t73+9.115599744691193d+1*t25*t
     9         27*t29*t73-9.115599744691193d+1*t24*t25*t26*t73+1.0634866
     :         368806391d+2*t22*t24*t26*t73-1.0634866368806391d+2*t21*t2
     ;         2*t23*t73+1.2154132992921589d+2*t19*t21*t23*t73-1.2154132
     <         992921589d+2*t18*t19*t20*t73+1.3673399617036788d+2*t16*t1
     =         8*t20*t73-1.3673399617036788d+2*t15*t16*t17*t73+1.5192666
     >         241151986d+2*t13*t15*t17*t73-1.5192666241151986d+2*t12*t1
     ?         3*t14*t73)-9.305257363490997d-1*t42*t5*(4.327658698487368
     @         d-5*gammaaa**2*t75/rhoa**6.333333333333333d+0-9.632497209
     1         735863d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)-1.24070
     2         09817987995d+0*rhoa**3.333333333333333d-1*t4*t42)*wght+Am
     3         at(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t47*t72*(4.327658698
     1         487368d-5*gammabb**2*t78/rhob**6.333333333333333d+0-9.632
     2         497209735863d-3*gammabb*t45/rhob**3.6666666666666664d+0)-
     3         9.305257363490997d-1*t46*t47*(-1.6711932865267184d+2*t50*
     4         t7*t76*t77+1.5192666241151986d+1*t40*t71*t76-1.5192666241
     5         151986d+1*t40*t49*t70*t76+3.038533248230397d+1*t37*t49*t7
     6         0*t76+1.6711932865267184d+2*t52*t53*t7*t76-3.038533248230
     7         397d+1*t37*t68*t69*t76+4.5577998723455965d+1*t34*t68*t69*
     8         t76-4.5577998723455965d+1*t34*t66*t67*t76+6.0770664964607
     9         95d+1*t31*t66*t67*t76-6.077066496460795d+1*t31*t64*t65*t7
     :         6+7.596333120575993d+1*t28*t64*t65*t76-7.596333120575993d
     ;         +1*t28*t62*t63*t76+9.115599744691193d+1*t25*t62*t63*t76-9
     <         .115599744691193d+1*t25*t60*t61*t76+1.0634866368806391d+2
     =         *t22*t60*t61*t76-1.0634866368806391d+2*t22*t58*t59*t76+1.
     >         2154132992921589d+2*t19*t58*t59*t76-1.2154132992921589d+2
     ?         *t19*t56*t57*t76+1.3673399617036788d+2*t16*t56*t57*t76-1.
     @         3673399617036788d+2*t16*t54*t55*t76+1.5192666241151986d+2
     1         *t13*t54*t55*t76-1.5192666241151986d+2*t13*t52*t53*t76)-1
     2         .2407009817987995d+0*rhob**3.333333333333333d-1*t46*t72)*
     3         wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t42*t
     1         5*(3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     2         aaa*t75/rhoa**5.333333333333333d+0)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t47*t
     1         72*(3.6121864536509485d-3*t43*t45-1.6228720119327628d-5*g
     2         ammabb*t78/rhob**5.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         t39*(-2.0d+0*t40*t9-4.0d+0*t37*t9)-2.2d+1*t10*t7*t74+t12*
     2         (-2.2d+1*t14*t7-2.0d+1*t13*t14)-2.0d+0*t40*t41+t36*(-4.0d
     3         +0*t37*t38-6.0d+0*t34*t38)+t33*(-6.0d+0*t34*t35-8.0d+0*t3
     4         1*t35)+t30*(-8.0d+0*t31*t32-1.0d+1*t28*t32)+t27*(-1.0d+1*
     5         t28*t29-1.2d+1*t25*t29)+t24*(-1.2d+1*t25*t26-1.4d+1*t22*t
     6         26)+t21*(-1.4d+1*t22*t23-1.6d+1*t19*t23)+t18*(-1.6d+1*t19
     7         *t20-1.8d+1*t16*t20)+t15*(-1.8d+1*t16*t17-2.0d+1*t13*t17)
     8         )*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t46*t47
     1         *(-2.2d+1*t50*t7*t77-2.0d+0*t40*t71+(-2.0d+0*t40*t49-4.0d
     2         +0*t37*t49)*t70+t52*(-2.2d+1*t53*t7-2.0d+1*t13*t53)+t68*(
     3         -4.0d+0*t37*t69-6.0d+0*t34*t69)+t66*(-6.0d+0*t34*t67-8.0d
     4         +0*t31*t67)+t64*(-8.0d+0*t31*t65-1.0d+1*t28*t65)+t62*(-1.
     5         0d+1*t28*t63-1.2d+1*t25*t63)+t60*(-1.2d+1*t25*t61-1.4d+1*
     6         t22*t61)+t58*(-1.4d+1*t22*t59-1.6d+1*t19*t59)+t56*(-1.6d+
     7         1*t19*t57-1.8d+1*t16*t57)+t54*(-1.8d+1*t16*t55-2.0d+1*t13
     8         *t55))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(11)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(10)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(9)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(8)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(7)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(6)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(5)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(4)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(3)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(2)
            t40 = 1/t10
            t41 = t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t30*t31*t
     1         32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t20+t15*t1
     2         6*t17+t12*t13*t14+param(1)
            t42 = rhoa**6.666666666666666d-1
            t43 = 1/t10**12
            t44 = 1/t2**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t41*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t4*t5*(-1.67119328652
     1         67184d+2*t42*t43*t6*t9-1.5192666241151986d+1*t38*t39*t42*
     2         t8+3.038533248230397d+1*t36*t38*t42*t8+1.6711932865267184
     3         d+2*t11*t13*t42*t6+1.5192666241151986d+1*t39*t40*t42-3.03
     4         8533248230397d+1*t35*t36*t37*t42+4.5577998723455965d+1*t3
     5         3*t35*t37*t42-4.5577998723455965d+1*t32*t33*t34*t42+6.077
     6         066496460795d+1*t30*t32*t34*t42-6.077066496460795d+1*t29*
     7         t30*t31*t42+7.596333120575993d+1*t27*t29*t31*t42-7.596333
     8         120575993d+1*t26*t27*t28*t42+9.115599744691193d+1*t24*t26
     9         *t28*t42-9.115599744691193d+1*t23*t24*t25*t42+1.063486636
     :         8806391d+2*t21*t23*t25*t42-1.0634866368806391d+2*t20*t21*
     ;         t22*t42+1.2154132992921589d+2*t18*t20*t22*t42-1.215413299
     <         2921589d+2*t17*t18*t19*t42+1.3673399617036788d+2*t15*t17*
     =         t19*t42-1.3673399617036788d+2*t14*t15*t16*t42+1.519266624
     >         1151986d+2*t12*t14*t16*t42-1.5192666241151986d+2*t11*t12*
     ?         t13*t42)*wght-9.305257363490997d-1*t41*(4.327658698487368
     @         d-5*gammaaa**2*t44/rhoa**6.333333333333333d+0-9.632497209
     1         735863d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)*t5*wght
     2         -1.2407009817987995d+0*rhoa**3.333333333333333d-1*t4*t41*
     3         wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t41*(
     1         3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gammaaa
     2         *t44/rhoa**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         -2.2d+1*t43*t6*t9+t38*(-2.0d+0*t39*t8-4.0d+0*t36*t8)+t11*
     2         (-2.2d+1*t13*t6-2.0d+1*t12*t13)-2.0d+0*t39*t40+t35*(-4.0d
     3         +0*t36*t37-6.0d+0*t33*t37)+t32*(-6.0d+0*t33*t34-8.0d+0*t3
     4         0*t34)+t29*(-8.0d+0*t30*t31-1.0d+1*t27*t31)+t26*(-1.0d+1*
     5         t27*t28-1.2d+1*t24*t28)+t23*(-1.2d+1*t24*t25-1.4d+1*t21*t
     6         25)+t20*(-1.4d+1*t21*t22-1.6d+1*t18*t22)+t17*(-1.6d+1*t18
     7         *t19-1.8d+1*t15*t19)+t14*(-1.8d+1*t15*t16-2.0d+1*t12*t16)
     8         )*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = param(12)
            t7 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t8 = t7-2.0d+0*taub
            t9 = t8**11
            t10 = 2.0d+0*taub+t7
            t11 = 1/t10**11
            t12 = param(11)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(10)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(9)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(8)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(7)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(6)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(5)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(4)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(3)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(2)
            t40 = 1/t10
            t41 = t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t30*t31*t
     1         32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t20+t15*t1
     2         6*t17+t12*t13*t14+param(1)
            t42 = rhob**6.666666666666666d-1
            t43 = 1/t10**12
            t44 = 1/t2**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t41*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t4*t5*(-1.67119328652
     1         67184d+2*t42*t43*t6*t9-1.5192666241151986d+1*t38*t39*t42*
     2         t8+3.038533248230397d+1*t36*t38*t42*t8+1.6711932865267184
     3         d+2*t11*t13*t42*t6+1.5192666241151986d+1*t39*t40*t42-3.03
     4         8533248230397d+1*t35*t36*t37*t42+4.5577998723455965d+1*t3
     5         3*t35*t37*t42-4.5577998723455965d+1*t32*t33*t34*t42+6.077
     6         066496460795d+1*t30*t32*t34*t42-6.077066496460795d+1*t29*
     7         t30*t31*t42+7.596333120575993d+1*t27*t29*t31*t42-7.596333
     8         120575993d+1*t26*t27*t28*t42+9.115599744691193d+1*t24*t26
     9         *t28*t42-9.115599744691193d+1*t23*t24*t25*t42+1.063486636
     :         8806391d+2*t21*t23*t25*t42-1.0634866368806391d+2*t20*t21*
     ;         t22*t42+1.2154132992921589d+2*t18*t20*t22*t42-1.215413299
     <         2921589d+2*t17*t18*t19*t42+1.3673399617036788d+2*t15*t17*
     =         t19*t42-1.3673399617036788d+2*t14*t15*t16*t42+1.519266624
     >         1151986d+2*t12*t14*t16*t42-1.5192666241151986d+2*t11*t12*
     ?         t13*t42)*wght-9.305257363490997d-1*t41*(4.327658698487368
     @         d-5*gammabb**2*t44/rhob**6.333333333333333d+0-9.632497209
     1         735863d-3*gammabb*t3/rhob**3.6666666666666664d+0)*t5*wght
     2         -1.2407009817987995d+0*rhob**3.333333333333333d-1*t4*t41*
     3         wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t41*(
     1         3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gammabb
     2         *t44/rhob**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t4*t5*(
     1         -2.2d+1*t43*t6*t9+t38*(-2.0d+0*t39*t8-4.0d+0*t36*t8)+t11*
     2         (-2.2d+1*t13*t6-2.0d+1*t12*t13)-2.0d+0*t39*t40+t35*(-4.0d
     3         +0*t36*t37-6.0d+0*t33*t37)+t32*(-6.0d+0*t33*t34-8.0d+0*t3
     4         0*t34)+t29*(-8.0d+0*t30*t31-1.0d+1*t27*t31)+t26*(-1.0d+1*
     5         t27*t28-1.2d+1*t24*t28)+t23*(-1.2d+1*t24*t25-1.4d+1*t21*t
     6         25)+t20*(-1.4d+1*t21*t22-1.6d+1*t18*t22)+t17*(-1.6d+1*t18
     7         *t19-1.8d+1*t15*t19)+t14*(-1.8d+1*t15*t16-2.0d+1*t12*t16)
     8         )*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_m06_2x functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_4} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_5} &=& {\it t_4}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_6} &=& 2.0\,\tau_\alpha+{\it t_4}\\\\ 
C>   {\it t_7} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{16}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_{19}} &=& {\it t_{18}}-2.0\,\tau_\beta\\\\ 
C>   {\it t_{20}} &=& 2.0\,\tau_\beta+{\it t_{18}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{22}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{23}} &=& {\it t_{22}}-2.0\,\tau_s\\\\ 
C>   {\it t_{24}} &=& 2.0\,\tau_s+{\it t_{22}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{17}}\,
C>    \sigma_{\beta\beta}}\over{0.00449276922095889\,{\it t_{17}}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left({{{\it t_{16}}
C>    \,{\it t_{19}}}\over{{\it t_{20}}}}+{{{\it t_{15}}\,{
C>    \it t_{19}}^2}\over{{\it t_{20}}^2}}+{{{\it t_{14}}\,{
C>    \it t_{19}}^3}\over{{\it t_{20}}^3}}+{{{\it t_{13}}\,{
C>    \it t_{19}}^4}\over{{\it t_{20}}^4}}+{{{\it t_{12}}\,{
C>    \it t_{19}}^5}\over{{\it t_{20}}^5}}+{{{\it t_{11}}\,{
C>    \it t_{19}}^6}\over{{\it t_{20}}^6}}+{{{\it t_{10}}\,{
C>    \it t_{19}}^7}\over{{\it t_{20}}^7}}+{{{\it t_9}\,{
C>    \it t_{19}}^8}\over{{\it t_{20}}^8}}+{{{\it t_8}\,{
C>    \it t_{19}}^9}\over{{\it t_{20}}^9}}+{{{\it t_7}\,{
C>    \it t_{19}}^{10}}\over{{\it t_{20}}^{10}}}+{{{\it t_3}\,{
C>    \it t_{19}}^{11}}\over{{\it t_{20}}^{11}}}+{\it t_2}\right)
C>    -0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.00449276922095889\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left({{{\it t_{16}}
C>    \,{\it t_5}}\over{{\it t_6}}}+{{{\it t_{15}}\,{
C>    \it t_5}^2}\over{{\it t_6}^2}}+{{{\it t_{14}}\,{
C>    \it t_5}^3}\over{{\it t_6}^3}}+{{{\it t_{13}}\,{
C>    \it t_5}^4}\over{{\it t_6}^4}}+{{{\it t_{12}}\,{
C>    \it t_5}^5}\over{{\it t_6}^5}}+{{{\it t_{11}}\,{
C>    \it t_5}^6}\over{{\it t_6}^6}}+{{{\it t_{10}}\,{
C>    \it t_5}^7}\over{{\it t_6}^7}}+{{{\it t_9}\,{\it t_5}^8}
C>    \over{{\it t_6}^8}}+{{{\it t_8}\,{\it t_5}^9}\over{{\it t_6}^9}}
C>    +{{{\it t_7}\,{\it t_5}^{10}}\over{{\it t_6}^{10}}}+{{{
C>    \it t_3}\,{\it t_5}^{11}}\over{{\it t_6}^{11}}}+{
C>    \it t_2}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{21}}\,\sigma_{ss}}
C>    \over{0.00449276922095889\,{\it t_{21}}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left({{{\it t_{16}}\,{\it t_{23}}}\over{{
C>    \it t_{24}}}}+{{{\it t_{15}}\,{\it t_{23}}^2}\over{{
C>    \it t_{24}}^2}}+{{{\it t_{14}}\,{\it t_{23}}^3}\over{{
C>    \it t_{24}}^3}}+{{{\it t_{13}}\,{\it t_{23}}^4}\over{{
C>    \it t_{24}}^4}}+{{{\it t_{12}}\,{\it t_{23}}^5}\over{{
C>    \it t_{24}}^5}}+{{{\it t_{11}}\,{\it t_{23}}^6}\over{{
C>    \it t_{24}}^6}}+{{{\it t_{10}}\,{\it t_{23}}^7}\over{{
C>    \it t_{24}}^7}}+{{{\it t_9}\,{\it t_{23}}^8}\over{{
C>    \it t_{24}}^8}}+{{{\it t_8}\,{\it t_{23}}^9}\over{{
C>    \it t_{24}}^9}}+{{{\it t_7}\,{\it t_{23}}^{10}}\over{{
C>    \it t_{24}}^{10}}}+{{{\it t_3}\,{\it t_{23}}^{11}}\over{{
C>    \it t_{24}}^{11}}}+{\it t_2}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 123, 194101 (2006)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2370993 "> 
C> 10.1063/1.2370993 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_m06_2x_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(11)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(10)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(9)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(8)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(7)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(6)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(5)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(4)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(3)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(2)
            t40 = 1/t10
            t41 = t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t30*t31*t
     1         32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t20+t15*t1
     2         6*t17+t12*t13*t14+param(1)
            t42 = rhoa**6.666666666666666d-1
            t43 = 1/t10**12
            t44 = -1.6711932865267184d+2*t42*t43*t6*t9-1.519266624115198
     1         6d+1*t38*t39*t42*t8+3.038533248230397d+1*t36*t38*t42*t8+1
     2         .6711932865267184d+2*t11*t13*t42*t6+1.5192666241151986d+1
     3         *t39*t40*t42-3.038533248230397d+1*t35*t36*t37*t42+4.55779
     4         98723455965d+1*t33*t35*t37*t42-4.5577998723455965d+1*t32*
     5         t33*t34*t42+6.077066496460795d+1*t30*t32*t34*t42-6.077066
     6         496460795d+1*t29*t30*t31*t42+7.596333120575993d+1*t27*t29
     7         *t31*t42-7.596333120575993d+1*t26*t27*t28*t42+9.115599744
     8         691193d+1*t24*t26*t28*t42-9.115599744691193d+1*t23*t24*t2
     9         5*t42+1.0634866368806391d+2*t21*t23*t25*t42-1.06348663688
     :         06391d+2*t20*t21*t22*t42+1.2154132992921589d+2*t18*t20*t2
     ;         2*t42-1.2154132992921589d+2*t17*t18*t19*t42+1.36733996170
     <         36788d+2*t15*t17*t19*t42-1.3673399617036788d+2*t14*t15*t1
     =         6*t42+1.5192666241151986d+2*t12*t14*t16*t42-1.51926662411
     >         51986d+2*t11*t12*t13*t42
            t45 = rhoa**3.333333333333333d-1
            t46 = gammaaa**2
            t47 = 1/t2**2
            t48 = 1/rhoa**6.333333333333333d+0
            t49 = 1/rhoa**3.6666666666666664d+0
            t50 = 4.327658698487368d-5*t46*t47*t48-9.632497209735863d-3*
     1         gammaaa*t3*t49
            t51 = 1/rhoa**5.333333333333333d+0
            t52 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t47*t51
            t53 = -2.2d+1*t43*t6*t9
            t54 = -2.2d+1*t13*t6-2.0d+1*t12*t13
            t55 = -1.8d+1*t15*t16-2.0d+1*t12*t16
            t56 = -1.6d+1*t18*t19-1.8d+1*t15*t19
            t57 = -1.4d+1*t21*t22-1.6d+1*t18*t22
            t58 = -1.2d+1*t24*t25-1.4d+1*t21*t25
            t59 = -1.0d+1*t27*t28-1.2d+1*t24*t28
            t60 = -8.0d+0*t30*t31-1.0d+1*t27*t31
            t61 = -6.0d+0*t33*t34-8.0d+0*t30*t34
            t62 = -4.0d+0*t36*t37-6.0d+0*t33*t37
            t63 = -2.0d+0*t39*t8-4.0d+0*t36*t8
            t64 = -2.0d+0*t39*t40
            t65 = t64+t38*t63+t35*t62+t32*t61+t29*t60+t26*t59+t23*t58+t2
     1         0*t57+t17*t56+t14*t55+t11*t54+t53
            t66 = 1/t10**13
            t67 = 1/t45
            t68 = 1/t2**3
            t69 = -2.0d+0*t38*t39*t8-4.0d+0*t36*t38*t8+t64-2.2d+1*t11*t1
     1         3*t6+t53-4.0d+0*t35*t36*t37-6.0d+0*t33*t35*t37-6.0d+0*t32
     2         *t33*t34-8.0d+0*t30*t32*t34-8.0d+0*t29*t30*t31-1.0d+1*t27
     3         *t29*t31-1.0d+1*t26*t27*t28-1.2d+1*t24*t26*t28-1.2d+1*t23
     4         *t24*t25-1.4d+1*t21*t23*t25-1.4d+1*t20*t21*t22-1.6d+1*t18
     5         *t20*t22-1.6d+1*t17*t18*t19-1.8d+1*t15*t17*t19-1.8d+1*t14
     6         *t15*t16-2.0d+1*t12*t14*t16-2.0d+1*t11*t12*t13
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t4*t41*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t41*t5*t50-9.3052573
     1         63490997d-1*t4*t44*t5-1.2407009817987995d+0*t4*t41*t45)*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t41*t
     1         5*t52*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         65*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(-1.114128
     1         8576844789d+2*t43*t6*t67*t9+3.046785819198518d+4*t5*t6*t6
     2         6*t9-1.0128444160767991d+1*t38*t39*t67*t8+2.0256888321535
     3         982d+1*t36*t38*t67*t8+4.6163421503007845d+2*t35*t39*t5*t8
     4         -1.846536860120314d+3*t35*t36*t5*t8+1.3849026450902355d+3
     5         *t33*t35*t5*t8+1.1141288576844789d+2*t11*t13*t6*t67+1.012
     6         8444160767991d+1*t39*t40*t67-2.0256888321535982d+1*t35*t3
     7         6*t37*t67+3.038533248230397d+1*t33*t35*t37*t67-3.03853324
     8         8230397d+1*t32*t33*t34*t67+4.0513776643071964d+1*t30*t32*
     9         t34*t67-4.0513776643071964d+1*t29*t30*t31*t67+5.064222080
     :         383994d+1*t27*t29*t31*t67-5.064222080383994d+1*t26*t27*t2
     ;         8*t67+6.077066496460795d+1*t24*t26*t28*t67-6.077066496460
     <         795d+1*t23*t24*t25*t67+7.089910912537593d+1*t21*t23*t25*t
     =         67-7.089910912537593d+1*t20*t21*t22*t67+8.102755328614393
     >         d+1*t18*t20*t22*t67-8.102755328614393d+1*t17*t18*t19*t67+
     ?         9.115599744691191d+1*t15*t17*t19*t67-9.115599744691191d+1
     @         *t14*t15*t16*t67+1.012844416076799d+2*t12*t14*t16*t67-1.0
     1         12844416076799d+2*t11*t12*t13*t67-5.585774001863949d+4*t1
     2         3*t43*t5*t6+2.5389881826654315d+4*t11*t16*t5*t6+2.5389881
     3         826654315d+4*t12*t13*t43*t5-4.6163421503007845d+2*t38*t39
     4         *t5+4.6163421503007845d+2*t36*t38*t5+1.3849026450902355d+
     5         3*t32*t36*t37*t5-4.1547079352707067d+3*t32*t33*t37*t5+2.7
     6         69805290180471d+3*t30*t32*t37*t5+2.769805290180471d+3*t29
     7         *t33*t34*t5-7.386147440481255d+3*t29*t30*t34*t5+4.6163421
     8         50300784d+3*t27*t29*t34*t5+4.616342150300784d+3*t26*t30*t
     9         31*t5-1.154085537575196d+4*t26*t27*t31*t5+6.9245132254511
     :         78d+3*t24*t26*t31*t5+6.924513225451178d+3*t23*t27*t28*t5-
     ;         1.661883174108283d+4*t23*t24*t28*t5+9.694318515631648d+3*
     <         t21*t23*t28*t5+9.694318515631648d+3*t20*t24*t25*t5-2.2620
     =         076536473846d+4*t20*t21*t25*t5+1.2925758020842196d+4*t18*
     >         t20*t25*t5+1.2925758020842196d+4*t17*t21*t22*t5-2.9544589
     ?         76192502d+4*t17*t18*t22*t5+1.6618831741082823d+4*t15*t17*
     @         t22*t5+1.6618831741082823d+4*t14*t18*t19*t5-3.73923714174
     1         36354d+4*t14*t15*t19*t5+2.0773539676353528d+4*t12*t14*t19
     2         *t5+2.0773539676353528d+4*t11*t15*t16*t5-4.61634215030078
     3         4d+4*t11*t12*t16*t5)-9.305257363490997d-1*t41*t5*(1.03696
     4         91626335495d-6*gammaaa**3*t68/rhoa**10-3.8948928286386303
     5         d-4*t46*t47/rhoa**7.333333333333333d+0+3.531915643569816d
     6         -2*gammaaa*t3/rhoa**4.666666666666667d+0)-1.8610514726981
     7         994d+0*t44*t5*t50-2.481401963597599d+0*t41*t45*t50-2.4814
     8         01963597599d+0*t4*t44*t45-4.1356699393293317d-1*t4*t41/t4
     9         2)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t41*t5*(-3.8886
     1         343598758105d-7*t46*t68/rhoa**9-9.632497209735863d-3*t3*t
     2         49+1.29829760954621d-4*gammaaa*t47*t48)-9.305257363490997
     3         d-1*t44*t5*t52-1.2407009817987995d+0*t41*t45*t52)*wght+Cm
     4         at2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t41*t5*(1.4582378849534292d-7*gammaaa*t68/rhoa**8-3.
     2         2457440238655255d-5*t47*t51)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(4.0108638
     1         876641245d+3*t42*t6*t66*t9+6.077066496460795d+1*t35*t39*t
     2         42*t8-1.8231199489382385d+2*t33*t35*t42*t8-3.342386573053
     3         437d+3*t11*t16*t42*t6+3.342386573053437d+3*t12*t13*t42*t4
     4         3-6.077066496460795d+1*t36*t38*t42+1.8231199489382385d+2*
     5         t32*t36*t37*t42-3.646239897876477d+2*t30*t32*t37*t42+3.64
     6         6239897876477d+2*t29*t33*t34*t42-6.077066496460795d+2*t27
     7         *t29*t34*t42+6.077066496460795d+2*t26*t30*t31*t42-9.11559
     8         9744691192d+2*t24*t26*t31*t42+9.115599744691192d+2*t23*t2
     9         7*t28*t42-1.276183964256767d+3*t21*t23*t28*t42+1.27618396
     :         4256767d+3*t20*t24*t25*t42-1.7015786190090226d+3*t18*t20*
     ;         t25*t42+1.7015786190090226d+3*t17*t21*t22*t42-2.187743938
     <         725886d+3*t15*t17*t22*t42+2.187743938725886d+3*t14*t18*t1
     =         9*t42-2.734679923407357d+3*t12*t14*t19*t42+2.734679923407
     >         357d+3*t11*t15*t16*t42)-9.305257363490997d-1*t5*t50*t69-1
     ?         .2407009817987995d+0*t4*t45*t69)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t52*t65*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(5.28d+2*t6*t66*t9+t35*(1.6d+1*t36*t8+2.4d+1*t33*t
     2         8-4.0d+0*t63)+t32*(-6.0d+0*t62+3.6d+1*t33*t37+4.8d+1*t30*
     3         t37)+t29*(-8.0d+0*t61+6.4d+1*t30*t34+8.0d+1*t27*t34)+t26*
     4         (-1.0d+1*t60+1.0d+2*t27*t31+1.2d+2*t24*t31)+t11*(4.4d+2*t
     5         16*t6-2.0d+1*t55+4.0d+2*t12*t16)+t43*(4.84d+2*t13*t6-2.2d
     6         +1*t54)+t23*(-1.2d+1*t59+1.44d+2*t24*t28+1.68d+2*t21*t28)
     7         +t20*(-1.4d+1*t58+1.96d+2*t21*t25+2.24d+2*t18*t25)+t17*(-
     8         1.6d+1*t57+2.56d+2*t18*t22+2.88d+2*t15*t22)+t14*(-1.8d+1*
     9         t56+3.24d+2*t15*t19+3.6d+2*t12*t19)+t38*(8.0d+0*t39+8.0d+
     :         0*t36))*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
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
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(1)
            t7 = param(12)
            t8 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(11)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(10)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(9)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(8)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(7)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(6)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(5)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(4)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(3)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(2)
            t41 = 1/t11
            t42 = t40*t41*t9+t10*t12*t7+t6+t37*t38*t39+t34*t35*t36+t31*t
     1         32*t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t21+t1
     2         6*t17*t18+t13*t14*t15
            t43 = 1/rhob**2.6666666666666666d+0
            t44 = 4.49276922095889d-3*gammabb*t43+1.0d+0
            t45 = 1/t44
            t46 = 3.6121864536509485d-3*gammabb*t43*t45+1.0d+0
            t47 = rhob**1.3333333333333333d+0
            t48 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t49 = t48-2.0d+0*taub
            t50 = t49**11
            t51 = 2.0d+0*taub+t48
            t52 = 1/t51**11
            t53 = t49**10
            t54 = 1/t51**10
            t55 = t49**9
            t56 = 1/t51**9
            t57 = t49**8
            t58 = 1/t51**8
            t59 = t49**7
            t60 = 1/t51**7
            t61 = t49**6
            t62 = 1/t51**6
            t63 = t49**5
            t64 = 1/t51**5
            t65 = t49**4
            t66 = 1/t51**4
            t67 = t49**3
            t68 = 1/t51**3
            t69 = t49**2
            t70 = 1/t51**2
            t71 = 1/t51
            t72 = t40*t49*t71+t37*t69*t70+t50*t52*t7+t34*t67*t68+t31*t65
     1         *t66+t28*t63*t64+t25*t61*t62+t22*t59*t60+t6+t19*t57*t58+t
     2         16*t55*t56+t13*t53*t54
            t73 = rhoa**6.666666666666666d-1
            t74 = 1/t11**12
            t75 = -1.5192666241151986d+1*t39*t40*t73*t9+3.03853324823039
     1         7d+1*t37*t39*t73*t9-1.6711932865267184d+2*t10*t7*t73*t74+
     2         1.6711932865267184d+2*t12*t14*t7*t73+1.5192666241151986d+
     3         1*t40*t41*t73-3.038533248230397d+1*t36*t37*t38*t73+4.5577
     4         998723455965d+1*t34*t36*t38*t73-4.5577998723455965d+1*t33
     5         *t34*t35*t73+6.077066496460795d+1*t31*t33*t35*t73-6.07706
     6         6496460795d+1*t30*t31*t32*t73+7.596333120575993d+1*t28*t3
     7         0*t32*t73-7.596333120575993d+1*t27*t28*t29*t73+9.11559974
     8         4691193d+1*t25*t27*t29*t73-9.115599744691193d+1*t24*t25*t
     9         26*t73+1.0634866368806391d+2*t22*t24*t26*t73-1.0634866368
     :         806391d+2*t21*t22*t23*t73+1.2154132992921589d+2*t19*t21*t
     ;         23*t73-1.2154132992921589d+2*t18*t19*t20*t73+1.3673399617
     <         036788d+2*t16*t18*t20*t73-1.3673399617036788d+2*t15*t16*t
     =         17*t73+1.5192666241151986d+2*t13*t15*t17*t73-1.5192666241
     >         151986d+2*t12*t13*t14*t73
            t76 = rhoa**3.333333333333333d-1
            t77 = gammaaa**2
            t78 = 1/t2**2
            t79 = 1/rhoa**6.333333333333333d+0
            t80 = 1/rhoa**3.6666666666666664d+0
            t81 = 4.327658698487368d-5*t77*t78*t79-9.632497209735863d-3*
     1         gammaaa*t3*t80
            t82 = rhob**6.666666666666666d-1
            t83 = 1/t51**12
            t84 = -1.6711932865267184d+2*t50*t7*t82*t83+1.51926662411519
     1         86d+1*t40*t71*t82-1.5192666241151986d+1*t40*t49*t70*t82+3
     2         .038533248230397d+1*t37*t49*t70*t82+1.6711932865267184d+2
     3         *t52*t53*t7*t82-3.038533248230397d+1*t37*t68*t69*t82+4.55
     4         77998723455965d+1*t34*t68*t69*t82-4.5577998723455965d+1*t
     5         34*t66*t67*t82+6.077066496460795d+1*t31*t66*t67*t82-6.077
     6         066496460795d+1*t31*t64*t65*t82+7.596333120575993d+1*t28*
     7         t64*t65*t82-7.596333120575993d+1*t28*t62*t63*t82+9.115599
     8         744691193d+1*t25*t62*t63*t82-9.115599744691193d+1*t25*t60
     9         *t61*t82+1.0634866368806391d+2*t22*t60*t61*t82-1.06348663
     :         68806391d+2*t22*t58*t59*t82+1.2154132992921589d+2*t19*t58
     ;         *t59*t82-1.2154132992921589d+2*t19*t56*t57*t82+1.36733996
     <         17036788d+2*t16*t56*t57*t82-1.3673399617036788d+2*t16*t54
     =         *t55*t82+1.5192666241151986d+2*t13*t54*t55*t82-1.51926662
     >         41151986d+2*t13*t52*t53*t82
            t85 = rhob**3.333333333333333d-1
            t86 = gammabb**2
            t87 = 1/t44**2
            t88 = 1/rhob**6.333333333333333d+0
            t89 = 1/rhob**3.6666666666666664d+0
            t90 = 4.327658698487368d-5*t86*t87*t88-9.632497209735863d-3*
     1         gammabb*t45*t89
            t91 = 1/rhoa**5.333333333333333d+0
            t92 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t78*t91
            t93 = 1/rhob**5.333333333333333d+0
            t94 = 3.6121864536509485d-3*t43*t45-1.6228720119327628d-5*ga
     1         mmabb*t87*t93
            t95 = -2.2d+1*t10*t7*t74
            t96 = -2.2d+1*t14*t7-2.0d+1*t13*t14
            t97 = -1.8d+1*t16*t17-2.0d+1*t13*t17
            t98 = -1.6d+1*t19*t20-1.8d+1*t16*t20
            t99 = -1.4d+1*t22*t23-1.6d+1*t19*t23
            t100 = -1.2d+1*t25*t26-1.4d+1*t22*t26
            t101 = -1.0d+1*t28*t29-1.2d+1*t25*t29
            t102 = -8.0d+0*t31*t32-1.0d+1*t28*t32
            t103 = -6.0d+0*t34*t35-8.0d+0*t31*t35
            t104 = -4.0d+0*t37*t38-6.0d+0*t34*t38
            t105 = -2.0d+0*t40*t9-4.0d+0*t37*t9
            t106 = -2.0d+0*t40*t41
            t107 = t21*t99+t18*t98+t15*t97+t12*t96+t95+t105*t39+t104*t36
     1         +t103*t33+t102*t30+t101*t27+t100*t24+t106
            t108 = -2.2d+1*t50*t7*t83
            t109 = -2.2d+1*t53*t7-2.0d+1*t13*t53
            t110 = -1.8d+1*t16*t55-2.0d+1*t13*t55
            t111 = -1.6d+1*t19*t57-1.8d+1*t16*t57
            t112 = -1.4d+1*t22*t59-1.6d+1*t19*t59
            t113 = -1.2d+1*t25*t61-1.4d+1*t22*t61
            t114 = -1.0d+1*t28*t63-1.2d+1*t25*t63
            t115 = -8.0d+0*t31*t65-1.0d+1*t28*t65
            t116 = -6.0d+0*t34*t67-8.0d+0*t31*t67
            t117 = -4.0d+0*t37*t69-6.0d+0*t34*t69
            t118 = -2.0d+0*t40*t49-4.0d+0*t37*t49
            t119 = -2.0d+0*t40*t71
            t120 = t118*t70+t117*t68+t116*t66+t115*t64+t114*t62+t113*t60
     1         +t112*t58+t111*t56+t110*t54+t109*t52+t119+t108
            t121 = 1/t11**13
            t122 = 1/t76
            t123 = 1/t2**3
            t124 = 1/t51**13
            t125 = 1/t85
            t126 = 1/t44**3
            t127 = t95-2.0d+0*t39*t40*t9-4.0d+0*t37*t39*t9-2.2d+1*t12*t1
     1         4*t7-4.0d+0*t36*t37*t38-6.0d+0*t34*t36*t38-6.0d+0*t33*t34
     2         *t35-8.0d+0*t31*t33*t35-8.0d+0*t30*t31*t32-1.0d+1*t28*t30
     3         *t32-1.0d+1*t27*t28*t29-1.2d+1*t25*t27*t29-1.2d+1*t24*t25
     4         *t26-1.4d+1*t22*t24*t26-1.4d+1*t21*t22*t23-1.6d+1*t19*t21
     5         *t23-1.6d+1*t18*t19*t20-1.8d+1*t16*t18*t20-1.8d+1*t15*t16
     6         *t17-2.0d+1*t13*t15*t17-2.0d+1*t12*t13*t14+t106
            t128 = -2.0d+0*t40*t49*t70-4.0d+0*t37*t49*t70-2.2d+1*t52*t53
     1         *t7-4.0d+0*t37*t68*t69-6.0d+0*t34*t68*t69-6.0d+0*t34*t66*
     2         t67-8.0d+0*t31*t66*t67-8.0d+0*t31*t64*t65-1.0d+1*t28*t64*
     3         t65-1.0d+1*t28*t62*t63-1.2d+1*t25*t62*t63-1.2d+1*t25*t60*
     4         t61-1.4d+1*t22*t60*t61-1.4d+1*t22*t58*t59-1.6d+1*t19*t58*
     5         t59-1.6d+1*t19*t56*t57-1.8d+1*t16*t56*t57-1.8d+1*t16*t54*
     6         t55-2.0d+1*t13*t54*t55-2.0d+1*t13*t52*t53+t119+t108
            t129 = 8.0d+0*t40+8.0d+0*t37
            fnc(iq) = (-9.305257363490997d-1*t46*t47*t72-9.3052573634909
     1         97d-1*t4*t42*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t42*t5*t81-1.2407009
     1         817987995d+0*t4*t42*t76-9.305257363490997d-1*t4*t5*t75)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t47*t72*t90-1.240700
     1         9817987995d+0*t46*t72*t85-9.305257363490997d-1*t46*t47*t8
     2         4)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t42*t
     1         5*t92*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t47*t
     1         72*t94*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t107*t4
     1         *t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t120*t4
     1         6*t47*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(4.6163421
     1         503007845d+2*t36*t40*t5*t9-1.846536860120314d+3*t36*t37*t
     2         5*t9+1.3849026450902355d+3*t34*t36*t5*t9-1.01284441607679
     3         91d+1*t122*t39*t40*t9+2.0256888321535982d+1*t122*t37*t39*
     4         t9-5.585774001863949d+4*t14*t5*t7*t74-1.1141288576844789d
     5         +2*t10*t122*t7*t74+2.5389881826654315d+4*t13*t14*t5*t74+2
     6         .5389881826654315d+4*t12*t17*t5*t7+3.046785819198518d+4*t
     7         10*t121*t5*t7+1.1141288576844789d+2*t12*t122*t14*t7-4.616
     8         3421503007845d+2*t39*t40*t5+4.6163421503007845d+2*t37*t39
     9         *t5+1.3849026450902355d+3*t33*t37*t38*t5-4.15470793527070
     :         67d+3*t33*t34*t38*t5+2.769805290180471d+3*t31*t33*t38*t5+
     ;         2.769805290180471d+3*t30*t34*t35*t5-7.386147440481255d+3*
     <         t30*t31*t35*t5+4.616342150300784d+3*t28*t30*t35*t5+4.6163
     =         42150300784d+3*t27*t31*t32*t5-1.154085537575196d+4*t27*t2
     >         8*t32*t5+6.924513225451178d+3*t25*t27*t32*t5+6.9245132254
     ?         51178d+3*t24*t28*t29*t5-1.661883174108283d+4*t24*t25*t29*
     @         t5+9.694318515631648d+3*t22*t24*t29*t5+9.694318515631648d
     1         +3*t21*t25*t26*t5-2.2620076536473846d+4*t21*t22*t26*t5+1.
     2         2925758020842196d+4*t19*t21*t26*t5+1.2925758020842196d+4*
     3         t18*t22*t23*t5-2.954458976192502d+4*t18*t19*t23*t5+1.6618
     4         831741082823d+4*t16*t18*t23*t5+1.6618831741082823d+4*t15*
     5         t19*t20*t5-3.7392371417436354d+4*t15*t16*t20*t5+2.0773539
     6         676353528d+4*t13*t15*t20*t5+2.0773539676353528d+4*t12*t16
     7         *t17*t5-4.616342150300784d+4*t12*t13*t17*t5+1.01284441607
     8         67991d+1*t122*t40*t41-2.0256888321535982d+1*t122*t36*t37*
     9         t38+3.038533248230397d+1*t122*t34*t36*t38-3.0385332482303
     :         97d+1*t122*t33*t34*t35+4.0513776643071964d+1*t122*t31*t33
     ;         *t35-4.0513776643071964d+1*t122*t30*t31*t32+5.06422208038
     <         3994d+1*t122*t28*t30*t32-5.064222080383994d+1*t122*t27*t2
     =         8*t29+6.077066496460795d+1*t122*t25*t27*t29-6.07706649646
     >         0795d+1*t122*t24*t25*t26+7.089910912537593d+1*t122*t22*t2
     ?         4*t26-7.089910912537593d+1*t122*t21*t22*t23+8.10275532861
     @         4393d+1*t122*t19*t21*t23-8.102755328614393d+1*t122*t18*t1
     1         9*t20+9.115599744691191d+1*t122*t16*t18*t20-9.11559974469
     2         1191d+1*t122*t15*t16*t17+1.012844416076799d+2*t122*t13*t1
     3         5*t17-1.012844416076799d+2*t12*t122*t13*t14)-2.4814019635
     4         97599d+0*t42*t76*t81-1.8610514726981994d+0*t5*t75*t81-9.3
     5         05257363490997d-1*t42*t5*(-3.8948928286386303d-4*t77*t78/
     6         rhoa**7.333333333333333d+0+3.531915643569816d-2*gammaaa*t
     7         3/rhoa**4.666666666666667d+0+1.0369691626335495d-6*gammaa
     8         a**3*t123/rhoa**10)-2.481401963597599d+0*t4*t75*t76-4.135
     9         6699393293317d-1*t4*t42/t73)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-2.481401963597599d+0*t72*t85*t90-1.86
     1         10514726981994d+0*t47*t84*t90-9.305257363490997d-1*t47*t7
     2         2*(-3.8948928286386303d-4*t86*t87/rhob**7.333333333333333
     3         d+0+3.531915643569816d-2*gammabb*t45/rhob**4.666666666666
     4         667d+0+1.0369691626335495d-6*gammabb**3*t126/rhob**10)-2.
     5         481401963597599d+0*t46*t84*t85-9.305257363490997d-1*t46*t
     6         47*(-5.585774001863949d+4*t47*t53*t7*t83-1.11412885768447
     7         89d+2*t125*t50*t7*t83+2.5389881826654315d+4*t13*t47*t53*t
     8         83+1.0128444160767991d+1*t125*t40*t71-1.0128444160767991d
     9         +1*t125*t40*t49*t70+2.0256888321535982d+1*t125*t37*t49*t7
     :         0-4.6163421503007845d+2*t40*t47*t70+4.6163421503007845d+2
     ;         *t37*t47*t70+2.5389881826654315d+4*t47*t52*t55*t7+1.11412
     <         88576844789d+2*t125*t52*t53*t7+3.046785819198518d+4*t124*
     =         t47*t50*t7-2.0256888321535982d+1*t125*t37*t68*t69+3.03853
     >         3248230397d+1*t125*t34*t68*t69+1.3849026450902355d+3*t37*
     ?         t47*t66*t69-4.1547079352707067d+3*t34*t47*t66*t69+2.76980
     @         5290180471d+3*t31*t47*t66*t69+4.6163421503007845d+2*t40*t
     1         47*t49*t68-1.846536860120314d+3*t37*t47*t49*t68+1.3849026
     2         450902355d+3*t34*t47*t49*t68-3.038533248230397d+1*t125*t3
     3         4*t66*t67+4.0513776643071964d+1*t125*t31*t66*t67+2.769805
     4         290180471d+3*t34*t47*t64*t67-7.386147440481255d+3*t31*t47
     5         *t64*t67+4.616342150300784d+3*t28*t47*t64*t67-4.051377664
     6         3071964d+1*t125*t31*t64*t65+5.064222080383994d+1*t125*t28
     7         *t64*t65+4.616342150300784d+3*t31*t47*t62*t65-1.154085537
     8         575196d+4*t28*t47*t62*t65+6.924513225451178d+3*t25*t47*t6
     9         2*t65-5.064222080383994d+1*t125*t28*t62*t63+6.07706649646
     :         0795d+1*t125*t25*t62*t63+6.924513225451178d+3*t28*t47*t60
     ;         *t63-1.661883174108283d+4*t25*t47*t60*t63+9.6943185156316
     <         48d+3*t22*t47*t60*t63-6.077066496460795d+1*t125*t25*t60*t
     =         61+7.089910912537593d+1*t125*t22*t60*t61+9.69431851563164
     >         8d+3*t25*t47*t58*t61-2.2620076536473846d+4*t22*t47*t58*t6
     ?         1+1.2925758020842196d+4*t19*t47*t58*t61-7.089910912537593
     @         d+1*t125*t22*t58*t59+8.102755328614393d+1*t125*t19*t58*t5
     1         9+1.2925758020842196d+4*t22*t47*t56*t59-2.954458976192502
     2         d+4*t19*t47*t56*t59+1.6618831741082823d+4*t16*t47*t56*t59
     3         -8.102755328614393d+1*t125*t19*t56*t57+9.115599744691191d
     4         +1*t125*t16*t56*t57+1.6618831741082823d+4*t19*t47*t54*t57
     5         -3.7392371417436354d+4*t16*t47*t54*t57+2.0773539676353528
     6         d+4*t13*t47*t54*t57-9.115599744691191d+1*t125*t16*t54*t55
     7         +1.012844416076799d+2*t125*t13*t54*t55+2.0773539676353528
     8         d+4*t16*t47*t52*t55-4.616342150300784d+4*t13*t47*t52*t55-
     9         1.012844416076799d+2*t125*t13*t52*t53)-4.1356699393293317
     :         d-1*t46*t72/t82)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817987995d+0*t42*t76*t92-9.
     1         305257363490997d-1*t5*t75*t92-9.305257363490997d-1*t42*t5
     2         *(-9.632497209735863d-3*t3*t80+1.29829760954621d-4*gammaa
     3         a*t78*t79-3.8886343598758105d-7*t123*t77/rhoa**9))*wght+C
     4         mat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.2407009817987995d+0*t72*t85*t94-9.
     1         305257363490997d-1*t47*t84*t94-9.305257363490997d-1*t47*t
     2         72*(-9.632497209735863d-3*t45*t89+1.29829760954621d-4*gam
     3         mabb*t87*t88-3.8886343598758105d-7*t126*t86/rhob**9))*wgh
     4         t+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t42*t5*(1.4582378849534292d-7*gammaaa*t123/rhoa**8-3
     2         .2457440238655255d-5*t78*t91)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t47*t72*(1.4582378849534292d-7*gammabb*t126/rhob**8-
     2         3.2457440238655255d-5*t87*t93)*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(6.0770664
     1         96460795d+1*t36*t40*t73*t9-1.8231199489382385d+2*t34*t36*
     2         t73*t9+3.342386573053437d+3*t13*t14*t73*t74-3.34238657305
     3         3437d+3*t12*t17*t7*t73+4.0108638876641245d+3*t10*t121*t7*
     4         t73-6.077066496460795d+1*t37*t39*t73+1.8231199489382385d+
     5         2*t33*t37*t38*t73-3.646239897876477d+2*t31*t33*t38*t73+3.
     6         646239897876477d+2*t30*t34*t35*t73-6.077066496460795d+2*t
     7         28*t30*t35*t73+6.077066496460795d+2*t27*t31*t32*t73-9.115
     8         599744691192d+2*t25*t27*t32*t73+9.115599744691192d+2*t24*
     9         t28*t29*t73-1.276183964256767d+3*t22*t24*t29*t73+1.276183
     :         964256767d+3*t21*t25*t26*t73-1.7015786190090226d+3*t19*t2
     ;         1*t26*t73+1.7015786190090226d+3*t18*t22*t23*t73-2.1877439
     <         38725886d+3*t16*t18*t23*t73+2.187743938725886d+3*t15*t19*
     =         t20*t73-2.734679923407357d+3*t13*t15*t20*t73+2.7346799234
     >         07357d+3*t12*t16*t17*t73)-9.305257363490997d-1*t127*t5*t8
     ?         1-1.2407009817987995d+0*t127*t4*t76)*wght+Mmat2(iq,D2_RA_
     @         TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t128*t47*t90-1.2
     1         407009817987995d+0*t128*t46*t85-9.305257363490997d-1*t46*
     2         t47*(3.342386573053437d+3*t13*t53*t82*t83-6.0770664964607
     3         95d+1*t37*t70*t82-3.342386573053437d+3*t52*t55*t7*t82+4.0
     4         108638876641245d+3*t124*t50*t7*t82+1.8231199489382385d+2*
     5         t37*t66*t69*t82-3.646239897876477d+2*t31*t66*t69*t82+6.07
     6         7066496460795d+1*t40*t49*t68*t82-1.8231199489382385d+2*t3
     7         4*t49*t68*t82+3.646239897876477d+2*t34*t64*t67*t82-6.0770
     8         66496460795d+2*t28*t64*t67*t82+6.077066496460795d+2*t31*t
     9         62*t65*t82-9.115599744691192d+2*t25*t62*t65*t82+9.1155997
     :         44691192d+2*t28*t60*t63*t82-1.276183964256767d+3*t22*t60*
     ;         t63*t82+1.276183964256767d+3*t25*t58*t61*t82-1.7015786190
     <         090226d+3*t19*t58*t61*t82+1.7015786190090226d+3*t22*t56*t
     =         59*t82-2.187743938725886d+3*t16*t56*t59*t82+2.18774393872
     >         5886d+3*t19*t54*t57*t82-2.734679923407357d+3*t13*t54*t57*
     ?         t82+2.734679923407357d+3*t16*t52*t55*t82))*wght+Mmat2(iq,
     @         D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t107*t5*t92*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t120*t47*t94*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(t18*(-1.6d+1*t99+2.56d+2*t19*t23+2.88d+2*t16*t23)
     2         +t15*(-1.8d+1*t98+3.24d+2*t16*t20+3.6d+2*t13*t20)+t12*(-2
     3         .0d+1*t97+4.4d+2*t17*t7+4.0d+2*t13*t17)+t74*(4.84d+2*t14*
     4         t7-2.2d+1*t96)+t36*(1.6d+1*t37*t9+2.4d+1*t34*t9-4.0d+0*t1
     5         05)+5.28d+2*t10*t121*t7+t129*t39+t33*(3.6d+1*t34*t38+4.8d
     6         +1*t31*t38-6.0d+0*t104)+t30*(6.4d+1*t31*t35+8.0d+1*t28*t3
     7         5-8.0d+0*t103)+t27*(1.0d+2*t28*t32+1.2d+2*t25*t32-1.0d+1*
     8         t102)+t24*(1.44d+2*t25*t29+1.68d+2*t22*t29-1.2d+1*t101)+t
     9         21*(1.96d+2*t22*t26+2.24d+2*t19*t26-1.4d+1*t100))*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t46*t47*((4.84d+2*t53*t7-2.2d+1*t109)*t83+t129*t70+t52*(
     2         4.4d+2*t55*t7+4.0d+2*t13*t55-2.0d+1*t110)+5.28d+2*t124*t5
     3         0*t7+t66*(3.6d+1*t34*t69+4.8d+1*t31*t69-6.0d+0*t117)+(1.6
     4         d+1*t37*t49+2.4d+1*t34*t49-4.0d+0*t118)*t68+t64*(6.4d+1*t
     5         31*t67+8.0d+1*t28*t67-8.0d+0*t116)+t62*(1.0d+2*t28*t65+1.
     6         2d+2*t25*t65-1.0d+1*t115)+t60*(1.44d+2*t25*t63+1.68d+2*t2
     7         2*t63-1.2d+1*t114)+t58*(1.96d+2*t22*t61+2.24d+2*t19*t61-1
     8         .4d+1*t113)+t56*(2.56d+2*t19*t59+2.88d+2*t16*t59-1.6d+1*t
     9         112)+t54*(3.24d+2*t16*t57+3.6d+2*t13*t57-1.8d+1*t111))*wg
     :         ht
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(11)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(10)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(9)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(8)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(7)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(6)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(5)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(4)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(3)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(2)
            t40 = 1/t10
            t41 = t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t30*t31*t
     1         32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t20+t15*t1
     2         6*t17+t12*t13*t14+param(1)
            t42 = rhoa**6.666666666666666d-1
            t43 = 1/t10**12
            t44 = -1.6711932865267184d+2*t42*t43*t6*t9-1.519266624115198
     1         6d+1*t38*t39*t42*t8+3.038533248230397d+1*t36*t38*t42*t8+1
     2         .6711932865267184d+2*t11*t13*t42*t6+1.5192666241151986d+1
     3         *t39*t40*t42-3.038533248230397d+1*t35*t36*t37*t42+4.55779
     4         98723455965d+1*t33*t35*t37*t42-4.5577998723455965d+1*t32*
     5         t33*t34*t42+6.077066496460795d+1*t30*t32*t34*t42-6.077066
     6         496460795d+1*t29*t30*t31*t42+7.596333120575993d+1*t27*t29
     7         *t31*t42-7.596333120575993d+1*t26*t27*t28*t42+9.115599744
     8         691193d+1*t24*t26*t28*t42-9.115599744691193d+1*t23*t24*t2
     9         5*t42+1.0634866368806391d+2*t21*t23*t25*t42-1.06348663688
     :         06391d+2*t20*t21*t22*t42+1.2154132992921589d+2*t18*t20*t2
     ;         2*t42-1.2154132992921589d+2*t17*t18*t19*t42+1.36733996170
     <         36788d+2*t15*t17*t19*t42-1.3673399617036788d+2*t14*t15*t1
     =         6*t42+1.5192666241151986d+2*t12*t14*t16*t42-1.51926662411
     >         51986d+2*t11*t12*t13*t42
            t45 = rhoa**3.333333333333333d-1
            t46 = gammaaa**2
            t47 = 1/t2**2
            t48 = 1/rhoa**6.333333333333333d+0
            t49 = 1/rhoa**3.6666666666666664d+0
            t50 = 4.327658698487368d-5*t46*t47*t48-9.632497209735863d-3*
     1         gammaaa*t3*t49
            t51 = 1/rhoa**5.333333333333333d+0
            t52 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t47*t51
            t53 = -2.2d+1*t43*t6*t9
            t54 = -2.2d+1*t13*t6-2.0d+1*t12*t13
            t55 = -1.8d+1*t15*t16-2.0d+1*t12*t16
            t56 = -1.6d+1*t18*t19-1.8d+1*t15*t19
            t57 = -1.4d+1*t21*t22-1.6d+1*t18*t22
            t58 = -1.2d+1*t24*t25-1.4d+1*t21*t25
            t59 = -1.0d+1*t27*t28-1.2d+1*t24*t28
            t60 = -8.0d+0*t30*t31-1.0d+1*t27*t31
            t61 = -6.0d+0*t33*t34-8.0d+0*t30*t34
            t62 = -4.0d+0*t36*t37-6.0d+0*t33*t37
            t63 = -2.0d+0*t39*t8-4.0d+0*t36*t8
            t64 = -2.0d+0*t39*t40
            t65 = t64+t38*t63+t35*t62+t32*t61+t29*t60+t26*t59+t23*t58+t2
     1         0*t57+t17*t56+t14*t55+t11*t54+t53
            t66 = 1/t10**13
            t67 = 1/t45
            t68 = 1/t2**3
            t69 = -2.0d+0*t38*t39*t8-4.0d+0*t36*t38*t8+t64-2.2d+1*t11*t1
     1         3*t6+t53-4.0d+0*t35*t36*t37-6.0d+0*t33*t35*t37-6.0d+0*t32
     2         *t33*t34-8.0d+0*t30*t32*t34-8.0d+0*t29*t30*t31-1.0d+1*t27
     3         *t29*t31-1.0d+1*t26*t27*t28-1.2d+1*t24*t26*t28-1.2d+1*t23
     4         *t24*t25-1.4d+1*t21*t23*t25-1.4d+1*t20*t21*t22-1.6d+1*t18
     5         *t20*t22-1.6d+1*t17*t18*t19-1.8d+1*t15*t17*t19-1.8d+1*t14
     6         *t15*t16-2.0d+1*t12*t14*t16-2.0d+1*t11*t12*t13
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t41*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t41*t5*t50*wght-9.305
     1         257363490997d-1*t4*t44*t5*wght-1.2407009817987995d+0*t4*t
     2         41*t45*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t41*t
     1         5*t52*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         65*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = -9.305257363490997d-1*t4*t5*(-1.1141288
     1         576844789d+2*t43*t6*t67*t9+3.046785819198518d+4*t5*t6*t66
     2         *t9-1.0128444160767991d+1*t38*t39*t67*t8+2.02568883215359
     3         82d+1*t36*t38*t67*t8+4.6163421503007845d+2*t35*t39*t5*t8-
     4         1.846536860120314d+3*t35*t36*t5*t8+1.3849026450902355d+3*
     5         t33*t35*t5*t8+1.1141288576844789d+2*t11*t13*t6*t67+1.0128
     6         444160767991d+1*t39*t40*t67-2.0256888321535982d+1*t35*t36
     7         *t37*t67+3.038533248230397d+1*t33*t35*t37*t67-3.038533248
     8         230397d+1*t32*t33*t34*t67+4.0513776643071964d+1*t30*t32*t
     9         34*t67-4.0513776643071964d+1*t29*t30*t31*t67+5.0642220803
     :         83994d+1*t27*t29*t31*t67-5.064222080383994d+1*t26*t27*t28
     ;         *t67+6.077066496460795d+1*t24*t26*t28*t67-6.0770664964607
     <         95d+1*t23*t24*t25*t67+7.089910912537593d+1*t21*t23*t25*t6
     =         7-7.089910912537593d+1*t20*t21*t22*t67+8.102755328614393d
     >         +1*t18*t20*t22*t67-8.102755328614393d+1*t17*t18*t19*t67+9
     ?         .115599744691191d+1*t15*t17*t19*t67-9.115599744691191d+1*
     @         t14*t15*t16*t67+1.012844416076799d+2*t12*t14*t16*t67-1.01
     1         2844416076799d+2*t11*t12*t13*t67-5.585774001863949d+4*t13
     2         *t43*t5*t6+2.5389881826654315d+4*t11*t16*t5*t6+2.53898818
     3         26654315d+4*t12*t13*t43*t5-4.6163421503007845d+2*t38*t39*
     4         t5+4.6163421503007845d+2*t36*t38*t5+1.3849026450902355d+3
     5         *t32*t36*t37*t5-4.1547079352707067d+3*t32*t33*t37*t5+2.76
     6         9805290180471d+3*t30*t32*t37*t5+2.769805290180471d+3*t29*
     7         t33*t34*t5-7.386147440481255d+3*t29*t30*t34*t5+4.61634215
     8         0300784d+3*t27*t29*t34*t5+4.616342150300784d+3*t26*t30*t3
     9         1*t5-1.154085537575196d+4*t26*t27*t31*t5+6.92451322545117
     :         8d+3*t24*t26*t31*t5+6.924513225451178d+3*t23*t27*t28*t5-1
     ;         .661883174108283d+4*t23*t24*t28*t5+9.694318515631648d+3*t
     <         21*t23*t28*t5+9.694318515631648d+3*t20*t24*t25*t5-2.26200
     =         76536473846d+4*t20*t21*t25*t5+1.2925758020842196d+4*t18*t
     >         20*t25*t5+1.2925758020842196d+4*t17*t21*t22*t5-2.95445897
     ?         6192502d+4*t17*t18*t22*t5+1.6618831741082823d+4*t15*t17*t
     @         22*t5+1.6618831741082823d+4*t14*t18*t19*t5-3.739237141743
     1         6354d+4*t14*t15*t19*t5+2.0773539676353528d+4*t12*t14*t19*
     2         t5+2.0773539676353528d+4*t11*t15*t16*t5-4.616342150300784
     3         d+4*t11*t12*t16*t5)*wght-9.305257363490997d-1*t41*t5*(1.0
     4         369691626335495d-6*gammaaa**3*t68/rhoa**10-3.894892828638
     5         6303d-4*t46*t47/rhoa**7.333333333333333d+0+3.531915643569
     6         816d-2*gammaaa*t3/rhoa**4.666666666666667d+0)*wght-1.8610
     7         514726981994d+0*t44*t5*t50*wght-2.481401963597599d+0*t41*
     8         t45*t50*wght-2.481401963597599d+0*t4*t44*t45*wght-4.13566
     9         99393293317d-1*t4*t41*wght/t42+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490997d-1*t41*t5*(-3.88863
     1         43598758105d-7*t46*t68/rhoa**9-9.632497209735863d-3*t3*t4
     2         9+1.29829760954621d-4*gammaaa*t47*t48)*wght-9.30525736349
     3         0997d-1*t44*t5*t52*wght-1.2407009817987995d+0*t41*t45*t52
     4         *wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t41*t5*(1.4582378849534292d-7*gammaaa*t68/rhoa**8-3.
     2         2457440238655255d-5*t47*t51)*wght
            Mmat2(iq,D2_RA_TA) = -9.305257363490997d-1*t4*t5*(4.01086388
     1         76641245d+3*t42*t6*t66*t9+6.077066496460795d+1*t35*t39*t4
     2         2*t8-1.8231199489382385d+2*t33*t35*t42*t8-3.3423865730534
     3         37d+3*t11*t16*t42*t6+3.342386573053437d+3*t12*t13*t42*t43
     4         -6.077066496460795d+1*t36*t38*t42+1.8231199489382385d+2*t
     5         32*t36*t37*t42-3.646239897876477d+2*t30*t32*t37*t42+3.646
     6         239897876477d+2*t29*t33*t34*t42-6.077066496460795d+2*t27*
     7         t29*t34*t42+6.077066496460795d+2*t26*t30*t31*t42-9.115599
     8         744691192d+2*t24*t26*t31*t42+9.115599744691192d+2*t23*t27
     9         *t28*t42-1.276183964256767d+3*t21*t23*t28*t42+1.276183964
     :         256767d+3*t20*t24*t25*t42-1.7015786190090226d+3*t18*t20*t
     ;         25*t42+1.7015786190090226d+3*t17*t21*t22*t42-2.1877439387
     <         25886d+3*t15*t17*t22*t42+2.187743938725886d+3*t14*t18*t19
     =         *t42-2.734679923407357d+3*t12*t14*t19*t42+2.7346799234073
     >         57d+3*t11*t15*t16*t42)*wght-9.305257363490997d-1*t5*t50*t
     ?         69*wght-1.2407009817987995d+0*t4*t45*t69*wght+Mmat2(iq,D2
     @         _RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t52*t65*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(5.28d+2*t6*t66*t9+t35*(1.6d+1*t36*t8+2.4d+1*t33*t
     2         8-4.0d+0*t63)+t32*(-6.0d+0*t62+3.6d+1*t33*t37+4.8d+1*t30*
     3         t37)+t29*(-8.0d+0*t61+6.4d+1*t30*t34+8.0d+1*t27*t34)+t26*
     4         (-1.0d+1*t60+1.0d+2*t27*t31+1.2d+2*t24*t31)+t11*(4.4d+2*t
     5         16*t6-2.0d+1*t55+4.0d+2*t12*t16)+t43*(4.84d+2*t13*t6-2.2d
     6         +1*t54)+t23*(-1.2d+1*t59+1.44d+2*t24*t28+1.68d+2*t21*t28)
     7         +t20*(-1.4d+1*t58+1.96d+2*t21*t25+2.24d+2*t18*t25)+t17*(-
     8         1.6d+1*t57+2.56d+2*t18*t22+2.88d+2*t15*t22)+t14*(-1.8d+1*
     9         t56+3.24d+2*t15*t19+3.6d+2*t12*t19)+t38*(8.0d+0*t39+8.0d+
     :         0*t36))*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = param(12)
            t7 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t8 = t7-2.0d+0*taub
            t9 = t8**11
            t10 = 2.0d+0*taub+t7
            t11 = 1/t10**11
            t12 = param(11)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(10)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(9)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(8)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(7)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(6)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(5)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(4)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(3)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(2)
            t40 = 1/t10
            t41 = t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t30*t31*t
     1         32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t20+t15*t1
     2         6*t17+t12*t13*t14+param(1)
            t42 = rhob**6.666666666666666d-1
            t43 = 1/t10**12
            t44 = -1.6711932865267184d+2*t42*t43*t6*t9-1.519266624115198
     1         6d+1*t38*t39*t42*t8+3.038533248230397d+1*t36*t38*t42*t8+1
     2         .6711932865267184d+2*t11*t13*t42*t6+1.5192666241151986d+1
     3         *t39*t40*t42-3.038533248230397d+1*t35*t36*t37*t42+4.55779
     4         98723455965d+1*t33*t35*t37*t42-4.5577998723455965d+1*t32*
     5         t33*t34*t42+6.077066496460795d+1*t30*t32*t34*t42-6.077066
     6         496460795d+1*t29*t30*t31*t42+7.596333120575993d+1*t27*t29
     7         *t31*t42-7.596333120575993d+1*t26*t27*t28*t42+9.115599744
     8         691193d+1*t24*t26*t28*t42-9.115599744691193d+1*t23*t24*t2
     9         5*t42+1.0634866368806391d+2*t21*t23*t25*t42-1.06348663688
     :         06391d+2*t20*t21*t22*t42+1.2154132992921589d+2*t18*t20*t2
     ;         2*t42-1.2154132992921589d+2*t17*t18*t19*t42+1.36733996170
     <         36788d+2*t15*t17*t19*t42-1.3673399617036788d+2*t14*t15*t1
     =         6*t42+1.5192666241151986d+2*t12*t14*t16*t42-1.51926662411
     >         51986d+2*t11*t12*t13*t42
            t45 = rhob**3.333333333333333d-1
            t46 = gammabb**2
            t47 = 1/t2**2
            t48 = 1/rhob**6.333333333333333d+0
            t49 = 1/rhob**3.6666666666666664d+0
            t50 = 4.327658698487368d-5*t46*t47*t48-9.632497209735863d-3*
     1         gammabb*t3*t49
            t51 = 1/rhob**5.333333333333333d+0
            t52 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         abb*t47*t51
            t53 = -2.2d+1*t43*t6*t9
            t54 = -2.2d+1*t13*t6-2.0d+1*t12*t13
            t55 = -1.8d+1*t15*t16-2.0d+1*t12*t16
            t56 = -1.6d+1*t18*t19-1.8d+1*t15*t19
            t57 = -1.4d+1*t21*t22-1.6d+1*t18*t22
            t58 = -1.2d+1*t24*t25-1.4d+1*t21*t25
            t59 = -1.0d+1*t27*t28-1.2d+1*t24*t28
            t60 = -8.0d+0*t30*t31-1.0d+1*t27*t31
            t61 = -6.0d+0*t33*t34-8.0d+0*t30*t34
            t62 = -4.0d+0*t36*t37-6.0d+0*t33*t37
            t63 = -2.0d+0*t39*t8-4.0d+0*t36*t8
            t64 = -2.0d+0*t39*t40
            t65 = t64+t38*t63+t35*t62+t32*t61+t29*t60+t26*t59+t23*t58+t2
     1         0*t57+t17*t56+t14*t55+t11*t54+t53
            t66 = 1/t10**13
            t67 = 1/t45
            t68 = 1/t2**3
            t69 = -2.0d+0*t38*t39*t8-4.0d+0*t36*t38*t8+t64-2.2d+1*t11*t1
     1         3*t6+t53-4.0d+0*t35*t36*t37-6.0d+0*t33*t35*t37-6.0d+0*t32
     2         *t33*t34-8.0d+0*t30*t32*t34-8.0d+0*t29*t30*t31-1.0d+1*t27
     3         *t29*t31-1.0d+1*t26*t27*t28-1.2d+1*t24*t26*t28-1.2d+1*t23
     4         *t24*t25-1.4d+1*t21*t23*t25-1.4d+1*t20*t21*t22-1.6d+1*t18
     5         *t20*t22-1.6d+1*t17*t18*t19-1.8d+1*t15*t17*t19-1.8d+1*t14
     6         *t15*t16-2.0d+1*t12*t14*t16-2.0d+1*t11*t12*t13
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t41*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t41*t5*t50*wght-9.305
     1         257363490997d-1*t4*t44*t5*wght-1.2407009817987995d+0*t4*t
     2         41*t45*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t41*t
     1         5*t52*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t4*t5*t
     1         65*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = -9.305257363490997d-1*t4*t5*(-1.1141288
     1         576844789d+2*t43*t6*t67*t9+3.046785819198518d+4*t5*t6*t66
     2         *t9-1.0128444160767991d+1*t38*t39*t67*t8+2.02568883215359
     3         82d+1*t36*t38*t67*t8+4.6163421503007845d+2*t35*t39*t5*t8-
     4         1.846536860120314d+3*t35*t36*t5*t8+1.3849026450902355d+3*
     5         t33*t35*t5*t8+1.1141288576844789d+2*t11*t13*t6*t67+1.0128
     6         444160767991d+1*t39*t40*t67-2.0256888321535982d+1*t35*t36
     7         *t37*t67+3.038533248230397d+1*t33*t35*t37*t67-3.038533248
     8         230397d+1*t32*t33*t34*t67+4.0513776643071964d+1*t30*t32*t
     9         34*t67-4.0513776643071964d+1*t29*t30*t31*t67+5.0642220803
     :         83994d+1*t27*t29*t31*t67-5.064222080383994d+1*t26*t27*t28
     ;         *t67+6.077066496460795d+1*t24*t26*t28*t67-6.0770664964607
     <         95d+1*t23*t24*t25*t67+7.089910912537593d+1*t21*t23*t25*t6
     =         7-7.089910912537593d+1*t20*t21*t22*t67+8.102755328614393d
     >         +1*t18*t20*t22*t67-8.102755328614393d+1*t17*t18*t19*t67+9
     ?         .115599744691191d+1*t15*t17*t19*t67-9.115599744691191d+1*
     @         t14*t15*t16*t67+1.012844416076799d+2*t12*t14*t16*t67-1.01
     1         2844416076799d+2*t11*t12*t13*t67-5.585774001863949d+4*t13
     2         *t43*t5*t6+2.5389881826654315d+4*t11*t16*t5*t6+2.53898818
     3         26654315d+4*t12*t13*t43*t5-4.6163421503007845d+2*t38*t39*
     4         t5+4.6163421503007845d+2*t36*t38*t5+1.3849026450902355d+3
     5         *t32*t36*t37*t5-4.1547079352707067d+3*t32*t33*t37*t5+2.76
     6         9805290180471d+3*t30*t32*t37*t5+2.769805290180471d+3*t29*
     7         t33*t34*t5-7.386147440481255d+3*t29*t30*t34*t5+4.61634215
     8         0300784d+3*t27*t29*t34*t5+4.616342150300784d+3*t26*t30*t3
     9         1*t5-1.154085537575196d+4*t26*t27*t31*t5+6.92451322545117
     :         8d+3*t24*t26*t31*t5+6.924513225451178d+3*t23*t27*t28*t5-1
     ;         .661883174108283d+4*t23*t24*t28*t5+9.694318515631648d+3*t
     <         21*t23*t28*t5+9.694318515631648d+3*t20*t24*t25*t5-2.26200
     =         76536473846d+4*t20*t21*t25*t5+1.2925758020842196d+4*t18*t
     >         20*t25*t5+1.2925758020842196d+4*t17*t21*t22*t5-2.95445897
     ?         6192502d+4*t17*t18*t22*t5+1.6618831741082823d+4*t15*t17*t
     @         22*t5+1.6618831741082823d+4*t14*t18*t19*t5-3.739237141743
     1         6354d+4*t14*t15*t19*t5+2.0773539676353528d+4*t12*t14*t19*
     2         t5+2.0773539676353528d+4*t11*t15*t16*t5-4.616342150300784
     3         d+4*t11*t12*t16*t5)*wght-9.305257363490997d-1*t41*t5*(1.0
     4         369691626335495d-6*gammabb**3*t68/rhob**10-3.894892828638
     5         6303d-4*t46*t47/rhob**7.333333333333333d+0+3.531915643569
     6         816d-2*gammabb*t3/rhob**4.666666666666667d+0)*wght-1.8610
     7         514726981994d+0*t44*t5*t50*wght-2.481401963597599d+0*t41*
     8         t45*t50*wght-2.481401963597599d+0*t4*t44*t45*wght-4.13566
     9         99393293317d-1*t4*t41*wght/t42+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490997d-1*t41*t5*(-3.88863
     1         43598758105d-7*t46*t68/rhob**9-9.632497209735863d-3*t3*t4
     2         9+1.29829760954621d-4*gammabb*t47*t48)*wght-9.30525736349
     3         0997d-1*t44*t5*t52*wght-1.2407009817987995d+0*t41*t45*t52
     4         *wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t41*t5*(1.4582378849534292d-7*gammabb*t68/rhob**8-3.
     2         2457440238655255d-5*t47*t51)*wght
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = -9.305257363490997d-1*t4*t5*(4.01086388
     1         76641245d+3*t42*t6*t66*t9+6.077066496460795d+1*t35*t39*t4
     2         2*t8-1.8231199489382385d+2*t33*t35*t42*t8-3.3423865730534
     3         37d+3*t11*t16*t42*t6+3.342386573053437d+3*t12*t13*t42*t43
     4         -6.077066496460795d+1*t36*t38*t42+1.8231199489382385d+2*t
     5         32*t36*t37*t42-3.646239897876477d+2*t30*t32*t37*t42+3.646
     6         239897876477d+2*t29*t33*t34*t42-6.077066496460795d+2*t27*
     7         t29*t34*t42+6.077066496460795d+2*t26*t30*t31*t42-9.115599
     8         744691192d+2*t24*t26*t31*t42+9.115599744691192d+2*t23*t27
     9         *t28*t42-1.276183964256767d+3*t21*t23*t28*t42+1.276183964
     :         256767d+3*t20*t24*t25*t42-1.7015786190090226d+3*t18*t20*t
     ;         25*t42+1.7015786190090226d+3*t17*t21*t22*t42-2.1877439387
     <         25886d+3*t15*t17*t22*t42+2.187743938725886d+3*t14*t18*t19
     =         *t42-2.734679923407357d+3*t12*t14*t19*t42+2.7346799234073
     >         57d+3*t11*t15*t16*t42)*wght-9.305257363490997d-1*t5*t50*t
     ?         69*wght-1.2407009817987995d+0*t4*t45*t69*wght+Mmat2(iq,D2
     @         _RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t5*t52*t65*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t4*t5*(5.28d+2*t6*t66*t9+t35*(1.6d+1*t36*t8+2.4d+1*t33*t
     2         8-4.0d+0*t63)+t32*(-6.0d+0*t62+3.6d+1*t33*t37+4.8d+1*t30*
     3         t37)+t29*(-8.0d+0*t61+6.4d+1*t30*t34+8.0d+1*t27*t34)+t26*
     4         (-1.0d+1*t60+1.0d+2*t27*t31+1.2d+2*t24*t31)+t11*(4.4d+2*t
     5         16*t6-2.0d+1*t55+4.0d+2*t12*t16)+t43*(4.84d+2*t13*t6-2.2d
     6         +1*t54)+t23*(-1.2d+1*t59+1.44d+2*t24*t28+1.68d+2*t21*t28)
     7         +t20*(-1.4d+1*t58+1.96d+2*t21*t25+2.24d+2*t18*t25)+t17*(-
     8         1.6d+1*t57+2.56d+2*t18*t22+2.88d+2*t15*t22)+t14*(-1.8d+1*
     9         t56+3.24d+2*t15*t19+3.6d+2*t12*t19)+t38*(8.0d+0*t39+8.0d+
     :         0*t36))*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_m06_2x functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_4} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_5} &=& {\it t_4}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_6} &=& 2.0\,\tau_\alpha+{\it t_4}\\\\ 
C>   {\it t_7} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{16}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_{19}} &=& {\it t_{18}}-2.0\,\tau_\beta\\\\ 
C>   {\it t_{20}} &=& 2.0\,\tau_\beta+{\it t_{18}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{22}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{23}} &=& {\it t_{22}}-2.0\,\tau_s\\\\ 
C>   {\it t_{24}} &=& 2.0\,\tau_s+{\it t_{22}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{17}}\,
C>    \sigma_{\beta\beta}}\over{0.00449276922095889\,{\it t_{17}}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left({{{\it t_{16}}
C>    \,{\it t_{19}}}\over{{\it t_{20}}}}+{{{\it t_{15}}\,{
C>    \it t_{19}}^2}\over{{\it t_{20}}^2}}+{{{\it t_{14}}\,{
C>    \it t_{19}}^3}\over{{\it t_{20}}^3}}+{{{\it t_{13}}\,{
C>    \it t_{19}}^4}\over{{\it t_{20}}^4}}+{{{\it t_{12}}\,{
C>    \it t_{19}}^5}\over{{\it t_{20}}^5}}+{{{\it t_{11}}\,{
C>    \it t_{19}}^6}\over{{\it t_{20}}^6}}+{{{\it t_{10}}\,{
C>    \it t_{19}}^7}\over{{\it t_{20}}^7}}+{{{\it t_9}\,{
C>    \it t_{19}}^8}\over{{\it t_{20}}^8}}+{{{\it t_8}\,{
C>    \it t_{19}}^9}\over{{\it t_{20}}^9}}+{{{\it t_7}\,{
C>    \it t_{19}}^{10}}\over{{\it t_{20}}^{10}}}+{{{\it t_3}\,{
C>    \it t_{19}}^{11}}\over{{\it t_{20}}^{11}}}+{\it t_2}\right)
C>    -0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.00449276922095889\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left({{{\it t_{16}}
C>    \,{\it t_5}}\over{{\it t_6}}}+{{{\it t_{15}}\,{
C>    \it t_5}^2}\over{{\it t_6}^2}}+{{{\it t_{14}}\,{
C>    \it t_5}^3}\over{{\it t_6}^3}}+{{{\it t_{13}}\,{
C>    \it t_5}^4}\over{{\it t_6}^4}}+{{{\it t_{12}}\,{
C>    \it t_5}^5}\over{{\it t_6}^5}}+{{{\it t_{11}}\,{
C>    \it t_5}^6}\over{{\it t_6}^6}}+{{{\it t_{10}}\,{
C>    \it t_5}^7}\over{{\it t_6}^7}}+{{{\it t_9}\,{\it t_5}^8}
C>    \over{{\it t_6}^8}}+{{{\it t_8}\,{\it t_5}^9}\over{{\it t_6}^9}}
C>    +{{{\it t_7}\,{\it t_5}^{10}}\over{{\it t_6}^{10}}}+{{{
C>    \it t_3}\,{\it t_5}^{11}}\over{{\it t_6}^{11}}}+{
C>    \it t_2}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{21}}\,\sigma_{ss}}
C>    \over{0.00449276922095889\,{\it t_{21}}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left({{{\it t_{16}}\,{\it t_{23}}}\over{{
C>    \it t_{24}}}}+{{{\it t_{15}}\,{\it t_{23}}^2}\over{{
C>    \it t_{24}}^2}}+{{{\it t_{14}}\,{\it t_{23}}^3}\over{{
C>    \it t_{24}}^3}}+{{{\it t_{13}}\,{\it t_{23}}^4}\over{{
C>    \it t_{24}}^4}}+{{{\it t_{12}}\,{\it t_{23}}^5}\over{{
C>    \it t_{24}}^5}}+{{{\it t_{11}}\,{\it t_{23}}^6}\over{{
C>    \it t_{24}}^6}}+{{{\it t_{10}}\,{\it t_{23}}^7}\over{{
C>    \it t_{24}}^7}}+{{{\it t_9}\,{\it t_{23}}^8}\over{{
C>    \it t_{24}}^8}}+{{{\it t_8}\,{\it t_{23}}^9}\over{{
C>    \it t_{24}}^9}}+{{{\it t_7}\,{\it t_{23}}^{10}}\over{{
C>    \it t_{24}}^{10}}}+{{{\it t_3}\,{\it t_{23}}^{11}}\over{{
C>    \it t_{24}}^{11}}}+{\it t_2}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 123, 194101 (2006)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2370993 "> 
C> 10.1063/1.2370993 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_m06_2x_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(11)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(10)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(9)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(8)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(7)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(6)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(5)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(4)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(3)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(2)
            t41 = 1/t11
            t42 = t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+t31*t32*
     1         t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t21+t16*t
     2         17*t18+t13*t14*t15+param(1)
            t43 = rhoa**6.666666666666666d-1
            t44 = 1/t11**12
            t45 = -1.5192666241151986d+1*t39*t40*t43*t9+3.03853324823039
     1         7d+1*t37*t39*t43*t9-1.6711932865267184d+2*t10*t43*t44*t6+
     2         1.6711932865267184d+2*t12*t14*t43*t6+1.5192666241151986d+
     3         1*t40*t41*t43-3.038533248230397d+1*t36*t37*t38*t43+4.5577
     4         998723455965d+1*t34*t36*t38*t43-4.5577998723455965d+1*t33
     5         *t34*t35*t43+6.077066496460795d+1*t31*t33*t35*t43-6.07706
     6         6496460795d+1*t30*t31*t32*t43+7.596333120575993d+1*t28*t3
     7         0*t32*t43-7.596333120575993d+1*t27*t28*t29*t43+9.11559974
     8         4691193d+1*t25*t27*t29*t43-9.115599744691193d+1*t24*t25*t
     9         26*t43+1.0634866368806391d+2*t22*t24*t26*t43-1.0634866368
     :         806391d+2*t21*t22*t23*t43+1.2154132992921589d+2*t19*t21*t
     ;         23*t43-1.2154132992921589d+2*t18*t19*t20*t43+1.3673399617
     <         036788d+2*t16*t18*t20*t43-1.3673399617036788d+2*t15*t16*t
     =         17*t43+1.5192666241151986d+2*t13*t15*t17*t43-1.5192666241
     >         151986d+2*t12*t13*t14*t43
            t46 = rhoa**3.333333333333333d-1
            t47 = gammaaa**2
            t48 = 1/t2**2
            t49 = 1/rhoa**6.333333333333333d+0
            t50 = 1/rhoa**3.6666666666666664d+0
            t51 = 4.327658698487368d-5*t47*t48*t49-9.632497209735863d-3*
     1         gammaaa*t3*t50
            t52 = 1/rhoa**5.333333333333333d+0
            t53 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t48*t52
            t54 = -2.2d+1*t10*t44*t6
            t55 = -2.2d+1*t14*t6-2.0d+1*t13*t14
            t56 = -1.8d+1*t16*t17-2.0d+1*t13*t17
            t57 = -1.6d+1*t19*t20-1.8d+1*t16*t20
            t58 = -1.4d+1*t22*t23-1.6d+1*t19*t23
            t59 = -1.2d+1*t25*t26-1.4d+1*t22*t26
            t60 = -1.0d+1*t28*t29-1.2d+1*t25*t29
            t61 = -8.0d+0*t31*t32-1.0d+1*t28*t32
            t62 = -6.0d+0*t34*t35-8.0d+0*t31*t35
            t63 = -4.0d+0*t37*t38-6.0d+0*t34*t38
            t64 = -2.0d+0*t40*t9-4.0d+0*t37*t9
            t65 = -2.0d+0*t40*t41
            t66 = t65+t39*t64+t36*t63+t33*t62+t30*t61+t27*t60+t24*t59+t2
     1         1*t58+t18*t57+t15*t56+t12*t55+t54
            t67 = 1/t11**13
            t68 = 1/t46
            t69 = -1.0128444160767991d+1*t39*t40*t68*t9+2.02568883215359
     1         82d+1*t37*t39*t68*t9+4.6163421503007845d+2*t36*t40*t5*t9-
     2         1.846536860120314d+3*t36*t37*t5*t9+1.3849026450902355d+3*
     3         t34*t36*t5*t9-1.1141288576844789d+2*t10*t44*t6*t68+1.1141
     4         288576844789d+2*t12*t14*t6*t68+1.0128444160767991d+1*t40*
     5         t41*t68-2.0256888321535982d+1*t36*t37*t38*t68+3.038533248
     6         230397d+1*t34*t36*t38*t68-3.038533248230397d+1*t33*t34*t3
     7         5*t68+4.0513776643071964d+1*t31*t33*t35*t68-4.05137766430
     8         71964d+1*t30*t31*t32*t68+5.064222080383994d+1*t28*t30*t32
     9         *t68-5.064222080383994d+1*t27*t28*t29*t68+6.0770664964607
     :         95d+1*t25*t27*t29*t68-6.077066496460795d+1*t24*t25*t26*t6
     ;         8+7.089910912537593d+1*t22*t24*t26*t68-7.089910912537593d
     <         +1*t21*t22*t23*t68+8.102755328614393d+1*t19*t21*t23*t68-8
     =         .102755328614393d+1*t18*t19*t20*t68+9.115599744691191d+1*
     >         t16*t18*t20*t68-9.115599744691191d+1*t15*t16*t17*t68+1.01
     ?         2844416076799d+2*t13*t15*t17*t68-1.012844416076799d+2*t12
     @         *t13*t14*t68+3.046785819198518d+4*t10*t5*t6*t67-5.5857740
     1         01863949d+4*t14*t44*t5*t6+2.5389881826654315d+4*t12*t17*t
     2         5*t6+2.5389881826654315d+4*t13*t14*t44*t5-4.6163421503007
     3         845d+2*t39*t40*t5+4.6163421503007845d+2*t37*t39*t5+1.3849
     4         026450902355d+3*t33*t37*t38*t5-4.1547079352707067d+3*t33*
     5         t34*t38*t5+2.769805290180471d+3*t31*t33*t38*t5+2.76980529
     6         0180471d+3*t30*t34*t35*t5-7.386147440481255d+3*t30*t31*t3
     7         5*t5+4.616342150300784d+3*t28*t30*t35*t5+4.61634215030078
     8         4d+3*t27*t31*t32*t5-1.154085537575196d+4*t27*t28*t32*t5+6
     9         .924513225451178d+3*t25*t27*t32*t5+6.924513225451178d+3*t
     :         24*t28*t29*t5-1.661883174108283d+4*t24*t25*t29*t5+9.69431
     ;         8515631648d+3*t22*t24*t29*t5+9.694318515631648d+3*t21*t25
     <         *t26*t5-2.2620076536473846d+4*t21*t22*t26*t5+1.2925758020
     =         842196d+4*t19*t21*t26*t5+1.2925758020842196d+4*t18*t22*t2
     >         3*t5-2.954458976192502d+4*t18*t19*t23*t5+1.66188317410828
     ?         23d+4*t16*t18*t23*t5+1.6618831741082823d+4*t15*t19*t20*t5
     @         -3.7392371417436354d+4*t15*t16*t20*t5+2.0773539676353528d
     1         +4*t13*t15*t20*t5+2.0773539676353528d+4*t12*t16*t17*t5-4.
     2         616342150300784d+4*t12*t13*t17*t5
            t70 = 1/t43
            t71 = gammaaa**3
            t72 = 1/t2**3
            t73 = 1/rhoa**10
            t74 = 1/rhoa**7.333333333333333d+0
            t75 = 1/rhoa**4.666666666666667d+0
            t76 = 3.531915643569816d-2*gammaaa*t3*t75-3.8948928286386303
     1         d-4*t47*t48*t74+1.0369691626335495d-6*t71*t72*t73
            t77 = 1/rhoa**9
            t78 = -3.8886343598758105d-7*t47*t72*t77-9.632497209735863d-
     1         3*t3*t50+1.29829760954621d-4*gammaaa*t48*t49
            t79 = 1/rhoa**8
            t80 = 1.4582378849534292d-7*gammaaa*t72*t79-3.24574402386552
     1         55d-5*t48*t52
            t81 = 6.077066496460795d+1*t36*t40*t43*t9-1.8231199489382385
     1         d+2*t34*t36*t43*t9+4.0108638876641245d+3*t10*t43*t6*t67-3
     2         .342386573053437d+3*t12*t17*t43*t6+3.342386573053437d+3*t
     3         13*t14*t43*t44-6.077066496460795d+1*t37*t39*t43+1.8231199
     4         489382385d+2*t33*t37*t38*t43-3.646239897876477d+2*t31*t33
     5         *t38*t43+3.646239897876477d+2*t30*t34*t35*t43-6.077066496
     6         460795d+2*t28*t30*t35*t43+6.077066496460795d+2*t27*t31*t3
     7         2*t43-9.115599744691192d+2*t25*t27*t32*t43+9.115599744691
     8         192d+2*t24*t28*t29*t43-1.276183964256767d+3*t22*t24*t29*t
     9         43+1.276183964256767d+3*t21*t25*t26*t43-1.701578619009022
     :         6d+3*t19*t21*t26*t43+1.7015786190090226d+3*t18*t22*t23*t4
     ;         3-2.187743938725886d+3*t16*t18*t23*t43+2.187743938725886d
     <         +3*t15*t19*t20*t43-2.734679923407357d+3*t13*t15*t20*t43+2
     =         .734679923407357d+3*t12*t16*t17*t43
            t82 = -2.0d+0*t39*t40*t9-4.0d+0*t37*t39*t9+t65-2.2d+1*t12*t1
     1         4*t6+t54-4.0d+0*t36*t37*t38-6.0d+0*t34*t36*t38-6.0d+0*t33
     2         *t34*t35-8.0d+0*t31*t33*t35-8.0d+0*t30*t31*t32-1.0d+1*t28
     3         *t30*t32-1.0d+1*t27*t28*t29-1.2d+1*t25*t27*t29-1.2d+1*t24
     4         *t25*t26-1.4d+1*t22*t24*t26-1.4d+1*t21*t22*t23-1.6d+1*t19
     5         *t21*t23-1.6d+1*t18*t19*t20-1.8d+1*t16*t18*t20-1.8d+1*t15
     6         *t16*t17-2.0d+1*t13*t15*t17-2.0d+1*t12*t13*t14
            t83 = 4.84d+2*t14*t6-2.2d+1*t55
            t84 = 4.0d+2*t13*t17
            t85 = 4.4d+2*t17*t6
            t86 = t85+t84-2.0d+1*t56
            t87 = 3.24d+2*t16*t20
            t88 = 3.6d+2*t13*t20
            t89 = t88+t87-1.8d+1*t57
            t90 = 2.56d+2*t19*t23
            t91 = 2.88d+2*t16*t23
            t92 = t91+t90-1.6d+1*t58
            t93 = 1.96d+2*t22*t26
            t94 = 2.24d+2*t19*t26
            t95 = t94+t93-1.4d+1*t59
            t96 = 1.44d+2*t25*t29
            t97 = 1.68d+2*t22*t29
            t98 = t97+t96-1.2d+1*t60
            t99 = 1.0d+2*t28*t32
            t100 = 1.2d+2*t25*t32
            t101 = t99-1.0d+1*t61+t100
            t102 = 6.4d+1*t31*t35
            t103 = 8.0d+1*t28*t35
            t104 = -8.0d+0*t62+t103+t102
            t105 = 3.6d+1*t34*t38
            t106 = 4.8d+1*t31*t38
            t107 = -6.0d+0*t63+t106+t105
            t108 = 1.6d+1*t37*t9
            t109 = 2.4d+1*t34*t9
            t110 = -4.0d+0*t64+t109+t108
            t111 = 8.0d+0*t37
            t112 = 8.0d+0*t40+t111
            t113 = rhoa**2
            t114 = 1/t11**14
            t115 = 1/t5
            t116 = 1/t2**4
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t4*t42*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t42*t5*t51-9.3052573
     1         63490997d-1*t4*t45*t5-1.2407009817987995d+0*t4*t42*t46)*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t42*t
     1         5*t53*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         66*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t42*t5*t76-4.135
     1         6699393293317d-1*t4*t42*t70-9.305257363490997d-1*t4*t5*t6
     2         9-1.8610514726981994d+0*t45*t5*t51-2.481401963597599d+0*t
     3         42*t46*t51-2.481401963597599d+0*t4*t45*t46)*wght+Amat2(iq
     4         ,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t42*t5*t78-9.30
     1         5257363490997d-1*t45*t5*t53-1.2407009817987995d+0*t42*t46
     2         *t53)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t42*t5*t80*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t5*t51*t82-1.240
     1         7009817987995d+0*t4*t46*t82-9.305257363490997d-1*t4*t5*t8
     2         1)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t53*t66*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(t24*t98+t21*t95+t18*t92+t15*t89+t12*t86+t44*t83+5
     2         .28d+2*t10*t6*t67+t112*t39+t110*t36+t107*t33+t104*t30+t10
     3         1*t27)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490997d-1*t4*t5*(9.2326
     1         84300601568d+2*t36*t40*t46*t9-3.693073720240627d+3*t36*t3
     2         7*t46*t9+2.769805290180471d+3*t34*t36*t46*t9+3.3761480535
     3         8933d+0*t115*t39*t40*t9-2.104036366334451d+4*t113*t33*t40
     4         *t9-6.75229610717866d+0*t115*t37*t39*t9+1.262421819800670
     5         4d+5*t113*t33*t37*t9-1.8936327297010064d+5*t113*t33*t34*t
     6         9+8.416145465337804d+4*t113*t31*t33*t9+6.093571638397036d
     7         +4*t10*t46*t6*t67+1.5275304019588112d+7*t113*t14*t6*t67-4
     8         .628880005935792d+6*t113*t13*t14*t67-1.1171548003727898d+
     9         5*t14*t44*t46*t6+5.077976365330863d+4*t12*t17*t46*t6-1.27
     :         29420016323426d+7*t113*t17*t44*t6+3.713762858948263d+1*t1
     ;         0*t115*t44*t6+3.4716600044518436d+6*t113*t12*t20*t6-3.713
     <         762858948263d+1*t115*t12*t14*t6-6.01754400771653d+6*t10*t
     =         113*t114*t6+5.077976365330862d+4*t13*t14*t44*t46-9.232684
     >         300601568d+2*t39*t40*t46+9.232684300601568d+2*t37*t39*t46
     ?         +2.769805290180471d+3*t33*t37*t38*t46-8.309415870541413d+
     @         3*t33*t34*t38*t46+5.539610580360942d+3*t31*t33*t38*t46+5.
     1         539610580360942d+3*t30*t34*t35*t46-1.4772294880962508d+4*
     2         t30*t31*t35*t46+9.232684300601568d+3*t28*t30*t35*t46+9.23
     3         2684300601568d+3*t27*t31*t32*t46-2.308171075150392d+4*t27
     4         *t28*t32*t46+1.3849026450902355d+4*t25*t27*t32*t46+1.3849
     5         026450902355d+4*t24*t28*t29*t46-3.323766348216566d+4*t24*
     6         t25*t29*t46+1.9388637031263292d+4*t22*t24*t29*t46+1.93886
     7         37031263292d+4*t21*t25*t26*t46-4.5240153072947686d+4*t21*
     8         t22*t26*t46+2.585151604168439d+4*t19*t21*t26*t46+2.585151
     9         604168439d+4*t18*t22*t23*t46-5.908917952385003d+4*t18*t19
     :         *t23*t46+3.3237663482165647d+4*t16*t18*t23*t46+3.32376634
     ;         82165647d+4*t15*t19*t20*t46-7.478474283487271d+4*t15*t16*
     <         t20*t46+4.154707935270705d+4*t13*t15*t20*t46+4.1547079352
     =         707056d+4*t12*t16*t17*t46-9.232684300601568d+4*t12*t13*t1
     >         7*t46-3.4716600044518436d+6*t113*t16*t17*t44+1.1572200014
     ?         839477d+7*t113*t13*t17*t44-3.37614805358933d+0*t115*t40*t
     @         41+2.104036366334451d+4*t113*t36*t40+6.75229610717866d+0*
     1         t115*t36*t37*t38-8.416145465337804d+4*t113*t30*t37*t38-1.
     2         0128444160767991d+1*t115*t34*t36*t38+3.787265459402013d+5
     3         *t113*t30*t34*t38-5.049687279202681d+5*t113*t30*t31*t38+2
     4         .1040363663344508d+5*t113*t28*t30*t38-4.208072732668902d+
     5         4*t113*t36*t37+2.104036366334451d+4*t113*t34*t36+1.012844
     6         4160767991d+1*t115*t33*t34*t35-2.104036366334451d+5*t113*
     7         t27*t34*t35-1.350459221435732d+1*t115*t31*t33*t35+8.41614
     8         5465337803d+5*t113*t27*t31*t35-1.0520181831672254d+6*t113
     9         *t27*t28*t35+4.208072732668902d+5*t113*t25*t27*t35+1.3504
     :         59221435732d+1*t115*t30*t31*t32-4.2080727326689016d+5*t11
     ;         3*t24*t31*t32-1.6880740267946648d+1*t115*t28*t30*t32+1.57
     <         8027274750838d+6*t113*t24*t28*t32-1.893632729701006d+6*t1
     =         13*t24*t25*t32+7.364127282170578d+5*t113*t22*t24*t32+1.68
     >         80740267946648d+1*t115*t27*t28*t29-7.364127282170578d+5*t
     ?         113*t21*t28*t29-2.0256888321535982d+1*t115*t25*t27*t29+2.
     @         6510858215814087d+6*t113*t21*t25*t29-3.0929334585116436d+
     1         6*t113*t21*t22*t29+1.1782603651472927d+6*t113*t19*t21*t29
     2         +2.0256888321535982d+1*t115*t24*t25*t26-1.178260365147292
     3         7d+6*t113*t18*t25*t26-2.3633036375125308d+1*t115*t22*t24*
     4         t26+4.123911278015524d+6*t113*t18*t22*t26-4.7130414605891
     5         71d+6*t113*t18*t19*t26+1.767390547720939d+6*t113*t16*t18*
     6         t26+2.3633036375125308d+1*t115*t21*t22*t23-1.767390547720
     7         9387d+6*t113*t15*t22*t23-2.700918442871464d+1*t115*t19*t2
     8         1*t23+6.059624735043219d+6*t113*t15*t19*t23-6.81707782692
     9         3622d+6*t113*t15*t16*t23+2.524843639601341d+6*t113*t13*t1
     :         5*t23+2.700918442871464d+1*t115*t18*t19*t20-2.52484363960
     ;         1341d+6*t113*t12*t19*t20-3.0385332482303967d+1*t115*t16*t
     <         18*t20+8.521347283654526d+6*t113*t12*t16*t20-9.4681636485
     =         05028d+6*t113*t12*t13*t20+3.0385332482303967d+1*t115*t15*
     >         t16*t17-3.3761480535893296d+1*t115*t13*t15*t17+3.37614805
     ?         35893296d+1*t115*t12*t13*t14)-2.791577209047299d+0*t45*t5
     @         *t76-3.7221029453963983d+0*t42*t46*t76-9.305257363490997d
     1         -1*t42*t5*(-1.9702414090037443d-5*t71*t72/rhoa**11+3.2794
     2         03591520427d-3*t47*t48/rhoa**8.333333333333334d+0-1.64822
     3         7300332581d-1*gammaaa*t3/rhoa**5.666666666666667d+0+3.727
     4         0905095708196d-8*gammaaa**4*t116/rhoa**1.3666666666666666
     5         d+1)+(-1.2407009817987995d+0*t42*t51-1.2407009817987995d+
     6         0*t4*t45)*t70+2.757113292886221d-1*t4*t42/t7-2.7915772090
     7         47299d+0*t5*t51*t69-3.7221029453963983d+0*t4*t46*t69-7.44
     8         4205890792797d+0*t45*t46*t51)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-1.8610514726981994d+0*t45*t5*t78-
     1         2.481401963597599d+0*t42*t46*t78-9.305257363490997d-1*t42
     2         *t5*(3.531915643569816d-2*t3*t75-9.376593846722628d-4*gam
     3         maaa*t48*t74+6.610678411788878d-6*t47*t72*t73-1.397658941
     4         0890572d-8*t116*t71/rhoa**1.2666666666666665d+1)-4.135669
     5         9393293317d-1*t42*t53*t70-9.305257363490997d-1*t5*t53*t69
     6         -2.481401963597599d+0*t45*t46*t53)*wght+Cmat3(iq,D3_RA_RA
     7         _GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t45*t5*t80-
     1         1.2407009817987995d+0*t42*t46*t80-9.305257363490997d-1*t4
     2         2*t5*(-1.944317179937905d-6*gammaaa*t72*t77+1.73106347939
     3         4947d-4*t48*t49+5.241221029083964d-9*t116*t47/rhoa**1.166
     4         6666666666665d+1))*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t42*t5*(4.374713654860287d-7*t72*t79-1.96545
     2         7885906487d-9*gammaaa*t116/rhoa**1.0666666666666666d+1)*w
     3         ght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Mmat3(iq,D3_RA_RA_TA) = (-9.305257363490997d-1*t4*t5*(4.0513
     1         776643071964d+1*t36*t40*t68*t9-1.2154132992921589d+2*t34*
     2         t36*t68*t9-2.769805290180471d+3*t33*t40*t5*t9+5.539610580
     3         360942d+3*t33*t37*t5*t9+8.309415870541416d+3*t33*t34*t5*t
     4         9-1.1079221160721883d+4*t31*t33*t5*t9+2.6739092584427493d
     5         +3*t10*t6*t67*t68-2.228257715368958d+3*t12*t17*t6*t68+2.2
     6         282577153689576d+3*t13*t14*t44*t68-4.0513776643071964d+1*
     7         t37*t39*t68+1.2154132992921589d+2*t33*t37*t38*t68-2.43082
     8         65985843178d+2*t31*t33*t38*t68+2.4308265985843178d+2*t30*
     9         t34*t35*t68-4.051377664307196d+2*t28*t30*t35*t68+4.051377
     :         6643071964d+2*t27*t31*t32*t68-6.077066496460795d+2*t25*t2
     ;         7*t32*t68+6.077066496460793d+2*t24*t28*t29*t68-8.50789309
     <         5045111d+2*t22*t24*t29*t68+8.507893095045113d+2*t21*t25*t
     =         26*t68-1.1343857460060149d+3*t19*t21*t26*t68+1.1343857460
     >         060149d+3*t18*t22*t23*t68-1.4584959591505905d+3*t16*t18*t
     ?         23*t68+1.4584959591505908d+3*t15*t19*t20*t68-1.8231199489
     @         382383d+3*t13*t15*t20*t68+1.8231199489382383d+3*t12*t16*t
     1         17*t68+6.702928802236738d+5*t14*t5*t6*t67-6.0935716383970
     2         35d+5*t13*t14*t5*t67+5.585774001863949d+5*t17*t44*t5*t6-4
     3         .5701787287977763d+5*t12*t20*t5*t6-7.921643129916147d+5*t
     4         10*t114*t5*t6-4.5701787287977763d+5*t16*t17*t44*t5+5.0779
     5         76365330863d+5*t13*t17*t44*t5+9.232684300601569d+2*t36*t4
     6         0*t5-1.1079221160721883d+4*t30*t37*t38*t5+1.6618831741082
     7         832d+4*t30*t34*t38*t5+2.2158442321443766d+4*t30*t31*t38*t
     8         5-2.7698052901804704d+4*t28*t30*t38*t5+1.846536860120314d
     9         +3*t36*t37*t5-2.769805290180471d+3*t34*t36*t5-2.769805290
     :         180471d+4*t27*t34*t35*t5+3.693073720240629d+4*t27*t31*t35
     ;         *t5+4.616342150300784d+4*t27*t28*t35*t5-5.539610580360942
     <         d+4*t25*t27*t35*t5-5.539610580360941d+4*t24*t31*t32*t5+6.
     =         924513225451174d+4*t24*t28*t32*t5+8.309415870541416d+4*t2
     >         4*t25*t32*t5-9.694318515631648d+4*t22*t24*t32*t5-9.694318
     ?         515631648d+4*t21*t28*t29*t5+1.1633182218757979d+5*t21*t25
     @         *t29*t5+1.3572045921884307d+5*t21*t22*t29*t5-1.5510909625
     1         010635d+5*t19*t21*t29*t5-1.5510909625010638d+5*t18*t25*t2
     2         6*t5+1.809606122917908d+5*t18*t22*t26*t5+2.06812128333475
     3         15d+5*t18*t19*t26*t5-2.3266364437515952d+5*t16*t18*t26*t5
     4         -2.3266364437515955d+5*t15*t22*t23*t5+2.659013078573252d+
     5         5*t15*t19*t23*t5+2.9913897133949086d+5*t15*t16*t23*t5-3.3
     6         237663482165647d+5*t13*t15*t23*t5-3.323766348216565d+5*t1
     7         2*t19*t20*t5+3.7392371417436354d+5*t12*t16*t20*t5+4.15470
     8         79352707056d+5*t12*t13*t20*t5)-9.305257363490997d-1*t5*t7
     9         6*t82-4.1356699393293317d-1*t4*t70*t82-2.481401963597599d
     :         +0*t46*t51*t82-1.8610514726981994d+0*t5*t51*t81-2.4814019
     ;         63597599d+0*t4*t46*t81)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_GAA_TA) = (-9.305257363490997d-1*t5*t78*t82-1
     1         .2407009817987995d+0*t46*t53*t82-9.305257363490997d-1*t5*
     2         t53*t81)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)-9.30525736
     1         3490997d-1*t5*t66*t80*wght
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)-9.305257363490
     1         997d-1*t4*t5*(t24*(-1.2d+1*(t99+t100)-1.44d+3*t25*t32-1.6
     2         8d+3*t22*t32-1.2d+1*t101)+t21*(-1.4d+1*t98-1.4d+1*(t97+t9
     3         6)-2.352d+3*t22*t29-2.688d+3*t19*t29)+t18*(-1.6d+1*t95-1.
     4         6d+1*(t94+t93)-3.584d+3*t19*t26-4.032d+3*t16*t26)+t15*(-1
     5         .8d+1*t92-1.8d+1*(t91+t90)-5.184d+3*t16*t23-5.76d+3*t13*t
     6         23)+t33*(-1.44d+2*t34*t9-1.92d+2*t31*t9-6.0d+0*t110-6.0d+
     7         0*(t109+t108))+t12*(-2.0d+1*t89-2.0d+1*(t88+t87)-7.92d+3*
     8         t20*t6-7.2d+3*t13*t20)+t44*(-2.2d+1*t86-2.2d+1*(t85+t84)-
     9         9.68d+3*t17*t6)+t67*(-2.4d+1*t83-1.1616d+4*t14*t6)-1.3728
     :         d+4*t10*t114*t6+t36*(-4.0d+0*(4.0d+0*t40+t111)-3.2d+1*t37
     ;         -4.8d+1*t34-4.0d+0*t112)+t30*(-3.84d+2*t31*t38-4.8d+2*t28
     <         *t38-8.0d+0*t107-8.0d+0*(t106+t105))+t27*(-8.0d+2*t28*t35
     =         -9.6d+2*t25*t35-1.0d+1*t104-1.0d+1*(t103+t102)))*wght
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
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
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(1)
            t7 = param(12)
            t8 = rhoa**1.6666666666666669d+0
            t9 = 9.115599744691192d+0*t8
            t10 = t9-2.0d+0*taua
            t11 = t10**11
            t12 = 2.0d+0*taua+t9
            t13 = 1/t12**11
            t14 = param(11)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(10)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(9)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(8)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(7)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(6)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(5)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(4)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(3)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(2)
            t42 = 1/t12
            t43 = t11*t13*t7+t6+t10*t41*t42+t38*t39*t40+t35*t36*t37+t32*
     1         t33*t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*t22+t
     2         17*t18*t19+t14*t15*t16
            t44 = 1/rhob**2.6666666666666666d+0
            t45 = 4.49276922095889d-3*gammabb*t44+1.0d+0
            t46 = 1/t45
            t47 = 3.6121864536509485d-3*gammabb*t44*t46+1.0d+0
            t48 = rhob**1.3333333333333333d+0
            t49 = rhob**1.6666666666666669d+0
            t50 = 9.115599744691192d+0*t49
            t51 = t50-2.0d+0*taub
            t52 = t51**11
            t53 = 2.0d+0*taub+t50
            t54 = 1/t53**11
            t55 = t51**10
            t56 = 1/t53**10
            t57 = t51**9
            t58 = 1/t53**9
            t59 = t51**8
            t60 = 1/t53**8
            t61 = t51**7
            t62 = 1/t53**7
            t63 = t51**6
            t64 = 1/t53**6
            t65 = t51**5
            t66 = 1/t53**5
            t67 = t51**4
            t68 = 1/t53**4
            t69 = t51**3
            t70 = 1/t53**3
            t71 = t51**2
            t72 = 1/t53**2
            t73 = 1/t53
            t74 = t41*t51*t73+t38*t71*t72+t35*t69*t70+t52*t54*t7+t32*t67
     1         *t68+t29*t65*t66+t26*t63*t64+t23*t61*t62+t20*t59*t60+t6+t
     2         17*t57*t58+t14*t55*t56
            t75 = rhoa**6.666666666666666d-1
            t76 = 1/t12**12
            t77 = -1.6711932865267184d+2*t11*t7*t75*t76+1.67119328652671
     1         84d+2*t13*t15*t7*t75+1.5192666241151986d+1*t41*t42*t75-1.
     2         5192666241151986d+1*t10*t40*t41*t75+3.038533248230397d+1*
     3         t10*t38*t40*t75-3.038533248230397d+1*t37*t38*t39*t75+4.55
     4         77998723455965d+1*t35*t37*t39*t75-4.5577998723455965d+1*t
     5         34*t35*t36*t75+6.077066496460795d+1*t32*t34*t36*t75-6.077
     6         066496460795d+1*t31*t32*t33*t75+7.596333120575993d+1*t29*
     7         t31*t33*t75-7.596333120575993d+1*t28*t29*t30*t75+9.115599
     8         744691193d+1*t26*t28*t30*t75-9.115599744691193d+1*t25*t26
     9         *t27*t75+1.0634866368806391d+2*t23*t25*t27*t75-1.06348663
     :         68806391d+2*t22*t23*t24*t75+1.2154132992921589d+2*t20*t22
     ;         *t24*t75-1.2154132992921589d+2*t19*t20*t21*t75+1.36733996
     <         17036788d+2*t17*t19*t21*t75-1.3673399617036788d+2*t16*t17
     =         *t18*t75+1.5192666241151986d+2*t14*t16*t18*t75-1.51926662
     >         41151986d+2*t13*t14*t15*t75
            t78 = rhoa**3.333333333333333d-1
            t79 = gammaaa**2
            t80 = 1/t2**2
            t81 = 1/rhoa**6.333333333333333d+0
            t82 = 1/rhoa**3.6666666666666664d+0
            t83 = 4.327658698487368d-5*t79*t80*t81-9.632497209735863d-3*
     1         gammaaa*t3*t82
            t84 = rhob**6.666666666666666d-1
            t85 = 1/t53**12
            t86 = -1.6711932865267184d+2*t52*t7*t84*t85+1.51926662411519
     1         86d+1*t41*t73*t84-1.5192666241151986d+1*t41*t51*t72*t84+3
     2         .038533248230397d+1*t38*t51*t72*t84-3.038533248230397d+1*
     3         t38*t70*t71*t84+4.5577998723455965d+1*t35*t70*t71*t84+1.6
     4         711932865267184d+2*t54*t55*t7*t84-4.5577998723455965d+1*t
     5         35*t68*t69*t84+6.077066496460795d+1*t32*t68*t69*t84-6.077
     6         066496460795d+1*t32*t66*t67*t84+7.596333120575993d+1*t29*
     7         t66*t67*t84-7.596333120575993d+1*t29*t64*t65*t84+9.115599
     8         744691193d+1*t26*t64*t65*t84-9.115599744691193d+1*t26*t62
     9         *t63*t84+1.0634866368806391d+2*t23*t62*t63*t84-1.06348663
     :         68806391d+2*t23*t60*t61*t84+1.2154132992921589d+2*t20*t60
     ;         *t61*t84-1.2154132992921589d+2*t20*t58*t59*t84+1.36733996
     <         17036788d+2*t17*t58*t59*t84-1.3673399617036788d+2*t17*t56
     =         *t57*t84+1.5192666241151986d+2*t14*t56*t57*t84-1.51926662
     >         41151986d+2*t14*t54*t55*t84
            t87 = rhob**3.333333333333333d-1
            t88 = gammabb**2
            t89 = 1/t45**2
            t90 = 1/rhob**6.333333333333333d+0
            t91 = 1/rhob**3.6666666666666664d+0
            t92 = 4.327658698487368d-5*t88*t89*t90-9.632497209735863d-3*
     1         gammabb*t46*t91
            t93 = 1/rhoa**5.333333333333333d+0
            t94 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t80*t93
            t95 = 1/rhob**5.333333333333333d+0
            t96 = 3.6121864536509485d-3*t44*t46-1.6228720119327628d-5*ga
     1         mmabb*t89*t95
            t97 = -2.2d+1*t11*t7*t76
            t98 = -2.2d+1*t15*t7-2.0d+1*t14*t15
            t99 = -1.8d+1*t17*t18-2.0d+1*t14*t18
            t100 = -1.6d+1*t20*t21-1.8d+1*t17*t21
            t101 = -1.4d+1*t23*t24-1.6d+1*t20*t24
            t102 = -1.2d+1*t26*t27-1.4d+1*t23*t27
            t103 = -1.0d+1*t29*t30-1.2d+1*t26*t30
            t104 = -8.0d+0*t32*t33-1.0d+1*t29*t33
            t105 = -6.0d+0*t35*t36-8.0d+0*t32*t36
            t106 = -4.0d+0*t38*t39-6.0d+0*t35*t39
            t107 = -2.0d+0*t10*t41-4.0d+0*t10*t38
            t108 = -2.0d+0*t41*t42
            t109 = t16*t99+t13*t98+t97+t107*t40+t106*t37+t105*t34+t104*t
     1         31+t103*t28+t102*t25+t101*t22+t100*t19+t108
            t110 = -2.2d+1*t52*t7*t85
            t111 = -2.2d+1*t55*t7-2.0d+1*t14*t55
            t112 = -1.8d+1*t17*t57-2.0d+1*t14*t57
            t113 = -1.6d+1*t20*t59-1.8d+1*t17*t59
            t114 = -1.4d+1*t23*t61-1.6d+1*t20*t61
            t115 = -1.2d+1*t26*t63-1.4d+1*t23*t63
            t116 = -1.0d+1*t29*t65-1.2d+1*t26*t65
            t117 = -8.0d+0*t32*t67-1.0d+1*t29*t67
            t118 = -6.0d+0*t35*t69-8.0d+0*t32*t69
            t119 = -4.0d+0*t38*t71-6.0d+0*t35*t71
            t120 = -2.0d+0*t41*t51-4.0d+0*t38*t51
            t121 = -2.0d+0*t41*t73
            t122 = t120*t72+t119*t70+t118*t68+t117*t66+t116*t64+t115*t62
     1         +t114*t60+t113*t58+t112*t56+t111*t54+t121+t110
            t123 = 1/t12**13
            t124 = 1/t78
            t125 = -5.585774001863949d+4*t15*t5*t7*t76-1.114128857684478
     1         9d+2*t11*t124*t7*t76+2.5389881826654315d+4*t14*t15*t5*t76
     2         +2.5389881826654315d+4*t13*t18*t5*t7+3.046785819198518d+4
     3         *t11*t123*t5*t7+1.1141288576844789d+2*t124*t13*t15*t7-4.6
     4         163421503007845d+2*t40*t41*t5+4.6163421503007845d+2*t10*t
     5         37*t41*t5+4.6163421503007845d+2*t38*t40*t5+1.384902645090
     6         2355d+3*t34*t38*t39*t5-4.1547079352707067d+3*t34*t35*t39*
     7         t5+2.769805290180471d+3*t32*t34*t39*t5-1.846536860120314d
     8         +3*t10*t37*t38*t5+1.3849026450902355d+3*t10*t35*t37*t5+2.
     9         769805290180471d+3*t31*t35*t36*t5-7.386147440481255d+3*t3
     :         1*t32*t36*t5+4.616342150300784d+3*t29*t31*t36*t5+4.616342
     ;         150300784d+3*t28*t32*t33*t5-1.154085537575196d+4*t28*t29*
     <         t33*t5+6.924513225451178d+3*t26*t28*t33*t5+6.924513225451
     =         178d+3*t25*t29*t30*t5-1.661883174108283d+4*t25*t26*t30*t5
     >         +9.694318515631648d+3*t23*t25*t30*t5+9.694318515631648d+3
     ?         *t22*t26*t27*t5-2.2620076536473846d+4*t22*t23*t27*t5+1.29
     @         25758020842196d+4*t20*t22*t27*t5+1.2925758020842196d+4*t1
     1         9*t23*t24*t5-2.954458976192502d+4*t19*t20*t24*t5+1.661883
     2         1741082823d+4*t17*t19*t24*t5+1.6618831741082823d+4*t16*t2
     3         0*t21*t5-3.7392371417436354d+4*t16*t17*t21*t5+2.077353967
     4         6353528d+4*t14*t16*t21*t5+2.0773539676353528d+4*t13*t17*t
     5         18*t5-4.616342150300784d+4*t13*t14*t18*t5+1.0128444160767
     6         991d+1*t124*t41*t42-1.0128444160767991d+1*t10*t124*t40*t4
     7         1+2.0256888321535982d+1*t10*t124*t38*t40-2.02568883215359
     8         82d+1*t124*t37*t38*t39+3.038533248230397d+1*t124*t35*t37*
     9         t39-3.038533248230397d+1*t124*t34*t35*t36+4.0513776643071
     :         964d+1*t124*t32*t34*t36-4.0513776643071964d+1*t124*t31*t3
     ;         2*t33+5.064222080383994d+1*t124*t29*t31*t33-5.06422208038
     <         3994d+1*t124*t28*t29*t30+6.077066496460795d+1*t124*t26*t2
     =         8*t30-6.077066496460795d+1*t124*t25*t26*t27+7.08991091253
     >         7593d+1*t124*t23*t25*t27-7.089910912537593d+1*t124*t22*t2
     ?         3*t24+8.102755328614393d+1*t124*t20*t22*t24-8.10275532861
     @         4393d+1*t124*t19*t20*t21+9.115599744691191d+1*t124*t17*t1
     1         9*t21-9.115599744691191d+1*t124*t16*t17*t18+1.01284441607
     2         6799d+2*t124*t14*t16*t18-1.012844416076799d+2*t124*t13*t1
     3         4*t15
            t126 = 1/t75
            t127 = gammaaa**3
            t128 = 1/t2**3
            t129 = 1/rhoa**10
            t130 = 1/rhoa**7.333333333333333d+0
            t131 = 1/rhoa**4.666666666666667d+0
            t132 = -3.8948928286386303d-4*t130*t79*t80+3.531915643569816
     1         d-2*gammaaa*t131*t3+1.0369691626335495d-6*t127*t128*t129
            t133 = 1/t53**13
            t134 = 1/t87
            t135 = -5.585774001863949d+4*t48*t55*t7*t85-1.11412885768447
     1         89d+2*t134*t52*t7*t85+2.5389881826654315d+4*t14*t48*t55*t
     2         85+1.0128444160767991d+1*t134*t41*t73-1.0128444160767991d
     3         +1*t134*t41*t51*t72+2.0256888321535982d+1*t134*t38*t51*t7
     4         2-4.6163421503007845d+2*t41*t48*t72+4.6163421503007845d+2
     5         *t38*t48*t72-2.0256888321535982d+1*t134*t38*t70*t71+3.038
     6         533248230397d+1*t134*t35*t70*t71+1.3849026450902355d+3*t3
     7         8*t48*t68*t71-4.1547079352707067d+3*t35*t48*t68*t71+2.769
     8         805290180471d+3*t32*t48*t68*t71+4.6163421503007845d+2*t41
     9         *t48*t51*t70-1.846536860120314d+3*t38*t48*t51*t70+1.38490
     :         26450902355d+3*t35*t48*t51*t70+2.5389881826654315d+4*t48*
     ;         t54*t57*t7+1.1141288576844789d+2*t134*t54*t55*t7+3.046785
     <         819198518d+4*t133*t48*t52*t7-3.038533248230397d+1*t134*t3
     =         5*t68*t69+4.0513776643071964d+1*t134*t32*t68*t69+2.769805
     >         290180471d+3*t35*t48*t66*t69-7.386147440481255d+3*t32*t48
     ?         *t66*t69+4.616342150300784d+3*t29*t48*t66*t69-4.051377664
     @         3071964d+1*t134*t32*t66*t67+5.064222080383994d+1*t134*t29
     1         *t66*t67+4.616342150300784d+3*t32*t48*t64*t67-1.154085537
     2         575196d+4*t29*t48*t64*t67+6.924513225451178d+3*t26*t48*t6
     3         4*t67-5.064222080383994d+1*t134*t29*t64*t65+6.07706649646
     4         0795d+1*t134*t26*t64*t65+6.924513225451178d+3*t29*t48*t62
     5         *t65-1.661883174108283d+4*t26*t48*t62*t65+9.6943185156316
     6         48d+3*t23*t48*t62*t65-6.077066496460795d+1*t134*t26*t62*t
     7         63+7.089910912537593d+1*t134*t23*t62*t63+9.69431851563164
     8         8d+3*t26*t48*t60*t63-2.2620076536473846d+4*t23*t48*t60*t6
     9         3+1.2925758020842196d+4*t20*t48*t60*t63-7.089910912537593
     :         d+1*t134*t23*t60*t61+8.102755328614393d+1*t134*t20*t60*t6
     ;         1+1.2925758020842196d+4*t23*t48*t58*t61-2.954458976192502
     <         d+4*t20*t48*t58*t61+1.6618831741082823d+4*t17*t48*t58*t61
     =         -8.102755328614393d+1*t134*t20*t58*t59+9.115599744691191d
     >         +1*t134*t17*t58*t59+1.6618831741082823d+4*t20*t48*t56*t59
     ?         -3.7392371417436354d+4*t17*t48*t56*t59+2.0773539676353528
     @         d+4*t14*t48*t56*t59-9.115599744691191d+1*t134*t17*t56*t57
     1         +1.012844416076799d+2*t134*t14*t56*t57+2.0773539676353528
     2         d+4*t17*t48*t54*t57-4.616342150300784d+4*t14*t48*t54*t57-
     3         1.012844416076799d+2*t134*t14*t54*t55
            t136 = 1/t84
            t137 = gammabb**3
            t138 = 1/t45**3
            t139 = 1/rhob**10
            t140 = 1/rhob**7.333333333333333d+0
            t141 = 1/rhob**4.666666666666667d+0
            t142 = -3.8948928286386303d-4*t140*t88*t89+3.531915643569816
     1         d-2*gammabb*t141*t46+1.0369691626335495d-6*t137*t138*t139
            t143 = 1/rhoa**9
            t144 = -9.632497209735863d-3*t3*t82+1.29829760954621d-4*gamm
     1         aaa*t80*t81-3.8886343598758105d-7*t128*t143*t79
            t145 = 1/rhob**9
            t146 = -9.632497209735863d-3*t46*t91+1.29829760954621d-4*gam
     1         mabb*t89*t90-3.8886343598758105d-7*t138*t145*t88
            t147 = 1/rhoa**8
            t148 = 1.4582378849534292d-7*gammaaa*t128*t147-3.24574402386
     1         55255d-5*t80*t93
            t149 = 1/rhob**8
            t150 = 1.4582378849534292d-7*gammabb*t138*t149-3.24574402386
     1         55255d-5*t89*t95
            t151 = 3.342386573053437d+3*t14*t15*t75*t76-3.34238657305343
     1         7d+3*t13*t18*t7*t75+4.0108638876641245d+3*t11*t123*t7*t75
     2         +6.077066496460795d+1*t10*t37*t41*t75-6.077066496460795d+
     3         1*t38*t40*t75+1.8231199489382385d+2*t34*t38*t39*t75-3.646
     4         239897876477d+2*t32*t34*t39*t75-1.8231199489382385d+2*t10
     5         *t35*t37*t75+3.646239897876477d+2*t31*t35*t36*t75-6.07706
     6         6496460795d+2*t29*t31*t36*t75+6.077066496460795d+2*t28*t3
     7         2*t33*t75-9.115599744691192d+2*t26*t28*t33*t75+9.11559974
     8         4691192d+2*t25*t29*t30*t75-1.276183964256767d+3*t23*t25*t
     9         30*t75+1.276183964256767d+3*t22*t26*t27*t75-1.70157861900
     :         90226d+3*t20*t22*t27*t75+1.7015786190090226d+3*t19*t23*t2
     ;         4*t75-2.187743938725886d+3*t17*t19*t24*t75+2.187743938725
     <         886d+3*t16*t20*t21*t75-2.734679923407357d+3*t14*t16*t21*t
     =         75+2.734679923407357d+3*t13*t17*t18*t75
            t152 = t97-2.2d+1*t13*t15*t7-2.0d+0*t10*t40*t41-4.0d+0*t10*t
     1         38*t40-4.0d+0*t37*t38*t39-6.0d+0*t35*t37*t39-6.0d+0*t34*t
     2         35*t36-8.0d+0*t32*t34*t36-8.0d+0*t31*t32*t33-1.0d+1*t29*t
     3         31*t33-1.0d+1*t28*t29*t30-1.2d+1*t26*t28*t30-1.2d+1*t25*t
     4         26*t27-1.4d+1*t23*t25*t27-1.4d+1*t22*t23*t24-1.6d+1*t20*t
     5         22*t24-1.6d+1*t19*t20*t21-1.8d+1*t17*t19*t21-1.8d+1*t16*t
     6         17*t18-2.0d+1*t14*t16*t18-2.0d+1*t13*t14*t15+t108
            t153 = 3.342386573053437d+3*t14*t55*t84*t85-6.07706649646079
     1         5d+1*t38*t72*t84+1.8231199489382385d+2*t38*t68*t71*t84-3.
     2         646239897876477d+2*t32*t68*t71*t84+6.077066496460795d+1*t
     3         41*t51*t70*t84-1.8231199489382385d+2*t35*t51*t70*t84-3.34
     4         2386573053437d+3*t54*t57*t7*t84+4.0108638876641245d+3*t13
     5         3*t52*t7*t84+3.646239897876477d+2*t35*t66*t69*t84-6.07706
     6         6496460795d+2*t29*t66*t69*t84+6.077066496460795d+2*t32*t6
     7         4*t67*t84-9.115599744691192d+2*t26*t64*t67*t84+9.11559974
     8         4691192d+2*t29*t62*t65*t84-1.276183964256767d+3*t23*t62*t
     9         65*t84+1.276183964256767d+3*t26*t60*t63*t84-1.70157861900
     :         90226d+3*t20*t60*t63*t84+1.7015786190090226d+3*t23*t58*t6
     ;         1*t84-2.187743938725886d+3*t17*t58*t61*t84+2.187743938725
     <         886d+3*t20*t56*t59*t84-2.734679923407357d+3*t14*t56*t59*t
     =         84+2.734679923407357d+3*t17*t54*t57*t84
            t154 = -2.0d+0*t41*t51*t72-4.0d+0*t38*t51*t72-4.0d+0*t38*t70
     1         *t71-6.0d+0*t35*t70*t71-2.2d+1*t54*t55*t7-6.0d+0*t35*t68*
     2         t69-8.0d+0*t32*t68*t69-8.0d+0*t32*t66*t67-1.0d+1*t29*t66*
     3         t67-1.0d+1*t29*t64*t65-1.2d+1*t26*t64*t65-1.2d+1*t26*t62*
     4         t63-1.4d+1*t23*t62*t63-1.4d+1*t23*t60*t61-1.6d+1*t20*t60*
     5         t61-1.6d+1*t20*t58*t59-1.8d+1*t17*t58*t59-1.8d+1*t17*t56*
     6         t57-2.0d+1*t14*t56*t57-2.0d+1*t14*t54*t55+t121+t110
            t155 = 4.84d+2*t15*t7-2.2d+1*t98
            t156 = 4.0d+2*t14*t18
            t157 = 4.4d+2*t18*t7
            t158 = -2.0d+1*t99+t157+t156
            t159 = 3.24d+2*t17*t21
            t160 = 3.6d+2*t14*t21
            t161 = t160+t159-1.8d+1*t100
            t162 = 2.56d+2*t20*t24
            t163 = 2.88d+2*t17*t24
            t164 = t163+t162-1.6d+1*t101
            t165 = 1.96d+2*t23*t27
            t166 = 2.24d+2*t20*t27
            t167 = t166+t165-1.4d+1*t102
            t168 = 1.44d+2*t26*t30
            t169 = 1.68d+2*t23*t30
            t170 = t169+t168-1.2d+1*t103
            t171 = 1.0d+2*t29*t33
            t172 = 1.2d+2*t26*t33
            t173 = t172+t171-1.0d+1*t104
            t174 = 6.4d+1*t32*t36
            t175 = 8.0d+1*t29*t36
            t176 = t175+t174-8.0d+0*t105
            t177 = 3.6d+1*t35*t39
            t178 = 4.8d+1*t32*t39
            t179 = t178+t177-6.0d+0*t106
            t180 = 1.6d+1*t10*t38
            t181 = 2.4d+1*t10*t35
            t182 = t181+t180-4.0d+0*t107
            t183 = 8.0d+0*t38
            t184 = 8.0d+0*t41+t183
            t185 = t155*t76+5.28d+2*t11*t123*t7+t184*t40+t182*t37+t179*t
     1         34+t176*t31+t173*t28+t170*t25+t167*t22+t164*t19+t16*t161+
     2         t13*t158
            t186 = 4.84d+2*t55*t7-2.2d+1*t111
            t187 = 4.0d+2*t14*t57
            t188 = 4.4d+2*t57*t7
            t189 = t188+t187-2.0d+1*t112
            t190 = 3.24d+2*t17*t59
            t191 = 3.6d+2*t14*t59
            t192 = t191+t190-1.8d+1*t113
            t193 = 2.56d+2*t20*t61
            t194 = 2.88d+2*t17*t61
            t195 = t194+t193-1.6d+1*t114
            t196 = 1.96d+2*t23*t63
            t197 = 2.24d+2*t20*t63
            t198 = t197+t196-1.4d+1*t115
            t199 = 1.44d+2*t26*t65
            t200 = 1.68d+2*t23*t65
            t201 = t200+t199-1.2d+1*t116
            t202 = 1.0d+2*t29*t67
            t203 = 1.2d+2*t26*t67
            t204 = t203+t202-1.0d+1*t117
            t205 = 6.4d+1*t32*t69
            t206 = 8.0d+1*t29*t69
            t207 = t206+t205-8.0d+0*t118
            t208 = 3.6d+1*t35*t71
            t209 = 4.8d+1*t32*t71
            t210 = t209+t208-6.0d+0*t119
            t211 = 1.6d+1*t38*t51
            t212 = 2.4d+1*t35*t51
            t213 = t212+t211-4.0d+0*t120
            t214 = t186*t85+t184*t72+t213*t70+5.28d+2*t133*t52*t7+t210*t
     1         68+t207*t66+t204*t64+t201*t62+t198*t60+t195*t58+t192*t56+
     2         t189*t54
            t215 = rhoa**2
            t216 = 1/t12**14
            t217 = 1/t5
            t218 = 1/t2**4
            t219 = rhob**2
            t220 = 1/t53**14
            t221 = 1/t48
            t222 = 1/t45**4
            t223 = -4.0d+0*(4.0d+0*t41+t183)-3.2d+1*t38-4.8d+1*t35-4.0d+
     1         0*t184
            fnc(iq) = (-9.305257363490997d-1*t47*t48*t74-9.3052573634909
     1         97d-1*t4*t43*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t43*t5*t83-1.2407009
     1         817987995d+0*t4*t43*t78-9.305257363490997d-1*t4*t5*t77)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t48*t74*t92-1.240700
     1         9817987995d+0*t47*t74*t87-9.305257363490997d-1*t47*t48*t8
     2         6)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t43*t
     1         5*t94*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t48*t
     1         74*t96*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t109*t4
     1         *t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t122*t4
     1         7*t48*wght
            Amat2(iq,D2_RA_RA) = (-2.481401963597599d+0*t43*t78*t83-1.86
     1         10514726981994d+0*t5*t77*t83-2.481401963597599d+0*t4*t77*
     2         t78-9.305257363490997d-1*t132*t43*t5-9.305257363490997d-1
     3         *t125*t4*t5-4.1356699393293317d-1*t126*t4*t43)*wght+Amat2
     4         (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-2.481401963597599d+0*t74*t87*t92-1.86
     1         10514726981994d+0*t48*t86*t92-2.481401963597599d+0*t47*t8
     2         6*t87-9.305257363490997d-1*t142*t48*t74-4.135669939329331
     3         7d-1*t136*t47*t74-9.305257363490997d-1*t135*t47*t48)*wght
     4         +Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817987995d+0*t43*t78*t94-9.
     1         305257363490997d-1*t5*t77*t94-9.305257363490997d-1*t144*t
     2         43*t5)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.2407009817987995d+0*t74*t87*t96-9.
     1         305257363490997d-1*t48*t86*t96-9.305257363490997d-1*t146*
     2         t48*t74)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t148*t43*t5*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t150*t48*t74*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t152*t5*t83-1.24
     1         07009817987995d+0*t152*t4*t78-9.305257363490997d-1*t151*t
     2         4*t5)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t154*t48*t92-1.2
     1         407009817987995d+0*t154*t47*t87-9.305257363490997d-1*t153
     2         *t47*t48)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t109*t5*t94*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t122*t48*t96*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t185*t4*t5*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t214*t47*t48*wght
            Amat3(iq,D3_RA_RA_RA) = (t126*(-1.2407009817987995d+0*t43*t8
     1         3-1.2407009817987995d+0*t4*t77)-7.444205890792797d+0*t77*
     2         t78*t83-2.791577209047299d+0*t125*t5*t83-9.30525736349099
     3         7d-1*t43*t5*(3.279403591520427d-3*t79*t80/rhoa**8.3333333
     4         33333334d+0-1.648227300332581d-1*gammaaa*t3/rhoa**5.66666
     5         6666666667d+0+3.7270905095708196d-8*gammaaa**4*t218/rhoa*
     6         *1.3666666666666666d+1-1.9702414090037443d-5*t127*t128/rh
     7         oa**11)+2.757113292886221d-1*t4*t43/t8-9.305257363490997d
     8         -1*t4*t5*(-1.1171548003727898d+5*t15*t7*t76*t78+5.0779763
     9         65330862d+4*t14*t15*t76*t78+5.077976365330863d+4*t13*t18*
     :         t7*t78+6.093571638397036d+4*t11*t123*t7*t78-9.23268430060
     ;         1568d+2*t40*t41*t78+9.232684300601568d+2*t10*t37*t41*t78+
     <         9.232684300601568d+2*t38*t40*t78+2.769805290180471d+3*t34
     =         *t38*t39*t78-8.309415870541413d+3*t34*t35*t39*t78+5.53961
     >         0580360942d+3*t32*t34*t39*t78-3.693073720240627d+3*t10*t3
     ?         7*t38*t78+2.769805290180471d+3*t10*t35*t37*t78+5.53961058
     @         0360942d+3*t31*t35*t36*t78-1.4772294880962508d+4*t31*t32*
     1         t36*t78+9.232684300601568d+3*t29*t31*t36*t78+9.2326843006
     2         01568d+3*t28*t32*t33*t78-2.308171075150392d+4*t28*t29*t33
     3         *t78+1.3849026450902355d+4*t26*t28*t33*t78+1.384902645090
     4         2355d+4*t25*t29*t30*t78-3.323766348216566d+4*t25*t26*t30*
     5         t78+1.9388637031263292d+4*t23*t25*t30*t78+1.9388637031263
     6         292d+4*t22*t26*t27*t78-4.5240153072947686d+4*t22*t23*t27*
     7         t78+2.585151604168439d+4*t20*t22*t27*t78+2.58515160416843
     8         9d+4*t19*t23*t24*t78-5.908917952385003d+4*t19*t20*t24*t78
     9         +3.3237663482165647d+4*t17*t19*t24*t78+3.3237663482165647
     :         d+4*t16*t20*t21*t78-7.478474283487271d+4*t16*t17*t21*t78+
     ;         4.154707935270705d+4*t14*t16*t21*t78+4.1547079352707056d+
     <         4*t13*t17*t18*t78-9.232684300601568d+4*t13*t14*t18*t78+3.
     =         713762858948263d+1*t11*t217*t7*t76-1.2729420016323426d+7*
     >         t18*t215*t7*t76-3.4716600044518436d+6*t17*t18*t215*t76+1.
     ?         1572200014839477d+7*t14*t18*t215*t76-3.713762858948263d+1
     @         *t13*t15*t217*t7-6.01754400771653d+6*t11*t215*t216*t7+3.4
     1         716600044518436d+6*t13*t21*t215*t7+1.5275304019588112d+7*
     2         t123*t15*t215*t7-3.37614805358933d+0*t217*t41*t42+3.37614
     3         805358933d+0*t10*t217*t40*t41+2.104036366334451d+4*t215*t
     4         37*t41-2.104036366334451d+4*t10*t215*t34*t41-6.7522961071
     5         7866d+0*t10*t217*t38*t40+6.75229610717866d+0*t217*t37*t38
     6         *t39-8.416145465337804d+4*t215*t31*t38*t39-1.012844416076
     7         7991d+1*t217*t35*t37*t39+3.787265459402013d+5*t215*t31*t3
     8         5*t39-5.049687279202681d+5*t215*t31*t32*t39+2.10403636633
     9         44508d+5*t215*t29*t31*t39-4.208072732668902d+4*t215*t37*t
     :         38+1.2624218198006704d+5*t10*t215*t34*t38+2.1040363663344
     ;         51d+4*t215*t35*t37+1.0128444160767991d+1*t217*t34*t35*t36
     <         -2.104036366334451d+5*t215*t28*t35*t36-1.350459221435732d
     =         +1*t217*t32*t34*t36+8.416145465337803d+5*t215*t28*t32*t36
     >         -1.0520181831672254d+6*t215*t28*t29*t36+4.208072732668902
     ?         d+5*t215*t26*t28*t36-1.8936327297010064d+5*t10*t215*t34*t
     @         35+8.416145465337804d+4*t10*t215*t32*t34+1.35045922143573
     1         2d+1*t217*t31*t32*t33-4.2080727326689016d+5*t215*t25*t32*
     2         t33-1.6880740267946648d+1*t217*t29*t31*t33+1.578027274750
     3         838d+6*t215*t25*t29*t33-1.893632729701006d+6*t215*t25*t26
     4         *t33+7.364127282170578d+5*t215*t23*t25*t33+1.688074026794
     5         6648d+1*t217*t28*t29*t30-7.364127282170578d+5*t215*t22*t2
     6         9*t30-2.0256888321535982d+1*t217*t26*t28*t30+2.6510858215
     7         814087d+6*t215*t22*t26*t30-3.0929334585116436d+6*t215*t22
     8         *t23*t30+1.1782603651472927d+6*t20*t215*t22*t30+2.0256888
     9         321535982d+1*t217*t25*t26*t27-1.1782603651472927d+6*t19*t
     :         215*t26*t27-2.3633036375125308d+1*t217*t23*t25*t27+4.1239
     ;         11278015524d+6*t19*t215*t23*t27-4.713041460589171d+6*t19*
     <         t20*t215*t27+1.767390547720939d+6*t17*t19*t215*t27+2.3633
     =         036375125308d+1*t217*t22*t23*t24-1.7673905477209387d+6*t1
     >         6*t215*t23*t24-2.700918442871464d+1*t20*t217*t22*t24+6.05
     ?         9624735043219d+6*t16*t20*t215*t24-6.817077826923622d+6*t1
     @         6*t17*t215*t24+2.524843639601341d+6*t14*t16*t215*t24+2.70
     1         0918442871464d+1*t19*t20*t21*t217-3.0385332482303967d+1*t
     2         17*t19*t21*t217+3.0385332482303967d+1*t16*t17*t18*t217-3.
     3         3761480535893296d+1*t14*t16*t18*t217+3.3761480535893296d+
     4         1*t13*t14*t15*t217-2.524843639601341d+6*t13*t20*t21*t215+
     5         8.521347283654526d+6*t13*t17*t21*t215-9.468163648505028d+
     6         6*t13*t14*t21*t215-4.628880005935792d+6*t123*t14*t15*t215
     7         )-3.7221029453963983d+0*t132*t43*t78-3.7221029453963983d+
     8         0*t125*t4*t78-2.791577209047299d+0*t132*t5*t77)*wght+Amat
     9         3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t136*(-1.2407009817987995d+0*t74*t9
     1         2-1.2407009817987995d+0*t47*t86)-7.444205890792797d+0*t86
     2         *t87*t92-2.791577209047299d+0*t135*t48*t92-9.305257363490
     3         997d-1*t48*t74*(3.279403591520427d-3*t88*t89/rhob**8.3333
     4         33333333334d+0-1.648227300332581d-1*gammabb*t46/rhob**5.6
     5         66666666666667d+0+3.7270905095708196d-8*gammabb**4*t222/r
     6         hob**1.3666666666666666d+1-1.9702414090037443d-5*t137*t13
     7         8/rhob**11)-9.305257363490997d-1*t47*t48*(-1.117154800372
     8         7898d+5*t55*t7*t85*t87+5.077976365330862d+4*t14*t55*t85*t
     9         87-9.232684300601568d+2*t41*t72*t87+9.232684300601568d+2*
     :         t38*t72*t87+2.769805290180471d+3*t38*t68*t71*t87-8.309415
     ;         870541413d+3*t35*t68*t71*t87+5.539610580360942d+3*t32*t68
     <         *t71*t87+9.232684300601568d+2*t41*t51*t70*t87-3.693073720
     =         240627d+3*t38*t51*t70*t87+2.769805290180471d+3*t35*t51*t7
     >         0*t87+5.077976365330863d+4*t54*t57*t7*t87+6.0935716383970
     ?         36d+4*t133*t52*t7*t87+5.539610580360942d+3*t35*t66*t69*t8
     @         7-1.4772294880962508d+4*t32*t66*t69*t87+9.232684300601568
     1         d+3*t29*t66*t69*t87+9.232684300601568d+3*t32*t64*t67*t87-
     2         2.308171075150392d+4*t29*t64*t67*t87+1.3849026450902355d+
     3         4*t26*t64*t67*t87+1.3849026450902355d+4*t29*t62*t65*t87-3
     4         .323766348216566d+4*t26*t62*t65*t87+1.9388637031263292d+4
     5         *t23*t62*t65*t87+1.9388637031263292d+4*t26*t60*t63*t87-4.
     6         5240153072947686d+4*t23*t60*t63*t87+2.585151604168439d+4*
     7         t20*t60*t63*t87+2.585151604168439d+4*t23*t58*t61*t87-5.90
     8         8917952385003d+4*t20*t58*t61*t87+3.3237663482165647d+4*t1
     9         7*t58*t61*t87+3.3237663482165647d+4*t20*t56*t59*t87-7.478
     :         474283487271d+4*t17*t56*t59*t87+4.154707935270705d+4*t14*
     ;         t56*t59*t87+4.1547079352707056d+4*t17*t54*t57*t87-9.23268
     <         4300601568d+4*t14*t54*t57*t87-1.2729420016323426d+7*t219*
     =         t57*t7*t85+3.713762858948263d+1*t221*t52*t7*t85-3.4716600
     >         044518436d+6*t17*t219*t57*t85+1.1572200014839477d+7*t14*t
     ?         219*t57*t85-3.37614805358933d+0*t221*t41*t73+3.3761480535
     @         8933d+0*t221*t41*t51*t72-6.75229610717866d+0*t221*t38*t51
     1         *t72+6.75229610717866d+0*t221*t38*t70*t71-1.0128444160767
     2         991d+1*t221*t35*t70*t71-8.416145465337804d+4*t219*t38*t66
     3         *t71+3.787265459402013d+5*t219*t35*t66*t71-5.049687279202
     4         681d+5*t219*t32*t66*t71+2.1040363663344508d+5*t219*t29*t6
     5         6*t71+2.104036366334451d+4*t219*t41*t70-4.208072732668902
     6         d+4*t219*t38*t70+2.104036366334451d+4*t219*t35*t70+3.4716
     7         600044518436d+6*t219*t54*t59*t7-3.713762858948263d+1*t221
     8         *t54*t55*t7+1.5275304019588112d+7*t133*t219*t55*t7-6.0175
     9         4400771653d+6*t219*t220*t52*t7+1.0128444160767991d+1*t221
     :         *t35*t68*t69-1.350459221435732d+1*t221*t32*t68*t69-2.1040
     ;         36366334451d+5*t219*t35*t64*t69+8.416145465337803d+5*t219
     <         *t32*t64*t69-1.0520181831672254d+6*t219*t29*t64*t69+4.208
     =         072732668902d+5*t219*t26*t64*t69-2.104036366334451d+4*t21
     >         9*t41*t51*t68+1.2624218198006704d+5*t219*t38*t51*t68-1.89
     ?         36327297010064d+5*t219*t35*t51*t68+8.416145465337804d+4*t
     @         219*t32*t51*t68+1.350459221435732d+1*t221*t32*t66*t67-1.6
     1         880740267946648d+1*t221*t29*t66*t67-4.2080727326689016d+5
     2         *t219*t32*t62*t67+1.578027274750838d+6*t219*t29*t62*t67-1
     3         .893632729701006d+6*t219*t26*t62*t67+7.364127282170578d+5
     4         *t219*t23*t62*t67+1.6880740267946648d+1*t221*t29*t64*t65-
     5         2.0256888321535982d+1*t221*t26*t64*t65-7.364127282170578d
     6         +5*t219*t29*t60*t65+2.6510858215814087d+6*t219*t26*t60*t6
     7         5-3.0929334585116436d+6*t219*t23*t60*t65+1.17826036514729
     8         27d+6*t20*t219*t60*t65+2.0256888321535982d+1*t221*t26*t62
     9         *t63-2.3633036375125308d+1*t221*t23*t62*t63-1.17826036514
     :         72927d+6*t219*t26*t58*t63+4.123911278015524d+6*t219*t23*t
     ;         58*t63-4.713041460589171d+6*t20*t219*t58*t63+1.7673905477
     <         20939d+6*t17*t219*t58*t63+2.3633036375125308d+1*t221*t23*
     =         t60*t61-2.700918442871464d+1*t20*t221*t60*t61-1.767390547
     >         7209387d+6*t219*t23*t56*t61+6.059624735043219d+6*t20*t219
     ?         *t56*t61-6.817077826923622d+6*t17*t219*t56*t61+2.52484363
     @         9601341d+6*t14*t219*t56*t61+2.700918442871464d+1*t20*t221
     1         *t58*t59-3.0385332482303967d+1*t17*t221*t58*t59-2.5248436
     2         39601341d+6*t20*t219*t54*t59+8.521347283654526d+6*t17*t21
     3         9*t54*t59-9.468163648505028d+6*t14*t219*t54*t59+3.0385332
     4         482303967d+1*t17*t221*t56*t57-3.3761480535893296d+1*t14*t
     5         221*t56*t57+3.3761480535893296d+1*t14*t221*t54*t55-4.6288
     6         80005935792d+6*t133*t14*t219*t55)-3.7221029453963983d+0*t
     7         142*t74*t87-3.7221029453963983d+0*t135*t47*t87-2.79157720
     8         9047299d+0*t142*t48*t86+2.757113292886221d-1*t47*t74/t49)
     9         *wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-2.481401963597599d+0*t77*t78*t94-
     1         9.305257363490997d-1*t125*t5*t94-4.1356699393293317d-1*t1
     2         26*t43*t94-9.305257363490997d-1*t43*t5*(-9.37659384672262
     3         8d-4*gammaaa*t130*t80+6.610678411788878d-6*t128*t129*t79+
     4         3.531915643569816d-2*t131*t3-1.3976589410890572d-8*t127*t
     5         218/rhoa**1.2666666666666665d+1)-2.481401963597599d+0*t14
     6         4*t43*t78-1.8610514726981994d+0*t144*t5*t77)*wght+Cmat3(i
     7         q,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-2.481401963597599d+0*t86*t87*t96-
     1         4.1356699393293317d-1*t136*t74*t96-9.305257363490997d-1*t
     2         135*t48*t96-9.305257363490997d-1*t48*t74*(-9.376593846722
     3         628d-4*gammabb*t140*t89+6.610678411788878d-6*t138*t139*t8
     4         8+3.531915643569816d-2*t141*t46-1.3976589410890572d-8*t13
     5         7*t222/rhob**1.2666666666666665d+1)-2.481401963597599d+0*
     6         t146*t74*t87-1.8610514726981994d+0*t146*t48*t86)*wght+Cma
     7         t3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t43*t5*(1.7
     1         31063479394947d-4*t80*t81+5.241221029083964d-9*t218*t79/r
     2         hoa**1.1666666666666665d+1-1.944317179937905d-6*gammaaa*t
     3         128*t143)-1.2407009817987995d+0*t148*t43*t78-9.3052573634
     4         90997d-1*t148*t5*t77)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = Cmat3(iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = Cmat3(iq,D3_RB_GAA_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = (-9.305257363490997d-1*t48*t74*(1.
     1         731063479394947d-4*t89*t90+5.241221029083964d-9*t222*t88/
     2         rhob**1.1666666666666665d+1-1.944317179937905d-6*gammabb*
     3         t138*t145)-1.2407009817987995d+0*t150*t74*t87-9.305257363
     4         490997d-1*t150*t48*t86)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*(4.374713654860287d-7*t128*t147-1.9654578859
     2         06487d-9*gammaaa*t218/rhoa**1.0666666666666666d+1)*t43*t5
     3         *wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490997d-1*(4.374713654860287d-7*t138*t149-1.9654578859
     2         06487d-9*gammabb*t222/rhob**1.0666666666666666d+1)*t48*t7
     3         4*wght
            Mmat3(iq,D3_RA_RA_TA) = (-2.481401963597599d+0*t152*t78*t83-
     1         1.8610514726981994d+0*t151*t5*t83-2.481401963597599d+0*t1
     2         51*t4*t78-9.305257363490997d-1*t4*t5*(5.585774001863949d+
     3         5*t18*t5*t7*t76-4.5701787287977763d+5*t17*t18*t5*t76+5.07
     4         7976365330863d+5*t14*t18*t5*t76+2.2282577153689576d+3*t12
     5         4*t14*t15*t76-7.921643129916147d+5*t11*t216*t5*t7-4.57017
     6         87287977763d+5*t13*t21*t5*t7+6.702928802236738d+5*t123*t1
     7         5*t5*t7-2.228257715368958d+3*t124*t13*t18*t7+2.6739092584
     8         427493d+3*t11*t123*t124*t7+9.232684300601569d+2*t37*t41*t
     9         5-2.769805290180471d+3*t10*t34*t41*t5-1.1079221160721883d
     :         +4*t31*t38*t39*t5+1.6618831741082832d+4*t31*t35*t39*t5+2.
     ;         2158442321443766d+4*t31*t32*t39*t5-2.7698052901804704d+4*
     <         t29*t31*t39*t5+1.846536860120314d+3*t37*t38*t5+5.53961058
     =         0360942d+3*t10*t34*t38*t5-2.769805290180471d+3*t35*t37*t5
     >         -2.769805290180471d+4*t28*t35*t36*t5+3.693073720240629d+4
     ?         *t28*t32*t36*t5+4.616342150300784d+4*t28*t29*t36*t5-5.539
     @         610580360942d+4*t26*t28*t36*t5+8.309415870541416d+3*t10*t
     1         34*t35*t5-1.1079221160721883d+4*t10*t32*t34*t5-5.53961058
     2         0360941d+4*t25*t32*t33*t5+6.924513225451174d+4*t25*t29*t3
     3         3*t5+8.309415870541416d+4*t25*t26*t33*t5-9.69431851563164
     4         8d+4*t23*t25*t33*t5-9.694318515631648d+4*t22*t29*t30*t5+1
     5         .1633182218757979d+5*t22*t26*t30*t5+1.3572045921884307d+5
     6         *t22*t23*t30*t5-1.5510909625010635d+5*t20*t22*t30*t5-1.55
     7         10909625010638d+5*t19*t26*t27*t5+1.809606122917908d+5*t19
     8         *t23*t27*t5+2.0681212833347515d+5*t19*t20*t27*t5-2.326636
     9         4437515952d+5*t17*t19*t27*t5-2.3266364437515955d+5*t16*t2
     :         3*t24*t5+2.659013078573252d+5*t16*t20*t24*t5+2.9913897133
     ;         949086d+5*t16*t17*t24*t5-3.3237663482165647d+5*t14*t16*t2
     <         4*t5-3.323766348216565d+5*t13*t20*t21*t5+3.73923714174363
     =         54d+5*t13*t17*t21*t5+4.1547079352707056d+5*t13*t14*t21*t5
     >         -6.093571638397035d+5*t123*t14*t15*t5+4.0513776643071964d
     ?         +1*t10*t124*t37*t41-4.0513776643071964d+1*t124*t38*t40+1.
     @         2154132992921589d+2*t124*t34*t38*t39-2.4308265985843178d+
     1         2*t124*t32*t34*t39-1.2154132992921589d+2*t10*t124*t35*t37
     2         +2.4308265985843178d+2*t124*t31*t35*t36-4.051377664307196
     3         d+2*t124*t29*t31*t36+4.0513776643071964d+2*t124*t28*t32*t
     4         33-6.077066496460795d+2*t124*t26*t28*t33+6.07706649646079
     5         3d+2*t124*t25*t29*t30-8.507893095045111d+2*t124*t23*t25*t
     6         30+8.507893095045113d+2*t124*t22*t26*t27-1.13438574600601
     7         49d+3*t124*t20*t22*t27+1.1343857460060149d+3*t124*t19*t23
     8         *t24-1.4584959591505905d+3*t124*t17*t19*t24+1.45849595915
     9         05908d+3*t124*t16*t20*t21-1.8231199489382383d+3*t124*t14*
     :         t16*t21+1.8231199489382383d+3*t124*t13*t17*t18)-9.3052573
     ;         63490997d-1*t132*t152*t5-4.1356699393293317d-1*t126*t152*
     <         t4)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (-2.481401963597599d+0*t154*t87*t92-
     1         1.8610514726981994d+0*t153*t48*t92-2.481401963597599d+0*t
     2         153*t47*t87-9.305257363490997d-1*t47*t48*(5.5857740018639
     3         49d+5*t48*t57*t7*t85-4.5701787287977763d+5*t17*t48*t57*t8
     4         5+5.077976365330863d+5*t14*t48*t57*t85+2.2282577153689576
     5         d+3*t134*t14*t55*t85-4.0513776643071964d+1*t134*t38*t72+1
     6         .2154132992921589d+2*t134*t38*t68*t71-2.4308265985843178d
     7         +2*t134*t32*t68*t71-1.1079221160721883d+4*t38*t48*t66*t71
     8         +1.6618831741082832d+4*t35*t48*t66*t71+2.2158442321443766
     9         d+4*t32*t48*t66*t71-2.7698052901804704d+4*t29*t48*t66*t71
     :         +4.0513776643071964d+1*t134*t41*t51*t70-1.215413299292158
     ;         9d+2*t134*t35*t51*t70+9.232684300601569d+2*t41*t48*t70+1.
     <         846536860120314d+3*t38*t48*t70-2.769805290180471d+3*t35*t
     =         48*t70-4.5701787287977763d+5*t48*t54*t59*t7-2.22825771536
     >         8958d+3*t134*t54*t57*t7+6.702928802236738d+5*t133*t48*t55
     ?         *t7-7.921643129916147d+5*t220*t48*t52*t7+2.67390925844274
     @         93d+3*t133*t134*t52*t7+2.4308265985843178d+2*t134*t35*t66
     1         *t69-4.051377664307196d+2*t134*t29*t66*t69-2.769805290180
     2         471d+4*t35*t48*t64*t69+3.693073720240629d+4*t32*t48*t64*t
     3         69+4.616342150300784d+4*t29*t48*t64*t69-5.539610580360942
     4         d+4*t26*t48*t64*t69-2.769805290180471d+3*t41*t48*t51*t68+
     5         5.539610580360942d+3*t38*t48*t51*t68+8.309415870541416d+3
     6         *t35*t48*t51*t68-1.1079221160721883d+4*t32*t48*t51*t68+4.
     7         0513776643071964d+2*t134*t32*t64*t67-6.077066496460795d+2
     8         *t134*t26*t64*t67-5.539610580360941d+4*t32*t48*t62*t67+6.
     9         924513225451174d+4*t29*t48*t62*t67+8.309415870541416d+4*t
     :         26*t48*t62*t67-9.694318515631648d+4*t23*t48*t62*t67+6.077
     ;         066496460793d+2*t134*t29*t62*t65-8.507893095045111d+2*t13
     <         4*t23*t62*t65-9.694318515631648d+4*t29*t48*t60*t65+1.1633
     =         182218757979d+5*t26*t48*t60*t65+1.3572045921884307d+5*t23
     >         *t48*t60*t65-1.5510909625010635d+5*t20*t48*t60*t65+8.5078
     ?         93095045113d+2*t134*t26*t60*t63-1.1343857460060149d+3*t13
     @         4*t20*t60*t63-1.5510909625010638d+5*t26*t48*t58*t63+1.809
     1         606122917908d+5*t23*t48*t58*t63+2.0681212833347515d+5*t20
     2         *t48*t58*t63-2.3266364437515952d+5*t17*t48*t58*t63+1.1343
     3         857460060149d+3*t134*t23*t58*t61-1.4584959591505905d+3*t1
     4         34*t17*t58*t61-2.3266364437515955d+5*t23*t48*t56*t61+2.65
     5         9013078573252d+5*t20*t48*t56*t61+2.9913897133949086d+5*t1
     6         7*t48*t56*t61-3.3237663482165647d+5*t14*t48*t56*t61+1.458
     7         4959591505908d+3*t134*t20*t56*t59-1.8231199489382383d+3*t
     8         134*t14*t56*t59-3.323766348216565d+5*t20*t48*t54*t59+3.73
     9         92371417436354d+5*t17*t48*t54*t59+4.1547079352707056d+5*t
     :         14*t48*t54*t59+1.8231199489382383d+3*t134*t17*t54*t57-6.0
     ;         93571638397035d+5*t133*t14*t48*t55)-9.305257363490997d-1*
     <         t142*t154*t48-4.1356699393293317d-1*t136*t154*t47)*wght+M
     =         mat3(iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RA_GAA_TA) = (-1.2407009817987995d+0*t152*t78*t9
     1         4-9.305257363490997d-1*t151*t5*t94-9.305257363490997d-1*t
     2         144*t152*t5)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_RB_GAA_TA) = Mmat3(iq,D3_RB_GAA_TA)
            Mmat3(iq,D3_RB_GAA_TB) = Mmat3(iq,D3_RB_GAA_TB)
            Mmat3(iq,D3_RB_GAB_TA) = Mmat3(iq,D3_RB_GAB_TA)
            Mmat3(iq,D3_RB_GAB_TB) = Mmat3(iq,D3_RB_GAB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = (-1.2407009817987995d+0*t154*t87*t9
     1         6-9.305257363490997d-1*t153*t48*t96-9.305257363490997d-1*
     2         t146*t154*t48)*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)-9.305257363490
     1         997d-1*t109*t148*t5*wght
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAA_GBB_TB) = Mmat3(iq,D3_GAA_GBB_TB)-9.30525736
     1         3490997d-1*t122*t150*t48*wght
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)-9.30525736
     1         3490997d-1*t185*t5*t94*wght
            Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
            Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
            Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)-9.3052573634
     1         90997d-1*t214*t48*t96*wght
            Mmat3(iq,D3_GAB_TA_TA) = Mmat3(iq,D3_GAB_TA_TA)-9.3052573634
     1         90997d-1*t4*t5*((-9.68d+3*t18*t7-2.2d+1*t158-2.2d+1*(t157
     2         +t156))*t76+t13*(-7.92d+3*t21*t7-7.2d+3*t14*t21-2.0d+1*t1
     3         61-2.0d+1*(t160+t159))+t123*(-1.1616d+4*t15*t7-2.4d+1*t15
     4         5)-1.3728d+4*t11*t216*t7+t31*(-3.84d+2*t32*t39-4.8d+2*t29
     5         *t39-8.0d+0*t179-8.0d+0*(t178+t177))+t223*t37+t28*(-8.0d+
     6         2*t29*t36-9.6d+2*t26*t36-1.0d+1*t176-1.0d+1*(t175+t174))+
     7         t34*(-1.44d+2*t10*t35-1.92d+2*t10*t32-6.0d+0*t182-6.0d+0*
     8         (t181+t180))+t25*(-1.44d+3*t26*t33-1.68d+3*t23*t33-1.2d+1
     9         *t173-1.2d+1*(t172+t171))+t22*(-2.352d+3*t23*t30-2.688d+3
     :         *t20*t30-1.4d+1*t170-1.4d+1*(t169+t168))+t19*(-3.584d+3*t
     ;         20*t27-4.032d+3*t17*t27-1.6d+1*t167-1.6d+1*(t166+t165))+t
     <         16*(-5.184d+3*t17*t24-5.76d+3*t14*t24-1.8d+1*t164-1.8d+1*
     =         (t163+t162)))*wght
            Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
            Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)-9.3052573634
     1         90997d-1*t47*t48*((-9.68d+3*t57*t7-2.2d+1*t189-2.2d+1*(t1
     2         88+t187))*t85+t66*(-3.84d+2*t32*t71-4.8d+2*t29*t71-8.0d+0
     3         *t210-8.0d+0*(t209+t208))+t223*t70+t54*(-7.92d+3*t59*t7-7
     4         .2d+3*t14*t59-2.0d+1*t192-2.0d+1*(t191+t190))+t133*(-1.16
     5         16d+4*t55*t7-2.4d+1*t186)-1.3728d+4*t220*t52*t7+t64*(-8.0
     6         d+2*t29*t69-9.6d+2*t26*t69-1.0d+1*t207-1.0d+1*(t206+t205)
     7         )+(-1.44d+2*t35*t51-1.92d+2*t32*t51-6.0d+0*t213-6.0d+0*(t
     8         212+t211))*t68+t62*(-1.44d+3*t26*t67-1.68d+3*t23*t67-1.2d
     9         +1*t204-1.2d+1*(t203+t202))+t60*(-2.352d+3*t23*t65-2.688d
     :         +3*t20*t65-1.4d+1*t201-1.4d+1*(t200+t199))+t58*(-3.584d+3
     ;         *t20*t63-4.032d+3*t17*t63-1.6d+1*t198-1.6d+1*(t197+t196))
     <         +t56*(-5.184d+3*t17*t61-5.76d+3*t14*t61-1.8d+1*t195-1.8d+
     =         1*(t194+t193)))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(11)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(10)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(9)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(8)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(7)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(6)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(5)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(4)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(3)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(2)
            t41 = 1/t11
            t42 = t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+t31*t32*
     1         t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t21+t16*t
     2         17*t18+t13*t14*t15+param(1)
            t43 = rhoa**6.666666666666666d-1
            t44 = 1/t11**12
            t45 = -1.5192666241151986d+1*t39*t40*t43*t9+3.03853324823039
     1         7d+1*t37*t39*t43*t9-1.6711932865267184d+2*t10*t43*t44*t6+
     2         1.6711932865267184d+2*t12*t14*t43*t6+1.5192666241151986d+
     3         1*t40*t41*t43-3.038533248230397d+1*t36*t37*t38*t43+4.5577
     4         998723455965d+1*t34*t36*t38*t43-4.5577998723455965d+1*t33
     5         *t34*t35*t43+6.077066496460795d+1*t31*t33*t35*t43-6.07706
     6         6496460795d+1*t30*t31*t32*t43+7.596333120575993d+1*t28*t3
     7         0*t32*t43-7.596333120575993d+1*t27*t28*t29*t43+9.11559974
     8         4691193d+1*t25*t27*t29*t43-9.115599744691193d+1*t24*t25*t
     9         26*t43+1.0634866368806391d+2*t22*t24*t26*t43-1.0634866368
     :         806391d+2*t21*t22*t23*t43+1.2154132992921589d+2*t19*t21*t
     ;         23*t43-1.2154132992921589d+2*t18*t19*t20*t43+1.3673399617
     <         036788d+2*t16*t18*t20*t43-1.3673399617036788d+2*t15*t16*t
     =         17*t43+1.5192666241151986d+2*t13*t15*t17*t43-1.5192666241
     >         151986d+2*t12*t13*t14*t43
            t46 = rhoa**3.333333333333333d-1
            t47 = gammaaa**2
            t48 = 1/t2**2
            t49 = 1/rhoa**6.333333333333333d+0
            t50 = 1/rhoa**3.6666666666666664d+0
            t51 = 4.327658698487368d-5*t47*t48*t49-9.632497209735863d-3*
     1         gammaaa*t3*t50
            t52 = 1/rhoa**5.333333333333333d+0
            t53 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t48*t52
            t54 = -2.2d+1*t10*t44*t6
            t55 = -2.2d+1*t14*t6-2.0d+1*t13*t14
            t56 = -1.8d+1*t16*t17-2.0d+1*t13*t17
            t57 = -1.6d+1*t19*t20-1.8d+1*t16*t20
            t58 = -1.4d+1*t22*t23-1.6d+1*t19*t23
            t59 = -1.2d+1*t25*t26-1.4d+1*t22*t26
            t60 = -1.0d+1*t28*t29-1.2d+1*t25*t29
            t61 = -8.0d+0*t31*t32-1.0d+1*t28*t32
            t62 = -6.0d+0*t34*t35-8.0d+0*t31*t35
            t63 = -4.0d+0*t37*t38-6.0d+0*t34*t38
            t64 = -2.0d+0*t40*t9-4.0d+0*t37*t9
            t65 = -2.0d+0*t40*t41
            t66 = t65+t39*t64+t36*t63+t33*t62+t30*t61+t27*t60+t24*t59+t2
     1         1*t58+t18*t57+t15*t56+t12*t55+t54
            t67 = 1/t11**13
            t68 = 1/t46
            t69 = -1.0128444160767991d+1*t39*t40*t68*t9+2.02568883215359
     1         82d+1*t37*t39*t68*t9+4.6163421503007845d+2*t36*t40*t5*t9-
     2         1.846536860120314d+3*t36*t37*t5*t9+1.3849026450902355d+3*
     3         t34*t36*t5*t9-1.1141288576844789d+2*t10*t44*t6*t68+1.1141
     4         288576844789d+2*t12*t14*t6*t68+1.0128444160767991d+1*t40*
     5         t41*t68-2.0256888321535982d+1*t36*t37*t38*t68+3.038533248
     6         230397d+1*t34*t36*t38*t68-3.038533248230397d+1*t33*t34*t3
     7         5*t68+4.0513776643071964d+1*t31*t33*t35*t68-4.05137766430
     8         71964d+1*t30*t31*t32*t68+5.064222080383994d+1*t28*t30*t32
     9         *t68-5.064222080383994d+1*t27*t28*t29*t68+6.0770664964607
     :         95d+1*t25*t27*t29*t68-6.077066496460795d+1*t24*t25*t26*t6
     ;         8+7.089910912537593d+1*t22*t24*t26*t68-7.089910912537593d
     <         +1*t21*t22*t23*t68+8.102755328614393d+1*t19*t21*t23*t68-8
     =         .102755328614393d+1*t18*t19*t20*t68+9.115599744691191d+1*
     >         t16*t18*t20*t68-9.115599744691191d+1*t15*t16*t17*t68+1.01
     ?         2844416076799d+2*t13*t15*t17*t68-1.012844416076799d+2*t12
     @         *t13*t14*t68+3.046785819198518d+4*t10*t5*t6*t67-5.5857740
     1         01863949d+4*t14*t44*t5*t6+2.5389881826654315d+4*t12*t17*t
     2         5*t6+2.5389881826654315d+4*t13*t14*t44*t5-4.6163421503007
     3         845d+2*t39*t40*t5+4.6163421503007845d+2*t37*t39*t5+1.3849
     4         026450902355d+3*t33*t37*t38*t5-4.1547079352707067d+3*t33*
     5         t34*t38*t5+2.769805290180471d+3*t31*t33*t38*t5+2.76980529
     6         0180471d+3*t30*t34*t35*t5-7.386147440481255d+3*t30*t31*t3
     7         5*t5+4.616342150300784d+3*t28*t30*t35*t5+4.61634215030078
     8         4d+3*t27*t31*t32*t5-1.154085537575196d+4*t27*t28*t32*t5+6
     9         .924513225451178d+3*t25*t27*t32*t5+6.924513225451178d+3*t
     :         24*t28*t29*t5-1.661883174108283d+4*t24*t25*t29*t5+9.69431
     ;         8515631648d+3*t22*t24*t29*t5+9.694318515631648d+3*t21*t25
     <         *t26*t5-2.2620076536473846d+4*t21*t22*t26*t5+1.2925758020
     =         842196d+4*t19*t21*t26*t5+1.2925758020842196d+4*t18*t22*t2
     >         3*t5-2.954458976192502d+4*t18*t19*t23*t5+1.66188317410828
     ?         23d+4*t16*t18*t23*t5+1.6618831741082823d+4*t15*t19*t20*t5
     @         -3.7392371417436354d+4*t15*t16*t20*t5+2.0773539676353528d
     1         +4*t13*t15*t20*t5+2.0773539676353528d+4*t12*t16*t17*t5-4.
     2         616342150300784d+4*t12*t13*t17*t5
            t70 = 1/t43
            t71 = gammaaa**3
            t72 = 1/t2**3
            t73 = 1/rhoa**10
            t74 = 1/rhoa**7.333333333333333d+0
            t75 = 1/rhoa**4.666666666666667d+0
            t76 = 3.531915643569816d-2*gammaaa*t3*t75-3.8948928286386303
     1         d-4*t47*t48*t74+1.0369691626335495d-6*t71*t72*t73
            t77 = 1/rhoa**9
            t78 = -3.8886343598758105d-7*t47*t72*t77-9.632497209735863d-
     1         3*t3*t50+1.29829760954621d-4*gammaaa*t48*t49
            t79 = 1/rhoa**8
            t80 = 1.4582378849534292d-7*gammaaa*t72*t79-3.24574402386552
     1         55d-5*t48*t52
            t81 = 6.077066496460795d+1*t36*t40*t43*t9-1.8231199489382385
     1         d+2*t34*t36*t43*t9+4.0108638876641245d+3*t10*t43*t6*t67-3
     2         .342386573053437d+3*t12*t17*t43*t6+3.342386573053437d+3*t
     3         13*t14*t43*t44-6.077066496460795d+1*t37*t39*t43+1.8231199
     4         489382385d+2*t33*t37*t38*t43-3.646239897876477d+2*t31*t33
     5         *t38*t43+3.646239897876477d+2*t30*t34*t35*t43-6.077066496
     6         460795d+2*t28*t30*t35*t43+6.077066496460795d+2*t27*t31*t3
     7         2*t43-9.115599744691192d+2*t25*t27*t32*t43+9.115599744691
     8         192d+2*t24*t28*t29*t43-1.276183964256767d+3*t22*t24*t29*t
     9         43+1.276183964256767d+3*t21*t25*t26*t43-1.701578619009022
     :         6d+3*t19*t21*t26*t43+1.7015786190090226d+3*t18*t22*t23*t4
     ;         3-2.187743938725886d+3*t16*t18*t23*t43+2.187743938725886d
     <         +3*t15*t19*t20*t43-2.734679923407357d+3*t13*t15*t20*t43+2
     =         .734679923407357d+3*t12*t16*t17*t43
            t82 = -2.0d+0*t39*t40*t9-4.0d+0*t37*t39*t9+t65-2.2d+1*t12*t1
     1         4*t6+t54-4.0d+0*t36*t37*t38-6.0d+0*t34*t36*t38-6.0d+0*t33
     2         *t34*t35-8.0d+0*t31*t33*t35-8.0d+0*t30*t31*t32-1.0d+1*t28
     3         *t30*t32-1.0d+1*t27*t28*t29-1.2d+1*t25*t27*t29-1.2d+1*t24
     4         *t25*t26-1.4d+1*t22*t24*t26-1.4d+1*t21*t22*t23-1.6d+1*t19
     5         *t21*t23-1.6d+1*t18*t19*t20-1.8d+1*t16*t18*t20-1.8d+1*t15
     6         *t16*t17-2.0d+1*t13*t15*t17-2.0d+1*t12*t13*t14
            t83 = 5.28d+2*t10*t6*t67
            t84 = 4.84d+2*t14*t6-2.2d+1*t55
            t85 = 4.0d+2*t13*t17
            t86 = 4.4d+2*t17*t6
            t87 = t86+t85-2.0d+1*t56
            t88 = 3.24d+2*t16*t20
            t89 = 3.6d+2*t13*t20
            t90 = t89+t88-1.8d+1*t57
            t91 = 2.56d+2*t19*t23
            t92 = 2.88d+2*t16*t23
            t93 = t92+t91-1.6d+1*t58
            t94 = 1.96d+2*t22*t26
            t95 = 2.24d+2*t19*t26
            t96 = t95+t94-1.4d+1*t59
            t97 = 1.44d+2*t25*t29
            t98 = 1.68d+2*t22*t29
            t99 = t98+t97-1.2d+1*t60
            t100 = 1.0d+2*t28*t32
            t101 = 1.2d+2*t25*t32
            t102 = -1.0d+1*t61+t101+t100
            t103 = 6.4d+1*t31*t35
            t104 = 8.0d+1*t28*t35
            t105 = -8.0d+0*t62+t104+t103
            t106 = 3.6d+1*t34*t38
            t107 = 4.8d+1*t31*t38
            t108 = -6.0d+0*t63+t107+t106
            t109 = 1.6d+1*t37*t9
            t110 = 2.4d+1*t34*t9
            t111 = -4.0d+0*t64+t110+t109
            t112 = 8.0d+0*t37
            t113 = 8.0d+0*t40+t112
            t114 = t24*t99+t21*t96+t18*t93+t15*t90+t12*t87+t44*t84+t83+t
     1         113*t39+t111*t36+t108*t33+t105*t30+t102*t27
            t115 = rhoa**2
            t116 = 1/t11**14
            t117 = 1/t5
            t118 = 1/t2**4
            t119 = 8.0d+0*t36*t40*t9+3.2d+1*t36*t37*t9+2.4d+1*t34*t36*t9
     1         +t83+9.68d+2*t14*t44*t6+4.4d+2*t12*t17*t6+4.4d+2*t13*t14*
     2         t44+8.0d+0*t39*t40+8.0d+0*t37*t39+2.4d+1*t33*t37*t38+7.2d
     3         +1*t33*t34*t38+4.8d+1*t31*t33*t38+4.8d+1*t30*t34*t35+1.28
     4         d+2*t30*t31*t35+8.0d+1*t28*t30*t35+8.0d+1*t27*t31*t32+2.0
     5         d+2*t27*t28*t32+1.2d+2*t25*t27*t32+1.2d+2*t24*t28*t29+2.8
     6         8d+2*t24*t25*t29+1.68d+2*t22*t24*t29+1.68d+2*t21*t25*t26+
     7         3.92d+2*t21*t22*t26+2.24d+2*t19*t21*t26+2.24d+2*t18*t22*t
     8         23+5.12d+2*t18*t19*t23+2.88d+2*t16*t18*t23+2.88d+2*t15*t1
     9         9*t20+6.48d+2*t15*t16*t20+3.6d+2*t13*t15*t20+3.6d+2*t12*t
     :         16*t17+8.0d+2*t12*t13*t17
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t42*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t42*t5*t51*wght-9.305
     1         257363490997d-1*t4*t45*t5*wght-1.2407009817987995d+0*t4*t
     2         42*t46*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t42*t
     1         5*t53*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         66*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = -9.305257363490997d-1*t42*t5*t76*wght-4
     1         .1356699393293317d-1*t4*t42*t70*wght-9.305257363490997d-1
     2         *t4*t5*t69*wght-1.8610514726981994d+0*t45*t5*t51*wght-2.4
     3         81401963597599d+0*t42*t46*t51*wght-2.481401963597599d+0*t
     4         4*t45*t46*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490997d-1*t42*t5*t78*wght-
     1         9.305257363490997d-1*t45*t5*t53*wght-1.2407009817987995d+
     2         0*t42*t46*t53*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t42*t5*t80*wght
            Mmat2(iq,D2_RA_TA) = -9.305257363490997d-1*t5*t51*t82*wght-1
     1         .2407009817987995d+0*t4*t46*t82*wght-9.305257363490997d-1
     2         *t4*t5*t81*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t53*t66*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t114*t4*t5*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = t70*(-1.2407009817987995d+0*t42*t51*
     1         wght-1.2407009817987995d+0*t4*t45*wght)-9.305257363490997
     2         d-1*t4*t5*(9.232684300601568d+2*t36*t40*t46*t9-3.69307372
     3         0240627d+3*t36*t37*t46*t9+2.769805290180471d+3*t34*t36*t4
     4         6*t9+3.37614805358933d+0*t117*t39*t40*t9-2.10403636633445
     5         1d+4*t115*t33*t40*t9-6.75229610717866d+0*t117*t37*t39*t9+
     6         1.2624218198006704d+5*t115*t33*t37*t9-1.8936327297010064d
     7         +5*t115*t33*t34*t9+8.416145465337804d+4*t115*t31*t33*t9+6
     8         .093571638397036d+4*t10*t46*t6*t67+1.5275304019588112d+7*
     9         t115*t14*t6*t67-4.628880005935792d+6*t115*t13*t14*t67-1.1
     :         171548003727898d+5*t14*t44*t46*t6+5.077976365330863d+4*t1
     ;         2*t17*t46*t6-1.2729420016323426d+7*t115*t17*t44*t6+3.7137
     <         62858948263d+1*t10*t117*t44*t6+3.4716600044518436d+6*t115
     =         *t12*t20*t6-3.713762858948263d+1*t117*t12*t14*t6-6.017544
     >         00771653d+6*t10*t115*t116*t6+5.077976365330862d+4*t13*t14
     ?         *t44*t46-9.232684300601568d+2*t39*t40*t46+9.2326843006015
     @         68d+2*t37*t39*t46+2.769805290180471d+3*t33*t37*t38*t46-8.
     1         309415870541413d+3*t33*t34*t38*t46+5.539610580360942d+3*t
     2         31*t33*t38*t46+5.539610580360942d+3*t30*t34*t35*t46-1.477
     3         2294880962508d+4*t30*t31*t35*t46+9.232684300601568d+3*t28
     4         *t30*t35*t46+9.232684300601568d+3*t27*t31*t32*t46-2.30817
     5         1075150392d+4*t27*t28*t32*t46+1.3849026450902355d+4*t25*t
     6         27*t32*t46+1.3849026450902355d+4*t24*t28*t29*t46-3.323766
     7         348216566d+4*t24*t25*t29*t46+1.9388637031263292d+4*t22*t2
     8         4*t29*t46+1.9388637031263292d+4*t21*t25*t26*t46-4.5240153
     9         072947686d+4*t21*t22*t26*t46+2.585151604168439d+4*t19*t21
     :         *t26*t46+2.585151604168439d+4*t18*t22*t23*t46-5.908917952
     ;         385003d+4*t18*t19*t23*t46+3.3237663482165647d+4*t16*t18*t
     <         23*t46+3.3237663482165647d+4*t15*t19*t20*t46-7.4784742834
     =         87271d+4*t15*t16*t20*t46+4.154707935270705d+4*t13*t15*t20
     >         *t46+4.1547079352707056d+4*t12*t16*t17*t46-9.232684300601
     ?         568d+4*t12*t13*t17*t46-3.4716600044518436d+6*t115*t16*t17
     @         *t44+1.1572200014839477d+7*t115*t13*t17*t44-3.37614805358
     1         933d+0*t117*t40*t41+2.104036366334451d+4*t115*t36*t40+6.7
     2         5229610717866d+0*t117*t36*t37*t38-8.416145465337804d+4*t1
     3         15*t30*t37*t38-1.0128444160767991d+1*t117*t34*t36*t38+3.7
     4         87265459402013d+5*t115*t30*t34*t38-5.049687279202681d+5*t
     5         115*t30*t31*t38+2.1040363663344508d+5*t115*t28*t30*t38-4.
     6         208072732668902d+4*t115*t36*t37+2.104036366334451d+4*t115
     7         *t34*t36+1.0128444160767991d+1*t117*t33*t34*t35-2.1040363
     8         66334451d+5*t115*t27*t34*t35-1.350459221435732d+1*t117*t3
     9         1*t33*t35+8.416145465337803d+5*t115*t27*t31*t35-1.0520181
     :         831672254d+6*t115*t27*t28*t35+4.208072732668902d+5*t115*t
     ;         25*t27*t35+1.350459221435732d+1*t117*t30*t31*t32-4.208072
     <         7326689016d+5*t115*t24*t31*t32-1.6880740267946648d+1*t117
     =         *t28*t30*t32+1.578027274750838d+6*t115*t24*t28*t32-1.8936
     >         32729701006d+6*t115*t24*t25*t32+7.364127282170578d+5*t115
     ?         *t22*t24*t32+1.6880740267946648d+1*t117*t27*t28*t29-7.364
     @         127282170578d+5*t115*t21*t28*t29-2.0256888321535982d+1*t1
     1         17*t25*t27*t29+2.6510858215814087d+6*t115*t21*t25*t29-3.0
     2         929334585116436d+6*t115*t21*t22*t29+1.1782603651472927d+6
     3         *t115*t19*t21*t29+2.0256888321535982d+1*t117*t24*t25*t26-
     4         1.1782603651472927d+6*t115*t18*t25*t26-2.3633036375125308
     5         d+1*t117*t22*t24*t26+4.123911278015524d+6*t115*t18*t22*t2
     6         6-4.713041460589171d+6*t115*t18*t19*t26+1.767390547720939
     7         d+6*t115*t16*t18*t26+2.3633036375125308d+1*t117*t21*t22*t
     8         23-1.7673905477209387d+6*t115*t15*t22*t23-2.7009184428714
     9         64d+1*t117*t19*t21*t23+6.059624735043219d+6*t115*t15*t19*
     :         t23-6.817077826923622d+6*t115*t15*t16*t23+2.5248436396013
     ;         41d+6*t115*t13*t15*t23+2.700918442871464d+1*t117*t18*t19*
     <         t20-2.524843639601341d+6*t115*t12*t19*t20-3.0385332482303
     =         967d+1*t117*t16*t18*t20+8.521347283654526d+6*t115*t12*t16
     >         *t20-9.468163648505028d+6*t115*t12*t13*t20+3.038533248230
     ?         3967d+1*t117*t15*t16*t17-3.3761480535893296d+1*t117*t13*t
     @         15*t17+3.3761480535893296d+1*t117*t12*t13*t14)*wght-2.791
     1         577209047299d+0*t45*t5*t76*wght-3.7221029453963983d+0*t42
     2         *t46*t76*wght-9.305257363490997d-1*t42*t5*(-1.97024140900
     3         37443d-5*t71*t72/rhoa**11+3.279403591520427d-3*t47*t48/rh
     4         oa**8.333333333333334d+0-1.648227300332581d-1*gammaaa*t3/
     5         rhoa**5.666666666666667d+0+3.7270905095708196d-8*gammaaa*
     6         *4*t118/rhoa**1.3666666666666666d+1)*wght+2.7571132928862
     7         21d-1*t4*t42*wght/t7-2.791577209047299d+0*t5*t51*t69*wght
     8         -3.7221029453963983d+0*t4*t46*t69*wght-7.444205890792797d
     9         +0*t45*t46*t51*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -1.8610514726981994d+0*t45*t5*t78*w
     1         ght-2.481401963597599d+0*t42*t46*t78*wght-9.3052573634909
     2         97d-1*t42*t5*(3.531915643569816d-2*t3*t75-9.3765938467226
     3         28d-4*gammaaa*t48*t74+6.610678411788878d-6*t47*t72*t73-1.
     4         3976589410890572d-8*t118*t71/rhoa**1.2666666666666665d+1)
     5         *wght-4.1356699393293317d-1*t42*t53*t70*wght-9.3052573634
     6         90997d-1*t5*t53*t69*wght-2.481401963597599d+0*t45*t46*t53
     7         *wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -9.305257363490997d-1*t45*t5*t80*w
     1         ght-1.2407009817987995d+0*t42*t46*t80*wght-9.305257363490
     2         997d-1*t42*t5*(-1.944317179937905d-6*gammaaa*t72*t77+1.73
     3         1063479394947d-4*t48*t49+5.241221029083964d-9*t118*t47/rh
     4         oa**1.1666666666666665d+1)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t42*t5*(4.374713654860287d-7*t72*t79-1.96545
     2         7885906487d-9*gammaaa*t118/rhoa**1.0666666666666666d+1)*w
     3         ght
            Mmat3(iq,D3_RA_RA_TA) = -9.305257363490997d-1*t4*t5*(4.05137
     1         76643071964d+1*t36*t40*t68*t9-1.2154132992921589d+2*t34*t
     2         36*t68*t9-2.769805290180471d+3*t33*t40*t5*t9+5.5396105803
     3         60942d+3*t33*t37*t5*t9+8.309415870541416d+3*t33*t34*t5*t9
     4         -1.1079221160721883d+4*t31*t33*t5*t9+2.6739092584427493d+
     5         3*t10*t6*t67*t68-2.228257715368958d+3*t12*t17*t6*t68+2.22
     6         82577153689576d+3*t13*t14*t44*t68-4.0513776643071964d+1*t
     7         37*t39*t68+1.2154132992921589d+2*t33*t37*t38*t68-2.430826
     8         5985843178d+2*t31*t33*t38*t68+2.4308265985843178d+2*t30*t
     9         34*t35*t68-4.051377664307196d+2*t28*t30*t35*t68+4.0513776
     :         643071964d+2*t27*t31*t32*t68-6.077066496460795d+2*t25*t27
     ;         *t32*t68+6.077066496460793d+2*t24*t28*t29*t68-8.507893095
     <         045111d+2*t22*t24*t29*t68+8.507893095045113d+2*t21*t25*t2
     =         6*t68-1.1343857460060149d+3*t19*t21*t26*t68+1.13438574600
     >         60149d+3*t18*t22*t23*t68-1.4584959591505905d+3*t16*t18*t2
     ?         3*t68+1.4584959591505908d+3*t15*t19*t20*t68-1.82311994893
     @         82383d+3*t13*t15*t20*t68+1.8231199489382383d+3*t12*t16*t1
     1         7*t68+6.702928802236738d+5*t14*t5*t6*t67-6.09357163839703
     2         5d+5*t13*t14*t5*t67+5.585774001863949d+5*t17*t44*t5*t6-4.
     3         5701787287977763d+5*t12*t20*t5*t6-7.921643129916147d+5*t1
     4         0*t116*t5*t6-4.5701787287977763d+5*t16*t17*t44*t5+5.07797
     5         6365330863d+5*t13*t17*t44*t5+9.232684300601569d+2*t36*t40
     6         *t5-1.1079221160721883d+4*t30*t37*t38*t5+1.66188317410828
     7         32d+4*t30*t34*t38*t5+2.2158442321443766d+4*t30*t31*t38*t5
     8         -2.7698052901804704d+4*t28*t30*t38*t5+1.846536860120314d+
     9         3*t36*t37*t5-2.769805290180471d+3*t34*t36*t5-2.7698052901
     :         80471d+4*t27*t34*t35*t5+3.693073720240629d+4*t27*t31*t35*
     ;         t5+4.616342150300784d+4*t27*t28*t35*t5-5.539610580360942d
     <         +4*t25*t27*t35*t5-5.539610580360941d+4*t24*t31*t32*t5+6.9
     =         24513225451174d+4*t24*t28*t32*t5+8.309415870541416d+4*t24
     >         *t25*t32*t5-9.694318515631648d+4*t22*t24*t32*t5-9.6943185
     ?         15631648d+4*t21*t28*t29*t5+1.1633182218757979d+5*t21*t25*
     @         t29*t5+1.3572045921884307d+5*t21*t22*t29*t5-1.55109096250
     1         10635d+5*t19*t21*t29*t5-1.5510909625010638d+5*t18*t25*t26
     2         *t5+1.809606122917908d+5*t18*t22*t26*t5+2.068121283334751
     3         5d+5*t18*t19*t26*t5-2.3266364437515952d+5*t16*t18*t26*t5-
     4         2.3266364437515955d+5*t15*t22*t23*t5+2.659013078573252d+5
     5         *t15*t19*t23*t5+2.9913897133949086d+5*t15*t16*t23*t5-3.32
     6         37663482165647d+5*t13*t15*t23*t5-3.323766348216565d+5*t12
     7         *t19*t20*t5+3.7392371417436354d+5*t12*t16*t20*t5+4.154707
     8         9352707056d+5*t12*t13*t20*t5)*wght-9.305257363490997d-1*t
     9         5*t76*t82*wght-4.1356699393293317d-1*t4*t70*t82*wght-2.48
     :         1401963597599d+0*t46*t51*t82*wght-1.8610514726981994d+0*t
     ;         5*t51*t81*wght-2.481401963597599d+0*t4*t46*t81*wght+Mmat3
     <         (iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_GAA_TA) = -9.305257363490997d-1*t5*t78*t82*wg
     1         ht-1.2407009817987995d+0*t46*t53*t82*wght-9.3052573634909
     2         97d-1*t5*t53*t81*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)-9.30525736
     1         3490997d-1*t5*t66*t80*wght
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_RA_TA_TA) = -9.305257363490997d-1*t4*t5*(-3.6462
     1         39897876477d+2*t33*t40*t43*t9-7.292479795752954d+2*t33*t3
     2         7*t43*t9+1.0938719693629433d+3*t33*t34*t43*t9+1.458495959
     3         1505908d+3*t31*t33*t43*t9-8.823900552861074d+4*t14*t43*t6
     4         *t67-8.021727775328248d+4*t13*t14*t43*t67+7.3532504607175
     5         61d+4*t17*t43*t44*t6+6.016295831496187d+4*t12*t20*t43*t6-
     6         1.0428246107926724d+5*t10*t116*t43*t6-6.016295831496187d+
     7         4*t16*t17*t43*t44-6.684773146106873d+4*t13*t17*t43*t44-1.
     8         2154132992921589d+2*t36*t40*t43-1.4584959591505908d+3*t30
     9         *t37*t38*t43-2.1877439387258865d+3*t30*t34*t38*t43+2.9169
     :         919183011817d+3*t30*t31*t38*t43+3.646239897876477d+3*t28*
     ;         t30*t38*t43+2.4308265985843178d+2*t36*t37*t43+3.646239897
     <         876477d+2*t34*t36*t43-3.646239897876477d+3*t27*t34*t35*t4
     =         3-4.8616531971686355d+3*t27*t31*t35*t43+6.077066496460795
     >         d+3*t27*t28*t35*t43+7.292479795752954d+3*t25*t27*t35*t43-
     ?         7.292479795752954d+3*t24*t31*t32*t43-9.115599744691192d+3
     @         *t24*t28*t32*t43+1.0938719693629431d+4*t24*t25*t32*t43+1.
     1         2761839642567668d+4*t22*t24*t32*t43-1.2761839642567668d+4
     2         *t21*t28*t29*t43-1.5314207571081204d+4*t21*t25*t29*t43+1.
     3         786657549959474d+4*t21*t22*t29*t43+2.0418943428108272d+4*
     4         t19*t21*t29*t43-2.0418943428108272d+4*t18*t25*t26*t43-2.3
     5         822100666126317d+4*t18*t22*t26*t43+2.722525790414436d+4*t
     6         18*t19*t26*t43+3.062841514216241d+4*t16*t18*t26*t43-3.062
     7         841514216241d+4*t15*t22*t23*t43-3.5003903019614174d+4*t15
     8         *t19*t23*t43+3.9379390897065947d+4*t15*t16*t23*t43+4.3754
     9         87877451772d+4*t13*t15*t23*t43-4.375487877451772d+4*t12*t
     :         19*t20*t43-4.9224238621332433d+4*t12*t16*t20*t43+5.469359
     ;         846814714d+4*t12*t13*t20*t43)*wght-9.305257363490997d-1*t
     <         119*t5*t51*wght-1.2407009817987995d+0*t119*t4*t46*wght+Mm
     =         at3(iq,D3_RA_TA_TA)
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)-9.3052573634
     1         90997d-1*t114*t5*t53*wght
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)-9.305257363490
     1         997d-1*t4*t5*(t21*(-1.4d+1*t99-1.4d+1*(t98+t97)-2.352d+3*
     2         t22*t29-2.688d+3*t19*t29)+t18*(-1.6d+1*t96-1.6d+1*(t95+t9
     3         4)-3.584d+3*t19*t26-4.032d+3*t16*t26)+t15*(-1.8d+1*t93-1.
     4         8d+1*(t92+t91)-5.184d+3*t16*t23-5.76d+3*t13*t23)+t12*(-2.
     5         0d+1*t90-2.0d+1*(t89+t88)-7.92d+3*t20*t6-7.2d+3*t13*t20)+
     6         t33*(-1.44d+2*t34*t9-1.92d+2*t31*t9-6.0d+0*t111-6.0d+0*(t
     7         110+t109))+t44*(-2.2d+1*t87-2.2d+1*(t86+t85)-9.68d+3*t17*
     8         t6)+t67*(-2.4d+1*t84-1.1616d+4*t14*t6)-1.3728d+4*t10*t116
     9         *t6+t36*(-4.0d+0*(4.0d+0*t40+t112)-3.2d+1*t37-4.8d+1*t34-
     :         4.0d+0*t113)+t30*(-3.84d+2*t31*t38-4.8d+2*t28*t38-8.0d+0*
     ;         t108-8.0d+0*(t107+t106))+t27*(-8.0d+2*t28*t35-9.6d+2*t25*
     <         t35-1.0d+1*t105-1.0d+1*(t104+t103))+t24*(-1.44d+3*t25*t32
     =         -1.68d+3*t22*t32-1.2d+1*t102-1.2d+1*(t101+t100)))*wght
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = param(12)
            t7 = rhob**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taub
            t10 = t9**11
            t11 = 2.0d+0*taub+t8
            t12 = 1/t11**11
            t13 = param(11)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(10)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(9)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(8)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(7)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(6)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(5)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(4)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(3)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(2)
            t41 = 1/t11
            t42 = t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+t31*t32*
     1         t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t21+t16*t
     2         17*t18+t13*t14*t15+param(1)
            t43 = rhob**6.666666666666666d-1
            t44 = 1/t11**12
            t45 = -1.5192666241151986d+1*t39*t40*t43*t9+3.03853324823039
     1         7d+1*t37*t39*t43*t9-1.6711932865267184d+2*t10*t43*t44*t6+
     2         1.6711932865267184d+2*t12*t14*t43*t6+1.5192666241151986d+
     3         1*t40*t41*t43-3.038533248230397d+1*t36*t37*t38*t43+4.5577
     4         998723455965d+1*t34*t36*t38*t43-4.5577998723455965d+1*t33
     5         *t34*t35*t43+6.077066496460795d+1*t31*t33*t35*t43-6.07706
     6         6496460795d+1*t30*t31*t32*t43+7.596333120575993d+1*t28*t3
     7         0*t32*t43-7.596333120575993d+1*t27*t28*t29*t43+9.11559974
     8         4691193d+1*t25*t27*t29*t43-9.115599744691193d+1*t24*t25*t
     9         26*t43+1.0634866368806391d+2*t22*t24*t26*t43-1.0634866368
     :         806391d+2*t21*t22*t23*t43+1.2154132992921589d+2*t19*t21*t
     ;         23*t43-1.2154132992921589d+2*t18*t19*t20*t43+1.3673399617
     <         036788d+2*t16*t18*t20*t43-1.3673399617036788d+2*t15*t16*t
     =         17*t43+1.5192666241151986d+2*t13*t15*t17*t43-1.5192666241
     >         151986d+2*t12*t13*t14*t43
            t46 = rhob**3.333333333333333d-1
            t47 = gammabb**2
            t48 = 1/t2**2
            t49 = 1/rhob**6.333333333333333d+0
            t50 = 1/rhob**3.6666666666666664d+0
            t51 = 4.327658698487368d-5*t47*t48*t49-9.632497209735863d-3*
     1         gammabb*t3*t50
            t52 = 1/rhob**5.333333333333333d+0
            t53 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         abb*t48*t52
            t54 = -2.2d+1*t10*t44*t6
            t55 = -2.2d+1*t14*t6-2.0d+1*t13*t14
            t56 = -1.8d+1*t16*t17-2.0d+1*t13*t17
            t57 = -1.6d+1*t19*t20-1.8d+1*t16*t20
            t58 = -1.4d+1*t22*t23-1.6d+1*t19*t23
            t59 = -1.2d+1*t25*t26-1.4d+1*t22*t26
            t60 = -1.0d+1*t28*t29-1.2d+1*t25*t29
            t61 = -8.0d+0*t31*t32-1.0d+1*t28*t32
            t62 = -6.0d+0*t34*t35-8.0d+0*t31*t35
            t63 = -4.0d+0*t37*t38-6.0d+0*t34*t38
            t64 = -2.0d+0*t40*t9-4.0d+0*t37*t9
            t65 = -2.0d+0*t40*t41
            t66 = t65+t39*t64+t36*t63+t33*t62+t30*t61+t27*t60+t24*t59+t2
     1         1*t58+t18*t57+t15*t56+t12*t55+t54
            t67 = 1/t11**13
            t68 = 1/t46
            t69 = -1.0128444160767991d+1*t39*t40*t68*t9+2.02568883215359
     1         82d+1*t37*t39*t68*t9+4.6163421503007845d+2*t36*t40*t5*t9-
     2         1.846536860120314d+3*t36*t37*t5*t9+1.3849026450902355d+3*
     3         t34*t36*t5*t9-1.1141288576844789d+2*t10*t44*t6*t68+1.1141
     4         288576844789d+2*t12*t14*t6*t68+1.0128444160767991d+1*t40*
     5         t41*t68-2.0256888321535982d+1*t36*t37*t38*t68+3.038533248
     6         230397d+1*t34*t36*t38*t68-3.038533248230397d+1*t33*t34*t3
     7         5*t68+4.0513776643071964d+1*t31*t33*t35*t68-4.05137766430
     8         71964d+1*t30*t31*t32*t68+5.064222080383994d+1*t28*t30*t32
     9         *t68-5.064222080383994d+1*t27*t28*t29*t68+6.0770664964607
     :         95d+1*t25*t27*t29*t68-6.077066496460795d+1*t24*t25*t26*t6
     ;         8+7.089910912537593d+1*t22*t24*t26*t68-7.089910912537593d
     <         +1*t21*t22*t23*t68+8.102755328614393d+1*t19*t21*t23*t68-8
     =         .102755328614393d+1*t18*t19*t20*t68+9.115599744691191d+1*
     >         t16*t18*t20*t68-9.115599744691191d+1*t15*t16*t17*t68+1.01
     ?         2844416076799d+2*t13*t15*t17*t68-1.012844416076799d+2*t12
     @         *t13*t14*t68+3.046785819198518d+4*t10*t5*t6*t67-5.5857740
     1         01863949d+4*t14*t44*t5*t6+2.5389881826654315d+4*t12*t17*t
     2         5*t6+2.5389881826654315d+4*t13*t14*t44*t5-4.6163421503007
     3         845d+2*t39*t40*t5+4.6163421503007845d+2*t37*t39*t5+1.3849
     4         026450902355d+3*t33*t37*t38*t5-4.1547079352707067d+3*t33*
     5         t34*t38*t5+2.769805290180471d+3*t31*t33*t38*t5+2.76980529
     6         0180471d+3*t30*t34*t35*t5-7.386147440481255d+3*t30*t31*t3
     7         5*t5+4.616342150300784d+3*t28*t30*t35*t5+4.61634215030078
     8         4d+3*t27*t31*t32*t5-1.154085537575196d+4*t27*t28*t32*t5+6
     9         .924513225451178d+3*t25*t27*t32*t5+6.924513225451178d+3*t
     :         24*t28*t29*t5-1.661883174108283d+4*t24*t25*t29*t5+9.69431
     ;         8515631648d+3*t22*t24*t29*t5+9.694318515631648d+3*t21*t25
     <         *t26*t5-2.2620076536473846d+4*t21*t22*t26*t5+1.2925758020
     =         842196d+4*t19*t21*t26*t5+1.2925758020842196d+4*t18*t22*t2
     >         3*t5-2.954458976192502d+4*t18*t19*t23*t5+1.66188317410828
     ?         23d+4*t16*t18*t23*t5+1.6618831741082823d+4*t15*t19*t20*t5
     @         -3.7392371417436354d+4*t15*t16*t20*t5+2.0773539676353528d
     1         +4*t13*t15*t20*t5+2.0773539676353528d+4*t12*t16*t17*t5-4.
     2         616342150300784d+4*t12*t13*t17*t5
            t70 = 1/t43
            t71 = gammabb**3
            t72 = 1/t2**3
            t73 = 1/rhob**10
            t74 = 1/rhob**7.333333333333333d+0
            t75 = 1/rhob**4.666666666666667d+0
            t76 = 3.531915643569816d-2*gammabb*t3*t75-3.8948928286386303
     1         d-4*t47*t48*t74+1.0369691626335495d-6*t71*t72*t73
            t77 = 1/rhob**9
            t78 = -3.8886343598758105d-7*t47*t72*t77-9.632497209735863d-
     1         3*t3*t50+1.29829760954621d-4*gammabb*t48*t49
            t79 = 1/rhob**8
            t80 = 1.4582378849534292d-7*gammabb*t72*t79-3.24574402386552
     1         55d-5*t48*t52
            t81 = 6.077066496460795d+1*t36*t40*t43*t9-1.8231199489382385
     1         d+2*t34*t36*t43*t9+4.0108638876641245d+3*t10*t43*t6*t67-3
     2         .342386573053437d+3*t12*t17*t43*t6+3.342386573053437d+3*t
     3         13*t14*t43*t44-6.077066496460795d+1*t37*t39*t43+1.8231199
     4         489382385d+2*t33*t37*t38*t43-3.646239897876477d+2*t31*t33
     5         *t38*t43+3.646239897876477d+2*t30*t34*t35*t43-6.077066496
     6         460795d+2*t28*t30*t35*t43+6.077066496460795d+2*t27*t31*t3
     7         2*t43-9.115599744691192d+2*t25*t27*t32*t43+9.115599744691
     8         192d+2*t24*t28*t29*t43-1.276183964256767d+3*t22*t24*t29*t
     9         43+1.276183964256767d+3*t21*t25*t26*t43-1.701578619009022
     :         6d+3*t19*t21*t26*t43+1.7015786190090226d+3*t18*t22*t23*t4
     ;         3-2.187743938725886d+3*t16*t18*t23*t43+2.187743938725886d
     <         +3*t15*t19*t20*t43-2.734679923407357d+3*t13*t15*t20*t43+2
     =         .734679923407357d+3*t12*t16*t17*t43
            t82 = -2.0d+0*t39*t40*t9-4.0d+0*t37*t39*t9+t65-2.2d+1*t12*t1
     1         4*t6+t54-4.0d+0*t36*t37*t38-6.0d+0*t34*t36*t38-6.0d+0*t33
     2         *t34*t35-8.0d+0*t31*t33*t35-8.0d+0*t30*t31*t32-1.0d+1*t28
     3         *t30*t32-1.0d+1*t27*t28*t29-1.2d+1*t25*t27*t29-1.2d+1*t24
     4         *t25*t26-1.4d+1*t22*t24*t26-1.4d+1*t21*t22*t23-1.6d+1*t19
     5         *t21*t23-1.6d+1*t18*t19*t20-1.8d+1*t16*t18*t20-1.8d+1*t15
     6         *t16*t17-2.0d+1*t13*t15*t17-2.0d+1*t12*t13*t14
            t83 = 5.28d+2*t10*t6*t67
            t84 = 4.84d+2*t14*t6-2.2d+1*t55
            t85 = 4.0d+2*t13*t17
            t86 = 4.4d+2*t17*t6
            t87 = t86+t85-2.0d+1*t56
            t88 = 3.24d+2*t16*t20
            t89 = 3.6d+2*t13*t20
            t90 = t89+t88-1.8d+1*t57
            t91 = 2.56d+2*t19*t23
            t92 = 2.88d+2*t16*t23
            t93 = t92+t91-1.6d+1*t58
            t94 = 1.96d+2*t22*t26
            t95 = 2.24d+2*t19*t26
            t96 = t95+t94-1.4d+1*t59
            t97 = 1.44d+2*t25*t29
            t98 = 1.68d+2*t22*t29
            t99 = t98+t97-1.2d+1*t60
            t100 = 1.0d+2*t28*t32
            t101 = 1.2d+2*t25*t32
            t102 = -1.0d+1*t61+t101+t100
            t103 = 6.4d+1*t31*t35
            t104 = 8.0d+1*t28*t35
            t105 = -8.0d+0*t62+t104+t103
            t106 = 3.6d+1*t34*t38
            t107 = 4.8d+1*t31*t38
            t108 = -6.0d+0*t63+t107+t106
            t109 = 1.6d+1*t37*t9
            t110 = 2.4d+1*t34*t9
            t111 = -4.0d+0*t64+t110+t109
            t112 = 8.0d+0*t37
            t113 = 8.0d+0*t40+t112
            t114 = t24*t99+t21*t96+t18*t93+t15*t90+t12*t87+t44*t84+t83+t
     1         113*t39+t111*t36+t108*t33+t105*t30+t102*t27
            t115 = rhob**2
            t116 = 1/t11**14
            t117 = 1/t5
            t118 = 1/t2**4
            t119 = 8.0d+0*t36*t40*t9+3.2d+1*t36*t37*t9+2.4d+1*t34*t36*t9
     1         +t83+9.68d+2*t14*t44*t6+4.4d+2*t12*t17*t6+4.4d+2*t13*t14*
     2         t44+8.0d+0*t39*t40+8.0d+0*t37*t39+2.4d+1*t33*t37*t38+7.2d
     3         +1*t33*t34*t38+4.8d+1*t31*t33*t38+4.8d+1*t30*t34*t35+1.28
     4         d+2*t30*t31*t35+8.0d+1*t28*t30*t35+8.0d+1*t27*t31*t32+2.0
     5         d+2*t27*t28*t32+1.2d+2*t25*t27*t32+1.2d+2*t24*t28*t29+2.8
     6         8d+2*t24*t25*t29+1.68d+2*t22*t24*t29+1.68d+2*t21*t25*t26+
     7         3.92d+2*t21*t22*t26+2.24d+2*t19*t21*t26+2.24d+2*t18*t22*t
     8         23+5.12d+2*t18*t19*t23+2.88d+2*t16*t18*t23+2.88d+2*t15*t1
     9         9*t20+6.48d+2*t15*t16*t20+3.6d+2*t13*t15*t20+3.6d+2*t12*t
     :         16*t17+8.0d+2*t12*t13*t17
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t42*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t42*t5*t51*wght-9.305
     1         257363490997d-1*t4*t45*t5*wght-1.2407009817987995d+0*t4*t
     2         42*t46*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t42*t
     1         5*t53*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t4*t5*t
     1         66*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = -9.305257363490997d-1*t42*t5*t76*wght-4
     1         .1356699393293317d-1*t4*t42*t70*wght-9.305257363490997d-1
     2         *t4*t5*t69*wght-1.8610514726981994d+0*t45*t5*t51*wght-2.4
     3         81401963597599d+0*t42*t46*t51*wght-2.481401963597599d+0*t
     4         4*t45*t46*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490997d-1*t42*t5*t78*wght-
     1         9.305257363490997d-1*t45*t5*t53*wght-1.2407009817987995d+
     2         0*t42*t46*t53*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t42*t5*t80*wght
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = -9.305257363490997d-1*t5*t51*t82*wght-1
     1         .2407009817987995d+0*t4*t46*t82*wght-9.305257363490997d-1
     2         *t4*t5*t81*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t5*t53*t66*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t114*t4*t5*wght
            Amat3(iq,D3_RB_RB_RB) = t70*(-1.2407009817987995d+0*t42*t51*
     1         wght-1.2407009817987995d+0*t4*t45*wght)-9.305257363490997
     2         d-1*t4*t5*(9.232684300601568d+2*t36*t40*t46*t9-3.69307372
     3         0240627d+3*t36*t37*t46*t9+2.769805290180471d+3*t34*t36*t4
     4         6*t9+3.37614805358933d+0*t117*t39*t40*t9-2.10403636633445
     5         1d+4*t115*t33*t40*t9-6.75229610717866d+0*t117*t37*t39*t9+
     6         1.2624218198006704d+5*t115*t33*t37*t9-1.8936327297010064d
     7         +5*t115*t33*t34*t9+8.416145465337804d+4*t115*t31*t33*t9+6
     8         .093571638397036d+4*t10*t46*t6*t67+1.5275304019588112d+7*
     9         t115*t14*t6*t67-4.628880005935792d+6*t115*t13*t14*t67-1.1
     :         171548003727898d+5*t14*t44*t46*t6+5.077976365330863d+4*t1
     ;         2*t17*t46*t6-1.2729420016323426d+7*t115*t17*t44*t6+3.7137
     <         62858948263d+1*t10*t117*t44*t6+3.4716600044518436d+6*t115
     =         *t12*t20*t6-3.713762858948263d+1*t117*t12*t14*t6-6.017544
     >         00771653d+6*t10*t115*t116*t6+5.077976365330862d+4*t13*t14
     ?         *t44*t46-9.232684300601568d+2*t39*t40*t46+9.2326843006015
     @         68d+2*t37*t39*t46+2.769805290180471d+3*t33*t37*t38*t46-8.
     1         309415870541413d+3*t33*t34*t38*t46+5.539610580360942d+3*t
     2         31*t33*t38*t46+5.539610580360942d+3*t30*t34*t35*t46-1.477
     3         2294880962508d+4*t30*t31*t35*t46+9.232684300601568d+3*t28
     4         *t30*t35*t46+9.232684300601568d+3*t27*t31*t32*t46-2.30817
     5         1075150392d+4*t27*t28*t32*t46+1.3849026450902355d+4*t25*t
     6         27*t32*t46+1.3849026450902355d+4*t24*t28*t29*t46-3.323766
     7         348216566d+4*t24*t25*t29*t46+1.9388637031263292d+4*t22*t2
     8         4*t29*t46+1.9388637031263292d+4*t21*t25*t26*t46-4.5240153
     9         072947686d+4*t21*t22*t26*t46+2.585151604168439d+4*t19*t21
     :         *t26*t46+2.585151604168439d+4*t18*t22*t23*t46-5.908917952
     ;         385003d+4*t18*t19*t23*t46+3.3237663482165647d+4*t16*t18*t
     <         23*t46+3.3237663482165647d+4*t15*t19*t20*t46-7.4784742834
     =         87271d+4*t15*t16*t20*t46+4.154707935270705d+4*t13*t15*t20
     >         *t46+4.1547079352707056d+4*t12*t16*t17*t46-9.232684300601
     ?         568d+4*t12*t13*t17*t46-3.4716600044518436d+6*t115*t16*t17
     @         *t44+1.1572200014839477d+7*t115*t13*t17*t44-3.37614805358
     1         933d+0*t117*t40*t41+2.104036366334451d+4*t115*t36*t40+6.7
     2         5229610717866d+0*t117*t36*t37*t38-8.416145465337804d+4*t1
     3         15*t30*t37*t38-1.0128444160767991d+1*t117*t34*t36*t38+3.7
     4         87265459402013d+5*t115*t30*t34*t38-5.049687279202681d+5*t
     5         115*t30*t31*t38+2.1040363663344508d+5*t115*t28*t30*t38-4.
     6         208072732668902d+4*t115*t36*t37+2.104036366334451d+4*t115
     7         *t34*t36+1.0128444160767991d+1*t117*t33*t34*t35-2.1040363
     8         66334451d+5*t115*t27*t34*t35-1.350459221435732d+1*t117*t3
     9         1*t33*t35+8.416145465337803d+5*t115*t27*t31*t35-1.0520181
     :         831672254d+6*t115*t27*t28*t35+4.208072732668902d+5*t115*t
     ;         25*t27*t35+1.350459221435732d+1*t117*t30*t31*t32-4.208072
     <         7326689016d+5*t115*t24*t31*t32-1.6880740267946648d+1*t117
     =         *t28*t30*t32+1.578027274750838d+6*t115*t24*t28*t32-1.8936
     >         32729701006d+6*t115*t24*t25*t32+7.364127282170578d+5*t115
     ?         *t22*t24*t32+1.6880740267946648d+1*t117*t27*t28*t29-7.364
     @         127282170578d+5*t115*t21*t28*t29-2.0256888321535982d+1*t1
     1         17*t25*t27*t29+2.6510858215814087d+6*t115*t21*t25*t29-3.0
     2         929334585116436d+6*t115*t21*t22*t29+1.1782603651472927d+6
     3         *t115*t19*t21*t29+2.0256888321535982d+1*t117*t24*t25*t26-
     4         1.1782603651472927d+6*t115*t18*t25*t26-2.3633036375125308
     5         d+1*t117*t22*t24*t26+4.123911278015524d+6*t115*t18*t22*t2
     6         6-4.713041460589171d+6*t115*t18*t19*t26+1.767390547720939
     7         d+6*t115*t16*t18*t26+2.3633036375125308d+1*t117*t21*t22*t
     8         23-1.7673905477209387d+6*t115*t15*t22*t23-2.7009184428714
     9         64d+1*t117*t19*t21*t23+6.059624735043219d+6*t115*t15*t19*
     :         t23-6.817077826923622d+6*t115*t15*t16*t23+2.5248436396013
     ;         41d+6*t115*t13*t15*t23+2.700918442871464d+1*t117*t18*t19*
     <         t20-2.524843639601341d+6*t115*t12*t19*t20-3.0385332482303
     =         967d+1*t117*t16*t18*t20+8.521347283654526d+6*t115*t12*t16
     >         *t20-9.468163648505028d+6*t115*t12*t13*t20+3.038533248230
     ?         3967d+1*t117*t15*t16*t17-3.3761480535893296d+1*t117*t13*t
     @         15*t17+3.3761480535893296d+1*t117*t12*t13*t14)*wght-2.791
     1         577209047299d+0*t45*t5*t76*wght-3.7221029453963983d+0*t42
     2         *t46*t76*wght-9.305257363490997d-1*t42*t5*(-1.97024140900
     3         37443d-5*t71*t72/rhob**11+3.279403591520427d-3*t47*t48/rh
     4         ob**8.333333333333334d+0-1.648227300332581d-1*gammabb*t3/
     5         rhob**5.666666666666667d+0+3.7270905095708196d-8*gammabb*
     6         *4*t118/rhob**1.3666666666666666d+1)*wght+2.7571132928862
     7         21d-1*t4*t42*wght/t7-2.791577209047299d+0*t5*t51*t69*wght
     8         -3.7221029453963983d+0*t4*t46*t69*wght-7.444205890792797d
     9         +0*t45*t46*t51*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -1.8610514726981994d+0*t45*t5*t78*w
     1         ght-2.481401963597599d+0*t42*t46*t78*wght-9.3052573634909
     2         97d-1*t42*t5*(3.531915643569816d-2*t3*t75-9.3765938467226
     3         28d-4*gammabb*t48*t74+6.610678411788878d-6*t47*t72*t73-1.
     4         3976589410890572d-8*t118*t71/rhob**1.2666666666666665d+1)
     5         *wght-4.1356699393293317d-1*t42*t53*t70*wght-9.3052573634
     6         90997d-1*t5*t53*t69*wght-2.481401963597599d+0*t45*t46*t53
     7         *wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -9.305257363490997d-1*t45*t5*t80*w
     1         ght-1.2407009817987995d+0*t42*t46*t80*wght-9.305257363490
     2         997d-1*t42*t5*(-1.944317179937905d-6*gammabb*t72*t77+1.73
     3         1063479394947d-4*t48*t49+5.241221029083964d-9*t118*t47/rh
     4         ob**1.1666666666666665d+1)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490997d-1*t42*t5*(4.374713654860287d-7*t72*t79-1.96545
     2         7885906487d-9*gammabb*t118/rhob**1.0666666666666666d+1)*w
     3         ght
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = -9.305257363490997d-1*t4*t5*(4.05137
     1         76643071964d+1*t36*t40*t68*t9-1.2154132992921589d+2*t34*t
     2         36*t68*t9-2.769805290180471d+3*t33*t40*t5*t9+5.5396105803
     3         60942d+3*t33*t37*t5*t9+8.309415870541416d+3*t33*t34*t5*t9
     4         -1.1079221160721883d+4*t31*t33*t5*t9+2.6739092584427493d+
     5         3*t10*t6*t67*t68-2.228257715368958d+3*t12*t17*t6*t68+2.22
     6         82577153689576d+3*t13*t14*t44*t68-4.0513776643071964d+1*t
     7         37*t39*t68+1.2154132992921589d+2*t33*t37*t38*t68-2.430826
     8         5985843178d+2*t31*t33*t38*t68+2.4308265985843178d+2*t30*t
     9         34*t35*t68-4.051377664307196d+2*t28*t30*t35*t68+4.0513776
     :         643071964d+2*t27*t31*t32*t68-6.077066496460795d+2*t25*t27
     ;         *t32*t68+6.077066496460793d+2*t24*t28*t29*t68-8.507893095
     <         045111d+2*t22*t24*t29*t68+8.507893095045113d+2*t21*t25*t2
     =         6*t68-1.1343857460060149d+3*t19*t21*t26*t68+1.13438574600
     >         60149d+3*t18*t22*t23*t68-1.4584959591505905d+3*t16*t18*t2
     ?         3*t68+1.4584959591505908d+3*t15*t19*t20*t68-1.82311994893
     @         82383d+3*t13*t15*t20*t68+1.8231199489382383d+3*t12*t16*t1
     1         7*t68+6.702928802236738d+5*t14*t5*t6*t67-6.09357163839703
     2         5d+5*t13*t14*t5*t67+5.585774001863949d+5*t17*t44*t5*t6-4.
     3         5701787287977763d+5*t12*t20*t5*t6-7.921643129916147d+5*t1
     4         0*t116*t5*t6-4.5701787287977763d+5*t16*t17*t44*t5+5.07797
     5         6365330863d+5*t13*t17*t44*t5+9.232684300601569d+2*t36*t40
     6         *t5-1.1079221160721883d+4*t30*t37*t38*t5+1.66188317410828
     7         32d+4*t30*t34*t38*t5+2.2158442321443766d+4*t30*t31*t38*t5
     8         -2.7698052901804704d+4*t28*t30*t38*t5+1.846536860120314d+
     9         3*t36*t37*t5-2.769805290180471d+3*t34*t36*t5-2.7698052901
     :         80471d+4*t27*t34*t35*t5+3.693073720240629d+4*t27*t31*t35*
     ;         t5+4.616342150300784d+4*t27*t28*t35*t5-5.539610580360942d
     <         +4*t25*t27*t35*t5-5.539610580360941d+4*t24*t31*t32*t5+6.9
     =         24513225451174d+4*t24*t28*t32*t5+8.309415870541416d+4*t24
     >         *t25*t32*t5-9.694318515631648d+4*t22*t24*t32*t5-9.6943185
     ?         15631648d+4*t21*t28*t29*t5+1.1633182218757979d+5*t21*t25*
     @         t29*t5+1.3572045921884307d+5*t21*t22*t29*t5-1.55109096250
     1         10635d+5*t19*t21*t29*t5-1.5510909625010638d+5*t18*t25*t26
     2         *t5+1.809606122917908d+5*t18*t22*t26*t5+2.068121283334751
     3         5d+5*t18*t19*t26*t5-2.3266364437515952d+5*t16*t18*t26*t5-
     4         2.3266364437515955d+5*t15*t22*t23*t5+2.659013078573252d+5
     5         *t15*t19*t23*t5+2.9913897133949086d+5*t15*t16*t23*t5-3.32
     6         37663482165647d+5*t13*t15*t23*t5-3.323766348216565d+5*t12
     7         *t19*t20*t5+3.7392371417436354d+5*t12*t16*t20*t5+4.154707
     8         9352707056d+5*t12*t13*t20*t5)*wght-9.305257363490997d-1*t
     9         5*t76*t82*wght-4.1356699393293317d-1*t4*t70*t82*wght-2.48
     :         1401963597599d+0*t46*t51*t82*wght-1.8610514726981994d+0*t
     ;         5*t51*t81*wght-2.481401963597599d+0*t4*t46*t81*wght+Mmat3
     <         (iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = -9.305257363490997d-1*t5*t78*t82*wg
     1         ht-1.2407009817987995d+0*t46*t53*t82*wght-9.3052573634909
     2         97d-1*t5*t53*t81*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)-9.30525736
     1         3490997d-1*t5*t66*t80*wght
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = -9.305257363490997d-1*t4*t5*(-3.6462
     1         39897876477d+2*t33*t40*t43*t9-7.292479795752954d+2*t33*t3
     2         7*t43*t9+1.0938719693629433d+3*t33*t34*t43*t9+1.458495959
     3         1505908d+3*t31*t33*t43*t9-8.823900552861074d+4*t14*t43*t6
     4         *t67-8.021727775328248d+4*t13*t14*t43*t67+7.3532504607175
     5         61d+4*t17*t43*t44*t6+6.016295831496187d+4*t12*t20*t43*t6-
     6         1.0428246107926724d+5*t10*t116*t43*t6-6.016295831496187d+
     7         4*t16*t17*t43*t44-6.684773146106873d+4*t13*t17*t43*t44-1.
     8         2154132992921589d+2*t36*t40*t43-1.4584959591505908d+3*t30
     9         *t37*t38*t43-2.1877439387258865d+3*t30*t34*t38*t43+2.9169
     :         919183011817d+3*t30*t31*t38*t43+3.646239897876477d+3*t28*
     ;         t30*t38*t43+2.4308265985843178d+2*t36*t37*t43+3.646239897
     <         876477d+2*t34*t36*t43-3.646239897876477d+3*t27*t34*t35*t4
     =         3-4.8616531971686355d+3*t27*t31*t35*t43+6.077066496460795
     >         d+3*t27*t28*t35*t43+7.292479795752954d+3*t25*t27*t35*t43-
     ?         7.292479795752954d+3*t24*t31*t32*t43-9.115599744691192d+3
     @         *t24*t28*t32*t43+1.0938719693629431d+4*t24*t25*t32*t43+1.
     1         2761839642567668d+4*t22*t24*t32*t43-1.2761839642567668d+4
     2         *t21*t28*t29*t43-1.5314207571081204d+4*t21*t25*t29*t43+1.
     3         786657549959474d+4*t21*t22*t29*t43+2.0418943428108272d+4*
     4         t19*t21*t29*t43-2.0418943428108272d+4*t18*t25*t26*t43-2.3
     5         822100666126317d+4*t18*t22*t26*t43+2.722525790414436d+4*t
     6         18*t19*t26*t43+3.062841514216241d+4*t16*t18*t26*t43-3.062
     7         841514216241d+4*t15*t22*t23*t43-3.5003903019614174d+4*t15
     8         *t19*t23*t43+3.9379390897065947d+4*t15*t16*t23*t43+4.3754
     9         87877451772d+4*t13*t15*t23*t43-4.375487877451772d+4*t12*t
     :         19*t20*t43-4.9224238621332433d+4*t12*t16*t20*t43+5.469359
     ;         846814714d+4*t12*t13*t20*t43)*wght-9.305257363490997d-1*t
     <         119*t5*t51*wght-1.2407009817987995d+0*t119*t4*t46*wght+Mm
     =         at3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
            Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)-9.3052573634
     1         90997d-1*t114*t5*t53*wght
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)-9.305257363490
     1         997d-1*t4*t5*(t21*(-1.4d+1*t99-1.4d+1*(t98+t97)-2.352d+3*
     2         t22*t29-2.688d+3*t19*t29)+t18*(-1.6d+1*t96-1.6d+1*(t95+t9
     3         4)-3.584d+3*t19*t26-4.032d+3*t16*t26)+t15*(-1.8d+1*t93-1.
     4         8d+1*(t92+t91)-5.184d+3*t16*t23-5.76d+3*t13*t23)+t12*(-2.
     5         0d+1*t90-2.0d+1*(t89+t88)-7.92d+3*t20*t6-7.2d+3*t13*t20)+
     6         t33*(-1.44d+2*t34*t9-1.92d+2*t31*t9-6.0d+0*t111-6.0d+0*(t
     7         110+t109))+t44*(-2.2d+1*t87-2.2d+1*(t86+t85)-9.68d+3*t17*
     8         t6)+t67*(-2.4d+1*t84-1.1616d+4*t14*t6)-1.3728d+4*t10*t116
     9         *t6+t36*(-4.0d+0*(4.0d+0*t40+t112)-3.2d+1*t37-4.8d+1*t34-
     :         4.0d+0*t113)+t30*(-3.84d+2*t31*t38-4.8d+2*t28*t38-8.0d+0*
     ;         t108-8.0d+0*(t107+t106))+t27*(-8.0d+2*t28*t35-9.6d+2*t25*
     <         t35-1.0d+1*t105-1.0d+1*(t104+t103))+t24*(-1.44d+3*t25*t32
     =         -1.68d+3*t22*t32-1.2d+1*t102-1.2d+1*(t101+t100)))*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
