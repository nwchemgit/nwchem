C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_m05.F 
C> The nwxcm_x_m05 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_m05 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_3} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_4} &=& {\it t_3}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_5} &=& 2.0\,\tau_\alpha+{\it t_3}\\\\ 
C>   {\it t_6} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}-2.0\,\tau_\beta\\\\ 
C>   {\it t_{19}} &=& 2.0\,\tau_\beta+{\it t_{17}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{21}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}-2.0\,\tau_s\\\\ 
C>   {\it t_{23}} &=& 2.0\,\tau_s+{\it t_{21}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{16}}\,
C>    \sigma_{\beta\beta}}\over{0.00449276922095889\,{\it t_{16}}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left({{{
C>    \it t_{15}}\,{\it t_{18}}}\over{{\it t_{19}}}}+{{{\it t_{14}}
C>    \,{\it t_{18}}^2}\over{{\it t_{19}}^2}}+{{{\it t_{13}}\,{
C>    \it t_{18}}^3}\over{{\it t_{19}}^3}}+{{{\it t_{12}}\,{
C>    \it t_{18}}^4}\over{{\it t_{19}}^4}}+{{{\it t_{11}}\,{
C>    \it t_{18}}^5}\over{{\it t_{19}}^5}}+{{{\it t_{10}}\,{
C>    \it t_{18}}^6}\over{{\it t_{19}}^6}}+{{{\it t_9}\,{
C>    \it t_{18}}^7}\over{{\it t_{19}}^7}}+{{{\it t_8}\,{
C>    \it t_{18}}^8}\over{{\it t_{19}}^8}}+{{{\it t_7}\,{
C>    \it t_{18}}^9}\over{{\it t_{19}}^9}}+{{{\it t_6}\,{
C>    \it t_{18}}^{10}}\over{{\it t_{19}}^{10}}}+{{{\it t_2}\,{
C>    \it t_{18}}^{11}}\over{{\it t_{19}}^{11}}}\right)+1.0\right)
C>    -0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.00449276922095889\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,
C>    \left({{{\it t_{15}}\,{\it t_4}}\over{{\it t_5}}}+{{{
C>    \it t_{14}}\,{\it t_4}^2}\over{{\it t_5}^2}}+{{{\it t_{13}}
C>    \,{\it t_4}^3}\over{{\it t_5}^3}}+{{{\it t_{12}}\,{
C>    \it t_4}^4}\over{{\it t_5}^4}}+{{{\it t_{11}}\,{
C>    \it t_4}^5}\over{{\it t_5}^5}}+{{{\it t_{10}}\,{
C>    \it t_4}^6}\over{{\it t_5}^6}}+{{{\it t_9}\,{\it t_4}^7}
C>    \over{{\it t_5}^7}}+{{{\it t_8}\,{\it t_4}^8}\over{{\it t_5}^8}}
C>    +{{{\it t_7}\,{\it t_4}^9}\over{{\it t_5}^9}}+{{{\it t_6}\,{
C>    \it t_4}^{10}}\over{{\it t_5}^{10}}}+{{{\it t_2}\,{
C>    \it t_4}^{11}}\over{{\it t_5}^{11}}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{20}}\,\sigma_{ss}}
C>    \over{0.00449276922095889\,{\it t_{20}}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{
C>    \it t_{22}}}\over{{\it t_{23}}}}+{{{\it t_{14}}\,{
C>    \it t_{22}}^2}\over{{\it t_{23}}^2}}+{{{\it t_{13}}\,{
C>    \it t_{22}}^3}\over{{\it t_{23}}^3}}+{{{\it t_{12}}\,{
C>    \it t_{22}}^4}\over{{\it t_{23}}^4}}+{{{\it t_{11}}\,{
C>    \it t_{22}}^5}\over{{\it t_{23}}^5}}+{{{\it t_{10}}\,{
C>    \it t_{22}}^6}\over{{\it t_{23}}^6}}+{{{\it t_9}\,{
C>    \it t_{22}}^7}\over{{\it t_{23}}^7}}+{{{\it t_8}\,{
C>    \it t_{22}}^8}\over{{\it t_{23}}^8}}+{{{\it t_7}\,{
C>    \it t_{22}}^9}\over{{\it t_{23}}^9}}+{{{\it t_6}\,{
C>    \it t_{22}}^{10}}\over{{\it t_{23}}^{10}}}+{{{\it t_2}\,{
C>    \it t_{22}}^{11}}\over{{\it t_{23}}^{11}}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, NE Schultz, DG Truhlar, J.Chem.Phys. 123
C>    , 161103 (2005)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2126975 "> 
C> 10.1063/1.2126975 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_m05(param,tol_rho,ipol,nq,wght, 
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
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
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
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
            t42 = 1/rhob**2.6666666666666666d+0
            t43 = 4.49276922095889d-3*gammabb*t42+1.0d+0
            t44 = 1/t43
            t45 = 3.6121864536509485d-3*gammabb*t42*t44+1.0d+0
            t46 = rhob**1.3333333333333333d+0
            t47 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t48 = t47-2.0d+0*taub
            t49 = t48**11
            t50 = 2.0d+0*taub+t47
            t51 = 1/t50**11
            t52 = t48**10
            t53 = 1/t50**10
            t54 = t48**9
            t55 = 1/t50**9
            t56 = t48**8
            t57 = 1/t50**8
            t58 = t48**7
            t59 = 1/t50**7
            t60 = t48**6
            t61 = 1/t50**6
            t62 = t48**5
            t63 = 1/t50**5
            t64 = t48**4
            t65 = 1/t50**4
            t66 = t48**3
            t67 = 1/t50**3
            t68 = t48**2
            t69 = 1/t50**2
            t70 = 1/t50
            t71 = 1.0d+0*(t39*t48*t70+t36*t68*t69+t33*t66*t67+t30*t64*t6
     1         5+t27*t62*t63+t24*t60*t61+t49*t51*t6+t21*t58*t59+t18*t56*
     2         t57+t15*t54*t55+t12*t52*t53)+1.0d+0
            t72 = rhoa**6.666666666666666d-1
            t73 = 1/t10**12
            t74 = 1/t2**2
            t75 = rhob**6.666666666666666d-1
            t76 = 1/t50**12
            t77 = 1/t43**2
            fnc(iq) = (-9.305257363490997d-1*t45*t46*t71-9.3052573634909
     1         97d-1*t4*t41*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-1.6711932865
     1         267184d+2*t6*t72*t73*t9-1.5192666241151986d+1*t38*t39*t72
     2         *t8+3.038533248230397d+1*t36*t38*t72*t8+1.671193286526718
     3         4d+2*t11*t13*t6*t72+1.5192666241151986d+1*t39*t40*t72-3.0
     4         38533248230397d+1*t35*t36*t37*t72+4.5577998723455965d+1*t
     5         33*t35*t37*t72-4.5577998723455965d+1*t32*t33*t34*t72+6.07
     6         7066496460795d+1*t30*t32*t34*t72-6.077066496460795d+1*t29
     7         *t30*t31*t72+7.596333120575993d+1*t27*t29*t31*t72-7.59633
     8         3120575993d+1*t26*t27*t28*t72+9.115599744691193d+1*t24*t2
     9         6*t28*t72-9.115599744691193d+1*t23*t24*t25*t72+1.06348663
     :         68806391d+2*t21*t23*t25*t72-1.0634866368806391d+2*t20*t21
     ;         *t22*t72+1.2154132992921589d+2*t18*t20*t22*t72-1.21541329
     <         92921589d+2*t17*t18*t19*t72+1.3673399617036788d+2*t15*t17
     =         *t19*t72-1.3673399617036788d+2*t14*t15*t16*t72+1.51926662
     >         41151986d+2*t12*t14*t16*t72-1.5192666241151986d+2*t11*t12
     ?         *t13*t72)-9.305257363490997d-1*t41*t5*(4.327658698487368d
     @         -5*gammaaa**2*t74/rhoa**6.333333333333333d+0-9.6324972097
     1         35863d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)-1.240700
     2         9817987995d+0*rhoa**3.333333333333333d-1*t4*t41)*wght+Ama
     3         t(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t46*t71*(4.327658698
     1         487368d-5*gammabb**2*t77/rhob**6.333333333333333d+0-9.632
     2         497209735863d-3*gammabb*t44/rhob**3.6666666666666664d+0)-
     3         9.305257363490997d-1*t45*t46*(-1.6711932865267184d+2*t49*
     4         t6*t75*t76+1.5192666241151986d+1*t39*t70*t75-1.5192666241
     5         151986d+1*t39*t48*t69*t75+3.038533248230397d+1*t36*t48*t6
     6         9*t75-3.038533248230397d+1*t36*t67*t68*t75+4.557799872345
     7         5965d+1*t33*t67*t68*t75-4.5577998723455965d+1*t33*t65*t66
     8         *t75+6.077066496460795d+1*t30*t65*t66*t75-6.0770664964607
     9         95d+1*t30*t63*t64*t75+7.596333120575993d+1*t27*t63*t64*t7
     :         5-7.596333120575993d+1*t27*t61*t62*t75+9.115599744691193d
     ;         +1*t24*t61*t62*t75-9.115599744691193d+1*t24*t59*t60*t75+1
     <         .0634866368806391d+2*t21*t59*t60*t75+1.6711932865267184d+
     =         2*t51*t52*t6*t75-1.0634866368806391d+2*t21*t57*t58*t75+1.
     >         2154132992921589d+2*t18*t57*t58*t75-1.2154132992921589d+2
     ?         *t18*t55*t56*t75+1.3673399617036788d+2*t15*t55*t56*t75-1.
     @         3673399617036788d+2*t15*t53*t54*t75+1.5192666241151986d+2
     1         *t12*t53*t54*t75-1.5192666241151986d+2*t12*t51*t52*t75)-1
     2         .2407009817987995d+0*rhob**3.333333333333333d-1*t45*t71)*
     3         wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t41*t
     1         5*(3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     2         aaa*t74/rhoa**5.333333333333333d+0)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t46*t
     1         71*(3.6121864536509485d-3*t42*t44-1.6228720119327628d-5*g
     2         ammabb*t77/rhob**5.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         -2.2d+1*t6*t73*t9+t38*(-2.0d+0*t39*t8-4.0d+0*t36*t8)+t11*
     2         (-2.2d+1*t13*t6-2.0d+1*t12*t13)-2.0d+0*t39*t40+t35*(-4.0d
     3         +0*t36*t37-6.0d+0*t33*t37)+t32*(-6.0d+0*t33*t34-8.0d+0*t3
     4         0*t34)+t29*(-8.0d+0*t30*t31-1.0d+1*t27*t31)+t26*(-1.0d+1*
     5         t27*t28-1.2d+1*t24*t28)+t23*(-1.2d+1*t24*t25-1.4d+1*t21*t
     6         25)+t20*(-1.4d+1*t21*t22-1.6d+1*t18*t22)+t17*(-1.6d+1*t18
     7         *t19-1.8d+1*t15*t19)+t14*(-1.8d+1*t15*t16-2.0d+1*t12*t16)
     8         )*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t45*t46
     1         *(-2.2d+1*t49*t6*t76-2.0d+0*t39*t70+(-2.0d+0*t39*t48-4.0d
     2         +0*t36*t48)*t69+t67*(-4.0d+0*t36*t68-6.0d+0*t33*t68)+t65*
     3         (-6.0d+0*t33*t66-8.0d+0*t30*t66)+t63*(-8.0d+0*t30*t64-1.0
     4         d+1*t27*t64)+t61*(-1.0d+1*t27*t62-1.2d+1*t24*t62)+t59*(-1
     5         .2d+1*t24*t60-1.4d+1*t21*t60)+t51*(-2.2d+1*t52*t6-2.0d+1*
     6         t12*t52)+t57*(-1.4d+1*t21*t58-1.6d+1*t18*t58)+t55*(-1.6d+
     7         1*t18*t56-1.8d+1*t15*t56)+t53*(-1.8d+1*t15*t54-2.0d+1*t12
     8         *t54))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
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
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t8 = t7-2.0d+0*taub
            t9 = t8**11
            t10 = 2.0d+0*taub+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
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
C> \brief Evaluate the nwxcm_x_m05 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_3} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_4} &=& {\it t_3}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_5} &=& 2.0\,\tau_\alpha+{\it t_3}\\\\ 
C>   {\it t_6} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}-2.0\,\tau_\beta\\\\ 
C>   {\it t_{19}} &=& 2.0\,\tau_\beta+{\it t_{17}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{21}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}-2.0\,\tau_s\\\\ 
C>   {\it t_{23}} &=& 2.0\,\tau_s+{\it t_{21}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{16}}\,
C>    \sigma_{\beta\beta}}\over{0.00449276922095889\,{\it t_{16}}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left({{{
C>    \it t_{15}}\,{\it t_{18}}}\over{{\it t_{19}}}}+{{{\it t_{14}}
C>    \,{\it t_{18}}^2}\over{{\it t_{19}}^2}}+{{{\it t_{13}}\,{
C>    \it t_{18}}^3}\over{{\it t_{19}}^3}}+{{{\it t_{12}}\,{
C>    \it t_{18}}^4}\over{{\it t_{19}}^4}}+{{{\it t_{11}}\,{
C>    \it t_{18}}^5}\over{{\it t_{19}}^5}}+{{{\it t_{10}}\,{
C>    \it t_{18}}^6}\over{{\it t_{19}}^6}}+{{{\it t_9}\,{
C>    \it t_{18}}^7}\over{{\it t_{19}}^7}}+{{{\it t_8}\,{
C>    \it t_{18}}^8}\over{{\it t_{19}}^8}}+{{{\it t_7}\,{
C>    \it t_{18}}^9}\over{{\it t_{19}}^9}}+{{{\it t_6}\,{
C>    \it t_{18}}^{10}}\over{{\it t_{19}}^{10}}}+{{{\it t_2}\,{
C>    \it t_{18}}^{11}}\over{{\it t_{19}}^{11}}}\right)+1.0\right)
C>    -0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.00449276922095889\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,
C>    \left({{{\it t_{15}}\,{\it t_4}}\over{{\it t_5}}}+{{{
C>    \it t_{14}}\,{\it t_4}^2}\over{{\it t_5}^2}}+{{{\it t_{13}}
C>    \,{\it t_4}^3}\over{{\it t_5}^3}}+{{{\it t_{12}}\,{
C>    \it t_4}^4}\over{{\it t_5}^4}}+{{{\it t_{11}}\,{
C>    \it t_4}^5}\over{{\it t_5}^5}}+{{{\it t_{10}}\,{
C>    \it t_4}^6}\over{{\it t_5}^6}}+{{{\it t_9}\,{\it t_4}^7}
C>    \over{{\it t_5}^7}}+{{{\it t_8}\,{\it t_4}^8}\over{{\it t_5}^8}}
C>    +{{{\it t_7}\,{\it t_4}^9}\over{{\it t_5}^9}}+{{{\it t_6}\,{
C>    \it t_4}^{10}}\over{{\it t_5}^{10}}}+{{{\it t_2}\,{
C>    \it t_4}^{11}}\over{{\it t_5}^{11}}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{20}}\,\sigma_{ss}}
C>    \over{0.00449276922095889\,{\it t_{20}}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{
C>    \it t_{22}}}\over{{\it t_{23}}}}+{{{\it t_{14}}\,{
C>    \it t_{22}}^2}\over{{\it t_{23}}^2}}+{{{\it t_{13}}\,{
C>    \it t_{22}}^3}\over{{\it t_{23}}^3}}+{{{\it t_{12}}\,{
C>    \it t_{22}}^4}\over{{\it t_{23}}^4}}+{{{\it t_{11}}\,{
C>    \it t_{22}}^5}\over{{\it t_{23}}^5}}+{{{\it t_{10}}\,{
C>    \it t_{22}}^6}\over{{\it t_{23}}^6}}+{{{\it t_9}\,{
C>    \it t_{22}}^7}\over{{\it t_{23}}^7}}+{{{\it t_8}\,{
C>    \it t_{22}}^8}\over{{\it t_{23}}^8}}+{{{\it t_7}\,{
C>    \it t_{22}}^9}\over{{\it t_{23}}^9}}+{{{\it t_6}\,{
C>    \it t_{22}}^{10}}\over{{\it t_{23}}^{10}}}+{{{\it t_2}\,{
C>    \it t_{22}}^{11}}\over{{\it t_{23}}^{11}}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, NE Schultz, DG Truhlar, J.Chem.Phys. 123
C>    , 161103 (2005)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2126975 "> 
C> 10.1063/1.2126975 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_m05_d2(param,tol_rho,ipol,nq,wght, 
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
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
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
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
            t42 = 1/rhob**2.6666666666666666d+0
            t43 = 4.49276922095889d-3*gammabb*t42+1.0d+0
            t44 = 1/t43
            t45 = 3.6121864536509485d-3*gammabb*t42*t44+1.0d+0
            t46 = rhob**1.3333333333333333d+0
            t47 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t48 = t47-2.0d+0*taub
            t49 = t48**11
            t50 = 2.0d+0*taub+t47
            t51 = 1/t50**11
            t52 = t48**10
            t53 = 1/t50**10
            t54 = t48**9
            t55 = 1/t50**9
            t56 = t48**8
            t57 = 1/t50**8
            t58 = t48**7
            t59 = 1/t50**7
            t60 = t48**6
            t61 = 1/t50**6
            t62 = t48**5
            t63 = 1/t50**5
            t64 = t48**4
            t65 = 1/t50**4
            t66 = t48**3
            t67 = 1/t50**3
            t68 = t48**2
            t69 = 1/t50**2
            t70 = 1/t50
            t71 = 1.0d+0*(t39*t48*t70+t36*t68*t69+t33*t66*t67+t30*t64*t6
     1         5+t27*t62*t63+t24*t60*t61+t49*t51*t6+t21*t58*t59+t18*t56*
     2         t57+t15*t54*t55+t12*t52*t53)+1.0d+0
            t72 = rhoa**6.666666666666666d-1
            t73 = 1/t10**12
            t74 = -1.6711932865267184d+2*t6*t72*t73*t9-1.519266624115198
     1         6d+1*t38*t39*t72*t8+3.038533248230397d+1*t36*t38*t72*t8+1
     2         .6711932865267184d+2*t11*t13*t6*t72+1.5192666241151986d+1
     3         *t39*t40*t72-3.038533248230397d+1*t35*t36*t37*t72+4.55779
     4         98723455965d+1*t33*t35*t37*t72-4.5577998723455965d+1*t32*
     5         t33*t34*t72+6.077066496460795d+1*t30*t32*t34*t72-6.077066
     6         496460795d+1*t29*t30*t31*t72+7.596333120575993d+1*t27*t29
     7         *t31*t72-7.596333120575993d+1*t26*t27*t28*t72+9.115599744
     8         691193d+1*t24*t26*t28*t72-9.115599744691193d+1*t23*t24*t2
     9         5*t72+1.0634866368806391d+2*t21*t23*t25*t72-1.06348663688
     :         06391d+2*t20*t21*t22*t72+1.2154132992921589d+2*t18*t20*t2
     ;         2*t72-1.2154132992921589d+2*t17*t18*t19*t72+1.36733996170
     <         36788d+2*t15*t17*t19*t72-1.3673399617036788d+2*t14*t15*t1
     =         6*t72+1.5192666241151986d+2*t12*t14*t16*t72-1.51926662411
     >         51986d+2*t11*t12*t13*t72
            t75 = rhoa**3.333333333333333d-1
            t76 = gammaaa**2
            t77 = 1/t2**2
            t78 = 1/rhoa**6.333333333333333d+0
            t79 = 1/rhoa**3.6666666666666664d+0
            t80 = 4.327658698487368d-5*t76*t77*t78-9.632497209735863d-3*
     1         gammaaa*t3*t79
            t81 = rhob**6.666666666666666d-1
            t82 = 1/t50**12
            t83 = -1.6711932865267184d+2*t49*t6*t81*t82+1.51926662411519
     1         86d+1*t39*t70*t81-1.5192666241151986d+1*t39*t48*t69*t81+3
     2         .038533248230397d+1*t36*t48*t69*t81-3.038533248230397d+1*
     3         t36*t67*t68*t81+4.5577998723455965d+1*t33*t67*t68*t81-4.5
     4         577998723455965d+1*t33*t65*t66*t81+6.077066496460795d+1*t
     5         30*t65*t66*t81-6.077066496460795d+1*t30*t63*t64*t81+7.596
     6         333120575993d+1*t27*t63*t64*t81-7.596333120575993d+1*t27*
     7         t61*t62*t81+9.115599744691193d+1*t24*t61*t62*t81-9.115599
     8         744691193d+1*t24*t59*t60*t81+1.0634866368806391d+2*t21*t5
     9         9*t60*t81+1.6711932865267184d+2*t51*t52*t6*t81-1.06348663
     :         68806391d+2*t21*t57*t58*t81+1.2154132992921589d+2*t18*t57
     ;         *t58*t81-1.2154132992921589d+2*t18*t55*t56*t81+1.36733996
     <         17036788d+2*t15*t55*t56*t81-1.3673399617036788d+2*t15*t53
     =         *t54*t81+1.5192666241151986d+2*t12*t53*t54*t81-1.51926662
     >         41151986d+2*t12*t51*t52*t81
            t84 = rhob**3.333333333333333d-1
            t85 = gammabb**2
            t86 = 1/t43**2
            t87 = 1/rhob**6.333333333333333d+0
            t88 = 1/rhob**3.6666666666666664d+0
            t89 = 4.327658698487368d-5*t85*t86*t87-9.632497209735863d-3*
     1         gammabb*t44*t88
            t90 = 1/rhoa**5.333333333333333d+0
            t91 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t77*t90
            t92 = 1/rhob**5.333333333333333d+0
            t93 = 3.6121864536509485d-3*t42*t44-1.6228720119327628d-5*ga
     1         mmabb*t86*t92
            t94 = -2.2d+1*t6*t73*t9
            t95 = -2.2d+1*t13*t6-2.0d+1*t12*t13
            t96 = -1.8d+1*t15*t16-2.0d+1*t12*t16
            t97 = -1.6d+1*t18*t19-1.8d+1*t15*t19
            t98 = -1.4d+1*t21*t22-1.6d+1*t18*t22
            t99 = -1.2d+1*t24*t25-1.4d+1*t21*t25
            t100 = -1.0d+1*t27*t28-1.2d+1*t24*t28
            t101 = -8.0d+0*t30*t31-1.0d+1*t27*t31
            t102 = -6.0d+0*t33*t34-8.0d+0*t30*t34
            t103 = -4.0d+0*t36*t37-6.0d+0*t33*t37
            t104 = -2.0d+0*t39*t8-4.0d+0*t36*t8
            t105 = -2.0d+0*t39*t40
            t106 = t23*t99+t20*t98+t17*t97+t14*t96+t11*t95+t94+t104*t38+
     1         t103*t35+t102*t32+t101*t29+t100*t26+t105
            t107 = -2.2d+1*t49*t6*t82
            t108 = -2.2d+1*t52*t6-2.0d+1*t12*t52
            t109 = -1.8d+1*t15*t54-2.0d+1*t12*t54
            t110 = -1.6d+1*t18*t56-1.8d+1*t15*t56
            t111 = -1.4d+1*t21*t58-1.6d+1*t18*t58
            t112 = -1.2d+1*t24*t60-1.4d+1*t21*t60
            t113 = -1.0d+1*t27*t62-1.2d+1*t24*t62
            t114 = -8.0d+0*t30*t64-1.0d+1*t27*t64
            t115 = -6.0d+0*t33*t66-8.0d+0*t30*t66
            t116 = -4.0d+0*t36*t68-6.0d+0*t33*t68
            t117 = -2.0d+0*t39*t48-4.0d+0*t36*t48
            t118 = -2.0d+0*t39*t70
            t119 = t117*t69+t116*t67+t115*t65+t114*t63+t113*t61+t112*t59
     1         +t111*t57+t110*t55+t109*t53+t108*t51+t118+t107
            t120 = 1/t10**13
            t121 = 1/t75
            t122 = 1/t2**3
            t123 = 1/t50**13
            t124 = 1/t84
            t125 = 1/t43**3
            t126 = t94-2.0d+0*t38*t39*t8-4.0d+0*t36*t38*t8-2.2d+1*t11*t1
     1         3*t6-4.0d+0*t35*t36*t37-6.0d+0*t33*t35*t37-6.0d+0*t32*t33
     2         *t34-8.0d+0*t30*t32*t34-8.0d+0*t29*t30*t31-1.0d+1*t27*t29
     3         *t31-1.0d+1*t26*t27*t28-1.2d+1*t24*t26*t28-1.2d+1*t23*t24
     4         *t25-1.4d+1*t21*t23*t25-1.4d+1*t20*t21*t22-1.6d+1*t18*t20
     5         *t22-1.6d+1*t17*t18*t19-1.8d+1*t15*t17*t19-1.8d+1*t14*t15
     6         *t16-2.0d+1*t12*t14*t16-2.0d+1*t11*t12*t13+t105
            t127 = -2.0d+0*t39*t48*t69-4.0d+0*t36*t48*t69-4.0d+0*t36*t67
     1         *t68-6.0d+0*t33*t67*t68-6.0d+0*t33*t65*t66-8.0d+0*t30*t65
     2         *t66-8.0d+0*t30*t63*t64-1.0d+1*t27*t63*t64-1.0d+1*t27*t61
     3         *t62-1.2d+1*t24*t61*t62-1.2d+1*t24*t59*t60-1.4d+1*t21*t59
     4         *t60-2.2d+1*t51*t52*t6-1.4d+1*t21*t57*t58-1.6d+1*t18*t57*
     5         t58-1.6d+1*t18*t55*t56-1.8d+1*t15*t55*t56-1.8d+1*t15*t53*
     6         t54-2.0d+1*t12*t53*t54-2.0d+1*t12*t51*t52+t118+t107
            t128 = 8.0d+0*t39+8.0d+0*t36
            fnc(iq) = (-9.305257363490997d-1*t45*t46*t71-9.3052573634909
     1         97d-1*t4*t41*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t41*t5*t80-1.2407009
     1         817987995d+0*t4*t41*t75-9.305257363490997d-1*t4*t5*t74)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t46*t71*t89-1.240700
     1         9817987995d+0*t45*t71*t84-9.305257363490997d-1*t45*t46*t8
     2         3)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t41*t
     1         5*t91*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t46*t
     1         71*t93*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t106*t4
     1         *t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t119*t4
     1         5*t46*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(-1.114128
     1         8576844789d+2*t121*t6*t73*t9+3.046785819198518d+4*t120*t5
     2         *t6*t9+4.6163421503007845d+2*t35*t39*t5*t8-1.846536860120
     3         314d+3*t35*t36*t5*t8+1.3849026450902355d+3*t33*t35*t5*t8-
     4         1.0128444160767991d+1*t121*t38*t39*t8+2.0256888321535982d
     5         +1*t121*t36*t38*t8-5.585774001863949d+4*t13*t5*t6*t73+2.5
     6         389881826654315d+4*t12*t13*t5*t73+2.5389881826654315d+4*t
     7         11*t16*t5*t6+1.1141288576844789d+2*t11*t121*t13*t6-4.6163
     8         421503007845d+2*t38*t39*t5+4.6163421503007845d+2*t36*t38*
     9         t5+1.3849026450902355d+3*t32*t36*t37*t5-4.154707935270706
     :         7d+3*t32*t33*t37*t5+2.769805290180471d+3*t30*t32*t37*t5+2
     ;         .769805290180471d+3*t29*t33*t34*t5-7.386147440481255d+3*t
     <         29*t30*t34*t5+4.616342150300784d+3*t27*t29*t34*t5+4.61634
     =         2150300784d+3*t26*t30*t31*t5-1.154085537575196d+4*t26*t27
     >         *t31*t5+6.924513225451178d+3*t24*t26*t31*t5+6.92451322545
     ?         1178d+3*t23*t27*t28*t5-1.661883174108283d+4*t23*t24*t28*t
     @         5+9.694318515631648d+3*t21*t23*t28*t5+9.694318515631648d+
     1         3*t20*t24*t25*t5-2.2620076536473846d+4*t20*t21*t25*t5+1.2
     2         925758020842196d+4*t18*t20*t25*t5+1.2925758020842196d+4*t
     3         17*t21*t22*t5-2.954458976192502d+4*t17*t18*t22*t5+1.66188
     4         31741082823d+4*t15*t17*t22*t5+1.6618831741082823d+4*t14*t
     5         18*t19*t5-3.7392371417436354d+4*t14*t15*t19*t5+2.07735396
     6         76353528d+4*t12*t14*t19*t5+2.0773539676353528d+4*t11*t15*
     7         t16*t5-4.616342150300784d+4*t11*t12*t16*t5+1.012844416076
     8         7991d+1*t121*t39*t40-2.0256888321535982d+1*t121*t35*t36*t
     9         37+3.038533248230397d+1*t121*t33*t35*t37-3.03853324823039
     :         7d+1*t121*t32*t33*t34+4.0513776643071964d+1*t121*t30*t32*
     ;         t34-4.0513776643071964d+1*t121*t29*t30*t31+5.064222080383
     <         994d+1*t121*t27*t29*t31-5.064222080383994d+1*t121*t26*t27
     =         *t28+6.077066496460795d+1*t121*t24*t26*t28-6.077066496460
     >         795d+1*t121*t23*t24*t25+7.089910912537593d+1*t121*t21*t23
     ?         *t25-7.089910912537593d+1*t121*t20*t21*t22+8.102755328614
     @         393d+1*t121*t18*t20*t22-8.102755328614393d+1*t121*t17*t18
     1         *t19+9.115599744691191d+1*t121*t15*t17*t19-9.115599744691
     2         191d+1*t121*t14*t15*t16+1.012844416076799d+2*t12*t121*t14
     3         *t16-1.012844416076799d+2*t11*t12*t121*t13)-2.48140196359
     4         7599d+0*t41*t75*t80-1.8610514726981994d+0*t5*t74*t80-9.30
     5         5257363490997d-1*t41*t5*(-3.8948928286386303d-4*t76*t77/r
     6         hoa**7.333333333333333d+0+3.531915643569816d-2*gammaaa*t3
     7         /rhoa**4.666666666666667d+0+1.0369691626335495d-6*gammaaa
     8         **3*t122/rhoa**10)-2.481401963597599d+0*t4*t74*t75-4.1356
     9         699393293317d-1*t4*t41/t72)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-2.481401963597599d+0*t71*t84*t89-1.86
     1         10514726981994d+0*t46*t83*t89-9.305257363490997d-1*t46*t7
     2         1*(-3.8948928286386303d-4*t85*t86/rhob**7.333333333333333
     3         d+0+3.531915643569816d-2*gammabb*t44/rhob**4.666666666666
     4         667d+0+1.0369691626335495d-6*gammabb**3*t125/rhob**10)-2.
     5         481401963597599d+0*t45*t83*t84-9.305257363490997d-1*t45*t
     6         46*(-5.585774001863949d+4*t46*t52*t6*t82-1.11412885768447
     7         89d+2*t124*t49*t6*t82+2.5389881826654315d+4*t12*t46*t52*t
     8         82+1.0128444160767991d+1*t124*t39*t70-1.0128444160767991d
     9         +1*t124*t39*t48*t69+2.0256888321535982d+1*t124*t36*t48*t6
     :         9-4.6163421503007845d+2*t39*t46*t69+4.6163421503007845d+2
     ;         *t36*t46*t69-2.0256888321535982d+1*t124*t36*t67*t68+3.038
     <         533248230397d+1*t124*t33*t67*t68+1.3849026450902355d+3*t3
     =         6*t46*t65*t68-4.1547079352707067d+3*t33*t46*t65*t68+2.769
     >         805290180471d+3*t30*t46*t65*t68+4.6163421503007845d+2*t39
     ?         *t46*t48*t67-1.846536860120314d+3*t36*t46*t48*t67+1.38490
     @         26450902355d+3*t33*t46*t48*t67-3.038533248230397d+1*t124*
     1         t33*t65*t66+4.0513776643071964d+1*t124*t30*t65*t66+2.7698
     2         05290180471d+3*t33*t46*t63*t66-7.386147440481255d+3*t30*t
     3         46*t63*t66+4.616342150300784d+3*t27*t46*t63*t66-4.0513776
     4         643071964d+1*t124*t30*t63*t64+5.064222080383994d+1*t124*t
     5         27*t63*t64+4.616342150300784d+3*t30*t46*t61*t64-1.1540855
     6         37575196d+4*t27*t46*t61*t64+6.924513225451178d+3*t24*t46*
     7         t61*t64-5.064222080383994d+1*t124*t27*t61*t62+6.077066496
     8         460795d+1*t124*t24*t61*t62+6.924513225451178d+3*t27*t46*t
     9         59*t62-1.661883174108283d+4*t24*t46*t59*t62+9.69431851563
     :         1648d+3*t21*t46*t59*t62-6.077066496460795d+1*t124*t24*t59
     ;         *t60+7.089910912537593d+1*t124*t21*t59*t60+9.694318515631
     <         648d+3*t24*t46*t57*t60-2.2620076536473846d+4*t21*t46*t57*
     =         t60+1.2925758020842196d+4*t18*t46*t57*t60+2.5389881826654
     >         315d+4*t46*t51*t54*t6+1.1141288576844789d+2*t124*t51*t52*
     ?         t6+3.046785819198518d+4*t123*t46*t49*t6-7.089910912537593
     @         d+1*t124*t21*t57*t58+8.102755328614393d+1*t124*t18*t57*t5
     1         8+1.2925758020842196d+4*t21*t46*t55*t58-2.954458976192502
     2         d+4*t18*t46*t55*t58+1.6618831741082823d+4*t15*t46*t55*t58
     3         -8.102755328614393d+1*t124*t18*t55*t56+9.115599744691191d
     4         +1*t124*t15*t55*t56+1.6618831741082823d+4*t18*t46*t53*t56
     5         -3.7392371417436354d+4*t15*t46*t53*t56+2.0773539676353528
     6         d+4*t12*t46*t53*t56-9.115599744691191d+1*t124*t15*t53*t54
     7         +1.012844416076799d+2*t12*t124*t53*t54+2.0773539676353528
     8         d+4*t15*t46*t51*t54-4.616342150300784d+4*t12*t46*t51*t54-
     9         1.012844416076799d+2*t12*t124*t51*t52)-4.1356699393293317
     :         d-1*t45*t71/t81)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817987995d+0*t41*t75*t91-9.
     1         305257363490997d-1*t5*t74*t91-9.305257363490997d-1*t41*t5
     2         *(-9.632497209735863d-3*t3*t79+1.29829760954621d-4*gammaa
     3         a*t77*t78-3.8886343598758105d-7*t122*t76/rhoa**9))*wght+C
     4         mat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.2407009817987995d+0*t71*t84*t93-9.
     1         305257363490997d-1*t46*t83*t93-9.305257363490997d-1*t46*t
     2         71*(-9.632497209735863d-3*t44*t88+1.29829760954621d-4*gam
     3         mabb*t86*t87-3.8886343598758105d-7*t125*t85/rhob**9))*wgh
     4         t+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t41*t5*(1.4582378849534292d-7*gammaaa*t122/rhoa**8-3
     2         .2457440238655255d-5*t77*t90)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t46*t71*(1.4582378849534292d-7*gammabb*t125/rhob**8-
     2         3.2457440238655255d-5*t86*t92)*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(4.0108638
     1         876641245d+3*t120*t6*t72*t9+6.077066496460795d+1*t35*t39*
     2         t72*t8-1.8231199489382385d+2*t33*t35*t72*t8+3.34238657305
     3         3437d+3*t12*t13*t72*t73-3.342386573053437d+3*t11*t16*t6*t
     4         72-6.077066496460795d+1*t36*t38*t72+1.8231199489382385d+2
     5         *t32*t36*t37*t72-3.646239897876477d+2*t30*t32*t37*t72+3.6
     6         46239897876477d+2*t29*t33*t34*t72-6.077066496460795d+2*t2
     7         7*t29*t34*t72+6.077066496460795d+2*t26*t30*t31*t72-9.1155
     8         99744691192d+2*t24*t26*t31*t72+9.115599744691192d+2*t23*t
     9         27*t28*t72-1.276183964256767d+3*t21*t23*t28*t72+1.2761839
     :         64256767d+3*t20*t24*t25*t72-1.7015786190090226d+3*t18*t20
     ;         *t25*t72+1.7015786190090226d+3*t17*t21*t22*t72-2.18774393
     <         8725886d+3*t15*t17*t22*t72+2.187743938725886d+3*t14*t18*t
     =         19*t72-2.734679923407357d+3*t12*t14*t19*t72+2.73467992340
     >         7357d+3*t11*t15*t16*t72)-9.305257363490997d-1*t126*t5*t80
     ?         -1.2407009817987995d+0*t126*t4*t75)*wght+Mmat2(iq,D2_RA_T
     @         A)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t127*t46*t89-1.2
     1         407009817987995d+0*t127*t45*t84-9.305257363490997d-1*t45*
     2         t46*(3.342386573053437d+3*t12*t52*t81*t82-6.0770664964607
     3         95d+1*t36*t69*t81+1.8231199489382385d+2*t36*t65*t68*t81-3
     4         .646239897876477d+2*t30*t65*t68*t81+6.077066496460795d+1*
     5         t39*t48*t67*t81-1.8231199489382385d+2*t33*t48*t67*t81+3.6
     6         46239897876477d+2*t33*t63*t66*t81-6.077066496460795d+2*t2
     7         7*t63*t66*t81+6.077066496460795d+2*t30*t61*t64*t81-9.1155
     8         99744691192d+2*t24*t61*t64*t81+9.115599744691192d+2*t27*t
     9         59*t62*t81-1.276183964256767d+3*t21*t59*t62*t81+1.2761839
     :         64256767d+3*t24*t57*t60*t81-1.7015786190090226d+3*t18*t57
     ;         *t60*t81-3.342386573053437d+3*t51*t54*t6*t81+4.0108638876
     <         641245d+3*t123*t49*t6*t81+1.7015786190090226d+3*t21*t55*t
     =         58*t81-2.187743938725886d+3*t15*t55*t58*t81+2.18774393872
     >         5886d+3*t18*t53*t56*t81-2.734679923407357d+3*t12*t53*t56*
     ?         t81+2.734679923407357d+3*t15*t51*t54*t81))*wght+Mmat2(iq,
     @         D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t106*t5*t91*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t119*t46*t93*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(t20*(-1.4d+1*t99+1.96d+2*t21*t25+2.24d+2*t18*t25)
     2         +t17*(-1.6d+1*t98+2.56d+2*t18*t22+2.88d+2*t15*t22)+t14*(-
     3         1.8d+1*t97+3.24d+2*t15*t19+3.6d+2*t12*t19)+t11*(-2.0d+1*t
     4         96+4.4d+2*t16*t6+4.0d+2*t12*t16)+t73*(4.84d+2*t13*t6-2.2d
     5         +1*t95)+5.28d+2*t120*t6*t9+t35*(1.6d+1*t36*t8+2.4d+1*t33*
     6         t8-4.0d+0*t104)+t128*t38+t32*(3.6d+1*t33*t37+4.8d+1*t30*t
     7         37-6.0d+0*t103)+t29*(6.4d+1*t30*t34+8.0d+1*t27*t34-8.0d+0
     8         *t102)+t26*(1.0d+2*t27*t31+1.2d+2*t24*t31-1.0d+1*t101)+t2
     9         3*(1.44d+2*t24*t28+1.68d+2*t21*t28-1.2d+1*t100))*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t45*t46*((4.84d+2*t52*t6-2.2d+1*t108)*t82+t128*t69+t65*(
     2         3.6d+1*t33*t68+4.8d+1*t30*t68-6.0d+0*t116)+(1.6d+1*t36*t4
     3         8+2.4d+1*t33*t48-4.0d+0*t117)*t67+t63*(6.4d+1*t30*t66+8.0
     4         d+1*t27*t66-8.0d+0*t115)+t61*(1.0d+2*t27*t64+1.2d+2*t24*t
     5         64-1.0d+1*t114)+t59*(1.44d+2*t24*t62+1.68d+2*t21*t62-1.2d
     6         +1*t113)+t57*(1.96d+2*t21*t60+2.24d+2*t18*t60-1.4d+1*t112
     7         )+t51*(4.4d+2*t54*t6+4.0d+2*t12*t54-2.0d+1*t109)+5.28d+2*
     8         t123*t49*t6+t55*(2.56d+2*t18*t58+2.88d+2*t15*t58-1.6d+1*t
     9         111)+t53*(3.24d+2*t15*t56+3.6d+2*t12*t56-1.8d+1*t110))*wg
     :         ht
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t8 = t7-2.0d+0*taua
            t9 = t8**11
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
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
            t6 = param(11)
            t7 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t8 = t7-2.0d+0*taub
            t9 = t8**11
            t10 = 2.0d+0*taub+t7
            t11 = 1/t10**11
            t12 = param(10)
            t13 = t8**10
            t14 = 1/t10**10
            t15 = param(9)
            t16 = t8**9
            t17 = 1/t10**9
            t18 = param(8)
            t19 = t8**8
            t20 = 1/t10**8
            t21 = param(7)
            t22 = t8**7
            t23 = 1/t10**7
            t24 = param(6)
            t25 = t8**6
            t26 = 1/t10**6
            t27 = param(5)
            t28 = t8**5
            t29 = 1/t10**5
            t30 = param(4)
            t31 = t8**4
            t32 = 1/t10**4
            t33 = param(3)
            t34 = t8**3
            t35 = 1/t10**3
            t36 = param(2)
            t37 = t8**2
            t38 = 1/t10**2
            t39 = param(1)
            t40 = 1/t10
            t41 = 1.0d+0*(t11*t6*t9+t39*t40*t8+t36*t37*t38+t33*t34*t35+t
     1         30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t18*t19*t2
     2         0+t15*t16*t17+t12*t13*t14)+1.0d+0
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
C> \brief Evaluate the nwxcm_x_m05 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_3} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_4} &=& {\it t_3}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_5} &=& 2.0\,\tau_\alpha+{\it t_3}\\\\ 
C>   {\it t_6} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}-2.0\,\tau_\beta\\\\ 
C>   {\it t_{19}} &=& 2.0\,\tau_\beta+{\it t_{17}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{21}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}-2.0\,\tau_s\\\\ 
C>   {\it t_{23}} &=& 2.0\,\tau_s+{\it t_{21}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{16}}\,
C>    \sigma_{\beta\beta}}\over{0.00449276922095889\,{\it t_{16}}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left({{{
C>    \it t_{15}}\,{\it t_{18}}}\over{{\it t_{19}}}}+{{{\it t_{14}}
C>    \,{\it t_{18}}^2}\over{{\it t_{19}}^2}}+{{{\it t_{13}}\,{
C>    \it t_{18}}^3}\over{{\it t_{19}}^3}}+{{{\it t_{12}}\,{
C>    \it t_{18}}^4}\over{{\it t_{19}}^4}}+{{{\it t_{11}}\,{
C>    \it t_{18}}^5}\over{{\it t_{19}}^5}}+{{{\it t_{10}}\,{
C>    \it t_{18}}^6}\over{{\it t_{19}}^6}}+{{{\it t_9}\,{
C>    \it t_{18}}^7}\over{{\it t_{19}}^7}}+{{{\it t_8}\,{
C>    \it t_{18}}^8}\over{{\it t_{19}}^8}}+{{{\it t_7}\,{
C>    \it t_{18}}^9}\over{{\it t_{19}}^9}}+{{{\it t_6}\,{
C>    \it t_{18}}^{10}}\over{{\it t_{19}}^{10}}}+{{{\it t_2}\,{
C>    \it t_{18}}^{11}}\over{{\it t_{19}}^{11}}}\right)+1.0\right)
C>    -0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.00449276922095889\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,
C>    \left({{{\it t_{15}}\,{\it t_4}}\over{{\it t_5}}}+{{{
C>    \it t_{14}}\,{\it t_4}^2}\over{{\it t_5}^2}}+{{{\it t_{13}}
C>    \,{\it t_4}^3}\over{{\it t_5}^3}}+{{{\it t_{12}}\,{
C>    \it t_4}^4}\over{{\it t_5}^4}}+{{{\it t_{11}}\,{
C>    \it t_4}^5}\over{{\it t_5}^5}}+{{{\it t_{10}}\,{
C>    \it t_4}^6}\over{{\it t_5}^6}}+{{{\it t_9}\,{\it t_4}^7}
C>    \over{{\it t_5}^7}}+{{{\it t_8}\,{\it t_4}^8}\over{{\it t_5}^8}}
C>    +{{{\it t_7}\,{\it t_4}^9}\over{{\it t_5}^9}}+{{{\it t_6}\,{
C>    \it t_4}^{10}}\over{{\it t_5}^{10}}}+{{{\it t_2}\,{
C>    \it t_4}^{11}}\over{{\it t_5}^{11}}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003612186453650949\,{\it t_{20}}\,\sigma_{ss}}
C>    \over{0.00449276922095889\,{\it t_{20}}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{
C>    \it t_{22}}}\over{{\it t_{23}}}}+{{{\it t_{14}}\,{
C>    \it t_{22}}^2}\over{{\it t_{23}}^2}}+{{{\it t_{13}}\,{
C>    \it t_{22}}^3}\over{{\it t_{23}}^3}}+{{{\it t_{12}}\,{
C>    \it t_{22}}^4}\over{{\it t_{23}}^4}}+{{{\it t_{11}}\,{
C>    \it t_{22}}^5}\over{{\it t_{23}}^5}}+{{{\it t_{10}}\,{
C>    \it t_{22}}^6}\over{{\it t_{23}}^6}}+{{{\it t_9}\,{
C>    \it t_{22}}^7}\over{{\it t_{23}}^7}}+{{{\it t_8}\,{
C>    \it t_{22}}^8}\over{{\it t_{23}}^8}}+{{{\it t_7}\,{
C>    \it t_{22}}^9}\over{{\it t_{23}}^9}}+{{{\it t_6}\,{
C>    \it t_{22}}^{10}}\over{{\it t_{23}}^{10}}}+{{{\it t_2}\,{
C>    \it t_{22}}^{11}}\over{{\it t_{23}}^{11}}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, NE Schultz, DG Truhlar, J.Chem.Phys. 123
C>    , 161103 (2005)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2126975 "> 
C> 10.1063/1.2126975 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_m05_d3(param,tol_rho,ipol,nq,wght, 
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
            t6 = param(11)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(10)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(9)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(8)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(7)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(6)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(5)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(4)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(3)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(2)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(1)
            t41 = 1/t11
            t42 = 1.0d+0*(t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+
     1         t31*t32*t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t
     2         21+t16*t17*t18+t13*t14*t15)+1.0d+0
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
            t6 = param(11)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(10)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(9)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(8)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(7)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(6)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(5)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(4)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(3)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(2)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(1)
            t41 = 1/t11
            t42 = 1.0d+0*(t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+
     1         t31*t32*t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t
     2         21+t16*t17*t18+t13*t14*t15)+1.0d+0
            t43 = 1/rhob**2.6666666666666666d+0
            t44 = 4.49276922095889d-3*gammabb*t43+1.0d+0
            t45 = 1/t44
            t46 = 3.6121864536509485d-3*gammabb*t43*t45+1.0d+0
            t47 = rhob**1.3333333333333333d+0
            t48 = rhob**1.6666666666666669d+0
            t49 = 9.115599744691192d+0*t48
            t50 = t49-2.0d+0*taub
            t51 = t50**11
            t52 = 2.0d+0*taub+t49
            t53 = 1/t52**11
            t54 = t50**10
            t55 = 1/t52**10
            t56 = t50**9
            t57 = 1/t52**9
            t58 = t50**8
            t59 = 1/t52**8
            t60 = t50**7
            t61 = 1/t52**7
            t62 = t50**6
            t63 = 1/t52**6
            t64 = t50**5
            t65 = 1/t52**5
            t66 = t50**4
            t67 = 1/t52**4
            t68 = t50**3
            t69 = 1/t52**3
            t70 = t50**2
            t71 = 1/t52**2
            t72 = 1/t52
            t73 = 1.0d+0*(t40*t50*t72+t37*t70*t71+t34*t68*t69+t31*t66*t6
     1         7+t28*t64*t65+t25*t62*t63+t22*t60*t61+t51*t53*t6+t19*t58*
     2         t59+t16*t56*t57+t13*t54*t55)+1.0d+0
            t74 = rhoa**6.666666666666666d-1
            t75 = 1/t11**12
            t76 = -1.5192666241151986d+1*t39*t40*t74*t9+3.03853324823039
     1         7d+1*t37*t39*t74*t9-1.6711932865267184d+2*t10*t6*t74*t75+
     2         1.6711932865267184d+2*t12*t14*t6*t74+1.5192666241151986d+
     3         1*t40*t41*t74-3.038533248230397d+1*t36*t37*t38*t74+4.5577
     4         998723455965d+1*t34*t36*t38*t74-4.5577998723455965d+1*t33
     5         *t34*t35*t74+6.077066496460795d+1*t31*t33*t35*t74-6.07706
     6         6496460795d+1*t30*t31*t32*t74+7.596333120575993d+1*t28*t3
     7         0*t32*t74-7.596333120575993d+1*t27*t28*t29*t74+9.11559974
     8         4691193d+1*t25*t27*t29*t74-9.115599744691193d+1*t24*t25*t
     9         26*t74+1.0634866368806391d+2*t22*t24*t26*t74-1.0634866368
     :         806391d+2*t21*t22*t23*t74+1.2154132992921589d+2*t19*t21*t
     ;         23*t74-1.2154132992921589d+2*t18*t19*t20*t74+1.3673399617
     <         036788d+2*t16*t18*t20*t74-1.3673399617036788d+2*t15*t16*t
     =         17*t74+1.5192666241151986d+2*t13*t15*t17*t74-1.5192666241
     >         151986d+2*t12*t13*t14*t74
            t77 = rhoa**3.333333333333333d-1
            t78 = gammaaa**2
            t79 = 1/t2**2
            t80 = 1/rhoa**6.333333333333333d+0
            t81 = 1/rhoa**3.6666666666666664d+0
            t82 = 4.327658698487368d-5*t78*t79*t80-9.632497209735863d-3*
     1         gammaaa*t3*t81
            t83 = rhob**6.666666666666666d-1
            t84 = 1/t52**12
            t85 = -1.6711932865267184d+2*t51*t6*t83*t84+1.51926662411519
     1         86d+1*t40*t72*t83-1.5192666241151986d+1*t40*t50*t71*t83+3
     2         .038533248230397d+1*t37*t50*t71*t83-3.038533248230397d+1*
     3         t37*t69*t70*t83+4.5577998723455965d+1*t34*t69*t70*t83-4.5
     4         577998723455965d+1*t34*t67*t68*t83+6.077066496460795d+1*t
     5         31*t67*t68*t83-6.077066496460795d+1*t31*t65*t66*t83+7.596
     6         333120575993d+1*t28*t65*t66*t83-7.596333120575993d+1*t28*
     7         t63*t64*t83+9.115599744691193d+1*t25*t63*t64*t83-9.115599
     8         744691193d+1*t25*t61*t62*t83+1.0634866368806391d+2*t22*t6
     9         1*t62*t83-1.0634866368806391d+2*t22*t59*t60*t83+1.2154132
     :         992921589d+2*t19*t59*t60*t83+1.6711932865267184d+2*t53*t5
     ;         4*t6*t83-1.2154132992921589d+2*t19*t57*t58*t83+1.36733996
     <         17036788d+2*t16*t57*t58*t83-1.3673399617036788d+2*t16*t55
     =         *t56*t83+1.5192666241151986d+2*t13*t55*t56*t83-1.51926662
     >         41151986d+2*t13*t53*t54*t83
            t86 = rhob**3.333333333333333d-1
            t87 = gammabb**2
            t88 = 1/t44**2
            t89 = 1/rhob**6.333333333333333d+0
            t90 = 1/rhob**3.6666666666666664d+0
            t91 = 4.327658698487368d-5*t87*t88*t89-9.632497209735863d-3*
     1         gammabb*t45*t90
            t92 = 1/rhoa**5.333333333333333d+0
            t93 = 3.6121864536509485d-3*t1*t3-1.6228720119327628d-5*gamm
     1         aaa*t79*t92
            t94 = 1/rhob**5.333333333333333d+0
            t95 = 3.6121864536509485d-3*t43*t45-1.6228720119327628d-5*ga
     1         mmabb*t88*t94
            t96 = -2.2d+1*t10*t6*t75
            t97 = -2.2d+1*t14*t6-2.0d+1*t13*t14
            t98 = -1.8d+1*t16*t17-2.0d+1*t13*t17
            t99 = -1.6d+1*t19*t20-1.8d+1*t16*t20
            t100 = -1.4d+1*t22*t23-1.6d+1*t19*t23
            t101 = -1.2d+1*t25*t26-1.4d+1*t22*t26
            t102 = -1.0d+1*t28*t29-1.2d+1*t25*t29
            t103 = -8.0d+0*t31*t32-1.0d+1*t28*t32
            t104 = -6.0d+0*t34*t35-8.0d+0*t31*t35
            t105 = -4.0d+0*t37*t38-6.0d+0*t34*t38
            t106 = -2.0d+0*t40*t9-4.0d+0*t37*t9
            t107 = -2.0d+0*t40*t41
            t108 = t18*t99+t15*t98+t12*t97+t96+t106*t39+t105*t36+t104*t3
     1         3+t103*t30+t102*t27+t101*t24+t100*t21+t107
            t109 = -2.2d+1*t51*t6*t84
            t110 = -2.2d+1*t54*t6-2.0d+1*t13*t54
            t111 = -1.8d+1*t16*t56-2.0d+1*t13*t56
            t112 = -1.6d+1*t19*t58-1.8d+1*t16*t58
            t113 = -1.4d+1*t22*t60-1.6d+1*t19*t60
            t114 = -1.2d+1*t25*t62-1.4d+1*t22*t62
            t115 = -1.0d+1*t28*t64-1.2d+1*t25*t64
            t116 = -8.0d+0*t31*t66-1.0d+1*t28*t66
            t117 = -6.0d+0*t34*t68-8.0d+0*t31*t68
            t118 = -4.0d+0*t37*t70-6.0d+0*t34*t70
            t119 = -2.0d+0*t40*t50-4.0d+0*t37*t50
            t120 = -2.0d+0*t40*t72
            t121 = t119*t71+t118*t69+t117*t67+t116*t65+t115*t63+t114*t61
     1         +t113*t59+t112*t57+t111*t55+t110*t53+t120+t109
            t122 = 1/t11**13
            t123 = 1/t77
            t124 = 4.6163421503007845d+2*t36*t40*t5*t9-1.846536860120314
     1         d+3*t36*t37*t5*t9+1.3849026450902355d+3*t34*t36*t5*t9-1.0
     2         128444160767991d+1*t123*t39*t40*t9+2.0256888321535982d+1*
     3         t123*t37*t39*t9-5.585774001863949d+4*t14*t5*t6*t75-1.1141
     4         288576844789d+2*t10*t123*t6*t75+2.5389881826654315d+4*t13
     5         *t14*t5*t75+2.5389881826654315d+4*t12*t17*t5*t6+3.0467858
     6         19198518d+4*t10*t122*t5*t6+1.1141288576844789d+2*t12*t123
     7         *t14*t6-4.6163421503007845d+2*t39*t40*t5+4.61634215030078
     8         45d+2*t37*t39*t5+1.3849026450902355d+3*t33*t37*t38*t5-4.1
     9         547079352707067d+3*t33*t34*t38*t5+2.769805290180471d+3*t3
     :         1*t33*t38*t5+2.769805290180471d+3*t30*t34*t35*t5-7.386147
     ;         440481255d+3*t30*t31*t35*t5+4.616342150300784d+3*t28*t30*
     <         t35*t5+4.616342150300784d+3*t27*t31*t32*t5-1.154085537575
     =         196d+4*t27*t28*t32*t5+6.924513225451178d+3*t25*t27*t32*t5
     >         +6.924513225451178d+3*t24*t28*t29*t5-1.661883174108283d+4
     ?         *t24*t25*t29*t5+9.694318515631648d+3*t22*t24*t29*t5+9.694
     @         318515631648d+3*t21*t25*t26*t5-2.2620076536473846d+4*t21*
     1         t22*t26*t5+1.2925758020842196d+4*t19*t21*t26*t5+1.2925758
     2         020842196d+4*t18*t22*t23*t5-2.954458976192502d+4*t18*t19*
     3         t23*t5+1.6618831741082823d+4*t16*t18*t23*t5+1.66188317410
     4         82823d+4*t15*t19*t20*t5-3.7392371417436354d+4*t15*t16*t20
     5         *t5+2.0773539676353528d+4*t13*t15*t20*t5+2.07735396763535
     6         28d+4*t12*t16*t17*t5-4.616342150300784d+4*t12*t13*t17*t5+
     7         1.0128444160767991d+1*t123*t40*t41-2.0256888321535982d+1*
     8         t123*t36*t37*t38+3.038533248230397d+1*t123*t34*t36*t38-3.
     9         038533248230397d+1*t123*t33*t34*t35+4.0513776643071964d+1
     :         *t123*t31*t33*t35-4.0513776643071964d+1*t123*t30*t31*t32+
     ;         5.064222080383994d+1*t123*t28*t30*t32-5.064222080383994d+
     <         1*t123*t27*t28*t29+6.077066496460795d+1*t123*t25*t27*t29-
     =         6.077066496460795d+1*t123*t24*t25*t26+7.089910912537593d+
     >         1*t123*t22*t24*t26-7.089910912537593d+1*t123*t21*t22*t23+
     ?         8.102755328614393d+1*t123*t19*t21*t23-8.102755328614393d+
     @         1*t123*t18*t19*t20+9.115599744691191d+1*t123*t16*t18*t20-
     1         9.115599744691191d+1*t123*t15*t16*t17+1.012844416076799d+
     2         2*t123*t13*t15*t17-1.012844416076799d+2*t12*t123*t13*t14
            t125 = 1/t74
            t126 = gammaaa**3
            t127 = 1/t2**3
            t128 = 1/rhoa**10
            t129 = 1/rhoa**7.333333333333333d+0
            t130 = 1/rhoa**4.666666666666667d+0
            t131 = -3.8948928286386303d-4*t129*t78*t79+3.531915643569816
     1         d-2*gammaaa*t130*t3+1.0369691626335495d-6*t126*t127*t128
            t132 = 1/t52**13
            t133 = 1/t86
            t134 = -5.585774001863949d+4*t47*t54*t6*t84-1.11412885768447
     1         89d+2*t133*t51*t6*t84+2.5389881826654315d+4*t13*t47*t54*t
     2         84+1.0128444160767991d+1*t133*t40*t72-1.0128444160767991d
     3         +1*t133*t40*t50*t71+2.0256888321535982d+1*t133*t37*t50*t7
     4         1-4.6163421503007845d+2*t40*t47*t71+4.6163421503007845d+2
     5         *t37*t47*t71-2.0256888321535982d+1*t133*t37*t69*t70+3.038
     6         533248230397d+1*t133*t34*t69*t70+1.3849026450902355d+3*t3
     7         7*t47*t67*t70-4.1547079352707067d+3*t34*t47*t67*t70+2.769
     8         805290180471d+3*t31*t47*t67*t70+4.6163421503007845d+2*t40
     9         *t47*t50*t69-1.846536860120314d+3*t37*t47*t50*t69+1.38490
     :         26450902355d+3*t34*t47*t50*t69-3.038533248230397d+1*t133*
     ;         t34*t67*t68+4.0513776643071964d+1*t133*t31*t67*t68+2.7698
     <         05290180471d+3*t34*t47*t65*t68-7.386147440481255d+3*t31*t
     =         47*t65*t68+4.616342150300784d+3*t28*t47*t65*t68-4.0513776
     >         643071964d+1*t133*t31*t65*t66+5.064222080383994d+1*t133*t
     ?         28*t65*t66+4.616342150300784d+3*t31*t47*t63*t66-1.1540855
     @         37575196d+4*t28*t47*t63*t66+6.924513225451178d+3*t25*t47*
     1         t63*t66-5.064222080383994d+1*t133*t28*t63*t64+6.077066496
     2         460795d+1*t133*t25*t63*t64+6.924513225451178d+3*t28*t47*t
     3         61*t64-1.661883174108283d+4*t25*t47*t61*t64+9.69431851563
     4         1648d+3*t22*t47*t61*t64-6.077066496460795d+1*t133*t25*t61
     5         *t62+7.089910912537593d+1*t133*t22*t61*t62+9.694318515631
     6         648d+3*t25*t47*t59*t62-2.2620076536473846d+4*t22*t47*t59*
     7         t62+1.2925758020842196d+4*t19*t47*t59*t62-7.0899109125375
     8         93d+1*t133*t22*t59*t60+8.102755328614393d+1*t133*t19*t59*
     9         t60+1.2925758020842196d+4*t22*t47*t57*t60-2.9544589761925
     :         02d+4*t19*t47*t57*t60+1.6618831741082823d+4*t16*t47*t57*t
     ;         60+2.5389881826654315d+4*t47*t53*t56*t6+1.114128857684478
     <         9d+2*t133*t53*t54*t6+3.046785819198518d+4*t132*t47*t51*t6
     =         -8.102755328614393d+1*t133*t19*t57*t58+9.115599744691191d
     >         +1*t133*t16*t57*t58+1.6618831741082823d+4*t19*t47*t55*t58
     ?         -3.7392371417436354d+4*t16*t47*t55*t58+2.0773539676353528
     @         d+4*t13*t47*t55*t58-9.115599744691191d+1*t133*t16*t55*t56
     1         +1.012844416076799d+2*t13*t133*t55*t56+2.0773539676353528
     2         d+4*t16*t47*t53*t56-4.616342150300784d+4*t13*t47*t53*t56-
     3         1.012844416076799d+2*t13*t133*t53*t54
            t135 = 1/t83
            t136 = gammabb**3
            t137 = 1/t44**3
            t138 = 1/rhob**10
            t139 = 1/rhob**7.333333333333333d+0
            t140 = 1/rhob**4.666666666666667d+0
            t141 = -3.8948928286386303d-4*t139*t87*t88+3.531915643569816
     1         d-2*gammabb*t140*t45+1.0369691626335495d-6*t136*t137*t138
            t142 = 1/rhoa**9
            t143 = -9.632497209735863d-3*t3*t81+1.29829760954621d-4*gamm
     1         aaa*t79*t80-3.8886343598758105d-7*t127*t142*t78
            t144 = 1/rhob**9
            t145 = -9.632497209735863d-3*t45*t90+1.29829760954621d-4*gam
     1         mabb*t88*t89-3.8886343598758105d-7*t137*t144*t87
            t146 = 1/rhoa**8
            t147 = 1.4582378849534292d-7*gammaaa*t127*t146-3.24574402386
     1         55255d-5*t79*t92
            t148 = 1/rhob**8
            t149 = 1.4582378849534292d-7*gammabb*t137*t148-3.24574402386
     1         55255d-5*t88*t94
            t150 = 6.077066496460795d+1*t36*t40*t74*t9-1.823119948938238
     1         5d+2*t34*t36*t74*t9+3.342386573053437d+3*t13*t14*t74*t75-
     2         3.342386573053437d+3*t12*t17*t6*t74+4.0108638876641245d+3
     3         *t10*t122*t6*t74-6.077066496460795d+1*t37*t39*t74+1.82311
     4         99489382385d+2*t33*t37*t38*t74-3.646239897876477d+2*t31*t
     5         33*t38*t74+3.646239897876477d+2*t30*t34*t35*t74-6.0770664
     6         96460795d+2*t28*t30*t35*t74+6.077066496460795d+2*t27*t31*
     7         t32*t74-9.115599744691192d+2*t25*t27*t32*t74+9.1155997446
     8         91192d+2*t24*t28*t29*t74-1.276183964256767d+3*t22*t24*t29
     9         *t74+1.276183964256767d+3*t21*t25*t26*t74-1.7015786190090
     :         226d+3*t19*t21*t26*t74+1.7015786190090226d+3*t18*t22*t23*
     ;         t74-2.187743938725886d+3*t16*t18*t23*t74+2.18774393872588
     <         6d+3*t15*t19*t20*t74-2.734679923407357d+3*t13*t15*t20*t74
     =         +2.734679923407357d+3*t12*t16*t17*t74
            t151 = t96-2.0d+0*t39*t40*t9-4.0d+0*t37*t39*t9-2.2d+1*t12*t1
     1         4*t6-4.0d+0*t36*t37*t38-6.0d+0*t34*t36*t38-6.0d+0*t33*t34
     2         *t35-8.0d+0*t31*t33*t35-8.0d+0*t30*t31*t32-1.0d+1*t28*t30
     3         *t32-1.0d+1*t27*t28*t29-1.2d+1*t25*t27*t29-1.2d+1*t24*t25
     4         *t26-1.4d+1*t22*t24*t26-1.4d+1*t21*t22*t23-1.6d+1*t19*t21
     5         *t23-1.6d+1*t18*t19*t20-1.8d+1*t16*t18*t20-1.8d+1*t15*t16
     6         *t17-2.0d+1*t13*t15*t17-2.0d+1*t12*t13*t14+t107
            t152 = 3.342386573053437d+3*t13*t54*t83*t84-6.07706649646079
     1         5d+1*t37*t71*t83+1.8231199489382385d+2*t37*t67*t70*t83-3.
     2         646239897876477d+2*t31*t67*t70*t83+6.077066496460795d+1*t
     3         40*t50*t69*t83-1.8231199489382385d+2*t34*t50*t69*t83+3.64
     4         6239897876477d+2*t34*t65*t68*t83-6.077066496460795d+2*t28
     5         *t65*t68*t83+6.077066496460795d+2*t31*t63*t66*t83-9.11559
     6         9744691192d+2*t25*t63*t66*t83+9.115599744691192d+2*t28*t6
     7         1*t64*t83-1.276183964256767d+3*t22*t61*t64*t83+1.27618396
     8         4256767d+3*t25*t59*t62*t83-1.7015786190090226d+3*t19*t59*
     9         t62*t83+1.7015786190090226d+3*t22*t57*t60*t83-2.187743938
     :         725886d+3*t16*t57*t60*t83-3.342386573053437d+3*t53*t56*t6
     ;         *t83+4.0108638876641245d+3*t132*t51*t6*t83+2.187743938725
     <         886d+3*t19*t55*t58*t83-2.734679923407357d+3*t13*t55*t58*t
     =         83+2.734679923407357d+3*t16*t53*t56*t83
            t153 = -2.0d+0*t40*t50*t71-4.0d+0*t37*t50*t71-4.0d+0*t37*t69
     1         *t70-6.0d+0*t34*t69*t70-6.0d+0*t34*t67*t68-8.0d+0*t31*t67
     2         *t68-8.0d+0*t31*t65*t66-1.0d+1*t28*t65*t66-1.0d+1*t28*t63
     3         *t64-1.2d+1*t25*t63*t64-1.2d+1*t25*t61*t62-1.4d+1*t22*t61
     4         *t62-1.4d+1*t22*t59*t60-1.6d+1*t19*t59*t60-2.2d+1*t53*t54
     5         *t6-1.6d+1*t19*t57*t58-1.8d+1*t16*t57*t58-1.8d+1*t16*t55*
     6         t56-2.0d+1*t13*t55*t56-2.0d+1*t13*t53*t54+t120+t109
            t154 = 4.84d+2*t14*t6-2.2d+1*t97
            t155 = 4.0d+2*t13*t17
            t156 = 4.4d+2*t17*t6
            t157 = -2.0d+1*t98+t156+t155
            t158 = 3.24d+2*t16*t20
            t159 = 3.6d+2*t13*t20
            t160 = -1.8d+1*t99+t159+t158
            t161 = 2.56d+2*t19*t23
            t162 = 2.88d+2*t16*t23
            t163 = t162+t161-1.6d+1*t100
            t164 = 1.96d+2*t22*t26
            t165 = 2.24d+2*t19*t26
            t166 = t165+t164-1.4d+1*t101
            t167 = 1.44d+2*t25*t29
            t168 = 1.68d+2*t22*t29
            t169 = t168+t167-1.2d+1*t102
            t170 = 1.0d+2*t28*t32
            t171 = 1.2d+2*t25*t32
            t172 = t171+t170-1.0d+1*t103
            t173 = 6.4d+1*t31*t35
            t174 = 8.0d+1*t28*t35
            t175 = t174+t173-8.0d+0*t104
            t176 = 3.6d+1*t34*t38
            t177 = 4.8d+1*t31*t38
            t178 = t177+t176-6.0d+0*t105
            t179 = 1.6d+1*t37*t9
            t180 = 2.4d+1*t34*t9
            t181 = t180+t179-4.0d+0*t106
            t182 = 8.0d+0*t37
            t183 = 8.0d+0*t40+t182
            t184 = t154*t75+5.28d+2*t10*t122*t6+t183*t39+t181*t36+t178*t
     1         33+t175*t30+t172*t27+t169*t24+t166*t21+t163*t18+t15*t160+
     2         t12*t157
            t185 = 4.84d+2*t54*t6-2.2d+1*t110
            t186 = 4.0d+2*t13*t56
            t187 = 4.4d+2*t56*t6
            t188 = t187+t186-2.0d+1*t111
            t189 = 3.24d+2*t16*t58
            t190 = 3.6d+2*t13*t58
            t191 = t190+t189-1.8d+1*t112
            t192 = 2.56d+2*t19*t60
            t193 = 2.88d+2*t16*t60
            t194 = t193+t192-1.6d+1*t113
            t195 = 1.96d+2*t22*t62
            t196 = 2.24d+2*t19*t62
            t197 = t196+t195-1.4d+1*t114
            t198 = 1.44d+2*t25*t64
            t199 = 1.68d+2*t22*t64
            t200 = t199+t198-1.2d+1*t115
            t201 = 1.0d+2*t28*t66
            t202 = 1.2d+2*t25*t66
            t203 = t202+t201-1.0d+1*t116
            t204 = 6.4d+1*t31*t68
            t205 = 8.0d+1*t28*t68
            t206 = t205+t204-8.0d+0*t117
            t207 = 3.6d+1*t34*t70
            t208 = 4.8d+1*t31*t70
            t209 = t208+t207-6.0d+0*t118
            t210 = 1.6d+1*t37*t50
            t211 = 2.4d+1*t34*t50
            t212 = t211+t210-4.0d+0*t119
            t213 = t185*t84+t183*t71+t212*t69+t209*t67+t206*t65+t203*t63
     1         +t200*t61+5.28d+2*t132*t51*t6+t197*t59+t194*t57+t191*t55+
     2         t188*t53
            t214 = rhoa**2
            t215 = 1/t11**14
            t216 = 1/t5
            t217 = 1/t2**4
            t218 = rhob**2
            t219 = 1/t52**14
            t220 = 1/t47
            t221 = 1/t44**4
            t222 = -4.0d+0*(4.0d+0*t40+t182)-3.2d+1*t37-4.8d+1*t34-4.0d+
     1         0*t183
            fnc(iq) = (-9.305257363490997d-1*t46*t47*t73-9.3052573634909
     1         97d-1*t4*t42*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t42*t5*t82-1.2407009
     1         817987995d+0*t4*t42*t77-9.305257363490997d-1*t4*t5*t76)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t47*t73*t91-1.240700
     1         9817987995d+0*t46*t73*t86-9.305257363490997d-1*t46*t47*t8
     2         5)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t42*t
     1         5*t93*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t47*t
     1         73*t95*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t108*t4
     1         *t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t121*t4
     1         6*t47*wght
            Amat2(iq,D2_RA_RA) = (-2.481401963597599d+0*t42*t77*t82-1.86
     1         10514726981994d+0*t5*t76*t82-2.481401963597599d+0*t4*t76*
     2         t77-9.305257363490997d-1*t131*t42*t5-9.305257363490997d-1
     3         *t124*t4*t5-4.1356699393293317d-1*t125*t4*t42)*wght+Amat2
     4         (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-2.481401963597599d+0*t73*t86*t91-1.86
     1         10514726981994d+0*t47*t85*t91-2.481401963597599d+0*t46*t8
     2         5*t86-9.305257363490997d-1*t141*t47*t73-4.135669939329331
     3         7d-1*t135*t46*t73-9.305257363490997d-1*t134*t46*t47)*wght
     4         +Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817987995d+0*t42*t77*t93-9.
     1         305257363490997d-1*t5*t76*t93-9.305257363490997d-1*t143*t
     2         42*t5)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.2407009817987995d+0*t73*t86*t95-9.
     1         305257363490997d-1*t47*t85*t95-9.305257363490997d-1*t145*
     2         t47*t73)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t147*t42*t5*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t149*t47*t73*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t151*t5*t82-1.24
     1         07009817987995d+0*t151*t4*t77-9.305257363490997d-1*t150*t
     2         4*t5)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t153*t47*t91-1.2
     1         407009817987995d+0*t153*t46*t86-9.305257363490997d-1*t152
     2         *t46*t47)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t108*t5*t93*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t121*t47*t95*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t184*t4*t5*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t213*t46*t47*wght
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490997d-1*t4*t5*(9.2326
     1         84300601568d+2*t36*t40*t77*t9-3.693073720240627d+3*t36*t3
     2         7*t77*t9+2.769805290180471d+3*t34*t36*t77*t9+3.3761480535
     3         8933d+0*t216*t39*t40*t9-2.104036366334451d+4*t214*t33*t40
     4         *t9-6.75229610717866d+0*t216*t37*t39*t9+1.262421819800670
     5         4d+5*t214*t33*t37*t9-1.8936327297010064d+5*t214*t33*t34*t
     6         9+8.416145465337804d+4*t214*t31*t33*t9-1.1171548003727898
     7         d+5*t14*t6*t75*t77+5.077976365330862d+4*t13*t14*t75*t77+5
     8         .077976365330863d+4*t12*t17*t6*t77+6.093571638397036d+4*t
     9         10*t122*t6*t77-9.232684300601568d+2*t39*t40*t77+9.2326843
     :         00601568d+2*t37*t39*t77+2.769805290180471d+3*t33*t37*t38*
     ;         t77-8.309415870541413d+3*t33*t34*t38*t77+5.53961058036094
     <         2d+3*t31*t33*t38*t77+5.539610580360942d+3*t30*t34*t35*t77
     =         -1.4772294880962508d+4*t30*t31*t35*t77+9.232684300601568d
     >         +3*t28*t30*t35*t77+9.232684300601568d+3*t27*t31*t32*t77-2
     ?         .308171075150392d+4*t27*t28*t32*t77+1.3849026450902355d+4
     @         *t25*t27*t32*t77+1.3849026450902355d+4*t24*t28*t29*t77-3.
     1         323766348216566d+4*t24*t25*t29*t77+1.9388637031263292d+4*
     2         t22*t24*t29*t77+1.9388637031263292d+4*t21*t25*t26*t77-4.5
     3         240153072947686d+4*t21*t22*t26*t77+2.585151604168439d+4*t
     4         19*t21*t26*t77+2.585151604168439d+4*t18*t22*t23*t77-5.908
     5         917952385003d+4*t18*t19*t23*t77+3.3237663482165647d+4*t16
     6         *t18*t23*t77+3.3237663482165647d+4*t15*t19*t20*t77-7.4784
     7         74283487271d+4*t15*t16*t20*t77+4.154707935270705d+4*t13*t
     8         15*t20*t77+4.1547079352707056d+4*t12*t16*t17*t77-9.232684
     9         300601568d+4*t12*t13*t17*t77+3.713762858948263d+1*t10*t21
     :         6*t6*t75-1.2729420016323426d+7*t17*t214*t6*t75-3.47166000
     ;         44518436d+6*t16*t17*t214*t75+1.1572200014839477d+7*t13*t1
     <         7*t214*t75-3.713762858948263d+1*t12*t14*t216*t6-6.0175440
     =         0771653d+6*t10*t214*t215*t6+3.4716600044518436d+6*t12*t20
     >         *t214*t6+1.5275304019588112d+7*t122*t14*t214*t6-3.3761480
     ?         5358933d+0*t216*t40*t41+2.104036366334451d+4*t214*t36*t40
     @         +6.75229610717866d+0*t216*t36*t37*t38-8.416145465337804d+
     1         4*t214*t30*t37*t38-1.0128444160767991d+1*t216*t34*t36*t38
     2         +3.787265459402013d+5*t214*t30*t34*t38-5.049687279202681d
     3         +5*t214*t30*t31*t38+2.1040363663344508d+5*t214*t28*t30*t3
     4         8-4.208072732668902d+4*t214*t36*t37+2.104036366334451d+4*
     5         t214*t34*t36+1.0128444160767991d+1*t216*t33*t34*t35-2.104
     6         036366334451d+5*t214*t27*t34*t35-1.350459221435732d+1*t21
     7         6*t31*t33*t35+8.416145465337803d+5*t214*t27*t31*t35-1.052
     8         0181831672254d+6*t214*t27*t28*t35+4.208072732668902d+5*t2
     9         14*t25*t27*t35+1.350459221435732d+1*t216*t30*t31*t32-4.20
     :         80727326689016d+5*t214*t24*t31*t32-1.6880740267946648d+1*
     ;         t216*t28*t30*t32+1.578027274750838d+6*t214*t24*t28*t32-1.
     <         893632729701006d+6*t214*t24*t25*t32+7.364127282170578d+5*
     =         t214*t22*t24*t32+1.6880740267946648d+1*t216*t27*t28*t29-7
     >         .364127282170578d+5*t21*t214*t28*t29-2.0256888321535982d+
     ?         1*t216*t25*t27*t29+2.6510858215814087d+6*t21*t214*t25*t29
     @         -3.0929334585116436d+6*t21*t214*t22*t29+1.178260365147292
     1         7d+6*t19*t21*t214*t29+2.0256888321535982d+1*t216*t24*t25*
     2         t26-1.1782603651472927d+6*t18*t214*t25*t26-2.363303637512
     3         5308d+1*t216*t22*t24*t26+4.123911278015524d+6*t18*t214*t2
     4         2*t26-4.713041460589171d+6*t18*t19*t214*t26+1.76739054772
     5         0939d+6*t16*t18*t214*t26+2.3633036375125308d+1*t21*t216*t
     6         22*t23-1.7673905477209387d+6*t15*t214*t22*t23-2.700918442
     7         871464d+1*t19*t21*t216*t23+6.059624735043219d+6*t15*t19*t
     8         214*t23-6.817077826923622d+6*t15*t16*t214*t23+2.524843639
     9         601341d+6*t13*t15*t214*t23+2.700918442871464d+1*t18*t19*t
     :         20*t216-3.0385332482303967d+1*t16*t18*t20*t216+3.03853324
     ;         82303967d+1*t15*t16*t17*t216-3.3761480535893296d+1*t13*t1
     <         5*t17*t216+3.3761480535893296d+1*t12*t13*t14*t216-2.52484
     =         3639601341d+6*t12*t19*t20*t214+8.521347283654526d+6*t12*t
     >         16*t20*t214-9.468163648505028d+6*t12*t13*t20*t214-4.62888
     ?         0005935792d+6*t122*t13*t14*t214)+t125*(-1.240700981798799
     @         5d+0*t42*t82-1.2407009817987995d+0*t4*t76)-7.444205890792
     1         797d+0*t76*t77*t82-2.791577209047299d+0*t124*t5*t82-9.305
     2         257363490997d-1*t42*t5*(3.279403591520427d-3*t78*t79/rhoa
     3         **8.333333333333334d+0-1.648227300332581d-1*gammaaa*t3/rh
     4         oa**5.666666666666667d+0+3.7270905095708196d-8*gammaaa**4
     5         *t217/rhoa**1.3666666666666666d+1-1.9702414090037443d-5*t
     6         126*t127/rhoa**11)-3.7221029453963983d+0*t131*t42*t77-3.7
     7         221029453963983d+0*t124*t4*t77-2.791577209047299d+0*t131*
     8         t5*t76+2.757113292886221d-1*t4*t42/t7)*wght+Amat3(iq,D3_R
     9         A_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t135*(-1.2407009817987995d+0*t73*t9
     1         1-1.2407009817987995d+0*t46*t85)-7.444205890792797d+0*t85
     2         *t86*t91-2.791577209047299d+0*t134*t47*t91-9.305257363490
     3         997d-1*t47*t73*(3.279403591520427d-3*t87*t88/rhob**8.3333
     4         33333333334d+0-1.648227300332581d-1*gammabb*t45/rhob**5.6
     5         66666666666667d+0+3.7270905095708196d-8*gammabb**4*t221/r
     6         hob**1.3666666666666666d+1-1.9702414090037443d-5*t136*t13
     7         7/rhob**11)-9.305257363490997d-1*t46*t47*(-1.117154800372
     8         7898d+5*t54*t6*t84*t86+5.077976365330862d+4*t13*t54*t84*t
     9         86-9.232684300601568d+2*t40*t71*t86+9.232684300601568d+2*
     :         t37*t71*t86+2.769805290180471d+3*t37*t67*t70*t86-8.309415
     ;         870541413d+3*t34*t67*t70*t86+5.539610580360942d+3*t31*t67
     <         *t70*t86+9.232684300601568d+2*t40*t50*t69*t86-3.693073720
     =         240627d+3*t37*t50*t69*t86+2.769805290180471d+3*t34*t50*t6
     >         9*t86+5.539610580360942d+3*t34*t65*t68*t86-1.477229488096
     ?         2508d+4*t31*t65*t68*t86+9.232684300601568d+3*t28*t65*t68*
     @         t86+9.232684300601568d+3*t31*t63*t66*t86-2.30817107515039
     1         2d+4*t28*t63*t66*t86+1.3849026450902355d+4*t25*t63*t66*t8
     2         6+1.3849026450902355d+4*t28*t61*t64*t86-3.323766348216566
     3         d+4*t25*t61*t64*t86+1.9388637031263292d+4*t22*t61*t64*t86
     4         +1.9388637031263292d+4*t25*t59*t62*t86-4.5240153072947686
     5         d+4*t22*t59*t62*t86+2.585151604168439d+4*t19*t59*t62*t86+
     6         2.585151604168439d+4*t22*t57*t60*t86-5.908917952385003d+4
     7         *t19*t57*t60*t86+3.3237663482165647d+4*t16*t57*t60*t86+5.
     8         077976365330863d+4*t53*t56*t6*t86+6.093571638397036d+4*t1
     9         32*t51*t6*t86+3.3237663482165647d+4*t19*t55*t58*t86-7.478
     :         474283487271d+4*t16*t55*t58*t86+4.154707935270705d+4*t13*
     ;         t55*t58*t86+4.1547079352707056d+4*t16*t53*t56*t86-9.23268
     <         4300601568d+4*t13*t53*t56*t86-1.2729420016323426d+7*t218*
     =         t56*t6*t84+3.713762858948263d+1*t220*t51*t6*t84-3.4716600
     >         044518436d+6*t16*t218*t56*t84+1.1572200014839477d+7*t13*t
     ?         218*t56*t84-3.37614805358933d+0*t220*t40*t72+3.3761480535
     @         8933d+0*t220*t40*t50*t71-6.75229610717866d+0*t220*t37*t50
     1         *t71+6.75229610717866d+0*t220*t37*t69*t70-1.0128444160767
     2         991d+1*t220*t34*t69*t70-8.416145465337804d+4*t218*t37*t65
     3         *t70+3.787265459402013d+5*t218*t34*t65*t70-5.049687279202
     4         681d+5*t218*t31*t65*t70+2.1040363663344508d+5*t218*t28*t6
     5         5*t70+2.104036366334451d+4*t218*t40*t69-4.208072732668902
     6         d+4*t218*t37*t69+2.104036366334451d+4*t218*t34*t69+1.0128
     7         444160767991d+1*t220*t34*t67*t68-1.350459221435732d+1*t22
     8         0*t31*t67*t68-2.104036366334451d+5*t218*t34*t63*t68+8.416
     9         145465337803d+5*t218*t31*t63*t68-1.0520181831672254d+6*t2
     :         18*t28*t63*t68+4.208072732668902d+5*t218*t25*t63*t68-2.10
     ;         4036366334451d+4*t218*t40*t50*t67+1.2624218198006704d+5*t
     <         218*t37*t50*t67-1.8936327297010064d+5*t218*t34*t50*t67+8.
     =         416145465337804d+4*t218*t31*t50*t67+1.350459221435732d+1*
     >         t220*t31*t65*t66-1.6880740267946648d+1*t220*t28*t65*t66-4
     ?         .2080727326689016d+5*t218*t31*t61*t66+1.578027274750838d+
     @         6*t218*t28*t61*t66-1.893632729701006d+6*t218*t25*t61*t66+
     1         7.364127282170578d+5*t218*t22*t61*t66+1.6880740267946648d
     2         +1*t220*t28*t63*t64-2.0256888321535982d+1*t220*t25*t63*t6
     3         4-7.364127282170578d+5*t218*t28*t59*t64+2.651085821581408
     4         7d+6*t218*t25*t59*t64-3.0929334585116436d+6*t218*t22*t59*
     5         t64+1.1782603651472927d+6*t19*t218*t59*t64+2.025688832153
     6         5982d+1*t220*t25*t61*t62-2.3633036375125308d+1*t22*t220*t
     7         61*t62-1.1782603651472927d+6*t218*t25*t57*t62+4.123911278
     8         015524d+6*t218*t22*t57*t62-4.713041460589171d+6*t19*t218*
     9         t57*t62+1.767390547720939d+6*t16*t218*t57*t62+2.363303637
     :         5125308d+1*t22*t220*t59*t60-2.700918442871464d+1*t19*t220
     ;         *t59*t60-1.7673905477209387d+6*t218*t22*t55*t60+6.0596247
     <         35043219d+6*t19*t218*t55*t60-6.817077826923622d+6*t16*t21
     =         8*t55*t60+2.524843639601341d+6*t13*t218*t55*t60+3.4716600
     >         044518436d+6*t218*t53*t58*t6-3.713762858948263d+1*t220*t5
     ?         3*t54*t6+1.5275304019588112d+7*t132*t218*t54*t6-6.0175440
     @         0771653d+6*t218*t219*t51*t6+2.700918442871464d+1*t19*t220
     1         *t57*t58-3.0385332482303967d+1*t16*t220*t57*t58-2.5248436
     2         39601341d+6*t19*t218*t53*t58+8.521347283654526d+6*t16*t21
     3         8*t53*t58-9.468163648505028d+6*t13*t218*t53*t58+3.0385332
     4         482303967d+1*t16*t220*t55*t56-3.3761480535893296d+1*t13*t
     5         220*t55*t56+3.3761480535893296d+1*t13*t220*t53*t54-4.6288
     6         80005935792d+6*t13*t132*t218*t54)-3.7221029453963983d+0*t
     7         141*t73*t86-3.7221029453963983d+0*t134*t46*t86-2.79157720
     8         9047299d+0*t141*t47*t85+2.757113292886221d-1*t46*t73/t48)
     9         *wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-2.481401963597599d+0*t76*t77*t93-
     1         9.305257363490997d-1*t124*t5*t93-4.1356699393293317d-1*t1
     2         25*t42*t93-9.305257363490997d-1*t42*t5*(-9.37659384672262
     3         8d-4*gammaaa*t129*t79+6.610678411788878d-6*t127*t128*t78+
     4         3.531915643569816d-2*t130*t3-1.3976589410890572d-8*t126*t
     5         217/rhoa**1.2666666666666665d+1)-2.481401963597599d+0*t14
     6         3*t42*t77-1.8610514726981994d+0*t143*t5*t76)*wght+Cmat3(i
     7         q,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-2.481401963597599d+0*t85*t86*t95-
     1         4.1356699393293317d-1*t135*t73*t95-9.305257363490997d-1*t
     2         134*t47*t95-9.305257363490997d-1*t47*t73*(-9.376593846722
     3         628d-4*gammabb*t139*t88+6.610678411788878d-6*t137*t138*t8
     4         7+3.531915643569816d-2*t140*t45-1.3976589410890572d-8*t13
     5         6*t221/rhob**1.2666666666666665d+1)-2.481401963597599d+0*
     6         t145*t73*t86-1.8610514726981994d+0*t145*t47*t85)*wght+Cma
     7         t3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t42*t5*(1.7
     1         31063479394947d-4*t79*t80+5.241221029083964d-9*t217*t78/r
     2         hoa**1.1666666666666665d+1-1.944317179937905d-6*gammaaa*t
     3         127*t142)-1.2407009817987995d+0*t147*t42*t77-9.3052573634
     4         90997d-1*t147*t5*t76)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-9.305257363490997d-1*t47*t73*(1.
     1         731063479394947d-4*t88*t89+5.241221029083964d-9*t221*t87/
     2         rhob**1.1666666666666665d+1-1.944317179937905d-6*gammabb*
     3         t137*t144)-1.2407009817987995d+0*t149*t73*t86-9.305257363
     4         490997d-1*t149*t47*t85)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*(4.374713654860287d-7*t127*t146-1.9654578859
     2         06487d-9*gammaaa*t217/rhoa**1.0666666666666666d+1)*t42*t5
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
     1         363490997d-1*(4.374713654860287d-7*t137*t148-1.9654578859
     2         06487d-9*gammabb*t221/rhob**1.0666666666666666d+1)*t47*t7
     3         3*wght
            Mmat3(iq,D3_RA_RA_TA) = (-9.305257363490997d-1*t4*t5*(-2.769
     1         805290180471d+3*t33*t40*t5*t9+5.539610580360942d+3*t33*t3
     2         7*t5*t9+8.309415870541416d+3*t33*t34*t5*t9-1.107922116072
     3         1883d+4*t31*t33*t5*t9+4.0513776643071964d+1*t123*t36*t40*
     4         t9-1.2154132992921589d+2*t123*t34*t36*t9+5.58577400186394
     5         9d+5*t17*t5*t6*t75-4.5701787287977763d+5*t16*t17*t5*t75+5
     6         .077976365330863d+5*t13*t17*t5*t75+2.2282577153689576d+3*
     7         t123*t13*t14*t75-7.921643129916147d+5*t10*t215*t5*t6-4.57
     8         01787287977763d+5*t12*t20*t5*t6+6.702928802236738d+5*t122
     9         *t14*t5*t6-2.228257715368958d+3*t12*t123*t17*t6+2.6739092
     :         584427493d+3*t10*t122*t123*t6+9.232684300601569d+2*t36*t4
     ;         0*t5-1.1079221160721883d+4*t30*t37*t38*t5+1.6618831741082
     <         832d+4*t30*t34*t38*t5+2.2158442321443766d+4*t30*t31*t38*t
     =         5-2.7698052901804704d+4*t28*t30*t38*t5+1.846536860120314d
     >         +3*t36*t37*t5-2.769805290180471d+3*t34*t36*t5-2.769805290
     ?         180471d+4*t27*t34*t35*t5+3.693073720240629d+4*t27*t31*t35
     @         *t5+4.616342150300784d+4*t27*t28*t35*t5-5.539610580360942
     1         d+4*t25*t27*t35*t5-5.539610580360941d+4*t24*t31*t32*t5+6.
     2         924513225451174d+4*t24*t28*t32*t5+8.309415870541416d+4*t2
     3         4*t25*t32*t5-9.694318515631648d+4*t22*t24*t32*t5-9.694318
     4         515631648d+4*t21*t28*t29*t5+1.1633182218757979d+5*t21*t25
     5         *t29*t5+1.3572045921884307d+5*t21*t22*t29*t5-1.5510909625
     6         010635d+5*t19*t21*t29*t5-1.5510909625010638d+5*t18*t25*t2
     7         6*t5+1.809606122917908d+5*t18*t22*t26*t5+2.06812128333475
     8         15d+5*t18*t19*t26*t5-2.3266364437515952d+5*t16*t18*t26*t5
     9         -2.3266364437515955d+5*t15*t22*t23*t5+2.659013078573252d+
     :         5*t15*t19*t23*t5+2.9913897133949086d+5*t15*t16*t23*t5-3.3
     ;         237663482165647d+5*t13*t15*t23*t5-3.323766348216565d+5*t1
     <         2*t19*t20*t5+3.7392371417436354d+5*t12*t16*t20*t5+4.15470
     =         79352707056d+5*t12*t13*t20*t5-6.093571638397035d+5*t122*t
     >         13*t14*t5-4.0513776643071964d+1*t123*t37*t39+1.2154132992
     ?         921589d+2*t123*t33*t37*t38-2.4308265985843178d+2*t123*t31
     @         *t33*t38+2.4308265985843178d+2*t123*t30*t34*t35-4.0513776
     1         64307196d+2*t123*t28*t30*t35+4.0513776643071964d+2*t123*t
     2         27*t31*t32-6.077066496460795d+2*t123*t25*t27*t32+6.077066
     3         496460793d+2*t123*t24*t28*t29-8.507893095045111d+2*t123*t
     4         22*t24*t29+8.507893095045113d+2*t123*t21*t25*t26-1.134385
     5         7460060149d+3*t123*t19*t21*t26+1.1343857460060149d+3*t123
     6         *t18*t22*t23-1.4584959591505905d+3*t123*t16*t18*t23+1.458
     7         4959591505908d+3*t123*t15*t19*t20-1.8231199489382383d+3*t
     8         123*t13*t15*t20+1.8231199489382383d+3*t12*t123*t16*t17)-2
     9         .481401963597599d+0*t151*t77*t82-1.8610514726981994d+0*t1
     :         50*t5*t82-2.481401963597599d+0*t150*t4*t77-9.305257363490
     ;         997d-1*t131*t151*t5-4.1356699393293317d-1*t125*t151*t4)*w
     <         ght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (-2.481401963597599d+0*t153*t86*t91-
     1         1.8610514726981994d+0*t152*t47*t91-2.481401963597599d+0*t
     2         152*t46*t86-9.305257363490997d-1*t46*t47*(5.5857740018639
     3         49d+5*t47*t56*t6*t84-4.5701787287977763d+5*t16*t47*t56*t8
     4         4+5.077976365330863d+5*t13*t47*t56*t84+2.2282577153689576
     5         d+3*t13*t133*t54*t84-4.0513776643071964d+1*t133*t37*t71+1
     6         .2154132992921589d+2*t133*t37*t67*t70-2.4308265985843178d
     7         +2*t133*t31*t67*t70-1.1079221160721883d+4*t37*t47*t65*t70
     8         +1.6618831741082832d+4*t34*t47*t65*t70+2.2158442321443766
     9         d+4*t31*t47*t65*t70-2.7698052901804704d+4*t28*t47*t65*t70
     :         +4.0513776643071964d+1*t133*t40*t50*t69-1.215413299292158
     ;         9d+2*t133*t34*t50*t69+9.232684300601569d+2*t40*t47*t69+1.
     <         846536860120314d+3*t37*t47*t69-2.769805290180471d+3*t34*t
     =         47*t69+2.4308265985843178d+2*t133*t34*t65*t68-4.051377664
     >         307196d+2*t133*t28*t65*t68-2.769805290180471d+4*t34*t47*t
     ?         63*t68+3.693073720240629d+4*t31*t47*t63*t68+4.61634215030
     @         0784d+4*t28*t47*t63*t68-5.539610580360942d+4*t25*t47*t63*
     1         t68-2.769805290180471d+3*t40*t47*t50*t67+5.53961058036094
     2         2d+3*t37*t47*t50*t67+8.309415870541416d+3*t34*t47*t50*t67
     3         -1.1079221160721883d+4*t31*t47*t50*t67+4.0513776643071964
     4         d+2*t133*t31*t63*t66-6.077066496460795d+2*t133*t25*t63*t6
     5         6-5.539610580360941d+4*t31*t47*t61*t66+6.924513225451174d
     6         +4*t28*t47*t61*t66+8.309415870541416d+4*t25*t47*t61*t66-9
     7         .694318515631648d+4*t22*t47*t61*t66+6.077066496460793d+2*
     8         t133*t28*t61*t64-8.507893095045111d+2*t133*t22*t61*t64-9.
     9         694318515631648d+4*t28*t47*t59*t64+1.1633182218757979d+5*
     :         t25*t47*t59*t64+1.3572045921884307d+5*t22*t47*t59*t64-1.5
     ;         510909625010635d+5*t19*t47*t59*t64+8.507893095045113d+2*t
     <         133*t25*t59*t62-1.1343857460060149d+3*t133*t19*t59*t62-1.
     =         5510909625010638d+5*t25*t47*t57*t62+1.809606122917908d+5*
     >         t22*t47*t57*t62+2.0681212833347515d+5*t19*t47*t57*t62-2.3
     ?         266364437515952d+5*t16*t47*t57*t62+1.1343857460060149d+3*
     @         t133*t22*t57*t60-1.4584959591505905d+3*t133*t16*t57*t60-2
     1         .3266364437515955d+5*t22*t47*t55*t60+2.659013078573252d+5
     2         *t19*t47*t55*t60+2.9913897133949086d+5*t16*t47*t55*t60-3.
     3         3237663482165647d+5*t13*t47*t55*t60-4.5701787287977763d+5
     4         *t47*t53*t58*t6-2.228257715368958d+3*t133*t53*t56*t6+6.70
     5         2928802236738d+5*t132*t47*t54*t6-7.921643129916147d+5*t21
     6         9*t47*t51*t6+2.6739092584427493d+3*t132*t133*t51*t6+1.458
     7         4959591505908d+3*t133*t19*t55*t58-1.8231199489382383d+3*t
     8         13*t133*t55*t58-3.323766348216565d+5*t19*t47*t53*t58+3.73
     9         92371417436354d+5*t16*t47*t53*t58+4.1547079352707056d+5*t
     :         13*t47*t53*t58+1.8231199489382383d+3*t133*t16*t53*t56-6.0
     ;         93571638397035d+5*t13*t132*t47*t54)-9.305257363490997d-1*
     <         t141*t153*t47-4.1356699393293317d-1*t135*t153*t46)*wght+M
     =         mat3(iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RA_GAA_TA) = (-1.2407009817987995d+0*t151*t77*t9
     1         3-9.305257363490997d-1*t150*t5*t93-9.305257363490997d-1*t
     2         143*t151*t5)*wght+Mmat3(iq,D3_RA_GAA_TA)
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
            Mmat3(iq,D3_RB_GBB_TB) = (-1.2407009817987995d+0*t153*t86*t9
     1         5-9.305257363490997d-1*t152*t47*t95-9.305257363490997d-1*
     2         t145*t153*t47)*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)-9.305257363490
     1         997d-1*t108*t147*t5*wght
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
     1         3490997d-1*t121*t149*t47*wght
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)-9.30525736
     1         3490997d-1*t184*t5*t93*wght
            Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
            Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
            Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)-9.3052573634
     1         90997d-1*t213*t47*t95*wght
            Mmat3(iq,D3_GAB_TA_TA) = Mmat3(iq,D3_GAB_TA_TA)-9.3052573634
     1         90997d-1*t4*t5*(t33*(-1.44d+2*t34*t9-1.92d+2*t31*t9-6.0d+
     2         0*t181-6.0d+0*(t180+t179))+(-9.68d+3*t17*t6-2.2d+1*t157-2
     3         .2d+1*(t156+t155))*t75+t12*(-7.92d+3*t20*t6-7.2d+3*t13*t2
     4         0-2.0d+1*t160-2.0d+1*(t159+t158))+t122*(-1.1616d+4*t14*t6
     5         -2.4d+1*t154)-1.3728d+4*t10*t215*t6+t30*(-3.84d+2*t31*t38
     6         -4.8d+2*t28*t38-8.0d+0*t178-8.0d+0*(t177+t176))+t222*t36+
     7         t27*(-8.0d+2*t28*t35-9.6d+2*t25*t35-1.0d+1*t175-1.0d+1*(t
     8         174+t173))+t24*(-1.44d+3*t25*t32-1.68d+3*t22*t32-1.2d+1*t
     9         172-1.2d+1*(t171+t170))+t21*(-2.352d+3*t22*t29-2.688d+3*t
     :         19*t29-1.4d+1*t169-1.4d+1*(t168+t167))+t18*(-3.584d+3*t19
     ;         *t26-4.032d+3*t16*t26-1.6d+1*t166-1.6d+1*(t165+t164))+t15
     <         *(-5.184d+3*t16*t23-5.76d+3*t13*t23-1.8d+1*t163-1.8d+1*(t
     =         162+t161)))*wght
            Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
            Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)-9.3052573634
     1         90997d-1*t46*t47*((-9.68d+3*t56*t6-2.2d+1*t188-2.2d+1*(t1
     2         87+t186))*t84+t65*(-3.84d+2*t31*t70-4.8d+2*t28*t70-8.0d+0
     3         *t209-8.0d+0*(t208+t207))+t222*t69+t63*(-8.0d+2*t28*t68-9
     4         .6d+2*t25*t68-1.0d+1*t206-1.0d+1*(t205+t204))+(-1.44d+2*t
     5         34*t50-1.92d+2*t31*t50-6.0d+0*t212-6.0d+0*(t211+t210))*t6
     6         7+t61*(-1.44d+3*t25*t66-1.68d+3*t22*t66-1.2d+1*t203-1.2d+
     7         1*(t202+t201))+t59*(-2.352d+3*t22*t64-2.688d+3*t19*t64-1.
     8         4d+1*t200-1.4d+1*(t199+t198))+t57*(-3.584d+3*t19*t62-4.03
     9         2d+3*t16*t62-1.6d+1*t197-1.6d+1*(t196+t195))+t55*(-5.184d
     :         +3*t16*t60-5.76d+3*t13*t60-1.8d+1*t194-1.8d+1*(t193+t192)
     ;         )+t53*(-7.92d+3*t58*t6-7.2d+3*t13*t58-2.0d+1*t191-2.0d+1*
     <         (t190+t189))+t132*(-1.1616d+4*t54*t6-2.4d+1*t185)-1.3728d
     =         +4*t219*t51*t6)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.49276922095889d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509485d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(11)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taua
            t10 = t9**11
            t11 = 2.0d+0*taua+t8
            t12 = 1/t11**11
            t13 = param(10)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(9)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(8)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(7)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(6)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(5)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(4)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(3)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(2)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(1)
            t41 = 1/t11
            t42 = 1.0d+0*(t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+
     1         t31*t32*t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t
     2         21+t16*t17*t18+t13*t14*t15)+1.0d+0
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
            t6 = param(11)
            t7 = rhob**1.6666666666666669d+0
            t8 = 9.115599744691192d+0*t7
            t9 = t8-2.0d+0*taub
            t10 = t9**11
            t11 = 2.0d+0*taub+t8
            t12 = 1/t11**11
            t13 = param(10)
            t14 = t9**10
            t15 = 1/t11**10
            t16 = param(9)
            t17 = t9**9
            t18 = 1/t11**9
            t19 = param(8)
            t20 = t9**8
            t21 = 1/t11**8
            t22 = param(7)
            t23 = t9**7
            t24 = 1/t11**7
            t25 = param(6)
            t26 = t9**6
            t27 = 1/t11**6
            t28 = param(5)
            t29 = t9**5
            t30 = 1/t11**5
            t31 = param(4)
            t32 = t9**4
            t33 = 1/t11**4
            t34 = param(3)
            t35 = t9**3
            t36 = 1/t11**3
            t37 = param(2)
            t38 = t9**2
            t39 = 1/t11**2
            t40 = param(1)
            t41 = 1/t11
            t42 = 1.0d+0*(t40*t41*t9+t10*t12*t6+t37*t38*t39+t34*t35*t36+
     1         t31*t32*t33+t28*t29*t30+t25*t26*t27+t22*t23*t24+t19*t20*t
     2         21+t16*t17*t18+t13*t14*t15)+1.0d+0
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
