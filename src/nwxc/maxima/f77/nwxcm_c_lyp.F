C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_lyp.F 
C> The nwxcm_c_lyp functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_lyp functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_1}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.349\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^{{{11}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 0.2533\,{\it t_3}\\\\ 
C>   {\it t_7} &=& e^ {- {\it t_6} }\\\\ 
C>   {\it t_8} &=& 0.349\,{\it t_3}\,{\it t_4}\\\\ 
C>   {\it t_9} &=& {\it t_8}+{\it t_6}-11.0\\\\ 
C>   {\it t_{10}} &=& {\it t_8}+{\it t_6}\\\\ 
C>   {\it t_{11}} &=& -3.0\,{\it t_{10}}\\\\ 
C>   f &=& 1.0\,\left(-0.006491760000000001\,{\it t_5}\,{\it t_4}
C>    \,\left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left({
C>    \it t_{11}}-1.0\,\rho_\beta\,{\it t_2}\,{\it t_9}+1.0\right)
C>    -\rho_\alpha^2\right)\,{\it t_7}\,\sigma_{\beta\beta}
C>    -0.006491760000000001\,{\it t_5}\,{\it t_4}\,
C>    \left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left(47.0
C>    -7.0\,{\it t_{10}}\right)-1.333333333333333\,{
C>    \it t_1}^2\right)\,{\it t_7}\,\sigma_{\alpha\beta}
C>    -0.006491760000000001\,{\it t_5}\,{\it t_4}\,
C>    \left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left({
C>    \it t_{11}}-1.0\,\rho_\alpha\,{\it t_2}\,{\it t_9}+1.0\right)
C>    -\rho_\beta^2\right)\,{\it t_7}\,\sigma_{\alpha\alpha}
C>    -0.2367051431943866\,\rho_\alpha\,\rho_\beta\,{\it t_5}
C>    \,\left(\rho_\beta^{{{8}\over{3}}}+\rho_\alpha^{{{8}
C>    \over{3}}}\right)\,{\it t_4}\,{\it t_7}-0.19672\,\rho_\alpha
C>    \,\rho_\beta\,{\it t_2}\,{\it t_4}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] C Lee, W Yang, RG Parr, Phys.Rev.B 37, 785 (1988)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.37.785 "> 
C> 10.1103/PhysRevB.37.785 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_lyp(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Cmat) 
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
      double precision rgamma(nq,*) !< [Input] The norm of the density 
                                    !< gradients 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma 
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
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**3.333333333333333d-1
            t2 = 1/(2.7700148356845083d-1*t1+1.0d+0)
            t3 = 2.010443432317725d-1*t1
            t4 = exp(-t3)
            t5 = rhoa**3.6666666666666664d+0
            t6 = 1/t5
            t7 = rhoa**2
            t8 = -5.333333333333332d+0*t7
            t9 = 2.7700148356845083d-1*t1*t2
            t10 = t9+t3
            t11 = -t7
            t12 = 7.937005259840998d-1
            t13 = 3.49d-1*t1*t12+1.0d+0
            t14 = 1/t13
            t15 = 1/t13**2
            t16 = 7.874506561842959d-2
            t17 = 2.533d-1*t1*t12
            t18 = 3.49d-1*t1*t12*t14
            t19 = t18+t17
            t20 = 4.7d+1-7.0d+0*t19
            t21 = 3.968502629920499d-1
            t22 = 1/rhoa**1.3333333333333333d+0
            t23 = -1.1633333333333332d-1*t14*t21*t22-8.443333333333334d-
     1         2*t21*t22+1.2788303649853786d-2*t15/rhoa**1.6666666666666
     2         669d+0
            t24 = exp(-t17)
            t25 = t18+t17-1.1d+1
            t26 = -5.0d-1*t25-3.0d+0*t19+1.0d+0
            t27 = 1.111111111111111d-1*rhoa*t26
            t28 = -3.5d+0*t23
            t29 = 1/rhoa
            t30 = 1/rhoa**5
            t31 = t8+1.111111111111111d-1*t20*t7
            t32 = 1.111111111111111d-1*t26*t7+t11
            t33 = rhoa**4.666666666666667d+0
            fnc(iq) = 1.0d+0*(-1.0223881343581931d-3*gammaaa*t2*t4*t6*(1
     1         .111111111111111d-1*t7*(-5.0d-1*(t9+t3-1.1d+1)-3.0d+0*t10
     2         +1.0d+0)+t11)-5.111940671790965d-4*gammaaa*t2*t4*t6*(t8+1
     3         .111111111111111d-1*(4.7d+1-7.0d+0*t10)*t7)-3.72787240661
     4         23486d-2*rhoa*t2*t4-9.836d-2*rhoa*t2)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(t14*t16*t6*(-6.491760000000001d-3*g
     1         ammaaa*t24*(1.111111111111111d-1*(2.5d-1*t25*t29+t28)*t7+
     2         t27-2*rhoa)-6.491760000000001d-3*gammaaa*t24*(1.111111111
     3         111111d-1*(t28-2.5d-1*t25*t29)*t7+t27)-6.491760000000001d
     4         -3*gammaaa*t24*(-7.777777777777777d-1*t23*t7+1.1111111111
     5         11111d-1*rhoa*t20-5.333333333333332d+0*rhoa)-1.1046240015
     6         73804d+0*t24*t5)+3.937253280921478d-2*t14*(1.735837716758
     7         835d+0*t24*t33+4.760624000000001d-2*gammaaa*t24*t32+2.380
     8         3120000000005d-2*gammaaa*t24*t31)/t33+3.125d-2*t14*t30*(-
     9         3.9971608514092083d-2*t24*t33-1.0962418720000001d-3*gamma
     :         aa*t24*t32-5.481209360000001d-4*gammaaa*t24*t31)+3.125d-2
     ;         *t15*t30*(-5.507339664989395d-2*t24*t33-1.51041616d-3*gam
     <         maaa*t24*t32-7.5520808d-4*gammaaa*t24*t31)-4.540977653965
     =         4695d-3*t1*t15-4.918d-2*t14)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.491760000000001d-3*t14*t
     1         16*t24*t32*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)-6.491760000000001d-3*t14*t
     1         16*t24*t31*t6*wght
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob+rhoa
            t2 = 1/t1
            t3 = 1/t1**3.333333333333333d-1
            t4 = 3.49d-1*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/t1**3.6666666666666664d+0
            t7 = rhoa**2.6666666666666666d+0
            t8 = rhob**2.6666666666666666d+0
            t9 = t8+t7
            t10 = 2.533d-1*t3
            t11 = exp(-t10)
            t12 = t1**2
            t13 = 3.49d-1*t3*t5
            t14 = t13+t10
            t15 = 4.7d+1-7.0d+0*t14
            t16 = 1.111111111111111d-1*rhoa*rhob*t15-1.333333333333333d+
     1         0*t12
            t17 = t13+t10-1.1d+1
            t18 = -3.0d+0*t14
            t19 = -1.0d+0*rhoa*t17*t2+t18+1.0d+0
            t20 = 1.111111111111111d-1*rhoa*rhob*t19-rhob**2
            t21 = -1.0d+0*rhob*t17*t2+t18+1.0d+0
            t22 = 1.111111111111111d-1*rhoa*rhob*t21-rhoa**2
            t23 = 1/t4**2
            t24 = -2.288509333333333d-2*rhoa*rhob*t23/t1**2.333333333333
     1         3334d+0
            t25 = 1/t12
            t26 = 1.9672d-1*rhoa*rhob*t25*t5
            t27 = -2.666666666666666d+0*t1
            t28 = 1/t1**1.3333333333333333d+0
            t29 = -1.1633333333333332d-1*t28*t5-8.443333333333334d-2*t28
     1         +4.060033333333332d-2*t23/t1**1.6666666666666669d+0
            t30 = -7.777777777777777d-1*rhoa*rhob*t29
            t31 = -3.0d+0*t29
            t32 = -1.0d+0*rhoa*t2*t29
            t33 = 1.0d+0*rhoa*t17*t25
            t34 = -1.0d+0*t17*t2
            t35 = -1.0d+0*rhob*t2*t29
            t36 = 1.0d+0*rhob*t17*t25
            t37 = 1/t1**5
            t38 = t23*t37*(-2.7536698324946973d-2*rhoa*rhob*t11*t9-7.552
     1         0808d-4*gammabb*t11*t22-7.5520808d-4*gammaaa*t11*t20-7.55
     2         20808d-4*gammaab*t11*t16)
            t39 = t37*t5*(-1.9985804257046041d-2*rhoa*rhob*t11*t9-5.4812
     1         09360000001d-4*gammabb*t11*t22-5.481209360000001d-4*gamma
     2         aa*t11*t20-5.481209360000001d-4*gammaab*t11*t16)
            t40 = t5*(8.679188583794175d-1*rhoa*rhob*t11*t9+2.3803120000
     1         000005d-2*gammabb*t11*t22+2.3803120000000005d-2*gammaaa*t
     2         11*t20+2.3803120000000005d-2*gammaab*t11*t16)/t1**4.66666
     3         6666666667d+0
            fnc(iq) = 1.0d+0*(-2.367051431943866d-1*rhoa*rhob*t11*t5*t6*
     1         t9-6.491760000000001d-3*gammabb*t11*t22*t5*t6-6.491760000
     2         000001d-3*gammaaa*t11*t20*t5*t6-6.491760000000001d-3*gamm
     3         aab*t11*t16*t5*t6-1.9672d-1*rhoa*rhob*t2*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(t5*t6*(-2.367051431943866d-1*rhob*t
     1         11*t9-6.312137151850309d-1*rhob*t11*t7-6.491760000000001d
     2         -3*gammabb*t11*(1.111111111111111d-1*rhoa*rhob*(t36+t35+t
     3         31)+1.111111111111111d-1*rhob*t21-2*rhoa)-6.4917600000000
     4         01d-3*gammaaa*t11*(1.111111111111111d-1*rhoa*rhob*(t34+t3
     5         3+t32+t31)+1.111111111111111d-1*rhob*t19)-6.4917600000000
     6         01d-3*gammaab*t11*(t30+t27+1.111111111111111d-1*rhob*t15)
     7         )-1.9672d-1*rhob*t2*t5+t40+t39+t38+t26+t24)*wght+Amat(iq,
     8         D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(t5*t6*(-2.367051431943866d-1*rhoa*t
     1         11*t9-6.312137151850309d-1*rhoa*t11*t8-6.491760000000001d
     2         -3*gammabb*t11*(1.111111111111111d-1*rhoa*rhob*(t36+t35+t
     3         34+t31)+1.111111111111111d-1*rhoa*t21)-6.491760000000001d
     4         -3*gammaaa*t11*(1.111111111111111d-1*rhoa*rhob*(t33+t32+t
     5         31)+1.111111111111111d-1*rhoa*t19-2*rhob)-6.4917600000000
     6         01d-3*gammaab*t11*(t30+t27+1.111111111111111d-1*rhoa*t15)
     7         )-1.9672d-1*rhoa*t2*t5+t40+t39+t38+t26+t24)*wght+Amat(iq,
     8         D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.491760000000001d-3*t11*t
     1         20*t5*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)-6.491760000000001d-3*t11*t
     1         16*t5*t6*wght
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.491760000000001d-3*t11*t
     1         22*t5*t6*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_lyp functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_1}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.349\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^{{{11}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 0.2533\,{\it t_3}\\\\ 
C>   {\it t_7} &=& e^ {- {\it t_6} }\\\\ 
C>   {\it t_8} &=& 0.349\,{\it t_3}\,{\it t_4}\\\\ 
C>   {\it t_9} &=& {\it t_8}+{\it t_6}-11.0\\\\ 
C>   {\it t_{10}} &=& {\it t_8}+{\it t_6}\\\\ 
C>   {\it t_{11}} &=& -3.0\,{\it t_{10}}\\\\ 
C>   f &=& 1.0\,\left(-0.006491760000000001\,{\it t_5}\,{\it t_4}
C>    \,\left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left({
C>    \it t_{11}}-1.0\,\rho_\beta\,{\it t_2}\,{\it t_9}+1.0\right)
C>    -\rho_\alpha^2\right)\,{\it t_7}\,\sigma_{\beta\beta}
C>    -0.006491760000000001\,{\it t_5}\,{\it t_4}\,
C>    \left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left(47.0
C>    -7.0\,{\it t_{10}}\right)-1.333333333333333\,{
C>    \it t_1}^2\right)\,{\it t_7}\,\sigma_{\alpha\beta}
C>    -0.006491760000000001\,{\it t_5}\,{\it t_4}\,
C>    \left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left({
C>    \it t_{11}}-1.0\,\rho_\alpha\,{\it t_2}\,{\it t_9}+1.0\right)
C>    -\rho_\beta^2\right)\,{\it t_7}\,\sigma_{\alpha\alpha}
C>    -0.2367051431943866\,\rho_\alpha\,\rho_\beta\,{\it t_5}
C>    \,\left(\rho_\beta^{{{8}\over{3}}}+\rho_\alpha^{{{8}
C>    \over{3}}}\right)\,{\it t_4}\,{\it t_7}-0.19672\,\rho_\alpha
C>    \,\rho_\beta\,{\it t_2}\,{\it t_4}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] C Lee, W Yang, RG Parr, Phys.Rev.B 37, 785 (1988)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.37.785 "> 
C> 10.1103/PhysRevB.37.785 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_lyp_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Cmat,Cmat2) 
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
      double precision rgamma(nq,*) !< [Input] The norm of the density 
                                    !< gradients 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma 
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
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**3.333333333333333d-1
            t2 = 1/(2.7700148356845083d-1*t1+1.0d+0)
            t3 = 2.010443432317725d-1*t1
            t4 = exp(-t3)
            t5 = rhoa**3.6666666666666664d+0
            t6 = 1/t5
            t7 = rhoa**2
            t8 = -5.333333333333332d+0*t7
            t9 = 2.7700148356845083d-1*t1*t2
            t10 = t9+t3
            t11 = -t7
            t12 = 7.937005259840998d-1
            t13 = 3.49d-1*t1*t12+1.0d+0
            t14 = 1/t13
            t15 = 1.9842513149602492d-1
            t16 = 1/t13**2
            t17 = 7.874506561842959d-2
            t18 = 2.533d-1*t1*t12
            t19 = 3.49d-1*t1*t12*t14
            t20 = t19+t18
            t21 = 4.7d+1-7.0d+0*t20
            t22 = 1/rhoa**1.6666666666666669d+0
            t23 = 3.968502629920499d-1
            t24 = 1/rhoa**1.3333333333333333d+0
            t25 = -1.1633333333333332d-1*t14*t23*t24-8.443333333333334d-
     1         2*t23*t24+1.2788303649853786d-2*t16*t22
            t26 = -7.777777777777777d-1*t25*t7+1.111111111111111d-1*rhoa
     1         *t21-5.333333333333332d+0*rhoa
            t27 = exp(-t18)
            t28 = t19+t18-1.1d+1
            t29 = -5.0d-1*t28-3.0d+0*t20+1.0d+0
            t30 = 1.111111111111111d-1*rhoa*t29
            t31 = -3.5d+0*t25
            t32 = 1/rhoa
            t33 = t31-2.5d-1*t28*t32
            t34 = 1.111111111111111d-1*t33*t7+t30
            t35 = 2.5d-1*t28*t32+t31
            t36 = 1.111111111111111d-1*t35*t7+t30-2*rhoa
            t37 = -1.104624001573804d+0*t27*t5-6.491760000000001d-3*gamm
     1         aaa*t27*t36-6.491760000000001d-3*gammaaa*t27*t34-6.491760
     2         000000001d-3*gammaaa*t26*t27
            t38 = 1/rhoa**5
            t39 = t8+1.111111111111111d-1*t21*t7
            t40 = 1.111111111111111d-1*t29*t7+t11
            t41 = rhoa**4.666666666666667d+0
            t42 = -5.507339664989395d-2*t27*t41-1.51041616d-3*gammaaa*t2
     1         7*t40-7.5520808d-4*gammaaa*t27*t39
            t43 = -3.9971608514092083d-2*t27*t41-1.0962418720000001d-3*g
     1         ammaaa*t27*t40-5.481209360000001d-4*gammaaa*t27*t39
            t44 = 3.937253280921477d-2
            t45 = 1/t41
            t46 = 1.735837716758835d+0*t27*t41+4.760624000000001d-2*gamm
     1         aaa*t27*t40+2.3803120000000005d-2*gammaaa*t27*t39
            t47 = 1/t13**3
            t48 = -5.324598382222221d-3*t17*t22*t47
            t49 = 7.568296089942449d-3*t16*t24
            t50 = -4.577018666666666d-2*t15*t16*t24
            t51 = -1.5555555555555553d+0*rhoa*t25
            t52 = rhoa**2.6666666666666666d+0
            t53 = rhoa**2.3333333333333334d+0
            t54 = 1/t53
            t55 = 1.551111111111111d-1*t14*t15*t54+1.1257777777777778d-1
     1         *t15*t54-1.2788303649853786d-2*t16/t52+1.1807930277777773
     2         d-3*t47/rhoa**3
            t56 = -7.777777777777777d-1*t55*t7
            t57 = -5.481209360000001d-4*gammaaa*t23*t24*t26*t27
            t58 = -3.5d+0*t55
            t59 = 1/t7
            t60 = -5.481209360000001d-4*gammaaa*t23*t24*t27*t34
            t61 = -5.481209360000001d-4*gammaaa*t23*t24*t27*t36
            t62 = -5.329547801878944d-2*t23*t27*t53
            t63 = -1.9985804257046041d-2*t12*t27*t53
            t64 = -3.6666666666666664d+0*t14*t37*t44*t45
            t65 = rhoa**3.3333333333333337d+0
            t66 = 3.125d-2*t14*t38*(-1.6874680727699207d-3*t12*t27*t65-9
     1         .326708653288152d-2*t27*t5-9.255935539253335d-5*gammaaa*t
     2         23*t24*t27*t40-4.6279677696266674d-5*gammaaa*t23*t24*t27*
     3         t39-5.481209360000001d-4*gammaaa*t27*t36-5.48120936000000
     4         1d-4*gammaaa*t27*t34-5.481209360000001d-4*gammaaa*t26*t27
     5         )
            t67 = t14*t44*t45*(7.328128227583548d-2*t12*t27*t65+4.050288
     1         005770614d+0*t27*t5+4.0195535306666674d-3*gammaaa*t23*t24
     2         *t27*t40+2.0097767653333337d-3*gammaaa*t23*t24*t27*t39+2.
     3         3803120000000005d-2*gammaaa*t27*t36+2.3803120000000005d-2
     4         *gammaaa*t27*t34+2.3803120000000005d-2*gammaaa*t26*t27)
            t68 = 1.2401570718501566d-2
            t69 = 1/rhoa**6.333333333333333d+0
            t70 = 2.3266666666666663d-1*t42*t47*t68*t69
            t71 = 1.1633333333333332d-1*t16*t43*t68*t69
            t72 = 1/rhoa**6
            t73 = -7.8125d-2*t14*t43*t72
            t74 = -2.3333333333333334d+0*t14*t44*t46/rhoa**5.66666666666
     1         6667d+0
            t75 = 3.125d-2*t16*t38*(-2.3250152285696893d-3*t12*t27*t65-1
     1         .2850459218308585d-1*t27*t5-1.2752947110933333d-4*gammaaa
     2         *t23*t24*t27*t40-6.376473555466666d-5*gammaaa*t23*t24*t27
     3         *t39+1.1633333333333332d-1*t37-7.5520808d-4*gammaaa*t27*t
     4         36-7.5520808d-4*gammaaa*t27*t34-7.5520808d-4*gammaaa*t26*
     5         t27)
            t76 = 1.5625d-2*t16*(1.1633333333333332d-1*t46-5*t42)*t72
            t77 = -2.36002525d-5*t16*t27*t38*t40
            t78 = -1.7128779250000004d-5*t14*t27*t38*t40
            t79 = 2.3803120000000005d-2*t14*t27*t40*t44*t45
            fnc(iq) = 1.0d+0*(-1.0223881343581931d-3*gammaaa*t2*t4*t6*(1
     1         .111111111111111d-1*t7*(-5.0d-1*(t9+t3-1.1d+1)-3.0d+0*t10
     2         +1.0d+0)+t11)-5.111940671790965d-4*gammaaa*t2*t4*t6*(t8+1
     3         .111111111111111d-1*(4.7d+1-7.0d+0*t10)*t7)-3.72787240661
     4         23486d-2*rhoa*t2*t4-9.836d-2*rhoa*t2)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(t14*t17*t37*t6+t14*t44*t45*t46+3.12
     1         5d-2*t14*t38*t43+3.125d-2*t16*t38*t42-2.288509333333333d-
     2         2*t1*t15*t16-4.918d-2*t14)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.491760000000001d-3*t14*t
     1         17*t27*t40*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)-6.491760000000001d-3*t14*t
     1         17*t27*t39*t6*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t76+t75+t74+t73+t71+t70+t14*t17
     1         *t6*(-6.491760000000001d-3*gammaaa*t27*(1.111111111111111
     2         d-1*(2.5d-1*t28*t59+t58-5.0d-1*t25*t32)*t7+2.222222222222
     3         222d-1*rhoa*t33)-6.491760000000001d-3*gammaaa*t27*(1.1111
     4         11111111111d-1*(-2.5d-1*t28*t59+t58+5.0d-1*t25*t32)*t7+2.
     5         222222222222222d-1*rhoa*t35-2)+t63+t62+t61+t60+t57-6.4917
     6         60000000001d-3*gammaaa*t27*(t56+t51-2.666666666666666d+0)
     7         -2.3144502890117796d+0*t27*t52)+t67+t66+t64+t50+t49+t48+4
     8         .918d-2*t14*t32)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t76+t75+t74+t73+t71+t70+t14*t17
     1         *t6*(-1.298352d-2*gammaaa*t27*(-3.8888888888888884d-1*t55
     2         *t7+1.111111111111111d-1*rhoa*t35+1.111111111111111d-1*rh
     3         oa*t33+1.111111111111111d-1*t29)+t63+t62+t61+t60+t57-6.49
     4         1760000000001d-3*gammaaa*t27*(t56+t51+1.111111111111111d-
     5         1*t21-2.666666666666666d+0)-1.735837716758835d+0*t27*t52)
     6         +t67+t66+t64+t50+t49+t48-4.918d-2*t14*t32)*wght+Amat2(iq,
     7         D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(t79+t78+t77-6.491760000000001d
     1         -3*t14*t17*t27*t34*t6)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(-6.491760000000001d-3*t14*t17*
     1         t26*t27*t6+2.3803120000000005d-2*t14*t27*t39*t44*t45-2.36
     2         002525d-5*t16*t27*t38*t39-1.7128779250000004d-5*t14*t27*t
     3         38*t39)*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = 1.0d+0*(t79+t78+t77-6.491760000000001d
     1         -3*t14*t17*t27*t36*t6)*wght+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob+rhoa
            t2 = 1/t1
            t3 = 1/t1**3.333333333333333d-1
            t4 = 3.49d-1*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/t1**3.6666666666666664d+0
            t7 = rhoa**2.6666666666666666d+0
            t8 = rhob**2.6666666666666666d+0
            t9 = t8+t7
            t10 = 2.533d-1*t3
            t11 = exp(-t10)
            t12 = t1**2
            t13 = 3.49d-1*t3*t5
            t14 = t13+t10
            t15 = 4.7d+1-7.0d+0*t14
            t16 = 1.111111111111111d-1*rhoa*rhob*t15-1.333333333333333d+
     1         0*t12
            t17 = t13+t10-1.1d+1
            t18 = -3.0d+0*t14
            t19 = -1.0d+0*rhoa*t17*t2+t18+1.0d+0
            t20 = 1.111111111111111d-1*rhoa*rhob*t19-rhob**2
            t21 = -1.0d+0*rhob*t17*t2+t18+1.0d+0
            t22 = 1.111111111111111d-1*rhoa*rhob*t21-rhoa**2
            t23 = 1/t1**2.3333333333333334d+0
            t24 = 1/t4**2
            t25 = -2.288509333333333d-2*rhoa*rhob*t23*t24
            t26 = 1/t12
            t27 = 1.9672d-1*rhoa*rhob*t26*t5
            t28 = -2.666666666666666d+0*t1
            t29 = 1/t1**1.3333333333333333d+0
            t30 = -1.1633333333333332d-1*t29*t5-8.443333333333334d-2*t29
     1         +4.060033333333332d-2*t24/t1**1.6666666666666669d+0
            t31 = -7.777777777777777d-1*rhoa*rhob*t30
            t32 = t31+t28+1.111111111111111d-1*rhob*t15
            t33 = -3.0d+0*t30
            t34 = -1.0d+0*rhoa*t2*t30
            t35 = 1.0d+0*rhoa*t17*t26
            t36 = -1.0d+0*t17*t2
            t37 = t36+t35+t34+t33
            t38 = 1.111111111111111d-1*rhoa*rhob*t37+1.111111111111111d-
     1         1*rhob*t19
            t39 = -1.0d+0*rhob*t2*t30
            t40 = 1.0d+0*rhob*t17*t26
            t41 = t40+t39+t33
            t42 = 1.111111111111111d-1*rhoa*rhob*t41+1.111111111111111d-
     1         1*rhob*t21-2*rhoa
            t43 = -2.367051431943866d-1*rhob*t11*t9-6.312137151850309d-1
     1         *rhob*t11*t7-6.491760000000001d-3*gammabb*t11*t42-6.49176
     2         0000000001d-3*gammaaa*t11*t38-6.491760000000001d-3*gammaa
     3         b*t11*t32
            t44 = 1/t1**5
            t45 = -2.7536698324946973d-2*rhoa*rhob*t11*t9-7.5520808d-4*g
     1         ammabb*t11*t22-7.5520808d-4*gammaaa*t11*t20-7.5520808d-4*
     2         gammaab*t11*t16
            t46 = t24*t44*t45
            t47 = -1.9985804257046041d-2*rhoa*rhob*t11*t9-5.481209360000
     1         001d-4*gammabb*t11*t22-5.481209360000001d-4*gammaaa*t11*t
     2         20-5.481209360000001d-4*gammaab*t11*t16
            t48 = t44*t47*t5
            t49 = 1/t1**4.666666666666667d+0
            t50 = 8.679188583794175d-1*rhoa*rhob*t11*t9+2.38031200000000
     1         05d-2*gammabb*t11*t22+2.3803120000000005d-2*gammaaa*t11*t
     2         20+2.3803120000000005d-2*gammaab*t11*t16
            t51 = t49*t5*t50
            t52 = t31+t28+1.111111111111111d-1*rhoa*t15
            t53 = t35+t34+t33
            t54 = 1.111111111111111d-1*rhoa*rhob*t53+1.111111111111111d-
     1         1*rhoa*t19-2*rhob
            t55 = t40+t39+t36+t33
            t56 = 1.111111111111111d-1*rhoa*rhob*t55+1.111111111111111d-
     1         1*rhoa*t21
            t57 = -2.367051431943866d-1*rhoa*t11*t9-6.312137151850309d-1
     1         *rhoa*t11*t8-6.491760000000001d-3*gammabb*t11*t56-6.49176
     2         0000000001d-3*gammaaa*t11*t54-6.491760000000001d-3*gammaa
     3         b*t11*t52
            t58 = 1/t4**3
            t59 = -5.324598382222221d-3*rhoa*rhob*t58*t6
            t60 = 7.628364444444444d-2*rhoa*rhob*t24/t1**3.3333333333333
     1         337d+0
            t61 = 1/t1**3
            t62 = -3.9344d-1*rhoa*rhob*t5*t61
            t63 = -3.6666666666666664d+0*t43*t49*t5
            t64 = -1.9985804257046041d-2*rhob*t11*t29*t9
            t65 = -5.329547801878944d-2*rhob*t11*t29*t7
            t66 = 9.446344222222218d-3*t58*t61+1.551111111111111d-1*t23*
     1         t5-8.120066666666665d-2*t24/t1**2.6666666666666666d+0+1.1
     2         257777777777778d-1*t23
            t67 = -7.777777777777777d-1*rhoa*rhob*t66
            t68 = -3.0d+0*t66
            t69 = -1.0d+0*rhob*t2*t66
            t70 = 2.0d+0*rhob*t26*t30
            t71 = -2.0d+0*rhob*t17*t61
            t72 = -1.0d+0*rhoa*t2*t66
            t73 = 2.0d+0*rhoa*t26*t30
            t74 = -2.0d+0*t2*t30
            t75 = -2.0d+0*rhoa*t17*t61
            t76 = 2.0d+0*t17*t26
            t77 = -5.481209360000001d-4*gammaab*t11*t29*t32
            t78 = -5.481209360000001d-4*gammaaa*t11*t29*t38
            t79 = -5.481209360000001d-4*gammabb*t11*t29*t42
            t80 = 1/t1**6.333333333333333d+0
            t81 = 2.3266666666666663d-1*t45*t58*t80
            t82 = 1.1633333333333332d-1*t24*t47*t80
            t83 = 1/t1**6
            t84 = -5*t47*t5*t83
            t85 = -4.666666666666667d+0*t5*t50/t1**5.666666666666667d+0
            t86 = -1.6874680727699207d-3*rhoa*rhob*t11*t29*t9
            t87 = -4.6279677696266674d-5*gammaab*t11*t16*t29
            t88 = -4.6279677696266674d-5*gammaaa*t11*t20*t29
            t89 = -4.6279677696266674d-5*gammabb*t11*t22*t29
            t90 = 7.328128227583548d-2*rhoa*rhob*t11*t29*t9
            t91 = 2.0097767653333337d-3*gammaab*t11*t16*t29
            t92 = 2.0097767653333337d-3*gammaaa*t11*t20*t29
            t93 = 2.0097767653333337d-3*gammabb*t11*t22*t29
            t94 = -2.3250152285696893d-3*rhoa*rhob*t11*t29*t9
            t95 = -6.376473555466666d-5*gammaab*t11*t16*t29
            t96 = -6.376473555466666d-5*gammaaa*t11*t20*t29
            t97 = -6.376473555466666d-5*gammabb*t11*t22*t29
            t98 = 1.1633333333333332d-1*t43
            t99 = t24*(1.1633333333333332d-1*t50-5*t45)*t83
            t100 = -1.0d+0*t2*t30
            t101 = 1.0d+0*t17*t26
            t102 = t44*t5*(-1.9985804257046041d-2*rhoa*t11*t9+t89+t88+t8
     1         7+t86-5.329547801878943d-2*rhoa*t11*t8-5.481209360000001d
     2         -4*gammabb*t11*t56-5.481209360000001d-4*gammaaa*t11*t54-5
     3         .481209360000001d-4*gammaab*t11*t52)
            t103 = t49*t5*(t93+t92+t91+t90+8.679188583794175d-1*rhoa*t11
     1         *t9+2.3144502890117796d+0*rhoa*t11*t8+2.3803120000000005d
     2         -2*gammabb*t11*t56+2.3803120000000005d-2*gammaaa*t11*t54+
     3         2.3803120000000005d-2*gammaab*t11*t52)
            t104 = -7.343119553319191d-2*rhoa*t11*t8
            t105 = -2.7536698324946973d-2*rhoa*t11*t9
            t106 = -7.5520808d-4*gammaab*t11*t52
            t107 = -7.5520808d-4*gammaaa*t11*t54
            t108 = -7.5520808d-4*gammabb*t11*t56
            t109 = -7.5520808d-4*t11*t20*t24*t44
            t110 = -5.481209360000001d-4*t11*t20*t44*t5
            t111 = 2.3803120000000005d-2*t11*t20*t49*t5
            t112 = -7.5520808d-4*t11*t16*t24*t44
            t113 = -5.481209360000001d-4*t11*t16*t44*t5
            t114 = 2.3803120000000005d-2*t11*t16*t49*t5
            t115 = -7.5520808d-4*t11*t22*t24*t44
            t116 = -5.481209360000001d-4*t11*t22*t44*t5
            t117 = 2.3803120000000005d-2*t11*t22*t49*t5
            fnc(iq) = 1.0d+0*(-2.367051431943866d-1*rhoa*rhob*t11*t5*t6*
     1         t9-6.491760000000001d-3*gammabb*t11*t22*t5*t6-6.491760000
     2         000001d-3*gammaaa*t11*t20*t5*t6-6.491760000000001d-3*gamm
     3         aab*t11*t16*t5*t6-1.9672d-1*rhoa*rhob*t2*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(t43*t5*t6+t51-1.9672d-1*rhob*t2*t5+
     1         t48+t46+t27+t25)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(t5*t57*t6+t51-1.9672d-1*rhoa*t2*t5+
     1         t48+t46+t27+t25)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.491760000000001d-3*t11*t
     1         20*t5*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)-6.491760000000001d-3*t11*t
     1         16*t5*t6*wght
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.491760000000001d-3*t11*t
     1         22*t5*t6*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t99+t24*t44*(t98+t97+t96+t95+t9
     1         4-2.7536698324946973d-2*rhob*t11*t9-7.343119553319191d-2*
     2         rhob*t11*t7-7.5520808d-4*gammabb*t11*t42-7.5520808d-4*gam
     3         maaa*t11*t38-7.5520808d-4*gammaab*t11*t32)+t49*t5*(t93+t9
     4         2+t91+t90+8.679188583794175d-1*rhob*t11*t9+2.314450289011
     5         7796d+0*rhob*t11*t7+2.3803120000000005d-2*gammabb*t11*t42
     6         +2.3803120000000005d-2*gammaaa*t11*t38+2.3803120000000005
     7         d-2*gammaab*t11*t32)+t44*t5*(-1.9985804257046041d-2*rhob*
     8         t11*t9+t89+t88+t87+t86-5.329547801878943d-2*rhob*t11*t7-5
     9         .481209360000001d-4*gammabb*t11*t42-5.481209360000001d-4*
     :         gammaaa*t11*t38-5.481209360000001d-4*gammaab*t11*t32)+t85
     ;         +t84+t82+t81+t5*t6*(t79+t78+t77-6.491760000000001d-3*gamm
     <         aaa*t11*(1.111111111111111d-1*rhoa*rhob*(t76+t75+t74+t73+
     =         t72+t68)+2.222222222222222d-1*rhob*t37)-6.491760000000001
     >         d-3*gammabb*t11*(1.111111111111111d-1*rhoa*rhob*(t71+t70+
     ?         t69+t68)+2.222222222222222d-1*rhob*t41-2)-6.4917600000000
     @         01d-3*gammaab*t11*(t67-1.5555555555555553d+0*rhob*t30-2.6
     1         66666666666666d+0)+t65+t64-2.3144502890117796d+0*rhoa**1.
     2         6666666666666669d+0*rhob*t11)+t63+t62+t60+t59+3.9344d-1*r
     3         hob*t26*t5-4.577018666666666d-2*rhob*t23*t24)*wght+Amat2(
     4         iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t99+t24*t44*(t98+t97+t96+t95+t9
     1         4+t108+t107+t106+t105+t104)+t5*t6*(-2.367051431943866d-1*
     2         t11*t9-6.312137151850309d-1*t11*t8+t79+t78+t77-6.49176000
     3         0000001d-3*gammaaa*t11*(1.111111111111111d-1*rhoa*rhob*(t
     4         75+t73+t72+t68+t101+t100)+1.111111111111111d-1*rhob*t53+1
     5         .111111111111111d-1*rhoa*t37+1.111111111111111d-1*t19)-6.
     6         491760000000001d-3*gammabb*t11*(1.111111111111111d-1*rhoa
     7         *rhob*(t71+t70+t69+t68+t101+t100)+1.111111111111111d-1*rh
     8         ob*t55+1.111111111111111d-1*rhoa*t41+1.111111111111111d-1
     9         *t21)-6.312137151850309d-1*t11*t7-6.491760000000001d-3*ga
     :         mmaab*t11*(t67-7.777777777777777d-1*rhob*t30-7.7777777777
     ;         77777d-1*rhoa*t30+1.111111111111111d-1*t15-2.666666666666
     <         666d+0)+t65+t64)+t85+t84+t82+t81+t63+t62+t60+t59+(1.9672d
     =         -1*rhob+1.9672d-1*rhoa)*t26*t5-1.9672d-1*t2*t5+(-2.288509
     >         333333333d-2*rhob-2.288509333333333d-2*rhoa)*t23*t24+t103
     ?         +t102)*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(t99+t24*t44*(t97+t96+t95+t94+1.
     1         1633333333333332d-1*t57+t108+t107+t106+t105+t104)+t5*t6*(
     2         -1.9985804257046041d-2*rhoa*t11*t29*t9-5.329547801878944d
     3         -2*rhoa*t11*t29*t8-6.491760000000001d-3*gammabb*t11*(1.11
     4         1111111111111d-1*rhoa*rhob*(t76+t74+t71+t70+t69+t68)+2.22
     5         2222222222222d-1*rhoa*t55)-6.491760000000001d-3*gammaaa*t
     6         11*(1.111111111111111d-1*rhoa*rhob*(t75+t73+t72+t68)+2.22
     7         2222222222222d-1*rhoa*t53-2)-6.491760000000001d-3*gammaab
     8         *t11*(t67-1.5555555555555553d+0*rhoa*t30-2.66666666666666
     9         6d+0)-5.481209360000001d-4*gammabb*t11*t29*t56-5.48120936
     :         0000001d-4*gammaaa*t11*t29*t54-5.481209360000001d-4*gamma
     ;         ab*t11*t29*t52-2.3144502890117796d+0*rhoa*rhob**1.6666666
     <         666666669d+0*t11)+t85+t84+t82+t81+t62+t60+t59-3.666666666
     =         6666664d+0*t49*t5*t57+3.9344d-1*rhoa*t26*t5-4.57701866666
     >         6666d-2*rhoa*t23*t24+t103+t102)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-6.491760000000001d-3*t11*t38*
     1         t5*t6+t111+t110+t109)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(-6.491760000000001d-3*t11*t32*
     1         t5*t6+t114+t113+t112)*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = 1.0d+0*(-6.491760000000001d-3*t11*t42*
     1         t5*t6+t117+t116+t115)*wght+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = 1.0d+0*(-6.491760000000001d-3*t11*t5*t
     1         54*t6+t111+t110+t109)*wght+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = 1.0d+0*(-6.491760000000001d-3*t11*t5*t
     1         52*t6+t114+t113+t112)*wght+Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(-6.491760000000001d-3*t11*t5*t
     1         56*t6+t117+t116+t115)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_lyp functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_1}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.349\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^{{{11}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 0.2533\,{\it t_3}\\\\ 
C>   {\it t_7} &=& e^ {- {\it t_6} }\\\\ 
C>   {\it t_8} &=& 0.349\,{\it t_3}\,{\it t_4}\\\\ 
C>   {\it t_9} &=& {\it t_8}+{\it t_6}-11.0\\\\ 
C>   {\it t_{10}} &=& {\it t_8}+{\it t_6}\\\\ 
C>   {\it t_{11}} &=& -3.0\,{\it t_{10}}\\\\ 
C>   f &=& 1.0\,\left(-0.006491760000000001\,{\it t_5}\,{\it t_4}
C>    \,\left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left({
C>    \it t_{11}}-1.0\,\rho_\beta\,{\it t_2}\,{\it t_9}+1.0\right)
C>    -\rho_\alpha^2\right)\,{\it t_7}\,\sigma_{\beta\beta}
C>    -0.006491760000000001\,{\it t_5}\,{\it t_4}\,
C>    \left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left(47.0
C>    -7.0\,{\it t_{10}}\right)-1.333333333333333\,{
C>    \it t_1}^2\right)\,{\it t_7}\,\sigma_{\alpha\beta}
C>    -0.006491760000000001\,{\it t_5}\,{\it t_4}\,
C>    \left(0.1111111111111111\,\rho_\alpha\,\rho_\beta\,\left({
C>    \it t_{11}}-1.0\,\rho_\alpha\,{\it t_2}\,{\it t_9}+1.0\right)
C>    -\rho_\beta^2\right)\,{\it t_7}\,\sigma_{\alpha\alpha}
C>    -0.2367051431943866\,\rho_\alpha\,\rho_\beta\,{\it t_5}
C>    \,\left(\rho_\beta^{{{8}\over{3}}}+\rho_\alpha^{{{8}
C>    \over{3}}}\right)\,{\it t_4}\,{\it t_7}-0.19672\,\rho_\alpha
C>    \,\rho_\beta\,{\it t_2}\,{\it t_4}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] C Lee, W Yang, RG Parr, Phys.Rev.B 37, 785 (1988)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.37.785 "> 
C> 10.1103/PhysRevB.37.785 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_lyp_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3) 
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
      double precision rgamma(nq,*) !< [Input] The norm of the density 
                                    !< gradients 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma 
                                    !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho 
      double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma 
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
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**3.333333333333333d-1
            t2 = 1/(2.7700148356845083d-1*t1+1.0d+0)
            t3 = 2.010443432317725d-1*t1
            t4 = exp(-t3)
            t5 = rhoa**3.6666666666666664d+0
            t6 = 1/t5
            t7 = rhoa**2
            t8 = -5.333333333333332d+0*t7
            t9 = 2.7700148356845083d-1*t1*t2
            t10 = t9+t3
            t11 = -t7
            t12 = 7.937005259840998d-1
            t13 = 3.49d-1*t1*t12+1.0d+0
            t14 = 1/t13
            t15 = 1.9842513149602492d-1
            t16 = 1/t13**2
            t17 = 7.874506561842959d-2
            t18 = 2.533d-1*t1*t12
            t19 = 3.49d-1*t1*t12*t14
            t20 = t19+t18
            t21 = 4.7d+1-7.0d+0*t20
            t22 = 3.149802624737183d-1
            t23 = rhoa**1.6666666666666669d+0
            t24 = 1/t23
            t25 = 3.968502629920499d-1
            t26 = rhoa**1.3333333333333333d+0
            t27 = 1/t26
            t28 = -1.1633333333333332d-1*t14*t25*t27-8.443333333333334d-
     1         2*t25*t27+4.060033333333332d-2*t16*t22*t24
            t29 = -7.777777777777777d-1*t28*t7+1.111111111111111d-1*rhoa
     1         *t21-5.333333333333332d+0*rhoa
            t30 = exp(-t18)
            t31 = t19+t18-1.1d+1
            t32 = -5.0d-1*t31-3.0d+0*t20+1.0d+0
            t33 = 1.111111111111111d-1*rhoa*t32
            t34 = -3.5d+0*t28
            t35 = 1/rhoa
            t36 = t34-2.5d-1*t31*t35
            t37 = 1.111111111111111d-1*t36*t7+t33
            t38 = 2.5d-1*t31*t35+t34
            t39 = 1.111111111111111d-1*t38*t7+t33-2*rhoa
            t40 = -1.104624001573804d+0*t30*t5-6.491760000000001d-3*gamm
     1         aaa*t30*t39-6.491760000000001d-3*gammaaa*t30*t37-6.491760
     2         000000001d-3*gammaaa*t29*t30
            t41 = 1/rhoa**5
            t42 = t8+1.111111111111111d-1*t21*t7
            t43 = 1.111111111111111d-1*t32*t7+t11
            t44 = rhoa**4.666666666666667d+0
            t45 = -5.507339664989395d-2*t30*t44-1.51041616d-3*gammaaa*t3
     1         0*t43-7.5520808d-4*gammaaa*t30*t42
            t46 = -3.9971608514092083d-2*t30*t44-1.0962418720000001d-3*g
     1         ammaaa*t30*t43-5.481209360000001d-4*gammaaa*t30*t42
            t47 = 3.937253280921477d-2
            t48 = 1/t44
            t49 = 1.735837716758835d+0*t30*t44+4.760624000000001d-2*gamm
     1         aaa*t30*t43+2.3803120000000005d-2*gammaaa*t30*t42
            t50 = 1/t13**3
            t51 = -5.324598382222221d-3*t17*t24*t50
            t52 = 9.921256574801247d-2
            t53 = 7.628364444444444d-2*t16*t27*t52
            t54 = -4.577018666666666d-2*t15*t16*t27
            t55 = -1.5555555555555553d+0*rhoa*t28
            t56 = 1/rhoa**3
            t57 = 1.5749013123685918d-1
            t58 = rhoa**2.6666666666666666d+0
            t59 = 1/t58
            t60 = rhoa**2.3333333333333334d+0
            t61 = 1/t60
            t62 = 1.551111111111111d-1*t14*t15*t61+1.1257777777777778d-1
     1         *t15*t61-8.120066666666665d-2*t16*t57*t59+1.1807930277777
     2         773d-3*t50*t56
            t63 = -7.777777777777777d-1*t62*t7
            t64 = t63+t55-2.666666666666666d+0
            t65 = -5.481209360000001d-4*gammaaa*t25*t27*t29*t30
            t66 = -3.5d+0*t62
            t67 = 1/t7
            t68 = 2.5d-1*t31*t67+t66-5.0d-1*t28*t35
            t69 = 1.111111111111111d-1*t68*t7+2.222222222222222d-1*rhoa*
     1         t36
            t70 = -2.5d-1*t31*t67+t66+5.0d-1*t28*t35
            t71 = 1.111111111111111d-1*t7*t70+2.222222222222222d-1*rhoa*
     1         t38-2
            t72 = -5.481209360000001d-4*gammaaa*t25*t27*t30*t37
            t73 = -5.481209360000001d-4*gammaaa*t25*t27*t30*t39
            t74 = -5.329547801878944d-2*t25*t30*t60
            t75 = -1.9985804257046041d-2*t12*t30*t60
            t76 = t75+t74+t73+t72-6.491760000000001d-3*gammaaa*t30*t71-6
     1         .491760000000001d-3*gammaaa*t30*t69+t65-6.491760000000001
     2         d-3*gammaaa*t30*t64-2.3144502890117796d+0*t30*t58
            t77 = -3.6666666666666664d+0*t14*t40*t47*t48
            t78 = rhoa**3.3333333333333337d+0
            t79 = -1.6874680727699207d-3*t12*t30*t78-9.326708653288152d-
     1         2*t30*t5-9.255935539253335d-5*gammaaa*t25*t27*t30*t43-4.6
     2         279677696266674d-5*gammaaa*t25*t27*t30*t42-5.481209360000
     3         001d-4*gammaaa*t30*t39-5.481209360000001d-4*gammaaa*t30*t
     4         37-5.481209360000001d-4*gammaaa*t29*t30
            t80 = 3.125d-2*t14*t41*t79
            t81 = 7.328128227583548d-2*t12*t30*t78+4.050288005770614d+0*
     1         t30*t5+4.0195535306666674d-3*gammaaa*t25*t27*t30*t43+2.00
     2         97767653333337d-3*gammaaa*t25*t27*t30*t42+2.3803120000000
     3         005d-2*gammaaa*t30*t39+2.3803120000000005d-2*gammaaa*t30*
     4         t37+2.3803120000000005d-2*gammaaa*t29*t30
            t82 = t14*t47*t48*t81
            t83 = 1.2401570718501566d-2
            t84 = 1/rhoa**6.333333333333333d+0
            t85 = 2.3266666666666663d-1*t45*t50*t83*t84
            t86 = 1.1633333333333332d-1*t16*t46*t83*t84
            t87 = 1/rhoa**6
            t88 = -7.8125d-2*t14*t46*t87
            t89 = 1/rhoa**5.666666666666667d+0
            t90 = -2.3333333333333334d+0*t14*t47*t49*t89
            t91 = -7.5520808d-4*gammaaa*t29*t30
            t92 = -7.5520808d-4*gammaaa*t30*t37
            t93 = -7.5520808d-4*gammaaa*t30*t39
            t94 = -6.376473555466666d-5*gammaaa*t25*t27*t30*t42
            t95 = -1.2752947110933333d-4*gammaaa*t25*t27*t30*t43
            t96 = -2.3250152285696893d-3*t12*t30*t78
            t97 = -1.2850459218308585d-1*t30*t5
            t98 = t97+t96+t95+t94+t93+t92+t91+1.1633333333333332d-1*t40
            t99 = 3.125d-2*t16*t41*t98
            t100 = 1.1633333333333332d-1*t49-5*t45
            t101 = 1.5625d-2*t100*t16*t87
            t102 = t63+t55+1.111111111111111d-1*t21-2.666666666666666d+0
            t103 = -3.8888888888888884d-1*t62*t7+1.111111111111111d-1*rh
     1         oa*t38+1.111111111111111d-1*rhoa*t36+1.111111111111111d-1
     2         *t32
            t104 = t75+t74+t73+t72+t65-1.735837716758835d+0*t30*t58-1.29
     1         8352d-2*gammaaa*t103*t30-6.491760000000001d-3*gammaaa*t10
     2         2*t30
            t105 = -2.36002525d-5*t16*t30*t41*t43
            t106 = -1.7128779250000004d-5*t14*t30*t41*t43
            t107 = 2.3803120000000005d-2*t14*t30*t43*t47*t48
            t108 = 1/t13**4
            t109 = -5.80714011061111d-5*t108*t56
            t110 = 3.7272188675555545d-2*t47*t50*t59
            t111 = -1.5973795146666664d-2*t17*t50*t59
            t112 = 4.960628287400625d-2
            t113 = -3.0004900148148145d-1*t112*t16*t61
            t114 = 2.288509333333333d-1*t16*t52*t61
            t115 = -4.4999148607197886d-3*rhoa*t30*t57
            t116 = -1.6874680727699207d-3*rhoa*t22*t30
            t117 = 7.106063735838591d-2*t15*t26*t30
            t118 = -2.333333333333333d+0*rhoa*t62
            t119 = 1/t78
            t120 = 2.3457970370370365d-1*t16*t17*t6-3.619259259259259d-1
     1         *t119*t14*t52-2.626814814814815d-1*t119*t52-2.95198256944
     2         44435d-3*t50/rhoa**4+3.296774133555554d-3*t108*t112/rhoa*
     3         *4.333333333333333d+0
            t121 = -7.777777777777777d-1*t120*t7
            t122 = -4.6279677696266674d-5*gammaaa*t29*t30*t57*t59
            t123 = 7.308279146666667d-4*gammaaa*t15*t29*t30*t61
            t124 = -3.5d+0*t120
            t125 = -4.6279677696266674d-5*gammaaa*t30*t37*t57*t59
            t126 = 7.308279146666667d-4*gammaaa*t15*t30*t37*t61
            t127 = -4.6279677696266674d-5*gammaaa*t30*t39*t57*t59
            t128 = 7.308279146666667d-4*gammaaa*t15*t30*t39*t61
            t129 = -1.4247855427754028d-4*t22*t30*t7
            t130 = -9.255935539253335d-5*gammaaa*t25*t27*t29*t30
            t131 = -9.255935539253335d-5*gammaaa*t25*t27*t30*t37
            t132 = -9.255935539253335d-5*gammaaa*t25*t27*t30*t39
            t133 = -3.907547453488116d-6*gammaaa*t30*t42*t57*t59
            t134 = 6.170623692835556d-5*gammaaa*t15*t30*t42*t61
            t135 = -7.815094906976232d-6*gammaaa*t30*t43*t57*t59
            t136 = 1.2341247385671113d-4*gammaaa*t15*t30*t43*t61
            t137 = -6.749872291079682d-3*t25*t30*t60
            t138 = -3.3749361455398413d-3*t12*t30*t60
            t139 = 6.187382933489709d-3*t22*t30*t7
            t140 = 4.0195535306666674d-3*gammaaa*t25*t27*t29*t30
            t141 = 4.0195535306666674d-3*gammaaa*t25*t27*t30*t37
            t142 = 4.0195535306666674d-3*gammaaa*t25*t27*t30*t39
            t143 = 1.6969215155297782d-4*gammaaa*t30*t42*t57*t59
            t144 = -2.679702353777778d-3*gammaaa*t15*t30*t42*t61
            t145 = 3.3938430310595563d-4*gammaaa*t30*t43*t57*t59
            t146 = -5.359404707555556d-3*gammaaa*t15*t30*t43*t61
            t147 = 2.9312512910334193d-1*t25*t30*t60
            t148 = 1.4656256455167097d-1*t12*t30*t60
            t149 = 8.555555555555555d+0*t14*t40*t47*t89
            t150 = 2.3266666666666663d-1*t16*t79*t83*t84
            t151 = -1.5625d-1*t14*t79*t87
            t152 = -2.3333333333333334d+0*t14*t17*t81*t89
            t153 = 4.921566601151848d-3
            t154 = 1/rhoa**7.666666666666667d+0
            t155 = 8.120066666666664d-2*t108*t153*t154*t45
            t156 = 2.706688888888888d-2*t153*t154*t46*t50
            t157 = 6.200785359250782d-3
            t158 = 1/rhoa**7.333333333333333d+0
            t159 = -1.3184444444444443d+0*t157*t158*t16*t46
            t160 = 1/rhoa**7
            t161 = 2.34375d-1*t14*t160*t46
            t162 = 1.9686266404607394d-2
            t163 = 1.322222222222222d+1*t14*t162*t49/rhoa**6.66666666666
     1         6667d+0
            t164 = -1.9630878579890076d-4*t22*t30*t7
            t165 = -1.2752947110933333d-4*gammaaa*t25*t27*t29*t30
            t166 = -1.2752947110933333d-4*gammaaa*t25*t27*t30*t37
            t167 = -1.2752947110933333d-4*gammaaa*t25*t27*t30*t39
            t168 = -5.383869171999022d-6*gammaaa*t30*t42*t57*t59
            t169 = 8.501964740622221d-5*gammaaa*t15*t30*t42*t61
            t170 = -1.0767738343998043d-5*gammaaa*t30*t43*t57*t59
            t171 = 1.7003929481244442d-4*gammaaa*t15*t30*t43*t61
            t172 = -4.6500304571393786d-3*t12*t30*t60
            t173 = t97+t96+t95+t94+t93+t92+t91
            t174 = 1.5625d-2*t16*t87*(-5*t98+2.3266666666666663d-1*t81-4
     1         .2655555555555547d-1*t40-5*t173)
            t175 = t50*t83*t84*(2.3266666666666663d-1*t98+2.326666666666
     1         6663d-1*t173)
            t176 = 7.8125d-3*t16*t160*(-5.428888888888889d-1*t49-6*t100)
            t177 = t157*t158*(2.3266666666666663d-1*t100-1.4735555555555
     1         552d+0*t45)*t50
            t178 = -7.777777777777777d-1*rhoa*t62
            t179 = -1.757117466133333d-4*t30*t43*t50*t83*t84
            t180 = -6.376473555466666d-5*t16*t30*t43*t83*t84
            t181 = 1.0226776083333333d-4*t16*t30*t43*t87
            t182 = 4.282194812500001d-5*t14*t30*t43*t87
            t183 = -1.110812266666667d-1*t14*t162*t30*t43*t89
            t184 = -5.481209360000001d-4*t25*t27*t30*t37
            t185 = -6.376473555466666d-5*t25*t27*t30*t43
            t186 = -4.6279677696266674d-5*t25*t27*t30*t43
            t187 = 2.0097767653333337d-3*t25*t27*t30*t43
            t188 = -1.757117466133333d-4*t30*t42*t50*t83*t84
            t189 = -6.376473555466666d-5*t16*t30*t42*t83*t84
            t190 = 1.0226776083333333d-4*t16*t30*t42*t87
            t191 = 4.282194812500001d-5*t14*t30*t42*t87
            t192 = -1.110812266666667d-1*t14*t162*t30*t42*t89
            t193 = -5.481209360000001d-4*t25*t27*t29*t30
            t194 = 3.125d-2*t16*t41*(-6.376473555466666d-5*t25*t27*t30*t
     1         42-1.51041616d-3*t29*t30)
            t195 = 3.125d-2*t14*t41*(-4.6279677696266674d-5*t25*t27*t30*
     1         t42-5.481209360000001d-4*t29*t30)
            t196 = t14*(2.0097767653333337d-3*t25*t27*t30*t42+4.76062400
     1         0000001d-2*t29*t30)*t47*t48
            t197 = 3.125d-2*t14*(t186-5.481209360000001d-4*t30*t39)*t41
            fnc(iq) = 1.0d+0*(-1.0223881343581931d-3*gammaaa*t2*t4*t6*(1
     1         .111111111111111d-1*t7*(-5.0d-1*(t9+t3-1.1d+1)-3.0d+0*t10
     2         +1.0d+0)+t11)-5.111940671790965d-4*gammaaa*t2*t4*t6*(t8+1
     3         .111111111111111d-1*(4.7d+1-7.0d+0*t10)*t7)-3.72787240661
     4         23486d-2*rhoa*t2*t4-9.836d-2*rhoa*t2)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(t14*t17*t40*t6+t14*t47*t48*t49+3.12
     1         5d-2*t14*t41*t46+3.125d-2*t16*t41*t45-2.288509333333333d-
     2         2*t1*t15*t16-4.918d-2*t14)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.491760000000001d-3*t14*t
     1         17*t30*t43*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)-6.491760000000001d-3*t14*t
     1         17*t30*t42*t6*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t99+t90+t88+t86+t85+t82+t80+t77
     1         +t14*t17*t6*t76+t54+t53+t51+4.918d-2*t14*t35+t101)*wght+A
     2         mat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t99+t90+t88+t86+t85+t82+t80+t77
     1         +t104*t14*t17*t6+t54+t53+t51-4.918d-2*t14*t35+t101)*wght+
     2         Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-6.491760000000001d-3*t14*t17*
     1         t30*t37*t6+t107+t106+t105)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(-6.491760000000001d-3*t14*t17*
     1         t29*t30*t6+2.3803120000000005d-2*t14*t30*t42*t47*t48-2.36
     2         002525d-5*t16*t30*t41*t42-1.7128779250000004d-5*t14*t30*t
     3         41*t42)*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = 1.0d+0*(-6.491760000000001d-3*t14*t17*
     1         t30*t39*t6+t107+t106+t105)*wght+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-3.6666666666666664d+0*t14*t
     1         17*t48*t76+3.125d-2*t16*t41*(2.3266666666666663d-1*t76-7.
     2         5520808d-4*gammaaa*t30*t71-7.5520808d-4*gammaaa*t30*t69-7
     3         .5520808d-4*gammaaa*t30*t64-9.300060914278757d-3*t25*t30*
     4         t60-2.69247716955037d-1*t30*t58+t172+t171+t170+t169+t168+
     5         t167+t166+t165+t164)+t14*t17*t6*(-1.0962418720000001d-3*g
     6         ammaaa*t25*t27*t30*t71-6.491760000000001d-3*gammaaa*t30*(
     7         3.333333333333333d-1*rhoa*t70+1.111111111111111d-1*(-7.5d
     8         -1*t28*t67+7.5d-1*t35*t62+3.75d-1*t31*t56+t124)*t7)-6.491
     9         760000000001d-3*gammaaa*t30*(1.111111111111111d-1*(7.5d-1
     :         *t28*t67-7.5d-1*t35*t62-3.75d-1*t31*t56+t124)*t7+3.333333
     ;         333333333d-1*rhoa*t68)-1.0962418720000001d-3*gammaaa*t25*
     <         t27*t30*t69-1.0962418720000001d-3*gammaaa*t25*t27*t30*t64
     =         -3.6418576646172784d-1*t25*t26*t30-3.857417148352966d+0*t
     >         23*t30-6.491760000000001d-3*gammaaa*(t121+t118)*t30+t128+
     ?         t127+t126+t125+t123+t122+t117+t116+t115)+t14*t47*t48*(2.3
     @         803120000000005d-2*gammaaa*t30*t71+2.3803120000000005d-2*
     1         gammaaa*t30*t69+2.3803120000000005d-2*gammaaa*t30*t64+8.4
     2         86317726376527d+0*t30*t58+t148+t147+t146+t145+t144+t143+t
     3         142+t141+t140+t139)+3.125d-2*t14*t41*(-5.481209360000001d
     4         -4*gammaaa*t30*t71-5.481209360000001d-4*gammaaa*t30*t69-5
     5         .481209360000001d-4*gammaaa*t30*t64-1.9541675273556125d-1
     6         *t30*t58+t138+t137+t136+t135+t134+t133+t132+t131+t130+t12
     7         9)-7.377d-2*t14*t67+t177+t176+t175+t174+t163+t161+t159+t1
     8         56+t155+t152+t151+t150+t149+t114+t113+t111+t110+t109)*wgh
     9         t+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 1.0d+0*(3.125d-2*t16*t41*(1.16333333
     1         33333332d-1*t76-9.300060914278756d-3*t25*t30*t60-2.019357
     2         8771627776d-1*t30*t58-1.51041616d-3*gammaaa*t103*t30-7.55
     3         20808d-4*gammaaa*t102*t30+t172+t171+t170+t169+t168+t167+t
     4         166+t165+t164+1.1633333333333332d-1*t104)+3.3333333333333
     5         33d-1*t14*t47*t48*(-11*t76-11*t104)+t14*t17*t6*(-5.481209
     6         360000001d-4*gammaaa*t25*t27*t30*t71-6.491760000000001d-3
     7         *gammaaa*t30*(1.111111111111111d-1*rhoa*t70+1.11111111111
     8         1111d-1*(-2.5d-1*t28*t67+2.5d-1*t35*t62+1.25d-1*t31*t56+t
     9         124)*t7+2.222222222222222d-1*t38+t178)-6.491760000000001d
     :         -3*gammaaa*t30*(1.111111111111111d-1*(2.5d-1*t28*t67-2.5d
     ;         -1*t35*t62-1.25d-1*t31*t56+t124)*t7+1.111111111111111d-1*
     <         rhoa*t68+2.222222222222222d-1*t36+t178)-5.481209360000001
     =         d-4*gammaaa*t25*t27*t30*t69-5.481209360000001d-4*gammaaa*
     >         t25*t27*t30*t64-6.491760000000001d-3*gammaaa*(-1.55555555
     ?         55555553d+0*t28+t121+t118)*t30-1.0962418720000001d-3*gamm
     @         aaa*t103*t25*t27*t30-5.481209360000001d-4*gammaaa*t102*t2
     1         5*t27*t30-2.7535996976374544d-1*t25*t26*t30-1.99858042570
     2         46041d-2*t12*t26*t30-2.3144502890117796d+0*t23*t30+t128+t
     3         127+t126+t125+t123+t122+t117+t116+t115)+2.459d-2*t14*t67-
     4         4.577018666666666d-2*t15*t16*t61+t14*t47*t48*(6.364738294
     5         782395d+0*t30*t58+4.760624000000001d-2*gammaaa*t103*t30+2
     6         .3803120000000005d-2*gammaaa*t102*t30+t148+t147+t146+t145
     7         +t144+t143+t142+t141+t140+t139)+3.125d-2*t14*t41*(-1.4656
     8         256455167094d-1*t30*t58-1.0962418720000001d-3*gammaaa*t10
     9         3*t30-5.481209360000001d-4*gammaaa*t102*t30+t138+t137+t13
     :         6+t135+t134+t133+t132+t131+t130+t129)+t177+t176+t175+t174
     ;         +t163+t161+t159+t156+t155+t152+t151+t150+t149+t114+t113+t
     <         111+t110+t109)*wght+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(t14*t17*t6*(t184-6.49176000
     1         0000001d-3*t30*t69)+t14*(4.760624000000001d-2*t30*t37+t18
     2         7)*t47*t48+3.125d-2*t14*(t186-5.481209360000001d-4*t30*t3
     3         7)*t41+3.125d-2*t16*(t185-1.51041616d-3*t30*t37)*t41+t183
     4         +t182+t181+t180+t179)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = 1.0d+0*(t14*t17*t6*(t193-6.49176000
     1         0000001d-3*t30*t64)+t196+t195+t194+t192+t191+t190+t189+t1
     2         88)*wght+Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = 1.0d+0*(t14*t17*t6*(-6.491760000000
     1         001d-3*t30*t71-5.481209360000001d-4*t25*t27*t30*t39)+t14*
     2         (4.760624000000001d-2*t30*t39+t187)*t47*t48+3.125d-2*t16*
     3         (t185-1.51041616d-3*t30*t39)*t41+t197+t183+t182+t181+t180
     4         +t179)*wght+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = 1.0d+0*(t14*t17*(t184-6.49176000000
     1         0001d-3*t103*t30)*t6+t14*(2.3803120000000005d-2*t30*t39+2
     2         .3803120000000005d-2*t30*t37+t187)*t47*t48+3.125d-2*t16*(
     3         -7.5520808d-4*t30*t39-7.5520808d-4*t30*t37+t185)*t41+t197
     4         +t183+t182+t181+t180+t179)*wght+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = 1.0d+0*(t14*t17*(t193-6.49176000000
     1         0001d-3*t102*t30)*t6+t196+t195+t194+t192+t191+t190+t189+t
     2         188)*wght+Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob+rhoa
            t2 = 1/t1
            t3 = 1/t1**3.333333333333333d-1
            t4 = 3.49d-1*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/t1**3.6666666666666664d+0
            t7 = rhoa**2.6666666666666666d+0
            t8 = rhob**2.6666666666666666d+0
            t9 = t8+t7
            t10 = 2.533d-1*t3
            t11 = exp(-t10)
            t12 = t1**2
            t13 = 3.49d-1*t3*t5
            t14 = t13+t10
            t15 = 4.7d+1-7.0d+0*t14
            t16 = 1.111111111111111d-1*rhoa*rhob*t15-1.333333333333333d+
     1         0*t12
            t17 = t13+t10-1.1d+1
            t18 = -3.0d+0*t14
            t19 = -1.0d+0*rhoa*t17*t2+t18+1.0d+0
            t20 = 1.111111111111111d-1*rhoa*rhob*t19-rhob**2
            t21 = -1.0d+0*rhob*t17*t2+t18+1.0d+0
            t22 = 1.111111111111111d-1*rhoa*rhob*t21-rhoa**2
            t23 = 1/t1**2.3333333333333334d+0
            t24 = 1/t4**2
            t25 = -2.288509333333333d-2*rhoa*rhob*t23*t24
            t26 = 1/t12
            t27 = 1.9672d-1*rhoa*rhob*t26*t5
            t28 = -2.666666666666666d+0*t1
            t29 = 1/t1**1.3333333333333333d+0
            t30 = -1.1633333333333332d-1*t29*t5-8.443333333333334d-2*t29
     1         +4.060033333333332d-2*t24/t1**1.6666666666666669d+0
            t31 = -7.777777777777777d-1*rhoa*rhob*t30
            t32 = t31+t28+1.111111111111111d-1*rhob*t15
            t33 = -3.0d+0*t30
            t34 = -1.0d+0*rhoa*t2*t30
            t35 = 1.0d+0*rhoa*t17*t26
            t36 = -1.0d+0*t17*t2
            t37 = t36+t35+t34+t33
            t38 = 1.111111111111111d-1*rhoa*rhob*t37+1.111111111111111d-
     1         1*rhob*t19
            t39 = -1.0d+0*rhob*t2*t30
            t40 = 1.0d+0*rhob*t17*t26
            t41 = t40+t39+t33
            t42 = 1.111111111111111d-1*rhoa*rhob*t41+1.111111111111111d-
     1         1*rhob*t21-2*rhoa
            t43 = -2.367051431943866d-1*rhob*t11*t9-6.312137151850309d-1
     1         *rhob*t11*t7-6.491760000000001d-3*gammabb*t11*t42-6.49176
     2         0000000001d-3*gammaaa*t11*t38-6.491760000000001d-3*gammaa
     3         b*t11*t32
            t44 = 1/t1**5
            t45 = -2.7536698324946973d-2*rhoa*rhob*t11*t9-7.5520808d-4*g
     1         ammabb*t11*t22-7.5520808d-4*gammaaa*t11*t20-7.5520808d-4*
     2         gammaab*t11*t16
            t46 = t24*t44*t45
            t47 = -1.9985804257046041d-2*rhoa*rhob*t11*t9-5.481209360000
     1         001d-4*gammabb*t11*t22-5.481209360000001d-4*gammaaa*t11*t
     2         20-5.481209360000001d-4*gammaab*t11*t16
            t48 = t44*t47*t5
            t49 = 1/t1**4.666666666666667d+0
            t50 = 8.679188583794175d-1*rhoa*rhob*t11*t9+2.38031200000000
     1         05d-2*gammabb*t11*t22+2.3803120000000005d-2*gammaaa*t11*t
     2         20+2.3803120000000005d-2*gammaab*t11*t16
            t51 = t49*t5*t50
            t52 = t31+t28+1.111111111111111d-1*rhoa*t15
            t53 = t35+t34+t33
            t54 = 1.111111111111111d-1*rhoa*rhob*t53+1.111111111111111d-
     1         1*rhoa*t19-2*rhob
            t55 = t40+t39+t36+t33
            t56 = 1.111111111111111d-1*rhoa*rhob*t55+1.111111111111111d-
     1         1*rhoa*t21
            t57 = -2.367051431943866d-1*rhoa*t11*t9-6.312137151850309d-1
     1         *rhoa*t11*t8-6.491760000000001d-3*gammabb*t11*t56-6.49176
     2         0000000001d-3*gammaaa*t11*t54-6.491760000000001d-3*gammaa
     3         b*t11*t52
            t58 = 1/t4**3
            t59 = -5.324598382222221d-3*rhoa*rhob*t58*t6
            t60 = 1/t1**3.3333333333333337d+0
            t61 = 7.628364444444444d-2*rhoa*rhob*t24*t60
            t62 = 1/t1**3
            t63 = -3.9344d-1*rhoa*rhob*t5*t62
            t64 = -3.6666666666666664d+0*t43*t49*t5
            t65 = rhoa**1.6666666666666669d+0
            t66 = -1.9985804257046041d-2*rhob*t11*t29*t9
            t67 = -5.329547801878944d-2*rhob*t11*t29*t7
            t68 = 1/t1**2.6666666666666666d+0
            t69 = -8.120066666666665d-2*t24*t68+9.446344222222218d-3*t58
     1         *t62+1.551111111111111d-1*t23*t5+1.1257777777777778d-1*t2
     2         3
            t70 = -7.777777777777777d-1*rhoa*rhob*t69
            t71 = t70-1.5555555555555553d+0*rhob*t30-2.666666666666666d+
     1         0
            t72 = -3.0d+0*t69
            t73 = -1.0d+0*rhob*t2*t69
            t74 = 2.0d+0*rhob*t26*t30
            t75 = -2.0d+0*rhob*t17*t62
            t76 = t75+t74+t73+t72
            t77 = 1.111111111111111d-1*rhoa*rhob*t76+2.222222222222222d-
     1         1*rhob*t41-2
            t78 = -1.0d+0*rhoa*t2*t69
            t79 = 2.0d+0*rhoa*t26*t30
            t80 = -2.0d+0*t2*t30
            t81 = -2.0d+0*rhoa*t17*t62
            t82 = 2.0d+0*t17*t26
            t83 = t82+t81+t80+t79+t78+t72
            t84 = 1.111111111111111d-1*rhoa*rhob*t83+2.222222222222222d-
     1         1*rhob*t37
            t85 = -5.481209360000001d-4*gammaab*t11*t29*t32
            t86 = -5.481209360000001d-4*gammaaa*t11*t29*t38
            t87 = -5.481209360000001d-4*gammabb*t11*t29*t42
            t88 = t87+t86+t85-6.491760000000001d-3*gammaaa*t11*t84-6.491
     1         760000000001d-3*gammabb*t11*t77-6.491760000000001d-3*gamm
     2         aab*t11*t71+t67+t66-2.3144502890117796d+0*rhob*t11*t65
            t89 = 1/t1**6.333333333333333d+0
            t90 = 2.3266666666666663d-1*t45*t58*t89
            t91 = 1.1633333333333332d-1*t24*t47*t89
            t92 = 1/t1**6
            t93 = -5*t47*t5*t92
            t94 = 1/t1**5.666666666666667d+0
            t95 = -4.666666666666667d+0*t5*t50*t94
            t96 = -1.6874680727699207d-3*rhoa*rhob*t11*t29*t9
            t97 = -4.6279677696266674d-5*gammaab*t11*t16*t29
            t98 = -4.6279677696266674d-5*gammaaa*t11*t20*t29
            t99 = -4.6279677696266674d-5*gammabb*t11*t22*t29
            t100 = t99+t98+t97+t96-1.9985804257046041d-2*rhob*t11*t9-5.3
     1         29547801878943d-2*rhob*t11*t7-5.481209360000001d-4*gammab
     2         b*t11*t42-5.481209360000001d-4*gammaaa*t11*t38-5.48120936
     3         0000001d-4*gammaab*t11*t32
            t101 = 7.328128227583548d-2*rhoa*rhob*t11*t29*t9
            t102 = 2.0097767653333337d-3*gammaab*t11*t16*t29
            t103 = 2.0097767653333337d-3*gammaaa*t11*t20*t29
            t104 = 2.0097767653333337d-3*gammabb*t11*t22*t29
            t105 = 8.679188583794175d-1*rhob*t11*t9+2.3144502890117796d+
     1         0*rhob*t11*t7+2.3803120000000005d-2*gammabb*t11*t42+2.380
     2         3120000000005d-2*gammaaa*t11*t38+2.3803120000000005d-2*ga
     3         mmaab*t11*t32+t104+t103+t102+t101
            t106 = -7.343119553319191d-2*rhob*t11*t7
            t107 = -2.7536698324946973d-2*rhob*t11*t9
            t108 = -2.3250152285696893d-3*rhoa*rhob*t11*t29*t9
            t109 = -7.5520808d-4*gammaab*t11*t32
            t110 = -6.376473555466666d-5*gammaab*t11*t16*t29
            t111 = -7.5520808d-4*gammaaa*t11*t38
            t112 = -6.376473555466666d-5*gammaaa*t11*t20*t29
            t113 = -7.5520808d-4*gammabb*t11*t42
            t114 = -6.376473555466666d-5*gammabb*t11*t22*t29
            t115 = 1.1633333333333332d-1*t43
            t116 = t115+t114+t113+t112+t111+t110+t109+t108+t107+t106
            t117 = 1.1633333333333332d-1*t50-5*t45
            t118 = t117*t24*t92
            t119 = -2.288509333333333d-2*rhob-2.288509333333333d-2*rhoa
            t120 = 1.9672d-1*rhob+1.9672d-1*rhoa
            t121 = t70-7.777777777777777d-1*rhob*t30-7.777777777777777d-
     1         1*rhoa*t30+1.111111111111111d-1*t15-2.666666666666666d+0
            t122 = -1.0d+0*t2*t30
            t123 = 1.0d+0*t17*t26
            t124 = t81+t79+t78+t72+t123+t122
            t125 = 1.111111111111111d-1*rhob*t53+1.111111111111111d-1*rh
     1         oa*t37+1.111111111111111d-1*t19+1.111111111111111d-1*rhoa
     2         *rhob*t124
            t126 = t75+t74+t73+t72+t123+t122
            t127 = 1.111111111111111d-1*rhob*t55+1.111111111111111d-1*rh
     1         oa*t41+1.111111111111111d-1*t21+1.111111111111111d-1*rhoa
     2         *rhob*t126
            t128 = -2.367051431943866d-1*t11*t9+t87+t86+t85-6.3121371518
     1         50309d-1*t11*t8-6.312137151850309d-1*t11*t7+t67+t66-6.491
     2         760000000001d-3*gammabb*t11*t127-6.491760000000001d-3*gam
     3         maaa*t11*t125-6.491760000000001d-3*gammaab*t11*t121
            t129 = t99+t98+t97+t96-1.9985804257046041d-2*rhoa*t11*t9-5.3
     1         29547801878943d-2*rhoa*t11*t8-5.481209360000001d-4*gammab
     2         b*t11*t56-5.481209360000001d-4*gammaaa*t11*t54-5.48120936
     3         0000001d-4*gammaab*t11*t52
            t130 = t129*t44*t5
            t131 = 8.679188583794175d-1*rhoa*t11*t9+2.3144502890117796d+
     1         0*rhoa*t11*t8+2.3803120000000005d-2*gammabb*t11*t56+2.380
     2         3120000000005d-2*gammaaa*t11*t54+2.3803120000000005d-2*ga
     3         mmaab*t11*t52+t104+t103+t102+t101
            t132 = t131*t49*t5
            t133 = -7.343119553319191d-2*rhoa*t11*t8
            t134 = -2.7536698324946973d-2*rhoa*t11*t9
            t135 = -7.5520808d-4*gammaab*t11*t52
            t136 = -7.5520808d-4*gammaaa*t11*t54
            t137 = -7.5520808d-4*gammabb*t11*t56
            t138 = t137+t136+t135+t134+t133+t115+t114+t112+t110+t108
            t139 = rhob**1.6666666666666669d+0
            t140 = t70-1.5555555555555553d+0*rhoa*t30-2.666666666666666d
     1         +0
            t141 = t81+t79+t78+t72
            t142 = 2.222222222222222d-1*rhoa*t53+1.111111111111111d-1*rh
     1         oa*rhob*t141-2
            t143 = t82+t80+t75+t74+t73+t72
            t144 = 2.222222222222222d-1*rhoa*t55+1.111111111111111d-1*rh
     1         oa*rhob*t143
            t145 = -1.9985804257046041d-2*rhoa*t11*t29*t9-5.329547801878
     1         944d-2*rhoa*t11*t29*t8-5.481209360000001d-4*gammabb*t11*t
     2         29*t56-5.481209360000001d-4*gammaaa*t11*t29*t54-5.4812093
     3         60000001d-4*gammaab*t11*t29*t52-6.491760000000001d-3*gamm
     4         abb*t11*t144-6.491760000000001d-3*gammaaa*t11*t142-6.4917
     5         60000000001d-3*gammaab*t11*t140-2.3144502890117796d+0*rho
     6         a*t11*t139
            t146 = 1.1633333333333332d-1*t57+t137+t136+t135+t134+t133+t1
     1         14+t112+t110+t108
            t147 = -7.5520808d-4*t11*t20*t24*t44
            t148 = -5.481209360000001d-4*t11*t20*t44*t5
            t149 = 2.3803120000000005d-2*t11*t20*t49*t5
            t150 = -7.5520808d-4*t11*t16*t24*t44
            t151 = -5.481209360000001d-4*t11*t16*t44*t5
            t152 = 2.3803120000000005d-2*t11*t16*t49*t5
            t153 = -7.5520808d-4*t11*t22*t24*t44
            t154 = -5.481209360000001d-4*t11*t22*t44*t5
            t155 = 2.3803120000000005d-2*t11*t22*t49*t5
            t156 = 1/t4**4
            t157 = -1.858284835395555d-3*rhoa*rhob*t156*t44
            t158 = 3.7272188675555545d-2*rhoa*rhob*t49*t58
            t159 = 1/t1**4.333333333333333d+0
            t160 = -3.0004900148148145d-1*rhoa*rhob*t159*t24
            t161 = 1/t1**4
            t162 = 1.18032d+0*rhoa*rhob*t161*t5
            t163 = 1.711111111111111d+1*t43*t5*t94
            t164 = -1.6874680727699207d-3*rhob*t11*t68*t9
            t165 = 2.6647739009394716d-2*rhob*t11*t23*t9
            t166 = -4.4999148607197886d-3*rhob*t11*t68*t7
            t167 = 7.106063735838591d-2*rhob*t11*t23*t7
            t168 = -3.619259259259259d-1*t5*t60-2.626814814814815d-1*t60
     1         +2.3457970370370365d-1*t24*t6-4.7231721111111097d-2*t161*
     2         t58+3.296774133555554d-3*t156*t159
            t169 = -7.777777777777777d-1*rhoa*rhob*t168
            t170 = -3.0d+0*t168
            t171 = -1.0d+0*rhob*t168*t2
            t172 = 3.0d+0*rhob*t26*t69
            t173 = -6.0d+0*rhob*t30*t62
            t174 = 6.0d+0*rhob*t161*t17
            t175 = -1.0d+0*rhoa*t168*t2
            t176 = 3.0d+0*rhoa*t26*t69
            t177 = -3.0d+0*t2*t69
            t178 = -6.0d+0*rhoa*t30*t62
            t179 = 6.0d+0*t26*t30
            t180 = 6.0d+0*rhoa*t161*t17
            t181 = -6.0d+0*t17*t62
            t182 = -4.6279677696266674d-5*gammaab*t11*t32*t68
            t183 = 7.308279146666667d-4*gammaab*t11*t23*t32
            t184 = -4.6279677696266674d-5*gammaaa*t11*t38*t68
            t185 = 7.308279146666667d-4*gammaaa*t11*t23*t38
            t186 = -4.6279677696266674d-5*gammabb*t11*t42*t68
            t187 = 7.308279146666667d-4*gammabb*t11*t23*t42
            t188 = 1/t1**7.666666666666667d+0
            t189 = 8.120066666666664d-2*t156*t188*t45
            t190 = 2.706688888888888d-2*t188*t47*t58
            t191 = 1/t1**7.333333333333333d+0
            t192 = -1.3184444444444443d+0*t191*t24*t47
            t193 = 1/t1**7
            t194 = 30*t193*t47*t5
            t195 = 2.644444444444444d+1*t5*t50/t1**6.666666666666667d+0
            t196 = 6.187382933489709d-3*rhoa*rhob*t11*t68*t9
            t197 = -9.770837636778064d-2*rhoa*rhob*t11*t23*t9
            t198 = 1.6969215155297782d-4*gammaab*t11*t16*t68
            t199 = -2.679702353777778d-3*gammaab*t11*t16*t23
            t200 = 1.6969215155297782d-4*gammaaa*t11*t20*t68
            t201 = -2.679702353777778d-3*gammaaa*t11*t20*t23
            t202 = 1.6969215155297782d-4*gammabb*t11*t22*t68
            t203 = -2.679702353777778d-3*gammabb*t11*t22*t23
            t204 = -1.4247855427754028d-4*rhoa*rhob*t11*t68*t9
            t205 = 2.249957430359894d-3*rhoa*rhob*t11*t23*t9
            t206 = -3.907547453488116d-6*gammaab*t11*t16*t68
            t207 = 6.170623692835556d-5*gammaab*t11*t16*t23
            t208 = -3.907547453488116d-6*gammaaa*t11*t20*t68
            t209 = 6.170623692835556d-5*gammaaa*t11*t20*t23
            t210 = -3.907547453488116d-6*gammabb*t11*t22*t68
            t211 = 6.170623692835556d-5*gammabb*t11*t22*t23
            t212 = -1.9630878579890076d-4*rhoa*rhob*t11*t68*t9
            t213 = 3.100020304759586d-3*rhoa*rhob*t11*t23*t9
            t214 = -5.383869171999022d-6*gammaab*t11*t16*t68
            t215 = 8.501964740622221d-5*gammaab*t11*t16*t23
            t216 = -5.383869171999022d-6*gammaaa*t11*t20*t68
            t217 = 8.501964740622221d-5*gammaaa*t11*t20*t23
            t218 = -5.383869171999022d-6*gammabb*t11*t22*t68
            t219 = 8.501964740622221d-5*gammabb*t11*t22*t23
            t220 = -4.2655555555555547d-1*t43
            t221 = t114+t113+t112+t111+t110+t109+t108+t107+t106
            t222 = -5*t116
            t223 = 2.3266666666666663d-1*t116
            t224 = t193*t24*(-5.428888888888889d-1*t50-6*t117)
            t225 = t191*(2.3266666666666663d-1*t117-1.4735555555555552d+
     1         0*t45)*t58
            t226 = -5.324598382222221d-3*rhoa
            t227 = 7.628364444444444d-2*rhoa
            t228 = -4.577018666666666d-2*t23*t24
            t229 = -3.9344d-1*rhoa
            t230 = 3.9344d-1*t26*t5
            t231 = -1.5555555555555553d+0*t30
            t232 = -1.0d+0*t2*t69
            t233 = 2.0d+0*t26*t30
            t234 = -2.0d+0*t17*t62
            t235 = -2.0d+0*t2*t69
            t236 = 4.0d+0*t26*t30
            t237 = -4.0d+0*t17*t62
            t238 = t137+t136+t135+t134+t133+t114+t112+t110+t108
            t239 = -5*t238
            t240 = 2.3266666666666663d-1*t238
            t241 = t49*t5*(1.4656256455167097d-1*rhoa*t11*t29*t9+3.90833
     1         50547112256d-1*rhoa*t11*t29*t8+4.0195535306666674d-3*gamm
     2         abb*t11*t29*t56+4.0195535306666674d-3*gammaaa*t11*t29*t54
     3         +4.0195535306666674d-3*gammaab*t11*t29*t52+t203+t202+t201
     4         +t200+t199+t198+t197+t196+2.3803120000000005d-2*gammabb*t
     5         11*t144+2.3803120000000005d-2*gammaaa*t11*t142+2.38031200
     6         00000005d-2*gammaab*t11*t140+8.486317726376527d+0*rhoa*t1
     7         1*t139)
            t242 = t44*t5*(-3.3749361455398413d-3*rhoa*t11*t29*t9-8.9998
     1         29721439576d-3*rhoa*t11*t29*t8-9.255935539253335d-5*gamma
     2         bb*t11*t29*t56-9.255935539253335d-5*gammaaa*t11*t29*t54-9
     3         .255935539253335d-5*gammaab*t11*t29*t52+t211+t210+t209+t2
     4         08+t207+t206+t205+t204-5.481209360000001d-4*gammabb*t11*t
     5         144-5.481209360000001d-4*gammaaa*t11*t142-5.4812093600000
     6         01d-4*gammaab*t11*t140-1.9541675273556125d-1*rhoa*t11*t13
     7         9)
            t243 = 2.3266666666666663d-1*t129*t24*t89
            t244 = -10*t129*t5*t92
            t245 = -9.333333333333333d+0*t131*t5*t94
            t246 = -2.69247716955037d-1*rhoa*t11*t139
            t247 = -4.6500304571393786d-3*rhoa*t11*t29*t9
            t248 = -1.2400081219038343d-2*rhoa*t11*t29*t8
            t249 = -7.5520808d-4*gammaab*t11*t140
            t250 = -7.5520808d-4*gammaaa*t11*t142
            t251 = -7.5520808d-4*gammabb*t11*t144
            t252 = -1.2752947110933333d-4*gammaab*t11*t29*t52
            t253 = -1.2752947110933333d-4*gammaaa*t11*t29*t54
            t254 = -1.2752947110933333d-4*gammabb*t11*t29*t56
            t255 = 2.3266666666666663d-1*t131
            t256 = -1.757117466133333d-4*t11*t20*t58*t89
            t257 = -6.376473555466666d-5*t11*t20*t24*t89
            t258 = 6.545136693333333d-3*t11*t20*t24*t92
            t259 = 2.7406046800000006d-3*t11*t20*t5*t92
            t260 = -1.110812266666667d-1*t11*t20*t5*t94
            t261 = -5.481209360000001d-4*t11*t29*t38
            t262 = -6.376473555466666d-5*t11*t20*t29
            t263 = -4.6279677696266674d-5*t11*t20*t29
            t264 = 2.0097767653333337d-3*t11*t20*t29
            t265 = -1.757117466133333d-4*t11*t16*t58*t89
            t266 = -6.376473555466666d-5*t11*t16*t24*t89
            t267 = 6.545136693333333d-3*t11*t16*t24*t92
            t268 = 2.7406046800000006d-3*t11*t16*t5*t92
            t269 = -1.110812266666667d-1*t11*t16*t5*t94
            t270 = -5.481209360000001d-4*t11*t29*t32
            t271 = -6.376473555466666d-5*t11*t16*t29
            t272 = -4.6279677696266674d-5*t11*t16*t29
            t273 = 2.0097767653333337d-3*t11*t16*t29
            t274 = -1.757117466133333d-4*t11*t22*t58*t89
            t275 = -6.376473555466666d-5*t11*t22*t24*t89
            t276 = 6.545136693333333d-3*t11*t22*t24*t92
            t277 = 2.7406046800000006d-3*t11*t22*t5*t92
            t278 = -1.110812266666667d-1*t11*t22*t5*t94
            t279 = -5.481209360000001d-4*t11*t29*t42
            t280 = -6.376473555466666d-5*t11*t22*t29
            t281 = -4.6279677696266674d-5*t11*t22*t29
            t282 = 2.0097767653333337d-3*t11*t22*t29
            t283 = t44*t5*(t263-5.481209360000001d-4*t11*t54)
            t284 = t44*t5*(t272-5.481209360000001d-4*t11*t52)
            t285 = t44*t5*(t281-5.481209360000001d-4*t11*t56)
            fnc(iq) = 1.0d+0*(-2.367051431943866d-1*rhoa*rhob*t11*t5*t6*
     1         t9-6.491760000000001d-3*gammabb*t11*t22*t5*t6-6.491760000
     2         000001d-3*gammaaa*t11*t20*t5*t6-6.491760000000001d-3*gamm
     3         aab*t11*t16*t5*t6-1.9672d-1*rhoa*rhob*t2*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(t43*t5*t6+t51-1.9672d-1*rhob*t2*t5+
     1         t48+t46+t27+t25)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(t5*t57*t6+t51-1.9672d-1*rhoa*t2*t5+
     1         t48+t46+t27+t25)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.491760000000001d-3*t11*t
     1         20*t5*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)-6.491760000000001d-3*t11*t
     1         16*t5*t6*wght
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.491760000000001d-3*t11*t
     1         22*t5*t6*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(t95+t93+t91+t90+t5*t6*t88+t64+t
     1         63+t61+t59+t105*t49*t5+t100*t44*t5+3.9344d-1*rhob*t26*t5+
     2         t116*t24*t44-4.577018666666666d-2*rhob*t23*t24+t118)*wght
     3         +Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.0d+0*(t95+t93+t91+t90+t64+t63+t61+t12
     1         8*t5*t6+t59+t120*t26*t5-1.9672d-1*t2*t5+t138*t24*t44+t119
     2         *t23*t24+t132+t130+t118)*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(t95+t93+t91+t90+t63+t61+t145*t5
     1         *t6+t59-3.6666666666666664d+0*t49*t5*t57+3.9344d-1*rhoa*t
     2         26*t5+t146*t24*t44-4.577018666666666d-2*rhoa*t23*t24+t132
     3         +t130+t118)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-6.491760000000001d-3*t11*t38*
     1         t5*t6+t149+t148+t147)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 1.0d+0*(-6.491760000000001d-3*t11*t32*
     1         t5*t6+t152+t151+t150)*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = 1.0d+0*(-6.491760000000001d-3*t11*t42*
     1         t5*t6+t155+t154+t153)*wght+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = 1.0d+0*(-6.491760000000001d-3*t11*t5*t
     1         54*t6+t149+t148+t147)*wght+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = 1.0d+0*(-6.491760000000001d-3*t11*t5*t
     1         52*t6+t152+t151+t150)*wght+Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(-6.491760000000001d-3*t11*t5*t
     1         56*t6+t155+t154+t153)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-9.333333333333333d+0*t105*t
     1         5*t94-10*t100*t5*t92+(t222-5*t221+t220+2.3266666666666663
     2         d-1*t105)*t24*t92+t49*t5*(1.4656256455167097d-1*rhob*t11*
     3         t29*t9+2.3803120000000005d-2*gammaaa*t11*t84+2.3803120000
     4         000005d-2*gammabb*t11*t77+2.3803120000000005d-2*gammaab*t
     5         11*t71+3.9083350547112256d-1*rhob*t11*t29*t7+8.4863177263
     6         76527d+0*rhob*t11*t65+4.0195535306666674d-3*gammabb*t11*t
     7         29*t42+4.0195535306666674d-3*gammaaa*t11*t29*t38+4.019553
     8         5306666674d-3*gammaab*t11*t29*t32+t203+t202+t201+t200+t19
     9         9+t198+t197+t196)+t44*t5*(-3.3749361455398413d-3*rhob*t11
     :         *t29*t9-5.481209360000001d-4*gammaaa*t11*t84-5.4812093600
     ;         00001d-4*gammabb*t11*t77-5.481209360000001d-4*gammaab*t11
     <         *t71-8.999829721439576d-3*rhob*t11*t29*t7-1.9541675273556
     =         125d-1*rhob*t11*t65-9.255935539253335d-5*gammabb*t11*t29*
     >         t42-9.255935539253335d-5*gammaaa*t11*t29*t38-9.2559355392
     ?         53335d-5*gammaab*t11*t29*t32+t211+t210+t209+t208+t207+t20
     @         6+t205+t204)+t24*t44*(-4.6500304571393786d-3*rhob*t11*t29
     1         *t9+2.3266666666666663d-1*t88-7.5520808d-4*gammaaa*t11*t8
     2         4-7.5520808d-4*gammabb*t11*t77-7.5520808d-4*gammaab*t11*t
     3         71-1.2400081219038343d-2*rhob*t11*t29*t7-2.69247716955037
     4         d-1*rhob*t11*t65-1.2752947110933333d-4*gammabb*t11*t29*t4
     5         2-1.2752947110933333d-4*gammaaa*t11*t29*t38-1.27529471109
     6         33333d-4*gammaab*t11*t29*t32+t219+t218+t217+t216+t215+t21
     7         4+t213+t212)+(t223+2.3266666666666663d-1*t221)*t58*t89+2.
     8         3266666666666663d-1*t100*t24*t89-7.333333333333333d+0*t49
     9         *t5*t88+t5*t6*(-1.0962418720000001d-3*gammaaa*t11*t29*t84
     :         -6.491760000000001d-3*gammaaa*t11*(3.333333333333333d-1*r
     ;         hob*t83+1.111111111111111d-1*rhoa*rhob*(t181+t180+t179+t1
     <         78+t177+t176+t175+t170))-1.0962418720000001d-3*gammabb*t1
     =         1*t29*t77-6.491760000000001d-3*gammabb*t11*(3.33333333333
     >         3333d-1*rhob*t76+1.111111111111111d-1*rhoa*rhob*(t174+t17
     ?         3+t172+t171+t170))-1.0962418720000001d-3*gammaab*t11*t29*
     @         t71-6.491760000000001d-3*gammaab*t11*(t169-2.333333333333
     1         333d+0*rhob*t69)-3.9083350547112256d-1*rhob*t11*t29*t65+t
     2         187+t186+t185+t184+t183+t182+t167+t166+t165+t164-3.857417
     3         148352966d+0*rhoa**6.666666666666666d-1*rhob*t11)-1.18032
     4         d+0*rhob*t5*t62+2.288509333333333d-1*rhob*t24*t60-1.59737
     5         95146666664d-2*rhob*t58*t6+t225+t224+t195+t194+t192+t190+
     6         t189+t163+t162+t160+t158+t157)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 1.0d+0*(3.333333333333333d-1*(-14*t1
     1         31-14*t105)*t5*t94+(-5*t129-5*t100)*t5*t92+(t239+t222+t22
     2         0+1.1633333333333332d-1*t131+1.1633333333333332d-1*t105)*
     3         t24*t92+t49*t5*(7.328128227583548d-2*rhob*t11*t29*t9+7.32
     4         8128227583548d-2*rhoa*t11*t29*t9+8.679188583794175d-1*t11
     5         *t9+1.9541675273556128d-1*rhoa*t11*t29*t8+2.3144502890117
     6         796d+0*t11*t8+1.9541675273556128d-1*rhob*t11*t29*t7+2.314
     7         4502890117796d+0*t11*t7+2.0097767653333337d-3*gammabb*t11
     8         *t29*t56+2.0097767653333337d-3*gammaaa*t11*t29*t54+2.0097
     9         767653333337d-3*gammaab*t11*t29*t52+2.0097767653333337d-3
     :         *gammabb*t11*t29*t42+2.0097767653333337d-3*gammaaa*t11*t2
     ;         9*t38+2.0097767653333337d-3*gammaab*t11*t29*t32+t203+t202
     <         +t201+t200+t199+t198+t197+t196+2.3803120000000005d-2*gamm
     =         abb*t11*t127+2.3803120000000005d-2*gammaaa*t11*t125+2.380
     >         3120000000005d-2*gammaab*t11*t121)+t44*t5*(-1.68746807276
     ?         99207d-3*rhob*t11*t29*t9-1.6874680727699207d-3*rhoa*t11*t
     @         29*t9-1.9985804257046041d-2*t11*t9-4.499914860719788d-3*r
     1         hoa*t11*t29*t8-5.329547801878943d-2*t11*t8-4.499914860719
     2         788d-3*rhob*t11*t29*t7-5.329547801878943d-2*t11*t7-4.6279
     3         677696266674d-5*gammabb*t11*t29*t56-4.6279677696266674d-5
     4         *gammaaa*t11*t29*t54-4.6279677696266674d-5*gammaab*t11*t2
     5         9*t52-4.6279677696266674d-5*gammabb*t11*t29*t42-4.6279677
     6         696266674d-5*gammaaa*t11*t29*t38-4.6279677696266674d-5*ga
     7         mmaab*t11*t29*t32+t211+t210+t209+t208+t207+t206+t205+t204
     8         -5.481209360000001d-4*gammabb*t11*t127-5.481209360000001d
     9         -4*gammaaa*t11*t125-5.481209360000001d-4*gammaab*t11*t121
     :         )+t24*t44*(-2.3250152285696893d-3*rhob*t11*t29*t9-2.32501
     ;         52285696893d-3*rhoa*t11*t29*t9-2.7536698324946973d-2*t11*
     <         t9+1.1633333333333332d-1*t88-6.200040609519172d-3*rhoa*t1
     =         1*t29*t8-7.343119553319191d-2*t11*t8-6.200040609519171d-3
     >         *rhob*t11*t29*t7-7.343119553319191d-2*t11*t7-6.3764735554
     ?         66666d-5*gammabb*t11*t29*t56-6.376473555466666d-5*gammaaa
     @         *t11*t29*t54-6.376473555466666d-5*gammaab*t11*t29*t52-6.3
     1         76473555466666d-5*gammabb*t11*t29*t42-6.376473555466666d-
     2         5*gammaaa*t11*t29*t38-6.376473555466666d-5*gammaab*t11*t2
     3         9*t32+t219+t218+t217+t216+t215+t214+t213+t212+1.163333333
     4         3333332d-1*t128-7.5520808d-4*gammabb*t11*t127-7.5520808d-
     5         4*gammaaa*t11*t125-7.5520808d-4*gammaab*t11*t121)+t5*t6*(
     6         -1.9985804257046041d-2*t11*t29*t9-5.481209360000001d-4*ga
     7         mmaaa*t11*t29*t84-6.491760000000001d-3*gammaaa*t11*(1.111
     8         111111111111d-1*rhoa*t83+2.222222222222222d-1*t37+1.11111
     9         1111111111d-1*rhoa*rhob*(t237+t236+t235+t180+t178+t176+t1
     :         75+t170)+2.222222222222222d-1*rhob*t124)-5.32954780187894
     ;         3d-2*t11*t29*t8-5.481209360000001d-4*gammabb*t11*t29*t77-
     <         6.491760000000001d-3*gammabb*t11*(1.111111111111111d-1*rh
     =         oa*t76+2.222222222222222d-1*t41+1.111111111111111d-1*rhoa
     >         *rhob*(t234+t233+t232+t174+t173+t172+t171+t170)+2.2222222
     ?         22222222d-1*rhob*t126)-5.481209360000001d-4*gammaab*t11*t
     @         29*t71-5.329547801878944d-2*t11*t29*t7-6.491760000000001d
     1         -3*gammaab*t11*(-1.5555555555555553d+0*rhob*t69-7.7777777
     2         77777777d-1*rhoa*t69+t231+t169)-1.9541675273556128d-1*rho
     3         b*t11*t29*t65-2.3144502890117796d+0*t11*t65-5.48120936000
     4         0001d-4*gammabb*t11*t127*t29-5.481209360000001d-4*gammaaa
     5         *t11*t125*t29-5.481209360000001d-4*gammaab*t11*t121*t29+t
     6         187+t186+t185+t184+t183+t182+t167+t166+t165+t164)+(t240+t
     7         223)*t58*t89+(1.1633333333333332d-1*t129+1.16333333333333
     8         32d-1*t100)*t24*t89+3.333333333333333d-1*t49*t5*(-11*t88-
     9         11*t128)+(t229-7.8688d-1*rhob)*t5*t62+(t227+1.52567288888
     :         88887d-1*rhob)*t24*t60+(t226-1.0649196764444441d-2*rhob)*
     ;         t58*t6+t230+t228+t225+t224+t195+t194+t192+t190+t189+t163+
     <         t162+t160+t158+t157)*wght+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = 1.0d+0*(t24*(t255+t239+t220-5*t138)*
     1         t92+t5*t6*(-3.9971608514092083d-2*t11*t29*t9-1.0659095603
     2         757887d-1*t11*t29*t8-1.0659095603757889d-1*t11*t29*t7-6.4
     3         91760000000001d-3*gammaab*t11*(-7.777777777777777d-1*rhob
     4         *t69-1.5555555555555553d+0*rhoa*t69+t231+t169)-6.49176000
     5         0000001d-3*gammabb*t11*(2.222222222222222d-1*t55+1.111111
     6         111111111d-1*rhoa*rhob*(t237+t236+t235+t174+t173+t172+t17
     7         1+t170)+1.111111111111111d-1*rhob*t143+2.222222222222222d
     8         -1*rhoa*t126)-6.491760000000001d-3*gammaaa*t11*(2.2222222
     9         22222222d-1*t53+1.111111111111111d-1*rhoa*rhob*(t234+t233
     :         +t232+t180+t178+t176+t175+t170)+1.111111111111111d-1*rhob
     ;         *t141+2.222222222222222d-1*rhoa*t124)-1.0962418720000001d
     <         -3*gammabb*t11*t127*t29-1.0962418720000001d-3*gammaaa*t11
     =         *t125*t29-1.0962418720000001d-3*gammaab*t11*t121*t29+t187
     >         +t186+t185+t184+t183+t182+t167+t166+t165+t164-2.314450289
     ?         0117796d+0*t11*t139)+(t240+2.3266666666666663d-1*t138)*t5
     @         8*t89+(t229-2*t120)*t5*t62+(t227+1.1633333333333332d-1*t1
     1         20)*t24*t60-2.3333333333333334d+0*t119*t24*t60+(t226+2.32
     2         66666666666663d-1*t119)*t58*t6-7.333333333333333d+0*t128*
     3         t49*t5+t24*(t254+t253+t252+t251+t250+t249+t248+t247+t246+
     4         t219+t218+t217+t216+t215+t214+t213+t212+2.326666666666666
     5         3d-1*t128)*t44+t245+t244+t243+t242+t241+t230+t228+t225+t2
     6         24+t195+t194+t192+t190+t189+t163+t162+t160+t158+t157)*wgh
     7         t+Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(1.711111111111111d+1*t5*t57*
     1         t94+t24*(-4.2655555555555547d-1*t57+t255+t239-5*t146)*t92
     2         +t5*t6*(-1.6874680727699207d-3*rhoa*t11*t68*t9+2.66477390
     3         09394716d-2*rhoa*t11*t23*t9-4.4999148607197886d-3*rhoa*t1
     4         1*t68*t8+7.106063735838591d-2*rhoa*t11*t23*t8-6.491760000
     5         000001d-3*gammaab*t11*(t169-2.333333333333333d+0*rhoa*t69
     6         )-4.6279677696266674d-5*gammabb*t11*t56*t68-4.62796776962
     7         66674d-5*gammaaa*t11*t54*t68-4.6279677696266674d-5*gammaa
     8         b*t11*t52*t68+7.308279146666667d-4*gammabb*t11*t23*t56+7.
     9         308279146666667d-4*gammaaa*t11*t23*t54+7.308279146666667d
     :         -4*gammaab*t11*t23*t52-1.0962418720000001d-3*gammabb*t11*
     ;         t144*t29-1.0962418720000001d-3*gammaaa*t11*t142*t29-1.096
     <         2418720000001d-3*gammaab*t11*t140*t29-3.9083350547112256d
     =         -1*rhoa*t11*t139*t29-6.491760000000001d-3*gammabb*t11*(1.
     >         111111111111111d-1*rhoa*rhob*(t181+t179+t177+t174+t173+t1
     ?         72+t171+t170)+3.333333333333333d-1*rhoa*t143)-6.491760000
     @         000001d-3*gammaaa*t11*(1.111111111111111d-1*rhoa*rhob*(t1
     1         80+t178+t176+t175+t170)+3.333333333333333d-1*rhoa*t141)-3
     2         .857417148352966d+0*rhoa*rhob**6.666666666666666d-1*t11)+
     3         (t240+2.3266666666666663d-1*t146)*t58*t89-1.18032d+0*rhoa
     4         *t5*t62+2.288509333333333d-1*rhoa*t24*t60-1.5973795146666
     5         664d-2*rhoa*t58*t6-7.333333333333333d+0*t145*t49*t5+t24*(
     6         t254+t253+t252+t251+t250+t249+t248+t247+t246+t219+t218+t2
     7         17+t216+t215+t214+t213+t212+2.3266666666666663d-1*t145)*t
     8         44+t245+t244+t243+t242+t241+t225+t224+t195+t194+t192+t190
     9         +t189+t162+t160+t158+t157)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(t5*t6*(t261-6.4917600000000
     1         01d-3*t11*t84)+(4.760624000000001d-2*t11*t38+t264)*t49*t5
     2         +(t263-5.481209360000001d-4*t11*t38)*t44*t5+t24*(t262-1.5
     3         1041616d-3*t11*t38)*t44+t260+t259+t258+t257+t256)*wght+Cm
     4         at3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = 1.0d+0*(t5*t6*(t270-6.4917600000000
     1         01d-3*t11*t71)+(4.760624000000001d-2*t11*t32+t273)*t49*t5
     2         +(t272-5.481209360000001d-4*t11*t32)*t44*t5+t24*(t271-1.5
     3         1041616d-3*t11*t32)*t44+t269+t268+t267+t266+t265)*wght+Cm
     4         at3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = 1.0d+0*(t5*t6*(t279-6.4917600000000
     1         01d-3*t11*t77)+(4.760624000000001d-2*t11*t42+t282)*t49*t5
     2         +(t281-5.481209360000001d-4*t11*t42)*t44*t5+t24*(t280-1.5
     3         1041616d-3*t11*t42)*t44+t278+t277+t276+t275+t274)*wght+Cm
     4         at3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = 1.0d+0*((t261-6.491760000000001d-3*
     1         t11*t125)*t5*t6+t49*t5*(2.3803120000000005d-2*t11*t54+2.3
     2         803120000000005d-2*t11*t38+t264)+t24*t44*(-7.5520808d-4*t
     3         11*t54-7.5520808d-4*t11*t38+t262)+t283+t260+t259+t258+t25
     4         7+t256)*wght+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = 1.0d+0*((t270-6.491760000000001d-3*
     1         t11*t121)*t5*t6+t49*t5*(2.3803120000000005d-2*t11*t52+2.3
     2         803120000000005d-2*t11*t32+t273)+t24*t44*(-7.5520808d-4*t
     3         11*t52-7.5520808d-4*t11*t32+t271)+t284+t269+t268+t267+t26
     4         6+t265)*wght+Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = 1.0d+0*((t279-6.491760000000001d-3*
     1         t11*t127)*t5*t6+t49*t5*(2.3803120000000005d-2*t11*t56+2.3
     2         803120000000005d-2*t11*t42+t282)+t24*t44*(-7.5520808d-4*t
     3         11*t56-7.5520808d-4*t11*t42+t280)+t285+t278+t277+t276+t27
     4         5+t274)*wght+Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = 1.0d+0*(t5*(-5.481209360000001d-4*t
     1         11*t29*t54-6.491760000000001d-3*t11*t142)*t6+t49*t5*(4.76
     2         0624000000001d-2*t11*t54+t264)+t24*t44*(t262-1.51041616d-
     3         3*t11*t54)+t283+t260+t259+t258+t257+t256)*wght+Cmat3(iq,D
     4         3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = 1.0d+0*(t5*(-5.481209360000001d-4*t
     1         11*t29*t52-6.491760000000001d-3*t11*t140)*t6+t49*t5*(4.76
     2         0624000000001d-2*t11*t52+t273)+t24*t44*(t271-1.51041616d-
     3         3*t11*t52)+t284+t269+t268+t267+t266+t265)*wght+Cmat3(iq,D
     4         3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*(t5*(-5.481209360000001d-4*t
     1         11*t29*t56-6.491760000000001d-3*t11*t144)*t6+t49*t5*(4.76
     2         0624000000001d-2*t11*t56+t282)+t24*t44*(t280-1.51041616d-
     3         3*t11*t56)+t285+t278+t277+t276+t275+t274)*wght+Cmat3(iq,D
     4         3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
