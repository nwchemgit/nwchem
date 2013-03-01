C> \ingroup nwxc
C> @{
C>
C> \file nwxc_eval.F
C> Routines that evaluate the current density functional
C>
C> The routines in this file evaluate the current density functional.
C> They use a variety of other routines that implement specific
C> functional terms. However, those other routines are not part of the
C> API.
C>
C> @}
C
C> \ingroup nwxc_api
C> @{
C>
C> \brief Evaluate the exchange-correlation energy and the 1st order
C> partial derivatives
C>
C> Whenever density functionals are used we need the 
C> exchange-correlation energy as well as the 1st order derivatives. 
C> Hence there seemed little point in implementing a routine for the
C> exchange-correlation energy only. 
C>
C> A given density functional may consist of a variety of terms. The
C> evaluation of the functional simply loops over all the terms and
C> adds each contribution to the outputs. Therefore we first of all
C> zero all the output variables, and subsequently add all the 
C> required terms.
C>
C> Experience has shown that for these kinds of routines it is 
C> essential to be precise about the arguments are. Hence we define
C> them here is detail. First of all the inputs are:
C>
C> Rho (density)
C> -------------
C> 
C> Rho is the electron density. In closed shell cases the total 
C> electron density is used, in open shell cases the \f$\alpha\f$-electron
C> and \f$\beta\f$-electron densities are used. The density is defined as
C> \f{eqnarray*}{
C>   \rho^\alpha(\vec{r}) 
C>   &=& \sum_{i \in \{occupied\}}{\;}^\alpha\phi^*_i(\vec{r})\;^\alpha\phi_i(\vec{r}) \\\\
C>   &=& \sum_{\mu\nu}\chi_\mu(\vec{r})D^\alpha_{\mu\nu}\chi_\nu(\vec{r}) \\\\
C>   \rho^\beta(\vec{r}) 
C>   &=& \sum_{i \in \{occupied\}}{\;}^\beta\phi^*_i(\vec{r})\;^\beta\phi_i(\vec{r}) \\\\
C>   &=& \sum_{\mu\nu}\chi_\mu(\vec{r})D^\beta_{\mu\nu}\chi_\nu(\vec{r}) \\\\
C>   \rho^t(\vec{r}) &=& \rho^\alpha(\vec{r})+\rho^\beta(\vec{r})
C> \f}
C> where \f$ D^\alpha \f$ and \f$ D^\beta \f$ are the \f$\alpha\f$- and 
C> \f$\beta\f$-electron density matrices, the functions \f$^\alpha\phi\f$ and
C> \f$^\beta\phi\f$ are spin-orbitals.
C>
C> Gamma (density gradient)
C> ------------------------
C>
C> Gamma is the norm of the density gradient squared. In closed shell cases
C> norm of the gradient of the total electron density squared is used. In open
C> shell cases the various spin components are used. The different components
C> are defined as:
C> \f{eqnarray*}{
C>   \gamma^{\alpha\alpha}(\vec{r})
C>   &=& \nabla\rho^\alpha(\vec{r})\cdot\nabla\rho^\alpha(\vec{r}) \\\\
C>   \gamma^{\alpha\beta}(\vec{r})
C>   &=& \nabla\rho^\alpha(\vec{r})\cdot\nabla\rho^\beta(\vec{r}) \\\\
C>   \gamma^{\beta\beta}(\vec{r})
C>   &=& \nabla\rho^\beta(\vec{r})\cdot\nabla\rho^\beta(\vec{r}) \\\\
C>   \gamma^{tt}(\vec{r})
C>   &=& \gamma^{\alpha\alpha}(\vec{r})+\gamma^{\beta\beta}(\vec{r})
C>    +  2\gamma^{\alpha\beta}(\vec{r})
C> \f}
C> In the above the gradient of the density for one spin channel is of course
C> \f{eqnarray*}{
C>    \nabla\rho^\sigma(\vec{r})
C>    &=& \sum_{\mu\nu}\left\{
C>    \left(\nabla\chi_\mu(\vec{r})\right)D^\sigma_{\mu\nu}\chi_\nu(\vec{r})+
C>    \chi_\mu(\vec{r})D^\sigma_{\mu\nu}\left(\nabla\chi_\nu(\vec{r})\right)
C>    \right\}
C> \f}
C>
C> Tau (kinetic energy density)
C> ----------------------------
C>
C> Tau is the kinetic energy density which encapsulates 2nd derivative
C> information about the electron density. This quantity is defined
C> as:
C> \f{eqnarray*}{
C>   \tau^\sigma(\vec{r}) 
C>   &=& \frac{1}{2}\sum_{i\in\{occupied\}}|\nabla\phi_i(\vec{r})|^2
C> \f}
C>
C> Next the outputs are:
C>
C> The functional
C> --------------
C>
C> The exchange-correlation energy at every point, in general terms defined
C> as:
C> \f{eqnarray*}{
C>   f &=& 
C>   f\left(\rho^\alpha,\rho^\beta,\gamma^{\alpha\alpha},\gamma^{\alpha\beta},
C>          \gamma^{\beta\beta}\right)
C> \f}
C>
C> The derivatives of \f$f\f$ with respect to the density
C> ------------------------------------------------------
C>
C> The derivative of the exchange-correlation energy with respect to the
C> electron density (`dfdr`), defined as:
C> \f{eqnarray*}{
C>   \frac{\mathrm{d}f}{\mathrm{d}\rho^\alpha} \\\\
C>   \frac{\mathrm{d}f}{\mathrm{d}\rho^\beta}
C> \f}
C>
C> The derivatives \f$f\f$ with respect to the density gradient
C> ------------------------------------------------------------
C>
C> The routine evaluates the following derivatives related to the density
C> gradients (`dfdg`):
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}f}{\mathrm{d}\gamma^{\alpha\alpha}} \\\\
C>    \frac{\mathrm{d}f}{\mathrm{d}\gamma^{\alpha\beta}} \\\\
C>    \frac{\mathrm{d}f}{\mathrm{d}\gamma^{\beta\beta}}
C> \f}
C>
C> Conversions
C> -----------
C>
C> This routine assumes that the functional is defined in terms of one 
C> particular set of quantities. In practice there are a number of different
C> conventions that may be used and that are all equally valid. For example
C> we have described the functional in terms of the \f$\alpha\f$- and 
C> \f$\beta\f$-electron density. However, the functional could also be defined
C> in terms of the total- and spin-densities. There are a number of these
C> equally valid representations that can be used. This leads to the question
C> quantities from one representation can be expressed in terms of those of 
C> another representation. Here we will present a number of these conversions
C> that we have encountered and worked out. 
C>
C> Converting between the norm of gradient squared and the norm
C> ------------------------------------------------------------
C>
C> In the Gaussian basis set DFT code in NWChem the functional is differentiated
C> with respect to the norm of the gradient squared. In the planewave codes
C> by contrast the norm of the gradient is used. Hence we need to be able to
C> convert the derivatives from one representation to another. Some of these
C> conversions are simple, for example:
C> \f{eqnarray*}{
C>   \frac{\partial f}{\partial |\gamma^{\alpha\alpha}|}
C>   &=& 
C>   \frac{\partial f}{\partial \gamma^{\alpha\alpha}}
C>   \frac{\partial \gamma^{\alpha\alpha}}{\partial |\gamma^{\alpha\alpha}|} \\\\
C>   \frac{\partial f}{\partial |\gamma^{\alpha\alpha}|}
C>   &=& 
C>   \frac{\partial f}{\partial \gamma^{\alpha\alpha}}
C>   \frac{\partial |\gamma^{\alpha\alpha}|^2}{\partial |\gamma^{\alpha\alpha}|} \\\\
C>   \frac{\partial f}{\partial |\gamma^{\alpha\alpha}|}
C>   &=& 2
C>   \frac{\partial f}{\partial \gamma^{\alpha\alpha}}
C>   |\gamma^{\alpha\alpha}|
C> \f}
C> Similarly we have for the \f$\gamma^{\beta\beta}\f$ component:
C> \f{eqnarray*}{
C>   \frac{\partial f}{\partial |\gamma^{\beta\beta}|}
C>   &=& 2
C>   \frac{\partial f}{\partial \gamma^{\beta\beta}}
C>   |\gamma^{\beta\beta}|
C> \f}
C> To work out how to convert
C> between these representations we start by considering the functional as
C> a functional of the density matrices. Next we express the functional in
C> terms of the desired quantities
C> \f{eqnarray*}{
C>   f\left(D^\alpha,D^\beta\right)
C>   &=& f\left(\gamma^{\alpha\alpha}(D^\alpha),\gamma^{\beta\beta}(D^\beta),
C>             \gamma^{\alpha\beta}(D^\alpha,D^\beta)\right) \\\\\
C>   &=& f\left(|\gamma^{\alpha\alpha}(D^\alpha)|,
C>              |\gamma^{\beta\beta}(D^\beta)|,
C>              |\gamma^{tt}(D^\alpha,D^\beta)|\right) \\\\\
C> \f}
C> Because of the equalities we know that the derivatives with respect to a
C> particular property also have to equate to eachother as, for example
C> \f{eqnarray*}{
C>   \frac{\mathrm{d} f\left(\gamma^{\alpha\alpha},
C>         \gamma^{\beta\beta},
C>         \gamma^{\alpha\beta}\right)}{\mathrm{d}\gamma^{\alpha\alpha}}
C>   &=&
C>   \frac{\mathrm{d} f\left(|\gamma^{\alpha\alpha}|,
C>         |\gamma^{\beta\beta}|,
C>         |\gamma^{tt}|\right)}{\mathrm{d}\gamma^{\alpha\alpha}} \\\\
C>   \frac{\mathrm{d} f\left(\gamma^{\alpha\alpha},
C>         \gamma^{\beta\beta},
C>         \gamma^{\alpha\beta}\right)}{\mathrm{d}\gamma^{\beta\beta}}
C>   &=&
C>   \frac{\mathrm{d} f\left(|\gamma^{\alpha\alpha}|,
C>         |\gamma^{\beta\beta}|,
C>         |\gamma^{tt}|\right)}{\mathrm{d}\gamma^{\beta\beta}} \\\\
C>   \frac{\mathrm{d} f\left(\gamma^{\alpha\alpha},
C>         \gamma^{\beta\beta},
C>         \gamma^{\alpha\beta}\right)}{\mathrm{d}\gamma^{\alpha\beta}}
C>   &=&
C>   \frac{\mathrm{d} f\left(|\gamma^{\alpha\alpha}|,
C>         |\gamma^{\beta\beta}|,
C>         |\gamma^{tt}|\right)}{\mathrm{d}\gamma^{\alpha\beta}} \\\\
C> \f}
C> On the left hand sides we have the quantities that this routine
C> `nwxc_eval_df` calculates. Here we are particularly interested in the last
C> equation.
C> Working out the right hand sides we get:
C> \f{eqnarray*}{
C>   \frac{\mathrm{d} f\left(\gamma^{\alpha\alpha},
C>         \gamma^{\beta\beta},
C>         \gamma^{\alpha\beta}\right)}{\mathrm{d}\gamma^{\alpha\beta}}
C>   &=&
C>   \frac{\mathrm{d} f\left(|\gamma^{\alpha\alpha}|,
C>         |\gamma^{\beta\beta}|,
C>         |\gamma^{tt}|\right)}{\mathrm{d}|\gamma^{tt}|}
C>   \frac{\mathrm{d}|\gamma^{tt}|}{\mathrm{d}\gamma^{\alpha\beta}}
C> \f}
C> Using the definition
C> \f$|\gamma^{tt}|=(\gamma^{\alpha\alpha}+\gamma^{\beta\beta}+2\gamma^{\alpha\beta})^{1/2}\f$
C> we can simplify this expression further to
C> \f{eqnarray*}{
C>   \frac{\mathrm{d} f\left(\gamma^{\alpha\alpha},
C>         \gamma^{\beta\beta},
C>         \gamma^{\alpha\beta}\right)}{\mathrm{d}\gamma^{\alpha\beta}}
C>   &=&
C>   \frac{\mathrm{d} f\left(|\gamma^{\alpha\alpha}|,
C>         |\gamma^{\beta\beta}|,
C>         |\gamma^{tt}|\right)}{\mathrm{d}|\gamma^{tt}|}
C>   (\gamma^{tt})^{-1/2} \\\\
C>   \frac{\mathrm{d} f\left(\gamma^{\alpha\alpha},
C>         \gamma^{\beta\beta},
C>         \gamma^{\alpha\beta}\right)}{\mathrm{d}\gamma^{\alpha\beta}}
C>   |\gamma^{tt}|
C>   &=&
C>   \frac{\mathrm{d} f\left(|\gamma^{\alpha\alpha}|,
C>         |\gamma^{\beta\beta}|,
C>         |\gamma^{tt}|\right)}{\mathrm{d}|\gamma^{tt}|}
C> \f}
C> Interestingly this relationship does not depend on \f$\gamma^{\alpha\alpha}\f$
C> or \f$\gamma^{\beta\beta}\f$ even though these quantities are contained
C> in the norm of the total density gradient.
C>
       subroutine nwxc_eval_df(ipol,nq,rho,gamma,tau,f,dfdr,dfdg,dfdt)
       implicit none
c#include "errquit.fh"
#include "nwxcP.fh"
#include "nwxc_param.fh"
       integer ipol !< [Input] The number of spin channels
       integer nq   !< [Input] The number of points
C
       double precision rho(nq,ipol)     !< [Input] Density
       double precision gamma(nq,ipol+1) !< [Input] |Density gradient|^2
       double precision tau(nq,ipol)     !< [Input] Kinetic energy 
                                         !< density
C
       double precision f(nq)           !< [Output] Energy
       double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
       double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
       double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
C
       integer iq !< Counter over points
       integer ip !< Counter over spin channels
       integer it !< Counter of functional terms
C
       logical nwxc_is_gga  !< Is the functional a GGA
       logical nwxc_is_mgga !< Is the functional a meta-GGA
C
       integer nd1r(2) !< The number of partial derivatives wrt rho as
                       !< a function of ipol
       integer nd1g(2) !< The number of partial derivatives wrt gamma as
                       !< a function of ipol
       integer nd1t(2) !< The number of partial derivatives wrt tau as
                       !< a function of ipol
C
       integer max_param
       parameter (max_param = 10)
       double precision param(max_param)
C
       data nd1r / D1_RA,  D1_RB  /
       data nd1g / D1_GAA, D1_GBB /
       data nd1t / D1_TA,  D1_TB  /
C
       do iq = 1, nq
         f(iq) = 0.0d0
       enddo
       do ip = 1, nd1r(ipol)
         do iq = 1, nq
           dfdr(iq,ip) = 0.0d0
         enddo
       enddo
       if (nwxc_is_gga()) then
         do ip = 1, nd1g(ipol)
           do iq = 1, nq
             dfdg(iq,ip) = 0.0d0
           enddo
         enddo
       endif
       if (nwxc_is_mgga()) then
         do ip = 1, nd1t(ipol)
           do iq = 1, nq
             dfdt(iq,ip) = 0.0d0
           enddo
         enddo
       endif
C
       do it = 1, nwxc_num_f
         select case (nwxc_ids(it))
           case (NWXCP_X_SLATER)
             call nwxc_x_dirac(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,f,dfdr)
           case (NWXCP_X_B88)
             call nwxc_x_b88(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f,dfdr,dfdg)
           case (NWXCP_X_BNL)
             param(1) = nwxc_cam_gamma
             call nwxc_x_bnl(param,nwxc_rho_tol,ipol,nq,
     +                       nwxc_wghts(it),rho,f,dfdr)
           case (NWXCP_X_CAMPBE)
             param(1) = 0.8040d0
             param(2) = nwxc_cam_alpha
             param(3) = nwxc_cam_beta
             param(4) = nwxc_cam_gamma
             call nwxc_x_campbe(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdg)
           case (NWXCP_X_PBE)
             param(1) = 0.8040d0
             call nwxc_x_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f,dfdr,dfdg)
           case (NWXCP_X_REVPBE)
             param(1) = 1.245d0
             call nwxc_x_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f,dfdr,dfdg)
           case (NWXCP_X_RPBE)
             call nwxc_x_rpbe(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,gamma,f,dfdr,dfdg)
           case (NWXCP_X_WPBE)
             param(1) = nwxc_cam_gamma
             call nwxc_x_wpbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,gamma,f,dfdr,dfdg)
           case (NWXCP_X_VS98)
             call nwxc_x_vs98(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f,dfdr,dfdg,dfdt)
           case (NWXCP_C_B95)
             param(1) = 0.0031d0
             param(2) = 0.038d0
             call nwxc_c_b95(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f,dfdr,dfdg,dfdt)
           case (NWXCP_C_PW6B95)
             param(1) = 0.00262d0
             param(2) = 0.03668d0
             call nwxc_c_b95(param,nwxc_rho_tol,ipol,nq,
     +                       nwxc_wghts(it),rho,gamma,tau,
     +                       f,dfdr,dfdg,dfdt)
           case (NWXCP_C_PWB6K)
             param(1) = 0.00353d0
             param(2) = 0.04120d0
             call nwxc_c_b95(param,nwxc_rho_tol,ipol,nq,
     +                       nwxc_wghts(it),rho,gamma,tau,
     +                       f,dfdr,dfdg,dfdt)
           case (NWXCP_C_PW91LDA)
             call nwxc_c_pw91lda(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,f,dfdr)
           case (NWXCP_C_PBE)
             call nwxc_c_pbe(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f,dfdr,dfdg)
           case (NWXCP_C_VS98)
             call nwxc_c_vs98(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,gamma,tau,f,dfdr,dfdg,dfdt)
           case (NWXCP_C_VWN1)
             call nwxc_c_vwn1(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,f,dfdr)
           case (NWXCP_C_VWN1_RPA)
             call nwxc_c_vwn1_rpa(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                            rho,f,dfdr)
           case (NWXCP_C_VWN2)
             call nwxc_c_vwn2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,f,dfdr)
           case (NWXCP_C_VWN3)
             call nwxc_c_vwn3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,f,dfdr)
           case (NWXCP_C_VWN4)
             call nwxc_c_vwn4(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,f,dfdr)
           case (NWXCP_C_VWN5)
             call nwxc_c_vwn5(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,f,dfdr)
           case default
c            call errquit("nwxc_eval_df: invalid functional",
c    +                    nwxc_ids(it),UERR)
             write(*,*)"nwxc_eval_df: invalid functional",nwxc_ids(it)
         end select
       enddo
cDEBUG
c      write(*,*)'func',rho(50,1),gamma(50,1),f(50),
c    + dfdr(50,1),dfdg(50,1)
cDEBUG
C
       end
C> @}
