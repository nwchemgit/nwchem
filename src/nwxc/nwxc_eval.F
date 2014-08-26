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
C>          \gamma^{\beta\beta},\tau^\alpha,\tau^\beta\right)
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
C> The derivatives of \f$f\f$ with respect to the density gradient
C> ---------------------------------------------------------------
C>
C> The routine evaluates the following derivatives related to the density
C> gradients (`dfdg`):
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}f}{\mathrm{d}\gamma^{\alpha\alpha}} \\\\
C>    \frac{\mathrm{d}f}{\mathrm{d}\gamma^{\alpha\beta}} \\\\
C>    \frac{\mathrm{d}f}{\mathrm{d}\gamma^{\beta\beta}}
C> \f}
C>
C> The derivatives of \f$f\f$ with respect to the kinetic energy density
C> ---------------------------------------------------------------------
C>
C> The routine evaluates the following derivatives related to the kinetic
C> energy density (`dfdt`):
C> \f{eqnarray*}{
C>   \frac{\mathrm{d}f}{\mathrm{d}\tau^{\alpha}} \\\\
C>   \frac{\mathrm{d}f}{\mathrm{d}\tau^{\beta}}
C> \f}
C> where we define the kinetic energy density as 
C> \f{eqnarray*}{
C>   \tau_\sigma &=& \sum_i \langle\phi_i^\sigma|\nabla^2|\phi_i^\sigma\rangle
C> \f}
C> I.e. we assume the factor \f$1/2\f$ to be absorbed into the functional.
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
C> how quantities from one representation can be expressed in terms of those of 
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
C> \brief The functional and 1st order partial derivative evaluation
C>
C> In actual fact this routine only sets up the memory and then calls the
C> driver routine to drive the actual functional evaluation. The lack of 
C> a Fortran90 typecast drives this code structure.
C>
      subroutine nwxc_eval_df(ipol,nq,rho,gamma,tau,f,dfdr,dfdg,dfdt)
#include "nwad.fh"
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
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
c
      type(nwad_dble) :: memory_test_nwad(2)
      double precision   memory_test_dble(2)
      integer length
      integer l_rho, k_rho !< Index for density
      integer l_gam, k_gam !< Index for gradient
      integer l_tau, k_tau !< Index for kinetic energy density
      integer l_fun, k_fun !< Index for functional
      integer ipl          !< Hack version of ipol
c
c     Work out how many double precision floating point words are needed
c     to represent a single nwad_dble
c
      ipl = 2
      length = (loc(memory_test_dble(2))-loc(memory_test_dble(1)))
      length = (loc(memory_test_nwad(2))-loc(memory_test_nwad(1))
     +          +length-1)/length
      if(.not.ma_push_get(mt_dbl,length*nq*ipl,"nwad rho",l_rho,k_rho))
     +  call errquit("nwxc_eval_df: cannot allocate nwad rho",
     +    ma_sizeof(mt_dbl,length*nq*ipl,mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*(ipl+1),"nwad gamma",
     +                    l_gam,k_gam))
     +  call errquit("nwxc_eval_df: cannot allocate nwad gamma",
     +    ma_sizeof(mt_dbl,length*nq*(ipl+1),mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*ipl,"nwad tau",l_tau,k_tau))
     +  call errquit("nwxc_eval_df: cannot allocate nwad tau",
     +    ma_sizeof(mt_dbl,length*nq*ipl,mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq,"nwad fun",l_fun,k_fun))
     +  call errquit("nwxc_eval_df: cannot allocate nwad fun",
     +    ma_sizeof(mt_dbl,length*nq,mt_byte),MA_ERR)
c
      call nwxc_eval_df_driver(ipol,nq,rho,gamma,tau,dbl_mb(k_rho),
     +     dbl_mb(k_gam),dbl_mb(k_tau),dbl_mb(k_fun),f,
     +     dfdr,dfdg,dfdt)
c
      if(.not.ma_pop_stack(l_fun)) call errquit(
     +  "nwxc_eval_df: cannot deallocate nwad fun",0,MA_ERR)
      if(.not.ma_pop_stack(l_tau)) call errquit(
     +  "nwxc_eval_df: cannot deallocate nwad tau",0,MA_ERR)
      if(.not.ma_pop_stack(l_gam)) call errquit(
     +  "nwxc_eval_df: cannot deallocate nwad gam",0,MA_ERR)
      if(.not.ma_pop_stack(l_rho)) call errquit(
     +  "nwxc_eval_df: cannot deallocate nwad rho",0,MA_ERR)
      end
C>
C> @}
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Driver routine for the functional and first derivative evaluation
C>
C> This driver routine initializes the active variables and loops over the
C> functional evaluations. Afterwards it interpolates the partial derivatives
C> and returns the values of the functional and its 1st order partial
C> derivatives. Obviously, interpolating the 1st order partial derivatives is
C> particularly simply, but this scheme applies to all orders of derivatives.
C> Due to different data types being required the memory for this routine
C> needs to be allocated in the calling routine. Calling this subroutine
C> effectively implements a typecast.
C>
C> One key point is that the implementation of the functionals assumes that
C> the kinetic energy density is 
C> \f{eqnarray*}{
C>   \tau_\sigma &=& \sum_i \langle\phi_i^\sigma|\nabla^2|\phi_i^\sigma\rangle
C> \f}
C> The NWChem code that generates the kinetic energy density computes the
C> proper kinetic energy 
C> \f{eqnarray*}{
C>   \tau_\sigma &=& \frac{1}{2}
C>               \sum_i \langle\phi_i^\sigma|\nabla^2|\phi_i^\sigma\rangle
C> \f}
C> Hence we need to convert this quantity to what the functionals expect
C> when we initialize the active variables (otherwise the conversion factor
C> will affect the derivatives).
C>
      subroutine nwxc_eval_df_driver(ipol,nq,rho,gamma,tau,
     +           nwad_rho,nwad_gam,nwad_tau,nwad_f,f,
     +           dfdr,dfdg,dfdt)
      use nwad1
      implicit none
#include "nwxc_param.fh"
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of points
C
      double precision rho(nq,ipol) !< [Input] Density
      double precision gamma(nq,3)  !< [Input] |Density gradient|^2
      double precision tau(nq,ipol) !< [Input] Kinetic energy 
                                    !< density
C
      double precision f(nq)           !< [Output] Energy
      double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
      double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
      double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
c
      type(nwad_dble) :: nwad_rho(nq,2)
      type(nwad_dble) :: nwad_gam(nq,3)
      type(nwad_dble) :: nwad_tau(nq,2)
      type(nwad_dble) :: nwad_f(nq)
C
      logical  nwxc_is_lda  !< Is the functional a GGA
      logical  nwxc_is_gga  !< Is the functional a GGA
      logical  nwxc_is_mgga !< Is the functional a meta-GGA
      external nwxc_is_gga  !< Is the functional a GGA
      external nwxc_is_mgga !< Is the functional a meta-GGA
c
      integer iq, ip
c
      if (ipol.eq.1) then
        do ip = 1, 2
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
          enddo
        enddo
      else  ! ipol.eq.1
        do ip = 1, 2
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(rho(iq,ip))
          enddo
        enddo
      endif ! ipol.eq.1
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        if (ipol.eq.1) then
          do ip = 1, 3
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(0.25d0*gamma(iq,G_TT))
            enddo
          enddo
        else  ! ipol.eq.1
          do ip = 1, 3
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(gamma(iq,ip))
            enddo
          enddo
        endif ! ipol.eq.1
      endif
      if (nwxc_is_mgga()) then
        if (ipol.eq.1) then
          do ip = 1, 2
            do iq = 1, nq
              nwad_tau(iq,ip) = inactive(0.5d0*tau(iq,T_T))
            enddo
          enddo
        else  ! ipol.eq.1
          do ip = 1, 2
            do iq = 1, nq
              nwad_tau(iq,ip) = inactive(1.0d0*tau(iq,ip))
            enddo
          enddo
        endif ! ipol.eq.1
      endif
c
c     Do derivatives wrt rho
c
      do ip = 1, ipol
c       do iq = 1, nq
c         nwad_f(iq) = 0.0d0
c       enddo
        if (ipol.eq.1) then
          do iq = 1, nq
            nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T))
          enddo
        else  ! ipol.eq.1
          do iq = 1, nq
            nwad_rho(iq,ip) = active(rho(iq,ip))
          enddo
        endif ! ipol.eq.1
        call nwxc_eval_df_doit(2,nq,nwad_rho,nwad_gam,
     +       nwad_tau,nwad_f,dfdr,dfdg,dfdt)
        do iq = 1, nq
          dfdr(iq,ip) = inter_d1_dx(nwad_f(iq))
        enddo
        if (ipol.eq.1) then
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
          enddo
        else  ! ipol.eq.1
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(rho(iq,ip))
          enddo
        endif ! ipol.eq.1
      enddo ! ip = 1, ipol
c
c     Do derivatives wrt gamma
c
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        do ip = 1, ipol+1
c         do iq = 1, nq
c           nwad_f(iq) = 0.0d0
c         enddo
          if (ipol.eq.1) then
            do iq = 1, nq
              nwad_gam(iq,ip) = active(0.25d0*gamma(iq,G_TT))
            enddo
          else  ! ipol.eq.1
            do iq = 1, nq
              nwad_gam(iq,ip) = active(gamma(iq,ip))
            enddo
          endif ! ipol.eq.1
          call nwxc_eval_df_doit(2,nq,nwad_rho,nwad_gam,
     +         nwad_tau,nwad_f,dfdr,dfdg,dfdt)
          do iq = 1, nq
            dfdg(iq,ip) = inter_d1_dx(nwad_f(iq))
          enddo
          if (ipol.eq.1) then
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(0.25d0*gamma(iq,G_TT))
            enddo
          else  ! ipol.eq.1
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(gamma(iq,ip))
            enddo
          endif ! ipol.eq.1
        enddo ! ip = 1, ipol
      endif
c
c     Do derivatives wrt tau
c
      if (nwxc_is_mgga()) then
        do ip = 1, ipol
c         do iq = 1, nq
c           nwad_f(iq) = 0.0d0
c         enddo
          if (ipol.eq.1) then
            do iq = 1, nq
              nwad_tau(iq,ip) = active(0.5d0*tau(iq,T_T))
            enddo
          else  ! ipol.eq.1
            do iq = 1, nq
              nwad_tau(iq,ip) = active(1.0d0*tau(iq,ip))
            enddo
          endif ! ipol.eq.1
          call nwxc_eval_df_doit(2,nq,nwad_rho,nwad_gam,
     +         nwad_tau,nwad_f,dfdr,dfdg,dfdt)
          do iq = 1, nq
            dfdt(iq,ip) = inter_d1_dx(nwad_f(iq))
          enddo
          if (ipol.eq.1) then
            do iq = 1, nq
              nwad_tau(iq,ip) = inactive(0.5d0*tau(iq,T_T))
            enddo
          else  ! ipol.eq.1
            do iq = 1, nq
              nwad_tau(iq,ip) = inactive(1.0d0*tau(iq,ip))
            enddo
          endif ! ipol.eq.1
        enddo ! ip = 1, ipol
      endif
c
c     Finally get the functional value
c
      do iq = 1, nq
        f(iq) = value(nwad_f(iq))
      enddo
C     
      end
C>
C> \brief Evaluate the functional and its first order derivatives
C>
C> Using automatic differentiation the functional and its first order
C> derivatives are calculated for one direction. The loop over directions
C> and the correct interpolation of the partial derivatives has to be handled
C> in a driver routine.
C>
      subroutine nwxc_eval_df_doit(ipol,nq,rho,gamma,tau,f,
     +                             dfdr,dfdg,dfdt)
      use nwad1
      implicit none
#include "intf_nwxc_c_p91.fh"
#include "intf_nwxc_c_pw91lda.fh"
c
#include "intf_nwxc_x_p91.fh"
c
#include "nwxcP.fh"
#include "nwxc_param.fh"
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of points
C
      type(nwad_dble)::rho(nq,ipol)              !< [Input] Density
      type(nwad_dble)::gamma(nq,ipol*(ipol+1)/2) !< [Input] |Density gradient|^2
      type(nwad_dble)::tau(nq,ipol)              !< [Input] Kinetic energy 
                                                 !< density
C
      type(nwad_dble)::f(nq)           !< [Output] Energy
      double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
      double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
      double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
C
      integer iq !< Counter over points
      integer ip !< Counter over spin channels
      integer it !< Counter of functional terms
      integer n  !< Counter
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
      external nwxc_k_dirac
      external nwxc_k_becke88
      external nwxc_k_pbe96
c
      integer max_param
      parameter (max_param = 50)
      double precision param(max_param)
      double precision Axlsda
      parameter (Axlsda = -0.9305257363491d0 )
C
      data nd1r / D1_RA,  D1_RB  /
      data nd1g / D1_GAA, D1_GBB /
      data nd1t / D1_TA,  D1_TB  /
C
      do iq = 1, nq
        f(iq) = 0.0d0
      enddo
c     do ip = 1, nd1r(ipol)
c       do iq = 1, nq
c         dfdr(iq,ip) = 0.0d0
c       enddo
c     enddo
c     if (nwxc_is_gga().or.nwxc_is_mgga()) then
c       do ip = 1, nd1g(ipol)
c         do iq = 1, nq
c           dfdg(iq,ip) = 0.0d0
c         enddo
c       enddo
c     endif
c     if (nwxc_is_mgga()) then
c       do ip = 1, nd1t(ipol)
c         do iq = 1, nq
c           dfdt(iq,ip) = 0.0d0
c         enddo
c       enddo
c     endif
cDEBUG
c     if (ipol.eq.2.and.nq.eq.3) then
c       write(*,*)"point 1"
c       call print(6,rho(1,1))
c       call print(6,rho(1,2))
c       call print(6,gamma(1,1))
c       call print(6,gamma(1,2))
c       call print(6,gamma(1,3))
c       call print(6,tau(1,1))
c       call print(6,tau(1,2))
c       write(*,*)
c       write(*,*)"point 2"
c       call print(6,rho(2,1))
c       call print(6,rho(2,2))
c       call print(6,gamma(2,1))
c       call print(6,gamma(2,2))
c       call print(6,gamma(2,3))
c       call print(6,tau(2,1))
c       call print(6,tau(2,2))
c       write(*,*)
c       write(*,*)"point 3"
c       call print(6,rho(3,1))
c       call print(6,rho(3,2))
c       call print(6,gamma(3,1))
c       call print(6,gamma(3,2))
c       call print(6,gamma(3,3))
c       call print(6,tau(3,1))
c       call print(6,tau(3,2))
c       write(*,*)
c     endif
cDEBUG
C
      do it = 1, nwxc_num_f
        select case (nwxc_ids(it))
          case (NWXCP_X_SLATER)
            call nwxc_x_dirac(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                        rho,f)
          case (NWXCP_X_B86B)
            call nwxc_x_b86b(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_B88)
            call nwxc_x_b88(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B97)
            param(1) = 2.0d0
            param(2) = 0.80940d+00
            param(3) = 0.50730d+00
            param(4) = 0.74810d+00
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B97_1)
            param(1) = 2.0d0
            param(2) = 0.789518d+00
            param(3) = 0.573805d+00
            param(4) = 0.660975d+00
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B97_2)
            param(1) = 2.0d0
            param(2) = 0.827642D+00
            param(3) = 0.478400D-01
            param(4) = 0.176125D+01
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B97_3)
            param(1) = 4.0d0
            param(2) = 7.334648D-01
            param(3) = 2.925270D-01
            param(4) = 3.338789D+00
            param(5) =-1.051158D+01
            param(6) = 1.060907D+01
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B97_D)
            param(1) = 2.0d0
            param(2) = 1.086620d+0
            param(3) =-0.521270d+00
            param(4) = 3.254290d+00
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B97_G)
            param(1) = 2.0d0
            param(2) = 1.1068d0
            param(3) =-0.8765d0
            param(4) = 4.2639d0
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_B98)
            param(1) = 2.0d0
            param(2) = 0.790194d00
            param(3) = 0.400271d00
            param(4) = 0.832857d00
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_BNL)
            param(1) = nwxc_cam_gamma
            call nwxc_x_bnl(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,f)
          case (NWXCP_X_CAMB88)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camb88(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMLSD)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camlsd(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,f)
          case (NWXCP_X_CAMPBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            call nwxc_x_campbe(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMREVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            call nwxc_x_campbe(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMRPBE)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camrpbe(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMS12G)
            param(1) = 1.03323556d0
            param(2) = 0.757d0
            param(3) = 0.00417251d0
            param(4) = 0.00115216d0
            param(5) = 0.00706184d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            call nwxc_x_cams12(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMS12H)
            param(1) = 1.02149642d0
            param(2) = 0.757d0
            param(3) = 0.00825905d0
            param(4) = 0.00235804d0
            param(5) = 0.00654977d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            call nwxc_x_cams12(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_DLDF)
            call nwxc_x_dldf(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f)
          case (NWXCP_X_FT97)
            call nwxc_x_ft97(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_GILL)
            call nwxc_x_gill(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_HCTH)
            param(1) = 4.0d0
            param(2) = 0.109320d+01
            param(3) =-0.744056d+00
            param(4) = 0.559920d+01
            param(5) =-0.678549d+01
            param(6) = 0.449357d+01
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_HCTH120)
            param(1) = 4.0d0
            param(2) = 1.09163d0
            param(3) =-0.74720d0
            param(4) = 5.07830d0
            param(5) =-4.10750d0
            param(6) = 1.17170d0
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_HCTH147)
            param(1) = 4.0d0
            param(2) = 1.09025d0
            param(3) =-0.79920d0
            param(4) = 5.57210d0
            param(5) =-5.86760d0
            param(6) = 3.04540d0
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_HCTH407)
            param(1) = 4.0d0
            param(2) = 1.08184d0
            param(3) =-0.5183d0
            param(4) = 3.4256d0
            param(5) =-2.6290d0
            param(6) = 2.2886d0
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_HCTH407P)
            param(1) = 4.0d0
            param(2) = 1.08018D0
            param(3) =-0.4117D0
            param(4) = 2.4368D0
            param(5) = 1.3890D0
            param(6) =-1.3529D0
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_HCTH_A)
            param(1) = 2.0d0
            param(2) = 0.109878d+01
            param(3) =-0.251173d+01
            param(4) = 0.156233d-01
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_HCTHP14)
            param(1) = 4.0d0
            param(2) = 0.103161d+01
            param(3) =-0.360781d+00
            param(4) = 0.351994d+01
            param(5) =-0.495944d+01
            param(6) = 0.241165d+01
            call nwxc_x_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_M05)
            param( 1) =    0.08151d0
            param( 2) =   -0.43956d0
            param( 3) =   -3.22422d0
            param( 4) =    2.01819d0
            param( 5) =    8.79431d0
            param( 6) =   -0.00295d0
            param( 7) =    9.82029d0
            param( 8) =   -4.82351d0
            param( 9) =  -48.17574d0
            param(10) =    3.64802d0
            param(11) =   34.02248d0
            call nwxc_x_m05(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_X_M05_2X)
            param( 1) =   -0.56833d0
            param( 2) =   -1.30057d0
            param( 3) =    5.50070d0
            param( 4) =    9.06402d0
            param( 5) =  -32.21075d0
            param( 6) =  -23.73298d0
            param( 7) =   70.22996d0
            param( 8) =   29.88614d0
            param( 9) =  -60.25778d0
            param(10) =  -13.22205d0
            param(11) =   15.23694d0
            call nwxc_x_m05(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_X_M06)
            param( 1) =  1.422057D-01*Axlsda
            param( 2) =  7.370319D-04*Axlsda
            param( 3) = -1.601373D-02*Axlsda
            param( 4) =  0.000000D+00
            param( 5) =  0.000000D+00
            param( 6) =  0.000000D+00
            param( 7) =  5.877943D-01
            param( 8) = -1.371776D-01
            param( 9) =  2.682367D-01
            param(10) = -2.515898D+00
            param(11) = -2.978892D+00
            param(12) =  8.710679D+00
            param(13) =  1.688195D+01
            param(14) = -4.489724D+00
            param(15) = -3.299983D+01
            param(16) = -1.449050D+01
            param(17) =  2.043747D+01
            param(18) =  1.256504D+01
            call nwxc_x_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_X_M06_HF)
            param( 1) = -1.179732D-01*Axlsda
            param( 2) = -2.500000D-03*Axlsda
            param( 3) = -1.180065D-02*Axlsda
            param( 4) =  0.000000D+00
            param( 5) =  0.000000D+00
            param( 6) =  0.000000D+00
            param( 7) =  1.179732D-01
            param( 8) = -1.066708D+00
            param( 9) = -1.462405D-01
            param(10) =  7.481848D+00
            param(11) =  3.776679D+00
            param(12) = -4.436118D+01
            param(13) = -1.830962D+01
            param(14) =  1.003903D+02
            param(15) =  3.864360D+01
            param(16) = -9.806018D+01
            param(17) = -2.557716D+01
            param(18) =  3.590404D+01
            call nwxc_x_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_X_M06_L)
            param( 1) =  6.012244D-01*Axlsda
            param( 2) =  4.748822D-03*Axlsda
            param( 3) = -8.635108D-03*Axlsda
            param( 4) = -9.308062D-06*Axlsda
            param( 5) =  4.482811D-05*Axlsda
            param( 6) =  0.000000D+00
            param( 7) =  3.987756D-01
            param( 8) =  2.548219D-01
            param( 9) =  3.923994D-01
            param(10) = -2.103655D+00
            param(11) = -6.302147D+00
            param(12) =  1.097615D+01
            param(13) =  3.097273D+01
            param(14) = -2.318489D+01
            param(15) = -5.673480D+01
            param(16) =  2.160364D+01
            param(17) =  3.421814D+01
            param(18) = -9.049762D+00
            call nwxc_x_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_X_M06_2X)
            param( 1) =  4.600000D-01
            param( 2) = -2.206052D-01
            param( 3) = -9.431788D-02
            param( 4) =  2.164494D+00
            param( 5) = -2.556466D+00
            param( 6) = -1.422133D+01
            param( 7) =  1.555044D+01
            param( 8) =  3.598078D+01
            param( 9) = -2.722754D+01
            param(10) = -3.924093D+01
            param(11) =  1.522808D+01
            param(12) =  1.522227D+01
            call nwxc_x_m06_2x(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_X_M08_HX)
c           parameters A
            param(01) =  1.3340172D+00
            param(02) = -9.4751087D+00
            param(03) = -1.2541893D+01
            param(04) =  9.1369974D+00
            param(05) =  3.4717204D+01
            param(06) =  5.8831807D+01
            param(07) =  7.1369574D+01
            param(08) =  2.3312961D+01
            param(09) =  4.8314679D+00
            param(10) = -6.5044167D+00
            param(11) = -1.4058265D+01
            param(12) =  1.2880570D+01
c           parameters B
            param(13) = -8.5631823D-01
            param(14) =  9.2810354D+00
            param(15) =  1.2260749D+01
            param(16) = -5.5189665D+00
            param(17) = -3.5534989D+01
            param(18) = -8.2049996D+01
            param(19) = -6.8586558D+01
            param(20) =  3.6085694D+01
            param(21) = -9.3740983D+00
            param(22) = -5.9731688D+01
            param(23) =  1.6587868D+01
            param(24) =  1.3993203D+01
c           parameters C and D
            do n = 25, 48
              param(n) = 0.0d0
            enddo
            call nwxc_x_m08(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_X_M08_SO)
c           parameters A
            param(01) = -3.4888428D-01
            param(02) = -5.8157416D+00
            param(03) =  3.7550810D+01
            param(04) =  6.3727406D+01
            param(05) = -5.3742313D+01
            param(06) = -9.8595529D+01
            param(07) =  1.6282216D+01
            param(08) =  1.7513468D+01
            param(09) = -6.7627553D+00
            param(10) =  1.1106658D+01
            param(11) =  1.5663545D+00
            param(12) =  8.7603470D+00
c           parameters B
            param(13) =  7.8098428D-01
            param(14) =  5.4538178D+00
            param(15) = -3.7853348D+01
            param(16) = -6.2295080D+01
            param(17) =  4.6713254D+01
            param(18) =  8.7321376D+01
            param(19) =  1.6053446D+01
            param(20) =  2.0126920D+01
            param(21) = -4.0343695D+01
            param(22) = -5.8577565D+01
            param(23) =  2.0890272D+01
            param(24) =  1.0946903D+01
c           parameters C and D
            do n = 25, 48
              param(n) = 0.0d0
            enddo
            call nwxc_x_m08(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_X_M11)
c           parameters A
            param(01) = -0.18399900D+00
            param(02) = -1.39046703D+01
            param(03) =  1.18206837D+01
            param(04) =  3.10098465D+01
            param(05) = -5.19625696D+01
            param(06) =  1.55750312D+01
            param(07) = -6.94775730D+00
            param(08) = -1.58465014D+02
            param(09) = -1.48447565D+00
            param(10) =  5.51042124D+01
            param(11) = -1.34714184D+01
            param(12) =  0.000000D+00
c           parameters B
            param(13) =  0.75599900D+00
            param(14) =  1.37137944D+01
            param(15) = -1.27998304D+01
            param(16) = -2.93428814D+01
            param(17) =  5.91075674D+01
            param(18) = -2.27604866D+01
            param(19) = -1.02769340D+01
            param(20) =  1.64752731D+02
            param(21) =  1.85349258D+01
            param(22) = -5.56825639D+01
            param(23) =  7.47980859D+00
            param(24) =  0.000000D+00
c           parameters C and D
            do n = 25, 48
              param(n) = 0.0d0
            enddo
            call nwxc_x_m11(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_X_M11_L)
c           parameters A
            param(01) =  8.121131D-01
            param(02) =  1.738124D+01
            param(03) =  1.154007D+00
            param(04) =  6.869556D+01
            param(05) =  1.016864D+02
            param(06) = -5.887467D+00
            param(07) =  4.517409D+01
            param(08) = -2.773149D+00
            param(09) = -2.617211D+01
            param(10) =  0.000000D+00
            param(11) =  0.000000D+00
            param(12) =  0.000000D+00
c           parameters B
            param(13) =  1.878869D-01
            param(14) = -1.653877D+01
            param(15) =  6.755753D-01
            param(16) = -7.567572D+01
            param(17) = -1.040272D+02
            param(18) =  1.831853D+01
            param(19) = -5.573352D+01
            param(20) = -3.520210D+00
            param(21) =  3.724276D+01
            param(22) =  0.000000D+00
            param(23) =  0.000000D+00
            param(24) =  0.000000D+00
c           parameters C
            param(25) = -4.386615D-01
            param(26) = -1.214016D+02
            param(27) = -1.393573D+02
            param(28) = -2.046649D+00
            param(29) =  2.804098D+01
            param(30) = -1.312258D+01
            param(31) = -6.361819D+00
            param(32) = -8.055758D-01
            param(33) =  3.736551D+00
            param(34) =  0.000000D+00
            param(35) =  0.000000D+00
            param(36) =  0.000000D+00
c           parameters D
            param(37) =  1.438662D+00
            param(38) =  1.209465D+02
            param(39) =  1.328252D+02
            param(40) =  1.296355D+01
            param(41) =  5.854866D+00
            param(42) = -3.378162D+00
            param(43) = -4.423393D+01
            param(44) =  6.844475D+00
            param(45) =  1.949541D+01
            param(46) =  0.000000D+00
            param(47) =  0.000000D+00
            param(48) =  0.000000D+00
            call nwxc_x_m11(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_X_MPW91)
            param(1) = 3.72d0
            param(2) = 0.00426D0
            call nwxc_x_pw91(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_OPT)
            call nwxc_x_opt(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_PW6B95)
            param(1) = 0.00538d0
            param(2) = 1.7382d0
            param(3) = 3.8901d0
            call nwxc_x_pw6(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_PWB6K)
            param(1) = 0.00539d0
            param(2) = 1.7077d0
            param(3) = 4.0876d0
            call nwxc_x_pw6(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_PW91)
            param(1) = 4.0d0
            param(2) = 0.0042D0
            call nwxc_x_pw91(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_PBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            call nwxc_x_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_PBESOL)
            param(1) = 0.8040d0
            param(2) = 10.0d0/81.0d0
            call nwxc_x_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_REVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            call nwxc_x_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_X_RPBE)
            call nwxc_x_rpbe(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_PKZB)
            call nwxc_x_pkzb99(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,tau,f)
          case (NWXCP_X_S12G)
            param(1) = 1.03842032d0
            param(2) = 0.757d0
            param(3) = 0.00403198d0
            param(4) = 0.00104596d0
            param(5) = 0.00594635d0
            call nwxc_x_s12(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_S12H)
            param(1) = 1.02543951d0
            param(2) = 0.757d0
            param(3) = 0.00761554d0
            param(4) = 0.00211063d0
            param(5) = 0.00604672d0
            call nwxc_x_s12(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SOGGA)
            param(1)  =  0.5d0
            param(2)  =  0.276d0
            param(3)  =  0.0d0
            param(4)  =  0.0d0
            param(5)  =  0.0d0
            param(6)  =  0.0d0
            param(7)  =  0.5d0
            param(8)  =  0.276d0
            param(9)  =  0.0d0
            param(10) =  0.0d0
            param(11) =  0.0d0
            param(12) =  0.0d0
            call nwxc_x_sogga(param,nwxc_rho_tol,ipol,nq,
     +                        nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SOGGA11)
            param(1)  =  0.5d0
            param(2)  = -2.95535d0
            param(3)  =  15.7974d0
            param(4)  = -91.1804d0
            param(5)  =  96.2030d0
            param(6)  =  0.18683d0
            param(7)  =  0.50000d0
            param(8)  =  3.50743d0
            param(9)  = -12.9523d0
            param(10) =  49.7870d0
            param(11) = -33.2545d0
            param(12) = -11.1396d0
            call nwxc_x_sogga(param,nwxc_rho_tol,ipol,nq,
     +                        nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SOGGA11_X)
            param(1)  =  2.99250d-01
            param(2)  =  3.21638d+00
            param(3)  = -3.55605d+00
            param(4)  =  7.65852d+00
            param(5)  = -1.12830d+01
            param(6)  =  5.25813d+00
            param(7)  =  2.99250d-01
            param(8)  = -2.88595d+00
            param(9)  =  3.23617d+00
            param(10) = -2.45393d+00
            param(11) = -3.75495d+00
            param(12) =  3.96613d+00
            call nwxc_x_sogga(param,nwxc_rho_tol,ipol,nq,
     +                        nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SSB_D)
            call nwxc_x_ssbD_1(nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_TPSS)
            call nwxc_x_tpss03(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,tau,f)
          case (NWXCP_X_WPBE)
            param(1) = nwxc_cam_gamma
            call nwxc_x_wpbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_X_VS98)
            param(1) =  -9.800683d-01
            param(2) =  -3.556788d-03
            param(3) =   6.250326d-03
            param(4) =  -2.354518d-05
            param(5) =  -1.282732d-04
            param(6) =   3.574822d-04
            call nwxc_x_vs98(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f)
          case (NWXCP_C_B95)
            param(1) = 0.0031d0
            param(2) = 0.038d0
            call nwxc_c_b95(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_B97)
            param(1) = 2.0d0
            param(2) = 0.17370d+00
            param(3) = 0.94540d+00
            param(4) = 0.23487d+01
            param(5) = 0.74710d+00
            param(6) =-0.24868d+01
            param(7) =-0.45961d+01
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_B97_1)
            param(1) = 2.0d0
            param(2) = 0.820011d-01
            param(3) = 0.955689d+00
            param(4) = 0.271681d+01
            param(5) = 0.788552d+00
            param(6) =-0.287103d+01
            param(7) =-0.547869d+01
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_B97_2)
            param(1) = 2.0d0
            param(2) = 0.585808D+00
            param(3) = 0.999849D+00
            param(4) =-0.691682D+00
            param(5) = 0.140626D+01
            param(6) = 0.394796D+00
            param(7) =-0.744060D+01
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_B97_3)
            param(1)  = 4.0d0
            param(2)  = 5.623649D-01
            param(3)  = 1.133830D+00
            param(4)  =-1.322980D+00
            param(5)  =-2.811967D+00
            param(6)  = 6.359191D+00
            param(7)  = 7.431302D+00
            param(8)  =-7.464002D+00
            param(9)  =-1.969342D+00
            param(10) = 1.827082D+00
            param(11) =-1.174423D+01
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_B97_D)
            param(1) = 2.0d0
            param(2) = 0.22340d+00
            param(3) = 0.690410d+00
            param(4) =-1.562080d+00
            param(5) = 6.302700d00
            param(6) = 1.942930d+0
            param(7) =-14.97120d+00
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_B97_G)
            param(1) = 2.0d0
            param(2) = 0.4883d0
            param(3) = 0.7961d0
            param(4) =-2.117d0
            param(5) = 5.7060d0
            param(6) = 2.3235d0
            param(7) =-14.9820d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_B98)
            param(1) = 2.0d0
            param(2) =-0.120163d00
            param(3) = 0.934715d00
            param(4) = 2.82332d0
            param(5) = 1.14105d0
            param(6) =-2.59412d0
            param(7) =-5.33398d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_DLDF)
            call nwxc_c_dldf(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f)
          case (NWXCP_C_FT97)
            call nwxc_c_ft97(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_C_HCTH)
            param(1)  = 4.0d0
            param(2)  = 0.222601d0
            param(3)  = 0.729974d0
            param(4)  =-3.38622d-002
            param(5)  = 3.352870d0
            param(6)  =-1.25170d-002
            param(7)  =-11.543d0
            param(8)  =-0.802496d0
            param(9)  = 8.085640d0
            param(10) = 1.553960d0
            param(11) =-4.478570d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_HCTH120)
            param(1)  = 4.0d0
            param(2)  = 0.48951d0
            param(3)  = 0.51473d0
            param(4)  =-0.26070d0
            param(5)  = 6.92980d0
            param(6)  = 0.43290d0
            param(7)  =-24.7070d0
            param(8)  =-1.99250d0
            param(9)  = 23.1100d0
            param(10) = 2.48530d0
            param(11) =-11.3230d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_HCTH147)
            param(1)  = 4.0d0
            param(2)  = 0.56258d0
            param(3)  = 0.54235d0
            param(4)  =-1.71000d-002
            param(5)  = 7.01460d0
            param(6)  =-1.30640d0
            param(7)  =-28.3820d0
            param(8)  = 1.05750d0
            param(9)  = 35.0330d0
            param(10) = 0.88540d0
            param(11) =-20.4280d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_HCTH407)
            param(1)  = 4.0d0
            param(2)  = 1.18777d0
            param(3)  = 0.58908d0
            param(4)  =-2.40290d0
            param(5)  = 4.42370d0
            param(6)  = 5.61740d0
            param(7)  =-19.2220d0
            param(8)  =-9.17920d0
            param(9)  = 42.5720d0
            param(10) = 6.24800d0
            param(11) =-42.0050d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_HCTH407P)
            param(1)  = 4.0d0
            param(2)  = 0.80302d0
            param(3)  = 0.73604d0
            param(4)  =-1.04790d0
            param(5)  = 3.02700d0
            param(6)  = 4.98070d0
            param(7)  =-10.0750d0
            param(8)  =-12.8900d0
            param(9)  = 20.6110d0
            param(10) = 9.64460d0
            param(11) =-29.4180d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_HCTH_A)
            param(1)  = 4.0d0
            param(2)  = 1.36823d-002
            param(3)  = 0.836897d0
            param(4)  = 0.268920d0
            param(5)  = 1.720510d0
            param(6)  =-0.550769d0
            param(7)  =-2.784980d0
            param(8)  = 1.039470d0
            param(9)  =-4.575040d0
            param(10) = 0.000000d0
            param(11) = 0.000000d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_HCTHP14)
            param(1)  = 4.0d0
            param(2)  = 2.82414d0     
            param(3)  = 8.21827d-002
            param(4)  = 3.18843d-002
            param(5)  = 4.56466d0     
            param(6)  =-1.78512d0     
            param(7)  =-13.5529d0     
            param(8)  = 2.39795d0     
            param(9)  = 13.3820d0     
            param(10) =-0.876909d0     
            param(11) =-3.174930d0
            call nwxc_c_b97(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_M05)
            param( 1) =   1.00000d0
            param( 2) =   3.78569d0
            param( 3) = -14.15261d0
            param( 4) =  -7.46589d0
            param( 5) =  17.94491d0
            param( 6) =   1.00000d0
            param( 7) =   3.77344d0
            param( 8) = -26.04463d0
            param( 9) =  30.69913d0
            param(10) =  -9.22695d0
            call nwxc_c_m05(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_M05_2X)
            param( 1) =   1.00000d0
            param( 2) =   1.09297d0
            param( 3) =  -3.79171d0
            param( 4) =   2.82810d0
            param( 5) = -10.58909d0
            param( 6) =   1.00000d0
            param( 7) =  -3.05430d0
            param( 8) =   7.61854d0
            param( 9) =   1.47665d0
            param(10) = -11.92365d0
            call nwxc_c_m05(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_M06)
c
            param( 1) =  -2.741539D+00
            param( 2) =  -6.720113D-01
            param( 3) =  -7.932688D-02
            param( 4) =   1.918681D-03
            param( 5) =  -2.032902D-03
            param( 6) =   0.000000D+00
            param( 7) =   4.905945D-01
            param( 8) =  -1.437348D-01
            param( 9) =   2.357824D-01
            param(10) =   1.871015D-03
            param(11) =  -3.788963D-03
            param(12) =   0.000000D+00
c
            param(13) =   3.741539D+00
            param(14) =   2.187098D+02
            param(15) =  -4.531252D+02
            param(16) =   2.936479D+02
            param(17) =  -6.287470D+01
            param(18) =   5.094055D-01
            param(19) =  -1.491085D+00
            param(20) =   1.723922D+01
            param(21) =  -3.859018D+01
            param(22) =   2.845044D+01
c
            call nwxc_c_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_M06_HF)
c
            param( 1) =  -6.746338D-01
            param( 2) =  -1.534002D-01
            param( 3) =  -9.021521D-02
            param( 4) =  -1.292037D-03
            param( 5) =  -2.352983D-04
            param( 6) =   0.000000D+00
            param( 7) =   8.976746D-01
            param( 8) =  -2.345830D-01
            param( 9) =   2.368173D-01
            param(10) =  -9.913890D-04
            param(11) =  -1.146165D-02
            param(12) =   0.000000D+00
c
            param(13) =   1.674634D+00
            param(14) =   5.732017D+01
            param(15) =   5.955416D+01
            param(16) =  -2.311007D+02
            param(17) =   1.255199D+02
            param(18) =   1.023254D-01
            param(19) =  -2.453783D+00
            param(20) =   2.913180D+01
            param(21) =  -3.494358D+01
            param(22) =   2.315955D+01
c
            call nwxc_c_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_M06_L)
c
            param( 1) =   3.957626D-01
            param( 2) =  -5.614546D-01
            param( 3) =   1.403963D-02
            param( 4) =   9.831442D-04
            param( 5) =  -3.577176D-03
            param( 6) =   0.000000D+00
            param( 7) =   4.650534D-01
            param( 8) =   1.617589D-01
            param( 9) =   1.833657D-01
            param(10) =   4.692100D-04
            param(11) =  -4.990573D-03
            param(12) =   0.000000D+00
c
            param(13) =   6.042374D-01
            param(14) =   1.776783D+02
            param(15) =  -2.513252D+02
            param(16) =   7.635173D+01
            param(17) =  -1.255699D+01
            param(18) =   5.349466D-01
            param(19) =   5.396620D-01
            param(20) =  -3.161217D+01
            param(21) =   5.149592D+01
            param(22) =  -2.919613D+01
c
            call nwxc_c_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_M06_2X)
c
            param( 1) =   1.166404D-01
            param( 2) =  -9.120847D-02
            param( 3) =  -6.726189D-02
            param( 4) =   6.720580D-05
            param( 5) =   8.448011D-04
            param( 6) =   0.000000D+00
            param( 7) =   6.902145D-01
            param( 8) =   9.847204D-02
            param( 9) =   2.214797D-01
            param(10) =  -1.968264D-03
            param(11) =  -6.775479D-03
            param(12) =   0.000000D+00
c
            param(13) =   8.833596D-01
            param(14) =   3.357972D+01
            param(15) =  -7.043548D+01
            param(16) =   4.978271D+01
            param(17) =  -1.852891D+01
            param(18) =   3.097855D-01
            param(19) =  -5.528642D+00
            param(20) =   1.347420D+01
            param(21) =  -3.213623D+01
            param(22) =   2.846742D+01
c
            call nwxc_c_m06(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,tau,f)
          case (NWXCP_C_M08_HX)
c           parameters A
            param(1)  =    1.0000000D+00
            param(2)  =   -4.0661387D-01
            param(3)  =   -3.3232530D+00
            param(4)  =    1.5540980D+00
            param(5)  =    4.4248033D+01
            param(6)  =   -8.4351930D+01
            param(7)  =   -1.1955581D+02
            param(8)  =    3.9147081D+02
            param(9)  =    1.8363851D+02
            param(10) =   -6.3268223D+02
            param(11) =   -1.1297403D+02
            param(12) =    3.3629312D+02
c           parameters B
            param(13) =    1.3812334D+00
            param(14) =   -2.4683806D+00
            param(15) =   -1.1901501D+01
            param(16) =   -5.4112667D+01
            param(17) =    1.0055846D+01
            param(18) =    1.4800687D+02
            param(19) =    1.1561420D+02
            param(20) =    2.5591815D+02
            param(21) =    2.1320772D+02
            param(22) =   -4.8412067D+02
            param(23) =   -4.3430813D+02
            param(24) =    5.6627964D+01
            call nwxc_c_m11(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M08_SO)
c           parameters A
            param(1)  =   1.0000000D+00
            param(2)  =   0.0000000D+00
            param(3)  =  -3.9980886D+00
            param(4)  =   1.2982340D+01
            param(5)  =   1.0117507D+02
            param(6)  =  -8.9541984D+01
            param(7)  =  -3.5640242D+02
            param(8)  =   2.0698803D+02
            param(9)  =   4.6037780D+02
            param(10) =  -2.4510559D+02
            param(11) = -1.9638425D+02
            param(12) =  1.1881459D+02
c           parameters B
            param(13) =   1.0000000D+00
            param(14) =  -4.4117403D+00
            param(15) =  -6.4128622D+00
            param(16) =   4.7583635D+01
            param(17) =   1.8630053D+02
            param(18) =  -1.2800784D+02
            param(19) =  -5.5385258D+02
            param(20) =   1.3873727D+02
            param(21) =   4.1646537D+02
            param(22) =  -2.6626577D+02
            param(23) =   5.6676300D+01
            param(24) =   3.1673746D+02
            call nwxc_c_m11(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M11)
c           parameters A
            param(1)  =  1.0000000D+00
            param(2)  =  0.0000000D+00
            param(3)  = -3.8933250D+00
            param(4)  = -2.1688455D+00
            param(5)  =  9.3497200D+00
            param(6)  = -1.9845140D+01
            param(7)  =  2.3455253D+00
            param(8)  =  7.9246513D+01
            param(9)  =  9.6042757D+00
            param(10) = -6.7856719D+01
            param(11) = -9.1841067D+00
            param(12) =  0.0000000D+00
c           parameters B
            param(13) =  7.2239798D-01
            param(14) =  4.3730564D-01
            param(15) = -1.6088809D+01
            param(16) = -6.5542437D+01
            param(17) =  3.2057230D+01
            param(18) =  1.8617888D+02
            param(19) =  2.0483468D+01
            param(20) = -7.0853739D+01
            param(21) =  4.4483915D+01
            param(22) = -9.4484747D+01
            param(23) = -1.1459868D+02
            param(24) =  0.0000000D+00
            call nwxc_c_m11(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M11_L)
c           parameters A
            param(1)  =  1.000000D+00
            param(2)  =  0.000000D+00
            param(3)  =  2.750880D+00
            param(4)  = -1.562287D+01
            param(5)  =  9.363381D+00
            param(6)  =  2.141024D+01
            param(7)  = -1.424975D+01
            param(8)  = -1.134712D+01
            param(9)  =  1.022365D+01
            param(10) =  0.000000D+00
            param(11) =  0.000000D+00
            param(12) =  0.000000D+00
c           parameters B
            param(13) =  1.000000D+00
            param(14) = -9.082060D+00
            param(15) =  6.134682D+00
            param(16) = -1.333216D+01
            param(17) = -1.464115D+01
            param(18) =  1.713143D+01
            param(19) =  2.480738D+00
            param(20) = -1.007036D+01
            param(21) = -1.117521D-01
            param(22) =  0.000000D+00
            param(23) =  0.000000D+00
            param(24) =  0.000000D+00
            call nwxc_c_m11(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_OP)
            param(1) = 2.3670d0
            call nwxc_c_op(nwxc_k_becke88,param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_OPT)
            call nwxc_c_opt(nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PW6B95)
            param(1) = 0.00262d0
            param(2) = 0.03668d0
            call nwxc_c_b95(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_PWB6K)
            param(1) = 0.00353d0
            param(2) = 0.04120d0
            call nwxc_c_b95(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_PW91LDA)
            call nwxc_c_pw91lda(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_LYP)
            call nwxc_c_lyp(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_PZ81)
            call nwxc_c_perdew81(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,f)
          case (NWXCP_C_P86)
            call nwxc_c_perdew86(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,gamma,f)
          case (NWXCP_C_P91)
            call nwxc_c_p91(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_P91_VWN5)
            call nwxc_c_p91_vwn5(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,gamma,f)
          case (NWXCP_C_PBE)
            param(1) = 0.066724550603149d0
            call nwxc_c_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_PBESOL)
            param(1) = 0.046d0
            call nwxc_c_pbe(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                      rho,gamma,f)
          case (NWXCP_C_PKZB)
            call nwxc_c_pkzb99(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,tau,f)
          case (NWXCP_C_SPBE)
            call nwxc_c_spbe(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_C_SOGGA11)
            param(1)  =  5.00000d-01
            param(2)  = -4.62334D+00
            param(3)  =  8.00410D+00
            param(4)  = -130.226D+00
            param(5)  =  38.2685D+00
            param(6)  =  69.5599D+00
            param(7)  =  5.00000d-01
            param(8)  =  3.62334D+00
            param(9)  =  9.36393D+00
            param(10) =  34.5114D+00
            param(11) = -18.5684D+00
            param(12) = -0.16519D+00
            call nwxc_c_sogga(param,nwxc_rho_tol,ipol,nq,
     +                        nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_SOGGA11_X)
            param(1)  =  5.00000d-01
            param(2)  =  7.82439d+01
            param(3)  =  2.57211d+01
            param(4)  = -1.38830d+01
            param(5)  = -9.87375d+00
            param(6)  = -1.41357d+01
            param(7)  =  5.00000d-01
            param(8)  = -7.92439d+01
            param(9)  =  1.63725d+01
            param(10) =  2.08129d+00
            param(11) =  7.50769d+00
            param(12) = -1.01861d+01
            call nwxc_c_sogga(param,nwxc_rho_tol,ipol,nq,
     +                        nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_TPSS)
            call nwxc_c_tpss03(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,tau,f)
          case (NWXCP_C_VS98)
            param(1)  =  7.035010d-01
            param(2)  =  7.694574d-03
            param(3)  =  5.152765d-02
            param(4)  =  3.394308d-05
            param(5)  = -1.269420d-03
            param(6)  =  1.296118d-03
            param(7)  =  3.270912d-01
            param(8)  = -3.228915d-02
            param(9)  = -2.942406d-02
            param(10) =  2.134222d-03
            param(11) = -5.451559d-03
            param(12) =  1.577575d-02
            call nwxc_c_vs98(param,nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f)
          case (NWXCP_C_VWN1)
            call nwxc_c_vwn1(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,f)
          case (NWXCP_C_VWN1_RPA)
            call nwxc_c_vwn1_rpa(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,f)
          case (NWXCP_C_VWN2)
            call nwxc_c_vwn2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,f)
          case (NWXCP_C_VWN3)
            call nwxc_c_vwn3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,f)
          case (NWXCP_C_VWN4)
            call nwxc_c_vwn4(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,f)
          case (NWXCP_C_VWN5)
            call nwxc_c_vwn5(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,f)
          case (NWXCP_BOP)
            param(1) = 2.3670d0
            call nwxc_c_op(nwxc_k_becke88,param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_KT1)
            call nwxc_xc_kt1(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_PBEOP)
            param(1) = 2.3789d0
            call nwxc_c_op(nwxc_k_pbe96,param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_SOP)
            param(1) = 2.5654d0
            call nwxc_c_op(nwxc_k_dirac,param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
          case default
            if (nwxc_oroot) then
              write(*,*)"nwxc_eval_df: invalid functional",nwxc_ids(it)
            endif
            call nwxc_printP()
            call errquit("nwxc_eval_df: invalid functional",
     +                    nwxc_ids(it),0)
        end select
      enddo
C
      end
C>
C> @}
C>
C> \ingroup nwad_api
C> @{
C>
C> \brief The functional, 1st and 2nd order partial derivative evaluation
C>
C> In actual fact this routine only sets up the memory and then calls the
C> driver routine to drive the actual functional evaluation. The lack of 
C> a Fortran90 typecast drives this code structure.
C>
      subroutine nwxc_eval_df2(ipol,nq,rho,gamma,tau,f,dfdr,dfdr2,
     +                         dfdg,dfdg2,dfdt,dfdt2)
      use nwad2
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
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
      double precision dfdr2(nq,NCOL_AMAT2) !< [Output] 2nd derivative wrt rho
      double precision dfdg2(nq,NCOL_CMAT2) !< [Output] 2nd derivative wrt gamma
      double precision dfdt2(nq,NCOL_MMAT2) !< [Output] 2nd derivative wrt tau
c
      type(nwad_dble) :: memory_test_nwad(2)
      double precision   memory_test_dble(2)
      integer length
      integer l_rho, k_rho !< Index for density
      integer l_gam, k_gam !< Index for gradient
      integer l_tau, k_tau !< Index for kinetic energy density
      integer l_fun, k_fun !< Index for functional
      integer l_scr, k_scr !< Index for scratch
      integer nvar         !< The number of variables in the functional
                           !< LDA: nvar=ipol, GGA: nvar=2*ipol+1, 
                           !< MGGA: nvar=3*ipol+1
      integer ndrv         !< The number of partial derivatives
      integer ipl          !< Hack version of ipol
C
      logical nwxc_is_gga  !< Is the functional a GGA
      logical nwxc_is_mgga !< Is the functional a meta-GGA
c
c     Work out nvar for this functional
c
      ipl = 2
      nvar = 2
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        nvar = nvar + 3
        if (nwxc_is_mgga()) then
          nvar = nvar + 2
        endif
      endif
      ndrv = nvar*(nvar+1)/2
c
c     Work out how many double precision floating point words are needed
c     to represent a single nwad_dble
c
      length = (loc(memory_test_dble(2))-loc(memory_test_dble(1)))
      length = (loc(memory_test_nwad(2))-loc(memory_test_nwad(1))
     +          +length-1)/length
      if(.not.ma_push_get(mt_dbl,length*nq*ipl,"nwad rho",l_rho,k_rho))
     +  call errquit("nwxc_eval_df2: cannot allocate nwad rho",
     +    ma_sizeof(mt_dbl,length*nq*ipl,mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*(ipl+1),"nwad gamma",
     +                    l_gam,k_gam))
     +  call errquit("nwxc_eval_df2: cannot allocate nwad gamma",
     +    ma_sizeof(mt_dbl,length*nq*(ipl+1),mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*ipl,"nwad tau",l_tau,k_tau))
     +  call errquit("nwxc_eval_df2: cannot allocate nwad tau",
     +    ma_sizeof(mt_dbl,length*nq*ipl,mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*ndrv,"nwad fun",l_fun,k_fun))
     +  call errquit("nwxc_eval_df2: cannot allocate nwad fun",
     +    ma_sizeof(mt_dbl,length*nq*ndrv,mt_byte),MA_ERR)
c
      call nwxc_eval_df2_driver(ipol,nq,length,rho,gamma,tau,
     +     dbl_mb(k_rho),dbl_mb(k_gam),dbl_mb(k_tau),dbl_mb(k_fun),
     +     nvar,ndrv,f,dfdr,dfdr2,dfdg,dfdg2,dfdt,dfdt2)
c
      if(.not.ma_pop_stack(l_fun)) call errquit(
     +  "nwxc_eval_df2: cannot deallocate nwad fun",0,MA_ERR)
      if(.not.ma_pop_stack(l_tau)) call errquit(
     +  "nwxc_eval_df2: cannot deallocate nwad tau",0,MA_ERR)
      if(.not.ma_pop_stack(l_gam)) call errquit(
     +  "nwxc_eval_df2: cannot deallocate nwad gam",0,MA_ERR)
      if(.not.ma_pop_stack(l_rho)) call errquit(
     +  "nwxc_eval_df2: cannot deallocate nwad rho",0,MA_ERR)
      end
C>
C> @}
C>
C> \ingroup nwad_priv
C> @{
C>
C> \brief Driver routine for the functional, 1st and 2nd derivative evaluation
C>
C> This driver routine initializes the active variables and loops over the
C> functional evaluations. Afterwards it interpolates the partial derivatives
C> and returns the values of the functional and its 1st order partial
C> derivatives. Obviously, interpolating the 1st order partial derivatives is
C> particularly simply, but this scheme applies to all orders of derivatives.
C> Due to different data types being required the memory for this routine
C> needs to be allocated in the calling routine. Calling this subroutine
C> effectively implements a typecast.
C>
C> An issue that comes into play at the 2nd derivative level relates to the
C> cross terms in closed shell calculations. In closed shell calculations
C> we exploit that the \f$\alpha\f$ and \f$\beta\f$ spin densities have the
C> same values. As a result we know that for every partial derivative there
C> is another partial derivative with the \f$\alpha\f$ and \f$\beta\f$ labels
C> interchanged that has the same value. I.e. 
C> \f$\frac{\partial^2f}{\partial\rho_\alpha^2} = 
C> \frac{\partial^2f}{\partial\rho_\beta^2}\f$. However, there are also cross
C> terms such as \f$\frac{\partial^2f}{\partial\rho_\alpha\partial\rho_\beta}\f$
C> that are new at this level of differentiation and that cannot be expressed
C> in terms of anything else. Hence for closed shell systems we need to 
C> construct the same number of partial derivatives as for the open shell 
C> systems, but we can exploit the additional symmetry to reduce the number
C> of partial derivatives we actually have to calculate. The partial derivatives
C> that are the same due to the spin interchange symmetry can simply be copied.
C> The loop structure has been adapted for the closed shell case to exploit
C> this.
C>
      subroutine nwxc_eval_df2_driver(ipol,nq,ndbl,rho,gamma,tau,
     +           nwad_rho,nwad_gam,nwad_tau,nwad_f,nvar,ndrv,
     +           f,dfdr,dfdr2,dfdg,dfdg2,dfdt,dfdt2)
      use nwad2
      implicit none
#include "errquit.fh"
#include "nwxc_param.fh"
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of points
      integer ndbl !< [Input] The number doubles per nwad_dble
      integer nvar !< [Input] The number of independent variables
      integer ndrv !< [Input] The number of 2nd order partial derivatives
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
      double precision dfdr2(nq,NCOL_AMAT2) !< [Output] 2nd derivative wrt rho
      double precision dfdg2(nq,NCOL_CMAT2) !< [Output] 2nd derivative wrt gamma
      double precision dfdt2(nq,NCOL_MMAT2) !< [Output] 2nd derivative wrt tau
c
      type(nwad_dble) :: nwad_rho(nq,2)
      type(nwad_dble) :: nwad_gam(nq,3)
      type(nwad_dble) :: nwad_tau(nq,2)
      type(nwad_dble) :: nwad_f(nq,ndrv)
      type(nwad_dble) :: nwad_scr(nq)
C
      logical nwxc_is_gga  !< Is the functional a GGA
      logical nwxc_is_mgga !< Is the functional a meta-GGA
c
      integer iq, ip, id, ir
      integer ip2, ir2, ipr, id2
c
      do ip = 1, 2
        if (ipol.eq.1) then
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
          enddo
        else  ! ipol.eq.1
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(rho(iq,ip))
          enddo
        endif ! ipol.eq.1
      enddo
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        if (ipol.eq.1) then
          do ip = 1, 3
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(0.25d0*gamma(iq,G_TT))
            enddo
          enddo
        else  ! ipol.eq.1
          do ip = 1, 3
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(gamma(iq,ip))
            enddo
          enddo
        endif ! ipol.eq.1
        if (nwxc_is_mgga()) then
          do ip = 1, 2
            if (ipol.eq.1) then
              do iq = 1, nq
                nwad_tau(iq,ip) = inactive(0.5d0*tau(iq,T_T))
              enddo
            else  ! ipol.eq.1
              do iq = 1, nq
                nwad_tau(iq,ip) = inactive(1.0d0*tau(iq,ip))
              enddo
            endif ! ipol.eq.1
          enddo
        endif
      endif
c
c     Do closed shell and open shell separately
c
      if (ipol.eq.1) then
c
c       Do closed shell
c
c       First do d2f/dx2
c
        do ip = 1, nvar
          id = ip*(ip+1)/2
          if (ip.eq.1) then
            do iq = 1, nq
              nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T),2)
            enddo
          elseif (ip.eq.2) then
c           do nothing, we copy the rho alpha-alpha result
            id2 = 1*(1+1)/2
          elseif (ip.le.4) then
            do iq = 1, nq
              nwad_gam(iq,ip-2) = active(0.25d0*gamma(iq,G_TT),2)
            enddo
          elseif (ip.eq.5) then
c           do nothing, we copy the gamma alpha-alpha result
            id2 = 3*(3+1)/2
          elseif (ip.eq.6) then
            do iq = 1, nq
              nwad_tau(iq,ip-5) = active(0.5d0*tau(iq,T_T),2)
            enddo
          elseif (ip.eq.7) then
c           do nothing, we copy the tau alpha-alpha result
            id2 = 6*(6+1)/2
          endif
          if (ip.eq.1.or.ip.eq.3.or.ip.eq.4.or.ip.eq.6) then
            call nwxc_eval_df2_doit(2,nq,nwad_rho,nwad_gam,
     +           nwad_tau,nwad_f(1,id),dfdr,dfdr2,dfdg,dfdg2,
     +           dfdt,dfdt2)
          else
            call dcopy(nq*ndbl,nwad_f(1,id2),1,nwad_f(1,id),1)
          endif
          if (ip.eq.1) then
            do iq = 1, nq
              nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
            enddo
          elseif (ip.eq.2) then
c           do nothing
          elseif (ip.le.4) then
            do iq = 1, nq
              nwad_gam(iq,ip-2) = inactive(0.25d0*gamma(iq,G_TT))
            enddo
          elseif (ip.eq.5) then
c           do nothing
          elseif (ip.eq.6) then
            do iq = 1, nq
              nwad_tau(iq,ip-5) = inactive(0.5d0*tau(iq,T_T))
            enddo
          elseif (ip.eq.7) then
c           do nothing
          endif
        enddo
c
c       Second do d2f/dxdy
c
        do ip = 1, nvar
          do ir = 1, ip-1
            id = ip*(ip-1)/2+ir
            if (ip.le.2) then
              do iq = 1, nq
                nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T))
              enddo
            elseif (ip.le.5) then
              do iq = 1, nq
                nwad_gam(iq,ip-2) = active(0.25d0*gamma(iq,G_TT))
              enddo
            elseif (ip.le.7) then
              do iq = 1, nq
                nwad_tau(iq,ip-5) = active(0.5d0*tau(iq,T_T))
              enddo
            endif
            if (ir.le.2) then
              do iq = 1, nq
                nwad_rho(iq,ir) = active(0.5d0*rho(iq,R_T))
              enddo
            elseif (ir.le.5) then
              do iq = 1, nq
                nwad_gam(iq,ir-2) = active(0.25d0*gamma(iq,G_TT))
              enddo
            elseif (ir.le.7) then
              do iq = 1, nq
                nwad_tau(iq,ir-5) = active(0.5d0*tau(iq,T_T))
              enddo
            endif
            if ((ip.eq.2.and.ir.eq.1).or.(ip.eq.3.and.ir.le.2).or.
     +          (ip.eq.4.and.(ir.eq.1.or.ir.eq.3)).or.
     +          (ip.eq.5.and.ir.eq.3).or.(ip.eq.6).or.
     +          (ip.eq.7.and.ir.eq.6)) then
              call nwxc_eval_df2_doit(2,nq,nwad_rho,nwad_gam,
     +             nwad_tau,nwad_f(1,id),dfdr,dfdr2,dfdg,dfdg2,
     +             dfdt,dfdt2)
            else
              if (ip.eq.4.and.ir.eq.2) then
c               d2_ra_gab -> d2_rb_gab
                id2 = 4*(4-1)/2+1
              elseif (ip.eq.5.and.ir.eq.1) then
c               d2_rb_gaa -> d2_ra_gbb
                id2 = 3*(3-1)/2+2
              elseif (ip.eq.5.and.ir.eq.2) then
c               d2_ra_gaa -> d2_rb_gbb
                id2 = 3*(3-1)/2+1
              elseif (ip.eq.5.and.ir.eq.4) then
c               d2_gaa_gab -> d2_gab_gbb
                id2 = 4*(4-1)/2+3
              elseif (ip.eq.7.and.ir.eq.1) then
c               d2_rb_ta -> d2_ra_tb
                id2 = 6*(6-1)/2+2
              elseif (ip.eq.7.and.ir.eq.2) then
c               d2_ra_ta -> d2_rb_tb
                id2 = 6*(6-1)/2+1
              elseif (ip.eq.7.and.ir.eq.3) then
c               d2_gbb_ta -> d2_gaa_tb
                id2 = 6*(6-1)/2+5
              elseif (ip.eq.7.and.ir.eq.4) then
c               d2_gab_ta -> d2_gab_tb
                id2 = 6*(6-1)/2+4
              elseif (ip.eq.7.and.ir.eq.5) then
c               d2_gaa_ta -> d2_gbb_tb
                id2 = 6*(6-1)/2+3
              endif
              call dcopy(nq*ndbl,nwad_f(1,id2),1,nwad_f(1,id),1)
            endif
            if (ip.le.2) then
              do iq = 1, nq
                nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
              enddo
            elseif (ip.le.5) then
              do iq = 1, nq
                nwad_gam(iq,ip-2) = inactive(0.25d0*gamma(iq,G_TT))
              enddo
            elseif (ip.le.7) then
              do iq = 1, nq
                nwad_tau(iq,ip-5) = inactive(0.5d0*tau(iq,T_T))
              enddo
            endif
            if (ir.le.2) then
              do iq = 1, nq
                nwad_rho(iq,ir) = inactive(0.5d0*rho(iq,R_T))
              enddo
            elseif (ir.le.5) then
              do iq = 1, nq
                nwad_gam(iq,ir-2) = inactive(0.25d0*gamma(iq,G_TT))
              enddo
            elseif (ir.le.7) then
              do iq = 1, nq
                nwad_tau(iq,ir-5) = inactive(0.5d0*tau(iq,T_T))
              enddo
            endif
          enddo ! ir
        enddo ! ip
c
c       Done with all functional evaluations, now interpolate the results
c
c       First get the functional itself
c
c       do iq = 1, nq
c         f(iq) = value(nwad_f(iq,1))
c       enddo
c
c       Interpolate the first order partial derivatives
c
c       do ip = 1, nvar
c         id = ip*(ip+1)/2
c         if (ip.eq.1) then
c           do iq = 1, nq
c             dfdr(iq,ip) = inter_d2_dx(nwad_f(iq,id))
c           enddo
c         elseif (ip.eq.2) then
c         elseif (ip.le.4) then
c           do iq = 1, nq
c             dfdg(iq,ip-2) = inter_d2_dx(nwad_f(iq,id))
c           enddo
c         elseif (ip.eq.6) then
c           do iq = 1, nq
c             dfdt(iq,ip-5) = inter_d2_dx(nwad_f(iq,id))
c           enddo
c         endif
c       enddo
c
c       Interpolate the second order partial derivatives
c
c       First do d2fdx2
c
c       do ip = 1, nvar
c         id = ip*(ip+1)/2
c         if (ip.eq.1) then
c           do iq = 1, nq
c             dfdr2(iq,D2_RA_RA) = inter_d2_dx2(nwad_f(iq,id))
c           enddo
c         elseif (ip.eq.2) then
c         elseif (ip.eq.3) then
c           do iq = 1, nq
c             dfdg2(iq,D2_GAA_GAA) = inter_d2_dx2(nwad_f(iq,id))
c           enddo
c         elseif (ip.eq.4) then
c           do iq = 1, nq
c             dfdg2(iq,D2_GAB_GAB) = inter_d2_dx2(nwad_f(iq,id))
c           enddo
c         elseif (ip.eq.6) then
c           do iq = 1, nq
c             dfdt2(iq,D2_TA_TA) = inter_d2_dx2(nwad_f(iq,id))
c           enddo
c         endif
c       enddo
c
c       Second do d2fdxy
c
c       do ip = 1, nvar
c         do ir = 1, ip-1
c           ip2 = ip*(ip+1)/2
c           ir2 = ir*(ir+1)/2
c           ipr = ip*(ip-1)/2+ir
c           id2 = 0
c           if (ip.eq.2.and.ir.eq.1) then
c             id2 = D2_RA_RB
c           elseif (ip.eq.3.and.ir.eq.1) then
c             id2 = D2_RA_GAA
c           elseif (ip.eq.3.and.ir.eq.2) then
c             id2 = D2_RB_GAA
c           elseif (ip.eq.4.and.ir.eq.1) then
c             id2 = D2_RA_GAB
c           elseif (ip.eq.4.and.ir.eq.3) then
c             id2 = D2_GAA_GAB
c           elseif (ip.eq.5.and.ir.eq.3) then
c             id2 = D2_GAA_GBB
c           elseif (ip.eq.6.and.ir.eq.1) then
c             id2 = D2_RA_TA
c           elseif (ip.eq.6.and.ir.eq.2) then
c             id2 = D2_RB_TA
c           elseif (ip.eq.6.and.ir.eq.3) then
c             id2 = D2_GAA_TA
c           elseif (ip.eq.6.and.ir.eq.4) then
c             id2 = D2_GAB_TA
c           elseif (ip.eq.6.and.ir.eq.5) then
c             id2 = D2_GBB_TA
c           elseif (ip.eq.7.and.ir.eq.6) then
c             id2 = D2_TA_TB
c           endif
c           if (id2.ne.0) then
c             if (ip.le.2) then
c               do iq = 1, nq
c                 dfdr2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
c    +                            nwad_f(iq,ipr),nwad_f(iq,ir2))
c               enddo
c             elseif (ip.le.5) then
c               do iq = 1, nq
c                 dfdg2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
c    +                            nwad_f(iq,ipr),nwad_f(iq,ir2))
c               enddo
c             else
c               do iq = 1, nq
c                 dfdt2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
c    +                            nwad_f(iq,ipr),nwad_f(iq,ir2))
c               enddo
c             endif
c           endif ! id2.ne.0
c         enddo ! ir
c       enddo ! ip
c
      elseif (ipol.eq.2) then
c
c        Do open shell
c
c        First do d2f/dx2
c
         do ip = 1, nvar
           id = ip*(ip+1)/2
           if (ip.le.2) then
             do iq = 1, nq
               nwad_rho(iq,ip) = active(rho(iq,ip),2)
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               nwad_gam(iq,ip-2) = active(gamma(iq,ip-2),2)
             enddo
           else
             do iq = 1, nq
               nwad_tau(iq,ip-5) = active(1.0d0*tau(iq,ip-5),2)
             enddo
           endif
           call nwxc_eval_df2_doit(ipol,nq,nwad_rho,nwad_gam,nwad_tau,
     +                             nwad_f(1,id),dfdr,dfdr2,dfdg,dfdg2,
     +                             dfdt,dfdt2)
           if (ip.le.2) then
             do iq = 1, nq
               nwad_rho(iq,ip) = inactive(rho(iq,ip))
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               nwad_gam(iq,ip-2) = inactive(gamma(iq,ip-2))
             enddo
           else
             do iq = 1, nq
               nwad_tau(iq,ip-5) = inactive(1.0d0*tau(iq,ip-5))
             enddo
           endif
         enddo
c
c        Second do d2f/dxdy
c
         do ip = 1, nvar
           do ir = 1, ip-1
             id = ip*(ip-1)/2+ir
             if (ip.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = active(rho(iq,ip))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ip-2) = active(gamma(iq,ip-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ip-5) = active(1.0d0*tau(iq,ip-5))
               enddo
             endif
             if (ir.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = active(rho(iq,ir))
               enddo
             elseif (ir.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ir-2) = active(gamma(iq,ir-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ir-5) = active(1.0d0*tau(iq,ir-5))
               enddo
             endif
             call nwxc_eval_df2_doit(ipol,nq,nwad_rho,nwad_gam,nwad_tau,
     +                               nwad_f(1,id),dfdr,dfdr2,dfdg,dfdg2,
     +                               dfdt,dfdt2)
             if (ip.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = inactive(rho(iq,ip))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ip-2) = inactive(gamma(iq,ip-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ip-5) = inactive(1.0d0*tau(iq,ip-5))
               enddo
             endif
             if (ir.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = inactive(rho(iq,ir))
               enddo
             elseif (ir.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ir-2) = inactive(gamma(iq,ir-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ir-5) = inactive(1.0d0*tau(iq,ir-5))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Done with all functional evaluations, now interpolate the results
c
c        First get the functional itself
c
c        do iq = 1, nq
c          f(iq) = value(nwad_f(iq,1))
c        enddo
c
c        Interpolate the first order partial derivatives
c
c        do ip = 1, nvar
c          id = ip*(ip+1)/2
c          if (ip.le.2) then
c            do iq = 1, nq
c              dfdr(iq,ip) = inter_d2_dx(nwad_f(iq,id))
c            enddo
c          elseif (ip.le.5) then
c            do iq = 1, nq
c              dfdg(iq,ip-2) = inter_d2_dx(nwad_f(iq,id))
c            enddo
c          else
c            do iq = 1, nq
c              dfdt(iq,ip-5) = inter_d2_dx(nwad_f(iq,id))
c            enddo
c          endif
c        enddo
c
c        Interpolate the second order partial derivatives
c
c        First do d2fdx2
c
c        do ip = 1, nvar
c          id = ip*(ip+1)/2
c          id2 = D2_RA_RA
c          if (ip.eq.2) then
c            id2 = D2_RB_RB
c          elseif (ip.eq.3) then
c            id2 = D2_GAA_GAA
c          elseif (ip.eq.4) then
c            id2 = D2_GAB_GAB
c          elseif (ip.eq.5) then
c            id2 = D2_GBB_GBB
c          elseif (ip.eq.6) then
c            id2 = D2_TA_TA
c          elseif (ip.eq.7) then
c            id2 = D2_TB_TB
c          endif
c          if (ip.le.2) then
c            do iq = 1, nq
c              dfdr2(iq,id2) = inter_d2_dx2(nwad_f(iq,id))
c            enddo
c          elseif (ip.le.5) then
c            do iq = 1, nq
c              dfdg2(iq,id2) = inter_d2_dx2(nwad_f(iq,id))
c            enddo
c          else
c            do iq = 1, nq
c              dfdt2(iq,id2) = inter_d2_dx2(nwad_f(iq,id))
c            enddo
c          endif
c        enddo
c
c        Second do d2fdxy
c
c        do ip = 1, nvar
c          do ir = 1, ip-1
c            ip2 = ip*(ip+1)/2
c            ir2 = ir*(ir+1)/2
c            ipr = ip*(ip-1)/2+ir
c            id2 = 0
c            if (ip.eq.2.and.ir.eq.1) then
c              id2 = D2_RA_RB
c            elseif (ip.eq.3.and.ir.eq.1) then
c              id2 = D2_RA_GAA
c            elseif (ip.eq.3.and.ir.eq.2) then
c              id2 = D2_RB_GAA
c            elseif (ip.eq.4.and.ir.eq.1) then
c              id2 = D2_RA_GAB
c            elseif (ip.eq.4.and.ir.eq.2) then
c              id2 = D2_RB_GAB
c            elseif (ip.eq.4.and.ir.eq.3) then
c              id2 = D2_GAA_GAB
c            elseif (ip.eq.5.and.ir.eq.1) then
c              id2 = D2_RA_GBB
c            elseif (ip.eq.5.and.ir.eq.2) then
c              id2 = D2_RB_GBB
c            elseif (ip.eq.5.and.ir.eq.3) then
c              id2 = D2_GAA_GBB
c            elseif (ip.eq.5.and.ir.eq.4) then
c              id2 = D2_GAB_GBB
c            elseif (ip.eq.6.and.ir.eq.1) then
c              id2 = D2_RA_TA
c            elseif (ip.eq.6.and.ir.eq.2) then
c              id2 = D2_RB_TA
c            elseif (ip.eq.6.and.ir.eq.3) then
c              id2 = D2_GAA_TA
c            elseif (ip.eq.6.and.ir.eq.4) then
c              id2 = D2_GAB_TA
c            elseif (ip.eq.6.and.ir.eq.5) then
c              id2 = D2_GBB_TA
c            elseif (ip.eq.7.and.ir.eq.1) then
c              id2 = D2_RA_TB
c            elseif (ip.eq.7.and.ir.eq.2) then
c              id2 = D2_RB_TB
c            elseif (ip.eq.7.and.ir.eq.3) then
c              id2 = D2_GAA_TB
c            elseif (ip.eq.7.and.ir.eq.4) then
c              id2 = D2_GAB_TB
c            elseif (ip.eq.7.and.ir.eq.5) then
c              id2 = D2_GBB_TB
c            elseif (ip.eq.7.and.ir.eq.6) then
c              id2 = D2_TA_TB
c            endif
c            if (ip.le.2) then
c              do iq = 1, nq
c                dfdr2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
c    +                           nwad_f(iq,ipr),nwad_f(iq,ir2))
c              enddo
c            elseif (ip.le.5) then
c              do iq = 1, nq
c                dfdg2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
c    +                           nwad_f(iq,ipr),nwad_f(iq,ir2))
c              enddo
c            else
c              do iq = 1, nq
c                dfdt2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
c    +                           nwad_f(iq,ipr),nwad_f(iq,ir2))
c              enddo
c            endif
c          enddo ! ir
c        enddo ! ip
c
      else
        call errquit("nwxc_eval_df2_driver: insane ipol",ipol,UERR)
      endif
c
c
c     Done with all functional evaluations, now interpolate the results
c
c     First get the functional itself
c
      do iq = 1, nq
        f(iq) = value(nwad_f(iq,1))
      enddo
c
c     Interpolate the first order partial derivatives
c
      do ip = 1, nvar
        id = ip*(ip+1)/2
        if (ip.le.2) then
          do iq = 1, nq
            dfdr(iq,ip) = inter_d2_dx(nwad_f(iq,id))
          enddo
        elseif (ip.le.5) then
          do iq = 1, nq
            dfdg(iq,ip-2) = inter_d2_dx(nwad_f(iq,id))
          enddo
        else
          do iq = 1, nq
            dfdt(iq,ip-5) = inter_d2_dx(nwad_f(iq,id))
          enddo
        endif
      enddo
c
c     Interpolate the second order partial derivatives
c
c     First do d2fdx2
c
      do ip = 1, nvar
        id = ip*(ip+1)/2
        id2 = D2_RA_RA
        if (ip.eq.2) then
          id2 = D2_RB_RB
        elseif (ip.eq.3) then
          id2 = D2_GAA_GAA
        elseif (ip.eq.4) then
          id2 = D2_GAB_GAB
        elseif (ip.eq.5) then
          id2 = D2_GBB_GBB
        elseif (ip.eq.6) then
          id2 = D2_TA_TA
        elseif (ip.eq.7) then
          id2 = D2_TB_TB
        endif
        if (ip.le.2) then
          do iq = 1, nq
            dfdr2(iq,id2) = inter_d2_dx2(nwad_f(iq,id))
          enddo
        elseif (ip.le.5) then
          do iq = 1, nq
            dfdg2(iq,id2) = inter_d2_dx2(nwad_f(iq,id))
          enddo
        else
          do iq = 1, nq
            dfdt2(iq,id2) = inter_d2_dx2(nwad_f(iq,id))
          enddo
        endif
      enddo
c
c     Second do d2fdxy
c
      do ip = 1, nvar
        do ir = 1, ip-1
          ip2 = ip*(ip+1)/2
          ir2 = ir*(ir+1)/2
          ipr = ip*(ip-1)/2+ir
          id2 = 0
          if (ip.eq.2.and.ir.eq.1) then
            id2 = D2_RA_RB
          elseif (ip.eq.3.and.ir.eq.1) then
            id2 = D2_RA_GAA
          elseif (ip.eq.3.and.ir.eq.2) then
            id2 = D2_RB_GAA
          elseif (ip.eq.4.and.ir.eq.1) then
            id2 = D2_RA_GAB
          elseif (ip.eq.4.and.ir.eq.2) then
            id2 = D2_RB_GAB
          elseif (ip.eq.4.and.ir.eq.3) then
            id2 = D2_GAA_GAB
          elseif (ip.eq.5.and.ir.eq.1) then
            id2 = D2_RA_GBB
          elseif (ip.eq.5.and.ir.eq.2) then
            id2 = D2_RB_GBB
          elseif (ip.eq.5.and.ir.eq.3) then
            id2 = D2_GAA_GBB
          elseif (ip.eq.5.and.ir.eq.4) then
            id2 = D2_GAB_GBB
          elseif (ip.eq.6.and.ir.eq.1) then
            id2 = D2_RA_TA
          elseif (ip.eq.6.and.ir.eq.2) then
            id2 = D2_RB_TA
          elseif (ip.eq.6.and.ir.eq.3) then
            id2 = D2_GAA_TA
          elseif (ip.eq.6.and.ir.eq.4) then
            id2 = D2_GAB_TA
          elseif (ip.eq.6.and.ir.eq.5) then
            id2 = D2_GBB_TA
          elseif (ip.eq.7.and.ir.eq.1) then
            id2 = D2_RA_TB
          elseif (ip.eq.7.and.ir.eq.2) then
            id2 = D2_RB_TB
          elseif (ip.eq.7.and.ir.eq.3) then
            id2 = D2_GAA_TB
          elseif (ip.eq.7.and.ir.eq.4) then
            id2 = D2_GAB_TB
          elseif (ip.eq.7.and.ir.eq.5) then
            id2 = D2_GBB_TB
          elseif (ip.eq.7.and.ir.eq.6) then
            id2 = D2_TA_TB
          endif
          if (ip.le.2) then
            do iq = 1, nq
              dfdr2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
     +                        nwad_f(iq,ipr),nwad_f(iq,ir2))
            enddo
          elseif (ip.le.5) then
            do iq = 1, nq
              dfdg2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
     +                        nwad_f(iq,ipr),nwad_f(iq,ir2))
            enddo
          else
            do iq = 1, nq
              dfdt2(iq,id2) = inter_d2_dxy(nwad_f(iq,ip2),
     +                        nwad_f(iq,ipr),nwad_f(iq,ir2))
            enddo
          endif
        enddo ! ir
      enddo ! ip
C      
      end
C>
C> \brief Evaluate the exchange-correlation energy and its 1st and 2nd
C> partial derivatives
C>
      subroutine nwxc_eval_df2_doit(ipol,nq,rho,gamma,tau,f,dfdr,dfdr2,
     +                         dfdg,dfdg2,dfdt,dfdt2)
      use nwad2
      implicit none
#include "nwxcP.fh"
#include "nwxc_param.fh"
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of points
C
      type(nwad_dble)::rho(nq,ipol)     !< [Input] Density
      type(nwad_dble)::gamma(nq,ipol+1) !< [Input] |Density gradient|^2
      type(nwad_dble)::tau(nq,ipol)     !< [Input] Kinetic energy 
                                        !< density
C
      type(nwad_dble)::f(nq)           !< [Output] Energy
      double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
      double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
      double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
c
      double precision dfdr2(nq,*)     !< [Output] 2nd derivative wrt
                                       !< rho
      double precision dfdg2(nq,*)     !< [Output] 2nd derivative wrt
                                       !< rho and gamma, and gamma
      double precision dfdt2(nq,*)     !< [Output] 2nd derivative wrt
                                       !< rho and tau, gamma and tau,
                                       !< and tau
C
      integer iq !< Counter over points
      integer ip !< Counter over spin channels
      integer it !< Counter of functional terms
      integer n  !< Counter
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
      integer nd2r(2) !< The number of 2nd partial derivatives wrt rho
                      !< as a function of ipol
      integer nd2g(2) !< The number of 2nd partial derivatives wrt
                      !< gamma as a function of ipol
      integer nd2t(2) !< The number of 2nd partial derivatives wrt tau
                      !< as a function of ipol
C
      external nwxc_k_dirac_d2
      external nwxc_k_becke88_d2
      external nwxc_k_pbe96_d2
c
      integer max_param
      parameter (max_param = 50)
      double precision param(max_param)
      double precision Axlsda
      parameter (Axlsda = -0.9305257363491d0 )
C
      data nd1r / D1_RA,  D1_RB  /
      data nd1g / D1_GAA, D1_GBB /
      data nd1t / D1_TA,  D1_TB  /
C
      data nd2r / D2_RA_RA,   D2_RB_RB  /
c      data nd2g / D2_GAA_GAA, D2_GBB_GBB /
      data nd2g / D2_GAA_GBB, D2_GBB_GBB /
      data nd2t / D2_TA_TA,   D2_TB_TB  /
C
      do iq = 1, nq
        f(iq) = 0.0d0
      enddo
c      do ip = 1, nd1r(ipol)
c        do iq = 1, nq
c          dfdr(iq,ip) = 0.0d0
c        enddo
c      enddo
c      do ip = 1, nd2r(ipol)
c        do iq = 1, nq
c          dfdr2(iq,ip) = 0.0d0
c        enddo
c      enddo
c      if (nwxc_is_gga().or.nwxc_is_mgga()) then
c        do ip = 1, nd1g(ipol)
c          do iq = 1, nq
c            dfdg(iq,ip) = 0.0d0
c          enddo
c        enddo
c        do ip = 1, nd2g(ipol)
c          do iq = 1, nq
c            dfdg2(iq,ip) = 0.0d0
c          enddo
c        enddo
c      endif
c      if (nwxc_is_mgga()) then
c        do ip = 1, nd1t(ipol)
c          do iq = 1, nq
c            dfdt(iq,ip) = 0.0d0
c          enddo
c        enddo
c        do ip = 1, nd2t(ipol)
c          do iq = 1, nq
c            dfdt2(iq,ip) = 0.0d0
c          enddo
c        enddo
c      endif
C
      do it = 1, nwxc_num_f
        select case (nwxc_ids(it))
          case (NWXCP_X_SLATER)
            call nwxc_x_dirac_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,f)
          case (NWXCP_X_B86b)
            call nwxc_x_b86b_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_B88)
            call nwxc_x_b88_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,f)
          case (NWXCP_X_B97)
            param(1) = 2.0d0
            param(2) = 0.80940d+00
            param(3) = 0.50730d+00
            param(4) = 0.74810d+00
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_B97_1)
            param(1) = 2.0d0
            param(2) = 0.789518d+00
            param(3) = 0.573805d+00
            param(4) = 0.660975d+00
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_B97_2)
            param(1) = 2.0d0
            param(2) = 0.827642D+00
            param(3) = 0.478400D-01
            param(4) = 0.176125D+01
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_B97_3)
            param(1) = 4.0d0
            param(2) = 7.334648D-01
            param(3) = 2.925270D-01
            param(4) = 3.338789D+00
            param(5) =-1.051158D+01
            param(6) = 1.060907D+01
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_B97_D)
            param(1) = 2.0d0
            param(2) = 1.086620d+0
            param(3) =-0.521270d+00
            param(4) = 3.254290d+00
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_B97_G)
            param(1) = 2.0d0
            param(2) = 1.1068d0
            param(3) =-0.8765d0
            param(4) = 4.2639d0
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_B98)
            param(1) = 2.0d0
            param(2) = 0.790194d00
            param(3) = 0.400271d00
            param(4) = 0.832857d00
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_BNL)
            param(1) = nwxc_cam_gamma
            call nwxc_x_bnl_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,f)
          case (NWXCP_X_CAMB88)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camb88_d2(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMLSD)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camlsd_d2(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,f)
          case (NWXCP_X_CAMPBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            call nwxc_x_campbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMREVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            call nwxc_x_campbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMRPBE)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camrpbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                             nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMS12G)
            param(1) = 1.03323556d0
            param(2) = 0.757d0
            param(3) = 0.00417251d0
            param(4) = 0.00115216d0
            param(5) = 0.00706184d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            call nwxc_x_cams12_d2(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMS12H)
            param(1) = 1.02149642d0
            param(2) = 0.757d0
            param(3) = 0.00825905d0
            param(4) = 0.00235804d0
            param(5) = 0.00654977d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            call nwxc_x_cams12_d2(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_DLDF)
             call nwxc_x_dldf_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +            rho,gamma,tau,f)
          case (NWXCP_X_FT97)
            call nwxc_x_ft97_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_GILL)
            call nwxc_x_gill_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_HCTH)
            param(1) = 4.0d0
            param(2) = 0.109320d+01
            param(3) =-0.744056d+00
            param(4) = 0.559920d+01
            param(5) =-0.678549d+01
            param(6) = 0.449357d+01
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH120)
            param(1) = 4.0d0
            param(2) = 1.09163d0
            param(3) =-0.74720d0
            param(4) = 5.07830d0
            param(5) =-4.10750d0
            param(6) = 1.17170d0
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH147)
            param(1) = 4.0d0
            param(2) = 1.09025d0
            param(3) =-0.79920d0
            param(4) = 5.57210d0
            param(5) =-5.86760d0
            param(6) = 3.04540d0
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH407)
            param(1) = 4.0d0
            param(2) = 1.08184d0
            param(3) =-0.5183d0
            param(4) = 3.4256d0
            param(5) =-2.6290d0
            param(6) = 2.2886d0
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH407P)
            param(1) = 4.0d0
            param(2) = 1.08018D0
            param(3) =-0.4117D0
            param(4) = 2.4368D0
            param(5) = 1.3890D0
            param(6) =-1.3529D0
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH_A)
            param(1) = 2.0d0
            param(2) = 0.109878d+01
            param(3) =-0.251173d+01
            param(4) = 0.156233d-01
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTHP14)
            param(1) = 4.0d0
            param(2) = 0.103161d+01
            param(3) =-0.360781d+00
            param(4) = 0.351994d+01
            param(5) =-0.495944d+01
            param(6) = 0.241165d+01
            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_M05)
             param( 1) =    0.08151d0
             param( 2) =   -0.43956d0
             param( 3) =   -3.22422d0
             param( 4) =    2.01819d0
             param( 5) =    8.79431d0
             param( 6) =   -0.00295d0
             param( 7) =    9.82029d0
             param( 8) =   -4.82351d0
             param( 9) =  -48.17574d0
             param(10) =    3.64802d0
             param(11) =   34.02248d0
             call nwxc_x_m05_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M05_2X)
             param( 1) =   -0.56833d0
             param( 2) =   -1.30057d0
             param( 3) =    5.50070d0
             param( 4) =    9.06402d0
             param( 5) =  -32.21075d0
             param( 6) =  -23.73298d0
             param( 7) =   70.22996d0
             param( 8) =   29.88614d0
             param( 9) =  -60.25778d0
             param(10) =  -13.22205d0
             param(11) =   15.23694d0
             call nwxc_x_m05_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06)
             param( 1) =  1.422057D-01*Axlsda
             param( 2) =  7.370319D-04*Axlsda
             param( 3) = -1.601373D-02*Axlsda
             param( 4) =  0.000000D+00
             param( 5) =  0.000000D+00
             param( 6) =  0.000000D+00
             param( 7) =  5.877943D-01
             param( 8) = -1.371776D-01
             param( 9) =  2.682367D-01
             param(10) = -2.515898D+00
             param(11) = -2.978892D+00
             param(12) =  8.710679D+00
             param(13) =  1.688195D+01
             param(14) = -4.489724D+00
             param(15) = -3.299983D+01
             param(16) = -1.449050D+01
             param(17) =  2.043747D+01
             param(18) =  1.256504D+01
             call nwxc_x_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06_HF)
             param( 1) = -1.179732D-01*Axlsda
             param( 2) = -2.500000D-03*Axlsda
             param( 3) = -1.180065D-02*Axlsda
             param( 4) =  0.000000D+00
             param( 5) =  0.000000D+00
             param( 6) =  0.000000D+00
             param( 7) =  1.179732D-01
             param( 8) = -1.066708D+00
             param( 9) = -1.462405D-01
             param(10) =  7.481848D+00
             param(11) =  3.776679D+00
             param(12) = -4.436118D+01
             param(13) = -1.830962D+01
             param(14) =  1.003903D+02
             param(15) =  3.864360D+01
             param(16) = -9.806018D+01
             param(17) = -2.557716D+01
             param(18) =  3.590404D+01
             call nwxc_x_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06_L)
             param( 1) =  6.012244D-01*Axlsda
             param( 2) =  4.748822D-03*Axlsda
             param( 3) = -8.635108D-03*Axlsda
             param( 4) = -9.308062D-06*Axlsda
             param( 5) =  4.482811D-05*Axlsda
             param( 6) =  0.000000D+00
             param( 7) =  3.987756D-01
             param( 8) =  2.548219D-01
             param( 9) =  3.923994D-01
             param(10) = -2.103655D+00
             param(11) = -6.302147D+00
             param(12) =  1.097615D+01
             param(13) =  3.097273D+01
             param(14) = -2.318489D+01
             param(15) = -5.673480D+01
             param(16) =  2.160364D+01
             param(17) =  3.421814D+01
             param(18) = -9.049762D+00
             call nwxc_x_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06_2X)
             param( 1) =  4.600000D-01
             param( 2) = -2.206052D-01
             param( 3) = -9.431788D-02
             param( 4) =  2.164494D+00
             param( 5) = -2.556466D+00
             param( 6) = -1.422133D+01
             param( 7) =  1.555044D+01
             param( 8) =  3.598078D+01
             param( 9) = -2.722754D+01
             param(10) = -3.924093D+01
             param(11) =  1.522808D+01
             param(12) =  1.522227D+01
             call nwxc_x_m06_2x_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M08_HX)
c            parameters A
             param(01) =  1.3340172D+00
             param(02) = -9.4751087D+00
             param(03) = -1.2541893D+01
             param(04) =  9.1369974D+00
             param(05) =  3.4717204D+01
             param(06) =  5.8831807D+01
             param(07) =  7.1369574D+01
             param(08) =  2.3312961D+01
             param(09) =  4.8314679D+00
             param(10) = -6.5044167D+00
             param(11) = -1.4058265D+01
             param(12) =  1.2880570D+01
c            parameters B
             param(13) = -8.5631823D-01
             param(14) =  9.2810354D+00
             param(15) =  1.2260749D+01
             param(16) = -5.5189665D+00
             param(17) = -3.5534989D+01
             param(18) = -8.2049996D+01
             param(19) = -6.8586558D+01
             param(20) =  3.6085694D+01
             param(21) = -9.3740983D+00
             param(22) = -5.9731688D+01
             param(23) =  1.6587868D+01
             param(24) =  1.3993203D+01
c            parameters C and D
             do n = 25, 48
               param(n) = 0.0d0
             enddo
             call nwxc_x_m08_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M08_SO)
c            parameters A
             param(01) = -3.4888428D-01
             param(02) = -5.8157416D+00
             param(03) =  3.7550810D+01
             param(04) =  6.3727406D+01
             param(05) = -5.3742313D+01
             param(06) = -9.8595529D+01
             param(07) =  1.6282216D+01
             param(08) =  1.7513468D+01
             param(09) = -6.7627553D+00
             param(10) =  1.1106658D+01
             param(11) =  1.5663545D+00
             param(12) =  8.7603470D+00
c            parameters B
             param(13) =  7.8098428D-01
             param(14) =  5.4538178D+00
             param(15) = -3.7853348D+01
             param(16) = -6.2295080D+01
             param(17) =  4.6713254D+01
             param(18) =  8.7321376D+01
             param(19) =  1.6053446D+01
             param(20) =  2.0126920D+01
             param(21) = -4.0343695D+01
             param(22) = -5.8577565D+01
             param(23) =  2.0890272D+01
             param(24) =  1.0946903D+01
c            parameters C and D
             do n = 25, 48
               param(n) = 0.0d0
             enddo
             call nwxc_x_m08_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M11)
c            parameters A
             param(01) = -0.18399900D+00
             param(02) = -1.39046703D+01
             param(03) =  1.18206837D+01
             param(04) =  3.10098465D+01
             param(05) = -5.19625696D+01
             param(06) =  1.55750312D+01
             param(07) = -6.94775730D+00
             param(08) = -1.58465014D+02
             param(09) = -1.48447565D+00
             param(10) =  5.51042124D+01
             param(11) = -1.34714184D+01
             param(12) =  0.000000D+00
c            parameters B
             param(13) =  0.75599900D+00
             param(14) =  1.37137944D+01
             param(15) = -1.27998304D+01
             param(16) = -2.93428814D+01
             param(17) =  5.91075674D+01
             param(18) = -2.27604866D+01
             param(19) = -1.02769340D+01
             param(20) =  1.64752731D+02
             param(21) =  1.85349258D+01
             param(22) = -5.56825639D+01
             param(23) =  7.47980859D+00
             param(24) =  0.000000D+00
c            parameters C and D
             do n = 25, 48
               param(n) = 0.0d0
             enddo
             call nwxc_x_m11_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M11_L)
c            parameters A
             param(01) =  8.121131D-01
             param(02) =  1.738124D+01
             param(03) =  1.154007D+00
             param(04) =  6.869556D+01
             param(05) =  1.016864D+02
             param(06) = -5.887467D+00
             param(07) =  4.517409D+01
             param(08) = -2.773149D+00
             param(09) = -2.617211D+01
             param(10) =  0.000000D+00
             param(11) =  0.000000D+00
             param(12) =  0.000000D+00
c            parameters B
             param(13) =  1.878869D-01
             param(14) = -1.653877D+01
             param(15) =  6.755753D-01
             param(16) = -7.567572D+01
             param(17) = -1.040272D+02
             param(18) =  1.831853D+01
             param(19) = -5.573352D+01
             param(20) = -3.520210D+00
             param(21) =  3.724276D+01
             param(22) =  0.000000D+00
             param(23) =  0.000000D+00
             param(24) =  0.000000D+00
c            parameters C
             param(25) = -4.386615D-01
             param(26) = -1.214016D+02
             param(27) = -1.393573D+02
             param(28) = -2.046649D+00
             param(29) =  2.804098D+01
             param(30) = -1.312258D+01
             param(31) = -6.361819D+00
             param(32) = -8.055758D-01
             param(33) =  3.736551D+00
             param(34) =  0.000000D+00
             param(35) =  0.000000D+00
             param(36) =  0.000000D+00
c            parameters D
             param(37) =  1.438662D+00
             param(38) =  1.209465D+02
             param(39) =  1.328252D+02
             param(40) =  1.296355D+01
             param(41) =  5.854866D+00
             param(42) = -3.378162D+00
             param(43) = -4.423393D+01
             param(44) =  6.844475D+00
             param(45) =  1.949541D+01
             param(46) =  0.000000D+00
             param(47) =  0.000000D+00
             param(48) =  0.000000D+00
             call nwxc_x_m11_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_X_MPW91)
            param(1) = 3.72d0
            param(2) = 0.00426D0
            call nwxc_x_pw91_d2(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_OPT)
             call nwxc_x_opt_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_PW6B95)
            param(1) = 0.00538d0
            param(2) = 1.7382d0
            param(3) = 3.8901d0
            call nwxc_x_pw6_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_PWB6K)
            param(1) = 0.00539d0
            param(2) = 1.7077d0
            param(3) = 4.0876d0
            call nwxc_x_pw6_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_PW91)
            param(1) = 4.0d0
            param(2) = 0.0042D0
            call nwxc_x_pw91_d2(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_PBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            call nwxc_x_pbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_PBESOL)
            param(1) = 0.8040d0
            param(2) = 10.0d0/81.0d0
            call nwxc_x_pbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_REVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            call nwxc_x_pbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_RPBE)
            call nwxc_x_rpbe_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_PKZB)
            call nwxc_x_pkzb99_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,tau,f)
          case (NWXCP_X_S12G)
            param(1) = 1.03842032d0
            param(2) = 0.757d0
            param(3) = 0.00403198d0
            param(4) = 0.00104596d0
            param(5) = 0.00594635d0
            call nwxc_x_s12_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_S12H)
            param(1) = 1.02543951d0
            param(2) = 0.757d0
            param(3) = 0.00761554d0
            param(4) = 0.00211063d0
            param(5) = 0.00604672d0
            call nwxc_x_s12_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SOGGA)
            param(1)  =  0.5d0
            param(2)  =  0.276d0
            param(3)  =  0.0d0
            param(4)  =  0.0d0
            param(5)  =  0.0d0
            param(6)  =  0.0d0
            param(7)  =  0.5d0
            param(8)  =  0.276d0
            param(9)  =  0.0d0
            param(10) =  0.0d0
            param(11) =  0.0d0
            param(12) =  0.0d0
            call nwxc_x_sogga_d2(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SOGGA11)
            param(1)  =  0.5d0
            param(2)  = -2.95535d0
            param(3)  =  15.7974d0
            param(4)  = -91.1804d0
            param(5)  =  96.2030d0
            param(6)  =  0.18683d0
            param(7)  =  0.50000d0
            param(8)  =  3.50743d0
            param(9)  = -12.9523d0
            param(10) =  49.7870d0
            param(11) = -33.2545d0
            param(12) = -11.1396d0
            call nwxc_x_sogga_d2(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SOGGA11_X)
            param(1)  =  2.99250d-01
            param(2)  =  3.21638d+00
            param(3)  = -3.55605d+00
            param(4)  =  7.65852d+00
            param(5)  = -1.12830d+01
            param(6)  =  5.25813d+00
            param(7)  =  2.99250d-01
            param(8)  = -2.88595d+00
            param(9)  =  3.23617d+00
            param(10) = -2.45393d+00
            param(11) = -3.75495d+00
            param(12) =  3.96613d+00
            call nwxc_x_sogga_d2(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_SSB_D)
            call nwxc_x_ssbD_1_d2(nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_TPSS)
            call nwxc_x_tpss03_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                            rho,gamma,tau,f)
          case (NWXCP_X_WPBE)
            param(1) = nwxc_cam_gamma
            call nwxc_x_wpbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_VS98)
             param(1) =  -9.800683d-01
             param(2) =  -3.556788d-03
             param(3) =   6.250326d-03
             param(4) =  -2.354518d-05
             param(5) =  -1.282732d-04
             param(6) =   3.574822d-04
             call nwxc_x_vs98_d2(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_B95)
            param(1) = 0.0031d0
            param(2) = 0.038d0
            call nwxc_c_b95_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_B97)
            param(1) = 2.0d0
            param(2) = 0.17370d+00
            param(3) = 0.94540d+00
            param(4) = 0.23487d+01
            param(5) = 0.74710d+00
            param(6) =-0.24868d+01
            param(7) =-0.45961d+01
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_B97_1)
            param(1) = 2.0d0
            param(2) = 0.820011d-01
            param(3) = 0.955689d+00
            param(4) = 0.271681d+01
            param(5) = 0.788552d+00
            param(6) =-0.287103d+01
            param(7) =-0.547869d+01
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_B97_2)
            param(1) = 2.0d0
            param(2) = 0.585808D+00
            param(3) = 0.999849D+00
            param(4) =-0.691682D+00
            param(5) = 0.140626D+01
            param(6) = 0.394796D+00
            param(7) =-0.744060D+01
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_B97_3)
            param(1)  = 4.0d0
            param(2)  = 5.623649D-01
            param(3)  = 1.133830D+00
            param(4)  =-1.322980D+00
            param(5)  =-2.811967D+00
            param(6)  = 6.359191D+00
            param(7)  = 7.431302D+00
            param(8)  =-7.464002D+00
            param(9)  =-1.969342D+00
            param(10) = 1.827082D+00
            param(11) =-1.174423D+01
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_B97_D)
            param(1) = 2.0d0
            param(2) = 0.22340d+00
            param(3) = 0.690410d+00
            param(4) =-1.562080d+00
            param(5) = 6.302700d00
            param(6) = 1.942930d+0
            param(7) =-14.97120d+00
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_B97_G)
            param(1) = 2.0d0
            param(2) = 0.4883d0
            param(3) = 0.7961d0
            param(4) =-2.117d0
            param(5) = 5.7060d0
            param(6) = 2.3235d0
            param(7) =-14.9820d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_B98)
            param(1) = 2.0d0
            param(2) =-0.120163d00
            param(3) = 0.934715d00
            param(4) = 2.82332d0
            param(5) = 1.14105d0
            param(6) =-2.59412d0
            param(7) =-5.33398d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_DLDF)
            call nwxc_c_dldf_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,tau,f)
          case (NWXCP_C_FT97)
            call nwxc_c_ft97_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f)
          case (NWXCP_C_HCTH)
            param(1)  = 4.0d0
            param(2)  = 0.222601d0
            param(3)  = 0.729974d0
            param(4)  =-3.38622d-002
            param(5)  = 3.352870d0
            param(6)  =-1.25170d-002
            param(7)  =-11.543d0
            param(8)  =-0.802496d0
            param(9)  = 8.085640d0
            param(10) = 1.553960d0
            param(11) =-4.478570d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_HCTH120)
            param(1)  = 4.0d0
            param(2)  = 0.48951d0
            param(3)  = 0.51473d0
            param(4)  =-0.26070d0
            param(5)  = 6.92980d0
            param(6)  = 0.43290d0
            param(7)  =-24.7070d0
            param(8)  =-1.99250d0
            param(9)  = 23.1100d0
            param(10) = 2.48530d0
            param(11) =-11.3230d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_HCTH147)
            param(1)  = 4.0d0
            param(2)  = 0.56258d0
            param(3)  = 0.54235d0
            param(4)  =-1.71000d-002
            param(5)  = 7.01460d0
            param(6)  =-1.30640d0
            param(7)  =-28.3820d0
            param(8)  = 1.05750d0
            param(9)  = 35.0330d0
            param(10) = 0.88540d0
            param(11) =-20.4280d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_HCTH407)
            param(1)  = 4.0d0
            param(2)  = 1.18777d0
            param(3)  = 0.58908d0
            param(4)  =-2.40290d0
            param(5)  = 4.42370d0
            param(6)  = 5.61740d0
            param(7)  =-19.2220d0
            param(8)  =-9.17920d0
            param(9)  = 42.5720d0
            param(10) = 6.24800d0
            param(11) =-42.0050d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_HCTH407P)
            param(1)  = 4.0d0
            param(2)  = 0.80302d0
            param(3)  = 0.73604d0
            param(4)  =-1.04790d0
            param(5)  = 3.02700d0
            param(6)  = 4.98070d0
            param(7)  =-10.0750d0
            param(8)  =-12.8900d0
            param(9)  = 20.6110d0
            param(10) = 9.64460d0
            param(11) =-29.4180d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_HCTH_A)
            param(1)  = 4.0d0
            param(2)  = 1.36823d-002
            param(3)  = 0.836897d0
            param(4)  = 0.268920d0
            param(5)  = 1.720510d0
            param(6)  =-0.550769d0
            param(7)  =-2.784980d0
            param(8)  = 1.039470d0
            param(9)  =-4.575040d0
            param(10) = 0.000000d0
            param(11) = 0.000000d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_HCTHP14)
            param(1)  = 4.0d0
            param(2)  = 2.82414d0     
            param(3)  = 8.21827d-002
            param(4)  = 3.18843d-002
            param(5)  = 4.56466d0     
            param(6)  =-1.78512d0     
            param(7)  =-13.5529d0     
            param(8)  = 2.39795d0     
            param(9)  = 13.3820d0     
            param(10) =-0.876909d0     
            param(11) =-3.174930d0
            call nwxc_c_b97_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_M05)
            param( 1) =   1.00000d0
            param( 2) =   3.78569d0
            param( 3) = -14.15261d0
            param( 4) =  -7.46589d0
            param( 5) =  17.94491d0
            param( 6) =   1.00000d0
            param( 7) =   3.77344d0
            param( 8) = -26.04463d0
            param( 9) =  30.69913d0
            param(10) =  -9.22695d0
            call nwxc_c_m05_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M05_2X)
            param( 1) =   1.00000d0
            param( 2) =   1.09297d0
            param( 3) =  -3.79171d0
            param( 4) =   2.82810d0
            param( 5) = -10.58909d0
            param( 6) =   1.00000d0
            param( 7) =  -3.05430d0
            param( 8) =   7.61854d0
            param( 9) =   1.47665d0
            param(10) = -11.92365d0
            call nwxc_c_m05_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M06)
c
            param( 1) =  -2.741539D+00
            param( 2) =  -6.720113D-01
            param( 3) =  -7.932688D-02
            param( 4) =   1.918681D-03
            param( 5) =  -2.032902D-03
            param( 6) =   0.000000D+00
            param( 7) =   4.905945D-01
            param( 8) =  -1.437348D-01
            param( 9) =   2.357824D-01
            param(10) =   1.871015D-03
            param(11) =  -3.788963D-03
            param(12) =   0.000000D+00
c
            param(13) =   3.741539D+00
            param(14) =   2.187098D+02
            param(15) =  -4.531252D+02
            param(16) =   2.936479D+02
            param(17) =  -6.287470D+01
            param(18) =   5.094055D-01
            param(19) =  -1.491085D+00
            param(20) =   1.723922D+01
            param(21) =  -3.859018D+01
            param(22) =   2.845044D+01
c
            call nwxc_c_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M06_HF)
c
            param( 1) =  -6.746338D-01
            param( 2) =  -1.534002D-01
            param( 3) =  -9.021521D-02
            param( 4) =  -1.292037D-03
            param( 5) =  -2.352983D-04
            param( 6) =   0.000000D+00
            param( 7) =   8.976746D-01
            param( 8) =  -2.345830D-01
            param( 9) =   2.368173D-01
            param(10) =  -9.913890D-04
            param(11) =  -1.146165D-02
            param(12) =   0.000000D+00
c
            param(13) =   1.674634D+00
            param(14) =   5.732017D+01
            param(15) =   5.955416D+01
            param(16) =  -2.311007D+02
            param(17) =   1.255199D+02
            param(18) =   1.023254D-01
            param(19) =  -2.453783D+00
            param(20) =   2.913180D+01
            param(21) =  -3.494358D+01
            param(22) =   2.315955D+01
c
            call nwxc_c_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M06_L)
c
            param( 1) =   3.957626D-01
            param( 2) =  -5.614546D-01
            param( 3) =   1.403963D-02
            param( 4) =   9.831442D-04
            param( 5) =  -3.577176D-03
            param( 6) =   0.000000D+00
            param( 7) =   4.650534D-01
            param( 8) =   1.617589D-01
            param( 9) =   1.833657D-01
            param(10) =   4.692100D-04
            param(11) =  -4.990573D-03
            param(12) =   0.000000D+00
c
            param(13) =   6.042374D-01
            param(14) =   1.776783D+02
            param(15) =  -2.513252D+02
            param(16) =   7.635173D+01
            param(17) =  -1.255699D+01
            param(18) =   5.349466D-01
            param(19) =   5.396620D-01
            param(20) =  -3.161217D+01
            param(21) =   5.149592D+01
            param(22) =  -2.919613D+01
c
            call nwxc_c_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M06_2X)
c
            param( 1) =   1.166404D-01
            param( 2) =  -9.120847D-02
            param( 3) =  -6.726189D-02
            param( 4) =   6.720580D-05
            param( 5) =   8.448011D-04
            param( 6) =   0.000000D+00
            param( 7) =   6.902145D-01
            param( 8) =   9.847204D-02
            param( 9) =   2.214797D-01
            param(10) =  -1.968264D-03
            param(11) =  -6.775479D-03
            param(12) =   0.000000D+00
c
            param(13) =   8.833596D-01
            param(14) =   3.357972D+01
            param(15) =  -7.043548D+01
            param(16) =   4.978271D+01
            param(17) =  -1.852891D+01
            param(18) =   3.097855D-01
            param(19) =  -5.528642D+00
            param(20) =   1.347420D+01
            param(21) =  -3.213623D+01
            param(22) =   2.846742D+01
c
            call nwxc_c_m06_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M08_HX)
c           parameters A
            param(1)  =    1.0000000D+00
            param(2)  =   -4.0661387D-01
            param(3)  =   -3.3232530D+00
            param(4)  =    1.5540980D+00
            param(5)  =    4.4248033D+01
            param(6)  =   -8.4351930D+01
            param(7)  =   -1.1955581D+02
            param(8)  =    3.9147081D+02
            param(9)  =    1.8363851D+02
            param(10) =   -6.3268223D+02
            param(11) =   -1.1297403D+02
            param(12) =    3.3629312D+02
c           parameters B
            param(13) =    1.3812334D+00
            param(14) =   -2.4683806D+00
            param(15) =   -1.1901501D+01
            param(16) =   -5.4112667D+01
            param(17) =    1.0055846D+01
            param(18) =    1.4800687D+02
            param(19) =    1.1561420D+02
            param(20) =    2.5591815D+02
            param(21) =    2.1320772D+02
            param(22) =   -4.8412067D+02
            param(23) =   -4.3430813D+02
            param(24) =    5.6627964D+01
            call nwxc_c_m11_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M08_SO)
c           parameters A
            param(1)  =   1.0000000D+00
            param(2)  =   0.0000000D+00
            param(3)  =  -3.9980886D+00
            param(4)  =   1.2982340D+01
            param(5)  =   1.0117507D+02
            param(6)  =  -8.9541984D+01
            param(7)  =  -3.5640242D+02
            param(8)  =   2.0698803D+02
            param(9)  =   4.6037780D+02
            param(10) =  -2.4510559D+02
            param(11) = -1.9638425D+02
            param(12) =  1.1881459D+02
c           parameters B
            param(13) =   1.0000000D+00
            param(14) =  -4.4117403D+00
            param(15) =  -6.4128622D+00
            param(16) =   4.7583635D+01
            param(17) =   1.8630053D+02
            param(18) =  -1.2800784D+02
            param(19) =  -5.5385258D+02
            param(20) =   1.3873727D+02
            param(21) =   4.1646537D+02
            param(22) =  -2.6626577D+02
            param(23) =   5.6676300D+01
            param(24) =   3.1673746D+02
            call nwxc_c_m11_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M11)
c            parameters A
            param(1)  =  1.0000000D+00
            param(2)  =  0.0000000D+00
            param(3)  = -3.8933250D+00
            param(4)  = -2.1688455D+00
            param(5)  =  9.3497200D+00
            param(6)  = -1.9845140D+01
            param(7)  =  2.3455253D+00
            param(8)  =  7.9246513D+01
            param(9)  =  9.6042757D+00
            param(10) = -6.7856719D+01
            param(11) = -9.1841067D+00
            param(12) =  0.0000000D+00
c           parameters B
            param(13) =  7.2239798D-01
            param(14) =  4.3730564D-01
            param(15) = -1.6088809D+01
            param(16) = -6.5542437D+01
            param(17) =  3.2057230D+01
            param(18) =  1.8617888D+02
            param(19) =  2.0483468D+01
            param(20) = -7.0853739D+01
            param(21) =  4.4483915D+01
            param(22) = -9.4484747D+01
            param(23) = -1.1459868D+02
            param(24) =  0.0000000D+00
            call nwxc_c_m11_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_M11_L)
c           parameters A
            param(1)  =  1.000000D+00
            param(2)  =  0.000000D+00
            param(3)  =  2.750880D+00
            param(4)  = -1.562287D+01
            param(5)  =  9.363381D+00
            param(6)  =  2.141024D+01
            param(7)  = -1.424975D+01
            param(8)  = -1.134712D+01
            param(9)  =  1.022365D+01
            param(10) =  0.000000D+00
            param(11) =  0.000000D+00
            param(12) =  0.000000D+00
c           parameters B
            param(13) =  1.000000D+00
            param(14) = -9.082060D+00
            param(15) =  6.134682D+00
            param(16) = -1.333216D+01
            param(17) = -1.464115D+01
            param(18) =  1.713143D+01
            param(19) =  2.480738D+00
            param(20) = -1.007036D+01
            param(21) = -1.117521D-01
            param(22) =  0.000000D+00
            param(23) =  0.000000D+00
            param(24) =  0.000000D+00
            call nwxc_c_m11_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_OP)
            param(1) = 2.3670d0
            call nwxc_c_op_d2(nwxc_k_becke88_d2,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_OPT)
            call nwxc_c_opt_d2(nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PW6B95)
            param(1) = 0.00262d0
            param(2) = 0.03668d0
            call nwxc_c_b95_d2(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_PWB6K)
            param(1) = 0.00353d0
            param(2) = 0.04120d0
            call nwxc_c_b95_d2(param,nwxc_rho_tol,ipol,nq,
     +                      nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_PW91LDA)
            call nwxc_c_pw91lda_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                             rho,f)
          case (NWXCP_C_LYP)
            call nwxc_c_lyp_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,f)
          case (NWXCP_C_PZ81)
            call nwxc_c_perdew81_d2(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),rho,f)
          case (NWXCP_C_P86)
            call nwxc_c_perdew86_d2(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_P91)
            call nwxc_c_p91_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,f)
          case (NWXCP_C_P91_VWN5)
            call nwxc_c_p91_vwn5_d2(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),rho,gamma,f,dfdr)
          case (NWXCP_C_PBE)
            param(1) = 0.066724550603149d0
            call nwxc_c_pbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PBESOL)
            param(1) = 0.046d0
            call nwxc_c_pbe_d2(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PKZB)
            call nwxc_c_pkzb99_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,tau,f)
          case (NWXCP_C_SPBE)
            call nwxc_c_spbe_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_C_SOGGA11)
            param(1)  =  5.00000d-01
            param(2)  = -4.62334D+00
            param(3)  =  8.00410D+00
            param(4)  = -130.226D+00
            param(5)  =  38.2685D+00
            param(6)  =  69.5599D+00
            param(7)  =  5.00000d-01
            param(8)  =  3.62334D+00
            param(9)  =  9.36393D+00
            param(10) =  34.5114D+00
            param(11) = -18.5684D+00
            param(12) = -0.16519D+00
            call nwxc_c_sogga_d2(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_SOGGA11_X)
            param(1)  =  5.00000d-01
            param(2)  =  7.82439d+01
            param(3)  =  2.57211d+01
            param(4)  = -1.38830d+01
            param(5)  = -9.87375d+00
            param(6)  = -1.41357d+01
            param(7)  =  5.00000d-01
            param(8)  = -7.92439d+01
            param(9)  =  1.63725d+01
            param(10) =  2.08129d+00
            param(11) =  7.50769d+00
            param(12) = -1.01861d+01
            call nwxc_c_sogga_d2(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_TPSS)
            call nwxc_c_tpss03_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                            rho,gamma,tau,f)
          case (NWXCP_C_VS98)
            param(1)  =  7.035010d-01
            param(2)  =  7.694574d-03
            param(3)  =  5.152765d-02
            param(4)  =  3.394308d-05
            param(5)  = -1.269420d-03
            param(6)  =  1.296118d-03
            param(7)  =  3.270912d-01
            param(8)  = -3.228915d-02
            param(9)  = -2.942406d-02
            param(10) =  2.134222d-03
            param(11) = -5.451559d-03
            param(12) =  1.577575d-02
            call nwxc_c_vs98_d2(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_VWN1)
            call nwxc_c_vwn1_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN1_RPA)
            call nwxc_c_vwn1_rpa_d2(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),
     +                              rho,f)
          case (NWXCP_C_VWN2)
            call nwxc_c_vwn2_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN3)
            call nwxc_c_vwn3_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN4)
            call nwxc_c_vwn4_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN5)
            call nwxc_c_vwn5_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_BOP)
            param(1) = 2.3670d0
            call nwxc_c_op_d2(nwxc_k_becke88_d2,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_KT1)
            call nwxc_xc_kt1_d2(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_PBEOP)
            param(1) = 2.3789d0
            call nwxc_c_op_d2(nwxc_k_pbe96_d2,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_SOP)
            param(1) = 2.5654d0
            call nwxc_c_op_d2(nwxc_k_dirac_d2,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case default
            if (nwxc_oroot) then
              write(*,*)"nwxc_eval_df2: invalid functional",
     +                  nwxc_ids(it)
            endif
            call nwxc_printP()
            call errquit("nwxc_eval_df2: invalid functional",
     +                   nwxc_ids(it),0)
        end select
      enddo
C
      end
C>
C> @}
C>
C> \ingroup nwad_api
C> @{
C>
C> \brief The functional, 1st, 2nd, and 3rd order partial derivative evaluation
C>
C> In actual fact this routine only sets up the memory and then calls the
C> driver routine to drive the actual functional evaluation. The lack of 
C> a Fortran90 typecast drives this code structure.
C>
      subroutine nwxc_eval_df3(ipol,nq,rho,gamma,tau,f,dfdr,dfdr2,dfdr3,
     +                         dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
      use nwad3
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
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
      double precision dfdr2(nq,NCOL_AMAT2) !< [Output] 2nd derivative wrt rho
      double precision dfdg2(nq,NCOL_CMAT2) !< [Output] 2nd derivative wrt gamma
      double precision dfdt2(nq,NCOL_MMAT2) !< [Output] 2nd derivative wrt tau
      double precision dfdr3(nq,NCOL_AMAT3) !< [Output] 3rd derivative wrt rho
      double precision dfdg3(nq,NCOL_CMAT3) !< [Output] 3rd derivative wrt gamma
      double precision dfdt3(nq,NCOL_MMAT3) !< [Output] 3rd derivative wrt tau
c
      type(nwad_dble) :: memory_test_nwad(2)
      double precision   memory_test_dble(2)
      integer length
      integer l_rho, k_rho !< Index for density
      integer l_gam, k_gam !< Index for gradient
      integer l_tau, k_tau !< Index for kinetic energy density
      integer l_fun, k_fun !< Index for functional
      integer l_scr, k_scr !< Index for functional
      integer nvar         !< The number of variables in the functional
                           !< LDA: nvar=ipol, GGA: nvar=2*ipol+1, 
                           !< MGGA: nvar=3*ipol+1
      integer ndrv         !< The number of partial derivatives
      integer ipl          !< Hack version of ipol
C
      logical nwxc_is_gga  !< Is the functional a GGA
      logical nwxc_is_mgga !< Is the functional a meta-GGA
c
c     Work out nvar for this functional
c
      ipl = 2
      nvar = ipol
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        nvar = nvar + ipol+1
        if (nwxc_is_mgga()) then
          nvar = nvar + ipol
        endif
      endif
      ndrv = nvar*(nvar+1)*(nvar+2)/6
c
c     Work out how many double precision floating point words are needed
c     to represent a single nwad_dble
c
      length = (loc(memory_test_dble(2))-loc(memory_test_dble(1)))
      length = (loc(memory_test_nwad(2))-loc(memory_test_nwad(1))
     +          +length-1)/length
      if(.not.ma_push_get(mt_dbl,length*nq*ipl,"nwad rho",l_rho,k_rho))
     +  call errquit("nwxc_eval_df3: cannot allocate nwad rho",
     +    ma_sizeof(mt_dbl,length*nq*ipl,mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*(ipl+1),"nwad gamma",
     +                    l_gam,k_gam))
     +  call errquit("nwxc_eval_df3: cannot allocate nwad gamma",
     +    ma_sizeof(mt_dbl,length*nq*(ipl+1),mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*ipl,"nwad tau",l_tau,k_tau))
     +  call errquit("nwxc_eval_df3: cannot allocate nwad tau",
     +    ma_sizeof(mt_dbl,length*nq*ipl,mt_byte),MA_ERR)
      if(.not.ma_push_get(mt_dbl,length*nq*ndrv,"nwad fun",l_fun,k_fun))
     +  call errquit("nwxc_eval_df3: cannot allocate nwad fun",
     +    ma_sizeof(mt_dbl,length*nq*ndrv,mt_byte),MA_ERR)
c
      call nwxc_eval_df3_driver(ipol,nq,rho,gamma,tau,dbl_mb(k_rho),
     +     dbl_mb(k_gam),dbl_mb(k_tau),dbl_mb(k_fun),
     +     nvar,ndrv,
     +     f,dfdr,dfdr2,dfdr3,dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
c
      if(.not.ma_pop_stack(l_fun)) call errquit(
     +  "nwxc_eval_df3: cannot deallocate nwad fun",0,MA_ERR)
      if(.not.ma_pop_stack(l_tau)) call errquit(
     +  "nwxc_eval_df3: cannot deallocate nwad tau",0,MA_ERR)
      if(.not.ma_pop_stack(l_gam)) call errquit(
     +  "nwxc_eval_df3: cannot deallocate nwad gam",0,MA_ERR)
      if(.not.ma_pop_stack(l_rho)) call errquit(
     +  "nwxc_eval_df3: cannot deallocate nwad rho",0,MA_ERR)
      end
C>
C> @}
C>
C> \ingroup nwad_priv
C> @{
C>
C> \brief Driver routine for the functional, 1st, 2nd and 3rd derivative
C> evaluation
C>
C> This driver routine initializes the active variables and loops over the
C> functional evaluations. Afterwards it interpolates the partial derivatives
C> and returns the values of the functional and its 1st order partial
C> derivatives. Obviously, interpolating the 1st order partial derivatives is
C> particularly simply, but this scheme applies to all orders of derivatives.
C> Due to different data types being required the memory for this routine
C> needs to be allocated in the calling routine. Calling this subroutine
C> effectively implements a typecast.
C>
      subroutine nwxc_eval_df3_driver(ipol,nq,rho,gamma,tau,
     +           nwad_rho,nwad_gam,nwad_tau,nwad_f,nvar,ndrv,
     +           f,dfdr,dfdr2,dfdr3,dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
      use nwad3
      implicit none
#include "errquit.fh"
#include "nwxc_param.fh"
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of points
      integer nvar !< [Input] The number of independent variables
      integer ndrv !< [Input] The number of 2nd order partial derivatives
C
      double precision rho(nq,ipol) !< [Input] Density
      double precision gamma(nq,3)  !< [Input] |Density gradient|^2
      double precision tau(nq,ipol) !< [Input] Kinetic energy 
                                    !< density
C
      double precision f(nq)           !< [Output] Energy
      double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
      double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
      double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
      double precision dfdr2(nq,NCOL_AMAT2) !< [Output] 2nd derivative wrt rho
      double precision dfdg2(nq,NCOL_CMAT2) !< [Output] 2nd derivative wrt gamma
      double precision dfdt2(nq,NCOL_MMAT2) !< [Output] 2nd derivative wrt tau
      double precision dfdr3(nq,NCOL_AMAT3) !< [Output] 3rd derivative wrt rho
      double precision dfdg3(nq,NCOL_CMAT3) !< [Output] 3rd derivative wrt gamma
      double precision dfdt3(nq,NCOL_MMAT3) !< [Output] 3rd derivative wrt tau
c
      type(nwad_dble) :: nwad_rho(nq,2)
      type(nwad_dble) :: nwad_gam(nq,3)
      type(nwad_dble) :: nwad_tau(nq,2)
      type(nwad_dble) :: nwad_f(nq,ndrv)
C
      logical nwxc_is_gga  !< Is the functional a GGA
      logical nwxc_is_mgga !< Is the functional a meta-GGA
c
      integer iq, ip, id, ir, it
      integer ip2, ir2, ipr, id2
      integer ip3, ir3, it3, ip2r, ip2t, ipr2, ir2t, ipt2, irt2, iprt
c
      if (ipol.eq.1) then
        do ip = 1, 2
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
          enddo
        enddo
      else  ! ipol.eq.1
        do ip = 1, ipol
          do iq = 1, nq
            nwad_rho(iq,ip) = inactive(rho(iq,ip))
          enddo
        enddo
      endif ! ipol.eq.1
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        if (ipol.eq.1) then
          do ip = 1, 3
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(0.25d0*gamma(iq,G_TT))
            enddo
          enddo
        else  ! ipol.eq.1
          do ip = 1, ipol+1
            do iq = 1, nq
              nwad_gam(iq,ip) = inactive(gamma(iq,ip))
            enddo
          enddo
        endif ! ipol.eq.1
        if (nwxc_is_mgga()) then
          if (ipol.eq.1) then
            do ip = 1, 2
              do iq = 1, nq
                nwad_tau(iq,ip) = inactive(1.0d0*tau(iq,T_T))
              enddo
            enddo
          else
            do ip = 1, ipol
              do iq = 1, nq
                nwad_tau(iq,ip) = inactive(2.0d0*tau(iq,ip))
              enddo
            enddo
          endif
        endif
      endif
c
c     Do closed shell and open shell separately
c
      if (ipol.eq.1) then
c
c       Do closed shell
c
c       First do d3f/dx3
c
        do ip = 1, nvar
          id = ip*(ip+1)*(ip+2)/6
          if (ip.eq.1) then
            do iq = 1, nq
              nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T),3)
            enddo
          elseif (ip.le.3) then
            do iq = 1, nq
              nwad_gam(iq,ip-1) = active(0.25d0*gamma(iq,G_TT),3)
            enddo
          else
            do iq = 1, nq
              nwad_tau(iq,ip-3) = active(1.0d0*tau(iq,T_T),3)
            enddo
          endif
          call nwxc_eval_df3_doit(2,nq,nwad_rho,nwad_gam,
     +         nwad_tau,nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +         dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
          if (ip.eq.1) then
            do iq = 1, nq
               nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
            enddo
          elseif (ip.le.3) then
            do iq = 1, nq
              nwad_gam(iq,ip-1) = inactive(0.25d0*gamma(iq,G_TT))
            enddo
          else
            do iq = 1, nq
              nwad_tau(iq,ip-3) = inactive(1.0d0*tau(iq,T_T))
            enddo
          endif
        enddo
c
c       Second do d3f/dx2dy
c
        do ip = 1, nvar
          do ir = 1, ip-1
            id = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
            if (ip.eq.1) then
              do iq = 1, nq
                nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T),2)
              enddo
            elseif (ip.le.3) then
              do iq = 1, nq
                nwad_gam(iq,ip-1) = active(0.25d0*gamma(iq,G_TT),2)
              enddo
            elseif (ip.eq.4) then
              do iq = 1, nq
                nwad_tau(iq,ip-3) = active(1.0d0*tau(iq,T_T),2)
              enddo
            endif
            if (ir.eq.1) then
              do iq = 1, nq
                nwad_rho(iq,ir) = active(0.5d0*rho(iq,R_T))
              enddo
            elseif (ir.le.3) then
              do iq = 1, nq
                nwad_gam(iq,ir-1) = active(0.25d0*gamma(iq,G_TT))
              enddo
            elseif (ir.eq.4) then
              do iq = 1, nq
                nwad_tau(iq,ir-3) = active(1.0d0*tau(iq,T_T))
              enddo
            endif
            call nwxc_eval_df3_doit(2,nq,nwad_rho,nwad_gam,
     +           nwad_tau,nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +           dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
             if (ip.eq.1) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
               enddo
             elseif (ip.le.3) then
               do iq = 1, nq
                 nwad_gam(iq,ip-1) = inactive(0.25d0*gamma(iq,G_TT))
               enddo
             elseif (ip.eq.4) then
               do iq = 1, nq
                 nwad_tau(iq,ip-3) = inactive(1.0d0*tau(iq,T_T))
               enddo
             endif
             if (ir.eq.1) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = inactive(0.5d0*rho(iq,R_T))
               enddo
             elseif (ir.le.3) then
               do iq = 1, nq
                 nwad_gam(iq,ir-1) = inactive(0.25d0*gamma(iq,G_TT))
               enddo
             elseif (ir.eq.4) then
               do iq = 1, nq
                 nwad_tau(iq,ir-3) = inactive(1.0d0*tau(iq,T_T))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Third do d3f/dxdy2
c
         do ip = 1, nvar
           do ir = 1, ip-1
             id = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             if (ip.eq.1) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T))
               enddo
             elseif (ip.le.3) then
               do iq = 1, nq
                 nwad_gam(iq,ip-1) = active(0.25d0*gamma(iq,G_TT))
               enddo
             elseif (ip.eq.4) then
               do iq = 1, nq
                 nwad_tau(iq,ip-3) = active(1.0d0*tau(iq,T_T))
               enddo
             endif
             if (ir.eq.1) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = active(0.5d0*rho(iq,R_T),2)
               enddo
             elseif (ir.le.3) then
               do iq = 1, nq
                 nwad_gam(iq,ir-1) = active(0.25d0*gamma(iq,G_TT),2)
               enddo
             elseif (ir.eq.4) then
               do iq = 1, nq
                 nwad_tau(iq,ir-3) = active(1.0d0*tau(iq,T_T),2)
               enddo
             endif
             call nwxc_eval_df3_doit(2,nq,nwad_rho,nwad_gam,
     +            nwad_tau,nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +            dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
             if (ip.eq.1) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
               enddo
             elseif (ip.le.3) then
               do iq = 1, nq
                 nwad_gam(iq,ip-1) = inactive(0.25d0*gamma(iq,G_TT))
               enddo
             elseif (ip.eq.4) then
               do iq = 1, nq
                 nwad_tau(iq,ip-3) = inactive(1.0d0*tau(iq,T_T))
               enddo
             endif
             if (ir.eq.1) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = inactive(0.5d0*rho(iq,R_T))
               enddo
             elseif (ir.le.3) then
               do iq = 1, nq
                 nwad_gam(iq,ir-1) = inactive(0.25d0*gamma(iq,G_TT))
               enddo
             elseif (ir.eq.4) then
               do iq = 1, nq
                 nwad_tau(iq,ir-3) = inactive(1.0d0*tau(iq,T_T))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Fourth do d3f/dxdydz
c
         do ip = 1, nvar
           do ir = 1, ip-1
             do it = 1, ir-1
               id = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+it
               if (ip.eq.1) then
                 do iq = 1, nq
                   nwad_rho(iq,ip) = active(0.5d0*rho(iq,R_T))
                 enddo
               elseif (ip.le.3) then
                 do iq = 1, nq
                   nwad_gam(iq,ip-1) = active(0.25d0*gamma(iq,G_TT))
                 enddo
               elseif (ip.eq.4) then
                 do iq = 1, nq
                   nwad_tau(iq,ip-3) = active(1.0d0*tau(iq,T_T))
                 enddo
               endif
               if (ir.eq.1) then
                 do iq = 1, nq
                   nwad_rho(iq,ir) = active(0.5d0*rho(iq,R_T))
                 enddo
               elseif (ir.le.3) then
                 do iq = 1, nq
                   nwad_gam(iq,ir-1) = active(0.25d0*gamma(iq,G_TT))
                 enddo
               elseif (ir.eq.4) then
                 do iq = 1, nq
                   nwad_tau(iq,ir-3) = active(1.0d0*tau(iq,T_T))
                 enddo
               endif
               if (it.eq.1) then
                 do iq = 1, nq
                   nwad_rho(iq,it) = active(0.5d0*rho(iq,R_T))
                 enddo
               elseif (it.le.3) then
                 do iq = 1, nq
                   nwad_gam(iq,it-1) = active(0.25d0*gamma(iq,G_TT))
                 enddo
               elseif (it.eq.4) then
                 do iq = 1, nq
                   nwad_tau(iq,it-3) = active(1.0d0*tau(iq,T_T))
                 enddo
               endif
               call nwxc_eval_df3_doit(2,nq,nwad_rho,
     +              nwad_gam,nwad_tau,nwad_f(1,id),
     +              dfdr,dfdr2,dfdr3,dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
               if (ip.eq.1) then
                 do iq = 1, nq
                   nwad_rho(iq,ip) = inactive(0.5d0*rho(iq,R_T))
                 enddo
               elseif (ip.le.3) then
                 do iq = 1, nq
                   nwad_gam(iq,ip-1) = inactive(0.25d0*gamma(iq,G_TT))
                 enddo
               elseif (ip.eq.4) then
                 do iq = 1, nq
                   nwad_tau(iq,ip-3) = inactive(1.0d0*tau(iq,T_T))
                 enddo
               endif
               if (ir.eq.1) then
                 do iq = 1, nq
                   nwad_rho(iq,ir) = inactive(0.5d0*rho(iq,R_T))
                 enddo
               elseif (ir.le.3) then
                 do iq = 1, nq
                   nwad_gam(iq,ir-1) = inactive(0.25d0*gamma(iq,G_TT))
                 enddo
               elseif (ir.eq.4) then
                 do iq = 1, nq
                   nwad_tau(iq,ir-3) = inactive(1.0d0*tau(iq,T_T))
                 enddo
               endif
               if (it.eq.1) then
                 do iq = 1, nq
                   nwad_rho(iq,it) = inactive(0.5d0*rho(iq,R_T))
                 enddo
               elseif (it.le.3) then
                 do iq = 1, nq
                   nwad_gam(iq,it-1) = inactive(0.25d0*gamma(iq,G_TT))
                 enddo
               elseif (it.eq.4) then
                 do iq = 1, nq
                   nwad_tau(iq,it-3) = inactive(1.0d0*tau(iq,T_T))
                 enddo
               endif
             enddo ! it
           enddo ! ir
         enddo ! ip
c
c        Done with all functional evaluations, now interpolate the results
c
c        First get the functional itself
c
         do iq = 1, nq
           f(iq) = value(nwad_f(iq,1))
         enddo
c
c        Interpolate the first order partial derivatives
c
         do ip = 1, nvar
           id = ip*(ip+1)*(ip+2)/6
           if (ip.eq.1) then
             do iq = 1, nq
               dfdr(iq,ip) = inter_d3_dx(nwad_f(iq,id))
             enddo
           elseif (ip.le.3) then
             do iq = 1, nq
               dfdg(iq,ip-1) = inter_d3_dx(nwad_f(iq,id))
             enddo
           elseif (ip.eq.4) then
             do iq = 1, nq
               dfdt(iq,ip-3) = inter_d3_dx(nwad_f(iq,id))
             enddo
           endif
         enddo
c
c        Interpolate the second order partial derivatives
c
c        First do d2fdx2
c
         do ip = 1, nvar
           id = ip*(ip+1)*(ip+2)/6
           if (ip.eq.1) then
             do iq = 1, nq
               dfdr2(iq,D2_RA_RA) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           elseif (ip.eq.2) then
             do iq = 1, nq
               dfdg2(iq,D2_GAA_GAA) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           elseif (ip.eq.3) then
             do iq = 1, nq
               dfdg2(iq,D2_GAB_GAB) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           elseif (ip.eq.4) then
             do iq = 1, nq
               dfdt2(iq,D2_TA_TA) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           endif
         enddo
c
c        Second do d2fdxy
c
         do ip = 1, nvar
           do ir = 1, ip-1
             ip3  = ip*(ip+1)*(ip+2)/6
             ir3  = ir*(ir+1)*(ir+2)/6
             ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             id2  = D2_RA_RA
             if (ip.eq.2.and.ir.eq.1) then
               id2 = D2_RA_GAA
             elseif (ip.eq.3.and.ir.eq.1) then
               id2 = D2_RA_GAB
             elseif (ip.eq.3.and.ir.eq.2) then
               id2 = D2_GAA_GAB
             elseif (ip.eq.4.and.ir.eq.1) then
               id2 = D2_RA_TA
             elseif (ip.eq.4.and.ir.eq.2) then
               id2 = D2_GAA_TA
             elseif (ip.eq.4.and.ir.eq.3) then
               id2 = D2_GAB_TA
             endif
             if (ip.le.3) then
               do iq = 1, nq
                 dfdg2(iq,id2) = inter_d3_dxy(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             else
               do iq = 1, nq
                 dfdt2(iq,id2) = inter_d3_dxy(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Interpolate the third order partial derivatives
c
c        First do d3fdx3
c
         do ip = 1, nvar
           id = ip*(ip+1)*(ip+2)/6
           if (ip.eq.1) then
             do iq = 1, nq
               dfdr3(iq,D3_RA_RA_RA) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           elseif (ip.eq.2) then
             do iq = 1, nq
               dfdg3(iq,D3_GAA_GAA_GAA) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           elseif (ip.eq.3) then
             do iq = 1, nq
               dfdg3(iq,D3_GAB_GAB_GAB) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           elseif (ip.eq.4) then
             do iq = 1, nq
               dfdt3(iq,D3_TA_TA_TA) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           endif
         enddo
c
c        Second do d3fdx2y
c
         do ip = 1, nvar
           do ir = 1, ip-1
             ip3  = ip*(ip+1)*(ip+2)/6
             ir3  = ir*(ir+1)*(ir+2)/6
             ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             id2  = D3_RA_RA_RA
             if (ip.eq.2.and.ir.eq.1) then
               id2 = D3_RA_GAA_GAA
             elseif (ip.eq.3.and.ir.eq.1) then
               id2 = D3_RA_GAB_GAB
             elseif (ip.eq.3.and.ir.eq.2) then
               id2 = D3_GAA_GAB_GAB
             elseif (ip.eq.4.and.ir.eq.1) then
               id2 = D3_RA_TA_TA
             elseif (ip.eq.4.and.ir.eq.2) then
               id2 = D3_GAA_TA_TA
             elseif (ip.eq.4.and.ir.eq.3) then
               id2 = D3_GAB_TA_TA
             endif
             if (ip.le.3) then
               do iq = 1, nq
                 dfdg3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             else
               do iq = 1, nq
                 dfdt3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Third do d3fdxy2
c
         do ip = 1, nvar
           do ir = 1, ip-1
             ip3  = ip*(ip+1)*(ip+2)/6
             ir3  = ir*(ir+1)*(ir+2)/6
             ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             id2  = D3_RA_RA_RA
             if (ip.eq.2.and.ir.eq.1) then
               id2 = D3_RA_RA_GAA
             elseif (ip.eq.3.and.ir.eq.1) then
               id2 = D3_RA_RA_GAB
             elseif (ip.eq.3.and.ir.eq.2) then
               id2 = D3_GAA_GAA_GAB
             elseif (ip.eq.4.and.ir.eq.1) then
               id2 = D3_RA_RA_TA
             elseif (ip.eq.4.and.ir.eq.2) then
               id2 = D3_GAA_GAA_TA
             elseif (ip.eq.4.and.ir.eq.3) then
               id2 = D3_GAB_GAB_TA
             endif
             if (ip.le.3) then
               do iq = 1, nq
                 dfdg3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ir3),
     +                           nwad_f(iq,ipr2),nwad_f(iq,ip2r),
     +                           nwad_f(iq,ip3))
               enddo
             else
               do iq = 1, nq
                 dfdt3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ir3),
     +                           nwad_f(iq,ipr2),nwad_f(iq,ip2r),
     +                           nwad_f(iq,ip3))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Fourth do d3fdxyz
c
         do ip = 1, nvar
           do ir = 1, ip-1
             do it = 1, ir-1
               ip3  = ip*(ip+1)*(ip+2)/6
               ir3  = ir*(ir+1)*(ir+2)/6
               it3  = it*(it+1)*(it+2)/6
               ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
               ip2t = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+it
               ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
               ir2t = (ir-1)*ir*(ir+1)/6+ir*(ir-1)/2+it
               ipt2 = (ip-1)*ip*(ip+1)/6+it*(it-1)/2+it
               irt2 = (ir-1)*ir*(ir+1)/6+it*(it-1)/2+it
               iprt = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+it
               id2  = D3_RA_RA_RA
               if (ip.eq.3.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_GAA_GAB
               elseif (ip.eq.4.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_GAA_TA
               elseif (ip.eq.4.and.ir.eq.3.and.it.eq.1) then
                 id2 = D3_RA_GAB_TA
               elseif (ip.eq.4.and.ir.eq.3.and.it.eq.2) then
                 id2 = D3_GAA_GAB_TA
               endif
               if (ip.le.3) then
                 do iq = 1, nq
                   dfdg3(iq,id2)
     +             = inter_d3_dxyz(nwad_f(iq,ip3),nwad_f(iq,ir3),
     +                             nwad_f(iq,it3),nwad_f(iq,ip2r),
     +                             nwad_f(iq,ip2t),nwad_f(iq,ipr2),
     +                             nwad_f(iq,ir2t),nwad_f(iq,ipt2),
     +                             nwad_f(iq,irt2),nwad_f(iq,iprt))
                 enddo
               else
                 do iq = 1, nq
                   dfdt3(iq,id2)
     +             = inter_d3_dxyz(nwad_f(iq,ip3),nwad_f(iq,ir3),
     +                             nwad_f(iq,it3),nwad_f(iq,ip2r),
     +                             nwad_f(iq,ip2t),nwad_f(iq,ipr2),
     +                             nwad_f(iq,ir2t),nwad_f(iq,ipt2),
     +                             nwad_f(iq,irt2),nwad_f(iq,iprt))
                 enddo
               endif
             enddo ! it
           enddo ! ir
         enddo ! ip
c
       elseif (ipol.eq.2) then
c
c        Do open shell
c
c        First do d3f/dx3
c
         do ip = 1, nvar
           id = ip*(ip+1)/2
           if (ip.le.2) then
             do iq = 1, nq
               nwad_rho(iq,ip) = active(rho(iq,ip),3)
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               nwad_gam(iq,ip-2) = active(gamma(iq,ip-2),3)
             enddo
           else
             do iq = 1, nq
               nwad_tau(iq,ip-5) = active(2.0d0*tau(iq,ip-5),3)
             enddo
           endif
           call nwxc_eval_df3_doit(ipol,nq,nwad_rho,nwad_gam,nwad_tau,
     +                             nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +                             dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
           if (ip.le.2) then
             do iq = 1, nq
               nwad_rho(iq,ip) = inactive(rho(iq,ip))
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               nwad_gam(iq,ip-2) = inactive(gamma(iq,ip-2))
             enddo
           else
             do iq = 1, nq
               nwad_tau(iq,ip-5) = inactive(2.0d0*tau(iq,ip-5))
             enddo
           endif
         enddo
c
c        Second do d3f/dx2dy
c
         do ip = 1, nvar
           do ir = 1, ip-1
             id = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             if (ip.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = active(rho(iq,ip),2)
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ip-2) = active(gamma(iq,ip-2),2)
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ip-5) = active(2.0d0*tau(iq,ip-5),2)
               enddo
             endif
             if (ir.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = active(rho(iq,ir))
               enddo
             elseif (ir.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ir-2) = active(gamma(iq,ir-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ir-5) = active(2.0d0*tau(iq,ir-5))
               enddo
             endif
             call nwxc_eval_df3_doit(ipol,nq,nwad_rho,nwad_gam,nwad_tau,
     +                               nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +                               dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
             if (ip.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = inactive(rho(iq,ip))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ip-2) = inactive(gamma(iq,ip-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ip-5) = inactive(2.0d0*tau(iq,ip-5))
               enddo
             endif
             if (ir.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = inactive(rho(iq,ir))
               enddo
             elseif (ir.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ir-2) = inactive(gamma(iq,ir-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ir-5) = inactive(2.0d0*tau(iq,ir-5))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Third do d3f/dxdy2
c
         do ip = 1, nvar
           do ir = 1, ip-1
             id = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             if (ip.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = active(rho(iq,ip))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ip-2) = active(gamma(iq,ip-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ip-5) = active(2.0d0*tau(iq,ip-5))
               enddo
             endif
             if (ir.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = active(rho(iq,ir),2)
               enddo
             elseif (ir.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ir-2) = active(gamma(iq,ir-2),2)
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ir-5) = active(2.0d0*tau(iq,ir-5),2)
               enddo
             endif
             call nwxc_eval_df3_doit(ipol,nq,nwad_rho,nwad_gam,nwad_tau,
     +                               nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +                               dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
             if (ip.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ip) = inactive(rho(iq,ip))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ip-2) = inactive(gamma(iq,ip-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ip-5) = inactive(2.0d0*tau(iq,ip-5))
               enddo
             endif
             if (ir.le.2) then
               do iq = 1, nq
                 nwad_rho(iq,ir) = inactive(rho(iq,ir))
               enddo
             elseif (ir.le.5) then
               do iq = 1, nq
                 nwad_gam(iq,ir-2) = inactive(gamma(iq,ir-2))
               enddo
             else
               do iq = 1, nq
                 nwad_tau(iq,ir-5) = inactive(2.0d0*tau(iq,ir-5))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Fourth do d3f/dxdydz
c
         do ip = 1, nvar
           do ir = 1, ip-1
             do it = 1, ir-1
               id = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+it
               if (ip.le.2) then
                 do iq = 1, nq
                   nwad_rho(iq,ip) = active(rho(iq,ip))
                 enddo
               elseif (ip.le.5) then
                 do iq = 1, nq
                   nwad_gam(iq,ip-2) = active(gamma(iq,ip-2))
                 enddo
               else
                 do iq = 1, nq
                   nwad_tau(iq,ip-5) = active(2.0d0*tau(iq,ip-5))
                 enddo
               endif
               if (ir.le.2) then
                 do iq = 1, nq
                   nwad_rho(iq,ir) = active(rho(iq,ir))
                 enddo
               elseif (ir.le.5) then
                 do iq = 1, nq
                   nwad_gam(iq,ir-2) = active(gamma(iq,ir-2))
                 enddo
               else
                 do iq = 1, nq
                   nwad_tau(iq,ir-5) = active(2.0d0*tau(iq,ir-5))
                 enddo
               endif
               if (it.le.2) then
                 do iq = 1, nq
                   nwad_rho(iq,it) = active(rho(iq,it))
                 enddo
               elseif (it.le.5) then
                 do iq = 1, nq
                   nwad_gam(iq,it-2) = active(gamma(iq,it-2))
                 enddo
               else
                 do iq = 1, nq
                   nwad_tau(iq,it-5) = active(2.0d0*tau(iq,it-5))
                 enddo
               endif
               call nwxc_eval_df3_doit(ipol,nq,nwad_rho,nwad_gam,
     +              nwad_tau,nwad_f(1,id),dfdr,dfdr2,dfdr3,
     +              dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
               if (ip.le.2) then
                 do iq = 1, nq
                   nwad_rho(iq,ip) = inactive(rho(iq,ip))
                 enddo
               elseif (ip.le.5) then
                 do iq = 1, nq
                   nwad_gam(iq,ip-2) = inactive(gamma(iq,ip-2))
                 enddo
               else
                 do iq = 1, nq
                   nwad_tau(iq,ip-5) = inactive(2.0d0*tau(iq,ip-5))
                 enddo
               endif
               if (ir.le.2) then
                 do iq = 1, nq
                   nwad_rho(iq,ir) = inactive(rho(iq,ir))
                 enddo
               elseif (ir.le.5) then
                 do iq = 1, nq
                   nwad_gam(iq,ir-2) = inactive(gamma(iq,ir-2))
                 enddo
               else
                 do iq = 1, nq
                   nwad_tau(iq,ir-5) = inactive(2.0d0*tau(iq,ir-5))
                 enddo
               endif
               if (it.le.2) then
                 do iq = 1, nq
                   nwad_rho(iq,it) = inactive(rho(iq,it))
                 enddo
               elseif (it.le.5) then
                 do iq = 1, nq
                   nwad_gam(iq,it-2) = inactive(gamma(iq,it-2))
                 enddo
               else
                 do iq = 1, nq
                   nwad_tau(iq,it-5) = inactive(2.0d0*tau(iq,it-5))
                 enddo
               endif
             enddo ! it
           enddo ! ir
         enddo ! ip
c
c        Done with all functional evaluations, now interpolate the results
c
c        First get the functional itself
c
         do iq = 1, nq
           f(iq) = value(nwad_f(iq,1))
         enddo
c
c        Interpolate the first order partial derivatives
c
         do ip = 1, nvar
           id = ip*(ip+1)*(ip+2)/6
           if (ip.le.2) then
             do iq = 1, nq
               dfdr(iq,ip) = inter_d3_dx(nwad_f(iq,id))
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               dfdg(iq,ip-2) = inter_d3_dx(nwad_f(iq,id))
             enddo
           else
             do iq = 1, nq
               dfdt(iq,ip-5) = inter_d3_dx(nwad_f(iq,id))
             enddo
           endif
         enddo
c
c        Interpolate the second order partial derivatives
c
c        First do d2fdx2
c
         do ip = 1, nvar
           id = ip*(ip+1)*(ip+2)/6
           id2 = D2_RA_RA
           if (ip.eq.2) then
             id2 = D2_RB_RB
           elseif (ip.eq.3) then
             id2 = D2_GAA_GAA
           elseif (ip.eq.4) then
             id2 = D2_GAB_GAB
           elseif (ip.eq.5) then
             id2 = D2_GBB_GBB
           elseif (ip.eq.6) then
             id2 = D2_TA_TA
           elseif (ip.eq.7) then
             id2 = D2_TB_TB
           endif
           if (ip.le.2) then
             do iq = 1, nq
               dfdr2(iq,id2) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               dfdg2(iq,id2) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           else
             do iq = 1, nq
               dfdt2(iq,id2) = inter_d3_dx2(nwad_f(iq,id))
             enddo
           endif
         enddo
c
c        Second do d2fdxy
c
         do ip = 1, nvar
           do ir = 1, ip-1
             ip3  = ip*(ip+1)*(ip+2)/6
             ir3  = ir*(ir+1)*(ip+2)/6
             ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             id2 = D2_RA_RA
             if (ip.eq.2.and.ir.eq.1) then
               id2 = D2_RA_RB
             elseif (ip.eq.3.and.ir.eq.1) then
               id2 = D2_RA_GAA
             elseif (ip.eq.3.and.ir.eq.2) then
               id2 = D2_RB_GAA
             elseif (ip.eq.4.and.ir.eq.1) then
               id2 = D2_RA_GAB
             elseif (ip.eq.4.and.ir.eq.2) then
               id2 = D2_RB_GAB
             elseif (ip.eq.4.and.ir.eq.3) then
               id2 = D2_GAA_GAB
             elseif (ip.eq.5.and.ir.eq.1) then
               id2 = D2_RA_GBB
             elseif (ip.eq.5.and.ir.eq.2) then
               id2 = D2_RB_GBB
             elseif (ip.eq.5.and.ir.eq.3) then
               id2 = D2_GAA_GBB
             elseif (ip.eq.5.and.ir.eq.4) then
               id2 = D2_GAB_GBB
             elseif (ip.eq.6.and.ir.eq.1) then
               id2 = D2_RA_TA
             elseif (ip.eq.6.and.ir.eq.2) then
               id2 = D2_RB_TA
             elseif (ip.eq.6.and.ir.eq.3) then
               id2 = D2_GAA_TA
             elseif (ip.eq.6.and.ir.eq.4) then
               id2 = D2_GAB_TA
             elseif (ip.eq.6.and.ir.eq.5) then
               id2 = D2_GBB_TA
             elseif (ip.eq.7.and.ir.eq.1) then
               id2 = D2_RA_TB
             elseif (ip.eq.7.and.ir.eq.2) then
               id2 = D2_RB_TB
             elseif (ip.eq.7.and.ir.eq.3) then
               id2 = D2_GAA_TB
             elseif (ip.eq.7.and.ir.eq.4) then
               id2 = D2_GAB_TB
             elseif (ip.eq.7.and.ir.eq.5) then
               id2 = D2_GBB_TB
             elseif (ip.eq.7.and.ir.eq.6) then
               id2 = D2_TA_TB
             endif
             if (ip.le.2) then
               do iq = 1, nq
                 dfdr2(iq,id2) = inter_d3_dxy(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 dfdg2(iq,id2) = inter_d3_dxy(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             else
               do iq = 1, nq
                 dfdt2(iq,id2) = inter_d3_dxy(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Interpolate the third order partial derivatives
c
c        First do d3fdx3
c
         do ip = 1, nvar
           id = ip*(ip+1)*(ip+2)/6
           id2 = D3_RA_RA_RA
           if (ip.eq.2) then
             id2 = D3_RB_RB_RB
           elseif (ip.eq.3) then
             id2 = D3_GAA_GAA_GAA
           elseif (ip.eq.4) then
             id2 = D3_GAB_GAB_GAB
           elseif (ip.eq.5) then
             id2 = D3_GBB_GBB_GBB
           elseif (ip.eq.6) then
             id2 = D3_TA_TA_TA
           elseif (ip.eq.7) then
             id2 = D3_TB_TB_TB
           endif
           if (ip.le.2) then
             do iq = 1, nq
               dfdr3(iq,id2) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           elseif (ip.le.5) then
             do iq = 1, nq
               dfdg3(iq,id2) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           else
             do iq = 1, nq
               dfdt3(iq,id2) = inter_d3_dx3(nwad_f(iq,id))
             enddo
           endif
         enddo
c
c        Second do d3fdx2y
c
         do ip = 1, nvar
           do ir = 1, ip-1
             ip3  = ip*(ip+1)*(ip+2)/6
             ir3  = ir*(ir+1)*(ip+2)/6
             ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             id2 = D3_RA_RA_RA
             if (ip.eq.2.and.ir.eq.1) then
               id2 = D3_RA_RB_RB
             elseif (ip.eq.3.and.ir.eq.1) then
               id2 = D3_RA_GAA_GAA
             elseif (ip.eq.3.and.ir.eq.2) then
               id2 = D3_RB_GAA_GAA
             elseif (ip.eq.4.and.ir.eq.1) then
               id2 = D3_RA_GAB_GAB
             elseif (ip.eq.4.and.ir.eq.2) then
               id2 = D3_RB_GAB_GAB
             elseif (ip.eq.4.and.ir.eq.3) then
               id2 = D3_GAA_GAB_GAB
             elseif (ip.eq.5.and.ir.eq.1) then
               id2 = D3_RA_GBB_GBB
             elseif (ip.eq.5.and.ir.eq.2) then
               id2 = D3_RB_GBB_GBB
             elseif (ip.eq.5.and.ir.eq.3) then
               id2 = D3_GAA_GBB_GBB
             elseif (ip.eq.5.and.ir.eq.4) then
               id2 = D3_GAB_GBB_GBB
             elseif (ip.eq.6.and.ir.eq.1) then
               id2 = D3_RA_TA_TA
             elseif (ip.eq.6.and.ir.eq.2) then
               id2 = D3_RB_TA_TA
             elseif (ip.eq.6.and.ir.eq.3) then
               id2 = D3_GAA_TA_TA
             elseif (ip.eq.6.and.ir.eq.4) then
               id2 = D3_GAB_TA_TA
             elseif (ip.eq.6.and.ir.eq.5) then
               id2 = D3_GBB_TA_TA
             elseif (ip.eq.7.and.ir.eq.1) then
               id2 = D3_RA_TB_TB
             elseif (ip.eq.7.and.ir.eq.2) then
               id2 = D3_RB_TB_TB
             elseif (ip.eq.7.and.ir.eq.3) then
               id2 = D3_GAA_TB_TB
             elseif (ip.eq.7.and.ir.eq.4) then
               id2 = D3_GAB_TB_TB
             elseif (ip.eq.7.and.ir.eq.5) then
               id2 = D3_GBB_TB_TB
             elseif (ip.eq.7.and.ir.eq.6) then
               id2 = D3_TA_TB_TB
             endif
             if (ip.le.2) then
               do iq = 1, nq
                 dfdr3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 dfdg3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             else
               do iq = 1, nq
                 dfdt3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ip3),
     +                           nwad_f(iq,ip2r),nwad_f(iq,ipr2),
     +                           nwad_f(iq,ir3))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Third do d3fdxy2
c
         do ip = 1, nvar
           do ir = 1, ip-1
             ip3  = ip*(ip+1)*(ip+2)/6
             ir3  = ir*(ir+1)*(ip+2)/6
             ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
             ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
             id2 = D3_RA_RA_RA
             if (ip.eq.2.and.ir.eq.1) then
               id2 = D3_RA_RA_RB
             elseif (ip.eq.3.and.ir.eq.1) then
               id2 = D3_RA_RA_GAA
             elseif (ip.eq.3.and.ir.eq.2) then
               id2 = D3_RB_RB_GAA
             elseif (ip.eq.4.and.ir.eq.1) then
               id2 = D3_RA_RA_GAB
             elseif (ip.eq.4.and.ir.eq.2) then
               id2 = D3_RB_RB_GAB
             elseif (ip.eq.4.and.ir.eq.3) then
               id2 = D3_GAA_GAA_GAB
             elseif (ip.eq.5.and.ir.eq.1) then
               id2 = D3_RA_RA_GBB
             elseif (ip.eq.5.and.ir.eq.2) then
               id2 = D3_RB_RB_GBB
             elseif (ip.eq.5.and.ir.eq.3) then
               id2 = D3_GAA_GAA_GBB
             elseif (ip.eq.5.and.ir.eq.4) then
               id2 = D3_GAB_GAB_GBB
             elseif (ip.eq.6.and.ir.eq.1) then
               id2 = D3_RA_RA_TA
             elseif (ip.eq.6.and.ir.eq.2) then
               id2 = D3_RB_RB_TA
             elseif (ip.eq.6.and.ir.eq.3) then
               id2 = D3_GAA_GAA_TA
             elseif (ip.eq.6.and.ir.eq.4) then
               id2 = D3_GAB_GAB_TA
             elseif (ip.eq.6.and.ir.eq.5) then
               id2 = D3_GBB_GBB_TA
             elseif (ip.eq.7.and.ir.eq.1) then
               id2 = D3_RA_RA_TB
             elseif (ip.eq.7.and.ir.eq.2) then
               id2 = D3_RB_RB_TB
             elseif (ip.eq.7.and.ir.eq.3) then
               id2 = D3_GAA_GAA_TB
             elseif (ip.eq.7.and.ir.eq.4) then
               id2 = D3_GAB_GAB_TB
             elseif (ip.eq.7.and.ir.eq.5) then
               id2 = D3_GBB_GBB_TB
             elseif (ip.eq.7.and.ir.eq.6) then
               id2 = D3_TA_TA_TB
             endif
             if (ip.le.2) then
               do iq = 1, nq
                 dfdr3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ir3),
     +                           nwad_f(iq,ipr2),nwad_f(iq,ip2r),
     +                           nwad_f(iq,ip3))
               enddo
             elseif (ip.le.5) then
               do iq = 1, nq
                 dfdg3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ir3),
     +                           nwad_f(iq,ipr2),nwad_f(iq,ip2r),
     +                           nwad_f(iq,ip3))
               enddo
             else
               do iq = 1, nq
                 dfdt3(iq,id2) = inter_d3_dx2y(nwad_f(iq,ir3),
     +                           nwad_f(iq,ipr2),nwad_f(iq,ip2r),
     +                           nwad_f(iq,ip3))
               enddo
             endif
           enddo ! ir
         enddo ! ip
c
c        Fourth do d3fdxyz
c
         do ip = 1, nvar
           do ir = 1, ip-1
             do it = 1, ir-1
               ip3  = ip*(ip+1)*(ip+2)/6
               ir3  = ir*(ir+1)*(ip+2)/6
               it3  = it*(it+1)*(it+2)/6
               ip2r = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+ir
               ip2t = (ip-1)*ip*(ip+1)/6+ip*(ip-1)/2+it
               ipr2 = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+ir
               ir2t = (ir-1)*ir*(ir+1)/6+ir*(ir-1)/2+it
               ipt2 = (ip-1)*ip*(ip+1)/6+it*(it-1)/2+it
               irt2 = (ir-1)*ir*(ir+1)/6+it*(it-1)/2+it
               iprt = (ip-1)*ip*(ip+1)/6+ir*(ir-1)/2+it
               id2 = D3_RA_RA_RA
               if (ip.eq.3.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_RB_GAA
               elseif (ip.eq.4.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_RB_GAB
               elseif (ip.eq.4.and.ir.eq.3.and.it.eq.1) then
                 id2 = D3_RA_GAA_GAB
               elseif (ip.eq.4.and.ir.eq.3.and.it.eq.2) then
                 id2 = D3_RB_GAA_GAB
               elseif (ip.eq.5.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_RB_GBB
               elseif (ip.eq.5.and.ir.eq.3.and.it.eq.1) then
                 id2 = D3_RA_GAA_GBB
               elseif (ip.eq.5.and.ir.eq.3.and.it.eq.2) then
                 id2 = D3_RB_GAA_GBB
               elseif (ip.eq.5.and.ir.eq.4.and.it.eq.1) then
                 id2 = D3_RA_GAB_GBB
               elseif (ip.eq.5.and.ir.eq.4.and.it.eq.2) then
                 id2 = D3_RB_GAB_GBB
               elseif (ip.eq.5.and.ir.eq.4.and.it.eq.3) then
                 id2 = D3_GAA_GAB_GBB
               elseif (ip.eq.6.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_RB_TA
               elseif (ip.eq.6.and.ir.eq.3.and.it.eq.1) then
                 id2 = D3_RA_GAA_TA
               elseif (ip.eq.6.and.ir.eq.3.and.it.eq.2) then
                 id2 = D3_RB_GAA_TA
               elseif (ip.eq.6.and.ir.eq.4.and.it.eq.1) then
                 id2 = D3_RA_GAB_TA
               elseif (ip.eq.6.and.ir.eq.4.and.it.eq.2) then
                 id2 = D3_RB_GAB_TA
               elseif (ip.eq.6.and.ir.eq.4.and.it.eq.3) then
                 id2 = D3_GAA_GAB_TA
               elseif (ip.eq.6.and.ir.eq.5.and.it.eq.1) then
                 id2 = D3_RA_GBB_TA
               elseif (ip.eq.6.and.ir.eq.5.and.it.eq.2) then
                 id2 = D3_RB_GBB_TA
               elseif (ip.eq.6.and.ir.eq.5.and.it.eq.3) then
                 id2 = D3_GAA_GBB_TA
               elseif (ip.eq.6.and.ir.eq.5.and.it.eq.4) then
                 id2 = D3_GAB_GBB_TA
               elseif (ip.eq.7.and.ir.eq.2.and.it.eq.1) then
                 id2 = D3_RA_RB_TB
               elseif (ip.eq.7.and.ir.eq.3.and.it.eq.1) then
                 id2 = D3_RA_GAA_TB
               elseif (ip.eq.7.and.ir.eq.3.and.it.eq.2) then
                 id2 = D3_RB_GAA_TB
               elseif (ip.eq.7.and.ir.eq.4.and.it.eq.1) then
                 id2 = D3_RA_GAB_TB
               elseif (ip.eq.7.and.ir.eq.4.and.it.eq.2) then
                 id2 = D3_RB_GAB_TB
               elseif (ip.eq.7.and.ir.eq.4.and.it.eq.3) then
                 id2 = D3_GAA_GAB_TB
               elseif (ip.eq.7.and.ir.eq.5.and.it.eq.1) then
                 id2 = D3_RA_GBB_TB
               elseif (ip.eq.7.and.ir.eq.5.and.it.eq.2) then
                 id2 = D3_RB_GBB_TB
               elseif (ip.eq.7.and.ir.eq.5.and.it.eq.3) then
                 id2 = D3_GAA_GBB_TB
               elseif (ip.eq.7.and.ir.eq.5.and.it.eq.4) then
                 id2 = D3_GAB_GBB_TB
               elseif (ip.eq.7.and.ir.eq.6.and.it.eq.1) then
                 id2 = D3_RA_TA_TB
               elseif (ip.eq.7.and.ir.eq.6.and.it.eq.2) then
                 id2 = D3_RB_TA_TB
               elseif (ip.eq.7.and.ir.eq.6.and.it.eq.3) then
                 id2 = D3_GAA_TA_TB
               elseif (ip.eq.7.and.ir.eq.6.and.it.eq.4) then
                 id2 = D3_GAB_TA_TB
               elseif (ip.eq.7.and.ir.eq.6.and.it.eq.5) then
                 id2 = D3_GBB_TA_TB
               endif
               if (ip.le.2) then
                 do iq = 1, nq
                   dfdr3(iq,id2)
     +             = inter_d3_dxyz(nwad_f(iq,ip3),nwad_f(iq,ir3),
     +               nwad_f(iq,it3),nwad_f(iq,ip2r),nwad_f(iq,ip2t),
     +               nwad_f(iq,ipr2),nwad_f(iq,ir2t),nwad_f(iq,ipt2),
     +               nwad_f(iq,irt2),nwad_f(iq,iprt))
                 enddo
               elseif (ip.le.5) then
                 do iq = 1, nq
                   dfdg3(iq,id2)
     +             = inter_d3_dxyz(nwad_f(iq,ip3),nwad_f(iq,ir3),
     +               nwad_f(iq,it3),nwad_f(iq,ip2r),nwad_f(iq,ip2t),
     +               nwad_f(iq,ipr2),nwad_f(iq,ir2t),nwad_f(iq,ipt2),
     +               nwad_f(iq,irt2),nwad_f(iq,iprt))
                 enddo
               else
                 do iq = 1, nq
                   dfdt3(iq,id2)
     +             = inter_d3_dxyz(nwad_f(iq,ip3),nwad_f(iq,ir3),
     +               nwad_f(iq,it3),nwad_f(iq,ip2r),nwad_f(iq,ip2t),
     +               nwad_f(iq,ipr2),nwad_f(iq,ir2t),nwad_f(iq,ipt2),
     +               nwad_f(iq,irt2),nwad_f(iq,iprt))
                 enddo
               endif
             enddo ! it
           enddo ! ir
         enddo ! ip
c
       else
         call errquit("nwxc_eval_df3_driver: insane ipol",ipol,UERR)
       endif
C      
       end
C>
C> \brief Evaluate the exchange-correlation energy and its 1st, 2nd
C> and 3rd partial derivatives
C>
      subroutine nwxc_eval_df3_doit(ipol,nq,rho,gamma,tau,f,
     +                               dfdr,dfdr2,dfdr3,
     +                               dfdg,dfdg2,dfdg3,
     +                               dfdt,dfdt2,dfdt3)
      use nwad3
      implicit none
#include "nwxcP.fh"
#include "nwxc_param.fh"
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of points
C
      type(nwad_dble)::rho(nq,ipol)     !< [Input] Density
      type(nwad_dble)::gamma(nq,ipol+1) !< [Input] |Density gradient|^2
      type(nwad_dble)::tau(nq,ipol)     !< [Input] Kinetic energy 
                                         !< density
C
      type(nwad_dble)::f(nq)           !< [Output] Energy
      double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
      double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
      double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
c
      double precision dfdr2(nq,*)     !< [Output] 2nd derivative wrt
                                       !< rho
      double precision dfdg2(nq,*)     !< [Output] 2nd derivative wrt
                                       !< rho and gamma, and gamma
      double precision dfdt2(nq,*)     !< [Output] 2nd derivative wrt
                                       !< rho and tau, gamma and tau,
                                       !< and tau
c
      double precision dfdr3(nq,*)     !< [Output] 3rd derivative wrt
                                       !< rho
      double precision dfdg3(nq,*)     !< [Output] 3rd derivative wrt
                                       !< rho and gamma, and gamma
      double precision dfdt3(nq,*)     !< [Output] 3rd derivative wrt
                                       !< rho, gamm and tau
C
      integer iq !< Counter over points
      integer ip !< Counter over spin channels
      integer it !< Counter of functional terms
      integer n  !< Counter
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
      integer nd2r(2) !< The number of 2nd partial derivatives wrt rho
                      !< as a function of ipol
      integer nd2g(2) !< The number of 2nd partial derivatives wrt
                      !< gamma as a function of ipol
      integer nd2t(2) !< The number of 2nd partial derivatives wrt tau
                      !< as a function of ipol
      integer nd3r(2) !< The number of 3rd partial derivatives wrt rho
                      !< as a function of ipol
      integer nd3g(2) !< The number of 3rd partial derivatives wrt
                      !< gamma as a function of ipol
      integer nd3t(2) !< The number of 3rd partial derivatives wrt rho
                      !< and tau as a function of ipol
      integer nd3gt(2)!< The number of 3rd partial derivatives wrt rho,
                      !< gamma and tau as a function of ipol
C
      external nwxc_k_dirac_d3
      external nwxc_k_becke88_d3
      external nwxc_k_pbe96_d3
c
      integer max_param
      parameter (max_param = 50)
      double precision param(max_param)
      double precision Axlsda
      parameter (Axlsda = -0.9305257363491d0 )
C
      data nd1r / D1_RA,  D1_RB  /
      data nd1g / D1_GAA, D1_GBB /
      data nd1t / D1_TA,  D1_TB  /
C
      data nd2r / D2_RA_RA,   D2_RB_RB   /
c      data nd2g / D2_GAA_GAA, D2_GBB_GBB /
      data nd2g / D2_GAA_GBB, D2_GBB_GBB /
      data nd2t / D2_TA_TA,   D2_TB_TB   /
C
      data nd3r  / D3_RA_RA_RA,    D3_RB_RB_RB    /
      data nd3g  / D3_GAA_GBB_GBB, D3_GBB_GBB_GBB /
      data nd3t  / D3_TA_TA_TA,    D3_TB_TB_TB    /
      data nd3gt / D3_GAA_TB_TB,   D3_GBB_TB_TB   /
C
      do iq = 1, nq
        f(iq) = 0.0d0
      enddo
c      do ip = 1, nd1r(ipol)
c        do iq = 1, nq
c          dfdr(iq,ip) = 0.0d0
c        enddo
c      enddo
c      do ip = 1, nd2r(ipol)
c        do iq = 1, nq
c          dfdr2(iq,ip) = 0.0d0
c        enddo
c      enddo
c      do ip = 1, nd3r(ipol)
c        do iq = 1, nq
c          dfdr3(iq,ip) = 0.0d0
c        enddo
c      enddo
c      if (nwxc_is_gga().or.nwxc_is_mgga()) then
c        do ip = 1, nd1g(ipol)
c          do iq = 1, nq
c            dfdg(iq,ip) = 0.0d0
c          enddo
c        enddo
c        do ip = 1, nd2g(ipol)
c          do iq = 1, nq
c            dfdg2(iq,ip) = 0.0d0
c          enddo
c        enddo
c        do ip = 1, nd3g(ipol)
c          do iq = 1, nq
c            dfdg3(iq,ip) = 0.0d0
c          enddo
c        enddo
c      endif
c      if (nwxc_is_mgga()) then
c        do ip = 1, nd1t(ipol)
c          do iq = 1, nq
c            dfdt(iq,ip) = 0.0d0
c          enddo
c        enddo
c        do ip = 1, nd2t(ipol)
c          do iq = 1, nq
c            dfdt2(iq,ip) = 0.0d0
c          enddo
c        enddo
c        do ip = 1, nd3t(ipol)
c          do iq = 1, nq
c            dfdt3(iq,ip) = 0.0d0
c          enddo
c        enddo
c      endif
C
      do it = 1, nwxc_num_f
        select case (nwxc_ids(it))
          case (NWXCP_X_SLATER)
            call nwxc_x_dirac_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,f)
          case (NWXCP_X_B86b)
            call nwxc_x_b86b_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_B88)
            call nwxc_x_b88_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,f)
c          case (NWXCP_X_B97)
c            param(1) = 2.0d0
c            param(2) = 0.80940d+00
c            param(3) = 0.50730d+00
c            param(4) = 0.74810d+00
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
c          case (NWXCP_X_B97_1)
c            param(1) = 2.0d0
c            param(2) = 0.789518d+00
c            param(3) = 0.573805d+00
c            param(4) = 0.660975d+00
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
c          case (NWXCP_X_B97_2)
c            param(1) = 2.0d0
c            param(2) = 0.827642D+00
c            param(3) = 0.478400D-01
c            param(4) = 0.176125D+01
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
c          case (NWXCP_X_B97_3)
c            param(1) = 4.0d0
c            param(2) = 7.334648D-01
c            param(3) = 2.925270D-01
c            param(4) = 3.338789D+00
c            param(5) =-1.051158D+01
c            param(6) = 1.060907D+01
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
c          case (NWXCP_X_B97_D)
c            param(1) = 2.0d0
c            param(2) = 1.086620d+0
c            param(3) =-0.521270d+00
c            param(4) = 3.254290d+00
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
c          case (NWXCP_X_B97_G)
c            param(1) = 2.0d0
c            param(2) = 1.1068d0
c            param(3) =-0.8765d0
c            param(4) = 4.2639d0
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
c          case (NWXCP_X_B98)
c            param(1) = 2.0d0
c            param(2) = 0.790194d00
c            param(3) = 0.400271d00
c            param(4) = 0.832857d00
c            call nwxc_x_b97_d2(param,nwxc_rho_tol,ipol,nq,
c    +                          nwxc_wghts(it),rho,gamma,f,dfdr,dfdr2,
c    +                          dfdg,dfdg2)
          case (NWXCP_X_BNL)
            param(1) = nwxc_cam_gamma
            call nwxc_x_bnl_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,f)
          case (NWXCP_X_CAMB88)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camb88_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMLSD)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camlsd_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,f)
          case (NWXCP_X_CAMPBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            call nwxc_x_campbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMREVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            call nwxc_x_campbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMRPBE)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            call nwxc_x_camrpbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                             nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMS12G)
            param(1) = 1.03323556d0
            param(2) = 0.757d0
            param(3) = 0.00417251d0
            param(4) = 0.00115216d0
            param(5) = 0.00706184d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            call nwxc_x_cams12_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_CAMS12H)
            param(1) = 1.02149642d0
            param(2) = 0.757d0
            param(3) = 0.00825905d0
            param(4) = 0.00235804d0
            param(5) = 0.00654977d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            call nwxc_x_cams12_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_DLDF)
            call nwxc_x_dldf_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,tau,f)
          case (NWXCP_X_FT97)
            call nwxc_x_ft97_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_GILL)
            call nwxc_x_gill_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_X_HCTH)
            param(1) = 4.0d0
            param(2) = 0.109320d+01
            param(3) =-0.744056d+00
            param(4) = 0.559920d+01
            param(5) =-0.678549d+01
            param(6) = 0.449357d+01
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH120)
            param(1) = 4.0d0
            param(2) = 1.09163d0
            param(3) =-0.74720d0
            param(4) = 5.07830d0
            param(5) =-4.10750d0
            param(6) = 1.17170d0
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH147)
            param(1) = 4.0d0
            param(2) = 1.09025d0
            param(3) =-0.79920d0
            param(4) = 5.57210d0
            param(5) =-5.86760d0
            param(6) = 3.04540d0
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH407)
            param(1) = 4.0d0
            param(2) = 1.08184d0
            param(3) =-0.5183d0
            param(4) = 3.4256d0
            param(5) =-2.6290d0
            param(6) = 2.2886d0
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH407P)
            param(1) = 4.0d0
            param(2) = 1.08018D0
            param(3) =-0.4117D0
            param(4) = 2.4368D0
            param(5) = 1.3890D0
            param(6) =-1.3529D0
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTH_A)
            param(1) = 2.0d0
            param(2) = 0.109878d+01
            param(3) =-0.251173d+01
            param(4) = 0.156233d-01
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_X_HCTHP14)
            param(1) = 4.0d0
            param(2) = 0.103161d+01
            param(3) =-0.360781d+00
            param(4) = 0.351994d+01
            param(5) =-0.495944d+01
            param(6) = 0.241165d+01
            call nwxc_x_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_M05)
             param( 1) =    0.08151d0
             param( 2) =   -0.43956d0
             param( 3) =   -3.22422d0
             param( 4) =    2.01819d0
             param( 5) =    8.79431d0
             param( 6) =   -0.00295d0
             param( 7) =    9.82029d0
             param( 8) =   -4.82351d0
             param( 9) =  -48.17574d0
             param(10) =    3.64802d0
             param(11) =   34.02248d0
             call nwxc_x_m05_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M05_2X)
             param( 1) =   -0.56833d0
             param( 2) =   -1.30057d0
             param( 3) =    5.50070d0
             param( 4) =    9.06402d0
             param( 5) =  -32.21075d0
             param( 6) =  -23.73298d0
             param( 7) =   70.22996d0
             param( 8) =   29.88614d0
             param( 9) =  -60.25778d0
             param(10) =  -13.22205d0
             param(11) =   15.23694d0
             call nwxc_x_m05_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06)
             param( 1) =  1.422057D-01*Axlsda
             param( 2) =  7.370319D-04*Axlsda
             param( 3) = -1.601373D-02*Axlsda
             param( 4) =  0.000000D+00
             param( 5) =  0.000000D+00
             param( 6) =  0.000000D+00
             param( 7) =  5.877943D-01
             param( 8) = -1.371776D-01
             param( 9) =  2.682367D-01
             param(10) = -2.515898D+00
             param(11) = -2.978892D+00
             param(12) =  8.710679D+00
             param(13) =  1.688195D+01
             param(14) = -4.489724D+00
             param(15) = -3.299983D+01
             param(16) = -1.449050D+01
             param(17) =  2.043747D+01
             param(18) =  1.256504D+01
             call nwxc_x_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06_HF)
             param( 1) = -1.179732D-01*Axlsda
             param( 2) = -2.500000D-03*Axlsda
             param( 3) = -1.180065D-02*Axlsda
             param( 4) =  0.000000D+00
             param( 5) =  0.000000D+00
             param( 6) =  0.000000D+00
             param( 7) =  1.179732D-01
             param( 8) = -1.066708D+00
             param( 9) = -1.462405D-01
             param(10) =  7.481848D+00
             param(11) =  3.776679D+00
             param(12) = -4.436118D+01
             param(13) = -1.830962D+01
             param(14) =  1.003903D+02
             param(15) =  3.864360D+01
             param(16) = -9.806018D+01
             param(17) = -2.557716D+01
             param(18) =  3.590404D+01
             call nwxc_x_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06_L)
             param( 1) =  6.012244D-01*Axlsda
             param( 2) =  4.748822D-03*Axlsda
             param( 3) = -8.635108D-03*Axlsda
             param( 4) = -9.308062D-06*Axlsda
             param( 5) =  4.482811D-05*Axlsda
             param( 6) =  0.000000D+00
             param( 7) =  3.987756D-01
             param( 8) =  2.548219D-01
             param( 9) =  3.923994D-01
             param(10) = -2.103655D+00
             param(11) = -6.302147D+00
             param(12) =  1.097615D+01
             param(13) =  3.097273D+01
             param(14) = -2.318489D+01
             param(15) = -5.673480D+01
             param(16) =  2.160364D+01
             param(17) =  3.421814D+01
             param(18) = -9.049762D+00
             call nwxc_x_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M06_2X)
             param( 1) =  4.600000D-01
             param( 2) = -2.206052D-01
             param( 3) = -9.431788D-02
             param( 4) =  2.164494D+00
             param( 5) = -2.556466D+00
             param( 6) = -1.422133D+01
             param( 7) =  1.555044D+01
             param( 8) =  3.598078D+01
             param( 9) = -2.722754D+01
             param(10) = -3.924093D+01
             param(11) =  1.522808D+01
             param(12) =  1.522227D+01
             call nwxc_x_m06_2x_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M08_HX)
c            parameters A
             param(01) =  1.3340172D+00
             param(02) = -9.4751087D+00
             param(03) = -1.2541893D+01
             param(04) =  9.1369974D+00
             param(05) =  3.4717204D+01
             param(06) =  5.8831807D+01
             param(07) =  7.1369574D+01
             param(08) =  2.3312961D+01
             param(09) =  4.8314679D+00
             param(10) = -6.5044167D+00
             param(11) = -1.4058265D+01
             param(12) =  1.2880570D+01
c            parameters B
             param(13) = -8.5631823D-01
             param(14) =  9.2810354D+00
             param(15) =  1.2260749D+01
             param(16) = -5.5189665D+00
             param(17) = -3.5534989D+01
             param(18) = -8.2049996D+01
             param(19) = -6.8586558D+01
             param(20) =  3.6085694D+01
             param(21) = -9.3740983D+00
             param(22) = -5.9731688D+01
             param(23) =  1.6587868D+01
             param(24) =  1.3993203D+01
c            parameters C and D
             do n = 25, 48
               param(n) = 0.0d0
             enddo
             call nwxc_x_m08_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M08_SO)
c            parameters A
             param(01) = -3.4888428D-01
             param(02) = -5.8157416D+00
             param(03) =  3.7550810D+01
             param(04) =  6.3727406D+01
             param(05) = -5.3742313D+01
             param(06) = -9.8595529D+01
             param(07) =  1.6282216D+01
             param(08) =  1.7513468D+01
             param(09) = -6.7627553D+00
             param(10) =  1.1106658D+01
             param(11) =  1.5663545D+00
             param(12) =  8.7603470D+00
c            parameters B
             param(13) =  7.8098428D-01
             param(14) =  5.4538178D+00
             param(15) = -3.7853348D+01
             param(16) = -6.2295080D+01
             param(17) =  4.6713254D+01
             param(18) =  8.7321376D+01
             param(19) =  1.6053446D+01
             param(20) =  2.0126920D+01
             param(21) = -4.0343695D+01
             param(22) = -5.8577565D+01
             param(23) =  2.0890272D+01
             param(24) =  1.0946903D+01
c            parameters C and D
             do n = 25, 48
               param(n) = 0.0d0
             enddo
             call nwxc_x_m08_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M11)
c            parameters A
             param(01) = -0.18399900D+00
             param(02) = -1.39046703D+01
             param(03) =  1.18206837D+01
             param(04) =  3.10098465D+01
             param(05) = -5.19625696D+01
             param(06) =  1.55750312D+01
             param(07) = -6.94775730D+00
             param(08) = -1.58465014D+02
             param(09) = -1.48447565D+00
             param(10) =  5.51042124D+01
             param(11) = -1.34714184D+01
             param(12) =  0.000000D+00
c            parameters B
             param(13) =  0.75599900D+00
             param(14) =  1.37137944D+01
             param(15) = -1.27998304D+01
             param(16) = -2.93428814D+01
             param(17) =  5.91075674D+01
             param(18) = -2.27604866D+01
             param(19) = -1.02769340D+01
             param(20) =  1.64752731D+02
             param(21) =  1.85349258D+01
             param(22) = -5.56825639D+01
             param(23) =  7.47980859D+00
             param(24) =  0.000000D+00
c            parameters C and D
             do n = 25, 48
               param(n) = 0.0d0
             enddo
             call nwxc_x_m11_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_M11_L)
c            parameters A
             param(01) =  8.121131D-01
             param(02) =  1.738124D+01
             param(03) =  1.154007D+00
             param(04) =  6.869556D+01
             param(05) =  1.016864D+02
             param(06) = -5.887467D+00
             param(07) =  4.517409D+01
             param(08) = -2.773149D+00
             param(09) = -2.617211D+01
             param(10) =  0.000000D+00
             param(11) =  0.000000D+00
             param(12) =  0.000000D+00
c            parameters B
             param(13) =  1.878869D-01
             param(14) = -1.653877D+01
             param(15) =  6.755753D-01
             param(16) = -7.567572D+01
             param(17) = -1.040272D+02
             param(18) =  1.831853D+01
             param(19) = -5.573352D+01
             param(20) = -3.520210D+00
             param(21) =  3.724276D+01
             param(22) =  0.000000D+00
             param(23) =  0.000000D+00
             param(24) =  0.000000D+00
c            parameters C
             param(25) = -4.386615D-01
             param(26) = -1.214016D+02
             param(27) = -1.393573D+02
             param(28) = -2.046649D+00
             param(29) =  2.804098D+01
             param(30) = -1.312258D+01
             param(31) = -6.361819D+00
             param(32) = -8.055758D-01
             param(33) =  3.736551D+00
             param(34) =  0.000000D+00
             param(35) =  0.000000D+00
             param(36) =  0.000000D+00
c            parameters D
             param(37) =  1.438662D+00
             param(38) =  1.209465D+02
             param(39) =  1.328252D+02
             param(40) =  1.296355D+01
             param(41) =  5.854866D+00
             param(42) = -3.378162D+00
             param(43) = -4.423393D+01
             param(44) =  6.844475D+00
             param(45) =  1.949541D+01
             param(46) =  0.000000D+00
             param(47) =  0.000000D+00
             param(48) =  0.000000D+00
             call nwxc_x_m11_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_X_MPW91)
             param(1) = 3.72d0
             param(2) = 0.00426D0
             call nwxc_x_pw91_d3(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_OPT)
             call nwxc_x_opt_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
           case (NWXCP_X_PW6B95)
             param(1) = 0.00538d0
             param(2) = 1.7382d0
             param(3) = 3.8901d0
             call nwxc_x_pw6_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_PWB6K)
             param(1) = 0.00539d0
             param(2) = 1.7077d0
             param(3) = 4.0876d0
             call nwxc_x_pw6_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_PW91)
             param(1) = 4.0d0
             param(2) = 0.0042D0
             call nwxc_x_pw91_d3(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_PBE)
             param(1) = 0.8040d0
             param(2) = 0.2195149727645171d0
             call nwxc_x_pbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_PBESOL)
             param(1) = 0.8040d0
             param(2) = 10.0d0/81.0d0
             call nwxc_x_pbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_REVPBE)
             param(1) = 1.245d0
             param(2) = 0.2195149727645171d0
             call nwxc_x_pbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_RPBE)
             call nwxc_x_rpbe_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,gamma,f)
           case (NWXCP_X_PKZB)
             call nwxc_x_pkzb99_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                             rho,gamma,tau,f)
           case (NWXCP_X_S12G)
             param(1) = 1.03842032d0
             param(2) = 0.757d0
             param(3) = 0.00403198d0
             param(4) = 0.00104596d0
             param(5) = 0.00594635d0
             call nwxc_x_s12_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_S12H)
             param(1) = 1.02543951d0
             param(2) = 0.757d0
             param(3) = 0.00761554d0
             param(4) = 0.00211063d0
             param(5) = 0.00604672d0
             call nwxc_x_s12_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_SOGGA)
             param(1)  =  0.5d0
             param(2)  =  0.276d0
             param(3)  =  0.0d0
             param(4)  =  0.0d0
             param(5)  =  0.0d0
             param(6)  =  0.0d0
             param(7)  =  0.5d0
             param(8)  =  0.276d0
             param(9)  =  0.0d0
             param(10) =  0.0d0
             param(11) =  0.0d0
             param(12) =  0.0d0
             call nwxc_x_sogga_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_SOGGA11)
             param(1)  =  0.5d0
             param(2)  = -2.95535d0
             param(3)  =  15.7974d0
             param(4)  = -91.1804d0
             param(5)  =  96.2030d0
             param(6)  =  0.18683d0
             param(7)  =  0.50000d0
             param(8)  =  3.50743d0
             param(9)  = -12.9523d0
             param(10) =  49.7870d0
             param(11) = -33.2545d0
             param(12) = -11.1396d0
             call nwxc_x_sogga_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_SOGGA11_X)
             param(1)  =  2.99250d-01
             param(2)  =  3.21638d+00
             param(3)  = -3.55605d+00
             param(4)  =  7.65852d+00
             param(5)  = -1.12830d+01
             param(6)  =  5.25813d+00
             param(7)  =  2.99250d-01
             param(8)  = -2.88595d+00
             param(9)  =  3.23617d+00
             param(10) = -2.45393d+00
             param(11) = -3.75495d+00
             param(12) =  3.96613d+00
             call nwxc_x_sogga_d3(param,nwxc_rho_tol,ipol,nq,
     +                            nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_SSB_D)
             call nwxc_x_ssbD_1_d3(nwxc_rho_tol,ipol,nq,
     +                             nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_TPSS)
             call nwxc_x_tpss03_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                             rho,gamma,tau,f)
          case (NWXCP_X_WPBE)
            param(1) = nwxc_cam_gamma
            call nwxc_x_wpbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_X_VS98)
             param(1) =  -9.800683d-01
             param(2) =  -3.556788d-03
             param(3) =   6.250326d-03
             param(4) =  -2.354518d-05
             param(5) =  -1.282732d-04
             param(6) =   3.574822d-04
             call nwxc_x_vs98_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_B95)
             param(1) = 0.0031d0
             param(2) = 0.038d0
             call nwxc_c_b95_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_B97)
             param(1) = 2.0d0
             param(2) = 0.17370d+00
             param(3) = 0.94540d+00
             param(4) = 0.23487d+01
             param(5) = 0.74710d+00
             param(6) =-0.24868d+01
             param(7) =-0.45961d+01
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_B97_1)
             param(1) = 2.0d0
             param(2) = 0.820011d-01
             param(3) = 0.955689d+00
             param(4) = 0.271681d+01
             param(5) = 0.788552d+00
             param(6) =-0.287103d+01
             param(7) =-0.547869d+01
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_B97_2)
             param(1) = 2.0d0
             param(2) = 0.585808D+00
             param(3) = 0.999849D+00
             param(4) =-0.691682D+00
             param(5) = 0.140626D+01
             param(6) = 0.394796D+00
             param(7) =-0.744060D+01
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_B97_3)
             param(1)  = 4.0d0
             param(2)  = 5.623649D-01
             param(3)  = 1.133830D+00
             param(4)  =-1.322980D+00
             param(5)  =-2.811967D+00
             param(6)  = 6.359191D+00
             param(7)  = 7.431302D+00
             param(8)  =-7.464002D+00
             param(9)  =-1.969342D+00
             param(10) = 1.827082D+00
             param(11) =-1.174423D+01
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_B97_D)
             param(1) = 2.0d0
             param(2) = 0.22340d+00
             param(3) = 0.690410d+00
             param(4) =-1.562080d+00
             param(5) = 6.302700d00
             param(6) = 1.942930d+0
             param(7) =-14.97120d+00
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_B97_G)
             param(1) = 2.0d0
             param(2) = 0.4883d0
             param(3) = 0.7961d0
             param(4) =-2.117d0
             param(5) = 5.7060d0
             param(6) = 2.3235d0
             param(7) =-14.9820d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_B98)
             param(1) = 2.0d0
             param(2) =-0.120163d00
             param(3) = 0.934715d00
             param(4) = 2.82332d0
             param(5) = 1.14105d0
             param(6) =-2.59412d0
             param(7) =-5.33398d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_DLDF)
             call nwxc_c_dldf_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,gamma,tau,f)
           case (NWXCP_C_FT97)
             call nwxc_c_ft97_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,gamma,f)
           case (NWXCP_C_HCTH)
             param(1)  = 4.0d0
             param(2)  = 0.222601d0
             param(3)  = 0.729974d0
             param(4)  =-3.38622d-002
             param(5)  = 3.352870d0
             param(6)  =-1.25170d-002
             param(7)  =-11.543d0
             param(8)  =-0.802496d0
             param(9)  = 8.085640d0
             param(10) = 1.553960d0
             param(11) =-4.478570d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_HCTH120)
             param(1)  = 4.0d0
             param(2)  = 0.48951d0
             param(3)  = 0.51473d0
             param(4)  =-0.26070d0
             param(5)  = 6.92980d0
             param(6)  = 0.43290d0
             param(7)  =-24.7070d0
             param(8)  =-1.99250d0
             param(9)  = 23.1100d0
             param(10) = 2.48530d0
             param(11) =-11.3230d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_HCTH147)
             param(1)  = 4.0d0
             param(2)  = 0.56258d0
             param(3)  = 0.54235d0
             param(4)  =-1.71000d-002
             param(5)  = 7.01460d0
             param(6)  =-1.30640d0
             param(7)  =-28.3820d0
             param(8)  = 1.05750d0
             param(9)  = 35.0330d0
             param(10) = 0.88540d0
             param(11) =-20.4280d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_HCTH407)
             param(1)  = 4.0d0
             param(2)  = 1.18777d0
             param(3)  = 0.58908d0
             param(4)  =-2.40290d0
             param(5)  = 4.42370d0
             param(6)  = 5.61740d0
             param(7)  =-19.2220d0
             param(8)  =-9.17920d0
             param(9)  = 42.5720d0
             param(10) = 6.24800d0
             param(11) =-42.0050d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_HCTH407P)
             param(1)  = 4.0d0
             param(2)  = 0.80302d0
             param(3)  = 0.73604d0
             param(4)  =-1.04790d0
             param(5)  = 3.02700d0
             param(6)  = 4.98070d0
             param(7)  =-10.0750d0
             param(8)  =-12.8900d0
             param(9)  = 20.6110d0
             param(10) = 9.64460d0
             param(11) =-29.4180d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_HCTH_A)
             param(1)  = 4.0d0
             param(2)  = 1.36823d-002
             param(3)  = 0.836897d0
             param(4)  = 0.268920d0
             param(5)  = 1.720510d0
             param(6)  =-0.550769d0
             param(7)  =-2.784980d0
             param(8)  = 1.039470d0
             param(9)  =-4.575040d0
             param(10) = 0.000000d0
             param(11) = 0.000000d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_HCTHP14)
             param(1)  = 4.0d0
             param(2)  = 2.82414d0     
             param(3)  = 8.21827d-002
             param(4)  = 3.18843d-002
             param(5)  = 4.56466d0     
             param(6)  =-1.78512d0     
             param(7)  =-13.5529d0     
             param(8)  = 2.39795d0     
             param(9)  = 13.3820d0     
             param(10) =-0.876909d0     
             param(11) =-3.174930d0
             call nwxc_c_b97_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_M05)
             param( 1) =   1.00000d0
             param( 2) =   3.78569d0
             param( 3) = -14.15261d0
             param( 4) =  -7.46589d0
             param( 5) =  17.94491d0
             param( 6) =   1.00000d0
             param( 7) =   3.77344d0
             param( 8) = -26.04463d0
             param( 9) =  30.69913d0
             param(10) =  -9.22695d0
             call nwxc_c_m05_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M05_2X)
             param( 1) =   1.00000d0
             param( 2) =   1.09297d0
             param( 3) =  -3.79171d0
             param( 4) =   2.82810d0
             param( 5) = -10.58909d0
             param( 6) =   1.00000d0
             param( 7) =  -3.05430d0
             param( 8) =   7.61854d0
             param( 9) =   1.47665d0
             param(10) = -11.92365d0
             call nwxc_c_m05_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M06)
 
             param( 1) =  -2.741539D+00
             param( 2) =  -6.720113D-01
             param( 3) =  -7.932688D-02
             param( 4) =   1.918681D-03
             param( 5) =  -2.032902D-03
             param( 6) =   0.000000D+00
             param( 7) =   4.905945D-01
             param( 8) =  -1.437348D-01
             param( 9) =   2.357824D-01
             param(10) =   1.871015D-03
             param(11) =  -3.788963D-03
             param(12) =   0.000000D+00
c
             param(13) =   3.741539D+00
             param(14) =   2.187098D+02
             param(15) =  -4.531252D+02
             param(16) =   2.936479D+02
             param(17) =  -6.287470D+01
             param(18) =   5.094055D-01
             param(19) =  -1.491085D+00
             param(20) =   1.723922D+01
             param(21) =  -3.859018D+01
             param(22) =   2.845044D+01
c
             call nwxc_c_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M06_HF)
c
             param( 1) =  -6.746338D-01
             param( 2) =  -1.534002D-01
             param( 3) =  -9.021521D-02
             param( 4) =  -1.292037D-03
             param( 5) =  -2.352983D-04
             param( 6) =   0.000000D+00
             param( 7) =   8.976746D-01
             param( 8) =  -2.345830D-01
             param( 9) =   2.368173D-01
             param(10) =  -9.913890D-04
             param(11) =  -1.146165D-02
             param(12) =   0.000000D+00
c
             param(13) =   1.674634D+00
             param(14) =   5.732017D+01
             param(15) =   5.955416D+01
             param(16) =  -2.311007D+02
             param(17) =   1.255199D+02
             param(18) =   1.023254D-01
             param(19) =  -2.453783D+00
             param(20) =   2.913180D+01
             param(21) =  -3.494358D+01
             param(22) =   2.315955D+01
c
             call nwxc_c_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M06_L)
c
             param( 1) =   3.957626D-01
             param( 2) =  -5.614546D-01
             param( 3) =   1.403963D-02
             param( 4) =   9.831442D-04
             param( 5) =  -3.577176D-03
             param( 6) =   0.000000D+00
             param( 7) =   4.650534D-01
             param( 8) =   1.617589D-01
             param( 9) =   1.833657D-01
             param(10) =   4.692100D-04
             param(11) =  -4.990573D-03
             param(12) =   0.000000D+00
c
             param(13) =   6.042374D-01
             param(14) =   1.776783D+02
             param(15) =  -2.513252D+02
             param(16) =   7.635173D+01
             param(17) =  -1.255699D+01
             param(18) =   5.349466D-01
             param(19) =   5.396620D-01
             param(20) =  -3.161217D+01
             param(21) =   5.149592D+01
             param(22) =  -2.919613D+01
c
             call nwxc_c_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M06_2X)
c
             param( 1) =   1.166404D-01
             param( 2) =  -9.120847D-02
             param( 3) =  -6.726189D-02
             param( 4) =   6.720580D-05
             param( 5) =   8.448011D-04
             param( 6) =   0.000000D+00
             param( 7) =   6.902145D-01
             param( 8) =   9.847204D-02
             param( 9) =   2.214797D-01
             param(10) =  -1.968264D-03
             param(11) =  -6.775479D-03
             param(12) =   0.000000D+00
c
             param(13) =   8.833596D-01
             param(14) =   3.357972D+01
             param(15) =  -7.043548D+01
             param(16) =   4.978271D+01
             param(17) =  -1.852891D+01
             param(18) =   3.097855D-01
             param(19) =  -5.528642D+00
             param(20) =   1.347420D+01
             param(21) =  -3.213623D+01
             param(22) =   2.846742D+01
c
             call nwxc_c_m06_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M08_HX)
c            parameters A
             param(1)  =    1.0000000D+00
             param(2)  =   -4.0661387D-01
             param(3)  =   -3.3232530D+00
             param(4)  =    1.5540980D+00
             param(5)  =    4.4248033D+01
             param(6)  =   -8.4351930D+01
             param(7)  =   -1.1955581D+02
             param(8)  =    3.9147081D+02
             param(9)  =    1.8363851D+02
             param(10) =   -6.3268223D+02
             param(11) =   -1.1297403D+02
             param(12) =    3.3629312D+02
c            parameters B
             param(13) =    1.3812334D+00
             param(14) =   -2.4683806D+00
             param(15) =   -1.1901501D+01
             param(16) =   -5.4112667D+01
             param(17) =    1.0055846D+01
             param(18) =    1.4800687D+02
             param(19) =    1.1561420D+02
             param(20) =    2.5591815D+02
             param(21) =    2.1320772D+02
             param(22) =   -4.8412067D+02
             param(23) =   -4.3430813D+02
             param(24) =    5.6627964D+01
             call nwxc_c_m11_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M08_SO)
c            parameters A
             param(1)  =   1.0000000D+00
             param(2)  =   0.0000000D+00
             param(3)  =  -3.9980886D+00
             param(4)  =   1.2982340D+01
             param(5)  =   1.0117507D+02
             param(6)  =  -8.9541984D+01
             param(7)  =  -3.5640242D+02
             param(8)  =   2.0698803D+02
             param(9)  =   4.6037780D+02
             param(10) =  -2.4510559D+02
             param(11) = -1.9638425D+02
             param(12) =  1.1881459D+02
c            parameters B
             param(13) =   1.0000000D+00
             param(14) =  -4.4117403D+00
             param(15) =  -6.4128622D+00
             param(16) =   4.7583635D+01
             param(17) =   1.8630053D+02
             param(18) =  -1.2800784D+02
             param(19) =  -5.5385258D+02
             param(20) =   1.3873727D+02
             param(21) =   4.1646537D+02
             param(22) =  -2.6626577D+02
             param(23) =   5.6676300D+01
             param(24) =   3.1673746D+02
             call nwxc_c_m11_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M11)
c            parameters A
             param(1)  =  1.0000000D+00
             param(2)  =  0.0000000D+00
             param(3)  = -3.8933250D+00
             param(4)  = -2.1688455D+00
             param(5)  =  9.3497200D+00
             param(6)  = -1.9845140D+01
             param(7)  =  2.3455253D+00
             param(8)  =  7.9246513D+01
             param(9)  =  9.6042757D+00
             param(10) = -6.7856719D+01
             param(11) = -9.1841067D+00
             param(12) =  0.0000000D+00
c            parameters B
             param(13) =  7.2239798D-01
             param(14) =  4.3730564D-01
             param(15) = -1.6088809D+01
             param(16) = -6.5542437D+01
             param(17) =  3.2057230D+01
             param(18) =  1.8617888D+02
             param(19) =  2.0483468D+01
             param(20) = -7.0853739D+01
             param(21) =  4.4483915D+01
             param(22) = -9.4484747D+01
             param(23) = -1.1459868D+02
             param(24) =  0.0000000D+00
             call nwxc_c_m11_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_M11_L)
c            parameters A
             param(1)  =  1.000000D+00
             param(2)  =  0.000000D+00
             param(3)  =  2.750880D+00
             param(4)  = -1.562287D+01
             param(5)  =  9.363381D+00
             param(6)  =  2.141024D+01
             param(7)  = -1.424975D+01
             param(8)  = -1.134712D+01
             param(9)  =  1.022365D+01
             param(10) =  0.000000D+00
             param(11) =  0.000000D+00
             param(12) =  0.000000D+00
c            parameters B
             param(13) =  1.000000D+00
             param(14) = -9.082060D+00
             param(15) =  6.134682D+00
             param(16) = -1.333216D+01
             param(17) = -1.464115D+01
             param(18) =  1.713143D+01
             param(19) =  2.480738D+00
             param(20) = -1.007036D+01
             param(21) = -1.117521D-01
             param(22) =  0.000000D+00
             param(23) =  0.000000D+00
             param(24) =  0.000000D+00
             call nwxc_c_m11_d3(param,nwxc_rho_tol,ipol,nq,
     +            nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_OP)
             param(1) = 2.3670d0
             call nwxc_c_op_d3(nwxc_k_becke88_d3,param,nwxc_rho_tol,
     +            ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_OPT)
            call nwxc_c_opt_d3(nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,f)
           case (NWXCP_C_PW6B95)
             param(1) = 0.00262d0
             param(2) = 0.03668d0
             call nwxc_c_b95_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,tau,f)
           case (NWXCP_C_PWB6K)
             param(1) = 0.00353d0
             param(2) = 0.04120d0
             call nwxc_c_b95_d3(param,nwxc_rho_tol,ipol,nq,
     +                          nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_PW91LDA)
            call nwxc_c_pw91lda_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                             rho,f)
          case (NWXCP_C_LYP)
            call nwxc_c_lyp_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,f)
          case (NWXCP_C_PZ81)
            call nwxc_c_perdew81_d3(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),rho,f)
          case (NWXCP_C_P86)
            call nwxc_c_perdew86_d3(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_P91)
            call nwxc_c_p91_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,gamma,f)
          case (NWXCP_C_P91_VWN5)
            call nwxc_c_p91_vwn5_d3(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PBE)
            param(1) = 0.066724550603149d0
            call nwxc_c_pbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PBESOL)
            param(1) = 0.046d0
            call nwxc_c_pbe_d3(param,nwxc_rho_tol,ipol,nq,
     +                         nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_PKZB)
            call nwxc_c_pkzb99_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                            rho,gamma,tau,f)
          case (NWXCP_C_SPBE)
            call nwxc_c_spbe_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_C_SOGGA11)
            param(1)  =  5.00000d-01
            param(2)  = -4.62334D+00
            param(3)  =  8.00410D+00
            param(4)  = -130.226D+00
            param(5)  =  38.2685D+00
            param(6)  =  69.5599D+00
            param(7)  =  5.00000d-01
            param(8)  =  3.62334D+00
            param(9)  =  9.36393D+00
            param(10) =  34.5114D+00
            param(11) = -18.5684D+00
            param(12) = -0.16519D+00
            call nwxc_c_sogga_d3(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_SOGGA11_X)
            param(1)  =  5.00000d-01
            param(2)  =  7.82439d+01
            param(3)  =  2.57211d+01
            param(4)  = -1.38830d+01
            param(5)  = -9.87375d+00
            param(6)  = -1.41357d+01
            param(7)  =  5.00000d-01
            param(8)  = -7.92439d+01
            param(9)  =  1.63725d+01
            param(10) =  2.08129d+00
            param(11) =  7.50769d+00
            param(12) = -1.01861d+01
            call nwxc_c_sogga_d3(param,nwxc_rho_tol,ipol,nq,
     +                           nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_C_TPSS)
            call nwxc_c_tpss03_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                            rho,gamma,tau,f)
          case (NWXCP_C_VS98)
            param(1)  =  7.035010d-01
            param(2)  =  7.694574d-03
            param(3)  =  5.152765d-02
            param(4)  =  3.394308d-05
            param(5)  = -1.269420d-03
            param(6)  =  1.296118d-03
            param(7)  =  3.270912d-01
            param(8)  = -3.228915d-02
            param(9)  = -2.942406d-02
            param(10) =  2.134222d-03
            param(11) = -5.451559d-03
            param(12) =  1.577575d-02
            call nwxc_c_vs98_d3(param,nwxc_rho_tol,ipol,nq,
     +           nwxc_wghts(it),rho,gamma,tau,f)
          case (NWXCP_C_VWN1)
            call nwxc_c_vwn1_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN1_RPA)
            call nwxc_c_vwn1_rpa_d3(nwxc_rho_tol,ipol,nq,
     +                              nwxc_wghts(it),
     +                              rho,f)
          case (NWXCP_C_VWN2)
            call nwxc_c_vwn2_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN3)
            call nwxc_c_vwn3_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN4)
            call nwxc_c_vwn4_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_C_VWN5)
            call nwxc_c_vwn5_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,f)
          case (NWXCP_BOP)
            param(1) = 2.3670d0
            call nwxc_c_op_d3(nwxc_k_becke88_d3,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_KT1)
            call nwxc_xc_kt1_d3(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                          rho,gamma,f)
          case (NWXCP_PBEOP)
            param(1) = 2.3789d0
            call nwxc_c_op_d3(nwxc_k_pbe96_d3,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case (NWXCP_SOP)
            param(1) = 2.5654d0
            call nwxc_c_op_d3(nwxc_k_dirac_d3,param,nwxc_rho_tol,
     +           ipol,nq,nwxc_wghts(it),rho,gamma,f)
          case default
            if (nwxc_oroot) then
              write(*,*)"nwxc_eval_df3: invalid functional",
     +                  nwxc_ids(it)
            endif
            call nwxc_printP()
            call errquit("nwxc_eval_df3: invalid functional",
     +                   nwxc_ids(it),0)
        end select
      enddo
C
      end
C> @}
c $Id$
