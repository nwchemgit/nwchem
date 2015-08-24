C> \ingroup nwxc
C> @{
C>
C> \file nwxc_eval_print.F
C> Routines that print the current density functional as a Maxima
C> expression
C>
C> The routines in this file evaluate the current density functional.
C> Using a modified version of the automatic differentiation library
C> they also write the expressions they evaluate in the language
C> of the Maxima symbolic algebra package. The Maxima expression is
C> written to a text file so that Maxima can do further processing
C> on it.
C>
C> @}
C
C> \ingroup nwxc_api
C> @{
C>
C> \brief Evaluate the exchange-correlation energy and the 1st order
C> partial derivatives and print the energy expression
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
#define NWAD_PRINT Y
#define NWAD_DEG 1
      subroutine nwxc_eval_df_print(ipol,nq,rho,gamma,tau,f,
     &                              dfdr,dfdg,dfdt)
      use nwadp1
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
      integer iq           !< Counter over elements
c
c     Work out how many double precision floating point words are needed
c     to represent a single nwad_dble
c
      ipl = 2
      memory_test_nwad(1)%d0   = 0.0d0
      memory_test_nwad(1)%ivar = 0
      memory_test_nwad(1)%nvar = 0
      memory_test_nwad(2)%d0   = 0.0d0
      memory_test_nwad(2)%ivar = 0
      memory_test_nwad(2)%nvar = 0
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
c     The Maxima printing capabilities depend on the initialization
c     of the NWAD_DBLE variables. Of course using the MA library
c     bypasses the regular Fortran instantiation of the variables.
c     Therefore, we must explicitly initialize the NWAD_DBLEs that
c     live on the MA memory segments.
c
      do iq = 0, nq-1
        call dcopy(length,memory_test_nwad(1),1,
     +                    dbl_mb(k_fun+length*iq),1)
      enddo
      do iq = 0, nq*ipl-1
        call dcopy(length,memory_test_nwad(1),1,
     +                    dbl_mb(k_rho+length*iq),1)
        call dcopy(length,memory_test_nwad(1),1,
     +                    dbl_mb(k_tau+length*iq),1)
      enddo
      do iq = 0, nq*(ipl+1)-1
        call dcopy(length,memory_test_nwad(1),1,
     +                    dbl_mb(k_gam+length*iq),1)
      enddo
c
      call nwxc_eval_df_driver_print(ipol,nq,rho,gamma,tau,
     +     dbl_mb(k_rho),dbl_mb(k_gam),dbl_mb(k_tau),dbl_mb(k_fun),f,
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
      subroutine nwxc_eval_df_driver_print(ipol,nq,rho,gamma,tau,
     +           nwad_rho,nwad_gam,nwad_tau,nwad_f,f,
     +           dfdr,dfdg,dfdt)
      use nwadp1
      implicit none
#include "nwxc_param.fh"
#include "errquit.fh"
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
      logical  nobeta       !< Is the beta spin channel not present (zero)?
      logical  nwxc_is_lda  !< Is the functional a GGA
      logical  nwxc_is_gga  !< Is the functional a GGA
      logical  nwxc_is_mgga !< Is the functional a meta-GGA
      external nwxc_is_gga  !< Is the functional a GGA
      external nwxc_is_mgga !< Is the functional a meta-GGA
c
      integer ix1           !< the variable rank
      integer iv1           !< the variable label
      integer nvar          !< the number of variables
      double precision val  !< the functional (derivative) value
      integer iq, ip
c
      if (ipol.eq.1) then
        do iq = 1, nq
          nwad_rho(iq,R_T) = set_rho_t(rho(iq,R_T))
        enddo
      else  ! ipol.eq.1
        do iq = 1, nq
          nwad_rho(iq,R_A) = set_rho_a(rho(iq,R_A))
        enddo
        do iq = 1, nq
          if (rho(iq,R_B).gt.0.0d0) then
            nwad_rho(iq,R_B) = set_rho_b(rho(iq,R_B))
          else
            nwad_rho(iq,R_B) = 0.0d0
          endif
        enddo
      endif ! ipol.eq.1
      if (nwxc_is_gga().or.nwxc_is_mgga()) then
        if (ipol.eq.1) then
          do iq = 1, nq
            nwad_gam(iq,G_TT) = set_gamma_tt(gamma(iq,G_TT))
          enddo
        else  ! ipol.eq.1
          do iq = 1, nq
            nwad_gam(iq,G_AA) = set_gamma_aa(gamma(iq,G_AA))
          enddo
          do iq = 1, nq
            if (gamma(iq,G_AB).gt.0.0d0) then
              nwad_gam(iq,G_AB) = set_gamma_ab(gamma(iq,G_AB))
            else
              nwad_gam(iq,G_AB) = 0.0d0
            endif
          enddo
          do iq = 1, nq
            if (gamma(iq,G_AB).gt.0.0d0) then
               nwad_gam(iq,G_BB) = set_gamma_bb(gamma(iq,G_BB))
            else
              nwad_gam(iq,G_BB) = 0.0d0
            endif
          enddo
        endif ! ipol.eq.1
      endif
      if (nwxc_is_mgga()) then
        if (ipol.eq.1) then
          do iq = 1, nq
            nwad_tau(iq,T_T) = set_tau_t(tau(iq,T_T))
          enddo
        else  ! ipol.eq.1
          do iq = 1, nq
            nwad_tau(iq,T_A) = set_tau_a(tau(iq,T_A))
          enddo
          do iq = 1, nq
            nwad_tau(iq,T_B) = set_tau_b(tau(iq,T_B))
          enddo
        endif ! ipol.eq.1
      endif
c
      call nwxc_eval_df_doit_print(ipol,nq,nwad_rho,nwad_gam,
     +     nwad_tau,nwad_f)
c
 10   format("f(zeta,rhoa,rhob,rho,sigmaaa,sigmaab,sigmabb,sigma,taua,ta
     +ub,tau):=",$)
      if (value(nwad_rho(1,2)).gt.0.0d0) then
        write(nwad_iunit,10)
        nobeta = .false.
      else
        write(nwad_iunit,'("G(rhoa,sigmaaa,taua,taub):=",$)')
        nobeta = .true.
      endif
      call nwad_print_entity(nwad_f(1))
      write(nwad_iunit,'(";")')
c
c     Now unpack the results
c
      do iq = 1, nq
        f(iq) = get_val(nwad_f(iq))
        nvar = get_nvar(nwad_f(iq))
        do ix1 = 1, nvar
          call get_d1(nwad_f(iq),ix1,val,iv1)
          select case (iv1)
            case (1)
              dfdr(iq,D1_RA) = val
            case (2)
              dfdr(iq,D1_RB) = val
              if (nobeta) then
                write(nwad_iunit,'("illegal argument rhob")')
                call errquit("nwxc_df_driver_print: illegal variable",
     +                       iv1,UERR)
              endif
            case (3)
              dfdg(iq,D1_GAA) = val
            case (4)
              dfdg(iq,D1_GAB) = val
              if (nobeta) then
                write(nwad_iunit,'("illegal argument gamma_ab")')
                call errquit("nwxc_df_driver_print: illegal variable",
     +                       iv1,UERR)
              endif
            case (5)
              dfdg(iq,D1_GBB) = val
              if (nobeta) then
                write(nwad_iunit,'("illegal argument gamma_bb")')
                call errquit("nwxc_df_driver_print: illegal variable",
     +                       iv1,UERR)
              endif
            case (6)
              dfdt(iq,D1_TA) = val
            case (7)
              dfdt(iq,D1_TB) = val
            case default
              call print(6,nwad_f(iq))
              call errquit("nwxc_df_driver_print: illegal variable",
     +                     iv1,UERR)
          end select
        enddo
      enddo
c
      end
C>
C> \brief Evaluate the functional and its first order derivatives
C>
C> Using automatic differentiation the functional and its first order
C> derivatives are calculated for one direction. The loop over directions
C> and the correct interpolation of the partial derivatives has to be handled
C> in a driver routine.
C>
      subroutine nwxc_eval_df_doit_print(ipol,nq,rho,gamma,tau,f)
      use nwadp1
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
!     double precision dfdr(nq,ipol)   !< [Output] Derivative wrt rho
!     double precision dfdg(nq,ipol+1) !< [Output] Derivative wrt gamma
!     double precision dfdt(nq,ipol)   !< [Output] Derivative wrt tau
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
      external nwxc_k_dirac_p
      external nwxc_k_becke88_p
      external nwxc_k_pbe96_p
c
      integer max_param
      parameter (max_param = 50)
      type(nwad_dble)::param(max_param)
      double precision weight
      character*(9) name
      double precision Axlsda
      parameter (Axlsda = -0.9305257363491d0 )
C
      data nd1r / D1_RA,  D1_RB  /
      data nd1g / D1_GAA, D1_GBB /
      data nd1t / D1_TA,  D1_TB  /
C
      do iq = 1, max_param
        write(name,'("param(",i2.2,")")')iq
        call set_name(param(iq),name)
      enddo
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
      next_term: do it = 1, nwxc_num_f
c
c       If the next term has a tiny weight skip it!
c
        if (abs(nwxc_wghts(it)).lt.1.0d-10) cycle next_term
c
        select case (nwxc_ids(it))
          case (NWXCP_X_SLATER)
            weight = nwxc_wghts(it)
            call nwxc_x_dirac_p(nwxc_rho_tol,ipol,nq,weight,
     +                          rho,f)
          case (NWXCP_X_B86B)
            weight = nwxc_wghts(it)
            call nwxc_x_b86b_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_B88)
            weight = nwxc_wghts(it)
            call nwxc_x_b88_p(nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B97)
            param(1) = 2.0d0
            param(2) = 0.80940d+00
            param(3) = 0.50730d+00
            param(4) = 0.74810d+00
            weight   = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B97_1)
            param(1) = 2.0d0
            param(2) = 0.789518d+00
            param(3) = 0.573805d+00
            param(4) = 0.660975d+00
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B97_2)
            param(1) = 2.0d0
            param(2) = 0.827642D+00
            param(3) = 0.478400D-01
            param(4) = 0.176125D+01
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B97_3)
            param(1) = 4.0d0
            param(2) = 7.334648D-01
            param(3) = 2.925270D-01
            param(4) = 3.338789D+00
            param(5) =-1.051158D+01
            param(6) = 1.060907D+01
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B97_D)
            param(1) = 2.0d0
            param(2) = 1.086620d+0
            param(3) =-0.521270d+00
            param(4) = 3.254290d+00
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B97_G)
            param(1) = 2.0d0
            param(2) = 1.1068d0
            param(3) =-0.8765d0
            param(4) = 4.2639d0
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_B98)
            param(1) = 2.0d0
            param(2) = 0.790194d00
            param(3) = 0.400271d00
            param(4) = 0.832857d00
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_BNL)
            param(1) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_bnl_p(param,nwxc_rho_tol,ipol,nq,
     +                        weight,rho,f)
          case (NWXCP_X_CAMB88)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_camb88_p(param,nwxc_rho_tol,ipol,nq,
     +                           weight,rho,gamma,f)
          case (NWXCP_X_CAMLSD)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_camlsd_p(param,nwxc_rho_tol,ipol,nq,
     +                           weight,rho,f)
          case (NWXCP_X_CAMPBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_campbe_p(param,nwxc_rho_tol,ipol,nq,
     +                           weight,rho,gamma,f)
          case (NWXCP_X_CAMREVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            param(3) = nwxc_cam_alpha
            param(4) = nwxc_cam_beta
            param(5) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_campbe_p(param,nwxc_rho_tol,ipol,nq,
     +                           weight,rho,gamma,f)
          case (NWXCP_X_CAMRPBE)
            param(1) = nwxc_cam_alpha
            param(2) = nwxc_cam_beta
            param(3) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_camrpbe_p(param,nwxc_rho_tol,ipol,nq,
     +                            weight,rho,gamma,f)
          case (NWXCP_X_CAMS12G)
            param(1) = 1.03323556d0
            param(2) = 0.757d0
            param(3) = 0.00417251d0
            param(4) = 0.00115216d0
            param(5) = 0.00706184d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_cams12_p(param,nwxc_rho_tol,ipol,nq,
     +                           weight,rho,gamma,f)
          case (NWXCP_X_CAMS12H)
            param(1) = 1.02149642d0
            param(2) = 0.757d0
            param(3) = 0.00825905d0
            param(4) = 0.00235804d0
            param(5) = 0.00654977d0
            param(6) = nwxc_cam_alpha
            param(7) = nwxc_cam_beta
            param(8) = nwxc_cam_gamma
            weight = nwxc_wghts(it)
            call nwxc_x_cams12_p(param,nwxc_rho_tol,ipol,nq,
     +                           weight,rho,gamma,f)
          case (NWXCP_X_DLDF)
            weight = nwxc_wghts(it)
            call nwxc_x_dldf_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,tau,f)
          case (NWXCP_X_FT97)
            weight = nwxc_wghts(it)
            call nwxc_x_ft97_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_GILL)
            weight = nwxc_wghts(it)
            call nwxc_x_gill_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_HCTH)
            param(1) = 4.0d0
            param(2) = 0.109320d+01
            param(3) =-0.744056d+00
            param(4) = 0.559920d+01
            param(5) =-0.678549d+01
            param(6) = 0.449357d+01
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_HCTH120)
            param(1) = 4.0d0
            param(2) = 1.09163d0
            param(3) =-0.74720d0
            param(4) = 5.07830d0
            param(5) =-4.10750d0
            param(6) = 1.17170d0
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_HCTH147)
            param(1) = 4.0d0
            param(2) = 1.09025d0
            param(3) =-0.79920d0
            param(4) = 5.57210d0
            param(5) =-5.86760d0
            param(6) = 3.04540d0
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_HCTH407)
            param(1) = 4.0d0
            param(2) = 1.08184d0
            param(3) =-0.5183d0
            param(4) = 3.4256d0
            param(5) =-2.6290d0
            param(6) = 2.2886d0
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_HCTH407P)
            param(1) = 4.0d0
            param(2) = 1.08018D0
            param(3) =-0.4117D0
            param(4) = 2.4368D0
            param(5) = 1.3890D0
            param(6) =-1.3529D0
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_HCTH_A)
            param(1) = 2.0d0
            param(2) = 0.109878d+01
            param(3) =-0.251173d+01
            param(4) = 0.156233d-01
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_HCTHP14)
            param(1) = 4.0d0
            param(2) = 0.103161d+01
            param(3) =-0.360781d+00
            param(4) = 0.351994d+01
            param(5) =-0.495944d+01
            param(6) = 0.241165d+01
            weight = nwxc_wghts(it)
            call nwxc_x_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight = nwxc_wghts(it)
            call nwxc_x_m05_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight = nwxc_wghts(it)
            call nwxc_x_m05_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight = nwxc_wghts(it)
            call nwxc_x_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m06_2x_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m08_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m08_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m11_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_m11_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
          case (NWXCP_X_MPW91)
            param(1) = 3.72d0
            param(2) = 0.00426D0
            weight    = nwxc_wghts(it)
            call nwxc_x_pw91_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_OPT)
            weight    = nwxc_wghts(it)
            call nwxc_x_opt_p(nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_PW6B95)
            param(1) = 0.00538d0
            param(2) = 1.7382d0
            param(3) = 3.8901d0
            weight    = nwxc_wghts(it)
            call nwxc_x_pw6_p(param,nwxc_rho_tol,ipol,nq,
     +                        weight,rho,gamma,f)
          case (NWXCP_X_PWB6K)
            param(1) = 0.00539d0
            param(2) = 1.7077d0
            param(3) = 4.0876d0
            weight    = nwxc_wghts(it)
            call nwxc_x_pw6_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_PW91)
            param(1) = 4.0d0
            param(2) = 0.0042D0
            weight    = nwxc_wghts(it)
            call nwxc_x_pw91_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_PBE)
            param(1) = 0.8040d0
            param(2) = 0.2195149727645171d0
            weight    = nwxc_wghts(it)
            call nwxc_x_pbe_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_PBESOL)
            param(1) = 0.8040d0
            param(2) = 10.0d0/81.0d0
            weight    = nwxc_wghts(it)
            call nwxc_x_pbe_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_REVPBE)
            param(1) = 1.245d0
            param(2) = 0.2195149727645171d0
            weight    = nwxc_wghts(it)
            call nwxc_x_pbe_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_X_RPBE)
            weight    = nwxc_wghts(it)
            call nwxc_x_rpbe_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_PKZB)
            weight    = nwxc_wghts(it)
            call nwxc_x_pkzb99_p(nwxc_rho_tol,ipol,nq,weight,
     +                           rho,gamma,tau,f)
          case (NWXCP_X_S12G)
            param(1) = 1.03842032d0
            param(2) = 0.757d0
            param(3) = 0.00403198d0
            param(4) = 0.00104596d0
            param(5) = 0.00594635d0
            weight    = nwxc_wghts(it)
            call nwxc_x_s12_p(param,nwxc_rho_tol,ipol,nq,
     +                        weight,rho,gamma,f)
          case (NWXCP_X_S12H)
            param(1) = 1.02543951d0
            param(2) = 0.757d0
            param(3) = 0.00761554d0
            param(4) = 0.00211063d0
            param(5) = 0.00604672d0
            weight    = nwxc_wghts(it)
            call nwxc_x_s12_p(param,nwxc_rho_tol,ipol,nq,
     +                        weight,rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_sogga_p(param,nwxc_rho_tol,ipol,nq,
     +                          weight,rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_sogga_p(param,nwxc_rho_tol,ipol,nq,
     +                          weight,rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_x_sogga_p(param,nwxc_rho_tol,ipol,nq,
     +                          weight,rho,gamma,f)
          case (NWXCP_X_SSB_D)
            weight    = nwxc_wghts(it)
            call nwxc_x_ssbD_1_p(nwxc_rho_tol,ipol,nq,
     +                           weight,rho,gamma,f)
          case (NWXCP_X_TPSS)
            weight    = nwxc_wghts(it)
            call nwxc_x_tpss03_p(nwxc_rho_tol,ipol,nq,weight,
     +                           rho,gamma,tau,f)
          case (NWXCP_X_WPBE)
            param(1) = nwxc_cam_gamma
            weight    = nwxc_wghts(it)
            call nwxc_x_wpbe_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_X_VS98)
            param(1) =  -9.800683d-01
            param(2) =  -3.556788d-03
            param(3) =   6.250326d-03
            param(4) =  -2.354518d-05
            param(5) =  -1.282732d-04
            param(6) =   3.574822d-04
            weight    = nwxc_wghts(it)
            call nwxc_x_vs98_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,tau,f)
          case (NWXCP_C_B95)
            param(1) = 0.0031d0
            param(2) = 0.038d0
            weight    = nwxc_wghts(it)
            call nwxc_c_b95_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
          case (NWXCP_C_B97)
            param(1) = 2.0d0
            param(2) = 0.17370d+00
            param(3) = 0.94540d+00
            param(4) = 0.23487d+01
            param(5) = 0.74710d+00
            param(6) =-0.24868d+01
            param(7) =-0.45961d+01
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_B97_1)
            param(1) = 2.0d0
            param(2) = 0.820011d-01
            param(3) = 0.955689d+00
            param(4) = 0.271681d+01
            param(5) = 0.788552d+00
            param(6) =-0.287103d+01
            param(7) =-0.547869d+01
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_B97_2)
            param(1) = 2.0d0
            param(2) = 0.585808D+00
            param(3) = 0.999849D+00
            param(4) =-0.691682D+00
            param(5) = 0.140626D+01
            param(6) = 0.394796D+00
            param(7) =-0.744060D+01
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_B97_D)
            param(1) = 2.0d0
            param(2) = 0.22340d+00
            param(3) = 0.690410d+00
            param(4) =-1.562080d+00
            param(5) = 6.302700d00
            param(6) = 1.942930d+0
            param(7) =-14.97120d+00
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_B97_G)
            param(1) = 2.0d0
            param(2) = 0.4883d0
            param(3) = 0.7961d0
            param(4) =-2.117d0
            param(5) = 5.7060d0
            param(6) = 2.3235d0
            param(7) =-14.9820d0
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_B98)
            param(1) = 2.0d0
            param(2) =-0.120163d00
            param(3) = 0.934715d00
            param(4) = 2.82332d0
            param(5) = 1.14105d0
            param(6) =-2.59412d0
            param(7) =-5.33398d0
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_DLDF)
            weight    = nwxc_wghts(it)
            call nwxc_c_dldf_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,tau,f)
          case (NWXCP_C_FT97)
            weight    = nwxc_wghts(it)
            call nwxc_c_ft97_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_b97_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m05_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m05_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m06_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m11_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m11_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m11_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_m11_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,tau,f)
          case (NWXCP_C_MPBE)
            param(1) = 0.066724550603149d0
            weight   = nwxc_wghts(it)
            call nwxc_c_mpbe_p(param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,f)
          case (NWXCP_C_OP)
            param(1) = 2.3670d0
            weight    = nwxc_wghts(it)
            call nwxc_c_op_p(nwxc_k_becke88_p,param,nwxc_rho_tol,ipol,
     +           nq,weight,rho,gamma,f)
          case (NWXCP_C_OPT)
            weight    = nwxc_wghts(it)
            call nwxc_c_opt_p(nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,f)
          case (NWXCP_C_PW6B95)
            param(1) = 0.00262d0
            param(2) = 0.03668d0
            weight    = nwxc_wghts(it)
            call nwxc_c_b95_p(param,nwxc_rho_tol,ipol,nq,
     +                        weight,rho,gamma,tau,f)
          case (NWXCP_C_PWB6K)
            param(1) = 0.00353d0
            param(2) = 0.04120d0
            weight    = nwxc_wghts(it)
            call nwxc_c_b95_p(param,nwxc_rho_tol,ipol,nq,
     +                        weight,rho,gamma,tau,f)
          case (NWXCP_C_PW91LDA)
            weight    = nwxc_wghts(it)
            call nwxc_c_pw91lda_p(nwxc_rho_tol,ipol,nq,weight,
     +                            rho,f)
          case (NWXCP_C_LYP)
            weight    = nwxc_wghts(it)
            call nwxc_c_lyp_p(nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_PZ81)
            weight    = nwxc_wghts(it)
            call nwxc_c_perdew81_p(nwxc_rho_tol,ipol,nq,weight,
     +                             rho,f)
          case (NWXCP_C_P86)
            weight    = nwxc_wghts(it)
            call nwxc_c_perdew86_p(nwxc_rho_tol,ipol,nq,weight,
     +                             rho,gamma,f)
          case (NWXCP_C_P91)
            weight    = nwxc_wghts(it)
            call nwxc_c_p91_p(nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_P91_VWN5)
            weight    = nwxc_wghts(it)
            call nwxc_c_p91_vwn5_p(nwxc_rho_tol,ipol,nq,weight,
     +                             rho,gamma,f)
          case (NWXCP_C_PBE)
            param(1) = 0.066724550603149d0
            weight    = nwxc_wghts(it)
            call nwxc_c_pbe_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_PBESOL)
            param(1) = 0.046d0
            weight    = nwxc_wghts(it)
            call nwxc_c_pbe_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                        rho,gamma,f)
          case (NWXCP_C_PKZB)
            param(1) = 0.066724550603149d0
            weight    = nwxc_wghts(it)
            call nwxc_c_pkzb99_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                           rho,gamma,tau,f)
          case (NWXCP_C_SPBE)
            weight    = nwxc_wghts(it)
            call nwxc_c_spbe_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_sogga_p(param,nwxc_rho_tol,ipol,nq,
     +                          weight,rho,gamma,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_sogga_p(param,nwxc_rho_tol,ipol,nq,
     +                          weight,rho,gamma,f)
          case (NWXCP_C_TPSS)
            weight    = nwxc_wghts(it)
            param(1) = 0.066724550603149d0
            call nwxc_c_tpss03_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                           rho,gamma,tau,f)
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
            weight    = nwxc_wghts(it)
            call nwxc_c_vs98_p(param,nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,tau,f)
          case (NWXCP_C_VWN1)
            weight    = nwxc_wghts(it)
            call nwxc_c_vwn1_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,f)
          case (NWXCP_C_VWN1_RPA)
            weight    = nwxc_wghts(it)
            call nwxc_c_vwn1_rpa_p(nwxc_rho_tol,ipol,nq,weight,
     +                             rho,f)
          case (NWXCP_C_VWN2)
            weight    = nwxc_wghts(it)
            call nwxc_c_vwn2_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,f)
          case (NWXCP_C_VWN3)
            weight    = nwxc_wghts(it)
            call nwxc_c_vwn3_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,f)
          case (NWXCP_C_VWN4)
            weight    = nwxc_wghts(it)
            call nwxc_c_vwn4_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,f)
          case (NWXCP_C_VWN5)
            weight    = nwxc_wghts(it)
            call nwxc_c_vwn5_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,f)
          case (NWXCP_BOP)
            param(1)  = 2.3670d0
            weight    = nwxc_wghts(it)
            call nwxc_c_op_p(nwxc_k_becke88_p,param,nwxc_rho_tol,ipol,
     +           nq,weight,rho,gamma,f)
          case (NWXCP_KT1)
            weight    = nwxc_wghts(it)
            call nwxc_xc_kt1_p(nwxc_rho_tol,ipol,nq,weight,
     +                         rho,gamma,f)
          case (NWXCP_PBEOP)
            param(1)  = 2.3789d0
            weight    = nwxc_wghts(it)
            call nwxc_c_op_p(nwxc_k_pbe96_p,param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,f)
          case (NWXCP_SOP)
            param(1)  = 2.5654d0
            weight    = nwxc_wghts(it)
            call nwxc_c_op_p(nwxc_k_dirac_p,param,nwxc_rho_tol,ipol,nq,
     +           weight,rho,gamma,f)
          case default
            if (nwxc_oroot) then
              write(*,*)"nwxc_eval_df: invalid functional",nwxc_ids(it)
            endif
            call nwxc_printP()
            call errquit("nwxc_eval_df: invalid functional",
     +                    nwxc_ids(it),0)
        end select
      enddo next_term
C
      end
C>
C> @}
