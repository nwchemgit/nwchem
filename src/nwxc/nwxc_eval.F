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
       subroutine nwxc_eval_df(ipol,nq,rho,gamma,tau,f,dfdr,dfdg,dfdt)
       implicit none
#include "nwxcP.fh"
       integer ipol !< [Input] The number of spin channels
       integer nq   !< [Input] The number of points
C
       double precision rho(nq,ipol)     !< [Input] Density
       double precision gamma(nq,ipol+1) !< [Input] |Density gradient|
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
       integer nd1r(2) !< The number of partial derivatives wrt rho as
                       !< a function of ipol
       integer nd1g(2) !< The number of partial derivatives wrt gamma as
                       !< a function of ipol
       integer nd1t(2) !< The number of partial derivatives wrt tau as
                       !< a function of ipol
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
           case (NWXC_X_SLATER)
             call nwxc_x_dirac(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                         rho,f,dfdr)
           case (NWXC_X_PBE)
             call nwxc_x_pbe(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f,dfdr,dfdg)
           case (NWXC_C_PW91LDA)
             call nwxc_c_pw91lda(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                           rho,f,dfdr)
           case (NWXC_C_PBE)
             call nwxc_c_pbe(nwxc_rho_tol,ipol,nq,nwxc_wghts(it),
     +                       rho,gamma,f,dfdr,dfdg)
           case default
             call errquit("nwxc_eval_df: invalid functional",
     +                    nwxc_ids(it),UERR)
         end select
       enddo
C
       end
C> @}
