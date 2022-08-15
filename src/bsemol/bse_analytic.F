      subroutine bse_analytic(pars)

      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bse.fh"

      type(bse_params_t) :: pars
      character(*), parameter :: pname = 'bse_analytic: '

      integer lSize,isp
      integer l_omega, k_omega
      integer l_dia,k_dia
      integer lW,kW
      integer maxpoles,maxnpoles,totnpoles,nri

      ! synchronize parallel processes
      call ga_sync

      !initialize useful variables
      nri = pars%nri
      maxpoles = maxval(pars%mynpoles(1:pars%ipol))
      maxnpoles = maxval(pars%npoles(1:pars%ipol))
      totnpoles = maxnpoles*pars%ipol

      ! allocation
      if(.not.ma_push_get(mt_dbl,totnpoles,'omega',l_omega,k_omega))
     &  call errquit(pname//'failed to allocate omega',0,MA_ERR) 
      if(.not.ma_push_get(mt_dbl,totnpoles,'dia',l_dia,k_dia))
     &  call errquit(pname//'failed to allocate Delta_ia',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nri**2,'W',lW,kW))
     &  call errquit(pname//'failed to allocate W',0,MA_ERR)

      ! obtain "old" eigenvalue difference to build W
      do isp=1,pars%ipol
         lSize = (isp-1)*pars%npoles(1)
         call gw_get_eia(dbl_mb(pars%k_mf_evals+(isp-1)*pars%nmo),
     &                   dbl_mb(k_dia+lSize),pars%nocc(isp),
     &                   pars%nvir(isp))    
      enddo 

      ! obtain screened Coulomb matrix
      call bse_buildw(pars,dbl_mb(k_dia),dbl_mb(kW),pars%nmo,nri,
     &                maxpoles,pars%ipol)

      ! transform ERIs using screened Coulomb matrix
      call bse_wmn(pars,dbl_mb(kW),pars%nmo,nri,pars%ipol)
      if(.not.ma_chop_stack(lW))
     &  call errquit(pname//'failed to chop stack',0,MA_ERR) 

      ! obtain "GW" eigenvalue difference to build Casida matrix
      do isp=1,pars%ipol
         lSize = (isp-1)*pars%npoles(1)
         call bse_get_eia(dbl_mb(pars%k_gw_evals+(isp-1)*pars%nmo),
     &                   dbl_mb(k_dia+lSize),pars%nocc(isp),
     &                   pars%nvir(isp),pars%ovlo(isp),
     %                   pars%ovhi(isp))    
      enddo 

      ! Build and diagonalize Casida-like matrix
      if (pars%tda) then
        call bse_analytic_tda(pars,dbl_mb(k_dia),dbl_mb(k_omega),
     &                        maxpoles,totnpoles,pars%singlet)
      else
        call bse_analytic_cis(pars,dbl_mb(k_dia),dbl_mb(k_omega),
     &                        maxpoles,totnpoles,pars%singlet)
      endif

      ! Compute oscillator strengths
      call bse_oscstr(pars,dbl_mb(k_omega),totnpoles)

      ! Destroy GAs
      if (.not.ga_destroy(pars%g_apb))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR)

      ! chop stack
      if(.not.ma_chop_stack(l_omega))
     $  call errquit(pname//'failed to chop stack',0,MA_ERR) 

      end subroutine bse_analytic
