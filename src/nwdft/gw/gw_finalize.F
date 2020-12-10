      subroutine gw_finalize(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "gw.fh"

      type(gw_params_t)  :: pars
c
      character(len=13), parameter :: pname = 'gw_finalize: '
c
      integer i
c
c     -----------------------------------------------------------------
c
      do i=1,ipol
        if (.not.ga_destroy(pars%g_movecs(i)))
     $    call errquit(pname//'failed to destroy g_movecs',0,GA_ERR)
      enddo
      if (.not.ga_destroy(pars%g_erioo))
     $  call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
      if (.not.ga_destroy(pars%g_eriov))
     $  call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
      if (pars%nstates.gt.pars%nocc(1)) then
        if (.not.ga_destroy(pars%g_erivv))
     $    call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
      endif
      if (.not.ma_chop_stack(pars%l_evals))
     $  call errquit(pname//'failed to chop stack',0,MA_ERR)
c
      call schwarz_tidy()
      call int_terminate()
      call fock_2e_tidy(pars%rtdb)

      if (.not.geom_destroy(geom))
     $  call errquit(pname//'geom_destroy failed',0,GEOM_ERR)
      if (.not.bas_destroy(ao_bas_han))
     $  call errquit(pname//'bas_destroy failed',0,BASIS_ERR)
      if (.not.bas_destroy(cd_bas_han))
     $  call errquit(pname//'bas_destroy failed',0,BASIS_ERR)
      return
      end
