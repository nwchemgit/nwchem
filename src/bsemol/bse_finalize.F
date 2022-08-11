      subroutine bse_finalize(pars)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "cdft.fh"
#include "bse.fh"

       type(bse_params_t) :: pars
       character(*),parameter :: pname = 'bse_finalize: '

       integer i
       logical ok

       do i=1,ipol
        if (.not.ga_destroy(pars%g_movecs(i)))
     $    call errquit(pname//'failed to destroy g_movecs',0,GA_ERR)
c        if (pars%npoles(i).gt.0) then
c        if (.not.ga_destroy(pars%g_eriov(i)))
c     $    call errquit(pname//'failed to destroy g_moeris',0,GA_ERR)
c        endif
       enddo

       ok = ma_chop_stack(pars%l_mf_evals)
       if (.not.ok)
     &   call errquit(pname//'could not chop stack',0,MA_ERR) 

       call schwarz_tidy()
       call int_terminate()
       call fock_2e_tidy(pars%rtdb)

      if (.not.geom_destroy(geom))
     $  call errquit(pname//'geom_destroy failed',0,GEOM_ERR)
      if (.not.bas_destroy(ao_bas_han))
     $  call errquit(pname//'bas_destroy failed',0,BASIS_ERR)
      if (.not.bas_destroy(cd_bas_han))
     $  call errquit(pname//'bas_destroy failed',0,BASIS_ERR)


       end subroutine
