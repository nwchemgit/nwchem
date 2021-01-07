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
#include "stdio.fh"
#include "util.fh"      
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
      if (.not.ga_destroy(pars%idx_oo))
     $  call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
      if (.not.ga_destroy(pars%idx_ov))
     $  call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
      if (pars%nstates.gt.pars%nocc(1)) then
        if (.not.ga_destroy(pars%g_erivv))
     $    call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
        if (.not.ga_destroy(pars%idx_vv))
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

      pars%timer%total = util_wallsec() - pars%timer%total

      if (pars%me.eq.0) then
        write(luout,*)
        write(luout,*) "    GW Timing Statistics"
        write(luout,*) "  ------------------------"
        write(luout,'("   Init.   : ",F12.1,"s")') pars%timer%init
        write(luout,'("   ERIs    : ",F12.1,"s")') pars%timer%eris
        write(luout,'("   Distr.  : ",F12.1,"s")') pars%timer%distr
        write(luout,'("   Vxc     : ",F12.1,"s")') pars%timer%vxc
        if (pars%cdgw) then
          write(luout,'("   Vmo     : ",F12.1,"s")') pars%timer%vmo
          write(luout,'("   iW      : ",F12.1,"s")') pars%timer%iw
          write(luout,'("   I_n     : ",F12.1,"s")') pars%timer%in
          write(luout,'("   R_n     : ",F12.1,"s")') pars%timer%rn
        else
          write(luout,'("   Sigma_x : ",F12.1,"s")') pars%timer%sigmax
          write(luout,'("   RPA     : ",F12.1,"s")') pars%timer%rpa
          write(luout,'("   Wmn     : ",F12.1,"s")') pars%timer%wm
          write(luout,'("   Sigma_c : ",F12.1,"s")') pars%timer%sigmac
        endif
        write(luout,  '("   GW total: ",F12.1,"s")') pars%timer%total
        write(luout,*)
      endif

      return
      end
