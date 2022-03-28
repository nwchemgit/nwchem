      subroutine gw_finalize(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
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
        if (pars%npoles(i).gt.0) then
        if (.not.ga_destroy(pars%g_eriov(i)))
     $    call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
        endif
        if (pars%nstates(i).eq.0) cycle
        if (.not.ga_destroy(pars%g_erioo(i)))
     $    call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
        if (pars%nstates(i).gt.pars%nocc(i)) then
          if (.not.ga_destroy(pars%g_erivv(i)))
     $      call errquit(pname//'failed to destroy g_moeris',0,GA_ERR) 
        endif
      enddo

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
  900   format("  Exact exchange energy   = ",F25.12)
 1000   format(" RPA correlation energy   = ",F25.12)
 1010   format("  GM correlation energy   = ",F25.12)
 1020   format("       Total RPA energy   = ",F25.12)
 1030   format("        Total GM energy   = ",F25.12)
        write(luout,*)
        if (pars%rpa) then
          write(luout,900) pars%ex
          write(luout,1000) pars%ecrpa
          write(luout,1010) pars%ecgm
          write(luout,*)
          write(luout,1020) pars%escf-pars%exscf-pars%ecscf+
     $                      pars%ex+pars%ecrpa
          write(luout,1030) pars%escf-pars%exscf-pars%ecscf+
     $                      pars%ex+pars%ecgm
          write(luout,*)
        endif
        write(luout,*) "    GW Timing Statistics"
        write(luout,*) "  ------------------------"
        write(luout,'("   Init.   : ",F12.1,"s")') pars%timer%init
        write(luout,'("     ERIs  : ",F12.1,"s")') pars%timer%eris
        write(luout,'("     Distr.: ",F12.1,"s")') pars%timer%distr
        write(luout,'("     Vxc   : ",F12.1,"s")') pars%timer%vxc
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
