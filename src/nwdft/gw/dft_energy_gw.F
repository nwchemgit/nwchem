      logical function dft_energy_gw(rtdb)
      implicit none
#include "rtdb.fh"
      integer rtdb
c
      logical nwdft
      external nwdft
c
      logical gw_dft
      external gw dft
c
      logical status
c
      status = rtdb_cput(rtdb,'dft:theory', 1, 'dft')
      status = nwdft(rtdb)
      if (status) status = gw_dft(rtdb)
      dft_energy_gw = status
      call grid_cleanup(.false.)
c
      end function
c
c     -----------------------------------------------------------------
c
      logical function gw_dft(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "case.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "gw.fh"

      integer, intent(in) :: rtdb
      character(*), parameter :: pname = "gw_dft: "

      type(gw_params_t) :: params
      character*80 theory
      integer me
c
c     -----------------------------------------------------------------
c
      gw_dft = .false.
c
c     Determine the theory we are using
c
      if (.not. rtdb_cget(rtdb, "dft:theory", 1, theory))
     $  call errquit(pname//"theory not specified",0,RTDB_ERR)
c
c     Initialize GW calculation
c
      params%rtdb = rtdb
      params%me = ga_nodeid()
c
      if (params%me.eq.0) then
        write(luout,*) " "
        write(luout,'(22X,A)') "***********************"
        write(luout,'(22X,A)') "***      GW@DFT     ***"
        write(luout,'(22X,A)') "***********************"
        write(luout,*) " "
        call util_flush(luout)
      endif
      call gw_init(params)
c
      if (theory.eq."dft") then
        
        if (params%cdgw) then
          call gw_cdgw(params)
        else
          if (ipol.eq.1) then
            call gw_run(params)
          else
            call errquit(pname//'Analytic open-shell not implemented',
     $                   0,0)
          endif
        endif
c
      elseif (theory.eq.'sodft') then
        call errquit(pname//"spin-orbit not implemented",0,0)
      else
        call errquit(pname//"invalid theory: "//trim(theory),0,0)
      endif
c
      call gw_finalize(params)
      gw_dft = .true.
c
c     -----------------------------------------------------------------
c
      end function
