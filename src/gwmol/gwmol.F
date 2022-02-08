      logical function gwmol(rtdb)
      implicit none
#include "rtdb.fh"
      integer rtdb
c
      logical nwdft
      external nwdft
c
      logical gw_doit
      external gw doit
c
      logical status
c
      status = rtdb_cput(rtdb,'dft:theory', 1, 'dft')
      if (status) status = gw_doit(rtdb)
      gwmol = status
      call grid_cleanup(.false.)
c
      end function
c
c     -----------------------------------------------------------------
c
      logical function gw_doit(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "gw.fh"

      integer, intent(in) :: rtdb
      character(*), parameter :: pname = "gw_doit: "

      type(gw_params_t) :: params
      character*80 theory
      integer me
c
c     -----------------------------------------------------------------
c
c     Initialize
c
      gw_doit = .false.
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
      call gw_init(params)
c
      if (theory.eq."dft") then
        
        if (params%cdgw) then
          call gw_cdgw(params)
        else
          call gw_analytic(params)
        endif
c
      elseif (theory.eq.'sodft') then
        call errquit(pname//"spin-orbit not implemented",0,0)
      else
        call errquit(pname//"invalid theory: "//trim(theory),0,0)
      endif
c
      call gw_finalize(params)
      gw_doit = .true.
c
c     -----------------------------------------------------------------
c
      end function
