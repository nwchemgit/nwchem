      logical function bsemol(rtdb)
      implicit none
#include "rtdb.fh"
      integer rtdb
c
      logical bse_doit
      external bse doit
c
      logical status
c
      status = rtdb_cput(rtdb,'dft:theory', 1, 'dft')
      if (status) status = bse_doit(rtdb)
      bsemol = status
      call grid_cleanup(.false.)
c
      end function
c
c     -----------------------------------------------------------------
c
      logical function bse_doit(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "bse.fh"

      integer, intent(in) :: rtdb
      character(*), parameter :: pname = "bse_doit: "

      type(bse_params_t) :: params
      character*80 theory
      integer me
c
c     -----------------------------------------------------------------
c
c     Initialize
c
      bse_doit = .false.
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
      call bse_init(params)
c
      if (theory.eq."dft") then
        
        if (params%davidson) then
          call bse_davidson(params)
        elseif (params%lanczos) then
          call bse_lanczos(params)
        else
          call bse_analytic(params)
        endif
c
      elseif (theory.eq.'sodft') then
        call errquit(pname//"spin-orbit not implemented",0,0)
      else
        call errquit(pname//"invalid theory: "//trim(theory),0,0)
      endif
c
      call bse_finalize(params)
      bse_doit = .true.
c
c     -----------------------------------------------------------------
c
      end function
