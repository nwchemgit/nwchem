      subroutine int_terminate()
c $Id: int_term.F,v 1.2 1996-07-11 22:21:16 d3e129 Exp $
      implicit none
#include "apiP.fh"
c
      if (init_int .eq. 0) then
        write(6,*)' warning: int_init has not been called '
      endif
      mem_2e4c  = 0
      init_int  = 0
      init_intd = 0
      mem_1e    = 0
      isz_2e4c  = 0
      isz_1e    = 0
      mem_1e_h1 = 0
      isz_2e3c  = 0
      isz_2e2c  = 0
      mem_2e3c  = 0
      mem_2e2c  = 0
      isz_3ov   = 0
      mem_3ov   = 0
c      
      call sp_terminate()
      call texas_terminate()
      end
      subroutine sp_terminate()
      implicit none
#include "apispP.fh"
#include "mafdecls.fh"
c
      if (no_sp_conts) return
c
      if (.not.MA_free_heap(hszspa))
     &       call errquit('sp_terminate: free of hszspa failed',911)
      if (.not.MA_free_heap(hszspb))
     &       call errquit('sp_terminate: free of hszspb failed',911)
      if (.not.MA_free_heap(hszspc))
     &       call errquit('sp_terminate: free of hszspc failed',911)
      if (.not.MA_free_heap(hszspd))
     &       call errquit('sp_terminate: free of hszspd failed',911)
      if (.not.MA_free_heap(hspscat))
     &       call errquit('sp_terminate: free of hspscat failed',911)
c
      end
      subroutine intd_terminate()
      implicit none
#include "apiP.fh"
c
      if (init_intd.eq.0) then
        write(6,*)' warning: intd_init has not been called '
      endif
      mem_2e4c  = 0
      init_int  = 0
      init_intd = 0
      mem_1e    = 0
      isz_2e4c  = 0
      isz_1e    = 0
      mem_1e_h1 = 0
      isz_2e3c  = 0
      isz_2e2c  = 0
      mem_2e3c  = 0
      mem_2e2c  = 0
      isz_3ov   = 0
      mem_3ov   = 0
c.. parameters for int derivative (hf2d)
      intd_memthresh = 0
      numd_tot       = 0
      numd_okay      = 0
      numd_red       = 0
c
      call sp_terminate()
c
      end

