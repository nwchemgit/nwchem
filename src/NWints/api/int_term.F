c:tex-
c:tex-%API Initialization and Termination Routines
c:tex-\subsection{int\_terminate}
c:tex-This is the main termination routine for integrals.
c:tex-After this call the INT-API is ready for re-initialization.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_terminate()
c:tex-\end{verbatim}
c:tex-No formal arguments
c:tex-
c $Id: int_term.F,v 1.9 1997-06-06 21:07:34 pg481 Exp $
      implicit none
#include "apiP.fh"
#include "candoP.fh"
c::functions
      logical spcart_terminate 
      external spcart_terminate 
c
      if (init_int .eq. 0) then
        write(6,*)' warning: int_init has not been called '
      endif
      call int_mem_zero()
c      
      call sp_terminate()
* def u=f d=f -> f.and.!f -> f -> e = t
* no txs u=t d=f -> t.and.!f -> t -> e = f
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_terminate()
      endif
      if (.not.spcart_terminate()) then
        call errquit('int_terminate: spcart_terminate failed',0)
      endif
      if (any_ecp) then
        call int_ecp_terminate()
        ecp_bsh = 0
      endif
      end
      subroutine sp_terminate()
*_not_used?:      implicit none
*_not_used?:#include "apispP.fh"
*_not_used?:#include "mafdecls.fh"
*_not_used?:c
*_not_used?:      if (no_sp_conts) return
*_not_used?:c
*_not_used?:      write(6,*)':::: de-allocating sp blocks of memory'
*_not_used?:      if (.not.MA_free_heap(hszspa))
*_not_used?:     &       call errquit('sp_terminate: free of hszspa failed',911)
*_not_used?:      write(6,*)'::done:: de-allocating sp blocks of memory'
*_not_used?:      if (.not.MA_free_heap(hszspb))
*_not_used?:     &       call errquit('sp_terminate: free of hszspb failed',911)
*_not_used?:      write(6,*)'::done:: de-allocating sp blocks of memory'
*_not_used?:      if (.not.MA_free_heap(hszspc))
*_not_used?:     &       call errquit('sp_terminate: free of hszspc failed',911)
*_not_used?:      write(6,*)'::done:: de-allocating sp blocks of memory'
*_not_used?:      if (.not.MA_free_heap(hszspd))
*_not_used?:     &       call errquit('sp_terminate: free of hszspd failed',911)
*_not_used?:      write(6,*)'::done:: de-allocating sp blocks of memory'
*_not_used?:      if (.not.MA_free_heap(hspscat))
*_not_used?:     &       call errquit('sp_terminate: free of hspscat failed',911)
*_not_used?:      write(6,*)'::done:: de-allocating sp blocks of memory'
*_not_used?:c
      end
c:tex-
c:tex-%API Initialization and Termination Routines
c:tex-\subsection{intd\_terminate}
c:tex-This is the main termination routine for integral
c:tex-derivatives.
c:tex-After this call the INT-API is ready for re-initialization.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intd_terminate()
c:tex-\end{verbatim}
c:tex-No formal arguments
c:tex-
      implicit none
#include "apiP.fh"
#include "candoP.fh"
c::functions
      logical spcart_terminate 
      external spcart_terminate 
c
      if (init_intd.eq.0) then
        write(6,*)' warning: intd_init has not been called '
      endif
      call int_mem_zero()
c.. parameters for int derivative (hf2d)
      intd_memthresh = 0
      numd_tot       = 0
      numd_okay      = 0
      numd_red       = 0
c
      call sp_terminate()
c
      if (.not.spcart_terminate()) then
        call errquit('int_terminate: spcart_terminate failed',0)
      endif
c
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_terminate()
      endif
      call int_app_unset_no_texas(rtdbIused)
c
      if (any_ecp) then
        call int_ecp_terminate()
        ecp_bsh = 0
      endif
c
      end



