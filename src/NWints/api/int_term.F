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
c $Id: int_term.F,v 1.11 1998-09-22 16:14:51 d3e129 Exp $
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
      if (.not.spcart_terminate()) then
        call errquit('intd_terminate: spcart_terminate failed',0)
      endif
c
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_terminate()
      endif
*:rak:this should not be here      call int_app_unset_no_texas(rtdbIused)
c
      if (any_ecp) then
        call int_ecp_terminate()
        ecp_bsh = 0
      endif
c
      end
c:tex-
c:tex-%API Initialization and Termination Routines
c:tex-\subsection{intdd\_terminate}
c:tex-This is the main termination routine for integral
c:tex-second derivatives.
c:tex-After this call the INT-API is ready for re-initialization.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intdd_terminate()
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
      if (init_intdd.eq.0) then
        write(6,*)' warning: intdd_init has not been called '
      endif
      call int_mem_zero()
      intd_memthresh = 0
c
      if (.not.spcart_terminate()) then
        call errquit('intdd_terminate: spcart_terminate failed',0)
      endif
c
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_terminate()
      endif
c
      if (any_ecp) then
        call int_ecp_terminate()
        ecp_bsh = 0
      endif
c
      end



