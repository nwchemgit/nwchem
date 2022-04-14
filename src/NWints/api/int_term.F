c $Id$
C> \ingroup nwint
C> @{
C>
C> \brief The main integral termination routine
C>
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
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "apiP.fh"
#include "candoP.fh"
#include "stdio.fh"
c::functions
      logical spcart_terminate 
      external spcart_terminate 
      logical hnd_terminate 
      external hnd_terminate 
c
      if (init_int .eq. 0) then
        write(luout,*)' warning: int_init has not been called '
      endif
      call int_mem_zero()
c      
* def u=f d=f -> f.and.!f -> f -> e = t
* no txs u=t d=f -> t.and.!f -> t -> e = f
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_terminate()
      endif
      if (.not.(user_cando_sim.and.(.not.def_cando_sim))) then
        call nwcsim_terminate()
      endif
      if (.not.spcart_terminate()) then
        call errquit('int_terminate: spcart_terminate failed',0, 
     &             INT_ERR)
      endif
      if (.not.hnd_terminate()) then
        call errquit('int_terminate: hnd_terminate failed', 555,
     &            INT_ERR)
      endif
      if (any_ecp) then
        call int_ecp_terminate()
        ecp_bsh = 0
      endif
c
c   Terminate relativistic integrals. This is an unconditional call
c   because the call to int_rel_init is also unconditional, and the
c   flag init_rel_int must be set to .false.
c
      if (dyall_mod_dir) call int_rel_term
c
c   Terminate Douglas-Kroll integrals. The fitting basis set will
c   be destroyed and the handle put to -1.
c
      if (doug_kroll) then
       if (.not. bas_destroy(dk_fitbas)) call errquit(
     &   'int_terminate: failed to destroy Douglas-Kroll fitting basis'
     &   ,0, INT_ERR)
       dk_fitbas = -1
      endif
c
c     terminate uncontracted basis used in the X2C code
      if (exact_2c) then
        if (.not. bas_destroy(x2c_bas)) call errquit(
     &    'int_terminate: failed to destroy x2c uncontracted basis'
     &    ,0, INT_ERR)
        x2c_bas = -1
      endif
c
      end
C>
C> \brief The derivative integral termination routine
C>
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
#include "bas.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "candoP.fh"
c::functions
      logical spcart_terminate 
      external spcart_terminate 
      logical hnd_terminate 
      external hnd_terminate 
c
      if (init_intd.eq.0) then
        write(luout,*)' warning: intd_init has not been called '
      endif
      call int_mem_zero()
c.. parameters for int derivative (hf2d)
      intd_memthresh = 0
      numd_tot       = 0
      numd_okay      = 0
      numd_red       = 0
c
      if (.not.spcart_terminate()) then
        call errquit('intd_terminate: spcart_terminate failed',0,
     &         INT_ERR)
      endif
c
      if (.not.hnd_terminate()) then
        call errquit('intd_terminate: hnd_terminate failed', 555,
     &         INT_ERR)
      endif
c
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_terminate()
      endif
      if (.not.(user_cando_sim.and.(.not.def_cando_sim))) then
        call nwcsim_terminate()
      endif
*:rak:this should not be here      call int_app_unset_no_texas(rtdbIused)
c
      if (any_ecp) then
        call int_ecp_terminate()
        ecp_bsh = 0
      endif
c
      if (dyall_mod_dir) call int_rel_term
c
c   Terminate Douglas-Kroll integrals. The fitting basis set will
c   be destroyed and the handle put to -1.
c
      if (doug_kroll) then
       if (.not. bas_destroy(dk_fitbas)) call errquit(
     &   'int_terminate: failed to destroy Douglas-Kroll fitting basis'
     &   ,0, INT_ERR)
       dk_fitbas = -1
      endif
c
      end
C>
C> \brief The 2nd derivative integral termination routine
C>
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
#include "errquit.fh"
#include "candoP.fh"
#include "stdio.fh"
c::functions
      logical spcart_terminate 
      external spcart_terminate 
      logical hnd_terminate 
      external hnd_terminate 
c
      if (init_intdd.eq.0) then
        write(luout,*)' warning: intdd_init has not been called '
      endif
      call int_mem_zero()
      intd_memthresh = 0
c
      if (.not.spcart_terminate()) then
        call errquit('intdd_terminate: spcart_terminate failed',0,
     &             INT_ERR)
      endif
c
      if (.not.hnd_terminate()) then
        call errquit('intd_terminate: hnd_terminate failed', 555,
     &          INT_ERR)
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
C> @}
