c $Id: int_mem.F,v 1.9 1997-05-29 00:36:09 d3e129 Exp $
*
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for both one electron and two electron (4 center) 
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem(max1e, maxg, mscratch_1e, mscratch_2e)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer max1e       ! [output] max 1e buffer size
      integer maxg        ! [output] max 2e4c buffer size
      integer mscratch_1e ! [output] max scr for 1e ints
      integer mscratch_2e ! [output] max scr for 2e ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem')) then
        max1e       = isz_1e
        maxg        = isz_2e4c
        mscratch_1e = mem_1e
        mscratch_2e = mem_2e4c
      else
        call errquit('int_mem: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_1e}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for one electron (2 center)
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_1e(max1e, mscratch_1e)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer max1e       ! [output] max 1e buffer size
      integer mscratch_1e ! [output] max scr for 1e ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_1e')) then
        max1e       = isz_1e
        mscratch_1e = mem_1e
      else
        call errquit('int_mem_1e: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_2e4c}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for two electron (4 center)
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_2e4c(maxg, mscratch_2e)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg        ! [output] max 2e4c buffer size
      integer mscratch_2e ! [output] max scr for 2e ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_2e4c')) then
        maxg        = isz_2e4c
        mscratch_2e = mem_2e4c
      else
        call errquit('int_mem_2e4c: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_h1}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for one electron hamiltonian 
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_h1(maxh1, mscratch_h1)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxh1       ! [output] max buffer size for H1 ints
      integer mscratch_h1 ! [output] max scr size for H1 ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_h1')) then
        maxh1       = isz_1e
        mscratch_h1 = mem_1e_h1
      else
        call errquit('int_mem_h1: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_2eec}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for two electron (3 center)
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_2e3c(maxg, mscratch_2e3c)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg           ! [output] max buf size for 2e3c ints
      integer mscratch_2e3c  ! [output] max scr size for 2e3c ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_2e3c')) then
        maxg          = isz_2e3c
        mscratch_2e3c = mem_2e3c
      else
        call errquit('int_mem_2e3c: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_2e2c}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for two electron (2 center)
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_2e2c(maxg, mscratch_2e2c)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg          ! [output] max buf size for 2e2c ints
      integer mscratch_2e2c ! [output] max scr size for 2e2c ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_2e2c')) then
        maxg          = isz_2e2c
        mscratch_2e2c = mem_2e2c
      else
        call errquit('int_mem_2e2c: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_3ov}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for one electron (3 center)
c:tex-integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_3ov(maxbuf, mscratch_3ov)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxbuf       ! [output] max buf size for 3 center ov
      integer mscratch_3ov ! [output] max scr size for 3 center ov
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_3ov')) then
        maxbuf       = isz_3ov
        mscratch_3ov = mem_3ov
      else
        call errquit('int_mem_3ov: int_init was not called' ,0)
      endif
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_mem\_print}
c:tex-This routine prints the maximum buffer and scratch array
c:tex-sizes for all known ``classes'' of integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_print()
c:tex-\end{verbatim}
      implicit none
#include "stdio.fh"
#include "apiP.fh"
c
      write(luout,*)' memory parameters for integral api:'
      write(luout,10000)
c
      if (init_int.eq.1)
     &    write(luout,*)' int_init was called '
      if (init_intd.eq.1)
     &    write(luout,*)' intd_init was called '

      write(luout,*)' mem_2e4c  : ',mem_2e4c
      write(luout,*)' mem_1e    : ',mem_1e
      write(luout,*)' isz_2e4c  : ',isz_2e4c
      write(luout,*)' isz_1e    : ',isz_1e
      write(luout,*)' mem_1e_h1 : ',mem_1e_h1
      write(luout,*)' isz_2e3c  : ',isz_2e3c
      write(luout,*)' isz_2e2c  : ',isz_2e2c
      write(luout,*)' mem_2e3c  : ',mem_2e3c
      write(luout,*)' mem_2e2c  : ',mem_2e2c
      write(luout,*)' isz_3ov   : ',isz_3ov
      write(luout,*)' mem_3ov   : ',mem_3ov
      write(luout,*)' iszb_2e4c : ',iszb_2e4c
      write(luout,*)' memb_2e4c : ',memb_2e4c
10000 format(//)
      end
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{intb\_mem\_2e4c}
c:tex-This routine returns the maximum buffer and scratch array
c:tex-sizes for two electron (4 center)
c:tex-integrals from the blocking interface.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intb_mem_2e4c(maxg, mscratch_2e)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg        ! [output] max buf size for blocked 2e4c ints
      integer mscratch_2e ! [output] max scr size for blocked 2e4c ints
c:tex-\end{verbatim}
c
      if (int_chk_init('intb_mem_2e4c')) then
        maxg        = iszb_2e4c
        mscratch_2e = memb_2e4c
      else
        call errquit('intb_mem_2e4c: int_init was not called' ,0)
      endif
      end
