c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Retrieve maximum buffer and scratch array sizes
C>
C> This routine returns the maximum buffer and scratch array
C> sizes for both one electron and two electron (4 center) 
C> integrals.
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer max1e       !< [Output] max 1e buffer size
      integer maxg        !< [Output] max 2e4c buffer size
      integer mscratch_1e !< [Output] max scr for 1e ints
      integer mscratch_2e !< [Output] max scr for 2e ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem')) then
        max1e       = isz_1e
        maxg        = isz_2e4c
        mscratch_1e = mem_1e
        mscratch_2e = mem_2e4c
      else
        call errquit('int_mem: int_init was not called' ,0, INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array size
C> for 1-electron integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer max1e       !< [Output] max 1e buffer size
      integer mscratch_1e !< [Output] max scr for 1e ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_1e')) then
        max1e       = isz_1e
        mscratch_1e = mem_1e
      else
        call errquit('int_mem_1e: int_init was not called' ,0, INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array size
C> for 2-electron integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg        !< [Output] max 2e4c buffer size
      integer mscratch_2e !< [Output] max scr for 2e ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_2e4c')) then
        maxg        = isz_2e4c
        mscratch_2e = mem_2e4c
      else
        call errquit('int_mem_2e4c: int_init was not called',0, INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array size
C> for 1-electron Hamiltonian integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxh1       !< [Output] max buffer size for H1 ints
      integer mscratch_h1 !< [Output] max scr size for H1 ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_h1')) then
        maxh1       = isz_1e
        mscratch_h1 = mem_1e_h1
      else
        call errquit('int_mem_h1: int_init was not called',0, INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array size
C> for 2-electron 3-center integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg           !< [Output] max buf size for 2e3c ints
      integer mscratch_2e3c  !< [Output] max scr size for 2e3c ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_2e3c')) then
        maxg          = isz_2e3c
        mscratch_2e3c = mem_2e3c
      else
        call errquit('int_mem_2e3c: int_init was not called',0, INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array size
C> for 2-electron 2-center integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg          !< [Output] max buf size for 2e2c ints
      integer mscratch_2e2c !< [Output] max scr size for 2e2c ints
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_2e2c')) then
        maxg          = isz_2e2c
        mscratch_2e2c = mem_2e2c
      else
        call errquit('int_mem_2e2c: int_init was not called',0, INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array size
C> for 1-electron 3-center integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxbuf       !< [Output] max buf size for 3 center ov
      integer mscratch_3ov !< [Output] max scr size for 3 center ov
c:tex-\end{verbatim}
c
      if (int_chk_init('int_mem_3ov')) then
        maxbuf       = isz_3ov
        mscratch_3ov = mem_3ov
      else
        call errquit('int_mem_3ov: int_init was not called' ,0, INT_ERR)
      endif
      end
C>
C> \brief Print buffer and scratch array sizes
C>
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
      integer ii, int_kbsize
      int_kbsize(ii) = ((ii*8)/1024)
c
      write(luout,*)' memory parameters for integral api:'
      write(luout,10000)
c
      if (init_int.eq.1)
     &    write(luout,*)' int_init was called '
      if (init_intd.eq.1)
     &    write(luout,*)' intd_init was called '

      write(luout,10001)' mem_2e4c  : ',mem_2e4c,int_kbsize(mem_2e4c)
      write(luout,10001)' mem_1e    : ',mem_1e,int_kbsize(mem_1e)
      write(luout,10001)' isz_2e4c  : ',isz_2e4c,int_kbsize(isz_2e4c)
      write(luout,10001)' isz_1e    : ',isz_1e,int_kbsize(isz_1e)
      write(luout,10001)' mem_1e_h1 : ',mem_1e_h1,int_kbsize(mem_1e_h1)
      write(luout,10001)' isz_2e3c  : ',isz_2e3c,int_kbsize(isz_2e3c)
      write(luout,10001)' isz_2e2c  : ',isz_2e2c,int_kbsize(isz_2e2c)
      write(luout,10001)' mem_2e3c  : ',mem_2e3c,int_kbsize(mem_2e3c)
      write(luout,10001)' mem_2e2c  : ',mem_2e2c,int_kbsize(mem_2e2c)
      write(luout,10001)' isz_3ov   : ',isz_3ov,int_kbsize(isz_3ov)
      write(luout,10001)' mem_3ov   : ',mem_3ov,int_kbsize(mem_3ov)
      write(luout,10001)' iszb_2e4c : ',iszb_2e4c,int_kbsize(iszb_2e4c)
      write(luout,10001)' memb_2e4c : ',memb_2e4c,int_kbsize(memb_2e4c)
10001 format(a,i12,1x,'64-bit-words',i20,1x,'Kbytes')
10000 format(//)
      end
C>
C> \brief Retrieve maximum buffer and scratch array sizes for the
C> blocking 2-electron integrals
C>
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
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer maxg        !< [Output] max buf size for blocked 2e4c ints
      integer mscratch_2e !< [Output] max scr size for blocked 2e4c ints
c:tex-\end{verbatim}
c
      if (int_chk_init('intb_mem_2e4c')) then
        maxg        = iszb_2e4c
        mscratch_2e = memb_2e4c
      else
        call errquit('intb_mem_2e4c: int_init was not called' ,0,
     &       INT_ERR)
      endif
      end
C>
C> \brief Retrieve the maximum buffer and scratch array sizes for 
C> the dipole/multipole integrals
C>
      subroutine int_mem_dipole(maxbuf, maxscr,
     &    basisin_i, basisin_j, lval_max)
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxbuf    !< [Output] maximum buffer size for dipole routines
      integer maxscr    !< [Output] maximum scratch size for dipole routines
      integer basisin_i !< [Input] basis set for part of multipole ints
      integer basisin_j !< [Input] other basis set for part of multipole ints
      integer lval_max  !< [Input] maximum order multipole used
*
      integer maxang_i
      integer maxang_j
*
      if (.not.int_chk_init('int_mem_dipole')) then
        call errquit('int_mem_dipole: int_init was not called' ,0,
     &       INT_ERR)
      endif
      maxscr = mem_3ov  ! should always work for dipole integrals
      if (.not.bas_high_angular(basisin_i,maxang_i))
     &    call errquit
     &    ('int_mem_dipole: bas_high_angular failed',911, INT_ERR)
      if (.not.bas_high_angular(basisin_j,maxang_j))
     &    call errquit
     &    ('int_mem_dipole: bas_high_angular failed',911, INT_ERR)
      maxbuf =          ((lval_max+1)*(lval_max+2)*(lval_max+3))/6
      maxbuf = maxbuf * ((maxang_i+1)*(maxang_i+2))/2
      maxbuf = maxbuf * ((maxang_j+1)*(maxang_j+2))/2
      end
C>
C> \brief Retrieve the maximum buffer and scratch array sizes for 
C> the dipole/multipole integrals
C>
C> This routine is identical to `int_mem_dipole`.
C>
      subroutine int_mem_dipole_l(maxbuf, maxscr,
     &    basisin_i, basisin_j, lval)
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxbuf    !< [Output] maximum buffer size for dipole routines
      integer maxscr    !< [Output] maximum scratch size for dipole routines
      integer basisin_i !< [Input] basis set for part of multipole ints
      integer basisin_j !< [Input] other basis set for part of multipole ints
      integer lval      !< [Input] order of multipole for the given call
*
      integer maxang_i
      integer maxang_j
*
      if (.not.int_chk_init('int_mem_dipole_l')) then
        call errquit('int_mem_dipole_l: int_init was not called' ,0,
     &       INT_ERR)
      endif
      maxscr = mem_3ov  ! should always work for dipole integrals
      if (.not.bas_high_angular(basisin_i,maxang_i))
     &    call errquit
     &    ('int_mem_dipole: bas_high_angular failed',911, INT_ERR)
      if (.not.bas_high_angular(basisin_j,maxang_j))
     &    call errquit
     &    ('int_mem_dipole: bas_high_angular failed',911, INT_ERR)
      maxbuf =          ((lval+1)*(lval+2)*(lval+3))/6
      maxbuf = maxbuf * ((maxang_i+1)*(maxang_i+2))/2
      maxbuf = maxbuf * ((maxang_j+1)*(maxang_j+2))/2
      end
C>
C> \brief Retrieve the maximum buffer and scratch array sizes for 
C> the dipole/multipole derivative integrals
C>
      subroutine intd_mem_dipole(maxbuf, maxscr,
     &    basisin_i, basisin_j, lval_max)
      implicit none
* computes memory for derivative dipole integrals
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxbuf    !< [Output] maximum buffer size for dipole routines
      integer maxscr    !< [Output] maximum scratch size for dipole routines
      integer basisin_i !< [Input] basis set for part of multipole ints
      integer basisin_j !< [Input] other basis set for part of multipole ints
      integer lval_max  !< [Input] maximum order multipole used
*
      if (.not.int_chk_init('intd_mem_dipole')) then
        call errquit('intd_mem_dipole: int_init was not called' ,0,
     &       INT_ERR)
      endif
      call int_mem_dipole(maxbuf,maxscr,basisin_i,basisin_j,lval_max)
      maxbuf = maxbuf * 9  ! (block,xyz,3) 3=atom1, atom2, multipole center
      end
C>
C> \brief Retrieve the maximum buffer and scratch array sizes for 
C> the dipole/multipole derivative integrals
C>
C> This routine is identical to `intd_mem_dipole`.
C>
      subroutine intd_mem_dipole_l(maxbuf, maxscr,
     &    basisin_i, basisin_j, lval)
      implicit none
#include "apiP.fh"
#include "errquit.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxbuf    !< [Output] maximum buffer size for dipole routines
      integer maxscr    !< [Output] maximum scratch size for dipole routines
      integer basisin_i !< [Input] basis set for part of multipole ints
      integer basisin_j !< [Input] other basis set for part of multipole ints
      integer lval      !< [Input] multipole order for the given call
*
      if (.not.int_chk_init('intd_mem_dipole_l')) then
        call errquit('intd_mem_dipole_l: int_init was not called' ,0,
     &       INT_ERR)
      endif
      call int_mem_dipole_l(maxbuf,maxscr,basisin_i,basisin_j,lval)
      maxbuf = maxbuf * 9  ! (block,xyz,3) 3=atom1, atom2, multipole center
      end
C> @}
