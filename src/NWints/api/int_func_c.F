c $Id: int_func_c.F,v 1.7 2003-10-17 22:53:41 carlfahl Exp $
*
c:tex-% this is an internal API routine
c:tex-\subsection{int\_nospherical\_check}
c:tex-This routine stubs out routines that are not ready for 
c:tex-spherical basis functions by forcing an error condition.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_nospherical_check(basisin,tine)
c:tex-\end{verbatim}
c
c routine to stop if basis is spherical
c
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer basisin     ! [input] basis set handle
      character*(*) tine  ! [input] routine identifier
c:tex-\end{verbatim}
c
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      if (bas_spherical(basis)) then
        write(luout,*)' basis set handle ',basisin
        write(luout,*)' basis set name   ',
     &      bs_name(basis)(1:len_bs_name(basis))
        write(luout,*)' routine identified by :',tine
        write(luout,*)' has been called with a spherical basis set'
        write(luout,*)' this particular routine does not have this ',
     &      'functionality yet'
        call errquit(' int_nospherical_check: error ',911, INT_ERR)
      endif
      end
c:tex-% this is an internal API routine
c:tex-\subsection{int\_nogencont\_check}
c:tex-This routine stubs out routines that are not ready for 
c:tex-general contraction basis functions by forcing an error 
c:tex-condition.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_nogencont_check(basisin,tine)
c:tex-\end{verbatim}
c
c routine to stop if basis has any general contractions
c
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer basisin     ! [input] basis set handle
      character*(*) tine  ! [input] routine identifier
c:tex-\end{verbatim}
c
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      if (bas_any_gc(basis)) then
        write(luout,*)' basis set handle ',basisin
        write(luout,*)' basis set name   ',
     &      bs_name(basis)(1:len_bs_name(basis))
        write(luout,*)' routine identified by :',tine
        write(luout,*)' has been called with a general contracted',
     &      ' basis set'
        write(luout,*)' this particular routine does not have this ',
     &      'functionality yet'
        call errquit(' int_nogencont_check: error ',911, INT_ERR)
      endif
      end
c:tex-% this is an internal API routine
c:tex-\subsection{int\_nospshell\_check}
c:tex-This routine stubs out routines that are not ready for 
c:tex-sp shells type basis functions by forcing an error 
c:tex-condition.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_nospshell_check(basisin,tine)
c:tex-\end{verbatim}
c
c routine to stop if basis/shell has sp shells
c
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer basisin     ! [input] basis set handle
      character*(*) tine  ! [input] routine identifier
c:tex-\end{verbatim}
c
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      if (bas_any_sp_shell(basis)) then
        write(luout,*)' basis set handle ',basisin
        write(luout,*)' basis set name   ',
     &      bs_name(basis)(1:len_bs_name(basis))
        write(luout,*)' routine identified by :',tine
        write(luout,*)' has been called with a  basis set that ',
     &      'has sp shells'
        write(luout,*)' this particular routine does not have this ',
     &      'functionality yet'
        call errquit(' int_nospshell_check: error ',911, INT_ERR)
      endif
      end
c:tex-% this is an internal API routine
c:tex-\subsection{int\_bothsp\_gc\_check}
c:tex-This routine checks to see if the basis sets used
c:tex-have both sp shells/type basis functions and general
c:tex-contractions.  The 2e integral codes cannot handle this.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_bothsp_gc_check(basesin,nbas,tine)
c:tex-\end{verbatim}
c
c routine to stop if basis/shell has sp shells
c
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer nbas           ! [input] number of basis sets
      integer basesin(nbas)  ! [input] basis set handles
      character*(*) tine     ! [input] routine identifier
c:tex-\end{verbatim}

      integer bas1, basis1
      integer bas2, basis2
      logical badcombo
*
      do bas1 = 1,nbas
        do bas2 = 1,bas1
          basis1 = basesin(bas1) + BASIS_HANDLE_OFFSET
          basis2 = basesin(bas2) + BASIS_HANDLE_OFFSET
          badcombo = bas_any_sp_shell(basis1).and.bas_any_gc(basis2)
          if (badcombo) then
            write(luout,*)' calling routine: ',tine
            write(luout,*)
     &          ' basis set combination sp shells with general ',
     &          'contractions'
            write(luout,*)' basis 1 handle=',basesin(bas1),
     &          ' name:',bs_name(basis1)(1:len_bs_name(basis1))
            write(luout,*)'         any sp shells: ',
     &          bas_any_sp_shell(basis1)
            write(luout,*)'         any gen. con.: ',
     &          bas_any_gc(basis1)
            write(luout,*)' basis 2 handle=',basesin(bas2),
     &          ' name:',bs_name(basis2)(1:len_bs_name(basis2))
            write(luout,*)'         any sp shells: ',
     &          bas_any_sp_shell(basis2)
            write(luout,*)'         any gen. con.: ',
     &          bas_any_gc(basis2)
            call errquit('int_bothsp_gc_check: fatal error',911,
     &       INT_ERR)
          endif
        enddo
      enddo
      end
