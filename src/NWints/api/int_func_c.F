c $Id: int_func_c.F,v 1.1 1997-04-23 22:47:42 d3e129 Exp $
*
c:tex-% this is an internal API routine
c:tex-\subsubsection{int\_nospherical\_check}
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
#include "basP.fh"
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
        call errquit(' int_nospherical_check: error ',911)
      endif
      end
c:tex-% this is an internal API routine
c:tex-\subsubsection{int\_nogencont\_check}
c:tex-This routine stubs out routines that are not ready for 
c:tex-general contraction basis functions by forcing an error 
c:tex-condition.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_nogencont_check(basisin,shell,tine)
c:tex-\end{verbatim}
c
c routine to stop if basis/shell has general contractions
c
      implicit none
#include "stdio.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c:tex-\begin{verbatim}
      integer basisin     ! [input] basis set handle
      integer shell       ! [input] lexical contraction index
      character*(*) tine  ! [input] routine identifier
c:tex-\end{verbatim}
c
      integer basis
      logical okay
      integer myucont
      integer my_gen
c
#include "bas_ibs_sfn.fh"
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      myucont = (sf_ibs_cn2ucn(shell,basis))
      my_gen  = infbs_cont(Cont_Ngen,myucont,basis)
      okay = my_gen.lt.2
      if (.not.okay) then
        write(luout,*)' basis set handle ',basisin
        write(luout,*)' basis set name   ',
     &      bs_name(basis)(1:len_bs_name(basis))
        write(luout,*)' lexical shell index is         :',shell
        write(luout,*)' unique basis set shell index is:',myucont
        write(luout,*)' routine identified by :',tine
        write(luout,*)' has been called with a general contracted',
     &      ' basis set'
        write(luout,*)' this particular routine does not have this ',
     &      'functionality yet'
        call errquit(' int_nogencont_check: error ',911)
      endif
      end
c:tex-% this is an internal API routine
c:tex-\subsubsection{int\_nospshell\_check}
c:tex-This routine stubs out routines that are not ready for 
c:tex-sp shells type basis functions by forcing an error 
c:tex-condition.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_nospshell_check(basisin,shell,tine)
c:tex-\end{verbatim}
c
c routine to stop if basis/shell has sp shells
c
      implicit none
#include "stdio.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c:tex-\begin{verbatim}
      integer basisin     ! [input] basis set handle
      integer shell       ! [input] lexical contraction index
      character*(*) tine  ! [input] routine identifier
c:tex-\end{verbatim}
c
      integer basis
      logical okay
      integer myucont
      integer my_type
c
#include "bas_ibs_sfn.fh"
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      myucont = (sf_ibs_cn2ucn(shell,basis))
      my_type  = infbs_cont(Cont_Type,myucont,basis)
      okay = my_type.ge.0
      if (.not.okay) then
        write(luout,*)' basis set handle ',basisin
        write(luout,*)' basis set name   ',
     &      bs_name(basis)(1:len_bs_name(basis))
        write(luout,*)' lexical shell index is         :',shell
        write(luout,*)' unique basis set shell index is:',myucont
        write(luout,*)' routine identified by :',tine
        write(luout,*)' has been called with a general contracted',
     &      ' basis set'
        write(luout,*)' this particular routine does not have this ',
     &      'functionality yet'
        call errquit(' int_nospshell_check: error ',911)
      endif
      end
