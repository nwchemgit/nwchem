c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Print a block of 1-electron integrals
C>
C> Print a generic block of 1-electron integrals including labels.
C> The routine support two print levels:
C>
C> * 0 - Print the integrals and labels
C>
C> * 1 - Print shell information in addition to the integrals and labels
C>
c:tex-% This is part of the standard API
c:tex-\subsection{int\_pgen1e}
c:tex-This routine prints a generic one electron block of integrals.
c:tex-This requires the labels be generated and input to this routine.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_pgen1e(msg, i_basis, ish, j_basis, jsh, 
     &       ilab, jlab, l1e, Gen1e, print_level)
c:tex-\end{verbatim}
      implicit none
c
c generic one electron printing routine
c  print_level = 0 print only lables and integrals
c              = 1 print shell information etc. 
c
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "geobasmapP.fh"
c
c::external subroutines used
c errquit
c::function
      logical int_chk_sh
      external int_chk_sh
c::passed
c:tex-\begin{verbatim}
      character*(*) msg           !< [Input] informational message
      integer i_basis             !< [Input] basis set handle for i shell
      integer j_basis             !< [Input] basis set handle for j shell
      integer ish                 !< [Input] i shell
      integer jsh                 !< [Input] j shell 
      integer l1e                 !< [Input] number of integrals and labels
      integer ilab(l1e)           !< [Input] i basis functions labels
      integer jlab(l1e)           !< [Input] j basis functions labels
      double precision Gen1e(l1e) !< [Input] integrals to be printed
      integer print_level         !< [Input] print level to be used
c. . . . . . . . . . . . . .! print_level = 0 print lables and integrals
c. . . . . . . . . . . . . .!             = 1 also print shell info.
c:tex-\end{verbatim}
c
c::local
      integer ii
      integer stdout
c
      stdout = 6
c
      if (print_level.ge.1) then
        write(stdout,*)'pgen1e:',msg
        write(stdout,10000)i_basis,j_basis,ish,jsh
      endif
c
      do 00100 ii = 1,l1e
        write(stdout,10001)ilab(ii),jlab(ii),Gen1e(ii)
00100 continue
c
10000 format(1x,'<ibasis:jbasis>-<',i4,':',i4,'>',5x,
     &       '<ish:jsh>-<',i3,':',i3,'>')
10001 format('i = ',i4,' j = ',i4,5x,'1eval = ',1pd20.10)
c
      end
C> @}
