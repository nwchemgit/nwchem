      subroutine dft_grad_top(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
c
      integer rtdb
c
      logical nwdft, grad_dft
      external nwdft, grad_dft
c
      if (.not. nwdft(rtdb))
     $     call errquit('dft_gradient: nwdft failed', 0, UNKNOWN_ERR)
c      write (6,*) 'dft_grad_top: before grad_dft'

      if (.not. grad_dft(rtdb))
     $     call errquit('dft_gradient: dft_gradients failed', 0,
     &       UNKNOWN_ERR)
c
      end
