      subroutine dft_grad_top(rtdb)
*
* $Id: dft_grad_top.F,v 1.8 1997-11-04 10:12:11 d3e129 Exp $
*
      implicit none
c
      integer rtdb
c
      logical nwdft, grad_dft
      external nwdft, grad_dft
c
      if (.not. nwdft(rtdb))
     $     call errquit('dft_gradient: nwdft failed', 0)
c      write (6,*) 'dft_grad_top: before grad_dft'

      if (.not. grad_dft(rtdb))
     $     call errquit('dft_gradient: dft_gradients failed', 0)
c
      end
