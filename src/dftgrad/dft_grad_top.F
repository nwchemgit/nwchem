      subroutine dft_grad_top(rtdb)
      implicit none
c
      integer rtdb
c
      logical nwdft, grad_dft
      external nwdft, grad_dft
c
      if (.not. nwdft(rtdb))
     $     call errquit('dft_gradient: nwdft failed', 0)
      if (.not. grad_dft(rtdb))
     $     call errquit('dft_gradient: dft_gradients failed', 0)
c
      end
