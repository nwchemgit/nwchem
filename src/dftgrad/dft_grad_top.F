      subroutine dft_grad_top(rtdb)
      implicit real*8 (a-h,o-z) 
c
      integer rtdb
#include<quad.fh>
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
