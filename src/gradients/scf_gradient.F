      subroutine scf_gradient(rtdb)
      implicit none
c
      integer rtdb
c
      logical scf, gradients
      external scf, gradients
c
      if (.not. scf(rtdb))
     $     call errquit('scf_gradient: scf failed', 0)
      if (.not. gradients(rtdb))
     $     call errquit('scf_gradient: gradients failed', 0)
c
      end
