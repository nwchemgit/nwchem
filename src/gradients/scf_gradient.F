      subroutine scf_gradient(rtdb)
C$Id: scf_gradient.F,v 1.2 1995-02-02 23:14:54 d3g681 Exp $
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
