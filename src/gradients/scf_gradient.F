      subroutine scf_gradient(rtdb)
c$Id: scf_gradient.F,v 1.3 1995-02-03 01:01:55 vg038 Exp $
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
