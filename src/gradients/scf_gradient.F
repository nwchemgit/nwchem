      logical function scf_gradient(rtdb)
c$Id: scf_gradient.F,v 1.4 1996-07-09 22:31:33 vg038 Exp $
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
      scf_gradient = .true.
      end
