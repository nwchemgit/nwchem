      logical function mcscf_gradient(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb
c     
      logical mcscf, gradients
      external mcscf, gradients
c
      if (.not. mcscf(rtdb))
     $     call errquit('mcscf_gradient: mcscf energy failed', 0,
     &       CALC_ERR)
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'mcscf')
      if (.not. gradients(rtdb))
     $     call errquit('mcscf_gradient: gradients failed', 0,
     &       CALC_ERR)
      call util_print_pop
c
      mcscf_gradient = .true.
c
      end
      logical function scf_gradient(rtdb)
c$Id$
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      logical scf, gradients
      external scf, gradients
c
      if (.not. scf(rtdb))
     $     call errquit('scf_gradient: scf energy failed', 0,
     &       CALC_ERR)
      call util_print_push
      call util_print_rtdb_load(rtdb,'scf')
      if (.not. gradients(rtdb))
     $     call errquit('scf_gradient: gradients failed', 0,
     &       CALC_ERR)
      call util_print_pop
c
      scf_gradient = .true.
c
      end
