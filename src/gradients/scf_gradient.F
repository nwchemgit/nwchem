      logical function scf_gradient(rtdb)
c$Id: scf_gradient.F,v 1.7 2000-01-05 00:25:49 windus Exp $
      implicit none
#include "rtdb.fh"
c
      integer rtdb
c
      character*32 theory
      logical mcscf, scf, gradients
      external mcscf, scf, gradients
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('(mc)scf gradients: no theory found',0)
c
      if (theory .eq. 'scf') then
         if (.not. scf(rtdb))
     $        call errquit('scf_gradient: scf energy failed', 0)
         call util_print_push
         call util_print_rtdb_load(rtdb,'scf')
      else if (theory .eq. 'mcscf') then
         if (.not. mcscf(rtdb))
     $        call errquit('mcscf_gradient: mcscf energy failed', 0)
         call util_print_push
         call util_print_rtdb_load(rtdb,'mcscf')
      else
         call errquit('(mc)scf gradients: unknown theory',0)
      endif
c
      if (.not. gradients(rtdb))
     $     call errquit('(mc)scf_gradient: gradients failed', 0)
c
      call util_print_pop
c
      scf_gradient = .true.
c
      end
