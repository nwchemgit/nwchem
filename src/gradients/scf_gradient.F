      logical function scf_gradient(rtdb)
c$Id: scf_gradient.F,v 1.5 1997-02-07 22:06:15 d3g681 Exp $
      implicit none
#include "rtdb.fh"
c
      integer rtdb
c
      character*20 theory
      logical mcscf, scf, gradients
      external mcscf, scf, gradients
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('(mc)scf gradients: no theory found',0)
c
      if (theory .eq. 'scf') then
         if (.not. scf(rtdb))
     $        call errquit('scf_gradient: scf energy failed', 0)
      else if (theory .eq. 'mcscf') then
         if (.not. mcscf(rtdb))
     $        call errquit('mcscf_gradient: mcscf energy failed', 0)
      else
         call errquit('(mc)scf gradients: unknown theory',0)
      endif
c
      if (.not. gradients(rtdb))
     $     call errquit('(mc)scf_gradient: gradients failed', 0)
c
      scf_gradient = .true.
c
      end
