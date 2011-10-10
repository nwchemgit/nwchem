      subroutine dft_get_conv_info(rtdb)
C$Id$
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "cscf.fh"
c     
      integer rtdb
c
      double precision mone
      parameter (mone = -1.0d0)
c     
c     if (.not.xc_gotxc()) then
c       if (.not. rtdb_get(rtdb, 'scf:maxiter', MT_INT, 1, maxiter))
c    &  then
c          if (scftype .ne. 'UHF') then
c             maxiter = 30
c          else
c             maxiter = 30
c          endif
c       endif
c       if (.not.rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, gnorm_tol)) 
c    $       gnorm_tol = 1.d-4
c
c       Ensure that the default integral selection is sufficient
c       for the request accuracy of the SCF.  However, allow user
c       override.
c
c       if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) 
c    $       tol2e = min(1.0d-7,gnorm_tol*1d-2)
c
c       if (rtdb_get(rtdb, 'scf:level shift info', MT_DBL, 6,shifts))
c    &  then
c          ouser_changed_conv = .true.
c       else
c          call dfill(6, mone, shifts, 1)
c       endif
c       if (rtdb_get(rtdb, 'scf:full hessian switch', MT_DBL, 1, 
c    $       nr_gswitch)) then
c          ouser_changed_conv = .true.
c       else 
c          nr_gswitch = 0.1d0
c       endif
c     endif
c
c     Apply defaults
c
      if (shifts(1) .eq. -1.0d0) shifts(1) = 5.0d0
      if (shifts(2) .eq. -1.0d0) shifts(2) = 0.5d0
      if (shifts(3) .eq. -1.0d0) shifts(3) = 0.0d0
      if (shifts(4) .eq. -1.0d0) shifts(4) = 0.0d0
      if (shifts(5) .eq. -1.0d0) shifts(5) = 0.0d0
      if (shifts(6) .eq. -1.0d0) shifts(6) = 0.0d0
c
      end
