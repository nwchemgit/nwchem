      subroutine riscf_tidy ( geom, basis, riscf_basis )
C$Id: riscf_tidy.F,v 1.3 1996-02-07 00:41:31 vg038 Exp $
      implicit none

#include "geom.fh"
#include "bas.fh"
#ifdef RISCF_USE_DA
#include "da.fh"
#endif
c
c     this is just a wrapper for del_3c_ints, to free ALL arrays of 3c 
c     integrals
c
c     and destroys basis and geometry (and DAs)
c
      integer geom, basis, riscf_basis

      integer istat
      logical status, riscf_del_3c_ints

      status = riscf_del_3c_ints(-1)
      if (.not.status)
     $     call errquit('riscf_tidy: failed to free 3c ints', 0)

C     status = bas_destroy ( basis )
C     if (.not.status)
C     $     call errquit('riscf_tidy: failed to destroy ao basis', 0)

      status = bas_destroy ( riscf_basis )
      if (.not.status)
     $     call errquit('riscf_tidy: failed to destroy ri basis', 0)

C     status = geom_destroy ( geom )
C     if (.not.status)
C     $     call errquit('riscf_tidy: failed to destroy geometry', 0)

#ifdef RISCF_USE_DA
      istat = da_terminate()
#endif

      return
      end
