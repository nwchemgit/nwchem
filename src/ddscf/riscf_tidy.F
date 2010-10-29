      subroutine riscf_tidy ( geom, basis, riscf_basis )
C$Id$
      implicit none
#include "errquit.fh"

#include "geom.fh"
#include "bas.fh"
#ifdef RISCF_USE_DRA
#include "dra.fh"
#endif
c
c     this is just a wrapper for del_3c_ints, to free ALL arrays of 3c 
c     integrals
c
c     and destroys basis and geometry (and DAs)
c
      integer geom, basis, riscf_basis

#ifdef RISCF_USE_DRA
      integer istat
#endif
      logical status, riscf_del_3c_ints

      status = riscf_del_3c_ints(-1)
      if (.not.status)
     $     call errquit('riscf_tidy: failed to free 3c ints', 0,
     &       UNKNOWN_ERR)

C     status = bas_destroy ( basis )
C     if (.not.status)
C     $     call errquit('riscf_tidy: failed to destroy ao basis', 0)

      status = bas_destroy ( riscf_basis )
      if (.not.status)
     $     call errquit('riscf_tidy: failed to destroy ri basis', 0,
     &       BASIS_ERR)

C     status = geom_destroy ( geom )
C     if (.not.status)
C     $     call errquit('riscf_tidy: failed to destroy geometry', 0)

#ifdef RISCF_USE_DRA
      istat = dra_terminate()
#endif

      return
      end
