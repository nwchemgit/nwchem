      subroutine ccsd_get_param(rtdb, maxit,convi,iprt)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
c     
      integer rtdb, maxit, convi, iprt
c
      double precision mone
      parameter (mone = -1.0d0)
c     
      if (.not. rtdb_get(rtdb, 'ccsd:maxiter',
     $     MT_INT, 1, maxit)) maxit = 20
      if (.not. rtdb_get(rtdb, 'ccsd:thresh',
     $     MT_INT, 1, convi)) convi = 8
      if (.not. rtdb_get(rtdb, 'ccsd:print',
     $     MT_INT, 1, iprt)) iprt = 0
c
c
      end
