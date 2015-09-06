      subroutine matvecReal(rtdb, nlen, g_vec, g_rhs, xyz, ldbg)
c   ==================================================================
c   Performs the operation A*x = b, where x and b are real.
c
c   In  : rtdb      - Handle to the RTDB
c         nlen      - Length of the work array
c   ==================================================================
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input Variables
      integer rtdb
      integer nlen
      integer g_vec
      integer g_rhs
      logical ldbg
      double precision xyz(3,nDIM)
c
c     Variables from dimqm common blocks
c      integer nDIM, nDIMTypes, nOrder, g_DIM_diag
c
c     Local Variables
      integer ntype, mtype
      integer m, n, ms, me, ns, ne, id
      integer mLo, mUp, mOnNode, nodes
      double precision time
      double precision sPol(nDIMTypes)
      double precision vec(nOrder)
      integer tx(nDIM) ! Type Index
      double precision diagS(nOrder)
      double precision fld(3)
      logical stat
      integer lo, hi
      call ga_sync()
c
c      Determine node ID
      id = ga_nodeid()
      if(ldbg .and. id.eq.0) write(LuOut,*)
     $   "Start MatvecReal routine"
      time = util_timer()
c
c      Zero right-hand side
      call ga_zero(g_rhs)
c
c      Determine chunk to be done on this node
      call pphilo(id, nDIM, mLo, mUp, mOnNode)
c
c      Pull values from RTDB and GAs
      stat = rtdb_get(rtdb, 'dimpar:sPol', mt_dbl, nDIMTypes, sPol)
      stat = rtdb_get(rtdb, 'dimpar:typeindex', mt_int, nDIM, tx)
      call ga_get(g_DIM_diag, 1, nOrder, 1, 1, diagS, 1)
      call ga_get(g_vec, 1, nOrder, 1, 1, vec, 1)
      call ga_init_fence()
c
c      Loop over atoms on this node
      do m = mLo, mUp
        mtype = tx(m)
        fld = ZERO
        ms = 3 * (m - 1) + 1
        me = 3 * (m - 1) + 3
        do n = 1, nDIM
          ntype = tx(n)
          ns = 3 * (n - 1) + 1
          ne = 3 * (n - 1) + 3
          fld(1:3) = fld(1:3) + atomPairField(m, n,
     $                                     xyz(:,m),     ! Coordinates to atom m
     $                                     xyz(:,n),     ! Coordinates to atom n
     $                                     sPol(mtype),  ! Static polarizability of atom m
     $                                     sPol(ntype),  ! Static polarizability of atom n
     $                                     diagS(ns:ne), ! Diagonal elements
     $                                     vec(ns:ne),   ! Vector x in A*x = b
     $                                     ldbg)         ! debug
        end do
        call ga_acc(g_rhs, ms, me, 1, 1, fld, 1, ONE)
      end do
      call ga_fence()
      time = util_timer() - time
      if(ldbg .and. id.eq.0) write(LuOut,*)
     $   "End MatvecReal Routine"

      contains 
      function atomPairField(m, n, mxyz, nxyz, mpol, npol,
     $                        diag, vec, debug) result (fld)
      implicit none
#include "stdio.fh"
#include "dimqm_constants.fh"
      integer, intent(in) :: m, n
      double precision, intent(in) :: mxyz(3), nxyz(3) ! Coordinates of the atoms
      double precision, intent(in) :: mpol, npol ! Static polarizability of the atoms
      double precision, intent(in) :: diag(3) ! Diagonal elements
      double precision, intent(in) :: vec(3)  ! Vector being multiplied
      logical, intent(in) :: debug
c
c     Output
      double precision fld(3)
c
c     Local variables
      double precision r(3)
      double precision t2(3,3)
      double precision dist, invdist
      interface
        function t2r(r, dist, invdist, a1, a2, pol1, pol2)
          integer, intent(in) :: a1, a2
          double precision, intent(in) :: dist, invdist
          double precision, intent(in) :: pol1, pol2
          double precision t2r(3,3)
          double precision, intent(in) :: r(3)
        end function t2r
      end interface
c
c   Diagonal element
      if(m .eq. n) then
        fld(1:3) = vec(1:3) * diag(1:3)
      else
C   Off diagonal
        r(1:3) = nxyz(1:3) - mxyz(1:3)
        dist = SQRT(DOT_PRODUCT(r, r))
        invdist = ONE / dist
        t2 = t2r(r, dist, invdist, m, n, mpol, npol)
        fld(1:3) = - MATMUL(t2, vec)
      end if
c      write(LuOut,*) "Fld:", fld
      end function atomPairField
      end subroutine matvecReal
    
