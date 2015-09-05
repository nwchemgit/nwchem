      subroutine dimqm_EqmN(rtdb, geom, xyz, eqmn)
c      use constants
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "geom.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "msgids.fh"
c
c      Input Variables
      integer rtdb
      integer geom
      double precision xyz(3, nDIM)
      double precision eqmn(3, nDIM)
c
c      Local Variables
      integer mLo, mUp, mOnNode, id
      integer nNuc
      integer l_qmxyz, k_qmxyz
      integer l_qmchrg, k_qmchrg
      integer m, n
c
      id = ga_nodeid()
      if(id.eq.0 .and. ldebug) write(luout,*) "Start EqmN"
c
c      Pull in QM geometry data
      if(.not. geom_ncent(geom, nNuc))
     $  call errquit("dimqm_EqmN: geom_ncent", 1, GEOM_ERR)
c
c      Allocate memory for coords and charges
      if(.not.ma_push_get(mt_dbl,3*nNuc,'QM xyz',l_qmxyz,k_qmxyz))
     $  call errquit("dimqm_EqmN malloc k_qmxyz failed", 1, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nNuc,'QM chrg',l_qmchrg,k_qmchrg))
     $  call errquit("dimqm_EqmN malloc k_qmchrg failed", 1, MA_ERR)
c
c      Populate the arrays
      if(.not.geom_cart_coords_get(geom, dbl_mb(k_qmxyz)))
     $  call errquit("dimqm_EqmN coords get failed", 1, GEOM_ERR)
      if(.not.geom_cart_get_charges(geom,nNuc,dbl_mb(k_qmchrg)))
     $  call errquit("dimqm_EqmN charges get failed", 1, GEOM_ERR)
c
c      Determine load balance
      call pphilo(id, nDIM, mLo, mUp, mOnNode)
c
      do m = mLo, mUp ! Loop over DIM atoms on this node
        do n = 1, nNuc ! Loop over all QM atoms
c
c         Calculate potential
          call dimqm_fldpot(xyz(:,m), dbl_mb(k_qmxyz+3*(n-1)),
     $                      dbl_mb(k_qmchrg+n-1), eqmn(:,m))
        end do
      end do
c      write(luout,*) "EqmN:"
c      write(luout,*) eqmn
      call ga_dgop(msg_dimqm, eqmn, 3*nDIM, '+')
      call ga_sync()
c
c      Deallocate memory
      if(.not.ma_chop_stack(l_qmxyz))
     $  call errquit('dimqm_EqmN chop failed', 1, MA_ERR)
c
      if(id.eq.0 .and. ldebug) write(luout,*) "End EqmN"
      end subroutine dimqm_EqmN
