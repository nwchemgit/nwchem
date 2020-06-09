      subroutine dimqm_broadcast()
c
c     Broadcasts all common variables collected on input so all processes have the correct value
c     
      implicit none
#include "dimqm.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer mitob1,mdtob1
c     Logicals
      call ga_brdcst(msg_dimqm, ldimqm, 1, 0)
!      call ga_brdcst(msg_dimqm, lcpim, 1, 0)
!      call ga_brdcst(msg_dimqm, lpim, 1, 0)
!      call ga_brdcst(msg_dimqm, lnopol, 1, 0)
!      call ga_brdcst(msg_dimqm, lnochar, 1, 0)
!      call ga_brdcst(msg_dimqm, lnocross, 1, 0)
      call ga_brdcst(msg_dimqm, lfreq, 1, 0)
!      call ga_brdcst(msg_dimqm, lfullgrid, 1, 0)
      call ga_brdcst(msg_dimqm, lclfld, 1, 0)
      call ga_brdcst(msg_dimqm, lrsp, 1, 0)
      call ga_brdcst(msg_dimqm, lefield, 1, 0)
!      call ga_brdcst(msg_dimqm, lnograd, 1, 0)
      call ga_brdcst(msg_dimqm, ldebug, 1, 0)
      call ga_brdcst(msg_dimqm, latmdip, 1, 0)
!      call ga_brdcst(msg_dimqm, ltypes, 1, 0)
      call ga_brdcst(msg_dimqm, lfirst, 1, 0)
      call ga_brdcst(msg_dimqm, lcmplx, 1, 0)
      call ga_brdcst(msg_dimqm, dimqm_noseed, 1, 0)
      call ga_brdcst(msg_dimqm, dimqm_go, 1, 0)
      mitob1= MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst(msg_dimqm, alg, mitob1, 0)
      call ga_brdcst(msg_dimqm, scrnType, mitob1, 0)
      call ga_brdcst(msg_dimqm, nDIM, mitob1, 0)
      call ga_brdcst(msg_dimqm, nOrder, mitob1, 0)
      call ga_brdcst(msg_dimqm, nDIMTypes, mitob1, 0)
      mdtob1= MA_sizeof(MT_DBL,1,MT_BYTE)
      call ga_brdcst(msg_dimqm, dimtol, mdtob1, 0)
      call ga_brdcst(msg_dimqm, scrnFactor, mdtob1, 0)
      call ga_brdcst(msg_dimqm, totcharge, mdtob1, 0)
      call ga_brdcst(msg_dimqm, lagrange, mdtob1, 0)
      call ga_brdcst(msg_dimqm, lagrange, mdtob1, 0)
      call ga_brdcst(msg_dimqm, finiteFieldX, mdtob1, 0)
      call ga_brdcst(msg_dimqm, finiteFieldY, mdtob1, 0)
      call ga_brdcst(msg_dimqm, finiteFieldZ, mdtob1, 0)
      call ga_sync()
      
      end subroutine dimqm_broadcast
