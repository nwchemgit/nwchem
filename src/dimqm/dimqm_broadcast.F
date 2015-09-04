      subroutine dimqm_broadcast()
c
c     Broadcasts all common variables collected on input so all processes have the correct value
c     
      implicit none
#include "dimqm.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
c
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

      call ga_brdcst(msg_dimqm, alg, mitob(1), 0)
      call ga_brdcst(msg_dimqm, scrnType, mitob(1), 0)
      call ga_brdcst(msg_dimqm, nDIM, mitob(1), 0)
      call ga_brdcst(msg_dimqm, nOrder, mitob(1), 0)
      call ga_brdcst(msg_dimqm, nDIMTypes, mitob(1), 0)

      call ga_brdcst(msg_dimqm, dimtol, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, scrnFactor, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, totcharge, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, lagrange, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, lagrange, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, finiteFieldX, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, finiteFieldY, mdtob(1), 0)
      call ga_brdcst(msg_dimqm, finiteFieldZ, mdtob(1), 0)
      call ga_sync()
      
      end subroutine dimqm_broadcast
