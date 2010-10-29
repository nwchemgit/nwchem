      subroutine bas_all_ce_info(basis, cnr, bfr)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
      integer basis, cnr(2,*), bfr(2,*)
      integer ncent, icent, geom
      logical status
c
      status = bas_geom(basis, geom)
      status = status .and.  geom_ncent(geom, ncent)
      if (.not. status) call errquit('bas_all_ce_info: handles?',0,
     &       BASIS_ERR)
c
      do icent = 1, ncent
         status = bas_ce2cnr(basis, icent, cnr(1,icent), cnr(2,icent))
         status = bas_ce2bfr(basis, icent, bfr(1,icent), bfr(2,icent))
      enddo
      if (.not. status) call errquit('bas_all_ce_info: cn/bfr?',0,
     &       BASIS_ERR)
c     
      end
      subroutine bas_all_cn_info(basis, info, bfr)
      implicit none
#include "errquit.fh"
#include "bas.fh"
c
      integer basis, info(4,*), bfr(2,*)
      integer ncn, icn
      logical status
c
      if (.not. bas_numcont(basis,ncn)) 
     $      call errquit('bas_all_cn_info: handles?',0,
     &       BASIS_ERR)
      status = .true.
      do icn = 1, ncn
         status = status .and. bas_continfo(basis, icn,
     $           info(1,icn), info(2,icn), info(3,icn), info(4,icn))
         status = status .and. bas_cn2bfr(basis, icn, 
     $           bfr(1,icn), bfr(2,icn))
      enddo
      if (.not. status) call errquit('bas_all_cn_info: bad ?',0,
     &       BASIS_ERR)
c
      end
