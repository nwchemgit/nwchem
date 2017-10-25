c     
c     $Id: dft_scf.F 28382 2016-07-21 01:23:33Z niri $
c     
      subroutine dimqm_used(ldimqm_out)
      implicit none
      logical ldimqm_out
#include "dimqm.fh"
      ldimqm_out=ldimqm
      return
      end
      subroutine dimqm_getgo(dimqm_go_out)
      implicit none
      logical dimqm_go_out
#include "dimqm.fh"
      dimqm_go_out=dimqm_go
      return
      end
      subroutine dimqm_setgo(dimqm_go_in)
      implicit none
      logical dimqm_go_in
#include "dimqm.fh"
      dimqm_go=dimqm_go_in
      return
      end
      subroutine dimqm_setlfirst(lfirst_in)
      implicit none
      logical lfirst_in
#include "dimqm.fh"
      lfirst=lfirst_in
      return
      end
      subroutine dimqm_getenergy(edimqm_out)
      implicit none
      double precision edimqm_out
#include "dimqm.fh"
      edimqm_out=edimqm
      return
      end
      subroutine dimqm_getlrsp(lrsp_out)
      implicit none
      logical lrsp_out
#include "dimqm.fh"
      lrsp_out=lrsp
      return
      end
      subroutine dimqm_getlclfld(lclfld_out)
      implicit none
      logical lclfld_out
#include "dimqm.fh"
      lclfld_out=lclfld
      return
      end
c      subroutine dimqm_getgdipeli(g_dipel_i_out)
c      implicit none
c      integer g_dipel_i_out
c#include "dimqm.fh"
c      g_dipel_i_out=g_dipel_i
c      return
c      end
      subroutine dimqm_init(rtdb)
      implicit none
#include "dimqm.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb ! Run time database
c     Some DIM/QM global variables are needed in the rtdb to
      if (.not. rtdb_get(rtdb, 'dimqm:ldimqm', mt_log, 1, ldimqm))
     A     ldimqm=.false.
      if (.not. rtdb_get(rtdb, 'dimqm:lrsp', mt_log, 1, lrsp))
     A     lrsp=.false.
      if (.not. rtdb_get(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
     A     lfirst=.false.
      if (.not. rtdb_get(rtdb, 'dimpar:ndim', mt_int, 1, ndim))
     A     ndim=0
      if (.not. rtdb_get(rtdb,'dimpar:ndimtypes', mt_int, 1, ndimtypes))
     A     ndimtypes=0
      if (.not. rtdb_get(rtdb, 'dimpar:ndim', mt_int, 1, ndim))
     A     ndim=0
      if (.not. rtdb_get(rtdb, 'dimpar:ndim', mt_int, 1, ndim))
     A     ndim=0
      if (.not. rtdb_get(rtdb,'dimpar:norder', mt_int, 1, norder))
     A     norder=0
      if (.not. rtdb_get(rtdb,'dimqm:alg', mt_int, 1, alg))
     A     alg=0
      if (.not. rtdb_get(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
     $  dimtol=0d0
      return
      end
