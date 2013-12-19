      subroutine et_fmo_getdata(rtdb, frag, g_s, g_vecs, ma_evls,
     &   ma_occ)
      implicit none
c
c  Returns overlap matrix (g_s), MO vectors (g_vecs),
c  MO energies (ma_evls) and MO occupancies (ma_occ)
c  from the frag rtdb record and movecs file
c
c  ma_evls stack is allocated first
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "apiP.fh"
#include "stdio.fh"
#include "cscf.fh"

c
      integer rtdb
      character*(*) frag
      character*255 str
      integer g_s, g_vecs
      integer ma_evls(2), ma_occ(2) ! mhandle 1st index, index 2nd index
      double precision dum
      integer i, j, k
c
      logical et_movecs_read
      external et_movecs_read
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c
c -------------------------------
c load geometry and symmetry info
c -------------------------------
c
      if (.not. geom_create(geom, frag))
     &     call errquit('et_calc: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, frag))
     &     call errquit('et_calc: no geometry for frag', 0, RTDB_ERR)
c
c ----------------------------------------
c load the basis set and get info about it
c ----------------------------------------
c
      if (.not. bas_create(basis, 'ao basis'))
     &     call errquit('et_calc: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     &     call errquit('et_calc: no ao basis set', 0, RTDB_ERR)
c
      if (.not.bas_numbf(basis,nbf))
     &     call errquit('bas_numbf failed',0, BASIS_ERR)
c
      str = frag//'_s'
      g_s = ga_create_atom_blocked(geom, basis, str)
      call ga_zero(g_s)

      str = frag//'_mos'
      g_vecs = ga_create_atom_blocked(geom, basis, str)
      call ga_zero(g_vecs)
c
c initialize the API 
c
      call int_init(rtdb, 1, basis)
      call int_acc_std()
      call int_normalize(rtdb,basis)
c
c  -------------------------------
c  calculate the overlap integrals
c  -------------------------------
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)

      if (util_print("et",print_debug)) then
         call ga_print(g_s)
      endif
c
      str = frag//'_evls'
      if (.not.ma_push_get(mt_dbl,nbf,str,ma_evls(1),ma_evls(2)))
     & call errquit('ma_push_get failed _evls',ma_evls(2), MA_ERR)
c
      str = frag//'_occ'
      if(.not.ma_push_get(mt_dbl,nbf,str,ma_occ(1),ma_occ(2)))
     &   call errquit('ma_push_get failed _occ',ma_occ(2), MA_ERR)
c
      str = frag//'.movecs'
      if(.not.et_movecs_read(str,1,dbl_mb(ma_occ(2)),
     &   dbl_mb(ma_evls(2)),g_vecs,dum))
     & call errquit('et_movecs_read failed for frag mos',0, DISK_ERR)
c
c     destroying
c
c
      if (.not. geom_destroy(geom))
     &  call errquit('et_calc: geom destroy failed', 0, GEOM_ERR)
c
      if (.not. bas_destroy(basis))
     &  call errquit('et_calc: basis destroy failed',0, BASIS_ERR)
c
      call int_terminate()
      return
      end
c $Id$
