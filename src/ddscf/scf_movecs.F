      subroutine scf_movecs_read
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
      logical movecs_read
      external movecs_read
c
      if (.not. movecs_read(movecs_in, 1,
     $     dbl_mb(k_occ), dbl_mb(k_eval),
     $     g_movecs))
     $     call errquit('scf_movecs_read failed ',0, INPUT_ERR)
c
      if (scftype .eq. 'UHF') then
         if (.not. movecs_read(movecs_in, 2,
     $        dbl_mb(k_occ+nbf), dbl_mb(k_eval+nbf),
     $        g_movecs(2))) then
            if (ga_nodeid() .eq. 0 .and.
     $           util_print('information', print_default)) then
               write(6,*) '  Duplicating RHF/ROHF vectors for UHF'
            endif
            call ga_copy(g_movecs(1), g_movecs(2))
            call dcopy(nbf,dbl_mb(k_occ),1,dbl_mb(k_occ+nbf),1)
            call dcopy(nbf,dbl_mb(k_eval),1,dbl_mb(k_eval+nbf),1)
         endif
      endif
c
      end
      subroutine scf_movecs_write(rtdb)
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "mafdecls.fh"
#include "cscf.fh"
      logical movecs_write
      external movecs_write
      integer rtdb              ! [input]
c
      integer nset, nmo2(2)
c
      nmo2(1) = nmo
      if (scftype .eq. 'UHF') then
         nset = 2
         nmo2(2) = nmo
      else
         nset = 1
      endif
c
      if (.not. movecs_write(rtdb, basis, movecs_out, 
     $     'scf', title, 
     $     nbf, nset, nmo2,
     $     dbl_mb(k_occ), nbf, dbl_mb(k_eval),
     $     nbf, g_movecs)) call errquit
     $     ('scf_movecs_write: movecs_write failed', 0, DISK_ERR)
c
      end
      subroutine scf_movecs_lock
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer l_e, k_e, l_o, k_o
      integer g_old
      integer iset, nset
      logical movecs_read
      external movecs_read
c
      if (scftype .eq. 'UHF') then
         nset = 2
      else
         nset = 1
      endif
c
c     Allocate space for the old vectors
c
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'scf_movecs_lock:old',
     $     nmo, 0, g_old)) call errquit
     $     ('scf_movecs_lock: ga_create of old failed',nmo*nbf, GA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'movecs_lock:e',l_e, k_e))
     $     call errquit('movecs_lock: insufficient memory',nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'movecs_lock:e',l_o, k_o))
     $     call errquit('movecs_lock: insufficient memory',nbf, MA_ERR)
c
      do iset = 1, nset
         if (.not. movecs_read(movecs_in, iset,
     $        dbl_mb(k_o), dbl_mb(k_e), g_old))
     $        call errquit('scf_movecs_lock: read failed',0, DISK_ERR)
c
         call movecs_lock(geom, basis, g_old, g_movecs(iset), 
     $        dbl_mb(k_eval+nbf*(iset-1)), dbl_mb(k_occ+nbf*(iset-1)))
      enddo
c
      if (.not. ma_pop_stack(l_o))call errquit('scf_movecs_lock: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_e))call errquit('scf_movecs_lock: ma?',0,
     &       MA_ERR)
      if (.not. ga_destroy(g_old))call errquit('scf_movecs_lock: ga?',0,
     &       GA_ERR)
c
      end
