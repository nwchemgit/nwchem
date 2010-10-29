      subroutine tddft_finalize(geom,ao_bas_han,ipol,g_movecs,
     1  l_evl,k_evl,l_occ,k_occ,l_irs,k_irs)
c
c $Id$
c
c Destroy handles for clean exit of TDDFT module.
c Written by So Hirata, Feb 2002.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      integer g_movecs(2)      ! GA handle for MO vectors
      integer l_evl,k_evl      ! MA handles for eigenvalues (AB)
      integer l_occ,k_occ      ! MA handles for occupancies (AB)
      integer l_irs(2),k_irs(2)! MA handles for occupancies (AB)
      integer ipol             ! =1 (RDFT); =2 (UDFT)
      integer i
c --------
c Geometry
c --------
      if (.not.geom_destroy(geom))
     1  call errquit('tddft_finalize: geom_destroy failed',0, GEOM_ERR)
c ---------
c Basis set
c ---------
      if (.not.bas_destroy(ao_bas_han))
     1  call errquit('tddft_finalize: bas_destroy failed',0, BASIS_ERR)
c ----------
c MO vectors
c ----------
      do i=1,ipol
        if (.not.ga_destroy(g_movecs(i))) call errquit
     1    ('tddft_finalize: failed to destroy g_movecs',0, GA_ERR)
      enddo
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_irs(i))) call errquit
     1    ('tddft_finalize: failed to deallocate irs',0, MA_ERR)
      enddo
      if (.not.ma_pop_stack(l_occ)) call errquit
     1  ('tddft_finalize: failed to deallocate occ',0, MA_ERR)
      if (.not.ma_pop_stack(l_evl)) call errquit
     1  ('tddft_finalize: failed to deallocate evl',0, MA_ERR)
c
c Return
c
      return
      end
