c
c  Create and save ROHF Lagrangian
c  for gradient module
c
      subroutine rohf_lagr
C$Id: rohf_lagr.F,v 1.1 1995-11-02 01:52:16 rg240 Exp $
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "cscf.fh"
#include "crohf.fh"
c
      integer g_lagr
      character*255 lagr_out
      integer nset
      logical movecs_write
      external movecs_write
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_lagr: u', 1, 1,
     $     g_lagr)) call errquit('rohf_lagr: ga_create failed', 0)
      call ga_zero(g_lagr)
      call ga_copy_patch('n', crohf_g_fcv, 1, nclosed, 1, nbf,
     $                        g_lagr, 1, nclosed, 1, nbf, )
      if (nopen.gt.0) then
        call ga_copy_patch('n',
     $           crohf_g_fpv, nclosed+1, nclosed+nopen, 1, nbf,
     $           g_lagr, nclosed+1, nclosed+nopen, 1, nbf )
        call ga_dscal_patch( g_lagr, nclosed+1, nclosed+nopen,
     $                       1, nbf, 0.5d0)
      endif
      nset = 1
c
c  Prepare to dump out
c
      call util_file_prefix('lagr', lagr_out)
      if (.not. movecs_write(lagr_out, title, trans,
     $                       nbf, nset, nbf,
     $                       dbl_mb(k_occ), nbf, dbl_mb(k_eval),
     $                       nbf, g_lagr))
     $     call errquit('rohf_lagr: movecs_write failed', 0)
      if (.not. ga_destroy(g_lagr))
     $     call errquit('rohf_lagr: destroy', 0)
      return
      end
