c
c  Create and save ROHF Lagrangian
c  for gradient module
c
      subroutine rohf_lagr
C$Id: rohf_lagr.F,v 1.5 1996-07-17 16:29:31 d3g681 Exp $
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
      integer nocc
      logical movecs_write
      external movecs_write
c
      nocc = nclosed + nopen
      if (.not. ga_create(MT_DBL, nocc, nocc, 'rohf_lagr: lagr', 1, 1,
     $     g_lagr)) call errquit('rohf_lagr: ga_create failed', 0)
      call ga_zero(g_lagr)
      call ga_copy_patch('n', crohf_g_fcv, 1, nclosed, 1, nocc,
     $                        g_lagr, 1, nclosed, 1, nocc )
      if (nopen.gt.0) then
        call ga_copy_patch('n',
     $           crohf_g_fpv, nclosed+1, nclosed+nopen, 1, nocc,
     $           g_lagr, nclosed+1, nclosed+nopen, 1, nocc )
        call ga_dscal_patch( g_lagr, nclosed+1, nclosed+nopen,
     $                       1, nocc, 0.5d0)
      endif
c$$$      CALL GA_PRINT(G_LAGR)
      call ga_screen(g_lagr,1.0d-13)
c
c  Dump out to disk
c
      nset = 1
      call util_file_name('lagr',.true.,.false.,lagr_out)
      call file_write_ga( lagr_out, g_lagr )
      if (.not. ga_destroy(g_lagr))
     $     call errquit('rohf_lagr: destroy', 0)
      return
      end
