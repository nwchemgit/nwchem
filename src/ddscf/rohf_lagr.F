c
c  Create and save ROHF Lagrangian
c  for gradient module
c
      subroutine rohf_lagr
C$Id$
      implicit none
#include "errquit.fh"
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
      logical movecs_write, file_write_ga
      external movecs_write, file_write_ga
c
      nocc = nclosed + nopen
*ga:1:0
      if (.not. ga_create(MT_DBL, nocc, nocc, 'rohf_lagr: lagr', 0, 0,
     $     g_lagr)) call errquit('rohf_lagr: ga_create failed', 0,
     &       GA_ERR)
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
      if (.not. file_write_ga( lagr_out, g_lagr ))
     $     call errquit('rohf: failed writing lagrangian',0, GA_ERR)
      if (.not. ga_destroy(g_lagr))
     $     call errquit('rohf_lagr: destroy', 0, GA_ERR)
      return
      end
      subroutine rohf_dumpfock(nopen,g_oexch,g_cfock)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nopen,g_oexch,g_cfock
      character*255 cfock_out,oexch_out
      logical file_write_ga
      external file_write_ga
c
      call util_file_name('cfock',.true.,.false.,cfock_out)
      call util_file_name('oexch',.true.,.false.,oexch_out)
      if (.not. file_write_ga( cfock_out, g_cfock ))
     $     call errquit('rohf: failed writing cfock',0, GA_ERR)
      if(nopen.ne.0) then
      if (.not. file_write_ga( oexch_out, g_oexch ))
     $     call errquit('rohf: failed writing oexch',0, GA_ERR)
      endif
      return
      end
      subroutine rohf_readfock(nopen,g_oexch,g_cfock)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nopen,g_oexch,g_cfock
      character*255 cfock_out,oexch_out
      logical file_read_ga
      external file_read_ga
c
      call util_file_name('cfock',.true.,.false.,cfock_out)
      call util_file_name('oexch',.true.,.false.,oexch_out)
      if (.not. file_read_ga( cfock_out, g_cfock ))
     $     call errquit('rohf: failed reading cfock',0, GA_ERR)
      if(nopen.ne.0) then
      if (.not. file_read_ga( oexch_out, g_oexch ))
     $     call errquit('rohf: failed reading oexch',0, GA_ERR)
      endif
      return
      end
