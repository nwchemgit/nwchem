C
C $Id: hess_restart.F,v 1.1 2000-08-10 03:00:06 windus Exp $
C
      subroutine hess_restart(rtdb, restr)
c
c     Write information to the database for restart and also put the
c     hessian and fock derivatives out to disk.  Before putting the
c     hessian out to disk, we will do a global summation to make sure
c     that we get all of the contributions out to disk.  Dividing by
c     the number of processors will prepare the hessian for the next
c     contribution.
c
      implicit none
c
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer rtdb              ! [input] Run-time database handle
      integer restr             ! [input] level of calculation completed
                                ! 1 = 1 e-; 2 = 2 e-
c
      call ga_dgop(msg_hess_exy,dbl_mb(k_exy),n3xyz2,'+')
      if (ga_nodeid().eq.0) then
        write(6,*) 'The Hessian after contributions: ',restr
        call flush(6)
        call hnd_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
      endif                                                                     
c
c  Write hessian out to disk
c
      call hess_wrt()
c
c  Scale for next step
c
      call dscal(n3xyz2,1.0d+00/ga_nnodes(),dbl_mb(k_exy),1)
c
c  Write out fock derivatives to disk
c
      call sx_write(g_rhs,n3xyz,nbf,fx_dra_handle)
c
c  Write out restr to the database
c
      if (.not. rtdb_put(rtdb, 'hessian:restart',mt_int,1,restr))
     $  call errquit('hess_restart: could not save restr to rtdb',
     $  555)
c
      return
      end
