C
C $Id: hess_restart.F,v 1.4 2001-05-23 17:44:43 windus Exp $
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
#include "stdio.fh"
#include "util.fh"
c
      integer rtdb              ! [input] Run-time database handle
      integer restr             ! [input] level of calculation completed
                                ! 1 = 1 e-; 2 = 2 e-
c
      logical oprint
      integer i
c
      call ga_dgop(msg_hess_exy,dbl_mb(k_exy),n3xyz2,'+')
c
c  Print out if requested
c
      oprint = util_print('hess_follow',print_high)
      if (oprint.and.(ga_nodeid().eq.0)) then
        write(LuOut,*)
        call util_flush(LuOut)
        write(LuOut,*) 'The Hessian after contributions:',restr
        call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
        write(LuOut,*)
        write(LuOut,*) 'The Fock derivatives :'
      endif
      if (oprint) then
        do i = 1, n3xyz*ndens
          call ga_print(g_rhs(i))
        enddo
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
      call sx_write(g_rhs,n3xyz*ndens,nbf,fx_dra_handle)
c
c  Write out restr to the database
c
      if (.not. rtdb_put(rtdb, 'hessian:restart',mt_int,1,restr))
     $  call errquit('hess_restart: could not save restr to rtdb',
     $  555)
c
      return
      end
