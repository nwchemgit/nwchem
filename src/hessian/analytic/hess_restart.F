C
C $Id$
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
#include "errquit.fh"
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
      integer ind, idens, irhs
      integer ilo(3), ihi(3)
      double precision dnrm
c
      logical oprint, olprint
      logical ocphfprint, olcphfprint
c
      call ga_dgop(msg_hess_exy,dbl_mb(k_exy),n3xyz2,'+')
c
c     Print out if requested
c
      oprint = util_print('hess_follow',print_high)
      olprint = (oprint.and.(ga_nodeid().eq.0))
      ocphfprint = util_print('hess_cont',print_debug)
      olcphfprint = (oprint.and.(ga_nodeid().eq.0))
c
      if (olprint) then
        write(LuOut,*)
        write(LuOut,*) 'The Hessian after contributions:',restr
        call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
      endif
      if (ocphfprint) then
        if (olcphfprint) then
          write(LuOut,*)
          write(LuOut,*) 'The Fock derivatives :'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do irhs = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_rhs,ilo,ihi,1)
          enddo
        enddo
      endif
      if (oprint) then
        if (olprint) then
          write(LuOut,*)
          write(LuOut,*) 'The Fock derivatives :'
        endif
        ind = 0
        do idens = 1, ndens
          do irhs = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_normf_patch(g_rhs,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": dFock =",i2,i4,f24.8)')
     +        'hess_restart',
     +        idens,irhs,dnrm
            endif
          enddo
        enddo
      endif
c
c  Write hessian out to disk
c
      if(ga_nodeid().eq.0) then
         call hess_wrt()
      endif
c
c  Scale for next step
c
      call dscal(n3xyz2,1.0d+00/ga_nnodes(),dbl_mb(k_exy),1)
c
c  Write out fock derivatives to disk
c
c      call sx_write(g_rhs,n3xyz*ndens,nbf,fx_dra_handle)
c
c  Write out restr to the database
c
c      if (.not. rtdb_put(rtdb, 'hess:restart',mt_int,1,restr))
c     $  call errquit('hess_restart: could not save restr to rtdb',
c     $  555, RTDB_ERR)
c
      return
      end
C
      subroutine hess_getrestart(irestart)
c
c     Get hessian and fock derivatives from disk.  After getting the
c     hessian from disk, we will divide the hessian by the number of
c     processors to prepare for the next contribution.
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "tcgmsg.fh"
#include "hess_info.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer irestart              ! [input] restart index
c
      integer i, ind, idens, irhs
      integer ilo(3), ihi(3)
      double precision dnrm
      logical oprint, olprint, ocphfprint, olcphfprint
c
      oprint = util_print('hess_follow',print_high)
      olprint = (oprint.and.(ga_nodeid().eq.0))
      ocphfprint = util_print('hess_cont',print_debug)
      olcphfprint = (oprint.and.(ga_nodeid().eq.0))
c
      if (ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*) 'Getting restart information'
        write(LuOut,*)
      endif
c
c  Read hessian and dipole derivative from disk and broadcast
c
      if(ga_nodeid().eq.0) then
         call hess_read(ncent, dbl_mb(k_exy))
         call ddip_read(ncent, dbl_mb(k_dder))
      endif
      call ga_brdcst(msg_hess_exy+MSGDBL, dbl_mb(k_exy), 
     &               mdtob(n3xyz2), 0)
      call ga_brdcst(msg_hess_ddip+MSGDBL, dbl_mb(k_dder),
     &               mdtob(3*n3xyz),0)
c
c  Read Fock derivatives from disk
c
      call sx_read(g_rhs,n3xyz*ndens,nbf,fx_dra_handle)
c
c  Print out if requested
c
      if (olprint) then
        write(LuOut,*)
        write(LuOut,*) 'The Hessian at restart:'
        call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
      endif
      if (ocphfprint) then
        if (olcphfprint) then
          write(LuOut,*)
          write(LuOut,*) 'The Fock derivatives :'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do irhs = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_rhs,ilo,ihi,1)
          enddo
        enddo
      endif
      if (oprint) then
        if (olprint) then
          write(LuOut,*)
          write(LuOut,*) 'The Fock derivatives :'
        endif
        ind = 0
        do idens = 1, ndens
          do irhs = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_normf_patch(g_rhs,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": dFock =",i2,i4,f24.8)')
     +        'hess_restart',
     +        idens,irhs,dnrm
            endif
          enddo
        enddo
      endif
c
c  Scale for next step
c
      call dscal(n3xyz2,1.0d+00/ga_nnodes(),dbl_mb(k_exy),1)
      call dscal(n3xyz*3,1.0d+00/ga_nnodes(),dbl_mb(k_dder),1)
c
c  Get rid of density allocation if we are done with the 2e- contribution
c
cng      if (irestart.ge.2) then
        do i = 1, ndens
          if (.not.ga_destroy(g_dens(i)))
     &     call errquit('hess_getrestart: failed to delete density',i,
     &       GA_ERR)
        enddo
        if (.not.ga_destroy(g_rhs_xc))
     &    call errquit('hess_getrestart: failed to delete g_rhs_xc',
     &                  555, GA_ERR)
         if (.not.ga_destroy(g_wdens))
     &   call errquit('hess_getrestart: failed to delete wdensity',
     &                 555, GA_ERR)
        call schwarz_tidy()
        call intdd_terminate()
cng      endif
c
      return
      end
