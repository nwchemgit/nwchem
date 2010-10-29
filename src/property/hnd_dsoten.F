      subroutine hnd_dsoten(rtdb,geom,basis,g_dens,ndens,pairlist,
     &                      translate,ipairs,unique,i_pert,i_resp,
     &                      tensor,coords,nbf,ndso)
c $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "apiP.fh"
c
      integer rtdb                             ! [input] rtdb handle
      integer basis                            ! [input] basis handle
      integer geom                             ! [input] geometry handle
      integer ndens                            ! [input] # of density matrices
      integer g_dens(ndens)                    ! [input] density
      integer ipairs                           ! [input] number of spin-spin pairs
      integer pairlist(2*ipairs)               ! [input] list of the pairs
      integer translate(2*ipairs)              ! [input] translation to unique list
      integer i_pert, i_resp                   ! [input] # of unique responding and perturbing atoms
      integer unique(i_pert+i_resp)            ! [input] list of unique atoms
      double precision tensor(3,3,5,ipairs)    ! [output] spin-spin tensor, one for each term
      double precision coords(3,i_pert+i_resp) ! [input] coordinates of unique atoms
      integer nbf                              ! [input] basis function info
      double precision ndso                    ! [input] dso prefactor
c
      integer ix, iy
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer g_fca, l_pr, k_pr
      integer i, j, iat, jat, k
      double precision tol2e
c
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
c
      logical     oskel
      double precision valuea, valueb
      data tol2e   /1.0d-16/
c
      oskel = .false.
c
c     Integral initialization
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 9*ipairs
      if (.not.nga_create(MT_DBL,3,ahi,'fca matrix',alo,g_fca)) call
     &    errquit('hnd_dso: nga_create failed g_fca',0,GA_ERR)
      call ga_zero(g_fca)
c
c     Generate coordinate list for DSO integrals
c
      if (.not. ma_push_get(mt_dbl,ipairs*6,'coords',l_pr,k_pr))
     &   call errquit('hnd_dso: ma_push_get failed k_pr',0,MA_ERR)
      do i = 1, ipairs
       j = (i-1)*3
       k = j+3*ipairs
       dbl_mb(k_pr+j  ) = coords(1,translate(i       ))
       dbl_mb(k_pr+j+1) = coords(2,translate(i       ))
       dbl_mb(k_pr+j+2) = coords(3,translate(i       ))
       dbl_mb(k_pr+k  ) = coords(1,translate(i+ipairs)+i_pert)
       dbl_mb(k_pr+k+1) = coords(2,translate(i+ipairs)+i_pert)
       dbl_mb(k_pr+k+2) = coords(3,translate(i+ipairs)+i_pert)
      enddo
      call int_giao_1ega(basis,basis,g_fca,'dso',dbl_mb(k_pr),ipairs,
     &                   oskel)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 0
      bhi(3) = 0
      do i = 1, ipairs
         iat = pairlist(i)
         jat = pairlist(i+ipairs)
         do ix = 1, 3
            do iy = 1, 3
               blo(3) = blo(3)+1
               bhi(3) = bhi(3)+1
          valuea=nga_ddot_patch(g_dens(1),'n',alo,ahi,g_fca,'n',blo,bhi)
          valueb=nga_ddot_patch(g_dens(2),'n',alo,ahi,g_fca,'n',blo,bhi)
            tensor(ix,iy,5,i)=(valuea+valueb)*ndso
c           if(ga_nodeid().eq.0) write(*,'(A,5i3,3F12.4)')
c    &                 'DSO iat jat ix iy', i,
c    &                pairlist(i),pairlist(i+ipairs),ix,iy,
c    &                tensor(ix,iy,5,i),valuea,valueb
            enddo
         enddo
c        if(ga_nodeid().eq.0)print*,''
      enddo
c
      call ga_sync()
c
c     Clean up all remaining memory
c
      if (.not.ga_destroy(g_fca)) call 
     &    errquit('hnd_dso: ga_destroy failed g_fca',0,GA_ERR)
c
      if (.not.ma_pop_stack(l_pr)) call
     &    errquit('hnd_dso: ma_pop_stack failed l_pr',0,MA_ERR)
c
      call schwarz_tidy()
      call int_terminate()
c
      return
      end
