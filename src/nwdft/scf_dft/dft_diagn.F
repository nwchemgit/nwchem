      subroutine dft_diagn(levelshifting,
     ,     nmo,ncanorg, g_xij,
     .     g_fock,g_s,g_movecs,g_tmp,g_svecs,
     ,     evals,svals,noc,homo,lumo,toll_s,
     ,     oprint_intermediate_fock,oprint_time,oprint_conv_details)
c     
c     $Id: dft_diagn.F,v 1.1 2001-10-05 09:15:59 edo Exp $
c     
      implicit none
#include "cdft.fh"
#include "global.fh"
#include "mafdecls.fh"
      logical levelshifting ! [in]
      integer g_tmp   ! [in]
      integer g_fock   ! [in]
      integer g_s   ! [in]
      integer g_movecs   ! [in]
      integer g_svecs   ! [in]
      integer noc ! [in]
      double precision evals(*)   ! [out]
      double precision svals(*)   ! [in]
      double precision toll_s   ! [in]
      double precision homo,lumo ! [out]
      logical oprint_intermediate_fock,oprint_time,
     ,     oprint_conv_details ![in]
      integer nmo ! [out]
c
      integer me,j,nproc,g_fockt,g_xij,g_movecst
      integer  ga_create_atom_blocked,info
      external ga_create_atom_blocked
      integer i,ncanorg,g_vecstt
c     
c     Level shifting is implemented here (similarity 
c     transformation before standard eigensolver).  Note,
c     levelshifting is appropriate once a transformation
c     is available which makes the resulting Fock matrix 
c     diagonally dominant, e.g., in an approximate MO basis.  
c     
      me=ga_nodeid()
      nproc=ga_nnodes()

c
c     get linear dep
c
      nmo =0 
      do i=1,nbf_ao
         if(svals(i).ge.toll_s) nmo=nmo+1
      enddo
      if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'movecstransf', 
     &     nbf, 0, g_movecst)) 
     . call errquit(' ddftdiagn: ga_create failed ',1)
      if (.not. ga_create(mt_dbl, nmo, nmo, 'fock transf', 
     &     nbf, 0, g_fockt)) 
     . call errquit(' ddftdiagn: ga_create failed ',2)
c
c     get xij through canonical orthog
c
      if(ncanorg.eq.0) then
         call dft_canorg(toll_s, nmo, svals, g_svecs, g_xij)
         ncanorg=1
       endif
c
c     F' = X+ * F * X
c
      call ga_dgemm('N', 'N', nbf_ao, nmo, nbf_ao, 1d0, 
     &     g_fock, g_xij, 0d0, g_tmp)
      call ga_dgemm('T', 'N', nmo, nmo, nbf_ao, 1d0, 
     &     g_xij, g_tmp, 0d0, g_fockt)
      

      if (levelshifting)then
c     
c     Build upper diagonal matrix. in MO basis
c     
         do j = noc+1+me, nmo, nproc
            call ga_put(g_fockt, j, j, j, j, rlshift, 1)
         enddo
      else
         rlshift = 0.0
      endif
c     
c     Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
      if (oprint_intermediate_fock)then     
         call ga_print(g_fockt)
      endif
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(10) = '   ga_dgemm   '
      call tcllct(t_pgm(1,10),tp4,tp3)
      call dft_second(tp3)
#endif
      if (.not. ga_create(mt_dbl, nmo, nmo, 'vecs transf1', 
     &     nmo, 0, g_vecstt)) 
     . call errquit(' ddftdiagn: ga_create failed ',2)
      if (oprint_time)
     &     call dft_tstamp(' calling diag ')
#if defined(PARALLEL_DIAG)
      call ga_diag_std(
     .     g_fockt,  g_vecstt, evals)
#else 
      call ga_diag_std_seq(g_fockt, g_vecstt, evals)
#endif
      if (oprint_time)
     &     call dft_tstamp(' called diag ')
c     
c     Check HOMO/LUMO gap.
c     
      if(noc.gt.0) then
         homo = evals(noc)
         lumo = evals(noc+1)
      else
         homo=-100.0d0
         lumo=100.0d0
      endif
c
c     backtransform  C(k) = X * C'(k)
c
      call ga_zero(g_movecst)
      call ga_dgemm('N', 'N', nbf_ao, nmo, nmo, 1d0, 
     &     g_xij, g_vecstt, 0d0, g_movecst)
      if (.not. ga_destroy(g_vecstt)) call errquit
     &     ('dft_scf: could not destroy g_scr', 0)

c     
c     Keep orbital ordering according to principle
c     of maximum overlap with previous iteration.
c     
      if (lmaxov)then
         if(nmo.ne.nbf_ao) call errquit(
     '        ' mxovl not ready for lindep ',0)
         call dft_mxovl(ao_bas_han, evals,
     &        g_movecst, g_movecs, g_s, g_fock,
     ,        noc,oprint_conv_details,
     ,        homo,lumo)
      endif
c     
c     put new movecs
c     
      call ga_copy(g_movecst, g_movecs)
         if (.not. ga_destroy(g_fockt)) call errquit
     &        ('dft_scf: could not destroy g_scr', 0)
         if (.not. ga_destroy(g_movecst)) call errquit
     &        ('dft_scf: could not destroy g_scr', 0)
      
      return
      end
