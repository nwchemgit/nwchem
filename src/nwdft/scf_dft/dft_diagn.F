#define USE_PDSYEVR 1
      subroutine dft_diagn(levelshifting,
     ,     nmo,ncanorg,
     .     g_fock,g_s,g_movecs,g_tmp,g_svecs,
     ,     evals,svals,nocc,homo,lumo,toll_s,
     ,     oprint_intermediate_fock,oprint_time,oprint_conv_details)
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "cdft.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "dftps.fh"
#include "util.fh"
#ifdef SCALAPACK
#include "dft_scaleig.fh"
#endif
      logical levelshifting ! [in]
      integer g_tmp   ! [in]
      integer g_fock   ! [in]
      integer g_s   ! [in]
      integer g_movecs   ! [in]
      integer g_svecs   ! [in]
      integer nocc ! [in]
      double precision evals(*)   ! [out]
      double precision svals(*)   ! [in]
      double precision toll_s   ! [in]
      double precision homo,lumo ! [out]
      logical oprint_intermediate_fock,oprint_time,
     ,     oprint_conv_details ![in]
      integer nmo ! [out]
c
      integer me,j,nproc,g_fockt,g_xij,g_movecst
      integer ncanorg,g_vecstt
      double precision fjj
      integer dft_lindep
      external dft_lindep
c
c     timing
c
      double precision time1,time2
c     
c     Level shifting is implemented here (similarity 
c     transformation before standard eigensolver).  Note,
c     levelshifting is appropriate once a transformation
c     is available which makes the resulting Fock matrix 
c     diagonally dominant, e.g., in an approximate MO basis.  
c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      nmo=dft_lindep(nbf_ao,toll_s,svals,oprint_conv_details,6)
      if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'movecstransf', 
     &     nbf, 0, g_movecst)) 
     . call errquit(' ddftdiagn: ga_create failed ',1, GA_ERR)
      if (.not. ga_create(mt_dbl, nmo, nmo, 'fock transf', 
     &     nbf, 0, g_fockt)) 
     . call errquit(' ddftdiagn: ga_create failed ',2, GA_ERR)
c
c     get xij through canonical orthog
c
#ifdef SCALAPACK
      if(yev.eq.'m') ncanorg=0
#endif
      if(ncanorg.eq.0) then
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'Xij canorg', 
     &        nbf, 0, g_xij)) 
     .        call errquit(' ddftdiagn: ga_create failed ',0, GA_ERR)
         call dft_canorg(nmo, svals, g_svecs, g_xij)
      else
         g_xij=g_movecs
       endif
c
c     F' = X+ * F * X
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N', 'N', nbf_ao, nmo, nbf_ao, 1d0, 
     &     g_fock, g_xij, 0d0, g_tmp)
      call ga_dgemm('T', 'N', nmo, nmo, nbf_ao, 1d0, 
     &     g_xij, g_tmp, 0d0, g_fockt)
      if (odftps) call pstat_off(ps_dgemm)
      

      if (levelshifting.and.ncanorg.ne.0)then
c     
c     Build upper diagonal matrix. in MO basis
c     
         do j = nocc+1+me, nmo, nproc
            call ga_get(g_fockt, j, j, j, j, fjj, 1)
            fjj=fjj+rlshift
            call ga_put(g_fockt, j, j, j, j, fjj, 1)
         enddo
      else
         rlshift = 0.0
      endif
c     
c     Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
      if(util_print('dft timings', print_high))
     &      time1=util_cpusec()   ! start diagonalization time
      if (oprint_intermediate_fock)then     
         call ga_print(g_fockt)
      endif
      if (.not. ga_create(mt_dbl, nmo, nmo, 'vecs transf1', 
     &     nmo, 0, g_vecstt)) 
     . call errquit(' ddftdiagn: ga_create failed ',2, GA_ERR)
      if (oprint_time)
     &     call dft_tstamp(' calling diag ')

      if (odftps) call pstat_on(ps_diag)
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
      call dft_scaleig(g_fockt, g_vecstt, evals,
     P     nocc) 
#else
      call ga_diag_std(
     .     g_fockt,  g_vecstt, evals)
#endif
#else 
      call ga_diag_std_seq(g_fockt, g_vecstt, evals)
#endif
cedo      call ga_print(g_fockt)
cedo      call ga_print(g_vecstt)
      if (odftps) call pstat_off(ps_diag)
      if (oprint_time)
     &     call dft_tstamp(' called diag ')
      time2=util_cpusec()   ! end diagonalization time
      if(util_print('dft timings', print_high))
     &      time2=util_cpusec()   ! end diagonalization time
c
c     print diagonalization time
      if(util_print('dft timings', print_high)) then
       if (me.eq.0) then
         write(*,"(4x,'Diagonalization Time:',F13.1,'s')") time2-time1
       endif
      end if
c     
c     Check HOMO/LUMO gap.
c     
      if(nocc.gt.0) then
         homo = evals(nocc)
         lumo = evals(nocc+1)
      else
         homo=-100.0d0
         lumo=100.0d0
      endif
c
c     backtransform  C(k) = X * C'(k)
c
      call ga_zero(g_movecst)
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N', 'N', nbf_ao, nmo, nmo, 1d0, 
     &     g_xij, g_vecstt, 0d0, g_movecst)
      if (odftps) call pstat_off(ps_dgemm)
      if (.not. ga_destroy(g_vecstt)) call errquit
     &     ('dft_scf: could not destroy g_scr', 0, GA_ERR)

c     
c     Keep orbital ordering according to principle
c     of maximum overlap with previous iteration.
c     
      if (lmaxov)then
         call dft_mxovl(ao_bas_han, evals,
     &        g_movecst, g_movecs, g_s, g_fock,
     ,        nocc,oprint_conv_details,
     ,        homo,lumo)
      endif
c     
c     put new movecs
c     
      call ga_copy(g_movecst, g_movecs)
         if (.not. ga_destroy(g_fockt)) call errquit
     &        ('dft_scf: could not destroy g_scr', 0, GA_ERR)
         if (.not. ga_destroy(g_movecst)) call errquit
     &        ('dft_scf: could not destroy g_scr', 0, GA_ERR)
         if(ncanorg.eq.0) then
            ncanorg=1
            if (.not. ga_destroy(g_xij)) call errquit
     &           ('dft_scf: could not destroy g_scr', 0, GA_ERR)
         endif
      
      return
      end
      subroutine dft_remshft(levelshifting,ipol,nbf_ao,noc,
     ,     rlshift,eval)
      implicit none
      logical levelshifting
      integer ipol
      integer nbf_ao
      double precision eval(nbf_ao,*),rlshift
      integer noc(2)
c
      integer i,j
c
      if(levelshifting) then
         do i=1,ipol
            do j = noc(i)+1, nbf_ao
               eval(j,i)=eval(j,i)- rlshift
            enddo
         enddo
      endif
      return
      end
#ifdef SCALAPACK
      subroutine dft_scaleig(g_s,g_svecs,evals,
     A      nocc)
      implicit none
#include "dft_scaleig.fh"
#include "global.fh"
      integer g_s,g_svecs
      double precision evals(*)
      integer nocc,i
c
      integer nmoout
c
      if(yev.eq.'v') then
         call ga_pdsyev(g_s, g_svecs, evals,nb)
      elseif(yev.eq.'x') then
         call ga_pdsyevx(g_s, g_svecs, evals,nb)
      elseif(yev.eq.'d') then
         call ga_pdsyevd(g_s, g_svecs, evals,nb)
#ifdef USE_PDSYEVR
      elseif(yev.eq.'r'.or.yev.eq.'m') then
         nmoout=0
         if(yev.eq.'m') nmoout=nocc
         call ga_pdsyevr(g_s, g_svecs, evals,nb,nmoout)
         if(nmoout.ne.0) 
c     patching non computed eigenvectors with identity
     d        call dft_idtyvecs(g_svecs,nmoout+1)
#endif
      else
         call errquit(' wrong pdsyev value',0,0)
      endif
      return
      end
      subroutine dft_scaleig_so(g_s,g_svecs,evals,
     A      nocc)
      implicit none
#include "dft_scaleig.fh"
      integer g_s(*),g_svecs(*)
      double precision evals(*)
      integer nocc
c
      integer nmoout
c
      if(yev.eq.'d') then
         call ga_pzheevd(g_s, g_svecs, evals,nb)
#ifdef USE_PDSYEVR
      elseif(yev.eq.'r'.or.yev.eq.'m') then
         nmoout=0
         if(yev.eq.'m') nmoout=nocc
         call ga_pzheevr(g_s, g_svecs, evals,nb,nmoout)
         if(nmoout.ne.0)  then
c     patching non computed eigenvectors with identity
             call dft_idtyvecs(g_svecs(1),nmoout+1)
             call dft_idtyvecs(g_svecs(2),nmoout+1)
          endif
#endif
      else
         call errquit(' wrong pzheevs value',0,0)
      endif
      return
      end
#endif
      subroutine dft_idtyvecs(g_a,istart)
      implicit none
#include "global.fh"
      integer g_a
      integer istart
c
      integer alo(2),ahi(2)
      integer dimm,i
      integer typpe

      call ga_inquire(g_a, typpe, dimm, dimm)
      if(istart.gt.dimm) return
      alo(1)=1
      ahi(1)=dimm
      alo(2)=istart
      ahi(2)=dimm
      call nga_zero_patch(g_a, alo, ahi)
      call ga_distribution(g_a,
     D     ga_nodeid(), alo(1), ahi(1), alo(2), ahi(2))
      if (alo(1).gt.0 .and. alo(1).le.ahi(1)) then
         if(ahi(1).ge.istart) then
            do i=max(alo(1),istart),min(ahi(1),dimm)
               if(i.ge.alo(2).and.i.le.ahi(2))
     P              call ga_put(g_a,i,i,i,i,1d0,0)
            enddo
         endif
      endif
      return
      end
