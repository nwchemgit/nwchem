cedo#define USE_PDSYEVR 1
cedo#define USE_ELPA 1
      subroutine dft_diagn(levelshifting,
     ,     nmo,ncanorg,
     .     g_fock,g_s,g_movecs,g_tmp,g_svecs,
     ,     evals,svals,nocc,homo,lumo,toll_s,
c
c MN solvation models -->
c
c     ,     oprint_intermediate_fock,oprint_time,oprint_conv_details)
     ,     oprint_intermediate_fock,oprint_time,oprint_conv_details,
     $ istep_cosmo_vem)
c
c <-- MN solvation models
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
c MN solvation models -->
c
      integer istep_cosmo_vem
      integer g_fockt_tmp
c
      if (istep_cosmo_vem.eq.2) levelshifting=.false.
c
c <-- MN solvation models
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
c MN solvation models -->
c
      if (istep_cosmo_vem.eq.2) then
c replace off-diagonal elements of the MO Fock matrix with zeros to skip
c actual diagonalization; evals will contain the expectation values <Phi(i)|F|Phi(i)> 
c rather than correct eigenvalues
       if (.not. ga_create(mt_dbl, nmo, nmo, 'fock transf tmp',
     $     nbf, 0, g_fockt_tmp))
     . call errquit(' ddftdiagn: ga_create failed g_fockt_tmp',
     $  2, GA_ERR)
       call ga_zero(g_fockt_tmp)
       do j = 1+me, nmo, nproc
            call ga_get(g_fockt, j, j, j, j, fjj, 1)
            call ga_put(g_fockt_tmp, j, j, j, j, fjj, 1)
       enddo
       call ga_copy(g_fockt_tmp,g_fockt)
       if (.not. ga_destroy(g_fockt_tmp)) call errquit
     $        ('dft_diagn: could not destroy g_fockt_tmp', 0, GA_ERR)
      endif
c
c <-- MN solvation models
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
     &     ('dft_diagn: could not destroy g_scr', 0, GA_ERR)

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
     &        ('dft_diagn: could not destroy g_scr', 0, GA_ERR)
         if (.not. ga_destroy(g_movecst)) call errquit
     &        ('dft_diagn: could not destroy g_scr', 0, GA_ERR)
         if(ncanorg.eq.0) then
          ncanorg=1
          if (.not. ga_destroy(g_xij)) call errquit
     &       ('dft_diagn: could not destroy g_scr', 0, GA_ERR)
         endif
      
      return
      end
c
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
#ifdef USE_ELPA
      elseif(yev.eq.'e') then
         nmoout=0
         call ga_evp_real(g_s, g_svecs, evals,nb,nmoout)
      elseif(yev.eq.'l') then
         nmoout=nocc
         call ga_evp_real(g_s, g_svecs, evals,nb,nmoout)
         call dft_idtyvecs(g_svecs,nmoout+1)
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
#ifdef USE_GAPZHEEVD
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
#else
         call errquit(' dftscaleigso: not working',0,0)
#endif
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
c     
      subroutine dft_eig_so(nbf_mo, g_fockso, g_moso, evals)
c
      implicit none
c
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "zora.fh" 
#include "case.fh"
c     
c     so
      integer nbf_mo 
      integer g_fockso(2), g_moso(2)
      double precision evals(*)
c
      integer la, ia            ! complex*16 a(nbf_ao, nbf_ao)
      integer llwork 
      integer lwork, iwork      ! complex*16 work(3) 
      integer lrwork, irwork    ! double precision rwork
      integer lbuff, ibuff
      integer info
      integer i,j,i1
c
c     ma for complex diagonalizer 
      if (.not.MA_Push_Get(MT_DCpl,nbf_mo*nbf_mo,'cpl a',la,ia))
     &     call errquit('dft_scf: cannot allocate cpl a',0, MA_ERR)
      llwork = max(1, 2*nbf_mo-1)
      if (.not.MA_Push_Get(MT_DCpl,llwork,'cpl work',lwork,iwork))
     &     call errquit('dft_scf: cannot allocate cpl work',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,max(1,3*nbf_mo-2),'w.s',lrwork,
     &     irwork))
     &     call errquit('dft_scf: cannot allocate w.s',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_mo,'buff',lbuff,ibuff))
     &     call errquit('dft_scf:push_get failed', 13, MA_ERR)
c
c     call diagonalizer
      do i = 1, nbf_mo 
         do j = 1, nbf_mo
            DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))=dcmplx(0.0, 0.0)
         enddo
      enddo
      do i = 1, nbf_mo 
         call ga_get(g_fockso(1), 1,i, i,i, dbl_mb(ibuff),1)
         do j=1,i 
            DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))=
     &           dcmplx(dbl_mb(ibuff+j-1),0d0)
         enddo 
         call ga_get(g_fockso(2), 1,i, i,i, dbl_mb(ibuff),1)
         do j=1,i 
            DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))=
     $               DCpl_mb(ia+(nbf_mo)*(i-1)+(j-1))
     $           +dcmplx(0d0,dbl_mb(ibuff+j-1))
         enddo 
      enddo
c
      call zheev( 'V', 'U', nbf_mo, DCpl_mb(ia), nbf_mo, evals,
     $            DCpl_mb(iwork), LLWORK, Dbl_mb(irwork), INFO )
c
      if(ga_nodeid().eq.0) then
       do i = 1, nbf_mo
         do j = 1, nbf_mo 
            dbl_mb(ibuff+j-1)=0.0d0
            dbl_mb(ibuff+j-1)=dble(DCpl_mb(ia+nbf_mo*(i-1)+(j-1)))
         enddo 
         i1=i
         call ga_put(g_moso(1),1,nbf_mo,i1,i1,dbl_mb(ibuff),1)
         do j = 1, nbf_mo 
            dbl_mb(ibuff+j-1)=0.0d0
            dbl_mb(ibuff+j-1)=
     $             dimag(dcmplx(DCpl_mb(ia+nbf_mo*(i-1)+(j-1))))
         enddo
         i1=i 
         call ga_put(g_moso(2),1,nbf_mo,i1,i1,dbl_mb(ibuff),1)
       enddo
      endif
c
      call ga_sync()
c
c     pop memory
      if (.not.ma_pop_stack(lbuff))
     &   call errquit('dft_eig_so: cannot pop lbuff',0, MA_ERR)
      if (.not.ma_pop_stack(lrwork))
     &   call errquit('dft_eig_so: cannot pop lrwork',0, MA_ERR)
      if (.not.ma_pop_stack(lwork))
     &   call errquit('dft_eig_so: cannot pop lwork',0, MA_ERR)
      if (.not.ma_pop_stack(la))
     &   call errquit('dft_eig_so: cannot pop la',0, MA_ERR)
c     
      return
      end
