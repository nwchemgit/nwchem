      subroutine tddft_contract(rtdb,geom,ao_bas_han,
     1  g_td,g_pd,g_qd,nbf_ao,ipol,tol2e,tda,oskel,kfac,
     2  lhashf,triplet)
c
c $Id$
c
c Contract trial density matrices with 2-e integrals.
c (1) A trial density matrix will be decomposed into 
c a sum of symmetric and anti-symmetric contributions.
c (2) For pure DFT, only symmetric density matrix will
c suffice.  For HF or hybrid HF DFT, we must consider
c anti-symmetric part.
c (3) For (not TDA), A+B needs only symmetric and A-B
c needs only anti-symmetric, because of cancellation,
c while for (TDA), A is a sum of symmetric and anti-
c symmetric.
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "case.fh"
c
      integer rtdb             ! Runtime database
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer g_td(2)          ! GA handle for trial density matrix
      integer g_pd(2)          ! GA handle for A+B product density matrix
      integer g_qd(2)          ! GA handle for A-B product density matrix
      double precision tol2e   ! 2-e integral tolerance
      logical tda              ! True if Tamm-Dancoff approximation
      logical oskel            ! True if use symmetry in shell_fock_build
      double precision jfac    ! Prefactor for Coulomb
      double precision kfac    ! Prefactor for HF exchange
      logical lhashf           ! True if XC has HF contribution
      logical triplet          ! True if restricted triplet
c
      integer idim(3)          ! Dimension array for NGA create
      integer ichunk(3)        ! Chunk array for NGA create
      integer alo(2),ahi(2)    ! Regions for patch
      integer blo(3),bhi(3)    ! Regions for patch
      integer g_sym(2)         ! GA handles for (anti-)symmetrized density matrix
      integer g_ntrl,g_nprd    ! NGA handles for work density and product density matrices
      double precision jf(4)   ! Prefactor for HF Coulomb as input for shell_fock_build
      double precision kf(4)   ! Prefactor for HF exchange as input for shell_fock_build
      integer i
      integer ndu,nfu
      integer g_temp           ! NGA handles for temp array
      double precision scale_j, scale_k
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c -------------------------------
c Get scale_j and scale_k factors
c -------------------------------
      scale_j = 1.d0
      if (.not.rtdb_get(rtdb,'tddft:scale_j',mt_dbl,1,scale_j))
     &   scale_j = 1.d0
c
      scale_k = 1.d0
      if (.not.rtdb_get(rtdb,'tddft:scale_k',mt_dbl,1,scale_k))
     &   scale_k = 1.d0
c
c --------------------
c Allocate work arrays
c --------------------
c
      do i=1,ipol
        g_sym(i)=ga_create_atom_blocked(geom,ao_bas_han,
     1    '(Anti-)symmetrized density matrix')
      enddo
      idim(1)=ipol
      idim(2)=nbf_ao
      idim(3)=nbf_ao
      ichunk(1)=ipol
      ichunk(2)=-1
      ichunk(3)=-1
      if (.not.nga_create(mt_dbl,3,idim,'Work trial',
     1  ichunk,g_ntrl)) call errquit
     2  ('tddft_contract: failed to create g_ntrl',0, GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'Work product',
     1  ichunk,g_nprd)) call errquit
     2  ('tddft_contract: failed to create g_nprd',0, GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'Temp array',
     1  ichunk,g_temp)) call errquit
     2  ('tddft_contract: failed to create g_work',0, GA_ERR)
c
c ---------------------------
c Symmetrize density matrices
c ---------------------------
c
      do i=1,ipol
        call ga_copy(g_td(i),g_sym(i))
        call ga_symmetrize(g_sym(i))
        if (util_print('contract',print_debug))
     1    call ga_print(g_sym(i))
      enddo
      if (ipol.eq.1) call ga_dscal(g_sym(1),2.0d0)
      do i=1,ipol
        alo(1)=1
        ahi(1)=nbf_ao
        alo(2)=1
        ahi(2)=nbf_ao
        blo(1)=i
        bhi(1)=i
        blo(2)=1
        bhi(2)=nbf_ao
        blo(3)=1
        bhi(3)=nbf_ao
        call nga_copy_patch('N',g_sym(i),alo,ahi,g_ntrl,blo,bhi)
      enddo
c
c ---------------------
c Call shell_fock_build
c ---------------------
c
      call ga_zero(g_nprd)
      if ((ipol.eq.1).and.(.not.triplet)) then
        jf(1)=1.0d0*scale_j
        kf(1)=-kfac*0.5d0*scale_k
        if (.not.cam_exch) then            ! regular calculations
          call shell_fock_build(geom,ao_bas_han,0,1,
     &      jf,kf,tol2e,g_ntrl,g_nprd,.true.)
        else                                  ! attenuated calculations
          call shell_fock_build_cam(geom,ao_bas_han,0,1,
     &      jf,kf,tol2e,g_ntrl,g_nprd,.true.)
        end if  ! cam_exch
      else if ((ipol.eq.1).and.(triplet)) then
        jf(1)=0.0d0*scale_j
        kf(1)=-kfac*0.5d0*scale_k
        if (.not.cam_exch) then            ! regular calculations
          call shell_fock_build(geom,ao_bas_han,0,1,
     1      jf,kf,tol2e,g_ntrl,g_nprd,.true.)
        else                                  ! attenuated calculations
          call shell_fock_build_cam(geom,ao_bas_han,0,1,
     1      jf,kf,tol2e,g_ntrl,g_nprd,.true.)
        end if  ! cam_exch
      else if (ipol.eq.2) then
        ndu=2
        nfu=2
        jf(1)=1.0d0*scale_j
        kf(1)=kfac*scale_k
        jf(2)=1.0d0*scale_j
        kf(2)=kfac*scale_k
        call uhf_fock_setup(g_ntrl,g_nprd,jf,kf,ndu,
     1       nfu,nbf_ao,.false.)
        if (.not.cam_exch) then     ! regular calculations
          call shell_fock_build(geom,ao_bas_han,0,ndu,
     1       jf,kf,tol2e,g_ntrl,g_nprd,.true.)
        else
          call shell_fock_build_cam(geom,ao_bas_han,0,ndu,
     1       jf,kf,tol2e,g_ntrl,g_nprd,.true.)
        end if  ! cam_exch
        call uhf_fock_finish(g_ntrl,g_nprd,jf,kf,ndu,
     1       nfu,nbf_ao,.false.)
      endif  ! ipol
c
      do i=1,ipol
        alo(1)=1
        ahi(1)=nbf_ao
        alo(2)=1
        ahi(2)=nbf_ao
        blo(1)=i
        bhi(1)=i
        blo(2)=1
        bhi(2)=nbf_ao
        blo(3)=1
        bhi(3)=nbf_ao
        call nga_copy_patch('N',g_nprd,blo,bhi,g_pd(i),alo,ahi)
      enddo
      do i=1,ipol
        if (util_print('contract',print_debug))
     1    call ga_print(g_pd(i))
      enddo
c
c --------------------------------------------------
c If (HF exchange), anti-symmetrize density matrices
c --------------------------------------------------
c
      if (lhashf) then
        if (ipol.eq.1) call ga_dscal(g_sym,0.5d0)
        do i=1,ipol
          call ga_add(1.0d0,g_td(i),-1.0d0,g_sym(i),g_sym(i)) 
          if (util_print('contract',print_debug))
     1      call ga_print(g_sym(i))
        enddo
        if (ipol.eq.1) call ga_dscal(g_sym(1),2.0d0)
        do i=1,ipol
          alo(1)=1
          ahi(1)=nbf_ao
          alo(2)=1
          ahi(2)=nbf_ao
          blo(1)=i
          bhi(1)=i
          blo(2)=1
          bhi(2)=nbf_ao
          blo(3)=1
          bhi(3)=nbf_ao
          call nga_copy_patch('N',g_sym(i),alo,ahi,g_ntrl,blo,bhi)
        enddo
c
c ---------------------
c Call shell_fock_build
c ---------------------
c
        call ga_zero(g_nprd)
        if (ipol.eq.1) then
         jf(1)=0.0d0*scale_j
         kf(1)=-kfac*0.5d0*scale_k
         if (.not.cam_exch) then  ! regular calculations
           call shell_fock_build(geom,ao_bas_han,0,1,
     &       jf,kf,tol2e,g_ntrl,g_nprd,.false.)
         else                        ! attenuated calculations
           call shell_fock_build_cam(geom,ao_bas_han,0,1,
     &       jf,kf,tol2e,g_ntrl,g_nprd,.false.)
         end if  ! cam_exch
        else
          ndu=2
          nfu=2
          jf(1)=0.0d0*scale_j
          kf(1)=kfac*scale_k
          jf(2)=0.0d0*scale_j
          kf(2)=kfac*scale_k
          call uhf_fock_setup(g_ntrl,g_nprd,jf,kf,ndu,
     1        nfu,nbf_ao,.false.)
          if (.not.cam_exch) then  ! regular calculations
            call shell_fock_build(geom,ao_bas_han,0,ndu,
     1        jf,kf,tol2e,g_ntrl,g_nprd,.false.)
          else
            call shell_fock_build_cam(geom,ao_bas_han,0,ndu,
     1        jf,kf,tol2e,g_ntrl,g_nprd,.false.)
          end if ! cam_exch
          call uhf_fock_finish(g_ntrl,g_nprd,jf,kf,ndu,
     1        nfu,nbf_ao,.false.)
        endif  ! ipol
c
        do i=1,ipol
          alo(1)=1
          ahi(1)=nbf_ao
          alo(2)=1
          ahi(2)=nbf_ao
          blo(1)=i
          bhi(1)=i
          blo(2)=1
          bhi(2)=nbf_ao
          blo(3)=1
          bhi(3)=nbf_ao
          if (.not.tda) then
            call nga_copy_patch('N',g_nprd,blo,bhi,g_qd(i),alo,ahi)
          else
            call nga_add_patch(1.0d0,g_nprd,blo,bhi,
     1        1.0d0,g_pd(i),alo,ahi,g_pd(i),alo,ahi)
          endif
        enddo
        if (.not.tda) then
          do i=1,ipol
            if (util_print('contract',print_debug))
     1        call ga_print(g_qd(i))
          enddo
        else
          do i=1,ipol
            if (util_print('contract',print_debug))
     1        call ga_print(g_pd(i))
          enddo
        endif
      else if (.not.tda) then
        do i=1,ipol
          call ga_zero(g_qd(i))
        enddo
      endif
c ---------------------------------------------------------------
c If (not TDA), double g_pd to account for B matrix contribution
c to A+B and there is no contribution from anti-symmetric density
c matrix to A+B.  Similarly, A-B is twice g_qd.
c ---------------------------------------------------------------
      if (.not.tda) then
        do i=1,ipol
          call ga_dscal(g_pd(i),2.0d0)
          call ga_dscal(g_qd(i),2.0d0)
        enddo
      endif
c
c ------
c Return
c ------
c
      do i=1,ipol
        if (.not.ga_destroy(g_sym(i))) call errquit
     1    ('tddft_contract: failed to deallocate g_sym',0, GA_ERR)
      enddo
      if (.not.ga_destroy(g_ntrl)) call errquit
     1  ('tddft_contract: failed to deallocate g_ntrl',0, GA_ERR)
      if (.not.ga_destroy(g_nprd)) call errquit
     1  ('tddft_contract: failed to deallocate g_nprd',0, GA_ERR)
      if (.not.ga_destroy(g_temp)) call errquit
     1  ('tddft_contract: failed to deallocate g_temp',0, GA_ERR)
c
      return
      end
c
      subroutine tddft_nga_cont(rtdb,geom,ao_bas_han,
     1  g_ntrl,g_nprd_p,g_nprd_m,nbf_ao,ipol,tol2e,tda,oskel,kfac,
     2  lhashf,triplet,nvectors,iwhich)
c
c $Id: tddft_contract.F 26045 2014-08-26 04:09:46Z niri $
c
c Contract trial density matrices with 2-e integrals.
c Multiple trial density matrices are processed simultaneously,
c so that the number of direct integral evaluations will be 
c drastically reduced.
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "case.fh"
c
      integer rtdb             ! Runtime database
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer g_ntrl           ! NGA handle for trial density matrix
      integer g_nprd_p         ! NGA handle for A+B product density matrix
      integer g_nprd_m         ! NGA handle for A-B product density matrix
      double precision tol2e   ! 2-e integral tolerance
      logical tda              ! True if Tamm-Dancoff approximation
      logical oskel            ! True if use symmetry in shell_fock_build
      double precision kfac    ! Prefactor for HF exchange
      double precision jfac    ! Prefactor for Coulomb exchange
      logical lhashf           ! True if XC has HF contribution
      logical triplet          ! True if restricted triplet
      integer nvectors         ! Number of vectors to be processed
c
      integer idim(3)          ! Dimension array for NGA create
      integer ichunk(3)        ! Chunk array for NGA create
      integer alo(3),ahi(3)    ! Regions for patch
      integer blo(3),bhi(3)    ! Regions for patch
      integer g_nsym           ! NGA handles for (anti-)symmetrized density matrix
      integer g_work           ! NGA handles for work array
      integer l_jf,k_jf        ! Prefactor for HF Coulomb as input for shell_fock_build
      integer l_kf,k_kf        ! Prefactor for HF exchange as input for shell_fock_build
      integer n
      integer ndu,nfu
      integer g_temp           ! NGA handles for temp array
c
      double precision scale_j, scale_k
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer iwhich
c
c     NOTE: 
c     iwhich = 1 and 2 are set in tddft_lanczos.F and 1 in tddft_lanczos_tda.F
c     iwhich = 0 elsewhere
c
c -------------------------------
c Get scale_j and scale_k factors
c -------------------------------
      scale_j = 1.d0
      if (.not.rtdb_get(rtdb,'tddft:scale_j',mt_dbl,1,scale_j))
     &   scale_j = 1.d0
c
      scale_k = 1.d0
      if (.not.rtdb_get(rtdb,'tddft:scale_k',mt_dbl,1,scale_k))
     &   scale_k = 1.d0
c
c --------------------
c Allocate work arrays
c --------------------
c
      if (ipol.eq.1) then
        if (.not.ma_push_get(mt_dbl,nvectors,'J factors',
     1    l_jf,k_jf)) call errquit
     2    ('tddft_nga_cont: failed to allocate jf',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,nvectors,'K factors',
     1    l_kf,k_kf)) call errquit
     2    ('tddft_nga_cont: failed to allocate kf',0, MA_ERR)
      else
        if (.not.ma_push_get(mt_dbl,4*nvectors,'J factors',
     1    l_jf,k_jf)) call errquit
     2    ('tddft_nga_cont: failed to allocate jf',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,4*nvectors,'K factors',
     1    l_kf,k_kf)) call errquit
     2    ('tddft_nga_cont: failed to allocate kf',0, MA_ERR)
      endif
      idim(1)=ipol*nvectors
      idim(2)=nbf_ao
      idim(3)=nbf_ao
      ichunk(1)=ipol*nvectors
      ichunk(2)=-1
      ichunk(3)=-1
      if (.not.nga_create(mt_dbl,3,idim,'Work sym',ichunk,g_nsym)) 
     1 call errquit('tddft_nga_cont: failed to create g_nsym',0, GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'temp array',ichunk,g_temp)) 
     1 call errquit('tddft_nga_cont: failed to create g_temp',0, GA_ERR)
      if (lhashf.and.tda) then
       if (.not.nga_create(mt_dbl,3,idim,'Work',ichunk,g_work)) 
     1 call errquit('tddft_nga_cont: failed to create g_work',0, GA_ERR)
      endif
c
c ---------------------------
c Symmetrize density matrices
c ---------------------------
c
      alo(1)=1
      ahi(1)=ipol*nvectors
      alo(2)=1
      ahi(2)=nbf_ao
      alo(3)=1
      ahi(3)=nbf_ao
      blo(1)=1
      bhi(1)=ipol*nvectors
      blo(2)=1
      bhi(2)=nbf_ao
      blo(3)=1
      bhi(3)=nbf_ao
      call nga_copy_patch('N',g_ntrl,alo,ahi,g_nsym,blo,bhi)
      call ga_symmetrize(g_nsym)
      if (ipol.eq.1) call ga_dscal(g_nsym,2.0d0)
      if (util_print('contract',print_debug))
     1  call ga_print(g_nsym)
c
c ---------------------
c Call shell_fock_build
c ---------------------
c

      if((iwhich.eq.1).or.(iwhich.eq.0).or.tda) then

      call ga_zero(g_nprd_p)
      if ((ipol.eq.1).and.(.not.triplet)) then
       do n=1,nvectors
          dbl_mb(k_jf+n-1)=1.0d0*scale_j
          dbl_mb(k_kf+n-1)=-kfac*0.5d0*scale_k
       enddo
       if(.not.cam_exch) then ! normal calculations
c        write(6,*)'HERE'
        call shell_fock_build(geom,ao_bas_han,0,nvectors,
     1    dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_p,.true.)
       else   ! attenuated calculations
c        write(6,*)'NOT HERE'

        call shell_fock_build_cam(geom,ao_bas_han,0,nvectors,
     1    dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_p,.true.)
       end if  ! cam_exch
      else if ((ipol.eq.1).and.(triplet)) then
       do n=1,nvectors
          dbl_mb(k_jf+n-1)=0.0d0*scale_j
          dbl_mb(k_kf+n-1)=-kfac*0.5d0*scale_k
       enddo
       if (.not.cam_exch) then   ! normal calculations
        call shell_fock_build(geom,ao_bas_han,0,nvectors,
     1    dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_p,.true.)
       else   ! attenuated calculations
        call shell_fock_build_cam(geom,ao_bas_han,0,nvectors,
     1    dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_p,.true.)
       end if  ! cam_exch
      else if (ipol.eq.2) then
        do n=1,2*nvectors
          dbl_mb(k_jf+n-1)=1.0d0*scale_j
          dbl_mb(k_kf+n-1)=kfac*scale_k
        enddo
        ndu=2*nvectors
        nfu=2*nvectors
        call uhf_fock_setup(g_nsym,g_nprd_p,dbl_mb(k_jf),dbl_mb(k_kf),
     1    ndu,nfu,nbf_ao,.false.)
        if (.not.cam_exch) then   ! normal calculations
          call shell_fock_build(geom,ao_bas_han,0,ndu,
     1    dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_p,.true.)
        else
          call shell_fock_build_cam(geom,ao_bas_han,0,ndu,
     1    dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_p,.true.)
        end if ! cam_exch
        call uhf_fock_finish(g_nsym,g_nprd_p,dbl_mb(k_jf),dbl_mb(k_kf),
     1    ndu,nfu,nbf_ao,.false.)
      endif
      if (util_print('contract',print_debug))
     1  call ga_print(g_nprd_p)
       
      endif
c
c --------------------------------------------------
c If (HF exchange), anti-symmetrize density matrices
c --------------------------------------------------
c
      if((iwhich.eq.2).or.(iwhich.eq.0).or.tda) then

      if (lhashf) then
        if (ipol.eq.1) call ga_dscal(g_nsym,0.5d0)
        alo(1)=1
        ahi(1)=ipol*nvectors
        alo(2)=1
        ahi(2)=nbf_ao
        alo(3)=1
        ahi(3)=nbf_ao
        blo(1)=1
        bhi(1)=ipol*nvectors
        blo(2)=1
        bhi(2)=nbf_ao
        blo(3)=1
        bhi(3)=nbf_ao
        call nga_add_patch(1.0d0,g_ntrl,alo,ahi,-1.0d0,g_nsym,blo,bhi,
     1    g_nsym,blo,bhi) 
        if (ipol.eq.1) call ga_dscal(g_nsym,2.0d0)
        if (util_print('contract',print_debug))
     1    call ga_print(g_nsym)
c
c ---------------------
c Call shell_fock_build
c ---------------------
c
        if (tda) then
          call ga_zero(g_work)
        else
          call ga_zero(g_nprd_m)
        endif
        if (ipol.eq.1) then
          do n=1,nvectors
            dbl_mb(k_jf+n-1)=0.0d0*scale_j
            dbl_mb(k_kf+n-1)=-kfac*0.5d0*scale_k
          enddo
          if (tda) then
            if (.not.cam_exch) then ! normal calculations
             call shell_fock_build(geom,ao_bas_han,0,nvectors,
     1          dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_work,.false.)
            else  ! attenuated calculations
             call shell_fock_build_cam(geom,ao_bas_han,0,nvectors,
     1          dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_work,.false.)
            end if  ! cam_exch
            alo(1)=1
            ahi(1)=ipol*nvectors
            alo(2)=1
            ahi(2)=nbf_ao
            alo(3)=1
            ahi(3)=nbf_ao
            blo(1)=1
            bhi(1)=ipol*nvectors
            blo(2)=1
            bhi(2)=nbf_ao
            blo(3)=1
            bhi(3)=nbf_ao
            call nga_add_patch(1.0d0,g_work,alo,ahi,1.0d0,
     1        g_nprd_p,blo,bhi,g_nprd_p,blo,bhi) 
          else
            if (.not.cam_exch) then ! normal calculations
              call shell_fock_build(geom,ao_bas_han,0,nvectors,
     1          dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_m,.false.)
            else  ! attenuated calculations
              call shell_fock_build_cam(geom,ao_bas_han,0,nvectors,
     1          dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_m,.false.)
            end if  ! cam_exch
          endif ! tda
        else
          do n=1,2*nvectors
            dbl_mb(k_jf+n-1)=0.0d0*scale_j
            dbl_mb(k_kf+n-1)=kfac*scale_k
          enddo
          ndu=2*nvectors
          nfu=2*nvectors
          if (tda) then
           call uhf_fock_setup(g_nsym,g_work,
     1        dbl_mb(k_jf),dbl_mb(k_kf),ndu,nfu,nbf_ao,.false.)
           if (.not.cam_exch) then
              call shell_fock_build(geom,ao_bas_han,0,ndu,
     1        dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_work,.false.)
           else
              call shell_fock_build_cam(geom,ao_bas_han,0,ndu,
     1        dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_work,.false.)
           end if ! cam_exch
           call uhf_fock_finish(g_nsym,g_work,
     1        dbl_mb(k_jf),dbl_mb(k_kf),ndu,nfu,nbf_ao,.false.)
c
           alo(1)=1
           ahi(1)=ipol*nvectors
           alo(2)=1
           ahi(2)=nbf_ao
           alo(3)=1
           ahi(3)=nbf_ao
           blo(1)=1
           bhi(1)=ipol*nvectors
           blo(2)=1
           bhi(2)=nbf_ao
           blo(3)=1
           bhi(3)=nbf_ao
           call nga_add_patch(1.0d0,g_work,alo,ahi,1.0d0,
     1        g_nprd_p,blo,bhi,g_nprd_p,blo,bhi) 
          else
           call uhf_fock_setup(g_nsym,g_nprd_m,
     1        dbl_mb(k_jf),dbl_mb(k_kf),ndu,nfu,nbf_ao,.false.)
           if (.not.cam_exch) then
             call shell_fock_build(geom,ao_bas_han,0,ndu,
     1        dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_m,.false.)
           else
             call shell_fock_build_cam(geom,ao_bas_han,0,ndu,
     1        dbl_mb(k_jf),dbl_mb(k_kf),tol2e,g_nsym,g_nprd_m,.false.)
           end if ! cam_exch
           call uhf_fock_finish(g_nsym,g_nprd_m,
     1        dbl_mb(k_jf),dbl_mb(k_kf),ndu,nfu,nbf_ao,.false.)
          endif ! tda
        endif
        if (tda) then
          if (util_print('contract',print_debug))
     1      call ga_print(g_work)
          if (util_print('contract',print_debug))
     1      call ga_print(g_nprd_p)
        else
          if (util_print('contract',print_debug))
     1      call ga_print(g_nprd_m)
        endif
      else if (.not.tda) then
        call ga_zero(g_nprd_m)
      endif

      endif
c ---------------------------------------------------------------
c If (not TDA), double g_pd to account for B matrix contribution
c to A+B and there is no contribution from anti-symmetric density
c matrix to A+B.  Similarly, A-B is twice g_qd.
c ---------------------------------------------------------------
      if (.not.tda) then
       if((iwhich.eq.1).or.(iwhich.eq.0))
     1   call ga_dscal(g_nprd_p,2.0d0)
       if((iwhich.eq.2).or.(iwhich.eq.0))
     1   call ga_dscal(g_nprd_m,2.0d0)
      endif

c      if (.not.tda) then
c        call ga_dscal(g_nprd_p,2.0d0)
c        call ga_dscal(g_nprd_m,2.0d0)
c      endif
c
c ------
c Return
c ------
c
      if (.not.ma_pop_stack(l_kf)) call errquit
     1  ('tddft_nga_cont: failed to deallocate kf',0, MA_ERR)
      if (.not.ma_pop_stack(l_jf)) call errquit
     1  ('tddft_nga_cont: failed to deallocate jf',0, MA_ERR)
      if (.not.ga_destroy(g_nsym)) call errquit
     1  ('tddft_nga_cont: failed to deallocate g_nsym',0, GA_ERR)
      if (.not.ga_destroy(g_temp)) call errquit
     1  ('tddft_nga_cont: failed to deallocate g_temp',0, GA_ERR)
      if (lhashf.and.tda) then
        if (.not.ga_destroy(g_work)) call errquit
     1    ('tddft_nga_cont: failed to deallocate g_work',0, GA_ERR)
      endif
      return
      end
