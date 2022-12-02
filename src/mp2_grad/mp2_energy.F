      subroutine mp2_energy(rtdb,geom,basis,movecs,task)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "stdio.fh"
#include "sym.fh"
#include "cmp2ps.fh"
#include "msgids.fh"
      character*8 task
      integer num_oa_sym(0:7),num_va_sym(0:7)    
      integer num_ob_sym(0:7),num_vb_sym(0:7)
      integer sym_lo_oa(0:7)    ! First MO in each sym
      integer sym_hi_oa(0:7)    ! Last MO in each sym
      integer sym_lo_va(0:7)
      integer sym_hi_va(0:7)
      integer sym_lo_ob(0:7)    ! First MO in each sym
      integer sym_hi_ob(0:7)    ! Last MO in each sym
      integer sym_lo_vb(0:7)
      integer sym_hi_vb(0:7)
c
      integer num_foa_sym(0:7),num_fva_sym(0:7) ! Ditto for frozen
      integer num_fob_sym(0:7),num_fvb_sym(0:7) 
      integer sym_lo_foa(0:7)
      integer sym_hi_foa(0:7)
      integer sym_lo_fva(0:7)
      integer sym_hi_fva(0:7)
      integer sym_lo_fob(0:7)
      integer sym_hi_fob(0:7)
      integer sym_lo_fvb(0:7)
      integer sym_hi_fvb(0:7)
      integer nmixed
      integer rtdb,basis,geom
      integer g_vecs_a,g_vecs_b
      integer nocc,noc(2)
      character*4 scftype
      character*20 scftypevecs
      integer i,nbf_vecs,nsets,nmos(2),nbf, nmo
      integer l_eval_a,k_eval_a,k_eval_b
      integer noa_frozen,nva_frozen,nob_frozen,nvb_frozen
      integer noa_lo,noa_hi,nob_lo,nob_hi,nva_lo,nva_hi,nvb_lo,nvb_hi
      integer noa_active,nob_active,nva_active,nvb_active
      integer k_irs_a,k_irs_b,l_irs_a,l_irs_b
      integer k_occ_a,k_occ_b,l_occ_a
      integer pass
      integer noa,nob,nva,nvb
      double precision tol2e,tol2e_fock,backtol,escf,e2,e_total,norm
      double precision esing, etrip
      integer ga_create_atom_blocked
      logical otrace, omatrix, oskel, osym
      external ga_create_atom_blocked
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
      character*(nw_max_path_len) movecs
      character*255 title_vecs,basis_vecs
      integer maxireps,maxops
      parameter(maxireps=20,maxops=120)
      integer k_work,l_work
      integer k_map_a,l_map_a
      integer k_map_b,l_map_b
      intrinsic nint
      logical status,oblk
      integer g_p_ab_a,g_p_ab_b,g_p_ij_a,g_p_ij_b,g_lai_a,g_lai_b
      integer g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b,g_w_ai_a,g_w_ai_b
      integer g_liy_a, g_liy_b
c     
      integer g_rhs, g_vecs_swap(2)
      integer grad_length
      integer g_p_tot_a,g_p_tot_b
      integer g_w_tot_a,g_w_tot_b
      logical transform
      logical cphf
      external cphf
      logical file_write_ga,file_read_ga
      external file_write_ga,file_read_ga

      character*(nw_max_path_len) cphf_sol
      integer natoms, l_grad, k_grad, g_tmpa, g_tmpb
c
      logical dftmp2            ! use DFT reference orbitals and energy
      double precision mp2fac   ! coefficient on MP2 contribution to total energy
      double precision ess,  fss, eos, fos ! Spin Component Scaled (SCS) MP2
      logical scs                          ! SCS gradient
      logical doprint
c
c     Check for dftmp2
c
      if (.not. rtdb_get(rtdb,'dft:dftmp2',mt_log,1,dftmp2))
     $     dftmp2 = .false.
c     
      if (.not.rtdb_get(rtdb, 'mp2:skeleton',MT_LOG, 1, oskel))
     $     oskel = sym_number_ops(geom) .gt. 0
      if (.not.rtdb_get(rtdb, 'mp2:osym',MT_LOG, 1, osym))
     $     osym = oskel .and. sym_abelian_group(geom)
c
      otrace = util_print('matrix trace',print_debug)
      omatrix = util_print('matrix debug', print_debug)
c     
c     get some space for the derivatives
c     
      status = geom_ncent ( geom, natoms )
      if(.not.ma_push_get(mt_dbl, 3*natoms,'mp2 derivatives',
     $     l_grad,k_grad))      
     $     call errquit('mp2_energy ma problem: derivatives',1, MA_ERR)
      call dfill ( 3*natoms, 0.D0, dbl_mb(k_grad), 1 )
c     
*     only blocking currently works
*      if(.not.rtdb_get(rtdb,'fourindex:aoblock',mt_log,oblk,1))
*     $     oblk=.false.
      oblk = .true.
c     
c     Get no. of MOs ... nmo != nbf if linear dependence present
c     
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('mp2_grad: could not get nbf',0, BASIS_ERR)
      if(.not.movecs_read_header(movecs,title_vecs,basis_vecs,
     $     scftypevecs, nbf_vecs,nsets,nmos,2))
     $     call errquit('mp2_grad: problem with movecs_read_header',1,
     &       DISK_ERR)
      if(nbf.ne.nbf_vecs)call errquit('mp2: MO vectors mismatch?',0,
     &       INPUT_ERR)
      nmo = nmos(1)
c     
C      if (rtdb_get(rtdb,'mp2:frozen core:freeze by atoms',mt_log, 1, 
C     $     status)) then
C         if (.not. geom_num_core(geom,noa_frozen)) 
C     $        call errquit('semi_direct_mp2: geom_num_core?',0)
C      else if (rtdb_get(rtdb, 'mp2:frozen core', MT_INT, 1,
C     $        noa_frozen)) then
C      else
C         noa_frozen = 0
C      endif

       if(.not. geom_freeze(rtdb,geom,'mp2',noa_frozen)) then
            call errquit('mp2_grad: unable to get frozen core',0,
     &       GEOM_ERR)
       endif
c
      if(.not.rtdb_put(rtdb,'mp2:number frozen core',MT_INT,1,
     $      noa_frozen))
     $     call errquit('mp2_grad: rtdb failure for frozen core',555,
     &       RTDB_ERR)
c
      if (.not.rtdb_get(rtdb, 'mp2:frozen virtual', MT_INT, 1,
     $     nva_frozen)) nva_frozen = 0
      if ((task.eq.'gradient') .and. (nva_frozen.gt.0)) call errquit
     $     ('mp2_gradient: unable to freeze virtuals for gradient',0,
     &       RTDB_ERR)
c
      nob_frozen = noa_frozen
      nvb_frozen = nva_frozen
c     
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
     $     call errquit('mp2:scftype absent',1, RTDB_ERR)
      if(scftype.eq.'UHF')then
        if (dftmp2) then
          if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc)) 
     1        call errquit('mp2: rtdb_get of noc failed',0,RTDB_ERR)
          noa=noc(1)
          nob=noc(2)
        else
          if(.not.rtdb_get(rtdb,'scf:nalpha',mt_int,1,noa))
     $        call errquit('no alphas',1, RTDB_ERR)
          if(.not.rtdb_get(rtdb,'scf:nbeta',mt_int,1,nob))
     $        call errquit('no betas',1, RTDB_ERR)
        endif
      else if(scftype.eq.'RHF')then
         if(.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nocc))
     $        call errquit('no occ',1, RTDB_ERR)
         noa=nocc
         nob=nocc
      else
         if(ga_nodeid().eq.0)  write(6,*)  ' SCFTYPE = ',SCFTYPE
         call errquit('mp2_grad: unknown SCFTYPE ', 0, INPUT_ERR)
      end if
      nva=nmo-noa
      nvb=nmo-nob
c     
c     Get the MO vectors, swap and then symmetry block them
c     
*ga:1:0
      if(.not.ga_create(mt_dbl,nbf,nmo,'alpha_vecs',nbf,0,
     $     g_vecs_a))call errquit('mp2_energy: ga_create failed',0,
     &       GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nbf,nmo,'beta_vecs',nbf,0,g_vecs_b))
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)
c     
      if(.not.ma_push_get(mt_dbl,nbf*2,'occ alpha',
     $     l_occ_a,k_occ_a))      
     $     call errquit('mp2_energy ma problem with alpha occ',1,
     &       MA_ERR)
      k_occ_b = k_occ_a + nbf   ! So that are contiguous for swap. Ugh.
c     
      if(.not.ma_push_get(mt_dbl,2*nbf,'evals',
     $     l_eval_a,k_eval_a))
     $     call errquit('problem with ma for evals',1, MA_ERR)
      k_eval_b = k_eval_a + nbf ! So that are contiguous for swap. Ugh.
c     
      if(.not.movecs_read(movecs,1,dbl_mb(k_occ_a),
     $     dbl_mb(k_eval_a),g_vecs_a))
     $     call errquit('could not read alpha vecs',1, DISK_ERR)

      if(scftype.eq.'UHF')then
         if(.not.movecs_read(movecs,2,dbl_mb(k_occ_b),
     $        dbl_mb(k_eval_b),g_vecs_b))
     $        call errquit('could not read beta vecs',1, DISK_ERR)
      endif
c
      g_vecs_swap(1) = g_vecs_a
      g_vecs_swap(2) = g_vecs_b
      call movecs_swap(rtdb, 'mp2', scftype, g_vecs_swap, 
     $     dbl_mb(k_occ_a), dbl_mb(k_eval_a))
c
      if (scftype .ne. 'UHF') then
         call ga_copy(g_vecs_a,g_vecs_b) ! RHF just duplicate vectors
         call dcopy(nmo, dbl_mb(k_eval_a), 1, dbl_mb(k_eval_b), 1)
      end if
c     
      if(.not.ma_push_get(mt_int,nmo,'irs_a',l_irs_a,k_irs_a))
     $     call errquit('scf_movecs_sym_adapt: memory?',nmo, MA_ERR)
c     
      if(.not.ma_push_get(mt_int,nmo,'irs_b',l_irs_b,k_irs_b))
     $     call errquit('scf_movecs_sym_adapt: memory?',nmo, MA_ERR)
c     
      call int_init (rtdb, 1, basis) ! Needed for overlap
c     
      if (osym) then
         call sym_movecs_adapt(basis,1d-8,g_vecs_a,int_mb(k_irs_a),
     $        nmixed)
         call sym_movecs_adapt(basis,1d-8,g_vecs_b,int_mb(k_irs_b),
     $        nmixed)      
      else
         do i = 0, nmo-1
            int_mb(k_irs_a + i) = 1
            int_mb(k_irs_b + i) = 1
         enddo
      endif
c
      if (omatrix) then
         call ga_print(g_vecs_a)
         call ga_print(g_vecs_b)
      endif
c     
      call int_terminate()
c     
c     Change irreps to start counting at 0
c     
      do i=1,nmo
         int_mb(k_irs_a+i-1)=int_mb(k_irs_a+i-1)-1
         int_mb(k_irs_b+i-1)=int_mb(k_irs_b+i-1)-1
      end do
c     
      noa_lo=noa_frozen+1
      noa_hi=noa
      noa_active=noa_hi-noa_lo+1
      nva_lo=noa+1
      nva_hi=nmo-nva_frozen
      nva_active=nva_hi-nva_lo+1
c     
      nob_lo=nob_frozen+1
      nob_hi=nob
      nob_active=nob_hi-nob_lo+1
      nvb_lo=nob+1
      nvb_hi=nmo-nvb_frozen
      nvb_active=nvb_hi-nvb_lo+1
c     
c     map_a/b hold reordering vector from MO file to symmetry blocked
c     
      if(.not.ma_push_get(mt_int,nmo,'map_a',l_map_a,k_map_a))
     $     call errquit('could not alloc map_a',1, MA_ERR)
      if(.not.ma_push_get(mt_int,nmo,'map_b',l_map_b,k_map_b))
     $     call errquit('could not alloc map_b',1, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbf,'work',l_work,k_work)) ! Holds an MO
     $     call errquit('could not alloc work',1, MA_ERR)
c
      if (noa_frozen .gt. 0)
     $     call moints_vecs_sym_sort(g_vecs_a,nbf,1,noa_lo-1,
     $     int_mb(k_irs_a),int_mb(k_map_a),dbl_mb(k_work),num_foa_sym,
     $     sym_lo_foa,sym_hi_foa)
      if (nob_frozen .gt. 0) 
     $     call moints_vecs_sym_sort(g_vecs_b,nbf,1,nob_lo-1,
     $     int_mb(k_irs_b),int_mb(k_map_b),dbl_mb(k_work),num_fob_sym,
     $     sym_lo_fob,sym_hi_fob)
c
      call moints_vecs_sym_sort(g_vecs_a,nbf,noa_lo,noa_hi,
     $     int_mb(k_irs_a),int_mb(k_map_a),dbl_mb(k_work),num_oa_sym,
     $     sym_lo_oa,sym_hi_oa)
      call moints_vecs_sym_sort(g_vecs_b,nbf,nob_lo,nob_hi,
     $     int_mb(k_irs_b),int_mb(k_map_b),dbl_mb(k_work),num_ob_sym,
     $     sym_lo_ob,sym_hi_ob)
c     
      call moints_vecs_sym_sort(g_vecs_a,nbf,nva_lo,nva_hi,
     $     int_mb(k_irs_a),int_mb(k_map_a),dbl_mb(k_work),num_va_sym,
     $     sym_lo_va,sym_hi_va)
      call moints_vecs_sym_sort(g_vecs_b,nbf,nvb_lo,nvb_hi,
     $     int_mb(k_irs_b),int_mb(k_map_b),dbl_mb(k_work),num_vb_sym,
     $     sym_lo_vb,sym_hi_vb)
c
      if (nva_frozen .gt. 0)
     $     call moints_vecs_sym_sort(g_vecs_a,nbf,nva_hi+1,nmo,
     $     int_mb(k_irs_a),int_mb(k_map_a),dbl_mb(k_work),num_fva_sym,
     $     sym_lo_fva,sym_hi_fva)
      if (nvb_frozen .gt. 0)
     $     call moints_vecs_sym_sort(g_vecs_b,nbf,nvb_hi+1,nmo,
     $     int_mb(k_irs_b),int_mb(k_map_b),dbl_mb(k_work),num_fvb_sym,
     $     sym_lo_fvb,sym_hi_fvb)
c
      if (omatrix) then
         call ga_print(g_vecs_a)
         call ga_print(g_vecs_b)
      endif
c     
c     Also reorder the evals
c     
      call dcopy(nbf,dbl_mb(k_eval_a),1,dbl_mb(k_work),1)
      do i=1,nva_hi
         dbl_mb(k_eval_a-1+int_mb(k_map_a+i-1))=dbl_mb(k_work+i-1)
      end do
c     
      call dcopy(nmo,dbl_mb(k_eval_b),1,dbl_mb(k_work),1)
      do i=1,nvb_hi
         dbl_mb(k_eval_b-1+int_mb(k_map_b+i-1))=dbl_mb(k_work+i-1)
      end do
c
*      write(6,*) ' reordering vector ', (int_mb(k_map_a+i-1),i=1,nmo)
*      write(6,*) ' reordering vector ', (int_mb(k_map_b+i-1),i=1,nmo)
*      call output(dbl_mb(k_eval_a),1,nmo,1,1,nmo,1,1)
*      call output(dbl_mb(k_eval_b),1,nmo,1,1,nmo,1,1)
c     
      if (.not. ma_pop_stack(l_work)) call errquit
     $     ('mp2_grad: ma failed popping workspace(1)', 0, MA_ERR)
c     
      if (.not. rtdb_get(rtdb, 'mp2:aotol2e', mt_dbl, 1, tol2e))
     $     tol2e = 1d-9  ! For screening energy integrals
      if (.not. rtdb_get(rtdb, 'mp2:aotol2e fock', mt_dbl, 1, 
     $     tol2e_fock)) tol2e_fock = 1d-9 ! For overall accuracy of Fock builds
      if (.not. rtdb_get(rtdb, 'mp2:backtol', mt_dbl, 1, backtol))
     $     backtol = 1d-9  ! For screening derivative integrals
c
      doprint=util_print('information',print_low).and.ga_nodeid().eq.0
      if(doprint) then
        write(6,123) nbf, nmo, noa_frozen, nva_frozen, 
     $        noa_active, nob_active, nva_active, nvb_active,
     $        osym, oskel, tol2e, tol2e_fock, backtol
 123     format(
     $        '  Basis functions       = ', i6/
     $        '  Molecular orbitals    = ', i6/
     $        '  Frozen core           = ', i6/
     $        '  Frozen virtuals       = ', i6/
     $        '  Active alpha occupied = ', i6/
     $        '  Active beta occupied  = ', i6/
     $        '  Active alpha virtual  = ', i6/
     $        '  Active beta virtual   = ', i6/
     $        '  Use MO symmetry       = ', l1/
     $        '  Use skeleton AO sym   = ', l1//
     $        '  AO/Fock/Back tols     = ', 1p,3d9.1)
         call util_flush(6)
      endif
c
c     Print number of alpha and beta electrons to the ecce.out file
c
      call ecce_print1('number of frozen orbitals',mt_int,
     $                  noa_frozen,1)
      call ecce_print1('number of active alpha electrons',mt_int,
     $                  noa_active,1)
      call ecce_print1('number of active beta electrons',mt_int,
     $                  nob_active,1)
c
c     RJH DEBUG
c
*      call testtwopdm(rtdb, g_vecs_a, noa, nmo, nbf)
c     
      e2=0.0d0
      esing = 0.0d0
      etrip = 0.0d0

c     flag and scaling factors for SCS-MP2

      if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $  fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
      if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $  fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
      if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $  scs = .false.       ! flag for SCS gradient
      if (doprint.and.scs) then
         write(6,'(''  SCS-MP2               = '', l1)') scs
         call util_flush(6)
      endif
c
c     this protect against cases where a single electron is present
c
      if(noa_active+nob_active.lt.2) goto 1984
c     
c...  A note on notation: arrays for four index quantities of pure spin
c...  denoted with a _a on the end. Mixed spin case with a _b. On first
c...  call to grad_uhf _a arrays are all alpha,_b are <alpha beta||alpha beta>.
c...  On second call to grad_uhf vice versa is the case. 
c...  Two index quantities of alpha spin are denoted with a _a, 
c...  and beta spin with _b. With this in mind, naming of arrays is obvious
c...  when Chem. Phys. Letts. 166 (1990) 275 is consulted.
c     
      if(task.eq.'gradient')then
*ga:1:0
         if(.not.ga_create(mt_dbl,noa,noa,'g_p_ij_a',
     $        noa,0,g_p_ij_a))call errquit(
     $        'could not galloc g_p_ij_a',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nob,nob,'g_p_ij_b',
     $        nob,0,g_p_ij_b))call errquit(
     $        'could not galloc g_p_ij_b',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nva,nva,'g_p_ab_a',
     $        nva,0,g_p_ab_a))call errquit(
     $        'could not galloc g_p_ab_a',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nvb,nvb,'g_p_ab_b',
     $        nvb,0,g_p_ab_b))call errquit(
     $        'could not galloc g_p_ab_b',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,noa,noa,'g_w_ij_a',
     $        noa,0,g_w_ij_a))call errquit(
     $        'could not galloc g_w_ij_a',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nob,nob,'g_w_ij_b',
     $        nob,0,g_w_ij_b))call errquit(
     $        'could not galloc g_w_ij_b',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nva,nva,'g_w_ab_a',
     $        nva,0,g_w_ab_a))call errquit(
     $        'could not galloc g_w_ab_a',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nvb,nvb,'g_w_ab_b',
     $        nvb,0,g_w_ab_b))call errquit(
     $        'could not galloc g_w_ab_b',1, GA_ERR)
c
*ga:1:0
         if(.not.ga_create(mt_dbl,noa,nva,'lai_a',noa,0,g_lai_a))
     $        call errquit('could not alloc ga for lai_a',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nob,nvb,'lai_b',nob,0,g_lai_b))
     $        call errquit('could not alloc ga for lai_b',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,noa,nva,'wai_a',noa,0,g_w_ai_a))
     $        call errquit('could not alloc ga for wai_a',1, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nob,nvb,'wai_b',nob,0,g_w_ai_b))
     $        call errquit('could not alloc ga for wai_b',1, GA_ERR)
c
         if (noa_frozen .gt. 0) then
*ga:1:0
            if(.not.ga_create(mt_dbl,noa-noa_frozen,noa_frozen,'lyi_a',
     $           0,0,g_liy_a))
     $           call errquit('could not alloc ga for lyi_a',1, GA_ERR)
*ga:1:0
            if(.not.ga_create(mt_dbl,nob-nob_frozen,nob_frozen,'lyi_b',
     $           0,0,g_liy_b))
     $           call errquit('could not alloc ga for lyi_b',1, GA_ERR)
c
            call ga_zero(g_liy_a)
            call ga_zero(g_liy_b)
         endif
c     
         call ga_zero(g_p_ij_a)
         call ga_zero(g_p_ab_a)
         call ga_zero(g_p_ij_b)
         call ga_zero(g_p_ab_b)
         call ga_zero(g_w_ij_a)
         call ga_zero(g_w_ij_b)
         call ga_zero(g_w_ab_a)
         call ga_zero(g_w_ab_b)
         call ga_zero(g_w_ai_a)
         call ga_zero(g_w_ai_b)
         call ga_zero(g_lai_a)
         call ga_zero(g_lai_b)
      end if
c     
      pass=1
      call mp2_grad_uhf(rtdb,g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     int_mb(k_irs_a),sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,
     $     num_oa_sym,num_va_sym,g_vecs_b,nob_lo,nob_hi,nvb_lo,
     $     nvb_hi,int_mb(k_irs_b),sym_lo_ob,sym_hi_ob,sym_lo_vb,
     $     sym_hi_vb,num_ob_sym,num_vb_sym,

     $     sym_lo_foa,sym_hi_foa,num_foa_sym,
     $     sym_lo_fob,sym_hi_fob,num_fob_sym,

     $     dbl_mb(k_eval_a),
     $     dbl_mb(k_eval_b),basis,geom,tol2e,e2,esing,etrip,
     $     nbf,nmo,task,
     $     pass,oblk,oskel,osym,
     $     g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,g_liy_a,g_liy_b,
     $     g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b,
     $     g_w_ai_a,g_w_ai_b,g_lai_a,g_lai_b,
     $     dbl_mb(k_grad),scftype)

      ess= etrip / 3.0d0
 
      if(scftype.eq.'UHF') then
         pass=2
         esing = 0.0d0
         etrip = 0.0d0
         call mp2_grad_uhf(rtdb,g_vecs_b,nob_lo,nob_hi,
     $        nvb_lo,nvb_hi,int_mb(k_irs_b),sym_lo_ob,
     $        sym_hi_ob,sym_lo_vb,sym_hi_vb,num_ob_sym,
     $        num_vb_sym,g_vecs_a,noa_lo,noa_hi,nva_lo,
     $        nva_hi,int_mb(k_irs_a),sym_lo_oa,sym_hi_oa,sym_lo_va,
     $        sym_hi_va,num_oa_sym,num_va_sym,

     $        sym_lo_fob,sym_hi_fob,num_fob_sym,
     $        sym_lo_foa,sym_hi_foa,num_foa_sym,

     $        dbl_mb(k_eval_b),
     $        dbl_mb(k_eval_a),basis,geom,tol2e,e2,esing,etrip,
     $        nbf,nmo,task,
     $        pass,oblk,oskel,osym,
     $        g_p_ij_b,g_p_ij_a,g_p_ab_b,g_p_ab_a,g_liy_b,g_liy_a,
     $        g_w_ij_b,g_w_ij_a,g_w_ab_b,g_w_ab_a,
     $        g_w_ai_b,g_w_ai_a,g_lai_b,g_lai_a,
     $        dbl_mb(k_grad),scftype)
         ess = ess + etrip / 3.0d0
         eos = e2 - ess
         esing = 0.0d0
         etrip = 0.0d0
      else
c     
c     For RHF basically just need to add the mixed and pure
c     spin pieces together and then set alpha=beta.  The
c     energy and non-separable gradient need doubling.
c     
         e2 = e2 * 2.0d0
         ess= ess * 2.0d0
         eos = e2 - ess
c     
         if (task .eq. 'gradient') then
            call dscal(3*natoms, 2.0d0, dbl_mb(k_grad), 1)
            call ga_dadd(1.0d0, g_p_ij_a, 1.0d0, g_p_ij_b, g_p_ij_a)
            call ga_copy(g_p_ij_a, g_p_ij_b)
            call ga_dadd(1.0d0, g_p_ab_a, 1.0d0, g_p_ab_b, g_p_ab_a)
            call ga_copy(g_p_ab_a, g_p_ab_b)
            call ga_dadd(1.0d0, g_w_ij_a, 1.0d0, g_w_ij_b, g_w_ij_a)
            call ga_copy(g_w_ij_a, g_w_ij_b)
            call ga_dadd(1.0d0, g_w_ab_a, 1.0d0, g_w_ab_b, g_w_ab_a)
            call ga_copy(g_w_ab_a, g_w_ab_b)
            call ga_dadd(1.0d0, g_w_ai_a, 1.0d0, g_w_ai_b, g_w_ai_a)
            call ga_copy(g_w_ai_a, g_w_ai_b)
            call ga_dadd(1.0d0, g_lai_a, 1.0d0, g_lai_b, g_lai_a)
            call ga_copy(g_lai_a, g_lai_b)
            if (noa_frozen .gt. 0) then
               call ga_dadd(1.0d0, g_liy_a, 1.0d0, g_liy_b, g_liy_a)
               call ga_copy(g_liy_a, g_liy_b)
            endif
         end if
      end if
 1984 continue
c     
c     Tally up the total energy
c
      if (dftmp2) then
        if(.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,escf))
     $       call errquit('could not get dft energy',1, RTDB_ERR)
        if (.not.rtdb_get(rtdb,'dft:mp2fac',mt_dbl,1,mp2fac)) then
          mp2fac = 1.0d0
          if (ga_nodeid().eq.0) then
            write(LuOut,*) 'dft:pt2fac not found, using 1.0'
          endif
        endif
      else
        if(.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,escf))
     $       call errquit('could not get scf energy',1, RTDB_ERR)
      endif
c
        if (dftmp2) then
          if (scs) then
            ! DHDF with SCS scaling
            e2      = fss*ess+fos*eos
            e_total = escf + mp2fac*e2
      if (doprint) 
     W      write(6,996) escf, ess, fss, eos, fos, e2, e_total
 996        format(//,10x,43('-'),
     $              /,10x,'DFT energy                 ',f24.12,
     $              /,10x,'Same spin pairs            ',f24.12,
     $              /,10x,'Same spin scaling factor   ',f24.12,
     $              /,10x,'Opposite spin pairs        ',f24.12,
     $              /,10x,'Opposite spin scaling fact.',f24.12,
     $              /,10x,'SCS-MP2 correlation energy ',f24.12,
     $              /,10x,'Total DFT+SCS-MP2 energy   ',f24.12,
     $              /,10x,43('-'),/)
          else
            ! DHDF with no SCS scaling
            e_total = escf + e2*mp2fac
            if(doprint) write(6,998) escf, e2, e_total
 998        format(//,10x,43('-'),
     $              /,10x,'DFT energy             ',f24.12,
     $              /,10x,'Unscaled MP2 energy    ',f24.12,
     $              /,10x,'Total DFT+MP2 energy   ',f24.12,
     $              /,10x,43('-'),/)
          endif
        else ! HF-based MP2
          e_total = escf + e2
          if(doprint)write(6,999) escf, e2, esing, etrip, e_total
 999      format(//,10x,43('-'),
     $            /,10x,'SCF energy         ',f24.12,
     $            /,10x,'Correlation energy ',f24.12,
     $            /,10x,'Singlet pairs      ',f24.12,
     $            /,10x,'Triplet pairs      ',f24.12,
     $            /,10x,'Total MP2 energy   ',f24.12,
     $            /,10x,43('-'),/)

          ! there should really be an if (scs) here...

          ! SCS-MP2

          ! why are we not setting e2 and e_total for SCS?
          ! (same for esing and etrip)
          ! these are written to RTDB and ECCE...
          if(doprint)write(6,997) ess, fss, eos, fos, fss*ess+fos*eos,
     $                 escf+fss*ess+fos*eos
 997      format(//,10x,51('-'),
     $            /,20x,'Spin Component Scaled (SCS) MP2'
     $            /,10x,'Same spin pairs            ',f24.12,
     $            /,10x,'Same spin scaling factor   ',f24.12,
     $            /,10x,'Opposite spin pairs        ',f24.12,
     $            /,10x,'Opposite spin scaling fact.',f24.12,
     $            /,10x,'SCS-MP2 correlation energy ',f24.12,
     $            /,10x,'Total SCS-MP2 energy       ',f24.12,
     $            /,10x,51('-'),/)
        endif
        if(doprint) call util_flush(6)

      if(scs) then
c     store scs energies for gradient
         e_total=escf+fss*ess+fos*eos
         e2=fss*ess+fos*eos
      endif
c
      call ecce_print1('total energy', mt_dbl, e_total, 1)
      call ecce_print1('correlation energy', mt_dbl, e2, 1)
      if(.not. rtdb_put(rtdb,'mp2:correlation energy',mt_dbl,1,e2))
     $     call errquit('semi_dir_mp2:failed storing energy in rtdb',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'mp2:energy',mt_dbl,1,e_total))
     $     call errquit('semi_dir_mp2:failed storing energy in rtdb',0,
     &       RTDB_ERR)
      if (dftmp2) then
        if(.not. rtdb_put(rtdb,'dft:energy',mt_dbl,1,e_total))
     1       call errquit('semi_dir_mp2:failed storing energy in rtdb',
     2                    0,RTDB_ERR)
      endif
      if(.not. rtdb_put(rtdb,'mp2:esing',mt_dbl,1,esing))
     $     call errquit('semi_dir_mp2:failed storing energy in rtdb',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'mp2:etrip',mt_dbl,1,etrip))
     $     call errquit('semi_dir_mp2:failed storing energy in rtdb',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'mp2:escf',mt_dbl,1,escf))
     $     call errquit('semi_dir_mp2:failed storing energy in rtdb',0,
     &       RTDB_ERR)
c     
      if(task.eq.'gradient') then
         if (omatrix) then
            call ga_print(g_p_ab_a)
            call ga_print(g_p_ij_a)
            call ga_print(g_w_ab_a)
            call ga_print(g_w_ij_a)
            call ga_print(g_w_ai_a)
            call ga_print(g_lai_a)
            if (noa_frozen .gt. 0) call ga_print(g_liy_a)
            call ga_print(g_p_ab_b)
            call ga_print(g_p_ij_b)
            call ga_print(g_w_ab_b)
            call ga_print(g_w_ij_b)
            call ga_print(g_w_ai_b)
            call ga_print(g_lai_b)
            if (nob_frozen .gt. 0) call ga_print(g_liy_b)
         end if
         if (otrace) then
            call mp2_print_trace_norm(g_w_ij_a,'w_ij_a',noa)
            call mp2_print_trace_norm(g_w_ab_a,'w_ab_a',nva)
            call mp2_print_trace_norm(g_w_ij_b,'w_ij_b',nob)
            call mp2_print_trace_norm(g_w_ab_b,'w_ab_b',nvb)
            norm = sqrt(ga_ddot(g_lai_a,g_lai_a))
            if (ga_nodeid().eq.0) write(6,771) 'Norm of lai_a', norm
            norm = sqrt(ga_ddot(g_lai_b,g_lai_b))
            if (ga_nodeid().eq.0) write(6,771) 'Norm of lai_b', norm
            norm = sqrt(ga_ddot(g_w_ai_a,g_w_ai_a))
            if (ga_nodeid().eq.0) write(6,771) 'Norm of w_ai_a', norm
            norm = sqrt(ga_ddot(g_w_ai_b,g_w_ai_b))
            if (ga_nodeid().eq.0) write(6,771) 'Norm of w_ai_b', norm
 771        format(1x,a,1x,f20.12)
            if (ga_nodeid().eq.0) call util_flush(6)
         end if
         call ga_sync()
c     
         if (noa_frozen .gt. 0) then
            call mp2_make_piy(noa_lo, noa_hi, g_p_ij_a, g_liy_a, 
     $           dbl_mb(k_eval_a))
            call mp2_make_piy(nob_lo, nob_hi, g_p_ij_b, g_liy_b, 
     $           dbl_mb(k_eval_b))
            call mp2_make_wiy(noa_lo, noa_hi, g_w_ij_a, g_liy_a, 
     $           dbl_mb(k_eval_a))
            call mp2_make_wiy(nob_lo, nob_hi, g_w_ij_b, g_liy_b, 
     $           dbl_mb(k_eval_b))
            if (omatrix) then
               call ga_print(g_p_ij_a)
               call ga_print(g_p_ij_b)
               call ga_print(g_w_ij_a)
               call ga_print(g_w_ij_b)
            endif
         endif
         call pstat_on(ps_laifock)
         if (doprint) then
            write(6,111) ' Start laifock', util_wallsec()
            call util_flush(6)
         end if
         call mp2_lai_fock_uhf_prepar(
     $        g_p_ab_a,g_p_ab_b,g_p_ij_a,g_p_ij_b,
     $        noa,nob,nva,nvb,g_tmpa, g_tmpb,nmo)
         call mp2_lai_fock_uhf(geom,basis,
     $        g_tmpa,g_tmpb,g_vecs_a,g_vecs_b,
     $        noa,nob,nva,nvb,g_lai_a,g_lai_b,rtdb,tol2e_fock)
         call mp2_lai_fock_uhf_tidy(g_tmpa, g_tmpb)
         call pstat_off(ps_laifock)
         if (doprint) then
            write(6,111) ' Done laifock', util_wallsec()
            call util_flush(6)
         end if
c
         if (omatrix) then
            call ga_print(g_lai_a)
            call ga_print(g_lai_b)
         endif
c
         transform=.false.
c     
         if (scftype .eq. 'UHF') then
            grad_length=(noa*nva)+(nob*nvb)
         else
            grad_length=noa*nva
         end if
c     
*ga:1:0
         if(.not.ga_create(mt_dbl,grad_length,1,'mp2g: rhs',
     $        0,0,g_rhs))call errquit('mp2: ga_rhs',0, GA_ERR)
c     
         call mp2_make_cphf_rhs(scftype, transform,basis,nmo,nbf,
     $        noa,nva,nob,nvb,int_mb(k_irs_a),int_mb(k_irs_b),
     $        int_mb(k_map_a),int_mb(k_map_b),
     $        num_oa_sym,num_ob_sym,num_va_sym,num_vb_sym,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,
     $        sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $        g_vecs_a,g_vecs_b,g_lai_a,g_lai_b,g_rhs)
c     
         if (omatrix) call ga_print(g_rhs)
         if (otrace) then
            norm = sqrt(ga_ddot(g_rhs,g_rhs))
            if (ga_nodeid().eq.0) write(6,771) 'Norm of RHS', norm
         end if
c     
         call pstat_on(ps_cphf)
         if(.not.cphf(rtdb))call errquit('could not do cphf',1,
     &       CALC_ERR)
         call pstat_off(ps_cphf)
c     
c...  the g_rhs will now hold the solution         
c     
         call cphf_fname('cphf_sol',cphf_sol)
         if(.not.file_read_ga(cphf_sol,g_rhs))
     $        call errquit('could not read cphf_rhs',1, DISK_ERR)
c     
         if (omatrix) call ga_print(g_rhs)
         if (otrace) then
            norm = sqrt(ga_ddot(g_rhs,g_rhs))
            if (ga_nodeid().eq.0) write(6,771) 'Norm of soln', norm
         end if
c     
*ga:1:0
         if(.not.ga_create(mt_dbl,nmo,nmo,'tot alpha dens',nmo,0,
     $        g_p_tot_a))call errquit(
     $        'mp2_energy: ga_create failed',0, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nmo,nmo,'tot beta dens',nmo,0,
     $        g_p_tot_b))call errquit(
     $        'mp2_energy: ga_create failed',0, GA_ERR)
c
         if (doprint) then
            write(6,111) ' Start make_tot_dens', util_wallsec()
            call util_flush(6)
         end if
         
         call mp2_make_tot_dens(g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,
     $        g_rhs,g_p_tot_a,g_p_tot_b,noa,nva,nob,nvb,nbf,nmo,
     $        int_mb(k_map_a), int_mb(k_map_b),
     $        g_vecs_a,g_vecs_b,basis,scftype)
         if (doprint) then
            write(6,111) ' Done make_tot_dens', util_wallsec()
            call util_flush(6)
         end if
c
         if (otrace) then
            call mp2_print_trace_norm(g_p_tot_a,'g_p_tot_a',nmo)
            call mp2_print_trace_norm(g_p_tot_b,'g_p_tot_b',nmo)
         endif
         if (omatrix) then
            call ga_print(g_p_tot_a)
            call ga_print(g_p_tot_b)
         endif
c     
         call mp2_make_dip_mom(rtdb, basis, nbf, nmo, noa, nob,
     $        g_p_tot_a, g_p_tot_b, g_vecs_a, g_vecs_b)
c
         call mp2_make_nos(rtdb, basis, nbf, nmo, noa, nob,
     $        g_p_tot_a, g_p_tot_b, g_vecs_a, g_vecs_b,
     $        scftype, noa_frozen, title_vecs)
c     
         call pstat_on(ps_wijfock)
         call mp2_wij_fock_uhf(rtdb, geom,basis, tol2e_fock,
     $        g_p_tot_a,g_p_tot_b,noa,nob,g_vecs_a,g_vecs_b,
     $        g_w_ij_a,g_w_ij_b)
         call pstat_off(ps_wijfock)
c     
         if (omatrix) then
            call ga_print(g_w_ij_a)
            call ga_print(g_w_ij_b)
         end if
         if (otrace) then
            call mp2_print_trace_norm(g_w_ij_a,'w_ij_a',noa)
            call mp2_print_trace_norm(g_w_ij_b,'w_ij_b',nob)
         end if
c     
*ga:1:0
         if(.not.ga_create(mt_dbl,nmo,nmo,'tot alpha wdens',nmo,0,
     $        g_w_tot_a))call errquit(
     $        'mp2_energy: ga_create failed',0, GA_ERR)
*ga:1:0
         if(.not.ga_create(mt_dbl,nmo,nmo,'tot beta wdens',nmo,0,
     $        g_w_tot_b))call errquit(
     $        'mp2_energy: ga_create failed',0, GA_ERR)
         call mp2_make_tot_wdens(g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b,
     $        g_w_ai_a,g_w_ai_b,g_w_tot_a,g_w_tot_b,
     $        noa_lo,noa_hi,nva,nob_lo,nob_hi,
     $        nvb,nbf,nmo,g_vecs_a,g_vecs_b,
     $        dbl_mb(k_eval_a),dbl_mb(k_eval_b),
     $        g_p_tot_a,g_p_tot_b)
c     
         call ga_dgop(1,dbl_mb(k_grad),3*natoms,'+')
c     
         if(ga_nodeid().eq.0)then
            if (util_print('nonsepgrad',print_high)) then
               write(6,*)
               write(6,*) ' Non-separable gradient'
               write(6,*)
               write(6,*) '       X              Y              Z'
               write (luout,1000) (dbl_mb(k_grad+i),i=0,3*natoms-1)
 1000          format (1000(3(F14.8,1X)/))
               call util_flush(6)
            end if
            do i=1, 3*natoms
               dbl_mb(k_grad+i-1) = -dbl_mb(k_grad+i-1)
            end do
            status=rtdb_parallel(.false.)
            if(.not.rtdb_put(rtdb,'mp2:nonseparable gradient',mt_dbl,
     $           3*natoms,dbl_mb(k_grad)))
     $           call errquit(
     $           'semi_dir_mp2:failed to store gradient in rtdb',0,
     &       RTDB_ERR)
            status=rtdb_parallel(.true.)
         end if
c
         if (noa_frozen .gt. 0) then
            if (.not. ga_destroy(g_liy_a)) call errquit(
     $           'could not destroy Liy a',0, GA_ERR)
            if (.not. ga_destroy(g_liy_b)) call errquit(
     $           'could not destroy Liy b',0, GA_ERR)
         endif
         if(.not.ga_destroy(g_rhs))call errquit(
     $        'could not destroy g_rhs',1, GA_ERR)
         if(.not.ga_destroy(g_p_ij_a))call errquit(
     $        'could not destroy p_ij_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_p_ij_b))call errquit(
     $        'could not destroy p_ij_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_p_ab_a))call errquit(
     $        'could not destroy p_ab_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_p_ab_b))call errquit(
     $        'could not destroy p_ab_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_p_tot_a))call errquit(
     $        'could not destroy p_tot_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_p_tot_b))call errquit(
     $        'could not destroy p_tot_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_ij_a))call errquit(
     $        'could not destroy g_w_ij_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_ij_b))call errquit(
     $        'could not destroy g_w_ij_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_ab_a))call errquit(
     $        'could not destroy g_w_ab_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_ab_b))call errquit(
     $        'could not destroy g_w_ab_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_ai_a))call errquit(
     $        'could not destroy g_w_ai_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_ai_b))call errquit(
     $        'could not destroy g_w_ai_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_tot_a))call errquit(
     $        'could not destroy g_w_tot_a handle',1, GA_ERR)
         if(.not.ga_destroy(g_w_tot_b))call errquit(
     $        'could not destroy g_w_tot_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_lai_b))call errquit(
     $        'could not destroy lai_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_lai_a))call errquit(
     $        'could not destroy lai_a handle',1, GA_ERR)
      end if
c     
      if(.not.ma_chop_stack(l_grad))call errquit(
     $     'could not chop stack in mp2_energy',1, MA_ERR)
      if(.not.ga_destroy(g_vecs_a))call errquit(
     $     'could not destroy g_vecs_a',1, GA_ERR)
      if(.not.ga_destroy(g_vecs_b))
     $     call errquit('could not destroy g_vecs_b',1, GA_ERR)
c...  
c...  call to Herbert's stuff to do the separable peice
c...  
      if (task .eq. 'gradient') then
         call pstat_on(ps_sep)
         call grad_force(rtdb,basis,geom) 
         call pstat_off(ps_sep)
c
c     Clean up the files (CPHF RHS and solution, densities)
c
         call cphf_fname('cphf_sol',cphf_sol)
         if (ga_nodeid() .eq. 0) call util_file_unlink(cphf_sol)
         call cphf_fname('cphf_rhs',cphf_sol)
         if (ga_nodeid() .eq. 0) call util_file_unlink(cphf_sol)
         call util_file_name('p_a_ao',.true.,.true.,cphf_sol)
         if (ga_nodeid() .eq. 0) call util_file_unlink(cphf_sol)
         call util_file_name('p_b_ao',.true.,.true.,cphf_sol)
         if (ga_nodeid() .eq. 0) call util_file_unlink(cphf_sol)
         call util_file_name('w_a_ao',.true.,.true.,cphf_sol)
         if (ga_nodeid() .eq. 0) call util_file_unlink(cphf_sol)
         call util_file_name('w_b_ao',.true.,.true.,cphf_sol)
         if (ga_nodeid() .eq. 0) call util_file_unlink(cphf_sol)
      endif
c     
      call ga_sync
 111  format(a,' at ', f10.1,'s')
c     
      end
      subroutine mp2_print_trace_norm(g, text, n)
      implicit none
#include "global.fh"
      integer g                 ! [input] GA handle
      character*(*) text        ! [input] info to print
      integer n                 ! [input] dimension of ga
c     
      double precision sum, xx, norm
      integer i
c     
      norm = sqrt(ga_ddot(g,g))
      if (ga_nodeid() .eq. 0) then
         sum = 0.0d0
         do i = 1, n
            call ga_get(g,i,i,i,i,xx,1)
            sum = sum + xx
         end do
         write(6,1) text,sum,norm
 1       format(1x,a,' trace=',f20.12,' norm=',f20.12)
         call util_flush(6)
      end if
c     
      end
