      subroutine mp2_energy(rtdb,geom,basis,movecs,task)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "stdio.fh"
#include "sym.fh"
      character*8 task
      integer num_oa_sym(0:7),num_va_sym(0:7)    
      integer num_ob_sym(0:7),num_vb_sym(0:7)
      integer sym_lo_oa(0:7)    ! [output] First MO in each sym
      integer sym_hi_oa(0:7)    ! [output] Last MO in each sym
      integer sym_lo_va(0:7)
      integer sym_hi_va(0:7)
      integer sym_lo_ob(0:7)    ! [output] First MO in each sym
      integer sym_hi_ob(0:7)    ! [output] Last MO in each sym
      integer sym_lo_vb(0:7)
      integer sym_hi_vb(0:7)
      integer nmixed
      integer rtdb,basis,geom
      integer g_vecs_a,g_vecs_b,nalpha,nbeta
      integer nocc
      character*4 scftype
      integer i,nbf_vecs,nsets,nmo,nbf
      integer l_eval_a,k_eval_a,l_eval_b,k_eval_b
      integer noa_frozen,nva_frozen,nob_frozen,nvb_frozen
      integer noa_lo,noa_hi,nob_lo,nob_hi,nva_lo,nva_hi,nvb_lo,nvb_hi
      integer noa_correl,nob_correl
      integer nmo2_a,nmo2_b
      integer k_irs_a,k_irs_b,l_irs_a,l_irs_b
      integer k_occ_a,k_occ_b,l_occ_a,l_occ_b
      integer pass
      integer noa,nob,nva,nvb
      integer noper_a,noper_b
      double precision tol2e,e2
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
      character*255 movecs
      character*255 title_vecs,basis_vecs
      integer maxireps,maxops
      parameter(maxireps=20,maxops=120)
      integer k_junk,l_junk
      integer k_work,l_work
      integer k_map_a,l_map_a
      integer k_map_b,l_map_b
      intrinsic nint
      logical status,oblk
      integer k_p_ij_a,k_p_ij_b
      integer l_p_ij_a,l_p_ij_b
      integer k_p_ab_a,k_p_ab_b
      integer l_p_ab_a,l_p_ab_b
      integer g_p_ab_a,g_p_ab_b,g_p_ij_a,g_p_ij_b,g_lai_a,g_lai_b
      integer g_p_ai_a,g_p_ai_b,g_dens_tot_a,g_dens_tot_b
c     
      integer g_rhs
      integer grad_length
      integer g_p_tot_a,g_p_tot_b
      logical transform
      logical cphf
      external cphf
      logical file_write_ga,file_read_ga
      external file_write_ga,file_read_ga
      character*256 cphf_sol
      integer cphf_sol_u
      parameter(cphf_sol_u=18)
      integer g_x_a,g_x_b,g_y_a,g_y_b,g_z_a,g_z_b,g_scratch
c
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
     $     call errquit('scf:scftype wrong',1)
      status=bas_numbf(basis,nbf)
      if(.not.status) call errquit('mp2_grad: could not get nbf',0)
      if(scftype.eq.'UHF')then
         if(.not.rtdb_get(rtdb,'scf:nalpha',mt_int,1,nalpha))
     $        call errquit('no alphas',1)
         if(.not.rtdb_get(rtdb,'scf:nbeta',mt_int,1,nbeta))
     $        call errquit('no betas',1)
         noa=nalpha
         nva=nbf-noa
         nob=nbeta
         nvb=nbf-nob
      endif
      if(scftype.eq.'RHF')then
         if(.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nocc))
     $        call errquit('no occ',1)
         noa=nocc
         nob=nocc
         nva=nbf-noa
         nvb=nbf-nob
      endif
c     
      if(.not.ga_create(mt_dbl,nbf,nbf,'alpha_vecs',nbf,1,
     $     g_vecs_a)) call errquit('mp2_energy: ga_create failed',0)
      if(scftype.eq.'UHF')then
         if(.not.ga_create(mt_dbl,nbf,nbf,'beta_vecs',nbf,1,g_vecs_b))
     $        call errquit('mp2_energy: ga_create failed',0)
      endif
      if(.not.ma_push_get(mt_dbl,nbf,'occ alpha',
     $     l_occ_a,k_occ_a))      
     $     call errquit('mp2_energy ma problem with alpha occ',1)
      if(scftype.eq.'UHF')then
         if(.not.ma_push_get(mt_dbl,nbf,'occ beta',
     $        l_occ_b,k_occ_b))
     $        call errquit('mp2_energy ma problem with beta occ',1)
      endif
      if(.not.movecs_read_header(movecs,title_vecs,basis_vecs,
     $     nbf_vecs,nsets,nmo,2))
     $     call errquit('problem with movecs_read_header',1)
      if(.not.ma_push_get(mt_dbl,nbf,'alpha evals',
     $     l_eval_a,k_eval_a))
     $     call errquit('problem with ma for a evals',1)
      if(scftype.eq.'UHF')then
         if(.not.ma_push_get(mt_dbl,nbf,'beta evals',
     $        l_eval_b,k_eval_b))
     $        call errquit('problem with ma for b evals',1)
      endif
      if(nbf.eq.nbf_vecs)then
         if(.not.movecs_read(movecs,1,dbl_mb(k_occ_a),
     $        dbl_mb(k_eval_a),g_vecs_a))
     $        call errquit('could not read alpha vecs',1)
         if(scftype.eq.'UHF')then
            if(.not.movecs_read(movecs,2,dbl_mb(k_occ_b),
     $           dbl_mb(k_eval_b),g_vecs_b))
     $           call errquit('could not read beta vecs',1)
         endif
      endif
      if(.not.rtdb_get(rtdb,'mp2_grad:alpha_occ_frozen',mt_int,1,
     $     noa_frozen))then
         if(.not.rtdb_get(rtdb,'mp2_grad:no_alpha_occ_frozen',
     $        mt_int,1,noa_frozen))noa_frozen=0
      endif
      if(.not.rtdb_get(rtdb,'mp2_grad:alpha_vir_frozen',mt_int,1,
     $     nva_frozen))then
         if(.not.rtdb_get(rtdb,'mp2_grad:no_alpha_vir_frozen',
     $        mt_int,1,nva_frozen))nva_frozen=0
      endif
      if(scftype.eq.'UHF')then
         if(.not.rtdb_get(rtdb,'mp2_grad:beta_occ_frozen',mt_int,1,
     $        nob_frozen))then
            if(.not.rtdb_get(rtdb,'mp2_grad:no_beta_occ_frozen',
     $           mt_int,1,nob_frozen))nob_frozen=0
         endif
         if(.not.rtdb_get(rtdb,'mp2_grad:beta_vir_frozen',mt_int,1,
     $        nvb_frozen))then
            if(.not.rtdb_get(rtdb,'mp2_grad:no_beta_vir_frozen',
     $           mt_int,1,nvb_frozen))nvb_frozen=0
         endif
      endif
      if(.not.rtdb_get(rtdb,'fourindex:aoblock',mt_log,oblk,1))
     $     oblk=.false.
      if(.not.ma_push_get(mt_int,nbf,'irs_a',l_irs_a,k_irs_a))
     $     call errquit('scf_movecs_sym_adapt: memory?',nbf)
      if(scftype.eq.'UHF')then
         if(.not.ma_push_get(mt_int,nbf,'irs_b',l_irs_b,k_irs_b))
     $        call errquit('scf_movecs_sym_adapt: memory?',nbf)
      endif
      call sym_movecs_adapt(basis,1d-8,g_vecs_a,int_mb(k_irs_a),
     $     nmixed)
      if(scftype.eq.'UHF')then
         call sym_movecs_adapt(basis,1d-8,g_vecs_b,int_mb(k_irs_b),
     $        nmixed)      
      endif
      do i=1,nbf
         int_mb(k_irs_a+i-1)=int_mb(k_irs_a+i-1)-1
         if(scftype.eq.'UHF')
     $        int_mb(k_irs_b+i-1)=int_mb(k_irs_b+i-1)-1
      enddo
      noa_lo=noa_frozen+1
      noa_hi=noa
      noa_correl=noa_hi-noa_lo+1
      nva_lo=noa+1
      nva_hi=nbf-nva_frozen
      nmo2_a=nva_hi-nva_lo+1
      noper_a=noa_hi-noa_lo+1
      if(scftype.eq.'UHF')then
         nob_lo=nob_frozen+1
         nob_hi=nob
         nob_correl=nob_hi-noa_lo+1
         nvb_lo=nob+1
         nvb_hi=nbf-nvb_frozen
         nmo2_b=nvb_hi-nvb_lo+1
         noper_b=nob_hi-nob_lo+1
      endif
c     
      if(.not.ma_push_get(mt_dbl,nbf,'work',l_work,k_work))
     $     call errquit('could not alloc work',1)
      if(.not.ma_push_get(mt_int,nbf,'map_a',l_map_a,k_map_a))
     $     call errquit('could not alloc map_a',1)
      if(scftype.eq.'UHF')then
         if(.not.ma_push_get(mt_int,nbf,'map_b',l_map_b,k_map_b))
     $        call errquit('could not alloc map_b',1)
      endif
      call moints_vecs_sym_sort(g_vecs_a,nbf,noa_lo,noa_hi,
     $     int_mb(k_irs_a),int_mb(k_map_a),dbl_mb(k_work),num_oa_sym,
     $     sym_lo_oa,sym_hi_oa)
      if(scftype.eq.'UHF')
     $     call moints_vecs_sym_sort(g_vecs_b,nbf,nob_lo,nob_hi,
     $     int_mb(k_irs_b),int_mb(k_map_b),dbl_mb(k_work),num_ob_sym,
     $     sym_lo_ob,sym_hi_ob)
      call moints_vecs_sym_sort(g_vecs_a,nbf,nva_lo,nva_hi,
     $     int_mb(k_irs_a),int_mb(k_map_a),dbl_mb(k_work),num_va_sym,
     $     sym_lo_va,sym_hi_va)
      if(scftype.eq.'UHF')
     $     call moints_vecs_sym_sort(g_vecs_b,nbf,nvb_lo,nvb_hi,
     $     int_mb(k_irs_b),int_mb(k_map_b),dbl_mb(k_work),num_vb_sym,
     $     sym_lo_vb,sym_hi_vb)
      if(.not.ma_push_get(mt_dbl,nbf,'junk',l_junk,k_junk))
     $     call errquit('could not alloc junk',1)
      call dcopy(nbf,dbl_mb(k_eval_a),1,dbl_mb(k_junk),1)
      do i=noa_lo,nva_hi
         dbl_mb(k_eval_a-1+int_mb(k_map_a+i-1))=dbl_mb(k_junk+i-1)
      enddo
      if(scftype.eq.'UHF')then
         call dcopy(nbf,dbl_mb(k_eval_b),1,dbl_mb(k_junk),1)
         do i=nob_lo,nvb_hi
            dbl_mb(k_eval_b-1+int_mb(k_map_b+i-1))=dbl_mb(k_junk+i-1)
         enddo
      endif
c     
      tol2e=1d-12
c     
      e2=0.0d0
      if(scftype.eq.'RHF')then
         call mp2_grad_rhf(g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $        int_mb(k_irs_a),sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,
     $        num_oa_sym,num_va_sym,dbl_mb(k_eval_a),basis,geom,tol2e,
     $        e2,nbf,task)
      endif
c     
c... A note on notation: arrays for four index quantities of pure spin
c... denoted with a _a on the end. Mixed spin case with a _b. On first
c... call to grad_uhf _a arrays are all alpha,_b are <alpha beta||alpha beta>.
c... On second call to grad_uhf vice versa is the case. 
c... Two index quantities of alpha spin are denoted with a _a, 
c... and beta spin with _b. With this in mind, naming of arrays is obvious
c... when Chem. Phys. Letts. 166 (1990) 275 is consulted.
c
      if(.not.ma_push_get(mt_dbl,noa*noa,'p_ij_a',
     $     l_p_ij_a,k_p_ij_a))call errquit(
     $     'could not alloc p_ij_a',1)
      if(.not.ma_push_get(mt_dbl,nva*nva,'p_ab_a',
     $     l_p_ab_a,k_p_ab_a))call errquit(
     $     'could not alloc p_ab_a',1)
      if(.not.ma_push_get(mt_dbl,nob*nob,'p_ij_b',
     $     l_p_ij_b,k_p_ij_b))call errquit(
     $     'could not alloc p_ij_b',1)
      if(.not.ma_push_get(mt_dbl,nvb*nvb,'p_ab_b',
     $     l_p_ab_b,k_p_ab_b))call errquit(
     $     'could not alloc p_ab_b',1)
      if(.not.ga_create(mt_dbl,noa,noa,'g_p_ij_a',
     $     noa,1,g_p_ij_a))call errquit(
     $     'could not alloc g_p_ij_a',1)
      if(.not.ga_create(mt_dbl,nob,nob,'g_p_ij_b',
     $     nob,1,g_p_ij_b))call errquit(
     $     'could not alloc g_p_ij_b',1)
      if(.not.ga_create(mt_dbl,nva,nva,'g_p_ab_a',
     $     nva,1,g_p_ab_a))call errquit(
     $     'could not alloc g_p_ab_a',1)
      if(.not.ga_create(mt_dbl,nvb,nvb,'g_p_ab_b',
     $     noa,1,g_p_ab_b))call errquit(
     $     'could not alloc g_p_ab_b',1)
      if(.not.ga_create(mt_dbl,noa,nva,'lai_a',noa,1,g_lai_a))
     $     call errquit('could not alloc ga for lai_a',1)
      if(.not.ga_create(mt_dbl,nob,nvb,'lai_b',nob,1,g_lai_b))
     $     call errquit('could not alloc ga for lai_b',1)
      call ga_zero(g_p_ij_a)
      call ga_zero(g_p_ab_a)
      call ga_zero(g_p_ij_b)
      call ga_zero(g_p_ab_b)
      call ga_zero(g_lai_a)
      call ga_zero(g_lai_b)
      if(scftype.eq.'UHF')then
         pass=1
         call mp2_grad_uhf(g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $        int_mb(k_irs_a),sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,
     $        num_oa_sym,num_va_sym,g_vecs_b,nob_lo,nob_hi,nvb_lo,
     $        nvb_hi,int_mb(k_irs_b),sym_lo_ob,sym_hi_ob,sym_lo_vb,
     $        sym_hi_vb,num_ob_sym,num_vb_sym,dbl_mb(k_eval_a),
     $        dbl_mb(k_eval_b),basis,geom,tol2e,e2,nbf,task,
     $        pass,oblk,
     $        dbl_mb(k_p_ij_a),dbl_mb(k_p_ij_b),dbl_mb(k_p_ab_a),
     $        dbl_mb(k_p_ab_b),g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,
     $        g_lai_a,g_lai_b)
         pass=2
         call mp2_grad_uhf(g_vecs_b,nob_lo,nob_hi,nvb_lo,nvb_hi,
     $        int_mb(k_irs_b),sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $        num_ob_sym,num_vb_sym,g_vecs_a,noa_lo,noa_hi,nva_lo,
     $        nva_hi,int_mb(k_irs_a),sym_lo_oa,sym_hi_oa,sym_lo_va,
     $        sym_hi_va,num_oa_sym,num_va_sym,dbl_mb(k_eval_b),
     $        dbl_mb(k_eval_a),basis,geom,tol2e,e2,nbf,task,
     $        pass,oblk,
     $        dbl_mb(k_p_ij_b),dbl_mb(k_p_ij_a),dbl_mb(k_p_ab_b),
     $        dbl_mb(k_p_ab_a),g_p_ij_b,g_p_ij_a,g_p_ab_b,g_p_ab_a,
     $        g_lai_b,g_lai_a)
      endif
      if(.not.rtdb_put(rtdb,'semi_dir_mp2:energy',mt_dbl,1,e2))
     $     call errquit('could not stash the energy',1)
c     
      if(task.eq.'gradient')call mp2_lai_fock_uhf(geom,basis,
     $     g_p_ab_a,g_p_ab_b,g_p_ij_a,g_p_ij_b,g_vecs_a,g_vecs_b,
     $     noa,nob,nva,nvb,g_lai_a,g_lai_b,rtdb,tol2e)
c
      if(task.eq.'gradient')then
         transform=.false.
      grad_length=(noa*nva)+(nob*nvb)
      if(.not.ga_create(mt_dbl,grad_length,1,'uhf_nr_solve: grad',
     $     1,1,g_rhs))call errquit('mp2: ga_grad',0)
         call mp2_make_cphf_rhs(transform,basis,nmo,nbf,
     $        noa,nva,nob,nvb,g_vecs_a,g_vecs_b,g_lai_a,g_lai_b,g_rhs)
      endif
      if(.not.cphf(rtdb))call errquit('could not do cphf',1)
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_read_ga(cphf_sol,g_rhs))
     $     call errquit('could not read cphf_rhs',1)
      if(.not.ga_create(mt_dbl,nbf,nbf,'tot alpha dens',nbf,1,
     $     g_p_tot_a))call errquit('mp2_energy: ga_create failed',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,'tot beta dens',nbf,1,
     $     g_p_tot_b))call errquit('mp2_energy: ga_create failed',0)
      call mp2_make_tot_dens(g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,
     $     g_rhs,g_p_tot_a,g_p_tot_b,noa,nva,nob,nvb,nbf,
     $     g_vecs_a,g_vecs_b,geom,basis,rtdb)
      call mp2_make_dip_mom(g_p_tot_a,g_p_tot_b,basis,rtdb,nbf,
     $     g_vecs_a,g_vecs_b)
c     
      if(ga_nodeid().eq.0)write(6,*) 'MP2 Correction = ', e2
c     
      if(.not.ga_destroy(g_rhs))call errquit(
     $     'could not destroy g_rhs',1)
      if(.not.ga_destroy(g_p_ij_a))call errquit(
     $     'could not destroy p_ij_a handle',1)
      if(.not.ga_destroy(g_p_ij_b))call errquit(
     $     'could not destroy p_ij_b handle',1)
      if(.not.ga_destroy(g_p_ab_a))call errquit(
     $     'could not destroy p_ab_a handle',1)
      if(.not.ga_destroy(g_p_ab_b))call errquit(
     $     'could not destroy p_ab_b handle',1)
      if(.not.ga_destroy(g_p_tot_a))call errquit(
     $     'could not destroy p_tot_a handle',1)
      if(.not.ga_destroy(g_p_tot_b))call errquit(
     $     'could not destroy p_tot_b handle',1)
      if(.not.ga_destroy(g_lai_b))call errquit(
     $     'could not destroy lai_b handle',1)
      if(.not.ga_destroy(g_lai_a))call errquit(
     $     'could not destroy lai_a handle',1)
c     
      if(.not.ma_chop_stack(l_occ_a))call errquit(
     $     'could not chop stack in mp2_energy',1)
      if(.not.ga_destroy(g_vecs_a))call errquit(
     $     'could not destroy g_vecs_a',1)
      if(scftype.eq.'UHF')then
         if(.not.ga_destroy(g_vecs_b))
     $        call errquit('could not destroy g_vecs_b',1)
      endif
      end
