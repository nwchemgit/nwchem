      subroutine cphf_solve(rtdb)
*
* $Id: cphf_solve.F,v 1.31 2001-06-28 03:19:10 zzhang Exp $
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
c!!! BGJ
#include "bgj.fh"
c!!! BGJ
      integer rtdb              ! [input] database handle
      double precision eone,etwo,enrep,energy
      integer g_grad, g_work, g_evecs
      integer max_precond_iter,max_precond_sub
      double precision precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      logical ga_iter_diag
      external ga_iter_diag
      logical oprint_solve
c$$$  logical ga_it_ls_mrhs
c$$$  external ga_it_ls_mrhs
      external rohf_hessv_precond,rohf_hessv
      external uhf_precond,uhf_hessv
      external tddft_precond
      logical oprint, oprint_conv
c     
      character*(nw_max_path_len) cphf_rhs
      character*(nw_max_path_len) cphf_sol
      logical file_write_ga,file_read_ga
      external file_write_ga,file_read_ga
      integer vlen, nvirt, nmos(2), nbf_vecs, nsets
      integer bases(2), nbases
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      double precision residual
      logical movecs_read_header, cphf_poliz
      external movecs_read_header, cphf_poliz
      logical int_normalize
      external int_normalize
c
      logical hf_job
      integer nroot, leval, ieval, iroot
      logical cphf_tddft, cphf_triplet
      cphf_tddft = .true. 
      cphf_triplet = .false. ! true for triplet, false for singlet
c
      hf_job = .not. bgj_have_xc()      
c
c     Print levels assume that CPHF is being invoked from within
c     another application.  By default be informative.
c     At low print, just indicate convergence or failure.
c
      oprint = util_print('information',print_default) .and.
     $     ga_nodeid() .eq. 0
      oprint_conv = util_print('information',print_low) .and.
     $     ga_nodeid() .eq. 0
      oprint_solve = util_print('solve',print_default)
c
c     Must use the same vectors as the converged SCF
c
      if (hf_job) then
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('cphf: no SCF MO vectors',0)
      else
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('cphf_solve: no DFT MO vectors',0)
c !!! Hack to set up scftype properly for later code !!! ? need this?
         if (bgj_restricted()) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
      endif
c      write(*,*)'*** cphf_solve: scftype ',scftype
c     
      call scf_init(rtdb)
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('cphf_solve: int_normalize failed',911)
      nbases = 1
      bases(1) = basis
      if (bgj_have_j_fit()) then
        nbases = 2
        bases(2) = bgj_CD_bas_han()
        write(*,*) __FILE__,': bases',bases
      endif
      call int_init(rtdb,nbases,bases)
      call schwarz_init(geom, basis)
      oskel=sym_number_ops(geom).gt.0
      call int_1e_cache_ga(basis,oskel)
c
c     Get the MO vectors
c
      if(.not.movecs_read_header(movecs_in,title_vecs,basis_vecs,
     $     scftype_vecs, nbf_vecs,nsets,nmos,2))
     $     call errquit('cphf: problem with movecs_read_header',1)
      if(nbf.ne.nbf_vecs)call errquit('cphf: MO vectors mismatch?',0)
      nmo = nmos(1)
c
      if (nmo .ne. nbf) then
         if (.not. ga_destroy(g_movecs(1))) call errquit('cphf: ga',1)
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs a',
     $        0, 0, g_movecs(1))) call errquit('scf_v_g: tmp MOs', 0)
         if (scftype .eq. 'UHF') then
            if (.not. ga_destroy(g_movecs(2)))call errquit('cphf: ga',2)
*ga:1:0
            if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs b',
     $           0, 0, g_movecs(2))) call errquit('scf_v_g: tmp MOs', 0)
         endif
      endif
c         
      call scf_movecs_read
c
      if (scftype.eq.'UHF') then 
         vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         nvirt = nmo - nclosed - nopen
         vlen = nclosed*(nopen+nvirt) + nopen*nvirt
      else
         call errquit('cphf: unknown SCF type',0)
      endif
c
      lshift=0.0d0
      if (.not. rtdb_get(rtdb, 'cphf:thresh', mt_dbl, 1, 
     $     precond_tol)) precond_tol=1d-4
      if (.not. rtdb_get(rtdb, 'cphf:maxiter', mt_int, 1,
     $     max_precond_iter)) max_precond_iter = 100
c
      pflg=2
      max_precond_sub = 20
      if (oprint) then
         write(6,170) scftype
 170     format('    scftype = ',a)
         if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
            write(6,171) nclosed, nopen
 171        format('    nclosed = ', i4/
     $             '      nopen = ', i4)
         else
            write(6,172) nalpha, nbeta
 172        format('     nalpha = ', i4/
     $             '      nbeta = ', i4)
         endif
         write(6,173) vlen, precond_tol, max_precond_iter
 173     format('  variables = ', i8/
     $          '  tolerance = ', 1p,d9.1,0p/
     $          '  max. iter = ', i4/)
      endif
c
*ga:1:0
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve: rhs',0,1,g_grad))
     $     call errquit('cphf_solve: ga_cre grad',vlen)
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)         
      if(.not.file_read_ga(cphf_rhs, g_grad))
     $     call errquit('cphf_solve: could not read cphf_rhs',1)
c
*ga:1:0
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve: grad',
     $     0,0,g_work))call errquit('cphf_solve: ga_cre grad',0)
c
c     Compute necessary Fock matrices
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 1))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0)
c !!! BGJ
      if(scftype .eq. 'UHF')then
         call uhf_init(rtdb)
         if(.not. rtdb_put(rtdb,'cphf_solve:cphf_uhf',MT_LOG,1,.true.))
     $        call errquit('cphf_solve: rtdb_put of calc_type failed',0)
         call uhf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
         if(.not.rtdb_delete(rtdb,'cphf_solve:cphf_uhf'))
     $        call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      else
         call rohf_init(rtdb)
         call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work,
     ,        rtdb)
      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',1)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
c !!! BGJ
c
      residual = sqrt(ga_ddot(g_work, g_work))
      if (residual .gt. 1d-4) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' SCF residual = ', residual
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('cphf_solve:the available MOs do not satisfy'//
     $        ' the SCF equations',0)
      endif
c
      call ga_zero(g_work)      ! Essential to get default initial guess
c
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'cphf_solve:cphf_tddft',
     &     MT_log,1,cphf_tddft))
     $     call errquit('cphf_solve: rtdb_put of cphf_tddft failed',0)
      if(.not. rtdb_put(rtdb,'cphf_solve:cphf_triplet',
     &     MT_log,1,cphf_triplet))
     $     call errquit('cphf_solve: rtdb_put of cphf_triplet failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0)
c !!! BGJ
      if (scftype .eq. 'UHF') then
         if(cphf_tddft)then
            nroot = 1
            if(.not.ga_create(mt_dbl,vlen,nroot,'cphf_solve: g_evecs',
     $           0,0,g_evecs))call errquit('cphf_solve: ga_cre evecs',0)
            call ga_zero(g_evecs)
            do iroot = 1, nroot
               call tddft_guess(g_grad, g_work, iroot)
               call ga_copy_patch('n', g_work, 1, vlen, 1, 1,
     &              g_evecs, 1, vlen, iroot, iroot)
            enddo
            if (.not.MA_Push_Get(mt_dbl, nroot, 'Eval', lEval,iEval))
     &           call errquit('cphf_solve: cannot allocate Eval',0)
            if(.not.ga_iter_diag(vlen, nroot, max_precond_iter,
     $           max_precond_sub, precond_tol, 
     $           tddft_precond, uhf_hessv, oprint_solve, 0.d0, 
     $           g_evecs, dbl_mb(ieval), precond_acc, precond_iter))
     $           call errquit('cphf did not converge',precond_iter)
            if(.not.ma_pop_stack(leval))
     &           call errquit("can not pop leval?", 0)
            if(.not.ga_destroy(g_evecs))call errquit
     &           ("ga_iter_diag: can not destroy g_evecs", 0)
         else
            if(.not.ga_iter_lsolve(vlen, max_precond_iter,
     $           max_precond_sub, precond_tol, 
     $           uhf_precond, uhf_hessv, .false., oprint_solve, 
     $           g_grad, g_work, precond_acc, precond_iter))
     $           call errquit('cphf did not converge',precond_iter)
         endif
      else 
         if(cphf_tddft)then
            nroot = 1
            if(.not.ga_create(mt_dbl,vlen,nroot,'cphf_solve: g_evecs',
     $           0,0,g_evecs))call errquit('cphf_solve: ga_cre evecs',0)
            call ga_zero(g_evecs)
            do iroot = 1, nroot
               call tddft_guess(g_grad, g_work, iroot)
               call ga_copy_patch('n', g_work, 1, vlen, 1, 1,
     &              g_evecs, 1, vlen, iroot, iroot)
            enddo
            if (.not.MA_Push_Get(mt_dbl, nroot, 'Eval', lEval,iEval))
     &           call errquit('cphf_solve: cannot allocate Eval',0)
            if(.not.ga_iter_diag(vlen, nroot, max_precond_iter,
     $           max_precond_sub, precond_tol, 
     $           tddft_precond, rohf_hessv, oprint_solve, 0.d0, 
     $           g_evecs, dbl_mb(ieval), precond_acc, precond_iter))
     $           call errquit('cphf did not converge',precond_iter)
            if(.not.ma_pop_stack(leval))
     &           call errquit("can not pop leval?", 0)
            if(.not.ga_destroy(g_evecs))call errquit
     &           ("ga_iter_diag: can not destroy g_evecs", 0)
         else
            if(.not.ga_iter_lsolve(vlen, max_precond_iter,
     $           max_precond_sub, precond_tol, 
     $           rohf_hessv_precond, rohf_hessv, .false., oprint_solve, 
     $           g_grad, g_work, precond_acc,precond_iter))
     $           call errquit('cphf did not converge', precond_iter)
         endif
      endif
c     !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0)
c !!! BGJ
c
      if (oprint_conv) then
         write(6,174) precond_acc, precond_iter
 174     format(/' CPHF converged to ',1p,d9.1,0p,' in ', i4, ' iters'/)
         call util_flush(6)
      endif
c
c$$$  if(.not.ga_it_ls_mrhs(cuhf_vlen,1,max_iter,
c$$$  $        max_sub,accuracy, 
c$$$  $        uhf_precond,uhf_hessv,oprint, 
c$$$  $        g_grad,g_work,precond_acc,precond_iter))
c$$$  $        call errquit('cphf did not converge',1)
c$$$  pol=ga_ddot(g_work,g_grad)
c$$$  if(ga_nodeid().eq.0)then
c$$$  write(6,911) pol
c$$$  911        format(//,15x,40('-'),
c$$$  $           /,15x,'UHF Polarizability is:',f16.6,
c$$$  $           /,15x,40('-'),/)
c$$$  endif
c
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_write_ga(cphf_sol,g_work))
     $     call errquit('cphf_solve: could not write cphf_sol',1)
      if(.not.ga_destroy(g_work))call errquit(
     $     'cphf_solve: could not destroy work array',1)
      if(.not.ga_destroy(g_grad))call errquit(
     $     'cphf_solve: could not destroy gradient handle',1)
c
      call int_1e_uncache_ga()
      call schwarz_tidy()
      call int_terminate()
      if(scftype .eq. 'UHF') then
         call uhf_tidy(rtdb)
      else if(scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         call rohf_tidy
      endif
      if(cphf_tddft)call hess_tidy()
      call scf_tidy(rtdb)
c
      return
c
c     if (.not. cphf_poliz(rtdb))     ! Never executed.
c    $     call errquit(' cphf: not possible!',0)
c
      end

      subroutine tddft_guess(g_grad, g_work, flag)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "cscf.fh"
      integer g_grad, g_work, flag
      integer vlen, i, j, ioff, nvira, l_init, k_init
      double precision ea, ei, eia
c
      call ga_zero(g_grad)
      call ga_zero(g_work)
c
c     initial guess of g_soln
c
      if (scftype.eq.'UHF') then 
         vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         vlen = nclosed*(nmo-nclosed)
      else
         call errquit('tddft_guess: unknown SCF type',0)
      endif

      if(.not.ma_push_get(mt_dbl,vlen,'initial guess',l_init,
     &     k_init))
     &     call errquit('cannot alloctate l_init',0)      
      call dfill(vlen, 0.0d0, dbl_mb(k_init), 1)
      if(scftype .eq. 'UHF')then
         nvira=nmo-nalpha
c         ioff = nvira*(nalpha-2)
         ioff = nvira*(nalpha-1)
         if(flag .eq. 1)then
            dbl_mb(k_init+ ioff) = 1.0d0 ! homo to lumo
            dbl_mb(k_init + ioff + nvira*nalpha) = -1.0d0 ! beta
         else if (flag.eq.2) then
            dbl_mb(k_init+ ioff + nvira*nalpha ) = 1.0d0
         else if (flag.eq.3) then
            dbl_mb(k_init+ ioff + 1) = 1.0d0
         else if (flag.eq.4) then
            dbl_mb(k_init+ ioff + nvira*nalpha +1) = 1.0d0
         endif
      else if(scftype .eq. 'RHF') then
         ioff = (nmo-nclosed)*(nclosed-1)
         if(flag .eq. 1)then
            dbl_mb(k_init+ ioff) = 1.0d0 ! homo to lumo
         else if (flag.eq.2) then
            dbl_mb(k_init+ ioff + 1) = -1.0d0
         else if (flag.eq.3) then
            dbl_mb(k_init+ ioff + 2) = 1.0d0
         else if (flag.eq.4) then
            dbl_mb(k_init+ ioff +3) = -1.0d0
         endif
      else
         call errquit('tddft_guess: unknown SCF type',0)
      endif
c     
c     scale the initial guess by epsilon(i,a)=epsilon(i) - epsilon(a)
c     
      if(scftype .eq. 'UHF')then 
         do i = 1, nalpha
            ei = dbl_mb(k_eval+(i-1))
            do j = 1, nbf-nalpha
               ea = dbl_mb(k_eval+(nalpha+j-1))
               eia = ea - ei 
               dbl_mb(k_init + ((i-1)*(nbf-nalpha)+j-1))=
     &              dbl_mb(k_init + ((i-1)*(nbf-nalpha)+j-1)) /
     &              sqrt(eia)
            enddo
         enddo
         ioff = nalpha*(nmo-nalpha)
         do i = 1, nbeta
            ei = dbl_mb(k_eval+nbf+(i-1))
            do j = 1, nbf-nbeta
               ea = dbl_mb(k_eval+nbf+(nbeta+j-1))
               eia = ea - ei 
               dbl_mb(k_init + ioff + ((i-1)*(nbf-nbeta)+j-1))=
     &              dbl_mb(k_init + ioff + ((i-1)*(nbf-nalpha)+j-1)) / 
     &              sqrt(eia)
            enddo
         enddo
      else if(scftype.eq.'RHF')then
         do i = 1, nclosed
            ei = dbl_mb(k_eval+(i-1))
            do j = 1, nbf-nclosed
               ea = dbl_mb(k_eval+(nalpha+j-1))
               eia = ea - ei 
               dbl_mb(k_init + ((i-1)*(nbf-nclosed)+j-1))=
     &              dbl_mb(k_init + ((i-1)*(nbf-nclosed)+j-1)) /
     &              sqrt(eia)
            enddo
         enddo
      else
         call errquit('tddft_guess: unknown SCF type',0)
      endif
      call ga_put(g_work,1,vlen,1,1,dbl_mb(k_init),1)
      if(.not.ma_pop_stack(l_init))
     &     call errquit("cphf_guess: cannot pop stack initial guess", 0) 
c     call ga_print(g_work)
      
      return 
      end

      subroutine tddft_scale(g_x, g_grad, flag)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "cscf.fh"
      integer g_grad, g_x, g_a, g_tmp, g_tmp2
      integer vlen, gtype, nvec, i, j, ioff, nvira, l_init, k_init
      double precision ea, ei, eia
      logical flag              ! flag = .false. : scale by sqrt(ea-ei) only
                                ! flag = .true.  : add the orbital energy and scale
      if(scftype.eq.'UHF')then
         vlen = nalpha*(nmo-nalpha)+nbeta*(nmo-nbeta)
      else if(scftype.eq.'RHF')then
         vlen=nclosed*(nmo-nclosed)
      else
         call errquit("unknown scftype?", 0)
      endif
      call ga_inquire(g_x,gtype,vlen,nvec)
      if(.not.ma_push_get(mt_dbl,vlen,'tmp_vec',l_init,
     &     k_init))
     &     call errquit('cannot alloctate l_init',0)      
      call dfill(vlen, 0.0d0, dbl_mb(k_init), 1)
      if(scftype.eq.'UHF')then
         nvira=nmo-nalpha
         ioff = nvira*(nalpha-1)
c     
c     form diagonal matrix of (sqrt of) eigenvalue differences 
c     
         do i = 1, nalpha
            ei = dbl_mb(k_eval+(i-1))
            do j = 1, nbf-nalpha
               ea = dbl_mb(k_eval+(nalpha+j-1))
               eia = ea - ei 
               dbl_mb(k_init + ((i-1)*(nbf-nalpha)+j-1)) = sqrt(eia)
            enddo
         enddo
         ioff = nalpha*(nmo-nalpha)
         do i = 1, nbeta
            ei = dbl_mb(k_eval+nbf+(i-1))
            do j = 1, nbf-nbeta
               ea = dbl_mb(k_eval+nbf+(nbeta+j-1))
               eia = ea - ei 
               dbl_mb(k_init + ioff + ((i-1)*(nbf-nbeta)+j-1))
     &              = sqrt(eia)
            enddo
         enddo
      else if (scftype.eq.'RHF')then
c     
c     form diagonal matrix of (sqrt of) eigenvalue differences 
c     
         do i = 1, nclosed
            ei = dbl_mb(k_eval+(i-1))
            do j = 1, nbf-nclosed
               ea = dbl_mb(k_eval+(nclosed+j-1))
               eia = ea - ei 
               dbl_mb(k_init + ((i-1)*(nbf-nclosed)+j-1)) = sqrt(eia)
            enddo
         enddo
      else
         call errquit("unknown scftype?",0)
      endif
c     
      if (.not. ga_create(MT_DBL, vlen, vlen, 'tddft_precond: diag',
     $     1, 0, g_a)) call errquit('diag', 0)
      if (.not. ga_create(MT_DBL, vlen, nvec, 'tddft_precond: tmp',
     $     1, 0, g_tmp)) call errquit('diag', 0)
      if (.not. ga_create(MT_DBL, vlen, vlen, 'tddft_precond: tmp',
     $     1, 0, g_tmp2)) call errquit('diag', 0)
      call ga_zero(g_a)
      call ga_zero(g_tmp)
c     
      call ga_copy_todiag(g_a,dbl_mb(k_init))
      if(flag)then
c     
c     form the diagonal eia matrix g_tmp = g_a*g_a and then add the
c     matrix vector product g_tmp * g_x to g_grad.
c     
         if(scftype.eq.'UHF')then
            call ga_dgemm('n', 'n', vlen, vlen, vlen, 1.0d0, g_a, g_a, 
     &           0.0d0, g_tmp2)
         else if(scftype.eq.'RHF')then
            call ga_dgemm('n', 'n', vlen, vlen, vlen, 1.0d0, g_a, g_a, 
     &           0.0d0, g_tmp2)
         endif
         call ga_dgemm('n', 'n', vlen, nvec, vlen, 1.0d0, g_tmp2, g_x, 
     &        1.0d0, g_grad)
      endif
      call ga_copy(g_grad, g_tmp)
      call ga_dgemm('n', 'n', vlen, nvec, vlen, 1.0d0, g_a, g_tmp, 
     &     0.0d0, g_grad) 
      if(.not.ga_destroy(g_a))
     &     call errquit("tddft_scale: destroy g_a?" ,0)
      if(.not.ga_destroy(g_tmp))
     &     call errquit("tddft_scale: destroy g_tmp?" ,0)
      if(.not.ga_destroy(g_tmp2))
     &     call errquit("tddft_scale: destroy g_tmp2?" ,0)
      if(.not.ma_pop_stack(l_init))
     &     call errquit("cphf_guess: cannot pop stack initial guess", 0) 
      return 
      end

