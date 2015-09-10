      subroutine tce_grad_force(rtdb,model,basis,geom,
     1       d_1pdm, k_1pdm_offset, size_1pdm,
     1       d_1pwdm,k_1pwdm_offset,size_1pwdm,
     1       d_2pdm, k_2pdm_offset, size_2pdm,
     1       noab,nvab,k_range,k_offset,k_spin,
     1       k_movecs_sorted,restricted)
c     $Id$
C     calculate energy gradients with respect to nuclear coordinates
C------------------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "apiP.fh"
#include "inp.fh"

      integer rtdb, basis,geom
c
      character*10 model
      integer d_1pdm,  k_1pdm_offset,  size_1pdm
      integer d_1pwdm, k_1pwdm_offset, size_1pwdm
      integer d_2pdm,  k_2pdm_offset,  size_2pdm
      integer noab,nvab,k_range,k_offset,k_movecs_sorted,k_spin
      integer k_2pdm_ao,l_2pdm_ao,k,l
      logical restricted
      integer g_1pdm,  g_1pwdm
      integer d_2pdm_ao,l_2pdm_ao_offset,k_2pdm_ao_offset,size_2pdm_ao
      integer idim,jdim,kdim,ldim,i_ind,j_ind,k_ind,l_ind,i1,j1,k1,l1
      integer ib,jb,kb,lb
      character*255 filename
      double precision res
c     ====================
c     atomic orbital tiles
c     ====================
      integer atpart2
      integer nalength2(200)
      integer a2length2(200)
c
      integer g1b
      integer g1
      integer l_a,k_a,l_dim_shb,k_dim_shb
      integer l_mo_tmp,k_mo_tmp,dima
      
c
      double precision cpu_tim(2), wall_tim(2)
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_dens(6),        ! density matrices(up to 6)
     $        g_wdens,          ! energy weighted density
     $        g_eigen_diag,     ! eigenvalue diagonal matrix or lagrangian
     $        g_force           ! forces on atoms(cartesian)
c
      character*255 movecs_in
c
      integer nproc, nat, nbf, nsh,max_sh_bf, max_at_bf, nopen,
     &        nclosed, 
     $        ndens, nbf2, nocc, max1e, max2e, mscratch_1e, mscratch_2e,
     $        lbuf, lscratch, lsqa, lsqatom, i, j, iat
c
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap, 
     $        l_bfglo, l_bfghi, l_labels, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap, 
     $        k_bfglo, k_bfghi, k_labels, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi
      integer k_shbflo, k_shbfhi
c      
      integer lforce, l_force, k_force, l_evals, k_evals, l_occ, k_occ,
     $        l_act, k_act, 
     $        l_buf, k_buf, l_scr, k_scr, l_dens, k_dens, 
     1        l_wdens, k_wdens,
     $        k_frc_nuc, k_frc_kin, k_frc_2el, k_frc_wgh, k_frc_cd,
     $        l_frc_nuc, l_frc_kin, l_frc_2el, l_frc_wgh, l_frc_cd,
     $        l_frc_mp2, k_frc_mp2

c     pdm2d is now used by all methods.  pdm2/a/b/c only by MCSCF
      integer l_pdm2, l_pdm2a, l_pdm2b, l_pdm2c, l_pdm2d ! MCSCF 2-pdm
      integer k_pdm2, k_pdm2a, k_pdm2b, k_pdm2c, k_pdm2d ! MCSCF 2-pdm
      integer l_coeff, k_coeff  ! MCSCF local copy of Active space MOs

      integer nactive, nshblocks

      integer blen              ! bf-blocking size for shell ordering
      integer maxblen           ! maximum value for blen = max no. bf in group
      integer maxsh             ! maximum no. of shells in a group
      integer maxq              ! max quartets in a request
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)
c
      integer 
     $     lh_ij,  ld_ij, lh_kl, ld_kl, lh_ik, ld_ik, 
     $     lh_jl,  ld_jl, lh_il, ld_il,
     $     lh_jk,  ld_jk, 
     $     lh_ij2, ld_ij2, lh_kl2, ld_kl2, lh_ik2, ld_ik2,
     $     lh_jl2, ld_jl2, lh_il2, ld_il2, lh_jk2, ld_jk2,
     $     lh_ij3, ld_ij3, lh_kl3, ld_kl3, lh_ik3, ld_ik3,
     $     lh_jl3, ld_jl3, lh_il3, ld_il3, lh_jk3, ld_jk3,
     $     lh_ij4, ld_ij4, lh_kl4, ld_kl4, lh_ik4, ld_ik4,
     $     lh_jl4, ld_jl4, lh_il4, ld_il4, lh_jk4, ld_jk4

      double precision crd(3),  ! atomic coordinates
     $     tol2e, q

      double precision grad_norm, grad_max

      external grad_norm, grad_max

      character*16 tag
      character*32 theory
      character*8  scftype
      character*32 rtdb_string

c    
c     Stuff for Douglas-Kroll
c
      logical  dkdiv_energy, num_grad, file_write_ga  ! Three data blocks needed
      external dkdiv_energy, num_grad, file_write_ga  ! for the possible addition
      integer  l_dkdivfrc,   k_dkdivfrc               ! of Douglas-Kroll contribution
      character*(nw_max_path_len) dkdiv_name
      character*128 mytheory
c
*     Stuff for DFT

      integer ipol, numfunc, noc(2)
      integer ld_ik5,ld_jl5,ld_il5,ld_jk5,
     $        ld_ik6,ld_jl6,ld_il6,ld_jk6,
     $        ld_ik7,ld_jl7,ld_il7,ld_jk7,
     $        ld_ik8,ld_jl8,ld_il8,ld_jk8

c
c     cosmo flag
c
      logical odbug
      logical ocosmo
      logical osome

      parameter (numfunc = 40)
      logical oskel, omp2, odft, ocdfit, status,frac_occ
      double precision xfac(numfunc), jfac, kfac
c->pdfan
      logical otce
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
c<-pdfan

      nproc = ga_nnodes()

      jfac = 1.0d0              ! Only changed for DFT
      kfac = 1.0d0
C     get information about basis set
      if (.not. geom_ncent(geom,nat))
     $     call errquit('tce_grad_force: could not get natoms',0, 
     2     GEOM_ERR)
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('tce_grad_force: could not get nbf',0, 
     1     BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('tce_grad_force: could not get nsh',0, 
     1     BASIS_ERR)
c     
c     Atom blocking now only used for the 1-e integrals (history)
c     
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('tce_grad_force: could not get max_at_bf',0, 
     1     BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('tce_grad_force: could not get max_sh_bf',0, 
     1     BASIS_ERR)

      blen = min(nbf,maxblen,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets

      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('tce_grad_force: failed getting theory',0, 
     1     RTDB_ERR)

C     get SCF MO vectors for density

      odft = .false.
      if (theory .eq. 'mcscf') then
         if (.not. rtdb_cget(rtdb, 'mcscf:input vectors', 1, movecs_in))
     $      call errquit('tce_grad_force: MCSCF MO vectors not defined',
     1      0,RTDB_ERR)
         scftype = 'mcscf'
      else if (theory .eq. 'dft') then
         odft = .true.
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('tce_grad_force: DFT MO vectors not defined',
     1        0,RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'dft:cdfit', mt_log, 1, ocdfit))
     $        ocdfit = .false.
         if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     $        call errquit('tce_grad_force: xfac not in rtdb ', 0,
     &        RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     $        ipol = 1
         if (ipol .eq. 1) then
            scftype = 'rhf'
         else
            scftype = 'uhf'
         endif
         jfac = 1.0d0
         if (ocdfit) jfac = 0.0d0
         kfac = xfac(1)
      else
         if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype))
     $        call errquit('tce_grad_force: scftype not defined',0,
     &        RTDB_ERR)
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('tce_grad_force: SCF MO vectors not defined',
     1        0,RTDB_ERR)
      end if

C     get information about type of calculation

      if ((theory .eq. 'mp2') 
     $     .or.(theory .eq. 'semi_dir_mp2')
     $     .or.(theory .eq. 'direct_mp2')) then
         omp2 = .true.
      else
         omp2 = .false.
      end if
c
c->pdfan
      if (theory.eq.'tce') then
        otce = .true.
      else
        otce = .false.
      endif
c<-pdfan
c
C     scftype: MCSCF, RHF, ROHF or UHF

      if (scftype .eq. 'mcscf') then
         if (.not. rtdb_get(rtdb, 'mcscf:nclosed', mt_int, 1, nclosed))
     $        nclosed = 0
         if (.not. rtdb_get(rtdb, 'mcscf:nact', mt_int, 1, nopen))
     $        call errquit('tce_grad_force: no mcscf active orbitals?',
     1        0,RTDB_ERR)
      else if (odft) then
         if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &        call errquit('tce_grad_force: rtdb_get of noc failed', 0,
     &       RTDB_ERR)
         if (ipol .eq. 1) then
            nclosed = noc(1)
            nopen   = 0
         else
            nclosed = min(noc(1),noc(2))
            nopen   = abs(noc(1) - noc(2))
         endif
      else
         if (.not. rtdb_get(rtdb, 'scf:nopen', mt_int, 1, nopen))
     $        nopen = 0
         if (.not. rtdb_get(rtdb, 'scf:nclosed', mt_int, 1, nclosed))
     $        call errquit('tce_grad_force: get of nclosed? ', 0, 
     1        RTDB_ERR)
         
C     is it MP2?
c         if (.not. omp2) then
          if((.not.omp2).and.(.not.otce)) then
            if (ga_nodeid() .eq. 0) then
               if (util_print('information', print_medium)) then
                  if (odft) then
                     write(luout,101) ' dft ', ipol
  101                format(/'  wavefunction    =   ', a, ' ipol =',i2/)
                  else
                     write(luout,1) scftype
    1                format(/'  wavefunction    =   ', a/)
                  endif
                  call util_flush(luout)
               end if
            end if
         end if
c         
      end if
c
c->pengdong.fan
c tce
      if (otce) then
         if (ga_nodeid().eq.0) then
            if (util_print('information',print_medium)) then
               if (model.eq.'ccsd') then
                   write(luout,11) scftype,'CCSD'
   11              format(/'  wavefunction    =   TCE ', a,a/)
               else if(model.eq.'ccsdt') then
                   write(luout,11) scftype,'CCSDT'
               else if(model.eq.'ccsdtq') then
                   write(luout,11) scftype,'CCSDTQ'
               else if(model.eq.'mbpt2') then
                   write(luout,11) scftype,'MP2'
               else if(model.eq.'mbpt3') then
                   write(luout,11) scftype,'MP3'
               else if(model.eq.'mbpt4') then
                   write(luout,11) scftype,'MP4'
               else if(model.eq.'cis') then
                   write(luout,11) scftype,'CIS'
               else if(model.eq.'cisd') then
                   write(luout,11) scftype,'CISD'
               else if(model.eq.'cisdt') then
                   write(luout,11) scftype,'CISDT'
               else if(model.eq.'cisdtq') then
                   write(luout,11) scftype,'CISDTQ'               
               else
                  call errquit('tce_grad_force: not implement yet',109,
     1                 capmis_err)
               endif
            endif
         endif
      endif
c 
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if (omp2) then
         if (scftype .eq. 'UHF') then
            nbf2 = 2 * nbf
            ndens = 5
         else if (scftype .eq. 'RHF') then
            nbf2 = nbf
            ndens = 3
         else                   ! ROHF
            nbf2 = nbf
            ndens = 3
            call errquit('tce_grad_force:no ROMP2 gradients yet', 110, 
     1      CAPMIS_ERR)
         end if
      else                      ! SCF
         if (scftype .eq. 'UHF') then
            nbf2 = 2 * nbf
            ndens = 2
         else if (scftype .eq. 'RHF') then
            nbf2 = nbf
            ndens = 1
         else if (scftype .eq. 'MCSCF') then
            nbf2 = nbf
            ndens = 1
         else                   ! ROHF
            nbf2 = nbf
            ndens = 3
         end if
      end if

C     allocate and initialize global and local memory

C     forces on atoms(3xnat)
*ga:1:0
      if (.not. ga_create(mt_dbl, 3, nat, 'forces', 3, 0, g_force))
     $     call errquit('tce_grad_force: failed to create force GA',0,
     &       GA_ERR)
      call ga_zero(g_force)

C     local replication(separate for the different pieces)
      lforce = nat * 3
      if (.not. ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     $     call errquit('tce_grad_force:could not allocate l_force',
     1     lforce,MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_nuc, 
     $     k_frc_nuc)) call errquit
     &    ('tce_grad_force:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_kin, 
     $     k_frc_kin)) call errquit
     &    ('tce_grad_force:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_wgh, 
     $     k_frc_wgh)) call errquit
     &    ('tce_grad_force:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_2el, 
     $     k_frc_2el)) call errquit
     &    ('tce_grad_force:could not allocate l_force',lforce, MA_ERR)
      if (omp2) then
         if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_mp2, 
     $        k_frc_mp2)) call errquit
     &      ('tce_grad_force:could not allocate l_force',lforce, MA_ERR)
      end if
      if (odft) then
         if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_cd, 
     $        k_frc_cd)) call errquit
     &      ('tce_grad_force:could not allocate l_force',lforce, MA_ERR)
      end if

C     global density
c      do i=1, ndens
c         g_dens(i) = ga_create_atom_blocked(geom, basis, 
c     $        'density matrix')
c      end do
       g_1pdm = ga_create_atom_blocked(geom, basis,
     $        'density matrix')
c     
      if (scftype .eq. 'MCSCF') then
         if (.not. ma_push_get(mt_dbl, nopen**4, 'pdm2', 
     $        l_pdm2 , k_pdm2 )) call errquit
     $        ('tce_grad_force:failed allocating pdm2',nopen**4, MA_ERR)
         if (.not. ma_push_get(mt_dbl, (blen**1)*(nopen**3), 'pdm2a',
     $        l_pdm2a, k_pdm2a)) call errquit
     $        ('tce_grad_force:failed allocating pdm2a',
     1        (blen**1)*(nopen**3),
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, (blen**2)*(nopen**2), 'pdm2b',
     $        l_pdm2b, k_pdm2b)) call errquit
     $        ('tce_grad_force:failed allocating pdm2b',
     1        (blen**2)*(nopen**2),
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, (blen**3)*(nopen**1), 'pdm2c',
     $        l_pdm2c, k_pdm2c)) call errquit
     $        ('tce_grad_force:failed allocating pdm2c',
     1        (blen**3)*(nopen**1),
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, nopen*nbf, 'coeff',
     $        l_coeff, k_coeff)) call errquit
     $        ('tce_grad_force:failed allocating coeff',nopen*nbf,
     1        MA_ERR)
      else
         k_pdm2  = 1            ! To avoid SEGV
         k_pdm2a = 1
         k_pdm2b = 1
         k_pdm2c = 1
         k_coeff = 1
      end if
c     
      if (.not. ma_push_get(mt_dbl, blen**4, 'pdm2d',
     $     l_pdm2d, k_pdm2d)) call errquit
     $     ('tce_grad_force:failed allocating pdm2d',blen**4, MA_ERR)

C     lookup table and list of active atoms
      if (.not. ma_push_get(MT_LOG,nat,'active atoms',l_act,k_act))
     $     call errquit('tce_grad_force: could not allocate l_act',nat,
     1     MA_ERR)

C     symmetry or not
      if (.not. rtdb_get(rtdb, 'gradients:use symmetry', mt_log, 1,
     $     oskel)) then
         if (.not. rtdb_get(rtdb, 'scf:skeleton', mt_log, 1,
     $        oskel)) then
            oskel = sym_number_ops(geom) .gt. 0
         end if
      end if
c     
      if (oskel .and. ga_nodeid().eq.0 .and. 
     $     util_print('information', print_default)) then
         write(luout,*) ' Using symmetry'
         call util_flush(luout)
      end if
c     
C     energy weighted density(NxN)
c      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')
      g_1pwdm = ga_create_atom_blocked(geom, basis, 'weighted density')

C     eigenvalue diagonal matrix for forming energy weighted density
      nocc = nopen + nclosed
*ga:1:0
      if (.not. ga_create(mt_dbl, nocc, nocc,'eigen_diag', 0, nocc, 
     $     g_eigen_diag)) call errquit('tce_grad_force: ga diags?',
     1     nocc*nocc,
     &     GA_ERR)
      call ga_zero(g_eigen_diag)

C     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('tce_grad_force: could not allocate l_evals',
     1     nbf2,
     &     MA_ERR)
C     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf2,'occ. numbers', l_occ, k_occ))
     $     call errquit('tce_grad_force: could not allocate l_occ',
     1     nbf2,
     &     MA_ERR)
c     
      call grad_active_atoms(rtdb, nat, log_mb(k_act), nactive)
c     
c
c     check if smearing is on
c
      if (.not.rtdb_get(rtdb, 'dft:fractional_occup', mt_log, 1,
     & frac_occ  )) frac_occ=.false.
c
c->pengdong.fan
c      call tce_grad_dens(geom, basis, g_dens, g_wdens, g_eigen_diag, 
c     $     dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2, 
c     $     nopen, nclosed, scftype, movecs_in, omp2, 
c     $     dbl_mb(k_pdm2), dbl_mb(k_coeff),frac_occ,.false.)
c
c       call ga_put(g_dens,1,nbf,1,nbf,dbl_mb(k_1pdm_ao),nbf)
c       call ga_put(g_wdens,1,nbf,1,nbf,dbl_mb(k_1pwdm_ao),nbf)
        call ao_tiles(atpart2,nalength2,a2length2,30)
c        call print_1pdm_mo(d_1pdm,k_1pdm_offset,size_1pdm)
        call btrans1(d_1pdm,k_1pdm_offset,g_1pdm,atpart2,nalength2)
        call ga_symmetrize(g_1pdm)
c        call ga_print(g_1pdm)
        call dscalfile(-1.d0,d_1pwdm,size_1pwdm)
c        call print_1pdm_mo(d_1pwdm,k_1pwdm_offset,size_1pwdm)
        call btrans1(d_1pwdm,k_1pwdm_offset,g_1pwdm,atpart2,nalength2)
        call ga_symmetrize(g_1pwdm)
c->debug
c       write(6,*) 'debug g_wdens'
c       call ga_print(g_wdens)
c       write(6,*) 'end of debug g_wdens'
c<-denug
c
c<-pengdong.fan
c            
C     free temporary arrays
      if (.not. ga_destroy(g_eigen_diag))
     $     call errquit('tce_grad_force: could not destroy g_eigen_diag'
     1     ,1,
     &     GA_ERR)

      if (.not.ma_pop_stack(l_occ))
     $     call errquit('tce_grad_force:ma free occ',1, MA_ERR)
      if (.not.ma_pop_stack(l_evals))
     $     call errquit('tce_grad_force:ma free eval',1, MA_ERR)

C     initialize for integral gradients
c->d3p975
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
      call int_terminate()
      call intd_init(rtdb, 1, basis)
c<-d3p975
      call dfill(lforce, 0.0D0, dbl_mb(k_force), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_nuc), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_kin), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_wgh), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el), 1)
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e) ! blocking algorithm
      max2e = max(max2e,1296*100)            ! 100 D quartets 
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)
C     one-electron contribution
C     buffers for one electron integral derivatives

      if (.not. ma_push_get(mt_dbl,lbuf,'deriv buffer',l_buf,k_buf))
     $     call errquit('tce_grad_force:could not allocate buffer',lbuf,
     &       MA_ERR)
*      call ma_summarize_allocated_blocks()
      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch', 
     $     l_scr, k_scr))
     1     call errquit('tce_grad_force: scratch alloc failed',
     $     lscratch, MA_ERR)
      if (.not. ma_push_get(mt_int,lbuf/3,'labels',l_labels,k_labels))
     $     call errquit('tce_grad_force: could not allocate labels',
     1     lbuf/3,
     &     MA_ERR)
c     
C     local density matrix block 
      lsqatom = max_at_bf * max_at_bf
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_dens,k_dens))
     $     call errquit('tce_grad_force:could not allocate l_dens',
     1     lsqatom,
     &     MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_w_density',l_wdens, 
     $     k_wdens)) call errquit
     &    ('tce_grad_force:could not allocate l_wdens',lsqatom, MA_ERR)
c     
      cpu_tim(1)  = util_cpusec()
      wall_tim(1) = util_wallsec()
c      call tce_grad1(dbl_mb(k_buf),lbuf,dbl_mb(k_scr),lscratch,
c     $     dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
c     $     dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh), g_force,
c     $     g_dens, g_wdens, basis, geom, nproc, nat, max_at_bf,
c     $     rtdb, oskel)
c
cc AJL/Begin/SPIN ECPs
cc Added ndens to grad1
cc Work is still needed to separate g_dens for alpha and beta
cc so spin-polarised will not work at this time. It throws
cc an error in tce_grad1. To be looked at later...
c
      call tce_grad1(dbl_mb(k_buf),lbuf,dbl_mb(k_scr),lscratch,
     $     dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
     $     dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh), g_force,
     $     g_1pdm, g_1pwdm, basis, geom, nproc, nat, max_at_bf,
     $     rtdb, oskel, ndens)
c
cc AJL/End
c
c
c     ----- cosmo contribution -----
c
      odbug=.false.
      odbug=odbug.and.ga_nodeid().eq.0
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
         if(odbug) then
            write(luout,*) '-cosmo- ... found in -gradients-',
     $                     ocosmo,ga_nodeid()
         endif
         if(ocosmo) then
            if(odbug) then
               osome=.true.
            else
               osome=.false.
            endif
            osome=osome.and.ga_nodeid().eq.0
            if(odbug) then
               write(luout,*) '-cosmo- ... found and .true. ',
     $                        ocosmo,ga_nodeid()
               write(luout,*) 'calling -cosder- ...'
            endif
c
c     .....
c
c            call grad_hnd_cos( dbl_mb(k_buf), lbuf, dbl_mb(k_scr), 
c     $           lscratch,
c     $           dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
c     $           dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh), g_force,
c     $           g_dens, g_wdens, basis, geom, nproc, nat, max_at_bf,
c     $           rtdb, oskel)
            call grad_hnd_cos( dbl_mb(k_buf), lbuf, dbl_mb(k_scr),
     $           lscratch,
     $           dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
     $           dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh), g_force,
     $           g_1pdm, g_1pwdm, basis, geom, nproc, nat, max_at_bf,
     $           rtdb, oskel)
         else
            if(odbug) then
               write(luout,*) '-cosmo- ... found but .false.',
     $                        ocosmo,ga_nodeid()
            endif
         endif
      else
         if(odbug) then
            write(luout,*) '-cosmo- not found in -gradients-'
         endif
      endif
      call ga_sync()
c
c     Check if we have to deal with Douglas-Kroll type integrals. If so,
c     we have to add the numerical derivative of the integral type
c     (Douglas-Kroll - Non-Relativistic). In this way the non-relativistic
c     kinetic and potential energy integral derivatives are computed
c     analytically and only the relativistic "Douglas-Kroll" correction
c     is done numerically. Numerical errors are reduced in this way and
c     at the same time the construction of the complex Douglas-Kroll
c     integral derivatives is avoided.
c
      if (doug_kroll) then
c
         call intd_terminate()
c
c        Write density to file for use by Douglas-Kroll energy routine
c
         call util_file_name('dkdiv_density',.true.,.true.,dkdiv_name)
c         if (.not. file_write_ga(dkdiv_name,g_dens(1)))
c     $       call errquit('tce_grad_force:file_write_ga density failed',
c     1       0,DISK_ERR)
c
         if (.not. file_write_ga(dkdiv_name,g_1pdm))
     $       call errquit('tce_grad_force:file_write_ga density failed',
     1       0,DISK_ERR)
c
c
c        Do numerical gradient for (Douglas-Kroll - Non-relativistic)
c
         if (ga_nodeid() .eq. 0) then
            call util_print_centered(luout,'Douglas-Kroll gradients',
     $                               40,.true.)
            write(luout,105)
  105       format(/
     $           12x,
     $           'The Douglas-Kroll correction to the gradient is '/
     $           12x,
     $           'computed numerically and added to the analytical '/
     $           12x,
     $           'non-relativistic gradients.'/)
         endif
         mytheory = 'Douglas-Kroll correction to'
         if (.not. num_grad(rtdb,dkdiv_energy,'grad_dk:dkdiv_energy',
     $       'grad_dk:dkdiv_frc', mytheory, .false.))
     $       call errquit('tce_grad_force: task_num_grad failed',0,
     &       UNKNOWN_ERR)
c
c        Add total Douglas-Kroll - Non-relativistic contribution to 
c        array on node 0. There is a global sum later on.
c
         if (ga_nodeid() .eq. 0) then
            status = rtdb_parallel(.false.)
            if (.not. ma_push_get(mt_dbl, 3*nat, 'dkdiv_frc',
     $          l_dkdivfrc, k_dkdivfrc)) call errquit
     $          ('tce_grad_force: ma_push_get dkdiv_frc failed',0,
     1          RTDB_ERR)
            if (.not. rtdb_get(rtdb, 'grad_dk:dkdiv_frc', mt_dbl,
     $          3*nat, dbl_mb(k_dkdivfrc))) call errquit
     $          ('tce_grad_force: could not read grad_dk:dkdiv_frc',0,
     &          RTDB_ERR)
            do i = 0, 3*nat-1
               dbl_mb(k_force+i) = dbl_mb(k_force+i) +
     $                             dbl_mb(k_dkdivfrc+i)
            enddo
            if (.not. ma_pop_stack(l_dkdivfrc))
     $         call errquit('tce_grad_force: ma_pop_stack failed',0,
     1         MA_ERR)
            status = rtdb_parallel(.true.)
         endif
         call util_file_unlink(dkdiv_name)
         call intd_init(rtdb, 1, basis)
      endif
c
      call ga_sync
c
      cpu_tim(1)  = util_cpusec() - cpu_tim(1)
      wall_tim(1) = util_wallsec() - wall_tim(1)
c     
      if (.not. ma_pop_stack(l_wdens)) 
     1     call errquit('tce_grad_force: MA?',0,
     &     MA_ERR)
      if (.not. ma_pop_stack(l_dens)) 
     1     call errquit('tce_grad_force: MA?',0,
     &     MA_ERR)
c      if (.not. ga_destroy(g_wdens)) 
c     1     call errquit('tce_grad_force: GA?',0,
c     &     GA_ERR)
      if (.not. ga_destroy(g_1pwdm))
     1     call errquit('tce_grad_force: GA?',0,
     &     GA_ERR)


C     two-electron contribution
C     allocate arrays for two-electron integral stuff(or rename old ones)
C     blocks of density matrix:

      lsqa = blen*blen

      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij', lh_ij, ld_ij))
     $     call errquit('tce_grad_force:could not allocate ld_ij',lsqa, 
     1     MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl', lh_kl, ld_kl))
     $     call errquit('tce_grad_force:could not allocate ld_kl',lsqa, 
     1     MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik', lh_ik, ld_ik))
     $     call errquit('tce_grad_force:could not allocate ld_ik',lsqa, 
     1     MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl', lh_jl, ld_jl))
     $     call errquit('tce_grad_force:could not allocate ld_jl',lsqa, 
     1     MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il', lh_il, ld_il))
     $     call errquit('tce_grad_force:could not allocate ld_il',lsqa, 
     1     MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk', lh_jk, ld_jk))
     $     call errquit('tce_grad_force:could not allocate ld_jk',lsqa, 
     1     MA_ERR)
c     
      if (scftype .eq. 'UHF' .or. 
     $     scftype .eq. 'ROHF' .or. omp2) then ! UHF or ROHF or MP2
         if (scftype .eq. 'ROHF' .or. omp2) then
            if (.not. ma_push_get(mt_dbl,lsqa,'ld_ij2',lh_ij2,ld_ij2))
     $           call errquit('tce_grad_force:could not allocate ld_ij2'
     1           ,lsqa,
     &           MA_ERR)
            if (.not. ma_push_get(mt_dbl,lsqa,'ld_kl2',lh_kl2,ld_kl2))
     $           call errquit('tce_grad_force:could not allocate ld_kl2'
     1           ,lsqa,
     &           MA_ERR)
         end if          
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik2', lh_ik2, ld_ik2))
     $        call errquit('tce_grad_force:could not allocate ld_ik2',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl2', lh_jl2, ld_jl2))
     $        call errquit('tce_grad_force:could not allocate ld_jl2',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il2', lh_il2, ld_il2))
     $        call errquit('tce_grad_force:could not allocate ld_il2',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk2', lh_jk2, ld_jk2))
     $        call errquit('tce_grad_force:could not allocate ld_jk2',
     1        lsqa,
     &        MA_ERR)
      end if

      if (scftype .eq. 'UHF' .and. omp2) then
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij3', lh_ij3, ld_ij3))
     $        call errquit('tce_grad_force:could not allocate ld_ij3',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl3', lh_kl3, ld_kl3))
     $        call errquit('tce_grad_force:could not allocate ld_jl3',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik3', lh_ik3, ld_ik3))
     $        call errquit('tce_grad_force:could not allocate ld_ik3',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl3', lh_jl3, ld_jl3))
     $        call errquit('tce_grad_force:could not allocate ld_jl3',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il3', lh_il3, ld_il3))
     $        call errquit('tce_grad_force:could not allocate ld_il3',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk3', lh_jk3, ld_jk3))
     $        call errquit('tce_grad_force:could not allocate ld_jk3',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij4', lh_ij4, ld_ij4))
     $        call errquit('tce_grad_force:could not allocate ld_ij4',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl4', lh_kl4, ld_kl4))
     $        call errquit('tce_grad_force:could not allocate ld_kl4',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik4', lh_ik4, ld_ik4))
     $        call errquit('tce_grad_force:could not allocate ld_ik4',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl4', lh_jl4, ld_jl4))
     $        call errquit('tce_grad_force:could not allocate ld_jl4',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il4', lh_il4, ld_il4))
     $        call errquit('tce_grad_force:could not allocate ld_il4',
     1        lsqa,
     &        MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk4', lh_jk4, ld_jk4))
     $        call errquit('tce_grad_force:could not allocate ld_jk4',
     1        lsqa,
     &        MA_ERR)
      end if
c
c     set SO ma handles to zero to keep compiler quier
c     all this stuff is used on SO only
c
      ld_ik5=0
      ld_jl5=0
      ld_il5=0
      ld_jk5=0
      ld_ik6=0
      ld_jl6=0
      ld_il6=0
      ld_jk6=0
      ld_ik7=0
      ld_jl7=0
      ld_il7=0
      ld_jk7=0
      ld_ik8=0
      ld_jl8=0
      ld_il8=0
      ld_jk8=0

C     define threshold for Schwarz screening(same as in SCF)
      if (.not. rtdb_get(rtdb, 'scf:tol2e', mt_dbl, 1, tol2e)) then
         if (rtdb_get(rtdb,'scf:thresh',mt_dbl,1,tol2e)) then
            tol2e = min(1d-7,tol2e * 1d-2)
         else
            tol2e = 1.0d-7
         end if
      end if
c     
c     Block the shells for Texas
c     
      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $     call errquit('tce_grad_force:could not allocate shmap',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $     call errquit('tce_grad_force:could not allocate blo',nsh,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $     call errquit('tce_grad_force:could not allocate bhi',nsh,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $     call errquit('tce_grad_force:could not allocate bflo',nsh,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $     call errquit('tce_grad_force:could not allocate bfhi',nsh,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $     call errquit('tce_grad_force:could not allocate blo',nbf,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $     call errquit('tce_grad_force:could not allocate bhi',nbf,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $     call errquit('tce_grad_force:could not allocate bfmap',nbf,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $     call errquit('tce_grad_force:could not allocate rbfmap',nbf,
     1     MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $     call errquit('tce_grad_force:could not allocate bftoat',nbf,
     1     MA_ERR)

      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap), 
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))
c-> debug
c      do i=1,nshblocks
c         print *,'debug bfglo=',i,int_mb(k_bfglo+i-1)
c         print *,'debug bfghi=',i,int_mb(k_bfghi+i-1)
c      enddo
c<- debug
      do i = 1, nbf
         if (.not. bas_bf2ce(basis, i, iat)) call errquit('bf2ce',i,
     &       BASIS_ERR)
         int_mb(k_bftoat+i-1) = iat
      end do
c
c     ==================================================================
c
c     BEWARE:  below here all AO indices have been put in texas order!
c
c     ==================================================================
c     
c     Now reorder all of the arrays with AO indices according to the map
c     
c     ==============
c     reorder movecs
c     ==============
c
      if (.not. ma_push_get(mt_dbl,nbf,'mo_tmp',l_mo_tmp,k_mo_tmp))
     1   call errquit('tce_grad_force: could not allocate mo_tmp',
     1   01,ma_err)
c
      do g1b=1,noab+nvab
         do g1=1,int_mb(k_range+g1b-1)
         do i=1,nbf
            dbl_mb(k_mo_tmp+int_mb(k_rbfmap+i-1)-1)=
     1      dbl_mb(k_movecs_sorted+(i-1)
     1      +(int_mb(k_offset+g1b-1)+g1-1)*nbf)
         enddo
         do i=1,nbf
            dbl_mb(k_movecs_sorted+(i-1)
     1      +(int_mb(k_offset+g1b-1)+g1-1)*nbf)=dbl_mb(k_mo_tmp+i-1)
         enddo
         enddo
      enddo
c      
      if (.not. ma_pop_stack(l_mo_tmp))
     1   call errquit('tce_grad_force: could not release l_mo_tmp',
     1   02, ma_err)
c
c      next=nxtask(-nprocs,1)  
c      call  ga_sync( )
c    ========================================
c    REMOVE HF REFERENCE CONTRIBUTION TO 1PDM
c    ========================================
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c      
      do g1b=1,noab
         if(next.eq.count) then
c
         if ((.not.restricted).or.(int_mb(k_spin+g1b-1).ne.2)) then
            dima=int_mb(k_range+g1b-1)*int_mb(k_range+g1b-1)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('ccsd_gradients: ma problem',2,ma_err)
            call get_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            do g1=1,int_mb(k_range+g1b-1)
               dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))=
     1         dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))-1.d0
            enddo
            call put_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            if (.not.ma_pop_stack(l_a))
     1      call errquit('ccsd_gradients: ma problem',3,ma_err)
         endif
         next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c     ======================================
c     allocate l_list and l_q4 for tce_grad2
c     ======================================
      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('tce_grad_force:could not allocate list',4*maxq,
     1     MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('tce_grad_force:could not allocate q4',maxq,
     1     MA_ERR)
c     ===============================================================
c     do the back transformation on d_2pdm with the reordered mo vecs
c     and the nblocked atomic tiles
c     ===============================================================
c      if(.not.ma_push_get(mt_dbl,nbf*nbf*nbf*nbf, '2pdm_ao',
c     1     l_2pdm_ao, k_2pdm_ao))
c     $     call errquit('tce_grad_force:could not allocate k_2pdm_ao',
c     1     ma_err)
c      do i=1,nbf*nbf*nbf*nbf
c         dbl_mb(k_2pdm_ao+i-1)=0.d0
c      enddo
      if(.not.ma_push_get(mt_int,nshblocks, 'dim_shb',
     1     l_dim_shb, k_dim_shb))
     $     call errquit('tce_grad_force:could not allocate dim_shb',
     1     ma_err)
c
      do i=1,nshblocks
         int_mb(k_dim_shb+i-1) =
     1   int_mb(k_bfghi+i-1)-int_mb(k_bfglo+i-1)+1
c         print *,'ishblock',i,int_mb(k_bfglo+i-1),int_mb(k_bfghi+i-1)
      enddo
c      call ga_sync( )
c
      call offset_btrans2_i0(l_2pdm_ao_offset,k_2pdm_ao_offset,
     1     size_2pdm_ao,nshblocks,int_mb(k_dim_shb))
      call tce_filename('2pdm_ao',filename)
      call createfile(filename,d_2pdm_ao,size_2pdm_ao)
c     =======================
c     DEBUG PRINT OUT 2PDM MO
c     =======================
c      call print_2pdm_mo(d_2pdm,k_2pdm_offset,size_2pdm)
c     ==========================
c     TRANSFORM 2PDM ON MO TO AO
c     ==========================
c
c      call btrans2(d_2pdm,k_2pdm_offset,d_1pdm,k_1pdm_offset,
c     1           k_2pdm_ao,atpart2,nalength2)
c
c      call ga_sync( )
c      next=nxtask(-nprocs,1)
      call btrans2(d_2pdm,k_2pdm_offset,d_1pdm,k_1pdm_offset,
     1             d_2pdm_ao,k_2pdm_ao_offset,size_2pdm_ao,
     1             nshblocks,int_mb(k_dim_shb),int_mb(k_bfglo),
     1             int_mb(k_bfghi))
c      call ga_sync( )
c      print *,'after btrans2'
c      call print_2pdm_ao(d_2pdm_ao,k_2pdm_ao_offset,nshblocks,
c     1     int_mb(k_bfglo),int_mb(k_bfghi))
c      call ga_sync( )

c      print *,'symmetrization'
c
      call sym_2pdm_ao(d_2pdm_ao,k_2pdm_ao_offset,nshblocks,
     1     int_mb(k_bfglo),int_mb(k_bfghi))
c      call ga_sync( )
c
c      if(.not.ma_pop_stack(l_dim_shb))
c     1  call errquit('tce_grad_force:',3,ma_err)      

c      call sym_2pdm_ao(d_2pdm_ao,k_2pdm_ao_offset,nshblocks,
c     1     int_mb(k_bfglo),int_mb(k_bfghi))
c     =========================================================
c     READ 2PDM AO FROM FILE
c     =========================================================
c     ======================
c     READ d_i0 to k_2pdm_ao
c     ======================
c      i_ind=0
c      do ib=1,nshblocks
c         idim = int_mb(k_dim_shb+ib-1)
c      j_ind=0
c      do jb=1,nshblocks
c         jdim = int_mb(k_dim_shb+jb-1)
c      k_ind=0
c      do kb=1,nshblocks
c         kdim = int_mb(k_dim_shb+kb-1)
c      l_ind=0
c      do lb=1,nshblocks
c         ldim = int_mb(k_dim_shb+lb-1)
c         dima=idim*jdim*kdim*ldim
c         if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
c     1   call errquit('btrans2: ma problem',3,ma_err)
c         call get_hash_block(d_2pdm_ao,dbl_mb(k_a),dima,
c     1        int_mb(k_2pdm_ao_offset),
c     1        (lb-1)+(kb-1)*nshblocks+(jb-1)*nshblocks*nshblocks+
c     1        (ib-1)*nshblocks*nshblocks*nshblocks)
c         do i=1,idim
c            i1=i_ind+i
c         do j=1,jdim
c            j1=j_ind+j
c         do k=1,kdim
c            k1=k_ind+k
c         do l=1,ldim
c            l1=l_ind+l
c            dbl_mb(k_2pdm_ao+(l1-1)+(j1-1)*nbf+(k1-1)*nbf*nbf
c     1      +(i1-1)*nbf*nbf*nbf) =
c     1      dbl_mb(k_2pdm_ao+(l1-1)+(j1-1)*nbf+(k1-1)*nbf*nbf
c     1      +(i1-1)*nbf*nbf*nbf) +
c     1      dbl_mb(k_a+(l-1)+(k-1)*ldim+(j-1)*ldim*kdim+
c     1      (i-1)*ldim*kdim*jdim)
c         enddo
c         enddo
c         enddo
c         enddo
c         if(.not.ma_pop_stack(l_a))
c     1   call errquit('btrans2: ma problem',4,ma_err)
c         l_ind=l_ind+ldim
c      enddo
c      k_ind=k_ind+kdim
c      enddo
c      j_ind=j_ind+jdim
c      enddo
c      i_ind=i_ind+idim
c      enddo
c     ====================================================================
c      do i=1,nbf
c      do j=1,nbf
c      do k=1,nbf
c      do l=1,nbf
c          write(6,'(a,4i5,a,f20.16)') 'FINAL 2PDM AO',i,j,k,l,'=',
c     1     dbl_mb(k_2pdm_ao+(l-1)+(k-1)*nbf+(j-1)*nbf*nbf
c     1       +(i-1)*nbf*nbf*nbf)
c      enddo
c      enddo
c      enddo
c      enddo
c
c      do i=1,nbf
c      do j=1,nbf
c      do k=1,nbf
c      do l=1,nbf
c         res=1.d0/8.d0*(
c     1       dbl_mb(k_2pdm_ao+(l-1)+(k-1)*nbf+(j-1)*nbf*nbf
c     1       +(i-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(l-1)+(k-1)*nbf+(i-1)*nbf*nbf
c     1       +(j-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(k-1)+(l-1)*nbf+(j-1)*nbf*nbf
c     1       +(i-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(k-1)+(l-1)*nbf+(i-1)*nbf*nbf
c     1       +(j-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(j-1)+(i-1)*nbf+(l-1)*nbf*nbf
c     1       +(k-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(i-1)+(j-1)*nbf+(l-1)*nbf*nbf
c     1       +(k-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(j-1)+(i-1)*nbf+(k-1)*nbf*nbf
c     1       +(l-1)*nbf*nbf*nbf)+
c     1       dbl_mb(k_2pdm_ao+(i-1)+(j-1)*nbf+(k-1)*nbf*nbf
c     1       +(l-1)*nbf*nbf*nbf))
c
c         dbl_mb(k_2pdm_ao+(l-1)+(k-1)*nbf+(j-1)*nbf*nbf
c     1       +(i-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(l-1)+(k-1)*nbf+(i-1)*nbf*nbf
c     1       +(j-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(k-1)+(l-1)*nbf+(j-1)*nbf*nbf
c     1       +(i-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(k-1)+(l-1)*nbf+(i-1)*nbf*nbf
c     1       +(j-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(j-1)+(i-1)*nbf+(l-1)*nbf*nbf
c     1       +(k-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(i-1)+(j-1)*nbf+(l-1)*nbf*nbf
c     1       +(k-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(j-1)+(i-1)*nbf+(k-1)*nbf*nbf
c     1       +(l-1)*nbf*nbf*nbf) = res
c         dbl_mb(k_2pdm_ao+(i-1)+(j-1)*nbf+(k-1)*nbf*nbf
c     1       +(l-1)*nbf*nbf*nbf) = res
c      enddo
c      enddo
c      enddo
c      enddo
c
c      do i=1,nbf
c      do j=1,nbf
c      do k=1,nbf
c      do l=1,nbf
c          write(6,'(a,4i5,a,f20.16)') 'FINAL 2PDM AO',i,j,k,l,'=',
c     1     dbl_mb(k_2pdm_ao+(l-1)+(k-1)*nbf+(j-1)*nbf*nbf
c     1       +(i-1)*nbf*nbf*nbf)
c      enddo
c      enddo
c      enddo
c      enddo
c
c     =========================================================
      
c      if(.not.ma_pop_stack(l_2pdm_ao))
c     1  call errquit('tce_grad_force:',3,ma_err)
c
c      call btrans2(d_2pdm,k_2pdm_offset,size_2pdm,
c     1             d_1pdm,k_1pdm_offset,size_1pdm,
c     1             d_2pdm_ao,k_2pdm_ao_offset,size_2pdm_ao,
c     1             nshblocks,
c     1             int_mb(k_bfglo),int_mb(k_bfghi))
c      print *,'BTRANS2 IS FINE!'
c     ==================
c     SYMMETRIZE 2PDM AO
c     ==================
c      call sym_2pdm_ao(d_2pdm_ao,k_2pdm_ao_offset,nshblocks,
c     1     int_mb(k_bfglo),int_mb(k_bfghi))
c     ==========================
c     DEBUG PRINT OUT 2PDM ON AO
c     ==========================
c      call ga_sync( )
c      call print_2pdm_ao(d_2pdm_ao,k_2pdm_ao_offset,nshblocks,
c     1     int_mb(k_bfglo),int_mb(k_bfghi))
c      call ga_sync( )
c     
c     
      cpu_tim(2)  = util_cpusec()
      wall_tim(2) = util_wallsec()
      if (jfac.ne.0.0d0 .or. kfac.ne.0.0d0) then
c->debug
c         do i=1,blen**4
c            dbl_mb(k_pdm2d+i-1)=0.d0
c         enddo
c<-debug
c         call tce_grad2(k_2pdm_ao,
c     $     dbl_mb(ld_ij), dbl_mb(ld_kl), dbl_mb(ld_ik),
c     $     dbl_mb(ld_jl), dbl_mb(ld_il), dbl_mb(ld_jk), 
c     $     dbl_mb(ld_ij2),dbl_mb(ld_kl2),dbl_mb(ld_ik2),
c     $     dbl_mb(ld_jl2),dbl_mb(ld_il2),dbl_mb(ld_jk2),
c     $     dbl_mb(ld_ij3),dbl_mb(ld_kl3),dbl_mb(ld_ik3),
c     $     dbl_mb(ld_jl3),dbl_mb(ld_il3),dbl_mb(ld_jk3),
c     $     dbl_mb(ld_ij4),dbl_mb(ld_kl4),dbl_mb(ld_ik4),
c     $     dbl_mb(ld_jl4),dbl_mb(ld_il4),dbl_mb(ld_jk4),
c     $     dbl_mb(ld_ik5),dbl_mb(ld_jl5),dbl_mb(ld_il5),
c     .     dbl_mb(ld_jk5),
c     $     dbl_mb(ld_ik6),dbl_mb(ld_jl6),dbl_mb(ld_il6),
c     .     dbl_mb(ld_jk6),
c     $     dbl_mb(ld_ik7),dbl_mb(ld_jl7),dbl_mb(ld_il7),
c     .     dbl_mb(ld_jk7),
c     $     dbl_mb(ld_ik8),dbl_mb(ld_jl8),dbl_mb(ld_il8),
c     .     dbl_mb(ld_jk8),
c     $     dbl_mb(k_frc_2el), g_dens, g_force, blen,
c     $     geom, basis, nproc, nat, 
c     $     lscratch, dbl_mb(k_scr), lbuf/12, dbl_mb(k_buf), 
c     $     int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
c     $     tol2e, nsh, 
c     $     log_mb(k_act), oskel, scftype, omp2, nopen, nbf,
c     $     dbl_mb(k_pdm2),dbl_mb(k_pdm2a),dbl_mb(k_pdm2b), ! MCSCF
c     $     dbl_mb(k_pdm2c),dbl_mb(k_pdm2d), dbl_mb(k_coeff), ! MCSCF
c     $     nshblocks,
c     $     int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi), 
c     $     int_mb(k_bfglo), int_mb(k_bfghi), 
c     $     int_mb(k_bfmap), int_mb(k_rbfmap), 
c     $     int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi),
c     $     jfac, kfac,.false.)
c         print *,'BEFORE TCE_GRAD2'
         call tce_grad2(d_2pdm_ao,k_2pdm_ao_offset,
     $     dbl_mb(ld_ij), dbl_mb(ld_kl), dbl_mb(ld_ik),
     $     dbl_mb(ld_jl), dbl_mb(ld_il), dbl_mb(ld_jk),
     $     dbl_mb(ld_ij2),dbl_mb(ld_kl2),dbl_mb(ld_ik2),
     $     dbl_mb(ld_jl2),dbl_mb(ld_il2),dbl_mb(ld_jk2),
     $     dbl_mb(ld_ij3),dbl_mb(ld_kl3),dbl_mb(ld_ik3),
     $     dbl_mb(ld_jl3),dbl_mb(ld_il3),dbl_mb(ld_jk3),
     $     dbl_mb(ld_ij4),dbl_mb(ld_kl4),dbl_mb(ld_ik4),
     $     dbl_mb(ld_jl4),dbl_mb(ld_il4),dbl_mb(ld_jk4),
     $     dbl_mb(ld_ik5),dbl_mb(ld_jl5),dbl_mb(ld_il5),
     .     dbl_mb(ld_jk5),
     $     dbl_mb(ld_ik6),dbl_mb(ld_jl6),dbl_mb(ld_il6),
     .     dbl_mb(ld_jk6),
     $     dbl_mb(ld_ik7),dbl_mb(ld_jl7),dbl_mb(ld_il7),
     .     dbl_mb(ld_jk7),
     $     dbl_mb(ld_ik8),dbl_mb(ld_jl8),dbl_mb(ld_il8),
     .     dbl_mb(ld_jk8),
     $     dbl_mb(k_frc_2el), g_1pdm, g_force, blen,
     $     geom, basis, nproc, nat,
     $     lscratch, dbl_mb(k_scr), lbuf/12, dbl_mb(k_buf),
     $     int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $     tol2e, nsh,
     $     log_mb(k_act), oskel, scftype, omp2, nopen, nbf,
     $     dbl_mb(k_pdm2),dbl_mb(k_pdm2a),dbl_mb(k_pdm2b), ! MCSCF
     $     dbl_mb(k_pdm2c),dbl_mb(k_pdm2d), dbl_mb(k_coeff), ! MCSCF
     $     nshblocks,
     $     int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi),
     $     int_mb(k_bfglo), int_mb(k_bfghi),
     $     int_mb(k_bfmap), int_mb(k_rbfmap),
     $     int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi),
     $     jfac, kfac,.false.)
c         print *,'AFTER TCE_GRAD2'
c
      endif
      cpu_tim(2)  = util_cpusec()  - cpu_tim(2)
      wall_tim(2) = util_wallsec() - wall_tim(2)
c
      call deletefile(d_2pdm_ao)
      if (.not.ma_pop_stack(l_2pdm_ao_offset))
     1   call errquit('tce_grad_force: ma problem',1,ma_err) 
      if (.not.ma_pop_stack(l_dim_shb))
     1   call errquit('tce_grad_force:',2,ma_err)
c      if(.not.ma_pop_stack(l_2pdm_ao))
c     1  call errquit('tce_grad_force:',3,ma_err)
c
c     
c     terminate integrals
c     
      call schwarz_tidy()
      call intd_terminate()

C     
      call ga_sync()
      call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_nuc), 3*nat, '+')
      call ga_dgop(msg_grad_wgh, dbl_mb(k_frc_wgh), 3*nat, '+')
      call ga_dgop(msg_grad_kin, dbl_mb(k_frc_kin), 3*nat, '+')
      call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_2el), 3*nat, '+')
      call ga_sync()
c     
      if (ga_nodeid() .eq. 0) then

         status = rtdb_parallel(.false.)
c->debug
         call sym_grad_symmetrize(geom, dbl_mb(k_frc_2el))
         call sym_grad_symmetrize(geom, dbl_mb(k_frc_nuc))
         call sym_grad_symmetrize(geom, dbl_mb(k_frc_wgh))
         call sym_grad_symmetrize(geom, dbl_mb(k_frc_kin))
c<-debug
         do i=0, 3*nat-1
            dbl_mb(k_force+i) = dbl_mb(k_force+i) +
     $           dbl_mb(k_frc_2el+i) + dbl_mb(k_frc_nuc+i) +
     $           dbl_mb(k_frc_wgh+i) + dbl_mb(k_frc_kin+i)
c->debug
c           write(6,'(a,i4,f20.16)') 'k_frc_2el',i,dbl_mb(k_frc_2el+i)
c           write(6,'(a,i4,f20.16)') 'k_frc_nuc',i,dbl_mb(k_frc_nuc+i)
c           write(6,'(a,i4,f20.16)') 'k_frc_wgh',i,dbl_mb(k_frc_wgh+i)
c           write(6,'(a,i4,f20.16)') 'k_frc_kin',i,dbl_mb(k_frc_kin+i)
c           write(6,'(a,i4,f20.16)') 'k_force', i, dbl_mb(k_force+i)
c<-debug
         end do

         if (odft) then
            if (.not. rtdb_get(rtdb, 'dft:cd+xc gradient', 
     $           mt_dbl, 3*nat, dbl_mb(k_frc_cd))) call errquit
     $           ('tce_grad_force: no dft cd+xc gradient',110, RTDB_ERR)
            do i=0, 3*nat-1
               dbl_mb(k_force+i) = dbl_mb(k_force+i) + 
     $              dbl_mb(k_frc_cd+i) 
            end do
         endif

         if (omp2) then
            if (.not. rtdb_get(rtdb, 'mp2:nonseparable gradient', 
     $           mt_dbl, 3*nat, dbl_mb(k_frc_mp2))) call errquit
     $           ('tce_grad_force: no nonseparable gradient',110,
     1           RTDB_ERR)
            do i=0, 3*nat-1
               dbl_mb(k_force+i) = dbl_mb(k_force+i) + 
     $              dbl_mb(k_frc_mp2+i) 
            end do
         end if

C     zero force contributions on inactive atoms
         call zero_forces(dbl_mb(k_force), log_mb(k_act), nat)
         
C     symmetrize
         if (oskel) then
            call sym_grad_symmetrize(geom, dbl_mb(k_force))
         end if

         if (.not. omp2) then
            if (scftype .eq. 'MCSCF') then
               theory = 'MCSCF'
            else if (omp2) then
               if (scftype .eq. 'UHF') then
                  theory = 'MP2(UHF)'
               else
                  theory = 'MP2(RHF)'
               end if
            else if (odft) then
               theory = 'DFT'
            else
               theory = scftype
            end if
         end if

         if (util_print('forces', print_high)) then
            write(luout,2200) 'nuclear repulsion gradient',' ',' ',
     $           ((dbl_mb(k_frc_nuc+i-1+3*(j-1)),i=1,3),j=1,nat)
            write(luout,2200) 'weighted density gradient',' ',' ',
     $           ((dbl_mb(k_frc_wgh+i-1+3*(j-1)),i=1,3),j=1,nat)
            write(luout,2200) 'kinetic energy gradient',' ',' ',
     $           ((dbl_mb(k_frc_kin+i-1+3*(j-1)),i=1,3),j=1,nat)
            write(luout,2200) '2-electron gradient',' ',' ',
     $           ((dbl_mb(k_frc_2el+i-1+3*(j-1)),i=1,3),j=1,nat)
            if (omp2) write(luout,2200) theory(1:inp_strlen(theory)),
     $           ' ', 'non-separable gradient',
     $           ((dbl_mb(k_frc_mp2+i-1+3*(j-1)),i=1,3),j=1,nat)
            if (odft) write(luout,2200) theory(1:inp_strlen(theory)),
     $           ' ', 'CD+XC gradient',
     $           ((dbl_mb(k_frc_cd+i-1+3*(j-1)),i=1,3),j=1,nat)
            write(luout,2200) 'total ',
     $           theory(1:inp_strlen(theory)), ' gradient',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)
 2200       format(A,A,A/,1000(3(1x,F12.6),/))
            call util_flush(luout)
         end if

C     store in rtdb
         if (omp2) then
            rtdb_string = 'mp2:gradient'
         else if (scftype.eq. 'MCSCF') then
            rtdb_string = 'mcscf:gradient'
         else if (odft) then
            rtdb_string = 'dft:gradient'
c->pdfan
         else if (otce) then
            rtdb_string = 'tce:gradient'
c<-pdfan
         else
            rtdb_string = 'scf:gradient'
         end if
         if (.not. rtdb_put(rtdb, rtdb_string, mt_dbl, 3*nat, 
     $        dbl_mb(k_force)))call errquit
     $        ('tce_grad_force: could not store gradients',1, RTDB_ERR)
c     
         status = rtdb_parallel(.true.)
c     
      end if
c     
      call ga_sync
c     

C     default: print the total forces 
      if (ga_nodeid() .eq. 0 
     $     .and. util_print('forces', print_low)) then
c         write(luout,1000) theory(1:inp_strlen(theory)),
c     $        'x','y','z','x','y','z'
         write(luout,1000) model,
     $        'x','y','z','x','y','z'
c         if (model.eq.'ccsd') then
c             write(luout,1000) 'CCSD',
c     $        'x','y','z','x','y','z'
c         else if(model.eq.'ccsdt') then
c             write(luout,1000) 'CCSDT',
c     $        'x','y','z','x','y','z'
c         else if(model.eq.'ccsdtq') then
c             write(luout,1000) 'CCSDTQ',
c     $        'x','y','z','x','y','z'
c         endif
c
         do 30, i=1, nat
            if (.not. geom_cent_get(geom, i, tag, crd, q)) call errquit
     $           ('tce_grad_force: geometry corrupt?',0, GEOM_ERR)
            write(luout,2000) i, tag,(crd(j),j=1,3),
     $           (dbl_mb(k_force+3*(i-1)+j),j=0,2)
   30    continue
         write(luout,*)
 1000    format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $        'coordinates',
     $        24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000    format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
         call util_flush(luout)
      end if

C     print timing information
      if (ga_nodeid().eq.0 .and.
     $     util_print('timing', print_default)) then
         write(luout,03000)cpu_tim,wall_tim
03000    format(17x,40('-'),/,
     &        17x,'|  Time  |  1-e(secs)   |  2-e(secs)   |',/,
     &        17x,40('-'),/,
     &        17x,'|  CPU   |',f11.2,3x,'|',f11.2,3x,'|',/,
     &        17x,40('-'),/,
     &        17x,'|  WALL  |',f11.2,3x,'|',f11.2,3x,'|',/,
     &        17x,40('-'))
         call util_flush(luout)
      end if
c     
      call ga_sync()

C     free memory

      if (.not. ma_verify_allocator_stuff())
     $     call errquit('tce_grad_force: ma corrupt',0, MA_ERR)

      if (.not. ma_chop_stack(l_force)) 
     $     call errquit('tce_grad_force: failed chopping MA stack',0,
     1     MA_ERR)
c     
      if (.not. ga_destroy(g_force)) 
     1   call errquit('tce_grad_force: GA?',0, GA_ERR)
c     
c      do i=1, ndens
c         if (.not. ga_destroy(g_dens(i))) call errquit 
c     $        ('tce_grad_force:error destroying density', 1, GA_ERR)
c      end do
       if (.not. ga_destroy(g_1pdm)) call errquit
     $    ('tce_grad_force:error destroying density', 1, GA_ERR)


c     
      end
c      double precision function grad_norm(frc, nat)
c      implicit none
c
c      integer nat
c      double precision frc(3*nat)
c
c      integer i
c      double precision g2
c
c      g2 = 0
c      do i=1, 3*nat
c        g2 = g2 + frc(i)*frc(i)
c      end do
c
c      grad_norm = sqrt(g2)
c
c      return
c      end
c      double precision function grad_max(frc, nat)
c      implicit none
c
c      integer nat
c      double precision frc(3,nat)
c
c      integer i
c      double precision norm, mxnorm
c
c      mxnorm = 0
c      do i=1, nat
c        norm = frc(1,i)*frc(1,i) + frc(2,i)*frc(2,i) + frc(3,i)*frc(3,i) 
c        mxnorm = max(norm,mxnorm)
c      end do
c
c      grad_max = sqrt(mxnorm)
c
c      return
c      end
c      subroutine grad_active_atoms(rtdb, natoms, oactive, nactive)
c      implicit none
c#include "errquit.fh"
c#include "rtdb.fh"
c#include "mafdecls.fh"
c#include "util.fh"
c#include "global.fh"
cc
c      integer rtdb              ! [input]
c      integer natoms            ! [input]
c      logical oactive(natoms)   ! [output]
c      integer nactive           ! [output]
cc
c      integer ma_type, l_actlist, k_actlist, i
cc
cc**** list of active atoms ****
c      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type, 
c     $     nactive, l_actlist)) then
c         if (.not. ma_get_index(l_actlist, k_actlist))
c     $        call errquit('grad_act_at: ma_get_inded failed',l_actlist,
c     &       MA_ERR)
c         if (nactive.le.0 .or. nactive.gt.natoms)
c     $        call errquit
c     $        ('grad_act_at: invalid number of active atoms',
c     $        natoms*10000+nactive, INPUT_ERR)
cC     fill lookup table
c         do i=1, natoms
c            oactive(i) = .false.
c         end do
c         do i=1, nactive
c            if (int_mb(k_actlist+i-1).le.0 .or. 
c     $           int_mb(k_actlist+i-1).gt.natoms) call errquit
c     $           ('grad_act_at: invalid active atom', INPUT_ERR,
c     $           int_mb(k_actlist+i-1))
c            oactive(int_mb(k_actlist+i-1)) = .true.
c         end do
c         if (.not. ma_free_heap(l_actlist)) call errquit
c     $        ('grad_act_at: free of actlist?',0, MA_ERR)
cc
cc**** list of inactive atoms ****
c      else if (rtdb_ma_get(rtdb, 'geometry:inactlist', ma_type, 
c     $     nactive, l_actlist)) then
c         if (.not. ma_get_index(l_actlist, k_actlist))
c     $        call errquit('grad_act_at: ma_get_inded failed',l_actlist,
c     &       MA_ERR)
c         if (nactive.le.0 .or. nactive.gt.natoms)
c     $        call errquit
c     $        ('grad_act_at: invalid number of inactive atoms',
c     $        natoms*10000+nactive, INPUT_ERR)
cC     fill lookup table
c         do i=1, natoms
c            oactive(i) = .true.
c         end do
c         do i=1, nactive
c            if (int_mb(k_actlist+i-1).le.0 .or. 
c     $           int_mb(k_actlist+i-1).gt.natoms) call errquit
c     $           ('grad_act_at: invalid active atom', INPUT_ERR,
c     $           int_mb(k_actlist+i-1))
c            oactive(int_mb(k_actlist+i-1)) = .false.
c         end do
cc
c         if (util_print('active atoms', print_low) .and. 
c     $       ga_nodeid().eq.0)   then
c           write(6,1) (int_mb(k_actlist+i-1),i=1,nactive)
c 1         format('  gradient active atoms ',15i4)
c*           write(6,2) (oactive(i),i=1,natoms)
c* 2         format('  oactive ',30l2)
c	 end if
cc
c         if (.not. ma_free_heap(l_actlist)) call errquit
c     $        ('grad_act_at: free of actlist?',0, MA_ERR)
cc
c      else
c         do i = 1, natoms
c            oactive(i) = .true.
c         end do
c         nactive = natoms
c      end if
cc
c      end
c
c      subroutine grad_shorder( basis, nsh, nbf, maxsh, blen,
c     $     nshblocks, shglo, shghi, shmap,
c     $     bfmap, rbfmap, bfglo, bfghi, shbflo, shbfhi)
c      implicit none
c#include "errquit.fh"
c#include "global.fh"
c#include "util.fh"
c#include "bas.fh"
c      integer basis             ! [input] basis handle
c      integer nsh, nbf          ! [input] no. of shells and basis functions
c      integer maxsh             ! [input]  Max. no. of shells in block
c      integer blen              ! [input]  Max no. of basis func in block
c      integer nshblocks         ! [output] Number of groups formed
c      integer shglo(*)          ! [output] First new shell in group
c      integer shghi(*)          ! [output] Last new shell in group
c      integer shmap(nsh)        ! [output] Shell new  --> orig
c      integer bfmap(nbf)        ! [output] Basis new  --> orig
c      integer rbfmap(nbf)       ! [output] Basis orig --> new
c      integer bfglo(*)          ! [output] First new basis func in group
c      integer bfghi(*)          ! [output] Last new basis func in group
c      integer shbflo(nsh)       ! [output] First new bf in new shell
c      integer shbfhi(nsh)       ! [output] Last new bf in new shell
cc     
c      logical oreorder          ! if true then actually reorder shells
c      logical oprint
cc
c      integer iuat, iatsh, itype, inprim, ingen, spch, iat, i, j
c      integer shlo, shhi, tmp, iold, bflo, bfhi, shdim, prevtype
c      integer shlen, bflen, newbf
c      double precision thetype
cc
c      data oreorder/.true./
cc
c      thetype(iuat,iatsh,itype,inprim,ingen) = 
c     $     iuat-1 + 32.0d0*(iatsh-1 + 32.0d0*(inprim-1 + 32.0d0*(
c     $     ingen-1 + 32.0d0*(itype+1))))
cc
c      oprint = util_print('reorder', print_debug) .and. ga_nodeid().eq.0
c      if (oreorder) then
cc
cc     Order shells by type (Texas 93 with suborder by uniq atom for Texas 95)
cc
c      do i = 1, nsh
c         if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
c     $        call errquit('fock_pairs: basis corrupt ',0,
c     &       BASIS_ERR)
c         if (.not. bas_cn2uce(basis,i,iuat))
c     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
c         if (.not. bas_cn2ce(basis,i,iat))
c     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
c         if (.not. bas_ce2cnr(basis,iat,shlo,shhi))
c     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
c         iatsh = i-shlo+1
cc     reorder L,s,p (-1,0,1) into order by complexity s,p,L (-1,0,1)
c         if (itype.eq.-1) then
c            itype = 1
c         else if (itype.eq.0) then
c            itype = -1
c         else if (itype.eq.1) then
c            itype = 0
c         end if
c         bfmap(i) = thetype(iuat,iatsh,itype,inprim,ingen) ! bfmap workspace
c         shmap(i) = i
c         do j = 1, i-1          ! Dumb nsh^2 sort into increaing order
c            if (bfmap(j).gt.bfmap(i)) then 
c               tmp = bfmap(j)
c               bfmap(j) = bfmap(i)
c               bfmap(i) = tmp
c               tmp = shmap(j)
c               shmap(j) = shmap(i)
c               shmap(i) = tmp
c            end if
c         end do
c      end do
c      else
c         do i = 1, nsh
c            shmap(i) = i
c            bfmap(i) = 1
c         end do
c      end if
cc
c      if (oprint) then
c         write(6,1) (shmap(i),i=1,nsh)
c 1       format('  shmap ', 12(i3,3x))
c         write(6,11)(bfmap(i),i=1,nsh)
c 11      format('  type  ', 12i9)
c      end if
cc
cc     Now have shmap (new to old shell map) and in bfmap the type of
cc     each new shell.  Perform the blocking subject to constraints 
cc     from blen and maxsh.
cc
c      iold = shmap(1)
c      if (.not. bas_cn2bfr(basis,iold,bflo,bfhi))
c     $     call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
c      shdim = bfhi - bflo + 1   ! Size of first new shell
cc
c      nshblocks = 1
c      shglo(1) = 1
c      bfglo(1) = 1
c      prevtype = bfmap(1)
c      shlen = 1
c      bflen = shdim
c      do i = 2, nsh
c         iold = shmap(i)
c         if (.not. bas_cn2bfr(basis,iold,bflo,bfhi))
c     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
c         shdim = bfhi - bflo + 1
c         itype = bfmap(i)
c         if (itype.ne.prevtype .or. shlen.eq.maxsh .or. 
c     $        (bflen+shdim).gt.blen) then
c            shghi(nshblocks) = shglo(nshblocks)+shlen-1
c            bfghi(nshblocks) = bfglo(nshblocks)+bflen-1
c            nshblocks = nshblocks + 1
c            shlen = 0
c            bflen = 0
c            shglo(nshblocks) = i
c            bfglo(nshblocks) = bfghi(nshblocks-1)+1
c            prevtype = itype
c         end if
c         shlen = shlen + 1
c         bflen = bflen + shdim
c      end do
c      shghi(nshblocks) = nsh
c      bfghi(nshblocks) = nbf
c      if (oprint) then
c         write(6,2) nshblocks, blen, maxsh
c 2       format(' nshblocks ', i3, '   blen ',i2,'  maxsh ',i2)
c         write(6,3) (shglo(i),i=1,nshblocks)
c 3       format('  shglo ', 12(i3,3x))
c         write(6,4) (shghi(i),i=1,nshblocks)
c 4       format('  shghi ', 12(i3,3x))
c         write(6,5) (bfglo(i),i=1,nshblocks)
c 5       format('  bfglo ', 12(i3,3x))
c         write(6,6) (bfghi(i),i=1,nshblocks)
c 6       format('  bfghi ', 12(i3,3x))
c      end if
cc
cc     Now have shmap, shglo/hi and bfglo/hi, and nshblocks
cc
cc     Just need to make the remaining maps
cc     bfmap(i) = new to old basis functions
cc     rbfmap(i) = old to new basis functions
cc
c      newbf = 0
c      do i = 1, nsh
c         iold = shmap(i)
c         if (.not. bas_cn2bfr(basis,iold,bflo,bfhi))
c     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
c         shbflo(i) = newbf + 1
c         do j = bflo,bfhi
c            newbf = newbf + 1
c            rbfmap(j) = newbf
c            bfmap(newbf) = j
c         end do
c         shbfhi(i) = newbf
c      end do
c      if (newbf .ne. nbf) call errquit('newbf is bad',0, BASIS_ERR)
cc
c      if (oprint) then
c         write(6,7) (bfmap(i),i=1,nbf)
c 7       format('  bfmap ', 12(i3,3x))
c         write(6,8) (rbfmap(i),i=1,nbf)
c 8       format(' rbfmap ', 12(i3,3x))
c         write(6,9) (shbflo(i),i=1,nsh)
c 9       format(' shbflo ', 12(i3,3x))
c         write(6,10) (shbfhi(i),i=1,nsh)
c 10      format(' shbfhi ', 12(i3,3x))
c      end if
cc
c      end
c      subroutine zero_forces ( force, oactive, nat )
c
c      double precision force, zero
c      logical oactive
c      integer nat
c
c      parameter ( zero = 0.D0 )
c
c      dimension force(3,nat)
c      dimension oactive(nat)
c
c      do i=1, nat
c        if ( .not.oactive(i) ) then
c          force(1,i) = zero
c          force(2,i) = zero
c          force(3,i) = zero
c        end if
c      end do
c
c      return
cc      end
