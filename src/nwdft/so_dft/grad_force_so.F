      subroutine grad_force_so(rtdb, basis, geom)
c     $Id$
C     calculate energy gradients with respect to nuclear coordinates
C------------------------------------------------------------------------------
C     ___                 ___                         ___
C     dE      \         dh(i,j)   \             d(mn|op)      \          dS(i,j)
C     -- = 2   > D(i,j) ------- +  > P(m,n,o,p) -------- - 2   > Dw(i,j) -------
C     dA      /           dA      /                dA         /            dA
C     ---                 ---                         ---
C     i,j                ijkl                         i,j
C     
C     
C     . dV(nuc-nuc)
C     + ----------
C     .    dA
C------------------------------------------------------------------------------
C     RHF:(H.Schlegel & S. Wolfe, JCP 63, p3632)
C     .                              1
C     P(i,j,k,l) = [2 D(i,j)D(k,l) - -(D(i,k)D(j,l) + D(i,l)D(j,k) ]
C     .                              2
C------------------------------------------------------------------------------
C     UHF:(M. Dupuis & H.F. King, JCP 68, p3998)
C     .                                1 
C     P(i,j,k,l) = [2 D+(i,j)D+(k,l) - -(D+(i,k)D+(j,l) + D-(i,k)D-(j,l) 
C     .                                2
C     + D+(i,l)D+(j,k) + D-(i,l)D-(j,k)) ]
C     D == D+
C------------------------------------------------------------------------------
C     ROHF:(Y. Yamaguchi, Y. Osamura, J.D. Goddard, H.F. Schaefer:
C     'A New Dimension to Quantum Chemistry', Oxford University Press '94, p74)
C     
C     -         I    J               I    J       I    J
C     P(ijkl) =  > [ a(IJ)D(ij)D(kl) + b(IJ)(D(ik)D(jl) + D(il)D(jk)) ]
C     -
C     IJ
C     I,J = open, closed
C------------------------------------------------------------------------------
C     UMP2:(M.J. Frisch, M. Head-Gordon, J.A. Pople, CPL 166, Nr. 3, p275,
C     explicit spin formalism(alpha/beta) from Robert
C     .           1
C     P(ijkl) = [ -((Da(ij)+Db(ij))*(Pa(kl)+Pb(kl))  
C     .           4
C     +(Pa(ij)+Pb(ij))*(Da(kl)+Db(kl)))
C     1
C     - -(Da(il)*Pa(jk) + Db(il)*Pb(jk) + Pa(il)*Da(jk) + Pb(il)*Db(jk)
C     8
C     + Da(jl)*Pa(ik) + Db(jl)*Pb(ik) + Pa(jl)*Da(ik) + Pb(jl)*Db(ik)
C     
C     ) ]
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
#include "xc.fh"
#include "case.fh"

      integer rtdb, basis, geom
      double precision cpu_tim(2), wall_tim(2)
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_dens(8),        ! density matrices(up to 8)
     $     g_wdens,             ! energy weighted density
     $     g_eigen_diag,        ! eigenvalue diagonal matrix or lagrangian
     $     g_force              ! forces on atoms(cartesian)
c
      character*255 movecs_in
c
      integer nproc, nat, nbf, nsh,max_sh_bf, max_at_bf, nopen,
     &     nclosed, 
     $     ndens, nbf2, nocc, max1e, max2e, mscratch_1e, mscratch_2e,
     $     lbuf, lscratch, lsqa, lsqatom, i, j, iat
c
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap, 
     $     l_bfglo, l_bfghi, l_labels, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap, 
     $     k_bfglo, k_bfghi, k_labels, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi
      integer k_shbflo, k_shbfhi
c      
      integer lforce, l_force, k_force, l_evals, k_evals, l_occ, k_occ,
     $     l_act, k_act, 
     $     l_buf, k_buf, l_scr, k_scr, l_dens, k_dens, l_wdens, k_wdens,
     $     k_frc_nuc, k_frc_kin, k_frc_2el, k_frc_2el_j, k_frc_2el_k, 
     $     k_frc_wgh, k_frc_cd,
     $     l_frc_nuc, l_frc_kin, l_frc_2el, l_frc_2el_j, l_frc_2el_k, 
     $     l_frc_wgh, l_frc_cd, l_frc_mp2, k_frc_mp2, l_frc_so, 
     $     k_frc_so,l_densx, k_densx, l_densy, k_densy, l_densz, k_densz

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
     $     lh_jl4, ld_jl4, lh_il4, ld_il4, lh_jk4, ld_jk4,
     $     lh_ik5, ld_ik5,
     $     lh_jl5, ld_jl5, lh_il5, ld_il5, lh_jk5, ld_jk5,
     $     lh_ik6, ld_ik6,
     $     lh_jl6, ld_jl6, lh_il6, ld_il6, lh_jk6, ld_jk6,
     $     lh_ik7, ld_ik7,
     $     lh_jl7, ld_jl7, lh_il7, ld_il7, lh_jk7, ld_jk7,
     $     lh_ik8, ld_ik8,
     $     lh_jl8, ld_jl8, lh_il8, ld_il8, lh_jk8, ld_jk8

      double precision crd(3),  ! atomic coordinates
     $     tol2e, q

      double precision grad_norm, grad_max
      external grad_norm, grad_max

      character*16 tag
      character*20 theory
      character*8 scftype
      character*32 rtdb_string

c    
c     Stuff for Douglas-Kroll
c
      logical num_grad, file_write_ga  ! Three data blocks needed
      external  num_grad, file_write_ga ! for the possible addition
c
*     Stuff for DFT

      integer ipol, noc(2)
      logical oskel, omp2, odft, ocdfit, status, exso 
      logical use_nwxc, lhashfxc
      double precision xfac(numfunc), jfac, kfac
c
      nproc = ga_nnodes()
c
      lh_ij2 = 0
      ld_ij2 = 0
      lh_kl2 = 0
      ld_kl2 = 0
      lh_ik2 = 0
      ld_ik2 = 0
      lh_jl2 = 0
      ld_jl2 = 0
      lh_il2 = 0
      ld_il2 = 0
      lh_jk2 = 0
      ld_jk2 = 0
      lh_ij3 = 0
      ld_ij3 = 0
      lh_kl3 = 0
      ld_kl3 = 0
      lh_ik3 = 0
      ld_ik3 = 0
      lh_jl3 = 0
      ld_jl3 = 0
      lh_il3 = 0
      ld_il3 = 0
      lh_jk3 = 0
      ld_jk3 = 0
      lh_ij4 = 0
      ld_ij4 = 0
      lh_kl4 = 0
      ld_kl4 = 0
      lh_ik4 = 0
      ld_ik4 = 0
      lh_jl4 = 0
      ld_jl4 = 0
      lh_il4 = 0
      ld_il4 = 0
      lh_jk4 = 0
      ld_jk4 = 0
      lh_ik5 = 0
      ld_ik5 = 0
      lh_jl5 = 0
      ld_jl5 = 0
      lh_il5 = 0
      ld_il5 = 0
      lh_jk5 = 0
      ld_jk5 = 0
      lh_ik6 = 0
      ld_ik6 = 0
      lh_jl6 = 0
      ld_jl6 = 0
      lh_il6 = 0
      ld_il6 = 0
      lh_jk6 = 0
      ld_jk6 = 0
      lh_ik7 = 0
      ld_ik7 = 0
      lh_jl7 = 0
      ld_jl7 = 0
      lh_il7 = 0
      ld_il7 = 0
      lh_jk7 = 0
      ld_jk7 = 0
      lh_ik8 = 0
      ld_ik8 = 0
      lh_jl8 = 0
      ld_jl8 = 0
      lh_il8 = 0
      ld_il8 = 0
      lh_jk8 = 0
      ld_jk8 = 0
c
      jfac = 1.0d0              ! Only changed for DFT
      kfac = 1.0d0
      exso=.false.

C     get information about basis set
      if (.not. geom_ncent(geom,nat))
     $     call errquit('grad: could not get natoms',0, GEOM_ERR)
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('grad: could not get nbf',0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('grad: could not get nsh',0, BASIS_ERR)
c     
c     Atom blocking now only used for the 1-e integrals (history)
c     
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('grad: could not get max_at_bf',0, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('grad: could not get max_sh_bf',0, BASIS_ERR)

      blen = min(nbf,maxblen,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets

      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('grad: failed getting theory',0, RTDB_ERR)

C     get SCF MO vectors for density

      odft = .false.
      if (theory .eq. 'mcscf') then
         if (.not.rtdb_cget(rtdb, 'mcscf:input vectors',1,movecs_in))
     $        call errquit('gradients: MCSCF MO vectors not defined',0,
     &       RTDB_ERR)
         scftype = 'mcscf'
      else if (theory .eq. 'dft') then
         odft = .true.
         if (.not. rtdb_cget(rtdb, 'dft:input vectors',1,movecs_in))
     $        call errquit('gradients: DFT MO vectors not defined',0,
     &       RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'dft:cdfit', mt_log, 1, ocdfit))
     $        ocdfit = .false.
         if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     $        call errquit('gradients: xfac not in rtdb ', 0,
     &       RTDB_ERR)
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
      else if (theory .eq. 'sodft') then
         odft = .true.
         if (.not. rtdb_cget(rtdb, 'dft:input vectors',1,movecs_in))
     $        call errquit('gradients: DFT MO vectors not defined',0,
     &       RTDB_ERR)
         call util_file_name_resolve(movecs_in,.false.)
         if (.not. rtdb_get(rtdb, 'dft:cdfit', mt_log, 1, ocdfit))
     $        ocdfit = .false.
         if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     $        call errquit('gradients: xfac not in rtdb ', 0,
     &       RTDB_ERR)
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
     $        call errquit('gradients: scftype not defined',0,
     &       RTDB_ERR)
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('gradients: SCF MO vectors not defined',0,
     &       RTDB_ERR)
      end if
      use_nwxc = util_module_avail("nwxc")
      if (use_nwxc) then
        call nwxc_rtdb_load(rtdb,"dft",use_nwxc)
      endif
      if (use_nwxc) then
        call nwxc_getvals("nwxc_is_on",use_nwxc)
      endif
      if (use_nwxc) then
        call nwxc_getvals("nwxc_has_hfx",lhashfxc)
        if (lhashfxc) then
          call nwxc_getwght("nwxc_wght_hfx",xfac(1))
          kfac = xfac(1)
        endif
      endif

C     get information about type of calculation

      if ((theory .eq. 'mp2') 
     $     .or.(theory .eq. 'semi_dir_mp2')
     $     .or.(theory .eq. 'direct_mp2')) then
         omp2 = .true.
      else
         omp2 = .false.
      end if

C     scftype: MCSCF, RHF, ROHF or UHF

      if (scftype .eq. 'mcscf') then
         if (.not. rtdb_get(rtdb, 'mcscf:nclosed',mt_int,1, nclosed))
     $        nclosed = 0
         if (.not. rtdb_get(rtdb, 'mcscf:nact', mt_int, 1, nopen))
     $        call errquit('gradients: no mcscf active orbitals?',0,
     &       RTDB_ERR)
      else if (odft) then
         if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &        call errquit('gradients: rtdb_get of noc failed', 0,
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
     $        call errquit('gradients: get of nclosed? ', 0,
     &       RTDB_ERR)
         
C     is it MP2?
         if (.not. omp2) then
            if (ga_nodeid() .eq. 0) then
               if (util_print('information', print_medium)) then
                  if (odft) then
                     write(luout,101) ' dft ', ipol
 101                 format(/'  wavefunction    =   ', a, ' ipol =',i2/)
                  else
                     write(luout,1) scftype
 1                   format(/'  wavefunction    =   ', a/)
                  endif
                  call util_flush(luout)
               end if
            end if
         end if
         
      end if
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
            call errquit('grad;no ROMP2 gradients yet', 110,
     &       CAPMIS_ERR)
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
      if (.not. ga_create(mt_dbl, 3, nat, 'forces', 3, 0, g_force))
     $     call errquit('gradients: failed to create force GA',0,
     &       GA_ERR)
      call ga_zero(g_force)

C     local replication(separate for the different pieces)
      lforce = nat * 3
      if (.not. ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     $     call errquit('grad:could not allocate l_force',lforce,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_nuc, 
     $     k_frc_nuc)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_kin, 
     $     k_frc_kin)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_wgh, 
     $     k_frc_wgh)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_so, 
     $     k_frc_so)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
c
      if (.not.cam_exch) then
       if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_2el, 
     $     k_frc_2el)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      else
       if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_2el_j, 
     $     k_frc_2el_j)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
       if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_2el_k, 
     $     k_frc_2el_k)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      end if ! cam_exch
c
      if (omp2) then
         if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_mp2, 
     $        k_frc_mp2)) call errquit
     &      ('grad:could not allocate l_force',lforce, MA_ERR)
      end if
      if (odft) then
         if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_cd, 
     $        k_frc_cd)) call errquit
     &      ('grad:could not allocate l_force',lforce, MA_ERR)
      end if

C     global density
cso   
      ndens = 8                 !alpha, beta, x, y and z densities
                                !Re(Daa), Re(Dbb), 
                                !Im(Daa), Im(Dbb), Re(Dab), Im(Dba)
      do i=1, ndens
         status = ga_create(mt_dbl,nbf,nbf,'density matrix', 0, 0, 
     &        g_dens(i))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' grad_force_so: ', 'density matrix'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('grad_force_so: ga_create ', 0, GA_ERR)
         endif
      end do
c     
      if (scftype .eq. 'MCSCF') then
         if (.not. ma_push_get(mt_dbl, nopen**4, 'pdm2', 
     $        l_pdm2 , k_pdm2 )) call errquit
     $        ('gradient:failed allocating pdm2',nopen**4, MA_ERR)
         if (.not. ma_push_get(mt_dbl, (blen**1)*(nopen**3), 'pdm2a',
     $        l_pdm2a, k_pdm2a)) call errquit
     $        ('gradient:failed allocating pdm2a',(blen**1)*(nopen**3),
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, (blen**2)*(nopen**2), 'pdm2b',
     $        l_pdm2b, k_pdm2b)) call errquit
     $        ('gradient:failed allocating pdm2b',(blen**2)*(nopen**2),
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, (blen**3)*(nopen**1), 'pdm2c',
     $        l_pdm2c, k_pdm2c)) call errquit
     $        ('gradient:failed allocating pdm2c',(blen**3)*(nopen**1),
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, nopen*nbf, 'coeff',
     $        l_coeff, k_coeff)) call errquit
     $        ('gradient:failed allocating coeff',nopen*nbf, MA_ERR)
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
     $     ('gradient:failed allocating pdm2d',blen**4, MA_ERR)

C     lookup table and list of active atoms
      if (.not. ma_push_get(MT_LOG,nat,'active atoms',l_act,k_act))
     $     call errquit('grad: could not allocate l_act',nat, MA_ERR)

C     symmetry or not
      if (.not. rtdb_get(rtdb, 'gradients:use symmetry',mt_log,1,
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
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')

C     eigenvalue diagonal matrix for forming energy weighted density
cso      nocc = nopen + nclosed
      nocc = nopen + 2*nclosed
      if (.not. ga_create(mt_dbl, nocc, nocc,'eigen_diag', 1, nocc, 
     $     g_eigen_diag)) call errquit('gradients: ga diags?',nocc*nocc,
     &       GA_ERR)
      call ga_zero(g_eigen_diag)

C     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('grad: could not allocate l_evals',nbf2, MA_ERR)
C     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf2,'occ. numbers', l_occ, k_occ))
     $     call errquit('grad: could not allocate l_occ',nbf2, MA_ERR)
c     
      call grad_active_atoms(rtdb, nat, log_mb(k_act), nactive)
c
c     
c     call grad_dens_so to get the total (weighted) density
c     
      call grad_dens_so(g_dens, g_wdens, g_eigen_diag, 
     $     dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2, 
     $     nopen, nclosed, scftype, movecs_in)

C     free temporary arrays
      if (.not. ga_destroy(g_eigen_diag))
     $     call errquit('grad: could not destroy g_eigen_diag',1,
     &       GA_ERR)

      if (.not.ma_pop_stack(l_occ))
     $     call errquit('grad:ma free occ',1, MA_ERR)
      if (.not.ma_pop_stack(l_evals))
     $     call errquit('grad:ma free eval',1, MA_ERR)

C     initialize for integral gradients
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
      call int_terminate()
      call intd_init(rtdb, 1, basis)
      call dfill(lforce, 0.0D0, dbl_mb(k_force), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_nuc), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_kin), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_wgh), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_so), 1)
      if (.not.cam_exch) then
         call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el), 1)
      else
         call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el_j), 1)
         call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el_k), 1)
      end if ! cam_exch
c
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e) ! blocking algorithm
      max2e = max2e/12
      max1e = max1e/12
      max2e = max(max2e,1296*100)          ! 100 D quartets 
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)

C     one-electron contribution
C     buffers for one electron integral derivatives

      if (.not. ma_push_get(mt_dbl,12*lbuf,'deriv buffer',l_buf,k_buf))
     $     call errquit('grad:could not allocate buffer',12*lbuf,
     &       MA_ERR)
*      call ma_summarize_allocated_blocks()
      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch', 
     $     l_scr, k_scr))call errquit('grad: scratch alloc failed',
     $     lscratch, MA_ERR)
      if (.not. ma_push_get(mt_int,4*lbuf,'labels',l_labels,k_labels))
     $     call errquit('grad: could not allocate labels',4*lbuf,
     &       MA_ERR)
c     
C     local density matrix block 
      lsqatom = max_at_bf * max_at_bf
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_dens,k_dens))
     $     call errquit('grad:could not allocate l_dens',lsqatom,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_densx,k_densx))
     $     call errquit('grad:could not allocate l_dens',lsqatom,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_densy,k_densy))
     $     call errquit('grad:could not allocate l_dens',lsqatom,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_densz,k_densz))
     $     call errquit('grad:could not allocate l_dens',lsqatom,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_w_density',l_wdens, 
     $     k_wdens)) call errquit
     &    ('grad;could not allocate l_wdens',lsqatom, MA_ERR)
c     
      cpu_tim(1)  = util_cpusec()
      wall_tim(1) = util_wallsec()
c      goto 111 !skip 1-e contribution 
      call grad1_so( dbl_mb(k_buf), lbuf, dbl_mb(k_scr), lscratch,
     $     dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
     $     dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh),
     $     g_dens, g_wdens, basis, geom, nproc, nat, max_at_bf,
     $     oskel, 
     &     dbl_mb(k_frc_so), 
     &     dbl_mb(k_densx), dbl_mb(k_densy), dbl_mb(k_densz)) 
 111  continue
      cpu_tim(1)  = util_cpusec() - cpu_tim(1)
      wall_tim(1) = util_wallsec() - wall_tim(1)
c     
      if (.not. ma_pop_stack(l_wdens)) call errquit('grad: MA?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_densz)) call errquit('grad: MA?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_densy)) call errquit('grad: MA?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_densx)) call errquit('grad: MA?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_dens)) call errquit('grad: MA?',0,
     &       MA_ERR)
      if (.not. ga_destroy(g_wdens)) call errquit('grad: GA?',0,
     &       GA_ERR)
C
C     two-electron contribution
C
C     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('grad: could not allocate l_evals',nbf2, MA_ERR)
C     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf2,'occ. numbers', l_occ, k_occ))
     $     call errquit('grad: could not allocate l_occ',nbf2, MA_ERR)
c
      if(theory .eq. 'sodft' .and. xfac(1) .ne. 0)
     &     call grad_dens_so2(g_dens,  
     $     dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2, 
     $     nopen, nclosed, movecs_in)
c
      if (.not.ma_pop_stack(l_occ))
     $     call errquit('grad:ma free occ',1, MA_ERR)
      if (.not.ma_pop_stack(l_evals))
     $     call errquit('grad:ma free eval',1, MA_ERR)
c
C     allocate arrays for two-electron integral stuff(or rename old ones)
C     blocks of density matrix:
c
      lsqa = blen*blen

      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij', lh_ij, ld_ij))
     $     call errquit('grad:could not allocate ld_ij',lsqa, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl', lh_kl, ld_kl))
     $     call errquit('grad:could not allocate ld_kl',lsqa, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik', lh_ik, ld_ik))
     $     call errquit('grad:could not allocate ld_ik',lsqa, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl', lh_jl, ld_jl))
     $     call errquit('grad:could not allocate ld_jl',lsqa, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il', lh_il, ld_il))
     $     call errquit('grad:could not allocate ld_il',lsqa, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk', lh_jk, ld_jk))
     $     call errquit('grad:could not allocate ld_jk',lsqa, MA_ERR)
c     
      if (scftype .eq. 'UHF' .or. 
     $     scftype .eq. 'ROHF' .or. omp2) then ! UHF or ROHF or MP2
         if (scftype .eq. 'ROHF' .or. omp2) then
            if (.not. ma_push_get(mt_dbl,lsqa,'ld_ij2',lh_ij2,ld_ij2))
     $           call errquit('grad:could not allocate ld_ij2',lsqa,
     &       MA_ERR)
            if (.not. ma_push_get(mt_dbl,lsqa,'ld_kl2',lh_kl2,ld_kl2))
     $           call errquit('grad:could not allocate ld_kl2',lsqa,
     &       MA_ERR)
         end if          
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik2', lh_ik2, ld_ik2))
     $        call errquit('grad:could not allocate ld_ik2',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl2', lh_jl2, ld_jl2))
     $        call errquit('grad:could not allocate ld_jl2',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il2', lh_il2, ld_il2))
     $        call errquit('grad:could not allocate ld_il2',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk2', lh_jk2, ld_jk2))
     $        call errquit('grad:could not allocate ld_jk2',lsqa,
     &       MA_ERR)
      end if

      if (scftype .eq. 'UHF' .and. omp2) then
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij3', lh_ij3, ld_ij3))
     $        call errquit('grad:could not allocate ld_ij3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl3', lh_kl3, ld_kl3))
     $        call errquit('grad:could not allocate ld_jl3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik3', lh_ik3, ld_ik3))
     $        call errquit('grad:could not allocate ld_ik3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl3', lh_jl3, ld_jl3))
     $        call errquit('grad:could not allocate ld_jl3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il3', lh_il3, ld_il3))
     $        call errquit('grad:could not allocate ld_il3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk3', lh_jk3, ld_jk3))
     $        call errquit('grad:could not allocate ld_jk3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ij4', lh_ij4, ld_ij4))
     $        call errquit('grad:could not allocate ld_ij4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_kl4', lh_kl4, ld_kl4))
     $        call errquit('grad:could not allocate ld_kl4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik4', lh_ik4, ld_ik4))
     $        call errquit('grad:could not allocate ld_ik4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl4', lh_jl4, ld_jl4))
     $        call errquit('grad:could not allocate ld_jl4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il4', lh_il4, ld_il4))
     $        call errquit('grad:could not allocate ld_il4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk4', lh_jk4, ld_jk4))
     $        call errquit('grad:could not allocate ld_jk4',lsqa,
     &       MA_ERR)
      end if

      if (theory .eq. 'sodft' .and. kfac .ne. 0) then
         exso = .true.
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik3', lh_ik3, ld_ik3))
     $        call errquit('grad:could not allocate ld_ik3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl3', lh_jl3, ld_jl3))
     $        call errquit('grad:could not allocate ld_jl3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il3', lh_il3, ld_il3))
     $        call errquit('grad:could not allocate ld_il3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk3', lh_jk3, ld_jk3))
     $        call errquit('grad:could not allocate ld_jk3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik4', lh_ik4, ld_ik4))
     $        call errquit('grad:could not allocate ld_ik4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl4', lh_jl4, ld_jl4))
     $        call errquit('grad:could not allocate ld_jl4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il4', lh_il4, ld_il4))
     $        call errquit('grad:could not allocate ld_il4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk4', lh_jk4, ld_jk4))
     $        call errquit('grad:could not allocate ld_jk4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik5', lh_ik5, ld_ik5))
     $        call errquit('grad:could not allocate ld_ik3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl5', lh_jl5, ld_jl5))
     $        call errquit('grad:could not allocate ld_jl3',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il5', lh_il5, ld_il5))
     $        call errquit('grad:could not allocate ld_il5',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk5', lh_jk5, ld_jk5))
     $        call errquit('grad:could not allocate ld_jk5',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik6', lh_ik6, ld_ik6))
     $        call errquit('grad:could not allocate ld_ik4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl6', lh_jl6, ld_jl6))
     $        call errquit('grad:could not allocate ld_jl4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il6', lh_il6, ld_il6))
     $        call errquit('grad:could not allocate ld_il4',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk6', lh_jk6, ld_jk6))
     $        call errquit('grad:could not allocate ld_jk6',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik7', lh_ik7, ld_ik7))
     $        call errquit('grad:could not allocate ld_ik7',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl7', lh_jl7, ld_jl7))
     $        call errquit('grad:could not allocate ld_jl7',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il7', lh_il7, ld_il7))
     $        call errquit('grad:could not allocate ld_il7',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk7', lh_jk7, ld_jk7))
     $        call errquit('grad:could not allocate ld_jk7',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_ik8', lh_ik8, ld_ik8))
     $        call errquit('grad:could not allocate ld_ik8',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jl8', lh_jl8, ld_jl8))
     $        call errquit('grad:could not allocate ld_jl8',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_il8', lh_il8, ld_il8))
     $        call errquit('grad:could not allocate ld_il8',lsqa,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, lsqa, 'ld_jk8', lh_jk8, ld_jk8))
     $        call errquit('grad:could not allocate ld_jk8',lsqa,
     &       MA_ERR)
      end if

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
     $     call errquit('grad:could not allocate shmap',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $     call errquit('grad:could not allocate blo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $     call errquit('grad:could not allocate bhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $     call errquit('grad:could not allocate bflo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $     call errquit('grad:could not allocate bfhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $     call errquit('grad:could not allocate blo',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $     call errquit('grad:could not allocate bhi',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $     call errquit('grad:could not allocate bfmap',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $     call errquit('grad:could not allocate rbfmap',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $     call errquit('grad:could not allocate bftoat',nbf, MA_ERR)

      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap), 
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))

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
      do i = 1, ndens
         call ga_reorder(g_dens(i), .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
      end do
c
      if (scftype .eq. 'MCSCF') then 
         call matrix_reorder(nbf, nopen,
     $        dbl_mb(k_coeff), .true., int_mb(k_rbfmap),
     $        .false., int_mb(k_rbfmap))
      end if
c

      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('grad:could not allocate list',4*maxq, MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('grad:could not allocate q4',maxq, MA_ERR)
c
c     
      cpu_tim(2)  = util_cpusec()
      wall_tim(2) = util_wallsec()
      if (jfac.ne.0.0d0 .or. kfac.ne.0.0d0) then
       if (.not.cam_exch) then
         call grad2(
     $     dbl_mb(ld_ij), dbl_mb(ld_kl), dbl_mb(ld_ik),
     $     dbl_mb(ld_jl), dbl_mb(ld_il), dbl_mb(ld_jk), 
     $     dbl_mb(ld_ij2),dbl_mb(ld_kl2),dbl_mb(ld_ik2),
     $     dbl_mb(ld_jl2),dbl_mb(ld_il2),dbl_mb(ld_jk2),
     $     dbl_mb(ld_ij3),dbl_mb(ld_kl3),dbl_mb(ld_ik3),
     $     dbl_mb(ld_jl3),dbl_mb(ld_il3),dbl_mb(ld_jk3),
     $     dbl_mb(ld_ij4),dbl_mb(ld_kl4),dbl_mb(ld_ik4),
     $     dbl_mb(ld_jl4),dbl_mb(ld_il4),dbl_mb(ld_jk4),
     $     dbl_mb(ld_ik5),
     $     dbl_mb(ld_jl5),dbl_mb(ld_il5),dbl_mb(ld_jk5),     
     $     dbl_mb(ld_ik6),
     $     dbl_mb(ld_jl6),dbl_mb(ld_il6),dbl_mb(ld_jk6),     
     $     dbl_mb(ld_ik7),
     $     dbl_mb(ld_jl7),dbl_mb(ld_il7),dbl_mb(ld_jk7),     
     $     dbl_mb(ld_ik8),
     $     dbl_mb(ld_jl8),dbl_mb(ld_il8),dbl_mb(ld_jk8),     
     $     dbl_mb(k_frc_2el), g_dens, g_force, blen,
     $     geom, basis, nproc, nat, 
     $     lscratch, dbl_mb(k_scr), lbuf, dbl_mb(k_buf), 
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
     $     jfac, kfac, exso)
      else
         call case_setflags(.false.)  ! for the J part
         call grad2(
     $     dbl_mb(ld_ij), dbl_mb(ld_kl), dbl_mb(ld_ik),
     $     dbl_mb(ld_jl), dbl_mb(ld_il), dbl_mb(ld_jk), 
     $     dbl_mb(ld_ij2),dbl_mb(ld_kl2),dbl_mb(ld_ik2),
     $     dbl_mb(ld_jl2),dbl_mb(ld_il2),dbl_mb(ld_jk2),
     $     dbl_mb(ld_ij3),dbl_mb(ld_kl3),dbl_mb(ld_ik3),
     $     dbl_mb(ld_jl3),dbl_mb(ld_il3),dbl_mb(ld_jk3),
     $     dbl_mb(ld_ij4),dbl_mb(ld_kl4),dbl_mb(ld_ik4),
     $     dbl_mb(ld_jl4),dbl_mb(ld_il4),dbl_mb(ld_jk4),
     $     dbl_mb(ld_ik5),
     $     dbl_mb(ld_jl5),dbl_mb(ld_il5),dbl_mb(ld_jk5),     
     $     dbl_mb(ld_ik6),
     $     dbl_mb(ld_jl6),dbl_mb(ld_il6),dbl_mb(ld_jk6),     
     $     dbl_mb(ld_ik7),
     $     dbl_mb(ld_jl7),dbl_mb(ld_il7),dbl_mb(ld_jk7),     
     $     dbl_mb(ld_ik8),
     $     dbl_mb(ld_jl8),dbl_mb(ld_il8),dbl_mb(ld_jk8),     
     $     dbl_mb(k_frc_2el_j), g_dens, g_force, blen,
     $     geom, basis, nproc, nat, 
     $     lscratch, dbl_mb(k_scr), lbuf, dbl_mb(k_buf), 
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
     $     jfac, 0.d0, .false.)
c
         call case_setflags(.true.)  ! for the K part
         call grad2(
     $     dbl_mb(ld_ij), dbl_mb(ld_kl), dbl_mb(ld_ik),
     $     dbl_mb(ld_jl), dbl_mb(ld_il), dbl_mb(ld_jk), 
     $     dbl_mb(ld_ij2),dbl_mb(ld_kl2),dbl_mb(ld_ik2),
     $     dbl_mb(ld_jl2),dbl_mb(ld_il2),dbl_mb(ld_jk2),
     $     dbl_mb(ld_ij3),dbl_mb(ld_kl3),dbl_mb(ld_ik3),
     $     dbl_mb(ld_jl3),dbl_mb(ld_il3),dbl_mb(ld_jk3),
     $     dbl_mb(ld_ij4),dbl_mb(ld_kl4),dbl_mb(ld_ik4),
     $     dbl_mb(ld_jl4),dbl_mb(ld_il4),dbl_mb(ld_jk4),
     $     dbl_mb(ld_ik5),
     $     dbl_mb(ld_jl5),dbl_mb(ld_il5),dbl_mb(ld_jk5),     
     $     dbl_mb(ld_ik6),
     $     dbl_mb(ld_jl6),dbl_mb(ld_il6),dbl_mb(ld_jk6),     
     $     dbl_mb(ld_ik7),
     $     dbl_mb(ld_jl7),dbl_mb(ld_il7),dbl_mb(ld_jk7),     
     $     dbl_mb(ld_ik8),
     $     dbl_mb(ld_jl8),dbl_mb(ld_il8),dbl_mb(ld_jk8),     
     $     dbl_mb(k_frc_2el_k), g_dens, g_force, blen,
     $     geom, basis, nproc, nat, 
     $     lscratch, dbl_mb(k_scr), lbuf, dbl_mb(k_buf), 
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
     $     0.d0, kfac, exso)
         call case_setflags(.false.)
       end if ! cam_exch
      end if
c
      cpu_tim(2)  = util_cpusec()  - cpu_tim(2)
      wall_tim(2) = util_wallsec() - wall_tim(2)
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
      call ga_dgop(msg_grad_so, dbl_mb(k_frc_so), 3*nat, '+')
      if (.not.cam_exch) then
        call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_2el), 3*nat, '+')
      else
        call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_2el_j), 3*nat, '+')
        call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_2el_k), 3*nat, '+')
      end if ! cam_exch
c
      call ga_sync()
c     
      if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
         do i=0, 3*nat-1
          if (.not.cam_exch) then
            dbl_mb(k_force+i) = 
     $           dbl_mb(k_frc_2el+i) + dbl_mb(k_frc_nuc+i) +
     $           dbl_mb(k_frc_wgh+i) + dbl_mb(k_frc_kin+i) + 
     $           dbl_mb(k_frc_so+i)
          else
            dbl_mb(k_force+i) = 
     $           dbl_mb(k_frc_2el_j+i) + dbl_mb(k_frc_2el_k+i) +
     $           dbl_mb(k_frc_nuc+i) +
     $           dbl_mb(k_frc_wgh+i) + dbl_mb(k_frc_kin+i) + 
     $           dbl_mb(k_frc_so+i)
          end if ! cam_exch
         end do

         if (odft) then
            if (.not. rtdb_get(rtdb, 'dft:cd+xc gradient', 
     $           mt_dbl, 3*nat, dbl_mb(k_frc_cd))) call errquit
     $           ('grad_force: no dft cd+xc gradient',110, RTDB_ERR)
            do i=0, 3*nat-1
               dbl_mb(k_force+i) = dbl_mb(k_force+i) + 
     $              dbl_mb(k_frc_cd+i) 
            end do
         endif

         if (omp2) then
            if (.not. rtdb_get(rtdb, 'mp2:nonseparable gradient', 
     $           mt_dbl, 3*nat, dbl_mb(k_frc_mp2))) call errquit
     $           ('grad_force: no nonseparable gradient',110, RTDB_ERR)
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
               theory = theory
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
c
            if (.not.cam_exch) then
              write(luout,2200) '2-electron gradient',' ',' ',
     $           ((dbl_mb(k_frc_2el+i-1+3*(j-1)),i=1,3),j=1,nat)
            else
              write(luout,2200) '2-electron gradient-J',' ',' ',
     $           ((dbl_mb(k_frc_2el_j+i-1+3*(j-1)),i=1,3),j=1,nat)
              write(luout,2200) '2-electron gradient-K',' ',' ',
     $           ((dbl_mb(k_frc_2el_k+i-1+3*(j-1)),i=1,3),j=1,nat)
            end if ! cam_exch
c
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
            if(theory.eq.'dft')rtdb_string = 'dft:gradient'
            if(theory.eq.'sodft')rtdb_string = 'sodft:gradient'
         else
            rtdb_string = 'scf:gradient'
         end if
         if (.not. rtdb_put(rtdb, rtdb_string, mt_dbl, 3*nat, 
     $        dbl_mb(k_force)))call errquit
     $        ('gradients: could not store gradients',1, RTDB_ERR)
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
         write(luout,1000) theory(1:inp_strlen(theory)),
     $        'x','y','z','x','y','z'
         do 30, i=1, nat
            if (.not. geom_cent_get(geom, i, tag, crd, q)) call errquit
     $           ('gradients: geometry corrupt?',0, GEOM_ERR)
            write(luout,2000) i, tag,(crd(j),j=1,3),
     $           (dbl_mb(k_force+3*(i-1)+j),j=0,2)
 30      continue
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
     $     call errquit('grad: ma corrupt',0, MA_ERR)

      if (.not. ma_chop_stack(l_force)) 
     $     call errquit('grad: failed chopping MA stack',0, MA_ERR)
c     
      if (.not. ga_destroy(g_force)) call errquit('grad: GA?',0, GA_ERR)
c     
      do i=1, ndens
         if (.not. ga_destroy(g_dens(i))) call errquit 
     $        ('error destroying density', 1, GA_ERR)
      end do
c     
      end
