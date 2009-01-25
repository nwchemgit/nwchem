C
C $Id: hess_init.F,v 1.63 2009-01-25 02:27:01 niri Exp $
C
      subroutine hess_init(rtdb)
c
c  This routine initializes many of the hessian variables and grabs memory
c  that will be in use during the whole calculation.
c  A lot of this could use a consolidation of the SCF and DFT energy
c  initialization routine.
c
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "bgj.fh"
#include "dra.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "sym.fh"
#include "util.fh"
#include "xc.fh"
c
c     integer rtdb, mult, g_vecs(2)
      integer rtdb, mult
      integer dftnoc(2), nbsize
      logical oprint, omp2
      logical cdfit, frac_occ
      integer nbf_cd, nbf2, mxnocc
      integer g_eigen_diag
      integer l_evals, k_evals, l_occ, k_occ, k_pdm2, k_coeff
      integer dims(3), chunk(3)
      character*256 sx_name
      character*256 fx_name
      character*255 movecs_in, title_vecs, basis_vecs
      character*20 scftype_vecs
      integer nbf_vecs, nsets, nmos(2)
      double precision zero, one, onem
      double precision xfac(numfunc)
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical int_normalize,int_norm_2c, movecs_read_header
      external int_normalize,int_norm_2c, movecs_read_header
c
c  print out follow information
c
      oprint = util_print('hess_follow',print_high)
      oprint = oprint .and. (ga_nodeid().eq.0)
c
c  get geometry and basis set handles and associated useful information
c
      if (oprint) write(LuOut,*) 'HESS: set up geometry and basis'
      if (.not. geom_create(geom, 'geometry'))
     *  call errquit('hess_init: geom_create problem', 555, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     *  call errquit('hess_init: no geometry ', 555, RTDB_ERR)
      if (.not. geom_ncent(geom, ncent))
     *  call errquit('hess_init: unable to get ncent ', 555, GEOM_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     *  call errquit('hess_init: bas_create problem', 555, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     *  call errquit('hess_init: no ao basis set', 555, BASIS_ERR)
      if (.not. bas_numbf(basis, nbf)) 
     *  call errquit ('hess_init: basis info', 555, BASIS_ERR)
      if (.not. bas_create(cd_basis, 'cd basis'))
     &   call errquit('hess_init: bas_create?',102, BASIS_ERR)
      cdfit=bas_rtdb_load(rtdb, geom, cd_basis, 'cd basis')
      if (cdfit)then
         if (.not. bas_numbf(cd_basis, nbf_cd)) 
     &      call errquit('hess_init: basis set error:', 86, BASIS_ERR)
      else
        if (.not. bas_destroy(cd_basis))
     &     call errquit('rdinput:not able to destroy CD_bas:',86,
     &       BASIS_ERR)
        cd_basis=0
      endif
c
c  Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hess_init: theory not specified',555, RTDB_ERR)
      if (theory.eq.'hyb') theory = 'dft'
      if (theory.eq.'scf') theory = 'hf'
c
c  initialize the second derivative integrals
c  making sure that the schwarz stuff is initialized
c
      if (oprint) write(LuOut,*) 'HESS: initializing integrals'
      nbasis=1
      basis_list(1)=basis
      if(cdfit) then
         nbasis=2
         basis_list(2)=cd_basis
      endif
      call int_init(rtdb, nbasis, basis_list)
      call schwarz_init(geom, basis)
      call int_terminate()
      call intdd_init(rtdb,nbasis,basis_list)
      if(.not.int_normalize(rtdb,basis))
     &        call errquit('hess_init: int_normalize failed', 208,
     &       INT_ERR)
      if(cdfit) then
         if (.not. int_norm_2c(rtdb,cd_basis))
     &        call errquit('hess_init: int_norm_2c failed', 208,
     &       INT_ERR)
      endif
c
c Are we using symmetry?  Assume not!
c
c     if (.not. rtdb_get(rtdb, 'scf:skeleton', mt_log, 1, skelsym)) then
c        skelsym = sym_number_ops(geom) .gt. 0
c     end if
      skelsym = .false.
c
c  get local memory for the nuclear hessian Exy and dipole derivative
c  tensor
c
      if (oprint) write(LuOut,*) 'HESS: initializing hessian'
      n3xyz = 3 * ncent
      n3xyz2 = n3xyz * n3xyz
      if (.not.ma_push_get(MT_DBL,n3xyz2,'hessian: exy',l_exy,k_exy))
     *     call errquit('hessian_exy: cannot allocate',555, MA_ERR)
      call dfill(n3xyz2,zero,dbl_mb(k_exy),1)
      if (.not.ma_push_get(MT_DBL,n3xyz2,'hessian: exy',l_hess_xc,
     *     k_hess_xc))
     *     call errquit('hessian_hess_xc: cannot allocate',555, MA_ERR)
      call dfill(n3xyz2,zero,dbl_mb(k_hess_xc),1)
      if (.not.ma_push_get(MT_DBL,3*n3xyz,'hessian: dder',l_dder,
     *     k_dder))
     *     call errquit('hessian_dder: cannot allocate',555, MA_ERR)
      call dfill(n3xyz*3,zero,dbl_mb(k_dder),1)
c
c  Initialization for each of the different types of wavefunctions.
c
      if (oprint) write(LuOut,*) 'HESS: getting wavefunction info'
      jfac(1) = 1.0d0
      kfac(1) = -0.5d0
      if (theory .eq. 'hf') then
c
c   Read number of mo's off of the movecs file
c
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $    call errquit('hess_init: failed to get input vectors', 555,
     &       RTDB_ERR)
        if (.not.movecs_read_header(movecs_in,title_vecs,basis_vecs,
     $    scftype_vecs, nbf_vecs,nsets,nmos,2))
     $    call errquit('hess_init:problem with movecs_read_header',1,
     &       DISK_ERR)
        if (nbf.ne.nbf_vecs) 
     $    call errquit('hessint: nbf mismatch?',0, INPUT_ERR)
        nmo = nmos(1)
c        if (.not. rtdb_get(rtdb, 'scf:nmo', MT_INT, 1, nmo))
c     *    call errquit('hess_init: failed to get nmo', 555)
        if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     *     call errquit('hess_init: failed to get nclosed',0, RTDB_ERR)
        nalpha = nclosed
        nbeta = 0
        if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     *    nopen = 0
c       nopen = nopen + 1
        if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
          if (nopen .eq. 0) then
            scftype = 'RHF'
          else if (scftype .ne. 'UHF') then
            scftype = 'ROHF'
          endif
        endif
        if (scftype .eq. 'UHF') then
          if (.not. rtdb_get(rtdb, 'scf:nalpha', MT_INT, 1, nalpha))
     *      call errquit('hess_init: failed to get nalpha',555,
     &       RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'scf:nbeta', MT_INT, 1, nbeta))
     *      call errquit('hess_init: failed to get nbeta',555,
     &       RTDB_ERR)
          jfac(2) = 1.0d0
          kfac(1) = 1.0d0
          kfac(2) = 1.0d0
        endif
      elseif (theory .eq. 'dft') then
        if (.not. rtdb_get(rtdb, 'dft:mult', MT_INT, 1,mult))
     *    mult = 1
        nopen = mult - 1
c
c     Get information about fitting and scale factors
c
        if (.not. rtdb_get(rtdb, 'dft:cdfit', mt_log, 1, ocdfit))
     $        ocdfit = .false.
        if(oprint) then
           write(luout,*) ' ocdfit ',ocdfit
        endif
        if (ocdfit) jfac(1) = 0.0d0 
        if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     $        call errquit('gradients: xfac not in rtdb ', 0, RTDB_ERR)
        kfac(1) = xfac(1)*(-0.5d0)
c
c  I don't like the way I am getting the number of mo's, but I can't find
c  another way to get this from the DFT.
c
c
c   Read number of mo's off of the movecs file
c
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $    call errquit('hess_init: failed to get input vectors', 555,
     &       RTDB_ERR)
        if (.not.movecs_read_header(movecs_in,title_vecs,basis_vecs,
     $    scftype_vecs, nbf_vecs,nsets,nmos,2))
     $    call errquit('hess_init:problem with movecs_read_header',1,
     &       DISK_ERR)
        if (nbf.ne.nbf_vecs) 
     $    call errquit('hessint: nbf mismatch?',0, INPUT_ERR)
        nmo = nmos(1)
c        if (.not. bas_numbf(basis, nmo))
c     *     call errquit('hess_init: error getting nmo:', 555)
        if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, dftnoc))
     *   call errquit('hess_init: rtdb_get of dftnoc failed', 555,
     &       RTDB_ERR)
        if (bgj_restricted()) then
          nclosed = dftnoc(1)
          nalpha = nclosed
          nbeta = 0
          scftype = 'RHF'
        else
          nalpha = dftnoc(1)
          nbeta  = dftnoc(2)
          nclosed = nbeta    ! assumes that there are more alphas than betas
          scftype = 'UHF'
          if (ocdfit) then
            jfac(2) = 0.0d0
          else
            jfac(2) = 1.0d0
          endif
          kfac(1) = 1.0d0*xfac(1)
          kfac(2) = 1.0d0*xfac(1)
          if (.not. rtdb_put(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     *         call errquit('scfinit:rtdbput nopen failed',nopen,
     &       RTDB_ERR)
        endif

      endif
c
      if (scftype.eq.'UHF') then
         if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, scftype))
     &        call errquit("hess_init: cannot put scftype", 0,
     &       RTDB_ERR)
         nocc(1) = nalpha
         nocc(2) = nbeta
         nvirt(1) = nmo-nalpha
         nvirt(2) = nmo-nbeta
         vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         nocc(1) = nclosed
         nvirt(1) = nmo - nclosed - nopen
         vlen = nclosed*(nopen+nvirt(1)) + nopen*nvirt(1)
      else
         call errquit('hess_init: unknown SCF type',0, INPUT_ERR)
      endif
c
      if (.not. ma_push_get(MT_LOG, ncent, 'active atoms',
     &   l_act, k_act))
     &   call errquit('grad: could not allocate l_act',1, MA_ERR)
      call grad_active_atoms(rtdb, ncent, log_mb(k_act), nactive)
c
c  Get the density and weighted density.  Need to set up some
c  temporary memory to go into grad_dens.
c
      if (oprint) write(LuOut,*) 'HESS: getting ga_memory'
      if (scftype .eq. 'UHF') then
        g_dens(1) = ga_create_atom_blocked
     *            (geom, basis, 'hessian:a dens')
        g_dens(2) = ga_create_atom_blocked
     *            (geom, basis, 'hessian:b dens')
        ndens = 2
      else
        g_dens(1) = ga_create_atom_blocked
     *            (geom, basis, 'hessian:dens')
        ndens = 1
      endif
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')
      nbf2 = nbf
      if (scftype .eq. 'UHF') nbf2 = nbf + nbf
      mxnocc = nopen + nclosed
c
c     eigenvalue diagonal matrix for forming energy weighted density
c
      if (.not. ga_create(MT_DBL, mxnocc, mxnocc,'eigen_diag', 0,
     $  mxnocc, g_eigen_diag)) call errquit('onedd_cont: ga diags?',
     $  mxnocc*mxnocc, GA_ERR)
      call ga_zero(g_eigen_diag)
C     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('onedd_cont: could not allocate l_evals',nbf2,
     &       MA_ERR)
C     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf2,'occ. numbers', l_occ, k_occ))
     $     call errquit('onedd_cont: could not allocate l_occ',nbf2,
     &       MA_ERR)
c
c     check if smearing is on
c
      if (.not.rtdb_get(rtdb, 'dft:fractional_occup', mt_log, 1,
     & frac_occ  )) frac_occ=.false.
      omp2 = .false.
      k_pdm2 = 1  ! dummy variables to pass to grad_dens
      k_coeff = 1
c
c get MO vectors
c
      if (theory.eq.'dft') then
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     *    call errquit('onedd_cont: DFT MO vectors not defined',555,
     &       RTDB_ERR)
      else
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $    call errquit('onedd_cont: SCF MO vectors not defined',0,
     &       RTDB_ERR)
      endif
c
c get density and weighted density
c
      if (oprint) write(LuOut,*) 'HESS: getting densities'
      call grad_dens(geom, basis, g_dens, g_wdens, g_eigen_diag,
     $     dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2,
     $     nopen, nclosed, scftype, movecs_in, omp2,
     $     dbl_mb(k_pdm2), dbl_mb(k_coeff),frac_occ,.true.)
c
c     free temporary arrays
c
      if (.not. ga_destroy(g_eigen_diag))
     $  call errquit('onedd_cont: could not destroy g_eigen_diag',1,
     &       GA_ERR)
 
      if (.not.ma_pop_stack(l_occ))
     $     call errquit('onedd_cont:ma free occ',1, MA_ERR)
      if (.not.ma_pop_stack(l_evals))
     $     call errquit('onedd_cont:ma free eval',1, MA_ERR)
c
c Get the Fock derivative matrices (rhs), note that we get a seperate
c one for the DFT code so that we can debug for the moment
c
      dims(1) = ndens*n3xyz
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'CPKS RHS',chunk, 
     &     g_rhs))
     &   call errquit('hess_init: could not allocate g_rhs',555, GA_ERR)
      call ga_zero(g_rhs)
      if (.not. nga_create (MT_DBL, 3, dims, 'DFT CPKS RHS',chunk, 
     &     g_rhs_xc))
     &   call errquit('hess_init: could not allocate g_rhs_xc',555,
     &       GA_ERR)
      call ga_zero(g_rhs_xc)
c
c  Initialize the Disk Resident Array (DRA) interface which will be needed
c  to hold the overlap derivatives.
c
      if (oprint) write(LuOut,*) 'HESS: getting dra_handles'
      call util_file_name('sxmatrices',.true.,.false.,sx_name)
      call util_file_name('fxmatrices',.true.,.false.,fx_name)
      nbsize = ma_sizeof(mt_dbl,ndens*ncent*3*nbf*nbf,mt_byte)
c     if (dra_init(1,nbsize,nbsize,-1).ne.0) call errquit
c    &    ('hess_init: Failed to initialize DRA interface',555)
      if (dra_init(2,nbsize,2*nbsize,-1).ne.0) call errquit
     &    ('hess_init: Failed to initialize DRA interface',555,
     &       INT_ERR)
      if (dra_create(mt_dbl,ndens*nbf*ncent*3,nbf,'Sx matrices',
     &    sx_name,
     &    DRA_RW,nbf,nbf,sx_dra_handle).ne.0) call errquit
     &   ('hess_init: Failed to create sx DRA',555, INT_ERR)
      if (dra_create(mt_dbl,ndens*nbf*ncent*3,nbf,'Fx matrices',
     &    fx_name,
     &    DRA_RW,nbf,nbf,fx_dra_handle).ne.0) call errquit
     &   ('hess_init: Failed to create fx DRA',555, INPUT_ERR)


      return
      end
c
c
      subroutine hess_energyinit(rtdb)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
c  This routine does a few checks to make sure that the energy will
c  be calculated without symmetry and with reasonable thresholds so
c  the rest of the Hessian calculation will work properly.  Note that
c  these values should be reset after the calculation is done.
c
      integer rtdb
      integer geomnew, filesize, memsize
      double precision ithresh,tol_rho,tol_rho_rd,e_conv,
     .     job_grid_acc
      logical xc_got2nd      
      external xc_got2nd
      logical grid_clinit,gstatus
      external grid_clinit
      character*15 gridtype
c     
c  Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hess_energyinit: theory not specified',555,
     &       INPUT_ERR)
      if (theory.eq.'hyb') theory = 'dft'
c
c  Check to see if the current geometry has symmetry.  If so, create
c  a new geometry without symmetry and put onto the rtdb.  When symmetry
c  is implemented into the code, this can go away.
c
      if (.not. geom_create(geomold, 'geometry'))
     *  call errquit('hess_init: geom_create problem', 555, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geomold, 'geometry'))
     *  call errquit('hess_init: no geometry ', 555, GEOM_ERR)
      if (.not. geom_makec1(geomold, geomnew))
     $  call errquit('hess_init: failure to decrease symmetry',555,
     &       GEOM_ERR)
      if (.not.geom_rtdb_store(rtdb,geomnew,'geometry'))
     &      call errquit(' hess_init: geom_rtdb_store failed',555,
     &       RTDB_ERR)
      if (.not.geom_destroy(geomnew))
     &      call errquit(' hess_init: geom_destroy failed',555,
     &       GEOM_ERR)
c
c  Need to check to make sure that the convergence tolerances for
c  the wavefunction are OK.  Use the user defined one if available,
c  otherwise check to make sure that the threshold gets set to a
c  "reasonable" value.
c
      if (rtdb_get(rtdb, 'hessian:thresh', MT_DBL,1, ithresh)) then
          if (.not. rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (.not. rtdb_put(rtdb, 'scf:thresh', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        if (theory.eq.'scf') then
        else
          if (.not. rtdb_get(rtdb, 'dft:g_conv', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (.not. rtdb_put(rtdb, 'dft:g_conv', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        endif
      else
        ithresh = 1.0d-6
          if (.not. rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (tthresh.lt.ithresh) ithresh = tthresh
          if (.not. rtdb_put(rtdb, 'scf:thresh', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        if (theory.eq.'scf') then
        else
        ithresh = 1.0d-6
          if (.not. rtdb_get(rtdb, 'dft:g_conv', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (tthresh.lt.ithresh) ithresh = tthresh
          if (.not. rtdb_put(rtdb, 'dft:g_conv', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        endif
      endif
      if(theory.eq.'dft') then
         if(.not.xc_got2nd()) call errquit(
     '        'analytic 2nds not ready for these XC functionals',0,
     &       CAPMIS_ERR)
c
c     set grid to fine or leave the user to set it to xfine
c
         if (rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype)) then
            if(gridtype.eq.'coarse'.or.gridtype.eq.'xcoarse'.
     $           or.gridtype.eq.'medium') gridtype='fine'
         else
            gridtype='fine'
         endif
 
         if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     $        call errquit('hessinit:rtdbcput failed',0, RTDB_ERR)
c
c     set grid_acc consistent with fine (at least), ie 1d-7
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc)) then
          job_grid_acc = 1d-7
      else
         job_grid_acc = min(job_grid_acc,1d-7)
      endif
      if (.not. rtdb_put(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))call errquit('hessint:rtdbput failed',1,
     &       RTDB_ERR)
c
c     weight derivatives bit turned off
c
      if (.not. rtdb_put(rtdb, 'dft:weight derivatives', mt_log, 1, 
     & .false.))call errquit('hessint:rtdbput failed',2, RTDB_ERR)
      
c
c     set tol_rho to at least 1d-11 (it is 1d-10 otherwise)
c
         tol_rho=1d-11
         if (rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho_rd))
     .        tol_rho=min(tol_rho,tol_rho_rd)
         if (.not.rtdb_put(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho)) 
     .        call errquit(' hessinit: rtdbput failed ',11,0)
c     
c     kill symmetry for XC grid build
c     
         if (.not. rtdb_put(rtdb, 'dft:gridsym',
     $        MT_LOG, 1, .false.)) call errquit(
     $        ' task_hessian: rtdb_put failed',0, RTDB_ERR)
c
c     get xfac for doing xc_gotxc right on restart
c     but ignore the status since there are not xfac 
c     on rtdb when we are in start mode
c
         gstatus=grid_clinit(rtdb)
         call grid_cleanup(.false.)
c
c     set deltaE at least to 1d-7 
c
      if (.not.rtdb_get(rtdb,'dft:e_conv',mt_dbl,1,e_conv))then
        e_conv=1d-7
      else
        e_conv=min(e_conv,1d-7)
      endif
      if (.not. rtdb_put(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))
     &   call errquit('hess_init: rtdb_put failed', 100,
     &       RTDB_ERR)

       call dft_inpanae(rtdb)
      endif
c
c     make the Fock build go direct
c
      filesize = -1
      memsize  = -1
      if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, filesize))
     $  call errquit('hess_energyinit: failed to store filesize',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $  call errquit('hess_energyinit: failed to store filesize',0,
     &       RTDB_ERR)
c
      return
      end
c
      logical function hess_check(rtdb)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer rtdb      ! [input] rtdb handle
c
      integer geom1, basis1
      integer ipol, mult
      logical cdfit, dmd_rel, dk_rel, zora_rel
      logical xc_chktau,xc_getxcfac,cam_exch
      external xc_chktau,xc_getxcfac
c
      hess_check = .false.
      cdfit = .false.
      ipol = 1
c
c  Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hess_check: theory not specified',555,
     &       RTDB_ERR)
      if (theory.eq.'hyb') theory = 'dft'
c
c  Check for relativistic all-electron methods
c
      if (.not. rtdb_get(rtdb, 'dyall_mod_dir', MT_LOG, 1, dmd_rel))
     $   dmd_rel = .false.
      if (.not. rtdb_get(rtdb, 'doug_kroll', MT_LOG, 1, dk_rel))
     $   dk_rel = .false.
      if (.not. rtdb_get(rtdb, 'zora', MT_LOG, 1, zora_rel))
     $   zora_rel = .false.
c 
      if (dmd_rel.or.dk_rel.or.zora_rel) then
         if (ga_nodeid().eq.0)
     $    write(LuOut,*) 'Relativistic all-electron method detected'
         return
      endif
c
c  Check for ROHF in SCF
c
      if (theory.eq.'scf') then
        if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     *    nopen = 0
        if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
          if (nopen .eq. 0) then
            scftype = 'RHF'
          else if (scftype .ne. 'UHF') then
            scftype = 'ROHF'
          endif
        endif
        if (scftype.eq.'ROHF') then
           if (ga_nodeid().eq.0)
     *       write(LuOut,*) 'ROHF detected'
           return
        endif
      endif
c
c  Check for UHF in DFT
c
      if (theory.eq.'dft') then
        if (.not. rtdb_get(rtdb, 'dft:nopen', MT_INT, 1, nopen))
     $    nopen = 0
        if (.not. rtdb_get(rtdb, 'dft:mult', MT_INT, 1, mult))
     $    mult = 1
        if (.not. rtdb_get(rtdb, 'dft:ipol', MT_INT, 1, ipol))
     $    ipol = 1
        if (ipol.gt.1.or.nopen.gt.0.or.mult.gt.1) then
          if (ga_nodeid().eq.0)
     $      write(LuOut,*) 'Open shell DFT detected'
          return
        endif
      endif
c
c  Check for CD fit
c
      if (theory.eq.'dft') then
         if (.not. geom_create(geom1, 'geometry'))
     $        call errquit('hess_check: geom_create problem', 555,
     &       GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom1, 'geometry'))
     $        call errquit('hess_check: no geometry ', 555,
     &       RTDB_ERR)
         if (.not. bas_create(basis1, 'cd basis'))
     $        call errquit('hess_check: bas_create?',102, BASIS_ERR)
         cdfit=bas_rtdb_load(rtdb, geom1, basis1, 'cd basis')
         if (.not. geom_destroy(geom1))
     $        call errquit('hess_check:not able to destroy geometry:',
     $        geom1, GEOM_ERR)
         if (.not. bas_destroy(basis1))
     $        call errquit('hess_check:not able to destroy CD_bas:',
     $                   basis1, BASIS_ERR)
      endif
c
c     check for metaGGAs
c
      if(theory.eq.'dft') then
         if(xc_getxcfac(rtdb)) then
            if(xc_chktau()) then
               if (ga_nodeid().eq.0)
     $              write(LuOut,*) 'MetaGGA xc detected'
               return
            endif
         endif
      endif
c
c     check for attenuated exchange-correlation forms
c
      cam_exch = .false.
      if(theory.eq.'dft') then
        if (.not.rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1,
     &   cam_exch)) cam_exch=.false.
         if (cam_exch) then
             if (ga_nodeid().eq.0) 
     &              write(LuOut,*) 'CAM/LC xc detected'
             return
         end if
      end if
c
c     check for cdfit
c      
      if (cdfit) then
        if (ga_nodeid().eq.0)
     $    write(LuOut,*) 'CD fit detected'
        return
      endif
c
      hess_check = .true.
c
      return
      end
