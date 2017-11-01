C
      subroutine hess_init(rtdb)
c
c  This routine initializes many of the hessian variables and grabs memory
c  that will be in use during the whole calculation.
c  A lot of this could use a consolidation of the SCF and DFT energy
c  initialization routine.
c
      implicit none
c
#include "errquit.fh"
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
#include "case.fh"
#include "inp.fh"
c
      integer rtdb, mult
      integer dftnoc(2), nbsize
      logical oprint, omp2, odft
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
      logical use_nwxc,lhashfxc
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical int_normalize,int_norm_2c, movecs_read_header
      external int_normalize,int_norm_2c, movecs_read_header
      character *1024 ccc
      CHARACTER*3   crcr

      crcr = CHAR(13)//CHAR(11)//CHAR(0) 
c  print out follow information
c
      oprint = util_print('hess_follow',print_high)
      oprint = oprint .and. (ga_nodeid().eq.0)
c
c  get geometry and basis set handles and associated useful information
c
      odft = .false.
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
c
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
        odft = .true.
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
     $        call errquit('hess_init: xfac not in rtdb ', 0, RTDB_ERR)
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
          endif
        endif
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
c     if oactive business is on, bail out and provide recipe for restart
c
      if(nactive.ne.ncent) then
         ccc= crcr//' WARNING '//crcr
     c        //' Analytical Hessian code not ready to be'//crcr
     c        //' used with an active atom list or '//crcr
     c        //' with constraints for optimization.'//crcr
     c        //' If you want to continue this calculation and'//crcr
     c        //' get a hessian with negative frequencies,'//crcr
     c        //' please add the following line before the '//crcr
     c        //' line task dft frequencies'//crcr
     c        // ' '//crcr
     c        //' unset geometry:actlist '//crcr
     c        //' unset cons:* '//crcr
     c        //'  '//crcr
     c        //'  '//crcr
         
         call errquit(ccc(1:inp_strlen(ccc)),
     C        0,CAPMIS_ERR)
      endif
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
     $     nopen, nclosed, scftype, movecs_in, omp2, odft,
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
      if (.not. nga_create (MT_DBL, 3, dims, 'DFT CPKS RHS',chunk, 
     &     g_rhs_sx))
     &   call errquit('hess_init: could not allocate g_rhs_xc',555,
     &       GA_ERR)
      call ga_zero(g_rhs_sx)
      if (.not. nga_create (MT_DBL, 3, dims, 'DFT CPKS RHS',chunk, 
     &     g_rhs_fx))
     &   call errquit('hess_init: could not allocate g_rhs_xc',555,
     &       GA_ERR)
      call ga_zero(g_rhs_fx)
c
c  Initialize the Disk Resident Array (DRA) interface which will be needed
c  to hold the overlap derivatives.
c
c      if (oprint) write(LuOut,*) 'HESS: getting dra_handles'
c      call util_file_name('sxmatrices',.true.,.false.,sx_name)
c      call util_file_name('fxmatrices',.true.,.false.,fx_name)
c      nbsize = ma_sizeof(mt_dbl,ndens*ncent*3*nbf*nbf,mt_byte)
c     if (dra_init(1,nbsize,nbsize,-1).ne.0) call errquit
c    &    ('hess_init: Failed to initialize DRA interface',555)
c      if (dra_init(2,nbsize,2*nbsize,-1).ne.0) call errquit
c     &    ('hess_init: Failed to initialize DRA interface',555,
c     &       INT_ERR)
c      if (dra_create(mt_dbl,ndens*nbf*ncent*3,nbf,'Sx matrices',
c     &    sx_name,
c     &    DRA_RW,nbf,nbf,sx_dra_handle).ne.0) call errquit
c     &   ('hess_init: Failed to create sx DRA',555, INT_ERR)
c      if (dra_create(mt_dbl,ndens*nbf*ncent*3,nbf,'Fx matrices',
c     &    fx_name,
c     &    DRA_RW,nbf,nbf,fx_dra_handle).ne.0) call errquit
c     &   ('hess_init: Failed to create fx DRA',555, INPUT_ERR)


      return
      end
c $Id$
