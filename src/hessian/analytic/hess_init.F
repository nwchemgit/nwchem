C
C $Id: hess_init.F,v 1.12 2000-08-14 20:27:00 windus Exp $
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
#include "bas.fh"
#include "bgj.fh"
#include "dra.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "sym.fh"
c
c     integer rtdb, mult, g_vecs(2)
      integer rtdb, mult
      integer i, j, dftnoc(2), numfunc, nbsize
      parameter(numfunc = 20) ! took this from gradient code, may need to be
                              ! changed in the future
      logical status
      character*256 sx_name
      character*256 fx_name
      double precision zero, one, onem
      double precision xfac(numfunc)
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical int_normalize
      external int_normalize
c
c  get geometry and basis set handles and associated useful information
c
      if (.not. geom_create(geom, 'geometry'))
     *  call errquit('hess_init: geom_create problem', 555)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     *  call errquit('hess_init: no geometry ', 555)
      if (.not. geom_ncent(geom, ncent))
     *  call errquit('hess_init: unable to get ncent ', 555) 
      if (.not. bas_create(basis, 'ao basis'))
     *  call errquit('hess_init: bas_create problem', 555)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     *  call errquit('hess_init: no ao basis set', 555)
      if (.not. bas_numbf(basis, nbf)) 
     *  call errquit ('scf_init: basis info', 555)
c
c  initialize the second derivative integrals
c  making sure that the schwarz stuff is initialized
c
      write(6,*) 'before the call to intdd_init'
      call flush(6)
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
      call int_terminate()
      call intdd_init(rtdb,1,basis)
      write(6,*) 'after the call to intdd_init'
      call flush(6)
      status = int_normalize(rtdb,basis)
      write(6,*) 'after the call to int_normalize'
      call flush(6)
c
c Are we using symmetry?
c
      if (.not. rtdb_get(rtdb, 'scf:skeleton', mt_log, 1, skelsym)) then
         skelsym = sym_number_ops(geom) .gt. 0
      end if
c
c  get local memory for the nuclear hessian Exy
c
      n3xyz = 3 * ncent
      n3xyz2 = n3xyz * n3xyz
      if (.not.ma_push_get(MT_DBL,n3xyz2,'hessian: exy',l_exy,k_exy))
     *     call errquit('hessian_exy: cannot allocate',555)
      call dfill(n3xyz2,zero,dbl_mb(k_exy),1)
      if (.not.ma_push_get(MT_DBL,n3xyz2,'hessian: exy',l_hess_xc,
     *     k_hess_xc))
     *     call errquit('hessian_hess_xc: cannot allocate',555)
      call dfill(n3xyz2,zero,dbl_mb(k_hess_xc),1)
      write(6,*) 'after the calls to hessian ma_push_get'
      call flush(6)
c
c  Initialization for each of the different types of wavefunctions.
c
c    Need to fix jfac and kfac for UHF
c
      jfac(1) = 1.0d0   ! Only changed for DFT
      kfac(1) = 1.0d0
      if (theory .eq. 'hf') then
        if (.not. rtdb_get(rtdb, 'scf:nmo', MT_INT, 1, nmo))
     *    call errquit('hess_init: failed to get nmo', 555)
        if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     *     call errquit('hess_init: failed to get nclosed',0)
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
     *      call errquit('hess_init: failed to get nalpha',555)
          if (.not. rtdb_get(rtdb, 'scf:nbeta', MT_INT, 1, nbeta))
     *      call errquit('hess_init: failed to get nbeta',555)
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
        if (ocdfit) jfac(1) = 0.0d0 
        if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     $        call errquit('gradients: xfac not in rtdb ', 0)
        kfac(1) = xfac(1)
c
c  I don't like the way I am getting the number of mo's, but I can't find
c  another way to get this from the DFT.
c
        if (.not. bas_numbf(basis, nmo))
     *     call errquit('hess_init: error getting nmo:', 555)
        if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, dftnoc))
     *   call errquit('hess_init: rtdb_get of dftnoc failed', 555)
        if (bgj_restricted()) then
          nclosed = dftnoc(1)
          nalpha = nclosed
          nbeta = 0
          scftype = 'RHF'
        else
          nalpha = dftnoc(1)
          nbeta  = dftnoc(2) !Is this in the correct order?
          nclosed = nbeta
          scftype = 'UHF'
        endif

      endif
c
      if (scftype.eq.'UHF') then
         vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         nvirt = nmo - nclosed - nopen
         vlen = nclosed*(nopen+nvirt) + nopen*nvirt
      else
         call errquit('hess_init: unknown SCF type',0)
      endif
c
      if (.not. ma_push_get(MT_LOG, ncent, 'active atoms',
     &   l_act, k_act))
     &   call errquit('grad: could not allocate l_act',1)
      call grad_active_atoms(rtdb, ncent, log_mb(k_act), nactive)
      write(6,*) 'after the call to grad_active_atoms'
      call flush(6)
c
c  Get the density and global fock derivative matrices (Qs or RHS).
c
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
c     
      j = 0
      do i = 1, ncent
         j = j + 1
         g_rhs(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a x')
         call ga_zero(g_rhs(j))
         if (scftype .eq. 'UHF') then
            j = j + 1
            g_rhs(j) = ga_create_atom_blocked
     &           (geom, basis, 'CPKS RHS test b x')
            call ga_zero(g_rhs(j))
         endif
         j = j + 1
         g_rhs(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a y')
         call ga_zero(g_rhs(j))
         if (scftype .eq. 'UHF') then
            j = j + 1
            g_rhs(j) = ga_create_atom_blocked
     &           (geom, basis, 'CPKS RHS test b y')
            call ga_zero(g_rhs(j))
         endif
         j = j + 1
         g_rhs(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a z')
         call ga_zero(g_rhs(j))
         if (scftype .eq. 'UHF') then
            j = j + 1
            g_rhs(j) = ga_create_atom_blocked
     &           (geom, basis, 'CPKS RHS test b z')
            call ga_zero(g_rhs(j))
         endif
      enddo
      j = 0
      do i = 1, ncent
         j = j + 1
         g_rhs_xc(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a x')
         call ga_zero(g_rhs_xc(j))
         if (scftype .eq. 'UHF') then
            j = j + 1
            g_rhs_xc(j) = ga_create_atom_blocked
     &           (geom, basis, 'CPKS RHS test b x')
            call ga_zero(g_rhs_xc(j))
         endif
         j = j + 1
         g_rhs_xc(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a y')
         call ga_zero(g_rhs_xc(j))
         if (scftype .eq. 'UHF') then
            j = j + 1
            g_rhs_xc(j) = ga_create_atom_blocked
     &           (geom, basis, 'CPKS RHS test b y')
            call ga_zero(g_rhs_xc(j))
         endif
         j = j + 1
         g_rhs_xc(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a z')
         call ga_zero(g_rhs_xc(j))
         if (scftype .eq. 'UHF') then
            j = j + 1
            g_rhs_xc(j) = ga_create_atom_blocked
     &           (geom, basis, 'CPKS RHS test b z')
            call ga_zero(g_rhs_xc(j))
         endif
      enddo
c
c  Initialize the Disk Resident Array (DRA) interface which will be needed
c  to hold the overlap derivatives.
c
      call util_file_name('sxmatrices',.true.,.false.,sx_name)
      call util_file_name('fxmatrices',.true.,.false.,fx_name)
      nbsize = ma_sizeof(mt_dbl,ncent*3*nbf*nbf,mt_byte)
c     if (dra_init(1,nbsize,nbsize,-1).ne.0) call errquit
c    &    ('hess_init: Failed to initialize DRA interface',555)
      if (dra_init(2,nbsize,2*nbsize,-1).ne.0) call errquit
     &    ('hess_init: Failed to initialize DRA interface',555)
      if (dra_create(mt_dbl,nbf*ncent*3,nbf,'Sx matrices',sx_name,
     &    DRA_RW,nbf,nbf,sx_dra_handle).ne.0) call errquit
     &   ('hess_init: Failed to create sx DRA',555)
      if (dra_create(mt_dbl,nbf*ncent*3,nbf,'Fx matrices',fx_name,
     &    DRA_RW,nbf,nbf,fx_dra_handle).ne.0) call errquit
     &   ('hess_init: Failed to create fx DRA',555)

      return
      end
