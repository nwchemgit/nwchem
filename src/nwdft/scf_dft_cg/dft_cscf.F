      subroutine dft_cscf_init(rtdb,dft_conv_tol,dft_charge,dft_oskel,
     &           dft_oadapt,dft_maxiter,dft_geom,dft_basis,dft_nopen,
     &           dft_nbf,dft_noc,dft_movecs_in,dft_movecs_out,
     &           dft_movecs_guess,dft_title)
      implicit none
c
c     The quasi-Newton-Raphson DFT solver is a combination of the
c     Hartree-Fock solver with a bit of DFT (fock_xc) added in. As
c     a result the functionality is controlled by parameters that
c     traditionally are part of the SCF energy expression as well as
c     parameters that belong to DFT. However, due to their independent
c     history there are also parameters that appear in both such as
c     the number of basis functions (nbf).
c
c     The consequence is that the common blocks in cdft.fh and cscf.fh
c     have name clashes, yet their settings need to be consistent for
c     the code to work correctly. Hence this subroutine that takes
c     variables from the DFT code and assign their values to variables
c     in the SCF common blocks. This way consistent settings can be
c     assured without variable name clashes.
c
c     In addition to setting these variables this subroutine may set
c     variables that are SCF specific and have no meaning to the DFT
c     code and are not relevant to execution path decisions. These
c     values may be obtained from the runtime database.
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "cscf.fh"
c
      integer rtdb                   ! handle of the runtime database
      double precision dft_conv_tol  ! DFT version of gnorm_tol
      double precision dft_charge    ! DFT version of charge
      logical dft_oskel              ! DFT version of oskel
      logical dft_oadapt             ! DFT version of oadapt
      integer dft_maxiter            ! DFT version of maxiter
      integer dft_geom               ! DFT version of geom
      integer dft_basis              ! DFT version of basis
      integer dft_nopen              ! DFT version of nopen
      integer dft_nbf                ! DFT version of nbf
      integer dft_noc(2)             ! DFT version of nalpha and nbeta
c     integer dft_movecs(2)          ! DFT version of MO vectors handles
      character*(*) dft_movecs_in    ! DFT version of movecs_in
      character*(*) dft_movecs_out   ! DFT version of movecs_out
      character*(*) dft_movecs_guess ! DFT version of movecs_guess
      character*(*) dft_title        ! DFT version of title
c
c     local stuff
c
c     integer dft_noc(2) ! number of occupied orbitals
c
c     Pass data from arguments to cscf.fh
c
      gnorm_tol = dft_conv_tol
      charge = dft_charge
      oskel = dft_oskel
      oadapt = dft_oadapt
      maxiter = dft_maxiter
      geom = dft_geom
      basis = dft_basis
      nopen = dft_nopen
      nbf = dft_nbf
      nmo = dft_nbf
      nalpha = dft_noc(1)
      nbeta = dft_noc(2)
      movecs_in = dft_movecs_in
      movecs_out = dft_movecs_out
      movecs_guess = dft_movecs_guess
      title = dft_title
c
c     Initialise settings that make little sense in quasi-Newton-Raphson
c     DFT code (mostly because we would use the DFT_SCF code for those).
c
      ododiag = .false.
      odiis = .false.
c
c     Work out remaining stuff
c
      if (nopen.eq.0) then
        nelec = 2*nalpha
        nclosed = nalpha
      else
        nelec = nalpha + nbeta
        nclosed = min(nalpha,nbeta)
      endif
c
c     Get stuff from RTDB (for SCF driver specific things)
c
c     if (rtdb_get(rtdb, 'dft:level shift info', MT_DBL, 6,shifts)) then
c       ouser_changed_conv = .true.
c     else
c       call dfill(6, -1.0d0, shifts, 1)
c     endif
c
c     if (rtdb_get(rtdb, 'dft:full hessian switch', MT_DBL, 1,
c    $     nr_gswitch)) then
c       ouser_changed_conv = .true.
c     else
c       nr_gswitch = 0.1d0
c     endif
c
c     if (.not. rtdb_get(rtdb, 'dft:tol2e', MT_DBL, 1, tol2e))
c    $     tol2e = min(1.0d-7,gnorm_tol*1d-2)
c
c     if (rtdb_get(rtdb, 'dft:ri', MT_INT, 1, nriscf)) then
c       if (.not. bas_create(riscf_basis, 'riscf basis'))
c    $       call errquit('dft_cscf_int: bas_create?', 0, BASIS_ERR)
c       if (.not. bas_rtdb_load(rtdb, geom, riscf_basis, 'riscf basis'))
c    $       call errquit('dft_cscf_int: no riscf basis set', 0,
c    $                    RTDB_ERR)
c     else
c       nriscf = 0
c       riscf_basis = 0
c     endif
c
c     if (.not.rtdb_get(rtdb, 'scf:lock',MT_LOG, 1, olock)) then
c        olock = .false.
c     endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dft_cscf_scftype(dft_scftype)
      implicit none
#include "cscf.fh"
      character*(*) dft_scftype
      scftype = dft_scftype
      end
c
c-----------------------------------------------------------------------
c
      logical function dft_cscf_equals_scftype(dft_scftype)
      implicit none
#include "cscf.fh"
      character*(*) dft_scftype
      dft_cscf_equals_scftype = dft_scftype.eq.scftype
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dft_cscf_movecs(dft_movecs)
      implicit none
c
c     Store the GA handles for the MO-vectors as defined in the DFT
c     code into the SCF common block. This way they can be used in
c     the quasi-Newton-Raphson solver as well.
c
#include "cscf.fh"
      integer dft_movecs(2)
      g_movecs(1) = dft_movecs(1)
      g_movecs(2) = dft_movecs(2)
      end
c
c-----------------------------------------------------------------------
c
      subroutine dft_cscf_movecs_inout(dft_movecs_in,dft_movecs_out)
      implicit none
c
c     Store the locations of the guess vector and final vectors
c
#include "cscf.fh"
      character*(*) dft_movecs_in, dft_movecs_out
      movecs_in  = dft_movecs_in
      movecs_out = dft_movecs_out
      end
c $Id$
