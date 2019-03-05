      subroutine aoresponse_driver (rtdb, basis, geom)

c     $Id$

c     ==================================================================

c     main driver routine for "aoresponse", nwchem implementation
c
c     purpose: perform frequency dependent SCF response calculation
c
c     This subroutine was
c     initially adapted from the property module hnd_giaox.F and
c     the subroutine "apt1w" from the original aoresponse module
c     which was interfaced with the ADF code (www.scm.com)

c     input:   rtdb  - nwchem run time database 
c              basis - basis set handle
c              geom  - geometry handle

c     The "aoresponse" software is a code package being developed
c     by the research group of Jochen Autschbach, State University
c     of New York at Buffalo (SUNY Buffalo). For questions please
c     contact us by email at jochena@buffalo.edu. 
c
c     some references to the methods being used here are:
c     JPC-A 110 (2006), 2461
c     JCP 123 (2005), 114103
c     JCP 122 (2005), 224115
c     JCP 122 (2005), 074105
c     and a number of submitted papers as well as the
c     work by D. P. Santry and coworkers from the 1970s. Check out
c    http://www.nsm.buffalo.edu/~jochena/publications/publications.html
c     for follow- papers by our group.

c     note: most of the ADF aoresponse functionality is not yet
c     available in NWChem (July 2006, update Jan 2011)

c     ==================================================================

      implicit none

#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "util.fh"
c#include "dimqm.fh"


      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
c
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer l_dia, k_dia, l_para, k_para
      integer l_xyz, k_xyz, l_zan, k_zan
      integer l_tmp, k_tmp
      integer g_dens(3), g_s10, g_h01, g_h11, g_d1, g_rhs, g_fock, g_u
      integer g_vectors(2), geomnew, i, j, ij, g_xc(3)
      double precision atn, tol2e, val, isotr, aniso, tolfr
      double precision jfac(3),kfac(3),a(6),axs(3,3),eig(3),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
c
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
c
      logical     oskel, status
      double precision ppm
      data ppm     /26.62566914d+00/ 
      data tol2e   /1.0d-10/
      data tolfr   /1.0d-50/

      integer nfreq, l_fr, k_fr, l_deps, k_deps, nocv
      data nfreq/0/  ! number of frequencies, static if zero
      integer response_order
      data response_order/1/ ! 1 for linear, 2 for quadratic, etc

      double precision dummy(3)
      double precision origin(3)
      data origin/0d0,0d0,0d0/

      integer nat
      parameter (nat=1)

      double precision betare(3,3), betaim(3,3)
      double precision alfare(3,3), alfaim(3,3)

c     GA handles for one-electron integrals and temp arrays
      integer g_dipel, g_smat0, g_temp, g_dipmag, g_smat1,
     &   g_sket1, g_quadel

c     GA handle for magnetic field perturbed MO vectors:
      integer g_vecB1

      integer idir, ipert ! field directions
      integer nbfx ! later set to nbf*(nbf+1)/2
      integer ntmp, npol

      double precision sum, sumim, zero, one, two
      parameter (zero = 0d0, one = 1d0, two = 2d0)
      double precision gamwidth
      double precision ga_trace_diag
      external ga_trace_diag

      logical debug, lgiao, lvelocity, lmagpert, lstatic, lifetime,
     $        lbeta, lgprime, lquad
      logical oprint

      double precision alfare2(3,3), alfaim2(3,3)
      double precision auXnm   ! Conversion factor from a.u. to nm
      parameter (auXnm=45.563353d0)
      double precision alphaiso     ! Isotropic dipole polarizability
      double precision alphaani     ! Anisotropic dipole polarizability
      double precision alphaevr(3)  ! Dipole polarizability tensor eigenvalues (real)
      double precision alphaevi(3)  ! Dipole polarizability tensor eigenvalues (imag)
      double precision aeigvec(3,3) ! dummy dgeev eigenvectors
      double precision awork(9)     ! dgeev workspace
      double precision alpha_rtdb(12) ! omega, XX, XY, XZ, YY, YZ, ZZ, LL, MM, NN, Iso, Ani
      logical lclfld,use_dimqm,ldimqm

      integer aeiginfo              ! dgeev information
cccc#endif

c     ==================================================================

      debug = .false. .and. ga_nodeid().eq.0 ! .true. during development

      oprint = util_print('information', print_low)
      oprint = oprint .and. (ga_nodeid().eq.0)
      use_dimqm=util_module_avail("dimqm")

      npol = 1 ! we assume closed shell
      
c     ----------------------
c     start-up: print header
c     ----------------------
      
      if (oprint) then
c jbecca
c        write (LuOut,*) 'Entering AOResponse driver routine'
        write (LuOut,*) 'Entering old AOResponse driver routine'
c jbecca end
      endif
      
      if (oprint) then
        write (luout,9000)
        write (luout,9010) 
      end if
      
c     -------------------------------------------
c     set defaults for options and read from RTDB
c     -------------------------------------------

      lgiao     = .false. 
      lmagpert  = .false.
      lvelocity = .false.
      lifetime  = .false.
      lbeta     = .false.
      lgprime   = .false.
      lquad     = .false.
      gamwidth  = 0d0

      ntmp = 1
      if (.not.rtdb_get(rtdb, 'aoresponse:velocity', mt_int, 1,ntmp)) 
     &    call errquit('aor_driver: velocity failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lvelocity = .true.

      ntmp = 1
      if (.not.rtdb_get(rtdb,'aoresponse:magneticpert', mt_int, 1,ntmp))
     &    call errquit('aor_driver: magneticpert failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lmagpert = .true.

      ntmp = 1
      if (.not.rtdb_get(rtdb, 'aoresponse:giao', mt_int, 1,ntmp)) 
     &    call errquit('aor_driver: giao failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lgiao = .true.

      ntmp = 1
      if (.not.rtdb_get(rtdb, 'aoresponse:orbeta', mt_int, 1,ntmp)) 
     &    call errquit('aor_driver: lbeta failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lbeta = .true.

      ntmp = 1
      if (.not.rtdb_get(rtdb, 'aoresponse:gprime', mt_int, 1,ntmp)) 
     &   call errquit('aor_driver: lgprime failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lgprime = .true.         

      ntmp = 1
      if (.not.rtdb_get(rtdb, 'aoresponse:quad', mt_int, 1,ntmp)) 
     &    call errquit('aor_driver: lquad failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lquad = .true.

      ntmp = 1
      if (.not.rtdb_get(rtdb, 'aoresponse:damping', mt_int, 1,ntmp)) 
     &    call errquit('aor_driver: damping failed ',0, RTDB_ERR)
      if (ntmp.eq.0) lifetime = .true.

      if (lifetime) then
        if (.not.rtdb_get(rtdb, 'aoresponse:gamwidth',
     &     mt_dbl, 1,gamwidth)) 
     &     call errquit('aor_driver: gamwidth failed ',0, RTDB_ERR)
      endif        

c     note: if lgiao.eq.true. we assume that the user is
c     interested to get origin independent optical rotation which
c     we have to calculate from the magnetically perturbed density
c     matrix with the code as set up. 
c     (NOT TESTED. DON't USE. Use LBETA functionality instead).

c     if lvelocity.eq.true. we assume that NO giaos are to be used.
c     in that case the perturbing operator is imaginary no matter if
c     we use the electric or the magnetic field because it's either
c     i nabla or i/2 (r x nabla)

      if (lvelocity .and. lgiao) then
        if (oprint) then
          write (LuOut,*) 
     &       '*** WARNING: Velocity and GIAO requested => GIAO disabled'
        endif
        lgiao = .false.
      endif
      if (lgiao .and. .not.lbeta) lmagpert = .true.

c     also check if lvelocity and lbeta are both true:

      if (lvelocity .and. lbeta) then
        if (oprint) then
          write (LuOut,*) 
     &     '*** WARNING: Velocity and ORBETA requested => BETA disabled'
        endif
        lbeta = .false.
      endif

c     for the time being, if we have GIAOs set, the code
c     must go through the "beta" static B-field routine

c$$$      if (lgiao) then 
c$$$        lbeta = .true.
c$$$        lmagpert = .false.
c$$$        if (oprint) then
c$$$          write (LuOut,*) 
c$$$     &       '*** NOTE: enabled ORBETA because of GIAO keyword ***'
c$$$        end if
c$$$      end if

c     finite lifetime not yet available with the "beta" code

      if (lbeta .and. lifetime) then
        if (oprint) then
          write (LuOut,*) 
     &       '*** WARNING: LIFETIME: ORBETA or GIAO disabled ***'
          lbeta = .false.
          lgiao = .false.
        endif
      end if

c     -----------------------
c     general initializations (adapted from hnd_giaox.F)
c     -----------------------

c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb (oskel = .false.)
     
      oskel = .false.  ! temporarily set to true to enable sym.

c     If DFT get part of the exact exchange defined

      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()
c
c     Integral initialization
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)

c     read number of frequencies from run-time database

      status = rtdb_parallel(.true.)
      if (.not.rtdb_get(rtdb, 'aoresponse:nfreq', mt_int, 1,nfreq)) 
     &   then  
        call errquit('aor_driver: nfreq failed ',0, RTDB_ERR)
      endif
      if (debug) write (luout,*) 'nfreq',nfreq

c     at present we will tie the response order to the number
c     of frequencies that the user inputs, and calculate the 
c     response only for one set of frequencies. i.e. we calculate
c     alpha(-sum_i \omega_i; omega_1, omega_2, omega_3, ...).

      lstatic = .false.
      if (nfreq .eq. 0) then
        response_order = 1
        lstatic = .true.
      elseif(nfreq.gt.0) then
        response_order = nfreq
      else
        call errquit ('aor_driver: nfreq .lt. zero',nfreq,RTDB_ERR)
      end if

c     allocate memory for frequencies list and read values 
c     into dbl_mb(k_fr). If nfreq = 0, assign a single frequency = 0

      if (nfreq.ne.0) then
c       dynamic response: get frequencies
        if (.not. ma_push_get(mt_dbl,nfreq,'freqs',l_fr,k_fr)) 
     &    call errquit('aor_driver: ma_push_get failed k_fr',
     &     l_fr,MA_ERR)
        if (.not. rtdb_get(rtdb, 'aoresponse:freqlist', mt_dbl, nfreq,
     &     dbl_mb(k_fr)))
     &     call errquit('aor_driver: freqlist failed', 0, RTDB_ERR)
      else
c       static response: set all frequencies to zero
        if (.not. ma_push_get(mt_dbl,response_order,'freqs',l_fr,k_fr)) 
     &     call errquit('aor_driver: ma_push_get failed k_fr',0,MA_ERR) 
        do i=1, response_order
          dbl_mb(k_fr + (i-1)) = 0d0
        end do
      endif                     ! nfreq.ne.0

c     exit if response order is larger than one (TEMPORARY)
      if (response_order.gt.1)
     &   call errquit ('aor_driver: order>1',response_order,INPUT_ERR)

c     force nfreq = 1 if response order equals 1 (no freq. input)
      if (response_order.eq.1 .and. nfreq.ne.1) nfreq = 1

c     -------------------------------------------------------
c     if we have lbeta and lgiao set, calculate the magnetic
c     field perturbed MO vectors upfront and store in g_vecB1
c     We also compute required GIAO AO matrices used later
c     -------------------------------------------------------

      if (lbeta) then

        if (oprint) write (luout,*) 'Computing some B-field data ...'
c       retrieve nbf = number of basis functions:
        if (.not. bas_numbf(basis,nbf)) call
     &     errquit('aor_driver: could not get nbf',0, BASIS_ERR)

        alo(1) = nbf
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nbf
        ahi(2) = nbf
        ahi(3) = 3
        if (.not.nga_create(MT_DBL,3,ahi,'aor vecB1',alo,g_vecB1))
     &     call 
     &     errquit('aor_driver: nga_create failed vecB1',0,GA_ERR)
        call ga_zero(g_vecB1)

        if (lgiao) then
c         also compute
c         first-order S-ket matrix (right-hand side S-derivative)
          if (.not.nga_create(MT_DBL,3,ahi,'s-ket1 matrix',alo,g_sket1))
     &       call 
     &       errquit('aor_driver: nga_create failed g_sket1',0,GA_ERR)
          call ga_zero(g_sket1)
          call giao_b1_movecs(rtdb,basis,geom,npol, g_vecB1, g_sket1)
        else
          call fiao_b1_movecs(rtdb,basis,geom,npol, g_vecB1)
        end if
        if (oprint) write (luout,*) '... done. Reinitializing integrals'

c       the [gf]iao_b1_movecs routines terminate the integrals.
c       need to re-initialize:

        call int_init(rtdb,1,basis)
        call schwarz_init(geom,basis)
        call hnd_giao_init(basis,1)
        call scf_get_fock_param(rtdb,tol2e)

      end if ! lbeta 

c     ---------------------------------------------------------
c     write to output a summary of what we are going to do next
c     ---------------------------------------------------------

      if (oprint) write (luout,8888)
 8888 format(
     1 /,10x,47(1h-),/,
     2 10x,'Solving response equations for perturbing field',/,
     3 10x,47(1h-),/)

      if (oprint) then
        if (nfreq.ne.0) then
          write (luout,'(1x,a,1x,i4)') 'number of frequencies: ',nfreq
          if (nfreq.gt.1) then
            write (luout,'(1x,a/(e15.7))') 'frequencies in a.u.:', 
     &         (dbl_mb(k_fr+i),i=0,nfreq-1)
          else
            write (luout,'(1x,a,e15.7)') 'frequency in a.u.:', 
     &         dbl_mb(k_fr)
          endif
        else 
          write (luout,*) 'STATIC response calculation'
        endif                   ! nfreq ?= 0

        if (lmagpert) then
          write (luout,'(1x,a)') 'Perturbing field: magnetic'
        else
          write (luout,'(1x,a)') 'Perturbing field: electric'
        end if
        if (lgiao) then
          write (luout,'(1x,a)') 'Using GIAO formalism'
          write (luout,'(1x,a)')
     &       '*** WARNING, GIAO FUNCTIONALITY NOT WELL TESTED ***'
        end if
        if (lquad) then
          write (luout,'(1x,a)') 'Adding quadrupole terms to OptRot'
        end if
        if (lvelocity) then
          write (luout,'(1x,a)') 'Using Dipole Velocity Gauge'
        else
          write (luout,'(1x,a)') 'Using Dipole Length Gauge'
        endif
        if (lifetime) then
          write (luout,*)
          write (luout,'(1x,a,f10.6,1x,a)')
     &       'Applying DAMPING constant of',gamwidth,'a.u.'
        endif
        write (luout,*)
      endif                     ! oprint

c     -------------------------------------------------------
c     assemble some data neded for response calculation:
c     MO coeffs, orbital energies, density matrix, occupation
c     numbers, basis set info, matrix elements, etc.
c     -------------------------------------------------------

      if (.not. bas_numbf(basis,nbf)) call
     &   errquit('aor_driver: could not get nbf',0, BASIS_ERR)
c     allocate dbl_mb(k_occ) = fractional occupation numbers:
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('aor_driver: ma_push_get failed k_occ',0,MA_ERR)
c     allocate dbl_mb(leval) = orbital energies:
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('aor_driver: ma_push_get failed k_eval',0,MA_ERR)

c     the following call is somewhat equivalent to the ADF aoresponse
c     call to fcprep:
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,g_vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)

c     Get Unperturbed Density Matrix pointer. The following call also 
c     contains a call to hnd_prp_vec_read inside hnd_prp_get_dens:
      if (.not.lbeta .or. debug) then
        call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &     nclosed,nopen,nvirt)
      end if

      if (oprint .and. debug)
     &   write (luout,*) 'nclosed,nopen,nvirt',nclosed(1),
     &   nopen(1), nvirt(1)
      
c     Exit gracefully if we have an open-shell system:
      if (scftyp.eq.'UHF') then
        if (oprint) write(luout,9020)
        goto 7000
c       Note: ndens = 1 means closed shell
c             ndens = 3 open shell, then g_dens has 3 components
c                       with 1=alpha, 2=beta, 3=total
      endif
      
c     perform some sanity checks with the orbital occupations:
      
c     for RHF system there should be no open shell orbitals
      if (nopen(1).ne.0) call
     &   errquit('aor_driver: nopen.ne.0',0, CALC_ERR)
c     in this case, nclosed and nvirt should add up to nmo
      if ((nclosed(1)+nvirt(1)).ne.nmo) call
     &   errquit('aor_driver: wrong no. of orbitals',0, CALC_ERR)
      
c     --------------------------------------------
c     calculate one-electron integrals in AO basis
c     --------------------------------------------

      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3

c     electron dipole moment integrals:
      if (.not.nga_create(MT_DBL,3,ahi,'e-dipole',alo,g_dipel)) call 
     &    errquit('aor_driver: nga_create failed g_dipel',0,GA_ERR)
      call ga_zero(g_dipel)

c     note: origin has been set to (0,0,0) for multipole integs.
      if (lvelocity) then
        call int_giao_1ega(basis,basis,g_dipel,'velocity',origin,
     &     nat,oskel)
c       invert sign of velocity integrals for compatibility 
c       w/ dipole-length optical rotation
        call ga_scale (g_dipel, -1d0) 
      else
        call int_mpole_1ega(basis,basis,g_dipel,'dipole',origin,
     &     oskel)
      endif
c
c     DIM/QM JEM
c     Add in local field perturbations to the dipole matrix
      if (use_dimqm) then
         call dimqm_used(ldimqm)
         call dimqm_getlclfld(lclfld)
      else
         ldimqm=.false.
         lclfld=.false.
      endif
      write(luout,*) "OMEGA: ", dbl_mb(k_fr+response_order-1)
      write(luout,*) "DIMQM: ", ldimqm
      if (ldimqm) then
        call dimqm_prep(rtdb, dbl_mb(k_fr+response_order-1))
        if (lclfld) then
          call dimqm_lclfld(g_dipel, dbl_mb(k_fr+response_order-1))
        end if
      end if

c     magnetic dipole moment integrals (r x nabla), times 0.5:
      
      if (.not.lbeta .or. debug) then
        if (.not.nga_create(MT_DBL,3,ahi,'mag-dipole',alo,g_dipmag))
     &     call errquit
     &     ('aor_driver: nga_create failed g_dipmag',0,GA_ERR)
        call ga_zero(g_dipmag)
        
        if (.not.lgiao) then
          call int_giao_1ega(basis,basis,g_dipmag,'angmom',origin,
     &       nat,oskel)
          call ga_scale (g_dipmag, 0.5d0)
        end if
      end if

c     zeroth-order basis-overlap integrals:
      if (.not.nga_create(MT_DBL,2,ahi,'smat0',alo,g_smat0)) call 
     &    errquit('aor_driver: nga_create failed g_smat0',0,GA_ERR)
      call ga_zero(g_smat0)
      call int_1e_ga(basis,basis,g_smat0,'overlap',oskel)

      if (lgiao .and. .not.lbeta) then

c       -----------------------------------------------------
c       compute additional GIAO integrals if needed.
c       if lbeta is set we have already computed what we need
c       -----------------------------------------------------
        
c       magnetic-field perturbed overlap integrals, first order
        if (.not.nga_create(MT_DBL,3,ahi,'smat1 matrix',alo,g_smat1))
     &     call 
     &     errquit('aor_driver: nga_create failed g_smat1',0,GA_ERR)
        call ga_zero(g_smat1)
        call int_giao_1ega(basis,basis,g_smat1,'s10',origin,
     &     nat,oskel)
        
c       first-order S-ket matrix (right-hand side S-derivative)
        if (.not.nga_create(MT_DBL,3,ahi,'s-ket1 matrix',alo,g_sket1))
     &     call 
     &     errquit('aor_driver: nga_create failed g_sket1',0,GA_ERR)
        call ga_zero(g_sket1)
        call int_giao_1ega(basis,basis,g_sket1,'srxRb',origin,
     &     nat,oskel)
        
c       magnetic-field derivatives of <mu | h | nu> in GIAO basis,
c       use g_dipmag array for that. h is the one-electron part of
c       the Fock-operator
c       

          call ga_zero(g_dipmag)
          call int_giao_1ega(basis,basis,g_dipmag,'l10',dummy,
     &       nat,oskel)
          call int_giao_1ega(basis,basis,g_dipmag,'tv10',dummy,
     &       nat,oskel)    
 

      endif                     ! lgiao     

c     electron quadrupole moment integrals
      if (lquad) then
        alo(1) = nbf
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nbf
        ahi(2) = nbf
        ahi(3) = 6 
        if (.not.nga_create(MT_DBL,3,ahi,'e-quad',alo,g_quadel)) call 
     &     errquit('aor_driver: nga_create failed g_quadel',0,GA_ERR)
        call ga_zero(g_quadel)
        call int_mpole_1ega(basis,basis,g_quadel,'quadrupole',origin,
     &     oskel)
      end if


c     --------------
c     some debugging (may not work after 01/11 code updates)
c     --------------

      if (debug) then

c       debug: multiply density matrix with overlap matrix
c       and extract diagonal elements => N(electrons)

        if (.not.nga_create(MT_DBL,2,ahi,'temp',alo,g_temp)) call 
     &     errquit('aor_driver: nga_create failed g_temp',0,GA_ERR)

        alo(1) = 1              ! range of A, first dimension 
        ahi(1) = nbf
        alo(2) = 1              ! range of A, second dimension
        ahi(2) = nbf
        blo(1) = 1              ! range of B, first dimension, etc.
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf 
        clo(1) = 1
        chi(1) = nbf
        clo(2) = 1
        chi(2) = nbf
        call nga_matmul_patch('n','n',1d0,0d0,
     &     g_dens(1),alo,ahi,
     &     g_smat0,blo,bhi,
     &     g_temp,clo,chi)
        sum = ga_trace_diag(g_temp)
        if (oprint)
     &     write (luout,'(1x,a,3x,e15.7)') 'tr(P*S): ',sum

c       debug: calculate electronic dipole moments

c       calculate C = A*B with nga_matmul_patch
        alo(1) = 1              ! range of A, first dimension 
        ahi(1) = nbf
        alo(2) = 1              ! range of A, second dimension
        ahi(2) = nbf
        blo(1) = 1              ! range of B, first dimension, etc.
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf 
        clo(1) = 1
        chi(1) = nbf
        clo(2) = 1
        chi(2) = nbf

        do idir = 1,3

          blo(3) = idir         ! B is here the dipole matrix with
          bhi(3) = idir         ! dimension (nbf, nbf, 3)

          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_dens(1),alo,ahi,
     &       g_dipel,blo,bhi,
     &       g_temp,clo,chi)
          sum = ga_trace_diag(g_temp)
          if (oprint)
     &       write (luout,'(1x,a,i3,e15.7)') 'tr(P0*D): ',idir,sum

          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_dens(1),alo,ahi,
     &       g_dipmag,blo,bhi,
     &       g_temp,clo,chi)
          sum = ga_trace_diag(g_temp)
          if (oprint)
     &       write (luout,'(1x,a,i3,e15.7)') 'tr(P0*M): ',idir,sum
c
c         debug Smat1 and S-ket matrices
c
          if (lgiao) then
            call nga_copy_patch('t',g_sket1,blo,bhi,g_temp,clo,chi)
            call nga_add_patch(1d0,g_sket1,blo,bhi,-1d0,g_temp,clo,chi,
     &         g_temp,clo,chi)
c            if (oprint) then
c              call nga_print_patch(g_smat1,blo,bhi,1)
c              call nga_print_patch(g_temp,clo,chi,1)
c            endif
            call nga_add_patch(1d0,g_smat1,blo,bhi,-1d0,g_temp,clo,chi,
     &         g_temp,clo,chi)
            if (oprint) then
              write (6,*)
     &         '(S-ket - transpose(S-ket)) - S10 in AO bas idir =',idir
              call ga_print(g_temp)
            endif
            
          endif                 ! lgiao
          
        enddo                   ! idir = 1,3
        
        if (.not.ga_destroy(g_temp)) call 
     &     errquit('aor_driver: ga_destroy failed g_temp',0,GA_ERR)
        
      end if                    ! debug

c     ------------------------------------------------
c     All necessary input data has now been assembled.
c     Allocate memory for temporary arrays and
c     call aoresponse routine to carry out the CPKS
c     and determine the response property:
c     ------------------------------------------------

c     Here is what we need as input to do the response. GA indicates
c     that we don't actually have the array here but instead we have
c     global arrays objects. MA indicates a memory allocator array.

c     number of frequencies       :   nfreq
c     response order              :   response_order
c     frequencies list            :   dbl_mb(k_fr) (MA)
c     unperturbed overlap matrix  :   g_smat0 (GA)
c     electric dipole integrals   :   g_dipel (GA)
c     magnetic dipole integrals   :   g_dipmag (GA)
c     unperturbed MO coeffs       :   g_vectors(2) (GA)
c     MO energies                 :   dbl_mb(k_eval) (MA)
c     MO occupations              :   dbl_mb(k_occ) (MA)
c     unperturbed density matrix  :   g_dens(3) (GA) 
c     number of basis functions   :   nbf    
c     number of MOs               :   nmo
c     number of closed shell MOs  :   nclosed(2) = occ MOs
c     number of open shell MOs    :   nopen(2) should be zero
c     number of empty MOs         :   nvir(2) = virt MOs

      nbfx = nbf*(nbf+1)/2 ! not needed at present

c     allocate temporary array to hold orbital energy difference
c     denominators:
      nocv = nclosed(1)*nvirt(1)
      if (.not. ma_push_get(mt_dbl,nocv,'deps_tmp',l_deps,k_deps)) 
     &   call errquit('aor_driver: ma_push_get failed k_deps',
     &   l_deps,MA_ERR)

      do idir = 1,3
        do ipert = 1,3
          alfare(idir,ipert) = zero ! polarizability
          alfaim(idir,ipert) = zero
          betare(idir,ipert) = zero ! optical rotation
          betaim(idir,ipert) = zero
        enddo
      enddo
      
      if (.not.lbeta) then
c     default call: compute dynamic polarizability and
c     optical rotation G' tensor
         if (debug) write (luout,*) 'calling aoresponse'
         call aoresponse (rtdb, basis, geom, nfreq, response_order,
     &        dbl_mb(k_fr),
     &        g_smat0, g_smat1, g_sket1,
     &        g_dipel, g_dipmag, g_vectors, dbl_mb(k_eval),
     &        dbl_mb(k_occ), g_dens, nbf, nbfx, nmo, nclosed, nvirt,
     &        dbl_mb(k_deps), lgprime, lgiao, lvelocity, lmagpert,
     &        lifetime, gamwidth, alfare, alfaim, betare, betaim)
         if (debug) write (luout,*) 'exiting aoresponse'
         
      else ! lbeta?

c     alternate call: compute directly optical rotation beta tensor
c     using frequency-dependent generalization of Amos' formula

         if (debug) write (luout,*) 'calling aoresponse_beta'

c        original aor_beta where the magnetically perturbed 
c        vectors are are generated within the routine:
c         call aoresponse_beta (rtdb,basis,geom, nfreq, response_order,
c     &        dbl_mb(k_fr),
c     &        g_smat0, g_smat1, g_sket1,
c     &        g_dipel, g_quadel, g_dipmag, g_vectors, dbl_mb(k_eval),
c     &        dbl_mb(k_occ), g_dens, nbf, nbfx, nmo, nclosed, nvirt,
c     &        dbl_mb(k_deps), lgiao, lquad, lvelocity, lmagpert,
c     &        lifetime, gamwidth, alfare, alfaim, betare, betaim)

c        new aor_beta where we use the g_vecB1 data
c        calculated earlier:
         call aoresponse_beta2 (rtdb, basis, geom, nfreq,
     &      response_order,
     &      dbl_mb(k_fr),
     &      g_smat0, g_sket1, g_vecB1,
     &      g_dipel, g_quadel, g_vectors, 
     &      dbl_mb(k_occ), nbf, nmo, nclosed, nvirt,
     &      dbl_mb(k_deps), lgiao, lquad, lvelocity, lmagpert,
     &      lifetime, gamwidth, alfare, alfaim, betare, betaim)
         if (debug) write (luout,*) 'exiting aoresponse_beta'

      end if ! lbeta

c     ------------
c     print result
c     ------------

      if (oprint) then
        write (LuOut,*)
     &     'Electric Dipole Response Matrix (nonzero elements):'
        if (lmagpert) then
          if (lstatic) then
            write (LuOut,*) 'Optical rotation tensor G`'
          else
            write (LuOut,*) 'Optical rotation tensor Beta'
          endif
        endif
        call output(alfare,1,3,1,3,3,3,1)
      endif
c
          do idir = 1,3
            do ipert = 1,3
              alfare2(idir,ipert) = alfare(idir,ipert)
              alfaim2(idir,ipert) = alfaim(idir,ipert)
            enddo
          enddo
c
c         diagonalize
          call hnd_diag(alfare2,alphaevr,3,.false.,.false.)  ! replaces dgeev call
c
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0 ! isotropic
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(alphaani/2.d0)  ! anisotropic
c
          if (oprint) then
            if (dabs(dbl_mb(k_fr)).gt.tolfr) then
              write(LuOut,9434) "DFT Linear Response",
     1              dbl_mb(k_fr),auXnm/dbl_mb(k_fr),
     2              alfare(1,1),alfare(1,2),alfare(1,3),
     3              alfare(2,1),alfare(2,2),alfare(2,3),
     4              alfare(3,1),alfare(3,2),alfare(3,3),
     5              alphaevr(1),alphaevr(2),alphaevr(3),
     6              alphaiso,alphaani
            else
              write(LuOut,9435) "DFT Linear Response",
     1              dbl_mb(k_fr),
     2              alfare(1,1),alfare(1,2),alfare(1,3),
     3              alfare(2,1),alfare(2,2),alfare(2,3),
     4              alfare(3,1),alfare(3,2),alfare(3,3),
     5              alphaevr(1),alphaevr(2),alphaevr(3),
     6              alphaiso,alphaani
            endif
          endif
          call util_flush(LuOut)
c
c        alpha_rtdb = [omega, XX, XY, XZ, YY, YZ, ZZ, LL, MM, NN, Iso, Ani]
c
          alpha_rtdb(1)  = dbl_mb(k_fr)
          alpha_rtdb(2)  = alfare(1,1)
          alpha_rtdb(3)  = alfare(1,2)
          alpha_rtdb(4)  = alfare(1,3)
          alpha_rtdb(5)  = alfare(2,2)
          alpha_rtdb(6)  = alfare(2,3)
          alpha_rtdb(7)  = alfare(3,3)
          alpha_rtdb(8)  = alphaevr(1)
          alpha_rtdb(9)  = alphaevr(2)
          alpha_rtdb(10) = alphaevr(3)
          alpha_rtdb(11) = alphaiso
          alpha_rtdb(12) = alphaani
c
c         Write polarizability information to the RTDB
c
          if (.not. rtdb_put(rtdb, 'aoresponse:alpha', mt_dbl,12,
     1                       alpha_rtdb))
     2       call errquit('aoresponse_driver: rtdb_put failed',12,
     3                    RTDB_ERR)
c
ccc#endif
c
c
c
      if (oprint) then
        if (lifetime) then
          write (luout,*)
          write (luout,*) 'Imaginary part:'
          call output(alfaim,1,3,1,3,3,3,1)
        endif
        write (LuOut,*)
        sum = 0d0
        sumim = 0d0
        do idir = 1,3
          sum = sum + alfare(idir,idir)
          sumim = sumim + alfaim(idir,idir)
        enddo
        write(luout,*) 'average: ',sum/3d0,' + I ',sumim/3d0
        write (luout,*)
        write (LuOut,*)
     &     'Magnetic Dipole Response Matrix (nonzero elements):'
        if (.not.lmagpert) then
          if ((lstatic .and. .not.lbeta) .or. lgprime) then
            write (LuOut,*) 'Optical rotation tensor -G'' '
          else
            write (LuOut,*) 'Optical rotation tensor Beta'
          endif
        endif
        if (lstatic .and. .not.lbeta) write(luout,*)
     &  '*** static G'' is incorrect. use GPRIME input for small omega'
        call output(betare,1,3,1,3,3,3,1)
        if (lifetime) then
          write (luout,*)
          write (luout,*) 'Imaginary part:'
          call output(betaim,1,3,1,3,3,3,1)
        endif
        sum = 0d0
        sumim = 0d0
        do idir = 1,3
          sum = sum + betare(idir,idir)
          sumim = sumim + betaim(idir,idir)
        enddo
        write(luout,*) 'average: ',sum/3d0, ' + I ',sumim/3d0
        write (LuOut,*)
      endif
c
c     DIM/QM JEM
c     Calculate the DIM polarizability
      if(ldimqm) then
        call dimqm_polarizability(rtdb, dbl_mb(k_fr))
      endif

c
c     == Raman bits ==
c
      if (.not. rtdb_put(rtdb, 'raman:alfare ', mt_dbl, 9, alfare))   ! real alpha
     &  call errquit('aoresponse:failed to write alfare', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'raman:alfaim ', mt_dbl, 9, alfaim))   ! imag alpha
     &  call errquit('aoresponse:failed to write alfaim', 0, RTDB_ERR)
c
c     ----------------------------------------
c     clean up, print exit statment and return
c     ----------------------------------------

 7000 continue ! jump here in case of abnormal termination

c     deallocate global arrays

      if (.not.ga_destroy(g_vectors)) call 
     &    errquit('aor_driver: ga_destroy failed g_vectors',0,GA_ERR)
      if (.not.lbeta .or. debug) then
        do i=1,ndens
          if (.not.ga_destroy(g_dens(i))) call 
     &       errquit('aor_driver: ga_destroy failed g_dens',0,GA_ERR)
        enddo
      end if
      if (.not.ga_destroy(g_dipel)) call 
     &    errquit('aor_driver: ga_destroy failed g_dipel',0,GA_ERR)

      if (.not.lbeta .or. debug) then
        if (.not.ga_destroy(g_dipmag)) call 
     &     errquit('aor_driver: ga_destroy failed g_dipmag',0,GA_ERR)
      end if

      if (.not.ga_destroy(g_smat0)) call 
     &    errquit('aor_driver: ga_destroy failed g_smat0',0,GA_ERR)

      if (lgiao  .and. .not.lbeta) then
        if (.not.ga_destroy(g_smat1)) call 
     &     errquit('aor_driver: ga_destroy failed g_smat1',0,GA_ERR)
      end if

      if (lgiao) then
        if (.not.ga_destroy(g_sket1)) call 
     &     errquit('aor_driver: ga_destroy failed g_sket1',0,GA_ERR)
      endif

      if (lbeta) then
        if (.not.ga_destroy(g_vecB1)) call 
     &     errquit('aor_driver: ga_destroy failed g_vecB1',0,GA_ERR)
      end if

      if (lquad) then
        if (.not.ga_destroy(g_quadel)) call 
     &     errquit('aor_driver: ga_destroy failed g_quadel',0,GA_ERR)
      end if

      call ga_sync()


c     ------------------------------------------------------------
c     deallocate all MA arrays allocated up to and including the
c     frequencies list which was the first array that we allocated
c     ------------------------------------------------------------

      if (.not.ma_chop_stack(l_fr)) call
     &   errquit('aor_driver: ma_chop_stack failed k_fr',l_fr,MA_ERR)

c      call schwarz_tidy()
c      call int_terminate()
      
      if (oprint) then
        write (LuOut,*) 'Exiting AOResponse driver routine'
      endif

c     ==================================================================

      return

 9000 format(//80(1h-)//,
     &   10x,16(1h*),/,10x,'*** RESPONSE ***',/,10x,16(1h*)/)
 9010 format(' Response module for NWChem and dynamic CPKS solver'/
     &       ' developed by J. Autschbach and coworkers, SUNY Buffalo'/
     &       ' The methodology used in this program is described in '/
     &       ' J. Chem. Phys. 123 (2005), 114103'/
     &       ' J. Chem. Phys. 122 (2005), 224115'/
     &       ' J. Chem. Phys. 122 (2005), 074105'/
     &       ' Comp. Lett. 3 (2007), 131-150 (contact JA for a copy)'/
     &       ' Please cite this work in publications based on results'/
     &       ' obtained with this code. Thank you!'/)
 9020 format(/,10x,'Open shell calculations not yet implemented',
     &      ' for the Response module. Exiting')
 9434 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength = ',f15.7,' / nm',/
     2  1x,'             X              Y              Z',/
     3  1x,'-----------------------------------------------',/
     4  1x,'X ',3f15.7,/
     5  1x,'Y ',3f15.7,/
     6  1x,'Z ',3f15.7,/
     3  1x,'-----------------------------------------------',/
     6  1x,'Eigenvalues = ',3f15.7,/
     6  1x,'Isotropic   = ',1f15.7,/
     6  1x,'Anisotropic = ',1f15.7,/
     7  1x,'-----------------------------------------------')
 9435 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength =        Infinity / nm',/
     2  1x,'             X              Y              Z',/
     3  1x,'-----------------------------------------------',/
     4  1x,'X ',3f15.7,/
     5  1x,'Y ',3f15.7,/
     6  1x,'Z ',3f15.7,/
     3  1x,'-----------------------------------------------',/
     6  1x,'Eigenvalues = ',3f15.7,/
     6  1x,'Isotropic   = ',1f15.7,/
     6  1x,'Anisotropic = ',1f15.7,/
     7  1x,'-----------------------------------------------')
      end
