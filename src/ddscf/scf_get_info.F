      block data cscfdata
      implicit none
#include "cscf.fh"
c
      data oinitialized /.false./
      data g_movecs /0,0/
      data k_occ /10000000/
      data k_eval /10000000/
c
      end
      subroutine scf_init(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
#include "cosmo.fh"
c     
      integer rtdb              ! database handle
c     
      double precision nuclear_charge
      character*255 name
      integer len_occ
      external cscfdata ! For T3D linker
c
      logical osome
c
      logical hf_job
      character*30 tag
      character*255 theory
      integer mult
      logical xc_gotxc,xc_gothfx
      external xc_gotxc,xc_gothfx
c
      hf_job = (.not. xc_gotxc()).and.(.not. xc_gothfx())
c
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
c
c     load geometry and symmetry info
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0, RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif
      if (.not.rtdb_get(rtdb, 'scf:adapt',MT_LOG, 1, oadapt)) then
         oadapt = sym_number_ops(geom) .gt. 0
      endif
      if (.not.rtdb_get(rtdb, 'scf:lock',MT_LOG, 1, olock)) then
         olock = .false.
      endif
c     
c     load the basis set and get info about it
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $        call errquit('scf_init: no ao basis set', 0, RTDB_ERR)
c
c     For debug ... call int_init and do the 2-e
c
      if (util_print('texas init debug',print_never)) then
        call int_init(rtdb, 1, basis)
        write(6,*) ' DONE INIT'
*        call schwarz_init(geom, basis)
*        call schwarz_tidy()
        call int_terminate
      endif
c
      if (.not. bas_name(basis, name, trans))
     $     call errquit('scf_init: bas_name?', 0, BASIS_ERR)
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
c
c     Is RI approximation to be used? If so get fitting basis set.
c
      if (rtdb_get(rtdb, 'scf:ri', MT_INT, 1, nriscf)) then
        if (.not. bas_create(riscf_basis, 'riscf basis'))
     $       call errquit('scf_init: bas_create?', 0, BASIS_ERR)
        if (.not. bas_rtdb_load(rtdb, geom, riscf_basis, 'riscf basis'))
     $       call errquit('scf_init: no riscf basis set', 0, RTDB_ERR)
      else
        nriscf = 0
        riscf_basis = 0
      endif
c
c     Figure input/output MO vectors
c     
      if (hf_job) then
         tag = 'scf:input vectors'
      else
         tag = 'dft:input vectors'
      endif
      if (.not. rtdb_cget(rtdb, tag, 1, movecs_in))
     $     movecs_in = 'atomic'
      if (hf_job) then
         tag = 'scf:output vectors'
      else
         tag = 'dft:output vectors'
      endif
      if (.not. rtdb_cget(rtdb, tag, 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
     $        movecs_in.eq.'project' .or. movecs_in.eq.'fragment'.or.
     $        movecs_in.eq.'molden'
     $          .or.movecs_in.eq.'rotate') then
            call util_file_name('movecs', .false.,.false.,movecs_out)
         else
            movecs_out = movecs_in
         endif
      endif
c
c     Resolve names of MO files to full paths defaulting to the
c     permanent directory
c
      if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
     $        movecs_in.eq.'project' .or. movecs_in.eq.'fragment'.or.
     $        movecs_in.eq.'molden'
     $          .or.movecs_in.eq.'rotate') then
         continue
      else
         call util_file_name_resolve(movecs_in, .false.)
      endif
      call util_file_name_resolve(movecs_out, .false.)
c     
c     Figure out the number of electrons from the required total
c     charge and the sum of nuclear charges
c     
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0, GEOM_ERR)
      nelec = nint(nuclear_charge - charge)
      if (nelec .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec, INPUT_ERR)
      if (abs(nuclear_charge - charge - dble(nelec)) .gt. 1d-8)
     $     call errquit('scf: non-integral no. of electrons ?', 0,
     &       INPUT_ERR)
c     
c     Determine no. of open and closed shells ... default is to run closed
c     shell unless told otherwise
c     
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     +     call errquit('task: no task input for theory?',0, INPUT_ERR)
c     if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
c     $     nopen = 0
      if(theory .eq. 'dft')then
         if (.not. rtdb_get(rtdb, 'dft:mult', MT_INT, 1,mult))
     *        mult = 1
         nopen = mult - 1
      endif
      if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     nopen = 0
      if (nopen .gt. nelec) call errquit
     $     ('scf: nopen > nelec ', nopen, INPUT_ERR)
      if (mod(nelec-nopen,2) .ne. 0) call errquit
     $     ('scf: no. of closed-shell electrons is not even!',nopen,
     &       INPUT_ERR)
      nclosed = (nelec-nopen) / 2
c     
      if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
         if (nopen .eq. 0) then
            scftype = 'RHF'
         else
            scftype = 'ROHF'
         endif
      endif
c
      call inp_ucase(scftype)
c
c     Take care of holes in the input routines
c
      if (scftype.eq.'RHF' .and. nopen.gt.0) then
         scftype = 'ROHF'
      else if (scftype.eq.'ROHF' .and. nopen.eq.0) then
         scftype = 'RHF'
      endif
c
      if ( scftype.ne.'ROHF' .and. scftype.ne.'RHF' .and.
     $     scftype.ne.'UHF' ) call errquit
     $     ('scf: only ROHF, RHF, and UHF currently supported', 0,
     &       INPUT_ERR)
c
c     Dump lagrangian?  Yes by default now since if the SCF 
c     has converged for an energy it will not be rerun for the gradient.
c
      if (.not.rtdb_get(rtdb, 'scf:lagrangian',MT_LOG, 1, olagr))
     $   olagr = scftype.eq.'ROHF'
c
      nalpha = nclosed + nopen
      nbeta  = nclosed
c
c   DIIS toggle
c
      if (.not.rtdb_get(rtdb, 'scf:diis',MT_LOG, 1, odiis))
     $   odiis = .false.
c
      call ga_sync()
c
c     For now set NMO = NBF, however this may change later when the 
c     linear dependency analysis is done just before the starting guess
c
      nmo = nbf
c
c     Store the derived info in the database for other wavefunction
c     modules and/or restart to access
c
      if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, scftype)) 
     $     call errquit('scf_init: put of scftyp failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     call errquit('scf_init: put of nopen failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit('scf_init: put of nclosed failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:nelec', MT_INT, 1, nelec))
     $     call errquit('scf_init: put of nelec failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:nmo', MT_INT, 1, nmo))
     $     call errquit('scf_init: put of nmo failed',0, RTDB_ERR)
      if (scftype .eq. 'UHF') then
         if (.not. rtdb_put(rtdb, 'scf:nalpha', MT_INT, 1, nalpha))
     $        call errquit('scf_init: put of nalpha failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'scf:nbeta', MT_INT, 1, nbeta))
     $        call errquit('scf_init: put of nbeta failed',0, RTDB_ERR)
      endif
c     
c     Allocate persistent local and global arrays ... these may
c     be reallocated later when the dependency analysis is done
c
      if (scftype .eq. 'UHF') then
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_init: alpha MOs',
     $        32, 32, g_movecs)) call errquit('scf_init: alpha MOs', 0,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_init: beta MOs',
     $        32, 32, g_movecs(2))) call errquit('scf_init: beta MOs',0,
     &       GA_ERR)
      else
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_init: MOs',
     $        32, 32, g_movecs)) call errquit('scf_init: MOs', 0,
     &       GA_ERR)
      endif
c     
      len_occ = nmo
      if (scftype .eq. 'UHF') len_occ = nbf * 2
      if (.not. ma_push_get(mt_dbl, len_occ, 'scf_init: mo evals',
     $     l_eval, k_eval)) call errquit
     $     ('scf_init: insufficient memory?', len_occ, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, len_occ, 'scf_init: mo occ',
     $     l_occ, k_occ)) call errquit
     $     ('scf_init: insufficient memory?', len_occ, MA_ERR)
c
      if (.not. ma_push_get(mt_int, len_occ, 'scf_init: mo irs',
     $     l_irs, k_irs)) call errquit
     $     ('scf_init: insufficient memory?', len_occ, MA_ERR)
c
      call ifill(len_occ, 1, int_mb(k_irs), 1) ! In case not adapting
c
c     Fill in the SCF convergence info
c     
      call scf_get_conv_info(rtdb)
      call scf_get_fock_param(rtdb, tol2e)
c
c     ----- cosmo initialization ----
c
      cosmo_last = .false.
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
         if(cosmo_on) then
            osome = util_print('information', print_low)
c
            call cosmo_initialize(rtdb,geom,basis,osome)
c
c           Turn cosmo on, we want to run the calculation
c           Start with gas_phase run
c
            cosmo_last = .true.
            cosmo_on = .true.
            if(.not.rtdb_get(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
     >         cosmo_phase = 1
         endif
      endif
      oinitialized = .true.
c     
      end
      subroutine scf_get_conv_info(rtdb)
C$Id$
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "cscf.fh"
c     
      integer rtdb
c
      double precision mone
      parameter (mone = -1.0d0)
c     
      ouser_changed_conv = .false.
c
      if (.not. rtdb_get(rtdb, 'scf:maxiter', MT_INT, 1, maxiter)) then
         if (scftype .ne. 'UHF') then
            maxiter = 30
         else
            maxiter = 30
         endif
      endif
      if (.not.rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, gnorm_tol)) 
     $     gnorm_tol = 1.d-4
c
c     Ensure that the default integral selection is sufficient
c     for the request accuracy of the SCF.  However, allow user override.
c
      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) 
     $     tol2e = min(1.0d-7,gnorm_tol*1d-2)
c
      if (rtdb_get(rtdb, 'scf:level shift info', MT_DBL, 6,shifts))then
         ouser_changed_conv = .true.
      else
         call dfill(6, mone, shifts, 1)
      endif
      if (rtdb_get(rtdb, 'scf:full hessian switch', MT_DBL, 1, 
     $     nr_gswitch)) then
         ouser_changed_conv = .true.
      else 
         nr_gswitch = 0.1d0
      endif
c
c     Apply defaults
c
      if (shifts(1) .eq. -1.0d0) shifts(1) = 5.0d0
      if (shifts(2) .eq. -1.0d0) shifts(2) = 0.5d0
      if (shifts(3) .eq. -1.0d0) shifts(3) = 0.0d0
      if (shifts(4) .eq. -1.0d0) shifts(4) = 0.0d0
      if (shifts(5) .eq. -1.0d0) shifts(5) = 0.0d0
      if (shifts(6) .eq. -1.0d0) shifts(6) = 0.0d0
c
      end
      subroutine scf_get_fock_param(rtdb, tol2e)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "cfock.fh"
#include "global.fh"
#include "util.fh"
      integer rtdb
      double precision tol2e
c
      logical orestart
      logical int2e_file_open
      external int2e_file_open
c
c     Load balancing information
c
      if (.not. rtdb_get(rtdb, 'fock:task_bf', MT_INT, 1, task_bf))
     $     task_bf = -1
c
c     Dentolmax constrains dentol=tol2e/denmax which is used
c     to screen integrals before screening integrals*density
c     against tol2e.  Making dentolmax smaller will, up to a
c     point, make the hessian-vector products more stable,
c     though slower.  In strongly convergent cases dentolmax
c     could be increased to perhaps 1d-3 with some speed gain.
c     Dentolmax is still used even if density
c     screening is turned off since otherwise the hessian 
c     vector products would be disastrously expensive.
c
c     C4H10 shows worse convergence with dentolmax > 1d-5
c     C60 shows worsse convergence with dentolmax > 1d-6
c
c     Some poor convergence in the CPHF and line search can
c     be attributed to dentolmax being too high when preconditioning
c     with Hessian-vector products.  I think there must be a problem
c     in the screening algorithm somewhere ... but where?
c
      if (.not. rtdb_get(rtdb, 'fock:dentolmax', mt_dbl, 1, dentolmax))
     $     dentolmax = 1d-6
c
c     odensityscreen enables/disables density screening inside
c     the fock build.  Turning it off can increase stability but will
c     slow things down.
c
      if (.not. rtdb_get(rtdb, 'fock:densityscreen', mt_log, 1, 
     $     odensityscreen)) odensityscreen = .true.
c
c     Integral caching/file information
c
#if defined(CRAYXT) || defined(BGP) || defined(BGQ) || defined(NOIO)
c     direct by default to avoid io
      if (.not. rtdb_get(rtdb, 'int2e:filesize',
     $     MT_INT, 1, filesize)) filesize = -1
      if (.not. rtdb_get(rtdb, 'int2e:memsize',
     $     MT_INT, 1, memsize)) memsize = -1
#else
      if (.not. rtdb_get(rtdb, 'int2e:filesize',
     $     MT_INT, 1, filesize)) filesize = 0
      if (.not. rtdb_get(rtdb, 'int2e:memsize',
     $     MT_INT, 1, memsize)) memsize = 0
#endif
      if (.not. rtdb_get(rtdb, 'int2e:restart',
     $     MT_LOG, 1, orestart)) orestart = .false.
c
c     The opening routine will put the .pid on the integral filename
c     (hence even tho' parallel file, open as sequential)
c
      if (.not. rtdb_cget(rtdb, 'int2e:filename', 1, int2efilename)) 
     $     call util_file_name('aoints',.true.,.false.,int2efilename)
c
      oreadfile = .false.
      owritefile= .false.
c
      if (filesize.gt.0 .or. memsize.gt.0) then
         if (.not. int2e_file_open(int2efilename, memsize, filesize,
     $        tol2e*0.1d0, orestart)) then
c     Probably memsize/filesize too small to be useful ... just force direct
            filesize = -1
            memsize = -1
         else
c     Managed to open a file or memory cache
            if (orestart) then
               oreadfile = .true.
            else
               owritefile = .true.
            endif
         endif
      endif
c
c     Blocking integral interface
c
      if (.not. rtdb_get(rtdb, 'fock:maxquartet',
     $     MT_INT, 1, maxquartet)) maxquartet = 10000
      if (.not. rtdb_get(rtdb, 'fock:maxpair',
     $     MT_INT, 1, maxquartet)) maxpair = sqrt(dble(maxquartet))
      if (.not. rtdb_get(rtdb, 'fock:maxeri',
     $     MT_INT, 1, maxeri)) maxeri = 1296*100 ! 100 d(6) quartets
      if (.not. rtdb_get(rtdb, 'fock:maxscr',
     $     MT_INT, 1, maxscr)) maxscr = 0 ! <=0 implies use texas estimate
      if (.not. rtdb_get(rtdb, 'fock:intacc',
     $     MT_DBL, 1, intacc)) intacc = 0.0d0 ! Default is variable
c
      if (.not. rtdb_get(rtdb, 'fock:replicated',
     $     mt_log, 1, oreplicated)) oreplicated = .true.
      if (.not. rtdb_get(rtdb, 'fock:noblock',
     $     mt_log, 1, onoblock)) onoblock = .false.
c
c  Use labels by default
c
      if (.not. rtdb_get( rtdb, 'fock:uselabels',  MT_LOG, 1,
     $                    oerilabel)) oerilabel = .true.
c
      if (ga_nodeid().eq.0 .and.
     $     util_print('fock param',print_high)) then
         write(6,2) task_bf, maxquartet, maxeri, maxscr, intacc,
     $        odensityscreen, dentolmax
 2       format(/
     $        ' Setting fock-build task_bf   :', i8/
     $        '                    maxquartet:', i8/
     $        '                    maxeri    :', i8/
     $        '                    maxscr    :', i8/
     $        '                    intacc    :', 1p,d8.1/
     $        '                    denscreen : ',l1/
     $        '                    dentol    :', d8.1/)
         call util_flush(6)
      endif
c
      end
      subroutine fock_force_direct(rtdb)
      implicit none
#include "cfock.fh"
      integer rtdb
c
c     Call after scf_get_fock_param to force a fock build
c     to be direct
c
      call fock_2e_tidy(rtdb)
c
      filesize = -1
      memsize  = -1
c
      end
      subroutine scf_tidy(rtdb)
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "cfock.fh"
#include "rtdb.fh"
#include "cosmo.fh"
      integer rtdb
c
      logical status
c
      if (oinitialized) then
         if (.not. geom_destroy(geom)) call errquit
     $        ('scf_tidy: geom destroy failed', 0, GEOM_ERR)
         if (.not. bas_destroy(basis)) call errquit
     $        ('scf_tidy: basis destroy failed',0, BASIS_ERR)
         status = ma_pop_stack(l_irs)
         status = ma_pop_stack(l_occ) .and. status
         status = ma_pop_stack(l_eval) .and. status
         if (.not. status) call errquit
     $        ('scf_tidy: failed to free irs/occupation/evals',0,
     &       MA_ERR)
         if (.not. ga_destroy(g_movecs)) call errquit
     $        ('scf_tidy: failed to free movecs',0, GA_ERR)
         if (scftype .eq. 'UHF') then
            if (.not. ga_destroy(g_movecs(2))) call errquit
     $           ('scf_tidy: failed to free beta movecs',0, GA_ERR)
         endif
         oinitialized = .false.
      endif
c
      call fock_2e_tidy(rtdb)
c
c     ----- cosmo cleanup and reset -----
c
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
         if(cosmo_on) then
           call cosmo_tidy(rtdb)
           cosmo_on = .false.
           cosmo_phase = 1
         endif
      endif
c
      end
      subroutine fock_2e_tidy(rtdb)
      implicit none
#include "errquit.fh"
#include "cfock.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
c
      logical int2e_file_close
      external int2e_file_close
      logical okeep
c
      if (.not. rtdb_get(rtdb, 'int2e:keep',
     $     MT_LOG, 1, okeep)) okeep = .false.
      if (.not. int2e_file_close(okeep))
     $     call errquit('scf_tidy: closing aoints?', 0, INT_ERR)
c
      end
