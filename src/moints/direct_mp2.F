      logical function direct_mp2(rtdb)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "moints.fh"
c
      integer rtdb                     ! [input] database handle
      integer geom, basis
      character*255 title
      character*255 movecs_in
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      integer nbf_vecs, nsets
      integer k_occ, l_occ, k_evals, l_evals, g_vecs
      integer k_epair, l_epair
      integer npairs
      integer nbf, nsh, maxbfsh
      integer nelec,  nocc, nopen, nocc_correl
      integer occ_frozen, occ_lo, occ_hi
      integer vir_frozen, vir_lo, vir_hi, nmoarray(2), nmo, nmo2
      integer noper, oseg_lo, oseg_hi, npass, blen, nocct
      integer g_coul, g_exch, g_epair
      integer i, j, ij, ii
      logical oprint, osym, status, omemprint
      logical ohalf, oblk, odisk, oaoreuse, oaokeep
      logical dftmp2            ! use DFT reference orbitals and energy
      double precision mp2fac ! coefficient on MP2 contribution to total energy
      double precision emp2, e2pair
      double precision tpass, ttotal, scf_energy
      double precision tol2e
      character*7 algo
      character*4 scftype
      double precision fss, fos ! Spin Component Scaled (SCS) MP2
      logical scs               ! SCS-MP2 flag
c
      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ga_create_JKblocked
      external ga_create_JKblocked
      logical scf, dft_energy
      external scf, dft_energy
      logical int_normalize
      external int_normalize
c
      data g_coul,g_exch/-1,-1/
      data ohalf/.false./
c
      if (.not. rtdb_get(rtdb,'dft:dftmp2',mt_log,1,dftmp2))
     $     dftmp2 = .false.
c
c     Check if mp2 is requested with dft orbitals
      if (dftmp2) then
        if (.not.dft_energy(rtdb)) then
           direct_mp2 = .false.
           return
        endif
      else
        if (.not. scf(rtdb)) then
           direct_mp2 = .false.
           return
        endif
      endif
c
      direct_mp2 = .true.
c
      ttotal = util_cpusec()
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call ecce_print_module_entry('mp2')
c
c Print settings
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'mp2')
      oprint = util_print('information',print_default)
      omemprint = util_print('memory',print_high)
*c
*c Memory settings ... RJH ... don't set these here ... leave under user control
*c
*      status = ma_set_hard_fail(.true.)
*      status = ma_set_auto_verify(.true.)
c
c Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $      title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $      call errquit('direct_mp2: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $      call errquit('direct_mp2: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $      call errquit('direct_mp2: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $      call errquit('direct_mp2: no ao basis ', 0, RTDB_ERR)
      if (.not. rtdb_cget(rtdb, 'input mo vector file', 1,
     $      movecs_in)) call util_file_name
     $      ('movecs',.false.,.false.,movecs_in)
      if (.not. bas_numbf(basis, nbf))
     $      call errquit('direct_mp2: cannot get basis info',0,
     &      BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $      call errquit('direct_mp2: cannot get basis info',0,
     &      BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,maxbfsh))
     $      call errquit('direct_mp2: cannot get basis info',0,
     &      BASIS_ERR)
      if (dftmp2) then
        if (.not. rtdb_get(rtdb, 'dft:energy', MT_DBL, 1, scf_energy))
     $        scf_energy = 0.0d0
        if (.not.rtdb_get(rtdb,'dft:mp2fac',mt_dbl,1,mp2fac)) then
          mp2fac = 1.0d0
          if (GA_NodeID().eq.0) then
            write(6,*) 'dft:mp2fac not found, using 1.0'
          endif
        endif
      else
        if (.not. rtdb_get(rtdb, 'scf:energy', MT_DBL, 1, scf_energy))
     $        scf_energy = 0.0d0
      endif
      if (.not. rtdb_get(rtdb, 'mp2:aotol2e', MT_DBL, 1, tol2e))
     $     tol2e = 1d-9
c
c Get algorithm choice
c
      if (.not. rtdb_cget(rtdb, 'fourindex:method ', 1, algo))
     $      algo = 'twofold'
c
c Check disk caching toggles
c
      if (.not.rtdb_get(rtdb,'fourindex:aodisk', MT_LOG,
     $                 1, odisk)) odisk = .false.
      if (.not.rtdb_get(rtdb,'fourindex:aoreuse', MT_LOG,
     $                 1, oaoreuse)) oaoreuse = .false.
      if (.not.rtdb_get(rtdb,'fourindex:aokeep', MT_LOG,
     $                 1, oaokeep)) oaokeep = .false.
c
c Reconcile number of electrons and occupied/virtual orbitals
c
      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nocc))
     $      call errquit('direct_mp2: get of nocc failed',0, RTDB_ERR)
c
c     looking at scftype is more general than nopen
c
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
     $     call errquit('scf:scftype missing',1, RTDB_ERR)
      if (scftype.ne.'RHF') 
     $      call errquit('direct_mp2: scftype not RHF',911,INPUT_ERR)
#ifdef UNUSED
      if (dftmp2) then
        if (.not. rtdb_get(rtdb, 'dft:nopen', MT_INT, 1, nopen))
     1        call errquit('direct_mp2: get of nopen failed',0,
     2                     RTDB_ERR)
        if (nopen.ne.0) 
     1        call errquit('direct_mp2: nopen is not zero',nopen,
     2        INPUT_ERR)
      else
        if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     1        call errquit('direct_mp2: get of nopen failed',0,
     2                     RTDB_ERR)
        if (nopen.ne.0) 
     1        call errquit('direct_mp2: nopen is not zero',nopen,
     2        INPUT_ERR)
      endif
#endif
c
      nelec = 2*nocc
      status = ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ)
c
c     Look for MP2 specific MO vectors.  If not there use previous SCF.  
c     If not there use default.
c
      if (dftmp2) then
        if (.not.rtdb_cget(rtdb,'mp2:input vectors',1,movecs_in)) then
          if (.not.rtdb_cget(rtdb,'dft:output vectors',1,movecs_in))
     $          call util_file_name('movecs',.false.,.false.,movecs_in)
        endif
      else
        if (.not.rtdb_cget(rtdb,'mp2:input vectors',1,movecs_in)) then
          if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_in))
     $          call util_file_name('movecs',.false.,.false.,movecs_in)
        endif
      endif
c
c     Read in MO coefficients and apply swaps as appropriate.
c
      if (.not. movecs_read_header(movecs_in, title_vecs,
     $    basis_vecs, scftype_vecs, nbf_vecs, nsets, nmoarray, 2))
     $    call errquit('direct_mp2: cannot read MO vectors header',0,
     &      INPUT_ERR)
      if (nsets .ne. 1) call errquit('direct_mp2: rhf only',0,
     &      INPUT_ERR)
      nmo = nmoarray(1)
c
      if (nbf.ne.nbf_vecs)  call errquit(
     $    'direct_mp2: number of basis fun does not match input',0,
     &      BASIS_ERR)
c
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'Molecular orbital', 
     $      nbf, 0, g_vecs)) call errquit
     $      ('direct_mp2: failed to create GA for MOs',0,
     &      GA_ERR)
c
      status = MA_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals)
      if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $                    dbl_mb(k_evals), g_vecs)) call
     $      errquit('direct_mp2: cannot read MO vectors',0, DISK_ERR)
      call movecs_swap(rtdb, 'mp2', 'RHF', g_vecs, dbl_mb(k_occ), 
     $      dbl_mb(k_evals))
c
c Get frozen core and virtual offset
c
      if(.not. geom_freeze(rtdb, geom, 'mp2', occ_frozen)) then
         call errquit('direct_mp2: unable to get frozen core',0,
     &      GEOM_ERR)
      endif

      if (.not.rtdb_get(rtdb, 'mp2:frozen virtual', MT_INT, 1,
     $                 vir_frozen)) vir_frozen = 0
c
c Symmetry toggle
c
      osym = (sym_number_ops(geom) .gt. 0) .and. 
     $      sym_abelian_group(geom)
c
c Set index ranges
c
      occ_lo = occ_frozen + 1
      occ_hi = nocc
      nocc_correl = occ_hi - occ_lo + 1
      vir_lo = nocc + 1
      vir_hi = nmo - vir_frozen
      nmo2 = vir_hi - vir_lo + 1
c
c Get Block-length
c
      if (.not.rtdb_get(rtdb,'fourindex:block',MT_INT,
     $                 1,blen)) blen = 32
      blen = max(blen,maxbfsh)
      if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,
     $                 1,oblk)) oblk = .true.
c
c Initialise integrals and screening
c
      if (.not.int_normalize(rtdb,basis))
     &     call errquit('direct_mp2: int_normalize failed',0, INT_ERR)
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
c
c Maximum number of operator matrices we can allocate
c
      nocct = (nocc_correl*(nocc_correl+1))/2
      noper = moints_numoper( basis, algo, nbf, nocc,
     $                      nmo2, blen, omemprint )
      if (noper.gt.nocct) noper = nocct
c
c Allocate operator matrix
c
      if ((algo.eq.'twofold').and.(.not.(ohalf))) then
        g_exch = ga_create_JKblocked(noper,nbf,nmo2,'X oper')
      else
        g_exch = ga_create_JKblocked(noper,nbf,nbf,'X oper')
      endif
c
c Allocate pair energy array
c
*ga:1:0
      if (.not.ga_create(MT_DBL,1,((nocc_correl*(nocc_correl+1))/2),
     $      'epair',0,0,g_epair))
     $      call errquit('direct_mp2:cannot create epair array',0,
     &      GA_ERR)
c
c Dummy loop to count number of passes
c
      npass = 0
      oseg_hi = occ_lo - 1
 322   continue
      npass = npass + 1
      oseg_lo = oseg_hi + 1
      oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
      if (oseg_hi.lt.occ_hi) goto 322
c
c  Enabled disk caching *only* if multipassing required
c
C      odisk = odisk.and.(npass.gt.1)
      if (.not.moints_aodisk_init( odisk, oaoreuse ))
     $    call errquit('direct_mp2: cannot initiliaze ao cache',0,
     &      DISK_ERR)

c     flag and scaling factors for SCS-MP2

      if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $  fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
      if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $  fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
      if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $  scs = .false.       ! flag for SCS-MP2
c
c Print info
c
      if ((ga_nodeid().eq.0).and.(oprint)) then
        call util_print_centered(6,'Direct MP2', 30, .true. )
        write(6,966) nbf, nsh, blen, occ_lo, occ_hi,
     $              vir_lo, vir_hi, movecs_in, noper,npass
 966     format(
     $       10x,'Basis functions:',20x,i5,/,
     $       10x,'Shells:',29x,i5,/,
     $       10x,'Block length:',23x,i5,/,
     $       10x,'Active occupied range:',8x,i4,1x,'-',i5,/,
     $       10x,'Active virtual range: ',8x,i4,1x,'-',i5,/,
     $       10x,'MO coefficients read from:',2x,a25,/,
     $       10x,'Operator matrices in core:',10x,i5,/,
     $       10x,'AO passes required:',17x,i5)
        if (oblk) then
          write(6,881)
 881      format(10x,'Use AO Integral blocking')
        else
          write(6,882)
 882      format(10x,'Use quartet at a time')
        endif
        if (odisk) then
          write(6,971)
 971      format(10x,'AO integral disk caching enabled')
        endif
        if (algo.eq.'sixfold') then
          write(6,975) 
 975      format(10x,'Sixfold algorithm used',/)
        else
          write(6,976) 
 976      format(10x,'Twofold algorithm used',/)
        endif
        if (scs) then
         write(6,'(10x,''SCS-MP2:'',32x,l1)')scs
         write(6,'(10x,''Same spin scaling factor:'',6x,f10.5)')fss
         write(6,'(10x,''Opposite spin scaling factor:'',2x,f10.5/)')fos
        endif
        call util_flush(6)
      endif
c
c Clear 4-index statistics
c
      call moints_stats_clear()
c
c Loop over segments of occupied orbitals
c
      npass = 0
      oseg_hi = occ_lo - 1
 22    continue
      tpass = util_cpusec()
      npass = npass + 1
      oseg_lo = oseg_hi + 1
      oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
c
c J & K build routines: two algorithm choices
c
      if (algo.eq.'sixfold') then
        call moints_build_6x(basis, osym, occ_lo,
     $                     oseg_lo, oseg_hi,
     $                     vir_lo, vir_hi, g_vecs,
     $                     g_coul, .false., g_exch, .true.,
     $                     blen, oblk )
      else
        call moints_build_2x(basis, ohalf, osym, occ_lo,
     $                     oseg_lo, oseg_hi,
     $                     vir_lo, vir_hi, g_vecs,
     $                     g_coul, .false., g_exch, .true.,
     $                     blen, oblk )
      endif


c
c Evaluate and store pair energies
c
      call moints_e2pair( rtdb, nocc, occ_lo, oseg_lo, oseg_hi, 
     $                  vir_lo, vir_hi, vir_lo, vir_hi,
     $                  g_exch, dbl_mb(k_evals), g_epair )
c
c Print progress of this pass
c
      tpass = util_cpusec() - tpass
      if ((ga_nodeid().eq.0).and.
     $     (util_print('progress reports',print_default))) then
        write(6,772) npass, oseg_lo, oseg_hi, tpass
 772     format(5x,'Pass:',i5,5x,'Index range:',i5,2x,'-',i5,
     $         5x,'Time:',f10.2)
        call util_flush(6)
      endif
c
c End multipass loop
c
      if (oseg_hi.lt.occ_hi) goto 22
c$$$           
c$$$      CALL GA_ZERO(G_EXCH)
c$$$      TPASS = UTIL_CPUSEC()
c$$$      CALL MOINTS_BUILD_2X(BASIS, OHALF, OSYM, OCC_LO,
c$$$     $                   OSEG_LO, OSEG_HI,
c$$$     $                   VIR_LO, VIR_HI, G_VECS,
c$$$     $                   G_COUL, .FALSE., G_EXCH, .TRUE.,
c$$$     $                   BLEN, OBLK )
c$$$      CALL MOINTS_E2PAIR( NOCC, OCC_LO, OSEG_LO, OSEG_HI, 
c$$$     $                  VIR_LO, VIR_HI, VIR_LO, VIR_HI,
c$$$     $                  G_EXCH, DBL_MB(K_EVALS), G_EPAIR )
c$$$      TPASS = UTIL_CPUSEC() - TPASS
c$$$      WRITE(6,771) TPASS
c$$$ 771   FORMAT(' Second Pass   Time:',f10.2)


      call moints_aodisk_close( oaokeep )
c
c Print statistics
c
      if (util_print('statistics',print_high))
     $      call moints_stats_print(algo)
c
c Gather pair energies and evaluate mp2 correction
c
      if (ga_nodeid().eq.0) then
        npairs = (nocc_correl*(nocc_correl+1))/2
        status = ma_push_get(MT_DBL, npairs, 'Pair energies',
     $                     l_epair, k_epair)
        ii = 0
        emp2 = 0.0d0
        do i=occ_lo,nocc
          do j=occ_lo,i
            ij = ((i-occ_frozen)*(i-occ_frozen-1))/2 + j - occ_frozen
            call ga_get(g_epair,1,1,ij,ij,e2pair,1)
            dbl_mb(k_epair+ii) = e2pair
            emp2 = emp2 + e2pair
            ii = ii + 1
          enddo
        enddo
c
c Print total and pair energies
c
        if (util_print('pair energies',print_debug)) then
          write(6,970)
 970      format(//,15x,'MP2 Pair Energies',/)
          write(6,970) (dbl_mb(k_epair+ii),ii=0,npairs-1)
 980      format(8f10.5)
        endif
        if (oprint) then
          if (dftmp2) then ! DHDF
            if (scs) then
              write(6,996) scf_energy,emp2,scf_energy+mp2fac*emp2
 996          format(//,10x,43('-'),
     $                /,10x,'DFT energy               ',f24.12,
     $                /,10x,'Unscaled SCS-MP2 energy  ',f24.12,
     $                /,10x,'Total DFT+SCS-MP2 energy ',f24.12,
     $                /,10x,43('-'),/)
            else
              write(6,998) scf_energy,emp2,scf_energy+mp2fac*emp2
 998          format(//,10x,43('-'),
     $                /,10x,'DFT energy             ',f24.12,
     $                /,10x,'Unscaled MP2 energy    ',f24.12,
     $                /,10x,'Total DFT+MP2 energy   ',f24.12,
     $                /,10x,43('-'),/)
            endif
          else ! HF-based MP2
            if (scs) then
              write(6,997) scf_energy, emp2, scf_energy+emp2
 997          format(//,10x,47('-'),
     $                /,10x,'SCF energy             ',f24.12,
     $                /,10x,'SCS Correlation energy ',f24.12,
     $                /,10x,'Total SCS-MP2 energy   ',f24.12,
     $                /,10x,47('-'),/)
            else
              write(6,999) scf_energy, emp2, scf_energy+emp2
 999          format(//,10x,43('-'),
     $                /,10x,'SCF energy         ',f24.12,
     $                /,10x,'Correlation energy ',f24.12,
     $                /,10x,'Total MP2 energy   ',f24.12,
     $                /,10x,43('-'),/)
            endif
          endif
        endif
        status = ma_pop_stack(l_epair)
        call util_flush(6)
      endif
      call ga_sync()
      call ga_brdcst(msg_direct_mp2_emp2, emp2, mdtob(1), 0 )
      ttotal = util_cpusec() - ttotal
      if ((ga_nodeid().eq.0).and.(oprint)) then
        write(6,823) ttotal
 823     format(/,10x,'Total MP2 time:',5x,f10.2)
        call util_flush(6)
      endif
c
      call ecce_print1('correlation energy', mt_dbl, emp2, 1)
c
      emp2 = scf_energy + emp2
      if (.not. rtdb_put(rtdb, 'direct_mp2:energy', mt_dbl, 1, emp2))
     $     call errquit('direct_mp2: failed storing energy in rtdb', 0,
     &      RTDB_ERR)
      if (dftmp2) then
        if (.not. rtdb_put(rtdb,'dft:energy',mt_dbl,1,emp2))
     1       call errquit('direct_mp2: failed storing energy in rtdb',
     2                    1,RTDB_ERR)
      endif
c
      call ecce_print1('total energy', mt_dbl, emp2, 1)
c
c Cleanup
c
      status = ga_destroy(g_epair)
      status = ga_destroy(g_exch)
      status = ma_pop_stack(l_evals)
      call schwarz_tidy()
      status = status.and.ma_pop_stack(l_occ)
      if (.not.status) call errquit('direct_mp2: cannot pop stack?',0,
     &      MA_ERR)
      status = ga_destroy(g_vecs)
      status = geom_destroy(geom)
      call int_terminate()
      if (.not.bas_destroy(basis))
     $   call errquit('direct_mp2: cannot free basis handle',0,
     &      BASIS_ERR)
c
      call util_print_pop()
c
      call ecce_print_module_exit('mp2','ok')
      return
      end


