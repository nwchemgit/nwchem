      logical function moints_full(rtdb)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "util.fh"
#include "inp.fh"
#include "cmointsbuf.fh"
c     
      integer rtdb              ! [input] database handle
c     
c     Rather inefficient use of the partial 4-index transformation
c     to do the complete 4-index transformation.
c
c     Currently collect all integrals and write out a single moints
c     file in various formats, however it would be easy to write
c     out parallel files.
c
c     Ascii format
c
c         title             ! from movecs file
c         norbs,nfzc        ! no. of active MOs, no. of frozen core
c         (map(i)=1,norbs)  ! map from current to old MO order
c         (sym(i)=1,norbs)  ! symmetry of orbitals (0-7 suitable for ieor)
c         (eval(i)=1,norbs) ! eigenvalues of orbitals
c         (occ(i)=1,norbs)  ! occupation numbers of orbitals from movecs
c         g i j k l         ! two-electron integral (one per line)
c         h i j 0 0         ! one-electron integral (one per line)
c         e 0 0 0 0         ! effective nuclear repulsion energy
c         0 -1 -1 -1 -1     ! EOF
c
c     Binary format (variable record length file)
c
c         title             ! c*80
c         norbs,nfzc,buflen ! 
c         efzc              ! effective nuclear repulsion energy
c         (map(i)=1,norbs)  !
c         (sym(i)=1,norbs)  !
c         (eval(i)=1,norbs) !
c         (occ(i)=1,norbs)  !
c         (h(ij),ij=1,norbs*(norbs+1)/2)
c         (g(1:buflen),ijkl(1:buflen),ninrec,pad) ! See moints_write_buf
c         (record with ninrec=0) ! EOF
c
c     COLUMBUS format
c
c         coming soon
c
c     FCIDUMP format
c
c     same as ASCII,
c     but no header
c     and no last line 0. -1 -1 -1 -1
c      
      character*255 movecs_in, moints_file
      integer moints_unit
      parameter (moints_unit=1)
      integer geom, basis
      character*255 title_vecs, basis_vecs, title
      character*20 scftype_vecs
      integer nbf_vecs, nmo(2), nsets
      integer k_occ, l_occ, k_evals, l_evals, g_vecs, k_g, l_g
      integer k_irs, l_irs, k_map, l_map, l_work, k_work
      integer g_mofock, g_aofock, g_aodens, g_tmp, g_hcore
      integer i, nmixed
      integer nbf, nsh, lo, hi
      integer nact, nacttri
      integer occ_frozen, vir_frozen
      integer noper, oseg_lo, oseg_hi, npass, blen, max_noper
      integer g_coul, nsym(8), sym_lo(8), sym_hi(8)
      logical osym, oprint,  ohalf, omemprint, osymblock, ocontam, oblk
      double precision tpass, efzc, eone, etwo
      double precision tol2e, tol_out, tol_sym
      character*7 algo
      character*10 intfmt
c     
      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ga_create_JKblocked
      external ga_create_JKblocked
      integer moints_numoper
      external moints_numoper
      integer moints_occrange
      external moints_occrange
      logical int_normalize
      external int_normalize
c     
      logical cckohn_ci
c
      data g_coul/-1/
      data ohalf/.true./
      data tol2e/1.d-12/
c
      if (.not. rtdb_get(rtdb,'selci:cckohn_ci',mt_log,1,cckohn_ci))
     $     cckohn_ci=.false.
      moints_full = .true.
c     
c     Print settings
c     
      call util_print_push()
      call util_print_rtdb_load(rtdb,'moints')
      oprint = util_print('information',print_default) .and.
     $     ga_nodeid().eq.0
      omemprint = util_print('memory',print_high)
c     
c     Basis and geometry
c     
      if (.not. geom_create(geom,'geometry'))
     $     call errquit('moints: geom failed',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     $     call errquit('moints: geometry not found',0,
     &       RTDB_ERR)
      if(cckohn_ci)then
         if (.not. bas_create(basis, 'all'))
     $        call errquit('moints: basis all failed',0, BASIS_ERR)
         if (.not. bas_rtdb_load(rtdb, geom, basis, 'all'))
     $        call errquit('moints: basis all not found',0,
     &       RTDB_ERR)
      else
         if (.not. bas_create(basis, 'ao basis'))
     $        call errquit('moints: basis failed',0, BASIS_ERR)
         if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $        call errquit('moints: ao basis not found',0, RTDB_ERR)
      endif
c
      if (.not. bas_numbf(basis, nbf)) 
     $     call errquit('moints: failed to get nbf', 0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('moints: failed to get nbf', 0, BASIS_ERR)
c     
c     Load MO vectors
c     
      if (.not.rtdb_cget(rtdb,'fourindex:input vectors',1,movecs_in))
     $     then
         if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_in))
     $        call util_file_name('movecs',.false.,.false.,movecs_in)
      end if
c
      if(cckohn_ci)then
         if (.not.rtdb_cget(rtdb,'cckohn:movecs',1,movecs_in))
     $        then
            call util_file_name
     &           ('kohn_movecs',.false.,.false.,movecs_in)
         end if
      endif
c     
      if (.not. movecs_read_header(movecs_in, title_vecs,
     $     basis_vecs, scftype_vecs, nbf_vecs, nsets, nmo, 2))
     $     call errquit('moints: cannot read MO vectors header',0,
     &       DISK_ERR)
      if (nbf.ne.nbf_vecs) call errquit
     $     ('moints: number of basis fun does not match input',0,
     &       INPUT_ERR)
*     if (nsets .ne. 1) call errquit('moints: cannot do UHF?',0)
c     
      if (.not. ma_push_get(MT_DBL, nbf, 'Evals', l_work, k_work))
     $     call errquit('moints: failed to alloc work ',nbf, MA_ERR)
      if (.not. ma_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals))
     $     call errquit('moints: failed to alloc evals ',nbf, MA_ERR)
      if (.not. ma_push_get(MT_DBL, nbf, 'Occs', l_occ, k_occ))
     $     call errquit('moints: failed to alloc occs ',nbf, MA_ERR)
      if (.not. ma_push_get(MT_INT, nmo(1), 'Irs', l_irs, k_irs))
     $     call errquit('moints: failed to alloc irs ',nmo(1), MA_ERR)
      if (.not. ma_push_get(MT_INT, nmo(1), 'Irs', l_map, k_map))
     $     call errquit('moints: failed to alloc map ',nmo(1), MA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nmo, 'MOs',
     $     0, 0, g_vecs)) call errquit('moints: ga failed', 0, GA_ERR)
      if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $     dbl_mb(k_evals), g_vecs))
     $     call errquit('moints: cannot read MO vectors',0,
     &       DISK_ERR)
c
c     Apply swaps of movecs
c
      call movecs_swap(rtdb,'fourindex', 'ROHF', 
     $     g_vecs,dbl_mb(k_occ),dbl_mb(k_evals))
c     
c     Get other parameters from data base
c     ... frozen core & virtual offsets, block
c     
      if (.not.rtdb_get(rtdb,'fourindex:occ_frozen',MT_INT,1,
     $     occ_frozen)) occ_frozen = 0
      if (.not.rtdb_get(rtdb,'fourindex:vir_frozen',MT_INT,1,
     $     vir_frozen)) vir_frozen = 0
      if (.not.rtdb_get(rtdb,'fourindex:block',MT_INT,1,
     $     blen)) blen = 10
      if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,1,
     $     oblk)) oblk = .false.
      if (.not. rtdb_get(rtdb,'fourindex:osym',MT_LOG,1,
     $     osym)) osym = (sym_number_ops(geom) .gt. 0) .and.
     $     sym_abelian_group(geom)
      if (.not. rtdb_get(rtdb,'fourindex:ignore contamination',MT_LOG,1,
     $     ocontam)) ocontam = .false.
      if (.not. rtdb_get(rtdb,'fourindex:tol', MT_DBL,1,
     $     tol_out)) tol_out = 1.0d-14
      tol_sym = max(tol_out*0.01d0,1d-16)
      if (.not. rtdb_get(rtdb,'fourindex:osymblock',MT_LOG,1,
     $     osymblock)) osymblock = .false.
      if (.not. rtdb_cget(rtdb, 'fourindex:method ', 1, algo))
     $     algo = 'twofold'
      if (.not.rtdb_cget(rtdb,'fourindex:moints file',1,moints_file))
     $     call util_file_name('moints',.true.,.false.,moints_file)
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. rtdb_cget(rtdb,'fourindex:file format',1, intfmt))
     $     intfmt = 'nwchem'
      if (intfmt .eq. 'ascii') then
      else if (intfmt .eq. 'fcidump') then
      else if (intfmt .eq. 'nwchem') then
      else if (intfmt .eq. 'columbus') then
         osymblock = .true.
      else
         call errquit('moints: bad intfmt',0, INPUT_ERR)
      end if
c
c     Set index ranges
c     
      lo   = occ_frozen + 1
      hi   = nbf - vir_frozen
      nact = hi - lo + 1
c
c     Initialise integrals and screening
c     
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('moints_full: int_normalize failed',911,
     &       INT_ERR)
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom,basis)
c     
c     Maximum number of operator matrices we can allocate
c     
      nacttri = (nact*(nact+1))/2
      max_noper = moints_numoper( basis, algo, nbf,
     $     nact, nmo, blen, omemprint )
c     
      noper = min(max_noper, nacttri)
**      noper = 20
c
c     Dummy loop to count number of passes
c
       npass = 0
       oseg_hi = lo - 1
 322   continue
       npass = npass + 1
       oseg_lo = oseg_hi + 1
       oseg_hi = moints_occrange( noper, lo, oseg_lo, hi )
       if (oseg_hi.lt.hi) goto 322
c     
c     Allocate operator matrices
c     
      g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb oper')
      if (.not. ma_push_get(MT_DBL, nact**2, 'g buffer', l_g, k_g))
     $     call errquit('moints: failed to alloc g ',nacttri, MA_ERR)
c     
c     Print info
c     
      if (oprint) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $        30,.true.)
         write(6,*)
         call util_print_centered(6,title,30,.false.)
         write(6,*)
         if (algo.eq.'sixfold') then
            write(6,975) 
 975        format(10x,'Sixfold algorithm used')
         else
            write(6,976) 
 976        format(10x,'Twofold algorithm used')
         end if
         write(6,966) nbf,nsh,nmo(1),nact,lo,hi,blen,noper,npass,
     $        osym, osymblock,intfmt,
     $        moints_file(1:inp_strlen(moints_file))
 966     format(10x,'Number of basis functions:',10x,i5,/,
     $        10x,'Number of shells:',19x,i5,/,
     $        10x,'Number of MOs:   ',19x,i5,/,
     $        10x,'Number of active MOs:',15x,i5,/,
     $        10x,'Active MO range: ',19x,i5,'-',i3,/,
     $        10x,'Block length:',23x,i5,/,
     $        10x,'Number of operator matrices in core:',i5,/,
     $        10x,'Number of passes:',19x,i5,/,
     $        10x,'Use symmetry:    ',23x,l1,/,
     $        10x,'Symmetry block MOs:',21x,l1,/,
     $        10x,'File format:     ',19x,a10,/,
     $        10x,'MO integrals file:',18x,a)
         write(6,9661) movecs_in(1:inp_strlen(movecs_in)),
     $        title_vecs(1:inp_strlen(title_vecs))
 9661    format(10x,'MO coefficients read from:',10x,a,//,10x,a/)
         write(6,*)
         call util_flush(6)
c
         call moints_open(moints_unit, moints_file, intfmt)
      end if
c
c     Symmetrize the MOs and reblock as requested
c
      if (osym .or. osymblock) then
*         call ga_summarize(0)
         call sym_movecs_adapt(basis,tol_sym,g_vecs,int_mb(k_irs),
     $        nmixed)
         if (nmixed .ne. 0) then
            if (.not. ocontam) call errquit
     $           ('moints: MOs are symmetry contaminated',nmixed,
     I           INPUT_ERR)
c     If forced to ignore contamination must orthog and 
c     repeat projection to be sure
            if (oprint) then
               write(6,911) nmixed
 911           format(/' Symmetry contamination of',i3,
     $              ' vectors ignored ... ',/,
     $                 '          repeating projection',
     $              ' and orthogonalization'/)
               call util_flush(6)
            end if
            call ga_orthog_mos(basis, g_vecs)
            call sym_movecs_adapt(basis,tol_sym,g_vecs,int_mb(k_irs),
     $           nmixed)
            if (nmixed .gt. 0) then
               write(6,912) nmixed
 912           format(/' Still have symmetry contamination of',i3,
     $              ' vectors ... cannot proceed')
               call util_flush(6)
               call errquit('moints: symmetry contamination',0,
     &       INPUT_ERR)
            end if
         end if
         if (oprint) then
            write(6,*) '          Symmetry of input active orbitals '
            write(6,*)
            write(6,222) (int_mb(k_irs+i-1),i=lo,hi)
 222        format(40i2)
            write(6,*)
            call util_flush(6)
         end if
         if (intfmt .eq. 'fcidump') then
ccc            
c            call fcidump_symm(iunit,isym,norb,orbsym)
            call fcidump_symm(moints_unit,
     c           hi-lo+1,sym_number_ops(geom)+1,
     c           int_mb(k_irs+lo-1))
         endif
         
         do i = 1, nmo(1)
            int_mb(k_irs+i-1) = int_mb(k_irs+i-1) - 1
         end do
      else
         call ifill(nmo(1),0,int_mb(k_irs),1)
         if (intfmt .eq. 'fcidump') then
            call fcidump_symm(moints_unit,
     c           hi-lo+1,1,
     c           int_mb(k_irs+lo-1))
         endif
      end if
c
c     Keep track of where the MOs will end up
c
      do i = 1,nmo(1)
         int_mb(k_map+i-1) = i
      end do
      if (osymblock) then
         call moints_vecs_sym_sort(g_vecs, nbf, lo, hi,
     $        int_mb(k_irs), int_mb(k_map), dbl_mb(k_work),
     $        nsym, sym_lo, sym_hi)
         if (oprint) then
            write(6,31) nsym
 31         format(/10x, 'Output MOs are blocked by symmetry: ', 8i3)
            call util_flush(6)
         end if
      end if
c
c     Offset so active MOs are counted from one and form the inverse map
c     which is actually what we want.
c
      do i = lo,hi              ! Offset to be in the active range
         int_mb(k_map+i-1) = int_mb(k_map+i-1) - lo + 1
         dbl_mb(k_work-1+int_mb(k_map+i-1)) = i
      end do
      if (oprint) then
         write(6,*)
         write(6,*) '          Map from input to output MO order'
         write(6,*) 
         if (nbf .lt. 100) then
            write(6,223) (i,int_mb(k_map+i-1),i=lo,hi)
 223        format(10(i2,'->',i2,2x))
         else
            write(6,224) (i,int_mb(k_map+i-1),i=lo,hi)
 224        format(9(i3,'->',i3,2x))
         end if
      end if
      do i = 1, nact
         int_mb(k_map+i-1) = nint(dbl_mb(k_work+i-1))
      end do
      if (oprint) then
         write(6,*)
         write(6,*) '          Map from output to input MO order'
         write(6,*) 
         if (nbf .lt. 100) then
            write(6,223) (i,int_mb(k_map+i-1),i=1,nact)
         else
            write(6,224) (i,int_mb(k_map+i-1),i=1,nact)
         end if
         call util_flush(6)
      end if
c
c     Reorder the occupation, eigenvalues to match the active orbitals
c
      call dcopy(nmo(1), dbl_mb(k_occ), 1, dbl_mb(k_work), 1)
      call gather(nact, dbl_mb(k_occ), dbl_mb(k_work), int_mb(k_map))
      call dcopy(nmo(1), dbl_mb(k_evals), 1, dbl_mb(k_work), 1)
      call gather(nact, dbl_mb(k_evals), dbl_mb(k_work), int_mb(k_map))
c
c     Two index transformation of effective core hamiltonian
c
      tpass = -tcgtime()
      if (.not. geom_nuc_rep_energy(geom, efzc)) call errquit
     $     ('moints: geom_nuc_rep_energy failed',0, GEOM_ERR)
      g_hcore = ga_create_atom_blocked(geom, basis, 'T+V')
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore,'kinetic', osym)
      call int_1e_ga(basis, basis, g_hcore,'potential', osym)
      if (util_print('AO T+V',print_debug)) call ga_print(g_hcore)
c
      if (occ_frozen .gt. 0) then
         g_aofock = ga_create_atom_blocked(geom, basis, 'ao fock')
         g_aodens = ga_create_atom_blocked(geom, basis, 'density')
         call ga_zero(g_aofock)
         call ga_dgemm('n','t',nbf, nbf, occ_frozen,
     $        2.0d0, g_vecs, g_vecs, 0.0d0, g_aodens)
         if (util_print('FZC AO density',print_debug)) 
     $        call ga_print(g_aodens)
c
         call scf_get_fock_param(rtdb, tol2e)
         call rhf_fock_2e(geom, basis, g_aodens, g_aofock, tol2e,
     $        .true., .true., osym )
         call scf_tidy(rtdb)
c
         eone = ga_ddot(g_hcore,g_aodens)
         etwo = 0.5d0*ga_ddot(g_aofock,g_aodens)
         efzc = efzc + eone + etwo
         call ga_dadd(1.0d0, g_hcore, 1.0d0, g_aofock, g_hcore)
         if (.not. ga_destroy(g_aodens)) call errquit
     $        ('moints: ga_destroy failed', 0, GA_ERR)
         if (.not. ga_destroy(g_aofock)) call errquit
     $        ('moints: ga_destroy failed', 0, GA_ERR)
      end if
      tpass = tpass + tcgtime()
c
      if (ga_nodeid() .eq. 0) then
         write(6,712) tpass, efzc
 712     format(/10x,'Time for two-index transform:',5x,f8.1,'s'/
     $           10x,'Effective frozen-core energy:',3x,f19.10/)
         call util_flush(6)
      end if
c
      if (osym) call sym_symmetrize(geom, basis, .false., g_hcore)
      if (util_print('AO FZC fock',print_debug)) 
     $     call ga_print(g_hcore)
c
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nmo(1), 'tmp', 0, 0, g_tmp)) 
     $     call errquit('moints: ga create of tmp failed', nmo(1)*nbf,
     &       GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nmo(1), nmo(1), 'MO fock', 
     $     0, 0, g_mofock)) 
     $     call errquit('moints: ga create of tmp failed', nmo(1)**2,
     &       GA_ERR)
c
      if (util_print('vectors',print_debug)) call ga_print(g_vecs)
      call two_index_transf(g_hcore, g_vecs, g_vecs, g_tmp, g_mofock)
      if (util_print('MO FZC fock',print_debug)) call ga_print(g_mofock)
      if (ga_nodeid() .eq. 0) then
         call moints_header_output(moints_unit,
     $        title_vecs, nact, occ_frozen, buflen, nsym,
     $        int_mb(k_map), int_mb(k_irs+occ_frozen),
     $        dbl_mb(k_evals), dbl_mb(k_occ), efzc, intfmt)
         call moints_one_output(moints_unit,
     $     g_mofock, lo, hi, tol_out,osym, int_mb(k_irs), intfmt)
      end if
c
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('moints: ga_destroy failed', 0, GA_ERR)
      if (.not. ga_destroy(g_mofock)) call errquit
     $     ('moints: ga_destroy failed', 0, GA_ERR)
      if (.not. ga_destroy(g_hcore)) call errquit
     $     ('moints: ga_destroy failed', 0, GA_ERR)
c
c     Clear 4-index statistics
c     
      call moints_stats_clear()
c     
c     Loop over segments of occupied orbitals
c     
      npass = 0
      oseg_hi = lo - 1
 22   continue
      tpass = tcgtime()
      npass = npass + 1
      oseg_lo = oseg_hi + 1
      oseg_hi = moints_occrange( noper, lo, oseg_lo, hi )
c     
c     J & K build routines
c     
      if (algo.eq.'sixfold') then
         call moints_build_6x( basis, osym, lo,
     $        oseg_lo, oseg_hi,
     $        lo, hi, g_vecs,
     $        g_coul, .true., -1, .false., blen, oblk )
      else
         call moints_build_2x( basis,  ohalf, osym, lo,
     $        oseg_lo, oseg_hi,
     $        lo, hi, g_vecs,
     $        g_coul, .true., -1, .false., blen, oblk )
      end if
c     
c     Write the MO integrals to disk
c     
      if (ga_nodeid() .eq. 0) call moints_two_output(moints_unit,
     $     g_coul, lo, hi, oseg_lo, oseg_hi, tol_out,
     $     dbl_mb(k_g), osym, int_mb(k_irs), intfmt)
c     
      tpass = tcgtime() - tpass
      if ((ga_nodeid().eq.0).and.
     $     (util_print('progress reports',print_default)))
     $     write(6,772) npass, oseg_lo, oseg_hi, tpass
 772  format(5x,'Pass:',i5,5x,'Index range:',i5,2x,'-',i5,
     $     5x,'Time:',f10.2)
      if (oseg_hi .lt. hi) goto 22
c     
c     End multipass loop
c     
      if (util_print('statistics',print_high))
     $     call moints_stats_print(algo)
c     
      if (ga_nodeid().eq.0) call moints_close(moints_unit, intfmt)
c
      call ga_sync()
c     
c     Tidy up
c     
      call schwarz_tidy()
      call int_terminate()
      if (.not. bas_destroy(basis))
     $     call errquit('moints: basis corrupt',0, BASIS_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit('moints: geom corrupt',0, GEOM_ERR)
      if (.not. ga_destroy(g_coul)) call errquit
     $     ('moints: ga destroy failed ',0, GA_ERR)
      if (.not. ma_pop_stack(l_g))
     $     call errquit('moints: ma failed',0, MA_ERR)
      if (.not. ma_pop_stack(l_map))
     $     call errquit('moints: ma failed',0, MA_ERR)
      if (.not. ma_pop_stack(l_irs))
     $     call errquit('moints: ma failed',0, MA_ERR)
      if (.not. ma_pop_stack(l_occ))
     $     call errquit('moints: ma failed',0, MA_ERR)
      if (.not. ma_pop_stack(l_evals))
     $     call errquit('moints: ma failed',0, MA_ERR)
      if (.not. ma_pop_stack(l_work))
     $     call errquit('moints: ma failed',0, MA_ERR)
      if (.not. ga_destroy(g_vecs))
     $     call errquit('moints: ga_destroy failed',0, GA_ERR)
c     
      call util_print_pop()
c

      end
      subroutine moints_two_output(moints_unit,
     $     g_coul,lo, hi, oseg_lo, oseg_hi,
     $     tol2e, g, osym, irs, intfmt)
      implicit none
#include "errquit.fh"
c
      integer moints_unit
      integer g_coul
      integer lo, hi
      integer oseg_lo, oseg_hi
      double precision tol2e, g(lo:hi,lo:hi)
      logical osym
      integer irs(*)
      character*(*) intfmt
c
c     buffer for nwchem binary format
c
      integer buflen
      parameter (buflen=1024)   ! MUST MATCH MAIN PROGRAM
      double precision buf(buflen)
      integer labels(4,buflen)
      integer n_in_buf
c
      integer i, j, k, l, ltop, ij, lo1, itype
      integer nact2, ijkl, nact
#include "bitops.fh"
c
      lo1 = lo - 1
      nact = hi-lo+1
      nact2 = nact**2
c
      if (intfmt .eq. 'ascii') then
         itype = 1
      else if (intfmt .eq. 'fcidump') then
         itype = 1
      else if (intfmt .eq. 'nwchem') then
         n_in_buf = 0
         itype = 2
      else if (intfmt .eq. 'columbus') then
         itype = 3
      else
         call errquit('moints_two: bad intfmt',0, INPUT_ERR)
      end if
c
      ij = 1
      do i = oseg_lo, oseg_hi
         do j = lo, i
            call ga_get(g_coul, 1, nact2, ij, ij, g, 1)
            do k = lo, i
               ltop = k
               if (k .eq. i) ltop = j
               do l = lo, ltop
                  if (abs(g(l,k)) .gt. tol2e) then
                     if (osym) then
                        ijkl = ieor(ieor(irs(i),irs(j)),
     $                       ieor(irs(k),irs(l)))
                        if (ijkl .ne. 0) then
                           write(6,*) i, j, k, l
                           write(6,*) irs(i), irs(j), irs(k), irs(l)
                           write(6,*) g(l,k)
                           call errquit
     $                          ('moints: integral zero by symmetry',0,
     &       GEOM_ERR)
                        end if
                     end if
                     if (itype .eq. 1) then
                        write(moints_unit,1) g(l,k), 
     $                       i-lo1, j-lo1, k-lo1, l-lo1
 1                      format(e23.16,4i4)
cfcidump -0.4759344611440753E+00                        
                     else if (itype .eq. 2) then
                        n_in_buf = n_in_buf + 1
                        buf(n_in_buf) = g(l,k)
                        labels(1,n_in_buf) = i - lo1
                        labels(2,n_in_buf) = j - lo1
                        labels(3,n_in_buf) = k - lo1
                        labels(4,n_in_buf) = l - lo1
                        if (n_in_buf .eq. buflen) then
                           call moints_wrt_buf(moints_unit,
     $                          n_in_buf, buf, labels)
                           n_in_buf = 0
                        end if
                     else
                        call errquit('moints: bad itype', itype,
     &       INPUT_ERR)
                     end if
                  end if
               end do
            end do
            ij = ij + 1
         end do
      end do
c
      if (itype.eq.2 .and. n_in_buf.gt.0) 
     $     call moints_wrt_buf(moints_unit, n_in_buf, buf, labels)
c
      end
      subroutine moints_one_output(moints_unit,
     $     g_mofock, lo, hi, tol_out, osym, irs, intfmt)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer moints_unit
      integer g_mofock
      integer lo, hi
      double precision tol_out
      logical osym
      integer irs(*)
      character*(*) intfmt
c
      double precision buf(1024)
      integer i, j, lo1
      integer k_buf, l_buf, nact, nacttri
c
c     First check that symmetry is obeyed
c
      if (osym) then
         do i = lo, hi
            call ga_get(g_mofock, 1, hi, i, i, buf, 1)
            do j = lo, i
               if (irs(i).ne.irs(j)) then
                  if (abs(buf(j)) .gt. tol_out) then
                     write(6,*) i, j, irs(i), irs(j), buf(j)
                     call errquit('moints: integral zero by sym',0,
     &       GEOM_ERR)
                  else
                     buf(j) = 0.0d0
                  end if
               end if
            end do
            call ga_put(g_mofock, 1, hi, i, i, buf, 1)
         end do
      end if
c
      if (intfmt .eq. 'ascii'.or.intfmt .eq. 'fcidump') then
         lo1 = lo - 1
         do i = lo, hi
            call ga_get(g_mofock, 1, hi, i, i, buf, 1)
            do j = lo, i
               if (abs(buf(j)) .gt. tol_out) then
                  write(moints_unit,1) buf(j), i-lo1, j-lo1, 0, 0
 1                format(e23.16,4i4)
               end if
            end do
         end do
      else if (intfmt .eq. 'nwchem') then
c
c     Get triangle local and write out
c
         nact = hi - lo + 1
         nacttri = nact*(nact+1)/2
         if (.not. ma_push_get(mt_dbl, nacttri, 'moints_one',
     $        l_buf, k_buf)) call errquit('moints: ma ', nacttri,
     &       MA_ERR)
         lo1 = k_buf
         do i = lo, hi
            call ga_get(g_mofock, lo, i, i, i, dbl_mb(lo1), 1)
            lo1 = lo1 + i - lo + 1
         end do
         write(moints_unit) (dbl_mb(k_buf+i-1),i=1,nacttri)
         if (.not. ma_pop_stack(l_buf))
     $        call errquit('moints: ma corrupt',0, MA_ERR)
      else if (intfmt .eq. 'columbus') then
c
c     Get triangle local and write out
c
         nact = hi - lo + 1
         nacttri = nact*(nact+1)/2
         if (.not. ma_push_get(mt_dbl, nacttri, 'moints_one',
     $        l_buf, k_buf)) call errquit('moints: ma ', nacttri,
     &       MA_ERR)
         lo1 = k_buf
         do i = lo, hi
            call ga_get(g_mofock, lo, i, i, i, dbl_mb(lo1), 1)
            lo1 = lo1 + i - lo + 1
         end do
         write(moints_unit) (dbl_mb(k_buf+i-1),i=1,nacttri)
         if (.not. ma_pop_stack(l_buf))
     $        call errquit('moints: ma corrupt',0, MA_ERR)
      else
         call errquit('moints_one: bad integral format',0, INPUT_ERR)

      end if
c
      end
      subroutine moints_close(moints_unit, intfmt)
      implicit none
#include "errquit.fh"
      integer moints_unit
      character*(*) intfmt
      double precision buf(1)   ! Not used
      integer labels(4)         ! Not used
c
      if (intfmt .eq. 'ascii') then
         write(moints_unit,1) 0.0d0, -1, -1, -1, -1
 1       format(e23.16,4i4)
         close(moints_unit, status='keep', err=1000)
      else if (intfmt .eq. 'nwchem') then
         call moints_wrt_buf(moints_unit, 0, buf, labels)
         close(moints_unit, status='keep', err=1000)
      else if (intfmt .eq. 'columbus'.or.
     c        intfmt .eq. 'fcidump') then
         close(moints_unit, status='keep', err=1000)
      else
         call errquit('moints: unknown integral format',0, INPUT_ERR)
      end if
c
      return
c
 1000 call errquit('moints_close: close failed',0, DISK_ERR)
c
      end
      subroutine moints_open(moints_unit, moints_file, intfmt)
      implicit none
#include "errquit.fh"
      integer moints_unit
      character*(*) moints_file
      character*(*) intfmt
c
      if (intfmt .eq. 'ascii') then
         open(moints_unit,file=moints_file,
     $        form='formatted',status='unknown',
     $        err=1000)
      elseif (intfmt .eq. 'fcidump') then
         open(moints_unit,file=moints_file,
     $        form='formatted',status='old',
     $        position='append',
     $        err=1000)
      else if (intfmt .eq. 'nwchem') then
         open(moints_unit,file=moints_file,
     $        form='unformatted',status='unknown',
     $        err=1000)
      else if (intfmt .eq. 'columbus') then
         open(moints_unit,file=moints_file,
     $        form='unformatted',status='unknown',
     $        err=1000)
      else
         call errquit('moints: unknown file format',0, INPUT_ERR)
      end if
c
      return
c
 1000 write(6,*) ' moints filename=', moints_file
      call errquit('moints_open: failed to open file',0, DISK_ERR)
c
      end
      subroutine moints_vecs_sym_sort(g_vecs, nbf, lo, hi,
     $     irs, map, work, nsym, sym_lo, sym_hi)
      implicit none
#include "errquit.fh"
#include "global.fh"
      integer g_vecs            ! [input] GA handle to vectors
      integer nbf               ! [input]
      integer lo, hi            ! [input] Range of vectors to sort
      integer irs(*)            ! [input/output] Symmetry of orbitals
      integer map(*)            ! [input/output] Mapping indices
      integer work(*)           ! [scratch] Length nbf
      integer nsym(0:7)         ! [output] No. of MOs in each sym in range
      integer sym_lo(0:7)       ! [output] First MO in each sym
      integer sym_hi(0:7)       ! [output] Last MO in each sym
c
c     !! irs(*) assumed to be in the range 0,..,7
c
c     sort the MO vectors in the range (lo:hi) by symmetry.
c
c     map(lo:hi) returns the map from input to output orbitals
c     
c     irs(lo:hi) are permuted to reflect the output order
c
c     nsym(sym)=sym_hi(sym)-sym_lo(sym)+1 return no. and range of orbitals
c     in a given symmetry (sym=0,...,7)
c
c     If there are no orbtials in a symmetry then you will get
c     sym_hi(sym)=sym_lo(sym)-1 which will form a null do-loop
c

c     generate map
c     
      integer i, isym, n, g_tmp
c
      n = lo-1
      do isym = 0, 7
         nsym(isym) = 0
         sym_lo(isym)=n+1
         do i = lo, hi
            if (irs(i).eq.isym) then
               n = n + 1
               nsym(isym) = nsym(isym) + 1
               map(i) = n
               work(map(i)) = irs(i) ! Reorder the actual symmetries
            end if
         end do
         sym_hi(isym)=sym_lo(isym)+nsym(isym)-1
      end do
      if (n .ne. hi) call errquit('moints: bad map?',0, INPUT_ERR)
      do i = lo, hi
         irs(i) = work(i)
      end do
c
c     Sort the vectors
c
      if (.not. ga_duplicate(g_vecs, g_tmp, 'sort buf'))
     $     call errquit('moints: no GA room for sort buf',0, GA_ERR)
      call ga_copy(g_vecs, g_tmp)
c
      do i = lo+ga_nodeid(), hi, ga_nnodes()
         call ga_get(g_tmp, 1, nbf, i, i, work, nbf)
         call ga_put(g_vecs,1, nbf, map(i), map(i), work, nbf)
      end do
c
      if (.not. ga_destroy(g_tmp)) call errquit('moints: GA destroy?',0,
     &       GA_ERR)
c
      end
      subroutine moints_header_output(moints_unit, title,
     $     nact, occ_frozen, buflen, nsym,
     $     map, irs, evals, occ, efzc, intfmt)
      implicit none
#include "errquit.fh"
      integer moints_unit
      character*(*) title
      integer nact
      integer occ_frozen
      integer buflen
      integer nsym(8)
      integer map(1:nact)
      integer irs(1:nact)
      double precision evals(1:nact)
      double precision occ(1:nact)
      double precision efzc
      character*(*) intfmt
c
      character*80 fmt
c
c     Junk for columbus
c
      character*4 slabel(255)
      character*8 bfnlab(255)
      integer numsym
c
      integer ifmt1,l1rec,n1max,ifmt2,l2rec,n2max,info,ierr,i
      common/cmointssifs/info(5)
c
      if (intfmt .eq. 'ascii') then
         write(moints_unit,'(a)') title(1:80)
         write(moints_unit,'(2i6,5x,a)') nact, occ_frozen, '! nmo nfzc'
         fmt = ' '
         write(fmt,1) nact
 1       format('(',i4,'i4,5x,''! map'')')
         write(moints_unit,fmt) map
         fmt = ' '
         write(fmt,2) nact
 2       format('(',i4,'i4,5x,''! irs'')')
         write(moints_unit,fmt) irs
         fmt = ' '
         write(fmt,3) nact
 3       format('(1p,',i4,'d16.8,5x,''! evals'')')
         write(moints_unit,fmt) evals
         fmt = ' '
         write(fmt,4) nact
 4       format('(1p,',i4,'d16.8,5x,''! occ'')')
         write(moints_unit,fmt) occ
         write(moints_unit,5) efzc, 0, 0, 0, 0
 5       format(e23.16,4i4,5x,'! effective nuclear-repulsion energy')
      else if (intfmt.eq.'fcidump') then
         write(moints_unit,6) efzc, 0, 0, 0, 0
 6       format(e23.16,4i4,5x)
c     nothing else to do here. fcidump_header()
      else if (intfmt .eq. 'nwchem') then
         write(moints_unit) title(1:80)
         write(moints_unit) nact, occ_frozen, buflen
         write(moints_unit) efzc
         write(moints_unit) map
         write(moints_unit) irs
         write(moints_unit) evals
         write(moints_unit) occ
      else if (intfmt .eq. 'columbus') then
c
c     Initialize SIFS record-sizes etc.
c
         call sifcfg( 1, -1, nact, 0, ifmt1, l1rec, n1max, ierr)
         if ( ierr .ne. 0 ) call errquit('moitns: sifcfg1 ', ierr,
     &       INPUT_ERR)
         call sifcfg( 2, -1, nact, 0, ifmt2, l2rec, n2max, ierr)
         if ( ierr .ne. 0 ) call errquit('moitns: sifcfg2 ', ierr,
     &       INPUT_ERR)
c
c     Fill up the header and stuff it out there
c
         info(1) = 1
         info(2) = l1rec
         info(3) = n1max
         info(4) = l2rec
         info(5) = n2max
c     
         numsym = 0             ! Figure out no. of symmetry types
         do i = 1, 8
            if (nsym(i).ne.0) numsym = i
         end do
         if (numsym .le. 2) then
            continue
         else if (numsym .le. 4) then
            numsym = 4
         else if (numsym .le. 8) then
            numsym = 8
         else
            call errquit('moints_columbus symmetry nonsense',numsym,
     &       INPUT_ERR)
         end if
c     
         do i = 1,numsym
            write(slabel(i),'(a3,i1)')'sym',i
         end do
         do i = 1,nact
            write(bfnlab(i),'(a4,i4.4)')'MO  ',i
         end do
c     
         call sifwh(moints_unit, 1, numsym, nact,
     &        5, 1, 0, title,
     &        nsym, slabel, info, bfnlab,
     &        -1, efzc, 0, 0, ierr )
         if( ierr.ne.0 ) call errquit('moints: sifwh failed', ierr,
     &       INPUT_ERR)
      else
         call errquit('moints header: unknown intfmt',0,
     &       INPUT_ERR)
      end if
      end
      subroutine moints_wrt_buf(moints_unit, n, buf, labels)
      implicit none
#include "cmointsbuf.fh"
#include "mafdecls.fh"
      integer moints_unit
      integer n,k
      double precision buf(*)
      integer labels(*)
c
c     SWRITE CALL ASSUMES 32 or 64 BIT INTEGER
c
      k = ma_sizeof(mt_int, 1, mt_byte)
      call dcopy(n, buf, 1, fbuf, 1)
      if (n.gt.0) then
         if (k.eq.4) then
            call util_pack_8(4*n, plab, labels)
         else
            call util_pack_16(4*n, plab, labels)
         endif
      endif
      n_in_rec = n
      call swrite(moints_unit, fbuf, buflen+buflen*k/8+1)
c
      end
      subroutine moints_read_buf(moints_unit, n, buf, labels)
      implicit none
#include "cmointsbuf.fh"
#include "mafdecls.fh"
#include "moints.fh"
      integer moints_unit
      integer n,k
      double precision buf(*)
      integer labels(*)
c
c     SWRITE CALL ASSUMES 32 or 64 BIT INTEGER
c
      k = ma_sizeof(mt_int, 1, mt_byte)
      call sread(moints_unit, fbuf, buflen+buflen*k/8+1)
      n = n_in_rec
      if (n.gt.0) then
         if (k.eq.4) then
             call util_unpack_8(4*buflen, plab, labels)
         else
             call util_unpack_16(4*buflen, plab, labels)
         endif
         call dcopy(n, fbuf, 1, buf, 1)
      endif
c
      end
