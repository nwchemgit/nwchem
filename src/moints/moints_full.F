      subroutine moints_full(rtdb)
C     $Id: moints_full.F,v 1.3 1996-02-01 01:33:55 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "context.fh"
#include "sym.fh"
#include "util.fh"
#include "inp.fh"
c     
      integer rtdb              ! [input] database handle
c     
c     Rather inefficient use of the partial 4-index transformation
c     to do the complete 4-index transformation.
c     
      character*255 movecs_in, moints_file
      integer moints_unit
      parameter (moints_unit=1)
      integer geom, basis
      character*255 title_vecs, basis_vecs
      integer nbf_vecs, nmo(2), nsets
      integer k_occ, l_occ, k_evals, l_evals, g_vecs, k_g, l_g
      integer k_irs, l_irs, k_map, l_map
      integer g_mofock, g_aofock, g_aodens, g_tmp, g_hcore
      integer i, nmixed
      integer nbf, nsh, lo, hi
      integer nact, nacttri
      integer occ_frozen, vir_frozen
      integer noper, oseg_lo, oseg_hi, npass, blen, max_noper
      integer g_coul, nsym(8), sym_lo(8), sym_hi(8)
      logical osym, oprint,  ohalf, omemprint, osymblock
      double precision tpass, efzc, eone, etwo
      double precision tol2e, tol_out
      character*7 algo
c     
      logical movecs_read_header, movecs_read, movecs_write
      external movecs_read_header, movecs_read, movecs_write
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ga_create_JKblocked
      external ga_create_JKblocked
      integer moints_numoper
      external moints_numoper
      integer moints_occrange
      external moints_occrange
c     
      data g_coul/-1/
      data ohalf/.true./
      data tol2e/1.d-12/
      data tol_out/1.d-9/
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
     $     call errquit('moints: geom failed',0)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     $     call errquit('moints: geometry not found',0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('moints: basis failed',0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('moints: ao basis not found',0)
c     
      if (.not. bas_numbf(basis, nbf)) 
     $     call errquit('moints: failed to get nbf', 0)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('moints: failed to get nbf', 0)
c     
c     Load MO vectors
c     
      if (.not.rtdb_cget(rtdb,'fourindex:input vectors',1,movecs_in))
     $     then
         if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_in))
     $        call util_file_prefix('movecs',movecs_in)
      endif
      if (.not. movecs_read_header(movecs_in, title_vecs,
     $     basis_vecs, nbf_vecs, nsets, nmo, 2))
     $     call errquit('moints: cannot read MO vectors header',0)
      if (nbf.ne.nbf_vecs) call errquit
     $     ('moints: number of basis fun does not match input',0)
      if (nsets .ne. 1) call errquit('moints: cannot do UHF?',0)
c     
      if (.not. ma_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals))
     $     call errquit('moints: failed to alloc evals ',nbf)
      if (.not. ma_push_get(MT_DBL, nbf, 'Occs', l_occ, k_occ))
     $     call errquit('moints: failed to alloc occs ',nbf)
      if (.not. ma_push_get(MT_INT, nmo, 'Irs', l_irs, k_irs))
     $     call errquit('moints: failed to alloc irs ',nmo)
      if (.not. ma_push_get(MT_INT, nmo, 'Irs', l_map, k_map))
     $     call errquit('moints: failed to alloc map ',nmo)
      if (.not. ga_create(MT_DBL, nbf, nmo, 'MOs',
     $     1, 1, g_vecs)) call errquit('rohf_nr_solve: ga_cre grad', 0)
      if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $     dbl_mb(k_evals), g_vecs))
     $     call errquit('moints: cannot read MO vectors',0)
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
      if (.not. rtdb_get(rtdb,'fourindex:osym',MT_LOG,1,
     $     osym)) osym = (sym_number_ops(geom) .gt. 0) 
      if (.not. rtdb_get(rtdb,'fourindex:osymblock',MT_LOG,1,
     $     osymblock)) osymblock = .false.
      if (.not. rtdb_cget(rtdb, 'fourindex:method ', 1, algo))
     $     algo = 'twofold'
      if (.not.rtdb_cget(rtdb,'fourindex:moints file',1,moints_file))
     $     call util_file_prefix('moints',moints_file)
c
c     Set index ranges
c     
      lo   = occ_frozen + 1
      hi   = nbf - vir_frozen
      nact = hi - lo + 1
c
c     Initialise integrals and screening
c     
      call int_init(rtdb, 1, basis)
      call int_normalize(basis)
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
     $     call errquit('moints: failed to alloc g ',nacttri)
c     
c     Print info
c     
      if (oprint) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $        30,.true.)
         write(6,*)
         if (algo.eq.'sixfold') then
            write(6,975) 
 975        format(10x,'Sixfold algorithm used')
         else
            write(6,976) 
 976        format(10x,'Twofold algorithm used')
         endif
         write(6,966) nbf,nsh,nmo(1),nact,lo,hi,blen,noper,npass,
     $        osym, osymblock,
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
     $        10x,'MO integrals file:',18x,a)
         write(6,9661) movecs_in(1:inp_strlen(movecs_in)),
     $        title_vecs(1:inp_strlen(title_vecs))
 9661    format(10x,'MO coefficients read from:',10x,a,//,10x,a/)
         write(6,*)
         call util_flush(6)
c
         call moints_open(moints_unit, moints_file)
      endif
c
c     Symmetrize the MOs and reblock as requested
c
      if (osym .or. osymblock) then
         call sym_movecs_adapt(basis,1.d-9,g_vecs,int_mb(k_irs),nmixed)
         if (nmixed .ne. 0) call errquit
     $        ('moints: MOs are symmetry contaminated',0)
         if (oprint) then
            write(6,*) '          Symmetry of input active orbitals '
            write(6,*)
            write(6,222) (int_mb(k_irs+i-1),i=lo,hi)
 222        format(40i2)
            write(6,*)
            call util_flush(6)
         endif
         do i = 1, nmo(1)
            int_mb(k_irs+i-1) = int_mb(k_irs+i-1) - 1
         enddo
      endif
c
c     Keep track of where the MOs will end up
c
      do i = lo, hi
         int_mb(k_map+i-1) = i
      enddo
      if (osymblock) then
         call moints_vecs_sym_sort(g_vecs, nbf, lo, hi,
     $        int_mb(k_irs), int_mb(k_map), dbl_mb(k_occ),
     $        nsym, sym_lo, sym_hi)
         if (oprint) then
            write(6,31) nsym
 31         format(/10x, 'Output MOs are blocked by symmetry: ', 8i3)
            call util_flush(6)
         endif
      else
      endif
      if (oprint) then
         write(6,*)
         write(6,*) '          Map from input to output MO order'
         write(6,*) 
         if (nbf .lt. 100) then
            write(6,223) (i,int_mb(k_map+i-1)-lo+1,i=lo,hi)
 223        format(10(i2,'->',i2,2x))
         else
            write(6,224) (i,int_mb(k_map+i-1)-lo+1,i=lo,hi)
 224        format(9(i3,'->',i3,2x))
         endif
         call util_flush(6)
      endif
c
c     Two index transformation of effective core hamiltonian
c
      tpass = -tcgtime()
      if (.not. geom_nuc_rep_energy(geom, efzc)) call errquit
     $     ('moints: geom_nuc_rep_energy failed',0)
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
     $        ('moints: ga_destroy failed', 0)
         if (.not. ga_destroy(g_aofock)) call errquit
     $        ('moints: ga_destroy failed', 0)
      endif
      tpass = tpass + tcgtime()
c
      if (ga_nodeid() .eq. 0) then
         write(6,712) tpass, efzc
 712     format(/10x,'Time for two-index transform:',5x,f8.1,'s'/
     $           10x,'Effective frozen-core energy:',3x,f19.10/)
         call util_flush(6)
      endif
c
      if (osym) call sym_symmetrize(geom, basis, .false., g_hcore)
      if (util_print('AO FZC fock',print_debug)) 
     $     call ga_print(g_hcore)
c
      if (.not. ga_create(MT_DBL, nbf, nmo(1), 'tmp', 1, 1, g_tmp)) 
     $     call errquit('moints: ga create of tmp failed', nmo(1)*nbf)
      if (.not. ga_create(MT_DBL, nmo(1), nmo(1), 'MO fock', 
     $     1, 1, g_mofock)) 
     $     call errquit('moints: ga create of tmp failed', nmo(1)**2)
c
      if (util_print('vectors',print_debug)) call ga_print(g_vecs)
      call two_index_transf(g_hcore, g_vecs, g_vecs, g_tmp, g_mofock)
      if (util_print('MO FZC fock',print_debug)) call ga_print(g_mofock)
      if (ga_nodeid() .eq. 0) then
         call moints_efzc_output(moints_unit, efzc)
         call moints_one_output(moints_unit,
     $     g_mofock, lo, hi, tol_out,osym, int_mb(k_irs))
      endif
c
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('moints: ga_destroy failed', 0)
      if (.not. ga_destroy(g_mofock)) call errquit
     $     ('moints: ga_destroy failed', 0)
      if (.not. ga_destroy(g_hcore)) call errquit
     $     ('moints: ga_destroy failed', 0)
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
     $        g_coul, .true., -1, .false., blen )
      else
         call moints_build_2x( basis,  ohalf, osym, lo,
     $        oseg_lo, oseg_hi,
     $        lo, hi, g_vecs,
     $        g_coul, .true., -1, .false., blen )
      endif
c     
c     Write the MO integrals to disk ... right now just the
c     old ascii full-ci format
c     
      if (ga_nodeid() .eq. 0) call moints_two_output(moints_unit,
     $     g_coul, lo, hi, oseg_lo, oseg_hi, tol_out,
     $     dbl_mb(k_g), osym, int_mb(k_irs))
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
      if (ga_nodeid().eq.0) call moints_close(moints_unit)
c
      call ga_sync()
c     
c     Tidy up
c     
      call schwarz_tidy()
      call int_terminate()
      if (.not. bas_destroy(basis))
     $     call errquit('moints: basis corrupt',0)
      if (.not. geom_destroy(geom))
     $     call errquit('moints: geom corrupt',0)
      if (.not. ga_destroy(g_coul)) call errquit
     $     ('moints: ga destroy failed ',0)
      if (.not. ma_pop_stack(l_g))
     $     call errquit('moints: ma failed',0)
      if (.not. ma_pop_stack(l_map))
     $     call errquit('moints: ma failed',0)
      if (.not. ma_pop_stack(l_irs))
     $     call errquit('moints: ma failed',0)
      if (.not. ma_pop_stack(l_occ))
     $     call errquit('moints: ma failed',0)
      if (.not. ma_pop_stack(l_evals))
     $     call errquit('moints: ma failed',0)
      if (.not. ga_destroy(g_vecs))
     $     call errquit('moints: ga_destroy failed',0)
c     
      call util_print_pop()
c

      end
      subroutine moints_two_output(moints_unit,
     $     g_coul,lo, hi, oseg_lo, oseg_hi,
     $     tol2e, g, osym, irs)
      implicit none
c
      integer moints_unit
      integer g_coul
      integer lo, hi
      integer oseg_lo, oseg_hi
      double precision tol2e, g(lo:hi,lo:hi)
      logical osym
      integer irs(*)
c
      integer i, j, k, l, ltop, ij, lo1
      integer nact2, ijkl, nact
#include "bitops.fh"
c
      lo1 = lo - 1
      nact = hi-lo+1
      nact2 = nact**2
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
     $                          ('moints: integral zero by symmetry',0)
                        endif
                     endif
                     write(moints_unit,1) g(l,k), 
     $                    i-lo1, j-lo1, k-lo1, l-lo1
 1                   format(f18.12,4i4)
                  endif
               enddo
            enddo
            ij = ij + 1
         enddo
      enddo
c
      end
      subroutine moints_one_output(moints_unit,
     $     g_mofock, lo, hi, tol_out, osym, irs)
      implicit none
c
      integer moints_unit
      integer g_mofock
      integer lo, hi
      double precision tol_out
      logical osym
      integer irs(*)
c
      double precision buf(1024)
      integer i, j, lo1
c
      lo1 = lo - 1
      do i = lo, hi
         call ga_get(g_mofock, 1, hi, i, i, buf, 1)
         do j = lo, i
            if (abs(buf(j)) .gt. tol_out) then
               if (osym) then
                  if (irs(i).ne.irs(j)) then
                     write(6,*) i, j, irs(i), irs(j), buf(j)
                     call errquit('moints: integral zero by symmetry',0)
                  endif
               endif
               write(moints_unit,1) buf(j), i-lo1, j-lo1, 0, 0
 1             format(f18.12,4i4)
            endif
         enddo
      enddo
c
      end
      subroutine moints_efzc_output(moints_unit, efzc)
      implicit none
      integer moints_unit
      double precision efzc
c
      write(moints_unit,1) efzc, 0, 0, 0, 0
 1    format(f18.12,4i4)
c
      end
      subroutine moints_close(moints_unit)
      implicit none
      integer moints_unit
c
      write(moints_unit,1) 0.0d0, -1, -1, -1, -1
 1    format(f18.12,4i4)
c
      close(moints_unit, status='keep', err=1000)
      return
c
 1000 call errquit('moints_close: close failed',0)
c
      end
      subroutine moints_open(moints_unit, moints_file)
      implicit none
      integer moints_unit
      character*(*) moints_file
c
      open(moints_unit,file=moints_file,
     $     form='formatted',status='unknown',
     $     err=1000)
c
      return
c
 1000 write(6,*) ' moints filename=', moints_file
      call errquit('moints_open: failed to open file',0)
c
      end
      subroutine moints_vecs_sym_sort(g_vecs, nbf, lo, hi,
     $     irs, map, work, nsym, sym_lo, sym_hi)
      implicit none
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
c     map(lo:hi) returns the map from input to output orbtials
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
            endif
         enddo
         sym_hi(isym)=sym_lo(isym)+nsym(isym)-1
      enddo
      if (n .ne. hi) call errquit('moints: bad map?',0)
      do i = lo, hi
         irs(i) = work(i)
      enddo
c
c     Sort the vectors
c
      if (.not. ga_duplicate(g_vecs, g_tmp, 'sort buf'))
     $     call errquit('moints: no GA room for sort buf',0)
      call ga_copy(g_vecs, g_tmp)
c
      do i = lo+ga_nodeid(), hi, ga_nnodes()
         call ga_get(g_tmp, 1, nbf, i, i, work, 1)
         call ga_put(g_vecs,1, nbf, map(i), map(i), work, 1)
      enddo
c
      if (.not. ga_destroy(g_tmp)) call errquit('moints: GA destroy?',0)
c
      write(6,*) ' nsym ', nsym
      write(6,*) '   lo ', sym_lo
      write(6,*) '   hi ', sym_hi
c
      call ga_print(g_vecs)

c
      end
