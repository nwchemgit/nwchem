      logical function riscf_init (rtdb, geom, ao_basis, ri_basis, 
     $        tol2e)
C$Id$

C     calculate three center ERIs, transform them with V**(-1/2)
C     (V-approximation: V(i,j)=(i|j))
C
C     parameters:
C     geom:      geometry handle [input]
C     ao_basis:  basis set handle (AO) [input]
C     ri_basis:  basis set handle (fit) [input]
C     tol2e:     integral threshold [input]
C 
      implicit none
#include "errquit.fh"

#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "stdio.fh"
#include "util.fh"
#ifdef RISCF_USE_DRA
#include "dra.fh"
#endif

      integer rtdb, geom, ao_basis, ri_basis, g_Vsqrt, g_indx,
     $        bases(2), bfsquare
      double precision tol2e, MinEval
      logical status, odisk

      logical riscf_put_3c_ints, Mk_Fit_Xf, riscf_use_disk

      integer k_T_tmp, l_T_tmp, k_V_row, l_V_row, nsh, nsh_F, npair, 
     $     nsqhalf, nbf, Nff, lbuf, l_buf, k_buf, l_scr, k_scr, lscr,
     $     max_bf_ao, max_bf_ri, g_three, nstore
C     g_three is either GA or DA handle for transformed integrals
C     (depending on odisk)

#ifdef RISCF_USE_DRA
      integer istat
      character*64 filename
#endif
C
      logical int_normalize, xc_gotxc
      external int_normalize, xc_gotxc

C     get basis set information
      status = bas_numbf (ao_basis, nbf)
      if (.not.status) call errquit ('RI-SCF: can''t get nbf', 110,
     &       BASIS_ERR)
      status = bas_numbf (ri_basis, Nff)
      if (.not.status) call errquit ('RI-SCF: can''t get Nff', 110,
     &       BASIS_ERR)
      status = bas_nbf_cn_max (ao_basis, max_bf_ao)
      if (.not.status) call errquit ('RI-SCF: can''t get max_bf_ao', 11,
     &       BASIS_ERR)
      bfsquare = max_bf_ao * max_bf_ao
      status = bas_nbf_cn_max (ri_basis, max_bf_ri)

C     initialize for integral calculation
      bases(1) = ao_basis
      bases(2) = ri_basis
      if (.not.int_normalize (rtdb,ri_basis))
     &    call errquit('RI-SCF: int_normalize failed', 12, INT_ERR)
      call int_init (rtdb, 2, bases)
      call schwarz_init(geom, ao_basis)

C  get memory from integral API
      call int_mem_2e3c (lbuf, lscr)

C     calculate V**(-1/2) with Daves routine
      MinEval = 1.D-8
      if (oscfps) call pstat_on(ps_vhalf)
      if ( util_print('information',print_medium)) then
        if ( ga_nodeid().eq.0 ) then
          call util_print_centered (LUOUT,'creating V**-1/2', 20, 
     $         .false.)
        endif
      endif
      status = Mk_Fit_Xf('V', 'Y', ri_basis, MinEval, g_Vsqrt)
      if (oscfps) call pstat_off(ps_vhalf)
      if (.not. status) call errquit ('error creating V**-1/2',110,
     &       UNKNOWN_ERR)

C     temporary arrays (local and global)
C     integral buffer
      status = ma_alloc_get (MT_DBL, lbuf, 'RI buffer', l_buf, k_buf)
      if (.not.status) 
     $     call errquit ('RI-SCF: can''t get buffer space',110, MA_ERR)
C     scratch space
      lscr = max (lscr, bfsquare*Nff) !reuse for transformed 3c-ints
      status = ma_alloc_get (MT_DBL, lscr, 'RI scratch', l_scr, 
     $                       k_scr)
      if (.not.status) 
     $     call errquit ('RI-SCF: can''t get scratch space',110, MA_ERR)
C     (ij|t) for one block of ij
      status = ma_alloc_get (MT_DBL, bfsquare*Nff, 'RI T',
     $                       l_T_tmp, k_T_tmp)
      if (.not.status) 
     $     call errquit ('RI-SCF: can''t get temporary space',110,
     &       MA_ERR)
C     local block of V**(-1/2)
      status = ma_alloc_get (MT_DBL, max_bf_ri*Nff,'RI V_row',
     $                       l_V_row, k_V_row)
      if (.not.status) call errquit ('RI-SCF: can''t get V_row',110,
     &       MA_ERR)

C     get and fill indexing array, determine space needed for 3c-ERIs
      if (oscfps) call pstat_on(ps_index)
      status = bas_check_handle(ao_basis,'RI: ao_basis: invalid handle')
      status = bas_check_handle(ri_basis,'RI: ri_basis: invalid handle')
      status = bas_numcont (ao_basis, nsh)
      status = bas_numcont (ri_basis, nsh_F)
      nsqhalf = (nsh*(nsh+1))/2
*ga:1:0
      status = ga_create(MT_INT, nsh, nsh,'3c index', nsh, 0, g_indx)
*ga:1:0
*     status = ga_create(MT_INT, nsh, nsh,'3c index', 0, nsh, g_indx)
      if (.not.status) call errquit ('RI-SCF: can''t get indx arr.',110,
     &       GA_ERR)
C     does this work for integer arrays?
      call ga_zero (g_indx)
      call riscf_make_index (g_indx, npair,nsh,nsqhalf,tol2e,ao_basis)
      if (oscfps) call pstat_off(ps_index)

C     disk-based method requested?
      if (xc_gotxc()) then
        status = rtdb_get(rtdb, 'dft:ri_store',MT_INT, 1, nstore)
      else
        status = rtdb_get(rtdb, 'scf:ri_store',MT_INT, 1, nstore)
      endif
      if ( .not. status ) call errquit
     $     ('no store option for RI in rtdb', 110, RTDB_ERR)
C     nstore = 1 : automatic choice based on memory size
C       "    = 2 : in-core
C       "    = 3 : disk based
      if ( nstore .eq. 2 ) then
        odisk = .false.
      else if ( nstore .eq. 3 ) then
        odisk = .true.
      else if ( nstore .eq. 1 ) then
        odisk = riscf_use_disk ( ao_basis, npair, nbf, Nff )
      else
        call errquit ('no store option for RI in rtdb', 110, RTDB_ERR)
      endif

 100  continue  ! jump back here after error message
      if ( odisk ) then
#ifdef RISCF_USE_DRA
        if ( ga_nodeid().eq.0 ) then
          if ( util_print('storage',print_medium) ) then
            call util_print_centered(LUOUT,'using disk based algorithm', 
     $           20, .false.)
          endif
        endif

C       create unique file name
        call util_file_prefix ('3c-eri', filename )
C       transformed 3-center integrals on disk
        istat = dra_init ( 10, 1.D10, 1.D10, 1.D7 ) 
C         (max_arrays, max_array_size, total_disk_space, max_memory)
        if ( istat .ne. 0 ) 
     $       call errquit ('RI-SCF: error initializing CHEMIO', istat ,
     &       UNKNOWN_ERR)

        istat = dra_create ( MT_DBL, npair, Nff, '3c-eri', filename, 
     $       DRA_RW, 1, npair, g_three )
        if ( istat .ne. 0 ) 
     $       call errquit ('RI-SCF: can''t create global 3c-ERIs',istat,
     &       UNKNOWN_ERR)
#else
        if ( ga_nodeid().eq.0 ) then
          if ( util_print('storage',print_low) ) then
            write (LUOUT,*)
            call util_print_centered (LUOUT,'DISK ARRAYS NOT AVAILABLE', 
     $           20, .true. )
            call util_print_centered (LUOUT,'attempting in-core method', 
     $         20, .false. )
            write (LUOUT,*)
          endif
        endif
C       switch to in-core and jump back
        odisk = .false.
        goto 100
#endif
      else      ! in-core
        if ( ga_nodeid().eq.0 ) then
          if ( util_print('storage',print_low) ) then
            call util_print_centered (LUOUT,'using in-core algorithm', 
     $           20, .false.)
            write (LUOUT,*)
          endif
        endif
*ga:1:0
        status = ga_create( MT_DBL, npair, Nff, '3c-eri', npair, 0, 
     $                      g_three )
        if (.not.status) 
     $       call errquit ('RI-SCF: can''t allocate 3c-eri array',110,
     &       GA_ERR)
        call ga_zero (g_three)
      endif

      if ( ga_nodeid().eq.0 ) then
        if ( util_print('information',print_medium) ) then
          call util_print_centered (LUOUT,
     $         'calculating 3-center integrals', 20, .false.)
        endif
      endif

      if (oscfps) call pstat_on(ps_transform)
      call riscf_trans_int (g_three, ao_basis, nsh, nbf, nsqhalf,
     $     bfsquare,
     $     npair, ri_basis, nsh_F, Nff, max_bf_ri, g_Vsqrt, tol2e, 
     $     g_indx, 
     $     dbl_mb(k_V_row), dbl_mb(k_T_tmp), dbl_mb(k_buf), lbuf, 
     $     dbl_mb(k_scr), lscr, odisk)
      if (oscfps) call pstat_off(ps_transform)

      if ( ga_nodeid().eq.0 ) then
        if ( util_print('information',print_medium) ) then
          call util_print_centered (LUOUT,'3-center integrals done'
     $         , 20, .false.)
          write (LUOUT,*)
        endif
      endif

      status = riscf_put_3c_ints (ao_basis, ri_basis, g_three, g_indx, 
     $     npair, odisk)
      if (.not.status) 
     $     call errquit ('RI-SCF: can''t put 3c-ERIS in list',110,
     &       UNKNOWN_ERR)

C     free the temporary arrays
      status = ma_free_heap ( l_V_row )
      status = status .and. ma_free_heap ( l_T_tmp )
      status = status .and. ma_free_heap ( l_scr )
      status = status .and. ma_free_heap ( l_buf )
      status = status .and. ga_destroy ( g_Vsqrt )
      if (.not.status)
     $     call errquit ('RI-SCF: problem freeing arrays',110, MEM_ERR)

C     destroy integral setup for expansion basis
      call int_terminate()
      call schwarz_tidy()

C     make ri the default for mixed methods
      call do_riscf (.true.)

      riscf_init = status

      return

      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine riscf_make_index (g_indx, npair, nsh, nsqhalf, tol2e,
     $     ao_basis)

      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "msgids.fh"

      integer g_indx, npair, nsh, nsqhalf, ao_basis
      double precision tol2e

      integer i, j, ifirst, ilast, jfirst, jlast, nfi, nfj, icount, 
     $     isum, idum, me, nproc, k_indx, l_indx, len
      logical status

      me = ga_nodeid()
      nproc = ga_nnodes()

      status = ma_alloc_get (MT_INT, nsh,'RI indices', l_indx, k_indx)
      if (.not. status) call errquit ('couldn''t get local index array',
     $     110, MA_ERR)
      npair = 0
      do i=1, nsh
        if ( mod(i,nproc) .eq. me ) then
          call ifill (i, 0, int_mb(k_indx), 1)
          do j=1, i
            if (schwarz_shell (i,j) .gt. tol2e) then
              status = bas_cn2bfr(ao_basis, i, ifirst, ilast)
              nfi = ilast - ifirst + 1
              if ( i.ne.j ) then
                status = bas_cn2bfr(ao_basis, j, jfirst, jlast)
                nfj = jlast - jfirst + 1
                int_mb(k_indx+j-1) = nfi * nfj
              else
                int_mb(k_indx+j-1) = nfi*(nfi+1)/2
              endif
              npair = npair + int_mb(k_indx+j-1)
            endif
          enddo
          call ga_put (g_indx, 1, i, i, i, int_mb(k_indx), i)
C         only one of them is probably necessary
C         check distribution (ga_create above)
        endif
      enddo
      call ga_sync()

C     sum up
C     is there a way to do it in parallel?
      icount = 0
      isum = 1
      if (me .eq. 0) then
        do i=1, nsh
          call ga_get (g_indx, 1, i, i, i, int_mb(k_indx), i)
          do j=0, i-1
            if ( int_mb(k_indx+j) .gt. 0 ) then
              icount = icount+1
              idum = int_mb(k_indx+j)
              int_mb(k_indx+j) = isum !starting point in compressed (i,j)
              isum = isum + idum
            endif
          enddo
          call ga_put (g_indx, 1, i, i, i, int_mb(k_indx), i)
          if ( i .gt. 1 ) then
            call ga_put (g_indx, i, i, 1, i-1, int_mb(k_indx), 1)
          endif
        enddo
        npair = isum - 1
      endif
      len = MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst (Msg_3c_index, npair, len, 0)

      status = ma_free_heap (l_indx)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function riscf_use_disk ( basis, npair, nbf, Nff )
C     find out if there is enough memory to hold the 3-center integrals in
C     core (=>false) or if disk based algorithm should be used (=>true)

      implicit none
#include "errquit.fh"

#include "msgids.fh"
#include "global.fh"
#include "bas.fh"

      integer basis, npair, nbf, Nff

      integer mfree, mneeded, nnodes, mxbf
      logical status

      nnodes = ga_nnodes()

      mfree = ga_memory_avail() / 8

      if ( mfree .le. 0 ) then
C       free memory cannot be determined
C       we are conservative and choose the disk based algorithm
        riscf_use_disk = .false.
        return
      endif

C     global arrays in riscf, assuming n_t = 1
      mfree = mfree - 2 * nbf * nbf / nnodes

C     local arrays from MA
      if ( ga_uses_ma() ) then
C       l_three, l_sthree
        status = bas_nbf_cn_max ( basis, mxbf )
        if ( .not. status ) 
     $     call errquit ('riscf_use_disk: error getting basis info',110,
     &       BASIS_ERR)
        mfree = mfree - nbf * mxbf

C       scratch arrays for ga_dgemm: 3 * block size
        mfree = mfree - 3 * 256 * 256
      endif

C     assume evenly distributed array: free space = min * nnodes
      call ga_igop ( msg_riscf_free, mfree, 1, 'min' )
      mfree = mfree * nnodes

C     size of the 3-center integral array
      mneeded = npair * Nff

C     don't use all the memory 
C     (there are several O(nbf) arrrays unaccounted for)

      if ( dble(mneeded) .gt. ( 0.8D0 * mfree ) ) then
        riscf_use_disk = .true.
      else
        riscf_use_disk = .false.
      endif

      end
