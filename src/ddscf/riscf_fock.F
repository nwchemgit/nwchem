      subroutine riscf_fock_2e(geom1, ao_basis, ri_basis, nfock,
     $     fcoul, fexch, tol2e1, oskeleton, g_dens, g_fock)
C$Id$
C     calculate the coulomb and exchange contribution to the Fock matrix 
C     from the density matrix and the transformed three-center ERIs

      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "tcgmsg.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "cscf.fh"
#include "cfock.fh"
#include "stdio.fh"
#include "util.fh"
#ifdef RISCF_USE_DRA
#include "dra.fh"
#endif

      integer ao_basis, ri_basis, g_dens, g_fock, g_three, g_indx, 
     $     g_sthree, g_T, g_tr, t, riscf_nr_nu, geom1, g_x1, g_x2, nfock
C     g_three can be GA or DA handle !!!

      logical status, ocoul, oexch, oskeleton, riscf_get_3c_ints
      double precision fcoul, fexch, tol2e1

      dimension g_dens(nfock), g_fock(nfock), fcoul(nfock), fexch(nfock)

      integer Nff, nsh, lsq, lhalf, npair, n_t, nt, nocc, in, il,
     $     l_indx, k_indx, l_three, k_three, l_sthree, k_sthree,
     $     ifirst, ifirst1, ifirst2, ilast, ilast1, ilast2, ifk, mxbf, 
     $     mxints, nloop, nrest, nproc, nt_dim, ic, last_t, g_three_sq

#ifdef RISCF_USE_DRA
      integer req, g_three_io, first_t, num_t, istat, iarray, 
     $     iarray2, next_t
      dimension g_three_io(2), req(2)
#endif

      double precision mone, zero, one, sq2, W, factor

      logical odisk

      character*64 message

      parameter (zero=0.0D0, one=1.0D0, mone=-1.D0)

      sq2 = sqrt(2.D0)
      if (oskeleton) then
        if (ga_nodeid().eq.0) then
          if ( util_print('information',print_low) ) then
            write(luout,*)
            call util_print_centered(luout,
     $           'Symmetry not implemented in RI', 25, .true.)
          endif
        endif
      endif

C     get basis set information
      status = bas_numbf (ao_basis, nbf)
      status = bas_numbf (ri_basis, Nff)
      status = bas_numcont(ao_basis, nsh)
      lsq = nbf * nbf
      lhalf = (nbf * (nbf + 1)) / 2

C     get 3center integrals and indices
      status = riscf_get_3c_ints (ao_basis, ri_basis, g_three, g_indx, 
     $                npair, odisk)
      if (.not. status) then
        call errquit ('RISCF not initialized', 110, UNKNOWN_ERR)
      endif

C     in case of in-core method pass the integral array to riscf_square
      g_three_sq = g_three

C     get number of expansion functions to be processed at the same time
      if ( scftype .eq. 'ROHF' .or. scftype .eq. 'RHF') then
        nocc = nclosed + nopen
      else if ( scftype .eq. 'UHF' ) then
        nocc = MAX ( nalpha, nbeta )
      endif
      n_t = riscf_nr_nu ( nbf, nsh, nocc, npair, ao_basis, dovecs, 
     $     odisk )
      n_t = min ( n_t, Nff )
      if (ga_nodeid().eq.0) then
        if ( util_print('information',print_high) ) then
          message = 'processing XXXXX expansion functions in one step'
          write ( message(12:16),'(I5)') n_t
          write(luout,*)
          call util_print_centered(luout, message, 25, .true.)
        endif
      endif

C     temporary arrays (local and global)
      status = bas_nbf_cn_max (ao_basis, mxbf)
      if (.not.status) call errquit ('problem with bas_nbf_cn_max', 110,
     &       BASIS_ERR)
C     indices (local)
      status = ma_alloc_get (MT_INT, nsh, 'RI indices', l_indx, k_indx)
      mxints = nbf * mxbf
C     compressed integrals
      status = ma_alloc_get (MT_DBL, mxints, 'RI 3c integrals compr.',
     $                       l_three, k_three)
C     block of square matrix
      status = ma_alloc_get (MT_DBL, mxints, 'RI 3c integrals full',
     $                       l_sthree, k_sthree)
C     full 3c square matrix (ij|t) (range of t fixed)
*ga:1:0
      status = ga_create(MT_DBL, n_t*nbf, nbf,'3c-eri square', 0, nbf, 
     $     g_sthree)
      if (.not.status) 
     $     call errquit ('riscf_fock: error creating 3c square', 110,
     &       GA_ERR)
C
      if (dovecs) then
C       Ken Dyalls MO-vector method

        if ( scftype .ne. 'ROHF' .and. scftype .ne. 'RHF' )
     $        call errquit ('Dyall method only for ROHF', 110,
     &       UNKNOWN_ERR)
*ga:1:0
        status = ga_create(MT_DBL, n_t*nbf, nocc, 'Dyall tmp 1',0, 
     $       nocc, g_x1)
        if (.not.status) 
     $       call errquit ('riscf_fock: error creating g_x1', 110,
     &       GA_ERR)
*ga:1:0
        status = ga_create(MT_DBL, nbf, n_t*nocc, 'Dyall tmp 2',0, 
     $       nocc, g_x2)
        if (.not.status) 
     $       call errquit ('riscf_fock: error creating g_x2', 110,
     &       GA_ERR)
      else
C       intermediate array T
*ga:1:0
        status = ga_create(MT_DBL, n_t*nbf, nbf,'3c-eri tmp', 0, nbf,
     $       g_T)
        if (.not.status) 
     $       call errquit ('riscf_fock: error creating g_T', 110,
     &       GA_ERR)
C       scratch array for transposition of T-blocks
*ga:1:0
        status = ga_create(MT_DBL, nbf, nbf,'scr', 0, nbf, g_tr)
        if (.not.status) 
     $       call errquit ('riscf_fock: error creating g_tr', 110,
     &       GA_ERR)
      endif

      nproc = ga_nnodes()

      nloop = Nff / n_t
      nrest = MOD(Nff,n_t)
      if ( nrest.eq.0 ) then
C       last loop is not shorter
        nrest = n_t
      else
C       last (partial) loop not covered by Nff / n_t
        nloop = nloop + 1
      endif

      nt = n_t ! # t values to be processed in one loop

      last_t = min ( nt, Nff )

      if ( odisk ) then
#ifdef RISCF_USE_DRA
C       global arrays to store part of three-center integrals read from da
C       two of them to read asynchronously
*ga:1:0
        status = ga_create(MT_DBL, npair, n_t,'3c-eri #1', 0, nbf, 
     $       g_three_io(1) )
        if (.not.status) 
     $       call errquit ('riscf_fock: error creating io array 1', 110,
     &       GA_ERR)
        status = ga_create(MT_DBL, npair, n_t,'3c-eri #2', 0, nbf, 
     $       g_three_io(2) )
        if (.not.status) 
     $       call errquit ('riscf_fock: error creating io array 2', 110,
     &       GA_ERR)

C       and read the first block
        if (oscfps) call pstat_on(ps_riscf_read)
        istat = dra_read_section (.false., g_three_io(1), 
     $                  1, npair, 1, last_t,
     $       g_three,   1, npair, 1, last_t, req(1) )
        if (oscfps) call pstat_off(ps_riscf_read)

#endif
      endif

      ic = 0
      do il = 1, nloop
        if ( il.eq.nloop ) nt = nrest
        nt_dim = nt * nbf
C       multiple fock matrices
        do ifk = 1, nfock
          ocoul = (fcoul(ifk).ne.0.D0)
          oexch = (fexch(ifk).ne.0.D0)

C         create full square matrices of 3center ERIs for this t range
          if (oscfps) call pstat_on(ps_extract)

          t = (il-1)*n_t

          if ( odisk ) then
#ifdef RISCF_USE_DRA
            next_t = t + n_t
            iarray = mod(ic,2) + 1
C           initialize next read
            if ( il .ne. nloop ) then
              first_t = next_t+1
              last_t = min ( first_t+nt-1, Nff )
              num_t = last_t - first_t + 1
              ic = ic + 1
              iarray2 = mod(ic,2) + 1
              if (oscfps) call pstat_on(ps_riscf_read)
              istat = dra_read_section (.false., 
     $             g_three_io(iarray2), 1, npair, 1, num_t,
     $             g_three,  1,  npair,     first_t, last_t, 
     $             req(iarray2) )
              if (oscfps) call pstat_off(ps_riscf_read)
            endif

C           finish this read
            if (oscfps) call pstat_on(ps_riscf_read)
            istat = dra_wait ( req(iarray) )
            if (oscfps) call pstat_off(ps_riscf_read)

            g_three_sq = g_three_io(iarray)
#endif
          endif

          call riscf_square (ao_basis, g_three_sq, g_sthree,
     $         dbl_mb(k_three), dbl_mb(k_sthree), g_indx, 
     $         int_mb(k_indx), t, nt, nbf, nsh, mxbf, nproc, odisk)

          if (oscfps) call pstat_off(ps_extract)
          if (oexch) then
C           exchange contribution to Fock matrix
            if (oscfps) call pstat_on(ps_exchange)
            if (dovecs) then
              if (scftype.eq.'RHF') then
                factor = 1.D0
              else
                factor = 1.D0
              endif
              
C             X := three x C
              call ga_matmul_patch ('N','N', factor, zero,
     $             g_sthree, 1, nt_dim, 1, nbf,
     $             g_vecs  , 1, nbf,    1, nocc,
     $             g_x1,     1, nt_dim, 1, nocc )
C             bring X in shape for second dgemm
              ifirst1 = 1
              ifirst2 = 1
              do in = 1, nt
                ilast1 = ifirst1 + nbf - 1
                ilast2 = ifirst2 + nocc - 1
                call ga_copy_patch ('N', g_x1, ifirst1, ilast1, 1, nocc,
     $                                   g_x2, 1, nbf, ifirst2, ilast2)
                ifirst1 = ilast1 + 1
                ifirst2 = ilast2 + 1
              enddo

C             K := K + X x X(T)
              call ga_matmul_patch ('N','T', mone, one,
     $             g_x2, 1, nbf, 1, nt*nocc,
     $             g_x2, 1, nt*nocc, 1, nbf,
     $             g_fock(ifk), 1, nbf, 1, nbf)

            else !not dovecs => use density

C             T := three x D
              call ga_matmul_patch ('N','N', one, zero,
     $             g_sthree,   1, nt_dim, 1, nbf,
     $             g_dens(ifk),1,    nbf, 1, nbf,
     $             g_T,        1, nt_dim, 1, nbf )
              
C             transpose blocks of T
              ifirst = 1
              do in = 1, nt
                ilast = ifirst + nbf - 1
                if (oscfps) call pstat_on(ps_copy)
                call ga_copy_patch ('N', g_T, ifirst, ilast, 1, nbf, 
     $               g_tr,     1,   nbf, 1, nbf)
                if (oscfps) call pstat_off(ps_copy)
                if (oscfps) call pstat_on(ps_transpose)
                call ga_copy_patch ('T', g_tr,     1,   nbf, 1, nbf,
     $               g_T, ifirst, ilast, 1, nbf) 
                if (oscfps) call pstat_off(ps_transpose)
                ifirst = ifirst + nbf
              enddo
              
C             K := K + T x three
              call ga_matmul_patch ('T','N', fexch(ifk), one,
     $             g_T,        1,     nbf, 1, nt_dim,
     $             g_sthree,   1, nt_dim, 1, nbf,
     $             g_fock(ifk),1,    nbf, 1, nbf)
              
            endif !dovecs or not

            if (oscfps) call pstat_off(ps_exchange)
          endif !exchange

          if (ocoul) then
C           Coulomb contribution
            if (oscfps) call pstat_on(ps_coulomb)
            do in = 1, nt
C             W = 2 * sum_kl (m|kl) * D_kl
              ifirst = (in-1) * nbf + 1
              ilast = ifirst + nbf - 1
              W = fcoul(ifk) * ga_ddot_patch ( 
     $             g_sthree,   'N', ifirst, ilast, 1, nbf,
     $             g_dens(ifk),'N',      1,   nbf, 1, nbf)
C             J_ij = J_ij + (m|ij) * W
              call ga_dadd_patch ( W, g_sthree,   ifirst, ilast, 1, nbf,
     $             one, g_fock(ifk),     1,   nbf, 1, nbf,
     $                  g_fock(ifk),     1,   nbf, 1, nbf)
            enddo
            if (oscfps) call pstat_off(ps_coulomb)
          endif ! exchange
        enddo   ! nfock
      enddo     ! nloop

C     free the temporary arrays
      status = ma_free_heap ( l_sthree )
      status = ma_free_heap ( l_three )
      status = ma_free_heap ( l_indx )
      if (dovecs) then
        status = ga_destroy ( g_x1)
        status = ga_destroy ( g_x2 )
      else
        status = ga_destroy ( g_tr )
        status = ga_destroy ( g_T )
      endif
      status = ga_destroy ( g_sthree )
#ifdef RISCF_USE_DRA
      if ( odisk ) then
        status = ga_destroy ( g_three_io(2) )
        status = ga_destroy ( g_three_io(1) )
      endif
#endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function riscf_nr_nu ( nbf, nsh, nocc, npair, basis, 
     $     dovecs, odisk )
C     calculate number of expansion functions to be processed in one ga_dgemm
C     from the available memory
      implicit none
#include "errquit.fh"
#include "msgids.fh"
#include "global.fh"
#include "bas.fh"

      integer nbf, nsh, nocc, npair, basis
      logical dovecs, odisk

      integer mfree, ndim, nnodes, mxbf
      logical status

      nnodes = ga_nnodes()

      mfree = ga_memory_avail() / 8

      if ( mfree .lt. 0 ) then
C       memory cannot be determined
        riscf_nr_nu = 10
        return
      endif

      if ( .not. dovecs ) then
C       g_tr
        mfree = mfree - ( nbf * nbf / nnodes )
      endif

C     local arrays from MA
      if ( ga_uses_ma() ) then
C       l_indx
        mfree = mfree - nsh
C       l_three, l_sthree
        status = bas_nbf_cn_max ( basis, mxbf )
        if ( .not. status ) 
     $     call errquit ('riscf_use_disk: error getting basis info',110,
     &       BASIS_ERR)
        mfree = mfree - 2 * nbf * mxbf

C       scratch arrays for ga_dgemm: 3 * block size
        mfree = mfree - 3 * 256 * 256
      endif

C     assume evenly distributed array: free space = min * nnodes
      call ga_igop ( msg_riscf_mem, mfree, 1, 'min' )
      mfree = mfree * nnodes

C     arrays dimensioned with the result of this function
C     g_sthree (and an internal copy in ga_ddot_patch)
      ndim = 2 * nbf * nbf

      if ( dovecs ) then
C       g_x1, g_x2
        ndim = ndim + ( 2 * nbf * nocc )
      else
C       g_T
        ndim = ndim + ( nbf * nbf )
      endif
      if ( odisk ) then
C       g_three_io(1:2)
        ndim = ndim + ( 2 * npair )
      endif

C     don't be too greedy, leave 20%
      mfree = mfree * 8 / 10

      riscf_nr_nu = mfree / ndim

      return
      end
