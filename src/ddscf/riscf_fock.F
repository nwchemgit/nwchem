      subroutine riscf_fock_2e(geom, ao_basis, ri_basis, nfock,
     $     fcoul, fexch, tol2e, oskeleton, g_dens, g_fock)
C$Id: riscf_fock.F,v 1.2 1995-12-16 06:19:46 d3g681 Exp $
C     calculate the coulomb and exchange contribution to the Fock matrix 
C     from the density matrix and the transformed three-center ERIs

      implicit double precision (a-h,o-z)

#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "pstat.fh"
#include "cscfps.fh"

      integer rtdb
      integer geom, ao_basis, ri_basis, g_dens, g_fock, g_three, g_indx, 
     $     g_sthree, g_T, g_tr, t, riscf_nr_nu
      logical status, ocoul, oexch, oskeleton, riscf_get_3c_ints, 
     $     riscf_init

      dimension g_dens(nfock), g_fock(nfock), fcoul(nfock), fexch(nfock)

      parameter (zero=0.0D0, one=1.0D0)

      if (oskeleton) then
        if (ga_nodeid().eq.0) then
          write(6,*)
          call util_print_centered(6,
     $         'Symmetry not implemented in RI', 25, .true.)
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
     $                npair)
      if (.not. status) then
        status = riscf_init (rtdb, geom, ao_basis, ri_basis, tol2e)
        status = riscf_get_3c_ints (ao_basis, ri_basis, g_three, g_indx,
     $                npair)
      endif

C     get number of expansion functions to be processed at the same time
      n_t = riscf_nr_nu()

C     temporary arrays (local and global)
      status = bas_nbf_cn_max (ao_basis, mxbf)
      if (.not.status) call errquit ('problem with bas_nbf_cn_max', 110)
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
      status = ga_create (MT_DBL, n_t*nbf, nbf,'3c-eri square', 1, nbf, 
     $     g_sthree)
C     intermediate array T
      status = ga_create (MT_DBL, n_t*nbf, nbf,'3c-eri tmp', 1, nbf,
     $     g_T)
C     scratch array for transposition of T-blocks
      status = ga_create (MT_DBL, nbf, nbf,'scr', 1, nbf, g_tr)

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

      do il = 1, nloop
        if ( il.eq.nloop ) nt = nrest
        nt_dim = nt * nbf
C       multiple fock matrices
        do ifk = 1, nfock
          ocoul = (fcoul(ifk).ne.0.D0)
          oexch = (fexch(ifk).ne.0.D0)

C         create full square matrix of 3center ERIs for this t range
          if (oscfps) call pstat_on(ps_extract)
c         
          t = (il-1)*n_t
          do in = 1, nt
            t = t+1
            call riscf_square (ao_basis, g_three, g_sthree,
     $           dbl_mb(k_three), dbl_mb(k_sthree), g_indx, 
     $           int_mb(k_indx), t, in, nbf, nsh, mxbf, nproc)
          enddo

          if (oscfps) call pstat_off(ps_extract)
          if (oexch) then
C           exchange contribution to Fock matrix
            if (oscfps) call pstat_on(ps_exchange)
C           T := three x D
            call ga_matmul_patch ('N','N', one, zero,
     $           g_sthree,   1, nt_dim, 1, nbf,
     $           g_dens(ifk),1,    nbf, 1, nbf,
     $           g_T,        1, nt_dim, 1, nbf )

C           transpose blocks of T
            ifirst = 1
            do in = 1, nt
              ilast = ifirst + nbf - 1
              if (oscfps) call pstat_on(ps_copy)
              call ga_copy_patch ('N', g_T, ifirst, ilast, 1, nbf, 
     $                                 g_tr,     1,   nbf, 1, nbf)
              if (oscfps) call pstat_off(ps_copy)
              if (oscfps) call pstat_on(ps_transpose)
              call ga_copy_patch ('T', g_tr,     1,   nbf, 1, nbf,
     $                                 g_T, ifirst, ilast, 1, nbf) 
              if (oscfps) call pstat_off(ps_transpose)
              ifirst = ifirst + nbf
            enddo

C           K := K + T x three
            call ga_matmul_patch ('T','N', fexch(ifk), one,
     $           g_T,        1,     nbf, 1, nt_dim,
     $           g_sthree,   1, nt_dim, 1, nbf,
     $           g_fock(ifk),1,    nbf, 1, nbf)

*           Call ga_dgemm ('T','N', nbf, nbf, nt_dim, fexch(ifk), 
*    $           g_sthree, g_T, one, g_fock(ifk))
            if (oscfps) call pstat_off(ps_exchange)
          endif
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
          endif
        enddo
      enddo

C     free the temporary arrays
      status = ma_free_heap ( l_sthree )
      status = ma_free_heap ( l_three )
      status = ma_free_heap ( l_indx )
      status = ga_destroy ( g_tr )
      status = ga_destroy ( g_T )
      status = ga_destroy ( g_sthree )

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function riscf_nr_nu()
C     calculate number of expansion functions to be processed in one ga_dgemm
C     from the available memory
C     at the moment a dummy (hoping for better memory control routines)
#include "mafdecls.fh"
#include "global.fh"

      logical lma, llim

      llim = ga_memory_limited()
      lma  = ga_uses_ma()

      riscf_nr_nu = 10

      return
      end
