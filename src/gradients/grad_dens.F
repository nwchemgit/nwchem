      subroutine grad_dens (geom, basis, g_dens, g_wdens, g_eigen_diag, 
     $     evals, occ, ndens, nbf, nbf2, nopen, nclosed, scftype, 
     $     movecs, omp2)

c$Id: grad_dens.F,v 1.6 1996-10-21 00:39:35 vg038 Exp $

C     calculate the density and energy-weighted density matrices from the
C     movecs (and Lagrangian in ROHF case) on disk

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"

      integer geom, basis, g_dens, g_wdens, g_eigen_diag, ndens, 
     $     nbf, nbf2, nopen, nclosed
      double precision evals, occ
      character*4 scftype
      character*(*) movecs
      logical omp2

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read, file_read_ga
      external movecs_read, file_read_ga

      character*255 fname ! filename for lagrangian (ROHF) 
                          ! or second-order contributions (MP2)

      dimension evals(nbf2), occ(nbf2), g_dens(ndens)

C     local variables
      integer nocc, l_ind, k_ind, l_ind2, k_ind2, g_vecs, g_vecs2, i
      double precision zero, one, mone, two, mtwo
      logical status

      parameter ( mtwo = -2.D0 )
      parameter ( mone = -1.D0 )
      parameter ( zero =  0.D0 )
      parameter ( one  =  1.D0 )
      parameter ( two  =  2.D0 )

C  get MO vectors from file (taken from Roberts scf_movecs_read)

      g_vecs = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      if ( scftype.eq. 'UHF' ) then
        g_vecs2 = ga_create_atom_blocked (geom, basis, 'uhf mo vectors')
      endif

      status = movecs_read (movecs, 1, occ, evals, g_vecs)
      if (.not.status) then
        call errquit('grad_dens: could not read mo vectors', 110)
      endif

      if (scftype .eq. 'UHF') then
        if (.not. movecs_read(movecs, 2, occ(nbf+1), evals(nbf+1),
     $       g_vecs2)) then
          if ( util_print('information',print_low)) then
            if (ga_nodeid() .eq. 0) then
              write(luout,*) '  Duplicating RHF/ROHF vectors for UHF'
            endif
          endif
          call ga_copy(g_vecs, g_vecs2)
          call dcopy (nbf, occ, 1, occ(nbf+1), 1)
          call dcopy (nbf, evals, 1, evals(nbf+1), 1)
        endif
      endif
      
      if (scftype .eq. 'RHF') then

C       energy-weighted density

C       create eigenvalue diagonal matrix
        nocc = nclosed
        if (ga_nodeid() .eq. 0) then
          status = MA_alloc_get(MT_INT, nocc, 'indices', l_ind, k_ind)
          status = MA_alloc_get(MT_INT, nocc, 'indices', l_ind2, k_ind2)
          do i=1, nocc
            int_mb( k_ind+i-1 ) = i
            int_mb( k_ind2+i-1 ) = i
          enddo
          call ga_scatter ( g_eigen_diag, evals, int_mb(k_ind),
     $           int_mb(k_ind2), nocc )
          status = MA_free_heap ( l_ind )
          status = MA_free_heap ( l_ind2 )
        endif

C       multiply to weighted density
        call ga_matmul_patch ('n','n', two, zero,
     $                         g_vecs,       1, nbf,  1, nocc,
     $                         g_eigen_diag, 1, nocc, 1, nocc,
     $                         g_dens(1),    1, nbf,  1, nocc )
        call ga_matmul_patch ('n','t', one, zero,
     $                         g_dens(1),    1, nbf,  1, nocc,
     $                         g_vecs,       1, nocc, 1, nbf,
     $                         g_wdens,      1, nbf,  1, nbf )

        if ( omp2 ) then
C         read 2nd order correction to weighted density from disk
C         use g_dens(3) as buffer
          call util_file_name ('w_a_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(3) )
          call ga_add ( one, g_wdens, mtwo, g_dens(3), g_wdens )
        endif

C       density matrix
        call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $       g_vecs,    1, nbf, 1, nclosed,
     $       g_vecs,    1, nclosed, 1, nbf,
     $       g_dens(1), 1, nbf, 1, nbf)
        call ga_symmetrize(g_dens(1))

        if ( omp2 ) then
C         g_dens(1) = D(HF+(2))     for 1-el. contribution
C         g_dens(2) = D(HF)         for 2-el. contr.
C         g_dens(3) = D(HF+2(2))        "

          call ga_copy ( g_dens(1), g_dens(2) )
C         read and add MP2 contribution
          call util_file_name ('p_a_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(3) )
          call ga_scale ( g_dens(3), two )
          call ga_dadd ( one, g_dens(3), one, g_dens(1), g_dens(1) )
          call ga_dadd ( one, g_dens(3), one, g_dens(1), g_dens(3) )
        endif ! mp2

      else if (scftype .eq. 'ROHF') then

        nocc = nopen + nclosed
C       read in Lagrangian matrix ('eigenvalue matrix' in ROHF case)
        call util_file_name('lagr', .true.,.false., fname )
        status = file_read_ga( fname, g_eigen_diag )
        if (.not.status) then
          call errquit('grad_dens: could not read Lagrangian', 110)
        endif

C       multiply to weighted density
        call ga_matmul_patch ('n','n', two, zero,
     $                         g_vecs,       1, nbf,  1, nocc,
     $                         g_eigen_diag, 1, nocc, 1, nocc,
     $                         g_dens(1),    1, nbf,  1, nocc )
        call ga_matmul_patch ('n','t', one, zero,
     $                         g_dens(1),    1, nbf,  1, nocc,
     $                         g_vecs,       1, nocc, 1, nbf,
     $                         g_wdens,      1, nbf,  1, nbf )

C       density matrices:
C       g_dens(1): D(closed) + 1/2 D(open)
C       g_dens(2): D(closed)
C       g_dens(3): D(open)
C       density contains a factor 2 compared to most literature
C       (same as in ddscf (I think...))

        call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $       g_vecs,    1, nbf,     1, nclosed,
     $       g_vecs,    1, nclosed, 1, nbf,
     $       g_dens(1), 1, nbf,     1, nbf )

        call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $       g_vecs,    1,         nbf,  nclosed+1, nocc,
     $       g_vecs,    nclosed+1, nocc, 1,         nbf,
     $       g_dens(3), 1,         nbf,  1,         nbf )

        call ga_copy ( g_dens(1), g_dens(2) )
        call ga_dadd (one, g_dens(1), 0.5D0, g_dens(3), g_dens(1))
        call ga_symmetrize(g_dens(1))
        call ga_symmetrize(g_dens(2))
        call ga_symmetrize(g_dens(3))

      else  ! UHF

C       energy-weighted density

C       create eigenvalue diagonal matrix
        nocc = nclosed + nopen  ! == nalpha
        if (ga_nodeid() .eq. 0) then
          status = MA_alloc_get(MT_INT, nocc, 'indices', l_ind, k_ind)
          status = MA_alloc_get(MT_INT, nocc, 'indices', l_ind2, k_ind2)
          do i=1, nocc
            int_mb( k_ind+i-1 ) = i
            int_mb( k_ind2+i-1 ) = i
          enddo
          call ga_scatter ( g_eigen_diag, evals, int_mb(k_ind),
     $           int_mb(k_ind2), nocc )
        endif

C       multiply to weighted density
        call ga_matmul_patch ('n','n', one, zero,
     $                         g_vecs,       1, nbf,  1, nocc,
     $                         g_eigen_diag, 1, nocc, 1, nocc,
     $                         g_dens(1),    1, nbf,  1, nocc )
        call ga_matmul_patch ('n','t', one, zero,
     $                         g_dens(1),    1, nbf,  1, nocc,
     $                         g_vecs,       1, nocc, 1, nbf,
     $                         g_wdens,      1, nbf,  1, nbf )

        nocc = nclosed   ! nbeta
        if (ga_nodeid() .eq. 0) then
          do i=1, nocc
            int_mb( k_ind+i-1 ) = i
            int_mb( k_ind2+i-1 ) = i
          enddo
          call ga_scatter ( g_eigen_diag, evals(nbf+1), int_mb(k_ind),
     $           int_mb(k_ind2), nocc )
          status = MA_free_heap ( l_ind )
          status = MA_free_heap ( l_ind2 )
        endif

C       multiply to weighted density and add to spin up part
        call ga_matmul_patch ('n','n', one, zero,
     $                         g_vecs2,      1, nbf,  1, nocc,
     $                         g_eigen_diag, 1, nocc, 1, nocc,
     $                         g_dens(1),    1, nbf,  1, nocc )
        call ga_matmul_patch ('n','t', one, one,
     $                         g_dens(1),    1, nbf,  1, nocc,
     $                         g_vecs2,      1, nocc, 1, nbf,
     $                         g_wdens,      1, nbf,  1, nbf )
        call ga_symmetrize(g_wdens)

        if ( omp2 ) then
C         read 2nd order correction to weighted density from disk
C         use g_dens(3) as buffer
          call util_file_name ('w_a_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(3) )
          call ga_add ( one, g_wdens, mone, g_dens(3), g_wdens )
          call util_file_name ('w_b_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(3) )
          call ga_add ( one, g_wdens, mone, g_dens(3), g_wdens )
        endif

C       density matrix
C       D_up and D_down (or and b)
        call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $       g_vecs,    1, nbf,           1, nopen+nclosed,
     $       g_vecs,    1, nopen+nclosed, 1, nbf,
     $       g_dens(1), 1, nbf,           1, nbf)
        call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $       g_vecs2,      1, nbf,     1, nclosed,
     $       g_vecs2,      1, nclosed, 1, nbf,
     $       g_vecs,       1, nbf,     1, nbf)

        if ( .not. omp2 ) then
C         D+ and D- (Formalism from Dupuis & King, JCP 68(9), p4000
          call ga_dadd ( one, g_dens(1), mone, g_vecs, g_dens(2) )
          call ga_dadd ( one, g_dens(1), one,  g_vecs, g_dens(1) )
          call ga_symmetrize(g_dens(1))
          call ga_symmetrize(g_dens(2))
        else   ! UMP2
C         g_dens(1) = D(a+b) (HF+(2))   for 1-el. contribution
C         g_dens(2) = D(a) (HF)         for 2-el. contr.
C         g_dens(3) = D(b) (HF)           "
C         g_dens(4) = D(a) (HF+2(2))      "
C         g_dens(5) = D(b) (HF+2(2))      "

C         move HF matrices to g_dens(2:3)
          call ga_copy ( g_dens(1), g_dens(2) )
          call ga_copy ( g_vecs, g_dens(3) )

          call ga_dadd ( one, g_dens(1), one,  g_dens(3), g_dens(1) )
C         read and add MP2 contribution for spin up
          call util_file_name ('p_a_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(4) )
          call ga_dadd ( one, g_dens(4), one, g_dens(1), g_dens(1) )
          call ga_dadd ( two, g_dens(4), one, g_dens(2), g_dens(4) )
C         spin down
          call util_file_name ('p_b_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(5) )
          call ga_dadd ( one, g_dens(5), one, g_dens(1), g_dens(1) )
          call ga_dadd ( two, g_dens(5), one, g_dens(3), g_dens(5) )
          call ga_symmetrize(g_dens(1))
          call ga_symmetrize(g_dens(2))
          call ga_symmetrize(g_dens(3))
          call ga_symmetrize(g_dens(4))
          call ga_symmetrize(g_dens(5))
        endif   !  UMP2
c       
      endif

      status = ga_destroy ( g_vecs )
      if ( scftype .eq. 'UHF' ) then
        status = ga_destroy ( g_vecs2 )
      endif

      return
      end
