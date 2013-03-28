      subroutine grad_dens (geom, basis, g_dens, g_wdens, g_eigen_diag, 
     $     evals, occ, ndens, nbf, nbf2, nopen, nclosed, scftype, 
     $     movecs, omp2, odft, pdm2, coeff, lfrac_occ, osep)

c     $Id$

C     calculate the density and energy-weighted density matrices from the
C     movecs (and Lagrangian in ROHF case) on disk

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "eaf.fh"

      integer ndens,geom, basis, g_dens(ndens), g_wdens, g_eigen_diag,
     $     nbf, nbf2, nopen, nclosed
      double precision evals(*), occ(*), pdm2(*), coeff(*)
      character*(*) scftype
      character*(*) movecs
      double precision offset   ! EAF file offset
      character*255 pdmfilename ! EAF file name
      integer pdmfile           ! EAF file handle
      integer l_buf, k_buf

      logical omp2
      logical odft              ! [input] doing DFT (ROKS needs to know)
      logical osep              ! [input] get seperate alpha and beta density
                                ! for UHF

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read, file_read_ga
      external movecs_read, file_read_ga

      character*255 fname       ! filename for lagrangian (ROHF) 
                                ! or second-order contributions (MP2)

C     local variables
      integer nocc, l_ind, k_ind, l_ind2, k_ind2, g_vecs, g_vecs2, i,
     .     g_tmp,g_tmp2, itmpm, ltmpm, nocca,noccb,mxnoc
      double precision zero, one, mone, two, mtwo
      logical status,lfrac_occ
      double precision toll
      double precision docca,doccb

      parameter ( mtwo = -2.D0 )
      parameter ( mone = -1.D0 )
      parameter ( zero =  0.D0 )
      parameter ( one  =  1.D0 )
      parameter ( two  =  2.D0 )
      parameter ( toll =  1.d-9)
c
C     get MO vectors from file (taken from Roberts scf_movecs_read)

      docca = 0.0d0
      doccb = 0.0d0
      g_vecs = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      if (.not. movecs_read (movecs, 1, occ, evals, g_vecs)) 
     $     call errquit('grad_dens: could not read mo vectors', 110,
     &       DISK_ERR)
      do i = 1, nbf
        docca = docca + occ(i)
      enddo
      if (scftype .eq. 'UHF') then
         g_vecs2 = ga_create_atom_blocked(geom, basis, 'uhf mo vectors')
         if (.not. movecs_read(movecs, 2, occ(nbf+1), evals(nbf+1),
     $        g_vecs2)) call errquit
     $        ('grad_dens: could not read beta UHF MOs',0, DISK_ERR)
         do i = 1, nbf
           doccb = doccb + occ(nbf+i)
         enddo
      end if
c
      if (scftype .eq. 'RHF') then
c
C        energy-weighted density
c
C        create eigenvalue diagonal matrix
         nocc = nclosed
         if(lfrac_occ) then
c
c           fractional occupation, therefore check new nocs
c
            if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &         call errquit('dftgforce: failed to alloc tmpm',0, MA_ERR)
c
            g_tmp = ga_create_atom_blocked(geom, basis,
     &              'frac vecs')
            call ga_zero(g_tmp)
            do i = ga_nodeid()+1, nbf, ga_nnodes()
              call get_col(g_vecs, nbf, i, DBL_MB(itmpm))
              call dscal(nbf, occ(i)*0.5d0, DBL_MB(itmpm), 1)
              call put_col(g_tmp, nbf, i, DBL_MB(itmpm))
           enddo
           do i=1,nbf 
              if(occ(i).ge.toll) nocc=i
           enddo
           if (.not.ma_pop_stack(ltmpm))
     &        call errquit('dftg_force: cannot pop stack',0, MA_ERR)
           if(nocc.gt.nclosed) then
             if(.not.ga_destroy ( g_eigen_diag )) call errquit(
     (            'grad_dens: cannot destroy eigendiag',11, GA_ERR)
             if (.not. ga_create(mt_dbl, nocc, nocc,'eigen_diag', 
     ,            1, nocc, g_eigen_diag)) 
     ,         call errquit('gradients: ga diags?',nocc*nocc, GA_ERR)
             call ga_zero(g_eigen_diag)
           endif
           nclosed=nocc
         else
           g_tmp=g_vecs
         endif
#if 0
         if (ga_nodeid() .eq. 0) then
            status = MA_alloc_get(MT_INT, nocc, 'indices', l_ind, k_ind)
            status = MA_alloc_get(MT_INT, nocc, 'indices',l_ind2,k_ind2)
            do i=1, nocc
               int_mb( k_ind+i-1 ) = i
               int_mb( k_ind2+i-1 ) = i
            end do
            call ga_scatter ( g_eigen_diag, evals, int_mb(k_ind),
     $           int_mb(k_ind2), nocc )
            status = MA_free_heap ( l_ind )
            status = MA_free_heap ( l_ind2 )
         end if
#else
         call grad_diagb(g_eigen_diag, evals, nocc)
#endif
c
C        multiply to weighted density
c
         call ga_matmul_patch ('n','n', two, zero,
     $        g_tmp,       1, nbf,  1, nocc,
     $        g_eigen_diag, 1, nocc, 1, nocc,
     $        g_dens(1),    1, nbf,  1, nocc )
         call ga_matmul_patch ('n','t', one, zero,
     $        g_dens(1),    1, nbf,  1, nocc,
     $        g_vecs,       1, nocc, 1, nbf,
     $        g_wdens,      1, nbf,  1, nbf )

         if ( omp2 ) then
c
C           read 2nd order correction to weighted density from disk
C           use g_dens(3) as buffer
c
            call util_file_name ('w_a_ao', .true., .true., fname )
            status = file_read_ga ( fname, g_dens(3) )
            call ga_add ( one, g_wdens, mtwo, g_dens(3), g_wdens )
         end if
c
C        density matrix
c
         call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $        g_tmp,    1, nbf, 1, nocc,
     $        g_vecs,    1, nocc, 1, nbf,
     $        g_dens(1), 1, nbf, 1, nbf)
         call ga_symmetrize(g_dens(1))
c
         if ( omp2 ) then
c
C           g_dens(1) = D(HF+(2))     for 1-el. contribution
C           g_dens(2) = D(HF)         for 2-el. contr.
C           g_dens(3) = D(HF+2(2))        "

            call ga_copy ( g_dens(1), g_dens(2) )
c
C           read and add MP2 contribution
c
            call util_file_name ('p_a_ao', .true., .true., fname )
            status = file_read_ga ( fname, g_dens(3) )
            call ga_scale ( g_dens(3), two )
            call ga_dadd ( one, g_dens(3), one, g_dens(1), g_dens(1) )
            call ga_dadd ( one, g_dens(3), one, g_dens(1), g_dens(3) )
         end if                  ! mp2

      else if (scftype .eq. 'ROHF') then

         nocc = nopen + nclosed
C        read in Lagrangian matrix ('eigenvalue matrix' in ROHF case)
         call util_file_name('lagr', .true.,.false., fname )
         status = file_read_ga( fname, g_eigen_diag )
         if (.not.status) then
            call errquit('grad_dens: could not read Lagrangian', 110,
     &       GA_ERR)
         end if
c
C        multiply to weighted density
c
         call ga_matmul_patch ('n','n', two, zero,
     $        g_vecs,       1, nbf,  1, nocc,
     $        g_eigen_diag, 1, nocc, 1, nocc,
     $        g_dens(1),    1, nbf,  1, nocc )
         call ga_matmul_patch ('n','t', one, zero,
     $        g_dens(1),    1, nbf,  1, nocc,
     $        g_vecs,       1, nocc, 1, nbf,
     $        g_wdens,      1, nbf,  1, nbf )
c
         if (.not.odft) then
c
C          density matrices:
C          g_dens(1): D(closed) + 1/2 D(open)
C          g_dens(2): D(closed)
C          g_dens(3): D(open)
C          density contains a factor 2 compared to most literature
C          (same as in ddscf (I think...))
c
           call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $          g_vecs,    1, nbf,     1, nclosed,
     $          g_vecs,    1, nclosed, 1, nbf,
     $          g_dens(1), 1, nbf,     1, nbf )

           call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $          g_vecs,    1,         nbf,  nclosed+1, nocc,
     $          g_vecs,    nclosed+1, nocc, 1,         nbf,
     $          g_dens(3), 1,         nbf,  1,         nbf )

           call ga_copy ( g_dens(1), g_dens(2) )
           call ga_dadd (one, g_dens(1), 0.5D0, g_dens(3), g_dens(1))
           call ga_symmetrize(g_dens(1))
           call ga_symmetrize(g_dens(2))
           call ga_symmetrize(g_dens(3))
c
         else
c
C          density matrices:
C          g_dens(1): D(alpha)
C          g_dens(2): D(beta)
c          if (.not.osep) then:
C          g_dens(1): D(alpha)+D(beta)
C          g_dens(2): D(alpha)-D(beta)
c
           call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $          g_vecs,    1, nbf,     1, nocc,
     $          g_vecs,    1, nocc,    1, nbf,
     $          g_dens(1), 1, nbf,     1, nbf )

           call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $          g_vecs,    1,   nbf,     1,   nclosed,
     $          g_vecs,    1,   nclosed, 1,   nbf,
     $          g_dens(2), 1,   nbf,     1,   nbf )

           if (.not. osep) then
             call ga_dadd ( one, g_dens(1),  one, g_dens(2), g_dens(1) )
             call ga_dadd ( one, g_dens(1), mtwo, g_dens(2), g_dens(2) )
           endif

           call ga_symmetrize(g_dens(1))
           call ga_symmetrize(g_dens(2))
c
         endif
c
      else if (scftype .eq. 'UHF') then
c
C        energy-weighted density
c
C        create eigenvalue diagonal matrix
         if (docca.ge.doccb) then
           nocca = nclosed + nopen ! == nalpha
           noccb = nclosed         ! nbeta
         else
           nocca = nclosed         ! == nalpha
           noccb = nclosed + nopen ! nbeta
         endif
         if(lfrac_occ ) then
c
c          fractional occupation, therefore check new nocs
c
           if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &        call errquit('dftgforce: failed to alloc tmpm',0, MA_ERR)
c
           g_tmp = ga_create_atom_blocked(geom, basis,
     &             'frac vecs')
           call ga_zero(g_tmp)
           do i = ga_nodeid()+1, nbf, ga_nnodes()
             call get_col(g_vecs, nbf, i, DBL_MB(itmpm))
             call dscal(nbf, occ(i), DBL_MB(itmpm), 1)
             call put_col(g_tmp, nbf, i, DBL_MB(itmpm))
           enddo
           do i=1,nbf 
             if(occ(i).ge.toll) nocca=i
           enddo
           do i=nbf+1,nbf*2 
             if(occ(i).ge.toll) noccb=i-nbf
           enddo
           if (.not.ma_pop_stack(ltmpm))
     &       call errquit('dftg_force: cannot pop stack',0, MA_ERR)
         else
           g_tmp=g_vecs
         endif
c
         if(nocca.gt.nopen+nclosed) then
           if(.not.ga_destroy ( g_eigen_diag )) call errquit(
     &        'grad_dens: cannot destroy eigendiag',11, GA_ERR)
           if (.not. ga_create(mt_dbl, nocca, nocca,'eigen_diag', 
     ,         1, nocca, g_eigen_diag)) 
     ,       call errquit('gradients: ga diags?',nocca*nocca, GA_ERR)
           call ga_zero(g_eigen_diag)
        endif
#if 0
        if (ga_nodeid() .eq. 0) then
           mxnoc=max(nocca,noccb)
           status=MA_alloc_get(MT_INT,mxnoc, 'indices', l_ind, k_ind)
           status=MA_alloc_get(MT_INT,mxnoc, 'indices',l_ind2,k_ind2)
           do i=1, nocca
              int_mb( k_ind+i-1 ) = i
              int_mb( k_ind2+i-1 ) = i
           end do
           call ga_scatter ( g_eigen_diag, evals, int_mb(k_ind),
     $          int_mb(k_ind2), nocca)
        end if
#else
        call grad_diagb(g_eigen_diag, evals, nocca)
#endif
c
C       multiply to weighted density
        call ga_matmul_patch ('n','n', one, zero,
     $       g_tmp,       1, nbf,  1, nocca,
     $       g_eigen_diag, 1, nocca, 1, nocca,
     $       g_dens(1),    1, nbf,  1, nocca )
        call ga_matmul_patch ('n','t', one, zero,
     $       g_dens(1),    1, nbf,  1, nocca,
     $       g_vecs,       1, nocca, 1, nbf,
     $       g_wdens,      1, nbf,  1, nbf )

        if(lfrac_occ) then
c
c         fractional occupation, therefore check new nocs
c
          if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &       call errquit('dftgforce: failed to alloc tmpm',0, MA_ERR)
c
          g_tmp2 = ga_create_atom_blocked(geom, basis,
     &             'frac vecs')
          call ga_zero(g_tmp2)
          do i = ga_nodeid()+1, nbf, ga_nnodes()
              call get_col(g_vecs2, nbf, i, DBL_MB(itmpm))
              call dscal(nbf, occ(i+nbf), DBL_MB(itmpm), 1)
              call put_col(g_tmp2, nbf, i, DBL_MB(itmpm))
          enddo
          if (.not.ma_pop_stack(ltmpm))
     &        call errquit('dftg_force: cannot pop stack',0, MA_ERR)
          if(noccb.gt.nclosed) then
             if(.not.ga_destroy ( g_eigen_diag )) call errquit(
     &            'grad_dens: cannot destroy eigendiag',11, GA_ERR)
             if (.not. ga_create(mt_dbl, noccb, noccb,'eigen_diag', 
     &            1, noccb, g_eigen_diag)) 
     &            call errquit('gradients: ga diags?',noccb*noccb,
     &                         GA_ERR)
             call ga_zero(g_eigen_diag)
          endif
        else
           g_tmp2=g_vecs2
        endif
#if 0
        if (ga_nodeid() .eq. 0) then
           do i=1, noccb
              int_mb( k_ind+i-1 ) = i
              int_mb( k_ind2+i-1 ) = i
           end do
           call ga_scatter ( g_eigen_diag, evals(nbf+1), int_mb(k_ind),
     $          int_mb(k_ind2), noccb )
           status = MA_free_heap ( l_ind )
           status = MA_free_heap ( l_ind2 )
        end if
#else
        call grad_diagb(g_eigen_diag, evals(nbf+1), noccb)
#endif

C       multiply to weighted density and add to spin up part
        call ga_matmul_patch ('n','n', one, zero,
     $       g_tmp2,      1, nbf,  1, noccb,
     $       g_eigen_diag, 1, noccb, 1, noccb,
     $       g_dens(1),    1, nbf,  1, noccb )
        call ga_matmul_patch ('n','t', one, one,
     $       g_dens(1),    1, nbf,  1, noccb,
     $       g_vecs2,      1, noccb, 1, nbf,
     $       g_wdens,      1, nbf,  1, nbf )
        call ga_symmetrize(g_wdens)

        if ( omp2 ) then
c
C         read 2nd order correction to weighted density from disk
C         use g_dens(3) as buffer
c
          call util_file_name ('w_a_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(3) )
          call ga_add ( one, g_wdens, mone, g_dens(3), g_wdens )
          call util_file_name ('w_b_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(3) )
          call ga_add ( one, g_wdens, mone, g_dens(3), g_wdens )
        end if
c
C       density matrix
C       D_up and D_down (or and b)
        call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $       g_tmp,    1, nbf,           1, nocca,
     $       g_vecs,    1, nocca, 1, nbf,
     $       g_dens(1), 1, nbf,           1, nbf)
        call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $       g_tmp2,      1, nbf,     1, noccb,
     $       g_vecs2,      1, noccb, 1, nbf,
     $       g_vecs,       1, nbf,     1, nbf)
c
        if ( .not. omp2 ) then
c
C         D+ and D- (Formalism from Dupuis & King, JCP 68(9), p4000
c
          if (.not. osep) then
            call ga_dadd ( one, g_dens(1), mone, g_vecs, g_dens(2) )
            call ga_dadd ( one, g_dens(1), one,  g_vecs, g_dens(1) )
          else
            call ga_copy (g_vecs,g_dens(2))
          endif
          call ga_symmetrize(g_dens(1))
          call ga_symmetrize(g_dens(2))
        else                   ! UMP2
c
C         g_dens(1) = D(a+b) (HF+(2))   for 1-el. contribution
C         g_dens(2) = D(a) (HF)         for 2-el. contr.
C         g_dens(3) = D(b) (HF)           "
C         g_dens(4) = D(a) (HF+2(2))      "
C         g_dens(5) = D(b) (HF+2(2))      "
c
C         move HF matrices to g_dens(2:3)
          call ga_copy ( g_dens(1), g_dens(2) )
          call ga_copy ( g_vecs, g_dens(3) )

          call ga_dadd ( one, g_dens(1), one,  g_dens(3), g_dens(1) )
c
C         read and add MP2 contribution for spin up
c
          call util_file_name ('p_a_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(4) )
          call ga_dadd ( one, g_dens(4), one, g_dens(1), g_dens(1) )
          call ga_dadd ( two, g_dens(4), one, g_dens(2), g_dens(4) )
c
C         spin down
c
          call util_file_name ('p_b_ao', .true., .true., fname )
          status = file_read_ga ( fname, g_dens(5) )
          call ga_dadd ( one, g_dens(5), one, g_dens(1), g_dens(1) )
          call ga_dadd ( two, g_dens(5), one, g_dens(3), g_dens(5) )
          call ga_symmetrize(g_dens(1))
          call ga_symmetrize(g_dens(2))
          call ga_symmetrize(g_dens(3))
          call ga_symmetrize(g_dens(4))
          call ga_symmetrize(g_dens(5))
        end if                  !  UMP2
c
      else if (scftype .eq. 'MCSCF') then
c
         if (ga_nodeid().eq.0) then
            call util_file_name('mcpdm',.true.,.false.,pdmfilename)
            if (eaf_open(pdmfilename, eaf_r, pdmfile) .ne. 0)
     $           call errquit('mcscf grad: failed opening 2pdm file',0,
     &                        GA_ERR)
            offset = 0.0d0

            if (.not. ma_push_get(mt_dbl,nbf*nbf,'buffer',l_buf,k_buf))
     $           call errquit('mcscf grad: memory problem',nbf**2,
     &                        MA_ERR)

            if (eaf_read(pdmfile, offset, dbl_mb(k_buf), 8*nbf**2).ne.0)
     $        call errquit('mcscf grad: failed reading lagrangian',0,
     &                     DISK_ERR)
            call ga_put(g_wdens, 1, nbf, 1, nbf, dbl_mb(k_buf), nbf)
            offset = offset + 8*nbf**2

            if (eaf_read(pdmfile, offset, dbl_mb(k_buf), 8*nbf**2).ne.0)
     $        call errquit('mcscf grad: failed reading 1-density',0,
     &                     DISK_ERR)
            call ga_put(g_dens(1), 1, nbf, 1, nbf, dbl_mb(k_buf), nbf)
            offset = offset + 8*nbf**2

            if (eaf_read(pdmfile, offset, pdm2, 8*nopen**4).ne.0)
     $           call errquit('mcscf grad: failed reading 2pdm ',0,
     &                        DISK_ERR)

            call ga_get(g_vecs, 1, nbf, nclosed+1, nclosed+nopen,
     $           coeff, nbf)

            if (eaf_close(pdmfile).ne.0)
     $           call errquit('mcscf: failed closing 2pdm file',0,
     &                        DISK_ERR)

*            call MA_summarize_allocated_blocks()
            
            if (.not. ma_pop_stack(l_buf)) 
     $           call errquit('gradients: ma corrupt buf',0, MA_ERR)

c
         end if
         call ga_sync()
c
         call ga_brdcst(33,pdm2, 8*nopen**4,0)
         call ga_brdcst(34,coeff,8*nopen*nbf,0)

         call ga_sync()
c
         call ga_symmetrize(g_wdens)
         call ga_symmetrize(g_dens(1))

      end if
c     
      status = ga_destroy ( g_vecs )
      if(lfrac_occ) then
        status = ga_destroy ( g_tmp )
        if ( scftype .eq. 'UHF' ) then
         status = ga_destroy ( g_tmp2 )
        end if
      endif
      if ( scftype .eq. 'UHF' ) then
         status = ga_destroy ( g_vecs2 )
      end if
c
      return
      end
c
      subroutine grad_diagb(g_a, evals, nocc)
      implicit none
#include "global.fh"
      integer g_a
      double precision evals(*)
      integer nocc
c
      integer k
      integer ilo, ihi, jlo, jhi
c
      call ga_distribution(g_a,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and.jlo.gt.0.and.
     A    ilo.le.ihi.and.jlo.le.ihi) then
         do k=1,nocc
            if(k.ge.ilo.and.k.le.ihi.and.k.ge.jlo.and.k.le.jhi)
     P           call ga_put(g_a,k,k,k,k,evals(k),1)
         enddo
      endif
      return
      end
