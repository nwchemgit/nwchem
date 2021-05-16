c     
c     
c***********************************************************************
c
c-->  Construct the density matrix.
c
c***********************************************************************
c
      subroutine dft_densm_so(g_dens,g_vecs,nbf,noc)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_vecs(2)
      integer nbf, nbf2
      integer noc
c
      double precision zero, one, mone
      parameter (zero=0.d0, one=1.d0, mone=-1.d0)
c
      nbf2 = 2*nbf 
c
      call ga_dgemm('n', 't', nbf2, nbf2, noc, one,  g_vecs(1),
     $  g_vecs(1), zero, g_dens(1))
      call ga_dgemm('n', 't', nbf2, nbf2, noc, one,  g_vecs(2),
     $  g_vecs(2), one,  g_dens(1)) 
      call ga_dgemm('n', 't', nbf2, nbf2, noc, one,  g_vecs(1),
     $  g_vecs(2), zero, g_dens(2))
      call ga_dgemm('n', 't', nbf2, nbf2, noc, mone, g_vecs(2),
     $  g_vecs(1), one,  g_dens(2)) 
c
      return
      end 

      subroutine ga_fock_sf(g_fock, g_fs, nbf)
C$Id$
      implicit none

      integer g_fock
      integer g_fs
      integer nbf, npone, ntwo 
c
      npone=nbf+1 
      ntwo=2*nbf
      call ga_dadd_patch( 1.d0, g_fock, 1, nbf, 1, nbf,
     &                    1.d0, g_fs,1, nbf, 1, nbf,
     &                          g_fs,1, nbf, 1, nbf)
      call ga_dadd_patch(1.d0, g_fock, 1, nbf, 1, nbf,
     &                   1.d0, g_fs,npone, ntwo,npone, ntwo,
     &                         g_fs,npone, ntwo,npone, ntwo)
      return
      end

      subroutine ga_fock_so(g_so, g_fs, nbf)
C$Id$
      implicit none

      integer g_so(3)
      integer g_fs(2)
      integer nbf, obf, tbf
c
      obf = nbf+1
      tbf = 2*nbf
      call ga_dadd_patch(1.d0,  g_so(1), 1, nbf, 1, nbf,
     &                   1.d0,  g_fs(2), 1, nbf, 1, nbf,
     &                          g_fs(2), 1, nbf, 1, nbf)
      call ga_dadd_patch(-1.d0, g_so(1), 1, nbf, 1, nbf,
     &                    1.d0, g_fs(2), obf, tbf, obf, tbf,
     &                          g_fs(2), obf, tbf, obf, tbf)

      call ga_dadd_patch(1.d0,  g_so(2), 1, nbf, 1, nbf,
     &                   1.d0,  g_fs(1), 1, nbf, obf, tbf,
     &                          g_fs(1), 1, nbf, obf, tbf)
      call ga_dadd_patch(-1.d0, g_so(2), 1, nbf, 1, nbf,
     &                   1.d0,  g_fs(1), obf, tbf, 1, nbf,
     &                          g_fs(1), obf, tbf, 1, nbf)

      call ga_dadd_patch(1.d0,  g_so(3), 1, nbf, 1, nbf,
     &                   1.d0,  g_fs(2), 1, nbf, obf, tbf,
     &                          g_fs(2), 1, nbf, obf, tbf)
      call ga_dadd_patch(1.d0,  g_so(3), 1, nbf, 1, nbf,
     &                   1.d0,  g_fs(2), obf, tbf, 1, nbf,
     &                          g_fs(2), obf, tbf, 1, nbf)

      return
      end

      subroutine ga_dens_so(g_tmp, g_a, nbf, op)
C$Id$
      implicit none

      integer g_tmp
      integer g_a(2)
      integer nbf, npone, ntwice
      character*1 op
c
      npone=nbf+1 
      ntwice=2*nbf
      call ga_zero(g_tmp)
      if (op .eq. 'z') then
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                     -1.d0, g_a(2), 1, nbf, 1, nbf,
     &                            g_tmp, 1, nbf, 1, nbf)
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                      1.d0, g_a(2),npone,ntwice,npone,ntwice, 
     &                            g_tmp, 1, nbf, 1, nbf)

      else if (op .eq. 'x') then 
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                     -1.d0, g_a(2), 1, nbf, npone,ntwice,
     &                            g_tmp, 1, nbf, 1, nbf)
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                     -1.d0, g_a(2),npone,ntwice,1,nbf,
     &                            g_tmp, 1, nbf, 1, nbf)

      else if (op .eq. 'y') then 
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                      1.d0, g_a(1), 1, nbf, npone,ntwice,
     &                            g_tmp, 1, nbf, 1, nbf)
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                     -1.d0, g_a(1),npone,ntwice,1,nbf,
     &                            g_tmp, 1, nbf, 1, nbf)

      else if (op. eq. 's') then 
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                      1.d0, g_a(1), 1, nbf, 1, nbf,
     &                            g_tmp, 1, nbf, 1, nbf)
         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
     &                      1.d0, g_a(1),npone,ntwice,npone,ntwice, 
     &                            g_tmp, 1, nbf, 1, nbf)
      else
         call errquit('ga_dens_so: op arg not valid',0,0)
      endif
      return
      end
 
      subroutine ga_dens_sf(g_dens, g_a, nbf)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_a(2)
      integer nbf, npone, ntwice
c
      npone=nbf+1 
      ntwice=2*nbf
c    
      call ga_copy_patch('n',g_a(1), 1, nbf, 1, nbf,
     &                         g_dens(1), 1, nbf, 1, nbf)
      call ga_copy_patch('n',g_a(1),npone,ntwice,npone,ntwice, 
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end
 
      subroutine dft_densm_sowght(g_dens,g_vecs,g_eigen_diag,nbf,noc)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_vecs(2)
      integer g_eigen_diag 
      integer nbf, nbf2
      integer noc
c
      double precision zero, one, mone
      parameter (zero=0.d0, one=1.d0, mone=-1.d0)
c
      nbf2 = 2*nbf 
cso   only real density matrix is needed for weighted density matrix
cso   
cso   real movecs contribution 
cso
         call ga_matmul_patch ('n','n', one, zero,
     $        g_vecs(1),       1, nbf2,  1, noc,
     $        g_eigen_diag, 1, noc, 1, noc,
     $        g_dens(2),    1, nbf2,  1, noc )
         call ga_matmul_patch ('n','t', one, zero,
     $        g_dens(2),    1, nbf2,  1, noc,
     $        g_vecs(1),       1, noc, 1, nbf2,
     $        g_dens(1),      1, nbf2,  1, nbf2 )
 
cso   
cso   imaginary movecs contribution 
cso
         call ga_matmul_patch ('n','n', one, zero,
     $        g_vecs(2),       1, nbf2,  1, noc,
     $        g_eigen_diag, 1, noc, 1, noc,
     $        g_dens(2),    1, nbf2,  1, noc )
         call ga_matmul_patch ('n','t', one, one,
     $        g_dens(2),    1, nbf2,  1, noc,
     $        g_vecs(2),       1, noc, 1, nbf2,
     $        g_dens(1),      1, nbf2,  1, nbf2 )

cso      call ga_symmetrize(g_dens(1))
cso      call ga_symmetrize(g_dens(2))
      return
      end 
c
      subroutine ga_dens_aabbr(g_dens, g_densso, nbf)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_densso(2)
      integer nbf, npone, ntwice
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c
      npone=nbf+1
      ntwice=2*nbf
c
c     == alpha (real) ==
      call ga_dadd_patch(0.d0, g_dens(1), 1, nbf, 1, nbf, 
     &                   1.d0, g_densso(1),  1, nbf, 1, nbf,
     &                         g_dens(1), 1, nbf, 1, nbf)
c     == beta (real) ==
      call ga_dadd_patch(0.d0, g_dens(2), 1, nbf, 1, nbf, 
     &                   1.d0, g_densso(1), npone,ntwice,npone,ntwice, 
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end
c
      subroutine ga_dens_aabbi(g_dens, g_densso, nbf)
C     $Id$
      implicit none
      
      integer g_dens(2)
      integer g_densso(2)
      integer nbf, npone, ntwice
c     
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c
      npone=nbf+1 
      ntwice=2*nbf
c     
c     == alpha (imag) ==
      call ga_dadd_patch(0.d0, g_dens(1), 1, nbf, 1, nbf, 
     &                   1.d0, g_densso(2),1, nbf, 1, nbf,
     &                         g_dens(1), 1, nbf, 1, nbf)
c     == beta (imag) ==
      call ga_dadd_patch(0.d0, g_dens(2), 1, nbf, 1, nbf, 
     &                   1.d0, g_densso(2),npone,ntwice,npone,ntwice, 
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end
c
      subroutine ga_dens_abbar(g_dens, g_a, nbf)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_a(2)
      integer nbf, npone, ntwice
c
      npone=nbf+1 
      ntwice=2*nbf
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c
      call ga_dadd_patch(1.d0, g_dens(1), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(1),    1, nbf, npone,ntwice,
     &                         g_dens(1),    1, nbf, 1, nbf)
      call ga_dadd_patch(1.d0, g_dens(2), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(1),    npone, ntwice, 1, nbf, 
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end

      subroutine ga_dens_abbai(g_dens, g_a, nbf)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_a(2)
      integer nbf, npone, ntwice
c
      npone=nbf+1 
      ntwice=2*nbf
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c
      call ga_dadd_patch(1.d0, g_dens(1), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(2),    1,nbf,  npone,ntwice, 
     &                         g_dens(1),    1, nbf, 1, nbf)
      call ga_dadd_patch(1.d0, g_dens(2), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(2),    npone,ntwice, 1,nbf,
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end

      subroutine ga_dens_abri(g_dens, g_a, nbf)
C$Id$
      implicit none

      integer g_dens(2)
      integer g_a(2)
      integer nbf, npone, ntwice
c
      npone=nbf+1 
      ntwice=2*nbf
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c
      call ga_dadd_patch(1.d0, g_dens(1), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(1),    1, nbf, npone,ntwice, 
     &                         g_dens(1), 1, nbf, 1, nbf)
      call ga_dadd_patch(1.d0, g_dens(2), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(2),    1, nbf, npone,ntwice, 
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end


      subroutine dftv_densm_so(dens,vecs,nbf_mo,noc)
C$Id$
      implicit none

      integer nbf_mo
      double precision dens(nbf_mo,nbf_mo,2)
      double precision vecs(nbf_mo,2)
      integer noc
c
      integer j
c

      call dgemm('n', 't', nbf_mo, nbf_mo, 1, 1d0, vecs(1,1), nbf_mo,
     $  vecs(1,1), nbf_mo,0d0, dens(1,1,1),nbf_mo)
      call dgemm('n', 't', nbf_mo, nbf_mo, 1, 1d0, vecs(1,2), nbf_mo,
     $  vecs(1,2), nbf_mo, 1d0,  dens(1,1,1), nbf_mo) 
      call dgemm('n', 't', nbf_mo, nbf_mo, 1, 1d0, vecs(1,1), nbf_mo,
     $  vecs(1,2), nbf_mo, 0d0, dens(1,1,2), nbf_mo)
      call dgemm('n', 't', nbf_mo, nbf_mo, 1,-1d0, vecs(1,2),nbf_mo,
     $  vecs(1,1), nbf_mo, 1d0,  dens(1,1,2), nbf_mo) 
      return
      end 
      subroutine vec_dens_so(vtmp, a, nbf, op)
C$Id$
      implicit none

      integer nbf
      double precision vtmp(nbf,nbf)
      integer a(2*nbf,2*nbf,2)
      character*1 op
c
      integer j
c      call ga_zero(g_tmp)
      call dcopy(nbf*nbf,0d0,0,vtmp,1)

      if (op .eq. 'z') then
         do j=1,nbf
            call daxpy(nbf,-1d0,a(1    ,j,2),1,vtmp(1,j),1)
         enddo
         do j=1,nbf
            call daxpy(nbf,1d0,a(1+nbf,j+nbf,2),1,vtmp(1,j),1)
         enddo

      else if (op .eq. 'x') then 
         do j=1,nbf
            call daxpy(nbf,-1d0,a(1,j+nbf,2),1,vtmp(1,j),1)
         enddo
         do j=1,nbf
            call daxpy(nbf,-1d0,a(1+nbf,j,2),1,vtmp(1,j),1)
         enddo

      else if (op .eq. 'y') then 
         do j=1,nbf
            call daxpy(nbf,1d0,a(1,j+nbf,1),1,vtmp(1,j),1)
         enddo
         do j=1,nbf
            call daxpy(nbf,-1d0,a(1+nbf,j,1),1,vtmp(1,j),1)
         enddo

      else if (op. eq. 's') then 
         do j=1,nbf
            call daxpy(nbf,1d0,a(1,j,1),1,vtmp(1,j),1)
         enddo
c         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
c     &                      1.d0, g_a(1), 1, nbf, 1, nbf,
c     &                            g_tmp, 1, nbf, 1, nbf)
         do j=1,nbf
            call daxpy(nbf,1d0,a(1+nbf,j+nbf,1),1,vtmp(1,j),1)
         enddo
c         call ga_dadd_patch(1.d0, g_tmp, 1, nbf, 1, nbf, 
c     &                      1.d0, g_a(1),1+nbf,2*nbf,1+nbf,2*nbf, 
c     &                            g_tmp, 1, nbf, 1, nbf)
      else
         call errquit(' vec_dens_so: op arg not valid',0,0)
      endif
      return
      end
