c     
c     
c***********************************************************************
c
c-->  Construct the density matrix.
c
c***********************************************************************
c
      subroutine dft_densm_so(g_dens,g_vecs,nbf,noc)
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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

      call ga_dgemm('n', 't', nbf2, nbf2, noc, one,  g_vecs(1),
     $  g_vecs(1), zero, g_dens(1))
      call ga_dgemm('n', 't', nbf2, nbf2, noc, one,  g_vecs(2),
     $  g_vecs(2), one,  g_dens(1)) 
      call ga_dgemm('n', 't', nbf2, nbf2, noc, one,  g_vecs(1),
     $  g_vecs(2), zero, g_dens(2))
      call ga_dgemm('n', 't', nbf2, nbf2, noc, mone, g_vecs(2),
     $  g_vecs(1), one,  g_dens(2)) 

cso      call ga_symmetrize(g_dens(1))
cso      call ga_symmetrize(g_dens(2))
      return
      end 

      subroutine ga_fock_sf(g_fock, g_fs, nbf)
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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
         stop
      endif
      return
      end
 
      subroutine ga_dens_sf(g_dens, g_a, nbf)
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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

      subroutine ga_dens_aabbi(g_dens, g_a, nbf)
C     $Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
      implicit none
      
      integer g_dens(2)
      integer g_a(2)
      integer nbf, npone, ntwice
c     
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c
      npone=nbf+1 
      ntwice=2*nbf
c     
      call ga_dadd_patch(0.d0, g_dens(1), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(2),    1, nbf, 1, nbf,
     &                         g_dens(1), 1, nbf, 1, nbf)
      call ga_dadd_patch(0.d0, g_dens(2), 1, nbf, 1, nbf, 
     &                   1.d0, g_a(2),    npone,ntwice,npone,ntwice, 
     &                         g_dens(2), 1, nbf, 1, nbf)
      return
      end

      subroutine ga_dens_abbar(g_dens, g_a, nbf)
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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
C$Id: dft_densmso.F,v 1.7 2007-07-02 17:57:53 niri Exp $
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

