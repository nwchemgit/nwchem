      subroutine util_diag_std(g_fock,g_evecs,evals)
      implicit none
      integer g_fock            ! matrix to diagonalize
      integer g_evecs           ! eigenvectors
      double precision evals(*) ! eigenvalues
c
      integer nb
      integer mout
      mout=0
      nb=0
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
#ifdef USE_ELPA
      call ga_evp_real(g_fock, g_evecs, evals,nb, mout)
#else      
#ifdef USE_PDSYEVR
      call ga_pdsyevr(g_fock, g_evecs, evals, nb, mout)
#else      
      call ga_pdsyevd(g_fock, g_evecs, evals, nb)
#endif
#endif      
#else
      call ga_diag_std(g_fock, g_evecs, evals)
#endif
#else 
      call ga_diag_std_seq(g_fock, g_evecs, evals)
#endif
      return
      end
