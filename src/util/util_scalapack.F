      subroutine util_diag_std(g_fock,g_evecs,evals)
      implicit none
#include "global.fh"
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
C#ifdef USE_ELPA
C      call ga_evp_real(g_fock, g_evecs, evals,nb, mout)
C#else
#ifdef USE_PDSYEVR
      if(ga_nnodes().gt.1) then
         call ga_pdsyevr(g_fock, g_evecs, evals, nb, mout)
      else
         call ga_pdsyevd(g_fock, g_evecs, evals, nb)
      endif
#else
      call ga_pdsyevd(g_fock, g_evecs, evals, nb)
#endif
C#endif
#else
      call ga_diag_std(g_fock, g_evecs, evals)
#endif
#else
      call ga_diag_std_seq(g_fock, g_evecs, evals)
#endif
      return
      end
