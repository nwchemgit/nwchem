      Subroutine dft_swapab(g_dens, g_movecs,nbf_ao,nocc,
     ,        evals,occ,g_tmp,ktmp)

C$Id$
      Implicit none
      integer nbf_ao              ! [input] number of ao basis functions
      integer g_dens(2)           ! dens       [input/out] 
      integer g_movecs(2)         ! eigenvecs [input/output]
      integer g_tmp            ! scratch space
      double precision evals(nbf_ao,2)! eigenvalues  [input/output]
      double precision occ(nbf_ao,2)! occup nos  [input/output]
      double precision ktmp(*)! [scratch]
      integer nocc(2)
c
c     local
c
      integer scsc
c
c     DM
c
      call ga_copy(g_dens(2),g_tmp)
      call ga_copy(g_dens(1),g_dens(2))
      call ga_copy(g_tmp,g_dens(1))
c
c     MOs
c
      call ga_copy(g_movecs(2),g_tmp)
      call ga_copy(g_movecs(1),g_movecs(2))
      call ga_copy(g_tmp,g_movecs(1))
c
c     evals
c
      call dcopy(nbf_ao,evals,1,ktmp,1)
      call dcopy(nbf_ao,evals(1,2),1,evals(1,1),1)
      call dcopy(nbf_ao,ktmp,1,evals(1,2),1)
c
c     evals
c
      call dcopy(nbf_ao,occ,1,ktmp,1)
      call dcopy(nbf_ao,occ(1,2),1,occ(1,1),1)
      call dcopy(nbf_ao,ktmp,1,occ(1,2),1)
c
      scsc=nocc(1)
      nocc(1)=nocc(2)
      nocc(2)=scsc
      return
      end
