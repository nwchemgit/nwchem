c***********************************************************************
c
c-->  Construct the density matrix.
c
c***********************************************************************
c
      subroutine dft_densm(g_dens,g_vecs,nbf,ipol,noc)
C$Id: dft_densm.F,v 1.4 1997-06-11 18:16:50 d3h449 Exp $
      implicit none

      integer g_dens
      integer g_vecs
      integer nbf
      integer ipol
      integer noc
c
      double precision zero,fact
      parameter (zero=0.d0)
c
      fact=2.d0/dble(ipol)
      call ga_dgemm('n', 't', nbf, nbf, noc, fact, g_vecs,
     $  g_vecs, zero, g_dens)
      return
      end
