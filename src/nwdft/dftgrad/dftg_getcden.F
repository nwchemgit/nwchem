      subroutine dftg_get_cden( d_ij, d_kl, ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l,
     &           iab4f, iab4l,
     $           g_dens )
C
c**** adapted from the scf gradients routine get_den.  
C     get atomic blocks from density matrix and make sure no block is
C     transferred twice
C
*
* $Id$
*
      implicit none
      
      integer ldim
      integer iab1f, iab1l, iab2f, iab2l, iab3f, iab3l,
     &           iab4f, iab4l

      double precision d_ij(ldim,ldim),d_kl(ldim,ldim)
      integer g_dens

      logical eij, ejk, ekl, eil, eik, ejl
c
      integer li,lj,lk,ll
c

      eij = ( iab1f.eq.iab2f )
      ejk = ( iab2f.eq.iab3f )
      ekl = ( iab3f.eq.iab4f )
      eil = ( iab1f.eq.iab4f )
      eik = ( iab1f.eq.iab3f )
      ejl = ( iab2f.eq.iab4f )

      li = iab1l - iab1f + 1
      lj = iab2l - iab2f + 1
      lk = iab3l - iab3f + 1
      ll = iab4l - iab4f + 1

C     i,j
      call ga_get (g_dens, iab1f, iab1l, iab2f, iab2l, d_ij,ldim)

C     k,l
      if ( eik.and.ejl ) then
        call grad_copy ( d_ij, d_kl, li, lj, ldim )
      else if ( eil.and.ejk ) then
        call grad_transp ( d_ij, d_kl, li, lj, ldim )
      else
        call ga_get (g_dens, iab3f, iab3l, iab4f, iab4l, d_kl,ldim)
      endif

      return
      end
