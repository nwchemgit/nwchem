      subroutine dftg_get_xden( d_ik, d_jl, d_il, d_jk, ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
c  $Id$
C
C     get atomic blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none
      
      integer ldim
      integer iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l
      integer g_dens

      double precision d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim)
c
      integer li,lj,lk,ll
c

      logical eij, ejk, ekl, eil, eik, ejl

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

C     i,k
        call ga_get (g_dens, iab1f, iab1l, iab3f, iab3l, d_ik,ldim)

C     j,l
      if ( eij.and.ekl ) then
        call grad_copy ( d_ik, d_jl, li, lk, ldim )
      else if ( ejk.and.eil ) then
        call grad_transp ( d_ik, d_jl, li, lk, ldim )
      else
        call ga_get (g_dens, iab2f, iab2l, iab4f, iab4l, d_jl,ldim)
      endif

C     i,l
      if ( eij ) then
        call grad_copy ( d_jl, d_il, lj, ll, ldim )
      else if ( ekl ) then
        call grad_copy ( d_ik, d_il, li, lk, ldim )
      else
        call ga_get (g_dens, iab1f, iab1l, iab4f, iab4l, d_il,ldim)
      endif

C     j,k
      if ( ekl ) then
        call grad_copy ( d_jl, d_jk, lj, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_ik, d_jk, li, lk, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_il, d_jk, li, ll, ldim )
      else if ( ejl.and.eik ) then
        call grad_transp ( d_il, d_jk, li, ll, ldim )
      else
        call ga_get (g_dens, iab2f, iab2l, iab3f, iab3l, d_jk,ldim)
      endif

      return
      end
