      subroutine grad_get_dens( d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
c$Id: grad_getdens.F,v 1.1 1995-12-13 01:35:16 d3g681 Exp $
C
C     get atomic blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, ldim, 
     $     iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $     li, lj, lk, ll

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk

      dimension d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim)

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

C     i,j
      call ga_get (g_dens, iab1f, iab1l, iab2f, iab2l, d_ij, ldim)

C     k,l
      if ( eik.and.ejl ) then
        call grad_copy ( d_ij, d_kl, li, lj, ldim )
      else if ( eil.and.ejk ) then
        call grad_transp ( d_ij, d_kl, li, lj, ldim )
      else
        call ga_get (g_dens, iab3f, iab3l, iab4f, iab4l, d_kl, ldim)
      endif

C     i,k
      if ( ejk ) then
        call grad_copy ( d_ij, d_ik, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_kl, d_ik, lk, ll, ldim )
      else
        call ga_get (g_dens, iab1f, iab1l, iab3f, iab3l, d_ik, ldim)
      endif

C     j,l
      if ( ejk ) then
        call grad_copy ( d_kl, d_jl, lk, ll, ldim )
      else if ( eij.and.ejl ) then
        call grad_copy ( d_ij, d_jl, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_ij, d_jl, li, lj, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_ik, d_jl, li, lk, ldim )
      else if ( ejk.and.eil ) then
        call grad_transp ( d_ik, d_jl, li, lk, ldim )
      else
        call ga_get (g_dens, iab2f, iab2l, iab4f, iab4l, d_jl, ldim)
      endif

C     i,l
      if ( ejl ) then
        call grad_copy ( d_ij, d_il, li, lj, ldim )
      else if ( eik ) then
        call grad_copy ( d_kl, d_il, lk, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_jl, d_il, lj, ll, ldim )
      else if ( ekl ) then
        call grad_copy ( d_ik, d_il, li, lk, ldim )
      else
        call ga_get (g_dens, iab1f, iab1l, iab4f, iab4l, d_il, ldim)
      endif

C     j,k
      if ( ekl ) then
        call grad_copy ( d_jl, d_jk, lj, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_ik, d_jk, li, lk, ldim )
      else if ( eik ) then
        call grad_transp ( d_ij, d_jk, li, lj, ldim )
      else if ( ejl ) then
        call grad_transp ( d_kl, d_jk, lk, ll, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_il, d_jk, li, ll, ldim )
      else if ( ejl.and.eik ) then
        call grad_transp ( d_il, d_jk, li, ll, ldim )
      else
        call ga_get (g_dens, iab2f, iab2l, iab3f, iab3l, d_jk, ldim)
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_rohf( d_ij, d_kl, d_ik, d_jl, d_il, d_jk,  
     $           d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
C
C     get atomic blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, ldim, 
     $     iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $     li, lj, lk, ll

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $          d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2

      dimension d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim),
     $          d_ij2(ldim,ldim),d_kl2(ldim,ldim),d_ik2(ldim,ldim),
     $          d_jl2(ldim,ldim),d_il2(ldim,ldim),d_jk2(ldim,ldim)
      dimension g_dens(3)

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

C     i,j
      call ga_get (g_dens(2), iab1f, iab1l, iab2f, iab2l, d_ij, ldim)
      call ga_get (g_dens(3), iab1f, iab1l, iab2f, iab2l, d_ij2, ldim)

C     k,l
      if ( eik.and.ejl ) then
        call grad_copy ( d_ij, d_kl, li, lj, ldim )
        call grad_copy ( d_ij2, d_kl2, li, lj, ldim )
      else if ( eil.and.ejk ) then
        call grad_transp ( d_ij, d_kl, li, lj, ldim )
        call grad_transp ( d_ij2, d_kl2, li, lj, ldim )
      else
        call ga_get (g_dens(2), iab3f, iab3l, iab4f, iab4l, d_kl, ldim)
        call ga_get (g_dens(3), iab3f, iab3l, iab4f, iab4l, d_kl2, ldim)
      endif

C     i,k
      if ( ejk ) then
        call grad_copy ( d_ij, d_ik, li, lj, ldim )
        call grad_copy ( d_ij2, d_ik2, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_kl, d_ik, lk, ll, ldim )
        call grad_transp ( d_kl2, d_ik2, lk, ll, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab3f, iab3l, d_ik, ldim)
        call ga_get (g_dens(3), iab1f, iab1l, iab3f, iab3l, d_ik2, ldim)
      endif

C     j,l
      if ( ejk ) then
        call grad_copy ( d_kl, d_jl, lk, ll, ldim )
        call grad_copy ( d_kl2, d_jl2, lk, ll, ldim )
      else if ( eij.and.ejl ) then
        call grad_copy ( d_ij, d_jl, li, lj, ldim )
        call grad_copy ( d_ij2, d_jl2, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_ij, d_jl, li, lj, ldim )
        call grad_transp ( d_ij2, d_jl2, li, lj, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_ik, d_jl, li, lk, ldim )
        call grad_copy ( d_ik2, d_jl2, li, lk, ldim )
      else if ( ejk.and.eil ) then
        call grad_transp ( d_ik, d_jl, li, lk, ldim )
        call grad_transp ( d_ik2, d_jl2, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab4f, iab4l, d_jl, ldim)
        call ga_get (g_dens(3), iab2f, iab2l, iab4f, iab4l, d_jl2, ldim)
      endif

C     i,l
      if ( ejl ) then
        call grad_copy ( d_ij, d_il, li, lj, ldim )
        call grad_copy ( d_ij2, d_il2, li, lj, ldim )
      else if ( eik ) then
        call grad_copy ( d_kl, d_il, lk, ll, ldim )
        call grad_copy ( d_kl2, d_il2, lk, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_jl, d_il, lj, ll, ldim )
        call grad_copy ( d_jl2, d_il2, lj, ll, ldim )
      else if ( ekl ) then
        call grad_copy ( d_ik, d_il, li, lk, ldim )
        call grad_copy ( d_ik2, d_il2, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab4f, iab4l, d_il, ldim)
        call ga_get (g_dens(3), iab1f, iab1l, iab4f, iab4l, d_il2, ldim)
      endif

C     j,k
      if ( ekl ) then
        call grad_copy ( d_jl, d_jk, lj, ll, ldim )
        call grad_copy ( d_jl2, d_jk2, lj, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_ik, d_jk, li, lk, ldim )
        call grad_copy ( d_ik2, d_jk2, li, lk, ldim )
      else if ( eik ) then
        call grad_transp ( d_ij, d_jk, li, lj, ldim )
        call grad_transp ( d_ij2, d_jk2, li, lj, ldim )
      else if ( ejl ) then
        call grad_transp ( d_kl, d_jk, lk, ll, ldim )
        call grad_transp ( d_kl2, d_jk2, lk, ll, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_il, d_jk, li, ll, ldim )
        call grad_copy ( d_il2, d_jk2, li, ll, ldim )
      else if ( ejl.and.eik ) then
        call grad_transp ( d_il, d_jk, li, ll, ldim )
        call grad_transp ( d_il2, d_jk2, li, ll, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab3f, iab3l, d_jk, ldim)
        call ga_get (g_dens(3), iab2f, iab2l, iab3f, iab3l, d_jk2, ldim)
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_uhf( d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           d_ik2, d_jl2, d_il2, d_jk2, ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
C
C     get atomic blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, ldim, 
     $     iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $     li, lj, lk, ll

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     d_ik2, d_jl2, d_il2, d_jk2

      dimension d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim),
     $          d_ik2(ldim,ldim),d_jl2(ldim,ldim),d_il2(ldim,ldim),
     $          d_jk2(ldim,ldim)
      dimension g_dens(2)

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

C     i,j
      call ga_get (g_dens(1), iab1f, iab1l, iab2f, iab2l, d_ij, ldim)

C     k,l
      if ( eik.and.ejl ) then
        call grad_copy ( d_ij, d_kl, li, lj, ldim )
      else if ( eil.and.ejk ) then
        call grad_transp ( d_ij, d_kl, li, lj, ldim )
      else
        call ga_get (g_dens(1), iab3f, iab3l, iab4f, iab4l, d_kl, ldim)
      endif

C     i,k
      if ( ejk ) then
        call grad_copy ( d_ij, d_ik, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_kl, d_ik, lk, ll, ldim )
      else
        call ga_get (g_dens(1), iab1f, iab1l, iab3f, iab3l, d_ik, ldim)
      endif
      call ga_get (g_dens(2), iab1f, iab1l, iab3f, iab3l, d_ik2, ldim)

C     j,l
      if ( ejk ) then
        call grad_copy ( d_kl, d_jl, lk, ll, ldim )
      else if ( eij.and.ejl ) then
        call grad_copy ( d_ij, d_jl, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_ij, d_jl, li, lj, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_ik, d_jl, li, lk, ldim )
      else if ( ejk.and.eil ) then
        call grad_transp ( d_ik, d_jl, li, lk, ldim )
      else
        call ga_get (g_dens(1), iab2f, iab2l, iab4f, iab4l, d_jl, ldim)
      endif


      if ( eij .and. ekl ) then
        call grad_copy ( d_ik2, d_jl2, li, lk, ldim )
      else if ( eil .and. ejk ) then
        call grad_transp ( d_ik2, d_jl2, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab4f, iab4l, d_jl2, ldim)
      endif

C     i,l
      if ( ejl ) then
        call grad_copy ( d_ij, d_il, li, lj, ldim )
      else if ( eik ) then
        call grad_copy ( d_kl, d_il, lk, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_jl, d_il, lj, ll, ldim )
      else if ( ekl ) then
        call grad_copy ( d_ik, d_il, li, lk, ldim )
      else
        call ga_get (g_dens(1), iab1f, iab1l, iab4f, iab4l, d_il, ldim)
      endif

      if ( ekl ) then
        call grad_copy ( d_ik2, d_il2, li, lk, ldim )
      else if ( eij ) then
        call grad_copy ( d_jl2, d_il2, lj, ll, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab4f, iab4l, d_il2, ldim)
      endif

C     j,k
      if ( ekl ) then
        call grad_copy ( d_jl, d_jk, lj, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_ik, d_jk, li, lk, ldim )
      else if ( eik ) then
        call grad_transp ( d_ij, d_jk, li, lj, ldim )
      else if ( ejl ) then
        call grad_transp ( d_kl, d_jk, lk, ll, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_il, d_jk, li, ll, ldim )
      else if ( ejl.and.eik ) then
        call grad_transp ( d_il, d_jk, li, ll, ldim )
      else
        call ga_get (g_dens(1), iab2f, iab2l, iab3f, iab3l, d_jk, ldim)
      endif

      if ( eij ) then
        call grad_copy ( d_ik2, d_jk2, li, lk, ldim )
      else if ( ekl ) then
        call grad_copy ( d_jl2, d_jk2, lj, ll, ldim )
      else if ( eik .and. ejl ) then
        call grad_transp ( d_il2, d_jk2, li, ll, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab3f, iab3l, d_jk2, ldim)
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_copy ( d_1, d_2, l1, l2, ld )
      implicit none
      integer l1, l2, ld, i, j
      double precision d_1, d_2
      dimension d_1(ld,ld), d_2(ld,ld)
      do 20, j=1, l2
        do 10, i=1, l1
          d_2(i,j) = d_1(i,j)
 10     continue
 20   continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_transp ( d_1, d_2, l1, l2, ld )
      implicit none
      integer l1, l2, ld, i, j
      double precision d_1, d_2
      dimension d_1(ld,ld), d_2(ld,ld)
      do 20, j=1, l2
        do 10, i=1, l1
          d_2(j,i) = d_1(i,j)
 10     continue
 20   continue
      return
      end
