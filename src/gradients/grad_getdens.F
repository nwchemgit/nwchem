      subroutine grad_get_dens( d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
c$Id: grad_getdens.F,v 1.3 1996-10-21 00:39:36 vg038 Exp $
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_rmp2 ( 
     $           d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           p_ij, p_kl, p_ik, p_jl, p_il, p_jk, 
     $           ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
C
C     get atomic blocks from density matrices and make sure no block is
C     transferred twice
C
C     (a,b) denotes spin, d is the HF part and p contains (twice) the 
C     second order correction
C
      implicit none

#include "global.fh"

      integer g_dens, ldim, 
     $     iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $     li, lj, lk, ll

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk,
     $                 p_ij, p_kl, p_ik, p_jl, p_il, p_jk

      dimension g_dens(3)
      dimension d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim),
     $          p_ij(ldim,ldim),p_kl(ldim,ldim),p_ik(ldim,ldim),
     $          p_jl(ldim,ldim),p_il(ldim,ldim),p_jk(ldim,ldim)

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
      call ga_get (g_dens(3), iab1f, iab1l, iab2f, iab2l, p_ij, ldim)

C     k,l
      if ( eik.and.ejl ) then
        call grad_copy ( d_ij, d_kl, li, lj, ldim )
        call grad_copy ( p_ij, p_kl, li, lj, ldim )
      else if ( eil.and.ejk ) then
        call grad_transp ( d_ij, d_kl, li, lj, ldim )
        call grad_transp ( p_ij, p_kl, li, lj, ldim )
      else
        call ga_get (g_dens(2), iab3f, iab3l, iab4f, iab4l, d_kl, ldim)
        call ga_get (g_dens(3), iab3f, iab3l, iab4f, iab4l, p_kl, ldim)
      endif

C     i,k
      if ( ejk ) then
        call grad_copy ( d_ij, d_ik, li, lj, ldim )
        call grad_copy ( p_ij, p_ik, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_kl, d_ik, lk, ll, ldim )
        call grad_transp ( p_kl, p_ik, lk, ll, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab3f, iab3l, d_ik, ldim)
        call ga_get (g_dens(3), iab1f, iab1l, iab3f, iab3l, p_ik, ldim)
      endif

C     j,l
      if ( ejk ) then
        call grad_copy ( d_kl, d_jl, lk, ll, ldim )
        call grad_copy ( p_kl, p_jl, lk, ll, ldim )
      else if ( eij.and.ejl ) then
        call grad_copy ( d_ij, d_jl, li, lj, ldim )
        call grad_copy ( p_ij, p_jl, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( d_ij, d_jl, li, lj, ldim )
        call grad_transp ( p_ij, p_jl, li, lj, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_ik, d_jl, li, lk, ldim )
        call grad_copy ( p_ik, p_jl, li, lk, ldim )
      else if ( ejk.and.eil ) then
        call grad_transp ( d_ik, d_jl, li, lk, ldim )
        call grad_transp ( p_ik, p_jl, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab4f, iab4l, d_jl, ldim)
        call ga_get (g_dens(3), iab2f, iab2l, iab4f, iab4l, p_jl, ldim)
      endif

C     i,l
      if ( ejl ) then
        call grad_copy ( d_ij, d_il, li, lj, ldim )
        call grad_copy ( p_ij, p_il, li, lj, ldim )
      else if ( eik ) then
        call grad_copy ( d_kl, d_il, lk, ll, ldim )
        call grad_copy ( p_kl, p_il, lk, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_jl, d_il, lj, ll, ldim )
        call grad_copy ( p_jl, p_il, lj, ll, ldim )
      else if ( ekl ) then
        call grad_copy ( d_ik, d_il, li, lk, ldim )
        call grad_copy ( p_ik, p_il, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab4f, iab4l, d_il, ldim)
        call ga_get (g_dens(3), iab1f, iab1l, iab4f, iab4l, p_il, ldim)
      endif

C     j,k
      if ( ekl ) then
        call grad_copy ( d_jl, d_jk, lj, ll, ldim )
        call grad_copy ( p_jl, p_jk, lj, ll, ldim )
      else if ( eij ) then
        call grad_copy ( d_ik, d_jk, li, lk, ldim )
        call grad_copy ( p_ik, p_jk, li, lk, ldim )
      else if ( eik ) then
        call grad_transp ( d_ij, d_jk, li, lj, ldim )
        call grad_transp ( p_ij, p_jk, li, lj, ldim )
      else if ( ejl ) then
        call grad_transp ( d_kl, d_jk, lk, ll, ldim )
        call grad_transp ( p_kl, p_jk, lk, ll, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( d_il, d_jk, li, ll, ldim )
        call grad_copy ( p_il, p_jk, li, ll, ldim )
      else if ( ejl.and.eik ) then
        call grad_transp ( d_il, d_jk, li, ll, ldim )
        call grad_transp ( p_il, p_jk, li, ll, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab3f, iab3l, d_jk, ldim)
        call ga_get (g_dens(3), iab2f, iab2l, iab3f, iab3l, p_jk, ldim)
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_ump2 ( 
     $           da_ij, da_kl, da_ik, da_jl, da_il, da_jk, 
     $           db_ij, db_kl, db_ik, db_jl, db_il, db_jk, 
     $           pa_ij, pa_kl, pa_ik, pa_jl, pa_il, pa_jk, 
     $           pb_ij, pb_kl, pb_ik, pb_jl, pb_il, pb_jk, 
     $           ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )
C
C     get atomic blocks from density matrices and make sure no block is
C     transferred twice
C
C     (a,b) denotes spin, d is the HF part and p contains (twice) the 
C     second order correction
C
      implicit none

#include "global.fh"

      integer g_dens, ldim, 
     $     iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $     li, lj, lk, ll

      double precision da_ij, da_kl, da_ik, da_jl, da_il, da_jk,
     $                 db_ij, db_kl, db_ik, db_jl, db_il, db_jk,
     $                 pa_ij, pa_kl, pa_ik, pa_jl, pa_il, pa_jk,
     $                 pb_ij, pb_kl, pb_ik, pb_jl, pb_il, pb_jk

      dimension g_dens(5)
      dimension da_ij(ldim,ldim),da_kl(ldim,ldim),da_ik(ldim,ldim),
     $          da_jl(ldim,ldim),da_il(ldim,ldim),da_jk(ldim,ldim),
     $          db_ij(ldim,ldim),db_kl(ldim,ldim),db_ik(ldim,ldim),
     $          db_jl(ldim,ldim),db_il(ldim,ldim),db_jk(ldim,ldim),
     $          pa_ij(ldim,ldim),pa_kl(ldim,ldim),pa_ik(ldim,ldim),
     $          pa_jl(ldim,ldim),pa_il(ldim,ldim),pa_jk(ldim,ldim),
     $          pb_ij(ldim,ldim),pb_kl(ldim,ldim),pb_ik(ldim,ldim),
     $          pb_jl(ldim,ldim),pb_il(ldim,ldim),pb_jk(ldim,ldim)

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
      call ga_get (g_dens(2), iab1f, iab1l, iab2f, iab2l, da_ij, ldim)
      call ga_get (g_dens(3), iab1f, iab1l, iab2f, iab2l, db_ij, ldim)
      call ga_get (g_dens(4), iab1f, iab1l, iab2f, iab2l, pa_ij, ldim)
      call ga_get (g_dens(5), iab1f, iab1l, iab2f, iab2l, pb_ij, ldim)

C     k,l
      if ( eik.and.ejl ) then
        call grad_copy ( da_ij, da_kl, li, lj, ldim )
        call grad_copy ( db_ij, db_kl, li, lj, ldim )
        call grad_copy ( pa_ij, pa_kl, li, lj, ldim )
        call grad_copy ( pb_ij, pb_kl, li, lj, ldim )
      else if ( eil.and.ejk ) then
        call grad_transp ( da_ij, da_kl, li, lj, ldim )
        call grad_transp ( db_ij, db_kl, li, lj, ldim )
        call grad_transp ( pa_ij, pa_kl, li, lj, ldim )
        call grad_transp ( pb_ij, pb_kl, li, lj, ldim )
      else
        call ga_get (g_dens(2), iab3f, iab3l, iab4f, iab4l, da_kl, ldim)
        call ga_get (g_dens(3), iab3f, iab3l, iab4f, iab4l, db_kl, ldim)
        call ga_get (g_dens(4), iab3f, iab3l, iab4f, iab4l, pa_kl, ldim)
        call ga_get (g_dens(5), iab3f, iab3l, iab4f, iab4l, pb_kl, ldim)
      endif

C     i,k
      if ( ejk ) then
        call grad_copy ( da_ij, da_ik, li, lj, ldim )
        call grad_copy ( db_ij, db_ik, li, lj, ldim )
        call grad_copy ( pa_ij, pa_ik, li, lj, ldim )
        call grad_copy ( pb_ij, pb_ik, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( da_kl, da_ik, lk, ll, ldim )
        call grad_transp ( db_kl, db_ik, lk, ll, ldim )
        call grad_transp ( pa_kl, pa_ik, lk, ll, ldim )
        call grad_transp ( pb_kl, pb_ik, lk, ll, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab3f, iab3l, da_ik, ldim)
        call ga_get (g_dens(3), iab1f, iab1l, iab3f, iab3l, db_ik, ldim)
        call ga_get (g_dens(4), iab1f, iab1l, iab3f, iab3l, pa_ik, ldim)
        call ga_get (g_dens(5), iab1f, iab1l, iab3f, iab3l, pb_ik, ldim)
      endif

C     j,l
      if ( ejk ) then
        call grad_copy ( da_kl, da_jl, lk, ll, ldim )
        call grad_copy ( db_kl, db_jl, lk, ll, ldim )
        call grad_copy ( pa_kl, pa_jl, lk, ll, ldim )
        call grad_copy ( pb_kl, pb_jl, lk, ll, ldim )
      else if ( eij.and.ejl ) then
        call grad_copy ( da_ij, da_jl, li, lj, ldim )
        call grad_copy ( db_ij, db_jl, li, lj, ldim )
        call grad_copy ( pa_ij, pa_jl, li, lj, ldim )
        call grad_copy ( pb_ij, pb_jl, li, lj, ldim )
      else if ( eil ) then
        call grad_transp ( da_ij, da_jl, li, lj, ldim )
        call grad_transp ( db_ij, db_jl, li, lj, ldim )
        call grad_transp ( pa_ij, pa_jl, li, lj, ldim )
        call grad_transp ( pb_ij, pb_jl, li, lj, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( da_ik, da_jl, li, lk, ldim )
        call grad_copy ( db_ik, db_jl, li, lk, ldim )
        call grad_copy ( pa_ik, pa_jl, li, lk, ldim )
        call grad_copy ( pb_ik, pb_jl, li, lk, ldim )
      else if ( ejk.and.eil ) then
        call grad_transp ( da_ik, da_jl, li, lk, ldim )
        call grad_transp ( db_ik, db_jl, li, lk, ldim )
        call grad_transp ( pa_ik, pa_jl, li, lk, ldim )
        call grad_transp ( pb_ik, pb_jl, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab4f, iab4l, da_jl, ldim)
        call ga_get (g_dens(3), iab2f, iab2l, iab4f, iab4l, db_jl, ldim)
        call ga_get (g_dens(4), iab2f, iab2l, iab4f, iab4l, pa_jl, ldim)
        call ga_get (g_dens(5), iab2f, iab2l, iab4f, iab4l, pb_jl, ldim)
      endif

C     i,l
      if ( ejl ) then
        call grad_copy ( da_ij, da_il, li, lj, ldim )
        call grad_copy ( db_ij, db_il, li, lj, ldim )
        call grad_copy ( pa_ij, pa_il, li, lj, ldim )
        call grad_copy ( pb_ij, pb_il, li, lj, ldim )
      else if ( eik ) then
        call grad_copy ( da_kl, da_il, lk, ll, ldim )
        call grad_copy ( db_kl, db_il, lk, ll, ldim )
        call grad_copy ( pa_kl, pa_il, lk, ll, ldim )
        call grad_copy ( pb_kl, pb_il, lk, ll, ldim )
      else if ( eij ) then
        call grad_copy ( da_jl, da_il, lj, ll, ldim )
        call grad_copy ( db_jl, db_il, lj, ll, ldim )
        call grad_copy ( pa_jl, pa_il, lj, ll, ldim )
        call grad_copy ( pb_jl, pb_il, lj, ll, ldim )
      else if ( ekl ) then
        call grad_copy ( da_ik, da_il, li, lk, ldim )
        call grad_copy ( db_ik, db_il, li, lk, ldim )
        call grad_copy ( pa_ik, pa_il, li, lk, ldim )
        call grad_copy ( pb_ik, pb_il, li, lk, ldim )
      else
        call ga_get (g_dens(2), iab1f, iab1l, iab4f, iab4l, da_il, ldim)
        call ga_get (g_dens(3), iab1f, iab1l, iab4f, iab4l, db_il, ldim)
        call ga_get (g_dens(4), iab1f, iab1l, iab4f, iab4l, pa_il, ldim)
        call ga_get (g_dens(5), iab1f, iab1l, iab4f, iab4l, pb_il, ldim)
      endif

C     j,k
      if ( ekl ) then
        call grad_copy ( da_jl, da_jk, lj, ll, ldim )
        call grad_copy ( db_jl, db_jk, lj, ll, ldim )
        call grad_copy ( pa_jl, pa_jk, lj, ll, ldim )
        call grad_copy ( pb_jl, pb_jk, lj, ll, ldim )
      else if ( eij ) then
        call grad_copy ( da_ik, da_jk, li, lk, ldim )
        call grad_copy ( db_ik, db_jk, li, lk, ldim )
        call grad_copy ( pa_ik, pa_jk, li, lk, ldim )
        call grad_copy ( pb_ik, pb_jk, li, lk, ldim )
      else if ( eik ) then
        call grad_transp ( da_ij, da_jk, li, lj, ldim )
        call grad_transp ( db_ij, db_jk, li, lj, ldim )
        call grad_transp ( pa_ij, pa_jk, li, lj, ldim )
        call grad_transp ( pb_ij, pb_jk, li, lj, ldim )
      else if ( ejl ) then
        call grad_transp ( da_kl, da_jk, lk, ll, ldim )
        call grad_transp ( db_kl, db_jk, lk, ll, ldim )
        call grad_transp ( pa_kl, pa_jk, lk, ll, ldim )
        call grad_transp ( pb_kl, pb_jk, lk, ll, ldim )
      else if ( eij.and.ekl ) then
        call grad_copy ( da_il, da_jk, li, ll, ldim )
        call grad_copy ( db_il, db_jk, li, ll, ldim )
        call grad_copy ( pa_il, pa_jk, li, ll, ldim )
        call grad_copy ( pb_il, pb_jk, li, ll, ldim )
      else if ( ejl.and.eik ) then
        call grad_transp ( da_il, da_jk, li, ll, ldim )
        call grad_transp ( db_il, db_jk, li, ll, ldim )
        call grad_transp ( pa_il, pa_jk, li, ll, ldim )
        call grad_transp ( pb_il, pb_jk, li, ll, ldim )
      else
        call ga_get (g_dens(2), iab2f, iab2l, iab3f, iab3l, da_jk, ldim)
        call ga_get (g_dens(3), iab2f, iab2l, iab3f, iab3l, db_jk, ldim)
        call ga_get (g_dens(4), iab2f, iab2l, iab3f, iab3l, pa_jk, ldim)
        call ga_get (g_dens(5), iab2f, iab2l, iab3f, iab3l, pb_jk, ldim)
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
