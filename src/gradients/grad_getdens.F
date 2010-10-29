      subroutine grad_get_dens_all (omp2,ouhf,orohf,
     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk,
     $     d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $     d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $     d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4, 
     $     d_ik5, d_jl5, d_il5, d_jk5,
     $     d_ik6, d_jl6, d_il6, d_jk6, 
     $     d_ik7, d_jl7, d_il7, d_jk7, 
     $     d_ik8, d_jl8, d_il8, d_jk8, 
     $     blen,
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     ilo_prev, jlo_prev, klo_prev,
     $     g_dens,exso)
c
c     Wrapper for all of the grad_get_dens routines.  Note that
c     some of the arrays passed into here are not allocated
c     if they are not being used at the given level of theory.
c     The addresses will be valid but DON'T write to them.
c
      implicit none
#include "errquit.fh"
#include "global.fh"
      logical omp2,ouhf,orohf
      integer g_dens(*), blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     ilo_prev, jlo_prev, klo_prev
      double precision 
     $     d_ij(*),  d_kl(*),  d_ik(*),  d_jl(*),  d_il(*),  d_jk(*),
     $     d_ij2(*), d_kl2(*), d_ik2(*), d_jl2(*), d_il2(*), d_jk2(*), 
     $     d_ij3(*), d_kl3(*), d_ik3(*), d_jl3(*), d_il3(*), d_jk3(*), 
     $     d_ij4(*), d_kl4(*), d_ik4(*), d_jl4(*), d_il4(*), d_jk4(*),
     $     d_ik5(*), d_jl5(*), d_il5(*), d_jk5(*),
     $     d_ik6(*), d_jl6(*), d_il6(*), d_jk6(*), 
     $     d_ik7(*), d_jl7(*), d_il7(*), d_jk7(*), 
     $     d_ik8(*), d_jl8(*), d_il8(*), d_jk8(*)
c
      logical ogetij, ogetik, ogetjk, exso
c
c     Reuse data we already have ... l is inner loop so cannot save this.
c
      ogetij = ilo.ne.ilo_prev .or. jlo.ne.jlo_prev
      ogetik = ilo.ne.ilo_prev .or. klo.ne.klo_prev
      ogetjk = jlo.ne.jlo_prev .or. klo.ne.klo_prev
c
      ilo_prev = ilo
      jlo_prev = jlo
      klo_prev = klo

c
      if (.not. omp2) then
         if (ouhf) then
            if(.not.exso)
     $           call grad_get_dens_uhf(d_ij, d_kl, d_ik, 
     $           d_jl, d_il, d_jk,
     $           d_ik2, d_jl2, d_il2, d_jk2, blen, 
     $           ilo, ihi, jlo, jhi, 
     $           klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
            if(exso)
     $           call grad_get_dens_exso(
     &           d_ij, d_kl, 
     &           d_ik,  d_jl,  d_il,  d_jk,
     $           d_ik2, d_jl2, d_il2, d_jk2, 
     &           d_ik3, d_jl3, d_il3, d_jk3, 
     &           d_ik4, d_jl4, d_il4, d_jk4,
     &           d_ik5, d_jl5, d_il5, d_jk5, 
     &           d_ik6, d_jl6, d_il6, d_jk6,
     &           d_ik7, d_jl7, d_il7, d_jk7, 
     &           d_ik8, d_jl8, d_il8, d_jk8,
     &           blen, 
     $           ilo, ihi, jlo, jhi, 
     $           klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
         else if (orohf) then
            call grad_get_dens_rohf(d_ij, d_kl, d_ik, 
     $           d_jl, d_il, d_jk,
     $           d_ij2, d_kl2, d_ik2, 
     $           d_jl2, d_il2, d_jk2, blen, 
     $           ilo, ihi, jlo, jhi, 
     $           klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
         else
            call grad_get_dens(d_ij, d_kl, d_ik, 
     $           d_jl, d_il, d_jk, blen,
     $           ilo, ihi, jlo, jhi, 
     $           klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens(1))
         end if
      else                      ! MP2
         if (ouhf) then
            call grad_get_dens_ump2(d_ij, d_kl, d_ik, 
     $           d_jl, d_il, d_jk,
     $           d_ij2, d_kl2, d_ik2, 
     $           d_jl2, d_il2, d_jk2, 
     $           d_ij3, d_kl3, d_ik3, 
     $           d_jl3, d_il3, d_jk3, 
     $           d_ij4, d_kl4, d_ik4, 
     $           d_jl4, d_il4, d_jk4, 
     $           blen,
     $           ilo, ihi, jlo, jhi, 
     $           klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
         else if (orohf) then
            call errquit('no ROMP2 gradients yet', 0, CAPMIS_ERR)
         else                   ! rmp2
            call grad_get_dens_rmp2(d_ij, d_kl, d_ik, 
     $           d_jl, d_il, d_jk,
     $           d_ij2, d_kl2, d_ik2, 
     $           d_jl2, d_il2, d_jk2, 
     $           blen,
     $           ilo, ihi, jlo, jhi, 
     $           klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
         end if
      end if
c
      end
      subroutine grad_get_dens(d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           blen,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
c$Id$
C
C     get blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     li, lj, lk, ll
      logical ogetij, ogetik, ogetjk

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk

      dimension d_ij(blen,blen),d_kl(blen,blen),d_ik(blen,blen),
     $          d_jl(blen,blen),d_il(blen,blen),d_jk(blen,blen)

      logical eij, ejk, ekl, eil, eik, ejl

      eij = (ilo.eq.jlo)
      ejk = (jlo.eq.klo)
      ekl = (klo.eq.llo)
      eil = (ilo.eq.llo)
      eik = (ilo.eq.klo)
      ejl = (jlo.eq.llo)

      li = ihi - ilo + 1
      lj = jhi - jlo + 1
      lk = khi - klo + 1
      ll = lhi - llo + 1

C     i,j
      if (ogetij) call ga_get (g_dens, ilo, ihi, jlo, jhi, d_ij, blen)

C     k,l
      if (eik.and.ejl) then
        call grad_copy (d_ij, d_kl, li, lj, blen)
      else if (eil.and.ejk) then
        call grad_transp (d_ij, d_kl, li, lj, blen)
      else
        call ga_get (g_dens, klo, khi, llo, lhi, d_kl, blen)
      endif

C     i,k
      if (ogetik) then
         if (ejk) then
            call grad_copy (d_ij, d_ik, li, lj, blen)
         else if (eil) then
            call grad_transp (d_kl, d_ik, lk, ll, blen)
         else
            call ga_get (g_dens, ilo, ihi, klo, khi, d_ik, blen)
         endif
      endif

C     j,l
      if (ejk) then
        call grad_copy (d_kl, d_jl, lk, ll, blen)
      else if (eij.and.ejl) then
        call grad_copy (d_ij, d_jl, li, lj, blen)
      else if (eil) then
        call grad_transp (d_ij, d_jl, li, lj, blen)
      else if (eij.and.ekl) then
        call grad_copy (d_ik, d_jl, li, lk, blen)
      else if (ejk.and.eil) then
        call grad_transp (d_ik, d_jl, li, lk, blen)
      else
        call ga_get (g_dens, jlo, jhi, llo, lhi, d_jl, blen)
      endif

C     i,l
      if (ejl) then
        call grad_copy (d_ij, d_il, li, lj, blen)
      else if (eik) then
        call grad_copy (d_kl, d_il, lk, ll, blen)
      else if (eij) then
        call grad_copy (d_jl, d_il, lj, ll, blen)
      else if (ekl) then
        call grad_copy (d_ik, d_il, li, lk, blen)
      else
        call ga_get (g_dens, ilo, ihi, llo, lhi, d_il, blen)
      endif

C     j,k
      if (ogetjk) then
         if (ekl) then
            call grad_copy (d_jl, d_jk, lj, ll, blen)
         else if (eij) then
            call grad_copy (d_ik, d_jk, li, lk, blen)
         else if (eik) then
            call grad_transp (d_ij, d_jk, li, lj, blen)
         else if (ejl) then
            call grad_transp (d_kl, d_jk, lk, ll, blen)
         else if (eij.and.ekl) then
            call grad_copy (d_il, d_jk, li, ll, blen)
         else if (ejl.and.eik) then
            call grad_transp (d_il, d_jk, li, ll, blen)
         else
            call ga_get (g_dens, jlo, jhi, klo, khi, d_jk, blen)
         endif
      endif
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_rohf(d_ij, d_kl, d_ik, d_jl, d_il, d_jk,  
     $           d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, blen,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
C
C     get blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     li, lj, lk, ll
      logical ogetij, ogetik, ogetjk

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $          d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2

      dimension d_ij(blen,blen),d_kl(blen,blen),d_ik(blen,blen),
     $          d_jl(blen,blen),d_il(blen,blen),d_jk(blen,blen),
     $          d_ij2(blen,blen),d_kl2(blen,blen),d_ik2(blen,blen),
     $          d_jl2(blen,blen),d_il2(blen,blen),d_jk2(blen,blen)
      dimension g_dens(3)

      logical eij, ejk, ekl, eil, eik, ejl

      eij = (ilo.eq.jlo)
      ejk = (jlo.eq.klo)
      ekl = (klo.eq.llo)
      eil = (ilo.eq.llo)
      eik = (ilo.eq.klo)
      ejl = (jlo.eq.llo)

      li = ihi - ilo + 1
      lj = jhi - jlo + 1
      lk = khi - klo + 1
      ll = lhi - llo + 1

C     i,j
      if (ogetij) then
         call ga_get (g_dens(2), ilo, ihi, jlo, jhi, d_ij, blen)
         call ga_get (g_dens(3), ilo, ihi, jlo, jhi, d_ij2, blen)
      endif

C     k,l
      if (eik.and.ejl) then
        call grad_copy (d_ij, d_kl, li, lj, blen)
        call grad_copy (d_ij2, d_kl2, li, lj, blen)
      else if (eil.and.ejk) then
        call grad_transp (d_ij, d_kl, li, lj, blen)
        call grad_transp (d_ij2, d_kl2, li, lj, blen)
      else
        call ga_get (g_dens(2), klo, khi, llo, lhi, d_kl, blen)
        call ga_get (g_dens(3), klo, khi, llo, lhi, d_kl2, blen)
      endif

C     i,k
      if (ogetik) then
         if (ejk) then
            call grad_copy (d_ij, d_ik, li, lj, blen)
            call grad_copy (d_ij2, d_ik2, li, lj, blen)
         else if (eil) then
            call grad_transp (d_kl, d_ik, lk, ll, blen)
            call grad_transp (d_kl2, d_ik2, lk, ll, blen)
         else 
            call ga_get (g_dens(2), ilo, ihi, klo, khi, d_ik, blen)
            call ga_get (g_dens(3), ilo, ihi, klo, khi, d_ik2, blen)
         endif
      endif

C     j,l
      if (ejk) then
        call grad_copy (d_kl, d_jl, lk, ll, blen)
        call grad_copy (d_kl2, d_jl2, lk, ll, blen)
      else if (eij.and.ejl) then
        call grad_copy (d_ij, d_jl, li, lj, blen)
        call grad_copy (d_ij2, d_jl2, li, lj, blen)
      else if (eil) then
        call grad_transp (d_ij, d_jl, li, lj, blen)
        call grad_transp (d_ij2, d_jl2, li, lj, blen)
      else if (eij.and.ekl) then
        call grad_copy (d_ik, d_jl, li, lk, blen)
        call grad_copy (d_ik2, d_jl2, li, lk, blen)
      else if (ejk.and.eil) then
        call grad_transp (d_ik, d_jl, li, lk, blen)
        call grad_transp (d_ik2, d_jl2, li, lk, blen)
      else
        call ga_get (g_dens(2), jlo, jhi, llo, lhi, d_jl, blen)
        call ga_get (g_dens(3), jlo, jhi, llo, lhi, d_jl2, blen)
      endif

C     i,l
      if (ejl) then
        call grad_copy (d_ij, d_il, li, lj, blen)
        call grad_copy (d_ij2, d_il2, li, lj, blen)
      else if (eik) then
        call grad_copy (d_kl, d_il, lk, ll, blen)
        call grad_copy (d_kl2, d_il2, lk, ll, blen)
      else if (eij) then
        call grad_copy (d_jl, d_il, lj, ll, blen)
        call grad_copy (d_jl2, d_il2, lj, ll, blen)
      else if (ekl) then
        call grad_copy (d_ik, d_il, li, lk, blen)
        call grad_copy (d_ik2, d_il2, li, lk, blen)
      else
        call ga_get (g_dens(2), ilo, ihi, llo, lhi, d_il, blen)
        call ga_get (g_dens(3), ilo, ihi, llo, lhi, d_il2, blen)
      endif

C     j,k
      if (ogetjk) then
         if (ekl) then
            call grad_copy (d_jl, d_jk, lj, ll, blen)
            call grad_copy (d_jl2, d_jk2, lj, ll, blen)
         else if (eij) then
            call grad_copy (d_ik, d_jk, li, lk, blen)
            call grad_copy (d_ik2, d_jk2, li, lk, blen)
         else if (eik) then
            call grad_transp (d_ij, d_jk, li, lj, blen)
            call grad_transp (d_ij2, d_jk2, li, lj, blen)
         else if (ejl) then
            call grad_transp (d_kl, d_jk, lk, ll, blen)
            call grad_transp (d_kl2, d_jk2, lk, ll, blen)
         else if (eij.and.ekl) then
            call grad_copy (d_il, d_jk, li, ll, blen)
            call grad_copy (d_il2, d_jk2, li, ll, blen)
         else if (ejl.and.eik) then
            call grad_transp (d_il, d_jk, li, ll, blen)
            call grad_transp (d_il2, d_jk2, li, ll, blen)
         else
            call ga_get (g_dens(2), jlo, jhi, klo, khi, d_jk, blen)
            call ga_get (g_dens(3), jlo, jhi, klo, khi, d_jk2, blen)
         endif
      endif
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_uhf(d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           d_ik2, d_jl2, d_il2, d_jk2, blen,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
C
C     get blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     li, lj, lk, ll
      logical ogetij, ogetik, ogetjk

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     d_ik2, d_jl2, d_il2, d_jk2

      dimension d_ij(blen,blen),d_kl(blen,blen),d_ik(blen,blen),
     $          d_jl(blen,blen),d_il(blen,blen),d_jk(blen,blen),
     $          d_ik2(blen,blen),d_jl2(blen,blen),d_il2(blen,blen),
     $          d_jk2(blen,blen)
      dimension g_dens(2)

      logical eij, ejk, ekl, eil, eik, ejl

      eij = (ilo.eq.jlo)
      ejk = (jlo.eq.klo)
      ekl = (klo.eq.llo)
      eil = (ilo.eq.llo)
      eik = (ilo.eq.klo)
      ejl = (jlo.eq.llo)

      li = ihi - ilo + 1
      lj = jhi - jlo + 1
      lk = khi - klo + 1
      ll = lhi - llo + 1

C     i,j
      if (ogetij) call ga_get(g_dens(1), ilo, ihi, jlo, jhi, d_ij, blen)

C     k,l
      if (eik.and.ejl) then
        call grad_copy (d_ij, d_kl, li, lj, blen)
      else if (eil.and.ejk) then
        call grad_transp (d_ij, d_kl, li, lj, blen)
      else
        call ga_get (g_dens(1), klo, khi, llo, lhi, d_kl, blen)
      endif

C     i,k
      if (ogetik) then
         if (ejk) then
            call grad_copy (d_ij, d_ik, li, lj, blen)
         else if (eil) then
            call grad_transp (d_kl, d_ik, lk, ll, blen)
         else
            call ga_get (g_dens(1), ilo, ihi, klo, khi, d_ik, blen)
         endif
         call ga_get (g_dens(2), ilo, ihi, klo, khi, d_ik2, blen)
      endif

C     j,l
      if (ejk) then
        call grad_copy (d_kl, d_jl, lk, ll, blen)
      else if (eij.and.ejl) then
        call grad_copy (d_ij, d_jl, li, lj, blen)
      else if (eil) then
        call grad_transp (d_ij, d_jl, li, lj, blen)
      else if (eij.and.ekl) then
        call grad_copy (d_ik, d_jl, li, lk, blen)
      else if (ejk.and.eil) then
        call grad_transp (d_ik, d_jl, li, lk, blen)
      else
        call ga_get (g_dens(1), jlo, jhi, llo, lhi, d_jl, blen)
      endif

      if (eij .and. ekl) then
        call grad_copy (d_ik2, d_jl2, li, lk, blen)
      else if (eil .and. ejk) then
        call grad_transp (d_ik2, d_jl2, li, lk, blen)
      else
        call ga_get (g_dens(2), jlo, jhi, llo, lhi, d_jl2, blen)
      endif

C     i,l
      if (ejl) then
        call grad_copy (d_ij, d_il, li, lj, blen)
      else if (eik) then
        call grad_copy (d_kl, d_il, lk, ll, blen)
      else if (eij) then
        call grad_copy (d_jl, d_il, lj, ll, blen)
      else if (ekl) then
        call grad_copy (d_ik, d_il, li, lk, blen)
      else
        call ga_get (g_dens(1), ilo, ihi, llo, lhi, d_il, blen)
      endif

      if (ekl) then
        call grad_copy (d_ik2, d_il2, li, lk, blen)
      else if (eij) then
        call grad_copy (d_jl2, d_il2, lj, ll, blen)
      else
        call ga_get (g_dens(2), ilo, ihi, llo, lhi, d_il2, blen)
      endif

C     j,k
      if (ogetjk) then
         if (ekl) then
            call grad_copy (d_jl, d_jk, lj, ll, blen)
         else if (eij) then
            call grad_copy (d_ik, d_jk, li, lk, blen)
         else if (eik) then
            call grad_transp (d_ij, d_jk, li, lj, blen)
         else if (ejl) then
            call grad_transp (d_kl, d_jk, lk, ll, blen)
         else if (eij.and.ekl) then
            call grad_copy (d_il, d_jk, li, ll, blen)
         else if (ejl.and.eik) then
            call grad_transp (d_il, d_jk, li, ll, blen)
         else
            call ga_get (g_dens(1), jlo, jhi, klo, khi, d_jk, blen)
         endif

         if (eij) then
            call grad_copy (d_ik2, d_jk2, li, lk, blen)
         else if (ekl) then
            call grad_copy (d_jl2, d_jk2, lj, ll, blen)
         else if (eik .and. ejl) then
            call grad_transp (d_il2, d_jk2, li, ll, blen)
         else
            call ga_get (g_dens(2), jlo, jhi, klo, khi, d_jk2, blen)
         endif
      endif
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_exso(d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           d_ik2, d_jl2, d_il2, d_jk2, 
     $           d_ik3, d_jl3, d_il3, d_jk3, 
     $           d_ik4, d_jl4, d_il4, d_jk4, 
     $           d_ik5, d_jl5, d_il5, d_jk5, 
     $           d_ik6, d_jl6, d_il6, d_jk6, 
     $           d_ik7, d_jl7, d_il7, d_jk7, 
     $           d_ik8, d_jl8, d_il8, d_jk8, 
     &           blen,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
C
C     get blocks from density matrix and make sure no block is
C     transferred twice
C
      implicit none

#include "global.fh"

      integer g_dens, blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     li, lj, lk, ll
      logical ogetij, ogetik, ogetjk

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     d_ik2, d_jl2, d_il2, d_jk2,
     $     d_ik3, d_jl3, d_il3, d_jk3,
     $     d_ik4, d_jl4, d_il4, d_jk4,
     $     d_ik5, d_jl5, d_il5, d_jk5,
     $     d_ik6, d_jl6, d_il6, d_jk6,
     $     d_ik7, d_jl7, d_il7, d_jk7,
     $     d_ik8, d_jl8, d_il8, d_jk8

      dimension d_ij(blen,blen),d_kl(blen,blen),d_ik(blen,blen),
     $          d_jl(blen,blen),d_il(blen,blen),d_jk(blen,blen),
     $          d_ik2(blen,blen),d_jl2(blen,blen),d_il2(blen,blen),
     $          d_jk2(blen,blen),
     $          d_ik3(blen,blen),d_jl3(blen,blen),d_il3(blen,blen),
     $          d_jk3(blen,blen),
     $          d_ik4(blen,blen),d_jl4(blen,blen),d_il4(blen,blen),
     $          d_jk4(blen,blen),
     $          d_ik5(blen,blen),d_jl5(blen,blen),d_il5(blen,blen),
     $          d_jk5(blen,blen),
     $          d_ik6(blen,blen),d_jl6(blen,blen),d_il6(blen,blen),
     $          d_jk6(blen,blen),
     $          d_ik7(blen,blen),d_jl7(blen,blen),d_il7(blen,blen),
     $          d_jk7(blen,blen),
     $          d_ik8(blen,blen),d_jl8(blen,blen),d_il8(blen,blen),
     $          d_jk8(blen,blen)
      dimension g_dens(8)

      logical eij, ejk, ekl, eil, eik, ejl

      eij = (ilo.eq.jlo)
      ejk = (jlo.eq.klo)
      ekl = (klo.eq.llo)
      eil = (ilo.eq.llo)
      eik = (ilo.eq.klo)
      ejl = (jlo.eq.llo)

      li = ihi - ilo + 1
      lj = jhi - jlo + 1
      lk = khi - klo + 1
      ll = lhi - llo + 1

C     i,j
      if (ogetij)call ga_get(g_dens(1), ilo, ihi, jlo, jhi, d_ij, blen)

C     k,l
      if (eik.and.ejl) then
        call grad_copy (d_ij, d_kl, li, lj, blen)
      else if (eil.and.ejk) then
        call grad_transp (d_ij, d_kl, li, lj, blen)
      else
        call ga_get (g_dens(1), klo, khi, llo, lhi, d_kl, blen)
      endif

C     i,k
      if (ogetik) then
c     dens(1) and dens(2)
         if (ejk) then
            call grad_copy (d_ij, d_ik, li, lj, blen)
         else if (eil) then
            call grad_transp (d_kl, d_ik, lk, ll, blen)
         else
            call ga_get (g_dens(1), ilo, ihi, klo, khi, d_ik, blen)
         endif
         call ga_get (g_dens(2), ilo, ihi, klo, khi, d_ik2, blen)
c     dens(3) and dens(4)
         call ga_get (g_dens(3), ilo, ihi, klo, khi, d_ik3, blen)
         call ga_get (g_dens(4), ilo, ihi, klo, khi, d_ik4, blen)
c     dens(5) and dens(6)
         call ga_get (g_dens(5), ilo, ihi, klo, khi, d_ik5, blen)
         call ga_get (g_dens(6), ilo, ihi, klo, khi, d_ik6, blen)
c     dens(7) and dens(8)
         call ga_get (g_dens(7), ilo, ihi, klo, khi, d_ik7, blen)
         call ga_get (g_dens(8), ilo, ihi, klo, khi, d_ik8, blen)
      endif

C     j,l
c     dens(1) and dens(2)
      if (ejk) then
         call grad_copy (d_kl, d_jl, lk, ll, blen)
      else if (eij.and.ejl) then
         call grad_copy (d_ij, d_jl, li, lj, blen)
      else if (eil) then
         call grad_transp (d_ij, d_jl, li, lj, blen)
      else if (eij.and.ekl) then
         call grad_copy (d_ik, d_jl, li, lk, blen)
      else if (ejk.and.eil) then
         call grad_transp (d_ik, d_jl, li, lk, blen)
      else
         call ga_get (g_dens(1), jlo, jhi, llo, lhi, d_jl, blen)
      endif

      if (eij .and. ekl) then
        call grad_copy (d_ik2, d_jl2, li, lk, blen)
      else if (eil .and. ejk) then
        call grad_transp (d_ik2, d_jl2, li, lk, blen)
      else
        call ga_get (g_dens(2), jlo, jhi, llo, lhi, d_jl2, blen)
      endif
c
c     dens(3) and dens(4)
c
      call ga_get (g_dens(3), jlo, jhi, llo, lhi, d_jl3, blen)
      call ga_get (g_dens(4), jlo, jhi, llo, lhi, d_jl4, blen)
c
c     dens(5) and dens(6)
c
      call ga_get (g_dens(5), jlo, jhi, llo, lhi, d_jl5, blen)
      call ga_get (g_dens(6), jlo, jhi, llo, lhi, d_jl6, blen)
c
c     dens(7) and dens(8)
c
      call ga_get (g_dens(7), jlo, jhi, llo, lhi, d_jl7, blen)
      call ga_get (g_dens(8), jlo, jhi, llo, lhi, d_jl8, blen)


C     i,l
c
c     dens(1) and dens(2)
c
      if (ejl) then
        call grad_copy (d_ij, d_il, li, lj, blen)
      else if (eik) then
        call grad_copy (d_kl, d_il, lk, ll, blen)
      else if (eij) then
        call grad_copy (d_jl, d_il, lj, ll, blen)
      else if (ekl) then
        call grad_copy (d_ik, d_il, li, lk, blen)
      else
        call ga_get (g_dens(1), ilo, ihi, llo, lhi, d_il, blen)
      endif

      if (ekl) then
        call grad_copy (d_ik2, d_il2, li, lk, blen)
      else if (eij) then
        call grad_copy (d_jl2, d_il2, lj, ll, blen)
      else
        call ga_get (g_dens(2), ilo, ihi, llo, lhi, d_il2, blen)
      endif
c
c     dens(3) and dnes(4)
c
      call ga_get (g_dens(3), ilo, ihi, llo, lhi, d_il3, blen)
      call ga_get (g_dens(4), ilo, ihi, llo, lhi, d_il4, blen)
c
c     dens(5) and dnes(6)
c     
      call ga_get (g_dens(5), ilo, ihi, llo, lhi, d_il5, blen)
      call ga_get (g_dens(6), ilo, ihi, llo, lhi, d_il6, blen)
c     
c     dens(7) and dnes(8)
c     
      call ga_get (g_dens(7), ilo, ihi, llo, lhi, d_il7, blen)
      call ga_get (g_dens(8), ilo, ihi, llo, lhi, d_il8, blen)

C     j,k
      if (ogetjk) then
c
c     dens(1) and dens(2)
c
         if (ekl) then
            call grad_copy (d_jl, d_jk, lj, ll, blen)
         else if (eij) then
            call grad_copy (d_ik, d_jk, li, lk, blen)
         else if (eik) then
            call grad_transp (d_ij, d_jk, li, lj, blen)
         else if (ejl) then
            call grad_transp (d_kl, d_jk, lk, ll, blen)
         else if (eij.and.ekl) then
            call grad_copy (d_il, d_jk, li, ll, blen)
         else if (ejl.and.eik) then
            call grad_transp (d_il, d_jk, li, ll, blen)
         else
            call ga_get (g_dens(1), jlo, jhi, klo, khi, d_jk, blen)
         endif
c
         if (eij) then
            call grad_copy (d_ik2, d_jk2, li, lk, blen)
         else if (ekl) then
            call grad_copy (d_jl2, d_jk2, lj, ll, blen)
         else if (eik .and. ejl) then
            call grad_transp (d_il2, d_jk2, li, ll, blen)
         else
            call ga_get (g_dens(2), jlo, jhi, klo, khi, d_jk2, blen)
         endif
c
c     dens(3) and dens(4)
c
            call ga_get (g_dens(3), jlo, jhi, klo, khi, d_jk3, blen)
            call ga_get (g_dens(4), jlo, jhi, klo, khi, d_jk4, blen)
c
c     dens(5) and dens(6)
c
            call ga_get (g_dens(5), jlo, jhi, klo, khi, d_jk5, blen)
            call ga_get (g_dens(6), jlo, jhi, klo, khi, d_jk6, blen)
c
c     dens(7) and dens(8)
c
            call ga_get (g_dens(7), jlo, jhi, klo, khi, d_jk7, blen)
            call ga_get (g_dens(8), jlo, jhi, klo, khi, d_jk8, blen)
      endif
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_rmp2 (
     $           d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $           p_ij, p_kl, p_ik, p_jl, p_il, p_jk, 
     $           blen,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
C
C     get blocks from density matrices and make sure no block is
C     transferred twice
C
C     (a,b) denotes spin, d is the HF part and p contains (twice) the 
C     second order correction
C
      implicit none

#include "global.fh"

      integer g_dens, blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     li, lj, lk, ll

      double precision d_ij, d_kl, d_ik, d_jl, d_il, d_jk,
     $                 p_ij, p_kl, p_ik, p_jl, p_il, p_jk
      logical ogetij, ogetik, ogetjk

      dimension g_dens(3)
      dimension d_ij(blen,blen),d_kl(blen,blen),d_ik(blen,blen),
     $          d_jl(blen,blen),d_il(blen,blen),d_jk(blen,blen),
     $          p_ij(blen,blen),p_kl(blen,blen),p_ik(blen,blen),
     $          p_jl(blen,blen),p_il(blen,blen),p_jk(blen,blen)

      logical eij, ejk, ekl, eil, eik, ejl

      eij = (ilo.eq.jlo)
      ejk = (jlo.eq.klo)
      ekl = (klo.eq.llo)
      eil = (ilo.eq.llo)
      eik = (ilo.eq.klo)
      ejl = (jlo.eq.llo)

      li = ihi - ilo + 1
      lj = jhi - jlo + 1
      lk = khi - klo + 1
      ll = lhi - llo + 1

C     i,j
      if (ogetij) then
         call ga_get (g_dens(2), ilo, ihi, jlo, jhi, d_ij, blen)
         call ga_get (g_dens(3), ilo, ihi, jlo, jhi, p_ij, blen)
      endif

C     k,l
      if (eik.and.ejl) then
        call grad_copy (d_ij, d_kl, li, lj, blen)
        call grad_copy (p_ij, p_kl, li, lj, blen)
      else if (eil.and.ejk) then
        call grad_transp (d_ij, d_kl, li, lj, blen)
        call grad_transp (p_ij, p_kl, li, lj, blen)
      else
        call ga_get (g_dens(2), klo, khi, llo, lhi, d_kl, blen)
        call ga_get (g_dens(3), klo, khi, llo, lhi, p_kl, blen)
      endif

C     i,k
      if (ogetik) then
         if (ejk) then
            call grad_copy (d_ij, d_ik, li, lj, blen)
            call grad_copy (p_ij, p_ik, li, lj, blen)
         else if (eil) then
            call grad_transp (d_kl, d_ik, lk, ll, blen)
            call grad_transp (p_kl, p_ik, lk, ll, blen)
         else
            call ga_get (g_dens(2), ilo, ihi, klo, khi, d_ik, blen)
            call ga_get (g_dens(3), ilo, ihi, klo, khi, p_ik, blen)
         endif
      endif

C     j,l
      if (ejk) then
         call grad_copy (d_kl, d_jl, lk, ll, blen)
         call grad_copy (p_kl, p_jl, lk, ll, blen)
      else if (eij.and.ejl) then
         call grad_copy (d_ij, d_jl, li, lj, blen)
         call grad_copy (p_ij, p_jl, li, lj, blen)
      else if (eil) then
         call grad_transp (d_ij, d_jl, li, lj, blen)
         call grad_transp (p_ij, p_jl, li, lj, blen)
      else if (eij.and.ekl) then
         call grad_copy (d_ik, d_jl, li, lk, blen)
         call grad_copy (p_ik, p_jl, li, lk, blen)
      else if (ejk.and.eil) then
         call grad_transp (d_ik, d_jl, li, lk, blen)
         call grad_transp (p_ik, p_jl, li, lk, blen)
      else
         call ga_get (g_dens(2), jlo, jhi, llo, lhi, d_jl, blen)
         call ga_get (g_dens(3), jlo, jhi, llo, lhi, p_jl, blen)
      endif

C     i,l
      if (ejl) then
        call grad_copy (d_ij, d_il, li, lj, blen)
        call grad_copy (p_ij, p_il, li, lj, blen)
      else if (eik) then
        call grad_copy (d_kl, d_il, lk, ll, blen)
        call grad_copy (p_kl, p_il, lk, ll, blen)
      else if (eij) then
        call grad_copy (d_jl, d_il, lj, ll, blen)
        call grad_copy (p_jl, p_il, lj, ll, blen)
      else if (ekl) then
        call grad_copy (d_ik, d_il, li, lk, blen)
        call grad_copy (p_ik, p_il, li, lk, blen)
      else
        call ga_get (g_dens(2), ilo, ihi, llo, lhi, d_il, blen)
        call ga_get (g_dens(3), ilo, ihi, llo, lhi, p_il, blen)
      endif

C     j,k
      if (ogetjk) then
         if (ekl) then
            call grad_copy (d_jl, d_jk, lj, ll, blen)
            call grad_copy (p_jl, p_jk, lj, ll, blen)
         else if (eij) then
            call grad_copy (d_ik, d_jk, li, lk, blen)
            call grad_copy (p_ik, p_jk, li, lk, blen)
         else if (eik) then
            call grad_transp (d_ij, d_jk, li, lj, blen)
            call grad_transp (p_ij, p_jk, li, lj, blen)
         else if (ejl) then
            call grad_transp (d_kl, d_jk, lk, ll, blen)
            call grad_transp (p_kl, p_jk, lk, ll, blen)
         else if (eij.and.ekl) then
            call grad_copy (d_il, d_jk, li, ll, blen)
            call grad_copy (p_il, p_jk, li, ll, blen)
         else if (ejl.and.eik) then
            call grad_transp (d_il, d_jk, li, ll, blen)
            call grad_transp (p_il, p_jk, li, ll, blen)
         else
            call ga_get (g_dens(2), jlo, jhi, klo, khi, d_jk, blen)
            call ga_get (g_dens(3), jlo, jhi, klo, khi, p_jk, blen)
         endif
      endif
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_get_dens_ump2 (
     $           da_ij, da_kl, da_ik, da_jl, da_il, da_jk, 
     $           db_ij, db_kl, db_ik, db_jl, db_il, db_jk, 
     $           pa_ij, pa_kl, pa_ik, pa_jl, pa_il, pa_jk, 
     $           pb_ij, pb_kl, pb_ik, pb_jl, pb_il, pb_jk, 
     $           blen,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           ogetij, ogetik, ogetjk,
     $           g_dens)
C
C     get blocks from density matrices and make sure no block is
C     transferred twice
C
C     (a,b) denotes spin, d is the HF part and p contains (twice) the 
C     second order correction
C
      implicit none

#include "global.fh"

      integer g_dens, blen, 
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     li, lj, lk, ll
      logical ogetij, ogetik, ogetjk

      double precision da_ij, da_kl, da_ik, da_jl, da_il, da_jk,
     $                 db_ij, db_kl, db_ik, db_jl, db_il, db_jk,
     $                 pa_ij, pa_kl, pa_ik, pa_jl, pa_il, pa_jk,
     $                 pb_ij, pb_kl, pb_ik, pb_jl, pb_il, pb_jk

      dimension g_dens(5)
      dimension da_ij(blen,blen),da_kl(blen,blen),da_ik(blen,blen),
     $          da_jl(blen,blen),da_il(blen,blen),da_jk(blen,blen),
     $          db_ij(blen,blen),db_kl(blen,blen),db_ik(blen,blen),
     $          db_jl(blen,blen),db_il(blen,blen),db_jk(blen,blen),
     $          pa_ij(blen,blen),pa_kl(blen,blen),pa_ik(blen,blen),
     $          pa_jl(blen,blen),pa_il(blen,blen),pa_jk(blen,blen),
     $          pb_ij(blen,blen),pb_kl(blen,blen),pb_ik(blen,blen),
     $          pb_jl(blen,blen),pb_il(blen,blen),pb_jk(blen,blen)

      logical eij, ejk, ekl, eil, eik, ejl

      eij = (ilo.eq.jlo)
      ejk = (jlo.eq.klo)
      ekl = (klo.eq.llo)
      eil = (ilo.eq.llo)
      eik = (ilo.eq.klo)
      ejl = (jlo.eq.llo)

      li = ihi - ilo + 1
      lj = jhi - jlo + 1
      lk = khi - klo + 1
      ll = lhi - llo + 1

C     i,j
      if (ogetij) then
         call ga_get (g_dens(2), ilo, ihi, jlo, jhi, da_ij, blen)
         call ga_get (g_dens(3), ilo, ihi, jlo, jhi, db_ij, blen)
         call ga_get (g_dens(4), ilo, ihi, jlo, jhi, pa_ij, blen)
         call ga_get (g_dens(5), ilo, ihi, jlo, jhi, pb_ij, blen)
      endif

C     k,l
      if (eik.and.ejl) then
        call grad_copy (da_ij, da_kl, li, lj, blen)
        call grad_copy (db_ij, db_kl, li, lj, blen)
        call grad_copy (pa_ij, pa_kl, li, lj, blen)
        call grad_copy (pb_ij, pb_kl, li, lj, blen)
      else if (eil.and.ejk) then
        call grad_transp (da_ij, da_kl, li, lj, blen)
        call grad_transp (db_ij, db_kl, li, lj, blen)
        call grad_transp (pa_ij, pa_kl, li, lj, blen)
        call grad_transp (pb_ij, pb_kl, li, lj, blen)
      else
        call ga_get (g_dens(2), klo, khi, llo, lhi, da_kl, blen)
        call ga_get (g_dens(3), klo, khi, llo, lhi, db_kl, blen)
        call ga_get (g_dens(4), klo, khi, llo, lhi, pa_kl, blen)
        call ga_get (g_dens(5), klo, khi, llo, lhi, pb_kl, blen)
      endif

C     i,k
      if (ogetik) then
         if (ejk) then
            call grad_copy (da_ij, da_ik, li, lj, blen)
            call grad_copy (db_ij, db_ik, li, lj, blen)
            call grad_copy (pa_ij, pa_ik, li, lj, blen)
            call grad_copy (pb_ij, pb_ik, li, lj, blen)
         else if (eil) then
            call grad_transp (da_kl, da_ik, lk, ll, blen)
            call grad_transp (db_kl, db_ik, lk, ll, blen)
            call grad_transp (pa_kl, pa_ik, lk, ll, blen)
            call grad_transp (pb_kl, pb_ik, lk, ll, blen)
         else
            call ga_get (g_dens(2), ilo, ihi, klo, khi, da_ik, blen)
            call ga_get (g_dens(3), ilo, ihi, klo, khi, db_ik, blen)
            call ga_get (g_dens(4), ilo, ihi, klo, khi, pa_ik, blen)
            call ga_get (g_dens(5), ilo, ihi, klo, khi, pb_ik, blen)
         endif
      endif

C     j,l
      if (ejk) then
        call grad_copy (da_kl, da_jl, lk, ll, blen)
        call grad_copy (db_kl, db_jl, lk, ll, blen)
        call grad_copy (pa_kl, pa_jl, lk, ll, blen)
        call grad_copy (pb_kl, pb_jl, lk, ll, blen)
      else if (eij.and.ejl) then
        call grad_copy (da_ij, da_jl, li, lj, blen)
        call grad_copy (db_ij, db_jl, li, lj, blen)
        call grad_copy (pa_ij, pa_jl, li, lj, blen)
        call grad_copy (pb_ij, pb_jl, li, lj, blen)
      else if (eil) then
        call grad_transp (da_ij, da_jl, li, lj, blen)
        call grad_transp (db_ij, db_jl, li, lj, blen)
        call grad_transp (pa_ij, pa_jl, li, lj, blen)
        call grad_transp (pb_ij, pb_jl, li, lj, blen)
      else if (eij.and.ekl) then
        call grad_copy (da_ik, da_jl, li, lk, blen)
        call grad_copy (db_ik, db_jl, li, lk, blen)
        call grad_copy (pa_ik, pa_jl, li, lk, blen)
        call grad_copy (pb_ik, pb_jl, li, lk, blen)
      else if (ejk.and.eil) then
        call grad_transp (da_ik, da_jl, li, lk, blen)
        call grad_transp (db_ik, db_jl, li, lk, blen)
        call grad_transp (pa_ik, pa_jl, li, lk, blen)
        call grad_transp (pb_ik, pb_jl, li, lk, blen)
      else
        call ga_get (g_dens(2), jlo, jhi, llo, lhi, da_jl, blen)
        call ga_get (g_dens(3), jlo, jhi, llo, lhi, db_jl, blen)
        call ga_get (g_dens(4), jlo, jhi, llo, lhi, pa_jl, blen)
        call ga_get (g_dens(5), jlo, jhi, llo, lhi, pb_jl, blen)
      endif

C     i,l
      if (ejl) then
        call grad_copy (da_ij, da_il, li, lj, blen)
        call grad_copy (db_ij, db_il, li, lj, blen)
        call grad_copy (pa_ij, pa_il, li, lj, blen)
        call grad_copy (pb_ij, pb_il, li, lj, blen)
      else if (eik) then
        call grad_copy (da_kl, da_il, lk, ll, blen)
        call grad_copy (db_kl, db_il, lk, ll, blen)
        call grad_copy (pa_kl, pa_il, lk, ll, blen)
        call grad_copy (pb_kl, pb_il, lk, ll, blen)
      else if (eij) then
        call grad_copy (da_jl, da_il, lj, ll, blen)
        call grad_copy (db_jl, db_il, lj, ll, blen)
        call grad_copy (pa_jl, pa_il, lj, ll, blen)
        call grad_copy (pb_jl, pb_il, lj, ll, blen)
      else if (ekl) then
        call grad_copy (da_ik, da_il, li, lk, blen)
        call grad_copy (db_ik, db_il, li, lk, blen)
        call grad_copy (pa_ik, pa_il, li, lk, blen)
        call grad_copy (pb_ik, pb_il, li, lk, blen)
      else
        call ga_get (g_dens(2), ilo, ihi, llo, lhi, da_il, blen)
        call ga_get (g_dens(3), ilo, ihi, llo, lhi, db_il, blen)
        call ga_get (g_dens(4), ilo, ihi, llo, lhi, pa_il, blen)
        call ga_get (g_dens(5), ilo, ihi, llo, lhi, pb_il, blen)
      endif

C     j,k
      if (ogetjk) then
         if (ekl) then
            call grad_copy (da_jl, da_jk, lj, ll, blen)
            call grad_copy (db_jl, db_jk, lj, ll, blen)
            call grad_copy (pa_jl, pa_jk, lj, ll, blen)
            call grad_copy (pb_jl, pb_jk, lj, ll, blen)
         else if (eij) then
            call grad_copy (da_ik, da_jk, li, lk, blen)
            call grad_copy (db_ik, db_jk, li, lk, blen)
            call grad_copy (pa_ik, pa_jk, li, lk, blen)
            call grad_copy (pb_ik, pb_jk, li, lk, blen)
         else if (eik) then
            call grad_transp (da_ij, da_jk, li, lj, blen)
            call grad_transp (db_ij, db_jk, li, lj, blen)
            call grad_transp (pa_ij, pa_jk, li, lj, blen)
            call grad_transp (pb_ij, pb_jk, li, lj, blen)
         else if (ejl) then
            call grad_transp (da_kl, da_jk, lk, ll, blen)
            call grad_transp (db_kl, db_jk, lk, ll, blen)
            call grad_transp (pa_kl, pa_jk, lk, ll, blen)
            call grad_transp (pb_kl, pb_jk, lk, ll, blen)
         else if (eij.and.ekl) then
            call grad_copy (da_il, da_jk, li, ll, blen)
            call grad_copy (db_il, db_jk, li, ll, blen)
            call grad_copy (pa_il, pa_jk, li, ll, blen)
            call grad_copy (pb_il, pb_jk, li, ll, blen)
         else if (ejl.and.eik) then
            call grad_transp (da_il, da_jk, li, ll, blen)
            call grad_transp (db_il, db_jk, li, ll, blen)
            call grad_transp (pa_il, pa_jk, li, ll, blen)
            call grad_transp (pb_il, pb_jk, li, ll, blen)
         else
            call ga_get (g_dens(2), jlo, jhi, klo, khi, da_jk, blen)
            call ga_get (g_dens(3), jlo, jhi, klo, khi, db_jk, blen)
            call ga_get (g_dens(4), jlo, jhi, klo, khi, pa_jk, blen)
            call ga_get (g_dens(5), jlo, jhi, klo, khi, pb_jk, blen)
         endif
      endif
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_copy (d_1, d_2, l1, l2, ld)
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
      subroutine grad_transp (d_1, d_2, l1, l2, ld)
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
