      subroutine fock_init_cmul(idim,jdim,ldim)
      implicit none
#include "cfockmul.fh"
      integer idim, jdim, ldim
c
c     Initialize array to eliminate integer multiplication in
c     Fock-build kernels
c
      integer i
c
      do i = 0,255
         lm(i) = ldim*i
         jm(i) = jdim*i
         im(i) = idim*i
      enddo
c
      end
      subroutine fock_2e_label(nfock, tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      integer nfock
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(neri)
      double precision dij(nfock,0:*)
      double precision dik(nfock,0:*)
      double precision dli(nfock,0:*)
      double precision djk(nfock,0:*)
      double precision dlj(nfock,0:*)
      double precision dlk(nfock,0:*)
      double precision fij(nfock,0:*)
      double precision fik(nfock,0:*)
      double precision fli(nfock,0:*)
      double precision fjk(nfock,0:*)
      double precision flj(nfock,0:*)
      double precision flk(nfock,0:*)
c     
      integer i, j, k, l, ind, v
      integer idim, jdim, ldim
      integer ij, lk, ik, li, lj, jk
      double precision g
c
      if (nfock .eq. 1) then
         call fock_2e_1_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else if (nfock .eq. 2) then
         call fock_2e_2_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else if (nfock .eq. 3) then
         call fock_2e_3_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else if (nfock .eq. 4) then
         call fock_2e_4_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else
         idim = ihi - ilo + 1
         jdim = jhi - jlo + 1
         ldim = lhi - llo + 1
         do ind = 1, neri
            g = eri(ind)
            if (abs(g) .gt. tol2e) then
               i = ilab(ind) - ilo
               j = jlab(ind) - jlo
               k = klab(ind) - klo
               l = llab(ind) - llo
c     
*     ij = i + j*idim
*     lk = l + k*ldim
*     ik = i + k*idim
*     li = l + i*ldim
*     lj = l + j*ldim
*     jk = j + k*jdim
               ij = i + im(j)
               lk = l + lm(k)
               ik = i + im(k)
               li = l + lm(i)
               lj = l + lm(j)
               jk = j + jm(k)
               do v = 1, nfock
                  fij(v,ij) = fij(v,ij) + g*dlk(v,lk)
                  flk(v,lk) = flk(v,lk) + g*dij(v,ij)
                  fik(v,ik) = fik(v,ik) + g*dlj(v,lj)
                  fli(v,li) = fli(v,li) + g*djk(v,jk)
                  flj(v,lj) = flj(v,lj) + g*dik(v,ik)
                  fjk(v,jk) = fjk(v,jk) + g*dli(v,li)
               end do
            end if
         end do
      endif
c
      end
      subroutine fock_2e_1_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(neri)
      double precision dij(0:*)
      double precision dik(0:*)
      double precision dli(0:*)
      double precision djk(0:*)
      double precision dlj(0:*)
      double precision dlk(0:*)
      double precision fij(0:*)
      double precision fik(0:*)
      double precision fli(0:*)
      double precision fjk(0:*)
      double precision flj(0:*)
      double precision flk(0:*)
c     
      integer i, j, k, l, ind
      double precision g
      integer idim, jdim, ldim
      integer ij, lk, ik, li, lj, jk
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
      ldim = lhi - llo + 1
c     
      do ind = 1, neri
         g = eri(ind)
         i = ilab(ind) - ilo
         j = jlab(ind) - jlo
         k = klab(ind) - klo
         l = llab(ind) - llo
c     
*      ij = i + j*idim
*      lk = l + k*ldim
*      ik = i + k*idim
*      li = l + i*ldim
*      lj = l + j*ldim
*      jk = j + k*jdim
         ij = i + im(j)
         lk = l + lm(k)
         ik = i + im(k)
         li = l + lm(i)
         lj = l + lm(j)
         jk = j + jm(k)
         fij(ij) = fij(ij) + g*dlk(lk)
         flk(lk) = flk(lk) + g*dij(ij)
         fik(ik) = fik(ik) + g*dlj(lj)
         fli(li) = fli(li) + g*djk(jk)
         flj(lj) = flj(lj) + g*dik(ik)
         fjk(jk) = fjk(jk) + g*dli(li)
      end do
c     
      end
      subroutine fock_2e_2_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(neri)
      double precision dij(2,0:*)
      double precision dik(2,0:*)
      double precision dli(2,0:*)
      double precision djk(2,0:*)
      double precision dlj(2,0:*)
      double precision dlk(2,0:*)
      double precision fij(2,0:*)
      double precision fik(2,0:*)
      double precision fli(2,0:*)
      double precision fjk(2,0:*)
      double precision flj(2,0:*)
      double precision flk(2,0:*)
c     
      integer i, j, k, l, ind
      double precision g
      integer idim, jdim, ldim
      integer ij, lk, ik, li, lj, jk
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
      ldim = lhi - llo + 1
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind) - ilo
            j = jlab(ind) - jlo
            k = klab(ind) - klo
            l = llab(ind) - llo
c
*            ij = i + j*idim
*            lk = l + k*ldim
*            ik = i + k*idim
*            li = l + i*ldim
*            lj = l + j*ldim
*            jk = j + k*jdim
            ij = i + im(j)
            lk = l + lm(k)
            ik = i + im(k)
            li = l + lm(i)
            lj = l + lm(j)
            jk = j + jm(k)
            fij(1,ij) = fij(1,ij) + g*dlk(1,lk)
            flk(1,lk) = flk(1,lk) + g*dij(1,ij)
            fik(1,ik) = fik(1,ik) + g*dlj(1,lj)
            fli(1,li) = fli(1,li) + g*djk(1,jk)
            flj(1,lj) = flj(1,lj) + g*dik(1,ik)
            fjk(1,jk) = fjk(1,jk) + g*dli(1,li)
            fij(2,ij) = fij(2,ij) + g*dlk(2,lk)
            flk(2,lk) = flk(2,lk) + g*dij(2,ij)
            fik(2,ik) = fik(2,ik) + g*dlj(2,lj)
            fli(2,li) = fli(2,li) + g*djk(2,jk)
            flj(2,lj) = flj(2,lj) + g*dik(2,ik)
            fjk(2,jk) = fjk(2,jk) + g*dli(2,li)
         end if
      end do
c
      end
      subroutine fock_2e_3_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(neri)
      double precision dij(3,0:*)
      double precision dik(3,0:*)
      double precision dli(3,0:*)
      double precision djk(3,0:*)
      double precision dlj(3,0:*)
      double precision dlk(3,0:*)
      double precision fij(3,0:*)
      double precision fik(3,0:*)
      double precision fli(3,0:*)
      double precision fjk(3,0:*)
      double precision flj(3,0:*)
      double precision flk(3,0:*)
c     
      integer i, j, k, l, ind
      double precision g
      integer idim, jdim, ldim
      integer ij, lk, ik, li, lj, jk
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
      ldim = lhi - llo + 1
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind) - ilo
            j = jlab(ind) - jlo
            k = klab(ind) - klo
            l = llab(ind) - llo
c
*            ij = i + j*idim
*            lk = l + k*ldim
*            ik = i + k*idim
*            li = l + i*ldim
*            lj = l + j*ldim
*            jk = j + k*jdim
            ij = i + im(j)
            lk = l + lm(k)
            ik = i + im(k)
            li = l + lm(i)
            lj = l + lm(j)
            jk = j + jm(k)
            fij(1,ij) = fij(1,ij) + g*dlk(1,lk)
            flk(1,lk) = flk(1,lk) + g*dij(1,ij)
            fik(1,ik) = fik(1,ik) + g*dlj(1,lj)
            fli(1,li) = fli(1,li) + g*djk(1,jk)
            flj(1,lj) = flj(1,lj) + g*dik(1,ik)
            fjk(1,jk) = fjk(1,jk) + g*dli(1,li)
            fij(2,ij) = fij(2,ij) + g*dlk(2,lk)
            flk(2,lk) = flk(2,lk) + g*dij(2,ij)
            fik(2,ik) = fik(2,ik) + g*dlj(2,lj)
            fli(2,li) = fli(2,li) + g*djk(2,jk)
            flj(2,lj) = flj(2,lj) + g*dik(2,ik)
            fjk(2,jk) = fjk(2,jk) + g*dli(2,li)
            fij(3,ij) = fij(3,ij) + g*dlk(3,lk)
            flk(3,lk) = flk(3,lk) + g*dij(3,ij)
            fik(3,ik) = fik(3,ik) + g*dlj(3,lj)
            fli(3,li) = fli(3,li) + g*djk(3,jk)
            flj(3,lj) = flj(3,lj) + g*dik(3,ik)
            fjk(3,jk) = fjk(3,jk) + g*dli(3,li)
         end if
      end do
c
      end
      subroutine fock_2e_4_label(tol2e, neri, 
     $     ilab, jlab, klab, llab, eri,
     $     ilo, ihi, jlo, jhi,
     $     klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(neri)
      double precision dij(4,0:*)
      double precision dik(4,0:*)
      double precision dli(4,0:*)
      double precision djk(4,0:*)
      double precision dlj(4,0:*)
      double precision dlk(4,0:*)
      double precision fij(4,0:*)
      double precision fik(4,0:*)
      double precision fli(4,0:*)
      double precision fjk(4,0:*)
      double precision flj(4,0:*)
      double precision flk(4,0:*)
c     
      integer i, j, k, l, ind
      double precision g
      integer idim, jdim, ldim
      integer ij, lk, ik, li, lj, jk
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
      ldim = lhi - llo + 1
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind) - ilo
            j = jlab(ind) - jlo
            k = klab(ind) - klo
            l = llab(ind) - llo
c
*            ij = i + j*idim
*            lk = l + k*ldim
*            ik = i + k*idim
*            li = l + i*ldim
*            lj = l + j*ldim
*            jk = j + k*jdim
            ij = i + im(j)
            lk = l + lm(k)
            ik = i + im(k)
            li = l + lm(i)
            lj = l + lm(j)
            jk = j + jm(k)
            fij(1,ij) = fij(1,ij) + g*dlk(1,lk)
            flk(1,lk) = flk(1,lk) + g*dij(1,ij)
            fik(1,ik) = fik(1,ik) + g*dlj(1,lj)
            fli(1,li) = fli(1,li) + g*djk(1,jk)
            flj(1,lj) = flj(1,lj) + g*dik(1,ik)
            fjk(1,jk) = fjk(1,jk) + g*dli(1,li)
            fij(2,ij) = fij(2,ij) + g*dlk(2,lk)
            flk(2,lk) = flk(2,lk) + g*dij(2,ij)
            fik(2,ik) = fik(2,ik) + g*dlj(2,lj)
            fli(2,li) = fli(2,li) + g*djk(2,jk)
            flj(2,lj) = flj(2,lj) + g*dik(2,ik)
            fjk(2,jk) = fjk(2,jk) + g*dli(2,li)
            fij(3,ij) = fij(3,ij) + g*dlk(3,lk)
            flk(3,lk) = flk(3,lk) + g*dij(3,ij)
            fik(3,ik) = fik(3,ik) + g*dlj(3,lj)
            fli(3,li) = fli(3,li) + g*djk(3,jk)
            flj(3,lj) = flj(3,lj) + g*dik(3,ik)
            fjk(3,jk) = fjk(3,jk) + g*dli(3,li)
            fij(4,ij) = fij(4,ij) + g*dlk(4,lk)
            flk(4,lk) = flk(4,lk) + g*dij(4,ij)
            fik(4,ik) = fik(4,ik) + g*dlj(4,lj)
            fli(4,li) = fli(4,li) + g*djk(4,jk)
            flj(4,lj) = flj(4,lj) + g*dik(4,ik)
            fjk(4,jk) = fjk(4,jk) + g*dli(4,li)
         end if
      end do
c
      end
c!!
c!!   READ THIS !!!!!!!!!!!!
c!!
c!!     Mod routines differ from above routines in that
c!!
c!!     1) have labels(4,neri) rather than i(neri), j(neri), ...
c!!     2) labels are offsets within the range, rather than
c!!        full labels (i.e., there is no need to subtract ilo)
c!!     3) No need to pass ranges since im/jm/lm are used
c!!
      subroutine fock_2e_mod_label(nfock, tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      integer nfock
      double precision tol2e
      integer neri
      integer labels(*)
      double precision eri(neri)
      double precision dij(nfock,0:*)
      double precision dik(nfock,0:*)
      double precision dli(nfock,0:*)
      double precision djk(nfock,0:*)
      double precision dlj(nfock,0:*)
      double precision dlk(nfock,0:*)
      double precision fij(nfock,0:*)
      double precision fik(nfock,0:*)
      double precision fli(nfock,0:*)
      double precision fjk(nfock,0:*)
      double precision flj(nfock,0:*)
      double precision flk(nfock,0:*)
c     
      integer i, j, k, l, ind, v, ind4
      integer ij, lk, ik, li, lj, jk
      double precision g
c
      if (nfock .eq. 1) then
         call fock_2e_mod_1_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else if (nfock .eq. 2) then
         call fock_2e_mod_2_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else if (nfock .eq. 3) then
         call fock_2e_mod_3_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else if (nfock .eq. 4) then
         call fock_2e_mod_4_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
      else
         ind4 = 1
         do ind = 1, neri
            g = eri(ind)
            if (abs(g) .gt. tol2e) then
               i = labels(  ind4)
               j = labels(1+ind4)
               k = labels(2+ind4)
               l = labels(3+ind4)
c     
               ij = i + im(j)
               lk = l + lm(k)
               ik = i + im(k)
               li = l + lm(i)
               lj = l + lm(j)
               jk = j + jm(k)
               do v = 1, nfock
                  fij(v,ij) = fij(v,ij) + g*dlk(v,lk)
                  flk(v,lk) = flk(v,lk) + g*dij(v,ij)
                  fik(v,ik) = fik(v,ik) + g*dlj(v,lj)
                  fli(v,li) = fli(v,li) + g*djk(v,jk)
                  flj(v,lj) = flj(v,lj) + g*dik(v,ik)
                  fjk(v,jk) = fjk(v,jk) + g*dli(v,li)
               end do
            end if
            ind4 = ind4 + 4
         end do
      endif
c
      end
      subroutine fock_2e_mod_1_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer labels(*)
      double precision eri(neri)
      double precision dij(0:*)
      double precision dik(0:*)
      double precision dli(0:*)
      double precision djk(0:*)
      double precision dlj(0:*)
      double precision dlk(0:*)
      double precision fij(0:*)
      double precision fik(0:*)
      double precision fli(0:*)
      double precision fjk(0:*)
      double precision flj(0:*)
      double precision flk(0:*)
c     
      integer i, j, k, l, ind, ind4
      double precision g
      integer ij, lk, ik, li, lj, jk
c
      ind4 = 1
      do ind = 1, neri
         g = eri(ind)
         i = labels(  ind4)
         j = labels(1+ind4)
         k = labels(2+ind4)
         l = labels(3+ind4)
         ind4 = ind4 + 4
c     
         ij = i + im(j)
         lk = l + lm(k)
         ik = i + im(k)
         li = l + lm(i)
         lj = l + lm(j)
         jk = j + jm(k)
         fij(ij) = fij(ij) + g*dlk(lk)
         flk(lk) = flk(lk) + g*dij(ij)
         fik(ik) = fik(ik) + g*dlj(lj)
         fli(li) = fli(li) + g*djk(jk)
         flj(lj) = flj(lj) + g*dik(ik)
         fjk(jk) = fjk(jk) + g*dli(li)
      end do
c     
      end
      subroutine fock_2e_mod_2_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer labels(*)
      double precision eri(neri)
      double precision dij(2,0:*)
      double precision dik(2,0:*)
      double precision dli(2,0:*)
      double precision djk(2,0:*)
      double precision dlj(2,0:*)
      double precision dlk(2,0:*)
      double precision fij(2,0:*)
      double precision fik(2,0:*)
      double precision fli(2,0:*)
      double precision fjk(2,0:*)
      double precision flj(2,0:*)
      double precision flk(2,0:*)
c     
      integer i, j, k, l, ind, ind4
      double precision g
      integer ij, lk, ik, li, lj, jk
c
      ind4 = 1
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(  ind4)
            j = labels(1+ind4)
            k = labels(2+ind4)
            l = labels(3+ind4)
c
            ij = i + im(j)
            lk = l + lm(k)
            ik = i + im(k)
            li = l + lm(i)
            lj = l + lm(j)
            jk = j + jm(k)
            fij(1,ij) = fij(1,ij) + g*dlk(1,lk)
            flk(1,lk) = flk(1,lk) + g*dij(1,ij)
            fik(1,ik) = fik(1,ik) + g*dlj(1,lj)
            fli(1,li) = fli(1,li) + g*djk(1,jk)
            flj(1,lj) = flj(1,lj) + g*dik(1,ik)
            fjk(1,jk) = fjk(1,jk) + g*dli(1,li)
            fij(2,ij) = fij(2,ij) + g*dlk(2,lk)
            flk(2,lk) = flk(2,lk) + g*dij(2,ij)
            fik(2,ik) = fik(2,ik) + g*dlj(2,lj)
            fli(2,li) = fli(2,li) + g*djk(2,jk)
            flj(2,lj) = flj(2,lj) + g*dik(2,ik)
            fjk(2,jk) = fjk(2,jk) + g*dli(2,li)
         end if
         ind4 = ind4 + 4
      end do
c
      end
      subroutine fock_2e_mod_3_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer labels(*)
      double precision eri(neri)
      double precision dij(3,0:*)
      double precision dik(3,0:*)
      double precision dli(3,0:*)
      double precision djk(3,0:*)
      double precision dlj(3,0:*)
      double precision dlk(3,0:*)
      double precision fij(3,0:*)
      double precision fik(3,0:*)
      double precision fli(3,0:*)
      double precision fjk(3,0:*)
      double precision flj(3,0:*)
      double precision flk(3,0:*)
c     
      integer i, j, k, l, ind, ind4
      double precision g
      integer ij, lk, ik, li, lj, jk
c
      ind4 = 1
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(  ind4)
            j = labels(1+ind4)
            k = labels(2+ind4)
            l = labels(3+ind4)
c
            ij = i + im(j)
            lk = l + lm(k)
            ik = i + im(k)
            li = l + lm(i)
            lj = l + lm(j)
            jk = j + jm(k)
            fij(1,ij) = fij(1,ij) + g*dlk(1,lk)
            flk(1,lk) = flk(1,lk) + g*dij(1,ij)
            fik(1,ik) = fik(1,ik) + g*dlj(1,lj)
            fli(1,li) = fli(1,li) + g*djk(1,jk)
            flj(1,lj) = flj(1,lj) + g*dik(1,ik)
            fjk(1,jk) = fjk(1,jk) + g*dli(1,li)
            fij(2,ij) = fij(2,ij) + g*dlk(2,lk)
            flk(2,lk) = flk(2,lk) + g*dij(2,ij)
            fik(2,ik) = fik(2,ik) + g*dlj(2,lj)
            fli(2,li) = fli(2,li) + g*djk(2,jk)
            flj(2,lj) = flj(2,lj) + g*dik(2,ik)
            fjk(2,jk) = fjk(2,jk) + g*dli(2,li)
            fij(3,ij) = fij(3,ij) + g*dlk(3,lk)
            flk(3,lk) = flk(3,lk) + g*dij(3,ij)
            fik(3,ik) = fik(3,ik) + g*dlj(3,lj)
            fli(3,li) = fli(3,li) + g*djk(3,jk)
            flj(3,lj) = flj(3,lj) + g*dik(3,ik)
            fjk(3,jk) = fjk(3,jk) + g*dli(3,li)
         end if
         ind4 = ind4 + 4
      end do
c
      end
      subroutine fock_2e_mod_4_label(tol2e, neri, 
     $     labels, eri,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_lab.F,v 1.6 1996-02-13 19:57:51 rg240 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.  It is ASSUMED that
c     all integrals in the input record will contribute only to the
c     six blocks of the given fock matrix
c     
      double precision tol2e
      integer neri
      integer labels(*)
      double precision eri(neri)
      double precision dij(4,0:*)
      double precision dik(4,0:*)
      double precision dli(4,0:*)
      double precision djk(4,0:*)
      double precision dlj(4,0:*)
      double precision dlk(4,0:*)
      double precision fij(4,0:*)
      double precision fik(4,0:*)
      double precision fli(4,0:*)
      double precision fjk(4,0:*)
      double precision flj(4,0:*)
      double precision flk(4,0:*)
c     
      integer i, j, k, l, ind, ind4
      double precision g
      integer ij, lk, ik, li, lj, jk
c
      ind4 = 1
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(  ind4)
            j = labels(1+ind4)
            k = labels(2+ind4)
            l = labels(3+ind4)
c
            ij = i + im(j)
            lk = l + lm(k)
            ik = i + im(k)
            li = l + lm(i)
            lj = l + lm(j)
            jk = j + jm(k)
            fij(1,ij) = fij(1,ij) + g*dlk(1,lk)
            flk(1,lk) = flk(1,lk) + g*dij(1,ij)
            fik(1,ik) = fik(1,ik) + g*dlj(1,lj)
            fli(1,li) = fli(1,li) + g*djk(1,jk)
            flj(1,lj) = flj(1,lj) + g*dik(1,ik)
            fjk(1,jk) = fjk(1,jk) + g*dli(1,li)
            fij(2,ij) = fij(2,ij) + g*dlk(2,lk)
            flk(2,lk) = flk(2,lk) + g*dij(2,ij)
            fik(2,ik) = fik(2,ik) + g*dlj(2,lj)
            fli(2,li) = fli(2,li) + g*djk(2,jk)
            flj(2,lj) = flj(2,lj) + g*dik(2,ik)
            fjk(2,jk) = fjk(2,jk) + g*dli(2,li)
            fij(3,ij) = fij(3,ij) + g*dlk(3,lk)
            flk(3,lk) = flk(3,lk) + g*dij(3,ij)
            fik(3,ik) = fik(3,ik) + g*dlj(3,lj)
            fli(3,li) = fli(3,li) + g*djk(3,jk)
            flj(3,lj) = flj(3,lj) + g*dik(3,ik)
            fjk(3,jk) = fjk(3,jk) + g*dli(3,li)
            fij(4,ij) = fij(4,ij) + g*dlk(4,lk)
            flk(4,lk) = flk(4,lk) + g*dij(4,ij)
            fik(4,ik) = fik(4,ik) + g*dlj(4,lj)
            fli(4,li) = fli(4,li) + g*djk(4,jk)
            flj(4,lj) = flj(4,lj) + g*dik(4,ik)
            fjk(4,jk) = fjk(4,jk) + g*dli(4,li)
         end if
         ind4 = ind4 + 4
      end do
c
      end
