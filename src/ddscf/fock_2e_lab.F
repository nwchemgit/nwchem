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
      do i = 0,maxim1
         lm(i) = ldim*i
         jm(i) = jdim*i
         im(i) = idim*i
      enddo
c
      end
      subroutine fock_init_rep_cmul(nbf)
      implicit none
#include "cfockmul.fh"
      integer nbf
c
c     Initialize array to eliminate integer multiplication in
c     Fock-build kernels ... replicated version
c
      integer i
c
      do i = 1,nbf
         im1(i) = (i-1)*nbf
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
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
      subroutine fock_2e_rep_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, ilab, jlab, klab, llab, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.
c     
c     The _rep_ routines assume that each process has the entire 
c     fock matrix and that we are not adding into distinct patches
c     of the matrix.  We also must add in the j/k scaling factors.
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      double precision eri(neri)
      double precision fock(nfock,nbf*nbf), dens(nfock,nbf*nbf)
      double precision jfac(nfock), kfac(nfock)
c     
      integer i, j, k, l, ind, v
      integer ij, lk, ik, li, lj, jk
      double precision g, gk, gj
c
      if (nfock .eq. 1) then
         call fock_2e_rep_1_label(nfock, nbf, jfac, kfac, tol2e,
     $        neri, ilab, jlab, klab, llab, eri, dens, fock)
         return
      else if (nfock .eq. 2) then
         call fock_2e_rep_2_label(nfock, nbf, jfac, kfac, tol2e,
     $        neri, ilab, jlab, klab, llab, eri, dens, fock)
         return
      elseif (nfock .eq. 3) then
         call fock_2e_rep_3_label(nfock, nbf, jfac, kfac, tol2e,
     $        neri, ilab, jlab, klab, llab, eri, dens, fock)
         return
      else if (nfock .eq. 4) then
         call fock_2e_rep_4_label(nfock, nbf, jfac, kfac, tol2e,
     $        neri, ilab, jlab, klab, llab, eri, dens, fock)
         return
      endif
c     
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind)
            j = jlab(ind)
            k = klab(ind)
            l = llab(ind)
*debug
*            if (i.le.0 .or. i.gt.nbf) call errquit(' i bad', i)
*            if (j.le.0 .or. j.gt.nbf) call errquit(' j bad', j)
*            if (l.le.0 .or. l.gt.nbf) call errquit(' l bad', l)
*            if (k.le.0 .or. k.gt.nbf) call errquit(' k bad', k)
*
            ij = i + im1(j)
            lj = l + im1(j)
            lk = l + im1(k)
            ik = i + im1(k)
            jk = j + im1(k)
            li = l + im1(i)
c
            do v = 1, nfock
               gj = g*jfac(v)
               gk = g*kfac(v)
               fock(v,ij) = fock(v,ij) + gj*dens(v,lk)
               fock(v,lk) = fock(v,lk) + gj*dens(v,ij)
               fock(v,ik) = fock(v,ik) + gk*dens(v,lj)
               fock(v,li) = fock(v,li) + gk*dens(v,jk)
               fock(v,lj) = fock(v,lj) + gk*dens(v,ik)
               fock(v,jk) = fock(v,jk) + gk*dens(v,li)
            end do
         end if
      end do
c     
      end
      subroutine fock_2e_rep_1_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, ilab, jlab, klab, llab, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.
c     
c     The _rep_ routines assume that each process has the entire 
c     fock matrix and that we are not adding into distinct patches
c     of the matrix.  We also must add in the j/k scaling factors.
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      double precision eri(neri)
      double precision fock(*), dens(*)
      double precision jfac, kfac
c     
      integer i, j, k, l, ind
      integer ij, kl, ik, il, jl, jk
      double precision g, gk, gj, fij, fkl, fik, fil, fjl, fjk
c     
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            gj = g*jfac
            gk = g*kfac
c
            i = ilab(ind)
            j = jlab(ind)
            k = klab(ind)
            l = llab(ind)
c
c
c     Logic behind this choice is that indices vary mostly with
c     i>=j and k>=l with kl varying fastest ... improved cache hits?
c     If change this mapping must change overwrite logic below
c
            ij = im1(i) + j
            kl = im1(k) + l
            if (ij .eq. kl) gj = gj + gj
c
c     Using overwrite logic ... separates stores/loads, and enables
c     partial pipelining of */+ pairs ... at least in principle
c
            fij = fock(ij) + gj*dens(kl)
            fkl = fock(kl) + gj*dens(ij)
            if (i .eq. j) gk = gk + gk
            fock(ij) = fij
            fock(kl) = fkl
c
            if (k .eq. l) gk = gk + gk
c
            ik = im1(i) + k
            il = im1(i) + l
            jl = im1(j) + l
            jk = im1(j) + k
c
            fik = fock(ik) + gk*dens(jl)
            fil = fock(il) + gk*dens(jk)
            fjk = fock(jk) + gk*dens(il)
            fjl = fock(jl) + gk*dens(ik)
            fock(ik) = fik
            fock(il) = fil
            fock(jl) = fjl
            fock(jk) = fjk
c
c     Code without overwrite logic
c
c$$$            fock(ij) = fock(ij) + gj*dens(kl)
c$$$            fock(kl) = fock(kl) + gj*dens(ij)
c$$$            fock(ik) = fock(ik) + gk*dens(jl)
c$$$            fock(il) = fock(il) + gk*dens(jk)
c$$$            fock(jl) = fock(jl) + gk*dens(ik)
c$$$            fock(jk) = fock(jk) + gk*dens(il)
         end if
      end do
c     
      end
      subroutine fock_2e_rep_2_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, ilab, jlab, klab, llab, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.
c     
c     The _rep_ routines assume that each process has the entire 
c     fock matrix and that we are not adding into distinct patches
c     of the matrix.  We also must add in the j/k scaling factors.
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      double precision eri(neri)
      double precision fock(nfock,nbf*nbf), dens(nfock,nbf*nbf)
      double precision jfac(nfock), kfac(nfock)
c     
      integer i, j, k, l, ind
      integer ij, lk, ik, li, lj, jk
      double precision g, gk1, gj1, gk2, gj2
c     
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind)
            j = jlab(ind)
            k = klab(ind)
            l = llab(ind)
c
            ij = i + im1(j)
            lj = l + im1(j)
            lk = l + im1(k)
            ik = i + im1(k)
            jk = j + im1(k)
            li = l + im1(i)
c
            gj1 = g*jfac(1)
            gk1 = g*kfac(1)
            fock(1,ij) = fock(1,ij) + gj1*dens(1,lk)
            fock(1,lk) = fock(1,lk) + gj1*dens(1,ij)
            fock(1,ik) = fock(1,ik) + gk1*dens(1,lj)
            fock(1,li) = fock(1,li) + gk1*dens(1,jk)
            fock(1,lj) = fock(1,lj) + gk1*dens(1,ik)
            fock(1,jk) = fock(1,jk) + gk1*dens(1,li)
c
            gj2 = g*jfac(2)
            gk2 = g*kfac(2)
            fock(2,ij) = fock(2,ij) + gj2*dens(2,lk)
            fock(2,lk) = fock(2,lk) + gj2*dens(2,ij)
            fock(2,ik) = fock(2,ik) + gk2*dens(2,lj)
            fock(2,li) = fock(2,li) + gk2*dens(2,jk)
            fock(2,lj) = fock(2,lj) + gk2*dens(2,ik)
            fock(2,jk) = fock(2,jk) + gk2*dens(2,li)
         end if
      end do
c     
      end
      subroutine fock_2e_rep_3_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, ilab, jlab, klab, llab, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.
c     
c     The _rep_ routines assume that each process has the entire 
c     fock matrix and that we are not adding into distinct patches
c     of the matrix.  We also must add in the j/k scaling factors.
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      double precision eri(neri)
      double precision fock(nfock,nbf*nbf), dens(nfock,nbf*nbf)
      double precision jfac(nfock), kfac(nfock)
c     
      integer i, j, k, l, ind
      integer ij, lk, ik, li, lj, jk
      double precision g, gk1, gj1, gk2, gj2, gk3, gj3
c     
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind)
            j = jlab(ind)
            k = klab(ind)
            l = llab(ind)
c
            ij = i + im1(j)
            lj = l + im1(j)
            lk = l + im1(k)
            ik = i + im1(k)
            jk = j + im1(k)
            li = l + im1(i)
c
            gj1 = g*jfac(1)
            gk1 = g*kfac(1)
            fock(1,ij) = fock(1,ij) + gj1*dens(1,lk)
            fock(1,lk) = fock(1,lk) + gj1*dens(1,ij)
            fock(1,ik) = fock(1,ik) + gk1*dens(1,lj)
            fock(1,li) = fock(1,li) + gk1*dens(1,jk)
            fock(1,lj) = fock(1,lj) + gk1*dens(1,ik)
            fock(1,jk) = fock(1,jk) + gk1*dens(1,li)
c
            gj2 = g*jfac(2)
            gk2 = g*kfac(2)
            fock(2,ij) = fock(2,ij) + gj2*dens(2,lk)
            fock(2,lk) = fock(2,lk) + gj2*dens(2,ij)
            fock(2,ik) = fock(2,ik) + gk2*dens(2,lj)
            fock(2,li) = fock(2,li) + gk2*dens(2,jk)
            fock(2,lj) = fock(2,lj) + gk2*dens(2,ik)
            fock(2,jk) = fock(2,jk) + gk2*dens(2,li)
c
            gj3 = g*jfac(3)
            gk3 = g*kfac(3)
            fock(3,ij) = fock(3,ij) + gj3*dens(3,lk)
            fock(3,lk) = fock(3,lk) + gj3*dens(3,ij)
            fock(3,ik) = fock(3,ik) + gk3*dens(3,lj)
            fock(3,li) = fock(3,li) + gk3*dens(3,jk)
            fock(3,lj) = fock(3,lj) + gk3*dens(3,ik)
            fock(3,jk) = fock(3,jk) + gk3*dens(3,li)
         end if
      end do
c     
      end
      subroutine fock_2e_rep_4_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, ilab, jlab, klab, llab, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.
c     
c     The _rep_ routines assume that each process has the entire 
c     fock matrix and that we are not adding into distinct patches
c     of the matrix.  We also must add in the j/k scaling factors.
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer ilab(neri),jlab(neri),klab(neri),llab(neri)
      double precision eri(neri)
      double precision fock(nfock,nbf*nbf), dens(nfock,nbf*nbf)
      double precision jfac(nfock), kfac(nfock)
c     
      integer i, j, k, l, ind
      integer ij, lk, ik, li, lj, jk
      double precision g, gk1, gj1, gk2, gj2, gk3, gj3, gk4, gj4
c     
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = ilab(ind)
            j = jlab(ind)
            k = klab(ind)
            l = llab(ind)
c
            ij = i + im1(j)
            lj = l + im1(j)
            lk = l + im1(k)
            ik = i + im1(k)
            jk = j + im1(k)
            li = l + im1(i)
c
            gj1 = g*jfac(1)
            gk1 = g*kfac(1)
            fock(1,ij) = fock(1,ij) + gj1*dens(1,lk)
            fock(1,lk) = fock(1,lk) + gj1*dens(1,ij)
            fock(1,ik) = fock(1,ik) + gk1*dens(1,lj)
            fock(1,li) = fock(1,li) + gk1*dens(1,jk)
            fock(1,lj) = fock(1,lj) + gk1*dens(1,ik)
            fock(1,jk) = fock(1,jk) + gk1*dens(1,li)
c
            gj2 = g*jfac(2)
            gk2 = g*kfac(2)
            fock(2,ij) = fock(2,ij) + gj2*dens(2,lk)
            fock(2,lk) = fock(2,lk) + gj2*dens(2,ij)
            fock(2,ik) = fock(2,ik) + gk2*dens(2,lj)
            fock(2,li) = fock(2,li) + gk2*dens(2,jk)
            fock(2,lj) = fock(2,lj) + gk2*dens(2,ik)
            fock(2,jk) = fock(2,jk) + gk2*dens(2,li)
c
            gj3 = g*jfac(3)
            gk3 = g*kfac(3)
            fock(3,ij) = fock(3,ij) + gj3*dens(3,lk)
            fock(3,lk) = fock(3,lk) + gj3*dens(3,ij)
            fock(3,ik) = fock(3,ik) + gk3*dens(3,lj)
            fock(3,li) = fock(3,li) + gk3*dens(3,jk)
            fock(3,lj) = fock(3,lj) + gk3*dens(3,ik)
            fock(3,jk) = fock(3,jk) + gk3*dens(3,li)
c
            gj4 = g*jfac(4)
            gk4 = g*kfac(4)
            fock(4,ij) = fock(4,ij) + gj4*dens(4,lk)
            fock(4,lk) = fock(4,lk) + gj4*dens(4,ij)
            fock(4,ik) = fock(4,ik) + gk4*dens(4,lj)
            fock(4,li) = fock(4,li) + gk4*dens(4,jk)
            fock(4,lj) = fock(4,lj) + gk4*dens(4,ik)
            fock(4,jk) = fock(4,jk) + gk4*dens(4,li)
         end if
      end do
c     
      end
      subroutine fock_2e_rep_mod_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
c     Add scaled labelled integrals into fock matrix blocks ... 
c     symmetrization will be needed later.
c
c     MOD routines again use label(4,*) and indices are
c     offset from zero.
c     
c     The _rep_ routines assume that each process has the entire 
c     fock matrix and that we are not adding into distinct patches
c     of the matrix.  We also must add in the j/k scaling factors.
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer labels(4,neri)
      double precision eri(neri)
      double precision fock(nfock,0:nbf*nbf-1), dens(nfock,0:nbf*nbf-1)
      double precision jfac(nfock), kfac(nfock)
c     
      integer i, j, k, l, ind, v
      integer ij, lk, ik, li, lj, jk
      double precision g, gk, gj
c
      if (nfock .eq. 1) then
         call fock_2e_rep_mod_1_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
         return
      else if (nfock .eq. 2) then
         call fock_2e_rep_mod_2_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
         return
      else if (nfock .eq. 4) then
         call fock_2e_rep_mod_4_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
         return
      endif
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
*     debug
*            if (i.lt.0 .or. i.ge.nbf) call errquit(' i bad', i)
*            if (j.lt.0 .or. j.ge.nbf) call errquit(' j bad', j)
*            if (l.lt.0 .or. l.ge.nbf) call errquit(' l bad', l)
*            if (k.lt.0 .or. k.ge.nbf) call errquit(' k bad', k)

            ij = i + im(j)      ! im NOT im1 since go from 0
            lj = l + im(j)
            lk = l + im(k)
            ik = i + im(k)
            jk = j + im(k)
            li = l + im(i)
c
            do v = 1, nfock
               gj = g*jfac(v)
               gk = g*kfac(v)
               fock(v,ij) = fock(v,ij) + gj*dens(v,lk)
               fock(v,lk) = fock(v,lk) + gj*dens(v,ij)
               fock(v,ik) = fock(v,ik) + gk*dens(v,lj)
               fock(v,li) = fock(v,li) + gk*dens(v,jk)
               fock(v,lj) = fock(v,lj) + gk*dens(v,ik)
               fock(v,jk) = fock(v,jk) + gk*dens(v,li)
            end do
         end if
      end do
c     
      end
      subroutine fock_2e_rep_mod_1_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer labels(4,neri)
      double precision eri(neri)
      double precision fock(0:*), dens(0:*)
      double precision jfac, kfac
c     
      integer i, j, k, l, ind
      integer ij, ik, jk, kl, il, jl
      double precision g, gk, gj, fij, fkl, fik, fil, fjl, fjk
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            gj = g*jfac
            gk = g*kfac
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
*     debug
*            if (i.lt.0 .or. i.ge.nbf) call errquit(' i bad', i)
*            if (j.lt.0 .or. j.ge.nbf) call errquit(' j bad', j)
*            if (l.lt.0 .or. l.ge.nbf) call errquit(' l bad', l)
*            if (k.lt.0 .or. k.ge.nbf) call errquit(' k bad', k)
c
c     Logic behind this choice is that indices vary mostly with
c     i>=j and k>=l with kl varying fastest ... improved cache hits?
c     If change this mapping must change overwrite logic below
c
            ij = im(i) + j      ! NOT im1 since go from 0
            kl = im(k) + l
            if (ij .eq. kl) gj = gj + gj
c
c     Using overwrite logic ... separates stores/loads, and enables
c     partial pipelining of */+ pairs ... at least in principle
c
            fij = fock(ij) + gj*dens(kl)
            fkl = fock(kl) + gj*dens(ij)
            if (i .eq. j) gk = gk + gk
            fock(ij) = fij
            fock(kl) = fkl
c
            if (k .eq. l) gk = gk + gk
c
            ik = im(i) + k
            il = im(i) + l
            jl = im(j) + l
            jk = im(j) + k
c
            fik = fock(ik) + gk*dens(jl)
            fil = fock(il) + gk*dens(jk)
            fjk = fock(jk) + gk*dens(il)
            fjl = fock(jl) + gk*dens(ik)
            fock(ik) = fik
            fock(il) = fil
            fock(jl) = fjl
            fock(jk) = fjk
c
c     Code without overwrite logic
c
c$$$            fock(ij) = fock(ij) + gj*dens(kl)
c$$$            fock(kl) = fock(kl) + gj*dens(ij)
c$$$            fock(ik) = fock(ik) + gk*dens(jl)
c$$$            fock(il) = fock(il) + gk*dens(jk)
c$$$            fock(jl) = fock(jl) + gk*dens(ik)
c$$$            fock(jk) = fock(jk) + gk*dens(il)
         end if
      end do
c     
      end
      subroutine fock_2e_rep_mod_2_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer labels(4,neri)
      double precision eri(neri)
      double precision fock(2,0:*), dens(2,0:*)
      double precision jfac(*), kfac(*)
c     
      integer i, j, k, l, ind
      integer ij, ik, jk, kl, il, jl
      double precision g, gk1, gj1, gj2, gk2
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
c
            ij = im(i) + j      ! NOT im1 since go from 0
            ik = im(i) + k
            il = im(i) + l
            jl = im(j) + l
            jk = im(j) + k
            kl = im(k) + l
c
            gj1 = g*jfac(1)
            gk1 = g*kfac(1)
            fock(1,ij) = fock(1,ij) + gj1*dens(1,kl)
            fock(1,kl) = fock(1,kl) + gj1*dens(1,ij)
            fock(1,ik) = fock(1,ik) + gk1*dens(1,jl)
            fock(1,il) = fock(1,il) + gk1*dens(1,jk)
            fock(1,jl) = fock(1,jl) + gk1*dens(1,ik)
            fock(1,jk) = fock(1,jk) + gk1*dens(1,il)
c
            gj2 = g*jfac(2)
            gk2 = g*kfac(2)
            fock(2,ij) = fock(2,ij) + gj2*dens(2,kl)
            fock(2,kl) = fock(2,kl) + gj2*dens(2,ij)
            fock(2,ik) = fock(2,ik) + gk2*dens(2,jl)
            fock(2,il) = fock(2,il) + gk2*dens(2,jk)
            fock(2,jl) = fock(2,jl) + gk2*dens(2,ik)
            fock(2,jk) = fock(2,jk) + gk2*dens(2,il)
         end if
      end do
c     
      end
      subroutine fock_2e_rep_mod_4_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.7 1996-05-14 17:59:57 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer labels(4,neri)
      double precision eri(neri)
      double precision fock(4,0:*), dens(4,0:*)
      double precision jfac(*), kfac(*)
c     
      integer i, j, k, l, ind
      integer ij, ik, jk, kl, il, jl
      double precision g, gk1, gj1, gj2, gk2, gj3, gk3, gj4, gk4
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
c
            ij = im(i) + j      ! NOT im1 since go from 0
            ik = im(i) + k
            il = im(i) + l
            jl = im(j) + l
            jk = im(j) + k
            kl = im(k) + l
c
            gj1 = g*jfac(1)
            gk1 = g*kfac(1)
            fock(1,ij) = fock(1,ij) + gj1*dens(1,kl)
            fock(1,kl) = fock(1,kl) + gj1*dens(1,ij)
            fock(1,ik) = fock(1,ik) + gk1*dens(1,jl)
            fock(1,il) = fock(1,il) + gk1*dens(1,jk)
            fock(1,jl) = fock(1,jl) + gk1*dens(1,ik)
            fock(1,jk) = fock(1,jk) + gk1*dens(1,il)
c
            gj2 = g*jfac(2)
            gk2 = g*kfac(2)
            fock(2,ij) = fock(2,ij) + gj2*dens(2,kl)
            fock(2,kl) = fock(2,kl) + gj2*dens(2,ij)
            fock(2,ik) = fock(2,ik) + gk2*dens(2,jl)
            fock(2,il) = fock(2,il) + gk2*dens(2,jk)
            fock(2,jl) = fock(2,jl) + gk2*dens(2,ik)
            fock(2,jk) = fock(2,jk) + gk2*dens(2,il)
c
            gj3 = g*jfac(3)
            gk3 = g*kfac(3)
            fock(3,ij) = fock(3,ij) + gj3*dens(3,kl)
            fock(3,kl) = fock(3,kl) + gj3*dens(3,ij)
            fock(3,ik) = fock(3,ik) + gk3*dens(3,jl)
            fock(3,il) = fock(3,il) + gk3*dens(3,jk)
            fock(3,jl) = fock(3,jl) + gk3*dens(3,ik)
            fock(3,jk) = fock(3,jk) + gk3*dens(3,il)
c
            gj4 = g*jfac(4)
            gk4 = g*kfac(4)
            fock(4,ij) = fock(4,ij) + gj4*dens(4,kl)
            fock(4,kl) = fock(4,kl) + gj4*dens(4,ij)
            fock(4,ik) = fock(4,ik) + gk4*dens(4,jl)
            fock(4,il) = fock(4,il) + gk4*dens(4,jk)
            fock(4,jl) = fock(4,jl) + gk4*dens(4,ik)
            fock(4,jk) = fock(4,jk) + gk4*dens(4,il)
         end if
      end do
c     
      end
