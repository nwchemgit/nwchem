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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
*      write(6,*) ' f2el ', neri, nfock
*      write(6,1) (eri(i),ilab(i),jlab(i),klab(i),llab(i),
*     $     i=1,min(neri,9))
*1     format(d16.8,4i5)
c
      if (nfock .eq. 1) then
         call fock_2e_rep_1_label(nfock, nbf, jfac, kfac, tol2e,
     $        neri, ilab, jlab, klab, llab, eri, dens, fock)
         return
      else if (nfock .eq. 2) then
         call fock_2e_rep_2_label(nfock, nbf, jfac, kfac, tol2e,
     $        neri, ilab, jlab, klab, llab, eri, dens, fock)
         return
      else if (nfock .eq. 3) then
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
      else if (nfock .eq. 3) then
         call fock_2e_rep_mod_3_label(nfock, nbf, jfac, kfac, tol2e,
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
#ifdef SOLARIS
*     Overwrite logic seems slower on UltraSparc
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            gj = g*jfac
            gk = g*kfac
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
            ij = im(i) + j
            kl = im(k) + l
            ik = im(i) + k
            il = im(i) + l
            jl = im(j) + l
            jk = im(j) + k
            fock(ij) = fock(ij) + gj*dens(kl)
            fock(kl) = fock(kl) + gj*dens(ij)
            fock(ik) = fock(ik) + gk*dens(jl)
            fock(il) = fock(il) + gk*dens(jk)
            fock(jl) = fock(jl) + gk*dens(ik)
            fock(jk) = fock(jk) + gk*dens(il)
         end if
      end do
#else
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            gj = g*jfac
            gk = g*kfac
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
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
            if (k .eq. l) gk = gk + gk
            ik = im(i) + k
            il = im(i) + l
            jl = im(j) + l
            jk = im(j) + k
            fik = fock(ik) + gk*dens(jl)
            fil = fock(il) + gk*dens(jk)
            fjk = fock(jk) + gk*dens(il)
            fjl = fock(jl) + gk*dens(ik)
            fock(ik) = fik
            fock(il) = fil
            fock(jl) = fjl
            fock(jk) = fjk
         end if
      end do
#endif
c     
      end
c$$$      subroutine fock_2e_rep_mod_1_label(nfock, nbf, jfac, kfac, tol2e,
c$$$     $     neri, labels, eri, dens, fock)
c$$$c     
c$$$c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
c$$$c     
c$$$      implicit none
c$$$#include "cfockmul.fh"
c$$$c     
c$$$      integer nfock, nbf
c$$$      double precision tol2e
c$$$      integer neri
c$$$      integer labels(4,neri)
c$$$      double precision eri(neri)
c$$$      double precision fock(0:*), dens(0:*)
c$$$      double precision jfac, kfac
c$$$c     
c$$$      integer i, j, k, l, ind
c$$$      integer ij, ik, jk, kl, il, jl
c$$$      double precision g, gk, gj, fij, fkl, fik, fil, fjl, fjk
c$$$c
c$$$      do ind = 1, neri
c$$$         g = eri(ind)
c$$$         if (abs(g) .gt. tol2e) then
c$$$            gj = g*jfac
c$$$            gk = g*kfac
c$$$            i = labels(1,ind)
c$$$            j = labels(2,ind)
c$$$            k = labels(3,ind)
c$$$            l = labels(4,ind)
c$$$*     debug
c$$$*            if (i.lt.0 .or. i.ge.nbf) call errquit(' i bad', i)
c$$$*            if (j.lt.0 .or. j.ge.nbf) call errquit(' j bad', j)
c$$$*            if (l.lt.0 .or. l.ge.nbf) call errquit(' l bad', l)
c$$$*            if (k.lt.0 .or. k.ge.nbf) call errquit(' k bad', k)
c$$$c
c$$$c     Logic behind this choice is that indices vary mostly with
c$$$c     i>=j and k>=l with kl varying fastest ... improved cache hits?
c$$$c     If change this mapping must change overwrite logic below
c$$$c
c$$$            ij = im(i) + j      ! NOT im1 since go from 0
c$$$            kl = im(k) + l
c$$$            if (ij .eq. kl) gj = gj + gj
c$$$c
c$$$c     Using overwrite logic ... separates stores/loads, and enables
c$$$c     partial pipelining of */+ pairs ... at least in principle
c$$$c
c$$$            fij = fock(ij) + gj*dens(kl)
c$$$            fkl = fock(kl) + gj*dens(ij)
c$$$            if (i .eq. j) gk = gk + gk
c$$$            fock(ij) = fij
c$$$            fock(kl) = fkl
c$$$c
c$$$            if (k .eq. l) gk = gk + gk
c$$$c
c$$$            ik = im(i) + k
c$$$            il = im(i) + l
c$$$            jl = im(j) + l
c$$$            jk = im(j) + k
c$$$
c$$$            fik = fock(ik) + gk*dens(jl)
c$$$            fil = fock(il) + gk*dens(jk)
c$$$            fjk = fock(jk) + gk*dens(il)
c$$$            fjl = fock(jl) + gk*dens(ik)
c$$$            fock(ik) = fik
c$$$            fock(il) = fil
c$$$            fock(jl) = fjl
c$$$            fock(jk) = fjk
c$$$c
c$$$c     Code without overwrite logic
c$$$c
c$$$c$$$            fock(ij) = fock(ij) + gj*dens(kl)
c$$$c$$$            fock(kl) = fock(kl) + gj*dens(ij)
c$$$c$$$            fock(ik) = fock(ik) + gk*dens(jl)
c$$$c$$$            fock(il) = fock(il) + gk*dens(jk)
c$$$c$$$            fock(jl) = fock(jl) + gk*dens(ik)
c$$$c$$$            fock(jk) = fock(jk) + gk*dens(il)
c$$$         end if
c$$$      end do
c$$$c     
c$$$      end
      subroutine fock_2e_rep_mod_2_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
      subroutine fock_2e_rep_mod_3_label(nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
c     
      implicit none
#include "cfockmul.fh"
c     
      integer nfock, nbf
      double precision tol2e
      integer neri
      integer labels(4,neri)
      double precision eri(neri)
      double precision fock(3,0:*), dens(3,0:*)
      double precision jfac(*), kfac(*)
c     
      integer i, j, k, l, ind
      integer ij, ik, jk, kl, il, jl
      double precision g, gk1, gj1, gj2, gk2, gj3, gk3
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
         end if
      end do
c     
      end
      subroutine fock_2e_rep_mod_4_label
     $     (nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
      double precision jfac(nfock), kfac(nfock)
c     
      integer i, j, k, l, ind
      integer ij, ik, jk, kl, il, jl
      double precision g, gj, gk, gk1, gj1, gj2, gk2, gj3, gk3, gj4, gk4
      double precision 
     $     fock1ij, fock1kl, fock2ij, fock2kl, fock3ij, fock3kl, 
     $     fock4ij, fock4kl, fock1ik, fock1il, fock1jl, fock1jk,
     $     fock2ik, fock2il, fock2jl, fock2jk, fock3ik, fock3il,
     $     fock3jl, fock3jk, fock4ik, fock4il, fock4jl, fock4jk
c     
      if (
     $     jfac(1).eq.1.0d0 .and. jfac(2).eq.0.0d0 .and. 
     $     jfac(3).eq.1.0d0 .and. jfac(4).eq.0.0d0 .and. 
     $     kfac(1).eq.0.0d0 .and. kfac(2).eq.0.5d0 .and. 
     $     kfac(3).eq.0.0d0 .and. kfac(4).eq.0.5d0 ) then ! UHF special
         do ind = 1, neri
            g = eri(ind)
            if (abs(g) .gt. tol2e) then
               i = labels(1,ind)
               j = labels(2,ind)
               k = labels(3,ind)
               l = labels(4,ind)
c     
               ij = (im(i) + j)*4
               kl = (im(k) + l)*4
c     
               gj = g
               if (ij .eq. kl) gj = gj + gj
c     
               fock1ij = fock(  ij) + gj*dens(  kl)
               fock3ij = fock(2+ij) + gj*dens(2+kl)
               fock1kl = fock(  kl) + gj*dens(  ij)
               fock3kl = fock(2+kl) + gj*dens(2+ij)
c     
               gk = g * 0.5d0
c
               fock(  ij) = fock1ij
               fock(2+ij) = fock3ij
               fock(  kl) = fock1kl
               fock(2+kl) = fock3kl
c     
               if (i .eq. j) gk = gk + gk
c     
               ik = (im(i) + k)*4 + 1
               jl = (im(j) + l)*4 + 1
c
               if (k .eq. l) gk = gk + gk
c
               il = (im(i) + l)*4 + 1
               jk = (im(j) + k)*4 + 1
c     
               fock2ik = fock(  ik) + gk*dens(  jl)
               fock4ik = fock(2+ik) + gk*dens(2+jl)
               fock2il = fock(  il) + gk*dens(  jk)
               fock4il = fock(2+il) + gk*dens(2+jk)
               fock2jl = fock(  jl) + gk*dens(  ik)
               fock4jl = fock(2+jl) + gk*dens(2+ik)
               fock2jk = fock(  jk) + gk*dens(  il)
               fock4jk = fock(2+jk) + gk*dens(2+il)
c     
               fock(  ik) = fock2ik
               fock(2+ik) = fock4ik
               fock(  il) = fock2il
               fock(2+il) = fock4il
               fock(  jl) = fock2jl
               fock(2+jl) = fock4jl
               fock(  jk) = fock2jk
               fock(2+jk) = fock4jk
c     
            end if
         end do

      else

         do ind = 1, neri
            g = eri(ind)
            if (abs(g) .gt. tol2e) then
               i = labels(1,ind)
               j = labels(2,ind)
               k = labels(3,ind)
               l = labels(4,ind)
c     
               ij = (im(i) + j)*4 
               kl = (im(k) + l)*4
c     
               gj = g
               if (ij .eq. kl) gj = gj + gj
               gj1 = gj*jfac(1)
               gj2 = gj*jfac(2)
               gj3 = gj*jfac(3)
               gj4 = gj*jfac(4)
c     
               fock1ij = fock(  ij) + gj1*dens(  kl)
               fock2ij = fock(1+ij) + gj2*dens(1+kl)
               fock3ij = fock(2+ij) + gj3*dens(2+kl)
               fock4ij = fock(3+ij) + gj4*dens(3+kl)
               fock1kl = fock(  kl) + gj1*dens(  ij)
               fock2kl = fock(1+kl) + gj2*dens(1+ij)
               fock3kl = fock(2+kl) + gj3*dens(2+ij)
               fock4kl = fock(3+kl) + gj4*dens(3+ij)
c     
               fock(  ij) = fock1ij
               fock(1+ij) = fock2ij
               fock(2+ij) = fock3ij
               fock(3+ij) = fock4ij
               fock(  kl) = fock1kl
               fock(1+kl) = fock2kl
               fock(2+kl) = fock3kl
               fock(3+kl) = fock4kl
c     
c     
               gk = g
               if (i .eq. j) gk = gk + gk
               if (k .eq. l) gk = gk + gk
               gk1 = gk*kfac(1)
               gk2 = gk*kfac(2)
               gk3 = gk*kfac(3)
               gk4 = gk*kfac(4)
c     
               ik = (im(i) + k)*4
               il = (im(i) + l)*4
               jl = (im(j) + l)*4
               jk = (im(j) + k)*4
c     
               fock1ik = fock(  ik) + gk1*dens(  jl)
               fock2ik = fock(1+ik) + gk2*dens(1+jl)
               fock3ik = fock(2+ik) + gk3*dens(2+jl)
               fock4ik = fock(3+ik) + gk4*dens(3+jl)
               fock1il = fock(  il) + gk1*dens(  jk)
               fock2il = fock(1+il) + gk2*dens(1+jk)
               fock3il = fock(2+il) + gk3*dens(2+jk)
               fock4il = fock(3+il) + gk4*dens(3+jk)
               fock1jl = fock(  jl) + gk1*dens(  ik)
               fock2jl = fock(1+jl) + gk2*dens(1+ik)
               fock3jl = fock(2+jl) + gk3*dens(2+ik)
               fock4jl = fock(3+jl) + gk4*dens(3+ik)
               fock1jk = fock(  jk) + gk1*dens(  il)
               fock2jk = fock(1+jk) + gk2*dens(1+il)
               fock3jk = fock(2+jk) + gk3*dens(2+il)
               fock4jk = fock(3+jk) + gk4*dens(3+il)
c     
               fock(  ik) = fock1ik
               fock(1+ik) = fock2ik
               fock(2+ik) = fock3ik
               fock(3+ik) = fock4ik
               fock(  il) = fock1il
               fock(1+il) = fock2il
               fock(2+il) = fock3il
               fock(3+il) = fock4il
               fock(  jl) = fock1jl
               fock(1+jl) = fock2jl
               fock(2+jl) = fock3jl
               fock(3+jl) = fock4jl
               fock(  jk) = fock1jk
               fock(1+jk) = fock2jk
               fock(2+jk) = fock3jk
               fock(3+jk) = fock4jk
c     
            end if
         end do
      endif
c     
      end
      subroutine fock_2e_rep_mod_4_label_new_ok
     $     (nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
      double precision g, g2, gk1, gj1, gj2, gk2, gj3, gk3, gj4, gk4
      double precision 
     $     fock1ij, fock1kl, fock2ij, fock2kl, fock3ij, fock3kl, 
     $     fock4ij, fock4kl, fock1ik, fock1il, fock1jl, fock1jk,
     $     fock2ik, fock2il, fock2jl, fock2jk, fock3ik, fock3il,
     $     fock3jl, fock3jk, fock4ik, fock4il, fock4jl, fock4jk
c
      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            i = labels(1,ind)
            j = labels(2,ind)
            k = labels(3,ind)
            l = labels(4,ind)
c
c
            ij = im(i) + j      ! NOT im1 since go from 0
            kl = im(k) + l
c
            if (ij .eq. kl) then
               g2 = g + g
               gj1 = g2*jfac(1)
               gj2 = g2*jfac(2)
               gj3 = g2*jfac(3)
               gj4 = g2*jfac(4)
            else
               gj1 = g*jfac(1)
               gj2 = g*jfac(2)
               gj3 = g*jfac(3)
               gj4 = g*jfac(4)
            endif
c
            fock1ij = fock(1,ij) + gj1*dens(1,kl)
            fock2ij = fock(2,ij) + gj2*dens(2,kl)
            fock3ij = fock(3,ij) + gj3*dens(3,kl)
            fock4ij = fock(4,ij) + gj4*dens(4,kl)
            fock1kl = fock(1,kl) + gj1*dens(1,ij)
            fock2kl = fock(2,kl) + gj2*dens(2,ij)
            fock3kl = fock(3,kl) + gj3*dens(3,ij)
            fock4kl = fock(4,kl) + gj4*dens(4,ij)
c
            fock(1,ij) = fock1ij
            fock(2,ij) = fock2ij
            fock(3,ij) = fock3ij
            fock(4,ij) = fock4ij
            fock(1,kl) = fock1kl
            fock(2,kl) = fock2kl
            fock(3,kl) = fock3kl
            fock(4,kl) = fock4kl
c
c     Working non-overwrite code for exchange
c
c$$$            ik = im(i) + k
c$$$            il = im(i) + l
c$$$            jl = im(j) + l
c$$$            jk = im(j) + k
c$$$c
c$$$            gk1 = g*kfac(1)
c$$$            fock(1,ik) = fock(1,ik) + gk1*dens(1,jl)
c$$$            fock(1,il) = fock(1,il) + gk1*dens(1,jk)
c$$$            fock(1,jl) = fock(1,jl) + gk1*dens(1,ik)
c$$$            fock(1,jk) = fock(1,jk) + gk1*dens(1,il)
c$$$c
c$$$            gk2 = g*kfac(2)
c$$$            fock(2,ik) = fock(2,ik) + gk2*dens(2,jl)
c$$$            fock(2,il) = fock(2,il) + gk2*dens(2,jk)
c$$$            fock(2,jl) = fock(2,jl) + gk2*dens(2,ik)
c$$$            fock(2,jk) = fock(2,jk) + gk2*dens(2,il)
c$$$c
c$$$            gk3 = g*kfac(3)
c$$$            fock(3,ik) = fock(3,ik) + gk3*dens(3,jl)
c$$$            fock(3,il) = fock(3,il) + gk3*dens(3,jk)
c$$$            fock(3,jl) = fock(3,jl) + gk3*dens(3,ik)
c$$$            fock(3,jk) = fock(3,jk) + gk3*dens(3,il)
c$$$c
c$$$            gk4 = g*kfac(4)
c$$$            fock(4,ik) = fock(4,ik) + gk4*dens(4,jl)
c$$$            fock(4,il) = fock(4,il) + gk4*dens(4,jk)
c$$$            fock(4,jl) = fock(4,jl) + gk4*dens(4,ik)
c$$$            fock(4,jk) = fock(4,jk) + gk4*dens(4,il)
c
c
c
            if (i .eq. j) g = g + g
            if (k .eq. l) then
               g2 = g + g
               gk1 = g2*kfac(1)
               gk2 = g2*kfac(2)
               gk3 = g2*kfac(3)
               gk4 = g2*kfac(4)
            else
               gk1 = g*kfac(1)
               gk2 = g*kfac(2)
               gk3 = g*kfac(3)
               gk4 = g*kfac(4)
            endif
c
            ik = im(i) + k
            il = im(i) + l
            jl = im(j) + l
            jk = im(j) + k
c
            fock1ik = fock(1,ik) + gk1*dens(1,jl)
            fock2ik = fock(2,ik) + gk2*dens(2,jl)
            fock3ik = fock(3,ik) + gk3*dens(3,jl)
            fock4ik = fock(4,ik) + gk4*dens(4,jl)
            fock1il = fock(1,il) + gk1*dens(1,jk)
            fock2il = fock(2,il) + gk2*dens(2,jk)
            fock3il = fock(3,il) + gk3*dens(3,jk)
            fock4il = fock(4,il) + gk4*dens(4,jk)
            fock1jl = fock(1,jl) + gk1*dens(1,ik)
            fock2jl = fock(2,jl) + gk2*dens(2,ik)
            fock3jl = fock(3,jl) + gk3*dens(3,ik)
            fock4jl = fock(4,jl) + gk4*dens(4,ik)
            fock1jk = fock(1,jk) + gk1*dens(1,il)
            fock2jk = fock(2,jk) + gk2*dens(2,il)
            fock3jk = fock(3,jk) + gk3*dens(3,il)
            fock4jk = fock(4,jk) + gk4*dens(4,il)
c
            fock(1,ik) = fock1ik
            fock(2,ik) = fock2ik
            fock(3,ik) = fock3ik
            fock(4,ik) = fock4ik
            fock(1,il) = fock1il
            fock(2,il) = fock2il
            fock(3,il) = fock3il
            fock(4,il) = fock4il
            fock(1,jl) = fock1jl
            fock(2,jl) = fock2jl
            fock(3,jl) = fock3jl
            fock(4,jl) = fock4jl
            fock(1,jk) = fock1jk
            fock(2,jk) = fock2jk
            fock(3,jk) = fock3jk
            fock(4,jk) = fock4jk
c
         end if
      end do
c     
      end
      subroutine fock_2e_rep_mod_4_label_old
     $     (nfock, nbf, jfac, kfac, tol2e,
     $     neri, labels, eri, dens, fock)
c     
c     $Id: fock_2e_lab.F,v 1.11 1999-08-07 23:33:08 d3g681 Exp $
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
