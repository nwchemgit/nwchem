      subroutine fock_2e_shell_to_labels(tol2e,
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     ilab, jlab, klab, llab,
     $     eri, scr,
     $     neri)
c
c     $Id$
c
c     Given a shell block of integrals from nwints, 
c     compress out the zeroes, constrain for canonical
c     indices and return labels
c
      implicit none
      double precision tol2e
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer ilab(*), jlab(*), klab(*), llab(*)
      double precision scr(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
      double precision eri(*)
      integer neri
c
      double precision g
      integer i, j, k, l
      integer jtop, ktop, ltop
      logical oij, okl, oikjl
c     
c     Assume identical index ranges ==> identical blocks
c     
      oij = ilo.eq.jlo
      okl = klo.eq.llo
      oikjl = (ilo.eq.klo).and.(jlo.eq.llo)
c
c     Move integrals into scratch space
c
      neri = (ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*(lhi-llo+1)
      call dcopy(neri, eri, 1, scr, 1)
c     
c     Compute labels and compress out zeroes
c
      neri = 0
      do i = ilo, ihi
         jtop = jhi
         if (oij) jtop = i
         do j = jlo, jtop
            ktop = khi
            if (oikjl) ktop = i
            do k = klo, ktop
               ltop = lhi
               if (okl) ltop = k
               if (oikjl .and. k.eq.i) ltop = j
               do l = llo, ltop
                  g = scr(l,k,j,i)
                  if (abs(g) .gt. tol2e) then
                     neri = neri + 1
                     eri(neri)  = g
                     ilab(neri) = i
                     jlab(neri) = j
                     klab(neri) = k
                     llab(neri) = l
*                     write(6,1) i, j, k, l, g
* 1                   format(4i5,f20.6)
                  endif
               enddo
            enddo
         enddo
      enddo
c
      end
