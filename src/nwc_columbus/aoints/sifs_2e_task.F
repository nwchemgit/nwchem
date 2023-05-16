cgk be sure int_acc_set is called properly.

      subroutine sifs_2e_task( geom, basis, oskel, iatlo, jatlo, 
     &     katlo, latlo, iathi, jathi, kathi, lathi, ijk_prev,
     $     tmp, tol2e, aoint2, ninfo, info, ifmt, ibvtyp, values, 
     &     ibitv, buffer, clabs, ilabs, jlabs, klabs, llabs, eris )
      implicit none
#include "errquit.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "cscfps.fh"
#include "sym.fh"
c
c
      integer geom
      integer basis
      logical oskel
      integer iatlo, jatlo, katlo, latlo
      integer iathi, jathi, kathi, lathi
      integer ijk_prev(3,2)
      double precision tmp(*)
      double precision tol2e
      integer ilabs(*), jlabs(*), klabs(*), llabs(*)
      double precision eris(*)
c      
c     Given a task which comprises computing the contributions
c     arising from a quartet of blocks of atoms, do everything necessary
c
      double precision sij, sijkl, q4, smax
      logical otest, otestij
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer iat, jat, kat, lat, lattop, kattop, jattop
      logical oij, okl, oikjl
      integer junk
      logical int2e_set_bf_range
      external int2e_set_bf_range
      integer clabs(4,*)
      integer aoint2, ninfo, ifmt, ibvtyp, ibitv
      integer info(ninfo)
      double precision values(*), buffer(*)
 
cgk debug
*     write(*,*)'gk: entered sifs_2e_task'
*     write(*,*)'gk: oskel = ', oskel
cgk end
c
c     With caching get density and flush necessary fock blocks
c     and return magnitude of density blocks for screening
c
      if (.not. bas_ce2bfr(basis, iatlo, ilo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, iathi, junk, ihi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
c
      if (.not. bas_ce2bfr(basis, jatlo, jlo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, jathi, junk, jhi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
c
      if (.not. bas_ce2bfr(basis, katlo, klo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, kathi, junk, khi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
c
      if (.not. bas_ce2bfr(basis, latlo, llo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, lathi, junk, lhi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0, BASIS_ERR)
c
      if (owritefile) then
         if (.not. int2e_set_bf_range(
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi))
     $        call errquit('should not happen', 0, UNKNOWN_ERR)
      end if
c     
      q4 = 1.0d0
      oij   = iatlo.eq.jatlo
      okl   = katlo.eq.latlo
      oikjl = (iatlo.eq.katlo) .and. (jatlo.eq.latlo)
c
      do iat = iatlo, iathi
         jattop = jathi
         if (oij) jattop = iat
         do jat = jatlo, jattop
            kattop = kathi
            if (oikjl) kattop = iat
            do kat = katlo, kattop
               lattop = lathi
               if (okl) lattop = kat
               if (oikjl .and. kat.eq.iat) lattop = jat
               do lat = latlo,lattop
                  if (oskel) otest = sym_atom_quartet
     $                 (geom, iat, jat, kat, lat, q4)
*                 if (otest) then
cgk debug
*      write(*,'(a,4i3)')'gk: iat,jat,kat,lat:  ',iat,jat,kat,lat
cgk end
                     call int_2e_sifs_b(basis, tol2e, q4, iat, jat, 
     &                    kat, lat, ilo, jlo, klo, llo, ihi, jhi, khi, 
     &                    lhi, aoint2, ninfo, info, ifmt, ibvtyp, 
     &                    values, ibitv, buffer, clabs,
     &                    ilabs, jlabs, klabs, llabs, eris)
*                 end if
               end do
            end do
         end do
      end do
cgk debug
*     write(*,*)'gk: leaving sifs_2e_task'
cgk end
c     
      end
