      subroutine fock_2e_b( basis, nfock, tol2e, q4_atom,
     $     iat, jat, kat, lat, 
     $     idim, jdim, kdim, ldim,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_b.F,v 1.2 1995-03-31 01:35:11 d3g681 Exp $
c     
      implicit none
#include "mafdecls.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cscfps.fh"
c     
c     
      integer basis
      integer nfock
      double precision tol2e
      double precision q4_atom
      integer iat, jat, kat, lat
      integer idim, jdim, kdim, ldim
      double precision dij(nfock,idim,jdim), dik(nfock,idim,kdim)
      double precision dli(nfock,ldim,idim), djk(nfock,jdim,kdim)
      double precision dlj(nfock,ldim,jdim), dlk(nfock,ldim,kdim)
      double precision fij(nfock,idim,jdim), fik(nfock,idim,kdim)
      double precision fli(nfock,ldim,idim), fjk(nfock,jdim,kdim)
      double precision flj(nfock,ldim,jdim), flk(nfock,ldim,kdim)
c     
      logical oij, okl, oikjl, status
      integer ish, jsh, ksh, lsh, mem2, max2e
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer iatbflo, iatbfhi, jatbflo, jatbfhi
      integer katbflo, katbfhi, latbflo, latbfhi
      integer ishbflo, ishbfhi, jshbflo, jshbfhi
      integer kshbflo, kshbfhi, lshbflo, lshbfhi
      integer jshtop, kshtop, lshtop
      integer l_g, l_scr, k_g, k_scr
      double precision sij, smax, q4_shell
#include "itri.fh"
c     
c     iat, jat, kat, lat externally loop in canonical order over
c     atoms.  Compute the appropriate canonical contributions
c     to the fock matrix
c     
c     scratch space needed for 2-electron routines
c     
c     l_scr ... integral scratch space
c     l_g   ... buffer for eris
c     
      call int_mem_2e4c(max2e, mem2)
      status = .true.
      status = status .and.
     $     ma_push_get(MT_DBL, max2e, 'fock_2e: buf', l_g, k_g)
      status = status .and.
     $     ma_push_get(MT_DBL, mem2, 'fock_2e: scr', l_scr, k_scr)
      if (.not. status) call errquit('rhf_fock_2e_b: not enuf memory',0)
c     
c     Permutation symmetry boolean
c     
      oij = iat .eq. jat
      okl = kat .eq. lat
      oikjl = iat.eq.kat .and. jat.eq.lat
c     
c     Shell and function ranges for each atom
c     
      status = .true.
      status = status .and. bas_ce2cnr(basis, iat, ishlo, ishhi)
      status = status .and. bas_ce2cnr(basis, jat, jshlo, jshhi)
      status = status .and. bas_ce2cnr(basis, kat, kshlo, kshhi)
      status = status .and. bas_ce2cnr(basis, lat, lshlo, lshhi)
      status = status .and. bas_ce2bfr(basis, iat, iatbflo, iatbfhi)
      status = status .and. bas_ce2bfr(basis, jat, jatbflo, jatbfhi)
      status = status .and. bas_ce2bfr(basis, kat, katbflo, katbfhi)
      status = status .and. bas_ce2bfr(basis, lat, latbflo, latbfhi)
      if (.not. status) call errquit
     $     ('rhf_fock_2e_b: bas_ce to cnr/bfr failed', 0)
c     
c     loop over shells on each atom
c     
      smax = schwarz_max()
      do ish = ishlo, ishhi
         if (.not. bas_cn2bfr(basis, ish, ishbflo, ishbfhi))
     $        call errquit('rhf_fock_2e_b: bas_cn2bfr', ish)
         jshtop = jshhi
         if (oij) jshtop = ish
         do jsh = jshlo, jshtop
            sij = schwarz_shell(ish,jsh)
            if (sij*smax.ge.tol2e) then
               if (.not. bas_cn2bfr(basis, jsh, jshbflo, jshbfhi))
     $              call errquit('rhf_fock_2e_b: bas_cn2bfr', jsh)
               kshtop = kshhi
               if (oikjl) kshtop = ish
               do ksh = kshlo, kshtop
                  if (.not. bas_cn2bfr(basis, ksh, kshbflo, kshbfhi))
     $                 call errquit('rhf_fock_2e_b: bas_cn2bfr', ksh)
                  lshtop = lshhi
                  if (okl) lshtop = ksh
                  if (oikjl .and. ksh.eq.ish) lshtop = jsh
                  do lsh = lshlo, lshtop
                     if (sij*schwarz_shell(ksh,lsh)
     $                    .ge. tol2e) then
                        if (.not. bas_cn2bfr(basis, lsh,
     $                       lshbflo, lshbfhi)) call errquit
     $                       ('rhf_fock_2e_b: bas_cn2bfr', lsh)
                        q4_shell = q4_atom
c     
c     Compute the shell block of integrals
c     
                        if (oscfps) call pstat_on(ps_int_2e)
                        call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                       mem2, dbl_mb(k_scr), max2e, dbl_mb(k_g))
                        if (oscfps) call pstat_off(ps_int_2e)
c     
c     Prescale ERIs
c     
                        if (oscfps) call pstat_on(ps_fock_add)
                        call fock_eri_scale( ishbflo, ishbfhi,
     $                       jshbflo, jshbfhi,
     $                       kshbflo, kshbfhi,
     $                       lshbflo, lshbfhi,
     $                       dbl_mb(k_g) )
c     
c     Add shell-block ERI into Fock components
c     
                        call fock_2e_c( nfock, tol2e, q4_shell,
     $                       dbl_mb(k_g),
     $                       ishbflo, ishbfhi, jshbflo, jshbfhi,
     $                       kshbflo, kshbfhi, lshbflo, lshbfhi,
     $                       iatbflo, iatbfhi, jatbflo, jatbfhi,
     $                       katbflo, katbfhi, latbflo, latbfhi,
     $                       dij, dik, dli, djk, dlj, dlk,
     $                       fij, fik, fli, fjk, flj, flk )
                        if (oscfps) call pstat_off(ps_fock_add)
                        iscreen(5) = iscreen(5) + 1
                     else
                        iscreen(4) = iscreen(4) + 1
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
c     
      status = ma_pop_stack(l_scr)
      status = status .and. ma_pop_stack(l_g)
      if (.not. status) call errquit
     $     ('rhf_fock_2e_b: pop stack failed?', 0)
c     
      end
