      subroutine fock_2e_b( basis, nfock, tol2e, q4_atom,
     $     iat, jat, kat, lat, 
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_b.F,v 1.3 1995-05-22 06:56:11 d3g681 Exp $
c     
      implicit none
#include "mafdecls.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "cfock.fh"
#include "util.fh"
c     
c     
      integer basis
      integer nfock
      double precision tol2e
      double precision q4_atom
      integer iat, jat, kat, lat
      double precision dij(nfock,*), dik(nfock,*)
      double precision dli(nfock,*), djk(nfock,*)
      double precision dlj(nfock,*), dlk(nfock,*)
      double precision fij(nfock,*), fik(nfock,*)
      double precision fli(nfock,*), fjk(nfock,*)
      double precision flj(nfock,*), flk(nfock,*)
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
      integer l_i, l_j, l_k, l_l
      integer k_i, k_j, k_k, k_l
      double precision sij, smax, q4_shell
      integer neri
      logical int2e_file_write, int2e_set_bf_range
      external int2e_file_write, int2e_set_bf_range
c     
c     iat, jat, kat, lat externally loop in canonical order over
c     atoms.  Compute the appropriate canonical contributions
c     to the fock matrix
c     
c     scratch space needed for 2-electron routines
c     
c     l_scr ... integral scratch space
c     l_g   ... buffer for eris
c     l_i, j, k, l ... integral labels
c     
      call int_mem_2e4c(max2e, mem2)
      mem2 = max(mem2,max2e)
      status = .true.
      status = status .and.
     $     ma_push_get(MT_DBL, max2e, 'fock_2e: buf', l_g, k_g)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'fock_2e: buf', l_i, k_i)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'fock_2e: buf', l_j, k_j)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'fock_2e: buf', l_k, k_k)
      status = status .and.
     $     ma_push_get(MT_INT, max2e, 'fock_2e: buf', l_l, k_l)
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
      if (owritefile) then
         if (.not. int2e_set_bf_range(
     $        iatbflo, iatbfhi, jatbflo, jatbfhi,
     $        katbflo, katbfhi, latbflo, latbfhi))
     $        call errquit('should not ahppen', 0)
      end if
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
#ifdef NOLABELS
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
#else
c
c     New code with labels
c
c     
c     Compute the shell block of integrals with labels
c     
                        if (oscfps) call pstat_on(ps_int_2e)
                        call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                       mem2, dbl_mb(k_scr), max2e, dbl_mb(k_g))
                        call fock_2e_shell_to_labels(tol2e,
     $                       ishbflo, ishbfhi, jshbflo, jshbfhi,
     $                       kshbflo, kshbfhi, lshbflo, lshbfhi,
     $                       int_mb(k_i), int_mb(k_j),
     $                       int_mb(k_k), int_mb(k_l), dbl_mb(k_g),
     $                       dbl_mb(k_scr),  neri)
                        if (oscfps) call pstat_off(ps_int_2e)
c     
c     Prescale ERIs
c     
                        if (oscfps) call pstat_on(ps_fock_add)
                        call fock_eri_label_scale(neri, q4_shell,
     $                       int_mb(k_i), int_mb(k_j),
     $                       int_mb(k_k), int_mb(k_l), dbl_mb(k_g))
c
c     Cache integrals if desired
c
                        if (owritefile) then
                           if (.not. int2e_file_write(neri,
     $                       int_mb(k_i), int_mb(k_j),
     $                       int_mb(k_k), int_mb(k_l), dbl_mb(k_g)))
     $                          call errquit('Arrgggghh',0)
                           if(util_print('int2efile write',print_never))
     $                          call print_integ_list(neri, 
     $                          int_mb(k_i), int_mb(k_j),
     $                          int_mb(k_k), int_mb(k_l), dbl_mb(k_g))
                        end if
c     
c     Add shell-block ERI into Fock components
c     
                        call fock_2e_label( nfock, tol2e,
     $                       neri,
     $                       int_mb(k_i), int_mb(k_j),
     $                       int_mb(k_k), int_mb(k_l), dbl_mb(k_g),
     $                       iatbflo, iatbfhi, jatbflo, jatbfhi,
     $                       katbflo, katbfhi, latbflo, latbfhi,
     $                       dij, dik, dli, djk, dlj, dlk,
     $                       fij, fik, fli, fjk, flj, flk )
                        if (oscfps) call pstat_off(ps_fock_add)
#endif

                        iscreen(5) = iscreen(5) + 1
                     else
                        iscreen(4) = iscreen(4) + 1
                     end if
                  end do
               end do
            end if
         end do
      end do
c     
      status = ma_pop_stack(l_scr)
      status = status .and. ma_pop_stack(l_l)
      status = status .and. ma_pop_stack(l_k)
      status = status .and. ma_pop_stack(l_j)
      status = status .and. ma_pop_stack(l_i)
      status = status .and. ma_pop_stack(l_g)
      if (.not. status) call errquit
     $     ('rhf_fock_2e_b: pop stack failed?', 0)
c     
      end
