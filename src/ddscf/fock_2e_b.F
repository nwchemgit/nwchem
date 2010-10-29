      subroutine fock_2e_b( basis, nfock, sijkl, tol2e,
     $     dentol, q4,
     $     iat, jat, kat, lat, 
     $     ilo, jlo, klo, llo,
     $     ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id$
c
c     Called by fock_2e_task once density matrix blocks have 
c     been fetched.
c     
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "util.fh"
#include "cscfps.fh"
#include "stdio.fh"
c
c     WARNING ... this piece of code will soon retire forever.
c     
      integer basis
      integer nfock
      double precision tol2e, dentol
      double precision sijkl
      double precision q4
      integer iat, jat, kat, lat
      integer ilo, jlo, klo, llo ! Range for which have density/fock
      integer ihi, jhi, khi, lhi
      double precision dij(nfock,*), dik(nfock,*)
      double precision dli(nfock,*), djk(nfock,*)
      double precision dlj(nfock,*), dlk(nfock,*)
      double precision fij(nfock,*), fik(nfock,*)
      double precision fli(nfock,*), fjk(nfock,*)
      double precision flj(nfock,*), flk(nfock,*)
c     
c
      double precision ncomputed
      common/ctemp/ncomputed
c
      logical oij, okl, oikjl, status
      integer ish, jsh, ksh, lsh, mem2, max2e
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer ibflo, jbflo, kbflo, lbflo, ibfhi, jbfhi, kbfhi, lbfhi
      integer ishbflo, ishbfhi, jshbflo, jshbfhi
      integer kshbflo, kshbfhi, lshbflo, lshbfhi
      integer jshtop, kshtop, lshtop
      integer l_g, l_scr, k_g, k_scr
      integer l_i, l_j, l_k, l_l
      integer k_i, k_j, k_k, k_l
c
      double precision sij, smax, denmax, dtol2e
      integer neri
      logical int2e_file_write
      external int2e_file_write
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
c     Shell and function ranges for each atom
c     
      status = .true.
      status = status .and. bas_ce2cnr(basis, iat, ishlo, ishhi)
      status = status .and. bas_ce2cnr(basis, jat, jshlo, jshhi)
      status = status .and. bas_ce2cnr(basis, kat, kshlo, kshhi)
      status = status .and. bas_ce2cnr(basis, lat, lshlo, lshhi)
      if (.not. owritefile) then
         status = status .and. bas_ce2bfr(basis, iat, ibflo, ibfhi)
         status = status .and. bas_ce2bfr(basis, jat, jbflo, jbfhi)
         status = status .and. bas_ce2bfr(basis, kat, kbflo, kbfhi)
         status = status .and. bas_ce2bfr(basis, lat, lbflo, lbfhi)
      endif
      if (.not. status) call errquit
     $     ('fock_2e_b: bas_ce to cnr/bfr failed', 0, BASIS_ERR)
c
      if (owritefile) then
         dtol2e = min(tol2e,dentol)
      else
c     
c     Density screening
c
         call fock_density_screen(nfock,
     $        ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $        ibflo, jbflo, kbflo, lbflo, ibfhi, jbfhi, kbfhi, lbfhi,
     $        dij, dik, dli, djk, dlj, dlk, denmax)
         if (sijkl*denmax .lt. tol2e) return
         dtol2e = dentol
      endif
c
      iscreen(3) = iscreen(3) + 1
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
      if (.not. status) call errquit('fock_2e_b: not enuf memory',0,
     &       MA_ERR)
c     
c     Permutation symmetry boolean
c     
      oij = iat .eq. jat
      okl = kat .eq. lat
      oikjl = iat.eq.kat .and. jat.eq.lat
c
c     loop over shells on each atom
c     
      smax = schwarz_max()
      do ish = ishlo, ishhi
         if (.not. bas_cn2bfr(basis, ish, ishbflo, ishbfhi))
     $        call errquit('fock_2e_b: bas_cn2bfr', ish, BASIS_ERR)
         jshtop = jshhi
         if (oij) jshtop = ish
         do jsh = jshlo, jshtop
            sij = schwarz_shell(ish,jsh)
            if (sij*smax.ge.dtol2e) then
               if (.not. bas_cn2bfr(basis, jsh, jshbflo, jshbfhi))
     $            call errquit('fock_2e_b: bas_cn2bfr',jsh,
     $                         BASIS_ERR)
               kshtop = kshhi
               if (oikjl) kshtop = ish
               do ksh = kshlo, kshtop
                  if (.not. bas_cn2bfr(basis, ksh, kshbflo, kshbfhi))
     $               call errquit('fock_2e_b: bas_cn2bfr',ksh,
     $                            BASIS_ERR)
                  lshtop = lshhi
                  if (okl) lshtop = ksh
                  if (oikjl .and. ksh.eq.ish) lshtop = jsh
                  do lsh = lshlo, lshtop
                     if (sij*schwarz_shell(ksh,lsh)
     $                    .ge. dtol2e) then
                        if (.not.bas_cn2bfr(basis,lsh,lshbflo,lshbfhi)) 
     $                     call errquit('fock_2e_b: bas_cn2bfr',lsh,
     $                                  BASIS_ERR)
c     
c     Compute the shell block of integrals with labels
c     
                        call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                       mem2, dbl_mb(k_scr), max2e, dbl_mb(k_g))
                        if (oerilabel) then
c
c     New code with labels
c        - Prescale ERIs
c        - Add shell-block ERI into Fock components
c
                          call fock_2e_shell_to_labels( dtol2e,
     $                            ishbflo, ishbfhi, jshbflo, jshbfhi,
     $                            kshbflo, kshbfhi, lshbflo, lshbfhi,
     $                            int_mb(k_i), int_mb(k_j),
     $                            int_mb(k_k), int_mb(k_l),
     $                            dbl_mb(k_g), dbl_mb(k_scr),  neri)
       ncomputed = ncomputed + neri
                          call fock_eri_label_scale(neri, q4,
     $                                    int_mb(k_i), int_mb(k_j),
     $                                    int_mb(k_k), int_mb(k_l),
     $                                    dbl_mb(k_g))
                          call fock_2e_label( nfock, dtol2e, neri,
     $                         int_mb(k_i), int_mb(k_j),
     $                         int_mb(k_k), int_mb(k_l), dbl_mb(k_g),
     $                         ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $                         dij, dik, dli, djk, dlj, dlk,
     $                         fij, fik, fli, fjk, flj, flk )
c
c     Cache integrals if desired
c
                          if (owritefile) then
                            if (.not. int2e_file_write(neri,
     $                         int_mb(k_i), int_mb(k_j),
     $                         int_mb(k_k), int_mb(k_l), dbl_mb(k_g)))
     .                            then
                  write(luout,*) 
                  write(luout,*) ' CHANGE TO DIRECT SCF '
                  write(luout,*) 
                  call util_flush(luout)
                               call errquit(
     .        ' fock_2e_b: write failed with neri=',neri,
     &       UNKNOWN_ERR)
                            endif
                        if (util_print('int2efile write',print_never))
     $                          call print_integ_list(neri, 
     $                          int_mb(k_i), int_mb(k_j),
     $                          int_mb(k_k), int_mb(k_l), dbl_mb(k_g))
                          end if
                        else
c
c    Old code with no labels
c        - Prescale ERIs
c        - Add shell-block ERI into Fock components 
c          (special code for 1 fock matrix)
c
                          call fock_eri_scale( ishbflo, ishbfhi,
     $                                         jshbflo, jshbfhi,
     $                                         kshbflo, kshbfhi,
     $                                         lshbflo, lshbfhi,
     $                                         dbl_mb(k_g) )
                          if (nfock.eq.1) then
                            call fock_2e_c1( dtol2e, q4, dbl_mb(k_g),
     $                            ishbflo, ishbfhi, jshbflo, jshbfhi,
     $                            kshbflo, kshbfhi, lshbflo, lshbfhi,
     $                            ilo, ihi, jlo, jhi,
     $                            klo, khi, llo, lhi,
     $                            dij, dik, dli, djk, dlj, dlk,
     $                            fij, fik, fli, fjk, flj, flk )
                          else
                            call fock_2e_c( nfock, dtol2e, q4,
     $                            dbl_mb(k_g),
     $                            ishbflo, ishbfhi, jshbflo, jshbfhi,
     $                            kshbflo, kshbfhi, lshbflo, lshbfhi,
     $                            ilo, ihi, jlo, jhi,
     $                            klo, khi, llo, lhi,
     $                            dij, dik, dli, djk, dlj, dlk,
     $                            fij, fik, fli, fjk, flj, flk )
                          endif
                        endif
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
     $     ('fock_2e_b: pop stack failed?', 0, MA_ERR)
c     
      end
