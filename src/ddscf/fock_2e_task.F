      subroutine fock_2e_task(
     $     geom, basis, oskel,
     $     iatlo, jatlo, katlo, latlo,
     $     iathi, jathi, kathi, lathi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, tol2e, dentol)
      implicit none
#include "errquit.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "cscfps.fh"
#include "sym.fh"
c
c     $Id$
c
c     Given an block of atomic quartets, fetch the necessary blocks
c     of the density matrices, call fock_2e_b to add in
c     in the fock matrix contribution and then accumulate
c     the contributions.
c
      integer geom
      integer basis
      logical oskel
      integer iatlo, jatlo, katlo, latlo
      integer iathi, jathi, kathi, lathi
      integer ijk_prev(3,2)
      integer nfock
      integer vg_dens(nfock), vg_fock(nfock)
      double precision jfac(nfock), kfac(nfock)
      double precision dij(nfock,*), dik(nfock,*), dli(nfock,*),
     $     djk(nfock,*), dlj(nfock,*), dlk(nfock,*)
      double precision fij(nfock,*), fik(nfock,*), fli(nfock,*),
     $     fjk(nfock,*), flj(nfock,*), flk(nfock,*)
      double precision tmp(*)
      double precision tol2e, dentol
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
      call fock_init_cmul(ihi-ilo+1,jhi-jlo+1,lhi-llo+1)
c
      call fock_2e_cache_dens_fock(
     $     ilo, jlo, klo, llo,
     $     ihi, jhi, khi, lhi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp)
c
      q4 = 1.0d0
      oij   = iatlo.eq.jatlo
      okl   = katlo.eq.latlo
      oikjl = (iatlo.eq.katlo) .and. (jatlo.eq.latlo)
c
      if (oscfps) call pstat_on(ps_fock_add)
      if (onoblock) then
         smax = schwarz_max()
         do iat = iatlo, iathi
            jattop = jathi
            if (oij) jattop = iat
            do jat = jatlo, jattop
               sij = schwarz_atom(iat,jat)
               otestij = (sij*smax) .gt. dentol
               if (otestij) then
                  kattop = kathi
                  if (oikjl) kattop = iat
                  do kat = katlo, kattop
                     lattop = lathi
                     if (okl) lattop = kat
                     if (oikjl .and. kat.eq.iat) lattop = jat
                     do lat = latlo,lattop
                        sijkl = sij*schwarz_atom(kat,lat)
                        otest = sijkl .gt. dentol
                        if (otest .and. oskel) otest = sym_atom_quartet
     $                       (geom, iat, jat, kat, lat, q4)
                        if (otest) then
                           call fock_2e_b(basis, nfock, sijkl, tol2e,
     $                          dentol, q4, iat, jat, kat, lat, 
     $                          ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $                          dij, dik, dli, djk, dlj, dlk,
     $                          fij, fik, fli, fjk, flj, flk)
                        end if
                     end do
                  end do
               end if
            end do
         end do
      else
         call fock_2e_texas(
     $        geom, basis, oskel,
     $        iatlo, jatlo, katlo, latlo,
     $        iathi, jathi, kathi, lathi,
     $        ijk_prev,
     $        nfock, vg_dens, vg_fock,
     $        jfac, kfac,
     $        ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $        dij, dik, dli, djk, dlj, dlk,
     $        fij, fik, fli, fjk, flj, flk,
     $        tmp, tol2e, dentol, oij, okl, oikjl)
      end if
      if (oscfps) call pstat_off(ps_fock_add)
c     
c     Update F blocks
c     
      if (oscfps) call pstat_on(ps_gaops)
      call fock_upd_blk(nfock, vg_fock,
     $     llo, lhi, ilo, ihi, kfac, fli, tmp)
      call fock_upd_blk(nfock, vg_fock,
     $     llo, lhi, jlo, jhi, kfac, flj, tmp)
      call fock_upd_blk(nfock, vg_fock,
     $     llo, lhi, klo, khi, jfac, flk, tmp)
      if (oscfps) call pstat_off(ps_gaops)
c     
      end
      subroutine fock_2e_texas(
     $     geom, basis, oskel,
     $     iatlo, jatlo, katlo, latlo,
     $     iathi, jathi, kathi, lathi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, tol2e, dentol, oij, okl, oikjl)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      integer geom
      integer basis
      logical oskel
      integer iatlo, jatlo, katlo, latlo
      integer iathi, jathi, kathi, lathi
      integer ijk_prev(3,2)
      integer nfock
      integer ilo, jlo, klo, llo, ihi, jhi, khi, lhi
c
      integer vg_dens(nfock), vg_fock(nfock)
      double precision jfac(nfock), kfac(nfock)
      double precision dij(nfock,*), dik(nfock,*), dli(nfock,*),
     $     djk(nfock,*), dlj(nfock,*), dlk(nfock,*)
      double precision fij(nfock,*), fik(nfock,*), fli(nfock,*),
     $     fjk(nfock,*), flj(nfock,*), flk(nfock,*)
      double precision tmp(*)
      double precision tol2e, dentol
      double precision smax, sijkl
      logical oij, okl, oikjl
      integer iat, jat, kat, lat, nsh
c
c     These for sorting shells within each atom pair block
c
      integer max_pair_t, max_sh, max_at
      parameter (max_at = 10)
      parameter (max_pair_t = 100000, max_sh = 64) ! Should dynamically alloc
      double precision skl(max_sh*max_sh) ! Should be just max shell block
      double precision sij(max_sh*max_sh) ! Should be just max shell block
      integer ijlist(2,max_sh*max_sh), ij
      integer kllist(2,max_sh*max_sh), kl
      integer nt_ij
      integer nt_kl
      integer k_ijinfo,l_ijinfo,k_klinfo,l_klinfo
c
      integer max_q, max_p
      parameter (max_q = 10000)
      integer ijkl(max_q,4)
      double precision q4(max_q)
      double precision rdij(max_at,max_at), rdik(max_at,max_at), 
     $     rdli(max_at,max_at),
     $     rdjk(max_at,max_at), rdlj(max_at,max_at), rdlk(max_at,max_at)
c
      integer ncent, ncent_unique
      integer i
      integer ia,ja,ka,la
      integer ijt, klt
      integer ish, jsh, ksh, lsh
      integer nq
c
c     These need to be refined
c
      double precision q4_at(max_at,max_at,max_at,max_at)
c
      integer maxcn
      parameter (maxcn = nw_max_shells)
      integer cn_to_ce(maxcn)
c
      logical status
      double precision q2, qq4, sijat, sklat
c
c     Form q4 factors over atoms and mapping from shells to atoms
c     ... very crude now for debug purposes
c
      if ((iathi-iatlo+1).gt.max_at) call errquit('too many ats',0,
     &       UNKNOWN_ERR)
      if ((jathi-jatlo+1).gt.max_at) call errquit('too many ats',0,
     &       UNKNOWN_ERR)
      if ((kathi-katlo+1).gt.max_at) call errquit('too many ats',0,
     &       UNKNOWN_ERR)
      if ((lathi-latlo+1).gt.max_at) call errquit('too many ats',0,
     &       UNKNOWN_ERR)
c
      smax = schwarz_max()
c
      call fock_den_sum(basis, nfock, iatlo, iathi, jatlo, jathi,
     $     dij, ilo, ihi, jlo, jhi, rdij, max_at, dentol/smax)
      call fock_den_sum(basis, nfock, iatlo, iathi, katlo, kathi,
     $     dik, ilo, ihi, klo, khi, rdik, max_at, dentol/smax)
      call fock_den_sum(basis, nfock, latlo, lathi, iatlo, iathi,
     $     dli, llo, lhi, ilo, ihi, rdli, max_at, dentol/smax)
      call fock_den_sum(basis, nfock, jatlo, jathi, katlo, kathi,
     $     djk, jlo, jhi, klo, khi, rdjk, max_at, dentol/smax)
      call fock_den_sum(basis, nfock, latlo, lathi, jatlo, jathi,
     $     dlj, llo, lhi, jlo, jhi, rdlj, max_at, dentol/smax)
      call fock_den_sum(basis, nfock, latlo, lathi, katlo, kathi,
     $     dlk, llo, lhi, klo, khi, rdlk, max_at, dentol/smax)
c
      if (oskel) then
         do iat = iatlo, iathi
            ia = iat-iatlo+1
            do jat = jatlo, jathi
               ja = jat-jatlo+1
               sijat = schwarz_atom(iat,jat)
               do kat = katlo, kathi
                  ka = kat-katlo+1
                  do lat = latlo, lathi
                     la = lat-latlo+1
                     sklat = schwarz_atom(kat,lat)
                     if (sklat*sijat .gt. dentol) then
                        status = sym_atom_quartet(geom, 
     $                       iat, jat, kat, lat, q4_at(ia,ja,ka,la))
                     else
                        q4_at(ia,ja,ka,la) = 0.0d0
                     end if
                  end do
               end do
            end do
         end do
      end if
c
*      write(6,*) ' tol2e, dentol ', tol2e, dentol
c
      if (.not. bas_numcont(basis, nsh)) call errquit('bas?',0,
     &       BASIS_ERR)
      if (nsh .gt. maxcn) call errquit('maxcn?',0, UNKNOWN_ERR)
      do i = 1, nsh
         if (.not. bas_cn2ce(basis, i, cn_to_ce(i)))
     $        call errquit('bas?',0, BASIS_ERR)
      end do
c
c     Form list of sorted interacting shell pairs
c
      if (.not. geom_ncent(geom,ncent)) call errquit('ftx:geom?',0,
     &       GEOM_ERR)
      if (.not. geom_ncent_unique(geom,ncent_unique)) 
     $     call errquit('ftx:geom?',0, GEOM_ERR)
      max_p = int(dble(max_q*ncent)/dble(ncent_unique))
      if (oikjl) max_p = max_p*2
      max_p = sqrt(dble(max_p))
c      
      if(.not.ma_alloc_get(MT_int, 2*max_pair_t,'ijinfo', 
     .     l_ijinfo, k_ijinfo)) call errquit(
     .     ' fock2etask: pushget failed ',2*max_pair_t, MA_ERR)
      call fock_pairs(basis, max_p, oij, iatlo, iathi, jatlo, jathi, 
     $     dentol, ijlist, sij, int_mb(k_ijinfo), max_pair_t, nt_ij)
      if(.not.ma_alloc_get(MT_int, 2*max_pair_t,'klinfo', 
     .     l_klinfo, k_klinfo)) call errquit(
     .     ' fock2etask: pushget failed ',2*max_pair_t, MA_ERR)
      call fock_pairs(basis, max_p, okl, katlo, kathi, latlo, lathi, 
     $     dentol, kllist, skl, int_mb(k_klinfo), max_pair_t, nt_kl)
c
c     Loop over pairs of pairs (have i>=j and k>=l already)
c
      nq = 0
      do ijt = 1,nt_ij
*         write(6,*) ' nij ', ijinfo(2,ijt)-ijinfo(1,ijt)+1
         do klt = 1, nt_kl
*            write(6,*) '      nkl ', klinfo(2,klt)-klinfo(1,klt)+1
            do ij = int_mb(k_ijinfo+(ijt-1)*2), 
     ,              int_mb(k_ijinfo+(ijt-1)*2+1) 
               ish = ijlist(1,ij)
               jsh = ijlist(2,ij)
               iat = cn_to_ce(ish)-iatlo+1
               jat = cn_to_ce(jsh)-jatlo+1
               q2 = 1.0d0
               if (oij .and. ish.eq.jsh) q2 = q2*0.5d0
               do kl = int_mb(k_klinfo+(klt-1)*2), 
     ,              int_mb(k_klinfo+(klt-1)*2+1) 
                  qq4 = q2
                  ksh = kllist(1,kl)
                  lsh = kllist(2,kl)
                  kat = cn_to_ce(ksh)-katlo+1
                  lat = cn_to_ce(lsh)-latlo+1
c     ikjl diagonal
                  if (oikjl) then
                     if (ish.lt.ksh) then
                        goto 300 ! Next kl
                     else if (ish.eq.ksh) then
                        if (jsh.lt.lsh) then
                           goto 300 ! Next kl
                        else if (jsh.eq.lsh) then
                           qq4 = qq4 * 0.5d0
                        end if
                     end if
                  end if
c     symmetry test
                  if (oskel) then
                     qq4 = qq4*q4_at(iat,jat,kat,lat)
                     if (qq4 .eq. 0.0d0) goto 300 ! Next kl
                  end if
c     sparsity test
                  sijkl = sij(ij)*skl(kl)
                  if (sijkl .lt. dentol) goto 300 ! Next kl
c     density screening
                  if (odensityscreen .and. (.not. owritefile)) then
                     sijkl=sijkl*max(
     $                    rdij(iat,jat),rdik(iat,kat),rdli(lat,iat),
     $                    rdjk(jat,kat),rdlj(lat,jat),rdlk(lat,kat))
                     if (sijkl .lt. tol2e) goto 300 ! Next kl
                  end if
c     kl diagonal
                  if (okl .and. ksh.eq.lsh) qq4 = qq4*0.5d0
c     do it!
                  if (nq .eq. max_q) then
                     call fock_txs(basis, nfock, vg_dens, vg_fock,
     $                    ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $                    dij, dik, dli, djk, dlj, dlk,
     $                    fij, fik, fli, fjk, flj, flk,
     $                    dentol, nq, ijkl(1,1), ijkl(1,2),
     $                    ijkl(1,3),ijkl(1,4), q4)
                     nq = 0
                  end if
                  nq = nq + 1
                  ijkl(nq,1) = ish
                  ijkl(nq,2) = jsh
                  ijkl(nq,3) = ksh
                  ijkl(nq,4) = lsh
                  q4(nq)     = qq4
 300           end do            ! Next KL pair
            end do
            call fock_txs(basis, nfock, vg_dens, vg_fock,
     $           ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $           dij, dik, dli, djk, dlj, dlk,
     $           fij, fik, fli, fjk, flj, flk,
     $           dentol, nq, ijkl(1,1), ijkl(1,2),
     $           ijkl(1,3),ijkl(1,4), q4)
            nq = 0
         end do
      end do
      if (.not. ma_free_heap(l_ijinfo))
     $     call errquit('fock2etask: cannot free iji',l_ijinfo, MA_ERR)
      if (.not. ma_free_heap(l_klinfo))
     $     call errquit('fock2etask: cannot free kli',l_klinfo, MA_ERR)
c     
      end
      subroutine fock_txs(basis, nfock, vg_dens, vg_fock,
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tol2e, nq, iq, jq, kq, lq, q4)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cfock.fh"
#include "util.fh"
#include "stdio.fh"
c
      double precision ncomputed
      common /ctemp/ ncomputed
c
      integer basis
      integer nfock
      integer ilo, jlo, klo, llo, ihi, jhi, khi, lhi
      integer vg_dens(nfock), vg_fock(nfock)
      double precision dij(nfock,*), dik(nfock,*), dli(nfock,*),
     $     djk(nfock,*), dlj(nfock,*), dlk(nfock,*)
      double precision fij(nfock,*), fik(nfock,*), fli(nfock,*),
     $     fjk(nfock,*), flj(nfock,*), flk(nfock,*)
      double precision tol2e
      integer nq
      integer iq(*), jq(*), kq(*), lq(*)
      double precision q4(*)
c
      integer lenscr, leneri, lenlab, maxg
      integer l_eri, k_eri, l_lab, k_lab, l_scr, k_scr
      integer k_i, k_j, k_k, k_l
c
      integer neri
      double precision block_eff
      logical more
c
      logical int2e_file_write, intb_2e4c, intb_init4c
      external int2e_file_write, intb_2e4c, intb_init4c
c
      if (nq .le. 0) return
*      write(6,*) ' fock_txs : ', nq
c
      if (maxscr .eq. 0) then
*1 quartet         call int_mem_2e4c(maxg, lenscr)
         call intb_mem_2e4c(maxg, lenscr) ! blocking algorithm
      else
         lenscr = maxscr
      end if
      leneri = max(maxeri,maxg)
      lenlab = 4*leneri
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri))
     $     call errquit('fock_txs: no room for eri',leneri, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr))
     $     call errquit('fock_txs: no room for scr',lenscr, MA_ERR)
      if (.not. ma_push_get(mt_int, lenlab, 'lab', l_lab, k_lab))
     $     call errquit('fock_txs: no room for lab',lenlab, MA_ERR)
c
      k_i = k_lab
      k_j = k_i   + leneri
      k_k = k_j   + leneri
      k_l = k_k   + leneri
c
      if (.not. intb_init4c(basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., lenscr, dbl_mb(k_scr), leneri, block_eff))
     $     call errquit('intb_init?',nq, INT_ERR)
c     
 50   more = intb_2e4c(basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., tol2e, .false.,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri), leneri, neri, lenscr, dbl_mb(k_scr))
c
      if (neri .gt. 0) then
c
         ncomputed = ncomputed + neri
c     
c     Cache integrals if desired
c     
         if (owritefile) then
            if (.not. int2e_file_write(neri,
     $           int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $           dbl_mb(k_eri))) then
                  write(luout,*) 
                  write(luout,*) ' CHANGE TO DIRECT SCF '
                  write(luout,*) 
               call errquit(
     .        ' fock_2e_task: write failed with neri=',neri, INT_ERR)
            endif
         end if
c     
         call fock_2e_label(nfock, tol2e, neri,
     $        int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $        dbl_mb(k_eri),
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $        dij, dik, dli, djk, dlj, dlk,
     $        fij, fik, fli, fjk, flj, flk )
      endif
c     
      if (more) goto 50
c
 100  if (.not. ma_pop_stack(l_lab))
     $     call errquit('fock_txs: cannot free lab',lenlab, MA_ERR)
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('fock_txs: cannot free scr',lenlab, MA_ERR)
      if (.not. ma_pop_stack(l_eri))
     $     call errquit('fock_txs: cannot free eri',lenlab, MA_ERR)
c
      end
      subroutine fock_den_sum(basis, nfock,
     $     iatlo, iathi, jatlo, jathi,
     $     dij, ilo, ihi, jlo, jhi, rdij, ldim, tol)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
c
      integer basis
      integer nfock
      double precision tol
      integer iatlo, iathi, jatlo, jathi
      integer ilo, ihi, jlo, jhi
      double precision dij(nfock,ilo:ihi,jlo:jhi)
      integer ldim
      double precision rdij(ldim,*)
c
      integer iat, jat, i, j, ifock
      integer ibflo, ibfhi, jbflo, jbfhi
      double precision sum
c
      do iat = iatlo, iathi
         if (.not. bas_ce2bfr(basis,iat,ibflo,ibfhi))
     $        call errquit('fock_den_sum :bas_ce2bfr iat', iat,
     &       BASIS_ERR)
         do jat = jatlo, jathi
            if (.not. bas_ce2bfr(basis,jat,jbflo,jbfhi))
     $           call errquit('fock_den_sum :bas_ce2bfr jat', jat,
     &       BASIS_ERR)
            sum = 0.0d0
            if (schwarz_atom(iat,jat).ge.tol) then
               if (nfock .eq. 1) then
                  do j = jlo, jhi
                     do i = ilo, ihi
                        sum = sum + abs(dij(1,i,j))
                     end do
                  end do
               else
                  do j = jlo, jhi
                     do i = ilo, ihi
                        do ifock = 1, nfock
                           sum = sum + abs(dij(ifock,i,j))
                        end do
                     end do
                  end do
               end if
            end if
            rdij(iat-iatlo+1,jat-jatlo+1) = sum
         end do
      end do
c
      end
