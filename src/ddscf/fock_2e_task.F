      subroutine fock_2e_task(
     $     geom, basis, oskel,
     $     iatlo, jatlo, katlo, latlo,
     $     iathi, jathi, kathi, lathi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, tol2e, den_tol)
      implicit none
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "cscfps.fh"
#include "sym.fh"
c
c     $Id: fock_2e_task.F,v 1.11 1996-01-04 21:28:59 d3g681 Exp $
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
      double precision tol2e, den_tol
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
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, iathi, junk, ihi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
c
      if (.not. bas_ce2bfr(basis, jatlo, jlo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, jathi, junk, jhi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
c
      if (.not. bas_ce2bfr(basis, katlo, klo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, kathi, junk, khi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
c
      if (.not. bas_ce2bfr(basis, latlo, llo, junk))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, lathi, junk, lhi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
c
      if (owritefile) then
         if (.not. int2e_set_bf_range(
     $        ilo, ihi, jlo, jhi, klo, khi, llo, lhi))
     $        call errquit('should not happen', 0)
      end if
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
      smax = schwarz_max()
      oij   = iatlo.eq.jatlo
      okl   = katlo.eq.latlo
      oikjl = (iatlo.eq.katlo) .and. (jatlo.eq.latlo)
c
      if (oscfps) call pstat_on(ps_fock_add)
#ifdef DDSCF_TEXAS      
      call fock_2e_texas(
     $     geom, basis, oskel,
     $     iatlo, jatlo, katlo, latlo,
     $     iathi, jathi, kathi, lathi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac,
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, tol2e, den_tol, smax, oij, okl, oikjl)
#else
      do iat = iatlo, iathi
         jattop = jathi
         if (oij) jattop = iat
         do jat = jatlo, jattop
            sij = schwarz_atom(iat,jat)
            otestij = (sij*smax) .gt. tol2e
            if (otestij) then
               kattop = kathi
               if (oikjl) kattop = iat
               do kat = katlo, kattop
                  lattop = lathi
                  if (okl) lattop = kat
                  if (oikjl .and. kat.eq.iat) lattop = jat
                  do lat = latlo,lattop
                     sijkl = sij*schwarz_atom(kat,lat)
                     otest = sijkl .gt. tol2e
                     if (otest .and. oskel) otest = sym_atom_quartet
     $                    (geom, iat, jat, kat, lat, q4)
                     if (otest) then
                        call fock_2e_b(basis, nfock, sijkl, den_tol, q4,
     $                       iat, jat, kat, lat, 
     $                       ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $                       dij, dik, dli, djk, dlj, dlk,
     $                       fij, fik, fli, fjk, flj, flk)
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
#endif
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
#ifdef DDSCF_TEXAS
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
     $     tmp, tol2e, den_tol, smax, oij, okl, oikjl)
      implicit none
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "cscfps.fh"
#include "sym.fh"
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
      double precision tol2e, den_tol
      double precision smax
      logical oij, okl, oikjl
      integer iat, jat, kat, lat, nsh
c
      integer max_sh_t, max_sh
      parameter (max_sh_t = 6, max_sh = 255)  ! i functions, 
      integer i_max, i_sh_r(2,-1:max_sh_t), i_sh(max_sh)
      integer j_max, j_sh_r(2,-1:max_sh_t), j_sh(max_sh)
      integer k_max, k_sh_r(2,-1:max_sh_t), k_sh(max_sh)
      integer l_max, l_sh_r(2,-1:max_sh_t), l_sh(max_sh)
c
      integer max_q, max_at
      parameter (max_q = 100001, max_at=10) ! For debug only
      integer ijkl(max_q,4)
      double precision q4(max_q), q4_at(max_at,max_at,max_at,max_at)
c
      integer lscr, leri
      parameter (lscr=1000000,leri=1000000)
      integer labels(leri,4)
      double precision scr(lscr), eri(leri), block_eff
c
      integer i, j, k, l
      integer it, jt, kt, lt
      integer ish, jsh, ksh, lsh
      integer jtop, ktop, ltop
      integer nq, neri
c
      integer cn_to_ce(1000)
c
      logical more, intb_2e4c, intb_init4c
      external intb_2e4c, intb_init4c
c
      logical fock_shell_logic, odoit
      external fock_shell_logic
      double precision qq4
c
c     Form q4 factors over atoms and mapping from shells to atoms
c     ... very crude now for debug purposes
c
      do iat = iatlo, iathi
         do jat = jatlo, jathi
            do kat = katlo, kathi
               do lat = latlo, lathi
                  odoit = sym_atom_quartet(geom, iat, jat, kat, lat,
     $                 q4_at(iat-iatlo+1,jat-jatlo+1,
     $                 kat-katlo+1,lat-latlo+1))
                  write(6,*) iat,jat,kat,lat,
     $                 q4_at(iat-iatlo+1,jat-jatlo+1,
     $                 kat-katlo+1,lat-latlo+1)
               enddo
            enddo
         enddo
      enddo
c
      if (.not. bas_numcont(basis, nsh)) call errquit('bas?',0)
      do i = 1, nsh
         if (.not. bas_cn2ce(basis, i, cn_to_ce(i)))
     $        call errquit('bas?',0)
      enddo
      write(6,*) ' Center map ', (cn_to_ce(i),i=1,nsh)
c
c
c     Form list of shells in each atom block, sorted by shell type. 
c
      call fock_2e_sh_list(basis, iatlo, iathi, i_max, i_sh_r, i_sh)
      call fock_2e_sh_list(basis, jatlo, jathi, j_max, j_sh_r, j_sh)
      call fock_2e_sh_list(basis, katlo, kathi, k_max, k_sh_r, k_sh)
      call fock_2e_sh_list(basis, latlo, lathi, l_max, l_sh_r, l_sh)
c
c     loop thru blocks of shell types and then over shells themselves
c
      do it = -1, max_sh_t
         if (i_sh_r(1,it).gt.i_sh_r(2,it)) goto 10
         do jt = -1, max_sh_t
            if (j_sh_r(1,jt).gt.j_sh_r(2,jt)) goto 20
            do kt = -1, max_sh_t
               if (k_sh_r(1,kt).gt.k_sh_r(2,kt)) goto 30
               do lt = -1, max_sh_t
                  if (l_sh_r(1,lt).gt.l_sh_r(2,lt)) goto 40
c
                  nq = 0
                  do i = i_sh_r(1,it),i_sh_r(2,it)
                     ish = i_sh(i)
                     iat = cn_to_ce(ish)-iatlo+1
                     do j = j_sh_r(1,jt),j_sh_r(2,jt)
                        jsh = j_sh(j)
                        jat = cn_to_ce(jsh)-jatlo+1
                        if ((.not. oij) .or. (ish.ge.jsh)) then
                           do k = k_sh_r(1,kt),k_sh_r(2,kt)
                              ksh = k_sh(k)
                              kat = cn_to_ce(ksh)-katlo+1
                              do l = l_sh_r(1,lt),l_sh_r(2,lt)
                                 lsh = l_sh(l)
                                 lat = cn_to_ce(lsh)-latlo+1
                                 odoit = fock_shell_logic(
     $                                oij, okl, oikjl,
     $                                ish, jsh, ksh, lsh, qq4)
c     Temp. discard factors of 2/4/8
                                 qq4 = qq4*q4_at(iat,jat,kat,lat)
                                 odoit = odoit .and. (qq4.gt.0.0d0)
                                 if (odoit) then
                                    nq = nq + 1
                                    ijkl(nq,1) = ish
                                    ijkl(nq,2) = jsh
                                    ijkl(nq,3) = ksh
                                    ijkl(nq,4) = lsh
                                    q4(nq)     = qq4
                                 endif
                              enddo
                           enddo
                        endif
                     enddo
                  enddo
c     This will be a REALLY wierd canonical shell order
c$$$                  nq = 0
c$$$                  do i = i_sh_r(1,it),i_sh_r(2,it)
c$$$                     ish = i_sh(i)
c$$$                     jtop = j_sh_r(2,jt)
c$$$                     if (oij) jtop = i
c$$$                     do j = j_sh_r(1,jt),jtop
c$$$                        jsh = j_sh(j)
c$$$                        ktop = k_sh_r(2,kt)
c$$$                        if (oikjl) ktop = i
c$$$                        do k = k_sh_r(1,kt),ktop
c$$$                           ksh = k_sh(k)
c$$$                           ltop = l_sh_r(2,lt)
c$$$                           if (okl) ltop = k
c$$$                           if (oikjl .and. k.eq.i) ltop = j
c$$$                           do l = l_sh_r(1,lt),ltop
c$$$                              lsh = l_sh(l)
c$$$                              nq = nq + 1
c$$$                              ijkl(nq,1) = ish
c$$$                              ijkl(nq,2) = jsh
c$$$                              ijkl(nq,3) = ksh
c$$$                              ijkl(nq,4) = lsh
c$$$                              q4(nq)     = 1.0d0
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c
                  if (nq .gt. 0) then
                     write(6,11)  it, jt, kt, lt, nq 
 11                  format( ' it, jt, kt, lt, nq ',4i4,2x,i4)
c
                     if (.not. intb_init4c(
     $                    basis, ijkl(1,1), ijkl(1,2),
     $                    basis, ijkl(1,3), ijkl(1,4),
     $                    nq, lscr, scr, leri, block_eff))
     $                    call errquit('intb_init?',nq)
c
 50                  more = intb_2e4c(
     $                    basis, ijkl(1,1), ijkl(1,2),
     $                    basis, ijkl(1,3), ijkl(1,4),
     $                    nq, q4, .true., tol2e, .false.,
     $                    labels(1,1), labels(1,2),
     $                    labels(1,3), labels(1,4),
     $                    eri, leri, neri, lscr, scr)
c
*                     call print_integ_list(neri, 
*     $                    labels(1,1), labels(1,2),
*     $                    labels(1,3), labels(1,4),
*     $                    eri)
c
c     Only need scale with nwchem blocking
c
*                     call fock_eri_label_scale(neri, 1.0d0,
*     $                    labels(1,1), labels(1,2),
*     $                    labels(1,3), labels(1,4), eri)
c
                     call fock_2e_label(nfock, tol2e, neri,
     $                    labels(1,1), labels(1,2),
     $                    labels(1,3), labels(1,4),
     $                    eri, 
     $                    ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $                    dij, dik, dli, djk, dlj, dlk,
     $                    fij, fik, fli, fjk, flj, flk )
c
                     if (more) goto 50
c
                  endif
 40            enddo
 30         enddo
 20      enddo
 10   enddo
c                              
      end
      subroutine fock_2e_sh_list(basis, iatlo, iathi,
     $     i_max, i_sh_r, i_sh)
      implicit none
#include "bas.fh"
      integer basis
      integer iatlo, iathi
      integer i_max
      integer i_sh_r(2,-1:*)
      integer i_sh(*)
c
      integer max_sh_t, max_sh
      parameter (max_sh_t = 6, max_sh = 255)  ! i functions
      integer ilo, ihi, junk, i
      integer type(max_sh), nprim(max_sh), ngen(max_sh), spch
      integer nt(-1:max_sh_t)
c
      integer i, ii, jj, tmp
      logical otype, ongen
c
      if (.not. bas_ce2cnr(basis, iatlo, ilo, junk)) call errquit
     $     ('fock_2e_sh_list: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, iathi, junk, ihi)) call errquit

c
c     Sort by type.  If the types are equal subsort by ngen.
c     If type and ngen are equal subsort by nprim.
c
      call ifill(max_sh_t+2,0,nt,1)
      do i = ilo, ihi
         ii = i - ilo + 1
         if (.not. bas_continfo(basis, i, type(ii), nprim(ii),
     $        ngen(ii), spch))
     $        call errquit('fock_2e_sh_list: basis screwed ',0)
         i_sh(ii) = i
         nt(type(ii)) = nt(type(ii)) + 1
         do jj = 1, ii-1
c
            otype = type(ii).eq.type(jj)
            ongen = (ngen(ii).eq.ngen(jj)) .and. otype
            if ( (type(ii).lt.type(jj)) .or.
     $           (otype .and. (ngen(ii).lt.ngen(jj))) .or.
     $           (ongen .and. (nprim(ii).lt.nprim(jj)))) then
               tmp = type(ii)
               type(ii) = type(jj)
               type(jj) = tmp
               tmp = ngen(ii)
               ngen(ii) = ngen(jj)
               ngen(jj) = tmp
               tmp = nprim(ii)
               nprim(ii) = nprim(jj)
               nprim(jj) = tmp
               tmp = i_sh(ii)
               i_sh(ii) = i_sh(jj)
               i_sh(jj) = tmp
            endif
         enddo
      enddo
c
c     Now form the ranges for each type ... for now just ignore
c     the subsorting by ngen/nprim but we may come back to this later.
c     If a certain angular momentum is not present then the index range 
c     will form a null do loop.
c
      i_sh_r(1,-1) = 1
      do i = -1, max_sh_t
         if (i.ne.-1) i_sh_r(1,i) = i_sh_r(2,i-1) + 1
         i_sh_r(2,i) = i_sh_r(1,i) + nt(i) - 1
      enddo
c
*      write(6,*) ' Final sorted info for shell range ', ilo, ihi
*      do i = -1, max_sh_t
*         if (i_sh_r(2,i) .ge. i_sh_r(1,i)) then
*            write(6,*) i, ' sorted shell range ',
*     $           i_sh_r(1,i),i_sh_r(2,i)
*            write(6,*) (i_sh(ii),ii=i_sh_r(1,i),i_sh_r(2,i))
*         endif
*      enddo
c
      end
      logical function fock_shell_logic(oij, okl, oikjl,
     $     ish, jsh, ksh, lsh, qq4)
      implicit none
      logical oij, okl, oikjl
      integer ish, jsh, ksh, lsh
      double precision qq4
c
      logical odoit
c
      odoit = .true.
      qq4 = 1.0d0
      if (oij) then
         if (ish.lt.jsh) then
            odoit = .false.
         else if (ish.eq.jsh) then
            qq4 = qq4 * 0.5d0
         endif
      endif
      if (okl) then
         if (ksh.lt.lsh) then
            odoit = .false.
         else if (ksh.eq.lsh) then
            qq4 = qq4 * 0.5d0
         endif
      endif
      if (oikjl) then
         if (ish.lt.ksh) then
            odoit = .false.
         else if (ish.eq.ksh) then
            if (jsh.lt.lsh) then
               odoit = .false.
            else if (jsh.eq.lsh) then
               qq4 = qq4 * 0.5d0
            endif
         endif
      endif
c
      fock_shell_logic = odoit
c
      end
#endif
