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
c     $Id: fock_2e_task.F,v 1.16 1996-01-22 18:33:31 d3g681 Exp $
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
     $     tmp, tol2e, den_tol, oij, okl, oikjl)
#else
      smax = schwarz_max()
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
     $     tmp, tol2e, den_tol, oij, okl, oikjl)
      implicit none
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
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
      double precision smax, sijkl
      logical oij, okl, oikjl
      integer iat, jat, kat, lat, nsh
c
c     These for sorting shells within each atom pair block
c
      integer max_pair_t, max_sh
      parameter (max_pair_t = 1024, max_sh = 255) ! allow for each shell diff
      double precision skl(max_sh*max_sh) ! Should be just max shell block
      double precision sij(max_sh*max_sh) ! Should be just max shell block
      integer ijlist(2,max_sh*max_sh), nij, ij
      integer kllist(2,max_sh*max_sh), nkl, kl
      integer ijinfo(2,max_pair_t), nt_ij
      integer klinfo(2,max_pair_t), nt_kl
c
      integer max_q, max_at, max_p
      parameter (max_q = 24001, max_at=20) ! For debug only
      integer ijkl(max_q,4)
      double precision q4(max_q), q4_at(max_at,max_at,max_at,max_at)
      double precision denmax(max_at,max_at,max_at,max_at)
c
      integer ncent, ncent_unique
      integer i, j, k, l
      integer ia,ja,ka,la
      integer it, jt, kt, lt, ijt, klt
      integer ish, jsh, ksh, lsh
      integer ibflo, jbflo, kbflo, lbflo, ibfhi, jbfhi, kbfhi, lbfhi
      integer nq
c
      integer cn_to_ce(1000)
c
      logical fock_shell_logic, odoit, status
      external fock_shell_logic
      double precision q2, qq4, test
c
c     Form q4 factors over atoms and mapping from shells to atoms
c     ... very crude now for debug purposes
c
      if ((iathi-iatlo+1).gt.max_at) call errquit('too many ats',0)
      if ((jathi-jatlo+1).gt.max_at) call errquit('too many ats',0)
      if ((kathi-katlo+1).gt.max_at) call errquit('too many ats',0)
      if ((lathi-latlo+1).gt.max_at) call errquit('too many ats',0)
      do iat = iatlo, iathi
         do jat = jatlo, jathi
            do kat = katlo, kathi
               do lat = latlo, lathi
                  ia = iat-iatlo+1
                  ja = jat-jatlo+1
                  ka = kat-katlo+1
                  la = lat-latlo+1
                  odoit = sym_atom_quartet(geom, iat, jat, kat, lat,
     $                 q4_at(ia,ja,ka,la))
                  if (.not. owritefile) then
                     status = bas_ce2bfr(basis,iat,ibflo,ibfhi)
                     status = bas_ce2bfr(basis,jat,jbflo,jbfhi)
                     status = bas_ce2bfr(basis,kat,kbflo,kbfhi)
                     status = bas_ce2bfr(basis,lat,lbflo,lbfhi)
                     call fock_density_screen(nfock,
     $                    ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $                    ibflo, jbflo, kbflo, lbflo,
     $                    ibfhi, jbfhi, kbfhi, lbfhi,
     $                    dij, dik, dli, djk, dlj, dlk,
     $                    denmax(ia,ja,ka,la))
                  endif
               enddo
            enddo
         enddo
      enddo
      smax = schwarz_max()
c
      if (.not. bas_numcont(basis, nsh)) call errquit('bas?',0)
      do i = 1, nsh
         if (.not. bas_cn2ce(basis, i, cn_to_ce(i)))
     $        call errquit('bas?',0)
      enddo
c
c     Form list of sorted interacting shell pairs
c
      if (.not. geom_ncent(geom,ncent)) call errquit('ftx:geom?',0)
      if (.not. geom_ncent_unique(geom,ncent_unique)) 
     $     call errquit('ftx:geom?',0)
      max_p = dble(max_q*ncent)/dble(ncent_unique)
*      if (oikjl) max_p = max_p*2
      max_p = sqrt(dble(max_p))
*      max_p = 109
c      
      call fock_pairs(basis, max_p, oij, iatlo, iathi, jatlo, jathi, 
     $     tol2e, ijlist, sij, ijinfo, nt_ij)
      call fock_pairs(basis, max_p, okl, katlo, kathi, latlo, lathi, 
     $     tol2e, kllist, skl, klinfo, nt_kl)
c
c     Loop over pairs of pairs (have i>=j and k>=l already)
c
      nq = 0
      do ijt = 1,nt_ij
         write(6,*) ' nij ', ijinfo(2,ijt)-ijinfo(1,ijt)+1
         do klt = 1, nt_kl
            write(6,*) '      nkl ', klinfo(2,klt)-klinfo(1,klt)+1
            do ij = ijinfo(1,ijt),ijinfo(2,ijt)
               ish = ijlist(1,ij)
               jsh = ijlist(2,ij)
               iat = cn_to_ce(ish)-iatlo+1
               jat = cn_to_ce(jsh)-jatlo+1
               if (cn_to_ce(ish).lt.iatlo .or. 
     $              cn_to_ce(ish).gt.iathi) call errquit('bad iat',0)
               if (cn_to_ce(jsh).lt.jatlo .or. 
     $              cn_to_ce(jsh).gt.jathi) call errquit('bad jat',0)
               q2 = 1.0d0
               if (oij .and. ish.eq.jsh) q2 = q2*0.5d0
               do kl = klinfo(1,klt),klinfo(2,klt)
                  qq4 = q2
                  ksh = kllist(1,kl)
                  lsh = kllist(2,kl)
                  kat = cn_to_ce(ksh)-katlo+1
                  lat = cn_to_ce(lsh)-latlo+1
                  if (cn_to_ce(ksh).lt.katlo .or. 
     $                 cn_to_ce(ksh).gt.kathi)
     $                 call errquit('bad kat',0)
                  if (cn_to_ce(lsh).lt.latlo .or. 
     $                 cn_to_ce(lsh).gt.lathi)
     $                 call errquit('bad lat',0)
c     ikjl diagonal
                  if (oikjl) then
                     if (ish.lt.ksh) then
                        goto 300 ! Next kl
                     else if (ish.eq.ksh) then
                        if (jsh.lt.lsh) then
                           goto 300 ! Next kl
                        else if (jsh.eq.lsh) then
                           qq4 = qq4 * 0.5d0
                        endif
                     endif
                  endif
c     symmetry test
                  if (oskel) then
                     qq4 = qq4*q4_at(iat,jat,kat,lat)
                     if (qq4 .eq. 0.0d0) goto 300 ! Next kl
                  endif
c     sparsity test
                  sijkl = sij(ij)*skl(kl)
                  if (sijkl .lt. tol2e) goto 300 ! Next kl
c     density screening
                  if (.not. owritefile) then
                     sijkl=sijkl*denmax(iat,jat,kat,lat)
                     if (sijkl .lt. den_tol) goto 300 ! Next kl
                  endif
c     kl diagonal
                  if (okl .and. ksh.eq.lsh) qq4 = qq4*0.5d0
c     do it!
                  if (nq .eq. max_q) then
                     call fock_txs(basis, nfock, vg_dens, vg_fock,
     $                    ilo, jlo, klo, llo, ihi, jhi,khi,lhi,
     $                    dij, dik, dli, djk, dlj, dlk,
     $                    fij, fik, fli, fjk, flj, flk,
     $                    tol2e, nq, ijkl(1,1),ijkl(1,2),
     $                    ijkl(1,3),ijkl(1,4), q4)
                     nq = 0
                  endif
                  nq = nq + 1
                  ijkl(nq,1) = ish
                  ijkl(nq,2) = jsh
                  ijkl(nq,3) = ksh
                  ijkl(nq,4) = lsh
                  q4(nq)     = qq4
 300           enddo            ! Next KL pair
            enddo
            call fock_txs(basis, nfock, vg_dens, vg_fock,
     $           ilo, jlo, klo, llo, ihi, jhi,khi,lhi,
     $           dij, dik, dli, djk, dlj, dlk,
     $           fij, fik, fli, fjk, flj, flk,
     $           tol2e, nq, ijkl(1,1),ijkl(1,2),
     $           ijkl(1,3),ijkl(1,4), q4)
            nq = 0
         enddo
      enddo
c     
      end
      subroutine fock_pairs(basis, max_p, oij, 
     $     iatlo, iathi, jatlo, jathi, 
     $     tol2e, ijlist, sij, ijinfo, nt_ij)
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis
      integer max_p
      logical oij
      integer iatlo, iathi, jatlo, jathi
      integer ijlist(2,*), ijinfo(2,*), nt_ij
      double precision sij(*)
c
      integer junk
      integer nij, ij, ijt, i, j, jtop, prev, spch
      integer jlo, jhi, jnprim, jngen, jt, jtype
      integer ilo, ihi, inprim, ingen, it, itype
      double precision smax, test, tol2e
      integer ijtype(255*255), n_in_b, old_nt
c
#include "itri.fh"
c
      if (.not. bas_ce2cnr(basis, iatlo, ilo, junk)) call errquit
     $     ('fock_pairs: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, iathi, junk, ihi)) call errquit
     $     ('fock_pairs: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, jatlo, jlo, junk)) call errquit
     $     ('fock_pairs: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, jathi, junk, jhi)) call errquit
     $     ('fock_pairs: basis screwed ',0)
c
      smax = schwarz_max()
c
      nt_ij = 0
      nij = 0
      do i = ilo, ihi
         if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
     $        call errquit('fock_pairs: basis screwed ',0)
         it = (1+itype)*1000+inprim*10+ingen
         jtop = jhi
         if (oij) jtop = i
         do j = jlo, jtop
            test = schwarz_shell(i,j)
            if (test*smax .gt. tol2e) then
               if (.not. bas_continfo(basis,j,jtype,jnprim,jngen,spch))
     $              call errquit('fock_pairs: basis screwed ',0)
               jt = (1+jtype)*1000+jnprim*10+jngen
               nij           = nij + 1
               ijlist(1,nij) = i
               ijlist(2,nij) = j
               ijtype(nij)   = itri(it,jt)
            endif
         enddo
      enddo
      if (nij .eq. 0) return
c
      call fock_pair_sort(nij,ijlist,ijtype)
c
c     Split blocks by type only
c
c$$$      nt_ij = 0
c$$$      prev = -99999999
c$$$      do ij = 1, nij
c$$$         if (ijtype(ij) .ne. prev) then
c$$$            if (nt_ij .ne. 0) ijinfo(2,nt_ij) = ij-1
c$$$            nt_ij = nt_ij + 1
c$$$            ijinfo(1,nt_ij) = ij
c$$$            prev = ijtype(ij)
c$$$         endif
c$$$      enddo
c$$$      ijinfo(2,nt_ij) = nij
c
c     Split blocks by size only
c
c$$$      nt_ij = 0
c$$$      do ij = 1, nij, max_p
c$$$         nt_ij = nt_ij + 1
c$$$         ijinfo(1,nt_ij) = ij
c$$$         ijinfo(2,nt_ij) = min(nij,ij + max_p - 1)
c$$$      enddo
c
c     Split pair blocks on type boundaries with regard to size
c     (first do type only then combine)
c
      nt_ij = 0
      prev = -99999999
      do ij = 1, nij
         if (ijtype(ij) .ne. prev) then
            if (nt_ij .ne. 0) ijinfo(2,nt_ij) = ij-1
            nt_ij = nt_ij + 1
            ijinfo(1,nt_ij) = ij
            prev = ijtype(ij)
         endif
      enddo
      ijinfo(2,nt_ij) = nij
c
      old_nt = nt_ij
      nt_ij = 1
      n_in_b = (ijinfo(2,1)-ijinfo(1,1)+1)
      do ijt = 2, old_nt
         n_in_b = n_in_b + (ijinfo(2,ijt)-ijinfo(1,ijt)+1)
         write(6,*) ijt, n_in_b, max_p
         if (n_in_b .le. max_p) then
            ijinfo(2,nt_ij) = ijinfo(2,ijt)
            write(6,*) ' setting end   ', nt_ij, ijt
         else
            nt_ij = nt_ij + 1
            write(6,*) ' setting start ', nt_ij, ijt
            ijinfo(1,nt_ij) = ijinfo(1,ijt)
            ijinfo(2,nt_ij) = ijinfo(2,ijt)
            n_in_b = (ijinfo(2,ijt)-ijinfo(1,ijt)+1)
         endif
      enddo
      ijinfo(2,nt_ij) = nij
c     
c     Gather schwarz info again
c
      do ij = 1, nij
         i = ijlist(1,ij)
         j = ijlist(2,ij)
         sij(ij) = schwarz_shell(i,j)
      enddo
c
      do ijt = 1, nt_ij
         write(6,*) ' ij type ', ijt, ijinfo(1,ijt),ijinfo(2,ijt)
c$$$         do ij = ijinfo(1,ijt),ijinfo(2,ijt)
c$$$            i = ijlist(1,ij)
c$$$            j = ijlist(2,ij)
c$$$            if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
c$$$     $           call errquit('fock_pairs: basis screwed ',0)
c$$$            if (.not. bas_continfo(basis,j,jtype,jnprim,jngen,spch))
c$$$     $           call errquit('fock_pairs: basis screwed ',0)
c$$$            it = (1+itype)*1000+inprim*10+ingen
c$$$            jt = (1+jtype)*1000+jnprim*10+jngen
c$$$            if (it .ge. jt) then
c$$$               write(6,11) ij,i,j,itype,jtype,inprim,jnprim,ingen,jngen
c$$$            else
c$$$               write(6,11) ij,j,i,jtype,itype,jnprim,inprim,jngen,ingen
c$$$            endif
c$$$ 11         format(1x,i6,2x,2i3,2x,2i3,2x,2i3,2x,2i3)
c$$$         enddo
      enddo
c
      end
c
      subroutine fock_pair_sort(n,ijlist,ijtype)
      implicit integer (a-z)
      dimension ijlist(2,n), ijtype(n)
c     
      if (n .eq. 1) return
c     
      l=n/2+1
      ir=n
 10   continue
      if(l.gt.1)then
         l=l-1
         rra=ijtype(l)
         rrb1=ijlist(1,l)
         rrb2=ijlist(2,l)
      else
         rra=ijtype(ir)
         rrb1=ijlist(1,ir)
         rrb2=ijlist(2,ir)
         ijtype(ir)=ijtype(1)
         ijlist(1,ir)=ijlist(1,1)
         ijlist(2,ir)=ijlist(2,1)
         ir=ir-1
         if(ir.eq.1)then
            ijtype(1)=rra
            ijlist(1,1)=rrb1
            ijlist(2,1)=rrb2
            return
         endif
      endif
      i=l
      j=l+l
 20   if(j.le.ir)then
         if(j.lt.ir)then
            if(ijtype(j).lt.ijtype(j+1))j=j+1
         endif
         if(rra.lt.ijtype(j))then
            ijtype(i)=ijtype(j)
            ijlist(1,i)=ijlist(1,j)
            ijlist(2,i)=ijlist(2,j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         go to 20
      endif
      ijtype(i)=rra
      ijlist(1,i)=rrb1
      ijlist(2,i)=rrb2
      go to 10
      end
c
      subroutine fock_2e_texas_old(
     $     geom, basis, oskel,
     $     iatlo, jatlo, katlo, latlo,
     $     iathi, jathi, kathi, lathi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, tol2e, den_tol, oij, okl, oikjl)
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
      double precision smax, sijkl
      logical oij, okl, oikjl
      integer iat, jat, kat, lat, nsh
c
c     These for sorting shells within each atom block
c
      integer max_sh_t, max_sh
      parameter (max_sh_t = 255, max_sh = 255) ! allow for each shell diff
      integer i_nt, i_sh_r(2,max_sh_t), i_sh(max_sh)
      integer j_nt, j_sh_r(2,max_sh_t), j_sh(max_sh)
      integer k_nt, k_sh_r(2,max_sh_t), k_sh(max_sh)
      integer l_nt, l_sh_r(2,max_sh_t), l_sh(max_sh)
      double precision skl(max_sh*max_sh) ! Should be just max shell block
      double precision sij(max_sh*max_sh) ! Should be just max shell block
c
      integer max_q, max_at
      parameter (max_q = 12001, max_at=20) ! For debug only
      integer ijkl(max_q,4)
      double precision q4(max_q), q4_at(max_at,max_at,max_at,max_at)
      double precision denmax(max_at,max_at,max_at,max_at)
c
      integer i, j, k, l
      integer ia,ja,ka,la
      integer it, jt, kt, lt
      integer ish, jsh, ksh, lsh
      integer ibflo, jbflo, kbflo, lbflo, ibfhi, jbfhi, kbfhi, lbfhi
      integer nq
c
      integer ijlist(2,255*255), nij, ij
      integer kllist(2,255*255), nkl, kl
c
      integer cn_to_ce(1000)
c
      logical fock_shell_logic, odoit, status
      external fock_shell_logic
      double precision q2, qq4, test
c
c     Form q4 factors over atoms and mapping from shells to atoms
c     ... very crude now for debug purposes
c
      if ((iathi-iatlo+1).gt.max_at) call errquit('too many ats',0)
      if ((jathi-jatlo+1).gt.max_at) call errquit('too many ats',0)
      if ((kathi-katlo+1).gt.max_at) call errquit('too many ats',0)
      if ((lathi-latlo+1).gt.max_at) call errquit('too many ats',0)
      do iat = iatlo, iathi
         do jat = jatlo, jathi
            do kat = katlo, kathi
               do lat = latlo, lathi
                  ia = iat-iatlo+1
                  ja = jat-jatlo+1
                  ka = kat-katlo+1
                  la = lat-latlo+1
                  odoit = sym_atom_quartet(geom, iat, jat, kat, lat,
     $                 q4_at(ia,ja,ka,la))
                  if (.not. owritefile) then
                     status = bas_ce2bfr(basis,iat,ibflo,ibfhi)
                     status = bas_ce2bfr(basis,jat,jbflo,jbfhi)
                     status = bas_ce2bfr(basis,kat,kbflo,kbfhi)
                     status = bas_ce2bfr(basis,lat,lbflo,lbfhi)
                     call fock_density_screen(nfock,
     $                    ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $                    ibflo, jbflo, kbflo, lbflo,
     $                    ibfhi, jbfhi, kbfhi, lbfhi,
     $                    dij, dik, dli, djk, dlj, dlk,
     $                    denmax(ia,ja,ka,la))
                  endif
               enddo
            enddo
         enddo
      enddo
      smax = schwarz_max()
c
      if (.not. bas_numcont(basis, nsh)) call errquit('bas?',0)
      do i = 1, nsh
         if (.not. bas_cn2ce(basis, i, cn_to_ce(i)))
     $        call errquit('bas?',0)
      enddo
c
c     Form list of shells in each atom block, sorted by shell type. 
c
      call fock_2e_sh_list(basis,iatlo,iathi,i_nt,i_sh_r,i_sh)
      call fock_2e_sh_list(basis,jatlo,jathi,j_nt,j_sh_r,j_sh)
      call fock_2e_sh_list(basis,katlo,kathi,k_nt,k_sh_r,k_sh)
      call fock_2e_sh_list(basis,latlo,lathi,l_nt,l_sh_r,l_sh)
c
c     Loop thru shell types and then over shells of each type
c
      nq = 0
      do it = 1, i_nt
         do jt = 1, j_nt
c     
c     Gather ij screening information. Unfortunately, we cannot test
c     on symmetry of the ij pairs since some quartets will not be 
c     canonical (i.e., kl can be > ij)
c
            nij = 0
            do i = i_sh_r(1,it),i_sh_r(2,it)
               ish = i_sh(i)
               do j = j_sh_r(1,jt),j_sh_r(2,jt)
                  jsh = j_sh(j)
                  test = schwarz_shell(ish,jsh)
                  if ((.not. (oij .and. ish.lt.jsh)) .and.
     $                 (test*smax .gt. tol2e)) then
                     nij           = nij + 1
                     ijlist(1,nij) = ish
                     ijlist(2,nij) = jsh
                     sij(nij)      = test
                  endif
               enddo
            enddo
            if (nij .eq. 0) goto 1010
            write(6,*) ' nij ', nij
c
            do kt = 1, k_nt
               do lt = 1, l_nt
c     
c     gather kl screening information 
c
                  nkl = 0
                  do k = k_sh_r(1,kt),k_sh_r(2,kt)
                     ksh = k_sh(k)
                     do l = l_sh_r(1,lt),l_sh_r(2,lt)
                        lsh = l_sh(l)
                        test = schwarz_shell(ksh,lsh)
                        if ((.not. (okl .and. ksh.lt.lsh)) .and.
     $                       (test*smax .gt. tol2e)) then
                           nkl           = nkl + 1
                           kllist(1,nkl) = ksh
                           kllist(2,nkl) = lsh
                           skl(nkl)      = test
                        endif
                     enddo
                  enddo
                  if (nkl .eq. 0) goto 1000
c     
c     Loop over pairs of pairs (have i>=j and k>=l already)
c     
                  do ij = 1, nij
                     ish = ijlist(1,ij)
                     jsh = ijlist(2,ij)
                     iat = cn_to_ce(ish)-iatlo+1
                     jat = cn_to_ce(jsh)-jatlo+1
                     q2 = 1.0d0
                     if (oij .and. ish.eq.jsh) q2 = q2*0.5d0
                     do kl = 1, nkl
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
                              endif
                           endif
                        endif
c     symmetry test
                        if (oskel) then
                           qq4 = qq4*q4_at(iat,jat,kat,lat)
                           if (qq4 .eq. 0.0d0) goto 300 ! Next kl
                        endif
c     sparsity test
                        sijkl = sij(ij)*skl(kl)
                        if (sijkl .lt. tol2e) goto 300 ! Next kl
c     density screening
                        if (.not. owritefile) then
                           sijkl=sijkl*denmax(iat,jat,kat,lat)
                           if (sijkl .lt. den_tol) goto 300 ! Next kl
                        endif
c     kl diagonal
                        if (okl .and. ksh.eq.lsh) qq4 = qq4*0.5d0
c     do it!
                        if (nq .eq. max_q) then
                           call fock_txs(basis, nfock, vg_dens, vg_fock,
     $                          ilo, jlo, klo, llo, ihi, jhi,khi,lhi,
     $                          dij, dik, dli, djk, dlj, dlk,
     $                          fij, fik, fli, fjk, flj, flk,
     $                          tol2e, nq, ijkl(1,1),ijkl(1,2),
     $                          ijkl(1,3),ijkl(1,4), q4)
                           nq = 0
                        endif
                        nq = nq + 1
                        ijkl(nq,1) = ish
                        ijkl(nq,2) = jsh
                        ijkl(nq,3) = ksh
                        ijkl(nq,4) = lsh
                        q4(nq)     = qq4
 300                 enddo      ! Next KL pair
                  enddo
c
 1000             continue      ! Jump here if no KL pairs
c
               enddo
            enddo
 1010             continue      ! Jump here if no IJ pairs
         enddo
      enddo
c                              
      call fock_txs(basis, nfock, vg_dens, vg_fock,
     $     ilo, jlo, klo, llo, ihi, jhi,khi,lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tol2e, nq, ijkl(1,1),ijkl(1,2),
     $     ijkl(1,3),ijkl(1,4), q4)
c     
      end
      subroutine fock_txs(basis, nfock, vg_dens, vg_fock,
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tol2e, nq, iq, jq, kq, lq, q4)
      implicit none
#include "mafdecls.fh"
#include "cfock.fh"
c
      double precision ncomputed
      common /crap/ ncomputed
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
      integer lenscr, leneri, lenlab
      parameter (lenscr=3200000, leneri=200003) ! Eventually dynamic
      parameter (lenlab=4*leneri)
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
      write(6,*) ' fock_txs : ', nq
c
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri))
     $     call errquit('fock_txs: no room for eri',leneri)
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr))
     $     call errquit('fock_txs: no room for scr',lenscr)
      if (.not. ma_push_get(mt_int, lenlab, 'lab', l_lab, k_lab))
     $     call errquit('fock_txs: no room for lab',lenlab)
c
      k_i = k_lab
      k_j = k_i   + leneri
      k_k = k_j   + leneri
      k_l = k_k   + leneri
c
      if (.not. intb_init4c(basis, iq, jq, basis, kq, lq,
     $     nq, lenscr, dbl_mb(k_scr), leneri, block_eff))
     $     call errquit('intb_init?',nq)
c     
 50   more = intb_2e4c(basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., tol2e, .false.,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri), leneri, neri, lenscr, dbl_mb(k_scr))
c
      if (neri .eq. 0) goto 100 ! Free memory before going home
c
      ncomputed = ncomputed + neri
c
c     Cache integrals if desired
c
      if (owritefile) then
         if (.not. int2e_file_write(neri,
     $        int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $        dbl_mb(k_eri))) call errquit('Arrgggghh',0)
      end if
c     
      call fock_2e_label(nfock, tol2e, neri,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri),
     $     ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk )
c     
      if (more) goto 50
c
 100  if (.not. ma_pop_stack(l_lab))
     $     call errquit('fock_txs: cannot free lab',lenlab)
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('fock_txs: cannot free scr',lenlab)
      if (.not. ma_pop_stack(l_eri))
     $     call errquit('fock_txs: cannot free eri',lenlab)
c
      end
      subroutine fock_2e_sh_list(basis, iatlo, iathi,
     $     i_nt, i_sh_r, i_sh)
      implicit none
#include "bas.fh"
      integer basis
      integer iatlo, iathi
      integer i_nt
      integer i_sh_r(2,*)
      integer i_sh(*)
c
      integer max_sh
      parameter (max_sh = 255)
      integer type(max_sh), nprim(max_sh), ngen(max_sh), spch
      integer ilo, ihi, junk
c
      integer i, ii, jj, tmp, ind
      logical otype, ongen
c
      if (.not. bas_ce2cnr(basis, iatlo, ilo, junk)) call errquit
     $     ('fock_2e_sh_list: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, iathi, junk, ihi)) call errquit
     $     ('fock_2e_sh_list: basis screwed ',0)
c
c     Sort by type.  If the types are equal subsort by ngen.
c     If type and ngen are equal subsort by nprim.
c
      do i = ilo, ihi
         ii = i - ilo + 1
         if (.not. bas_continfo(basis, i, type(ii), nprim(ii),
     $        ngen(ii), spch))
     $        call errquit('fock_2e_sh_list: basis screwed ',0)
         i_sh(ii) = i
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
c     Now form the ranges for each distinct shell block
c
      i_nt = 0
      ii = 5551212
      do i = 1, ihi-ilo+1
         ind = nprim(i)+100*(ngen(i)+10*type(i))
         if (ind .ne. ii) then
            if (i_nt.ne.0) i_sh_r(2,i_nt) = i-1
            i_nt = i_nt + 1
            i_sh_r(1,i_nt) = i
            ii = ind
         endif
      enddo
      if (i_nt.ne.0) i_sh_r(2,i_nt) = ihi-ilo+1
c
*      write(6,*) ' Final sorted info for shell range ', ilo, ihi
*      do i = 1, i_nt
*         write(6,*) ' sorted shell range ', i,
*     $        i_sh_r(1,i),i_sh_r(2,i)
*         write(6,31) (i_sh(ii),type(ii),nprim(ii),ngen(ii),
*     $        ii=i_sh_r(1,i),i_sh_r(2,i))
* 31      format(1x,i5,5x,3i5)
*      enddo
c
      end
c$$$      logical function fock_shell_logic(oij, okl, oikjl,
c$$$     $     ish, jsh, ksh, lsh, qq4)
c$$$      implicit none
c$$$      logical oij, okl, oikjl
c$$$      integer ish, jsh, ksh, lsh
c$$$      double precision qq4
c$$$c
c$$$c     Not used, but it makes explicit the tests and scaling
c$$$c
c$$$      qq4 = 1.0d0
c$$$      if (oij) then
c$$$         if (ish.lt.jsh) then
c$$$            fock_shell_logic = .false.
c$$$            return
c$$$         else if (ish.eq.jsh) then
c$$$            qq4 = qq4 * 0.5d0
c$$$         endif
c$$$      endif
c$$$      if (okl) then
c$$$         if (ksh.lt.lsh) then
c$$$            fock_shell_logic = .false.
c$$$            return
c$$$         else if (ksh.eq.lsh) then
c$$$            qq4 = qq4 * 0.5d0
c$$$         endif
c$$$      endif
c$$$      if (oikjl) then
c$$$         if (ish.lt.ksh) then
c$$$            fock_shell_logic = .false.
c$$$            return
c$$$         else if (ish.eq.ksh) then
c$$$            if (jsh.lt.lsh) then
c$$$               fock_shell_logic = .false.
c$$$               return
c$$$            else if (jsh.eq.lsh) then
c$$$               qq4 = qq4 * 0.5d0
c$$$            endif
c$$$         endif
c$$$      endif
c$$$c
c$$$      fock_shell_logic = .true.
c$$$c
c$$$      end
#endif
