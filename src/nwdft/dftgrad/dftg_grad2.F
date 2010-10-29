      subroutine dftg_grad2(d_ija, d_kla, frc, g_dens, ldim,
     &                      basis, nproc, nat, scr, lscr, buf, 
     &                      lbuf, tol2e, oactive, nactive,
     &                      geom, s34, oskel)
c    
c     exchange density matrices removed from call list for DFT
c     coulomb-only term.
c
c     Turned into wrapper to call generalized function
c
c     BGJ - 9/98
c    
*
* $Id$
*
      implicit none
c      
      integer g_dens(2), basis, nactive, geom
      integer ldim
      integer nat
      integer nproc
      integer lbuf
      integer lscr
      double precision tol2e, s34(nat)
      logical oskel, oactive(nat)
c
      double precision scr(lscr), buf(lbuf), frc(3, nat)
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
c
c     Dummy arguments
c
      double precision hess
      integer g_dj
c
      call dftg_grad2_gen(d_ija, d_kla, frc, hess, g_dj, g_dens,
     &                    ldim, basis, nproc, nat, scr, lscr, buf, 
     &                    lbuf, tol2e, oactive, nactive,
     &                    geom, s34, oskel, 1)
c
      return
      end
c
c     Generalized from dftg_grad2 to calculate J derivative matrices
c     and J hessian as well as J forces
c
c     BGJ - 9/98
c
      subroutine dftg_grad2_gen(d_ija, d_kla, frc, hess, g_dj, g_dens,
     &                      ldim, basis, nproc, nat, scr, lscr, buf, 
     &                      lbuf, tol2e, oactive, nactive,
     &                      geom, s34, oskel, calc_type)
c    
c     exchange density matrices removed from call list for DFT
c     coulomb-only term.
c    
*
* $Id$
*
      implicit none
#include "errquit.fh"
c      
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "tcgmsg.fh"
#include "sym.fh"
c
      integer g_dens(2), basis, nactive, geom
      integer ldim
      integer nat
      integer nproc
      integer lbuf
      integer lscr
      integer g_dj(3, nat)
      double precision tol2e,  s34(nat), s12, q2, q4, scale
      logical osym, oskel, oschwarz12, diagb, diagk, ijij
      integer calc_type
c
      logical status, oij, okl, oikjl, oactive(nat)
c
      double precision scr(lscr), buf(lbuf), frc(3, nat),
     &     hess(3, nat, 3, nat)
      integer idatom(4)
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
c
c      integer iproc,iat3,iat1,iat2,iat4,ihi4,ijatom,next,
      integer iproc,iat3,iat1,iat2,iat4,ihi4,ijklshell,next,
     &     iab1f,iab1l,iab2f,iab2l,iab3f,iab3l,iab4f,iab4l,
     &     iac2f,iac2l,iac3f,iac3l,iac4f,iac4l,
     &     ish1,ish2,ish3,ish4,iac1l,itop3,itop4,
     &     if1,il1,if2,il2,if3,il3,if4,il4,iac1f,
     &     id1,id2,id3,id4,ig1,im1,ig2,im2,ig3,im3,ig4,im4,
     &     nshbfp,nshbfq,Nintegrals
c
      logical do_force, do_cpks_r, do_hess
      integer l_dfij, k_dfij, l_dfkl, k_dfkl
      integer nxtask
      external nxtask
c
      do_force = calc_type .eq. 1
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_force .or. do_cpks_r .or. do_hess)) then
         call errquit('dft_grad2_gen: illegal calculation type',1,
     &       INPUT_ERR)
      endif
      if (do_cpks_r) then
c
c     Need scratch space for locally formed J mat derivative blocks
c
         if (.not.MA_push_get(MT_DBL, ldim*ldim, 'dfij',
     &        l_dfij, k_dfij)) call
     &        errquit('dftg_grad2_gen: cannot allocate dfij',1, MA_ERR)
         if (.not.MA_push_get(MT_DBL, ldim*ldim, 'dfkl',
     &        l_dfkl, k_dfkl)) call
     &        errquit('dftg_grad2_gen: cannot allocate dfkl',1, MA_ERR)
      endif
c
      iproc = ga_nodeid()
      q4 = 1.0d0                ! If not using symmetry
c
c      ijatom = -1
      ijklshell = -1
      next = nxtask(nproc, 1)
c
      do 95 iat3 = nat, 1, -1
         do iat4 = 1, iat3
            s34(iat4) = schwarz_atom(iat3, iat4) ! Precompute
         enddo
         do 90 iat1 = nat, iat3, -1   ! Want canonical atom list for symmetry
            do 85 iat2 = iat1, 1, -1
c
               s12 = schwarz_atom(iat1, iat2)
               oschwarz12 = (s12*schwarz_max().ge.tol2e)
c     
               if (oschwarz12 .and. oskel) then
                  oschwarz12 = sym_atom_pair(geom, iat1, iat2, q2)
               endif
c     
               if (oschwarz12) then
                  ihi4 = iat3
                  if (iat1.eq.iat3) ihi4 = iat2
                  do 80 iat4 = 1, ihi4
c
c                    include only terms with at least one active center
c
                     if (.not.oactive(iat1) .and. 
     &                   .not.oactive(iat2) .and.
     &                   .not.oactive(iat3) .and. 
     &                   .not.oactive(iat4))goto 1010
c
c                    check screening
c
                     if (s12 * s34(iat4) .lt. tol2e)goto 1010
c
c                    check symmetry
c
                     if (oskel) then
                        osym = sym_atom_quartet(geom,
     &                         iat1, iat2, iat3, iat4, q4)
                        if (.not. osym) goto 1010
                     endif
c
c                    check for functions on each center
c
                     status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0,
     &       BASIS_ERR)
                     status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0,
     &       BASIS_ERR)
                     status = bas_ce2bfr(basis,iat3,iab3f,iab3l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0,
     &       BASIS_ERR)
                     status = bas_ce2bfr(basis,iat4,iab4f,iab4l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0,
     &       BASIS_ERR)
                     if (iab1f.le.0 .or. 
     &                   iab2f.le.0 .or. 
     &                   iab3f.le.0 .or.
     &                   iab4f.le.0)goto 1010
c                              
                     status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0,
     &       BASIS_ERR)
                     status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0,
     &       BASIS_ERR)
                     status = bas_ce2cnr(basis,iat3,iac3f,iac3l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0,
     &       BASIS_ERR)
                     status = bas_ce2cnr(basis,iat4,iac4f,iac4l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0,
     &       BASIS_ERR)
c
c                     ijatom = ijatom + 1
c
c                     if (next .eq. ijatom) then
c
c                       get density matrix blocks
c               
c                       modified for coulomb-only and UHF/RHF.  For exact exchange
c                       need to use get_dens calls (for UHF, 2 calls)
c               
                        call dftg_get_cden (d_ija, d_kla, ldim,
     &                                      iab1f, iab1l, iab2f, 
     &                                      iab2l, iab3f, iab3l, 
     &                                      iab4f, iab4l, g_dens(1))
c    
                        oij = iat1.eq.iat2
                        okl = iat3.eq.iat4
                        oikjl = iat1.eq.iat3 .and. iat2.eq.iat4
                        do 75 ish1 = iac1f, iac1l
                           status = bas_cn2bfr(basis,ish1,if1,il1)
                           ig1 = if1
                           im1 = il1
                           if1 = if1 - iab1f + 1
                           il1 = il1 - iab1f + 1
                           id1 = il1 - if1   + 1
                           if (oij) iac2l = ish1
                           do 70 ish2 = iac2f, iac2l
                              status = bas_cn2bfr(basis,ish2,if2,il2)
                              ig2 = if2
                              im2 = il2
                              if2 = if2 - iab2f + 1
                              il2 = il2 - iab2f + 1
                              id2 = il2 - if2   + 1
                              if (schwarz_shell(ish1,ish2) * 
     &                            schwarz_max()*q4 .ge. tol2e) then
                                 itop3 = iac3l
                                 if (iat1.eq.iat3 .and. 
     &                               iat2.eq.iat4) itop3 = ish1
                                 do 65 ish3 = iac3f, itop3
                                    status = bas_cn2bfr
     &                                       (basis,ish3,if3,il3)
                                    ig3 = if3
                                    im3 = il3
                                    if3 = if3 - iab3f + 1
                                    il3 = il3 - iab3f + 1
                                    id3 = il3 - if3   + 1
                                    itop4 = iac4l
                                    if (okl) itop4 = ish3
                                    if (oikjl .and. 
     &                                 ish1.eq.ish3) itop4 = ish2
                                    do 60 ish4 = iac4f, itop4
                                       if (schwarz_shell(ish1,ish2)*
     &                                     schwarz_shell(ish3,ish4)*q4
     &                                     .ge. tol2e) then
c
                                          ijklshell = ijklshell + 1
c
                                          if (next .eq. ijklshell) then
c
                                          status = bas_cn2bfr
     &                                             (basis,ish4,if4,il4)
                                          ig4 = if4
                                          im4 = il4
                                          if4 = if4 - iab4f + 1
                                          il4 = il4 - iab4f + 1
                                          id4 = il4 - if4   + 1
c
c                                         Determine the number of ERIs in this block of integrals.
c
                                          nshbfp = id1*id2
                                          nshbfq = id3*id4
                                          Nintegrals = nshbfp*nshbfq
c
c                                         2el. integral derivatives
c
                                          if (do_force.or.do_cpks_r)
     &                                       then
                                             call intd_2e4c(basis, ish1, 
     &                                            ish2, basis, ish3,
     &                                            ish4, lscr, scr, lbuf,
     &                                            buf, idatom)
                                          else if (do_hess) then
                                             call intdd_2e4c(basis,
     &                                            ish1, ish2, basis,
     &                                            ish3, ish4, lscr, scr,
     &                                            lbuf, buf, idatom)
                                          endif
c
                                          if (do_force .or. do_hess)
     &                                       then
                                             scale = 1.0d0
                                             if (.not.oij .or. 
     &                                          ish1.ne.ish2) 
     &                                          scale = scale*2.0d0
                                             if (.not.okl .or. 
     &                                          ish3.ne.ish4) 
     &                                          scale = scale*2.0d0
                                             if (.not.oikjl .or. 
     &                                          ish1.ne.ish3 .or. 
     &                                          ish2.ne.ish4) 
     &                                          scale = scale*2.0d0
                                          else if (do_cpks_r) then
                                             diagb = oij .and. 
     &                                          ish1.eq.ish2
                                             diagk = okl .and. 
     &                                          ish3.eq.ish4 
                                             ijij = oikjl .and.
     &                                          ish1.eq.ish3 .and.
     &                                          ish2.eq.ish4
                                          endif
c
c - indent - region redefine
c---------------------------------------------------------------------
c
      if (do_force) then
c
         call dint_block_jfrc(d_ija, d_kla, ldim, frc, buf, q4*scale,
     &        idatom, if1, il1, if2, il2, if3, il3, if4, il4)
c
      else if (do_cpks_r) then
c
         call dint_block_jdmat(d_ija(if1,if2), d_kla(if3,if4), ldim,
     &        g_dj, dbl_mb(k_dfij), dbl_mb(k_dfkl), buf, q4,
     &        diagb, diagk, ijij, idatom, ig1, im1, ig2, im2,
     &        ig3, im3, ig4, im4)
c
      else if (do_hess) then
c
         call d2int_block_jhess(d_ija, d_kla, ldim, hess, nat, buf,
     &        q4*scale, idatom, if1, il1, if2, il2, if3, il3, if4, il4)
c
      endif
c
      next = nxtask(nproc, 1)
c
      endif     ! next
c
c---------------------------------------------------------------------
c - indent - region pick back up
c
                                       endif   ! schwarz(1234)
   60                               continue   ! ish4
   65                            continue      ! ish3
                              endif            ! schwarz(12)
   70                      continue            ! ish2
   75                   continue               ! ish1
c
c                        next = nxtask(nproc, 1)
c
c                     endif                     ! next
c
 1010                continue                  ! jump here for inactive/zero
c
   80             continue                     ! iat4
               endif                           ! schwarz(12)
   85       continue                           ! iat2
   90    continue                              ! iat1
   95 continue                                 ! iat3
      next = nxtask(-nproc, 1)
c
      if (do_cpks_r) then
         if (.not.MA_pop_stack(l_dfkl))
     &        call errquit('dftg_grad2_gen: cannot pop dfkl',1, MA_ERR)
         if (.not.MA_pop_stack(l_dfij))
     &        call errquit('dftg_grad2_gen: cannot pop dfij',1, MA_ERR)
      endif
c
      return
      end
