      subroutine dftg_grad2(d_ija, d_kla, frc, g_dens, ldim,
     &                      basis, nproc, nat, scr, lscr, buf, 
     &                      lbuf, tol2e, oactive, nactive,
     &                      geom, s34, oskel)
c    
c     exchange density matrices removed from call list for DFT
c     coulomb-only term.
c    
*
* $Id: dftg_grad2.F,v 1.6 1997-11-04 10:12:15 d3e129 Exp $
*
      implicit none
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
      double precision tol2e,  s34(nat), s12, q2, q4, scale
      logical osym, oskel, oadapt, oschwarz12
c
      logical status, oij, okl, oikjl, oactive(nat)
c
      double precision scr(lscr), buf(lbuf), frc(3, nat)
      integer idatom(4)
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
c
      integer iproc,iat3,iat1,iat2,iat4,ihi4,ijatom,next,
     &     iab1f,iab1l,iab2f,iab2l,iab3f,iab3l,iab4f,iab4l,
     &     iac2f,iac2l,iac3f,iac3l,iac4f,iac4l,
     &     ish1,ish2,ish3,ish4,iac1l,itop3,itop4,
     &     if1,il1,if2,il2,if3,il3,if4,il4,iac1f,
     &     id1,id2,id3,id4,
     &     nshbfp,nshbfq,Nintegrals,icount,i,icart,ip1,ip2,ip3,ip4
      double precision dE,den2
c
      iproc = ga_nodeid()
      q4 = 1.0d0                ! If not using symmetry
c
      ijatom = -1
      next = nxtval(nproc)
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
     &                  ('dftg_grad2: error with bas_ce2bfr',0)
                     status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0)
                     status = bas_ce2bfr(basis,iat3,iab3f,iab3l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0)
                     status = bas_ce2bfr(basis,iat4,iab4f,iab4l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2bfr',0)
                     if (iab1f.le.0 .or. 
     &                   iab2f.le.0 .or. 
     &                   iab3f.le.0 .or.
     &                   iab4f.le.0)goto 1010
c                              
                     status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0)
                     status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0)
                     status = bas_ce2cnr(basis,iat3,iac3f,iac3l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0)
                     status = bas_ce2cnr(basis,iat4,iac4f,iac4l)
                     if (.not.status)call errquit 
     &                  ('dftg_grad2: error with bas_ce2cnr',0)
c
                     ijatom = ijatom + 1
c
                     if (next .eq. ijatom) then
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
                           if1 = if1 - iab1f + 1
                           il1 = il1 - iab1f + 1
                           id1 = il1 - if1   + 1
                           if (oij) iac2l = ish1
                           do 70 ish2 = iac2f, iac2l
                              status = bas_cn2bfr(basis,ish2,if2,il2)
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
                                          status = bas_cn2bfr
     &                                             (basis,ish4,if4,il4)
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
                                          call intd_2e4c(basis, ish1, 
     &                                         ish2, basis, ish3, ish4,
     &                                         lscr, scr, lbuf, buf, 
     &                                         idatom)
c
                                          scale = 1.0d0
                                          if (.not.oij .or. 
     &                                       ish1.ne.ish2) 
     &                                       scale = scale*2.0d0
                                          if (.not.okl .or. 
     &                                       ish3.ne.ish4) 
     &                                       scale = scale*2.0d0
                                          if (.not.oikjl .or. 
     &                                       ish1.ne.ish3 .or. 
     &                                       ish2.ne.ish4) 
     &                                       scale = scale*2.0d0
c
c - indent - region redefine
c---------------------------------------------------------------------
c
      icount = 1
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               dE = 0.D0
               do ip1 = if1, il1
                  do ip2 = if2, il2
                     do ip3 = if3, il3
                        do ip4 = if4, il4
c                   
c                          hacked for DFT; exchange contributions deleted
c                   
c                          den2 = 0.5D0   * (d_ij(ip1,ip2) * d_kl(ip3,ip4))
c                   &           - 0.125D0 * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
c                   &                     +  d_il(ip1,ip4) * d_jk(ip2,ip3))
c
                           den2 = 0.5D0 *
     &                            (d_ija(ip1,ip2) * d_kla(ip3,ip4))
c
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
c
                        enddo 
                     enddo
                  enddo
               enddo
c
               frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     &                                dE * q4 * scale
c
            enddo
         else
            icount = icount + 3 * Nintegrals
         endif
      enddo 
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
                        next = nxtval(nproc)
c
                     endif                     ! next
c
 1010                continue                  ! jump here for inactive/zero
c
   80             continue                     ! iat4
               endif                           ! schwarz(12)
   85       continue                           ! iat2
   90    continue                              ! iat1
   95 continue                                 ! iat3
      next = nxtval(-nproc)
      return
      end
