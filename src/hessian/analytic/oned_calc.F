c
c $Id: oned_calc.F,v 1.2 2000-04-26 21:49:16 windus Exp $
c
      subroutine oned_calc(nxyz,g_rhs,Ibuf,lbuf,Iscr,lscr,
     *      skelsym,doS,doT,doV,geom,basis,nat)
c
      implicit none
c
#include "bas.fh"
#include "global.fh"
#include "sym.fh"
c
      integer nat                  ! [input] the number of atoms
      integer nxyz                 ! [input] 3*number of atoms
      integer g_rhs(nxyz)          ! [output] handles for the RHSs
      integer lbuf                 ! [input] length of the integral buffer
      double precision Ibuf(lbuf)  ! [input] Buffer for integrals
      integer lscr                 ! [input] length of the scratch buffer
      double precision Iscr(lscr)  ! [input] scratch for integral calculation
      logical skelsym              ! [input] skeleton symmetry?
      logical doS                  ! [input] do overlap contribution?
      logical doT                  ! [input] do kinetic contribution?
      logical doV                  ! [input] do potential contribution?
      integer geom                 ! [input] the current geometry handle
      integer basis                ! [input] the current basis handle
c
      integer task_size, ijatom, next, iat1, iat2, attemp1
      integer iac1f, iac1l, iac2f, iac2l, idatom(2)
      integer iab1f, iab1l, iab2f, iab2l, ish1, ish2, if1, il1, if2, il2
      integer ic, icart1, icart2, ip1, ip2, i, j, at1, at2
      integer nproc, atx, nint
      logical status
      double precision qfac, de
c
      integer nxtask
      external nxtask
c
      nproc = ga_nnodes()
      task_size = 1
      ijatom = -1
      next = nxtask(nproc,task_size)
c
      do 90 iat1 = 1, nat
        do 80 iat2 = 1, nat
c       do 80 iat2 = 1, iat1

          ijatom = ijatom + 1
          if ( ijatom .eq. next ) then

            status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
            status = bas_ce2bfr(basis,iat2,iab2f,iab2l)

            if (iab1f.le.0 .or. iab2f.le.0) then
c
c     At least one center has no functions on it ... next atom
c
              goto 1010
            endif

            if (skelsym) then
              if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
     *              goto 1010
            else
              qfac = 1.0d0
            endif

            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)

            do 70 ish1 = iac1f, iac1l
c             if ( iat1.eq.iat2 ) iac2l = ish1
              do 60 ish2 = iac2f, iac2l

c
c               shell block in atomic (D/Dw)-matrix block
c
                status = bas_cn2bfr(basis,ish1,if1,il1)
c               if1 = if1 - iab1f + 1
c               il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
c               if2 = if2 - iab2f + 1
c               il2 = il2 - iab2f + 1

                nint = (il1 - if1 + 1) * (il2 - if2 + 1)
                write(6,*) 'nint is ',nint
c
C     dS
c
                if (doS) then
                call intd_1eov(basis,ish1,basis,ish2,lscr,Iscr,
     &               lbuf,Ibuf,idatom)

                if ( idatom(1) .ge. 1 ) then
c               idatom(1).ge.0 <=> idatom(2).ge.0 (no check necessary)
                  qfac = qfac*2
                  call dscal(lbuf,qfac,Ibuf,1)
                  qfac = qfac/2
                  ic = 1
                  do 27 icart1 = 1, 3 
                    call ga_acc(g_rhs(idatom(1)),if1,il1,if2,il2,
     *               Ibuf(ic),il1-if1+1,1.0)
                    call ga_acc(g_rhs(idatom(2)),if1,il1,if2,il2,
     *               Ibuf(ic),il1-if1+1,-1.0)
c put in the ga_acc for idatoms 1 and 2
 27               continue
                endif   ! if idatom
                endif   ! if doS
c
C     dH
c
                if (doV.and.doT) then
                 call intd_1eh1(basis,ish1,basis,ish2,lscr,Iscr,
     *                lbuf,Ibuf)
c                if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) qfac=qfac*2
c                call dscal(lbuf,qfac,Ibuf,1)
c                if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) qfac=qfac/2
                 ic = 1
                 do atx = 1, nxyz
                   write(6,*) 'if1,il1,if2,il2, ',if1,il1,if2,il2
                   write(6,*) 'iab1f, iab1l, iab2f, iab2l', iab1f,
     *              iab1l, iab2f, iab2l
                   write(6,*) 'Ibuf(ic) ',Ibuf(ic), ic
c                  call ga_acc(g_rhs(atx),iab1f,iab1l,iab2f,iab2l,
c    *              Ibuf(ic),iab1l-iab1f+1,1.0d0)
                   call ga_acc(g_rhs(atx),if1,il1,if2,il2,
     *              Ibuf(ic),il1-if1+1,1.0d0)
c                  if (.not.((iab1f.eq.iab2f).and.
c    *               (iab1l.eq.iab2l))) call ga_acc(g_rhs(atx),
c    *               iab2f,iab2l,iab1f,iab1l,Ibuf(ic),nint,1.0d0)
                    write(6,*) 'for g_rhs(atx) ',atx
                    call ga_print(g_rhs(atx))
                    ic = ic + nint
c put in the ga_accs
                 enddo
c  first print out the rhs values (should be all zeroes)!
c     write(6,*) 'g_rhs in oned_calc'
c     do i = 1, nxyz
c       call ga_print(g_rhs(i))
c     enddo 
                endif   ! doV and doT

 60           continue
 70         continue
 1010       continue

            next = nxtask(nproc,task_size)
          endif  ! if my task

 80     continue
 90   continue
      next = nxtask(-nproc,task_size)
c
      end
