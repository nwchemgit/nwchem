c
c $Id$
c
      subroutine oned_calc(nxyz,g_rhs,Ibuf,lbuf,Iscr,lscr,
     *      skelsym,doS,doT,doV,geom,basis,nat,ecp_channel)
c
      implicit none
c
#include "bas.fh"
#include "global.fh"
#include "sym.fh"
c
      integer nat                  ! [input] the number of atoms
      integer nxyz                 ! [input] 3*number of atoms
      integer g_rhs                ! [output] handles for the RHSs
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
      integer task_size, ijatom, next, iat1, iat2
      integer iac1f, iac1l, iac2f, iac2l, idatom(2)
      integer iab1f, iab1l, iab2f, iab2l, ish1, ish2, if1, il1, if2, il2
      integer ic, icart1, i
      integer nproc, atx, nint
      integer lo(3), hi(3), ld(2)
      logical status
      double precision qfac
c
      integer nxtask
      external nxtask
cc AJL/Begin/SPIN ECPs
      integer ecp_channel
cc AJL/End
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
                status = bas_cn2bfr(basis,ish2,if2,il2)

                nint = (il1 - if1 + 1) * (il2 - if2 + 1)
c               write(6,*) 'nint is ',nint

                lo(2) = if2
                hi(2) = il2
                lo(3) = if1
                hi(3) = il1
                ld(1) = 1
                ld(2) = il2 - if2 + 1
c
C     dS
c
                if (doS) then
                call intd_1eov(basis,ish1,basis,ish2,lscr,Iscr,
     &               lbuf,Ibuf,idatom)

                if ( idatom(1) .ge. 1 ) then
c               write(6,*) 'integral buffer for ',idatom(1),
c    *              idatom(2)
c
                  ic = 1
                  do i = 1, 2
                  do icart1 = 1, 3 
                    lo(1) = (idatom(i)-1) * 3 + icart1
                    hi(1) = (idatom(i)-1) * 3 + icart1
cedo                    call nga_acc(g_rhs, lo, hi, Ibuf(ic), ld, 1.0d0)
                    call nga_put(g_rhs, lo, hi, Ibuf(ic), ld)
                    ic = ic + nint
                  enddo
                  enddo
                endif   ! if idatom
                endif   ! if doS
c
C     dH
c
                if (doV.and.doT) then
c
cc AJL/Begin/SPIN ECPs
c
                 if (ecp_channel.eq.1) then
                   call intd_1eh1(basis,ish1,basis,ish2,lscr,Iscr,
     *                  lbuf,Ibuf)
                 else
                   call intd_1eh1_beta(basis,ish1,basis,ish2,lscr,Iscr,
     *                  lbuf,Ibuf)
                 end if
c
cc AJL/End
c
                 ic = 1
                 do atx = 1, nxyz
                   lo(1) = atx
                   hi(1) = atx
!                   call nga_acc(g_rhs, lo, hi, Ibuf(ic), ld, 1.0d0)
                   call nga_put(g_rhs, lo, hi, Ibuf(ic), ld)
                   ic = ic + nint
                 enddo
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
      return
      end
