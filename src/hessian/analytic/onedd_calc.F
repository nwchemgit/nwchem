c
c $Id$
c
      subroutine onedd_calc(Hxy,nxyz,Ibuf,lbuf,Iscr,lscr,
     *      skelsym,doS,doT,doV,doECP,
     *      g_wdens,g_dens,geom,basis,nat,ndens,dens,wdens,matb,
     *      ecp_channels)
c
      implicit none
c
#include "bas.fh"
#include "global.fh"
#include "sym.fh"
c
      integer nat                  ! [input] the number of atoms
      double precision Hxy(3,nat,3,nat) ! [output] the hessian matrix
      integer nxyz                 ! [input] 3*number of atoms
      integer lbuf                 ! [input] length of the integral buffer
      double precision Ibuf(lbuf)  ! [input] Buffer for integrals
      integer lscr                 ! [input] length of the scratch buffer
      double precision Iscr(lscr)  ! [input] scratch for integral calculation
      logical skelsym              ! [input] skeleton symmetry?
      logical doS                  ! [input] do overlap contribution?
      logical doT                  ! [input] do kinetic contribution?
      logical doV                  ! [input] do potential contribution?
      logical doECP                ! [input] do ECP contribution?
      integer ndens                ! [input] number of density matrices
      integer g_wdens              ! [input] the weighted density handle
      integer g_dens(ndens)        ! [input] the density handle
      integer geom                 ! [input] the current geometry handle
      integer basis                ! [input] the current basis handle
      integer matb                 ! [input] max. atom buffer space needed
      double precision dens(matb,matb,ndens)  ! [input] scratch for density blocks
      double precision wdens(matb,matb) ! [input] scratch for wdensity blocks
      
c
      integer task_size, ijatom, next, iat1, iat2, attemp1
      integer iac1f, iac1l, iac2f, iac2l, idatom(2)
      integer iab1f, iab1l, iab2f, iab2l, ish1, ish2, if1, il1, if2, il2
      integer ic, icart1, icart2, ip1, ip2, i, j, at1, at2
      integer nproc
      logical status
      double precision qfac, de
c
      integer nxtask
      external nxtask
c
cc AJL/Begin/SPIN ECPs
      integer ecp_channels
cc Buffer for Beta Channel ECP integrals
      double precision Ibuf_beta(lbuf)
cc AJL/End
c
      nproc = ga_nnodes()
      task_size = 1
      ijatom = -1
      next = nxtask(nproc,task_size)
c
      do 90 iat1 = 1, nat
        do 80 iat2 = 1, iat1

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

            call ga_get (g_dens(1), iab1f,iab1l,iab2f,iab2l,
     *                   dens,matb)
            if (ndens.eq.2) then
               call ga_get (g_dens(2), iab1f,iab1l,iab2f,iab2l,
     *                      dens(1,1,2),matb)
c
cc AJL/Begin/SPIN ECPs
cc We only want to sum this if ecp_channels are one
cc Otherwise we are going to apply different ddH to each density channel
               if (ecp_channels.eq.1) then
                 call daxpy(matb*matb,1.0d+00,dens(1,1,2),1,dens,1)
               endif
cc AJL/End
c
            endif
            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,matb)

            do 70 ish1 = iac1f, iac1l
              if ( iat1.eq.iat2 ) iac2l = ish1
              do 60 ish2 = iac2f, iac2l
c               write(6,*) 'ish1,ish2 ',ish1,ish2
c
c               shell block in atomic (D/Dw)-matrix block
c
                status = bas_cn2bfr(basis,ish1,if1,il1)
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1

C               overlap derivatives

                if (doS) then
                call intdd_1eov(basis,ish1,basis,ish2,lscr,Iscr,
     &               lbuf,Ibuf,idatom)

C     w x ddS

                if ( idatom(1) .ge. 1 ) then
C               idatom(1).ge.0 <=> idatom(2).ge.0 (no check necessary)
c                 at1 = idatom(1)
c                 at2 = idatom(2)
                  ic = 1
                  do 29 i = 1, 2
                  do 28 j = 1, 2   ! for both atoms
                    at1 = idatom(i) ! switch atom numbers for the second half
                    at2 = idatom(j)
                    do 27 icart2 = 1, 3 
                    do 26 icart1 = 1, 3
                      de = 0.D0           ! contract wdensity and integrals
                      do 22 ip1 = if1, il1
                        do 20 ip2 = if2, il2
                          dE = dE - wdens(ip1,ip2) * Ibuf(ic)
                          ic = ic + 1
 20                     continue
 22                   continue
                      dE = dE * qfac
                      dE = dE + dE
                      Hxy(icart1,at1,icart2,at2) = 
     *                 Hxy(icart1,at1,icart2,at2) + dE
 26                 continue
 27                 continue
 28               continue
 29               continue
                endif   ! if idatom
                endif   ! if doS
c
C               1el. derivatives
c
                if (doV.and.doT) then
                call intdd_1eh1(basis,ish1,basis,ish2,lscr,Iscr,
     *               lbuf,Ibuf)
c
C     D x ddH
c
c
c    Handle (i|d2V|j) first.  
c
                ic = 1
                do at1 = 1, nat
                  do icart2 = 1, 3
                    do icart1 = 1, 3
                      dE = 0.0D0
                      do ip1 = if1, il1
                        do ip2 = if2, il2
                          dE = dE + dens(ip1,ip2,1)*Ibuf(ic)
c
cc AJL/Begin/SPIN ECPs
cc We haven't summed the local density so need to add this
                          if (ecp_channels.gt.1) then
                            dE = dE + dens(ip1,ip2,2)*Ibuf(ic)
                          end if
cc AJL/End
c
                          ic = ic + 1
                        enddo
                      enddo
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                      dE = dE * qfac
                      Hxy(icart1,at1,icart2,at1) =
     *                 Hxy(icart1,at1,icart2,at1) + dE
                    enddo
                  enddo
                enddo
c
c    Handle (di|dV|j), (i|dV|dj) 
c
                do 50 attemp1 = 1, 2
                if (attemp1.eq.1) then
                  at1 = iat1
                else
                  at1 = iat2
                endif 
                do 50 at2 = 1, nat
                  do 40 icart2 = 1, 3
                  do 40 icart1 = 1, 3
                    dE = 0.D0
                    do 31 ip1 = if1, il1
                      do 30 ip2 = if2, il2
                        dE = dE + dens(ip1,ip2,1) * Ibuf(ic)
c
cc AJL/Begin/SPIN ECPs
cc We haven't summed the local density so need to add this
                        if (ecp_channels.gt.1) then
                          dE = dE + dens(ip1,ip2,2)*Ibuf(ic)
                        end if
cc AJL/End
c
c                       write(6,*) 'dens and int ',dens(ip1,ip2,1),
c    *                     Ibuf(ic)
                        ic = ic + 1
 30                   continue
 31                 continue
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                    dE = dE * qfac
c                   write(6,*) 'at1,at2,icart1,icart2,dE ',at1,at2,
c    *               icart1,icart2,dE
                    Hxy(icart1,at1,icart2,at2) = 
     *               Hxy(icart1,at1,icart2,at2) + dE
                    Hxy(icart2,at2,icart1,at1) = 
     *               Hxy(icart2,at2,icart1,at1) + dE
 40               continue
 50             continue
c
c do (di2|V+T|j), (di|V+T|dj), and (i|V+T|dj2)
c
                do i = 1, 2
                do j = 1, i
                  if (i.eq.1) then
                    at2 = iat1
                  else
                    at2 = iat2
                  endif
                  if (j.eq.1) then
                    at1 = iat1
                  else
                    at1 = iat2
                  endif
                  do icart2 = 1, 3
                  do icart1 = 1, 3
                    dE = 0.D0
                    do ip1 = if1, il1
                      do ip2 = if2, il2
                        dE = dE + dens(ip1,ip2,1) * Ibuf(ic)
c
cc AJL/Begin/SPIN ECPs
cc We haven't summed the local density so need to add this
                        if (ecp_channels.gt.1) then
                          dE = dE + dens(ip1,ip2,2)*Ibuf(ic)
                        end if
cc AJL/End
c
                        ic = ic + 1
                      enddo
                    enddo
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE+dE
                    dE = dE * qfac
                    if (i.eq.j) then
                      Hxy(icart1,at1,icart2,at2) =
     *                 Hxy(icart1,at1,icart2,at2) + dE
                    else
                      Hxy(icart1,at1,icart2,at2) =
     *                 Hxy(icart1,at1,icart2,at2) + dE
                      Hxy(icart2,at2,icart1,at1) = 
     *                 Hxy(icart2,at2,icart1,at1) + dE
                    endif
                  enddo
                  enddo  
                enddo
                enddo
                endif   ! doV and doT
c
c If there are ECPs we need the derivatives for those too
c
                if (doECP) then
                call intdd_1eecp1(basis,ish1,basis,ish2,
     &               lscr,Iscr,
     &               lbuf,Ibuf)
c
cc AJL/Begin/SPIN ECPs
                if (ecp_channels.gt.1) then
                  call intdd_1eecp1_beta(basis,ish1,basis,ish2,
     &                 lscr,Iscr,
     &                 lbuf,Ibuf_beta)
                end if
cc AJL/End
c
                ic = 0
                do at2 = 1, nat
                do at1 = 1, at2
                  do icart2 = 1, 3
                  do icart1 = 1, 3
                    dE = 0.D0
c                   do ip2 = if2, il2
c                     do ip1 = if1, il1
                    do ip1 = if1, il1
                      do ip2 = if2, il2
                        ic = ic + 1
                        dE = dE + dens(ip1,ip2,1) * Ibuf(ic)
c
cc AJL/Begin/SPIN ECPs
cc Now we compute the Hecp_beta x Dens_beta
                        if (ecp_channels.gt.1) then
                          dE = dE + dens(ip1,ip2,2) * Ibuf_beta(ic)
                        end if
cc AJL/End
c
                      enddo
                    enddo
                    if (iat1.ne.iat2 .or. ish1.ne.ish2) dE = dE + dE
                    dE = dE * qfac
                    if (at1.eq.at2) then
                      Hxy(icart1,at1,icart2,at2) =
     *                 Hxy(icart1,at1,icart2,at2) + dE 
                    else 
                      Hxy(icart1,at1,icart2,at2) =
     *                 Hxy(icart1,at1,icart2,at2) + dE 
                      Hxy(icart2,at2,icart1,at1) =
     *                 Hxy(icart2,at2,icart1,at1) + dE 
                    endif
                  enddo
                  enddo
                enddo
                enddo
                endif

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
