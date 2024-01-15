      subroutine grad1 ( H, lbuf, scr, lscr, dens, wdens, frc_nuc,
     $                   frc_kin, frc_wgh, g_force,
     $                   g_dens, g_wdens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb, oskel )
c$Id: grad1.F,v 1.17 2008/04/22 22:00:19 niri Exp $

C     one electron contribution to RHF, ROHF and UHF gradients
C     now also UMP2

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "bq_params.fh"

C-------------------------parameters--------------------------------
      integer lbuf, lscr,
     $     g_dens,        ! density matrix (summed if ROHF, UHF)
     $     g_wdens,       ! weighted density (Lagrangian)
     $     g_force,       ! global force array
     $     basis, geom, nproc, nat, max_at_bf, rtdb

      double precision H, ! integral derivatives
     $     scr, 
     $     dens,          ! local density block
     $     wdens,         ! local weighted density block
     $     frc_nuc, frc_kin, frc_wgh   ! forces arrays

      dimension H ( lbuf ), frc_nuc(3, nat), frc_kin(3, nat),
     $          frc_wgh(3, nat), scr(lscr),
     $          dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf)

      logical oskel   ! symmetry?

C-------------------------local variables--------------------------

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2

      double precision crd1, crd2 ! atomic coordinates
      dimension crd1(3), crd2(3)

      integer idatom
      dimension idatom(2)

      double precision dE, dx, dy, dz, qfac, fact, q1, q2

      logical status, pointforce,dobq

      character*16 name

      integer bq_ncent
      integer i_qbq,i_cbq
      double precision r12
c
      integer nxtask, task_size
      external nxtask

      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes

      pointforce = geom_include_bqbq(geom)
      dobq = geom_extbq_on()
      call hf_print_set(1)
cgk debug
*      write(*,*)'gk: in grad1, dobq = ', dobq
*      write(*,*)'gk: in grad1, pointforce = ', pointforce
cgk end

      ijatom = -1
      next = nxtask(nproc,task_size)
      do 90, iat1 = 1, nat
        do 80, iat2 = 1, iat1

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

            if (oskel) then
               if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
     $              goto 1010
            else
               qfac = 1.0d0
            endif

            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)

            call ga_get (g_dens, iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)

            do 70, ish1 = iac1f, iac1l
              if ( iat1.eq.iat2 ) iac2l = ish1
              do 60, ish2 = iac2f, iac2l

C               shell block in atomic (D/Dw)-matrix block
                status = bas_cn2bfr(basis,ish1,if1,il1)
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1

                nint = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )

C               overlap derivatives
                call intd_1eov(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H,idatom)

C     Dw x S

                if ( idatom(1) .ge. 1 ) then
C               idatom(1).ge.0 <=> idatom(2).ge.0 (no check necessary)
                  ic = 1
                  do 28, icart = 1, 3
                    de = 0.D0
                    do 22, ip1 = if1, il1
                      do 20, ip2 = if2, il2
                        dE = dE + wdens(ip1,ip2) * H(ic)
                        ic = ic + 1
 20                   continue
 22                 continue
                    dE = dE * qfac
                    frc_wgh(icart,idatom(1)) = frc_wgh(icart,idatom(1)) 
     $                                      - dE - dE
                    frc_wgh(icart,idatom(2)) = frc_wgh(icart,idatom(2)) 
     $                                      + dE + dE
 28               continue
                endif

C               1el. derivatives
                if(.not.dobq) then
                call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)
                else
                call intd_1epot(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)
                end if
                  

C     D x H

                ic=1
                do 50, iat3 = 1, nat
                  do 40, icart = 1, 3
                    dE = 0.D0
                    do 31, ip1 = if1, il1
                      do 30, ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
                        ic = ic + 1
 30                   continue
 31                 continue
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                    dE = dE * qfac
                    frc_kin(icart,iat3) = frc_kin(icart,iat3) + dE
 40               continue
 50             continue

 60           continue
 70         continue

 1010       continue

C     Vnn 

            if ( iat1 .NE. iat2 ) then
              if (iab1f.ne.0 .or. iab2f.ne.0 .or. pointforce ) then
C               no forces between point charges (for John Nicholas)
                status = geom_cent_get (geom, iat1, name, crd1, q1)
                status = geom_cent_get (geom, iat2, name, crd2, q2)
                dx = crd2(1) - crd1(1)
                dy = crd2(2) - crd1(2)
                dz = crd2(3) - crd1(3)
                fact = q1 * q2 / SQRT ( dx*dx + dy*dy + dz*dz ) **3
                dE = dx * fact
                frc_nuc(1,iat1) = frc_nuc(1,iat1) + dE
                frc_nuc(1,iat2) = frc_nuc(1,iat2) - dE
                dE = dy * fact
                frc_nuc(2,iat1) = frc_nuc(2,iat1) + dE
                frc_nuc(2,iat2) = frc_nuc(2,iat2) - dE
                dE = dz * fact
                frc_nuc(3,iat1) = frc_nuc(3,iat1) + dE
                frc_nuc(3,iat2) = frc_nuc(3,iat2) - dE
              endif
            endif

            next = nxtask(nproc,task_size)
          endif

 80     continue
 90   continue
      next = nxtask(-nproc,task_size)

c     
c     continue with forces if there are external charges (MV)
c     --------------------------------------------------
      if(.not.dobq) return
cgk debug
*      write(*,*)'gk:  continuing with forces for external charges'
cgk end

      bq_ncent = geom_extbq_ncenter()
      i_cbq = geom_extbq_coord()
      i_qbq = geom_extbq_charge()

      ijatom = -1
      next = nxtask(nproc,task_size)
      do 91, iat1 = 1, nat
        do 81, iat2 = 1, bq_ncent

          ijatom = ijatom + 1
          if ( ijatom .eq. next ) then
            status = geom_cent_get (geom, iat1, name, crd1, q1)

            q2      = dbl_mb(i_qbq+iat2-1)
            crd2(1) = dbl_mb(i_cbq+(iat2-1)*3) 
            crd2(2) = dbl_mb(i_cbq+(iat2-1)*3+1)
            crd2(3) = dbl_mb(i_cbq+(iat2-1)*3+2)
            
            dx = crd2(1) - crd1(1)
            dy = crd2(2) - crd1(2)
            dz = crd2(3) - crd1(3)

c           == check for small distances (the regions could overlap) ==
            r12 = sqrt(dx*dx + dy*dy + dz*dz) 
            if (r12.gt.bq_smalldist) then
              fact = q1 * q2 / SQRT ( dx*dx + dy*dy + dz*dz ) **3
              dE = dx * fact
              frc_nuc(1,iat1) = frc_nuc(1,iat1) + dE
              dE = dy * fact
              frc_nuc(2,iat1) = frc_nuc(2,iat1) + dE
              dE = dz * fact
              frc_nuc(3,iat1) = frc_nuc(3,iat1) + dE
            end if
c
            next = nxtask(nproc,task_size)
          endif
 81     continue
 91   continue
      next = nxtask(-nproc,task_size)

cgk debug
*      write(*,*)'gk: leaving grad1'
cgk end

      return
      end
      subroutine grad_hnd_cos ( H, lbuf, scr, lscr, 
     $                   dens, wdens, frc_nuc,
     $                   frc_kin, frc_wgh, g_force,
     $                   g_dens, g_wdens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb, oskel )
c$Id: grad1.F,v 1.17 2008/04/22 22:00:19 niri Exp $

C     one electron contribution to RHF, ROHF and UHF gradients
C     now also UMP2

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "stdio.fh"

C-------------------------parameters--------------------------------
      integer lbuf, lscr,
     $     g_dens,        ! density matrix (summed if ROHF, UHF)
     $     g_wdens,       ! weighted density (Lagrangian)
     $     g_force,       ! global force array
     $     basis, geom, nproc, nat, max_at_bf, rtdb

      double precision H, ! integral derivatives
     $     scr, 
     $     dens,          ! local density block
     $     wdens,         ! local weighted density block
     $     frc_nuc, frc_kin, frc_wgh   ! forces arrays

      dimension H ( lbuf ), frc_nuc(3, nat), frc_kin(3, nat),
     $          frc_wgh(3, nat), scr(lscr),
     $          dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf)

      logical oskel   ! symmetry?

C-------------------------local variables--------------------------

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2

      double precision dE, qfac

      logical status, pointforce

      integer nxtask, task_size
      external nxtask

c     ---- -cosmo- gradient term -----
      logical odbug

      odbug=.true.
      if(odbug) then
         write(Luout,*) 'in -grad1_hnd_cos- ...'
      endif
c

      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes

      pointforce = geom_include_bqbq(geom)

      call hf_print_set(1)

      ijatom = -1
      next = nxtask(nproc,task_size)
      do 90, iat1 = 1, nat
        do 80, iat2 = 1, iat1

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

            if (oskel) then
               if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
     $              goto 1010
            else
               qfac = 1.0d0
            endif

            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)

            call ga_get (g_dens, iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)

            do 70, ish1 = iac1f, iac1l
              if ( iat1.eq.iat2 ) iac2l = ish1
              do 60, ish2 = iac2f, iac2l

C               shell block in atomic (D/Dw)-matrix block
                status = bas_cn2bfr(basis,ish1,if1,il1)
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1

                nint = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )

                ic=1
                do iat3 = 1, nat
                  do icart = 1, 3
                    do ip1 = if1, il1
                      do ip2 = if2, il2
                        H(ic)=0.D0
                        ic = ic + 1
                      enddo
                    enddo
                  enddo
                enddo

C               1el. -cosmo- derivatives
c-              call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
c-   &               lbuf,H)

C     D x H

                ic=1
                do 50, iat3 = 1, nat
                  do 40, icart = 1, 3
                    dE = 0.D0
                    do 31, ip1 = if1, il1
                      do 30, ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
                        ic = ic + 1
 30                   continue
 31                 continue
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                    dE = dE * qfac
                    frc_kin(icart,iat3) = frc_kin(icart,iat3) + dE
 40               continue
 50             continue

 60           continue
 70         continue

 1010       continue

            next = nxtask(nproc,task_size)
          endif

 80     continue
 90   continue
      next = nxtask(-nproc,task_size)

      return
      end
