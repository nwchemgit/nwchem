      subroutine grad1 ( H, lbuf, scr, lscr, dens, wdens, frc, 
     $                   g_dens, g_wdens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb )

      implicit real*8 (a-h,o-z)

#include "mafdecls.h"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"

      dimension frc ( 3, nat ), idatom(2), scr(lscr)
      dimension H ( lbuf ), 
     $          dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf),
     $          crd1(3), crd2(3)

      integer g_dens, g_wdens, basis, geom
      integer rtdb

      logical status, pointforce, point

      character*16 name

      iproc = nodeid()
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      pointforce = rtdb_get(rtdb,'gradients:point', MT_LOG, 1,point) 
      pointforce = (pointforce .and. point)

      call hf_print_set(1)

      ijatom = -1
      next = nxtval(nproc)
      do 90, iat1 = 1, nat
        do 80, iat2 = 1, iat1

          ijatom = ijatom + 1
C         if (mod(ijatom, nproc) .eq. iproc) then
          if ( ijatom .eq. next ) then

            status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
            status = bas_ce2bfr(basis,iat2,iab2f,iab2l)

            if (iab1f.le.0 .or. iab2f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c
              goto 1010
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
C               idatom(1).ge.0 <=> idatom(1).ge.0 (no check necessary)
                  ic = 1
                  do 28, icart = 1, 3
                    de = 0.D0
                    do 22, ip1 = if1, il1
                      do 20, ip2 = if2, il2
                        dE = dE + wdens(ip1,ip2) * H(ic)
                        ic = ic + 1
 20                   continue
 22                 continue
                    frc(icart,idatom(1)) = frc(icart,idatom(1)) 
     $                                     - dE - dE
                    frc(icart,idatom(2)) = frc(icart,idatom(2)) 
     $                                     + dE + dE
 28               continue
                endif

C               1el. derivatives
                call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)

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
                    frc(icart,iat3) = frc(icart,iat3) + dE
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
                frc(1,iat1) = frc(1,iat1) + dE
                frc(1,iat2) = frc(1,iat2) - dE
                dE = dy * fact
                frc(2,iat1) = frc(2,iat1) + dE
                frc(2,iat2) = frc(2,iat2) - dE
                dE = dz * fact
                frc(3,iat1) = frc(3,iat1) + dE
                frc(3,iat2) = frc(3,iat2) - dE
              endif
            endif

            next = nxtval(nproc)
          endif

 80     continue
 90   continue
      next = nxtval(-nproc)

      return
      end
