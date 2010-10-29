      subroutine grad1_so
     $     ( H, lbuf, scr, lscr, dens, wdens, frc_nuc,
     $     frc_kin, frc_wgh, 
     $     g_dens, g_wdens, basis, geom, nproc, nat, 
     $     max_at_bf, oskel, 
     &     frc_so, densx, densy, densz)
c$Id$

C     one electron contribution to RHF, ROHF and UHF gradients
C     now also UMP2

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"

C-------------------------parameters--------------------------------
      integer lbuf, lscr,
     $     g_dens(6),        ! density matrix (summed if ROHF, UHF)
     $     g_wdens,       ! weighted density (Lagrangian)
     $     basis, geom, nproc, nat, max_at_bf

      double precision H, ! integral derivatives
     $     scr, 
     $     dens, densx, densy, densz, ! local density block
     $     wdens,               ! local weighted density block
     $     frc_nuc, frc_kin, frc_wgh, frc_so ! forces arrays

      dimension H ( lbuf ), frc_nuc(3, nat), frc_kin(3, nat),
     $     frc_wgh(3, nat), frc_so(3, nat), scr(lscr),
     $     dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf),
     $     densx(max_at_bf,max_at_bf), densy(max_at_bf,max_at_bf),
     $     densz(max_at_bf,max_at_bf) 

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

      logical status, pointforce

      character*16 name
      integer nxtask, task_size
      external nxtask

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

            call ga_get (g_dens(1), 
     &           iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
            call ga_get(g_wdens,
     &           iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)
            call ga_get (g_dens(3), 
     &           iab1f,iab1l,iab2f,iab2l,densz,max_at_bf) 
            call ga_get (g_dens(4), 
     &           iab1f,iab1l,iab2f,iab2l,densy,max_at_bf)
            call ga_get (g_dens(5), 
     &           iab1f,iab1l,iab2f,iab2l,densx,max_at_bf)
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
                    dE = dE * qfac
                    frc_kin(icart,iat3) = frc_kin(icart,iat3) + dE
 40               continue
 50             continue

C     1el. so. derivatives
                call intd_1eso(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)
C     Dso x Hso
                
                ic=1
                do 150, iat3 = 1, nat
                   do 140, icart = 1, 3
c     z componet 
                      dE = 0.D0
                      do 131, ip1 = if1, il1
                         do 128, ip2 = if2, il2
                            dE = dE - densz(ip1,ip2)*H(ic)*0.5d0
                            ic = ic + 1
 128                     continue
 131                  continue
                      if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                      dE = dE * qfac
                      frc_so(icart,iat3) = frc_so(icart,iat3) + dE
c     y componet 
                      dE = 0.D0
                      do 230, ip1 = if1, il1
                         do 231, ip2 = if2, il2
                            dE = dE - densy(ip1,ip2)*H(ic)*0.5d0
                            ic = ic + 1
 231                     continue
 230                  continue
                      if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                      dE = dE * qfac
                      frc_so(icart,iat3) = frc_so(icart,iat3) + dE
c     x component 
                      dE = 0.D0
                      do 250, ip1 = if1, il1
                         do 251, ip2 = if2, il2
                            dE = dE - densx(ip1,ip2)*H(ic)*0.5d0
                            ic = ic + 1
 251                     continue
 250                  continue
                      if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                      dE = dE * qfac
                      frc_so(icart,iat3) = frc_so(icart,iat3) + dE
 140               continue
 150            continue
                
 60          continue
 70       continue
          
 1010     continue

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
c      write(*,'("forces",9f10.5)')((frc_nuc(i,j),i=1,3),j=1,nat)
c      write(*,'("forces",9f10.5)')((frc_kin(i,j),i=1,3),j=1,nat)
c      write(*,'("forces",9f10.5)')((frc_wgh(i,j),i=1,3),j=1,nat)
c      write(*,'("forces",9f10.5)')((frc_so(i,j),i=1,3),j=1,nat)
 
      return
      end
