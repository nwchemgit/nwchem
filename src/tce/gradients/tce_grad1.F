      subroutine tce_grad1 ( H, lbuf, scr, lscr, dens, wdens, frc_nuc,
     $                   frc_kin, frc_wgh, g_force,
     $                   g_dens, g_wdens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb, oskel, ndens )
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
#include "stdio.fh"

C-------------------------parameters--------------------------------
      integer lbuf, lscr,
     $     g_dens,        ! density matrix (summed if ROHF, UHF)
     $     g_wdens,       ! weighted density (Lagrangian)
     $     g_force,       ! global force array
     $     basis, geom, nproc, nat, max_at_bf, rtdb,
cc AJL/Begin/SPIN ECPs
     $     ndens          ! Number of density matrices in g_dens
cc AJL/End

      double precision H, ! integral derivatives
     $     scr, 
     $     dens,          ! local density block
     $     wdens,         ! local weighted density block
     $     frc_nuc, frc_kin, frc_wgh   ! forces arrays

      dimension H ( lbuf ), frc_nuc(3, nat), frc_kin(3, nat),
     $          frc_wgh(3, nat), scr(lscr),
     $          dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf),
cc AJL/Begin/SPIN ECPs
     $          g_dens(ndens)
cc AJL/End

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

      integer nxtask, task_size
      external nxtask

cc AJL/Begin/SPIN ECPs
      integer ecp_channels
      integer iecp
      double precision H_beta (lbuf)
      double precision dens_beta (max_at_bf,max_at_bf)

cc    Read this value from rtdb vvvv
      if (.not.rtdb_get(rtdb, 'dft:spin_polarised_ecps',
     &  mt_int, 1, ecp_channels)) then
        ecp_channels = 1
        H_beta = 0
      endif
cc AJL: With spin-polarised ECPs Hcore will be spin dependent
cc            See Szabo and Ostlund pg. 215
cc      So we need to separate out the densities
cc        
      if (ecp_channels.gt.1) then
c     
c       Restore alpha and beta densities to calculate spin-polarised
c       derivatives 
c
cc Currently, this will not work, as g_dens is formulated only in the spin-paired
cc manner, and not separable. Therefore, for now I will just put an
cc an error. In the long term we need to be able to supply the separable
cc alpha and beta densities             
      
        write(LuOut,*)'-------------- WARNING ---------------'
        write(LuOut,*)'Spin-polarised ECPs have not been implemented'
        write(LuOut,*)'for TCE gradients. Any results received will be'
        write(LuOut,*)'calculated assuming all density is in the alpha'
        write(LuOut,*)'channel, and therefore must be assumed INCORRECT'
        write(LuOut,*)''
        
        stop

c Note: all of the framework is already in place, for future reference
c we just need to pass in the alpha and beta densities combined in
c g_dens(1), and the beta densities in g_dens(2), and comment out
c the above error.
c
        call ga_dadd(1d0, g_dens(1), -1d0, g_dens(2), g_dens(1))
      end if
cc AJL/End

      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c

      pointforce = geom_include_bqbq(geom)
      dobq = geom_extbq_on()
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
c
cc AJL/Begin/SPIN ECPs
cc            call ga_get(g_dens,iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
            do 75, iecp = 1, ecp_channels
              if (iecp.eq.1) then
                call ga_get(g_dens(iecp), iab1f,iab1l,
     $               iab2f,iab2l,dens,max_at_bf)
              else
                call ga_get(g_dens(iecp), iab1f,iab1l,
     $               iab2f,iab2l,dens_beta,max_at_bf)
              endif
 75         continue
cc g_wdens is not dependent on spin, so can leave this
            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)
cc AJL/End
c
c->debug
c            call ga_print(g_dens)
c             call ga_print(g_wdens)
c<-debug
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
c                write(6,*) 'Dw x S'
                if ( idatom(1) .ge. 1 ) then
C               idatom(1).ge.0 <=> idatom(2).ge.0 (no check necessary)
                  ic = 1
                  do 28, icart = 1, 3
                    de = 0.D0
                    do 22, ip1 = if1, il1
                      do 20, ip2 = if2, il2
                        dE = dE + wdens(ip1,ip2) * H(ic)
                        ic = ic + 1
   20                 continue
   22               continue
                    dE = dE * qfac
                    frc_wgh(icart,idatom(1)) = frc_wgh(icart,idatom(1)) 
     $                                      - dE - dE
                    frc_wgh(icart,idatom(2)) = frc_wgh(icart,idatom(2)) 
     $                                      + dE + dE
c                    write(6,'(a,i4,i4,f20.16)') 'icart,idatom',icart,
c     1                 idatom(1), frc_wgh(icart,idatom(1))
c                    write(6,'(a,i4,i4,f20.16)') 'icart,idatom',icart,
c     1                 idatom(2), frc_wgh(icart,idatom(2))

   28             continue
                endif

C               1el. derivatives
                if(.not.dobq) then
                call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)
                else
                call intd_1epot(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)
                end if
c
cc AJL/Begin/SPIN ECPs
cc AJL: With spin-polarised ECPs Hcore will be spin dependent
cc            See Szabo and Ostlund pg. 215
                if (ecp_channels.gt.1) then
C               1el. derivatives
                  if(.not.dobq) then
cc For now this will do, but this could be more efficiently done
                    call intd_1eh1_beta(basis,ish1,basis,ish2,lscr,scr,
     &                   lbuf,H_beta)
                  else
                    call intd_1epot_beta(basis,ish1,basis,ish2,lscr,scr,
     &                   lbuf,H_beta)
                  end if
                end if
cc AJL/End                  
c
C     D x H
c                write(6,*) 'D x H'
                ic=1
                do 50, iat3 = 1, nat
                  do 40, icart = 1, 3
                    dE = 0.D0
                    do 31, ip1 = if1, il1
                      do 30, ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
cc AJL/Begin
                        if (ecp_channels.gt.1) then
                          dE = dE + dens_beta(ip1,ip2) * H_beta(ic)
                        end if
cc AJL/End
c                        write(6,'(a,i4,i4,f20.16)') 'dens',ip1,ip2,
c     1                        dens(ip1,ip2)
c                        write(6,'(a,i4,f20.16)') 'H',ic,H(ic)
                        ic = ic + 1
   30                 continue
   31               continue
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                    dE = dE * qfac
                    frc_kin(icart,iat3) = frc_kin(icart,iat3) + dE
c                    write(6,'(a,i4,i4,f20.16)') 'icart,iat3',icart,iat3,
c     1                   frc_kin(icart,iat3)
   40             continue
   50           continue

   60         continue
   70       continue

 1010       continue

C     Vnn 
c            write(6,*) 'Vnn'
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
c                write(6,'(a,i4,f20.16)') '1, iat1',iat1,
c     1                frc_nuc(1,iat1)
c                write(6,'(a,i4,f20.16)') '1, iat2',iat2,
c     1                frc_nuc(1,iat2)
                dE = dy * fact
                frc_nuc(2,iat1) = frc_nuc(2,iat1) + dE
                frc_nuc(2,iat2) = frc_nuc(2,iat2) - dE
c                write(6,'(a,i4,f20.16)') '2, iat1',iat1,
c     1                frc_nuc(2,iat1)
c                write(6,'(a,i4,f20.16)') '2, iat2',iat2,
c     1                frc_nuc(2,iat2)
                dE = dz * fact
                frc_nuc(3,iat1) = frc_nuc(3,iat1) + dE
                frc_nuc(3,iat2) = frc_nuc(3,iat2) - dE
c                write(6,'(a,i4,f20.16)') '3, iat1',iat1,
c     1                frc_nuc(3,iat1)
c                write(6,'(a,i4,f20.16)') '3, iat2',iat2,
c     1                frc_nuc(3,iat2)
              endif
            endif

            next = nxtask(nproc,task_size)
          endif

   80   continue
   90 continue
      next = nxtask(-nproc,task_size)

c     
c     continue with forces if there are external charges (MV)
c     --------------------------------------------------
      if(.not.dobq) return

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
            fact = q1 * q2 / SQRT ( dx*dx + dy*dy + dz*dz ) **3
            dE = dx * fact
            frc_nuc(1,iat1) = frc_nuc(1,iat1) + dE
            dE = dy * fact
            frc_nuc(2,iat1) = frc_nuc(2,iat1) + dE
            dE = dz * fact
            frc_nuc(3,iat1) = frc_nuc(3,iat1) + dE

            next = nxtask(nproc,task_size)
          endif
   81   continue
   91 continue
      next = nxtask(-nproc,task_size)


      return
      end
c      subroutine grad_hnd_cos ( H, lbuf, scr, lscr, 
c     $                   dens, wdens, frc_nuc,
c     $                   frc_kin, frc_wgh, g_force,
c     $                   g_dens, g_wdens, basis, geom, nproc, nat, 
c     $                   max_at_bf, rtdb, oskel )
cc$Id$
c
cC     one electron contribution to RHF, ROHF and UHF gradients
cC     now also UMP2
c
c      implicit none
c
c#include "mafdecls.fh"
c#include "global.fh"
c#include "geom.fh"
c#include "bas.fh"
c#include "rtdb.fh"
c#include "sym.fh"
c#include "stdio.fh"
c
cC-------------------------parameters--------------------------------
c      integer lbuf, lscr,
c     $     g_dens,        ! density matrix (summed if ROHF, UHF)
c     $     g_wdens,       ! weighted density (Lagrangian)
c     $     g_force,       ! global force array
c     $     basis, geom, nproc, nat, max_at_bf, rtdb
c
c      double precision H, ! integral derivatives
c     $     scr, 
c     $     dens,          ! local density block
c     $     wdens,         ! local weighted density block
c     $     frc_nuc, frc_kin, frc_wgh   ! forces arrays
c
c      dimension H ( lbuf ), frc_nuc(3, nat), frc_kin(3, nat),
c     $          frc_wgh(3, nat), scr(lscr),
c     $          dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf)
c
c      logical oskel   ! symmetry?
c
cC-------------------------local variables--------------------------
c
c      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
c     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
c     $     if1, il1, if2, il2,
c     $     icart, ic, nint, ip1, ip2
c
c      double precision dE, qfac
c
c      logical status, pointforce
c
c      integer nxtask, task_size
c      external nxtask
c
cc     ---- -cosmo- gradient term -----
c      logical odbug
c
c      odbug=.true.
c      if(odbug) then
c         write(Luout,*) 'in -grad1_hnd_cos- ...'
c      endif
cc
c
c      task_size = 1
c      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c      pointforce = geom_include_bqbq(geom)
c
c      call hf_print_set(1)
c
c      ijatom = -1
c      next = nxtask(nproc,task_size)
c      do 90, iat1 = 1, nat
c        do 80, iat2 = 1, iat1
c
c          ijatom = ijatom + 1
c          if ( ijatom .eq. next ) then
c
c            status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
c            status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
c
c            if (iab1f.le.0 .or. iab2f.le.0) then
cc     
cc     At least one center has no functions on it ... next atom
cc
c              goto 1010
c            endif
c
c            if (oskel) then
c               if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
c     $              goto 1010
c            else
c               qfac = 1.0d0
c            endif
c
c            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
c            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
c
c            call ga_get (g_dens, iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
c            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)
c
c            do 70, ish1 = iac1f, iac1l
c              if ( iat1.eq.iat2 ) iac2l = ish1
c              do 60, ish2 = iac2f, iac2l
c
cC               shell block in atomic (D/Dw)-matrix block
c                status = bas_cn2bfr(basis,ish1,if1,il1)
c                if1 = if1 - iab1f + 1
c                il1 = il1 - iab1f + 1
c                status = bas_cn2bfr(basis,ish2,if2,il2)
c                if2 = if2 - iab2f + 1
c                il2 = il2 - iab2f + 1
c
c                nint = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )
c
c                ic=1
c                do iat3 = 1, nat
c                  do icart = 1, 3
c                    do ip1 = if1, il1
c                      do ip2 = if2, il2
c                        H(ic)=0.D0
c                        ic = ic + 1
c                      enddo
c                    enddo
c                  enddo
c                enddo
c
cC               1el. -cosmo- derivatives
cc-              call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
cc-   &               lbuf,H)
c
cC     D x H
c
c                ic=1
c                do 50, iat3 = 1, nat
c                  do 40, icart = 1, 3
c                    dE = 0.D0
c                    do 31, ip1 = if1, il1
c                      do 30, ip2 = if2, il2
c                        dE = dE + dens(ip1,ip2) * H(ic)
c                        ic = ic + 1
c   30                 continue
c   31               continue
c                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
c                    dE = dE * qfac
c                    frc_kin(icart,iat3) = frc_kin(icart,iat3) + dE
c   40             continue
c   50           continue
c
c   60         continue
c   70       continue
c
c 1010       continue
c
c            next = nxtask(nproc,task_size)
c          endif
c
c   80   continue
c   90 continue
c      next = nxtask(-nproc,task_size)
c
c      return
c      end
