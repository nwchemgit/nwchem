      subroutine grad1 ( H, lbuf, scr, lscr, dens, wdens, frc_nuc,
     $                   frc_kin, frc_wgh, g_force,
     $                   g_dens, g_wdens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb, oskel )
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
     &                 lbuf,H)
                else
                  call intd_1epot(basis,ish1,basis,ish2,lscr,scr,
     &                 lbuf,H)
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


      return
      end
c
C> \brief calculate the gradient terms due to the interaction with the 
C> COSMO charges
      subroutine grad_hnd_cos ( H, lbuf, scr, lscr, 
     $                   dens, frc_cos_nucq, frc_cos_elq,
     $                   frc_cos_qq, 
     $                   g_dens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb, oskel )
c$Id$

C     COSMO one electron contribution to RHF, ROHF and UHF gradients
C     now also UMP2 ??? unlikely as that requires solutions to the
c     CPHF equation???
c
c     Terms included in this subroutine are:
c     1. Electron - COSMO charge interactions
c     2. Nuclear - COSMO charge interactions
c     3. COSMO charge - COSMO charge interactions
c
c     Terms NOT included are:
c     1. All regular QM derivatives

      implicit none

#include "nwc_const.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "stdio.fh"
#include "prop.fh"

C-------------------------parameters--------------------------------
      integer g_dens !< [Input] the total electron density matrix GA
      integer basis  !< [Input] the basis set handle
      integer geom   !< [Input] the geometry handle
      integer rtdb   !< [Input] the RTDB handle
      integer lbuf   !< [Input] the length of the integral buffer
      integer lscr,  !< [Input] the length of the scratch space
     $     nproc, nat, max_at_bf

      double precision frc_cos_nucq(3,nat) !< [Output] the forces due
                                           !< nuclear-COSMO charge 
                                           !< interaction
      double precision frc_cos_elq(3,nat)  !< [Output] the forces due
                                           !< electron-COSMO charge 
                                           !< interaction
      double precision frc_cos_qq(3,nat)   !< [Output] the forces due
                                           !< COSMO charge-COSMO charge 
                                           !< interaction
      double precision H(lbuf)   !< [Scratch] the derivative integrals
      double precision scr(lscr) !< [Scratch] scratch space
      double precision dens(max_at_bf,max_at_bf) !< [Scratch] local
                                                 !< density block

      logical oskel   ! symmetry?
c
      double precision dielec, screen, rsolv
      common/hnd_cospar/dielec, screen ,rsolv

C-------------------------local variables--------------------------

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2
      integer im1, im2, nprim, ngen, sphcart, ityp1, ityp2
      integer ich1, ich2

      integer nefc        ! the number of COSMO charges
      integer l_efciat    ! the handle of the COSMO charge-atom map
      integer k_efciat    ! the index of the COSMO charge-atom map
      integer nefcl       ! the number of COSMO charge for a given atom
      integer iefc        ! counter over COSMO charges

      double precision dE, qfac, fact, dx, dy, dz, rr
      double precision invscreen

      logical status, pointforce

      integer nxtask, task_size
      external nxtask

c     ---- -cosmo- gradient term -----
      logical odbug

      odbug=.false.
      if(odbug) then
         write(Luout,*) 'in -grad1_hnd_cos- ...'
      endif
c

      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
      pointforce = geom_include_bqbq(geom)
      if (.not.geom_create(cosmo_geom_efc,'cosmo efc geom'))
     $   call errquit("grad_hnd_cos: geom_create on cosmo failed",
     $                0,GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,cosmo_geom_efc,
     $                        'cosmo efc geom'))
     $   call errquit("grad_hnd_cos: geom_rtdb_load on cosmo failed",
     $                0,RTDB_ERR)
      if (.not.geom_ncent(cosmo_geom_efc,nefc))
     $  call errquit("grad_hnd_cos: could not establish number of "//
     $               "COSMO charges",cosmo_geom_efc,GEOM_ERR)
      if (.not.ma_push_get(MT_INT,nefc,"efciat",l_efciat,k_efciat))
     $  call errquit("grad_hnd_cos: could not allocate efciat",
     $               ma_sizeof(MT_BYTE,nefc,MT_INT),MA_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efciat',mt_int,nefc,
     $                 int_mb(k_efciat)))
     $   call errquit('grad_hnd_cos: rtdb get failed for iatefc',915,
     $                rtdb_err)

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
            qfac = -qfac

            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)

            call ga_get (g_dens, iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)

            do 70, ish1 = iac1f, iac1l
              if ( iat1.eq.iat2 ) iac2l = ish1
              do 60, ish2 = iac2f, iac2l

C               shell block in atomic (D/Dw)-matrix block
                status = bas_cn2bfr(basis,ish1,if1,il1)
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
c
c               Work out the number of Cartesian basis functions
c               The integrals are evaluated in the Cartesian basis set
c               and then transformed to spherical harmonics. So the
c               buffer size depends on the number of Cartesian functions
c
                status = bas_continfo(basis,ish1,ityp1,nprim,ngen,
     +                                sphcart)
                if (sphcart.eq.1.and.ityp1.ge.2) then
                  im1 = if1 + (ityp1+1)*(ityp1+2)/2 - 1
                else
                  im1 = il1
                endif
                status = bas_cn2bfr(basis,ish2,if2,il2)
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1
c
c               Same Cartesian vs spherical harmonic catastrophy as
c               for ish1.
c
                status = bas_continfo(basis,ish2,ityp2,nprim,ngen,
     +                                sphcart)
                if (sphcart.eq.1.and.ityp2.ge.2) then
                  im2 = if2 + (ityp2+1)*(ityp2+2)/2 - 1
                else
                  im2 = il2
                endif

                nint = ( im1 - if1 + 1 ) * ( im2 - if2 + 1 )

                do iefc = 1, nefc

                  ic=1
                  do iat3 = 1, 3 ! centers A, B, and C
                    do icart = 1, 3
                      do ip1 = if1, im1
                        do ip2 = if2, im2
                          H(ic)=0.0D0
                          ic = ic + 1
                        enddo
                      enddo
                    enddo
                  enddo

C                 1el. -cosmo- derivatives
c                 Currently calculated on for every COSMO charge
c                 separately.
                  call intd_1epot_cosmo(basis,ish1,basis,ish2,lscr,scr,
     &                 lbuf,H,coords(1,iefc,cosmo_geom_efc),
     &                 charge(iefc,cosmo_geom_efc),1)

C     D x H
c
c                 Do center A (associated with ish1)
c
                  ic=1
                  do icart = 1, 3
                    dE = 0.D0
                    do ip1 = if1, il1
                      do ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
                        ic = ic + 1
                      enddo
                    enddo
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                    dE = dE * qfac
                    frc_cos_elq(icart,iat1)
     &              = frc_cos_elq(icart,iat1) - dE
                  enddo
c
c                 Do center B (associated with ish2)
c
                  do icart = 1, 3
                    dE = 0.D0
                    do ip1 = if1, il1
                      do ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
                        ic = ic + 1
                      enddo
                    enddo
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) then
                      dE = dE + dE
                    else
                      dE = 0.0d0
                    endif
                    dE = dE * qfac
                    frc_cos_elq(icart,iat2)
     &              = frc_cos_elq(icart,iat2) - dE
                  enddo
c
c                 Do center C, i.e. the Cosmo charge
c
                  do icart = 1, 3
                    dE = 0.D0
                    do ip1 = if1, il1
                      do ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
                        ic = ic + 1
                      enddo
                    enddo
                    if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
                    dE = dE * qfac
                    frc_cos_elq(icart,int_mb(k_efciat+iefc-1))
     &              = frc_cos_elq(icart,int_mb(k_efciat+iefc-1)) - dE
                  enddo

                enddo

 60           continue
 70         continue

 1010       continue

            next = nxtask(nproc,task_size)
          endif

 80     continue
 90   continue
      next = nxtask(-nproc,task_size)
c
      if (ga_nodeid().eq.0) then
c
c       Do the Nuclear - Cosmo charge part
c
        do iat1 = 1, nat
          do ich1 = 1, nefc
            dx = coords(1,iat1,geom) - coords(1,ich1,cosmo_geom_efc)
            dy = coords(2,iat1,geom) - coords(2,ich1,cosmo_geom_efc)
            dz = coords(3,iat1,geom) - coords(3,ich1,cosmo_geom_efc)
            fact = -charge(iat1,geom)*charge(ich1,cosmo_geom_efc) / 
     $             sqrt(dx*dx+dy*dy+dz*dz)**3
            dE = dx * fact
            frc_cos_nucq(1,iat1) = frc_cos_nucq(1,iat1) + dE
            frc_cos_nucq(1,int_mb(k_efciat+ich1-1)) 
     $      = frc_cos_nucq(1,int_mb(k_efciat+ich1-1)) - dE
            dE = dy * fact
            frc_cos_nucq(2,iat1) = frc_cos_nucq(2,iat1) + dE
            frc_cos_nucq(2,int_mb(k_efciat+ich1-1)) 
     $      = frc_cos_nucq(2,int_mb(k_efciat+ich1-1)) - dE
            dE = dz * fact
            frc_cos_nucq(3,iat1) = frc_cos_nucq(3,iat1) + dE
            frc_cos_nucq(3,int_mb(k_efciat+ich1-1)) 
     $      = frc_cos_nucq(3,int_mb(k_efciat+ich1-1)) - dE
          enddo
        enddo
c
c       Do cosmo charge - cosmo charge interaction
c
        invscreen = 1.0d0/(1.0d0*screen)
        do ich1 = 1, nefc
          do ich2 = 1, ich1-1
            if (ich1.ne.ich2) then
              dx = coords(1,ich1,cosmo_geom_efc)
     $           - coords(1,ich2,cosmo_geom_efc)
              dy = coords(2,ich1,cosmo_geom_efc)
     $           - coords(2,ich2,cosmo_geom_efc)
              dz = coords(3,ich1,cosmo_geom_efc)
     $           - coords(3,ich2,cosmo_geom_efc)
              rr = sqrt(dx*dx+dy*dy+dz*dz)
              if (rr.lt.1.0d-10) then
                fact = 0.0d0
              else
                fact = -invscreen*charge(ich1,cosmo_geom_efc)
     $               * charge(ich2,cosmo_geom_efc) / 
     $                 rr**3
              endif
              dE = dx * fact
              frc_cos_qq(1,int_mb(k_efciat+ich1-1)) 
     $        = frc_cos_qq(1,int_mb(k_efciat+ich1-1)) + dE
              frc_cos_qq(1,int_mb(k_efciat+ich2-1)) 
     $        = frc_cos_qq(1,int_mb(k_efciat+ich2-1)) - dE
              dE = dy * fact
              frc_cos_qq(2,int_mb(k_efciat+ich1-1)) 
     $        = frc_cos_qq(2,int_mb(k_efciat+ich1-1)) + dE
              frc_cos_qq(2,int_mb(k_efciat+ich2-1)) 
     $        = frc_cos_qq(2,int_mb(k_efciat+ich2-1)) - dE
              dE = dz * fact
              frc_cos_qq(3,int_mb(k_efciat+ich1-1)) 
     $        = frc_cos_qq(3,int_mb(k_efciat+ich1-1)) + dE
              frc_cos_qq(3,int_mb(k_efciat+ich2-1)) 
     $        = frc_cos_qq(3,int_mb(k_efciat+ich2-1)) - dE
            endif
          enddo
        enddo
      endif
c
      if (.not.ma_pop_stack(l_efciat))
     $   call errquit("grad_hnd_cos: could not deallocate l_efciat",
     $                0,MA_ERR)
      if (.not.geom_destroy(cosmo_geom_efc))
     $   call errquit("grad_hnd_cos: geom_destroy on cosmo failed",
     $                0,GEOM_ERR)

      return
      end
