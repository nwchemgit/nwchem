#define NO_BQGEM 1
      subroutine grad1 ( H, lbuf, scr, lscr, dens, wdens, frc_nuc,
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
#include "bq_params.fh"

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
     $          dens(max_at_bf,max_at_bf), 
     $          wdens(max_at_bf,max_at_bf),
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
      double precision r12
c
      integer nxtask, task_size
      external nxtask

cc AJL/Begin/SPIN ECPs
      integer ecp_channels
      integer iecp
      double precision H_beta (lbuf)
      double precision dens_beta (max_at_bf,max_at_bf)
#ifdef NO_BQGEM
#include "inp.fh"
      character*2 bqchar
#endif

cc    Read this value from rtdb vvvv
      if (.not.rtdb_get(rtdb, 'dft:spin_polarised_ecps',
     &  mt_int, 1, ecp_channels)) then
        ecp_channels = 1
      endif 

cc AJL: With spin-polarised ECPs Hcore will be spin dependent
cc            See Szabo and Ostlund pg. 215
cc      So we need to separate out the densities
cc        
c      if (ecp_channels.gt.1) then
c     
c       Restore alpha and beta densities to calculate spin-polarised
c       derivatives 
c
c        call ga_print(g_dens(1))
c        call ga_print(g_dens(2))
c        call ga_dadd(1d0, g_dens(1), -1d0, g_dens(2), g_dens(1))
c        call ga_print(g_dens(1))
c        call ga_print(g_dens(2))
c      end if   
cc AJL/End

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
c
cc AJL/Begin/SPIN ECPs
c            call ga_get(g_dens,iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
            do 75, iecp = 1, ecp_channels
              if (iecp.eq.1) then
                call ga_get(g_dens(iecp),iab1f,iab1l,
     $               iab2f,iab2l,dens,max_at_bf)
              else
                call ga_get(g_dens(iecp),iab1f,iab1l,
     $               iab2f,iab2l,dens_beta,max_at_bf)
              endif
 75         continue
cc Recombine g_dens, as it is not used again
c            if (ecp_channels.gt.1) then
c              call ga_dadd(1d0, g_dens(1), 1d0, g_dens(2), g_dens(1))
c            end if
cc g_wdens is not dependent on spin, so can leave this
            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)
cc AJL/End
c
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
c
cc AJL/Begin/SPIN ECPs
cc With spin-polarised ECPs Hcore will be spin dependent
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

                ic=1
                do 50, iat3 = 1, nat
                  do 40, icart = 1, 3
                    dE = 0.D0
                    do 31, ip1 = if1, il1
                      do 30, ip2 = if2, il2
                        dE = dE + dens(ip1,ip2) * H(ic)
cc AJL/Begin/SPIN ECPs
                        if (ecp_channels.gt.1) then
                          dE = dE + dens_beta(ip1,ip2) * H_beta(ic)
                        end if
cc AJL/End                 
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
#ifdef NO_BQGEM
         status = geom_cent_get (geom, iat1, name, crd1, q1)
c         bqchar(1:1)=byte_mb(i_t+16*(iat-1))
c         bqchar(2:2)=byte_mb(i_t+16*(iat-1)+1)
         if (.not.inp_compare(.false.,'bq', name(1:2))) then
#endif
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
#ifdef NO_BQGEM
      endif
#endif
 91   continue
      next = nxtask(-nproc,task_size)


      return
      end
c
C> \brief calculate the gradient terms due to the interaction with the 
C> COSMO charges
C>
C> Evaluate the gradient contributions from the COSMO embedding. The
C> original part is from Klamt and Sch&uuml;&uuml;rmann [1] 
C> (see Eqs.(13-16)). The derivatives of matrix \f$A\f$ have been 
C> modified by York and Karplus [2] (see Eqs.(73-76)) to obtain smooth
C> potential energy surfaces. York and Karplus also modified matrix
C> \f$B\f$ which is easy to do in their classical force field code.
C> In an ab-initio code this not so easy to do and as it is not 
C> required to eliminate singularities the original expression from [1]
C> for \f$B\f$ is used here.
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
C> [2] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="https://doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
