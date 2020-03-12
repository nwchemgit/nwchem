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
      subroutine grad_hnd_cos ( H, lbuf, scr, lscr, 
     $                   dens, frc_cos_nucq, frc_cos_elq,
     $                   frc_cos_qq, 
     $                   g_dens, basis, geom, nproc, nat, 
     $                   max_at_bf, rtdb, oskel )
c$Id: grad1.F 27488 2015-09-10 08:44:11Z alogsdail $

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
#include "cosmoP.fh"
#include "cosmo_params.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bq.fh"
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
      integer l_rad       ! the handle of the atom radii
      integer k_rad       ! the index of the atom radii
      integer nefcl       ! the number of COSMO charge for a given atom
      integer iefc        ! counter over COSMO charges
      integer iefc_c      ! memory index for COSMO charge coordinates
      integer iefc_q      ! memory index for COSMO charges

      double precision dE, qfac, fact, dx, dy, dz, rr
      double precision invscreen, bohr, pi
      double precision zeta1, zeta2, zeta12
      parameter (bohr=0.529177249d0)

      logical status, pointforce

      double precision util_erf, cosff, cosdff
      external         util_erf, cosff, cosdff

      integer nxtask, task_size
      external nxtask
c
      double precision rin, rout, alphai, xyzff
      integer iat
      parameter (alphai = 0.5d0)
      rin(iat)=dbl_mb(k_rad-1+iat)
     &        *(1.0d0-alphai*gammas*sqrt(0.25d0**minbem))
      rout(iat)=dbl_mb(k_rad-1+iat)
     &         *(1.0d0+(1.0d0-alphai)*gammas*sqrt(0.25d0**minbem))

c     ---- -cosmo- gradient term -----
      logical odbug

      pi = acos(-1.0d0)
      odbug=.false.
      if(odbug) then
         write(Luout,*) 'in -grad1_hnd_cos- ...'
      endif
c
      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
      pointforce = geom_include_bqbq(geom)
      if (.not.bq_create('cosmo efc bq',cosmo_bq_efc))
     $   call errquit("grad_hnd_cos: bq_create on cosmo failed",
     $                0,GEOM_ERR)
      if (.not.bq_rtdb_load(rtdb,cosmo_bq_efc))
     $   call errquit('grad_hnd_cos: rtdb load failed for Bq',916,
     $                rtdb_err)
      if (.not.bq_ncenter(cosmo_bq_efc,nefc))
     $   call errquit('grad_hnd_cos: could not retrieve nefc',917,
     $                GEOM_ERR)
      if (.not.bq_index_coord(cosmo_bq_efc,iefc_c))
     $   call errquit('grad_hnd_cos: could not get coordinate index Bq',
     $                cosmo_bq_efc,MA_ERR)
      if (.not.bq_index_charge(cosmo_bq_efc,iefc_q))
     $   call errquit('grad_hnd_cos: could not get charge index Bq',
     $                cosmo_bq_efc,MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nat,"rad",l_rad,k_rad))
     $  call errquit("grad_hnd_cos: could not allocate rad",
     $               ma_sizeof(MT_BYTE,nat,MT_DBL),MA_ERR)
      call cosmo_def_radii(rtdb,geom,nat,dbl_mb(k_rad))
      status = rtdb_get(rtdb,'cosmo:radius',mt_dbl,nat,
     $                  dbl_mb(k_rad))
      do iat1=0,nat-1
        dbl_mb(k_rad+iat1) = dbl_mb(k_rad+iat1)/bohr
      enddo
c
      if (.not.ma_push_get(MT_INT,nefc,"efciat",l_efciat,k_efciat))
     $  call errquit("grad_hnd_cos: could not allocate efciat",
     $               ma_sizeof(MT_BYTE,nefc,MT_INT),MA_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efciat',mt_int,nefc,
     $                 int_mb(k_efciat)))
     $   call errquit('grad_hnd_cos: rtdb get failed for iatefc',915,
     $                rtdb_err)
c
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
     &                 lbuf,H,dbl_mb(iefc_c+3*(iefc-1)),
     &                 dbl_mb(iefc_q+iefc-1),1)

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
c                 Do center C, i.e. the Cosmo charge (associated with
c                 the atom stored in int_mb(k_efciat))
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
c       - 1. The derivative of matrix B (i.e. the Coulomb interaction
c            between the nuclear charge and the surface charge).
c       - 2. The derivative due to the change in the switching
c            function (see [2] Eq.(74)). This only applies for the
c            York-Karplus model.
c
        invscreen = 1.0d0/(1.0d0*screen)
        do ich1 = 1, nefc
          if (do_cosmo_model.eq.DO_COSMO_YK) then
            zeta1 = zeta*sqrt(ptspatm)
     $            / (dbl_mb(k_rad+int_mb(k_efciat+ich1-1)-1)
     $               * sqrt(2.0d0*pi))
            xyzff = 1.0d0
            do iat1 = 1, nat
              if (iat1.ne.int_mb(k_efciat+ich1-1)) then
                dx = coords(1,iat1,geom) - dbl_mb(0+3*(ich1-1)+iefc_c)
                dy = coords(2,iat1,geom) - dbl_mb(1+3*(ich1-1)+iefc_c)
                dz = coords(3,iat1,geom) - dbl_mb(2+3*(ich1-1)+iefc_c)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                if (rr.lt.rout(iat1)) then
                  xyzff = xyzff
     $                  * cosff((rr-rin(iat1))/(rout(iat1)-rin(iat1)))
                endif
              endif
            enddo
          endif
          do iat1 = 1, nat
            dx = coords(1,iat1,geom) - dbl_mb(0+3*(ich1-1)+iefc_c)
            dy = coords(2,iat1,geom) - dbl_mb(1+3*(ich1-1)+iefc_c)
            dz = coords(3,iat1,geom) - dbl_mb(2+3*(ich1-1)+iefc_c)
            rr = sqrt(dx*dx+dy*dy+dz*dz)
c
c           - term 1.
c
            fact = -charge(iat1,geom)*dbl_mb(iefc_q+ich1-1) / 
     $              rr**3
c
c           - term 2.
c
            if (do_cosmo_model.eq.DO_COSMO_YK) then
              if (iat1.ne.int_mb(k_efciat+ich1-1)) then
                fact = fact - 2.0d0*zeta1*sqrt(2.0d0/pi)*invscreen
     $               * dbl_mb(iefc_q+ich1-1)**2
     $               * cosdff((rr-rin(iat1))/(rout(iat1)-rin(iat1)))
     $               / (rr*xyzff*(rout(iat1)-rin(iat1))
     $                  *cosff((rr-rin(iat1))/(rout(iat1)-rin(iat1))))
              endif
            endif
c
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
          if (do_cosmo_model.eq.DO_COSMO_YK) then
            zeta1 = zeta*sqrt(ptspatm)
     $            / (dbl_mb(k_rad+int_mb(k_efciat+ich1-1)-1)
     $               * sqrt(2.0d0*pi))
          endif
          do ich2 = 1, ich1-1
            if (do_cosmo_model.eq.DO_COSMO_YK) then
              zeta2 = zeta*sqrt(ptspatm)
     $              / (dbl_mb(k_rad+int_mb(k_efciat+ich2-1)-1)
     $                 * sqrt(2.0d0*pi))
              zeta12 = zeta1*zeta2/sqrt(zeta1**2+zeta2**2)
            endif
            if (ich1.ne.ich2) then
              dx = dbl_mb(0+3*(ich1-1)+iefc_c)
     $           - dbl_mb(0+3*(ich2-1)+iefc_c)
              dy = dbl_mb(1+3*(ich1-1)+iefc_c)
     $           - dbl_mb(1+3*(ich2-1)+iefc_c)
              dz = dbl_mb(2+3*(ich1-1)+iefc_c)
     $           - dbl_mb(2+3*(ich2-1)+iefc_c)
              rr = sqrt(dx*dx+dy*dy+dz*dz)
              if (rr.lt.1.0d-6) then
                fact = 0.0d0
              else
                if (do_cosmo_model.eq.DO_COSMO_KS) then
                  fact = -invscreen*dbl_mb(iefc_q+ich1-1)
     $                 * dbl_mb(iefc_q+ich2-1)/(rr**3)
                else if (do_cosmo_model.eq.DO_COSMO_YK) then
                  fact = +invscreen*dbl_mb(iefc_q+ich1-1)
     $                 * dbl_mb(iefc_q+ich2-1) 
     $                 * (2.0d0*zeta12/sqrt(pi)/(rr**2)
     $                    * exp(-(zeta12*rr)**2)
     $                    - util_erf(zeta12*rr)/(rr**3))
                else
                  call errquit("grad_hnd_cos: panic",
     $                         do_cosmo_model,UERR)
                endif
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
      if (.not.ma_pop_stack(l_rad))
     $   call errquit("grad_hnd_cos: could not deallocate l_rad",
     $                0,MA_ERR)
      if (.not.bq_destroy(cosmo_bq_efc))
     $   call errquit("grad_hnd_cos: bq_destroy on cosmo failed",
     $                0,GEOM_ERR)

      return
      end
