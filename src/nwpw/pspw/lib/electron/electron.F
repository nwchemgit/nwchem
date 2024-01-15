*
* $Id$
*


*     ***********************************
*     *					*
*     *		electron_init		*
*     *					*
*     ***********************************
      subroutine electron_init()
      implicit none

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"


*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter


*     **** local variables ****
      logical value
      integer n2ft3d

*     **** external functions ****
      logical  ion_chargeexist,ion_mmexist,psp_pawexist
      external ion_chargeexist,ion_mmexist,psp_pawexist
      logical  control_Efield_on
      external control_Efield_on
      integer  psi_ispin,psi_ne,psi_neq,control_version
      external psi_ispin,psi_ne,psi_neq,control_version

!$OMP MASTER
      counter = 0
!$OMP END MASTER

      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)
      neq(1) = psi_neq(1)
      neq(2) = psi_neq(2)
      field_exist = ion_chargeexist().or.
     >              ion_mmexist()    .or.
     >              control_Efield_on()
      confine_exist = .false.

*     **** get nfft3d, and n2ft3d ****
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

c      paw_exist = psp_pawexist()

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'Hpsi_k',Hpsi_k(2),Hpsi_k(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d*(neq(1)+neq(2)),
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'vl2',vl(2),vl(1))

      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,
     >                     'vl_lr',vl_lr(2),vl_lr(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,
     >                     'v_field',v_field(2),v_field(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'vall',vall(2),vall(1))


      if (control_version().eq.3) then
        value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'vc',vc(2),vc(1))
      end if

      if (control_version().eq.4) then
        value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,
     >                     'vc',vc(2),vc(1))
      end if

      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'xcp',xcp(2),xcp(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'xce',xce(2),xce(1))
      if (.not. value) 
     >  call errquit('electron_init: out of heap memory',0, MA_ERR)
c      call dcopy(n2ft3d*(neq(1)+neq(2)),0.0d0,0,dbl_mb(psi_r(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(vall(1)),1)
c      call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)
c      call dcopy(n2ft3d,0.0d0,0,dbl_mb(vl_lr(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(xcp(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(xce(1)),1)
      call Parallel_shared_vector_zero(.false.,(neq(1)+neq(2))*n2ft3d,
     >                                 dbl_mb(psi_r(1)))
      call Parallel_shared_vector_zero(.false.,2*n2ft3d,dbl_mb(vall(1)))
      call Parallel_shared_vector_zero(.false.,n2ft3d,
     >                                 dbl_mb(v_field(1)))
      call Parallel_shared_vector_zero(.false.,n2ft3d,dbl_mb(vl_lr(1)))
      call Parallel_shared_vector_zero(.false.,2*n2ft3d,dbl_mb(xcp(1)))
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,dbl_mb(xce(1)))

*     *** initializer symmetrizer ***
      call rho_symmetrizer_init()

      return
      end


*     ***********************************
*     *					*
*     *		electron_finalize       *
*     *					*
*     ***********************************
      subroutine electron_finalize()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      logical value

*     **** free heap  memory ****
      call rho_symmetrizer_end()

      value = BA_free_heap(Hpsi_k(2))
      value = value.and.
     >        BA_free_heap(psi_r(2))
      value = value.and.
     >        BA_free_heap(vl(2))
      value = value.and.
     >        BA_free_heap(vl_lr(2))
      value = value.and.
     >        BA_free_heap(v_field(2))
      value = value.and.
     >        BA_free_heap(vall(2))
      value = value.and.
     >        BA_free_heap(vc(2))
      value = value.and.
     >        BA_free_heap(xcp(2))
      value = value.and.
     >        BA_free_heap(xce(2))
      if (.not. value) 
     >  call errquit('electron_init: error freeing heap memory',0,
     &       MA_ERR)

      return
      end

*     ***********************************
*     *					*
*     *		electron_count		*
*     *					*
*     ***********************************
      integer function electron_count()
      implicit none

*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter

      electron_count = counter
      return
      end


*     ***********************************
*     *                                 *
*     *         electron_run0           *
*     *                                 *
*     ***********************************
      subroutine electron_run0(psi_k)
      implicit none
      complex*16 psi_k(*)
     
*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter
     
!$OMP MASTER                  
      counter = counter+1     
!$OMP END MASTER
      
      call electron_gen_psi_r(psi_k)
      call electron_gen_Hvallpsi_k(psi_k)

      return
      end



*     ***********************************
*     *					*
*     *		electron_run		*
*     *					*
*     ***********************************
      subroutine electron_run(psi_k,dn,dng,dnall,fractional,occ)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)
      logical    fractional
      real*8     occ(*)

*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter

!$OMP MASTER
      counter = counter+1
!$OMP END MASTER

      call electron_gen_psi_r(psi_k)
      call electron_gen_densities(psi_k,dn,dng,dnall,fractional,occ)
      call electron_gen_scf_potentials(dn,dng,dnall)
      call electron_gen_Hpsi_k(psi_k)

      return
      end

*     ***********************************
*     *					*
*     *		electron_genrho		*
*     *					*
*     ***********************************
      subroutine electron_genrho(psi_k,dn,fractional,occ)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      logical    fractional
      real*8     occ(*)

      call electron_gen_psi_r(psi_k)
      call electron_gen_density(psi_k,dn,fractional,occ)
      return
      end


*     ***********************************
*     *					*
*     *		electron_run_orb	*
*     *					*
*     ***********************************
      subroutine electron_run_orb(i,psi_k)
      implicit none
      integer    i
      complex*16 psi_k(*)

      call electron_gen_psi_r_orb(i,psi_k)
      call electron_gen_Hpsi_k_orb(i,psi_k)

      return
      end

*     ***********************************
*     *					*
*     *		electron_sd_update  	*
*     *					*
*     ***********************************
      subroutine electron_sd_update(psi1,psi2,dte)
      implicit none
      complex*16 psi1(*),psi2(*)
      real*8     dte

#include "bafdecls.fh"
#include "electron_common.fh"



c      call ke_Precondition(npack1,(ne(1)+ne(2)),psi1,dcpl_mb(Hpsi_k(1)))
c     call electron_sd_subupdate(npack1,(ne(1)+ne(2)),
c    >                           psi1,psi2,dcpl_mb(Hpsi_k(1)),
c    >                           dte)
      call dcopy(2*(npack1)*(neq(1)+neq(2)),psi1,1,psi2,1)
      call daxpy(2*(npack1)*(neq(1)+neq(2)),
     >           (-dte),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           psi2,1)

      return
      end
*     ***********************************
*     *                                 *
*     *         electron_cpmd_update    *
*     *                                 *
*     ***********************************
      subroutine electron_cpmd_update(psi0,psi1,psi2,hml,dte)
      implicit none
      complex*16 psi0(*),psi1(*),psi2(*)
      real*8     hml(*)
      real*8     dte

#include "bafdecls.fh"
#include "electron_common.fh"

      !**** psi2 = 2*psi1 - psi0 + dte*Hpsi ****
      !**** - note that Hpsi = minus the gradient in electron ****

      !**** rotate Hpsi ****
c      call electron_gen_Hpsi_k(psi1)
      call Dneall_fmf_Multiply(0,dcpl_mb(Hpsi_k(1)),npack1,
     >                         hml,(-dte),
     >                         psi2,0.0d0)
      call daxpy(2*(npack1)*(neq(1)+neq(2)),
     >           1.0d0,
     >           psi1,1,
     >           psi2,1)

c      call daxpy(2*(npack1)*(neq(1)+neq(2)),
c     >           -1.0d0,
c     >           psi0,1,
c     >           psi2,1)
c      call daxpy(2*(npack1)*(neq(1)+neq(2)),
c     >           (-dte),
c     >           dcpl_mb(Hpsi_k(1)),1,
c     >           psi2,1)

      return
      end



c      subroutine electron_sd_subupdate(npack1,nn,
c     >                                 psi1,psi2,Hpsi,dte)
c      implicit none
c      integer    npack1,nn
c      complex*16 psi1(npack1,nn)
c      complex*16 psi2(npack1,nn)
c      complex*16 Hpsi(npack1,nn)
c      real*8     dte
c
c      integer n
c
c*     ************************************
c*     **** do a steepest descent step ****
c*     ************************************
c      do n=1,nn
c        call Pack_c_SMul(1,(-dte),Hpsi(1,n),psi2(1,n))
c        call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
c      end do
c
c      return
c      end


 
*     ***********************************
*     *					*
*     *		electron_energy		*
*     *					*
*     ***********************************
      real*8 function electron_energy(psi_k,dn,dng,dnall,fractional,occ)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)
      logical fractional
      real*8  occ(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer n2ft3d
      integer ii,ms,n1(2),n2(2),nx,ny,nz
      real*8  sum,eorbit,ehartr,exc,pxc,exc2,pxc2,dv,eion_core
      real*8  ehartree_atom,ecmp_cmp,ecmp_pw,exc_atom,pxc_atom,eke_core
      real*8  total_energy

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2


*     **** external functions *****
      logical  pspw_SIC,pspw_SIC_relaxed,psp_U_psputerm,psp_pawexist
      logical  pspw_HFX,pspw_HFX_relaxed,meta_found,nwpw_meta_gga_on
      integer  control_version
      real*8   lattice_omega,coulomb_e,electron_ehartree2
      real*8   nwpw_meta_gga_pxc,psp_kinetic_core,psp_ion_core
      real*8   psp_hartree_atom,psp_hartree_cmp_cmp
      real*8   psp_hartree_cmp_pw
      external pspw_SIC,pspw_SIC_relaxed,psp_U_psputerm,psp_pawexist
      external pspw_HFX,pspw_HFX_relaxed,meta_found,nwpw_meta_gga_on
      external control_version
      external lattice_omega,coulomb_e,electron_ehartree2
      external nwpw_meta_gga_pxc,psp_kinetic_core,psp_ion_core
      external psp_hartree_atom,psp_hartree_cmp_cmp
      external psp_hartree_cmp_pw
      logical  pspw_V_APC_on
      external pspw_V_APC_on
      integer  pspw_Efield_type
      external pspw_Efield_type
      real*8   dipole_Efield_e,dipole_Efield_p
      external dipole_Efield_e,dipole_Efield_p
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
     

      dv = lattice_omega()/dble(nx*ny*nz)

      n2ft3d = 2*nfft3d
      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)


*     *** get orbital energies ****
!$OMP MASTER
      e1 = 0.0d0
      e2 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      eorbit = 0.0d0
      if (fractional) then
         do ms=1,ispin
         do ii=n1(ms),n2(ms)
           call Pack_cc_idot(1,
     >                       psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       e1)
!$OMP MASTER
           !eorbit = eorbit + sum*occ(ii)
           e2 = e2 + e1*occ(ii)
!$OMP END MASTER
         end do
         end do
      else
         do ms=1,ispin
         do ii=n1(ms),n2(ms)
           call Pack_cc_idot(1,
     >                       psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       e1)
!$OMP MASTER
           !eorbit = eorbit + sum
           e2 = e2 + e1
!$OMP END MASTER
         end do
         end do
      end if
!$OMP BARRIER
      call Parallel_SumAll(e2)
      eorbit = e2
      if (ispin.eq.1) eorbit = eorbit+eorbit


           
*     **** get coulomb energy ****
      ehartr = 0.0d0
      if (control_version().eq.3) ehartr = coulomb_e(dng)
      if (control_version().eq.4) ehartr = electron_ehartree2(dn)
        


*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dnall(1),
     >                dbl_mb(xce(1)),
     >                e1)
      call D3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     .                 e2)
      exc = e1
      pxc = e2
      if (ispin.eq.1) then
         exc= exc + exc 
         pxc= pxc + pxc 
      else
         call D3dB_rr_dot(1,dnall(1+n2ft3d),
     >                    dbl_mb(xce(1)),
     >                    e1)
         call D3dB_rr_dot(1,dn(1+n2ft3d),
     >                    dbl_mb(xcp(1)+n2ft3d),
     >                    e2)
         exc2 = e1
         pxc2 = e2
         exc = exc + exc2
         pxc = pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv

*     **** meta_GGA energy ****
      if (nwpw_meta_gga_on()) then
         pxc = pxc + nwpw_meta_gga_pxc(ispin,neq,psi_k)
      end if

      total_energy = eorbit + exc - ehartr - pxc



*     **** PAW ee terms ****
      if (psp_pawexist()) then
         eke_core      = psp_kinetic_core()
         eion_core     = psp_ion_core()
         ehartree_atom = psp_hartree_atom(ispin,neq,psi_k)
         ecmp_cmp      = psp_hartree_cmp_cmp(ispin)
         ecmp_pw       = psp_hartree_cmp_pw(ispin,dng,dn)
         call psp_xc_atom(ispin,neq,psi_k,exc_atom,pxc_atom)

         total_energy = total_energy + exc_atom - pxc_atom
     >                - ehartree_atom - ecmp_cmp - ecmp_pw
c     >                + eke_core + eion_core
      end if

*     **** SIC corrections ****
      if (pspw_SIC()) then
      if (pspw_SIC_relaxed()) then
         call pspw_energy_SIC(ispin,dbl_mb(psi_r(1)),
     >                        ehsic,
     >                        phsic,
     >                        exsic,
     >                        pxsic)
         total_energy = total_energy + ehsic + exsic - phsic - pxsic
      end if
      end if

*     **** HFX energy ****
      if (pspw_HFX()) then
      if (pspw_HFX_relaxed()) then
         call pspw_energy_HFX(ispin,dbl_mb(psi_r(1)),
     >                        ehfx,
     >                        phfx)
         total_energy = total_energy + ehfx - phfx
      end if
      end if

*     **** DFT+U energy ****
      if (psp_U_psputerm()) then
         call psp_U_psputerm_energy(edftu,pdftu)
         total_energy = total_energy + edftu - pdftu
      end if


*     **** metadynamics energy ****
      if (meta_found()) then
         call meta_energypotential(ispin,neq,psi_k,emeta,pmeta)
         total_energy = total_energy + emeta - pmeta
      end if

*     **** APC energy ****
      if (pspw_V_APC_on()) then
         call pspw_E_APC(eapc,papc)
         total_energy = total_energy + eapc - papc
      end if

*     **** get pspw_charge and pspw_Efield energies ****
      if ((field_exist).and.(pspw_Efield_type().eq.0)) then
         total_energy=total_energy+dipole_Efield_e()-dipole_Efield_p()
      end if


*     **** total energy ****
      electron_energy = total_energy

      return
      end

*     ***********************************
*     *					*
*     *		electron_energy0	*
*     *					*
*     ***********************************
      real*8 function electron_energy0(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer ii,ms,n1(2),n2(2)
      real*8  eorbit,total_energy

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2

      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)

*     *** get orbital energies ****
!$OMP MASTER
      e1 = 0.0d0
      e2 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      !eorbit = 0.0d0
      do ms=1,ispin
      do ii=n1(ms),n2(ms)
        call Pack_cc_idot(1,psi_k(1+(ii-1)*npack1),
     >                    dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                    e1)
!$OMP MASTER
        e2 = e2 + e1
!$OMP END MASTER
      end do
      end do
!$OMP BARRIER
      call Parallel_SumAll(e2)
      eorbit = e2
      if (ispin.eq.1) eorbit = eorbit+eorbit

      total_energy = eorbit !*** add other energies here ? ***

*     **** total energy ****
      electron_energy0 = total_energy

      return
      end



*     ***********************************
*     *                                 *
*     *    electron_energy_remainder    *
*     *                                 *
*     ***********************************
      real*8 function electron_energy_remainder(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer n2ft3d
      integer ii,ms,n1(2),n2(2),nx,ny,nz
      real*8  sum,eorbit,ehartr,exc,pxc,exc2,pxc2,dv,eion_core
      real*8  ehartree_atom,ecmp_cmp,ecmp_pw,exc_atom,pxc_atom,eke_core
      real*8  total_energy

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2


*     **** external functions *****
      logical  pspw_SIC,pspw_SIC_relaxed,psp_U_psputerm,psp_pawexist
      logical  pspw_HFX,pspw_HFX_relaxed,meta_found,nwpw_meta_gga_on
      integer  control_version
      real*8   lattice_omega,coulomb_e,electron_ehartree2
      real*8   nwpw_meta_gga_pxc,psp_kinetic_core,psp_ion_core
      real*8   psp_hartree_atom,psp_hartree_cmp_cmp
      real*8   psp_hartree_cmp_pw
      external pspw_SIC,pspw_SIC_relaxed,psp_U_psputerm,psp_pawexist
      external pspw_HFX,pspw_HFX_relaxed,meta_found,nwpw_meta_gga_on
      external control_version
      external lattice_omega,coulomb_e,electron_ehartree2
      external nwpw_meta_gga_pxc,psp_kinetic_core,psp_ion_core
      external psp_hartree_atom,psp_hartree_cmp_cmp
      external psp_hartree_cmp_pw
      logical  pspw_V_APC_on
      external pspw_V_APC_on
      integer  pspw_Efield_type
      external pspw_Efield_type
      real*8   dipole_Efield_e,dipole_Efield_p
      external dipole_Efield_e,dipole_Efield_p

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)


      dv = lattice_omega()/dble(nx*ny*nz)

      n2ft3d = 2*nfft3d
      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)


*     *** get orbital energies - zero out eorbit ****
      eorbit = 0.0d0 

*     **** get coulomb energy ****
      ehartr = 0.0d0
      if (control_version().eq.3) ehartr = coulomb_e(dng)
      if (control_version().eq.4) ehartr = electron_ehartree2(dn)

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dnall(1),
     >                dbl_mb(xce(1)),
     >                e1)
      call D3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     .                 e2)
      exc = e1
      pxc = e2
      if (ispin.eq.1) then
         exc= exc + exc
         pxc= pxc + pxc
      else
         call D3dB_rr_dot(1,dnall(1+n2ft3d),
     >                    dbl_mb(xce(1)),
     >                    e1)
         call D3dB_rr_dot(1,dn(1+n2ft3d),
     >                    dbl_mb(xcp(1)+n2ft3d),
     >                    e2)
         exc2 = e1
         pxc2 = e2
         exc = exc + exc2
         pxc = pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv

*     **** meta_GGA energy ****
      if (nwpw_meta_gga_on()) then
         pxc = pxc + nwpw_meta_gga_pxc(ispin,neq,psi_k)
      end if

      total_energy = eorbit + exc - ehartr - pxc



*     **** PAW ee terms ****
      if (psp_pawexist()) then
         eke_core      = psp_kinetic_core()
         eion_core     = psp_ion_core()
         ehartree_atom = psp_hartree_atom(ispin,neq,psi_k)
         ecmp_cmp      = psp_hartree_cmp_cmp(ispin)
         ecmp_pw       = psp_hartree_cmp_pw(ispin,dng,dn)
         call psp_xc_atom(ispin,neq,psi_k,exc_atom,pxc_atom)

         total_energy = total_energy + exc_atom - pxc_atom
     >                - ehartree_atom - ecmp_cmp - ecmp_pw
c     >                + eke_core + eion_core
      end if

*     **** SIC corrections ****
      if (pspw_SIC()) then
      if (pspw_SIC_relaxed()) then
         call pspw_energy_SIC(ispin,dbl_mb(psi_r(1)),
     >                        ehsic,
     >                        phsic,
     >                        exsic,
     >                        pxsic)
         total_energy = total_energy + ehsic + exsic - phsic - pxsic
      end if
      end if

*     **** HFX energy ****
      if (pspw_HFX()) then
      if (pspw_HFX_relaxed()) then
         call pspw_energy_HFX(ispin,dbl_mb(psi_r(1)),
     >                        ehfx,
     >                        phfx)
         total_energy = total_energy + ehfx - phfx
      end if
      end if

*     **** DFT+U energy ****
      if (psp_U_psputerm()) then
         call psp_U_psputerm_energy(edftu,pdftu)
         total_energy = total_energy + edftu - pdftu
      end if


*     **** metadynamics energy ****
      if (meta_found()) then
         call meta_energypotential(ispin,neq,psi_k,emeta,pmeta)
         total_energy = total_energy + emeta - pmeta
      end if

*     **** APC energy ****
      if (pspw_V_APC_on()) then
         call pspw_E_APC(eapc,papc)
         total_energy = total_energy + eapc - papc
      end if

*     **** get pspw_charge and pspw_Efield energies ****
      if ((field_exist).and.(pspw_Efield_type().eq.0)) then
         total_energy=total_energy+dipole_Efield_e()-dipole_Efield_p()
      end if


*     **** total energy ****
      electron_energy_remainder = total_energy

      return
      end




*     ***********************************
*     *                                 *
*     *         electron_eorbit_noocc   *
*     *                                 *
*     ***********************************
      real*8 function electron_eorbit_noocc(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"

*     **** local variables ****
      integer ii,ms,n1(2),n2(2)
      real*8  sum,eorbit

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2


      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)

*     *** get orbital energies ****
!$OMP MASTER
      e1 = 0.0d0
      e2 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      !eorbit = 0.0d0
      do ms=1,ispin
      do ii=n1(ms),n2(ms)
        call Pack_cc_idot(1,psi_k(1+(ii-1)*npack1),
     >                    dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                    e1)
!$OMP MASTER
        e2 = e2 + e1
!$OMP END MASTER
      end do
      end do
!$OMP BARRIER
      call Parallel_SumAll(e2)
      eorbit = e2
      if (ispin.eq.1) eorbit = eorbit+eorbit

      electron_eorbit_noocc = eorbit
      return
      end



*     ***********************************
*     *					*
*     *		electron_eorbit		*
*     *					*
*     ***********************************
      real*8 function electron_eorbit(psi_k,fractional,occ)
      implicit none
      complex*16 psi_k(*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer ii,ms,n1(2),n2(2)
      real*8  sum,eorbit

      common /eelectron_ejtmp/ sum,eorbit

      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)


*     *** get orbital energies ****
      eorbit = 0.0d0
      if (fractional) then
         do ms=1,ispin
         do ii=n1(ms),n2(ms)
           call Pack_cc_idot(1,psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       sum)
!$OMP MASTER
           eorbit = eorbit + sum*occ(ii)
!$OMP END MASTER
         end do
         end do
      else
         do ms=1,ispin
         do ii=n1(ms),n2(ms)
           call Pack_cc_idot(1,psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       sum)
!$OMP MASTER
           eorbit = eorbit + sum
!$OMP END MASTER
         end do
         end do
      end if
!$OMP MASTER
      if (ispin.eq.1) eorbit = eorbit+eorbit
!$OMP END MASTER
!$BARRIER
      call Parallel_SumAll(eorbit)

      electron_eorbit = eorbit 

      return
      end
 

*     ***********************************
*     *					*
*     *		electron_ehartree	*
*     *					*
*     ***********************************
      real*8 function electron_ehartree(dng)
      implicit none
      complex*16 dng(*)


*     **** external functions ****
      real*8   coulomb_e
      external coulomb_e

      electron_ehartree = coulomb_e(dng)

      return
      end

*     ***********************************
*     *					*
*     *		electron_ehartree2	*
*     *					*
*     ***********************************
      real*8 function electron_ehartree2(dn)
      implicit none
      real*8     dn(*)

#include "bafdecls.fh"
#include "electron_common.fh"

*     **** local variables ****
      real*8 ehartr, ehart1,ehart2,dv
      integer nx,ny,nz

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2

*     ***** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)

      call D3dB_rr_dot(1,dn(1),
     >                   dcpl_mb(vc(1)),
     >                   e1)
      call D3dB_rr_dot(1,dn(1+(ispin-1)*2*nfft3d),
     >                   dcpl_mb(vc(1)),
     >                   e2)
      ehartr = 0.5d0*(e1+e2)*dv

      electron_ehartree2 = ehartr

      return
      end

*     ***********************************
*     *					*
*     *		electron_exc		*
*     *					*
*     ***********************************
      real*8 function electron_exc(dnall)
      implicit none
      real*8 dnall(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer nx,ny,nz
      real*8  exc,exc2,dv

      common /eelectron_ejtmp/ exc,exc2

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      
      dv = lattice_omega()/dble(nx*ny*nz)

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dnall,
     >                 dbl_mb(xce(1)),
     >                 exc)
      if (ispin.eq.1) then
!$OMP MASTER
         exc= exc + exc 
!$OMP END MASTER
      else
         call D3dB_rr_dot(1,dnall(1+2*nfft3d),
     >                    dbl_mb(xce(1)),
     >                    exc2)
!$OMP MASTER
         exc= exc + exc2
!$OMP END MASTER
      end if
!$OMP MASTER
      exc = exc*dv
!$OMP END MASTER
!$OMP BARRIER
        
      electron_exc =  exc 

      return
      end


*     ***********************************
*     *					*
*     *		electron_pxc		*
*     *					*
*     ***********************************
      real*8 function electron_pxc(dn)
      implicit none
      real*8 dn(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer nx,ny,nz
      real*8  pxc,pxc2,dv

      common /eelectron_ejtmp/ pxc,pxc2

*     **** external functions *****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      
      dv = lattice_omega()/dble(nx*ny*nz)

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     >                 pxc)
      if (ispin.eq.1) then
!$OMP MASTER
         pxc= pxc + pxc 
!$OMP END MASTER
      else
         call D3dB_rr_dot(1,dn(1+2*nfft3d),
     >                    dbl_mb(xcp(1)+2*nfft3d),
     >                    pxc2)
!$OMP MASTER
         pxc= pxc + pxc2
!$OMP END MASTER
      end if
!$OMP MASTER
      pxc = pxc*dv
!$OMP END MASTER
!$OMP BARRIER
        
      electron_pxc =  pxc

      return
      end

*     ***********************************
*     *					*
*     *		electron_pxc_rho	*
*     *					*
*     ***********************************
      real*8 function electron_pxc_rho(rho)
      implicit none
      real*8 rho(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer nx,ny,nz
      real*8  pxc,pxc2,dv

*     **** external functions *****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      
      dv = lattice_omega()/dble(nx*ny*nz)

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,rho,
     >                 dbl_mb(xcp(1)),
     >                 pxc)
      if (ispin.eq.1) then
         pxc = pxc + pxc
      else
         call D3dB_rr_dot(1,rho,
     >                    dbl_mb(xcp(1)+2*nfft3d),
     >                    pxc2)
         pxc = (pxc + pxc2)
      end if
      pxc = pxc*dv
          
      electron_pxc_rho =  pxc
      return
      end


*     ***********************************
*     *                                 *
*     *         electron_xcp_ptr        *
*     *                                 *
*     ***********************************
      integer function electron_xcp_ptr()
      implicit none

#include "electron_common.fh"

      electron_xcp_ptr = xcp(1)
      return
      end


*     ***********************************
*     *                                 *
*     *         electron_psi_r_ptr      *
*     *                                 *
*     ***********************************
      integer function electron_psi_r_ptr()
      implicit none

#include "electron_common.fh"

      electron_psi_r_ptr = psi_r(1)
      return
      end

*     ***********************************
*     *                                 *
*     *         electron_Hpsi_k_ptr      *
*     *                                 *
*     ***********************************
      integer function electron_Hpsi_k_ptr()
      implicit none

#include "electron_common.fh"

      electron_Hpsi_k_ptr = Hpsi_k(1)
      return
      end



*     ***********************************
*     *					*
*     *	    electron_SIC_energies	*
*     *					*
*     ***********************************
      
      subroutine electron_SIC_energies(ehsic0,phsic0,exsic0,pxsic0)
      implicit none
      real*8 ehsic0,phsic0
      real*8 exsic0,pxsic0

#include "bafdecls.fh"
#include "electron_common.fh"

      logical  pspw_SIC_relaxed
      external pspw_SIC_relaxed
      
      if (.not.pspw_SIC_relaxed()) then
         call pspw_energy_SIC(ispin,dbl_mb(psi_r(1)),
     >                        ehsic,
     >                        phsic,
     >                        exsic,
     >                        pxsic)
         phsic = 0.0d0
         pxsic = 0.0d0
      end if

      ehsic0 = ehsic
      exsic0 = exsic
      phsic0 = phsic
      pxsic0 = pxsic
      return
      end



*     ***********************************
*     *                                 *
*     *      electron_SIC_stress        *
*     *                                 *
*     ***********************************
      subroutine electron_SIC_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "bafdecls.fh"
#include "electron_common.fh"


      call pspw_SIC_euv(ispin,dbl_mb(psi_r(1)),stress)
      return
      end


*     ***********************************
*     *                                 *
*     *      electron_HFX_stress        *
*     *                                 *
*     ***********************************
      subroutine electron_HFX_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "bafdecls.fh"
#include "electron_common.fh"


      call pspw_energy_euv_HFX(ispin,dbl_mb(psi_r(1)),stress)
      return
      end



*     ***********************************
*     *					*
*     *	    electron_HFX_energies	*
*     *					*
*     ***********************************
      
      subroutine electron_HFX_energies(ehfx0,phfx0)
      implicit none
      real*8 ehfx0,phfx0

#include "bafdecls.fh"
#include "electron_common.fh"

      logical  pspw_HFX_relaxed
      external pspw_HFX_relaxed
      
      if (.not.pspw_HFX_relaxed()) then
         call pspw_energy_HFX(ispin,dbl_mb(psi_r(1)),
     >                        ehfx,
     >                        phfx)
         phfx = 0.0d0
      end if

      ehfx0 = ehfx
      phfx0 = phfx
      return
      end


*     ***********************************
*     *                                 *
*     *     electron_U_energies         *
*     *                                 *
*     ***********************************

      subroutine electron_U_energies(edftu0,pdftu0)
      implicit none
      real*8 edftu0,pdftu0

#include "bafdecls.fh"
#include "electron_common.fh"


      edftu0 = edftu
      pdftu0 = pdftu
      return
      end


*     ***********************************
*     *                                 *
*     *     electron_meta_energies      *
*     *                                 *
*     ***********************************

      subroutine electron_meta_energies(emeta0,pmeta0)
      implicit none
      real*8 emeta0,pmeta0

#include "bafdecls.fh"
#include "electron_common.fh"

      emeta0 = emeta
      pmeta0 = pmeta
      return
      end


*     ***********************************
*     *                                 *
*     *      electron_apc_energies      *
*     *                                 *
*     ***********************************

      subroutine electron_apc_energies(eapc0,papc0)
      implicit none
      real*8 eapc0,papc0

#include "bafdecls.fh"
#include "electron_common.fh"

      eapc0 = eapc
      papc0 = papc
      return
      end


*     ***********************************
*     *					*
*     *		electron_get_Hpsi_k	*
*     *					*
*     ***********************************
      subroutine electron_get_Hpsi_k(Hpsi_k_new)
      implicit none
      complex*16 Hpsi_k_new(*)

#include "bafdecls.fh"
#include "electron_common.fh"

c      call dcopy(2*npack1*(neq(1)+neq(2)),
c     >           dcpl_mb(Hpsi_k(1)),1,
c     >           Hpsi_k_new,1)
      call Parallel_shared_vector_copy(.true.,
     >             2*npack1*(neq(1)+neq(2)),
     >             dcpl_mb(Hpsi_k(1)),
     >             Hpsi_k_new)
      return
      end



*     ***************************
*     *				*
*     *	   electron_ispin	*
*     *				*
*     ***************************
      integer function electron_ispin()
      implicit none

#include "electron_common.fh"

      electron_ispin = ispin
      return
      end


*     ***************************
*     *				*
*     *	     electron_ne	*
*     *				*
*     ***************************
      integer function electron_ne(ms)
      implicit none
      integer ms

#include "electron_common.fh"

      electron_ne = ne(ms)
      return
      end

*     ***************************
*     *				*
*     *	     electron_neq	*
*     *				*
*     ***************************
      integer function electron_neq(ms)
      implicit none
      integer ms

#include "electron_common.fh"

      electron_neq = neq(ms)
      return
      end


*     ***********************************
*     *					*
*     *	    electron_get_Tgradient 	*
*     *					*
*     ***********************************

      subroutine electron_get_Tgradient(psi_k,hml,THpsi_k)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     ***** local variables ****
      integer ms,n,shift1,shift2
     
      call  Dneall_fmf_Multiply(0,
     >                           psi_k,npack1,
     >                           hml,    1.0d0,
     >                           THpsi_k,0.0d0)
      call DAXPY_OMP(2*npack1*(neq(1)+neq(2)),
     >           (-1.0d0),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           THpsi_k,1)
      return
      end


*     ***********************************
*     *					*
*     *	    electron_gen_Tangent 	*
*     *					*
*     ***********************************

      subroutine electron_gen_Tangent(psi_k,hml,THpsi_k)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"

*     ***** local variables ****
c      integer ms,n,shift1,shift2
     
      call Dneall_fmf_Multiply(0,psi_k,npack1,
     >                            hml,1.0d0,
     >                            THpsi_k,-1.0d0)

c      do ms=1,ispin
c         n     = ne(ms)
c         if (n.le.0) go to 30
c         shift1 = 1 + (ms-1)*ne(1)*npack1
c         shift2 = 1 + (ms-1)*ne(1)*ne(1)
c         call DGEMM('N','N',2*npack1,n,n,
c     >             (1.0d0),
c     >             psi_k(shift1),  2*npack1,
c     >             hml(shift2),    n,
c     >             (-1.0d0),
c     >             THpsi_k(shift1),2*npack1)
c   30    continue
c      end do

      return
      end


*     ***********************************
*     *					*
*     *	    electron_get_Gradient 	*
*     *					*
*     ***********************************

      subroutine electron_get_Gradient(THpsi_k)
      implicit none
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"


c      call dcopy(2*npack1*(neq(1)+neq(2)),
c     >           dcpl_mb(Hpsi_k(1)),1,
c     >           THpsi_k,1)      
      call Parallel_shared_vector_copy(.true.,
     >             2*npack1*(neq(1)+neq(2)),
     >             dcpl_mb(Hpsi_k(1)),
     >             THpsi_k)      
      return
      end


*     ***********************************
*     *					*
*     *	    electron_get_gradient_orb 	*
*     *					*
*     ***********************************

      subroutine electron_get_gradient_orb(i,Horb)
      implicit none
      integer i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "electron_common.fh"

      call Pack_c_Copy(1,dcpl_mb(Hpsi_k(1)+(i-1)*npack1),Horb) 
      
      return
      end



*     ***********************************
*     *					*
*     *	    electron_get_TMgradient 	*
*     *					*
*     ***********************************

      subroutine electron_get_TMgradient(psi_k,THpsi_k)
      implicit none
#include "errquit.fh"
      complex*16 psi_k(*)
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     ***** local variables ****
      logical value
      integer ms,n,n1(2),shift,mhml(2)


      n1(1) = 1
      n1(2) = ne(1)+1

      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'mhml',mhml(2),mhml(1))
      if (.not. value)
     >   call errquit('electron_get_Tgradient: push stack',0, MA_ERR)



*     **** generate M*H|psi> *****
      call Grsm_gg_Copy(npack1,(neq(1)+neq(2)),
     >                  dcpl_mb(Hpsi_k(1)),
     >                  THpsi_k)

      call ke_Precondition(npack1,(neq(1)+neq(2)),
     >                  psi_k,
     >                  THpsi_k)

*     **** generate mhml = <psi|M*H|psi> ****
      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         n     = ne(ms)
         call Grsm_ggm_dot(npack1,n,
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
     >                     dbl_mb(mhml(1)+shift))
      end do

*     **** mhml = -mhml ****
      call DSCAL_OMP(2*ne(1)*ne(1),(-1.0d0),dbl_mb(mhml(1)),1)


*     **** generate TMG = M*H|psi> - |psi>*mhml ****
      do ms=1,ispin
            shift = (ms-1)*ne(1)*ne(1)
            n     = ne(ms)
            call Grsm_gmg_daxpy(npack1,n,
     >                        psi_k(1+(ms-1)*ne(1)*npack1),
     >                        dbl_mb(mhml(1)+shift),
     >                        THpsi_k(1+(ms-1)*ne(1)*npack1))
      end do

      call Grsm_gg_dScale1(npack1,(neq(1)+neq(2)),
     >                    (-1.0d0),
     >                    THpsi_k)


      value = BA_pop_stack(mhml(2))
      if (.not. value)
     > call errquit('electron_get_Tradient: popping stack',0, MA_ERR)

      return
      end


*     ***************************
*     *				*
*     *	    electron_gen_hml 	*
*     *				*
*     ***************************

      subroutine electron_gen_hml(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "electron_common.fh"


c*     **** local variables ****
c      integer ms,n,n1(2),shift
     
      call Dneall_ffm_sym_Multiply(0,psi_k,
     >                                dcpl_mb(Hpsi_k(1)),npack1,
     >                                hml)

c      n1(1) = 1
c      n1(2) = ne(1) + 1

c      do ms=1,ispin
c         shift = (ms-1)*ne(1)*ne(1)
c         n     = ne(ms)
c         if (n.le.0) go to 30
cc        call Grsm_ggm_sym_dot(npack1,n,
cc    >                     psi_k(1+(ms-1)*ne(1)*npack1),
cc    >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
cc    >                     hml(shift+1))
c         call Pack_ccm_sym_dot(1,n,
c     >                     psi_k(1+(ms-1)*ne(1)*npack1),
c     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
c     >                     hml(shift+1))
c  30     continue
c      end do

      return
      end




*     ***************************
*     *                         *
*     *     electron_gen_hml_g  *
*     *                         *
*     ***************************

      subroutine electron_gen_hml_g(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "electron_common.fh"


c*     **** local variables ****
c      integer ms,n,n1(2),shift
     
      call Dneall_ffm_Multiply(0,psi_k,
     >                           dcpl_mb(Hpsi_k(1)),npack1,
     >                           hml)

c      n1(1) = 1
c      n1(2) = ne(1) + 1

c      do ms=1,ispin
c         n     = ne(ms)
c         if (n.le.0) go to 30
c         shift = (ms-1)*ne(1)*ne(1)
c         call Pack_ccm_dot(1,n,
c     >                     psi_k(1+(ms-1)*ne(1)*npack1),
c     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
c     >                     hml(shift+1))
c   30    continue
c      end do
      return
      end


*     ***********************************
*     *				        *
*     *	    electron_gen_hmlt     	*
*     *				        *
*     ***********************************

      subroutine electron_gen_hmlt(psi_k,hmlt)
      implicit none
      complex*16 psi_k(*)
      real*8     hmlt(*)

#include "bafdecls.fh"
#include "electron_common.fh"


c*     **** local variables ****
c      integer ms,n,n1(2),shift
     
      call Dneall_ffm_Multiply(0,dcpl_mb(Hpsi_k(1)),
     >                            psi_k,npack1,
     >                            hmlt)

c      n1(1) = 1
c      n1(2) = ne(1) + 1

c      do ms=1,ispin
c         n     = ne(ms)
c         if (n.le.0) go to 30
c         shift = (ms-1)*ne(1)*ne(1)
cc        call Pack_ccm_sym_dot(1,n,
cc    >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
cc    >                     psi_k(1+(ms-1)*ne(1)*npack1),
cc    >                     hmlt(shift+1))
c         call Pack_ccm_dot(1,n,
c     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
c     >                     psi_k(1+(ms-1)*ne(1)*npack1),
c     >                     hmlt(shift+1))
c   30    continue
c      end do


      return
      end



*     **************************************
*     *				           *
*     *	    electron_gen_psiTangenthml 	   *
*     *				           *
*     **************************************

      subroutine electron_gen_psiTangenthml(psi_k,THpsi_k,hml)
      implicit none
      complex*16 psi_k(*)
      complex*16 THpsi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "electron_common.fh"


c*     **** local variables ****
c      integer ms,n,shift
     
      call Dneall_ffm_sym_Multiply(0,psi_k,
     >                                THpsi_k,npack1,
     >                                hml)

c      do ms=1,ispin
c         n     = ne(ms)
c         if (n.le.0) go to 30
c         shift = (ms-1)*ne(1)*ne(1)
c         call Pack_ccm_sym_dot(1,n,
c     >                     psi_k(1+(ms-1)*ne(1)*npack1),
c     >                     THpsi_k(1+(ms-1)*ne(1)*npack1),
c     >                     hml(shift+1))
c   30    continue
c      end do

      return
      end






**************************************************************************
**************************************************************************
*******    routines below this line are for internal use only    *********
**************************************************************************
**************************************************************************

*     ***********************************
*     *					*
*     *		electron_gen_Hpsi_k	*
*     *					*
*     ***********************************

      subroutine electron_gen_Hpsi_k(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"
cccc#include "frac_occ.fh"


*     **** local variables ****
      logical move,fractional
      integer n
      real*8  fion(3,1)

*     **** external functions ****
      integer  control_version
      external control_version

      move = .false.
      fractional = .false.
*     ******************
*     **** get Hpsi ****
*     ******************
      if (control_version().eq.3)
     >  call psi_H(ispin,neq,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(vc(1)),
     >             dbl_mb(xcp(1)),
     >             dcpl_mb(Hpsi_k(1)),
     >             move,
     >             fion,
     >             fractional,fion)

      if (control_version().eq.4)
     >  call psi_Hv4(ispin,neq,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dcpl_mb(vl(1)),
     >             dbl_mb(vl_lr(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(vc(1)),
     >             dbl_mb(xcp(1)),
     >             dcpl_mb(Hpsi_k(1)),
     >             move,     
     >             fion,
     >             fractional,fion)

      call Grsm_gg_dScale1(npack1,(neq(1)+neq(2)),(-1.0d0),
     >                     dcpl_mb(Hpsi_k(1)))

      return
      end


*     ***********************************
*     *                                 *
*     *     electron_gen_Hvallpsi_k     *
*     *                                 *
*     ***********************************

      subroutine electron_gen_Hvallpsi_k(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      logical move,fractional
      integer n
      real*8  fion(3,1)

*     **** external functions ****
      integer  control_version
      external control_version

      move = .false.
      fractional = .false.
*     ******************
*     **** get Hpsi ****
*     ******************
      if (control_version().eq.3)
     >  call psi_H_vall(ispin,neq,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dbl_mb(vall(1)),
     >             dcpl_mb(Hpsi_k(1)))

      if (control_version().eq.4)
     >  call psi_Hv4_vall(ispin,neq,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dbl_mb(vall(1)),
     >             dcpl_mb(Hpsi_k(1)))

      call Grsm_gg_dScale1(npack1,(neq(1)+neq(2)),(-1.0d0),
     >                     dcpl_mb(Hpsi_k(1)))

      return
      end




*     ***********************************
*     *					*
*     *     electron_gen_Hpsi_k_orb	*
*     *					*
*     ***********************************

      subroutine electron_gen_Hpsi_k_orb(n,psi_k)
      implicit none
      integer n
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer ms,index1,index1r,index2

*     **** external functions ****
      integer  control_version
      external control_version

     
      if (n.le.neq(1)) then
        ms=1
      else
        ms=2
      end if 
      index1  = (n-1)*nfft3d  
      index1r = 2*index1
      index2  = (n-1)*npack1

*     ******************
*     **** get Hpsi ****
*     ******************
      if (control_version().eq.3)
     >  call psi_Horb(.true.,ispin,ms,
     >             dbl_mb(psi_r(1)),
     >             dbl_mb(vall(1)),
     >             psi_k(index2+1),
     >             dbl_mb(  psi_r(1)+index1r),
     >             dcpl_mb(Hpsi_k(1)+index2))

      if (control_version().eq.4)
     >  call psi_Horbv4(.true.,ispin,ms,
     >             dbl_mb(psi_r(1)),
     >             dbl_mb(vall(1)),
     >             psi_k(index2+1),
     >             dbl_mb(  psi_r(1)+index1r),
     >             dcpl_mb(Hpsi_k(1)+index2))

c      call Pack_c_SMul(1,(-1.0d0),
c     >                 dcpl_mb(Hpsi_k(1)+index2),
c     >                 dcpl_mb(Hpsi_k(1)+index2))
      call Pack_c_SMul1(1,(-1.0d0),
     >                 dcpl_mb(Hpsi_k(1)+index2))

      return
      end

*     ***********************************
*     *					*
*     *   electron_get_gradient_virtual *
*     *					*
*     ***********************************

      subroutine electron_get_gradient_virtual(ms,orb,Horb)
      implicit none
      integer    ms
      complex*16 orb(*)
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value,hfxon
      integer n2ft3d
      integer tmp_r(2)

*     **** external functions ****
      integer  control_version
      external control_version
      logical  pspw_HFX_virtual
      external pspw_HFX_virtual


      hfxon = pspw_HFX_virtual()
      n2ft3d = 2*nfft3d

      value = BA_push_get(mt_dbl,(n2ft3d),'tmp_r',tmp_r(2),tmp_r(1))
      if (.not. value)
     >   call errquit('electron_get_gradient_virtual: push stack',0,
     &       MA_ERR)


      call Pack_c_Copy(1,orb,dbl_mb(tmp_r(1)))
      call Pack_c_unpack(1,  dbl_mb(tmp_r(1)))
      call D3dB_cr_pfft3b(1,1,dbl_mb(tmp_r(1))) 

*     **** get Hpsi ****
      if (control_version().eq.3)
     >  call psi_Horb_replicated(hfxon,ispin,ms,
     >             dbl_mb(psi_r(1)),
     >             dbl_mb(vall(1)),
     >             orb,
     >             dbl_mb(tmp_r(1)),
     >             Horb)

      if (control_version().eq.4)
     >  call psi_Horbv4_replicated(hfxon,ispin,ms,
     >             dbl_mb(psi_r(1)),
     >             dbl_mb(vall(1)),
     >             orb,
     >             dbl_mb(tmp_r(1)),
     >             Horb)

      call Pack_c_SMul1(1,(-1.0d0),Horb)

      value = BA_pop_stack(tmp_r(2))
      if (.not. value) call errquit(
     >     'electron_get_gradient_virtual: poping stack',1, MA_ERR)


      return
      end

*     ***********************************
*     *                                 *
*     *    electron_get_H0psi_k_orb     *
*     *                                 *
*     ***********************************
      subroutine electron_get_H0psi_k_orb(orb,Horb)
      implicit none
      complex*16 orb(*)
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer n2ft3d
      integer tmp_r(2)

*     **** external functions ****
      integer  control_version
      external control_version

      n2ft3d = 2*nfft3d

      value = BA_push_get(mt_dbl,n2ft3d,'tmp_r',tmp_r(2),tmp_r(1))
      if (.not.value)
     >   call errquit('electron_get_H0psi_k_orb: push stack',0,MA_ERR)


      call Pack_c_Copy(1,orb,dbl_mb(tmp_r(1)))
      call Pack_c_unpack(1,  dbl_mb(tmp_r(1)))
      call D3dB_cr_pfft3b(1,1,dbl_mb(tmp_r(1)))

*     **** get Hpsi ****
      call psi_H0orb_replicated(
     >             dbl_mb(vall(1)),
     >             orb,
     >             dbl_mb(tmp_r(1)),
     >             Horb)

      call Pack_c_SMul1(1,(-1.0d0),Horb)

      value = BA_pop_stack(tmp_r(2))
      if (.not.value) call errquit(
     >     'electron_get_H0psi_k_orb: poping stack',1,MA_ERR)

      return
      end



*     ***************************
*     *				*
*     *	    electron_gen_psi_r	*
*     *				*
*     ***************************

      subroutine electron_gen_psi_r(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"

*     **** local variables ****
      integer n,nemax,n2ft3d
      integer i,taskid
      call Parallel_taskid(taskid)


*     ***** generate compensation charge ****
      n2ft3d = 2*nfft3d
      nemax = neq(1) + neq(2)

c     call Grsm_gg_Copy(npack1,nemax,psi_k,dbl_mb(psi_r(1)))
!$OMP DO private(n)
      do n=1,nemax
         call Pack_c_Copy0(1,psi_k(1+(n-1)*npack1),
     >                    dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      end do
!$OMP END DO

      call Grsm_gh_fftb(nfft3d,nemax,dbl_mb(psi_r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemax,dbl_mb(psi_r(1)))  !*** probably not neeeded!

*     ****  generate tau functions ****
      call nwpw_meta_gga_gen_tau(ispin,neq,psi_k)

      return
      end

*     ***********************************
*     *                                 *
*     *      electron_gen_psi_r_orb     *
*     *                                 *
*     ***********************************

      subroutine electron_gen_psi_r_orb(n,psi_k)
      implicit none
      integer    n
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "electron_common.fh"

*     **** local variables ****
      integer n2ft3d

      n2ft3d = 2*nfft3d

      call Pack_c_Copy(1,psi_k(1+(n-1)*npack1),
     >                        dbl_mb(psi_r(1)+(n-1)*n2ft3d))
     
      call Pack_c_unpack(1,   dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      !call D3dB_cr_fft3b(1,   dbl_mb(psi_r(1)+(n-1)*n2ft3d)) 
      call D3dB_cr_pfft3b(1,1,dbl_mb(psi_r(1)+(n-1)*n2ft3d)) 
      call D3dB_r_Zero_Ends(1,dbl_mb(psi_r(1)+(n-1)*n2ft3d))

      return
      end

*     ***************************
*     *				*
*     *	 electron_gen_density	*
*     *				*
*     ***************************

      subroutine electron_gen_density(psi_k,dn,fractional,occ)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "electron_common.fh"
cccc#include "frac_occ.fh"

     
*     **** local variables ****
      integer i
      integer ms,n2ft3d
      integer n,n1(2),n2(2)
      real*8  scal2,wf
      integer tmp1(2)
      logical value

*     ***** external functions *****
      logical  psp_semicore
      real*8   lattice_omega
      external psp_semicore
      external lattice_omega


      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)

      n2ft3d = 2*nfft3d
      scal2 = 1.0d0/lattice_omega()


*     *********************
*     **** generate dn ****
*     *********************
c      call dcopy(2*n2ft3d,0.0d0,0,dn,1)
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,dn)
      if (fractional) then
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            wf = occ(n)
            do i=1,n2ft3d
               dn(i+(ms-1)*n2ft3d)
     >            = dn(i+(ms-1)*n2ft3d)
     >            + wf*scal2*(dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d)**2)
            end do
         end do
         call D3dB_r_Zero_Ends(1,dn(1+(ms-1)*n2ft3d))
         call D1dB_Vector_SumAll(n2ft3d,dn(1+(ms-1)*n2ft3d))
      end do
      
      else
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            do i=1,n2ft3d
               dn(i+(ms-1)*n2ft3d) 
     >            = dn(i+(ms-1)*n2ft3d) 
     >            + scal2*(dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d)**2)
            end do
         end do
         call D3dB_r_Zero_Ends(1,dn(1+(ms-1)*n2ft3d))
         call D1dB_Vector_SumAll(n2ft3d,dn(1+(ms-1)*n2ft3d))
      end do
      end if
      call rho_symmetrizer_dn(ispin,n2ft3d,dn)

      return
      end



*     ***************************
*     *				*
*     *	 electron_gen_densities	*
*     *				*
*     ***************************

      subroutine electron_gen_densities(psi_k,dn,dng,dnall,
     >                                  fractional,occ)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)
      logical fractional
      real*8     occ(*)

#include "bafdecls.fh"
#include "electron_common.fh"
cccccccccc#include "frac_occ.fh"

     
*     **** local variables ****
      integer i
      integer ms,n2ft3d
      integer n,n1(2),n2(2)
      real*8  scal2,wf
      integer tmp1(2)
      logical value

*     ***** external functions *****
      logical  psp_semicore
      real*8   lattice_omega
      external psp_semicore
      external lattice_omega


      n1(1) = 1
      n1(2) = neq(1) + 1
      n2(1) = neq(1)
      n2(2) = neq(1) + neq(2)

      n2ft3d = 2*nfft3d
      scal2 = 1.0d0/lattice_omega()


*     *********************
*     **** generate dn ****
*     *********************

      !call dcopy(2*n2ft3d,0.0d0,0,dn,1)
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,dn)
      if (fractional) then
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            wf = occ(n)
!$OMP DO private(i)
            do i=1,n2ft3d
               dn(i+(ms-1)*n2ft3d)
     >            = dn(i+(ms-1)*n2ft3d)
     >            + wf*scal2*(dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d)**2)
            end do
!$OMP END DO
         end do
         call D3dB_r_Zero_Ends(1,dn(1+(ms-1)*n2ft3d)) !*** probably not needed!
         call D1dB_Vector_SumAll(n2ft3d,dn(1+(ms-1)*n2ft3d))
      end do
      else
      do ms=1,ispin
         do n=n1(ms),n2(ms)
!$OMP DO private(i)
            do i=1,n2ft3d
               dn(i+(ms-1)*n2ft3d) 
     >            = dn(i+(ms-1)*n2ft3d) 
     >            + scal2*(dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d)**2)
            end do
!$OMP END DO
         end do
         call D3dB_r_Zero_Ends(1,dn(1+(ms-1)*n2ft3d))   !*** probably not needed!
         call D1dB_Vector_SumAll(n2ft3d,dn(1+(ms-1)*n2ft3d))
      end do
      end if
      call rho_symmetrizer_dn(ispin,n2ft3d,dn)


*     **** generate dng and dnall ****
      call electron_gen_dng_dnall(dn,dng,dnall)

      return
      end



*     ***************************
*     *				*
*     *	 electron_gen_dng_dnall	*
*     *				*
*     ***************************

      subroutine electron_gen_dng_dnall(dn,dng,dnall)
      implicit none
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"

     
*     **** local variables ****
      integer i
      integer ms,nx,ny,nz,n2ft3d
      integer n,n1(2),n2(2)
      real*8  scal1
      integer tmp1(2)
      logical value

*     ***** external functions *****
      logical  psp_semicore
      external psp_semicore

      n2ft3d = 2*nfft3d
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
        
*     **********************
*     **** generate dng ****
*     **********************
      value = BA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_dng_dnall: out of stack memory',0, MA_ERR)

      call D3dB_rr_Sum(1,dn,dn(1+(ispin-1)*n2ft3d),dbl_mb(tmp1(1)))
c      call D3dB_r_SMul(1,scal1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
      call D3dB_r_SMul1(1,scal1,dbl_mb(tmp1(1)))
      !call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
      call D3dB_rc_pfft3f(1,0,dbl_mb(tmp1(1)))
      call Pack_c_pack(0,dbl_mb(tmp1(1)))
      call Pack_c_Copy(0,dbl_mb(tmp1(1)),dng)

*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (psp_semicore(0)) then
           call semicore_density(dbl_mb(tmp1(1)))
           call D3dB_r_SMul1(1,0.5d0,dbl_mb(tmp1(1)))
        else
           call D3dB_r_Zero(1,dbl_mb(tmp1(1)))
        end if
        do ms=1,ispin
          call D3dB_rr_Sum(1,dn(1+(ms-1)*n2ft3d),
     >                     dbl_mb(tmp1(1)),
     >                     dnall(1+(ms-1)*n2ft3d))
        end do


      value = BA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >     'electron_gen_dng_dnall: poping stack',1, MA_ERR)


      return
      end


*     ***********************************
*     *			 		*
*     *        electron_gen_vall        *
*     *					*
*     ***********************************

      subroutine electron_gen_vall()
      implicit none
   

#include "bafdecls.fh"
#include "electron_common.fh"
     
*     **** local variables ****
      integer ms
      real*8 scal2
      
*     **** external functions ****
      integer  control_version
      real*8   lattice_omega
      external control_version
      external lattice_omega
      

      scal2 = 1.0d0/lattice_omega()

      if (control_version().eq.3) then

*       **** add up k-space potentials, vall = scal2*vl + vc  ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))

c        call Pack_cc_Sum(0,dbl_mb(vall(1)),
c     >                     dcpl_mb(vc(1)),
c     >                     dbl_mb(vall(1)))
        call Pack_cc_Sum2(0,dcpl_mb(vc(1)),dbl_mb(vall(1)))
      
*       **** fourier transform k-space potentials ****
        call Pack_c_unpack(0,dbl_mb(vall(1)))
        !call D3dB_cr_fft3b(1,dbl_mb(vall(1)))
        call D3dB_cr_pfft3b(1,0,dbl_mb(vall(1)))

*       **** add v_field to vall ****
c        if (field_exist)
c     >    call D3dB_rr_Sum(1,dbl_mb(vall(1)),
c     >                       dbl_mb(v_field(1)),
c     >                       dbl_mb(vall(1)))
        if (field_exist)
     >    call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),dbl_mb(vall(1)))

      else

*       **** add up k-space potentials, vall = scal2*vsr_l    ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))
      
*        **** fourier transform k-space potentials ****
         call Pack_c_unpack(0,dbl_mb(vall(1)))
         !call D3dB_cr_fft3b(1,dbl_mb(vall(1)))
         call D3dB_cr_pfft3b(1,0,dbl_mb(vall(1)))

*        **** add vlr_l, vc and v_field to vall ****
c         call D3dB_rr_Sum(1,dbl_mb(vall(1)),
c     >                      dbl_mb(vl_lr(1)),
c     >                      dbl_mb(vall(1)))
c         call D3dB_rr_Sum(1,dbl_mb(vall(1)),
c     >                      dcpl_mb(vc(1)),   
c     >                      dbl_mb(vall(1)))
c         if (field_exist)
c     >     call D3dB_rr_Sum(1,dbl_mb(vall(1)),
c     >                        dbl_mb(v_field(1)),
c     >                        dbl_mb(vall(1)))
         call D3dB_rr_Sum2(1,dbl_mb(vl_lr(1)),dbl_mb(vall(1)))
         call D3dB_rr_Sum2(1,dcpl_mb(vc(1)),dbl_mb(vall(1)))
         if (field_exist)
     >     call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),dbl_mb(vall(1)))

         if (confine_exist)
     >     call D3dB_rr_Sum2(1,dbl_mb(v_confine(1)),dbl_mb(vall(1)))

      end if

*     **** add xcp to vall ****
c      do ms=ispin,1,-1      
c        call D3dB_rr_Sum(1,dbl_mb(vall(1)),
c     >                   dbl_mb(xcp(1)+(ms-1)*2*nfft3d),
c     >                   dbl_mb(vall(1)+(ms-1)*2*nfft3d))
c        call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)+(ms-1)*2*nfft3d))
c      end do

      if (ispin.eq.2) then
        call D3dB_rr_Sum(1,dbl_mb(vall(1)),
     >                   dbl_mb(xcp(1) +2*nfft3d),
     >                   dbl_mb(vall(1)+2*nfft3d))
        call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)+2*nfft3d))
      end if
      call D3dB_rr_Sum2(1,dbl_mb(xcp(1)),dbl_mb(vall(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)))

      return
      end


*     ***********************************
*     *                                 *
*     *        electron_gen_vall0       *
*     *                                 *
*     ***********************************

*     Generates non self-consistent potentials only.

      subroutine electron_gen_vall0()
      implicit none

#include "bafdecls.fh"
#include "electron_common.fh"

*     **** local variables ****
      integer ms
      real*8 scal2

*     **** external functions ****
      integer  control_version
      real*8   lattice_omega
      external control_version
      external lattice_omega


      scal2 = 1.0d0/lattice_omega()


      if (control_version().eq.3) then
*       **** add up k-space potentials, vall = scal2*vl  ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))

*       **** fourier transform k-space potentials ****
        call Pack_c_unpack(0,dbl_mb(vall(1)))
        call D3dB_cr_pfft3b(1,0,dbl_mb(vall(1)))

        if (field_exist)
     >    call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),dbl_mb(vall(1)))

      else
*       **** add up k-space potentials, vall = scal2*vsr_l    ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))

*        **** fourier transform k-space potentials ****
         call Pack_c_unpack(0,dbl_mb(vall(1)))
         call D3dB_cr_pfft3b(1,0,dbl_mb(vall(1)))

         call D3dB_rr_Sum2(1,dbl_mb(vl_lr(1)),dbl_mb(vall(1)))
         if (field_exist)
     >     call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),dbl_mb(vall(1)))

         if (confine_exist)
     >     call D3dB_rr_Sum2(1,dbl_mb(v_confine(1)),dbl_mb(vall(1)))

      end if
      if (ispin.eq.2) then
         call D3dB_r_Copy(1,dbl_mb(vall(1)),dbl_mb(vall(1)+2*nfft3d))
      end if

      return 
      end 





*     ***********************************
*     *			 		*
*     *     electron_add_oep_to_vall	*
*     *					*
*     ***********************************

      subroutine electron_add_oep_to_vall(dn)
      implicit none
#include "errquit.fh"
      real*8 dn(*)

#include "bafdecls.fh"
#include "electron_common.fh"

      logical value
      integer ms,n2ft3d,v_oep(2)

      call D3dB_n2ft3d(1,n2ft3d)
      value = BA_push_get(mt_dbl,(2*n2ft3d),'V_OEP',v_oep(2),v_oep(1))
      if (.not. value)
     >  call errquit('electron_add_oep_to_vall:out of stack memory',0,0)
      

      call pspw_potential_SIC_OEP(ispin,ne,
     >                            dn,
     >                            dbl_mb(psi_r(1)),
     >                            dbl_mb(v_oep(1)))

*     **** add v_oep to vall ****
      do ms=1,ispin
c        call D3dB_rr_Sum(1,dbl_mb( vall(1)+(ms-1)*n2ft3d),
c     >                     dbl_mb(v_oep(1)+(ms-1)*n2ft3d),
c     >                     dbl_mb( vall(1)+(ms-1)*n2ft3d))
        call D3dB_rr_Sum2(1,dbl_mb(v_oep(1)+(ms-1)*n2ft3d),
     >                      dbl_mb( vall(1)+(ms-1)*n2ft3d))
        call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)+(ms-1)*n2ft3d))
      end do

      value = BA_pop_stack(v_oep(2))
      if (.not. value)
     >  call errquit(
     >  'electron_add_oep_to_vall:popping stack memory',1,0)

      return
      end



*     ***********************************
*     *			 		*
*     *        electron_get_vall	*
*     *					*
*     ***********************************

      subroutine electron_get_vall(vall_out)
      implicit none
      real*8 vall_out(*)

#include "bafdecls.fh"
#include "electron_common.fh"

      !call dcopy(4*nfft3d,dbl_mb(vall(1)),1,vall_out,1)
      call Parallel_shared_vector_copy(.true.,
     >             4*nfft3d,dbl_mb(vall(1)),vall_out)
      return
      end


*     ***********************************
*     *			 		*
*     *        electron_set_vall	*
*     *					*
*     ***********************************

      subroutine electron_set_vall(vall_in)
      implicit none
      real*8 vall_in(*)

#include "bafdecls.fh"
#include "electron_common.fh"

      !call dcopy(4*nfft3d,vall_in,1,dbl_mb(vall(1)),1)
      call Parallel_shared_vector_copy(.true.,
     >             4*nfft3d,vall_in,dbl_mb(vall(1)))
      return
      end


*     ***********************************
*     *			 		*
*     *   electron_gen_scf_potentials	*
*     *					*
*     ***********************************

      subroutine electron_gen_scf_potentials(dn,dng,dnall)
      implicit none
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "electron_common.fh"


*     ***** local variables ****
      logical value
      integer n2ft3d,gga
      integer tmp1(2),r_grid(2)
      real*8  gmma,ftmp(3)

*     **** external functions ****
      integer  control_gga,control_version
      external control_gga,control_version
      real*8   control_attenuation
      external control_attenuation
      logical  nwpw_cosmo1_on,nwpw_cosmo2_on,pspw_V_APC_on
      external nwpw_cosmo1_on,nwpw_cosmo2_on,pspw_V_APC_on

      n2ft3d = 2*nfft3d
      gga = control_gga()

      
      if (control_version().eq.3) then
            call coulomb_v(dng,dcpl_mb(vc(1)))
      end if

      if (control_version().eq.4)  then
         value = BA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >   'electron_gen_scf_potentials: out of stack memory',0,MA_ERR)
      
         call D3dB_rr_Sum(1,dn(1),dn(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(tmp1(1)))

         call coulomb2_v(dbl_mb(tmp1(1)),dcpl_mb(vc(1)))

         if (nwpw_cosmo1_on()) then
            if (.not.BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',
     >                           r_grid(2),r_grid(1)))
     >         call errquit("electron_gen_scf_potentials: push stack",
     >                      1,MA_ERR)
            call lattice_r_grid(dbl_mb(r_grid(1)))
            call v_local_set_cosmo_BQ(dbl_mb(r_grid(1)),
     >                                dbl_mb(tmp1(1)),dng)
            call electron_gen_vl_potential()
            if (.not.BA_pop_stack(r_grid(2)))
     >         call errquit("electron_gen_scf_potentials:pop stack",
     >                      1,MA_ERR)
         end if

c         if (pspw_V_APC_on()) then
c            call electron_gen_vl_potential()
c            call pspw_V_APC(ispin,ne,dng,dcpl_mb(vl(1)),eapc,papc)
c         end if

         value = BA_pop_stack(tmp1(2)) 
         if (.not. value) call errquit(
     >   'electron_gen_scf_potentials: error popping stack memory',0,
     >       MA_ERR)
      end if

      if (pspw_V_APC_on()) then
         call electron_gen_vl_potential()
         call pspw_V_APC(ispin,ne,dng,dcpl_mb(vl(1)),eapc,papc,
     >                   .false.,ftmp)
      end if

*    **** xc potential ****
      call v_bwexc_all(gga,n2ft3d,ispin,dnall,
     >                 dbl_mb(xcp(1)),dbl_mb(xce(1)))


      return
      end


*     ***********************************
*     *                                 *
*     *   electron_zero_scf_potentials  *
*     *                                 *
*     ***********************************

      subroutine electron_zero_scf_potentials()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "electron_common.fh"

*     ***** local variables ****
      integer n2ft3d

*     **** external functions ****
      integer  control_version
      external control_version

      n2ft3d = 2*nfft3d

*     **** zero out scf potentials ****
      if (control_version().eq.3) then
         call Parallel_shared_vector_zero(.true., 
     >                                    2*npack0,dcpl_mb(vc(1)))
      else
         call Parallel_shared_vector_zero(.true.,n2ft3d,dcpl_mb(vc(1)))
      end if

      call Parallel_shared_vector_zero(.true.,2*n2ft3d,dbl_mb(xce(1)))
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,dbl_mb(xcp(1)))

      return
      end



*     ***********************************
*     *			 		*
*     *   electron_gen_vl_potential 	*
*     *					*
*     ***********************************

      subroutine electron_gen_vl_potential()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      logical move,value
      integer n2ft3d
      integer tmp1(2)
      integer tmp2(2)
      integer r_grid(2)

*     **** external functions *****
      logical  pspw_charge_found,pspw_Efield_found
      external pspw_charge_found,pspw_Efield_found
      integer  control_version
      external control_version

      field_exist = pspw_charge_found().or.pspw_Efield_found()
      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        BA_push_get(mt_dbl,(3),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call 
     >   errquit('electron_gen_vl_potential: out of stack memory',0,
     &       MA_ERR)

      move = .false.
      call v_local(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(tmp1(1)),
     >               dbl_mb(tmp2(1)))

      value = BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >  'electron_gen_vl_potential: error popping stack memory',0,
     &       MA_ERR)

*     **** generate real-space fields ****
      if ((control_version().eq.4).or.field_exist) then

         value = BA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                       r_grid(2),r_grid(1)) 

         call lattice_r_grid(dbl_mb(r_grid(1)))

*        **** generate long-range psp potential ****
         if (control_version().eq.4) then
            call v_lr_local(dbl_mb(r_grid(1)),
     >                      dbl_mb(vl_lr(1)))
         end if

*        **** zero out v_field ****
         !call dcopy(2*nfft3d,0.0d0,0,dbl_mb(v_field(1)),1)
         call Parallel_shared_vector_zero(.true.,
     >               2*nfft3d,dbl_mb(v_field(1)))

*        **** generate other realspace fields ****
         if (field_exist) then
           field_exist = .true.
           n2ft3d = 2*nfft3d
*          **** generate charge potential ****
           call pspw_charge_Generate_V(n2ft3d,
     >                                 dbl_mb(r_grid(1)),
     >                                 dbl_mb(v_field(1)))
*          **** generate Efield potential ****
           call pspw_Efield_Generate_V(n2ft3d,
     >                                    dbl_mb(r_grid(1)),
     >                                    dbl_mb(v_field(1)))
         end if

         value = BA_pop_stack(r_grid(2))
         if (.not. value) call errquit(
     >   'electron_gen_vl_potential: error popping stack memory',0,
     &       MA_ERR)

      end if
 


      return
      end



*     ***********************************
*     *			 		*
*     *   electron_psi_vl_ave	 	*
*     *					*
*     ***********************************

      real*8 function electron_psi_vl_ave(psi1,dn)
      implicit none
      complex*16 psi1(*)
      real*8     dn(*)

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2)
      integer nx,ny,nz,n2ft3d,np
      real*8 elocal,sum,scal1,scal2,dv
      integer tmp1(2),tmp2(2)

      common /eelectron_ejtmp/ sum,elocal


*     **** external functions ***
      integer  control_version
      real*8   lattice_omega
      external control_version
      external lattice_omega


      call Parallel_np(np)

      n2ft3d = 2*nfft3d
      value = BA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'electron_psi_vl_ave: out of stack memory',0, MA_ERR)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1) + 1
      n2(2) = neq(1) + neq(2)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1/scal2

!$OMP MASTER
      elocal = 0.0d0
!$OMP END MASTER

*     **** average Kohn-Sham v_local energy ****
      call Pack_c_Copy(0,dcpl_mb(vl(1)),dbl_mb(tmp1(1)))
      call Pack_c_unpack(0,dbl_mb(tmp1(1)))
      !call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
      call D3dB_cr_pfft3b(1,0,dbl_mb(tmp1(1)))
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            call D3dB_rr_Mul(1,
     >                       dbl_mb(tmp1(1)),
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)))

c           call D3dB_rc_fft3f(1,dbl_mb(tmp2(1)))
c           call Pack_c_pack(1,dbl_mb(tmp2(1)))
c           call Pack_cc_dot(1,psi1(1+(n-1)*npack1),
c    >                         dbl_mb(tmp2(1)),
c    >                         sum)
            
            call D3dB_rr_idot(1,
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)),
     >                       sum)

!$OMP MASTER
            elocal = elocal + sum*scal1*scal2
!$OMP END MASTER
         end do
      end do
!$OMP BARRIER
      if (np.gt.1) call Parallel_SumAll(elocal)
!$OMP MASTER
      if (ispin.eq.1) elocal = 2.0d0*elocal
!$OMP END MASTER

*     *** add in long range part of psp ****
      if (control_version().eq.4) then
       call D3dB_rr_dot(1,dn(1),dbl_mb(vl_lr(1)),sum)
!$OMP MASTER
       elocal = elocal + sum*dv
!$OMP END MASTER
       call D3dB_rr_dot(1,dn(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(vl_lr(1)),sum)
!$OMP MASTER
       elocal = elocal + sum*dv
!$OMP END MASTER

      end if

*     **** add in other real-space fields ****
      if (field_exist) then
       call D3dB_rr_dot(1,dn(1),dbl_mb(v_field(1)),sum)
!$OMP MASTER
       elocal = elocal + sum*dv
!$OMP END MASTER
       call D3dB_rr_dot(1,dn(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(v_field(1)),sum)
!$OMP MASTER
       elocal = elocal + sum*dv
!$OMP END MASTER
      end if
!$OMP BARRIER


*     ***** ncmp*Vl+ncmp_smooth*vlpaw terms ****
c      if (paw_exist) then
c         call nwpw_compcharge_gen_dn_cmp2(ispin,
c     >                                    dbl_mb(tmp1(1)),
c     >                                    dbl_mb(tmp2(1)))
c         call Pack_cc_dot(0,
c     >                    dbl_mb(tmp1(1)),
c     >                    dcpl_mb(vl(1)),
c     >                    sum)
c         elocal = elocal + sum
c         call Pack_cc_dot(0,
c     >                    dbl_mb(tmp2(1)),
c     >                    dcpl_mb(vlpaw(1)),
c     >                    sum)
c         elocal = elocal + sum
c
c         if ((control_version().eq.4).or.(field_exist)) then
c            call Pack_c_unpack(0,dbl_mb(tmp1(1)))
c            call D3dB_cr_pfft3b(1,0,dbl_mb(tmp1(1)))
c            call Pack_c_unpack(0,dbl_mb(tmp2(1)))
c            call D3dB_cr_pfft3b(1,0,dbl_mb(tmp2(1)))
c         end if
c         if (control_version().eq.4) then
c            call D3dB_rr_dot(1,
c     >                       dbl_mb(tmp1(1)),
c     >                       dbl_mb(vl_lr(1)),
c     >                       sum)
c            elocal = elocal + sum*dv
c            call D3dB_rr_dot(1,
c     >                       dbl_mb(tmp2(1)),
c     >                       dbl_mb(vl_lr_paw(1)),
c     >                       sum)
c            elocal = elocal + sum*dv
c         end if
c         if (field_exist) then
c            call D3dB_rr_dot(1,dbl_mb(tmp1(1)),
c     >                         dbl_mb(v_field(1)),sum)
c            elocal = elocal + sum*dv
c         end if
c      end if
 
      value = BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >           'electron_psi_vl_ave: error popping stack memory',0,
     &       MA_ERR)

      electron_psi_vl_ave = elocal
      return
      end



*     ***********************************
*     *			 		*
*     *   electron_psi_vnl_ave	 	*
*     *					*
*     ***********************************

      real*8 function electron_psi_vnl_ave(psi1,fractional,occ)
      implicit none
      complex*16 psi1(*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"

      real*8   E_vnonlocal
      external E_vnonlocal

c*     **** local variables ****
c      logical value
c      integer i,n,ms,n1(2),n2(2),np
c      integer nee(2)
c      integer n2ft3d
c      real*8 enlocal,sum
c      integer tmp1(2),tmp2(2)
c
c
c      call Parallel_np(np)
c
c      n2ft3d = 2*nfft3d
c      value = BA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
c      value = value.and.
c     >        BA_push_get(mt_dbl,(n2ft3d),'tmp2',tmp2(2),tmp2(1))
c      if (.not. value) call errquit(
c     >            'electron_psi_vl_ave: out of stack memory',0, MA_ERR)
c
c      n1(1) = 1
c      n2(1) = neq(1)
c      n1(2) = neq(1) + 1
c      n2(2) = neq(1) + neq(2)
c
c
c
c*     **** average Kohn-Sham v_nonlocal energy ****
c      nee(1) = 1
c      nee(2) = 0
c      enlocal = 0.0d0
c      do ms=1,ispin
c         do n=n1(ms),n2(ms)
c            call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
c            call v_nonlocal(ispin,nee,psi1(1+(n-1)*npack1),
c     >                      dbl_mb(tmp1(1)),
c     >                      .false.,dbl_mb(tmp2(1)),fractional,occ)
c            call Pack_cc_idot(1,psi1(1+(n-1)*npack1),
c     >                         dbl_mb(tmp1(1)),
c     >                         sum)
c            if (fractional) then
c               call Dneall_qton(n,i)
c               sum=sum*occ(i)
c            end if
c            enlocal = enlocal - sum
c         end do
c      end do
c      if (np.gt.1) call Parallel_SumAll(enlocal)
c      if (ispin.eq.1) enlocal = enlocal+enlocal
c
c 
c      value = BA_pop_stack(tmp2(2))
c      value = value.and.
c     >        BA_pop_stack(tmp1(2))
c      if (.not. value) call errquit(
c     >           'electron_psi_vl_ave: error popping stack memory',0,
c     &       MA_ERR)
c
c      electron_psi_vnl_ave = enlocal

      electron_psi_vnl_ave = E_vnonlocal(ispin,neq,fractional,occ)
      return
      end

*     ***********************************
*     *			 		*
*     *   electron_psi_v_field_ave	*
*     *					*
*     ***********************************

      real*8 function electron_psi_v_field_ave(psi1,dn)
      implicit none
      complex*16 psi1(*)
      real*8     dn(*)

#include "bafdecls.fh"
#include "electron_common.fh"


*     **** local variables ****
      integer nx,ny,nz,n2ft3d
      real*8 elocal,sum,scal1,scal2,dv

*     **** external functions ***
      real*8   lattice_omega
      external lattice_omega


      n2ft3d = 2*nfft3d

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1/scal2


      elocal = 0.0d0

*     **** add in other real-space fields ****
      if (field_exist) then
       call D3dB_rr_dot(1,dn(1),dbl_mb(v_field(1)),sum)
       elocal = elocal + sum*dv
       call D3dB_rr_dot(1,dn(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(v_field(1)),sum)
       elocal = elocal + sum*dv
      end if
      
      electron_psi_v_field_ave = elocal
      return
      end



*     ***********************************
*     *			 		*
*     *   electron_semicoreforce 	*
*     *					*
*     ***********************************

      subroutine electron_semicoreforce(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "electron_common.fh"

     
      call semicore_xc_F(ispin,dbl_mb(xcp(1)),fion)

      return
      end

*     ***********************************
*     *			 		*
*     *   electron_create_confine	*
*     *					*
*     ***********************************
      subroutine electron_create_confine()
      implicit none

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"

      !*** local variables ****
      logical value
      integer n2ft3d,r_grid(2)

      n2ft3d = 2*nfft3d
      confine_exist = .true.
      value = BA_alloc_get(mt_dbl,n2ft3d,
     >                     'v_confine',v_confine(2),v_confine(1))
      value = value.and.BA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                       r_grid(2),r_grid(1))
      if (.not.value) 
     >  call errquit('electron_create_confine: error allocate heap',0,
     >               MA_ERR)

      call Parallel_shared_vector_zero(.true.,n2ft3d,
     >                                 dbl_mb(v_confine(1)))
      call lattice_r_grid(dbl_mb(r_grid(1)))
      call electron_atom_Generate_V(1,n2ft3d,
     >                          dbl_mb(r_grid(1)),
     >                          dbl_mb(v_confine(1)))
      value = BA_pop_stack(r_grid(2))
      if (.not. value) call errquit(
     >   'electron_create_confine: error popping stack memory',0,
     >       MA_ERR)
      return
      end

*     ***********************************
*     *                                 *
*     *   electron_destroy_confine      *
*     *                                 *
*     ***********************************
      subroutine electron_destroy_confine()
      implicit none

#include "bafdecls.fh"
#include "electron_common.fh"
#include "errquit.fh"

      !*** local variables ****
      logical value

      confine_exist = .false.
      if (.not.BA_free_heap(v_confine(2)))
     >  call errquit('electron_destroy_confine: error freeing heap',0,
     >               MA_ERR)
      return
      end



*     **********************************
*     *                                *
*     *     electron_atom_Generate_V   *
*     *                                *
*     **********************************

      subroutine electron_atom_Generate_V(ctype,n2ft3d,rgrid,Vqm)
      implicit none
      integer ctype
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)

#include "stdio.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      integer ii,ia,k
      real*8  x1,y1,z1,q1,r,rs,re
      real*8  r0,dr0,beta

*     **** external functions ****
      integer  ion_katm_qm,ion_nion_qm,ion_nkatm_qm
      external ion_katm_qm,ion_nion_qm,ion_nkatm_qm
      real*8   ion_rion
      external ion_rion
      character*4 ion_atom_qm
      external    ion_atom_qm

      if (ion_nion_qm().gt.0) then
         call Parallel_taskid(taskid)

*        **** point charges ****
         if (ctype.eq.1) then
            if (taskid.eq.MASTER) then
               write(luout,*)
               write(luout,*) " - Confining virtual space"  
               write(luout,*) 
     >         " - Parameters: Symbol      Ve      rs      re"
            end if
            do ia = 1,ion_nkatm_qm()
               call control_fetch_confine(ion_atom_qm(ia),q1,rs,re)
               if (taskid.eq.MASTER) 
     >         write(luout,'(15x,A7,F8.3,F8.3,F8.3)') 
     >         ion_atom_qm(ia),q1,rs,re
            end do
            if (taskid.eq.MASTER) write(luout,*)

           do ii=1,ion_nion_qm()
              ia = ion_katm_qm(ii)
              call control_fetch_confine(ion_atom_qm(ia),q1,rs,re)
              r0  = 0.5d0*(rs+re)
              dr0 = (re-rs)
              beta = dlog(999.0d0)/dr0

              x1 = ion_rion(1,ii)
              y1 = ion_rion(2,ii)
              z1 = ion_rion(3,ii)
!$OMP DO
              do k=1,n2ft3d
                r = (rgrid(1,k)-x1)**2
     >            + (rgrid(2,k)-y1)**2
     >            + (rgrid(3,k)-z1)**2
                r = dsqrt(r)
                Vqm(k) = Vqm(k)
     >                 + q1*(1.0d0-1.0d0/(1.0d0+dexp(beta*(r-r0))))
              end do
!$OMP END DO
           end do
         end if
      end if
      return
      end





c*     ***********************************
c*     *			 		*
c*     *     electron_dn_cmp_coulomb     *
c*     *					*
c*     ***********************************
c
c      real*8 function electron_dn_cmp_coulomb()
c      implicit none
c#include "bafdecls.fh"
c#include "errquit.fh"
c#include "electron_common.fh"
c
c*     **** local variables ****
c      real*8 E
c
c      E = 0.0d0
c      if (paw_exist) then
c         call  nwpw_compcharge_gen_dn_cmp2(ispin,
c     >                                     dcpl_mb(dng_cmp(1)),
c     >                                     dcpl_mb(dng_cmp_smooth(1)))
c      end if
c
c      electron_dn_cmp_coulomb = E
c      return
c      end

