*
* $Id$
*


*     ***********************************
*     *                                 *
*     *         paw_electron_init       *
*     *                                 *
*     ***********************************
      subroutine paw_electron_init()
      implicit none

#include "bafdecls.fh"
#include "paw_electron.fh"
#include "errquit.fh"


*     **** paw_electron_counter common block ****
      integer counter
      common / paw_electron_counter / counter


*     **** local variables ****
      logical value
      integer n2ft3d

*     **** external functions ****
      logical  ion_chargeexist,ion_mmexist
      integer  paw_psi_ispin,paw_psi_ne,control_version
      external ion_chargeexist,ion_mmexist
      external paw_psi_ispin,paw_psi_ne,control_version

      counter = 0

      ispin = paw_psi_ispin()
      ne(1) = paw_psi_ne(1)
      ne(2) = paw_psi_ne(2)
      !field_exist = ion_chargeexist().or.ion_mmexist()

*     **** get nfft3d, and n2ft3d ****
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'Hpsi_k',Hpsi_k(2),Hpsi_k(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d*(ne(1)+ne(2)),
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'vl2',vl(2),vl(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'vall',vall(2),vall(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'vc',vc(2),vc(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,
     >                     'vcomp',vcomp(2),vcomp(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,
     >                     'vh',vh(2),vh(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'dn_cmp_smooth',
     >                     dn_cmp_smooth(2),
     >                     dn_cmp_smooth(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'xcp',xcp(2),xcp(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'xce',xce(2),xce(1))
      if (.not. value) 
     >  call errquit('paw_electron_init: out of heap memory',0, MA_ERR)


      return
      end


*     ***********************************
*     *                                 *
*     *         paw_electron_finalize   *
*     *                                 *
*     ***********************************
      subroutine paw_electron_finalize()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "paw_electron.fh"


*     **** local variables ****
      logical value

*     **** free heap  memory ****
      value = BA_free_heap(Hpsi_k(2))
      value = value.and.
     >        BA_free_heap(psi_r(2))
      value = value.and.
     >        BA_free_heap(vl(2))
      value = value.and.
     >        BA_free_heap(vall(2))
      value = value.and.
     >        BA_free_heap(vc(2))
      value = value.and.
     >        BA_free_heap(vcomp(2))
      value = value.and.
     >        BA_free_heap(vh(2))
      value = value.and.
     >        BA_free_heap(dn_cmp_smooth(2))
      value = value.and.
     >        BA_free_heap(xcp(2))
      value = value.and.
     >        BA_free_heap(xce(2))
      if (.not. value) 
     >  call errquit('paw_electron_init: error freeing heap memory',0,
     &       MA_ERR)

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_electron_count              *
*     *                                 *
*     ***********************************
      integer function paw_electron_count()
      implicit none

*     **** paw_electron_counter common block ****
      integer counter
      common / paw_electron_counter / counter

      paw_electron_count = counter
      return
      end


*     ***********************************
*     *                                 *
*     *         paw_electron_run                *
*     *                                 *
*     ***********************************
      subroutine paw_electron_run(psi_k,dn,dng)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)

#include "bafdecls.fh"

      logical value

*     **** paw_electron_counter common block ****
      integer counter
      common / paw_electron_counter / counter

      counter = counter+1


*     **** set overlaps ****
      call paw_ovlp_coeff_set(psi_k)
      call paw_ovlp_weights_set()

*     **** generate comp charge potential ****
      call paw_comp_charge_update()
      call paw_pot_comp_solve()


*     **** generate realspace psi ***
      call paw_electron_gen_psi_r(psi_k)

*     **** generate densities ****
      call paw_electron_gen_densities(psi_k,dn,dng)

*     **** generate hartree, compensation, and xc potentials ****
      call paw_electron_gen_scf_potentials(dn,dng)


*     **** generate Hpsi **** 
      call paw_electron_gen_Hpsi_k(psi_k)

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_electron_run_orb    *
*     *                                 *
*     ***********************************
      subroutine paw_electron_run_orb(i,psi_k)
      implicit none
      integer    i
      complex*16 psi_k(*)

      call paw_electron_gen_psi_r_orb(i,psi_k)
c      call paw_electron_gen_Hpsi_k_orb(i,psi_k)

      return
      end

*     ***********************************
*     *                                 *
*     *     paw_electron_sd_update      *
*     *                                 *
*     ***********************************
      subroutine paw_electron_sd_update(psi1,psi2,dte)
      implicit none
      complex*16 psi1(*),psi2(*)
      real*8     dte

#include "bafdecls.fh"
#include "paw_electron.fh"



      call dcopy(2*(npack1)*(ne(1)+ne(2)),psi1,1,psi2,1)
      call daxpy(2*(npack1)*(ne(1)+ne(2)),
     >           (-dte),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           psi2,1)

      return
      end


 
*     ***********************************
*     *                                 *
*     *         paw_electron_energy     *
*     *                                 *
*     ***********************************
      real*8 function paw_electron_energy(psi_k,dn,dng)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "paw_energy_kin_atom.fh"
#include "paw_energy_vloc_atom.fh"
#include "paw_energy_ion_atom.fh"
#include "paw_energy_core_atom.fh"
#include "paw_energy_hartree_atom.fh"
#include "paw_xc.fh"
#include "paw_electron.fh"


*     **** local variables ****
      integer n2ft3d
      integer ii,ms,n1(2),n2(2),nx,ny,nz
      real*8  sum,ehartr,exc,exc2,dv,eke,evloc
      real*8  total_energy,occ(1)

*     **** external functions *****
      integer  control_version
      real*8   lattice_omega,paw_electron_ehartree
      real*8   paw_mult_energy_atom_comp
      real*8   paw_mult_energy_atom_self
      real*8   paw_mult_energy_atom_mult

      external control_version
      external lattice_omega,paw_electron_ehartree
      external paw_mult_energy_atom_comp
      external paw_mult_energy_atom_self
      external paw_mult_energy_atom_mult
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
     

      dv = lattice_omega()/dble(nx*ny*nz)

      n2ft3d = 2*nfft3d


*     **** get kinetic energies ****
      call ke_ave(ispin,ne,psi_k,eke,.false.,occ)
      eke = eke + paw_energy_kin_atom()
c      write(*,*) "eke:",eke

      
*     **** get coulomb energy ****
      ehartr = paw_electron_ehartree(dng)
     >       + paw_energy_hartree_atom()
     >       + paw_energy_ion_atom()
     >       + paw_energy_core_atom()
     >       + paw_energy_core_ion_atom()
     >       + paw_energy_core_self_atom()
     >       + paw_mult_energy_atom_self()
     >       + paw_mult_energy_atom_comp()
     >       + paw_mult_energy_atom_mult()
c      write(*,*) "ehartr:",ehartr
c      write(*,*) "    pw:",paw_electron_ehartree(dng)
c      write(*,*) "  atom:",  paw_energy_hartree_atom()
c     >       + paw_energy_ion_atom()
c     >       + paw_energy_core_atom()
c     >       + paw_energy_core_ion_atom()
c     >       + paw_energy_core_self_atom()
c     >       + paw_mult_energy_atom_self()
c     >       + paw_mult_energy_atom_comp()
c      write(*,*) "  mult:",paw_mult_energy_atom_mult()



*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dn(1),dbl_mb(xce(1)),exc)
      if (ispin.eq.1) then
         exc= exc + exc
      else
         call D3dB_rr_dot(1,dn(1+n2ft3d),dbl_mb(xce(1)),exc2)
         exc= exc + exc2
      end if
      exc = exc*dv + paw_energy_xc_atom()      
c      write(*,*) "exc:",exc

*     **** get local potential energy ****
      call Pack_cc_dot(0,dng,dcpl_mb(vl(1)),evloc)
      evloc = evloc + paw_energy_vloc_atom()
c      write(*,*) "evloc:",evloc


          
      total_energy = eke + ehartr + exc + evloc

c      write(*,*) "total_energy:",total_energy

*     **** total energy ****
      paw_electron_energy = total_energy

      return
      end




*     ***********************************
*     *                                 *
*     *         paw_electron_eorbit     *
*     *                                 *
*     ***********************************
      real*8 function paw_electron_eorbit(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      integer ii,ms,n1(2),n2(2)
      real*8  sum,eorbit

      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)


*     *** get orbital energies ****
      eorbit = 0.0d0
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
           call Pack_cc_idot(1,psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       sum)

           eorbit = eorbit + sum

         end do
      end do
      call D3dB_SumAll(eorbit)
      if (ispin.eq.1) eorbit = eorbit+eorbit

           
      paw_electron_eorbit = eorbit 

      return
      end
 

*     ***********************************
*     *                                 *
*     *      paw_electron_ehartree      *
*     *                                 *
*     ***********************************
      real*8 function paw_electron_ehartree(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "paw_electron.fh"
#include "errquit.fh"

      !**** local varialbles ****
      logical value
      integer tmp1(2),tmp2(2)
      real*8 ehartr
      

      !**** external functions ****
      real*8   lattice_omega
      external lattice_omega

     
      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >'paw_electron_ehartree: out of stack memory',0,MA_ERR)
 



*     **** get coulomb energy ****
      call Pack_cc_Sum(0,
     >                   dng,
     >                   dbl_mb(dn_cmp_smooth(1)),
     >                   dcpl_mb(tmp1(1)))
      call Pack_c_Copy(0,
     >                   dcpl_mb(vcomp(1)),
     >                   dcpl_mb(tmp2(1)))
      call Pack_cc_daxpy(0,0.5d0,
     >                   dcpl_mb(vh(1)),
     >                   dcpl_mb(tmp2(1)))
      call Pack_cc_dot(0,
     >                  dcpl_mb(tmp1(1)),
     >                  dcpl_mb(tmp2(1)),
     >                  ehartr)

c      call Pack_cc_dot(0,
c     >                  dbl_mb(dn_cmp_smooth(1)),
c     >                  dcpl_mb(tmp2(1)),
c     >                  ehartr)

c      call Pack_cc_dot(0,
c     >                  dng,
c     >                  dcpl_mb(vcomp(1)),
c     >                  ehartr)
      !ehartr = 0.5d0*ehartr

      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >'paw_electron_ehartree: error popping stack memory',0,
     > MA_ERR)


      paw_electron_ehartree = ehartr*lattice_omega()

      return
      end



*     ***********************************
*     *                                 *
*     *         paw_electron_exc        *
*     *                                 *
*     ***********************************
      real*8 function paw_electron_exc(dn)
      implicit none
      real*8 dn(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      integer nx,ny,nz
      real*8  exc,exc2,dv

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dn,
     >                 dbl_mb(xce(1)),
     >                 exc)
      if (ispin.eq.1) then
         exc= exc + exc 
      else
         call D3dB_rr_dot(1,dn(1+2*nfft3d),
     >                    dbl_mb(xce(1)),
     >                    exc2)
         exc= exc + exc2
      end if
      exc = exc*dv
        
  
      paw_electron_exc =  exc 

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_electron_pxc        *
*     *                                 *
*     ***********************************
      real*8 function paw_electron_pxc(dn)
      implicit none
      real*8 dn(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


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
      call D3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     >                 pxc)
      if (ispin.eq.1) then
         pxc= pxc + pxc 
      else
         call D3dB_rr_dot(1,dn(1+2*nfft3d),
     >                    dbl_mb(xcp(1)+2*nfft3d),
     >                    pxc2)
         pxc= pxc + pxc2
      end if
      pxc = pxc*dv
        
  
      paw_electron_pxc =  pxc

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_electron_pxc_rho    *
*     *                                 *
*     ***********************************
      real*8 function paw_electron_pxc_rho(rho)
      implicit none
      real*8 rho(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


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
          
      paw_electron_pxc_rho =  pxc
      return
      end





*     ***********************************
*     *                                 *
*     *         paw_electron_get_Hpsi_k *
*     *                                 *
*     ***********************************
      subroutine paw_electron_get_Hpsi_k(Hpsi_k_new)
      implicit none
      complex*16 Hpsi_k_new(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


      call dcopy(2*npack1*(ne(1)+ne(2)),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           Hpsi_k_new,1)
      return
      end



*     ***************************
*     *                         *
*     *    paw_electron_ispin   *
*     *                         *
*     ***************************
      integer function paw_electron_ispin()
      implicit none

#include "paw_electron.fh"

      paw_electron_ispin = ispin
      return
      end


*     ***************************
*     *                         *
*     *      paw_electron_ne    *
*     *                         *
*     ***************************
      integer function paw_electron_ne(ms)
      implicit none
      integer ms

#include "paw_electron.fh"

      paw_electron_ne = ne(ms)
      return
      end


*     ***********************************
*     *                                 *
*     *     paw_electron_get_Tgradient  *
*     *                                 *
*     ***********************************

      subroutine paw_electron_get_Tgradient(psi_k,hml,THpsi_k)

      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     ***** local variables ****
      integer ms,n,shift1,shift2
     


*     **** no preconditioning ***
      do ms=1,ispin
            n     = ne(ms)
            if (n.le.0) go to 30
            shift1 = 1 + (ms-1)*ne(1)*npack1
            shift2 = 1 + (ms-1)*ne(1)*ne(1)

            call dgemm('N','N',2*npack1,n,n,
     >                (1.0d0),
     >                psi_k(shift1),  2*npack1,
     >                hml(shift2),    n,
     >                (0.0d0),
     >                THpsi_k(shift1),2*npack1)
  30        continue
      end do
      call daxpy(2*npack1*(ne(1)+ne(2)),
     >           (-1.0d0),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           THpsi_k,1)
      return
      end


*     ***********************************
*     *                                 *
*     *     paw_electron_gen_Tangent    *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_Tangent(psi_k,hml,THpsi_k)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

*     ***** local variables ****
      integer ms,n,shift1,shift2
     
      do ms=1,ispin
         n     = ne(ms)
         if (n.le.0) go to 30
         shift1 = 1 + (ms-1)*ne(1)*npack1
         shift2 = 1 + (ms-1)*ne(1)*ne(1)
         call dgemm('N','N',2*npack1,n,n,
     >             (-1.0d0),
     >             psi_k(shift1),  2*npack1,
     >             hml(shift2),    n,
     >             (1.0d0),
     >             THpsi_k(shift1),2*npack1)
   30    continue
      end do
      return
      end


*     ***********************************
*     *                                 *
*     *     paw_electron_get_Gradient   *
*     *                                 *
*     ***********************************

      subroutine paw_electron_get_Gradient(THpsi_k)
      implicit none
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


      call dcopy(2*npack1*(ne(1)+ne(2)),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           THpsi_k,1)      
      return
      end


*     *******************************************
*     *                                         *
*     *     paw_electron_get_gradient_orb       *
*     *                                         *
*     *******************************************

      subroutine paw_electron_get_gradient_orb(i,Horb)
      implicit none
      integer i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

      call Pack_c_Copy(1,dcpl_mb(Hpsi_k(1)+(i-1)*npack1),Horb) 
      
      return
      end



*     *******************************************
*     *                                         *
*     *     paw_electron_get_TMgradient         *
*     *                                         *
*     *******************************************

      subroutine paw_electron_get_TMgradient(psi_k,THpsi_k)
      implicit none
#include "errquit.fh"
      complex*16 psi_k(*)
      complex*16 THpsi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     ***** local variables ****
      logical value
      integer ms,n,n1(2),shift,mhml(2)
     

      n1(1) = 1
      n1(2) = ne(1)+1

      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'mhml',mhml(2),mhml(1))
      if (.not. value)
     >call errquit('paw_electron_get_Tgradient: push stack',0, MA_ERR)



*     **** generate M*H|psi> *****
      call Grsm_gg_Copy(npack1,(ne(1)+ne(2)),
     >                  dcpl_mb(Hpsi_k(1)),
     >                  THpsi_k)
      call ke_precondition(npack1,(ne(1)+ne(2)),
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
      call dscal(2*ne(1)*ne(1),(-1.0d0),dbl_mb(mhml(1)),1)


*     **** generate TMG = M*H|psi> - |psi>*mhml ****
      do ms=1,ispin
            shift = (ms-1)*ne(1)*ne(1)
            n     = ne(ms)
            call Grsm_gmg_daxpy(npack1,n,
     >                        psi_k(1+(ms-1)*ne(1)*npack1),
     >                        dbl_mb(mhml(1)+shift),
     >                        THpsi_k(1+(ms-1)*ne(1)*npack1))
      end do

      call Grsm_gg_dScale(npack1,(ne(1)+ne(2)),
     >                    (-1.0d0),
     >                    THpsi_k,
     >                    THpsi_k)


      value = BA_pop_stack(mhml(2))
      if (.not. value)
     >call errquit('paw_electron_get_Tradient: popping stack',0,MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *     paw_electron_gen_hml        *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_hml(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

*     **** local variables ****
      integer ms,n,shift

      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         n     = ne(ms)
         if (n.le.0) go to 30
         call Pack_ccm_sym_dot(1,n,
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
     >                     hml(shift+1))
  30     continue
      end do
      return
      end

*     ***********************************
*     *                                 *
*     *     paw_electron_gen_hml_S      *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_hml_S(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

*     **** local variables ****
      integer ms,n,shift

      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         n     = ne(ms)
         if (n.le.0) go to 30
         call paw_overlap_matrix_gen(n,n,
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
     >                     hml(shift+1))
  30     continue
      end do
      return
      end




*     ***********************************
*     *                                 *
*     *     paw_electron_gen_hml_g      *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_hml_g(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      integer ms,n,n1(2),shift


      n1(1) = 1
      n1(2) = ne(1) + 1


      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         n     = ne(ms)
         call Pack_ccm_dot(1,n,
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
     >                     hml(shift+1))
      end do


      return
      end


*     ***********************************
*     *                                 *
*     *     paw_electron_gen_hmlt       *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_hmlt(psi_k,hmlt)
      implicit none
      complex*16 psi_k(*)
      real*8     hmlt(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      integer ms,n,n1(2),shift
     

      n1(1) = 1
      n1(2) = ne(1) + 1


      do ms=1,ispin
         n     = ne(ms)
         if (n.le.0) go to 30
         shift = (ms-1)*ne(1)*ne(1)
         call Pack_ccm_dot(1,n,
     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     hmlt(shift+1))
   30    continue
      end do


      return
      end



*     **********************************************
*     *                                            *
*     *     paw_electron_gen_psiTangenthml_S       *
*     *                                            *
*     **********************************************

      subroutine paw_elecpsitang(psi_k,THpsi_k,hml)
      implicit none
      complex*16 psi_k(*)
      complex*16 THpsi_k(*)
      real*8     hml(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      integer ms,n,shift

      do ms=1,ispin
         n     = ne(ms)
         if (n.le.0) go to 30
         shift = (ms-1)*ne(1)*ne(1)
         call paw_overlap_matrix_gen(n,n,
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     THpsi_k(1+(ms-1)*ne(1)*npack1),
     >                     hml(shift+1))
   30    continue
      end do

      return
      end






**************************************************************************
**************************************************************************
*******    routines below this line are for internal use only    *********
**************************************************************************
**************************************************************************

*     *******************************************
*     *                                         *
*     *         paw_electron_gen_Hpsi_k         *
*     *                                         *
*     *******************************************

      subroutine paw_electron_gen_Hpsi_k(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      logical move
      integer n
      real*8  fion(3,1)

*     **** external functions ****
      integer  control_version
      external control_version

      move = .false.

*     ******************
*     **** get Hpsi ****
*     ******************
        call paw_psi_H(ispin,ne,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dcpl_mb(vl(1)),
     >             dcpl_mb(vc(1)),
     >             dbl_mb(xcp(1)),
     >             dcpl_mb(Hpsi_k(1)),
     >             move,
     >             fion)


      call Grsm_gg_dscale(npack1,(ne(1)+ne(2)),(-1.0d0),
     >                    dcpl_mb(Hpsi_k(1)),
     >                    dcpl_mb(Hpsi_k(1)))

      return
      end



c*     ***********************************
c*     *                                        *
c*     *     paw_electron_gen_Hpsi_k_orb        *
c*     *                                        *
c*     ***********************************
c
c      subroutine paw_electron_gen_Hpsi_k_orb(n,psi_k)
c      implicit none
c      integer n
c      complex*16 psi_k(*)
c
c#include "bafdecls.fh"
c#include "paw_electron.fh"
c
c
c*     **** local variables ****
c      integer nne(2)
c
c*     **** external functions ****
c      integer  control_version
c      external control_version
c
c
c    
c      if (n.le.ne(1)) then
c        nne(1) = 1
c        nne(2) = 0
c      else
c        nne(1) = 0
c        nne(2) = 1
c      end if 
c
c*     ******************
c*     **** get Hpsi ****
c*     ******************
c      call paw_psi_Horb(ispin,nne,psi_k(1+(n-1)*npack1),
c     >             dbl_mb(psi_r(1)+(n-1)*2*nfft3d),
c     >             dbl_mb(vall(1)),
c     >             dcpl_mb(Hpsi_k(1)+(n-1)*npack1))
c
c
c      call Pack_c_SMul(1,(-1.0d0),
c     >                 dcpl_mb(Hpsi_k(1)+(n-1)*npack1),
c     >                 dcpl_mb(Hpsi_k(1)+(n-1)*npack1))
c
c      return
c      end

c*     ***********************************
c*     *                                        *
c*     *   paw_electron_get_gradient_virtual *
c*     *                                        *
c*     ***********************************
c
c      subroutine paw_electron_get_gradient_virtual(ms,orb,Horb)
c      implicit none
c#include "errquit.fh"
c      integer    ms
c      complex*16 orb(*)
c      complex*16 Horb(*)
c
c#include "bafdecls.fh"
c#include "paw_electron.fh"
c
c
c*     **** local variables ****
c      logical value
c      integer nne(2),n2ft3d
c      integer tmp_r(2)
c
c*     **** external functions ****
c      integer  control_version
c      external control_version
c
c      if (ms.eq.1) then
c        nne(1) = 1
c        nne(2) = 0
c      else
c        nne(1) = 0
c        nne(2) = 1
c      end if 
c
c      n2ft3d = 2*nfft3d
c
c      value = BA_push_get(mt_dbl,(n2ft3d),'tmp_r',tmp_r(2),tmp_r(1))
c      if (.not. value)
c     >   call errquit('paw_electron_get_gradient_virtual: push stack',0,
c     &       MA_ERR)
c
c
c      call Pack_c_Copy(1,orb,dbl_mb(tmp_r(1)))
c      call Pack_c_unpack(1,  dbl_mb(tmp_r(1)))
c      call D3dB_cr_fft3b(1,  dbl_mb(tmp_r(1))) 
c
c
c*     **** get Hpsi ****
c      if (control_version().eq.3)
c     >  call psi_Horb(ispin,nne,orb,
c     >             dbl_mb(tmp_r(1)),
c     >             dbl_mb(vall(1)),
c     >             Horb)
c
c      if (control_version().eq.4)
c     >  call psi_Horbv4(ispin,nne,orb,
c     >             dbl_mb(tmp_r(1)),
c     >             dbl_mb(vall(1)),
c     >             Horb)
c
c      call Pack_c_SMul(1,(-1.0d0),Horb,Horb)
c
c      value = BA_pop_stack(tmp_r(2))
c      if (.not. value) call errquit(
c     >     'paw_electron_get_gradient_virtual: poping stack',1, MA_ERR)
c
c
c      return
c      end



*     ***********************************
*     *                                 *
*     *     paw_electron_gen_psi_r      *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_psi_r(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

*     **** local variables ****
      integer n,nemax,n2ft3d


      n2ft3d = 2*nfft3d
      nemax = ne(1) + ne(2)

c     call Grsm_gg_Copy(npack1,nemax,psi_k,dbl_mb(psi_r(1)))
      do n=1,nemax
         call Pack_c_Copy(1,psi_k(1+(n-1)*npack1),
     >                    dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      end do

      call Grsm_gh_fftb(nfft3d,nemax,dbl_mb(psi_r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemax,dbl_mb(psi_r(1)))

      return
      end

*     ***********************************
*     *                                 *
*     *     paw_electron_gen_psi_r_orb  *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_psi_r_orb(n,psi_k)
      implicit none
      integer    n
      complex*16 psi_k(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

*     **** local variables ****
      integer n2ft3d

      n2ft3d = 2*nfft3d

      call Pack_c_Copy(1,psi_k(1+(n-1)*npack1),
     >                        dbl_mb(psi_r(1)+(n-1)*n2ft3d))
     
      call Pack_c_unpack(1,   dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      call D3dB_cr_fft3b(1,   dbl_mb(psi_r(1)+(n-1)*n2ft3d)) 
      call D3dB_r_Zero_Ends(1,dbl_mb(psi_r(1)+(n-1)*n2ft3d))

      return
      end


*     ***********************************
*     *                                 *
*     *    paw_electron_gen_densities   *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_densities(psi_k,dn,dng)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "paw_electron.fh"

     
*     **** local variables ****
      integer i
      integer ms,n2ft3d
      integer n,n1(2),n2(2)
      real*8  scal2
      integer tmp1(2)
      logical value

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega


      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)

      n2ft3d = 2*nfft3d
      scal2 = 1.0d0/lattice_omega()


*     *********************
*     **** generate dn ****
*     *********************
      call dcopy(2*n2ft3d,0.0d0,0,dn,1)

      do ms=1,ispin
         do n=n1(ms),n2(ms)
            do i=1,n2ft3d
               dn(i+(ms-1)*n2ft3d) 
     >            = dn(i+(ms-1)*n2ft3d) 
     >            + scal2*(dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d)**2)
            end do
         end do
         call D3dB_r_Zero_Ends(1,dn(1+(ms-1)*n2ft3d))
      end do

*     **** generate dng  ****
      call paw_electron_gen_dng(dn,dng)

      return
      end



*     ***********************************
*     *                                 *
*     *      paw_electron_gen_dng       *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_dng(dn,dng)
      implicit none
      real*8     dn(*)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "paw_electron.fh"
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
     >     'paw_electron_gen_dng: out of stack memory',0, MA_ERR)

      call D3dB_rr_Sum(1,dn,dn(1+(ispin-1)*n2ft3d),dbl_mb(tmp1(1)))
c      call D3dB_r_SMul(1,scal1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
      call D3dB_r_SMul1(1,scal1,dbl_mb(tmp1(1)))
      call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
      call Pack_c_pack(0,dbl_mb(tmp1(1)))
      call Pack_c_Copy(0,dbl_mb(tmp1(1)),dng)


      value = BA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >     'paw_electron_gen_dng: poping stack',1, MA_ERR)


      return
      end


*     ***********************************
*     *                                 *
*     *        paw_electron_gen_vall    *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_vall()
      implicit none
   

#include "bafdecls.fh"
#include "paw_electron.fh"
     
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
        call Pack_cc_Sum2(0,dcpl_mb(vc(1)),
     >                     dbl_mb(vall(1)))
      
*       **** fourier transform k-space potentials ****
        call Pack_c_unpack(0,dbl_mb(vall(1)))
        call D3dB_cr_fft3b(1,dbl_mb(vall(1)))

*       **** add v_field to vall ****
c        if (field_exist)
c     >    call D3dB_rr_Sum(1,dbl_mb(vall(1)),
c     >                       dbl_mb(v_field(1)),
c     >                       dbl_mb(vall(1)))
        if (field_exist)
     >    call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),
     >                        dbl_mb(vall(1)))

      else

*       **** add up k-space potentials, vall = scal2*vsr_l    ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))
      
*        **** fourier transform k-space potentials ****
         call Pack_c_unpack(0,dbl_mb(vall(1)))
         call D3dB_cr_fft3b(1,dbl_mb(vall(1)))

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
         call D3dB_rr_Sum2(1,dbl_mb(vl_lr(1)),
     >                       dbl_mb(vall(1)))
         call D3dB_rr_Sum2(1,dcpl_mb(vc(1)),   
     >                       dbl_mb(vall(1)))
         if (field_exist)
     >     call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),
     >                         dbl_mb(vall(1)))

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
     >                   dbl_mb(xcp(1)+2*nfft3d),
     >                   dbl_mb(vall(1)+2*nfft3d))
        call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)+2*nfft3d))
      end if
      call D3dB_rr_Sum2(1,dbl_mb(xcp(1)),dbl_mb(vall(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)))

      return
      end



*     ***********************************
*     *                                 *
*     *        paw_electron_get_vall    *
*     *                                 *
*     ***********************************

      subroutine paw_electron_get_vall(vall_out)
      implicit none
      real*8 vall_out

#include "bafdecls.fh"
#include "paw_electron.fh"

      call dcopy(4*nfft3d,dbl_mb(vall(1)),1,vall_out,1)
      return
      end


*     ***********************************
*     *                                 *
*     *        paw_electron_set_vall    *
*     *                                 *
*     ***********************************

      subroutine paw_electron_set_vall(vall_in)
      implicit none
      real*8 vall_in

#include "bafdecls.fh"
#include "paw_electron.fh"

      call dcopy(4*nfft3d,vall_in,1,dbl_mb(vall(1)),1)
      return
      end


*     *******************************************
*     *                                         *
*     *   paw_electron_gen_scf_potentials       *
*     *                                         *
*     *******************************************

      subroutine paw_electron_gen_scf_potentials(dn,dng)
      implicit none
      real*8     dn(*)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "paw_electron.fh"
#include "nwpwxc.fh"
#include "util.fh"
#include "errquit.fh"


*     ***** local variables ****
      logical value
      logical use_lda
      integer n2ft3d,nx,ny,nz,gga
      integer tmp1(2),tmp2(2)
      real*8  vzero,scal1

*     **** external functions ****
      integer  control_gga,control_version
      external control_gga,control_version

     
      n2ft3d = 2*nfft3d
      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >'paw_electron_gen_scf_potentials: out of stack memory',0,MA_ERR)
 

*     ************************************
*     **** generate coulomb potential ****
*     ************************************
      call paw_pot_hartree_solve()

*     *** atomic portion ***

      call paw_mult_dn_cmp_get(dcpl_mb(tmp1(1)),
     >                         dbl_mb(dn_cmp_smooth(1)))

      if (control_version().eq.3)  then
           call Pack_cc_Sub(0,
     >                      dcpl_mb(tmp1(1)),
     >                      dbl_mb(dn_cmp_smooth(1)),
     >                      dcpl_mb(tmp2(1))) !** tmp2 = dn_cmp - dn_cmp_smooth
           call Pack_cc_Sum(0,
     >                      dbl_mb(dn_cmp_smooth(1)),
     >                      dng,
     >                      dcpl_mb(tmp1(1))) !** tmp1 = dng+dn_cmp_smooth **

           !**** vh *****
           call coulomb_v(dcpl_mb(tmp1(1)),
     >                    dcpl_mb(vh(1)))


           !**** vcmp *****
           call coulomb_v(dcpl_mb(tmp2(1)),
     >                    dcpl_mb(vcomp(1)))
           call paw_mult_vzero(vzero)
           call Pack_c_setzero(0,vzero,dcpl_mb(vcomp(1)))

           call paw_mult_coeff_set(dcpl_mb(vh(1)),dcpl_mb(vcomp(1)))

           call Pack_cc_Sum(0,
     >                      dcpl_mb(vh(1)),
     >                      dcpl_mb(vcomp(1)),
     >                      dcpl_mb(vc(1)))

      end if

      if (control_version().eq.4) then

           call D3dB_nx(1,nx)
           call D3dB_ny(1,ny)
           call D3dB_nz(1,nz)
           scal1 = 1.0d0/dble(nx*ny*nz)
           call Pack_cc_Sub(0,
     >                      dcpl_mb(tmp1(1)),
     >                      dbl_mb(dn_cmp_smooth(1)),
     >                      dcpl_mb(tmp2(1))) !** tmp2 = dn_cmp - dn_cmp_smooth
c           call Pack_c_Copy(0,
c     >                      dcpl_mb(tmp1(1)),
c     >                      dcpl_mb(tmp2(1))) !** tmp2 = dn_cmp
           call Pack_cc_Sum(0,
     >                      dbl_mb(dn_cmp_smooth(1)),
     >                      dng,
     >                      dcpl_mb(tmp1(1))) !** tmp1 = dng+dn_cmp_smooth **


           call Pack_c_unpack(0,dcpl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dcpl_mb(tmp1(1)))

           call coulomb2_v(dcpl_mb(tmp1(1)),dcpl_mb(vh(1)))

           call D3dB_rc_fft3f(1,dcpl_mb(vh(1)))
           call D3dB_r_SMul1(1,scal1,dcpl_mb(vh(1)))
           call Pack_c_pack(0,dcpl_mb(vh(1)))


           call Pack_c_unpack(0,dcpl_mb(tmp2(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp2(1)))
           call D3dB_r_Zero_Ends(1,dcpl_mb(tmp2(1)))

           call coulomb2_v(dcpl_mb(tmp2(1)),dcpl_mb(vcomp(1)))

           call D3dB_rc_fft3f(1,dcpl_mb(vcomp(1)))
           call D3dB_r_SMul1(1,scal1,dcpl_mb(vcomp(1)))
           call Pack_c_pack(0,dcpl_mb(vcomp(1)))

           call paw_mult_vzero(vzero)
           call Pack_c_setzero(0,vzero,dcpl_mb(vcomp(1)))

           call paw_mult_coeff_set(dcpl_mb(vh(1)),dcpl_mb(vcomp(1)))

           call Pack_cc_Sum(0,
     >                      dcpl_mb(vh(1)),
     >                      dcpl_mb(vcomp(1)),
     >                      dcpl_mb(vc(1)))
      end if


*     *************************************************
*     **** generate exchange-correlation potential ****
*     *************************************************

*     **** local portion ****
c      call paw_density_solve()
      call paw_xc_solve()
      

*     **** plane wave portion ****
      gga = control_gga()
      
      use_lda = ((.not.nwpwxc_is_on().and.gga.eq.0).or.
     &             (nwpwxc_is_on().and.nwpwxc_is_lda()))
      
      if (use_lda) then
        call vxc(n2ft3d,ispin,dn,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)),
     >                    dcpl_mb(tmp1(1)))

      else
        call v_bwexc(gga,n2ft3d,ispin,dn,
     >                    1.0d0,1.0d0,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)))
      end if

      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >'paw_electron_gen_scf_potentials: error popping stack memory',0,
     > MA_ERR)


      return
      end




*     ***********************************
*     *                                 *
*     *   paw_electron_gen_vl_potential *
*     *                                 *
*     ***********************************

      subroutine paw_electron_gen_vl_potential()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "paw_electron.fh"


*     **** local variables ****
      logical move,value
      integer n2ft3d
      integer tmp1(2)
      integer tmp2(2)

*     **** external functions *****
      integer  control_version
      external control_version

      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        BA_push_get(mt_dbl,(3),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call 
     >   errquit('paw_electron_gen_vl_potential: out of stack memory',0,
     &       MA_ERR)

      move = .false.
      call paw_vloc(dcpl_mb(vl(1)),
     >              move,
     >              dcpl_mb(tmp1(1)),
     >              dbl_mb(tmp2(1)))

      value = BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >  'paw_electron_gen_vl_potential: error popping stack memory',0,
     &       MA_ERR)

      return
      end



*     ***********************************
*     *                                 *
*     *   paw_electron_psi_vl_ave       *
*     *                                 *
*     ***********************************

      real*8 function paw_electron_psi_vl_ave(psi1,dn)
      implicit none
      complex*16 psi1(*)
      real*8     dn(*)

#include "bafdecls.fh"
#include "paw_electron.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2)
      integer nx,ny,nz,n2ft3d,np
      real*8 elocal,sum,scal1,scal2,dv
      integer tmp1(2),tmp2(2)

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
     >'paw_electron_psi_vl_ave: out of stack memory',0, MA_ERR)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1/scal2


*     **** average Kohn-Sham v_local energy ****
       call Pack_c_Copy(0,dcpl_mb(vl(1)),dbl_mb(tmp1(1)))
      call Pack_c_unpack(0,dbl_mb(tmp1(1)))
      call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
      elocal = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            call D3dB_rr_Mul(1,
     >                       dbl_mb(tmp1(1)),
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)))

            
            call D3dB_rr_idot(1,
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)),
     >                       sum)

            elocal = elocal + sum*scal1*scal2
         end do
      end do
      if (np.gt.1) call D3dB_SumAll(elocal)
      if (ispin.eq.1) elocal = 2.0d0*elocal

 
      value = BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >'paw_electron_psi_vl_ave: error popping stack memory',0,MA_ERR)

      paw_electron_psi_vl_ave = elocal
      return
      end

*     ***********************************
*     *                                 *
*     *       paw_electron_force        *
*     *                                 *
*     ***********************************

      subroutine paw_electron_force(psi1,dng1,fion)
      implicit none
      complex*16 psi1(*)
      complex*16 dng1(*)
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "paw_electron.fh"

      !**** local variables ****
      logical value
      integer tmp(2),nemax


      nemax = ne(1)+ne(2)
      call paw_vloc_force(dng1,fion)
      call paw_mult_pw_force(dcpl_mb(vh(1)),
     >                        dcpl_mb(vcomp(1)),
     >                        fion)
      call paw_force_solve(psi1,fion)

      value = BA_push_get(mt_dbl,(8*nemax*nemax),
     >                     'tmp',tmp(2),tmp(1))
      if (.not. value) 
     > call errquit('paw_electron_force:out of stack memory',0,0)


      call paw_electron_gen_hml(psi1,dbl_mb(tmp(1)))

      call paw_force_constraint(dbl_mb(tmp(1)),fion)

      value = BA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('paw_electron_force:popping stack memory',0,0)



       return
       end


