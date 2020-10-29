*
* $Id$
*


*     **********************************
*     *                                *
*     *      c_electron_iptr_psir      *
*     *                                *
*     **********************************
      integer function c_electron_iptr_psir()
      implicit none

#include "bafdecls.fh"
#include "c_electron_common.fh"

      integer iptr
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      iptr = cpsi_data_get_chnk(psi_r_tag,1)
      c_electron_iptr_psir = iptr
      return
      end


*     ***********************************
*     *					*
*     *		c_electron_init		*
*     *					*
*     ***********************************
      subroutine c_electron_init()
      implicit none

#include "bafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer neall

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne
      integer  cpsi_data_alloc,Pneb_ispinq,Pneb_nbrillq
      logical  cpsi_spin_orbit,control_hfxon_virtual
      external cpsi_ispin,cpsi_ne
      external cpsi_data_alloc,Pneb_ispinq,Pneb_nbrillq
      external cpsi_spin_orbit,control_hfxon_virtual

      counter = 0

      ispin  = cpsi_ispin()
      ne(1)  = cpsi_ne(1)
      ne(2)  = cpsi_ne(2)
      spin_orbit= cpsi_spin_orbit()
      hfxon_virtual = control_hfxon_virtual()

      ispinq  = Pneb_ispinq()
      call Pneb_neq(neq)
      nbrillq = Pneb_nbrillq()

*     **** get nfft3d, and npack0 and nbrill****
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      call C3dB_nfft3d(1,nfft3d)

      neall  = neq(1)+neq(2)
      
*     **** allocate memory ****
c      value = BA_alloc_get(mt_dcpl,nbrill*npack1*(ne(1)+ne(2)),
c     >                     'Hpsi_k',Hpsi_k(2),Hpsi_k(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,nbrill*nfft3d*(ne(1)+ne(2)),
c     >                     'psi_r',psi_r(2),psi_r(1))

      Hpsi_k_tag = cpsi_data_alloc(nbrillq,neall,2*npack1)
      psi_r_tag  = cpsi_data_alloc(nbrillq,neall,2*nfft3d)

      value = BA_alloc_get(mt_dcpl,npack0,'vl2',vl(2),vl(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,'vc',vc(2),vc(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,2*nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nfft3d,'xcp',xcp(2),xcp(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nfft3d,'xce',xce(2),xce(1))
      if (.not. value) 
     >  call errquit('c_electron_init: out of heap memory',0, MA_ERR)

*     *** initializer symmetrizer ***
      call c_rho_symmetrizer_init()

      return
      end


*     ***********************************
*     *					*
*     *		c_electron_finalize     *
*     *					*
*     ***********************************
      subroutine c_electron_finalize()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value

*     **** free heap  memory ****
c      value = BA_free_heap(Hpsi_k(2))
c      value = value.and.
c     >        BA_free_heap(psi_r(2))


*     **** free heap  memory ****
      call c_rho_symmetrizer_end()

      call cpsi_data_dealloc(Hpsi_k_tag)
      call cpsi_data_dealloc(psi_r_tag)

      value =           BA_free_heap(vl(2))
      value = value.and.BA_free_heap(vall(2))
      value = value.and.BA_free_heap(vc(2))
      value = value.and.BA_free_heap(xcp(2))
      value = value.and.BA_free_heap(xce(2))
      if (.not. value) 
     >  call errquit(
     >   'c_electron_finalize: error freeing heap memory',0, MA_ERR)
      return
      end

*     ***********************************
*     *					*
*     *		c_electron_cout		*
*     *					*
*     ***********************************
      integer function c_electron_count()
      implicit none

#include "c_electron_common.fh"

      c_electron_count = counter
      return
      end

*     ***********************************
*     *					*
*     *		c_electron_run_noscf	*
*     *					*
*     ***********************************
      subroutine c_electron_run_noscf(psi_k_tag,dn,dng,dnall)
      implicit none
c      complex*16 psi_k(*)
      integer    psi_k_tag
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "c_electron_common.fh"

      counter = counter+1

      call c_electron_gen_psi_r(psi_k_tag)
c      call c_electron_gen_densities(dn,dng,dnall)
      call c_electron_gen_Hpsi_k_vall(psi_k_tag)
      return
      end

*     ***********************************
*     *					*
*     *		c_electron_run		*
*     *					*
*     ***********************************
      subroutine c_electron_run(psi_k_tag,dn,dng,dnall)
      implicit none
c      complex*16 psi_k(*)
      integer    psi_k_tag
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "c_electron_common.fh"

!$OMP MASTER
      counter = counter+1
!$OMP END MASTER

      call c_electron_gen_psi_r(psi_k_tag)
      call c_electron_gen_densities(dn,dng,dnall)
      call c_electron_gen_scf_potentials(dn,dng,dnall)
      call c_electron_gen_Hpsi_k(psi_k_tag)

      return
      end

*     ***********************************
*     *                                 *
*     *         c_electron_run_orb      *
*     *                                 *
*     ***********************************
      subroutine c_electron_run_orb(nb,i,psi_k_tag)
      implicit none
      integer    nb,i
c      complex*16 psi_k(*)
      integer    psi_k_tag

      call c_electron_gen_psi_r_orb(nb, i,psi_k_tag)
      call c_electron_gen_Hpsi_k_orb(nb,i,psi_k_tag)

      return
      end


*     ***********************************
*     *                                 *
*     *    c_electron_get_gradient_orb  *
*     *                                 *
*     ***********************************

      subroutine c_electron_get_gradient_orb(nb,i,Horb)
      implicit none
      integer nb,i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer Hpsi_shift

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      Hpsi_shift = cpsi_data_get_ptr(Hpsi_k_tag,nb,i)
      call Cram_c_Copy(nb,dbl_mb(Hpsi_shift),Horb)

      if (spin_orbit) then
        Hpsi_shift = cpsi_data_get_ptr(Hpsi_k_tag,nb,i+neq(1))
        call Cram_c_Copy(nb,dbl_mb(Hpsi_shift),Horb(npack1+1))
      end if
      return
      end

*     **********************************************
*     *                                            *
*     *         c_electron_get_gradient_virtual    *
*     *                                            *
*     **********************************************

      subroutine c_electron_get_gradient_virtual(nb,ms,orb,Horb)
      implicit none
      integer    nb,ms
      complex*16 orb(*),Horb(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"

*     **** local variables ****
      integer tmp_r(2)

      if (.not.BA_push_get(mt_dcpl,(nfft3d),'tmp_r',tmp_r(2),tmp_r(1)))
     >   call errquit('c_electron_get_gradient_virtual: push stack',0,
     &       MA_ERR)

      call Cram_c_Copy(nb,orb,dcpl_mb(tmp_r(1)))
      call Cram_c_unpack(nb,dcpl_mb(tmp_r(1)))
      call C3dB_cr_pfft3b(1,nb,dcpl_mb(tmp_r(1)))

*     ******************
*     **** get Hpsi ****
*     ******************
      call cpsi_Horb(hfxon_virtual,nb,ispin,ms,psi_r_tag,
     >              dcpl_mb(vall(1)),
     >              orb,
     >              dcpl_mb(tmp_r(1)),
     >              Horb)
      call Cram_c_SMul1(nb,(-1.0d0),Horb)

c      if (spin_orbit) then
c       call Cram_c_SMul1(nb,(-1.0d0),dbl_mb(hpsi_shift))
c      end if


      if (.not.BA_pop_stack(tmp_r(2)))
     >   call errquit('c_electron_get_gradient_virtual: poping stack',1,
     >                MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		c_electron_sd_update  	*
*     *					*
*     ***********************************
*     ***********************************
      subroutine c_electron_sd_update(psi1_tag,psi2_tag,dte)
      implicit none
      integer    psi1_tag,psi2_tag
      real*8     dte

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer nb,psi1_shift,psi2_shift,Hpsi_shift

*     ***** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call cpsi_data_update(psi2_tag)
      do nb=1,nbrillq
         psi1_shift = cpsi_data_get_chnk(  psi1_tag,nb)
         psi2_shift = cpsi_data_get_chnk(  psi2_tag,nb)
         Hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nb)
         call c_electron_sd_subupdate(nb,npack1,(neq(1)+neq(2)),
     >                           dbl_mb(psi1_shift),
     >                           dbl_mb(psi2_shift),
     >                           dbl_mb(Hpsi_shift),
     >                           dte)
      end do
      call cpsi_data_noupdate(psi2_tag)
      return
      end
cccccccccccc
      subroutine c_electron_sd_subupdate(nb,nfft3d,nn,
     >              psi1,psi2,Hpsi,dte)
      implicit none
      integer    nb,nfft3d,nn
      complex*16 psi1(nfft3d,nn)
      complex*16 psi2(nfft3d,nn)
      complex*16 Hpsi(nfft3d,nn)
      real*8     dte

      integer n
*     ************************************
*     **** do a steepest descent step ****
*     ************************************
      do n=1,nn
        call Cram_c_SMul(nb,(-dte),Hpsi(1,n),psi2(1,n))
c        call Cram_cc_Sum(nb,psi2(1,n),psi1(1,n),psi2(1,n))
        call Cram_cc_Sum2(nb,psi1(1,n),psi2(1,n))
      end do

      return
      end

cccccccccccccccccc


*     ***********************************
*     *					*
*     *		c_electron_energy	*
*     *					*
*     ***********************************
      real*8 function c_electron_energy(psi_k_tag,dn,dng,dnall)
      implicit none
      integer    psi_k_tag
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      integer nbq,ii,nx,ny,nz,neall
      integer psi_shift,Hpsi_shift,nshift1,nshift2,occ_shift,occ_tag
      real*8  sum,eorbit,ehartr,exc,pxc,exc2,pxc2,dv,weight
      real*8  total_energy

*     **** external functions *****
      logical  band_HFX,band_HFX_relaxed
      integer  cpsi_data_get_chnk,cpsi_data_nsize,Pneb_convert_nb
      integer  cpsi_data_get_next
      real*8   lattice_omega,c_coulomb_e,brillioun_weight
      external band_HFX,band_HFX_relaxed
      external cpsi_data_get_chnk,cpsi_data_nsize,Pneb_convert_nb
      external cpsi_data_get_next
      external lattice_omega,c_coulomb_e,brillioun_weight

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)

      neall = neq(1) + neq(2)
      nshift1 = cpsi_data_nsize(psi_k_tag)
      nshift2 = cpsi_data_nsize(Hpsi_k_tag)
      occ_tag = cpsi_data_get_next(psi_k_tag)

*     *** get orbital energies ****
      eorbit = 0.0d0
      do nbq=1,nbrillq
        !weight = brillioun_weight(Pneb_convert_nb(nbq))
        weight = brillioun_weight(nbq)
        psi_shift  = cpsi_data_get_chnk( psi_k_tag,nbq)
        Hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nbq)

        if (occ_tag.gt.0) then
          occ_shift = cpsi_data_get_chnk(occ_tag,nbq)
          do ii=1,neall
           call Cram_cc_idot(nbq,
     >                       dbl_mb(psi_shift),
     >                       dbl_mb(Hpsi_shift),
     >                       sum)
           eorbit = eorbit + sum*weight*dbl_mb(occ_shift)
           psi_shift  = psi_shift  + nshift1
           Hpsi_shift = Hpsi_shift + nshift2
           occ_shift  = occ_shift + 1
          end do
        else
          do ii=1,neall
           call Cram_cc_idot(nbq,
     >                       dbl_mb(psi_shift),
     >                       dbl_mb(Hpsi_shift),
     >                       sum)
           eorbit = eorbit + sum*weight
           psi_shift  = psi_shift  + nshift1
           Hpsi_shift = Hpsi_shift + nshift2
          end do
        end if

      end do
      call Parallel_SumAll(eorbit)
      if (ispin.eq.1) eorbit = eorbit+eorbit

           
*     **** get coulomb energy ****
      ehartr = c_coulomb_e(dng)
              
*     **** get exchange-correlation energy ****
      call C3dB_rr_dot(1,dnall(1),
     >                dbl_mb(xce(1)),
     >                exc)
      call C3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     .                 pxc)
      if (ispin.eq.1) then
         exc= exc + exc 
         pxc= pxc + pxc 
      else
         call C3dB_rr_dot(1,dnall(1+nfft3d),
     >                    dbl_mb(xce(1)),
     >                    exc2)
         call C3dB_rr_dot(1,dn(1+nfft3d),
     >                    dbl_mb(xcp(1)+nfft3d),
     >                    pxc2)
         exc= exc + exc2
         pxc= pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv

      total_energy = eorbit + exc - ehartr - pxc


*     **** HFX energy ****
      if (band_HFX()) then
      if (band_HFX_relaxed()) then
         call band_energy_HFX(ispin,psi_r_tag,
     >                        ehfx,
     >                        phfx)
         total_energy = total_energy + ehfx - phfx
      end if
      end if

*     **** total energy ****
      c_electron_energy = total_energy
      return
      end


*     ***********************************
*     *					*
*     *		c_electron_eorbit	*
*     *					*
*     ***********************************
      real*8 function c_electron_eorbit(psi_k_tag)
      implicit none
      integer psi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer nbq,ii,neall
      integer psi_shift,Hpsi_shift,nshift1,nshift2,occ_shift,occ_tag
      real*8  sum,eorbit,weight

*     **** external functions *****
      integer  cpsi_data_get_chnk,cpsi_data_nsize,Pneb_convert_nb
      integer  cpsi_data_get_next
      real*8   brillioun_weight
      external cpsi_data_get_chnk,cpsi_data_nsize,Pneb_convert_nb
      external cpsi_data_get_next
      external brillioun_weight

      neall = neq(1) + neq(2)
      nshift1 = cpsi_data_nsize(psi_k_tag)
      nshift2 = cpsi_data_nsize(Hpsi_k_tag)
      occ_tag = cpsi_data_get_next(psi_k_tag)

*     *** get orbital energies ****
      eorbit = 0.0d0
      do nbq=1,nbrillq
        !weight = brillioun_weight(Pneb_convert_nb(nbq))
        weight = brillioun_weight(nbq)
        psi_shift  = cpsi_data_get_chnk( psi_k_tag,nbq)
        Hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nbq)
        if (occ_tag.gt.0) then
          occ_shift = cpsi_data_get_chnk(occ_tag,nbq)
          do ii=1,neall
           call Cram_cc_idot(nbq,
     >                       dbl_mb(psi_shift),
     >                       dbl_mb(Hpsi_shift),
     >                       sum)
           eorbit = eorbit + sum*weight*dbl_mb(occ_shift)
           psi_shift  = psi_shift  + nshift1
           Hpsi_shift = Hpsi_shift + nshift2
           occ_shift  = occ_shift + 1
          end do
        else
          do ii=1,neall
           call Cram_cc_idot(nbq,
     >                       dbl_mb(psi_shift),
     >                       dbl_mb(Hpsi_shift),
     >                       sum)
           eorbit = eorbit + sum*weight
           psi_shift  = psi_shift  + nshift1
           Hpsi_shift = Hpsi_shift + nshift2
          end do
        end if

      end do
      call Parallel_SumAll(eorbit)
      if (ispin.eq.1) eorbit = eorbit+eorbit

      c_electron_eorbit = eorbit 
      return
      end
 

*     ***********************************
*     *					*
*     *		c_electron_ehartree	*
*     *					*
*     ***********************************
      real*8 function c_electron_ehartree(dng)
      implicit none
      complex*16 dng(*)


*     **** external functions ****
      real*8   c_coulomb_e
      external c_coulomb_e

      c_electron_ehartree = c_coulomb_e(dng)
      return
      end


*     ***********************************
*     *					*
*     *		c_electron_exc		*
*     *					*
*     ***********************************
      real*8 function c_electron_exc(dnall)
      implicit none
      real*8 dnall(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      integer nx,ny,nz
      real*8  exc,exc2,dv

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)

*     **** get exchange-correlation energy ****
      call C3dB_rr_dot(1,dnall,
     >                 dbl_mb(xce(1)),
     >                 exc)
      if (ispin.eq.1) then
         exc= exc + exc 
      else
         call C3dB_rr_dot(1,dnall(1+nfft3d),
     >                    dbl_mb(xce(1)),
     >                    exc2)
         exc= exc + exc2
      end if
      exc = exc*dv
        
      c_electron_exc =  exc 
      return
      end


*     ***********************************
*     *					*
*     *		c_electron_pxc		*
*     *					*
*     ***********************************
      real*8 function c_electron_pxc(dn)
      implicit none
      real*8 dn(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer nx,ny,nz
      real*8  pxc,pxc2,dv

*     **** external functions *****
      real*8   lattice_omega
      external lattice_omega

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** get exchange-correlation energy ****
      call C3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     >                 pxc)
      if (ispin.eq.1) then
         pxc= pxc + pxc 
      else
         call C3dB_rr_dot(1,dn(1+nfft3d),
     >                    dbl_mb(xcp(1)+nfft3d),
     >                    pxc2)
         pxc= pxc + pxc2
      end if
      pxc = pxc*dv
       
      c_electron_pxc =  pxc
      return
      end


 

*     ***********************************
*     *                                 *
*     *         c_electron_pxc_rho      *
*     *                                 *
*     ***********************************
      real*8 function c_electron_pxc_rho(rho)
      implicit none
      real*8 rho(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer nx,ny,nz
      real*8  pxc,pxc2,dv

*     **** external functions *****
      real*8   lattice_omega
      external lattice_omega

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** get exchange-correlation energy ****
      call C3dB_rr_dot(1,rho,
     >                 dbl_mb(xcp(1)),
     >                 pxc)
      if (ispin.eq.1) then
         pxc= pxc + pxc
      else
         call C3dB_rr_dot(1,rho,
     >                    dbl_mb(xcp(1)+nfft3d),
     >                    pxc2)
         pxc= pxc + pxc2
      end if
      pxc = pxc*dv

      c_electron_pxc_rho =  pxc
      return
      end


*     ***********************************
*     *                                 *
*     *     c_electron_HFX_energies      *
*     *                                 *
*     ***********************************

      subroutine c_electron_HFX_energies(ehfx0,phfx0)
      implicit none
      real*8 ehfx0,phfx0

#include "bafdecls.fh"
#include "c_electron_common.fh"

      logical  band_HFX_relaxed
      external band_HFX_relaxed
      
      if (.not.band_HFX_relaxed()) then
         call band_energy_HFX(ispin,psi_r_tag,
     >                        ehfx,
     >                        phfx)
         phfx = 0.0d0
      end if

      ehfx0 = ehfx
      phfx0 = phfx
      return
      end




*     ***********************************
*     *					*
*     *		c_electron_get_Hpsi_k	*
*     *					*
*     ***********************************
      subroutine c_electron_get_Hpsi_k(Hpsi_k_new_tag)
      implicit none
      integer  Hpsi_k_new_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

      call cpsi_data_copyall(Hpsi_k_tag,Hpsi_k_new_tag)
      return
      end



*     ***************************
*     *				*
*     *	   c_electron_ispin	*
*     *				*
*     ***************************
      integer function c_electron_ispin()
      implicit none

#include "c_electron_common.fh"

      c_electron_ispin = ispin
      return
      end


*     ***************************
*     *				*
*     *	     c_electron_ne	*
*     *				*
*     ***************************
      integer function c_electron_ne(ms)
      implicit none
      integer ms

#include "c_electron_common.fh"

      c_electron_ne = ne(ms)
      return
      end


*     ***********************************
*     *					*
*     *	    c_electron_get_Tgradient 	*
*     *					*
*     ***********************************

      subroutine c_electron_get_Tgradient(psi_k_tag,hml_tag,THpsi_k_tag)
      implicit none
      integer    psi_k_tag
      integer    hml_tag
      integer    THpsi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      integer nb,nsize,psi_shift,hpsi_shift,tpsi_shift,hml_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      nsize = 2*(neq(1)+neq(2))*npack1

      call cpsi_data_update(THpsi_k_tag)
      do nb=1,nbrillq
         psi_shift  = cpsi_data_get_chnk(psi_k_tag,nb)
         hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nb)
         tpsi_shift = cpsi_data_get_chnk(THpsi_k_tag,nb)
         hml_shift  = cpsi_data_get_chnk(hml_tag,nb)

         call Pneb_fwf_Multiply(0,nb,
     >                       one,
     >                       dbl_mb(psi_shift),npack1,
     >                       dbl_mb(hml_shift),
     >                       zero,
     >                       dbl_mb(tpsi_shift))
         call c_electron_ff_sub(nsize,
     >                          dbl_mb(hpsi_shift),
     >                          dbl_mb(tpsi_shift))
c         call BGrsm_fff_Sub(dcpl_mb(Hpsi_k(1)),
c     >                   THpsi_k,
c     >                   THpsi_k)
c
c         call BGrsm_ff_dScale(-1.0d0,
c     >                     THpsi_k,
c     >                     THpsi_k)

      end do
      call cpsi_data_noupdate(THpsi_k_tag)
      return
      end

      subroutine c_electron_ff_sub(nsize,A,B)
      implicit none
      integer nsize
      real*8 A(*)
      real*8 B(*)
      integer i
      do i=1,nsize
         B(i) = B(i) - A(i)
      end do
      return
      end

*     ***************************
*     *				*
*     *	    c_electron_gen_hml 	*
*     *				*
*     ***************************
      subroutine c_electron_gen_hml(psi_k_tag,hml_tag)
      implicit none
      integer psi_k_tag
      integer hml_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer nb,psi_shift,hpsi_shift,hml_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk
 
      call cpsi_data_update(hml_tag)
      do nb=1,nbrillq
         psi_shift  = cpsi_data_get_chnk(psi_k_tag,nb)
         hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nb)
         hml_shift  = cpsi_data_get_chnk(hml_tag,nb)
         call Pneb_ffw_hermit_Multiply(0,nb,
     >                  dbl_mb(psi_shift),
     >                  dbl_mb(hpsi_shift),npack1,
     >                  dbl_mb(hml_shift))
      end do
      call cpsi_data_noupdate(hml_tag)

      return
      end

**************************************************************************
**************************************************************************
*******    routines below this line are for internal use only    *********
**************************************************************************
**************************************************************************

*     ***********************************
*     *					*
*     *		c_electron_gen_Hpsi_k	*
*     *					*
*     ***********************************

      subroutine c_electron_gen_Hpsi_k(psi_k_tag)
      implicit none
      integer  psi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical move
      integer nb,hpsi_shift,nsize
      real*8  fion(3,1)

*     **** external functions ****
      integer  control_version,cpsi_data_get_chnk
      external control_version,cpsi_data_get_chnk

      move = .false.
*     ******************
*     **** get Hpsi ****
*     ******************
      call cpsi_H(ispin,neq,
     >           psi_k_tag,
     >           psi_r_tag,
     >           dcpl_mb(vl(1)),
     >           dcpl_mb(vc(1)),
     >           dbl_mb(xcp(1)),
     >           Hpsi_k_tag,
     >           move,
     >           fion)

      nsize = 2*(neq(1)+neq(2))*npack1
      do nb=1,nbrillq
        hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nb)
        call dscal(nsize,(-1.0d0),dbl_mb(hpsi_shift),1)
      end do
      return
      end

*     ***********************************
*     *                                 *
*     *   c_electron_gen_Hpsi_k_vall    *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_Hpsi_k_vall(psi_k_tag)
      implicit none
      integer psi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer nb,hpsi_shift,nsize

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

*     ******************
*     **** get Hpsi ****
*     ******************
      call cpsi_H_vall(ispin,ne,
     >           psi_k_tag,
     >           psi_r_tag,
     >           dcpl_mb(vall(1)),
     >           Hpsi_k_tag)

      nsize = 2*(neq(1)+neq(2))*npack1
      do nb=1,nbrillq
        hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nb)
        call dscal(nsize,(-1.0d0),dbl_mb(hpsi_shift),1)
      end do
      return
      end


*     ***********************************
*     *                                 *
*     *     c_electron_gen_Hpsi_k_orb   *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_Hpsi_k_orb(nb,n,psi_k_tag)
      implicit none
      integer nb,n
      integer psi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      integer ms
      integer psir_shift,psik_shift,hpsi_shift

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      if (n.le.ne(1)) then
        ms=1
      else
        ms=2
      end if
      psir_shift = cpsi_data_get_ptr( psi_r_tag,nb,n)
      psik_shift = cpsi_data_get_ptr( psi_k_tag,nb,n)
      hpsi_shift = cpsi_data_get_ptr(Hpsi_k_tag,nb,n)

      call cpsi_data_update(Hpsi_k_tag)
*     ******************
*     **** get Hpsi ****
*     ******************
      call cpsi_Horb(.true.,nb,ispin,ms,
     >              psi_r_tag,
     >              dcpl_mb(vall(1)),
     >              dbl_mb(psik_shift),
     >              dbl_mb(psir_shift),
     >              dbl_mb(hpsi_shift))

c      call Cram_c_SMul(nb,(-1.0d0),
c     >                 dbl_mb(hpsi_shift),
c     >                 dbl_mb(hpsi_shift))
      call Cram_c_SMul1(nb,(-1.0d0),dbl_mb(hpsi_shift))

      if (spin_orbit) then
         hpsi_shift = cpsi_data_get_ptr(Hpsi_k_tag,nb,n+neq(1))
c         call Cram_c_SMul(nb,(-1.0d0),
c     >                 dbl_mb(hpsi_shift),
c     >                 dbl_mb(hpsi_shift))
         call Cram_c_SMul1(nb,(-1.0d0),dbl_mb(hpsi_shift))
      end if               
      call cpsi_data_noupdate(Hpsi_k_tag)
      return
      end






*     *******************************
*     *				    *
*     *	    c_electron_gen_psi_r    *
*     *				    *
*     *******************************

      subroutine c_electron_gen_psi_r(psi_k_tag)
      implicit none
      integer psi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

c     **** local variables ****
      logical done
      integer k,neall,nb,nbq
      integer ki,kj,nbi,nbj
      integer shiftA,nshiftA
      integer shiftB,nshiftB

c     **** external functions ***
      logical  C3dB_cr_pfft3_queue_filled
      integer  Pneb_convert_nb,cpsi_data_get_chnk,cpsi_data_get_next
      external C3dB_cr_pfft3_queue_filled
      external Pneb_convert_nb,cpsi_data_get_chnk,cpsi_data_get_next

c      call BGrsm_fh_fftb(psi_k,
c     >                   dcpl_mb(psi_r(1)))

      neall   = neq(1)+neq(2)
      nshiftA = 2*npack1
      nshiftB = 2*nfft3d

c      do nbq=1,nbrillq
c         nb = Pneb_convert_nb(nbq)
c         shiftA = cpsi_data_get_chnk(psi_k_tag,nbq)
c         shiftB = cpsi_data_get_chnk(psi_r_tag,nbq)
c         do k=1,neall
c            call Cram_c_Copy(nb,dbl_mb(shiftA),dbl_mb(shiftB))
c            call Cram_c_unpack(nb,dbl_mb(shiftB))
c            !call C3dB_cr_fft3b(1,dbl_mb(shiftB))
c            call C3dB_cr_pfft3b(1,nb,dbl_mb(shiftB))
c            shiftA = shiftA + nshiftA
c            shiftB = shiftB + nshiftB
c         end do
c      end do


      nbi = 1
      ki  = 1 
      shiftA =  cpsi_data_get_chnk(psi_k_tag,nbi)

      nbj = 1
      kj  = 1 
      shiftB =  cpsi_data_get_chnk(psi_r_tag,nbj)

c      done = .false.
      done = (nbrillq.eq.0)
      do while (.not.done)

         if (nbi.le.nbrillq) then
            call C3dB_cr_pfft3b_queuein(nbi,dbl_mb(shiftA))
            shiftA = shiftA + nshiftA
            ki = ki+1
            if (ki.gt.neall) then
               nbi = nbi + 1
               ki  = 1
               if (nbi.le.nbrillq) then
                  shiftA =  cpsi_data_get_chnk(psi_k_tag,nbi)
               endif
            end if
         end if

         if ((C3dB_cr_pfft3_queue_filled()).or.(nbi.gt.nbrillq)) then
            call C3dB_cr_pfft3b_queueout(nbj,dbl_mb(shiftB))
            shiftB = shiftB + nshiftB
            kj = kj+1
            if (kj.gt.neall) then
               nbj = nbj + 1
               kj  = 1
               if (nbj.le.nbrillq) then
                  shiftB =  cpsi_data_get_chnk(psi_r_tag,nbj)
               endif
            end if
         endif
         done = ((nbi.gt.nbrillq).and.(nbj.gt.nbrillq))
      end do

*     **** set the occupations ****
      call cpsi_data_set_next(psi_r_tag,cpsi_data_get_next(psi_k_tag))

      return
      end


*     ***********************************
*     *                                 *
*     *     c_electron_gen_psi_r_orb    *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_psi_r_orb(nb,n,psi_k_tag)
      implicit none
      integer    nb,n
      integer    psi_k_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer shiftk,shiftr

*     ***** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      call cpsi_data_update(psi_r_tag)
      shiftk = cpsi_data_get_ptr(psi_k_tag,nb,n)
      shiftr = cpsi_data_get_ptr(psi_r_tag,nb,n)

      call Cram_c_Copy(nb,dbl_mb(shiftk),dbl_mb(shiftr))

      call Cram_c_unpack(nb, dbl_mb(shiftr))
      !call C3dB_cr_fft3b(1,dbl_mb(shiftr))
      call C3dB_cr_pfft3b(1,nb,dbl_mb(shiftr))

      if (spin_orbit) then
        shiftk = cpsi_data_get_ptr(psi_k_tag,nb,n+neq(1))
        shiftr = cpsi_data_get_ptr(psi_r_tag,nb,n+neq(1))
        call Cram_c_Copy(nb,dbl_mb(shiftk),dbl_mb(shiftr))
        call Cram_c_unpack(nb, dbl_mb(shiftr))
        !call C3dB_cr_fft3b(1,dbl_mb(shiftr))
        call C3dB_cr_pfft3b(1,nb,dbl_mb(shiftr))
      end if                                   
      call cpsi_data_noupdate(psi_r_tag)
      return
      end


*     ***********************************
*     *				        *
*     *	    c_electron_gen_densities	*
*     *				        *
*     ***********************************

      subroutine c_electron_gen_densities(dn,dng,dnall)
      implicit none
#include "bafdecls.fh"
#include "errquit.fh"
#include "c_electron_common.fh"

      real*8     dn(nfft3d,*)
      complex*16 dng(*)
      real*8     dnall(nfft3d,*)

     
*     **** local variables ****
      integer ms,nx,ny,nz,neall,n,n1(2),n2(2)
      integer shiftA,nshiftA,shifto,nbq,nb,occ_tag
      real*8  scal1,scal2,weight
      integer tmp1(2)
      logical value

*     ***** external functions *****
      logical  cpsp_semicore
      integer  cpsi_data_get_chnk,Pneb_convert_nb,cpsi_data_get_next
      real*8   lattice_omega,brillioun_weight,cpsi_data_get_ptr
      external cpsp_semicore
      external cpsi_data_get_chnk,Pneb_convert_nb,cpsi_data_get_next
      external lattice_omega,brillioun_weight,cpsi_data_get_ptr

      neall = neq(1)+neq(2)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()

*     *********************
*     **** generate dn ****
*     *********************
c      call C3dB_r_Zero(1,dn)
c      if (ispin.gt.1) call C3dB_r_Zero(1,dn(1,ispin))
      call dcopy(ispin*nfft3d,0.0d0,0,dn,1)

c      call BGrsm_hr_SumSqr(dcpl_mb(psi_r(1)),
c     >                     dn)
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1)+1
      n2(2) = neq(1)+neq(2)
      nshiftA = 2*nfft3d

      occ_tag = cpsi_data_get_next(psi_r_tag)
      if (occ_tag.gt.0) then
       do nbq=1,nbrillq
        weight = brillioun_weight(nbq)
        shiftA = cpsi_data_get_chnk(psi_r_tag,nbq)
        shifto = cpsi_data_get_chnk(occ_tag,nbq)
        do ms=1,ispin
        do n=n1(ms),n2(ms)
           call C3dB_cr_aSqrpy(1,dbl_mb(shifto)*weight,
     >                           dbl_mb(shiftA),
     >                           dn(1,ms))
           shiftA = shiftA + nshiftA
           shifto = shifto + 1
        end do
        end do
       end do
      else
       do nbq=1,nbrillq
        weight = brillioun_weight(nbq)
        shiftA = cpsi_data_get_chnk(psi_r_tag,nbq)
        do ms=1,ispin
        do n=n1(ms),n2(ms)
           call C3dB_cr_aSqrpy(1,weight,dbl_mb(shiftA),dn(1,ms))
           shiftA = shiftA + nshiftA
        end do
        end do
       end do
      end if


*     *** reduce over brillioun zone ***
      do ms=1,ispin
       call K1dB_Vector_SumAll(nfft3d,dn(1,ms))
      end do

      call C3dB_r_SMul1(1,scal2,dn)
      if (ispin.gt.1) 
     >  call C3dB_r_SMul1(1,scal2,dn(1,ispin))

      call c_rho_symmetrizer_dn(ispin,nfft3d,dn)


        
*     **********************
*     **** generate dng ****
*     **********************
      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_densities: out of stack memory',0, MA_ERR)

      call C3dB_rrc_Sum(1,dn,dn(1,ispin),dcpl_mb(tmp1(1)))
      call C3dB_rc_pfft3f(1,0,dcpl_mb(tmp1(1)))
      call C3dB_c_SMul1(1,scal1,dcpl_mb(tmp1(1)))
      call Cram_c_pack(0,dcpl_mb(tmp1(1)))
      call Cram_c_Copy(0,dcpl_mb(tmp1(1)),dng)


*     ********************************************************
*     **** generate dnall - used for semicore corrections ****
*     ********************************************************
      if (cpsp_semicore(0)) then
         call c_semicore_density(dcpl_mb(tmp1(1)))
c         call C3dB_r_SMul(1,0.5d0,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
         call C3dB_r_SMul1(1,0.5d0,dcpl_mb(tmp1(1)))
      else
         call C3dB_r_Zero(1,dcpl_mb(tmp1(1)))
      end if
      do ms=1,ispin
        call C3dB_rr_Sum(1,dn(1,ms),
     >                     dcpl_mb(tmp1(1)),
     >                     dnall(1,ms))
      end do

      value = BA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >     'c_electron_gen_densities: error popping stack memory',0,
     &       MA_ERR)
      return
      end


*     ***********************************
*     *                                 *
*     *     c_electron_gen_dng_dnall    *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_dng_dnall(dn,dng,dnall)
      implicit none
#include "bafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"

      real*8     dn(nfft3d,*)
      complex*16 dng(*)
      real*8     dnall(nfft3d,*)



*     **** local variables ****
      integer ms,nx,ny,nz,neall
      real*8  scal1
      integer tmp1(2)
      logical value

*     ***** external functions *****
      logical  cpsp_semicore
      real*8   lattice_omega
      external cpsp_semicore
      external lattice_omega

      neall = ne(1)+ne(2)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)


*     **********************
*     **** generate dng ****
*     **********************
      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
         if (.not. value) call errquit(
     >     'electron_gen_dng_dnall: out of stack memory',0, MA_ERR)

      call C3dB_rrc_Sum(1,dn,dn(1,ispin),dcpl_mb(tmp1(1)))
c      call C3dB_c_SMul(1,scal1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
      !call C3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
      call C3dB_rc_pfft3f(1,0,dcpl_mb(tmp1(1)))
      call C3dB_c_SMul1(1,scal1,dcpl_mb(tmp1(1)))
      call Cram_c_pack(0,dcpl_mb(tmp1(1)))
      call Cram_c_Copy(0,dcpl_mb(tmp1(1)),dng)


*     ********************************************************
*     **** generate dnall - used for semicore corrections ****
*     ********************************************************
      if (cpsp_semicore(0)) then
         call c_semicore_density(dcpl_mb(tmp1(1)))
         call C3dB_r_SMul1(1,0.5d0,dcpl_mb(tmp1(1)))
      else
         call C3dB_r_Zero(1,dcpl_mb(tmp1(1)))
      end if
      do ms=1,ispin
        call C3dB_rr_Sum(1,dn(1,ms),
     >                     dcpl_mb(tmp1(1)),
     >                     dnall(1,ms))
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >     'c_electron_gen_dng_dnall: error popping stack memory',0,
     &       MA_ERR)
      return
      end



*     ***********************************
*     *                                 *
*     *        c_electron_gen_vall0     *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_vall0()
      implicit none


#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer ms
      real*8 scal2

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      scal2 = 1.0d0/lattice_omega()

*     **** add up k-space potentials, vall = scal2*vl + vc  ****
      call Cram_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                         dcpl_mb(vall(1)))


*     **** fourier transform k-space potentials ****
      call Cram_c_unpack(0,dcpl_mb(vall(1)))
      !call C3dB_cr_fft3b(1,dcpl_mb(vall(1)))
      call C3dB_cr_pfft3b(1,0,dcpl_mb(vall(1)))


      if (ispin.eq.2)
     >   call C3dB_c_Copy(1,dcpl_mb(vall(1)),
     >                      dcpl_mb(vall(1)+nfft3d))

      return
      end




*     ***********************************
*     *                                 *
*     *        c_electron_gen_vall      *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_vall()
      implicit none


#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer ms
      real*8 scal2

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      scal2 = 1.0d0/lattice_omega()

*     **** add up k-space potentials, vall = scal2*vl + vc  ****
      call Cram_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                         dcpl_mb(vall(1)))
c      call Cram_cc_Sum(0,dcpl_mb(vall(1)),
c     >                   dcpl_mb(vc(1)),
c     >                   dcpl_mb(vall(1)))
      call Cram_cc_Sum2(0,dcpl_mb(vc(1)),
     >                    dcpl_mb(vall(1)))

*     **** fourier transform k-space potentials ****
      call Cram_c_unpack(0,dcpl_mb(vall(1)))
      !call C3dB_cr_fft3b(1,dcpl_mb(vall(1)))
      call C3dB_cr_pfft3b(1,0,dcpl_mb(vall(1)))

*     **** add xcp to vall ****
c      do ms=ispin,1,-1
c        call C3dB_rc_Sum(1,
c     >                   dbl_mb(xcp(1)+(ms-1)*nfft3d),
c     >                   dcpl_mb(vall(1)),
c     >                   dcpl_mb(vall(1)+(ms-1)*nfft3d))
c      end do
      if (ispin.eq.2)
     >   call C3dB_rc_Sum(1,dbl_mb(xcp(1)+nfft3d),
     >                      dcpl_mb(vall(1)),
     >                      dcpl_mb(vall(1)+nfft3d))
      call C3dB_rc_Sum2(1,dbl_mb(xcp(1)),dcpl_mb(vall(1)))

      return
      end


*     ***********************************
*     *                                 *
*     *        c_electron_get_vall      *
*     *                                 *
*     ***********************************

      subroutine c_electron_get_vall(vall_out)
      implicit none
      complex*16 vall_out(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

      call dcopy(4*nfft3d,dcpl_mb(vall(1)),1,vall_out,1)
      return
      end


*     ***********************************
*     *                                 *
*     *        c_electron_set_vall        *
*     *                                 *
*     ***********************************

      subroutine c_electron_set_vall(vall_in)
      implicit none
      complex*16 vall_in(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

      call dcopy(4*nfft3d,vall_in,1,dcpl_mb(vall(1)),1)
      return
      end





*     ***********************************
*     *			 		*
*     *   c_electron_gen_scf_potentials	*
*     *					*
*     ***********************************

      subroutine c_electron_gen_scf_potentials(dn,dng,dnall)
      implicit none
#include "bafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"

      real*8     dn(nfft3d,*)
      complex*16 dng(*)
      real*8     dnall(nfft3d,*)



*     ***** local variables ****
      logical value
      integer tmp1(2)
      integer trho(2),txcp(2),txce(2),n2ft3d
      integer ms,i,j,k,q,nx,ny,nz,nq,indx1,indx2,gga

*     **** external functions ****
      integer  control_gga,control_version
      external control_gga,control_version
      logical  nwpwxc_is_on,control_only_lda
      external nwpwxc_is_on,control_only_lda
     

      call c_coulomb_v(dng,dcpl_mb(vc(1)))

      gga = control_gga()

      !*** LDA calculation ****
      !if ((.not.nwpwxc_is_on()).and.(gga.eq.0)) then
      if (control_only_lda()) then
         value = BA_push_get(mt_dbl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',0, MA_ERR)
        call vxc(nfft3d,ispin,dnall,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)),
     >                    dbl_mb(tmp1(1)))
        value = BA_pop_stack(tmp1(2)) 
        if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: error popping stack memory',0,
     >       MA_ERR)


      !*** this code needs to be optimized!!! ****
      !*** GGA calculation ****
      !else if ((gga.ge.10).and.(gga.lt.200)).or. then
      else 
         call D3dB_n2ft3d(1,n2ft3d)
         value = BA_push_get(mt_dbl,(2*n2ft3d),'trho',trho(2),trho(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(2*n2ft3d),'txcp',txcp(2),txcp(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(2*n2ft3d),'txce',txce(2),txce(1))
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',1, MA_ERR)

         do ms=1,ispin
           call C3dB_D3dB_r_Copy(1,dnall(1,ms),
     >                             dbl_mb(trho(1)+(ms-1)*n2ft3d))
         end do

         call v_bwexc_all(gga,n2ft3d,ispin,
     >                    dbl_mb(trho(1)),
     >                    dbl_mb(txcp(1)),
     >                    dbl_mb(txce(1)))

         do ms=1,ispin
           call D3dB_C3dB_r_Copy(1,dbl_mb(txcp(1)+(ms-1)*n2ft3d),
     >                             dbl_mb(xcp(1) +(ms-1)*nfft3d))
         end do

         call D3dB_C3dB_r_Copy(1,dbl_mb(txce(1)),
     >                           dbl_mb(xce(1)))

         value =           BA_pop_stack(txce(2))
         value = value.and.BA_pop_stack(txcp(2))
         value = value.and.BA_pop_stack(trho(2))
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: error popping stack memory',1,
     &       MA_ERR)

      end if
      return
      end


*     ***********************************
*     *			 		*
*     *   c_electron_gen_vl_potential 	*
*     *					*
*     ***********************************
      subroutine c_electron_gen_vl_potential()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      logical move,value
      integer tmp1(2)
      integer tmp2(2)
     

      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.BA_push_get(mt_dbl,(3),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call 
     >   errquit('c_electron_gen_vl_potential: out of stack memory',0,
     >       MA_ERR)

      move = .false.
      call cpsp_v_local(dcpl_mb(vl(1)),
     >                  move,
     >                  dcpl_mb(tmp1(1)),
     >                  dbl_mb(tmp2(1)))

      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >  'c_electron_gen_vl_potential: error popping stack memory',0,
     >       MA_ERR)

      return
      end



*     ***********************************
*     *			 		*
*     *   c_electron_psi_vl_ave	 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vl_ave(psi1_tag)
      implicit none
      integer psi1_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value
      integer n,nb
      integer nx,ny,nz,np,neall,psi_shift,nshift
      real*8 elocal,sum,scal1,scal2,weight,scal
      integer tmp1(2),tmp2(2)
      integer occ_tag,occ_shift

*     **** external functions ***
      integer  cpsi_data_get_chnk,cpsi_data_nsize,Pneb_convert_nb
      integer  cpsi_data_get_next
      real*8   lattice_omega,brillioun_weight
      external cpsi_data_get_chnk,cpsi_data_nsize,Pneb_convert_nb
      external cpsi_data_get_next
      external lattice_omega,brillioun_weight

      call Parallel_np(np)

      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'c_electron_psi_vl_ave: out of stack memory',0,
     &       MA_ERR)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      neall = neq(1) + neq(2)
      nshift= cpsi_data_nsize(psi_r_tag)
      occ_tag = cpsi_data_get_next(psi1_tag)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      scal=scal1*scal2

*     **** average Kohn-Sham v_local energy ****
      call Cram_c_Copy(0,dcpl_mb(vl(1)),dcpl_mb(tmp1(1)))
      call Cram_c_unpack(0,dcpl_mb(tmp1(1)))
      !call C3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
      call C3dB_cr_pfft3b(1,0,dcpl_mb(tmp1(1)))


      if (occ_tag.gt.0) then
         elocal = 0.0d0
         do nb=1,nbrillq
            !weight = brillioun_weight(Pneb_convert_nb(nb))
            weight = brillioun_weight(nb)
            psi_shift = cpsi_data_get_chnk(psi_r_tag,nb)
            occ_shift = cpsi_data_get_chnk(occ_tag,nb)

            do n=1,neall
               call C3dB_bb_Mul(1,
     >                          dcpl_mb(tmp1(1)),
     >                          dbl_mb(psi_shift),
     >                          dcpl_mb(tmp2(1)))
               call C3dB_bb_idot(1,
     >                          dbl_mb(psi_shift),
     >                          dcpl_mb(tmp2(1)),
     >                          sum)
               elocal = elocal 
     >                + sum*scal*weight*dbl_mb(occ_shift)

               psi_shift = psi_shift + nshift
               occ_shift = occ_shift + 1
            end do
         end do
      else
         elocal = 0.0d0
         do nb=1,nbrillq
            !weight = brillioun_weight(Pneb_convert_nb(nb))
            weight = brillioun_weight(nb)
            psi_shift = cpsi_data_get_chnk(psi_r_tag,nb)
             
            do n=1,neall
               call C3dB_bb_Mul(1,
     >                          dcpl_mb(tmp1(1)),
     >                          dbl_mb(psi_shift),
     >                          dcpl_mb(tmp2(1)))
               call C3dB_bb_idot(1,
     >                          dbl_mb(psi_shift),
     >                          dcpl_mb(tmp2(1)),
     >                          sum)
               psi_shift = psi_shift + nshift
               elocal = elocal + sum*scal*weight
            end do
         end do
      end if
      if (np.gt.1) call Parallel_SumAll(elocal)
      if (ispin.eq.1) elocal = 2.0d0*elocal
 
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >   'c_electron_psi_vl_ave: error popping stack memory',0,MA_ERR)

      c_electron_psi_vl_ave = elocal
      return
      end


*     ***********************************
*     *					*
*     *   c_electron_psi_vnl_ave 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vnl_ave(psi1_tag)
      implicit none
      integer psi1_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      integer n,np,nb
      integer neall,nshift,psi1_shift,tmp1_shift,occ1_shift
      real*8 enlocal,sum,weight,fion(3)
      integer tmp1_tag,occ1_tag

*     **** external functions ****
      integer  cpsi_data_push_stack,cpsi_data_get_chnk
      integer  cpsi_data_get_next
      real*8   brillioun_weight
      external cpsi_data_push_stack,cpsi_data_get_chnk
      external cpsi_data_get_next
      external brillioun_weight

      call Parallel_np(np)

      neall = neq(1) + neq(2)
      nshift = 2*npack1

      occ1_tag = cpsi_data_get_next(psi1_tag)
      tmp1_tag = cpsi_data_push_stack(nbrillq,neall,nshift)

*     **** average Kohn-Sham v_nonlocal energy ****
      call BGrsm_f_Zero_tag(tmp1_tag)
      call cpsp_v_nonlocal(ispin,ne,
     >                     psi1_tag,
     >                     tmp1_tag,
     >                    .false.,fion)
      if (occ1_tag.gt.0) then
         enlocal = 0.0d0
         do nb=1,nbrillq
            weight = brillioun_weight(nb)
            psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
            tmp1_shift = cpsi_data_get_chnk(tmp1_tag,nb)
            occ1_shift = cpsi_data_get_chnk(occ1_tag,nb)
            do n=1,neall
               call Cram_cc_idot(nb,
     >                           dbl_mb(psi1_shift),
     >                           dbl_mb(tmp1_shift),
     >                           sum)
               enlocal = enlocal - sum*weight*dbl_mb(occ1_shift)
               psi1_shift = psi1_shift + nshift
               tmp1_shift = tmp1_shift + nshift
               occ1_shift = occ1_shift + 1
            end do
         end do
      else
         enlocal = 0.0d0
         do nb=1,nbrillq
            weight = brillioun_weight(nb)
            psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
            tmp1_shift = cpsi_data_get_chnk(tmp1_tag,nb)
            do n=1,neall
               call Cram_cc_idot(nb,
     >                           dbl_mb(psi1_shift),
     >                           dbl_mb(tmp1_shift),
     >                           sum)
               enlocal = enlocal - sum*weight
               psi1_shift = psi1_shift + nshift
               tmp1_shift = tmp1_shift + nshift
            end do
         end do
      end if
      if (np.gt.1) call Parallel_SumAll(enlocal)
      if (ispin.eq.1) enlocal = 2.0d0*enlocal
 
      call cpsi_data_pop_stack(tmp1_tag)

      c_electron_psi_vnl_ave = enlocal
      return
      end


*     ***********************************
*     *					*
*     *   c_electron_psi_vnlso_ave 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vnlso_ave(psi1_tag)
      implicit none
      integer psi1_tag

#include "bafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      integer n,np,nb
      integer neall,nshift,psi1_shift,tmp1_shift,occ1_shift
      real*8 enlocal,sum,weight,fion(3)
      integer tmp1_tag,occ1_tag

*     **** external functions ****
      integer  cpsi_data_push_stack,cpsi_data_get_chnk
      integer  cpsi_data_get_next
      real*8   brillioun_weight
      external cpsi_data_push_stack,cpsi_data_get_chnk
      external cpsi_data_get_next
      external brillioun_weight

      call Parallel_np(np)

      neall = neq(1) + neq(2)
      nshift = 2*npack1

      occ1_tag = cpsi_data_get_next(psi1_tag)
      tmp1_tag = cpsi_data_push_stack(nbrillq,neall,nshift)

*     **** average Kohn-Sham v_nonlocal energy ****
      call BGrsm_f_Zero_tag(tmp1_tag)
      call cpsp_v_spin_orbit(ispin,ne,
     >                     psi1_tag,
     >                     tmp1_tag,
     >                    .false.,fion)
      if (occ1_tag.gt.0) then
         enlocal = 0.0d0
         do nb=1,nbrillq
            weight = brillioun_weight(nb)
            psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
            tmp1_shift = cpsi_data_get_chnk(tmp1_tag,nb)
            occ1_shift = cpsi_data_get_chnk(occ1_tag,nb)
            do n=1,neall
               call Cram_cc_idot(nb,
     >                           dbl_mb(psi1_shift),
     >                           dbl_mb(tmp1_shift),
     >                           sum)
               enlocal = enlocal - sum*weight*dbl_mb(occ1_shift)
               psi1_shift = psi1_shift + nshift
               tmp1_shift = tmp1_shift + nshift
               occ1_shift = occ1_shift + 1
            end do
         end do
      else
         enlocal = 0.0d0
         do nb=1,nbrillq
            weight = brillioun_weight(nb)
            psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
            tmp1_shift = cpsi_data_get_chnk(tmp1_tag,nb)
            do n=1,neall
               call Cram_cc_idot(nb,
     >                           dbl_mb(psi1_shift),
     >                           dbl_mb(tmp1_shift),
     >                           sum)
               enlocal = enlocal - sum*weight
               psi1_shift = psi1_shift + nshift
               tmp1_shift = tmp1_shift + nshift
            end do
         end do
      end if
      if (np.gt.1) call C3dB_SumAll(enlocal)
      if (ispin.eq.1) enlocal = 2.0d0*enlocal
 
      call cpsi_data_pop_stack(tmp1_tag)

      c_electron_psi_vnlso_ave = enlocal
      return
      end



*     ***********************************
*     *			 		*
*     *   c_electron_semicoreforce 	*
*     *					*
*     ***********************************

      subroutine c_electron_semicoreforce(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

     
      call c_semicore_xc_F(ispin,dbl_mb(xcp(1)),fion)

      return
      end


*     ******************************************
*     *                                        *
*     *     c_electron_gen_weighted_density    *
*     *                                        *
*     ******************************************

      subroutine c_electron_gen_weighted_density(ms,weight,rho)
      implicit none
      integer ms
      real*8  weight(*)
      real*8  rho(*)

#include "bafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer psi_shift
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

*     **** generate weighted rho ****
      psi_shift = cpsi_data_get_ptr(psi_r_tag,1,1+(ms-1)*ne(1))
      call BGrsm_hr_aSqrpy(nfft3d,ne(ms),
     >                      weight,
     >                      dbl_mb(psi_shift),
     >                      rho)
      return
      end

*    *************************************
*    *  c_electron spin orbit
*    *************************************
      logical function c_electron_spin_orbit()
#include "c_electron_common.fh"
      c_electron_spin_orbit=spin_orbit 
      return
      end

