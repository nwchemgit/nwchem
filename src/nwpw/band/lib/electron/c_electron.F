*
* $Id: c_electron.F,v 1.4 2001-12-20 00:30:51 bylaska Exp $
*

*     ***********************************
*     *					*
*     *		c_electron_init		*
*     *					*
*     ***********************************
      subroutine c_electron_init()
      implicit none

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value

*     **** external functions ****
      integer  psi_ispin,psi_ne,brillioun_nbrillioun
      external psi_ispin,psi_ne,brillioun_nbrillioun

      counter = 0

      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)

*     **** get nfft3d, and npack0 and nbrill****
      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)
      nbrill = brillioun_nbrillioun()

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,nbrill*nfft3d*(ne(1)+ne(2)),
     >                     'Hpsi_k',Hpsi_k(2),Hpsi_k(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nbrill*nfft3d*(ne(1)+ne(2)),
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'vl2',vl(2),vl(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'vc',vc(2),vc(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'xcp',xcp(2),xcp(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'xce',xce(2),xce(1))
      if (.not. value) 
     >  call errquit('c_electron_init: out of heap memory',0)

      return
      end


*     ***********************************
*     *					*
*     *		c_electron_finalize       *
*     *					*
*     ***********************************
      subroutine c_electron_finalize()
      implicit none

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value

*     **** free heap  memory ****
      value = MA_free_heap(Hpsi_k(2))
      value = value.and.
     >        MA_free_heap(psi_r(2))
      value = value.and.
     >        MA_free_heap(vl(2))
      value = value.and.
     >        MA_free_heap(vc(2))
      value = value.and.
     >        MA_free_heap(xcp(2))
      value = value.and.
     >        MA_free_heap(xce(2))
      if (.not. value) 
     >  call errquit(
     >   'c_electron_finalize: error freeing heap memory',0)
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
*     *		c_electron_run		*
*     *					*
*     ***********************************
      subroutine c_electron_run(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "c_electron_common.fh"

      counter = counter+1

      call c_electron_gen_psi_r(psi_k)
      call c_electron_gen_densities(psi_k,dn,dng,dnall)
      call c_electron_gen_scf_potentials(dn,dng,dnall)
      call c_electron_gen_Hpsi_k(psi_k)
      return
      end


*     ***********************************
*     *					*
*     *		c_electron_sd_update  	*
*     *					*
*     ***********************************
      subroutine c_electron_sd_update(psi1,psi2,dte)
      implicit none
      complex*16 psi1(*),psi2(*)
      real*8     dte

#include "mafdecls.fh"
#include "c_electron_common.fh"


      call c_electron_sd_subupdate(nfft3d,(ne(1)+ne(2)),nbrill,
     >                           psi1,psi2,dcpl_mb(Hpsi_k(1)),
     >                           dte)

      return
      end


      subroutine c_electron_sd_subupdate(nfft3d,nn,nbrill,
     >                                 psi1,psi2,Hpsi,dte)
      implicit none
      integer    nfft3d,nn,nbrill
      complex*16 psi1(nfft3d,nn,nbrill)
      complex*16 psi2(nfft3d,nn,nbrill)
      complex*16 Hpsi(nfft3d,nn,nbrill)
      real*8     dte

      integer n,nb

*     ************************************
*     **** do a steepest descent step ****
*     ************************************
      do nb=1,nbrill
      do n=1,nn
        call Cram_c_SMul(nb,(-dte),Hpsi(1,n,nb),psi2(1,n,nb))
        call Cram_cc_Sum(nb,psi2(1,n,nb),psi1(1,n,nb),psi2(1,n,nb))
      end do
      end do

      return
      end

*     ***********************************
*     *					*
*     *		c_electron_energy		*
*     *					*
*     ***********************************
      real*8 function c_electron_energy(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      integer nb,ii,ms,n1(2),n2(2),nx,ny,nz,indx,neall
      real*8  sum,eorbit,ehartr,exc,pxc,exc2,pxc2,dv,weight

*     **** external functions *****
      real*8   lattice_omega,c_coulomb_e,brillioun_weight
      external lattice_omega,c_coulomb_e,brillioun_weight
      
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
     

      dv = lattice_omega()/dble(nx*ny*nz)

      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)


*     *** get orbital energies ****
      eorbit = 0.0d0
      do nb=1,nbrill
        weight = brillioun_weight(nb)
        do ms=1,ispin
           do ii=n1(ms),n2(ms)
             indx = (ii-1)*nfft3d + (nb-1)*nfft3d*neall
             call Cram_cc_idot(nb,
     >                         psi_k(1+indx),
     >                         dcpl_mb(Hpsi_k(1)+indx),
     >                         sum)
             eorbit = eorbit + sum*weight

           end do
        end do
      end do
      call C3dB_SumAll(eorbit)
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
          

*     **** total energy ****
      c_electron_energy = eorbit + exc - ehartr - pxc

      return
      end


*     ***********************************
*     *					*
*     *		c_electron_eorbit		*
*     *					*
*     ***********************************
      real*8 function c_electron_eorbit(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer ii,ms,n1(2),n2(2),nb,indx,neall
      real*8  sum,eorbit,weight

      real*8  brillioun_weight
      external brillioun_weight

      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)

*     *** get orbital energies ****
      eorbit = 0.0d0
      do nb=1,nbrill
        weight = brillioun_weight(nb)
        do ms=1,ispin
           do ii=n1(ms),n2(ms)
             indx = (ii-1)*nfft3d + (nb-1)*nfft3d*neall
             call Cram_cc_idot(nb,psi_k(1+indx),
     >                           dcpl_mb(Hpsi_k(1)+indx),
     >                           sum)
             eorbit = eorbit + sum*weight
           end do
        end do
      end do
      call C3dB_SumAll(eorbit)
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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
     >                    dbl_mb(xcp(1)+2*nfft3d),
     >                    pxc2)
         pxc= pxc + pxc2
      end if
      pxc = pxc*dv
       
      c_electron_pxc =  pxc
      return
      end


 

*     ***********************************
*     *					*
*     *		c_electron_get_Hpsi_k	*
*     *					*
*     ***********************************
      subroutine c_electron_get_Hpsi_k(Hpsi_k_new)
      implicit none
      complex*16 Hpsi_k_new(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

      call Stfl_gg_Copy(ispin,ne,nfft3d,(ne(1)+ne(2)),
     >                  dcpl_mb(Hpsi_k(1)),
     >                  Hpsi_k_new)
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

      subroutine c_electron_get_Tgradient(psi_k,hml,THpsi_k)
      implicit none
      complex*16 psi_k(*)
      complex*16     hml(*)
      complex*16 THpsi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     ***** local variables ****
      integer n,neall

       
      neall = ne(1)+ne(2)
      n     = ne(1)

*     **** no preconditioning ***
      call Stfl_gqg_Mul(ispin,ne,nfft3d,neall,nbrill,n,
     >                  psi_k,hml,THpsi_k)
      call Stfl_ggg_Sub(ispin,ne,nfft3d,neall,nbrill,
     >                  dcpl_mb(Hpsi_k(1)),
     >                  THpsi_k,
     >                  THpsi_k)
      call Stfl_gg_dScale(ispin,ne,nfft3d,neall,nbrill,
     >                    (-1.0d0),
     >                    THpsi_k,
     >                    THpsi_k)
      return
      end

*     ***************************
*     *				*
*     *	    c_electron_gen_hml 	*
*     *				*
*     ***************************

      subroutine c_electron_gen_hml(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      complex*16 hml(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer n,neall
     
      n     = ne(1)
      neall = ne(1)+ne(2)

      call Stfl_ggq_dot(ispin,ne,nfft3d,neall,nbrill,n,
     >                  psi_k,
     >                  dcpl_mb(Hpsi_k(1)),
     >                  hml)
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

      subroutine c_electron_gen_Hpsi_k(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical move
      real*8  fion(3,1)

*     **** external functions ****
      integer  control_version
      external control_version

      move = .false.
*     ******************
*     **** get Hpsi ****
*     ******************
      call cpsi_H(ispin,ne,
     >           psi_k,
     >           dcpl_mb(psi_r(1)),
     >           dcpl_mb(vl(1)),
     >           dcpl_mb(vc(1)),
     >           dbl_mb(xcp(1)),
     >           dcpl_mb(Hpsi_k(1)),
     >           move,
     >           fion)

      call Stfl_gg_dscale(ispin,ne,nfft3d,(ne(1)+ne(2)),nbrill,
     >                    (-1.0d0),
     >                    dcpl_mb(Hpsi_k(1)),
     >                    dcpl_mb(Hpsi_k(1)))
      return
      end

*     ***************************
*     *				*
*     *	    c_electron_gen_psi_r	*
*     *				*
*     ***************************

      subroutine c_electron_gen_psi_r(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      integer neall

      neall = ne(1)+ne(2)


      call Stfl_gg_Copy(ispin,ne,nfft3d,neall,nbrill,
     >                  psi_k,
     >                  dcpl_mb(psi_r(1)))
      call Stfl_gh_fftb(ispin,ne,nfft3d,neall,nbrill,dcpl_mb(psi_r(1)))

      return
      end


*     ***********************************
*     *				        *
*     *	 c_electron_gen_densities	*
*     *				        *
*     ***********************************

      subroutine c_electron_gen_densities(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

     
*     **** local variables ****
      integer ms,nx,ny,nz,neall
      real*8  scal1,scal2
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
      scal2 = 1.0d0/lattice_omega()


*     *********************
*     **** generate dn ****
*     *********************
c     call dcopy(2*nfft3d,0.0d0,0,dn,1)
      call Stfl_hr_SumSqr(ispin,ne,nfft3d,neall,nbrill,
     >                    dcpl_mb(psi_r(1)),
     >                    dn)
      call C3dB_r_SMul(1,scal2,dn,dn)
      if (ispin.gt.1) 
     >  call C3dB_r_SMul(1,scal2,dn(1+nfft3d),dn(1+nfft3d))


        
*     **********************
*     **** generate dng ****
*     **********************
      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_densities: out of stack memory',0)

      call C3dB_rrc_Sum(1,dn,dn(1+(ispin-1)*nfft3d),dcpl_mb(tmp1(1)))
      call C3dB_c_SMul(1,scal1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
      call C3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
      call Cram_c_pack(0,dcpl_mb(tmp1(1)))
      call Cram_c_Copy(0,dcpl_mb(tmp1(1)),dng)


*     ********************************************************
*     **** generate dnall - used for semicore corrections ****
*     ********************************************************
      if (cpsp_semicore(0)) then
         call c_semicore_density(dcpl_mb(tmp1(1)))
         call C3dB_r_SMul(1,0.5d0,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
      else
         call C3dB_r_Zero(1,dcpl_mb(tmp1(1)))
      end if
      do ms=1,ispin
        call C3dB_rr_Sum(1,dn(1+(ms-1)*nfft3d),
     >                     dcpl_mb(tmp1(1)),
     >                     dnall(1+(ms-1)*nfft3d))
      end do

      value = MA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >     'c_electron_gen_densities: error popping stack memory',0)
      return
      end



*     ***********************************
*     *			 		*
*     *   c_electron_gen_scf_potentials	*
*     *					*
*     ***********************************

      subroutine c_electron_gen_scf_potentials(dn,dng,dnall)

      implicit none
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     ***** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      integer  control_gga,control_version
      external control_gga,control_version
     

      call c_coulomb_v(dng,dcpl_mb(vc(1)))

      if (control_gga().eq.0) then
         value = MA_push_get(mt_dbl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',0)
        call vxc(nfft3d,ispin,dnall,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)),
     >                    dbl_mb(tmp1(1)))
        value = MA_pop_stack(tmp1(2)) 
        if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: error popping stack memory',0)

      else
         write(*,*) "gga not yet implemented"
c        call v_bwexc(n2ft3d,ispin,dnall,
c     >                    dbl_mb(xcp(1)),
c     >                    dbl_mb(xce(1)))
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

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      logical move,value
      integer tmp1(2)
      integer tmp2(2)
     

      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl,(3),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call 
     >   errquit('electron_gen_vl_potential: out of stack memory',0)

      move = .false.
      call cpsp_v_local(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(tmp1(1)),
     >               dbl_mb(tmp2(1)))

      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >  'c_electron_gen_vl_potential: error popping stack memory',0)

      return
      end



*     ***********************************
*     *			 		*
*     *   c_electron_psi_vl_ave	 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vl_ave(psi1)
      implicit none
      complex*16 psi1(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2),nb,indx
      integer nx,ny,nz,np,neall
      real*8 elocal,sum,scal1,scal2,weight
      integer tmp1(2),tmp2(2)

*     **** external functions ***
      real*8   lattice_omega,brillioun_weight
      external lattice_omega,brillioun_weight

      call Parallel_np(np)

      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'c_electron_psi_vl_ave: out of stack memory',0)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()


*     **** average Kohn-Sham v_local energy ****
      call Cram_c_Copy(0,dcpl_mb(vl(1)),dcpl_mb(tmp1(1)))
      call Cram_c_unpack(0,dcpl_mb(tmp1(1)))
      call C3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
      elocal = 0.0d0
      do nb=1,nbrill
         weight = brillioun_weight(nb)
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            indx = (n-1)*nfft3d + (nb-1)*nfft3d*neall
            call C3dB_cc_Mul(1,
     >                       dcpl_mb(tmp1(1)),
     >                       dcpl_mb(psi_r(1)+indx),
     >                       dcpl_mb(tmp2(1)))

            call C3dB_cc_idot(1,
     >                       dcpl_mb(psi_r(1)+indx),
     >                       dcpl_mb(tmp2(1)),
     >                       sum)

            elocal = elocal + sum*scal1*scal2*weight
         end do
         end do
      end do
      if (np.gt.1) call C3dB_SumAll(elocal)
      if (ispin.eq.1) elocal = 2.0d0*elocal
 
      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >           'electron_psi_vl_ave: error popping stack memory',0)

      c_electron_psi_vl_ave = elocal
      return
      end


*     ***********************************
*     *					*
*     *   c_electron_psi_vnl_ave	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vnl_ave(psi1)
      implicit none
      complex*16 psi1(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2),np,nb
      integer nee(2),neall,indx
      real*8 enlocal,sum,weight
      integer tmp1(2),tmp2(2)

*     **** external functions ****
      real*8   brillioun_weight
      external brillioun_weight

      call Parallel_np(np)

      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'electron_psi_vl_ave: out of stack memory',0)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)

*     **** average Kohn-Sham v_nonlocal energy ****
      nee(1) = 1
      nee(2) = 0
      enlocal = 0.0d0
      do nb=1,nbrill
         weight = brillioun_weight(nb)
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            indx = (n-1)*nfft3d + (nb-1)*nfft3d*neall
            call C3dB_c_Zero(1,dcpl_mb(tmp1(1)))
            call cpsp_v_nonlocal(ispin,nee,
     >                      psi1(1+indx),
     >                      dcpl_mb(tmp1(1)),
     >                      .false.,dcpl_mb(tmp2(1)))
            call Cram_cc_idot(nb,psi1(1+indx),
     >                         dcpl_mb(tmp1(1)),
     >                         sum)
            enlocal = enlocal - sum*weight
         end do
         end do
      end do
      if (np.gt.1) call C3dB_SumAll(enlocal)
      if (ispin.eq.1) enlocal = 2.0d0*enlocal

 
      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >           'electron_psi_vl_ave: error popping stack memory',0)

      c_electron_psi_vnl_ave = enlocal
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

#include "mafdecls.fh"
#include "c_electron_common.fh"

     
      call c_semicore_xc_F(ispin,dbl_mb(xcp(1)),fion)

      return
      end




