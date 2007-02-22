*
* $Id: c_electron.F,v 1.29 2007-02-22 18:40:01 d3p708 Exp $
*


*     **********************************
*     *                                *
*     *      c_electron_iptr_psir      *
*     *                                *
*     **********************************
      integer function c_electron_iptr_psir()
      implicit none
      integer i

#include "mafdecls.fh"
#include "c_electron_common.fh"

      integer iptr

      iptr = psi_r(1)
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

#include "mafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      logical  cpsi_spin_orbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_spin_orbit

      counter = 0

      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      spin_orbit= cpsi_spin_orbit()

*     **** get nfft3d, and npack0 and nbrill****
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      call C3dB_nfft3d(1,nfft3d)
      nbrill = brillioun_nbrillioun()

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,nbrill*npack1*(ne(1)+ne(2)),
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
     >        MA_alloc_get(mt_dcpl,2*nfft3d,
     >                     'vall',vall(2),vall(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'xcp',xcp(2),xcp(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'xce',xce(2),xce(1))
      if (.not. value) 
     >  call errquit('c_electron_init: out of heap memory',0, MA_ERR)

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
     >        MA_free_heap(vall(2))
      value = value.and.
     >        MA_free_heap(vc(2))
      value = value.and.
     >        MA_free_heap(xcp(2))
      value = value.and.
     >        MA_free_heap(xce(2))
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
*     *                                 *
*     *         c_electron_run_orb        *
*     *                                 *
*     ***********************************
      subroutine c_electron_run_orb(nb,i,psi_k)
      implicit none
      integer    nb,i
      complex*16 psi_k(*)

      call c_electron_gen_psi_r_orb(nb,i,psi_k)    
      call c_electron_gen_Hpsi_k_orb(nb,i,psi_k)

      return
      end


*     ***********************************
*     *                                 *
*     *    c_electron_get_gradient_orb  *
*     *                                 *
*     ***********************************

      subroutine c_electron_get_gradient_orb(nb,i,Horb)
      implicit none
      integer nb,i,indx
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

      indx=Hpsi_k(1)+(i-1)*npack1+(nb-1)*(ne(1)+ne(2))*npack1
      call Cram_c_Copy(nb,
     >                 dcpl_mb(indx),
     >                 Horb)
      if (spin_orbit) then
        indx=indx+ne(1)*npack1
        call Cram_c_Copy(nb,
     >                 dcpl_mb(indx),
     >                 Horb(npack1+1))
      end if
      return
      end


*     ***********************************
*     *					*
*     *		c_electron_sd_update  	*
*     *					*
*     ***********************************
*     ***********************************
      subroutine c_electron_sd_update(psi1,psi2,dte)
      implicit none
      complex*16 psi1(*),psi2(*)
      real*8     dte

#include "mafdecls.fh"
#include "c_electron_common.fh"


      call c_electron_sd_subupdate(npack1,(ne(1)+ne(2)),nbrill,
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

cccccccccccccccccc


*     ***********************************
*     *					*
*     *		c_electron_energy	*
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
             indx = (ii-1)*npack1 + (nb-1)*npack1*neall
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
*     *		c_electron_eorbit	*
*     *					*
*     ***********************************
      real*8 function c_electron_eorbit(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      real*8  eorbit

      call Pneb_ff_traceall(0,0,psi_k,dcpl_mb(Hpsi_k(1)),npack1,
     >eorbit)
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
*     *					*
*     *		c_electron_get_Hpsi_k	*
*     *					*
*     ***********************************
      subroutine c_electron_get_Hpsi_k(Hpsi_k_new)
      implicit none
      complex*16 Hpsi_k_new(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

      call BGrsm_ff_Copy(dcpl_mb(Hpsi_k(1)),
     >                   Hpsi_k_new)
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

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      call Pneb_fwf_Multiply(0,0,
     >                       one,
     >                       psi_k,npack1,
     >                       hml,
     >                       zero,
     >                       THpsi_k)

      call BGrsm_fff_sub(dcpl_mb(Hpsi_k(1)),
     >                   THpsi_k,
     >                   THpsi_k)

      call BGrsm_ff_dScale(-1.0d0,
     >                     THpsi_k,
     >                     THpsi_k)
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
 

      call Pneb_ffw_multiply(0,0,
     >                  psi_k,
     >                  dcpl_mb(Hpsi_k(1)),
     >                  npack1,
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

      call BGrsm_ff_dScale((-1.0d0),
     >                     dcpl_mb(Hpsi_k(1)),
     >                     dcpl_mb(Hpsi_k(1)))
      return
      end





*     ***********************************
*     *                                 *
*     *     c_electron_gen_Hpsi_k_orb   *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_Hpsi_k_orb(nb,n,psi_k)
      implicit none
      integer nb,n
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      integer ms,index1,index2

      if (n.le.ne(1)) then
        ms=1
      else
        ms=2
      end if
      index1  = (n-1)*nfft3d  + (nb-1)*(ne(1)+ne(2))*nfft3d
      index2  = (n-1)*npack1  + (nb-1)*(ne(1)+ne(2))*npack1

*     ******************
*     **** get Hpsi ****
*     ******************
      call cpsi_Horb(nb,ispin,ms,
     >              dcpl_mb(psi_r(1)),
     >              dcpl_mb(vall(1)),
     >              psi_k(index2+1),
     >              dcpl_mb(psi_r(1)+index1),
     >              dcpl_mb(Hpsi_k(1)+index2))

      call Cram_c_SMul(nb,(-1.0d0),
     >                 dcpl_mb(Hpsi_k(1)+index2),
     >                 dcpl_mb(Hpsi_k(1)+index2))

      if (spin_orbit) then
        index2=index2+ne(1)*npack1
        call Cram_c_SMul(nb,(-1.0d0),
     >                 dcpl_mb(Hpsi_k(1)+index2),
     >                 dcpl_mb(Hpsi_k(1)+index2))
      end if               
      return
      end






*     *******************************
*     *				    *
*     *	    c_electron_gen_psi_r    *
*     *				    *
*     *******************************

      subroutine c_electron_gen_psi_r(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

      call BGrsm_fh_fftb(psi_k,
     >                   dcpl_mb(psi_r(1)))
      return
      end


*     ***********************************
*     *                                 *
*     *     c_electron_gen_psi_r_orb    *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_psi_r_orb(nb,n,psi_k)
      implicit none
      integer    nb,n
      complex*16 psi_k(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer shiftk,shiftr

      shiftk =        1 + (n-1)*npack1 + (nb-1)*(ne(1)+ne(2))*npack1
      shiftr = psi_r(1) + (n-1)*nfft3d + (nb-1)*(ne(1)+ne(2))*nfft3d

      call Cram_c_Copy(nb,psi_k(shiftk),
     >                    dcpl_mb(shiftr))

      call Cram_c_unpack(nb, dcpl_mb(shiftr))
      call C3dB_cr_fft3b(1,  dcpl_mb(shiftr))

      if (spin_orbit) then
        shiftk=shiftk+ne(1)*npack1
        shiftr=shiftr+ne(1)*nfft3d
        call Cram_c_Copy(nb,psi_k(shiftk),
     >                    dcpl_mb(shiftr))
        call Cram_c_unpack(nb, dcpl_mb(shiftr))
        call C3dB_cr_fft3b(1,  dcpl_mb(shiftr))
      end if                                   
      return
      end


*     ***********************************
*     *				        *
*     *	    c_electron_gen_densities	*
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
#include "errquit.fh"

     
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
      call BGrsm_hr_SumSqr(dcpl_mb(psi_r(1)),
     >                     dn)
      call C3dB_r_SMul(1,scal2,dn,dn)
      if (ispin.gt.1) 
     >  call C3dB_r_SMul(1,scal2,dn(1+nfft3d),dn(1+nfft3d))


        
*     **********************
*     **** generate dng ****
*     **********************
      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_densities: out of stack memory',0, MA_ERR)

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
     >     'c_electron_gen_densities: error popping stack memory',0,
     &       MA_ERR)
      return
      end


*     ***********************************
*     *                                 *
*     *        c_electron_gen_vall      *
*     *                                 *
*     ***********************************

      subroutine c_electron_gen_vall()
      implicit none


#include "mafdecls.fh"
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
      call Cram_cc_Sum(0,dcpl_mb(vall(1)),
     >                   dcpl_mb(vc(1)),
     >                   dcpl_mb(vall(1)))

*     **** fourier transform k-space potentials ****
      call Cram_c_unpack(0,dcpl_mb(vall(1)))
      call C3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     **** add xcp to vall ****
      do ms=ispin,1,-1
        call C3dB_rc_Sum(1,
     >                   dbl_mb(xcp(1)+(ms-1)*nfft3d),
     >                   dcpl_mb(vall(1)),
     >                   dcpl_mb(vall(1)+(ms-1)*nfft3d))
      end do

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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"
#include "errquit.fh"


*     ***** local variables ****
      logical value
      integer tmp1(2)
      integer trho(2),txcp(2),txce(2),n2ft3d
      integer ms,i,j,k,q,nx,ny,nz,nq,indx1,indx2,gga

*     **** external functions ****
      integer  control_gga,control_version
      external control_gga,control_version
     

      call c_coulomb_v(dng,dcpl_mb(vc(1)))

      gga = control_gga()

      !*** LDA calculation ****
      if (gga.eq.0) then
         value = MA_push_get(mt_dbl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',0, MA_ERR)
        call vxc(nfft3d,ispin,dnall,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)),
     >                    dbl_mb(tmp1(1)))
        value = MA_pop_stack(tmp1(2)) 
        if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: error popping stack memory',0,
     &       MA_ERR)


      !*** this code needs to be optimized!!! ****
      !*** GGA calculation ****
      else
         !call D3dB_nx(1,nx)
         !call D3dB_ny(1,ny)
         !call D3dB_nz(1,nz)
         !call D3dB_nq(1,nq)
         !n2ft3d = (nx+2)*ny*nq
         call D3dB_n2ft3d(1,n2ft3d)
         value = MA_push_get(mt_dbl,(2*n2ft3d),'trho',trho(2),trho(1))
         value = value.and.
     >           MA_push_get(mt_dbl,(2*n2ft3d),'txcp',txcp(2),txcp(1))
         value = value.and.
     >           MA_push_get(mt_dbl,(2*n2ft3d),'txce',txce(2),txce(1))
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',1, MA_ERR)

         !call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(trho(1)),1)

c         do ms=1,ispin
c         do q=1,nq
c         do j=1,ny
c         do i=1,nx
c            indx1 = (i-1) 
c     >            + (j-1)*nx     
c     >            + (q-1)*nx*ny     
c     >            + (ms-1)*nfft3d
c            indx2 = (i-1) 
c     >            + (j-1)*(nx+2) 
c     >            + (q-1)*(nx+2)*ny 
c     >            + (ms-1)*n2ft3d
c            dbl_mb(trho(1)+indx2) = dnall(indx1+1)
c         end do
c         end do
c         end do
c         end do

         do ms=1,ispin
           call C3dB_D3dB_r_Copy(1,dnall(1+(ms-1)*nfft3d),
     >                             dbl_mb(trho(1)+(ms-1)*n2ft3d))
         end do

         call v_bwexc(gga,n2ft3d,ispin,
     >                    dbl_mb(trho(1)),
     >                    1.0d0,1.0d0,
     >                    dbl_mb(txcp(1)),
     >                    dbl_mb(txce(1)))

         do ms=1,ispin
           call D3dB_C3dB_r_Copy(1,dbl_mb(txcp(1)+(ms-1)*n2ft3d),
     >                             dbl_mb(xcp(1) +(ms-1)*nfft3d))
         end do
         call D3dB_C3dB_r_Copy(1,dbl_mb(txce(1)),
     >                           dbl_mb(xce(1)))

c         do ms=1,ispin
c         do q=1,nq
c         do j=1,ny
c         do i=1,nx
c            indx1 = (i-1) 
c     >            + (j-1)*nx     
c     >            + (q-1)*nx*ny     
c     >            + (ms-1)*nfft3d
c            indx2 = (i-1) 
c     >            + (j-1)*(nx+2) 
c     >            + (q-1)*(nx+2)*ny 
c     >            + (ms-1)*n2ft3d
c            dbl_mb(xcp(1)+indx1) = dbl_mb(txcp(1)+indx2)
c         end do
c         end do
c         end do
c         end do
c
c         do q=1,nq
c         do j=1,ny
c         do i=1,nx
c            indx1 = (i-1) 
c     >            + (j-1)*nx     
c     >            + (q-1)*nx*ny     
c            indx2 = (i-1) 
c     >            + (j-1)*(nx+2) 
c     >            + (q-1)*(nx+2)*ny 
c            dbl_mb(xce(1)+indx1) = dbl_mb(txce(1)+indx2)
c         end do
c         end do
c         end do

         value =           MA_pop_stack(txce(2))
         value = value.and.MA_pop_stack(txcp(2))
         value = value.and.MA_pop_stack(trho(2))
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
#include "errquit.fh"

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
     >   errquit('electron_gen_vl_potential: out of stack memory',0,
     &       MA_ERR)

      move = .false.
      call cpsp_v_local(dcpl_mb(vl(1)),
     >                  move,
     >                  dcpl_mb(tmp1(1)),
     >                  dbl_mb(tmp2(1)))

      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >  'c_electron_gen_vl_potential: error popping stack memory',0,
     &       MA_ERR)

      return
      end



*     ***********************************
*     *			 		*
*     *   c_electron_psi_vl_ave	 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vl_ave(psi1)
      implicit none
#include "errquit.fh"
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
     >            'c_electron_psi_vl_ave: out of stack memory',0,
     &       MA_ERR)

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
     >           'electron_psi_vl_ave: error popping stack memory',0,
     &       MA_ERR)

      c_electron_psi_vl_ave = elocal
      return
      end


*     ***********************************
*     *					*
*     *   c_electron_psi_vnl_ave 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vnl_ave(psi1)
      implicit none
#include "errquit.fh"
      complex*16 psi1(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2),np,nb
      integer neall,indx
      real*8 enlocal,sum,weight
      integer tmp1(2),tmp2(2)

*     **** external functions ****
      real*8   brillioun_weight
      external brillioun_weight

      call Parallel_np(np)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)

      value = MA_push_get(mt_dcpl,(npack1*neall*nbrill),
     >                    'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'electron_psi_vl_ave: out of stack memory',0,
     &       MA_ERR)


*     **** average Kohn-Sham v_nonlocal energy ****
      call BGrsm_f_Zero(dcpl_mb(tmp1(1)))
      call cpsp_v_nonlocal(ispin,ne,
     >                     psi1,
     >                     dcpl_mb(tmp1(1)),
     >                    .false.,dcpl_mb(tmp2(1)))
      enlocal = 0.0d0
      do nb=1,nbrill
         weight = brillioun_weight(nb)
         do ms=1,ispin
         do n=n1(ms),n2(ms)
         indx = (n-1)*npack1 + (nb-1)*npack1*neall

            call Cram_cc_idot(nb,psi1(1+indx),
     >                         dcpl_mb(tmp1(1)+indx),
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
     >           'electron_psi_vl_ave: error popping stack memory',0,
     &       MA_ERR)

      c_electron_psi_vnl_ave = enlocal
      return
      end



*     ***********************************
*     *					*
*     *   c_electron_psi_vnlso_ave 	*
*     *					*
*     ***********************************

      real*8 function c_electron_psi_vnlso_ave(psi1)
      implicit none
#include "errquit.fh"
      complex*16 psi1(*)

#include "mafdecls.fh"
#include "c_electron_common.fh"


*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2),np,nb
      integer neall,indx
      real*8 enlocal,sum,weight
      integer tmp1(2),tmp2(2)

*     **** external functions ****
      real*8   brillioun_weight
      external brillioun_weight

      if (.not.spin_orbit) then
         c_electron_psi_vnlso_ave=0.0d0
         return
      end if

      call Parallel_np(np)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)

      value = MA_push_get(mt_dcpl,(npack1*neall*nbrill),
     >                    'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'electron_psi_vl_ave: out of stack memory',0,
     &       MA_ERR)


*     **** average Kohn-Sham v_nonlocal energy ****
      call BGrsm_f_Zero(dcpl_mb(tmp1(1)))
      call cpsp_v_spin_orbit(ispin,ne,
     >                     psi1,
     >                     dcpl_mb(tmp1(1)),
     >                    .false.,dcpl_mb(tmp2(1)))

      enlocal = 0.0d0
      do nb=1,nbrill
         weight = brillioun_weight(nb)
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            indx = (n-1)*npack1 + (nb-1)*npack1*neall

            call Cram_cc_idot(nb,psi1(1+indx),
     >                         dcpl_mb(tmp1(1)+indx),
     >                         sum)
            enlocal = enlocal - sum*weight
         end do
         end do
      end do
      if (np.gt.1) call C3dB_SumAll(enlocal)
      

 
      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >           'electron_psi_vl_ave: error popping stack memory',0,
     &       MA_ERR)

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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
#include "c_electron_common.fh"

*     **** local variables ****
      integer ishift

*     **** generate weighted rho ****
      ishift = (ms-1)*ne(1)*nfft3d
      call BGrsm_hr_aSqrpy(nfft3d,ne(ms),
     >                      weight,
     >                      dcpl_mb(psi_r(1)+ishift),
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

