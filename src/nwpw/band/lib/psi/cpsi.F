*
* $Id: cpsi.F,v 1.1 2001-12-14 20:02:44 bylaska Exp $
*

*     ***************************
*     *				*
*     *		cpsi_1to2	*
*     *				*
*     ***************************
      subroutine cpsi_1to2()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call zcopy(nfft3d*(ne(1)+ne(2))*nbrillioun,
     >           dcpl_mb(psi1(1)),1,
     >           dcpl_mb(psi2(1)),1)

      return
      end



*     ***************************
*     *				*
*     *		cpsi_2to1	*
*     *				*
*     ***************************
      subroutine cpsi_2to1()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

 
      call zcopy(nfft3d*(ne(1)+ne(2))*nbrillioun,
     >           dcpl_mb(psi2(1)),1,
     >           dcpl_mb(psi1(1)),1)

      return
      end

*     ***************************
*     *				*
*     *		cpsi_check	*
*     *				*
*     ***************************
      subroutine cpsi_check()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"


      call c_OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end



*     ***************************
*     *				*
*     *		c_rho_2to1	*
*     *				*
*     ***************************
      subroutine c_rho_2to1()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      call dcopy(4*nfft3d,
     >           dbl_mb(rho2(1)),1,
     >           dbl_mb(rho1(1)),1)

      call dcopy(4*nfft3d,
     >           dbl_mb(rho2_all(1)),1,
     >           dbl_mb(rho1_all(1)),1)
      return
      end

*     ***************************
*     *				*
*     *		c_dng_2to1	*
*     *				*
*     ***************************
      subroutine c_dng_2to1()
      implicit none
 
#include "mafdecls.fh"
#include "cpsi_common.fh"

      call zcopy(npack0,
     >           dcpl_mb(dng2(1)),1,
     >           dcpl_mb(dng1(1)),1)

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_1toelectron		*
*     *					*
*     ***********************************
      subroutine cpsi_1toelectron()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      call c_electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                  dbl_mb(rho1_all(1)))

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_1energy		*
*     *					*
*     ***********************************
      real*8 function cpsi_1energy()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

      call c_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))
      cpsi_1energy = c_electron_energy(dcpl_mb(psi1(1)),
     >                               dbl_mb(rho1(1)),
     >                              dcpl_mb(dng1(1)),
     >                              dbl_mb(rho1_all(1)))

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_2energy		*
*     *					*
*     ***********************************
      real*8 function cpsi_2energy()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

      call c_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      cpsi_2energy = c_electron_energy(dcpl_mb(psi2(1)),
     >                               dbl_mb(rho2(1)),
     >                              dcpl_mb(dng2(1)),
     >                              dbl_mb(rho2_all(1)))

      return
      end



*     ***********************
*     *						*
*     *		cpsi_1eorbit		*
*     *						*
*     ***********************
      real*8 function cpsi_1eorbit()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_eorbit
      external c_electron_eorbit

      cpsi_1eorbit = c_electron_eorbit(dcpl_mb(psi1(1)))

      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1ke 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1ke()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ***
      real*8 ave

      
      call cke_ave(ispin,ne,dcpl_mb(psi1(1)),ave)
      cpsi_1ke = ave
      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1vl 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1vl()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_psi_vl_ave
      external c_electron_psi_vl_ave
 
      cpsi_1vl = c_electron_psi_vl_ave(dcpl_mb(psi1(1)),dbl_mb(rho1(1)))

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_1vnl 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1vnl()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_psi_vnl_ave
      external c_electron_psi_vnl_ave
 
      cpsi_1vnl = c_electron_psi_vnl_ave(dcpl_mb(psi1(1)))
      return
      end



*     ***********************
*     *						*
*     *		c_rho_1exc		*
*     *						*
*     ***********************
      real*8 function c_rho_1exc()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** external functions ****
      real*8   c_electron_exc
      external c_electron_exc

      c_rho_1exc = c_electron_exc(dbl_mb(rho1_all(1)))
      return
      end

*     *******************
*     *					*
*     *		c_rho_1pxc	*
*     *					*
*     *******************
      real*8 function c_rho_1pxc()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_pxc
      external c_electron_pxc

      c_rho_1pxc = c_electron_pxc(dbl_mb(rho1(1)))
      return
      end


*     *******************************
*     *								*
*     *		c_dng_1ehartree         *
*     *								*
*     *******************************
      real*8 function c_dng_1ehartree()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_ehartree
      external c_electron_ehartree

      c_dng_1ehartree = c_electron_ehartree(dcpl_mb(dng1(1)))
      return
      end



*     ***********************************
*     *					*
*     *		cpsi_2toelectron		*
*     *					*
*     ***********************************
      subroutine cpsi_2toelectron()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      call c_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                   dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))

      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1get_Tgradient	*
*     *					*
*     ***********************************
      subroutine cpsi_1get_Tgradient(THpsi)
      implicit none
      real*8 THpsi(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer tmp1(2)
 
      real*8 dsum
      external dsum

      value = MA_push_get(mt_dcpl,(2*ne(1)*ne(1))*nbrillioun,
     >                   'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0)


      call c_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))
      call c_electron_gen_hml(dcpl_mb(psi1(1)),
     >                       dcpl_mb(tmp1(1)))
      call c_electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dcpl_mb(tmp1(1)),
     >                            THpsi)
      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_1get_Tradient',0)

      return
      end



*     ***********************************
*     *					*
*     *		cpsi_2get_Tgradient	*
*     *					*
*     ***********************************
      subroutine cpsi_2get_Tgradient(THpsi)
      implicit none
      complex*16 THpsi(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     *** local variables ****
      logical value
      integer tmp1(2)

      value = MA_push_get(mt_dcpl,(2*ne(1)*ne(1)*nbrillioun),
     >                   'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0)


      call c_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      call c_electron_gen_hml(dcpl_mb(psi2(1)),
     >                       dcpl_mb(tmp1(1)))
      call c_electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dcpl_mb(tmp1(1)),
     >                             THpsi)
      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_1get_Tradient',0)

      return
      end



*     ***********************************
*     *					*
*     *	    cpsi_1geodesic_transport	*
*     *					*
*     ***********************************
      subroutine cpsi_1geodesic_transport(t,H0)
      implicit none
      real*8 t
      complex*16 H0(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_transport(t,dcpl_mb(psi1(1)),H0)

      return
      end


*     ***********************************
*     *									*
*     *	    cpsi_1geodesic_Gtransport	*
*     *									*
*     ***********************************
      subroutine cpsi_1geodesic_Gtransport(t,G0)
      implicit none
      real*8 t
      complex*16 G0(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_Gtransport(t,dcpl_mb(psi1(1)),G0)

      return
      end



*     *******************************
*     *								*
*     *		cpsi_geodesic_energy 	*
*     *								*
*     *******************************
      real*8 function cpsi_geodesic_energy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      real*8 e_new
*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

     
      call c_geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))
      call c_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      e_new =  c_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)))

      cpsi_geodesic_energy = e_new
      return
      end

*     ***********************************
*     *					*
*     *		cpsi_geodesic_denergy 	*
*     *					*
*     ***********************************
      real*8 function cpsi_geodesic_denergy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_eorbit
      external c_electron_eorbit


      call c_geodesic_transport(t,dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)))
      cpsi_geodesic_denergy = 2.0d0*c_electron_eorbit(dcpl_mb(psi2(1)))

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_geodesic_final 	*
*     *					*
*     ***********************************
      subroutine cpsi_geodesic_final(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "cpsi_common.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
c     real*8 sum1,sum2
     
      call Parallel_taskid(taskid)

      call c_geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))

c*     **** Ortho Check ****     
c      call Grsm_gg_trace(npack1,(ne(1)+ne(2)),
c     >                        dcpl_mb(psi2(1)),
c     >                        dcpl_mb(psi2(1)),
c     >                        sum2)
c  
c      sum1 = dble(ne(1) + ne(2))
c      if (dabs(sum2-sum1).gt.1.0d-10) then
c         if (taskid.eq.MASTER) 
c     >      write(*,*) "Waring: Gram-Schmidt being performed on psi2"
c         call Grsm_g_MakeOrtho(npack1,ne(1),dcpl_mb(psi2(1)))
c         if (ispin.gt.1) then
c           call Grsm_g_MakeOrtho(npack1,ne(2),
c     >                           dcpl_mb(psi2(1) + ne(1)*npack1))
c         end if
c      end if
          
      return
      end

*     ***********************************
*     *									*
*     *			cpsito2_sd_update	 	*
*     *									*
*     ***********************************
      subroutine cpsi1to2_sd_update(dte)
      implicit none
      real*8 dte

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer nemax,ierr
      integer lmd(2),tmp_L(2)

      call c_electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                  dbl_mb(rho1_all(1)))

*     **** do a steepest descent step ****
      call c_electron_sd_update(dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >						  dte)

*     **** lagrange multiplier corrections ****
      nemax = ne(1)+ne(2)

*     **** allocate MA local variables ****
      value = MA_push_get(mt_dbl,(8*nemax*nemax),
     >                    'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(2*nemax*nemax),
     >                    'lmd',lmd(2),lmd(1))

      call cpsi_lmbda(ispin,ne,nemax,nfft3d,
     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
     >                 dbl_mb(lmd(1)),
     >                 dbl_mb(tmp_L(1)),ierr)

      value = value.and.MA_pop_stack(lmd(2))
      value = value.and.MA_pop_stack(tmp_L(2))
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update: stack failure', 0)
      return
      end


*     ***********************************
*     *					*
*     *		psi_1force              *
*     *					*
*     ***********************************
      subroutine psi_1force(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
*     ***** psi common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer psi1(2)
      integer psi2(2)
      integer rho1(2)
      integer rho2(2)
      integer dng1(2)
      integer dng2(2)
      integer hml(2)
      integer eig(2)
      common / psi_block / psi1,psi2,rho1,rho2,dng1,dng2,
     >                     hml,eig,ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      logical value
      integer r_grid(2),tmp(2)

*     **** external functions ****
      integer  control_version
      external control_version

c     call electron_gen_psi_r(dcpl_mb(psi1(1)))
c     call electron_gen_densities(dcpl_mb(psi1(1)),
c    >                             dbl_mb(rho1(1)),
c    >                            dcpl_mb(dng1(1)))

      call f_vlocal(dcpl_mb(dng1(1)),fion)

      if (control_version().eq.4) then
          value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',
     >                        tmp(2),tmp(1))
          value = value.and.
     >            MA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                        r_grid(2),r_grid(1))
         if (.not. value) call errquit('out of stack memory',0)

          call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                       dbl_mb(rho1(1)+(ispin-1)*2*nfft3d),
     >                       dbl_mb(tmp(1)))
          call lattice_r_grid(dbl_mb(r_grid(1)))
          call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                         dbl_mb(tmp(1)),
     >                         fion)

          value = MA_pop_stack(r_grid(2))
          value = value.and.MA_pop_stack(tmp(2))
         if (.not. value) call errquit('error popping stack memory',0)
      end if

      call f_vnonlocal(ispin,
     >                 ne,
     >                 dcpl_mb(psi1(1)),
     >                 fion)
      return
      end





*     ***********************************
*     *					*
*     *		psi_1gen_hml      	*
*     *					*
*     ***********************************
      subroutine psi_1gen_hml()
      implicit none

#include "mafdecls.fh"
*     ***** psi common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer psi1(2)
      integer psi2(2)
      integer rho1(2)
      integer rho2(2)
      integer dng1(2)
      integer dng2(2)
      integer hml(2)
      integer eig(2)
      common / psi_block / psi1,psi2,rho1,rho2,dng1,dng2,
     >                     hml,eig,ne,ispin,nfft3d,npack1,npack0

      call electron_gen_hml(dcpl_mb(psi1(1)),dbl_mb(hml(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_eigenvalue    	*
*     *					*
*     ***********************************
      real*8  function psi_eigenvalue(ms,i)
      implicit none
      integer ms
      integer i

#include "mafdecls.fh"
*     ***** psi common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer psi1(2)
      integer psi2(2)
      integer rho1(2)
      integer rho2(2)
      integer dng1(2)
      integer dng2(2)
      integer hml(2)
      integer eig(2)
      common / psi_block / psi1,psi2,rho1,rho2,dng1,dng2,
     >                     hml,eig,ne,ispin,nfft3d,npack1,npack0

      psi_eigenvalue = dbl_mb(eig(1) + (i-1) + (ms-1)*ne(1))

      return
      end

*     ***********************************
*     *					*
*     *		psi_hml		   	*
*     *					*
*     ***********************************
      real*8  function psi_hml(ms,i,j)
      implicit none
      integer ms
      integer i,j

#include "mafdecls.fh"
*     ***** psi common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer psi1(2)
      integer psi2(2)
      integer rho1(2)
      integer rho2(2)
      integer dng1(2)
      integer dng2(2)
      integer hml(2)
      integer eig(2)
      common / psi_block / psi1,psi2,rho1,rho2,dng1,dng2,
     >                     hml,eig,ne,ispin,nfft3d,npack1,npack0

      psi_hml = dbl_mb(hml(1)-1 + i 
     >                          + (j-1)*ne(ms) 
     >                          + (ms-1)*ne(1)*ne(1)) 

      return
      end


*     ***********************************
*     *					*
*     *		psi_spin_density  	*
*     *					*
*     ***********************************
      subroutine psi_spin_density(en)
      implicit none
      real*8 en(2)

#include "mafdecls.fh"
*     ***** psi common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer psi1(2)
      integer psi2(2)
      integer rho1(2)
      integer rho2(2)
      integer dng1(2)
      integer dng2(2)
      integer hml(2)
      integer eig(2)
      common / psi_block / psi1,psi2,rho1,rho2,dng1,dng2,
     >                     hml,eig,ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      integer ms,nx,ny,nz,n2ft3d
      real*8  scale,sumall

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nfft3d(1,n2ft3d)
      n2ft3d = 2*n2ft3d
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scale = lattice_omega()/dble(nx*ny*nz)

*     **** check total number of electrons ****
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(rho1(1)+(ms-1)*n2ft3d),sumall)
         en(ms) = sumall*scale
      end do
      
      return
      end

*     ***********************************
*     *					*
*     *		psi_spin2     	        *
*     *					*
*     ***********************************
      subroutine psi_spin2(Sab)
      implicit none
      real*8 Sab

#include "mafdecls.fh"
*     ***** psi common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer psi1(2)
      integer psi2(2)
      integer rho1(2)
      integer rho2(2)
      integer dng1(2)
      integer dng2(2)
      integer hml(2)
      integer eig(2)
      common / psi_block / psi1,psi2,rho1,rho2,dng1,dng2,
     >                     hml,eig,ne,ispin,nfft3d,npack1,npack0

      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi1(1)),Sab)
           
      return
      end







*     ***************************
*     *							*
*     *		cpsi_1rotate2      	*
*     *							*
*     ***************************
      subroutine cpsi_1rotate2()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables *****
      integer ms,index,i,j,shift1,shift2,nb

      call dcopy(2*nfft3d*(ne(1)+ne(2)),0.0d0,0,dcpl_mb(psi2(1)),1)

      do nb=1,nbrillioun
      do ms=1,ispin
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)
         do j=1,ne(ms)
           do i=1,ne(ms)
              index = (i-1) + (j-1)*ne(ms) + shift2
             
              call Cram_cc_zaxpy(nb,dbl_mb(hml(1)+index),
     >                           dcpl_mb(psi1(1)+(i-1+shift1)*nfft3d),
     >                           dcpl_mb(psi2(1)+(j-1+shift1)*nfft3d)) 

           end do
         end do
      end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *		cpsi_diagonalize_hml	*
*     *					*
*     ***********************************
      subroutine cpsi_diagonalize_hml()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"


*     ***** local variables ****
      logical value
      integer ms,shift1,shift2,ierr,nb,nemax
      integer work(2),rwork(2),lwork


      lwork = 2*ne(1)*ne(1)
      value = MA_push_get(mt_dbl,lwork,
     >                    'rwork',rwork(2),rwork(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,lwork,
     >                    'work',work(2),work(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_diagonalize_hml',0)


      nemax = ne(1)+ne(2)

*     ***** diagonalize the hamiltonian matrix *****
      call dcopy(2*nbrillioun*(ne(1)+ne(2)),0.0d0,0,dbl_mb(eig(1)),1)
      do nb=1,nbrillioun
      do ms=1,ispin
         shift1 = (ms-1)*ne(1)      + (nb-1)*2*nemax
         shift2 = (ms-1)*ne(1)*ne(1)+ (nb-1)*2*nemax*nemax
         call ZHEEV('V','L',
     >              ne(ms),
     >              dcpl_mb(hml(1)+shift2),
     >              ne(ms),
     >              dbl_mb(eig(1)+shift1),
     >              dcpl_mb(work(1)),
     >              lwork,
     >              dbl_mb(rwork(1)),
     >              ierr)

      end do
      end do

      
      value = MA_pop_stack(work(2))
      value = value.and.MA_pop_stack(rwork(2))
      if (.not. value) 
     > call errquit('error popping stack in cpsi_diagonalize_hml',0)


      return
      end



*     ***********************
*     *						*
*     *		cpsi_error		*
*     *						*
*     ***********************
      real*8 function cpsi_error()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      logical value
      integer k,n,nb,npack,indx
      real*8  error,sum,size
      integer tmp1(2)

      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_error',0)


      error = 0.0d0
      size =  dble(ne(1)+ne(2))
      do nb=1,nbrillioun
        call Cram_npack(nb,npack)
        do n=1, (ne(1)+ne(2))
          indx = (n-1)*nfft3d + (nb-1)*nfft3d*(ne(1)+ne(2))
          do k=1,npack
             dcpl_mb(tmp1(1)+k-1) = dcpl_mb(psi2(1)+k-1+indx) 
     >                            - dcpl_mb(psi1(1)+k-1+indx)
          end do
          call Cram_cc_dot(nb,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)

          error = error + sum
        end do
      end do
      error = dsqrt(error)/size
      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_error',0)


      cpsi_error = error
      return
      end

*     ***********************
*     *						*
*     *		c_rho_error		*
*     *						*
*     ***********************
      real*8 function c_rho_error()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      logical value
      integer k,nx,ny,nz
      real*8  error,scale
      integer tmp1(2)

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega

      value = MA_push_get(mt_dbl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in rho_error',0)


      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scale = lattice_omega()

      scale = (scale)/dble(nx*ny*nz)
*     scale = (scale)/dble(nx*ny*nz)
*     scale = (scale*scale)

      do k=1,(nfft3d)
         dbl_mb(tmp1(1)+k-1) = (dbl_mb(rho2(1)+k-1)
     >                         -dbl_mb(rho1(1)+k-1)) 
         dbl_mb(tmp1(1)+k-1) = dbl_mb(tmp1(1)+k-1) 
     >                      + (dbl_mb(rho2(1)+k-1+(ispin-1)*(nfft3d))
     >                        -dbl_mb(rho1(1)+k-1+(ispin-1)*(nfft3d))) 
      end do
      call C3dB_rr_dot(1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)),error)
      error = error*scale
*     error = dsqrt(error)

      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in rho_error',0)

      c_rho_error = error
      return
      end


*     ***********************
*     *						*
*     *		cpsi_ispin		*
*     *						*
*     ***********************
      integer function cpsi_ispin()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_ispin = ispin
      return
      end


*     *******************
*     *					*
*     *		cpsi_ne		*
*     *					*
*     *******************
      integer function cpsi_ne(ms)
      implicit none
      integer ms

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_ne = ne(ms)
      return
      end

*     ***************************
*     *				*
*     *		cpsi_initialize 	*
*     *				*
*     ***************************

      logical function cpsi_initialize()
      implicit none 

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value
      integer nemax
      real*8 sum1,sum2

*     **** external functions ****
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun


*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne(ispin,ne)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      nbrillioun = brillioun_nbrillioun()
      nemax = ne(1)+ne(2)

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2))*nbrillioun,
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2))*nbrillioun,
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1',rho1(2),rho1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2',rho2(2),rho2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'dng1',dng1(2),dng1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'dng2',dng2(2),dng2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(2*nemax*nemax*nbrillioun),
     >                    'hml',hml(2),hml(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax*nbrillioun),
     >                    'eig',eig(2),eig(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) call errquit('out of heap memory',0)

*     *****  read initial wavefunctions into psi1  ****
      call cpsi_read(ispin,ne,nbrillioun,dcpl_mb(psi1(1)))

c*     **** Ortho Check ****
c      call Grsm_gg_trace(npack1,(ne(1)+ne(2)),
c     >                        dcpl_mb(psi1(1)),
c     >                        dcpl_mb(psi1(1)),
c     >                        sum2)
c
c      sum1 = dble(ne(1) + ne(2))
c      if (dabs(sum2-sum1).gt.1.0d-10) then
c         call Parallel_taskid(taskid)
c         call Grsm_g_MakeOrtho(npack1,ne(1),dcpl_mb(psi1(1)))
c         if (ispin.gt.1) then
c           call Grsm_g_MakeOrtho(npack1,ne(2),
c     >                           dcpl_mb(psi1(1) + ne(1)*npack1))
c         end if
c         call Grsm_gg_trace(npack1,(ne(1)+ne(2)),
c     >                        dcpl_mb(psi1(1)),
c     >                        dcpl_mb(psi1(1)),
c     >                        sum2)
c         if (taskid.eq.MASTER)
c     >    write(*,*) "Waring - Gram-Schmidt being performed on psi:",
c     >               dabs(sum2-sum1)
c
c
c      end if



      cpsi_initialize = value
      return
      end



*     ***************************
*     *							*
*     *		cpsi_finalize		*
*     *							*
*     ***************************

      logical function cpsi_finalize()
      implicit none 

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      logical value

*     ***** write psi1 wavefunctions ****
      call cpsi_write(ispin,ne,nbrillioun,dcpl_mb(psi1(1)))
    
      value = MA_free_heap(eig(2))
      value = MA_free_heap(hml(2))
      value = MA_free_heap(dng2(2))
      value = MA_free_heap(dng1(2))
      value = MA_free_heap(rho2(2))
      value = MA_free_heap(rho1(2))
      value = MA_free_heap(psi2(2))
      value = MA_free_heap(psi1(2))
      value = MA_free_heap(rho2_all(2))
      value = MA_free_heap(rho1_all(2))

      cpsi_finalize = value
      return
      end




