*
* $Id: cpsi.F,v 1.7 2002-01-07 03:04:58 bylaska Exp $
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
*     *		            		*
*     *		cpsi_check      	*
*     *				            *
*     ***************************
      subroutine cpsi_check()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"


      call c_OrthoCheck(ispin,ne,nbrillioun,dcpl_mb(psi1(1)))  
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

      call dcopy(2*nfft3d,
     >           dbl_mb(rho2(1)),1,
     >           dbl_mb(rho1(1)),1)

      call dcopy(2*nfft3d,
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
      subroutine cpsi_1get_Tgradient(THpsi,E0)
      implicit none
      complex*16 THpsi(*)
      real*8     E0

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer tmp1(2)
 
*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

      value = MA_push_get(mt_dcpl,(2*ne(1)*ne(1))*nbrillioun,
     >                   'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0)


      call c_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))

      E0 = c_electron_energy(dcpl_mb(psi1(1)),
     >                       dbl_mb(rho1(1)),
     >                       dcpl_mb(dng1(1)),
     >                       dbl_mb(rho1_all(1)))

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
      subroutine cpsi_2get_Tgradient(THpsi,E0)
      implicit none
      complex*16 THpsi(*)
      real*8     E0

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     *** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy


      value = MA_push_get(mt_dcpl,(2*ne(1)*ne(1)*nbrillioun),
     >                   'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0)


      call c_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))

      E0 = c_electron_energy(dcpl_mb(psi2(1)),
     >                       dbl_mb(rho2(1)),
     >                       dcpl_mb(dng2(1)),
     >                       dbl_mb(rho2_all(1)))

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

      real*8 sum
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
      integer n,neall,ierr
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
      neall = ne(1)+ne(2)
      n     = ne(1)

*     **** allocate MA local variables ****
      value = MA_push_get(mt_dcpl,(8*n*n),
     >                    'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n*nbrillioun),
     >                    'lmd',lmd(2),lmd(1))

      call cpsi_lmbda(ispin,ne,nfft3d,neall,nbrillioun,n,
     >                 dcpl_mb(psi1(1)),
     >                 dcpl_mb(psi2(1)),
     >                 dte,
     >                 dcpl_mb(lmd(1)),
     >                 dcpl_mb(tmp_L(1)),ierr)

      value = value.and.MA_pop_stack(lmd(2))
      value = value.and.MA_pop_stack(tmp_L(2))
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update: stack failure', 0)
      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1force              *
*     *					*
*     ***********************************
      subroutine cpsi_1force(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     *** what about semicore???****
c      call csp_f_vlocal(dcpl_mb(dng1(1)),fion)
      call cpsp_f_vnonlocal(ispin,
     >                 ne,
     >                 dcpl_mb(psi1(1)),
     >                 fion)
      return
      end





*     ***********************************
*     *					*
*     *		cpsi_1gen_hml      	*
*     *					*
*     ***********************************
      subroutine cpsi_1gen_hml()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call c_electron_gen_hml(dcpl_mb(psi1(1)),dcpl_mb(hml(1)))

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_eigenvalue    	*
*     *					*
*     ***********************************
      real*8  function cpsi_eigenvalue(nb,ms,i)
      implicit none
      integer nb,ms
      integer i

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_eigenvalue = dbl_mb(eig(1) 
     >                        + (i-1) 
     >                        + (ms-1)*ne(1)
     >                        + (nb-1)*(ne(1)+ne(2)))

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_hml		   	*
*     *					*
*     ***********************************
      complex*16  function cpsi_hml(nb,ms,i,j)
      implicit none
      integer nb,ms
      integer i,j

#include "mafdecls.fh"
#include "cpsi_common.fh"


      cpsi_hml = dcpl_mb(hml(1)-1 + i 
     >                          + (j-1)*ne(ms) 
     >                          + (ms-1)*ne(1)*ne(1)
     >                          + (nb-1)*2*ne(1)*ne(1)) 

      return
      end


*     ***********************************
*     *					*
*     *		cpsi_spin_density  	*
*     *					*
*     ***********************************
      subroutine cpsi_spin_density(en)
      implicit none
      real*8 en(2)

#include "mafdecls.fh"
#include "cpsi_common.fh"


*     **** local variables ****
      integer ms,nx,ny,nz
      real*8  scale,sumall

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scale = lattice_omega()/dble(nx*ny*nz)

*     **** check total number of electrons ****
      do ms =1,ispin
         call C3dB_r_dsum(1,dbl_mb(rho1(1)+(ms-1)*nfft3d),sumall)
         en(ms) = sumall*scale
      end do
      
      return
      end

c*     ***********************************
c*     *					*
c*     *		psi_spin2     	        *
c*     *					*
c*     ***********************************
c      subroutine psi_spin2(Sab)
c      implicit none
c      real*8 Sab
c
c#include "mafdecls.fh"
c#include "cpsi_common.fh"
c
c      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi1(1)),Sab)
c           
c      return
c      end







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
      integer n,neall

      neall = ne(1) +ne(2)
      n     = ne(1)

      call Stfl_gqg_Mul(ispin,ne,nfft3d,neall,nbrillioun,n,
     >                  dcpl_mb(psi1(1)),
     >                  dcpl_mb(hml(1)),
     >                  dcpl_mb(psi2(1)))


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
      integer neall,n

      neall = ne(1)+ne(2)
      n     = ne(1)

*     ***** diagonalize the hamiltonian matrix *****
      call dcopy(nbrillioun*neall,0.0d0,0,dbl_mb(eig(1)),1)
      call Stfl_q_diag(ispin,ne,neall,nbrillioun,n,
     >     dbl_mb(eig(1)),
     >     dcpl_mb(hml(1)))
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
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value
      integer neall,n

*     **** external functions ****
      logical  pspw_reformat_c_wvfnc
      integer  brillioun_nbrillioun,psi_get_version
      external pspw_reformat_c_wvfnc
      external brillioun_nbrillioun,psi_get_version

*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne(ispin,ne)

      nbrillioun = brillioun_nbrillioun()
      if (psi_get_version().eq.3) then
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
          value= pspw_reformat_c_wvfnc(nbrillioun)
        end if
      end if

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      nbrillioun = brillioun_nbrillioun()
      neall = ne(1)+ne(2)
      n     = ne(1)

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,nfft3d*neall*nbrillioun,
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*neall*nbrillioun,
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho1',rho1(2),rho1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho2',rho2(2),rho2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'dng1',dng1(2),dng1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'dng2',dng2(2),dng2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(2*n*n*nbrillioun),
     >                    'hml',hml(2),hml(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(neall*nbrillioun),
     >                    'eig',eig(2),eig(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) call errquit('out of heap memory',0)

*     *****  read initial wavefunctions into psi1  ****
      call cpsi_read(ispin,ne,nbrillioun,dcpl_mb(psi1(1)))


*     **** Ortho Check ****
      call Stfl_g_OrthoCheckMake(ispin,ne,nfft3d,neall,nbrillioun,
     >                        dcpl_mb(psi1(1)))

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




