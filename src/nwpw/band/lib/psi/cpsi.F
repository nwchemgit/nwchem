* * $Id: cpsi.F,v 1.47 2007-09-05 20:23:12 bylaska Exp $ *


*     ***************************
*     *                         *
*     *      cpsi_iptr_psi      *
*     *                         *
*     ***************************
      integer function cpsi_iptr_psi(i)
      implicit none
      integer i

#include "mafdecls.fh"
#include "cpsi_common.fh"

      integer iptr

      if (i.eq.2) then
        iptr = psi2(1)
      else
        iptr = psi1(1)
      end if
      cpsi_iptr_psi = iptr
      return
      end



*     ***************************
*     *                         *
*     *      cpsi_iptr_dn       *
*     *                         *
*     ***************************
      integer function cpsi_iptr_dn(i)
      implicit none
      integer i

#include "mafdecls.fh"
#include "cpsi_common.fh"

      integer iptr

      if (i.eq.2) then
        iptr = rho2(1)
      else
        iptr = rho1(1)
      end if
      cpsi_iptr_dn = iptr
      return
      end

*     ***************************
*     *				*
*     *		cpsi_1to2	*
*     *				*
*     ***************************
      subroutine cpsi_1to2()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call BGrsm_ff_Copy(dcpl_mb(psi1(1)),dcpl_mb(psi2(1)))
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
      call BGrsm_ff_Copy(dcpl_mb(psi2(1)),dcpl_mb(psi1(1)))
      return
      end

*     ***************************
*     *		            	*
*     *		cpsi_check     	*
*     *			        *
*     ***************************
      subroutine cpsi_check()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call Pneb_orthoCheckMake(.true.,0,0,npack1,dcpl_mb(psi1(1)))

      return
      end

*     ***************************
*     *		            	*
*     *	    cpsi_check_indx    	*
*     *			        *
*     ***************************
      subroutine cpsi_check_indx(indx)
      implicit none
      integer indx

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call Pneb_orthoCheck(indx,0,0,npack1,dcpl_mb(psi1(1)))
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
*     *                         *
*     *         c_rho_1to2      *
*     *                         *
*     ***************************
      subroutine c_rho_1to2()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      call dcopy(2*nfft3d,
     >           dbl_mb(rho1(1)),1,
     >           dbl_mb(rho2(1)),1)

      call dcopy(2*nfft3d,
     >           dbl_mb(rho1_all(1)),1,
     >           dbl_mb(rho2_all(1)),1)
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
*     *		cpsi_1toelectron	*
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



*     ***************************
*     *				*
*     *		cpsi_1eorbit	*
*     *				*
*     ***************************
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
 
      cpsi_1vl = c_electron_psi_vl_ave(dcpl_mb(psi1(1)))

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


*     ***********************************
*     *					*
*     *		cpsi_1vnlso 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1vnlso()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_psi_vnlso_ave
      external c_electron_psi_vnlso_ave
 
      cpsi_1vnlso = c_electron_psi_vnlso_ave(dcpl_mb(psi1(1)))
      return
      end


*     *******************************
*     *				    *
*     *		c_rho_1exc	    *
*     *				    *
*     *******************************
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

*     ***************************
*     *			        *
*     *		c_rho_1pxc	*
*     *				*
*     ***************************
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


*     ***********************************
*     *					*
*     *		c_dng_1ehartree         *
*     *					*
*     ***********************************
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
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer tmp1(2)
 
*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      real*8   c_electron_energy
      external Pneb_w_push_get,Pneb_w_pop_stack
      external c_electron_energy

      if (.not.Pneb_w_push_get(0,0,tmp1))
     >  call errquit("psi_1_get_Tgradient:out of stack memory",0,MA_ERR)


      call c_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))

      E0 = c_electron_energy(dcpl_mb(psi1(1)),
     >                       dbl_mb(rho1(1)),
     >                       dcpl_mb(dng1(1)),
     >                       dbl_mb(rho1_all(1)))

      call c_electron_gen_hml(dcpl_mb(psi1(1)),
     >                       dcpl_mb(tmp1(1)))

      call c_electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dcpl_mb(tmp1(1)),
     >                             THpsi)
      

      if (.not.Pneb_w_pop_stack(tmp1))
     > call errquit("psi_1_get_Tgradient:error popping stack",0,MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		cpsi_2get_Tgradient	*
*     *					*
*     ***********************************
      subroutine cpsi_2get_Tgradient(option,THpsi,E0)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     E0

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"
#include "errquit.fh"

*     *** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy


      value = MA_push_get(mt_dcpl,(2*ne(1)*ne(1)*nbrillioun),
     >                   'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0,
     &       MA_ERR)


      if (option.le.1) then
      call c_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                   dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      end if

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
     > call errquit('error popping stack memory in psi_1get_Tradient',0,
     &       MA_ERR)

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
*     *					*
*     *	    cpsi_1geodesic_Gtransport	*
*     *					*
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



*     ***********************************
*     *					*
*     *		cpsi_geodesic_energy 	*
*     *					*
*     ***********************************
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
     >                      dcpl_mb(psi2(1)))


      call c_electron_run(dcpl_mb(psi2(1)),
     >                    dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                    dbl_mb(rho2_all(1)))


      e_new =  c_electron_energy(dcpl_mb(psi2(1)),
     >                           dbl_mb(rho2(1)),
     >                           dcpl_mb(dng2(1)),
     >                           dbl_mb(rho2_all(1)))


      cpsi_geodesic_energy = e_new
      return
      end

*     ***************************************
*     *					    *
*     *		cpsi_geodesic_denergy 	    *
*     *					    *
*     ***************************************
      real*8 function cpsi_geodesic_denergy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_eorbit
      external c_electron_eorbit


      call c_geodesic_transport(t,dcpl_mb(psi1(1)),
     >                            dcpl_mb(psi2(1)))

      cpsi_geodesic_denergy = 
     >   2.0d0*c_electron_eorbit(dcpl_mb(psi2(1)))

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

      return
      end

*     ***********************************
*     *					*
*     *		cpsi1to2_sd_update	*
*     *					*
*     ***********************************
      subroutine cpsi1to2_sd_update(dte)
      implicit none
#include "errquit.fh"
      real*8 dte

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer n,neall,ierr
      integer lmd(2),tmp_L(2)


      call c_electron_run(dcpl_mb(psi1(1)),
     >                     dbl_mb(rho1(1)),
     >                    dcpl_mb(dng1(1)),
     >                     dbl_mb(rho1_all(1)))

*     **** do a steepest descent step ****

      call c_electron_sd_update(dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)),
     >	         		dte)

*     **** lagrange multiplier corrections ****
      neall = ne(1)+ne(2)
      n     = ne(1)

*     **** allocate MA local variables ****
      value = MA_push_get(mt_dcpl,(8*n*n),
     >                    'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(4*n*n*nbrillioun),
     >                    'lmd',lmd(2),lmd(1))


      call cpsi_lmbda(ispin,ne,npack1,neall,nbrillioun,n,
     >                 dcpl_mb(psi1(1)),
     >                 dcpl_mb(psi2(1)),
     >                 dte,
     >                 dcpl_mb(lmd(1)),
     >                 dcpl_mb(tmp_L(1)),ierr)

      value = value.and.MA_pop_stack(lmd(2))
      value = value.and.MA_pop_stack(tmp_L(2))
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update: stack failure', 0, MA_ERR)

      return
      end



*     ***********************************
*     *                                 *
*     *    cpsi1to2_sd_update_noscf     *
*     *                                 *
*     ***********************************
      subroutine cpsi1to2_sd_update_noscf(dte)
      implicit none
      real*8 dte

#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer n,neall,ierr
      integer lmd(2),tmp_L(2)


      call c_electron_run_noscf(dcpl_mb(psi1(1)),
     >                          dbl_mb(rho1(1)),
     >                          dcpl_mb(dng1(1)),
     >                          dbl_mb(rho1_all(1)))

*     **** do a steepest descent step ****

      call c_electron_sd_update(dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)),
     >                          dte)

*     **** lagrange multiplier corrections ****
      neall = ne(1)+ne(2)
      n     = ne(1)

*     **** allocate MA local variables ****
      value = MA_push_get(mt_dcpl,(8*n*n),
     >                    'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(4*n*n*nbrillioun),
     >                    'lmd',lmd(2),lmd(1))
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update_noscf: stack failure',0,MA_ERR)


      call cpsi_lmbda(ispin,ne,npack1,neall,nbrillioun,n,
     >                 dcpl_mb(psi1(1)),
     >                 dcpl_mb(psi2(1)),
     >                 dte,
     >                 dcpl_mb(lmd(1)),
     >                 dcpl_mb(tmp_L(1)),ierr)

      value = value.and.MA_pop_stack(lmd(2))
      value = value.and.MA_pop_stack(tmp_L(2))
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update_noscf: stack failure',0,MA_ERR)

      return
      end



*     ***************************
*     *				*
*     *		cpsi_1force     *
*     *				*
*     ***************************

      subroutine cpsi_1force(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     *** what about semicore???****
      call cpsp_f_vlocal(dcpl_mb(dng1(1)),fion)
      call cpsp_f_vnonlocal(ispin,
     >                 ne,
     >                 dcpl_mb(psi1(1)),
     >                 fion)
      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_1ke_stress         *
*     *                                 *
*     ***********************************
      subroutine cpsi_1ke_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call cke_euv(ispin,ne,dcpl_mb(psi1(1)),stress)
      return
      end

*     ***********************************
*     *                                 *
*     *         cpsi_1coulomb_stress     *
*     *                                 *
*     ***********************************
      subroutine cpsi_1coulomb_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call c_coulomb_euv(dcpl_mb(dng1(1)),stress)
      return
      end


*     ***********************************
*     *                                 *
*     *         c_rho_1exc_stress       *
*     *                                 *
*     ***********************************
      subroutine c_rho_1exc_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"
#include "errquit.fh"



*     ***** local variables ****
      logical value
      integer u,v,gga
      real*8 exc,pxc
      real*8 pi,scal,hm(3,3),tstress(3,3)
      integer trho(2),n2ft3d
      integer ms,i,j,q,nx,ny,nz,nq,indx1,indx2

*

*     **** external functions ****
      integer  control_gga
      real*8   c_rho_1exc,c_rho_1pxc,lattice_unitg,lattice_omega
      external control_gga
      external c_rho_1exc,c_rho_1pxc,lattice_unitg,lattice_omega

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     **** LDA part ****
      exc = c_rho_1exc()
      pxc = c_rho_1pxc()
      do v=1,3
      do u=1,3
         stress(u,v) = (exc-pxc)*hm(u,v)
      end do
      end do

*     **** PBE96 GGA part ****
      gga = control_gga()
      if ((gga.ge.10).and.(gga.lt.100)) then
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         call D3dB_nq(1,nq)
         n2ft3d = (nx+2)*ny*nq
         value = MA_push_get(mt_dbl,(2*n2ft3d),'trho',trho(2),trho(1))
         if (.not. value) call errquit(
     >     'c_rho_1exc_stress: out of stack memory',1, MA_ERR)


         call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(trho(1)),1)
         do ms=1,ispin
         do q=1,nq
         do j=1,ny
         do i=1,nx
            indx1 = (i-1)
     >            + (j-1)*nx
     >            + (q-1)*nx*ny
     >            + (ms-1)*nfft3d
            indx2 = (i-1)
     >            + (j-1)*(nx+2)
     >            + (q-1)*(nx+2)*ny
     >            + (ms-1)*n2ft3d
            dbl_mb(trho(1)+indx2) = dbl_mb(rho1_all(1)+indx1)

         end do
         end do
         end do
         end do
         
         call v_bwexc_euv(gga,n2ft3d,ispin,dbl_mb(trho(1)),
     >                  1.0d0,1.0d0,tstress)
         do v=1,3
         do u=1,3
          stress(u,v) = stress(u,v) + tstress(u,v)
         end do
         end do

         value = MA_pop_stack(trho(2))
         if (.not. value) call errquit(
     >     'c_rho_1exc_stress: error popping stack memory',1,
     >       MA_ERR)
      end if

      return
      end

*     ***********************************
*     *                                 *
*     *       c_rho_1semicore_stress    *
*     *                                 *
*     ***********************************
      subroutine c_rho_1semicore_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** not finished ****
      call c_semicore_euv(stress)

      return
      end



*     ***********************************
*     *                                 *
*     *         c_dng_1vlocal_stress    *
*     *                                 *
*     ***********************************

      subroutine c_dng_1vlocal_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "cpsi_common.fh"

      call cpsp_v_local_euv(dcpl_mb(dng1(1)),stress)

      return
      end

*     ***********************************
*     *                                 *
*     *        cpsi_1vnonlocal_stress   *
*     *                                 *
*     ***********************************
      subroutine cpsi_1vnonlocal_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "cpsi_common.fh"



*     ***** local variables ****
      integer u,v
      real*8 evnl
      real*8 pi,scal,hm(3,3)

*     **** external functions ****
      real*8   cpsi_1vnl,lattice_unitg
      external cpsi_1vnl,lattice_unitg

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

      call cpsp_v_nonlocal_euv_2(ispin,ne,dcpl_mb(psi1(1)),stress)
      evnl = cpsi_1vnl()
      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) - evnl*hm(u,v)
      end do
      end do

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
      integer indx

      if (spin_orbit) then
        indx=eig(1)+i-1+(nb-1)*ne(1)
      else
        indx=eig(1)+i-1+(ms-1)*ne(1)+(nb-1)*(ne(1)+ne(2))
      end if
      cpsi_eigenvalue = dbl_mb(indx) 

      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_1gen_svector       *
*     *                                 *
*     ***********************************
      subroutine cpsi_1gen_svector()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      logical value
      integer nb
      logical  BGrsm_f_Spins
      external BGrsm_f_Spins

      if (spin_orbit) then
         do nb=1,nbrillq
         value=BGrsm_f_Spins(nb,
     >                       dcpl_mb(psi1(1)),
     >                       dbl_mb(svec(1)+(nb-1)*3*ne(1)))
         end do
      else
        call dcopy(3*ne(1)*nbrillq,0.0d0,0,dbl_mb(svec(1)),1)
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         cpsi_svector            *
*     *                                 *
*     ***********************************
      real*8  function cpsi_svector(nb,i,xyz)
      implicit none
      integer nb,ms
      integer i,xyz

#include "mafdecls.fh"
#include "cpsi_common.fh"

      integer indx
      real*8  w

      if (spin_orbit) then
         indx=svec(1) + (xyz-1) 
     >                + 3*(i-1)
     >                + (nb-1)*3*ne(1)
         w = dbl_mb(indx)
      else
         w = 0.0d0
      end if
      cpsi_svector = w
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
      integer indx
      complex*16 Pneb_w_value
      external Pneb_w_value
      cpsi_hml=Pneb_w_value(ms,nb,i,j,dcpl_mb(hml(1)))
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

*     ***************************
*     *				*
*     *	      cpsi_spin2        *
*     *				*
*     ***************************
      subroutine cpsi_spin2(Sab)
      implicit none
      real*8 Sab
    
#include "mafdecls.fh"
#include "cpsi_common.fh"

      call Calculate_cpsi_spin2(nbrillioun,ispin,ne,npack1,
     >                          dcpl_mb(psi1(1)),Sab)

      return
      end


*     ***************************
*     *				*
*     *	      cpsi_1rotate2     *
*     *				*
*     ***************************
      subroutine cpsi_1rotate2()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables *****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      call Pneb_fwf_Multiply(0,0,
     >                       one,
     >                       dcpl_mb(psi1(1)),npack1,
     >                       dcpl_mb(hml(1)),
     >                       zero,
     >                       dcpl_mb(psi2(1)))

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
      integer neall

      neall = neq(1)+neq(2)

      call dcopy(nbrillq*neall,0.0d0,0,dbl_mb(eig(1)),1)

*     ***** diagonalize the hamiltonian matrix *****
      call Pneb_w_diag(0,0,dbl_mb(eig(1)),dcpl_mb(hml(1)))
      return
      end



*     ***************************
*     *				*
*     *		cpsi_error	*
*     *				*
*     ***************************
      real*8 function cpsi_error()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      logical value
      integer k,n,nb,npack,indx
      real*8  error,sum,size
      integer tmp1(2)

      value = MA_push_get(mt_dcpl,(npack1),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_error',0, MA_ERR)


      error = 0.0d0
      size =  dble(ne(1)+ne(2))
      do nb=1,nbrillioun
        call Cram_npack(nb,npack)
        do n=1, (ne(1)+ne(2))
          indx = (n-1)*npack1 + (nb-1)*npack1*(ne(1)+ne(2))
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
     > call errquit('error popping stack memory in psi_error',0, MA_ERR)

      cpsi_error = error
      return
      end

*     ***************************
*     *			        *
*     *		c_rho_error	*
*     *			        *
*     ***************************
      real*8 function c_rho_error()
      implicit none
#include "errquit.fh"

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
     >   call errquit('out of stack memory in rho_error',0, MA_ERR)


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
     > call errquit('error popping stack memory in rho_error',0, MA_ERR)

      c_rho_error = error
      return
      end

*     ***************************
*     *                         *
*     *         cpsi_nbrill     *
*     *                         *
*     ***************************
      integer function cpsi_nbrill()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_nbrill = nbrillioun
      return
      end



*     ***************************
*     *                         *
*     *         cpsi_nbrillq     *
*     *                         *
*     ***************************
      integer function cpsi_nbrillq()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_nbrillq = nbrillq
      return
      end



*     ***************************
*     *				*
*     *		cpsi_ispin      *
*     *			        *
*     ***************************
      integer function cpsi_ispin()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_ispin = ispin
      return
      end

*     ***************************
*     *                         *
*     *         cpsi_ispinq      *
*     *                         *
*     ***************************
      integer function cpsi_ispinq()
      implicit none

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_ispinq = ispinq
      return
      end


*     **********************
*     *			   *
*     *		cpsi_ne	   *
*     *			   *
*     **********************
      integer function cpsi_ne(ms)
      implicit none
      integer ms

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_ne = ne(ms)
      return
      end


*     **********************
*     *                    *
*     *         cpsi_neq   *
*     *                    *
*     **********************
      integer function cpsi_neq(ms)
      implicit none
      integer ms

#include "mafdecls.fh"
#include "cpsi_common.fh"

      cpsi_neq = neq(ms)
      return
      end



*     ***************************
*     *				*
*     *	    cpsi_initialize 	*
*     *				*
*     ***************************

      logical function cpsi_initialize(ortho)
      implicit none 
      logical ortho

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer neall,n,vers
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*30 filename
      character*30 control_input_psi
      external     control_input_psi
      logical  c_wvfnc_expander,control_spin_orbit
      external c_wvfnc_expander,control_spin_orbit
  
*     **** external functions ****
      logical  pspw_reformat_c_wvfnc,Pneb_w_allocate
      integer  brillioun_nbrillioun,psi_get_version
      integer  Pneb_ispinq,Pneb_nbrillq
      double precision util_random
      external pspw_reformat_c_wvfnc,Pneb_w_allocate
      external brillioun_nbrillioun,psi_get_version
      external Pneb_ispinq,Pneb_nbrillq
      external         util_random


*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne(ispin,ne)  !** only ne is used 
      nbrillioun = brillioun_nbrillioun()

      call Pneb_neq(neq)
      ispinq     = Pneb_ispinq()
      nbrillq    = Pneb_nbrillq()

      vers = psi_get_version()
      if ((vers.eq.3).or.(vers.eq.4)) then
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
          value= pspw_reformat_c_wvfnc(1)
        end if
      end if

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)

      neall = neq(1)+neq(2)
      n     = ne(1)
       
*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,npack1*neall*nbrillq,
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*neall*nbrillq,
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
c      value = value.and.
c     >        MA_alloc_get(mt_dcpl,(2*n*n*nbrillioun),
c     >                    'hml',hml(2),hml(1))
      value = value.and.Pneb_w_allocate(0,0,hml)
      value = value.and.
     >        MA_alloc_get(mt_dbl,(neall*nbrillq),
     >                    'eig',eig(2),eig(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*ne(1)*nbrillq),
     >                    'svec',svec(2),svec(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) 
     >  call errquit('cpsi_initialize:out of heap memory',0, MA_ERR)

c     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nogrid',
     >                  mt_log,1,psi_nogrid))
     >   psi_nogrid = .true.

      if (psi_nogrid) then

        call psi_get_header(hversion,hnfft,hunita,hispin,hne)

        if ( (hnfft(1).ne.control_ngrid(1)) .or.
     >       (hnfft(2).ne.control_ngrid(2)) .or.
     >       (hnfft(3).ne.control_ngrid(3)) ) then

        hnfft(1) = control_ngrid(1)
        hnfft(2) = control_ngrid(2)
        hnfft(3) = control_ngrid(3)
        call Parallel_taskid(taskid)
        value = rtdb_parallel(.false.)
        if (taskid.eq.MASTER) then

          filename =  control_input_psi()

          ind = index(filename,' ') - 1
          if (.not. rtdb_cput(rtdb,'c_xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_cput(rtdb,'c_xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_put(rtdb,'c_xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'c_wvfnc_expander_input: rtdb_put failed', 0, RTDB_ERR)

          write(*,*)
          write(*,*) "Grid is being converted:"
          write(*,*) "------------------------"
          write(*,*)
          write(*,*) "To turn off automatic grid conversion:"
          write(*,*)
          write(*,*) "set nwpw:psi_nogrid .false."
          write(*,*)
          value = c_wvfnc_expander(rtdb)

        end if
        value = rtdb_parallel(.true.)

        end if

      end if


*     *****  read initial wavefunctions into psi1  ****
      call cpsi_read(spin_orbit,ispin,ne,nbrillioun,dcpl_mb(psi1(1)))


*     **** Ortho Check ****
      call Pneb_orthoCheckMake(.true.,0,0,npack1,
     >                        dcpl_mb(psi1(1)))

      cpsi_initialize = value
      return
      end



*     ***************************
*     *				*
*     *	     cpsi_finalize	*
*     *				*
*     ***************************

      logical function cpsi_finalize(wpsi)
      implicit none 
      logical wpsi

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      logical value

*     **** external functions ****
      logical  Pneb_w_free
      external Pneb_w_free

*     ***** write psi1 wavefunctions ****
      if (wpsi) then
        call cpsi_write(spin_orbit,ispin,ne,nbrillioun,dcpl_mb(psi1(1)))
      end if
    
      value = MA_free_heap(eig(2))
      value = value.and.MA_free_heap(svec(2))
      value = value.and.Pneb_w_free(hml)
      value = value.and.MA_free_heap(dng2(2))
      value = value.and.MA_free_heap(dng1(2))
      value = value.and.MA_free_heap(rho2(2))
      value = value.and.MA_free_heap(rho1(2))
      value = value.and.MA_free_heap(psi2(2))
      value = value.and.MA_free_heap(psi1(2))
      value = value.and.MA_free_heap(rho2_all(2))
      value = value.and.MA_free_heap(rho1_all(2))

      cpsi_finalize = value
      return
      end

*     ***************************
*     *				*
*     *	     cpsi_band_finalize	*
*     *				*
*     ***************************

      logical function cpsi_band_finalize(wpsi)
      implicit none 
      logical wpsi

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      logical value

*     **** external functions ****
      logical  Pneb_w_free
      external Pneb_w_free

*     ***** write psi1 wavefunctions ****

      call cpsi_write2(spin_orbit,ispin,ne,nbrillioun,
     >      dcpl_mb(psi1(1)),wpsi)
    
      value = MA_free_heap(eig(2))
      value = value.and.Pneb_w_free(hml)
      value = value.and.MA_free_heap(dng2(2))
      value = value.and.MA_free_heap(dng1(2))
      value = value.and.MA_free_heap(rho2(2))
      value = value.and.MA_free_heap(rho1(2))
      value = value.and.MA_free_heap(psi2(2))
      value = value.and.MA_free_heap(psi1(2))
      value = value.and.MA_free_heap(rho2_all(2))
      value = value.and.MA_free_heap(rho1_all(2))

      cpsi_band_finalize = value
      return
      end


*     ***************************
*     *                         *
*     *      cpsi_dospsi_write  *
*     *                         *
*     ***************************

      subroutine cpsi_dospsi_write(k)
      implicit none
      integer k

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      character*30 filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

*     ***** write psi1 wavefunctions ****
      filename = 'dos'//c_index_name(k)//'.movecs'
      call cpsi_write_filename(filename,spin_orbit,ispin,ne,nbrillioun,
     >                         dcpl_mb(psi1(1)))

      return
      end



*     ***************************
*     *                         *
*     *      cpsi_dospsi_read   *
*     *                         *
*     ***************************

      subroutine cpsi_dospsi_read(k)
      implicit none
      integer k

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      character*30 filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

*     ***** write psi1 wavefunctions ****
      filename = 'dos'//c_index_name(k)//'.movecs'
      call cpsi_read_filename(filename,spin_orbit,ispin,ne,nbrillioun,
     >                         dcpl_mb(psi1(1)))

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function cpsi_spin_orbit()
#include "cpsi_common.fh"
      cpsi_spin_orbit=spin_orbit
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine cpsi_spins_output(nbb,outfile)
c      implicit none
c      integer outfile,nbb
c
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "cpsi_common.fh"
c
c      integer svec(2)
c      integer i,indx
c      logical value
c      real*8 ei,eiv,EV
c
cc    ** external **
c      integer  brillioun_nbrillioun
c      real*8   cpsi_eigenvalue,cpsi_svector
c      real*8   brillioun_k,brillioun_ks,brillioun_weight
c      external brillioun_nbrillioun
c      external cpsi_eigenvalue,cpsi_svector
c      external brillioun_k,brillioun_ks,brillioun_weight
c 
c      EV=27.2116d0
c      write(outfile,*)" E(H)               (E eV)        Spin(Sz,Sy,Sx)"
c      do i=0,ne(1)-1
c         ei=cpsi_eigenvalue(nbb,1,ne(1)-i)
c         eiv=ei*EV
c         write(outfile,111)
c     >      ei,eiv,  
c     >      cpsi_svector(nbb,i,3),
c     >      cpsi_svector(nbb,i,2),
c     >      cpsi_svector(nbb,i,1)
c      end do
c      call flush(outfile)
c      if (.not.value) then
c       call errquit("trouble popping stack",0,0)  
c      end if
c 108  format(/' Brillouin zone point ',i3,
c     >       /' weight = ',f10.6,
c     >       /' k      = <',3f8.3,'> . b= <b1,b2,b3> ',
c     >       /'        = <',3f8.3,'>') 
c 111  format(D16.7,' (',F8.3,' eV) (',F8.3,' zhat + ',F8.3,' yhat +',
c     > F8.3,' xhat )') 
c      return
c      end            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cpsi_output21(outfile,nb,w,k)
      implicit none
      integer outfile,nb,k
      real*8 w(*)
#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
      integer i,indx,ii
      logical value
      real*8 ei,eiv,EV
c    ** external **
      logical BGrsm_f_Spins 
      integer  brillioun_nbrillioun
      real*8 cpsi_eigenvalue
      real*8   brillioun_k,brillioun_ks,brillioun_weight
      external brillioun_nbrillioun
      external cpsi_eigenvalue
      external brillioun_k,brillioun_ks,brillioun_weight
      external BGrsm_f_Spins
 
      EV=27.2116
      do i=0,ne(1)-1
         ii=ne(1)-i
         ei=cpsi_eigenvalue(nb,1,ii)
         eiv=ei*EV
         write(outfile,111)
     >      ei,eiv,w(ii+(k-1)*ne(1))  
      end do
      call flush(outfile)
 108  format(/' Brillouin zone point ',i3,
     >       /' weight = ',f10.6,
     >       /' k      = <',3f8.3,'> . b= <b1,b2,b3> ',
     >       /'        = <',3f8.3,'>') 
 111  FORMAT(1(E18.7,' (',F8.3,'eV)  dplot weight=',F8.3))
      return
      end            

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cpsi_output22(outfile,nb)
      implicit none
      integer outfile,nb
#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
      integer i,indx,ii,iend
      real*8 ei,eiv,EV
c    ** external **
      real*8 cpsi_eigenvalue
      external cpsi_eigenvalue
 
      iend=ne(1)-1
      EV=27.2116
      do i=0,iend
         ii=ne(1)-i
         ei=cpsi_eigenvalue(nb,1,ii)
         eiv=ei*EV
         write(outfile,111)
     >      ei,eiv  
      end do
      call flush(outfile)
 108  format(/' Brillouin zone point ',i3,
     >       /' weight = ',f10.6,
     >       /' k      = <',3f8.3,'> . b= <b1,b2,b3> ',
     >       /'        = <',3f8.3,'>') 
 111  FORMAT(E18.7,' (',F8.3,'eV)')
      return
      end            

