* * $Id$ *


*     ***************************
*     *                         *
*     *      cpsi_iptr_psi      *
*     *                         *
*     ***************************
      integer function cpsi_iptr_psi(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "cpsi_common.fh"

      integer iptr

      if (i.eq.2) then
        iptr = psi2_tag
      else
        iptr = psi1_tag
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      !call BGrsm_ff_Copy(dcpl_mb(psi1(1)),dcpl_mb(psi2(1)))
      call cpsi_data_copyall(psi1_tag,psi2_tag)
      return
      end
*     ***************************
*     *				*
*     *		cpsi_2to1	*
*     *				*
*     ***************************
      subroutine cpsi_2to1()
      implicit none
#include "bafdecls.fh"
#include "cpsi_common.fh"

      !call BGrsm_ff_Copy(dcpl_mb(psi2(1)),dcpl_mb(psi1(1)))
      call cpsi_data_copyall(psi2_tag,psi1_tag)
      return
      end

*     ***************************
*     *				*
*     *		cpsi_swap12	*
*     *				*
*     ***************************
      subroutine cpsi_swap12()
      implicit none
#include "bafdecls.fh"
#include "cpsi_common.fh"
      integer itag
      itag = psi1_tag
      psi1_tag = psi2_tag
      psi2_tag = itag
      return
      end



*     ***************************
*     *		            	*
*     *		cpsi_check     	*
*     *			        *
*     ***************************
      subroutine cpsi_check()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call Pneb_orthoCheckMake_tag(.true.,0,0,npack1,psi1_tag)

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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call Pneb_orthoCheck_tag(indx,0,0,npack1,psi1_tag)
      return
      end



*     ***************************
*     *				*
*     *		c_rho_2to1	*
*     *				*
*     ***************************
      subroutine c_rho_2to1()
      implicit none

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
*     *                         *
*     *         c_rho_swap12    *
*     *                         *
*     ***************************
      subroutine c_rho_swap12()
      implicit none
#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"
      integer itag
      itag    = rho1(1)
      rho1(1) = rho2(1)
      rho2(1) = itag

      itag        = rho1_all(1)
      rho1_all(1) = rho2_all(1)
      rho2_all(1) = itag
      return
      end


*     ***************************
*     *				*
*     *		c_dng_2to1	*
*     *				*
*     ***************************
      subroutine c_dng_2to1()
      implicit none
 
#include "bafdecls.fh"
#include "cpsi_common.fh"

      call zcopy(npack0,
     >           dcpl_mb(dng2(1)),1,
     >           dcpl_mb(dng1(1)),1)

      return
      end

*     ***************************
*     *                         *
*     *         c_dng_swap12    *
*     *                         *
*     ***************************
      subroutine c_dng_swap12()
      implicit none
#include "bafdecls.fh"
#include "cpsi_common.fh"
      integer itag
      itag    = dng1(1)
      dng1(1) = dng2(1)
      dng2(1) = itag
      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1toelectron	*
*     *					*
*     ***********************************
      subroutine cpsi_1toelectron()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      call c_electron_run(psi1_tag,
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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

      call c_electron_run(psi1_tag,
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))
      cpsi_1energy = c_electron_energy(psi1_tag,
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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

      call c_electron_run(psi2_tag,
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      cpsi_2energy = c_electron_energy(psi2_tag,
     >                               dbl_mb(rho2(1)),
     >                              dcpl_mb(dng2(1)),
     >                               dbl_mb(rho2_all(1)))

      return
      end



*     ***************************
*     *				*
*     *		cpsi_1eorbit	*
*     *				*
*     ***************************
      real*8 function cpsi_1eorbit()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_eorbit
      external c_electron_eorbit

      cpsi_1eorbit = c_electron_eorbit(psi1_tag)

      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1ke 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1ke()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ***
      real*8 ave

      
      call cke_ave(ispin,neq,psi1_tag,ave)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_psi_vl_ave
      external c_electron_psi_vl_ave
 
      cpsi_1vl = c_electron_psi_vl_ave(psi1_tag)

      return
      end

*     ***********************************
*     *					*
*     *		cpsi_1vnl 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1vnl()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_psi_vnl_ave
      external c_electron_psi_vnl_ave
 
      cpsi_1vnl = c_electron_psi_vnl_ave(psi1_tag)
      return
      end


*     ***********************************
*     *					*
*     *		cpsi_1vnlso 		*
*     *					*
*     ***********************************
      real*8 function cpsi_1vnlso()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_psi_vnlso_ave
      external c_electron_psi_vnlso_ave
 
      cpsi_1vnlso = c_electron_psi_vnlso_ave(psi1_tag)
      return
      end


*     *******************************
*     *				    *
*     *		c_rho_1exc	    *
*     *				    *
*     *******************************
      real*8 function c_rho_1exc()
      implicit none

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      call c_electron_run(psi2_tag,
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
      subroutine cpsi_1get_Tgradient(THpsi_tag,E0)
      implicit none
      integer    THpsi_tag
      real*8     E0

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      integer tmp1_tag
 
*     **** external functions ****
      integer  cpsi_data_push_stack,Pneb_w_size
      real*8   c_electron_energy
      external cpsi_data_push_stack,Pneb_w_size
      external c_electron_energy

      tmp1_tag = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,1))
c      if (.not.Pneb_w_push_get(0,0,tmp1))
c     >  call errquit("psi_1_get_Tgradient:out of stack memory",0,MA_ERR)


      call c_electron_run(psi1_tag,
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))
      E0 = c_electron_energy(psi1_tag,
     >                       dbl_mb(rho1(1)),
     >                       dcpl_mb(dng1(1)),
     >                       dbl_mb(rho1_all(1)))
      call c_electron_gen_hml(psi1_tag,
     >                        tmp1_tag)

      call c_electron_get_Tgradient(psi1_tag,
     >                             tmp1_tag,
     >                             THpsi_tag)
      
      call cpsi_data_pop_stack(tmp1_tag)
c      if (.not.Pneb_w_pop_stack(tmp1))
c     > call errquit("psi_1_get_Tgradient:error popping stack",0,MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		cpsi_2get_Tgradient	*
*     *					*
*     ***********************************
      subroutine cpsi_2get_Tgradient(option,THpsi_tag,E0)
      implicit none
      integer    option
      integer    THpsi_tag
      real*8     E0

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"
#include "errquit.fh"

*     *** local variables ****
      integer tmp1_tag

*     **** external functions ****
      integer  cpsi_data_push_stack,Pneb_w_size
      real*8   c_electron_energy
      external cpsi_data_push_stack,Pneb_w_size
      external c_electron_energy


      tmp1_tag = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,1))
c      value = BA_push_get(mt_dcpl,(2*ne(1)*ne(1)*nbrillioun),
c     >                   'tmp1',tmp1(2),tmp1(1))
c      if (.not. value) 
c     >   call errquit('out of stack memory in psi_1get_Tradient',0,
c     &       MA_ERR)


      if (option.le.1) then
      call c_electron_run(psi2_tag,
     >                   dbl_mb(rho2(1)),
     >                   dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      end if

      E0 = c_electron_energy(psi2_tag,
     >                       dbl_mb(rho2(1)),
     >                       dcpl_mb(dng2(1)),
     >                       dbl_mb(rho2_all(1)))



      call c_electron_gen_hml(psi2_tag,
     >                       tmp1_tag)
      call c_electron_get_Tgradient(psi2_tag,
     >                             tmp1_tag,
     >                             THpsi_tag)
      
      call cpsi_data_pop_stack(tmp1_tag)
c      value = BA_pop_stack(tmp1(2))
c      if (.not. value) 
c     > call errquit('error popping stack memory in psi_1get_Tradient',0,
c     &       MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *	    cpsi_1geodesic_transport	*
*     *					*
*     ***********************************
      subroutine cpsi_1geodesic_transport(t,H0_tag)
      implicit none
      real*8 t
      integer H0_tag

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_transport(t,psi1_tag,H0_tag)

      return
      end


*     ***********************************
*     *					*
*     *	    cpsi_1geodesic_Gtransport	*
*     *					*
*     ***********************************
      subroutine cpsi_1geodesic_Gtransport(t,G0_tag)
      implicit none
      real*8 t
      integer G0_tag

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_Gtransport(t,psi1_tag,G0_tag)

      return
      end


*     ****************************************
*     *                                      *
*     *     cpsi_1geodesic_transport_junk    *
*     *                                      *
*     ****************************************
      subroutine cpsi_1geodesic_transport_junk(t,H0)
      implicit none
      real*8 t
      complex*16 H0(*)

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_transport_junk(t,psi1_tag,H0)

      return
      end



*     ****************************************
*     *                                      *
*     *     cpsi_1geodesic_Gtransport_junk   *
*     *                                      *
*     ****************************************
      subroutine cpsi_1geodesic_Gtransport_junk(t,G0)
      implicit none
      real*8 t
      complex*16 G0(*)

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_Gtransport_junk(t,psi1_tag,G0)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

      real*8 e_new

*     **** external functions ****
      real*8   c_electron_energy
      external c_electron_energy

      call c_geodesic_get(t,psi1_tag,
     >                      psi2_tag)


      call c_electron_run(psi2_tag,
     >                    dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                    dbl_mb(rho2_all(1)))


      e_new =  c_electron_energy(psi2_tag,
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** external functions ****
      real*8   c_electron_eorbit
      external c_electron_eorbit

      call c_geodesic_transport(t,psi1_tag,psi2_tag)
      cpsi_geodesic_denergy = 2.0d0*c_electron_eorbit(psi2_tag)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call c_geodesic_get(t,psi1_tag,psi2_tag)
      return
      end

*     ***********************************
*     *					*
*     *		cpsi1to2_sd_update	*
*     *					*
*     ***********************************
      subroutine cpsi1to2_sd_update(dte)
      implicit none
      real*8 dte

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      logical value
      integer ms,nb,ierr
      integer tmp_L(2)
      integer psi1_shift,psi2_shift
      real*8 tsum
      complex*16 tttc(16,16)
      integer i,j

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_push_get_block,Pneb_w_pop_stack
      integer  cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_push_get_block,Pneb_w_pop_stack
      external cpsi_data_get_chnk

      call c_electron_run(psi1_tag,
     >                     dbl_mb(rho1(1)),
     >                    dcpl_mb(dng1(1)),
     >                     dbl_mb(rho1_all(1)))

*     **** do a steepest descent step ****
      call c_electron_sd_update(psi1_tag,psi2_tag,dte)

*     **** lagrange multiplier corrections ****

*     **** allocate MA local variables ****
      if (.not.Pneb_w_push_get_block(0,1,8,tmp_L))
     >  call errquit('psi1to2_sd_update: stack failure',0,MA_ERR)

      do nb=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
        psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
        do ms=1,ispinq
        if (neq(ms).gt.0) then
           !psi1_shift = cpsi_data_get_ptr(psi1_tag,nb,1+(ms-1)*neq(1))
           !psi2_shift = cpsi_data_get_ptr(psi2_tag,nb,1+(ms-1)*neq(1))
           call cpsi_lmbda(ms,nb,npack1,
     >                 dbl_mb(psi1_shift),
     >                 dbl_mb(psi2_shift),
     >                 dte,
     >                 dcpl_mb(tmp_L(1)),ierr)

        end if
        end do
      end do

      if (.not.Pneb_w_pop_stack(tmp_L))
     >  call errquit('psi1to2_sd_update: stack failure',0,MA_ERR)

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

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      integer ms,nb,ierr
      integer tmp_L(2)
      integer psi1_shift,psi2_shift

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_push_get_block,Pneb_w_pop_stack
      integer  cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_push_get_block,Pneb_w_pop_stack
      external cpsi_data_get_chnk


      call c_electron_run_noscf(psi1_tag,
     >                          dbl_mb(rho1(1)),
     >                          dcpl_mb(dng1(1)),
     >                          dbl_mb(rho1_all(1)))

*     **** do a steepest descent step ****
      call c_electron_sd_update(psi1_tag,psi2_tag,dte)

*     **** lagrange multiplier corrections ****

*     **** allocate MA local variables ****
      if (.not.Pneb_w_push_get_block(0,1,8,tmp_L))
     >     call errquit(
     >          'psi1to2_sd_update_noscf: stack failure',0,MA_ERR)

      do nb=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
        psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
        do ms=1,ispinq
         !psi1_shift = cpsi_data_get_ptr(psi1_tag,nb,1+(ms-1)*neq(1))
         !psi2_shift = cpsi_data_get_ptr(psi2_tag,nb,1+(ms-1)*neq(1))
         call cpsi_lmbda(ms,nb,npack1,
     >                 dbl_mb(psi1_shift),
     >                 dbl_mb(psi2_shift),
     >                 dte,
     >                 dcpl_mb(tmp_L(1)),ierr)
      end do
      end do

      if (.not.Pneb_w_pop_stack(tmp_L))
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call cpsp_f_vlocal(dcpl_mb(dng1(1)),fion)
      call cpsp_f_vnonlocal(ispin,ne,psi1_tag,fion)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call cke_euv(ispin,neq,psi1_tag,stress)
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
      logical  control_only_lda
      external control_only_lda

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
      !if ((gga.ge.10).and.(gga.lt.100)) then
      if (.not.control_only_lda()) then
         call D3dB_n2ft3d(1,n2ft3d)
         if (.not.
     >      BA_push_get(mt_dbl,(ispin*n2ft3d),'trho',trho(2),trho(1)))
     >      call errquit(
     >      'c_rho_1exc_stress: out of stack memory',1,MA_ERR)

         do ms=1,ispin
            call C3dB_D3dB_r_Copy(1,dbl_mb(rho1_all(1)+(ms-1)*nfft3d),
     >                              dbl_mb(trho(1) +(ms-1)*n2ft3d))
         end do
         call v_bwexc_euv(gga,n2ft3d,ispin,dbl_mb(trho(1)),
     >                  1.0d0,1.0d0,tstress)
         do v=1,3
         do u=1,3
          stress(u,v) = stress(u,v) + tstress(u,v)
         end do
         end do

         if (.not.BA_pop_stack(trho(2)))
     >      call errquit(
     >     'c_rho_1exc_stress: error popping stack memory',1,MA_ERR)
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

      call cpsp_v_nonlocal_euv_2(ispin,neq,psi1_tag,stress)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      call c_electron_gen_hml(psi1_tag,hml_tag)

      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_scal_hml           *
*     *                                 *
*     ***********************************
      subroutine cpsi_scal_hml(alpha)
      implicit none
      real*8 alpha

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer nb,hml_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk
     
      do nb=1,nbrillq
         hml_shift = cpsi_data_get_chnk(hml_tag,nb)
         call Pneb_w_scal(0,nb,alpha,dbl_mb(hml_shift))
      end do
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,eig_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      eig_shift = cpsi_data_get_chnk(eig_tag,nb)
      if (spin_orbit) then
        indx=eig_shift+i-1
      else
        indx=eig_shift+i-1+(ms-1)*ne(1)
      end if
      cpsi_eigenvalue = dbl_mb(indx) 
      return
      end



*     ***********************************
*     *                                 *
*     *         cpsi_virtual            *
*     *                                 *
*     ***********************************
      real*8  function cpsi_virtual(nb,ms,i)
      implicit none
      integer nb,ms
      integer i

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,eig_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      eig_shift = cpsi_data_get_chnk(eig_excited_tag,nb)
      if (spin_orbit) then
        indx=eig_shift+i-1
      else
        indx=eig_shift+i-1+(ms-1)*ne(1)
      end if
      cpsi_virtual = dbl_mb(indx)
      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_occupation         *
*     *                                 *
*     ***********************************
      real*8  function cpsi_occupation(nb,ms,i)
      implicit none
      integer nb,ms
      integer i

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,occ_shift,occ_tag
      real*8  occupation

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_get_next
      external cpsi_data_get_chnk,cpsi_data_get_next

      occ_tag   = cpsi_data_get_next(psi1_tag)
      if (occ_tag.ge.0) then
         occ_shift = cpsi_data_get_chnk(occ_tag,nb)
         if (spin_orbit) then
        indx=occ_shift+i-1
         else
        indx=occ_shift+i-1+(ms-1)*ne(1)
         end if
         occupation = dbl_mb(indx)
      else
         occupation = 1.0d0
      end if
      cpsi_occupation = occupation
      return
      end


*     ***********************************
*     *                                 *
*     *      cpsi_1define_occupation    *
*     *                                 *
*     ***********************************
      subroutine cpsi_1define_occupation(initial_alpha,use_hml)
      implicit none
      real*8  initial_alpha
      logical use_hml

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer it,itmax
      parameter (itmax=50)

      integer ms,nb,n,shift1,shift2,occ1_tag,ndiff
      real*8 e,x,kT,f,g,wb,alpha,pi
      real*8 ZZ,Z(2),Zlower,Zmid,Zupper,elower,emid,eupper
      real*8 flower,fmid,fupper

*     **** external functions ****
      integer  cpsi_data_get_next,cpsi_data_get_chnk
      integer  control_multiplicity
      real*8   brillioun_weight
      real*8   control_TotalCharge,ion_TotalCharge_qm
      real*8   cpsi_occ_distribution,control_fractional_alpha
      external cpsi_data_get_next,cpsi_data_get_chnk
      external control_multiplicity
      external brillioun_weight
      external control_TotalCharge,ion_TotalCharge_qm
      external cpsi_occ_distribution,control_fractional_alpha

      if (use_hml) then
         do ms=1,ispin
            do nb=1,nbrillq
               shift1 = cpsi_data_get_chnk(eig_tag, nb) + (ms-1)*ne(1)
               shift2 = cpsi_data_get_chnk(hml_tag,nb)
     >                + (ms-1)*ne(1)*ne(1)
               do n=1,ne(ms)
                  dbl_mb(shift1+n-1) = dbl_mb(shift2+(n-1)+(n-1)*ne(ms))
               end do
            end do
         end do
      end if

      smearfermi(1)      = 0.0d0
      smearfermi(2)      = 0.0d0
      smearcorrection = 0.0d0
      occ1_tag = cpsi_data_get_next(psi1_tag)
      if (occ1_tag.gt.0) then
         if (initial_alpha.lt.0.0d0) then
            alpha = control_fractional_alpha()
         else
            alpha = initial_alpha
         end if
         kT    = smearkT
         ZZ  = ion_TotalCharge_qm() - control_TotalCharge()
         if (ispin.eq.2) then
            ndiff = control_multiplicity() - 1
            Z(1) = 0.5d0*(ZZ+ndiff)
            Z(2) = 0.5d0*(ZZ-ndiff)
         else
            Z(1) = 0.5d0*ZZ
            Z(2) = 0.0d0
         end if
       
         pi    = 4.0d0*datan(1.0d0)
         !if (initial) alpha = 1.0d0
         if (smeartype.le.0) alpha = 0.0d0

*        **** outer loop over spins ****
         smearcorrection = 0.0d0
         do ms=1,ispin

*           **** find eupper and elower ****
            elower =  9.9d12
            eupper = -9.9d12
            do nb=1,nbrillq
               shift1 = cpsi_data_get_chnk(eig_tag, nb) + (ms-1)*ne(1)
               do n=1,ne(ms)
                 e       = dbl_mb(shift1)
                 if (e.lt.elower) elower = e
                 if (e.gt.eupper) eupper = e
                 shift1  = shift1 + 1
               end do
            end do
            call K1dB_MinAll(elower)
            call K1dB_MaxAll(eupper)

*           **** find fermi level ****
  10        Zlower = 0.0d0
            Zupper = 0.0d0
            do nb=1,nbrillq
               shift1 = cpsi_data_get_chnk(eig_tag, nb) + (ms-1)*ne(1)
               wb = brillioun_weight(nb)
               do n=1,ne(ms)
                 e = dbl_mb(shift1)
                 Zlower = Zlower 
     >            + wb*cpsi_occ_distribution(smeartype,(e-elower)/kT)
                 Zupper = Zupper 
     >            + wb*cpsi_occ_distribution(smeartype,(e-eupper)/kT)
                 shift1  = shift1 + 1
               end do
            end do
            call K1dB_SumAll(Zlower)
            call K1dB_SumAll(Zupper)
            flower = Zlower - Z(ms)
            fupper = Zupper - Z(ms)

            if (flower*fupper.ge.0.0d0) 
     >       call errquit(
     >            'cpsi_1define_occupation:Fermi energy not found',ms,0)

            it = 0
  20        it = it + 1
            emid = 0.5d0*(elower + eupper)
            Zmid = 0.0d0
            do nb=1,nbrillq
               shift1 = cpsi_data_get_chnk(eig_tag, nb) + (ms-1)*ne(1)
               wb = brillioun_weight(nb)
               do n=1,ne(ms)
                 e = dbl_mb(shift1)
                 Zmid = Zmid 
     >                + wb*cpsi_occ_distribution(smeartype,(e-emid)/kT)
                 shift1  = shift1 + 1
               end do
            end do
            call K1dB_SumAll(Zmid)
            fmid = Zmid - Z(ms)
            if (fmid.lt.0.0d0) then
               flower = fmid
               elower = emid
            else 
               fupper = fmid
               eupper = emid
            end if
            if ( (dabs(fmid)     .gt.1.0d-11) .and.
     >           ((eupper-elower).gt.1.0d-11) .and.
     >           (it.lt.itmax))   goto 20
   
            smearfermi(ms) = emid

*           **** determine filling and correction ****
            do nb=1,nbrillq
               shift1 = cpsi_data_get_chnk(eig_tag, nb) + (ms-1)*ne(1)
               shift2 = cpsi_data_get_chnk(occ1_tag,nb) + (ms-1)*ne(1)
               wb = brillioun_weight(nb)
               do n=1,ne(ms)
                 e = dbl_mb(shift1)
                 x = (e - smearfermi(ms))/kT
                 f = cpsi_occ_distribution(smeartype,x)

                 dbl_mb(shift2) = (1.0d0-alpha)*dbl_mb(shift2) + alpha*f

                 if (smeartype.eq.1) then
                    if (  (dbl_mb(shift2)       .gt.1.0d-6) .and. 
     >                  ( (1.0d0-dbl_mb(shift2)).gt.1.0d-6)) then
                   smearcorrection = smearcorrection  
     >              + kT*wb*( dbl_mb(shift2)*log(dbl_mb(shift2)) 
     >              + (1.0d0-dbl_mb(shift2))*log(1.0d0-dbl_mb(shift2)) )
                    end if
                 else if (smeartype.eq.2) then
                   smearcorrection 
     >                 = smearcorrection 
     >                 - kT*wb*dexp(-x*x)/(4.0d0*dsqrt(pi))
                 else if (smeartype.eq.4) then
                   smearcorrection
     >                 = smearcorrection
     >                 - kT*wb*dexp(-(x+dsqrt(0.5d0))*(x+dsqrt(0.5d0)))
     >                 * (1.0d0 + dsqrt(2.0d0) * x)
     >                 / (2.0d0 * dsqrt(pi))
                 end if
                 shift1  = shift1 + 1
                 shift2  = shift2 + 1
               end do
            end do

         end do !** ms***
         if (ms.eq.1) smearcorrection=smearcorrection+smearcorrection
         call K1dB_SumAll(smearcorrection)

      end if

      return
      end


c  set nwpw:fractional_smeartype 1 #0-none, 1-Fermi-Dirac, 2-Gaussian, 3-Hermite
c                                   4-Marzari-Vanderbilt

      real*8 function cpsi_occ_distribution(smeartype,e)
      implicit none
      integer smeartype
      real*8 e
      real*8 f

*     **** external functions ****
      real*8   util_erfc
      external util_erfc

      if (smeartype.eq.1) then
         if (e.gt.30.0d0) then
           f = 0.0d0
         else if (e.lt.(-30.0d0)) then
           f = 1.0d0
         else
           f = 1.0d0/(1.0d0+dexp(e))
         end if
      else if (smeartype.eq.2) then
         f = 0.5d0*util_erfc(e)
      else if (smeartype.eq.4) then
         f = dexp(-(e + dsqrt(0.5d0)) * (e + dsqrt(0.5d0)))
     >     * dsqrt(0.125d0 / datan(1.0d0))
     >     + 0.5d0 * util_erfc(e + dsqrt(0.5d0))
      else 
         if (e.gt.0.0d0) then
           f = 0.0d0
         else
           f = 1.0d0
         end if
      end if
      cpsi_occ_distribution = f 
      return
      end

      real*8 function cpsi_smearfermi(ms)
      implicit none
      integer ms
#include "cpsi_common.fh"
      cpsi_smearfermi = smearfermi(ms)
      return
      end
      real*8 function cpsi_smearcorrection()
      implicit none
#include "cpsi_common.fh"
      cpsi_smearcorrection = smearcorrection
      return
      end





*     ***********************************
*     *                                 *
*     *         cpsi_1gen_svector       *
*     *                                 *
*     ***********************************
      subroutine cpsi_1gen_svector()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      logical value
      integer nb,psi1_shift,svec_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      if (spin_orbit) then
         call Pneb_f_SOSpins_tag(psi1_tag,svec_tag)
c         do nb=1,nbrillq
c            psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
c            svec_shift = cpsi_data_get_chnk(svec_tag,nb)
c            call BGrsm_f_Spins(nb,dbl_mb(psi1_shift),        
c     >                            dbl_mb(svec_shift))
c         end do
      else
         call cpsi_data_update(svec_tag)
         do nb=1,nbrillq
            svec_shift = cpsi_data_get_chnk(svec_tag,nb)
            call dcopy(3*ne(1),0.0d0,0,dbl_mb(svec_shift),1)
         end do
         call cpsi_data_noupdate(svec_tag)
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
      integer nb
      integer i,xyz

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,svec_shift
      real*8  w

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      w=0.0d0
      if (spin_orbit) then
         svec_shift = cpsi_data_get_chnk(svec_tag,nb)
         indx=svec_shift + (xyz-1) 
     >                   + 3*(i-1)
         w = dbl_mb(indx)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer hml_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk
      complex*16 Pneb_w_value
      external   Pneb_w_value

      hml_shift = cpsi_data_get_chnk(hml_tag,nb)
      cpsi_hml=Pneb_w_value(ms,1,i,j,dbl_mb(hml_shift))
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

#include "bafdecls.fh"
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
    
#include "cpsi_common.fh"

      call Calculate_cpsi_spin2(nbrillq,ispin,ne,npack1,psi1_tag,Sab)
      return
      end


*     ***************************
*     *				*
*     *	      cpsi_1rotate2     *
*     *				*
*     ***************************
      subroutine cpsi_1rotate2()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables *****
      integer nb,hml_shift,psi1_shift,psi2_shift
      complex*16 zero,one

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)

      call cpsi_data_update(psi2_tag)
      do nb=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
        psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
        hml_shift  = cpsi_data_get_chnk(hml_tag,nb)
        call Pneb_fwf_Multiply(0,nb,
     >                       one,
     >                       dbl_mb(psi1_shift),npack1,
     >                       dbl_mb(hml_shift),
     >                       zero,
     >                       dbl_mb(psi2_shift))

      end do
      call cpsi_data_noupdate(psi2_tag)

      return
      end


*     ***********************************
*     *					*
*     *		cpsi_diagonalize_hml	*
*     *					*
*     ***********************************
      subroutine cpsi_diagonalize_hml()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      integer nb,hml_shift,eig_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call cpsi_data_update(hml_tag)
      call cpsi_data_update(eig_tag)
      do nb=1,nbrillq
        hml_shift = cpsi_data_get_chnk(hml_tag,nb)
        eig_shift = cpsi_data_get_chnk(eig_tag,nb)
        call Pneb_w_diag(0,nb,dbl_mb(eig_shift),dbl_mb(hml_shift))
      end do
      call cpsi_data_noupdate(hml_tag)
      call cpsi_data_noupdate(eig_tag)

      return
      end



*     ***************************
*     *				*
*     *		cpsi_error	*
*     *				*
*     ***************************
      real*8 function cpsi_error()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      integer n,nb
      integer psi1_shift,psi2_shift,neall
      real*8  error,sum,size
      integer tmp1(2)

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      if (.not.BA_push_get(mt_dcpl,(npack1),'tmp1',tmp1(2),tmp1(1)))
     >   call errquit('out of stack memory in psi_error',0,MA_ERR)

      error = 0.0d0
      neall = neq(1)+neq(2)
      size =  dble(ne(1)+ne(2))
      do nb=1,nbrillq
      do n=1,neall
         psi1_shift = cpsi_data_get_ptr(psi1_tag,nb,n)
         psi2_shift = cpsi_data_get_ptr(psi2_tag,nb,n)
         call Cram_cc_Sub(nb,dbl_mb(psi2_shift),
     >                       dbl_mb(psi1_shift),
     >                       dcpl_mb(tmp1(1)))
         call Cram_cc_dot(nb,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)
         error = error + sum
      end do
      end do
      error = dsqrt(error)/size
      
      if (.not.BA_pop_stack(tmp1(2)))
     > call errquit('error popping stack memory in psi_error',0,MA_ERR)

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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     ***** local variables ****
      logical value
      integer k,nx,ny,nz
      real*8  error,scale
      integer tmp1(2)

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega

      value = BA_push_get(mt_dbl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
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

      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      cpsi_nbrill = nbrillioun
      return
      end



*     ***************************
*     *                         *
*     *         cpsi_nbrillq    *
*     *                         *
*     ***************************
      integer function cpsi_nbrillq()
      implicit none

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "cpsi_common.fh"

      cpsi_neq = neq(ms)
      return
      end


*     ******************************
*     *                            *
*     *         cpsi_ne_excited    *
*     *                            *
*     ******************************
      integer function cpsi_ne_excited(ms)
      implicit none
      integer ms

#include "bafdecls.fh"
#include "cpsi_common.fh"

      cpsi_ne_excited = ne_excited(ms)
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,psi_nogrid,expand_psi
      integer neall,n,vers
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*50 filename
      character*50 control_input_psi
      external     control_input_psi
      logical  c_wvfnc_expander,control_spin_orbit
      external c_wvfnc_expander,control_spin_orbit
  
*     **** external functions ****
      logical  pspw_reformat_c_wvfnc
      integer  brillioun_nbrillioun,psi_get_version
      integer  Pneb_ispinq,Pneb_nbrillq,Pneb_w_size
      integer  cpsi_data_alloc,control_fractional_smeartype
      double precision util_random,control_fractional_kT
      external pspw_reformat_c_wvfnc
      external brillioun_nbrillioun,psi_get_version
      external Pneb_ispinq,Pneb_nbrillq,Pneb_w_size
      external cpsi_data_alloc,control_fractional_smeartype
      external util_random,control_fractional_kT

      integer  cpsi_data_get_ptr,l
      external cpsi_data_get_ptr
      real*8   w

      real*8   cpsi_1ke
      external cpsi_1ke



      expand_psi = .false.

*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)  !** only ne and occupation are used 
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
        expand_psi = .true.
      end if

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)

      neall = neq(1)+neq(2)
      n     = ne(1)
       
*     **** allocate memory ****

c      value = BA_alloc_get(mt_dcpl,npack1*neall*nbrillq,
c     >                     'psi2',psi2(2),psi2(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,npack1*neall*nbrillq,
c     >                     'psi1',psi1(2),psi1(1))
      psi2_tag = cpsi_data_alloc(nbrillq,neall,2*npack1)
      psi1_tag = cpsi_data_alloc(nbrillq,neall,2*npack1)

*     *** fractional orbitals ***
      if (smearoccupation.gt.0) then
        call cpsi_data_set_next(psi2_tag,
     >                 cpsi_data_alloc(nbrillq,neall,1))
        call cpsi_data_set_next(psi1_tag,
     >                 cpsi_data_alloc(nbrillq,neall,1))
        smeartype = control_fractional_smeartype()
        smearkT   = control_fractional_kT()
      end if

      value = BA_alloc_get(mt_dbl,2*nfft3d,'rho1',rho1(2),rho1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho2',rho2(2),rho2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'dng1',dng1(2),dng1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'dng2',dng2(2),dng2(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,(2*n*n*nbrillioun),
c     >                    'hml',hml(2),hml(1))
c      value = value.and.Pneb_w_allocate(0,0,hml)
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,(neall*nbrillq),
c     >                    'eig',eig(2),eig(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,(3*ne(1)*nbrillq),
c     >                    'svec',svec(2),svec(1))
      hml_tag = cpsi_data_alloc(nbrillq,1,2*Pneb_w_size(0,1))
      eig_tag = cpsi_data_alloc(nbrillq,neall,1)
      svec_tag =cpsi_data_alloc(nbrillq,neq(1),3)

      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) 
     >  call errquit('cpsi_initialize:out of heap memory',0, MA_ERR)

c     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nogrid',
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
        value = btdb_parallel(.false.)
        if (taskid.eq.MASTER) then

          filename =  control_input_psi()

          ind = index(filename,' ') - 1
          if (.not. btdb_cput(rtdb,'c_xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_cput(rtdb,'c_xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_put(rtdb,'c_xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_put failed', 0, RTDB_ERR)

          write(*,*)
          write(*,*) "Grid is being converted:"
          write(*,*) "------------------------"
          write(*,*)
          write(*,*) "To turn off automatic grid conversion:"
          write(*,*)
          write(*,*) "set nwpw:psi_nogrid .false."
          write(*,*)
          value = c_wvfnc_expander(rtdb)
          !expand_psi = .true.

        end if
        expand_psi = .true.
        value = btdb_parallel(.true.)

        end if

      end if

*     *****  read initial wavefunctions into psi1  ****
      call cpsi_read(spin_orbit,ispin,ne,nbrillioun,psi1_tag)

*     **** Ortho Check ****
      call Pneb_orthoCheckMake_tag(.true.,0,0,npack1,
     >                             psi1_tag)

      ne_excited(1) = 0
      ne_excited(2) = 0

      !cpsi_initialize = value
      cpsi_initialize = expand_psi
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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      logical value
      integer next1,next2

*     **** external functions ****
      logical  Pneb_w_free
      integer  cpsi_data_get_next
      external Pneb_w_free
      external cpsi_data_get_next

*     ***** write psi1 wavefunctions ****
      if (wpsi) then
        call cpsi_write(spin_orbit,ispin,ne,nbrillioun,psi1_tag)
      end if
    
c      value = BA_free_heap(eig(2))
c      value = value.and.BA_free_heap(svec(2))
c      value = value.and.Pneb_w_free(hml)
c      value = value.and.BA_free_heap(psi2(2))
c      value = value.and.BA_free_heap(psi1(2))
      next1 = cpsi_data_get_next(psi1_tag)
      if (next1.ge.0) call cpsi_data_dealloc(next1)
      call cpsi_data_dealloc(psi1_tag)

      next2 = cpsi_data_get_next(psi2_tag)
      if (next2.ge.0) call cpsi_data_dealloc(next2)
      call cpsi_data_dealloc(psi2_tag)

      call cpsi_data_dealloc(eig_tag)
      call cpsi_data_dealloc(svec_tag)
      call cpsi_data_dealloc(hml_tag)

      value =           BA_free_heap(dng2(2))
      value = value.and.BA_free_heap(dng1(2))
      value = value.and.BA_free_heap(rho2(2))
      value = value.and.BA_free_heap(rho1(2))
      value = value.and.BA_free_heap(rho2_all(2))
      value = value.and.BA_free_heap(rho1_all(2))


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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      logical value

*     **** external functions ****
      logical  Pneb_w_free
      external Pneb_w_free

*     ***** write psi1 wavefunctions ****
      call cpsi_write2(spin_orbit,ispin,ne,nbrillioun,psi1_tag,wpsi)
    
c      value = BA_free_heap(eig(2))
c      value = value.and.Pneb_w_free(hml)
c      value = value.and.BA_free_heap(psi2(2))
c      value = value.and.BA_free_heap(psi1(2))
      call cpsi_data_dealloc(eig_tag)
      call cpsi_data_dealloc(hml_tag)
      call cpsi_data_dealloc(svec_tag)
      call cpsi_data_dealloc(psi1_tag)
      call cpsi_data_dealloc(psi2_tag)
      value =           BA_free_heap(dng2(2))
      value = value.and.BA_free_heap(dng1(2))
      value = value.and.BA_free_heap(rho2(2))
      value = value.and.BA_free_heap(rho1(2))
      value = value.and.BA_free_heap(rho2_all(2))
      value = value.and.BA_free_heap(rho1_all(2))

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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      character*50 filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

*     ***** write psi1 wavefunctions ****
      filename = 'dos'//c_index_name(k)//'.movecs'
      call cpsi_write_filename(filename,spin_orbit,ispin,ne,nbrillioun,
     >                         psi1_tag)

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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


*     **** local variables ****
      character*50 filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

*     ***** write psi1 wavefunctions ****
      filename = 'dos'//c_index_name(k)//'.movecs'
      call cpsi_read_filename(filename,spin_orbit,ispin,ne,nbrillioun,
     >                         psi1_tag)

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
c#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer i,ii
      real*8 ei,eiv,EV

c     **** external functions ****
      integer  brillioun_nbrillioun
      real*8 cpsi_eigenvalue
      real*8   brillioun_k,brillioun_ks,brillioun_weight
      external brillioun_nbrillioun
      external cpsi_eigenvalue
      external brillioun_k,brillioun_ks,brillioun_weight
 
      EV=27.2116d0
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
#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"
      integer i,ii,iend
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



*     ***********************************
*     *					*
*     *		cpsi_eigenvalue_brdcst 	*
*     *					*
*     ***********************************
      real*8  function cpsi_eigenvalue_brdcst(nb,ms,i)
      implicit none
      integer nb,ms
      integer i

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,eig_shift,nbq,taskid_k,p
      real*8  e

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)
      
      e = 0.0d0
      if (p.eq.taskid_k) then
         eig_shift = cpsi_data_get_chnk(eig_tag,nbq)
         if (spin_orbit) then
           indx=eig_shift+i-1
         else
           indx=eig_shift+i-1+(ms-1)*ne(1)
         end if
         e = dbl_mb(indx) 
      end if
      call K1dB_SumAll(e)

      cpsi_eigenvalue_brdcst = e
      return
      end




*     ***********************************
*     *                                 *
*     *         cpsi_virtual_brdcst     *
*     *                                 *
*     ***********************************
      real*8  function cpsi_virtual_brdcst(nb,ms,i)
      implicit none
      integer nb,ms
      integer i

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,eig_shift,nbq,taskid_k,p
      real*8  e

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      e = 0.0d0
      if (p.eq.taskid_k) then
         eig_shift = cpsi_data_get_chnk(eig_excited_tag,nbq)
         if (spin_orbit) then
           indx=eig_shift+i-1
         else
           indx=eig_shift+i-1+(ms-1)*ne_excited(1)
         end if
         e = dbl_mb(indx)
      end if
      call K1dB_SumAll(e)

      cpsi_virtual_brdcst = e
      return
      end



*     ***********************************
*     *                                 *
*     *        cpsi_eig_brdcst_tag      *
*     *                                 *
*     ***********************************
      real*8 function cpsi_eig_brdcst_tag(ne,spin_orbit,eig_tag,nb,ms,i)
      implicit none
      integer ne(2)
      logical spin_orbit
      integer eig_tag
      integer nb,ms
      integer i

#include "bafdecls.fh"

*     **** local variables ****
      integer indx,eig_shift,nbq,taskid_k,p
      real*8  e

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      e = 0.0d0
      if (p.eq.taskid_k) then
         eig_shift = cpsi_data_get_chnk(eig_tag,nbq)
         if (spin_orbit) then
           indx=eig_shift+i-1
         else
           indx=eig_shift+i-1+(ms-1)*ne(1)
         end if
         e = dbl_mb(indx)
      end if
      call K1dB_SumAll(e)

      cpsi_eig_brdcst_tag = e
      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_occupation_brdcst  *
*     *                                 *
*     ***********************************
      real*8  function cpsi_occupation_brdcst(nb,ms,i)
      implicit none
      integer nb,ms
      integer i

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,occ_shift,occ_tag,nbq,taskid_k,p
      real*8  occupation

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_get_next
      external cpsi_data_get_chnk,cpsi_data_get_next

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      occupation = 0.0d0
      if (p.eq.taskid_k) then
         occ_tag   = cpsi_data_get_next(psi1_tag)
         if (occ_tag.ge.0) then
            occ_shift = cpsi_data_get_chnk(occ_tag,nbq)
            if (spin_orbit) then
               indx=occ_shift+i-1
            else
               indx=occ_shift+i-1+(ms-1)*ne(1)
            end if
            occupation = dbl_mb(indx)
         else
            occupation = 1.0d0
         end if
      end if
      call K1dB_SumAll(occupation)

      cpsi_occupation_brdcst = occupation
      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_occ_brdcst_tag     *
*     *                                 *
*     ***********************************
      real*8 function cpsi_occ_brdcst_tag(ne,spin_orbit,occ_tag,nb,ms,i)
      implicit none
      integer ne(2)
      logical spin_orbit
      integer occ_tag
      integer nb,ms
      integer i

#include "bafdecls.fh"

*     **** local variables ****
      integer indx,occ_shift,nbq,taskid_k,p
      real*8  occupation

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_get_next
      external cpsi_data_get_chnk,cpsi_data_get_next

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      occupation = 0.0d0
      if (p.eq.taskid_k) then
         if (occ_tag.ge.0) then
            occ_shift = cpsi_data_get_chnk(occ_tag,nbq)
            if (spin_orbit) then
               indx=occ_shift+i-1
            else
               indx=occ_shift+i-1+(ms-1)*ne(1)
            end if
            occupation = dbl_mb(indx)
         else
            occupation = 1.0d0
         end if
      end if
      call K1dB_SumAll(occupation)

      cpsi_occ_brdcst_tag = occupation
      return
      end


*     ***********************************
*     *                                 *
*     *         cpsi_svector_brdcst     *
*     *                                 *
*     ***********************************
      real*8  function cpsi_svector_brdcst(nb,i,xyz)
      implicit none
      integer nb
      integer i,xyz

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer indx,svec_shift,nbq,taskid_k,p
      real*8  w

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      w=0.0d0
      if (p.eq.taskid_k) then
         if (spin_orbit) then
            svec_shift = cpsi_data_get_chnk(svec_tag,nbq)
            indx=svec_shift + (xyz-1) + 3*(i-1)
            w = dbl_mb(indx)
          end if
      end if
      call K1dB_SumAll(w)

      cpsi_svector_brdcst = w
      return
      end



*     ***********************************
*     *                                 *
*     *         cpsi_sv_brdcst_tag     *
*     *                                 *
*     ***********************************
      real*8  function cpsi_sv_brdcst_tag(spin_orbit,svec_tag,nb,i,xyz)
      implicit none
      logical spin_orbit
      integer svec_tag
      integer nb
      integer i,xyz

#include "bafdecls.fh"

*     **** local variables ****
      integer indx,svec_shift,nbq,taskid_k,p
      real*8  w

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      w=0.0d0
      if (p.eq.taskid_k) then
         if (spin_orbit) then
            svec_shift = cpsi_data_get_chnk(svec_tag,nbq)
            indx=svec_shift + (xyz-1) + 3*(i-1)
            w = dbl_mb(indx)
          end if
      end if
      call K1dB_SumAll(w)

      cpsi_sv_brdcst_tag = w
      return
      end




*     ***************************
*     *				*
*     *	    ecpsi_initialize 	*
*     *				*
*     ***************************

      logical function ecpsi_initialize()
      implicit none 
      logical ortho

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer neall,n,vers,ispin0
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*50 filename
      character*50 control_input_epsi
      external     control_input_epsi
      logical  c_wvfnc_expander
      external c_wvfnc_expander
  
*     **** external functions ****
      logical  pspw_reformat_c_wvfnc
      integer  brillioun_nbrillioun,psi_get_version
      integer  Pneb_ispinq,Pneb_nbrillq,Pneb_w_size
      integer  cpsi_data_alloc
      double precision util_random,control_fractional_kT
      external pspw_reformat_c_wvfnc
      external brillioun_nbrillioun,psi_get_version
      external cpsi_data_alloc

      integer  cpsi_data_get_ptr,l
      external cpsi_data_get_ptr
      real*8   w

      value =.true.

*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne_excited(ispin0,ne_excited)
      neall  = ne_excited(1)  + ne_excited(2)

       
*     **** allocate memory ****
      psi2_excited_tag = cpsi_data_alloc(nbrillq,neall,2*npack1)
      psi1_excited_tag = cpsi_data_alloc(nbrillq,neall,2*npack1)
      eig_excited_tag  = cpsi_data_alloc(nbrillq,neall,1)


c     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nogrid',
     >                  mt_log,1,psi_nogrid))
     >   psi_nogrid = .true.

      if (psi_nogrid) then

        filename = control_input_epsi()
        call psi_get_header_filename(filename,
     >                      hversion,hnfft,hunita,hispin,hne)

        if ( (hnfft(1).ne.control_ngrid(1)) .or.
     >       (hnfft(2).ne.control_ngrid(2)) .or.
     >       (hnfft(3).ne.control_ngrid(3)) ) then

        hnfft(1) = control_ngrid(1)
        hnfft(2) = control_ngrid(2)
        hnfft(3) = control_ngrid(3)
        call Parallel_taskid(taskid)
        value = btdb_parallel(.false.)
        if (taskid.eq.MASTER) then

          filename =  control_input_epsi()

          ind = index(filename,' ') - 1
          if (.not. btdb_cput(rtdb,'c_xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_cput(rtdb,'c_xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_put(rtdb,'c_xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_put failed', 0, RTDB_ERR)

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
        value = btdb_parallel(.true.)

        end if

      end if

*     *****  read initial wavefunctions into psi1  ****
      call ecpsi_read(spin_orbit,ispin,ne_excited,
     >                nbrillioun,psi1_excited_tag)


      ecpsi_initialize = value
      return
      end



*     ***************************
*     *				*
*     *	     ecpsi_finalize	*
*     *				*
*     ***************************

      logical function ecpsi_finalize(wpsi)
      implicit none 
      logical wpsi

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value
      integer next1,next2

*     **** external functions ****
      integer  cpsi_data_get_next
      external cpsi_data_get_next

*     ***** write psi1 wavefunctions ****
      if (wpsi) then
        call ecpsi_write(spin_orbit,ispin,ne_excited,
     >                   nbrillioun,psi1_excited_tag)
      end if
    
      next1 = cpsi_data_get_next(psi1_excited_tag)
      if (next1.ge.0) call cpsi_data_dealloc(next1)
      call cpsi_data_dealloc(psi1_excited_tag)

      next2 = cpsi_data_get_next(psi2_excited_tag)
      if (next2.ge.0) call cpsi_data_dealloc(next2)
      call cpsi_data_dealloc(psi2_excited_tag)

      call cpsi_data_dealloc(eig_excited_tag)

      value = .true.
      ecpsi_finalize = value
      return
      end

 
*     ***************************
*     *				*
*     *	  ecpsi_generate_kb_H   *
*     *				*
*     ***************************
      subroutine ecpsi_generate_kb_H(dk,nb,H)
      implicit none
      real*8     dk(3)
      integer    nb
      complex*16 H(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value 
      integer ms,i,j,neall(2),isize
      integer nbq,p,taskid_k,indx,indxt
      integer Gx(2),Gy(2),Gz(2),tmp(2),psii,psij
      real*8  kx,ky,kz
      complex*16 vv(3),ww

*     **** external functions ****
      integer  c_G_indx,cpsi_data_get_ptr
      real*8   brillioun_k
      external c_G_indx,cpsi_data_get_ptr
      external brillioun_k

      neall(1) = ne(1)+ne_excited(1)
      neall(2) = ne(2)+ne_excited(2)

      isize = neall(1)*neall(1) + neall(2)*neall(2)
      call dcopy(2*isize,0.0d0,0,H,1)

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)
      if (p.eq.taskid_k) then
         value = BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
         value = value.and.
     >           BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
         value = value.and.
     >           BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
         value = value.and.
     >           BA_push_get(mt_dbl, 2*npack1,'tmp',tmp(2),tmp(1))
         if (.not. value) 
     >      call errquit('ecpsi_generate_kb_H:pushing stack',1,MA_ERR)

*        **** define Gx,Gy and Gz in packed space ****
         kx = brillioun_k(1,nbq)
         ky = brillioun_k(2,nbq)
         kz = brillioun_k(3,nbq)
         call C3dB_t_Copy(1,dbl_mb(c_G_indx(1)),dbl_mb(Gx(1)))
         call C3dB_t_Copy(1,dbl_mb(c_G_indx(2)),dbl_mb(Gy(1)))
         call C3dB_t_Copy(1,dbl_mb(c_G_indx(3)),dbl_mb(Gz(1)))
         do i=1,nfft3d
            dbl_mb(Gx(1)+i-1) = dbl_mb(Gx(1)+i-1) + kx
            dbl_mb(Gy(1)+i-1) = dbl_mb(Gy(1)+i-1) + ky
            dbl_mb(Gz(1)+i-1) = dbl_mb(Gz(1)+i-1) + kz
         end do
         call Cram_r_pack(nbq,dbl_mb(Gx(1)))
         call Cram_r_pack(nbq,dbl_mb(Gy(1)))
         call Cram_r_pack(nbq,dbl_mb(Gz(1)))


         do ms=1,ispin
         do j=1,neall(ms)
            if (j.le.ne(ms)) then
               psij = cpsi_data_get_ptr(psi1_tag,nbq,j)
            else
               psij = cpsi_data_get_ptr(psi1_excited_tag,nbq,j-ne(ms))
            end if

            do i=j,neall(ms)
               indx  = i+(j-1)*neall(ms)+(ms-1)*neall(1)**2
               indxt = j+(i-1)*neall(ms)+(ms-1)*neall(1)**2
               if (i.le.ne(ms)) then
                 psii = cpsi_data_get_ptr(psi1_tag,nbq,i)
               else
                 psii = cpsi_data_get_ptr(psi1_excited_tag,nbq,i-ne(ms))
               end if
               call Cram_rc_Mul(nbq,
     >                          dbl_mb(Gx(1)),
     >                          dbl_mb(psij),
     >                          dbl_mb(tmp(1)))
               call Cram_cc_zdot(nbq,
     >                           dbl_mb(psii),
     >                           dbl_mb(tmp(1)),
     >                           vv(1))
               call Cram_rc_Mul(nbq,
     >                          dbl_mb(Gy(1)),
     >                          dbl_mb(psij),
     >                          dbl_mb(tmp(1)))
               call Cram_cc_zdot(nbq,
     >                           dbl_mb(psii),
     >                           dbl_mb(tmp(1)),
     >                           vv(2))
               call Cram_rc_Mul(nbq,
     >                          dbl_mb(Gz(1)),
     >                          dbl_mb(psij),
     >                          dbl_mb(tmp(1)))
               call Cram_cc_zdot(nbq,
     >                           dbl_mb(psii),
     >                           dbl_mb(tmp(1)),
     >                           vv(3))
               H(indx) = -dk(1)*vv(1) - dk(2)*vv(2) - dk(3)*vv(3)
               if (i.ne.j) H(indxt) = H(indx)
            end do
         end do
         end do
         value =           BA_pop_stack(tmp(2))
         value = value.and.BA_pop_stack(Gz(2))
         value = value.and.BA_pop_stack(Gy(2))
         value = value.and.BA_pop_stack(Gx(2))
         if (.not. value) 
     >      call errquit('ecpsi_generate_kb_H:popping stack',1,MA_ERR)
      end if
      isize = neall(1)*neall(1) + neall(2)*neall(2)
      call K1dB_Vector_SumAll(2*isize,H)

      return
      end




*     ***********************
*     *                     *
*     *     cpsi_Mulliken   *
*     *                     * 
*     ***********************

      subroutine cpsi_Mulliken(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "cpsi_common.fh"


*     **** Lubin Water Analysis ****

*     **** Atom Analysis ****

*     **** Mulliken Analysis ****
      call band_analysis(0,rtdb,ispin,ne,nbrillioun,nbrillq,
     >                   psi1_tag,eig_tag)
      
      call band_gen_APC(ispin,ne,dcpl_mb(dng1(1)))
      call band_print_APC(luout)
      call band_gen_Efield(rtdb,ispin,dcpl_mb(dng1(1)))
     

      return
      end



*     ***********************************
*     *                                 *
*     *     cpsi_projected_dos_weights  *
*     *                                 * 
*     ***********************************

      subroutine cpsi_projected_dos_weights(rtdb,dosgrid,k,pweight,
     >                                      lmax,norbs)
      implicit none
      integer rtdb,dosgrid(3),k
      real*8  pweight(*)
      integer lmax,norbs

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value
      integer l,n,p,nb,nbq,indx,nnmax,taskid_k
      integer weight(2),tmp(2)

      call Parallel3d_taskid_k(taskid_k)

      nnmax = (ne(1)+ne(2))*nbrillq*(4+norbs)
      value = BA_push_get(mt_dbl,(nnmax),'eweight',weight(2),weight(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(ne(1)+ne(2)),'tmp',tmp(2),tmp(1))
      if (.not. value)
     >  call errquit(
     >         'cpsi_projected_dos_weighs:out of stack memory',0,MA_ERR)
      call dcopy(nnmax,0.0d0,0,dbl_mb(weight(1)),1)


*     **** Mulliken Analysis ****
      call band_projecteddos_analysis(0,rtdb,ispin,ne,nbrillq,
     >                                psi1_tag,dbl_mb(weight(1)),
     >                                lmax,norbs)

      do l=0,lmax+norbs
      do nb=1,nbrillioun
         call K1dB_ktoqp(nb,nbq,p)
         if (p.eq.taskid_k) then
            indx = (nbq-1)*(ne(1)+ne(2)) + l*(ne(1)+ne(2))*nbrillq
            call dcopy(ne(1)+ne(2),
     >                 dbl_mb(weight(1)+indx),1,
     >                 dbl_mb(tmp(1)),1)
         else
            call dcopy(ne(1)+ne(2),0.0d0,0,dbl_mb(tmp(1)),1)
         end if
         call K1dB_Vector_SumAll(ne(1)+ne(2),dbl_mb(tmp(1)))

         do n=1,ne(1)+ne(2)
            indx = 1 + (k+nb-2) 
     >         + (n-1)*dosgrid(1)*dosgrid(2)*dosgrid(3)
     >         +     l*dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))
            pweight(indx) = dbl_mb(tmp(1)+n-1)
         end do
      end do
      end do

      value =           BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit(
     >       'cpsi_projected_dos_weights:error freeing stack',0,MA_ERR)
      return
      end


*     ***********************
*     *                     *
*     *    ecpsi_Mulliken   *
*     *                     *
*     ***********************

      subroutine ecpsi_Mulliken(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

      call band_analysis(1,rtdb,ispin,ne_excited,nbrillioun,nbrillq,
     >                   psi1_excited_tag,eig_excited_tag)

      return
      end



*     ***********************
*     *                     *
*     *     cpsi_DOS        *
*     *                     *
*     ***********************

      subroutine cpsi_DOS(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"


*     **** local variables ****
      logical value
      integer npoints,ii,ms,nbq
      integer weight(2),dstates(2),nemax,eig_shift
      real*8 emin,emax,alpha
      character*255 filename

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk
      real*8   brillioun_weight
      external brillioun_weight

      nemax = ne(1)
      value = BA_push_get(mt_dbl,(nemax),'weight',weight(2),weight(1))
      if (.not. value)
     >  call errquit('cpsi_dos:out of stack memory',0, MA_ERR)
      call dcopy(nemax,1.0d0,0,dbl_mb(weight(1)),1)


      if (.not.btdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.btdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if
      value = BA_push_get(mt_dbl,(2*npoints),
     >                    'dstates',dstates(2),dstates(1))
      if (.not. value)
     >  call errquit('cpsi_dos:out of stack memory',0,MA_ERR)


      if (.not.btdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do nbq=1,nbrillq
            eig_shift = cpsi_data_get_chnk(eig_tag,nbq)
            do ii=1,ne(1)+ne(2)
               if (dbl_mb(eig_shift+ii-1).lt.emin) 
     >            emin = dbl_mb(eig_shift+ii-1)
            end do 
         end do
         call K1dB_MinAll(emin)
         emin = emin - 0.1d0
      end if

      if (.not.btdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do nbq=1,nbrillq
            eig_shift = cpsi_data_get_chnk(eig_tag,nbq)
            do ii=1,ne(1)+ne(2)
               if (dbl_mb(eig_shift+ii-1).gt.emax) 
     >            emax = dbl_mb(eig_shift+ii-1)
            end do
         end do
         call K1dB_MaxAll(emin)
         emax = emax + 0.1d0
      end if

*     **** generate DENSITY OF STATES *****
      call dcopy(2*npoints,0.0d0,0,dbl_mb(dstates(1)),1)
      do nbq=1,nbrillq
         eig_shift = cpsi_data_get_chnk(eig_tag,nbq)
         do ms=1,ispin
            call adddensityofstates(
     >                dbl_mb(eig_shift+(ms-1)*ne(1)),
     >                dbl_mb(weight(1)),ne(ms),
     >                (3-2*ms)*brillioun_weight(nbq),
     >                alpha,npoints,emin,emax,
     >                dbl_mb(dstates(1)+(ms-1)*npoints))
         end do
      end do
      call K1dB_Vector_SumAll(2*npoints,dbl_mb(dstates(1)))

      if (ispin.eq.1) then
        filename = "smear_dos_both"
        call writedensityofstates(filename,
     >                     1.0d0,alpha,npoints,emin,emax,
     >                     dbl_mb(dstates(1)))
      end if

      if (ispin.eq.2) then
        filename = "smear_dos_alpha"
        call writedensityofstates(filename,
     >                     1.0d0,alpha,npoints,emin,emax,
     >                     dbl_mb(dstates(1)))
        filename = "smear_dos_beta"
        call writedensityofstates(filename,
     >               -1.0d0,alpha,npoints,emin,emax,
     >                dbl_mb(dstates(1)+npoints))

        filename = "smear_dos_both"
        call daxpy(npoints,-1.0d0,dbl_mb(dstates(1)+npoints),1,
     >                            dbl_mb(dstates(1)),1)
        call writedensityofstates(filename,
     >                     1.0d0,alpha,npoints,emin,emax,
     >                     dbl_mb(dstates(1)))
      end if

      value =           BA_pop_stack(dstates(2))
      value = value.and.BA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit('cpsi_dos: error freeing stack',0, MA_ERR)

      return
      end


*     ***********************
*     *                     *
*     *     ecpsi_DOS       *
*     *                     *
*     ***********************

      subroutine ecpsi_DOS(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"


*     **** local variables ****
      logical value
      integer npoints,ii,ms,nbq
      integer weight(2),dstates(2),nemax,eig_shift
      real*8 emin,emax,alpha
      character*255 filename

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk
      real*8   brillioun_weight
      external brillioun_weight

      nemax = ne_excited(1)
      value = BA_push_get(mt_dbl,(nemax),'weight',weight(2),weight(1))
      if (.not. value)
     >  call errquit('ecpsi_dos:out of stack memory',0, MA_ERR)
      call dcopy(nemax,1.0d0,0,dbl_mb(weight(1)),1)


      if (.not.btdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.btdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if
      value = BA_push_get(mt_dbl,(2*npoints),
     >                    'dstates',dstates(2),dstates(1))
      if (.not. value)
     >  call errquit('ecpsi_dos:out of stack memory',0,MA_ERR)


      if (.not.btdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do nbq=1,nbrillq
            eig_shift = cpsi_data_get_chnk(eig_excited_tag,nbq)
            do ii=1,ne_excited(1)+ne_excited(2)
               if (dbl_mb(eig_shift+ii-1).lt.emin) 
     >            emin = dbl_mb(eig_shift+ii-1)
            end do 
         end do
         call K1dB_MinAll(emin)
         emin = emin - 0.1d0
      end if

      if (.not.btdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do nbq=1,nbrillq
            eig_shift = cpsi_data_get_chnk(eig_excited_tag,nbq)
            do ii=1,ne_excited(1)+ne_excited(2)
               if (dbl_mb(eig_shift+ii-1).gt.emax) 
     >            emax = dbl_mb(eig_shift+ii-1)
            end do
         end do
         call K1dB_MaxAll(emin)
         emax = emax + 0.1d0
      end if

*     **** generate DENSITY OF STATES *****
      call dcopy(2*npoints,0.0d0,0,dbl_mb(dstates(1)),1)
      do nbq=1,nbrillq
         eig_shift = cpsi_data_get_chnk(eig_excited_tag,nbq)
         do ms=1,ispin
            call adddensityofstates(
     >                dbl_mb(eig_shift+(ms-1)*ne_excited(1)),
     >                dbl_mb(weight(1)),ne_excited(ms),
     >                (3-2*ms)*brillioun_weight(nbq),
     >                alpha,npoints,emin,emax,
     >                dbl_mb(dstates(1)+(ms-1)*npoints))
         end do
      end do
      call K1dB_Vector_SumAll(2*npoints,dbl_mb(dstates(1)))

      if (ispin.eq.1) then
        filename = "smear_vdos_both"
        call writedensityofstates(filename,
     >                     1.0d0,alpha,npoints,emin,emax,
     >                     dbl_mb(dstates(1)))
      end if
      if (ispin.eq.2) then
        filename = "smear_vdos_alpha"
        call writedensityofstates(filename,
     >                     1.0d0,alpha,npoints,emin,emax,
     >                     dbl_mb(dstates(1)))
        filename = "smear_vdos_beta"
        call writedensityofstates(filename,
     >               -1.0d0,alpha,npoints,emin,emax,
     >                dbl_mb(dstates(1)+npoints))
        filename = "smear_vdos_both"
        call daxpy(npoints,-1.0d0,dbl_mb(dstates(1)+npoints),1,
     >                            dbl_mb(dstates(1)),1)
        call writedensityofstates(filename,
     >                     1.0d0,alpha,npoints,emin,emax,
     >                     dbl_mb(dstates(1)))
      end if

      value =           BA_pop_stack(dstates(2))
      value = value.and.BA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit('ecpsi_dos: error freeing stack',0, MA_ERR)
      return
      end


