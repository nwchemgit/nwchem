
*     ***************************
*     *				*
*     *		psi_1to2	*
*     *				*
*     ***************************
      subroutine psi_1to2()
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

   
      call zcopy(npack1*(ne(1)+ne(2)),
     >           dcpl_mb(psi1(1)),1,
     >           dcpl_mb(psi2(1)),1)

      return
      end



*     ***************************
*     *				*
*     *		psi_2to1	*
*     *				*
*     ***************************
      subroutine psi_2to1()
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
 
      call zcopy(npack1*(ne(1)+ne(2)),
     >           dcpl_mb(psi2(1)),1,
     >           dcpl_mb(psi1(1)),1)

c      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end

*     ***************************
*     *				*
*     *		psi_check	*
*     *				*
*     ***************************
      subroutine psi_check()
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
 

      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end



*     ***************************
*     *				*
*     *		rho_2to1	*
*     *				*
*     ***************************
      subroutine rho_2to1()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

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
*     *		dng_2to1	*
*     *				*
*     ***************************
      subroutine dng_2to1()
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
 
      call zcopy(npack0,
     >           dcpl_mb(dng2(1)),1,
     >           dcpl_mb(dng1(1)),1)

      return
      end




*     ***********************************
*     *					*
*     *		psi_1toelectron		*
*     *					*
*     ***********************************
      subroutine psi_1toelectron()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all
 
      call electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                  dbl_mb(rho1_all(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_1energy		*
*     *					*
*     ***********************************
      real*8 function psi_1energy()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy

      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))
      psi_1energy = electron_energy(dcpl_mb(psi1(1)),
     >                               dbl_mb(rho1(1)),
     >                              dcpl_mb(dng1(1)),
     >                              dbl_mb(rho1_all(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_2energy		*
*     *					*
*     ***********************************
      real*8 function psi_2energy()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy

      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      psi_2energy = electron_energy(dcpl_mb(psi2(1)),
     >                               dbl_mb(rho2(1)),
     >                              dcpl_mb(dng2(1)),
     >                              dbl_mb(rho2_all(1)))

      return
      end



*     ***********************************
*     *					*
*     *		psi_1eorbit		*
*     *					*
*     ***********************************
      real*8 function psi_1eorbit()
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

*     **** external functions ****
      real*8   electron_eorbit
      external electron_eorbit

      psi_1eorbit = electron_eorbit(dcpl_mb(psi1(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_1ke 		*
*     *					*
*     ***********************************
      real*8 function psi_1ke()
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

*     **** local variables ***
      real*8 ave

*     **** external functions ****
      real*8   electron_eorbit
      external electron_eorbit
      
      call ke_ave(ispin,ne,dcpl_mb(psi1(1)),ave)
      psi_1ke = ave

      return
      end


*     ***********************************
*     *					*
*     *		psi_1vl 		*
*     *					*
*     ***********************************
      real*8 function psi_1vl()
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


*     **** external functions ****
      real*8   electron_psi_vl_ave
      external electron_psi_vl_ave
 
      psi_1vl = electron_psi_vl_ave(dcpl_mb(psi1(1)),dbl_mb(rho1(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_1vnl 		*
*     *					*
*     ***********************************
      real*8 function psi_1vnl()
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


*     **** external functions ****
      real*8   electron_psi_vnl_ave
      external electron_psi_vnl_ave
 
      psi_1vnl = electron_psi_vnl_ave(dcpl_mb(psi1(1)))

      return
      end

*     ***********************************
*     *									*
*     *		rho_1Fmm					*
*     *									*
*     ***********************************
      subroutine rho_1Fmm(F_mm)
      implicit none
      real*8 F_mm(*)

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
      integer n2ft3d,nx,ny,nz
      integer r_grid(2),rho(2)
      real*8  dv

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

*     **** Initializationsr ****
      call D3dB_n2ft3d(1,n2ft3d) 
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** Push memory ****
      value = MA_push_get(mt_dbl,(3*n2ft3d),'r_grid',
     >                       r_grid(2),r_grid(1))
      value = value.and. 
     >        MA_push_get(mt_dbl,(3*n2ft3d),'rho',
     >                       rho(2),rho(1))
      if (.not. value) call errquit(
     >     'rho_1Fmm: out of stack memory',0)


*     **** Get r_grid and rho ****
      call lattice_r_grid(dbl_mb(r_grid(1)))
      call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                   dbl_mb(rho1(1)+(ispin-1)*n2ft3d),
     >                   dbl_mb(rho(1)))

*     **** Now calculate F_mm ****
      call pspw_qmmm_Fmm(n2ft3d,dbl_mb(r_grid(1)),
     >                            dbl_mb(rho(1)),
     >                            dv,F_mm)

*     **** Pop memory ****
      value =           MA_pop_stack(rho(2))
      value = value.and.MA_pop_stack(r_grid(2))
      if (.not. value) call errquit(
     >     'rho_1Fmm: error popping stack memory',0)

      return
      end


*     ***********************************
*     *					*
*     *		rho_1exc		*
*     *					*
*     ***********************************
      real*8 function rho_1exc()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_exc
      external electron_exc

      rho_1exc = electron_exc(dbl_mb(rho1_all(1)))
      return
      end

*     ***********************************
*     *					*
*     *		rho_1pxc		*
*     *					*
*     ***********************************
      real*8 function rho_1pxc()
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

*     **** external functions ****
      real*8   electron_pxc
      external electron_pxc

      rho_1pxc = electron_pxc(dbl_mb(rho1(1)))
      return
      end


*     ***********************************
*     *					*
*     *		dng_1ehartree           *
*     *					*
*     ***********************************
      real*8 function dng_1ehartree()
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

*     **** external functions ****
      integer  control_version
      real*8   electron_ehartree,electron_ehartree2
      external control_version
      external electron_ehartree,electron_ehartree2

*     **** local variables *****
      real*8 eh

      eh = 0.0d0
      if (control_version().eq.3) 
     >    eh = electron_ehartree(dcpl_mb(dng1(1)))

      if (control_version().eq.4) 
     >    eh = electron_ehartree2(dbl_mb(rho1(1)))

      dng_1ehartree = eh
      return
      end



*     ***********************************
*     *					*
*     *		psi_2toelectron		*
*     *					*
*     ***********************************
      subroutine psi_2toelectron()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                   dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_1get_Tgradient	*
*     *					*
*     ***********************************
      subroutine psi_1get_Tgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      logical value
      integer tmp1(2)
 
*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      value = MA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0)


      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))

      Eout =  electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)))


      call electron_gen_hml(dcpl_mb(psi1(1)),
     >                       dbl_mb(tmp1(1)))
      call electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)
      

      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_1get_Tradient',0)

      return
      end



*     ***********************************
*     *					*
*     *		psi_2get_Tgradient	*
*     *					*
*     ***********************************
      subroutine psi_2get_Tgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     *** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      value = MA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_1get_Tradient',0)

      if (option.le.1) then
        call electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                     dbl_mb(rho2_all(1)))
      end if

      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)))

      call electron_gen_hml(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))
      call electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dbl_mb(tmp1(1)),
     >                             THpsi)
      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_1get_Tradient',0)

      return
      end


*     ***********************************
*     *					*
*     *		psi_1get_TMgradient	*
*     *					*
*     ***********************************
      subroutine psi_1get_TMgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8     Eout

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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)))

      Eout =  electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)))

      call electron_get_TMgradient(dcpl_mb(psi1(1)),
     >                            THpsi)

      return
      end



*     ***********************************
*     *					*
*     *		psi_2get_TMgradient	*
*     *					*
*     ***********************************
      subroutine psi_2get_TMgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8     Eout

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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))

      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)))

      call electron_get_TMgradient(dcpl_mb(psi2(1)),
     >                             THpsi)
      
      return
      end


*     ***********************************
*     *					*
*     *	    psi_1geodesic_transport	*
*     *					*
*     ***********************************
      subroutine psi_1geodesic_transport(t,H0)
      implicit none
      real*8 t
      complex*16 H0(*)

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


      call geodesic_transport(t,dcpl_mb(psi1(1)),H0)

      return
      end


*     ***********************************
*     *					*
*     *	    psi_1geodesic_Gtransport	*
*     *					*
*     ***********************************
      subroutine psi_1geodesic_Gtransport(t,G0)
      implicit none
      real*8 t
      complex*16 G0(*)

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

      call geodesic_Gtransport(t,dcpl_mb(psi1(1)),G0)

      return
      end



*     ***********************************
*     *					*
*     *		psi_geodesic_energy 	*
*     *					*
*     ***********************************
      real*8 function psi_geodesic_energy(t)
      implicit none
      real*8 t

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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      real*8 e_new
*     **** external functions ****
      real*8   electron_energy
      external electron_energy

     
      call geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))
      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)))
      e_new =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)))

c      write(*,*) "psi_geodesic_energy:",t,e_new
      psi_geodesic_energy = e_new
      return
      end

*     ***********************************
*     *					*
*     *		psi_geodesic_denergy 	*
*     *					*
*     ***********************************
      real*8 function psi_geodesic_denergy(t)
      implicit none
      real*8 t

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

*     **** external functions ****
      real*8   electron_eorbit
      external electron_eorbit


      call geodesic_transport(t,dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)))
      psi_geodesic_denergy =  2.0d0*electron_eorbit(dcpl_mb(psi2(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_geodesic_final 	*
*     *					*
*     ***********************************
      subroutine psi_geodesic_final(t)
      implicit none
      real*8 t

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

      integer taskid,MASTER
      parameter (MASTER=0)
c     real*8 sum1,sum2
     
      call Parallel_taskid(taskid)

      call geodesic_get(t,dcpl_mb(psi1(1)),
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
*     *			psito2_sd_update	 	*
*     *									*
*     ***********************************
      subroutine psi1to2_sd_update(dte)
      implicit none
      real*8 dte

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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      logical value
      integer nemax,ierr
      integer lmd(2),tmp_L(2)

      call electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                  dbl_mb(rho1_all(1)))

*     **** do a steepest descent step ****
      call electron_sd_update(dcpl_mb(psi1(1)),
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

      call psi_lmbda(ispin,ne,nemax,npack1,
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
*     *									*
*     *		psi_1ke_stress              *
*     *									*
*     ***********************************
      subroutine psi_1ke_stress(stress)
      implicit none
      real*8 stress(3,3)

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

      call ke_euv(ispin,ne,dcpl_mb(psi1(1)),stress)
      return
      end

*     ***********************************
*     *									*
*     *		psi_1coulomb_stress         *
*     *									*
*     ***********************************
      subroutine psi_1coulomb_stress(stress)
      implicit none
      real*8 stress(3,3)

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

      call coulomb_euv(dcpl_mb(dng1(1)),stress)
      return
      end

*     ***********************************
*     *									*
*     *		rho_1exc_stress 	        *
*     *									*
*     ***********************************
      subroutine rho_1exc_stress(stress)
      implicit none
      real*8 stress(3,3)


*     ***** local variables ****
      integer u,v
      real*8 exc,pxc
      real*8 pi,scal,hm(3,3)

*     **** external functions ****
      real*8   rho_1exc,rho_1pxc,lattice_unitg,lattice_omega
      external rho_1exc,rho_1pxc,lattice_unitg,lattice_omega

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     **** LDA part ****
      exc = rho_1exc()
      pxc = rho_1pxc()
      do v=1,3
      do u=1,3
         stress(u,v) = (exc-pxc)*hm(u,v)
      end do
      end do
      write(*,*) "hm(1,1):",hm(1,1),1.0d0/hm(1,1)
      write(*,*) "exc:",exc,pxc
      write(*,*) "D:",stress(1,1)

*     **** GGA part ****
*     **** not finished ****

      return
      end

*     ***********************************
*     *									*
*     *		dng_1vlocal_stress          *
*     *									*
*     ***********************************

      subroutine dng_1vlocal_stress(stress)
      implicit none
      real*8 stress(3,3)


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

      call v_local_euv(dcpl_mb(dng1(1)),stress)

      return
      end

*     ***********************************
*     *									*
*     *		psi_1vnonlocal_stress       *
*     *									*
*     ***********************************
      subroutine psi_1vnonlocal_stress(stress)
      implicit none
      real*8 stress(3,3)

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


*     ***** local variables ****
      integer u,v
      real*8 evnl
      real*8 pi,scal,hm(3,3)

*     **** external functions ****
      real*8   psi_1vnl,lattice_unitg
      external psi_1vnl,lattice_unitg

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)      
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

      call v_nonlocal_euv_2(ispin,ne,dcpl_mb(psi1(1)),stress)
      evnl = psi_1vnl()
      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) - evnl*hm(u,v)
      end do
      end do

      return
      end

*     ***********************************
*     *					*
*     *		psi_1Orb_Analysis       *
*     *					*
*     ***********************************
      subroutine psi_1Orb_Analysis(iunit)
      implicit none
      integer iunit

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

c      call Orb_Analysis(iunit,ispin,ne,dcpl_mb(psi1(1)))
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







*     ***********************************
*     *					*
*     *		psi_1rotate2       	*
*     *					*
*     ***********************************
      subroutine psi_1rotate2()
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

*     ***** local variables *****
      integer ms,index,i,j,shift1,shift2

      call dcopy(2*npack1*(ne(1)+ne(2)),0.0d0,0,dcpl_mb(psi2(1)),1)
      do ms=1,ispin
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)
         do j=1,ne(ms)
           do i=1,ne(ms)
              index = (i-1) + (j-1)*ne(ms) + shift2
             
c              call D3dB_cc_daxpy(1,dbl_mb(hml(1)+index),
c     >                           dcpl_mb(psi1(1)+(i-1+shift1)*nfft3d),
c     >                           dcpl_mb(psi2(1)+(j-1+shift1)*nfft3d)) 
              call Pack_cc_daxpy(1,dbl_mb(hml(1)+index),
     >                           dcpl_mb(psi1(1)+(i-1+shift1)*npack1),
     >                           dcpl_mb(psi2(1)+(j-1+shift1)*npack1)) 

           end do
         end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *		psi_diagonalize_hml	*
*     *					*
*     ***********************************
      subroutine psi_diagonalize_hml()
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


*     ***** local variables ****
      logical value
      integer ms,shift1,shift2,ierr
      integer tmp1(2)


      value = MA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_diagonalize_hml',0)


*     ***** diagonalize the hamiltonian matrix *****
      call dcopy((ne(1)+ne(2)),0.0d0,0,dbl_mb(eig(1)),1)
      do ms=1,ispin
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)
         call eigen(ne(ms),ne(ms),
     >              dbl_mb(hml(1)+shift2),
     >              dbl_mb(eig(1)+shift1),
     >              dbl_mb(tmp1(1)),ierr)
      end do

      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack in psi_diagonalize_hml',0)


      return
      end



*     ***************************
*     *				*
*     *		psi_error	*
*     *				*
*     ***************************
      real*8 function psi_error()
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

*     ***** local variables ****
      logical value
      integer k,n
      real*8  error,sum,size
      integer tmp1(2)

      value = MA_push_get(mt_dcpl,(npack1),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_error',0)


      error = 0.0d0
      size =  dble(ne(1)+ne(2))
      do n=1, (ne(1)+ne(2))
         do k=1,npack1
            dcpl_mb(tmp1(1)+k-1) = dcpl_mb(psi2(1)+k-1+(n-1)*npack1) 
     >                           - dcpl_mb(psi1(1)+k-1+(n-1)*npack1)
         end do
c         call D3dB_cc_dot(1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)
         call Pack_cc_dot(1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)

         error = error + sum
      end do
      error = dsqrt(error)/size
      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_error',0)


      psi_error = error
      return
      end

*     ***************************
*     *				*
*     *		rho_error	*
*     *				*
*     ***************************
      real*8 function rho_error()
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

*     ***** local variables ****
      logical value
      integer k,nx,ny,nz
      real*8  error,scale
      integer tmp1(2)

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega

      value = MA_push_get(mt_dbl,(2*nfft3d),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in rho_error',0)


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scale = lattice_omega()

      scale = (scale)/dble(nx*ny*nz)
*     scale = (scale)/dble(nx*ny*nz)
*     scale = (scale*scale)

      do k=1,(2*nfft3d)
         dbl_mb(tmp1(1)+k-1) = (dbl_mb(rho2(1)+k-1)
     >                         -dbl_mb(rho1(1)+k-1)) 
         dbl_mb(tmp1(1)+k-1) = dbl_mb(tmp1(1)+k-1) 
     >                      + (dbl_mb(rho2(1)+k-1+(ispin-1)*(2*nfft3d))
     >                        -dbl_mb(rho1(1)+k-1+(ispin-1)*(2*nfft3d))) 
      end do
      call D3dB_rr_dot(1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)),error)
      error = error*scale
*     error = dsqrt(error)

      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in rho_error',0)

     

      rho_error = error
      return
      end


*     ***************************
*     *				*
*     *		rho_dipole	*
*     *				*
*     ***************************
      subroutine rho_dipole(dipole)
      implicit none
      real*8 dipole(3)

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

      call Calculate_Dipole(ispin,ne,2*nfft3d,dbl_mb(rho1(1)),dipole)
      return
      end


*     ***************************
*     *				*
*     *		psi_ispin	*
*     *				*
*     ***************************
      integer function psi_ispin()
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

      psi_ispin = ispin
      return
      end


*     ***************************
*     *				*
*     *		psi_ne		*
*     *				*
*     ***************************
      integer function psi_ne(ms)
      implicit none
      integer ms

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

      psi_ne = ne(ms)
      return
      end

*     ***************************
*     *				*
*     *		psi_initialize 	*
*     *				*
*     ***************************

      logical function psi_initialize()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value
      integer nemax
      real*8 sum1,sum2

*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne(ispin,ne)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      nemax = ne(1)+ne(2)

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
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
     >        MA_alloc_get(mt_dbl,(2*nemax*nemax),'hml',hml(2),hml(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) call errquit('out of heap memory',0)

*     *****  read initial wavefunctions into psi1  ****
      call psi_read(ispin,ne,dcpl_mb(psi1(1)))

*     **** Ortho Check ****
      call Grsm_gg_trace(npack1,(ne(1)+ne(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi1(1)),
     >                        sum2)

      sum1 = dble(ne(1) + ne(2))
      if (dabs(sum2-sum1).gt.1.0d-10) then
         call Parallel_taskid(taskid)
         call Grsm_g_MakeOrtho(npack1,ne(1),dcpl_mb(psi1(1)))
         if (ispin.gt.1) then
           call Grsm_g_MakeOrtho(npack1,ne(2),
     >                           dcpl_mb(psi1(1) + ne(1)*npack1))
         end if
         call Grsm_gg_trace(npack1,(ne(1)+ne(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi1(1)),
     >                        sum2)
         if (taskid.eq.MASTER)
     >    write(*,*) "Waring - Gram-Schmidt being performed on psi:",
     >               dabs(sum2-sum1)


      end if



      psi_initialize = value
      return
      end



*     ***************************
*     *				*
*     *		psi_finalize	*
*     *				*
*     ***************************

      logical function psi_finalize()
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

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all


*     **** local variables ****
      logical value

*     ***** write psi1 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi1(1)))
    
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

      psi_finalize = value
      return
      end



*     ***********************
*     *						*
*     *		psi_Mulliken	*
*     *						*
*     ***********************

      subroutine psi_Mulliken(rtdb)
      implicit none 
      integer rtdb

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

*     **** Lubin Water Analysis ***
      call pspw_Lubin_water_analysis(rtdb,ispin,ne,2*nfft3d,
     >                                 dbl_mb(rho1(1)))

      call pspw_analysis(rtdb,ispin,ne,dcpl_mb(psi2(1)),
     >                                   dbl_mb(eig(1)))

      call pspw_gen_APC(ispin,ne,dcpl_mb(dng1(1)))
      call pspw_print_APC(6)

      return
      end

