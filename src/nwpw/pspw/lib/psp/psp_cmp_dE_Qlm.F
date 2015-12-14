c     *****************************************
c     *                                       *
c     *         psp_cmp_dE_Qlm_pw             *
c     *                                       *
c     *****************************************
c    0.5*ncmp*ncmp + ncmp*n
      subroutine psp_cmp_dE_Qlm_pw(ispin,move,fion)
      implicit none
      integer ispin
      logical move
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     **** local variables ****
      logical periodic,value
      integer npack0,nfft3d,nx,ny,nz
      integer dng_cmp(2),vcmp(2)
      real*8  dv,scal1,omega

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      integer  control_version
      external control_version

      periodic = (control_version().eq.3)

      !****************
      !*** periodic ***
      !****************
      if (periodic) then
         omega = lattice_omega()
         call Pack_npack(0,npack0)
         value = BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                    dng_cmp(2),dng_cmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'vcmp',vcmp(2),vcmp(1))
         if (.not.value)
     >     call errquit('psp_cmp_dE_Qlm_pw:out stack',0,MA_ERR)

         !*** generate vcmp ****
         call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
         call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp(1)))
         call Pack_c_Copy(0,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))

         !*** generate vcmp = vcmp+vc ****
         call Pack_cc_Sum2(0,
     >                     dcpl_mb(vc_tmp(1)),
     >                     dcpl_mb(vcmp(1)))

         !*** generate dEcmp/dQlm ****
         call Pack_c_SMul1(0,omega,dcpl_mb(vcmp(1)))
         call nwpw_compcharge_gen_dE_Qlm(1,ispin,
     >                                   dcpl_mb(vcmp(1)),
     >                                   dcpl_mb(vcmp(1)),
     >                                   move,fion)

         value =           BA_pop_stack(vcmp(2))
         value = value.and.BA_pop_stack(dng_cmp(2))
         if (.not.value)
     >      call errquit('psp_cmp_dE_Qlm_pw:popping stack',1,MA_ERR)


      !*****************
      !*** aperiodic ***
      !*****************
      else
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         dv = lattice_omega()*scal1

         call D3dB_nfft3d(1,nfft3d)
         value = BA_push_get(mt_dcpl,nfft3d,'dng_cmp',
     >                    dng_cmp(2),dng_cmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'vcmp',vcmp(2),vcmp(1))
         if (.not.value)
     >     call errquit('psp_cmp_dE_Qlm_pw:out stack',2,MA_ERR)
         call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
         call Pack_c_unpack(0,dcpl_mb(dng_cmp(1)))
         call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp(1)))

         !*** generate vcmp ****
         call coulomb2_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp(1)))
         call D3dB_r_Copy(0,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))

         !*** generate vcmp = vcmp+vc ****
         call D3dB_rr_Sum2(0,
     >                     dcpl_mb(vc_tmp(1)),
     >                     dcpl_mb(vcmp(1)))
         call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp(1)))
         call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp(1)))
         call Pack_c_pack(0,dcpl_mb(vcmp(1)))

         !*** generate dEcmp/dQlm ***
         call nwpw_compcharge_gen_dE_Qlm(1,ispin,
     >                                   dcpl_mb(vcmp(1)),
     >                                   dcpl_mb(vcmp(1)),
     >                                   move,fion)

         value =           BA_pop_stack(vcmp(2))
         value = value.and.BA_pop_stack(dng_cmp(2))
         if (.not.value)
     >      call errquit('psp_cmp_dE_Qlm_pw:popping stack',3,MA_ERR)
      end if

      return
      end 

c     *****************************************
c     *                                       *
c     *         psp_cmp_dE_Qlm                *
c     *                                       *
c     *****************************************

      subroutine psp_cmp_dE_Qlm(ispin,move,fion)
      implicit none
      integer ispin
      logical move
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     **** local variables ****
      logical periodic,value
      integer npack0,nfft3d,nx,ny,nz
      integer dng_cmp(2),dng_cmp_smooth(2),vcmp(2),vcmp_smooth(2)
      real*8  dv,scal1,omega

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      integer  control_version  
      external control_version  

      periodic = (control_version().eq.3)

      !****************
      !*** periodic ***
      !****************
      if (periodic) then

         omega = lattice_omega()
         call Pack_npack(0,npack0)
         value = BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                       dng_cmp(2),dng_cmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                       dng_cmp_smooth(2),dng_cmp_smooth(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'vcmp',
     >                       vcmp(2),vcmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'vcmp_smooth',
     >                       vcmp_smooth(2),vcmp_smooth(1))
         if (.not.value) 
     >      call errquit('psp_cmp_dE_Qlm:out stack',0,MA_ERR)

         call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                    dcpl_mb(dng_cmp(1)),
     >                                    dcpl_mb(dng_cmp_smooth(1)))

         !*** compute hartree potential of ntilde + ncmp_tilde ***
         !*** compute hartree potential ncmp   - ncmp_tilde ***
         call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp(1)))
         call coulomb_v(dcpl_mb(dng_cmp_smooth(1)),
     >                  dcpl_mb(vcmp_smooth(1)))
         call Pack_c_Copy(0,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))

         call Pack_cc_Sub2(0,
     >                     dcpl_mb(vcmp_smooth(1)),
     >                     dcpl_mb(vcmp(1)))

         call Pack_cc_Sum2(0,dcpl_mb(vc_tmp(1)),
     >                       dcpl_mb(vcmp_smooth(1)))

         call Pack_c_SMul1(0,omega,dcpl_mb(vcmp(1)))
         call Pack_c_SMul1(0,omega,dcpl_mb(vcmp_smooth(1)))

         call nwpw_compcharge_gen_dE_Qlm(3,ispin,
     >                                   dcpl_mb(vcmp_smooth(1)),
     >                                   dcpl_mb(vcmp(1)),
     >                                   move,fion)

         value =           BA_pop_stack(vcmp_smooth(2))
         value = value.and.BA_pop_stack(vcmp(2))
         value = value.and.BA_pop_stack(dng_cmp_smooth(2))
         value = value.and.BA_pop_stack(dng_cmp(2))
         if (.not.value) 
     >     call errquit('psp_cmp_dE_Qlm:pop stack',1,MA_ERR)


      !******************
      !*** aperiodic ****
      !******************
      else
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         dv = lattice_omega()*scal1

         call D3dB_nfft3d(1,nfft3d)
         value = BA_push_get(mt_dcpl,nfft3d,'dng_cmp',
     >                       dng_cmp(2),dng_cmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'dng_cmp_smooth',
     >                       dng_cmp_smooth(2),dng_cmp_smooth(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'vcmp',
     >                       vcmp(2),vcmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'vcmp_smooth',
     >                       vcmp_smooth(2),vcmp_smooth(1))
         if (.not.value)
     >      call errquit('psp_cmp_dE_Qlm2:out stack',2,MA_ERR)

         call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                    dcpl_mb(dng_cmp(1)),
     >                                    dcpl_mb(dng_cmp_smooth(1)))

         !*** dng_cmp(G),dng_cmp_smooth(G) --> dng_cmp(r),dng_cmp_smooth(r) ***
         call Pack_c_unpack(0,dcpl_mb(dng_cmp(1)))
         call Pack_c_unpack(0,dcpl_mb(dng_cmp_smooth(1)))
         call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp(1)))
         call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp_smooth(1)))

         !*** generate vcmp and vcmp_smooth ***
         call coulomb2_v(dcpl_mb(dng_cmp(1)),
     >                   dcpl_mb(vcmp(1)))
         call coulomb2_v(dcpl_mb(dng_cmp_smooth(1)),
     >                   dcpl_mb(vcmp_smooth(1)))
         call D3dB_r_Copy(1,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))


         !*** vcmp        = vcmp-vcmp_smooth ***
         !*** vcmp_smooth = vcmp_smooth + vc ***
         call D3dB_rr_Sub2(1,
     >                     dcpl_mb(vcmp_smooth(1)),
     >                     dcpl_mb(vcmp(1)))
         call D3dB_rr_Sum2(1,
     >                     dcpl_mb(vc_tmp(1)),
     >                     dcpl_mb(vcmp_smooth(1)))

         call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp(1)))
         call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp_smooth(1)))
         call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp(1)))
         call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp_smooth(1)))
         call Pack_c_pack(0,dcpl_mb(vcmp(1)))
         call Pack_c_pack(0,dcpl_mb(vcmp_smooth(1)))

         call nwpw_compcharge_gen_dE_Qlm(3,ispin,
     >                                   dcpl_mb(vcmp_smooth(1)),
     >                                   dcpl_mb(vcmp(1)),move,fion)

         value =           BA_pop_stack(vcmp_smooth(2))
         value = value.and.BA_pop_stack(vcmp(2))
         value = value.and.BA_pop_stack(dng_cmp_smooth(2))
         value = value.and.BA_pop_stack(dng_cmp(2))
         if (.not.value)
     >      call errquit('psp_cmp_dE_Qlm:pop stack',3,MA_ERR)

      end if
         
      call nwpw_compcharge_gen_dEmult_Qlm(ispin)
      call nwpw_compcharge_add_dEmult_Qlm(ispin)

      return
      end


