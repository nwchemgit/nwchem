*
* $Id$
*

*     ***********************************
*     *					*
*     *	 	  psp_init  		*
*     *					*
*     ***********************************

      subroutine psp_init()      
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal

      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer npack1,npack0,nion,i,j,ii,n2ft3d,l,m,lm
      integer tmp_rgrid(2)
      real*8  r2

*     **** external functions *****
      logical  control_pspspin,control_use_grid_cmp,FMM_fmm
      external control_pspspin,control_use_grid_cmp,FMM_fmm
      integer  ion_nkatm,ion_nion,control_pspnuterms,FMM_lmax
      external ion_nkatm,ion_nion,control_pspnuterms,FMM_lmax
      integer  control_version
      external control_version
 
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      npsp = ion_nkatm()

c      value = BA_alloc_get(mt_dcpl,(nmax_max*lmmax_max*npack1),
c     >                    'prjtmp',prjtmp(2),prjtmp(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,(nmax_max*lmmax_max*npsp*npack1),
c     >                    'vnl',vnl(2),vnl(1))
c     value = value.and.
c    >         BA_alloc_get(mt_dbl,(lmmax_max*npsp),
c    >                    'vnlnrm',vnlnrm(2),vnlnrm(1))

      value = BA_alloc_get(mt_dbl,(npsp*npack0),'vl',vl(2),vl(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp*npack0),
     >                     'vlpaw',vlpaw(2),vlpaw(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'vnl',vnl(2),vnl(1))
c      value = value.and.
c     >      BA_alloc_get(mt_dbl,(nmax_max*nmax_max*(lmax_max+1)*npsp),
c     >                    'Gijl',Gijl(2),Gijl(1))
      value = value.and.BA_alloc_get(mt_int,npsp,'Gijl',Gijl(2),Gijl(1))

      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'nprj',nprj(2),nprj(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'n_projector',n_projector(2),n_projector(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'l_projector',l_projector(2),l_projector(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'm_projector',m_projector(2),m_projector(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'b_projector',b_projector(2),b_projector(1))

      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp),'zv',zv(2),zv(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp),'amass',amass(2),amass(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp*(lmax_max+1)),'rc',rc(2),rc(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'lmmax',lmmax(2),lmmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'lmax',lmax(2),lmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'locp',locp(2),locp(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'nmax',nmax(2),nmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),
     >                     'psp_type',psp_type(2),psp_type(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp),'rlocal',rlocal(2),rlocal(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'r3_matrix',
     >                     r3_matrix(2),r3_matrix(1))

*     ***** PAW variables *****
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'hartree_matrix',
     >                     hartree_matrix(2),hartree_matrix(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'comp_charge_matrix',
     >                     comp_charge_matrix(2),comp_charge_matrix(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'comp_pot_matrix',
     >                     comp_pot_matrix(2),comp_pot_matrix(1))
      value = value.and.
     > BA_alloc_get(mt_dbl,(npsp),'log_amesh',log_amesh(2),log_amesh(1))
      value = value.and.
     >  BA_alloc_get(mt_dbl,(npsp),'r1',r1(2),r1(1))
      value = value.and.
     >  BA_alloc_get(mt_dbl,(npsp),'rmax',rmax(2),rmax(1))
      value = value.and.
     >  BA_alloc_get(mt_dbl,(npsp),'sigma',sigma(2),sigma(1))
      value = value.and.
     >  BA_alloc_get(mt_dbl,(npsp),'zion',zion(2),zion(1))
      value = value.and.
     >  BA_alloc_get(mt_dbl,(npsp),'core_kin',core_kin(2),core_kin(1))
      value = value.and.
     >  BA_alloc_get(mt_dbl,(npsp),'core_ion',core_ion(2),core_ion(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'n1dgrid',n1dgrid(2),n1dgrid(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,
     >                     'n1dbasis',n1dbasis(2),n1dbasis(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'nae',nae(2),nae(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'nps',nps(2),nps(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'lps',lps(2),lps(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'icut',icut(2),icut(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'eig',eig(2),eig(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'phi_ae',phi_ae(2),phi_ae(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'dphi_ae',dphi_ae(2),dphi_ae(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'phi_ps',phi_ps(2),phi_ps(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'dphi_ps',dphi_ps(2),dphi_ps(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'core_ae',core_ae(2),core_ae(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'core_ps',core_ps(2),core_ps(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'core_ae_prime',
     >                    core_ae_prime(2),core_ae_prime(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'core_ps_prime',
     >                     core_ps_prime(2),core_ps_prime(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'rgrid',rgrid(2),rgrid(1))

*     **** setup pspspin structure - used for generating antiferromagnetic structures ****
      pspspin = control_pspspin()
      if (pspspin) then
         nion  = ion_nion()
         value = value.and.
     >   BA_alloc_get(mt_log,nion,'pspspin_upions',
     >                pspspin_upions(2),pspspin_upions(1))
         value = value.and.
     >   BA_alloc_get(mt_log,nion,'pspspin_downions',
     >                pspspin_downions(2),pspspin_downions(1))
         value = value.and.
     >   BA_alloc_get(mt_int,nion,'pspspin_upl',
     >                pspspin_upl(2),pspspin_upl(1))
         value = value.and.
     >   BA_alloc_get(mt_int,nion,'pspspin_downl',
     >                pspspin_downl(2),pspspin_downl(1))
         value = value.and.
     >   BA_alloc_get(mt_int,nion,'pspspin_upm',
     >                pspspin_upm(2),pspspin_upm(1))
         value = value.and.
     >   BA_alloc_get(mt_int,nion,'pspspin_downm',
     >                pspspin_downm(2),pspspin_downm(1))
         value = value.and.
     >   BA_alloc_get(mt_dbl,nion,'pspspin_upscale',
     >                pspspin_upscale(2),pspspin_upscale(1))
         value = value.and.
     >   BA_alloc_get(mt_dbl,nion,'pspspin_downscale',
     >                pspspin_downscale(2),pspspin_downscale(1))
         if (.not. value) 
     >   call errquit('psp_init:out of heap memory',0, MA_ERR)

         do ii=1,nion
            log_mb(pspspin_upions(1)+ii-1)    = .false.
            log_mb(pspspin_downions(1)+ii-1)  = .false.
            dbl_mb(pspspin_upscale(1)+ii-1)   = 1.0d0
            dbl_mb(pspspin_downscale(1)+ii-1) = 1.0d0
            int_mb(pspspin_upl(1)+ii-1)       = 0
            int_mb(pspspin_downl(1)+ii-1)     = 0
            int_mb(pspspin_upm(1)+ii-1)       = 0
            int_mb(pspspin_downm(1)+ii-1)     = 0
         end do
         call control_set_pspspin(nion,
     >                            dbl_mb(pspspin_upscale(1)),
     >                            dbl_mb(pspspin_downscale(1)),
     >                            int_mb(pspspin_upl(1)),
     >                            int_mb(pspspin_downl(1)),
     >                            int_mb(pspspin_upm(1)),
     >                            int_mb(pspspin_downm(1)),
     >                            log_mb(pspspin_upions(1)),
     >                            log_mb(pspspin_downions(1)))
      end if

      if (.not. value) 
     >   call errquit('psp_init:out of heap memory',0, MA_ERR)

      call dcopy(npsp*npack0,          0.0d0,0,dbl_mb(vl(1)), 1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(zv(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(amass(1)),1)
      call dcopy(npsp*(lmax_max+1),    0.0d0,0,dbl_mb(rc(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(rlocal(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(core_kin(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(core_ion(1)),1)

*     **** allocate semicore data ****
      call semicore_init()

*     *** set pawexist and use_grid_cmp ****
      pawexist = .false.
      use_grid_cmp = control_use_grid_cmp()


*     **** set FMM variables ****
      psp_fmm = ((control_version().eq.4).and.FMM_fmm())
      if (psp_fmm) then
         psp_fmm_lmax = FMM_lmax()
         call D3dB_n2ft3d(1,n2ft3d)

         value=BA_alloc_get(mt_dbl,n2ft3d*(psp_fmm_lmax+1)**2,
     >          'psp_fmm_rTlm',psp_fmm_rTlm(2),psp_fmm_rTlm(1))
         value=value.and.
     >         BA_alloc_get(mt_dbl,(psp_fmm_lmax+1)**2,
     >          'psp_fmm_Llm',psp_fmm_Llm(2),psp_fmm_Llm(1))
         value=value.and.
     >         BA_alloc_get(mt_dbl,(psp_fmm_lmax+1)**2,
     >          'psp_fmm_Mlm',psp_fmm_Mlm(2),psp_fmm_Mlm(1))
         if (.not.value) 
     >      call errquit('psp_init:out of heap memory',0, MA_ERR)

         value=BA_push_get(mt_dbl,3*n2ft3d,
     >          'tmp_rgrid',tmp_rgrid(2),tmp_rgrid(1))
         if (.not.value) 
     >      call errquit('psp_init:out of stack memory',0,MA_ERR)
         call lattice_r_grid(dbl_mb(tmp_rgrid(1)))
         psp_fmm_rmax2 = 0.0d0
         do i=1,n2ft3d
            r2 = dbl_mb(tmp_rgrid(1)+3*(i-1)  )**2
     >         + dbl_mb(tmp_rgrid(1)+3*(i-1+1))**2
     >         + dbl_mb(tmp_rgrid(1)+3*(i-1)+2)**2
            if (r2.gt.psp_fmm_rmax2) psp_fmm_rmax2 = r2
         end do
         call Parallel_MaxAll(psp_fmm_rmax2)

         lm = 0
         do l=0,psp_fmm_lmax
         do m=-l,l
            call Tesseral3_rgrid_lm_rl(l,m,
     >                              n2ft3d,dbl_mb(rgrid(1)),
     >                              dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d))
            lm = lm + 1
         end do
         end do

         if (.not.BA_pop_stack(tmp_rgrid(2)))
     >      call errquit('psp_init:popping stack',0,MA_ERR)
      end if

      return
      end

*     ****************************************************
*     *                                                  *
*     *                psp_proj_init                     *
*     *                                                  *
*     ****************************************************
      subroutine psp_proj_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     ***** local variables ****
      logical value,periodic
      integer npack1,npack0,n2ft3d

*     ***** external functions *****
      integer  psp_nprj_max,control_nprj_mult,control_version
      external psp_nprj_max,control_nprj_mult,control_version

      periodic = (control_version().eq.3)

      if (periodic) then
         call Pack_npack(0,npack0)
      else
         call D3dB_n2ft3d(1,n2ft3d)
         npack0 = n2ft3d/2
      end if

      call Pack_npack(1,npack1)
      nprj_mult= control_nprj_mult()
      nprj_max = psp_nprj_max()
      value = BA_alloc_get(mt_dcpl,npack1*nprj_max*nprj_mult,
     >                     'prjtmp',prjtmp(2),prjtmp(1))
      if (pawexist) then
         value = value.and.
     >            BA_alloc_get(mt_dbl,2*nprj_max*nprj_max,
     >                        'wtmp',wtmp(2),wtmp(1))
         value = value.and.
     >            BA_alloc_get(mt_dcpl,npack0,'vc_tmp',
     >                         vc_tmp(2),vc_tmp(1))
         value = value.and.
     >            BA_alloc_get(mt_dcpl,npack0,'vcmptmp',
     >                         vcmp_tmp(2),vcmp_tmp(1))
      end if
      if (.not. value)
     >   call errquit('psp_proj_init:out of heap memory',0,MA_ERR)
      return
      end


      subroutine psp_set_vc(vcin)
      implicit none
      complex*16 vcin(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"

      logical  periodic
      integer  control_version
      external control_version

      periodic = (control_version().eq.3)

      if (periodic) then
         call Pack_c_Copy(0,vcin,dcpl_mb(vc_tmp(1))) 
      else
         call D3dB_r_Copy(1,vcin,dcpl_mb(vc_tmp(1))) 
      end if

      return
      end 


      subroutine psp_get_vcmp(vcmpout)
      implicit none
      complex*16 vcmpout(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"

      logical  periodic
      integer  control_version
      external control_version

      periodic = (control_version().eq.3)

      if (periodic) then
         call Pack_c_Copy(0,dcpl_mb(vcmp_tmp(1)),vcmpout) 
      else
         call D3dB_r_Copy(1,dcpl_mb(vcmp_tmp(1)),vcmpout) 
      end if

      return
      end 


      subroutine psp_add_vcmp(vall)
      implicit none
      complex*16 vall(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"

      logical  periodic
      integer  control_version
      external control_version

      periodic = (control_version().eq.3)

      if (periodic) then
         call Pack_cc_Sum2(0,dcpl_mb(vcmp_tmp(1)),vall)
      else
         call D3dB_rr_Sum2(1,dcpl_mb(vcmp_tmp(1)),vall)
      end if
      return
      end




*     ***********************************
*     *					*
*     *	 	  psp_end  		*
*     *					*
*     ***********************************

      subroutine psp_end()      
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal


      logical value
      integer ia

*     **** external functions ****

*     **** deallocate FMM data ****
      if (psp_fmm) then
         value =           BA_free_heap(psp_fmm_rTlm(2))
         value = value.and.BA_free_heap(psp_fmm_Llm(2))
         value = value.and.BA_free_heap(psp_fmm_Mlm(2))
         if (.not.value)
     >     call errquit('psp_end:error freeing fmm heap',0,MA_ERR)
      end if

*     **** deallocate semicore data ****
      call semicore_end()

*     **** deallocate nonlocal pspw data ****
      do ia=1,npsp
        if (int_mb(nprj(1)+ia-1).gt.0) then
           call psi_data_dealloc(int_mb(vnl(1)+ia-1))
           call psi_data_dealloc(int_mb(Gijl(1)+ia-1))
        end if
        if (int_mb(psp_type(1)+ia-1).eq.4) then
           call psi_data_dealloc(int_mb(hartree_matrix(1)+ia-1))
           call psi_data_dealloc(int_mb(comp_charge_matrix(1)+ia-1))
           call psi_data_dealloc(int_mb(comp_pot_matrix(1)+ia-1))
           call psi_data_dealloc(int_mb(eig(1)+ia-1))
           call psi_data_dealloc(int_mb(phi_ae(1)+ia-1))
           call psi_data_dealloc(int_mb(dphi_ae(1)+ia-1))
           call psi_data_dealloc(int_mb(phi_ps(1)+ia-1))
           call psi_data_dealloc(int_mb(dphi_ps(1)+ia-1))
           call psi_data_dealloc(int_mb(core_ae(1)+ia-1))
           call psi_data_dealloc(int_mb(core_ps(1)+ia-1))
           call psi_data_dealloc(int_mb(core_ae_prime(1)+ia-1))
           call psi_data_dealloc(int_mb(core_ps_prime(1)+ia-1))
           call psi_data_dealloc(int_mb(rgrid(1)+ia-1))
        end if
        if (int_mb(psp_type(1)+ia-1).eq.9) then
           call psi_data_dealloc(int_mb(r3_matrix(1)+ia-1))
        end if
      end do

      value = BA_free_heap(prjtmp(2))
      if (pawexist) then
         call psp_paw_end()
         value = value.and.BA_free_heap(wtmp(2))
         value = value.and.BA_free_heap(vc_tmp(2))
         value = value.and.BA_free_heap(vcmp_tmp(2))
      end if
      value = value.and.BA_free_heap(vl(2))
      value = value.and.BA_free_heap(vnl(2))
c     value = value.and.BA_free_heap(vnlnrm(2))
      value = value.and.BA_free_heap(Gijl(2))
      value = value.and.BA_free_heap(nprj(2))
      value = value.and.BA_free_heap(n_projector(2))
      value = value.and.BA_free_heap(l_projector(2))
      value = value.and.BA_free_heap(m_projector(2))
      value = value.and.BA_free_heap(b_projector(2))
      value = value.and.BA_free_heap(zv(2))
      value = value.and.BA_free_heap(amass(2))
      value = value.and.BA_free_heap(rc(2))
      value = value.and.BA_free_heap(lmmax(2))
      value = value.and.BA_free_heap(lmax(2))
      value = value.and.BA_free_heap(locp(2))
      value = value.and.BA_free_heap(nmax(2))
      value = value.and.BA_free_heap(psp_type(2))
      value = value.and.BA_free_heap(rlocal(2))

*     **** PAW variables ****
      if (pawexist) then
         call nwpw_compcharge_end()
         call nwpw_xc_end()
      end if
      value = value.and.BA_free_heap(vlpaw(2))
      value = value.and.BA_free_heap(r3_matrix(2))
      value = value.and.BA_free_heap(hartree_matrix(2))
      value = value.and.BA_free_heap(comp_charge_matrix(2))
      value = value.and.BA_free_heap(comp_pot_matrix(2))
      value = value.and.BA_free_heap(log_amesh(2))
      value = value.and.BA_free_heap(r1(2))
      value = value.and.BA_free_heap(rmax(2))
      value = value.and.BA_free_heap(sigma(2))
      value = value.and.BA_free_heap(zion(2))
      value = value.and.BA_free_heap(core_kin(2))
      value = value.and.BA_free_heap(core_ion(2))
      value = value.and.BA_free_heap(n1dgrid(2))
      value = value.and.BA_free_heap(n1dbasis(2))
      value = value.and.BA_free_heap(nae(2))
      value = value.and.BA_free_heap(nps(2))
      value = value.and.BA_free_heap(lps(2))
      value = value.and.BA_free_heap(icut(2))
      value = value.and.BA_free_heap(eig(2))
      value = value.and.BA_free_heap(phi_ae(2))
      value = value.and.BA_free_heap(dphi_ae(2))
      value = value.and.BA_free_heap(phi_ps(2))
      value = value.and.BA_free_heap(dphi_ps(2))
      value = value.and.BA_free_heap(core_ae(2))
      value = value.and.BA_free_heap(core_ps(2))
      value = value.and.BA_free_heap(core_ae_prime(2))
      value = value.and.BA_free_heap(core_ps_prime(2))
      value = value.and.BA_free_heap(rgrid(2))

      if (pspspin) then
         value = value.and.BA_free_heap(pspspin_upions(2))
         value = value.and.BA_free_heap(pspspin_downions(2))
         value = value.and.BA_free_heap(pspspin_upm(2))
         value = value.and.BA_free_heap(pspspin_downm(2))
         value = value.and.BA_free_heap(pspspin_upl(2))
         value = value.and.BA_free_heap(pspspin_downl(2))
         value = value.and.BA_free_heap(pspspin_upscale(2))
         value = value.and.BA_free_heap(pspspin_downscale(2))
      end if


      if (.not. value) 
     >  call errquit('psp_end:error freeing heap memory',0,MA_ERR)

*     **** deallocate prj_indx ****
      call psp_prj_indx_end()

      return
      end


*     ****************************************************
*     *                                                  *
*     *             psp_prj_indx_init                    *
*     *                                                  *
*     ****************************************************
*
*     This routine sets up the prj_indx indexes:
*           shift_prj_indx,ii_prj_indx,ia_prj_indx,sd_function_prj_indx
*
      subroutine psp_prj_indx_init()
      implicit none

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,sd_function
      integer k,l,ii,ia,nproj,l_prj,m_prj,count1,count2,shift,nion

*     **** external functions ****
      integer  ion_nion,ion_katm,psi_data_get_ptr
      external ion_nion,ion_katm,psi_data_get_ptr

      nion = ion_nion()
      n_prj_indx    = 0
      nion_prj_indx = 0
      do ii=1,nion
        ia=ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then
           nion_prj_indx = nion_prj_indx + 1

        do l=1,nproj
           n_prj_indx = n_prj_indx + 1
        end do
        end if !** nproj>0 **
      end do !** ii **

      value = BA_alloc_get(mt_int,nion_prj_indx,'ii_prj_indx',
     >                     ii_prj_indx(2),ii_prj_indx(1))
      value = value.and.
     >        BA_alloc_get(mt_int,nion_prj_indx,'ia_prj_indx',
     >                     ia_prj_indx(2),ia_prj_indx(1))
      value = value.and.
     >        BA_alloc_get(mt_int,nion_prj_indx,'nproj_prj_indx',
     >                     nproj_prj_indx(2),nproj_prj_indx(1))
      value = value.and.
     >        BA_alloc_get(mt_int,nion_prj_indx,'swstart_prj_indx',
     >                     swstart_prj_indx(2),swstart_prj_indx(1))

      value = BA_alloc_get(mt_int,n_prj_indx,'shift_prj_indx',
     >                     shift_prj_indx(2),shift_prj_indx(1))
      value = BA_alloc_get(mt_int,n_prj_indx,'l_prj_prj_indx',
     >                     l_prj_prj_indx(2),l_prj_prj_indx(1))
      value = BA_alloc_get(mt_int,n_prj_indx,'m_prj_prj_indx',
     >                     m_prj_prj_indx(2),m_prj_prj_indx(1))
      value = value.and.
     >        BA_alloc_get(mt_log,n_prj_indx,'sd_function_prj_indx',
     >              sd_function_prj_indx(2),sd_function_prj_indx(1))
      if (.not.value)
     >  call errquit('psp_prj_indx_init: out of memory',0, MA_ERR)
      swaset = .false.
      swann  = 0

      count1 = 0
      count2 = 0
      do ii=1,nion
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)
        if (nproj.gt.0) then
           int_mb(ii_prj_indx(1)+count1)      = ii
           int_mb(ia_prj_indx(1)+count1)      = ia
           int_mb(nproj_prj_indx(1)+count1)   = nproj
           int_mb(swstart_prj_indx(1)+count1) = count2
           count1 = count1 + 1

           do l=1,nproj
              shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                  + (ia-1)*(nmax_max*lmmax_max))
              m_prj = int_mb(m_projector(1)+(l-1)
     >                                  + (ia-1)*(nmax_max*lmmax_max))
#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

              int_mb(shift_prj_indx(1)+count2)       = shift
              int_mb(l_prj_prj_indx(1)+count2)       = l_prj
              int_mb(m_prj_prj_indx(1)+count2)       = m_prj
              log_mb(sd_function_prj_indx(1)+count2) = sd_function

              count2 = count2 + 1
           end do
        end if !** nproj>0 **
      end do !** ii **

      return 
      end


*     ****************************************************
*     *                                                  *
*     *         psp_prj_indx_alloc_sw1a_sw2a             *
*     *                                                  *
*     ****************************************************
      subroutine psp_prj_indx_alloc_sw1a_sw2a(nn)
      implicit none
      integer nn

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     **** local variables ****
      logical ok
      integer nthr,nthr0
      integer  Parallel_maxthreads,Parallel_nthreads
      external Parallel_maxthreads,Parallel_nthreads

      nthr  = Parallel_maxthreads()
      nthr0 = Parallel_nthreads()
      if (nthr0.gt.nthr) nthr = nthr0

      if (swaset.and.(nn.gt.swann)) then
         ok =        BA_free_heap(sw1a(2))
         ok = ok.and.BA_free_heap(sw2a(2))
         if (.not.ok)
     >     call errquit('psp_prj_indx_alloc_sw1a_sw2a:freeing heap',
     >                  0,MA_ERR)
!$OMP MASTER
         swaset = .false.
!$OMP END MASTER
!$OMP BARRIER
      end if

      if (.not.swaset) then
         ok = BA_alloc_get(mt_dbl,nn*n_prj_indx,'sw1a',sw1a(2),sw1a(1))
         ok = ok.and.
     >        BA_alloc_get(mt_dbl,nn*n_prj_indx,'sw2a',sw2a(2),sw2a(1))
         ok = ok.and.
     >        BA_alloc_get(mt_dbl,nthr*nn*n_prj_indx,
     >                     'sw1t',sw1t(2),sw1t(1))
         if (.not.ok)
     >     call errquit('psp_prj_indx_alloc_sw1a_sw2a: out of memory',
     >                  0,MA_ERR)
!$OMP MASTER
         swaset = .true.
         swann  = nn
!$OMP END MASTER
!$OMP BARRIER
      end if

      return
      end



*     ****************************************************
*     *                                                  *
*     *             psp_prj_indx_end                     *
*     *                                                  *
*     ****************************************************
      subroutine psp_prj_indx_end()
      implicit none

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

      value =           BA_free_heap(ii_prj_indx(2))
      value = value.and.BA_free_heap(ia_prj_indx(2))
      value = value.and.BA_free_heap(nproj_prj_indx(2))
      value = value.and.BA_free_heap(swstart_prj_indx(2))
      value = value.and.BA_free_heap(shift_prj_indx(2))
      value = value.and.BA_free_heap(l_prj_prj_indx(2))
      value = value.and.BA_free_heap(m_prj_prj_indx(2))
      value = value.and.BA_free_heap(sd_function_prj_indx(2))
      if (swaset) then
         value = value.and.BA_free_heap(sw1a(2))
         value = value.and.BA_free_heap(sw2a(2))
         value = value.and.BA_free_heap(sw1t(2))
      end if
      if (.not.value)
     >  call errquit('psp_prj_indx_end:error freeing heap',0,MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *	 	   psp_rlocal		*
*     *					*
*     ***********************************
      real*8 function psp_rlocal(ia)
      implicit none
      integer ia

#include "bafdecls.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal

      psp_rlocal = dbl_mb(rlocal(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_r1               *
*     *					*
*     ***********************************
      real*8 function psp_r1(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_r1 = dbl_mb(r1(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_rmax             *
*     *					*
*     ***********************************
      real*8 function psp_rmax(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_rmax = dbl_mb(rmax(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_sigma            *
*     *					*
*     ***********************************
      real*8 function psp_sigma(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_sigma = dbl_mb(sigma(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	       psp_sphere_radius        *
*     *					*
*     ***********************************
      real*8 function psp_sphere_radius(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      integer  indx
      integer  psi_data_get_chnk
      external psi_data_get_chnk

      indx = psi_data_get_chnk(int_mb(rgrid(1)+ia-1))
      psp_sphere_radius = dbl_mb(indx + int_mb(icut(1)+ia-1)-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_zv		*
*     *					*
*     ***********************************
      real*8 function psp_zv(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_zv = dbl_mb(zv(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   psp_zion		*
*     *					*
*     ***********************************
      real*8 function psp_zion(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_zion = dbl_mb(zion(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *            psp_zv_ptr           *
*     *                                 *
*     ***********************************
      integer function psp_zv_ptr()
      implicit none

#include "psp.fh"

      psp_zv_ptr = zv(1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_amass		*
*     *					*
*     ***********************************
      real*8 function psp_amass(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_amass = dbl_mb(amass(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_rc		*
*     *					*
*     ***********************************
      real*8 function psp_rc(i,ia)
      implicit none
      integer i,ia

#include "bafdecls.fh"
#include "psp.fh"

c     psp_rc = rc(i,ia)
      psp_rc = dbl_mb(rc(1) + i + (lmax_max+1)*(ia-1))
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_atom		*
*     *					*
*     ***********************************
      character*2 function psp_atom(ia)
      implicit none
      integer  ia

#include "psp.fh"

      psp_atom = atom(ia)
      return
      end


*     ***********************************
*     *                                 *
*     *            psp_comment          *
*     *                                 *
*     ***********************************
      character*(*) function psp_comment(ia)
      implicit none
      integer  ia

#include "psp.fh"

      psp_comment = comment(ia)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_lmmax		*
*     *					*
*     ***********************************
      integer function psp_lmmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_lmmax = int_mb(lmmax(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *            psp_nprj             *
*     *                                 *
*     ***********************************
      integer function psp_nprj(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_nprj  = int_mb(nprj(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *            psp_nprj_max         *
*     *                                 *
*     ***********************************
      integer function psp_nprj_max()
      implicit none

#include "bafdecls.fh"
#include "psp.fh"

*     **** local variables ****
      integer ia,nprjmax,nprjtmp


      nprjmax = 0
      do ia=1,npsp
         nprjtmp = (int_mb(nprj(1)+ia-1))
         if (nprjtmp.gt.nprjmax) nprjmax = nprjtmp
      end do

      psp_nprj_max  = nprjmax
      return
      end


*     ***********************************
*     *                                 *
*     *            psp_psp_type         *
*     *                                 *
*     ***********************************

      integer function psp_psp_type(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_psp_type  = int_mb(psp_type(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   psp_lmax		*
*     *					*
*     ***********************************
      integer function psp_lmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_lmax = int_mb(lmax(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_nmax		*
*     *					*
*     ***********************************
      integer function psp_nmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_nmax = int_mb(nmax(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_locp		*
*     *					*
*     ***********************************
      integer function psp_locp(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_locp = int_mb(locp(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_n1dgrid          *
*     *					*
*     ***********************************
      integer function psp_n1dgrid(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_n1dgrid = int_mb(n1dgrid(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_n1dbasis         *
*     *					*
*     ***********************************
      integer function psp_n1dbasis(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_n1dbasis = int_mb(n1dbasis(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_icut             *
*     *					*
*     ***********************************
      integer function psp_icut(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_icut = int_mb(icut(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_nae              *
*     *					*
*     ***********************************
      integer function psp_nae(i,ia)
      implicit none
      integer  i,ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_nae =int_mb(nae(1)+(ia-1)*(nmax_max*lmmax_max)+i-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_nps              *
*     *					*
*     ***********************************
      integer function psp_nps(i,ia)
      implicit none
      integer  i,ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_nps =int_mb(nps(1)+(ia-1)*(nmax_max*lmmax_max)+i-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_lps              *
*     *					*
*     ***********************************
      integer function psp_lps(i,ia)
      implicit none
      integer  i,ia

#include "bafdecls.fh"
#include "psp.fh"

      psp_lps =int_mb(lps(1)+(ia-1)*(nmax_max*lmmax_max)+i-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_eig              *
*     *					*
*     ***********************************
      real*8 function psp_eig(i,ia)
      implicit none
      integer  i,ia

#include "bafdecls.fh"
#include "psp.fh"

      integer  indx
      integer  psi_data_get_chnk
      external psi_data_get_chnk

      indx = psi_data_get_chnk(int_mb(eig(1)+ia-1))
      psp_eig = dbl_mb(indx+i-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   psp_npsp		*
*     *					*
*     ***********************************
      integer function psp_npsp()
      implicit none

#include "psp.fh"

      psp_npsp = npsp
      return
      end



*     ***********************************
*     *					*
*     *	    v_local_set_cosmo_BQ  	*
*     *					*
*     ***********************************
      subroutine v_local_set_cosmo_BQ(r_grid,rho,dng)
      implicit none
      real*8     r_grid(3,*)
      real*8     rho(*)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid,taskid_j,tid,np,np_j,nthreads
      integer nion,nion_q
      integer i,j,ii,ia,nx,ny,nz
      integer n2ft3d,n2ft3d_map,npack0
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3)
      real*8  c,x,y,z,r,dv,scal1,sqrt_pi
      integer BQ_ptr

      real*8 sa,sb 
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** Error function parameters ****
      real*8  yerf,verf,xerf

c     real*8 c1,c2,c3,c4,c5,c6,yerf,fterf,verf
      real*8 c1,c2,c3,c4,c5,c6,fterf
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)


*     **** external functions ****
      integer  ion_nion,ion_katm,control_version
      external ion_nion,ion_katm,control_version
      real*8   ion_rion,psp_rlocal,lattice_omega,psp_zv
      external ion_rion,psp_rlocal,lattice_omega,psp_zv
      integer  ion_nion_q
      external ion_nion_q
      logical  control_fast_erf
      external control_fast_erf
      real*8   util_erf
      external util_erf
      integer  nwpw_cosmo_BQ_ptr
      external nwpw_cosmo_BQ_ptr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  ion_nkatm
      external ion_nkatm

      call nwpw_timing_start(5)
      call Parallel_np(np)
      call Parallel2d_np_j(np_j)
      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      tid = Parallel_threadid()
      nthreads = Parallel_nthreads()
      call Pack_npack(0,npack0)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      call D3dB_n2ft3d(1,n2ft3d)
      nion     = ion_nion()
      nion_q   = ion_nion_q()
      BQ_ptr   = nwpw_cosmo_BQ_ptr()
      periodic = (control_version().eq.3)

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      dv    = scal1*lattice_omega()


      call Parallel_shared_vector_zero(.true.,nion_q,dbl_mb(BQ_ptr))

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) 
     >  call errquit('v_local_BQ:out of stack memory',0,MA_ERR)

      do ii=taskid_j+nion+1,nion+nion_q,np_j
    
c          if (mod(ii-1,np_j).eq.taskid_j) then

             if (.not.periodic) then
                rxyz(1) = ion_rion(1,ii)
                rxyz(2) = ion_rion(2,ii)
                rxyz(3) = ion_rion(3,ii)
                call lattice_r1_to_frac(1,rxyz,fxyz)
                inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                   (dabs(fxyz(2)).le.0.4d0).and.
     >                   (dabs(fxyz(3)).le.0.4d0))
             else
                inside = .true.
             end if

             if (inside) then
              !ia=ion_katm(ii)
              ia=ion_nkatm()


*             **** structure factor and local pseudopotential ****
              call strfac_pack(0,ii,dcpl_mb(exi(1)))
              call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                          dcpl_mb(exi(1)),
     >                          dcpl_mb(vtmp(1)))


              
              call Pack_cc_idot(0,dng,dcpl_mb(vtmp(1)),
     >                            dbl_mb(BQ_ptr+ii-nion-1))
!$OMP MASTER
              dbl_mb(BQ_ptr+ii-nion-1) = -dbl_mb(BQ_ptr+ii-nion-1)
!$OMP END MASTER
            
              !write(*,*) "HERA,ii=",ii,BQ_ptr,dbl_mb(BQ_ptr+ii-nion-1)

             end if

c          end if

      end do

      value = .true.
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('v_local_BQ:popping stack',0,MA_ERR)

      if (control_fast_erf()) then

      do j=taskid_j+nion+1,nion+nion_q,np_j

c         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            c = 1.0d0/psp_rlocal(ia)

!$OMP MASTER
            sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sb)
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                 yerf = (1.0d0 - 1.0d0/yerf**4)
                 sb = sb - rho(i)*(yerf/r)*dv
               else
                 sb = sb - 2.0d0*rho(i)*c/sqrt_pi*dv
               end if
            end do
!$OMP END DO
!$OMP MASTER
           dbl_mb(BQ_ptr+j-nion-1) = dbl_mb(BQ_ptr+j-nion-1) + sb
!$OMP END MASTER
!$OMP BARRIER

c         end if

      end do

      else

      do j=taskid_j+nion+1,nion+nion_q,np_j

c         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            c = 1.0d0/psp_rlocal(ia)

!$OMP MASTER
            sb = 0.0d0
!$OMP END MASTER
!$OMP DO REDUCTION(+:sb)
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = util_erf(xerf)
                 sb = sb - rho(i)*(yerf/r)*dv
               else
                 sb = sb - rho(i)*2.0d0*c/sqrt_pi*dv
               end if
            end do
!$OMP END DO
!$OMP MASTER
           dbl_mb(BQ_ptr+j-nion-1) = dbl_mb(BQ_ptr+j-nion-1) + sb
           !write(*,*) "HERB,j=",j,BQ_ptr,dbl_mb(BQ_ptr+j-nion-1)
!$OMP END MASTER
!$OMP BARRIER

c         end if

      end do

      end if

      !**** ion-ion BQ contribution ****
      do j=tid+nion+1,nion+nion_q,nthreads
         do ii=taskid+1,nion,np
            ia=ion_katm(ii)
            x = ion_rion(1,j)-ion_rion(1,ii)
            y = ion_rion(2,j)-ion_rion(2,ii)
            z = ion_rion(3,j)-ion_rion(3,ii)
            r = dsqrt(x*x + y*y +z*z)
            sb = psp_zv(ia)/r
           
           dbl_mb(BQ_ptr+j-nion-1) = dbl_mb(BQ_ptr+j-nion-1) + sb
         end do
         !write(*,*) "HERC,j=",j,BQ_ptr,dbl_mb(BQ_ptr+j-nion-1)
      end do

      call Parallel_Vector_SumAll(nion_q,dbl_mb(BQ_ptr))

      call nwpw_cosmo_solve_q()

      call nwpw_timing_end(5)
      return 
      end






*     **********************************
*     *                                *
*     *        grad_v_lr_local         *
*     *                                *
*     **********************************
*
*     This routine calculates the gradient of the long-range part of the
*     local pseudopotential (used by version 4)
*
*     Entry -
*     Exit -
*
      subroutine grad_v_lr_local(r_grid,rho,fion)
      implicit none
      real*8  r_grid(3,*)
      real*8  rho(*)
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** Error function parameters ****
      real*8  yerf,verf

c     real*8 c1,c2,c3,c4,c5,c6,yerf,fterf,verf
      real*8 c1,c2,c3,c4,c5,c6,fterf
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer ftmp(2)
      integer taskid_j,np_j
      integer i,j,ia,np1,np2,np3,n2ft3d_map,nion,l,m,lm,n2ft3d,nion_qm
      real*8 x,y,z,q,c,r,sqrt_pi,dv,v,rx,ry,rz,fx,fy,fz
      real*8 lmbda0,lmbda

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      logical  control_fast_erf
      integer  ion_nion,ion_nion_qm,ion_katm
      real*8   lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf
      real*8   pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external control_fast_erf
      external ion_nion,ion_nion_qm,ion_katm
      external lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf
      external pspw_qmmm_lambda
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      call D3dB_n2ft3d(1,n2ft3d)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

*     **** constants ****
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      dv = lattice_omega()/dble(np1*np2*np3)

*     ***** allocate temporary space ****
      if (.not.BA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1)))
     > call errquit('grad_v_lr_local:out of stack memory',0, MA_ERR)

c      call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)
      call Parallel_shared_vector_zero(.true.,3*nion,dbl_mb(ftmp(1)))
      if (psp_fmm) then
c         call dcopy((psp_fmm_lmax+1)**2,0.0d0,0,
c     >              dbl_mb(psp_fmm_Mlm(1)),1)
         call Parallel_shared_vector_zero(.true.,
     >                                    (psp_fmm_lmax+1)**2,
     >                                    dbl_mb(psp_fmm_Mlm(1)))
         !lm = 0
         !do l=0,psp_fmm_lmax
         !do m=-l,l
!$OMP DO
         do lm=0,(psp_fmm_lmax+1)**2-1
            !lm = l*l+l+m
            do i=1,n2ft3d_map
               dbl_mb(psp_fmm_Mlm(1)+lm) = dbl_mb(psp_fmm_Mlm(1)+lm) 
     >            + dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d+i-1)*rho(i)
            end do
         end do
!$OMP END DO
            !lm = lm + 1
         !end do
         !end do
         call D3dB_Vector_SumAll((psp_fmm_lmax+1)**2,
     >                           dbl_mb(psp_fmm_Mlm(1)))
         call FMM_fion_Mlm(psp_fmm_lmax,nion,nion_qm,
     >                     int_mb(ion_katm_ptr()),
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(psp_zv_ptr()),
     >                     lmbda0,
     >                     psp_fmm_rmax2,
     >                     dbl_mb(psp_fmm_Mlm(1)),
     >                     dbl_mb(ftmp(1)))

      end if

      if (control_fast_erf()) then

!$OMP DO
      do j=1,nion
         if (mod(j-1,np_j).eq.taskid_j) then
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            ia=ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            do i=1,n2ft3d_map
               rx = x - r_grid(1,i)
               ry = y - r_grid(2,i)
               rz = z - r_grid(3,i)
               r  = dsqrt( rx**2 + ry**2 + rz**2)

               if (r .gt. 1.0d-8) then
                 yerf=r*c
                 fterf = (1.0d0
     >                 + yerf*(c1 + yerf*(c2
     >                 + yerf*(c3 + yerf*(c4
     >                 + yerf*(c5 + yerf*c6))))))**4
                 verf = (1.0d0 - 1.0d0/fterf**4)
c                 verf = util_erf(yerf)
                 v    = q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                    - verf)/r**3
               else
                 v = 0.0d0
               end if

               fx = fx + rho(i)*rx*v
               fy = fy + rho(i)*ry*v
               fz = fz + rho(i)*rz*v
            end do
            dbl_mb(ftmp(1)+3*(j-1))   = -fx*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+1) = -fy*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+2) = -fz*dv*lmbda
            end if
         end if
      end do
!$OMP END DO

      else

!$OMP DO
      do j=1,nion
       
         if (mod(j-1,np_j).eq.taskid_j) then
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            do i=1,n2ft3d_map
               rx = x - r_grid(1,i)
               ry = y - r_grid(2,i)
               rz = z - r_grid(3,i)
               r  = dsqrt( rx**2 + ry**2 + rz**2)

               if (r .gt. 1.0d-8) then
                 yerf=r*c
c                fterf = (1.0d0
c    >                 + yerf*(c1 + yerf*(c2
c    >                 + yerf*(c3 + yerf*(c4
c    >                 + yerf*(c5 + yerf*c6))))))**4
c                verf = (1.0d0 - 1.0d0/fterf**4)
                 verf = util_erf(yerf)
                 v    = q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                    - verf)/r**3
               else
                 v = 0.0d0
               end if

               fx = fx + rho(i)*rx*v
               fy = fy + rho(i)*ry*v
               fz = fz + rho(i)*rz*v
            end do
            dbl_mb(ftmp(1)+3*(j-1))   = -fx*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+1) = -fy*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+2) = -fz*dv*lmbda

*        fion(1,j) = fion(1,j) - ddot(n2ft3d,rho,1,gv(1,1),3)*dv
*        fion(2,j) = fion(2,j) - ddot(n2ft3d,rho,1,gv(2,1),3)*dv
*        fion(3,j) = fion(3,j) - ddot(n2ft3d,rho,1,gv(3,1),3)*dv
c         call D3dB_SumAll(fx)
c         call D3dB_SumAll(fy)
c         call D3dB_SumAll(fz)
c         fion(1,j) = fion(1,j) - fx*dv
c         fion(2,j) = fion(2,j) - fy*dv
c         fion(3,j) = fion(3,j) - fz*dv

            end if
         end if
      end do
!$OMP END DO

      end if

      call Parallel_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
c      call daxpy(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)
!$OMP DO
      do j=1,nion
         fion(1,j) = fion(1,j) + dbl_mb(ftmp(1)+3*(j-1))
         fion(2,j) = fion(2,j) + dbl_mb(ftmp(1)+3*(j-1)+1)
         fion(3,j) = fion(3,j) + dbl_mb(ftmp(1)+3*(j-1)+2)
      end do
!$OMP END DO

      if (.not.BA_pop_stack(ftmp(2)))
     > call errquit('grad_v_lr_local:popping stack',1,MA_ERR)

      call nwpw_timing_end(5)

      return
      end





*     ***********************************
*     *				        *
*     *	 	   v_lr_local  		*
*     *					*
*     ***********************************
*
*     This routine calculates the long-range part of the
*     local pseudopotential (used by version4)
*
      subroutine v_lr_local(r_grid,vlr_out)
      implicit none
      real*8     r_grid(3,*)
      real*8     vlr_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** Error function parameters ****
      real*8 xerf,yerf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,xerf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer taskid_j,np_j
      integer i,j,ia,n2ft3d,n2ft3d_map,l,m,lm,nion,nion_qm,nion_q
      real*8 x,y,z,q,c,r,sqrt_pi,lmbda0,lmbda

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      logical  control_fast_erf
      integer  ion_nion,ion_nion_qm,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf,pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external control_fast_erf
      external ion_nion,ion_nion_qm,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf,pspw_qmmm_lambda
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      integer  ion_nion_q
      external ion_nion_q
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_q  = ion_nion_q()

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      !call dcopy(n2ft3d,0.0d0,0,vlr_out,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,vlr_out)

      if (psp_fmm) then
         call FMM_rion_Llm(psp_fmm_lmax,ion_nion(),nion_qm,
     >                     int_mb(ion_katm_ptr()),
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(psp_zv_ptr()),
     >                     lmbda0,
     >                     psp_fmm_rmax2,
     >                     dbl_mb(psp_fmm_Llm(1)))
!$OMP DO
         do i=1,n2ft3d_map
            lm = 0
            do l=0,psp_fmm_lmax
            do m=-l,l
               vlr_out(i) = vlr_out(i) 
     >                   - dbl_mb(psp_fmm_Llm(1)+lm)
     >                    *dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d+i-1)
               lm = lm + 1
            end do
            end do
         end do
!$OMP END DO
      end if

      if (control_fast_erf()) then

      do j=1,nion+nion_q
      
         if (mod(j-1,np_j).eq.taskid_j) then
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            if (j.gt.nion) then
               q = -nwpw_cosmo_screen_qc0(j-nion)
            else
               q = -psp_zv(ia)
            end if
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                 yerf = (1.0d0 - 1.0d0/yerf**4)
c                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + lmbda*(q/r)*yerf
               else
                 vlr_out(i) = vlr_out(i) + lmbda*2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if
      
      end do

      else

      do j=1,nion+nion_q
 
         if (mod(j-1,np_j).eq.taskid_j) then
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            if (j.gt.nion) then
               q = -nwpw_cosmo_screen_qc0(j-nion)
            else
               q = -psp_zv(ia)
            end if
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + lmbda*(q/r)*yerf
c                vlr_out(i) = vlr_out(i) + (q/r)*erf(r*c)
               else
                 vlr_out(i) = vlr_out(i) + lmbda*2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if

      end do

      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_out)

      call nwpw_timing_end(5)


      return
      end





*     ***********************************
*     *					*
*     *	 	   v_local  		*
*     *					*
*     ***********************************

      subroutine v_local(vl_out,move,dng,fion)
      implicit none
      complex*16 vl_out(*)
      logical    move
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm,nion_q
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3),lmbda,lmbda0

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      integer  Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      real*8   ion_rion,pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      external ion_rion,pspw_qmmm_lambda
      integer  ion_nion_q
      external ion_nion_q
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_q  = ion_nion_q()
      periodic = (control_version().eq.3)

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) 
     >  call errquit('v_local:out of stack memory',0,MA_ERR)

*     **** define Gx,Gy and Gz in packed space ****
      if (move) then
         value = BA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))
         if (.not. value) 
     >   call errquit('v_local: out of stack memory',0, MA_ERR)
         G(1)  = Pack_G_indx(0,1)
         G(2)  = Pack_G_indx(0,2)
         G(3)  = Pack_G_indx(0,3)
         !call dcopy(3*nion,0.0d0,0,fion,1)
         call Parallel_shared_vector_zero(.false.,3*nion,fion)
      end if

      !call dcopy((2*npack0),0.0d0,0,vl_out,1)
      call Parallel_shared_vector_zero(.true.,2*npack0,vl_out)
      do ii=1,nion+nion_q
    
          if (mod(ii-1,np_j).eq.taskid_j) then
             if (ii.gt.nion_qm) then
                lmbda = lmbda0
             else
                lmbda = 1.0d0
             end if

             if (.not.periodic) then
                rxyz(1) = ion_rion(1,ii)
                rxyz(2) = ion_rion(2,ii)
                rxyz(3) = ion_rion(3,ii)
                call lattice_r1_to_frac(1,rxyz,fxyz)
                inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                   (dabs(fxyz(2)).le.0.4d0).and.
     >                   (dabs(fxyz(3)).le.0.4d0))
             else
                inside = .true.
             end if

             if (inside) then
              ia=ion_katm(ii)

*             **** structure factor and local pseudopotential ****
              call strfac_pack(0,ii,dcpl_mb(exi(1)))
       
*             **** add to local psp ****
              if (ii.gt.nion) then
                 call Pack_tc_aMul(0,
     >                     dabs(nwpw_cosmo_screen_qc0(ii-nion)),
     >                             dbl_mb(vl(1)+npack0*(ia-1)),
     >                             dcpl_mb(exi(1)),
     >                             dcpl_mb(vtmp(1)))
              else if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                 call Pack_ttcc_AddMul(0,
     >                                 dbl_mb(vl(1)+npack0*(ia-1)),
     >                                 dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                                 dcpl_mb(exi(1)),
     >                                 dcpl_mb(vtmp(1)))
              else
                 call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))
              end if
              call Pack_c_SMul1(0,lmbda,dcpl_mb(vtmp(1)))
              call Pack_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)

              if ((move).and.(ii.le.nion)) then
              
c#ifndef CRAY
c!DIR$ ivdep
c#endif
c                do i=1,npack0
c                   dbl_mb(xtmp(1)+i-1) 
c     >             = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
c     >              - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
c                end do
               call Pack_cct_iconjgMulb(0,
     >                                  dng,
     >                                  dcpl_mb(vtmp(1)),
     >                                  dbl_mb(xtmp(1)))
               
               call Pack_tt_idot(0,dbl_mb(G(1)),
     >                             dbl_mb(xtmp(1)),fion(1,ii))
               call Pack_tt_idot(0,dbl_mb(G(2)),
     >                             dbl_mb(xtmp(1)),fion(2,ii))
               call Pack_tt_idot(0,dbl_mb(G(3)),
     >                             dbl_mb(xtmp(1)),fion(3,ii))
              end if
             end if

          end if

      end do

      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(2*npack0,vl_out)
      end if
      if (move) call Parallel_Vector_SumAll(3*nion,fion)

      value = .true.
      if (move) value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('v_local:popping stack',0,MA_ERR)

      call nwpw_timing_end(5)

      return 
      end






*     ***********************************
*     *				        *
*     *	 	   v_lr_local_mm 	*
*     *					*
*     ***********************************
*
*     This routine calculates the long-range part of the
*     local pseudopotential (used by version4)
*
      subroutine v_lr_local_mm(r_grid,vlr_out)
      implicit none
      real*8     r_grid(3,*)
      real*8     vlr_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** Error function parameters ****
      real*8 xerf,yerf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,xerf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer taskid_j,np_j
      integer i,j,ia,n2ft3d,n2ft3d_map,l,m,lm,nion_qm,nion
      real*8 x,y,z,q,c,r,sqrt_pi

*     **** external functions ****
      logical  control_fast_erf
      integer  ion_nion,ion_nion_qm,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf
      external control_fast_erf
      external ion_nion,ion_nion_qm,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      !call dcopy(n2ft3d,0.0d0,0,vlr_out,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,vlr_out)

      if (psp_fmm) then
         call FMM_rion_Llm_mm(psp_fmm_lmax,ion_nion(),nion_qm,
     >                     int_mb(ion_katm_ptr()),
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(psp_zv_ptr()),
     >                     psp_fmm_rmax2,
     >                     dbl_mb(psp_fmm_Llm(1)))
!$OMP DO
         do i=1,n2ft3d_map
            lm = 0
            do l=0,psp_fmm_lmax
            do m=-l,l
               vlr_out(i) = vlr_out(i) 
     >                   - dbl_mb(psp_fmm_Llm(1)+lm)
     >                    *dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d+i-1)
               lm = lm + 1
            end do
            end do
         end do
!$OMP END DO
      end if

      if (control_fast_erf()) then

      do j=nion_qm+1,nion
      
         if (mod(j-1,np_j).eq.taskid_j) then

            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                 yerf = (1.0d0 - 1.0d0/yerf**4)
c                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + (q/r)*yerf
               else
                 vlr_out(i) = vlr_out(i) + 2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if
      
      end do

      else

      do j=nion_qm+1,nion
 
         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + (q/r)*yerf
c                vlr_out(i) = vlr_out(i) + (q/r)*erf(r*c)
               else
                 vlr_out(i) = vlr_out(i) + 2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if

      end do

      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_out)

      call nwpw_timing_end(5)


      return
      end





*     ***********************************
*     *					*
*     *	 	   v_local_mm  		*
*     *					*
*     ***********************************

      subroutine v_local_mm(vl_out)
      implicit none
      complex*16 vl_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3)

*     **** external functions ****
      integer  ion_nion,ion_nion_qm,ion_katm,control_version
      real*8   ion_rion
      external ion_nion,ion_nion_qm,ion_katm,control_version
      external ion_rion

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      nion_qm  = ion_nion_qm()
      periodic = (control_version().eq.3)


      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) 
     >  call errquit('v_local:out of stack memory',0,MA_ERR)


      !call dcopy((2*npack0),0.0d0,0,vl_out,1)
      call Parallel_shared_vector_zero(.true.,2*npack0,vl_out)
      do ii=nion_qm+1,nion
    
          if (mod(ii-1,np_j).eq.taskid_j) then

             if (.not.periodic) then
                rxyz(1) = ion_rion(1,ii)
                rxyz(2) = ion_rion(2,ii)
                rxyz(3) = ion_rion(3,ii)
                call lattice_r1_to_frac(1,rxyz,fxyz)
                inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                   (dabs(fxyz(2)).le.0.4d0).and.
     >                   (dabs(fxyz(3)).le.0.4d0))
             else
                inside = .true.
             end if

             if (inside) then
              ia=ion_katm(ii)

*             **** structure factor and local pseudopotential ****
              call strfac_pack(0,ii,dcpl_mb(exi(1)))
       
*             **** add to local psp ****
              if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                 call Pack_ttcc_AddMul(0,
     >                                 dbl_mb(vl(1)+npack0*(ia-1)),
     >                                 dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                                 dcpl_mb(exi(1)),
     >                                 dcpl_mb(vtmp(1)))
              else
                 call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))
              end if
              call Pack_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)

             end if

          end if

      end do

      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(2*npack0,vl_out)
      end if

      value = .true.
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('v_local_mm:popping stack',0,MA_ERR)

      call nwpw_timing_end(5)

      return 
      end




*     ***********************************
*     *				        *
*     *	 	   v_lr_local_qm	*
*     *					*
*     ***********************************
*
*     This routine calculates the long-range part of the
*     local pseudopotential (used by version4)
*
      subroutine v_lr_local_qm(r_grid,vlr_out)
      implicit none
      real*8     r_grid(3,*)
      real*8     vlr_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** Error function parameters ****
      real*8 xerf,yerf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,xerf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer taskid_j,np_j
      integer i,j,ia,n2ft3d,n2ft3d_map,l,m,lm,nion_qm
      real*8 x,y,z,q,c,r,sqrt_pi

*     **** external functions ****
      logical  control_fast_erf
      integer  ion_nion,ion_nion_qm,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf
      external control_fast_erf
      external ion_nion,ion_nion_qm,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      nion_qm = ion_nion_qm()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      !call dcopy(n2ft3d,0.0d0,0,vlr_out,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,vlr_out)

      if (control_fast_erf()) then

      do j=1,nion_qm
      
         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                 yerf = (1.0d0 - 1.0d0/yerf**4)
c                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + (q/r)*yerf
               else
                 vlr_out(i) = vlr_out(i) + 2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if
      
      end do

      else

      do j=1,nion_qm
 
         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + (q/r)*yerf
c                vlr_out(i) = vlr_out(i) + (q/r)*erf(r*c)
               else
                 vlr_out(i) = vlr_out(i) + 2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if

      end do

      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_out)

      call nwpw_timing_end(5)

      return
      end

*     ***********************************
*     *					*
*     *	 	   v_local_qm  		*
*     *					*
*     ***********************************

      subroutine v_local_qm(vl_out)
      implicit none
      complex*16 vl_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3)

*     **** external functions ****
      integer  ion_nion,ion_nion_qm,ion_katm,control_version
      real*8   ion_rion
      external ion_nion,ion_nion_qm,ion_katm,control_version
      external ion_rion

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion_qm  = ion_nion_qm()
      periodic = (control_version().eq.3)

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) 
     >  call errquit('v_local_qm:out of stack memory',0,MA_ERR)

      !call dcopy((2*npack0),0.0d0,0,vl_out,1)
      call Parallel_shared_vector_zero(.true.,2*npack0,vl_out)
      do ii=1,nion_qm
    
          if (mod(ii-1,np_j).eq.taskid_j) then

             if (.not.periodic) then
                rxyz(1) = ion_rion(1,ii)
                rxyz(2) = ion_rion(2,ii)
                rxyz(3) = ion_rion(3,ii)
                call lattice_r1_to_frac(1,rxyz,fxyz)
                inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                   (dabs(fxyz(2)).le.0.4d0).and.
     >                   (dabs(fxyz(3)).le.0.4d0))
             else
                inside = .true.
             end if

             if (inside) then
              ia=ion_katm(ii)

*             **** structure factor and local pseudopotential ****
              call strfac_pack(0,ii,dcpl_mb(exi(1)))
       
*             **** add to local psp ****
              if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                 call Pack_ttcc_AddMul(0,
     >                                 dbl_mb(vl(1)+npack0*(ia-1)),
     >                                 dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                                 dcpl_mb(exi(1)),
     >                                 dcpl_mb(vtmp(1)))
              else
                 call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))
              end if
              call Pack_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)

             end if

          end if

      end do

      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(2*npack0,vl_out)
      end if

      value = .true.
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('v_local_qm:popping stack',0,MA_ERR)

      call nwpw_timing_end(5)

      return 
      end



*     ***********************************
*     *				        *
*     *	        v_lr_local_cosmo        *
*     *					*
*     ***********************************
*
*     This routine calculates the long-range part of the
*     local pseudopotential (used by version4)
*
      subroutine v_lr_local_cosmo(r_grid,vlr_out)
      implicit none
      real*8     r_grid(3,*)
      real*8     vlr_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** Error function parameters ****
      real*8 xerf,yerf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,xerf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer taskid_j,np_j
      integer i,j,ia,n2ft3d,n2ft3d_map,l,m,lm,nion,nion_qm,nion_q
      real*8 x,y,z,q,c,r,sqrt_pi,lmbda0,lmbda

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      logical  control_fast_erf
      integer  ion_nion,ion_nion_qm,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf,pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external control_fast_erf
      external ion_nion,ion_nion_qm,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf,pspw_qmmm_lambda
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      integer  ion_nion_q
      external ion_nion_q
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_q  = ion_nion_q()

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      !call dcopy(n2ft3d,0.0d0,0,vlr_out,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,vlr_out)

      if (psp_fmm) then
         call FMM_rion_Llm(psp_fmm_lmax,ion_nion(),nion_qm,
     >                     int_mb(ion_katm_ptr()),
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(psp_zv_ptr()),
     >                     lmbda0,
     >                     psp_fmm_rmax2,
     >                     dbl_mb(psp_fmm_Llm(1)))
!$OMP DO
         do i=1,n2ft3d_map
            lm = 0
            do l=0,psp_fmm_lmax
            do m=-l,l
               vlr_out(i) = vlr_out(i) 
     >                   - dbl_mb(psp_fmm_Llm(1)+lm)
     >                    *dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d+i-1)
               lm = lm + 1
            end do
            end do
         end do
!$OMP END DO
      end if

      if (control_fast_erf()) then

      do j=nion+1,nion+nion_q
      
         if (mod(j-1,np_j).eq.taskid_j) then
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            if (j.gt.nion) then
               q = -nwpw_cosmo_screen_qc0(j-nion)
            else
               q = -psp_zv(ia)
            end if
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                 yerf = (1.0d0 - 1.0d0/yerf**4)
c                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + lmbda*(q/r)*yerf
               else
                 vlr_out(i) = vlr_out(i) + lmbda*2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if
      
      end do

      else

      do j=nion+1,nion+nion_q
 
         if (mod(j-1,np_j).eq.taskid_j) then
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            if (j.gt.nion) then
               q = -nwpw_cosmo_screen_qc0(j-nion)
            else
               q = -psp_zv(ia)
            end if
            c = 1.0d0/psp_rlocal(ia)

!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                 xerf=r*c
                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + lmbda*(q/r)*yerf
c                vlr_out(i) = vlr_out(i) + (q/r)*erf(r*c)
               else
                 vlr_out(i) = vlr_out(i) + lmbda*2.0d0*q*c/sqrt_pi
               end if
            end do
!$OMP END DO
            end if

         end if

      end do

      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_out)

      call nwpw_timing_end(5)


      return
      end





*     ***********************************
*     *					*
*     *	 	   v_local_cosmo	*
*     *					*
*     ***********************************

      subroutine v_local_cosmo(vl_out)
      implicit none
      complex*16 vl_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm,nion_q
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3),lmbda,lmbda0

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      integer  Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      real*8   ion_rion,pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      external ion_rion,pspw_qmmm_lambda
      integer  ion_nion_q
      external ion_nion_q
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_q  = ion_nion_q()
      periodic = (control_version().eq.3)

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) 
     >  call errquit('v_local:out of stack memory',0,MA_ERR)


      !call dcopy((2*npack0),0.0d0,0,vl_out,1)
      call Parallel_shared_vector_zero(.true.,2*npack0,vl_out)
      do ii=nion+1,nion+nion_q
    
          if (mod(ii-1,np_j).eq.taskid_j) then
             if (ii.gt.nion_qm) then
                lmbda = lmbda0
             else
                lmbda = 1.0d0
             end if

             if (.not.periodic) then
                rxyz(1) = ion_rion(1,ii)
                rxyz(2) = ion_rion(2,ii)
                rxyz(3) = ion_rion(3,ii)
                call lattice_r1_to_frac(1,rxyz,fxyz)
                inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                   (dabs(fxyz(2)).le.0.4d0).and.
     >                   (dabs(fxyz(3)).le.0.4d0))
             else
                inside = .true.
             end if

             if (inside) then
              ia=ion_katm(ii)

*             **** structure factor and local pseudopotential ****
              call strfac_pack(0,ii,dcpl_mb(exi(1)))
       
*             **** add to local psp ****
              if (ii.gt.nion) then
                 call Pack_tc_aMul(0,
     >                     dabs(nwpw_cosmo_screen_qc0(ii-nion)),
     >                             dbl_mb(vl(1)+npack0*(ia-1)),
     >                             dcpl_mb(exi(1)),
     >                             dcpl_mb(vtmp(1)))
              else if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                 call Pack_ttcc_AddMul(0,
     >                                 dbl_mb(vl(1)+npack0*(ia-1)),
     >                                 dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                                 dcpl_mb(exi(1)),
     >                                 dcpl_mb(vtmp(1)))
              else
                 call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))
              end if
              call Pack_c_SMul1(0,lmbda,dcpl_mb(vtmp(1)))
              call Pack_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)
             end if

          end if

      end do

      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(2*npack0,vl_out)
      end if

      value = .true.
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('v_local:popping stack',0,MA_ERR)

      call nwpw_timing_end(5)

      return 
      end






*     ***********************************
*     *                                 *
*     *         v_local_addto_vcmp      *
*     *                                 *
*     ***********************************

      subroutine v_local_addto_vcmp(vcmp_out)
      implicit none
      complex*16 vcmp_out(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm
      integer i,ii,ia
      integer exi(2),vtmp(2)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3),scal1,lmbda0,lmbda

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      integer  Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      real*8   ion_rion,lattice_omega,pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      external ion_rion,lattice_omega,pspw_qmmm_lambda

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      nion_qm  = ion_nion_qm()
      periodic = (control_version().eq.3)
      scal1 = 1.0d0/lattice_omega()

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value)
     >  call errquit('v_local_addto_vcmp:out of stack memory',0,MA_ERR)

      call dcopy((2*npack0),0.0d0,0,dcpl_mb(vtmp(1)),1)
      do ii=1,nion
       if (mod(ii-1,np_j).eq.taskid_j) then
          if (ii.gt.nion_qm) then
             lmbda = lmbda0
          else
             lmbda = 1.0d0
          end if

          if (.not.periodic) then
             rxyz(1) = ion_rion(1,ii)
             rxyz(2) = ion_rion(2,ii)
             rxyz(3) = ion_rion(3,ii)
             call lattice_r1_to_frac(1,rxyz,fxyz)
             inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                (dabs(fxyz(2)).le.0.4d0).and.
     >                (dabs(fxyz(3)).le.0.4d0))
          else
             inside = .true.
          end if

          if (inside) then
           ia=ion_katm(ii)

*          **** structure factor and local pseudopotential ****
           call strfac_pack(0,ii,dcpl_mb(exi(1)))

*          **** add to local psp ****
c           call Pack_tc_MulAdd(0,dbl_mb(vl(1)+npack0*(ia-1)),
c     >                      dcpl_mb(exi(1)),
c     >                      dcpl_mb(vtmp(1)))
           call Pack_tc_aMulAdd(0,lmbda,dbl_mb(vl(1)+npack0*(ia-1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(vtmp(1)))
      
          end if

       end if
      end do
      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(2*npack0,dcpl_mb(vtmp(1)))
      end if
      call Pack_cc_daxpy(0,scal1,dcpl_mb(vtmp(1)),vcmp_out)


      value = .true. 
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value)
     >  call errquit('v_local_addto_vcmp:popping stack',1,MA_ERR)

      call nwpw_timing_end(5)
      return
      end





*     ***********************************
*     *					*
*     *	 	   f_vlocal  		*
*     *					*
*     ***********************************

      subroutine f_vlocal(dng,fion)
      implicit none
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      logical value,periodic,inside
      integer taskid_j,np_j
      integer npack0,nion,nion_qm
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
c      integer Gx(2),Gy(2),Gz(2)
      real*8 rxyz(3),fxyz(3),lmbda0,lmbda

*     **** external functions ****
      logical  pspw_qmmm_lambda_flag
      integer  Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      real*8   ion_rion,pspw_qmmm_lambda
      external pspw_qmmm_lambda_flag
      external Pack_G_indx,ion_nion,ion_nion_qm,ion_katm,control_version
      external ion_rion,pspw_qmmm_lambda

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      nion_qm  = ion_nion_qm()
      periodic = (control_version().eq.3)

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

c     **** define Gx,Gy and Gz in packed space ****
      G(1)  = Pack_G_indx(0,1)
      G(2)  = Pack_G_indx(0,2)
      G(3)  = Pack_G_indx(0,3)
      !call dcopy(3*nion,0.0d0,0,fion,1)
      call Parallel_shared_vector_zero(.true.,3*nion,fion)

      do ii=1,nion
       
         if (mod(ii-1,np_j).eq.taskid_j) then

          if (ii.gt.nion_qm) then
             lmbda = lmbda0
          else
             lmbda = 1.0d0
          end if

          if (.not.periodic) then
             rxyz(1) = ion_rion(1,ii)
             rxyz(2) = ion_rion(2,ii)
             rxyz(3) = ion_rion(3,ii)
             call lattice_r1_to_frac(1,rxyz,fxyz)
             inside =((dabs(fxyz(1)).le.0.4d0).and.
     >                (dabs(fxyz(2)).le.0.4d0).and.
     >                (dabs(fxyz(3)).le.0.4d0))
          else 
             inside = .true.
          endif

          if (inside) then
           ia=ion_katm(ii)

*          **** structure factor and local pseudopotential ****
           call strfac_pack(0,ii,dcpl_mb(exi(1)))

*          **** add to local psp ****
           if ((int_mb(psp_type(1)+ia-1).eq.4)) then
              call Pack_ttcc_AddMul(0,
     >                              dbl_mb(vl(1)+npack0*(ia-1)),
     >                              dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))
           else
              call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(vtmp(1)))
           end if
           call Pack_c_SMul1(0,lmbda,dcpl_mb(vtmp(1)))

c#ifndef CRAY
c!DIR$ ivdep
c#endif
c            do i=1,npack0
c              dbl_mb(xtmp(1)+i-1) 
c     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
c     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
c           end do
           call Pack_cct_iconjgMulb(0,
     >                              dng,
     >                              dcpl_mb(vtmp(1)),
     >                              dbl_mb(xtmp(1)))
           call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),fion(1,ii))
           call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),fion(2,ii))
           call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),fion(3,ii))
          end if

         end if

      end do
      call Parallel_Vector_SumAll(3*nion,fion)
      !call write_force(nion,fion,"f_vlocal_check")

      value =           BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))

      call nwpw_timing_end(5)
      return 
      end


*     ************************************************
*     *                                              *
*     *              v_nonlocal_paw_op1              *
*     *                                              *
*     ************************************************
*
*  paw operations #1 - generate it's compcharge, add atomic coulomb, and add atomic xc potential
*
      subroutine v_nonlocal_paw_op1(ispin,ne,nproj,ii,ia,
     >                              sw1,sw2,ww,eh_atom)
      implicit none 
      integer ispin,ne(2),nproj,ii,ia
      real*8 sw1(*),sw2(*),ww(*)
      real*8 eh_atom

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

      integer  psi_data_get_chnk
      external psi_data_get_chnk

      call nwpw_timing_start(39)

*     **** paw atom - generate it's atomic density matrix ****
      call psp_gen_density_matrix(ispin,ne,nproj,sw1,ww)

*     **** paw atom - generate it's compcharge ***
      call nwpw_compcharge_gen_Qlm(ii,ia,ispin,nproj,ww)


*     **** atomic coulomb matrix - sw2 = sw2 + Vhatomijl*sw1  ****
      call nwpw_compcharge_coulomb_atom(ii,ia,ispin,ne,nproj,
     >                             ww,sw1,sw2,eh_atom)

*     **** xc matrix - sw2 = sw2 + Vxcijl*sw1 ******
      call nwpw_xc_solve(ii,ia,
     >  int_mb(n1dgrid(1)+ia-1),
     >  int_mb(n1dbasis(1)+ia-1),
     >  dbl_mb(psi_data_get_chnk(int_mb(phi_ae(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(phi_ps(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(dphi_ae(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(dphi_ps(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(core_ae(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(core_ps(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(core_ae_prime(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(core_ps_prime(1)+ia-1))),
     >  dbl_mb(psi_data_get_chnk(int_mb(rgrid(1)+ia-1))),
     >  dbl_mb(log_amesh(1)+ia-1),
     >  ispin,ne,nproj,sw1,sw2)

      call nwpw_timing_end(39)
      return
      end


*     ***********************************
*     *					*
*     *	       v_nonlocal               *
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 12-1-2013 to handle PAW operators.
*
*  To Do -  For very large numbers of atoms the code will need to be restructured
*           to distribute the sw1a and sw2a matrices over np_i.  Basically, if the orthogonalization matrices
*           need to be distributed then sw1a and sw2a will need to be distributed as well.  A simple algorithm
*           to do this will be to keep the loop structure the same but instead of the
*           call to D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1))), this will need to be, 
*           changed to D3dB_Vector_SumAll(nn*nproj,dbl_mb(sw1(1))) and then place sw1 --> sw1a(distributed)
*
      subroutine v_nonlocal(ispin,ne,psi1,psi2,move,fion,
     >                      fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer pcount,taskid_j,np_j
      integer G(3),npack0,npack1,nion,nu,mult_l,m,ms,lm
      integer i,j,ii,ia,l,n,nn,l1,ip,iip,iipmax,jp,swstart,jpl
      integer k,shift,l_prj,m_prj,nproj,Gijl_indx
      integer nx,ny,nz
      real*8  omega,scal,ff(3),scal1,dv
      complex*16 ctmp
      integer exi(2),xtmp(2),sum(2)
      integer dng_cmp(2),dng_cmp_smooth(2)
      integer vcmp(2),vcmp_smooth(2)
      logical value,sd_function,periodic
      real*8 vmm(50),eh_atom
      integer ld_ptr,r,tid,nthr,Mtid,Itid

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega,ddot
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega,ddot
      integer  nwpw_compcharge_mult_l,control_version
      external nwpw_compcharge_mult_l,control_version
      real*8   nwpw_compcharge_Qlm
      external nwpw_compcharge_Qlm
      integer  Parallel_nthreads,Parallel_threadid
      external Parallel_nthreads,Parallel_threadid

      call nwpw_timing_start(6) 

      periodic = (control_version().eq.3)
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)

      call Pack_npack(1,npack1)
      nthr = Parallel_nthreads()
      tid  = Parallel_threadid()
      r    = mod(npack1,nthr)
      if (tid.lt.r) then
         Itid = tid*(npack1/nthr+1) + 1
         Mtid = npack1/nthr + 1
      else
          Itid = r + tid*(npack1/nthr) + 1
          Mtid = npack1/nthr
      end if


*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)

      call psp_prj_indx_alloc_sw1a_sw2a(nn)
      value = BA_push_get(mt_dcpl,npack1,'exi',exi(2), exi(1))
      if (.not.value) 
     >  call errquit('v_nonlocal:out of stack',0, MA_ERR)

      if (move) then
       value = value.and.BA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
       value = value.and.BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
       if (.not. value) 
     >  call errquit('v_nonlocal:out of stack',1,MA_ERR)

       G(1)  = Pack_G_indx(1,1)
       G(2)  = Pack_G_indx(1,2)
       G(3)  = Pack_G_indx(1,3)
      end if

      omega = lattice_omega()
      scal = 1.0d0/(omega)


      call nwpw_timing_start(56) 
      jp = 0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)

*        **** structure factor and local pseudopotential ****
         call strfac_pack(1,ii,dcpl_mb(exi(1)))

         jpl = jp
         call nwpw_timing_start(58) 
         do l=1,nproj
            shift       = int_mb(shift_prj_indx(1)+jp)
            sd_function = log_mb(sd_function_prj_indx(1)+jp)
            jp = jp + 1

*           **** phase factor does not matter therefore ****
*           **** (-i)^l is the same as (i)^l in the     ****
*           **** Rayleigh scattering formula            ****

*           **** phase fact DOES matter for compensation charge!!!!     ****
*           **** assume that sign factor for proj is in kbpp formatting ****

*           *** current function is s or d ****
            if (sd_function) then
               call Pack_tc_Mul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)+(l-1)*npack1))

*           *** current function is p or f ****
            else
               call Pack_tc_iMul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)+(l-1)*npack1))

            end if
c            call Pack_cc_indot(1,nn,
c     >                       psi1,
c     >                       dcpl_mb(prjtmp(1)),
c     >                       dbl_mb(sw1a(1)+(jp-1)*nn))
c            call Pack_cc_indot_thread(1,nn,
c     >                   psi1,
c     >                   dcpl_mb(prjtmp(1)),
c     >                   dbl_mb(sw1t(1)+(jp-1)*nn+tid*nn*n_prj_indx))
         end do
         call nwpw_timing_end(58) 
         call nwpw_timing_start(59) 
         call Pack_cc_indot_thread(1,nn,nproj,
     >                   psi1,
     >                   dcpl_mb(prjtmp(1)),
     >                   dbl_mb(sw1t(1)+jpl*nn+tid*nn*n_prj_indx))
         call nwpw_timing_end(59) 
      end do
!$OMP SINGLE
      call dcopy(nn*n_prj_indx,0.0d0,0,dbl_mb(sw1a(1)),1)
!$OMP END SINGLE
!$OMP CRITICAL
      call daxpy(nn*n_prj_indx,1.0d0,
     >           dbl_mb(sw1t(1)+tid*nn*n_prj_indx),1,
     >           dbl_mb(sw1a(1)),1)
!$OMP END CRITICAL
      call D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1)))
      call nwpw_timing_end(56) 


*     **** Compute sw2  ****
      eh_atom = 0.0d0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)
         swstart     = int_mb(swstart_prj_indx(1)+ip-1)


*        **** sw2 = Gijl*sw1 ******
         Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
         call Multiply_Gijl_sw1(nn,
     >                          nproj,
     >                          int_mb(nmax(1)+ia-1),
     >                          int_mb(lmax(1)+ia-1),
     >                          int_mb(n_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(l_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(m_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          dbl_mb(Gijl_indx),
     >                          dbl_mb(sw1a(1)+swstart*nn),
     >                          dbl_mb(sw2a(1)+swstart*nn))

     
*        **** paw operations #1 - generate it's compcharge, add atomic coulomb, ****
*        **** and add atomic xc potential                                       ****
         if ((int_mb(psp_type(1)+ia-1).eq.4)) then
            call v_nonlocal_paw_op1(ispin,ne,nproj,ii,ia,
     >                              dbl_mb(sw1a(1)+swstart*nn),
     >                              dbl_mb(sw2a(1)+swstart*nn),
     >                              dbl_mb(wtmp(1)),
     >                              eh_atom)
         end if
      end do



*     **** paw operations #2 - generate vcmp,.... and Gaussian Multipole  ****
      if (pawexist) then
         call nwpw_timing_start(39)

         if (use_grid_cmp) then
            call psp_cmp_dE_Qlm_pw(ispin,move,fion)
         else
            call psp_cmp_dE_Qlm(ispin,move,fion)
         end if

*        **** update ncmp*Vloc ****
         if (move) call psp_dE_ncmp_vloc_Qlm(ispin,move,fion)
         call nwpw_compcharge_add_dElocal_Qlm(ispin)

         do ip=1,nion_prj_indx
            ia = int_mb(ia_prj_indx(1)+ip-1)
            if ((int_mb(psp_type(1)+ia-1).eq.4)) then
               ii          = int_mb(ii_prj_indx(1)+ip-1)
               nproj       = int_mb(nproj_prj_indx(1)+ip-1)
               swstart     = int_mb(swstart_prj_indx(1)+ip-1)
               call nwpw_compcharge_gen_sw2(ii,ia,ispin,ne,nproj,
     >                                   dbl_mb(sw1a(1)+swstart*nn),
     >                                   dbl_mb(sw2a(1)+swstart*nn))
            end if
         end do

         call nwpw_timing_end(39)
      end if


*     **** do Kleinman-Bylander Multiplication ****
      call dscal_omp(nn*n_prj_indx,scal,dbl_mb(sw2a(1)),1)

*     **** apply the sw2 to psi ****
      jp  = 0
      do iip=1,nion_prj_indx,nprj_mult 

         swstart = int_mb(swstart_prj_indx(1)+iip-1)
         l1      = 0
         iipmax  = (iip+nprj_mult-1)
         if (iipmax.gt.nion_prj_indx) iipmax = nion_prj_indx

         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

*           **** structure factor and local pseudopotential ****
            call strfac_pack(1,ii,dcpl_mb(exi(1)))

            do l=1,nproj
               shift       = int_mb(shift_prj_indx(1)+jp)
               l_prj       = int_mb(l_prj_prj_indx(1)+jp)
               m_prj       = int_mb(m_prj_prj_indx(1)+jp)
               sd_function = log_mb(sd_function_prj_indx(1)+jp)
               jp = jp + 1

*              **** phase factor does not matter therefore ****
*              **** (-i)^l is the same as (i)^l in the     ****
*              **** Rayleigh scattering formula            ****

*              *** current function is s or d ****
               if (sd_function) then
                  call Pack_tc_Mul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))

*              *** current function is p or f ****
               else
                  call Pack_tc_iMul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))
               end if

*              ***** scale (sw2a) psp by factor - used for generating antiferromagnetic structures ****
*              **** nwchem input: pspspin up/down scale l ion_numbers                              ****
               if (pspspin) then
!$OMP MASTER
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_upl(1)+ii-1)) .and.
     >           ((m_prj.ne.int_mb(pspspin_upm(1)+ii-1)).or.
     >            (int_mb(pspspin_upm(1)+ii-1).gt.999))) 
     >            call dscal(ne(1),dbl_mb(pspspin_upscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(swstart+l1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_downl(1)+ii-1)).and.
     >           ((m_prj.ne.int_mb(pspspin_downm(1)+ii-1)).or.
     >            (int_mb(pspspin_downm(1)+ii-1).gt.999)))
     >            call dscal(ne(2),dbl_mb(pspspin_downscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(swstart+l1)*nn+ne(1)),1)
!$OMP END MASTER
!$OMP BARRIER
               end if

               l1 = l1 + 1
            end do
         end do

      call nwpw_timing_start(57) 
c         call DGEMM_OMP('N','T',2*npack1,nn,l1,
c     >             (-1.0d0),
c     >             dcpl_mb(prjtmp(1)), 2*npack1,
c     >             dbl_mb(sw2a(1)+swstart*nn),   nn,
c     >             (1.0d0),
c     >             psi2,               2*npack1)
         call DGEMM('N','T',2*Mtid,nn,l1,
     >             (-1.0d0),
     >             dcpl_mb(prjtmp(1)+Itid-1), 2*npack1,
     >             dbl_mb(sw2a(1)+swstart*nn),   nn,
     >             (1.0d0),
     >             psi2(Itid),               2*npack1)
      call nwpw_timing_end(57) 

         if (move) then
           l1 = 0
           do ip=iip,iipmax
              ii    = int_mb(ii_prj_indx(1)+ip-1)
              ia    = int_mb(ia_prj_indx(1)+ip-1)
              nproj = int_mb(nproj_prj_indx(1)+ip-1)

              do l=1,nproj
                 do n=1,nn
                    call Pack_cct_iconjgMul(1,
     >                                 dcpl_mb(prjtmp(1)+l1*npack1),
     >                                 psi1(1+(n-1)*npack1),
     >                                 dbl_mb(xtmp(1)))
                   call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                               dbl_mb(sum(1)+3*(n-1)))
                   call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                               dbl_mb(sum(1)+1+3*(n-1)))
                   call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                               dbl_mb(sum(1)+2+3*(n-1)))
                 end do
                 call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

                 !**** fractional weighting ****
                 if (fractional) then
!$OMP DO
                  do n=1,nn
                   call Dneall_qton(n,i)
                   dbl_mb(sum(1)+3*(n-1))  
     >                =dbl_mb(sum(1)  +3*(n-1))*occ(i)
                   dbl_mb(sum(1)+1+3*(n-1))
     >                =dbl_mb(sum(1)+1+3*(n-1))*occ(i)
                   dbl_mb(sum(1)+2+3*(n-1))
     >                =dbl_mb(sum(1)+2+3*(n-1))*occ(i)
                  end do
!$OMP END DO
                 end if
 
!OMP BARRIER

                 ff(1) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                                dbl_mb(sum(1)),3)
                 ff(2) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                                dbl_mb(sum(1)+1),3)
                 ff(3) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                                dbl_mb(sum(1)+2),3)

                 call D1dB_Vector_SumAll(3,ff)
!$OMP MASTER
                 fion(1,ii) = fion(1,ii)  + ff(1)*(3-ispin)
                 fion(2,ii) = fion(2,ii)  + ff(2)*(3-ispin)
                 fion(3,ii) = fion(3,ii)  + ff(3)*(3-ispin)
!$OMP END MASTER
!$OMP BARRIER

                 l1 = l1 + 1
              end do !** l **
           end do !** ip **
         end if !** move **

      end do

      value = .true.
      if (move) then

*        *** add multipole force here ***
         if (pawexist)
     >    call nwpw_compcharge_F_multipole_all(ispin,dbl_mb(zv(1)),fion)

         value = value.and.BA_pop_stack(sum(2))
         value = value.and.BA_pop_stack(xtmp(2))
      end if

      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('v_nonlocal: popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)
      return 
      end






      subroutine Multiply_Gijl_sw1(nn,nprj,nmax,lmax,
     >                             n_prj,l_prj,m_prj,
     >                             G,
     >                             sw1,sw2)
      implicit none
      integer nn
      integer nprj,nmax,lmax
      integer n_prj(nprj)
      integer l_prj(nprj)
      integer m_prj(nprj)
      real*8  G(nmax,nmax,0:lmax)
      real*8  sw1(nn,nprj)
      real*8  sw2(nn,nprj)

      !**** local variables ****
      integer a,b,na,nb,la,lb,ma,mb,ab
      

      !call dcopy(nn*nprj,0.0d0,0,sw2,1)
      call Parallel_shared_vector_zero(.true.,nn*nprj,sw2)
!$OMP DO
      do ab=0,(nprj*nprj-1)
         a = mod(ab,nprj)
         b = (ab-a)/nprj
         a = a+1
         b = b+1
         la = l_prj(a)
         lb = l_prj(b)
         ma = m_prj(a)
         mb = m_prj(b)
         if ((la.eq.lb).and.(ma.eq.mb)) then
            na = n_prj(a)
            nb = n_prj(b)
            call daxpy(nn,G(nb,na,la),sw1(1,a),1,sw2(1,b),1)
         end if
      end do
!$OMP END DO

c      do b=1,nprj
c         lb = l_prj(b)
c         mb = m_prj(b)
c
c         do a=1,nprj
c            la = l_prj(a)
c            ma = m_prj(a)
c
c            if ((la.eq.lb).and.(ma.eq.mb)) then
c              na = n_prj(a)
c              nb = n_prj(b)
c              call daxpy(nn,G(nb,na,la),sw1(1,a),1,sw2(1,b),1)
c            end if
c
c         end do
c      end do
      return
      end

*     ***********************************
*     *					*
*     *	 	   E_vnonlocal  	*
*     *					*
*     ***********************************
      real*8 function E_vnonlocal(ispin,ne,fractional,occ)
      implicit none
      integer ispin,ne(2)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     **** local variables ****
      real*8 E

      real*8 sum1,sum2
      common /eelectron_ejtmp/ sum1,sum2

*     **** external functions ****
      real*8   E_vnonlocal_sub
      external E_vnonlocal_sub

      E = E_vnonlocal_sub(ne(1)+ne(2),n_prj_indx,
     >                    dbl_mb(sw1a(1)),dbl_mb(sw2a(1)),
     >                    fractional,occ)
!$OMP MASTER
      sum1 = E
!$OMP END MASTER
!$OMP BARRIER
      call D1dB_SumAll(sum1)
      E = sum1
      if (ispin.eq.1) E = E+E

      E_vnonlocal = E
      return
      end

      real*8 function E_vnonlocal_sub(nn,nprjs,sw1,sw2,fractional,occ)
      implicit none
      integer nn,nprjs
      real*8 sw1(nn,nprjs),sw2(nn,nprjs)
      logical fractional
      real*8 occ(*)

*     **** local variables ****
      integer i,ip,n
      real*8 sum,sum2

      common /eelectron_ejtmp/ sum,sum2


!$OMP MASTER
      sum = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      if (fractional) then
!$OMP DO reduction(+:sum)
         do ip=1,nprjs
            do n=1,nn
               call Dneall_qton(n,i)
               sum = sum + sw1(n,ip)*sw2(n,ip)*occ(i)
            end do
         end do
!$OMP END DO
      else
!$OMP DO reduction(+:sum)
         do ip=1,nprjs
            do n=1,nn
               sum = sum + sw1(n,ip)*sw2(n,ip)
            end do
         end do
!$OMP END DO
      end if

      E_vnonlocal_sub = sum
      return
      end


*     ***********************************
*     *					*
*     *	       f_vnonlocal              *
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 12-1-2013 to handle PAW operators.
*
*  To Do -  For very large numbers of atoms the code will need to be restructured
*           to distribute the sw1a and sw2a matrices over np_i.  Basically, if the orthogonalization matrices
*           need to be distributed then sw1a and sw2a will need to be distributed as well.  A simple algorithm
*           to do this will be to keep the loop structure the same but instead of the
*           call to D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1))), this will need to be, 
*           changed to D3dB_Vector_SumAll(nn*nproj,dbl_mb(sw1(1))) and then place sw1 --> sw1a(distributed)
*
      subroutine f_vnonlocal(ispin,ne,psi1,fion,fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer pcount,taskid_j,np_j
      integer G(3),npack0,npack1,nion,nu,mult_l,m,ms,lm
      integer i,j,ii,ia,l,n,nn,l1,ip,iip,iipmax,jp,swstart
      integer k,shift,l_prj,m_prj,nproj,Gijl_indx
      integer nx,ny,nz
      real*8  omega,scal,ff(3),scal1,dv
      complex*16 ctmp
      integer exi(2),xtmp(2),sum(2)
      integer dng_cmp(2),dng_cmp_smooth(2)
      integer vcmp(2),vcmp_smooth(2)
      logical value,sd_function,periodic
      real*8 vmm(50),eh_atom
      integer ld_ptr
      !real*8  fion2(3,50),fion3(3,50)

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega,ddot
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega,ddot
      integer  nwpw_compcharge_mult_l,control_version
      external nwpw_compcharge_mult_l,control_version
      real*8   nwpw_compcharge_Qlm
      external nwpw_compcharge_Qlm


      call nwpw_timing_start(6) 

      periodic = (control_version().eq.3)
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      !call write_force(nion,fion,"f_vnonlocal 0, flocal?")
      !call dcopy(3*nion,fion,1,fion2,1)

      call psp_prj_indx_alloc_sw1a_sw2a(nn)
      value = BA_push_get(mt_dcpl,npack1,'exi',exi(2), exi(1))
      if (.not.value) 
     >  call errquit('v_nonlocal:out of stack',0, MA_ERR)

      value = value.and.BA_push_get(mt_dbl,npack1,
     >                              'xtmp',xtmp(2),xtmp(1))
      value = value.and.BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) 
     > call errquit('v_nonlocal:out of stack',1,MA_ERR)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      omega = lattice_omega()
      scal = 1.0d0/(omega)

      jp = 0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)

*        **** structure factor and local pseudopotential ****
         call strfac_pack(1,ii,dcpl_mb(exi(1)))

         do l=1,nproj
            shift       = int_mb(shift_prj_indx(1)+jp)
            sd_function = log_mb(sd_function_prj_indx(1)+jp)
            jp = jp + 1

*           **** phase factor does not matter therefore ****
*           **** (-i)^l is the same as (i)^l in the     ****
*           **** Rayleigh scattering formula            ****

*           **** phase fact DOES matter for compensation charge!!!!     ****
*           **** assume that sign factor for proj is in kbpp formatting ****

*           *** current function is s or d ****
            if (sd_function) then
               call Pack_tc_Mul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)))

*           *** current function is p or f ****
            else
               call Pack_tc_iMul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)))

            end if
            call Pack_cc_indot(1,nn,
     >                       psi1,
     >                       dcpl_mb(prjtmp(1)),
     >                       dbl_mb(sw1a(1)+(jp-1)*nn))
         end do
      end do
      call D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1)))


*     **** Compute sw2  ****
      eh_atom = 0.0d0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)
         swstart     = int_mb(swstart_prj_indx(1)+ip-1)


*        **** sw2 = Gijl*sw1 ******
         Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
         call Multiply_Gijl_sw1(nn,
     >                          nproj,
     >                          int_mb(nmax(1)+ia-1),
     >                          int_mb(lmax(1)+ia-1),
     >                          int_mb(n_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(l_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(m_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          dbl_mb(Gijl_indx),
     >                          dbl_mb(sw1a(1)+swstart*nn),
     >                          dbl_mb(sw2a(1)+swstart*nn))

     
*        **** paw operations #1 - generate it's compcharge, add atomic coulomb, and add atomic xc potential ****
         if ((int_mb(psp_type(1)+ia-1).eq.4)) then
            call v_nonlocal_paw_op1(ispin,ne,nproj,ii,ia,
     >                              dbl_mb(sw1a(1)+swstart*nn),
     >                              dbl_mb(sw2a(1)+swstart*nn),
     >                              dbl_mb(wtmp(1)),
     >                              eh_atom)
         end if
      end do



*     **** paw operations #2 - generate vcmp,.... and Gaussian Multipole  ****
      if (pawexist) then
         call nwpw_timing_start(39)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         dv = lattice_omega()*scal1

         if (periodic) then
            call Pack_npack(0,npack0)
         else
            call D3dB_n2ft3d(1,npack0)
            npack0 = npack0/2
         end if


         if (use_grid_cmp) then


            value = BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                          dng_cmp(2),dng_cmp(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'vcmp',
     >                          vcmp(2),vcmp(1))
            if (.not.value) 
     >      call errquit('f_vnonlocal:out of stack',4, MA_ERR)

            call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))

            call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp_tmp(1)))
            call Pack_cc_Sum(0,
     >                       dcpl_mb(vc_tmp(1)),
     >                       dcpl_mb(vcmp_tmp(1)),
     >                       dcpl_mb(vcmp(1)))
            call nwpw_compcharge_gen_dE_Qlm(1,ispin,
     >                                   dcpl_mb(vcmp(1)),
     >                                   dcpl_mb(vcmp(1)),.true.,fion)

            value =           BA_pop_stack(vcmp(2))
            value = value.and.BA_pop_stack(dng_cmp(2))
            if (.not.value) 
     >      call errquit('f_vnonlocal:popping stack',4, MA_ERR)

         else

            value = BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                          dng_cmp(2),dng_cmp(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                          dng_cmp_smooth(2),dng_cmp_smooth(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'vcmp',
     >                          vcmp(2),vcmp(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'vcmp_smooth',
     >                          vcmp_smooth(2),vcmp_smooth(1))
            if (.not.value) 
     >        call errquit('v_nonlocal:out of stack',4, MA_ERR)

            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dcpl_mb(dng_cmp(1)),
     >                                       dcpl_mb(dng_cmp_smooth(1)))

            !*** compute hartree potential of ntilde + ncmp_tilde ***
            !*** compute hartree potential ncmp   - ncmp_tilde ***
            if (periodic) then
               call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp(1)))
               call Pack_c_Copy(0,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))
               call coulomb_v(dcpl_mb(dng_cmp_smooth(1)),
     >                        dcpl_mb(vcmp_smooth(1)))
               call Pack_cc_Sub2(0,
     >                           dcpl_mb(vcmp_smooth(1)),
     >                           dcpl_mb(vcmp(1)))
               call Pack_cc_Sum2(0,dcpl_mb(vc_tmp(1)),
     >                             dcpl_mb(vcmp_smooth(1)))
               call Pack_c_SMul1(0,omega,dcpl_mb(vcmp(1)))
               call Pack_c_SMul1(0,omega,dcpl_mb(vcmp_smooth(1)))
            else
               !*** dng_cmp(G),dng_cmp_smooth(G) --> dng_cmp(r),dng_cmp_smooth(r) ***
               call Pack_c_unpack(0,dcpl_mb(dng_cmp(1)))
               call Pack_c_unpack(0,dcpl_mb(dng_cmp_smooth(1)))
               call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp(1)))
               call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp_smooth(1)))

               !*** generate vcmp and vcmp_smooth ***
               call coulomb2_v(dcpl_mb(dng_cmp(1)),
     >                         dcpl_mb(vcmp(1)))
               call coulomb2_v(dcpl_mb(dng_cmp_smooth(1)),
     >                         dcpl_mb(vcmp_smooth(1)))
               call D3dB_r_Copy(1,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))

               !*** vcmp        = vcmp-vcmp_smooth ***
               !*** vcmp_smooth = vcmp_smooth + vc ***
               call D3dB_rr_Sub2(1,
     >                           dcpl_mb(vcmp_smooth(1)),
     >                           dcpl_mb(vcmp(1)))
               call D3dB_rr_Sum2(1,
     >                           dcpl_mb(vc_tmp(1)),
     >                           dcpl_mb(vcmp_smooth(1)))

               !*** vcmp(r),vcmp_smooth(r) --> vcmp(G),vcmp_smooth(G) ***
               !*** May want to change nwpw_compcharge_gen_dE_Qlm to remove these FFTs ***
               call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp(1)))
               call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp_smooth(1)))
               call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp(1)))
               call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp_smooth(1)))
               call Pack_c_pack(0,dcpl_mb(vcmp(1)))
               call Pack_c_pack(0,dcpl_mb(vcmp_smooth(1)))

            end if

            call nwpw_compcharge_gen_dE_Qlm(3,ispin,
     >                                   dcpl_mb(vcmp_smooth(1)),
     >                                   dcpl_mb(vcmp(1)),.true.,fion)
         
            call nwpw_compcharge_gen_dEmult_Qlm(ispin)
            call nwpw_compcharge_add_dEmult_Qlm(ispin)

            value =           BA_pop_stack(vcmp_smooth(2))
            value = value.and.BA_pop_stack(vcmp(2))
            value = value.and.BA_pop_stack(dng_cmp_smooth(2))
            value = value.and.BA_pop_stack(dng_cmp(2))
            if (.not.value) 
     >        call errquit('v_nonlocal:popping stack',5,MA_ERR)
         end if

*        **** update ncmp*Vloc ****
         call psp_dE_ncmp_vloc_Qlm(ispin,.true.,fion)
         call nwpw_compcharge_add_dElocal_Qlm(ispin)

         do ip=1,nion_prj_indx
            ia = int_mb(ia_prj_indx(1)+ip-1)
            if ((int_mb(psp_type(1)+ia-1).eq.4)) then
               ii      = int_mb(ii_prj_indx(1)+ip-1)
               nproj   = int_mb(nproj_prj_indx(1)+ip-1)
               swstart = int_mb(swstart_prj_indx(1)+ip-1)
               call nwpw_compcharge_gen_sw2(ii,ia,ispin,ne,nproj,
     >                                   dbl_mb(sw1a(1)+swstart*nn),
     >                                   dbl_mb(sw2a(1)+swstart*nn))
            end if
         end do
     
         call nwpw_timing_end(39)
      end if


*     **** do Kleinman-Bylander Multiplication ****
      call DSCAL_OMP(nn*n_prj_indx,scal,dbl_mb(sw2a(1)),1)

*     **** apply the sw2 to psi ****
      jp  = 0
      do iip=1,nion_prj_indx,nprj_mult 

         swstart = int_mb(swstart_prj_indx(1)+iip-1)
         l1      = 0
         iipmax  = (iip+nprj_mult-1)
         if (iipmax.gt.nion_prj_indx) iipmax = nion_prj_indx

         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

*           **** structure factor and local pseudopotential ****
            call strfac_pack(1,ii,dcpl_mb(exi(1)))

            do l=1,nproj
               shift       = int_mb(shift_prj_indx(1)+jp)
               l_prj       = int_mb(l_prj_prj_indx(1)+jp)
               m_prj       = int_mb(m_prj_prj_indx(1)+jp)
               sd_function = log_mb(sd_function_prj_indx(1)+jp)
               jp = jp + 1

*              **** phase factor does not matter therefore ****
*              **** (-i)^l is the same as (i)^l in the     ****
*              **** Rayleigh scattering formula            ****

*              *** current function is s or d ****
               if (sd_function) then
                  call Pack_tc_Mul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))

*              *** current function is p or f ****
               else
                  call Pack_tc_iMul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))
               end if

*              ***** scale (sw2a) psp by factor - used for generating antiferromagnetic structures ****
*              **** nwchem input: pspspin up/down scale l ion_numbers                              ****
               if (pspspin) then
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_upl(1)+ii-1)) .and.
     >           ((m_prj.ne.int_mb(pspspin_upm(1)+ii-1)).or.
     >            (int_mb(pspspin_upm(1)+ii-1).gt.999))) 
     >            call dscal(ne(1),dbl_mb(pspspin_upscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(swstart+l1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_downl(1)+ii-1)).and.
     >           ((m_prj.ne.int_mb(pspspin_downm(1)+ii-1)).or.
     >            (int_mb(pspspin_downm(1)+ii-1).gt.999)))
     >            call dscal(ne(2),dbl_mb(pspspin_downscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(swstart+l1)*nn+ne(1)),1)
               end if

               l1 = l1 + 1
            end do
         end do

         l1 = 0
         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

            do l=1,nproj
               do n=1,nn
                  call Pack_cct_iconjgMul(1,
     >                               dcpl_mb(prjtmp(1)+l1*npack1),
     >                               psi1(1+(n-1)*npack1),
     >                               dbl_mb(xtmp(1)))
                 call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+3*(n-1)))
                 call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+1+3*(n-1)))
                 call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+2+3*(n-1)))

               end do
               call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

               !**** fractional weighting ****
               if (fractional) then
!$OMP DO
                do n=1,nn
                 call Dneall_qton(n,i)
                 dbl_mb(sum(1)+3*(n-1))  
     >              =dbl_mb(sum(1)  +3*(n-1))*occ(i)
                 dbl_mb(sum(1)+1+3*(n-1))
     >              =dbl_mb(sum(1)+1+3*(n-1))*occ(i)
                 dbl_mb(sum(1)+2+3*(n-1))
     >              =dbl_mb(sum(1)+2+3*(n-1))*occ(i)
                end do
!$OMP END DO
               end if
 
!$OMP MASTER
               ff(1) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)),3)
               ff(2) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)+1),3)
               ff(3) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)+2),3)
!$OMP END MASTER
!$OMP BARRIER
               call D1dB_Vector_SumAll(3,ff)
!$OMP MASTER
               fion(1,ii) = fion(1,ii)  + ff(1)*(3-ispin)
               fion(2,ii) = fion(2,ii)  + ff(2)*(3-ispin)
               fion(3,ii) = fion(3,ii)  + ff(3)*(3-ispin)
!$OMP END MASTER

               l1 = l1 + 1
            end do !** l **
         end do !** ip **

      end do !** ii **

*     *** add multipole force here ***
      if (pawexist)
     >  call nwpw_compcharge_F_multipole_all(ispin,dbl_mb(zv(1)),fion)

      value =           BA_pop_stack(sum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('v_nonlocal: popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)

c      call f_vnonlocal_check(ispin,ne,psi1,fractional,occ)
      return 
      end




*     ********************************************
*     *                                          *
*     *             psp_read                     *
*     *                                          *
*     ********************************************
      subroutine psp_read(fname, comment,
     >                       psp_type,
     >                       version,
     >                       nfft,unita,
     >                       atom,amass,zv,lmmax,lmax,locp,nmax,
     >                       rc,
     >                       nprj,n_projector,l_projector,m_projector,
     >                       b_projector,
     >                       Gijl_tag,
     >                       rlocal,
     >                       nfft3d,npack1,npack0,
     >                       vl,vlpaw,vnl_tag,
     >                       semicore,rcore,ncore,
     >                       log_amesh,r1,rmax,sigma,zion,
     >                       n1dgrid,n1dbasis,nae,nps,lps,icut,eig_tag,
     >                       phi_ae_tag,dphi_ae_tag,
     >                       phi_ps_tag,dphi_ps_tag,
     >                       core_ae_tag,core_ps_tag,
     >                       core_ae_prime_tag,core_ps_prime_tag,
     >                       rgrid_tag,
     >                       core_kin_energy,core_ion_energy,
     >                       hartree_tag,comp_charge_tag,comp_pot_tag,
     >                       r3_tag,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*50 fname
      character*(*) comment
      integer psp_type
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      character*2 atom
      real*8 amass,zv
      integer lmmax
      integer lmax
      integer locp
      integer nmax
      real*8 rc(*)

      integer nprj,n_projector(*),l_projector(*),m_projector(*)
      integer b_projector(*)
      integer Gijl_tag

      real*8 rlocal
      integer nfft3d,npack1,npack0
      real*8 vl(*),vlpaw(*)
      integer vnl_tag
      logical semicore
      real*8  rcore
      real*8  ncore(*)
      
      real*8  log_amesh,r1,rmax,sigma,zion
      integer n1dgrid,n1dbasis,nae(*),nps(*),lps(*),icut,eig_tag
      integer phi_ae_tag,dphi_ae_tag,phi_ps_tag,dphi_ps_tag
      integer core_ae_tag,core_ps_tag
      integer core_ae_prime_tag,core_ps_prime_tag
      integer rgrid_tag
      real*8  core_kin_energy,core_ion_energy
      integer hartree_tag,comp_charge_tag,comp_pot_tag
      integer r3_tag

      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#include "bafdecls.fh"
#include "util.fh"

*    *** local variables ***
      logical pio
      integer MASTER,taskid,taskid_i,taskid_p,com_p
      parameter(MASTER=0)
      integer n,l
      integer msglen
      integer iatom(2),Gijl_indx,indx
      character*255 full_filename
      real*8 sum1

*     **** external functions ****
      integer  psi_data_alloc,psi_data_get_chnk
      external psi_data_alloc,psi_data_get_chnk
      real*8   lattice_omega
      external lattice_omega
      logical  control_parallel_io
      external control_parallel_io

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      rlocal = 0.0d0
      n1dgrid = 0
      n1dbasis = 0
      sigma = 0.0d0

      pio = control_parallel_io()
      if (pio) then
         taskid_p = taskid_i
         com_p    = 1
      else
         taskid_p = taskid
         com_p    = 0
      end if


*     **** open fname binary file ****
      if (taskid_p.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)

         call cread(5,comment,80)
         call iread(5,psp_type,1)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call cread(5,atom,2)
         call dread(5,amass,1)
         call dread(5,zv,1)
         call iread(5,lmax,1)
         call iread(5,locp,1)
         call iread(5,nmax,1)
         lmmax=(lmax+1)**2 - (2*locp+1) !* number of projectors *
         amass = amass*1822.89d0
         call dread(5,rc,lmax+1)

         call iread(5,nprj,1)
         if (nprj.gt.0) then
         call iread(5,n_projector,nprj)
         call iread(5,l_projector,nprj)
         call iread(5,m_projector,nprj)
         call iread(5,b_projector,nprj)
         end if
      end if

*     **** send header data to all processors ****
      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,version)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,psp_type)
      msglen = 3
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,nfft)
      msglen = 9
      call Parallela_Brdcst_values(com_p,MASTER,msglen,unita)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,iatom)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call Parallela_Brdcst_values(com_p,MASTER,msglen,amass)
      msglen = 1
      call Parallela_Brdcst_values(com_p,MASTER,msglen,zv)
      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,lmax)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,locp)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,nmax)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,nprj)
      lmmax=(lmax+1)**2 - (2*locp+1)
      if (psp_type.eq.4) then
         n1dbasis = locp
      else
         n1dbasis = lmax+1
      end if


      msglen=lmax+1
      call Parallela_Brdcst_values(com_p,MASTER,msglen,rc)
      msglen = 1

      msglen=nprj
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,n_projector)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,l_projector)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,m_projector)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,b_projector)


*     ***** read in Gijl's ****
      if (nprj.gt.0) then
         if (psp_type.eq.4) then
           n = 5 !*** paw atom  ***
         else
           n = 1 !*** pspw atom ***
         end if
         msglen    = nmax*nmax*(lmax+1)
         Gijl_tag  = psi_data_alloc(n,msglen)
         Gijl_indx = psi_data_get_chnk(Gijl_tag)
         if (taskid_p.eq.MASTER) then
         call dread(5,dbl_mb(Gijl_indx),n*msglen)      !** number of matrix elements = n*nmax*nmax*(lmax+1) **
         end if
         call Parallela_Brdcst_values(com_p,MASTER,n*msglen,
     >                                dbl_mb(Gijl_indx))
      end if

*     ***** read in rlocal and rcore  ****
      if (taskid_p.eq.MASTER) then
         if (version.eq.4) call dread(5,rlocal,1)
         call dread(5,rcore,1)
      end if
      msglen=1
      call Parallela_Brdcst_values(com_p,MASTER,msglen,rlocal)
      msglen=1
      call Parallela_Brdcst_values(com_p,MASTER,msglen,rcore)


*     **** determine semicore value ****
      if (rcore.gt.0.0d0) then
         semicore = .true.
      else
         semicore = .false.
      end if

*     ***** Miscellaneous paw energies and 1d wavefunctions****
      if (psp_type.eq.4) then

         !**** allocate and reading hartree_matrix(n1dbasis,n1dbasis,n1dbasis,n1dbasis,2*lmax+1) ****
         msglen=n1dbasis*n1dbasis*n1dbasis*n1dbasis
         hartree_tag  = psi_data_alloc(2*lmax+1,msglen)
         indx = psi_data_get_chnk(hartree_tag)
         if (taskid_p.eq.MASTER) then
            call dread(5,dbl_mb(indx),(2*lmax+1)*msglen)      
         end if
         call Parallela_Brdcst_values(com_p,MASTER,(2*lmax+1)*msglen,
     >                                dbl_mb(indx))

         !**** allocate and reading comp_charge_matrix(n1dbasis,n1dbasis,2*lmax+1) ****
         msglen=n1dbasis*n1dbasis
         comp_charge_tag  = psi_data_alloc(2*lmax+1,msglen)
         indx = psi_data_get_chnk(comp_charge_tag)
         if (taskid_p.eq.MASTER) then
            call dread(5,dbl_mb(indx),(2*lmax+1)*msglen)      
         end if
         call Parallela_Brdcst_values(com_p,MASTER,(2*lmax+1)*msglen,
     >                                dbl_mb(indx))

         !**** allocate and reading comp_pot_matrix(n1dbasis,n1dbasis,2*lmax+1) ****
         msglen=n1dbasis*n1dbasis
         comp_pot_tag  = psi_data_alloc(2*lmax+1,msglen)
         indx = psi_data_get_chnk(comp_pot_tag)
         if (taskid_p.eq.MASTER) then
            call dread(5,dbl_mb(indx),(2*lmax+1)*msglen)      
         end if
         call Parallela_Brdcst_values(com_p,MASTER,(2*lmax+1)*msglen,
     >                                dbl_mb(indx))


        if (taskid_p.eq.MASTER) call dread(5,core_kin_energy,1)
        if (taskid_p.eq.MASTER) call dread(5,core_ion_energy,1)
        msglen=1
        call Parallela_Brdcst_values(com_p,MASTER,msglen,
     >                               core_kin_energy)
        call Parallela_Brdcst_values(com_p,MASTER,msglen,
     >                               core_ion_energy)

        if (taskid_p.eq.MASTER) then
           call iread(5,n1dgrid,1)
           call iread(5,icut,1)
           call dread(5,log_amesh,1)
           call dread(5,r1,1)
           call dread(5,rmax,1)
           call dread(5,sigma,1)
           call dread(5,zion,1)
        end if
        msglen=1
        call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,n1dgrid)
        call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,icut)
        call Parallela_Brdcst_values(com_p,MASTER,msglen,log_amesh)
        call Parallela_Brdcst_values(com_p,MASTER,msglen,r1)
        call Parallela_Brdcst_values(com_p,MASTER,msglen,rmax)
        call Parallela_Brdcst_values(com_p,MASTER,msglen,sigma)
        call Parallela_Brdcst_values(com_p,MASTER,msglen,zion)


        eig_tag     = psi_data_alloc(n1dbasis,1)
        phi_ae_tag  = psi_data_alloc(n1dbasis,n1dgrid)
        dphi_ae_tag = psi_data_alloc(n1dbasis,n1dgrid)
        phi_ps_tag  = psi_data_alloc(n1dbasis,n1dgrid)
        dphi_ps_tag = psi_data_alloc(n1dbasis,n1dgrid)
        core_ae_tag = psi_data_alloc(1,n1dgrid)
        core_ps_tag = psi_data_alloc(1,n1dgrid)
        core_ae_prime_tag = psi_data_alloc(1,n1dgrid)
        core_ps_prime_tag = psi_data_alloc(1,n1dgrid)
        rgrid_tag   = psi_data_alloc(1,n1dgrid)
        
        indx = psi_data_get_chnk(eig_tag)
        if (taskid_p.eq.MASTER) call dread(5,dbl_mb(indx),n1dbasis)
        call Parallela_Brdcst_values(com_p,MASTER,n1dbasis,dbl_mb(indx))

        if (taskid_p.eq.MASTER) call iread(5,nae,n1dbasis)
        call Parallela_Brdcst_ivalues(com_p,MASTER,n1dbasis,nae)

        if (taskid_p.eq.MASTER) call iread(5,nps,n1dbasis)
        call Parallela_Brdcst_ivalues(com_p,MASTER,n1dbasis,nps)

        if (taskid_p.eq.MASTER) call iread(5,lps,n1dbasis)
        call Parallela_Brdcst_ivalues(com_p,MASTER,n1dbasis,lps)

        indx = psi_data_get_chnk(rgrid_tag)
        if (taskid_p.eq.MASTER) call dread(5,dbl_mb(indx),n1dgrid)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid,dbl_mb(indx))

        indx = psi_data_get_chnk(phi_ae_tag)
        if(taskid_p.eq.MASTER) 
     >     call dread(5,dbl_mb(indx),n1dgrid*n1dbasis)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid*n1dbasis,
     >                              dbl_mb(indx))

        indx = psi_data_get_chnk(dphi_ae_tag)
        if(taskid_p.eq.MASTER) 
     >     call dread(5,dbl_mb(indx),n1dgrid*n1dbasis)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid*n1dbasis,
     >                              dbl_mb(indx))

        indx = psi_data_get_chnk(phi_ps_tag)
        if(taskid_p.eq.MASTER) 
     >     call dread(5,dbl_mb(indx),n1dgrid*n1dbasis)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid*n1dbasis,
     >                              dbl_mb(indx))

        indx = psi_data_get_chnk(dphi_ps_tag)
        if(taskid_p.eq.MASTER) 
     >     call dread(5,dbl_mb(indx),n1dgrid*n1dbasis)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid*n1dbasis,
     >                              dbl_mb(indx))

        indx = psi_data_get_chnk(core_ae_tag)
        if (taskid_p.eq.MASTER) call dread(5,dbl_mb(indx),n1dgrid)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid,dbl_mb(indx))

        indx = psi_data_get_chnk(core_ps_tag)
        if (taskid_p.eq.MASTER) call dread(5,dbl_mb(indx),n1dgrid)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid,dbl_mb(indx))

        indx = psi_data_get_chnk(core_ae_prime_tag)
        if (taskid_p.eq.MASTER) call dread(5,dbl_mb(indx),n1dgrid)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid,dbl_mb(indx))

        indx = psi_data_get_chnk(core_ps_prime_tag)
        if (taskid_p.eq.MASTER) call dread(5,dbl_mb(indx),n1dgrid)
        call Parallela_Brdcst_values(com_p,MASTER,n1dgrid,dbl_mb(indx))

      end if

*     *** read in vl 3d block ***
      if (pio) then
         call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
      else
         call D3dB_t_read(1,5,tmp2,tmp,-1)
      end if
      call Pack_t_pack(0,tmp2)
      call Pack_t_Copy(0,tmp2,vl)

*     **** read in vlpaw 3d block - extra paw local potentials ****
      if (psp_type.eq.4) then
         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,vlpaw)
      end if


*     **** read in vnl 3d blocks ****
      if (nprj.gt.0) then
         vnl_tag = psi_data_alloc(nprj,npack1)

         do n=1,nprj
            if (pio) then
               call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
            else
               call D3dB_t_read(1,5,tmp2,tmp,-1)
            end if
            call Pack_t_pack(1,tmp2)
            !call Pack_t_Copy(1,tmp2,vnl(1,n))
            call psi_data_add(vnl_tag,n,tmp2)
         end do
      end if

*     **** read in semicore density block ****
      if (semicore) then
         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1))

         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+2*npack0))

         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if 
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+3*npack0))

         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+4*npack0))
      end if

*     **** read in r3_matrix if not paw(for now) ****
      if (psp_type.eq.9) then

         !**** allocate and reading r3_matrix(n1dbasis,n1dbasis) ****
         msglen=n1dbasis*n1dbasis
         r3_tag = psi_data_alloc(1,msglen)
         indx   = psi_data_get_chnk(r3_tag)
         if (taskid_p.eq.MASTER) then
            call dread(5,dbl_mb(indx),msglen)      
         end if
         call Parallela_Brdcst_values(com_p,MASTER,msglen,
     >                                dbl_mb(indx))

      end if

*     *** close fname binary file ***
      if (taskid_p.eq.MASTER) then
c       close(11)
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *					*
*     *	 	  psp_readall  		*
*     *					*
*     ***********************************

      subroutine psp_readall()
      implicit none
      
#include "bafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"

#include "psp.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore


*     **** local variables ****
      integer ngp(3),version,nfft3d,npack1,npack0
      integer ia,l
      real*8 unita(3,3)
      character*12 boundry
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box
      character*5  element
      character*50 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind,control_print
      integer      control_ngrid,ion_nion
      real*8       control_unita
      character*12 control_boundry
      character*4  ion_atom
      external     nwpw_filefind,control_print
      external     control_ngrid,ion_nion
      external     control_unita
      external     control_boundry
      external     ion_atom
      
      call nwpw_timing_start(50)

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call Parallel_taskid(taskid)

*     *** set semicore(0) *****
      log_mb(semicore(1)) = .false.

      value = BA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** read pseudopotentials ****
      do ia=1,npsp

*      **** define formatted psp name ****
       element = '     '
       element = ion_atom(ia)
       l = index(element,' ') - 1
       fname = element(1:l)//'.vpp'
        

       found = .false.
       do while (.not.found)

         if (nwpw_filefind(fname)) then
            call psp_read(fname, comment(ia),
     >                  int_mb(psp_type(1)+ia-1),
     >                  version,
     >                  ngp,unita,
     >                  atom(ia),
     >                  dbl_mb(amass(1)+ia-1),
     >                  dbl_mb(zv(1)+ia-1),
     >                  int_mb(lmmax(1)+ia-1),
     >                  int_mb(lmax(1)+ia-1),
     >                  int_mb(locp(1)+ia-1),
     >                  int_mb(nmax(1)+ia-1),
     >                  dbl_mb(rc(1) + (ia-1)*(lmax_max+1)),
     >                  int_mb(nprj(1)+ia-1),
     >                  int_mb(n_projector(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(l_projector(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(m_projector(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(b_projector(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(Gijl(1)+(ia-1)), 
     >                  dbl_mb(rlocal(1) + (ia-1)),
     >                  nfft3d,npack1,npack0,
     >                  dbl_mb(vl(1) + (ia-1)*npack0),
     >                  dbl_mb(vlpaw(1) + (ia-1)*npack0),
     >                  int_mb(vnl(1)+ (ia-1)),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(rcore(1)+ia-1),
     >                  dbl_mb(ncore(1)+ (ia-1)*npack0*5),
     >                  dbl_mb(log_amesh(1)+ia-1),
     >                  dbl_mb(r1(1)+ia-1),
     >                  dbl_mb(rmax(1)+ia-1),
     >                  dbl_mb(sigma(1)+ia-1),
     >                  dbl_mb(zion(1)+ia-1),
     >                  int_mb(n1dgrid(1)+ia-1),
     >                  int_mb(n1dbasis(1)+ia-1),
     >                  int_mb(nae(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(nps(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(lps(1) 
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(icut(1)+ia-1),
     >                  int_mb(eig(1)+ia-1),
     >                  int_mb(phi_ae(1)+ia-1),
     >                  int_mb(dphi_ae(1)+ia-1),
     >                  int_mb(phi_ps(1)+ia-1),
     >                  int_mb(dphi_ps(1)+ia-1),
     >                  int_mb(core_ae(1)+ia-1),
     >                  int_mb(core_ps(1)+ia-1),
     >                  int_mb(core_ae_prime(1)+ia-1),
     >                  int_mb(core_ps_prime(1)+ia-1),
     >                  int_mb(rgrid(1)+ia-1),
     >                  dbl_mb(core_kin(1)+ia-1),
     >                  dbl_mb(core_ion(1)+ia-1),
     >                  int_mb(hartree_matrix(1)+ia-1),
     >                  int_mb(comp_charge_matrix(1)+ia-1),
     >                  int_mb(comp_pot_matrix(1)+ia-1),
     >                  int_mb(r3_matrix(1)+ia-1),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)


*          **** set semicore(0) ****
           if (log_mb(semicore(1)+ia)) log_mb(semicore(1)) = .true.
           if (ierr.gt.0) go to 9000

*          **** set pawexist ****
           if (int_mb(psp_type(1)+ia-1).eq.4) pawexist = .true.

*          **************************************************************
*          ***** logic for finding out if psp is correctly formatted ****
*          **************************************************************
           correct_box = .true.
           boundry = control_boundry()
           l =index(boundry,' ') - 1
           if ( (ngp(1).ne.control_ngrid(1)) .or.
     >       (ngp(2).ne.control_ngrid(2)) .or.
     >       (ngp(3).ne.control_ngrid(3)) .or. 
     >       (unita(1,1).ne.control_unita(1,1)) .or.
     >       (unita(2,1).ne.control_unita(2,1)) .or.
     >       (unita(3,1).ne.control_unita(3,1)) .or.
     >       (unita(1,2).ne.control_unita(1,2)) .or.
     >       (unita(2,2).ne.control_unita(2,2)) .or.
     >       (unita(3,2).ne.control_unita(3,2)) .or.
     >       (unita(1,3).ne.control_unita(1,3)) .or.
     >       (unita(2,3).ne.control_unita(2,3)) .or.
     >       (unita(3,3).ne.control_unita(3,3)) .or.
     >       ((boundry(1:l).eq.'periodic').and.(version.ne.3)).or.
     >       ((boundry(1:l).eq.'aperiodic').and.(version.ne.4))) then
              correct_box = .false.
              if ((taskid.eq.MASTER).and.
     >           control_print(print_medium)) then
              write(luout,*) 
     >         "pseudopotential is not correctly formatted:",fname
              end if


*             *** deallocate memory ***
              if (int_mb(nprj(1)+ia-1).gt.0) then
              call psi_data_dealloc(int_mb(vnl(1)+ia-1))
              call psi_data_dealloc(int_mb(Gijl(1)+ia-1))
               if (int_mb(psp_type(1)+ia-1).eq.4) then
               call psi_data_dealloc(int_mb(hartree_matrix(1)+ia-1))
               call psi_data_dealloc(int_mb(comp_charge_matrix(1)+ia-1))
               call psi_data_dealloc(int_mb(comp_pot_matrix(1)+ia-1))
               call psi_data_dealloc(int_mb(eig(1)+ia-1))
               call psi_data_dealloc(int_mb(phi_ae(1)+ia-1))
               call psi_data_dealloc(int_mb(dphi_ae(1)+ia-1))
               call psi_data_dealloc(int_mb(phi_ps(1)+ia-1))
               call psi_data_dealloc(int_mb(dphi_ps(1)+ia-1))
               call psi_data_dealloc(int_mb(core_ae(1)+ia-1))
               call psi_data_dealloc(int_mb(core_ps(1)+ia-1))
               call psi_data_dealloc(int_mb(core_ae_prime(1)+ia-1))
               call psi_data_dealloc(int_mb(core_ps_prime(1)+ia-1))
               call psi_data_dealloc(int_mb(rgrid(1)+ia-1))
               end if
               if (int_mb(psp_type(1)+ia-1).eq.9) then
               call psi_data_dealloc(int_mb(r3_matrix(1)+ia-1))
               end if
              end if


           end if
           if (correct_box) found = .true.

         end if

*        **** generate formatted pseudopotential atom.vpp *****
         if (.not.found) then
             call psp_formatter_auto(ion_atom(ia))
         end if

       end do !***do while ****


      end do
 9000 value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not. value)
     > call errquit('psp_readall:error popping stack',0,MA_ERR)



*     **** done reading set nprj_max and prjtmp for nonlocal psp operator ****
      call psp_proj_init()

*     **** initialize prj_indx, used by v_nonlocal_new ****
      call psp_prj_indx_init()


      if (pawexist) then
         call psp_paw_init()
         call nwpw_compcharge_init(ion_nion(),npsp,
     >                             int_mb(nprj(1)),
     >                             int_mb(n1dbasis(1)),
     >                             int_mb(psp_type(1)),
     >                             int_mb(lmax(1)),
     >                             dbl_mb(sigma(1)),
     >                             nmax_max*lmmax_max,
     >                             int_mb(l_projector(1)),
     >                             int_mb(m_projector(1)),
     >                             int_mb(b_projector(1)),
     >                             int_mb(comp_charge_matrix(1)),
     >                             int_mb(hartree_matrix(1)))

         call nwpw_xc_init(ion_nion(),npsp,
     >                     int_mb(nprj(1)),
     >                     int_mb(n1dbasis(1)),
     >                     int_mb(n1dgrid(1)),
     >                     int_mb(psp_type(1)),
     >                     int_mb(lmax(1)),
     >                     nmax_max*lmmax_max,
     >                     int_mb(l_projector(1)),
     >                     int_mb(m_projector(1)),
     >                     int_mb(b_projector(1)))

      end if


      call nwpw_timing_end(50)
      return
      end


*     ***********************************
*     *					*
*     *	 	  psp_print 		*
*     *					*
*     ***********************************
      subroutine psp_print(ia)      
      implicit none
      integer ia

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "psp.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
      logical oprint
      integer i
      character*255 cmment

      logical       control_print,psp_semicore
      external      control_print,psp_semicore
      character     spdf_name
      external      spdf_name
      character*4   ion_atom
      external      ion_atom
      real*8        psp_rc,psp_rlocal,psp_rcore,psp_ncore,psp_zv
      external      psp_rc,psp_rlocal,psp_rcore,psp_ncore,psp_zv
      real*8        psp_zion,psp_sigma,psp_sphere_radius
      external      psp_zion,psp_sigma,psp_sphere_radius
      real*8        psp_r1,psp_rmax,psp_eig
      external      psp_r1,psp_rmax,psp_eig
      integer       psp_lmax,psp_locp,psp_psp_type,psp_nprj,psp_icut
      external      psp_lmax,psp_locp,psp_psp_type,psp_nprj,psp_icut
      integer       psp_nae,psp_nps,psp_lps,psp_n1dbasis,psp_n1dgrid
      external      psp_nae,psp_nps,psp_lps,psp_n1dbasis,psp_n1dgrid
      integer       control_version,inp_strlen
      external      control_version,inp_strlen
      
      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))

      if (oprint) then
         if (psp_psp_type(ia).ne.4) then
           write(luout,1150) ia,ion_atom(ia),psp_zv(ia),psp_lmax(ia)
           cmment = comment(ia)
           i = inp_strlen(cmment)
           write(luout,1157) cmment(1:i)
           write(luout,1158) psp_psp_type(ia)
           write(luout,1152) psp_lmax(ia)
           write(luout,1153) psp_locp(ia)
           write(luout,1154) psp_nprj(ia)
           if (control_version().eq.4) write(luout,1156) psp_rlocal(ia)
           if (psp_semicore(ia)) 
     >        write(luout,1155) psp_rcore(ia),psp_ncore(ia)
           write(luout,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
           write(luout,*)
         else 
            cmment = comment(ia)
            i = inp_strlen(cmment)
            write(luout,2141) ia,ion_atom(ia),psp_zv(ia),
     >                        psp_zion(ia)-psp_zv(ia)
            write(luout,1157) cmment(1:i)
            write(luout,1158) psp_psp_type(ia)
            write(luout,2142) psp_r1(ia)
            write(luout,2143) psp_rmax(ia)
            write(luout,2144) psp_n1dgrid(ia)
            write(luout,2145) psp_sphere_radius(ia),
     >                        psp_icut(ia),
     >                        psp_icut(ia)
            write(luout,2146) psp_sigma(ia)
            write(luout,2150) psp_nprj(ia)
            write(luout,2151)
            do i=1,psp_n1dbasis(ia)
              write(luout,2152) psp_nps(i,ia),
     >                          psp_nae(i,ia),
     >                          spdf_name(psp_lps(i,ia)),
     >                          psp_eig(i,ia),
     >                          2*psp_lps(i,ia)+1
            end do
            write(luout,*)
         end if
      end if

 1150 FORMAT(5X,I2,': ',A4,'valence charge:',F8.4,'  lmax=',I3)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i3)
 1153 FORMAT(12X,' local potential used           : ',i3)
 1154 FORMAT(12X,' number of non-local projections: ',i3)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1157 FORMAT(12X,' comment    : ',A)
 1158 FORMAT(12X,' pseudpotential type            : ',i3)

 2141 FORMAT(5X,I2,': ',A4,'valence charge:',F4.1,'  core charge:',F4.1)
 2142 FORMAT(12x,' loggrid parameter r0        :',E10.3)
 2143 FORMAT(12x,' loggrid parameter rmax      :',E10.3)
 2144 FORMAT(12x,' loggrid parameter npoints   :',I10)
 2145 FORMAT(12x,' augmentation sphere radius  :',F10.3,
     >           ' (',I5,' npoints',I5,' per task)')
 2146 FORMAT(12x,' compensation sigma          :',F10.3)
 2150 FORMAT(12x,' total number of projectors  :',I10)
 2151 FORMAT(12x,' n_ps (n) l          eig      #projector')
 2152 FORMAT(14X,I3,' (',I1,') ',A,F13.6,I16)

      return
      end

