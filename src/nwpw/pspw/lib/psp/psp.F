*
* $Id: psp.F,v 1.36 2009-02-07 03:50:57 bylaska Exp $
*


*     ***********************************
*     *					*
*     *	 	  psp_init  		*
*     *					*
*     ***********************************

      subroutine psp_init()      
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal



      integer npack1,npack0
      logical value

*     **** external functions *****
      integer  ion_nkatm
      external ion_nkatm
 

      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      npsp = ion_nkatm()

c      value = MA_alloc_get(mt_dcpl,(nmax_max*lmmax_max*npack1),
c     >                    'prjtmp',prjtmp(2),prjtmp(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dbl,(nmax_max*lmmax_max*npsp*npack1),
c     >                    'vnl',vnl(2),vnl(1))
c     value = value.and.
c    >         MA_alloc_get(mt_dbl,(lmmax_max*npsp),
c    >                    'vnlnrm',vnlnrm(2),vnlnrm(1))

      value = MA_alloc_get(mt_dbl,(npsp*npack0),'vl',vl(2),vl(1))
      value = value.and.MA_alloc_get(mt_int,npsp,'vnl',vnl(2),vnl(1))
c      value = value.and.
c     >      MA_alloc_get(mt_dbl,(nmax_max*nmax_max*(lmax_max+1)*npsp),
c     >                    'Gijl',Gijl(2),Gijl(1))
      value = value.and.MA_alloc_get(mt_int,npsp,'Gijl',Gijl(2),Gijl(1))

      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'nprj',nprj(2),nprj(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'n_projector',n_projector(2),n_projector(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'l_projector',l_projector(2),l_projector(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'm_projector',m_projector(2),m_projector(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'zv',zv(2),zv(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'amass',amass(2),amass(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp*(lmax_max+1)),'rc',rc(2),rc(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'lmmax',lmmax(2),lmmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'lmax',lmax(2),lmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'locp',locp(2),locp(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'nmax',nmax(2),nmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),
     >                     'psp_type',psp_type(2),psp_type(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'rlocal',rlocal(2),rlocal(1))

      value = value.and.
     >  MA_alloc_get(mt_dbl,(npsp),'core_kin',core_kin(2),core_kin(1))

      if (.not. value) 
     >   call errquit('psp_init:out of heap memory',0, MA_ERR)

      call dcopy(npsp*npack0,          0.0d0,0,dbl_mb(vl(1)), 1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(zv(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(amass(1)),1)
      call dcopy(npsp*(lmax_max+1),    0.0d0,0,dbl_mb(rc(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(rlocal(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(core_kin(1)),1)

*     **** allocate semicore data ****
      call semicore_init()

*     *** set pawexist ****
      pawexist = .false.

      return
      end

*     ****************************************************
*     *                                                  *
*     *                psp_proj_init                     *
*     *                                                  *
*     ****************************************************
      subroutine psp_proj_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     ***** local variables ****
      logical value
      integer npack1

*     ***** external functions *****
      integer  psp_nprj_max,control_nprj_mult
      external psp_nprj_max,control_nprj_mult


      call Pack_npack(1,npack1)
      nprj_mult= control_nprj_mult()
      nprj_max = psp_nprj_max()
      value = MA_alloc_get(mt_dcpl,npack1*nprj_max*nprj_mult,
     >                     'prjtmp',prjtmp(2),prjtmp(1))
      if (.not. value)
     >   call errquit('psp_proj_init:out of heap memory',0,MA_ERR)
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
#include "mafdecls.fh"
#include "psp.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal


      logical value
      integer ia

*     **** external functions ****

*     **** deallocate semicore data ****
      call semicore_end()

*     **** deallocate nonlocal pspw data ****
      do ia=1,npsp
        if (int_mb(nprj(1)+ia-1).gt.0) then
        call psi_data_dealloc(int_mb(vnl(1)+ia-1))
        call psi_data_dealloc(int_mb(Gijl(1)+ia-1))
        end if
      end do

      value = MA_free_heap(prjtmp(2))
      value = value.and.MA_free_heap(vl(2))
      value = value.and.MA_free_heap(vnl(2))
c     value = value.and.MA_free_heap(vnlnrm(2))
      value = value.and.MA_free_heap(Gijl(2))
      value = value.and.MA_free_heap(nprj(2))
      value = value.and.MA_free_heap(n_projector(2))
      value = value.and.MA_free_heap(l_projector(2))
      value = value.and.MA_free_heap(m_projector(2))
      value = value.and.MA_free_heap(zv(2))
      value = value.and.MA_free_heap(amass(2))
      value = value.and.MA_free_heap(rc(2))
      value = value.and.MA_free_heap(lmmax(2))
      value = value.and.MA_free_heap(lmax(2))
      value = value.and.MA_free_heap(locp(2))
      value = value.and.MA_free_heap(nmax(2))
      value = value.and.MA_free_heap(psp_type(2))
      value = value.and.MA_free_heap(rlocal(2))
      value = value.and.MA_free_heap(core_kin(2))
      if (.not. value) 
     >  call errquit('psp_end:error freeing heap memory',0,MA_ERR)

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

#include "mafdecls.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal


      psp_rlocal = dbl_mb(rlocal(1)+ia-1)
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

#include "mafdecls.fh"
#include "psp.fh"


      psp_zv = dbl_mb(zv(1)+ia-1)
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
#include "psp.fh"


      psp_lmax = int_mb(lmax(1)+ia-1)
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

#include "mafdecls.fh"
#include "psp.fh"


      psp_locp = int_mb(locp(1)+ia-1)
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
*     *                                 *
*     *            psp_pawexist         *
*     *                                 *
*     ***********************************

      logical function psp_pawexist()
      implicit none

#include "psp.fh"

      psp_pawexist = pawexist
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

#include "mafdecls.fh"
#include "errquit.fh"

*     **** Error function parameters ****
      real*8  yerf,verf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,fterf,verf
c     parameter (c1=0.07052307840d0,c2=0.04228201230d0)
c     parameter (c3=0.00927052720d0)
c     parameter (c4=0.00015201430d0,c5=0.00027656720d0)
c     parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer ftmp(2)
      integer taskid_j,np_j
      integer i,j,np1,np2,np3,n2ft3d,nion
      real*8 x,y,z,q,c,r,sqrt_pi,dv,v,rx,ry,rz,fx,fy,fz

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf
      external ion_nion,ion_katm
      external lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      nion = ion_nion()

*     **** constants ****
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      dv = lattice_omega()/dble(np1*np2*np3)

*     ***** allocate temporary space ****
      if (.not.MA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1)))
     > call errquit('grad_v_lr_local:out of stack memory',0, MA_ERR)

      call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)
      do j=1,nion
         if (mod(j-1,np_j).eq.taskid_j) then
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/psp_rlocal(ion_katm(j))
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            do i=1,n2ft3d
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
            dbl_mb(ftmp(1)+3*(j-1))   = -fx*dv
            dbl_mb(ftmp(1)+3*(j-1)+1) = -fy*dv
            dbl_mb(ftmp(1)+3*(j-1)+2) = -fz*dv

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
      end do

      call Parallel_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
      call daxpy(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)

      if (.not.MA_pop_stack(ftmp(2)))
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


*     **** Error function parameters ****
      real*8 xerf,yerf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,xerf
c     parameter (c1=0.07052307840d0,c2=0.04228201230d0)
c     parameter (c3=0.00927052720d0)
c     parameter (c4=0.00015201430d0,c5=0.00027656720d0)
c     parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer taskid_j,np_j
      integer i,j,n2ft3d
      real*8 x,y,z,q,c,r,sqrt_pi

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf
      external ion_nion,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      call dcopy(n2ft3d,0.0d0,0,vlr_out,1)

      do j=1,ion_nion()
         if (mod(j-1,np_j).eq.taskid_j) then
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/psp_rlocal(ion_katm(j))

            do i=1,n2ft3d
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r .gt. 1.0d-15) then
                 xerf=r*c
c                yerf = (1.0d0
c    >                 + xerf*(c1 + xerf*(c2
c    >                 + xerf*(c3 + xerf*(c4
c    >                 + xerf*(c5 + xerf*c6))))))**4
c                yerf = (1.0d0 - 1.0d0/yerf**4)
                 yerf = util_erf(xerf)
                 vlr_out(i) = vlr_out(i) + (q/r)*yerf
c                vlr_out(i) = vlr_out(i) + (q/r)*erf(r*c)
               else
                 vlr_out(i) = vlr_out(i) + 2.0d0*q*c/sqrt_pi
               end if
            end do

         end if
      end do
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d,vlr_out)

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

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3)

*     **** external functions ****
      integer  Pack_G_indx,ion_nion,ion_katm,control_version
      real*8   ion_rion
      external Pack_G_indx,ion_nion,ion_katm,control_version
      external ion_rion

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      periodic = (control_version().eq.3)

      value = MA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) 
     >  call errquit('v_local:out of stack memory',0,MA_ERR)

*     **** define Gx,Gy and Gz in packed space ****
      if (move) then
         value = MA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))
         if (.not. value) 
     >   call errquit('v_local: out of stack memory',0, MA_ERR)
         G(1)  = Pack_G_indx(0,1)
         G(2)  = Pack_G_indx(0,2)
         G(3)  = Pack_G_indx(0,3)
         call dcopy(3*nion,0.0d0,0,fion,1)
      end if

      call dcopy((2*npack0),0.0d0,0,vl_out,1)
      do ii=1,nion
       if (mod(ii-1,np_j).eq.taskid_j) then

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
           call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(vtmp(1)))
           call Pack_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)

           if (move) then
#ifndef CRAY
!DIR$ ivdep
#endif
             do i=1,npack0
                dbl_mb(xtmp(1)+i-1) 
     >          = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >           - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
             end do
            call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),fion(1,ii))
            call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),fion(2,ii))
            call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),fion(3,ii))
           end if
          end if

       end if
      end do
      if (np_j.gt.1) call D1dB_Vector_SumAll(2*npack0,vl_out)
      if (move) call Parallel_Vector_SumAll(3*nion,fion)

      value = .true.
      if (move) value = value.and.MA_pop_stack(xtmp(2))
      value = value.and.MA_pop_stack(vtmp(2))
      value = value.and.MA_pop_stack(exi(2))

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

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      logical value,periodic,inside
      integer taskid_j,np_j
      integer npack0,nion
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
c      integer Gx(2),Gy(2),Gz(2)
      real*8 rxyz(3),fxyz(3)

*     **** external functions ****
      integer  Pack_G_indx,ion_nion,ion_katm,control_version
      real*8   ion_rion
      external Pack_G_indx,ion_nion,ion_katm,control_version
      external ion_rion

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      periodic = (control_version().eq.3)

      value = MA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

c     **** define Gx,Gy and Gz in packed space ****
      G(1)  = Pack_G_indx(0,1)
      G(2)  = Pack_G_indx(0,2)
      G(3)  = Pack_G_indx(0,3)
      call dcopy(3*nion,0.0d0,0,fion,1)

      do ii=1,nion
         if (mod(ii-1,np_j).eq.taskid_j) then

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
           call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(vtmp(1)))

#ifndef CRAY
!DIR$ ivdep
#endif
            do i=1,npack0
              dbl_mb(xtmp(1)+i-1) 
     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
           end do
           call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),fion(1,ii))
           call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),fion(2,ii))
           call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),fion(3,ii))
          end if

         end if
      end do
      call Parallel_Vector_SumAll(3*nion,fion)

      value =           MA_pop_stack(xtmp(2))
      value = value.and.MA_pop_stack(vtmp(2))
      value = value.and.MA_pop_stack(exi(2))

      call nwpw_timing_end(5)
      return 
      end




*     ***********************************
*     *					*
*     *	 	   v_nonlocal  		*
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 5-13-2002 to improve
*         parallel efficiency.
*
      subroutine v_nonlocal(ispin,ne,psi1,psi2,move,fion,fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "mafdecls.fh"
#include "psp.fh"
#include "errquit.fh"
ccccccc#include "frac_occ.fh"


*     *** local variables ***
      integer G(3),npack1,nion
      integer i,ii,ia,l,n,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),xtmp(2),sw1(2),sw2(2),sum(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   MA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   MA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('v_nonlocal: out of stack',0, MA_ERR)

      if (move) then
       value = value.and.MA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
       value = value.and.MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
       if (.not. value) 
     >  call errquit('v_nonlocal:out of stack memory',1,MA_ERR)

       G(1)  = Pack_G_indx(1,1)
       G(2)  = Pack_G_indx(1,2)
       G(3)  = Pack_G_indx(1,3)
      end if

      omega = lattice_omega()
      scal = 1.0d0/(omega)

      do ii=1,nion
        ia=ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
c        call strfac(ii,dcpl_mb(exi(1)))
c        call Pack_c_pack(1,dcpl_mb(exi(1)))
        call strfac_pack(1,ii,dcpl_mb(exi(1)))


*       **** generate sw1's and projectors ****
        do l=1,nproj

           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))

           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****


*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
cc#ifndef CRAY
c!DIR$ ivdep
cc#endif
c              do k=0,npack1-1
c               dcpl_mb(prjtmp(1)+(l-1)*npack1 + k) 
c     >           = dbl_mb(k+shift)*dcpl_mb(exi(1)+k)
c              end do

*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
cc#ifndef CRAY
c!DIR$ ivdep
cc#endif
c              do k=0,npack1-1
c               dcpl_mb(prjtmp(1)+(l-1)*npack1 + k) 
c     >                = dbl_mb(k+shift)
c     >                 *dcmplx(-dimag(dcpl_mb(exi(1)+k)),
c     >                           dble(dcpl_mb(exi(1)+k)))
c              end do

           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
        end do
        call D3dB_Vector_SumAll((nn*nproj*nprj_mult),dbl_mb(sw1(1)))


*       **** sw2 = Gijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
        call Multiply_Gijl_sw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl_indx),
     >                         dbl_mb(sw1(1)),
     >                         dbl_mb(sw2(1)))
        
*       **** do Kleinman-Bylander Multiplication ****
        !scal = 1.0d0/(omega)
        call dscal(nn*int_mb(nprj(1)+ia-1)*nprj_mult,
     >             scal,dbl_mb(sw2(1)),1)
        call DGEMM('N','T',2*npack1,nn,int_mb(nprj(1)+ia-1)*nprj_mult,
     >             (-1.0d0),
     >             dcpl_mb(prjtmp(1)), 2*npack1,
     >             dbl_mb(sw2(1)),     nn,
     >             (1.0d0),
     >             psi2,               2*npack1)


        if (move) then
        do l=1,nproj

c          call DGER((2*npack1),nn,
c     >              (-1.0d0),
c     >              dcpl_mb(prjtmp(1)+(l-1)*npack1),   1,
c     >              dbl_mb(sw2(1)+(l-1)*nn),           1,
c     >              psi2,                    (2*npack1))


             do n=1,nn
                if (ispin.eq.1) 
     >            dbl_mb(sw2(1)+n-1+(l-1)*nn)
     >            =2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change

#ifndef CRAY
!DIR$ ivdep
#endif
                do i=1,npack1
                   ctmp = psi1(i+(n-1)*npack1)
     >                  *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1 + i-1))
                   dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                end do
                call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+1+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+2+3*(n-1)))

             end do

             call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

             !**** fractional weighting ****
             if (fractional) then
              do n=1,nn
               call Dneall_qton(n,i)
               dbl_mb(sum(1)+3*(n-1))=dbl_mb(sum(1)+3*(n-1))*occ(i)
               dbl_mb(sum(1)+1+3*(n-1))=dbl_mb(sum(1)+1+3*(n-1))*occ(i)
               dbl_mb(sum(1)+2+3*(n-1))=dbl_mb(sum(1)+2+3*(n-1))*occ(i)
              end do
             end if
  
             ff(1) = 0.0d0
             ff(2) = 0.0d0
             ff(3) = 0.0d0
             do n=1,nn
                ff(1) = ff(1) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(sum(1)+  3*(n-1))
                ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(sum(1)+1+3*(n-1))
                ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(sum(1)+2+3*(n-1))
             end do
             call D1dB_Vector_SumAll(3,ff)
             fion(1,ii) = fion(1,ii)  + ff(1)
             fion(2,ii) = fion(2,ii)  + ff(2)
             fion(3,ii) = fion(3,ii)  + ff(3)
       
        end do !** l **
        end if !** move **


        end if !** nproj>0 **
      end do !** ii **

      value = .true.
      if (move) then
      value = value.and.MA_pop_stack(sum(2))
      value = value.and.MA_pop_stack(xtmp(2))
      end if
      value = value.and.MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(exi(2))
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
      integer a,b,na,nb,la,lb,ma,mb
      

      call dcopy(nn*nprj,0.0d0,0,sw2,1)
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)

         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)

            if ((la.eq.lb).and.(ma.eq.mb)) then
              na = n_prj(a)
              nb = n_prj(b)
              call daxpy(nn,G(nb,na,la),sw1(1,a),1,sw2(1,b),1)
            end if

         end do
      end do
      return
      end


*     ***********************************
*     *					*
*     *	 	   f_vnonlocal 		*
*     *					*
*     ***********************************

      subroutine f_vnonlocal(ispin,ne,psi1,fion,fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "errquit.fh"
#include "mafdecls.fh"
#include "psp.fh"
cccccccccccc#include "frac_occ.fh"

*     *** local variables ***
      integer G(3),npack1,shift,Gijl_indx
      integer i,ii,ia,k,l,n,nn,l_prj,nproj
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),vtmp(2),xtmp(2),sw1(2),sw2(2),sum(2)
c      integer Gx(2),Gy(2),Gz(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx,psi_data_get_ptr
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx,psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)
      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, npack1,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >   MA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   MA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      value = value.and.
     >      MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) 
     >  call errquit('f_vnonlocal: out of stack memory',0, MA_ERR)

c     **** define Gx,Gy and Gz in packed space ****
      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      omega = lattice_omega()

      do ii=1,ion_nion()
        ia=ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

        do l=1,nproj

           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  +(ia-1)*(nmax_max*lmmax_max))

           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))
           end if

           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_Sumall((nn*nproj),dbl_mb(sw1(1)))
 

*       **** sw2 = Gijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
        call Multiply_Gijl_sw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl_indx),
     >                         dbl_mb(sw1(1)),
     >                         dbl_mb(sw2(1)))




*       **** do Kleinman-Bylander Multiplication ****
        do l=1,nproj


           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))

           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))
           end if


          scal = 1.0d0/(omega)
          call dscal(nn,scal,dbl_mb(sw2(1)+(l-1)*nn),1)


          do n=1,nn
             if (ispin.eq.1) 
     >         dbl_mb(sw2(1)+n-1+(l-1)*nn)
     >         =2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change

#ifndef CRAY
!DIR$ ivdep
#endif
             do i=1,npack1
                ctmp = psi1(i+(n-1)*npack1)
     >               *dconjg(dcpl_mb(vtmp(1)+i-1))
                dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
             end do
             call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                         dbl_mb(sum(1)+3*(n-1)))
             call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                         dbl_mb(sum(1)+1+3*(n-1)))
             call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                         dbl_mb(sum(1)+2+3*(n-1)))

          end do

          call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

          !**** fractional weighting ****
          if (fractional) then
           do n=1,nn
            call Dneall_qton(n,i)
            dbl_mb(sum(1)+3*(n-1))  =dbl_mb(sum(1)  +3*(n-1))*occ(i)
            dbl_mb(sum(1)+1+3*(n-1))=dbl_mb(sum(1)+1+3*(n-1))*occ(i)
            dbl_mb(sum(1)+2+3*(n-1))=dbl_mb(sum(1)+2+3*(n-1))*occ(i)
           end do
          end if
  
          ff(1) = 0.0d0
          ff(2) = 0.0d0
          ff(3) = 0.0d0
          do n=1,nn
             ff(1) = ff(1) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                            *dbl_mb(sum(1)+3*(n-1))
             ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                            *dbl_mb(sum(1)+1+3*(n-1))
             ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                            *dbl_mb(sum(1)+2+3*(n-1))
          end do
          call D1dB_Vector_SumAll(3,ff)
          fion(1,ii) = fion(1,ii) + ff(1)
          fion(2,ii) = fion(2,ii) + ff(2)
          fion(3,ii) = fion(3,ii) + ff(3)
       
        end do !** l **

        end if !** nproj>0 **

      end do !** ii **

      value = MA_pop_stack(sum(2))
      value = value.and.MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(xtmp(2))
      value = value.and.MA_pop_stack(vtmp(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('f_vnonlocal: popping stack memory',0, MA_ERR)


      call nwpw_timing_end(6)
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
     >                       Gijl_tag,
     >                       rlocal,
     >                       nfft3d,npack1,npack0,
     >                       vl,vnl_tag,
     >                       semicore,rcore,ncore,
     >                       core_kin_energy,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*20 fname
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
      integer Gijl_tag

      real*8 rlocal
      integer nfft3d,npack1,npack0
      real*8 vl(*)
      integer vnl_tag
      logical semicore
      real*8  rcore
      real*8  ncore(*)

      real*8  core_kin_energy

      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#include "mafdecls.fh"
#include "util.fh"

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer n,l
      integer msglen
      integer iatom(2),Gijl_indx
      character*255 full_filename
      real*8 sum1

*     **** external functions ****
      integer  psi_data_alloc,psi_data_get_chnk
      external psi_data_alloc,psi_data_get_chnk
      real*8   lattice_omega
      external lattice_omega

      call Parallel_taskid(taskid)
      rlocal = 0.0d0

*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
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
         end if
      end if
 
*     **** send header data to all processors ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      call Parallel_Brdcst_ivalues(MASTER,msglen,psp_type)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,iatom)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,amass)
      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,zv)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp)
      call Parallel_Brdcst_ivalues(MASTER,msglen,nmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,nprj)
      lmmax=(lmax+1)**2 - (2*locp+1)

      msglen=lmax+1
      call Parallel_Brdcst_values(MASTER,msglen,rc)
      msglen = 1

      msglen=nprj
      call Parallel_Brdcst_ivalues(MASTER,msglen,n_projector)
      call Parallel_Brdcst_ivalues(MASTER,msglen,l_projector)
      call Parallel_Brdcst_ivalues(MASTER,msglen,m_projector)


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
         if (taskid.eq.MASTER) then
         call dread(5,dbl_mb(Gijl_indx),n*msglen)      !** number of matrix elements = n*nmax*nmax*(lmax+1) **
         end if
         call Parallel_Brdcst_values(MASTER,n*msglen,dbl_mb(Gijl_indx))
      end if

*     ***** read in rlocal and rcore  ****
      if (taskid.eq.MASTER) then
         if (version.eq.4) call dread(5,rlocal,1)
         call dread(5,rcore,1)
      end if
      msglen=1
      call Parallel_Brdcst_values(MASTER,msglen,rlocal)
      msglen=1
      call Parallel_Brdcst_values(MASTER,msglen,rcore)


*     **** determine semicore value ****
      if (rcore.gt.0.0d0) then
         semicore = .true.
      else
         semicore = .false.
      end if

 
*     *** read in vl 3d block ***
      call D3dB_t_read(1,5,tmp2,tmp,-1)
      call Pack_t_pack(0,tmp2)
      call Pack_t_Copy(0,tmp2,vl)


*     **** read in vnl 3d blocks ****
      if (nprj.gt.0) then
         vnl_tag = psi_data_alloc(nprj,npack1)

         do n=1,nprj
            call D3dB_t_read(1,5,tmp2,tmp,-1)
            call Pack_t_pack(1,tmp2)
            !call Pack_t_Copy(1,tmp2,vnl(1,n))
            call psi_data_add(vnl_tag,n,tmp2)
         end do
      end if

*     **** read in semicore density block ****
      if (semicore) then
         call D3dB_t_read(1,5,tmp2,tmp,-1)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1))

         call D3dB_t_read(1,5,tmp2,tmp,-1)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+2*npack0))

         call D3dB_t_read(1,5,tmp2,tmp,-1)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+3*npack0))

         call D3dB_t_read(1,5,tmp2,tmp,-1)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+4*npack0))
      end if

*     ***** Miscellaneous paw energies ****
      if (psp_type.eq.4) then
        if (taskid.eq.MASTER) call dread(5,core_kin_energy,1)
        msglen=1
        call Parallel_Brdcst_values(MASTER,msglen,core_kin_energy)
      end if

*     *** close fname binary file ***
      if (taskid.eq.MASTER) then
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
      
#include "mafdecls.fh"
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
      character*20 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind,control_print
      integer      control_ngrid
      real*8       control_unita
      character*12 control_boundry
      character*4  ion_atom
      external     nwpw_filefind,control_print
      external     control_ngrid
      external     control_unita
      external     control_boundry
      external     ion_atom
      

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call Parallel_taskid(taskid)

*     *** set semicore(0) *****
      log_mb(semicore(1)) = .false.

      value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = MA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
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
     >                  int_mb(Gijl(1)+(ia-1)), 
     >                  dbl_mb(rlocal(1) + (ia-1)),
     >                  nfft3d,npack1,npack0,
     >                  dbl_mb(vl(1) + (ia-1)*npack0),
     >                  int_mb(vnl(1)+ (ia-1)),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(rcore(1)+ia-1),
     >                  dbl_mb(ncore(1)+ (ia-1)*npack0*5),
     >                  dbl_mb(core_kin(1)+ia-1),
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
 9000 value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value)
     > call errquit('psp_readall:error popping stack',0,MA_ERR)


*     **** done reading set nprj_max and prjtmp for nonlocal psp operator ****
      call psp_proj_init()

      return
      end

