*
* $Id$
*

*     ***********************************
*     *                                 *
*     *            psp_paw_init         *
*     *                                 *
*     ***********************************
      subroutine psp_paw_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      logical value
      integer ii,ia,iii,iia,i,j
      integer nx,ny,nz,nrho,nray,npack0
      real*8  unita(3,3),bmesh,log_bmesh,pi,rcut,rs,w,sigma_smooth
      integer Gray(2),vlray(2),tmpray(2),f(2),rho(2)

*     *** external functions ***
      integer  ion_nion_qm,ion_katm_qm,ion_nkatm_qm,kbpp_calc_nray
      external ion_nion_qm,ion_katm_qm,ion_nkatm_qm,kbpp_calc_nray
      real*8   lattice_unita,control_rcut
      external lattice_unita,control_rcut

      if (pawexist) then
         nion_paw = 0
         do ii=1,ion_nion_qm()
           ia = ion_katm_qm(ii)
           if (int_mb(psp_type(1)+ia-1).eq.4) then
              nion_paw = nion_paw + 1
           end if
         end do

*        **** allocate dummy variables ****
         if (.not.BA_alloc_get(mt_int,nion_paw,"ion_pawtoion",
     >                         ion_pawtoion(2),ion_pawtoion(1)))
     >      call errquit("psp_paw_init:allocate heap",0,MA_ERR)

         iii = 0
         do ii=1,ion_nion_qm()
            ia = ion_katm_qm(ii)
            if (int_mb(psp_type(1)+ia-1).eq.4) then
               iii = iii + 1
               int_mb(ion_pawtoion(1)+iii-1) = ii
            end if
         end do
      end if

      return
      end 


    


*     ***********************************
*     *                                 *
*     *            psp_paw_end          *
*     *                                 *
*     ***********************************
      subroutine psp_paw_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

      if (pawexist) then
        if (.not.BA_free_heap(ion_pawtoion(2)))
     >      call errquit("psp_paw_end:deallocate heap",0,MA_ERR)
      end if

      return
      end




*     ***********************************
*     *				        *
*     *	    v_lr_local_seperate_paw     *
*     *					*
*     ***********************************
*
*     This routine calculates the long-range part of the
*     local pseudopotential (used by version4)
*
      subroutine v_lr_local_seperate_paw(r_grid,vlr_paw,vlr_notpaw)
      implicit none
      real*8     r_grid(3,*)
      real*8     vlr_paw(*)
      real*8     vlr_notpaw(*)

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

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      nion_qm = ion_nion_qm()

      lmbda0 = 1.0d0
      if (pspw_qmmm_lambda_flag()) lmbda0 = pspw_qmmm_lambda()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      !call dcopy(n2ft3d,0.0d0,0,vlr_paw,1)
      !call dcopy(n2ft3d,0.0d0,0,vlr_notpaw,1)
      call Parallel_shared_vector_zero(.false.,n2ft3d,vlr_paw)
      call Parallel_shared_vector_zero(.true.,n2ft3d,vlr_notpaw)

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
               vlr_notpaw(i) = vlr_notpaw(i) 
     >                       - dbl_mb(psp_fmm_Llm(1)+lm)
     >                        *dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d+i-1)
               lm = lm + 1
            end do
            end do
         end do
!$OMP END DO
      end if

      if (control_fast_erf()) then

      do j=1,ion_nion()
      
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
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)

            if ((int_mb(psp_type(1)+ia-1).eq.4)) then
!$OMP DO
               do i=1,n2ft3d_map
                  r = dsqrt( (r_grid(1,i)-x)**2
     >                     + (r_grid(2,i)-y)**2
     >                     + (r_grid(3,i)-z)**2)
                  if (r.gt.1.0d-15) then
                    xerf=r*c
                    yerf = (1.0d0
     >                    + xerf*(c1 + xerf*(c2
     >                    + xerf*(c3 + xerf*(c4
     >                    + xerf*(c5 + xerf*c6))))))**4
                    yerf = (1.0d0 - 1.0d0/yerf**4)
                    vlr_paw(i) = vlr_paw(i) + (q/r)*yerf
                  else
                    vlr_paw(i) = vlr_paw(i) + 2.0d0*q*c/sqrt_pi
                  end if
               end do
!$OMP END DO

            else
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
!$OMP DO
               do i=1,n2ft3d_map
                  r = dsqrt( (r_grid(1,i)-x)**2
     >                     + (r_grid(2,i)-y)**2
     >                     + (r_grid(3,i)-z)**2)
                  if (r.gt.1.0d-15) then
                    xerf=r*c
                    yerf = (1.0d0
     >                    + xerf*(c1 + xerf*(c2
     >                    + xerf*(c3 + xerf*(c4
     >                    + xerf*(c5 + xerf*c6))))))**4
                    yerf = (1.0d0 - 1.0d0/yerf**4)
                    vlr_notpaw(i) = vlr_notpaw(i) + lmbda*(q/r)*yerf
                  else
                    vlr_notpaw(i) = vlr_notpaw(i) 
     >                            + lmbda*2.0d0*q*c/sqrt_pi
                  end if
               end do
!$OMP END DO
            end if
            end if

         end if
      end do

      else

      do j=1,ion_nion()
 
         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if

            if ((int_mb(psp_type(1)+ia-1).eq.4)) then
!$OMP DO
               do i=1,n2ft3d_map
                  r = dsqrt( (r_grid(1,i)-x)**2
     >                     + (r_grid(2,i)-y)**2
     >                     + (r_grid(3,i)-z)**2)
                  if (r.gt.1.0d-15) then
                    xerf=r*c
                    yerf = util_erf(xerf)
                    vlr_paw(i) = vlr_paw(i) + (q/r)*yerf
                  else
                    vlr_paw(i) = vlr_paw(i) + 2.0d0*q*c/sqrt_pi
                  end if
               end do
!$OMP END DO
            else
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
!$OMP DO
               do i=1,n2ft3d_map
                  r = dsqrt( (r_grid(1,i)-x)**2
     >                     + (r_grid(2,i)-y)**2
     >                     + (r_grid(3,i)-z)**2)
                  if (r.gt.1.0d-15) then
                    xerf=r*c
                    yerf = util_erf(xerf)
                    vlr_notpaw(i) = vlr_notpaw(i) + lmbda*(q/r)*yerf
                  else
                    vlr_notpaw(i) = vlr_notpaw(i) 
     >                            + lmbda*2.0d0*q*c/sqrt_pi
                  end if
               end do
!$OMP END DO
            end if
            end if
         end if
      end do

      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_paw)
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_notpaw)

      call nwpw_timing_end(5)

      return
      end


*     **********************************
*     *                                *
*     *     f_lr_local_seperate_paw    *
*     *                                *
*     **********************************
*
*     This routine calculates the gradient of the long-range part of the
*     local pseudopotential (used by version 4)
*
*     Entry -
*     Exit -
*
      subroutine f_lr_local_seperate_paw(r_grid,rho_paw,rho_notpaw,fion)
      implicit none
      real*8  r_grid(3,*)
      real*8  rho_paw(*)
      real*8  rho_notpaw(*)
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
      logical ispawv
      integer ftmp(2)
      integer taskid_j,np_j
      integer i,j,ia,np1,np2,np3,n2ft3d_map,nion,l,m,lm,n2ft3d,nion_qm
      real*8 x,y,z,q,c,r,sqrt_pi,dv,v,rx,ry,rz,fx,fy,fz,lmbda0,lmbda

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
      nion = ion_nion()
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
     > call errquit('f_lr_local_seperate_paw:out of stack',0,MA_ERR)

      !call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)
      call Parallel_shared_vector_zero(.true.,3*nion,dbl_mb(ftmp(1)))
      if (psp_fmm) then
c         call dcopy((psp_fmm_lmax+1)**2,0.0d0,0,
c     >              dbl_mb(psp_fmm_Mlm(1)),1)
         call Parallel_shared_vector_zero(.true.,(psp_fmm_lmax+1)**2,
     >                                    dbl_mb(psp_fmm_Mlm(1)),1)
         lm = 0
         do l=0,psp_fmm_lmax
         do m=-l,l
!$OMP DO
            do i=1,n2ft3d_map
               dbl_mb(psp_fmm_Mlm(1)+lm) = dbl_mb(psp_fmm_Mlm(1)+lm) 
     >            + dbl_mb(psp_fmm_rTlm(1)+lm*n2ft3d+i-1)*rho_notpaw(i)
            end do
!$OMP END DO
            lm = lm + 1
         end do
         end do
         call D3dB_Vector_SumAll((psp_fmm_lmax+1)**2,
     >                           dbl_mb(psp_fmm_Mlm(1)))
         call FMM_fion_Mlm(psp_fmm_lmax,ion_nion(),nion_qm,
     >                     int_mb(ion_katm_ptr()),
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(psp_zv_ptr()),
     >                     lmbda0,
     >                     psp_fmm_rmax2,
     >                     dbl_mb(psp_fmm_Mlm(1)),
     >                     dbl_mb(ftmp(1)))

      end if

      if (control_fast_erf()) then

!$OMP MASTER
      do j=1,nion
         if (mod(j-1,np_j).eq.taskid_j) then
            ia=ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            ispawv = (int_mb(psp_type(1)+ia-1).eq.4)
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

               if (ispawv) then
                  fx = fx + rho_paw(i)*rx*v
                  fy = fy + rho_paw(i)*ry*v
                  fz = fz + rho_paw(i)*rz*v
               else
                  fx = fx + rho_notpaw(i)*rx*v
                  fy = fy + rho_notpaw(i)*ry*v
                  fz = fz + rho_notpaw(i)*rz*v
               end if
            end do
            dbl_mb(ftmp(1)+3*(j-1))   = -fx*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+1) = -fy*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+2) = -fz*dv*lmbda
            end if
         end if
      end do
!$OMP END MASTER
!$OMP BARRIER

      else

!$OMP MASTER
      do j=1,nion
       
         if (mod(j-1,np_j).eq.taskid_j) then
            ia= ion_katm(j)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            if (j.gt.nion_qm) then
               lmbda = lmbda0
            else
               lmbda = 1.0d0
            end if
            if ((.not.psp_fmm).or.((x*x+y*y+z*z).le.psp_fmm_rmax2)) then
            q = -psp_zv(ia)
            c = 1.0d0/psp_rlocal(ia)
            ispawv = (int_mb(psp_type(1)+ia-1).eq.4)
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

               if (ispawv) then
                  fx = fx + rho_paw(i)*rx*v
                  fy = fy + rho_paw(i)*ry*v
                  fz = fz + rho_paw(i)*rz*v
               else
                  fx = fx + rho_notpaw(i)*rx*v
                  fy = fy + rho_notpaw(i)*ry*v
                  fz = fz + rho_notpaw(i)*rz*v
               end if
            end do
            dbl_mb(ftmp(1)+3*(j-1))   = -fx*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+1) = -fy*dv*lmbda
            dbl_mb(ftmp(1)+3*(j-1)+2) = -fz*dv*lmbda

            end if
         end if
      end do
!$OMP END MASTER
!$OMP BARRIER

      end if

      call Parallel_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
      call daxpy_omp(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)

      if (.not.BA_pop_stack(ftmp(2)))
     > call errquit('f_lr_local_seperate_paw:popping stack',1,MA_ERR)

      call nwpw_timing_end(5)

      return
      end



*     ***********************************
*     *					*
*     *	     v_local_seperate_paw  	*
*     *					*
*     ***********************************

      subroutine v_local_seperate_paw(vl_paw,vl_notpaw)
      implicit none
      complex*16 vl_paw(*)
      complex*16 vl_notpaw(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm
      integer i,ii,ia
      integer exi(2)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3),lmbda0,lmbda

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
      if (.not.value) 
     > call errquit('v_local_seperate_paw:out of stack',0,MA_ERR)

      !call dcopy((2*npack0),0.0d0,0,vl_paw,1)
      !call dcopy((2*npack0),0.0d0,0,vl_notpaw,1)
      call Parallel_shared_vector_zero(.false.,2*npack0,vl_paw)
      call Parallel_shared_vector_zero(.true., 2*npack0,vl_notpaw)

      do ii=1,nion
    
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
                if (ii.gt.nion_qm) then
                   lmbda = lmbda0
                else
                   lmbda = 1.0d0
                end if

*               **** structure factor and local pseudopotential ****
                call strfac_pack(0,ii,dcpl_mb(exi(1)))
       
*               **** add to local psp ****
                if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                   call Pack_tc_MulAdd(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                                dcpl_mb(exi(1)),
     >                                vl_paw)
                else
c                   call Pack_tc_MulAdd(0,dbl_mb(vl(1)+npack0*(ia-1)),
c     >                                dcpl_mb(exi(1)),
c     >                                vl_notpaw)
                   call Pack_tc_aMulAdd(0,lmbda,
     >                                dbl_mb(vl(1)+npack0*(ia-1)),
     >                                dcpl_mb(exi(1)),
     >                                vl_notpaw)
                end if
             end if

          end if

      end do

      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(2*npack0,vl_paw)
         call D1dB_Vector_SumAll(2*npack0,vl_notpaw)
      end if

      value = BA_pop_stack(exi(2))
      if (.not.value) 
     >  call errquit('v_local_seperate_paw:popping stack',0,MA_ERR)

      call nwpw_timing_end(5)
      return 
      end

*     ***********************************
*     *                                 *
*     *      f_local_seperate_paw       *
*     *                                 *
*     ***********************************

      subroutine f_local_seperate_paw(dng_paw,dng_notpaw,fion)
      implicit none
      complex*16 dng_paw(*)
      complex*16 dng_notpaw(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion,nion_qm
      integer i,ii,ia
      integer exi(2),ftmp(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3),lmbda0,lmbda

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

      value = BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,npack0,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1))
      if (.not.value)
     > call errquit('f_local_seperate_paw:out of stack',0,MA_ERR)

      !call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)
      call Parallel_shared_vector_zero(.true.,3*nion,dbl_mb(ftmp(1)))
      G(1) = Pack_G_indx(0,1)
      G(2) = Pack_G_indx(0,2)
      G(3) = Pack_G_indx(0,3)


      do ii=1,nion

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
                if (ii.gt.nion_qm) then
                   lmbda = lmbda0
                else
                   lmbda = 1.0d0
                end if
                ia=ion_katm(ii)
*               **** structure factor and local pseudopotential ****
                call strfac_pack(0,ii,dcpl_mb(exi(1)))

*               **** add to local psp ****
                if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                   call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(vtmp(1)))
                   call Pack_cct_iconjgMulb(0,dng_paw,
     >                                 dcpl_mb(vtmp(1)),
     >                                 dbl_mb(xtmp(1)))
                else
                   call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(vtmp(1)))
                   call Pack_c_SMul1(0,lmbda,dcpl_mb(vtmp(1)))
                   call Pack_cct_iconjgMulb(0,dng_notpaw,
     >                                 dcpl_mb(vtmp(1)),
     >                                 dbl_mb(xtmp(1)))
                end if
                call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                              dbl_mb(ftmp(1)+3*(ii-1)))
                call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                              dbl_mb(ftmp(1)+3*(ii-1)+1))
                call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                              dbl_mb(ftmp(1)+3*(ii-1)+2))
             end if

          end if

      end do
      call Parallel_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
      call daxpy_omp(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)


      value =           BA_pop_stack(ftmp(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value)
     >  call errquit('f_local_seperate_paw:popping stack',0,MA_ERR)

      call nwpw_timing_end(5)
      return
      end





*     **********************************
*     *                                *
*     *       grad_v_lr_local_paw      *
*     *                                *
*     **********************************
*
*     This routine calculates the gradient of the long-range part of the
*     local pseudopotential (used by version 4)
*
*     Entry -
*     Exit -
*
      subroutine grad_v_lr_local_paw(r_grid,rho,fion)
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
      integer i,j,jj,np1,np2,np3,n2ft3d_map,nion,l,m,lm,n2ft3d
      real*8 x,y,z,q,c,r,sqrt_pi,dv,v,rx,ry,rz,fx,fy,fz

*     **** external functions ****
      logical  control_fast_erf
      integer  ion_nion,ion_katm
      real*8   lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf
      external control_fast_erf
      external ion_nion,ion_katm
      external lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d_map(1,n2ft3d_map)
      call D3dB_n2ft3d(1,n2ft3d)
      nion = ion_nion()

*     **** constants ****
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      dv = lattice_omega()/dble(np1*np2*np3)

*     ***** allocate temporary space ****
      if (.not.BA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1)))
     > call errquit('grad_v_lr_local:out of stack memory',0, MA_ERR)


      call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)

      if (control_fast_erf()) then

      do jj=1,nion_paw
         if (mod(jj-1,np_j).eq.taskid_j) then
            j = int_mb(ion_pawtoion(1)+jj-1)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/psp_rlocal(ion_katm(j))
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
            dbl_mb(ftmp(1)+3*(j-1))   = -fx*dv
            dbl_mb(ftmp(1)+3*(j-1)+1) = -fy*dv
            dbl_mb(ftmp(1)+3*(j-1)+2) = -fz*dv
            
         end if
      end do

      else

      do jj=1,nion_paw
         if (mod(jj-1,np_j).eq.taskid_j) then
            j = int_mb(ion_pawtoion(1)+jj-1)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)

            
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/psp_rlocal(ion_katm(j))
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

      end if

      call Parallel_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
      call daxpy(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)

      if (.not.BA_pop_stack(ftmp(2)))
     > call errquit('paw_grad_v_lr_local:popping stack',1,MA_ERR)

      call nwpw_timing_end(5)

      return
      end


*     ***********************************
*     *				        *
*     *	 	  v_lr_local_paw  	*
*     *					*
*     ***********************************
*
*     This routine calculates the long-range part of the
*     local pseudopotential (used by version4)
*
      subroutine v_lr_local_paw(r_grid,vlr_out)
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
      integer i,j,jj,n2ft3d,n2ft3d_map,l,m,lm
      real*8 x,y,z,q,c,r,sqrt_pi

*     **** external functions ****
      logical  control_fast_erf
      integer  ion_nion,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf
      external control_fast_erf
      external ion_nion,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf
      integer  ion_katm_ptr,ion_rion_ptr,psp_zv_ptr
      external ion_katm_ptr,ion_rion_ptr,psp_zv_ptr

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      call dcopy(n2ft3d,0.0d0,0,vlr_out,1)


      if (control_fast_erf()) then

      do jj=1,nion_paw
         if (mod(jj-1,np_j).eq.taskid_j) then
            j = int_mb(ion_pawtoion(1)+jj-1)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/psp_rlocal(ion_katm(j))

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

         end if
      end do

      else

      do jj=1,nion_paw
         if (mod(jj-1,np_j).eq.taskid_j) then
            j = int_mb(ion_pawtoion(1)+jj-1)
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
          
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/psp_rlocal(ion_katm(j))

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

         end if
      end do

      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_out)

      call nwpw_timing_end(5)

      return
      end




*     ***********************************
*     *					*
*     *	 	v_locals_paw		*
*     *					*
*     ***********************************

      subroutine v_locals_paw(vlpaw_out,move,dng,fion)
      implicit none
      complex*16 vlpaw_out(*)
      logical    move
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      integer taskid_j,np_j
      integer npack0,nion
      integer i,ii,ia,iii
      integer exi(2),vtmp(2),xtmp(2),G(3)
      logical value,periodic,inside
      real*8  rxyz(3),fxyz(3)

*     **** external functions ****
      integer  Pack_G_indx,ion_nion,ion_katm,control_version
      real*8   ion_rion
      external Pack_G_indx,ion_nion,ion_katm,control_version
      external ion_rion

      if (pawexist)  then

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      periodic = (control_version().eq.3)

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      if (.not. value) 
     >  call errquit('v_locals_paw:out of stack memory',0,MA_ERR)

*     **** define Gx,Gy and Gz in packed space ****
      if (move) then
         value = BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
         value = value.and.
     >           BA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))
         if (.not. value) 
     >   call errquit('v_locals_paw: out of stack memory',0, MA_ERR)
         G(1)  = Pack_G_indx(0,1)
         G(2)  = Pack_G_indx(0,2)
         G(3)  = Pack_G_indx(0,3)
         call dcopy(3*nion,0.0d0,0,fion,1)
      end if

      call dcopy((4*npack0),0.0d0,0,vlpaw_out,1)
      do iii=1,nion_paw
       if (mod(iii-1,np_j).eq.taskid_j) then
          ii = int_mb(ion_pawtoion(1)+iii-1)

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
           call Pack_tc_MulAdd(0,dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                      dcpl_mb(exi(1)),
     >                      vlpaw_out)
           call Pack_tc_MulAdd(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                      dcpl_mb(exi(1)),
     >                      vlpaw_out(1+npack0))

           if (move) then
              call Pack_ttcc_AddMul(0,dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                              dbl_mb(vl(1)+npack0*(ia-1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))
            call Pack_cct_iconjgMulb(0,dng,
     >                                 dcpl_mb(vtmp(1)),
     >                                 dbl_mb(xtmp(1)))
            call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),fion(1,ii))
            call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),fion(2,ii))
            call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),fion(3,ii))
           end if
          end if

       end if
      end do
      if (np_j.gt.1) then
         call D1dB_Vector_SumAll(4*npack0,vlpaw_out)
      end if
      if (move) call Parallel_Vector_SumAll(3*nion,fion)

      value = .true.
      if (move) then
         value = value.and.BA_pop_stack(xtmp(2))
         value = value.and.BA_pop_stack(vtmp(2))
      end if
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >   call errquit('v_locals_paw:popping stack',0, MA_ERR)

      call nwpw_timing_end(5)
      end if
      return 
      end




*     ***********************************
*     *					*
*     *	 	   f_vlocals_paw	*
*     *					*
*     ***********************************
      subroutine f_vlocal_paw(dng,fion)
      implicit none
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"


*     *** local variables ***
      logical value,periodic,inside
      integer taskid_j,np_j
      integer npack0,nion
      integer i,ii,ia,iii
      integer exi(2),vtmp(2),xtmp(2),G(3)
c      integer Gx(2),Gy(2),Gz(2)
      real*8 rxyz(3),fxyz(3)

*     **** external functions ****
      integer  Pack_G_indx,ion_nion,ion_katm,control_version
      real*8   ion_rion
      external Pack_G_indx,ion_nion,ion_katm,control_version
      external ion_rion

      if (pawexist)  then
      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call Pack_npack(0,npack0)
      nion     = ion_nion()
      periodic = (control_version().eq.3)

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
      call dcopy(3*nion,0.0d0,0,fion,1)

      do iii=1,nion_paw
        if (mod(iii-1,np_j).eq.taskid_j) then
          ii = int_mb(ion_pawtoion(1)+iii-1)

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
           call Pack_tc_Mul(0,dbl_mb(vlpaw(1)+npack0*(ia-1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(vtmp(1)))

c#ifndef CRAY
c!DIR$ ivdep
c#endif
c            do i=1,npack0
c              dbl_mb(xtmp(1)+i-1) 
c     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
c     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
c           end do
           call Pack_cct_iconjgMulb(0,dng,
     >                                dcpl_mb(vtmp(1)),
     >                                dbl_mb(xtmp(1)))
           call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),fion(1,ii))
           call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),fion(2,ii))
           call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),fion(3,ii))
          end if

        end if
      end do
      call Parallel_Vector_SumAll(3*nion,fion)

      value =           BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) call errquit('popping stack',0, MA_ERR)

      call nwpw_timing_end(5)
      end if
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

*     ***********************************
*     *                                 *
*     *        psp_paw_use_grid_cmp     *
*     *                                 *
*     ***********************************
      logical function psp_paw_use_grid_cmp()
      implicit none

#include "psp.fh"

      psp_paw_use_grid_cmp = use_grid_cmp
      return
      end

*     ***********************************
*     *                                 *
*     *         psp_paw_mult_l_max      *
*     *                                 *
*     ***********************************
      integer function psp_paw_mult_l_max()
      implicit none

      integer  nwpw_xc_mult_l_max,mult_l
      external nwpw_xc_mult_l_max

      mult_l = nwpw_xc_mult_l_max()
      psp_paw_mult_l_max = mult_l
      return
      end

*     ***********************************
*     *                                 *
*     *         psp_overlap_S           *
*     *                                 *
*     ***********************************

*    This routine computes the paw overlap S operator to psi1
*      psi2 = S*psi1
*
      subroutine psp_overlap_S(ispin,neq,psi1,psi2)
      implicit none
      integer    ispin,neq(2)
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer npack1,ij
      integer ii,ia,iii,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega


      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nn = neq(1)+neq(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))

      if (.not.value) 
     >  call errquit('psp_overlap_S: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal = 1.0d0/(omega)
      scalsqr = scal*scal

      !call dcopy(2*npack1*nn,psi1,1,psi2,1)
      call Parallel_shared_vector_copy(.true.,2*npack1*nn,psi1,psi2)
      do iii=1,nion_paw
        ii = int_mb(ion_pawtoion(1)+iii-1)
        ia = ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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

*           **** phase fact DOES matter for compensation charge!!!!     ****
*           **** assume that sign factor for proj is in kbpp formatting ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1, dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
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
        !nproj = int_mb(nprj(1)+ia-1)
        !call dscal(nn*int_mb(nprj(1)+ia-1),scal,dbl_mb(sw2(1)),1)
        call DSCAL_OMP(nn*nproj,scal,dbl_mb(sw2(1)),1)

        call DGEMM_OMP('N','T',2*npack1,nn,nproj,
     >             (1.0d0),
     >             dcpl_mb(prjtmp(1)), 2*npack1,
     >             dbl_mb(sw2(1)),     nn,
     >             (1.0d0),
     >             psi2,               2*npack1)
        end if !** nproj>0 **
      end do !** iii **

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap_S: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)
      return 
      end




*     ***********************************
*     *                                 *
*     *    psp_add_paw_extra_overlap1   *
*     *                                 *
*     ***********************************

*    This routine adds the extra part of the paw overlap S operator to psi1
*      S = S + Smatrix, where Smatrix = <psi1| S-I|psi1>
*
      subroutine psp_add_paw_extra_overlap1(mb,psi1,S)
      implicit none
      integer    mb
      complex*16 psi1(*)
      real*8     S(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn,ms,shifts,shiftsw
      integer k,shift,l_prj,nproj,Gijl_indx,neq(2),ishift
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6) 


cc*     **** S = transpose(psi)*psi ****
cc      call Dneall_ffm_sym_Multiply(mb,psi1,psi1,npack1,S)


      if (pawexist) then

      call Pack_npack(1,npack1)
      call Dneall_neq(neq)
      if (mb.eq.0) then
         nn = neq(1) + neq(2)
         ishift = 1
      else
         nn = neq(mb)
         ishift = 1 + (mb-1)*neq(1)*npack1
      end if

*     **** allocate local memory ****
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_paw_extra_overlap1:out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal

      do iii=1,nion_paw

        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1(ishift),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))

*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
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

*       **** S = S + sw1*transpose(sw2) ****
        call Dneall_m_add_sw1sw2(mb,nproj,scal,
     >                           dbl_mb(sw1(1)),
     >                           dbl_mb(sw2(1)),S)

c*       *** routine needs to be parallelized over orbitals ****
c        do ms=1,ispin
c          shifts  = 1+(ms-1)*ne(1)*ne(1)
c          shiftsw =   (ms-1)*ne(1)
c          call DGEMM('N','T',
c     >              ne(ms),ne(ms),int_mb(nprj(1)+ia-1),
c     >              (scal),
c     >              dbl_mb(sw1(1)+shiftsw), nn,
c     >              dbl_mb(sw2(1)+shiftsw), nn,
c     >              (1.0d0),
c     >              S(shifts), ne(ms))
c        end do
        

        end if !** nproj>0 **
      end do !** iii **

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     >   call errquit('psp_paw_extra_overlap1: popping stack',3,MA_ERR)
      end if
      call nwpw_timing_end(6)
      return 
      end


*     ***********************************
*     *                                 *
*     *    psp_add_paw_extra_overlap    *
*     *                                 *
*     ***********************************

*    This routine adds the extra part of the paw overlap S operator to psi1
*      S = S + Smatrix where Smatrix = <psi1|S-I|psi2>
*
      subroutine psp_add_paw_extra_overlap2(mb,psi1,psi2,S)
      implicit none
      integer    mb
      complex*16 psi1(*)
      complex*16 psi2(*)
      real*8     S(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn,ms,shifts,shiftsw,ishift
      integer k,shift,l_prj,nproj,Gijl_indx,neq(2)
      real*8  omega,scal,scalsqr
      integer exi(2),sw0(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6) 

c*     **** S = transpose(psi)*psi ****
c      call Dneall_ffm_Multiply(mb,psi1,psi2,npack1,S)


      if (pawexist) then

      call Pack_npack(1,npack1)
      call Dneall_neq(neq)
      if (mb.eq.0) then
         nn = neq(1) + neq(2)
         ishift = 1
      else
         nn = neq(mb)
         ishift = 1+ (mb-1)*neq(1)*npack1
      end if

*     **** allocate local memory ****
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw0',sw0(2),sw0(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     > call errquit('psp_add_paw_extra_overlap2:out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal

      do iii=1,nion_paw

        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                        psi1(ishift),
     >                        dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                        dbl_mb(sw0(1)+(l-1)*nn))
           call Pack_cc_indot(1,nn,
     >                        psi2(ishift),
     >                        dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                        dbl_mb(sw2(1)+(l-1)*nn))
        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw0(1)))
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw2(1)))

*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
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
     >                         dbl_mb(sw0(1)),
     >                         dbl_mb(sw1(1)))

*       **** S = S + sw1*transpose(sw2) ****
        call Dneall_m_add_sw1sw2(mb,nproj,scal,
     >                           dbl_mb(sw1(1)),
     >                           dbl_mb(sw2(1)),S)

c*       *** routine needs to be parallelized over orbitals ****
c        do ms=1,ispin
c          shifts  = 1+(ms-1)*ne(1)*ne(1)
c          shiftsw =   (ms-1)*ne(1)
c          call DGEMM('N','T',
c     >              ne(ms),ne(ms),int_mb(nprj(1)+ia-1),
c     >              (scal),
c     >              dbl_mb(sw1(1)+shiftsw), nn,
c     >              dbl_mb(sw2(1)+shiftsw), nn,
c     >              (1.0d0),
c     >              S(shifts), ne(ms))
c        end do
        

        end if !** nproj>0 **
      end do !** iii **

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(sw0(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('psp_add_paw_extra_overlap2:popping stack',3,MA_ERR)
      end if
      call nwpw_timing_end(6)
      return 
      end




*     ***********************************
*     *                                 *
*     *       psp_paw_overlap_fion      *
*     *                                 *
*     ***********************************

*    This routine adds the extra part of the paw overlap S operator to psi1
*      S = S + Smatrix where Smatrix = <psi1|S-I|psi2>
*
      subroutine psp_paw_overlap_fion(ispin,lmbda,psi1,fion)
      implicit none
      integer    ispin
      real*8     lmbda(*)
      complex*16 psi1(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn,ms,shifts,shiftsw
      integer k,shift,l_prj,nproj,Gijl_indx,neq(2)
      real*8  omega,scal,scalsqr,ff(3)
      integer exi(2),sw1(2),sw2(2),sw1x(2),sw1y(2),sw1z(2)
      integer S12x(2),S12y(2),S12z(2),Gx,Gy,Gz
      logical value,sd_function

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      integer  ion_katm,Pack_G_indx
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external Dneall_m_push_get,Dneall_m_pop_stack
      external ion_katm,Pack_G_indx
      external psi_data_get_ptr
      external lattice_omega


      call nwpw_timing_start(6) 
      if (pawexist) then

      call Pack_npack(1,npack1)
      call Dneall_neq(neq)
      nn = neq(1) + neq(2)


*     **** allocate local memory ****
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1x',sw1x(2),sw1x(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1y',sw1y(2),sw1y(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1z',sw1z(2),sw1z(1))
      value = value.and.Dneall_m_push_get(0,S12x)
      value = value.and.Dneall_m_push_get(0,S12y)
      value = value.and.Dneall_m_push_get(0,S12z)
      if (.not.value) 
     > call errquit('psp_paw_overlap_fion:out of stack',0,MA_ERR)

      Gx = Pack_G_indx(1,1)
      Gy = Pack_G_indx(1,2)
      Gz = Pack_G_indx(1,3)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal

      do iii=1,nion_paw

        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))
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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                        psi1,
     >                        dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                        dbl_mb(sw1(1)+(l-1)*nn))

           call Pack_conjg_tcc_indot(1,nn,
     >                        dbl_mb(Gx),
     >                        psi1,
     >                        dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                        dbl_mb(sw1x(1)+(l-1)*nn))
           call Pack_conjg_tcc_indot(1,nn,
     >                        dbl_mb(Gy),
     >                        psi1,
     >                        dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                        dbl_mb(sw1y(1)+(l-1)*nn))
           call Pack_conjg_tcc_indot(1,nn,
     >                        dbl_mb(Gz),
     >                        psi1,
     >                        dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                        dbl_mb(sw1z(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1x(1)))
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1y(1)))
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1z(1)))

*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
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

        if (ispin.eq.1) call dscal_omp(nn*nproj,2.0d0,dbl_mb(sw2(1)),1)

*       **** Sx = Sx + sw1x*transpose(sw2) ****
        call Dneall_m_zero(0,dbl_mb(S12x(1)))
        call Dneall_m_add_sw1sw2(0,nproj,scal,
     >                           dbl_mb(sw1x(1)),
     >                           dbl_mb(sw2(1)),
     >                           dbl_mb(S12x(1)))
*       **** Sy = Sy + sw1y*transpose(sw2) ****
        call Dneall_m_zero(0,dbl_mb(S12y(1)))
        call Dneall_m_add_sw1sw2(0,nproj,scal,
     >                           dbl_mb(sw1y(1)),
     >                           dbl_mb(sw2(1)),
     >                           dbl_mb(S12y(1)))
*       **** Sz = Sz + sw1z*transpose(sw2) ****
        call Dneall_m_zero(0,dbl_mb(S12z(1)))
        call Dneall_m_add_sw1sw2(0,nproj,scal,
     >                           dbl_mb(sw1z(1)),
     >                           dbl_mb(sw2(1)),
     >                           dbl_mb(S12z(1)))
!$OMP MASTER
        call Dneall_mm_sum(0,lmbda,dbl_mb(S12x(1)),ff(1))
        call Dneall_mm_sum(0,lmbda,dbl_mb(S12y(1)),ff(2))
        call Dneall_mm_sum(0,lmbda,dbl_mb(S12z(1)),ff(3))

        fion(1,ii) = fion(1,ii)  - 2.0d0*ff(1)
        fion(2,ii) = fion(2,ii)  - 2.0d0*ff(2)
        fion(3,ii) = fion(3,ii)  - 2.0d0*ff(3)
!$OMP END MASTER
!$OMP BARRIER

        end if !** nproj>0 **

      end do !** iii **

      value =           Dneall_m_pop_stack(S12z)
      value = value.and.Dneall_m_pop_stack(S12y)
      value = value.and.Dneall_m_pop_stack(S12x)
      value = value.and.BA_pop_stack(sw1z(2))
      value = value.and.BA_pop_stack(sw1y(2))
      value = value.and.BA_pop_stack(sw1x(2))
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('psp_paw_overlap_fion:popping stack',3,MA_ERR)
     
      end if
      call nwpw_timing_end(6)


      return 
      end





c*     ***********************************
c*     *					*
c*     *	       psp_overlap_orb	        *
c*     *					*
c*     ***********************************
c
c*    This routine computes the paw overlap S operator to psi1
c*      psi2 = S*psi1
c*
c      subroutine psp_overlap_orb(n,psi1,S)
c      implicit none
c      integer    n
c      complex*16 psi1(*)
c      real*8     S(*)
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c#include "psp.fh"
c
c*     *** local variables ***
c      integer npack1
c      integer ii,ia,iii,l
c      integer k,shift,l_prj,nproj,Gijl_indx
c      real*8  omega,scal,scalsqr
c      integer exi(2),sw1(2),sw2(2)
c      logical value,sd_function
c
c*     **** external functions ****
c      integer  ion_katm
c      integer  psi_data_get_ptr
c      real*8   lattice_omega
c      external ion_katm
c      external psi_data_get_ptr
c      external lattice_omega
c
c      if (pawexist) then
c      call nwpw_timing_start(6) 
c
c*     **** S = transpose(psi)*psi ****
c      call Pack_ccm_sym_dot(1,n,psi1,psi1,S)
c
c
c*     **** allocate local memory ****
c      call Pack_npack(1,npack1)
c      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
c      value = value.and.
c     >        BA_push_get(mt_dbl,n*nprj_max,'sw1',sw1(2),sw1(1))
c      value = value.and.
c     >        BA_push_get(mt_dbl,n*nprj_max,'sw2',sw2(2),sw2(1))
c      if (.not.value) 
c     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)
c
c      omega = lattice_omega()
c      scal    = 1.0d0/(omega)
c      scalsqr = scal*scal
c
c
c      do iii=1,nion_paw
c        ii    = int_mb(ion_pawtoion(1)+iii-1)
c        ia    = ion_katm(ii)
c        nproj = int_mb(nprj(1)+ia-1)
c
c        if (nproj.gt.0) then
c
c*       **** structure factor and local pseudopotential ****
c        call strfac_pack(1,ii,dcpl_mb(exi(1)))
c
c*       **** generate sw1's and projectors ****
c        do l=1,nproj
c
c           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
c           l_prj = int_mb(l_projector(1)+(l-1) 
c     >                                  + (ia-1)*(nmax_max*lmmax_max))
c           !sd_function = .not.and(l_prj,1)
c#ifdef GCC4
c           k = iand(l_prj,1)
c#else
c           k = and(l_prj,1)
c#endif
c           sd_function = (k.eq.0)
c
c
c*          **** phase factor does not matter therefore ****
c*          **** (-i)^l is the same as (i)^l in the     ****
c*          **** Rayleigh scattering formula            ****
c
c*          *** current function is s or d ****
c           if (sd_function) then
c              call Pack_tc_Mul(1,dbl_mb(shift),
c     >                           dcpl_mb(exi(1)),
c     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
c*          *** current function is p or f ****
c           else
c              call Pack_tc_iMul(1,dbl_mb(shift),
c     >                           dcpl_mb(exi(1)),
c     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
c           end if
c           call Pack_cc_indot(1,n,
c     >                      psi1,
c     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
c     >                      dbl_mb(sw1(1)+(l-1)*n))
c
c        end do
c        call D3dB_Vector_SumAll((n*nproj),dbl_mb(sw1(1)))
c
c
c*       **** sw2 = Sijl*sw1 ******
c        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
c        call Multiply_Gijl_sw1(n,
c     >                         nproj,
c     >                         int_mb(nmax(1)+ia-1),
c     >                         int_mb(lmax(1)+ia-1),
c     >                         int_mb(n_projector(1)
c     >                                + (ia-1)*(nmax_max*lmmax_max)),
c     >                         int_mb(l_projector(1)
c     >                                + (ia-1)*(nmax_max*lmmax_max)),
c     >                         int_mb(m_projector(1)
c     >                                + (ia-1)*(nmax_max*lmmax_max)),
c     >                         dbl_mb(Gijl_indx),
c     >                         dbl_mb(sw1(1)),
c     >                         dbl_mb(sw2(1)))
c
c*       *** routine needs to be parallelized over orbitals ****
c*       **** S = S + sw1*transpose(sw2) ****
c        call DGEMM('N','T',n,n,int_mb(nprj(1)+ia-1),
c     >              (scal),
c     >              dbl_mb(sw1(1)), n,
c     >              dbl_mb(sw2(1)), n,
c     >              (1.0d0),
c     >              S, n)
c
c        end if !** nproj>0 **
c      end do !** ii **
c
c      value =           BA_pop_stack(sw2(2))
c      value = value.and.BA_pop_stack(sw1(2))
c      value = value.and.BA_pop_stack(exi(2))
c      if (.not.value) call errquit('psp_overlap_orb: popping stack',3,
c     &       MA_ERR)
c      call nwpw_timing_end(6)
c      end if
c      return 
c      end
c



*     ***********************************
*     *					*
*     *	        psp_kinetic_core   	*
*     *					*
*     ***********************************
*
*    This routine returns the paw kinetic energy for the core density
*
      real*8 function psp_kinetic_core()
      implicit none

#include "bafdecls.fh"
#include "psp.fh"

*     *** local variables ***
      integer ii,ia,iii
      real*8  ecore

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      ecore = 0.0d0
      if (pawexist) then
         do iii=1,nion_paw
           ii    = int_mb(ion_pawtoion(1)+iii-1)
           ia    = ion_katm(ii)
           ecore = ecore + dbl_mb(core_kin(1)+ia-1)
         end do
      end if

      psp_kinetic_core = ecore
      return 
      end

*     ***********************************
*     *                                 *
*     *         psp_ion_core            *
*     *                                 *
*     ***********************************
*
*    This routine returns the paw ion-core energy 
*
      real*8 function psp_ion_core()
      implicit none

#include "bafdecls.fh"
#include "psp.fh"

*     *** local variables ***
      integer ii,ia,iii
      real*8  ecore

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      ecore = 0.0d0
      if (pawexist) then
         do iii=1,nion_paw
            ii    = int_mb(ion_pawtoion(1)+iii-1)
            ia    = ion_katm(ii)
            ecore = ecore + dbl_mb(core_ion(1)+ia-1)
         end do
      end if

      psp_ion_core = ecore
      return
      end


*     ***********************************
*     *					*
*     *	        psp_kinetic_atom	*
*     *					*
*     ***********************************

*    This routine computes the paw atomic kinetic energy
*
      real*8 function psp_kinetic_atom(ispin,ne,psi1)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

      real*8  kinetic_atom

*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega

      kinetic_atom = 0.0d0
      if (pawexist) then

      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do iii=1,nion_paw
        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** sw2 = Tijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),3)
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

*       **** keatom = transpose(sw1)*sw2) ****
        do l=0,(nn*nproj-1)
         kinetic_atom = kinetic_atom+dbl_mb(sw1(1)+l)*dbl_mb(sw2(1)+l)
        end do

        end if !** nproj>0 **
      end do !** ii **

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap_orb: popping stack',3,
     &       MA_ERR)

      if (ispin.eq.1) kinetic_atom = kinetic_atom+kinetic_atom
      kinetic_atom = kinetic_atom*scal

      end if
      call D1dB_SumAll(kinetic_atom)
      psp_kinetic_atom = kinetic_atom
      return 
      end


*     ***********************************
*     *					*
*     *	     psp_valence_core_atom	*
*     *					*
*     ***********************************
*
*    This routine computes the paw atomic valence_core energy
*
      real*8 function psp_valence_core_atom(ispin,ne,psi1)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

      real*8  valence_core_atom

*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega

      valence_core_atom = 0.0d0
      if (pawexist) then

      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do iii=1,nion_paw
        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** sw2 = Vcoreijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),5)
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

*       **** keatom = transpose(sw1)*sw2) ****
        do l=0,(nn*nproj-1)
           valence_core_atom = valence_core_atom
     >                       + dbl_mb(sw1(1)+l)*dbl_mb(sw2(1)+l)
        end do

        end if !** nproj>0 **
      end do !** ii **

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap_orb: popping stack',3,
     &       MA_ERR)

      if (ispin.eq.1) 
     >   valence_core_atom = valence_core_atom+valence_core_atom
      valence_core_atom = valence_core_atom*scal
      end if
      call D1dB_SumAll(valence_core_atom)

      psp_valence_core_atom = valence_core_atom
      return 
      end





*     ***********************************
*     *					*
*     *	        psp_vloc_atom		*
*     *					*
*     ***********************************

*    This routine computes the paw atomic local psp energy
*
      real*8 function psp_vloc_atom(ispin,ne,psi1)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

      real*8  vloc_atom,vloc_atom0
      common /vloc_atom_tmp/ vloc_atom0

      integer i,j,n,li,lj,mi,mj,ni,nj
      real*8 w,ee

*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega

      vloc_atom = 0.0d0
      if (pawexist) then

      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do iii=1,nion_paw
        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))

        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),4)
c        ee = 0.0d0
c        do i=1,nproj
c        do j=1,nproj
c           w = 0.0d0
c           do n=1,nn
c              w = w + dbl_mb(sw1(1)+n-1+nn*(i-1))
c     >               *dbl_mb(sw1(1)+n-1+nn*(j-1))
c           end do 
c           ni = int_mb(n_projector(1)+(i-1)
c     >                 + (ia-1)*(nmax_max*lmmax_max)) 
c           li = int_mb(l_projector(1)+(i-1)
c     >                 + (ia-1)*(nmax_max*lmmax_max)) 
c           mi = int_mb(m_projector(1)+(i-1)
c     >                 + (ia-1)*(nmax_max*lmmax_max)) 
c           nj = int_mb(n_projector(1)+(j-1)
c     >                 + (ia-1)*(nmax_max*lmmax_max)) 
c           lj = int_mb(l_projector(1)+(j-1)
c     >                 + (ia-1)*(nmax_max*lmmax_max)) 
c           mj = int_mb(l_projector(1)+(j-1)
c     >                 + (ia-1)*(nmax_max*lmmax_max)) 
c
c          if ((li.eq.lj).and.(mi.eq.mj)) then
c           ee = ee + w*scal*dbl_mb(Gijl_indx+(ni-1)
c     >             +(nj-1)*int_mb(nmax(1)+ia-1)
c     >             + li*int_mb(nmax(1)+ia-1)**2) 
c        
c          end if
c
c        end do
c        end do


*       **** sw2 = Tijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),4)
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

*       **** vloc_atom = transpose(sw1)*sw2) ****
!$OMP MASTER
        vloc_atom0 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:vloc_atom0)
        do l=0,(nn*nproj-1)
         vloc_atom0 = vloc_atom0+dbl_mb(sw1(1)+l)*dbl_mb(sw2(1)+l)
        end do
!$OMP END DO
        vloc_atom = vloc_atom + vloc_atom0

        end if !** nproj>0 **
      end do !** ii **

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_vloc_atom:pop stack',3,MA_ERR)

      if (ispin.eq.1) vloc_atom = vloc_atom+vloc_atom
      vloc_atom = vloc_atom*scal
      end if

      call D1dB_SumAll(vloc_atom)

      psp_vloc_atom = vloc_atom
      return 
      end

*     *************************************************
*     *                                               *
*     *              psp_ncmp_vloc                    *
*     *                                               *
*     *************************************************
*
*    This routine calulates the hartree energy of the compensation charge density.
*
*                  /                     /                           / 
*   E_ncmp_vloc =  | ncmp(r)*Vl1(r) dr = | ncmp_smooth(r)*Vl1(r)dr + | (ncmp-ncmp_smooth)*Vl1(r) dr
*                  /                     /                           /
*
*                  /
*                = | ncmp_smooth(r)*Vl1(r) dr                   
*                  /
*
*                  /
*                + | (ncmp(r)-ncmp_smooth(r))*(Vl1(r)-Vl2(r)) dr 
*                  /
*
*                  /
*                + | (ncmp(r)-ncmp_smooth(r))*Vl2(r) dr          
*                  /
*-----------------------------------------------------------------------------------------------------
*  
*                  /
*   E_ncmp_vloc  = | (ncmp(r)*Vl2(r) + ncmp_smooth(r)*(Vl1(r)-Vl2(r))) dr - plane-wave integrals
*                  /
*
*                  /
*                + | (ncmp(r)-ncmp_smooth(r))*(Vl1(r)-Vl2(r)) dr  - Gaussian two-center integrals
*                  /
*  
*-----------------------------------------------------------------------------------------------------

      real*8 function psp_ncmp_vloc(ispin)
      implicit none
      integer ispin

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** local variables ****
      logical ok,periodic,move
      integer npack0,dng_cmp(2),dng_cmp_smooth(2),vl1(2),vl_notpaw(2)
      integer n2ft3d,rho_cmp(2),rho_cmp_smooth(2),vl2(2),r_grid(2)
      integer nx,ny,nz
      real*8 eh,e1,e2,eh0,eh1,eh2,scal1,dv
      real*8 dum1(3),dum2(3)

*     **** external functions ****
      integer  control_version
      external control_version
      real*8   nwpw_compcharge_E_multipole_zv_ee,lattice_omega
      external nwpw_compcharge_E_multipole_zv_ee,lattice_omega


      eh = 0.0d0
      if (pawexist) then

      periodic = (control_version().eq.3)

*     *************************************
*     **** Periodic Boundary Condtions ****
*     *************************************
      if (periodic) then
         call Pack_npack(0,npack0)
         ok =        BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                           dng_cmp(2),dng_cmp(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                           dng_cmp_smooth(2),dng_cmp_smooth(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vl1',
     >                           vl1(2),vl1(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vl_notpaw',
     >                           vl_notpaw(2),vl_notpaw(1))
         if (.not.ok)
     >      call errquit('psp_ncmp_vloc:out of stack',0,MA_ERR)

*        **** Using pw grid only ***
         if (use_grid_cmp) then
            call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
            move = .false.
            call v_local(dcpl_mb(vl1(1)),
     >               move,
     >               dum1,dum2)
            call Pack_cc_dot(0,dcpl_mb(dng_cmp(1)),
     >                         dcpl_mb(vl1(1)),eh)

*        **** Using gaussian two-electron integrals and pw grid ***
         else
            eh = nwpw_compcharge_E_multipole_zv_ee(ispin,dbl_mb(zv(1)))  !*** Gaussian integrals

            call v_local_seperate_paw(dcpl_mb(vl1(1)),
     >                                dcpl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_v_cmp_smooth(dbl_mb(zv(1)),
     >                                            dcpl_mb(dng_cmp(1)))
            call Pack_cc_Sub2(0,dcpl_mb(dng_cmp(1)),
     >                          dcpl_mb(vl1(1)))
            call Pack_cc_Sum2(0,dcpl_mb(dng_cmp(1)),
     >                          dcpl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dcpl_mb(dng_cmp(1)),
     >                                       dcpl_mb(dng_cmp_smooth(1)))
            call Pack_cc_dot(0,dcpl_mb(dng_cmp(1)),
     >                         dcpl_mb(vl_notpaw(1)),e1)

            call Pack_cc_dot(0,dcpl_mb(dng_cmp_smooth(1)),
     >                         dcpl_mb(vl1(1)),e2)
            eh = eh + (e1+e2)
         end if
         ok =        BA_pop_stack(vl_notpaw(2))
         ok = ok.and.BA_pop_stack(vl1(2))
         ok = ok.and.BA_pop_stack(dng_cmp_smooth(2))
         ok = ok.and.BA_pop_stack(dng_cmp(2))
         if (.not.ok)
     >      call errquit('psp_vloc_residual:pop stack',1,MA_ERR)

*     **************************************
*     **** APeriodic Boundary Condtions ****
*     **************************************
      else
         call D3dB_n2ft3d(1,n2ft3d)
         ok =         BA_push_get(mt_dbl,n2ft3d,'rho_cmp',
     >                            rho_cmp(2),rho_cmp(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'rho_cmp_smooth',
     >                           rho_cmp_smooth(2),rho_cmp_smooth(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl1',
     >                           vl1(2),vl1(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl_notpaw',
     >                           vl_notpaw(2),vl_notpaw(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl2',
     >                           vl2(2),vl2(1))
         ok = ok.and.BA_push_get(mt_dbl,3*n2ft3d,'rgrid_cmp',
     >                           r_grid(2),r_grid(1))
         if (.not.ok)
     >      call errquit('psp_ncmp_vloc:out of stack',2,MA_ERR)
         call lattice_r_grid(dbl_mb(r_grid(1)))
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/(nx*ny*nz)
         dv = scal1*lattice_omega()

*        **** Using pw grid only ***
         if (use_grid_cmp) then
            call nwpw_compcharge_gen_dn_cmp(ispin,dbl_mb(rho_cmp(1)))
     >                                     
            move = .false.
            call v_local(dbl_mb(vl1(1)),
     >               move,
     >               dum1,dum2)
            call Pack_cc_dot(0,dbl_mb(rho_cmp(1)),
     >                         dbl_mb(vl1(1)),eh)

            call Pack_c_unpack(0,dbl_mb(rho_cmp(1)))
            call D3dB_cr_fft3b(1,dbl_mb(rho_cmp(1)))

            call v_lr_local(dbl_mb(r_grid(1)),
     >                      dbl_mb(vl1(1)))
            call D3dB_rr_dot(1,dbl_mb(rho_cmp(1)),
     >                         dbl_mb(vl1(1)),e1)
            eh = eh + e1*dv

*        **** Using gaussian two-electron integrals and pw grid ***
         else
            eh = nwpw_compcharge_E_multipole_zv_ee(ispin,dbl_mb(zv(1)))

*           **** long-range terms ****
            call v_lr_local_seperate_paw(dbl_mb(r_grid(1)),
     >                                   dbl_mb(vl1(1)),
     >                                   dbl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_vlr_cmp_smooth(dbl_mb(zv(1)),
     >                                           dbl_mb(r_grid(1)),
     >                                           dbl_mb(rho_cmp(1)))

            call D3dB_rr_Sum2(1,dbl_mb(rho_cmp(1)),dbl_mb(vl_notpaw(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(vl_notpaw(1)))
            call Pack_c_pack(0,dbl_mb(vl_notpaw(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(vl_notpaw(1)))

            call D3dB_rr_Sub2(1,dbl_mb(rho_cmp(1)),dbl_mb(vl1(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(vl1(1)))
            call Pack_c_pack(0,dbl_mb(vl1(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(vl1(1)))

*           **** short-range terms ****
            call v_local_seperate_paw(dbl_mb(rho_cmp_smooth(1)),
     >                                dbl_mb(rho_cmp(1)))
            call Pack_cc_Sum2(0,
     >                        dbl_mb(rho_cmp(1)),
     >                        dbl_mb(vl_notpaw(1)))
            call Pack_cc_Sum2(0,
     >                        dbl_mb(rho_cmp_smooth(1)),
     >                        dbl_mb(vl1(1)))

            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dbl_mb(rho_cmp(1)),
     >                                       dbl_mb(rho_cmp_smooth(1)))
            call Pack_cc_dot(0,dbl_mb(rho_cmp(1)),
     >                         dbl_mb(vl_notpaw(1)),e1)
            call Pack_cc_dot(0,dbl_mb(rho_cmp_smooth(1)),
     >                         dbl_mb(vl1(1)),e2)
            eh = eh + (e1+e2)

         end if
         ok =        BA_pop_stack(r_grid(2))
         ok = ok.and.BA_pop_stack(vl2(2))
         ok = ok.and.BA_pop_stack(vl_notpaw(2))
         ok = ok.and.BA_pop_stack(vl1(2))
         ok = ok.and.BA_pop_stack(rho_cmp_smooth(2))
         ok = ok.and.BA_pop_stack(rho_cmp(2))
         if (.not.ok)
     >      call errquit('psp_vloc_ncmp_vloc:pop stack',3,MA_ERR)

      end if


      end if

      psp_ncmp_vloc = eh
      return
      end

*     *************************************************
*     *                                               *
*     *              psp_dE_ncmp_vloc_Qlm             *
*     *                                               *
*     *************************************************

      subroutine psp_dE_ncmp_vloc_Qlm(ispin,move,fion)
      implicit none
      integer ispin
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** local variables ****
      logical ok,periodic
      integer npack0,dng_cmp(2),dng_cmp_smooth(2),vl1(2),vl_notpaw(2)
      integer n2ft3d,rho_cmp(2),rho_cmp_smooth(2),vl2(2),r_grid(2)
      integer nx,ny,nz
      real*8 eh,e1,e2,eh0,eh1,eh2,scal1,dv
      real*8 dum1(3),dum2(3)

*     **** external functions ****
      integer  control_version
      external control_version
      real*8   nwpw_compcharge_E_multipole_zv_ee,lattice_omega
      external nwpw_compcharge_E_multipole_zv_ee,lattice_omega


      eh = 0.0d0
      if (pawexist) then

      periodic = (control_version().eq.3)

*     *************************************
*     **** Periodic Boundary Condtions ****
*     *************************************
      if (periodic) then
         call Pack_npack(0,npack0)
         ok =        BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                           dng_cmp(2),dng_cmp(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                           dng_cmp_smooth(2),dng_cmp_smooth(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vl1',
     >                           vl1(2),vl1(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vl_notpaw',
     >                           vl_notpaw(2),vl_notpaw(1))
         if (.not.ok)
     >      call errquit('psp_ncmp_vloc:out of stack',0,MA_ERR)

*        **** Using pw grid only ***
         if (use_grid_cmp) then
c            call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
            call v_local(dcpl_mb(vl1(1)),
     >               .false.,
     >               dum1,dum2)
c            call Pack_cc_dot(0,dcpl_mb(dng_cmp(1)),
c     >                         dcpl_mb(vl1(1)),eh)

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_pw(ispin,
     >                                               dcpl_mb(vl1(1)),
     >                                               move,fion)

*        **** Using gaussian two-electron integrals and pw grid ***
         else
            call v_local_seperate_paw(dcpl_mb(vl1(1)),
     >                                dcpl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_v_cmp_smooth(dbl_mb(zv(1)),
     >                                            dcpl_mb(dng_cmp(1)))

            call Pack_cc_Sub2(0,dcpl_mb(dng_cmp(1)),
     >                          dcpl_mb(vl1(1)))
            call Pack_cc_Sum2(0,dcpl_mb(dng_cmp(1)),
     >                          dcpl_mb(vl_notpaw(1)))

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm(ispin,
     >                                            dbl_mb(zv(1)),
     >                                            dcpl_mb(vl_notpaw(1)),
     >                                            dcpl_mb(vl1(1)),
     >                                            move,fion)
            if (move) then
               call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dcpl_mb(dng_cmp(1)),
     >                                       dcpl_mb(dng_cmp_smooth(1)))
               call f_local_seperate_paw(dcpl_mb(dng_cmp_smooth(1)),
     >                                dcpl_mb(dng_cmp(1)),
     >                                fion)
               call Pack_cc_Sub2(0,dcpl_mb(dng_cmp_smooth(1)),
     >                             dcpl_mb(dng_cmp(1)))
               call nwpw_compcharge_gen_f_cmp_smooth(dbl_mb(zv(1)),
     >                                            dcpl_mb(dng_cmp(1)),
     >                                            fion)
            end if

         end if
         ok =        BA_pop_stack(vl_notpaw(2))
         ok = ok.and.BA_pop_stack(vl1(2))
         ok = ok.and.BA_pop_stack(dng_cmp_smooth(2))
         ok = ok.and.BA_pop_stack(dng_cmp(2))
         if (.not.ok)
     >      call errquit('psp_vloc_residual:pop stack',1,MA_ERR)

*     **************************************
*     **** APeriodic Boundary Condtions ****
*     **************************************
      else
         call D3dB_n2ft3d(1,n2ft3d)
         ok =         BA_push_get(mt_dbl,n2ft3d,'rho_cmp',
     >                            rho_cmp(2),rho_cmp(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'rho_cmp_smooth',
     >                           rho_cmp_smooth(2),rho_cmp_smooth(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl1',
     >                           vl1(2),vl1(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl_notpaw',
     >                           vl_notpaw(2),vl_notpaw(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl2',
     >                           vl2(2),vl2(1))
         ok = ok.and.BA_push_get(mt_dbl,3*n2ft3d,'rgrid_cmp',
     >                           r_grid(2),r_grid(1))
         if (.not.ok)
     >      call errquit('psp_ncmp_vloc:out of stack',2,MA_ERR)
         call lattice_r_grid(dbl_mb(r_grid(1)))
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/(nx*ny*nz)
         dv = scal1*lattice_omega()

*        **** Using pw grid only ***
         if (use_grid_cmp) then
            call v_local(dbl_mb(vl1(1)),
     >               .false.,
     >               dum1,dum2)
            call v_lr_local(dbl_mb(r_grid(1)),
     >                      dbl_mb(rho_cmp(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(rho_cmp(1)))
            call Pack_c_pack(0,dbl_mb(rho_cmp(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(rho_cmp(1)))
            call Pack_cc_Sum2(0,dbl_mb(rho_cmp(1)),dbl_mb(vl1(1)))

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_pw(ispin,
     >                                               dbl_mb(vl1(1)),
     >                                               move,fion)

*        **** Using gaussian two-electron integrals and pw grid ***
         else

*           **** long-range terms ****
            call v_lr_local_seperate_paw(dbl_mb(r_grid(1)),
     >                                   dbl_mb(vl1(1)),
     >                                   dbl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_vlr_cmp_smooth(dbl_mb(zv(1)),
     >                                           dbl_mb(r_grid(1)),
     >                                           dbl_mb(rho_cmp(1)))
            call D3dB_rr_Sum2(1,dbl_mb(rho_cmp(1)),dbl_mb(vl_notpaw(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(vl_notpaw(1)))
            call Pack_c_pack(0,dbl_mb(vl_notpaw(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(vl_notpaw(1)))

            call D3dB_rr_Sub2(1,dbl_mb(rho_cmp(1)),dbl_mb(vl1(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(vl1(1)))
            call Pack_c_pack(0,dbl_mb(vl1(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(vl1(1)))

*           **** short-range terms ****
            call v_local_seperate_paw(dbl_mb(rho_cmp_smooth(1)),
     >                                dbl_mb(rho_cmp(1)))
            call Pack_cc_Sum2(0,
     >                        dbl_mb(rho_cmp(1)),
     >                        dbl_mb(vl_notpaw(1)))
            call Pack_cc_Sum2(0,
     >                        dbl_mb(rho_cmp_smooth(1)),
     >                        dbl_mb(vl1(1)))

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm(ispin,
     >                                            dbl_mb(zv(1)),
     >                                            dbl_mb(vl_notpaw(1)),
     >                                            dbl_mb(vl1(1)),
     >                                            move,fion)
            if (move) then
               call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dbl_mb(rho_cmp(1)),
     >                                       dbl_mb(rho_cmp_smooth(1)))
               call f_local_seperate_paw(dbl_mb(rho_cmp_smooth(1)),
     >                                   dbl_mb(rho_cmp(1)),
     >                                   fion)

               call Pack_c_unpack(0,dbl_mb(rho_cmp(1)))
               call D3dB_cr_pfft3b(1,0,dbl_mb(rho_cmp(1)))

               call Pack_c_unpack(0,dbl_mb(rho_cmp_smooth(1)))
               call D3dB_cr_pfft3b(1,0,dbl_mb(rho_cmp_smooth(1)))

               call f_lr_local_seperate_paw(dbl_mb(r_grid(1)),
     >                                      dbl_mb(rho_cmp_smooth(1)),
     >                                      dbl_mb(rho_cmp(1)),
     >                                      fion)
               call D3dB_rr_Sub2(1,dbl_mb(rho_cmp_smooth(1)),
     >                             dbl_mb(rho_cmp(1)))
               call nwpw_compcharge_gen_f_lr_cmp_smooth(dbl_mb(zv(1)),
     >                                               dbl_mb(r_grid(1)),
     >                                               dbl_mb(rho_cmp(1)),
     >                                               fion)
            end if
         end if
         ok =        BA_pop_stack(r_grid(2))
         ok = ok.and.BA_pop_stack(vl2(2))
         ok = ok.and.BA_pop_stack(vl_notpaw(2))
         ok = ok.and.BA_pop_stack(vl1(2))
         ok = ok.and.BA_pop_stack(rho_cmp_smooth(2))
         ok = ok.and.BA_pop_stack(rho_cmp(2))
         if (.not.ok)
     >      call errquit('psp_vloc_ncmp_vloc:pop stack',3,MA_ERR)

      end if


      end if

      return
      end


*     ***********************************
*     *					*
*     *	        psp_xc_atom		*
*     *					*
*     ***********************************

*    This routine computes the paw atomic exc and pxc psp energies
*
      subroutine psp_xc_atom(ispin,ne,psi1,exc,pxc)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 exc,pxc

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

      real*8 pxc0
      common /pxc_atom_pdx0/ pxc0


*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_chnk,psi_data_get_ptr
      real*8   lattice_omega,nwpw_xc_energy_atom
      external ion_katm
      external psi_data_get_chnk,psi_data_get_ptr
      external lattice_omega,nwpw_xc_energy_atom


      exc = 0.0d0
      pxc = 0.0d0
!$OMP MASTER
      pxc0 = 0.0d0
!$OMP END MASTER

      if (pawexist) then
      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)


      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do iii=1,nion_paw
        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** sw2 = sw2 + Vxcijl*sw1 ******
        !call dcopy(nn*nproj,0.0d0,0,dbl_mb(sw2(1)),1)
        call Parallel_shared_vector_zero(.true.,nn*nproj,dbl_mb(sw2(1)))
        call nwpw_xc_solve(ii,ia,
     >     int_mb(n1dgrid(1)+ia-1),
     >     int_mb(n1dbasis(1)+ia-1),
     >     dbl_mb(psi_data_get_chnk(int_mb(phi_ae(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(phi_ps(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(dphi_ae(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(dphi_ps(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(core_ae(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(core_ps(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(core_ae_prime(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(core_ps_prime(1)+ia-1))),
     >     dbl_mb(psi_data_get_chnk(int_mb(rgrid(1)+ia-1))),
     >     dbl_mb(log_amesh(1)+ia-1),
     >     ispin,ne,int_mb(nprj(1)+ia-1),dbl_mb(sw1(1)),dbl_mb(sw2(1)))

*          **** pxc = transpose(sw1)*sw2) ****
!$OMP MASTER
           do l=0,(nn*nproj-1)
            pxc0 = pxc0 + dbl_mb(sw1(1)+l)*dbl_mb(sw2(1)+l)
           end do
!$OMP END MASTER

        end if !** nproj>0 **
      end do !** ii **
!$OMP MASTER
      pxc0 = pxc0*scal 
      if (ispin.eq.1) pxc0 = pxc0 + pxc0
!$OMP END MASTER
!$OMP BARRIER

      value =           BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_xc_atom: popping stack',3,
     &       MA_ERR)

      call D1dB_SumAll(pxc0)
      pxc = pxc0
      exc = nwpw_xc_energy_atom()
      end if
      return 
      end




*     *******************************************************
*     *                                                     *
*     *                 psp_rholm_solve                     *
*     *                                                     *
*     *******************************************************

      subroutine psp_rholm_solve(ispin,ne,nproj,sw1,
     >                            l_prj,m_prj,projtobasis,
     >                            n1dgrid,n1dbasis,
     >                            rgrid,phi_ae,phi_ps,
     >                            lmax,lmax2,
     >                            rholm_ae,rholm_ps)
      implicit none
      integer ispin,ne(2),nproj
      real*8  sw1(ne(1)+ne(2),nproj)
      integer l_prj(*), m_prj(*),projtobasis(*)

      integer n1dgrid,n1dbasis
      real*8  rgrid(n1dgrid)
      real*8  phi_ae(n1dgrid,n1dbasis)
      real*8  phi_ps(n1dgrid,n1dbasis)
      integer lmax,lmax2
      real*8 rholm_ae(n1dgrid,ispin,lmax2)
      real*8 rholm_ps(n1dgrid,ispin,lmax2)

*     ***** local variables *****
      integer i,j,l,m,lm,ms,n,ig,n1(2),n2(2)
      real*8  wij,taunt

*     ***** external functions *****
      real*8   taunt_coeff
      external taunt_coeff

      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)

      do i=1,nproj
         do j=1,nproj

*           **** generate overlap matrix wij(ms) = Sum(n=1,ne(ms)) <psi(n)|prj(i)> * <prj(j)*psi(n)> ****
            do ms=1,ispin
               wij = 0.0
               do n=n1(ms),n2(ms)
                  wij = wij + sw1(n,i)*sw1(n,j)
               end do

               do ig=1,n1dgrid
                  rholm_ae(ig,ms,1) = wij
     >                               *phi_ae(ig,projtobasis(i))
     >                               *phi_ae(ig,projtobasis(j))
     >                               /rgrid(ig)**2
                  rholm_ps(ig,ms,1) = wij
     >                               *phi_ps(ig,projtobasis(i))
     >                               *phi_ps(ig,projtobasis(j))
     >                               /rgrid(ig)**2
               end do
            end do

            lm = 2
            do l=1,lmax
               do m=-l,l
c                  taunt = taunt_coeff(l,m,
c     >                                l_prj(j),m_prj(j),
c     >                                l_prj(i),m_prj(i))
                  do ms=1,ispin
                     do ig=1,n1dgrid
                        rholm_ae(ig,ms,lm) = taunt*rholm_ae(ig,ms,1)
                        rholm_ps(ig,ms,lm) = taunt*rholm_ps(ig,ms,1)
                     end do
                  end do
                 lm = lm + 1
               end do
            end do
c            taunt = taunt_coeff(0,0,
c     >                          l_prj(j),m_prj(j),
c     >                          l_prj(i),m_prj(i))
            do ms=1,ispin
               do ig=1,n1dgrid
                  rholm_ae(ig,ms,1) = taunt*rholm_ae(ig,ms,1)
                  rholm_ps(ig,ms,1) = taunt*rholm_ps(ig,ms,1)
               end do
            end do

         end do
      end do
      return
      end


*     ***********************************
*     *					*
*     *	        psp_qlm_atom		*
*     *					*
*     ***********************************

*    This routine computes the multipole expansion
*
      subroutine psp_qlm_atom(ispin,ne,psi1)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1
      integer ii,ia,iii,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function


*     **** external functions ****
      integer  ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_katm
      external psi_data_get_ptr
      external lattice_omega

      if (pawexist) then

      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      if (.not.value) 
     >  call errquit('psp_qlm_atom: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal

      do iii=1,nion_paw
        ii    = int_mb(ion_pawtoion(1)+iii-1)
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** paw atom - generate it's atomic density matrix ****
        call psp_gen_density_matrix(ispin,ne,nproj,
     >                              dbl_mb(sw1(1)),
     >                              dbl_mb(wtmp(1)))

*       **** paw atom - generate it's compcharge ***
        call nwpw_compcharge_gen_Qlm(ii,ia,ispin,nproj,
     >                               dbl_mb(wtmp(1)))



        end if !** nproj>0 **
      end do !** ii **

      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_qlm_atom: popping stack',3,
     >       MA_ERR)

      end if
      return 
      end


*     *************************************************
*     *                                               *
*     *              psp_hartree_cmp_cmp              *
*     *                                               *
*     *************************************************
*
*    This routine calulates the hartree energy of the compensation charge density.
*
*                         //
*    E_cmp_cmp   = 0.5 * || (ncmp(r))*(ncmp(r'))
*                        || --------------------  dr dr'
*                       //       |r-r'|
*
*
*
*    using either a Fourier grid (use_grid_cmp=.true.)
*    or a combinations of Fourier grids and two electron two center Gaussian Coulomb integrals (use_grid_cmp=.false.).
*
      real*8 function psp_hartree_cmp_cmp(ispin)
      implicit none
      integer ispin

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** local variables ****
      logical ok,periodic
      integer npack0,dng_cmp(2),dng_cmp_smooth(2),vcmp_smooth(2)
      integer n2ft3d,rho_cmp(2),rho_cmp_smooth(2)
      integer nx,ny,nz
      real*8 eh,e1,e2,eh0,eh1,eh2,scal1,dv
      common /psp_cmp_eh12/ eh,e1,e2
      
*     **** external functions ****
      integer  control_version
      external control_version
      real*8   nwpw_compcharge_E_multipole_zv_ee,lattice_omega
      external nwpw_compcharge_E_multipole_zv_ee,lattice_omega
      real*8   nwpw_compcharge_E_multipole_zv
      external nwpw_compcharge_E_multipole_zv
      real*8   nwpw_compcharge_E_multipole_zv_zv
      external nwpw_compcharge_E_multipole_zv_zv
      real*8   nwpw_compcharge_E_multipole
      external nwpw_compcharge_E_multipole

      if (pawexist) then

      periodic = (control_version().eq.3)

*     *************************************
*     **** Periodic Boundary Condtions ****
*     *************************************
      if (periodic) then
         call Pack_npack(0,npack0)
         ok =        BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                           dng_cmp(2),dng_cmp(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vcmp_smooth',
     >                           vcmp_smooth(2),vcmp_smooth(1))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_cmp:out of stack',0,MA_ERR)

*        **** Using pw grid ***
         if (use_grid_cmp) then
            call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
            call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp_smooth(1)))
            call Pack_cc_dot(0,dcpl_mb(dng_cmp(1)),
     >                         dcpl_mb(vcmp_smooth(1)),e1)
            eh = 0.5d0*e1*lattice_omega()

*        **** Using gaussian Multipole energies ****
         else

*           **** multipole energy ****
            eh0 = nwpw_compcharge_E_multipole(ispin)

            ok = BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                       dng_cmp_smooth(2),dng_cmp_smooth(1))
            if (.not.ok)
     >         call errquit('psp_hartree_cmp_cmp:out of stack',1,MA_ERR)

            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dcpl_mb(dng_cmp(1)),
     >                                       dcpl_mb(dng_cmp_smooth(1)))
            call coulomb_v(dcpl_mb(dng_cmp_smooth(1)),
     >                     dcpl_mb(vcmp_smooth(1)))
            call Pack_cc_dot(0,dcpl_mb(dng_cmp(1)),
     >                         dcpl_mb(vcmp_smooth(1)),e1)
            call Pack_cc_dot(0,dcpl_mb(dng_cmp_smooth(1)),
     >                         dcpl_mb(vcmp_smooth(1)),e2)

!$OMP MASTER
*           **** add cmp energies ****
            eh = eh0 + (e1-0.5d0*e2)*lattice_omega()
!$OMP END MASTER


            ok = ok.and.BA_pop_stack(dng_cmp_smooth(2))
            if (.not.ok)
     >        call errquit('psp_hartree_cmp_cmp:popping stack',2,MA_ERR)

         end if

         ok =        BA_pop_stack(vcmp_smooth(2))
         ok = ok.and.BA_pop_stack(dng_cmp(2))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_cmp:popping stack',3,MA_ERR)


*     **************************************
*     **** APeriodic Boundary Condtions ****
*     **************************************
      else
         call D3dB_n2ft3d(1,n2ft3d)
         ok =         BA_push_get(mt_dbl,n2ft3d,'rho_cmp',
     >                            rho_cmp(2),rho_cmp(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vcmp_smooth',
     >                           vcmp_smooth(2),vcmp_smooth(1))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_cmp:out of stack',4,MA_ERR)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/(nx*ny*nz)
         dv = scal1*lattice_omega()

*        **** Using pw grid ***
         if (use_grid_cmp) then
            call nwpw_compcharge_gen_dn_cmp(ispin,dbl_mb(rho_cmp(1)))

            call Pack_c_unpack(0,dbl_mb(rho_cmp(1)))
            call D3dB_cr_fft3b(1,dbl_mb(rho_cmp(1)))
            call coulomb2_v(dbl_mb(rho_cmp(1)),
     >                      dbl_mb(vcmp_smooth(1)))
            call D3dB_rr_dot(1,dbl_mb(rho_cmp(1)),
     >                         dbl_mb(vcmp_smooth(1)),e1)
!$OMP MASTER
            eh = 0.5d0*e1*dv
!$OMP END MASTER

*        **** Using gaussian Multipole energies ****
         else

*           **** multipole energy ****
            eh0  = nwpw_compcharge_E_multipole(ispin)
c            eh0 = nwpw_compcharge_E_multipole_zv(ispin,dbl_mb(zv(1)))
c            eh2 = nwpw_compcharge_E_multipole_zv_zv(ispin,dbl_mb(zv(1)))
c            eh1 = nwpw_compcharge_E_multipole_zv_ee(ispin,dbl_mb(zv(1)))


            ok = BA_push_get(mt_dbl,n2ft3d,'rho_cmp_smooth',
     >                       rho_cmp_smooth(2),rho_cmp_smooth(1))
            if (.not.ok)
     >         call errquit('psp_hartree_cmp_cmp:out of stack',5,MA_ERR)

            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dbl_mb(rho_cmp(1)),
     >                                       dbl_mb(rho_cmp_smooth(1)))

            call Pack_c_unpack(0,dbl_mb(rho_cmp(1)))
            call D3dB_cr_fft3b(1,dbl_mb(rho_cmp(1)))
            call Pack_c_unpack(0,dbl_mb(rho_cmp_smooth(1)))
            call D3dB_cr_fft3b(1,dbl_mb(rho_cmp_smooth(1)))

            call coulomb2_v(dbl_mb(rho_cmp_smooth(1)),
     >                      dbl_mb(vcmp_smooth(1)))
            call D3dB_rr_dot(1,dbl_mb(rho_cmp(1)),
     >                         dbl_mb(vcmp_smooth(1)),e1)
            call D3dB_rr_dot(1,dbl_mb(rho_cmp_smooth(1)),
     >                         dbl_mb(vcmp_smooth(1)),e2)

!$OMP MASTER
            eh = eh + (e1-0.5d0*e2)*dv
!$OMP END MASTER

            ok = BA_pop_stack(rho_cmp_smooth(2))
            if (.not.ok)
     >        call errquit('psp_hartree_cmp_cmp:popping stack',2,MA_ERR)
         end if

         ok =        BA_pop_stack(vcmp_smooth(2))
         ok = ok.and.BA_pop_stack(rho_cmp(2))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_cmp:popping stack',6,MA_ERR)
      end if

      end if
!$OMP BARRIER
      psp_hartree_cmp_cmp = eh
      return
      end

*     *************************************************
*     *                                               *
*     *              psp_hartree_cmp_pw               *
*     *                                               *
*     *************************************************
*
*    This routine calulates the hartree energy of the compensation charge density.
*
*                         //
*    E_cmp_smooth     =  || (ncmp(r))*(npw(r'))
*                        || -----------------------  dr dr'
*                       //         |r-r'|
*
*    using either a Fourier grid 
*
      real*8 function psp_hartree_cmp_pw(ispin,dng,rho)
      implicit none
      integer ispin
      complex*16 dng(*)
      real*8     rho(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** local variables ****
      logical ok,periodic
      integer npack0,dng_cmp(2),vcmp(2)
      integer n2ft3d,rho_cmp(2),nx,ny,nz
      real*8 eh,scal1,dv
      common /psp_cmp_eh/ eh

*     **** external functions ****
      integer  control_version
      external control_version
      real*8   lattice_omega
      external lattice_omega

      if (pawexist) then

      periodic = (control_version().eq.3)

*     ****************************************
*     ***** periodic boundary conditions *****
*     ****************************************
      if (periodic) then
         call Pack_npack(0,npack0)
         ok =BA_push_get(mt_dcpl,npack0,'dng_cmp',dng_cmp(2),dng_cmp(1))
         ok =ok.and.BA_push_get(mt_dcpl,npack0,'vcmp',vcmp(2),vcmp(1))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_pw:out of stack',0,MA_ERR)

         call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
         call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp(1)))
         call Pack_cc_dot(0,dng,dcpl_mb(vcmp(1)),eh)
!$OMP MASTER
         eh = eh*lattice_omega()
!$OMP END MASTER

         ok =        BA_pop_stack(vcmp(2))
         ok = ok.and.BA_pop_stack(dng_cmp(2))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_pw:popping stack',1,MA_ERR)


*     *****************************************
*     ***** aperiodic boundary conditions *****
*     *****************************************
      else
         call D3dB_n2ft3d(1,n2ft3d)
         ok =BA_push_get(mt_dbl,n2ft3d,'rho_cmp',rho_cmp(2),rho_cmp(1))
         ok =ok.and.BA_push_get(mt_dbl,n2ft3d,'vcmp',vcmp(2),vcmp(1))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_pw:out of stack',2,MA_ERR)
         call D3db_nx(1,nx)
         call D3db_ny(1,ny)
         call D3db_nz(1,nz)
         scal1 = 1.0d0/(nx*ny*nz)
         dv    = lattice_omega()*scal1

         call nwpw_compcharge_gen_dn_cmp(ispin,dbl_mb(rho_cmp(1)))

         call Pack_c_unpack(0,dbl_mb(rho_cmp(1)))
         call D3dB_cr_fft3b(1,dbl_mb(rho_cmp(1)))
         call coulomb2_v(dbl_mb(rho_cmp(1)),dbl_mb(vcmp(1)))
         call D3dB_rr_Sum(1,rho,rho(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(rho_cmp(1)))
         call D3dB_rr_dot(1,dbl_mb(vcmp(1)),dbl_mb(rho_cmp(1)),eh)
!$OMP MASTER
         eh = eh*dv
!$OMP END MASTER

         ok =        BA_pop_stack(vcmp(2))
         ok = ok.and.BA_pop_stack(rho_cmp(2))
         if (.not.ok)
     >      call errquit('psp_hartree_cmp_pw:popping stack',2,MA_ERR)
      end if

      end if
      psp_hartree_cmp_pw = eh
      return
      end


*     *************************************************
*     *                                               *
*     *              psp_hartree_atom                 *
*     *                                               *
*     *************************************************
*
*   This routine computes the sum of Watom = Sum(I=1,nionpaw) W1atom^I + W2atom^I+ W3atom^I where
*
*                      //
*    W1atom^I = 0.5 * || (n^I(r)*n^I(r') - ntilde^I(r)*ntilde^I(r'))
*                     || -------------------------------------------  dr dr'
*                    //                    |r-r'|
*
*                      //
*    W2atom^I =   -   || ntilde^I(r)*ncmp^I(r') 
*                     || ----------------------  dr dr'
*                    //         |r-r'|
*
*                      //
*    W3atom^I = -0.5* || ncmp^I(r)*ncmp^I(r') 
*                     || --------------------  dr dr'
*                    //         |r-r'|
*
      real*8 function psp_hartree_atom(ispin,neq,psi1)
      implicit none
      integer    ispin,neq(2)
      complex*16 psi1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** local variables ****
      logical value,sd_function
      integer ii,ia,iii,nproj,npack1,nn,k,l,shift,l_prj
      integer exi(2),sw1(2)
      real*8  Watom

*     ***** external functions ****
      integer  ion_katm,psi_data_get_ptr
      external ion_katm,psi_data_get_ptr
      real*8   nwpw_compcharge_coulomb_e_atom
      external nwpw_compcharge_coulomb_e_atom

      Watom = 0.0d0
      if (pawexist) then

      nn = neq(1)+neq(2)

*     **** allocate local memory ****
      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      if (.not.value) 
     >  call errquit('psp_hartree_atom: out of stack',0,MA_ERR)

      do iii=1,nion_paw
         ii = int_mb(ion_pawtoion(1)+iii-1)
         ia = ion_katm(ii)
         nproj = int_mb(nprj(1)+ia-1)
         if (nproj.gt.0) then

*           **** structure factor ****
            call strfac_pack(1,ii,dcpl_mb(exi(1)))

*           **** generate sw1's and projectors ****
            do l=1,nproj

               shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
               l_prj = int_mb(l_projector(1)+(l-1)
     >                                  + (ia-1)*(nmax_max*lmmax_max))
#ifdef GCC4
               k = iand(l_prj,1)
#else
               k = and(l_prj,1)
#endif
               sd_function = (k.eq.0)

*              **** phase factor does not matter therefore ****
*              **** (-i)^l is the same as (i)^l in the     ****
*              **** Rayleigh scattering formula            ****

*              *** current function is s or d ****
               if (sd_function) then
                  call Pack_tc_Mul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*              *** current function is p or f ****
               else
                  call Pack_tc_iMul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
               end if
               call Pack_cc_indot(1,nn,
     >                          psi1,
     >                          dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                          dbl_mb(sw1(1)+(l-1)*nn))

            end do
            call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))

*           **** paw atom - generate it's atomic density matrix ****
            call psp_gen_density_matrix(ispin,neq,nproj,
     >                                  dbl_mb(sw1(1)),
     >                                  dbl_mb(wtmp(1)))

c*           **** paw atom - generate it's compcharge ***
c            call nwpw_compcharge_gen_Qlm(ii,ia,ispin,nproj,
c     >                                   dbl_mb(wtmp(1)))


            Watom = Watom 
     >            + nwpw_compcharge_coulomb_e_atom(ii,ia,ispin,nproj,
     >                                             dbl_mb(wtmp(1)))
         end if
      end do

*     **** deallocate stack memory ****
      value =           BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     >  call errquit('psp_hartree_atom:popping stack',1,MA_ERR)

      end if
      !call D1dB_SumAll(Watom)


      psp_hartree_atom = Watom
      return
      end

*     ***********************************
*     *					*
*     *	    psp_gen_density_matrix      *
*     *					*
*     ***********************************
*
*   This routine computes the atomic spin density matrix from a atomic spin sw1.
*
      subroutine psp_gen_density_matrix(ispin,ne,nprj,sw1,wmatrix)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  wmatrix(nprj,nprj,ispin)

*     **** local variables ****
      integer i,j,ms,n,n1(2),n2(2)
      real*8  tmp,taskid_i,np_i,icount

      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_np_i(np_i)
      icount = 0

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

c     **************************************************************
c     **** this loop should be restructed to parallelize over n ****
c     **************************************************************
      !call dcopy(ispin*nprj*nprj,0.0d0,0,wmatrix,1)
      call Parallel_shared_vector_zero(.true.,ispin*nprj*nprj,wmatrix)
!$OMP DO 
      do j=1,nprj
         icount = j-1
         if (mod(icount,np_i).eq.taskid_i) then
            do ms=1,ispin
            do n=n1(ms),n2(ms)
               wmatrix(j,j,ms) = wmatrix(j,j,ms) + sw1(n,j)*sw1(n,j)
            end do
            end do
         end if
         icount = icount + 1
      !end do

      !do j=1,nprj
         do i=j+1,nprj
            if (mod(icount,np_i).eq.taskid_i) then
               do ms=1,ispin
               do n=n1(ms),n2(ms)
                  tmp = sw1(n,i)*sw1(n,j)
                  wmatrix(i,j,ms) = wmatrix(i,j,ms) + tmp
                  wmatrix(j,i,ms) = wmatrix(j,i,ms) + tmp
               end do
               end do
            end if
            icount = icount + 1
         end do
      end do
!$OMP END DO 
      !call D1dB_Vector_SumAll(ispin*nprj*nprj,wmatrix)
      call Parallel_Vector_SumAll(ispin*nprj*nprj,wmatrix)

      return
      end



*     ***********************************
*     *					*
*     *	        psp_efg_atoms		*
*     *					*
*     ***********************************
*    This routine computes the efg
*
      subroutine psp_efg_atoms(ispin,ne,psi1,efg)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 efg(3,3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1,nion
      integer ii,ia,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_nion,ion_katm
      integer  psi_data_get_chnk,psi_data_get_ptr
      real*8   lattice_omega,nwpw_xc_energy_atom
      external ion_nion,ion_katm
      external psi_data_get_chnk,psi_data_get_ptr
      external lattice_omega,nwpw_xc_energy_atom

      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      nion = ion_nion()
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      if (.not.value) 
     >  call errquit('psp_efg_atoms: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do ii=1,nion
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if ((int_mb(psp_type(1)+ia-1).ne.4).and.(nproj.gt.0)) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))

        call psp_efg_solve(ia,int_mb(lmax(1)+ia-1),
     >     int_mb(l_projector(1)+(ia-1)*(nmax_max*lmmax_max)),
     >     int_mb(m_projector(1)+(ia-1)*(nmax_max*lmmax_max)),
     >     dbl_mb(psi_data_get_chnk(int_mb(r3_matrix(1)+ia-1))),
     >     ispin,ne,int_mb(nprj(1)+ia-1),dbl_mb(sw1(1)),efg(1,1,ii))

        end if !** nproj>0 **
      end do !** ii **

      value =           BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_efg_atom: popping stack',3,
     &       MA_ERR)
      return 
      end

*     ********************************************************
*     *                                                      *
*     *                psp_efg_solve                         *
*     *                                                      *
*     ********************************************************
******************************************************************************************
*** Warning! This routine is not functioning. Need to rederive the efg for atoms...EJB ***
******************************************************************************************
      subroutine psp_efg_solve(ia,lmax,l_prj,m_prj,
     >                          r3_matrix,
     >                          ispin,ne,nprj,sw1,
     >                          efg)
      implicit none
      integer ia,lmax
      integer l_prj(*)
      integer m_prj(*)
      real*8 r3_matrix(0:lmax,0:lmax)
      integer ispin,ne(2),nprj
      real*8 sw1(ne(1)+ne(2),nprj)
      real*8 efg(3,3)

*     **** external functions ****
      real*8   nwpw_gaunt,lattice_omega
      external nwpw_gaunt,lattice_omega

*     **** local variables ****
      integer i,j,li,lj,mi,mj,l,m,lm,n
      real*8  coeflm(6),tmp,scal,pi,c1,c2,c3

      pi = 4.0d0*datan(1.0d0)

      scal = 1.0d0/lattice_omega()
      do lm=1,6
         coeflm(lm) = 0.0d0
      end do

cccc these formulas need to be redirived!!!
c      do j=1,nprj
c         lj=l_prj(j)
c         mj=m_prj(j)
c         do i=1,nprj
c            li=l_prj(i)
c            mi=m_prj(i)
c
c            tmp = 0.0d0
c            do n=1,(ne(1)+ne(2))
c               tmp = tmp + sw1(n,i)*sw1(n,j)
c            end do
c            tmp = tmp*scal*r3_matrix(li,lj)
c
c            lm = 1
c            do l=0,2,2
c               do m=-l,l
c                  coeflm(lm) = coeflm(lm) 
c     >                       + tmp*nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
c                  lm = lm + 1
c               end do
c            end do
c         end do
c      end do

      c1 = 2.0d0*dsqrt(pi)
      c2 = 6.0d0*dsqrt(pi/15.0d0)
      c3 = 2.0d0*dsqrt(pi/5.0d0)

      !*** 2*sqrt(pi)*(l=0,m=0) + 6*sqrt(pi/15)*(l=2,m=2) + 2*sqrt(pi/5)*(l=2,m=0) ****
      efg(1,1) = efg(1,1) + c1*coeflm(1) + c2*coeflm(6) + c3*coeflm(4)

      !*** 6*sqrt(pi/15)*(l=2,m=-2)***
      efg(2,1) = efg(2,1) + c1*coeflm(2)
      efg(1,2) = efg(1,2) + efg(2,1)

      !*** 6*sqrt(pi/15)*(l=2,m=1)***
      efg(3,1) = efg(3,1) + c1*coeflm(5)
      efg(1,3) = efg(1,3) + efg(3,1)

      !*** 2*sqrt(pi)*(l=0,m=0) - 6*sqrt(pi/15)*(l=2,m=2) + 2*sqrt(pi/5)*(l=2,m=0) ****
      efg(2,2) = efg(2,2) + c1*coeflm(1) -c2*coeflm(6) + c3*coeflm(4)

      !*** 6*sqrt(pi/15)*(l=2,m=-1)***
      efg(3,2) = efg(3,2) + c2*coeflm(3)
      efg(2,3) = efg(2,3) + efg(3,2)

      !*** 4*(l=2,m=0)+ 2*sqrt(pi)(l=0,m=0) ***
      efg(3,3) = efg(3,3) + 4.0d0*coeflm(4) + c1*coeflm(1)

      return
      end 
      





*     *************************************************
*     *                                               *
*     *              psp_dE_ncmp_vloc_Qlm_test        *
*     *                                               *
*     *************************************************

      subroutine psp_dE_ncmp_vloc_Qlm_test(ispin,move,fion)
      implicit none
      integer ispin
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     **** local variables ****
      logical ok,periodic
      integer npack0,dng_cmp(2),dng_cmp_smooth(2),vl1(2),vl_notpaw(2)
      integer n2ft3d,rho_cmp(2),rho_cmp_smooth(2),vl2(2),r_grid(2)
      integer nx,ny,nz
      real*8 eh,e1,e2,eh0,eh1,eh2,scal1,dv
      real*8 dum1(3),dum2(3)

*     **** external functions ****
      integer  control_version
      external control_version
      real*8   nwpw_compcharge_E_multipole_zv_ee,lattice_omega
      external nwpw_compcharge_E_multipole_zv_ee,lattice_omega


      eh = 0.0d0
      if (pawexist) then

      periodic = (control_version().eq.3)

*     *************************************
*     **** Periodic Boundary Condtions ****
*     *************************************
      if (periodic) then
         call Pack_npack(0,npack0)
         ok =        BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                           dng_cmp(2),dng_cmp(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                           dng_cmp_smooth(2),dng_cmp_smooth(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vl1',
     >                           vl1(2),vl1(1))
         ok = ok.and.BA_push_get(mt_dcpl,npack0,'vl_notpaw',
     >                           vl_notpaw(2),vl_notpaw(1))
         if (.not.ok)
     >      call errquit('psp_ncmp_vloc:out of stack',0,MA_ERR)

*        **** Using pw grid only ***
         if (use_grid_cmp) then
c            call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))
            call v_local(dcpl_mb(vl1(1)),
     >               .false.,
     >               dum1,dum2)
c            call Pack_cc_dot(0,dcpl_mb(dng_cmp(1)),
c     >                         dcpl_mb(vl1(1)),eh)

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_pw(ispin,
     >                                               dcpl_mb(vl1(1)),
     >                                               move,fion)

*        **** Using gaussian two-electron integrals and pw grid ***
         else
            call v_local_seperate_paw(dcpl_mb(vl1(1)),
     >                                dcpl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_v_cmp_smooth(dbl_mb(zv(1)),
     >                                            dcpl_mb(dng_cmp(1)))

            call Pack_cc_Sub2(0,dcpl_mb(dng_cmp(1)),
     >                          dcpl_mb(vl1(1)))
            call Pack_cc_Sum2(0,dcpl_mb(dng_cmp(1)),
     >                          dcpl_mb(vl_notpaw(1)))

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_test(
     >                                            ispin,dbl_mb(zv(1)),
     >                                            dcpl_mb(vl_notpaw(1)),
     >                                            dcpl_mb(vl1(1)),
     >                                            move,fion)
            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dcpl_mb(dng_cmp(1)),
     >                                       dcpl_mb(dng_cmp_smooth(1)))
            call f_local_seperate_paw(dcpl_mb(dng_cmp_smooth(1)),
     >                                dcpl_mb(dng_cmp(1)),
     >                                fion)
            call Pack_cc_Sub2(0,dcpl_mb(dng_cmp_smooth(1)),
     >                          dcpl_mb(dng_cmp(1)))
            call nwpw_compcharge_gen_f_cmp_smooth(dbl_mb(zv(1)),
     >                                            dcpl_mb(dng_cmp(1)),
     >                                            fion)
         end if
         ok =        BA_pop_stack(vl_notpaw(2))
         ok = ok.and.BA_pop_stack(vl1(2))
         ok = ok.and.BA_pop_stack(dng_cmp_smooth(2))
         ok = ok.and.BA_pop_stack(dng_cmp(2))
         if (.not.ok)
     >      call errquit('psp_vloc_residual:pop stack',1,MA_ERR)

*     **************************************
*     **** APeriodic Boundary Condtions ****
*     **************************************
      else
         call D3dB_n2ft3d(1,n2ft3d)
         ok =         BA_push_get(mt_dbl,n2ft3d,'rho_cmp',
     >                            rho_cmp(2),rho_cmp(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'rho_cmp_smooth',
     >                           rho_cmp_smooth(2),rho_cmp_smooth(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl1',
     >                           vl1(2),vl1(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl_notpaw',
     >                           vl_notpaw(2),vl_notpaw(1))
         ok = ok.and.BA_push_get(mt_dbl,n2ft3d,'vl2',
     >                           vl2(2),vl2(1))
         ok = ok.and.BA_push_get(mt_dbl,3*n2ft3d,'rgrid_cmp',
     >                           r_grid(2),r_grid(1))
         if (.not.ok)
     >      call errquit('psp_ncmp_vloc:out of stack',2,MA_ERR)
         call lattice_r_grid(dbl_mb(r_grid(1)))
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/(nx*ny*nz)
         dv = scal1*lattice_omega()

*        **** Using pw grid only ***
         if (use_grid_cmp) then
            call v_local(dbl_mb(vl1(1)),
     >               .false.,
     >               dum1,dum2)
            call v_lr_local(dbl_mb(r_grid(1)),
     >                      dbl_mb(rho_cmp(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(rho_cmp(1)))
            call Pack_c_pack(0,dbl_mb(rho_cmp(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(rho_cmp(1)))
            call Pack_cc_Sum2(0,dbl_mb(rho_cmp(1)),dbl_mb(vl1(1)))

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_pw(ispin,
     >                                               dbl_mb(vl1(1)),
     >                                               move,fion)

*        **** Using gaussian two-electron integrals and pw grid ***
         else

*           **** long-range terms ****
            call v_lr_local_seperate_paw(dbl_mb(r_grid(1)),
     >                                   dbl_mb(vl1(1)),
     >                                   dbl_mb(vl_notpaw(1)))
            call nwpw_compcharge_gen_vlr_cmp_smooth(dbl_mb(zv(1)),
     >                                           dbl_mb(r_grid(1)),
     >                                           dbl_mb(rho_cmp(1)))

            call D3dB_rr_Sum2(1,dbl_mb(rho_cmp(1)),dbl_mb(vl_notpaw(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(vl_notpaw(1)))
            call Pack_c_pack(0,dbl_mb(vl_notpaw(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(vl_notpaw(1)))

            call D3dB_rr_Sub2(1,dbl_mb(rho_cmp(1)),dbl_mb(vl1(1)))
            call D3dB_rc_pfft3f(1,0,dbl_mb(vl1(1)))
            call Pack_c_pack(0,dbl_mb(vl1(1)))
            call Pack_c_SMul1(0,dv,dbl_mb(vl1(1)))

*           **** short-range terms ****
            call v_local_seperate_paw(dbl_mb(rho_cmp_smooth(1)),
     >                                dbl_mb(rho_cmp(1)))
            call Pack_cc_Sum2(0,
     >                        dbl_mb(rho_cmp(1)),
     >                        dbl_mb(vl_notpaw(1)))
            call Pack_cc_Sum2(0,
     >                        dbl_mb(rho_cmp_smooth(1)),
     >                        dbl_mb(vl1(1)))

            call nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_test(
     >                                            ispin,dbl_mb(zv(1)),
     >                                            dbl_mb(vl_notpaw(1)),
     >                                            dbl_mb(vl1(1)),
     >                                            move,fion)

               call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dbl_mb(rho_cmp(1)),
     >                                       dbl_mb(rho_cmp_smooth(1)))
               call f_local_seperate_paw(dbl_mb(rho_cmp_smooth(1)),
     >                                   dbl_mb(rho_cmp(1)),
     >                                   fion)

               call Pack_c_unpack(0,dbl_mb(rho_cmp(1)))
               call D3dB_cr_pfft3b(1,0,dbl_mb(rho_cmp(1)))

               call Pack_c_unpack(0,dbl_mb(rho_cmp_smooth(1)))
               call D3dB_cr_pfft3b(1,0,dbl_mb(rho_cmp_smooth(1)))

               call f_lr_local_seperate_paw(dbl_mb(r_grid(1)),
     >                                      dbl_mb(rho_cmp_smooth(1)),
     >                                      dbl_mb(rho_cmp(1)),
     >                                      fion)

               call D3dB_rr_Sub2(1,dbl_mb(rho_cmp_smooth(1)),
     >                             dbl_mb(rho_cmp(1)))
               call nwpw_compcharge_gen_f_lr_cmp_smooth(dbl_mb(zv(1)),
     >                                               dbl_mb(r_grid(1)),
     >                                               dbl_mb(rho_cmp(1)),
     >                                               fion)
         end if
         ok =        BA_pop_stack(r_grid(2))
         ok = ok.and.BA_pop_stack(vl2(2))
         ok = ok.and.BA_pop_stack(vl_notpaw(2))
         ok = ok.and.BA_pop_stack(vl1(2))
         ok = ok.and.BA_pop_stack(rho_cmp_smooth(2))
         ok = ok.and.BA_pop_stack(rho_cmp(2))
         if (.not.ok)
     >      call errquit('psp_vloc_ncmp_vloc:pop stack',3,MA_ERR)

      end if


      end if

      return
      end
