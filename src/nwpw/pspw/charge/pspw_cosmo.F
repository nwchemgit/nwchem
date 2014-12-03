c
c $Id$
c
  
***********************************************************
*                                                         *
*                   PSPW-COSMO module                     *
*                                                         *
*                                                         *
*   -- developed by Eric J. Bylaska on September 1, 2012  *
*                                                         *
***********************************************************
*
*
*  

*     **********************************
*     *	                               *
*     *          pspw_cosmo_init       *
*     *                                *
*     **********************************
 
      subroutine pspw_cosmo_init()
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "pspw_cosmo.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value

*     **** external functions ****
      logical  control_cosmo
      external control_cosmo

      cosmo_found = control_cosmo()

      if (cosmo_found) then
         value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,'fcut',fcut(2),fcut(1))
         value = value.and.
     >        BA_alloc_get(mt_dbl,3*n2ft3d,'gradfcut',
     >                     gradfcut(2),gradfcut(1))
         value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,'lapfcut',
     >                     lapfcut(2),lapfcut(1))

         value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,'vcond',vcond(2),vcond(1))
         value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,'Avcond',Avcond(2),Avcond(1))

      end if

      return
      end



*     **********************************
*     *	                               *
*     *        pspw_cosmo_end          *
*     *                                *
*     **********************************
      subroutine pspw_cosmo_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_cosmo.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

      if (cosmo_found) then
         value =           BA_free_heap(fcut(2))
         value = value.and.BA_free_heap(vcond(2))
      end if

      return
      end


*     **********************************
*     *	                               *
*     *         pspw_cosmo_found       *
*     *                                *
*     **********************************
      logical function pspw_cosmo_found()
      implicit none

#include "pspw_cosmo.fh"
 

      pspw_cosmo_found = cosmo_found
      return
      end



*     **********************************
*     *	                               *
*     *     pspw_cosmo_nq_total        *
*     *                                *
*     **********************************
      integer function pspw_cosmo_nq_total()
      implicit none

#include "pspw_cosmo.fh"

      pspw_cosmo_nq_total = nq_total
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_cosmo_nq_ion          *
*     *                                *
*     **********************************
      integer function pspw_cosmo_nq_ion(ii)
      implicit none
      integer ii

#include "bafdecls.fh"
#include "pspw_cosmo.fh"
#include "errquit.fh"

      pspw_cosmo_nq_ion = int_mb(nq_ion(1)+ii-1)
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_cosmo                 *
*     *                                *
*     **********************************

      subroutine pspw_cosmo(n2ft3d,r_grid,vc,vcosmo)
      impicit none
      integer n2ft3d
      real*8 r_grid(3,n2ft3d),vc(n2ft3d),vcosmo(n2ft3d)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer grx(2),gry(2),grz(2)

*     **** external functions ****
      integer  G_indx
      external G_indx


      call pspw_cosmo_vions(n2ft3d,0.1d0,r_grid,vcosmo)
      call daxpy(n2ft3d,1.0d0,vc,1,vcosmo,1)

      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'fcut',fcut(2),grx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'vconductor',vconductor(2),vconductor(1))

      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grx',grx(2),grx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'gry',gry(2),gry(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grz',grz(2),grz(1))
      if (.not. value) 
     >   call errquit('pspw_cosmo:out of stack memory',0,MA_ERR)


      call D3dB_r_Zero_Ends(1,vcosmo)
      call D3dB_r_SMul1(1,scal1,vcosmo)
      call D3dB_rc_fft3f(1,vcosmo)
      call mask_C(0,vcosmo)

c     **** calculate  gradient vcosmo ****
      call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                   vcosmo,
     >                   dbl_mb(grx(1)))
      call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                   vcosmo,
     >                   dbl_mb(gry(1)))
      call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                   vcosmo,
     >                   dbl_mb(grz(1)))
      call D3dB_cr_fft3b(1,dbl_mb(grx(1)))
      call D3dB_cr_fft3b(1,dbl_mb(gry(1)))
      call D3dB_cr_fft3b(1,dbl_mb(grz(1)))




      return
      end



*     ***********************************
*     *				        *
*     *	        pspw_cosmo_vions        *
*     *					*
*     ***********************************
      subroutine pspw_cosmo_vions(n2ft3d,rcut,r_grid,vlr_out)
      implicit none
      integer n2ft3d
      real*8 rcut
      real*8 r_grid(3,*)
      real*8 vlr_out(*)

*     **** Error function parameters ****
      real*8 xerf,yerf

*     **** local variables ****
      integer taskid_j,np_j
      integer i,j
      real*8 x,y,z,q,c,r,sqrt_pi

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   ion_rion,psp_zv,util_erf
      external ion_nion,ion_katm
      external ion_rion,psp_zv,util_erf

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      call dcopy(n2ft3d,0.0d0,0,vlr_out,1)
      do j=1,ion_nion()
         if (mod(j-1,np_j).eq.taskid_j) then
            x = ion_rion(1,j)
            y = ion_rion(2,j)
            z = ion_rion(3,j)
            q = -psp_zv(ion_katm(j))
            c = 1.0d0/rcut
            do i=1,n2ft3d
               r = dsqrt( (r_grid(1,i)-x)**2
     >                  + (r_grid(2,i)-y)**2
     >                  + (r_grid(3,i)-z)**2)
               if (r .gt. 1.0d-15) then
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
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d,vlr_out)
      call nwpw_timing_end(5)
      return
      end

*     ***********************************
*     *				        *
*     *	        pspw_cosmo_gen_fcut     *
*     *					*
*     ***********************************
      subroutine pspw_cosmo_gen_fcut(rcut,nion,rcosmo,rion,
     >                               n2ft3d,r_grid,fcut)
      implicit none
      real*8  rcut
      integer nion
      real*8  rcosmo(*),rion(3,*)
      integer n2ft3d
      real*8  r_grid(3,n2ft3d),fcut(n2ft3d)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer i
      integer tmp(2)

*     **** external functions ****
      real*8   util_erf
      external util_erf

      if (.not.BA_push_get(mt_dbl, n2ft3d,'tmp',tmp(2),tmp(1)))
     >   call errquit('pspw_cosmo_gen_fcut:out of stack',0,MA_ERR)

      call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp(1)),1)
      do i=1,nion
         call pspw_cosmo_gen_fcut_ion(rcut,rcosmo(i),rion(1,i),
     >                                n2ft3d,r_grid,dbl_mb(tmp(1)))
      end do

      call dcopy(n2ft3d,0.0d0,0,fcut,1)
      do i=1,n2ft3d
         fcut(i) = util_erf(dbl_mb(tmp(1)+i-1)/rcut)
      end do

      if (.not.BA_pop_stack(tmp(2))) 
     >   call errquit('pspw_cosmo_gen_fcut:popping stack',1,MA_ERR)
      return
      end


*     ***********************************
*     *				        *
*     *	    pspw_cosmo_add_fcut_ion     *
*     *					*
*     ***********************************

      subroutine pspw_cosmo_add_fcut_ion(rcut,rsphere,rion,
     >                                   n2ft3d,r_grid,fcut_ion)
      implicit none
      real*8  rcut,rsphere,rcut,rion(3)
      integer n2ft3d
      real*8  r_grid(3,n2ft3d)
      real*8  fcut_ion(n2ft3d)

*     **** local variables ****
      integer i
      real*8 x,y,z,r,u

*     **** external functions ****
      real*8   util_erfc
      external util_erfc
      
      do i=1,n2ft3d
         x = rion(1) - r_grid(1,i)
         y = rion(2) - r_grid(2,i)
         z = rion(3) - r_grid(3,i)
         r = dsqrt(x*x + y*y + z*z)
         u = (r-rsphere)/rcut
         fcut_ion(i) = fcut_ion(i) + 0.5d0*util_erfc(u)
      end do
      return
      end
      
c*     ***********************************
c*     *				        *
c*     *	    pspw_cosmo_gen_dfcut_ion    *
c*     *					*
c*     ***********************************
c      
c      subroutine pspw_cosmo_gen_dfcut_ion(rcut,rsphere,rion,
c     >                                   n2ft3d,r_grid,dfcut_ion)
c      implicit none
c      real*8  rcut,rsphere,rion(3)
c      integer n2ft3d
c      real*8  r_grid(3,n2ft3d)
c      real*8  dfcut_ion(3,n2ft3d)
c      
c*     **** local variables ****
c      integer i
c      real*8 x,y,z,r,u,df,const
c
c*     **** external functions ****
c      real*8   util_erfc
c      external util_erfc
c
c      const = -1.0d0/dsqrt(4.0d0*datan(1.0d0))
c      const = const / rcut
c
c      call dcopy(3*n2ft3d,0.0d0,0,dfcut_ion,1)
c      do i=1,n2ft3d
c         x = rion(1) - r_grid(1,i)
c         y = rion(2) - r_grid(2,i)
c         z = rion(3) - r_grid(3,i)
c         r = dsqrt(x*x + y*y + z*z)
c         u = (r-rsphere)/rcut
c         df = const*dexp(-u*u)/r
c         
c         dfcut_ion(1,i) = df*x
c         dfcut_ion(1,i) = df*y
c         dfcut_ion(1,i) = df*z
c
c      end do
c      return
c      end

      subroutine pspw_cosmo_cg(n2ft3d,f1,f2xyz,f3,gamma)
      implicit none
      integer n2ft3d
      real*8 f1(n2ft3d),f2xyz(3,n2ft3d),f3(n2ft3d),b(n2ft3d)

      call pspw_cosmo_A(n2ft3d,f1,f2xyz,f3,gamma,gradp,lapp,Ap)
      call dcopy(n2ft3d,b,1,r,1)
      call daxpy(n2ft3d,-1.0d0,Ap,1,r,1)
      call dcopy(n2ft3d,r,1,p,1)
      call D3dB_rr_dot(1,r,r,rsold)

      done = .false.
      while (.not. done) do
         call pspw_cosmo_A(n2ft3d,f1,f2xyz,f3,p,gradp,lapp,Ap)
         call D3dB_rr_dot(1,p0,Ap0,x)
         alpha = rsold/x
         call daxpy(n2ft3d, alpha,p,1,r, 1)
         call daxpy(n2ft3d,-alpha,p,1,Ap,1)
         call D3dB_rr_dot(1,r,r,rsnew)
         done = (rsnew.lt.tol)
         if (.not.done) then
            call dscal(n2ft3d,(rsnew/rsold),p,1)
            call daxpy(n2ft3d,1.0d0,r,1,p,1)
            rsold = rsnew
         end if
      end do

      return
      end

*     ***********************************
*     *				        *
*     *	    pspw_cosmo_Agamma           *
*     *					*
*     ***********************************
      subroutine pspw_cosmo_Agamma(n2ft3d,tmpG,f1,f1x,f1y,f1z,lapf1,
     >                             tmp,tmp2,
     >                             gamma,Agamma)
      implicit none
      integer n2ft3d
      real*8 tmpG(n2ft3d/2,4)
      real*8 f1(n2ft3d)
      real*8 f1x(n2ft3d)
      real*8 f1y(n2ft3d)
      real*8 f1z(n2ft3d)
      real*8 lapf1(n2ft3d)
      real*8 tmp1(n2ft3d)
      real*8 tmp2(n2ft3d)
      real*8 gamma(n2ft3d)
      real*8 Agamma(n2ft3d)

      do i=1,n2ft3d
         Agamma(i) = lapf1(i)*gamma(i)
      end do

      call dcopy(n2ft3d,gamma,1,tmp1,1)
      call D3dB_r_Zero_Ends(1,tmp1)
      call D3dB_r_SMul1(1,scal1,tmp1)
      call D3dB_rc_fft3f(1,tmp1)
      call mask_C(0,tmp1)

c     **** calculate  gradient gamma ****
      call D3dB_ic_Mul(1,tmpG(1,1),tmp1,tmp2)
      call D3dB_cr_fft3b(1,tmp2)
      do i=1,n2ft3d
         Agamma(i) = 2.0d0*f1x(i)*tmp2(i)
      end do

      call D3dB_ic_Mul(1,tmpG(1,2),tmp1,tmp2)
      call D3dB_cr_fft3b(1,tmp2)
      do i=1,n2ft3d
         Agamma(i) = Agamma(i) + 2.0d0*f1y(i)*tmp2(i)
      end do

      call D3dB_ic_Mul(1,tmpG(1,3),tmp1,tmp2)
      call D3dB_cr_fft3b(1,tmp2)
      do i=1,n2ft3d
         Agamma(i) = Agamma(i) + 2.0d0*f1z(i)*tmp2(i)
      end do

      call D3dB_tc_Mul(1,tmpG(1,4),tmp1,tmp2)
      call D3dB_cr_fft3b(1,tmp2)
      do i=1,n2ft3d
         Agamma(i) = Agamma(i) + f1(i)*tmp2(i)
      end do

      return
      end 
