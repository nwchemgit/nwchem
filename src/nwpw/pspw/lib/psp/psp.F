#define TCGMSG

      subroutine psp_init()      
*
* $Id: psp.F,v 1.1 2002-01-14 21:06:28 bylaska Exp $
*
      implicit none

#include "mafdecls.fh"


*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal


*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      integer npack1,npack0
      logical value

*     **** external functions *****
      integer  ion_nkatm
      external ion_nkatm
 

      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      npsp = ion_nkatm()

      value = MA_alloc_get(mt_dbl,(npsp*npack0),
     >                    'vl',vl(2),vl(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(lmmax_max*npsp*npack1),
     >                    'vnl',vnl(2),vnl(1))
      value = value.and.
     >         MA_alloc_get(mt_dbl,(lmmax_max*npsp),
     >                    'vnlnrm',vnlnrm(2),vnlnrm(1))
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
     >        MA_alloc_get(mt_dbl,(npsp),'rlocal',rlocal(2),rlocal(1))

      if (.not. value) call errquit('out of heap memory',0)

      call dcopy(npsp*npack0,          0.0d0,0,dbl_mb(vl(1)), 1)
      call dcopy(lmmax_max*npsp*npack1,0.0d0,0,dbl_mb(vnl(1)),1)
      call dcopy(lmmax_max*npsp,       0.0d0,0,dbl_mb(vnlnrm(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(zv(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(amass(1)),1)
      call dcopy(npsp*(lmax_max+1),    0.0d0,0,dbl_mb(rc(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(rlocal(1)),1)

*     **** allocate semicore data ****
      call semicore_init()

      return
      end




      subroutine psp_end()      
      implicit none

#include "mafdecls.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal

*     **** common blocks ****
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      logical value

*     **** external functions ****

*     **** deallocate semicore data ****
      call semicore_end()

      value = MA_free_heap(vl(2))
      value = MA_free_heap(vnl(2))
      value = MA_free_heap(vnlnrm(2))
      value = MA_free_heap(zv(2))
      value = MA_free_heap(amass(2))
      value = MA_free_heap(rc(2))
      value = MA_free_heap(lmmax(2))
      value = MA_free_heap(lmax(2))
      value = MA_free_heap(locp(2))
      value = MA_free_heap(rlocal(2))

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


*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      psp_zv = dbl_mb(zv(1)+ia-1)
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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      psp_atom = atom(ia)
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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx),lmax(2),locp(2)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      psp_lmmax = int_mb(lmmax(1)+ia-1)
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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

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


*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      psp_npsp = npsp
      return
      end


c*     ***********************************
c*     *					*
c*     *	 	   v_local  		*
c*     *					*
c*     ***********************************
c
c      subroutine v_local(vl_out)
c      implicit none
c      complex*16 vl_out(*)
c
c#include "mafdecls.fh"
c
c*     **** common blocks ****
cc     real*8     vl(nfft3d,nkatmx)
cc     real*8     vnl(nfft3d,16,nkatmx)
cc     real*8     vnlnrm(16,nkatmx)
cc     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
cc     integer    lmmax(nkatmx)
c      integer    vl(2)
c      integer    vnl(2)
c      integer    vnlnrm(2)
c      integer    zv(2),amass(2),rc(2)
c      integer    lmmax(2),lmax(2)
c      integer    npsp
c   
c      integer nkatmx,lmmax_max,lmax_max
c      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
c      character*2 atom(nkatmx)
c      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
c     >                     lmmax,lmax,npsp,atom
c
c
c
c*     *** local variables ***
c      integer nfft3d
c      integer i,ii,ia,q
c      integer exi(2),vtmp(2)
c      logical value
c
c*     **** external functions ****
c      integer  ion_nion,ion_katm
c      external ion_nion,ion_katm
c
c      call D3dB_nfft3d(1,nfft3d)
c
c      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
c      value = MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
c
c      call dcopy((2*nfft3d),0.0d0,0,vl_out,1)
c      do ii=1,ion_nion()
c        ia=ion_katm(ii)
c
c*       **** structure factor and local pseudopotential ****
c        call strfac(ii,dcpl_mb(exi(1)))
c        call mask_C(0, dcpl_mb(exi(1)))
c        
c*       **** add to local psp ****
c        call D3dB_tc_MUL(1,dbl_mb(vl(1)+nfft3d*(ia-1)),
c     >                   dcpl_mb(exi(1)),
c     >                   dcpl_mb(vtmp(1)))
c        call D3dB_cc_Sum(1,vl_out,dcpl_mb(vtmp(1)),vl_out)
c
cc       if (move) then
cc       
cc          do i=1,nfft3d
cc             xtmp(i) = dimag(dng(i))*dble(vtmp(i))
cc    >                     - dble(dng(i))*dimag(vtmp(i))
cc          end do
cc          call D3dB_tt_dot(1,G(1,1),xtmp,fion(1,ii))
cc          call D3dB_tt_dot(1,G(1,2),xtmp,fion(2,ii))
cc          call D3dB_tt_dot(1,G(1,3),xtmp,fion(3,ii))
cc        
cc       end if
c        
c
c      end do
c      value = MA_pop_stack(vtmp(2))
c      value = MA_pop_stack(exi(2))
c
c
c      return 
c      end

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


*     **** Error function parameters ****
      real*8  yerf,verf
c     real*8 c1,c2,c3,c4,c5,c6,yerf,fterf,verf
c     parameter (c1=0.07052307840d0,c2=0.04228201230d0)
c     parameter (c3=0.00927052720d0)
c     parameter (c4=0.00015201430d0,c5=0.00027656720d0)
c     parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer i,j,np1,np2,np3,n2ft3d
      real*8 x,y,z,q,c,r,sqrt_pi,dv,v,rx,ry,rz,fx,fy,fz

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf
      external ion_nion,ion_katm
      external lattice_omega,ion_rion,psp_rlocal,psp_zv,util_erf

      call D3dB_n2ft3d(1,n2ft3d)

*     **** constants ****
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      call D3dB_nx(1,np1)
      call D3dB_nx(1,np2)
      call D3dB_nx(1,np3)
      dv = lattice_omega()/dble(np1*np2*np3)

      do j=1,ion_nion()
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
c             fterf = (1.0d0
c    >              + yerf*(c1 + yerf*(c2
c    >              + yerf*(c3 + yerf*(c4
c    >              + yerf*(c5 + yerf*c6))))))**4
c             verf = (1.0d0 - 1.0d0/fterf**4)
              verf = util_erf(yerf)
              v    = q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                 - verf)/r**3
            else
              v = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do


*        fion(1,j) = fion(1,j) - ddot(n2ft3d,rho,1,gv(1,1),3)*dv
*        fion(2,j) = fion(2,j) - ddot(n2ft3d,rho,1,gv(2,1),3)*dv
*        fion(3,j) = fion(3,j) - ddot(n2ft3d,rho,1,gv(3,1),3)*dv

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fion(1,j) = fion(1,j) - fx*dv
         fion(2,j) = fion(2,j) - fy*dv
         fion(3,j) = fion(3,j) - fz*dv
      end do

      return
      end





*     ***********************************
*     *				                 	*
*     *	 	   v_lr_local  				*
*     *									*
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
      integer i,j,n2ft3d
      real*8 x,y,z,q,c,r,sqrt_pi

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   ion_rion,psp_rlocal,psp_zv,util_erf
      external ion_nion,ion_katm
      external ion_rion,psp_rlocal,psp_zv,util_erf

      call D3dB_n2ft3d(1,n2ft3d)

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      call dcopy(n2ft3d,0.0d0,0,vlr_out,1)

      do j=1,ion_nion()
         x = ion_rion(1,j)
         y = ion_rion(2,j)
         z = ion_rion(3,j)
         q = -psp_zv(ion_katm(j))
         c = 1.0d0/psp_rlocal(ion_katm(j))

         do i=1,n2ft3d
            r = dsqrt( (r_grid(1,i)-x)**2
     >               + (r_grid(2,i)-y)**2
     >               + (r_grid(3,i)-z)**2)

            if (r .gt. 1.0d-15) then
              xerf=r*c
c             yerf = (1.0d0
c    >              + xerf*(c1 + xerf*(c2
c    >              + xerf*(c3 + xerf*(c4
c    >              + xerf*(c5 + xerf*c6))))))**4
c             yerf = (1.0d0 - 1.0d0/yerf**4)
              yerf = util_erf(xerf)
              vlr_out(i) = vlr_out(i) + (q/r)*yerf
c             vlr_out(i) = vlr_out(i) + (q/r)*erf(r*c)
            else
              vlr_out(i) = vlr_out(i) + 2.0d0*q*c/sqrt_pi
            end if

         end do
      end do

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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom



*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      call nwpw_timing_start(5)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))

      call dcopy((2*npack0),0.0d0,0,vl_out,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(0,dcpl_mb(exi(1)))
        
*       **** add to local psp ****
        call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))
        call Pack_cc_Sum(0,vl_out,dcpl_mb(vtmp(1)),vl_out)


        if (move) then
       
          do i=1,npack0
             dbl_mb(xtmp(1)+i-1) 
     >          = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >           - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
          end do
         call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
         call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
         call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))
   
       end if
        

      end do
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

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

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      call nwpw_timing_start(5)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(0,dcpl_mb(exi(1)))
        
*       **** add to local psp ****
        call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))

         do i=1,npack0
           dbl_mb(xtmp(1)+i-1) 
     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
        end do

        call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
        call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
        call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))
 
      end do
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

      call nwpw_timing_end(5)
      return 
      end




*     ***********************************
*     *					*
*     *	 	   v_nonlocal  		*
*     *					*
*     ***********************************

      subroutine v_nonlocal(ispin,ne,psi1,psi2,move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


*     *** local variables ***
      integer nfft3d,G(3),npack1
      integer i,ii,ia,l,n,nn
      real*8  omega
      complex*16 ctmp
      integer exi(2),vtmp(2),xtmp(2),sw1(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,G_indx
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn,'sw1',sw1(2),sw1(1))
      value = value.and.
     >      MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))

      omega = lattice_omega()

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(1,dcpl_mb(exi(1)))

        do l=1,int_mb(lmmax(1)+ia-1)

*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (is_sORd(l,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
              call Pack_tc_Mul(1,dbl_mb(vnl(1)+(l-1)*npack1
     >                           +(ia-1)*npack1*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Pack_itc_Mul(1,dbl_mb(vnl(1)+(l-1)*npack1
     >                           +(ia-1)*npack1*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
c             do i=1,npack1
c                dcpl_mb(vtmp(1)+i-1) 
c    >                     = dbl_mb(vnl(1) + (i-1)
c    >                                     + (l-1)*npack1
c    >                          + (ia-1)*npack1*lmmax_max)
c    >                * dcmplx(-dimag(dcpl_mb(exi(1)+i-1)),
c    >                          dble(dcpl_mb(exi(1)+i-1)))
c             end do
           end if


*         **** do kleinman-bylander multiplication ****
          call Pack_cc_ndot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)))

          do n=1,nn
          
             dbl_mb(sw1(1)+n-1) = dbl_mb(sw1(1)+n-1)
     >                      /(omega*dbl_mb(vnlnrm(1)
     >                         +(l-1)+(ia-1)*lmmax_max))
             call Pack_cc_daxpy(1,(-dbl_mb(sw1(1)+n-1)),
     >                          dcpl_mb(vtmp(1)),
     >                          psi2(1+(n-1)*npack1)) 

          end do

          if (move) then
             do n=1,nn
                if (ispin.eq.1) 
     >            dbl_mb(sw1(1)+n-1)=2.0d0*dbl_mb(sw1(1)+n-1)

                do i=1,npack1
                   ctmp = psi1(i+(n-1)*npack1)
     >                  *dconjg(dcpl_mb(vtmp(1)+i-1))
                   dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                end do
                call Pack_tt_idot(1,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+1+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+2+3*(n-1)))

             end do

             call D3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
             do n=1,nn
                fion(1,ii) = fion(1,ii) 
     >                     + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                            *dbl_mb(sum(1)+3*(n-1))
                fion(2,ii) = fion(2,ii) 
     >                     + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                            *dbl_mb(sum(1)+1+3*(n-1))
                fion(3,ii) = fion(3,ii) 
     >                     + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                            *dbl_mb(sum(1)+2+3*(n-1))
             end do
          end if
       
        
        end do

      end do

      value = MA_pop_stack(sum(2))
      value = MA_pop_stack(sw1(2))
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

      call nwpw_timing_end(6)
      return 
      end



*     ***********************************
*     *					*
*     *	 	   f_vnonlocal 		*
*     *					*
*     ***********************************

      subroutine f_vnonlocal(ispin,ne,psi1,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom



*     *** local variables ***
      integer nfft3d,G(3),npack1
      integer i,ii,ia,l,n,nn
      real*8  omega
      complex*16 ctmp
      integer exi(2),vtmp(2),xtmp(2),sw1(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,G_indx
      external lattice_omega

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn,'sw1',sw1(2),sw1(1))
      value = value.and.
     >      MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))

      omega = lattice_omega()

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(1,dcpl_mb(exi(1)))

        do l=1,int_mb(lmmax(1)+ia-1)

*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (is_sORd(l,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
              call Pack_tc_Mul(1,dbl_mb(vnl(1)+(l -1)*npack1
     >                           +(ia-1)*npack1*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              do i=1,npack1
                 dcpl_mb(vtmp(1)+i-1) 
     >                     = dbl_mb(vnl(1) + (i-1)
     >                                     + (l-1)*npack1
     >                          + (ia-1)*npack1*lmmax_max)
     >                * dcmplx(-dimag(dcpl_mb(exi(1)+i-1)),
     >                          dble(dcpl_mb(exi(1)+i-1)))
              end do
           end if

*         **** do kleinman-bylander multiplication ****
          call Pack_cc_ndot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)))

          do n=1,nn
             dbl_mb(sw1(1)+n-1) = dbl_mb(sw1(1)+n-1)
     >                      /(omega*dbl_mb(vnlnrm(1)
     >                         +(l-1)+(ia-1)*lmmax_max))
 
          end do

          
          do n=1,nn
             if (ispin.eq.1) 
     >         dbl_mb(sw1(1)+n-1)=2.0d0*dbl_mb(sw1(1)+n-1)

               do i=1,npack1
                  ctmp = psi1(i+(n-1)*npack1)
     >                 *dconjg(dcpl_mb(vtmp(1)+i-1))
                  dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
               end do

               call Pack_tt_idot(1,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+3*(n-1)))
               call Pack_tt_idot(1,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                           dbl_mb(sum(1)+1+3*(n-1)))
               call Pack_tt_idot(1,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                           dbl_mb(sum(1)+2+3*(n-1)))

          end do

          call D3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
          do n=1,nn
             fion(1,ii) = fion(1,ii) 
     >                  + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                         *dbl_mb(sum(1)+3*(n-1))
             fion(2,ii) = fion(2,ii) 
     >                  + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                         *dbl_mb(sum(1)+1+3*(n-1))
             fion(3,ii) = fion(3,ii) 
     >                  + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                         *dbl_mb(sum(1)+2+3*(n-1))
          end do
       
        
        end do

      end do

      value = MA_pop_stack(sum(2))
      value = MA_pop_stack(sw1(2))
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

      call nwpw_timing_end(6)
      return 
      end




      subroutine psp_read(fname,
     >                       version,
     >                       nfft,unita,
     >                       atom,amass,zv,lmmax,lmax,locp,
     >                       rc,vnlnrm,rlocal,
     >                       nfft3d,npack1,npack0,
     >                       vl,vnl,
     >                       semicore,rcore,ncore,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*20 fname
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      character*2 atom
      real*8 amass,zv
      integer lmmax
      integer lmax
      integer locp
      real*8 rc(*),vnlnrm(*)
      real*8 rlocal
      integer nfft3d,npack1,npack0
      real*8 vl(*)
      real*8 vnl(npack1,*)
      logical semicore
      real*8  rcore
      real*8  ncore(*)
      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer n,l
      integer msglen
      integer iatom(2)
      character*255 full_filename


      call Parallel_taskid(taskid)
      rlocal = 0.0d0
  
*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call cread(5,atom,2)
         call dread(5,amass,1)
         call dread(5,zv,1)
         call iread(5,lmax,1)
         call iread(5,locp,1)
         lmmax=(lmax+1)**2 - (2*locp+1)
         amass = amass*1822.89d0
         call dread(5,rc,lmax+1)
         call dread(5,vnlnrm,lmmax)
         if (version.eq.4) call dread(5,rlocal,1)
         call dread(5,rcore,1)
      end if


*     **** send header data to all processors ****
#ifdef MPI
      msglen = 1
      call MPI_BCAST(version,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 3
      call MPI_BCAST(nfft,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 9
      call MPI_BCAST(unita,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 2
      call MPI_BCAST(atom,msglen,MPI_CHARACTER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 1
      call MPI_BCAST(amass,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 1
      call MPI_BCAST(zv,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 1
      call MPI_BCAST(lmax,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      call MPI_BCAST(locp,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      lmmax=(lmax+1)**2 - (2*locp+1)

      msglen=lmax+1
      call MPI_BCAST(rc,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen=lmmax
      call MPI_BCAST(vnlnrm,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen=1
      call MPI_BCAST(rlocal,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen=1
      call MPI_BCAST(rcore,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,version,mitob(msglen),MASTER)
      msglen = 3
      call BRDCST(9+MSGINT,nfft,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call BRDCST(9+MSGCHR,iatom,mitob(msglen),MASTER)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call BRDCST(9+MSGDBL,amass,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGDBL,zv,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,lmax,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,locp,mitob(msglen),MASTER)
      lmmax=(lmax+1)**2 - (2*locp+1)

      msglen=lmax+1
      call BRDCST(9+MSGDBL,rc,mdtob(msglen),MASTER)
      msglen=lmmax
      call BRDCST(9+MSGDBL,vnlnrm,mdtob(msglen),MASTER)
      msglen=1
      call BRDCST(9+MSGDBL,rlocal,mdtob(msglen),MASTER)
      msglen=1
      call BRDCST(9+MSGDBL,rcore,mdtob(msglen),MASTER)

#endif

*     **** determine semicore value ****
      if (rcore.gt.0.0d0) then
         semicore = .true.
      else
         semicore = .false.
      end if


 
*     *** read in vl 3d block ***
      call D3dB_t_read(1,5,tmp2,tmp)
      call Pack_t_pack(0,tmp2)
      call Pack_t_Copy(0,tmp2,vl)

*     **** read in vnl 3d blocks ****
      do n=1,lmmax
         call D3dB_t_read(1,5,tmp2,tmp)
         call Pack_t_pack(1,tmp2)
         call Pack_t_Copy(1,tmp2,vnl(1,n))
      end do

*     **** read in semicore density block ****
      if (semicore) then
         call D3dB_t_read(1,5,tmp2,tmp)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1))

         call D3dB_t_read(1,5,tmp2,tmp)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+2*npack0))

         call D3dB_t_read(1,5,tmp2,tmp)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+3*npack0))

         call D3dB_t_read(1,5,tmp2,tmp)
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,ncore(1+4*npack0))
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
*     *									*
*     *	 	  psp_readall  				*
*     *									*
*     ***********************************

      subroutine psp_readall()
      implicit none
      
#include "mafdecls.fh"

*     **** version4 common block ****
      integer rlocal(2)
      common / version4 / rlocal

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** psp common block ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


*     **** local variables ****
      integer ngp(3),version,nfft3d,npack1,npack0
      integer ia,l
      real*8 unita(3,3)
      character*12 boundry
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box
      character*4  element
      character*20 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind
      integer      control_ngrid
      real*8       control_unita
      character*12 control_boundry
      character*2  ion_atom
      external     nwpw_filefind
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
      if (.not. value) call errquit('out of stack memory',0)

      value = MA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** read pseudopotentials ****
      do ia=1,npsp

*      **** define formatted psp name ****
       element = '    '
       element = ion_atom(ia)
       l = index(element,' ') - 1
       fname = element(1:l)//'.vpp'
        

       found = .false.
       do while (.not.found)

         if (nwpw_filefind(fname)) then
            call psp_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  atom(ia),
     >                  dbl_mb(amass(1)+ia-1),
     >                  dbl_mb(zv(1)+ia-1),
     >                  int_mb(lmmax(1)+ia-1),
     >                  int_mb(lmax(1)+ia-1),
     >                  int_mb(locp(1)+ia-1),
     >                  dbl_mb(rc(1) + (ia-1)*(lmax_max+1)),
     >                  dbl_mb(vnlnrm(1) + (ia-1)*lmmax_max),
     >                  dbl_mb(rlocal(1) + (ia-1)),
     >                  nfft3d,npack1,npack0,
     >                  dbl_mb(vl(1) + (ia-1)*npack0),
     >                  dbl_mb(vnl(1)+ (ia-1)*npack1*lmmax_max),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(rcore(1)+ia-1),
     >                  dbl_mb(ncore(1)+ (ia-1)*npack0*5),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)


*          **** set semicore(0) ****
           if (log_mb(semicore(1)+ia)) log_mb(semicore(1)) = .true.
           if (ierr.gt.0) go to 9000

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
              if (taskid.eq.MASTER) then
              write(6,*) "pseudopotential is not correctly formatted:",
     >                    fname
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
 9000 value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp(2))

      return
      end

