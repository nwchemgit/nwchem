*
* $Id: coulomb_screened.F,v 1.2 2002-01-11 17:55:33 bylaska Exp $
*
      subroutine coulomb_screened_init(rcut)
      implicit none 
      real*8 rcut

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     **** local variables ****
      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny,nxh,nyh
      real*8  fourpi,gg,rr,f
      logical value
      integer tmp1(2)


*     **** external functions ****
*     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

      call nwpw_timing_start(7)

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

*     **** allocate vc memory ****
      value = MA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) call errquit('out of heap memory',0)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)

      call Parallel_taskid(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      nxh=nx/2
      nyh=ny/2


*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call D3dB_ktoqp(1,k+1,qzero,pzero)
      zero = (qzero-1)*(nxh+1)*ny
     >     + j*(nxh+1)
     >     + i+1

*     ***** form Vc = 4*pi/G**2  *****
      rr     = rcut*rcut
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))
      do i = 1,nfft3d
         
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         f = 1.0d0 - dexp(-0.25d0*rr*gg)

         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            dbl_mb(tmp1(1)+i-1) = 0.0d0
         else
            dbl_mb(tmp1(1)+i-1) = -f*fourpi/gg
         end if
        
      end do
      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx))
      value = MA_pop_stack(tmp1(2))

      call nwpw_timing_end(7)


      return
      end

      subroutine coulomb_screened_end()
      implicit none
#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl
      logical value

      value = MA_free_heap(vc_hndl)
      return
      end


      subroutine coulomb_screened_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "mafdecls.fh"


*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

      call nwpw_timing_start(7)
      
      call Pack_tc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
     
      call nwpw_timing_end(7)

      return 
      end



      real*8 function coulomb_screened_e(dng)
      implicit none
      complex*16 dng(*)

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl


*     **** local variables ****
      integer npack0
      real*8 ec

c     real*8  tmp1(*)
      integer tmp1(2)
      logical value

*     **** external functions ****
      real*8 lattice_omega
      external lattice_omega

      
      call nwpw_timing_start(7)

      call Pack_npack(0,npack0)
      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)
      
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()

      value = MA_pop_stack(tmp1(2))
      
      call nwpw_timing_end(7)

      coulomb_screened_e = ec
      return 
      end

