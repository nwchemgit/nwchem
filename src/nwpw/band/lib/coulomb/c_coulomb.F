*
* $Id: c_coulomb.F,v 1.3 2003-10-24 18:56:58 marat Exp $
*
      subroutine c_coulomb_init()
      implicit none 

#include "mafdecls.fh"
#include "c_coulomb_common.fh"

*     **** local variables ****
      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny
      real*8  fourpi,gg
      logical value
      integer tmp1(2)


*     **** external functions ****
*     real*8 G(nfft3d,3)
      integer  c_G_indx
      external c_G_indx

      call nwpw_timing_start(7)
      
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      G(1) = c_G_indx(1)
      G(2) = c_G_indx(2)
      G(3) = c_G_indx(3)

*     **** allocate vc memory ****
      value = MA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) call errquit('out of heap memory',0,0)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0,0)

      call Parallel_taskid(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call C3dB_ktoqp(1,k+1,qzero,pzero)
      zero = (qzero-1)*(nx)*ny
     >     + j*(nx)
     >     + i+1


*     ***** form Vc = 4*pi/G**2  *****
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))
      do i = 1,nfft3d
         
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            dbl_mb(tmp1(1)+i-1) = 0.0d0
         else
            dbl_mb(tmp1(1)+i-1) = fourpi/gg
         end if
        
      end do
      call Cram_r_pack(0,dbl_mb(tmp1(1)))
      call Cram_r_Copy(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx))
      value = MA_pop_stack(tmp1(2))

      call nwpw_timing_end(7)


      return
      end

      subroutine c_coulomb_end()
      implicit none
#include "mafdecls.fh"
#include "c_coulomb_common.fh"

      logical value

      value = MA_free_heap(vc_hndl)
      return
      end


      subroutine c_coulomb_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "mafdecls.fh"
#include "c_coulomb_common.fh"

      call nwpw_timing_start(7)
     
      call Cram_rc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
      
      call nwpw_timing_end(7)

      return 
      end



      real*8 function c_coulomb_e(dng)
      implicit none
      complex*16 dng(*)

#include "mafdecls.fh"
#include "c_coulomb_common.fh"

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
     
      call Cram_npack(0,npack0)
      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0,0)
      
      call Cram_cr_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()


      value = MA_pop_stack(tmp1(2))
      
      call nwpw_timing_end(7)

      c_coulomb_e = ec
      return 
      end




