*
* $Id$
*
      subroutine c_coulomb_init()
      implicit none 

#include "bafdecls.fh"
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
      value = BA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) 
     > call errquit('c_coulomb_init:out of heap memory',0,0)

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('c_coulomb_init:out of stack memory',0,0)

      call Parallel3d_taskid_i(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
c     call C3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nx)*ny
c    >     + j*(nx)
c    >     + i+1
      call C3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)


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
      value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(7)


      return
      end

      subroutine c_coulomb_end()
      implicit none
#include "bafdecls.fh"
#include "c_coulomb_common.fh"

      logical value

      value = BA_free_heap(vc_hndl)
      return
      end


      subroutine c_coulomb_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "bafdecls.fh"
#include "c_coulomb_common.fh"

      call nwpw_timing_start(7)
     
      call Cram_rc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
      
      call nwpw_timing_end(7)

      return 
      end



      real*8 function c_coulomb_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
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
      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('c_coulomb_e:out of stack memory',0,0)
      
      call Cram_cr_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()


      value = BA_pop_stack(tmp1(2))
      
      call nwpw_timing_end(7)

      c_coulomb_e = ec
      return 
      end




      subroutine c_coulomb_euv(dng,euv)
      implicit none 
      complex*16 dng(*)
      real*8 euv(3,3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_coulomb_common.fh"


*     **** local variables ****
      integer npack0,nfft3d,G(2,3)
      integer i,j
      integer u,v,s
      logical value

      real*8 pi,fourpi,scal,ss,sum
      real*8 hm(3,3),Bus(3,3),ecoul 
      integer tmp1(2),tmp2(2)

*     **** external functions ****
      integer  c_G_indx
      external c_G_indx

      real*8   lattice_unitg,lattice_omega,c_coulomb_e
      external lattice_unitg,lattice_omega,c_coulomb_e


      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      


      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)

      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) 
     > call errquit('c_coulomb_euv:out of stack memory',0,MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) 
     > call errquit('c_coulomb_euv:out of stack memory',0,MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) 
     > call errquit('c_coulomb_euv:out of stack  memory',0,MA_ERR)

      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('c_coulomb_euv:out of stack memory',0,MA_ERR)

      value = BA_push_get(mt_dbl,npack0,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('c_coulomb_euv:out of stack memory',0,MA_ERR)

      call dcopy(nfft3d,dbl_mb(c_G_indx(1)),1,dbl_mb(G(1,1)),1)
      call dcopy(nfft3d,dbl_mb(c_G_indx(2)),1,dbl_mb(G(1,2)),1)
      call dcopy(nfft3d,dbl_mb(c_G_indx(3)),1,dbl_mb(G(1,3)),1)
      call Cram_r_pack(0,dbl_mb(G(1,1)))
      call Cram_r_pack(0,dbl_mb(G(1,2)))
      call Cram_r_pack(0,dbl_mb(G(1,3)))

*     **** tmp2(G) = (n(G)**2) * (4*pi/G**2)**2  ****
      call Cram_cr_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Cram_rr_Mul(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),
     >                                   dbl_mb(tmp2(1)))
c      call Cram_rr_Mul(0,dbl_mb(tmp2(1)),dbl_mb(vc_indx),
c     >                                   dbl_mb(tmp2(1)))
      call Cram_rr_Mul2(0,dbl_mb(vc_indx),dbl_mb(tmp2(1)))


*     **** Bus = Sum(G) (omega/4*pi)*tmp2(G)*Gu*Gs ****
      call dcopy(9,0.0d0,0,Bus,1)
      ss     = lattice_omega()/fourpi
      do u=1,3
      do s=u,3
        call Cram_rr_Mul(0,dbl_mb(G(1,u)),
     >                     dbl_mb(G(1,s)),
     >                     dbl_mb(tmp1(1)))
        call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)

        Bus(u,s) = ss*sum
      end do
      end do
      do u=1,3
      do s=u+1,3
         Bus(s,u) = Bus(u,s)
      end do
      end do

      ecoul = c_coulomb_e(dng)
      do v=1,3
      do u=1,3
         euv(u,v) = -ecoul*hm(u,v)
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do
 
      value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp1(2))
      value = BA_pop_stack(G(2,3)) 
      value = BA_pop_stack(G(2,2)) 
      value = BA_pop_stack(G(2,1))

      return
      end
