
*
* $Id$
*

      subroutine coulomb_init()
      implicit none 
#include "errquit.fh"

#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
c     real*8 vc(nfft3d)
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl

*     **** local variables ****
      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny,nxh,nyh
      real*8  fourpi,gg
      logical value
      integer tmp1(2)

*     **** external functions ****
*     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx
      double precision toll
      parameter (toll=1d-16)

      call nwpw_timing_start(7)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

*     **** allocate vc memory ****
      value = BA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      nxh=nx/2
      nyh=ny/2

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
c     call D3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nxh+1)*ny
c    >     + j*(nxh+1)
c    >     + i+1
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

*     ***** form Vc = 4*pi/G**2  *****
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))
      do i = 1,nfft3d
         
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if (((pzero.eq.taskid) .and. (i.eq.zero)).or.
     E    (abs(gg) .lt.toll))  then
            dbl_mb(tmp1(1)+i-1) = 0.0d0
         else
            dbl_mb(tmp1(1)+i-1) = fourpi/gg
         end if
        
      end do
      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx))
      value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(7)


      return
      end

      subroutine coulomb_end()
      implicit none
#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl
      logical value

      value = BA_free_heap(vc_hndl)
      return
      end


      subroutine coulomb_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "bafdecls.fh"


*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl

      call nwpw_timing_start(7)
      call Pack_tc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
      call nwpw_timing_end(7)

      return 
      end


      real*8 function coulomb_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
*     real*8 vc(nfft3d)
*     common / vc_block / vc
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     **** local variables ****
      integer npack0
      real*8 ec

c     real*8  tmp1(*)
      integer tmp1(2)
      logical value

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(7)
      call Pack_npack(0,npack0)
      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_ttp_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()

      value = BA_pop_stack(tmp1(2))
      call nwpw_timing_end(7)

      coulomb_e = ec
      return 
      end





      subroutine coulomb_euv(dng,euv)
*
* $Id$
*
      implicit none 
#include "errquit.fh"
      complex*16 dng(*)
      real*8 euv(3,3)

#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     **** local variables ****
      integer npack0,nfft3d,G(2,3)
      integer i,j
      integer u,v,s
      logical value

      real*8 pi,fourpi,scal,ss,sum
      real*8 hm(3,3),Bus(3,3),ecoul 
      integer tmp1(2),tmp2(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

      real*8   lattice_unitg,lattice_omega,coulomb_e
      external lattice_unitg,lattice_omega,coulomb_e



      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('out of stack  memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,npack0,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(G(1,1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(G(1,2)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(G(1,3)),1)
      call Pack_t_pack(0,dbl_mb(G(1,1)))
      call Pack_t_pack(0,dbl_mb(G(1,2)))
      call Pack_t_pack(0,dbl_mb(G(1,3)))

*     **** tmp2(G) = (n(G)**2) * (4*pi/G**2)**2  ****
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_Mul(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),
     >                                   dbl_mb(tmp2(1)))
c      call Pack_tt_Mul(0,dbl_mb(tmp2(1)),dbl_mb(vc_indx),
c     >                                   dbl_mb(tmp2(1)))
      call Pack_tt_Mul2(0,dbl_mb(vc_indx),dbl_mb(tmp2(1)))


*     **** Bus = Sum(G) (omega/4*pi)*tmp2(G)*Gu*Gs ****
      call dcopy(9,0.0d0,0,Bus,1)
      ss     = lattice_omega()/fourpi
      do u=1,3
      do s=u,3
        call Pack_tt_Mul(0,dbl_mb(G(1,u)),
     >                     dbl_mb(G(1,s)),
     >                     dbl_mb(tmp1(1)))
        call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)

        Bus(u,s) = ss*sum
      end do
      end do
      do u=1,3
      do s=u+1,3
         Bus(s,u) = Bus(u,s)
      end do
      end do

      ecoul = coulomb_e(dng)
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


*     **********************************************
*     *                                            *
*     *              coulomb_efg                   *
*     *                                            *
*     **********************************************
      subroutine coulomb_efg(dng,efg_smoothr,efg_smoothi)

      implicit none 

      complex*16 dng(*)
      real*8 efg_smoothr(3,3,*),efg_smoothi(3,3,*) ! real and complex parts

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl

*     **** local variables ****
      logical value
      integer npack0,nfft3d
      integer ii,k
      integer u,v
      integer Gtmp(2),G(3),tmp1(2),tmp2(2),exi(2)
      real*8 w,gg
      integer mu,nu,nion
      real*8 rdng,cdng,gvec(3),termgg
      real*8 pi,two_pi,four_pi,sqrt_pi
      real*8 phase,cosgr,singr
      real*8 r1,r2,r3
      complex*16 zsum

*     **** external functions ****
      integer  G_indx,ion_nion
      external G_indx,ion_nion
      real*8   lattice_omega,ion_rion
      external lattice_omega,ion_rion

      pi = 4.0d0*datan(1.0d0)
      sqrt_pi  = dsqrt(pi)
      two_pi   = 2.0d0*pi
      four_pi  = 4.0d0*pi

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      value = BA_push_get(mt_dbl,3*nfft3d,'Gtmp',Gtmp(2),Gtmp(1))
      G(1) = Gtmp(1)
      G(2) = Gtmp(1)+nfft3d
      G(3) = Gtmp(2)+nfft3d
      value = value.and.
     >        BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,npack0,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      if (.not. value) 
     > call errquit('coulomb_efg:out of stack memory',0, MA_ERR)

      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(G(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(G(2)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(G(3)),1)
      call Pack_t_pack(0,dbl_mb(G(1)))
      call Pack_t_pack(0,dbl_mb(G(2)))
      call Pack_t_pack(0,dbl_mb(G(3)))

*     **** initialize ***
      nion = ion_nion()
      do mu = 1,3
        do nu = 1,3
          do ii=1,nion
             efg_smoothr(mu,nu,ii) = 0.0d0
             efg_smoothi(mu,nu,ii) = 0.0d0
          end do ! ii
        end do ! nu
      end do ! mu

      !*** compute (-1/3)gg ***
       do k=1,npack0
          dbl_mb(tmp1(1)+k-1) = -( dbl_mb(G(1)+k-1)**2
     >                           + dbl_mb(G(2)+k-1)**2
     >                           + dbl_mb(G(3)+k-1)**2)/3.0d0
       end do ! k
c
       do ii=1,nion
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         do mu=1,3
            call Pack_ttt_dzaxpy(0,dbl_mb(G(mu)),
     >                             dbl_mb(G(mu)),
     >                             dbl_mb(tmp1(1)),
     >                             dbl_mb(tmp2(1)))
            call Pack_tt_Mul2(0,dbl_mb(vc_indx),dbl_mb(tmp2(1)))

            zsum = dcmplx(0.0d0,0.0d0) 
            do k=1,npack0
               zsum = zsum 
     >              - 2.0d0*dng(k)
     >                     *dcpl_mb(exi(1)+k-1)
     >                     *dbl_mb(tmp2(1)+k-1)
            end do
            efg_smoothr(mu,mu,ii) = dble(zsum)
            efg_smoothi(mu,mu,ii) = dimag(zsum)

            do nu=mu+1,3
               call Pack_tt_Mul(0,dbl_mb(G(mu)),
     >                            dbl_mb(G(nu)),
     >                            dbl_mb(tmp2(1)))
               call Pack_tt_Mul2(0,dbl_mb(vc_indx),dbl_mb(tmp2(1)))

               zsum = dcmplx(0.0d0,0.0d0) 
               do k=1,npack0
                  zsum = zsum - 2.0d0*dng(k)
     >                               *dcpl_mb(exi(1)+k-1)
     >                               *dbl_mb(tmp2(1)+k-1)
               end do
               efg_smoothr(mu,nu,ii) = dble(zsum)
               efg_smoothr(nu,mu,ii) = dble(zsum)
               efg_smoothi(mu,nu,ii) = dimag(zsum)
               efg_smoothi(nu,mu,ii) = dimag(zsum)

            end do  ! nu
         end do ! mu
      end do ! ii

c      do k=1,npack0
c        gvec(1) = dbl_mb(G(1)+k-1)
c        gvec(2) = dbl_mb(G(2)+k-1)
c        gvec(3) = dbl_mb(G(3)+k-1)
c        gg = gvec(1)*gvec(1) + gvec(2)*gvec(2) + gvec(3)*gvec(3)
c        if (gg .gt. 0.d0) then  ! eliminate g=o term
c          do ii=1,nion
c            r1 = ion_rion(1,ii)
c            r2 = ion_rion(2,ii)
c            r3 = ion_rion(3,ii)
c            phase = (gvec(1)*r1+gvec(2)*r2+gvec(3)*r3)
c            cosgr = cos(phase)  ! cos(G.R)
c            singr = sin(phase)  ! sin(G.R)
c            do mu=1,3
c             do nu=1,3
c                termgg = gvec(mu)*gvec(nu) 
c                if (mu == nu) termgg = termgg - gg/3.d0
c                rdng = real(dng(k))
c                cdng = aimag(dng(k))
c!
c                efg_smoothr(mu,nu,ii) = efg_smoothr(mu,nu,ii) -   ! real part
c     &            2.0d0*four_pi*termgg*(rdng*cosgr-cdng*singr)/gg 
c                efg_smoothi(mu,nu,ii) = efg_smoothi(mu,nu,ii) -   ! imag part
c     &            2.0d0*four_pi*termgg*(rdng*singr+cdng*cosgr)/gg 
c             end do ! nu
c            end do ! mu
c          end do ! ii
c        end if ! gg ne 0
c      end do ! k

      call D3dB_Vector_SumAll(9*nion,efg_smoothr)
      call D3dB_Vector_SumAll(9*nion,efg_smoothi)
!
      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and. BA_pop_stack(Gtmp(2))
      if (.not. value) call errquit('coulomb_efg:popstack',1,MA_ERR)
!
      return
      end
