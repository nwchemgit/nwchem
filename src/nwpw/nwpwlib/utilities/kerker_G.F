*
* $Id: kerker_G.F,v 1.2 2009-02-08 03:26:14 bylaska Exp $
*

      subroutine kerker_G_init(g0)
      implicit none 
      real*8 g0

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer npack0,nfft3d,G(3)
      integer i
      real*8  gg,gg0
      logical value
 
      integer tmp1(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

*     **** common block used for kerker_G.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / kerker_g_block / tg_indx,tg_hndl


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack0)
      G(1)= G_indx(1)
      G(2)= G_indx(2)
      G(3)= G_indx(3)

      value = MA_alloc_get(mt_dbl,npack0,
     >                     'tg',tg_hndl,tg_indx)
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      gg0 = g0*g0
      do i = 1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))
         dbl_mb(tmp1(1)+i-1) = gg/(gg+gg0)
      end do

      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(tg_indx))
      

      value = MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      return
      end

      subroutine kerker_G_end()
      implicit none 

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block used for kerker_G.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / kerker_G_block / tg_indx,tg_hndl

      logical value

      value = MA_free_heap(tg_hndl)
      if (.not. value) call errquit('error freeing heap',0, MA_ERR)
      return
      end


      subroutine kerker_G(ispin,V)
      implicit none
      integer    ispin
      real*8     V(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block used for kerker_G.F ****
      integer tg_indx,tg_hndl
      common / kerker_g_block / tg_indx,tg_hndl

*     **** local variables ****
      logical value
      integer nfft3d,n2ft3d,tmp1(2)
      integer ms,nx,ny,nz
      real*8  scal1


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

      
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d


      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      do ms=1,ispin
         call dcopy(n2ft3d,V((ms-1)*n2ft3d+1),1,dcpl_mb(tmp1(1)),1)
         call D3dB_r_Zero_Ends(1,dcpl_mb(tmp1(1)))
c         call D3dB_r_SMul(1,scal1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
         call D3dB_r_SMul1(1,scal1,dcpl_mb(tmp1(1)))
         call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))

c         call Pack_tc_Mul(0,dbl_mb(tg_indx),
c     >                      dcpl_mb(tmp1(1)),
c     >                      dcpl_mb(tmp1(1)))
         call Pack_tc_Mul2(0,dbl_mb(tg_indx),
     >                      dcpl_mb(tmp1(1)))

         call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
         call D3dB_r_Zero_Ends(1,dcpl_mb(tmp1(1)))
         call dcopy(n2ft3d,dcpl_mb(tmp1(1)),1,V((ms-1)*n2ft3d+1),1)
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      
      return 
      end

