*
* $Id$
*

      subroutine kerker_G_init()
      implicit none 

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer npack0,nfft3d,G(3)
      integer i,ccode
      real*8  gg,gg0,g0
      integer tmp1(2)

*     **** external functions ****
      real*8   control_kerker_g0
      external control_kerker_g0
      integer  G_indx,c_G_indx,control_code
      external G_indx,c_G_indx,control_code

*     **** common block used for kerker_G.F ****
c     real*8 tg(nfft3d)
      logical isband,dokerker
      integer tg_indx,tg_hndl
      common / kerker_g_block / tg_indx,tg_hndl,isband,dokerker

   
      g0 = control_kerker_g0()

!$OMP MASTER
      dokerker =  (g0.gt.0.0d0)
!$OMP END MASTER
!$OMP BARRIER

      if (dokerker) then
      
      ccode = control_code()
!$OMP MASTER
      isband =((ccode.eq.5).or.(ccode.eq.13).or.(ccode.eq.14))
!$OMP END MASTER
!$OMP BARRIER

      if (isband) then
         call C3dB_nfft3d(1,nfft3d)
         call Cram_npack(0,npack0)
         G(1)= c_G_indx(1)
         G(2)= c_G_indx(2)
         G(3)= c_G_indx(3)
      else
         call D3dB_nfft3d(1,nfft3d)
         call Pack_npack(0,npack0)
         G(1)= G_indx(1)
         G(2)= G_indx(2)
         G(3)= G_indx(3)
      end if

      if (.not.BA_alloc_get(mt_dbl,npack0,
     >                     'tg',tg_hndl,tg_indx))
     >   call errquit('kerker_G_init:out of heap memory',0,MA_ERR)

      if (.not.BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit('kerker_G_init:out of stack memory',0,MA_ERR)


      gg0 = g0*g0
!$OMP DO
      do i = 1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))
         dbl_mb(tmp1(1)+i-1) = gg/(gg+gg0)
      end do
!$OMP END DO

      if (isband) then
         call Cram_r_pack(0,dbl_mb(tmp1(1)))
         call Cram_r_Copy(0,dbl_mb(tmp1(1)),dbl_mb(tg_indx))
      else
         call Pack_t_pack(0,dbl_mb(tmp1(1)))
         call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(tg_indx))
      end if
      if (.not.BA_pop_stack(tmp1(2)))
     >   call errquit('kerker_G_init:popping stack memory',0,MA_ERR)

      end if
      return
      end

      subroutine kerker_G_end()
      implicit none 

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for kerker_G.F ****
      logical isband,dokerker
      integer tg_indx,tg_hndl
      common / kerker_G_block / tg_indx,tg_hndl,isband,dokerker

      if (dokerker) then
         if (.not.BA_free_heap(tg_hndl))
     >   call errquit('error freeing heap',0, MA_ERR)
      end if
      return
      end


      subroutine kerker_G(V)
      implicit none
      real*8     V(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for kerker_G.F ****
      logical isband,dokerker
      integer tg_indx,tg_hndl
      common / kerker_g_block / tg_indx,tg_hndl,isband,dokerker

*     **** local variables ****
      integer nfft3d,n2ft3d,tmp1(2)
      integer nx,ny,nz
      real*8  scal1

      if (dokerker) then

      if (isband) then
         call C3dB_nx(1,nx)
         call C3dB_ny(1,ny)
         call C3dB_nz(1,nz)
         call C3dB_nfft3d(1,nfft3d)
         n2ft3d = nfft3d
      else
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         call D3dB_nfft3d(1,nfft3d)
         n2ft3d = 2*nfft3d
      end if
      scal1 = 1.0d0/dble(nx*ny*nz)


      if (.not.BA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit('kerker_G: out of stack memory',0, MA_ERR)
      if (isband) then
         call C3dB_rc_SMul(1,scal1,V,dcpl_mb(tmp1(1)))
         call C3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
         call Cram_c_pack(0,dcpl_mb(tmp1(1)))
         call Cram_rc_Mul2(0,dbl_mb(tg_indx),dcpl_mb(tmp1(1)))
         call Cram_c_unpack(0,dcpl_mb(tmp1(1)))
         call C3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
         call C3dB_cr_real(1,dcpl_mb(tmp1(1)),V)
      else
         call D3dB_r_SMul(1,scal1,V,dcpl_mb(tmp1(1)))
         call D3dB_r_Zero_Ends(1,dcpl_mb(tmp1(1)))
         call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
         call Pack_c_pack(0,dcpl_mb(tmp1(1)))
         call Pack_tc_Mul2(0,dbl_mb(tg_indx),dcpl_mb(tmp1(1)))
         call Pack_c_unpack(0,dcpl_mb(tmp1(1)))
         call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
         call D3dB_r_Zero_Ends(1,dcpl_mb(tmp1(1)))
         !call dcopy(n2ft3d,dcpl_mb(tmp1(1)),1,V,1)
         call Parallel_shared_vector_copy(.true.,n2ft3d,
     >                dcpl_mb(tmp1(1)),V)
      end if
      if (.not.BA_pop_stack(tmp1(2)))
     >   call errquit('kerker_G: popping stack memory',0,MA_ERR)


      endif
      return 
      end

