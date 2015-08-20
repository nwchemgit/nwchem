*
* $Id: silvestrelli_overlap.F 26429 2014-12-03 21:31:11Z bylaska $
*


*     *********************************
*     *                               *
*     *     psi_dmatrix_localize      *
*     *                               *
*     *********************************
*
*   This routine calculates the indexes using a simple rank deficient qr
*
*
*     Entry - ispin: 
*             ne  : number of wavefunctions
*             psi_r: wavefunctions
*             psi_r2: tmp wavefunctions
*
*     Exit - 
*

      subroutine psi_dmatrix_localize(ispin,ne,neq,n2ft3d,psi_r,psi_r2)
      implicit none
      integer    ispin,ne(2),neq(2),n2ft3d
      real*8     psi_r(n2ft3d,neq(1)+neq(2))
      real*8     psi_r2(n2ft3d,neq(1)+neq(2))

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      real*8 x,y,z
      integer ms,i,neall,nn,ishift,taskid_i0,k0
      integer dm_indexs(2),rgrid(2)

*     **** external functions ****

      nn = 2*(ne(1)+ne(2))*ispin
      if (.not.BA_push_get(mt_int,nn,'dm_indexs',
     >                     dm_indexs(2),dm_indexs(1))) then
         call errquit('psi_dmatrix_localize:push stack',0,MA_ERR)
      end if
      if (.not.BA_push_get(mt_dbl,3*n2ft3d,'rgrid',rgrid(2),rgrid(1)))
     >   call errquit('psi_dmatrix_localize:push stack',1,MA_ERR)
      call lattice_r_grid(dbl_mb(rgrid(1)))

      neall = neq(1) + neq(2)
      call dcopy(n2ft3d*neall,psi_r,1,psi_r2,1)

      do ms=1,ispin
         ishift = 1 + (ms-1)*(neq(1))
         do i=1,ne(ms)
            call psi_dmatrix_maxcolumn(neq(ms),n2ft3d,psi_r2(1,ishift),
     >                        taskid_i0,k0)
            call psi_dmatrix_projectout(neq(ms),n2ft3d,psi_r2(1,ishift),
     >                                  taskid_i0,k0)
            nn = (i-1)*2 + (ms-1)*2*(neq(1)+neq(2))
            int_mb(dm_indexs(1)+nn)   = taskid_i0
            int_mb(dm_indexs(1)+nn+1) = k0
            write(*,*) "ms,n,taskid_i,k=",ms,i,taskid_i0,k0
            x = dbl_mb(rgrid(1)+3*(k0-1))
            y = dbl_mb(rgrid(1)+3*(k0-1)+1)
            z = dbl_mb(rgrid(1)+3*(k0-1)+2)
            write(*,'(A,3F10.3)') "x,y,z=",x,y,z
            write(*,*)
         end do
      end do

*     **** generate psi_r2 from columns of density matrix ****
*     **** compute overlap of psi_r2 ****
     

      if (.not.BA_pop_stack(rgrid(2)))
     >   call errquit('psi_dmatrix_localize:pop stack',0,MA_ERR)

      if (.not.BA_pop_stack(dm_indexs(2))) then
         call errquit('psi_dmatrix_localize:pop stack',0,MA_ERR)
      end if
      return
      end

*     ****************************************
*     *                                      *
*     *          psi_dmatrix_maxcolumn       *
*     *                                      *
*     ****************************************
      subroutine psi_dmatrix_maxcolumn(nn,n2ft3d,psi_r,taskid_i0,k0)
      implicit none
      integer nn,n2ft3d
      real*8  psi_r(n2ft3d,nn)
      integer taskid_i0,k0

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer np_i,taskid_i
      integer i,j,k,kk
      integer maxsum0(2)
      real*8  tsum,maxsum

      call Parallel2d_np_i(np_i) 
      call Parallel2d_taskid_i(taskid_i) 

      if (.not.BA_push_get(mt_dbl,np_i,'maxsum0',maxsum0(2),maxsum0(1)))
     > call errquit('psi_dmatrix_maxcolumn:push stack',0,MA_ERR)
      call dcopy(np_i,0.0d0,0,dbl_mb(maxsum0(1)),1)

      kk = 1
      maxsum = -9.9d99
      do k=1,n2ft3d
         tsum = 0.0d0
         do j=1,nn
            tsum = tsum +  psi_r(k,j)**2
         end do
         call D1dB_SumAll(tsum)
         if (tsum.gt.maxsum) then
            maxsum = tsum
            kk = k
         end if
      end do
      dbl_mb(maxsum0(1)+taskid_i) = maxsum
      write(*,*) "maxsum=",maxsum,kk,dbl_mb(maxsum0(1)+taskid_i)

      call D3dB_Vector_SumAll(np_i,dbl_mb(maxsum0(1)))

      tsum = -9.9d99
      do i = 1,np_i
        if (dbl_mb(maxsum0(1)+i-1).gt.tsum) then
           tsum = dbl_mb(maxsum0(1)+i-1)
           taskid_i0 = i-1
        end if
      end do

      if (taskid_i0.eq.taskid_i) then
         k0 = kk
      else
         k0 = 0
      end if
      call D3dB_ISumAll(k0)


      if (.not.BA_pop_stack(maxsum0(2)))
     > call errquit('psi_dmatrix_maxcolumn:pop stack',0,MA_ERR)
      return
      end 


*     ****************************************
*     *                                      *
*     *          psi_dmatrix_projectout      *
*     *                                      *
*     ****************************************
      subroutine psi_dmatrix_projectout(nn,n2ft3d,psi_r,
     >                                  taskid_i0,k0)
      implicit none
      integer nn,n2ft3d
      real*8  psi_r(n2ft3d,nn)
      integer taskid_i0,k0

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      real*8  uv,vv,ss
      integer taskid_i,i,k
      integer v0(2)
      

      call Parallel2d_taskid_i(taskid_i)

      if (.not.BA_push_get(mt_dbl,nn,'v0',v0(2),v0(1)))
     > call errquit('psi_dmaxtrix_projectout:push stack',0,MA_ERR)
      call dcopy(nn,0.0d0,0,dbl_mb(v0(1)),1)
      

c     **** collect the vector across tasks ****
      if (taskid_i.eq.taskid_i0) then
         do i=1,nn
            dbl_mb(v0(1)+i-1) = psi_r(k0,i)
         end do
      end if
      call D3dB_Vector_SumAll(nn,dbl_mb(v0(1)))
      vv = 0.0d0
      do i=1,nn
         vv = vv + dbl_mb(v0(1)+i-1)**2
      end do
      call D1dB_SumAll(vv)
      if (vv.gt.1.0d-11) then
         ss = 1.0d0/vv
      else
         ss = 0.0d0
      end if

c     **** project out vector ****
      do k=1,n2ft3d
         uv = 0.0d0
         do i=1,nn
            uv = uv + psi_r(k,i)*dbl_mb(v0(1)+i-1)
         end do
         call D1dB_SumAll(uv)
         do i=1,nn
            psi_r(k,i) = psi_r(k,i) - ss*uv*dbl_mb(v0(1)+i-1)
         end do
      end do

      if (.not.BA_pop_stack(v0(2)))
     > call errquit('psi_dmaxtrix_projectout:pop stack',0,MA_ERR)

      return
      end 

   

