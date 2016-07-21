*
* $Id: psi_dmatrix_localize.F 26429 2014-12-03 21:31:11Z bylaska $
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
*             neq  : number of wavefunctions
*             psi_r: wavefunctions
*
*     Exit - 
*             psi_r2: localized wavefunctions
*

      subroutine psi_dmatrix_localize(ispin,neq,n2ft3d,psi_r,psi_r2)
      implicit none
      integer    ispin,neq(2),n2ft3d
      real*8     psi_r(n2ft3d,neq(1)+neq(2))
      real*8     psi_r2(n2ft3d,neq(1)+neq(2))

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid_i,taskid_j
      real*8 x,y,z
      integer ms,i,j,ne(2),n1(2),n2(2),neall,nn,ishift
      integer i0,k0,j0,ii,jj,taskid_i0,taskid_j0
      integer dm_indexs(2)


      ne(1) = neq(1)
      ne(2) = neq(2)
      call D1dB_Vector_ISumAll(2,ne)

*     **** allocate memory ****
      nn = 2*(ne(1)+ne(2))*ispin
      if (.not.BA_push_get(mt_int,nn,'dm_indexs',
     >                     dm_indexs(2),dm_indexs(1))) then
         call errquit('psi_dmatrix_localize:push stack',0,MA_ERR)
      end if

      call psi_dmatrix_maxallcolumns(ispin,ne,neq,n2ft3d,
     >                                     psi_r,psi_r2,
     >                                     int_mb(dm_indexs(1)))
      call psi_dmatrix_rotate(ispin,ne,neq,n2ft3d,int_mb(dm_indexs(1)),
     >                              psi_r,psi_r2)

*     **** dealloate memory ****
      if (.not.BA_pop_stack(dm_indexs(2))) 
     >   call errquit('psi_dmatrix_localize:pop stack',2,MA_ERR)
     
      return
      end


*     *********************************
*     *                               *
*     *     psi_dmatrix_rotate        *
*     *                               *
*     *********************************
*
*   This routine give the location of the density matrix columns rotates psi_r to psi_r2
*
*
*     Entry - ispin: 
*             ne   : number of wavefunctions
*             neq  : distributed number of wavefunctions
*             n2ft3d : size of first dimension of psi_r and psi_r2
*             dm_indexes: location of density matrix columns
*             psi_r: wavefunctions
*
*     Exit - 
*             psi_r2: localized wavefunctions
*

      subroutine psi_dmatrix_rotate(ispin,ne,neq,n2ft3d,dm_indexes,
     >                              psi_r,psi_r2)
      implicit none
      integer    ispin,ne(2),neq(2),n2ft3d
      integer    dm_indexes(2,ispin,*)
      real*8     psi_r(n2ft3d,neq(1)+neq(2))
      real*8     psi_r2(n2ft3d,neq(1)+neq(2))

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid_i,taskid_j
      integer ms,i,j,n1(2),n2(2),ishift
      integer i0,k0,j0,ii,jj,taskid_i0,taskid_j0
      integer summ(2),msize,nfft3d

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      call nwpw_timing_start(54)

      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      call Dneall_m_size(0,msize)

*     **** allocate tempory matrix ****
      if (.not.Dneall_m_push_get(0,summ))
     >   call errquit('psi_dmatrix_rotate:push stack',2,MA_ERR)


*     **** generate summ matrix ****
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      call Dneall_m_zero(0,dbl_mb(summ(1)))
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
            i = ii-n1(ms)+1
            taskid_i0 = dm_indexes(1,ms,i)
            i0        = dm_indexes(2,ms,i)
            if (taskid_i0.eq.taskid_i) then
               do jj=n1(ms),n2(ms)
                  j = jj-n1(ms)+1
                  call Dneall_ntoqp(jj,j0,taskid_j0)
                  if (taskid_j0.eq.taskid_j) then
                     call Dneall_m_add_value(psi_r(i0,j0),0,
     >                                       ms,j,i,dbl_mb(summ(1)))
                  end if
               end do
            end if
         end do
      end do
      call Parallel_Vector_SumAll(msize,dbl_mb(summ(1)))

*     **** generate psi_r2 from columns of density matrix ****
      call Dneall_gmg_Multiply(0,psi_r,n2ft3d,
     >                         dbl_mb(summ(1)),1.0d0,
     >                         psi_r2,0.0d0)
c         
*     **** inverse dmatrix patch from columns of density matrix ****
      call Dneall_m_zero(0,dbl_mb(summ(1)))
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
            i = ii-n1(ms)+1
            taskid_i0 = dm_indexes(1,ms,i)
            i0        = dm_indexes(2,ms,i)
            if (taskid_i0.eq.taskid_i) then
               do jj=n1(ms),n2(ms)
                  j = jj-n1(ms)+1
                  call Dneall_ntoqp(jj,j0,taskid_j0)
                  if (taskid_j0.eq.taskid_j) then
                     call Dneall_m_add_value(psi_r2(i0,j0),0,
     >                                       ms,j,i,dbl_mb(summ(1)))
                  end if
               end do
            end if
         end do
      end do
      call Parallel_Vector_SumAll(msize,dbl_mb(summ(1)))

*     **** generate localized psi_r2 from cholesky ****
      call Dneall_m_cholesky(0,dbl_mb(summ(1)))
      call Dneall_mg_forwardsolve(0,dbl_mb(summ(1)),n2ft3d,psi_r2)

*     **** deallocate tempory matrix ****
      if (.not.Dneall_m_pop_stack(summ))
     >   call errquit('psi_dmatrix_rotate:pop stack',0,MA_ERR)

      call nwpw_timing_end(54)

      return
      end





*     *********************************
*     *                               *
*     *   psi_dmatrix_maxallcolumns   *
*     *                               *
*     *********************************
*
*   This routine calculates the indexes using a simple rank deficient qr
*
*
*     Entry - ispin: 
*             neq  : number of wavefunctions
*             psi_r: wavefunctions
*             psi_r2: tempory data
*
*     Exit - 
*             dm_indexes: location of density matrix columns
*             psi_r2: tempory data destroyed
*

      subroutine psi_dmatrix_maxallcolumns(ispin,ne,neq,n2ft3d,
     >                                     psi_r,psi_r2,
     >                                     dm_indexes)
      implicit none
      integer    ispin,ne(2),neq(2),n2ft3d
      real*8     psi_r(n2ft3d,neq(1)+neq(2))
      real*8     psi_r2(n2ft3d,neq(1)+neq(2))
      integer    dm_indexes(2,ispin,*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer ms,i,neall,ishift
      integer k0,taskid_i0

      call nwpw_timing_start(53)

      neall = neq(1) + neq(2)
      call dcopy(n2ft3d*neall,psi_r,1,psi_r2,1)

*     **** find the maximum columns of density matrix ****
      do ms=1,ispin
         ishift = 1 + (ms-1)*(neq(1))
         do i=1,ne(ms)
            call psi_dmatrix_maxcolumn(neq(ms),n2ft3d,psi_r2(1,ishift),
     >                        taskid_i0,k0)
            call psi_dmatrix_projectout(neq(ms),n2ft3d,psi_r2(1,ishift),
     >                                  taskid_i0,k0)
            dm_indexes(1,ms,i) = taskid_i0
            dm_indexes(2,ms,i) = k0
         end do
      end do

      call nwpw_timing_end(53)

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
      integer np_i,taskid_i,np_j
      integer i,j,k,kk
      integer maxsum0(2)
      real*8  tsum,maxsum

      call Parallel2d_np_i(np_i) 
      call Parallel2d_np_j(np_j) 
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
         if (np_j.gt.1) call D1dB_SumAll(tsum)
         if (tsum.gt.maxsum) then
            maxsum = tsum
            kk = k
         end if
      end do
      dbl_mb(maxsum0(1)+taskid_i) = maxsum
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
      integer taskid_i,i,k,np_j
      integer v0(2)

*     **** external functions ****
      real*8   ddot
      external ddot

      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_i(taskid_i)

      if (.not.BA_push_get(mt_dbl,nn,'v0',v0(2),v0(1)))
     > call errquit('psi_dmaxtrix_projectout:push stack',0,MA_ERR)
      call dcopy(nn,0.0d0,0,dbl_mb(v0(1)),1)
      

c     **** collect the vector across tasks ****
      if (taskid_i.eq.taskid_i0) then
!$OMP DO private(i)
         do i=1,nn
            dbl_mb(v0(1)+i-1) = psi_r(k0,i)
         end do
!$OMP END DO
      end if
      call D3dB_Vector_SumAll(nn,dbl_mb(v0(1)))
      vv = 0.0d0
      do i=1,nn
         vv = vv + dbl_mb(v0(1)+i-1)**2
      end do
!$OMP MASTER
      !vv = ddot(nn,dbl_mb(v0(1)),1,dbl_mb(v0(1)),1)
      if (np_j.gt.1) call D1dB_SumAll(vv)
      if (vv.gt.1.0d-11) then
         ss = 1.0d0/vv
      else
         ss = 0.0d0
      end if
!$OMP END MASTER
!$OMP BARRIER

c     **** project out vector ****
!$OMP DO private(k)
      do k=1,n2ft3d
         uv = 0.0d0
         do i=1,nn
            uv = uv + psi_r(k,i)*dbl_mb(v0(1)+i-1)
         end do
         !uv = ddot(nn,psi_r(k,1),n2ft3d,dbl_mb(v0(1)),1)

         if (np_j.gt.1) call D1dB_SumAll(uv)
         do i=1,nn
            psi_r(k,i) = psi_r(k,i) - ss*uv*dbl_mb(v0(1)+i-1)
         end do
         !call daxpy(nn,(-ss*uv),dbl_mb(v0(1)),1,psi_r(k,1),n2ft3d)
      end do
!$OMP END DO
      !do i=1,nn
      !   call D3dB_r_Zero_Ends(1,psi_r(1,i))
      !end do

      if (.not.BA_pop_stack(v0(2)))
     > call errquit('psi_dmaxtrix_projectout:pop stack',0,MA_ERR)

      return
      end 

   

