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
      integer dm_indexs(2),rgrid(2)
      integer summ(2),msize,nfft3d

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      ne(1) = neq(1)
      ne(2) = neq(2)
      call D1dB_Vector_ISumAll(2,ne)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      call Dneall_m_size(0,msize)

      nn = 2*(ne(1)+ne(2))*ispin
      if (.not.BA_push_get(mt_int,nn,'dm_indexs',
     >                     dm_indexs(2),dm_indexs(1))) then
         call errquit('psi_dmatrix_localize:push stack',0,MA_ERR)
      end if
      if (.not.BA_push_get(mt_dbl,3*n2ft3d,'rgrid',rgrid(2),rgrid(1)))
     >   call errquit('psi_dmatrix_localize:push stack',1,MA_ERR)
      call lattice_r_grid(dbl_mb(rgrid(1)))

      if (.not.Dneall_m_push_get(0,summ))
     >   call errquit('psi_dmatrix_localize:push stack',2,MA_ERR)

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
            nn = (i-1)*2 + (ms-1)*2*(ne(1)+ne(2))
            int_mb(dm_indexs(1)+nn)   = taskid_i0
            int_mb(dm_indexs(1)+nn+1) = k0
            x = dbl_mb(rgrid(1)+3*(k0-1))
            y = dbl_mb(rgrid(1)+3*(k0-1)+1)
            z = dbl_mb(rgrid(1)+3*(k0-1)+2)
            !write(*,*) "ms,i,taskid_i,k0,x,y,z=",ms,i,taskid_i0,k0,x,y,z
         end do
      end do

*     **** generate summ matrix ****
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      call Dneall_m_zero(0,dbl_mb(summ(1)))
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
            i = ii-n1(ms)+1
            nn = (i-1)*2 + (ms-1)*2*(ne(1)+ne(2))
            taskid_i0 = int_mb(dm_indexs(1)+nn)   
            i0        = int_mb(dm_indexs(1)+nn+1)   
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
            nn = (i-1)*2 + (ms-1)*2*(ne(1)+ne(2))
            taskid_i0 = int_mb(dm_indexs(1)+nn)   
            i0        = int_mb(dm_indexs(1)+nn+1)   
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
      !write(*,*) "dn patch="
      !do i=1,ne(1)
      !   write(*,*) (dbl_mb(summ(1)+(i-1) +(j-1)*ne(1)),j=1,ne(1))
      !end do
      !write(*,*)

*     **** generate localized psi_r2 from cholesky ****
      call Dneall_m_cholesky(0,dbl_mb(summ(1)))
      !write(*,*) "L="
      !do i=1,ne(1)
      !   write(*,*) (dbl_mb(summ(1)+(i-1) +(j-1)*ne(1)),j=1,ne(1))
      !end do
      !write(*,*)
      !write(*,*)
      !write(*,*) "n2ft3d=",n2ft3d
      call Dneall_mg_forwardsolve(0,dbl_mb(summ(1)),n2ft3d,psi_r2)
      !write(*,*) "out forward_solve"

      !call D3dB_nx(1,i)
      !call D3dB_ny(1,j)
      !call D3dB_nz(1,ii)
      !z = 1.0d0/dble(i*j*ii)
      !do j=1,ne(1)
      !do i=1,ne(1)
      !   call D3dB_rr_dot(1,psi_r2(1,i),psi_r2(1,j),x)
      !   call D3dB_rr_dot(1,psi_r(1,i),psi_r(1,j),y)
      !   write(*,*) "i,j,x,y=",i,j,x*z,y*z
      !end do
      !end do

*     **** compute overlap of psi_r2 ****
     

      if (.not.Dneall_m_pop_stack(summ))
     >   call errquit('psi_dmatrix_localize:pop stack',0,MA_ERR)
      if (.not.BA_pop_stack(rgrid(2)))
     >   call errquit('psi_dmatrix_localize:pop stack',1,MA_ERR)
      if (.not.BA_pop_stack(dm_indexs(2))) 
     >   call errquit('psi_dmatrix_localize:pop stack',2,MA_ERR)
     
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

   

