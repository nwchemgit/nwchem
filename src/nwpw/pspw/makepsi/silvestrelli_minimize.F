*
* $Id$
*


*     *********************************
*     *                               *
*     *      silvestrelli_minimize    *
*     *                               *
*     *********************************
*
*   This routine minimizes
*

      subroutine silvestrelli_minimize(maxiter,ms,n,rank,wts,
     >                                 X1,X2,X3,X4,X5,X6,A)
      implicit none
      integer maxiter,ms,n,rank
      real*8 wts(*)
      complex*16 X1(N,N),X2(N,N),X3(N,N),X4(N,N),X5(N,N),X6(N,N)
      complex*16 A(n,n)

#include "bafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,oprint
      integer i,j,k,sign,iter,nn

      integer At(2),Att(2)
      integer expplusA(2),expminusA(2),expminusAtotal(2)

      real*8     timestep,maxgrad,omega,omega1,crit
      real*8     wx,wy,wz

      complex*16 one,zero,scal

      real*8   control_wannier_timestep,Dneall_w_tracesqr
      real*8   Dneall_w_max
      external control_wannier_timestep,Dneall_w_tracesqr
      external Dneall_w_max

      logical  Dneall_w_push_get,Dneall_w_pop_stack
      external Dneall_w_push_get,Dneall_w_pop_stack


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)
      call Dneall_m_size(ms,nn)

*     **** allocate space from stack ****
      value =           Dneall_w_push_get(ms,At)
      value = value.and.Dneall_w_push_get(ms,Att)
      value = value.and.Dneall_w_push_get(ms,expplusA)
      value = value.and.Dneall_w_push_get(ms,expminusA)
      value = value.and.Dneall_w_push_get(ms,expminusAtotal)
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:out of stack memory',0,0)


      call Dneall_w_eye(ms,dcpl_mb(expminusAtotal(1)),one)

      maxgrad  = 9999.0d0
      timestep = control_wannier_timestep()
      iter = 0
      crit   = 9999.0d0
      omega1 = 0.0d0
      

      do while( (iter.lt.10)    .or.
     >          (((crit.gt.1.0d-7).or.(maxgrad.gt.1.0d-08))
     >           .and.(iter.lt.maxiter)))

         iter=iter+1

         do i=1,n
         do j=1,n
           A(i,j)= wts(1)*( 
     >             X1(j,i)* (dconjg(X1(j,j))-dconjg(X1(i,i)) ) 
     >           - dconjg(X1(i,j))*( X1(i,i)-X1(j,j))) 
     >           + wts(2)*(X2(j,i)*(dconjg(X2(j,j))-dconjg(X2(i,i))) 
     >           - dconjg(X2(i,j))*(X2(i,i)-X2(j,j))) 
     >           + wts(3)*(X3(j,i)*(dconjg(X3(j,j))-dconjg(X3(i,i))) 
     >           - dconjg(X3(i,j))*(X3(i,i)-X3(j,j)))
           if (rank.gt.3) then
            A(i,j)= A(i,j)+
     >             wts(4)*( X4(j,i)*(dconjg(X4(j,j))-dconjg(X4(i,i))) 
     >           - dconjg(X4(i,j))*(X4(i,i)-X4(j,j)))
            if (rank.gt.4) then 
             A(i,j)= A(i,j)+
     >             wts(5)*(X5(j,i)*(dconjg(X5(j,j))-dconjg(X5(i,i))) 
     >           - dconjg(X5(i,j))*(X5(i,i)-X5(j,j)))
             if (rank.eq.6) then
              A(i,j)=A(i,j)+ 
     >             wts(6)*(X6(j,i)*(dconjg(X6(j,j))-dconjg(X6(i,i))) 
     >           - dconjg(X6(i,j))*(X6(i,i)-X6(j,j)))
             end if
            end if
           end if      
         end do
         end do

         !A = conjg(A)
         !A = timestep * A
         do j=1,n
         do i=1,n
           A(i,j) = timestep*dconjg(A(i,j))
         end do
         end do

         call Dneall_w_eye(ms,dcpl_mb(At(1)),one)
         call Dneall_w_eye(ms,dcpl_mb(expplusA(1)),one)
         call Dneall_w_eye(ms,dcpl_mb(expminusA(1)),one)

         sign = 1             
         do k=1,12
            scal = dcmplx(1.0d0/dble(k),0.0d0)
            call Dneall_www_Multiply(ms,
     >                       dcpl_mb(At(1)),
     >                       A,scal,
     >                       dcpl_mb(Att(1)),zero)
            call Dneall_ww_copy(ms,dcpl_mb(Att(1)),
     >                             dcpl_mb(At(1)))

            sign = -sign
            scal = dcmplx(dble(sign),0.0d0)
            call zaxpy(nn,one, dcpl_mb(At(1)),1,
     >                          dcpl_mb(expplusA(1)),1)
            call zaxpy(nn,scal,dcpl_mb(At(1)),1,
     >                          dcpl_mb(expminusA(1)),1) 
         end do



         !X = matmul(expminusA , matmul(X , expplusA))
         !Y = matmul(expminusA , matmul(Y , expplusA))
         !Z = matmul(expminusA , matmul(Z , expplusA))

         call Dneall_www_Multiply(ms,
     >                       X1,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X1,zero)
         call Dneall_www_Multiply(ms,
     >                       X2,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X2,zero)
         call Dneall_www_Multiply(ms,
     >                       X3,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X3,zero)

         if (rank.gt.3) then
         call Dneall_www_Multiply(ms,
     >                       X4,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X4,zero)
         if (rank.gt.4) then
         call Dneall_www_Multiply(ms,
     >                       X5,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X5,zero)
         if (rank.gt.5) then
         call Dneall_www_Multiply(ms,
     >                       X6,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X6,zero)
         end if
         end if
         end if 

         !expminusAtotal = matmul(expminusA,expminusAtotal)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(expminusAtotal(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_ww_copy(ms,dcpl_mb(At(1)),
     >                          dcpl_mb(expminusAtotal(1)))

         maxgrad =  Dneall_w_max(ms,A)

         omega = 0.0d0
         omega = omega + wts(1)*Dneall_w_tracesqr(ms,X1)
         omega = omega + wts(2)*Dneall_w_tracesqr(ms,X2)
         omega = omega + wts(3)*Dneall_w_tracesqr(ms,X3)
         if (rank.gt.3) omega = omega + wts(4)*Dneall_w_tracesqr(ms,X4)
         if (rank.gt.4) omega = omega + wts(5)*Dneall_w_tracesqr(ms,X5)
         if (rank.gt.5) omega = omega + wts(6)*Dneall_w_tracesqr(ms,X6)

         if (omega.lt.omega1) timestep = 0.5d0*timestep

         crit = dabs(omega-omega1)/omega
         omega1 = omega

         if (oprint) then
            write(6,'(I5,4E15.6)') iter,omega,maxgrad,crit,timestep
         end if
      end do ! end do while

      if (oprint) then
        write(6,1000)
        write(6,1001) iter, omega, maxgrad,crit
 1000 FORMAT(//' fixed step steepest descent iteration results:')
 1001 FORMAT(5x,'iterations performed=',I12,
     >       /5x,'               omega=',F12.6,
     >       /5x,'               errors=',2E12.3)

      end if

      call zcopy(nn,dcpl_mb(expminusAtotal(1)),1,A,1)

*     **** pop memory ***
      value =           Dneall_w_pop_stack(expminusAtotal)
      value = value.and.Dneall_w_pop_stack(expminusA)
      value = value.and.Dneall_w_pop_stack(expplusA)
      value = value.and.Dneall_w_pop_stack(Att)
      value = value.and.Dneall_w_pop_stack(At)
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:popping stack memory',1,0)


      return
      end




*     *********************************
*     *                               *
*     *      silvestrelli_iterate     *
*     *                               *
*     *********************************
*
*   This routine minimizes
*

      subroutine silvestrelli_iterate(maxiter,ms,n,rank,wts,
     >                                 X1,X2,X3,X4,X5,X6,A)
      implicit none
      integer maxiter,ms,n,rank
      real*8 wts(*)
      complex*16 X1(N,N),X2(N,N),X3(N,N),X4(N,N),X5(N,N),X6(N,N)
      complex*16 A(n,n)

#include "bafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,oprint
      integer i,j,k,sign,iter,nn

      integer At(2),Att(2)
      integer expplusA(2),expminusA(2),expminusAtotal(2)

      real*8     timestep,maxgrad,omega,omega1,crit
      real*8     wx,wy,wz

      complex*16 one,zero,scal

      real*8   control_wannier_timestep,Dneall_w_tracesqr
      real*8   Dneall_w_max
      external control_wannier_timestep,Dneall_w_tracesqr
      external Dneall_w_max

      logical  Dneall_w_push_get,Dneall_w_pop_stack
      external Dneall_w_push_get,Dneall_w_pop_stack


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)
      oprint = .false.

      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)
      call Dneall_m_size(ms,nn)

*     **** allocate space from stack ****
      value =           Dneall_w_push_get(ms,At)
      value = value.and.Dneall_w_push_get(ms,Att)
      value = value.and.Dneall_w_push_get(ms,expplusA)
      value = value.and.Dneall_w_push_get(ms,expminusA)
      value = value.and.Dneall_w_push_get(ms,expminusAtotal)
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:out of stack memory',0,0)


      call Dneall_w_eye(ms,dcpl_mb(expminusAtotal(1)),one)

      maxgrad  = 9999.0d0
      timestep = control_wannier_timestep()
      iter = 0
      crit   = 9999.0d0
      omega1 = 0.0d0
      

      do while (((crit.gt.1.0d-7).or.(maxgrad.gt.1.0d-08))
     >           .and.(iter.lt.maxiter))

         iter=iter+1

         do i=1,n
         do j=1,n
           A(i,j)= wts(1)*( 
     >             X1(j,i)* (dconjg(X1(j,j))-dconjg(X1(i,i)) ) 
     >           - dconjg(X1(i,j))*( X1(i,i)-X1(j,j))) 
     >           + wts(2)*(X2(j,i)*(dconjg(X2(j,j))-dconjg(X2(i,i))) 
     >           - dconjg(X2(i,j))*(X2(i,i)-X2(j,j))) 
     >           + wts(3)*(X3(j,i)*(dconjg(X3(j,j))-dconjg(X3(i,i))) 
     >           - dconjg(X3(i,j))*(X3(i,i)-X3(j,j)))
           if (rank.gt.3) then
            A(i,j)= A(i,j)+
     >             wts(4)*( X4(j,i)*(dconjg(X4(j,j))-dconjg(X4(i,i))) 
     >           - dconjg(X4(i,j))*(X4(i,i)-X4(j,j)))
            if (rank.gt.4) then 
             A(i,j)= A(i,j)+
     >             wts(5)*(X5(j,i)*(dconjg(X5(j,j))-dconjg(X5(i,i))) 
     >           - dconjg(X5(i,j))*(X5(i,i)-X5(j,j)))
             if (rank.eq.6) then
              A(i,j)=A(i,j)+ 
     >             wts(6)*(X6(j,i)*(dconjg(X6(j,j))-dconjg(X6(i,i))) 
     >           - dconjg(X6(i,j))*(X6(i,i)-X6(j,j)))
             end if
            end if
           end if      
         end do
         end do

         !A = conjg(A)
         !A = timestep * A
         do j=1,n
         do i=1,n
           A(i,j) = timestep*dconjg(A(i,j))
         end do
         end do

         call Dneall_w_eye(ms,dcpl_mb(At(1)),one)
         call Dneall_w_eye(ms,dcpl_mb(expplusA(1)),one)
         call Dneall_w_eye(ms,dcpl_mb(expminusA(1)),one)

         sign = 1             
         do k=1,12
            scal = dcmplx(1.0d0/dble(k),0.0d0)
            call Dneall_www_Multiply(ms,
     >                       dcpl_mb(At(1)),
     >                       A,scal,
     >                       dcpl_mb(Att(1)),zero)
            call Dneall_ww_copy(ms,dcpl_mb(Att(1)),
     >                             dcpl_mb(At(1)))

            sign = -sign
            scal = dcmplx(dble(sign),0.0d0)
            call zaxpy(nn,one, dcpl_mb(At(1)),1,
     >                          dcpl_mb(expplusA(1)),1)
            call zaxpy(nn,scal,dcpl_mb(At(1)),1,
     >                          dcpl_mb(expminusA(1)),1) 
         end do



         !X = matmul(expminusA , matmul(X , expplusA))
         !Y = matmul(expminusA , matmul(Y , expplusA))
         !Z = matmul(expminusA , matmul(Z , expplusA))

         call Dneall_www_Multiply(ms,
     >                       X1,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X1,zero)
         call Dneall_www_Multiply(ms,
     >                       X2,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X2,zero)
         call Dneall_www_Multiply(ms,
     >                       X3,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X3,zero)

         if (rank.gt.3) then
         call Dneall_www_Multiply(ms,
     >                       X4,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X4,zero)
         if (rank.gt.4) then
         call Dneall_www_Multiply(ms,
     >                       X5,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X5,zero)
         if (rank.gt.5) then
         call Dneall_www_Multiply(ms,
     >                       X6,
     >                       dcpl_mb(expplusA(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(At(1)),one,
     >                       X6,zero)
         end if
         end if
         end if 

         !expminusAtotal = matmul(expminusA,expminusAtotal)
         call Dneall_www_Multiply(ms,
     >                       dcpl_mb(expminusA(1)),
     >                       dcpl_mb(expminusAtotal(1)),one,
     >                       dcpl_mb(At(1)),zero)
         call Dneall_ww_copy(ms,dcpl_mb(At(1)),
     >                          dcpl_mb(expminusAtotal(1)))

         maxgrad =  Dneall_w_max(ms,A)

         omega = 0.0d0
         omega = omega + wts(1)*Dneall_w_tracesqr(ms,X1)
         omega = omega + wts(2)*Dneall_w_tracesqr(ms,X2)
         omega = omega + wts(3)*Dneall_w_tracesqr(ms,X3)
         if (rank.gt.3) omega = omega + wts(4)*Dneall_w_tracesqr(ms,X4)
         if (rank.gt.4) omega = omega + wts(5)*Dneall_w_tracesqr(ms,X5)
         if (rank.gt.5) omega = omega + wts(6)*Dneall_w_tracesqr(ms,X6)

         if (omega.lt.omega1) timestep = 0.5d0*timestep

         crit = dabs(omega-omega1)/omega
         omega1 = omega

         if (oprint) then
            write(6,'(I5,4E15.6)') iter,omega,maxgrad,crit,timestep
         end if
      end do ! end do while

      if (oprint) then
        write(6,1000)
        write(6,1001) iter, omega, maxgrad,crit
 1000 FORMAT(//' fixed step steepest descent iteration results:')
 1001 FORMAT(5x,'iterations performed=',I12,
     >       /5x,'               omega=',F12.6,
     >       /5x,'               errors=',2E12.3)

      end if

      call zcopy(nn,dcpl_mb(expminusAtotal(1)),1,A,1)

*     **** pop memory ***
      value =           Dneall_w_pop_stack(expminusAtotal)
      value = value.and.Dneall_w_pop_stack(expminusA)
      value = value.and.Dneall_w_pop_stack(expplusA)
      value = value.and.Dneall_w_pop_stack(Att)
      value = value.and.Dneall_w_pop_stack(At)
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:popping stack memory',1,0)


      return
      end
