*
* $Id: silvestrelli_minimize.F,v 1.9 2007-12-27 20:18:45 d3p307 Exp $
*


*     *********************************
*     *                               *
*     *      silvestrelli_minimize    *
*     *                               *
*     *********************************
*
*   This routine minimizes
*

      subroutine silvestrelli_minimize(n,rank,wts,X1,X2,X3,X4,X5,X6,A)
      implicit none
      integer n,rank
      real*8 wts(*)
      complex*16 X1(N,N),X2(N,N),X3(N,N),X4(N,N),X5(N,N),X6(N,N)
      complex*16 A(n,n)

#include "mafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      integer MAXITER
      parameter (MAXITER=10000)

      logical value,oprint
      integer i,j,k,sign,iter

      integer At(2),Att(2)
      integer expplusA(2),expminusA(2),expminusAtotal(2)

      real*8     timestep,maxgrad,omega,omega1,crit
      real*8     wx,wy,wz

      complex*16 one,zero,scal


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)

*     **** allocate space from stack ****
      value = MA_push_get(mt_dcpl,(n*n),'At',
     >                    At(2),At(1))
      value = MA_push_get(mt_dcpl,(n*n),'Att',
     >                    Att(2),Att(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(n*n),'expplusA',
     >                    expplusA(2),expplusA(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(n*n),'expminusA',
     >                    expminusA(2),expminusA(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(n*n),'expminusAtotal',
     >                    expminusAtotal(2),expminusAtotal(1))
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:out of stack memory',0,0)


      call zcopy(n*n,zero,0,dcpl_mb(expminusAtotal(1)),1)
      do i=1,n
        dcpl_mb(expminusAtotal(1)+(i-1)+(i-1)*n) = one
      end do



      maxgrad  = 9999.0d0
      timestep =    3.0d-2
      iter = 0
      crit   = 9999.0d0
      omega1 = 0.0d0
      

      do while( (iter.lt.10)    .or.
     >          (crit.gt.1.0d-7).or.
     >          ((maxgrad .GT. 1.0d-08).and.(iter.lt.MAXITER)))

         iter=iter+1

         do i=1,n
         do j=1,n
           A(i,j)= wts(1)*( 
     >             X1(j,i)* (dconjg(X1(j,j))-dconjg(X1(i,i)) ) 
     >	         - dconjg(X1(i,j))*( X1(i,i)-X1(j,j))) 
     >           + wts(2)*(X2(j,i)*(dconjg(X2(j,j))-dconjg(X2(i,i))) 
     >           - dconjg(X2(i,j))*(X2(i,i)-X2(j,j))) 
     >           + wts(3)*(X3(j,i)*(dconjg(X3(j,j))-dconjg(X3(i,i))) 
     >           - dconjg(X3(i,j))*(X3(i,i)-X3(j,j)))
           if (rank.gt.3) then
           A(i,j)= A(i,j)+
     >             wts(4)*( X4(j,i)*(dconjg(X4(j,j))-dconjg(X4(i,i))) 
     >	         - dconjg(X4(i,j))*(X4(i,i)-X4(j,j)))

            if (rank.gt.4) then 
            A(i,j)= A(i,j)+
     >            wts(5)*(X5(j,i)*(dconjg(X5(j,j))-dconjg(X5(i,i))) 
     >           - dconjg(X5(i,j))*(X5(i,i)-X5(j,j)))
           
             if (rank.eq.6) then
             A(i,j)=A(i,j)+ 
     >            wts(6)*(X6(j,i)*(dconjg(X6(j,j))-dconjg(X6(i,i))) 
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

         call zcopy(n*n,zero,0,dcpl_mb(At(1)),1)
         call zcopy(n*n,zero,0,dcpl_mb(expplusA(1)),1)
         call zcopy(n*n,zero,0,dcpl_mb(expminusA(1)),1)
         do i=1,n
           dcpl_mb(At(1)       +(i-1)+(i-1)*n) = one
           dcpl_mb(expplusA(1) +(i-1)+(i-1)*n) = one
           dcpl_mb(expminusA(1)+(i-1)+(i-1)*n) = one
         end do

         sign = 1             
         do k=1,12
            scal = dcmplx(1.0d0/dble(k),0.0d0)
            call ZGEMM('N','N',n,n,n,scal,
     >                 dcpl_mb(At(1)),n,
     >                 A,n,
     >                 zero,
     >                 dcpl_mb(Att(1)),n)
            call zcopy(n*n,dcpl_mb(Att(1)),1,
     >                     dcpl_mb(At(1)),1)

            sign = -sign
            scal = dcmplx(dble(sign),0.0d0)
            call zaxpy(n*n,one, dcpl_mb(At(1)),1,
     >                          dcpl_mb(expplusA(1)),1)
            call zaxpy(n*n,scal,dcpl_mb(At(1)),1,
     >                          dcpl_mb(expminusA(1)),1) 
         end do



	 !X = matmul(expminusA , matmul(X , expplusA))
	 !Y = matmul(expminusA , matmul(Y , expplusA))
         !Z = matmul(expminusA , matmul(Z , expplusA))

         call ZGEMM('N','N',n,n,n,one,
     >              X1,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X1,n)

         call ZGEMM('N','N',n,n,n,one,
     >              X2,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X2,n)

         call ZGEMM('N','N',n,n,n,one,
     >              X3,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X3,n)

         if (rank.gt.3) then
         call ZGEMM('N','N',n,n,n,one,
     >              X4,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X4,n)
          if (rank.gt.4) then
         call ZGEMM('N','N',n,n,n,one,
     >              X5,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X5,n)
           if (rank.gt.5) then
         call ZGEMM('N','N',n,n,n,one,
     >              X6,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X6,n)
           end if
          end if
         end if 

         !expminusAtotal = matmul(expminusA,expminusAtotal)
	 call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(expminusAtotal(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call zcopy(n*n,dcpl_mb(At(1)),1,dcpl_mb(expminusAtotal(1)),1)       


         maxgrad=0.d0
         do i=1,n
         do j=1,n
           wx = dble(A(i,j))**2 + dimag(A(i,j))**2
           wx = dsqrt(wx)
           if(wx .GT. maxgrad) maxgrad=wx
         end do
         end do

         omega=0.d0
         do i=1,n
            wx = wts(1)*(dble(X1(i,i))**2 + dimag(X1(i,i))**2)
            wy = wts(2)*(dble(X2(i,i))**2 + dimag(X2(i,i))**2)
            wz = wts(3)*(dble(X3(i,i))**2 + dimag(X3(i,i))**2)
            omega = omega + (wx+wy+wz)
            if (rank.gt.3) then
             wx = wts(4)*(dble(X4(i,i))**2 + dimag(X4(i,i))**2)
             omega=omega+wx
            if (rank.gt.4) then
             wx = wts(5)*(dble(X5(i,i))**2 + dimag(X5(i,i))**2)
             omega=omega+wx
            if (rank.gt.5) then 
             wx = wts(6)*(dble(X6(i,i))**2 + dimag(X6(i,i))**2)
             omega = omega + wx
            end if
            end if 
            end if
         end do

         crit = dabs(omega-omega1)/omega
         omega1 = omega

         !if (oprint) then
         !write(6,*) iter, omega, maxgrad,crit
         !end if
      end do ! end do while

      if (oprint) then
        write(6,1000)
        write(6,1001) iter, omega, maxgrad,crit
 1000 FORMAT(//' fixed step steepest descent iteration results:')
 1001 FORMAT(5x,'iterations performed=',I12,
     >       /5x,'               omega=',F12.6,
     >       /5x,'               errors=',2E12.3)

      end if

      call zcopy(n*n,dcpl_mb(expminusAtotal(1)),1,A,1)

*     **** pop memory ***
      value =           MA_pop_stack(expminusAtotal(2))
      value = value.and.MA_pop_stack(expminusA(2))
      value = value.and.MA_pop_stack(expplusA(2))
      value = value.and.MA_pop_stack(Att(2))
      value = value.and.MA_pop_stack(At(2))  
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:popping stack memory',1,0)


      return
      end
