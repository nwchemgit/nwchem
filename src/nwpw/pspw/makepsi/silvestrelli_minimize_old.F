*
* $Id$
*


*     *********************************
*     *                               *
*     *      silvestrelli_minimize_old    *
*     *                               *
*     *********************************
*
*   This routine minimizes
*

      subroutine silvestrelli_minimize_old(n,X,Y,Z,A)
      implicit none
      integer n
      complex*16 X(n,n),Y(n,n),Z(n,n)
      complex*16 A(n,n)

#include "bafdecls.fh"

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
      value = BA_push_get(mt_dcpl,(n*n),'At',
     >                    At(2),At(1))
      value = BA_push_get(mt_dcpl,(n*n),'Att',
     >                    Att(2),Att(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(n*n),'expplusA',
     >                    expplusA(2),expplusA(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(n*n),'expminusA',
     >                    expminusA(2),expminusA(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(n*n),'expminusAtotal',
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
           A(i,j)= X(j,i)        *(dconjg(X(j,j))-dconjg(X(i,i))) 
     >	         - dconjg(X(i,j))*(X(i,i)-X(j,j)) 
     >           + Y(j,i)        *(dconjg(Y(j,j))-dconjg(Y(i,i))) 
     >           - dconjg(Y(i,j))*(Y(i,i)-Y(j,j)) 
     >           + Z(j,i)        *(dconjg(Z(j,j))-dconjg(Z(i,i))) 
     >           - dconjg(Z(i,j))*(Z(i,i)-Z(j,j))      
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
     >              X,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              X,n)

         call ZGEMM('N','N',n,n,n,one,
     >              Y,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              Y,n)

         call ZGEMM('N','N',n,n,n,one,
     >              Z,n,
     >              dcpl_mb(expplusA(1)),n,
     >              zero,
     >              dcpl_mb(At(1)),n)
         call ZGEMM('N','N',n,n,n,one,
     >              dcpl_mb(expminusA(1)),n,
     >              dcpl_mb(At(1)),n,
     >              zero,
     >              Z,n)



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
            wx = dble(X(i,i))**2 + dimag(X(i,i))**2
            wy = dble(Y(i,i))**2 + dimag(Y(i,i))**2
            wz = dble(Z(i,i))**2 + dimag(Z(i,i))**2
            omega = omega + (wx+wy+wz)
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
      value =           BA_pop_stack(expminusAtotal(2))
      value = value.and.BA_pop_stack(expminusA(2))
      value = value.and.BA_pop_stack(expplusA(2))
      value = value.and.BA_pop_stack(Att(2))
      value = value.and.BA_pop_stack(At(2))  
      if (.not. value) 
     >  call errquit('silvestrelli_minimize:popping stack memory',1,0)


      return
      end
