c
c$Id$
c
      Subroutine drdy_pmlqsp
C
C Calculate second-order term in Taylor series of x(s) at the saddle
C   point
C
      Implicit none
#include "errquit.fh"
C
#include "drdyP.fh"
      Integer i,j,jrank,neq
      Double Precision det,d3lx,tdx,sum1,sum2,vfv,vgvv
C
      Double Precision eps
C
      Data eps/1.d-9/
c
c Save initial value of x and make initial step 
c  on the product side
      d3lx = 0.1d00*dels
      do i=1,n3
         x(i)=xspsv(i)+vec0sv(i)*d3lx
      enddo
C Calculate the force constant matrix
      call drdy_pot2 (x,v,vec1sv,fold,hess,scr1,
     &      amass,natom,n3,n3tm,nbardm)
C Initial step on the reactant side
      do i=1,n3
         x(i)=xspsv(i)-vec0sv(i)*d3lx
      enddo
C Calculate the force constant matrix
      call drdy_pot2 (x,v,vec1sv,f,hess,scr1,amass,natom,n3,n3tm,nbardm)
c Calculate third derivatives of the energy with respect to x with 
C   and use it to calculate vector with vec0sv.
      tdx = 2.0d00*d3lx
      vgvv = 0.0d00
      vfv = 0.0d00
      do i = 1, n3
         sum1 = 0.0d00
         sum2 = 0.0d00
         do j = 1, n3
            sum1 = sum1 + vec0sv(j)*(fold(i,j)-f(i,j))/tdx
            sum2 = sum2 + fspsv(i,j)*vec0sv(j)
         enddo
         vec1sv(i) = -sum1
         vgvv = vgvv + vec0sv(i)*sum1
         vfv = vfv + vec0sv(i)*sum2
      enddo
c Construct matrix for linear equations
      do i = 1, n3
         do j = 1, n3
            u0(i,j) = fspsv(i,j)
         enddo
         u0(i,i) = u0(i,i) - 2.0d00*vfv
         u0(i,n3+1) = vec1sv(i) + vgvv*vec0sv(i)
      enddo
c
      neq = 1
      call drdy_mxlneq(u0,n3,n3tm,det,jrank,eps,mode,neq)
      if (jrank .lt. n3) then
         write (fu6,6000) jrank
         call errquit('Problem in pmlqsp.F',0, UNKNOWN_ERR)
      endif
c
      do i = 1, n3
         vec1sv(i) = u0(i,n3+1)
      enddo
c
      return
6000  Format(' Problem with mxlneq in subroutine pmlqsp, jrank = ',
     * i5)
      end
