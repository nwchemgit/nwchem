c
c$Id: pmlqsp.F,v 1.2 2000-07-27 16:09:20 bjohnson Exp $
c
      SUBROUTINE drdy_pmlqsp
C
C Calculate second-order term in Taylor series of x(s) at the saddle
C   point
C
      Implicit Double Precision (A-H,O-Z)
C
#include "drdyP.fh"
      integer fu6 
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
      call drdy_pot2 (x,v,vec1sv,f,hess,scr1,amass,natom,n3,n3tm,nbar)
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
      fu6 = 6
c
c!!! BGJ - The following call does not match the definition in mxnleq.F.
c!!! BGJ - I removed the last argument to allow it to link under WIN32.
c!!! BGJ - Someone who knows what this code should make sure this is
c!!! BGJ - the correct argument to delete.  It doesn't look like it based
c!!! BGJ - on the mxlneq code but I'm not sure what the proper fix is.
c!!!      call drdy_mxlneq(u0,n3,n3tm,det,jrank,eps,mode,-1,n3+1)
      call errquit('Bad call in pmlqsp.F - see code comments',0)
      call drdy_mxlneq(u0,n3,n3tm,det,jrank,eps,mode,-1)
      if (jrank .lt. n3) write (fu6,*) 'jrank = ', jrank
c
      do i = 1, n3
         vec1sv(i) = u0(i,n3+1)
      enddo
c
      return
      end
