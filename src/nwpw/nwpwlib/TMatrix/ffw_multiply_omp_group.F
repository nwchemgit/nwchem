c     ********************************************
c     *                                          *
c     *            zgemm2c_omp_group             *
c     *                                          *
c     ********************************************
c
c     Entry - group
c           - tid,nthr
c           - n,m,npack
c           - A,B
c           - thrmatrix: private memory
c           - C: assummed to be zero
c     Exit - C
c
      subroutine zgemm2c_omp_group(group,tid,nthr,n,m,npack2,
     >                             A,B,C,thrmatrix)
      implicit none
      integer group,tid,nthr,n,m,npack2
      complex*16 A(*),B(*),C(*),thrmatrix(*)

*     **** local variables ****
      integer r1,s1,shift1,npack2_block


*     **** define blocking ****
      s1 = npack2/nthr
      r1 = mod(npack2,nthr)
      shift1  = 1 + tid*s1 + min(tid,r1)
      npack2_block = s1 
      
      if (tid.lt.r1) npack2_block = npack2_block + 1

      call ZGEMM('C','N',n,m,npack2_block,
     >           2.0d0,
     >           A(shift1),npack2,
     >           B(shift1),npack2,
     >           0.0d0,
     >           thrmatrix,n)

      if (group.eq.1) then
!$OMP CRITICAL (group1)
         call zaxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (group1)

      else if (group.eq.2) then
!$OMP CRITICAL (ggroup2)
         call zaxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (ggroup2)

      else if (group.eq.3) then
!$OMP CRITICAL (hhhroup3)
         call zaxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (hhhroup3)
      end if

      return
      end
