c     ********************************************
c     *                                          *
c     *            dgemm2c_omp_group             *
c     *                                          *
c     ********************************************
c
c     Entry - group
c           - tid,nthr
c           - n,m,npack,nida2
c           - A,B
c           - thrmatrix: private memory
c           - C: assummed to be zero
c     Exit - C
c
      subroutine dgemm2c_omp_group(group,tid,nthr,n,m,npack2,nida2,
     >                             A,B,C,thrmatrix)
      implicit none
      integer group,tid,nthr,n,m,npack2,nida2
      real*8 A(*),B(*),C(*),thrmatrix(*)

*     **** local variables ****
      integer r1,s1,r2,s2,shift1,shift1a,npack2_block,nida2_block,i


*     **** define blocking ****
      s1 = npack2/nthr
      r1 = mod(npack2,nthr)
      s2 = nida2/nthr
      r2 = mod(nida2,nthr)
      shift1  = 1 + tid*s1 + min(tid,r1)
      shift1a = 1 + tid*s2 + min(tid,r2)
      npack2_block = s1 
      nida2_block = s2 
      if (tid.lt.r1) npack2_block = npack2_block + 1
      if (tid.lt.r2) nida2_block = nida2_block + 1

      call DGEMM('T','N',n,m,npack2_block,
     >           2.0d0,
     >           A(shift1),npack2,
     >           B(shift1),npack2,
     >           0.0d0,
     >           thrmatrix,n)
      if (nida2_block.gt.0) then
         call DGEMM('T','N',n,m,nida2_block,
     >              -1.0d0,
     >              A(shift1a),npack2,
     >              B(shift1a),npack2,
     >              1.0d0,
     >              thrmatrix,n)
      end if

      if (group.eq.1) then
!$OMP CRITICAL (group1)
         call daxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (group1)

      else if (group.eq.2) then
!$OMP CRITICAL (ggroup2)
         call daxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (ggroup2)

      else if (group.eq.3) then
!$OMP CRITICAL (hhhroup3)
         call daxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (hhhroup3)
      end if

      return
      end


c     ********************************************
c     *                                          *
c     *            dgemm2c0_omp_group            *
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
      subroutine dgemm2c0_omp_group(group,tid,nthr,n,m,npack2,
     >                             A,B,C,thrmatrix)
      implicit none
      integer group,tid,nthr,n,m,npack2
      real*8 A(*),B(*),C(*),thrmatrix(*)

*     **** local variables ****
      integer r1,s1,shift1,npack2_block,i


*     **** define blocking ****
      s1 = npack2/nthr
      r1 = mod(npack2,nthr)
      shift1  = 1 + tid*s1 + min(tid,r1)
      npack2_block = s1
      if (tid.lt.r1) npack2_block = npack2_block + 1

      call DGEMM('T','N',n,m,npack2_block,
     >           1.0d0,
     >           A(shift1),npack2,
     >           B(shift1),npack2,
     >           0.0d0,
     >           thrmatrix,n)

      if (group.eq.1) then
!$OMP CRITICAL (group1)
         call daxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (group1)

      else if (group.eq.2) then
!$OMP CRITICAL (ggroup2)
         call daxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (ggroup2)

      else if (group.eq.3) then
!$OMP CRITICAL (hhhroup3)
         call daxpy(n*m,1.0d0,thrmatrix,1,C,1)
!$OMP END CRITICAL (hhhroup3)
      end if

      return
      end
