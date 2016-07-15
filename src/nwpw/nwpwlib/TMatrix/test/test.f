      program testab
      implicit none
      integer npack,N
      parameter (npack=20000,N=256)
      real*8 A(npack,N)
      real*8 C(npack,N)
      real*8 B(N,N)
      real*8 t0,t1
      integer i,j,it
      integer Mtid,Itid,tid,r,nthr,M
      integer  omp_get_num_threads, omp_get_thread_num
      external omp_get_num_threads, omp_get_thread_num

      do j=1,N
      do i=1,N
         B(i,j) = 0.01d0*(i+j) + 1.3d0
      end do
      end do
      do j=1,N
      do i=1,npack
          A(i,j) = dcos(0.01d0*i)+dsin(0.23d0*j)
      end do
      end do
      
      call current_second(t0)
      write(*,*) "hera"

!$OMP PARALLEL private(it,tid,nthr,Itid,Mtid,r)
      tid  = omp_get_thread_num()
      nthr = omp_get_num_threads()

      do i=1,10
         write(*,*) "into dgemm"
         M = npack
         r    = mod(M,nthr)
         if (tid.lt.r) then
             Itid = tid*(M/nthr+1) + 1
             Mtid = M/nthr + 1
         else
             Itid = r + tid*(M/nthr) + 1
             Mtid = M/nthr
         end if

         call DGEMM('N','N',npack,N,N,
     >           1.0d0,
     >           A(Itid,1),npack,
     >           B,N,
     >           0.0d0,
     >           C(Itid,1),npack)

         write(*,*) "out dgemm"
      end do
!$OMP END PARALLEL
      call current_second(t1)
      write(*,*) "time=",t1-t0
      return
      end
     
