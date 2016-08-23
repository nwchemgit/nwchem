      program testab
      implicit none
      integer npack,N
      parameter (npack=12000,N=256)
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

!$OMP PARALLEL private(it,tid,nthr,Itid,Mtid,r) shared(A,B,C)
      tid  = omp_get_thread_num()
      nthr = omp_get_num_threads()
      write(*,*) "into dgemm",tid,nthr
      do it=1,10
         call  dgemmtest(npack,N,N,A,B,C)
      end do
!$OMP END PARALLEL
      call current_second(t1)
      write(*,*) "time=",t1-t0
      return
      end

      subroutine dgemmtest(M,N,K,A,B,C)
      implicit none
      integer M,N,K
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)
      integer r,tid,nthr,Itid,Mtid
      integer  omp_get_num_threads, omp_get_thread_num
      external omp_get_num_threads, omp_get_thread_num

      tid  = omp_get_thread_num()
      nthr = omp_get_num_threads()
      r    = mod(M,nthr)
      if (tid.lt.r) then
          Itid = tid*(M/nthr+1) + 1
          !Itid = tid*(M/nthr+1)*N + 1
          Mtid = M/nthr + 1
      else
          Itid = r + tid*(M/nthr) + 1
          !Itid = (r + tid*(M/nthr))*N + 1
          Mtid = M/nthr
      end if
      write(*,*) tid,"  Itid,Mtid,N,K=",Itid,Mtid,N,K
      call DGEMM('N','N',Mtid,N,K,
     >        1.0d0,
     >        A(Itid),M,
     >        B,K,
     >        0.0d0,
     >        C(Itid),M)
     
!$OMP BARRIER

       return
       end
