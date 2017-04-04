*
* $Id$
*

* Parallel.f
* Author - Eric Bylaska
*
*   These routines are to be used to keep track of the parallel message
* passing variables, as well as iniitialize and deinitialize the
* message passing routines.
*



*     *************************************
*     *                                   *
*     *        Parallel_Finalize          *
*     *                                   *
*     *************************************

      subroutine Parallel_Finalize()
      implicit none

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"


#ifdef USE_SUBGROUPS
#ifdef GA_DEV
#ifdef MPI
      integer myio,MASTER,mydepth
      parameter (MASTER=0)
      integer  util_sgroup_unset_io,util_sgroup_depth
      external util_sgroup_unset_io,util_sgroup_depth
      mydepth = util_sgroup_depth()
      if ((taskid.eq.MASTER).and.(mydepth.gt.0)) then
         myio = util_sgroup_unset_io()
      end if
#endif
#endif
#endif
      if (.not.BA_free_heap(procNd(2)))
     >  call errquit('Parralel_finalalize',1,MA_ERR)

      return
      end


*     *************************************
*     *                                   *
*     *        Parallel_np                *
*     *                                   *
*     *************************************

      subroutine Parallel_np(np_out)
      implicit none
      integer np_out

#include "Parallel.fh"

      np_out = np
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel_taskid            *
*     *                                   *
*     *************************************

      subroutine Parallel_taskid(task_out)
      implicit none
      integer task_out
      
#include "Parallel.fh"

      task_out = taskid
      return 
      end




*     *************************************
*     *                                   *
*     *        Parallel2d_np_i            *
*     *                                   *
*     *************************************
      subroutine Parallel2d_np_i(np_out)
      implicit none
      integer np_out

#include  "Parallel.fh"

      np_out = np_i
      return
      end

*     *************************************
*     *                                   *
*     *        Parallel2d_np_j            *
*     *                                   *
*     *************************************
      subroutine Parallel2d_np_j(np_out)
      implicit none
      integer np_out

#include  "Parallel.fh"

      np_out = np_j
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel2d_taskid_i        *
*     *                                   *
*     *************************************
      subroutine Parallel2d_taskid_i(taskid_out)
      implicit none
      integer taskid_out

#include "Parallel.fh"

      taskid_out = taskid_i
      return
      end 

*     *************************************
*     *                                   *
*     *        Parallel2d_taskid_j        *
*     *                                   *
*     *************************************
      subroutine Parallel2d_taskid_j(taskid_out)
      implicit none
      integer taskid_out

#include "Parallel.fh"

      taskid_out = taskid_j
      return
      end


*     *************************************
*     *                                   *
*     *     Parallel2d_convert_taskid_i   *
*     *                                   *
*     *************************************
      integer function Parallel2d_convert_taskid_i(i)
      implicit none
      integer i

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel2d_convert_taskid_i = int_mb(procNd(1)
     >                                     + i
     >                                     + taskid_j*np_i 
     >                                     + taskid_k*np_i*np_j)
      return
      end

*     *************************************
*     *                                   *
*     *     Parallel2d_convert_taskid_j   *
*     *                                   *
*     *************************************
      integer function Parallel2d_convert_taskid_j(j)
      implicit none
      integer j

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel2d_convert_taskid_j = int_mb(procNd(1)
     >                                    + taskid_i
     >                                    + j*np_i 
     >                                    + taskid_k*np_i*np_j)
      return
      end




*     *************************************
*     *                                   *
*     *   Parallel2d_convert_taskid_ij    *
*     *                                   *
*     *************************************
      integer function Parallel2d_convert_taskid_ij(i,j)
      implicit none
      integer i,j

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel2d_convert_taskid_ij = int_mb(procNd(1)+i+j*np_i 
     >                                    + taskid_k*np_i*np_j)
      return
      end



*     *************************************
*     *                                   *
*     *         Parallel2d_comm_i         *
*     *                                   *
*     *************************************
      integer function Parallel2d_comm_i()
      implicit none

#include "Parallel.fh"

      Parallel2d_comm_i = comm_i
      return
      end

*     *************************************
*     *                                   *
*     *         Parallel2d_comm_j         *
*     *                                   *
*     *************************************
      integer function Parallel2d_comm_j()
      implicit none

#include "Parallel.fh"

      Parallel2d_comm_j = comm_j
      return
      end




*     *************************************
*     *                                   *
*     *        Parallel3d_np_i            *
*     *                                   *
*     *************************************
      subroutine Parallel3d_np_i(np_out)
      implicit none
      integer np_out

#include  "Parallel.fh"

      np_out = np_i
      return
      end

*     *************************************
*     *                                   *
*     *        Parallel3d_np_j            *
*     *                                   *
*     *************************************
      subroutine Parallel3d_np_j(np_out)
      implicit none
      integer np_out

#include  "Parallel.fh"

      np_out = np_j
      return
      end

*     *************************************
*     *                                   *
*     *        Parallel3d_np_k            *
*     *                                   *
*     *************************************
      subroutine Parallel3d_np_k(np_out)
      implicit none
      integer np_out

#include  "Parallel.fh"

      np_out = np_k
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel3d_taskid_i        *
*     *                                   *
*     *************************************
      subroutine Parallel3d_taskid_i(taskid_out)
      implicit none
      integer taskid_out

#include "Parallel.fh"

      taskid_out = taskid_i
      return
      end 

*     *************************************
*     *                                   *
*     *        Parallel3d_taskid_j        *
*     *                                   *
*     *************************************
      subroutine Parallel3d_taskid_j(taskid_out)
      implicit none
      integer taskid_out

#include "Parallel.fh"

      taskid_out = taskid_j
      return
      end

*     *************************************
*     *                                   *
*     *        Parallel3d_taskid_k        *
*     *                                   *
*     *************************************
      subroutine Parallel3d_taskid_k(taskid_out)
      implicit none
      integer taskid_out

#include "Parallel.fh"

      taskid_out = taskid_k
      return
      end



*     *************************************
*     *                                   *
*     *     Parallel3d_convert_taskid_i   *
*     *                                   *
*     *************************************
      integer function Parallel3d_convert_taskid_i(i)
      implicit none
      integer i

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel3d_convert_taskid_i = int_mb(procNd(1)
     >                                    + i
     >                                    + taskid_j*np_i 
     >                                    + taskid_k*np_i*np_j)
      return
      end

*     *************************************
*     *                                   *
*     *     Parallel3d_convert_taskid_j   *
*     *                                   *
*     *************************************
      integer function Parallel3d_convert_taskid_j(j)
      implicit none
      integer j

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel3d_convert_taskid_j = int_mb(procNd(1)
     >                                    + taskid_i
     >                                    + j*np_i
     >                                    + taskid_k*np_i*np_j)
      return
      end

*     *************************************
*     *                                   *
*     *     Parallel3d_convert_taskid_k   *
*     *                                   *
*     *************************************
      integer function Parallel3d_convert_taskid_k(k)
      implicit none
      integer k

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel3d_convert_taskid_k = int_mb(procNd(1)
     >                                    + taskid_i
     >                                    + taskid_j*np_i
     >                                    + k*np_i*np_j)
      return
      end



*     *************************************
*     *                                   *
*     *   Parallel3d_convert_taskid_ijk   *
*     *                                   *
*     *************************************
      integer function Parallel3d_convert_taskid_ijk(i,j,k)
      implicit none
      integer i,j,k

#include "Parallel.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      Parallel3d_convert_taskid_ijk = int_mb(procNd(1)
     >                                      + i
     >                                      + j*np_i
     >                                      + k*np_i*np_j)

      return
      end



*     *************************************
*     *                                   *
*     *         Parallel3d_comm_i         *
*     *                                   *
*     *************************************
      integer function Parallel3d_comm_i()
      implicit none

#include "Parallel.fh"

      Parallel3d_comm_i = comm_i
      return
      end

*     *************************************
*     *                                   *
*     *         Parallel3d_comm_j         *
*     *                                   *
*     *************************************
      integer function Parallel3d_comm_j()
      implicit none

#include "Parallel.fh"

      Parallel3d_comm_j = comm_j
      return
      end

*     *************************************
*     *                                   *
*     *         Parallel3d_comm_k         *
*     *                                   *
*     *************************************
      integer function Parallel3d_comm_k()
      implicit none

#include "Parallel.fh"

      Parallel3d_comm_k = comm_k
      return
      end

*     *************************************
*     *                                   *
*     *         Parallel_comm_world       *
*     *                                   *
*     *************************************
      integer function Parallel_comm_world()
      implicit none

#include "Parallel.fh"

      Parallel_comm_world = comm_world
      return
      end





*     *************************************
*     *                                   *
*     *        Parallel_thread_init       *
*     *                                   *
*     *************************************
      subroutine Parallel_thread_init()
      implicit none
#include "Parallel.fh"

#ifdef USE_OPENMP
      integer  omp_get_thread_num,omp_get_num_threads
      external omp_get_thread_num,omp_get_num_threads
#endif

#ifdef USE_OPENMP
      threadid = omp_get_thread_num()
      nthreads = omp_get_num_threads()
#else
      threadid = 0
      nthreads = 1
#endif

      !write(*,*) "---threadid=",threadid
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel_threadid          *
*     *                                   *
*     *************************************
      integer function Parallel_threadid()
      implicit none
#include "Parallel.fh"
#ifdef USE_OPENMP
      integer  omp_get_thread_num
      external omp_get_thread_num
#endif

#ifdef USE_OPENMP
      Parallel_threadid = omp_get_thread_num()
#else
      Parallel_threadid = 0
#endif
      return
      end

*     *************************************
*     *                                   *
*     *        Parallel_nthreads          *
*     *                                   *
*     *************************************
      integer function Parallel_nthreads()
      implicit none
#include "Parallel.fh"
#ifdef USE_OPENMP
      integer  omp_get_num_threads
      external omp_get_num_threads
#endif

#ifdef USE_OPENMP
      Parallel_nthreads = omp_get_num_threads()
#else
      Parallel_nthreads = 1
#endif

      return
      end

*     *************************************
*     *                                   *
*     *        Parallel_maxthreads        *
*     *                                   *
*     *************************************
      integer function Parallel_maxthreads()
      implicit none
#include "Parallel.fh"
#ifdef USE_OPENMP
      integer  omp_get_max_threads
      external omp_get_max_threads
#endif

#ifdef USE_OPENMP
      Parallel_maxthreads = omp_get_max_threads()
#else
      Parallel_maxthreads = 1
#endif

      return
      end


*     *************************************
*     *                                   *
*     *    Parallel_shared_vector_copy    *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_vector_copy(barrier,n,a,b)
      implicit none
      logical barrier
      integer n
      real*8 a(*),b(*)

#ifdef USE_OPENMP
      integer nbig,tid,nthr,nsz,istart
      parameter (nbig=1000)

c     **** external functions ****
      integer  Parallel_nthreads,Parallel_threadid
      external Parallel_nthreads,Parallel_threadid

      nthr = Parallel_nthreads()
      if (nthr.eq.1) then
         call dcopy(n,a,1,b,1)
      else

         if (n.gt.nbig) then
            tid    = Parallel_threadid()
            istart = 1 + tid*(n/nthr)
            if (tid.eq.(nthr-1)) then
               nsz = n/nthr + mod(n,nthr)
            else
               nsz = n/nthr
            end if
            if (barrier) then
               call dcopy(nsz,a(istart),1,b(istart),1)
!$OMP BARRIER
            else
               call dcopy(nsz,a(istart),1,b(istart),1)
            end if

         else
            if (barrier) then
!$OMP MASTER
               call dcopy(n,a,1,b,1)
!$OMP END MASTER
!$OMP BARRIER
            else
!$OMP MASTER
               call dcopy(n,a,1,b,1)
!$OMP END MASTER
            end if
         end if
      end if
#else
      call dcopy(n,a,1,b,1)
#endif
      return
      end 



*     *************************************
*     *                                   *
*     *    Parallel_shared_vector_icopy   *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_vector_icopy(barrier,n,a,b)
      implicit none
      logical barrier
      integer n
      integer a(*),b(*)

#ifdef USE_OPENMP
      integer nbig,tid,nthr,nsz,istart
      parameter (nbig=1000)

c     **** external functions ****
      integer  Parallel_nthreads,Parallel_threadid
      external Parallel_nthreads,Parallel_threadid

      nthr = Parallel_nthreads()
      if (nthr.eq.1) then
         call icopy(n,a,1,b,1)
      else

         if (n.gt.nbig) then
            tid    = Parallel_threadid()
            istart = 1 + tid*(n/nthr)
            if (tid.eq.(nthr-1)) then
               nsz = n/nthr + mod(n,nthr)
            else
               nsz = n/nthr
            end if
            if (barrier) then
               call icopy(nsz,a(istart),1,b(istart),1)
!$OMP BARRIER
            else
               call icopy(nsz,a(istart),1,b(istart),1)
            end if

         else
            if (barrier) then
!$OMP MASTER
               call icopy(n,a,1,b,1)
!$OMP END MASTER
!$OMP BARRIER
            else
!$OMP MASTER
               call icopy(n,a,1,b,1)
!$OMP END MASTER
            end if
         end if
      end if
#else
      call icopy(n,a,1,b,1)
#endif
      return
      end




*     *************************************
*     *                                   *
*     *    Parallel_shared_vector_zero    *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_vector_zero(barrier,n,b)
      implicit none
      logical barrier
      integer n
      real*8 b(*)

#ifdef USE_OPENMP
      integer nbig,tid,nthr,nsz,istart
      parameter (nbig=1000)

c     **** external functions ****
      integer  Parallel_nthreads,Parallel_threadid
      external Parallel_nthreads,Parallel_threadid

      nthr = Parallel_nthreads()
      if (nthr.eq.1) then
         call dcopy(n,0.0d0,0,b,1)
      else

         if (n.gt.nbig) then
            tid    = Parallel_threadid()
            istart = 1 + tid*(n/nthr)
            if (tid.eq.(nthr-1)) then
               nsz = n/nthr + mod(n,nthr)
            else
               nsz = n/nthr
            end if
            if (barrier) then
               call dcopy(nsz,0.0d0,0,b(istart),1)
!$OMP BARRIER
            else
               call dcopy(nsz,0.0d0,0,b(istart),1)
            end if

         else
            if (barrier) then
!$OMP MASTER
               call dcopy(n,0.0d0,0,b,1)
!$OMP END MASTER
!$OMP BARRIER
            else
!$OMP MASTER
               call dcopy(n,0.0d0,0,b,1)
!$OMP END MASTER
            end if
         end if
      end if
#else
      call dcopy(n,0.0d0,0,b,1)
#endif
      return
      end 


*     *************************************
*     *                                   *
*     *    Parallel_shared_matrix_zero    *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_matrix_zero(barrier,n1,n2,ldc,c)
      implicit none
      logical barrier
      integer n1,n2,ldc
      real*8 c(ldc,*)

#ifdef USE_OPENMP
      integer nbig,tid,nthr,nsz,istart,i,j
      parameter (nbig=30)

c     **** external functions ****
      integer  Parallel_nthreads,Parallel_threadid
      external Parallel_nthreads,Parallel_threadid

      nthr = Parallel_nthreads()
      if (nthr.eq.1) then
         do j=1,n2
            call dcopy(n1,0.0d0,0,c(1,j),1)
         end do
      else
         if (n2.gt.nbig) then
            tid    = Parallel_threadid()
            if (barrier) then
               do j=1+tid,n2,nthr
                  call dcopy(n1,0.0d0,0,c(1,j),1)
               end do
!$OMP BARRIER
            else
               do j=1+tid,n2,nthr
                  call dcopy(n1,0.0d0,0,c(1,j),1)
               end do
            end if
         else
            if (barrier) then
!$OMP MASTER
               do j=1,n2
                  call dcopy(n1,0.0d0,0,c(1,j),1)
               end do
!$OMP END MASTER
!$OMP BARRIER
            else
!$OMP MASTER
               do j=1,n2
                  call dcopy(n1,0.0d0,0,c(1,j),1)
               end do
!$OMP END MASTER
            end if
         end if
      end if
#else
 
      integer j
      do j=1,n2
         call dcopy(n1,0.0d0,0,c(1,j),1)
      end do
#endif
      return
      end 




*     *************************************
*     *                                   *
*     *        Parallel_shared_copy       *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_copy(barrier,a,b)
      implicit none
      logical barrier
      real*8 a,b

#ifdef USE_OPENMP
c     **** external functions ****
      integer  Parallel_nthreads
      external Parallel_nthreads

      if (Parallel_nthreads().eq.1) then
         b = a
      else
         if (barrier) then
!$OMP MASTER
            b = a
!$OMP END MASTER
!$OMP BARRIER
         else
!$OMP MASTER
            b = a
!$OMP END MASTER
         end if
      end if
#else
      b = a
#endif
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel_shared_icopy      *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_icopy(barrier,a,b)
      implicit none
      logical barrier
      integer a,b

#ifdef USE_OPENMP
c     **** external functions ****
      integer  Parallel_nthreads
      external Parallel_nthreads

      if (Parallel_nthreads().eq.1) then
         b = a
      else
         if (barrier) then
!$OMP MASTER
            b = a
!$OMP END MASTER
!$OMP BARRIER
         else
!$OMP MASTER
            b = a
!$OMP END MASTER
         end if
      end if
#else
      b = a
#endif
      return
      end 

*     *************************************
*     *                                   *
*     *        Parallel_shared_iinc       *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_iinc(barrier,i,a)
      implicit none
      logical barrier
      integer i,a

#ifdef USE_OPENMP
c     **** external functions ****
      integer  Parallel_nthreads
      external Parallel_nthreads


      if (Parallel_nthreads().eq.1) then
         a = a + i
      else
         if (barrier) then
!$OMP MASTER
            a = a + i
!$OMP END MASTER
!$OMP BARRIER
         else
!$OMP MASTER
            a = a + i
!$OMP END MASTER
         end if
      end if
#else
      a = a + i
#endif
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel_shared_irot       *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_irot(barrier,imax,a)
      implicit none
      logical barrier
      integer imax,a

#ifdef USE_OPENMP
c     **** external functions ****
      integer  Parallel_nthreads
      external Parallel_nthreads

      !write(*,*) "irot nthreads = ",Parallel_nthreads()

      if (Parallel_nthreads().eq.1) then
         a = a + 1
         if (a.gt.imax) a = 1
      else
         if (barrier) then
!$OMP MASTER
            a = a + 1
            if (a.gt.imax) a = 1
!$OMP END MASTER
!$OMP BARRIER
         else
!$OMP MASTER
            a = a + 1
            if (a.gt.imax) a = 1
!$OMP END MASTER
         end if
      end if
#else
      a = a + 1
      if (a.gt.imax) a = 1
#endif
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel_shared_lshift     *
*     *                                   *
*     *************************************
      subroutine Parallel_shared_lshift(barrier,n,a)
      implicit none
      logical barrier
      integer n,a(*)
      integer i

#ifdef USE_OPENMP
c     **** external functions ****
      integer  Parallel_nthreads
      external Parallel_nthreads

      if (Parallel_nthreads().eq.1) then
         do i=1,n
           a(i) = a(i+1)
         end do
      else
         if (barrier) then
!$OMP MASTER
            do i=1,n
              a(i) = a(i+1)
            end do
!$OMP END MASTER
!$OMP BARRIER
         else
!$OMP MASTER
            do i=1,n
              a(i) = a(i+1)
            end do
!$OMP END MASTER
         end if
      end if
#else
      do i=1,n
        a(i) = a(i+1)
      end do
#endif
      return
      end 




