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

      integer function Parallel_threadid()
      implicit none
#include "Parallel.fh"

      Parallel_threadid = threadid
      return
      end

      integer function Parallel_nthreads()
      implicit none
#include "Parallel.fh"

      Parallel_nthreads = nthreads
      return
      end



