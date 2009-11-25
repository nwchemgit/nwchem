*
* $Id: Parallel.F,v 1.4 2008-04-21 22:20:45 bylaska Exp $
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
#include "mafdecls.fh"
#include "errquit.fh"

      if (.not.MA_free_heap(procNd(2)))
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
#include "mafdecls.fh"
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
#include "mafdecls.fh"
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
#include "mafdecls.fh"
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
#include "mafdecls.fh"
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
#include "mafdecls.fh"
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
#include "mafdecls.fh"
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
#include "mafdecls.fh"
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




