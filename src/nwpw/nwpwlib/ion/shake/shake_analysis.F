*
* $Id$
*

*     ***************************
*     *			        *
*     *	     shake_analysis     *
*     *				*
*     ***************************
      subroutine shake_analysis(rtdb)
      implicit none 
      integer   rtdb

#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value
      character*50  filename
      character*255 full_filename

*     **** external functions ****

      call ga_sync()

      call Parallel_taskid(taskid)


      call molecular_init(rtdb)
      value = btdb_parallel(.false.)
      if (taskid.eq.MASTER) then

         write(luout,*)
         write(luout,*) '===== Performing shake analysis ===='
         write(luout,*)


         if (.not.btdb_cget(rtdb,'cpmd:fei_filename',1,filename))
     >     call util_file_prefix('fei',filename)

         call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

         call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

         call molecular_force_average(full_filename)
         


      end if
      call ga_sync()
      value = btdb_parallel(.true.)
      call molecular_end()
   

      return
      end



