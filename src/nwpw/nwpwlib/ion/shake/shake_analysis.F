*
* $Id: shake_analysis.F,v 1.2 2006-12-28 01:56:52 bylaska Exp $
*

*     ***************************
*     *			        *
*     *	     shake_analysis     *
*     *				*
*     ***************************
      subroutine shake_analysis(rtdb)
      implicit none 
      integer   rtdb

#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value
      character*30  filename
      character*255 full_filename

*     **** external functions ****

      call ga_sync()

      call Parallel_taskid(taskid)


      value = rtdb_parallel(.false.)
      if (taskid.eq.MASTER) then

         write(luout,*)
         write(luout,*) '===== Performing shake analysis ===='
         write(luout,*)

         call molecular_init(rtdb)

         if (.not.rtdb_cget(rtdb,'cpmd:fei_filename',1,filename))
     >     call util_file_prefix('fei',filename)

         call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

         call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

         call molecular_force_average(full_filename)
         

         call molecular_end()

      end if
      call ga_sync()
      value = rtdb_parallel(.true.)
   

      return
      end



