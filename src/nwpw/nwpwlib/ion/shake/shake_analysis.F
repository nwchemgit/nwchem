*
* $Id: shake_analysis.F,v 1.3 2006-12-28 18:29:26 bylaska Exp $
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


      call molecular_init(rtdb)
      value = rtdb_parallel(.false.)
      if (taskid.eq.MASTER) then

         write(luout,*)
         write(luout,*) '===== Performing shake analysis ===='
         write(luout,*)


         if (.not.rtdb_cget(rtdb,'cpmd:fei_filename',1,filename))
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
      value = rtdb_parallel(.true.)
      call molecular_end()
   

      return
      end



