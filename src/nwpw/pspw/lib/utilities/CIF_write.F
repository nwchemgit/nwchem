*
* $Id$
*

*     *********************************
*     *                               *
*     *          CIF_write            *
*     *                               *
*     *********************************
      subroutine CIF_write(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"

      integer taskid,MASTER
      parameter (MASTER=0)

      integer       flen
      character*50 filename
      character*255 full_filename

      if (btdb_cget(rtdb,'pspw:cif_filename',1,filename)) then
         flen  = index(filename,' ') - 1
         filename = filename(1:flen)//'.cif'
     
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)

         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
            write(6,*) 
            write(6,*) 'Writing Crystallographic Information File:', 
     >                 full_filename
            open(unit=19,file=full_filename,form='formatted')
            call ion_Print_CIF(19)
            close(19)
         end if
      end if


      return
      end

