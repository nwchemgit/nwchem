c
c     $Id$
c

*     *********************************
*     *                               *
*     *          cpmd_properties      *
*     *                               *
*     *********************************
      subroutine cpmd_properties(rtdb)
      implicit none
      integer rtdb

#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,doit,oprint
      character*80 filename
      character*255 full_filename,full_filename2

*     **** external functions ****
      logical  control_Fei
      external control_Fei

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      if (.not.rtdb_get(rtdb,'nwpw:cpmd_properties',mt_log,1,doit))
     >   doit = .false.

      if (oprint) then
         write(luout,*) 
         write(luout,*) "================================"
         write(luout,*) "==== CPMD Property Analysis ===="
         write(luout,*) "================================"
         if (doit) then
            write(luout,*) 
            write(luout,*) "CPMD property analysis is on. ",
     >         "To turn analysis off set cpmd_properties, i.e."
            write(luout,*) "nwpw"
            write(luout,*) "   cpmd_properties off"
            write(luout,*) "end"
            write(luout,*)
         else
            write(luout,*)
            write(luout,*) "CPMD property analysis is off. ",
     >         "To turn analysis on set cpmd_properties, i.e."
            write(luout,*) "nwpw"
            write(luout,*) "   cpmd_properties on"
            write(luout,*) "end"
            write(luout,*)
         end if
      end if

      if (doit) then

      if (.not.rtdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename))
     >    call util_file_prefix('ion_motion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

      if (.not.rtdb_cget(rtdb,'cpmd:unfolded_xyz_filename',1,filename))
     >    call util_file_prefix('unfolded.xyz',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename2)

      call pair_dist(full_filename)
      call pair_distab(full_filename)
      call velocity_autocorrelation(full_filename)
      call generate_unfolded_xyz(full_filename,full_filename2)


*     **** perform shake force averaging if Fei set *****
      if (control_Fei()) then
         call shake_analysis(rtdb)
         call FixIon_analysis(rtdb)
      end if

      end if

      return
      end

