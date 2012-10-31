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
#include "errquit.fh"

*     **** local variables ****
      logical value
      character*80 filename
      character*255 full_filename,full_filename2

*     **** external functions ****
      logical  control_Fei
      external control_Fei

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

      return
      end

