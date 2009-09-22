c
c     $Id: cpmd_properties.F,v 1.7 2007-04-02 15:38:38 bylaska Exp $
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
      character*30 filename
      character*255 full_filename

*     **** external functions ****
      logical  control_Fei
      external control_Fei

      if (.not.rtdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename))
     >  call util_file_prefix('ion_motion',filename)

      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)


      call pair_dist(full_filename)
      call pair_distab(full_filename)
      call velocity_autocorrelation(full_filename)


*     **** perform shake force averaging if Fei set *****
      if (control_Fei()) then
         call shake_analysis(rtdb)
         call FixIon_analysis(rtdb)
      end if

      return
      end

