c     $Id: cpmd_properties.F,v 1.2 2001-10-26 20:13:25 edo Exp $
*     *********************************
*     *                               *
*     *          cpmd_properties      *
*     *                               *
*     *********************************
      subroutine cpmd_properties(rtdb)
      implicit none
      integer rtdb

#include "rtdb.fh"

      logical value
      character*30 filename
      character*255 full_filename

      value = rtdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename)
      if (.not. value)
     >    call errquit('cpmd_properties: cant get MOTION',0)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)



      call pair_dist(full_filename)
      call velocity_autocorrelation(full_filename)
      return
      end

