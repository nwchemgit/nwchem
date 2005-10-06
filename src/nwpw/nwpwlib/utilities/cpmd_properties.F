c
c     $Id: cpmd_properties.F,v 1.5 2005-10-06 01:04:34 bylaska Exp $
c

*     *********************************
*     *                               *
*     *          cpmd_properties      *
*     *                               *
*     *********************************
      subroutine cpmd_properties(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "rtdb.fh"

      logical value
      character*30 filename
      character*255 full_filename

      if (.not.rtdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename))
     >  call util_file_prefix('ion_motion',filename)

      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)



      call pair_dist(full_filename)
      call velocity_autocorrelation(full_filename)
      return
      end

