      subroutine stpr_filenames
*
* $Id: stpr_filenames.F,v 1.3 1997-10-31 20:50:35 d3e129 Exp $
*
      implicit none
#include "cstprfiles.fh"
c
c     initialize the common with filenames
c
      call util_file_name('stpr41',.false., .false.,FILE41)
      call util_file_name('stpr44',.true.,  .false.,FILE44)
      call util_file_name('stpr45',.true.,  .false.,FILE45)
      call util_file_name('hess',  .false., .false.,FILEHESS)
      call util_file_name('fdrst', .false.,.false.,FILEATR)
      call util_file_name('fd_ddipole',.false., .false.,FILEDDIPOLE)
c
      end

      
      
