C> \ingroup stpr_priv
C> @{
C>
C> \brief Initialize the file names
C>
C> The filenames used by the stepper module are stored in a common
C> block. This subroutine initializes the filename variables.
C>
      subroutine stpr_filenames
*
* $Id$
*
      implicit none
#include "cstprfiles.fh"
#include "util.fh"
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
c
c     initialize the common with filenames
c
      call util_file_name('stpr41',.false., .false.,FILE41)
      call util_file_name('stpr44',.true.,  .false.,FILE44)
      call util_file_name('stpr45',.true.,  .false.,FILE45)
      call util_file_name('hess',  .false., .false.,FILEHESS)
      call util_file_name('fdrst', .false., .false.,FILEATR)
      call util_file_name('fd_ddipole',.false., .false.,FILEDDIPOLE)
c
      if (util_print("files",print_debug).and.ga_nodeid().eq.0) then
        write(LuOut,*)"STEPPER filenames"
        write(LuOut,*)"================="
        write(LuOut,*)"file41  : ",FILE41(1:inp_strlen(FILE41))
        write(LuOut,*)"file44  : ",FILE44(1:inp_strlen(FILE44))
        write(LuOut,*)"file45  : ",FILE45(1:inp_strlen(FILE45))
        write(LuOut,*)"hessian : ",FILEHESS(1:inp_strlen(FILEHESS))
        write(LuOut,*)"restart : ",FILEATR(1:inp_strlen(FILEATR))
        write(LuOut,*)"ddipole : ",
     &                FILEDDIPOLE(1:inp_strlen(FILEDDIPOLE))
      endif
      end
C>
C> @}
      
      
