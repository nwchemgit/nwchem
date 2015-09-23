      Subroutine grid_setquad(rtdb)
c
c$Id$
c
      implicit none
      logical oprint, ogrdebug
c      
#include "rtdb.fh"
#include "cdft.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c
      integer rtdb
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      ogrdebug = util_print('griddebug', print_high)
c
c     Set up needed info by atom type (ntypes info; e.g., nuclear charge,
c     Bragg-Slater radii, etc)
c
      call grid_atom_type_info
cc AJL/Begin/FDE
      if (frozemb_fde) call grid_atom_type_info_fde
cc AJL/End
c
c     Determine desired integration accuracy and set 
c     corresponding quadrature definition.
c
      call grid_acc_def(rtdb)
c
c     Define the radial quadrature points and weights.
c
      call grid_setrad(rtdb)
c
c     Define the angular quadrature points and weights.
c
      call grid_setang(rtdb)
c
c     Define the rotation matrix that will be used to ensure that 
c     the quadrature grid is rotationally invariant.
c
      call grid_invariant
c
c     Compute/set any spatial weight parameters.
c
      call grid_setspac_params
c
c     Generate the complete grid and decompose it into blocks
c     which can define a computational task size on a file system (if desired).
c  
      call grid_gen_pts(rtdb)
      if(ogrdebug) write(LuOut,*) ' me ',ga_nodeid(),'  done genpts'
c
c     Update any req'd data in rtdb
c
      call grid_update_rtdb(rtdb)
      if(ogrdebug) write(LuOut,*) ' me ',ga_nodeid(),'  done uprt'
c
c     Memory management (as before)
c
      call grid_memory(rtdb)
      if(ogrdebug) write(LuOut,*) ' me ',ga_nodeid(),' done mem'
      call ga_sync()
      return
      end
