      Subroutine grid_setquad(rtdb)
c
c$Id: grid_setquad.F,v 1.7 2000-03-20 21:11:49 edo Exp $
c
      implicit none
      logical oprint, grid_on_file,ogrdebug
c      
#include "rtdb.fh"
#include "cdft.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
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
      if(ogrdebug) write(6,*) ' me ',ga_nodeid(),'  done genpts'
c
c     Update any req'd data in rtdb
c
      call grid_update_rtdb(rtdb)
      if(ogrdebug) write(6,*) ' me ',ga_nodeid(),'  done uprt'
c
c     Setup reduced matrices (eval gaussian on grid pts and reduce to 
c     atom centers) - only if points are on file system
c
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &   grid_on_file))then
         grid_on_file = .false.
      endif
cedo rchi stuff moved in grid_quadv
      call grid_rchi0
c
c     Memory management (as before)
c
      call grid_memory(rtdb)
      if(ogrdebug) write(6,*) ' me ',ga_nodeid(),' done mem'
      call ga_sync()
      return
      end
