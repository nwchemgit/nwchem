      Subroutine grid_setquad(rtdb)
c
c$Id: grid_setquad.F,v 1.3 1999-06-14 22:39:23 d3h449 Exp $
c
      implicit none
      logical oprint
c      
#include "rtdb.fh"
#include "cdft.fh"
#include "util.fh"
#include "global.fh"
c
      integer rtdb
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
c
c     Set up needed info by atom type (ntypes info; e.g., nuclear charge,
c     Bragg-Slater radii, etc)
c
      call grid_atom_type_info
      call ga_sync()
c
c     Determine desired integration accuracy and set 
c     corresponding quadrature definition.
c
      call grid_acc_def(rtdb)
      call ga_sync()
c
c     Define the radial quadrature points and weights.
c
      call grid_setrad(rtdb)
      call ga_sync()
c
c     Define the angular quadrature points and weights.
c
      call grid_setang(rtdb)
      call ga_sync()
c
c     Define the rotation matrix that will be used to ensure that 
c     the quadrature grid is rotationally invariant.
c
      call grid_invariant
      call ga_sync()
c
c     Compute/set any spatial weight parameters.
c
      call grid_setspac_params
      call ga_sync()
c
c     Generate the complete grid and decompose it into blocks
c     which can define the computational task size.
c  
      call grid_gen_pts(rtdb)
      call ga_sync()
c
c     Update any req'd data in rtdb
c
      call grid_update_rtdb(rtdb)
      call ga_sync()
c
c     Setup reduced matrices (eval gaussian on grid pts and reduce to atom centers) 
c
      call grid_rchi(rtdb)
      call ga_sync()
c
c     Memory management (as before)
c
      call grid_memory(rtdb)
      call ga_sync()
c
c     Try to see if you can read the quadrature file.
c
c      call grid_read_ck(rtdb)
c
c     Screening parameters are generated here.
c
c      call grid_prescreen_1
c      call grid_prescreen_2
c
c     Memory management must be determined either here or 
c     above the prescreening.
c
c      call grid_memory
      return
      end
