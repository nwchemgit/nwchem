      logical function esp(irtdb0)
c
c $Id: esp.F,v 1.1 1997-10-06 23:00:30 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
c
      integer irtdb0
c
c     me : number current node
c     np : number of nodes
c
      me=ga_nodeid()
      np=ga_nnodes()
      irtdb=irtdb0
c
      if(me.eq.0) then
      write(6,*)
      write(6,*)
      call util_print_centered(6,
     + 'NWChem Electrostatic Potential Fit Module',40,.true.)
      write(6,*)
      write(6,*)
      call ecce_print_module_entry('nwesp')
      endif
c
c     read input from rtdb
c
      call esp_rtdbin
c
c     allocate memory
c
      call esp_malloc
c
c     setup the esp grid
c
      call esp_grid(dbl_mb(ixptr),dbl_mb(iqptr),dbl_mb(igptr))
c
c     write input variables to output     
c
      call esp_wrtinp
c
c     calculate the electrostatic potential on the grid
c
      call esp_esp(dbl_mb(ieptr),dbl_mb(ioptr),dbl_mb(idptr),
     + dbl_mb(ipptr),dbl_mb(ixptr),dbl_mb(iqptr),dbl_mb(igptr),
     + dbl_mb(ivptr))
c
c     do the fitting
c
      call esp_fit(dbl_mb(igptr),dbl_mb(ivptr),dbl_mb(ixptr),
     + dbl_mb(iqptr),dbl_mb(iaptr),dbl_mb(ibptr),dbl_mb(icptr),
     + dbl_mb(iakptr),dbl_mb(iqkptr),int_mb(inptr))
c
c     calculate the electrical moments
c
      call esp_elemom(dbl_mb(ixptr),dbl_mb(iqptr))
c
c     write results
c
      call esp_out(byte_mb(itptr),dbl_mb(ixptr),dbl_mb(iqptr))
c
c     deallocate memory
c
      call esp_dalloc
c
      esp=.true.
      return
      end
