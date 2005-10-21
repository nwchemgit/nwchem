      logical function esp(irtdb0)
c
c $Id: esp.F,v 1.12 2005-10-21 19:22:39 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
c
      integer irtdb0
c
      call util_print_push
      call util_print_rtdb_load(irtdb,'esp')
c
c     me : number current node
c     np : number of nodes
c
      me=ga_nodeid()
      np=ga_nnodes()
      irtdb=irtdb0
      lfnout=6
      lfngrd=24
      lfner=25
      lfnplt=26
      lfnxyz=27
      lfnchg=28
c
c     read input from rtdb
c
      call esp_rtdbin
c
      if(util_print('output',print_medium)) then
      if(me.eq.0) then
      if(util_print('information',print_default)) then
      write(lfnout,*)
      write(lfnout,*)
      call util_print_centered(6,
     + 'NWChem Electrostatic Potential Fit Module',40,.true.)
      write(lfnout,*)
      write(lfnout,*)
      endif
      call ecce_print_module_entry('nwesp')
      endif
      endif
c
c     allocate memory
c
      call esp_malloc
c
c     setup the esp grid
c
      call esp_grid(dbl_mb(ixptr),dbl_mb(iqptr),dbl_mb(igptr),
     + dbl_mb(ivptr))
c
c     write input variables to output     
c
      call esp_wrtinp(dbl_mb(iqptr))
c
c     calculate the electrostatic potential on the grid
c
      if(newgrd.ne.0) then
      call esp_esp(dbl_mb(ieptr),dbl_mb(ioptr),dbl_mb(idptr),
     + dbl_mb(ixptr),dbl_mb(iqptr),dbl_mb(igptr),
     + dbl_mb(ivptr))
      endif
c
c     do the fitting
c
      call esp_fit(dbl_mb(igptr),dbl_mb(ivptr),dbl_mb(ixptr),
     + dbl_mb(iqptr),dbl_mb(iaptr),dbl_mb(ibptr),dbl_mb(icptr),
     + dbl_mb(iakptr),dbl_mb(iqkptr),int_mb(inptr),dbl_mb(iqfptr),
     + int_mb(ilptr),byte_mb(itptr))
c
c     calculate the electrical moments
c
      call esp_elemom(dbl_mb(ixptr),dbl_mb(iqfptr))
c
c     calculate rms deviation
c
      call esp_rms(dbl_mb(ixptr),dbl_mb(iqfptr),dbl_mb(igptr),
     + dbl_mb(ivptr))
c
c     write results
c
      call esp_out(byte_mb(itptr),dbl_mb(ixptr),dbl_mb(iqfptr),
     + dbl_mb(igptr),dbl_mb(ivptr))
c
c     electrostatic potential from partial charges on full grid
c
      call esp_plt(dbl_mb(ixptr),dbl_mb(iqfptr),dbl_mb(igptr),
     + byte_mb(itptr))
c
c     deallocate memory
c
      call esp_dalloc
c
      call util_print_pop
c
      esp=.true.
      return
      end
