c $Id$
      subroutine util_test(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer rtdb

      if (ga_nodeid() .eq. 0) then
        call util_print_centered
     1  (6,'Global Array & ParI/O Tests',40,.true.)
        write(6,*)
      endif

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing GA Basic Test"
      write(6,*) "------------------------"
      write(6,*) 
      endif
      
      call util_ga_test
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing GA Patch Test"
      write(6,*) "------------------------"
      write(6,*) 
      endif
      
      call util_patch_test
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing NGA Test"
      write(6,*) "-------------------"
      write(6,*) 
      endif
      
      call util_ndim_test
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing GA Matrix Algebra Test"
      write(6,*) "---------------------------------"
      write(6,*) 
      endif
      
      call util_test_eig
      call ga_sync()
      call util_test_lu
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing DRA Test"
      write(6,*) "-------------------"
      write(6,*) 
      endif
      
      call util_dra_test
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing SF Test"
      write(6,*) "------------------"
      write(6,*) 
      endif
      
      call util_sf_test
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing EAF Test"
      write(6,*) "-------------------"
      write(6,*) 
      endif
      
      call util_eaf_test
      call ga_sync()

      if (ga_nodeid() .eq. 0) then
      write(6,*) 
      write(6,*) "Commencing GA Performance Test"
      write(6,*) "------------------------------"
      write(6,*) 
      endif
      
      call util_perf_test
      call ga_sync()

      return
      end
