*
* $Id$
*

!  **************************************************
!  
!       name:
!  
!       purpose
!  
!       created:        7/7/2002
!  **************************************************
      logical function my_alloc (v_type, v_size, v_name, v)
 
 
#include "bafdecls.fh"
 
      character v_name*(*)
      integer v_type, v_size
      integer v(4)
 
!    *** local variables ***
 
      v(4) = v_type
      v(3) = v_size
      my_alloc = ma_alloc_get(v(4),v(3),v_name,v(2),v(1))
 
      end 
 
!  **************************************************
!  
!       name:
!  
!       purpose
!  
!       created:        7/7/2002
!  **************************************************
      logical function my_dealloc (v)
 
#include "bafdecls.fh" 

      integer v(4)
 
      my_dealloc = BA_free_heap(v(2))
 
 
      end 
 
      subroutine paw_ma_print(v)

      implicit none
      
      integer v(4)


      return
      end
