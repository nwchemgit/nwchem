      program test
*
* $Id$
*
      implicit none
#include "mafdecls.fh"
      double precision a(2,3)
      integer b(3,2)
      logical c(5,2)
      character*8 d(4)
      integer old, junk
c
      data a/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0/
      data b/6, 5, 4, 3, 2, 1/
      data c/.false., .false., .false., .false., .false., 
     $     .true., .true., .true., .true., .true./
      data d/'one', 'two', 'three', 'four'/
c
      call ecce_print_file_open('/tmp/ecce_output')
      call ecce_print_module_entry('Module A')
      call ecce_print1('A vector of 6 doubles', mt_dbl, a, 6)
      call ecce_print_module_entry('Module B')
      call ecce_print2('A 2x3 matrix of doubles', mt_dbl, a, 2, 2, 3)
      call ecce_print2('A 3x2 array of integers', mt_int, b, 3, 3, 2)
*      call ecce_print_module_exit('Not there', 'Oh dear')
      call ecce_print_module_exit('Module B', 'OK')
c
      call ecce_print_control(0, old)
      call ecce_print1('This should not be printed', mt_dbl, a, 6)
      call ecce_print_control(old, junk)
c
      call ecce_print1_char('A 4 element array of char*8', d, 4)
      call ecce_print2('A 5x2 array of logicals', mt_log, c, 5, 5, 2)
      call ecce_print_module_exit('Module A', 'Not so OK')
c
      call ecce_print_file_close()
c
      end

