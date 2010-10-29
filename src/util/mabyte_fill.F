      subroutine mabyte_fill(n,char2use,k_ma_array,inc_array)
*  $Id$
      implicit none
c
#include "mafdecls.fh"
c
      integer n            ! [input] number of bytes in ma_byte array
      character*1 char2use ! [input] fill character
      integer k_ma_array   ! [input] ma index if byte array in byte_mb
      integer inc_array    ! [input] increment in array
c
      integer i
c
      if (inc_array.eq.1) then
        do i = 0,(n-1)
          byte_mb(k_ma_array+i) = char2use
        enddo
      else
        do i = 0,((n-1)*inc_array),inc_array
          byte_mb(k_ma_array+i) = char2use
        enddo
      endif
      end
