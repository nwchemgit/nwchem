      logical function spcart_init(lmaxin)
* 
* initialization of spherical to cartesean tranformation array
*... stored on heap.  
*... stored up to lmax values
*...
*   spcart(iccart,icsp,l)  => spcart((lmax+1)*(lmax+2)/2,1:2*lmax+1,0:lmax)
*    lmax = 5 h functions  => size = 21*11*6 = 1386 
*   store array 34% sparse for simplicity.  1386 doubles is 11 Kbytes.
*
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
c::passed
      integer lmaxin   ! [input] init transformed values up to lmaxin
c::local
      integer size_sp2c         ! size of array
c
      size_sp2c = lmaxin+1
      size_sp2c = size_sp2c*(2*lmaxin+1)
      size_sp2c = size_sp2c*(((lmaxin+1)*(lmaxin+2))/2)
c
      if (sphcart_init.eq.SPHCART_INIT_VALUE) then
        if (lmaxin.gt.lmax_init) then
          if (.not.ma_free_heap(h_sp2c)) call errquit
     &        ('spcart_init: freeing old array failed',911)
          k_sp2c = 0
          h_sp2c = 0
        else
          spcart_init = .true.
          return   ! initialization already done to cover lmaxin
        endif
      endif
*
      spcart_init =
     &    ma_alloc_get(MT_DBL,size_sp2c,'sph 2 cart trans array',
     &    h_sp2c,k_sp2c)
      if (.not.spcart_init)  call errquit
     &    ('spcart_init: ma_alloc_get failed for size',size_sp2c)
c
      call xlmcoeff(lmaxin,dbl_mb(k_sp2c),.true.)  ! always normalize 
c
      sphcart_init = SPHCART_INIT_VALUE
      lmax_init = lmaxin
c
      end
      logical function spcart_terminate()
      implicit none
c
c terminates spcart data structure and initialization
c
      if (sphcart_init.eq.SPHCART_INIT_VALUE) then
        spcart_terminate = ma_free_heap(h_sp2c)
        k_sp2c = 0
        h_sp2c = 0
        sphcart_init = 0
      else
        spcart_terminate = .false.
      endif
      end

