
*     ***************************************
*     *										*
*     *			cpsp_formatter_auto			*
*     *										*
*     ***************************************
      subroutine cpsp_formatter_auto(atom)
      implicit none
      character*2  atom

#include "mafdecls.fh"

*     ***** local variables *****
      logical       value
      integer       i,j,l
      integer	    locp,lmax
      character*20  psp_filename,formatted_filename
      character*4   element
      real*8        rlocal
      integer       nbrillioun,kv(2),nb

*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version

*     ***** external functions *****
      integer      control_version
      integer      control_ngrid,brillioun_nbrillioun
      real*8       control_unita,brillioun_k
      logical      kbpp_band, pspw_filefind
      external     control_version
      external     control_ngrid,brillioun_nbrillioun
      external     control_unita,brillioun_k
      external     kbpp_band, pspw_filefind

      value  = .true.

*     **** get lmax,locp,and rlocal from psp library ****
      lmax  = -1
      locp  = -1
      rlocal = 1.0d0

*     **** define psp_filename and formatted_filename ****
      element = '    '
      element = atom
      l = index(element,' ') - 1
      psp_filename       = element(1:l)//'.psp'
      formatted_filename = element(1:l)//'.cpp'

*     **** generate one-dimensional pseudopotential file ****
      if (.not.pspw_filefind(psp_filename)) then
        call psp_generator_auto(atom,lmax,locp,rlocal)

*     **** get lmax,locp,and rlocal from rtdb ****
      else
         call control_pspparameters(atom,lmax,locp,rlocal)
      end if


*     **** cell information ****
      version = control_version()
      do i=1,3
         ngrid(i) = control_ngrid(i)
         do j=1,3
           unita(i,j) = control_unita(i,j)
         end do
      end do

*     **** brillioun zone information ****
      nbrillioun = brillioun_nbrillioun()
      value = MA_push_get(mt_dbl,3*nbrillioun,'kv',kv(2),kv(1))
      do nb=1,nbrillioun
         dbl_mb(kv(1)+3*(nb-1))   = brillioun_k(1,nb)
         dbl_mb(kv(1)+3*(nb-1)+1) = brillioun_k(2,nb)
         dbl_mb(kv(1)+3*(nb-1)+2) = brillioun_k(3,nb)
      end do


      value = value.and.
     >     kbpp_band(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,
     >            nbrillioun,dbl_mb(kv(1)))

      value = value.and.MA_pop_stack(kv(2))
*     *** exit if error found ***
      if(.not. value) call errquit('psp_formatter_auto: error',0)

      return
      end


