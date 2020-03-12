c
c $Id$
c

*     ***************************************
*     *				  	    *
*     *	     cpsp_stress_formatter_auto     *
*     *					    *
*     ***************************************

      subroutine cpsp_stress_formatter_auto(atom)
      implicit none
      character*4  atom

#include "bafdecls.fh"
#include "errquit.fh"


*     ***** local variables *****
      logical       value
      integer       i,j,l
      integer       locp,lmax
      character*50  psp_filename,formatted_filename
      character*5   element
      real*8        rlocal
      integer       nbrillioun,kv(2),nb


*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version,psp_type

*     ***** external functions *****
      integer      control_version
      integer      control_ngrid,brillioun_nbrillioun
      real*8       control_unita,brillioun_all_k
      logical      kbpp_band_stress, nwpw_filefind,kbpp_e_band_stress
      external     control_version
      external     control_ngrid,brillioun_nbrillioun
      external     control_unita,brillioun_all_k
      external     kbpp_band_stress, nwpw_filefind,kbpp_e_band_stress

      value  = .true.

*     **** get lmax,locp,and rlocal from psp library ****
      lmax  = -1
      locp  = -1
      rlocal = 1.0d0

*     **** define psp_filename and formatted_filename ****
      element = '     '
      element = atom
      l = index(element,' ') - 1
      psp_filename       = element(1:l)//'.psp'
      formatted_filename = element(1:l)//'.cpp2'

*     **** generate one-dimensional pseudopotential file ****
      if (.not.nwpw_filefind(psp_filename)) then
         call psp_generator_auto(1,atom,lmax,locp,rlocal)

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
      value = BA_push_get(mt_dbl,3*nbrillioun,'kv',kv(2),kv(1))
      do nb=1,nbrillioun
         dbl_mb(kv(1)+3*(nb-1))   = brillioun_all_k(1,nb)
         dbl_mb(kv(1)+3*(nb-1)+1) = brillioun_all_k(2,nb)
         dbl_mb(kv(1)+3*(nb-1)+2) = brillioun_all_k(3,nb)
      end do

      call get_psp_type(psp_filename,psp_type)

      if (psp_type.eq.0) then
      value = value.and.
     >       kbpp_band_stress(.false.,
     >              version,psp_filename,formatted_filename,ngrid,unita,
     >              locp,lmax,
     >              nbrillioun,dbl_mb(kv(1)))
      else if (psp_type.eq.1) then
        call errquit('cpsp_stress_formatter_auto: unknown psp_type',
     >              psp_type,INPUT_ERR)

      else if (psp_type.eq.2) then
      value = value.and.
     >       kbpp_e_band_stress(.false.,
     >              version,psp_filename,formatted_filename,ngrid,unita,
     >              locp,lmax,
     >              nbrillioun,dbl_mb(kv(1)))

      else if (psp_type.eq.7) then
        call errquit('cpsp_stress_formatter_auto: unknown psp_type',
     >              psp_type,INPUT_ERR)
      end if

      value = value.and.BA_pop_stack(kv(2))


*     *** exit if error found ***
      if(.not. value) call errquit('cpsp_stress_formatter_auto: error',
     >                             0,RTDB_ERR)

      return
      end


