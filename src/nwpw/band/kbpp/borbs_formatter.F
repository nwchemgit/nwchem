c
c $Id$
c

*     ***************************************
*     *                                     *
*     *		borbs_formatter_auto        *
*     *                                     *
*     ***************************************
      subroutine borbs_formatter_auto(atom,rcut,lmbda)
      implicit none
      character*4  atom
      real*8 rcut,lmbda


#include "bafdecls.fh"
#include "errquit.fh"

*     ***** local variables *****
      logical       value
      integer       i,j,l,nb
      integer	    locp,lmax,kv(2),nbrillioun
      character*50  psp_filename,formatted_filename
      character*4   element
      real*8        rlocal

*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version

*     ***** external functions *****
      integer      control_version,brillioun_nbrillioun
      integer      control_ngrid
      real*8       control_unita,brillioun_all_k
      logical      kbpp_band_orb, nwpw_filefind
      external     control_version,brillioun_nbrillioun
      external     control_ngrid
      external     control_unita,brillioun_all_k
      external     kbpp_band_orb, nwpw_filefind

      value  = .true.

*     **** get lmax,locp,and rlocal from psp library ****
      lmax  = -1
      locp  = -1
      rlocal = 1.0d0

*     **** define psp_filename and formatted_filename ****
      element = '    '
      element = atom
      l = index(element,' ') - 1
      psp_filename       = element(1:l)//'.psp1'
      formatted_filename = element(1:l)//'.borb'


*     **** generate one-dimensional pseudopotential file ****
      if (.not.nwpw_filefind(psp_filename)) then
        call psp_generator_auto1(atom,lmax,locp,rlocal)

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

      value = value.and.
     >     kbpp_band_orb(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,
     >            nbrillioun,dbl_mb(kv(1)),rcut,lmbda)

*     *** exit if error found ***
      value = value.and.BA_pop_stack(kv(2))
      if(.not. value) call errquit('borb_formatter_auto: error',0,
     >       CALC_ERR)

      return
      end


