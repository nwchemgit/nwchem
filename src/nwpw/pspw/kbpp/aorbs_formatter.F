c
c $Id$
c

*     ***********************************
*     *					*
*     *		aorbs_formatter_auto	*
*     *	  				*
*     ***********************************
      subroutine aorbs_formatter_auto(atom,rcut,lmbda)
      implicit none
      character*4  atom
      real*8 rcut,lmbda

#include "errquit.fh"

*     ***** local variables *****
      logical       value
      integer       i,j,l
      integer	    locp,lmax
      character*50  psp_filename,formatted_filename
      character*4   element
      real*8        rlocal

*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version

*     ***** external functions *****
      integer      control_version
      integer      control_ngrid
      real*8       control_unita
      logical      kbpp_orb, nwpw_filefind
      external     control_version
      external     control_ngrid
      external     control_unita
      external     kbpp_orb, nwpw_filefind

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
      formatted_filename = element(1:l)//'.aorb'


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

      value = value.and.
     >     kbpp_orb(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,rcut,lmbda)
*     *** exit if error found ***
      if(.not. value) call errquit('aorb_formatter_auto: error',0,
     &       CALC_ERR)

      return
      end


