
*     ***************************************
*     *					    *
*     *		cpsp_formatter_auto	    *
*     *					    *
*     ***************************************
      subroutine cpsp_formatter_auto(atom)
      implicit none
      character*4  atom

#include "mafdecls.fh"
#include "errquit.fh"

*     ***** local variables *****
      logical       value
      integer       i,j,l
      integer	    locp,lmax,psp_type
      character*20  psp_filename,formatted_filename
      character*5   element
      real*8        rlocal
      integer       nbrillioun,kv(2),nb

*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version

*     ***** external functions *****
      integer      control_version
      integer      control_ngrid,brillioun_nbrillioun
      real*8       control_unita,brillioun_k
      logical      kbpp_band,nwpw_filefind,hghpp_band,rkbpp_band
      external     control_version
      external     control_ngrid,brillioun_nbrillioun
      external     control_unita,brillioun_k
      external     kbpp_band,nwpw_filefind,hghpp_band
      external     rkbpp_band
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
      formatted_filename = element(1:l)//'.cpp'

      write(6,*)"CAINE A"
      call flush(6)
*     **** generate one-dimensional pseudopotential file ****
      if (.not.nwpw_filefind(psp_filename)) then
        write(6,*)"CAINE bxxx"
        call flush(6)
        call psp_generator_auto(1,atom,lmax,locp,rlocal)
        write(6,*)"CAINE bxx"
        call flush(6)

*     **** get lmax,locp,and rlocal from rtdb ****
      else
         call control_pspparameters(atom,lmax,locp,rlocal)
        call psp_generator_auto(1,atom,lmax,locp,rlocal)
        write(6,*)"CAINE bx"
        call flush(6)
      end if

      write(6,*)"CAINE c"
      call flush(6)

*     **** cell information ****
      version = control_version()
      do i=1,3
         ngrid(i) = control_ngrid(i)
         do j=1,3
           unita(i,j) = control_unita(i,j)
         end do
      end do

      write(6,*)"CAINE d"
      call flush(6)
*     **** brillioun zone information ****
      nbrillioun = brillioun_nbrillioun()
      value = MA_push_get(mt_dbl,3*nbrillioun,'kv',kv(2),kv(1))
      do nb=1,nbrillioun
         dbl_mb(kv(1)+3*(nb-1))   = brillioun_k(1,nb)
         dbl_mb(kv(1)+3*(nb-1)+1) = brillioun_k(2,nb)
         dbl_mb(kv(1)+3*(nb-1)+2) = brillioun_k(3,nb)
      end do

      write(6,*)"CAINE e"
      call flush(6)
      call get_psp_type(psp_filename,psp_type)

      write(6,*)"CAINE f"
      call flush(6)
      if (psp_type.eq.1) then
         value = value.and.
     >     hghpp_band(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            nbrillioun,dbl_mb(kv(1)))
      end if
      if (psp_type.eq.7) then
         value = value.and.
     >     rkbpp_band(.false.,version,psp_filename,formatted_filename,
     >           ngrid,unita,lmax,dbl_mb(kv(1)))
         call set_two_component_pseudopotential() 
         write(6,*)"RKBPP GENERATED"
         flush(6)
      end if 
      if ((psp_type.ne.7).and.(psp_type.ne.1)) then
         value = value.and.
     >     kbpp_band(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,
     >            nbrillioun,dbl_mb(kv(1)))
         write(6,*)"KBPP GENERATED"
         flush(6)
      end if

      value = value.and.MA_pop_stack(kv(2))
*     *** exit if error found ***
      if(.not. value) call errquit('psp_formatter_auto: error',0,
     &       MA_ERR)

      return
      end


