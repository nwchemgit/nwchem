c
c $Id$
c

*     ***************************************
*     *	  				    *
*     *		   psp_formatter	    *
*     *					    *
*     ***************************************

      logical function psp_formatter(rtdb)
      implicit none
#include "errquit.fh"
#include "bafdecls.fh"
#include "global.fh"
#include "btdb.fh"
c
      integer rtdb
      logical      value

      integer taskid,MASTER
      parameter (MASTER=0)

*     ***** rtdb variables *****
      real*8        rlocal
      integer       l
      integer	    locp,lmax
      character*50  cell_name,psp_filename,formatted_filename
      character*50  rtdb_unita
      character*50  rtdb_ngrid
      character*50  rtdb_boundry
      integer       ia,kn,kbpp_len,kbpp_number
      character*50  rtdb_name,kbpp_name

*     ***** simulation cell  variables *****
      character*50     boundry
      double precision unita(3,3)
      integer          ngrid(3)

*     ***** external functions *****
      logical  kbppv3,kbppv3_stress
      external kbppv3,kbppv3_stress

      taskid = ga_nodeid()
      value  = .true.
*     **** read kbpp_number ***
      if(.not. btdb_get(rtdb,'kbpp_number',mt_int,1,kbpp_number))
     >    call errquit(
     >       'psp_formatter: kbpp_number not read', 0, RTDB_ERR)


*     **** loop over kbpp_number ****
      do kn = 1,kbpp_number
*       **** define kbpp_name  - not very elegent and could break if ****
*       ****                     kbpp_number becomes very large      **** 
        ia        = ICHAR('a')
        kbpp_name = 'kbpp'//CHAR(kn-1+ia)
        kbpp_len  = index(kbpp_name,' ') - 1
        if (taskid.eq.MASTER) then
         write(*,*) 
         write(*,*) 'Formatting rtdb psp inputblock: ',kbpp_name
         write(*,*) '   ',kn,' of ',kbpp_number
         write(*,*) 
        end if



      rtdb_name = kbpp_name(1:kbpp_len)//':locp'
      value = btdb_get(rtdb,rtdb_name,
     >                  mt_int,1,locp)

      rtdb_name = kbpp_name(1:kbpp_len)//':lmax'
      value = value.and.btdb_get(rtdb,rtdb_name,
     >                  mt_int,1,lmax)

      rtdb_name = kbpp_name(1:kbpp_len)//':rlocal'
      value = value.and.btdb_get(rtdb,rtdb_name,
     >                  mt_dbl,1,rlocal)

      rtdb_name = kbpp_name(1:kbpp_len)//':cell_name'
      value = value.and.btdb_cget(rtdb,rtdb_name,
     >                  1,cell_name)

      rtdb_name = kbpp_name(1:kbpp_len)//':psp_filename'
      value = value.and.btdb_cget(rtdb,rtdb_name,
     >                  1,psp_filename)

      rtdb_name = kbpp_name(1:kbpp_len)//':formatted_filename'
      value = value.and.btdb_cget(rtdb,rtdb_name,
     >                  1,formatted_filename)

      l = index(cell_name,' ') - 1
      rtdb_unita   = cell_name(1:l)//':unita'
      rtdb_ngrid   = cell_name(1:l)//':ngrid'
      rtdb_boundry = cell_name(1:l)//':boundry'

      value = value.and.btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)
      value = value.and.btdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)
      value = value.and.btdb_cget(rtdb,rtdb_boundry,1,boundry)


      l = index(boundry,' ') - 1

*     ***** run Periodic code *****
      if (boundry(1:l).eq.'periodic') then
         value = value.and.
     >     kbppv3(.true.,3,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,rlocal)

*     ***** run APeriodic code *****
      else if (boundry(1:l).eq.'aperiodic') then
         value = value.and.
     >     kbppv3(.true.,4,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,rlocal)
      else
         write(*,*) "unknown boundry type: ",boundry(1:l)
      end if

*     *** exit if error found ***
      if(.not. value) call errquit('psp_formatter: error',0,
     &       UNKNOWN_ERR)

      end do
*     **** end kbpp_number loop ***

      psp_formatter = value
      return
      end


*     ***************************************
*     *                                     *
*     *		psp_formatter_auto          *
*     *                                     *
*     ***************************************
      subroutine psp_formatter_auto(atom)
      implicit none
      character*4  atom

#include "errquit.fh"

*     ***** local variables *****
      logical       value
      integer       i,j,l
      integer	    locp,lmax
      character*50  psp_filename,formatted_filename
      character*5   element
      real*8        rlocal

*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version,psp_type

*     ***** external functions *****
      integer      control_version
      integer      control_ngrid
      real*8       control_unita
      logical      kbppv3d,kbppv3e,hghppv1,pawppv1,nwpw_filefind
      external     control_version
      external     control_ngrid
      external     control_unita
      external     kbppv3d,kbppv3e,hghppv1,pawppv1,nwpw_filefind


*     **** get lmax,locp,and rlocal from psp library ****
      lmax  = -1
      locp  = -1
      rlocal = 1.0d0

*     **** define psp_filename and formatted_filename ****
      element = '     '
      element = atom
      l = index(element,' ') - 1
      psp_filename       = element(1:l)//'.psp'
      formatted_filename = element(1:l)//'.vpp'

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

      call get_psp_type(psp_filename,psp_type)

      if ((psp_type.eq.0).or.(psp_type.eq.9)) then
         value = kbppv3d(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,rlocal)
      else if (psp_type .eq.1) then
         value = hghppv1(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            rlocal)
      else if (psp_type.eq.2) then
         value = kbppv3e(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,rlocal)
      else if ((psp_type.eq.4).or.(psp_type.eq.6)) then
         value = pawppv1(.false.,
     >            version,psp_filename,formatted_filename,ngrid,unita,
     >            locp,lmax,rlocal)
      else 
         value = .false.
      end if


*     *** exit if error found ***
      if(.not. value) call errquit('psp_formatter_auto: error',0,
     >       UNKNOWN_ERR)

      return
      end


