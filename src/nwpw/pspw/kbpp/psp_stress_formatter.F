c
c $Id$
c

      logical function psp_stress_formatter(rtdb)
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
      character*50  cell_name,psp_filename,formatted_stress_filename
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
      logical  kbppv3_stress
      external kbppv3_stress

      taskid = ga_nodeid()
      value  = .true.
*     **** read kbpp_number ***
      if(.not. btdb_get(rtdb,'kbpp_number',mt_int,1,kbpp_number))
     >    call errquit(
     >       'psp_stress_formatter: kbpp_number not read', 0, RTDB_ERR)


*     **** loop over kbpp_number ****
      do kn = 1,kbpp_number
*       **** define kbpp_name  - not very elegent and could break if ****
*       ****                     kbpp_number becomes very large      **** 
        ia        = ICHAR('a')
        kbpp_name = 'kbpp'//CHAR(kn-1+ia)
        kbpp_len  = index(kbpp_name,' ') - 1
        if (taskid.eq.MASTER) then
         write(*,*) 
         write(*,*) 'Formatting rtdb psp_stress inputblock: ',kbpp_name
         write(*,*) '   ',kn,' of ',kbpp_number
         write(*,*) 
        end if


     

      rtdb_name = kbpp_name(1:kbpp_len)//':locp'
      value = btdb_get(rtdb,rtdb_name,
     >                  mt_int,1,locp)

      if(.not. value) call errquit('psp_stress_formatter: error',1,
     &       RTDB_ERR)
      
      rtdb_name = kbpp_name(1:kbpp_len)//':lmax'
      value = value.and.btdb_get(rtdb,rtdb_name,
     >                  mt_int,1,lmax)

      if(.not. value) call errquit('psp_stress_formatter: error',2,
     &       RTDB_ERR)
      
      rtdb_name = kbpp_name(1:kbpp_len)//':rlocal'
      value = value.and.btdb_get(rtdb,rtdb_name,
     >                  mt_dbl,1,rlocal)

      if(.not. value) call errquit('psp_stress_formatter: error',3,
     &       RTDB_ERR)
      rtdb_name = kbpp_name(1:kbpp_len)//':cell_name'
      value = value.and.btdb_cget(rtdb,rtdb_name,
     >                  1,cell_name)

      if(.not. value) call errquit('psp_stress_formatter: error',4,
     &       RTDB_ERR)
      rtdb_name = kbpp_name(1:kbpp_len)//':psp_filename'
      value = value.and.btdb_cget(rtdb,rtdb_name,
     >                  1,psp_filename)

      if(.not. value) call errquit('psp_stress_formatter: error',5,
     &       RTDB_ERR)
      rtdb_name = kbpp_name(1:kbpp_len)//':formatted_stress_filename'
      value = value.and.btdb_cget(rtdb,rtdb_name,
     >                  1,formatted_stress_filename)

      if(.not. value) call errquit('psp_stress_formatter: error',6,
     &       RTDB_ERR)

      l = index(cell_name,' ') - 1
      rtdb_unita   = cell_name(1:l)//':unita'
      rtdb_ngrid   = cell_name(1:l)//':ngrid'
      rtdb_boundry = cell_name(1:l)//':boundry'

      value = value.and.btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)
      value = value.and.btdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)
      value = value.and.btdb_cget(rtdb,rtdb_boundry,1,boundry)
      if(.not. value) call errquit('psp_stress_formatter: error',7,
     &       RTDB_ERR)


      l = index(boundry,' ') - 1

*     ***** run Periodic code *****
      if (boundry(1:l).eq.'periodic') then
         value = value.and.
     >     kbppv3_stress(.true.,3,psp_filename,
     >            formatted_stress_filename,
     >            ngrid,unita,locp,lmax,rlocal)

*     ***** run APeriodic code *****
      else 
        write(*,*) "Only periodic boundary conditions appropriate: ",
     >             boundry(1:l)
        value = .false.
      end if

*     *** exit if error found ***
      if(.not. value) call errquit('psp_stress_formatter: error',0,
     &       RTDB_ERR)

      end do
*     **** end kbpp_number loop ***

      psp_stress_formatter = value
      return
      end



*     ***************************************
*     *					    *
*     *		psp_stress_formatter_auto   *
*     *					    *
*     ***************************************
      subroutine psp_stress_formatter_auto(atom)
      implicit none
      character*4  atom

#include "errquit.fh"

*     ***** local variables *****
      logical       value
      integer       i,j,l,psp_type
      integer	    locp,lmax
      character*50  psp_filename,formatted_filename
      character*5   element
      real*8        rlocal

*     ***** simulation cell  variables *****
      double precision unita(3,3)
      integer          ngrid(3),version

*     ***** external functions *****
      integer      control_version
      integer      control_ngrid
      real*8       control_unita
      logical      kbppv3d_stress,nwpw_filefind,kbppv3e_stress
      external     control_version
      external     control_ngrid
      external     control_unita
      external     kbppv3d_stress,nwpw_filefind,kbppv3e_stress

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
      formatted_filename = element(1:l)//'.vpp2'

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



      if (version .eq.3) then
         call get_psp_type(psp_filename,psp_type)
         if (psp_type.eq.2) then
            value = value.and.
     >              kbppv3e_stress(.false.,
     >              version,psp_filename,formatted_filename,ngrid,unita,
     >              locp,lmax,rlocal)
         else
            value = value.and.
     >              kbppv3d_stress(.false.,
     >              version,psp_filename,formatted_filename,ngrid,unita,
     >              locp,lmax,rlocal)
         end if
      else
        write(*,*) 
     >  "Only periodic boundary conditions appropriate for stress."
        value = .false.
      end if

*     *** exit if error found ***
      if(.not. value) call errquit('psp_stress_formatter_auto: error',0,
     &       RTDB_ERR)

      return
      end


