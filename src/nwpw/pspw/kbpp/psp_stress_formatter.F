      logical function psp_stress_formatter(rtdb)
c
c $Id: psp_stress_formatter.F,v 1.1 2001-08-30 00:38:59 edo Exp $
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
c
      integer rtdb
      logical      value

      integer taskid,MASTER
      parameter (MASTER=0)

*     ***** rtdb variables *****
      real*8        rlocal
      integer       l
      integer	    locp,lmax
      character*25  cell_name,psp_filename,formatted_stress_filename
      character*30  rtdb_unita
      character*30  rtdb_ngrid
      character*30  rtdb_boundry
      integer       ia,kn,kbpp_len,kbpp_number
      character*50  rtdb_name,kbpp_name

*     ***** simulation cell  variables *****
      character*25     boundry
      double precision unita(3,3)
      integer          ngrid(3)

*     ***** external functions *****
      logical  kbppv3_stress
      external kbppv3_stress

      taskid = ga_nodeid()
      value  = .true.
*     **** read kbpp_number ***
      if(.not. rtdb_get(rtdb,'kbpp_number',mt_int,1,kbpp_number))
     >    call errquit(
     >       'psp_stress_formatter: kbpp_number not read', 0)


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
      value = rtdb_get(rtdb,rtdb_name,
     >                  mt_int,1,locp)

      if(.not. value) call errquit('psp_stress_formatter: error',1)
      
      rtdb_name = kbpp_name(1:kbpp_len)//':lmax'
      value = value.and.rtdb_get(rtdb,rtdb_name,
     >                  mt_int,1,lmax)

      if(.not. value) call errquit('psp_stress_formatter: error',2)
      
      rtdb_name = kbpp_name(1:kbpp_len)//':rlocal'
      value = value.and.rtdb_get(rtdb,rtdb_name,
     >                  mt_dbl,1,rlocal)

      if(.not. value) call errquit('psp_stress_formatter: error',3)
      rtdb_name = kbpp_name(1:kbpp_len)//':cell_name'
      value = value.and.rtdb_cget(rtdb,rtdb_name,
     >                  1,cell_name)

      if(.not. value) call errquit('psp_stress_formatter: error',4)
      rtdb_name = kbpp_name(1:kbpp_len)//':psp_filename'
      value = value.and.rtdb_cget(rtdb,rtdb_name,
     >                  1,psp_filename)

      if(.not. value) call errquit('psp_stress_formatter: error',5)
      rtdb_name = kbpp_name(1:kbpp_len)//':formatted_stress_filename'
      value = value.and.rtdb_cget(rtdb,rtdb_name,
     >                  1,formatted_stress_filename)

      if(.not. value) call errquit('psp_stress_formatter: error',6)

      l = index(cell_name,' ') - 1
      rtdb_unita   = cell_name(1:l)//':unita'
      rtdb_ngrid   = cell_name(1:l)//':ngrid'
      rtdb_boundry = cell_name(1:l)//':boundry'

      value = value.and.rtdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)
      value = value.and.rtdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)
      value = value.and.rtdb_cget(rtdb,rtdb_boundry,1,boundry)
      if(.not. value) call errquit('psp_stress_formatter: error',7)


      l = index(boundry,' ') - 1

*     ***** run Periodic code *****
      if (boundry(1:l).eq.'periodic') then
         value = value.and.
     >     kbppv3_stress(3,psp_filename,formatted_stress_filename,
     >            ngrid,unita,locp,lmax,rlocal)

*     ***** run APeriodic code *****
      else 
        write(*,*) "Only periodic boundary conditions appropriate: ",
     >             boundry(1:l)
        value = .false.
      end if

*     *** exit if error found ***
      if(.not. value) call errquit('psp_stress_formatter: error',0)

      end do
*     **** end kbpp_number loop ***

      psp_stress_formatter = value
      return
      end
