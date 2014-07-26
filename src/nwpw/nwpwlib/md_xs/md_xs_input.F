c $Id: md_xs_input.F 25011 2013-12-19 17:45:25Z d3y133 $

*     *********************************************
*     *                                           *
*     *              md_xs_input                  *
*     *                                           *
*     *********************************************
      subroutine md_xs_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"



      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 5)
      character*18 dirs(num_dirs)
      data dirs / 'spectroscopy',
     >            'absorber',
     >            'edge',
     >            'center',
     >            'end'/

c
      logical       value
      character*80  rtdb_name
      character*255 test,test2
      integer ind               ! Index of matched directive



 10   if (.not. inp_read()) 
     >  call errquit(
     >  'md_xs_input: inp_read failed',0,INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'md_xs_input: failed to read keyword',0,INPUT_ERR)

      if (.not. inp_match(num_dirs,.false.,test,dirs,ind))
     >  call errquit(
     >  'md_xs_input: unknown_directive',ind,INPUT_ERR)

      goto ( 100, 200,300, 400,
     >      9999) ind
      call errquit(
     >     'md_xs_input: unimplemented directive',ind,INPUT_ERR)

c
c  spectroscopy
c
 100  if (inp_a(test2)) then
         rtdb_name = 'ms_xs:spectroscopy'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing spectroscopy',100,RTDB_ERR)
      end if

      goto 10

c
c  absorbers
c
 200  if (inp_a(test2)) then
         rtdb_name = 'md_xs:absorber'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',200,RTDB_ERR)
      end if

      goto 10

c
c  edge
c
 300  if (inp_a(test2)) then
         rtdb_name = 'md_xs:edge'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',300, RTDB_ERR)
      end if

      goto 10

c
c  center
c
 400  if (inp_a(test2)) then
         rtdb_name = 'md_xs:center'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',400, RTDB_ERR)
      end if

      goto 10


 9999 continue
      
      return
      end
