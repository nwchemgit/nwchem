      subroutine psp_formatter_input(rtdb)
c
c $Id$
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 10)
      character*26 dirs(num_dirs)
      character*255 test

      character*50 rtdb_name
      logical      reset_kbpp_number
      integer      kbpp_number,ia,kbpp_len
      character*50 kbpp_name
      character*50 cell_name,psp_filename,formatted_filename
      character*50 formatted_stress_filename
      integer locp,lmax
      character l_read
      real*8    rlocal

      data dirs / 'cell_name:',
     >            'psp_filename:',
     >            'formatted_filename:',
     >            'locp:',
     >            'lmax:',
     >            'rlocal:',
     >            'r_local:',
     >            'reset',
     >            'formatted_stress_filename:',
     >            'end'/
c
c     
c
c
c
      reset_kbpp_number = .false.
      cell_name          = 'cell_default'
      psp_filename       = ' '
      formatted_filename = ' '
      formatted_stress_filename = ' '
      locp = -1
      lmax = -1
      rlocal = 1.0d0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'psp_formatter_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'psp_formatter_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'psp_formatter_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,  200, 300, 400, 500, 600, 600, 700, 800,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'psp_formatter_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c psp_filename
c
 200  if (.not. inp_a(psp_filename))
     >     call errquit(
     >          'psp_formatter_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c  formatted_filename
c 
 300  if (.not. inp_a(formatted_filename))
     >     call errquit(
     >          'psp_formatter_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10

c
c  locp
c 
 400  if (.not. inp_a(l_read))
     >     call errquit(
     >          'psp_formatter_input: failed to read locp', 0,
     &       INPUT_ERR)
      if ((l_read.eq.'s').or.(l_read.eq.'S')) locp = 0
      if ((l_read.eq.'0'))                    locp = 0
      if ((l_read.eq.'p').or.(l_read.eq.'P')) locp = 1
      if ((l_read.eq.'1'))                    locp = 1
      if ((l_read.eq.'d').or.(l_read.eq.'D')) locp = 2
      if ((l_read.eq.'2'))                    locp = 2
      if ((l_read.eq.'f').or.(l_read.eq.'F')) locp = 3
      if ((l_read.eq.'3'))                    locp = 3
      goto 10

c
c  lmax
c 
 500  if (.not. inp_a(l_read))
     >     call errquit(
     >          'psp_formatter_input: failed to read locp', 0,
     &       INPUT_ERR)
      if ((l_read.eq.'s').or.(l_read.eq.'S')) lmax = 0
      if ((l_read.eq.'0'))                    lmax = 0
      if ((l_read.eq.'p').or.(l_read.eq.'P')) lmax = 1
      if ((l_read.eq.'1'))                    lmax = 1
      if ((l_read.eq.'d').or.(l_read.eq.'D')) lmax = 2
      if ((l_read.eq.'2'))                    lmax = 2
      if ((l_read.eq.'f').or.(l_read.eq.'F')) lmax = 3
      if ((l_read.eq.'3'))                    lmax = 3
      goto 10

c
c  rlocal
c 
 600  if (.not. inp_f(rlocal))
     >     call errquit(
     >          'psp_formatter_input: failed to read rlocal', 0,
     &       INPUT_ERR)
      goto 10

c
c  reset (kbpp_number)
c 
 700  reset_kbpp_number = .true.
      goto 10

c
c  formatted_stress_filename
c 
 800  if (.not. inp_a(formatted_stress_filename))
     >     call errquit(
     >          'psp_formatter_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


*     ***** add simulation_cell to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('psp_formatter_input: nead a cell_name', 0,
     &       INPUT_ERR)
      if (psp_filename.eq.' ') 
     >  call errquit('psp_formatter_input: nead a psp_filename', 0,
     &       INPUT_ERR)
      if (formatted_filename.eq.' ') 
     >  call errquit(
     >       'psp_formatter_input: nead a formatted_filename', 0,
     &       INPUT_ERR)

      if (formatted_stress_filename.eq.' ') then
        ind = index(formatted_filename,' ') - 1
        formatted_stress_filename = formatted_filename(1:ind)//'2'
      end if

*     **** read kbpp_number ***
      if (.not. rtdb_get(rtdb,'kbpp_number',mt_int,1,kbpp_number)) 
     >    kbpp_number = 0
      kbpp_number = kbpp_number + 1
      if (reset_kbpp_number) kbpp_number = 1

*     **** define kbpp_name  - not very elegent and could break if ****
*     ****                     kbpp_number becomes very large      **** 
      ia        = ICHAR('a')
      kbpp_name = 'kbpp'//CHAR(kbpp_number-1+ia)
      kbpp_len  = index(kbpp_name,' ') - 1

*     **** output new kbpp_number ****
      if (.not. rtdb_put(rtdb,'kbpp_number',
     >                     mt_int,1,kbpp_number))
     >     call errquit(
     >     'psp_formatter_input: rtdb_put failed', 0, RTDB_ERR)


      rtdb_name = kbpp_name(1:kbpp_len)//':cell_name'
      ind = index(cell_name,' ') - 1
      if (.not. rtdb_cput(rtdb,rtdb_name,1,cell_name(1:ind)))
     >     call errquit(
     >     'psp_formatter_input: rtdb_cput failed', 0, RTDB_ERR)

      rtdb_name = kbpp_name(1:kbpp_len)//':psp_filename'
      ind = index(psp_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,rtdb_name,
     >                    1,psp_filename(1:ind)))
     >     call errquit(
     >     'psp_formatter_input: rtdb_cput failed', 0, RTDB_ERR)

      rtdb_name = kbpp_name(1:kbpp_len)//':formatted_filename'
      ind = index(formatted_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,rtdb_name,
     >                     1,formatted_filename(1:ind)))
     >     call errquit(
     >     'psp_formatter_input: rtdb_cput failed', 0, RTDB_ERR)
      
      rtdb_name = kbpp_name(1:kbpp_len)//':formatted_stress_filename'
      ind = index(formatted_stress_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,rtdb_name,
     >                     1,formatted_stress_filename(1:ind)))
     >     call errquit(
     >     'psp_formatter_input: rtdb_cput failed', 0, RTDB_ERR)

      rtdb_name = kbpp_name(1:kbpp_len)//':locp'
      if (.not. rtdb_put(rtdb,rtdb_name,
     >                     mt_int,1,locp))
     >     call errquit(
     >     'psp_formatter_input: rtdb_put failed', 0, RTDB_ERR)

      rtdb_name = kbpp_name(1:kbpp_len)//':lmax'
      if (.not. rtdb_put(rtdb,rtdb_name,
     >                     mt_int,1,lmax))
     >     call errquit(
     >     'psp_formatter_input: rtdb_put failed', 0, RTDB_ERR)

      rtdb_name = kbpp_name(1:kbpp_len)//':rlocal'
      if (.not. rtdb_put(rtdb,rtdb_name,
     >                     mt_dbl,1,rlocal))
     >     call errquit(
     >     'psp_formatter_input: rtdb_put failed', 0, RTDB_ERR)


      return
      end
