c
c $Id: pspw_dplot_input.F,v 1.6 2002-12-10 01:37:23 bylaska Exp $
c
      subroutine pspw_dplot_input(rtdb)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 8)
      character*22 dirs(num_dirs)
      data dirs / 'vectors',
     >            'density',
     >            'orbital',
     >            'position_tolerance',
     >            'elf',
     >            '2d_grid',
     >            '3d_grid',
     >            'end'/

      integer num_dnames          ! No. of density directives
      parameter (num_dnames = 8)
      character*22 dnames(num_dnames)
      data dnames / 'total',
     >              'diff',
     >              'alpha',
     >              'beta',
     >              'laplacian',
     >              'potential',
     >              'up',
     >              'down' /
      integer num_enames          ! No. of ELF directives
      parameter (num_enames = 5)
      character*22 enames(num_enames)
      data enames / 'restricted',
     >              'up',
     >              'down',
     >              'alpha',
     >              'beta'   /

      logical value
      integer ind               ! Index of matched directive
      integer       number,count,ia
      integer       name1_len,name2_len
      character*30  name1,name2
      character*30  filename
      character*30  wavefunction_filename
      character*255 test
      real*8        position_tolerance
      real*8        o(3),x(3),y(3),z(3)
      real*8        sizex(2),sizey(2),sizez(2)
c
c
*     *** initializations ****
      position_tolerance=0.001d0
      call util_file_prefix('movecs',wavefunction_filename)
      ia    = ICHAR('a')
      count = 0
      value = rtdb_put(rtdb,'pspw_dplot:count',mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'pspw_dplot_input: rtdb_put failed for count', 0)


 10   if (.not. inp_read()) 
     >     call errquit(
     >           'pspw_dplot_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >           'pspw_dplot_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'pspw_dplot_input: unknown directive', 0)


      goto ( 100,  200, 300, 400, 500, 600,700,
     >      9999) ind
      call errquit(
     >      'pspw_dplot_input: unimplemented directive', ind)


c     
c  vectors
c     
 100  if (.not. inp_a(wavefunction_filename))
     >    call errquit(
     >         'pspw_dplot_input: failed to read vector', 0)
      goto 10

c     
c  density
c     
*     *** read density option ****
 200  if (.not. inp_a(test))
     >   call errquit(
     >       'pspw_dplot_input: failed to read density keyword', 0)

*     *** density number ****
      if (.not.inp_match(num_dnames,.false.,test,dnames,number))
     >     number = 1
      if (number .eq. 7) number = 3 
      if (number .eq. 8) number = 4 
*     !*** number = 1 - total
*     !*** number = 2 - difference
*     !*** number = 3 - alpha
*     !*** number = 4 - beta
*     !*** number = 5 - laplacian
*     !*** number = 6 - potential
*     !*** number = 7 - ELF

*     **** make density numbers negative ****
      number = -number

*     *** read filename ****
      if (.not. inp_a(filename))
     >    call errquit(
     >         'pspw_dplot_input: failed to read density filename', 0)

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'pspw_dplot:filename'//CHAR(count-1+ia)
      name2 = 'pspw_dplot:number'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'pspw_dplot_input: rtdb_put failed for density', 0)

      goto 10


c     
c  orbital
c     
*     *** read orbital number ****
 300  if (.not. inp_i(number))
     >   call errquit(
     >       'pspw_dplot_input: failed to read orbital number', 0)

*     *** read filename ****
      if (.not. inp_a(filename))
     >    call errquit(
     >         'pspw_dplot_input: failed to read orbital filename', 0)

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'pspw_dplot:filename'//CHAR(count-1+ia)
      name2 = 'pspw_dplot:number'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'pspw_dplot_input: rtdb_put failed for orbital', 0)

      goto 10
c     
c position_tolerance
c     
*     **** read position_tolerance ****
 400  if (.not. inp_f(position_tolerance))
     >   call errquit(
     >       'pspw_dplot_input: failed to read position tolerance', 0)

      goto 10

c    
c ELF
c     
 500  if (.not. inp_a(test))
     >   call errquit(
     >       'pspw_dplot_input: failed to read ELF keyword', 0)

*     *** ELF number ****
      if (.not.inp_match(num_enames,.false.,test,enames,number))
     >     number = -7
      if (number .eq. 1) number = -7 
      if (number .eq. 2) number = -7
      if (number .eq. 3) number = -8
      if (number .eq. 4) number = -7
      if (number .eq. 5) number = -8
*     !*** number = -7 - ELF up
*     !*** number = -8 - ELF down

*      *** read ELF filename ****  
      if (.not. inp_a(filename))
     >   filename = test

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'pspw_dplot:filename'//CHAR(count-1+ia)
      name2 = 'pspw_dplot:number'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'pspw_dplot_input: rtdb_put failed for ELF', 0)


      goto 10

c     
c 2d_grid
c     
*     **** read o,x,sizex,y,sizey ****
 600  value =            inp_read()
      value =  value.and.inp_f(o(1))
      value =  value.and.inp_f(o(2))
      value =  value.and.inp_f(o(3))

      value =  value.and.inp_read()
      value =  value.and.inp_f(x(1))
      value =  value.and.inp_f(x(2))
      value =  value.and.inp_f(x(3))
      value =  value.and.inp_f(sizex(1))
      value =  value.and.inp_f(sizex(2))

      value =  value.and.inp_read()
      value =  value.and.inp_f(y(1))
      value =  value.and.inp_f(y(2))
      value =  value.and.inp_f(y(3))
      value =  value.and.inp_f(sizey(1))
      value =  value.and.inp_f(sizey(2))
      

      value = value.and.rtdb_put(rtdb,'pspw_dplot:2d_grid:o',
     >                           mt_dbl,3,o)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:2d_grid:x',
     >                           mt_dbl,3,x)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:2d_grid:y',
     >                           mt_dbl,3,y)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:2d_grid:sizex',
     >                  mt_dbl,2,sizex)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:2d_grid:sizey',
     >                  mt_dbl,2,sizey)

      if (.not. value)
     >   call errquit(
     >    'pspw_dplot_input: 2d_grid failed to read', 0)

      goto 10

c     
c 3d_grid
c     
*     **** read o,x,sizex,y,sizey,y,sizez ****
 700  value =            inp_read()
      value =  value.and.inp_f(o(1))
      value =  value.and.inp_f(o(2))
      value =  value.and.inp_f(o(3))

      value =  value.and.inp_read()
      value =  value.and.inp_f(x(1))
      value =  value.and.inp_f(x(2))
      value =  value.and.inp_f(x(3))
      value =  value.and.inp_f(sizex(1))
      value =  value.and.inp_f(sizex(2))

      value =  value.and.inp_read()
      value =  value.and.inp_f(y(1))
      value =  value.and.inp_f(y(2))
      value =  value.and.inp_f(y(3))
      value =  value.and.inp_f(sizey(1))
      value =  value.and.inp_f(sizey(2))

      value =  value.and.inp_read()
      value =  value.and.inp_f(z(1))
      value =  value.and.inp_f(z(2))
      value =  value.and.inp_f(z(3))
      value =  value.and.inp_f(sizez(1))
      value =  value.and.inp_f(sizez(2))

      value = value.and.rtdb_put(rtdb,'pspw_dplot:3d_grid:o',
     >                           mt_dbl,3,o)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:3d_grid:x',
     >                           mt_dbl,3,x)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:3d_grid:y',
     >                           mt_dbl,3,y)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:3d_grid:z',
     >                           mt_dbl,3,z)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:sizex',
     >                  mt_dbl,2,sizex)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:sizey',
     >                  mt_dbl,2,sizey)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:sizez',
     >                  mt_dbl,2,sizez)

      if (.not. value)
     >   call errquit(
     >    'pspw_dplot_input: 3d_grid failed to read', 0)

      goto 10

 9999 continue

      ind = index(wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'pspw_dplot:wavefunction_filename',
     >                    1,wavefunction_filename(1:ind))
      value = value.and.
     >         rtdb_put(rtdb,'pspw_dplot:position_tolerance',
     >                  mt_dbl,1,position_tolerance)
      if (.not.value)
     >     call errquit(
     >     'pspw_dplot_input: rtdb_put failed for vector', 0)

      return
      end
