c
c $Id$
c


*     *****************************************************
*     *                                                   *
*     *                   band_dplot_input                *
*     *                                                   *
*     *****************************************************

      subroutine band_dplot_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 10)
      character*22 dirs(num_dirs)
      data dirs / 'vectors',
     >            'density',
     >            'orbital',
     >            'position_tolerance',
     >            'elf',
     >            '2d_grid',
     >            '3d_grid',
     >            'origin',
     >            'limitxyz',
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


      integer num_fnames          ! No. of orbital directives
      parameter (num_fnames = 3)
      character*22 fnames(num_fnames)
      data fnames / 'density',
     >              'real',
     >              'complex' /

      logical value
      integer ind               ! Index of matched directive
      integer       number,count,ia,nb,onumber
      integer       name1_len,name2_len,name3_len,name4_len
      character*50  name1,name2,name3,name4
      character*50  filename
      character*50  wavefunction_filename
      character*255 test
      real*8        position_tolerance
      real*8        o(3),x(3),y(3),z(3)
      real*8        sizex(2),sizey(2),sizez(2),scal
      integer       nx,ny,nz
c
c
*     *** initializations ****
      position_tolerance=0.001d0
      call util_file_prefix('movecs',wavefunction_filename)
      ia    = ICHAR('a')
      count = 0
      value = rtdb_put(rtdb,'band_dplot:count',mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'band_dplot_input: rtdb_put failed for count', 0, RTDB_ERR)


 10   if (.not. inp_read()) 
     >     call errquit(
     >           'band_dplot_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'band_dplot_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'band_dplot_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,  200, 300, 400, 500, 600,700,800,900,
     >      9999) ind
      call errquit(
     >      'band_dplot_input: unimplemented directive', ind, INPUT_ERR)


c     
c  vectors
c     
 100  if (.not. inp_a(wavefunction_filename))
     >    call errquit(
     >         'band_dplot_input: failed to read vector', 0, INPUT_ERR)
      goto 10

c     
c  density
c     
*     *** read density option ****
 200  if (.not. inp_a(test))
     >   call errquit(
     >       'band_dplot_input: failed to read density keyword', 0,
     &       INPUT_ERR)

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
     >         'band_dplot_input: failed to read density filename', 0,
     &       INPUT_ERR)

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'band_dplot:filename'//CHAR(count-1+ia)
      name2 = 'band_dplot:number'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,'band_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'band_dplot_input: rtdb_put failed for density', 0,
     &       RTDB_ERR)

      goto 10


c     
c  orbital
c     
*     *** read orbital number ****
 300  if (.not. inp_a(test))
     >   call errquit(
     >       'band_dplot_input: failed to read orbital keyword', 0,
     &       INPUT_ERR)

*     *** orbital number ****
      if (.not.inp_match(num_fnames,.false.,test,fnames,onumber))
     >     onumber = 1
*     !*** onumber = 1 - density
*     !*** onumber = 2 - real
*     !*** onumber = 3 - complex

      if (.not. inp_i(number))
     >   call errquit(
     >       'band_dplot_input: failed to read orbital number', 0,
     &       INPUT_ERR)

*     *** read brillioun zone number ***
      if (.not. inp_i(nb)) nb=1

*     *** read filename ****
      if (.not. inp_a(filename))
     >    call errquit(
     >         'band_dplot_input: failed to read orbital filename', 0,
     &       INPUT_ERR)

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'band_dplot:filename'//CHAR(count-1+ia)
      name2 = 'band_dplot:number'//CHAR(count-1+ia)
      name3 = 'band_dplot:nb'//CHAR(count-1+ia)
      name4 = 'band_dplot:onumber'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1
      name3_len = index(name3,' ') - 1
      name4_len = index(name4,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,name3(1:name3_len),
     >                           mt_int,1,nb)
      value = value.and.rtdb_put(rtdb,name4(1:name4_len),
     >                           mt_int,1,onumber)
      value = value.and.rtdb_put(rtdb,'band_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'band_dplot_input: rtdb_put failed for orbital', 0, RTDB_ERR)

      goto 10
c     
c position_tolerance
c     
*     **** read position_tolerance ****
 400  if (.not. inp_f(position_tolerance))
     >   call errquit(
     >       'band_dplot_input: failed to read position tolerance', 0,
     &       INPUT_ERR)

      goto 10

c    
c ELF
c     
 500  if (.not. inp_a(test))
     >   call errquit(
     >       'band_dplot_input: failed to read ELF keyword', 0,
     &       INPUT_ERR)

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
      name1 = 'band_dplot:filename'//CHAR(count-1+ia)
      name2 = 'band_dplot:number'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,'band_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'band_dplot_input: rtdb_put failed for ELF', 0, RTDB_ERR)


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
      

      value = value.and.rtdb_put(rtdb,'band_dplot:2d_grid:o',
     >                           mt_dbl,3,o)
      value = value.and.rtdb_put(rtdb,'band_dplot:2d_grid:x',
     >                           mt_dbl,3,x)
      value = value.and.rtdb_put(rtdb,'band_dplot:2d_grid:y',
     >                           mt_dbl,3,y)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:2d_grid:sizex',
     >                  mt_dbl,2,sizex)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:2d_grid:sizey',
     >                  mt_dbl,2,sizey)

      if (.not. value)
     >   call errquit(
     >    'band_dplot_input: 2d_grid failed to read', 0, RTDB_ERR)

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
      value =  value.and.inp_i(nx)

      value =  value.and.inp_read()
      value =  value.and.inp_f(y(1))
      value =  value.and.inp_f(y(2))
      value =  value.and.inp_f(y(3))
      value =  value.and.inp_f(sizey(1))
      value =  value.and.inp_f(sizey(2))
      value =  value.and.inp_i(ny)

      value =  value.and.inp_read()
      value =  value.and.inp_f(z(1))
      value =  value.and.inp_f(z(2))
      value =  value.and.inp_f(z(3))
      value =  value.and.inp_f(sizez(1))
      value =  value.and.inp_f(sizez(2))
      value =  value.and.inp_i(nz)

      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:o',
     >                           mt_dbl,3,o)
      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:x',
     >                           mt_dbl,3,x)
      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:y',
     >                           mt_dbl,3,y)
      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:z',
     >                           mt_dbl,3,z)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:sizex',
     >                  mt_dbl,2,sizex)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:sizey',
     >                  mt_dbl,2,sizey)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:sizez',
     >                  mt_dbl,2,sizez)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:nx',
     >                  mt_int,1,nx)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:ny',
     >                  mt_int,1,ny)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:nz',
     >                  mt_int,1,nz)

      if (.not. value)
     >   call errquit(
     >    'band_dplot_input: 3d_grid failed to read', 0, RTDB_ERR)

      goto 10


c
c translate_origin
c
 800  value =            inp_f(o(1))
      value =  value.and.inp_f(o(2))
      value =  value.and.inp_f(o(3))

      value = value.and.rtdb_put(rtdb,'band_dplot:origin',
     >                           mt_dbl,3,o)

      if (.not. value)
     >   call errquit(
     >    'band_dplot_input: translate_origin failed to read',0,0)

      goto 10


c
c limitxyz
c
*     **** read o,x,sizex,y,sizey,y,sizez ****
 900  call get_scalefrominput(scal)

      value =  value.and.inp_read()
      value =  value.and.inp_f(sizex(1))
      value =  value.and.inp_f(sizex(2))
      value =  value.and.inp_i(nx)
      sizex(1) = scal*sizex(1)
      sizex(2) = scal*sizex(2)

      value =  value.and.inp_read()
      value =  value.and.inp_f(sizey(1))
      value =  value.and.inp_f(sizey(2))
      value =  value.and.inp_i(ny)
      sizey(1) = scal*sizey(1)
      sizey(2) = scal*sizey(2)

      value =  value.and.inp_read()
      value =  value.and.inp_f(sizez(1))
      value =  value.and.inp_f(sizez(2))
      value =  value.and.inp_i(nz)
      sizez(1) = scal*sizez(1)
      sizez(2) = scal*sizez(2)

      !*** set origin and axes ***
      o(1) = 0.0d0
      o(2) = 0.0d0
      o(3) = 0.0d0

      x(1) = 1.0d0
      x(2) = 0.0d0
      x(3) = 0.0d0

      y(1) = 0.0d0
      y(2) = 1.0d0
      y(3) = 0.0d0

      z(1) = 0.0d0
      z(2) = 0.0d0
      z(3) = 1.0d0

      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:o',
     >                           mt_dbl,3,o)
      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:x',
     >                           mt_dbl,3,x)
      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:y',
     >                           mt_dbl,3,y)
      value = value.and.rtdb_put(rtdb,'band_dplot:3d_grid:z',
     >                           mt_dbl,3,z)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:sizex',
     >                  mt_dbl,2,sizex)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:sizey',
     >                  mt_dbl,2,sizey)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:sizez',
     >                  mt_dbl,2,sizez)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:nx',
     >                  mt_int,1,nx)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:ny',
     >                  mt_int,1,ny)
      value = value.and.
     >        rtdb_put(rtdb,'band_dplot:3d_grid:nz',
     >                  mt_int,1,nz)

      if (.not. value)
     >   call errquit(
     >  'band_dplot_input: failed to read limitxyz',0,RTDB_ERR)

      goto 10






 9999 continue

      ind = index(wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'band_dplot:wavefunction_filename',
     >                    1,wavefunction_filename(1:ind))
      value = value.and.
     >         rtdb_put(rtdb,'band_dplot:position_tolerance',
     >                  mt_dbl,1,position_tolerance)
      if (.not.value)
     >     call errquit(
     >     'band_dplot_input: rtdb_put failed for vector', 0, RTDB_ERR)

      return
      end
