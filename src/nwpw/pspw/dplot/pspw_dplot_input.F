c
c $Id$
c


*     *****************************************************
*     *                                                   *
*     *                   pspw_dplot_input                *
*     *                                                   *
*     *****************************************************

      subroutine pspw_dplot_input(rtdb)
      implicit none
      integer rtdb
c
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
c
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 16)
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
     >            'ncell',
     >            'atom_truncate',
     >            '1d_grid',
     >            'orbital2',
     >            'density_matrix',
     >            'structure_factor',
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
      integer       number,number2,count,ia
      integer       name1_len,name2_len,name3_len,name4_len
      character*50  name1,name2,name3,name4
      character*50  filename
      character*50  wavefunction_filename
      character*255 test
      real*8        position_tolerance
      real*8        o(3),x(3),y(3),z(3)
      real*8        sizex(2),sizey(2),sizez(2),scal
      integer       j,jstart,jlast,jstride
      integer       nx,ny,nz,ncell(3),idx(2)
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
     >     'pspw_dplot_input: rtdb_put failed for count', 0, RTDB_ERR)


 10   if (.not. inp_read()) 
     >     call errquit(
     >           'pspw_dplot_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'pspw_dplot_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'pspw_dplot_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,  200, 300, 400, 500, 600,700,800,900,1000,1100,1200,
     >       350, 360, 370,
     >      9999) ind
      call errquit(
     >      'pspw_dplot_input: unimplemented directive', ind, INPUT_ERR)


c     
c  vectors
c     
 100  if (.not. inp_a(wavefunction_filename))
     >    call errquit(
     >         'pspw_dplot_input: failed to read vector', 0, INPUT_ERR)
      goto 10

c     
c  density
c     
*     *** read density option ****
 200  if (.not. inp_a(test))
     >   call errquit(
     >       'pspw_dplot_input: failed to read density keyword', 0,
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
     >         'pspw_dplot_input: failed to read density filename', 0,
     &       INPUT_ERR)

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
     >     'pspw_dplot_input: rtdb_put failed for density', 0,
     &       RTDB_ERR)

      goto 10


c     
c  orbital
c     
*     *** read orbital number ****
 300  if (.not. inp_i(number))
     >   call errquit(
     >       'pspw_dplot_input: failed to read orbital number', 0,
     &       INPUT_ERR)

*     *** read filename ****
      if (.not. inp_a(filename))
     >    call errquit(
     >         'pspw_dplot_input: failed to read orbital filename', 0,
     &       INPUT_ERR)

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
     >     'pspw_dplot_input: rtdb_put failed for orbital', 0, RTDB_ERR)

      goto 10


c     
c  orbital2
c     
*     *** read orbital number ****
 350  if (.not. inp_i(number))
     >   call errquit(
     >       'pspw_dplot_input: failed to read orbital number1',0,
     &       INPUT_ERR)
      if (.not. inp_i(number2))
     >   call errquit(
     >       'pspw_dplot_input: failed to read orbital number2',0,
     &       INPUT_ERR)

*     *** read filename ****
      if (.not. inp_a(filename))
     >    call errquit(
     >         'pspw_dplot_input: failed to read orbital2 filename',0,
     &       INPUT_ERR)

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'pspw_dplot:filename'//CHAR(count-1+ia)
      name2 = 'pspw_dplot:number1'//CHAR(count-1+ia)
      name3 = 'pspw_dplot:number2'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1
      name3_len = index(name3,' ') - 1

      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,name3(1:name3_len),
     >                           mt_int,1,number2)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >     call errquit(
     >     'pspw_dplot_input:rtdb_put failed for orbital2',0,RTDB_ERR)

      goto 10

c     
c  density_matrix ms x,y,z
c     
*     *** read ms,x,y,z ****
 360  if (.not. inp_i(number2))
     >   call errquit(
     >       'pspw_dplot_input: failed to read density_matrix ms',
     >       0,INPUT_ERR)
      if (.not. inp_f(x(1)))
     >  call errquit('pspw_dplot_input:failed to read density_matrix x',
     >               0,INPUT_ERR)
      if (.not. inp_f(x(2)))
     >  call errquit('pspw_dplot_input:failed to read density_matrix y',
     >               0,INPUT_ERR)
      if (.not. inp_f(x(3)))
     >  call errquit('pspw_dplot_input:failed to read density_matrix z',
     >                0,INPUT_ERR)

*     *** read filename ****
      if (.not. inp_a(filename))
     >    call errquit(
     >       'pspw_dplot_input: failed to read density_matrix filename',
     >       0,INPUT_ERR)

*     **** define name  - not very elegent and could break if ****
*     ****                      count becomes very large      ****
      count = count + 1
      name1 = 'pspw_dplot:filename'//CHAR(count-1+ia)
      name2 = 'pspw_dplot:number'//CHAR(count-1+ia)
      name3 = 'pspw_dplot:dmatrix_ms'//CHAR(count-1+ia)
      name4 = 'pspw_dplot:dmatrix_xyz'//CHAR(count-1+ia)
      name1_len = index(name1,' ') - 1
      name2_len = index(name2,' ') - 1
      name3_len = index(name3,' ') - 1
      name4_len = index(name4,' ') - 1

      number = -9
      ind = index(filename,' ') - 1
      value = rtdb_cput(rtdb,name1(1:name1_len),1,filename(1:ind))
      value = value.and.rtdb_put(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
      value = value.and.rtdb_put(rtdb,name3(1:name3_len),
     >                           mt_int,1,number2)
      value = value.and.rtdb_put(rtdb,name4(1:name4_len),
     >                           mt_dbl,3,x)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:count',
     >                           mt_int,1,count)
      if (.not.value)
     >  call errquit(
     >  'pspw_dplot_input:rtdb_put failed for densty_dmatix',0,RTDB_ERR)

      goto 10



c     
c  structure_factor
c     
*     *** read structure factor option ****
 370  if (.not. inp_a(filename))
     >   call errquit(
     >   'pspw_dplot_input: failed read structure_factor filename', 0,
     &       INPUT_ERR)

      number = -10
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
     >     'pspw_dplot_input: rtdb_put failed for structure_factor', 0,
     &       RTDB_ERR)

      goto 10




c     
c position_tolerance
c     
*     **** read position_tolerance ****
 400  if (.not. inp_f(position_tolerance))
     >   call errquit(
     >       'pspw_dplot_input: failed to read position tolerance', 0,
     &       INPUT_ERR)

      goto 10

c    
c ELF
c     
 500  if (.not. inp_a(test))
     >   call errquit(
     >       'pspw_dplot_input: failed to read ELF keyword', 0,
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
     >     'pspw_dplot_input: rtdb_put failed for ELF', 0, RTDB_ERR)


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
     >    'pspw_dplot_input: 2d_grid failed to read', 0, RTDB_ERR)

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
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:nx',
     >                  mt_int,1,nx)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:ny',
     >                  mt_int,1,ny)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:nz',
     >                  mt_int,1,nz)

      if (.not. value)
     >   call errquit(
     >    'pspw_dplot_input: 3d_grid failed to read', 0, RTDB_ERR)

      goto 10


c
c translate_origin
c
 800  value =            inp_f(o(1))
      value =  value.and.inp_f(o(2))
      value =  value.and.inp_f(o(3))

      value = value.and.rtdb_put(rtdb,'pspw_dplot:origin',
     >                           mt_dbl,3,o)

      if (.not. value)
     >   call errquit(
     >    'pspw_dplot_input: translate_origin failed to read',0,0)

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
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:nx',
     >                  mt_int,1,nx)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:ny',
     >                  mt_int,1,ny)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:3d_grid:nz',
     >                  mt_int,1,nz)

      if (.not. value)
     >   call errquit(
     >  'pspw_dplot_input: failed to write limitxyz',0,RTDB_ERR)

      goto 10


c
c ncell
c
 1000 if (.not.inp_i(ncell(1))) ncell(1) = 0
      if (.not.inp_i(ncell(2))) ncell(2) = 0
      if (.not.inp_i(ncell(3))) ncell(3) = 0
      if (.not.rtdb_put(rtdb,'pspw_dplot:ncell',mt_int,3,ncell))
     >   call errquit(
     >  'pspw_dplot_input: failed to write ncell',0,RTDB_ERR)

      goto 10

c
c atom_truncate
c
 1100 if (.not.BA_push_get(mt_int,nw_max_atom,'idx',idx(2),idx(1)))
     >   call errquit(
     >   'pspw_dplot_input:failed allocating idx',0,MA_ERR)

      nx = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(idx(1)+nx) = j
          nx = nx+1
         end do
      end do
      if (.not.rtdb_put(rtdb,'pspw_dplot:atom_truncate_size',
     >                  mt_int,1,nx))
     >  call errquit('pspw_dplot_input:failed write atom_truncate_size',
     >                0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'pspw_dplot:atom_truncate',
     >                  mt_int,nx,int_mb(idx(1))))
     >   call errquit('pspw_dplot_input:failed writing atom_truncate',
     >                0,RTDB_ERR)
      if (.not.BA_pop_stack(idx(2)))
     >   call errquit('pspw_dplot_input:failed deallocating idx',
     >                 0,MA_ERR)

      goto 10

c     
c 1d_grid
c     
*     **** read o,x,sizex,y,sizey ****
 1200 value =            inp_read()
      value =  value.and.inp_i(nx)

      value =  value.and.inp_read()
      value =  value.and.inp_f(o(1))
      value =  value.and.inp_f(o(2))
      value =  value.and.inp_f(o(3))

      value =  value.and.inp_read()
      value =  value.and.inp_f(x(1))
      value =  value.and.inp_f(x(2))
      value =  value.and.inp_f(x(3))


      value = value.and.rtdb_put(rtdb,'pspw_dplot:1d_grid:o',
     >                           mt_dbl,3,o)
      value = value.and.rtdb_put(rtdb,'pspw_dplot:1d_grid:x',
     >                           mt_dbl,3,x)
      value = value.and.
     >        rtdb_put(rtdb,'pspw_dplot:1d_grid:nx',
     >                  mt_int,1,nx)


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
     >     'pspw_dplot_input: rtdb_put failed for vector', 0, RTDB_ERR)

      return
      end
