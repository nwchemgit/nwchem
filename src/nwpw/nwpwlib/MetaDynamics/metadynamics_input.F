*     *********************************************
*     *                                           *
*     *           metadynamics_input              *
*     *                                           *
*     *********************************************
      subroutine metadynamics_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "errquit.fh"
c
c

      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 8)
      character*15 dirs(num_dirs)
      data dirs / 'update',
     >            'print_shift',
     >            'bond',
     >            'angle',
     >            'dihedral',
     >            'coord_number',
     >            'spline',
     >            'end'/

      integer num_dirs2          ! No. of known directives
      parameter (num_dirs2 = 12)
      character*10 dirs2(num_dirs2)
      data dirs2 / 'index1',
     >             'index2',
     >             'n',
     >             'm',
     >             'r0',
     >             'w',
     >             'sigma',
     >             'maxcoord',
     >             'range',
     >             'wab',
     >             'sigmaab',
     >             'sprik'/
c
      logical       value,spline
      character*80  rtdb_name
      character*255 test,dir_name,sdir_name
      character*4   elem1,elem2
      integer ind,ind2               ! Index of matched directive
      integer nka,nkb,nkc,nkd,nka_gs,nkb_gs,nkc_gs,nkd_gs
      integer j,jstart,jlast,jstride,n1,n2
      integer sindx(50),index_start1(2),index_start2(2)
      real*8  sdist(50)

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      nka   = 0
      nkb   = 0
      nkc   = 0
      nkd   = 0
      nka_gs  = 0
      nkb_gs  = 0
      nkc_gs  = 0
      nkd_gs  = 0
      spline = .false.

 10   if (.not. inp_read()) 
     >  call errquit(
     >  'metadyanmics_input: inp_read failed',0,INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'metadynamics_input: failed to read keyword',0,INPUT_ERR)
      if (.not. inp_match(num_dirs,.false.,test,dirs,ind))
     >  call errquit(
     >  'metadynamics_input: unknown directive',0,INPUT_ERR)

      goto ( 100, 150, 200, 300, 400, 500, 700,
     >      9999) ind
      call errquit(
     >     'metadynamics_input: unimplemented directive',ind,INPUT_ERR)

c
c  update
c
 100  if (inp_i(sindx(1))) then
         rtdb_name = 'metadynamics_update'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,sindx)
         if (.not. value) call errquit(
     >   'metadynamics_input:failed writing update',0, RTDB_ERR)
      end if

      goto 10

c
c  print_shift
c
 150  if (inp_i(sindx(1))) then
         rtdb_name = 'metadynamics_print_shift'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,sindx)
         if (.not. value) call errquit(
     >   'metadynamics_input:failed writing update',150, RTDB_ERR)
      end if

      goto 10


c
c  bond
c
 200  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_f(sdist(1))) then
      if (inp_f(sdist(2))) then
         if (spline) then
            nkb = nkb + 1
            rtdb_name = 'metadynamics_nbond'
            value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkb)
            rtdb_name = 
     >          'metadynamics_bond_index:'//c_index_name(nkb)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,2,sindx)
            rtdb_name = 
     >          'metadynamics_bond_ksigma:'//c_index_name(nkb)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
            if (.not. value) call errquit(
     >      'metadynamics_input:failed writing bond',0,RTDB_ERR)
         else
            nkb_gs = nkb_gs + 1
            rtdb_name = 'metadynamics_gs_nbond'
            value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkb_gs)
            rtdb_name = 
     >          'metadynamics_gs_bond_index:'//c_index_name(nkb_gs)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,2,sindx)
            rtdb_name = 
     >          'metadynamics_gs_bond_ksigma:'//c_index_name(nkb_gs)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
            if (.not. value) call errquit(
     >      'metadynamics_input:failed writing bond',1,RTDB_ERR)
         end if
      end if
      end if
      end if
      end if

      goto 10

c
c  angle
c
 300  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then
      if (inp_f(sdist(1))) then
      if (inp_f(sdist(2))) then
         if (spline) then
            nka = nka + 1
            rtdb_name = 'metadynamics_nangle'
            value = rtdb_put(rtdb,rtdb_name,mt_int,1,nka)
            rtdb_name = 
     >          'metadynamics_angle_index:'//c_index_name(nka)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,3,sindx)
            rtdb_name = 
     >          'metadynamics_angle_ksigma:'//c_index_name(nka)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
            if (.not. value) call errquit(
     >      'metadynamics_input:failed writing angle',0, RTDB_ERR)
         else
            nka_gs = nka_gs + 1
            rtdb_name = 'metadynamics_gs_nangle'
            value = rtdb_put(rtdb,rtdb_name,mt_int,1,nka_gs)
            rtdb_name = 
     >          'metadynamics_gs_angle_index:'//c_index_name(nka_gs)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,3,sindx)
            rtdb_name = 
     >          'metadynamics_gs_angle_ksigma:'//c_index_name(nka_gs)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
            if (.not. value) call errquit(
     >         'metadynamics_input:failed writing angle',1,RTDB_ERR)
         end if
      end if
      end if
      end if
      end if
      end if

      goto 10

c
c  dihedral
c
 400  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then
      if (inp_i(sindx(4))) then
      if (inp_f(sdist(1))) then
      if (inp_f(sdist(2))) then
         if (spline) then
            nkd = nkd + 1
            rtdb_name = 'metadynamics_ndihedral'
            value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkd)
            rtdb_name =
     >          'metadynamics_dihedral_index:'//c_index_name(nkd)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,4,sindx)
            rtdb_name = 
     >          'metadynamics_dihedral_ksigma:'//c_index_name(nkd)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
            if (.not. value) call errquit(
     >         'metadynamics_input:failed writing dihedral',0,RTDB_ERR)
         else
            nkd_gs = nkd_gs + 1
            rtdb_name = 'metadynamics_gs_ndihedral'
            value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkd_gs)
            rtdb_name =
     >          'metadynamics_gs_dihedral_index:'//c_index_name(nkd_gs)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,4,sindx)
            rtdb_name = 
     >          'metadynamics_gs_dihedral_ksigma:'//c_index_name(nkd_gs)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
            if (.not. value) call errquit(
     >         'metadynamics_input:failed writing dihedral',1,RTDB_ERR)
         end if
      end if
      end if
      end if
      end if
      end if
      end if

      goto 10

c
c  coord_number
c
 500  if (.not.spline) goto 600
      nkc = nkc + 1
      rtdb_name = 'metadynamics_ncoordnum'
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkc)
      sdist(1) = -1.0d0
      sdist(2) = -1.0d0
      sdist(3) = -1.0d0
      sdist(4) = -1.0d0
      sdist(5) = -1.0d0
      sdist(6) = -1.0d0
      sdist(7) = -1.0d0
      sdist(8) = -1.0d0
      sdist(9) = -1.0d0
      sdist(10) = -1.0d0
      sdist(11) = -1.0d0
      n1 = 0
      n2 = 0
 501  if (.not. inp_a(test)) goto 599
      if (.not. inp_match(num_dirs2,.false.,test,dirs2,ind2))
     >  call errquit(
     >  'metadynamics_input: unknown directive',0,INPUT_ERR)

      goto ( 510,520,530,540,550,560,570,580,590,591,592,595 ) ind2
      call errquit(
     >     'metadynamics_input: unimplemented directive',ind2,INPUT_ERR)
c
c index1
c
 510   if (.not.MA_push_get(mt_int,nw_max_atom,'indx_start1',
     >                    index_start1(2),index_start1(1)))
     >    call errquit(
     > 'metadynamics_input:failed allocating index_start1',0,MA_ERR)
      n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do
      rtdb_name='metadynamics_coordnum_indxsize1:'//c_index_name(nkc)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,n1)
      rtdb_name ='metadynamics_coordnum_indx1:'//c_index_name(nkc)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,n1,
     >                              int_mb(index_start1(1)))
      value = value.and.MA_pop_stack(index_start1(2))
      if (.not. value) call errquit(
     > 'metadynamics_input:failed deallocating index_start1',0,MA_ERR)

      goto 501
c
c index2 
c
 520   if (.not.MA_push_get(mt_int,nw_max_atom,'indx_start2',
     >                    index_start2(2),index_start2(1)))
     >    call errquit(
     > 'metadynamics_input:failed allocating index_start2',0,MA_ERR)
      n2 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start2(1)+n2) = j
          n2 = n2+1
         end do
      end do
      rtdb_name='metadynamics_coordnum_indxsize2:'//c_index_name(nkc)
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,n2)
      rtdb_name ='metadynamics_coordnum_indx2:'//c_index_name(nkc)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,n2,
     >                              int_mb(index_start2(1)))
      value = value.and.MA_pop_stack(index_start2(2))
      if (.not. value) call errquit(
     > 'metadynamics_input:failed deallocating index_start2',0,MA_ERR)

      goto 501
c
c n
c
 530  if (.not.inp_f(sdist(1))) sdist(1)=-1.0d0
      goto 501
c
c m
c
 540  if (.not.inp_f(sdist(2))) sdist(2)=-1.0d0
      goto 501
c
c r0
c
 550  if (.not.inp_f(sdist(3))) sdist(3)=-1.0d0
      goto 501
c
c w
c
 560  if (.not.inp_f(sdist(4))) sdist(4)=-1.0d0
      goto 501
c
c sigma
c
 570  if (.not.inp_f(sdist(5))) sdist(5)=-1.0d0
      goto 501
c
c maxcoord
c
 580  if (.not.inp_f(sdist(6))) sdist(6)=-1.0d0
      goto 501
c
c range
c
 590  if (.not.inp_f(sdist(7))) sdist(7)=-1.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=-1.0d0
      goto 501
c
c wab
c
 591  if (.not.inp_f(sdist(9))) sdist(9)=-1.0d0
      goto 501
c
c sigmaab
c
 592  if (.not.inp_f(sdist(10))) sdist(10)=-1.0d0
      goto 501
c
c sprik
c
 595  sdist(11)=1.0d0
      goto 501

 599  if ((n1.le.0).or.(n2.le.0)) then
         nkc = nkc - 1
         rtdb_name = 'metadynamics_ncoordnum'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkc)
      else
         rtdb_name ='metadynamics_coordnum_ksigma:'//c_index_name(nkc)
         if (.not.rtdb_put(rtdb,rtdb_name,mt_dbl,11,sdist))
     >    call errquit('metadynamics_input:failed writing coord_number',
     >                 0,RTDB_ERR)
      end if

      goto 10

c
c  coord_number_gs
c
 600  nkc_gs = nkc_gs + 1
      rtdb_name = 'metadynamics_gs_ncoordnum'
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkc_gs)
      sdist(1) = -1.0d0
      sdist(2) = -1.0d0
      sdist(3) = -1.0d0
      sdist(4) = -1.0d0
      sdist(5) = -1.0d0
      sdist(6) = -1.0d0
      sdist(7) = -1.0d0
      sdist(8) = -1.0d0
      sdist(9) = -1.0d0
      n1 = 0
      n2 = 0
 601  if (.not. inp_a(test)) goto 699
      if (.not. inp_match(num_dirs2,.false.,test,dirs2,ind2))
     >  call errquit(
     >  'metadynamics_input: unknown directive',0,INPUT_ERR)

      goto ( 610,620,630,640,650,660,670,680,690,695 ) ind2
      call errquit(
     >     'metadynamics_input: unimplemented directive',ind2,INPUT_ERR)
c
c index1
c
 610   if (.not.MA_push_get(mt_int,nw_max_atom,'indx_start1',
     >                    index_start1(2),index_start1(1)))
     >    call errquit(
     > 'metadynamics_input:failed allocating index_start1',0,MA_ERR)
      n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do
      rtdb_name=
     >'metadynamics_gs_coordnum_indxsize1:'//c_index_name(nkc_gs)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,n1)
      rtdb_name =
     >'metadynamics_gs_coordnum_indx1:'//c_index_name(nkc_gs)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,n1,
     >                              int_mb(index_start1(1)))
      value = value.and.MA_pop_stack(index_start1(2))
      if (.not. value) call errquit(
     > 'metadynamics_input:failed deallocating index_start1',0,MA_ERR)

      goto 601
c
c index2 
c
 620   if (.not.MA_push_get(mt_int,nw_max_atom,'indx_start2',
     >                    index_start2(2),index_start2(1)))
     >    call errquit(
     > 'metadynamics_input:failed allocating index_start2',0,MA_ERR)
      n2 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start2(1)+n2) = j
          n2 = n2+1
         end do
      end do
      rtdb_name=
     >'metadynamics_gs_coordnum_indxsize2:'//c_index_name(nkc_gs)
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,n2)
      rtdb_name=
     >'metadynamics_gs_coordnum_indx2:'//c_index_name(nkc_gs)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,n2,
     >                              int_mb(index_start2(1)))
      value = value.and.MA_pop_stack(index_start2(2))
      if (.not. value) call errquit(
     > 'metadynamics_input:failed deallocating index_start2',0,MA_ERR)

      goto 601
c
c n
c
 630  if (.not.inp_f(sdist(1))) sdist(1)=-1.0d0
      goto 601
c
c m
c
 640  if (.not.inp_f(sdist(2))) sdist(2)=-1.0d0
      goto 601
c
c r0
c
 650  if (.not.inp_f(sdist(3))) sdist(3)=-1.0d0
      goto 601
c
c w
c
 660  if (.not.inp_f(sdist(4))) sdist(4)=-1.0d0
      goto 601
c
c sigma
c
 670  if (.not.inp_f(sdist(5))) sdist(5)=-1.0d0
      goto 601
c
c maxcoord
c
 680  if (.not.inp_f(sdist(6))) sdist(6)=-1.0d0
      goto 601
c
c range
c
 690  if (.not.inp_f(sdist(7))) sdist(7)=-1.0d0
      if (.not.inp_f(sdist(8))) sdist(8)=-1.0d0
      goto 601
c
c sprik
c
 695  sdist(9)=1.0d0
      goto 601

 699  if ((n1.le.0).or.(n2.le.0)) then
         nkc_gs = nkc_gs - 1
         rtdb_name = 'metadynamics_gs_ncoordnum'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkc_gs)
      else
        rtdb_name=
     >  'metadynamics_gs_coordnum_ksigma:'//c_index_name(nkc_gs)
         if (.not.rtdb_put(rtdb,rtdb_name,mt_dbl,9,sdist))
     >    call errquit('metadynamics_input:failed writing coord_number',
     >                 0,RTDB_ERR)
      end if

      goto 10


c
c  spline
c
 700  spline = .true.
      if (inp_a(test)) then
         if (inp_compare(.false.,test,'on'))      spline = .true.
         if (inp_compare(.false.,test,'yes'))     spline = .true.
         if (inp_compare(.false.,test,'.true.'))  spline = .true.
         if (inp_compare(.false.,test,'T'))       spline = .true.
         if (inp_compare(.false.,test,'off'))     spline = .false.
         if (inp_compare(.false.,test,'no'))      spline = .false.
         if (inp_compare(.false.,test,'.false.')) spline = .false.
         if (inp_compare(.false.,test,'F'))       spline = .false.
      end if

      goto 10


 9999 continue

      return
      end
