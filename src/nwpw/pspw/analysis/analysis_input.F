      subroutine analysis_input(rtdb)
c
c $Id: analysis_input.F,v 1.1 2001-08-30 00:25:09 edo Exp $
c
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 3)
      integer lmax,lorb,i,nrho
      real*8  b0(0:5),u2(0:5),u3(0:5)
      real*8  x2,y2,x3,y3,dum
      character*2 atom
      logical     lorb_read
 
      character*30 dirs(num_dirs)
      character*255 test,id,full_filename

      data dirs / 'psp_filename:',
     >            'max_orbital:',
     >            'end'/

      lorb_read = .false.

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'analysis_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >      'analysis_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'analysis_input: unknown directive', 0)


      goto ( 100, 200,
     >      9999) ind
      call errquit(
     >   'analysis_input: unimplemented directive', ind)


c
c  psp_filename
c
 100  if (.not. inp_a(test))
     >     call errquit(
     >          'analysis_input: failed to read psp_filename', 0)
      call util_file_name_noprefix(test,.false.,.false.,
     >                    full_filename)

      ind = index(full_filename,' ') - 1
      open(82,FILE=full_filename,STATUS='OLD',FORM='FORMATTED')
        read(82,'(A2)') atom
        read(82,*) dum,dum,lmax
        read(82,*)
        read(82,*) nrho,dum
        do i=1,nrho
          read(82,*)
        end do
        read(82,*) 
        read(82,*) x2,(u2(i),i=0,lmax)
        read(82,*) x3,(u3(i),i=0,lmax)
      close(82)
      if (.not.lorb_read) lorb = lmax -1
      do i=0,lorb
          Y2=U2(i)/X2**(i+1)
          Y3=U3(i)/X3**(i+1)
          B0(i)=Y2-X2*(Y2-Y3)/(X2-X3)
      end do
      WRITE(6,*)
      WRITE(6,1300) test(1:ind)
      WRITE(6,1305) 'ATOM ','S','P','D','F'
      WRITE(6,1310) ATOM,(B0(i),i=0,lorb)
      WRITE(6,*)
      WRITE(6,*)
      call util_flush(6)
 1300 FORMAT('== ',A,' expansion coefficients =='/)
 1305 FORMAT(A5,6X,A,12X,A,12X,A,12X,A)
 1310 FORMAT(A2,' : ',4E13.5)

      id = 'analysis:lorb'//atom
      if (.not. rtdb_put(rtdb,id,mt_int,1,lorb))
     >        call errquit(
     >        'analysis_input: rtdb_put failed', 0)

      id = 'analysis:expansion'//atom
      if (.not. rtdb_put(rtdb,id,mt_dbl,(lorb+1),b0))
     >        call errquit(
     >        'analysis_input: rtdb_put failed', 0)

      lorb_read = .false.
      goto 10
 
c
c  lorb
c
 200  if (.not. inp_i(lorb))
     >        call errquit(
     >        'analysis_input: lorb read failed', 0)
      lorb_read = .true.
      goto 10

 9999 continue

      return
      end
