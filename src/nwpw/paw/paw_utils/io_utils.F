c
c $Id$
c

!*******************************************************
!
!   Name    :
!
!   Purpose :
!
!   Created :
!*******************************************************
      subroutine get_io_unit(fp)
      implicit none
      integer  fp
    
      !*** local variables ***
      integer  min_io_unit,max_io_unit 
      parameter(min_io_unit = 20,max_io_unit = 100)

      logical open_status,exist_status
      integer i,ierr

      DO i=min_io_unit,max_io_unit
      
        inquire(UNIT=i,
     >        OPENED = open_status,
     >        EXIST = exist_status,
     >        IOSTAT= ierr)
        if(.not. open_status .and.
     >    exist_status .and. ierr.eq.0) then
          fp = i
          return
        end if

      end do

      call errquit("cannot get free file unit",0,0)

      return
      end 
!*******************************************************
!
!   Name    :
!
!   Purpose :
!
!   Created :
!*******************************************************
      function assert_file(filename) 
      implicit none
      character*(*) filename
       
      logical assert_file !result

      !*** local variables ***
      logical  file_exist

      inquire(FILE=filename, EXIST=file_exist)
      assert_file = file_exist
      return
      end

!*******************************************************
!
!   Name    :
!
!   Purpose :
!
!   Created :
!*******************************************************
      subroutine open_file(FileName,mode,fp)
      implicit none
      character*(*) FileName
      character*(*) mode
      integer      fp

*     !**** local variables ****
      integer  l,ierr
      character*20 action, status,position
      character*255 junk

      l = index(mode,' ') -1
      if (l.lt.1) l = 1
      if (mode(1:l).eq.'r') then
          action   = 'READ'
          status   = 'OLD'
          position = 'REWIND'
      else if (mode(1:l).eq.'w') then
          action = 'WRITE'
          status = 'REPLACE'
          position = "REWIND"
      else if (mode(1:l).eq.'w+') then
          action   = 'WRITE'
          status   = 'OLD'
          position = "APPEND"
      else
          junk = "unknown parameter in open_file, "//mode(1:l)
          call errquit(junk,0,0)
      end if

      call get_io_unit(fp)
    
      l = index(FileName,' ') -1
      open(UNIT=fp,FILE=FileName(1:l),
     >    STATUS=status,
     >    FORM='formatted',
     >    IOSTAT=ierr)

      if (ierr.ne.0)  then
          junk = " opening "//FileName
          call errquit(junk,0,0)
      end if

      return
      end !SUBROUTINE open_file

