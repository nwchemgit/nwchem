*----------------------------------------------------------------------*
      subroutine fh_init(ioccun,noccun,iprint)
*----------------------------------------------------------------------*

      implicit none
      include "freeunits.inc"
      
      integer, intent(in) ::
     &     noccun,            ! number of units to be reserved
     &     ioccun(noccun),    ! unit numbers to be reserved
     &     iprint             ! print level

      integer :: iunit
      logical :: lopen

      ifrunit(1:6)       = 3 ! reserved
      ifrunit(7:mxpunit) = 0
      do iunit = 1, noccun
        ifrunit(ioccun(iunit)) = 2 ! reserved by user
      end do
      do iunit = 7, mxpunit
        if (ifrunit(iunit).eq.0) then
          inquire(unit=iunit,opened=lopen)
          if (lopen) then
            ifrunit(iunit) = -2
          end if
        end if
      end do

      if (iprint.gt.10) then
        call fh_prstat(6)
      end if

      end

*----------------------------------------------------------------------*
      subroutine fh_prstat(luout)
*----------------------------------------------------------------------*
      implicit none
      include "freeunits.inc"

      integer, intent(in) ::
     &     luout         ! output unit

      integer ::
     &     iunit
      character ::
     &     cstat*20

      write(6,'(2(/x,a))')
     &       'Status of unit numbers',
     &       '======================'
      do iunit = 1, mxpunit
        if (ifrunit(iunit).eq.0)  cstat = 'FREE                '
        if (ifrunit(iunit).eq.3)  cstat = 'RESERVED BY SYSTEM  '
        if (ifrunit(iunit).eq.2)  cstat = 'RESERVED BY USER    '
        if (ifrunit(iunit).eq.-2) cstat = 'OCCUPIED ON INIT    '
        if (ifrunit(iunit).eq.1)  cstat = 'RESERVED BY IGETUN  '
        if (ifrunit(iunit).eq.-1) cstat = 'OCCUPIED AFTER INIT '
        write(luout,'(x,i3,2x,a)') iunit, cstat
      end do

      end
*----------------------------------------------------------------------*
      integer function igetunit(iwish)
*----------------------------------------------------------------------*
      
      implicit none
      include "freeunits.inc"

      integer, intent(in) ::
     &     iwish

      integer :: iunit
      logical :: lopen

      ! trying to fulfil your wishes (if you asked for that)
      if (iwish.ge.7) then
        iunit = iwish
        if (ifrunit(iwish).eq.0) then
          inquire(unit=iunit,opened=lopen)
          if (lopen) then
            ifrunit(iunit) = -1
          else
            ifrunit(iunit) = 1
            igetunit = iunit
            return
          end if
        end if
      end if

      ! looking for the lowest free unit
      do iunit = 7, mxpunit
        if (ifrunit(iunit).eq.0) then
          inquire(unit=iunit,opened=lopen)
          if (lopen) then
            ifrunit(iunit) = -1
          else
            ifrunit(iunit) = 1
            igetunit = iunit
            return
          end if
        end if
      end do

      write(6,*) 'Fatal: run out of unit numbers in igetun'
      call fh_prstat(6)
      stop 'run out of unit numbers'

      end
*----------------------------------------------------------------------*
      subroutine relunit(iunit,cstat)
*----------------------------------------------------------------------*

      implicit none
      include "freeunits.inc"

      integer, intent(in) :: iunit
      character, intent(in) :: cstat*(*)

      logical :: lopen, lnamed
      integer :: len
      character :: filnam*256

      if (iunit.le.mxpunit.and.iunit.gt.0.and.ifrunit(iunit).eq.1) then
        inquire(unit=iunit,opened=lopen)
        if (lopen) close(iunit,status=cstat)
        ifrunit(iunit) = 0
      else
        write(6,*) 'relunit called for illegal unit ',
     &       iunit, ifrunit(iunit)
        inquire(unit=iunit,opened=lopen,named=lnamed,name=filnam)
        if (lopen) then
          write(6,*) ' the file was open'
        else
          write(6,*) ' the file was not open'
        end if
        if (lnamed) then
          len = len_trim(filnam)
          write(6,*) ' the file was named: "',filnam(1:len),'"'
        else
          write(6,*) ' the file was not named'
        end if
        stop 'illegal unit in relunit'
      end if
      
      end

*----------------------------------------------------------------------*
      integer function iopen_uus()
*----------------------------------------------------------------------*
* open  unnamed unformatted sequential

      implicit none
      include "freeunits.inc"

      integer :: iunit
      integer, external :: igetunit

      iunit = igetunit(-1)
      open(unit=iunit,form="unformatted",access="sequential")

      iopen_uus = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nus(filename)
*----------------------------------------------------------------------*
* open  named unformatted sequential

      implicit none
      include "freeunits.inc"

      character, intent(in) :: filename*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="sequential")

      iopen_nus = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nfs(filename)
*----------------------------------------------------------------------*
* open  named formatted sequential

      implicit none
      include "freeunits.inc"

      character, intent(in) :: filename*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="formatted",access="sequential")

      iopen_nfs = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nuss(filename,statstr)
*----------------------------------------------------------------------*
* open  named unformatted sequential in status 'statstr'

      implicit none
      include "freeunits.inc"

      character, intent(in) :: filename*(*), statstr*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="sequential",status=statstr)

      iopen_nuss = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nud(filename,lrec)
*----------------------------------------------------------------------*
* open  named unformatted direct access

      implicit none
      include "freeunits.inc"

      character, intent(in) :: filename*(*)
      integer, intent(in)   :: lrec

      integer :: iunit, len
      integer, external :: igetunit


      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="direct",recl=lrec)

      iopen_nud = iunit

      end

*----------------------------------------------------------------------*
      subroutine unit_info(lu)
*----------------------------------------------------------------------*
      implicit none
      include "freeunits.inc"
      
      integer, intent(in) :: lu

      logical :: lnamed, lopen
      integer :: len, lrec, nrec
      character :: str*256

      inquire(unit=lu,named=lnamed,opened=lopen,name=str)

      write(6,'(/x,78("="))')
      write(6,*) ' Information on unit ',lu
      if (lopen) then
        write(6,*) '  file is open'
        if (lnamed) then
          len = len_trim(str)
          write(6,*) '  filename: ',str(1:len)
        else
          write(6,*) '  unnamed'
        end if
        inquire(unit=lu,form=str)
        len = len_trim(str)
        write(6,*) '  format: ',str(1:len)
        inquire(unit=lu,access=str)
        len = len_trim(str)
        write(6,*) '  access: ',str(1:len)
        if (str(1:len).eq.'DIRECT'.or.str(1:len).eq.'direct') then
          inquire(unit=lu,recl=lrec,nextrec=nrec)
          write(6,*) '  record-length: ',lrec
          write(6,*) '  next record  : ',nrec
        end if
      else
        write(6,*) ' file is not open'
      end if
      write(6,'(x,78("="))')

      end
c $Id$
