      block data inp_data
C$Id$
      implicit none
#include "inpP.fh"
c
      data iread   /5/
      data iwrite  /6/
      data jrec    /-1/
      data jump    /0/
      data oswit   /.false. /
      data nerr    /999/
      data nline   /0/
      data noline  /0/
      data ierrpos /-1/
      data errmsg  /' '/
      data input_line /0/
      data xblnk /' '/
#if defined(CRAY) || defined(HPUX) || defined(WIN32) || defined(PSCALE) || defined(__FLANG) || ( __GNUC__ >= 4)
      data xtab  /'	'/      ! Tab ... no backslash necessary
#elif (defined(LINUX) || defined(MACX)) && !defined(PGLINUX) && !defined(XLFLINUX) &&!( __GNUC__ >= 4)
      data xtab  /9/            ! Tab ... g77 has trouble with escape sequence
#else
      data xtab  /'\	'/      ! Tab ... note backslash for cpp
#endif
      data xsplit/';'/
      data xcomm /'#'/
      data xback /'\\'/         ! Backslash ... note backslash for cpp
      data xquote/'"'/
c
      data save_level /0/            ! Depth of nesting in save/restore
      data include_level /0/
c
      end
c
C> \brief Returns the length of the contents of a string
C>
C> In Fortran the length of a string is determined by the amount of 
C> memory allocated for it. Often one just wants to know how long the
C> value of the string is without any trailing spaces. Which is what
C> this function returns.
C>
C> \return The length of the value of a string without any trailing 
C> spaces
c
      integer function inp_strlen(a)
      implicit none
#include "inpP.fh"
      character*(*) a !< [Input] The string
      integer i
      logical ois_ws
      intrinsic len
      character*1 xtest
      ois_ws(xtest) = (xtest.eq.xblnk .or. xtest.eq.xtab)
c
      do i = len(a),1,-1
         if (.not. ois_ws(a(i:i))) goto 10
      enddo
c
 10   inp_strlen = i
c
      end
