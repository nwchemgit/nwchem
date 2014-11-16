      integer function inp_strlen(a)
      implicit none
      character*(*) a !< [Input] The string
      integer i
      logical ois_ws
      intrinsic len
      character*1 xtest
      character*1 xblnk
      character*1 xtab
      parameter(xblnk = ' ')
      parameter(xtab = '	')
      ois_ws(xtest) = (xtest.eq.xblnk .or. xtest.eq.xtab)
c
      do i = len(a),1,-1
         if (.not. ois_ws(a(i:i))) goto 10
      enddo
c
 10   inp_strlen = i
c
      end

