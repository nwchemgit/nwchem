c
c $Id: util_io.F,v 1.1 2007-12-11 19:09:14 marat Exp $
c

      function util_get_io_unit(fn)

      implicit none
      integer fn
      logical util_get_io_unit
c 
      integer k
      logical ostatus
c
      do k=80,90
        INQUIRE(UNIT=k,OPENED=ostatus)
        ostatus = .not.ostatus
        if(ostatus) 
     >    INQUIRE(UNIT=k,EXIST=ostatus)
        if(ostatus) then
          fn = k
          util_get_io_unit = .true.
          return
        end if 
      end do
      util_get_io_unit = .false.
      return
      end


