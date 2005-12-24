      function util_io_unit(io_min,io_max,fn)

      implicit none
      integer io_min,io_max
      integer fn
      logical util_io_unit
c 
      integer k
      logical ostatus
c
      do k=io_min,io_max
        INQUIRE(UNIT=k,OPENED=ostatus)
        ostatus = .not.ostatus
        if(ostatus) 
     >    INQUIRE(UNIT=k,EXIST=ostatus)
        if(ostatus) then
          fn = k
          util_io_unit = .true.
          return
        end if 
      end do
      util_io_unit = .false.
      return
      end

