C> \brief Find an available file unit
C>
C> This function scans file unit numbers in the range given by 
C> `io_min` to `io_max`. If an available unit number is encountered
C> the function returns `.TRUE.` and the unit number if given in `fn`.
C> Otherwise the function returns `.FALSE.` and `fn` remains unchanged.
C>
C> \return `.TRUE.` if an available unit number was found and the
C> unit number is returned in `fn`, otherwise `.FALSE.` is returned
C>
      function util_io_unit(io_min,io_max,fn)

      implicit none
      integer io_min !< [Input] The minimum unit number
      integer io_max !< [Input] The maximum unit number
      integer fn     !< [Output] The first available unit number
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

c $Id$
