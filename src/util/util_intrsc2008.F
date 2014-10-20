
      integer function util_popcnt(i)
      implicit none
      integer i,j,s

      util_popcnt = count((/ (BTEST(i,j),j=0,bit_size(i)-i) /))
      return
      end

      integer function util_trailz(i)
      implicit none
      integer i
      logical done
      integer s,k

      s = 1
      done = .false.
      do while (.not.done)
         k = ishft(ishft(i,-s),s)
         if ((k.eq.i).and.(s.le.64)) then
            s = s + 1
         else
            done = .true.
         end if
      end do
      util_trailz = s-1
      return
      end 

      integer function util_leadz(i)
      implicit none
      integer i
      logical done
      integer s,k

      s = 1
      done = .false.
      do while (.not.done)
         k = ishft(ishft(i,s),-s)
         if ((k.eq.i).and.(s.le.64)) then
            s = s + 1
         else
            done = .true.
         end if
      end do
      util_leadz = s-1
      return
      end

