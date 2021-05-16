C>
C> \brief Implementation of POPCNT
C>
C> POPCNT is a Fortran 2008 intrinsic. This is an explicit 
C> implementation because not all compilers support this function
C> yet. The result returned is the number of bits of value 1 in 
C> argument I.
C>
      integer function util_popcnt(i)
      implicit none
      integer i !< [Input] The value I
      integer j
      util_popcnt = count((/ (BTEST(i,j),j=0,bit_size(i)-1) /))
      return
      end
C>
C> \brief Implementation of TRAILZ
C>
C> TRAILZ is a Fortran 2008 intrinsic. This is an explicit 
C> implementation because not all compilers support this function
C> yet. The result returned is the number of trailing 0 bits in 
C> argument I. If I equals 0 then the trailing number of 0 bits
C> equals the total number of bits.
C>
      integer function util_trailz(i)
      implicit none
      integer i !< [Input] The value I
      integer k

      if (i.eq.0) then
         k = bit_size(i)
      else
         k = 0
         do while (.not.btest(i,k))
            k = k + 1
         end do
      end if
      util_trailz = k
      return
      end
C>
C> \brief Implementation of LEADZ
C>
C> LEADZ is a Fortran 2008 intrinsic. This is an explicit 
C> implementation because not all compilers support this function
C> yet. The result returned is the number of leading 0 bits in 
C> argument I. If I equals 0 then the leading number of 0 bits
C> equals the total number of bits.
C>
      integer function util_leadz(i)
      implicit none
      integer i !< [Input] The value I
      integer k, l

      if (i.eq.0) then
         k = bit_size(i)
      else
         l = bit_size(i)-1
         do while (.not.btest(i,l))
            l = l - 1
         end do
         k = bit_size(i)-1-l
      end if
      util_leadz = k
      return
      end

