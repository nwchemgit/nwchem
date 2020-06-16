      integer function selci_popcnt32(i)
      integer i
      selci_popcnt32 = popcnt(i)
      end

      subroutine selci_vpop32(n, input, output)
      integer n
      integer input(*)
      integer output(*)
      integer i
      do i = 1,n
        output(i) = popcnt(input(i))
      enddo
      end

      integer function selci_leadz32(i)
      integer i
      integer function util_leadz
      selci_leadz32 = leadz(i)
      end
