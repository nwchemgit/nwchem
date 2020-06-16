      integer function selci_popcnt32(i)
      integer i
#ifdef NWCHEM_USE_FORTRAN2008
      selci_popcnt32 = popcnt(i)
#else
      selci_popcnt32 = util_popcnt(i)
#endif
      end

      subroutine selci_vpop32(n, input, output)
      integer n
      integer input(*)
      integer output(*)
      integer i
      do i = 1,n
#ifdef NWCHEM_USE_FORTRAN2008
        output(i) = popcnt(input(i))
#else
        output(i) = util_popcnt(input(i))
#endif
      enddo
      end

      integer function selci_leadz32(i)
      integer i
#ifdef NWCHEM_USE_FORTRAN2008
      selci_leadz32 = leadz(i)
#else
      selci_leadz32 = util_leadz(i)
#endif
      end
