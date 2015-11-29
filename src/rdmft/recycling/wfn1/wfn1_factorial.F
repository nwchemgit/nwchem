C> \ingroup wfn1
C> @{
C>
C> \brief Computes the factorial of a (small) integer
C>
C> Compute the factorial of an integer using a simple algorithm.
C> We need this number to establish how many permutations of a given
C> number of things are possible. The input integers must be 
C> non-negative.
C>
      integer function wfn1_factorial(nn)
      implicit none
      integer nn   !< [Input] The integer to compute the factorial off
      integer i    !< Counter
      integer ires !< The result
c
      ires = 1
      do i = 2, nn
        ires = ires * i
      enddo
c
      wfn1_factorial = ires
c
      return
      end
C>
C> @}

