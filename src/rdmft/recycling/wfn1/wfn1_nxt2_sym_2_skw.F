C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Transform a symmetric matrix into a skew symmetric matrix
C>
C> This is a rather trivial operation. We simply flip the sign on
C> the upper triangle and zero the diagonal.
C>
      subroutine wfn1_nxt2_sym_2_skw(nbf,q)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision q(nbf,nbf) !< [In/Output] The matrix
c
      integer i,j !< Counters
c
      do j = 1, nbf
        do i = 1, j-1
          q(i,j) = -q(i,j)
        enddo
      enddo
      do i = 1, nbf
        q(i,i) = 0.0d0
      enddo
c
      end
C> @}
