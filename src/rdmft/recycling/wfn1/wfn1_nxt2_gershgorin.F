C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Compute the Gershgorin radii of a matrix
C>
C> For a given matrix compute the Gershgorin radius [1] for each and
C> every diagonal element.
C>
C> ### References ###
C>
C> [1] S. Gershgorin, "Über die Abgrenzung der Eigenwerte einer Matrix."
C>     Izv. Akad. Nauk. USSR Otd. Fiz.-Mat. Nauk <b>6</b>, 749–754, 1931
C>     <a href="http://mi.mathnet.ru/eng/izv/y1931/i6/p749">
C>     http://mi.mathnet.ru/eng/izv/y1931/i6/p749</a>
C>
      subroutine wfn1_nxt2_gershgorin(nbf,a,r)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision a(nbf,nbf) !< [Input] The matrix
      double precision r(nbf) !< [Output] The Gershgorin radii
c
      integer ir,ic
c
      do ir = 1, nbf
        r(ir) = -abs(a(ir,ir))
      enddo
      do ic = 1, nbf
        do ir = 1, nbf
          r(ir) = r(ir) + abs(a(ir,ic))
        enddo
      enddo
c
      end
C>
C> @}
      
