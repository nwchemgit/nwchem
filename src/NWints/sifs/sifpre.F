      subroutine sifpre( nlist, nenrgy, energy, ietype )
c
c  print out the energy(*) array with character labels.
c
c  24-apr-92 cnvginf values added. -rls
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer  nlist,  nenrgy, ietype(*)
      real*8   energy(*)
c
      integer  i,      itypea, itypeb
      character*8      chrtyp
c
      do 10 i = 1, nenrgy
c
         itypea = ietype(i) / 1024
         itypeb = mod( ietype(i), 1024 )
         call siftyp( itypea, itypeb, chrtyp )
c
         if ( ietype(i) .ge. 0 ) then
c
c           # fcore from a 1-e hamiltonian array.
c
            write(nlist,6010) i, energy(i), ietype(i), 'fcore', chrtyp
c
         elseif ( itypea .eq. 0 ) then
c
c           # core energy value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'core', chrtyp
c
         elseif ( itypea .eq. -1 ) then
c
c           # total energy value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'total', chrtyp
c
         elseif ( itypea .eq. -2 ) then
c
c           # energy or wave function convergence value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'cnvginf', chrtyp
c
         else
c
c           # undefined energy type.
c
            write(nlist,6010) i, energy(i), ietype(i), 'unknown', chrtyp
c
         endif
10    continue
c
      return
6010  format(' energy(', i2, ')=', 1pe20.12, ', ietype=', i5,
     & ', ', a7, ' energy of type: ', a )
      end
