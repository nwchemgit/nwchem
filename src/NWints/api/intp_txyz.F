C> \ingroup nwint
C> @{
C>
C> \brief Generate periodic images of atomic centers
C>
C> See [1] for details.
C>
C> [1] JE Jaffe, AC Hess,
C>     <i>"Gaussian basis density functional theory for systems 
C>     periodic in two or three dimensions: Energy and forces"</i>,
C>    J.Chem.Phys. <b>105</b>, 10983-10998 (1996), DOI:
C>    <a href="https://doi.org/10.1063/1.472866">
C>    10.1063/1.472866</a>
C>
      subroutine intp_txyz(n_cent, n_geom, R, xyz_new)
c $Id$
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c::external subroutines used
c... dgemv
c::passed
      integer n_cent               !< [Input] center to be translated
      integer n_geom               !< [Input] geometry handle to use
      double precision R(3)        !< [Input] fractional translation vector
      double precision xyz_new(3)  !< [Output] new coordinates based on R translations
c::local
      double precision xyz_frac(3) !< fractional coordinates 
      integer ixyz
c
c.. zero local arrays
      call dcopy(3,0.0d00,0,xyz_new,1)
      call dcopy(3,0.0d00,0,xyz_frac,1)
c
c... compute n center fractional coordinates
c....    xyz_frac = amatrix_inv*coords(*,n_cent,n_geom)
c
      call dgemv('n',3,3,1.0d00,amatrix_inv(1,1,n_geom),3,
     &    coords(1,n_cent,n_geom),1,
     &    0.0d00,xyz_frac,1)
c
c... compute j translated fractional coordinates
c.........  xyz_frac = xyz_frac + R
c
      do 00100 ixyz=1,3
          xyz_frac(ixyz) = xyz_frac(ixyz) + R(ixyz)
00100 continue
c
c... computer j translated cart. coordinates.
c........ xyz_new = amatrix*xyz_frac
c
      call dgemv('n',3,3,1.0d00,amatrix(1,1,n_geom),3,
     &    xyz_frac,1,
     &    0.0d00,xyz_new,1)
c
      end
C> @}
