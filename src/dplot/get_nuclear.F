!>    @brief express the nuclear charges on a grid
!>
!>    Given the positions and nuclear charges of atoms express the 
!>    charge from the nucleii on a grid.
!>
!>    @param[in]  ngrid  The number of grid points
!>    @param[in]  natom  The number of atoms
!>    @param[in]  ri     The grid point positions
!>    @param[out] vl     The charge values on the grid
!>    @param[in]  z      The nuclear charges
!>    @param[in]  ra     The nuclear positions
!>    @param[in]  blur   The spread of the nuclear charges
!>
      subroutine get_nuclear(ngrid,natom,ri,vl,z,ra,blur)
      implicit none
#include "global.fh"
      integer,          intent(in)  :: ngrid
      integer,          intent(in)  :: natom
      double precision, intent(in)  :: ri(3,ngrid)
      double precision, intent(out) :: vl(ngrid)
      double precision, intent(in)  :: z(natom)
      double precision, intent(in)  :: ra(3,natom)
      double precision, intent(in)  :: blur
!
      double precision, external :: get_gaussian
!
      integer igrid, iatom
      integer nproc, iproc
!
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      vl = 0.0d0
      do igrid = 1, ngrid
        if (mod(igrid,nproc).eq.iproc) then
          do iatom = 1, natom
            vl(igrid) = vl(igrid)
     +                + get_gaussian(ra(1:3,iatom),ri(1:3,igrid),
     +                               blur,z(iatom))
          enddo
        endif
      enddo
      end
