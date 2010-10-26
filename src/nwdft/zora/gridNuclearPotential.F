c
c     == Calculate the nuclear potential on the grid ==
      subroutine gridNuclearPotential(geom,natoms,nqpts,qxyz,qwght,
     &   closegridpts,amat_nucl)
c
      implicit none
c
#include "geom.fh"
#include "stdio.fh"
#include "zora.fh"
c
      integer iatom,natoms,geom,igrid,nqpts
      double precision nucCharge, nucCoords(3)
      double precision qxyz(3,nqpts),qwght(nqpts)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision amat_nucl(nqpts)
      integer closegridpts(*)
c
c     == get the total nuclear potential on the grid ==
      do igrid = 1,nqpts   ! == loop over the grid points ==
       amat_nucl(igrid) = 0.d0  
       do iatom = 1,natoms ! == loop over the atoms ==
c
c     == get an atom (needs error handling) ==
        lSuccess = geom_cent_get(geom, iatom, tags, nucCoords, 
     &             nucCharge)
c
c     == distance from the grid points to the atom centers == 
        rx = nucCoords(1) - qxyz(1,igrid)
        ry = nucCoords(2) - qxyz(2,igrid)
        rz = nucCoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff) then ! == check cutoff ==
          amat_nucl(igrid) = amat_nucl(igrid) - nucCharge/dist
        else
          closegridpts(igrid) = igrid
        end if
       end do
      end do
c
      return
      end 
