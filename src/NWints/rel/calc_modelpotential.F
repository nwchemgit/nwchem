c
c     Calculate the model potential
c     
c Molecular density functional calculations in the regular relativistic approximation:
c Method, application to coinage metal diatomics, hydrides, fluorides and chlorides, and
c comparison with first-order relativistic calculations
c Christoph van Wüllen
c The Journal of Chemical Physics 109, 392 (1998); doi: 10.1063/1.476576
c
c Accurate and efficient treatment of two-electron contributions in 
c quasirelativistic high-order Douglas-Kroll density-functional calculations
c Christoph van Wüllen and Christine Michauk
c J. Chem. Phys. 123, 204113 (2005); https://doi.org/10.1063/1.2133731
c
c
      subroutine calc_modelpotential(geom,natoms,nqpts,qxyz,qwght,
     &                      closegridpts,gexpo,gcoef,amat_mp)
c
      implicit none
c
#include "geom.fh"
#include "stdio.fh"
#include "zora.fh"
c
      integer geom
      integer natoms
      integer nqpts
      double precision qxyz(3,nqpts),qwght(nqpts)
      integer closegridpts(*)
      double precision gexpo(natoms,50),gcoef(natoms,50)
      double precision amat_mp(nqpts)
c
      logical lSuccess
      integer i,j
      character*16 tags(natoms)
      double precision nucCharge, nucCoords(3)
      double precision rx,ry,rz,dist,pi,rtemp,rtemp2
c
      integer igrid,igau
      double precision atom_mp
c
      double precision util_erf
      external util_erf
c
      double precision dgami
      external dgami  ! Evaluates Incomplete Gamma Function P(a,x)
c
      pi=dacos(-1d0)
c
c     get the total zora potential of all atoms at grid points
      do igrid = 1,nqpts 
       amat_mp(igrid) = 0.d0
c
       do i = 1,natoms
        lSuccess = geom_cent_get(geom,i,tags,nucCoords,nucCharge)

c       distance from the grid points to the atom centers
        rx = nucCoords(1) - qxyz(1,igrid)
        ry = nucCoords(2) - qxyz(2,igrid)
        rz = nucCoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
c
c       tally up zora potential at grid point
        if (dist.gt.zoracutoff) then 
          atom_mp = 0.d0
          do j = 1,50
           rtemp = dsqrt(gexpo(i,j))*dist
           rtemp2 = rtemp*rtemp
c           atom_mp = atom_mp + gcoef(i,j)*util_erf(rtemp)
           atom_mp = atom_mp + gcoef(i,j)*dgami(0.5d0,rtemp2)/dsqrt(pi)
          end do ! j
          amat_mp(igrid) = amat_mp(igrid) + atom_mp/dist
        else
          closegridpts(igrid) = igrid
        end if
       end do ! i
      end do ! igrid
c
      return 
      end
