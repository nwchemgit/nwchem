c
c     Constrain molecular dynamics within a cavity
c
      subroutine qmd_cavity(natom,cavity_axis,cavity_spring,cent,r,g)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
c
      integer natom
      double precision cavity_axis(3)   ! a,b,c axes of the ellipsoid
      double precision cavity_spring    ! cavity spring
      double precision cent(3)          ! centroid of molecule
      double precision r(3,natom)       ! coordinates of molecule
      double precision g(3,natom)       ! gradients
c
      integer i,j
      double precision a2, b2, c2
      double precision d_r
      double precision cavity(3),rc(3)
      double precision d_atom, d_cavity, d_atom_cavity
      double precision cavity_force, normal(3), norm
c
c     define an cavity (if the axes lengths are the same, we have a sphere)
      a2= cavity_axis(1)*cavity_axis(1)
      b2= cavity_axis(2)*cavity_axis(2)
      c2= cavity_axis(3)*cavity_axis(3)
c 
c     loop over atoms and determine if the atoms need to be contained in the cavity
      do i=1,natom
c
        rc(1) = r(1,i)-cent(1)    ! with respect to the centroid
        rc(2) = r(2,i)-cent(2)
        rc(3) = r(3,i)-cent(3)

        d_r = (rc(1)*rc(1)/a2) + (rc(2)*rc(2)/b2) + (rc(3)*rc(3)/c2)
c      
        if (d_r > 1.d0) then      ! outside cavity
c         
c         calculate intercept with cavity wall
          cavity(1) = rc(1)/dsqrt(d_r)
          cavity(2) = rc(2)/dsqrt(d_r)
          cavity(3) = rc(3)/dsqrt(d_r)
c               
c         calculate distance of the atom and cavity wall from the center of the cavity
          d_atom = dsqrt(rc(1)**2 + rc(2)**2 + rc(3)**2)
          d_cavity = dsqrt(cavity(1)**2 + cavity(2)**2 + cavity(3)**2)
c         
c         calculate distance between cavity wall and atom
          d_atom_cavity = d_atom - d_cavity
c
c         calculate the magnitude of the reflection force from the cavity wall
          cavity_force = 0.5d0*cavity_spring*(d_atom_cavity**2) ! spring like bounce
c
c         calculate the normal to the cavity wall
          normal(1) = cavity(1)/a2 
          normal(2) = cavity(2)/b2 
          normal(3) = cavity(3)/c2 
c         
c         normalize the normal
          norm = normal(1)**2 + normal(2)**2 + normal(3)**2
          if (norm.ne.0.d0) then
            normal(1) = normal(1)/dsqrt(norm)
            normal(2) = normal(2)/dsqrt(norm)
            normal(3) = normal(3)/dsqrt(norm)
          end if
c         
c         add the cavity force (gradient) to the total gradient
          do j=1,3
            g(j,i) = g(j,i) - normal(j)*(-1.d0*cavity_force) ! -ve of force for gradient
          enddo
c         
        end if  ! outside barrier
c
      enddo  ! natom
c
      return
      end
