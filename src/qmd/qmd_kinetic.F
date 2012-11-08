c
c     calculate kinetic energy 
c
      subroutine qmd_kinetic(nat,mass,v,ekin)
c
      implicit none
c
      integer nat
      double precision mass(nat)
      double precision v(3,nat)
      double precision ekin
c
      integer i,j
c
c     total kinetic energy
      ekin = 0.d0
      do i = 1,nat
        do j = 1,3
            ekin = ekin + v(j,i)*v(j,i)
        end do
        ekin = mass(i)*ekin
      end do
      ekin = 0.5d0*ekin
c
      return
      end
