c
c     calculate kinetic energy 
c
      subroutine qmd_kinetic(nat,m,v,ekin)
c
      implicit none
c
      integer nat
      double precision m(nat)
      double precision v(3,nat)
      double precision ekin
c
      integer i
c
c     total kinetic energy
      ekin = 0.d0
      do i = 1,nat
       ekin = ekin +  m(i)*(v(1,i)*v(1,i)+v(2,i)*v(2,i)+v(3,i)*v(3,i))
      end do
      ekin = 0.5d0*ekin
c
      return
      end
