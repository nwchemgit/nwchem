*
*     $Id$
*

*     ***************************
*     *                         *
*     *     zero_cell_contraint *
*     *                         *
*     ***************************
*
*   This routine adjusts the cell gradients so
* that specified constraints are satisfied.
*
*   strategy==1 - a,b,c,alpha,beta,gamma constraints
*   strategy==2 - lattice vector constraints
*
      subroutine zero_cell_constraint(stress,lstress)
      implicit none
      real*8 stress(3,3)
      real*8 lstress(6)

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      integer i,u,v

*     **** external functions ****
      integer  control_optimize_cell_strategy
      integer  control_optimize_lattice
      integer  control_optimize_lattice_vectors
      external control_optimize_cell_strategy
      external control_optimize_lattice
      external control_optimize_lattice_vectors

      call Parallel_taskid(taskid)

*      **** apply lattice constraints ****
       if (control_optimize_cell_strategy().eq.1) then

*         *** dE/da ==0,dE/db==0 .... dE/dgamma==0 ***
          do i=1,6
            if (control_optimize_lattice(i).eq.0) then
              lstress(i) = 0.0d0
              if (taskid.eq.MASTER) write(*,1800) i
            end if
          end do


*      **** apply lattice_vector constraints ****
       else if (control_optimize_cell_strategy().eq.2) then
          do v=1,3
          do u=1,3
            if (control_optimize_lattice_vectors(u,v).eq.0)  then
              stress(u,v) = 0.0d0
              if (taskid.eq.MASTER) write(*,1801) u,v
            end if
          end do
          end do
       end if

      return
 1800 FORMAT('Lattice Constraint: gradient(',i1,')-->0')
 1801 FORMAT('Lattice Constraint: gradient(',i1,',',i1,')-->0')
      end

*     *******************************
*     *                             *
*     *     cell_lattice_gradient   *
*     *                             *
*     *******************************
*
*   Given the gradient, stress, wrt to lattice vectors
*  this routine calculates the gradient,lstress, wrt a,b,c,alpha,beta,gamma.
*
      subroutine cell_lattice_gradient(stress,lstress)
      implicit none
      real*8 stress(3,3)
      real*8 lstress(6)

      integer i,u,v
      real*8  lvector(6)

      call dcopy(6,0.0d0,0,lstress,1)
      do v=1,3
      do u=1,3
        call cell_lattice_lvector(u,v,lvector)
        do i=1,6
          lstress(i) = lstress(i) + stress(u,v)*lvector(i)
        end do
      end do
      end do

      return
      end


c*     *******************************
c*     *                             *
c*     *     cell_lattice_avector    *
c*     *                             *
c*     *******************************
c*   This routine calculates the jacobian,
c*         dl/da11, dl/da12, dl/da13
c*         dl/da21, dl/da22, dl/da23
c*         dl/da31, dl/da32, dl/da33
c*   where l = a,b,c,alpha,beta,gamma
c*
c      subroutine cell_lattice_avector(i,avector)
c      implicit none
c      integer i
c      real*8 avector(3,3)
c      return
c      end


*     *******************************
*     *                             *
*     *     cell_lattice_lvector    *
*     *                             *
*     *******************************
*   This routine calculates the jacobian,
*      dA_ij/da, dA_ij/db, dA_ij/dc, dA_ij/dalpha, dA_ij/dbeta, dA_ij/dgamma
*
      subroutine cell_lattice_lvector(i,j,lvector)
      implicit none
      integer i,j
      real*8 lvector(6)

*     *** local variables ***
      real*8   a,b,c,alpha,beta,gamma,pi

      real*8  dad11,dab2d11,dac2d11,dbotd11,dtopd11,dbetad11,dgammad11
      real*8  dad21,dab2d21,dac2d21,dbotd21,dtopd21,dbetad21,dgammad21
      real*8  dad31,dab2d31,dac2d31,dbotd31,dtopd31,dbetad31,dgammad31
      real*8  dbd12,dbc2d12,dab2d12,dbotd12,dtopd12,dalphad12,dgammad12
      real*8  dbd22,dbc2d22,dab2d22,dbotd22,dtopd22,dalphad22,dgammad22
      real*8  dbd32,dbc2d32,dab2d32,dbotd32,dtopd32,dalphad32,dgammad32
      real*8  dcd33,dbc2d33,dac2d33,dbotd33,dtopd33,dalphad33,dbetad33
      real*8  dcd13,dbc2d13,dac2d13,dbotd13,dtopd13,dalphad13,dbetad13
      real*8  dcd23,dbc2d23,dac2d23,dbotd23,dtopd23,dalphad23,dbetad23

      real*8  ac2,ab2,top,bc2,bot
      real*8  amat(3,3),vol,c1,c2,c3,s3

      integer ii, jj

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita

      call lattice_abc_abg(a,b,c,alpha,beta,gamma)
      pi = 4.0d0*datan(1.0d0)
      alpha = alpha*pi/180.0d0
      beta  = beta *pi/180.0d0
      gamma = gamma*pi/180.0d0

      c1=dcos(alpha)
      c2=dcos(beta)
      c3=dcos(gamma)
      s3=dsin(gamma)
      vol = (a*b*c)*sqrt(1-(c1**2+c2**2+c3**2)+(2*c1*c2*c3))

      do jj = 1, 3
       do ii = 1, 3
         amat(ii,jj) = lattice_unita(ii,jj)
       end do
      enddo

         ! cos(angle) = f, where f is complicated
         ! angle = arcos(f)
         ! d/dx angle = d/dx arccos(f), this what we want
         ! But, we also know that:
         ! d/dx arccos(f) = -1/sqrt(1-f**2) * df/dx
         ! Therefore:
         ! d/dx angle = -1/sqrt(1-f**2) * df/dx
         !            = -1/sqrt(1-cos(angle)**2) * df/dx
         !            = -1/sqrt(sin(angle)**2) * df/dx
         !            = (df/dx) * (-1/sin(angle))
         ! d(angle)/dx = -Cosecant(angle) * (df/dx)
         ! 
         ! f is always in the form f=(top/bot), thus
         ! d(angle)/dx =-1/sin(angle)*(1/bot**2)*(bot*(dtop/dx) - top*(dbot/dx))

      if (     (i.eq.1).and.(j.eq.1)) then
!        Changing 1,1 (part of "a") changes gamma, beta, and a. Not b,c,alpha
         lvector(2) = 0.0d0
         lvector(3) = 0.0d0
         lvector(4) = 0.0d0
         ! a = sqrt(amat(1,1)**2 + amat(2,1)**2 + amat(3,1)**2)
         dad11 = amat(1,1)/a
         lvector(1) = dad11 
         ac2 = (amat(1,3)-amat(1,1))**2 + (amat(2,3)-amat(2,1))**2 +
     >         (amat(3,3)-amat(3,1))**2 ! square of length of a to c vector
         dac2d11 = (amat(1,3)-amat(1,1))*(2.0d0)*(-1.0d0)
         ab2 = (amat(1,1)-amat(1,2))**2 + (amat(2,1)-amat(2,2))**2 +
     >         (amat(3,1)-amat(3,2))**2
         dab2d11 = (amat(1,1)-amat(1,2))*(2.0d0)*(1.0d0)

         ! beta = acos((c*c + a*a - ac2)/(2.0d0*c*a))
         bot=(2.0d0*c*a)
         dbotd11 = 2.0d0*c*dad11
         top = (c*c + a*a - ac2)
         dtopd11 = (0+2*a*dad11 - dac2d11)
         dbetad11 = (-1.0d0/sin(beta))*(1.0d0/bot**2)*
     >              (bot*dtopd11 - top*dbotd11)
         lvector(5) = dbetad11

         ! gamma = acos((a*a + b*b - ab2)/(2.0d0*a*b))
         bot = (2.0d0*a*b)
         dbotd11 = 2.0d0*b*dad11
         top = (a*a + b*b - ab2)
         dtopd11 = 2*a*dad11 + 0 - dab2d11
         dgammad11 = (-1.0d0/sin(gamma))*(1.0d0/bot**2)*
     >              (bot*dtopd11 - top*dbotd11)
         lvector(6) = dgammad11

      else if ((i.eq.2).and.(j.eq.1)) then
!        Changing 2,1 (part of "a") changes gamma, beta, and a. Not b,c,alpha
         lvector(2) = 0.0d0
         lvector(3) = 0.0d0
         lvector(4) = 0.0d0
         ! a = sqrt(amat(1,1)**2 + amat(2,1)**2 + amat(3,1)**2)
         dad21 = amat(2,1)/a
         lvector(1) = dad21 
         ac2 = (amat(1,3)-amat(1,1))**2 + (amat(2,3)-amat(2,1))**2 +
     >         (amat(3,3)-amat(3,1))**2 ! square of length of a to c vector
         dac2d21 = (amat(2,3)-amat(2,1))*(2.0d0)*(-1.0d0)
         ab2 = (amat(1,1)-amat(1,2))**2 + (amat(2,1)-amat(2,2))**2 +
     >         (amat(3,1)-amat(3,2))**2
         dab2d21 = (amat(2,1)-amat(2,2))*(2.0d0)*(1.0d0)

         ! beta = acos((c*c + a*a - ac2)/(2.0d0*c*a))
         bot=(2.0d0*c*a)
         dbotd21 = 2.0d0*c*dad21
         top = (c*c + a*a - ac2)
         dtopd21 = (0+2*a*dad21 - dac2d21)
         dbetad21 = (-1.0d0/sin(beta))*(1.0d0/bot**2)*
     >              (bot*dtopd21 - top*dbotd21)
         lvector(5) = dbetad21

         ! gamma = acos((a*a + b*b - ab2)/(2.0d0*a*b))
         bot = (2.0d0*a*b)
         dbotd21 = 2.0d0*b*dad21
         top = (a*a + b*b - ab2)
         dtopd21 = 2*a*dad21 + 0 - dab2d21
         dgammad21 = (-1.0d0/sin(gamma))*(1.0d0/bot**2)*
     >              (bot*dtopd21 - top*dbotd21)
         lvector(6) = dgammad21

      else if ((i.eq.3).and.(j.eq.1)) then
!        Changing 3,1 (part of "a") changes gamma, beta, and a. Not b,c,alpha
         lvector(2) = 0.0d0
         lvector(3) = 0.0d0
         lvector(4) = 0.0d0
         ! a = sqrt(amat(1,1)**2 + amat(2,1)**2 + amat(3,1)**2)
         dad31 = amat(3,1)/a
         lvector(1) = dad31
         ac2 = (amat(1,3)-amat(1,1))**2 + (amat(2,3)-amat(2,1))**2 +
     >         (amat(3,3)-amat(3,1))**2 ! square of length of a to c vector
         dac2d31 = (amat(3,3)-amat(3,1))*(2.0d0)*(-1.0d0)
         ab2 = (amat(1,1)-amat(1,2))**2 + (amat(2,1)-amat(2,2))**2 +
     >         (amat(3,1)-amat(3,2))**2
         dab2d31 = (amat(3,1)-amat(3,2))*(2.0d0)*(1.0d0)

         ! beta = acos((c*c + a*a - ac2)/(2.0d0*c*a))
         bot=(2.0d0*c*a)
         dbotd31 = 2.0d0*c*dad31
         top = (c*c + a*a - ac2)
         dtopd31 = (0+2*a*dad31 - dac2d31)
         dbetad31 = (-1.0d0/sin(beta))*(1.0d0/bot**2)*
     >              (bot*dtopd31 - top*dbotd31)
         lvector(5) = dbetad31

         ! gamma = acos((a*a + b*b - ab2)/(2.0d0*a*b))
         bot = (2.0d0*a*b)
         dbotd31 = 2.0d0*b*dad31
         top = (a*a + b*b - ab2)
         dtopd31 = 2*a*dad31 + 0 - dab2d31
         dgammad31 = (-1.0d0/sin(gamma))*(1.0d0/bot**2)*
     >              (bot*dtopd31 - top*dbotd31)
         lvector(6) = dgammad31

      else if ((i.eq.1).and.(j.eq.2)) then
!        Changing 1,2 (part of "b") changes gamma, alpha, and b. Not a,c,beta
         lvector(1) = 0.0d0
         lvector(3) = 0.0d0
         lvector(5) = 0.0d0
         !b = sqrt(amat(1,2)**2 + amat(2,2)**2 + amat(3,2)**2)
         dbd12 = amat(1,2)/b
         lvector(2) = dbd12

         bc2 = (amat(1,2)-amat(1,3))**2 + (amat(2,2)-amat(2,3))**2 +
     >         (amat(3,2)-amat(3,3))**2
         dbc2d12 = (amat(1,2)-amat(1,3))*(2.0d0)*(1.0d0)
         ab2 = (amat(1,1)-amat(1,2))**2 + (amat(2,1)-amat(2,2))**2 +
     >         (amat(3,1)-amat(3,2))**2
         dab2d12 = (amat(1,1)-amat(1,2))*(2.0d0)*(-1.0d0)

         ! alpha = acos((b*b + c*c - bc2)/(2.0d0*b*c))
         bot = (2.0d0*b*c)
         dbotd12 = 2.0d0*dbd12*c
         top = (b*b + c*c - bc2)
         dtopd12 = 2*b*dbd12 + 0 - dbc2d12
         dalphad12 = (-1.0d0/sin(alpha))*(1.0d0/bot**2)*
     >               (bot*dtopd12 - top*dbotd12)
         lvector(4) = dalphad12

         ! gamma = acos((a*a + b*b - ab2)/(2.0d0*a*b))
         bot = (2.0d0*a*b)
         dbotd12 = 2.0d0*a*dbd12
         top = (a*a + b*b - ab2)
         dtopd12 = 0 + 2*b*dbd12 - dab2d12
         dgammad12 = (-1.0d0/sin(gamma))*(1.0d0/bot**2)*
     >               (bot*dtopd12 - top*dbotd12)
         lvector(6) = dgammad12

      else if ((i.eq.2).and.(j.eq.2)) then
!        Changing 2,2 (part of "b") changes gamma, alpha, and b. Not a,c,beta
         lvector(1) = 0.0d0
         lvector(3) = 0.0d0
         lvector(5) = 0.0d0

         ! b = sqrt(amat(1,2)**2 + amat(2,2)**2 + amat(3,2)**2)
         dbd22 = amat(2,2)/b
         lvector(2) = dbd22

         bc2 = (amat(1,2)-amat(1,3))**2 + (amat(2,2)-amat(2,3))**2 +
     >         (amat(3,2)-amat(3,3))**2
         dbc2d22 = (amat(2,2)-amat(2,3))*(2.0d0)*(1.0d0)
         ab2 = (amat(1,1)-amat(1,2))**2 + (amat(2,1)-amat(2,2))**2 +
     >         (amat(3,1)-amat(3,2))**2
         dab2d22 = (amat(2,1)-amat(2,2))*(2.0d0)*(-1.0d0)

         ! alpha = acos((b*b + c*c - bc2)/(2.0d0*b*c))
         bot = (2.0d0*b*c)
         dbotd22 = 2.0d0*dbd22*c
         top = (b*b + c*c - bc2)
         dtopd22 = 2*b*dbd22 + 0 - dbc2d22
         dalphad22 = (-1.0d0/sin(alpha))*(1.0d0/bot**2)*
     >               (bot*dtopd22 - top*dbotd22)
         lvector(4) = dalphad22

         ! gamma = acos((a*a + b*b - ab2)/(2.0d0*a*b))
         bot = (2.0d0*a*b)
         dbotd22 = 2.0d0*a*dbd22
         top = (a*a + b*b - ab2)
         dtopd22 = 0 + 2*b*dbd22 - dab2d22
         dgammad22 = (-1.0d0/sin(gamma))*(1.0d0/bot**2)*
     >               (bot*dtopd22 - top*dbotd22)
         lvector(6) = dgammad22

      else if ((i.eq.3).and.(j.eq.2)) then
!        Changing 3,2 (part of "b") changes gamma, alpha, and b. Not a,c,beta
         lvector(1) = 0.0d0
         lvector(3) = 0.0d0
         lvector(5) = 0.0d0

         ! b = sqrt(amat(1,2)**2 + amat(2,2)**2 + amat(3,2)**2)
         dbd32 = amat(3,2)/b
         lvector(2) = dbd32

         bc2 = (amat(1,2)-amat(1,3))**2 + (amat(2,2)-amat(2,3))**2 +
     >         (amat(3,2)-amat(3,3))**2
         dbc2d32 = (amat(3,2)-amat(3,3))*(2.0d0)*(1.0d0)
         ab2 = (amat(1,1)-amat(1,2))**2 + (amat(2,1)-amat(2,2))**2 +
     >         (amat(3,1)-amat(3,2))**2
         dab2d32 = (amat(3,1)-amat(3,2))*(2.0d0)*(-1.0d0)

         ! alpha = acos((b*b + c*c - bc2)/(2.0d0*b*c))
         bot = (2.0d0*b*c)
         dbotd32 = 2.0d0*dbd32*c
         top = (b*b + c*c - bc2)
         dtopd32 = 2*b*dbd32 + 0 - dbc2d32
         dalphad32 = (-1.0d0/sin(alpha))*(1.0d0/bot**2)*
     >               (bot*dtopd32 - top*dbotd32)
         lvector(4) = dalphad32

         ! gamma = acos((a*a + b*b - ab2)/(2.0d0*a*b))
         bot = (2.0d0*a*b)
         dbotd32 = 2.0d0*a*dbd32
         top = (a*a + b*b - ab2)
         dtopd32 = 0 + 2*b*dbd32 - dab2d32
         dgammad32 = (-1.0d0/sin(gamma))*(1.0d0/bot**2)*
     >               (bot*dtopd32 - top*dbotd32)
         lvector(6) = dgammad32

      else if ((i.eq.1).and.(j.eq.3)) then
!        Changing 2,3 (part of "c") changes alpha,beta and c. Not a,b,gamma
         lvector(1) = 0.0d0
         lvector(2) = 0.0d0
         lvector(6) = 0.0d0
         !c = sqrt(amat(1,3)**2 + amat(2,3)**2 + amat(3,3)**2)
         dcd13 = amat(1,3)/c
         lvector(3) = dcd13

         bc2 = (amat(1,2)-amat(1,3))**2 + (amat(2,2)-amat(2,3))**2 +
     >         (amat(3,2)-amat(3,3))**2
         dbc2d13 = (amat(1,2)-amat(1,3))*(2.0d0)*(-1.0d0)
         ac2 = (amat(1,1)-amat(1,3))**2 + (amat(2,1)-amat(2,3))**2 +
     >         (amat(3,1)-amat(3,3))**2
         dac2d13 = (amat(1,1)-amat(1,3))*(2.0d0)*(-1.0d0)

         ! alpha = acos((b*b + c*c - bc2)/(2.0d0*b*c))
         bot = (2.0d0*b*c)
         dbotd13 = 2.0d0*b*dcd13
         top = (b*b + c*c - bc2)
         dtopd13 = 0 + 2*c*dcd13 - dbc2d13
         dalphad13 = (-1.0d0/sin(alpha))*(1.0d0/bot**2)*
     >               (bot*dtopd13 - top*dbotd13)
         lvector(4) = dalphad13

         ! beta = acos((c*c + a*a - ac2)/(2.0d0*c*a))
         bot=(2.0d0*c*a)
         dbotd13 = 2.0d0*dcd13*a
         top = (c*c + a*a - ac2)
         dtopd13 = (2*c*dcd13 + 0 - dac2d13)
         dbetad13 = (-1.0d0/sin(beta))*(1.0d0/bot**2)*
     >              (bot*dtopd13 - top*dbotd13)
         lvector(5) = dbetad13

      else if ((i.eq.2).and.(j.eq.3)) then
!        Changing 2,3 (part of "c") changes alpha,beta and c. Not a,b,gamma
         lvector(1) = 0.0d0
         lvector(2) = 0.0d0
         lvector(6) = 0.0d0
         !c = sqrt(amat(1,3)**2 + amat(2,3)**2 + amat(3,3)**2)
         dcd23 = amat(2,3)/c
         lvector(3) = dcd23

         bc2 = (amat(1,2)-amat(1,3))**2 + (amat(2,2)-amat(2,3))**2 +
     >         (amat(3,2)-amat(3,3))**2
         dbc2d23 = (amat(2,2)-amat(2,3))*(2.0d0)*(-1.0d0)
         ac2 = (amat(1,1)-amat(1,3))**2 + (amat(2,1)-amat(2,3))**2 +
     >         (amat(3,1)-amat(3,3))**2
         dac2d23 = (amat(2,1)-amat(2,3))*(2.0d0)*(-1.0d0)

         ! alpha = acos((b*b + c*c - bc2)/(2.0d0*b*c))
         bot = (2.0d0*b*c)
         dbotd23 = 2.0d0*b*dcd23
         top = (b*b + c*c - bc2)
         dtopd23 = 0 + 2*c*dcd23 - dbc2d23
         dalphad23 = (-1.0d0/sin(alpha))*(1.0d0/bot**2)*
     >               (bot*dtopd23 - top*dbotd23)
         lvector(4) = dalphad23

         ! beta = acos((c*c + a*a - ac2)/(2.0d0*c*a))
         bot=(2.0d0*c*a)
         dbotd23 = 2.0d0*dcd23*a
         top = (c*c + a*a - ac2)
         dtopd23 = (2*c*dcd23 + 0 - dac2d23)
         dbetad23 = (-1.0d0/sin(beta))*(1.0d0/bot**2)*
     >              (bot*dtopd23 - top*dbotd23)
         lvector(5) = dbetad23

      else if ((i.eq.3).and.(j.eq.3)) then
!        Changing 3,3 (part of "c") changes alpha,beta and c. Not a,b,gamma
         lvector(1) = 0.0d0
         lvector(2) = 0.0d0
         lvector(6) = 0.0d0
         !c = sqrt(amat(1,3)**2 + amat(2,3)**2 + amat(3,3)**2)
         dcd33 = amat(3,3)/c
         lvector(3) = dcd33

         bc2 = (amat(1,2)-amat(1,3))**2 + (amat(2,2)-amat(2,3))**2 +
     >         (amat(3,2)-amat(3,3))**2
         dbc2d33 = (amat(3,2)-amat(3,3))*(2.0d0)*(-1.0d0)
         ac2 = (amat(1,1)-amat(1,3))**2 + (amat(2,1)-amat(2,3))**2 +
     >         (amat(3,1)-amat(3,3))**2
         dac2d33 = (amat(3,1)-amat(3,3))*(2.0d0)*(-1.0d0)

         ! alpha = acos((b*b + c*c - bc2)/(2.0d0*b*c))
         bot = (2.0d0*b*c)
         dbotd33 = 2.0d0*b*dcd33
         top = (b*b + c*c - bc2)
         dtopd33 = 0 + 2*c*dcd33 - dbc2d33
         dalphad33 = (-1.0d0/sin(alpha))*(1.0d0/bot**2)*
     >               (bot*dtopd33 - top*dbotd33)
         lvector(4) = dalphad33

         ! beta = acos((c*c + a*a - ac2)/(2.0d0*c*a))
         bot=(2.0d0*c*a)
         dbotd33 = 2.0d0*dcd33*a
         top = (c*c + a*a - ac2)
         dtopd33 = (2*c*dcd33 + 0 - dac2d33)
         dbetad33 = (-1.0d0/sin(beta))*(1.0d0/bot**2)*
     >              (bot*dtopd33 - top*dbotd33)
         lvector(5) = dbetad33
      endif

      return
      end

     
