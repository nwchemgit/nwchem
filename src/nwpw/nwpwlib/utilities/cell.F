*
*     $Id: cell.F,v 1.2 2008-12-15 21:03:03 bylaska Exp $
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
      real*8  V,Va,Vb,Vc,Vd,Ve,Vf
      real*8  V2,V2a,V2b,V2c,V2d,V2e,V2f
      real*8  bot,bota,botb,botc,botd,bote,botf
      real*8  g11,g21,g31,g12,g22,g32,g13,g23,g33
      real*8  g11a,g21a,g31a,g12a,g22a,g32a,g13a,g23a,g33a
      real*8  g11b,g21b,g31b,g12b,g22b,g32b,g13b,g23b,g33b
      real*8  g11c,g21c,g31c,g12c,g22c,g32c,g13c,g23c,g33c
      real*8  g11d,g21d,g31d,g12d,g22d,g32d,g13d,g23d,g33d
      real*8  g11e,g21e,g31e,g12e,g22e,g32e,g13e,g23e,g33e
      real*8  g11f,g21f,g31f,g12f,g22f,g32f,g13f,g23f,g33f

      call lattice_abc_abg(a,b,c,alpha,beta,gamma)
      pi = 4.0d0*datan(1.0d0)
      alpha = alpha*pi/180.0d0
      beta  = beta *pi/180.0d0
      gamma = gamma*pi/180.0d0

      lvector(1) = 0.0d0
      lvector(2) = 0.0d0
      lvector(3) = 0.0d0
      lvector(4) = 0.0d0
      lvector(5) = 0.0d0
      lvector(6) = 0.0d0

      if ((i.eq.1).and.(j.eq.1)) then
         lvector(1) = dsin(gamma)
         lvector(6) = a*dcos(gamma)

      else if ((i.eq.2).and.(j.eq.1)) then
         lvector(1) = dcos(gamma)
         lvector(6) = -a*dsin(gamma)

      else if ((i.eq.2).and.(j.eq.2)) then
         lvector(2) = 1.0d0

      else if ((i.eq.1).and.(j.eq.3)) then
         lvector(3) = (dcos(beta)-dcos(alpha)*dcos(gamma))/dsin(gamma)
         lvector(4) = c*dsin(alpha)*dcos(gamma)/dsin(gamma)
         lvector(5) = -c*dsin(beta)/dsin(gamma)
         lvector(6) = c*cos(alpha) 
     >              + c*cos(alpha)*(cos(gamma)/sin(gamma))**2

      else if ((i.eq.2).and.(j.eq.3)) then
         lvector(3) = dcos(gamma)
         lvector(6) = -c*dsin(gamma)

      else if ((i.eq.3).and.(j.eq.3)) then
         g11 = a*a
         g11a = 2.0d0*a
         g11b = 0.0d0
         g11c = 0.0d0
         g11d = 0.0d0
         g11e = 0.0d0
         g11f = 0.0d0

         g21 = a*b*dcos(gamma)
         g21a = b*dcos(gamma)
         g21b = a*dcos(gamma)
         g21c = 0.0d0
         g21d = 0.0d0
         g21e = 0.0d0
         g21f = -a*b*dsin(gamma)
     
         g31 = a*c*dcos(beta)
         g31a = c*dcos(beta)
         g31b = 0.0d0
         g31c = a*dcos(beta)
         g31d = 0.0d0
         g31e = -a*c*dsin(beta)
         g31f = 0.0d0

         g12 = a*b*dcos(gamma)
         g12a = b*dcos(gamma)
         g12b = a*dcos(gamma)
         g12c = 0.0d0
         g12d = 0.0d0
         g12e = 0.0d0
         g12f = -a*b*dsin(gamma)

         g22 = b*b
         g22a = 0.0d0
         g22b = 2.0d0*b
         g22c = 0.0d0
         g22d = 0.0d0
         g22e = 0.0d0
         g22f = 0.0d0

         g32 = b*c*dcos(alpha)
         g32a = 0.0d0
         g32b = c*dcos(alpha)
         g32c = b*dcos(alpha)
         g32d = -b*c*dsin(alpha)
         g32e = 0.0d0
         g32f = 0.0d0

         g13 = a*b*dcos(beta)
         g13a = c*dcos(beta)
         g13b = 0.0d0
         g13c = a*dcos(beta)
         g13d = 0.0d0
         g13e = -a*c*dsin(beta)
         g13f = 0.0d0

         g23 = b*c*dcos(alpha)
         g23a = 0.0d0
         g23b = c*dcos(alpha)
         g23c = b*dcos(alpha)
         g23d = -b*c*dsin(alpha)
         g23e  = 0.0d0
         g23f = 0.0d0

         g33 = c*c
         g33a = 0.0d0
         g33b = 0.0d0
         g33c = 2.0d0*c
         g33d = 0.0d0
         g33e = 0.0d0
         g33f = 0.0d0

*                             |g11  g12  g13|
*        V = sqrt(V2),   V2 = |g21  g22  g23|
*                             |g31  g32  g33|
*
         V2 = g11*(g22*g33 - g23*g32) 
     >      + g12*(g23*g31 - g21*g33) 
     >      + g13*(g21*g32 - g22*g31)

         V2a = g11a*(g22*g33 - g23*g32) 
     >       + g11*(g22a*g33 + g22*g33a - g23a*g32 - g23*g32a) 
     >       + g12a*(g23*g31 - g21*g33) 
     >       + g12*(g23a*g31 + g23*g31a - g21a*g33 - g21*g33a) 
     >       + g13a*(g21*g32 - g22*g31)
     >       + g13*(g21a*g32 + g21*g32a - g22a*g31 - g22*g31a)

         V2b = g11b*(g22*g33 - g23*g32) 
     >       + g11*(g22b*g33 + g22*g33b - g23b*g32 - g23*g32b) 
     >       + g12b*(g23*g31 - g21*g33) 
     >       + g12*(g23b*g31 + g23*g31b - g21b*g33 - g21*g33b) 
     >       + g13b*(g21*g32 - g22*g31)
     >       + g13*(g21b*g32 + g21*g32b - g22b*g31 - g22*g31b)

         V2c = g11c*(g22*g33 - g23*g32) 
     >       + g11*(g22c*g33 + g22*g33c - g23c*g32 - g23*g32c) 
     >       + g12c*(g23*g31 - g21*g33) 
     >       + g12*(g23c*g31 + g23*g31c - g21c*g33 - g21*g33c) 
     >       + g13c*(g21*g32 - g22*g31)
     >       + g13*(g21c*g32 + g21*g32c - g22c*g31 - g22*g31c)

         V2d = g11d*(g22*g33 - g23*g32) 
     >       + g11*(g22d*g33 + g22*g33d - g23d*g32 - g23*g32d) 
     >       + g12d*(g23*g31 - g21*g33) 
     >       + g12*(g23d*g31 + g23*g31d - g21d*g33 - g21*g33d) 
     >       + g13d*(g21*g32 - g22*g31)
     >       + g13*(g21d*g32 + g21*g32d - g22d*g31 - g22*g31d)

         V2e = g11e*(g22*g33 - g23*g32) 
     >       + g11*(g22e*g33 + g22*g33e - g23e*g32 - g23*g32e) 
     >       + g12e*(g23*g31 - g21*g33) 
     >       + g12*(g23e*g31 + g23*g31e - g21e*g33 - g21*g33e) 
     >       + g13e*(g21*g32 - g22*g31)
     >       + g13*(g21e*g32 + g21*g32e - g22e*g31 - g22*g31e)

         V2f = g11f*(g22*g33 - g23*g32) 
     >       + g11*(g22f*g33 + g22*g33f - g23f*g32 - g23*g32f) 
     >       + g12f*(g23*g31 - g21*g33) 
     >       + g12*(g23f*g31 + g23*g31f - g21f*g33 - g21*g33f) 
     >       + g13f*(g21*g32 - g22*g31)
     >       + g13*(g21f*g32 + g21*g32f - g22f*g31 - g22*g31f)

         !a33 = V/(a*b*sin(gamma))
         !a33 = V/bot
         !a33a = Va/bot - V*bota/(bot**2)


         V = dsqrt(V2)
         Va   = (0.5d0*V2a/V)
         Vb   = (0.5d0*V2b/V)
         Vc   = (0.5d0*V2c/V)
         Vd   = (0.5d0*V2d/V)
         Ve   = (0.5d0*V2e/V)
         Vf   = (0.5d0*V2f/V)

         bot = a*b*dsin(gamma)
         bota = b*dsin(gamma)
         botb = a*dsin(gamma)
         botc = 0.0d0
         botd = 0.0d0
         bote = 0.0d0
         botf = a*b*dcos(gamma)


         lvector(1) = Va/bot - V*bota/(bot**2)
         lvector(2) = Vb/bot - V*botb/(bot**2)
         lvector(3) = Vc/bot - V*botc/(bot**2)
         lvector(4) = Vd/bot - V*botd/(bot**2)
         lvector(5) = Ve/bot - V*bote/(bot**2)
         lvector(6) = Vf/bot - V*botf/(bot**2)
      endif

      return
      end

     
