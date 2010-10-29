      subroutine sym_op_type(r, op_type, order, power, direction)
C$Id$
      implicit none
#include "errquit.fh"
      double precision r(3,3)
      character*(*) op_type
      integer order, power
      double precision direction(3)
c
c                     proper         improper
c     operation       rotn(theta)    rotn(theta)   mirror  inversion  identity
c     
c     trace           1+2cos(theta) -1+2cos(theta)    1      -3          3
c
c     determinant     1             -1               -1      -1          1
c
      double precision trace, deter, angle, twopi, pi, acos_arg
      double precision deter3, e, x(3), tmp(3), z(3), phase, dnorm
      double precision xnorm, sintheta
      integer m, n, i, j, k, repeat
      external deter3
c
      trace = r(1,1)+r(2,2)+r(3,3)
      deter = deter3(r)
      pi    = 4.0d0*atan(1.0d0)
      twopi = 2.0d0 * pi
      e     = exp(1.0d0)
      order = 0
      op_type = 'unknown'
c      
c     Determine type of operation
c
      if (abs(deter - 1.0d0) .lt. 1.0e-6) then
         if (abs(trace - 3.0d0) .lt. 1.0e-6) then
            op_type = 'identity'
         else
            op_type = 'rotation'
         endif
      else if (abs(deter + 1.0d0) .lt. 1.0e-6) then
         if (abs(trace+3.0d0) .lt. 1.0e-6) then
            op_type = 'inversion'
         else if (abs(trace - 1.0d0) .lt. 1.0e-6) then
            op_type = 'mirror'
         else
            op_type = 'improper rotation'
         endif
      else
         write(6,*) ' sym_op_type: inconsistency ', deter, trace
         call errquit('sym_op_type: deter/trace?', 0, GEOM_ERR)
      endif
c
c     Determine angles of rotations
c
cedo      write(0,*) ' ttt ',trace
      if (op_type.eq.'rotation' .or. op_type.eq.'improper rotation')then
         if (op_type.eq.'rotation') then
            acos_arg = 0.5d0*(trace - 1.0d0)
         else
            acos_arg = 0.5d0*(trace + 1.0d0)
         endif
         if (abs(acos_arg) .gt. (1.0d0+1d-6)) then
            write(6,*) ' sym_op_type: invalid trace for rotation', trace
            call errquit('sym_op_type: bad trace', 0, GEOM_ERR)
         else if (abs(acos_arg) .gt. 1.0d0) then
            acos_arg = sign(1.0d0,acos_arg)
         endif
         angle = acos(acos_arg) ! Angle should be 2*pi*m/n
cedo         write(0,*) ' aaangle ',angle,acos_arg
         do n = 2, 12
            do m = 1, n-1
               if (abs(twopi*dble(m)/dble(n) - angle) .lt.
     &            1.0d-6) goto 10
            enddo
         enddo
         write(6,*) ' sym_op_type: angle is not 2*pi*m/n ', angle
         call errquit('sym_op_type: bad', 0, GEOM_ERR)
 10      order = n
         power = m
cedo         write(0,*) ' mmmm ',m
c     
      endif
c
c     Determine directions by projecting out fixed point
c
      x(1) = pi                 ! This point is guaranteed generic
      x(2) = e
      x(3) = 1
      if (op_type .eq. 'mirror') then
         repeat = 1
         phase  = -1
      else if (op_type .eq. 'rotation') then
         repeat = order - 1
         phase  = 1
      else if (op_type .eq. 'improper rotation') then
         phase = -1
         if (mod(order,2) .eq. 0) then
            repeat = order - 1
         else
            repeat = 2*order - 1
         endif
      else
         repeat = 0
         phase  = 0
      endif
c
cedo      write(0,*) ' repeat ',repeat
      if (repeat .gt. 0) then
        direction(1) = x(1)
        direction(2) = x(2)
        direction(3) = x(3)
        do i = 1, repeat
           do j = 1, 3
              tmp(j) = 0.0d0
              do k = 1, 3
                 tmp(j) = tmp(j) + r(j,k)*x(k)
              enddo
           enddo
           do j = 1, 3
              x(j) = phase*tmp(j)
              direction(j) = direction(j) + x(j)
           enddo
        enddo
        do i = 1, 3
           if (abs(direction(i)).lt.1d-10) direction(i) = 0.0d0
        enddo
        dnorm = direction(1)**2+direction(2)**2+direction(3)**2
        dnorm = sqrt(dnorm)
        if (dnorm .lt. 1.0d-3) call errquit
     $       (' sym_op_type: zero norm for direction?',0, GEOM_ERR)
c
c     Establish convention for positive directions
c     a) make z>0
c     b) if (z==0) make x>0
c     c) if (z==x==0) make y>0
c
        if (direction(3).lt.0.0d0) then
           dnorm = -dnorm
        else if (direction(3).eq.0.0d0) then
           if (direction(1).lt.0.0d0) then
              dnorm = -dnorm
           else if (direction(1).eq.0.0d0) then
              if (direction(2).lt.0.0d0) then
                 dnorm = -dnorm
              endif
           endif
        endif
        direction(1) = direction(1)/dnorm
        direction(2) = direction(2)/dnorm
        direction(3) = direction(3)/dnorm
c
c     verify direction is preserved
c
        do j = 1, 3
           tmp(j) = 0.0d0
           do k = 1, 3
              tmp(j) = tmp(j) + r(j,k)*direction(k)
           enddo
           if (abs(tmp(j)*phase - direction(j)) .gt. 1.0e-6) 
     $          call errquit
     $          (' sym_op_type: direction is not correct? ',0, GEOM_ERR)
        enddo
      endif
c
c     Figure out if rotations are clockwise/counter clockwise
c
c     Apply rotation to unit vector ... cross product gives
c     sin(theta)*direction
c
      if ((op_type.eq.'rotation' .or. op_type.eq.'improper rotation')
     $     .and. order.gt.2)then
         x(1) = pi              ! This point is guaranteed generic
         x(2) = e
         x(3) = 1
         xnorm = sqrt(x(1)**2+x(2)**2+x(3)**2)
         x(1) = x(1)/xnorm
         x(2) = x(2)/xnorm
         x(3) = x(3)/xnorm
         do j = 1, 3
            tmp(j) = 0.0d0
            do k = 1, 3
               tmp(j) = tmp(j) + r(j,k)*x(k)
            enddo
         enddo
         call cross_product(x, tmp, z)
         sintheta=z(1)*direction(1)+z(2)*direction(2)+z(3)*direction(3)
         sintheta=sintheta*phase
         if (sintheta .lt. -1d-6) power = -power
ceod         write(0,*) ' sinth ',sintheta,power
      endif
c
c     Evaluation of power is mod(order) ... thus can end
c     up with seemingly even improper rotations
c
      if (op_type.eq.'improper rotation') then 
         if (mod(power,2).eq.0) then
            power = power + order
         elseif (power .gt. order) then
            power = 2*order - power
         endif
      endif
c
cedo      write(6,*) ' trace, deter, op_type, order, power ',
cedo     $     trace, deter, op_type, order, power 
c
      end
