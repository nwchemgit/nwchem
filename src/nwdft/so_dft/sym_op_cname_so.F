      subroutine sym_op_class_name_so(op_type, order, power, direction,
     $     group_name, class, rclass)
C$Id$
      implicit none
#include "inp.fh"
c     
      character*(*) op_type     ! [input]
      integer order, power      ! [input]
      double precision direction(3) ! [input]
      character*(*) group_name  ! [input]
      character*(*) class       ! [output]
      character*(*) rclass
      logical one_sided, two_sided  
c     
      integer abspow, rpow
c
c     Cn,   Sn,    Dn,    T,    O,    I 
c     Cnh,         Dnh,   Th,   Oh,   Ih 
c     Cnv,         Dnd,   Td,    
c     
c     Try to assign the class name of an operator given all the other info
c     
c     Powers of rotations in Cn, Cnh, Sn, T and Th  are all in separate 
c     classes (rotations about single-sided axis) 
c     For all the point groups, the normal point group 
c     elements are classified by the order and power of the rotations. 
c     For the above list of groups, the double elements belong to different 
c     seperate classes, which can also be designated by the order and power 
c     of the rotations. For the rest of the groups, the double elements 
c     belong to the same classes as the normal point groups. The double 
c     element CR, where C is a rotation of power and order, belongs to the 
c     same class as C' with (order-power) and order. Special cases, i.e., 
c     mirror planes and rotations by pi are handled seperately. 
c     
      one_sided = .false. 
      two_sided = .true. 
      if (power .lt. 0) then
         if (op_type.eq.'rotation' .or. mod(order,2).eq.0) then
            abspow = order+power
         else
            abspow = 2*order+power ! Odd improper rotation
         endif
      endif
      if (op_type.eq.'rotation' .or. mod(order,2).eq.0) then
         rpow = order-abspow
      else
         rpow = 2*order-abspow   ! Odd improper rotation
      endif
      if (group_name(1:1).eq.'S' .or. 
     $     group_name.eq.'T' .or. group_name.eq.'Th' .or.
     $     (group_name(1:1).eq.'C' .and.
     $     (group_name(3:3).eq.'h' .or. inp_strlen(group_name).eq.2)))
     $     then                 ! single-sided axis
         one_sided = .true.
         two_sided = .false. 
         if (op_type.eq.'rotation' .or. mod(order,2).eq.0) then
            rpow = order+abspow
         else
            rpow = 2*order+abspow ! Odd improper rotation
         endif
      endif
c
      class = ' '
      if (op_type .eq. 'identity') then
         class = 'e'
         rclass = 'e.r'
      else if (op_type .eq. 'inversion') then
         if (group_name.eq.'Ci' .or. group_name.eq.'S2') then
            class = 'i=s2'
            rclass = 'i=s2.r'
         else
            class = 'i'
            rclass = 'Ri'
         endif
      else if (op_type .eq. 'rotation') then 
         if (order .gt. 2) then 
            if(two_sided)then  
               if (power .gt. 1)then 
                  write(class,'(''c'',i1,''^'',i1)')order,power
                  write(rclass,'(''Rc'',i1,''^'',i1)')order,power
               else if(power .eq. 1)then 
                  write(class,'(''c'',i1)') order
                  write(rclass,'(''Rc'',i1)') order
               else if (power .eq. -1)then
                  write(class,'(''Rc'',i1)') order
                  write(rclass,'(''c'',i1)') order
               else if (power .lt. -1)then
                  write(class,'(''Rc'',i1,''^'',i1)')order,abs(power)
                  write(rclass,'(''c'',i1,''^'',i1)')order,abs(power)
               endif
            else if (one_sided)then 
               if (power .gt. 1)then 
                  write(class,'(''c'',i1,''^'',i1)')order,power
                  write(rclass,'(''Rc'',i1,''^'',i1)')order,power
               else if(power .eq. 1)then 
                  write(class,'(''c'',i1)') order
                  write(rclass,'(''Rc'',i1)') order
               else if(power .lt. 0)then 
                  write(class,'(''c'',i1,''^'',i1)')order,order+power
                  write(rclass,'(''Rc'',i1,''^'',i1)')order,order+power
               endif 
            endif
         else                   ! second order rotations
            write(class,'(''c2'',''^'',i1)')abs(power)
            if (group_name.eq.'D2' .or. group_name.eq.'D2h') then
               if (abs(direction(3)-1.0d0) .lt. 1.0e-6) then
                  class = 'c2z'
                  rclass = class
               else if (abs(direction(2)-1.0d0) .lt. 1.0e-6) then
                  class = 'c2y'
                  rclass = class
               else if (abs(direction(1)-1.0d0) .lt. 1.0e-6) then
                  class = 'c2x'
                  rclass = class
               endif
            else if (group_name.eq.'D4' .or. group_name.eq.'D4h') then
               if (abs(direction(3)-1.0d0) .lt. 1d-6) then
                  class = 'c2'
                  rclass = class
               else if ((abs(direction(1)-1.0d0) .lt. 1d-6) .or.
     $                 (abs(direction(2)-1.0d0) .lt. 1d-6)) then
                  class = 'c2'''
                  rclass = class 
               else
                  class = 'c2"'
                  rclass = class 
               endif
            else if (group_name.eq.'D2d' .or. group_name.eq.'D4d' .or. 
     $              group_name.eq.'D6d') then
               if (abs(direction(3)-1.0d0) .gt. 1.0e-6) then
                  class = 'c2'''
                  rclass = class
               endif
            else if (group_name(1:2).eq.'D6') then
c
c     C2' in plane of sigma.v ... -75, -15, 45 from X
c     
               if (abs(direction(3)-1.0d0).lt.1e-3) then
                  class = 'c2'
               else if (
     $                 abs(direction(2)-0.7071).lt.1d-3 .or.
     $                 abs(direction(2)+0.2588).lt.1d-3 .or.
     $                 abs(direction(2)+0.9659).lt.1d-3) then
                  class = 'c2'''
               else
                  class = 'c2"'
               endif
               rclass = class
            else if (group_name(1:1).eq.'O') then
c     
c     C2 normal to x/y/z axes, C2' imbetween
c     
               if (abs(direction(1)-1.0d0).le.1d-3 .or.
     $             abs(direction(2)-1.0d0).le.1d-3 .or.
     $             abs(direction(3)-1.0d0).le.1d-3) then
                  class = 'c2'
               else
                  class = 'c2'''
               endif
               rclass = class
            endif
         endif
      else if (op_type .eq. 'improper rotation') then 
         if(one_sided)then 
            if(power .gt. 1)then
               write(class,'("s",i1,"^",i1)')order,power
               write(rclass,'("Rs",i1,"^",i1)')order,power
            else if(power .eq. 1)then 
               write(class,'("s",i1)')order
               write(rclass,'("Rs",i1)')order
            else if(power .lt. 0)then 
               if(mod(order,2) .eq. 0)then 
                  write(class,'("s",i1,"^",i1)')order,order+power
                  write(rclass,'("Rs",i1,"^",i1)')order,order+power
               else 
                  write(class,'("s",i1,"^",i1)')order,2*order+power
                  write(rclass,'("Rs",i1,"^",i1)')order,2*order+power
               endif
            endif
         endif
         if(two_sided)then 
            if(power .gt. 1)then 
               write(class,'("s",i1,"^",i1)')order,power
               write(rclass,'("Rs",i1,"^",i1)')order,power
            else if(power .eq. 1)then 
               write(class,'("s",i1)')order
               write(rclass,'("Rs",i1)')order
            else if(power .eq. -1)then 
               write(class,'("Rs",i1)')order
               write(rclass,'("s",i1)')order
            else if(power .lt. -1)then 
               write(class,'("Rs",i1,"^",i1)')order,abs(power)
               write(rclass,'("s",i1,"^",i1)')order,abs(power)
            endif
         endif
      else if (op_type .eq. 'mirror') then
         if (group_name.eq.'C2v') then
            if (abs(direction(1)-1.0d0) .lt. 1.0d-6) then
               class = 'si.yz'
               rclass = class
            else  if (abs(direction(2)-1.0d0) .lt. 1.0d-6) then
               class = 'si.xz'
               rclass = class
            else
               write(6,*)' sym_op_class_name: c2v bad mirror'
            endif
         else if (group_name.eq.'D2h') then
            if (abs(direction(1)-1.0d0) .lt. 1.0d-6) then
               class = 'si.yz'
            else  if (abs(direction(2)-1.0d0) .lt. 1.0d-6) then
               class = 'si.xz'
            else  if (abs(direction(3)-1.0d0) .lt. 1.0d-6) then
               class = 'si.xy'
            else
               write(6,*)' sym_op_class_name: d2h bad mirror'
            endif
            rclass = class 
         else if (group_name(1:2) .eq. 'Ih') then
            class = 'si'
         else if (group_name(1:1) .eq. 'T' .or.
     $           group_name(3:3) .eq. 'd') then
            class = 'si.d'
         else if (group_name(1:2) .eq. 'Oh') then
c
c     sigmah normal to axes, sigmad between
c
            if (abs(direction(1)-1.0d0).le.1d-3 .or.
     $           abs(direction(2)-1.0d0).le.1d-3 .or.
     $           abs(direction(3)-1.0d0).le.1d-3) then
               class = 'si.h'
            else
               class = 'si.d'
            endif
         else if (abs(direction(3)-1.0d0) .lt. 1.0d-6) then
            class = 'si.h'
         else 
            class = 'si.v'
            if (group_name.eq.'C4v' .or. group_name.eq.'D4h') then
               if (abs(direction(1)-1.0d0).gt.1.0d-6 .and.
     $             abs(direction(2)-1.0d0).gt.1.0d-6) then
                  class = 'si.d'
               endif
            else if (group_name.eq.'C6v' .or. group_name.eq.'D6h') then
c
c     Adopt a convention that Sigma.v normals make angles of -45, 
c     +15 and 75 with the x-axis with the adopted direction convention
c
               if (
     $              abs(direction(2)+0.7071).lt.1d-3 .or.
     $              abs(direction(2)-0.2588).lt.1d-3 .or.
     $              abs(direction(2)-0.9659).lt.1d-3) then
                  class = 'si.v'
               else
                  class = 'si.d'
               endif
            endif
         endif
         rclass = class
      else
         write(6,*) ' sym_op_class_name: unknown op type ', op_type
      endif
c
      end



