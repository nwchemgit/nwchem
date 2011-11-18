      subroutine sym_op_class_name(op_type, order, power, direction,
     $     group_name, class)
C$Id$
      implicit none
#include "inp.fh"
c     
      character*(*) op_type     ! [input]
      integer order, power      ! [input]
      double precision direction(3) ! [input]
      character*(*) group_name  ! [input]
      character*(*) class       ! [output]
c     
      integer abspow
c     
c     Try to assign the class name of an operator given all the other info
c     
c     Powers of rotations in Cn, Cnh and Sn are all in separate classes
c     Also for T and Th
c     
      abspow = abs(power)
      if (group_name(1:1).eq.'S' .or. 
     $     group_name.eq.'T' .or. group_name.eq.'Th' .or.
     $     (group_name(1:1).eq.'C' .and.
     $     (group_name(3:3).eq.'h' .or. inp_strlen(group_name).eq.2)))
     $     then
         if (power .lt. 0) then
            if (op_type.eq.'rotation' .or. mod(order,2).eq.0) then
               abspow = order+power
            else
               abspow = 2*order+power ! Odd improper rotation
            endif
         endif
      endif
c
      class = ' '
      if (op_type .eq. 'identity') then
         class = 'e'
      else if (op_type .eq. 'inversion') then
         if (group_name.eq.'Ci' .or. group_name.eq.'S2') then
            class = 'i=s2'
         else
            class = 'i'
         endif
      else if (op_type .eq. 'rotation') then
         if (order .gt. 2) then
            if (abspow .gt. 1) then
               write(class,'(''c'',i1,''^'',i1)') order, abspow
            else
               write(class,'(''c'',i1)') order
            endif
         else
            class = 'c2'
            if (group_name.eq.'D2' .or. group_name.eq.'D2h') then
               if (abs(direction(3)-1.0d0) .lt. 1.0e-6) then
                  class = 'c2z'
               else if (abs(direction(2)-1.0d0) .lt. 1.0e-6) then
                  class = 'c2y'
               else if (abs(direction(1)-1.0d0) .lt. 1.0e-6) then
                  class = 'c2x'
               endif
            else if (group_name.eq.'D4' .or. group_name.eq.'D4h') then
               if (abs(direction(3)-1.0d0) .lt. 1d-6) then
                  class = 'c2'
               else if ((abs(direction(1)-1.0d0) .lt. 1d-6) .or.
     $                 (abs(direction(2)-1.0d0) .lt. 1d-6)) then
                  class = 'c2'''
               else
                  class = 'c2"'
               endif
            else if (group_name.eq.'D2d' .or. group_name.eq.'D4d' .or. 
     $              group_name.eq.'D6d') then
               if (abs(direction(3)-1.0d0) .gt. 1.0e-6) then
                  class = 'c2'''
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
            else if (group_name(1:2).eq.'D8') then
c
c     Added by Kurt on 16 NOV 2011
c
               if (abs(direction(3)-1.0d0).lt.1e-3) then
                  class = 'c2'
               else if (
     $                 abs(direction(1)-0.7071).lt.1d-3 .or.
     $                 abs(direction(1)-1.0000).lt.1d-3 .or.
     $                 abs(direction(2)-1.0000).lt.1d-3) then
                  class = 'c2'''
               else
                  class = 'c2"'
               endif
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
            endif
         endif
      else if (op_type .eq. 'improper rotation') then
         if (group_name(1:1) .eq. 'D') then
c
c     Additional mirror planes of Dnh force some ops into the same class
c     whereas for Cnh they are all in different classes.
c
            if (abspow .gt. order) abspow = 2*order - abspow
         endif
         if (abspow .gt. 1) then
            if (order .lt. 10) then
               write(class,'(''s'',i1,''^'',i1)') order, abspow
            else
               write(class,'(''s'',i2,''^'',i1)') order, abspow
            endif
         else
            if (order .lt. 10) then
               write(class,'(''s'',i1)') order
            else
               write(class,'(''s'',i2)') order
            endif
         endif
      else if (op_type .eq. 'mirror') then
         if (group_name.eq.'C2v') then
            if (abs(direction(1)-1.0d0) .lt. 1.0d-6) then
               class = 'si.yz'
            else  if (abs(direction(2)-1.0d0) .lt. 1.0d-6) then
               class = 'si.xz'
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

            else if (group_name.eq.'C8v' .or. group_name.eq.'D8h') then
c
c    Added by Kurt 16 Nov 2011
c
               if (
     $              abs(direction(1)-0.7071).lt.1d-3 .or.
     $              abs(direction(1)-1.0000).lt.1d-3 .or.
     $              abs(direction(2)-1.0000).lt.1d-3) then
                  class = 'si.d'
               else
                  class = 'si.v'
               endif
            endif
         endif
      else
         write(6,*) ' sym_op_class_name: unknown op type ', op_type
      endif
c
c      print*,op_type, order, power, direction(1:3),group_name,
c     &     class,acos(direction(1))*180/3.1415926535897932384626d0

      end
