      subroutine sym_print_ops(geom)
C$Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "inp.fh"
      integer geom
c
      double precision r(3,4), direction(3)
      integer nops, op, i, j, order, power
      integer maxireps, maxops
      parameter (maxireps=20, maxops=120)
      character*20 group_name, op_class_name, op_type
      integer op_class_index(maxops), lclass
c
c     Print out detailed info on the group operations
c
      call sym_op_classify(geom, op_class_index)
      call sym_group_name(geom,group_name)
      nops = sym_number_ops(geom)
c
      do op = 1, nops
         call sym_get_cart_op(geom, op, r)
         call sym_op_type(r, op_type, order, power, direction)
         call sym_op_class_name(op_type, order, power, direction,
     $     group_name, op_class_name)
         lclass = inp_strlen(op_class_name)
         if (op_type .eq. 'mirror') then
            write(6,3) op, op_class_name(1:lclass),(direction(i),i=1,3)
 3          format('     op ', i3, ' = ', a, ' = mirror normal to (',
     $           f5.2,',',f5.2,',',f5.2,')')
         else if (op_type .eq. 'inversion') then
            write(6,4) op, op_class_name(1:lclass)
 4          format('     op ', i3, ' = ', a, ' = inversion')
         else if (op_type .eq. 'rotation' .or.
     $            op_type .eq. 'improper rotation') then
            write(6,5) op, op_class_name(1:lclass),
     $           op_type(1:inp_strlen(op_type)), (direction(i),i=1,3)
 5          format('     op ', i3, ' = ', a, ' = ',a,
     $           ' about axis (',
     $           f5.2,',',f5.2,',',f5.2,')')
         else
            call errquit('sym_print_ops: unknown op?',op, INPUT_ERR)
         endif
c
         do i = 1, 3
            write(6,6) (r(i,j),j=1,4)
 6          format(10x,4f12.6)
         enddo
c
      enddo
c
      end
