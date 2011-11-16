      subroutine sym_op_classify(geom, op_class_index)
C$Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer geom              ! [input]
      integer op_class_index(*) ! [output] Class number of each operator
c     
c     Return an array that has for each operator the number of the
c     class to which it belongs.  This index makes the connection
c     between the operator and the character table.
c
      integer h, g, f, i, j, ir, count(20)
      double precision ff(3,4), hh(3,4), gg(3,4), t1(3,3), t2(3,3), 
     $     diff,gg1(3,3)
      integer maxop, maxireps
      parameter (maxop = 120, maxireps=20)
      character*20 op_type, group_name, op_class_name(maxop)
C      integer op_class_index(maxop)
c
      integer nop, nir,  nop_table, iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      integer order, power, npass
      double precision direction(3)
      logical oagain
c
c     Get character table information including class names
c
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table(group_name, nop_table, nir, iclass_dim,
     $     zir, zclass, chars)) 
     $     call errquit('sym_classify: no char table',0, GEOM_ERR)
c     
c     Verify number of operators match
c
      nop = sym_number_ops(geom)
      if ((nop+1).ne.nop_table) call errquit('sym_classify: ops?',0,
     &       GEOM_ERR)
c     
c     Derive class name of each operator just from the matrices
c
      if (util_print('sym_op_classify',print_never))
     $        write(6,*) 0, ' = E is implied'
      do h = 1, nop
         call sym_get_cart_op(geom, h, hh)
         call sym_op_type(hh, op_type, order, power, direction)
         call sym_op_class_name(op_type, order, power, direction,
     $     group_name, op_class_name(h))
         if (util_print('sym_op_classify',print_never))
     $        write(6,*) h, ' = ', op_class_name(h)
      enddo
c
c     Now having derived the class name from the operator themslves,
c     group the operators into classes using the conjugacy relationship,
c     verify that all ops in a class have the same class name,
c     match this name against the list of standard class names stored
c     for each group, and finally assign to each operator the class
c     index (i.e., the lexical index of its class in the standard order)
c     which is the link needed to make use of the character tables.
c
c     Mmmm ... not quite sure about how best to resolve this.  For at
c     least T and Th there is a problem about picking consistent
c     directions for rotations in sym_op_type ... since C3 and C3^2
c     or S6 and S6^2 are in different classes we need to be consistent.
c     Until I bother to figure out a more elegant way to pick the direction
c     simply relabel classes below if operators are found to be inconsistently
c     labelled ... run the algorithm a second time just to be sure.
c
      npass = 0
c
 1101 do i = 1, nir
         count(i) = 0
      enddo
      count(1) = 1
c
      do h = 1, nop
         op_class_index(h) = 0
      enddo
      oagain = .false.
      npass = npass + 1
c
      do h = 1, nop
         call sym_get_cart_op(geom, h, hh)
         if (op_class_index(h) .eq. 0) then
            if (.not. inp_match(nir, .true., op_class_name(h),
     $           zclass, op_class_index(h))) then
               write(6,*) ' unknown class ', h, op_class_name(h)
               call errquit('sym_clasify: bad class', 0,
     &       GEOM_ERR)
            endif
            count(op_class_index(h)) = count(op_class_index(h)) + 1
            do g = 1, nop
c     
c     Indentify f=gThg
c     
               call sym_get_cart_op(geom, g, gg)
               call dgemm('n','n',3,3,3,1.0d0,hh,3,gg,3,0.0d0,t1,3)
               call dcopy(9,gg,1,gg1,1)
               call dgemm('t','n',3,3,3,1.0d0,gg1,3,t1,3,0.0d0,t2,3)
c     
               do f= 1, nop
                  call sym_get_cart_op(geom, f, ff)
                  diff = 0.0d0
                  do i = 1, 3
                     do j = 1, 3
                        diff = diff + abs(t2(j,i)-ff(j,i))
                     enddo
                  enddo
                  if (diff.lt.1d-6) goto 5
               enddo
               call errquit('sym_classify: out of group', 0, GEOM_ERR)
c     
c     Have f=gThg ... identify class and check for consistency
c     
 5             if (op_class_index(f) .eq. 0) then
                  op_class_index(f) = op_class_index(h)
                  count(op_class_index(f))=count(op_class_index(f))+1
                  if (op_class_name(f) .ne. op_class_name(h)) then
                     if (group_name.eq.'T' .or. group_name.eq.'Th') then
c     See comments above about this switch
                        op_class_name(f) = op_class_name(h)
                        oagain = .true.
                     else
                        write(6,*) ' bad class: h ', h, op_class_name(h)
                        write(6,*) ' bad class: f ', f, op_class_name(f)
                        call errquit('sym_test: bad class', 0, GEOM_ERR)
                     endif
                  endif
               else if (op_class_index(f) .ne. op_class_index(h)) then
                     write(6,*) ' multiple class: h ', h, 
     $                 op_class_name(h)
                     write(6,*) ' multiple class: f ', f, 
     $                    op_class_name(f)
                  call errquit('sym_classify: multiple class', 0,
     &       GEOM_ERR)
               endif
            enddo
         endif
      enddo
      if (oagain .and. npass.eq.1) then
         goto 1101
      else if (oagain .and. npass.ne.1) then
         call errquit('sym_op_classify: patch for T/Th failed',0,
     &       GEOM_ERR)
      endif
c
      do ir = 1, nir
         if (count(ir).ne.iclass_dim(ir)) then
            if (ga_nodeid().eq.0) then
             write(6,*) "Name       Found        Expected"
             do i = 1, nir
              write(6,*) zclass(i),count(i),iclass_dim(i)
             enddo
            endif
            call errquit('sym_classify: classes wrong', ir,
     &       GEOM_ERR)
         endif
      enddo
c
c     If we get here all is consistent
c     
      end
