      subroutine sym_op_classify_so(geom, op_class_index, sym_ap_so)
C$Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "inp.fh"
#include "util.fh"
      integer geom              ! [input]
      integer op_class_index(*) ! [output] Class number of each operator
      complex*16 sym_ap_so(*)
c     
c     Return an array that has for each operator the number of the
c     class to which it belongs.  This index makes the connection
c     between the operator and the character table.
c
      integer h, g, f, i, j, ir, count(20)
      double precision ff(3,4), hh(3,4), gg(3,4), t1(3,3), t2(3,3), 
     $     diff
      integer maxop, maxireps
      parameter (maxop = 120, maxireps=20)
      character*20 op_type, group_name, op_class_name(maxop), 
     &     op_rclass_name(maxop)    
      integer op_R_index
c
      integer nop, nir,  nop_table, iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      integer order, power, npass
      double precision direction(3)
      logical oagain
      complex*16 sigma_hh(2,2), sigma_gg(2,2), sigma_ff(2,2),
     &     sigma_t1(2,2), sigma_t2(2,2)
      logical sym_char_table_so
      external sym_char_table_so
c     
c     Get character table information including class names
c
      call sym_group_name(geom,group_name)
      if (.not.sym_char_table_so(group_name,nop_table,nir,iclass_dim,
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
      do h = 1, nop
         call sym_get_cart_op(geom, h, hh)
         call sym_op_type(hh, op_type, order, power, direction)
         call sym_get_spin_op(geom,h,sym_ap_so(4*(h-1)+1))
         call sym_get_spin_op(geom,nop+h,sym_ap_so(4*(nop+h-1)+1))
         call sym_op_class_name_so(op_type, order, power, direction,
     $     group_name, op_class_name(h), op_rclass_name(h))
         if (util_print('sym_op_classify',print_never))
     $        write(6,*) h, ' = ', op_class_name(h)
         write(6,*) h, ' = ', op_class_name(h), op_rclass_name(h)
         op_class_name(nop+h)=op_rclass_name(h)
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
      if (.not.inp_match(nir,.true.,'R',zclass,op_R_index))
     +     call errquit('R must be one of the classes',0, INPUT_ERR)
      count(op_R_index)=1
c
      do h = 1, 2*nop
         op_class_index(h) = 0
      enddo
      oagain = .false.
      npass = npass + 1
c
      do h = 1, 2*nop ! double group has twice as many elements
         call sym_get_cart_op(geom, mod(h-1,nop)+1, hh)
         call sym_get_spin_op(geom, h, sigma_hh)
         if (op_class_index(h) .eq. 0) then
            if (.not. inp_match(nir, .true., op_class_name(h),
     $           zclass, op_class_index(h))) then
               write(6,*) ' unknown class ', h, op_class_name(h)
               call errquit('sym_clasify: bad class', 0, INPUT_ERR)
            endif
            count(op_class_index(h)) = count(op_class_index(h)) + 1
            do g = 1, 2*nop
c     
c     Indentify f=gThg
c     
               call sym_get_cart_op(geom, mod(g-1,nop)+1, gg)
               call dgemm('n','n',3,3,3,1.0d0,hh,3,gg,3,0.0d0,t1,3)
               call dgemm('t','n',3,3,3,1.0d0,gg,3,t1,3,0.0d0,t2,3)
c 
               call sym_get_spin_op(geom, g, sigma_gg)
               call zgemm('n','n',2,2,2,(1.0d0,0.0d0),sigma_hh,2,
     +              sigma_gg,2,(0.0d0,0.0d0),sigma_t1,2)
               call zgemm('c','n',2,2,2,(1.0d0,0.0d0),sigma_gg,2,
     +              sigma_t1,2,(0.0d0,0.0d0),sigma_t2,2)
c     
               do f= 1, 2*nop
                  call sym_get_cart_op(geom, mod(f-1,nop)+1, ff)
                  call sym_get_spin_op(geom, f, sigma_ff)
                  diff = 0.0d0
                  do i = 1, 3
                     do j = 1, 3
                        diff = diff + abs(t2(j,i)-ff(j,i))
                     enddo
                  enddo
                  do i = 1,2
                     do j = 1,2
                        diff=diff+abs(sigma_t2(j,i)-sigma_ff(j,i))
                     enddo
                  enddo
                  if (diff.lt.1d-6) goto 5
               enddo
c               write(*,*)"h,g,f,problem", h,g,f
c               call output(hh, 1, 3, 1, 3, 3, 3, 1)
c               call output(gg, 1, 3, 1, 3, 3, 3, 1)
c               call output(ff, 1, 3, 1, 3, 3, 3, 1)
               call errquit('sym_classify: out of group', 0,
     &       INPUT_ERR)
c     
c     Have f=gThg ... identify class and check for consistency
c     
 5             continue 
c               write(*,*)"h,g,f,hh,gg,gg", h,g,f
c               call output(hh, 1, 3, 1, 3, 3, 3, 1)
c               call output(gg, 1, 3, 1, 3, 3, 3, 1)
c               call output(ff, 1, 3, 1, 3, 3, 3, 1)
               if (op_class_index(f) .eq. 0) then
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
                        call errquit('sym_test: bad class', 0,
     &       INPUT_ERR)
                     endif
                  endif
               else if (op_class_index(f) .ne. op_class_index(h)) then
                     write(6,*) ' multiple class: h ', h, 
     $                 op_class_name(h)
                     write(6,*) ' multiple class: f ', f, 
     $                    op_class_name(f)
                  call errquit('sym_classify: multiple class', 0,
     &       INPUT_ERR)
               endif
            enddo
         endif
      enddo
      if (oagain .and. npass.eq.1) then
         goto 1101
      else if (oagain .and. npass.ne.1) then
         call errquit('sym_op_classify: patch for T/Th failed',0,
     &       INPUT_ERR)
      endif
c
      do ir = 1, nir
         if (count(ir).ne.iclass_dim(ir)) then
            call errquit('sym_classify: classes wrong', ir,
     &       INPUT_ERR)
         endif
      enddo
c
c     If we get here all is consistent
c     
      end

      subroutine sym_get_spin_op(geom,h,sigma) 
      implicit none 
#include "errquit.fh"
#include "sym.fh"
#include "inp.fh"
#include "util.fh"
      character*20 op_type
      double precision phi, direction(3), pi  
      complex*16 sigma_x(2,2), sigma_y(2,2), sigma_z(2,2), iden(2,2), 
     &     sigma(2,2),fact
      double precision hh(3,4), imatrix(3,3),hmatrix(3,3)
      integer h,geom,nop, order, power 
c
      data imatrix 
     +     /-1.0d0,0.0d0,0.0d0,
     +     0.0d0,-1.0d0,0.0d0,
     +     0.0d0,0.0d0,-1.0d0/
c
      sigma(1,1)=dcmplx(0.0d0,0.0d0)
      sigma(2,1)=dcmplx(0.0d0,0.0d0)
      sigma(1,2)=dcmplx(0.0d0,0.0d0)
      sigma(2,2)=dcmplx(0.0d0,0.0d0)
c
      iden(1,1)=dcmplx(1.0d0,0.0d0)
      iden(2,1)=dcmplx(0.0d0,0.0d0)
      iden(1,2)=dcmplx(0.0d0,0.0d0)
      iden(2,2)=dcmplx(1.0d0,0.0d0)
c
      sigma_x(1,1)=dcmplx(0.0d0,0.0d0)
      sigma_x(2,1)=dcmplx(1.0d0,0.0d0)
      sigma_x(1,2)=dcmplx(1.0d0,0.0d0)
      sigma_x(2,2)=dcmplx(0.0d0,0.0d0)
c
      sigma_y(1,1)=dcmplx(0.0d0, 0.0d0)
      sigma_y(2,1)=dcmplx(0.0d0, 1.0d0)
      sigma_y(1,2)=dcmplx(0.0d0,-1.0d0)
      sigma_y(2,2)=dcmplx(0.0d0, 0.0d0)
c
      sigma_z(1,1)=dcmplx(1.0d0,0.0d0)
      sigma_z(2,1)=dcmplx(0.0d0,0.0d0)
      sigma_z(1,2)=dcmplx(0.0d0,0.0d0)
      sigma_z(2,2)=dcmplx(-1.0d0,0.0d0)
c      
      pi = 4.0d0*atan(1.0d0)
      nop = sym_number_ops(geom)
      call sym_get_cart_op(geom, mod(h-1,nop)+1, hh)
      call sym_op_type(hh, op_type, order, power, direction)
      if(op_type.eq.'improper rotation'.or.op_type.eq.'mirror'.or.
     +     op_type.eq.'inversion')then 
         call dgemm('n','n',3,3,3,1.0d0,imatrix,3,hh,3,0.0d0,hmatrix,3)
         call dcopy(9,hmatrix,1,hh,1)
         call sym_op_type(hh, op_type, order, power, direction)
         if(op_type.ne.'rotation'.and.op_type.ne.'identity')
     +        call errquit('sym_get_spin_op: error in imatrix',0,
     &       INPUT_ERR)
      endif
c
c     Now the property of the operation h is found and this information 
c     can be used to construct the corresponding transformation matrix 
c     of this operator on the spin functions alpha and beta according to 
c
c     U(h)=cos(phi/2)I - i sin(phi/2)(nx sigma_x - ny sigma_y - nz sigma_z)  
c
c     where nx, ny and nz are the directions of the rotation and sigma_x, 
c     sigma_y and sigma_z are the pauli matrix. And phi is the rotation 
c     angle
c     
      if(op_type .eq. 'identity')then 
         power = 0 
         order = 1
      endif
      if(power .lt. 0)power=order+power
      phi = -pi*power/order 
      fact = dcmplx(cos(phi),0d0)
      call zaxpy(4, fact, iden,    1, sigma, 1)
      fact = dcmplx(0.0d0, sin(phi)*direction(1))
      call zaxpy(4, fact, sigma_x, 1, sigma, 1)
      fact = dcmplx(0.0d0, sin(phi)*direction(2))
      call zaxpy(4, fact, sigma_y, 1, sigma, 1)
      fact = dcmplx(0.0d0, sin(phi)*direction(3))
      call zaxpy(4, fact, sigma_z, 1, sigma, 1)
      if(h .gt. nop)call zscal(4,(-1.0d0,0.0d0),sigma,1)
      return
      end 




