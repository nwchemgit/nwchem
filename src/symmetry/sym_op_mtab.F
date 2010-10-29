      subroutine sym_op_mult_table(geom, table, ld)
C$Id$
      implicit none
#include "sym.fh"
#include "inp.fh"
#include "util.fh"
      integer geom              ! [input]
      integer ld
      integer table(ld,*)
c
c  !! THIS ROUTINE HAS NOT BEEN COMPILED OR TESTED !!
c     
c     Return in table the multiplication table for the operators
c     excluding the identity ... inside the table the identity 
c     is labelled as zero
c
      integer h, g, f, i, j
      double precision ff(3,4), hh(3,4), gg(3,4), t1(3,3), diff
c
      do h = 1, nop
         call sym_get_cart_op(geom, h, hh)
         do g = 1, nop
            call sym_get_cart_op(geom, g, gg)
            call dgemm('n','n',3,3,3,1.0d0,hh,3,gg,3,0.0d0,t1,3)
c     
            do f= 0, nop
               if (f .eq. 0) then
                  diff = check for identity
               else
                  call sym_get_cart_op(geom, f, ff)
                  diff = 0.0d0
                  do i = 1, 3
                     do j = 1, 3
                        diff = diff + abs(t1(j,i)-ff(j,i))
                     enddo
                  enddo
               endif
               if (diff.lt.1d-6) goto 5
            enddo
c
            error
c
 5          table(h,g) = f
c
         enddo
      enddo
c
      end

