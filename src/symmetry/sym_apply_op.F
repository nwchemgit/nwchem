      subroutine sym_apply_op(geom, opnum, r, r1)
C$Id$
      implicit none
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom, opnum
      double precision r(3), r1(3)
      integer i, j
c
c     Apply operator to vector of FRACTIONAL coordinates
c
      if (opnum .eq. -9999) then! -9999 means identity(back to box only - abusing this code)
       do i = 1, 3
        r1(i)=r(i)
       enddo
      else if (opnum .lt. 0) then ! Minus means no translations: needed for gradients
       do i = 1, 3
        r1(i)=0
        do j = 1, 3
          r1(i) = r1(i) + sym_ops((-opnum-1)*3+i,j,geom)*r(j)
        enddo
       enddo
      else ! normal stuff
       do i = 1, 3
        r1(i) = sym_ops((opnum-1)*3+i,4,geom)
        do j = 1, 3
          r1(i) = r1(i) + sym_ops((opnum-1)*3+i,j,geom)*r(j)
        enddo
       enddo
      endif
c     
c--> shift it into the home cell. This operation depends on sysytem type
c     
      if(isystype(geom).gt.0) then
 3000    continue
         do 380 i=1,isystype(geom)
            if(dabs(r1(i)).lt.1.0d-10 .or. 
     $           dabs((dabs(r1(i))-1.0d0)).lt.1.0d-10) then
               r1(i)=0.0d0
            endif
            if(r1(i).lt.0.0) then
               r1(i)=r1(i)+1.0d0
               goto 3000
            elseif (r1(i).ge.1.0d0) then
               r1(i)=r1(i)-1.0d0
               goto 3000
            endif
 380     continue
      endif
c
*      write(6,*) ' sym_apply_op r ', r
*      write(6,*) ' sym_apply_op op', opnum
*      write(6,*) ' sym_apply_op r1', r1
c
      end
