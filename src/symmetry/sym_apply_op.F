      subroutine sym_apply_op(geom, opnum, r, r1)
C$Id: sym_apply_op.F,v 1.8 1999-07-14 00:29:10 d3e129 Exp $
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
      do i = 1, 3
        r1(i) = sym_ops((opnum-1)*3+i,4,geom)
        do j = 1, 3
          r1(i) = r1(i) + sym_ops((opnum-1)*3+i,j,geom)*r(j)
        enddo
      enddo
c     
c--> shift it into the home cell. This operation depends on sysytem type
c--> in place for 3-d, 2-d and molecules now.
c     
      if(isystype(geom).eq.3) then
 3000    continue
         do 380 i=1,3
            if(dabs(r1(i)).lt.1.0d-10 .or. 
     $           dabs((dabs(r1(i))-1.0d0)).lt.1.0d-10) then
               r1(i)=0.0d0
            endif
            if(r1(i).lt.0.0) then
               r1(i)=r1(i)+1.0d0
            elseif (r1(i).ge.1.0d0) then
               r1(i)=r1(i)-1.0d0
            endif
            if(r1(i).lt.0.0d0.or.r1(i).ge.1.0d0) then
               goto 3000
            endif
 380     continue
      elseif(isystype(geom).eq.2) then ! surfaces
 3020    do 382 i=1,2
            if(dabs(r1(i)).lt.1.0d-10.or.
     $           dabs((dabs(r1(i))-1.0d0)).lt.1.0d-10) then
               r1(i)=0.0d0
            endif
            if(r1(i).lt.0.0d0) then
               r1(i)=r1(i)+1.0d0
            elseif (r1(i).ge.1.0d0) then
               r1(i)=r1(i)-1.0d0
            endif
            if(r1(i).lt.0.0d0.or.r1(i).ge.1.0d0) then
               goto 3020
            endif
 382     continue
      endif
c
*      write(6,*) ' sym_apply_op r ', r
*      write(6,*) ' sym_apply_op op', opnum
*      write(6,*) ' sym_apply_op r1', r1
c
      end
