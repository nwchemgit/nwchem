C$Id$

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    The matrix SYMOPS contains the matrix reps. of all group operators
c   except the identity. The variable NOPS holds the number of operators
c   in SYMOPS.
c
c***********************************************************************
      subroutine gensym_extra(itype,numgrp,numset,symops,nops,oprint,
     >                        group_name)
      implicit none
      integer maxops
      parameter(maxops=192)

      integer itype,numgrp,numset
      real*8 symops(maxops*3,4)
      integer nops
      character*(*) group_name
      logical oprint

*     **** local varialbles ****
      integer iop

      nops = 0

******* Extra Triclinic groups ******
c 
c    A1 
c 
      if (numgrp.eq.231) then
      group_name = 'A1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      
      end if
c 
c    B1 
c 
      if (numgrp.eq.232) then
      group_name = 'B1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      end if
c 
c    C1 
c 
      if (numgrp.eq.233) then
      group_name = 'C1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    F1 
c 
      if (numgrp.eq.234) then
      group_name = 'F1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    I1 
c 
      if (numgrp.eq.235) then
      group_name = 'I1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      end if

c 
c    A-1
c 
      if (numgrp.eq.236) then
      group_name = 'A-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      end if
c 
c    B-1
c 
      if (numgrp.eq.237) then
      group_name = 'B-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      end if
c 
c    C-1 
c 
      if (numgrp.eq.238) then
      group_name = 'C-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    F-1 
c 
      if (numgrp.eq.239) then
      group_name = 'F-1'
      nops = 7
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 3
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 4
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 5
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      iop = 6
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    I-1 
c 
      if (numgrp.eq.240) then
      group_name = 'I-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      end if


      return
      end

