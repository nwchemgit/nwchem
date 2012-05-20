      subroutine util_sorteps(n, ra, ind, eps)

      implicit none  

*     **** input/output variables ****

      integer n  
      integer ind(n)  
      real*8 ra(n)
      real*8 eps

*     **** local variables ****
      integer i, ir, j, l, iind  
      real*8 rra  

*     **** initialize index array ****
      do i = 1, n  
        ind (i) = i  
      enddo

*     **** nothing to order ****
      if (n.lt.2) return  

      l = n / 2 + 1  
      ir = n  

  10  continue 
  
      if ( l .gt. 1 ) then  
        l    = l - 1  
        rra  = ra (l)  
        iind = ind (l)  
      else  
       rra  = ra (ir)  
       iind = ind (ir)  
       ra (ir) = ra (1)  
       ind (ir) = ind (1)  
       ir = ir - 1  
       if ( ir .eq. 1 ) then  
          ra (1)  = rra  
          ind (1) = iind  
          return 
       endif
      endif
      i = l  
      j = l + l  
      do while ( j .le. ir )  
       if ( j .lt. ir ) then  
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then  
             if (ra(j).lt.ra(j+1)) j = j + 1
          else
             if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
       endif
       if ( abs(rra - ra(j)).ge.eps ) then  
          if (rra.lt.ra(j)) then
             ra (i) = ra (j)  
             ind (i) = ind (j)  
             i = j  
             j = j + j  
          else
             j = ir + 1  
          end if
       else
          if (iind.lt.ind (j) ) then
             ra (i) = ra (j)
             ind (i) = ind (j)
             i = j
             j = j + j
          else
             j = ir + 1
          endif
       end if
      enddo
      ra (i) = rra  
      ind (i) = iind  

      goto 10    

      return
      end 
c $Id$
