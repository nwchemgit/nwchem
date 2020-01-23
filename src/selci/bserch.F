C> \ingroup selci
C> @{
C>
C> \brief Insert a new configuration into the list
C>
C> The routine works on a list of items where each item
C> `nintpo` integers long. This routine inserts a new item
C> at the position given by `ipos`. To make space for the new
C> item all subsequent items are shuffled to the right,
C> then the new item is inserted at the specifiedd position.
C>
C> Next the number of elements on the list `n` is increased
C> before the routine returns.
C>
      subroutine selci_insert(n,item,list,nintpo,ipos)
      implicit none
c
#include "errquit.fh"
c
      integer n       !< [In/Output] The number of items on the list
      integer nintpo  !< [Input] The number of integers in an item
      integer ipos    !< [Input] The position for the new item
      integer list(*) !< [In/Output] The list of items
      integer item(*) !< [Input] The new item
*
* $Id$
*
c
c     insert item into the list before position ipos
c     each item is nintpo integers long
c
c     n is incremented before return
c
c     first shuffle the array to the right ... this can
c     be vectorized on some machines but NOT run concurrently
c     without modification
c
c     Local
c
      integer i, ilast, ifirst
c
      ilast = n*nintpo
      ifirst = (ipos-1)*nintpo
      do 10 i = ilast,ifirst+1,-1
         list(i+nintpo) = list(i)
10    continue
c
c     now insert the item
c
cvd$  nodepchk
cvd$  noconcur
      do 20 i = 1,nintpo
          list(i+ifirst) = item(i)
20    continue
c
      n = n + 1
c
      end
C>
C> \brief Find the location for a configuration in the configuration table
C>
C> Perform a binary search through an ordered list of items (each item
C> is `nintpo` integers long). The result is returned in `ipos`.
C> if `ipos` is positive then it represents the element in the list the
C> new item should precede, if `ipos` is negative it is already present
C> in the list.
C>
      subroutine selci_bserch(n,item,list,nintpo,ipos)
      implicit none
c
#include "errquit.fh"
c
      integer n              !< [Input] The length of the list
      integer nintpo         !< [Input] The number of integers per item
      integer list(nintpo,*) !< [Input] The list of items
      integer item(*)        !< [Input] The new item
      integer ipos           !< [Output] The position
                             !< - if ipos > 0: the new item should precede 
                             !<   list(1:nintpo,ipos)
                             !< - if ipos < 0: then item(1:nintpo) equals
                             !<   list(1:nintpo,-ipos)
      integer selci_icmp
      integer*4 isum1
c
c     binary search thru ordered list of items (each nintpo integers).
c     return in ipos:
c                    if +ve item should precede item at position ipos
c                    if -ve item is already present at position |ipos|
c
c     Local
c
      integer middle, left, iright
      integer ifist, ilast, i
c
      if (n.lt.0)  call errquit('bserch: n.lt.0 ',n, UNKNOWN_ERR)
c
      left = 1
      iright = n
 10   if ((iright-left).le.1) goto 50
      middle=(left+iright)/2
      isum1=selci_icmp(item,list(1,middle),nintpo)
      if (isum1) 20,30,40
c
 20   iright = middle
      goto 10
c
 30   ipos = -middle
      return
c
 40   left = middle
      goto 10
c
 50   isum1=selci_icmp(item,list(1,iright),nintpo)
      if (isum1) 60,70,80
c
c item<list(iright)
 60     isum1=selci_icmp(item,list(1,left), nintpo)
        if (isum1) 61,62,63
c
c item<list(left)
 61       ipos = left
          return
c
c item=list(left)
 62       ipos = -left
          return
c
c item>list(left) & item<list(iright)
 63       ipos = iright
          return
c
c item=list(iright)
 70     ipos = -iright
        return
c
c item>list(iright)
 80     ipos = iright + 1
        return
c
      end
C>
C> \brief Compare two configurations
C>
C> Compares two configurations stored in compressed form.
C>
C> \return Returns
C> - -1 if item1 < item2
C> -  0 if item1 = item2
C> -  1 if item1 > item2
C>
      integer function selci_icmp(item1,item2,n)
      implicit none
      integer n        !< [Input] The length of the items in number of integers
      integer item1(n) !< [Input] Item1
      integer item2(n) !< [Input] Item2
c
c     item1 and item2 are packed orbital occupations
c
c     icmp = -1 item1<item2, 0 item1=item2, 1 item1>item2
c
      integer i !< Counter
c
      do 10 i = 1,n
         if (item1(i).gt.item2(i)) then
            selci_icmp = 1
            return
         else if (item1(i).lt.item2(i)) then
            selci_icmp = -1
            return
         endif
10    continue
      selci_icmp = 0
c
      end
C>
C> @}
