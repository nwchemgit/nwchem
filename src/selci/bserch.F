      subroutine selci_insert(n,item,list,nintpo,ipos)
#include "errquit.fh"
*
* $Id: bserch.F,v 1.1 2003-10-17 22:33:53 carlfahl Exp $
*
      integer list(*),item(*)
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
      subroutine selci_bserch(n,item,list,nintpo,ipos)
      integer list(nintpo,*),item(*)
      integer selci_icmp
c
c     binary search thru ordered list of items (each nintpo integers).
c     return in ipos:
c                    if +ve item should precede item at position ipos
c                    if -ve item is already present at position |ipos|
c
      if (n.lt.0)  call errquit('bserch: n.lt.0 ',n, UNKNOWN_ERR)
c
      left = 1
      iright = n
 10   if ((iright-left).le.1) goto 50
      middle=(left+iright)/2
      if (selci_icmp(item,list(1,middle),nintpo)) 20,30,40
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
 50   if (selci_icmp(item,list(1,iright),nintpo)) 60,70,80
c
c item<list(iright)
 60     if (selci_icmp(item,list(1,left), nintpo)) 61,62,63
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
      integer function selci_icmp(item1,item2,n)
      dimension item1(*),item2(*)
c
c     item1 and item2 are packed orbital occupations
c
c     icmp = -1 item1<item2, 0 item1=item2, 1 item1>item2
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
