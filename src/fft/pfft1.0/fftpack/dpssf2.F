      subroutine dpssf2 (ido,l1,cc,ch,wa1)
*
* $Id: dpssf2.F,v 1.2 1997-11-04 10:07:22 d3e129 Exp $
*
      double precision cc(ido,2,l1), ch(ido,l1,2), wa1(1), ti2, tr2
c
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
  103    continue
  104 continue
c
      return
      end
