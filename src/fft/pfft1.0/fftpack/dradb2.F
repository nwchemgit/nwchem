#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradb2 (ido,l1,cc,ch,wa1)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,2,l1), ch(ido,l1,2), wa1(1), ti2, tr2
c      
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
c      
      if (ido - 2 .lt. 0) go to 107
      if (ido - 2 .eq. 0) go to 105
  102 continue
      idp2 = ido + 2
      do 104 k=1,l1
      do i = 1, (ido - 1)/2
         ch(i*2,k,1) = cc(i*2,1,k) + cc(idp2-(i+1)*2,2,k)
         tr2 = cc(i*2,1,k) - cc(idp2-(i+1)*2,2,k)
         ch(1+i*2,k,1) = cc(1+i*2,1,k) - cc(idp2-1-i*2,2,k)
         ti2 = cc(1+i*2,1,k) + cc(idp2-1-i*2,2,k)
         ch(i*2,k,2) = wa1(i*2-1)*tr2 - wa1(i*2)*ti2
         ch(1+i*2,k,2) = wa1(i*2-1)*ti2 + wa1(i*2)*tr2
      end do
  104 continue
c      
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
c      
  107 return
      end
