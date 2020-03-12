#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradf4 (ido,l1,cc,ch,wa1,wa2,wa3)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,l1,4), ch(ido,4,l1), wa1(1), wa2(1),
     1  wa3(1), ci2, ci3, ci4, cr2, cr3, cr4, hsqt2, ti1, ti2, ti3,
     2  ti4, tr1, tr2, tr3, tr4
      data hsqt2 /   .7071067811 8654752440 0844362104 85 d0 /
c      
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
c      
      if (ido - 2 .lt. 0) go to 107
      if (ido - 2 .eq. 0) go to 105
  102 continue
      idp2 = ido + 2
      do 104 k=1,l1
      do i = 1, ((ido - 1)/2)
         cr2 = wa1(i*2-1)*cc(i*2,k,2) + wa1(i*2)*cc(1+i*2,k,2)
         ci2 = wa1(i*2-1)*cc(1+i*2,k,2) - wa1(i*2)*cc(i*2,k,2)
         cr3 = wa2(i*2-1)*cc(i*2,k,3) + wa2(i*2)*cc(1+i*2,k,3)
         ci3 = wa2(i*2-1)*cc(1+i*2,k,3) - wa2(i*2)*cc(i*2,k,3)
         cr4 = wa3(i*2-1)*cc(i*2,k,4) + wa3(i*2)*cc(1+i*2,k,4)
         ci4 = wa3(i*2-1)*cc(1+i*2,k,4) - wa3(i*2)*cc(i*2,k,4)
         tr1 = cr2 + cr4
         tr4 = cr4 - cr2
         ti1 = ci2 + ci4
         ti4 = ci2 - ci4
         ti2 = cc(1+i*2,k,1) + ci3
         ti3 = cc(1+i*2,k,1) - ci3
         tr2 = cc(i*2,k,1) + cr3
         tr3 = cc(i*2,k,1) - cr3
         ch(i*2,1,k) = tr1 + tr2
         ch(idp2-(i+1)*2,4,k) = tr2 - tr1
         ch(1+i*2,1,k) = ti1 + ti2
         ch(idp2-1-i*2,4,k) = ti1 - ti2
         ch(i*2,3,k) = ti4 + tr3
         ch(idp2-(i+1)*2,2,k) = tr3 - ti4
         ch(1+i*2,3,k) = tr4 + ti3
         ch(idp2-1-i*2,2,k) = tr4 - ti3
      end do
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
c      
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
c      
  107 return
      end
