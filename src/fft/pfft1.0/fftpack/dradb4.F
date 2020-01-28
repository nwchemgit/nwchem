#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradb4 (ido,l1,cc,ch,wa1,wa2,wa3)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,4,l1), ch(ido,l1,4), wa1(1), wa2(1),
     1  wa3(1), ci2, ci3, ci4, cr2, cr3, cr4, sqrt2, ti1, ti2, ti3, ti4,
     2  tr1, tr2, tr3, tr4
      data sqrt2 /  1.414213562 3730950488 0168872420 970 d0 /
c      
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
c      
      if (ido - 2 .lt. 0) go to 107
      if (ido - 2 .eq. 0) go to 105
  102 continue
      idp2 = ido + 2
      do 104 k=1,l1
      do i = 1, ((ido - 1)/2)
         ti1 = cc(1+i*2,1,k) + cc(idp2-1-i*2,4,k)
         ti2 = cc(1+i*2,1,k) - cc(idp2-1-i*2,4,k)
         ti3 = cc(1+i*2,3,k) - cc(idp2-1-i*2,2,k)
         tr4 = cc(1+i*2,3,k) + cc(idp2-1-i*2,2,k)
         tr1 = cc(i*2,1,k) - cc(idp2-(i+1)*2,4,k)
         tr2 = cc(i*2,1,k) + cc(idp2-(i+1)*2,4,k)
         ti4 = cc(i*2,3,k) - cc(idp2-(i+1)*2,2,k)
         tr3 = cc(i*2,3,k) + cc(idp2-(i+1)*2,2,k)
         ch(i*2,k,1) = tr2 + tr3
         cr3 = tr2 - tr3
         ch(1+i*2,k,1) = ti2 + ti3
         ci3 = ti2 - ti3
         cr2 = tr1 - tr4
         cr4 = tr1 + tr4
         ci2 = ti1 + ti4
         ci4 = ti1 - ti4
         ch(i*2,k,2) = wa1(i*2-1)*cr2 - wa1(i*2)*ci2
         ch(1+i*2,k,2) = wa1(i*2-1)*ci2 + wa1(i*2)*cr2
         ch(i*2,k,3) = wa2(i*2-1)*cr3 - wa2(i*2)*ci3
         ch(1+i*2,k,3) = wa2(i*2-1)*ci3 + wa2(i*2)*cr3
         ch(i*2,k,4) = wa3(i*2-1)*cr4 - wa3(i*2)*ci4
         ch(1+i*2,k,4) = wa3(i*2-1)*ci4 + wa3(i*2)*cr4
      end do
  104 continue
      if (mod(ido,2) .eq. 1) return
c      
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
c      
  107 return
      end
