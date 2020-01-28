#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,l1,5), ch(ido,5,l1), wa1(1), wa2(1),
     1  wa3(1), wa4(1), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2,
     2  di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4,
     3  ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c      
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
c      
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
      do i = 1, ((ido - 1)/2)
         dr2 = wa1(i*2-1)*cc(i*2,k,2) + wa1(i*2)*cc(1+i*2,k,2)
         di2 = wa1(i*2-1)*cc(1+i*2,k,2) - wa1(i*2)*cc(i*2,k,2)
         dr3 = wa2(i*2-1)*cc(i*2,k,3) + wa2(i*2)*cc(1+i*2,k,3)
         di3 = wa2(i*2-1)*cc(1+i*2,k,3) - wa2(i*2)*cc(i*2,k,3)
         dr4 = wa3(i*2-1)*cc(i*2,k,4) + wa3(i*2)*cc(1+i*2,k,4)
         di4 = wa3(i*2-1)*cc(1+i*2,k,4) - wa3(i*2)*cc(i*2,k,4)
         dr5 = wa4(i*2-1)*cc(i*2,k,5) + wa4(i*2)*cc(1+i*2,k,5)
         di5 = wa4(i*2-1)*cc(1+i*2,k,5) - wa4(i*2)*cc(i*2,k,5)
         cr2 = dr2 + dr5
         ci5 = dr5 - dr2
         cr5 = di2 - di5
         ci2 = di2 + di5
         cr3 = dr3 + dr4
         ci4 = dr4 - dr3
         cr4 = di3 - di4
         ci3 = di3 + di4
         ch(i*2,1,k) = cc(i*2,k,1) + cr2 + cr3
         ch(1+i*2,1,k) = cc(1+i*2,k,1) + ci2 + ci3
         tr2 = cc(i*2,k,1) + tr11*cr2 + tr12*cr3
         ti2 = cc(1+i*2,k,1) + tr11*ci2 + tr12*ci3
         tr3 = cc(i*2,k,1) + tr12*cr2 + tr11*cr3
         ti3 = cc(1+i*2,k,1) + tr12*ci2 + tr11*ci3
         tr5 = ti11*cr5 + ti12*cr4
         ti5 = ti11*ci5 + ti12*ci4
         tr4 = ti12*cr5 - ti11*cr4
         ti4 = ti12*ci5 - ti11*ci4
         ch(i*2,3,k) = tr2 + tr5
         ch(idp2-(i+1)*2,2,k) = tr2 - tr5
         ch(1+i*2,3,k) = ti2 + ti5
         ch(idp2-1-i*2,2,k) = ti5 - ti2
         ch(i*2,5,k) = tr3 + tr4
         ch(idp2-(i+1)*2,4,k) = tr3 - tr4
         ch(1+i*2,5,k) = ti3 + ti4
         ch(idp2-1-i*2,4,k) = ti4 - ti3
      end do
  103 continue
c      
      return
      end
