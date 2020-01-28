#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,5,l1), ch(ido,l1,5), wa1(1), wa2(1),
     1  wa3(1), wa4(1), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     2  di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3,
     3  ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c      
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
c      
      idp2 = ido+2
      do 103 k=1,l1
      do i = 1, ((ido - 1)/2)
         ti5 = cc(1+i*2,3,k) + cc(idp2-1-i*2,2,k)
         ti2 = cc(1+i*2,3,k) - cc(idp2-1-i*2,2,k)
         ti4 = cc(1+i*2,5,k) + cc(idp2-1-i*2,4,k)
         ti3 = cc(1+i*2,5,k) - cc(idp2-1-i*2,4,k)
         tr5 = cc(i*2,3,k) - cc(idp2-(i+1)*2,2,k)
         tr2 = cc(i*2,3,k) + cc(idp2-(i+1)*2,2,k)
         tr4 = cc(i*2,5,k) - cc(idp2-(i+1)*2,4,k)
         tr3 = cc(i*2,5,k) + cc(idp2-(i+1)*2,4,k)
         ch(i*2,k,1) = cc(i*2,1,k) + tr2 + tr3
         ch(1+i*2,k,1) = cc(1+i*2,1,k) + ti2 + ti3
         cr2 = cc(i*2,1,k) + tr11*tr2 + tr12*tr3
         ci2 = cc(1+i*2,1,k) + tr11*ti2 + tr12*ti3
         cr3 = cc(i*2,1,k) + tr12*tr2 + tr11*tr3
         ci3 = cc(1+i*2,1,k) + tr12*ti2 + tr11*ti3
         cr5 = ti11*tr5 + ti12*tr4
         ci5 = ti11*ti5 + ti12*ti4
         cr4 = ti12*tr5 - ti11*tr4
         ci4 = ti12*ti5 - ti11*ti4
         dr3 = cr3 - ci4
         dr4 = cr3 + ci4
         di3 = ci3 + cr4
         di4 = ci3 - cr4
         dr5 = cr2 + ci5
         dr2 = cr2 - ci5
         di5 = ci2 - cr5
         di2 = ci2 + cr5
         ch(i*2,k,2) = wa1(i*2-1)*dr2 - wa1(i*2)*di2
         ch(1+i*2,k,2) = wa1(i*2-1)*di2 + wa1(i*2)*dr2
         ch(i*2,k,3) = wa2(i*2-1)*dr3 - wa2(i*2)*di3
         ch(1+i*2,k,3) = wa2(i*2-1)*di3 + wa2(i*2)*dr3
         ch(i*2,k,4) = wa3(i*2-1)*dr4 - wa3(i*2)*di4
         ch(1+i*2,k,4) = wa3(i*2-1)*di4 + wa3(i*2)*dr4
         ch(i*2,k,5) = wa4(i*2-1)*dr5 - wa4(i*2)*di5
         ch(1+i*2,k,5) = wa4(i*2-1)*di5 + wa4(i*2)*dr5
      end do
  103 continue
c      
      return
      end
