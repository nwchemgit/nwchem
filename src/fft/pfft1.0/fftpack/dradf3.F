#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradf3 (ido,l1,cc,ch,wa1,wa2)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,l1,3), ch(ido,3,l1), wa1(1), wa2(1),
     1  ci2, cr2, di2, di3, dr2, dr3, taui, taur, ti2, ti3, tr2, tr3
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c      
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
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
         cr2 = dr2 + dr3
         ci2 = di2 + di3
         ch(i*2,1,k) = cc(i*2,k,1) + cr2
         ch(1+i*2,1,k) = cc(1+i*2,k,1) + ci2
         tr2 = cc(i*2,k,1) + taur*cr2
         ti2 = cc(1+i*2,k,1) + taur*ci2
         tr3 = taui*(di2 - di3)
         ti3 = taui*(dr3 - dr2)
         ch(i*2,3,k) = tr2 + tr3
         ch(idp2-(i+1)*2,2,k) = tr2 - tr3
         ch(1+i*2,3,k) = ti2 + ti3
         ch(idp2-1-i*2,2,k) = ti3 - ti2
      end do
  103 continue
c      
      return
      end
