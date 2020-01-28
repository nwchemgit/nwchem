#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradb3 (ido,l1,cc,ch,wa1,wa2)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,3,l1), ch(ido,l1,3), wa1(1), wa2(1),
     1  ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c      
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
c      
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
      do i = 1, ((ido - 1)/2)
         tr2 = cc(i*2,3,k) + cc(idp2-(i+1)*2,2,k)
         cr2 = cc(i*2,1,k) + taur*tr2
         ch(i*2,k,1) = cc(i*2,1,k) + tr2
         ti2 = cc(1+i*2,3,k) - cc(idp2-1-i*2,2,k)
         ci2 = cc(1+i*2,1,k) + taur*ti2
         ch(1+i*2,k,1) = cc(1+i*2,1,k) + ti2
         cr3 = taui*(cc(i*2,3,k)-cc(idp2-(i+1)*2,2,k))
         ci3 = taui*(cc(1+i*2,3,k)+cc(idp2-1-i*2,2,k))
         dr2 = cr2 - ci3
         dr3 = cr2 + ci3
         di2 = ci2 + cr3
         di3 = ci2 - cr3
         ch(i*2,k,2) = wa1(i*2-1)*dr2 - wa1(i*2)*di2
         ch(1+i*2,k,2) = wa1(i*2-1)*di2 + wa1(i*2)*dr2
         ch(i*2,k,3) = wa2(i*2-1)*dr3 - wa2(i*2)*di3
         ch(1+i*2,k,3) = wa2(i*2-1)*di3 + wa2(i*2)*dr3
      end do
  103 continue
c      
      return
      end
