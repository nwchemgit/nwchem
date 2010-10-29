      subroutine argos_cafe_fw(iwfr,iwto,xw,fw,iwdt,iwatm,iwq,
     + lpbc,eww,vdw,chg,
     + mwb,nwb,nbp,ibnd,bnd,rbnd,mwh,nwh,nhp,iang,ang,rang,rub,
     + mwd,nwd,ndp,idih,dih,rdih,mwo,nwo,nop,iimp,dimp,rimp,
     + mwt,nwt,idwt,mwn,nwn,idwn)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer iwfr,iwto
      integer mwb,mwh,mwd,mwo,nbp,nhp,ndp,nop,mwt,mwn
      integer nwb,nwh,nwd,nwo,nwt,nwn
      real*8 xw(mwm,3,mwa),fw(mwm,3,mwa,2)
      integer iwdt(mwm),iwq(mwa),iwatm(mwa)
      logical lpbc
      real*8 vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      integer ibnd(mwb,3),iang(mwh,4),idih(mwd,5),iimp(mwo,5)
      real*8 bnd(mwb,nbp,mset),ang(mwh,nhp,mset)
      real*8 dih(mwd,ndp,mset),dimp(mwo,nop,mset)
      real*8 rbnd(mwb,2),rang(mwh,2),rub(mwh,2),rdih(mwd,2),rimp(mwo,2)
c
c      real*8 ca6(mat,mat,6),ca12(mat,mat,6)
c      real*8 cb6(mat,mat,6),cb12(mat,mat,6)
c      integer iwl(mwm,miw2),
c
       integer idwt(0:mwt,2),idwn(0:mwn,2)
c
c      real*8 cdwb(mwb,6),ddwb(mwb,6)
c      integer iwbs(mwb),idwb(mwb),jdwb(mwb),iwatm(mwa)
c      real*8 cdwh(mwh,6),ddwh(mwh,6)
c      integer idwh(mwh),jdwh(mwh),kdwh(mwh)
c      real*8 cdwd(mwd,6),ddwd(mwd,6),edwd(mwd,6)
c      integer idwd(mwd),jdwd(mwd),kdwd(mwd),ldwd(mwd)
c      real*8 cdwo(mwo,6),ddwo(mwo,6)
c      integer idwo(mwo),jdwo(mwo),kdwo(mwo),ldwo(mwo)
c      real*8 uwb(mwb),uwh(mwh),uwd(mwd),uwo(mwo)
c
      integer iwb,iwa,jwa,iwm,iwh,kwa,iwd,lwa,iwo,iwt,iwn
      real*8 bond,for,rwx1,rwx2,rwx3,rww,rwwi,dbond,dfor,dfw1,dfw2,dfw3
      real*8 angle,xwij1,xwij2,xwij3,xwkj1,xwkj2,xwkj3,rwij2,rwij2i
      real*8 rwkj2,rwkj2i,cphi,phi,dangle,sphi,rmul
      real*8 xwkl1,xwkl2,xwkl3,xwik1,xwik2,xwik3,xwjl1,xwjl2,xwjl3
      real*8 xm1,xm2,xm3,xn1,xn2,xn3,rm2i,rn2i,rmni,s,rpa
      real*8 xd1,xd2,xd3,xe1,xe2,xe3,dfwi1,dfwi2,dfwi3
      real*8 dfwj1,dfwj2,dfwj3,dfwk1,dfwk2,dfwk3,dfwl1,dfwl2,dfwl3
      real*8 danglep,c6p1,c12p1,c6p2,c12p2,qip1,qjp1,qip2,qjp2
      real*8 c6,c12,c6t,c12t,qit,qjt,cf6,cf12,qi,qj,q,qp1,qp2
      real*8 ep2l,ep3l,ep2q,ep3q,rxx,rxy,rxz,r2,r2i,r1i,r6i,dfw
      real*8 rwikji,sphii,qij,rwi,ferfc,fderfc,eww(mpe,2)
      real*8 etermq,eterml,eub
c
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
      c6t=zero
      c12t=zero
      qit=zero
      qjt=zero
      qp1=zero
      qp2=zero
c
      do 10 iwb=1,nwb
      if(iand(ibnd(iwb,3),icnstr).eq.0) then
      iwa=ibnd(iwb,1)
      jwa=ibnd(iwb,2)
      bond=bnd(iwb,1,iset)
      for=bnd(iwb,2,iset)
      rbnd(iwb,2)=zero
      do 20 iwm=iwfr,iwto
      rwx1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rwx2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rwx3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      rww=sqrt(rwx1**2+rwx2**2+rwx3**2)
      if(rww.lt.tiny) then
      rwwi=one
      else
      rwwi=one/rww
      endif
      dbond=rww-bond
      if(iand(iwdt(iwm),mdynam).eq.ldynam)
     + rbnd(iwb,2)=rbnd(iwb,2)+half*for*(rww-bond)**2
      dfor=for*dbond*rwwi
      dfw1=dfor*rwx1
      dfw2=dfor*rwx2
      dfw3=dfor*rwx3
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfw1
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)+dfw1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfw2
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)+dfw2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfw3
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)+dfw3
      if(ip2(6))
     + ep2(1)=ep2(1)+half*bnd(iwb,2,2)*(rww-bnd(iwb,1,2))**2
      if(ip3(6))
     + ep3(1)=ep3(1)+half*bnd(iwb,2,3)*(rww-bnd(iwb,1,3))**2
      if(ith(6)) then
      deriv(6,1)=deriv(6,1)+
     + dbond*(half*dbond*bnd(iwb,2,4)-for*bnd(iwb,1,4))
      endif
   20 continue
      eww(1,1)=eww(1,1)+rbnd(iwb,2)
      if(ip2(6)) ep2(1)=ep2(1)-rbnd(iwb,2)
      if(ip3(6)) ep3(1)=ep3(1)-rbnd(iwb,2)
      endif
      if(ipme.ne.0) then
      iwa=ibnd(iwb,1)
      jwa=ibnd(iwb,2)
      qij=chg(iwq(iwa),1,iset)*chg(iwq(jwa),1,iset)
      do 21 iwm=iwfr,iwto
      rwx1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rwx2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rwx3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      rww=sqrt(rwx1**2+rwx2**2+rwx3**2)
      rwi=one/rww
      ferfc=one-erfc(ealpha*rww)
      fderfc=-(ealpha*derfc(ealpha*rww))
      epmecw=epmecw-ferfc*qij*rwi
      eww(9,1)=eww(9,1)-ferfc*qij*rwi
      dfor=-(qij*rwi*rwi*(ferfc*rwi-fderfc))
      dfw1=dfor*rwx1
      dfw2=dfor*rwx2
      dfw3=dfor*rwx3
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfw1
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)+dfw1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfw2
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)+dfw2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfw3
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)+dfw3
      vpmeb(1)=vpmeb(1)+dfw1*rwx1
      vpmeb(2)=vpmeb(2)+dfw2*rwx1
      vpmeb(3)=vpmeb(3)+dfw3*rwx1
      vpmeb(4)=vpmeb(4)+dfw2*rwx2
      vpmeb(5)=vpmeb(5)+dfw3*rwx2
      vpmeb(6)=vpmeb(6)+dfw3*rwx3
   21 continue
      endif
   10 continue
      do 40 iwh=1,nwh
      iwa=iang(iwh,1)
      jwa=iang(iwh,2)
      kwa=iang(iwh,3)
      angle=ang(iwh,1,iset)
      for=ang(iwh,2,iset)
      rang(iwh,2)=zero
      do 50 iwm=iwfr,iwto
      xwij1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      xwij2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      xwij3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      xwkj1=xw(iwm,1,kwa)-xw(iwm,1,jwa)
      xwkj2=xw(iwm,2,kwa)-xw(iwm,2,jwa)
      xwkj3=xw(iwm,3,kwa)-xw(iwm,3,jwa)
      rwij2=xwij1**2+xwij2**2+xwij3**2
      rwkj2=xwkj1**2+xwkj2**2+xwkj3**2
      rwij2i=one/rwij2
      rwkj2i=one/rwkj2
      rwikji=one/sqrt(rwij2*rwkj2)
      cphi=rwikji*(xwij1*xwkj1+xwij2*xwkj2+xwij3*xwkj3)
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
      phi=acos(cphi)
      dangle=phi-angle
      if(iand(iwdt(iwm),mdynam).eq.ldynam)
     + rang(iwh,2)=rang(iwh,2)+half*for*dangle*dangle
      sphi=sin(phi)
      if(abs(sphi).lt.small) sphi=small
      dfor=for*dangle/sphi
      dfw1=dfor*(xwkj1*rwikji-xwij1*rwij2i*cphi)
      dfw2=dfor*(xwkj2*rwikji-xwij2*rwij2i*cphi)
      dfw3=dfor*(xwkj3*rwikji-xwij3*rwij2i*cphi)
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)+dfw1
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)-dfw1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)+dfw2
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)-dfw2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)+dfw3
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)-dfw3
      dfw1=dfor*(xwij1*rwikji-xwkj1*rwkj2i*cphi)
      dfw2=dfor*(xwij2*rwikji-xwkj2*rwkj2i*cphi)
      dfw3=dfor*(xwij3*rwikji-xwkj3*rwkj2i*cphi)
      fw(iwm,1,kwa,1)=fw(iwm,1,kwa,1)+dfw1
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)-dfw1
      fw(iwm,2,kwa,1)=fw(iwm,2,kwa,1)+dfw2
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)-dfw2
      fw(iwm,3,kwa,1)=fw(iwm,3,kwa,1)+dfw3
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)-dfw3
      if(ip2(8))
     + ep2(1)=ep2(1)+half*ang(iwh,2,2)*(phi-ang(iwh,1,2))**2
      if(ip3(8))
     + ep3(1)=ep3(1)+half*ang(iwh,2,3)*(phi-ang(iwh,1,3))**2
      if(ith(8)) then
      deriv(8,1)=deriv(8,1)+
     + dangle*(half*dangle*ang(iwh,2,4)-for*ang(iwh,1,4))
      endif
   50 continue
      eww(2,1)=eww(2,1)+rang(iwh,2)
      if(ip2(8)) ep2(1)=ep2(1)-rang(iwh,2)
      if(ip3(8)) ep3(1)=ep3(1)-rang(iwh,2)
      if(ipme.ne.0) then
      iwa=iang(iwh,1)
      jwa=iang(iwh,3)
      qij=chg(iwq(iwa),1,iset)*chg(iwq(jwa),1,iset)
      do 41 iwm=iwfr,iwto
      rwx1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rwx2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rwx3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      rww=sqrt(rwx1**2+rwx2**2+rwx3**2)
      rwi=one/rww
      ferfc=one-erfc(ealpha*rww)
      fderfc=-(ealpha*derfc(ealpha*rww))
      epmecw=epmecw-ferfc*qij*rwi
      eww(9,1)=eww(9,1)-ferfc*qij*rwi
      dfor=-(qij*rwi*rwi*(ferfc*rwi-fderfc))
      dfw1=dfor*rwx1
      dfw2=dfor*rwx2
      dfw3=dfor*rwx3
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfw1
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)+dfw1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfw2
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)+dfw2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfw3
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)+dfw3
      vpmeb(1)=vpmeb(1)+dfw1*rwx1
      vpmeb(2)=vpmeb(2)+dfw2*rwx1
      vpmeb(3)=vpmeb(3)+dfw3*rwx1
      vpmeb(4)=vpmeb(4)+dfw2*rwx2
      vpmeb(5)=vpmeb(5)+dfw3*rwx2
      vpmeb(6)=vpmeb(6)+dfw3*rwx3
   41 continue
      endif
   40 continue
      if(iffld.eq.2) then
      do 1140 iwh=1,nwh
      iwa=iang(iwh,1)
      kwa=iang(iwh,3)
      bond=ang(iwb,3,iset)
      for=ang(iwb,4,iset)
      eub=zero
      do 150 iwm=iwfr,iwto
      rwx1=xw(iwm,1,iwa)-xw(iwm,1,kwa)
      rwx2=xw(iwm,2,iwa)-xw(iwm,2,kwa)
      rwx3=xw(iwm,3,iwa)-xw(iwm,3,kwa)
      rww=sqrt(rwx1**2+rwx2**2+rwx3**2)
      if(rww.lt.tiny) then
      rwwi=one
      else
      rwwi=one/rww
      endif
      dbond=rww-bond
      if(iand(iwdt(iwm),mdynam).eq.ldynam)
     + eub=eub+half*for*(rww-bond)**2
      dfor=for*dbond*rwwi
      dfw1=dfor*rwx1
      dfw2=dfor*rwx2
      dfw3=dfor*rwx3
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfw1
      fw(iwm,1,kwa,1)=fw(iwm,1,kwa,1)+dfw1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfw2
      fw(iwm,2,kwa,1)=fw(iwm,2,kwa,1)+dfw2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfw3
      fw(iwm,3,kwa,1)=fw(iwm,3,kwa,1)+dfw3
      if(ip2(8))
     + ep2(1)=ep2(1)+half*ang(iwh,4,2)*(rww-ang(iwh,3,2))**2
      if(ip3(8))
     + ep3(1)=ep3(1)+half*ang(iwh,4,3)*(rww-ang(iwh,3,3))**2
      if(ith(8)) then
      deriv(8,1)=deriv(8,1)+
     + dbond*(half*dbond*ang(iwh,4,4)-for*ang(iwh,3,4))
      endif
  150 continue
      rub(iwh,2)=rub(iwh,2)+eub
      eww(13,1)=eww(13,1)+eub
      if(ip2(8)) ep2(1)=ep2(1)-eub
      if(ip3(8)) ep3(1)=ep3(1)-eub
 1140 continue
      endif
      do 70 iwd=1,nwd
      iwa=idih(iwd,1)
      jwa=idih(iwd,2)
      kwa=idih(iwd,3)
      lwa=idih(iwd,4)
      angle=dih(iwd,2,iset)
      for=dih(iwd,3,iset)
      rmul=dih(iwd,1,iset)
      rdih(iwd,2)=zero
      do 80 iwm=iwfr,iwto
      xwij1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      xwij2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      xwij3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      xwkj1=xw(iwm,1,kwa)-xw(iwm,1,jwa)
      xwkj2=xw(iwm,2,kwa)-xw(iwm,2,jwa)
      xwkj3=xw(iwm,3,kwa)-xw(iwm,3,jwa)
      xwkl1=xw(iwm,1,kwa)-xw(iwm,1,lwa)
      xwkl2=xw(iwm,2,kwa)-xw(iwm,2,lwa)
      xwkl3=xw(iwm,3,kwa)-xw(iwm,3,lwa)
      xwik1=xwij1-xwkj1
      xwik2=xwij2-xwkj2
      xwik3=xwij3-xwkj3
      xwjl1=xwkl1-xwkj1
      xwjl2=xwkl2-xwkj2
      xwjl3=xwkl3-xwkj3
      xm1=xwij2*xwkj3-xwkj2*xwij3
      xm2=xwij3*xwkj1-xwkj3*xwij1
      xm3=xwij1*xwkj2-xwkj1*xwij2
      xn1=xwkj2*xwkl3-xwkl2*xwkj3
      xn2=xwkj3*xwkl1-xwkl3*xwkj1
      xn3=xwkj1*xwkl2-xwkl1*xwkj2
      rm2i=one/(xm1**2+xm2**2+xm3**2)
      rn2i=one/(xn1**2+xn2**2+xn3**2)
      rmni=sqrt(rm2i*rn2i)
      cphi=(xm1*xn1+xm2*xn2+xm3*xn3)*rmni
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
      phi=acos(cphi)
      s=xwkj1*(xm2*xn3-xm3*xn2) +xwkj2*(xm3*xn1-xm1*xn3)
     + +xwkj3*(xm1*xn2-xm2*xn1)
      if(s.lt.zero) phi=-phi
      sphi=sin(phi)
      rpa=rmul*phi-angle
      if(iand(iwdt(iwm),mdynam).eq.ldynam)
     + rdih(iwd,2)=rdih(iwd,2)+for*(one+cos(rpa))
      dfor=(-for)*rmul*sin(rpa)
      if(ip2(8)) ep2(1)=ep2(1)+
     + dih(iwd,3,2)*(one+cos(dih(iwd,1,2)*phi-dih(iwd,2,2)))
      if(ip3(8)) ep3(1)=ep3(1)+
     + dih(iwd,3,3)*(one+cos(dih(iwd,1,3)*phi-dih(iwd,2,3)))
      if(abs(sphi).lt.small) sphi=sign(small,sphi)
      sphii=one/sphi
      xd1=(-dfor)*sphii*(rmni*xn1-cphi*rm2i*xm1)
      xe1=(-dfor)*sphii*(rmni*xm1-cphi*rn2i*xn1)
      xd2=(-dfor)*sphii*(rmni*xn2-cphi*rm2i*xm2)
      xe2=(-dfor)*sphii*(rmni*xm2-cphi*rn2i*xn2)
      xd3=(-dfor)*sphii*(rmni*xn3-cphi*rm2i*xm3)
      xe3=(-dfor)*sphii*(rmni*xm3-cphi*rn2i*xn3)
      dfwi1=xwkj2*xd3-xwkj3*xd2
      dfwi2=xwkj3*xd1-xwkj1*xd3
      dfwi3=xwkj1*xd2-xwkj2*xd1
      dfwj1=xwik2*xd3-xwik3*xd2-xwkl2*xe3+xwkl3*xe2
      dfwj2=xwik3*xd1-xwik1*xd3-xwkl3*xe1+xwkl1*xe3
      dfwj3=xwik1*xd2-xwik2*xd1-xwkl1*xe2+xwkl2*xe1
      dfwk1=xwjl2*xe3-xwjl3*xe2-xwij2*xd3+xwij3*xd2
      dfwk2=xwjl3*xe1-xwjl1*xe3-xwij3*xd1+xwij1*xd3
      dfwk3=xwjl1*xe2-xwjl2*xe1-xwij1*xd2+xwij2*xd1
      dfwl1=xwkj2*xe3-xwkj3*xe2
      dfwl2=xwkj3*xe1-xwkj1*xe3
      dfwl3=xwkj1*xe2-xwkj2*xe1
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfwi1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfwi2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfwi3
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)-dfwj1
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)-dfwj2
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)-dfwj3
      fw(iwm,1,kwa,1)=fw(iwm,1,kwa,1)-dfwk1
      fw(iwm,2,kwa,1)=fw(iwm,2,kwa,1)-dfwk2
      fw(iwm,3,kwa,1)=fw(iwm,3,kwa,1)-dfwk3
      fw(iwm,1,lwa,1)=fw(iwm,1,lwa,1)-dfwl1
      fw(iwm,2,lwa,1)=fw(iwm,2,lwa,1)-dfwl2
      fw(iwm,3,lwa,1)=fw(iwm,3,lwa,1)-dfwl3
      if(ith(9)) then
      deriv(9,1)=deriv(9,1)+(one+cos(rpa))*dih(iwd,3,4)
     + -for*sin(rpa)*(phi*dih(iwd,1,4)-dih(iwd,2,4))
      endif
   80 continue
      eww(3,1)=eww(3,1)+rdih(iwd,2)
      if(ip2(8)) ep2(1)=ep2(1)-rdih(iwd,2)
      if(ip3(8)) ep3(1)=ep3(1)-rdih(iwd,2)
   70 continue
      do 90 iwo=1,nwo
      iwa=iimp(iwo,1)
      jwa=iimp(iwo,2)
      kwa=iimp(iwo,3)
      lwa=iimp(iwo,4)
      angle=dimp(iwo,2,iset)
      for=dimp(iwo,3,iset)
      rimp(iwo,2)=zero
      do 100 iwm=iwfr,iwto
      xwij1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      xwij2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      xwij3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      xwkj1=xw(iwm,1,kwa)-xw(iwm,1,jwa)
      xwkj2=xw(iwm,2,kwa)-xw(iwm,2,jwa)
      xwkj3=xw(iwm,3,kwa)-xw(iwm,3,jwa)
      xwkl1=xw(iwm,1,kwa)-xw(iwm,1,lwa)
      xwkl2=xw(iwm,2,kwa)-xw(iwm,2,lwa)
      xwkl3=xw(iwm,3,kwa)-xw(iwm,3,lwa)
      xwik1=xwij1-xwkj1
      xwik2=xwij2-xwkj2
      xwik3=xwij3-xwkj3
      xwjl1=xwkl1-xwkj1
      xwjl2=xwkl2-xwkj2
      xwjl3=xwkl3-xwkj3
      xm1=xwij2*xwkj3-xwkj2*xwij3
      xm2=xwij3*xwkj1-xwkj3*xwij1
      xm3=xwij1*xwkj2-xwkj1*xwij2
      xn1=xwkj2*xwkl3-xwkl2*xwkj3
      xn2=xwkj3*xwkl1-xwkl3*xwkj1
      xn3=xwkj1*xwkl2-xwkl1*xwkj2
      rm2i=one/(xm1**2+xm2**2+xm3**2)
      rn2i=one/(xn1**2+xn2**2+xn3**2)
      rmni=sqrt(rm2i*rn2i)
      cphi=(xm1*xn1+xm2*xn2+xm3*xn3)
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
      phi=acos(cphi)
      s=xwkj1*(xm2*xn3-xm3*xn2) +xwkj2*(xm3*xn1-xm1*xn3)
     + +xwkj3*(xm1*xn2-xm2*xn1)
      if(s.lt.zero) phi=-phi
      sphi=sin(phi)
      dangle=(phi-angle)-nint((phi-angle)/twopi)*twopi
      dfor=for*dangle
      if(iand(iwdt(iwm),mdynam).eq.ldynam) rimp(iwo,2)=half*dfor*dangle
      if(ip2(9)) then
      danglep=(phi-dimp(iwo,2,2))-nint((phi-dimp(iwo,2,2))/twopi)*twopi
      ep2(1)=ep2(1)+half*dimp(iwo,3,2)*danglep**2
      endif
      if(ip3(9)) then
      danglep=(phi-dimp(iwo,2,3))-nint((phi-dimp(iwo,2,3))/twopi)*twopi
      ep3(1)=ep3(1)+half*dimp(iwo,3,3)*danglep**2
      endif
      if(abs(sphi).lt.small) sphi=sign(small,sphi)
      sphii=one/sphi
      xd1=(-dfor)*sphii*(rmni*xn1-cphi*rm2i*xm1)
      xe1=(-dfor)*sphii*(rmni*xm1-cphi*rn2i*xn1)
      xd2=(-dfor)*sphii*(rmni*xn2-cphi*rm2i*xm2)
      xe2=(-dfor)*sphii*(rmni*xm2-cphi*rn2i*xn2)
      xd3=(-dfor)*sphii*(rmni*xn3-cphi*rm2i*xm3)
      xe3=(-dfor)*sphii*(rmni*xm3-cphi*rn2i*xn3)
      dfwi1=xwkj2*xd3-xwkj3*xd2
      dfwi2=xwkj3*xd1-xwkj1*xd3
      dfwi3=xwkj1*xd2-xwkj2*xd1
      dfwj1=xwik2*xd3-xwik3*xd2-xwkl2*xe3+xwkl3*xe2
      dfwj2=xwik3*xd1-xwik1*xd3-xwkl3*xe1+xwkl1*xe3
      dfwj3=xwik1*xd2-xwik2*xd1-xwkl1*xe2+xwkl2*xe1
      dfwk1=xwjl2*xe3-xwjl3*xe2-xwij2*xd3+xwij3*xd2
      dfwk2=xwjl3*xe1-xwjl1*xe3-xwij3*xd1+xwij1*xd3
      dfwk3=xwjl1*xe2-xwjl2*xe1-xwij1*xd2+xwij2*xd1
      dfwl1=xwkj2*xe3-xwkj3*xe2
      dfwl2=xwkj3*xe1-xwkj1*xe3
      dfwl3=xwkj1*xe2-xwkj2*xe1
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfwi1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfwi2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfwi3
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)-dfwj1
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)-dfwj2
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)-dfwj3
      fw(iwm,1,kwa,1)=fw(iwm,1,kwa,1)-dfwk1
      fw(iwm,2,kwa,1)=fw(iwm,2,kwa,1)-dfwk2
      fw(iwm,3,kwa,1)=fw(iwm,3,kwa,1)-dfwk3
      fw(iwm,1,lwa,1)=fw(iwm,1,lwa,1)-dfwl1
      fw(iwm,2,lwa,1)=fw(iwm,2,lwa,1)-dfwl2
      fw(iwm,3,lwa,1)=fw(iwm,3,lwa,1)-dfwl3
      if(ith(10)) then
      deriv(10,1)=deriv(10,1)+
     + dangle*(half*dangle*dimp(iwo,3,4)-for*dimp(iwo,2,4))
      endif
  100 continue
      eww(4,1)=eww(4,1)+rimp(iwo,2)
      if(ip2(9)) ep2(1)=ep2(1)-rimp(iwo,2)
      if(ip3(9)) ep3(1)=ep3(1)-rimp(iwo,2)
   90 continue
      c6p1=zero
      c12p1=zero
      c6p2=zero
      c12p2=zero
      qip1=zero
      qjp1=zero
      qip2=zero
      qjp2=zero
      do 110 iwt=1,nwt
      iwa=idwt(iwt,1)
      jwa=idwt(iwt,2)
      c6=vdw(iwatm(iwa),iwatm(jwa),2,iset)
      c12=vdw(iwatm(iwa),iwatm(jwa),4,iset)
      if(ip2(2)) then
      c6p1=vdw(iwatm(iwa),iwatm(jwa),2,2)
      c12p1=vdw(iwatm(iwa),iwatm(jwa),4,2)
      endif
      if(ip3(2)) then
      c6p2=vdw(iwatm(iwa),iwatm(jwa),2,3)
      c12p2=vdw(iwatm(iwa),iwatm(jwa),4,3)
      endif
      if(ith(2).or.ith(4)) then
      c6t=vdw(iwatm(iwa),iwatm(jwa),2,4)
      c12t=vdw(iwatm(iwa),iwatm(jwa),4,4)
      qit=chg(iwq(iwa),1,4)*q14fac
      qjt=chg(iwq(jwa),1,4)
      endif
      cf6=six*c6
      cf12=twelve*c12
      qi=chg(iwq(iwa),1,iset)*q14fac
      qj=chg(iwq(jwa),1,iset)
      q=qi*qj
      if(ip2(4)) then
      qip1=chg(iwq(iwa),1,2)*q14fac
      qjp1=chg(iwq(jwa),1,2)
      qp1=qip1*qjp1
      endif
      if(ip3(4)) then
      qip2=chg(iwq(iwa),1,3)*q14fac
      qjp2=chg(iwq(jwa),1,3)
      qp2=qip2*qjp2
      endif
      ep2l=zero
      ep3l=zero
      ep2q=zero
      ep3q=zero
      do 120 iwm=iwfr,iwto
      rxx=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rxy=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rxz=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      r2=rxx*rxx+rxy*rxy+rxz*rxz
      r2i=one/r2
      r1i=sqrt(r2i)
      r6i=r2i*r2i*r2i
      eterml=(c12*r6i-c6)*r6i
      etermq=q*r1i
      if(iand(iwdt(iwm),mdynam).eq.ldynam) eww(5,1)=eww(5,1)+eterml
      if(iand(iwdt(iwm),mdynam).eq.ldynam) eww(6,1)=eww(6,1)+etermq
      if(ip2(2)) ep2l=ep2l-eterml+(c12p1*r6i-c6p1)*r6i
      if(ip3(2)) ep3l=ep3l-eterml+(c12p2*r6i-c6p2)*r6i
      if(ip2(4)) ep2q=ep2q-etermq+qp1*r1i
      if(ip3(4)) ep3q=ep3q-etermq+qp2*r1i
      dfw=((cf12*r6i-cf6)*r6i+q*r1i)*r2i
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)+dfw*rxx
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)+dfw*rxy
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)+dfw*rxz
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)-dfw*rxx
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)-dfw*rxy
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)-dfw*rxz
      if(ith(2)) then
      deriv(2,1)=deriv(2,1)+(c12t*r6i-c6t)*r6i
      endif
      if(ith(4)) then
      deriv(4,1)=deriv(4,1)+(qi*qjt+qj*qit)*r1i
      endif
  120 continue
      ep2(1)=ep2(1)+ep2l+ep2q
      ep3(1)=ep3(1)+ep3l+ep3q
      if(ipme.ne.0) then
      qij=(one-q14fac)*chg(iwq(iwa),1,iset)*chg(iwq(jwa),1,iset)
      do 111 iwm=iwfr,iwto
      rwx1=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rwx2=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rwx3=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      rww=sqrt(rwx1**2+rwx2**2+rwx3**2)
      rwi=one/rww
      ferfc=one-erfc(ealpha*rww)
      fderfc=-(ealpha*derfc(ealpha*rww))
      epmecw=epmecw-ferfc*qij*rwi
      eww(6,1)=eww(6,1)-ferfc*qij*rwi
      dfor=-(qij*rwi*rwi*(ferfc*rwi-fderfc))
      dfw1=dfor*rwx1
      dfw2=dfor*rwx2
      dfw3=dfor*rwx3
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)-dfw1
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)+dfw1
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)-dfw2
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)+dfw2
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)-dfw3
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)+dfw3
      vpmeb(1)=vpmeb(1)+dfw1*rwx1
      vpmeb(2)=vpmeb(2)+dfw2*rwx1
      vpmeb(3)=vpmeb(3)+dfw3*rwx1
      vpmeb(4)=vpmeb(4)+dfw2*rwx2
      vpmeb(5)=vpmeb(5)+dfw3*rwx2
      vpmeb(6)=vpmeb(6)+dfw3*rwx3
  111 continue
      endif
  110 continue
      do 130 iwn=1,nwn
      iwa=idwn(iwn,1)
      jwa=idwn(iwn,2)
      c6=vdw(iwatm(iwa),iwatm(jwa),1,iset)
      c12=vdw(iwatm(iwa),iwatm(jwa),3,iset)
      if(ip2(2)) then
      c6p1=vdw(iwatm(iwa),iwatm(jwa),1,2)
      c12p1=vdw(iwatm(iwa),iwatm(jwa),3,2)
      endif
      if(ip3(2)) then
      c6p2=vdw(iwatm(iwa),iwatm(jwa),1,3)
      c12p2=vdw(iwatm(iwa),iwatm(jwa),3,3)
      endif
      if(ith(2).or.ith(4)) then
      c6t=vdw(iwatm(iwa),iwatm(jwa),1,4)
      c12t=vdw(iwatm(iwa),iwatm(jwa),3,4)
      qit=chg(iwq(iwa),1,4)
      qjt=chg(iwq(jwa),1,4)
      endif
      cf6=six*c6
      cf12=twelve*c12
      qi=chg(iwq(iwa),1,iset)
      qj=chg(iwq(jwa),1,iset)
      q=qi*qj
      if(ip2(4)) then
      qip1=chg(iwq(iwa),1,2)
      qjp1=chg(iwq(jwa),1,2)
      qp1=qip1*qjp1
      endif
      if(ip3(4)) then
      qip2=chg(iwq(iwa),1,3)
      qjp2=chg(iwq(jwa),1,3)
      qp2=qip2*qjp2
      endif
      ep2l=zero
      ep3l=zero
      ep2q=zero
      ep3q=zero
      do 140 iwm=iwfr,iwto
      rxx=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rxy=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rxz=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      r2=rxx*rxx+rxy*rxy+rxz*rxz
      r2i=one/r2
      r1i=sqrt(r2i)
      r6i=r2i*r2i*r2i
      ferfc=one
      fderfc=zero
      if(ipme.ne.0) then
      ferfc=erfc(ealpha/r1i)
      fderfc=ealpha+derfc(ealpha/r1i)
      endif
      eterml=(c12*r6i-c6)*r6i
      etermq=ferfc*q*r1i
      if(iand(iwdt(iwm),mdynam).eq.ldynam) then
      eww(5,1)=eww(5,1)+eterml
      eww(6,1)=eww(6,1)+etermq
      endif
      if(ip2(2)) ep2l=ep2l-eterml+(c12p1*r6i-c6p1)*r6i
      if(ip3(2)) ep3l=ep3l-eterml+(c12p2*r6i-c6p2)*r6i
      if(ip2(4)) ep2q=ep2q-etermq+qp1*r1i
      if(ip3(4)) ep3q=ep3q-etermq+qp2*r1i
      dfw=((cf12*r6i-cf6)*r6i+q*(ferfc*r1i-fderfc))*r2i
      fw(iwm,1,iwa,1)=fw(iwm,1,iwa,1)+dfw*rxx
      fw(iwm,2,iwa,1)=fw(iwm,2,iwa,1)+dfw*rxy
      fw(iwm,3,iwa,1)=fw(iwm,3,iwa,1)+dfw*rxz
      fw(iwm,1,jwa,1)=fw(iwm,1,jwa,1)-dfw*rxx
      fw(iwm,2,jwa,1)=fw(iwm,2,jwa,1)-dfw*rxy
      fw(iwm,3,jwa,1)=fw(iwm,3,jwa,1)-dfw*rxz
      if(ith(2)) deriv(2,1)=deriv(2,1)+(c12t*r6i-c6t)*r6i
      if(ith(4)) deriv(4,1)=deriv(4,1)+(qi*qjt+qj*qit)*r1i
  140 continue
      ep2(1)=ep2(1)+ep2l+ep2q
      ep3(1)=ep3(1)+ep3l+ep3q
  130 continue
c
      return
      end
