      subroutine argos_cafe_fsw(xs,xsm,fs,zs,
     + isga,isat,isdt,ismf,isml,isss,isq1,
     + isfrom,nums,xw,xwm,fw,rtos,iwdt,lpbc,lpbcs,esw,esa,
     + vdw,chg,iwatm,iwq,iass,lswndx,lswjpt,lswin,lswj,
     + xi,xj,rwx,rwi1,rwi2,rwi6,rwc,f,fi,fj,facu,
     + rw,isal,isrx,list,dera)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
c
      real*8 xs(msa,3),xsm(msm,3),fs(msa,3,2)
      real*8 zs(msf,3,3,2),esw(msf,mpe,2)
      integer isga(msa),isat(msa),isdt(msa),ismf(msa)
      integer isml(msa),isss(msa),isq1(msa)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),fw(mwm,3,mwa,2),rtos(mwm)
      real*8 esa(nsa)
      integer iwdt(mwm)
      integer isfrom
      logical lpbc,lpbcs
c
      real*8 vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      integer iass(mat,mat),iwatm(mwa),iwq(mwa)
c
      real*8 xi(mscr,3),xj(mscr,3,mwa),rwx(mscr,3)
      real*8 rwi1(mscr),rwi2(mscr),rwi6(mscr),rw(mscr),rwc(mscr,3)
      real*8 f(mscr),fi(mscr,3,mwa),fj(mscr,3,mwa),facu(mscr)
      integer isal(mscr),isrx(mscr)
c
      integer lswj(*)
      integer nums,i
      integer lswndx(0:msa,2),lswjpt(nums,2),lswin(nums,2)
      integer list(0:msa)
c
      real*8 dera(6,nsatot)
c
      integer isatm,nswlen(2)
      integer isfr,iwm,ipsw,number,isa,ispm,isf,nax,ism
      integer ispj,ismn,lswptr,iwa,iax,iwatmi,ix,iy
      integer iwatyp
      real*8 c6,cf6,c12,cf12,sumen
      real*8 c64,c124,dercon,qj,qj4,derco1,derco2
      real*8 drvco1,drvco2,derco3,drvco3,c6p,c12p,etermq,eterml
      real*8 q,qwas,ferfc,fderfc
      real*8 boxi(3),dx,ri1,ri2,ri6,eq,eq0,e6,e8,e12,e14,ff,df
      real*8 xix,xiy,xiz,fact
      real*8 rx,ry,rz,er,p2qi,p3qi,dqi,qi,dd
      integer is,inum,jnum,jwm,iptr,jwa,iss,isg
      logical lid,ljd
      real*8 rtmp
c
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
      etermq=zero
c
      isfr=isfrom-1
      boxi(1)=one/box(1)
      boxi(2)=one/box(2)
      boxi(3)=one/box(3)
      if(npbtyp.eq.1.and.nbxtyp.eq.0.and.icntrl.eq.2) then
      do 101 ipsw=1,npsw
      do 102 is=1,nums
      isa=isfr+is
      ism=isml(isa)
      isf=ismf(isa)
      iss=0
      inum=lswjpt(is,ipsw)-1
      if(iand(isss(isa),6).eq.2) iss=-1
      if(iand(isss(isa),6).eq.4) iss=-1
c      write(*,'(a,4i5)') 'iss ',isa,isss(isa),iand(isss(isa),6),iss
      isg=isga(isa)
      lid=iand(isdt(isa),mdynam).eq.ldynam
      qi=chg(isq1(isa),1,iset)
      iptr=isat(isa)
      xix=xs(isa,1)
      xiy=xs(isa,2)
      xiz=xs(isa,3)
c
      if(lpbc) then
      do 104 ix=1,3
      do 103 jnum=1,lswin(is,ipsw)
      jwm=lswj(inum+jnum)
      dx=xsm(ism,ix)-xwm(jwm,ix)
      rwc(jwm,ix)=dx
      if(abs(dx).gt.boxh(ix)) then
      rwx(jwm,ix)=anint(dx*boxi(ix))*box(ix)
      rwc(jwm,ix)=dx-rwx(jwm,ix)
      else
      rwx(jwm,ix)=zero
      endif
  103 continue
  104 continue
      endif
c
      if(.not.ithint.and..not.ipert2.and..not.ipert3) then
      do 106 jwa=1,nwa
      q=qi*chg(iwq(jwa),1,iset)
      c6=vdw(iptr,iwatm(jwa),1,iset)
      c12=vdw(iptr,iwatm(jwa),3,iset)
      cf6=six*c6
      cf12=twelve*c12
      eq=zero
      eq0=zero
      e6=zero
      e12=zero
      if(ipme.eq.0) then
      do 107 jnum=1,lswin(is,ipsw)
      jwm=lswj(inum+jnum)
      ljd=iand(iwdt(jwm),mdynam).eq.ldynam
      fact=one
      if(.not.lid.or..not.ljd) fact=half
      rx=xix-xw(jwm,1,jwa)-rwx(jwm,1)
      ry=xiy-xw(jwm,2,jwa)-rwx(jwm,2)
      rz=xiz-xw(jwm,3,jwa)-rwx(jwm,3)
      ri2=one/(rx*rx+ry*ry+rz*rz)
      if(iss.gt.0) ri2=one/(one/ri2+shift0(1))
      if(iss.lt.0) ri2=one/(one/ri2+shift1(1))
      ri1=sqrt(ri2)
      ri6=ri2*ri2*ri2
      eq=eq+fact*ri1
      e6=e6+fact*ri6
      e12=e12+fact*ri6*ri6
      ff=(q*ri1+(cf12*ri6-cf6)*ri6)*ri2
      df=ff*rx
      fs(isa,1,ipsw)=fs(isa,1,ipsw)+df
      fw(jwm,1,jwa,ipsw)=fw(jwm,1,jwa,ipsw)-df
      zs(isf,1,1,ipsw)=zs(isf,1,1,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,1,ipsw)=zs(isf,2,1,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,1,ipsw)=zs(isf,3,1,ipsw)-half*df*rwc(jwm,3)
      zw(1,1,ipsw)=zw(1,1,ipsw)-half*df*rwc(jwm,1)
      zw(2,1,ipsw)=zw(2,1,ipsw)-half*df*rwc(jwm,2)
      zw(3,1,ipsw)=zw(3,1,ipsw)-half*df*rwc(jwm,3)
      df=ff*ry
      fs(isa,2,ipsw)=fs(isa,2,ipsw)+df
      fw(jwm,2,jwa,ipsw)=fw(jwm,2,jwa,ipsw)-df
      zs(isf,1,2,ipsw)=zs(isf,1,2,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,2,ipsw)=zs(isf,2,2,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,2,ipsw)=zs(isf,3,2,ipsw)-half*df*rwc(jwm,3)
      zw(1,2,ipsw)=zw(1,2,ipsw)-half*df*rwc(jwm,1)
      zw(2,2,ipsw)=zw(2,2,ipsw)-half*df*rwc(jwm,2)
      zw(3,2,ipsw)=zw(3,2,ipsw)-half*df*rwc(jwm,3)
      df=ff*rz
      fs(isa,3,ipsw)=fs(isa,3,ipsw)+df
      fw(jwm,3,jwa,ipsw)=fw(jwm,3,jwa,ipsw)-df
      zs(isf,1,3,ipsw)=zs(isf,1,3,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,3,ipsw)=zs(isf,2,3,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,3,ipsw)=zs(isf,3,3,ipsw)-half*df*rwc(jwm,3)
      zw(1,3,ipsw)=zw(1,3,ipsw)-half*df*rwc(jwm,1)
      zw(2,3,ipsw)=zw(2,3,ipsw)-half*df*rwc(jwm,2)
      zw(3,3,ipsw)=zw(3,3,ipsw)-half*df*rwc(jwm,3)
  107 continue
      esw(isf,5,ipsw)=esw(isf,5,ipsw)+c12*e12-c6*e6
      esw(isf,6,ipsw)=esw(isf,6,ipsw)+q*eq
      if(npener.ne.0) esa(isg)=esa(isg)+c12*e12-c6*e6+q*eq
      else
      do 108 jnum=1,lswin(is,ipsw)
      jwm=lswj(inum+jnum)
      ljd=iand(iwdt(jwm),mdynam).eq.ldynam
      fact=one
      if(.not.lid.or..not.ljd) fact=half
      rx=xix-xw(jwm,1,jwa)-rwx(jwm,1)
      ry=xiy-xw(jwm,2,jwa)-rwx(jwm,2)
      rz=xiz-xw(jwm,3,jwa)-rwx(jwm,3)
      ri2=one/(rx*rx+ry*ry+rz*rz)
      if(iss.gt.0) ri2=one/(one/ri2+shift0(1))
      if(iss.lt.0) ri2=one/(one/ri2+shift1(1))
      ri1=sqrt(ri2)
      ri6=ri2*ri2*ri2
      er=ealpha/ri1
      ferfc=erfc(er)
      fderfc=ealpha*derfc(er)
      eq0=eq0+fact*ri1
      eq=eq+fact*ri1*ferfc
      e6=e6+fact*ri6
      e12=e12+fact*ri6*ri6
      ff=(q*(ri1*ferfc-fderfc)+(cf12*ri6-cf6)*ri6)*ri2
      df=ff*rx
      fs(isa,1,ipsw)=fs(isa,1,ipsw)+df
      fw(jwm,1,jwa,ipsw)=fw(jwm,1,jwa,ipsw)-df
      zs(isf,1,1,ipsw)=zs(isf,1,1,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,1,ipsw)=zs(isf,2,1,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,1,ipsw)=zs(isf,3,1,ipsw)-half*df*rwc(jwm,3)
      zw(1,1,ipsw)=zw(1,1,ipsw)-half*df*rwc(jwm,1)
      zw(2,1,ipsw)=zw(2,1,ipsw)-half*df*rwc(jwm,2)
      zw(3,1,ipsw)=zw(3,1,ipsw)-half*df*rwc(jwm,3)
      df=ff*ry
      fs(isa,2,ipsw)=fs(isa,2,ipsw)+df
      fw(jwm,2,jwa,ipsw)=fw(jwm,2,jwa,ipsw)-df
      zs(isf,1,2,ipsw)=zs(isf,1,2,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,2,ipsw)=zs(isf,2,2,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,2,ipsw)=zs(isf,3,2,ipsw)-half*df*rwc(jwm,3)
      zw(1,2,ipsw)=zw(1,2,ipsw)-half*df*rwc(jwm,1)
      zw(2,2,ipsw)=zw(2,2,ipsw)-half*df*rwc(jwm,2)
      zw(3,2,ipsw)=zw(3,2,ipsw)-half*df*rwc(jwm,3)
      df=ff*rz
      fs(isa,3,ipsw)=fs(isa,3,ipsw)+df
      fw(jwm,3,jwa,ipsw)=fw(jwm,3,jwa,ipsw)-df
      zs(isf,1,3,ipsw)=zs(isf,1,3,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,3,ipsw)=zs(isf,2,3,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,3,ipsw)=zs(isf,3,3,ipsw)-half*df*rwc(jwm,3)
      zw(1,3,ipsw)=zw(1,3,ipsw)-half*df*rwc(jwm,1)
      zw(2,3,ipsw)=zw(2,3,ipsw)-half*df*rwc(jwm,2)
      zw(3,3,ipsw)=zw(3,3,ipsw)-half*df*rwc(jwm,3)
  108 continue
      esw(isf,5,ipsw)=esw(isf,5,ipsw)+c12*e12-c6*e6
      esw(isf,6,ipsw)=esw(isf,6,ipsw)+q*eq
      if(npener.ne.0) esa(isg)=esa(isg)+c12*e12-c6*e6+q*eq
      endif
  106 continue
      else
      dqi=zero
      p2qi=zero
      p3qi=zero
      if(ipert2) p2qi=chg(isq1(isa),1,2)
      if(ipert3) p3qi=chg(isq1(isa),1,3)
      if(ithint) dqi=chg(isq1(isa),1,4)
      do 116 jwa=1,nwa
      q=qi*chg(iwq(jwa),1,iset)
      c6=vdw(iptr,iwatm(jwa),1,iset)
      c12=vdw(iptr,iwatm(jwa),3,iset)
      cf6=six*c6
      cf12=twelve*c12
      eq=zero
      eq0=zero
      e6=zero
      e8=zero
      e12=zero
      e14=zero
      do 117 jnum=1,lswin(is,ipsw)
      jwm=lswj(inum+jnum)
      ljd=iand(iwdt(jwm),mdynam).eq.ldynam
      fact=one
      if(.not.lid.or..not.ljd) fact=half
      rx=xix-xw(jwm,1,jwa)-rwx(jwm,1)
      ry=xiy-xw(jwm,2,jwa)-rwx(jwm,2)
      rz=xiz-xw(jwm,3,jwa)-rwx(jwm,3)
      ri2=one/(rx*rx+ry*ry+rz*rz)
      if(iss.gt.0) ri2=one/(one/ri2+shift0(1))
      if(iss.lt.0) ri2=one/(one/ri2+shift1(1))
      ri1=sqrt(ri2)
      ri6=ri2*ri2*ri2
      if(ipme.eq.0) then
      eq=eq+fact*ri1
      ff=(q*ri1+(cf12*ri6-cf6)*ri6)*ri2
      else
      er=ealpha/ri1
      ferfc=erfc(er)
      fderfc=ealpha*derfc(er)
      eq=eq+fact*ri1*ferfc
      ff=(q*(ri1*ferfc-fderfc)+(cf12*ri6-cf6)*ri6)*ri2
      endif
      e6=e6+fact*ri6
      e12=e12+fact*ri6*ri6
      if(iss.ne.0) then
      e8=e8+fact*ri6*ri2
      e14=e14+fact*ri6*ri6*ri2
      endif
      df=ff*rx
      fs(isa,1,ipsw)=fs(isa,1,ipsw)+df
      fw(jwm,1,jwa,ipsw)=fw(jwm,1,jwa,ipsw)-df
      zs(isf,1,1,ipsw)=zs(isf,1,1,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,1,ipsw)=zs(isf,2,1,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,1,ipsw)=zs(isf,3,1,ipsw)-half*df*rwc(jwm,3)
      zw(1,1,ipsw)=zw(1,1,ipsw)-half*df*rwc(jwm,1)
      zw(2,1,ipsw)=zw(2,1,ipsw)-half*df*rwc(jwm,2)
      zw(3,1,ipsw)=zw(3,1,ipsw)-half*df*rwc(jwm,3)
      df=ff*ry
      fs(isa,2,ipsw)=fs(isa,2,ipsw)+df
      fw(jwm,2,jwa,ipsw)=fw(jwm,2,jwa,ipsw)-df
      zs(isf,1,2,ipsw)=zs(isf,1,2,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,2,ipsw)=zs(isf,2,2,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,2,ipsw)=zs(isf,3,2,ipsw)-half*df*rwc(jwm,3)
      zw(1,2,ipsw)=zw(1,2,ipsw)-half*df*rwc(jwm,1)
      zw(2,2,ipsw)=zw(2,2,ipsw)-half*df*rwc(jwm,2)
      zw(3,2,ipsw)=zw(3,2,ipsw)-half*df*rwc(jwm,3)
      df=ff*rz
      fs(isa,3,ipsw)=fs(isa,3,ipsw)+df
      fw(jwm,3,jwa,ipsw)=fw(jwm,3,jwa,ipsw)-df
      zs(isf,1,3,ipsw)=zs(isf,1,3,ipsw)-half*df*rwc(jwm,1)
      zs(isf,2,3,ipsw)=zs(isf,2,3,ipsw)-half*df*rwc(jwm,2)
      zs(isf,3,3,ipsw)=zs(isf,3,3,ipsw)-half*df*rwc(jwm,3)
      zw(1,3,ipsw)=zw(1,3,ipsw)-half*df*rwc(jwm,1)
      zw(2,3,ipsw)=zw(2,3,ipsw)-half*df*rwc(jwm,2)
      zw(3,3,ipsw)=zw(3,3,ipsw)-half*df*rwc(jwm,3)
  117 continue
      esw(isf,5,ipsw)=esw(isf,5,ipsw)+c12*e12-c6*e6
      esw(isf,6,ipsw)=esw(isf,6,ipsw)+q*eq
      if(npener.ne.0) esa(isg)=esa(isg)+c12*e12-c6*e6+q*eq
      if(ithint) then
      dd=half*(vdw(iptr,iwatm(jwa),3,4)*e12-vdw(iptr,iwatm(jwa),1,4)*c6)
      if(iss.eq.0) dd=dd+dd
      if(iss.gt.0) dd=dd+shift0(4)*(e8*vdw(iptr,iwatm(jwa),1,4)-
     + e14*vdw(iptr,iwatm(jwa),3,4))
      if(iss.lt.0) dd=dd+shift1(4)*(e8*vdw(iptr,iwatm(jwa),1,4)-
     + e14*vdw(iptr,iwatm(jwa),3,4))
      deriv(3,ipsw)=deriv(3,ipsw)+dd
      deriv(14,ipsw)=deriv(14,ipsw)+dd
      dd=half*(qi*chg(iwq(jwa),1,4)+chg(iwq(jwa),1,iset)*dqi)*eq0
      deriv(5,ipsw)=deriv(5,ipsw)+dd
      deriv(16,ipsw)=deriv(16,ipsw)+dd
      endif
      if(ipert2) then
      ep2(ipsw)=ep2(ipsw)+(p2qi*chg(iwq(jwa),1,2)-q)*eq+
     + (vdw(iptr,iwatm(jwa),3,2)-c12)*e12-
     + (vdw(iptr,iwatm(jwa),1,2)-c6)*e6
      endif
      if(ipert3) then
      ep3(ipsw)=ep3(ipsw)+(p3qi*chg(iwq(jwa),1,3)-q)*eq+
     + (vdw(iptr,iwatm(jwa),3,3)-c12)*e12-
     + (vdw(iptr,iwatm(jwa),1,3)-c6)*e6
      endif
  116 continue
      endif
c
  102 continue
  101 continue
      return
      endif
c
cx new stuff end
c     this subroutine evaluates the solute-solvent forces for nums
c     solute atoms starting from isfrom. the interacting solvent
c     molecules are determined from the pairlist.
c
      isfr=isfrom-1
c
      if(nrwrec.gt.0) then
      do 1 iwm=1,mwm
      rtos(iwm)=zero
    1 continue
      endif
c
c
      do 2 ipsw=1,lpsw
c
c     evaluate outer index array
c
      nswlen(ipsw)=0
      lswndx(0,ipsw)=0
      number=0
      do 3 isa=1,nums
      if(number+lswin(isa,ipsw).gt.mscr .or.
     + (ismf(isfr+isa).ne.ismf(isfr+isa-1).and.
     + number.gt.0)) then
      nswlen(ipsw)=nswlen(ipsw)+1
      lswndx(nswlen(ipsw),ipsw)=isa-1
      number=0
      endif
      number=number+lswin(isa,ipsw)
    3 continue
      if(number.gt.0) then
      nswlen(ipsw)=nswlen(ipsw)+1
      lswndx(nswlen(ipsw),ipsw)=nums
      endif
c
      do 4 ispm=1,nswlen(ipsw)
      isf=ismf(isfr+lswndx(ispm,ipsw))
      do 5 isa=0,nums
      list(isa)=0
    5 continue
      nax=0
c
      do 6 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      ispj=lswjpt(isa,ipsw)-1
      ism=isml(isfr+isa)
      if(lpbc.or.lpbcs.or.ism.eq.0) then
      do 7 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      rwc(nax+ismn,1)=xs(isfr+isa,1)-xwm(lswptr,1)
      rwc(nax+ismn,2)=xs(isfr+isa,2)-xwm(lswptr,2)
      rwc(nax+ismn,3)=xs(isfr+isa,3)-xwm(lswptr,3)
      isrx(nax+ismn)=0
    7 continue
      if(lpbc.or.lpbcs)
     + call argos_cafe_pbc(0,rwc,mscr,rwx,mscr,nax,1,lswin(isa,ipsw))
      endif
      if(ism.gt.0) then
      do 8 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      rwc(nax+ismn,1)=xsm(ism,1)-xwm(lswptr,1)
      rwc(nax+ismn,2)=xsm(ism,2)-xwm(lswptr,2)
      rwc(nax+ismn,3)=xsm(ism,3)-xwm(lswptr,3)
    8 continue
      endif
c
c      if(lssscl) then
c      isrst=iand(isss(isfr+isa),3)
c      isatm=isat(isfr+isa)
c      do 9 iwa=1,mwa
c      iasst=iass(isatm,iwatm(iwa))
c      if(iasst.le.0.or.iasst.ge.3.or.isrst.ne.iasst) isrst=0
c    9 continue
c      do 10 ismn=1,lswin(isa,ipsw)
c      isrx(nax+ismn)=isrst
c   10 continue
c      endif
c
c      write(*,'(4i5,2f12.6)')
c     + lssscl,isga(isa),isss(isfr+isa),iand(isss(isfr+isa),6),
c     + shift0(1),shift1(1)
      if(lssscl) then
      do 10 ismn=1,lswin(isa,ipsw)
c      isrx(nax+ismn)=isss(isfr+isa)
      if(iand(isss(isfr+isa),6).eq.2) isrx(nax+ismn)=-1
      if(iand(isss(isfr+isa),6).eq.4) isrx(nax+ismn)=1
   10 continue
      endif
c
      if(iand(isdt(isfr+isa),mdynam).eq.ldynam) then
      do 11 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xi(nax+ismn,1)=xs(isfr+isa,1)
      xi(nax+ismn,2)=xs(isfr+isa,2)
      xi(nax+ismn,3)=xs(isfr+isa,3)
      isal(nax+ismn)=isfr+isa
c      if(iand(iwdt(lswptr),mdynam).ne.ldynam) then
c      facu(nax+ismn)=half
c      else
      facu(nax+ismn)=one
c      endif
c      if(includ.eq.1) facu(nax+ismn)=one
   11 continue
      else
      do 12 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xi(nax+ismn,1)=xs(isfr+isa,1)
      xi(nax+ismn,2)=xs(isfr+isa,2)
      xi(nax+ismn,3)=xs(isfr+isa,3)
      isal(nax+ismn)=isfr+isa
      if(iand(iwdt(lswptr),mdynam).eq.ldynam) then
      facu(nax+ismn)=one
      else
      facu(nax+ismn)=zero
      endif
      if(includ.eq.1) facu(nax+ismn)=one
   12 continue
      endif
c
      if(.not.lpbc.and..not.lpbcs) then
      do 13 iwa=1,mwa
      do 14 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xj(nax+ismn,1,iwa)=xw(lswptr,1,iwa)
      xj(nax+ismn,2,iwa)=xw(lswptr,2,iwa)
      xj(nax+ismn,3,iwa)=xw(lswptr,3,iwa)
   14 continue
   13 continue
      else
      do 15 ismn=1,lswin(isa,ipsw)
      rwc(nax+ismn,1)=rwc(nax+ismn,1)-rwx(ismn,1)
      rwc(nax+ismn,2)=rwc(nax+ismn,2)-rwx(ismn,2)
      rwc(nax+ismn,3)=rwc(nax+ismn,3)-rwx(ismn,3)
   15 continue
      do 16 iwa=1,mwa
      do 17 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xj(nax+ismn,1,iwa)=xw(lswptr,1,iwa)+rwx(ismn,1)
      xj(nax+ismn,2,iwa)=xw(lswptr,2,iwa)+rwx(ismn,2)
      xj(nax+ismn,3,iwa)=xw(lswptr,3,iwa)+rwx(ismn,3)
   17 continue
   16 continue
      endif
c
      nax=nax+lswin(isa,ipsw)
      list(isa)=nax
    6 continue
c
      do 22 iax=1,nax
      fi(iax,1,1)=zero
      fi(iax,2,1)=zero
      fi(iax,3,1)=zero
   22 continue
      do 23 iwa=1,mwa
      do 24 iax=1,nax
      fj(iax,1,iwa)=zero
      fj(iax,2,iwa)=zero
      fj(iax,3,iwa)=zero
   24 continue
   23 continue
c      if(npener.ne.0) then
c      do 25 iax=1,nax
c      u(iax)=zero
c   25 continue
c      endif
      do 26 iwa=1,mwa
      do 27 iax=1,nax
      f(iax)=zero
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rtmp=rwi2(iax)
      if(isrx(iax).gt.0) rwi2(iax)=one/(one/rwi2(iax)+shift0(1))
      if(isrx(iax).lt.0) rwi2(iax)=one/(one/rwi2(iax)+shift1(1))
c      write(*,'(3i5,2f12.6)')
c     + isga(isal(iax)),isal(iax),isrx(iax),rtmp,rwi2(iax)
   27 continue
c
c     Lennard-Jones interactions
c
      iwatmi=iwatm(iwa)
      eterml=zero
      do 28 iax=1,nax
      isa=isal(iax)
      isatm=isat(isa)
      c6=vdw(isatm,iwatmi,1,iset)
      cf6=six*c6
      c12=vdw(isatm,iwatmi,3,iset)
      cf12=twelve*c12
      rwi6(iax)=rwi2(iax)*rwi2(iax)*rwi2(iax)
      rw(iax)=facu(iax)*(c12*rwi6(iax)-c6)*rwi6(iax)
      eterml=eterml+rw(iax)
      if(npener.ne.0) then
      esa(isga(isa))=esa(isga(isa))+half*rw(iax)
      endif
      f(iax)=f(iax)+(cf12*rwi6(iax)-cf6)*rwi6(iax)*rwi2(iax)
   28 continue
      esw(isf,5,ipsw)=esw(isf,5,ipsw)+eterml
c
c
c     electrostatic interactions
c
      qwas=chg(iwq(iwa),1,iset)
      if(abs(qwas).gt.small.or.ithint.or.
     + (ifstep-1.eq.((ifstep-1)/nfrdf)*nfrdf.and.ngrsw.gt.0)) then
      do 29 iax=1,nax
      rwi1(iax)=sqrt(rwi2(iax))
   29 continue
      endif
      etermq=zero
      if(abs(qwas).gt.small) then
      if(ipme.eq.0) then
      do 30 iax=1,nax
      isa=isal(iax)
      q=qwas*chg(isq1(isa),1,iset)
      rw(iax)=facu(iax)*q*rwi1(iax)
      etermq=etermq+rw(iax)
      if(npener.ne.0) then
      esa(isga(isa))=esa(isga(isa))+half*rw(iax)
      endif
      f(iax)=f(iax)+q*rwi1(iax)*rwi2(iax)
   30 continue
      else
      do 31 iax=1,nax
      isa=isal(iax)
      q=qwas*chg(isq1(isa),1,iset)
      ferfc=erfc(ealpha/rwi1(iax))
      fderfc=ealpha*derfc(ealpha/rwi1(iax))
      rw(iax)=facu(iax)*q*rwi1(iax)
      etermq=etermq+ferfc*rw(iax)
      if(npener.ne.0) then
      esa(isga(isa))=esa(isga(isa))+half*rw(iax)
      endif
      f(iax)=f(iax)+q*rwi2(iax)*(ferfc*rwi1(iax)-fderfc)
   31 continue
      endif
      esw(isf,6,ipsw)=esw(isf,6,ipsw)+etermq
      endif
c
c     reaction field contribution
c
      if(ireact.ne.0) then
      do 32 iax=1,nax
      isa=isal(iax)
      q=qwas*chg(isq1(isa),1,iset)
      rw(iax)=facu(iax)*q*rffsw/rwi2(iax)
      if(npener.ne.0) then
      esa(isga(isa))=esa(isga(isa))+half*q*rffsw/rwi2(iax)
      endif
      f(iax)=f(iax)-two*rffsw*q
   32 continue
      do 33 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      if(list(isa).gt.list(isa-1)) then
      sumen=zero
      do 34 iax=list(isa-1)+1,list(isa)
      sumen=sumen+rw(iax)
   34 continue
      endif
   33 continue
      endif
c
      do 35 iax=1,nax
      fi(iax,1,1)=fi(iax,1,1)+f(iax)*rwx(iax,1)
      fi(iax,2,1)=fi(iax,2,1)+f(iax)*rwx(iax,2)
      fi(iax,3,1)=fi(iax,3,1)+f(iax)*rwx(iax,3)
      fj(iax,1,iwa)=fj(iax,1,iwa)-f(iax)*rwx(iax,1)
      fj(iax,2,iwa)=fj(iax,2,iwa)-f(iax)*rwx(iax,2)
      fj(iax,3,iwa)=fj(iax,3,iwa)-f(iax)*rwx(iax,3)
   35 continue
      do 136 iy=1,3
      do 36 ix=1,3
      sumen=zero
      do 37 iax=1,nax
      sumen=sumen-half*f(iax)*rwx(iax,iy)*rwc(iax,ix)
   37 continue
      zs(isf,ix,iy,ipsw)=zs(isf,ix,iy,ipsw)+sumen
      zw(ix,iy,ipsw)=zw(ix,iy,ipsw)+sumen
   36 continue
  136 continue
c
c     Radial distribution functions
c
c      if(ifstep-1.eq.((ifstep-1)/nfrdf)*nfrdf .and. ngrsw.gt.0) then
c      do 38 igc=1,ngc
c      if(ngt(igc).eq.2) then
c      if(iagc(igc).eq.iwa) then
c      igr=igrc(igc)
c      do 39 iax=1,nax
c      if(isga(isal(iax)).eq.jagc(igc)) then
c      indx=int(one/(rwi1(iax)*drdf))
c      if(indx.gt.ngl) indx=ngl
c      rdf(indx,igr)=rdf(indx,igr)+rdfvol
c      endif
c   39 continue
c      endif
c      endif
c   38 continue
c      endif
c
c     Thermodynamic integration
c
      if(ithint) then
      if(ith(2).or.ith(14)) then
      if(.not.lssscl) then
      do 40 iax=1,nax
      isa=isal(iax)
      isatm=isat(isa)
      c64=vdw(isatm,iwatm(iwa),1,4)
      c124=vdw(isatm,iwatm(iwa),3,4)
      dercon=half*(c124*rwi6(iax)-c64)*rwi6(iax)
      deriv(3,ipsw)=deriv(3,ipsw)+dercon
      deriv(14,ipsw)=deriv(14,ipsw)+dercon
      if(npgdec.gt.1) dera(1,isga(isa))=dera(1,isga(isa))+dercon
   40 continue
      else
      do 41 iax=1,nax
      isa=isal(iax)
      isatm=isat(isa)
      c64=vdw(isatm,iwatm(iwa),1,4)
      c124=vdw(isatm,iwatm(iwa),3,4)
      dercon=half*(c124*rwi6(iax)-c64)*rwi6(iax)
      if(isrx(iax).gt.0) then
      c64=half*three*vdw(isatm,iwatm(iwa),1,iset)
      c124=three*vdw(isatm,iwatm(iwa),3,iset)
      dercon=dercon+shift0(4)*rwi2(iax)*rwi6(iax)*(c64-c124*rwi6(iax))
      elseif(isrx(iax).lt.0) then
      c64=half*three*vdw(isatm,iwatm(iwa),1,iset)
      c124=three*vdw(isatm,iwatm(iwa),3,iset)
      dercon=dercon+shift1(4)*rwi2(iax)*rwi6(iax)*(c64-c124*rwi6(iax))
      else
      c64=vdw(isatm,iwatm(iwa),1,4)
      c124=vdw(isatm,iwatm(iwa),3,4)
      dercon=half*(c124*rwi6(iax)-c64)*rwi6(iax)
      endif
      deriv(3,ipsw)=deriv(3,ipsw)+dercon
      deriv(14,ipsw)=deriv(14,ipsw)+dercon
      if(npgdec.gt.1) dera(1,isga(isa))=dera(1,isga(isa))+dercon
   41 continue
      endif
      endif
      if(ith(4).or.ith(16)) then
      qj=chg(iwq(iwa),1,iset)
      qj4=chg(iwq(iwa),1,4)
      derco1=zero
      derco2=zero
      if(ipme.eq.0) then
      if(.not.lssscl) then
      do 42 iax=1,nax
      isa=isal(iax)
      drvco1=qj*chg(isq1(isa),1,4)*rwi1(iax)
      derco1=derco1+drvco1
      drvco2=chg(isq1(isa),1,iset)*qj4*rwi1(iax)
      derco2=derco2+drvco2
      if(npgdec.gt.1)
     + dera(2,isga(isa))=dera(2,isga(isa))+half*(drvco1+drvco2)
   42 continue
      deriv(5,ipsw)=deriv(5,ipsw)+derco1
      deriv(16,ipsw)=deriv(16,ipsw)+derco2
      else
      derco3=zero
      do 43 iax=1,nax
      isa=isal(iax)
      drvco1=qj*chg(isq1(isa),1,4)*rwi1(iax)
      derco1=derco1+drvco1
      drvco2=chg(isq1(isa),1,iset)*qj4*rwi1(iax)
      derco2=derco2+drvco2
      drvco3=zero
      if(isrx(iax).gt.0) then
      drvco3=(-half)*shift0(4)*chg(isq1(isa),1,iset)*
     + qj*rwi1(iax)*rwi2(iax)
      elseif(isrx(iax).lt.0) then
      drvco3=(-half)*shift1(4)*chg(isq1(isa),1,iset)*
     + qj*rwi1(iax)*rwi2(iax)
      endif
      derco3=derco3+drvco3
      if(npgdec.gt.1) dera(2,isga(isa))=dera(2,isga(isa))+
     + half*(drvco1+drvco2+drvco3)
   43 continue
      deriv(5,ipsw)=deriv(5,ipsw)+derco1+half*derco3
      deriv(16,ipsw)=deriv(16,ipsw)+derco2+half*derco3
      endif
      else
      if(.not.lssscl) then
      do 142 iax=1,nax
      isa=isal(iax)
      drvco1=qj*chg(isq1(isa),1,4)*rwi1(iax)
      derco1=derco1+drvco1
      drvco2=chg(isq1(isa),1,iset)*qj4*rwi1(iax)
      derco2=derco2+drvco2
      if(npgdec.gt.1)
     + dera(2,isga(isa))=dera(2,isga(isa))+half*(drvco1+drvco2)
  142 continue
      deriv(5,ipsw)=deriv(5,ipsw)+derco1
      deriv(16,ipsw)=deriv(16,ipsw)+derco2
      else
      derco3=zero
      do 143 iax=1,nax
      isa=isal(iax)
      drvco1=qj*chg(isq1(isa),1,4)*rwi1(iax)
      derco1=derco1+drvco1
      drvco2=chg(isq1(isa),1,iset)*qj4*rwi1(iax)
      derco2=derco2+drvco2
      drvco3=zero
      if(isrx(iax).gt.0) then
      drvco3=(-half)*shift0(4)*chg(isq1(isa),1,iset)*
     + qj*rwi1(iax)*rwi2(iax)
      elseif(isrx(iax).lt.0) then
      drvco3=(-half)*shift1(4)*chg(isq1(isa),1,iset)*
     + qj*rwi1(iax)*rwi2(iax)
      endif
      derco3=derco3+drvco3
      if(npgdec.gt.1) dera(2,isga(isa))=dera(2,isga(isa))+
     + half*(drvco1+drvco2+drvco3)
  143 continue
      deriv(5,ipsw)=deriv(5,ipsw)+derco1+half*derco3
      deriv(16,ipsw)=deriv(16,ipsw)+derco2+half*derco3
      endif
      endif
      endif
      endif
c
c     Thermodynamic perturbation 1
c
      if(ipert2) then
      if(ip2(2).or.ip2(14)) then
      iwatyp=iwatm(iwa)
      if(.not.lssscl) then
      do 44 iax=1,nax
      isa=isal(iax)
      c6p=vdw(isat(isa),iwatyp,1,2)
      c12p=vdw(isat(isa),iwatyp,3,2)
      ep2(ipsw)=ep2(ipsw)+facu(iax)*(c12p*rwi6(iax)-c6p)*rwi6(iax)
   44 continue
      else
      do 45 iax=1,nax
      isa=isal(iax)
      c6p=vdw(isat(isa),iwatyp,1,2)
      c12p=vdw(isat(isa),iwatyp,3,2)
      if(isrx(iax).gt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift0(1)+shift0(2)))**3
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift1(1)+shift1(2)))**3
      else
      rwi6(iax)=rwi2(iax)**3
      endif
      ep2(ipsw)=ep2(ipsw)+facu(iax)*(c12p*rwi6(iax)-c6p)*rwi6(iax)
   45 continue
      endif
      ep2(ipsw)=ep2(ipsw)-eterml
      endif
      if(ip2(4).or.ip2(5).or.ip2(16).or.ip2(17)) then
      qj=chg(iwq(iwa),1,2)
      if(ipme.eq.0) then
      if(.not.lssscl) then
      do 46 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep2(ipsw)=ep2(ipsw)+facu(iax)*chg(isq1(isa),1,2)*qj*rwi1(iax)
   46 continue
      else
      do 47 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(2))
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(2))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep2(ipsw)=ep2(ipsw)+facu(iax)*chg(isq1(isa),1,2)*qj*rwi1(iax)
   47 continue
      endif
      else
      if(.not.lssscl) then
      do 146 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep2(ipsw)=ep2(ipsw)+facu(iax)*erfc(ealpha/rwi1(iax))*
     + chg(isq1(isa),1,2)*qj*rwi1(iax)
  146 continue
      else
      do 147 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(2))
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(2))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep2(ipsw)=ep2(ipsw)+facu(iax)*erfc(ealpha/rwi1(iax))*
     + chg(isq1(isa),1,2)*qj*rwi1(iax)
  147 continue
      endif
      endif
      ep2(ipsw)=ep2(ipsw)-etermq
      endif
      endif
c
c     Thermodynamic perturbation 2
c
      if(ipert3) then
      if(ip3(2).or.ip3(14)) then
      iwatyp=iwatm(iwa)
      if(.not.lssscl) then
      do 48 iax=1,nax
      isa=isal(iax)
      c6p=vdw(isat(isa),iwatyp,1,3)
      c12p=vdw(isat(isa),iwatyp,3,3)
      ep3(ipsw)=ep3(ipsw)+facu(iax)*(c12p*rwi6(iax)-c6p)*rwi6(iax)
   48 continue
      else
      do 49 iax=1,nax
      isa=isal(iax)
      c6p=vdw(isat(isa),iwatyp,1,3)
      c12p=vdw(isat(isa),iwatyp,3,3)
      if(isrx(iax).gt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift0(1)+shift0(3)))**3
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift1(1)+shift1(3)))**3
      else
      rwi6(iax)=rwi2(iax)**3
      endif
      ep3(ipsw)=ep3(ipsw)+facu(iax)*(c12p*rwi6(iax)-c6p)*rwi6(iax)
   49 continue
      endif
      ep3(ipsw)=ep3(ipsw)-eterml
      endif
      if(ip2(4).or.ip2(5).or.ip2(16).or.ip2(17)) then
      qj=chg(iwq(iwa),1,3)
      if(ipme.eq.0) then
      if(.not.lssscl) then
      do 50 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep3(ipsw)=ep3(ipsw)+facu(iax)*chg(isq1(isa),1,3)*qj*rwi1(iax)
   50 continue
      else
      do 51 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(3))
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(3))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep3(ipsw)=ep3(ipsw)+facu(iax)*chg(isq1(isa),1,3)*qj*rwi1(iax)
   51 continue
      endif
      else
      if(.not.lssscl) then
      do 150 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep3(ipsw)=ep3(ipsw)+facu(iax)*erfc(ealpha/rwi1(iax))*
     + chg(isq1(isa),1,3)*qj*rwi1(iax)
  150 continue
      else
      do 151 iax=1,nax
      isa=isal(iax)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(3))
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(3))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep3(ipsw)=ep3(ipsw)+facu(iax)*erfc(ealpha/rwi1(iax))*
     + chg(isq1(isa),1,3)*qj*rwi1(iax)
  151 continue
      endif
      endif
      ep3(ipsw)=ep3(ipsw)-etermq
      endif
      endif
   26 continue
c
      iax=0
      do 52 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      ispj=lswjpt(isa,ipsw)-1
      do 53 ismn=1,lswin(isa,ipsw)
      fs(isfr+isa,1,ipsw)=fs(isfr+isa,1,ipsw)+fi(iax+ismn,1,1)
      fs(isfr+isa,2,ipsw)=fs(isfr+isa,2,ipsw)+fi(iax+ismn,2,1)
      fs(isfr+isa,3,ipsw)=fs(isfr+isa,3,ipsw)+fi(iax+ismn,3,1)
   53 continue
      do 54 iwa=1,mwa
      do 55 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      fw(lswptr,1,iwa,ipsw)=fw(lswptr,1,iwa,ipsw)+fj(iax+ismn,1,iwa)
      fw(lswptr,2,iwa,ipsw)=fw(lswptr,2,iwa,ipsw)+fj(iax+ismn,2,iwa)
      fw(lswptr,3,iwa,ipsw)=fw(lswptr,3,iwa,ipsw)+fj(iax+ismn,3,iwa)
   55 continue
c
      if(nrwrec.gt.0) then
      do 56 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      if(rtos(lswptr).lt.rwi2(iax+ismn)) rtos(lswptr)=rwi2(iax+ismn)
   56 continue
      endif
   54 continue
c
c      if(npener.ne.0) then
c      do 57 ismn=1,lswin(isa,ipsw)
c      lswptr=lswj(ispj+ismn)
c      uwms(lswptr)=uwms(lswptr)+u(iax+ismn)
c   57 continue
c      endif
c
      iax=iax+lswin(isa,ipsw)
   52 continue
    4 continue
    2 continue
c
      return
      end
