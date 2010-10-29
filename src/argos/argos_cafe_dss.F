      subroutine argos_cafe_dss(xs,xsm,ps,psp,ismf,isml,isq2,isq3,
     + isfrom,nums,
     + lpbc,chg,lssndx,lssjpt,lssin,lssj,xi,xj,rwc,rwi1,rwi2,rwi6,
     + rwx,rw,fi,fj,f,isal,jsal,jmal,jfal,pl,pj)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 xs(msa,3),xsm(msm,3),ps(msa,3,2),psp(msa,3,2,2)
      integer ismf(msa),isml(msa),isq2(msa),isq3(msa)
      integer isfrom,nums
      logical lpbc
      real*8 chg(mqt,mqp,mset)
      integer lssndx(0:msa,2),lssjpt(nums,2),lssin(nums,2)
c
      real*8 xi(mscr,3),xj(mscr,3),rwx(mscr,3),rwi1(mscr)
      real*8 rwi2(mscr),rwi6(mscr),rwc(mscr,3),rw(mscr)
      real*8 f(mscr),fi(mscr,3),fj(mscr,3),pl(mscr,3),pj(mscr,3)
      integer isal(mscr),jsal(mscr),jmal(mscr),jfal(mscr)
      integer lssj(*)
c
      integer isa,jsa,jsf,ipset
      integer isfr,nsslen(2)
      integer ipss,number,isslen,nax,nax2,jsaptr
      integer jnum,lssptr,iax,ism,jsm
c
      real*8 ri1,ri2,ri3,rx,ry,rz,pix,piy,piz,pjx,pjy,pjz,rmi,rmj
      real*8 qai,pai,qaj,paj
c
#include "bitops.fh"
c
      real*8 qfaci
c
      qfaci=one/qfac
c
c     solute non-bonded pairs
c     =======================
c
      isfr=isfrom-1
c
c     loop over both short and long range pairlists
c
      do 1 ipss=1,npss
c
c     evaluate outer index array
c
      nsslen(ipss)=0
      lssndx(0,ipss)=0
      number=0
      do 2 isa=1,nums
      if(number+lssin(isa,ipss).gt.mscr .or.
     + (ismf(isfr+isa).ne.ismf(isfr+isa-1).and.number.gt.0)) then
      nsslen(ipss)=nsslen(ipss)+1
      lssndx(nsslen(ipss),ipss)=isa-1
      number=0
      endif
      number=number+lssin(isa,ipss)
    2 continue
      if(number.gt.0) then
      nsslen(ipss)=nsslen(ipss)+1
      lssndx(nsslen(ipss),ipss)=nums
      endif
c
c     loop over number of cycles to complete pairlists
c
      do 3 isslen=1,nsslen(ipss)
      nax=0
      ism=isml(isfr+lssndx(isslen,ipss))
c
c     collect coordinates into workarrays
c
c      if(npbtyp.eq.0) then
      do 4 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      ism=isml(isfr+isa)
      if(lpbc) then
      do 151 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      rwc(nax+jnum,1)=xs(isfr+isa,1)-xs(lssptr,1)
      rwc(nax+jnum,2)=xs(isfr+isa,2)-xs(lssptr,2)
      rwc(nax+jnum,3)=xs(isfr+isa,3)-xs(lssptr,3)
  151 continue
      call argos_cafe_pbc(0,rwc,mscr,rwx,mscr,nax,1,lssin(isa,ipss))
      endif
      do 5 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      jsf=ismf(lssptr)
      isal(nax+jnum)=isfr+isa
      jsal(nax+jnum)=lssptr
      jmal(nax+jnum)=jsf
      jsm=isml(lssptr)
      jfal(nax+jnum)=2
      if(ism.ne.jsm) jfal(nax+jnum)=3
      if(ism.gt.0) then
      if(jsm.gt.0) then
      rwc(nax+jnum,1)=xsm(ism,1)-xsm(jsm,1)
      rwc(nax+jnum,2)=xsm(ism,2)-xsm(jsm,2)
      rwc(nax+jnum,3)=xsm(ism,3)-xsm(jsm,3)
      else
      rwc(nax+jnum,1)=xsm(ism,1)-xs(lssptr,1)
      rwc(nax+jnum,2)=xsm(ism,2)-xs(lssptr,2)
      rwc(nax+jnum,3)=xsm(ism,3)-xs(lssptr,3)
      endif
      else
      if(jsm.gt.0) then
      rwc(nax+jnum,1)=xs(isfr+isa,1)-xsm(jsm,1)
      rwc(nax+jnum,2)=xs(isfr+isa,2)-xsm(jsm,2)
      rwc(nax+jnum,3)=xs(isfr+isa,3)-xsm(jsm,3)
      else
      rwc(nax+jnum,1)=xs(isfr+isa,1)-xs(lssptr,1)
      rwc(nax+jnum,2)=xs(isfr+isa,2)-xs(lssptr,2)
      rwc(nax+jnum,3)=xs(isfr+isa,3)-xs(lssptr,3)
      endif
      endif
    5 continue
      if(.not.lpbc) then
      do 7 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      xi(nax+jnum,1)=xs(isfr+isa,1)
      xi(nax+jnum,2)=xs(isfr+isa,2)
      xi(nax+jnum,3)=xs(isfr+isa,3)
      xj(nax+jnum,1)=xs(lssptr,1)
      xj(nax+jnum,2)=xs(lssptr,2)
      xj(nax+jnum,3)=xs(lssptr,3)
      pl(nax+jnum,1)=ps(isfr+isa,1,2)
      pl(nax+jnum,2)=ps(isfr+isa,2,2)
      pl(nax+jnum,3)=ps(isfr+isa,3,2)
      pj(nax+jnum,1)=ps(lssptr,1,2)
      pj(nax+jnum,2)=ps(lssptr,2,2)
      pj(nax+jnum,3)=ps(lssptr,3,2)
      isal(nax+jnum)=isfr+isa
      jsal(nax+jnum)=lssptr
    7 continue
      else
      do 8 jnum=1,lssin(isa,ipss)
      rwc(nax+jnum,1)=rwc(nax+jnum,1)-rwx(jnum,1)
      rwc(nax+jnum,2)=rwc(nax+jnum,2)-rwx(jnum,2)
      rwc(nax+jnum,3)=rwc(nax+jnum,3)-rwx(jnum,3)
      lssptr=lssj(jsaptr+jnum)
      xi(nax+jnum,1)=xs(isfr+isa,1)
      xi(nax+jnum,2)=xs(isfr+isa,2)
      xi(nax+jnum,3)=xs(isfr+isa,3)
      xj(nax+jnum,1)=xs(lssptr,1)+rwx(jnum,1)
      xj(nax+jnum,2)=xs(lssptr,2)+rwx(jnum,2)
      xj(nax+jnum,3)=xs(lssptr,3)+rwx(jnum,3)
      pl(nax+jnum,1)=ps(isfr+isa,1,2)
      pl(nax+jnum,2)=ps(isfr+isa,2,2)
      pl(nax+jnum,3)=ps(isfr+isa,3,2)
      pj(nax+jnum,1)=ps(lssptr,1,2)
      pj(nax+jnum,2)=ps(lssptr,2,2)
      pj(nax+jnum,3)=ps(lssptr,3,2)
      isal(nax+jnum)=isfr+isa
      jsal(nax+jnum)=lssptr
    8 continue
      endif
      nax=nax+lssin(isa,ipss)
    4 continue
c      else
c      do 8 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
c      jsaptr=lssjpt(isa,ipss)-1
c      ism=isl(isfr+isa,lsmol)
c      do 9 jnum=1,lssin(isa,ipss)
c      lssptr=lssj(jsaptr+jnum)
c      jsa=lssptr
c      rwx(jnum,1)=zero
c      rwx(jnum,2)=zero
c      rwx(jnum,3)=zero
c      if(abs(xs(isfr+isa,1)-xs(jsa,1)).gt.boxh(1)) then
c      rwx(jnum,1)=sign(box(1),xs(isfr+isa,1))
c      endif
c      if(abs(xs(isfr+isa,2)-xs(jsa,2)).gt.boxh(2)) then
c      rwx(jnum,2)=sign(box(2),xs(isfr+isa,2))
c      endif
c      if(npbtyp.eq.1) then
c      if(abs(xs(isfr+isa,3)-xs(jsa,3)).gt.boxh(3)) then
c      rwx(jnum,3)=sign(box(3),xs(isfr+isa,3))
c      endif
c      endif
c      jsf=isl(lssptr,lsfrc)
c      isal(nax+jnum)=isfr+isa
c      jsal(nax+jnum)=lssptr
c      jmal(nax+jnum)=jsf
c      jsm=isl(lssptr,lsmol)
c      jfal(nax+jnum)=2
c      if(ism.ne.jsm) jfal(nax+jnum)=3
c      if(ism.gt.0) then
c      if(jsm.gt.0) then
c      rwc(nax+jnum,1)=xsm(ism,1)-xsm(jsm,1)-rwx(jnum,1)
c      rwc(nax+jnum,2)=xsm(ism,2)-xsm(jsm,2)-rwx(jnum,2)
c      rwc(nax+jnum,3)=xsm(ism,3)-xsm(jsm,3)-rwx(jnum,3)
c      else
c      rwc(nax+jnum,1)=xsm(ism,1)-xs(lssptr,1)-rwx(jnum,1)
c      rwc(nax+jnum,2)=xsm(ism,2)-xs(lssptr,2)-rwx(jnum,2)
c      rwc(nax+jnum,3)=xsm(ism,3)-xs(lssptr,3)-rwx(jnum,3)
c      endif
c      else
c      if(jsm.gt.0) then
c      rwc(nax+jnum,1)=xs(isfr+isa,1)-xsm(jsm,1)-rwx(jnum,1)
c      rwc(nax+jnum,2)=xs(isfr+isa,2)-xsm(jsm,2)-rwx(jnum,2)
c      rwc(nax+jnum,3)=xs(isfr+isa,3)-xsm(jsm,3)-rwx(jnum,3)
c      else
c      rwc(nax+jnum,1)=xs(isfr+isa,1)-xs(lssptr,1)-rwx(jnum,1)
c      rwc(nax+jnum,2)=xs(isfr+isa,2)-xs(lssptr,2)-rwx(jnum,2)
c      rwc(nax+jnum,3)=xs(isfr+isa,3)-xs(lssptr,3)-rwx(jnum,3)
c      endif
c      endif
c    9 continue
c      do 11 jnum=1,lssin(isa,ipss)
c      lssptr=lssj(jsaptr+jnum)
c      xi(nax+jnum,1)=xs(isfr+isa,1)
c      xi(nax+jnum,2)=xs(isfr+isa,2)
c      xi(nax+jnum,3)=xs(isfr+isa,3)
c      xj(nax+jnum,1)=xs(lssptr,1)+rwx(jnum,1)
c      xj(nax+jnum,2)=xs(lssptr,2)+rwx(jnum,2)
c      xj(nax+jnum,3)=xs(lssptr,3)+rwx(jnum,3)
c      pl(nax+jnum,1)=ps(isfr+isa,1,2)
c      pl(nax+jnum,2)=ps(isfr+isa,2,2)
c      pl(nax+jnum,3)=ps(isfr+isa,3,2)
c      pj(nax+jnum,1)=ps(lssptr,1,2)
c      pj(nax+jnum,2)=ps(lssptr,2,2)
c      pj(nax+jnum,3)=ps(lssptr,3,2)
c   11 continue
c      nax=nax+lssin(isa,ipss)
c    8 continue
c      endif
c
      do 12 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jfal(iax).eq.2) then
      qai=qfaci*chg(isq2(isa),1,iset)
      pai=chg(isq2(isa),2,iset)
      qaj=qfaci*chg(isq2(jsa),1,iset)
      paj=chg(isq2(jsa),2,iset)
      else
      qai=qfaci*chg(isq3(isa),1,iset)
      pai=chg(isq3(isa),2,iset)
      qaj=qfaci*chg(isq3(jsa),1,iset)
      paj=chg(isq3(jsa),2,iset)
      endif
      rx=xj(iax,1)-xi(iax,1)
      ry=xj(iax,2)-xi(iax,2)
      rz=xj(iax,3)-xi(iax,3)
      ri2=one/(rx*rx+ry*ry+rz*rz)
      ri1=sqrt(ri2)
      ri3=ri1*ri2
      pix=pai*pl(iax,1)
      piy=pai*pl(iax,2)
      piz=pai*pl(iax,3)
      pjx=paj*pj(iax,1)
      pjy=paj*pj(iax,2)
      pjz=paj*pj(iax,3)
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      fi(iax,1)=((rmj-qaj)*rx-pjx)*ri3
      fi(iax,2)=((rmj-qaj)*ry-pjy)*ri3
      fi(iax,3)=((rmj-qaj)*rz-pjz)*ri3
      fj(iax,1)=((rmi+qai)*rx-pix)*ri3
      fj(iax,2)=((rmi+qai)*ry-piy)*ri3
      fj(iax,3)=((rmi+qai)*rz-piz)*ri3
   12 continue
c
c     accumulate fields into solute field arrays
c
      iax=0
      do 13 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      do 14 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      ps(isfr+isa,1,1)=ps(isfr+isa,1,1)+fi(iax+jnum,1)
      ps(isfr+isa,2,1)=ps(isfr+isa,2,1)+fi(iax+jnum,2)
      ps(isfr+isa,3,1)=ps(isfr+isa,3,1)+fi(iax+jnum,3)
      ps(lssptr,1,1)=ps(lssptr,1,1)+fj(iax+jnum,1)
      ps(lssptr,2,1)=ps(lssptr,2,1)+fj(iax+jnum,2)
      ps(lssptr,3,1)=ps(lssptr,3,1)+fj(iax+jnum,3)
   14 continue
      iax=iax+lssin(isa,ipss)
   13 continue
c
c     thermodynamic integration and perturbation
c
      do 15 ipset=2,3
      if((ipset.eq.2.and.ipert2).or.
     + (ipset.eq.3.and.ipert3)) then
c
      nax2=0
      do 16 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      do 17 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      pl(nax2+jnum,1)=psp(isfr+isa,1,ipset-1,2)
      pl(nax2+jnum,2)=psp(isfr+isa,2,ipset-1,2)
      pl(nax2+jnum,3)=psp(isfr+isa,3,ipset-1,2)
      pj(nax2+jnum,1)=psp(lssptr,1,ipset-1,2)
      pj(nax2+jnum,2)=psp(lssptr,2,ipset-1,2)
      pj(nax2+jnum,3)=psp(lssptr,3,ipset-1,2)
   17 continue
      nax2=nax2+lssin(isa,ipss)
   16 continue
c
      if(nax2.ne.nax) call md_abort('Error in dips',me)
c
      do 18 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jfal(iax).eq.2) then
      qai=qfaci*chg(isq2(isa),1,ipset)
      pai=chg(isq2(isa),2,ipset)
      qaj=qfaci*chg(isq2(jsa),1,ipset)
      paj=chg(isq2(jsa),2,ipset)
      else
      qai=qfaci*chg(isq3(isa),1,ipset)
      pai=chg(isq3(isa),2,ipset)
      qaj=qfaci*chg(isq3(jsa),1,ipset)
      paj=chg(isq3(jsa),2,ipset)
      endif
      rx=xj(iax,1)-xi(iax,1)
      ry=xj(iax,2)-xi(iax,2)
      rz=xj(iax,3)-xi(iax,3)
      ri2=one/(rx*rx+ry*ry+rz*rz)
      ri1=sqrt(ri2)
      ri3=ri1*ri2
      pix=pai*pl(iax,1)
      piy=pai*pl(iax,2)
      piz=pai*pl(iax,3)
      pjx=paj*pj(iax,1)
      pjy=paj*pj(iax,2)
      pjz=paj*pj(iax,3)
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      fi(iax,1)=((rmj-qaj)*rx-pjx)*ri3
      fi(iax,2)=((rmj-qaj)*ry-pjy)*ri3
      fi(iax,3)=((rmj-qaj)*rz-pjz)*ri3
      fj(iax,1)=((rmi+qai)*rx-pix)*ri3
      fj(iax,2)=((rmi+qai)*ry-piy)*ri3
      fj(iax,3)=((rmi+qai)*rz-piz)*ri3
   18 continue
c
c     accumulate fields into solute field arrays
c
      iax=0
      do 19 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      do 20 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      psp(isfr+isa,1,ipset-1,1)=psp(isfr+isa,1,ipset-1,1)+fi(iax+jnum,1)
      psp(isfr+isa,2,ipset-1,1)=psp(isfr+isa,2,ipset-1,1)+fi(iax+jnum,2)
      psp(isfr+isa,3,ipset-1,1)=psp(isfr+isa,3,ipset-1,1)+fi(iax+jnum,3)
      psp(lssptr,1,ipset-1,1)=psp(lssptr,1,ipset-1,1)+fj(iax+jnum,1)
      psp(lssptr,2,ipset-1,1)=psp(lssptr,2,ipset-1,1)+fj(iax+jnum,2)
      psp(lssptr,3,ipset-1,1)=psp(lssptr,3,ipset-1,1)+fj(iax+jnum,3)
   20 continue
      iax=iax+lssin(isa,ipss)
   19 continue
c
      endif
   15 continue
    3 continue
    1 continue
c
      return
      end
