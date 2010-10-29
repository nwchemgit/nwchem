      subroutine argos_cafe_dsw(xs,xsm,ps,psp,isdt,ismf,isml,
     + isq1,isfrom,nums,
     + xw,xwm,pw,pwp,iwq,lpbc,chg,lswndx,lswjpt,lswin,lswj,
     + rwc,xi,xj,rwx,pl,pj,fi,fj,isal)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 xs(msa,3),xsm(msm,3),ps(msa,3,2),psp(msa,3,2,2)
      integer isdt(msa),ismf(msa),isml(msa),isq1(msa)
      integer isfrom,nums
      real*8 xw(mwm,3,mwa),xwm(mwm,3),pw(mwm,3,mwa,2),pwp(mwm,3,mwa,2,2)
      integer iwq(mwa)
      logical lpbc
      real*8 chg(mqt,mqp,mset)
      integer lswndx(0:msa,2),lswjpt(nums,2),lswin(nums,2),lswj(*)
      real*8 rwc(mscr,3),xi(mscr,3),xj(mscr,3,mwa),rwx(mscr,3)
      real*8 pl(mscr,3),pj(mscr,3,mwa),fi(mscr,3),fj(mscr,3,mwa)
      integer isal(mscr)
c
      integer ispj,ism,lswptr,ipset,nswlen(2)
      integer isfr,ipsw,number,isa,jwa,ismn,ispm,iax,nax,nax2
      real*8 qai,qaj,pai,paj,rx,ry,rz,pix,piy,piz,pjx,pjy,pjz
      real*8 ri2,ri3,rmi,rmj
#include "bitops.fh"
c
      real*8 qfaci
c
      qfaci=one/qfac
c
c     this subroutine evaluates the solute-solvent forces for nums
c     solute atoms starting from isfrom. the interacting solvent
c     molecules are determined from the pairlist.
c
      isfr=isfrom-1
c
c     loop over short and long range pairs
c
      do 1 ipsw=1,npsw
c
c     evaluate outer index array
c
      nswlen(ipsw)=0
      lswndx(0,ipsw)=0
      number=0
      do 2 isa=1,nums
      if(number+lswin(isa,ipsw).gt.mscr .or.
     + (ismf(isfr+isa).ne.ismf(isfr+isa-1).and.number.gt.0)) then
      nswlen(ipsw)=nswlen(ipsw)+1
      lswndx(nswlen(ipsw),ipsw)=isa-1
      number=0
      endif
      number=number+lswin(isa,ipsw)
    2 continue
      if(number.gt.0) then
      nswlen(ipsw)=nswlen(ipsw)+1
      lswndx(nswlen(ipsw),ipsw)=nums
      endif
c
c     loop over number of cycles to complete pairlist
c
      do 3 ispm=1,nswlen(ipsw)
      nax=0
c
c     vacuo conditions
c
c      if(npbtyp.eq.0) then
c      if(.not.lpbc) then
      do 4 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      ispj=lswjpt(isa,ipsw)-1
      ism=isml(isfr+isa)
c
c     collect center of mass distance vectors
c
      if(lpbc.or.ism.eq.0) then
      do 6 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      rwc(nax+ismn,1)=xs(isfr+isa,1)-xwm(lswptr,1)
      rwc(nax+ismn,2)=xs(isfr+isa,2)-xwm(lswptr,2)
      rwc(nax+ismn,3)=xs(isfr+isa,3)-xwm(lswptr,3)
    6 continue
      if(lpbc)
     + call argos_cafe_pbc(0,rwc,mscr,rwx,mscr,nax,1,lswin(isa,ipsw))
      else
      do 5 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      rwc(nax+ismn,1)=xsm(ism,1)-xwm(lswptr,1)
      rwc(nax+ismn,2)=xsm(ism,2)-xwm(lswptr,2)
      rwc(nax+ismn,3)=xsm(ism,3)-xwm(lswptr,3)
    5 continue
      endif
c
c     collect solute coordinates and atomic polarization fields
c
c      if(iand(isdt(isfr+isa),mdynam).eq.ldynam) then
      do 7 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xi(nax+ismn,1)=xs(isfr+isa,1)
      xi(nax+ismn,2)=xs(isfr+isa,2)
      xi(nax+ismn,3)=xs(isfr+isa,3)
      pl(nax+ismn,1)=ps(isfr+isa,1,2)
      pl(nax+ismn,2)=ps(isfr+isa,2,2)
      pl(nax+ismn,3)=ps(isfr+isa,3,2)
      isal(nax+ismn)=isfr+isa
    7 continue
c      else
c      do 8 ismn=1,lswin(isa,ipsw)
c      lswptr=lswj(ispj+ismn)
c      xi(nax+ismn,1)=xs(isfr+isa,1)
c      xi(nax+ismn,2)=xs(isfr+isa,2)
c      xi(nax+ismn,3)=xs(isfr+isa,3)
c      pl(nax+ismn,1)=ps(isfr+isa,1,2)
c      pl(nax+ismn,2)=ps(isfr+isa,2,2)
c      pl(nax+ismn,3)=ps(isfr+isa,3,2)
c      isal(nax+ismn)=isfr+isa
c    8 continue
c      endif
c
c     collect solvent coordinates and atomic polarization fields
c
      do 8 jwa=1,mwa
      if(lpbc) then
      do 9 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xj(nax+ismn,1,jwa)=xw(lswptr,1,jwa)+rwx(ismn,1)
      xj(nax+ismn,2,jwa)=xw(lswptr,2,jwa)+rwx(ismn,2)
      xj(nax+ismn,3,jwa)=xw(lswptr,3,jwa)+rwx(ismn,3)
      pj(nax+ismn,1,jwa)=pw(lswptr,1,jwa,2)
      pj(nax+ismn,2,jwa)=pw(lswptr,2,jwa,2)
      pj(nax+ismn,3,jwa)=pw(lswptr,3,jwa,2)
    9 continue
      else
      do 10 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      xj(nax+ismn,1,jwa)=xw(lswptr,1,jwa)
      xj(nax+ismn,2,jwa)=xw(lswptr,2,jwa)
      xj(nax+ismn,3,jwa)=xw(lswptr,3,jwa)
      pj(nax+ismn,1,jwa)=pw(lswptr,1,jwa,2)
      pj(nax+ismn,2,jwa)=pw(lswptr,2,jwa,2)
      pj(nax+ismn,3,jwa)=pw(lswptr,3,jwa,2)
   10 continue
      endif
    8 continue
      nax=nax+lswin(isa,ipsw)
    4 continue
c      else
cc
cc     periodic boundary conditions
cc
c      do 11 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
c      ispj=lswjpt(isa,ipsw)-1
c      ism=isml(isfr+isa)
cc
cc     collect center of mass distance vectors
cc
c      do 12 ismn=1,lswin(isa,ipsw)
c      lswptr=lswj(ispj+ismn)
c      rwc(nax+ismn,1)=xs(isfr+isa,1)-xwm(lswptr,1)
c      rwc(nax+ismn,2)=xs(isfr+isa,2)-xwm(lswptr,2)
c      rwc(nax+ismn,3)=xs(isfr+isa,3)-xwm(lswptr,3)
c      rwx(ismn,1)=zero
c      rwx(ismn,2)=zero
c      rwx(ismn,3)=zero
c      if(abs(rwc(nax+ismn,1)).gt.boxh(1)) then
c      rwx(ismn,1)=sign(box(1),xs(isfr+isa,1))
c      endif
c      if(abs(rwc(nax+ismn,2)).gt.boxh(2)) then
c      rwx(ismn,2)=sign(box(2),xs(isfr+isa,2))
c      endif
c      if(npbtyp.eq.1) then
c      if(abs(rwc(nax+ismn,3)).gt.boxh(3)) then
c      rwx(ismn,3)=sign(box(3),xs(isfr+isa,3))
c      endif
c      endif
c      if(ism.gt.0) then
c      rwc(nax+ismn,1)=xsm(ism,1)-xwm(lswptr,1)-rwx(ismn,1)
c      rwc(nax+ismn,2)=xsm(ism,2)-xwm(lswptr,2)-rwx(ismn,2)
c      rwc(nax+ismn,3)=xsm(ism,3)-xwm(lswptr,3)-rwx(ismn,3)
c      else
c      rwc(nax+ismn,1)=xs(isfr+isa,1)-xwm(lswptr,1)-rwx(ismn,1)
c      rwc(nax+ismn,2)=xs(isfr+isa,2)-xwm(lswptr,2)-rwx(ismn,2)
c      rwc(nax+ismn,3)=xs(isfr+isa,3)-xwm(lswptr,3)-rwx(ismn,3)
c      endif
c   12 continue
cc
cc     collect solute coordinates and atomic polarization fields
cc
c      do 13 ismn=1,lswin(isa,ipsw)
c      lswptr=lswj(ispj+ismn)
c      xi(nax+ismn,1)=xs(isfr+isa,1)
c      xi(nax+ismn,2)=xs(isfr+isa,2)
c      xi(nax+ismn,3)=xs(isfr+isa,3)
c      pl(nax+ismn,1)=ps(isfr+isa,1,2)
c      pl(nax+ismn,2)=ps(isfr+isa,2,2)
c      pl(nax+ismn,3)=ps(isfr+isa,3,2)
c      isal(nax+ismn)=isfr+isa
c   13 continue
cc
cc     collect solvent coordinates and atomic polarization fields
cc
c      do 14 jwa=1,mwa
c      do 15 ismn=1,lswin(isa,ipsw)
c      lswptr=lswj(ispj+ismn)
c      xj(nax+ismn,1,jwa)=xw(lswptr,1,jwa)+rwx(ismn,1)
c      xj(nax+ismn,2,jwa)=xw(lswptr,2,jwa)+rwx(ismn,2)
c      xj(nax+ismn,3,jwa)=xw(lswptr,3,jwa)+rwx(ismn,3)
c      pj(nax+ismn,1,jwa)=pw(lswptr,1,jwa,2)
c      pj(nax+ismn,2,jwa)=pw(lswptr,2,jwa,2)
c      pj(nax+ismn,3,jwa)=pw(lswptr,3,jwa,2)
c   15 continue
c   14 continue
c      nax=nax+lswin(isa,ipsw)
c   11 continue
c      endif
c
c     zero temparary arays fi and fj
c
      do 16 iax=1,nax
      fi(iax,1)=zero
      fi(iax,2)=zero
      fi(iax,3)=zero
   16 continue
      do 17 jwa=1,mwa
      do 18 iax=1,nax
      fj(iax,1,jwa)=zero
      fj(iax,2,jwa)=zero
      fj(iax,3,jwa)=zero
   18 continue
   17 continue
c
c     loop over the number of atoms in a solvent molecule
c
c     calculated here is 4*pi*epsilon*field and not just the field
c     since the polarization is given in alpha/(4*pi*epsilon) in
c     stead of just alpha, the induced dipole is obtained by the
c     product of pwa and pw
c
      do 19 jwa=1,mwa
      qaj=qfaci*chg(iwq(jwa),1,iset)
      paj=chg(iwq(jwa),2,iset)
      do 20 iax=1,nax
      isa=isal(iax)
      qai=qfaci*chg(isq1(isa),1,iset)
      pai=chg(isq1(isa),2,iset)
      rx=xj(iax,1,jwa)-xi(iax,1)
      ry=xj(iax,2,jwa)-xi(iax,2)
      rz=xj(iax,3,jwa)-xi(iax,3)
      pix=pai*pl(iax,1)
      piy=pai*pl(iax,2)
      piz=pai*pl(iax,3)
      pjx=paj*pj(iax,1,jwa)
      pjy=paj*pj(iax,2,jwa)
      pjz=paj*pj(iax,3,jwa)
      ri2=one/(rx**2+ry**2+rz**2)
      ri3=sqrt(ri2)*ri2
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      fi(iax,1)=fi(iax,1)+((rmj-qaj)*rx-pjx)*ri3
      fi(iax,2)=fi(iax,2)+((rmj-qaj)*ry-pjy)*ri3
      fi(iax,3)=fi(iax,3)+((rmj-qaj)*rz-pjz)*ri3
      fj(iax,1,jwa)=fj(iax,1,jwa)+((rmi+qai)*rx-pix)*ri3
      fj(iax,2,jwa)=fj(iax,2,jwa)+((rmi+qai)*ry-piy)*ri3
      fj(iax,3,jwa)=fj(iax,3,jwa)+((rmi+qai)*rz-piz)*ri3
   20 continue
   19 continue
c
c     update the electric field arrays
c
      iax=0
      do 21 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      ispj=lswjpt(isa,ipsw)-1
      do 22 ismn=1,lswin(isa,ipsw)
      ps(isfr+isa,1,1)=ps(isfr+isa,1,1)+fi(iax+ismn,1)
      ps(isfr+isa,2,1)=ps(isfr+isa,2,1)+fi(iax+ismn,2)
      ps(isfr+isa,3,1)=ps(isfr+isa,3,1)+fi(iax+ismn,3)
   22 continue
      do 23 jwa=1,mwa
      do 24 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      pw(lswptr,1,jwa,1)=pw(lswptr,1,jwa,1)+fj(iax+ismn,1,jwa)
      pw(lswptr,2,jwa,1)=pw(lswptr,2,jwa,1)+fj(iax+ismn,2,jwa)
      pw(lswptr,3,jwa,1)=pw(lswptr,3,jwa,1)+fj(iax+ismn,3,jwa)
   24 continue
   23 continue
      iax=iax+lswin(isa,ipsw)
   21 continue
c
c     thermodynamic integration and perturbation
c
      do 30 ipset=2,3
      if((ipset.eq.2.and.ipert2).or.
     + (ipset.eq.3.and.ipert3)) then
c
      nax2=0
      do 31 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      ispj=lswjpt(isa,ipsw)-1
      do 32 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      pl(nax2+ismn,1)=psp(isfr+isa,1,2,ipset-1)
      pl(nax2+ismn,2)=psp(isfr+isa,2,2,ipset-1)
      pl(nax2+ismn,3)=psp(isfr+isa,3,2,ipset-1)
   32 continue
      do 33 jwa=1,mwa
      do 34 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      pj(nax2+ismn,1,jwa)=pwp(lswptr,1,jwa,2,ipset-1)
      pj(nax2+ismn,2,jwa)=pwp(lswptr,2,jwa,2,ipset-1)
      pj(nax2+ismn,3,jwa)=pwp(lswptr,3,jwa,2,ipset-1)
   34 continue
   33 continue
      nax2=nax2+lswin(isa,ipsw)
   31 continue
c
      if(nax.ne.nax2) call md_abort('Error in dipsw',me)
c
      do 41 iax=1,nax
      fi(iax,1)=zero
      fi(iax,2)=zero
      fi(iax,3)=zero
   41 continue
      do 42 jwa=1,mwa
      do 43 iax=1,nax
      fj(iax,1,jwa)=zero
      fj(iax,2,jwa)=zero
      fj(iax,3,jwa)=zero
   43 continue
   42 continue
c
      do 35 jwa=1,mwa
      qaj=qfaci*chg(iwq(jwa),1,ipset)
      paj=chg(iwq(jwa),2,ipset)
      do 36 iax=1,nax
      isa=isal(iax)
      qai=qfaci*chg(isq1(isa),1,ipset)
      pai=chg(isq1(isa),2,ipset)
      rx=xj(iax,1,jwa)-xi(iax,1)
      ry=xj(iax,2,jwa)-xi(iax,2)
      rz=xj(iax,3,jwa)-xi(iax,3)
      pix=pai*pl(iax,1)
      piy=pai*pl(iax,2)
      piz=pai*pl(iax,3)
      pjx=paj*pj(iax,1,jwa)
      pjy=paj*pj(iax,2,jwa)
      pjz=paj*pj(iax,3,jwa)
      ri2=one/(rx**2+ry**2+rz**2)
      ri3=sqrt(ri2)*ri2
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      fi(iax,1)=fi(iax,1)+((rmj-qaj)*rx-pjx)*ri3
      fi(iax,2)=fi(iax,2)+((rmj-qaj)*ry-pjy)*ri3
      fi(iax,3)=fi(iax,3)+((rmj-qaj)*rz-pjz)*ri3
      fj(iax,1,jwa)=fj(iax,1,jwa)+((rmi+qai)*rx-pix)*ri3
      fj(iax,2,jwa)=fj(iax,2,jwa)+((rmi+qai)*ry-piy)*ri3
      fj(iax,3,jwa)=fj(iax,3,jwa)+((rmi+qai)*rz-piz)*ri3
   36 continue
   35 continue
c
c     update the electric field arrays
c
      iax=0
      do 37 isa=lswndx(ispm-1,ipsw)+1,lswndx(ispm,ipsw)
      ispj=lswjpt(isa,ipsw)-1
      do 38 ismn=1,lswin(isa,ipsw)
      psp(isfr+isa,1,ipset-1,1)=psp(isfr+isa,1,ipset-1,1)+fi(iax+ismn,1)
      psp(isfr+isa,2,ipset-1,1)=psp(isfr+isa,2,ipset-1,1)+fi(iax+ismn,2)
      psp(isfr+isa,3,ipset-1,1)=psp(isfr+isa,3,ipset-1,1)+fi(iax+ismn,3)
   38 continue
      do 39 jwa=1,mwa
      do 40 ismn=1,lswin(isa,ipsw)
      lswptr=lswj(ispj+ismn)
      pwp(lswptr,1,jwa,ipset-1,1)=pwp(lswptr,1,jwa,ipset-1,1)+
     + fj(iax+ismn,1,jwa)
      pwp(lswptr,2,jwa,ipset-1,1)=pwp(lswptr,2,jwa,ipset-1,1)+
     + fj(iax+ismn,2,jwa)
      pwp(lswptr,3,jwa,ipset-1,1)=pwp(lswptr,3,jwa,ipset-1,1)+
     + fj(iax+ismn,3,jwa)
   40 continue
   39 continue
      iax=iax+lswin(isa,ipsw)
   37 continue
      endif
   30 continue
c
    3 continue
    1 continue
c
      return
      end
