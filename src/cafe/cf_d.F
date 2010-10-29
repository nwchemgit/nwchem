      subroutine cf_dww(xw,xwm,pw,pwp,iwfrom,nlocw,lpbc,chg,iwq,
     + lwwndx,lwwjpt,lwwin,lwwj,rwc,xi,xj,rwx,pl,pj,fi,fj)
c
c $Id$
c
      implicit none
c
#include "cf_common.fh"
c
      real*8 xw(mwm,3,mwa),xwm(mwm,3)
      real*8 pw(mwm,3,mwa,2),pwp(mwm,3,mwa,2,2)
      integer iwfrom,nlocw
      logical lpbc
      real*8 chg(mqt,mqp,mset)
      integer iwq(mwa)
      integer lwwndx(0:mwm,2),lwwin(nlocw,2),lwwjpt(nlocw,2),lwwj(*)
c
      real*8 rwc(mscr,3),rwx(mscr,3)
      real*8 xi(mscr,3,mwa),xj(mscr,3,mwa)
      real*8 pl(mscr,3,mwa),pj(mscr,3,mwa)
      real*8 fi(mscr,3,mwa),fj(mscr,3,mwa)
c
      integer ix,ipset,nax2,nwwlen(2)
      integer iwfr,ipww,number,iwm,iwpm,nax,iwpj,iwa,jwa,iax,iwmn,lwwptr
      real*8 pai,paj,qai,qaj
      real*8 rx,ry,rz,ri2,ri3,rmi,rmj,pix,piy,piz,pjx,pjy,pjz
c
      real*8 qfaci
c
      qfaci=one/qfac
c
c     calculation of solvent-solvent intermolecular energies and forces
c
      iwfr=iwfrom-1
c
c     loop over both short and long range parts
c
      do 1 ipww=1,npww
c
c     Evaluate the outer index array
c
      nwwlen(ipww)=0
      lwwndx(0,ipww)=0
      number=0
      do 2 iwm=1,nlocw
      if(number+lwwin(iwm,ipww).gt.mscr) then
      nwwlen(ipww)=nwwlen(ipww)+1
      lwwndx(nwwlen(ipww),ipww)=iwm-1
      number=0
      endif
      number=number+lwwin(iwm,ipww)
    2 continue
      if(number.gt.0) then
      nwwlen(ipww)=nwwlen(ipww)+1
      lwwndx(nwwlen(ipww),ipww)=nlocw
      endif
c
c     loop over number of cycles to complete pairlist
c
      do 3 iwpm=1,nwwlen(ipww)
      nax=0
c
c     collect coordinates into workarrays
c
      do 4 iwm=lwwndx(iwpm-1,ipww)+1,lwwndx(iwpm,ipww)
      iwpj=lwwjpt(iwm,ipww)-1
      do 5 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      rwc(nax+iwmn,1)=xwm(lwwptr,1)-xwm(iwfr+iwm,1)
      rwc(nax+iwmn,2)=xwm(lwwptr,2)-xwm(iwfr+iwm,2)
      rwc(nax+iwmn,3)=xwm(lwwptr,3)-xwm(iwfr+iwm,3)
    5 continue
c
      if(.not.lpbc) then
      do 6 iwa=1,mwa
      do 7 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      xi(nax+iwmn,1,iwa)=xw(iwfr+iwm,1,iwa)
      xi(nax+iwmn,2,iwa)=xw(iwfr+iwm,2,iwa)
      xi(nax+iwmn,3,iwa)=xw(iwfr+iwm,3,iwa)
      xj(nax+iwmn,1,iwa)=xw(lwwptr,1,iwa)
      xj(nax+iwmn,2,iwa)=xw(lwwptr,2,iwa)
      xj(nax+iwmn,3,iwa)=xw(lwwptr,3,iwa)
      pl(nax+iwmn,1,iwa)=pw(iwfr+iwm,1,iwa,2)
      pl(nax+iwmn,2,iwa)=pw(iwfr+iwm,2,iwa,2)
      pl(nax+iwmn,3,iwa)=pw(iwfr+iwm,3,iwa,2)
      pj(nax+iwmn,1,iwa)=pw(lwwptr,1,iwa,2)
      pj(nax+iwmn,2,iwa)=pw(lwwptr,2,iwa,2)
      pj(nax+iwmn,3,iwa)=pw(lwwptr,3,iwa,2)
    7 continue
    6 continue
      else
      call cf_pbc(0,rwc,mscr,rwx,mscr,nax,1,lwwin(iwm,ipww))
      do 9 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      rwc(nax+iwmn,1)=rwc(nax+iwmn,1)-rwx(iwmn,1)
      rwc(nax+iwmn,2)=rwc(nax+iwmn,2)-rwx(iwmn,2)
      rwc(nax+iwmn,3)=rwc(nax+iwmn,3)-rwx(iwmn,3)
    9 continue
c
      do 10 iwa=1,mwa
      do 11 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      xi(nax+iwmn,1,iwa)=xw(iwfr+iwm,1,iwa)
      xi(nax+iwmn,2,iwa)=xw(iwfr+iwm,2,iwa)
      xi(nax+iwmn,3,iwa)=xw(iwfr+iwm,3,iwa)
      xj(nax+iwmn,1,iwa)=xw(lwwptr,1,iwa)-rwx(iwmn,1)
      xj(nax+iwmn,2,iwa)=xw(lwwptr,2,iwa)-rwx(iwmn,2)
      xj(nax+iwmn,3,iwa)=xw(lwwptr,3,iwa)-rwx(iwmn,3)
      pl(nax+iwmn,1,iwa)=pw(iwfr+iwm,1,iwa,2)
      pl(nax+iwmn,2,iwa)=pw(iwfr+iwm,2,iwa,2)
      pl(nax+iwmn,3,iwa)=pw(iwfr+iwm,3,iwa,2)
      pj(nax+iwmn,1,iwa)=pw(lwwptr,1,iwa,2)
      pj(nax+iwmn,2,iwa)=pw(lwwptr,2,iwa,2)
      pj(nax+iwmn,3,iwa)=pw(lwwptr,3,iwa,2)
   11 continue
   10 continue
      endif
      nax=nax+lwwin(iwm,ipww)
    4 continue
c
c     zero temporary arrays fi and fj
c
      do 12 iwa=1,mwa
      do 13 ix=1,3
      do 14 iax=1,nax
      fi(iax,ix,iwa)=zero
      fj(iax,ix,iwa)=zero
   14 continue
   13 continue
   12 continue
c
c     loops over number of atoms in a solvent molecule
c
c     calculated here is 4*pi*epsilon*field and not just the field
c     since the polarization is given in alpha/(4*pi*epsilon) in
c     stead of just alpha, the induced dipole is obtained by the
c     product of pwa and pw
c
      do 20 iwa=1,mwa
      qai=qfaci*chg(iwq(iwa),1,iset)
      pai=chg(iwq(iwa),2,iset)
      do 21 jwa=1,mwa
      qaj=qfaci*chg(iwq(jwa),1,iset)
      paj=chg(iwq(jwa),2,iset)
      do 22 iax=1,nax
      rx=xj(iax,1,jwa)-xi(iax,1,iwa)
      ry=xj(iax,2,jwa)-xi(iax,2,iwa)
      rz=xj(iax,3,jwa)-xi(iax,3,iwa)
      pix=pai*pl(iax,1,iwa)
      piy=pai*pl(iax,2,iwa)
      piz=pai*pl(iax,3,iwa)
      pjx=paj*pj(iax,1,jwa)
      pjy=paj*pj(iax,2,jwa)
      pjz=paj*pj(iax,3,jwa)
      ri2=one/(rx**2+ry**2+rz**2)
      ri3=sqrt(ri2)*ri2
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      fi(iax,1,iwa)=fi(iax,1,iwa)+((rmj-qaj)*rx-pjx)*ri3
      fi(iax,2,iwa)=fi(iax,2,iwa)+((rmj-qaj)*ry-pjy)*ri3
      fi(iax,3,iwa)=fi(iax,3,iwa)+((rmj-qaj)*rz-pjz)*ri3
      fj(iax,1,jwa)=fj(iax,1,jwa)+((rmi+qai)*rx-pix)*ri3
      fj(iax,2,jwa)=fj(iax,2,jwa)+((rmi+qai)*ry-piy)*ri3
      fj(iax,3,jwa)=fj(iax,3,jwa)+((rmi+qai)*rz-piz)*ri3
   22 continue
   21 continue
   20 continue
c
c     Update the electric field arrays
c
      iax=0
      do 23 iwm=lwwndx(iwpm-1,ipww)+1,lwwndx(iwpm,ipww)
      iwpj=lwwjpt(iwm,ipww)-1
      do 24 iwa=1,mwa
      do 25 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      pw(iwfr+iwm,1,iwa,1)=pw(iwfr+iwm,1,iwa,1)+fi(iax+iwmn,1,iwa)
      pw(iwfr+iwm,2,iwa,1)=pw(iwfr+iwm,2,iwa,1)+fi(iax+iwmn,2,iwa)
      pw(iwfr+iwm,3,iwa,1)=pw(iwfr+iwm,3,iwa,1)+fi(iax+iwmn,3,iwa)
      pw(lwwptr,1,iwa,1)=pw(lwwptr,1,iwa,1)+fj(iax+iwmn,1,iwa)
      pw(lwwptr,2,iwa,1)=pw(lwwptr,2,iwa,1)+fj(iax+iwmn,2,iwa)
      pw(lwwptr,3,iwa,1)=pw(lwwptr,3,iwa,1)+fj(iax+iwmn,3,iwa)
   25 continue
   24 continue
      iax=iax+lwwin(iwm,ipww)
   23 continue
c
c     thermodynamic perturbation and integration
c
      do 30 ipset=2,3
      if((ipset.eq.2.and.ipert2).or.
     + (ipset.eq.3.and.ipert3)) then
c
      nax2=0
      do 31 iwm=lwwndx(iwpm-1,ipww)+1,lwwndx(iwpm,ipww)
      iwpj=lwwjpt(iwm,ipww)-1
      do 32 iwa=1,mwa
      do 33 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      pl(nax2+iwmn,1,iwa)=pwp(iwfr+iwm,1,iwa,ipset-1,2)
      pl(nax2+iwmn,2,iwa)=pwp(iwfr+iwm,2,iwa,ipset-1,2)
      pl(nax2+iwmn,3,iwa)=pwp(iwfr+iwm,3,iwa,ipset-1,2)
      pj(nax2+iwmn,1,iwa)=pwp(lwwptr,1,iwa,ipset-1,2)
      pj(nax2+iwmn,2,iwa)=pwp(lwwptr,2,iwa,ipset-1,2)
      pj(nax2+iwmn,3,iwa)=pwp(lwwptr,3,iwa,ipset-1,2)
   33 continue
   32 continue
      nax2=nax2+lwwin(iwm,ipww)
   31 continue
c
      if(nax.ne.nax2) call md_abort('Error in dipww',me)
c
c
      do 40 iwa=1,mwa
      do 41 ix=1,3
      do 42 iax=1,nax
      fi(iax,ix,iwa)=zero
      fj(iax,ix,iwa)=zero
   42 continue
   41 continue
   40 continue
c
      do 34 iwa=1,mwa
      qai=qfaci*chg(iwq(iwa),1,ipset)
      pai=chg(iwq(iwa),2,ipset)
      do 35 jwa=1,mwa
      qaj=qfaci*chg(iwq(jwa),1,ipset)
      paj=chg(iwq(jwa),2,ipset)
      do 36 iax=1,nax
      rx=xj(iax,1,jwa)-xi(iax,1,iwa)
      ry=xj(iax,2,jwa)-xi(iax,2,iwa)
      rz=xj(iax,3,jwa)-xi(iax,3,iwa)
      pix=pai*pl(iax,1,iwa)
      piy=pai*pl(iax,2,iwa)
      piz=pai*pl(iax,3,iwa)
      pjx=paj*pj(iax,1,jwa)
      pjy=paj*pj(iax,2,jwa)
      pjz=paj*pj(iax,3,jwa)
      ri2=one/(rx**2+ry**2+rz**2)
      ri3=sqrt(ri2)*ri2
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      fi(iax,1,iwa)=fi(iax,1,iwa)+((rmj-qaj)*rx-pjx)*ri3
      fi(iax,2,iwa)=fi(iax,2,iwa)+((rmj-qaj)*ry-pjy)*ri3
      fi(iax,3,iwa)=fi(iax,3,iwa)+((rmj-qaj)*rz-pjz)*ri3
      fj(iax,1,jwa)=fj(iax,1,jwa)+((rmi+qai)*rx-pix)*ri3
      fj(iax,2,jwa)=fj(iax,2,jwa)+((rmi+qai)*ry-piy)*ri3
      fj(iax,3,jwa)=fj(iax,3,jwa)+((rmi+qai)*rz-piz)*ri3
   36 continue
   35 continue
   34 continue
c
c     Update the electric field arrays
c
      iax=0
      do 37 iwm=lwwndx(iwpm-1,ipww)+1,lwwndx(iwpm,ipww)
      iwpj=lwwjpt(iwm,ipww)-1
      do 38 iwa=1,mwa
      do 39 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      pwp(iwfr+iwm,1,iwa,ipset-1,1)=pwp(iwfr+iwm,1,iwa,ipset-1,1)+
     + fi(iax+iwmn,1,iwa)
      pwp(iwfr+iwm,2,iwa,ipset-1,1)=pwp(iwfr+iwm,2,iwa,ipset-1,1)+
     + fi(iax+iwmn,2,iwa)
      pwp(iwfr+iwm,3,iwa,ipset-1,1)=pwp(iwfr+iwm,3,iwa,ipset-1,1)+
     + fi(iax+iwmn,3,iwa)
      pwp(lwwptr,1,iwa,ipset-1,1)=pwp(lwwptr,1,iwa,ipset-1,1)+
     + fj(iax+iwmn,1,iwa)
      pwp(lwwptr,2,iwa,ipset-1,1)=pwp(lwwptr,2,iwa,ipset-1,1)+
     + fj(iax+iwmn,2,iwa)
      pwp(lwwptr,3,iwa,ipset-1,1)=pwp(lwwptr,3,iwa,ipset-1,1)+
     + fj(iax+iwmn,3,iwa)
   39 continue
   38 continue
      iax=iax+lwwin(iwm,ipww)
   37 continue
      endif
   30 continue
c
c
    3 continue
    1 continue
c
      return
      end
      subroutine cf_dsw(xs,xsm,ps,psp,isdt,ismf,isml,isq1,isfrom,nums,
     + xw,xwm,pw,pwp,iwq,lpbc,chg,lswndx,lswjpt,lswin,lswj,
     + rwc,xi,xj,rwx,pl,pj,fi,fj,isal)
c
c $Id$
c
      implicit none
c
#include "cf_common.fh"
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
      if(lpbc) call cf_pbc(0,rwc,mscr,rwx,mscr,nax,1,lswin(isa,ipsw))
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
      subroutine cf_dss(xs,xsm,ps,psp,ismf,isml,isq2,isq3,isfrom,nums,
     + lpbc,chg,lssndx,lssjpt,lssin,lssj,xi,xj,rwc,rwi1,rwi2,rwi6,
     + rwx,rw,fi,fj,f,isal,jsal,jmal,jfal,pl,pj)
c
c $Id$
c
      implicit none
c
#include "cf_common.fh"
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
      call cf_pbc(0,rwc,mscr,rwx,mscr,nax,1,lssin(isa,ipss))
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
