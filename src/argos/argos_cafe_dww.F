      subroutine argos_cafe_dww(xw,xwm,pw,pwp,iwfrom,nlocw,lpbc,chg,iwq,
     + lwwndx,lwwjpt,lwwin,lwwj,rwc,xi,xj,rwx,pl,pj,fi,fj)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
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
      call argos_cafe_pbc(0,rwc,mscr,rwx,mscr,nax,1,lwwin(iwm,ipww))
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
