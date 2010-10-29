      subroutine argos_cafe_fpww(xw,xwm,fw,pw,pwp,idt,iwfrom,
     + nwloc,lpbc,eww,
     + vdw,chg,iwatm,iwq,lwwndx,lwwjpt,lwwin,lwwj,
     + xi,xj,rwx,rwi1,rwi2,rwi6,rwc,
     + f,fi,fj,facu,pl,pj)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
c
      real*8 xw(mwm,3,mwa),xwm(mwm,3),fw(mwm,3,mwa,2),eww(mpe,2)
      integer idt(mwm)
      integer iwfrom,nwloc
      logical lpbc
c
      real*8 vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      integer iwatm(mwa),iwq(mwa)
c
      real*8 xi(mscr,3,mwa),xj(mscr,3,mwa),rwx(mscr,3)
      real*8 rwi1(mscr),rwi2(mscr),rwi6(mscr),rwc(mscr,3)
      real*8 f(mscr),fi(mscr,3,mwa),fj(mscr,3,mwa)
c
      real*8 facu(mscr)
c     real*8 rdf(mgl,mgr)
c
      integer lwwj(*)
      integer lwwndx(0:mwm,2),lwwjpt(nwloc,2),lwwin(nwloc,2)
c
 
      real*8 pw(mwm,3,mwa,2),pwp(mwm,3,mwa,2,2)
      real*8 pl(mscr,3,mwa),pj(mscr,3,mwa)
c      integer nax2,ipset
      real*8 qai,qaj,pai,paj,pix,piy,piz,pjx,pjy,pjz
      real*8 ri3,rmi,rmj,fri,fmi,fmj,rmm,qfaci
      real*8 rx,ry,rz,ri1,ri2,ewwpsm,etermp
      real*8 ewwqsm
c
      integer iwfr,ipww,number,iwm,iwpm,nax
      integer iwmn,lwwptr,iwa,iax,jwa,iptr,jptr,iwpj
      real*8 ewwl6,ewwl12,q
      real*8 c64,c124,qi,qj,qi4,qj4,dercon
      real*8 c6p,c12p,qp,ep2tmp,ep3tmp
      real*8 c6,cf6,c12,cf12
c
      integer nwwlen(2)
      real*8 eterml,etermq
c
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
cx new stuff begin
c
cx new stuff end
c
c     calculation of solvent-solvent intermolecular energies and forces
c
c     subtract 1 from first molecule index for use as offset
c
      iwfr=iwfrom-1
c
c     loop over short and long range parts
c
      do 1 ipww=1,lpww
c
c     Evaluate the outer index array
c
      nwwlen(ipww)=0
      lwwndx(0,ipww)=0
      number=0
      do 2 iwm=1,nwloc
      if(number+lwwin(iwm,ipww).gt.mscr) then
      nwwlen(ipww)=nwwlen(ipww)+1
      lwwndx(nwwlen(ipww),ipww)=iwm-1
      number=0
      endif
      number=number+lwwin(iwm,ipww)
    2 continue
      if(number.gt.0) then
      nwwlen(ipww)=nwwlen(ipww)+1
      lwwndx(nwwlen(ipww),ipww)=nwloc
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
      rwc(nax+iwmn,1)=xwm(iwfr+iwm,1)-xwm(lwwptr,1)
      rwc(nax+iwmn,2)=xwm(iwfr+iwm,2)-xwm(lwwptr,2)
      rwc(nax+iwmn,3)=xwm(iwfr+iwm,3)-xwm(lwwptr,3)
      facu(nax+iwmn)=one
c      if( (iand(idt(iwm),mdynam).eq.ldynam.and.
c     + iand(idt(lwwptr),mdynam).ne.ldynam).or.
c     + (iand(idt(iwm),mdynam).ne.ldynam.and.
c     + iand(idt(lwwptr),mdynam).eq.ldynam) ) facu(nax+iwmn)=half
      if(iand(idt(iwm),mdynam).ne.ldynam.and.
     + iand(idt(lwwptr),mdynam).ne.ldynam) facu(nax+iwmn)=zero
      if(includ.eq.1) facu(nax+iwmn)=one
    5 continue
c
      do 6 iwa=1,mwa
      do 7 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      xi(nax+iwmn,1,iwa)=xw(iwfr+iwm,1,iwa)
      xi(nax+iwmn,2,iwa)=xw(iwfr+iwm,2,iwa)
      xi(nax+iwmn,3,iwa)=xw(iwfr+iwm,3,iwa)
      xj(nax+iwmn,1,iwa)=xw(lwwptr,1,iwa)
      xj(nax+iwmn,2,iwa)=xw(lwwptr,2,iwa)
      xj(nax+iwmn,3,iwa)=xw(lwwptr,3,iwa)
      pl(nax+iwmn,1,iwa)=pw(iwfr+iwm,1,iwa,1)
      pl(nax+iwmn,2,iwa)=pw(iwfr+iwm,2,iwa,1)
      pl(nax+iwmn,3,iwa)=pw(iwfr+iwm,3,iwa,1)
      pj(nax+iwmn,1,iwa)=pw(lwwptr,1,iwa,1)
      pj(nax+iwmn,2,iwa)=pw(lwwptr,2,iwa,1)
      pj(nax+iwmn,3,iwa)=pw(lwwptr,3,iwa,1)
    7 continue
    6 continue
      if(lpbc) then
      call argos_cafe_pbc(0,rwc,mscr,rwx,mscr,nax,1,lwwin(iwm,ipww))
      do 8 iwmn=1,lwwin(iwm,ipww)
      rwc(nax+iwmn,1)=rwc(nax+iwmn,1)-rwx(iwmn,1)
      rwc(nax+iwmn,2)=rwc(nax+iwmn,2)-rwx(iwmn,2)
      rwc(nax+iwmn,3)=rwc(nax+iwmn,3)-rwx(iwmn,3)
    8 continue
      do 9 iwa=1,mwa
      do 10 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      xj(nax+iwmn,1,iwa)=xj(nax+iwmn,1,iwa)+rwx(iwmn,1)
      xj(nax+iwmn,2,iwa)=xj(nax+iwmn,2,iwa)+rwx(iwmn,2)
      xj(nax+iwmn,3,iwa)=xj(nax+iwmn,3,iwa)+rwx(iwmn,3)
   10 continue
    9 continue
      endif
c
      nax=nax+lwwin(iwm,ipww)
    4 continue
c
c     initializations
c
c      if(npener.ne.0) then
c      do 12 iax=1,nax
c      u(iax)=zero
c   12 continue
c      endif
c
c     loops over number of atoms in a solvent molecule
c
      qfaci=one/qfac
      do 13 iwa=1,mwa
      qi=chg(iwq(iwa),1,iset)
      pai=chg(iwq(iwa),2,iset)
      qai=qfaci*qi
      do 14 jwa=1,mwa
      qj=chg(iwq(jwa),1,iset)
      q=qi*qj
      paj=chg(iwq(jwa),2,iset)
      qaj=qfaci*qj
c
      do 15 iax=1,nax
      f(iax)=zero
      rwx(iax,1)=xi(iax,1,iwa)-xj(iax,1,jwa)
      rwx(iax,2)=xi(iax,2,iwa)-xj(iax,2,jwa)
      rwx(iax,3)=xi(iax,3,iwa)-xj(iax,3,jwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
   15 continue
c
c
c     van der Waals contribution
c     --------------------------
c
      iptr=iwatm(iwa)
      jptr=iwatm(jwa)
      c6=vdw(iptr,jptr,1,iset)
      cf6=six*c6
      c12=vdw(iptr,jptr,3,iset)
      cf12=twelve*c12
c
      eterml=zero
      if(c6.ne.zero.or.c12.ne.zero) then
      ewwl6=zero
      ewwl12=zero
      do 20 iax=1,nax
      rwi6(iax)=rwi2(iax)*rwi2(iax)*rwi2(iax)
      ewwl6=ewwl6+facu(iax)*rwi6(iax)
      ewwl12=ewwl12+facu(iax)*rwi6(iax)*rwi6(iax)
      f(iax)=f(iax)+(cf12*rwi6(iax)-cf6)*rwi6(iax)*rwi2(iax)
   20 continue
      eterml=c12*ewwl12-c6*ewwl6
      eww(7,ipww)=eww(7,ipww)+eterml
      endif
c
c
c     force vectors
c     -------------
c
      if(iwa.eq.1) then
      do 22 iax=1,nax
      fj(iax,1,jwa)=(-f(iax))*rwx(iax,1)
      fj(iax,2,jwa)=(-f(iax))*rwx(iax,2)
      fj(iax,3,jwa)=(-f(iax))*rwx(iax,3)
   22 continue
      else
      do 23 iax=1,nax
      fj(iax,1,jwa)=fj(iax,1,jwa)-f(iax)*rwx(iax,1)
      fj(iax,2,jwa)=fj(iax,2,jwa)-f(iax)*rwx(iax,2)
      fj(iax,3,jwa)=fj(iax,3,jwa)-f(iax)*rwx(iax,3)
   23 continue
      endif
c
      if(jwa.eq.1) then
      do 24 iax=1,nax
      fi(iax,1,iwa)=f(iax)*rwx(iax,1)
      fi(iax,2,iwa)=f(iax)*rwx(iax,2)
      fi(iax,3,iwa)=f(iax)*rwx(iax,3)
   24 continue
      else
      do 25 iax=1,nax
      fi(iax,1,iwa)=fi(iax,1,iwa)+f(iax)*rwx(iax,1)
      fi(iax,2,iwa)=fi(iax,2,iwa)+f(iax)*rwx(iax,2)
      fi(iax,3,iwa)=fi(iax,3,iwa)+f(iax)*rwx(iax,3)
   25 continue
      endif
      do 26 iax=1,nax
      zw(1,1,ipww)=zw(1,1,ipww)-f(iax)*rwx(iax,1)*rwc(iax,1)
      zw(2,1,ipww)=zw(2,1,ipww)-f(iax)*rwx(iax,1)*rwc(iax,2)
      zw(3,1,ipww)=zw(3,1,ipww)-f(iax)*rwx(iax,1)*rwc(iax,3)
      zw(1,2,ipww)=zw(1,2,ipww)-f(iax)*rwx(iax,2)*rwc(iax,1)
      zw(2,2,ipww)=zw(2,2,ipww)-f(iax)*rwx(iax,2)*rwc(iax,2)
      zw(3,2,ipww)=zw(3,2,ipww)-f(iax)*rwx(iax,2)*rwc(iax,3)
      zw(1,3,ipww)=zw(1,3,ipww)-f(iax)*rwx(iax,3)*rwc(iax,1)
      zw(2,3,ipww)=zw(2,3,ipww)-f(iax)*rwx(iax,3)*rwc(iax,2)
      zw(3,3,ipww)=zw(3,3,ipww)-f(iax)*rwx(iax,3)*rwc(iax,3)
   26 continue
c
c
c     electrostatic and polarization contribution
c     -------------------------------------------
c
      ewwqsm=zero
      ewwpsm=zero
      do 117 iax=1,nax
      pix=pai*pl(iax,1,iwa)
      piy=pai*pl(iax,2,iwa)
      piz=pai*pl(iax,3,iwa)
      pjx=paj*pj(iax,1,jwa)
      pjy=paj*pj(iax,2,jwa)
      pjz=paj*pj(iax,3,jwa)
      rx=-rwx(iax,1)
      ry=-rwx(iax,2)
      rz=-rwx(iax,3)
      ri1=rwi1(iax)
      ri2=rwi2(iax)
      ri3=qfac*qfac*ri1*ri2
      rmi=three*(rx*pix+ry*piy+rz*piz)*ri2
      rmj=three*(rx*pjx+ry*pjy+rz*pjz)*ri2
      if(ipolt.eq.1) then
      fri=((-qai)*qaj+qai*rmj-qaj*rmi)*ri3
      fmi=(qaj)*ri3
      fmj=(-qai)*ri3
      else
      rmm=three*(pix*pjx+piy*pjy+piz*pjz)*ri2
      fri=((-qai)*qaj+qai*rmj-qaj*rmi+5.0*rmi*rmj/three-rmm)*ri3
      fmi=(qaj-rmj)*ri3
      fmj=((-qai)-rmi)*ri3
      endif
      fi(iax,1,iwa)=fi(iax,1,iwa)+fri*rx+fmi*pix+fmj*pjx
      fi(iax,2,iwa)=fi(iax,2,iwa)+fri*ry+fmi*piy+fmj*pjy
      fi(iax,3,iwa)=fi(iax,3,iwa)+fri*rz+fmi*piz+fmj*pjz
      fj(iax,1,jwa)=fj(iax,1,jwa)-(fri*rx+fmi*pix+fmj*pjx)
      fj(iax,2,jwa)=fj(iax,2,jwa)-(fri*ry+fmi*piy+fmj*pjy)
      fj(iax,3,jwa)=fj(iax,3,jwa)-(fri*rz+fmi*piz+fmj*pjz)
      zw(1,1,ipww)=zw(1,1,ipww)-(fri*rx+fmi*pix+fmj*pjx)*rwc(iax,1)
      zw(2,1,ipww)=zw(2,1,ipww)-(fri*rx+fmi*pix+fmj*pjx)*rwc(iax,2)
      zw(3,1,ipww)=zw(3,1,ipww)-(fri*rx+fmi*pix+fmj*pjx)*rwc(iax,3)
      zw(1,2,ipww)=zw(1,2,ipww)-(fri*ry+fmi*piy+fmj*pjy)*rwc(iax,1)
      zw(2,2,ipww)=zw(2,2,ipww)-(fri*ry+fmi*piy+fmj*pjy)*rwc(iax,2)
      zw(3,2,ipww)=zw(3,2,ipww)-(fri*ry+fmi*piy+fmj*pjy)*rwc(iax,3)
      zw(1,3,ipww)=zw(1,3,ipww)-(fri*rz+fmi*piz+fmj*pjz)*rwc(iax,1)
      zw(2,3,ipww)=zw(2,3,ipww)-(fri*rz+fmi*piz+fmj*pjz)*rwc(iax,2)
      zw(3,3,ipww)=zw(3,3,ipww)-(fri*rz+fmi*piz+fmj*pjz)*rwc(iax,3)
      ewwpsm=ewwpsm+facu(iax)*(qai*rmj-qaj*rmi)*ri1
      ewwqsm=ewwqsm+facu(iax)*ri1
  117 continue
      etermp=-qfac*qfac*ewwpsm/three
      eww(8,ipww)=eww(8,ipww)+etermp
      etermq=q*ewwqsm
      eww(8,ipww)=eww(8,ipww)+etermq
c
c     Radial distribution functions
c
c      if(ifstep-1.eq.((ifstep-1)/nfrdf)*nfrdf .and. ngrww.gt.0) then
c      do 27 igc=1,ngc
c      if(ngt(igc).eq.1) then
c      if(iagc(igc).eq.iwa .and. jagc(igc).eq.jwa) then
c      igr=igrc(igc)
c      do 28 iax=1,nax
c      indx=int(one/(rwi1(iax)*drdf))
c      if(indx.le.ngl) rdf(indx,igr)=rdf(indx,igr)+rdfvol
c   28 continue
c      endif
c      endif
c   27 continue
c      endif
c
c     Thermodynamic integration
c
      if(ithint) then
      if(ith(2)) then
      c64=vdw(iwatm(iwa),iwatm(jwa),1,4)
      c124=vdw(iwatm(iwa),iwatm(jwa),3,4)
      ewwl6=zero
      ewwl12=zero
      do 29 iax=1,nax
      ewwl6=ewwl6+facu(iax)*rwi6(iax)
      ewwl12=ewwl12+facu(iax)*rwi6(iax)*rwi6(iax)
   29 continue
      deriv(2,ipww)=deriv(2,ipww)+c124*ewwl12-c64*ewwl6
      endif
      if(ith(4)) then
      qi=chg(iwq(iwa),1,iset)
      qj=chg(iwq(jwa),1,iset)
      qi4=chg(iwq(iwa),1,4)
      qj4=chg(iwq(jwa),1,4)
      dercon=zero
      if(ipme.eq.0) then
      do 30 iax=1,nax
      dercon=dercon+rwi1(iax)
   30 continue
      else
      do 130 iax=1,nax
      dercon=dercon+rwi1(iax)
  130 continue
      endif
      deriv(4,ipww)=deriv(4,ipww)+(qi*qj4+qj*qi4)*dercon
      if(ireact.ne.0) then
      dercon=zero
      do 31 iax=1,nax
      dercon=dercon+one/rwi2(iax)
   31 continue
      deriv(4,ipww)=deriv(4,ipww)+(qi*qj4+qj*qi4)*rffww*dercon
      endif
      endif
      endif
c
c     Thermodynamic perturbation 1
c
      if(ipert2) then
      if(ip2(2)) then
      c6p=vdw(iwatm(iwa),iwatm(jwa),1,2)
      c12p=vdw(iwatm(iwa),iwatm(jwa),3,2)
      do 32 iax=1,nax
      ep2(ipww)=ep2(ipww)+facu(iax)*(c12p*rwi6(iax)-c6p)*rwi6(iax)
   32 continue
      ep2(ipww)=ep2(ipww)-eterml
      endif
      if(ip2(4).or.ip2(5)) then
      qp=chg(iwq(iwa),1,2)*chg(iwq(jwa),1,2)
      ep2tmp=zero
      do 33 iax=1,nax
      rwx(iax,1)=xi(iax,1,iwa)-xj(iax,1,jwa)
      rwx(iax,2)=xi(iax,2,iwa)-xj(iax,2,jwa)
      rwx(iax,3)=xi(iax,3,iwa)-xj(iax,3,jwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      if(ipme.eq.0) then
      ep2tmp=ep2tmp+facu(iax)*rwi1(iax)
      else
      ep2tmp=ep2tmp+facu(iax)*erfc(ealpha/rwi1(iax))*rwi1(iax)
      endif
   33 continue
      ep2(ipww)=ep2(ipww)+qp*ep2tmp-etermq
      if(ireact.ne.0) then
      ep2tmp=zero
      do 34 iax=1,nax
      ep2tmp=ep2tmp+facu(iax)/rwi2(iax)
   34 continue
      ep2(ipww)=ep2(ipww)+qp*rffww*ep2tmp
      endif
      endif
      endif
c
c     Thermodynamic perturbation 2
c
      if(ipert3) then
      if(ip3(2)) then
      c6p=vdw(iwatm(iwa),iwatm(jwa),1,3)
      c12p=vdw(iwatm(iwa),iwatm(jwa),3,3)
      do 35 iax=1,nax
      ep3(ipww)=ep3(ipww)+facu(iax)*(c12p*rwi6(iax)-c6p)*rwi6(iax)
   35 continue
      ep3(ipww)=ep3(ipww)-eterml
      endif
      if(ip2(4).or.ip2(5)) then
      qp=chg(iwatm(iwa),1,3)*chg(iwatm(jwa),1,3)
      ep3tmp=zero
      do 36 iax=1,nax
      rwx(iax,1)=xi(iax,1,iwa)-xj(iax,1,jwa)
      rwx(iax,2)=xi(iax,2,iwa)-xj(iax,2,jwa)
      rwx(iax,3)=xi(iax,3,iwa)-xj(iax,3,jwa)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      if(ipme.eq.0) then
      ep3tmp=ep3tmp+facu(iax)*rwi1(iax)
      else
      ep3tmp=ep3tmp+facu(iax)*erfc(ealpha/rwi1(iax))*rwi1(iax)
      endif
   36 continue
      ep3(ipww)=ep3(ipww)+qp*ep3tmp-etermq
      if(ireact.ne.0) then
      ep3tmp=zero
      do 37 iax=1,nax
      ep3tmp=ep3tmp+facu(iax)/rwi2(iax)
   37 continue
      ep3(ipww)=ep3(ipww)+qp*rffww*ep3tmp
      endif
      endif
      endif
   14 continue
   13 continue
c
c     Update force arrays
c
      iax=0
      do 38 iwm=lwwndx(iwpm-1,ipww)+1,lwwndx(iwpm,ipww)
      iwpj=lwwjpt(iwm,ipww)-1
      do 39 iwa=1,mwa
      do 40 iwmn=1,lwwin(iwm,ipww)
      lwwptr=lwwj(iwpj+iwmn)
      fw(iwfr+iwm,1,iwa,ipww)=fw(iwfr+iwm,1,iwa,ipww)+fi(iax+iwmn,1,iwa)
      fw(iwfr+iwm,2,iwa,ipww)=fw(iwfr+iwm,2,iwa,ipww)+fi(iax+iwmn,2,iwa)
      fw(iwfr+iwm,3,iwa,ipww)=fw(iwfr+iwm,3,iwa,ipww)+fi(iax+iwmn,3,iwa)
      fw(lwwptr,1,iwa,ipww)=fw(lwwptr,1,iwa,ipww)+fj(iax+iwmn,1,iwa)
      fw(lwwptr,2,iwa,ipww)=fw(lwwptr,2,iwa,ipww)+fj(iax+iwmn,2,iwa)
      fw(lwwptr,3,iwa,ipww)=fw(lwwptr,3,iwa,ipww)+fj(iax+iwmn,3,iwa)
   40 continue
   39 continue
c
c     update energy arrays if appropriate print option was set
c
c      if(npener.ne.0) then
c      do 41 iwmn=1,lwwin(iwm,ipww)
c      lwwptr=lwwj(iwpj+iwmn)
c      uwmw(iwfr+iwm)=uwmw(iwfr+iwm)+u(iax+iwmn)
c      uwmw(lwwptr)=uwmw(lwwptr)+u(iax+iwmn)
c   41 continue
c      endif
c
      iax=iax+lwwin(iwm,ipww)
   38 continue
    3 continue
c
    1 continue
c
      return
      end
