      subroutine argos_cafe_fpsw(xs,xsm,fs,zs,ps,psp,
     + isga,isat,isdt,ismf,isml,isss,isq1,
     + isfrom,nums,xw,xwm,fw,pw,pwp,rtos,iwdt,lpbc,lpbcs,esw,esa,
     + vdw,chg,iwatm,iwq,iass,lswndx,lswjpt,lswin,lswj,
     + xi,xj,rwx,rwi1,rwi2,rwi6,rwc,f,fi,fj,facu,
     + rw,isal,isrx,list,pl,pj)
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
      real*8 ps(msa,3,2),psp(msa,3,2,2)
      real*8 pw(mwm,3,mwa,2),pwp(mwm,3,mwa,2,2)
      real*8 pl(mscr,3),pj(mscr,3,mwa)
c
      integer isatm,nswlen(2)
      integer isfr,iwm,ipsw,number,isa,ispm,isf,nax,ism
      integer ispj,ismn,lswptr,iwa,iax,iwatmi,ix,iy
      integer iwatyp
      real*8 c6,cf6,c12,cf12,sumen
      real*8 c64,c124,dercon,qj,qj4,derco1,derco2
      real*8 drvco1,drvco2,derco3,drvco3,c6p,c12p,etermq,eterml
      real*8 qi,qai,qaj,pai,paj,pix,piy,piz,pjx,pjy,pjz
      real*8 rx,ry,rz,ri1,ri2,ri3,rmi,rmj,fri,fmi,fmj,rmm
      real*8 zxx,zxy,zxz,zyx,zyy,zyz,zzx,zzy,zzz
      real*8 eswqsm,eswpsm,qfaci
      real*8 rtmp
c
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
      etermq=zero
c
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
      qfaci=one/qfac
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
      pl(nax+ismn,1)=ps(isfr+isa,1,1)
      pl(nax+ismn,2)=ps(isfr+isa,2,1)
      pl(nax+ismn,3)=ps(isfr+isa,3,1)
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
      pl(nax+ismn,1)=ps(isfr+isa,1,1)
      pl(nax+ismn,2)=ps(isfr+isa,2,1)
      pl(nax+ismn,3)=ps(isfr+isa,3,1)
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
      pj(nax+ismn,1,iwa)=pw(lswptr,1,iwa,1)
      pj(nax+ismn,2,iwa)=pw(lswptr,2,iwa,1)
      pj(nax+ismn,3,iwa)=pw(lswptr,3,iwa,1)
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
      pj(nax+ismn,1,iwa)=pw(lswptr,1,iwa,1)
      pj(nax+ismn,2,iwa)=pw(lswptr,2,iwa,1)
      pj(nax+ismn,3,iwa)=pw(lswptr,3,iwa,1)
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
      qj=chg(iwq(iwa),1,iset)
      qaj=qfaci*qj
c      dqj=qwa(iwa,4)
c      dqaj=qfaci*dqj
      paj=chg(iwq(iwa),2,iset)
c      dpaj=pwa(iwa,4)
      eswqsm=zero
      eswpsm=zero
c      dswqsm=zero
c      dswqws=zero
c      dswqps=zero
c      dswpss=zero
c      dswpws=zero
      do 21 iax=1,nax
      isa=isal(iax)
      qi=chg(isq1(isa),1,iset)
c      dqi=qsa(isa,4,1)
      qai=qfaci*qi
c      dqai=qfaci*dqi
      pai=chg(isq1(isa),2,iset)
c      dpai=psa(isa,4)
      pix=pai*pl(iax,1)
      piy=pai*pl(iax,2)
      piz=pai*pl(iax,3)
      pjx=paj*pj(iax,1,iwa)
      pjy=paj*pj(iax,2,iwa)
      pjz=paj*pj(iax,3,iwa)
      rwx(iax,1)=xi(iax,1)-xj(iax,1,iwa)
      rwx(iax,2)=xi(iax,2)-xj(iax,2,iwa)
      rwx(iax,3)=xi(iax,3)-xj(iax,3,iwa)
      rx=-rwx(iax,1)
      ry=-rwx(iax,2)
      rz=-rwx(iax,3)
      rwi2(iax)=one/(rx**2+ry**2+rz**2)
      rwi1(iax)=sqrt(rwi2(iax))
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
      fi(iax,1,1)=fi(iax,1,1)+fri*rx+fmi*pix+fmj*pjx
      fi(iax,2,1)=fi(iax,2,1)+fri*ry+fmi*piy+fmj*pjy
      fi(iax,3,1)=fi(iax,3,1)+fri*rz+fmi*piz+fmj*pjz
      fj(iax,1,iwa)=fj(iax,1,iwa)-(fri*rx+fmi*pix+fmj*pjx)
      fj(iax,2,iwa)=fj(iax,2,iwa)-(fri*ry+fmi*piy+fmj*pjy)
      fj(iax,3,iwa)=fj(iax,3,iwa)-(fri*rz+fmi*piz+fmj*pjz)
      eswqsm=eswqsm+qi*facu(iax)*ri1
      eswpsm=eswpsm+facu(iax)*(qai*rmj-qaj*rmi)*ri1
c      if(ithint.ne.0) then
c      dpix=dpai*pl(iax,1)
c      dpiy=dpai*pl(iax,2)
c      dpiz=dpai*pl(iax,3)
c      dpjx=dpaj*pj(iax,1,jwa)
c      dpjy=dpaj*pj(iax,2,jwa)
c      dpjz=dpaj*pj(iax,3,jwa)
c      drmi=three*(rx*dpix+ry*dpiy+rz*dpiz)*ri2
c      drmj=three*(rx*dpjx+ry*dpjy+rz*dpjz)*ri2
c      dswqsm=dswqsm+dqi*facu(iax)*ri1
c      dswqws=dswqws+drmj*ri1
c      dswqps=dswqps+dqai*rmj*ri1
c      dswpss=dswpss-drmi*ri1
c      dswpws=dswpws+qai*drmj*ri1
c      endif
      zxx=(-half)*(fri*rx+fmi*pix+fmj*pjx)*rwc(iax,1)
      zyx=(-half)*(fri*rx+fmi*pix+fmj*pjx)*rwc(iax,2)
      zzx=(-half)*(fri*rx+fmi*pix+fmj*pjx)*rwc(iax,3)
      zxy=(-half)*(fri*ry+fmi*piy+fmj*pjy)*rwc(iax,1)
      zyy=(-half)*(fri*ry+fmi*piy+fmj*pjy)*rwc(iax,2)
      zzy=(-half)*(fri*ry+fmi*piy+fmj*pjy)*rwc(iax,3)
      zxz=(-half)*(fri*rz+fmi*piz+fmj*pjz)*rwc(iax,1)
      zyz=(-half)*(fri*rz+fmi*piz+fmj*pjz)*rwc(iax,2)
      zzz=(-half)*(fri*rz+fmi*piz+fmj*pjz)*rwc(iax,3)
      zw(1,1,ipsw)=zw(1,1,ipsw)+zxx
      zw(2,1,ipsw)=zw(2,1,ipsw)+zyx
      zw(3,1,ipsw)=zw(3,1,ipsw)+zzx
      zw(1,2,ipsw)=zw(1,2,ipsw)+zxy
      zw(2,2,ipsw)=zw(2,2,ipsw)+zyy
      zw(3,2,ipsw)=zw(3,2,ipsw)+zzy
      zw(1,3,ipsw)=zw(1,3,ipsw)+zxz
      zw(2,3,ipsw)=zw(2,3,ipsw)+zyz
      zw(3,3,ipsw)=zw(3,3,ipsw)+zzz
      zs(isf,1,1,ipsw)=zs(isf,1,1,ipsw)+zxx
      zs(isf,2,1,ipsw)=zs(isf,2,1,ipsw)+zyx
      zs(isf,3,1,ipsw)=zs(isf,3,1,ipsw)+zzx
      zs(isf,1,2,ipsw)=zs(isf,1,2,ipsw)+zxy
      zs(isf,2,2,ipsw)=zs(isf,2,2,ipsw)+zyy
      zs(isf,3,2,ipsw)=zs(isf,3,2,ipsw)+zzy
      zs(isf,1,3,ipsw)=zs(isf,1,3,ipsw)+zxz
      zs(isf,2,3,ipsw)=zs(isf,2,3,ipsw)+zyz
      zs(isf,3,3,ipsw)=zs(isf,3,3,ipsw)+zzz
   21 continue
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
