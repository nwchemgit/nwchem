      subroutine argos_cafe_fss(xs,xsm,fs,zs,
     + isga,isat,isdt,ismf,isml,isss,isq2,isq3,isgm,
     + isfrom,nums,lpbc,lpbcs,ess,fss,esa,
     + vdw,chg,iass,
     + lssndx,lssjpt,lssin,lssj,
     + xi,xj,rwx,rwi1,rwi2,rwi6,rwc,f,fi,fj,facu,
     + rw,isal,jsal,jmal,jfal,isrx,qsa2,qsa3,dera,lda,rda,uda,lseq)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 rtmp
      real*8 xs(msa,3),xsm(msm,3),fs(msa,3,2)
      real*8 zs(msf,3,3,2),ess(msf,msf,mpe,2)
      real*8 fss(msf,msf,3,2)
      real*8 esa(nsa)
      integer isga(msa),isat(msa),isdt(msa),ismf(msa)
      integer isml(msa),isss(msa),isq2(msa),isq3(msa)
      integer isgm(msa),lseq(mseq)
c
      real*8 vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      logical lpbc,lpbcs
c
      real*8 xi(mscr,3),xj(mscr,3),rwx(mscr,3),rwi1(mscr)
      real*8 rwi2(mscr),rwi6(mscr),rwc(mscr,3),rw(mscr)
      real*8 f(mscr),fi(mscr,3),fj(mscr,3)
      real*8 qsa2(mscr),qsa3(mscr)
      integer isal(mscr),jsal(mscr),jmal(mscr),jfal(mscr),isrx(mscr)
      integer lssj(*)
      real*8 facu(mscr)
      integer nums
      integer lssndx(0:msa,2),lssjpt(nums,2),lssin(nums,2)
      integer iass(mat,mat),nsslen(2)
c
      real*8 dera(6,nsatot)
      integer lda(16,*)
      real*8 rda(11,*),uda(4,*)
c
      integer isa,jsa,i,isf,jsf,ix
      integer isfr,isfrom,ism,jsm
      integer ipss,number,isslen,nax,jsaptr
      integer jnum,lssptr,iax
      real*8 dercon
c
      real*8 c6,c12,cf6,cf12
      real*8 c64,c124
      real*8 q14,sumen1,sumen2,sumen3
      real*8 etermq,eterml
      integer istt,jstt
      real*8 q,ferfc,fderfc
c
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
c
      if(nfhop.eq.0) then
      do 112 i=1,msa
      if(isq2(i).le.0.or.isq3(i).le.0.or.
     + isq2(i).gt.mqt.or.isq3(i).gt.mqt) goto 113
      qsa2(i)=chg(isq2(i),1,iset)
      qsa3(i)=chg(isq3(i),1,iset)
  112 continue
  113 continue
      else
      do 1112 i=1,msa
      if(isq2(i).le.0.or.isq3(i).le.0.or.
     + isq2(i).gt.mqt.or.isq3(i).gt.mqt) goto 1113
      qsa2(i)=chg(isq2(i),1,lseq(isgm(i)))
      qsa3(i)=chg(isq3(i),1,lseq(isgm(i)))
 1112 continue
 1113 continue
      endif
c
c     solute non-bonded interactions
c     ==============================
c
      isfr=isfrom-1
c
c     loop over short and long range pairlists
c
      do 11 ipss=1,lpss
c
c     evaluate outer index array
c
      nsslen(ipss)=0
      lssndx(0,ipss)=0
      number=0
      do 12 isa=1,nums
      if(number+lssin(isa,ipss).gt.mscr.or.
     + (ismf(isfr+isa).ne.ismf(isfr+isa-1).and.
     + number.gt.0)) then
      nsslen(ipss)=nsslen(ipss)+1
      lssndx(nsslen(ipss),ipss)=isa-1
      number=0
      endif
      number=number+lssin(isa,ipss)
   12 continue
      if(number.gt.0) then
      nsslen(ipss)=nsslen(ipss)+1
      lssndx(nsslen(ipss),ipss)=nums
      endif
c
c     loop over number of cycles to complete pairlists
c
      do 13 isslen=1,nsslen(ipss)
c
      etermq=zero
      eterml=zero
c
      nax=0
      isf=ismf(isfr+lssndx(isslen,ipss))
c
c     collect coordinates into workarrays
c
      do 14 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      ism=isml(isfr+isa)
      if(lpbc.or.lpbcs) then
      if(ipbtyp.eq.1) then
      do 15 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      rwc(nax+jnum,1)=xs(isfr+isa,1)-xs(lssptr,1)
      rwc(nax+jnum,2)=xs(isfr+isa,2)-xs(lssptr,2)
      rwc(nax+jnum,3)=xs(isfr+isa,3)-xs(lssptr,3)
      isrx(nax+jnum)=0
   15 continue
      else
      do 115 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      jsm=isml(lssptr)
      rwc(nax+jnum,1)=xsm(ism,1)-xsm(jsm,1)
      rwc(nax+jnum,2)=xsm(ism,2)-xsm(jsm,2)
      rwc(nax+jnum,3)=xsm(ism,3)-xsm(jsm,3)
      isrx(nax+jnum)=0
  115 continue
      endif
      call argos_cafe_pbc(0,rwc,mscr,rwx,mscr,nax,1,lssin(isa,ipss))
      endif
      do 16 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      jsf=ismf(lssptr)
      isal(nax+jnum)=isfr+isa
      jsal(nax+jnum)=lssptr
      jfal(nax+jnum)=jsf
      jmal(nax+jnum)=0
      jsm=isml(lssptr)
      if(ism.ne.jsm) jmal(nax+jnum)=1
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
c
      isrx(nax+jnum)=0
c
      if(lssscl) then
c
      istt=iand(isss(isfr+isa),48)
      jstt=iand(isss(lssptr),48)
      if(ism.ne.jsm) then
      if(istt.eq.16.or.jstt.eq.16) isrx(nax+jnum)=-1
      if(istt.eq.32.or.jstt.eq.32) isrx(nax+jnum)=1
      endif
c
      istt=iand(isss(isfr+isa),384)
      jstt=iand(isss(lssptr),384)
      if(istt.eq.128.or.jstt.eq.128) isrx(nax+jnum)=-2
      if(istt.eq.256.or.jstt.eq.256) isrx(nax+jnum)=2
c
      istt=iand(isss(isfr+isa),384)
      jstt=iand(isss(lssptr),384)
      if(istt.eq.128.and.jstt.eq.256) isrx(nax+jnum)=999
      if(istt.eq.256.and.jstt.eq.128) isrx(nax+jnum)=999
c
c      write(*,'(5i5)')
c     + isga(isfr+isa),isga(lssptr),istt,jstt,isrx(nax+jnum)
c
      endif
c
   16 continue
c
      do 17 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      facu(nax+jnum)=zero
      if(iand(isdt(isfr+isa),mdynam).eq.ldynam.or.
     + iand(isdt(lssptr),mdynam).eq.ldynam) facu(nax+jnum)=one
c      if((iand(isdt(isfr+isa),mdynam).eq.ldynam.and.
c     + iand(isdt(lssptr),mdynam).ne.ldynam) .or.
c     + (iand(isdt(isfr+isa),mdynam).ne.ldynam.and.
c     + iand(isdt(lssptr),mdynam).eq.ldynam)) facu(nax+jnum)=half
      if(includ.eq.1) facu(nax+jnum)=one
   17 continue
c
      if(.not.lpbc.and..not.lpbcs) then
      do 18 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      xi(nax+jnum,1)=xs(isfr+isa,1)
      xi(nax+jnum,2)=xs(isfr+isa,2)
      xi(nax+jnum,3)=xs(isfr+isa,3)
      xj(nax+jnum,1)=xs(lssptr,1)
      xj(nax+jnum,2)=xs(lssptr,2)
      xj(nax+jnum,3)=xs(lssptr,3)
      isal(nax+jnum)=isfr+isa
      jsal(nax+jnum)=lssptr
   18 continue
      else
      do 19 jnum=1,lssin(isa,ipss)
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
      isal(nax+jnum)=isfr+isa
      jsal(nax+jnum)=lssptr
   19 continue
      endif
c
      nax=nax+lssin(isa,ipss)
   14 continue
c
c
c     evaluate electrostatic energies and forces
c
c      etermq=zero
      if(.not.lssscl) then
      if(ipme.eq.0.or.isolvo.ne.0) then
      do 24 iax=1,nax
      f(iax)=zero
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      isa=isal(iax)
      jsa=jsal(iax)
c
      if(jfal(iax).ne.isf) then
      q=qsa2(isa)*qsa2(jsa)
      else
      q=qsa3(isa)*qsa3(jsa)
      endif
c
      rw(iax)=facu(iax)*q*rwi1(iax)
      f(iax)=f(iax)+q*rwi1(iax)*rwi2(iax)
      if(ireact.ne.0) then
      ess(isf,jfal(iax),6,ipss)=ess(isf,jfal(iax),6,ipss)+
     + facu(iax)*q*rffss/rwi2(iax)
      if(npener.ne.0) then
      esa(isga(isa))=esa(isga(isa))+half*facu(iax)*q*rffss/rwi2(iax)
      esa(isga(jsa))=esa(isga(jsa))+half*facu(iax)*q*rffss/rwi2(iax)
      endif
      f(iax)=f(iax)-two*q*rffss
      endif
cx      if(ihess.gt.0) then
cx      h(iax)=three*q*rwi1(iax)*rwi2(iax)*rwi2(iax)
cx      endif
   24 continue
      else
      do 25 iax=1,nax
      f(iax)=zero
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).eq.999) then
      rwi2(iax)=zero
      isrx(iax)=0
      endif
      rwi1(iax)=sqrt(rwi2(iax))
      isa=isal(iax)
      jsa=jsal(iax)
      if(jfal(iax).ne.isf) then
      q=qsa2(isa)*qsa2(jsa)
      else
      q=qsa3(isa)*qsa3(jsa)
      endif
c
      ferfc=erfc(ealpha/rwi1(iax))
      fderfc=ealpha*derfc(ealpha/rwi1(iax))
      rw(iax)=ferfc*facu(iax)*q*rwi1(iax)
      f(iax)=f(iax)+q*rwi2(iax)*(ferfc*rwi1(iax)-fderfc)
   25 continue
      endif
      else
      if(ipme.eq.0.or.isolvo.ne.0) then
      do 26 iax=1,nax
      f(iax)=zero
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).eq.999) then
      rwi2(iax)=zero
      isrx(iax)=0
      endif
      rtmp=rwi2(iax)
      if(isrx(iax).gt.0) then
      rwi2(iax)=one/(one/rwi2(iax)+shift0(1))
      elseif(isrx(iax).lt.0) then
      rwi2(iax)=one/(one/rwi2(iax)+shift1(1))
      endif
      rwi1(iax)=sqrt(rwi2(iax))
      isa=isal(iax)
      jsa=jsal(iax)
c      write(*,'(3i5,4f12.6)')
c     + isga(isa),isga(jsa),isrx(iax),shift0(1),shift1(1),rtmp,rwi2(iax)
      if(jfal(iax).ne.isf) then
      q=qsa2(isa)*qsa2(jsa)
      else
      q=qsa3(isa)*qsa3(jsa)
      endif
      rw(iax)=facu(iax)*q*rwi1(iax)
      f(iax)=f(iax)+q*rwi1(iax)*rwi2(iax)
      if(ireact.ne.0) then
      rw(iax)=rw(iax)+facu(iax)*q*rffss/rwi2(iax)
      f(iax)=f(iax)-two*q*rffss
      endif
   26 continue
      else
      do 126 iax=1,nax
      f(iax)=zero
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).eq.999) then
      rwi2(iax)=zero
      isrx(iax)=0
      endif
      if(isrx(iax).gt.0) then
      rwi2(iax)=one/(one/rwi2(iax)+shift0(1))
      elseif(isrx(iax).lt.0) then
      rwi2(iax)=one/(one/rwi2(iax)+shift1(1))
      endif
      rwi1(iax)=sqrt(rwi2(iax))
      isa=isal(iax)
      jsa=jsal(iax)
      if(jfal(iax).ne.isf) then
      q=qsa2(isa)*qsa2(jsa)
      else
      q=qsa3(isa)*qsa3(jsa)
      endif
      ferfc=erfc(ealpha/rwi1(iax))
      fderfc=ealpha*derfc(ealpha/rwi1(iax))
      rw(iax)=ferfc*facu(iax)*q*rwi1(iax)
      f(iax)=f(iax)+q*rwi2(iax)*(ferfc*rwi1(iax)-fderfc)
      if(ireact.ne.0) then
      rw(iax)=rw(iax)+facu(iax)*q*rffss/rwi2(iax)
      f(iax)=f(iax)-two*q*rffss
      endif
  126 continue
      endif
      endif
c
c     accumulate electrostatic energies per solute molecule
c
c      etermq=zero
      do 27 iax=1,nax
      if(npener.ne.0) then
      esa(isga(isal(iax)))=esa(isga(isal(iax)))+half*rw(iax)
      esa(isga(jsal(iax)))=esa(isga(jsal(iax)))+half*rw(iax)
      endif
      ess(isf,jfal(iax),6,ipss)=ess(isf,jfal(iax),6,ipss)+rw(iax)
      etermq=etermq+rw(iax)
   27 continue
c
c      do 27 jsf=1,msf
c      sumen=zero
c      do 28 iax=1,nax
c      if(jfal(iax).eq.jsf) sumen=sumen+rw(iax)
c      if(npener.ne.0) then
c      esa(isga(isal(iax)))=esa(isga(isal(iax)))+half*rw(iax)
c      esa(isga(jsal(iax)))=esa(isga(jsal(iax)))+half*rw(iax)
c      endif
c   28 continue
c      ess(isf,jsf,6,ipss)=ess(isf,jsf,6,ipss)+sumen
c      etermq=etermq+sumen
c   27 continue
c
c
c
c     Lennard-Jones energies and forces
c     =================================
c
      do 29 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      rwi6(iax)=rwi2(iax)*rwi2(iax)*rwi2(iax)
      c6=vdw(isat(isa),isat(jsa),1,iset)
      c12=vdw(isat(isa),isat(jsa),3,iset)
      cf6=six*c6
      cf12=twelve*c12
      rw(iax)=facu(iax)*(c12*rwi6(iax)-c6)*rwi6(iax)
      f(iax)=f(iax)+(cf12*rwi6(iax)-cf6)*rwi6(iax)*rwi2(iax)
cx      if(ihess.gt.0) then
cx      h(iax)=h(iax)+(forten*cf12*rwi6(iax)-eight*cf6)*rwi6(iax)*
cx     + rwi2(iax)*rwi2(iax)
cx      endif
   29 continue
c
c     accumulate Lennard-Jones energies per solute molecule
c
c      eterml=zero
c      etermq=zero
      do 30 iax=1,nax
      if(npener.ne.0) then
      esa(isga(isal(iax)))=esa(isga(isal(iax)))+half*rw(iax)
      esa(isga(jsal(iax)))=esa(isga(jsal(iax)))+half*rw(iax)
      endif
      ess(isf,jfal(iax),5,ipss)=ess(isf,jfal(iax),5,ipss)+rw(iax)
      eterml=eterml+rw(iax)
   30 continue
c
c      do 30 jsf=1,msf
c      sumen=zero
c      do 31 iax=1,nax
c      if(jfal(iax).eq.jsf) sumen=sumen+rw(iax)
c      if(npener.ne.0) then
c      esa(isga(isal(iax)))=esa(isga(isal(iax)))+half*rw(iax)
c      esa(isga(jsal(iax)))=esa(isga(jsal(iax)))+half*rw(iax)
c      endif
c   31 continue
c      ess(isf,jsf,5,ipss)=ess(isf,jsf,5,ipss)+sumen
c      eterml=eterml+sumen
c   30 continue
c
c     evaluate and accumulate the solute-solute virial contributions
c     allow virial contributions from interactions between a solute
c     molecule and its own image
c
      do 132 ix=1,3
      do 32 jsf=1,msf
      sumen1=zero
      sumen2=zero
      sumen3=zero
      do 33 iax=1,nax
cx      if(jfal(iax).eq.jsf.and.jmal(iax).eq.1) then
      if(jfal(iax).eq.jsf) then
      sumen1=sumen1-half*f(iax)*rwx(iax,1)*rwc(iax,ix)
      sumen2=sumen2-half*f(iax)*rwx(iax,2)*rwc(iax,ix)
      sumen3=sumen3-half*f(iax)*rwx(iax,3)*rwc(iax,ix)
      endif
   33 continue
      zs(isf,ix,1,ipss)=zs(isf,ix,1,ipss)+sumen1
      zs(jsf,ix,1,ipss)=zs(jsf,ix,1,ipss)+sumen1
      zs(isf,ix,2,ipss)=zs(isf,ix,2,ipss)+sumen2
      zs(jsf,ix,2,ipss)=zs(jsf,ix,2,ipss)+sumen2
      zs(isf,ix,3,ipss)=zs(isf,ix,3,ipss)+sumen3
      zs(jsf,ix,3,ipss)=zs(jsf,ix,3,ipss)+sumen3
   32 continue
  132 continue
c
c     generate radial distribution functions
c
c      if(ifstep-1.eq.((ifstep-1)/nfrdf)*nfrdf.and.ngrss.gt.0) then
c      do 34 iax=1,nax
c      isa=isal(iax)
c      jsa=jsal(iax)
c      do 35 igc=1,ngc
c      if(ngt(igc).eq.3) then
c      if((isga(isa).eq.iagc(igc).and.
c     + isga(jsa).eq.jagc(igc)).or.
c     + (isga(isa).eq.jagc(igc).and.
c     + isga(jsa).eq.iagc(igc))) then
c      igr=igrc(igc)
c      indx=int(one/(rwi1(iax)*drdf))
c      if(indx.gt.ngl) indx=ngl
c      rdf(indx,igr)=rdf(indx,igr)+rdfvol
c      endif
c      endif
c   35 continue
c   34 continue
c      endif
c
c     accumulate forces into solute force arrays
c
      nax=0
      do 36 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      do 37 jnum=1,lssin(isa,ipss)
      lssptr=lssj(jsaptr+jnum)
      fs(isfr+isa,1,ipss)=fs(isfr+isa,1,ipss)+
     + f(nax+jnum)*rwx(nax+jnum,1)
      fs(isfr+isa,2,ipss)=fs(isfr+isa,2,ipss)+
     + f(nax+jnum)*rwx(nax+jnum,2)
      fs(isfr+isa,3,ipss)=fs(isfr+isa,3,ipss)+
     + f(nax+jnum)*rwx(nax+jnum,3)
      fs(lssptr,1,ipss)=fs(lssptr,1,ipss)-f(nax+jnum)*rwx(nax+jnum,1)
      fs(lssptr,2,ipss)=fs(lssptr,2,ipss)-f(nax+jnum)*rwx(nax+jnum,2)
      fs(lssptr,3,ipss)=fs(lssptr,3,ipss)-f(nax+jnum)*rwx(nax+jnum,3)
      isf=ismf(isfr+isa)
      jsf=ismf(lssptr)
      fss(isf,jsf,1,ipss)=fss(isf,jsf,1,ipss)+
     + f(nax+jnum)*rwx(nax+jnum,1)
      fss(isf,jsf,2,ipss)=fss(isf,jsf,2,ipss)+
     + f(nax+jnum)*rwx(nax+jnum,2)
      fss(isf,jsf,3,ipss)=fss(isf,jsf,3,ipss)+
     + f(nax+jnum)*rwx(nax+jnum,3)
   37 continue
cx      if(ihess.gt.0) then
cx      do 137 jnum=1,lssin(isa,ipss)
cx      lssptr=lssj(jsaptr+jnum)
cx      hs(isfr+isa,1,ipss)=hs(isfr+isa,1,ipss)-f(nax+jnum)+
cx     + h(nax+jnum)**rwx(nax+jnum,1)*rwx(nax+jnum,1)
cx      hs(isfr+isa,2,ipss)=hs(isfr+isa,2,ipss)+
cx     + h(nax+jnum)*rwx(nax+jnum,1)*rwx(nax+jnum,2)
cx      hs(isfr+isa,3,ipss)=hs(isfr+isa,3,ipss)+
cx     + h(nax+jnum)*rwx(nax+jnum,1)*rwx(nax+jnum,3)
cx      hs(isfr+isa,4,ipss)=hs(isfr+isa,4,ipss)-f(nax+jnum)+
cx     + h(nax+jnum)*rwx(nax+jnum,2)*rwx(nax+jnum,2)
cx      hs(isfr+isa,5,ipss)=hs(isfr+isa,5,ipss)+
cx     + h(nax+jnum)*rwx(nax+jnum,2)*rwx(nax+jnum,3)
cx      hs(isfr+isa,6,ipss)=hs(isfr+isa,6,ipss)-f(nax+jnum)+
cx     + h(nax+jnum)*rwx(nax+jnum,3)*rwx(nax+jnum,3)
cx      hs(lssptr,1,ipss)=hs(lssptr,1,ipss)+f(nax+jnum)-
cx     + h(nax+jnum)*rwx(nax+jnum,1)*rwx(nax+jnum,1)
cx      hs(lssptr,2,ipss)=hs(lssptr,2,ipss)-
cx     + h(nax+jnum)*rwx(nax+jnum,1)*rwx(nax+jnum,2)
cx      hs(lssptr,3,ipss)=hs(lssptr,3,ipss)-
cx     + h(nax+jnum)*rwx(nax+jnum,1)*rwx(nax+jnum,3)
cx      hs(lssptr,4,ipss)=hs(lssptr,4,ipss)+f(nax+jnum)-
cx     + h(nax+jnum)*rwx(nax+jnum,2)*rwx(nax+jnum,2)
cx      hs(lssptr,5,ipss)=hs(lssptr,5,ipss)-
cx     + h(nax+jnum)*rwx(nax+jnum,2)*rwx(nax+jnum,3)
cx      hs(lssptr,6,ipss)=hs(lssptr,6,ipss)+f(nax+jnum)-
cx     + h(nax+jnum)*rwx(nax+jnum,3)*rwx(nax+jnum,3)
cx  137 continue
cx      endif
      nax=nax+lssin(isa,ipss)
   36 continue
c
c     thermodynamic integration
c
      if(ithint) then
      if(ith(14)) then
      nax=0
      do 38 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
c
      if(.not.lssscl) then
      do 39 jnum=1,lssin(isa,ipss)
      jsa=lssj(jsaptr+jnum)
      dercon=(vdw(isat(isfr+isa),isat(jsa),3,4)*rwi6(nax+jnum)
     + -vdw(isat(isfr+isa),isat(jsa),1,4))*rwi6(nax+jnum)
      deriv(15,ipss)=deriv(15,ipss)+dercon
      if(npgdec.gt.1) then
      dera(3,isga(isa))=dera(3,isga(isa))+half*dercon
      dera(3,isga(jsa))=dera(3,isga(jsa))+half*dercon
      endif
   39 continue
      else
      do 40 jnum=1,lssin(isa,ipss)
      jsa=lssj(jsaptr+jnum)
      dercon=(vdw(isat(isfr+isa),isat(jsa),3,4)*rwi6(nax+jnum)
     + -vdw(isat(isfr+isa),isat(jsa),1,4))*rwi6(nax+jnum)
      if(isrx(nax+jnum).gt.0) then
      c64=three*vdw(isat(isfr+isa),isat(jsa),1,iset)
      c124=six*vdw(isat(isfr+isa),isat(jsa),3,iset)
      dercon=dercon+shift0(4)*
     + rwi2(nax+jnum)*rwi6(nax+jnum)*(c64-c124*rwi6(nax+jnum))
      elseif(isrx(nax+jnum).lt.0) then
      c64=three*vdw(isat(isfr+isa),isat(jsa),1,iset)
      c124=six*vdw(isat(isfr+isa),isat(jsa),3,iset)
      dercon=dercon+shift1(4)*
     + rwi2(nax+jnum)*rwi6(nax+jnum)*(c64-c124*rwi6(nax+jnum))
      endif
      deriv(15,ipss)=deriv(15,ipss)+dercon
c      write(*,'(a,3i5,4f12.6)') 'gv ',
c     + isga(isfr+isa),isga(jsa),isrx(nax+jnum),shift0(4),shift1(4),
c     + dercon,deriv(15,ipss)
      if(npgdec.gt.1) then
      dera(3,isga(isfr+isa))=dera(3,isga(isfr+isa))+half*dercon
      dera(3,isga(jsa))=dera(3,isga(jsa))+half*dercon
      endif
   40 continue
      endif
c
      nax=nax+lssin(isa,ipss)
   38 continue
      endif
c
      if(ith(16)) then
      nax=0
      do 41 isa=lssndx(isslen-1,ipss)+1,lssndx(isslen,ipss)
      jsaptr=lssjpt(isa,ipss)-1
      ism=isml(isfr+isa)
      if(ipme.eq.0) then
      if(.not.lssscl) then
      do 42 jnum=1,lssin(isa,ipss)
      jsa=lssj(jsaptr+jnum)
      if(isml(jsa).ne.ism) then
      dercon=(chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,4)
     + +chg(isq2(isfr+isa),1,4)*chg(isq2(jsa),1,iset))
      else
      dercon=(chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,4)
     + +chg(isq3(isfr+isa),1,4)*chg(isq3(jsa),1,iset))
      endif
      deriv(17,ipss)=deriv(17,ipss)+dercon*rwi1(nax+jnum)
      if(ireact.ne.0) then
      deriv(17,ipss)=deriv(17,ipss)+dercon*rffss/rwi2(nax+jnum)
      endif
c      write(*,'(a,3i5,4f12.6)') 'gq ',
c     + isga(isfr+isa),isga(jsa),isrx(nax+jnum),shift0(4),shift1(4),
c     + dercon,deriv(17,ipss)
      if(npgdec.gt.1) then
      dera(4,isga(isfr+isa))=dera(4,isga(isfr+isa))+
     + half*dercon*rwi1(nax+jnum)
      dera(4,isga(jsa))=dera(4,isga(jsa))+half*dercon*rwi1(nax+jnum)
      endif
   42 continue
      else
      do 43 jnum=1,lssin(isa,ipss)
      jsa=lssj(jsaptr+jnum)
      if(isml(jsa).ne.ism) then
      dercon=(chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,4)
     + +chg(isq2(isfr+isa),1,4)*chg(isq2(jsa),1,iset))
      if(isrx(nax+jnum).gt.0) then
      dercon=dercon-half*shift0(4)*
     + chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,iset)*rwi2(nax+jnum)
      elseif(isrx(nax+jnum).lt.0) then
      dercon=dercon-half*shift1(4)*
     + chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,iset)*rwi2(nax+jnum)
      endif
      else
      dercon=(chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,4)
     + +chg(isq3(isfr+isa),1,4)*chg(isq3(jsa),1,iset))
      if(isrx(nax+jnum).gt.1) then
      dercon=dercon-half*shift0(4)*
     + chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,iset)*rwi2(nax+jnum)
      elseif(isrx(nax+jnum).lt.-1) then
      dercon=dercon-half*shift1(4)*
     + chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,iset)*rwi2(nax+jnum)
      endif
      endif
      deriv(17,ipss)=deriv(17,ipss)+dercon*rwi1(nax+jnum)
      if(ireact.ne.0) then
      deriv(17,ipss)=deriv(17,ipss)+dercon*rffss/rwi2(nax+jnum)
      endif
      if(npgdec.gt.1) then
      dera(4,isga(isfr+isa))=dera(4,isga(isfr+isa))+
     + half*dercon*rwi1(nax+jnum)
      dera(4,isga(jsa))=dera(4,isga(jsa))+half*dercon*rwi1(nax+jnum)
      endif
   43 continue
      endif
      else
      if(.not.lssscl) then
      do 142 jnum=1,lssin(isa,ipss)
      jsa=lssj(jsaptr+jnum)
      if(isml(jsa).ne.ism) then
      dercon=(chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,4)
     + +chg(isq2(isfr+isa),1,4)*chg(isq2(jsa),1,iset))
      else
      dercon=(chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,4)
     + +chg(isq3(isfr+isa),1,4)*chg(isq3(jsa),1,iset))
      endif
      deriv(17,ipss)=deriv(17,ipss)+dercon*rwi1(nax+jnum)
      if(ireact.ne.0) then
      deriv(17,ipss)=deriv(17,ipss)+dercon*rffss/rwi2(nax+jnum)
      endif
      if(npgdec.gt.1) then
      dera(4,isga(isfr+isa))=dera(4,isga(isfr+isa))+
     + half*dercon*rwi1(nax+jnum)
      dera(4,isga(jsa))=dera(4,isga(jsa))+half*dercon*rwi1(nax+jnum)
      endif
  142 continue
      else
      do 143 jnum=1,lssin(isa,ipss)
      jsa=lssj(jsaptr+jnum)
      if(isml(jsa).ne.ism) then
      dercon=(chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,4)
     + +chg(isq2(isfr+isa),1,4)*chg(isq2(jsa),1,iset))
      if(isrx(nax+jnum).gt.0) then
      dercon=dercon-half*shift0(4)*
     + chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,iset)*rwi2(nax+jnum)
      elseif(isrx(nax+jnum).lt.0) then
      dercon=dercon-half*shift1(4)*
     + chg(isq2(isfr+isa),1,iset)*chg(isq2(jsa),1,iset)*rwi2(nax+jnum)
      endif
      else
      dercon=(chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,4)
     + +chg(isq3(isfr+isa),1,4)*chg(isq3(jsa),1,iset))
      if(isrx(nax+jnum).gt.1) dercon=dercon-half*shift0(4)*
     + chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,iset)*rwi2(nax+jnum)
      if(isrx(nax+jnum).lt.-1) dercon=dercon-half*shift1(4)*
     + chg(isq3(isfr+isa),1,iset)*chg(isq3(jsa),1,iset)*rwi2(nax+jnum)
      endif
      deriv(17,ipss)=deriv(17,ipss)+dercon*rwi1(nax+jnum)
      if(ireact.ne.0) then
      deriv(17,ipss)=deriv(17,ipss)+dercon*rffss/rwi2(nax+jnum)
      endif
      if(npgdec.gt.1) then
      dera(4,isga(isfr+isa))=dera(4,isga(isfr+isa))+
     + half*dercon*rwi1(nax+jnum)
      dera(4,isga(jsa))=dera(4,isga(jsa))+half*dercon*rwi1(nax+jnum)
      endif
  143 continue
      endif
      endif
      nax=nax+lssin(isa,ipss)
   41 continue
      endif
      endif
c
c     thermodynamic perturbation 1
c
      if(ipert2) then
      if(ip2(14)) then
      if(.not.lssscl) then
      do 44 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      ep2(ipss)=ep2(ipss)
     + +facu(iax)*(vdw(isat(isa),isat(jsa),3,2)*rwi6(iax)
     + -vdw(isat(isa),isat(jsa),1,2))*rwi6(iax)
   44 continue
      else
      do 45 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      rwi6(iax)=rwi2(iax)**3
      if(isrx(iax).gt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift0(1)+shift0(2)))**3
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift1(1)+shift1(2)))**3
      endif
      ep2(ipss)=ep2(ipss)
     + +facu(iax)*(vdw(isat(isa),isat(jsa),3,2)*rwi6(iax)
     + -vdw(isat(isa),isat(jsa),1,2))*rwi6(iax)
   45 continue
      endif
      ep2(ipss)=ep2(ipss)-eterml
      endif
      if(ip2(16).or.ip2(17)) then
      if(ipme.eq.0) then
      if(.not.lssscl) then
      do 46 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,2)*chg(isq2(jsa),1,2)
      else
      q14=chg(isq3(isa),1,2)*chg(isq3(jsa),1,2)
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
   46 continue
      else
      do 47 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,2)*chg(isq2(jsa),1,2)
      istt=0
      else
      q14=chg(isq3(isa),1,2)*chg(isq3(jsa),1,2)
      istt=1
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(2))
      elseif(isrx(iax).lt.-istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(2))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
   47 continue
      endif
      else
      if(.not.lssscl) then
      do 146 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,2)*chg(isq2(jsa),1,2)*
     + erfc(ealpha/rwi1(iax))
      else
      q14=chg(isq3(isa),1,2)*chg(isq3(jsa),1,2)*
     + erfc(ealpha/rwi1(iax))
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
  146 continue
      else
      do 147 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,2)*chg(isq2(jsa),1,2)*
     + erfc(ealpha/rwi1(iax))
      istt=0
      else
      q14=chg(isq3(isa),1,2)*chg(isq3(jsa),1,2)*
     + erfc(ealpha/rwi1(iax))
      istt=1
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(2))
      elseif(isrx(iax).lt.-istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(2))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep2(ipss)=ep2(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
  147 continue
      endif
      endif
      ep2(ipss)=ep2(ipss)-etermq
      endif
      endif
c
c     thermodynamic perturbation 2
c
      if(ipert3) then
      if(ip3(14)) then
      if(.not.lssscl) then
      do 48 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      ep3(ipss)=ep3(ipss)
     + +facu(iax)*(vdw(isat(isa),isat(jsa),3,3)*rwi6(iax)
     + -vdw(isat(isa),isat(jsa),1,3))*rwi6(iax)
   48 continue
      else
      do 49 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      rwi6(iax)=rwi2(iax)**3
      if(isrx(iax).gt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift0(1)+shift0(3)))**3
      elseif(isrx(iax).lt.0) then
      rwi6(iax)=(one/(one/rwi2(iax)-shift1(1)+shift1(3)))**3
      endif
      ep3(ipss)=ep3(ipss)
     + +facu(iax)*(vdw(isat(isa),isat(jsa),3,3)*rwi6(iax)
     + -vdw(isat(isa),isat(jsa),1,3))*rwi6(iax)
   49 continue
      endif
      ep3(ipss)=ep3(ipss)-eterml
      endif
      if(ip2(16).or.ip2(17)) then
      if(ipme.eq.0) then
      if(.not.lssscl) then
      do 50 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,3)*chg(isq2(jsa),1,3)
      else
      q14=chg(isq3(isa),1,3)*chg(isq3(jsa),1,3)
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
   50 continue
      else
      do 51 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,3)*chg(isq2(jsa),1,3)
      istt=0
      else
      q14=chg(isq3(isa),1,3)*chg(isq3(jsa),1,3)
      istt=1
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(3))
      elseif(isrx(iax).lt.-istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(3))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
   51 continue
      endif
      else
      if(.not.lssscl) then
      do 150 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,3)*chg(isq2(jsa),1,3)*
     + erfc(ealpha/rwi1(iax))
      else
      q14=chg(isq3(isa),1,3)*chg(isq3(jsa),1,3)*
     + erfc(ealpha/rwi1(iax))
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi2(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      rwi1(iax)=sqrt(rwi2(iax))
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
  150 continue
      else
      do 151 iax=1,nax
      isa=isal(iax)
      jsa=jsal(iax)
      if(jmal(iax).ne.0) then
      q14=chg(isq2(isa),1,3)*chg(isq2(jsa),1,3)*
     + erfc(ealpha/rwi1(iax))
      istt=0
      else
      q14=chg(isq3(isa),1,3)*chg(isq3(jsa),1,3)*
     + erfc(ealpha/rwi1(iax))
      istt=1
      endif
      rwx(iax,1)=xi(iax,1)-xj(iax,1)
      rwx(iax,2)=xi(iax,2)-xj(iax,2)
      rwx(iax,3)=xi(iax,3)-xj(iax,3)
      rwi6(iax)=one/(rwx(iax,1)**2+rwx(iax,2)**2+rwx(iax,3)**2)
      if(isrx(iax).gt.istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift0(3))
      elseif(isrx(iax).lt.-istt) then
      rwi6(iax)=one/(one/rwi6(iax)+shift1(3))
      endif
      rwi1(iax)=sqrt(rwi6(iax))
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rwi1(iax)
      if(ireact.ne.0) then
      ep3(ipss)=ep3(ipss)+facu(iax)*q14*rffss/rwi2(iax)
      endif
  151 continue
      endif
      endif
      ep3(ipss)=ep3(ipss)-etermq
      endif
      endif
   13 continue
   11 continue
c
c     accumulate radial distribution function contributions from
c     the excluded pairlist
c
c      if(ifstep-1.eq.((ifstep-1)/nfrdf)*nfrdf.and.ngrss.gt.0) then
c      do 52 isx=1,nsx
c      isa=idsx(isx)
c      jsa=jdsx(isx)
c      do 53 igc=1,ngc
c      if(ngt(igc).eq.3) then
c      if((isa.eq.iagc(igc).and.jsa.eq.jagc(igc)).or.
c     + (isa.eq.iagc(igc).and.jsa.eq.jagc(igc))) then
c      igr=igrc(igc)
c      indx=int(sqrt((xs(isa,1)-xs(jsa,1))**2+(xs(isa,2)-xs(jsa,2))**2+
c     + (xs(isa,3)-xs(jsa,3))**2)/drdf)
c      if(indx.gt.ngl) indx=ngl
c      rdf(indx,igr)=rdf(indx,igr)+rdfvol
c      endif
c      endif
c   53 continue
c   52 continue
c      endif
c
c
      return
      end
