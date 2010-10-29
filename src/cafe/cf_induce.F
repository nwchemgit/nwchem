      subroutine induce(lself,lpbcs,xw,xwm,pw,pwp,iwdt,iwz,
     + iwfr,iwto,jwfr,jwto,xs,xsm,ps,psp,
     + isga,isat,isdt,isgr,ismf,isml,isss,isq1,isq2,isq3,isgm,ishop,isz,
     + isfr,isto,jsfr,jsto,lpbc,lstptr,lseq)
c
c     in     log : lself           : true for box self interactions
c     in     r*8 : xw(mwm,3,mwa)   : solvent coordinates
c     in     r*8 : xwm(mwm,3)      : solvent molecule center of mass coordinates
c     out    r*8 : rtos(mwm)       : solvent distances to closes solute atom
c     in     int : iwdt(mwm)       : solvent dynamics type
c     out    int : iwz(mwm)        : solvent boundary type
c     in     int : iwfr,iwto       : first and last solvent molecule i
c     in     int : jwfr,jwto       : first and last solvent molecule j
c     in     r*8 : xs(msa,3)       : solute atom coordinates
c     in     r*8 : xsm(msm,3)      : solute molecule center of mass coordinates
c     in     int : isga(msa)       : solute global atom number
c     in     int : isat(msa)       : solute atom type
c     in     int : isdt(msa)       : solute dynamics type
c     in     int : isgr(msa)       : solute charge group
c     in     int : ismf(msa)       : solute molecule fraction
c     in     int : isml(msa)       : solute molecule
c     in     int : isss(msa)       : solute separation shifted scaling type
c     in     int : isq1(msa)       : solute charge type 1
c     in     int : isq2(msa)       : solute charge type 2
c     in     int : isq3(msa)       : solute charge type 3
c     out    int : isz(msa)        : solute boundary type 
c     in     int : isfr,isto       : first and last solute atom i
c     in     int : jsfr,jsto       : first and last solute atom j
c     in     log : lpbc            : flag to consider periodic boundary conditions
c     in/out int : lstptr          : list pointer
c
c     dimensions nwm,nwa and nsa need to have been given by a call to cf_initx
c
c $Id$
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),xwm(mwm,3)
      real*8 pw(mwm,3,mwa,2),pwp(mwm,3,mwa,2,2)
      real*8 xs(msa,3),xsm(msm,3)
      real*8 ps(msa,3,2),psp(msa,3,2,2)
      integer iwdt(mwm),iwz(mwm),isz(msa)
      integer isga(msa),isat(msa),isdt(msa),isgr(msa),ismf(msa)
      integer isml(msa),isss(msa),isq1(msa),isq2(msa),isq3(msa)
      integer isgm(msa),ishop(msa)
      integer lseq(mseq)
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer lstptr
      logical lself,lpbc,lpbcs
c
      integer nwloc,nsloc,nwnon,nsnon,npairs
      integer lptr
c
      if(lself) then
      jwfr=iwfr
      jwto=iwto
      jsfr=isfr
      jsto=isto
      endif
c
      nwloc=iwto-iwfr+1
      if(iwfr.eq.0.or.iwto.lt.iwfr) nwloc=0
      nwnon=jwto-jwfr+1
      if(jwfr.eq.0.or.jwto.lt.jwfr) nwnon=0
      nsloc=isto-isfr+1
      if(isfr.eq.0.or.isto.lt.isfr) nsloc=0
      nsnon=jsto-jsfr+1
      if(jsfr.eq.0.or.jsto.lt.jsfr) nsnon=0
c
c     list [ list(1)=int_mb(lptr) ] with offset ndxp
c
c     list( 1) : index to pairlist for solvent-solvent
c     list( 2) : index to pairlist for solute-solvent
c     list( 3) : index to pairlist for solvent-solute
c     list( 4) : index to pairlist for solute-solute
c     list( 5) : index to pairlist for solute bonds
c     list( 6) : index to pairlist for solute angles
c     list( 7) : index to pairlist for solute dihedrals
c     list( 8) : index to pairlist for solute impropers
c     list( 9) : index to pairlist for solute third neighbors
c     list(10) : index to pairlist for solute excluded pairs
c
c
c     pairlists
c     ---------
c
      if(lpair) call pairs(lself,lpbcs,xw,xwm,iwdt,iwz,
     + iwfr,iwto,jwfr,jwto,xs,xsm,
     + isga,isat,isdt,isgr,isgm,ismf,isml,isss,isq1,isq2,isq3,ishop,isz,
     + isfr,isto,jsfr,jsto,lpbc,lstptr,lseq)
c
c     induced dipoles
c     ---------------
c
      ndxp=lstptr
      lptr=i_list+ndxp+24
c
      if(lself.and.nwloc.gt.0) then
c      call cf_fw(iwfr,iwto,xw,fw,iwdt,int_mb(i_iwa),int_mb(i_iwq),
c     + lpbc,dbl_mb(i_vdw),dbl_mb(i_chg),mbt(1),numb(1),mbp(1),
c     + int_mb(i_ibnd(1)),dbl_mb(i_bnd(1)),dbl_mb(i_rbnd(1)),
c     + mht(1),numh(1),mhp(1),
c     + int_mb(i_iang(1)),dbl_mb(i_ang(1)),dbl_mb(i_rang(1)),
c     + mdt(1),numd(1),mdp(1),
c     + int_mb(i_idih(1)),dbl_mb(i_dih(1)),dbl_mb(i_rdih(1)),
c     + mit(1),numi(1),mip(1),
c     + int_mb(i_iimp(1)),dbl_mb(i_imp(1)),dbl_mb(i_rimp(1)),
c     + mtt(1),numt(1),int_mb(i_itrd(1)),
c     + mxt(1),numx(1),int_mb(i_ixcl(1)))
      endif
c
c     solvent-solvent polarization
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
      call cf_dww(xw,xwm,pw,pwp,iwfr,nwloc,lpbc,
     + dbl_mb(i_chg),int_mb(i_iwq),
     + int_mb(i_s2i1),
     + int_mb(lptr+1),int_mb(lptr+1+2*nwloc),int_mb(lptr+1+4*nwloc),
     + dbl_mb(i_s3r1),dbl_mb(i_smr1),dbl_mb(i_smr2),dbl_mb(i_s3r2),
     + dbl_mb(i_smr3),dbl_mb(i_smr4),dbl_mb(i_smr5),dbl_mb(i_smr6))
      lptr=lptr+4*nwloc+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute-solvent polarization
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
      call cf_dsw(xs,xsm,ps,psp,isdt,ismf,isml,isq1,isfr,nsloc,
     + xw,xwm,pw,pwp,int_mb(i_iwq),lpbc,dbl_mb(i_chg),int_mb(i_s2i1),
     + int_mb(lptr+1),int_mb(lptr+1+2*nsloc),int_mb(lptr+1+4*nsloc),
     + dbl_mb(i_s3r1),dbl_mb(i_smr1),dbl_mb(i_smr2),dbl_mb(i_s3r2),
     + dbl_mb(i_smr3),dbl_mb(i_smr4),dbl_mb(i_smr5),dbl_mb(i_smr6),
     + int_mb(i_s1i1))
      lptr=lptr+4*nsloc+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solvent-solute forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
      call cf_dsw(xs,xsm,ps,psp,isdt,ismf,isml,isq1,jsfr,nsnon,
     + xw,xwm,pw,pwp,int_mb(i_iwq),lpbc,dbl_mb(i_chg),int_mb(i_s2i1),
     + int_mb(lptr+1),int_mb(lptr+1+2*nsloc),int_mb(lptr+1+4*nsloc),
     + dbl_mb(i_s3r1),dbl_mb(i_smr1),dbl_mb(i_smr2),dbl_mb(i_s3r2),
     + dbl_mb(i_smr3),dbl_mb(i_smr4),dbl_mb(i_smr5),dbl_mb(i_smr6),
     + int_mb(i_s1i1))
      lptr=lptr+4*nsnon+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute-solute forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
      call cf_dss(xs,xsm,ps,psp,ismf,isml,isq2,isq3,isfr,nsloc,
     + lpbc,dbl_mb(i_chg),int_mb(i_s2i1),
     + int_mb(lptr+1),int_mb(lptr+1+2*nsloc),int_mb(lptr+1+4*nsloc),
     + dbl_mb(i_smr1),dbl_mb(i_smr2),dbl_mb(i_s3r1),dbl_mb(i_s1r1),
     + dbl_mb(i_s1r2),dbl_mb(i_s1r3),dbl_mb(i_s3r2),dbl_mb(i_s1r6),
     + dbl_mb(i_smr3),dbl_mb(i_smr4),dbl_mb(i_s1r5),
     + int_mb(i_s1i1),int_mb(i_s1i2),int_mb(i_s1i3),int_mb(i_s1i4),
     + dbl_mb(i_smr3),dbl_mb(i_smr4))
      lptr=lptr+4*nsloc+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute bond forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
c      call cf_fsb(npairs,int_mb(lptr+1),mbt(2),mbp(2),
c     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
c     + max(isto,jsto),msa,isga,ismf,isdt,isq1,dbl_mb(i_chg),
c     + xs,fs,esb,esp,lpbc)
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute angle forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
c      call cf_fsh(npairs,int_mb(lptr+1),mht(2),mhp(2),
c     + int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)),
c     + max(isto,jsto),msa,isga,ismf,isdt,isq1,dbl_mb(i_chg),
c     + xs,fs,esh,esp,lpbc,lupden)
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute torsion forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
c      call cf_fsd(npairs,int_mb(lptr+1),mdt(2),mdp(2),
c     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)),
c     + max(isto,jsto),msa,isga,ismf,isdt,
c     + xs,fs,esd,esp,lpbc,lupden)
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute improper dihedral forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
c      call cf_fso(npairs,int_mb(lptr+1),mdt(2),mdp(2),
c     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)),
c     + max(isto,jsto),msa,isga,ismf,isdt,
c     + xs,fs,eso,esp,lpbc,lupden)
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute third neighbor forces
c
      npairs=int_mb(lptr)
      if(npairs.gt.0) then
c      call cf_fst(npairs,int_mb(lptr+1),mtt(2),int_mb(i_itrd(2)),
c     + dbl_mb(i_vdw),dbl_mb(i_chg),
c     + max(isto,jsto),msa,isat,isga,isml,isdt,isq2,
c     + xs,fs,esp,essl,essq,essr,lpbc)
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
      ndxp=lptr-i_list
c
      return
      end
