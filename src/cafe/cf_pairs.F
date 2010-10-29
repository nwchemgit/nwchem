      subroutine pairs(lself,lpbcs,xw,xwm,iwdt,iwz,
     + iwfr,iwto,jwfr,jwto,xs,xsm,
     + isga,isat,isdt,isgr,isgm,ismf,isml,isss,isq1,isq2,isq3,ishop,isz,
     + isfr,isto,jsfr,jsto,lpbc,lstptr,lseq)
c
c     in     log : lself           : true for box self interactions
c     in     r*8 : xw(mwm,3,mwa)   : solvent coordinates
c     in     r*8 : xwm(mwm,3)      : solvent molecule center of mass coordinates
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
      real*8 xs(msa,3),xsm(msm,3)
      integer iwdt(mwm),iwz(mwm),isz(msa)
      integer isga(msa),isat(msa),isdt(msa),isgr(msa),ismf(msa)
      integer isml(msa),isss(msa),isq1(msa),isq2(msa),isq3(msa)
      integer ishop(msa),isgm(msa)
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer lstptr
      logical lself,lpbc,lpbcs
      integer lseq(mseq)
c
      integer nwloc,nsloc,nwnon,nsnon,npairs,mpairs
      integer lptr,lptrn
      integer nconst
c
      if(lself) then
      jwfr=iwfr
      jwto=iwto
      jsfr=isfr
      jsto=isto
      endif
c
      lstptr=ndxp
      lptrn=i_list+ndxp
      lptr=i_list+ndxp+24
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
c     pairlists
c     ---------
c
c     solvent-solvent pairlist
c
      npairs=0
      mpairs=maxl-(lptr-i_list)-4*nwloc-1
      if(nwloc.gt.0.and.nwnon.gt.0) then
      call cf_lww(lself,lpbc,xwm,iwdt,iwfr,iwto,jwfr,jwto,nwloc,
     + mpairs,npairs,int_mb(lptr+1),int_mb(lptr+1+2*nwloc),
     + int_mb(lptr+1+4*nwloc),
     + int_mb(i_s2i1),dbl_mb(i_s3r1),dbl_mb(i_s1r1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn)=lptr
      if(npairs.gt.0) then
      lptr=lptr+4*nwloc+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute-solvent pairlist
c
      npairs=0
      mpairs=maxl-(lptr-i_list)-4*nsloc-1
      if(nsloc.gt.0.and.nwnon.gt.0) then
      call cf_lsw(lpbc,lpbcs,
     + xs,isdt,isgr,isfr,isto,xwm,iwdt,iwz,jwfr,jwto,
     + nsloc,mpairs,npairs,int_mb(lptr+1),int_mb(lptr+1+2*nsloc),
     + int_mb(lptr+1+4*nsloc),
     + int_mb(i_s2i1),dbl_mb(i_s3r1),dbl_mb(i_s1r1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+1)=lptr
      if(npairs.gt.0) then
      lptr=lptr+4*nsloc+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solvent-solute pairlist
c
      npairs=0
      mpairs=maxl-(lptr-i_list)-4*nsnon-1
      if(nsnon.gt.0.and.nwloc.gt.0.and..not.lself) then
      call cf_lsw(lpbc,lpbcs,
     + xs,isdt,isgr,jsfr,jsto,xwm,iwdt,iwz,iwfr,iwto,
     + nsnon,mpairs,npairs,int_mb(lptr+1),int_mb(lptr+1+2*nsnon),
     + int_mb(lptr+1+4*nsnon),
     + int_mb(i_s2i1),dbl_mb(i_s3r1),dbl_mb(i_s1r1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+2)=lptr
      if(npairs.gt.0) then
      lptr=lptr+4*nsnon+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute-solute pairlist
c
      npairs=0
      mpairs=maxl-(lptr-i_list)-4*nsloc-1
      if(nsloc.gt.0.and.nsnon.gt.0) then
      call cf_lss(lself,lpbc,lpbcs,
     + xs,isga,isdt,isgr,isgm,isss,ishop,isz,isfr,isto,
     + jsfr,jsto,
     + nsloc,mpairs,npairs,int_mb(lptr+1),int_mb(lptr+1+2*nsloc),
     + int_mb(lptr+1+4*nsloc),
     + int_mb(i_s2i1),dbl_mb(i_s3r1),dbl_mb(i_s1r1),
     + int_mb(i_itrd(2)),mtt(2),int_mb(i_ixcl(2)),mxt(2),
     + int_mb(i_lda),dbl_mb(i_rda),lseq,int_mb(i_lsthop),
     + int_mb(i_mprot))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+3)=lptr
      if(npairs.gt.0) then
      lptr=lptr+4*nsloc+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute bond list
c
      npairs=0
      mpairs=maxl-(lptr-i_list)
      if(nsloc.gt.0.and.(lself.or.nsnon.gt.0)) then
      call cf_lsb(lself,isga,isdt,isz,isfr,isto,jsfr,jsto,
     + int_mb(i_ibnd(2)),mbt(2),npairs,mpairs,nconst,
     + int_mb(lptr+1),int_mb(i_s2i1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+4)=lptr
      if(npairs.gt.0) then
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute angle list
c
      npairs=0
      mpairs=maxl-(lptr-i_list)
      if(nsloc.gt.0.and.(lself.or.nsnon.gt.0)) then
      call cf_lsh(lself,isga,isdt,isz,isfr,isto,jsfr,jsto,
     + int_mb(i_iang(2)),mht(2),npairs,mpairs,nconst,
     + int_mb(lptr+1),int_mb(i_s2i1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+5)=lptr
      if(npairs.gt.0) then
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute torsion list
c
      npairs=0
      mpairs=maxl-(lptr-i_list)
      if(nsloc.gt.0.and.(lself.or.nsnon.gt.0)) then
      call cf_lsd(lself,isga,isdt,isz,isfr,isto,jsfr,jsto,
     + int_mb(i_idih(2)),mdt(2),npairs,mpairs,nconst,
     + int_mb(lptr+1),int_mb(i_s2i1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+6)=lptr
      if(npairs.gt.0) then
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute improper torsion list
c
      npairs=0
      mpairs=maxl-(lptr-i_list)
      if(nsloc.gt.0.and.(lself.or.nsnon.gt.0)) then
      call cf_lso(lself,isga,isdt,isz,isfr,isto,jsfr,jsto,
     + int_mb(i_iimp(2)),mit(2),npairs,mpairs,nconst,
     + int_mb(lptr+1),int_mb(i_s2i1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+7)=lptr
      if(npairs.gt.0) then
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
c     solute third neighbor list
c
      npairs=0
      mpairs=maxl-(lptr-i_list)
      if(nsloc.gt.0.and.(lself.or.nsnon.gt.0)) then
      call cf_lst(lself,isga,isdt,isz,isfr,isto,jsfr,jsto,
     + int_mb(i_itrd(2)),mtt(2),npairs,mpairs,
     + int_mb(lptr+1),int_mb(i_s2i1))
      endif
      int_mb(lptr)=npairs
      int_mb(lptrn+8)=lptr
      if(npairs.gt.0) then
      lptr=lptr+1+npairs
      else
      lptr=lptr+1
      endif
c
      llist=.true.
      ndxp=lptr-i_list
c
      return
      end



