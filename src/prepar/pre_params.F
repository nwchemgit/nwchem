      logical function pre_params(lfnpar,lfnout,nparms,mparms,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + ffield,releps,q14fac,ignore,
     + latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,limp,kimp,rimp,mimp,nimp,
     + latt,lats,catt,patt,ratt,matt,natt,mats,nats,
     + lseq,cseq,mseq,nseq)
c
c $Id: pre_params.F,v 1.16 2003-07-17 19:38:57 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
      logical pre_ffield,pre_parcnv,pre_check,pre_dummy,pre_multip
      external pre_ffield,pre_parcnv,pre_check,pre_dummy,pre_multip
c
      integer lfnpar,lfnout,ignore,nparms,mparms
      character*255 dir_s,dir_x,dir_q,dir_u,dir_t,dir_c
      character*255 par_s,par_x,par_q,par_u,par_t,par_c
      character*80 ffield
      integer matm,natm
      integer latm(11,matm)
      character*6 catm(mparms,matm)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(nparms,2,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(nparms,2,mang)
      integer mdih,ndih
      integer ldih(6,mdih),nfdih(nparms,mdih),kdih(6,nparms,mdih)
      real*8 rdih(6,nparms,2,mdih)
      integer mimp,nimp
      integer limp(6,mimp),kimp(nparms,mimp)
      real*8 rimp(nparms,2,mimp)
      integer matt,natt,mats,nats
      integer latt(matt),lats(nparms,mats)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer mseq,nseq
      integer lseq(6,mseq)
      character*10 cseq(mseq)
c
      integer i,j,k,kk,lp,ld,icomb
      character*255 filnam
      real*8 releps,q14fac
c
c     make atom type list
c     -------------------
c
c     natt : number of atom types
c
      natt=0
      do 1 i=1,natm
      do 2 j=1,nparms
      do 3 k=1,natt
      if(catm(1+j,i).eq.catt(1,k)) goto 2
    3 continue
      natt=natt+1
      if(natt.gt.matt) call errquit('increase matt',9999)
      catt(1,natt)=catm(1+j,i)
      latt(natt)=0
    2 continue
    1 continue
c
      nats=0
      do 5 i=1,natm
      do 4 j=1,nats
      kk=0
      do 44 k=1,nparms
      if(catm(k+1,i).eq.catt(1,lats(k,j))) kk=kk+1
   44 continue
      if(kk.eq.nparms) then
      latm(3,i)=j
      goto 5
      endif
    4 continue
      nats=nats+1
      if(nats.gt.mats) call errquit('increase mats',9999)
      do 6 k=1,natt
      do 66 kk=1,nparms
      if(catm(kk+1,i).eq.catt(1,k)) lats(kk,nats)=k
   66 continue
    6 continue
      latm(3,i)=nats
    5 continue
c
c     substitute parameters parameter files
c     -------------------------------------
c
      icomb=0
      lp=index(par_s,' ')
      ld=index(dir_s,' ')-1
      if(ld.gt.0) then
      filnam=dir_s(1:ld)//par_s(1:lp)
      if(.not.pre_ffield(1,lfnpar,filnam,lfnout,q14fac,releps,
     + ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)) then
      endif
      endif
      lp=index(par_x,' ')
      ld=index(dir_x,' ')-1
      if(ld.gt.0) then
      filnam=dir_x(1:ld)//par_x(1:lp)
      if(.not.pre_ffield(2,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)) then
      endif
      endif
      lp=index(par_q,' ')
      ld=index(dir_q,' ')-1
      if(ld.gt.0) then
      filnam=dir_q(1:ld)//par_q(1:lp)
      if(.not.pre_ffield(2,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)) then
      endif
      endif
      lp=index(par_u,' ')
      ld=index(dir_u,' ')-1
      if(ld.gt.0) then
      filnam=dir_u(1:ld)//par_u(1:lp)
      if(.not.pre_ffield(3,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)) then
      endif
      endif
      lp=index(par_t,' ')
      ld=index(dir_t,' ')-1
      if(ld.gt.0) then
      filnam=dir_t(1:ld)//par_t(1:lp)
      if(.not.pre_ffield(4,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)) then
      endif
      endif
      lp=index(par_c,' ')
      ld=index(dir_c,' ')-1
      if(ld.gt.0) then
      filnam=dir_c(1:ld)//par_c(1:lp)
      if(.not.pre_ffield(4,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)) then
      endif
      endif
c
      if(util_print('topology',print_default)) then
      write(lfnout,1000)
 1000 format(' ')
      endif
c
c     convert non-bonded parameters to C6 and C12
c     -------------------------------------------
c
      if(.not.pre_parcnv(icomb,latt,patt,catt,matt,natt))
     + call errquit('pre_parcnv failed',9999)
c
c     copy bonded parameters for dummy atoms
c     --------------------------------------
c
      if(.not.pre_dummy(lfnout,latt,catt,patt,ratt,matt,natt,
     + latm,catm,matm,natm,lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms))
     + call errquit('pre_dummy failed',9999)
c
c     check if all required parameters have been found
c     ------------------------------------------------
c
      if(.not.pre_check(lfnout,ffield,nparms,mparms,ignore,
     + latt,catt,matt,natt,catm,latm,matm,natm,lbnd,mbnd,nbnd,
     + lang,mang,nang,ldih,mdih,ndih,limp,mimp,nimp,
     + lseq,cseq,mseq,nseq))
     + call errquit('pre_check failed',9999)
c
c     check for angles, dihedrals and impropers that involve more than 2 processors
c
      if(.not.pre_multip(lseq,mseq,nseq,latm,matm,natm,lang,mang,nang,
     + ldih,mdih,ndih,limp,mimp,nimp))
     + call errquit('pre_multip failed',9999)
c
      pre_params=.true.
      return
      end

