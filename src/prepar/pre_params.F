      logical function pre_params(lfnpar,lfnout,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + ffield,releps,q14fac,ignore,
     + latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,kdih,rdih,mdih,ndih,limp,rimp,mimp,nimp,
     + latt,lats,catt,patt,ratt,matt,natt,mats,nats)
c
c $Id: pre_params.F,v 1.15 2001-11-06 06:15:21 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
      logical pre_ffield,pre_parcnv,pre_check,pre_dummy
      external pre_ffield,pre_parcnv,pre_check,pre_dummy
c
      integer lfnpar,lfnout,ignore
      character*255 dir_s,dir_x,dir_q,dir_u,dir_t,dir_c
      character*255 par_s,par_x,par_q,par_u,par_t,par_c
      character*80 ffield
      integer matm,natm
      integer latm(9,matm)
      character*6 catm(4,matm)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(6,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(6,mang)
      integer mdih,ndih
      integer ldih(9,mdih),kdih(6,3,mdih)
      real*8 rdih(6,6,mdih)
      integer mimp,nimp
      integer limp(9,mimp)
      real*8 rimp(6,mimp)
      integer matt,natt,mats,nats
      integer latt(matt),lats(3,mats)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
c
      integer i,j,k,lp,ld,icomb
      character*255 filnam
      real*8 releps,q14fac
c
c     make atom type list
c     -------------------
c
      natt=0
      do 1 i=1,natm
      do 2 j=1,3
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
      if(catm(2,i).eq.catt(1,lats(1,j)).and.
     + catm(3,i).eq.catt(1,lats(2,j)).and.
     + catm(4,i).eq.catt(1,lats(3,j))) then
      latm(3,i)=j
      goto 5
      endif
    4 continue
      nats=nats+1
      if(nats.gt.mats) call errquit('increase mats',9999)
      do 6 k=1,natt
      if(catm(2,i).eq.catt(1,k)) lats(1,nats)=k
      if(catm(3,i).eq.catt(1,k)) lats(2,nats)=k
      if(catm(4,i).eq.catt(1,k)) lats(3,nats)=k
    6 continue
      latm(3,i)=nats
    5 continue
c
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
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)) then
      endif
      endif
      lp=index(par_x,' ')
      ld=index(dir_x,' ')-1
      if(ld.gt.0) then
      filnam=dir_x(1:ld)//par_x(1:lp)
      if(.not.pre_ffield(2,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)) then
      endif
      endif
      lp=index(par_q,' ')
      ld=index(dir_q,' ')-1
      if(ld.gt.0) then
      filnam=dir_q(1:ld)//par_q(1:lp)
      if(.not.pre_ffield(2,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)) then
      endif
      endif
      lp=index(par_u,' ')
      ld=index(dir_u,' ')-1
      if(ld.gt.0) then
      filnam=dir_u(1:ld)//par_u(1:lp)
      if(.not.pre_ffield(3,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)) then
      endif
      endif
      lp=index(par_t,' ')
      ld=index(dir_t,' ')-1
      if(ld.gt.0) then
      filnam=dir_t(1:ld)//par_t(1:lp)
      if(.not.pre_ffield(4,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)) then
      endif
      endif
      lp=index(par_c,' ')
      ld=index(dir_c,' ')-1
      if(ld.gt.0) then
      filnam=dir_c(1:ld)//par_c(1:lp)
      if(.not.pre_ffield(4,lfnpar,filnam,lfnout,q14fac,releps,
     +  ffield,icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)) then
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
     + ldih,kdih,rdih,mdih,ndih,limp,rimp,mimp,nimp))
     + call errquit('pre_dummy failed',9999)
c
c     check if all required parameters have been found
c     ------------------------------------------------
c
      if(.not.pre_check(lfnout,ffield,ignore,
     + latt,catt,matt,natt,catm,matm,natm,lbnd,mbnd,nbnd,
     + lang,mang,nang,ldih,mdih,ndih,limp,mimp,nimp))
     + call errquit('pre_check failed',9999)
c
      pre_params=.true.
      return
      end
