      subroutine argos_space_trvl(xw,vw,xwcr,gw,iwl,iwlp,
     + xs,vs,gs,isl,islp,
     + boxsiz,ibownr,ipl,ndx,itmp,rtmp,lenx)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      real*8 gw(mwm,3,mwa),gs(msa,3)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      integer lenx
      real*8 boxsiz(maxbox,3)
      integer ipl(mbox,mip2),ibownr(maxbox,3)
      integer ndx(lenx),itmp(lenx)
      real*8 rtmp(lenx)
      logical lrec(27)
c
      integer i,indexw,indexs,j,k,ibx,iby,ibz,ipx,ipy,ipz
      integer isbox,isnod,nrbox
      integer ilp,ihp,jlp,jhp
      integer il,ih,jl,jh,ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer iliw,ihiw,jliw,jhiw
      integer ilis,ihis,jlis,jhis
      integer iwm,iwstay,jwstay,lwstay,nwgo,nwgosm
      integer nwgtsm
      integer nwnew,nwstay,iwmloc,jwmloc,lwmloc,irw
      integer isa,jsa,isstay,jsstay,icsgm,ifsgm,ilsgm
      integer nsnew,nsstay,lsstay,isaloc,jsaloc,lsaloc,irs
      integer nsgo,jnode,iwfr,iwto,isfr,isto
      real*8 factor,xscx,xscy,xscz,boxi(3)
      integer itemps,nfold
      logical lend
      character*255 string
c
      boxi(1)=one/box(1)
      boxi(2)=one/box(2)
      boxi(3)=one/box(3)
      nfold=0
      lpbc9=.false.
c
      nwstay=0
c
c     order the solvent molecules
c
      if(nwmloc.gt.0) then
      do 1 i=1,nwmloc
      ndx(i)=i
    1 continue
      endif
      if(nwmloc.gt.1) then
      lwmloc=nwmloc/2+1
      irw=nwmloc
    2 continue
      if(lwmloc.gt.1) then
      lwmloc=lwmloc-1
      itemps=ndx(lwmloc)
      else
      itemps=ndx(irw)
      ndx(irw)=ndx(1)
      irw=irw-1
      if(irw.eq.1) then
      ndx(1)=itemps
      goto 3
      endif
      endif
      iwmloc=lwmloc
      jwmloc=lwmloc+lwmloc
    4 continue
      if(jwmloc.le.irw) then
      if(jwmloc.lt.irw) then
      if((iwl(ndx(jwmloc),lwnod).eq.iwl(ndx(jwmloc+1),lwnod).and.
     + iwl(ndx(jwmloc),lwbox).le.iwl(ndx(jwmloc+1),lwbox)).or.
     + ((iwl(ndx(jwmloc),lwnod).eq.me.or.
     + (iwl(ndx(jwmloc),lwnod).ne.me.and.
     + iwl(ndx(jwmloc),lwnod).le.iwl(ndx(jwmloc+1),lwnod))).and.
     + iwl(ndx(jwmloc+1),lwnod).ne.me)) jwmloc=jwmloc+1
      endif
      if((iwl(itemps,lwnod).eq.iwl(ndx(jwmloc),lwnod).and.
     + iwl(itemps,lwbox).le.iwl(ndx(jwmloc),lwbox)).or.
     + ((iwl(itemps,lwnod).eq.me.or. (iwl(itemps,lwnod).ne.me.and.
     + iwl(itemps,lwnod).le.iwl(ndx(jwmloc),lwnod))).and.
     + iwl(ndx(jwmloc),lwnod).ne.me)) then
      ndx(iwmloc)=ndx(jwmloc)
      iwmloc=jwmloc
      jwmloc=jwmloc+jwmloc
      else
      jwmloc=irw+1
      endif
      goto 4
      endif
      ndx(iwmloc)=itemps
      goto 2
    3 continue
c
      do 5 k=1,3
      do 8 i=1,nwmloc
      rtmp(i)=xwcr(i,k)
    8 continue
      do 9 i=1,nwmloc
      xwcr(i,k)=rtmp(ndx(i))
    9 continue
      do 10 j=1,nwa
      do 11 i=1,nwmloc
      rtmp(i)=xw(i,k,j)
   11 continue
      do 12 i=1,nwmloc
      xw(i,k,j)=rtmp(ndx(i))
   12 continue
      do 13 i=1,nwmloc
      rtmp(i)=vw(i,k,j)
   13 continue
      do 14 i=1,nwmloc
      vw(i,k,j)=rtmp(ndx(i))
   14 continue
      if(iguide.gt.0) then
      do 113 i=1,nwmloc
      rtmp(i)=gw(i,k,j)
  113 continue
      do 114 i=1,nwmloc
      gw(i,k,j)=rtmp(ndx(i))
  114 continue
      endif
   10 continue
    5 continue
      do 18 k=1,miw2
      do 19 i=1,nwmloc
      itmp(i)=iwl(i,k)
   19 continue
      do 20 i=1,nwmloc
      iwl(i,k)=itmp(ndx(i))
   20 continue
   18 continue
      endif
c
      if(nwmloc.gt.0) then
      do 21 iwm=1,nwmloc
      if(iwl(iwm,lwnod).eq.me) then
      nwstay=iwm
      else
c
c     check if moving atoms go to neighboring processor
c
      do 222 k=1,27
      if(iwl(iwm,lwnod).eq.neighb(k,1)) goto 223
  222 continue
      write(string,'(a,i4,a,i4,a,i4,9f6.2)')
     +  'argos_space_travel: solvent molecule ',
     + iwl(iwm,lwgmn),' moving to non-neighbor ',iwl(iwm,lwnod),
     + ' from ',me,((xw(iwm,i,j),i=1,3),j=1,3)
      call md_abort(string,me)
  223 continue
c
      endif
c
c     testcode
c
      if(iand(idebug,8).eq.8) then
      if(iwl(iwm,lwnod).ne.me) write(lfndbg,'(a,3i5)')
     +  'Travel w fnd ',me,iwl(iwm,lwnod),iwl(iwm,lwgmn)
      endif
c
c     end test code
c
   21 continue
      endif
c
c     order the solute atoms
c
c     isl(isa,lsbox) : box
c     isl(isa,lsnod) : node
c     isl(isa,lssgm) : segment
c
      nsstay=0
      if(nsaloc.gt.0) then
      do 22 i=1,nsaloc
      ndx(i)=i
   22 continue
      endif
c
      if(nsaloc.gt.1) then
      lsaloc=nsaloc/2+1
      irs=nsaloc
   23 continue
      if(lsaloc.gt.1) then
      lsaloc=lsaloc-1
      itemps=ndx(lsaloc)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 24
      endif
      endif
      isaloc=lsaloc
      jsaloc=lsaloc+lsaloc
   25 continue
      if(jsaloc.le.irs) then
      if(jsaloc.lt.irs) then
      if((isl(ndx(jsaloc),lsnod).eq.isl(ndx(jsaloc+1),lsnod).and.
     + (isl(ndx(jsaloc),lsbox).lt.isl(ndx(jsaloc+1),lsbox).or.
     + (isl(ndx(jsaloc),lsbox).eq.isl(ndx(jsaloc+1),lsbox).and.
     + isl(ndx(jsaloc),lssgm).le.isl(ndx(jsaloc+1),lssgm)))).or.
     + ((isl(ndx(jsaloc),lsnod).eq.me.or.
     + (isl(ndx(jsaloc),lsnod).ne.me.and.
     + isl(ndx(jsaloc),lsnod).le.isl(ndx(jsaloc+1),lsnod))).and.
     + isl(ndx(jsaloc+1),lsnod).ne.me)) jsaloc=jsaloc+1
      endif
      if((isl(itemps,lsnod).eq.isl(ndx(jsaloc),lsnod).and.
     + (isl(itemps,lsbox).lt.isl(ndx(jsaloc),lsbox).or.
     + (isl(itemps,lsbox).eq.isl(ndx(jsaloc),lsbox).and.
     + isl(itemps,lssgm).le.isl(ndx(jsaloc),lssgm)))).or.
     + ((isl(itemps,lsnod).eq.me.or. (isl(itemps,lsnod).ne.me.and.
     + isl(itemps,lsnod).le.isl(ndx(jsaloc),lsnod))).and.
     + isl(ndx(jsaloc),lsnod).ne.me)) then
      ndx(isaloc)=ndx(jsaloc)
      isaloc=jsaloc
      jsaloc=jsaloc+jsaloc
      else
      jsaloc=irs+1
      endif
      goto 25
      endif
      ndx(isaloc)=itemps
      goto 23
   24 continue
c
      do 26 k=1,3
      do 27 i=1,nsaloc
      rtmp(i)=xs(i,k)
   27 continue
      do 28 i=1,nsaloc
      xs(i,k)=rtmp(ndx(i))
   28 continue
      do 29 i=1,nsaloc
      rtmp(i)=vs(i,k)
   29 continue
      do 30 i=1,nsaloc
      vs(i,k)=rtmp(ndx(i))
   30 continue
      if(iguide.gt.0) then
      do 2129 i=1,nsaloc
      rtmp(i)=gs(i,k)
 2129 continue
      do 2130 i=1,nsaloc
      gs(i,k)=rtmp(ndx(i))
 2130 continue
      endif
   26 continue
      do 40 k=1,mis2
      do 41 i=1,nsaloc
      itmp(i)=isl(i,k)
   41 continue
      do 42 i=1,nsaloc
      isl(i,k)=itmp(ndx(i))
   42 continue
   40 continue
      endif
c
      if(nsa.gt.0) then
      do 43 isa=1,nsaloc
      if(isl(isa,lsnod).eq.me) then
      nsstay=isa
      else
c
c     check if moving atoms go to neighboring processor
c
      do 444 k=1,27
      if(isl(isa,lsnod).eq.neighb(k,1)) goto 445
  444 continue
      write(string,'(a,i4,a,i4,a,i4,3f6.2)')
     +  'argos_space_travel: solute segment ',
     + isl(isa,lssgm),' moving to non-neighbor ',isl(isa,lsnod),
     + ' from ',me,(xs(isa,i),i=1,3)
      call md_abort(string,me)
  445 continue
c
      endif
   43 continue
      endif
c
c     make packages ready for shipment
c
c     loop over all neighboring nodes
c
      call ga_distribution(ga_iw,me,iliw,ihiw,jliw,jhiw)
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_distribution(ga_is,me,ilis,ihis,jlis,jhis)
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
c
      indexw=0
      indexs=0
      nwgosm=0
c
      do 70 i=1,27
      jnode=neighb(i,1)
      if(jnode.ge.0.and.jnode.ne.me) then
c
c     for the solvent
c
      iwfr=0
      iwto=0
      do 71 iwm=nwstay+1,nwmloc
      if(iwl(iwm,lwnod).eq.jnode) then
      if(iwfr.eq.0) iwfr=iwm
      iwto=iwm
c
c     testcode
c
      if(iand(idebug,8).eq.8) then
      if(iwl(iwm,lwnod).ne.me) write(lfndbg,'(a,3i5)')
     +  'Travel w snd ',me,iwl(iwm,lwnod),iwl(iwm,lwgmn)
      endif
c
c     end test code
c
      endif
   71 continue
c
c     if molecules need to travel copy coordinates etc into global array
c
      nwgo=iwto-iwfr+1
      if(iwfr.eq.0) nwgo=0
      ipl(1,1)=0
      ipl(1,2)=0
c
      if(nwgo.gt.0) then
      nwgosm=nwgosm+nwgo
      il=iliw+indexw
      ih=il+nwgo-1
      if(npackw.eq.0) then
      call ga_put(ga_iw,il,ih,jliw,jhiw,iwl(iwfr,1),mwm)
      else
      call argos_space_packw(ih-il+1,iwl(iwfr,1),iwlp(iwfr,1))
      call ga_put(ga_iw,il,ih,jliw,jliw+npackw-1,iwlp(iwfr,1),mwm)
      endif
      il=ilw+indexw
      ih=il+nwgo-1
      call ga_put(ga_w,il,ih,jlw,jlw+3*mwa-1,xw(iwfr,1,1),mwm)
      call ga_put(ga_w,il,ih,jlw+3*mwa,jlw+6*mwa-1,vw(iwfr,1,1),mwm)
      call ga_put(ga_w,il,ih,jlw+6*mwa,jlw+6*mwa+2,xwcr(iwfr,1),mwm)
      if(iguide.gt.0) then
      call ga_put(ga_w,il,ih,jlw+6*mwa+3,jlw+9*mwa+2,gw(iwfr,1,1),mwm)
      endif
      ipl(1,1)=indexw+1
      ipl(1,2)=indexw+nwgo
      indexw=indexw+nwgo
      endif
c
c     for the solute
c
      isfr=0
      isto=0
      do 72 isa=nsstay+1,nsaloc
      if(isl(isa,lsnod).eq.jnode) then
      if(isfr.eq.0) isfr=isa
      isto=isa
      endif
   72 continue
      nsgo=isto-isfr+1
      if(isfr.eq.0) nsgo=0
      ipl(1,3)=0
      ipl(1,4)=0
      if(nsgo.gt.0) then
      il=ilis+indexs
      ih=il+nsgo-1
      if(npack.eq.0) then
      call ga_put(ga_is,il,ih,jlis,jhis,isl(isfr,1),msa)
      else
      call argos_space_pack(ih-il+1,isl(isfr,1),islp(isfr,1))
      call ga_put(ga_is,il,ih,jlis,jlis+npack-1,islp(isfr,1),msa)
      endif
      call ga_put(ga_s,il,ih,jls,jls+2,xs(isfr,1),msa)
      call ga_put(ga_s,il,ih,jls+3,jls+5,vs(isfr,1),msa)
      if(iguide.gt.0) then
      call ga_put(ga_s,il,ih,jls+6,jls+8,gs(isfr,1),msa)
      endif
      ipl(1,3)=indexs+1
      ipl(1,4)=indexs+nsgo
      indexs=indexs+nsgo
      endif
c
c     inform other node of number of molecules to get
c
      if(ipl(1,1).gt.0.or.ipl(1,3).gt.0) then
      call ga_distribution(ga_ip,jnode,ilp,ihp,jlp,jhp)
      ilp=ilp+2+i
      call ga_put(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      endif
      endif
      lrec(i)=.false.
   70 continue
c
      call ga_sync()
c
c     receive molecules from other nodes
c
      nwgtsm=0
c
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp+30,jlp,jhp,ipl,mbox)
c
      do 74 i=1,27
      jnode=neighb(i,2)
      if(jnode.ge.0.and.jnode.ne.me.and..not.lrec(i)) then
c
      iwfr=ipl(3+i,1)
      iwto=ipl(3+i,2)
      isfr=ipl(3+i,3)
      isto=ipl(3+i,4)
c
      nwnew=iwto-iwfr+1
      nsnew=isto-isfr+1
c
      if(iwfr.eq.0) nwnew=0
      if(isfr.eq.0) nsnew=0
c
      if(nwstay+nwnew.gt.mwm) then
      write(string,'(a,i7,a,i7)')
     + 'Travel: mwm needs increase with ',nwnew,' to ',nwstay+nwnew
      call md_abort(string,me)
      endif
      if(nsstay+nsnew.gt.msa) then
      write(string,'(a,i7,a,i7)')
     + 'Travel: msa needs increase with ',nsnew,' to ',nsstay+nsnew
      call md_abort(string,me)
      endif
c
      lrec(i)=.true.
c
      if(iwfr.gt.0) then
      nwgtsm=nwgtsm+nwnew
      iwto=ipl(3+i,2)
      call ga_distribution(ga_iw,jnode,iliw,ihiw,jliw,jhiw)
      call ga_distribution(ga_w,jnode,ilw,ihw,jlw,jhw)
c
c     get data for additional molecules
c
      il=iliw+iwfr-1
      ih=iliw+iwto-1
      if(npackw.eq.0) then
      call ga_get(ga_iw,il,ih,jliw,jhiw,iwl(nwstay+1,1),mwm)
      else
      call ga_get(ga_iw,il,ih,jliw,jliw+npackw-1,iwlp(nwstay+1,1),mwm)
      call argos_space_unpackw(ih-il+1,iwl(nwstay+1,1),iwlp(nwstay+1,1))
      endif
      call ga_get(ga_w,il,ih,jlw,jlw+3*mwa-1,xw(nwstay+1,1,1),mwm)
      call ga_get(ga_w,il,ih,jlw+3*mwa,jlw+6*mwa-1,vw(nwstay+1,1,1),mwm)
      call ga_get(ga_w,il,ih,jlw+6*mwa,jlw+6*mwa+2,xwcr(nwstay+1,1),mwm)
      if(iguide.gt.0) then
      call ga_get(ga_w,il,ih,jlw+6*mwa+3,jlw+9*mwa+2,
     + gw(nwstay+1,1,1),mwm)
      endif
c
c     testcode
c
      if(iand(idebug,8).eq.8) then
      write(lfndbg,'(a,3i5)')
     +  ('Travel w rcv ',me,jnode,iwl(nwstay+k,lwgmn),k=1,nwnew)
      endif
c
c     end test code
c
c
c     update number of local solvent molecules
c
      nwstay=nwstay+nwnew
c
      endif
c
c     for the solute
c
      if(isfr.gt.0) then
      call ga_distribution(ga_is,jnode,ilis,ihis,jlis,jhis)
      call ga_distribution(ga_s,jnode,ils,ihs,jls,jhs)
      il=ilis+isfr-1
      ih=ilis+isto-1
      jl=jlis
      jh=jhis
      if(npack.eq.0) then
      call ga_get(ga_is,il,ih,jlis,jhis,isl(nsstay+1,1),msa)
      else
      call ga_get(ga_is,il,ih,jlis,jlis+npack-1,islp(nsstay+1,1),msa)
      call argos_space_unpack(ih-il+1,isl(nsstay+1,1),islp(nsstay+1,1))
      endif
      call ga_get(ga_s,il,ih,jls,jls+2,xs(nsstay+1,1),msa)
      call ga_get(ga_s,il,ih,jls+3,jls+5,vs(nsstay+1,1),msa)
      if(iguide.gt.0) then
      call ga_get(ga_s,il,ih,jls+6,jls+8,gs(nsstay+1,1),msa)
      endif
c
      nsstay=nsstay+nsnew
      endif
c
      endif
c
c     reset the pointers to zero
c
      ipl(3+i,1)=0
      ipl(3+i,2)=0
      ipl(3+i,3)=0
      ipl(3+i,4)=0
c
   74 continue
c
c     reset ipl in global array
c
      call ga_put(ga_ip,ilp,ilp+30,jlp,jhp,ipl,mbox)
c
c     order the solvent molecules according to subbox and
c     store indices into ip
c
c     ip(1,1)    : number of boxes on this node
c     ip(1,2)    : number of solvent molecules on this node
c     ip(2,2)    : number of solute atoms on this node
c
c     ip(3+i,1)  : index for solvents to be moved to the i-th neighbor
c
c     ip(30+i,1) : number of i-th box on this node
c     ip(30+i,2) : index to first solvent in i-th box
c     ip(30+i,3) : index to lasst solvent in i-th box
c
      if(nwstay.gt.0.and.(nwgosm.gt.0.or.nwgtsm.gt.0)) then
      do 81 i=1,nwstay
      ndx(i)=i
   81 continue
      if(nwstay.gt.1) then
      lwstay=nwstay/2+1
      irw=nwstay
   82 continue
      if(lwstay.gt.1) then
      lwstay=lwstay-1
      itemps=ndx(lwstay)
      else
      itemps=ndx(irw)
      ndx(irw)=ndx(1)
      irw=irw-1
      if(irw.eq.1) then
      ndx(1)=itemps
      goto 83
      endif
      endif
      iwstay=lwstay
      jwstay=lwstay+lwstay
   84 continue
      if(jwstay.le.irw) then
      if(jwstay.lt.irw) then
      if(iwl(ndx(jwstay),lwbox).le.iwl(ndx(jwstay+1),lwbox))
     + jwstay=jwstay+1
      endif
      if(iwl(itemps,lwbox).le.iwl(ndx(jwstay),lwbox)) then
      ndx(iwstay)=ndx(jwstay)
      iwstay=jwstay
      jwstay=jwstay+jwstay
      else
      jwstay=irw+1
      endif
      goto 84
      endif
      ndx(iwstay)=itemps
      goto 82
   83 continue
c
      do 85 k=1,3
      do 88 i=1,nwstay
      rtmp(i)=xwcr(i,k)
   88 continue
      do 89 i=1,nwstay
      xwcr(i,k)=rtmp(ndx(i))
   89 continue
      do 90 j=1,nwa
      do 91 i=1,nwstay
      rtmp(i)=xw(i,k,j)
   91 continue
      do 92 i=1,nwstay
      xw(i,k,j)=rtmp(ndx(i))
   92 continue
      do 93 i=1,nwstay
      rtmp(i)=vw(i,k,j)
   93 continue
      do 94 i=1,nwstay
      vw(i,k,j)=rtmp(ndx(i))
   94 continue
      if(iguide.gt.0) then
      do 193 i=1,nwstay
      rtmp(i)=gw(i,k,j)
  193 continue
      do 194 i=1,nwstay
      gw(i,k,j)=rtmp(ndx(i))
  194 continue
      endif
   90 continue
   85 continue
      do 98 k=1,miw2
      do 99 i=1,nwstay
      itmp(i)=iwl(i,k)
   99 continue
      do 100 i=1,nwstay
      iwl(i,k)=itmp(ndx(i))
  100 continue
   98 continue
c
      endif
      endif
c
c     order the solute according to segment
c
      if(nsstay.gt.0) then
      do 122 i=1,nsstay
      ndx(i)=i
  122 continue
      if(nsstay.gt.1) then
      lsstay=nsstay/2+1
      irs=nsstay
  123 continue
      if(lsstay.gt.1) then
      lsstay=lsstay-1
      itemps=ndx(lsstay)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 124
      endif
      endif
      isstay=lsstay
      jsstay=lsstay+lsstay
  125 continue
      if(jsstay.le.irs) then
      if(jsstay.lt.irs) then
      if(isl(ndx(jsstay),lssgm).le.isl(ndx(jsstay+1),lssgm))
     + jsstay=jsstay+1
      endif
      if(isl(itemps,lssgm).le.isl(ndx(jsstay),lssgm)) then
      ndx(isstay)=ndx(jsstay)
      isstay=jsstay
      jsstay=jsstay+jsstay
      else
      jsstay=irs+1
      endif
      goto 125
      endif
      ndx(isstay)=itemps
      goto 123
  124 continue
      endif
c
c     for each segment : 1. determine box number
c                        2. assign box number to each atom
c                        3. when box not owned by node:
c                           a. assign box number
c                           b. assign correct node number
c
      goto 666
      icsgm=isl(ndx(1),lssgm)
      ifsgm=1
      ilsgm=1
      do 126 isa=2,nsstay+1
c
c     if isa is first atom of a new segment or very last atom
c
 
      if(isa.le.nsstay) then
      lend=isl(ndx(isa),lssgm).ne.icsgm
      else
      lend=.true.
      endif
      if(lend) then
      if(isa.gt.nsstay) ilsgm=nsstay
      if(ifsgm.gt.0.and.ilsgm.ge.ifsgm) then
      xscx=zero
      xscy=zero
      xscz=zero
      do 127 jsa=ifsgm,ilsgm
      xscx=xscx+xs(ndx(jsa),1)
      xscy=xscy+xs(ndx(jsa),2)
      xscz=xscz+xs(ndx(jsa),3)
  127 continue
      factor=one/dble(ilsgm-ifsgm+1)
      xscx=factor*xscx
      xscy=factor*xscy
      xscz=factor*xscz
      if(npbtyp.ne.0) then
      if(abs(xscx).gt.boxh(1)) then
      xscx=xscx-nint(xscx*boxi(1))*box(1)
      nfold=1
      endif
      if(abs(xscy).gt.boxh(2)) then
      xscy=xscy-nint(xscy*boxi(2))*box(2)
      nfold=1
      endif
      if(abs(xscz).gt.boxh(3)) then
      xscz=xscz-nint(xscz*boxi(3))*box(3)
      nfold=1
      endif
      endif
c
c     determine the box number
c
      ibx=0
      iby=0
      ibz=0
      do 128 i=1,nbx-1
      if(xscx+boxh(1).gt.boxsiz(i,1)) ibx=i
  128 continue
      do 129 i=1,nby-1
      if(xscy+boxh(2).gt.boxsiz(i,2)) iby=i
  129 continue
      do 1130 i=1,nbz-1
      if(xscz+boxh(3).gt.boxsiz(i,3)) ibz=i
 1130 continue
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
      isbox=(ibz*nby+iby)*nbx+ibx
      isnod=(ipz*npy+ipy)*npx+ipx
c
c     assign box and node numbers
c
      do 1131 jsa=ifsgm,ilsgm
      isl(ndx(jsa),lsbox)=isbox
      isl(ndx(jsa),lsnod)=isnod
 1131 continue
c
      endif
      if(isa.le.nsstay) icsgm=isl(ndx(isa),lssgm)
      ifsgm=isa
      else
      ilsgm=isa
      endif
  126 continue
  666 continue
c
c     order solute according to box, segment, charge group, atom number
c
      if(nsstay.gt.1) then
      lsstay=nsstay/2+1
      irs=nsstay
  132 continue
      if(lsstay.gt.1) then
      lsstay=lsstay-1
      itemps=ndx(lsstay)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 133
      endif
      endif
      isstay=lsstay
      jsstay=lsstay+lsstay
  134 continue
      if(jsstay.le.irs) then
      if(jsstay.lt.irs) then
      if(isl(ndx(jsstay),lsbox).lt.isl(ndx(jsstay+1),lsbox).or.
     + (isl(ndx(jsstay),lsbox).eq.isl(ndx(jsstay+1),lsbox).and.
     + (isl(ndx(jsstay),lssgm).lt.isl(ndx(jsstay+1),lssgm).or.
     + (isl(ndx(jsstay),lssgm).eq.isl(ndx(jsstay+1),lssgm).and.
     + (isl(ndx(jsstay),lsgrp).lt.isl(ndx(jsstay+1),lsgrp).or.
     + (isl(ndx(jsstay),lsgrp).eq.isl(ndx(jsstay+1),lsgrp).and.
     + isl(ndx(jsstay),lsgan).le.isl(ndx(jsstay+1),lsgan)))))))
     + jsstay=jsstay+1
      endif
      if(isl(itemps,lsbox).lt.isl(ndx(jsstay),lsbox).or.
     + (isl(itemps,lsbox).eq.isl(ndx(jsstay),lsbox).and.
     + (isl(itemps,lssgm).lt.isl(ndx(jsstay),lssgm).or.
     + (isl(itemps,lssgm).eq.isl(ndx(jsstay),lssgm).and.
     + (isl(itemps,lsgrp).lt.isl(ndx(jsstay),lsgrp).or.
     + (isl(itemps,lsgrp).eq.isl(ndx(jsstay),lsgrp).and.
     + isl(itemps,lsgan).le.isl(ndx(jsstay),lsgan))))))) then
      ndx(isstay)=ndx(jsstay)
      isstay=jsstay
      jsstay=jsstay+jsstay
      else
      jsstay=irs+1
      endif
      goto 134
      endif
      ndx(isstay)=itemps
      goto 132
  133 continue
      endif
c
      do 135 k=1,3
      do 136 i=1,nsstay
      rtmp(i)=xs(i,k)
  136 continue
      do 137 i=1,nsstay
      xs(i,k)=rtmp(ndx(i))
  137 continue
      do 138 i=1,nsstay
      rtmp(i)=vs(i,k)
  138 continue
      do 139 i=1,nsstay
      vs(i,k)=rtmp(ndx(i))
  139 continue
      if(iguide.gt.0) then
      do 1138 i=1,nsstay
      rtmp(i)=gs(i,k)
 1138 continue
      do 1139 i=1,nsstay
      gs(i,k)=rtmp(ndx(i))
 1139 continue
      endif
  135 continue
      do 149 k=1,mis2
      do 150 i=1,nsstay
      itmp(i)=isl(i,k)
  150 continue
      do 151 i=1,nsstay
      isl(i,k)=itmp(ndx(i))
  151 continue
  149 continue
c
      endif
c
      do 200 i=1,ipl(1,1)
      ipl(30+i,2)=0
      ipl(30+i,3)=0
      ipl(30+i,4)=0
      ipl(30+i,5)=0
  200 continue
c
      do 201 i=1,ipl(1,1)
      nrbox=ipl(30+i,1)
      if(nwstay.gt.0) then
      do 202 iwm=1,nwstay
      if(iwl(iwm,lwbox).eq.nrbox) then
      if(ipl(30+i,2).eq.0) ipl(30+i,2)=iwm
      ipl(30+i,3)=iwm
      endif
  202 continue
      endif
      if(nsstay.gt.0) then
      do 203 isa=1,nsstay
      if(isl(isa,lsbox).eq.nrbox) then
      if(ipl(30+i,4).eq.0) ipl(30+i,4)=isa
      ipl(30+i,5)=isa
      endif
  203 continue
      endif
  201 continue
c
      nwmloc=nwstay
      ipl(1,2)=nwmloc
      nsaloc=nsstay
      ipl(2,2)=nsaloc
c
      call ga_igop(msp_23,nfold,1,'+')
      lpbc9=nfold.gt.0
c
      return
      end
c $Id$
