      logical function pre_charge(irtdb,lfnout,source,jlo,ilo,ihi,jhi,
     + latm,catm,xatm,qatm,matm,natm,lbnd,mbnd,nbnd,maxscf,qscale)
c
c $Id: pre_charge.F,v 1.21 1999-07-09 14:46:03 d3j191 Exp $
c
      implicit none
c
#include "msgids.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
c
      logical pre_task
      external pre_task
c
c     conversion factor
c
c     cnm2au : from nm to bohr
c
      real*8 cnm2au
      parameter(cnm2au=1.88972598786d+01)
c
c     in  : latm(2,matm)      = atomic number
c           matm              = dimension atom list
c           lbnd(1:2,mbnd)    = bond indices
c           mbnd              = dimension bond list
c           nbnd              = length bond list
c
c     out : qatm(matm)        = guestimated atomic partial charges
c
      character*80 source
      integer irtdb,lfnout,jlo,ilo,ihi,jhi
      integer matm,natm,maxscf
      integer latm(5,matm)
      character*6 catm(3,matm)
      real*8 xatm(3,matm),qatm(matm),qscale
      integer mbnd,nbnd
      integer lbnd(2,mbnd)
c
      integer i,j,ia,ja,ka,it,jt,itmp,iab(20),nb,nox,nhx,itask,maxit
      integer igeom2,nout,icon(100),ncon,nconst,irest,ihfree,idump,n
      real*8 rcon,resa,resb,toler
      character*16 tag
      real*8 coord(3),chrge,qsum,charge
      integer newgrd
c
      newgrd=1
c
      do 1 i=1,matm
      qatm(i)=0.0d0
    1 continue
c
      charge=0.0d0
      do 15 i=ilo,ihi
      nb=0
      do 16 j=1,nbnd
      ia=lbnd(1,j)
      ja=lbnd(2,j)
      if(ia.eq.i) then
      nb=nb+1
      iab(nb)=ja
      else
      if(ja.eq.i) then
      nb=nb+1
      iab(nb)=ia
      endif
      endif
   16 continue
c
c     C with 3 bonds
c
      if(latm(2,i).eq.6.and.nb.eq.3) then
      nox=0
      do 17 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1) nox=nox+1
   17 continue
      if(nox.eq.2) charge=charge-1.0d0
      endif
c
c     N with 4 bonds
c
      if(latm(2,i).eq.7.and.nb.eq.4) then
      nhx=0
      do 18 j=1,nb
      if(latm(2,iab(j)).eq.1.and.latm(3,iab(j)).eq.1) nhx=nhx+1
   18 continue
      if(nhx.eq.3) charge=charge+1.0d0
      endif
c
c     P with 4 bonds
c
      if(latm(2,i).eq.15.and.nb.eq.4) then
      nox=0
      do 19 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1) nox=nox+1
   19 continue
      if(nox.eq.2) charge=charge-1.0d0
      if(nox.eq.3) charge=charge-2.0d0
      if(nox.eq.4) charge=charge-3.0d0
      endif
c
c
c     other
c
      if(latm(2,i).eq.3.and.nb.eq.0) charge=charge+1.0d0
      if(latm(2,i).eq.11.and.nb.eq.0) charge=charge+1.0d0
      if(latm(2,i).eq.19.and.nb.eq.0) charge=charge+1.0d0
      if(latm(2,i).eq.20.and.nb.eq.0) charge=charge+2.0d0
      if(latm(2,i).eq.9.and.nb.eq.0) charge=charge-1.0d0
      if(latm(2,i).eq.17.and.nb.eq.0) charge=charge-1.0d0
      if(latm(2,i).eq.35.and.nb.eq.0) charge=charge-1.0d0
      if(latm(2,i).eq.53.and.nb.eq.0) charge=charge-1.0d0
   15 continue
c
      if(source(1:8).eq.'geometry') then
      nout=natm
      if(util_print('sequence',print_default)) write(lfnout,1000) charge
 1000 format(' Fragment charge',t40,f12.6,/)
c      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,charge))
c     + call errquit('pre_charge: rtdb_put nconst failed',0)
      goto 100
      endif
c
      if(ihi-ilo+1.gt.maxscf) goto 200
c
      if(util_print('sequence',print_default)) write(lfnout,1000) charge
      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,charge))
     + call errquit('pre_charge: rtdb_put nconst failed',0)
      if(.not.geom_create(igeom2,'geometry'))
     + call errquit('pre_charge: geom_create 1 failed',9999)
      if(.not.geom_set_user_units(igeom2,'nanometer'))
     + call errquit('pre_charge: geom_set_user_units failed',9999)
      nout=0
      do 12 i=1,nbnd
      ia=lbnd(1,i)
      ja=lbnd(2,i)
      ka=0
      if((ia.lt.ilo.or.ia.gt.ihi).and.ja.ge.ilo.and.ja.le.ihi) ka=ja
      if((ja.lt.ilo.or.ja.gt.ihi).and.ia.ge.ilo.and.ia.le.ihi) ka=ia
      if(ka.gt.0) nout=nout+1
   12 continue
      nout=nout+ihi-ilo+1
      if(.not.geom_cart_set(igeom2,nout,tag,xatm,xatm))
     + call errquit('pre_charge: geom_cart_set failed',9999)
      do 10 i=ilo,ihi
      tag='                '
      if(latm(2,i).eq.1) then
      tag(1:1)='H'
      else
      if(catm(1,i)(1:1).eq.' ') then
      tag(1:1)=catm(1,i)(2:2)
      else
      tag(1:2)=catm(1,i)(1:2)
      endif
      endif
      chrge=dble(latm(2,i))
      coord(1)=cnm2au*xatm(1,i)
      coord(2)=cnm2au*xatm(2,i)
      coord(3)=cnm2au*xatm(3,i)
      if(.not.geom_cent_set(igeom2,i-ilo+1,tag,coord,chrge))
     + call errquit('pre_charge: geom_cent_set failed',9999)
   10 continue
      nout=ihi-ilo+1
      do 11 i=1,nbnd
      ia=lbnd(1,i)
      ja=lbnd(2,i)
      ka=0
      if((ia.lt.ilo.or.ia.gt.ihi).and.ja.ge.ilo.and.ja.le.ihi) ka=ia
      if((ja.lt.ilo.or.ja.gt.ihi).and.ia.ge.ilo.and.ia.le.ihi) ka=ja
      if(ka.gt.0) then
      nout=nout+1
      tag='H               '
      chrge=1.0d0
      coord(1)=cnm2au*xatm(1,ka)
      coord(2)=cnm2au*xatm(2,ka)
      coord(3)=cnm2au*xatm(3,ka)
      if(.not.geom_cent_set(igeom2,nout,tag,coord,chrge))
     + call errquit('pre_charge: geom_cent_set failed',9999)      
      endif
   11 continue
      if(.not.geom_rtdb_store(irtdb,igeom2,'geometry'))
     + call errquit('pre_charge: geom_rtdb_store failed',9999)
      if(.not.geom_destroy(igeom2))
     + call errquit('pre_charge: geom_destroy failed',9999)
c
      if(.not.rtdb_put(irtdb,'esp:newgrd',mt_int,1,newgrd))
     + call errquit('pre_charge: rtdb_put newgrd failed',0)
      nconst=0
      if(nout.gt.ihi-ilo+1) then
      ncon=ihi-ilo+1
      icon(1)=nout-ncon
      do 13 i=ncon+1,nout
      icon(i-ncon+1)=i
   13 continue
      ncon=nout-ncon+1
      nconst=1
      rcon=0.0d0
      endif
      if(.not.rtdb_put(irtdb,'esp:nconst',mt_int,1,nconst))
     + call errquit('pre_charge: rtdb_put nconst failed',0)
      if(nconst.gt.0) then
      if(.not.rtdb_put(irtdb,'esp:ncon',mt_int,1,ncon))
     + call errquit('pre_charge: rtdb_put ncon failed',0)
      if(.not.rtdb_put(irtdb,'esp:icon',mt_int,ncon,icon))
     + call errquit('pre_charge: rtdb_put icon failed',0)
      if(.not.rtdb_put(irtdb,'esp:rcon',mt_dbl,nconst,rcon))
     + call errquit('pre_charge: rtdb_put rcon failed',0)
      endif
  100 continue
      irest=2
      resa=0.001d0
      resb=0.1d0
      ihfree=0
      maxit=25
      toler=1.0d-04
      idump=0
      if(.not.rtdb_put(irtdb,'esp:irest',mt_int,1,irest))
     + call errquit('pre_charge: rtdb_put irest failed',0)
      if(.not.rtdb_put(irtdb,'esp:resa',mt_dbl,1,resa))
     + call errquit('pre_charge: rtdb_put resa failed',0)
      if(.not.rtdb_put(irtdb,'esp:resb',mt_dbl,1,resb))
     + call errquit('pre_charge: rtdb_put resb failed',0)
      if(.not.rtdb_put(irtdb,'esp:hfree',mt_int,1,ihfree))
     + call errquit('pre_charge: rtdb_put hfree failed',0)
      if(.not.rtdb_put(irtdb,'esp:maxit',mt_int,1,maxit))
     + call errquit('pre_charge: rtdb_put maxit failed',0)
      if(.not.rtdb_put(irtdb,'esp:toler',mt_dbl,1,toler))
     + call errquit('pre_charge: rtdb_put toler failed',0)
      if(.not.rtdb_put(irtdb,'esp:dump',mt_int,1,idump))
     + call errquit('pre_charge: rtdb_put dump failed',0)
c
      itask=1
      if(.not.pre_task(irtdb,itask,lfnout))
     + call errquit('pre_task failed',9999)
c
      if(.not.geom_create(igeom2,'pre:geometry'))
     + call errquit('pre_charge: geom_create 2 failed',0)
c
      if(rtdb_get(irtdb,'geometry:esp_fit:geometry:ncenter',mt_int,1,n))
     + then
      if(.not.geom_rtdb_load(irtdb,igeom2,'esp_fit:geometry'))
     + call errquit('pre_charge: geom_rtdb_load failed',0) 
      else
      call errquit('pre_charge: esp_fit geometry not found on rtdb',0)
      endif
c
      qsum=0.0d0
      do 9 i=1,ihi-ilo+1
      if(.not.geom_cent_get(igeom2,i,tag,coord,chrge))
     + call errquit('pre_charge: geom_cent_get failed',i)
      qatm(ilo+i-1)=chrge
      qsum=qsum+chrge
    9 continue
      if(.not.geom_destroy(igeom2))
     + call errquit('pre_charge: geom_destroy failed',9999)
      if(abs(qsum).lt.1.0d-5) then
      do 14 i=1,ihi-ilo+1
      qatm(ilo+i-1)=qscale*qatm(ilo+i-1)
   14 continue
      endif
      pre_charge=.true.
      return
c
  200 continue
c
      do 2 i=1,nbnd
      ia=lbnd(1,i)
      ja=lbnd(2,i)
      it=latm(2,ia)
      jt=latm(2,ja)
      if(it.ne.jt) then
c
c     smaller atom number in it/ia
c
      if(it.gt.jt) then
      itmp=ia
      ia=ja
      ja=itmp
      itmp=it
      it=jt
      jt=itmp
      endif
c
c     H - C : 0.05  : -0.05
c
      if(it.eq.1.and.jt.eq.6) then
      qatm(ia)=qatm(ia)+0.05d0
      qatm(ja)=qatm(ja)-0.05d0
      endif
c
c     H - N : 0.27  : -0.27
c
      if(it.eq.1.and.jt.eq.7) then
      qatm(ia)=qatm(ia)+0.27d0
      qatm(ja)=qatm(ja)-0.27d0
      endif
c
c     H - O : 0.19  : -0.19
c
      if(it.eq.1.and.jt.eq.8) then
      qatm(ia)=qatm(ia)+0.40d0
      qatm(ja)=qatm(ja)-0.40d0
      endif
c
c     H - S : 0.19  : -0.19
c
      if(it.eq.1.and.jt.eq.16) then
      qatm(ia)=qatm(ia)+0.19d0
      qatm(ja)=qatm(ja)-0.19d0
      endif
c
c     C - N : 0.07  : -0.07
c
      if(it.eq.6.and.jt.eq.7) then
      qatm(ia)=qatm(ia)+0.07d0
      qatm(ja)=qatm(ja)-0.07d0
      endif
c
c     C - O : 0.30  : -0.30
c
      if(it.eq.6.and.jt.eq.8.and.latm(3,ja).eq.2) then
      qatm(ia)=qatm(ia)+0.30d0
      qatm(ja)=qatm(ja)-0.30d0
      endif
c
c     C = O : 0.57  : -0.57
c
      if(it.eq.6.and.jt.eq.8.and.latm(3,ja).eq.1) then
      qatm(ia)=qatm(ia)+0.57d0
      qatm(ja)=qatm(ja)-0.57d0
      endif
c
c     C - S : 0.11  : -0.11
c
      if(it.eq.6.and.jt.eq.16) then
      qatm(ia)=qatm(ia)+0.11d0
      qatm(ja)=qatm(ja)-0.11d0
      endif
c
      endif
    2 continue
c
      do 3 i=1,natm
      nb=0
      do 4 j=1,nbnd
      ia=lbnd(1,j)
      ja=lbnd(2,j)
      if(ia.eq.i) then
      nb=nb+1
      iab(nb)=ja
      else
      if(ja.eq.i) then
      nb=nb+1
      iab(nb)=ia
      endif
      endif
    4 continue
c
c     -COO (-)
c
      if(latm(2,i).eq.6.and.nb.eq.3) then
      nox=0
      do 5 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1) nox=nox+1
    5 continue
      if(nox.eq.2) then
      do 6 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))-0.23
    6 continue
      qatm(i)=qatm(i)-0.54
      endif
      endif
c
c     -NH3 (+)
c
      if(latm(2,i).eq.7.and.nb.eq.4) then
      nhx=0
      do 7 j=1,nb
      if(latm(2,iab(j)).eq.1.and.latm(3,iab(j)).eq.1) nhx=nhx+1
    7 continue
      if(nhx.eq.3) then
      do 8 j=1,nb
      if(latm(2,iab(j)).eq.1.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))+0.07
    8 continue
      qatm(i)=qatm(i)+0.79
      endif
      endif
c
    3 continue
c
      pre_charge=.true.
      return
c
 9999 pre_charge=.false.
      return
      end
