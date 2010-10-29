      logical function argos_prepare_charge(irtdb,lfnout,lfnpar,
     + source,jlo,ilo,ihi,jhi,
     + latm,catm,xatm,qatm,matm,natm,lbnd,mbnd,nbnd,maxscf,qscale)
c
c $Id$
c
      implicit none
c
#include "msgids.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "argos_prepare_common.fh"
c
      logical argos_prepare_task
      external argos_prepare_task
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
      integer irtdb,lfnout,lfnpar,jlo,ilo,ihi,jhi
      integer matm,natm,maxscf
      integer latm(5,matm)
      character*6 catm(3,matm)
      real*8 xatm(3,matm),qatm(matm),qscale
      integer mbnd,nbnd
      integer lbnd(2,mbnd)
c
      integer i,j,ia,ja,ka,iab(20),nb,nox,nhx,itask,maxit
      integer igeom2,nout,icon(100),ncon,nconst,irest,ihfree,idump,n
      real*8 rcon,resa,resb,toler
      character*16 tag
      real*8 coord(3),chrge,qsum,charge
      integer newgrd
      integer i_qb,l_qb
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
c     + call md_abort('argos_prepare_charge: rtdb_put nconst failed',0)
      goto 100
      endif
c
      if(ihi-ilo+1.gt.maxscf) goto 200
c
      if(util_print('sequence',print_default)) write(lfnout,1000) charge
      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,charge))
     + call md_abort('argos_prepare_charge: rtdb_put nconst failed',0)
      if(.not.geom_create(igeom2,'geometry'))
     + call md_abort('argos_prepare_charge: geom_create 1 failed',9999)
      if(.not.geom_set_user_units(igeom2,'nanometer'))
     + call md_abort('argos_prepare_charge: geom_set_user_units failed',
     + 9999)
      nout=0
      do 12 i=1,nbnd
      ia=ilo-1+lbnd(1,i)
      ja=ilo-1+lbnd(2,i)
      ka=0
      if((ia.lt.ilo.or.ia.gt.ihi).and.ja.ge.ilo.and.ja.le.ihi) ka=ja
      if((ja.lt.ilo.or.ja.gt.ihi).and.ia.ge.ilo.and.ia.le.ihi) ka=ia
      if(ka.gt.0) nout=nout+1
   12 continue
      nout=nout+ihi-ilo+1
      if(.not.geom_cart_set(igeom2,nout,tag,xatm,xatm))
     + call md_abort('argos_prepare_charge: geom_cart_set failed',9999)
      do 10 i=1,ihi-ilo+1
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
      if(.not.geom_cent_set(igeom2,i,tag,coord,chrge))
     + call md_abort('argos_prepare_charge: geom_cent_set failed',9999)
   10 continue
      nout=ihi-ilo+1
      do 11 i=1,nbnd
      ia=ilo-1+lbnd(1,i)
      ja=ilo-1+lbnd(2,i)
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
     + call md_abort('argos_prepare_charge: geom_cent_set failed',9999)
      endif
   11 continue
      if(.not.geom_rtdb_store(irtdb,igeom2,'geometry'))
     + call md_abort('argos_prepare_charge: geom_rtdb_store failed',
     + 9999)
      if(.not.geom_destroy(igeom2))
     + call md_abort('argos_prepare_charge: geom_destroy failed',9999)
c
      if(.not.rtdb_put(irtdb,'esp:newgrd',mt_int,1,newgrd))
     + call md_abort('argos_prepare_charge: rtdb_put newgrd failed',0)
      if(rtdb_get(irtdb,'esp:nconst',mt_int,1,nconst)) goto 101
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
     + call md_abort('argos_prepare_charge: rtdb_put nconst failed',0)
      if(nconst.gt.0) then
      if(.not.rtdb_put(irtdb,'esp:ncon',mt_int,1,ncon))
     + call md_abort('argos_prepare_charge: rtdb_put ncon failed',0)
      if(.not.rtdb_put(irtdb,'esp:icon',mt_int,ncon,icon))
     + call md_abort('argos_prepare_charge: rtdb_put icon failed',0)
      if(.not.rtdb_put(irtdb,'esp:rcon',mt_dbl,nconst,rcon))
     + call md_abort('argos_prepare_charge: rtdb_put rcon failed',0)
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
     + call md_abort('argos_prepare_charge: rtdb_put irest failed',0)
      if(.not.rtdb_put(irtdb,'esp:resa',mt_dbl,1,resa))
     + call md_abort('argos_prepare_charge: rtdb_put resa failed',0)
      if(.not.rtdb_put(irtdb,'esp:resb',mt_dbl,1,resb))
     + call md_abort('argos_prepare_charge: rtdb_put resb failed',0)
      if(.not.rtdb_put(irtdb,'esp:hfree',mt_int,1,ihfree))
     + call md_abort('argos_prepare_charge: rtdb_put hfree failed',0)
      if(.not.rtdb_put(irtdb,'esp:maxit',mt_int,1,maxit))
     + call md_abort('argos_prepare_charge: rtdb_put maxit failed',0)
      if(.not.rtdb_put(irtdb,'esp:toler',mt_dbl,1,toler))
     + call md_abort('argos_prepare_charge: rtdb_put toler failed',0)
      if(.not.rtdb_put(irtdb,'esp:dump',mt_int,1,idump))
     + call md_abort('argos_prepare_charge: rtdb_put dump failed',0)
  101 continue
c
      itask=1
      if(.not.argos_prepare_task(irtdb,itask,lfnout))
     + call md_abort('argos_prepare_task failed',9999)
c
      if(.not.geom_create(igeom2,'pre:geometry'))
     + call md_abort('argos_prepare_charge: geom_create 2 failed',0)
c
      if(rtdb_get(irtdb,'geometry:esp_fit:geometry:ncenter',mt_int,1,n))
     + then
      if(.not.geom_rtdb_load(irtdb,igeom2,'esp_fit:geometry'))
     + call md_abort('argos_prepare_charge: geom_rtdb_load failed',0)
      else
      call md_abort('argos_prepare_charge: esp_fit rtdb geometry',0)
      endif
c
      qsum=0.0d0
      do 9 i=1,ihi-ilo+1
      if(.not.geom_cent_get(igeom2,i,tag,coord,chrge))
     + call md_abort('argos_prepare_charge: geom_cent_get failed',i)
      qatm(i)=chrge
      qsum=qsum+chrge
    9 continue
      if(.not.geom_destroy(igeom2))
     + call md_abort('argos_prepare_charge: geom_destroy failed',9999)
      if(abs(qsum).lt.1.0d-5) then
      do 14 i=1,ihi-ilo+1
      qatm(i)=qscale*qatm(i)
   14 continue
      endif
      argos_prepare_charge=.true.
      return
c
  200 continue
c
      if(.not.ma_push_get(mt_dbl,mbnd,'qb',l_qb,i_qb))
     + call md_abort('Failed to allocate qb',0)
c
      call argos_prepare_quick_charge(lfnout,lfnpar,
     + ilo,ihi,
     + latm,catm,qatm,matm,natm,
     + lbnd,dbl_mb(i_qb),mbnd,nbnd)
c
      if(.not.ma_pop_stack(l_qb))
     + call md_abort('Failed to deallocate qb',0)
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
     + qatm(iab(j))=qatm(iab(j))-0.33
    6 continue
      qatm(i)=qatm(i)-0.34
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
     + qatm(iab(j))=qatm(iab(j))+0.25
    8 continue
      qatm(i)=qatm(i)+0.25
      endif
      endif
c
c     -O-POO- (-)
c
      if(latm(2,i).eq.15.and.nb.eq.4) then
      nox=0
      do 20 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1) nox=nox+1
   20 continue
      if(nox.eq.2) then
      do 21 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))-0.33
   21 continue
      qatm(i)=qatm(i)-0.34
      endif
      if(nox.eq.3) then
      do 22 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))-0.50
   22 continue
      qatm(i)=qatm(i)-0.50
      endif
      endif
c
c     --SOOO (-)
c
      if(latm(2,i).eq.16.and.nb.eq.4) then
      nox=0
      do 23 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1) nox=nox+1
   23 continue
      if(nox.eq.3) then
      do 24 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))-0.25
   24 continue
      qatm(i)=qatm(i)-0.25
      endif
      endif
c
    3 continue
c
      argos_prepare_charge=.true.
      return
c
 9999 argos_prepare_charge=.false.
      return
      end
