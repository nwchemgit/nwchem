      subroutine hnd_giaspd(gint,gijkl,gnkl,gnm,b00,b01,b10,c00,d00,
     1                      f00,dij,dkl,abv,cv,rwv)
c
c $Id: hnd_giaspd.F,v 1.1 2004-07-21 16:24:35 bert Exp $
c
      implicit none
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
#include "hnd_rys.fh"
c
c     Input parameters
c
      double precision gint(*),gijkl(*),gnkl(*),gnm(*)
      double precision b00(*),b01(*),b10(*),c00(*),d00(*),f00(*)
      double precision dij(*),dkl(*)
      double precision abv(5,*),cv(18,*),rwv(2,*)
c
c     Local variabels
c
      double precision pi252
      data pi252   /34.986836655250d+00/
      integer nimax,njmax,nkmax,nlmax,nmax,mmax,ng,maxg
      logical last,testy
      integer klg,ijg,numg,i,j,k
      double precision dtol2,q4,dtol
      double precision db,bb,xb,yb,zb,xd,yd,zd,dxkl,dykl,dzkl,q4db
      double precision da,aa,xa,ya,za,aandb1,q4dbda,dum,rho,xx
      double precision xc,yc,zc,dxij,dyij,dzij
c
      numi = lit*(lit+1)/2
      numj = ljt*(ljt+1)/2
      numk = lkt*(lkt+1)/2
      numl = lmt*(lmt+1)/2
      nimax=lit+1
      njmax=ljt
      nkmax=lkt+1
      nlmax=lmt
      nmax =lit+ljt
      mmax =lkt+lmt
      nmaxs=nmax.eq.1
      nmaxp=nmax.le.2
      mmaxs=mmax.eq.1
      mmaxp=mmax.le.2
c
c     Doing integrals in subblocks to save memory
c
      maxg=maxxyz/nroots
c
      dtol=1.0d+01**(-itol)
      dtol2=dtol*dtol
      q4=pi252*1.0d0
c
c     ----- pair of k,l primitives -----
c
      ng=0
      do 10 klg = 1, nkl
         db=bcharg( 1,klg)
         bb=bcharg( 2,klg)
         xb=bcharg( 3,klg)
         yb=bcharg( 4,klg)
         zb=bcharg( 5,klg)
         xd=bcharg( 6,klg)
         yd=bcharg( 7,klg)
         zd=bcharg( 8,klg)
         dxkl=bcharg( 9,klg)
         dykl=bcharg(10,klg)
         dzkl=bcharg(11,klg)
         q4db=q4*db
c
c     ----- pair of i,j primitives -----
c
         do 20 ijg = 1, nij
            da=acharg( 1,ijg)
            aa=acharg( 2,ijg)
            xa=acharg( 3,ijg)
            ya=acharg( 4,ijg)
            za=acharg( 5,ijg)
            aandb1=1.0d0/(aa+bb)
            q4dbda=q4db*da
            dum   =q4dbda*q4dbda*aandb1
            if(dum.le.dtol2) go to 20
            q4dbda=q4dbda* sqrt(aandb1)
            rho   =aa*bb*aandb1
            xx=rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
c
            ng=ng+1
            abv(1,ng)=aa
            abv(2,ng)=bb
            abv(3,ng)=rho
            abv(4,ng)=q4dbda
            abv(5,ng)=xx
c
            xc=acharg( 6,ijg)
            yc=acharg( 7,ijg)
            zc=acharg( 8,ijg)
            dxij=acharg( 9,ijg)
            dyij=acharg(10,ijg)
            dzij=acharg(11,ijg)
c
            if(mmaxs) go to 21
            cv( 1,ng)=aa*(xa-xd)
            cv( 2,ng)=bb*(xb-xd)
            cv( 3,ng)=aa*(ya-yd)
            cv( 4,ng)=bb*(yb-yd)
            cv( 5,ng)=aa*(za-zd)
            cv( 6,ng)=bb*(zb-zd)
   21       if(nmaxs) go to 22
            cv( 7,ng)=aa*(xa-xc)
            cv( 8,ng)=bb*(xb-xc)
            cv( 9,ng)=aa*(ya-yc)
            cv(10,ng)=bb*(yb-yc)
            cv(11,ng)=aa*(za-zc)
            cv(12,ng)=bb*(zb-zc)
   22       continue
            cv(13,ng)=dxij
            cv(14,ng)=dyij
            cv(15,ng)=dzij
            cv(16,ng)=dxkl
            cv(17,ng)=dykl
            cv(18,ng)=dzkl
c
            last = (klg.eq.nkl).and.(ijg.eq.nij)
            if (ng.lt.maxg.and.(.not.last)) goto 20 
            numg=ng
            ng=0
            if (numg.eq.0) return 
c
c     ----- compute roots and weights for quadrature -----
c
            call hnd_jkwrys(rwv,abv,numg)
c
c     ----- compute coefficients for recursion formulae -----
c
            call hnd_jkbcdf(b00,b01,b10,c00,d00,f00,dij,dkl,
     1                      abv,cv,rwv,numg)
c
c     ----- compute -x- , -y- , -z- integrals ( 2 centers,2-d ) -----
c
            call hnd_giagnm(gnm,numg*nroots*3,nmax,mmax,
     1                      b00,b01,b10,c00,d00,f00)
c
c     ----- compute -x- , -y- , -z- integrals ( 4 centers,2-d ) -----
c
            call hnd_giaxyz(gijkl,gijkl,gnkl,gnkl,gnkl,gnm,
     1                      numg*nroots*3,nmax,mmax,nimax,njmax,
     2                      nkmax,nlmax,dij,dkl)
c
c     ----- -giao- specific -----
c
            call hnd_xyzgia(gijkl,numg*nroots,nimax,njmax,nkmax,nlmax)
c
c     ----- compute  - (i,j/k,l) -  integrals -----
c
            call hnd_spdgia(numg,nroots,gint,gijkl)
c
   20    continue
   10 continue
c
c     ----- ready for -giao- combinations -----
c
c     call hnd_giamak(gint,ijkl)
 
      return
      end
