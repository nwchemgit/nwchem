      subroutine qhop_asite(msa,xs,idon,idong,resd,istatd,iacc,
     + iaccg,resa,istata,xacc,doff,aoff)
c
c $Id: qhop_asite.F,v 1.1 2003-07-17 19:36:51 d3j191 Exp $
c
      implicit none
c
#include "qhop_common.fh"
c
c     declaration of arguments
c

      integer msa,idon,idong,resd,iacc,iaccg,resa,istatd,istata,doff,
     + aoff
      real*8 xacc(5),xs(msa,3)
c
c     declaration of variables
c
      integer typ,i,idr(2),idhv(2),inh,inhv,numh,nhv(2)
      real*8 dist,ang
      real*8 cost,sint,cosp,sinp,th,ph
      real*8 trsl,rds,fct
      real*8 xh2(3),xacc2(2),xdr(3),xr12(3),xar(3),xdrp(2),norm,r1(3),
     + sid1(3),sid2(3),nsid(2),xrh(3),nrh
      real*8 tol,pi
c
      parameter(tol=1E-10)
      parameter(pi=3.14159)      
c
      typ=arat(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2),2)  
      nhv(1)=arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),3)
      nhv(2)=arat(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2),3)
c
c     Calculating indexes and eq. dist and agle
c
      ang=aneq(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2))
      dist=deq(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2))
c
      idhv(1)=arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),1)+idon
      idhv(2)=arat(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2),1)+iacc
c
      idr(1)=idhv(2)+
     + arat(ptarat(ptseq(resa,1))+idong-
     + ptseq(resd,2)+arat(ptarat(ptseq(resd,1))
     + +idong-ptseq(resd,2),1),1)
      idr(2)=idhv(2)+
     + arat(ptarat(ptseq(resa,1))+idong-
     + ptseq(resd,2)+arat(ptarat(ptseq(resd,1))
     + +idong-ptseq(resd,2),1),2)
c
c     H in a "cone" position: typ=1 or OH-
c
      if(typ.eq.1.or.(typ.eq.3.and.istata.eq.3))then
c
      if(typ.eq.3)idr(1)=idhv(2)+1
      xdr(1)=xs(idhv(2),1)-xs(idr(1),1)
      xdr(2)=xs(idhv(2),2)-xs(idr(1),2)
      xdr(3)=xs(idhv(2),3)-xs(idr(1),3)
      call cf_pbc(1,xdr,1,xdr,1,0,1,1)
c     
      norm=sqrt(xdr(1)**2+xdr(2)**2+xdr(3)**2)
      trsl=-dist*cos(ang)+norm
      rds=dist*sin(pi-ang)
c  
      if(xdr(1).lt.tol.and.xdr(2).lt.tol)then
c   
c     no rotation needed
c            
      fct=sqrt(xh2(1)**2+xh2(2)**2)
      if(fct.gt.tol)then
      xacc(1)=xs(iacc,1)*(rds/fct)
      xacc(2)=xs(iacc,2)*(rds/fct)
      else
c
c     H projection on the circle plane coincides with the center of 
c     the circle, all points are equally distant
c     
      xacc(1)=rds
      xacc(2)=0.
      endif
      return
      else
      ph=atan(-(xdr(2)/(xdr(1))))
      endif
c
      sinp=sin(ph)
      cosp=cos(ph)
c
      xdrp(1)=xdr(1)*cos(ph)-xdr(2)*sin(ph)
      if((xdr(3)).lt.tol)then
      th=pi/2.
      else
      th=atan(xdrp(1)/xdr(3))
      endif
      xdrp(2)=sin(th)*xdrp(1)+cos(th)*xdr(3)
      if(xdrp(2).lt.0)th=th+pi
c      
      sint=sin(th)
      cost=cos(th)
c     
c     trasforming the coordinates of the hydrogen in the new system
c     
      xar(1)=xs(iacc,1)-xs(idr(1),1)
      xar(2)=xs(iacc,2)-xs(idr(1),2)
      xar(3)=xs(iacc,3)-xs(idr(1),3)
      call cf_pbc(1,xar,1,xar,1,0,1,1)
c
      xh2(1)=cosp*cost*xar(1)-cost*sinp*xar(2)-sint*xar(3)
      xh2(2)=sinp*xar(1)+cosp*xar(2)
      xh2(3)=sint*cosp*xar(1)-sint*sinp*xar(2)+cost*xar(3)
         
c
c     calculating the new hydrogen coordinates in transformed system
c      
      fct=sqrt(xh2(1)**2+xh2(2)**2)
      if(fct.gt.tol)then
      xacc2(1)=xh2(1)*(rds/fct)
      xacc2(2)=xh2(2)*(rds/fct)
      else
      xacc2(1)=rds
      xacc2(2)=0
      endif
c     
c     Retrasforming the coordinates into the original system 
c     
      xacc(1)=cost*cosp*xacc2(1)+sinp*xacc2(2)+sint*cosp
     +     *trsl+xs(idr(1),1)
      xacc(2)=-cost*sinp*xacc2(1)+cosp*xacc2(2)-sint*sinp
     +     *trsl+xs(idr(1),2)
      xacc(3)=-sint*xacc2(1)+cost*trsl+xs(idr(1),3)
c          
c     H in the "bisetrix" position 
c     
      else if(typ.eq.2.or.(typ.eq.3.and.istata.eq.2))then
c
      if(typ.eq.3)then
c
c     water: O is attacched to the two H
c
      idr(1)=idhv(2)+1
      idr(2)=idhv(2)+2
      endif
c
      xr12(1)=(xs(idr(1),1)+xs(idr(2),1))/2
      xr12(2)=(xs(idr(1),2)+xs(idr(2),2))/2
      xr12(3)=(xs(idr(1),3)+xs(idr(2),3))/2
c
      xrh(1)=xs(idhv(2),1)-xr12(1)
      xrh(2)=xs(idhv(2),2)-xr12(2)
      xrh(3)=xs(idhv(2),3)-xr12(3)
      call cf_pbc(1,xrh,1,xrh,1,0,1,1)
c
      nrh=sqrt(xrh(1)**2+xrh(2)**2+xrh(3)**2)
      xacc(1)=xs(idhv(2),1)+dist*(xrh(1)/nrh)
      xacc(2)=xs(idhv(2),2)+dist*(xrh(2)/nrh)
      xacc(3)=xs(idhv(2),3)+dist*(xrh(3)/nrh)
c
c     H2O: idr(1) and idr(2) are the indexes of the two non dummy H
c
      else       
c 
      endif
c
c     Calculating angle and donor-acceptor dist
c
      r1(1)=xs(idhv(1),1)-xs(idhv(2),1)
      r1(2)=xs(idhv(1),2)-xs(idhv(2),2)
      r1(3)=xs(idhv(1),3)-xs(idhv(2),3)
      call cf_pbc(1,r1,1,r1,1,0,1,1)
      xacc(4)=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
c
      sid1(1)=xs(idhv(1),1)-xs(idon,1)
      sid1(2)=xs(idhv(1),2)-xs(idon,2)
      sid1(3)=xs(idhv(1),3)-xs(idon,3)
      call cf_pbc(1,sid1,1,r1,1,0,1,1)
      sid2(1)=xs(idhv(2),1)-xs(idon,1)
      sid2(2)=xs(idhv(2),2)-xs(idon,2)
      sid2(3)=xs(idhv(2),3)-xs(idon,3)
      call cf_pbc(1,sid2,1,r1,1,0,1,1)
c
      nsid(1)=sqrt(sid1(1)**2+sid1(2)**2+sid1(3)**2)
      nsid(2)=sqrt(sid2(1)**2+sid2(2)**2+sid2(3)**2)
c
      xacc(5)=(sid1(1)*sid2(1)+sid1(2)*sid2(2)+sid1(3)*sid2(3))*
     +         (1/(nsid(1)*nsid(2)))
      xacc(5)=acos(xacc(5))
c
c     Shift for permutation of donor H
c
      numh=arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),4)
      inh=-arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),1)
      inhv=arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),3)
      doff=numh-inh-iand(((istatd-1)/(2**(inhv-1))+1),1)
      if(typ.eq.3)doff=numh-inh-(istatd-1)
c
c     Shift for permutation of acceptor dummy site
c     Currently if nhv>1,only one dummy is allowed
c
      numh=arat(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2),4)
      inh=-arat(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2),1)
      inhv=arat(ptarat(ptseq(resa,1))+iaccg-
     + ptseq(resa,2),3)
      aoff=inh-numh-1+iand(((istata-1)/(2**(inhv-1))+1),1)
      aoff=0
      if(typ.eq.3)aoff=numh-inh-(istata-2)
c      print*,numh,inh,istata,aoff
c
c      print 1021,typ,doff,aoff,xacc,numh,inh,istatd
 1021 format(3i5,5f8.3,3i5)
      return
      end

