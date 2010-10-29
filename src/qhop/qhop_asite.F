      subroutine qhop_asite(msa,xs,idon,idong,resd,istatd,iacc,
     + iaccg,resa,istata,xacc,doff,aoff,flada)
c
c $Id$
c
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
c
c     declaration of arguments
c
      integer msa,idon,idong,resd,iacc,iaccg,resa,istatd,istata,doff,
     + aoff,flada
      real*8 xacc(11),xs(msa,3)
c
      call qhop_asite2(msa,xs,idon,idong,resd,istatd,iacc,
     + iaccg,resa,istata,xacc,doff,aoff,flada,
     + int_mb(i_ptseq),int_mb(i_ptarat),int_mb(i_iarat),
     + dbl_mb(i_racs),dbl_mb(i_racs+mxar),int_mb(i_ptpar),
     + dbl_mb(i_par),dbl_mb(i_par+7*maxpar),dbl_mb(i_par+14*maxpar),
     + dbl_mb(i_par+17*maxpar),dbl_mb(i_par+23*maxpar),
     + dbl_mb(i_par+28*maxpar))
c
      return
      end
      subroutine qhop_asite2(msa,xs,idon,idong,resd,istatd,iacc,
     + iaccg,resa,istata,xacc,doff,aoff,flada,
     + ptseq,ptarat,arat,deq,aneq,ptpar,tunnel,defe12,zpef,
     + tdsgl,tstval,e12fxy)
c
c $Id$
c
      implicit none
c
#include "qhop_common.fh"
c
c     declaration of arguments
c
      integer msa,idon,idong,resd,iacc,iaccg,resa,istatd,istata,doff,
     + aoff,flada
      real*8 xacc(11),xs(msa,3)
c
      real*8 tunnel(maxpar,7),defe12(maxpar,7),zpef(maxpar,3),
     + tdsgl(maxpar,6),tstval(maxpar,5),e12fxy(maxpar,3)
      real*8 deq(mxar),aneq(mxar)
      integer ptseq(maxseq,2),ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
c
c     declaration of variables
c
      integer typ,i,idr(2),idhv(2),inh,inhv,numh,nhv(2)
      real*8 dist,ang,d2
      real*8 cost,sint,cosp,sinp,th,ph,alpha
      real*8 trsl,rds,fct
      real*8 xh2(3),xacc2(2),xdr(3),xr12(3),xar(3),xdrp(2),norm,r1(3),
     + sid1(3),sid2(3),nsid(2),xrh(3),nrh,r1pbc(3)
      real*8 tol,pi
c
      parameter(tol=1E-15)
      parameter(pi=3.14159265)      
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
c      write (*,98),'angdi ',ang,dist
c 98   format(a6,2f10.6)
c
      ang=1.910610
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
c     H projection on the circle plane coincides with the center of 
c     the circle, all points are equally distant
c     
      xacc(1)=rds+xs(idhv(2),1)
      xacc(2)=xs(idhv(2),2)
      xacc(3)=xs(idhv(2),3)+sign(trsl,xs(idhv(2),3)-xs(idr(1),3))
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
      else if(typ.eq.4)then      
c
c     if acceptor is water then fix exactly the 2 angles
c
      xdr(1)=xs(idhv(2)+1,1)-xs(idhv(2),1)
      xdr(2)=xs(idhv(2)+1,2)-xs(idhv(2),2)
      xdr(3)=xs(idhv(2)+1,3)-xs(idhv(2),3)     
      call cf_pbc(1,xdr,1,xdr,1,0,1,1)
c      write (*,91),'xdr1 ',xdr
      
      if(abs(xdr(1)).gt.tol)then                                  
         if(abs(xdr(2)).gt.tol)then
            th=atan(-xdr(1)/xdr(2))
c            write (*,92),'theta ',th,xdr(1),xdr(2)
 92         format(a6,3f10.6)
         else
            th=-sign(1.d0,xdr(1))*(pi/2)
c            write *,'thetanotol ',xdr(1)-tol,xdr(2)-tol
         endif
      else
         th=0                                              
      endif 
      xdrp(1)=-xdr(1)*sin(th)+xdr(2)*cos(th)                     
                                                              
      if(abs(xdrp(1)).gt.tol)then                               
         if(abs(xdr(3)).gt.tol)then                              
            ph=atan(-xdrp(1)/xdr(3)) 
c            write (*,92),'phi ',ph,xdrp(1),xdr(3)
         else                                                 
            ph=-sign(1.d0,xdrp(1))                      
c            write *,'phinotol ',xdrp(1)-tol,xdr(3)-tol
         endif                                                
      else                                                    
         ph=0                                                
      endif                                                   
      xdrp(2)=-xdrp(1)*sin(ph)+xdr(3)*cos(ph)
      if(xdrp(2).le.0)ph=ph+pi

      xdr(1)=xs(idhv(2)+2,1)-xs(idhv(2),1)
      xdr(2)=xs(idhv(2)+2,2)-xs(idhv(2),2)
      xdr(3)=xs(idhv(2)+2,3)-xs(idhv(2),3)
      call cf_pbc(1,xdr,1,xdr,1,0,1,1)

      call rot(th,xdr(1),xdr(2),xdr(3),3)
      call rot(ph,xdr(1),xdr(2),xdr(3),1)
c      write (*,91),'xdr2 ',xdr
      if(abs(xdr(2)).gt.tol)then
         if(abs(xdr(1)).gt.tol)then
            alpha=atan(xdr(2)/xdr(1))
c            write (*,92),'alpha ',alpha,xdr(2),xdr(1)
         else
            alpha=sign(1.d0,xdr(2))
         endif
      else
         alpha=0
      endif
      xdrp(1)=xdr(1)*cos(alpha)+xdr(2)*sin(alpha)
      if(xdrp(1).le.0)then
         if(alpha.le.0)then
            alpha=alpha+pi
         else
            alpha=alpha-pi
         endif
      endif
c
c     Found the angles, perform rotation
c
      xdr(1)=xs(idhv(1),1)-xs(idhv(2),1)
      xdr(2)=xs(idhv(1),2)-xs(idhv(2),2)
      xdr(3)=xs(idhv(1),3)-xs(idhv(2),3)
      call cf_pbc(1,xdr,1,xdr,1,0,1,1)

      call rot(th,xdr(1),xdr(2),xdr(3),3)
      call rot(ph,xdr(1),xdr(2),xdr(3),1)
      call rot(alpha,xdr(1),xdr(2),xdr(3),3)
      xr12(1)=xdr(1)
      xr12(2)=xdr(2)
      xr12(3)=xdr(3)
      xdr(1)=xs(idhv(2)+2,1)-xs(idhv(2),1)
      xdr(2)=xs(idhv(2)+2,2)-xs(idhv(2),2)
      xdr(3)=xs(idhv(2)+2,3)-xs(idhv(2),3)
      call cf_pbc(1,xdr,1,xdr,1,0,1,1)

      call rot(th,xdr(1),xdr(2),xdr(3),3)
      call rot(ph,xdr(1),xdr(2),xdr(3),1)
      call rot(alpha,xdr(1),xdr(2),xdr(3),3)
c      write (*,91),'xdr3 ',xdr
      d2=dist**2
      xacc(3)=dist*cos(ang)
      xacc(1)=(1/xdr(1))*xacc(3)*(0.1-xdr(3))
      xacc(2)=sqrt(d2-xacc(3)**2-xacc(1)**2)
c      write (*,93),'xacc1 ',xacc(1),xacc(2),xacc(3),d2,cos(ang)
c      write(*,94),'norm ',xacc(1)**2+xacc(2)**2+xacc(3)**2
      d2=abs(xr12(2)-xacc(2))-abs(xr12(2)+xacc(2))
      if(d2.ge.0)xacc(2)=-xacc(2)
      xdr(1)=xacc(1)
      xdr(2)=xacc(2)
      xdr(3)=xacc(3)
      call rot(-alpha,xdr(1),xdr(2),xdr(3),3)
      call rot(-ph,xdr(1),xdr(2),xdr(3),1)
      call rot(-th,xdr(1),xdr(2),xdr(3),3)
c      write (*,94),'norm ',xacc(1)**2+xacc(2)**2+xacc(3)**2
      xacc(1)=xdr(1)+xs(idhv(2),1)
      xacc(2)=xdr(2)+xs(idhv(2),2)
      xacc(3)=xdr(3)+xs(idhv(2),3)
c      write (*,93),'xacc2 ',xacc(1),xacc(2),xacc(3),d2,cos(ang)
 91   format(a5,3f10.6)
 93   format(a6,5f10.6)
 94   format(a5,f10.6)

      else if(typ.eq.5)then
         
         xar(1)=xs(idon,1)-xs(idhv(2),1)
         xar(2)=xs(idon,2)-xs(idhv(2),2)
         xar(3)=xs(idon,3)-xs(idhv(2),3)
         call cf_pbc(1,xar,1,xar,1,0,1,1)
         norm=sqrt(xar(1)**2+xar(2)**2+xar(3)**2)
         xacc(1)=xs(idhv(2),1)+dist*(xar(1)/norm)
         xacc(2)=xs(idhv(2),2)+dist*(xar(2)/norm)
         xacc(3)=xs(idhv(2),3)+dist*(xar(3)/norm)
      endif
c
c     Calculating angle and donor-acceptor dist
c
      r1(1)=xs(idhv(1),1)-xs(idhv(2),1)
      r1(2)=xs(idhv(1),2)-xs(idhv(2),2)
      r1(3)=xs(idhv(1),3)-xs(idhv(2),3)
      r1pbc(1)=r1(1)
      r1pbc(2)=r1(2)
      r1pbc(3)=r1(3)
      call cf_pbc(1,r1,1,r1,1,0,1,1)
      xacc(4)=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
c  
c     calculating midpoint of donor-acceptor
c
      if(abs(r1pbc(1)-r1(1)).lt.0.001) then
      xacc(6)=(xs(idhv(1),1)+xs(idhv(2),1))/2      
      xacc(9)=xacc(6)
      else
      xacc(6)=(xs(idhv(1),1)+(r1pbc(1)-r1(1))+xs(idhv(2),1))/2
      xacc(9)=xacc(6)-(r1pbc(1)-r1(1))
      endif
c
      if(abs(r1pbc(2)-r1(2)).lt.0.001) then
      xacc(7)=(xs(idhv(1),2)+xs(idhv(2),2))/2      
      xacc(10)=xacc(7)
      else
      xacc(7)=(xs(idhv(1),2)+(r1pbc(2)-r1(2))+xs(idhv(2),2))/2
      xacc(10)=xacc(7)-(r1pbc(2)-r1(2))
      endif
c
      if(abs(r1pbc(3)-r1(3)).lt.0.001) then
      xacc(8)=(xs(idhv(1),3)+xs(idhv(2),3))/2      
      xacc(11)=xacc(8)
      else
      xacc(8)=(xs(idhv(1),3)+(r1pbc(3)-r1(3))+xs(idhv(2),3))/2
      xacc(11)=xacc(8)-(r1pbc(3)-r1(3))
      endif
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
c
c     debug: print shift of proton
c
      r1(1)=xs(idon,1)-xacc(1)
      r1(2)=xs(idon,2)-xacc(2)
      r1(3)=xs(idon,3)-xacc(3)
      call cf_pbc(1,r1,1,r1,1,0,1,1)
c      dshop=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
c      print*,"mov ",resd,idong,resa,norm,xacc(4),xacc(5)
c     Shift for permutation of donor H
c
      numh=arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),4)
      inh=-arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),1)
      inhv=arat(ptarat(ptseq(resd,1))+idong-
     + ptseq(resd,2),3)
      if(numh.eq.1)then
      doff=0
      else if(numh.eq.2)then
      doff=numh-inh
c      doff=numh-inh-iand(((istatd-1)/(2**(inhv-1))+1),1)
      else
      doff=numh-inh-(istatd-1)
      endif
c      print*,'don ',resd,inh,doff,istatd
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
c      if(numh>2)aoff=numh-inh-(istata-2)
c     print*,'acc ',resa,inh,doff,istata
c
c      print 1021,typ,doff,aoff,xacc,numh,inh,istatd
c 1021 format(3i5,5f8.3,3i5)
c
c      flag to see if res are same type
c
      flada=0
      if(ptseq(resa,1).eq.ptseq(resd,1))flada=1
      return
      end

      subroutine rot(angle,x,y,z,axis)

      implicit none

      real*8 angle,x,y,z,x1,y1,z1
      integer axis

      x1=x
      y1=y
      z1=z

      if(axis.eq.3)then
         x=x1*cos(angle)+y1*sin(angle)
         y=-x1*sin(angle)+y1*cos(angle)
      else if(axis.eq.1) then
         y=y1*cos(angle)+z1*sin(angle)
         z=-y1*sin(angle)+z1*cos(angle)
      else

      endif
      end


