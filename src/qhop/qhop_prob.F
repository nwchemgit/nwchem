      subroutine qhop_prob(gd,rsd,prd,ga,rsa,pra,uda,dang,prob)
c
c $Id: qhop_prob.F,v 1.2 2003-08-07 01:33:55 d3j191 Exp $
c
c     
c     The hopping probability prob is calculated given the difference in 
c     energy between the reactants and the products states e12, and 
c     the distance between donor and acceptor r1.
c     The probability calculation depends on the position in the
c     (e12,r1) plane. 3 different regimes exist: for low e12,r1
c     probs are calculated with an analytical formula obtained
c     fitting several solutions of different Time Dependent Schroedinger 
c     equations (Helms, Lill, J Chem Phys 115, 7985 (2001)).
c     For high values of e12,r1 classical TST with tunnel and zero 
c     point energy corrections is used. In the intermediate regime 
c     probs are calculated via linear interpolation on a logarithmic
c     scale.
c     
      implicit none

c
#include "qhop_common.fh"
c
c     declaration of arguments
c
      real*8 uda(4),dang(2),prob
      integer gd,rsd,prd,ga,rsa,pra
     
c
c     declaration of local variables
c
      real*8 e12,deltae,zpe1
      real*8 tdsglpar1,tdsglpar2
      real*8 kabba,probright,e12right,tstpar1,tstpar2,e12left,
     + preldeltae
      integer i,ipar,index
      integer iseq(2),ihv(2),lhv(2)
c      
c     Searching the correct parameters set ipar
c      
c      print*,'entering prob'
      if(dang(2).le.2.0943)then
      prob=0
      return
      endif

      iseq(1)=ptseq(rsd,1) 
      iseq(2)=ptseq(rsa,1)
      ihv(1)=arat(ptarat(ptseq(rsd,1))+gd-ptseq(rsd,2),3)
      ihv(2)=arat(ptarat(ptseq(rsa,1))+ga-ptseq(rsa,2),3)     
c
      index=lshift(iseq(1),shbit(5))+lshift(ihv(1),shbit(4))+
     + lshift(prd,shbit(3))+lshift(iseq(2),shbit(2))+
     + lshift(ihv(2),shbit(1))+pra
c
      ipar=1
c      print 1026,ipar,index,ptpar(ipar),iseq(1),ihv(1),prd,iseq(2),
c     + ihv(2),pra
c 1026 format('prob ',9i5)
      do while(index.ne.ptpar(ipar))
      ipar=ipar+1 
      if(ipar.gt.maxpar)then
      write(*,1001)
 1001 format('Prob params not found')
      prob=0
      return
      endif
      enddo 
c
c     Calculating e12 
c
      dang(1)=dang(1)*10.
      e12=(uda(3)+uda(4)-uda(1)-uda(2))/4.184
      e12=e12+e12fxy(ipar,1)+e12fxy(ipar,2)*dang(1)+
     + e12fxy(ipar,3)*(dang(1)**2)
c      
c      lhv(1)=ld+arat(ptarat(ptseq(iseq(1),1)+gd-ptseq(iseq(1),2),1)
c      lhv(2)=la+arat(ptarat(ptseq(iseq(2),1)+ga-ptseq(iseq(2),2),1)
c
      tstpar1=tstval(ipar,1)+tstval(ipar,2)*300
      tstpar2=tstval(ipar,3)+tstval(ipar,4)*300
      e12right=tstpar1*(dang(1)-tstpar2)**2+tstval(ipar,5)
      
      if(e12.gt.e12right)then
c
c     TST regime
c        
      preldeltae=defe12(ipar,1)*(dang(1)-defe12(ipar,2))**2
     + +defe12(ipar,3)+defe12(ipar,4)*e12+(defe12(ipar,5)
     + *exp(-defe12(ipar,6)*(dang(1)-2.0))+defe12(ipar,7))*e12*e12
      zpe1=zpef(ipar,1)*exp(-zpef(ipar,2)*preldeltae)+zpef(ipar,3)
      deltae=preldeltae-zpe1
      if(e12.lt.0) then
      kabba=exp(tunnel(ipar,1)+(tunnel(ipar,2)+tunnel(ipar,3)*300
     + +tunnel(ipar,4)*300**2)*preldeltae+(tunnel(ipar,5)
     + +tunnel(ipar,6)*300+tunnel(ipar,7)*300**2)*((preldeltae)**2))
      else
      kabba=exp(tunnel(ipar,1)+(tunnel(ipar,2)+tunnel(ipar,3)*300+
     + tunnel(ipar,4)*300**2)*(preldeltae-e12)+(tunnel(ipar,5)+
     + tunnel(ipar,6)*300+tunnel(ipar,7)*300**2)*(preldeltae-e12)**2)
      endif
      prob=kabba*0.0750*exp(-deltae/(300*0.0019858))
c 
      else
c
      tdsglpar1=tdsgl(ipar,1)*exp(tdsgl(ipar,2)*(dang(1)-2.3))+
     +  tdsgl(ipar,3)
      tdsglpar2=tdsgl(ipar,4)*exp(tdsgl(ipar,5)*(dang(1)-2.3))+
     + tdsgl(ipar,6)
      e12left=(tdsglpar2+1.0986)/tdsglpar1
      if(e12.lt.e12left)then
c
c     low barrier regime
c
      prob=0.5+0.5*tanh(-tdsglpar1*e12+tdsglpar2)
      else
c
c     intermediate regime 
c       
      preldeltae=defe12(ipar,1)*(dang(1)-defe12(ipar,2))**2
     + +defe12(ipar,3)+defe12(ipar,4)*e12right
     + +(defe12(ipar,5)*exp(-defe12(ipar,6)*(dang(1)-2.0)
     + +defe12(ipar,7))*e12right*e12right)
      zpe1=zpef(ipar,1)*exp(-zpef(ipar,2)*preldeltae)+zpef(ipar,3)
      deltae=preldeltae-zpe1
      if(e12right.lt.0) then
      kabba=exp(tunnel(ipar,1)+(tunnel(ipar,2)+tunnel(ipar,3)*300+
     + tunnel(ipar,4)*300**2)*preldeltae+(tunnel(ipar,5)+
     + tunnel(ipar,6)*300+tunnel(ipar,7)*300**2)*preldeltae**2)
      else
      kabba=exp(tunnel(ipar,1)+(tunnel(ipar,2)+tunnel(ipar,3)*300+
     + tunnel(ipar,4)*300**2)*(preldeltae-e12right)+(tunnel(ipar,5)+
     + tunnel(ipar,6)*300+tunnel(ipar,7)*300**2)
     + *(preldeltae-e12right)**2)
      endif
c      
      probright=kabba*0.0750*exp(-deltae/(300*0.0019858))
      prob=0.1*10**(((log10(probright)+1.0)/(e12right-e12left))*(e12
     +        -e12left))
      endif
      endif 
c      
c      print 1101,iseq(1),iseq(2),prd,pra,ptpar(ipar),dang(1),e12,
c     + e12right,e12left,prob
c 1101 format('prob2 ',5i5,5f12.4)
      return
      end

