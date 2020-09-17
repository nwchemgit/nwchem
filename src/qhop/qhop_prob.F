      subroutine qhop_prob(gd,rsd,prd,ga,rsa,pra,uda,dang,prob,e120)
c     
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
c
c     declaration of arguments
c
      real*8 uda(4),dang(2),prob,e120
      integer gd,rsd,prd,ga,rsa,pra
c
      call qhop_prob2(gd,rsd,prd,ga,rsa,pra,uda,dang,prob,
     + int_mb(i_ptseq),int_mb(i_ptarat),int_mb(i_iarat),
     + dbl_mb(i_racs),dbl_mb(i_racs+mxar),int_mb(i_ptpar),
     + dbl_mb(i_par),dbl_mb(i_par+7*maxpar),dbl_mb(i_par+14*maxpar),
     + dbl_mb(i_par+17*maxpar),dbl_mb(i_par+23*maxpar),
     + dbl_mb(i_par+28*maxpar),e120)
c
      e120=e120*4.184
c
      return
      end
      subroutine qhop_prob2(gd,rsd,prd,ga,rsa,pra,uda,dang,prob,
     + ptseq,ptarat,arat,deq,aneq,ptpar,tunnel,defe12,zpef,
     + tdsgl,tstval,e12fxy,e120)
c
c $Id$
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
      real*8 uda(4),dang(2),prob,e120
      integer gd,rsd,prd,ga,rsa,pra
c
      real*8 tunnel(maxpar,7),defe12(maxpar,7),zpef(maxpar,3),
     + tdsgl(maxpar,6),tstval(maxpar,5),e12fxy(maxpar,3)
      real*8 deq(mxar),aneq(mxar)
      integer ptseq(maxseq,2),ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
     
c
c     declaration of local variables
c
      real*8 e12,deltae,zpe1
      real*8 tdsglpar1,tdsglpar2
      real*8 kabba,probright,e12right,tstpar1,tstpar2,e12left,
     + preldeltae
      real*8 dist,angle
      integer i,ipar,index
      integer iseq(2),ihv(2),lhv(2)
c
      e12 = 0.0d0
c      
c     Searching the correct parameters set ipar
c      
      dist=dang(1)*10.0d0
      angle=dang(2)
c      print*,'entering prob'
      if(angle.le.2.0943)then
      prob=0
      return
      endif

      iseq(1)=ptseq(rsd,1) 
      iseq(2)=ptseq(rsa,1)
      ihv(1)=arat(ptarat(ptseq(rsd,1))+gd-ptseq(rsd,2),3)
      ihv(2)=arat(ptarat(ptseq(rsa,1))+ga-ptseq(rsa,2),3)     
c
      index=ishft(iseq(1),shbit(5))+ishft(ihv(1),shbit(4))+
     + ishft(prd,shbit(3))+ishft(iseq(2),shbit(2))+
     + ishft(ihv(2),shbit(1))+pra
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
      e120=e12fxy(ipar,1)+e12fxy(ipar,2)*dist+e12fxy(ipar,3)*(dist**2)
      e12=e12+e120+(uda(3)+uda(4)-uda(1)-uda(2))/4.184
c      
c      lhv(1)=ld+arat(ptarat(ptseq(iseq(1),1)+gd-ptseq(iseq(1),2),1)
c      lhv(2)=la+arat(ptarat(ptseq(iseq(2),1)+ga-ptseq(iseq(2),2),1)
c
      tstpar1=tstval(ipar,1)+tstval(ipar,2)*300
      tstpar2=tstval(ipar,3)+tstval(ipar,4)*300
      e12right=tstpar1*(dist-tstpar2)**2+tstval(ipar,5)
      
      if(e12.gt.e12right)then
c
c     TST regime
c        
      preldeltae=defe12(ipar,1)*(dist-defe12(ipar,2))**2
     + +defe12(ipar,3)+defe12(ipar,4)*e12+(defe12(ipar,5)
     + *exp(-defe12(ipar,6)*(dist-2.0))+defe12(ipar,7))*e12*e12
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
      tdsglpar1=tdsgl(ipar,1)*exp(tdsgl(ipar,2)*(dist-2.3))+
     +  tdsgl(ipar,3)
      tdsglpar2=tdsgl(ipar,4)*exp(tdsgl(ipar,5)*(dist-2.3))+
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
      preldeltae=defe12(ipar,1)*(dist-defe12(ipar,2))**2
     + +defe12(ipar,3)+defe12(ipar,4)*e12right
     + +(defe12(ipar,5)*exp(-defe12(ipar,6)*(dist-2.0)
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
c      print 1101,iseq(1),iseq(2),prd,pra,ptpar(ipar),dist,e12,
c     + e12right,e12left,prob
c 1101 format('prob2 ',5i5,5f12.4)
      return
      end

