c $Id$
C New concept of neglecting integrals (kw 05/26/95)
c==================================================================
      subroutine cshneg(bl,ibl,kbl,nbls,ipres,iapb,icpd, estij,estkl,
     *                   npij,npkl,lcij,lckl,idx1,idx2 )
c---------------------------------------------------------------
c it is called from Calcint2 after pair-precalculations are complited.
c it is used to neglect part of cont. shell quartets at the contracted
c shells level. Ipres(i) is 0 or 1 and shows if i-quartet can be negl.
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
cccc  common /big/ bl(1)
      common /begin/ ijbegin,klbegin
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      dimension estij(npij,lcij),estkl(npkl,lckl)
      dimension ipres(nbls),idx1(nbls),idx2(nbls)
      dimension bl(*)
c---------------------------------------------------------------
      call getmem(npij,ioveij)
      call getmem(npkl,iovekl)
      call getmem(npij,iexpij)
      call getmem(npkl,iexpkl)
c---------------------------------------------------------------
c Now the overlap is : 
c
c  Sab=(a*b)**3/4 * (a+b)**-1 * exp{ -ab/(a+b) * Rab**2 }
c
      call finder(npij,lcij,estij,bl(iapb),ijpar1,ijbegin,
     *            bl(ioveij),bl(iexpij) )
c
c and do the same for pairs cd
c
      call finder(npkl,lckl,estkl,bl(icpd),klpar1,klbegin,
     *            bl(iovekl),bl(iexpkl) )
c---------------------------------------------------------------
c calculate upper bond for the Estimator 
c
      call uprest(nbls,idx1,idx2,ipres, 
     *            bl(ioveij),bl(iovekl),bl(iexpij),bl(iexpkl) )
c     
c on return the Ipres(ijkl) has a value of 0 if it is to be neglect.
c
c---------------------------------------------------------------
      call retmem(4)
c
      end
c==========================================
      subroutine finder(npij,lcij,estij,apb,ijpar1,ijbegin,oveij,expij)
      implicit real*8 (a-h,o-z)
      dimension apb(ijpar1,lcij)
      dimension estij(npij,lcij)
      dimension oveij(npij),expij(npij)
c
      do 10 ijc=1,npij
      ijpar=ijbegin-1+ijc
      coefmax=0.d0
      expomin=1000000.d0
        do 20 ijp=1,lcij
        apb1= apb(ijpar,ijp)
        coef1= estij(ijc,ijp)
        if(coef1.gt.coefmax) coefmax=coef1
        if(apb1.lt.expomin ) expomin=apb1
   20   continue
        oveij(ijc)=coefmax
        expij(ijc)=expomin
   10 continue
      end
c==================================================================
      subroutine uprest(nbls,idx1,idx2,ipres,oveij,ovekl,expij,expkl)
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
      dimension oveij(*),ovekl(*), expij(*),expkl(*)
      dimension ipres(nbls),idx1(nbls),idx2(nbls)
c
      thres=epsr
c
      do 10 ijkl=1,nbls
c it can be already 0 by symmetry 
      if(ipres(ijkl).eq.0) go to 10
      ijc=idx1(ijkl)
      klc=idx2(ijkl)
      overl=oveij(ijc)*ovekl(klc)
      esti1=overl/(expij(ijc)+expkl(klc))
      estim=esti1
      if(estim.le.thres) ipres(ijkl)=0
ctest
c     if(estim.le.thres) then
c        write(6,66) esti1,estim,thres
c  66 format('est1=',e15.5,' est=',e15.5,' thre=',e15.5)
c     endif
ctest
   10 continue
c
      end
c==================================================================
