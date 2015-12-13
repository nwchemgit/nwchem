      subroutine perfexcit(erank,eoper,jr,ir,
     &   ndoexcit,iexfrom,iexto,iexspin,signfact,espin,wtp,wth)
cjp
cjp this soubroutine converts excitation with respect to jr fermi vacuum
cjp to one with respect to ir fermi vacuum and keeps tract of the sign
cjp from possible permutation of 2ndQ operators
cjp output convention uses - signs for spins and operators sorted like a+ b+ j i
cjp no recognition of spin cases is done/neccessary i nthe current implementation
cjp

      implicit none
#include "tce.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "sym.fh"
#include "tce_mrcc.fh"
#include "global.fh"
#include "tce_main.fh"
      integer erank,eoper,jr,ir,i,j,k,is,signfact,erank2,eoper2,espin
      logical nodezero
      integer optyp,ndoexcit,iexfrom,iexto,iexspin,noper,isp,noper2
      integer noper2c,noper2a,eoper2c,eoper2a,k1,ii,k2,l
      integer nca, nab, eoper2rc, eoper2ra,wtp,wth
      dimension eoper(4*maxexcit),optyp(4*maxexcit)
      dimension eoper2(4*maxexcit),espin(4*maxexcit)
      dimension eoper2c(4*maxexcit),eoper2a(4*maxexcit)
      dimension eoper2rc(4*maxexcit), eoper2ra(4*maxexcit)
      dimension iexfrom(8),iexto(8),iexspin(8)
      nodezero = (ga_nodeid().eq.0)

c annihilators
!      dimension eoper(4*maxexcit),optyp(4*maxexcit),eoper2(4*maxexcit)
!      dimension iexfrom(maxexcit),iexto(maxexcit),iexspin(maxexcit)
      erank=0
      signfact=0
      wtp=0
      wth=0
cjp 1. convert input operators from jr-vacuum to absolute orbital numbers
c
       noper=2*ndoexcit
       do i=1,ndoexcit
         optyp(i)=1
         optyp(noper+1-i)=-1
         isp=iexspin(i)
         if(isp.eq.1) then
            eoper(i)=iexto(i)!invpnum(iexto(i),jr,isp)
            eoper(noper+1-i)=iexfrom(i)!invhnum(iexfrom(i),jr,isp)
         else
            eoper(i)=-iexto(i)!invpnum(iexto(i),jr,isp)
            eoper(noper+1-i)=-iexfrom(i)!invhnum(iexfrom(i),jr,isp)
         endif
       enddo
c     if(nodezero)write(100,*)'TESTSIMPLIFY', (eoper(i),i=1,noper),
c    1signfact
c     do i = 1, ndoexcit
c       do j = i+1,ndoexcit
c        if((iexspin(i).eq.iexspin(j)).and.(iexto(i).lt.iexto(j)))then
c        signfact=signfact+1
c        l=iexto(i)
c        iexto(i)=iexto(j)
c        iexto(j)=l
c        endif
c      enddo
c     enddo
c
c      do i = 1, ndoexcit
c       do j = i+1,ndoexcit
c         if((iexspin(i).eq.iexspin(j)).and.
c    1              (iexfrom(i).lt.iexfrom(j)))then
c          signfact=signfact+1
c          l=iexfrom(i)
c          iexfrom(i)=iexfrom(j)
c          iexfrom(j)=l
c        endif
c      enddo
c     enddo  

c      do i = 1, ndoexcit
c        do j = i+1, ndoexcit
c          if((eoper(i)*eoper(j).gt.0).and.
c    1        (abs(eoper(i)).gt.abs(eoper(j)))) then
c             k = eoper(i)
c             eoper(i) = eoper(j)
c             eoper(j) = k
c             signfact=signfact+1
c          endif
c        enddo
c      enddo
c
c      do i = ndoexcit+1, noper
c         do j = ndoexcit+2, noper
c         if ((eoper(i)*eoper(j).gt.0).and.
c    1      (abs(eoper(i)).lt.abs(eoper(j)))) then
c           signfact = signfact+1
c           l = eoper(i)
c           eoper(i) = eoper(j)
c           eoper(j) = l
c         endif
c         enddo
c      enddo
csb
c
cjp 2. combine them with operator string which defines jr reference as excitation from ir vacuum (in absolute orbital numbers already)
       do i=1,2*abs(ihefferank(jr,ir))
       noper=noper+1
       if(i.le.abs(ihefferank(jr,ir))) then
            optyp(noper)=1
       else
            optyp(noper)=-1
       endif
       eoper(noper)=mutualoper(i,jr,ir)
       enddo
cjp 3. simplify the operator string, keep the sign
c note .. optyp=0 used for paired operators
c assumed that no pairing inside the both parts is possible
       do i=1,2*ndoexcit
         do j=2*ndoexcit+1,noper
            if(eoper(i).eq.eoper(j) .and. optyp(i)*optyp(j).lt.0)then
c delete the operators, count how many not deleted operators are in between
              optyp(i)=0
              optyp(j)=0
              do k=i+1,j-1
                 if(optyp(k).ne.0) signfact=signfact+1
              enddo
            endif
         enddo
       enddo
c
       do i = 1, ndoexcit
          if(optyp(i).eq.0)wtp=wtp+1
       enddo
       do i = ndoexcit+1,2*ndoexcit
          if(optyp(i).eq.0) wth=wth+1
       enddo
c     if(nodezero)write(100,*)'TESTSIMPLIFY', (eoper(i),i=1,noper),
c    1wt
cjp 4. squeeze the array, sort it to the a+ b+ j i order, keep the sign
cjp we do NOT care to have same spin for i,a etc!
c k counts the number of annihilators we skip when picking up a creator
       k=0
       noper2=0
       noper2a=0
       noper2c=0
c creators
       do i=1,noper
          if(optyp(i).lt.0) k=k+1
          if(optyp(i).gt.0) then
c             noper2=noper2+1
             noper2c=noper2c+1
c             eoper2(noper2)=eoper(i)
             signfact=signfact+k
             eoper2c(noper2c) = eoper(i)
          endif
       enddo
       do i=1,noper
          if(optyp(i).lt.0) then
c              noper2=noper2+1
              noper2a = noper2a+1
c              eoper2(noper2)=eoper(i)
              eoper2a(noper2a)=eoper(i)
          endif
       enddo
c
c      if(nodezero)write(100,*)"BANIK",(eoper2c(i),i=1,noper2c),
c    1(eoper2a(i),i=1,noper2a)
c	ordering the creation operator sting as AA..BB..
c
      k1 = 0
      nca = 0
      do i = 1, noper2c
         if (eoper2c(i).lt.0) k1=k1+1
         if(eoper2c(i).gt.0) then
         nca = nca+1
         signfact = signfact+k1
         eoper2rc(nca) = eoper2c(i)
         endif
      enddo
      do i = 1, noper2c
        if (eoper2c(i).lt.0) then
        nca = nca+1
        eoper2rc(nca) = eoper2c(i)
        endif
      enddo
c
c 	ordering the annihilation opertator as BB..A..
c
      k1 = 0
      nab = 0
      do i = 1, noper2a
         if(eoper2a(i).gt.0) k1 = k1+1
         if (eoper2a(i).lt.0) then
            nab = nab + 1
            signfact = signfact + k1
            eoper2ra(nab) = eoper2a(i)
         endif
      enddo
      do i = 1, noper2a
         if (eoper2a(i).gt.0) then
           nab = nab + 1 
           eoper2ra(nab) = eoper2a(i)
         endif
      enddo
c
csb  the final orbital indices to build the blocks
c
       do i = 1, noper2c
          noper2 = noper2+1
          eoper2(noper2) = eoper2rc(i)
       enddo
       do i = 1, noper2a
          noper2=noper2+1
          eoper2(noper2) = eoper2ra(i)
       enddo
c        signfact=signfact+signfact1
c
      do i=1,noper2
         if(eoper2(i).gt.0) then
              espin(i)=1
         else
             espin(i) = 2
         endif
              eoper(i)=abs(eoper2(i))
      enddo

csb puting the creation operator in same spin as a<b<.... 
csb
c     do i = 1, noper2c
c        do j = i+1, noper2c
c           if((eoper2rc(i)*eoper2rc(j).gt.0).and.(abs(eoper2rc(i)).gt.
c    1         abs(eoper2rc(j)))) then
c           signfact = signfact+1
c           k = eoper2rc(i)
c           eoper2rc(i) = eoper2rc(j)
c           eoper2rc(j) = k
c           endif
c        enddo
c     enddo
csb
csb puting the annihilation operator in the same order of the creation operators
c
c     do i = 1, noper2a
c        do j = i+1, noper2a
c           if((eoper2ra(i)*eoper2ra(j).gt.0).and.
c    1      (abs(eoper2ra(i)).lt.abs(eoper2ra(j)))) then
c           signfact = signfact+1
c           k = eoper2ra(i)
c           eoper2ra(i) = eoper2ra(j)
c           eoper2ra(j) = k
c           endif
c        enddo
c     enddo

c
c     do i=1,noper2
c        if(eoper2(i).gt.0) then
c            eoper(i)=iphnum(eoper2(i),ir,1)
c        else
c            eoper(i)=-iphnum(-eoper2(i),ir,2)
c        endif
c     enddo
      erank=noper2
      return
      end
      subroutine signf(jr,ir,signfact)
csb
csb this soubroutine gives the sign between two reference function.
csb

      implicit none
#include "tce.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "sym.fh"
#include "tce_mrcc.fh"
#include "global.fh"
#include "tce_main.fh"
      integer jr,ir,i,j,k,signfact,espin
      integer i1,i2
      logical nodezero
      integer iocc0(maxorb,2)
      dimension espin(4*maxexcit)
      nodezero = (ga_nodeid().eq.0)

        do i=1,2
          do j = 1, nmo(i)
            iocc0(j,i)=iocc(j,jr,i)
          enddo
        enddo

       signfact=0

       do i=1,abs(ihefferank(ir,jr))

         j=abs(mutualoper(i,ir,jr))
         if(mutualoper(i,ir,jr).lt.0)then
           espin(j)=2
         else
           espin(j)=1
         endif
         k=abs(mutualoper(2*abs(ihefferank(ir,jr))-i+1,ir,jr))
csb  
csb    Count how many occupied orbitals are there in between
csb
         do i1=min(j,k)+1,max(j,k)-1
           if(iocc0(i1,espin(j)).eq.1)then
              signfact=signfact+1
           endif
         enddo
csb
csb	Interchanges the occupations of holes and particles
csb
           iocc0(j,espin(j)) = 1
           iocc0(k,espin(j)) = 0
       enddo
c
      return
      end

      integer function logtoint(x)
      implicit none
      integer i
      logical x
      if(x.eqv..true.) then
         logtoint=1
      else
         logtoint=0
      endif
      return
      end


