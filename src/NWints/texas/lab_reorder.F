      subroutine reordertxs(ncs,inx,iny,ncshell,ncfunct,datnuc,datbas)
c $Id$
      implicit none
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      integer ncs
      integer inx(12,*),iny(12,*)
      integer ncshell(ncs), ncfunct(*)
      double precision datnuc(5,*), datbas(13,*)
c
c     Local
c
      integer icf,ics,i,k,ii,iexch,itemp
      double precision nzi0,nzi1
      double precision exp0,exp1
      integer ics0,ics1
      integer ina0,ina1
      integer icb0,icb1
      integer ict0,ict1
      integer iclen0,iclen1
c------------------------------------
c iny - keeps the original basis set 
c------------------------------------
      do 100 ics=1,ncs
      do 200 i=1,12
       iny(i,ics)=inx(i,ics)
  200 continue
  100 continue
c------------------------------------
c    reorder the  basis set from the lowest to the highest ang.
c     momentum. within one ang. momentum, reorder them so that
c     the contraction depth decreases. Within one ang. mom. and
c     contr. depth, reorder them so that the basis functions
c     for the same atom follow each other.
c     reorder first the atoms so that the same symmetry type
c     for the atoms follow each other
c      maybe this latter is not needed 
c find the lowest ang. momentum, and within that the highest contr.
c length and within that the lowest atom number in the new
c atomic ordering       
      if (util_print("basis_order",print_debug).and.
     +    ga_nodeid().eq.0) then
        write(luout,*) 'nwchem basis set'
        do ics=1,ncs
          write(luout,666) ics,(inx(ii,ics),ii=1,5),
     +                         (inx(9+ii,ics),ii=1,3)
        enddo
      endif
c------------------------------------
 400  continue
      iexch=0
      do 500 i=1,ncs-1
        ics0=i
        ics1=i+1
c
c       square of type (1=s,2=p,3=sp,4=d,5=f,...) and add number of functions
c       in shell to distinguish between 5d and 6d shells, etc.
        ict0=inx(12,ics0)**2+inx(3,ics0)
        ict1=inx(12,ics1)**2+inx(3,ics1)
c
        ina0=inx(2,ics0)
        ina1=inx(2,ics1)
charge:
        if(ina0.eq.0) then
          nzi0=0
        else
          nzi0=datnuc(1,ina0)
        endif
        if(ina1.eq.0) then
          nzi1=0
        else
          nzi1=datnuc(1,ina1)
        endif
c
c contr. length :
        iclen0=inx(5,ics0)-inx(1,ics0)
        iclen1=inx(5,ics1)-inx(1,ics1)
c 
c first primitive exponent:
        icb0=inx(1,ics0)+1
        icb1=inx(1,ics1)+1
        exp0=datbas(1,icb0)
        exp1=datbas(1,icb1)
c
        if (ict1.lt.ict0) then
          iexch=1
          do 410 k=1,12
            itemp=inx(k,ics0)
            inx(k,ics0)=inx(k,ics1)
            inx(k,ics1)=itemp
 410      continue
        else if (ict1.eq.ict0.and.iclen0.lt.iclen1) then
          iexch=1
          do 420 k=1,12
            itemp=inx(k,ics0)
            inx(k,ics0)=inx(k,ics1)
            inx(k,ics1)=itemp
 420      continue
        else if(ict1.eq.ict0.and.iclen0.eq.iclen1
     *                      .and.nzi0.gt.nzi1) then
c
ctxs *                            .and.nzi1.gt.nzi0) then
ctxs :    from heavier to lighter .......................
cpnl *                            .and.nzi0.gt.nzi1) then
cpnl :    from lighter to heavier .......................
          iexch=1
          do 430 k=1,12
            itemp=inx(k,ics0)
            inx(k,ics0)=inx(k,ics1)
            inx(k,ics1)=itemp
 430      continue
        else if(ict1.eq.ict0.and.iclen0.eq.iclen1
     *                      .and.nzi0.eq.nzi1
     *                      .and.exp0.gt.exp1) then
ccc  *                      .and.exp0.lt.exp1) then
c.......  from BIG exp to SMALL exp
c
ccc  *                      .and.exp0.gt.exp1) then
c.......  from SMALL exp to BIG exp
c06 NEW
          iexch=1
          do 440 k=1,12
            itemp=inx(k,ics0)
            inx(k,ics0)=inx(k,ics1)
            inx(k,ics1)=itemp
 440      continue
        end if
cccccc  endif
 500  continue
c     now we have to re-generate the beginning-ending contraction
c     arrays (inx(11,i) and inx(10,i))
      icf=0
      do 700 i=1,ncs
        inx(11,i)=icf
c       end of contr= beginning+(shell-size)*(1+number of gen. contr.)
        inx(10,i)=icf+inx(3,i)*(1+inx(4,i))
        icf=inx(10,i)
 700  continue
      if (iexch.eq.1) go to 400
c------------------------------------
      if (util_print("basis_order",print_debug).and.
     +    ga_nodeid().eq.0) then
        write(luout,*) 'texas basis set'
        do ics=1,ncs
         write(luout,666) ics,(inx(ii,ics),ii=1,5),
     +                        (inx(9+ii,ics),ii=1,3)
 666     format(' ics=',i3,' ib=',i3,' nat=',i3,' shz=',i3,' gc=',i3,
     1          ' ie=',i3,' lcf=',i3,' fcf=',i3, ' type=',i2)
        enddo
      endif
c------------------------------------
c     write(6,*)'----------------------------------------'
c     write(6,*)'       reordered basis set'
c     write(6,*)'----------------------------------------'
c     write(6,*)' shell,center,itype,icont,igenc,ichar, exponent'
c     write(6,*)'----------------------------------------'
c     do ics=1,ncs
c       itype=inx(12,ics)                 ! cond 1
c       icont=inx(5,ics)-inx(1,ics)       ! cond 2
c       igenc=inx(4,ics)                  ! cond 3
c         iatom=inx(2,ics)
c       ichar=datnuc(1,iatom)             ! cond 4
c         icb=inx(1,ics)+1
c       expon=datbas(1,icb)               ! cond 5
c       write(6,777) ics,iatom,itype,icont,igenc,ichar,expon
c     enddo
c 777 format(6(i5,1x),f15.8)
c------------------------------------
c  set up a basis set relation between 
c  original (PNL) and re-ordered (TXS)
c  like this 
c            ncshell(ics_old)----->ics_new
c                      pnl          texas
c and 
c            ncfunct(icf_new)----->icf_old
c                      texas         pnl
c
      call new_old(ncs,inx,iny,ncshell,ncfunct) 
c------------------------------------
      if (util_print("basis_order",print_debug).and.
     +    ga_nodeid().eq.0) then
        write(luout,*) 'translation shell-nwchem -> shell-texas'
        do ics = 1, ncs
          write(luout,*)ics,ncshell(ics)
        enddo
      endif
      end
      subroutine new_old(ncs,inx,iny,ncshell,ncfunct) 
      implicit none
      integer ncs
      integer inx(12,*),iny(12,*)
      integer ncshell(*),ncfunct(*)
c
      integer iiii
      integer icsx,icsy
      integer icfx,icfy
      integer itypx,itypy
      integer iatox,iatoy
      integer icobx,icoby
      integer icoex,icoey
      integer igenx,igeny
      integer ilenx,ileny
c
      do icsx=1,ncs
        itypx=inx(12,icsx)
        iatox=inx(2,icsx)
        icobx=inx(1,icsx)
        icoex=inx(5,icsx)
        igenx=inx(4,icsx)
        ilenx=inx(10,icsx)-inx(11,icsx)
c    
        do icsy=1,ncs
          itypy=iny(12,icsy)
          iatoy=iny(2,icsy)
          icoby=iny(1,icsy)
          icoey=iny(5,icsy)
          igeny=iny(4,icsy)
          ileny=iny(10,icsy)-iny(11,icsy)
c
          if(itypy.eq.itypx) then
            if(icoby.eq.icobx) then
              if(icoey.eq.icoex) then
                if(igeny.eq.igenx) then
                  if(iatoy.eq.iatox) then
                    if(ileny.eq.ilenx) then
c
                      ncshell(icsy)=icsx 
c
                      do iiii = 1, ilenx
                        icfx=inx(11,icsx)+iiii
                        icfy=iny(11,icsy)+iiii
                        ncfunct(icfx)=icfy
                      enddo ! iiii
                      go to 30
c
                    endif
                  endif
                endif
              endif
            endif
          endif
        enddo ! icsy
   30   continue
      enddo ! icsx
c-----------------------------------------------------------
c contracted function mapping from TXS to PNL
c
c      write(8,*)' ncf_pnl=',iny(10,ncs),' ncf_txs=',inx(10,ncs)
c-----------------------------------------------------------
ctest
c     write(8,*) 'from re-order : shell_pnl -----> shell_txs '
c     do 100 ii=1,ncs
c     write(8,88) ii, ncshell(ii)
c 100 continue
c  88 format(20x,i5,8x,i5)
ctest
c     write(8,*) 'from re-order : funct_txs -----> funct_pnl '
c     do 110 icf=1,inx(10,ncs)
c     write(8,*)'       icf_txs=',icf,'  icf_pnl=',ncfunct(icf)
c 110 continue
c-----------------------------------------------------------
      end
