      subroutine int_hf1sp(
     &       xyzi,expi,coefi, i_nprim, i_ngen, Li,
     &       xyzj,expj,coefj, j_nprim, j_ngen, Lj,
     &       xyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &       dryrun,scr,lscr,msg)
c $Id: int_hf1sp.F,v 1.1 1994-06-07 00:29:58 d3e129 Exp $
      implicit none
#include "apispP.fh"
#include "mafdecls.fh"
c::passed
      integer i_nprim, i_ngen, Li
      integer j_nprim, j_ngen, Lj
      integer nat,lscr,lstv
      double precision xyzi(3),expi(i_nprim),coefi(i_nprim,i_ngen)
      double precision xyzj(3),expj(j_nprim),coefj(j_nprim,j_ngen)
      double precision xyz(3,nat),zan(nat),scr(lscr)
      double precision S(lstv),T(lstv),V(lstv)
      logical doS, doT, doV, canAB, dryrun
      character*(*) msg
c::local
      integer offset
c
      if (Li.gt.-1.and.Lj.gt.-1) then
        call hf1(
     &         xyzi,expi,coefi,i_nprim,i_ngen,Li,
     &         xyzj,expj,coefj,j_nprim,j_ngen,Lj,
     &         xyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
        return
c
      elseif (Li.eq.-1.and.Lj.eq.-1) then
c.............................................................. (SP|SP)
c________________________________ (S|S)
        call hf1(
     &         xyzi,expi,coefi,i_nprim,1,0,
     &         xyzj,expj,coefj,j_nprim,1,0,
     &         xyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
c________________________________ (S|P)
        offset = 1
        call hf1(
     &         xyzi,
     &         expi,coefi(1,1),
     &         i_nprim,1,0,
     &         xyzj,
     &         expj,coefj(1,2),
     &         j_nprim,1,1,
     &         xyz,zan,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
c________________________________ (P|S)
        offset = 1 + 3 
        call hf1(
     &         xyzi,
     &         expi,coefi(1,2),
     &         i_nprim,1,1,
     &         xyzj,
     &         expj,coefj(1,1),
     &         j_nprim,1,0,
     &         xyz,zan,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
c________________________________ (P|P)
        offset = 1 + 3 + 3
        call hf1(
     &         xyzi,
     &         expi,coefi(1,2),
     &         i_nprim,1,1,
     &         xyzj,
     &         expj,coefj(1,2),
     &         j_nprim,1,1,
     &         xyz,zan,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
        call spfix1(li,lj,int_mb(kspscat),
     &         2,int_mb(kszspa),2,int_mb(kszspb),
     &         scr,lscr,doS,doT,doV,S,T,V)
        return
      elseif (Li.eq.-1) then
c.............................................................. (SP|X)
c________________________________ (S|X)
        call hf1(
     &         xyzi,expi,coefi,i_nprim,1,0,
     &         xyzj,expj,coefj,j_nprim,j_ngen,Lj,
     &         xyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
c________________________________ (P|X)
        offset = 1*(Lj+1)*(Lj+2)/2
        call hf1(
     &         xyzi,
     &         expi,coefi(1,2),
     &         i_nprim,1,1,
     &         xyzj,
     &         expj,coefj,j_nprim,j_ngen,Lj,
     &         xyz,zan,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
        return
      elseif (Lj.eq.-1) then
c.............................................................. (X|SP)
c________________________________ (X|S)
        call hf1(
     &         xyzi,expi,coefi,i_nprim,i_ngen,Li,
     &         xyzj,expj,coefj,j_nprim,1,0,
     &         xyz,zan,nat,S,T,V,lstv,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
c________________________________ (X|P)
        offset = 1*(Li+1)*(Li+2)/2
        call hf1(
     &         xyzi,
     &         expi,coefi,i_nprim,i_ngen,Li,
     &         xyzj,
     &         expj,coefj(1,2),
     &         j_nprim,1,1,
     &         xyz,zan,nat,
     &         S(offset+1),T(offset+1),V(offset+1),
     &         (lstv-offset)
     &         ,doS,doT,doV,canAB,
     &         dryrun,scr,lscr)
c
        call spfix1(li,lj,int_mb(kspscat),
     &         1,int_mb(kszspa),2,int_mb(kszspb),
     &         scr,lscr,doS,doT,doV,S,T,V)
        return
      endif
      write(6,*)'int_hf1sp called by ',msg
      call errquit('int_hf1sp: not finished',911)
      stop ' okedokee'
c
      end
      subroutine spfix1(li,lj,iscat,
     &       nbi,iszi,nbj,iszj,
     &       scr,lscr,doS,doT,doV,S,T,V)
      implicit none
c
c
c::passed
      integer Li
      integer Lj
      integer lscr
      integer iscat(*)
      integer nbi,nbj,iszi(*),iszj(*)
      double precision scr(lscr)
      logical doS,doT,doV
      double precision S(*), T(*), V(*)
c
      integer sizei,sizej
      integer ii,jj,nint
c
      if (li.lt.0) then
        nbi = abs(li) + 1
        sizei =  0
        do 00100 ii=0,abs(li)
          iszi(ii+1) = (ii+1)*(ii+2)/2
          sizei = sizei + iszi(ii+1)
00100   continue
      else
        nbi = 1
        iszi(1) = (li+1)*(li+2)/2
        sizei = iszi(1)
      endif
c
      if (lj.lt.0) then
        nbj = abs(lj) + 1
        sizej = 0
        do 00101 jj=0,abs(lj)
          iszj(jj+1) = (jj+1)*(jj+2)/2
          sizej = sizej + iszj(jj+1)
00101   continue
      else
        nbj = 1
        iszj(1) = (lj+1)*(lj+2)/2 
        sizej = iszj(1)
      endif
c
      nint = sizei*sizej
      call spscat1(iscat,nbi,iszi,nbj,iszj)
c
      if (doS) then
        call dcopy(nint,S,1,scr,1)
        do 00200 ii = 1,nint
          S(ii) = scr(iscat(ii))
00200   continue
      endif
c
      if (doT) then
        call dcopy(nint,T,1,scr,1)
        do 00300 ii = 1,nint
          T(ii) = scr(iscat(ii))
00300   continue
      endif
c
      if (doV) then
        call dcopy(nint,V,1,scr,1)
        do 00400 ii = 1,nint
          V(ii) = scr(iscat(ii))
00400   continue
      endif
c
      end
      subroutine spscat1(iscat,na,sza,nb,szb)
      implicit none
c
c routine to determine scatter array for sp, spd, ... shells 
c evaluated left/right (a|b) blocked.  i.e., for an (sp|sp) 
c one electron integeral (s|s), (s|p), (p|s) (p|p)     
c
      integer iscat(*)    ! [output] scatter array
      integer na          ! [input] number of blocks in a "L"-block
      integer nb          ! [input] number of blocks in b "L"-block
      integer sza(na)     ! [input] size of each block in a ! (la+1)*(la+2)/2 for 
      integer szb(nb)     ! [input] size of each block in b ! (lb+1)*(lb+2)/2 for 
c
      integer iba,isa,ibb,isb
      integer jba,jsa,jbb,jsb
      integer counti,countj
c
      counti = 0
      do 00100 iba = 1,na
        do 00200 isa = 1,sza(iba)
          do 00300 ibb = 1,nb
            do 00400 isb = 1,szb(ibb)
              counti = counti + 1
              countj = 0
              do 00500 jba = 1,na
                do 00600 jbb = 1,nb
                  do 00700 jsa = 1,sza(jba)
                    do 00800 jsb = 1,szb(jbb)
                      countj = countj + 1
                      if (
     &                    iba.eq.jba.and.
     &                    ibb.eq.jbb.and.
     &                    isa.eq.jsa.and.
     &                    isb.eq.jsb
     &                   ) then
                        iscat(counti) = countj
                        goto 00401
                      endif
00800               continue
00700             continue
00600           continue
00500         continue
00401         continue
00400       continue
00300     continue
00200   continue
00100 continue
c
      end
