      subroutine sif2arr(aodens, d2, mx_d2, valbuf, nipv, labbuf, buf, 
     &   infocd, b2s, sminb, smaxb, shoff, nsp, nshell, dbg)

#include "errquit.fh"
      integer*4 aodens
      double precision d2(mx_d2)
      integer mx_d2
      double precision valbuf(*)
      integer nipv
      integer*4 labbuf(nipv,*)
      double precision buf(*)
      integer*4 infocd(*)
      integer b2s(*)
      integer sminb(*)
      integer smaxb(*)
      integer shoff(nsp,nsp)
      integer nsp
      integer nshell
      integer dbg(4,*)

      integer*4 nbuf, last
      integer i, j, k, l, ij, kl
      integer i0, j0, k0, l0
      integer s, t, u, v, st, uv
      integer s0, t0, u0, v0, st0, uv0
      integer is, js, ks, ls, a, b, c, d, ab, cd, abcd
      integer imin, jmin, kmin, lmin, imax, jmax
      integer slen, tlen, stlen, ulen, vlen
      integer ints, ijkl, ab1, cd1
      logical sort_err
      integer indx
      external indx
      integer idbg
      integer imn, jmn, kmn, lmn

cgk debug
*     write(*,*) 'gk: entered sif2arr'
cgk end

      do ints=1,mx_d2
       d2(ints)=0d0
      enddo
      
  100 continue
 
cgk debug
*     write(*,*)'gk: calling rd_d2bl'
cgk end
      call rd_d2bl(valbuf, labbuf, nbuf, buf, last, aodens, infocd)
cgk debug
*     write(*,*)'gk: back from rd_d2bl'
*     write(*,*)'gk:   nbuf = ', nbuf
cgk end
      
c     loop over integrals in buffer and place in shell block order

      do ints=1,nbuf

c      create mapping from [ij;kl] in AO basis order to [ab;cd] shell
c      blocked order

       sort_err=.false.
       i0=labbuf(1,ints)
       j0=labbuf(2,ints)
       k0=labbuf(3,ints)
       l0=labbuf(4,ints)

*      ij=((i-1)*i)/2+j
*      kl=((k-1)*k)/2+l

       u0=b2s(k0) ! shell for AO k
       v0=b2s(l0) ! shell for AO l
       uv0=(u0*(u0-1))/2+v0
       s0=b2s(i0) ! shell for AO i
       t0=b2s(j0) ! shell for AO j
       st0=(s0*(s0-1))/2+t0

*      write(26,*)'gk:# start index computation'
       if(uv0.gt.st0) then ! switch shell pairs
         i=k0
         j=l0
         k=i0
         l=j0
         s=u0
         t=v0
         u=s0
         v=t0
         st=uv0
         uv=st0
       else
         i=i0
         j=j0
         k=k0
         l=l0
         s=s0
         t=t0
         u=u0
         v=v0
         st=st0
         uv=uv0
       endif

       kmn=sminb(u)
       ks=k-kmn+1 ! min AO index in shell u
       lmn=sminb(v)
       ls=l-lmn+1 ! min AO index in shell v
       ulen=smaxb(u)-kmn+1
       vlen=smaxb(v)-lmn+1

c      construct the (cd) for indexing within the uv shell pair
*      if(u.eq.v) then
*       c=max(ks,ls)
*       d=min(ks,ls)
*       cd=((c-1)*c)/2+d ! linearized index within uv shell pair
*      else
*       c=ks
*       d=js
*       cd=(ls-1)*ulen+ks ! linearized index within uv shell pair
*      endif

       imn=sminb(s)
       is=i-imn+1 ! min index for shell s
       jmn=sminb(t)
       js=j-jmn+1 ! min index of shell t
       slen=smaxb(s)-imn+1
       tlen=smaxb(t)-jmn+1

       abcd=indx(is,js,ks,ls,slen,tlen,ulen,vlen,s,t,u,v,.false.,
     &  imn,jmn,kmn,lmn)

       ijkl=shoff(st,uv)+abcd ! sorted index
       if(dbg(1,ijkl).ne.i) sort_err=.true.
       if(dbg(2,ijkl).ne.j) sort_err=.true.
       if(dbg(3,ijkl).ne.k) sort_err=.true.
       if(dbg(4,ijkl).ne.l) sort_err=.true.

       d2(ijkl)=valbuf(ints)
cgk debug
*       write(26,2321) 's=',s,' t=',t,' u=',u,' v=', v
*       write(26,2321) 'sl=',slen,' tl=',tlen,' ul=',ulen,' vl=', vlen
 2321 format(4(a,i2))
*       write(26,2312) '(st)=', st, ' (uv)=',uv
*       if (st .lt. uv) write(26,*) 'st < uv'
 2312 format(2(a,i2))
*       write(26,2322) '[ab;cd] = [',is,',',js,';',ks,',',ls,']'
*       write(26,2312) 'ab=', ab, ' cd=',cd
 2322   format(4(a,i2),a1)
*       write(26,'(2(a,i6))') 'shoff=',shoff(st,uv), ' abcd=', abcd
*       write(26,'((a,i6))') 'ijkl=',ijkl
*       write(26,2322) '[',i,',',j,';',k,',',l,']'
*       write(26,2324) 'd2(',i,',',j,';',k,',',l,')=',d2(ijkl)
*       write(26,2322) 'debug(',dbg(1,ijkl),',',dbg(2,ijkl),';',
*    &   dbg(3,ijkl),',',dbg(4,ijkl),')'
 2324   format(4(a,i2),a,f20.12)
     
cgk end
       if(sort_err) then
*       write(26,*)'is = ', is
*       write(26,*)'s = ', s
*       write(26,*)'sminb(s) = ', sminb(s)
*       write(26,*)'is = ', is
         do idbg=1,nshell
           write(26,'(a,i2,a,i3)')'sminb(',idbg,') = ', sminb(idbg)
         enddo
        abcd=indx(slen,tlen,ulen,vlen,slen,tlen,ulen,vlen,s,t,u,v,
     &   .true., imn,jmn,kmn,lmn)
        abcd=indx(is,js,ks,ls,slen,tlen,ulen,vlen,s,t,u,v,.true.,
     &   imn,jmn, kmn,lmn)
        call errquit ('sif2arr: sort error', 0, GA_ERR)
       endif

      enddo

      if (last.eq.2) goto 1000
      go to 100     
 1000 continue
cgk debug
*     write(*,*)'gk: leaving sif2arr'
cgk end

      return
      end
