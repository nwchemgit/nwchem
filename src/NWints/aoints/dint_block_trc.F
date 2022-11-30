      subroutine dint_block_trc(blklen, d2blk, ilo, ihi, jlo, jhi,
     &  klo, khi, llo, lhi, nshell, geri, idatom, nat, frc_2el, fact,
     &  shoff, s, t, u, v, nbft, b2s, dbg)

      integer blklen
      double precision d2blk(*)
      integer  ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer nshell
      double precision geri(llo:lhi, klo:khi, jlo:jhi, ilo:ihi, 3, 4)
      integer idatom(4)
      integer nat
      double precision frc_2el(3,nat)
      double precision fact
      integer shoff
      integer s, t, u, v
      integer nbft
      integer b2s(nbft)
     
      integer slen, tlen, ulen, vlen, is, js, ks, ls
      integer st, uv, a, b, c, d, ab, cd, abcd, stlen, ab1, cd1
      integer i, j, k, l, ijkl
      integer icart, iat
      double precision dE
      integer dbg(4,*)
      integer i0, j0, k0, l0
      logical index_err
      character*2 deriv(3)
      data deriv/'dx','dy','dz'/
      double precision cntr

      st=(s*(s-1))/2+t
      uv=(u*(u-1))/2+v
   
      slen=ihi-ilo+1
      tlen=jhi-jlo+1
      ulen=khi-klo+1
      vlen=lhi-llo+1

*     do iat=1,nat
*      do icart=1,3
*        frc_2el(icart,iat)=0d0
*      enddo
*     enddo

*     write(6,2321) 'slen=',slen,' tlen=',tlen,' ulen=',ulen,' vlen=', 
*    &  vlen
 2321 format(4(a,i2))
*     write(6,2312) '(st)=', st, ' (uv)=',uv
*     write(6,'(a,4(i2,x))')'idatom: ', (idatom(iat), iat=1,4)

      do iat = 1, 4
       if(idatom(iat).ge.1) then
        do icart = 1, 3
*        write(*,'(a,i3)') '!! atom # = ', idatom(iat)
*        write(*,'(a,i3)') '$$ icart = ', icart
         dE = 0d0
         do i=ilo, ihi
          is=i-ilo+1
          do j=jlo, jhi
           js=j-jlo+1
*          s=b2s(a)
*          t=b2s(b)
*          st=nshell*(t-1)-((t-1)*t)/2+s
           if(s.eq.t) then
            a=max(is,js)
            b=min(is,js)
            ab=((a-1)*a)/2+b
            stlen=((slen+1)*slen)/2
           else
            a=is
            b=js
*           ab=(b-1)*slen+a
            ab=(a-1)*tlen+b
            stlen=slen*tlen
           endif
           do k=klo, khi
            ks=k-klo+1
            do l=llo, lhi
             ls=l-llo+1
             if(u.eq.v) then
              c=max(ks,ls)
              d=min(ks,ls)
              cd=((c-1)*c)/2+d 
             else
              c=ks
              d=ls
*             cd=(d-1)*ulen+c
              cd=(c-1)*vlen+d
             endif

             if(st.eq.uv) then
              ab1=max(ab,cd)
              cd1=min(ab,cd)
              abcd=stlen*(cd1-1)-((cd1-1)*cd1)/2+ab1
             else
              abcd=stlen*(cd-1)+ab 
             endif

             ijkl=shoff+abcd
cgk debug
        if(i.lt.j) then
          i0=j
          j0=i
        else
          i0=i
          j0=j
        endif
        if(k.lt.l) then
          k0=l
          l0=k
        else
          k0=k
          l0=l
        endif
        index_err=.false.
        if(dbg(1,ijkl).ne.i0) index_err=.true.
        if(dbg(2,ijkl).ne.j0) index_err=.true.
        if(dbg(3,ijkl).ne.k0) index_err=.true.
        if(dbg(4,ijkl).ne.l0) index_err=.true.
cgk I do not know how to switch bra and kets to make the error trap
cgk absolute
*       if(index_err) then
*         write(6,'(a)')'!!!WARNING Index error'
*         write(6,2322) 'normal ordered [',i0,',',j0,';',k0,',',l0,']'
*         write(6,2322) 'debug index[',dbg(1,ijkl),',',dbg(2,ijkl),';',
*    &     dbg(3,ijkl),',',dbg(4,ijkl),']'
*       endif
 2322   format(4(a,i2),a1)
*       write(6,2322) '[ab;cd] = [',a,',',b,';',c,',',d,']'
*       write(6,2312) 'ab=', ab, ' cd=',cd
 2312   format(2(a,i2))
*       write(6,'(a,i6)') ' abcd=', abcd
*       write(6,'((a,i6))') 'ijkl=',ijkl      
*       write(6,2324) 'd2(',i,',',j,';',k,',',l,')=',d2blk(ijkl)
 2324   format(4(a,i2),a,f20.12)
*       write(6,1200) idatom(iat),deriv(icart),'geri=',
*    &    geri(l,k,j,i,icart,iat)
 1200   format(i2,a2,a,f20.12)
cgk end
*            df=1d0
*            if(k.nq.l) df=df*2d0
*            if(i.nq.j) df=df*2d0
*            if(i.ne.j.and.k.ne.l) df=df*2d0
             cntr=d2blk(ijkl)*geri(l,k,j,i,icart,iat)
             dE = dE + cntr
       if(icart.eq.1.and.abs(cntr).gt.1d-10) then
*       write(27,2324) 'd2(',i,',',j,';',k,',',l,')=',d2blk(ijkl)
*       write(27,1200) idatom(iat),deriv(icart),'geri=',
*    &    geri(l,k,j,i,icart,iat)
*       write(27,*)'cntr = ', cntr
       endif
            enddo
           enddo
          enddo
         enddo
         frc_2el(icart,idatom(iat)) = frc_2el(icart,idatom(iat)) + 
     &    dE*fact
        enddo
       endif
      enddo

*     if (st.lt.uv) then
*       write(*,*)'st < uv'
*     endif
*     write(6,2200) '2 electron energy gradient',' ',' ',
*    $  ((frc_2el(i,j),i=1,3),j=1,nat)
*2200 format(A,A,A/,1000(3(1x,F20.10),/))

cgk debug
*     if(idatom(1).gt.0.or.idatom(2).gt.0.or.idatom(3).gt.0.or.
*    & idatom(4).gt.0) then
*       write(27,2222)'gk: shells: ', s, t, u, v
*2222   format(a,4i2)
*        WRITE(27, 4020) ((frc_2el(i,j),i=1,3),j=1,nat)
*4020  FORMAT (6F12.6)
*     endif
cgk end

      return
      end
