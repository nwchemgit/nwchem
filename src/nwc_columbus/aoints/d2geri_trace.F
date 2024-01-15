      subroutine d2geri_trace(mx_d2, d2, lscr, scr, lbuf, buf, nbft,
     &  b2s, nshell, sminb, smaxb, nsp, shoff, nat, frc_2el, ibas,
     &  dbg)
      
      integer mx_d2
      double precision d2(mx_d2)
      integer lscr 
      double precision scr(lscr)
      integer lbuf 
      double precision buf(lbuf)
      integer nbft
      double precision b2s(nbft)
      integer nshell
      integer sminb(nshell)
      integer smaxb(nshell)
      integer nsp
      integer shoff(nsp,nsp)
      integer nat
      double precision frc_2el(3,nat)
      integer ibas
      integer dbg(mx_d2)


      integer ish, jsh, ksh, lsh, ijsh, klsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer slen, tlen, ulen, vlen
      integer i, j, k, l
      integer shf, blklen
      integer idatom(4)
      double precision fact
 
cgk debug
*     write(*,*)'gk: in d2geri_trace'
*     open (unit=27,form='formatted')
cgk end

      do ksh=1,nshell
       klo=sminb(ksh)
       khi=smaxb(ksh)
       ulen=khi-klo+1
       do lsh=1,ksh
*       klsh=nshell*(lsh-1)-((lsh-1)*lsh)/2+ksh
        klsh=(ksh*(ksh-1))/2+lsh
        llo=sminb(lsh)
        lhi=smaxb(lsh)
        vlen=lhi-llo+1
        do ish=ksh,nshell
         ilo=sminb(ish)
         ihi=smaxb(ish)
         slen=ihi-ilo+1
         do jsh=1,ish
*         ijsh=nshell*(jsh-1)-((jsh-1)*jsh)/2+ish
          ijsh=(ish*(ish-1))/2+jsh
          jlo=sminb(jsh)
          jhi=smaxb(jsh)
          tlen=jhi-jlo+1

          if(ijsh.ge.klsh) then
cgk debug
*     write(*,*)'#####################################'
*     write(6,2312) 'ijsh=', ijsh, ' klsh=',klsh
*     write(*,2321) 'd2geri s=',ish,' t=',jsh,' u=',ksh,' v=', lsh
*     write(*,2321) 'ilo=',ilo,' ihi=',ihi,' jlo=',jlo,' jhi=',jhi
*     write(*,2321) 'klo=',klo,' khi=',khi,' llo=',llo,' lhi=',lhi
cgk end
           blklen=slen*tlen*ulen*vlen*4*3
           if(blklen.gt.lbuf) write(*,*)'ERROR in lbuf'
*          call dfill (lbuf, 0d0, buf, 1)
*          call dfill (lscr, 0d0, scr, 1)
           call intd_2e4c(ibas, ish, jsh, ibas, ksh, lsh, lscr, scr,
     &      lbuf, buf, idatom)

*          call print_dint_block(ilo, ihi, jlo, jhi, klo, khi,
*    &      llo, lhi, buf, idatom)

           shf=shoff(ijsh,klsh)
*     write(6,'(a,i5)') 'shoff=',shf
 2312   format(2(a,i2))
           fact=5d-1
*          if(ijsh.lt.klsh) fact=fact*2d0
           if(ish.ne.jsh) fact=fact*2d0
           if(ksh.ne.lsh) fact=fact*2d0
           if(ish.ne.ksh.or.jsh.ne.lsh) fact=fact*2d0

*     write(*,*)'gk: calling dint_block_trc'

           call dint_block_trc(blklen, d2, ilo, ihi, jlo, 
     &      jhi, klo, khi, llo, lhi, nshell, buf, idatom, nat, 
     &      frc_2el, fact, shf, ish, jsh, ksh, lsh, nbft, b2s, dbg)
 2321 format(4(a,i2))

*     write(*,*)'gk: back from dint_block_trc'

          endif
         enddo
        enddo
       enddo
      enddo

*     write(*,*)'gk: leaving d2geri_trace'
*     close(27)

      return
      end
