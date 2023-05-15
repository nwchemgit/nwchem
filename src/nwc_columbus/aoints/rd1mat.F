      subroutine rd1mat(ntape, lbuf, nipbuf, infomo, buf, ilab, val,
     & d1, nnmt, sym, map, nsym, nmpsy, kntin, itypea, itypeb, nbf)
c
c  read a 1-e density matrix
c
c  input:
c  ntape = input unit number.
c  lbuf = buffer length.
c  nipbuf = maximum number of integrals in a buffer.
c  infomo = info array
c  buf(*) = buffer for integrals and packed labels.
c  ilab(*)= unpacked labels.
c  val(*) = unpacked values.
c  sym(*) = symmetry of each orbital
c  nnmt = total size of symmetric, symmetry-packed, lower-triangular-
c         packed matrix.
c  map(*) = orbital index mapping vector.
c
c  output:
c  d1(*) = 1-e density matrix
c
      implicit none
c
      integer*4 lbuf, nipbuf, ntape, nnmt, nsym, nbf
      integer*4 itypea,  itypeb
      integer*4 nipv,   msame,  nmsame,     nomore
      parameter(nipv=2, msame=0,nmsame = 1, nomore = 2 )
      integer*4 last, nrec, ierr
      real*8 buf(lbuf),val(nipbuf),fcore
      double precision d1(nnmt)
      integer*4 ilab(nipv,nipbuf),sym(*),map(*)
      integer*4 nmpsy(*)
      integer*4 infomo(*),kntin(*)
      integer*4 iretbv,    symoff(1),lasta,lastb
      parameter(iretbv = 0)
      real*8    zero
      parameter(zero = 0.0d0)
      integer*4   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
      integer*4 i
      integer numd
c
cgk debug
*     write(*,*)'gk: calling sifr1x from rd1mat'
*     write(*,*)'gk: ntape = ', ntape
*     write(*,*)'gk: itypea = ', itypea
*     write(*,*)'gk: itypeb = ', itypeb
*     write(*,*)'gk: nipbuf = ', nipbuf
*     write(*,*)'gk: lbuf = ', lbuf
*     write(*,*)'gk: nsym = ', nsym
*     write(*,*)'gk: nnmt = ', nnmt
*     write(*,*)'gk: nbf = ', nbf
*     write(*,*)'gk: nmpsy(*) = ', (nmpsy(i),i=1,nsym)
*     write(*,*)'gk: infomo(1) = ', infomo(1)
*     write(*,*)'gk: infomo(2) = ', infomo(2)
*     write(*,*)'gk: infomo(3) = ', infomo(3)
*     write(*,*)'gk: infomo(4) = ', infomo(4)
*     write(*,*)'gk: infomo(5) = ', infomo(5)
*     write(*,*)'gk: symoff(1) = ', symoff(1)
cgk end
      fcore= zero
      symoff(1)=0
c
cgk debug
*     write(*,*)'gk:calling dfill'
*     write(*,*)'  nnmt = ', nnmt
*     call flush(6)
cgk end
      numd=nnmt
      call dfill(numd,0d0,d1,1)
cgk debug
*     write(*,*)'gk:back from dfill'
*     call flush(6)
cgk end
      do i=1,nbf
        map(i)=i
      enddo
cgk debug
*      write(*,*)'gk: calling sifr1x'
*      WRITE(*,*)"LB: itypeb=",itypeb
cgk end
      call sifr1x(
     1 ntape, infomo, itypea, itypeb,
     2 nsym,  nmpsy,  symoff, buf,
     3 val,   ilab,   map,    sym,
     4 nnmt,  d1,     fcore,  kntin,
     5 lasta,lastb,   last,   nrec,
     6 ierr)
cgk debug
*     write(*,*)'gk: back from sifr1x in rd1mat'
*     write(*,*)'gk:  nmpsy(1)=', nmpsy(1)
cgk end
      if ( ierr.ne. 0) then
        call bummer('rd1mat: from sifr1x, ierr=',ierr,faterr)
      elseif (fcore .ne. zero) then
          write(6,*)' frozen core orbitals are not yet allowed'
        call bummer('rd1mat: from sifr1x(), fcore=',fcore,wrnerr)
      endif
      return
      end
