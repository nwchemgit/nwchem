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
      integer lbuf, nipbuf, ntape, nnmt, nsym, nbf
      integer itypea,  itypeb
      integer nipv,   msame,  nmsame,     nomore
      parameter(nipv=2, msame=0,nmsame = 1, nomore = 2 )
      integer last, nrec, ierr
      real*8 buf(lbuf),val(nipbuf),fcore
      real*8 d1(nnmt)
      integer ilab(nipv,nipbuf),sym(*),map(*)
      integer nmpsy(*)
      integer infomo(*),kntin(*)
      integer iretbv,    symoff(1),lasta,lastb
      parameter(iretbv = 0)
      real*8    zero
      parameter(zero = 0.0d0)
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
      integer i
      integer numd
c
cgk debug
      !write(6,*)'gk: calling sifr1x from rd1mat'
      !write(6,*)'gk: ntape = ', ntape
      !write(6,*)'gk: itypea = ', itypea
      !write(6,*)'gk: itypeb = ', itypeb
      !write(6,*)'gk: nipbuf = ', nipbuf
      !write(6,*)'gk: lbuf = ', lbuf
      !write(6,*)'gk: nsym = ', nsym
      !write(6,*)'gk: nnmt = ', nnmt
      !write(6,*)'gk: nbf = ', nbf
      !write(6,*)'gk: nmpsy(*) = ', (nmpsy(i),i=1,nsym)
      !write(6,*)'gk: infomo(1) = ', infomo(1)
      !write(6,*)'gk: infomo(2) = ', infomo(2)
      !write(6,*)'gk: infomo(3) = ', infomo(3)
      !write(6,*)'gk: infomo(4) = ', infomo(4)
      !write(6,*)'gk: infomo(5) = ', infomo(5)
      !write(6,*)'gk: symoff(1) = ', symoff(1)
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
      !write(*,*)'gk: calling sifr1x'
      !WRITE(*,*)"LB: itypea=",itypea
      !WRITE(*,*)"LB: itypeb=",itypeb
cgk end
      call sifr1x(
     1 ntape, infomo, itypea, itypeb,
     2 nsym,  nmpsy,  symoff, buf,
     3 val,   ilab,   map,    sym,
     4 nnmt,  d1,     fcore,  kntin,
     5 lasta,lastb,   last,   nrec,
     6 ierr)
cgk debug
      !write(6,*)'gk: back from sifr1x in rd1mat'
      !WRITE(6,*)"LB ierr =",ierr
      !write(6,*)'gk:  nmpsy(1)=', nmpsy(1)
cgk end
      if ( ierr.ne. 0) then
        call bummer('rd1mat: from sifr1x, ierr=',ierr,faterr)
      elseif (fcore .ne. zero) then
          write(6,*)' frozen core orbitals are not yet allowed'
        call bummer('rd1mat: from sifr1x(), fcore=',fcore,wrnerr)
      endif
      return
      end
