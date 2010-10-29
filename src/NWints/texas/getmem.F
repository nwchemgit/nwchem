c $Id$
c
      subroutine getmem(amount,addr)
c  this subroutine reserves amount space in the common block 
c  bl, checks for max. memory, and returns an integer address
c  bl(addr) to bl(addr+amount-1) is reserved.
c  use: e.g. to reserve 3 matrices,  a(n,n), b(n,3*k),c(n,3*k),
c  and calculate c=a*b. let us assume that only c is needed later,
c  then give c the lower address!!!
c  in the main (calling) program:
c     call getmem(n*3*k,ic)
c     call getmem(n**2,ia)
c     call getmem(n*3*k,ib)
c     call calccmx(bl(ia),bl(ib),bl(ic),n,3*k)
c  ....
c     in the called program calccmx:
c     subroutine calccmx(a,b,c,n,k3)
c     implicit real*8 (a-h,o-z)
c     dimension a(n,n),b(n,k3),c(n,k3)
c     ... give values to a and b
c     call seta(a,n)
c     call setb(b,n,k3)
c     call mtrxmul(a,b,c,n,n,k3)
c     call retmem(2)
c      these last call frees the space taken up by a and b
      integer amount,addr
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c     nreq, maxmem and nadr(1) must be initialized to 0 
c     see retall 
c     lcore is the maximum available memory
      nreq=nreq+1
      addr=nadr(nreq)+1
      ixx=addr+amount-1
      if (ixx.gt.lcore) then
         write(6,*)
     *   ' from getmem: mem. needed=',ixx,' , mem. available=',lcore
         call nerror(1,'getmem',
     1  ' memory overflow : call no., amount requested :',nreq,amount)
         call nerror(1,'getmem',
     1  ' total amount of memory before and after call :',addr,ixx   )
cc>> 1  ' memory overflow request no., amount needed',nreq,ixx)
      end if
      if (ixx.gt.maxmem) maxmem=ixx
      if (nreq.lt.maxall) then
        nadr(nreq+1)=ixx
      else
        call nerror(2,'getmem','too many allocation requests',
     1  nreq,maxall)
      end if
c
c>>d     write(*,*) 'allocation request no=',nreq,'addr=',addr,
c>>d    1 'ceiling=',ixx
      end
c
      subroutine getmem_zero(bl,amount,addr)
      double precision bl(*)
c
c     RJH and HPP ... this routine is IDENTICAL to getmem EXCEPT
c     that it zeroes the data being allocated.  This is a workaround
c     for floating point exceptions that occur on the T3E/Dec ALPHA.
c
c  this subroutine reserves amount space in the common block 
c  bl, checks for max. memory, and returns an integer address
c  bl(addr) to bl(addr+amount-1) is reserved.
c  use: e.g. to reserve 3 matrices,  a(n,n), b(n,3*k),c(n,3*k),
c  and calculate c=a*b. let us assume that only c is needed later,
c  then give c the lower address!!!
c  in the main (calling) program:
c     call getmem(n*3*k,ic)
c     call getmem(n**2,ia)
c     call getmem(n*3*k,ib)
c     call calccmx(bl(ia),bl(ib),bl(ic),n,3*k)
c  ....
c     in the called program calccmx:
c     subroutine calccmx(a,b,c,n,k3)
c     implicit real*8 (a-h,o-z)
c     dimension a(n,n),b(n,k3),c(n,k3)
c     ... give values to a and b
c     call seta(a,n)
c     call setb(b,n,k3)
c     call mtrxmul(a,b,c,n,n,k3)
c     call retmem(2)
c      these last call frees the space taken up by a and b
      integer amount,addr
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c     nreq, maxmem and nadr(1) must be initialized to 0 
c     see retall 
c     lcore is the maximum available memory
      nreq=nreq+1
      addr=nadr(nreq)+1
      ixx=addr+amount-1
      if (ixx.gt.lcore) then
         write(6,*)
     *   ' from getmem: mem. needed=',ixx,' , mem. available=',lcore
         call nerror(1,'getmem',
     1  ' memory overflow : call no., amount requested :',nreq,amount)
         call nerror(1,'getmem',
     1  ' total amount of memory before and after call :',addr,ixx   )
cc>> 1  ' memory overflow request no., amount needed',nreq,ixx)
      end if
      if (ixx.gt.maxmem) maxmem=ixx
      if (nreq.lt.maxall) then
        nadr(nreq+1)=ixx
      else
        call nerror(2,'getmem','too many allocation requests',
     1  nreq,maxall)
      end if
c
#if defined(CRAY)  
      call dcopy(amount, 0.0d0, 0, bl(addr), 1)
#endif
c
c>>d     write(*,*) 'allocation request no=',nreq,'addr=',addr,
c>>d    1 'ceiling=',ixx
      end
c
      subroutine retmem(n)
c     removes the reservation for the n last occupied blocks
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      do 100 i=1,n
        if (nreq.gt.0) nreq=nreq-1
 100  continue
c>>d     write(*,*) 'memory deallocation, last valid request and ceiling=',
c>>d    1  nreq,nadr(nreq+1)
      end
c
      subroutine retall
c     resets the occupied memory to nothing
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      nreq=0
      maxmem=0
      nadr(1)=0
      mark=0
      end
c
      subroutine mmark
c  this subroutine puts down a mark in the memory management if called.
c  its purpose is to facilitate the return of unneded memory
c  retmem serves well if a few block have to be returned but sometimes
c  one does not know the exact number of blocks called. retall is
c  less practical in such  cases because it return everything.
c  mmark is used in conjunction with retmark. If retmark is called,
c  it returns memory allocated after the last call to mmark.
c  mmark and retmark may be used multiply like parentheses, i.e.
c     call mmark
c       call getmem ( ... array1..)
c       call getmem (.. array 2..)
c       call getmem (... array 3)
c       call mmark
c         call getmem(... array 4)
c         call getmem (...array 5)
c       call retmark !returns array locations 4 and 5 to the memory pool
c     call retmark ! returns arrays 1 ,2 and 3 
c
c  The maximum depth of these mmark-retmark parentheses is 100
c   
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      mark=mark+1
      if(mark.gt.100) call nerror(1,'mmark','too many marks, max=100',
     1 mark, nreq)
      marks(mark)=nreq
      end
c
      subroutine retmark
c  see the comments at mmark. Returns memory to the last mark.
c  do not use this in inner loops
c  if the memory has been released by using retmem, this routine
c  will not do anything
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      if(mark.gt.0) then
        newallo=nreq-marks(mark)
        if(newallo.gt.0) then
          call retmem(newallo)
c>>d         write(*,*) 'number of returned blocks=',newallo
        endif
        mark=mark-1
      endif
      end 
c
      subroutine memstat(nreque,nmark,lastadr,memtot)
c   this routine returns the memory status
c   it returns the current memory request number, i.e.
c   the number of allocated blocks, the number of current marks,
c   the last used address, and the total available memory
      parameter (maxall=1000)
      common /mmanag/nreq,maxmem,mark,marks(100),nadr(maxall)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      nreque=nreq
      nmark=mark
      lastadr=maxmem
      memtot=lcore
      end
c
      subroutine getint(amount,addr)
c  this subroutine reserves amount space in the common block 
c  /ctxs_index/, checks for max. memory, and returns an integer address
c  inx(addr) to inx(addr+amount-1) is reserved.
c  use: e.g. to reserve 3 integer matrices,  ii(n,n),jj(n,3*m),kk(n,3*m)
c  and calculate ii=jj*kk. let us assume that only kk is needed later,
c  then give kk the lower address!!!
c  in the main (calling) program:
c     call getint(n*3*m,ic)
c     call getint(n**2,ia)
c     call getint(n*3*m,ib)
c     call calccmx(inx(ia),inx(ib),inx(ic),n,3*m)
c  ....
c     in the called program calccmx:
c     subroutine calccmx(ii,jj,kk,n,m3)
c     implicit real*8 (a-h,o-z)
c     dimension ii(n,n),jj(n,m3),kk(n,m3)
c     ... give values to a and b
c     call seta(ii,n)
c     call setb(jj,n,k3)
c     call mtrxmul(ii,jj,kk,n,n,m3)
c     call retint(2)
c      these last call frees the space taken up by ii and jj
      integer amount,addr
      common /imanag/nreq,maxmem,nadr(400)
      common /ctxs_index/ maxsh,ifp,inx(1)
c     nreq, maxmem and nadr(1) must be initialized to 0 
c     see retinall 
c     maxsh is the maximum available memory
      nreq=nreq+1
      addr=nadr(nreq)+1
      ixx=addr+amount-1
      if (ixx.gt.maxsh) then
         write(*,*) 'integer overflow,needed',ixx,'available',maxsh
         call txs_error
      end if
      if (ixx.gt.maxmem) maxmem=ixx
      nadr(nreq+1)=ixx
c
      end
      subroutine retint(n)
c     removes the reservation for the last occupied block
      common /imanag/nreq,maxmem,nadr(400)
      do 100 i=1,n
        if (nreq.ge.0) nreq=nreq-1
 100  continue
      end
      subroutine retinall
c     resets the occupied memory to nothing
      common /imanag/nreq,maxmem,nadr(400)
      nreq=0
      maxmem=0
      nadr(1)=0
      end
