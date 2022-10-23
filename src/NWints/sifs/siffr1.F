      subroutine siffr1(
     & ninput,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  fcore,
     & ibitv,   ierr )
c
c  read a formatted 1-e integral record.
c  (see also routine siffw1().)
c
c  input:
c  ninput = input file unit number.
c  info(*) = info array for this file.
c
c  output:
c  nipv = number of integers per value returned.
c       = 1 return two orbital labels packed in each labels(*) entry.
c       = 2 return one orbital label in each labels(*) entry.
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format of the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values.
c  labels(1:nipv,1:num) = integral labels
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c  ierr = return code. 0 for normal return.
c
c  24-apr-92 nipv added to the input record. -rls
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  ninput, nipv,   num,   last,
     & itypea, itypeb, ifmt,   ibvtyp, ierr
      integer  info(*),        labels(*),     ibitv(*)
      real*8   values(*),      fcore
c
      integer  i, ioff, lab1
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
c
c     # read in the "dword" information.
c
      read(ninput,*,iostat=ierr)
     & num,    lab1,   ibvtyp, itypea,
     & itypeb, ifmt,   last,   nipv
      if ( ierr .ne. 0 ) then
         call bummer(' siffr1() dword, ierr=', ierr, wrnerr )
         return
      endif
c
      read(ninput,*,iostat=ierr) fcore
      if ( ierr .ne. 0 ) then
         call bummer(' siffr1() fcore, ierr=', ierr, wrnerr )
         return
      endif
c
c     # read the values and labels.
c
      if ( nipv .eq. 1 ) then
         do 10 i = 1, num
            read(ninput,*,iostat=ierr)
     &       values(i), labels(i)
            if ( ierr .ne. 0 ) then
               call bummer(' siffr1() nipv=1, ierr=', ierr, wrnerr )
               return
            endif
10       continue
      elseif ( nipv .eq. 2 ) then
         ioff = 0
         do 20 i = 1, num
            read(ninput,*,iostat=ierr)
     &       values(i), labels(ioff+1), labels(ioff+2)
            if ( ierr .ne. 0 ) then
               call bummer(' siffr1() nipv=2, ierr=', ierr, wrnerr )
               return
            endif
            ioff = ioff + 2
20       continue
      else
         ierr = -10
         return
      endif
c
c     # read the bit vector if present.
c
      if ( ibvtyp .ne. 0 ) then
         read(ninput,*,iostat=ierr) ( ibitv(i), i = 1, num )
         if ( ierr .ne. 0 ) then
            call bummer(' siffr1() ibitv, ierr=', ierr, wrnerr )
            return
         endif
      endif
c
      return
      end
