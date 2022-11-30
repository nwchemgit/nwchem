      subroutine sifew1(
     & aoints,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  fcore,
     & ibitv,   buffer,  nrec,    ierr )
c
c  encode and write a 1-e integral record.
c
c  input:
c  aoints = output file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 two orbital labels are packed in each labels(*) entry.
c       = 2 one orbital label is stored in each labels(*) entry.
c  num = actual number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n1max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c
c  output:
c  num = reset to the number of unwritten elements in values(*).
c        if (last.ne.0) then num is always zero on return. otherwise,
c        the calling program should not assume that this is set to zero.
c        note: this provision is to allow future data-dependent
c              value(*) and labels(*) packing methods.  use:
c                   numtot = numtot + num
c                   call (...num...)
c                   numtot = numtot - num
c              in the calling program to compute correctly the total
c              number of output values.  -rls
c  values(1:num) = elements that were not written on this call.
c  labels(1:nipv,1:num) = corresponding unwritten labels.
c  ibitv(1:num) = corresponding unwritten bit-vector elements.
c  buffer(1:l1rec) = packed  buffer.
c  nrec = updated record count.
c         note: the calling program should not assume that this is
c               always incremented by 1.
c  ierr = error return code.  0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, nipv,   num,    last,
     & itypea, itypeb, ifmt,   ibvtyp, nrec,   ierr
      integer  info(*),        labels(*),      ibitv(*)
      real*8   values(*),      buffer(*),      fcore
c
      integer  l1rec,  n1max
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
cgk debug
      integer i
cgk end
c
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
cgk debug
c      write(*,*)'gk: **** in sifew1'
c      write(*,*)'gk: aoints = ', aoints
c      write(*,*)'gk: l1rec = ', l1rec
c      write(*,*)'gk: n1max = ', n1max
c      write(*,*)'gk: num = ', num
c      write(*,*)'gk: fcore = ', fcore
c      do i = 1, num
c        write(*,*)'gk: i,values: ', i, values(i)
c      enddo
cgk end
c
c     # pack the buffer...
c
cgk debug
c      write(*,*)'gk:  calling sife1'
c      write(*,*)'LB: num=',num
cgk end
      call sife1(
     & info,   nipv,   num,   last,
     & itypea, itypeb, ifmt,  ibvtyp,
     & values, labels, fcore, ibitv,
     & buffer, ierr )
cgk debug
c      write(*,*)'gk:  back from sife1'
cgk end
      if ( ierr .ne. 0 ) return
c
c     # write to the output file...
c
      call seqwbf( aoints, buffer, l1rec )
cgk debug
c      write(*,*)'gk:  back from seqwbf'
cgk end
c     # seqwbf() does not return ierr.
      ierr = 0
c
c     # update nrec, reset num, and move any unwritten values..
c
      nrec = nrec + 1
      num  = 0
c
      return
      end
