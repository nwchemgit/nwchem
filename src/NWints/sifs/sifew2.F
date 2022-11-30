      subroutine sifew2(
     & aoint2,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  ibitv,
     & buffer,  iwait,   nrec,    reqnum,
     & ierr )
c
c  encode and write a 2-e integral record.
c
c  input:
c  aoint2 = output file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 four orbital labels are packed in each labels(*) entry.
c       = 2 two orbital labels are packed in each labels(*) entry.
c       = 4 one orbital label is stored in each labels(*) entry.
c  num = actual number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c  iwait   = asynchronous i/o parameter.
c          = 0  don't wait.  use asynch i/o and return without
c               waiting for i/o completion.  the calling program
c               must call sif2w8() before reusing the buffer or
c               before calling this routine again with the same
c               output buffer.
c          = 1  wait for i/o completion before returning.
c               no sif2w8() call is required in the calling program.
c               buffer(*) can be reused immediately upon return.
c
c  output:
c  num = reset to the number of unwritten elements in values(*).
c        if (last.ne.0) then num is always zero on return. otherwise,
c        the calling program should not assume that this is set to zero.
c        note: this provision is to allow future data-dependent
c              value(*) and labels(*) packing methods. use:
c                   numtot = numtot + num
c                   call (...num...)
c                   numtot = numtot - num
c        in the calling program to compute correctly the total
c        number of output values.  -rls
c  values(1:num) = elements that were not written on this call.
c  labels(1:nipv,1:num) = corresponding unwritten labels.
c  ibitv(1:num) = corresponding unwritten bit-vector elements.
c  buffer(1:l2rec) = packed  buffer.
c  nrec = updated record count.
c         note: the calling program should not assume that this is
c               always incremented by 1.
c  reqnum = i/o request number for the i/o operation associated with
c           this record.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoint2, nipv,   num,    last,   itypea, itypeb,
     & ifmt,   ibvtyp, iwait,  nrec,   reqnum, ierr
      integer  info(*),        labels(*),      ibitv(*)
      real*8   values(*),      buffer(*)
c
      integer  l2rec,  n2max
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l2rec = info(4)
      n2max = info(5)
c
c     # pack the buffer...
c
      call sife2(
     & info,   nipv,   num,   last,
     & itypea, itypeb, ifmt,  ibvtyp,
     & values, labels, ibitv, buffer,
     & ierr )
      if ( ierr .ne. 0 ) return
c
c     # write to the output file.
c
      call sifw2( aoint2, iwait, info, buffer, reqnum, ierr )
      if ( ierr .ne. 0 ) return
c
c     # update nrec, reset num, and move any uncopied values..
c
      nrec = nrec + 1
      num  = 0
c
      return
      end
