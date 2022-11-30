      subroutine sifw1x(
     & aoints,  info,    lstflg,  itypea,
     & itypeb,  ifmt,    mapout,  array,
     & nsym,    nbpsy,   symoff,  kntin,
     & buffer,  values,  labels,  fcore,
     & small,   kntout,  numtot,  nrec,
     & ierr )
c
c  write the designated 1-e integral array.
c
c  the output file must be correctly positioned on entry.
c
c  this is a no-frills routine to write generic 1-e integral arrays.
c
c  input:
c  aoints = output file unit number.
c  info(*) = info array for this file.
c  lstflg = flag value to write for the last record.
c  itypea, itypeb = generic and specific integral types.
c  ifmt = output format.
c  mapout(*) = bfn-to-output_bfn mapping vector.
c  array(*) = array to be output.
c             "standard" order of the symmetry blocks and elements
c             within a block is assumed.  see below for details.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  symoff(*) = symmetry block offsets used for nontotally symmetric
c              arrays. the elements are referenced as
c              ( nndxf(isym) + jsym ).  this allows compact
c              storage of sparse nonsymmetric operator matrices.
c              symoff(*) is not referenced for itypea=0 totally
c              symmetric operator arrays.
c  kntin(*) = nominally the number of elements in each symmetry block.
c             this routine only uses the fact that an element is
c             nonzero to determine whether to write out the block.
c             kntin(*) is referenced as (nndxf(isym)+jsym).   when
c             combined with symoff(*), this allows compact storage
c             of sparse nonsymmetric operator matrices in which the
c             blocks are stored arbitrarily in memory.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  fcore = frozen core contribution.
c  small = outout integral cutoff threshold.
c  nrec = initial aoints record count.
c
c  output:
c  array(*) = unmodified, but written to the output file.
c  kntout(1:nndx(nsym+1)) = number of output elements in each symmetry
c                           block.  referenced as (nndxf(isym)+jsym).
c  numtot = actual total number of values written to aoints.
c  nrec = updated record count.
c  ierr = error return code.  0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  15-aug-89  written by ron shepard.
c
      implicit none
c
c     # ibvtyp = 0; bit vectors are assumed to be associated with
c     #             individual integral records, and not with
c     #             the array elements. -rls
c
      integer   nipv,   msame,   ibvtyp
      parameter(nipv=2, msame=0, ibvtyp=0)
c
      integer  aoints, lstflg, itypea, itypeb, ifmt,
     & nsym,   numtot, nrec,   ierr
      integer  info(*),        mapout(*),      nbpsy(nsym),
     & symoff(*),      kntin(*),       labels(nipv,*), kntout(*)
      real*8   fcore,  small
      real*8   array(*),       buffer(*),      values(*)
c
      integer  i,      j,      ij,     ij0,    n1max,
     & num,    i2,     isym,   i1,     ijsym,  last,
     & jsym,   j1,     numtx,  skipd,  j2
      integer  idummy(1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
      nndxf(i) = ( i * (i-1) ) / 2
c
      ierr   = 0
      n1max  = info(3)
      num    = 0
      numtot = 0
      call izero( nndxf(nsym+1), kntout, 1 )
c
      if ( itypea.eq.0 ) then
c
c        # symmetric diagonal-symmetry-blocked input array corresponding
c        # to a totally symmetric operator.
c
         ij0 = 0
         i2  = 0
         do 130 isym = 1, nsym
            i1 = i2 + 1
            i2 = i2 + nbpsy(isym)
c
            ij  = ij0
            ij0 = ij0 + nndxf( nbpsy(isym) + 1 )
c
            ijsym = nndxf(isym+1)
c           # check to make sure this block is to be written.
            if ( kntin(ijsym) .gt. 0 ) then
c
               do 120 i = i1, i2
                  do 110 j = i1, i
                     ij = ij + 1
                     if ( abs(array(ij)).gt.small) then
                        if ( num.eq.n1max ) then
                           numtot = numtot + num
                           last   = msame
                           call sifew1(
     &                      aoints, info,   nipv,   num,
     &                      last,   itypea, itypeb, ifmt,
     &                      ibvtyp, values, labels, fcore,
     &                      idummy, buffer, nrec,   ierr )
                           if ( ierr .ne. 0 ) return
                           numtot = numtot - num
                        endif
                        kntout(ijsym) = kntout(ijsym) + 1
                        num           = num + 1
                        values(num)   = array(ij)
                        labels(1,num) = mapout(i)
                        labels(2,num) = mapout(j)
                     endif
110               continue
120            continue
            endif
130      continue
      elseif ( itypea.eq.1 .or. itypea.eq.2 ) then
c
c        # symmetric and antisymmetric arrays of nontotally symmetric
c        # operators are both stored the same.
c
         if ( itypea .eq. 1 ) then
c           # diagonal elements will be written.
            skipd = 0
         elseif ( itypea .eq. 2 ) then
c           # diagonal elements will be ignored.
            skipd = 1
         endif
c
         i2 = 0
         do 260 isym = 1, nsym
            i1 = i2 + 1
            i2 = i2 + nbpsy(isym)
c
            j2 = 0
            do 230 jsym = 1, (isym-1)
               j1 = j2 + 1
               j2 = j2 + nbpsy(jsym)
c
               ijsym = nndxf(isym) + jsym
               if ( kntin(ijsym) .gt. 0 ) then
c
                  ij = symoff(ijsym)
                  do 220 j = j1, j2
                     do 210 i = i1, i2
                        ij = ij + 1
                        if ( abs(array(ij)) .gt. small ) then
                           if ( num.eq.n1max ) then
                              numtot = numtot + num
                              last   = msame
                              call sifew1(
     &                         aoints, info,   nipv,   num,
     &                         last,   itypea, itypeb, ifmt,
     &                         ibvtyp, values, labels, fcore,
     &                         idummy, buffer, nrec,   ierr )
                              if ( ierr .ne. 0 ) return
                              numtot = numtot - num
                           endif
                           kntout(ijsym) = kntout(ijsym) + 1
                           num           = num + 1
                           values(num)   = array(ij)
                           labels(1,num) = mapout(i)
                           labels(2,num) = mapout(j)
                        endif
210                  continue
220               continue
               endif
230         continue
c
c           # diagonal symmetry block of a nontotally symmetric array.
c
c           # note that diagonal itypea=2 array elements are ignored.
c
            ijsym = nndxf(isym+1)
            if ( kntin(ijsym) .gt. 0 ) then
c
               do 250 i = (i1+skipd), i2
                  ij = symoff(ijsym) + nndxf(i - i1 + 1)
                  do 240 j = i1, (i-skipd)
                     ij = ij + 1
                     if ( abs(array(ij)) .gt. small ) then
                        if ( num.eq.n1max ) then
                           numtot = numtot + num
                           last   = msame
                           call sifew1(
     &                      aoints, info,   nipv,   num,
     &                      last,   itypea, itypeb, ifmt,
     &                      ibvtyp, values, labels, fcore,
     &                      idummy, buffer, nrec,   ierr )
                           if ( ierr .ne. 0 ) return
                           numtot = numtot - num
                        endif
                        kntout(ijsym) = kntout(ijsym) + 1
                        num           = num + 1
                        values(num)   = array(ij)
                        labels(1,num) = mapout(i)
                        labels(2,num) = mapout(j)
                     endif
240               continue
250            continue
            endif
260      continue
      else
         call bummer('sifw1x: unsupported itypea=',itypea,wrnerr)
         ierr = -1
         return
      endif
c
c     # dump the last buffer.
      numtot = numtot + num
      last=lstflg
      call sifew1(
     & aoints, info,   nipv,   num,
     & last,   itypea, itypeb, ifmt,
     & ibvtyp, values, labels, fcore,
     & idummy, buffer, nrec,   ierr )
      if ( ierr .ne. 0 ) return
c
c     check for consistency between kntout(*) and numtot.
      numtx = 0
      do 320 isym = 1, nsym
         do 310 jsym = 1, isym
            numtx = numtx + kntout( nndxf(isym) + jsym )
310      continue
320   continue
      if ( numtx .ne. numtot ) then
         call bummer('sifw1x: (numtx-numtot)=',(numtx-numtot),wrnerr)
         ierr = -2
         return
      endif
c
      return
      end
