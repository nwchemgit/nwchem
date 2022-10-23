      subroutine sifcfg(
     & itype,   lrecal,  nbft,    ibvtyp,
     & ifmt,    lrec,    nmax,    ierr )
c
c  return a set of consistent configuration parameters for a standard
c  integral file structure.
c
c  input:
c  itype = integral type.
c          1 for 1-e integrals,
c          2 for 2-e integrals.
c  lrecal = maximum buffer length to be allocated.
c         =-1 is a special case for default output values.
c  nbft = total number of basis functions.
c  ibvtyp = 0 if no bit vector is to be stored.
c         .ne.0 if a bit vector is to be stored in the record buffer.
c
c  notes: (1) ifmt, lrec, and nmax may eventually have meaning on input
c         in future versions of this routine.  for now the input values
c         are ignored.
c         (2) for extensibility, the input variables may be passed
c         into this routine as array elements, instead of scalars, in
c         future versions of this routine.
c
c  output:
c  ifmt = ifmt parameter.
c  lrec = actual buffer length to be written.
c  nmax = maximum number of values in each record.
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 for itype error.
c       = -2 for nbft error.
c       = -3 for lrec error.
c
c  08-oct-90 (columbus day) 1-e fcore change. chunk added. -rls
c  05-jul-89 ibvtyp added. -rls
c  01-aug-89  written by ron shepard.
c
      implicit none
c
      integer itype, lrecal, nbft, ibvtyp, ifmt, lrec, nmax, ierr
c
      integer nbig, chunk, nword, ibvfac, space
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # nbftmx = largest nbft consistent with label packing methods.
c     # lrecmn = minimum allowed record length.
c     # lrecmx = maximum allowed record length.  this should be
c     #          consistent with the record "dword" encoding.
c     # lrcinc = record length block size increments.
c     # ndeflt = default number of block increments. such that
c     #          the default record length = (lrecmn + ndeflt * lrcinc)
c
      integer    nbftmx
      parameter( nbftmx=2**16-1 )
      integer    lrecmn,     lrecmx,         lrcinc,      ndeflt
      parameter( lrecmn=2048, lrecmx=2**16,   lrcinc=2048,  ndeflt=1 )
cgk debug
      write(*,*)'gk: in sifcfg'
cgk end
c
c     # check the input parameters...
c
      ierr = 0
c
      if ( itype .ne. 1    .and.    itype .ne. 2 ) then
         ierr = -1
         return
      elseif ( nbft .le. 0    .or.    nbft .gt. nbftmx ) then
         ierr = -2
         return
      endif
c
c     # adjust lrec for efficient I/O.
c     # machine-dependent quirks, magic record lengths, etc. should
c     # be localized here and in the above parameter definitions..
c
      if ( lrecal .eq. -1 ) then
c        # return the default values for everything.
         lrec = lrecmn + (ndeflt * lrcinc)
      else
c        # want:     lrec = lrecmn + (n * lrcinc) <= lrecal
c        #           for the largest possible n.
         lrec = min( lrecmx, lrecal )
         lrec = lrecmn + ( ((lrec - lrecmn) / lrcinc ) * lrcinc )
      endif
      lrec = min( lrecmx, lrec )
      if ( lrec .lt. lrecmn ) then
         ierr = -3
         return
      endif
c
c     # for this record length, compute the number of values.
c
c     # 1-e records: dword // values // fcore // packed_labels // ibitv
c     # 2-e records: dword // values //          packed_labels // ibitv
c
      if ( nbft .lt. 2**8 ) then
c        # use standard 8-bit orbital label packing.
         ifmt=0
      else
c        # use standard 16-bit orbital label packing.
         ifmt=1
      endif
c
c     # compute a reasonable upper bound to nmax.
c     # chunk is used to constrain nmax so that ulab8() or ulab16() do
c     # not overwrite the labels(*) array.
c     #
c     # another choice would be to require the programmer to account
c     # for these overruns when allocating labels(*).  this would
c     # result in more efficient buffer(*) use, but the chunk choice
c     # is simpler.
c
c     # note that the programmer must explicitly account for ibitv(*)
c     # overruns; constraining chunk to 64 is too wasteful. -rls
c
      if ( itype .eq. 1 ) then
         if ( ifmt .eq. 0 ) then
c           # n + n/4 <= lrec
            nbig  = ( 4 * lrec ) / 5
            chunk = 4
         elseif ( ifmt .eq. 1 ) then
c           # n + n/2 <= lrec
            nbig  = ( 2 * lrec ) / 3
            chunk = 2
         endif
      elseif ( itype .eq. 2 ) then
         if ( ifmt .eq. 0 ) then
c           # n + n/2 <= lrec
            nbig  = ( 2 * lrec ) / 3
            chunk = 2
         elseif ( ifmt .eq. 1 ) then
c           # n + n <= lrec
            nbig  = lrec / 2
            chunk = 1
         endif
      endif
c
c     # round up to a higher multiple of chunk.
      nbig = ( ( nbig + 2*chunk ) / chunk ) * chunk
c
c     # account for dword and fcore.
      if ( itype .eq. 1 ) then
         nword = 2
      else
         nword = 1
      endif
c
c     # packed_bit_vector_space = ibvfac * ((n+63)/64)
      if ( ibvtyp .ne. 0 ) then
         ibvfac = 1
      else
         ibvfac = 0
      endif
c
c     # loop backwards from the upper bound until a valid nmax is found.
      do 10 nmax = nbig, 0, -chunk
         space = nword + nmax + ( ( nmax + chunk - 1) / chunk ) +
     &    ibvfac * ( ( nmax + 63 ) / 64 )
         if ( space .le. lrec ) goto 11
10    continue
11    continue
      if ( nmax .le. 0 ) then
         ierr = -4
         return
      endif
c
c     # ifmt, lrec, and nmax are all ok.
c
      return
      end
