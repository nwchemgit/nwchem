      subroutine sifr1x(
     & aoints,  info,    itypea,  itypeb,
     & nsym,    nbpsy,   symoff,  buffer,
     & values,  labels,  mapin,   symb,
     & lena,    array,   fcore,   kntin,
     & lasta,   lastb,   last,    nrec,
     & ierr )
c
c  read the designated 1-e integral array and accumulate the matrix
c  elements into array(*).
c
c  this is a no-frills routine to read generic integral arrays.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  itypea, itypeb = generic and specific integral type to be read.
c                   if itypea>=0, then search the entire 1-e integral
c                                 file for candidates.  the file is
c                                 left positioned at the end of the
c                                 1-e integral records.
c                   if itypea<0,  then read the next single array only.
c                                 the file is not repositioned before
c                                 the read.
c                                 the file is left positioned at the end
c                                 of this array.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  symoff(*) = symmetry block offsets used for nontotally symmetric
c              operators.  This includes itypea=1 and itypea=2 arrays.
c              the symmetry blocks within array(*) are offset by
c              symoff( nndxf(isym) + jsym). this allows compact
c              storage of sparse nonsymmetric operator matrices.
c              itypea=0 arrays are referenced in the normal manner,
c              and symoff(*) is not used.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  mapin(*) = input_ao-to-ao mapping vector.
c  lena = length of array(*).
c  array(1:lena) = initial values for the integral array.
c  fcore = initial value for the frozen core contribution.
c
c  output:
c  symb(1:nbft) = symmetry index of each basis function.
c  array(*) = designated 1-e integral array (only if ierr = 0).
c             symmetry block order is determined by the symoff(*) array.
c             diagonal symmetry blocks are stored lower-triangle packed
c             by rows.  off-diagonal symmetry blocks are stored in
c             standard column order.
c  fcore = updated with core contributions from the designated arrays.
c  kntin(*) = number of elements read in each symmetry block.  this
c              array is always referenced as ( nndxf(isym) + jsym ).
c  lasta, lastb = last integral type read and used.
c  last = last parameter from the last record read.
c  ierr = error return code.
c       =  0     for normal return.
c       = -1     eof was found on aoints.
c       = -2     unsatisfied search.
c       = -n*100 if n symmetry blocking errors were detected.
c       >  0     iostat error while reading aoints..
c
c  08-oct-90 (columbus day) 1-e fcore added. sifr1n() interface
c            used.  ierr added. -rls
c  01-jul-90 multiple matrix capability added. -rls
c  15-aug-89 written by ron shepard.
c
      implicit none
c
      integer   nipv,   nmsame,   nomore,   iretbv
      parameter(nipv=2, nmsame=1, nomore=2, iretbv=0)
c
      integer  aoints, itypea, itypeb, nsym,   lena,
     & lasta,  lastb,  last,   nrec,   ierr
      integer  info(*),        nbpsy(nsym),    labels(nipv,*),
     & mapin(*),       symb(*),        symoff(*),      kntin(*)
      real*8   fcore,  buffer(*),      values(*),      array(*)
c
c     # local:
      integer    btypmx
      parameter( btypmx=40 )
      integer btypes(0:btypmx)
      real*8 fcorex(1)
c
      if ( itypea .ge. 0 ) then
         if ( (itypeb .lt. 0) .or. (itypeb .gt. btypmx) ) then
c           # unsupported itypeb value.
            ierr = -3
            return
         else
c           # setup the btypes(*) array for this integral type.
            call izero( (itypeb+1), btypes, 1 )
            btypes(itypeb) = +1
         endif
      endif
c
      fcorex(1) = (0)
c
      call sifr1n(
     & aoints, info,   itypea, itypeb,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  symoff, mapin,
     & lena,   array,  fcorex, symb,
     & kntin,  lasta,  lastb,  last,
     & nrec,   ierr )

ctm
       if (ierr.eq.-4) ierr=0
ctm
c
c     # update the frozen core value.
      fcore = fcore + fcorex(1)
c
      return
      end
