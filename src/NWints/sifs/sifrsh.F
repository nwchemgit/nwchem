      subroutine sifrsh(
     & aoints,  info,    buffer,  values,
     & labels,  nsym,    nbpsy,   mapin,
     & nnbft,   s1h1,    score,   hcore,
     & symb,    ierr )
c
c  read the overlap and 1-e hamiltonian matrices.
c
c  this is a basic, no-frills, routine to read the 1-e integral arrays
c  necessary for energy calculations.
c
c  on exit, the integral file is positioned after the
c  last 1-e integral record.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  mapin(*) = input_ao-to-ao mapping vector.
c  nnbft = leading dimension of s1h1(*,1:2).
c
c  output:
c  s1h1(*) = the overlap s1(*) matrix is returned in s1h1(*,1) and
c            the total h1(*) matrix is returned in s1h1(*,2).
c            both are symmetry-blocked lower-triangle packed by rows.
c            all 1-e contributions on the integral file are summed
c            into this array.  consequently, the entries on the file
c            must be only the distinct array elements.
c  score = frozen core contribution. tr( s1 * dfc ) = nfrzct
c  hcore = frozen core contribution. tr( h1 * dfc ) = total_core_energy
c  symb(1:nbft) = symmetry index of each basis function
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 if no arrays were found on the integral file.
c       = -n if n symmetry blocking errors were detected.
c       >  0 for iostat error.
c
c  08-oct-90 (columbus day) 1-e fcore change.  sifr1n() interface
c            used. ierr added. -rls
c  04-oct-90 sifskh() call added. -rls
c  26-jul-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, nsym,   nnbft,  ierr
      integer  info(*),        nbpsy(nsym),    labels(2,*),
     & mapin(*),       symb(*)
      real*8   score,  hcore
      real*8   buffer(*),      values(*),      s1h1(nnbft,2)
c
c     # local...
      integer    itypea,   btypmx
      parameter( itypea=0, btypmx=6 )
c
      integer i,      nntot,  isym,   nrec,   last,   lastb,  lasta
      integer symtot(36), idummy(1), btypes(0:btypmx)
      real*8  fcore(2)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
c
c     # set the btypes(*) array:
c     #            0:s1, 1:t1, 2:v1, 3:vec, 4:vfc, 5:vref, 6:generic_h1
      data btypes/ 1,    2,    2,    2,     2,     2,      2          /
c
      nndxf(i) = (i * (i - 1)) / 2
c
c     # nntot = the actual number of elements in the arrays.
      nntot=0
      do 20 isym=1,nsym
         nntot = nntot + nndxf( nbpsy(isym) + 1 )
20    continue
c
      if ( nntot .gt. nnbft ) then
c        # inconsistent nnbft value.
         call bummer('sifrsh: (nntot-nnbft)=',(nntot-nnbft),wrnerr)
         ierr = -2
         return
      endif
c
c     # initialize the output arrays...
c
      call wzero( nntot, s1h1(1,1), 1 )
      call wzero( nntot, s1h1(1,2), 1 )
      fcore(1) = (0)
      fcore(2) = (0)
c
      call sifr1n(
     & aoints, info,   itypea, btypmx,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  idummy, mapin,
     & nnbft,  s1h1,   fcore,  symb,
     & symtot, lasta,  lastb,  last,
     & nrec,   ierr )

c      call plblks('s1',s1h1(1,1),nsym, nbpsy, 'MO',3,6)
c      call plblks('h1',s1h1(1,2),nsym, nbpsy, 'MO',3,6)

ctm
c      sifr1n (recently) returns ierr=-4 if some of the
c      requested integral types are not retrievable from
c      the aoints file; however this situation occurs
c      frequently for effective core potentials etc.
       if (ierr.eq.-4) ierr=0
ctm

c
c     # save the appropriate core values.
      score = fcore(1)
      hcore = fcore(2)
c
      return
      end
