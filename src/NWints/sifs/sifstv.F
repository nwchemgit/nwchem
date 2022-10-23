      subroutine sifstv(
     & aoints,  info,    buffer,  values,
     & labels,  nsym,    nbpsy,   mapin,
     & nnbft,   stv,     score,   tcore,
     & vcore,   symb,    ierr )
c
c  read the overlap, kinetic energy, and potential energy arrays.
c
c  this is a basic, no-frills, routine to read the 1-e integral arrays
c  necessary for energy and virial ratio calculations.
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
c  nnbft = leading dimension of stv(*,1:3).
c
c  output:
c  stv(*) = the overlap s1(*) matrix is returned in stv(*,1), the total
c           kinetic matrix T1_total(*) is returned in stv(*,2), and the
c           total 1-e potential matrix V1_total(*) is returned in
c           stv(*,3).
c           all arrays are returned symmetry-blocked lower-triangle
c           packed by rows.
c           all 1-e contributions on the integral file are summed
c           into this array.  consequently, the entries on the file
c           must be only the distinct array elements.
c  score = frozen core contribution. tr( s1 * dfc ) = nfrzct
c  tcore = frozen core contribution. tr( t1 * dfc )
c  vcore = frozen core contribution. tr( v1 * dfc )
c  symb(1:nbft) = symmetry index of each basis function
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 if no arrays were found on the integral file.
c       = -n if n symmetry blocking errors were detected.
c       >  0 for iostat error.
c
c  31-oct-90 sifstv() created from sifrsh(). -rls
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
      real*8   score,  tcore,  vcore
      real*8   buffer(*),      values(*),      stv(nnbft,3)
c
c     # local...
      integer    itypea,   btypmx
      parameter( itypea=0, btypmx=6 )
c
      integer i,      nntot,  isym,   nrec,   last,   lastb,  lasta
      integer symtot(36), idummy(1), btypes(0:btypmx)
      real*8  fcore(3)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
c
c     # set the btypes(*) array:
c     #            0:s1, 1:t1, 2:v1, 3:vec, 4:vfc, 5:vref, 6:generic_h1
      data btypes/ 1,    2,    3,    3,     3,     3,      3          /
c
c     # note that generic_h1(*) contributions are added to the total 1-e
c     # potential, and any corresponding fcore contribution is included
c     # into vcore.  This may not be best in all cases, but these
c     # terms must be included somewhere.  -rls
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
      call wzero( nntot, stv(1,1), 1 )
      call wzero( nntot, stv(1,2), 1 )
      call wzero( nntot, stv(1,3), 1 )
      fcore(1) = (0)
      fcore(2) = (0)
      fcore(3) = (0)
c
      call sifr1n(
     & aoints, info,   itypea, btypmx,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  idummy, mapin,
     & nnbft,  stv,    fcore,  symb,
     & symtot, lasta,  lastb,  last,
     & nrec,   ierr )
c
c     # save the appropriate core values.
      score = fcore(1)
      tcore = fcore(2)
      vcore = fcore(3)
c
      return
      end
