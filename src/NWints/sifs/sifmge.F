      subroutine sifmge( nenrgy, energy, ietype, nnew, enew, typnew )
c
c  merge two energy(*) vectors such that exact duplicates are not
c  repeated.
c
c  input: nenrgy    = number of current energy(*) entries.
c         energy(*) = current energy entries. (distinct)
c         ietype(*) = energy types.
c         nnew      = number of new candidates.
c         enew(*)   = new energy values. (distinct)
c         typnew(*) = new types.
c
c  output: nenrgy    = updated value. may be as large as
c                      (original_nenrgy + nnew).
c          energy(1:nenrgy) = updated values.
c          ietype(1:nenrgy) = updated types..
c
c  this routine should be called by programs that write integral
c  files and which read energy(*) arrays from several sources.  If
c  input energy(*) values from two different sources have some common
c  ancestor, then duplication of common energy(*) elements should be
c  eliminated.  this is not only desirable for total energ values, but
c  it is necessary for frozen-core energy contributions in order to
c  avoid overcounting.
c
c  note that tests for exact equality are used in the comparisons.
c  this means that nothing should be done to input energy(*) values
c  which are passed on to other programs other than copy operations.
c  in particular, the programmer should avoid adding and then
c  subtracting the nuclear repulsion energy from total energy values.
c  it is the original, unmodified, values which should be transferred.
c
c  18-oct-90 written by ron shepard.
c
      implicit none
c
      integer  nenrgy, ietype(*),      nnew,   typnew(*)
      real*8   energy(*),      enew(*)
c
      integer  nold,   inew,   iold
c
c     # current energy(*) values are assumed to be distinct (i.e. no
c     # repetitions of the same energy contribution) among themselves.
c
c     # new contributions must be distinct among themselves.
c
      nold = nenrgy
c
      do 20 inew = 1, nnew
         do 10 iold = 1, nold
            if ( (energy(iold) .eq. enew(inew)) .and.
     &       ( ietype(iold) .eq. typnew(inew)) ) goto 20
10       continue
c        # loop exit means a distinct new energy has been found.
         nenrgy = nenrgy + 1
         energy(nenrgy) = enew(inew)
         ietype(nenrgy) = typnew(inew)
20    continue
c
      return
      end
