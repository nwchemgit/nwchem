C> \ingroup wfn1
C> @{
C>
C> \brief Swap two vectors
C>
C> Swap two vectors (columns) in a matrix. The columns are indicated
C> by `io` for the occupied orbital and `iv` for the virtual orbital.
C> If `io` and `iv` are the same (i.e. `ne`) then the identity
C> permutation is applied (i.e. do nothing).
C>
      subroutine wfn1_swap_vectors(ne,nbf,io,iv,a)
      implicit none
c
#include "errquit.fh"
c
      integer ne  !< [Input] The number of electrons, i.e. occupied
                  !< orbitals
      integer nbf !< [Input] The number of orbitals
      integer io  !< [Input] The occupied orbital number
      integer iv  !< [Input] The virtual orbital number
c
      double precision a(nbf,nbf) !< [In/Output] On input the matrix,
                                  !< on output the same matrix but
                                  !< two columns swapped
      double precision tmp
c
      integer no !< Number of occupied orbitals
      integer nv !< Number of virtual orbitals
      integer ii !< Counter
c
      if (ne.eq.0.and.io.eq.0.and.iv.eq.0) then
c
c       As there are no electrons no sensible permutation is can be
c       formulated. The above combination is the only sensible thing
c       that can be asked for, i.e. perform no permutation between
c       occupied and unoccupied orbitals as there are no 
c       occupied states. So just return.
c
        return
      endif
      if (io.lt.1) then
        call errquit("wfn1_swap_vectors: occupied orbital number "
     &             //"too low",io,UERR)
      endif
      if (io.gt.ne) then
        call errquit("wfn1_swap_vectors: occupied orbital number "
     &             //"too high",io,UERR)
      endif
      if (iv.lt.ne) then
        call errquit("wfn1_swap_vectors: virtual orbital number "
     &             //"too low",iv,UERR)
      endif
      if (iv.gt.nbf) then
        call errquit("wfn1_swap_vectors: virtual orbital number "
     &             //"too high",iv,UERR)
      endif
      if (iv.eq.ne) then
        if (io.ne.ne) then
          call errquit("wfn1_swap_vectors: identity permutation "
     &               //"invalidly represented",io,UERR)
        endif
      endif
c
      no = ne
      nv = nbf - ne
      do ii = 1, nbf
        tmp      = a(ii,io)
        a(ii,io) = a(ii,iv)
        a(ii,iv) = tmp
      enddo
c
      end
C> @}
