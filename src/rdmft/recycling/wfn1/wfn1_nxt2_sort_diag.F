C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Sort the diagonal elements of the Fock matrix in 
C> natural orbital basis
C>
C> This routine sorts the diagonal elements of the Fock matrix in
C> the natural orbital basis in increasing order. When two elements
C> on the diagonal need to be interchanged this means 
C> 1. Swapping the Fock matrix rows
C> 2. Swapping the Fock matrix columns
C> 3. Swapping the natural orbital columns
C> 4. Swapping the occupation function rows
C>
      subroutine wfn1_nxt2_sort_diag(nbf,f,n,o)
      implicit none
      integer nbf                 !< [Input] The number of basis
                                  !< functions
      double precision f(nbf,nbf) !< [In/Output] The Fock matrix
      double precision n(nbf,nbf) !< [In/Output] The natural orbitals
      double precision o(nbf,nbf) !< [In/Output] The occupation functions
c
      integer ir  !< Counter over rows
      integer ic  !< Counter over columns
      integer icl !< Counter over columns local to a degenerate set
      integer is  !< The selected columns
      integer nd  !< The number of degenerate vectors
      integer ii  !< Counter
      integer jj  !< Counter
      integer id  !< Diagonal element of minimal value
c
      double precision ds !< The value of the selected element
      double precision dt !< Temporary
c
      do ii = 1, nbf-1
        ds = f(ii,ii)
        id = ii
        do jj = ii+1, nbf
          if (f(jj,jj).lt.ds) then
            id = jj
            ds = f(jj,jj)
          endif
        enddo
c
        if (ds.lt.f(ii,ii)) then
c
c         Swap elements ii and id
c
c         1. Swap rows of the Fock matrix
c
          do jj = 1, nbf
            dt = f(ii,jj)
            f(ii,jj) = f(id,jj)
            f(id,jj) = dt
          enddo
c
c         2. Swap columns of the Fock matrix
c
          do jj = 1, nbf
            dt = f(jj,ii)
            f(jj,ii) = f(jj,id)
            f(jj,id) = dt
          enddo
c
c         3. Swap columns of the natural orbitals
c
          do jj = 1, nbf
            dt = n(jj,ii)
            n(jj,ii) = n(jj,id)
            n(jj,id) = dt
          enddo
c
c         4. Swap rows of the occupation functions
c
          do jj = 1, nbf
            dt = o(ii,jj)
            o(ii,jj) = o(id,jj)
            o(id,jj) = dt
          enddo
        endif
      enddo
cDEBUG
c     write(*,*)'eigenvalues'
c     do ii = 1, nbf
c       write(*,*)ii,e(ii)
c     enddo
c     write(*,*)'eigenvectors'
c     call hess_hssout(u,nbf,nbf,nbf)
cDEBUG
c
      end
C>
C> @}
