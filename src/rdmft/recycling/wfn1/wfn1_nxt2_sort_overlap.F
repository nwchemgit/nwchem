C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Sort the overlap matrix
C>
C> Given the overlap matrix between the ground state orbitals and
C> the orbitals of an excited state we sort the columns. This problem
C> is tackled row by row:
C> 1. For each row i:
C>    a. Find j = iabsmax(row(i,i:n))
C>    b. Swap columns i and j of the overlap matrix
C>    c. Swap occupation functions i and j
C>
      subroutine wfn1_nxt2_sort_overlap(nbf,ov,o)
      implicit none
      integer nbf                  !< [Input] The number of basis
                                   !< functions
      double precision ov(nbf,nbf) !< [In/Output] The overlap matrix
      double precision o(nbf,nbf)  !< [In/Output] The occupation functions
c
      integer ii  !< Counter
      integer jj  !< Counter
      integer id  !< Diagonal element of minimal value
c
      double precision ds !< The value of the selected element
      double precision dt !< Temporary
c
      do ii = 1, nbf-1
        ds = abs(ov(ii,ii))
        id = ii
        do jj = ii+1, nbf
          if (abs(ov(ii,jj)).gt.ds) then
            id = jj
            ds = abs(ov(ii,jj))
          endif
        enddo
c
        if (ii.ne.id) then
c
c         Swap elements ii and id
c
c         1. Swap columns of the overlap matrix
c
          do jj = 1, nbf
            dt        = ov(jj,ii)
            ov(jj,ii) = ov(jj,id)
            ov(jj,id) = dt
          enddo
c
c         2. Swap columns of the occupation functions
c
          do jj = 1, nbf
            dt       = o(jj,ii)
            o(jj,ii) = o(jj,id)
            o(jj,id) = dt
          enddo
c
        endif
c
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
