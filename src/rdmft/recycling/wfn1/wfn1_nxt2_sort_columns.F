C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Sort the columns within a degenerate set of a unitary matrix
C> so that the maximum values appear on the diagonal
C>
C> Dsyev sorted the eigenvalues and vector so we need to reorder
C> them back into place. For every row find the absolute maximum
C> value then permute the columns to put that maximum value on the
C> diagonal.
C>
      subroutine wfn1_nxt2_sort_columns(nbf,u,e)
      implicit none
      integer nbf                 !< [Input] The number of basis
                                  !< functions
      double precision u(nbf,nbf) !< [In/Output] The (unitary) matrix
      double precision e(nbf)     !< [In/Output] The eigenvalues
c
      integer ir  !< Counter over rows
      integer ic  !< Counter over columns
      integer icl !< Counter over columns local to a degenerate set
      integer is  !< The selected columns
      integer nd  !< The number of degenerate vectors
      integer ii  !< Counter
c
      double precision ds !< The value of the selected element
      double precision dt !< Temporary
c
      ii = 1
      ic = 1
      do while (ic.lt.nbf)
        nd = 1
        do while((ic+nd.lt.nbf).and.(abs(e(ic+nd)-e(ic)).lt.0.001d0))
          nd = nd + 1
        enddo
        if (nd.gt.1) then

          do ir = ic, ic+nd-1
            is = ir
            ds = abs(u(ir,is))
            do icl = ir+1, nbf
              if (abs(u(ir,icl)).gt.ds) then
                is = icl
                ds = abs(u(ir,icl))
              endif
            enddo
            if (is.gt.ir) then
c
c             Swap the columns
c
              do ii = 1, nbf
                dt       = u(ii,ir)
                u(ii,ir) = u(ii,is)
                u(ii,is) = dt
              enddo
              dt    = e(ir)
              e(ir) = e(is)
              e(is) = dt
            endif
          enddo

        endif
        ic = ic + nd
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
