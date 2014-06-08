C> \ingroup util_ga
C> @{
C>
C> \brief Find the element with the absolute maximum value in 
C> a global array
C>
C> The subroutine finds the element with the absolute maximum
C> value in a global array and returns the position of that 
C> element as well as the absolute maximum value itself.
C> Internally it is assumed that the array to hold
C> the element position is `ga_max_dim` long. However, only
C> the first `ga_ndim(g_a)` elements will be accessed.
C>
      subroutine nga_idamax(g_a, indx, valmx)
C     $Id$
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "errquit.fh"
      integer g_a              !< [Input] The global array
      integer indx(ga_max_dim) !< [Output] The coordinates of the
                               !< array element with maximum value
      double precision valmx   !< [Output] The absolute maximum value
c     
      integer ilo(ga_max_dim), ihi(ga_max_dim)
      integer itype, ndim, dims(ga_max_dim)
      integer id, nproc, iproc, ii
      logical doit
      integer l_indx, k_indx ! local copies for all processors
      integer l_val , k_val  ! local copies for all processors
      integer nga_idamax_dbl_doit, nga_idamax_dcpl_doit
      external nga_idamax_dbl_doit, nga_idamax_dcpl_doit
c     
c     Each process first determines max of data local to it
c     
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      call ga_sync()
      call nga_inquire(g_a,itype,ndim,dims)
      call nga_distribution(g_a, iproc, ilo, ihi)
      do id = ndim+1, ga_max_dim
        ilo(id) = 0
        ihi(id) = 0
      enddo
c
c     Get some memory
c
      if (.not.ma_push_get(MT_INT,ga_max_dim*nproc,"indx",
     &                     l_indx,k_indx))
     &  call errquit("nga_idamax: allocate indx failed",
     &               ma_sizeof(MT_INT,ga_max_dim*nproc,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nproc,"val",l_val,k_val))
     &  call errquit("nga_idamax: allocate val failed",
     &               ma_sizeof(MT_DBL,nproc,MT_BYTE),MA_ERR)
c
      call ifill(ga_max_dim*nproc,0,int_mb(k_indx),1)
      call dfill(nproc,0.0d0,dbl_mb(k_val),1)
c
      doit = .true.
      do id = 1, ndim
        doit = doit .and. (ihi(id).gt.0) .and. (ihi(id).ge.ilo(id))
      enddo
      if (doit) then
        if (itype.eq.MT_DBL) then
          call nga_access_callback_release(
     $         g_a, ilo, ihi, nga_idamax_dbl_doit, ndim,
     $         int_mb(k_indx+iproc*ga_max_dim), dbl_mb(k_val+iproc),
     $         0, 0, 0,0)
        else if (itype.eq.MT_DCPL) then
          call nga_access_callback_release(
     $         g_a, ilo, ihi, nga_idamax_dcpl_doit, ndim,
     $         int_mb(k_indx+iproc*ga_max_dim), dbl_mb(k_val+iproc),
     $         0, 0, 0,0)
        else
          call errquit("nga_idamax: invalid data type",itype,UERR)
        endif
      endif
c     
c     Finally a global sum via message passing
c     
      call ga_sync()
      call ga_dgop(Msg_MxElt, dbl_mb(k_val), nproc, '+')
      call ga_igop(Msg_MxElt, int_mb(k_indx), ga_max_dim*nproc, '+')
      call ga_sync()
c
c     Now each processor finds the maximum value and the corresponding
c     coordinates
c     
      valmx = dbl_mb(k_val)
      do ii = 1, ndim
        indx(ii) = int_mb(k_indx+ii-1)
      enddo
      do id = 2, nproc
        if (dbl_mb(k_val+id-1).gt.valmx) then
          valmx = dbl_mb(k_val+id-1)
          do ii = 1, ndim
            indx(ii) = int_mb(k_indx+(id-1)*ga_max_dim+ii-1)
          enddo
        endif
      enddo
c
c     Clean up memory and return
c
      if (.not.ma_pop_stack(l_val))
     &  call errquit("nga_idamax: deallocate val failed",0,MA_ERR)
      if (.not.ma_pop_stack(l_indx))
     &  call errquit("nga_idamax: allocate indx failed",0,MA_ERR)
c
      end
C>
C> \brief Find the absolute maximum value in a block of a double 
C> precision GA
C>
C> Finds and returns the absolute maximum value in a block of a 
C> global array and returns this value as well as its position in
C> the global array.
C>
      integer function nga_idamax_dbl_doit(
     $     g_a, ilo, ihi, a, ld, ndim, indx, value, 
     $     idum4, idum5, idum6, idum7)
      implicit none
c
c     Double precision version
c
#include "global.fh"
#include "errquit.fh"
c
      integer g_a        !< [Input] The global array handle (not used)
      integer ndim       !< [Input] The number of dimensions in the GA
      integer ilo(ga_max_dim)  !< [Input] The lower GA block limits
      integer ihi(ga_max_dim)  !< [Input] The upper GA block limits
      integer ld(ga_max_dim)   !< [Input] The leading dimensions of
                               !< the block
      integer indx(ga_max_dim) !< [Output] The coordinates of the
                               !< maximum
c
      double precision a(*)  !< [Input] The block of data
      double precision value !< [Output] The absolute maximum value
c
      integer idum4  !< [Dummy]
      integer idum5  !< [Dummy]
      integer idum6  !< [Dummy]
      integer idum7  !< [Dummy]
c
      integer i1, i2, i3, i4, i5, i6, i7
      integer indx1, indx2, indx3, indx4, indx5, indx6, indx7
      integer ld1, ld2, ld3, ld4, ld5, ld6, ld7
c
      if (ndim.gt.7)
     &  call errquit("nga_idamax_dbl_doit: ndim out of range",
     &               ndim, UERR)
c
      ld1 = 1
      ld2 = ld1*ld(1)
      ld3 = ld2*ld(2)
      ld4 = ld3*ld(3)
      ld5 = ld4*ld(4)
      ld6 = ld5*ld(5)
      ld7 = ld6*ld(6)
      value = 0.0d0
      do i7 = 1, (ihi(7)-ilo(7)+1)
        indx7 = (i7-1)*ld7
        do i6 = 1, (ihi(6)-ilo(6)+1)
          indx6 = indx7 + (i6-1)*ld6
          do i5 = 1, (ihi(5)-ilo(5)+1)
            indx5 = indx6 + (i5-1)*ld5
            do i4 = 1, (ihi(4)-ilo(4)+1)
              indx4 = indx5 + (i4-1)*ld4
              do i3 = 1, (ihi(3)-ilo(3)+1)
                indx3 = indx4 + (i3-1)*ld3
                do i2 = 1, (ihi(2)-ilo(2)+1)
                  indx2 = indx3 + (i2-1)*ld2
                  do i1 = 1, (ihi(1)-ilo(1)+1)
                    indx1 = indx2 + i1
                    if (abs(a(indx1)).gt.value) then
                      value = abs(a(indx1))
                      indx(1) = ilo(1) + i1 - 1
                      indx(2) = ilo(2) + i2 - 1
                      indx(3) = ilo(3) + i3 - 1
                      indx(4) = ilo(4) + i4 - 1
                      indx(5) = ilo(5) + i5 - 1
                      indx(6) = ilo(6) + i6 - 1
                      indx(7) = ilo(7) + i7 - 1
                    endif
                  end do ! i1
                end do ! i2
              end do ! i3
            end do ! i4
          end do ! i5
        end do ! i6
      end do ! i7
c
      nga_idamax_dbl_doit = 0    ! no update necessary
c
      end
C>
C> \brief Find the absolute maximum value in a block of a double 
C> complex GA
C>
C> Finds and returns the absolute maximum value in a block of a 
C> global array and returns this value as well as its position in
C> the global array.
C>
      integer function nga_idamax_dcpl_doit(
     $     g_a, ilo, ihi, a, ld, ndim, indx, value, 
     $     idum4, idum5, idum6, idum7)
      implicit none
c
c     Double complex version
c
#include "global.fh"
#include "errquit.fh"
c
      integer g_a        !< [Input] The global array handle (not used)
      integer ndim       !< [Input] The number of dimensions in the GA
      integer ilo(ga_max_dim)  !< [Input] The lower GA block limits
      integer ihi(ga_max_dim)  !< [Input] The upper GA block limits
      integer ld(ga_max_dim)   !< [Input] The leading dimensions of
                               !< the block
      integer indx(ga_max_dim) !< [Output] The coordinates of the
                               !< maximum
c
      double complex   a(*)  !< [Input] The block of data
      double precision value !< [Output] The absolute maximum value
c
      integer idum4  !< [Dummy]
      integer idum5  !< [Dummy]
      integer idum6  !< [Dummy]
      integer idum7  !< [Dummy]
c
      integer i1, i2, i3, i4, i5, i6, i7
      integer indx1, indx2, indx3, indx4, indx5, indx6, indx7
      integer ld1, ld2, ld3, ld4, ld5, ld6, ld7
c
      if (ndim.gt.7)
     &  call errquit("nga_idamax_dbl_doit: ndim out of range",
     &               ndim, UERR)
c
      ld1 = 1
      ld2 = ld1*ld(1)
      ld3 = ld2*ld(2)
      ld4 = ld3*ld(3)
      ld5 = ld4*ld(4)
      ld6 = ld5*ld(5)
      ld7 = ld6*ld(6)
      value = 0.0d0
      do i7 = 1, (ihi(7)-ilo(7)+1)
        indx7 = (i7-1)*ld7
        do i6 = 1, (ihi(6)-ilo(6)+1)
          indx6 = indx7 + (i6-1)*ld6
          do i5 = 1, (ihi(5)-ilo(5)+1)
            indx5 = indx6 + (i5-1)*ld5
            do i4 = 1, (ihi(4)-ilo(4)+1)
              indx4 = indx5 + (i4-1)*ld4
              do i3 = 1, (ihi(3)-ilo(3)+1)
                indx3 = indx4 + (i3-1)*ld3
                do i2 = 1, (ihi(2)-ilo(2)+1)
                  indx2 = indx3 + (i2-1)*ld2
                  do i1 = 1, (ihi(1)-ilo(1)+1)
                    indx1 = indx2 + i1
                    if (abs(a(indx1)).gt.value) then
                      value = abs(a(indx1))
                      indx(1) = ilo(1) + i1 - 1
                      indx(2) = ilo(2) + i2 - 1
                      indx(3) = ilo(3) + i3 - 1
                      indx(4) = ilo(4) + i4 - 1
                      indx(5) = ilo(5) + i5 - 1
                      indx(6) = ilo(6) + i6 - 1
                      indx(7) = ilo(7) + i7 - 1
                    endif
                  end do ! i1
                end do ! i2
              end do ! i3
            end do ! i4
          end do ! i5
        end do ! i6
      end do ! i7
c
      nga_idamax_dcpl_doit = 0    ! no update necessary
c
      end
C>
C> @}
