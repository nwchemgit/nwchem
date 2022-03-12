      Subroutine gw_invdiag(rfp, ld, n, ecrpa, packed)
C
C     Compute the inverse of a matrix, possibly packed in RFP format,
C     by calling the appropriate eigendecomposition subroutines.
C
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c
      logical packed
      integer ld, n
      double precision rfp(ld,*), ecrpa

      integer g_tmp, g_tmp2, g_tmp3, g_diag
      integer iev,lev,idum,ll,ul,me,i,j,idx

      double precision det
      double precision, parameter :: toll = 1d-8

      me = ga_nodeid()

      if (.not.ga_create(mt_dbl, n, n, 'ga_temp', n, 0, g_tmp))
     $  call errquit('gw_invdiag: error creating ga_temp',0,GA_ERR)
      if (.not.ga_create(mt_dbl, n, n, 'ga_temp2', n, 0, g_tmp2))
     $  call errquit('gw_invdiag: error creating ga_temp2',0,GA_ERR)
      if (.not.nga_create(mt_dbl,1,n,'diagonal',n,g_diag))
     $  call errquit('gw_invdiag: error creating ga_diag',0,GA_ERR)
      if(.not.MA_Push_Get(MT_Dbl,n,'evals',lev,iev))
     &  call errquit('dft_invdiag: cannot allocate evals',0, MA_ERR)

      call ga_distribution(g_tmp,me,idum,idum,ll,ul)
      call ga_zero(g_tmp)

      ecrpa = 0d0

      !
      ! If matrix is in RFP format, we first need to unpack it.
      !
      if (packed) then

        ! N is odd
        if (ld.eq.n) then
          do i=1,n/2+1
            if (i.ge.ll .and. i.le.ul) then
              rfp(i,i) = 0.5d0*rfp(i,i)
              call ga_put(g_tmp,i,n,i,i,rfp(i,i),1)
              rfp(i,i) = 2d0*rfp(i,i)
            endif
            if (i.gt.1) then
              j = i + n/2
              if (j.ge.ll .and. j.le.ul) then
                rfp(i-1,i) = 0.5d0*rfp(i-1,i)
                call ga_put(g_tmp,n/2+2,j,j,j,rfp(1,i),1)
                rfp(i-1,i) = 2d0*rfp(i-1,i)
              endif
            endif
          enddo

        ! N is even
        else
          do i=1,n/2
            j = i + n/2
            rfp(i:i+1,i) = 0.5d0*rfp(i:i+1,i)
            if (i.ge.ll .and. i.le.ul) then
              call ga_put(g_tmp,i,n,i,i,rfp(i+1,i),1)
            endif
            if (j.ge.ll .and. j.le.ul) then
              call ga_put(g_tmp,n/2+1,j,j,j,rfp(1,i),1)
            endif
            rfp(i:i+1,i) = 2d0*rfp(i:i+1,i)
          enddo
        endif

        call ga_symmetrize(g_tmp)
        call ga_scale(g_tmp,2d0)

      else

        call ga_put(g_tmp,1,n,ll,ul,rfp(1,ll),ld)

      endif
c
      !
      ! Compute the contribution to the RPA correlation energy
      !
      call ga_get_diag(g_tmp,g_diag)
      if (me.eq.0) then
        call nga_access(g_diag,1,n,idx,idum)
        do i=0,n-1
         ecrpa = ecrpa - dbl_mb(idx+i)
        enddo
        ecrpa = ecrpa + dble(n)
        call nga_release(g_diag,1,n)
      endif

c
      !
      ! Perform the diagonalization
      !
      call ga_sync

#ifdef SCALAPACK
      call dft_scaleig(g_tmp, g_tmp2, DBl_MB(iev),
     P     0) 
#else
      call ga_diag_std(g_tmp,g_tmp2,DBl_MB(iev))
#endif

      !
      ! Remove small eigenvalues
      !
      do i=0,n-1
        if (dbl_mb(iev+i).gt.0d0) ecrpa = ecrpa + log(dbl_mb(iev+i))
        if(abs(DBl_MB(iev+i)).lt.toll) then
          if(me.eq.0) write(LuOut,*) ' GAFACT - singular eigenvalue',i
          call util_flush(LuOut)
          DBl_MB(iev+i)=0.d0
        else
          DBl_MB(iev+i)=1.d0/DBl_MB(iev+i)
        endif 
      enddo

      !
      ! Obtain the (pseudo) inverse of the matrix
      !
      call ga_copy(g_tmp2,g_tmp)
      call ga_access(g_tmp,1,n,ll,ul,idx,idum)
      do i=ll,ul
        call dscal(n,dbl_mb(iev+i-1),dbl_mb(idx+(i-ll)*n),1)
      enddo
      call ga_release_update(g_tmp,1,n,ll,ul)

      if (.not.ga_create(mt_dbl, n, n, 'ga_temp3', n, 0, g_tmp3))
     $  call errquit('gw_invdiag: error creating ga_temp2',0,GA_ERR)

      call ga_dgemm('N','T',n,n,n,1.d0,g_tmp,g_tmp2,0.d0,g_tmp3)

      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('dft_invdiag: could not destroy g_tmp2', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit
     &   ('dft_invdiag: could not destroy g_tmp2', 0, GA_ERR)
      if(.not.ma_pop_stack(lev))
     &  call errquit('dft_invdiag: cannot pop stack',0, MA_ERR)


      !
      ! Go back to RFP format
      !
      if (packed) then
        ! N is odd
        if (ld.eq.n) then
          do i=1,n/2+1
            call ga_get(g_tmp3,i,n,i,i,rfp(i,i),1)
            if (i.gt.1) then
              j = i + n/2
              call ga_get(g_tmp3,n/2+2,j,j,j,rfp(1,i),1)
            endif
          enddo
        ! N is even
        else
          do i=1,n/2
            j = i + n/2
            call ga_get(g_tmp3,i,n,i,i,rfp(i+1,i),1)
            call ga_get(g_tmp3,n/2+1,j,j,j,rfp(1,i),1)
          enddo
        endif
      else
        call dfill(n**2,0d0,rfp,1)
        call ga_get(g_tmp3,1,n,ll,ul,rfp(1,ll),ld)
      endif

      if (.not. ga_destroy(g_tmp3)) call errquit
     &   ('dft_invdiag: could not destroy g_tmp3', 0, GA_ERR)
      if (.not. ga_destroy(g_diag)) call errquit
     &   ('dft_invdiag: could not destroy g_diag', 0, GA_ERR)
      return
      end
