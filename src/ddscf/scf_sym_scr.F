      subroutine scf_sym_screen(nclosed, nopen, nvir, irs,
     $                          g_v, ivoff, oprint, maxelt)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer nclosed, nopen, nvir ! [input] No. of MOs in each space
      integer irs(*) ! [input] symmetry of each MO
      integer g_v    ! [input/output] rotation vector to screen
      integer ivoff  ! [input] Offset within g_v to start of vector
      logical oprint ! [input] If true print warnings
      double precision maxelt ! [output] Returns maximum forbidden element
c
c     Screen the orbital rotation/gradient vector for symmetry
c     forbidden rotations.  Zero them and return the largest.
c
c     Vector is assumed to hold (nclosed+nopen)*nvir + nclosed*nopen
c
      integer i, a, ioff
      double precision test
      integer l_buf, k_buf
      double precision tmp(2)
      integer nmixed 
      integer i0,i1,j0,j1
      integer i_min,i_max,len_i
c
      if (.not. ma_push_get(mt_dbl, nclosed+nopen+nvir, 
     $     'scf_sym_screen', l_buf, k_buf)) 
     $     call errquit('scf_sym_screen: ma?',nclosed+nopen+nvir,
     &       MA_ERR)
c
      call ga_sync
c
      maxelt = 0.0d0
      nmixed = 0
      call ga_distribution(g_v,
     .     ga_nodeid(), i0, i1, j0, j1)
c      do i = 1+ga_nodeid(),nclosed+nopen,ga_nnodes()
      if (i0.gt.0 .and. i0.le.i1.and.j0.gt.0.and.j0.le.j1) then
      do i = 1,nclosed+nopen
        if ((nclosed+nopen)*nvir .gt. 0) then
          ioff = 1+(i-1)*nvir + ivoff
          if(.not.(
     A         (ioff+nvir-1.lt.i0).or.
     B         (ioff.gt.i1))) then
          i_min=max(i0,ioff)
          i_max=min(i1,ioff+nvir-1)
          len_i=i_max-i_min+1
          if(len_i.lt.1) then
             write(6,431) ioff,nvir,i0,i1,i_min,i_max
 431         format(' ioff=',i7,' nvir= ',i7,' i0= ',i7,' i1= ',i7,
     M            ' i_min=',i7,' i_max= ',i7)
             call errquit('0len_i neg ',len_i,0)
          endif
          call ga_get(g_v, i_min, i_max, 1, 1, dbl_mb(k_buf+i_min-ioff), 
     C         len_i)
          do a=i_min-ioff+1,i_max-ioff
            if (irs(i).ne.irs(a+nclosed+nopen)) then
               test = abs(dbl_mb(k_buf-1+a))
               maxelt = max(maxelt,test)
               dbl_mb(k_buf-1+a) = 0.0d0
               if (test.gt.1d-6) nmixed = nmixed + 1
            endif
          enddo
          call ga_put(g_v, i_min, i_max, 1, 1, dbl_mb(k_buf+i_min-ioff),
     C         len_i)
        endif 
        endif ! if nc+nop
c
        if (i.le.nclosed .and. nopen*nclosed.gt.0) then
          ioff = 1+(nclosed+nopen)*nvir+(i-1)*nopen + ivoff
          if(.not.(
     A         (ioff+nopen-1.lt.i0).or.
     B         (ioff.gt.i1))) then
          i_min=max(i0,ioff)
          i_max=min(i1,ioff+nopen-1)
          len_i=i_max-i_min+1
          if(len_i.lt.1) then
             write(6,432) ioff,nopen,i0,i1,i_min,i_max
 432         format(' ioff=',i7,' nopen= ',i7,' i0= ',i7,' i1= ',i7,
     M            ' i_min=',i7,' i_max= ',i7)
             call errquit('1len_i neg ',len_i,0)
          endif
           if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
          call ga_get(g_v, i_min, i_max, 1, 1, dbl_mb(k_buf+i_min-ioff),
     C         len_i)
           if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
c          do a=1,nopen
         do a=i_min-ioff+1,i_max-ioff
             if (irs(i).ne.irs(a+nclosed)) then
               test = abs(dbl_mb(k_buf-1+a))
               maxelt = max(maxelt,test)
               dbl_mb(k_buf-1+a) = 0.0d0
               if (test.gt.1d-2) nmixed = nmixed + 1
            endif
          enddo
           if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
          call ga_put(g_v, i_min, i_max, 1, 1, dbl_mb(k_buf+i_min-ioff),
     C         len_i)
           if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
c          call ga_put(g_v, ioff, ioff+nopen-1, 1, 1, dbl_mb(k_buf), 1)
        endif
        endif
      enddo ! do i=1,nc+no
      endif
c
      if (.not. ma_pop_stack(l_buf))
     $     call errquit('scf_sym_screen: ma?',-1, MA_ERR)
      call ga_sync
c
      tmp(1) = nmixed
      tmp(2) = maxelt
      call ga_dgop(88,tmp,2,'max')
      nmixed = tmp(1)
      maxelt = tmp(2)
      if (nmixed.gt.0 .and. oprint .and. ga_nodeid().eq.0) then
        write(6,1) nmixed, maxelt
 1      format(/' !!',i8,' symmetry forbidden rotations with max = ',
     $         1p,d9.2/)
        call util_flush(6)
      endif
c
      end
      subroutine scf_sym_mat_screen(g_m, nmo, irs)
      implicit none
#include "global.fh"
      integer g_m               ! [input/output] MO GA to be screened.
      integer nmo               ! [input] NO. of MOs
      integer irs(nmo)          ! [input] Symmetry of each MO
c
c     Force to zero elements of the MO array that should be zero
c     by symmetry.
c
c     Tested but no longer used.
c
      integer i, j
      integer ilo, ihi, jlo, jhi
c
      call ga_sync()
      call ga_distribution(g_m, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         do i = ilo, ihi
            do j = jlo, jhi
               if (irs(i).ne.irs(j)) 
     $              call ga_put(g_m, i, i, j, j, 0.0d0, 1)
            enddo
         enddo
      endif
      call ga_sync()
c
      end

      
