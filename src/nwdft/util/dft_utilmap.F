C> \ingroup nwdft_xc
C> @{
C>
C> \file dft_utilmap.F
C> Various routines that build or use mappings between atoms, matrix
C> elements and what not.
C>
      subroutine build_maps(basis, cntoce, cntobfr, cetobfr, natoms, 
     &                      nshells)
*
* $Id$
*
      implicit none
#include "errquit.fh"
      integer basis, natoms, nshells
      integer cntoce(nshells), cntobfr(2,nshells), cetobfr(2,natoms)
c      
#include "bas.fh"
c
      integer ish, iat
c     
c     Build maps (for speed).
c     
      do ish = 1, nshells
         if (.not. bas_cn2ce(basis, ish, cntoce(ish)))
     &        call errquit('build_maps: bad basis', 0, BASIS_ERR)
      end do
c
      do ish = 1, nshells
         if (.not. bas_cn2bfr(basis,ish,cntobfr(1,ish),cntobfr(2,ish)))
     &        call errquit('build_maps: bad basis', 0, BASIS_ERR)
      end do
c
      do iat = 1, natoms
         if (.not. bas_ce2bfr(basis,iat,cetobfr(1,iat),cetobfr(2,iat)))
     &        call errquit('build_maps: bad basis', 0, BASIS_ERR)
      end do
c
      return
      end
      subroutine util_ga_mat_reduce(n, nr, map, g_a, n_a, r, op,
     &                              scr, mxmap,lfirst)
      implicit none
c     
#include "global.fh"
#include "msgids.fh"
      integer n                         ! Original size [input]
      integer nr                        ! Reduced size  [input]
      integer map(2,nr)                 ! map(1,*)=lo, map(2,*)=hi [input]
      integer n_a                       ! number of GAs to reduce
      integer g_a(n_a)                  ! Original GA handle(s) [input]
      double precision r(nr,nr,n_a)     ! Reduced matrix [output]
      integer mxmap                     ! max size of map vector
      double precision scr(mxmap,*) ! scratch
      character*(*) op                  ! Reduction operation
      logical lfirst
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax', 'rms' (extend as necessary)
c     
      integer ir, jr, k, ielems, jelems
      double precision sum
      integer me,nproc,itask
      double precision util_matops
      external util_matops
      integer ityp, ndim, dims(3)
      integer lo(3),hi(3),ld(2),ld12,offs
      integer distr_lo(3),distr_hi(3)
      integer get_lo(3),get_hi(3)
      logical l3d,getdone
c     
      call nga_inquire(g_a, ityp, ndim, dims)
      l3d=.false.
      if(ndim.gt.2) then
#ifdef DEBUG
         write(6,*) ' matreduce: ndim gt 2!, =',ndim
         write(6,*) ' na ',n_a, ' dims1 ',dims(1),dims(1)/n_a
         write(6,*) '  dims ',dims
#endif
         l3d=.true.
         offs=dims(1)/n_a
      endif
      me=ga_nodeid()
      nproc=ga_nnodes()
      itask=0
      if(lfirst)
     A     call dfill(n_a*nr*nr, 0.0d0, r, 1)
      call nga_distribution(g_a,
     .           ga_nodeid(), distr_lo, distr_hi)
      if (distr_lo(1).gt.0 .and. distr_lo(1).le.distr_hi(1)) then

      do k = 1, n_a
c
c         write(6,*) ' util_mat_reduce: input matrix '
c         call ga_print(g_a(k))
c     
         do jr = 1, nr
cedo            do ir = 1, nr
            do ir = 1, jr
c               itask=itask+1
c               if(mod(itask,nproc).eq.me) then
                  ielems = map(2,ir) - map(1,ir) + 1
                  jelems = map(2,jr) - map(1,jr) + 1
                  if (ielems.gt.0 .and. jelems.gt.0) then
                     getdone=.false.
                     if(l3d) then
                        ielems=ielems*offs
                        lo(1)=1+(k-1)*offs
                        hi(1)=lo(1)+offs-1
#if 0
                        lo(2)=map(1,ir)
                        hi(2)=map(2,ir)
                        lo(3)=map(1,jr)
                        hi(3)=map(2,jr)
#else
                        lo(2)=max(map(1,ir),distr_lo(1))
                        hi(2)=min(map(2,ir),distr_hi(1))
                        lo(3)=max(map(1,jr),distr_lo(2))
                        hi(3)=min(map(2,jr),distr_hi(2))
#endif
                        ld(1)=offs
                        ld(2)=mxmap
                        getdone=(get_lo(2).le.get_hi(2)).and.
     A                       (get_lo(3).le.get_hi(3))
                        if(getdone)
     G                   call nga_get(g_a, lo,hi,
     .                       scr, ld)
                        ld12=ld(1)*ld(2)
                     else
c                        call ga_get(g_a(k), map(1,ir), map(2,ir), 
c     .                       map(1,jr), map(2,jr), scr, mxmap)
                        get_lo(1)=max(map(1,ir),distr_lo(1))
                        get_hi(1)=min(map(2,ir),distr_hi(1))
                        get_lo(2)=max(map(1,jr),distr_lo(2))
                        get_hi(2)=min(map(2,jr),distr_hi(2))
                        getdone=(get_lo(1).le.get_hi(1)).and.
     A                       (get_lo(2).le.get_hi(2))
                        if(getdone)
     G                  call ga_get(g_a(k), get_lo(1), get_hi(1), 
     .                       get_lo(2), get_hi(2), scr, mxmap)
                        ld12=mxmap
                     endif
                     if(getdone) then
                     sum=util_matops(op,ielems,jelems,scr,ld12)
                     r(ir,jr,k) = max(r(ir,jr,k),sum)
c     
c     copy upper triangle of r(ij) to upper triangle
c     
                     if(ir.ne.jr) r(jr,ir,k)=r(ir,jr,k)
                     endif
                  endif
cc               endif
            end do
         end do
c         write(6,*) ' util_mat_reduce: reduced matrix '
c         call output(r(1,1,k), 1, nr, 1, nr, nr, nr, 1)
      enddo
      endif
c
c     global sum
c
c      write(6,*) ' length ',nr*nr*n_a
      call ga_dgop(msg_gop_rdens,r,nr*nr*n_a, 'absmax')
c     
      end
      double precision function util_matops(op,ielems,jelems,scr,ld)
      implicit none
#include "errquit.fh"
      character*(*) op
      integer ielems,jelems
      integer ld
      double precision scr(ld,*)
c
      integer i,j
      double precision sum
      sum=0d0
      if (op .eq. 'abssum') then
         do j = 1, jelems
            do i = 1, ielems
               sum = sum + abs(scr(i,j))
            end do
         end do
      else if (op .eq. 'absmax') then
         do j = 1, jelems
            do i = 1, ielems
               sum = max(sum, abs(scr(i,j)))
            end do
         end do
      else if (op .eq. 'rms') then
         do j = 1, jelems
            do i = 1, ielems
               sum = sum + scr(i,j)*scr(i,j)
            end do
         enddo
         sum = sqrt(sum)
      else
         call errquit('util_ga_mat_red: unknown op',0, UNKNOWN_ERR)
      end if
      util_matops=sum
      return
      end
      subroutine util_irreg_mat_reduce(n_row, n_col, nr_row, nr_col, 
     &                                 row_map, col_map, a, r, op)
      implicit none
#include "errquit.fh"
c     
      integer n_row                      ! Original row size [input]
      integer n_col                      ! Original col size [input]
      integer nr_row                     ! Reduced row size  [input]
      integer nr_col                     ! Reduced col size  [input]
      integer row_map(2,nr_row)          ! map(1,*)=lo, map(2,*)=hi [input]
      integer col_map(2,nr_col)          ! map(1,*)=lo, map(2,*)=hi [input]
      double precision a(n_row, n_col)   ! Original matrix [input]
      double precision r(nr_row, nr_col) ! Reduced matrix [output]
      character*(*) op                   ! Reduction operation
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax', 'rms' (extend as necessary)
c     
      integer ir, jr, i, j
      double precision sum
c     
      do jr = 1, nr_col
         do ir = 1, nr_row
            sum = 0.0d0
            if (op .eq. 'abssum') then
               do j = col_map(1,jr), col_map(2,jr)
                  do i = row_map(1,ir), row_map(2,ir)
                     if (i.ne.0.and.j.ne.0)sum = sum + abs(a(i,j))
                  end do
               end do
            else if (op .eq. 'absmax') then
               do j = col_map(1,jr), col_map(2,jr)
                  do i = row_map(1,ir), row_map(2,ir)
                     if (i.ne.0.and.j.ne.0)sum = max(sum, abs(a(i,j)))
                  end do
               end do
            else if (op .eq. 'rms') then
               do j = col_map(1,jr), col_map(2,jr)
                  do i = row_map(1,ir), row_map(2,ir)
                     if (i.ne.0.and.j.ne.0)sum = sum + a(i,j)*a(i,j)
                  end do
               enddo
               sum = sqrt(sum)
            else
               call errquit('util_irreg_mat_reduce: unknown op', 0,
     &       UNKNOWN_ERR)
            end if
            r(ir,jr) = max(r(ir,jr),sum)
         end do
      end do
c     
c      write(6,*) ' util_irreg_mat_reduce: input matrix '
c      call output(a, 1, n_row, 1, n_col, n_row, n_col, 1)
c
c      write(6,*) ' Row map begin: ',(row_map(1,ir),ir = 1,nr_row)
c      write(6,*) ' Row map end: ',(row_map(2,ir),ir = 1,nr_row)
c      write(6,*) ' Col map begin: ',(col_map(1,ir),ir = 1,nr_col)
c      write(6,*) ' Col map end: ',(col_map(2,ir),ir = 1,nr_col)
c     
c      write(6,*) ' util_irreg_mat_reduce: reduced matrix '
c      call output(r, 1, nr_row, 1, nr_col, nr_row, nr_col, 1)
c
      end
C>
C> \brief Compute the root-mean-square value of the basis functions
C> or basis function gradients for each atom
C>
C> Compute the root-mean-square value of basis functions or basis
C> function gradients for each atom. This data is usefull for screening
C> purposes and is a measure for the size of the basis functions on
C> a particular atom at the grid points.
C>
C> The exact quantity computed depends on the inputs for `chi` and
C> `nqsmall`. If `chi` represents basis functions then `rchi` will
C> be computed to be
C> \f{eqnarray*}{
C>   \mathrm{rchi}(i_{at})
C>   &=& \sqrt{\sum_{j\in\{\chi(i_{at})\}}\sum_q\chi_j^2(r_q)}
C> \f}
C> If `chi` represents the gradients of the basis functions then `rchi`
C> will be computed to be (the number of components in the gradient is
C> in practice encoded in `nqsmall`)
C> \f{eqnarray*}{
C>   \mathrm{rchi}(i_{at})
C>   &=& \sqrt{\sum_{j\in\{\chi(i_{at})\}}\sum_q\sum_{c\in\{x,y,z\}}
C>             \nabla_c\chi_j^2(r_q)}
C> \f}
C> where \f$\{\chi(i_{at})\}\f$ designates the set of all basis
C> functions on atom \f$i_{at}\f$.
C>
      subroutine util_rmsatbf(nqsmall, natoms,iniz,ifin,
     .     chi,rchi)
      implicit none
      integer nqsmall !< [Input] The number of data points 
                      !< - the number of grid points, `nq`, for basis
                      !< functions
                      !< - `3*nq` for gradients of basis functions
      integer natoms  !< [Input] The number of atoms
      integer iniz(*) !< [Input] The initial basis function on each
                      !< atom
      integer ifin(*) !< [Input] The final basis function on each atom
      double precision chi(nqsmall,*) !< [Input] The basis function
                                      !< (gradient) amplitudes
      double precision rchi(*) !< [Output] The root-mean-square of 
                               !< of all amplitudes on each atom
c
      integer q,iat,jf
      double precision sum,mxsum
      double precision,external :: ddot
c
      mxsum=0d0
      do iat=1,natoms
        if(iniz(iat).eq.0) then
          rchi(iat)=0d0
        else
           sum=0d0
           do jf=iniz(iat),ifin(iat)
              do q=1,nqsmall
                 sum=sum+chi(q,jf)*chi(q,jf)
              enddo
           enddo
           rchi(iat)=max(mxsum,sqrt(sum))
        endif
      enddo
      return
      end
C>
C> @}
