      subroutine build_maps(basis, cntoce, cntobfr, cetobfr, natoms, 
     &                      nshells)
      implicit none
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
     &        call errquit('build_maps: bad basis', 0)
      end do
c
      do ish = 1, nshells
         if (.not. bas_cn2bfr(basis,ish,cntobfr(1,ish),cntobfr(2,ish)))
     &        call errquit('build_maps: bad basis', 0)
      end do
c
      do iat = 1, natoms
         if (.not. bas_ce2bfr(basis,iat,cetobfr(1,iat),cetobfr(2,iat)))
     &        call errquit('build_maps: bad basis', 0)
      end do
c
      return
      end
      subroutine util_ga_mat_reduce(n, nr, map, g_a, n_a, r, op,
     &                              scr, mxmap)
      implicit none
c     
      integer n                         ! Original size [input]
      integer nr                        ! Reduced size  [input]
      integer map(2,nr)                 ! map(1,*)=lo, map(2,*)=hi [input]
      integer n_a                       ! number of GAs to reduce
      integer g_a(n_a)                  ! Original GA handle(s) [input]
      double precision r(nr,nr,n_a)     ! Reduced matrix [output]
      integer mxmap                     ! max size of map vector
      double precision scr(mxmap,mxmap) ! scratch
      character*(*) op                  ! Reduction operation
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax', 'rms' (extend as necessary)
c     
      integer ir, jr, i, j, k, ielems, jelems
      double precision sum
c     
      do k = 1, n_a
c
c         write(6,*) ' util_mat_reduce: input matrix '
c         call ga_print(g_a(k))
c     
         do jr = 1, nr
            do ir = 1, nr
               sum = 0.0d0
               ielems = map(2,ir) - map(1,ir) + 1
               jelems = map(2,jr) - map(1,jr) + 1
               if (ielems.gt.0 .and. jelems.gt.0) then
                  call ga_get(g_a(k), map(1,ir), map(2,ir), map(1,jr), 
     &                 map(2,jr), scr, mxmap)
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
                     call errquit('util_ga_mat_reduce: unknown op', 0)
                  end if
                  r(ir,jr,k) = max(r(ir,jr,k),sum)
               endif
            end do
         end do
c         write(6,*) ' util_mat_reduce: reduced matrix '
c         call output(r(1,1,k), 1, nr, 1, nr, nr, nr, 1)
      enddo
c     
      end
      subroutine util_irreg_mat_reduce(n_row, n_col, nr_row, nr_col, 
     &                                 row_map, col_map, a, r, op)
      implicit none
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
               call errquit('util_irreg_mat_reduce: unknown op', 0)
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
      subroutine util_rmsatbf(nqsmall, natoms,iniz,ifin,
     .     chi,rchi)
      implicit none
      integer nqsmall ! [input]
      integer natoms  ! [input]
      integer iniz(*) ! [input]
      integer ifin(*) ! [input]
      double precision chi(nqsmall,*) ! [input]
      double precision rchi(*) ! [output]
      
c
      integer q,iat,jf
      double precision sum,mxsum
c
      do iat=1,natoms
        if(iniz(iat).eq.0) then
          rchi(iat)=0d0
        else
          mxsum=0d0
          do jf=iniz(iat),ifin(iat)
            sum=0d0
            do q=1,nqsmall
              sum=sum+chi(q,jf)*chi(q,jf)
            enddo
            mxsum=max(mxsum,sqrt(sum))
          enddo
          rchi(iat)=mxsum
        endif
      enddo
      return
      end
