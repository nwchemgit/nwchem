      subroutine ga_orthog(g_vecs, g_over, ometric)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_vecs, g_over, ometric
c
      integer ulo, uhi, type, ni, nu, uchunk, ld
      integer l_tmp, k_tmp, k_s, l_s, l_over, k_over, l_w, k_w
c     
c     Redistribute the input matrix ... block the leading
c     dimension, leave second dimension undistributed
c     ... each process has a(ulo:uhi,1:ni).  If
c
      call ga_inquire(g_vecs, type, nu, ni)
      uchunk = max(64, nu/ga_nnodes())
      ulo  = ga_nodeid()*uchunk + 1
      uhi  = min(ulo + uchunk - 1, nu)
      ld   = uhi - ulo + 1
      if (ulo .gt. uhi) then
         ulo = 0
         uhi = -1
         ld  = 1
      end if
c
c     Allocate local workspace
c
      if (.not. ma_push_get(mt_dbl, ld*ni, 'ga_orthog: tmp', 
     $     l_tmp, k_tmp)) call errquit('ga_orthog: ma get failed', 
     $     ld*ni)
      if (.not. ma_push_get(mt_dbl, ld*nu, 'ga_orthog: over', 
     $     l_over, k_over)) call errquit('ga_orthog: ma get failed', 
     $     ld*ni)
      if (.not. ma_push_get(mt_dbl, ni, 'ga_orthog: sji', 
     $     l_s, k_s)) call errquit('ga_orthog: ma get failed', ni)
      if (.not. ma_push_get(mt_dbl, nu, 'ga_orthog: sji', 
     $     l_w, k_w)) call errquit('ga_orthog: ma get failed', ni)
c
      if (uhi .ge. ulo) then
         call ga_get(g_vecs, ulo, uhi, 1, ni, dbl_mb(k_tmp), ld)
         if (ometric) 
     $        call ga_get(g_over, 1, nu, ulo, uhi, dbl_mb(k_tmp), ld)
      endif
c
c     Do the deed
c
      call ga_sync()
      call ga_orthog_vecs(dbl_mb(k_tmp), ld, ni, ulo, uhi, dbl_mb(k_s),
     $     dbl_mb(k_over), nu, dbl_mb(k_w), ometric)
      call ga_sync()
c
c     Put results back
c
      if (uhi .ge. ulo)
     $     call ga_put(g_vecs, ulo, uhi, 1, ni, dbl_mb(k_tmp), ld)
c
c     Tidy up memory
c
      if (.not. ma_pop_stack(l_w)) call errquit
     $     ('ga_orthog: ma error popping w', 0)
      if (.not. ma_pop_stack(l_s)) call errquit
     $     ('ga_orthog: ma error popping s', 0)
      if (.not. ma_pop_stack(l_over)) call errquit
     $     ('ga_orthog: ma error popping over', 0)
      if (.not. ma_pop_stack(l_tmp)) call errquit
     $     ('ga_orthog: ma error popping tmp', 0)
c
      end
      subroutine ga_orthog_vecs(vecs, ld, ni, ulo, uhi, s, o, nu, w, 
     $     ometric)
      implicit none
c
      integer ld, ni, ulo, uhi, nu
      double precision 
     $     vecs(ulo:(ulo+ld-1),1:ni), ! [input/output] Local patch of vectors
     $     o(1:nu,ulo:(ulo+ld-1)),    ! [input] metric matrix if (ometric)
     $     s(ni),               ! [scratch] workspace for making overlap
     $     w(nu)                ! [scratch] workspace for making overlap
      logical ometric
c
      integer i, j, u, v
      double precision si, scale
c
c     Orthogonalize columns of a matrix distributed so that
c     each process has vecs(ulo:uhi,1:ni) ... uses global sums only.
c
      do i = 1, ni
 10      continue
c
c     If have a metric then first form overlap*vec(i)
c
       call dfill(i, 0.0d0, s, 1)
         do j = 1, i
            do u = ulo, uhi     ! May be a null range
               s(j) = s(j) + vecs(u,j)*vecs(u,i)
            enddo
         enddo
         call ga_dgop(11, s, i, '+')
         do j = 1, i-1
            do u = ulo, uhi     ! May be a null range
               vecs(u,i) = vecs(u,i) - vecs(u,j)*s(j)
            enddo
         enddo
c
         if (ometric) then
            call dfill(nu, 0.0d0, w, 1)
            do v = 1, nu
               do u = ulo, uhi
                  w(v) = w(v) + o(v,u)*vecs(u,i)
               enddo
            enddo
            call dgop(12, w, nu, '+') ! Actually only need ulo:hi
         else
            do u = ulo, uhi
               w(u) = vecs(u,i)
            enddo
         endif
c
         si = 0.0d0
         do u = ulo, uhi
            si = si + vecs(u,i)*w(u)
         enddo
         call ga_dgop(13, si, 1, '+')
         scale = 1.0d0/sqrt(si)
         do u = ulo, uhi
            vecs(u,i) = vecs(u,i) * scale
         enddo
c
c     If the vector norm changed a lot then repeat
c
         scale = abs(si/s(i))
         if (scale.lt.0.8 .or. scale.gt.1.2) then
            write(6,*) ' Doing again ', i
            goto 10
         endif
c
      enddo
c
      end



