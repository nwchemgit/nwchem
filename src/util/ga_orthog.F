      subroutine ga_orthog(g_vecs, g_over, ometric)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_vecs, g_over
      logical ometric
c
      integer ulo, uhi, type, ni, nu, uchunk, ld
      integer l_tmp, l_s, l_over, l_w
c     MAPOINTER k_tmp, k_s, k_over, k_w ! sign of things to come...
      integer k_tmp, k_s, k_over, k_w
c     
c     Redistribute the input matrix ... block the leading
c     dimension, leave second dimension undistributed
c     ... each process has a(ulo:uhi,1:ni).  If
c
c     nu = no. of basis functions
c     ni = no. of vectors
c
c     overlap should be (nu,nu)
c
      
      call ga_inquire(g_vecs, type, nu, ni)
      uchunk = max(8, (nu-1)/ga_nnodes()+1)
      ulo  = ga_nodeid()*uchunk + 1
      uhi  = min(ulo + uchunk - 1, nu)
      ld   = uhi - ulo + 1
      if (ulo .gt. uhi) then
         ulo = 0
         uhi = -1
         ld  = 1
      end if
*      write(6,*) ' node ulo uhi ', ga_nodeid(), ulo, uhi
*      call util_flush(6)
*      call ga_sync
*      call ga_print(g_vecs)
*      if (ometric) call ga_print(g_over)
c
c     Allocate local workspace
c
      if (.not. ma_push_get(mt_dbl, ld*ni, 'ga_orthog: tmp', 
     $     l_tmp, k_tmp)) call errquit('ga_orthog: ma get failed', 
     $     ld*ni, MA_ERR)
      if (ometric) then
         if (.not. ma_push_get(mt_dbl, ld*nu, 'ga_orthog: over', 
     $        l_over, k_over)) call errquit('ga_orthog: ma get failed', 
     $        ld*ni, MA_ERR)
      end if
      if (.not. ma_push_get(mt_dbl, ni, 'ga_orthog: sji', 
     $     l_s, k_s)) call errquit('ga_orthog: ma get failed', ni,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nu, 'ga_orthog: w', 
     $     l_w, k_w)) call errquit('ga_orthog: ma get failed', nu,
     &       MA_ERR)
c
      call ga_sync()
      if (uhi .ge. ulo) then
         call ga_get(g_vecs, ulo, uhi, 1, ni, dbl_mb(k_tmp), ld)
         if (ometric) 
     $        call ga_get(g_over, 1, nu, ulo, uhi, dbl_mb(k_over), nu)
      end if
c
c     Do the deed
c
      call ga_orthog_vecs(dbl_mb(k_tmp), ld, ni, ulo, uhi, dbl_mb(k_s),
     $     dbl_mb(k_over), nu, dbl_mb(k_w), ometric)
c
c     Put results back
c
      if (uhi .ge. ulo)
     $     call ga_put(g_vecs, ulo, uhi, 1, ni, dbl_mb(k_tmp), ld)
c
c     Tidy up memory
c
      if (.not. ma_pop_stack(l_w)) call errquit
     $     ('ga_orthog: ma error popping w', 0, MA_ERR)
      if (.not. ma_pop_stack(l_s)) call errquit
     $     ('ga_orthog: ma error popping s', 0, MA_ERR)
      if (ometric) then
         if (.not. ma_pop_stack(l_over)) call errquit
     $        ('ga_orthog: ma error popping over', 0, MA_ERR)
      end if
      if (.not. ma_pop_stack(l_tmp)) call errquit
     $     ('ga_orthog: ma error popping tmp', 0, MA_ERR)
c
      call ga_sync
c
      end
      subroutine ga_orthog_vecs(vecs, ld, ni, ulo, uhi, s, o, nu, w, 
     $     ometric)
      implicit none
#include "errquit.fh"
#include "msgids.fh"
c     
      integer ld, ni, ulo, uhi, nu
      double precision 
     $     vecs(ulo:(ulo+ld-1),1:ni), ! [input/output] Local patch of vectors
     $     o(1:nu,ulo:(ulo+ld-1)), ! [input] metric matrix if (ometric)
     $     s(ni),               ! [scratch] workspace for making overlap
     $     w(nu)                ! [scratch] workspace for making overlap
      logical ometric
c     
      integer i, j, u, npass
      double precision si, scale
c     
c     Orthogonalize columns of a matrix distributed so that
c     each process has vecs(ulo:uhi,1:ni) ... uses global sums only.
c     
      do i = 1, ni
         npass = 0
   10    npass = npass + 1
c     
c     If have a metric then first form overlap*vec(i)
c     
         if (ometric) then
            call dfill(nu, 0.0d0, w, 1)
            if ((uhi-ulo+1) .gt. 0) 
     $           call dgemv('n', nu, (uhi-ulo+1), 1.0d0, o, nu,
     $           vecs(ulo,i), 1, 0.0d0, w, 1)

            call ga_dgop(msg_orthog_1, w, nu, '+') ! Actually only need ulo:hi local
         else
            do u = ulo, uhi
               w(u) = vecs(u,i)
            end do
         end if
c     
c     Now form overlap between vector i and vectors 1...i
c     
         call dfill(i, 0.0d0, s, 1)
         if ((uhi-ulo+1) .gt. 0) 
     $        call dgemv('t', (uhi-ulo+1), i, 1.0d0,
     $        vecs(ulo,1), ld, w(ulo), 1, 0.0d0, s(1), 1)

         call ga_dgop(msg_orthog_2, s, i, '+')

c
c     Apply the rotation
c     
         if ((uhi-ulo+1).gt.0 .and. i.gt.1) then
            call dgemv('n', (uhi-ulo+1), i-1, -1.0d0,
     $           vecs(ulo,1), ld, s(1), 1, 1.0d0, vecs(ulo,i), 1)
         end if
c     
c     Renormalize vector i
c     
         si = s(i)
         do j = 1, i-1
            si = si - s(j)*s(j)
         end do
c     
c     If the vector norm changed a lot then repeat
c     
         if (i .gt. 1) then
            scale = si/s(i)
            if (scale .lt. 0.9d0) then
               if (npass .lt. 3) then
                  goto 10
               else
                  call errquit('ga_orthog: failed to orthog vector', i,
     &       GA_ERR)
               end if
            end if
         end if
c
         if (si .eq. 0.0d0) call errquit('ga_orthog: hard zero', i,
     &       0)
         if (si .lt. 0.0d0) then
            write(6,*) 'Negative value orthog: ',si
            call errquit('ga_orthog: negative', i, 0)
         endif
         scale = 1.0d0/sqrt(si)
         do u = ulo, uhi
            vecs(u,i) = vecs(u,i) * scale
         end do
c     
      end do
c
      end
