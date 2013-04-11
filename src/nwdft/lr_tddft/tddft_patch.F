      double precision function tga_ddot(g_a,g_b,k,l)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer g_a,g_b
      integer a,k,l
c
      integer me
      integer anlo,anhi,amlo,amhi,k_a,adim
      integer bnlo,bnhi,bmlo,bmhi,k_b,bdim
      integer m,oa,ob,alen, k_t1, l_t1, k_t2, l_t2
      double precision x
c
      me = ga_nodeid()
      x=0.0d0
      if (g_a.eq.g_b) then
         call ga_distribution(g_a,me,anlo,anhi,amlo,amhi)
         adim=anhi-anlo+1
         if (.not.ma_push_get(mt_dbl,adim,'Temp1',l_t1,k_t1)) 
     $      call errquit('tga_ddot:fail allocate t1',0, MA_ERR)
         if (k.eq.l) then
            call ga_get(g_a,anlo,anhi,k,k,dbl_mb(k_t1),adim)
            call tga_ddot_oper(x,dbl_mb(k_t1),dbl_mb(k_t1),adim)
         else
            if (.not.ma_push_get(mt_dbl,adim,'Temp2',l_t2,k_t2))   
     $         call errquit('tga_ddot:fail allocate t2',0, MA_ERR)
            call ga_get(g_a,anlo,anhi,l,l,dbl_mb(k_t2),adim)
            call tga_ddot_oper(x,dbl_mb(k_t1),dbl_mb(k_t2),adim)
            if (.not.ma_pop_stack(l_t2)) call errquit
     $         ('tga_ddot:failed to deallocate t2',0, MA_ERR)
         endif
         if (.not.ma_pop_stack(l_t1)) call errquit
     $      ('tga_ddot:failed to deallocate t1',0, MA_ERR)
      else
         call ga_distribution(g_a,me,anlo,anhi,amlo,amhi)
         call ga_distribution(g_b,me,bnlo,bnhi,bmlo,bmhi)
         adim=anhi-anlo+1
         bdim=bnhi-bnlo+1
         if (.not.ma_push_get(mt_dbl,adim,'Temp',l_t1,k_t1))
     $      call errquit('tga_ddot:fail allocate t1',0, MA_ERR)
         if (.not.ma_push_get(mt_dbl,bdim,'Temp',l_t2,k_t2))
     $      call errquit('tga_ddot:fail allocate t2',0, MA_ERR)
         call ga_get(g_a,anlo,anhi,k,k,dbl_mb(k_t1),adim)
         call ga_get(g_b,bnlo,bnhi,l,l,dbl_mb(k_t2),bdim)
         call tga_ddot_oper(x,dbl_mb(k_t1),dbl_mb(k_t2),adim)
         if (.not.ma_pop_stack(l_t2)) call errquit
     $      ('tga_ddot:failed to deallocate t2',0, MA_ERR)
         if (.not.ma_pop_stack(l_t1)) call errquit
     $      ('tga_ddot:failed to deallocate t1',0, MA_ERR)
      endif
      call ga_dgop(2123,x,1,'+')
      tga_ddot=x
      return
      end
c
      subroutine tga_ddot_oper(x,a,b,lda)
c
      implicit none
c
      double precision a(lda), b(lda), x
      integer lda, k, l, i
c
      do i=1,lda
         x=x+a(i)*b(i)
      enddo
      return
      end
c
      subroutine tga_add(x,g_a,g_b,k,l,m)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer g_a,g_b
      integer k,l,m
      double precision x
c
      integer me
      integer anlo,anhi,amlo,amhi,k_a,adim
      integer bnlo,bnhi,bmlo,bmhi,k_b,bdim
      integer n,chi,oa,ob,alen, k_t1, l_t1, k_t2, l_t2
c
      me = ga_nodeid()
      call ga_distribution(g_a,me,anlo,anhi,amlo,amhi)
      call ga_distribution(g_b,me,bnlo,bnhi,bmlo,bmhi)
      adim=anhi-anlo+1
      bdim=bnhi-bnlo+1
      if (.not.ma_push_get(mt_dbl,adim*m,'Temp1',l_t1,k_t1))
     $   call errquit('tga_ddot:fail allocate t1',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,bdim*m,'Temp2',l_t2,k_t2))
     $   call errquit('tga_ddot:fail allocate t1',0, MA_ERR)
      call ga_get(g_a,anlo,anhi,k,k+m-1,dbl_mb(k_t1),adim)
      call ga_get(g_b,bnlo,bnhi,l,l+m-1,dbl_mb(k_t2),bdim)
      call tga_add_oper(x,dbl_mb(k_t1),dbl_mb(k_t2),adim,m)
      call ga_put(g_b,bnlo,bnhi,l,l+m-1,dbl_mb(k_t2),bdim)
      if (.not.ma_pop_stack(l_t2)) call errquit
     $   ('tga_ddot:failed to deallocate t2',0, MA_ERR)
      if (.not.ma_pop_stack(l_t1)) call errquit
     $   ('tga_ddot:failed to deallocate t1',0, MA_ERR)
      return
      end
c
      subroutine tga_add_oper(x,a,b,lda,m)
c
      implicit none
c
      double precision a(lda,m), b(lda,m), x
      integer lda, m
      integer i, j
c
      do i=1,m
         do j=1,lda
            b(j,i)=b(j,i)+x*a(j,i)
         enddo
      enddo
      return
      end

