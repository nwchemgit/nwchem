       subroutine et_2e(g_pa,g_pb,omega2)
c
c      calculates the two-electron contribution to Vba
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
c
       integer g_pa,g_pb,g_pt,i,j,k,l
       double precision int2e1,int2e2,omega2
       double precision pa1,pa2,pb1,pb2,pt1,pt2
       double precision scr(mem_2e4c)
c
c      create array pt = pa + pb
c 
      if(.not.ga_create(mt_dbl,nbf,nbf,'pt',0,0,g_pt))
     $        call errquit('task_et ga_create failed', g_pt)
      call ga_dadd (1.0d0,g_pa,1.0d0,g_pb,g_pt)
c
      omega2 = 0.0d0
c
       DO i = 1,nbf
       DO j = 1,nbf
       DO k = 1,nbf
       DO l = 1,nbf
         call ga_get(g_pt, i, i, j, j, pt1, 1)
         call ga_get(g_pt, k, k, l, l, pt2, 1)
         call ga_get(g_pa, i, i, j, j, pa1, 1)
         call ga_get(g_pa, k, k, l, l, pa2, 1)
         call ga_get(g_pb, i, i, j, j, pb1, 1)
         call ga_get(g_pb, k, k, l, l, pb2, 1)
         call int_2e4c(basis,i,j,basis,k,l,mem_2e4c,scr,1,int2e1)
         call int_2e4c(basis,i,l,basis,k,j,mem_2e4c,scr,1,int2e2)
        omega2 = omega2 + ( pt1*pt2 )*int2e1
        omega2 = omega2 + ( -pa1*pa2-pb1*pb2 )*int2e2
       ENDDO
       ENDDO
       ENDDO
       ENDDO
c
c  Pt array is no longer needed.
      if (.not.ga_destroy(g_pt))
     &  call errquit('could not destroy Pt array',g_pt)
c
c
       return
       end
