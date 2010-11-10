c---------------- get_chi_centers_ga() ------------- START
      subroutine get_chi_centers_ga(g_chi_cntr,basis,
     &                              nbf,geom,mcenters)
      implicit none

#include "rtdb.fh" 
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "msgids.fh"
      integer basis,geom,nbf
      integer g_chi_cntr(3) ! OUTPUT
      double precision cnt(3),valZ
      integer ictr,ic1,ic2,icset
      integer l,nprim,ncontr,isphere,nshbf
      integer mcenters,i,n1
      integer iniz(mcenters),ifin(mcenters)
      character*16 at_tag  
      integer l_buf,k_buf
      integer lo1(3),hi1(3),ld(2)
      logical status     
c ----- allocate array to store centers ---- START  
       if(.not.MA_push_get(MT_DBL,3*nbf,'get_chi_centers_ga:buf',
     &                    l_buf,k_buf))
     $     call errquit('get_chi_centers_ga: ma failed',
     &                  3*nbf, MA_ERR) 
c ----- allocate array to store centers ---- END
        call ifill(mcenters,0,iniz,1)
        call ifill(mcenters,0,ifin,1)
        n1=0
        do ictr=1,mcenters
         if (.not.bas_ce2cnr(basis,ictr,ic1,ic2))
     &       call errquit('Exiting in get_chi_centers_ga.',
     &                    11, BASIS_ERR)
         do icset = ic1,ic2       
c ----- get info about current contraction set      
          if (.not. bas_continfo(basis,icset,l,nprim,
     &         ncontr,isphere))
     &         call errquit('Exiting in get_chi_centers_ga.',
     &                       5, BASIS_ERR)
          nshbf=ncontr*(((l+1)*(l+2))/2)
          if(isphere.eq.1) then
            nshbf=ncontr*(2*l+1)
          endif
          if (iniz(ictr).eq.0) iniz(ictr)=n1+1
          n1=n1+nshbf
         enddo ! end loop icset
         ifin(ictr)= n1
         status=geom_cent_get(geom,ictr,at_tag,
     &                        cnt,valZ)
         do i=iniz(ictr),ifin(ictr)
           dbl_mb(k_buf      +i-1)=cnt(1)
           dbl_mb(k_buf+nbf  +i-1)=cnt(2)
           dbl_mb(k_buf+2*nbf+i-1)=cnt(3)
         enddo ! end loop i 
        enddo ! end loop ictr
c ----- store in g_chi_cntr() --- START
c       dbl_mb() ---> g_chi_cntr()
        do i=1,3
         ld(1)=nbf
         lo1(1)=1
         hi1(1)=nbf
         lo1(2)=i
         hi1(3)=i
         call nga_put(g_chi_cntr(i),
     &                lo1,hi1,dbl_mb(k_buf+(i-1)*nbf),ld)
        enddo ! end-loop-i
c ----- store in g_chi_cntr() --- END
c --- Free memory 
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('get_chi_centers_ga: pop failed', 0, GA_ERR)
      return
      end
c---------------- get_chi_centers_ga() ------------- END
      subroutine get_3rdterm_R(g_N,     ! to be scaled
     &                         g_R,     ! scaling
     &                         ind_a,   ! from kab=123,231,312
     &                         ind_b,   ! from kab=123,231,312
     &                         g_tmp2,  ! scratch
     &                         g_N_scld)! output
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_N,g_R(3),g_M
      integer g_tmp2,g_N_scld
      integer ind_a,ind_b
      call ga_copy(g_N,g_tmp2)
      call ga_copy(g_N,g_N_scld)
      call ga_scale_cols(g_tmp2,g_R(ind_b))   ! R_{nu,b} g_N
      call ga_scale_rows(g_tmp2,g_R(ind_a))   ! R_{mu,a} [R_{nu,b} g_N] -> g_tmp2
      call ga_scale_cols(g_N_scld,g_R(ind_a)) ! R_{nu,a} g_N
      call ga_scale_rows(g_N_scld,g_R(ind_b)) ! R_{mu,b} [R_{nu,a} g_N] -> g_N_scld
      call ga_add(1.0d0,g_tmp2,-1.0d0,g_N_scld,g_N_scld)
      return
      end

      subroutine get_scld_A(g_A,  ! ga-arr to scale - OUT
     &                      g_R,  ! scaling arr
     &                      g_tmp)! scratch arr 
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_A,g_R
      integer g_tmp
      integer nbf
c     Purpose: Compute R_{nu} U_{munu} - R_{mu} U_{munu}
c              g_R ->  R_{mu}
c              g_A -> U_{munu}
      call ga_copy(g_A,g_tmp)
      call ga_scale_cols(g_A  ,g_R)
      call ga_scale_rows(g_tmp,g_R)
      call ga_add(1.0d0,g_A,-1.0d0,g_tmp,g_A)
      return
      end
