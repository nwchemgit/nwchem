*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_accy(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "global.fh"
      logical int_normalize
      external int_normalize
      integer rtdb
      integer geom, basis
      integer szbuf, maxg, maxscr, lmax, nshell
      integer size_grad, size_scr, size_labs
      integer h_grad_1, h_grad_2, h_labs, h_scr
      integer k_grad_1, k_grad_2, k_labs, k_scr
      logical status
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_accy: geom_create failed?',911,
     &       GEOM_ERR)
c
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_accy: geom_rtdb_load -ref failed',911, RTDB_ERR)
      
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &      ('bas_rtdb_load failed',911, RTDB_ERR)
c
      if (ga_nodeid().eq.0) then
        if (.not.geom_print(geom)) stop ' print error'
        if (.not.bas_print(basis)) stop ' print error'
        if (.not.gbs_map_print(basis)) stop ' gbs_map_print 2?'
      endif
      if (.not.int_normalize(rtdb,basis)) stop ' norm error 1'
      if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont fe'
      if (.not.bas_high_angular(basis,Lmax)) stop 'bas_ha fe'
c
      szbuf = (Lmax+1)*(Lmax+2)/2
      szbuf = szbuf**4
c
      call intd_init(rtdb,1,basis)
      call intb_mem_2e4c(maxg,maxscr)
      maxscr = maxscr + szbuf
      maxscr = maxscr + maxscr/5 + 1
      maxscr = max(51000,maxscr)
c
      size_grad = 12*szbuf ! 2 grad buffers
      size_scr  = maxscr
      size_labs = 4*szbuf ! for labels
c
      if (.not.ma_push_get(mt_dbl,size_grad,' grad buf 1 ',
     &    h_grad_1, k_grad_1)) stop ' push failed grad 1'
      if (.not.ma_push_get(mt_dbl,size_grad,' grad buf 2 ',
     &    h_grad_2, k_grad_2)) stop ' push failed grad 2'
      if (.not.ma_push_get(mt_dbl,size_scr,' scratch ',
     &    h_scr, k_scr)) stop ' push failed scratch '
      if (.not.ma_push_get(mt_int,size_labs,' labels ',
     &    h_labs, k_labs)) stop ' push failed labels'
      call task_accy_a(rtdb,geom,basis,
     &    size_grad,size_scr,size_labs,nshell,
     &    dbl_mb(k_grad_1),dbl_mb(k_grad_2),
     &    dbl_mb(k_scr),int_mb(k_labs))
c
      
c
      status = .true.
      status = status.and.ma_pop_stack(h_labs)
      status = status.and.ma_pop_stack(h_scr)
      status = status.and.ma_pop_stack(h_grad_1)
      status = status.and.ma_pop_stack(h_grad_2)
      task_accy = status
c
      end
C> @}
      subroutine task_accy_a(rtdb,geom,basis,
     &    lgrad,lscr,llab,nshell,
     &    grad1, grad2, scr, labs)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
c::functions
      logical intb_init4c
      logical intb_2e4c
      external intb_init4c
      external intb_2e4c
      integer rtdb,geom,basis
      integer lgrad, lscr, llab
      integer nshell
      double precision grad1(12,*)
      double precision grad2(12,*)
      double precision scr(lscr)
      integer labs(4,*)
c
      double precision q4
      double precision dummy
      integer me, nproc
      integer task_count
      integer i, j, k, l
      integer ish, jsh, ksh, lsh
      integer tish, tjsh, tksh, tlsh
      integer ilo, ihi, nbfi
      integer jlo, jhi, nbfj
      integer klo, khi, nbfk 
      integer llo, lhi, nbfl
      integer nint
c
      q4 = 1.0d00
      dummy = 1.0d00
      me    = ga_nodeid()
      nproc = ga_nnodes()
      task_count = (me-1)
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) stop 'cn2bfr i'        
        nbfi = ihi - ilo + 1
        do jsh = 1,nshell
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) stop 'cn2bfr j'        
          nbfj = jhi - jlo + 1
          do ksh = 1,nshell
            if (.not.bas_cn2bfr(basis,ksh,klo,khi)) stop 'cn2bfr k'        
            nbfk = khi - klo + 1
            do lsh = 1,nshell
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi)) stop 'cn2bfr l'        
              nbfl = lhi - llo + 1
              nint = nbfi*nbfj*nbfk*nbfl
              task_count = task_count + 1
              if (mod(task_count,nproc).eq.0) then
                tish = ish
                tjsh = jsh
                tksh = ksh
                tlsh = lsh
                call ifill((4*nint),0,labs,1)
                call dfill(nint,0.0d00,grad1,1)
                if (.not.intb_init4c(
     &              basis,tish,tjsh,
     &              basis,tksh,tlsh,1,
     &              q4,.false.,lscr,scr,nint,dummy)) stop 'error'
00100           continue
              endif
            enddo
          enddo
        enddo
      enddo
c
      end
