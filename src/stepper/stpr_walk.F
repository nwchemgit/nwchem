      logical function stpr_walk(rtdb)
c $Id: stpr_walk.F,v 1.4 1994-08-06 09:18:12 d3e129 Exp $
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "context.fh"
c
      integer rtdb    ! [input] run-time-data-base handle
c
c parameters read from rtdb
c scf:energy (later) rohf and uhf
c
      integer geom              ! geometry handle
      integer needed            ! core needed
      integer nat               ! number of atoms 
      character*255 name        ! name buffer for rtdb context
      character*255 title       ! title buffer
      double precision energy   ! energy (scf, rohf, uhf)
      logical converged         ! energy convergence flag
c
      integer h_grad, k_grad    ! MA handle/index of gradient
      integer h_coord, k_coord  ! MA handle/index of coords
      integer h_core, k_core    ! MA handle/index of scratch array for stepper
      integer h_chg, k_chg      ! MA handle/index of charge array
      character*16 stpr_tags(max_cent) ! local tags array
      logical status
      double precision convgg, convge
c
      stpr_walk = .false.
c.................................... set broadcast off
      status = rtdb_parallel(.false.) 
c
*rak:      status = rtdb_print(rtdb,.true.)
      if (ga_nodeid().eq.0)  then
        if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $      title = ' '
        call util_print_centered(6, 'NWCHEM STEPPER Module', 40, .true.)
        write(6,*)
        write(6,*)
        if (title .ne. ' ') then
          call util_print_centered(6, title, 40, .false.)
          write(6,*)
          write(6,*)
        endif
c      
c     Push context down to SCF later look for wavefunction!!
c
        if (.not. context_push('scf'))
     &      call errquit('stpr_walk: context_push <scf> failed', 911)
c
        if (.not. context_prefix('converged', name))
     &      call errquit
     &      ('stpr_walk: context_prefix failed', 911)
        if (.not. rtdb_get(rtdb, name, MT_LOG, 1, converged))
     &      call errquit
     &      ('stpr_walk: failed to read converged in rtdb', 911)
        if (.not. converged)
     &      write(6,*)' WARNING: calculation not converged '
c
        if (.not. context_prefix('energy', name))
     &      call errquit
     &      ('stpr_walk: context_prefix failed', 911)
        if (.not. rtdb_get(rtdb, name, MT_DBL, 1, energy))
     &      call errquit
     &      ('stpr_walk: failed to read energy in rtdb', 911)
c... remove scf context
        if (.not. context_pop('scf'))
     &      call errquit('stpr_walk: context_pop failed',0)
c
        if (.not. context_push('stepper'))
     &      call errquit('stpr_walk: context_push <scf> failed', 911)
c
        if (.not. geom_create(geom, 'geometry'))
     &      call errquit('stpr_walk: geom_create?', 911)
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &      call errquit('stpr_walk: no geometry ', 911)
c
*        if (.not. geom_print(geom))
*     &      call errquit('stpr_walk: geom_print ?',911)      
*        call util_flush(6)
        if (.not. geom_ncent(geom,nat))
     &      call errquit('stpr_walk: geom_ncent?',911)
c
        call stpr_cneed(needed,nat) ! determine core needed
        needed = (needed + 1023)*1024/1024 ! nearest 1024 block
c
c..... space for core array
        if (.not.MA_Push_Get(MT_DBL,needed,'stpr core',h_core,k_core))
     &      call errquit('stpr_walk: allocation for core failed?',911)
c..... space for coords
        if (.not.
     &      MA_Push_Get(MT_DBL,(3*nat),'stpr coords',h_coord,k_coord))
     &      call errquit
     &        ('stpr_walk: allocation for coords failed?',911)
c..... space for gradient
        if (.not.
     &      MA_Push_Get(MT_DBL,(3*nat),'stpr gradient',h_grad,k_grad))
     &      call errquit
     &        ('stpr_walk: allocation for gradient failed?',911)
c..... space for charge
        if (.not.
     &      MA_Push_Get(MT_DBL,nat,'stpr geom charge',h_chg,k_chg))
     &      call errquit
     &        ('stpr_walk: allocation for geom charge failed?',911)
c
c..... read gradients
        if (.not. rtdb_get(rtdb, 'scf:gradients', MT_DBL, (3*nat),
     &      dbl_mb(k_grad)))
     &      call errquit
     &        ('stpr_walk: reading gradients failed',911)
c
c..... get coordinates
        if (.not. geom_cart_get(geom,nat,
     &      stpr_tags,
     &      dbl_mb(k_coord),
     &      dbl_mb(k_chg)))
     &      call errquit
     &        ('stpr_walk: geom_cart_get failed ?',911)
c
c... take a step
c
        convgg = 1.0d-04
        convge = 1.0d-07
        call stpr_stepper(
     &      dbl_mb(k_core),
     &      needed,
     &      dbl_mb(k_grad),
     &      dbl_mb(k_coord),
     &      energy,nat,convge,convgg,stpr_walk,
     &      rtdb)
c
c..... set coordinates
        if (.not. geom_cart_set(geom,nat,
     &      stpr_tags,
     &      dbl_mb(k_coord),
     &      dbl_mb(k_chg)))
     &      call errquit
     &        ('stpr_walk: geom_cart_set failed ?',911)
c
c..... update rtdb with new coordinates        
c
*        if (.not. geom_print(geom))
*     &      call errquit('stpr_walk: geom_print ?',911)
        if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &      call errquit
     &      ('stpr_walk: geom_rtdb_store failed ?',911)
c
c... clean up
        if (.not.geom_destroy(geom))
     &      call errquit('stpr_walk: geom_destroy failed?',911)
        if (.not.MA_Pop_Stack(h_chg))
     &      call errquit('stpr_walk: ma_pop_stack h_chg failed',911)
        if (.not.MA_Pop_Stack(h_grad))
     &      call errquit('stpr_walk: ma_pop_stack h_grad failed',911)
        if (.not.MA_Pop_Stack(h_coord))
     &      call errquit('stpr_walk: ma_pop_stack h_coord failed',911)
        if (.not.MA_Pop_Stack(h_core))
     &      call errquit('stpr_walk: ma_pop_stack h_core failed',911)
c
c... remove stepper context
        if (.not. context_pop('stepper'))
     &      call errquit('stpr_walk: context_pop failed',0)
      endif
c
c.... restore parallel access to rtdb
        status = rtdb_parallel(.true.) 
c      
      call ga_sync()
c
      end
      
