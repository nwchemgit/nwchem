c
c     quantum molecular dynamics task handler
c
      function task_qmd(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c
      integer rtdb
c
      logical  task_qmd
c
      logical status
      double precision cpu,wall
c
      character*32 theory
      character*84 tag
      character*30 pname
c
      logical qmd_doit
      external qmd_doit
c
      pname = "task_qmd: "
c
c     start clock
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     invoke the qmd driver
      status = qmd_doit(rtdb)
c
c     stop clock
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (status) then
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit(pname//'failed storing cputime',0,RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit(pname//'failed storing walltime',0,RTDB_ERR)
      endif
c         
      task_qmd = status
c
      end
c
c     main qmd driver
      function qmd_doit(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c
      integer rtdb
c
      logical  qmd_doit
c
      integer geom
      integer l_c, k_c
      integer l_g, k_g
      integer nat
      integer istep,nstep
      double precision energy_qm
c
      logical status
c
      character*32 theory
      character*84 tag
      character*30 pname
c
      logical task_gradient_doit
      external task_gradient_doit
c
      pname = "qmd_doit: "
c
c     get number of steps
c

c     outer md loop
      do istep = 1,nstep

c       geometry data structure
        if (.not. geom_create(geom, 'geometry'))
     &   call errquit(pname//'geom_create?', 911, GEOM_ERR)
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &   call errquit(pname//'no geometry ', geom, RTDB_ERR)
c
c       get the theory level
        if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     &   call errquit(pname//'theory not specified',0, RTDB_ERR)
c
c       number of atoms
        if (.not. geom_ncent(geom,nat))
     $   call errquit(pname//'nat?',geom, GEOM_ERR)
c
c       coordinates of geometry
        if(.not.ma_push_get(mt_dbl,3*nat,'coords',l_c,k_c))
     $     call errquit(pname//'memory', 3*nat, MA_ERR)
        if (.not. geom_cart_coords_get(geom, dbl_mb(k_c)))
     $        call errquit(pname//'geom?',0, GEOM_ERR)
c
c       perform energy gradient calculation
        if(.not.ma_push_get(mt_dbl,3*nat,'grads',l_g,k_g))
     $     call errquit(pname//'memory',3*nat, MA_ERR)
        energy_qm = 0.d0
        status = task_gradient_doit(rtdb,theory,energy_qm,dbl_mb(k_g))
c
c       change gradient to forces


c
c       move the atoms...
c
c       apply constraints
c
c       temperature


c
c       store updated geometry
        if(.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &     call errquit(pname//'failed to store geometry',0, GEOM_ERR)


      end do  ! istep
c
c     clear local memory




      qmd_doit = status
c
      end
