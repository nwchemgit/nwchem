c
c $Id: qmmm_esp.F,v 1.1 2005-06-15 05:39:00 marat Exp $
c
c
      function task_qmmm_esp(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb
      logical task_qmmm_esp
c
      character*32 pname
      character*80 espfile
      integer lespfile
c
      integer i,k
      logical ignore
c    
      integer igeom
      integer i_c_qm,h_c_qm
      integer i_q_qm,h_q_qm
      integer i_t_qm,h_t_qm      
      integer nqm
c
      integer i_ifg
c
      logical esp
      external esp
c
      pname = "task_qmmm_esp"
c
      task_qmmm_esp = .true.
c
      lespfile = 33
c
      if(.not.esp(rtdb)) 
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(rtdb,igeom,"esp:geometry"))
     + call errquit('qmmm: Failed to load geom to rtdb',0, RTDB_ERR)

      if(.not. geom_ncent(igeom,nqm))
     +    call errquit('qmmm: geom_ncent',0,0)

      if(.not.ma_push_get(mt_dbl,3*nqm,'c_qm',h_c_qm,i_c_qm))
     + call errquit('qmmm: Failed to qmocate memory for c_qm',
     + 3*nqm, MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*nqm,'q_qm',h_q_qm,i_q_qm))
     + call errquit('qmmm: Failed to qmocate memory for q_qm',
     + 3*nqm, MA_ERR)
  
      if(.not.ma_push_get(mt_byte,16*nqm,'t_qm',h_t_qm,i_t_qm))
     + call errquit('qmmm: Failed to qmocate memory for q_qm',nqm,
     &       MA_ERR)

      if(.not.geom_cart_get(igeom,nqm,byte_mb(i_t_qm),
     + dbl_mb(i_c_qm),dbl_mb(i_q_qm)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

        if(.not.geom_print_pdb(igeom,6))
     +    call errquit('qmmm: Failed to print xyzq geometry',0,0)

      do i=1,nqm
       write(luout,'(4f20.10)') (dbl_mb(i_c_qm+(i-1)*3+k-1),k=1,3),
     >                        dbl_mb(i_q_qm+i-1)
      end do

      if (.not. rtdb_put(rtdb,'qmmm:esp',
     >                   mt_dbl,nqm,dbl_mb(i_q_qm)))
     &    call errquit(pname//'qmmm:esp',911,RTDB_ERR)

      i_ifg = qmmm_geom_get_i_ifg()

      if(ga_nodeid().eq.0) then
        call util_file_name('esp',.false.,.false.,espfile)
        open(unit=lespfile,file=espfile,
     >       form='formatted',status='unknown')
        do  i=1,nqm
          write(lespfile,*) int_mb(i_ifg+i-1),
     >                           dbl_mb(i_q_qm+i-1)
        end do
        close(unit=lespfile)
      end if


       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      if (.not.ma_pop_stack(h_t_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_q_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_c_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)

      end

      subroutine qmmm_esp_creat_grid(rtdb,n,c)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb
      integer n
      double precision c(3,n)

      character*32 pname
      character*3 theory
      character*255 vectors_in
      character*80 grdfil
c
      integer i,k
      logical ignore
      integer h_epot,i_epot
c    
      integer lfngrd
      integer igeom
c      integer i_qs
c      integer nqm,nlink
      double precision charge
c
      logical esp
      external esp
      logical property
      external property
c
      pname = "qmmm_bq_energy"
c
      lfngrd = 33
c
      if (.not. ma_push_get(mt_dbl,n,'qmmm epot',h_epot,i_epot))
     &    call errquit(pname//' ma failed',911,MA_ERR)
c
c     make grid coords available for property calcs 
c     -------------------------------------------------------- 
      if (.not. rtdb_put(rtdb,'prop:xyz',
     >                   mt_dbl,3*n,c))
     &    call errquit(pname//'prop:xyz failed',911,RTDB_ERR)
c
c     export movecs file to property to avoid doing energy again
c     ----------------------------------------------------------
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit(pname//'theory not specified',555,
     &       INPUT_ERR)

       if (theory.eq.'dft') then
          if (.not. rtdb_cget(rtdb,'dft:input vectors',1,vectors_in))
     &        call errquit(pname//'rtdb_cget failed', 100,
     &                     RTDB_ERR)
       elseif (theory.eq.'scf') then
          if (.not. rtdb_cget(rtdb,'scf:input vectors',1,vectors_in))
     &        call errquit(pname//'rtdb_cget failed', 100,
     &                     RTDB_ERR)
       endif

      if (.not. rtdb_cput(rtdb,'prop:vectors',1,vectors_in))
     &    call errquit(pname//'rtdb_cput failed', 100,
     &                 RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:esp',
     >                   mt_int,1,0))
     &    call errquit(pname//'prop: efield failed',911,RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:nocheck',
     >                   mt_log,1,.true.))
     &    call errquit(pname//'prop:nocheck failed',911,RTDB_ERR)
c
      if(.not.property(rtdb)) 
     &    call errquit(pname//'failed property call', 100,
     &                 RTDB_ERR)

c
c     get back electric field at bq positions
c     ---------------------------------------
      if (.not. rtdb_get(rtdb,'prop:epot_xyz',
     >                   mt_dbl,n,dbl_mb(i_epot)))
     &    call errquit(pname//'prop:efield failed',911,RTDB_ERR)
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))then
        charge = 0.0d0
      endif

      if(ga_nodeid().eq.0) then
        call util_file_name('grid',.false.,.false.,grdfil)
        open(unit=lfngrd,file=grdfil,form='formatted',status='unknown')
        write(lfngrd,'(i10,f20.10)') n,charge
        do  i=1,n
          write(lfngrd,'(4f20.10)') c(1,i),c(2,i),c(3,i),
     >                           dbl_mb(i_epot+i-1)
        end do 
        close(unit=lfngrd)
      end if
c 
c     clean up rtdb
c     ------------
      ignore = rtdb_delete(rtdb,"prop:epot_xyz")
      ignore = rtdb_delete(rtdb,"prop:xyz")
      ignore = rtdb_delete(rtdb,"prop:vectors")

      if (.not.ma_pop_stack(h_epot)) call errquit
     &   (pname//'ma_pop_stack of h_epot failed',911,MA_ERR)


      end


