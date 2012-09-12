c
c $Id$
c
C> \ingroup task
C> @{
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
#include "inp.fh"
c
      integer rtdb
      logical task_qmmm_esp
c
      character*32 pname
      character*255 espfile
      character*255 prefix
      character*255 movfil
      character*80 tag
c      character*80 movecs_tag
c      character*30 theory
      integer lespfile
      logical lrism
c
      integer i,k
      logical ignore
      logical omovecs
c    
      integer igeom
      integer i_c_qm,h_c_qm
      integer i_q_qm,h_q_qm
      integer i_t_qm,h_t_qm      
      integer nqm
c
      integer i_ifg
      integer itmp
      logical status
      logical readesp
      logical ofile
      logical bqgrid
c
      logical  geom_rtdb_ncent
      external geom_rtdb_ncent
c
      logical  qmmm_esp_load
      external qmmm_esp_load
c
      logical esp
      external esp
c
      pname = "task_qmmm_esp"
c
      if(qmmm_print_debug()) write(luout,*) "in "//pname
c
      task_qmmm_esp = .true.
c
c      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
c     + call errquit('task: no task input for theory?',0, RTDB_ERR)

c
c      call qmmm_esp_input_file(rtdb)
c
      if(.not.rtdb_get(rtdb,'qmmm:readesp',mt_log,1,readesp))
     +   readesp = .false.
c
      if(.not.rtdb_get(rtdb,'qmmm:bqgrid',mt_log,1,bqgrid))
     +   bqgrid = .false.
c      
c     read esp charges from external file if any
c     ------------------------------------------
      if(readesp) then
        if(.not.qmmm_esp_load(rtdb))
     >    call errquit('cannot load espfile',0,0)
        return
      end if
c
      lespfile = 33
c
      i_ifg = qmmm_geom_get_i_ifg()
c
      if(bqgrid) then
        call qmmm_esp_creat_bq_grid(rtdb)
      else
        if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)
      end if
c
c      if(theory.eq."dft") then
c        movecs_tag = 'dft:input vectors'
c      else if(theory.eq."scf") then
c        movecs_tag = "scf:input vectors"
c      else
c        movecs_tag = " "
c      end if
c      if(movecs_tag.ne." ") then
c        if(.not.rtdb_cget(rtdb,movecs_tag,1,movfil))
c     + call errquit(pname//'no movecs?',0, RTDB_ERR)
c        if(.not.rtdb_cput(rtdb,"esp:input vectors",1,movfil))
c     + call errquit(pname//'put movecs?',0, RTDB_ERR)
c      end if
cc
c
c      call util_file_name('movecs',.false.,.false.,movfil)
c      if(ga_nodeid().eq.0) then
c           inquire(file=movfil,exist=omovecs)
c           if(omovecs) then
c            write(*,*) "qmmm esp using movecs file ",
c     >                movfil(1:inp_strlen(movfil))
c           else
c            call errquit(pname//"default movecs file"
c     >                   //movfil(1:inp_strlen(movfil))
c     >                   //" not found",0,0)
c           end if
c       end if

c
      if(.not.esp(rtdb)) 
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

c
c     pick the best available charge geometry
c     ---------------------------------------
      tag = "cresp:geometry"
      status = geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      if(.not.status) then
        tag = "resp:geometry"
        status =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if
      if(.not.status) then
        tag = "cesp:geometry"
        status =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if
      if(.not.status) then
        tag = "esp:geometry"
        status =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if

      if(.not.status)
     + call errquit(pname,0, RTDB_ERR)

      status = geom_rtdb_load(rtdb,igeom,tag(1:inp_strlen(tag)))

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

c      if(ga_nodeid().eq.0) then
c        call util_print_centered(luout,"QM/MM esp charges",32,.true.)
c        do i=1,nqm
c         write(luout,'(4f20.10)') (dbl_mb(i_c_qm+(i-1)*3+k-1),k=1,3),
c     >                          dbl_mb(i_q_qm+i-1)
c        end do
c      end if

      ignore = rtdb_delete(rtdb,'qmmm:esp')
      if (.not. rtdb_put(rtdb,'qmmm:esp',
     >                   mt_dbl,nqm,dbl_mb(i_q_qm)))
     &    call errquit(pname//'qmmm:esp',911,RTDB_ERR)

      if (.not.rtdb_cget(rtdb,'qmmm:esp:save:file',1,espfile)) then
        call mm_restart_filename(prefix)
        i = index(prefix,".rst")-1
        espfile = prefix(1:i)//".esp"
      end if
      call util_file_name_resolve(espfile,.false.)
      if(ga_nodeid().eq.0) then
        open(unit=lespfile,file=espfile,
     >       form='formatted',status='unknown')
        do  i=1,nqm
          write(lespfile,*) int_mb(i_ifg+i-1),
     >                           dbl_mb(i_q_qm+i-1)
        end do
        close(unit=lespfile)
      end if

c
c     temporary rism hookup
c     to produce charge file
c     ---------------------
      if (.not.rtdb_get(rtdb,'qmmm:rism',mt_log,1,lrism)) lrism=.false.
      if(lrism) then
        call mm_restart_filename(prefix)
        i = index(prefix,".rst")-1
        espfile = prefix(1:i)//".xyzq"
      end if
      call util_file_name_resolve(espfile,.false.)
      if(ga_nodeid().eq.0) then
        open(unit=lespfile,file=espfile,
     >       form='formatted',status='unknown')
        do  i=1,nqm
         write(luout,*) (dbl_mb(i_c_qm+(i-1)*3+k-1),k=1,3),
     >                          dbl_mb(i_q_qm+i-1)
        end do
        close(unit=lespfile)
      end if
c     end of rism hookup

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      if (.not.ma_pop_stack(h_t_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_q_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_c_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)

      if(qmmm_print_debug()) write(luout,*) "out "//pname

      return

 134  call errquit(pname//'error reading espfile',0, 0)

      end
C> @}

      function qmmm_esp_load(rtdb)
c     reads esp charges from external file located
c     in perm directory
c     generates fatal error if file not found
c     or if there are read errors
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
c
      integer rtdb
      logical qmmm_esp_load
c
      character*32 pname
      character*255 espfile
      character*255 prefix
      integer lespfile
      integer i
      logical ignore
      integer i_q_qm,h_q_qm
      integer nfg
      integer i_ifg
      integer itmp
      logical status
      logical ofile
      logical oprint
c
      pname = "qmmm_esp_load"
c
      if(qmmm_print_debug()) write(luout,*) "in "//pname
c
      oprint = ga_nodeid().eq.0
      qmmm_esp_load = .false.
c
      if(.not.qmmm_get_io_unit(lespfile))
     >    call errquit(pname//"no free file unit",0,0)

c     
c     construct espfilename
c     ----------------------
      if(.not.rtdb_cget(rtdb,"qmmm:espfilename",1,espfile)) then
        call mm_restart_filename(prefix)
        i = index(prefix,".rst")-1
        espfile = prefix(1:i)//".esp"
        call util_file_name_resolve(espfile,.false.)
      end if
      call util_file_name_resolve(espfile,.false.)
c     
c     aborting if no espfile
c     ----------------------
      if(ga_nodeid().eq.0) then
        inquire(file=espfile,exist=ofile)
        if(.not.ofile) 
     >   call errquit(pname//'espfile not found'//espfile,0,0)
      end if
c
c     number of esp charges = full geometry
c     -------------------------------------
      nfg = qmmm_geom_get_nfg()
      i_ifg = qmmm_geom_get_i_ifg()
c
      if(.not.ma_push_get(mt_dbl,3*nfg,'q_qm',h_q_qm,i_q_qm))
     + call errquit(pname//'Failed to qmocate memory for q_qm',
     + 3*nfg, MA_ERR)
 
c
      if(oprint) 
     >  write(luout,*) "loading esp charges from external file"
      if(ofile.and.ga_nodeid().eq.0) then
        open(unit=lespfile,file=espfile,
     >       form='formatted',status='unknown')
        do  i=1,nfg
          read(lespfile,*,err=134,end=134) itmp,
     >                           dbl_mb(i_q_qm+i-1)
          if(itmp.ne.int_mb(i_ifg+i-1))
     >      call errquit(pname//'global index mismatch',0,0)
          if(oprint)
     >      write(luout,*) itmp,dbl_mb(i_q_qm+i-1)
        end do
        close(unit=lespfile)
      end if
c
c     push esp charges to rtdb
c     ------------------------
      ignore = rtdb_parallel(.false.)
      if(ga_nodeid().eq.0) then
        ignore = rtdb_delete(rtdb,'qmmm:esp')
        if (.not. rtdb_put(rtdb,'qmmm:esp',
     >                     mt_dbl,nfg,dbl_mb(i_q_qm)))
     &      call errquit(pname//'qmmm:esp',911,RTDB_ERR)
      end if
      ignore = rtdb_parallel(.true.)
      call ga_sync()
c
      if (.not.ma_pop_stack(h_q_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)

      if(qmmm_print_debug()) write(luout,*) "out "//pname

      qmmm_esp_load = .true.

      return

 134  call errquit(pname//'error reading espfile',0, 0)
      

      end

      function qmmm_esp_load_mix(rtdb,lambda)
c     reads esp charges from external file located
c     in perm directory
c     generates fatal error if file not found
c     or if there are read errors
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
c
      integer rtdb
      double precision lambda
      logical qmmm_esp_load_mix
c
      character*32 pname
      character*255 espfile
      character*255 prefix
      integer lespfile
      integer i
      logical ignore
      integer i_q_qm,h_q_qm
      integer i_q_qm0,h_q_qm0
      integer nfg
      integer i_ifg
      integer itmp
      logical status
      logical ofile
      logical oprint
c
      pname = "qmmm_esp_load"
c
      if(qmmm_print_debug()) write(luout,*) "in "//pname
c
      oprint = ga_nodeid().eq.0
      qmmm_esp_load_mix = .false.
c
      if(.not.qmmm_get_io_unit(lespfile))
     >    call errquit(pname//"no free file unit",0,0)

c     
c     construct espfilename
c     ----------------------
      if(.not.rtdb_cget(rtdb,"qmmm:espfilename",1,espfile)) then
        call mm_restart_filename(prefix)
        i = index(prefix,".rst")-1
        espfile = prefix(1:i)//".esp"
        call util_file_name_resolve(espfile,.false.)
      end if
      call util_file_name_resolve(espfile,.false.)
c     
c     aborting if no espfile
c     ----------------------
      if(ga_nodeid().eq.0) then
        inquire(file=espfile,exist=ofile)
        if(.not.ofile) 
     >   call errquit(pname//'espfile not found'//espfile,0,0)
      end if
c
c     number of esp charges = full geometry
c     -------------------------------------
      nfg = qmmm_geom_get_nfg()
      i_ifg = qmmm_geom_get_i_ifg()
c
      if(.not.ma_push_get(mt_dbl,nfg,'q_qm',h_q_qm,i_q_qm))
     + call errquit(pname//'Failed to qmocate memory for q_qm',
     + nfg, MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,nfg,'q_qm0',h_q_qm0,i_q_qm0))
     + call errquit(pname//'Failed to qmocate memory for q_qm0',
     + nfg, MA_ERR)
 
c
      if(oprint) 
     >  write(luout,*) "loading esp charges from external file"
      if(ofile.and.ga_nodeid().eq.0) then
        open(unit=lespfile,file=espfile,
     >       form='formatted',status='unknown')
        do  i=1,nfg
          read(lespfile,*,err=134,end=134) itmp,
     >                           dbl_mb(i_q_qm+i-1)
          if(itmp.ne.int_mb(i_ifg+i-1))
     >      call errquit(pname//'global index mismatch',0,0)
          if(oprint)
     >      write(luout,*) itmp,dbl_mb(i_q_qm+i-1)
        end do
        close(unit=lespfile)
      end if
c
c     push esp charges to rtdb
c     ------------------------
      ignore = rtdb_parallel(.false.)
      if(ga_nodeid().eq.0) then
        if(oprint) 
     >    write(luout,*) "scaling esp charges"
        if (.not. rtdb_get(rtdb,'qmmm:esp',
     >                     mt_dbl,nfg,dbl_mb(i_q_qm0)))
     &      call errquit(pname//'qmmm:esp',911,RTDB_ERR)
        do i=1,nfg
          dbl_mb(i_q_qm+i-1)=dbl_mb(i_q_qm0+i-1)+lambda*
     >                       (dbl_mb(i_q_qm+i-1)-dbl_mb(i_q_qm0+i-1))
          if(oprint)
     >      write(luout,*) dbl_mb(i_q_qm+i-1)
        end do
        if (.not. rtdb_put(rtdb,'qmmm:esp',
     >                     mt_dbl,nfg,dbl_mb(i_q_qm)))
     &      call errquit(pname//'qmmm:esp',911,RTDB_ERR)
        call util_flush(luout)
      end if
      ignore = rtdb_parallel(.true.)
      call ga_sync()
c
      if (.not.ma_pop_stack(h_q_qm0)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_q_qm)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)

      if(qmmm_print_debug()) write(luout,*) "out "//pname

      qmmm_esp_load_mix = .true.

      return

 134  call errquit(pname//'error reading espfile',0, 0)
      

      end

      subroutine qmmm_esp_creat_bq_grid(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "bq.fh"
c
      integer rtdb
c
      character*30 pname
      integer n
      integer i_cbq, bq_handle

      pname = "qmmm_esp_creat_bq_grid"
      if(.not.geom_extbq_on())
     >   call errquit(pname//'no active bqs',0,0) 
      n = geom_extbq_ncenter()
      i_cbq = geom_extbq_coord()
c
        if(.not.bq_get_active(bq_handle)) 
     &    call errquit(pname//'bq_handle',911,RTDB_ERR)
        if(.not.bq_deactivate(bq_handle)) call errquit(
     &        pname//'unable to deactivate bqs',0,0)

c
      call qmmm_esp_creat_grid(rtdb,n,dbl_mb(i_cbq))

      if(.not.bq_activate(bq_handle)) call errquit(
     &        pname//'unable to activate bqs',0,0)

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
      character*255 grdfil
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


       if (rtdb_cget(rtdb,'esp:input vectors',1,vectors_in)) then
        if (.not. rtdb_cput(rtdb,'prop:vectors',1,vectors_in))
     &      call errquit(pname//'rtdb_cput failed', 100,
     &                   RTDB_ERR)
       end if

      if (.not. rtdb_put(rtdb,'prop:esp',
     >                   mt_int,1,0))
     &    call errquit(pname//'prop: efield failed',911,RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:nocheck',
     >                   mt_log,1,.true.))
     &    call errquit(pname//'prop:nocheck failed',911,RTDB_ERR)
c
c
      if(.not.property(rtdb)) 
     &    call errquit(pname//'failed property call', 100,
     &                 RTDB_ERR)

c
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

      subroutine qmmm_esp_reset(rtdb)
      implicit none
#include "rtdb.fh"
c
      integer rtdb
      logical ignore

      ignore = rtdb_delete(rtdb,'qmmm:esp')

      end

      subroutine qmmm_bq_esp_grad(rtdb,nbq,qbq,cbq,g)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "qmmm.fh"
c
      integer rtdb
      integer nbq
      double precision qbq(nbq)
      double precision cbq(3,nbq)
      double precision g(3,*)
c
      character*32 pname
c
      integer i,k
c    
      integer nqm
      integer nlink
      integer igeom
c
      double precision zqm,zbq
      double precision xp, yp, zp, xn, yn, zn
      double precision rr,rr3
c
      integer i_c_esp,h_c_esp
      integer i_q_esp,h_q_esp
      integer i_cbq
      integer i_qbq
      logical master 
      logical esp
      external esp
      logical  task_qmmm_esp
      external task_qmmm_esp
c
      pname = "qmmm_bq_esp_energy"
c
      master = qmmm_master()
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
c
      if(.not.geom_rtdb_load(rtdb,igeom,"geometry"))
     + call errquit('qmmm: Failed to load geom to rtdb',0, RTDB_ERR)
c
      if(.not. geom_ncent(igeom,nqm))
     +    call errquit('qmmm: geom_ncent',0,0)
c
      if(.not.ma_push_get(mt_dbl,3*nqm,'c_esp',h_c_esp,i_c_esp))
     + call errquit('qmmm: Failed to qmocate memory for c_qm',
     + 3*nqm, MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,3*nqm,'q_esp',h_q_esp,i_q_esp))
     + call errquit('qmmm: Failed to qmocate memory for q_qm',
     + 3*nqm, MA_ERR)
c  
      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c_esp)))
     + call errquit(pname,0, GEOM_ERR)
c
       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)
c
111   continue
      if (.not. rtdb_get(rtdb,'qmmm:esp',
     >    mt_dbl,nqm,dbl_mb(i_q_esp))) then
       if(master) write(luout,*) "redoing esp charges"
       if(.not.task_qmmm_esp(rtdb)) 
     >   call errquit('failed qmmm esp',0, 0)
       goto 111
      else
       if(master) write(luout,*) "using old esp charges"
      end if
 
      do i=1,nqm
        xp = dbl_mb(i_c_esp  +3*(i-1))
        yp = dbl_mb(i_c_esp+1+3*(i-1))
        zp = dbl_mb(i_c_esp+2+3*(i-1))
        zqm = dbl_mb(i_q_esp+i-1)
        do k = 1,nbq
          zbq = qbq(k)
          xn  = cbq(1,k) - xp
          yn  = cbq(2,k) - yp
          zn  = cbq(3,k) - zp
          rr =  sqrt(xn*xn + yn*yn + zn*zn)
          rr3=rr*rr*rr
          g(1,k)= g(1,k) - zqm*zbq*xn/rr3
          g(2,k)= g(2,k) - zqm*zbq*yn/rr3
          g(3,k)= g(3,k) - zqm*zbq*zn/rr3
        end do
      end do

      if (.not.ma_pop_stack(h_q_esp)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_c_esp)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)


      end

      subroutine qmmm_bq_esp_energy(rtdb,nbq,qbq,cbq,e)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "qmmm.fh"
c
      integer rtdb
      integer nbq
      double precision qbq(nbq)
      double precision cbq(3,nbq)
      double precision e
c
      character*32 pname
c
      integer i,k
c    
      integer nqm
      integer nlink
      integer igeom
c
      double precision zqm,zbq
      double precision xp, yp, zp, xn, yn, zn
      double precision rr,rr3
c
      integer i_c_esp,h_c_esp
      integer i_q_esp,h_q_esp
      integer i_cbq
      integer i_qbq
      logical master
      logical esp
      external esp
      logical  task_qmmm_esp
      external task_qmmm_esp
c
      pname = "qmmm_bq_esp_energy"
c
      master = qmmm_master()
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
c
      if(.not.geom_rtdb_load(rtdb,igeom,"geometry"))
     + call errquit('qmmm: Failed to load geom to rtdb',0, RTDB_ERR)
c
      if(.not. geom_ncent(igeom,nqm))
     +    call errquit('qmmm: geom_ncent',0,0)
c
      if(.not.ma_push_get(mt_dbl,3*nqm,'c_esp',h_c_esp,i_c_esp))
     + call errquit('qmmm: Failed to qmocate memory for c_qm',
     + 3*nqm, MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,3*nqm,'q_esp',h_q_esp,i_q_esp))
     + call errquit('qmmm: Failed to qmocate memory for q_qm',
     + 3*nqm, MA_ERR)
c  
      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c_esp)))
     + call errquit(pname,0, GEOM_ERR)
c
       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)
c
111   continue
      if (.not. rtdb_get(rtdb,'qmmm:esp',
     >    mt_dbl,nqm,dbl_mb(i_q_esp))) then
       if(master) write(luout,*) "redoing esp charges"
       if(.not.task_qmmm_esp(rtdb)) 
     >   call errquit('failed qmmm esp',0, 0)
       goto 111
      else
       if(master) write(luout,*) "using old esp charges"
      end if
 
      e = 0.0d0
      do i=1,nqm
        xp = dbl_mb(i_c_esp  +3*(i-1))
        yp = dbl_mb(i_c_esp+1+3*(i-1))
        zp = dbl_mb(i_c_esp+2+3*(i-1))
        zqm = dbl_mb(i_q_esp+i-1)
        do k = 1,nbq
          zbq = qbq(k)
          xn  = cbq(1,k) - xp
          yn  = cbq(2,k) - yp
          zn  = cbq(3,k) - zp
          rr =  sqrt(xn*xn + yn*yn + zn*zn)
          e = e + zqm*zbq/rr
        end do
      end do

      if (.not.ma_pop_stack(h_q_esp)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)
      if (.not.ma_pop_stack(h_c_esp)) call errquit
     &   (pname//'ma_pop_stack failed',911,MA_ERR)


      end

      subroutine qmmm_esp_input_file(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
c     cannot get mm charges on qm atoms
c     since we zeroed them out in the prepare
c     also cannot do esp input several times in
c     succession for constraints it doest not
c     keep the memory
      integer rtdb
c
      integer i
      integer nfg
      integer nqm
      integer nlink
      integer ail
      integer aig
      integer aiq
      integer un
      integer i_q,h_q
      integer i_tmp,h_tmp,n_tmp
      integer i_grp,h_grp,n_grp
      integer i_grpd,h_grpd,n_grpd
      character*30 filename
      character*30 pname
c
      logical nw_inp_from_file
      external nw_inp_from_file
      logical nw_inp_from_unit
      external nw_inp_from_unit
 
      pname = "qmmm_esp_input_file"
      filename = "esp.nw"
      
      if(qmmm_print_debug()) write(*,*) "in ",pname
      nfg   =  qmmm_geom_get_nfg()
      nqm   =  qmmm_get_nqm()
      nlink =  qmmm_get_nlink()
      ail   =  qmmm_get_link_atom_index()
      aiq   =  qmmm_get_i_ai()
      aig   =  qmmm_geom_get_i_ifg()

      n_grp = nqm
     
      if(qmmm_master()) then
        if(.not.qmmm_get_io_unit(un)) 
     >     call errquit("cannot get file number",0,0)

        open(unit=un,status="scratch",form="formatted")
c        open(unit=un,status="unknown",form="formatted",file=filename)
      end if

      if(qmmm_master())
     >  write(un,899)

      if(qmmm_master())
     >  write(un,900)

      if(.not.ma_push_get(mt_int,nfg,'tmp qmmm index',h_tmp,i_tmp))
     + call errquit(pname//'memory allocation for tmp qmmm index ',
     + nfg, MA_ERR)

      if(.not.ma_push_get(mt_int,n_grp,'tmp qmmm group',h_grp,i_grp))
     + call errquit(pname//'memory allocation for tmp qmmm group ',
     + nfg, MA_ERR)

      if(.not.ma_push_get(mt_int,n_grp,'tmp qmmm group',h_grpd,i_grpd))
     + call errquit(pname//'memory allocation for tmp qmmm group ',
     + nfg, MA_ERR)

      if(.not.ma_push_get(mt_dbl,n_grp,'tmp qmmm charge',h_q,i_q))
     + call errquit(pname//'memory allocation for tmp qmmm group ',
     + nfg, MA_ERR)

      if(.not.qmmm_map(nfg,int_mb(aig),
     >                 nlink,int_mb(ail),
     >                 int_mb(i_tmp)))
     + call errquit(pname//'failed qmmm map ',
     + 0, MA_ERR)

      do i=1,nlink
      if(qmmm_master()) 
     >   write(un,901) int_mb(i_tmp+i-1)
      end do

      call mm_get_solute_charge_gen(nqm,
     >                       mm_quant,
     >                       int_mb(aiq),
     >                       dbl_mb(i_q))

      if(qmmm_master()) then
      write(*,*) "debug 0"
      do i=1,nqm
        write(*,*) int_mb(aiq+i-1),dbl_mb(i_q+i-1)
      end do
      end if
 
      call mm_get_solute_res(nqm,
     >                       mm_quant,
     >                       int_mb(aiq),
     >                       int_mb(i_grp))

      call qmmm_sort_unique(nqm,int_mb(aiq),n_grpd,int_mb(i_grpd))
      if(.not.ma_pop_stack(h_q))
     + call errquit(pname//'memory deallocation for tmp qmmm index ',
     + nfg, MA_ERR)

      if(.not.ma_pop_stack(h_grpd))
     + call errquit(pname//'memory deallocation for tmp qmmm index ',
     + nfg, MA_ERR)

      if(.not.ma_pop_stack(h_grp))
     + call errquit(pname//'memory deallocation for tmp qmmm index ',
     + nfg, MA_ERR)

      if(.not.ma_pop_stack(h_tmp))
     + call errquit(pname//'memory deallocation for tmp qmmm index ',
     + nfg, MA_ERR)

      if(qmmm_master())
     >  write(un,902)

899   FORMAT("esp")
900   FORMAT("recalculate")
901   FORMAT("constrain 0.0 ",I5)
902   FORMAT("end")

      call ga_sync()
c      if(.not.nw_inp_from_file(rtdb, filename))
      if(qmmm_master())
     >  rewind(un)
      if(.not.nw_inp_from_unit(rtdb, un))
     + call errquit(pname//'parsing input file',0,0)
      close(un)

      end

      logical function nw_inp_from_unit(rtdb,un)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb
      integer un
c
      logical input_parse, ignore
      external input_parse
c
      if (ga_nodeid() .eq. 0) then
         call inp_save_state()  ! save state of any inp unit
         call inp_init(un,6)
      endif
c
      ignore = input_parse(rtdb)
c
      if (ga_nodeid() .eq. 0) then
         call inp_restore_state() ! restore state of any inp unit
      endif
c
      nw_inp_from_unit = .true.
      return
c
      end

