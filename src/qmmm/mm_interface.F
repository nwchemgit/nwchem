c
c $Id$
c
      subroutine mm_init(irtdb0)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
c
      integer irtdb0
c
      integer includ
      character*32 pname

      pname = "mm_init"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call timer_init()
      call mm_set_irtdb(irtdb0)
c
c     include fixed atoms
c     otherwise prrroblems
c     -------------------
      includ=1
      if(.not.rtdb_put(irtdb0,'md:includ',mt_int,1,includ))
     + call errquit('mm_init: Failed to set includ',0,
     &       RTDB_ERR)
c
      call md_rdinp()
      call md_start()

c     this call is very important
c
c     for some strange reasons md_sp() has to be called twice
c     to stabilize bqzone flag for waters (MV)
      call md_sp()
      call md_sp()
      if(qmmm_print_debug()) write(*,*) "out "//pname
c
      return
      end

      subroutine mm_set_restart_file(irtdb0,filename)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
c
      integer irtdb0
      character*(*) filename
c
      character*32 pname
      integer j
c
      pname = "mm_set_restart_file"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
c     set the name of the new restart file
c     ------------------------------------
      j = index(filename,'.rst')-1
      if(.not.rtdb_cput(irtdb0,'md:project',1,filename(1:j)))
     + call errquit(pname//'Failed to set md:project',0,
     +       RTDB_ERR)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
c
      return
      end

      subroutine mm_set_top_file(filename)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
c
      character*(*) filename
c
      character*32 pname
      character*255 oldfile
      integer j
c
      pname = "mm_set_restart_file"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c      
      call mm_top_file(oldfile)     
c
c     copy new file over an old one
c     -----------------------------
      call util_file_name_resolve(filename,.false.)
      if(ga_nodeid().eq.0) then 
        write(*,*) "copying filename", filename
        call util_file_copy(filename(1:inp_strlen(filename)),
     >                      oldfile(1:inp_strlen(oldfile)))
      end if
      call ga_sync()
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
c
      return
      end

      subroutine mm_reload_rst(irtdb0,filename)
c     convenience routine to load
c     different restart file into qmmm simulations
c     irtdb0   - rtdb handle              [in]
c     filename - name of the restart file [in]
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
c
      integer irtdb0
      character*(*) filename
c
      character*32 pname
      integer j
      character*255 filename0
      character*255 filcmd,filcmd0
c
      pname = "mm_reload_rst"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
c     get old restart filename
c     ------------------------
      filename0 = ' '
      call mm_restart_filename(filename0)
      j = index(filename0,'.rst')-1
      filcmd0 = filename0(1:j)//".cmd"
c 
c     now clean up all prior md stuff
c     ---------------------------------
      call md_finish()
c
c     set the name of the new restart file
c     ------------------------------------
      j = index(filename,'.rst')-1
      if(.not.rtdb_cput(irtdb0,'md:project',1,filename(1:j)))
     + call errquit(pname//'Failed to set md:project',0,
     +       RTDB_ERR)
      filcmd = filename(1:j)//".cmd"
c
c     make copy of cmd file
c     ---------------------
      if(ga_nodeid().eq.0) then
        call util_file_copy(filcmd0(1:inp_strlen(filcmd0)),
     >                      filcmd(1:inp_strlen(filcmd)))
      end if
      call ga_sync()
c
c     now initialize md structures again
c     ----------------------------------
      call mm_init(irtdb0)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
c
      return
      end

      subroutine mm_reload_rst_default(irtdb0,filename)
c     convenience routine to load
c     restart file in perm directory into default restart file
c     irtdb0   - rtdb handle              [in]
c     filename - name of the restart file in the perm directory[in]
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
#include "inp.fh"
c
      integer irtdb0
      character*(*) filename
c
      character*32 pname
      character*255 filename0
      integer j
c
      pname = "mm_reload_rst_default"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
c     get default restart filename
c     ----------------------------
      filename0 = ' '
      call mm_restart_filename(filename0)
c
c     copy restart file to a default one
c     ---------------------------------
      call util_file_name_resolve(filename,.false.)
      if(ga_nodeid().eq.0) then 
        write(*,*) "copying filename", filename
        call util_file_copy(filename(1:inp_strlen(filename)),
     >                      filename0(1:inp_strlen(filename0)))
      end if
      call ga_sync()
      call mm_reload_rst(irtdb0,filename0)
c 
      if(qmmm_print_debug()) write(*,*) "out "//pname
c
      return
      end

      subroutine mm_write_restart()
      implicit none
c
#include "util.fh"
#include "qmmm.fh"
#include "inp.fh"
c
      character*255 rfile
      integer lfnrst
      character*32 pname
c
      pname = "mm_write_restart"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      rfile = " "
      call util_file_name("rst", .false., .false., rfile)

      lfnrst = 54
      call md_wrtrst(lfnrst,rfile,.false.)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_restart_filename_prefix(filename)
      implicit none
c
#include "util.fh"
#include "md_common.fh"
#include "inp.fh"
c
      character*(*) filename
      integer j

      j = index(filrst,'.rst')-1
      filename = filrst(1:j)
      end

      subroutine mm_restart_filename(filename)
      implicit none
c
#include "util.fh"
#include "md_common.fh"
#include "inp.fh"
c
      character*(*) filename

      filename = filrst
      end

      subroutine mm_trajectory_filename(filename)
      implicit none
c
#include "util.fh"
#include "md_common.fh"
#include "inp.fh"
c
      character*(*) filename

      filename = filtrj
      end

      subroutine mm_system_prefix(prefix)
      implicit none
c
#include "util.fh"
#include "md_common.fh"
#include "inp.fh"
c
      character*(*) prefix
      integer i

      i=index(filtop,'.top')-1
      prefix = filtop(1:i)
      end

      subroutine mm_top_file(filename)
      implicit none
c
#include "util.fh"
#include "md_common.fh"
#include "inp.fh"
c
      character*(*) filename

      filename = filtop
      end

      subroutine mm_set_qmmm(oqmmm)
      implicit none
c
#include "util.fh"
#include "md_common.fh"
#include "inp.fh"
c
      logical oqmmm
      lqmmm = oqmmm

      end

      subroutine mm_write_restart_named(rfile,lv)
      implicit none
c
#include "util.fh"
#include "qmmm.fh"
#include "inp.fh"
c
      character*(*) rfile
      logical lv
c
      integer lfnrst
      character*32 pname
c
      pname = "mm_write_restart_named"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      lfnrst = 54
      call md_wrtrst(lfnrst,rfile,lv)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_end()
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
      character*32 pname
c
      pname = "mm_end"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c      call md_wrtrst(lfnrst,rfile,.false.)
c      call md_wrtrst(lfnqrs,filqrs,.false.)
      call md_finish()
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_write_restart_default()
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
      character*32 pname
c
      pname = "mm_write_restart_default"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      call md_wrtrst(lfnrst,rfile,.false.)
      call md_wrtrst(lfnqrs,filqrs,.false.)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_write_restart_default_dynamics()
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
      character*32 pname
c
      pname = "mm_write_restart_default_dynamics"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      call md_wrtrst(lfnrst,rfile,.true.)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_rewind_trajectory()
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
      rewind(lfntrj)
      return
      end

      function mm_skip_frame(fn,n)
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
      integer n,fn
      logical mm_skip_frame
      character*32 pname
c
      external sp_skip
      logical  sp_skip
c
      pname = "mm_skip_frame"
c
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      mm_skip_frame = sp_skip(fn,n)
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      function mm_read_frame(fn,n)
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer n,fn
      logical mm_read_frame
      character*32 pname
c
      external sp_skip,sp_rdtrj
      logical  sp_skip,sp_rdtrj
c
      pname = "mm_read_frame"
c
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      if(n.gt.1) then
        mm_read_frame = sp_skip(fn,n-1)
        if(.not.mm_read_frame) return
      end if
      call ga_sync()
c 
c     note that lxw,..., are output flags indicating whether
c     corresponding quantites (coords, vel, ..)
c     were found/read from file
c     ------------------------------------------------------
      mm_read_frame = sp_rdtrj(fn,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs),nwmloc,nsaloc)
c
      call ga_brdcst(msg_qmmm_misc,mm_read_frame,
     >               ma_sizeof(mt_log,1,mt_byte),0)
      call ga_sync()

      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_write_frame(fn)
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
      integer fn
      character*32 pname
c
      pname = "mm_write_frame"
c
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
c     atom redistribution
c
      call sp_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)

c      call md_sp()
c
c     here we writing only solvent and solute coords
c     we could also write out velocities,forces, etc
c     ----------------------------------------------
      call sp_wrttrj(fn,.true.,.false.,.false.,.true.,.false.,.false.,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))

      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_print_info()
      implicit none
#include "qmmm.fh"
#include "md_common.fh"
#include "mafdecls.fh"
c
      call prp_print()
c
c     print energies
c
      call cf_print_energy(lfnout)
c
      call sp_printf(filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs),npener,dbl_mb(i_esa))
c
      return
      end

      subroutine mm_task_sp()
      implicit none
#include "qmmm.fh"
#include "md_common.fh"
      character*32 pname
      pname = "mm_task_sp"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
c      lpair = .false.
c      llong = .false.
      call md_sp_qmmm()
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_get_tot_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nqm
      character*32 pname
      pname = "mm_get_tot_nqm"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solute_tot_na_gen(nqm,mm_quant)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_tot_nact(nact)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nact,nacts,nactw
      character*32 pname
      pname = "mm_get_tot_nact"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solute_tot_nact(nacts)
      call mm_get_solvent_tot_nact(nactw)
      nact = nacts + nactw
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_tot_nact(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nqm
      character*32 pname
      pname = "mm_get_tot_nqm"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solute_tot_na_gen(nqm,mm_active)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_tot_nact(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nqm
      character*32 pname
      pname = "mm_get_tot_nqm"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solvent_tot_na_gen(nqm,mm_active_w)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_tot_nact_bqs(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nqm
      character*32 pname
      pname = "mm_get_tot_nqm"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solute_tot_na_gen(nqm,mm_active_bqs)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_tot_nlink(n)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer n
      character*32 pname
      pname = "mm_get_tot_nlink"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solute_tot_na_gen(n,mm_link)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nqm
      character*32 pname
      pname = "mm_get_nqm"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c

      call mm_get_solute_na_gen(nqm,mm_quant_link)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_tot_na_gen(na,mm_type)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "qmmm.fh"
c
      integer na
      logical mm_type
      external mm_type
      integer tmp(1)
      character*32 pname
      pname = "mm_get_solute_tot_na_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solute_na_gen(tmp(1),mm_type)
      call ga_igop(msg_qmmm_nqm,tmp,1,'+')
      na=tmp(1)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_na_gen(na,mm_type)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer i
      logical mm_type
      external mm_type

      integer na
      integer nsaloc
      character*32 pname

      pname = "mm_get_solute_na_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsaloc)
c
      na=0
      do i=1,nsaloc
        if(mm_type(i)) then
          na=na+1
        endif
      end do
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_tot_nbqw(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "qmmm.fh"
c
      integer nbq
      integer tmp(1)
      character*32 pname
      pname = "mm_get_tot_nbqw"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_nbqw(tmp(1))
      call ga_igop(msg_qmmm_nbqw,tmp,1,'+')
      nbq=tmp(1)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_tot_na_gen(na,mm_type)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "qmmm.fh"
c
      integer na
      logical mm_type
      external mm_type

      integer tmp(1)
      character*32 pname
      pname = "mm_get_solvent_na_gen"

      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_solvent_na_gen(tmp(1),mm_type)
      call ga_igop(msg_qmmm_nbqw,tmp,1,'+')
      na=tmp(1)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_na_gen(na,mm_type)
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer na
      logical mm_type
      external mm_type
c
      integer i
      integer mwa
      integer nwmloc
      character*32 pname
      pname = "mm_get_solvent_na_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_mwa(mwa)
      call mm_get_nwmloc(nwmloc)

      na=0
      do i=1,nwmloc
        if(mm_type(i)) na=na+mwa
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_nbqw(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nbq
c
      integer i
      integer mwa
      integer nwmloc
      character*32 pname
      pname = "mm_get_nbqw"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_mwa(mwa)
      call mm_get_nwmloc(nwmloc)

      nbq=0
      do i=1,nwmloc
        if(mm_bqwzone(i)) nbq=nbq+mwa
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_tot_nbqs(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nbq
      character*32 pname
      pname = "mm_get_tot_nbqs"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      call mm_get_solute_tot_na_gen(nbq,mm_bqzone)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_activate_bqszone()
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
c
      integer i
      integer nsaloc
      integer nwmloc
      integer psz
      character*32 pname
      pname = "mm_activate_bqszone"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsaloc)
      call mm_get_nwmloc(nwmloc)
      call mm_get_psz(psz)

      do i=1,nsaloc
          int_mb(psz+i-1)=1
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_activate_bqwzone()
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
c
      integer i
      integer nsaloc
      integer nwmloc
      integer pwz
      character*32 pname
      pname = "mm_activate_bqwzone"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsaloc)
      call mm_get_nwmloc(nwmloc)
      call mm_get_pwz(pwz)

      do i=1,nwmloc
          int_mb(pwz+i-1)=1
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_prune_bqzone()
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
c
      integer i
      integer nsaloc
      integer nwmloc
      integer bq_exclude
      character*32 pname
      pname = "mm_prune_bqzone"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsaloc)
      call mm_get_nwmloc(nwmloc)

      bq_exclude = qmmm_get_bq_exclude()

      if(bq_exclude.eq.qmmm_all_atoms) then
        do i=1,nsaloc
         call mm_clear_bqzone(i)
        end do
        do i=1,nwmloc
          call mm_clear_bqwzone(i)
        end do
      else if(bq_exclude.eq.qmmm_linkbond) then
        do i=1,nsaloc
          if(mm_lqsoft(i)) 
     >      call mm_clear_bqzone(i)
        end do
      else if(bq_exclude.eq.qmmm_linkbond_H) then
        do i=1,nsaloc
          if(mm_lqsoft(i).and.(mm_element(i).eq.1)) 
     >      call mm_clear_bqzone(i)
        end do
      end if
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_quant_link_ind(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_quant_link_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_quant_link,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_link_ind(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_link_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_link,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_bq_ind(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_bq_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_bqzone,ind)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_quant_ind(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_quant_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_quant,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_ind_act(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_ind_act"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_active,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_ind_act(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solvent_ind_act"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solvent_ind_gen(nt,mm_active_w,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_act_bqs_ind(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_act_bq_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_active_bqs,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_info_gen(nt,mm_type,ind,ires)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      logical mm_type
      integer ind(nt)
      integer ires(nt)
      external mm_type

c     local variables
      integer i
      integer ioff,ioq
      integer me
      integer np
      integer nqa
      integer l_ioff,i_ioff
      integer nsa
      integer psl
      character*32 pname
      pname = "mm_get_solute_ind_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call ifill(nt,0,ind,1)

      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)

      np=ga_nnodes()

      if(.not.ma_push_get(mt_int,np,'ioff',l_ioff,i_ioff))
     + call errquit('qmmm: Failed to allocate memory for ioff',np,
     &       MA_ERR)

      call ifill(np,0,int_mb(i_ioff),1)

c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
c      call mm_get_nqm(nqa)
      call mm_get_solute_na_gen(nqa,mm_type)
      do i=me+1,np
        int_mb(i_ioff+i-1)=nqa
      end do    
      call ga_igop(msg_qmmm_offset,int_mb(i_ioff),np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioff=int_mb(i_ioff+me-1)
      else
        ioff=0
      end if

      ioq=ioff
      do i=1,nsa
      if(mm_type(i)) then
         ioq=ioq+1
         ind(ioq)=int_mb(psl+i-1)
      endif
      end do
       
      call ga_igop(msg_qmmm_ind,ind,nt,'+')

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',np,
     &       MA_ERR)

      call qmmm_sort(nt,ind)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_ind_gen(nt,mm_type,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      logical mm_type
      integer ind(nt)
      external mm_type

c     local variables
      integer i
      integer ioff,ioq
      integer me
      integer np
      integer nqa
      integer l_ioff,i_ioff
      integer nsa
      integer psl
      character*32 pname
      pname = "mm_get_solute_ind_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call ifill(nt,0,ind,1)

      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)

      np=ga_nnodes()

      if(.not.ma_push_get(mt_int,np,'ioff',l_ioff,i_ioff))
     + call errquit('qmmm: Failed to allocate memory for ioff',np,
     &       MA_ERR)

      call ifill(np,0,int_mb(i_ioff),1)

c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
c      call mm_get_nqm(nqa)
      call mm_get_solute_na_gen(nqa,mm_type)
      do i=me+1,np
        int_mb(i_ioff+i-1)=nqa
      end do    
      call ga_igop(msg_qmmm_offset,int_mb(i_ioff),np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioff=int_mb(i_ioff+me-1)
      else
        ioff=0
      end if

      ioq=ioff
      do i=1,nsa
      if(mm_type(i)) then
         ioq=ioq+1
         ind(ioq)=int_mb(psl+i-1)
      endif
      end do
       
      call ga_igop(msg_qmmm_ind,ind,nt,'+')

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',np,
     &       MA_ERR)

      call qmmm_sort(nt,ind)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_ind_gen(nt,mm_type,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      logical mm_type
      integer ind(nt)
      external mm_type

c     local variables
      integer i
      integer ioff,ioq
      integer me
      integer np
      integer nqa
      integer l_ioff,i_ioff
      integer nwm
      integer mwa
      integer pwl
      character*32 pname
      pname = "mm_get_solvent_ind_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call ifill(nt,0,ind,1)

      call mm_get_nwmloc(nwm)
      call mm_get_mwa(mwa)

      call mm_get_pwl(pwl)

      np=ga_nnodes()

      if(.not.ma_push_get(mt_int,np,'ioff',l_ioff,i_ioff))
     + call errquit('qmmm: Failed to allocate memory for ioff',np,
     &       MA_ERR)

      call ifill(np,0,int_mb(i_ioff),1)

c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
      call mm_get_solvent_na_gen(nqa,mm_type)
c      call mm_get_nbqw(nqa)
c     number of water molecules = 
c     number of atoms/number of atoms per molecule 
      nqa=nqa/mwa
      do i=me+1,np
        int_mb(i_ioff+i-1)=nqa
      end do    
      call ga_igop(msg_qmmm_nbqw_off,int_mb(i_ioff),np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioff=int_mb(i_ioff+me-1)
      else
        ioff=0
      end if

      ioq=ioff
      do i=1,nwm
      if(mm_type(i)) then
         ioq=ioq+1
         ind(ioq)=int_mb(pwl+i-1)
      endif
      end do

      call ga_igop(msg_qmmm_ind_solvent,ind,nt,'+')

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',np,
     &       MA_ERR)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_ind_bq(nt,ind)
      implicit none
#include "mm_utils.fh"
#include "qmmm.fh"
      integer nt
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solvent_ind_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c

      call mm_get_solvent_ind_gen(nt,mm_bqwzone,ind)
    
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_coord_quant(nt,
     >                       link,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      logical link
      integer ai(nt)
      double precision c(3,nt)
      character*32 pname
      pname = "mm_get_solute_coord_quant"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      if(link) then
        call mm_get_solute_coord_gen(nt,
     >                         mm_quant_link,
     >                         ai,
     >                         c)
      else
        call mm_get_solute_coord_gen(nt,
     >                         mm_quant,
     >                         ai,
     >                         c)
      end if

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_geom(nt,
     >                       ai,
     >                       inum,
     >                       t,
     >                       c,
     >                       q,
     >                       m)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      integer inum(nt)
      character*16 t(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_geom_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)
      call dfill(nt,0.0d0,m,1)
      call ifill(nt,0,inum,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
        do j=1,nt
          if(int_mb(psl+i-1).eq.ai(j)) then
             do k=1,3
               c(k,j)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
             end do
             inum(j)=mm_element(i) 
             q(j) = dble(inum(j))
             if(.not.geom_atn_to_default_mass(inum(j),m(j)))
     >         call errquit(pname,0, GEOM_ERR)
             if(mm_bqzone(i).or.mm_link(i)) then
c               inum(j) = 0
               q(j) = mm_solute_charge(i)
             end if
          end if
        end do
      end do

      call ga_dgop(msg_qmmm_c,c,3*nt,'+')
      call ga_dgop(msg_qmmm_chg,q,nt,'+')
      call ga_igop(msg_qmmm_inum,inum,nt,'+')
      call ga_dgop(msg_qmmm_m,m,nt,'+')

      do i=1,nt
         call cf_num2tag(inum(i),t(i))
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_geom_gen(nt,
     >                       mm_type,
     >                       ai,
     >                       inum,
     >                       t,
     >                       c,
     >                       q)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      logical mm_type
      integer ai(nt)
      integer inum(nt)
      character*16 t(nt)
      double precision c(3,nt)
      double precision q(nt)
      external mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_geom_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)
      call ifill(nt,0,inum,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)


      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 do k=1,3
                   c(k,j)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
                 end do
                 inum(j)=mm_element(i)
                 q(j) = mm_solute_charge(i)
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_c,c,3*nt,'+')
      call ga_dgop(msg_qmmm_chg,q,nt,'+')
      call ga_igop(msg_qmmm_inum,inum,nt,'+')

      do i=1,nt
         call cf_num2tag(inum(i),t(i))
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_charge_gen(nt,
     >                       mm_type,
     >                       ai,
     >                       q)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      logical mm_type
      integer ai(nt)
      double precision q(nt)
      external mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_geom_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(nt,0.0d0,q,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)


      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 q(j) = mm_solute_charge(i)
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_chg,q,nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end
c
      subroutine mm_get_solute_cons_gen(nt,
     >                       mm_type,
     >                       ai,
     >                       act)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      logical act(nt)
      logical mm_type
      external mm_type

      integer i,j
      integer i_act,l_act


c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_cons_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(.not.ma_push_get(mt_int,nt,'i_act',l_act,i_act))
     + call errquit('qmmm: Failed to allocate memory for i_act',nt,
     &       MA_ERR)
      call ifill(nt,0,int_mb(i_act),1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)


      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 if(mm_fixed(i)) int_mb(i_act+j-1)=1
              end if
            end do
        end if
      end do

      call ga_igop(msg_qmmm_act,int_mb(i_act),nt,'+')


      do i=1,nt
         act(i)=act(i).and.(int_mb(i_act+i-1).eq.0)
      end do

      if(.not.ma_pop_stack(l_act))
     + call errquit('qmmm: Failed to deallocate memory for act',nt,
     &       MA_ERR)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_cons(nt,
     >                       ai,
     >                       act)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      logical act(nt)

      integer i,j
      integer i_act,l_act


c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_cons_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(.not.ma_push_get(mt_int,nt,'i_act',l_act,i_act))
     + call errquit('qmmm: Failed to allocate memory for i_act',nt,
     &       MA_ERR)
      call ifill(nt,0,int_mb(i_act),1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)


      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              if(mm_fixed(i)) int_mb(i_act+j-1)=1
           end if
         end do
      end do

      call ga_igop(msg_qmmm_act,int_mb(i_act),nt,'+')


      do i=1,nt
         act(i)=int_mb(i_act+i-1).eq.0
      end do

      if(.not.ma_pop_stack(l_act))
     + call errquit('qmmm: Failed to deallocate memory for act',nt,
     +       MA_ERR)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_fix_solute_region(nt,mm_type,ai)
      implicit none

#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      logical mm_type
      external mm_type
      integer nt
      integer ai(nt)

      integer i,j

c     mm local variables
      integer nsa
      integer psdt
      integer psl
      character*32 pname
      pname = "mm_fix_solute_region"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)

      do i=1,nsa
        if(mm_type(i)) then
        if(int_mb(psl+i-1).eq.ai(j)) then
        do j=1,nt
           int_mb(psdt+i-1)=ior(int_mb(psdt+i-1),lfixed)
        end do
        end if
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_free_solute_region(nt,mm_type,ai)
      implicit none

#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      logical mm_type
      external mm_type
      integer nt
      integer ai(nt)

      integer i,j

c     mm local variables
      integer nsa
      integer psdt
      integer psl
      character*32 pname
      pname = "mm_fix_solute_region"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)

      do i=1,nsa
        if(mm_type(i)) then
        if(int_mb(psl+i-1).eq.ai(j)) then
        do j=1,nt
           int_mb(psdt+i-1)=iand(int_mb(psdt+i-1),not(mfixed))
        end do
        end if
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_fix_solute(mm_type)
      implicit none

#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      logical mm_type
      external mm_type

      integer i

c     mm local variables
      integer nsa
      integer psdt
      character*32 pname
      pname = "mm_fix_solute"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)

      do i=1,nsa
        if(mm_type(i)) then
           int_mb(psdt+i-1)=ior(int_mb(psdt+i-1),lfixed)
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_free_solute(mm_type)
      implicit none

#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      logical mm_type
      external mm_type

      integer i

c     mm local variables
      integer nsa
      integer psdt
      character*32 pname
      pname = "mm_free_solute"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)

      do i=1,nsa
        if(mm_type(i)) then
           int_mb(psdt+i-1)=iand(int_mb(psdt+i-1),not(mfixed))
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_fix_solute_res(nr,resid)
      implicit none

#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      integer nr
      integer resid(nr)

c     mm local variables
      integer i,ir
      integer psg
      integer nsa
      integer psdt
      character*32 pname
      pname = "mm_fix_solute_res"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psg(psg)

      do i=1,nsa
        do ir=1,nr
        if(int_mb(psg+i-1).eq.resid(ir)) then
           int_mb(psdt+i-1)=ior(int_mb(psdt+i-1),lfixed)
        end if
        end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_free_solute_res(nr,resid)
      implicit none

#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      integer nr
      integer resid(nr)

c     mm local variables
      integer i,ir
      integer psg
      integer nsa
      integer psdt
      character*32 pname
      pname = "mm_free_solute_res"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psg(psg)

      do i=1,nsa
        do ir=1,nr
        if(int_mb(psg+i-1).eq.resid(ir)) then
           int_mb(psdt+i-1)=iand(int_mb(psdt+i-1),not(mfixed))
        else
           int_mb(psdt+i-1)=ior(int_mb(psdt+i-1),lfixed)
        end if
        end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_coord_gen(nt,
     >                       mm_type,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      logical mm_type
      integer ai(nt)
      double precision c(3,nt)
      external mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_coord_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(3*nt,0.0d0,c,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 do k=1,3
                   c(k,j)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
                 end do
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_c1,c,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_group(nt,
     >                       mm_type,
     >                       ai,
     >                       grp)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      integer grp(nt)
      external mm_type
      logical mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer psg
      character*32 pname
      pname = "mm_get_solute_group"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call ifill(nt,0,grp,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psg(psg)

      do i=1,nsa
        if(mm_type(i)) then
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              grp(j)=int_mb(psg+i-1)
           end if
         end do
        end if
      end do

      call ga_igop(msg_qmmm_c1,grp,nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_res(nt,
     >                       mm_type,
     >                       ai,
     >                       sgm)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      integer sgm(nt)
      external mm_type
      logical mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer psgm
      character*32 pname
      pname = "mm_get_solute_group"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call ifill(nt,0,sgm,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psgm(psgm)

      do i=1,nsa
        if(mm_type(i)) then
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              sgm(j)=int_mb(psgm+i-1)
           end if
         end do
        end if
      end do

      call ga_igop(msg_qmmm_c1,sgm,nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_coord(nt,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_coord_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(3*nt,0.0d0,c,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                c(k,j)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
              end do
           end if
         end do
      end do

      call ga_dgop(msg_qmmm_c1,c,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_coord_raw(nt,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_coord_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(3*nt,0.0d0,c,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                c(k,j)=dbl_mb(pxs+(i-1)+(k-1)*msa)
              end do
           end if
         end do
      end do

      call ga_dgop(msg_qmmm_c1,c,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_print_solute_coord()
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"


      integer i,k

c     mm local variables
      integer msa
      integer nsa
      integer pxs
      character*32 pname
      pname = "mm_print_solute_coord"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_pxs(pxs)
      
      if(qmmm_master()) then
        write(*,*) "mm_print_solute_coord"
        do i=1,nsa
            write(*,300) mm_solute_charge(i),
     >        (dbl_mb(pxs+(i-1)+(k-1)*msa)*
     >                    cnm2au*cau2ang,k=1,3)
        end do
      end if

300   FORMAT(4F12.6)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_geom(nt,
     >                       ai,
     >                       inum,
     >                       t,
     >                       c,
     >                       q,
     >                       mw)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision mw(nt)
      integer inum(nt)
      character*16 t(nt)
      double precision cf_wcharge
      external cf_wcharge

      integer i,j,k,m,j1

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pxw
      integer pwl
      character*32 pname
      pname = "mm_get_solvent_geom"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)
      call dfill(nt,0.0d0,mw,1)
      call ifill(nt,0,inum,1)

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pxw(pxw)
      call mm_get_pwl(pwl)

 
      do i=1,nwm
         do j=1,nt/mwa
           if(int_mb(pwl+i-1).eq.ai(j)) then
              do m=1,mwa
                j1=(j-1)*mwa+m
                do k=1,3
                  c(k,j1)=dbl_mb(pxw+(i-1)+(k-1)*mwm+(m-1)*3*mwm)
                  c(k,j1)=c(k,j1)*cnm2au
                end do
                inum(j1)=0
                q(j1) = cf_wcharge(m)
                mw(j1) = mm_ww(m)
              end do
           end if
         end do
      end do


      call ga_dgop(msg_qmmm_cbqw,c,3*nt,'+')
      call ga_dgop(msg_qmmm_qbqw,q,nt,'+')
      call ga_dgop(msg_qmmm_mw,mw,nt,'+')
      call ga_igop(msg_qmmm_inumw,inum,nt,'+')

      do i=1,nt
          t(i) = "Bqw"
c         call cf_num2tag(inum(i),t(i))
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_geom_bq(nt,
     >                       ai,
     >                       c,
     >                       q)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      logical mm_type
      integer ai(nt)
      double precision c(3,nt)
      double precision q(nt)
      external mm_type
      double precision cf_wcharge
      external cf_wcharge

      integer i,j,k,m,j1

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pxw
      integer pwl
      character*32 pname
      pname = "mm_get_solvent_geom"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pxw(pxw)
      call mm_get_pwl(pwl)

 
      do i=1,nwm
         do j=1,nt/mwa
           if(int_mb(pwl+i-1).eq.ai(j)) then
              do m=1,mwa
                j1=(j-1)*mwa+m
                do k=1,3
                  c(k,j1)=dbl_mb(pxw+(i-1)+(k-1)*mwm+(m-1)*3*mwm)
                  c(k,j1)=c(k,j1)*cnm2au
                end do
                q(j1) = cf_wcharge(m)
              end do
           end if
         end do
      end do


      call ga_dgop(msg_qmmm_cbqw,c,3*nt,'+')
      call ga_dgop(msg_qmmm_qbqw,q,nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_geom_bq(nt,
     >                       ai,
     >                       c,
     >                       q)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)
      double precision q(nt)

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_get_solute_geom_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)


      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                c(k,j)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
              end do
              q(j) = mm_solute_charge(i)
           end if
         end do
      end do

      call ga_dgop(msg_qmmm_c,c,3*nt,'+')
      call ga_dgop(msg_qmmm_chg,q,nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solvent_coord(nt,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k,m,j1

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pxw
      integer pwl
      character*32 pname
      pname = "mm_get_solvent_geom"
      if(qmmm_print_debug()) write(*,*) "in "//pname



      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pxw(pxw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
c        if(mm_bqwzone(i)) then
            do j=1,nt/mwa
              if(int_mb(pwl+i-1).eq.ai(j)) then
                 do m=1,mwa
                   j1=(j-1)*mwa+m
                   do k=1,3
                     dbl_mb(pxw+(i-1)+(k-1)*mwm+(m-1)*3*mwm)=c(k,j1)*
     >                                                       cau2nm
                   end do
                 end do
              end if
            end do
c        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_coord(nt,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k,m,j1

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pxw
      integer pwl
      character*32 pname
      pname = "mm_get_solvent_coord"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,c,1)

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pxw(pxw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
         do j=1,nt/mwa
           if(int_mb(pwl+i-1).eq.ai(j)) then
              do m=1,mwa
                j1=(j-1)*mwa+m
                do k=1,3
                  c(k,j1) = 
     >             dbl_mb(pxw+(i-1)+(k-1)*mwm+(m-1)*3*mwm)/cau2nm
                end do
              end do
           end if
         end do
      end do

      call ga_dgop(msg_qmmm_cw,c,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_cons(nt,
     >                       ai,
     >                       act)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      logical act(nt)

      integer i,j,m

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pxw
      integer pwl
      integer i_act,l_act

      character*32 pname
      pname = "mm_get_solvent_cons_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(.not.ma_push_get(mt_int,nt,'i_act',l_act,i_act))
     + call errquit('qmmm: Failed to allocate memory for i_act',nt,
     &       MA_ERR)
      call ifill(nt,0,int_mb(i_act),1)


      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pxw(pxw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
        if(mm_fixed_w(i)) then
            do j=1,nt/mwa
              if(int_mb(pwl+i-1).eq.ai(j)) then
                 do m=1,mwa
                 int_mb(i_act+(j-1)*mwa+m-1)=1
                 end do
              end if
            end do
        end if
      end do

      call ga_igop(msg_qmmm_wact,int_mb(i_act),nt,'+')

      do i=1,nt
         act(i)=act(i).and.(int_mb(i_act+i-1).eq.0)
      end do

      if(.not.ma_pop_stack(l_act))
     + call errquit('qmmm: Failed to deallocate memory for act',nt,
     &       MA_ERR)


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_fix_solvent()
      implicit none

#include "qmmm_params.fh"
#include "qmmm.fh"
#include "mafdecls.fh"
#include "errquit.fh"

c     mm local variables
      integer nwm
      integer pwdt
      integer i

      character*32 pname
      pname = "mm_get_solvent_cons_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_pwdt(pwdt)
      call mm_get_nwmloc(nwm)
 
      do i=1,nwm
        int_mb(pwdt+i-1)=ior(int_mb(pwdt+i-1),lfixed)
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_force(nt,
     >                       ai,
     >                       g)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision g(3,nt)

      integer i,j,k,m,j1,ioff

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pfw
      integer pwl
      character*32 pname
      pname = "mm_get_solvent_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,g,1)

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pfw(pfw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
c        if(mm_bqwzone(i)) then
            do j=1,nt/mwa
              if(int_mb(pwl+i-1).eq.ai(j)) then
                 do m=1,mwa
                   j1=(j-1)*mwa+m
                   do k=1,3
                     ioff=pfw+(i-1)+(k-1)*mwm+(m-1)*3*mwm
                     g(k,j1)=-dbl_mb(ioff)*cau2nm/cau2kj
                   end do
                 end do
              end if
            end do
c        end if
      end do

      call ga_dgop(msg_qmmm_forcew,g,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_rms_solvent_force_gen(mm_type,grms)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      logical mm_type
      external mm_type
      double precision grms

      integer i,k,m,ioff

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pfw
      integer ntot
      double precision g
      character*32 pname
      pname = "mm_get_rms_solvent_force_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pfw(pfw)
 
      ntot = 0
      grms = 0.0d0
      do i=1,nwm
        if(mm_type(i)) then
          ntot = ntot + 1
            do m=1,mwa
              do k=1,3
                ioff=pfw+(i-1)+(k-1)*mwm+(m-1)*3*mwm
                g=dbl_mb(ioff)
                grms = grms + g*g
              end do
            end do
       end if
      end do

      call ga_dgop(msg_qmmm_forcew,grms,1,'+')
      call ga_dgop(msg_qmmm_nqm,ntot,1,'+')
c
      grms = sqrt(grms/dble(ntot))*cau2nm/cau2kj
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solvent_force(nt,
     >                       ai,
     >                       g)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision g(3,nt)

      integer i,j,k,m,j1,ioff

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pfw
      integer pwl
      character*32 pname
      pname = "mm_set_solvent_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pfw(pfw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
        if(mm_bqwzone(i)) then
            do j=1,nt/mwa
              if(int_mb(pwl+i-1).eq.ai(j)) then
                 do m=1,mwa
                   j1=(j-1)*mwa+m
                   do k=1,3
                     ioff=pfw+(i-1)+(k-1)*mwm+(m-1)*3*mwm
                     dbl_mb(ioff)=-g(k,j1)*cau2kj/cau2nm
                   end do
                 end do
              end if
            end do
        end if
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_zero_solvent_force(nt,
     >                       ai,
     >                       map)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      logical map(nt)

      integer i,j,k,m,j1,ioff

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pfw
      integer pwl
      character*32 pname
      pname = "mm_set_solvent_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pfw(pfw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
         do j=1,nt/mwa
           if(int_mb(pwl+i-1).eq.ai(j)) then
              do m=1,mwa
                j1=(j-1)*mwa+m
                if(.not.map(j1)) then
                do k=1,3
                  ioff=pfw+(i-1)+(k-1)*mwm+(m-1)*3*mwm
                  dbl_mb(ioff)=0.0d0
                end do
                end if
              end do
           end if
         end do
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_print_solvent_force(un)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer un

      integer i,j,k,m,j1,ioff

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pfw
      integer pwl
      character*32 pname
      pname = "mm_print_solvent_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pfw(pfw)
      call mm_get_pwl(pwl)
 
      do i=1,nwm
       do m=1,mwa
         j1=(j-1)*mwa+m
         write(un,*)  
     +    (dbl_mb(pfw+(i-1)+(k-1)*mwm+(m-1)*3*mwm),k=1,3)
       end do
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_add_solvent_force(nt,
     >                       ai,
     >                       g)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision g(3,nt)

      integer i,j,k,m,j1,ioff

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pfw
      integer pwl
      character*32 pname
      pname = "mm_add_solvent_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pfw(pfw)
      call mm_get_pwl(pwl)

      do i=1,nwm
        if(mm_bqwzone(i)) then
            do j=1,nt/mwa
              if(int_mb(pwl+i-1).eq.ai(j)) then
                 do m=1,mwa
                   j1=(j-1)*mwa+m
                   do k=1,3
                     ioff=pfw+(i-1)+(k-1)*mwm+(m-1)*3*mwm
                     dbl_mb(ioff)=dbl_mb(ioff)-g(k,j1)*cau2kj/cau2nm
                   end do
                 end do
              end if
            end do
        end if
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_force_quant(nt,link,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "rtdb.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
      logical link
      character*32 pname
      pname = "mm_get_solute_force_quant"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      if(link) then
      call mm_get_solute_force_gen(nt,mm_quant_link,ai,g)
      else
      call mm_get_solute_force_gen(nt,mm_quant,ai,g)
      end if
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_force_quant(nt,link,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "rtdb.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
      logical link
      character*32 pname

      pname = " mm_set_solute_force_quant"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(qmmm_print_debug())
     >  write(*,*) "in "//pname

      if(link) then
      call mm_set_solute_force_gen(nt,mm_quant_link,ai,g)
      else
      call mm_set_solute_force_gen(nt,mm_quant,ai,g)
      end if

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_force_bq(nt,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "rtdb.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
      character*32 pname
      pname = "mm_get_solute_force_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_force_gen(nt,mm_bqzone,ai,g)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_force_bq(nt,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "rtdb.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
      character*32 pname
      pname = "mm_set_solute_force_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_set_solute_force_gen(nt,mm_bqzone,ai,g)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_force_gen(nt,mm_type,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer nt
      logical mm_type
      integer ai(nt)
      double precision  g(3,nt)
      external mm_type
c     local variables
      integer i,j,k
      integer msa
      integer nsa
      integer psdt
      integer psl
      integer pfs
      character*32 pname

      pname = "mm_set_solute_force_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 do k=1,3
                   dbl_mb(pfs+(i-1)+(k-1)*msa)=-g(k,j)*cau2kj/cau2nm
                 end do
              end if
            end do
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_zero_solute_force(nt,ai,map)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer nt
      integer ai(nt)
      logical map(nt)
c     local variables
      integer i,j,k
      integer msa
      integer nsa
      integer psdt
      integer psl
      integer pfs
      character*32 pname

      pname = "mm_set_solute_force_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      do i=1,nsa
        do j=1,nt
          if(int_mb(psl+i-1).eq.ai(j)) then
             if(.not.map(j)) then
             do k=1,3
               dbl_mb(pfs+(i-1)+(k-1)*msa)=0.0d0
             end do
             end if
          end if
        end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_add_solute_force(nt,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
c     local variables
      integer i,j,k
      integer msa
      integer nsa
      integer psdt
      integer psl
      integer pfs
      character*32 pname

      pname = "mm_add_solute_force_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                dbl_mb(pfs+(i-1)+(k-1)*msa)=
     >             dbl_mb(pfs+(i-1)+(k-1)*msa)-g(k,j)*cau2kj/cau2nm
              end do
           end if
         end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_force(nt,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
c     local variables
      integer i,j,k
      integer msa
      integer nsa
      integer psdt
      integer psl
      integer pfs
      character*32 pname

      pname = "mm_get_solute_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,g,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                g(k,j)=-dbl_mb(pfs+(i-1)+(k-1)*msa)*cau2nm/cau2kj
              end do
           end if
         end do
      end do

      call ga_dgop(msg_qmmm_force,g,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_force_raw(nt,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer nt
      integer ai(nt)
      double precision  g(3,nt)
c     local variables
      integer i,j,k
      integer msa
      integer nsa
      integer psdt
      integer psl
      integer pfs
      character*32 pname

      pname = "mm_get_solute_force"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,g,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                g(k,j)=-dbl_mb(pfs+(i-1)+(k-1)*msa)
              end do
           end if
         end do
      end do

      call ga_dgop(msg_qmmm_force,g,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_rms_solute_force_gen(mm_type,grms)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      logical mm_type
      double precision  grms
      external mm_type
c     local variables
      integer i,k
      integer msa
      integer nsa
      integer pfs
      integer ntot
      character*32 pname

      pname = "mm_get_rms_solute_force_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname
c
      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_pfs(pfs)
c
      ntot = 0
      grms = 0.0d0
      do i=1,nsa
        if(mm_type(i)) then
           ntot = ntot + 1
           do k=1,3
             grms=grms + dbl_mb(pfs+(i-1)+(k-1)*msa)*
     >                   dbl_mb(pfs+(i-1)+(k-1)*msa) 
           end do
        end if
      end do
c
      call ga_dgop(msg_qmmm_force,grms,1,'+')
      call ga_igop(msg_qmmm_nqm,ntot,1,'+')
c
      grms = sqrt(grms/dble(ntot))*cau2nm/cau2kj
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_force_gen(nt,mm_type,ai,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer nt
      logical mm_type
      integer ai(nt)
      double precision  g(3,nt)
      external mm_type
c     local variables
      integer i,j,k
      integer msa
      integer nsa
      integer psdt
      integer psl
      integer pfs
      character*32 pname

      pname = "mm_get_solute_force_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call dfill(3*nt,0.0d0,g,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 do k=1,3
                   g(k,j)=-dbl_mb(pfs+(i-1)+(k-1)*msa)*cau2nm/cau2kj
                 end do
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_force,g,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_coord_quant(nt,
     >                       link,
     >                       ai,
     >                       c)
      implicit none

#include "mm_utils.fh"
#include "qmmm.fh"

      integer nt
      logical link
      integer ai(nt)
      double precision c(3,nt)
      character*32 pname
      pname = "mm_set_solute_coord_quant"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(link) then
      call mm_set_solute_coord_gen(nt,
     >                       mm_quant_link,
     >                       ai,
     >                       c)

      else
      call mm_set_solute_coord_gen(nt,
     >                       mm_quant,
     >                       ai,
     >                       c)


      end if
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_coord_bq(nt,
     >                       ai,
     >                       c)
      implicit none

#include "mm_utils.fh"
#include "qmmm.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)
      character*32 pname
      pname = "mm_set_solute_coord_bq"

      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_set_solute_coord_gen(nt,
     >                       mm_bqzone,
     >                       ai,
     >                       c)

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_coord_gen(nt,
     >                       mm_type,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      logical mm_type
      integer ai(nt)
      double precision c(3,nt)
      external mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_set_solute_coord_gen"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 do k=1,3
                   dbl_mb(pxs+(i-1)+(k-1)*msa)=c(k,j)/cnm2au
                 end do
              end if
            end do
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_coord_raw(nt,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_set_solute_coord"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                dbl_mb(pxs+(i-1)+(k-1)*msa)=c(k,j)
              end do
           end if
         end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_set_solute_coord(nt,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_set_solute_coord"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                dbl_mb(pxs+(i-1)+(k-1)*msa)=c(k,j)/cnm2au
              end do
           end if
         end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_mix_solute_coord(nt,
     >                       lambda,
     >                       ai,
     >                       c)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      double precision lambda
      integer ai(nt)
      double precision c(3,nt)

      double precision x0,x1,x2
      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer pxs
      character*32 pname
      pname = "mm_set_solute_coord"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
         do j=1,nt
           if(int_mb(psl+i-1).eq.ai(j)) then
              do k=1,3
                x0 = dbl_mb(pxs+(i-1)+(k-1)*msa)
                x1 = c(k,j)/cnm2au
                x2=(x1-x0)*lambda+x0
                dbl_mb(pxs+(i-1)+(k-1)*msa)=x2
              end do
           end if
         end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_mix_solute_coord_file(aname,lambda)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) aname    
      double precision lambda
c     local variables
      integer ns
      integer i
      integer k
      logical title


      integer i_itmp,h_itmp
      integer i_ctmp,h_ctmp
      character*32 pname
      character*30 buf
      character*255 filename, dir

      integer fn
      logical master

      pname = "mm_set_solute_coord_file"
      master = qmmm_master()

      filename = aname(1:inp_strlen(aname))
      if(master)  
     +  call util_print_centered(luout,
     + "Loading external geometry from "//
     +   filename,40,.true.)
      dir      = ' '
      call util_directory_name(dir, .false., 0)
      filename = ' '
      write(filename,13) dir(1:inp_strlen(dir)), 
     $     aname(1:inp_strlen(aname))
 13      format(a,'/',a)
      if(master) then  
        if(.not.qmmm_get_io_unit(fn)) 
     >       call errquit("cannot get file number",0,0)
c
        open(fn,file=filename,form='formatted',status='old',
     $          err=133)
      end if
      
c     get number of atoms
      if(master) then
         read(fn,*) ns
      else
         ns = 0
      end if


      call ga_igop(msg_qmmm_nqm,ns,1,'+')
      call ga_sync()

      if(.not.ma_push_get(mt_int,ns,'itmp',h_itmp,i_itmp))
     + call errquit( pname//'Failed to allocate memory for itmp',
     + ns, MA_ERR)
      call ifill(ns,0,int_mb(i_itmp),1)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ns, MA_ERR)
      call dfill(3*ns,0.0d0,dbl_mb(i_ctmp),1)

c
c     read the coords and global index
c     --------------------------------
      if(master) then
c     testing for the presence of the title
        rewind(fn)
        title =.false.
        do i=1,ns+2
          read(fn,*,end=1) buf
        end do
        title = .true.
1       continue
c       title or not read coordinates now
        rewind(fn) 
        read(fn,*) buf
        if(title) 
     +     read(fn,*) buf

        do i=1,ns
         read(fn,*) buf, (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3),
     +               int_mb(i_itmp+i-1)
         write(luout,'(A,3F12.6,I5)') buf, 
     +              (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3),
     +               int_mb(i_itmp+i-1)

        end do
        close(fn)
        do i=1,3*ns
          dbl_mb(i_ctmp+i-1) = dbl_mb(i_ctmp+i-1)/cau2ang
        end do
      end if
c
      call ga_igop(msg_qmmm_ind,int_mb(i_itmp),ns,'+')
      call ga_dgop(msg_qmmm_c1,dbl_mb(i_ctmp),3*ns,'+')
c
      call ga_sync
c
      call mm_mix_solute_coord(ns,
     >                         lambda,
     >                         int_mb(i_itmp),
     >                         dbl_mb(i_ctmp))

c
      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_itmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)

      return

 133  call errquit(pname//'error opening/closing xyz file',0, GEOM_ERR)

      end

      subroutine mm_set_solute_coord_file(aname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) aname    
c     local variables
      integer ns
      integer i
      integer k
      logical title


      integer i_itmp,h_itmp
      integer i_ctmp,h_ctmp
      character*32 pname
      character*30 buf
      character*255 filename, dir

      integer fn
      logical master

      pname = "mm_set_solute_coord_file"
      master = qmmm_master()

      filename = aname(1:inp_strlen(aname))
      if(master)  
     +  call util_print_centered(luout,
     + "Loading external geometry from "//
     +   filename,40,.true.)
      dir      = ' '
      call util_directory_name(dir, .false., 0)
      filename = ' '
      write(filename,13) dir(1:inp_strlen(dir)), 
     $     aname(1:inp_strlen(aname))
 13      format(a,'/',a)
      if(master) then  
        if(.not.qmmm_get_io_unit(fn)) 
     >       call errquit("cannot get file number",0,0)
c
        open(fn,file=filename,form='formatted',status='old',
     $          err=133)
      end if
      
c     get number of atoms
      if(master) then
         read(fn,*) ns
      else
         ns = 0
      end if


      call ga_igop(msg_qmmm_nqm,ns,1,'+')
      call ga_sync()

      if(.not.ma_push_get(mt_int,ns,'itmp',h_itmp,i_itmp))
     + call errquit( pname//'Failed to allocate memory for itmp',
     + ns, MA_ERR)
      call ifill(ns,0,int_mb(i_itmp),1)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ns, MA_ERR)
      call dfill(3*ns,0.0d0,dbl_mb(i_ctmp),1)

c
c     read the coords and global index
c     --------------------------------
      if(master) then
c     testing for the presence of the title
        rewind(fn)
        title =.false.
        do i=1,ns+2
          read(fn,*,end=1) buf
        end do
        title = .true.
1       continue
c       title or not read coordinates now
        rewind(fn) 
        read(fn,*) buf
        if(title) 
     +     read(fn,*) buf

        do i=1,ns
         read(fn,*) buf, (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3),
     +               int_mb(i_itmp+i-1)
         write(luout,*) buf, (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3),
     +               int_mb(i_itmp+i-1)

        end do
        close(fn)
        do i=1,3*ns
          dbl_mb(i_ctmp+i-1) = dbl_mb(i_ctmp+i-1)/cau2ang
        end do
      end if
c
      call ga_igop(msg_qmmm_ind,int_mb(i_itmp),ns,'+')
      call ga_dgop(msg_qmmm_c1,dbl_mb(i_ctmp),3*ns,'+')
c
      call ga_sync
c
      call mm_set_solute_coord(ns,
     >                         int_mb(i_itmp),
     >                         dbl_mb(i_ctmp))

c
      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_itmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)

      return

 133  call errquit(pname//'error opening/closing xyz file',0, GEOM_ERR)

      end

      subroutine mm_create_trj_from_xyz(in_xyz,out_trj)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) in_xyz
      character*(*) out_trj
c     local variables
      integer ns
      integer i
      integer k
      logical title
      integer i_itmp,h_itmp
      integer i_ctmp,h_ctmp
      character*32 pname
      character*30 buf
      character*255 filename
      character*255 xyzfile
      character*255 trjfile
      character*255 message

      integer fn_xyz,fn_trj
      logical end_of_file
      logical master

      pname = "mm_create_trj_from_xyz"
      master = qmmm_master()
c
c     we assume that xyz file has a title
c     -----------------------------------
      title = .true.
c
      xyzfile = in_xyz(1:inp_strlen(in_xyz))
      call util_file_name_resolve(xyzfile, .false.)
c
      trjfile = out_trj(1:inp_strlen(out_trj))
      call util_file_name_resolve(trjfile, .false.)
c
      filename = in_xyz(1:inp_strlen(xyzfile))
      if(master) 
     + call util_print_centered(luout,
     + "Loading external xyz frames from "//
     +   filename,
     +   40,.true.)
c
      filename = out_trj(1:inp_strlen(out_trj))
      if(master) 
     + call util_print_centered(luout,
     + "Creating mm trajectory file "//
     +   filename,
     +   40,.true.)
c
c     prepare files for reading/writing
c     ---------------------------------
      if(master) then  
        if(.not.qmmm_get_io_unit(fn_xyz)) 
     >       call errquit("cannot get file number",0,0)
        filename = xyzfile
        open(fn_xyz,file=filename,form='formatted',status='old',
     $          err=133)
c
        if(.not.qmmm_get_io_unit(fn_trj)) 
     >       call errquit("cannot get file number",0,0)
        filename = trjfile
        open(fn_trj,file=filename,form='formatted',status='unknown',
     $          err=133)
c       write header for trj file
        call cf_trjhdr(fn_trj)
      end if
c
c     get number of atoms
c     ------------------
      message = " number of atoms "
      if(master) then
         read(fn_xyz,*,err=134) ns
         rewind(fn_xyz)
      end if
      call ga_brdcst(msg_qmmm_nqm,
     >     ns,
     >     ma_sizeof(mt_int,1,mt_byte),
     >     0)
      call ga_sync()

c
c     temporary stack memory
c     ----------------------
      if(.not.ma_push_get(mt_int,ns,'itmp',h_itmp,i_itmp))
     + call errquit( pname//'Failed to allocate memory for itmp',
     + ns, MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ns, MA_ERR)

c
c     read the coords and global index
c     --------------------------------
2     continue
      if(master) then
        end_of_file = .true.
        message = " number of atoms field"
        read(fn_xyz,*,err=134,end=3) buf
        end_of_file = .false.
        message = " title field"
        if(title) 
     +     read(fn_xyz,*,err=134,end=135) buf

        do i=1,ns
         read(fn_xyz,*,err=134,end=135) buf, 
     +        (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3),
     +        int_mb(i_itmp+i-1)
        end do
        call dscal(3*ns,1/cau2ang,dbl_mb(i_ctmp),1)
      end if
c
3     continue
c
c     check if end of file was found
c     ------------------------------
      call ga_brdcst(msg_qmmm_misc,end_of_file,
     >               ma_sizeof(mt_log,1,mt_byte),0)
      if(end_of_file) goto 4
c
c     send coords/index to others
c     ---------------------------
      call ga_brdcst(msg_qmmm_ind,
     >     int_mb(i_itmp),
     >     ns*ma_sizeof(mt_int,1,mt_byte),
     >     0)
      call ga_brdcst(msg_qmmm_c1,
     >     dbl_mb(i_ctmp),
     >     3*ns*ma_sizeof(mt_dbl,1,mt_byte),
     >     0)
      call ga_sync()
c
c     push coords to mm
c     -----------------
      call mm_set_solute_coord(ns,
     >                         int_mb(i_itmp),
     >                         dbl_mb(i_ctmp))
c
      call mm_write_frame(fn_trj)
      goto 2
c
 4    continue
c
      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_itmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)

      if(master) then
        close(fn_xyz)
        close(fn_trj)
      end if
      return

 133  call errquit(pname//'error opening/closing '//filename,0, 0)
 134  call errquit(pname//'error reading xyz file'//message,0, 0)
 135  call errquit(pname//'error end of file at'//message,0, 0)

      end

      subroutine mm_debug0()
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"



c     local variables
      integer i
      integer msa
      integer nsa
      integer psl
      integer psz
      character*32 pname
      pname = "mm_debug0"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psz(psz)

      do i=1,nsa
        write(56,*) int_mb(psl+i-1),
     >                int_mb(psz+i-1)

      end do

      call ga_sync()
      if(qmmm_print_debug()) write(*,*) "out "//pname

      end

      subroutine mm_debug()
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"



c     local variables
      integer i
      integer msa
      integer nsa
      integer psl
      integer psz
      character*32 pname
      pname = "mm_debug"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psz(psz)

      do i=1,nsa
        if(int_mb(psl+i-1).eq.639) then
          write(*,*) "found atom",int_mb(psl+i-1),
     >                int_mb(psz+i-1)

        end if
      end do

      call ga_sync()

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_print_system_pdb()
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"



c     local variables
      integer i,k
      integer msa
      integer nsa
      integer psl
      integer psz
      integer pxs
      integer psg
      double precision c(3)
      character*16 t
      character*32 pname
      pname = "mm_print_system"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psz(psz)
      call mm_get_pxs(pxs)
      call mm_get_psg(psg)

      if(ga_nodeid().eq.0) then
      do i=1,nsa
         do k=1,3
c           c(k)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
           c(k)=dbl_mb(pxs+(i-1)+(k-1)*msa)
         end do
         call cf_num2tag(mm_element(i),t)
         write(23,FMT=9000)
     >            int_mb(psl+i-1), 
     >            t,
     >            int_mb(psg+i-1),
c     >            c(1)*0.529177249d00,
c     >            c(2)*0.529177249d00,
c     >            c(3)*0.529177249d00,
     >            c(1),
     >            c(2),
     >            c(3),
     >            int_mb(psz+i-1),
     >            mm_quant(i),
     >            mm_fixed(i)

      end do
      end if

9000  FORMAT("ATOM",T7,I5,T13,A4,T23,I4,T31,
     >       F8.3,T39,F8.3,T47,F8.3,T55,I2,
     >       T60,L2,T62,L2)


      call ga_sync()

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_print_system()
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"



c     local variables
      integer i,k
      integer msa
      integer nsa
      integer psl
      integer psz
      integer pxs
      integer psg
      double precision c(3)
      character*16 t
      character*32 pname
      pname = "mm_print_system"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psz(psz)
      call mm_get_pxs(pxs)
      call mm_get_psg(psg)

      if(ga_nodeid().eq.0) then
      do i=1,nsa
         do k=1,3
c           c(k)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
           c(k)=dbl_mb(pxs+(i-1)+(k-1)*msa)
         end do
         call cf_num2tag(mm_element(i),t)
         write(23,FMT=9000)
     >            int_mb(psl+i-1), 
     >            t,
     >            int_mb(psg+i-1),
c     >            c(1)*0.529177249d00,
c     >            c(2)*0.529177249d00,
c     >            c(3)*0.529177249d00,
     >            c(1),
     >            c(2),
     >            c(3),
     >            int_mb(psz+i-1),
     >            mm_quant(i),
     >            mm_fixed(i)

      end do
      end if

9000  FORMAT("ATOM",T7,I5,T13,A4,T23,I4,T31,
     >       F12.6,4X,F12.6,4X,F12.6,4X,I2,
     >       4x,L2,4x,L2)


      call ga_sync()

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_print_system_file(fn)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"


      integer fn

c     local variables
      integer i,k
      integer msa
      integer nsa
      integer psl
      integer psz
      integer pxs
      integer psg
      double precision c(3)
      character*16 t
      character*32 pname
      pname = "mm_print_system"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psz(psz)
      call mm_get_pxs(pxs)
      call mm_get_psg(psg)

      if(ga_nodeid().eq.0) then
      do i=1,nsa
         do k=1,3
           c(k)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
         end do
         call cf_num2tag(mm_element(i),t)
         write(fn,FMT=9000)
     >            int_mb(psl+i-1), 
     >            t,
     >            int_mb(psg+i-1),
     >            c(1)*0.529177249d00,
     >            c(2)*0.529177249d00,
     >            c(3)*0.529177249d00,
     >            int_mb(psz+i-1),
     >            mm_quant(i),
     >            mm_fixed(i)

      end do
      end if

9000  FORMAT("ATOM",T7,I5,T13,A4,T23,I4,T31,
     >       F8.3,T39,F8.3,T47,F8.3,T55,I2,
     >       T60,L2,T62,L2)


      call ga_sync()

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_links(na,ai,nlinks,ilinks)
      implicit none
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mafdecls.fh"

      integer na
      integer ai(na)
      integer nlinks(na)
      integer ilinks(na,max_nb+1)
c
c     local variables:
      integer i,in
      integer psb
      integer msb
      integer ind1,ind2,ind0
      character*32 pname
      pname = "mm_links"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call cf_get_msb(msb)
      call cf_get_psb(psb)

      do in=1,na
       nlinks(in) = 0
      end do

      do i=1,msb
        ind1=int_mb(psb+i-1)
        ind2=int_mb(psb+msb+i-1)
        write(*,*) ind1,ind2
        do in=1,na
          ind0=ilinks(in,1)
          if(ind0.eq.ind1) then
             nlinks(in)=nlinks(in)+1
             ilinks(in,nlinks(in)+1)=ind2
          else if(ind0.eq.ind2) then
             nlinks(in)=nlinks(in)+1
             ilinks(in,nlinks(in)+1)=ind1
          end if
        end do
      end do
      write(*,*) "link atom list"

      do in=1,na
        write(*,*) nlinks(in), (ilinks(in,i),i=1,max_nb+1)
      end do
      
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
 
      end

      subroutine mm_count_links(nl,h_lb,i_lb)
      implicit none
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "errquit.fh"

      integer nl
      integer h_lb
      integer i_lb
c     local variables:
      integer i
      integer psb
      integer msb
      integer ind1,ind2
      integer ind_qm,ind_l
      integer iqm
      integer nqm
      integer h_qind,i_qind
      integer il,nl0
      integer h_lind,i_lind
      integer h_link,i_link
      character*32 pname
      pname = "mm_count_links"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_tot_nlink(nl0)
      if(qmmm_print_debug()) then
        write(*,*) pname,"number of link atoms",nl0
      end if
      if(nl0.eq.0)  then
        nl=nl0
        return
      end if

      call mm_get_tot_nqm(nqm)
      if(.not.ma_push_get(mt_int,nqm,'tmp qind',h_qind,i_qind))
     + call errquit(pname//'Failed to allocate stack',nqm,
     &       MA_ERR)

      call mm_get_solute_ind_gen(nqm,mm_quant,int_mb(i_qind))

      if(.not.ma_push_get(mt_int,nl0,'tmp lind',h_lind,i_lind))
     + call errquit(pname//'Failed to allocate stack',nl,
     &       MA_ERR)


      call mm_get_solute_ind_gen(nl0,mm_link,int_mb(i_lind))

      if(.not.ma_push_get(mt_int,nl0*max_nb,'tmp links',h_link,i_link))
     + call errquit(pname//'Failed to allocate stack',nl,
     &       MA_ERR)
      call ifill(nl0*max_nb,-1,int_mb(i_link),1)


      call cf_get_msb(msb)
      call cf_get_psb(psb)


      nl = 0
      do i=1,msb
        ind1=int_mb(psb+i-1)
        ind2=int_mb(psb+msb+i-1)
        do il=1,nl0
          ind_l = int_mb(i_lind+il-1)
          if((ind_l.eq.ind1) .or. (ind_l.eq.ind2)) then
            do iqm=1,nqm
              ind_qm = int_mb(i_qind+iqm-1)
              if((ind_qm.eq.ind1) .or. (ind_qm.eq.ind2)) then
                nl = nl+1
                int_mb(i_link+2*(nl-1))=ind_qm
                int_mb(i_link+2*nl-1)=ind_l
              end if
            end do
          end if 
        end do
      end do

      if(.not.ma_alloc_get(mt_int,2*nl,'link bond index',h_lb,i_lb))
     + call errquit(pname//'Failed to allocate heap',nl,
     &       MA_ERR)

      do i=1,nl
        int_mb(i_lb+i-1) = int_mb(i_link+2*(i-1))
        int_mb(i_lb+nl+i-1) = int_mb(i_link+2*i-1)
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname

      if(.not.ma_pop_stack(h_link))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_lind))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_qind))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)
      return
 
      end

      subroutine mm_count_qlinks(nps,h_ips,i_ips)
      implicit none
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "errquit.fh"

      integer nps
      integer h_ips
      integer i_ips
c     local variables:
      integer nl,i,j
      integer psb
      integer msb
      integer ind1,ind2
      integer ind_qm,ind_l
      integer iqm
      integer nqm
      integer h_qind,i_qind
      integer il,nl0
      integer h_lind,i_lind
      integer h_link,i_link
      logical duplicate
      character*32 pname
      pname = "mm_count_links"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_tot_nlink(nl0)
      if(qmmm_print_debug()) then
        write(*,*) pname,"number of link atoms",nl0
      end if
      if(nl0.eq.0) then
         nps = 0
         return
      end if

      call mm_get_tot_nqm(nqm)
      if(.not.ma_push_get(mt_int,nqm,'tmp qind',h_qind,i_qind))
     + call errquit(pname//'Failed to allocate stack',nqm,
     &       MA_ERR)

      call mm_get_solute_ind_gen(nqm,mm_quant,int_mb(i_qind))

      if(.not.ma_push_get(mt_int,nl0,'tmp lind',h_lind,i_lind))
     + call errquit(pname//'Failed to allocate stack',nl0,
     &       MA_ERR)

      call mm_get_solute_ind_gen(nl0,mm_link,int_mb(i_lind))

      if(.not.ma_push_get(mt_int,nl0*max_nb,'tmp links',h_link,i_link))
     + call errquit(pname//'Failed to allocate stack',nl,
     &       MA_ERR)
      call ifill(nqm,-1,int_mb(i_link),1)


      call cf_get_msb(msb)
      call cf_get_psb(psb)


      nl = 0
      do i=1,msb
        ind1=int_mb(psb+i-1)
        ind2=int_mb(psb+msb+i-1)
        do il=1,nl0
          ind_l = int_mb(i_lind+il-1)
          if((ind_l.eq.ind1) .or. (ind_l.eq.ind2)) then
            do iqm=1,nqm
              ind_qm = int_mb(i_qind+iqm-1)
              if((ind_qm.eq.ind1) .or. (ind_qm.eq.ind2)) then
                duplicate = .false.
                do j=1,nl
                  if(int_mb(i_link+j-1).eq.ind_qm) go to 10
                end do
                nl = nl+1
                int_mb(i_link+nl-1)=ind_qm
10              continue
              end if
            end do
          end if 
        end do
      end do
      
      nps = nl

      if(.not.ma_alloc_get(mt_int,nps,'pseudo atom index',h_ips,i_ips))
     + call errquit(pname//'Failed to allocate heap',nps,
     &       MA_ERR)

      do i=1,nps
        int_mb(i_ips+i-1) = int_mb(i_link+i-1)
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname

      if(.not.ma_pop_stack(h_link))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_lind))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_qind))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)

      return
 
      end

      subroutine mm_reset_solvent_bqzone(nt,
     >                       ai)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)

      integer i,j

c     mm local variables
      integer mwa
      integer mwm
      integer nwm
      integer pwz
      integer pwl
      character*32 pname
      pname = "mm_reset_solvent_bqzone"
      if(qmmm_print_debug()) write(*,*) "in "//pname


      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pwz(pwz)
      call mm_get_pwl(pwl)

 
      do i=1,nwm
        do j=1,nt/mwa
          if(int_mb(pwl+i-1).eq.ai(j)) then
             int_mb(pwz+i-1)=1
          end if
        end do
      end do


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_reset_solute_bqzone(nt,
     >                       ai)
      implicit none

#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "msgids.fh"
#include "mm_utils.fh"

      integer nt
      integer ai(nt)

      integer i,j

c     mm local variables
      integer msa
      integer nsa
      integer psl
      integer psz
      character*32 pname

      pname = "mm_reset_solute_bqzone"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_psz(psz)

      do i=1,nsa
        do j=1,nt
          if(int_mb(psl+i-1).eq.ai(j)) then
             int_mb(psz+i-1)=1
          end if
        end do
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

