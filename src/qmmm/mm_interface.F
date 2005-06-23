c
c $Id: mm_interface.F,v 1.56 2005-06-23 17:32:02 marat Exp $
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
      if(qmmm_print_debug()) write(*,*) "0 "//pname
      call md_start()
      if(qmmm_print_debug()) write(*,*) "1 "//pname
  

c     this call is very important
c      call mm_task_sp()
      call md_sp()
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
      integer nf
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

      subroutine mm_end()
      implicit none
#include "md_common.fh"
#include "qmmm.fh"
      character*32 pname
c
      pname = "mm_end"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      call md_wrtrst(lfnrst,rfile,.false.)
c      call md_wrtrst(lfnqrs,filqrs,.false.)
      call md_finish()
c
      if(qmmm_print_debug()) write(*,*) "out "//pname
      return
      end

      subroutine mm_print_info()
      implicit none
#include "qmmm.fh"
#include "md_common.fh"
#include "mafdecls.fh"
c
      call prp_proper(0,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
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
      call md_sp()
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

      subroutine mm_get_tot_nbqs(nbq,link)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "qmmm.fh"
c
      integer nbq
      logical link
      character*32 pname
      pname = "mm_get_tot_nbqs"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      if(link) then
      call mm_get_solute_tot_na_gen(nbq,mm_link_bqzone)
      else
      call mm_get_solute_tot_na_gen(nbq,mm_bqzone)
      end if
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
      integer bq_exclude
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
      integer bq_exclude
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

      subroutine mm_get_solute_bq_ind(nt,link,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "msgids.fh"
      integer nt
      logical link
      integer ind(nt)
      character*32 pname
      pname = "mm_get_solute_bq_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      if(link) then
      call mm_get_solute_ind_gen(nt,mm_link_bqzone,ind)
      else
      call mm_get_solute_ind_gen(nt,mm_bqzone,ind)
      end if

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

      subroutine mm_get_solute_act_ind(nt,ind)
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
      pname = "mm_get_solute_act_ind"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_ind_gen(nt,mm_active,ind)
      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_act_ind(nt,ind)
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
      pname = "mm_get_solute_act_ind"
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

      
c      do i=1,nt
c       write(85,*) ind(i)
c      end do
c      call qmmm_sort(nt/mwa,ind)
c      do i=1,nt
c       write(86,*) ind(i)
c      end do
c      stop

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

      subroutine mm_get_solute_geom_quant(nt,
     >                       ai,
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
      character*16 t(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)
c
      integer i
      integer i_num,h_num
      character*32 pname
      pname = "mm_get_solute_geom_quant"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(.not.ma_push_get(mt_int,nt,'inum',h_num,i_num))
     + call errquit('qmmm: Failed to allocate memory for inum',nt,
     &       MA_ERR)
 
      call mm_get_solute_geom_gen(nt,
     >                       mm_quant_link,
     >                       ai,
     >                       int_mb(i_num),
     >                       t,
     >                       c,
     >                       q)

      do i=1,nt
         if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
         q(i)=dble(int_mb(i_num+i-1))
      end do

      if(.not.ma_pop_stack(h_num))
     & call errquit('qmmm: Failed to deallocate stack',nt,
     &       MA_ERR)
 
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

      subroutine mm_get_solute_geom_bq(nt,
     >                       link,
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
      logical link
      integer ai(nt)
      double precision c(3,nt)
      double precision q(nt)
c
      integer i_num,h_num
      integer i_t,h_t
      character*32 pname
      pname = "mm_get_solute_geom_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      if(.not.ma_push_get(mt_int,nt,'inum',h_num,i_num))
     + call errquit('qmmm: Failed to allocate memory for inum',nt,
     &       MA_ERR)
      if(.not.ma_push_get(MT_BYTE, 16*nt, 'bq tag array',
     &      h_t, i_t) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nt, MA_ERR)

      if(link) then
      call mm_get_solute_geom_gen(nt,
     >                       mm_link_bqzone,
     >                       ai,
     >                       int_mb(i_num),
     >                       byte_mb(i_t),
     >                       c,
     >                       q)

      else
      call mm_get_solute_geom_gen(nt,
     >                       mm_bqzone,
     >                       ai,
     >                       int_mb(i_num),
     >                       byte_mb(i_t),
     >                       c,
     >                       q)
      end if

      if(.not.ma_pop_stack(h_t))
     & call errquit('qmmm: Failed to deallocate stack h_t',nt,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_num))
     & call errquit('qmmm: Failed to deallocate stack h_num',nt,
     &       MA_ERR)

 
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
               inum(j) = 0
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

      subroutine mm_get_solute_cons_quant(nt,
     >                       ai,
     >                       act)
      implicit none

#include "mm_utils.fh"
#include "qmmm.fh"


      integer nt
      integer ai(nt)
      logical act(nt)
      character*32 pname
      pname = "mm_get_solute_cons_quant"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_cons_gen(nt,
     >                       mm_quant,
     >                       ai,
     >                       act)


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_cons_bq(nt,
     >                       ai,
     >                       act)
      implicit none

#include "mm_utils.fh"
#include "qmmm.fh"


      integer nt
      integer ai(nt)
      logical act(nt)
      character*32 pname
      pname = "mm_get_solute_cons_bq"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_cons_gen(nt,
     >                       mm_bqzone,
     >                       ai,
     >                       act)


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solute_cons_quant_link(nt,
     >                       ai,
     >                       act)
      implicit none

#include "mm_utils.fh"
#include "qmmm.fh"

      integer nt
      integer ai(nt)
      logical act(nt)
      character*32 pname
      pname = "mm_get_solute_cons_quant_link"
      if(qmmm_print_debug()) write(*,*) "in "//pname
      
      call mm_get_solute_cons_gen(nt,
     >                       mm_quant_link,
     >                       ai,
     >                       act)


      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

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

      integer i,j,k
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
         act(i)=int_mb(i_act+i-1).eq.0
      end do

      if(.not.ma_pop_stack(l_act))
     + call errquit('qmmm: Failed to deallocate memory for act',nt,
     &       MA_ERR)

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

      integer i,j,k

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

      integer i,j,k
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
                 if(mm_fixed(i)) then
                   int_mb(i_act+j-1)=1
                 end if
              end if
            end do
      end do

      call ga_igop(msg_qmmm_act,int_mb(i_act),nt,'+')


      do i=1,nt
         act(i) = int_mb(i_act+i-1).eq.0
      end do

      if(.not.ma_pop_stack(l_act))
     + call errquit('qmmm: Failed to deallocate memory for act',nt,
     &       MA_ERR)

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
      logical mm_type
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
         call cf_num2tag(inum(i),t(i))
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
        if(mm_bqwzone(i)) then
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
        end if
      end do


      call ga_dgop(msg_qmmm_cbqw,c,3*nt,'+')
      call ga_dgop(msg_qmmm_qbqw,q,nt,'+')

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
        if(mm_bqwzone(i)) then
            do j=1,nt/mwa
              if(int_mb(pwl+i-1).eq.ai(j)) then
                 do m=1,mwa
                   j1=(j-1)*mwa+m
                   do k=1,3
                     c(k,j1) = 
     >                dbl_mb(pxw+(i-1)+(k-1)*mwm+(m-1)*3*mwm)/cau2nm
                   end do
                 end do
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_cw,c,3*nt,'+')

      if(qmmm_print_debug()) write(*,*) "out "//pname
      end

      subroutine mm_get_solvent_cons_bq(nt,
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

      integer i,j,k,m,j1

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
        if(mm_bqwzone(i).and.mm_fixed_w(i)) then
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
         act(i)=int_mb(i_act+i-1).eq.0
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
c
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

      subroutine mm_get_solute_coord_bq(nt,
     >                       ai,
     >                       c)
      implicit none

#include "mm_utils.fh"
#include "qmmm.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)
      character*32 pname
      pname = "mm_get_solute_coord_bq"

      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_solute_coord_gen(nt,
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
           c(k)=dbl_mb(pxs+(i-1)+(k-1)*msa)*cnm2au
         end do
         call cf_num2tag(mm_element(i),t)
         write(23,FMT=9000)
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
      integer i,in
      integer psb
      integer msb
      integer ind1,ind2,ind0
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
      integer nl,i,j,in
      integer psb
      integer msb
      integer ind1,ind2,ind0
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

