c
c $Id: mm_interface.F,v 1.18 2004-09-14 17:17:36 marat Exp $
c
      subroutine mm_init(irtdb0)
      implicit none
c
      integer irtdb0
c
      call timer_init()
      call mm_set_irtdb(irtdb0)
      call md_rdinp()
      call md_start()
      call mm_task_sp()
c
      return
      end

      subroutine mm_end()
      implicit none
#include "md_common.fh"
c
      call md_wrtrst(lfnqrs,filqrs,.false.)
      call md_finish()
c
      return
      end

      subroutine mm_task_sp()
      implicit none
c
      call md_sp()
c
      return
      end

      subroutine mm_get_tot_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
c
      integer nqm
c
      call mm_get_solute_tot_na_gen(nqm,mm_quant_link)
      end

      subroutine mm_get_tot_nlink(n)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
#include "mm_utils.fh"
c
      integer n
c
      call mm_get_solute_tot_na_gen(n,mm_link)
      end

      subroutine mm_get_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
c
      integer nqm
c

      call mm_get_solute_na_gen(nqm,mm_quant_link)

      end

      subroutine mm_get_solute_tot_na_gen(na,mm_type)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer na
      logical mm_type
      external mm_type
      integer tmp(1)
c
      call mm_get_solute_na_gen(tmp(1),mm_type)
      call ga_igop(msg_qmmm_nqm,tmp,1,'+')
      na=tmp(1)
      end

      subroutine mm_get_solute_na_gen(na,mm_type)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
c
      integer i
      logical mm_type
      external mm_type

      integer na
      integer nsaloc

      call mm_get_nsaloc(nsaloc)
c
      na=0
      do i=1,nsaloc
        if(mm_type(i)) then
          na=na+1
        endif
      end do
      end

      subroutine mm_get_tot_nbqw(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nbq
      integer tmp(1)
c
      call mm_get_nbqw(tmp(1))
      call ga_igop(msg_qmmm_nbqw,tmp,1,'+')
      nbq=tmp(1)
      end

      subroutine mm_get_nbqw(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
c
      integer nbq
c
      integer i
      integer mwa
      integer nwmloc

      call mm_get_mwa(mwa)
      call mm_get_nwmloc(nwmloc)

      nbq=0
      do i=1,nwmloc
        if(mm_bqwzone(i)) nbq=nbq+mwa
      end do

      end

      subroutine mm_get_tot_nbqs(nbq,link)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
c
      integer nbq
      logical link
      if(link) then
      call mm_get_solute_tot_na_gen(nbq,mm_link_bqzone)
      else
      call mm_get_solute_tot_na_gen(nbq,mm_bqzone)
      end if
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
      integer bq_exclude

      call mm_get_nsaloc(nsaloc)

      bq_exclude = qmmm_get_bq_exclude()
      do i=1,nsaloc
         if(mm_lqsoft(i)) then
           if((bq_exclude.eq.qmmm_hydrogen).and.(mm_element(i).eq.1))
     >     call mm_clear_bqzone(i)
           if(bq_exclude.eq.qmmm_all_atoms) call mm_clear_bqzone(i)
         endif
      end do

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
      
      call mm_get_solute_ind_gen(nt,mm_quant_link,ind)
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
      
      call mm_get_solute_ind_gen(nt,mm_link,ind)
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
      
      if(link) then
      call mm_get_solute_ind_gen(nt,mm_link_bqzone,ind)
      else
      call mm_get_solute_ind_gen(nt,mm_bqzone,ind)
      end if

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
      
      call mm_get_solute_ind_gen(nt,mm_quant,ind)
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
      integer pwl

      call ifill(nt,0,ind,1)

      call mm_get_nwmloc(nwm)

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
      call mm_get_nbqw(nqa)
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

      end

      subroutine mm_get_solvent_ind_bq(nt,ind)
      implicit none
#include "mm_utils.fh"
      integer nt
      integer ind(nt)

      call mm_get_solvent_ind_gen(nt,mm_bqwzone,ind)
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

      end

      subroutine mm_get_solvent_geom(nt,
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
                     dbl_mb(pxw+(i-1)+(k-1)*mwm+(m-1)*3*mwm)=c(k,j1)*
     >                                                       cau2nm
                   end do
                 end do
              end if
            end do
        end if
      end do

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


      call dfill(3*nt,0.0d0,g,1)

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
                     g(k,j1)=dbl_mb(ioff)*cau2nm/cau2kj
                   end do
                 end do
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_forcew,g,3*nt,'+')

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
      
      if(link) then
      call mm_get_solute_force_gen(nt,mm_quant_link,ai,g)
      else
      call mm_get_solute_force_gen(nt,mm_quant,ai,g)
      end if
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
      
      call mm_get_solute_force_gen(nt,mm_bqzone,ai,g)
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

      end
c
      subroutine mm_set_solute_coord_quant(nt,
     >                       link,
     >                       ai,
     >                       c)
      implicit none

#include "mm_utils.fh"

      integer nt
      logical link
      integer ai(nt)
      double precision c(3,nt)

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
      end

      subroutine mm_set_solute_coord_bq(nt,
     >                       ai,
     >                       c)
      implicit none

#include "mm_utils.fh"

      integer nt
      integer ai(nt)
      double precision c(3,nt)

      call mm_set_solute_coord_gen(nt,
     >                       mm_bqzone,
     >                       ai,
     >                       c)

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


      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      do i=1,nsa
        if(mm_type(i)) then
            do j=1,nt
              if(int_mb(psl+i-1).eq.ai(j)) then
                 do k=1,3
                   dbl_mb(pxs+(i-1)+(k-1)*msa)=c(k,j)*cau2nm
                 end do
              end if
            end do
        end if
      end do

      end

