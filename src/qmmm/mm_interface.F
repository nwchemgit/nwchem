c
c $Id: mm_interface.F,v 1.13 2004-09-03 20:07:56 marat Exp $
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
c
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
      integer i
      integer nqatom,nlatom
      integer nsaloc

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
      integer nqm
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

      subroutine mm_get_tot_nbqs(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "mm_utils.fh"
c
      integer nbq
      call mm_get_solute_tot_na_gen(nbq,mm_bqzone)
      end

      subroutine mm_get_nbqs(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
c
      integer nbq
c
      call mm_get_solute_na_gen(nbq,mm_bqzone)

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

      subroutine mm_get_geom_bqs(nt,
     >                       bqi,
     >                       c,
     >                       q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"

      integer nt
      integer bqi(*)
      double precision c(3,nt)
      double precision q(nt)


c     local variables
      integer i
      integer np,me
      integer nbqs
      integer l_off,i_off

      integer ioq
      integer nsa
      integer msa
      integer psz
      integer psl
      integer psq
      integer pxs
      integer psdt

      call dfill(3*nt,0,c,1)
      call dfill(nt,0,q,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psz(psz)
      call mm_get_pxs(pxs)
      call mm_get_psq(psq)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)

      np=ga_nnodes()

      if(.not.ma_push_get(mt_int,np,'off',l_off,i_off))
     + call errquit('qmmm: Failed to allocate memory for ioff',np,
     &       MA_ERR)

      do i=1,np
        int_mb(i_off+i-1)=0
      end do
c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
      call mm_get_nbqs(nbqs)
      do i=me+1,np
        int_mb(i_off+i-1)=nbqs
      end do    
      call ga_igop(msg_qmmm_nbqs_off,int_mb(i_off),np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioq=int_mb(i_off+me-1)
      else
        ioq=0
      end if

      write(*,*) "before mm_get_geom_bqs0"
     
      call mm_get_geom_bqs0(nt,
     &                   msa,
     &                   nsa,
     &                   ioq,
     &                   int_mb(psl),
     &                   int_mb(psdt),
     &                   int_mb(psz),
     &                   int_mb(psq),
     &                   dbl_mb(pxs),
     &                   bqi,
     &                   c,
     &                   q)

      write(*,*) "after mm_get_geom_bqs0"
 
      if(.not.ma_pop_stack(l_off))
     & call errquit('qmmm: Failed to deallocate stack',nt,
     &       MA_ERR)
      
      write(43,*) nbqs
      write(43,*) "  "
      do i=1,nbqs
       write(43,'("Bq  ",3(10X,F10.6),5X,"charge ",10X,F10.6)')
     >           c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do

      end

      subroutine mm_get_geom_bqs0(nt,
     &                   msa,
     &                   nsa,
     &                   ioff,
     &                   isl,
     &                   isdt,
     &                   isz,
     &                   isq,
     &                   xs,
     &                   bqi,
     &                   c,
     &                   q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
#include "msgids.fh"

      double precision cf_charge
      external cf_charge


      integer nt,msa,nsa
      integer ioff
      integer isl(msa)
      integer isdt(msa)
      integer isz(msa)
      integer isq(msa)
      double precision xs(msa,3)
      integer bqi(nt)
      double precision c(3,nt)
      double precision q(nt)


c     local variables
      integer i
      integer ioq



      ioq = ioff
      do i=1,nsa
      if(mm_bqzone(i)) then
         ioq=ioq+1
         c(1,ioq)=xs(i,1)*cnm2au
         c(2,ioq)=xs(i,2)*cnm2au
         c(3,ioq)=xs(i,3)*cnm2au
         q(ioq)=cf_charge(isq(i))
         bqi(ioq)=isl(i)
      endif
      end do

      call ga_dgop(msg_qmmm_cbqs,c,3*nt,'+')
      call ga_dgop(msg_qmmm_qbqs,q,nt,'+')

      end

      subroutine mm_get_geom_bqw(nt,
     >                       c,
     >                       q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "msgids.fh"
#include "qmmm.fh"

      integer nt
      double precision c(3,nt)
      double precision q(nt)


c     local variables
      integer i
      integer np,me
      integer nbqw
      integer l_off,i_off

      integer ioq
      integer mwa
      integer mwm
      integer nwm
      integer pwz
      integer pxw

      call dfill(3*nt,0,c,1)
      call dfill(nt,0,q,1)

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pwz(pwz)
      call mm_get_pxw(pxw)
      
      np=ga_nnodes()

      if(.not.ma_push_get(mt_int,np,'off',l_off,i_off))
     + call errquit('qmmm: Failed to allocate memory for ioff',np,
     &       MA_ERR)

      do i=1,np
        int_mb(i_off+i-1)=0
      end do
c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
      call mm_get_nbqw(nbqw)
      do i=me+1,np
        int_mb(i_off+i-1)=nbqw
      end do    
      call ga_igop(msg_qmmm_nbqw_off,int_mb(i_off),np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioq=int_mb(i_off+me-1)
      else
        ioq=0
      end if

      write(*,*) "before mm_get_geom_bqw0"
     
      call mm_get_geom_bqw0(nt,
     &                   mwm,
     &                   mwa,
     &                   nwm,
     &                   ioq,
     &                   int_mb(pwz),
     &                   dbl_mb(pxw),
     &                   c,
     &                   q)

      if(.not.ma_pop_stack(l_off))
     & call errquit('qmmm: Failed to deallocate stack',nt,
     &       MA_ERR)


      end

      subroutine mm_get_geom_bqw0(nt,
     &                   mwm,
     &                   mwa,
     &                   nwm,
     &                   ioff,
     &                   iwz,
     &                   xw,
     &                   c,
     &                   q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
#include "qmmm_params.fh"
#include "msgids.fh"

      double precision cf_wcharge
      external cf_wcharge


      integer nt,mwa,nwm,mwm
      integer ioff
      integer iwz(mwm)
      double precision xw(mwm,3,mwa)
      double precision c(3,nt)
      double precision q(nt)


c     local variables
      integer i,j
      integer ioq


      ioq = ioff
      do i=1,nwm
      if(mm_bqwzone(i)) then
         do j=1,mwa
         ioq=ioq+1
         c(1,ioq)=xw(i,1,j)*cnm2au
         c(2,ioq)=xw(i,2,j)*cnm2au
         c(3,ioq)=xw(i,3,j)*cnm2au
         q(ioq)=cf_wcharge(j)
          write(*,'(I5,3(10X,F10.6),5X,"charge ",10X,F10.6)')
     >     ioq,
     >     c(1,i),
     >      c(2,i),
     >      c(3,i), q(i)

         end do
      endif
      end do

      call ga_dgop(msg_qmmm_cbqw,c,3*nt,'+')
      call ga_dgop(msg_qmmm_qbqw,q,nt,'+')

      write(48,*) nt
      write(48,*) "  "
      do i=1,nt
       write(*,*) "i=",i
       write(48,'("Bqw  ",3(10X,F10.6),5X,"charge ",10X,F10.6)')
     >           c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do

      call util_flush(48)

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
      
      call mm_get_solute_ind_gen(nt,mm_bqzone,ind)
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
      integer i,j
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
      call mm_get_nqm(nqa)
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
      write(*,*) "me,ioff",me,ioff

      ioq=ioff
      do i=1,nsa
      if(mm_type(i)) then
         ioq=ioq+1
         ind(ioq)=int_mb(psl+i-1)
      endif
      end do

      call ga_igop(msg_qmmm_ind,ind,nt,'+')

      do i=0,np-1
        if(i.eq.me) then
        write(*,*) "qmmm index array on proc",i
        do j=1,nt
        write(*,*) ind(j)
        end do
        end if
      end do

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',np,
     &       MA_ERR)

      end

      subroutine mm_get_solute_geom_gen(nt,
     >                       mm_type,
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
      logical mm_type
      integer ai(nt)
      integer inum(nt)
      character*16 t(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)
      external mm_type

      integer i,j,k

c     mm local variables
      integer msa
      integer nsa
      integer psdt
      integer psq
      integer psl
      integer psat
      integer pxs

      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)
      call ifill(nt,0,inum,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psdt(psdt)
      call mm_get_psat(psat)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)
      call mm_get_psq(psq)

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
         if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do

      end

      subroutine mm_get_geom(nt,
     >                       ai,
     >                       inum,
     >                       t,
     >                       c,
     >                       q,
     >                       m)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"

      integer nt
      integer ai(nt)
      integer inum(nt)
      character*16 t(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)

      integer i

c     mm local variables
      integer msa
      integer nsaloc
      integer psdt
      integer psq
      integer psl
      integer psat
      integer pxs

      call dfill(3*nt,0.0d0,c,1)
      call dfill(nt,0.0d0,q,1)
      call dfill(nt,0.0d0,m,1)
      call ifill(nt,0.0d0,inum,1)

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psat(psat)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)
      call mm_get_psq(psq)


      call mm_get_geom0(nt,msa,nsaloc,
     &                   int_mb(psl),
     &                   int_mb(psdt),
     &                   int_mb(psat),
     &                   dbl_mb(pxs),
     &                   int_mb(psq),
     &                   ai,
     &                   inum,
     &                   c,
     &                   q,
     &                   m,
     &                   t)

      write(42,*) nt
      write(42,*) "qm geometry"
      do i=1,nt
       write(42,'(A4,3(10X,F10.6),5X,"charge ",10X,F10.6)')
     >           t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do

      do i=1,nt
       write(45,FMT=9000)
     >           i,t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(45)

      end

      subroutine mm_get_geom0(nt,msa,nsa,
     &                   isl,
     &                   isdt,
     &                   isat,
     &                   xs,
     &                   isq,
     &                   ai,
     &                   inum,
     &                   c,
     &                   q,
     &                   m,
     &                   t)
      implicit none
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
#include "inp.fh"
#include "mm_utils.fh"
      integer nt,msa,nsa
      integer isl(msa)
      integer isdt(msa)
      integer isat(msa)
      integer   isq(msa)
      double precision xs(msa,3)
      integer ai(nt)
      integer inum(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)
      character*16 t(nt)
c
      integer cf_element
      external cf_element
      double precision cf_charge
      external cf_charge

c     local variables
      logical islink
      integer atn
      integer i,j
      character*16 tag
      integer i_ndx,h_ndx
      double precision r(3)
      

      if(.not.ma_push_get(mt_int,nt,'ndx',h_ndx,i_ndx))
     + call errquit('qmmm: Failed to allocate memory for ndx',
     + 3*nt, MA_ERR)
      

      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        if(iand(isdt(i),lquant).ne.0 .or.
     &   islink) then
            do j=1,nt
              if(isl(i).eq.ai(j)) then
                 call mm_get_solute_coord(i,r)
                 write(*,*) "found match",isl(i),ai(j)
                 write(*,*) "xs=",xs(i,1),xs(i,2),xs(i,3)
                 write(*,*) "r=",r(1),r(2),r(3)
                 c(1,j)=xs(i,1)*cnm2au
                 c(2,j)=xs(i,2)*cnm2au
                 c(3,j)=xs(i,3)*cnm2au
                 atn=mm_element(i)
                 int_mb(i_ndx+j-1)=isdt(i)
                 if(islink) then
                   atn = 9
c                   if(atn.le.85) atn=85
c                   if(atn.le.53) atn=53
c                   if(atn.le.35) atn=35
c                   if(atn.le.17) atn=17
c                   if(atn.le.9)  atn=9
c                   if(atn.le.1)  atn=1
                   call cf_num2tag(atn,tag)
                   write(*,*) "found link atom",atn,tag
                 end if
                 inum(j)=atn
c                 q(j) = cf_charge(isq(i))
                 q(j) = mm_solute_charge(i)
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_c,c,3*nt,'+')
      call ga_dgop(msg_qmmm_chg,q,nt,'+')
      call ga_igop(msg_qmmm_inum,inum,nt,'+')
      call ga_igop(msg_qmmm_ndx,int_mb(i_ndx),nt,'+')

      do i=1,nt
         call cf_num2tag(inum(i),t(i))
         if(iand(int_mb(i_ndx+i-1),lqlink).ne.0)
     &      t(i)=t(i)(1:index(t(i),' ')-1)//'_L'
         if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do

      if(.not.ma_pop_stack(h_ndx))
     + call errquit('qmmm: Failed to deallocate memory for ndx',
     + 3*nt, MA_ERR)

      end

      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isl,isat,isdt,isq,isz,xs,fs,energy,lnkatm,nbq,lesp)
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      external task_gradient,esp
      logical task_gradient,esp
c
      integer cf_element
      external cf_element
c
      integer inrtdb
      logical lesp
      integer mwm,nwm,mwa,nwa,msa,nsa,lnkatm,nbq
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa),isl(msa)
      real*8 energy,xw(mwm,3,mwa),fw(mwm,3,mwa),xs(msa,3),fs(msa,3)
c
      integer i
      integer j
      integer nt
      integer i_ai
      integer ai
      logical islink
      logical isquant


      nt   = qmmm_get_nqm()
      i_ai = qmmm_get_i_ai()
c      i_g = qmmm_get_i_g()

      write(*,*) "in qmmm force"


      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        isquant=iand(isdt(i),lquant).ne.0
        if(isquant .or.  islink) then
            do j=1,nt
              ai = int_mb(i_ai+j-1)
              if(isl(i).eq.ai) then
                 write(*,*) "found match in qmmm force",isl(i),ai
                 if(isquant) then
                 write(*,*) "zeroing qmmm forces"
                 fs(i,1)=0.0d0
                 fs(i,2)=0.0d0
                 fs(i,3)=0.0d0
c                 else
c                 fs(i,1)=fs(i,1)-g(1,j)*(cau2kj/cau2nm) 
c                 fs(i,2)=fs(i,2)-g(2,j)*(cau2kj/cau2nm)
c                 fs(i,3)=fs(i,3)-g(3,j)*(cau2kj/cau2nm)
                 end if
              end if
            end do
        end if
      end do

      end
c
      subroutine mm_get_force(nt,ind,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
      integer nt
      integer ind(nt)
      double precision  g(3,nt)
c     local variables
      integer msa
      integer nsaloc
      integer psdt
      integer psl
      integer pfs

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      call mm_get_force0(nt,msa,nsaloc,
     &                   int_mb(psl),
     &                   int_mb(psdt),
     &                   dbl_mb(pfs),
     &                   ind,
     &                   g)

      end

      subroutine mm_get_force0(nt,msa,nsa,
     &                   isl,
     &                   isdt,
     &                   fs,
     &                   ai,
     &                   g)
      implicit none
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
#include "inp.fh"
      integer nt,msa,nsa
      integer isl(msa)
      integer isdt(msa)
      double precision fs(msa,3)
      integer ai(nt)
      double precision g(3,nt)
c     local variables
      logical islink
      integer i,j


      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        if(iand(isdt(i),lquant).ne.0 .or.
     &   islink) then
            do j=1,nt
              if(isl(i).eq.ai(j)) then
                 g(1,j)=-fs(i,1)
                 g(2,j)=-fs(i,2)
                 g(3,j)=-fs(i,3)
              end if
            end do
        end if
      end do

      end

      subroutine mm_set_geom(nt,ind,c)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer nt
      integer ind(nt)
      double precision c(3,nt)

c     local variables
      integer msa
      integer nsaloc
      integer psdt
      integer psl
      integer pxs

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)

      call mm_set_geom0(nt,msa,nsaloc,
     &                   int_mb(psl),
     &                   int_mb(psdt),
     &                   dbl_mb(pxs),
     &                   ind,
     &                   c)

      end

      subroutine mm_set_geom0(nt,msa,nsa,
     &                   isl,
     &                   isdt,
     &                   xs,
     &                   ai,
     &                   c)
      implicit none
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
#include "inp.fh"
      integer nt,msa,nsa
      integer isl(msa)
      integer isdt(msa)
      double precision xs(msa,3)
      integer ai(nt)
      double precision c(3,nt)
c     local variables
      logical islink
      integer i,j

      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        if(iand(isdt(i),lquant).ne.0 .or.
     &   islink) then
            do j=1,nt
              if(isl(i).eq.ai(j)) then
                 xs(i,1)=c(1,j)*cau2nm
                 xs(i,2)=c(2,j)*cau2nm
                 xs(i,3)=c(3,j)*cau2nm
              end if
            end do
        end if
      end do

      end

      subroutine mm_update_force(nt,ind,g)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "rtdb.fh"
      integer nt
      integer ind(nt)
      double precision  g(3,nt)
c     local variables
      integer msa
      integer nsaloc
      integer psdt
      integer psl
      integer pfs

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)
      call mm_get_pfs(pfs)

      call mm_update_force0(nt,msa,nsaloc,
     &                   int_mb(psl),
     &                   int_mb(psdt),
     &                   dbl_mb(pfs),
     &                   ind,
     &                   g)

      end

      subroutine mm_update_force0(nt,msa,nsa,
     &                   isl,
     &                   isdt,
     &                   fs,
     &                   ai,
     &                   g)
      implicit none
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
#include "inp.fh"
      integer nt,msa,nsa
      integer isl(msa)
      integer isdt(msa)
      double precision fs(msa,3)
      integer ai(nt)
      double precision g(3,nt)
c     local variables
      logical islink
      integer i,j


      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        if(iand(isdt(i),lquant).ne.0 .or.
     &   islink) then
            do j=1,nt
              if(isl(i).eq.ai(j)) then
                 fs(i,1)=fs(i,1)-g(1,j)
                 fs(i,2)=fs(i,2)-g(2,j)
                 fs(i,3)=fs(i,3)-g(3,j)
              end if
            end do
        end if
      end do

      end


