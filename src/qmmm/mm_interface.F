c
c $Id: mm_interface.F,v 1.4 2004-07-21 01:29:47 marat Exp $
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
c
      integer nqm
      integer tmp(1)
c
      call mm_get_nqm(tmp(1))
      call ga_igop(msg_qmmm_nat1,tmp,1,'+')
      nqm=tmp(1)
      end

      subroutine mm_get_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
c
      integer nqm
c
      integer i
      integer nqatom,nlatom
      integer msa
      integer nsaloc
      integer psdt
      integer lquant
      integer lqlink

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_lquant(lquant)
      call mm_get_lqlink(lqlink)
c
      nqatom=0
      nlatom=0
      do i=1,nsaloc
        if(iand(int_mb(psdt+i-1),lquant).ne.0) then
          nqatom=nqatom+1
        else if (and(int_mb(psdt+i-1),lqlink).ne.0) then
          nlatom=nlatom+1
        endif
      end do
      nqm = nqatom + nlatom
      end

      subroutine mm_get_nbqw(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
c
      integer nbq
c
      integer mwm
      integer mwa
      integer nwmloc
      integer pwz

      call mm_get_mwm(mwm)
      call mm_get_mwa(mwa)
      call mm_get_nwmloc(nwmloc)
      call mm_get_pwz(pwz)

      call mm_get_nbqw0(mwm,nwmloc,
     &        mwa,
     &        int_mb(pwz),
     &        nbq)
      end

      subroutine mm_get_nbqw0(mwm,nwm,mwa,
     &                     iwz,nbqw)
      implicit none
#include "util.fh"
#include "qmmm_common.fh"
c
      integer mwm,nwm,mwa
      integer iwz(mwm)
      integer nbqw
c
      integer i
c
      nbqw=0
      do i=1,nwm
        if(iwz(i).gt.0) nbqw=nbqw+mwa
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
c      call ga_igop(msg_qmmm_nat2,tmp,1,'+')
      nbq=tmp(1)
      end

      subroutine mm_get_nbqs(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_params.fh"
c
      integer nbq
c
      integer msa
      integer nsaloc
      integer psdt
      integer psz

      write(*,*) "mm_get_nbq"

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psz(psz)

      call mm_get_nbqs0(msa,nsaloc,
     &        int_mb(psdt),
     &        int_mb(psz),
     &        nbq)
      end

      subroutine mm_get_tot_nbqs(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nbq
      integer tmp(1)
c
      call mm_get_nbqs(tmp(1))
      call ga_igop(msg_qmmm_nat2,tmp,1,'+')
      nbq=tmp(1)
      end

      subroutine mm_get_nbqs0(msa,nsa,
     &                     isdt,isz,nbq)
      implicit none
#include "util.fh"
#include "qmmm_common.fh"
c
      integer msa,nsa
      integer isdt(msa)
      integer isz(msa)
      integer nbq
c
      integer i
c
      nbq=0
      do i=1,nsa
         write(*,*) "isdt(i),iand(isdt(i),lquant),isz(i)",
     >     isdt(i),iand(isdt(i),lquant),isz(i)
c        if(iand(isdt(i),lqsoft).ne.0) isz(i)=0
        if(iand(isdt(i),lquant).eq.0 .and. 
     >    iand(isdt(i),lqlink).eq.0) then
        if(isz(i).gt.0) nbq=nbq+1
        end if
      end do

      end

      subroutine mm_get_geom_bqs(nt,
     >                       c,
     >                       q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"

      integer nt
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
      integer psq
      integer pxs
      integer psdt

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsa)
      call mm_get_psz(psz)
      call mm_get_pxs(pxs)
      call mm_get_psq(psq)
      call mm_get_psdt(psdt)

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
c      call ga_igop(msg_qmmm_offset,loff,np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioq=int_mb(i_off+me-1)
      else
        ioq=0
      end if

      write(*,*) "before mm_get_geom_bqs0"
     
      call mm_get_geom_bqs0(nbqs,
     &                   msa,
     &                   nsa,
     &                   ioq,
     &                   int_mb(psdt),
     &                   int_mb(psz),
     &                   int_mb(psq),
     &                   dbl_mb(pxs),
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
     &                   isdt,
     &                   isz,
     &                   isq,
     &                   xs,
     &                   c,
     &                   q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"

      double precision cf_charge
      external cf_charge


      integer nt,msa,nsa
      integer ioff
      integer isdt(msa)
      integer isz(msa)
      integer isq(msa)
      double precision xs(msa,3)
      double precision c(3,nt)
      double precision q(nt)


c     local variables
      integer i
      integer ioq


      ioq = ioff
      do i=1,nsa
      if(isz(i).ne.0 .and. iand(isdt(i),lquant).eq.0 .and.
     >   iand(isdt(i),lqlink).eq.0) then
         ioq=ioq+1
         c(1,ioq)=xs(i,1)*cnm2au
         c(2,ioq)=xs(i,2)*cnm2au
         c(3,ioq)=xs(i,3)*cnm2au
         q(ioq)=cf_charge(isq(i))
      endif
      end do

      end

      subroutine mm_get_geom_bqw(nt,
     >                       c,
     >                       q)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
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

      call mm_get_mwa(mwa)
      call mm_get_mwm(mwm)
      call mm_get_nwmloc(nwm)
      call mm_get_pwz(pwz)
      call mm_get_pxw(pxw)
      
      write(*,*) "mwm,mwa,nwm",mwm,mwa,nwm

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
c      call ga_igop(msg_qmmm_offset,loff,np,'+')
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

      write(*,*) "after mm_get_geom_bqw0"
 
      write(*,*) " mm_get_geom_bqw0 1"
      
c      write(48,*) nbqw
c      write(48,*) "  "
c      do i=1,nbqw
c       write(*,*) "i=",i
c       write(48,'("Bqw  ",3(10X,F10.6),5X,"charge ",10X,F10.6)')
c     >           c(1,i)*0.529177249d00,
c     >            c(2,i)*0.529177249d00,
c     >            c(3,i)*0.529177249d00, q(i)
c
c      end do
c      write(*,*) " mm_get_geom_bqw0 1"

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
#include "qmmm_params.fh"

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
      if(iwz(i).gt.0) then
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

      subroutine mm_get_ind(nt,ind)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
      integer nt
      integer ind(nt)

c     local variables
      integer n
      integer l_ioff,i_ioff
      integer msa
      integer nsaloc
      integer psdt
      integer psl

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psl(psl)


      n=ga_nnodes()

      if(.not.ma_push_get(mt_int,n,'ioff',l_ioff,i_ioff))
     + call errquit('qmmm: Failed to allocate memory for ioff',n,
     &       MA_ERR)

      call mm_get_ind0(nt,n,msa,nsaloc,
     &                   int_mb(i_ioff),
     &                   int_mb(psdt),
     &                   int_mb(psl),
     &                   ind)

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',n,
     &       MA_ERR)

      end

      subroutine mm_get_ind0(nt,np,msa,nsa,loff,isdt,isl,ai)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "qmmm_params.fh"
      integer nt,np
      integer msa,nsa
      integer loff(np)
      integer ai(nt)
      integer isl(msa)
      integer isdt(msa)

      integer i
      integer me
      integer nqa
      integer ioff
      integer ioq

      do i=1,np
        loff(i)=0
      end do
c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
      call mm_get_nqm(nqa)
      do i=me+1,np
        loff(i)=nqa
      end do    
      call ga_igop(msg_qmmm_offset,loff,np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioff=loff(me)
      else
        ioff=0
      end if

      ioq=ioff
      do i=1,nsa
      if(iand(isdt(i),lquant).ne.0 .or.
     &   iand(isdt(i),lqlink).ne.0) then
         ioq=ioq+1
         ai(ioq)=isl(i)    
      endif
      end do
      do i=1,nt
        write(*,*) "i,ai",i,ai(i)
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

      do i=1,nsa
         islink=iand(isdt(i),lquant).ne.0
         j=cf_element(isat(i))
         call cf_num2tag(j,tag)
         write(*,*) i,tag,xs(i,1),xs(i,2),xs(i,3),isl(i),islink
      end do

      write(*,*) "quantum atoms"
      call dfill(3*nt,0.0d0,c,1)
      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        if(iand(isdt(i),lquant).ne.0 .or.
     &   islink) then
            do j=1,nt
              write(*,*) "j,nt",j,nt
              if(isl(i).eq.ai(j)) then
                 write(*,*) "found match",isl(i),ai(j)
                 write(*,*) "xs=",xs(i,1),xs(i,2),xs(i,3)
                 c(1,j)=xs(i,1)*cnm2au
                 c(2,j)=xs(i,2)*cnm2au
                 c(3,j)=xs(i,3)*cnm2au
                 atn=cf_element(isat(i))
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
                 q(j) = cf_charge(isq(i))
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

c      subroutine mm_create_bq_geom(irtdb)
c      implicit none
c
c#include "global.fh"
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "geom.fh"
c#include "qmmm.fh"
c
c      integer irtdb
cc     qm local variables
c      integer np
c      integer nbq
c      integer i_cbq,l_cbq
c      integer i_m,l_m
c      integer i_q,l_q
c      integer i_t,l_t
c      integer i_i,l_i
c      integer i_ioff,l_ioff
c      integer igeom
cc     mm local variables
c      integer msa
c      integer nsaloc
c      integer psdt
c      integer psat
c      integer pxs
c      integer psz
c      integer psq
c
cc
cc     get total number of bq atoms
c      call mm_get_tot_nbq(nbq)
c
cc     number of nodes
cc     --------------
c      np = ga_nnodes()
c
c      if(.not.ma_push_get(mt_dbl,3*nbq,'m',l_m,i_m))
c     + call errquit('qmmm: Failed to allocate for m',nbq, MA_ERR)
c      if(.not.ma_push_get(mt_dbl,nbq,'q',l_q,i_q))
c     + call errquit('qmmm: Failed to allocate for q',nbq, MA_ERR)
c      if(.not.ma_push_get(mt_byte,16*nbq,'t',l_t,i_t))
c     + call errquit('qmmm: Failed to allocate for t',nbq, MA_ERR)
c      if(.not.ma_push_get(mt_int,nbq,'i',l_i,i_i))
c     + call errquit('qmmm: Failed to allocate for i',nbq, MA_ERR)
c      if(.not.ma_push_get(mt_dbl,3*nbq,'cbq',l_cbq,i_cbq))
c     + call errquit('qmmm: Failed to allocate for m',3*nbq, MA_ERR)
c      if(.not.ma_push_get(mt_int,np,'ioff',l_ioff,i_ioff))
c     + call errquit('qmmm: Failed to allocate memory for ioff',np,
c     &       MA_ERR)
c
c 
c      call mm_get_msa(msa)
c      call mm_get_nsaloc(nsaloc)
c      call mm_get_psdt(psdt)
c      call mm_get_psat(psat)
c      call mm_get_pxs(pxs)
c      call mm_get_psz(psz)
c      call mm_get_psq(psq)
c
c
cc      call mm_get_bq_geom0(nt,np,msa,nsaloc,
cc     &                   int_mb(i_ioff),
cc     &                   int_mb(psdt),
cc     &                   int_mb(psat),
cc     &                   int_mb(psq),
cc     &                   int_mb(psz),
cc     &                   dbl_mb(pxs),
cc     &                   int_mb(i_i),
cc     &                   dbl_mb(i_c),
cc     &                   dbl_mb(i_q),
cc     &                   dbl_mb(i_m),
cc     &                   byte_mb(i_t))
cc
cc
cc       if(.not.geom_create(igeom,'geometry'))
cc     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
cc      if(.not.geom_cart_set(igeom,nt,byte_mb(i_t),
cc     + dbl_mb(i_c),dbl_mb(i_q)))
cc     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
cc      if(.not.geom_masses_set(igeom,nt,dbl_mb(i_m)))
cc     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
cc      call geom_compute_values(igeom)
cc
cc      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
cc     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)
cc
cc       if(.not.geom_destroy(igeom))
cc     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)
c
c
c      if(.not.ma_pop_stack(l_ioff))
c     & call errquit('qmmm: Failed to deallocate stack',nbq,
c     &       MA_ERR)
c      if(.not.ma_pop_stack(l_cbq))
c     & call errquit('qmmm: Failed to deallocate stack',nbq,
c     &       MA_ERR)
c      if(.not.ma_pop_stack(l_i))
c     & call errquit('qmmm: Failed to deallocate stack',nbq,
c     &       MA_ERR)
c      if(.not.ma_pop_stack(l_t))
c     & call errquit('qmmm: Failed to deallocate stack',nbq,
c     &       MA_ERR)
c      if(.not.ma_pop_stack(l_q))
c     & call errquit('qmmm: Failed to deallocate stack',nbq,
c     &       MA_ERR)
c      if(.not.ma_pop_stack(l_m))
c     & call errquit('qmmm: Failed to deallocate stack',nbq,
c     &       MA_ERR)
c
c
c
c      end
c
c
c      subroutine mm_get_bq_geom0(nt,np,msa,nsa,
c     &                   loff,
c     &                   isdt,
c     &                   isat,
c     &                   isq,
c     &                   isz,
c     &                   xs,
c     &                   inum,
c     &                   cbq,
c     &                   q,
c     &                   m,
c     &                   t)
c      implicit none
c#include "global.fh"
c#include "mafdecls.fh"
c#include "msgids.fh"
c#include "qmmm_params.fh"
c      integer nt,np
c      integer msa,nsa
c      integer loff(np)
c      integer isdt(msa)
c      integer isat(msa)
c      integer isz(msa)
c      integer xs(msa)
c      integer inum(nt)
c      double precision cbq(3,nt)
c      double precision q(nt)
c      double precision m(nt)
c      character*16 t(nt)
c
c      integer i
c      integer me
c      integer nbq
c      integer ioff
c      integer ioq
c
c      do i=1,np
c        loff(i)=0
c      end do
cc
cc     accumulate offsets
cc     ------------------
c      me = ga_nodeid()
c      call mm_get_nbq(nbq)
c      do i=me+1,np
c        loff(i)=nbq
c      end do    
c      call ga_igop(msg_qmmm_offset1,loff,np,'+')
cc
cc     set offset for a current node
cc     -----------------------------
c      if(me.gt.0) then
c        ioff=loff(me)
c      else
c        ioff=0
c      end if
c
c      ioq=ioff
c      do i=1,nsa
c      if(iand(isdt(i),lqsoft).ne.0) 
c         if(isz(i).gt.0) then
c         ioq=ioq+1
c         cbq(1,ioq)=xs(i,1)*cnm2au
c         cbq(2,ioq)=xs(i,2)*cnm2au
c         cbq(3,ioq)=xs(i,3)*cnm2au
c         inum(ioq)=cf_element(isat(i))
c         end if
c      endif
c      end do
c
c      call ga_dgop(msg_qmmm_cbq,cbq,3*nt,'+')
c      call ga_igop(msg_qmmm_inumbq,inum,nt,'+')
c
c      do i=1,nt
c         q(i) = cf_charge(isq(i))
c         t(i) = 'Bq' 
c         if(.not.geom_tag_to_default_mass(t(i),m(i)))
c     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
c      end do
c
c      end

