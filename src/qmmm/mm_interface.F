c
c $Id: mm_interface.F,v 1.2 2004-02-12 23:14:33 marat Exp $
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

      subroutine mm_get_tot_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nqm
      integer tmp(1)
c
      integer msa
      integer nsaloc
      integer psdt

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)

      call mm_get_nqm0(msa,nsaloc,
     &        int_mb(psdt),
     &        tmp(1))
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
      integer msa
      integer nsaloc
      integer psdt

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)

      call mm_get_nqm0(msa,nsaloc,
     &        int_mb(psdt),
     &        nqm)
      end

      subroutine mm_get_nqm0(msa,nsa,
     &                     isdt,nqm)
      implicit none
#include "util.fh"
#include "qmmm_common.fh"
c
      integer msa,nsa
      integer isdt(msa)
      integer nqm
c
      integer i
      integer nqatom,nlatom
c
      nqatom=0
      nlatom=0
      do i=1,nsa
        if(iand(isdt(i),lquant).ne.0) then
          nqatom=nqatom+1
        else if (and(isdt(i),lqlink).ne.0) then
          nlatom=nlatom+1
        endif
      end do
      nqm = nqatom + nlatom
      end

      subroutine mm_qm_mapping(nt,ai)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer nt
      integer ai(nt)
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

      call mm_qm_mapping0(nt,n,msa,nsaloc,
     &                   int_mb(i_ioff),
     &                   int_mb(psdt),
     &                   int_mb(psl),
     &                   ai)

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',n,
     &       MA_ERR)

      end

      subroutine mm_qm_mapping0(nt,np,msa,nsa,loff,isdt,isl,ai)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "qmmm_data.fh"
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

      subroutine mm_get_qm_geom(nt,ai,inum,c,q,m,t)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer nt
      integer ai(nt)
      integer inum(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)
      character*16 t(nt)
c     local variables
      integer msa
      integer nsaloc
      integer psdt
      integer psl
      integer psat
      integer pxs

      call mm_get_msa(msa)
      call mm_get_nsaloc(nsaloc)
      call mm_get_psdt(psdt)
      call mm_get_psat(psat)
      call mm_get_psl(psl)
      call mm_get_pxs(pxs)


      call mm_get_qm_geom0(nt,msa,nsaloc,
     &                   int_mb(psl),
     &                   int_mb(psdt),
     &                   int_mb(psat),
     &                   dbl_mb(pxs),
     &                   ai,
     &                   inum,
     &                   c,
     &                   q,
     &                   m,
     &                   t)

      end

      subroutine mm_get_qm_geom0(nt,msa,nsa,
     &                   isl,
     &                   isdt,
     &                   isat,
     &                   xs,
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
#include "qmmm_data.fh"
#include "msgids.fh"
#include "inp.fh"
      integer nt,msa,nsa
      integer isl(msa)
      integer isdt(msa)
      integer isat(msa)
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
              if(isl(i).eq.ai(j)) then
                 write(*,*) "found match",isl(i),ai(j)
                 write(*,*) "xs=",xs(i,1),xs(i,2),xs(i,3)
                 c(1,j)=xs(i,1)*cnm2au
                 c(2,j)=xs(i,2)*cnm2au
                 c(3,j)=xs(i,3)*cnm2au
                 atn=cf_element(isat(i))
                 if(islink) then
c                   if(atn.le.85) atn=85
c                   if(atn.le.53) atn=53
c                   if(atn.le.35) atn=35
c                   if(atn.le.17) atn=17
c                   if(atn.le.9)  atn=9
c                   if(atn.le.1)  atn=1
                    atn=1
                 end if
                 inum(j)=atn
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_c,c,3*nt,'+')
      call ga_igop(msg_qmmm_inum,inum,nt,'+')

      do i=1,nt
         q(i) = dble(inum(i))
         call cf_num2tag(inum(i),t(i))
         if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do

      end

