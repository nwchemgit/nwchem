c
c $Id$
c
      subroutine qmmm_links_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_link_data.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"

      integer irtdb
c
      integer i
      character*32 pname
      character*255 basis_string
      character*255 basis_name
      character*255 link_ecp
      integer nb
       

      pname = "qmmm_init links: "

      if(qmmm_print_debug())
     >   write(*,*) "in",pname

      call qmmm_find_links()

c
c     set link atoms ecp if any
c     ------------------------
      if (.not.rtdb_cget(irtdb,"qmmm:link_ecp",1,link_ecp)) 
     >    call errquit('qmmm_geom: failed to get qmmm:link_ecp',
     &       0,RTDB_ERR)
c
      if(link_ecp.eq."auto") then
        if(qmmm_halogen_link()) then
          if(qmmm_master()) then
            write(6,'(//)')
          end if
          basis_name = "qmmm_zhang_6-31G*_ecp"
          nb = inp_strlen(basis_name)
          basis_string='ecp; F_L library "'//basis_name(1:nb)//'"; end'
          call nw_inp_from_character(irtdb, basis_string)
        end if
      end if


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_find_links()
      implicit none
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "mm_utils.fh"
#include "errquit.fh"
#include "qmmm_link_data.fh"

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
      pname = "qmmm_init_hlinks"
      if(qmmm_print_debug()) write(*,*) "in "//pname

      call mm_get_tot_nlink(nl0)
      if(qmmm_print_debug()) then
        write(*,*) pname,"number of link atoms",nl0
      end if
      if(nl0.eq.0)  then
        nlink=nl0
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
     + call errquit(pname//'Failed to allocate stack',nl0,
     &       MA_ERR)
      call ifill(nl0*max_nb,-1,int_mb(i_link),1)


      call cf_get_msb(msb)
      call cf_get_psb(psb)


      nlink = 0
      do i=1,msb
        ind1=int_mb(psb+i-1)
        ind2=int_mb(psb+msb+i-1)
        do il=1,nl0
          ind_l = int_mb(i_lind+il-1)
          if((ind_l.eq.ind1) .or. (ind_l.eq.ind2)) then
            do iqm=1,nqm
              ind_qm = int_mb(i_qind+iqm-1)
              if((ind_qm.eq.ind1) .or. (ind_qm.eq.ind2)) then
                nlink = nlink+1
                int_mb(i_link+2*(nlink-1))=ind_qm
                int_mb(i_link+2*nlink-1)=ind_l
              end if
            end do
          end if 
        end do
      end do

      if(.not.ma_alloc_get(mt_int,2*nlink,'link bond index',h_lb,i_lb))
     + call errquit(pname//'Failed to allocate heap',nlink,
     &       MA_ERR)

      do i=1,nlink
        int_mb(i_lb+i-1) = int_mb(i_link+2*(i-1))
        int_mb(i_lb+nlink+i-1) = int_mb(i_link+2*i-1)
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname

      if(.not.ma_pop_stack(h_link))
     + call errquit(pname//'Failed to pop stack',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_lind))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_qind))
     + call errquit(pname//'Failed to pop stack',nqm,
     &       MA_ERR)

      return
 
      end

      subroutine qmmm_link_data_release()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_link_data.fh"

      if(nlink.ne.0) then
       if(.not.ma_free_heap(h_lb))
     &   call errquit('qmmm h_lb: Failed to deallocate heap',nlink,
     &       MA_ERR)
      end if

       return
      end

      function qmmm_get_nlink()
      implicit none
#include "qmmm_link_data.fh"
#include "qmmm_params.fh"


      integer qmmm_get_nlink

      if(link_atom.eq.qmmm_hydrogen) then
      qmmm_get_nlink = nlink
      else
      qmmm_get_nlink = 0
      end if

      end 

      function qmmm_get_nps()
      implicit none
#include "qmmm_link_data.fh"

      integer qmmm_get_nps

      qmmm_get_nps = nlink

      end 

      function qmmm_get_link_atom()
      implicit none
#include "qmmm_link_data.fh"

      integer qmmm_get_link_atom

      qmmm_get_link_atom = link_atom

      end 

      function qmmm_h_link()
      implicit none
#include "qmmm_link_data.fh"
#include "qmmm_params.fh"

      logical qmmm_h_link

      qmmm_h_link=link_atom.eq.qmmm_hydrogen

      end 

      function qmmm_halogen_link()
      implicit none
#include "qmmm_link_data.fh"
#include "qmmm_params.fh"

      logical qmmm_halogen_link

      qmmm_halogen_link=link_atom.eq.qmmm_halogen

      end 

      function qmmm_get_i_lb()
      implicit none
#include "qmmm_link_data.fh"
      integer qmmm_get_i_lb

      qmmm_get_i_lb = i_lb

      end 

      function qmmm_get_link_atom_index()
      implicit none
#include "qmmm_link_data.fh"
      integer qmmm_get_link_atom_index

      qmmm_get_link_atom_index = i_lb+nlink

      end 

      function qmmm_get_i_ips()
      implicit none
#include "qmmm_link_data.fh"
      integer qmmm_get_i_ips

      qmmm_get_i_ips = i_lb

      end 

      subroutine qmmm_links_adjust(n, ai, t, c, q)
      implicit none
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_link_data.fh"

      integer n        
      integer ai(n)
      character*16 t(n)    
      double precision c(3, n) 
      double precision q(n) 
 
c     local variables:
      integer i,j,k
      integer nqm
      double precision a
      double precision rq,rl,rh
      integer icl,icq
      character*32 pname

      pname = "qmmm_links_adjust"

      if(qmmm_print_debug()) write(*,*) "in "//pname

      a = 0.709d0

      do i=1,nlink
        icl = 0
        icq = 0
        do j=1,n
          if(ai(j).eq.int_mb(i_lb+i-1)) icq = j
          if(ai(j).eq.int_mb(i_lb+nlink+i-1)) icl = j
        end do
        if(qmmm_h_link()) then
         if(icq.ne.0 .and. icl.ne.0) then
           do k=1,3
              rq  = c(k,icq)
              rl  = c(k,icl)
              rh = rq - (rq-rl)*a
              c(k,icl) = rh
           end do
           q(icl) = 1.0d0
           t(icl) = "H_L"
         else
           write(*,*) "no links found"
         end if
        end if
        if(qmmm_halogen_link()) then
         if( icq.ne.0) then
           q(icq) = 9.0d0
           t(icq) = "F_L"
         end if
        end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname

      return
 
      end

      subroutine qmmm_links_adjust_forces(n, ai, g)
      implicit none
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_link_data.fh"

      integer n        
      integer ai(n)
      double precision g(3, n) 
 
c     local variables:
      integer i,j,k
      integer nqm
      double precision a
      double precision gq,gl,gh
      integer il,iq
      character*32 pname

      if(qmmm_halogen_link()) return

      pname = "qmmm_links_adjust_forces"

      if(qmmm_print_debug()) write(*,*) "in "//pname

      a = 0.709d0

      do i=1,nlink
        il = 0
        iq = 0
        do j=1,n
          if(ai(j).eq.int_mb(i_lb+i-1)) iq = j
          if(ai(j).eq.int_mb(i_lb+nlink+i-1)) il = j
        end do
         if(iq.ne.0 .and. il.ne.0) then
           do k=1,3
              gq  = g(k,iq)
              gl  = g(k,il)
              g(k,il) = a*gl
              g(k,iq) = gq+(1.0d0-a)*gl
           end do
         else
           write(*,*) "no links found"
         end if
      end do

      if(qmmm_print_debug()) write(*,*) "out "//pname

      return
 
      end

