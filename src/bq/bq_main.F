c $Id: bq_main.F,v 1.1 2005-04-01 03:10:53 marat Exp $

      function bq_init(rtdb)
      implicit none
#include "rtdb.fh"
#include "bq.fh"
#include "errquit.fh"

      integer rtdb
      logical bq_init
c
      integer handle
      character*255 namespace
      character*32 pname

      pname = "bq_init"
      
c
      bq_init = .true.
      if(.not. rtdb_cget(rtdb,"bq" , 1,namespace)) then
        bq_init = .false. 
        write(*,*) "no namespace"
        return
      end if
 
      if(.not.bq_create(namespace,handle)) then
         call errquit(pname//'failed to create bq',
     >                0,0)
      end if

      if(.not.bq_rtdb_load(rtdb,handle)) then
         call errquit(pname//'failed to load bq',
     >                0,RTDB_ERR)
      end if

      if(.not.bq_activate(handle)) then
         call errquit(pname//'failed to activate bq',
     >                0,RTDB_ERR)
      end if



      end

      function bq_add_nuc_rep_energy(geom,e)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"

      integer geom
      double precision e
      logical bq_add_nuc_rep_energy
c
      character*32 pname
      integer bq_handle
      integer bq_ncent
      integer ncent

      integer i,j
      integer i_cbq,i_qbq
      integer h_c,h_t,h_q
      integer i_c,i_t,i_q
      integer ioff,joff
      double precision rx,ry,rz,r
c
       logical  bq_ncenter
       external bq_ncenter

       logical  bq_index_coord
       external bq_index_coord

       logical  bq_index_charge
       external bq_index_charge

       logical  bq_destroy
       external bq_destroy

       logical  bq_get_active
       external bq_get_active



c   
      pname = "bq_add_nuc_rep_energy"
c
      bq_add_nuc_rep_energy = .true.
c
c     get bq information
c     ------------------
      if(.not.bq_get_active(bq_handle))
     >   call errquit(pname//':no active bq handle',0,0)
      if(.not.bq_ncenter(bq_handle,bq_ncent))
     >   call errquit(pname//':no bq centers',0,0)
      if(.not.bq_index_coord(bq_handle,i_cbq))
     >   call errquit(pname//':no bq coords',0,0)
      if(.not.bq_index_charge(bq_handle,i_qbq))
     >   call errquit(pname//':no bq coords',0,0)
c
      if(.not. geom_ncent(geom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(.not.ma_push_get(mt_dbl,3*ncent,'bq c',h_c,i_c))
     + call errquit(pname//': Failed to allocate memory for c',
     + 3*ncent, MA_ERR)
      if(.not.ma_push_get(mt_dbl,ncent,'bq q',h_q,i_q))
     + call errquit(pname//': Failed to allocate memory for q',ncent,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*ncent,'t',h_t,i_t))
     + call errquit(pname//': Failed to allocate memory for q',ncent,
     &       MA_ERR)

      if(.not.geom_cart_get(geom,ncent,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit(pname//':Failed to get geometry',0, GEOM_ERR)

      do i = 1,ncent
         ioff = i_c+(i-1)*3
         do j = 1, bq_ncent
           joff = i_cbq+(j-1)*3
           rx = dbl_mb(ioff)-dbl_mb(joff)
           rx = rx*rx
           ry = dbl_mb(ioff+1)-dbl_mb(joff+1)
           ry = ry*ry
           rz = dbl_mb(ioff+2)-dbl_mb(joff+2)
           rz = rz*rz
           r  = sqrt(rx+ry+rz)
           e = e + dbl_mb(i_q+i-1)*dbl_mb(i_qbq+j-1)/r
         end do
      end do
c
      if(.not.ma_pop_stack(h_t))
     & call errquit(pname//': Failed to deallocate stack t_all',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q))
     & call errquit(pname//': Failed to deallocate stack q_all',0,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit(pname//': Failed to deallocate stack c_all',0,
     &       MA_ERR)


      end 
