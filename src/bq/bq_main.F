#define NO_BQGEM 1
c $Id$
C> \ingroup bq
C> @{
C>
C> \brief Initialize the Bq instances
C>
C> Initialize the Bq instances with data from the RTDB. 
C> The routine looks for the "bq" tag on the RTDB and loads
C> the corresponding string as the name of a Bq instance. 
C> Subsequently a Bq instance is created, the data of the
C> named instance loaded and the instance activated for
C> embedding calculations.
C>
      function bq_init(rtdb)
      implicit none
#include "rtdb.fh"
#include "bq.fh"
#include "errquit.fh"
#include "global.fh"

      integer rtdb    !< [Input] The RTDB handle
      logical bq_init
c
      integer handle
      character*255 namespace
      character*32 pname

      pname = "bq_init"
      
c
c      write(*,*) "in "//pname
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
C>
C> \brief Clean up all Bq instances
C>
      function bq_end()
      implicit none
#include "rtdb.fh"
#include "bq.fh"
#include "errquit.fh"
      logical bq_end
      
      character*30 pname
      pname ="bq_end"
c      write(*,*) "in "//pname
      bq_end = bq_destroy_all()

      end
C>
C> \brief Calculate the nuclear the interaction between the Bq and
C> the nuclear charges
C>
C> \return Return .true. if successfull and .false. otherwise.
C>
      function bq_nuc_energy(rtdb,ebq)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"
#include "bq_params.fh"
#include "stdio.fh"
c
      integer rtdb         !< [Input] The RTDB handle
      double precision ebq !< [Output] The nuclear-bq interaction energy
      logical bq_nuc_energy
c
      character*32 pname
      integer bq_handle
      integer bq_ncent
      integer ncent
c
      integer i,j
      integer i_cbq,i_qbq
      integer h_c,h_t,h_q
      integer i_c,i_t,i_q
      integer ioff,joff
      double precision rx,ry,rz,r
      integer geom
      double precision qiqj
c
      logical  bq_ncenter
      external bq_ncenter
c
      logical  bq_index_coord
      external bq_index_coord
c
      logical  bq_index_charge
      external bq_index_charge
c
      logical  bq_destroy
      external bq_destroy
c
      logical  bq_get_active
      external bq_get_active
c   
      pname = "bq_nuc_rep_energy"
c
      bq_nuc_energy = .true.
c
c     get bq information
c     ------------------
      if(.not.bq_get_active(bq_handle)) then
         ebq = 0.0d0
         return
      end if
      if(.not.bq_ncenter(bq_handle,bq_ncent))
     >   call errquit(pname//':no bq centers',0,0)
      if(.not.bq_index_coord(bq_handle,i_cbq))
     >   call errquit(pname//':no bq coords',0,0)
      if(.not.bq_index_charge(bq_handle,i_qbq))
     >   call errquit(pname//':no bq coords',0,0)
c
      if(.not.geom_create(geom,'geometry'))
     + call errquit(pname//'Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(rtdb,geom,'geometry'))
     + call errquit(pname//'Failed to create geometry',0, GEOM_ERR)
c
      if(.not. geom_ncent(geom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)
c
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

      ebq = 0.0d0
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
           qiqj = dbl_mb(i_q+i-1)*dbl_mb(i_qbq+j-1)
           if (r.gt.bq_smalldist) then
             ebq = ebq + qiqj/r
           end if
         end do
      end do
c
      if(ga_nodeid().eq.0) then
        write(luout,*) "Bq nuclear interaction energy =", ebq 
      end if
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

       if(.not.geom_destroy(geom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)

      end 
C>
C> \brief Adds the Bq-nuclear interaction energy to the energy argument
C> passed
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      function bq_add_nuc_rep_energy(geom,e)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"
#include "bq_params.fh"
#include "stdio.fh"
c
      integer geom       !< [Input] The geometry handle for the nuclear
                         !< coordinates
      double precision e !< [In/Output] The energy
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
#ifdef NO_BQGEM
#include "inp.fh"
      character*2 bqchar
#endif
      double precision rx,ry,rz,r
      double precision ebq
      double precision qiqj
c
      logical  bq_ncenter
      external bq_ncenter
c
      logical  bq_index_coord
      external bq_index_coord
c
      logical  bq_index_charge
      external bq_index_charge
c
      logical  bq_destroy
      external bq_destroy
c
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

      ebq = 0.0d0
      do i = 1,ncent
         ioff = i_c+(i-1)*3
#ifdef NO_BQGEM
         bqchar(1:1)=byte_mb(i_t+16*(i-1))
         bqchar(2:2)=byte_mb(i_t+16*(i-1)+1)
         if (.not.inp_compare(.false.,'bq', bqchar) ) then
#endif
         do j = 1, bq_ncent
           joff = i_cbq+(j-1)*3
           rx = dbl_mb(ioff)-dbl_mb(joff)
           rx = rx*rx
           ry = dbl_mb(ioff+1)-dbl_mb(joff+1)
           ry = ry*ry
           rz = dbl_mb(ioff+2)-dbl_mb(joff+2)
           rz = rz*rz
           r  = sqrt(rx+ry+rz)
           qiqj = dbl_mb(i_q+i-1)*dbl_mb(i_qbq+j-1)
           if (r.gt.bq_smalldist) then
             ebq = ebq + qiqj/r
           end if
         end do
#ifdef NO_BQGEM
      endif
#endif
      end do
c
      if(ga_nodeid().eq.0) then
        write(luout,*) "Nuclear repulsion energy =", e 
        write(luout,*) "Bq nuclear interaction energy =", ebq 
      end if
c 
      e = e + ebq
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
C> @}
