
      subroutine mm_geom_init(rtdb,geomname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"
      
      integer rtdb
      character*(*) geomname

      character*255 aname
      character*32 pname
      character*30 operation
      logical orestrt
      integer geom
      integer i, j      
      logical ignore      
      integer i_c,h_c
      integer i_tag,h_tag
      integer i_m,h_m
      integer i_atn,h_atn
      integer i_dbq,h_dbq
      double precision scale
      logical geom_tag_to_atn
      logical geom_tag_to_charge
      external geom_tag_to_atn
      external geom_tag_to_charge
      

      pname = "mm_geom_init"     
c     write(*,*) pname 
       
      ignore = rtdb_cget(rtdb,'task:operation',1,operation)

c     deallocate all previous allocated arrays just in case
      call mm_geom_end()

      nact = nqm
      if(qmlink) nact = nqm + nlink

      nfg = nqm + nlink

      aux_geom = .false.
c     TP: refer to qmmm.F for operation = "hessian" and
c         operation = "freq". aux_geom = .true.
      if(operation.eq."optimize")  aux_geom =.true.
      if(operation.eq."neb")  aux_geom =.true.
      if(operation.eq."freq")  aux_geom =.true.
      if(operation.eq."hessian")  aux_geom =.true.
      if(operation.eq."saddle")  aux_geom =.true.

c     initialize indexing for a full qm geometry
      if(.not.ma_alloc_get(mt_int,nact,'mm act atom ind',
     &                     h_iact,i_iact))
     & call errquit(pname//'Failed to allocate heap',nact,
     &       MA_ERR)

      if(.not.ma_alloc_get(mt_int,nfg,'mm fullg ind',h_ifg,i_ifg))
     & call errquit(pname//'Failed to allocate heap',nfg,
     &       MA_ERR)

      call icopy(nqm,int_mb(i_iqm),1,int_mb(i_ifg),1)
      call icopy(nlink,int_mb(i_lb+nlink),1,int_mb(i_ifg+nqm),1)
       
      do i=1,nact
        int_mb(i_iact+i-1) = int_mb(i_iqml+i-1)
      end do

      if(operation.eq."neb") then
        if(.not.rtdb_get(rtdb,"mm:neb:restart",mt_log,1,orestrt))
     >     orestrt = .false.
        if(orestrt) return
      end if


      if(.not.ma_push_get(mt_dbl,3*nact,'c',h_c,i_c))
     & call errquit('mm: Failed to allocate memory for c',
     & 3*nact, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nact,'q',h_dbq,i_dbq))
     & call errquit('mm: Failed to allocate memory for q',nact,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nact,'m',h_m,i_m))
     & call errquit('mm: Failed to allocate memory for m',nact,
     &       MA_ERR)
      if(.not.ma_push_get(mt_int,nact,'inum',h_atn,i_atn))
     & call errquit('mm: Failed to allocate memory for atn',nact,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*nact,'t',h_tag,i_tag))
     & call errquit('mm: Failed to allocate memory for t',nact,
     &       MA_ERR)

c     assign coordinates for active atoms
      do i=1,nact
        dbl_mb(i_c+(i-1)*3)   = dbl_mb(i_rqml+3*(i-1))
        dbl_mb(i_c+(i-1)*3+1) = dbl_mb(i_rqml+3*(i-1)+1)
        dbl_mb(i_c+(i-1)*3+2) = dbl_mb(i_rqml+3*(i-1)+2)
        call mm_coords_tag_set(byte_mb(i_tqml+16*(i-1)),
     &                         byte_mb(i_tag+16*(i-1)))
      end do
      
      call util_convert_units("angstrom","au",scale)
      call dscal(3*nact, scale,dbl_mb(i_c),1)

      if(.not.geom_tag_to_charge(nact,byte_mb(i_tag),dbl_mb(i_dbq)))
     & call errquit('mm: Failed to get charge from tag',nact,MA_ERR)
      
      if(.not.geom_tag_to_atn(nact,byte_mb(i_tag),int_mb(i_atn)))
     & call errquit('mm: Failed to get atn from tag',nact,MA_ERR)

c     if nlink > 0 assign empirical charge to mmlink
      if(qmlink) then
        do i=1,nlink
          dbl_mb(i_dbq+nqm+i-1) = dbl_mb(i_lnkchg+i-1) 
        end do
      end if
 
      call mm_atn_get_mass(nact,int_mb(i_atn),dbl_mb(i_m))

c     ignore = rtdb_delete(rtdb,"geometry")
      ignore = rtdb_delete(rtdb,geomname(1:inp_strlen(geomname)))
 
c     if(.not.geom_create(geom,"geometry"))
      if(.not.geom_create(geom,geomname(1:inp_strlen(geomname))))
     & call errquit('mm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(geom,nact,byte_mb(i_tag),
     &                                dbl_mb(i_c),
     &                                dbl_mb(i_dbq)))
     & call errquit('mm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(geom,nact,dbl_mb(i_m)))
     & call errquit('mm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(geom)

      if(.not.geom_rtdb_store(rtdb,geom,
     &        geomname(1:inp_strlen(geomname))))
     & call errquit('mm: Failed to store geom to rtdb',0, RTDB_ERR)

      if(.not.geom_destroy(geom))
     & call errquit('mm: Failed to destroy geometry',0, GEOM_ERR)

      if(.not.ma_pop_stack(h_tag))
     & call errquit('mm: Failed to deallocate stack t_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_atn))
     & call errquit('mm: Failed to deallocate stack atn_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m))
     & call errquit('mm: Failed to deallocate stack m_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_dbq))
     & call errquit('mm: Failed to deallocate stack q_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit('mm: Failed to deallocate stack c_all',nact,
     &       MA_ERR)

      end 

      subroutine mm_geom_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_geom_data.fh"
      
      character*30 pname
      pname = "mm_geom_end"

      if(nfg.gt.0) then
        if (.not.ma_free_heap(h_ifg))    goto 911
        if (.not.ma_free_heap(h_iact))   goto 911
        nfg = 0
        nact = 0
        aux_geom = .false.
      end if

      return
911   call errquit("error "//trim(pname),0,-1)

      end

      subroutine mm_geom_create_full(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_geom_data.fh"

      character*32 pname
      integer rtdb
      integer i_c,h_c
      integer i_t,h_t
      integer i_m,h_m
      integer i_atn,h_atn
      integer i_dbq,h_dbq

      pname = "geom_create_full"
c     write(*,*) 'in ', pname

      if(.not.ma_push_get(mt_dbl,3*nfg,'c',h_c,i_c))
     & call errquit('mm: Failed to allocate memory for c',
     & 3*nfg, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nfg,'q',h_dbq,i_dbq))
     & call errquit('mm: Failed to allocate memory for q',nfg,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nfg,'m',h_m,i_m))
     & call errquit('mm: Failed to allocate memory for m',nfg,
     &       MA_ERR)
      if(.not.ma_push_get(mt_int,nfg,'inum',h_atn,i_atn))
     & call errquit('mm: Failed to allocate memory for atn',nfg,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*nfg,'t',h_t,i_t))
     & call errquit('mm: Failed to allocate memory for t',nfg,
     &       MA_ERR)
       
      call mm_get_geom(rtdb,nfg,int_mb(i_ifg),dbl_mb(i_c),
     &                 dbl_mb(i_dbq),dbl_mb(i_m),
     &                 int_mb(i_atn),byte_mb(i_t))

      if(.not.rtdb_cget(rtdb,'geometry',1,oldgeom))
     & oldgeom = ' '

      if(.not.rtdb_cput(rtdb,'geometry',1,'full_geom'))
     & call errquit(pname//' storing geom name to rtdb',0, RTDB_ERR)

c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t))
     & call errquit('mm: Failed to deallocate stack t_all',nfg,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_atn))
     & call errquit('mm: Failed to deallocate stack atn_all',nfg,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m))
     & call errquit('mm: Failed to deallocate stack m_all',nfg,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_dbq))
     & call errquit('mm: Failed to deallocate stack q_all',nfg,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit('mm: Failed to deallocate stack c_all',nfg,
     &       MA_ERR)


      end

      subroutine mm_get_geom(rtdb,nt,ifg,c,dbq,m,atn,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mm_bond_coords_data.fh"
#include "mm_coords_data.fh"
      
      integer rtdb
      integer geom
      integer nt
      integer ifg(nt)
      double precision c(3,nt)
      double precision dbq(nt)
      double precision m(nt)
      integer atn(nt)
      character*16 t(nt)
      
      integer i
      integer nat
      double precision scale
      logical status
      logical ignore
      character*30 message
      character*30 pname

      pname = "mm_get_geom"
      
c     if (.not. geom_create(geom, 'full_geom'))
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_create_geom',0, GEOM_ERR)

c     if (.not. geom_rtdb_load(rtdb, geom, 'full_geom'))
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_load_geom',0, RTDB_ERR)

      status=geom_ncent(geom,nat)

      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)

      if(.not.geom_cart_get2(geom,nat,t,c,dbq,atn))
     &    goto 911

      call util_convert_units("au","angstrom",scale)
      call dscal(3*nat, scale,c,1)

      if(.not.geom_destroy(geom)) 
     &    goto 911

      ignore = geom_rtdb_delete(rtdb,"full_geom")

      if(.not.geom_create(geom,"full_geom"))
     & call errquit('mm: Failed to create geometry',0, GEOM_ERR)

      call mm_links_adjust(nt,ifg,atn,t,c,dbq)


      call util_convert_units("angstrom","au",scale)
      call dscal(3*nt, scale,c,1)

      if(.not.geom_cart_set(geom,nt,t,c,dbq))
     & call errquit('mm: Failed to initialize geometry',0, GEOM_ERR)

      call mm_atn_get_mass(nt,atn,m) 

      if(.not.geom_masses_set(geom,nt,m))
     & call errquit('mm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(geom)

      if(.not.geom_rtdb_store(rtdb,geom,"full_geom"))
     & call errquit('mm: Failed to store geom to rtdb',0, RTDB_ERR)

      if(.not.geom_destroy(geom))
     & call errquit('mm: Failed to destroy geometry',0, GEOM_ERR)

      return
911   call errquit("error "//trim(message),0,-1)
      end 

      subroutine mm_atn_get_mass(n,atn,m)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
      integer n
      integer atn(n)
      double precision m(n)

c     local variables
      integer i
      character*32 pname
      
      pname = "mm_atn_get_mass"
      
      do i=1,n
        if(.not.geom_atn_to_default_mass(atn(i),m(i)))
     &    call errquit(pname,0, GEOM_ERR)

      end do

      end 

      subroutine mm_geom_push_active(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_geom_data.fh"
#include "mm_bond_coords_data.fh"
#include "mm_coords_data.fh"

      integer rtdb
     
      integer geom
      integer ncent
      integer i_ctmp,h_ctmp
      integer i, j
      integer nt
      double precision scale
      character*32 pname
       
      pname = "mm_geom_push_active"
c     write(*,*) 'in ', pname

      nt = nact

c     --------------------------------------
c     get active coordinates out of geometry
c     --------------------------------------
      if(.not.geom_create(geom,'geometry'))
     & call errquit('mm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(rtdb,geom,"geometry"))
     & call errquit('mm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(geom, ncent) ) 
     &    call errquit("mm:geom_ncent",0,0)

      if(.not.ma_push_get(mt_dbl,3*ncent,'ctmp',h_ctmp,i_ctmp))
     & call errquit( pname//'Failed to allocate memory for ctmp',
     & 3*ncent, MA_ERR)

      if(.not. geom_cart_coords_get(geom,dbl_mb(i_ctmp))) 
     &    call errquit("mm:geom_cart_coords_get",0,0)
 
 
      call util_convert_units("au","angstrom",scale)
      call dscal(3*ncent, scale,dbl_mb(i_ctmp),1)
      call mm_set_coord(nt,int_mb(i_iact),dbl_mb(i_ctmp))
      call mm_vdw_qmcoords_load(rtdb)     
      call mm_vdw14_qmcoords_load(rtdb)     
      call mm_bond_qmcoords_load(rtdb)     

      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('mm: 
     &              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

      if(.not.geom_destroy(geom))
     & call errquit('mm: Failed to destroy geometry',0, GEOM_ERR)

  
      end     

      subroutine mm_set_coord(nt,ai,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "mm_bond_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"

      integer rtdb
      integer nt
      integer ai(nt)
      double precision c(3,nt)

      integer i,j,k
      integer indx, iqm, imm
      character*32 pname 

      pname = "mm_set_coord"
c     write(*,*) 'in ', pname

      do i=1,nt
        do j=1,nqml
          indx = int_mb(i_iqml+j-1)
          if(ai(i).eq.indx) then
            do k = 1,3
              dbl_mb(i_rqml+3*(j-1)+k-1) = c(k,i)
            end do
          end if
        end do
      end do

c     -------------------
c     update i_rqm array
c     -------------------

      do i=1,nqm
        iqm = int_mb(i_iqm+i-1)
        do j=1,nqml
          indx = int_mb(i_iqml+j-1)
          if(iqm.eq.indx) then 
            do k = 1,3
              dbl_mb(i_rqm+3*(i-1)+k-1) = dbl_mb(i_rqml+3*(j-1)+k-1)
            end do
            exit
          end if
        end do
      end do

c     -------------------
c     update i_rmm array
c     -------------------
      if(qmlink) then 
        do i=1,nlink
          indx = int_mb(i_lb+nlink+i-1)
          do j=1,nmm
            imm = int_mb(i_imm+j-1)
            if(imm.eq.indx) then 
              do k = 1,3
                dbl_mb(i_rmm+3*(j-1)+k-1) = dbl_mb(i_rqml+
     >                                      (nqm+i-1)*3+k-1)
              end do
              exit 
            end if
          end do
        end do
      end if

      end 

      subroutine mm_geom_restore(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mm_geom_data.fh"
      integer rtdb
      logical ignore
      character*32 pname

      pname = "mm_restore_geom"
c     write(*,*) pname
      
      if(.not.aux_geom) return

      ignore = rtdb_delete(rtdb,'geometry')

      if(oldgeom.ne.' ') then
        if(.not.rtdb_cput(rtdb,'geometry',1,oldgeom))
     &    call errquit(pname//' storing geom to rtdb', 0, RTDB_ERR)
      end if
      
      end

      subroutine mm_geom_create_xyz_file(rtdb)
      implicit none
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"

      integer rtdb

      integer geom
      character*30 pname
      character*30 geomname
      character*30 operation
      character*50 filename
      character*255 full_filename
      integer ncent
      logical status
      logical ignore

      pname = "mm_geom_create_xyz_file"
      geomname = "geometry"
      

      if (ga_nodeid().eq.0) then
        call util_file_prefix('opt.xyz',filename)
        call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)
         
        if(.not.geom_create(geom,'geom_tmp'))
     >    call errquit('mm: Failed to create geometry',0, GEOM_ERR)
        if(.not.geom_rtdb_load(rtdb,geom,geomname))
     >    call errquit('mm: Failed to load geometry',0, GEOM_ERR)
        status =  geom_ncent(geom,ncent)
        open(88,file=full_filename,form='formatted')
        if (.not. geom_print_xyz(geom, 88))
     $    call errquit('mm:geom_print_xyz?',0, GEOM_ERR)
        close(88)

        if(.not.geom_destroy(geom))
     >    call errquit('mm: Failed to destroy geometry',0, GEOM_ERR)

      endif

      end

      subroutine mm_map_fixed_atoms(rtdb)
      implicit none
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"

      integer rtdb

      character*30 pname
      integer i,j,k
      integer nact0, nact1
      integer i_cons,h_cons
      integer i_act,h_act
      integer ma_type
      integer jj0,jj,ii,ii0
      integer nal,h_al,i_al
      integer h_am,i_am
      logical aflag
       
      pname = "mm_map_fixed_atoms"
c     write(*,*) pname

c     ---------------
c     get active list            
c     ---------------
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact0, h_cons)) then
        if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)

        if (.not.rtdb_delete(rtdb, 'geometry:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)
        do i=1,nact0
        end do
        if (.not.rtdb_put(rtdb, 'qmmm:actlist',
     >       mt_int,nact0,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      else
        nact0 = nact
        if(.not.ma_alloc_get(mt_int, nact0, 'qmmm actlist',
     &      h_cons, i_cons) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nact, MA_ERR)
        do i=1,nact0
           int_mb(i_cons+i-1) = i
        end do
      end if
      
c    ----------------------
c    create new active list
c    ----------------------
      if(.not.ma_push_get(MT_INT, nfg, 'qmmm fixed atom list',
     &      h_act, i_act) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nact, MA_ERR)

c     if(.not.ma_push_get(MT_LOG, ng, 'tmp qmmm act atom list',
c    &      h_am, i_am) ) call errquit(
c    &      pname//' unable to allocate stack',
c    &      ng, MA_ERR)

c     TP: call following subroutine to map active atom.
c     For now, assume all (qm) atoms are active.
     
c       call mm_cons_get_map(ng,int_mb(i_ig),
c    >                       log_mb(i_am))

c     --------------------------------------------
c     find total number (nal) of link atoms
c     that are bonded to active qm atoms
c     and store their global index in i_al
c     Note that link atoms are always indexed
c     by the global index of the classical atom 
c     --------------------------------------------
      if(nlink.ne.0) then

        if(.not.ma_push_get(MT_INT, nlink, 'qmmm tmp link',
     &        h_al, i_al) ) call errquit(
     &        pname//' unable to allocate stack',
     &        nact, MA_ERR)

        nal=0
        do i=1,nact0
c         get global index next active atom in the auxilary geometry
          ii0 = int_mb(i_cons+i-1)
          ii  = int_mb(i_iqml+ii0-1)
          aflag = .true.
          if(aflag) then
            do j=1,nlink
c             get index of qm atom bonded to a link atom
              jj  = int_mb(i_lb+j-1)
c             if qm atom bonded to a link atom is active
c             store corresponding index of the link atom
              if(aflag.and.(jj.eq.ii)) then
                nal=nal+1
                int_mb(i_al+nal-1)=int_mb(i_lb+nlink+j-1)
              end if
            end do 
          end if
        end do
      else
        nal = 0
      end if

c     ------------------------------------------------

c     ------------------------------------------------
c     construct active atom index
c     note that if classical boundary atom is
c     active the corresponding link atom is automatically
c     flagged as active in the first do loop because it 
c     carries the same global index. The second do loop takes
c     care of the case when qm boundary atom is active using
c     link index constructed above
c     ------------------------------------------------
      nact1=0
      do i=1,nfg
        ii = int_mb(i_iqml+i-1)
        do j=1,nact0
         jj0 = int_mb(i_cons+j-1)
         jj  = int_mb(i_iqml+jj0-1)
         aflag = .true.
         if(aflag.and.(ii.eq.jj)) then
            nact1 = nact1 + 1
            int_mb(i_act+nact1-1) = i
            goto 1
          end if
        end do
        do j=1,nal
         jj  = int_mb(i_al+j-1)
         if(ii.eq.jj) then
            nact1 = nact1 + 1
            int_mb(i_act+nact1-1) = i
            goto 1
          end if
        end do
1       continue
      end do

      if(nact1.ne.0) then
        if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >          mt_int,nact1,int_mb(i_act)))
     >         call errquit(pname,0,
     >         RTDB_ERR)
      end if

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)


      if (.not. ma_chop_stack(h_act) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)

      end

      subroutine mm_restore_fixed_atoms(rtdb)
      implicit none
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
      integer rtdb

      character*30 pname
      integer nact
      integer i_cons,h_cons
      integer ma_type
      logical ignore
     
      pname = "mm_restore_fixed_atoms"
c     write(*,*) pname

      ignore = rtdb_delete(rtdb, 'geometry:actlist')

      if (rtdb_ma_get(rtdb, 'qmmm:actlist', ma_type,
     $        nact, h_cons)) then
            if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)
      else
        return
      end if

      if (.not.rtdb_delete(rtdb, 'qmmm:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >       mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)
      
      end 
