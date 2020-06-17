c
c $Id$
c
      subroutine qmmm_geom_init(irtdb,ageom)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
#include "qmmm_params.fh"
      integer irtdb
      logical ageom
c     local variables
      integer nqm,nlink
      integer i_ind,i_lb
      integer mwa
      character*255 geom
      character*32 pname
      character*30 region(3)
      integer nregion
      character*84 tag

      pname = "qmmm_geom_init"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname


      nqm   = qmmm_get_nqm()
      nlink = qmmm_get_nlink()
      i_ind = qmmm_get_i_ai()
      i_lb  = qmmm_get_i_lb()
     
      call mm_get_mwa(mwa)

c
c     initialize indexing for a full qm geometry
c     ------------------------------------------
      nfg = nqm + nlink
      if(nfg.eq.0) goto 111
      if(.not.ma_alloc_get(mt_int,nfg,'qmmm fullg ind',h_ifg,i_ifg))
     + call errquit(pname//'Failed to allocate heap',nfg,
     &       MA_ERR)
      call icopy(nqm,int_mb(i_ind),1,int_mb(i_ifg),1)
      call icopy(nlink,int_mb(i_lb+nlink),1,int_mb(i_ifg+nqm),1)
c
c    allocate gradient array
c    -----------------------
      if(.not.ma_alloc_get(MT_DBL, 3*nfg, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*nfg, MA_ERR)
      call dfill(3*nfg,0,dbl_mb(i_g),1)
c
c    allocate esp charges
c    --------------------
      if(.not.ma_alloc_get(MT_DBL, 3*nfg, 'qmmm esp charges',
     &      h_qs, i_qs) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*nfg, MA_ERR)
      call dfill(3*nfg,0,dbl_mb(i_qs),1)
c
      if (rtdb_cget(irtdb,"qmmm:geom_file",1,geom)) then
        call mm_set_solute_coord_file(geom)
      end if
c
c     create active geometry if any
c     -----------------------------
111   continue
      if(ageom) then
c
c       process region definitions
c       --------------------------
        tag ="qmmm:region"
        if (.not.rtdb_get(irtdb,tag(1:inp_strlen(tag))//"_n",
     >                   mt_int,1,nregion)) 
     >        call errquit(pname//tag,0,RTDB_ERR)
        if(nregion.gt.0) then
          if (.not.rtdb_cget(irtdb,tag,nregion,region)) 
     >          call errquit(pname//tag,0,RTDB_ERR)
          call qmmm_cons_free_exclusive(region(1))
        end if
c
        call qmmm_geom_init_aux(irtdb,"geometry")
        call qmmm_geom_print(irtdb,"geometry","-active.pdb")
c        call qmmm_bq_data_update_active()

      else
        ngs = nfg
        ngw = 0
        ng =  nfg
        i_ig = i_ifg
        if(nfg.ne.0) then
          call qmmm_geom_create_full(irtdb)
          call qmmm_geom_print(irtdb,"geometry","-active.pdb")
        end if
      end if


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_init_aux(irtdb,aname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
#include "qmmm_params.fh"
      integer irtdb
      character*(*) aname
c     local variables
      integer i
      character*255 geom
      character*32 pname

      pname = "qmmm_geom_init_aux"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      if(aux_geom) call qmmm_geom_aux_end()

      call qmmm_cons_get_nacts(ngs)
      ngw = 0
      ng = ngs+ngw

      if(.not.ma_alloc_get(mt_int,ng,'qmmm geom index',h_ig,i_ig))
     + call errquit(pname//'Failed to allocate heap',ng,
     &       MA_ERR)

      call qmmm_cons_get_i_acts(int_mb(i_ig))

      if (rtdb_cget(irtdb,"qmmm:geomrtdb",1,geom)) then
        call qmmm_geom_push(irtdb,
     >                        ngs,
     >                        int_mb(i_ig),
     >                        ngw,
     >                        int_mb(i_ig+ngs),
     >                        geom)

        if(.not.rtdb_delete(irtdb,"qmmm:geomrtdb"))
     >     call errquit(pname//' qmmm:geomrtdb',0,
     &       RTDB_ERR)
      end if

      if(aname.ne."none")
     +  call qmmm_geom_create_active(irtdb,aname)

      if(qmmm_print_debug()) then
         write(*,*) "active solute index"
         do i=1,ngs
           write(*,*) int_mb(i_ig+i-1)
         end do
      end if

      aux_geom = .true.

      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_init_neb(irtdb)
c      initalizes geometry structures needed
c      for neb runs
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
#include "qmmm_params.fh"
#include "msgids.fh"
      integer irtdb
c     local variables
      character*255 filename
      character*255 filename0
      character*255 region
      character*255 prefix
      character*255 rstfile(2)
      character*255 geomname
      character*32 pname
      integer nbeads,i
      logical ofile
      logical custom_path

      pname = "qmmm_geom_init_neb"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

c     get neb region definition
c     -------------------------
      if (.not. rtdb_cget(irtdb,"qmmm:neb_region",
     +                    1,region)) then
       region="qmlink"
       if (.not. rtdb_cput(irtdb,"qmmm:neb_region",
     +                     1,region)) 
     +  call errquit(pname//'Failed neb region',
     + 1, RTDB_ERR)
      end if
c
        filename0 = ' '
        call mm_restart_filename(filename0)
c
c     activate atoms in neb region only
c     ---------------------------------
      call qmmm_cons_free_exclusive(region)
c
c     initialize auxiliary geometry
c     but do not create geometry object itself
c     ----------------------------------------
      call qmmm_geom_init_aux(irtdb,"none")
c
c     get rstfile names for the starting and ending
c     point of neb path
c     ---------------------------------------------
      if (rtdb_cget(irtdb,"qmmm:neb_path_limits",2,rstfile)) then 
        call qmmm_geom_create_rst(irtdb,"neb_start",rstfile(1))
        call qmmm_geom_create_rst(irtdb,"neb_end",rstfile(2))
        call mm_reload_rst(irtdb,filename0)
        custom_path = .false.
        if(.not.rtdb_put(irtdb,'neb:custom_path',mt_log,1,custom_path)) 
     >         call errquit(pname,0,RTDB_ERR)
      else
       custom_path = .true.
       if(.not.rtdb_put(irtdb,'neb:custom_path',mt_log,1,custom_path)) 
     >        call errquit(pname,0,RTDB_ERR)
       if (.not.rtdb_get(irtdb,'neb:nbeads',mt_int,1,nbeads)) then
          nbeads = 5
       end if
       if (.not.rtdb_cget(irtdb,"qmmm:neb_path",1,prefix))
     >   call mm_restart_filename_prefix(prefix)
       do i=1,nbeads
         filename = ' '
         write(filename,13) 
     >      prefix(1:inp_strlen(prefix)),i
         call util_file_name_resolve(filename,.false.)
         if(ga_nodeid().eq.0) then 
           inquire(file=filename,exist=ofile)
           if(ofile) then 
            write(*,*) "neb: found restart file",
     >                filename(1:inp_strlen(filename))
           call util_file_copy(filename(1:inp_strlen(filename)),
     >                         filename0(1:inp_strlen(filename0)))
           end if
         end if
         call ga_brdcst(msg_qmmm_misc, ofile,
     $        MA_sizeof(MT_INT,1,MT_BYTE), 0)
         call ga_sync()
         if(ofile) then
           geomname = ' '
           write(geomname,14) i 
           call qmmm_geom_create_rst
     >             (irtdb,
     >              geomname,
     >              filename0(1:inp_strlen(filename0)))
         else
           call errquit(pname//"cannot find"//
     >          filename(1:inp_strlen(filename)),0,0)
         end if
       end do
      end if

      if(qmmm_print_debug())
     >    write(*,*) "out",pname

 13      format(a,i3.3,'.rst')
 14      format('neb_bead_',i6.6,":geom")

      end

      subroutine qmmm_geom_create_rst(irtdb,geomname,rstfile)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
#include "qmmm_params.fh"
      integer irtdb
      character*(*) geomname
      character*(*) rstfile
c     local variables
      character*32 pname

      pname = "qmmm_geom_create_rst"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname
c
      call mm_reload_rst(irtdb,rstfile)
      call qmmm_geom_create_active(irtdb,geomname)

      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_aux_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_geom_data.fh"
      
      if(.not.aux_geom) return
       if(.not.ma_free_heap(h_ig))
     &   call errquit('Failed to deallocate heap',ng,
     &       MA_ERR)
       aux_geom = .false.
       ngs = nfg
       ngw = 0
       ng =  nfg
       i_ig = i_ifg

       return
      end

      subroutine qmmm_geom_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_geom_data.fh"
      
      if(aux_geom) call qmmm_geom_aux_end()

      if(nfg.ne.0) then
       if(.not.ma_free_heap(h_qs))
     &   call errquit('qmmm h_qs: Failed to deallocate heap',h_qs,
     &       MA_ERR)
       if(.not.ma_free_heap(h_g))
     &   call errquit('qmmm h_g: Failed to deallocate heap',h_g,
     &       MA_ERR)
       if(.not.ma_free_heap(h_ifg))
     &   call errquit('qmmm h_ifg: Failed to deallocate heap',h_ifg,
     &       MA_ERR)
       nfg = 0
      end if


       return
      end

      subroutine qmmm_geom_create(irtdb,ns,ais,nw,aiw,aname,link)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "stdio.fh"
#include "inp.fh"
#include "msgids.fh"
      integer irtdb
      integer ns
      integer ais(ns)
      integer nw
      integer aiw(nw)
      character*(*) aname
      logical link
c
      integer igeom
      character*32 pname
      character*255 buf 
      integer i_q,h_q
      integer i_c,h_c
      integer i_t,h_t
      integer i_m,h_m
      integer i_atn,h_atn
      integer n
      logical ignore
      integer geomlen
      logical oprint_low
      logical master
      character*255 filename
      character*16 tag
      integer j,ai0,fn
 
      pname = "qmmm_geom_create"
      if(qmmm_print_debug())
     > write(*,*) "in "//pname

      master = ga_nodeid().eq.0

      oprint_low = qmmm_print_low()
      geomlen = inp_strlen(aname)

      n = ns+nw
      if(.not.ma_push_get(mt_dbl,3*n,'c',h_c,i_c))
     + call errquit('qmmm: Failed to allocate memory for c',
     + 3*n, MA_ERR)
      if(.not.ma_push_get(mt_dbl,n,'q',h_q,i_q))
     + call errquit('qmmm: Failed to allocate memory for q',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,n,'m',h_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_int,n,'inum',h_atn,i_atn))
     + call errquit('qmmm: Failed to allocate memory for atn',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*n,'t',h_t,i_t))
     + call errquit('qmmm: Failed to allocate memory for t',n,
     &       MA_ERR)

      if(ns.ne.0) 
     >  call mm_get_solute_geom(ns,
     >                 ais,
     >                 int_mb(i_atn),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))


      if(nw.ne.0) 
     >  call mm_get_solvent_geom(nw,
     >                   aiw,
     >                   int_mb(i_atn+ns),
     >                   byte_mb(i_t+16*ns),
     >                   dbl_mb(i_c+3*ns),
     >                   dbl_mb(i_q+ns),
     >                   dbl_mb(i_m+ns))

c    ----------------------------
c    begin custom names from file
c    ----------------------------

      if (rtdb_cget(irtdb,"qmmm:rename",1,filename)) then

      if(master)  
     +  call util_print_centered(luout,
     + "Loading custom atom names from  "//
     +   filename,40,.true.)

      if(master) then
        if(.not.qmmm_get_io_unit(fn)) 
     >       call errquit("cannot get file number",0,0)
c
        open(fn,file=filename,form='formatted',status='old',
     $          err=133)
   
        fileloop: do 
        read(fn,*,err=133,end=20) tag,ai0
        write(*,*) "debug ", tag, ai0
        do j=1,ns
          if(ai0.eq.ais(j)) then
             call util_set_ma_char_arrays(16,tag,byte_mb(i_t+16*(j-1)))
             cycle fileloop
          end if
        end do
        end do fileloop
20      continue
        close(fn)
      end if

      call ga_brdcst(msg_qmmm_ind,
     >     byte_mb(i_t),
     >     16*n*ma_sizeof(mt_byte,1,mt_byte),
     >     0)

      end if
c    ----------------------------
c    end custom names from file
c    ----------------------------

 
      ignore = rtdb_delete(irtdb,aname(1:geomlen))

      if(.not.geom_create(igeom,aname(1:geomlen)))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(link)  
     +   call qmmm_links_adjust(ns,ais,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q))


      if(.not.geom_cart_set(igeom,n,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,n,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(.not.oprint_low.and.master) then
        buf = " "
        write(buf,*) "globally indexed xyz file for ",
     +    aname(1:geomlen)
        call util_print_centered(luout,
     +   buf,40,.true.)
        call qmmm_print_xyzi(luout,ns,ais,
     +      dbl_mb(i_c),byte_mb(i_t))
        if(.not.geom_print(igeom))
     +   call errquit('qmmm: Failed to print geom',0, RTDB_ERR)
      end if

      if(.not.geom_rtdb_store(irtdb,igeom,aname(1:geomlen)))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

      if(master) then
        call qmmm_print_xyzi(6,ns,ais,dbl_mb(i_c),byte_mb(i_t))
      end if

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t))
     & call errquit('qmmm: Failed to deallocate stack t_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_atn))
     & call errquit('qmmm: Failed to deallocate stack t_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m))
     & call errquit('qmmm: Failed to deallocate stack m_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q))
     & call errquit('qmmm: Failed to deallocate stack q_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit('qmmm: Failed to deallocate stack c_all',n,
     &       MA_ERR)

      if(qmmm_print_debug())
     > write(*,*) "out "//pname

      return
 133  call errquit('error reading file of names',0, 0)
 
      end

      subroutine qmmm_region_print(un,ns,nw,ais,aform)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
      integer un
      integer ns
      integer nw
      integer ais(ns+nw)
      character*(*) aform
c
      character*32 pname
      integer i_q,h_q
      integer i_c,h_c
      integer i_t,h_t
      integer i_m,h_m
      integer i_atn,h_atn
      integer n
 
      pname = "qmmm_region_print"
      if(qmmm_print_debug())
     > write(*,*) "in "//pname

      n = ns+nw
      if(n.eq.0) goto 999
      if(.not.ma_push_get(mt_dbl,3*n,'c',h_c,i_c))
     + call errquit('qmmm: Failed to allocate memory for c',
     + 3*n, MA_ERR)
      if(.not.ma_push_get(mt_dbl,n,'q',h_q,i_q))
     + call errquit('qmmm: Failed to allocate memory for q',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,n,'m',h_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_int,n,'inum',h_atn,i_atn))
     + call errquit('qmmm: Failed to allocate memory for atn',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*n,'t',h_t,i_t))
     + call errquit('qmmm: Failed to allocate memory for t',n,
     &       MA_ERR)

      if(ns.ne.0) 
     >   call mm_get_solute_geom(ns,
     >                 ais,
     >                 int_mb(i_atn),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))

      if(nw.ne.0) 
     >   call mm_get_solvent_geom(nw,
     >                 ais(ns+1),
     >                 int_mb(i_atn+ns),
     >                 byte_mb(i_t+16*ns),
     >                 dbl_mb(i_c+3*ns),
     >                 dbl_mb(i_q+ns),
     >                 dbl_mb(i_m+ns))


      if(qmmm_master()) then
        if(aform.eq."pdbi") then
        call qmmm_print_pdbi(un,n,
     >                      ais,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q),
     >                      byte_mb(i_t))
        else if (aform.eq."xyz") then
        call qmmm_print_xyz(un,n,
     >                      dbl_mb(i_c),
     >                      byte_mb(i_t))

        else if (aform.eq."xyzi") then
        call qmmm_print_xyzi(un,n,ais,
     >                      dbl_mb(i_c),
     >                      byte_mb(i_t))
        end if
      end if
c
c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t))
     & call errquit('qmmm: Failed to deallocate stack t_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_atn))
     & call errquit('qmmm: Failed to deallocate stack t_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m))
     & call errquit('qmmm: Failed to deallocate stack m_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q))
     & call errquit('qmmm: Failed to deallocate stack q_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit('qmmm: Failed to deallocate stack c_all',n,
     &       MA_ERR)

999   continue

      if(qmmm_print_debug())
     > write(*,*) "out "//pname
 
      end

      subroutine qmmm_geom_create_full(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
      integer irtdb
c     local variables
      character*32 pname
      logical update

      pname = "qmmm_geom_create_full"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      if(qmmm_print_debug())
     > write(*,*) "nfg=",nfg

      if(nfg.eq.0) return

      if (.not.rtdb_get(irtdb,"qmmm:update_geom",
     >                   mt_log,1,update)) 
     >   update = .true.
      
      if(update)  
     >  call qmmm_geom_create(irtdb,
     >                      nfg,
     >                      int_mb(i_ifg),
     >                      0,
     >                      0,
     >                      "full_geometry",.true.)


      if(.not.rtdb_cget(irtdb,'geometry',1,oldgeom))
     + oldgeom = ' '

      if(.not.rtdb_cput(irtdb,'geometry',1,'full_geometry'))
     + call errquit(pname//' storing geom name to rtdb',0, RTDB_ERR)

c      update = .false.
c      if (.not.rtdb_put(irtdb,"qmmm:update_geom",
c     >                   mt_log,1,update)) 
c     >   call errquit(pname//'qmmm:update_geom',0, RTDB_ERR)
     
      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_create_active(irtdb,aname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
      integer irtdb
      character*(*) aname
c     local variables
      character*32 pname

      pname = "qmmm_geom_create_active"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      call qmmm_geom_create(irtdb,
     >                      ngs,
     >                      int_mb(i_ig),
     >                      ngw,
     >                      int_mb(i_ig+ngs),
     >                      aname,.false.)


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_push_active(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
      integer irtdb
c     local variables
      character*32 pname

      pname = "qmmm_geom_push_active"

      if(.not.aux_geom) return

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      call qmmm_geom_push(irtdb,
     >                      ngs,
     >                      int_mb(i_ig),
     >                      ngw,
     >                      int_mb(i_ig+ngs),
     >                      "geometry")


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_push_gen(irtdb,aname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
      integer irtdb
      character*(*) aname
c     local variables
      character*32 pname

      pname = "qmmm_geom_push_active"

      if(.not.aux_geom) return

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      call qmmm_geom_push(irtdb,
     >                      ngs,
     >                      int_mb(i_ig),
     >                      ngw,
     >                      int_mb(i_ig+ngs),
     >                      aname)


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_restore(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_geom_data.fh"
#include "bq.fh"
      integer irtdb
c     local variables
      logical ignore
      character*32 pname
c
      pname = "qmmm_restore_qm_geom"

      if(.not.aux_geom) return

      if(qmmm_print_debug())
     > write(*,*) " in "//pname

      ignore = rtdb_delete(irtdb,'geometry')

      if(oldgeom.ne.' ') then
        if(.not.rtdb_cput(irtdb,'geometry',1,oldgeom))
     +    call errquit(pname//' storing geom name to rtdb',0, RTDB_ERR)
      end if

      if(qmmm_print_debug())
     > write(*,*) " out "//pname

      end

      subroutine qmmm_geom_push(irtdb,ns,ais,nw,aiw,aname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
      integer irtdb
      integer ns
      integer ais(ns)
      integer nw
      integer aiw(nw)
      character*(*) aname    
c     local variables
      integer ntot
      integer ncent


      integer igeom
      integer i_ctmp,h_ctmp
      character*32 pname

      pname = "qmmm_geom_push"

      ntot = ns+nw

c    --------------------------------------
c    get qm coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,aname))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(.not.ma_push_get(mt_dbl,3*ncent,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ncent, MA_ERR)

      if(.not. geom_cart_coords_get(igeom,dbl_mb(i_ctmp))) 
     >    call errquit("qmmm:geom_cart_coords_get",0,0)

      if(ncent.ne.ntot) 
     >    call errquit(pname//"mismatch in geometries",
     >                 ncent+10000*ntot,0)

      call mm_set_solute_coord(ns,
     >                         ais,
     >                         dbl_mb(i_ctmp))

      call mm_set_solvent_coord(nw,
     >                         aiw,
     >                         dbl_mb(i_ctmp+3*ns))


      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end

      subroutine qmmm_geom_print(irtdb,gname,fname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
      integer irtdb
      character*(*) gname    
      character*(*) fname    
c     local variables
      integer ncent
      integer igeom
      integer i_c,h_c,i_q,h_q,i_t,h_t
      character*32 pname

      pname = "qmmm_geom_push"

c    --------------------------------------
c    get qm coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(.not.ma_push_get(mt_dbl,3*ncent,'ctmp',h_c,i_c))
     + call errquit( pname//'Failed to allocate memory for c',
     + 3*ncent, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ncent,'qtmp',h_q,i_q))
     + call errquit('qmmm: Failed to allocate memory for q',ncent,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*ncent,'ttmp',h_t,i_t))
     + call errquit('qmmm: Failed to allocate memory for t',ncent,
     &       MA_ERR)

      if(.not.geom_cart_get(igeom,ncent,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to get geometry',0, GEOM_ERR)


      if(qmmm_master()) then
        call qmmm_print_pdb(ncent,fname,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q),
     >                       byte_mb(i_t))
      end if

      if(.not.ma_pop_stack(h_t))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_q))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_c))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end

      subroutine qmmm_geom_get_grad(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "qmmm_geom_data.fh"
      integer irtdb
c
      integer n,ns,nw
      integer i_igt,i_gt,h_gt

      if(aux_geom) then
        ns=ngs
        nw=ngw
        i_igt = i_ig
      else
        ns = ng
        nw = 0
        i_igt = i_ifg
      end if

      n = ns+nw
      
      if(.not.ma_push_get(MT_DBL, 3*n, 'qmmm grad array',
     &      h_gt, i_gt) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*n, MA_ERR)
      call dfill(3*n,0,dbl_mb(i_gt),1)


      if(ns.ne.0) 
     &  call mm_get_solute_force(ns,
     &                     int_mb(i_igt),
     &                     dbl_mb(i_gt))

      if(nw.ne.0) 
     &  call mm_get_solvent_force(nw,
     &                     int_mb(i_igt+ns),
     &                     dbl_mb(i_gt+3*ns))

      
      if(.not.rtdb_put(irtdb,'qmmm:gradient',mt_dbl,
     +                 3*n,dbl_mb(i_gt)))
     +   call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      if(.not.ma_pop_stack(h_gt))
     &      call errquit(
     &      'qmmm_data_alloc: unable to pop stack',
     &      3*n, MA_ERR)
 
      end

      subroutine qmmm_geom_print_current(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh" 
#include "geom.fh"
#include "global.fh"
      integer irtdb
      integer igeom

      character*23 pname

      pname = "qmmm_geom_print_current"
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,"geometry"))
     >    call errquit(pname//":load geometry",0,0)

c     print geometry
c     -------------
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,"Current QM/MM geometry", 36, .true.)
      if(.not.geom_print(igeom))
     +    call errquit('qmmm: Failed to print geometry',0)
      end if

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)



      end

      function qmmm_geom_get_nfg()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_nfg

      qmmm_geom_get_nfg = nfg

      end 

      function qmmm_geom_get_ng()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_ng

      qmmm_geom_get_ng = ng

      end 

      function qmmm_geom_get_ngs()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_ngs

      qmmm_geom_get_ngs = ngs

      end 

      function qmmm_geom_get_ngw()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_ngw

      qmmm_geom_get_ngw = ngw

      end 

      function qmmm_geom_get_i_ig()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_i_ig

      qmmm_geom_get_i_ig = i_ig

      end 

      function qmmm_geom_get_i_g()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_i_g

      qmmm_geom_get_i_g = i_g

      end 

      function qmmm_geom_get_i_ifg()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_i_ifg

      qmmm_geom_get_i_ifg = i_ifg

      end 

      function qmmm_aux_geom()
      implicit none
#include "qmmm_geom_data.fh"

      logical qmmm_aux_geom

      qmmm_aux_geom = aux_geom

      end 
