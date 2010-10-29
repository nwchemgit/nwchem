c
c $Id$
c
c
      subroutine qmmm_cons_init(irtdb)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
#include "qmmm_cons_data.fh"
      integer irtdb
c
      character*32 pname
      logical fixw
      logical fixs
      logical freeqmlink
      logical freeqm
      logical fixqm
      logical fixqmlink
c
      integer i
      character*26 ma_date
      integer ma_type
      integer nr
      integer h_ires,i_ires
c
c     get total number of active atoms and global index
c     as specified by mm input
c     ------------------------------------------------
      call mm_get_solute_tot_na_gen(nacts_qm,mm_active_qm)
      call mm_get_solute_tot_na_gen(nacts_link,mm_active_link)
      call mm_get_solute_tot_na_gen(nacts_mm,mm_active_classical)
      call mm_get_solute_tot_nact(nacts)
      call mm_get_solvent_tot_nact(nactw)

      if(.not.ma_alloc_get(MT_INT, nacts, 'qmmm acts',
     &      h_acts, i_acts) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nacts, MA_ERR)
      call ifill(nacts,-1,int_mb(i_acts),1)

      if(.not.ma_alloc_get(MT_INT, nactw, 'qmmm actw',
     &      h_actw, i_actw) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nacts, MA_ERR)
      call ifill(nactw,-1,int_mb(i_actw),1)
c
      if(nacts_qm.ne.0)  
     >  call mm_get_solute_ind_gen(nacts_qm,
     >                           mm_active_qm,
     >                           int_mb(i_acts))
c 
      if(nacts_link.ne.0)  
     >  call mm_get_solute_ind_gen(nacts_link,
     >                           mm_active_link,
     >                           int_mb(i_acts+nacts_qm))
c
      if(nacts_mm.ne.0)  
     >  call mm_get_solute_ind_gen(nacts_mm,
     >                    mm_active_classical,
     >                    int_mb(i_acts+nacts_qm+nacts_link))

      call mm_get_solvent_ind_act(nactw,int_mb(i_actw))

      if(.not.ma_alloc_get(MT_LOG, nacts, 'qmmm maps',
     &      h_maps, i_maps) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nacts, MA_ERR)

      if(.not.ma_alloc_get(MT_LOG, nactw, 'qmmm mapw',
     &      h_mapw, i_mapw) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nactw, MA_ERR)

      do i=1,nacts
        log_mb(i_maps+i-1) = .true.
      end do

      do i=1,nactw
        log_mb(i_mapw+i-1) = .true.
      end do

      if (.not. rtdb_get(irtdb,'qmmm:fixsolvent',mt_log,1,fixw))
     $     fixw=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixsolute',mt_log,1,fixs))
     $     fixs=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqmlink',mt_log,1,freeqmlink))
     $     freeqmlink=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqm',mt_log,1,freeqm))
     $     freeqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqm',mt_log,1,fixqm))
     $     fixqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqmlink',mt_log,1,fixqmlink))
     $     fixqmlink=.false.


      if(fixw) call qmmm_cons_set("fix","solvent")
      if(fixs) call qmmm_cons_set("fix","mm_solute")
      if(fixqm) call qmmm_cons_set("fix","qm")
      if(fixqmlink) call qmmm_cons_set("fix","qmlink")
      if(freeqmlink) call qmmm_cons_set("free","qmlink")
      if(freeqm) call qmmm_cons_set("free","qm")

      return
      end
c
      subroutine qmmm_cons_init0(irtdb)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
#include "qmmm_cons_data.fh"
      integer irtdb
c
      character*32 pname
      logical fixw
      logical fixs
      logical freeqmlink
      logical freeqm
      logical fixqm
      logical fixqmlink
c
      integer i
      character*26 ma_date
      integer ma_type
      integer nr
      integer h_ires,i_ires
c
c     get total number of active atoms and global index
c     as specified by mm input
c     ------------------------------------------------
      call mm_get_solute_tot_na_gen(nacts_qm,mm_active_qm)
      call mm_get_solute_tot_na_gen(nacts_link,mm_active_link)
      call mm_get_solute_tot_na_gen(nacts_mm,mm_active_classical)
      call mm_get_solute_tot_nact(nacts)
      call mm_get_solvent_tot_nact(nactw)

      if(.not.ma_alloc_get(MT_INT, nacts, 'qmmm acts',
     &      h_acts, i_acts) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nacts, MA_ERR)
      call ifill(nacts,-1,int_mb(i_acts),1)

      if(.not.ma_alloc_get(MT_INT, nactw, 'qmmm actw',
     &      h_actw, i_actw) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nacts, MA_ERR)
      call ifill(nactw,-1,int_mb(i_actw),1)
c
      if(nacts_qm.ne.0)  
     >  call mm_get_solute_ind_gen(nacts_qm,
     >                           mm_active_qm,
     >                           int_mb(i_acts))
c 
      if(nacts_link.ne.0)  
     >  call mm_get_solute_ind_gen(nacts_link,
     >                           mm_active_link,
     >                           int_mb(i_acts+nacts_qm))
c
      if(nacts_mm.ne.0)  
     >  call mm_get_solute_ind_gen(nacts_mm,
     >                    mm_active_classical,
     >                    int_mb(i_acts+nacts_qm+nacts_link))

      call mm_get_solvent_ind_act(nactw,int_mb(i_actw))

      if(.not.ma_alloc_get(MT_LOG, nacts, 'qmmm maps',
     &      h_maps, i_maps) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nacts, MA_ERR)

      if(.not.ma_alloc_get(MT_LOG, nactw, 'qmmm mapw',
     &      h_mapw, i_mapw) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nactw, MA_ERR)

      do i=1,nacts
        log_mb(i_maps+i-1) = .true.
      end do

      do i=1,nactw
        log_mb(i_mapw+i-1) = .true.
      end do


      if (.not. rtdb_get(irtdb,'qmmm:fixsolvent',mt_log,1,fixw))
     $     fixw=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixsolute',mt_log,1,fixs))
     $     fixs=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqmlink',mt_log,1,freeqmlink))
     $     freeqmlink=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqm',mt_log,1,freeqm))
     $     freeqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqm',mt_log,1,fixqm))
     $     fixqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqmlink',mt_log,1,fixqmlink))
     $     fixqmlink=.false.

      if(.not.rtdb_get_info(irtdb, 'qmmm:fixres', ma_type,
     $        nr, ma_date))   
     $   nr = 0
       
      if(nr.ne.0) then
        if(.not.ma_push_get(mt_int,nr,'ires',h_ires,i_ires))
     +   call errquit( pname//'Failed to allocate memory for ires',
     +   nr, MA_ERR)
        if (.not. rtdb_get(irtdb,'qmmm:fixres',
     +            mt_int,nr,int_mb(i_ires)))
     +   call errquit( pname//'Failed ires',0,0)
      end if
     

      if(fixw) then
        do i=1,nactw
          log_mb(i_mapw+i-1) = .false.
        end do
      end if
c
      if(fixs) then
        do i=nacts_qm+nacts_link+1,nacts
          log_mb(i_maps+i-1) = .false.
        end do
      end if
c
      if(fixqm) then
        do i=1,nacts_qm
          log_mb(i_maps+i-1) = .false.
        end do
      end if
c
      if(fixqmlink) then
        do i=nacts_qm+1,nacts_qm+nacts_link
          log_mb(i_maps+i-1) = .false.
        end do
      end if
c      if(fixs) call mm_fix_solute(mm_classical)
c      if(fixqm) then
c         write(*,*) "fixing quantum atoms"
c         call mm_fix_solute(mm_quant)
c      end if
c      if(fixqmlink) then
c         write(*,*) "fixing quantum+link atoms"
c         call mm_fix_solute(mm_quant_link)
c      end if
c
c      if(nr.ne.0) call mm_fix_solute_res(nr,int_mb(i_ires))
c
c      if(freeqmlink) call mm_free_solute(mm_quant_link)
c      if(freeqm) call mm_free_solute(mm_quant)

      if(nr.ne.0) then
        if(.not.ma_pop_stack(h_ires))
     &   call errquit(pname//'stack deallocation h_ires',nr,
     &         MA_ERR)
      end if
      return
      end
c
      subroutine qmmm_cons_get_map(nt,ai,map)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
#include "qmmm_cons_data.fh"
      integer nt
      integer ai(nt)
      logical map(nt)
c
      character*32 pname
      integer i,j

      do j=1,nt
         map(j) = .false.
      end do
      do i=1,nacts
        do j=1,nt
        if(int_mb(i_acts+i-1).eq.ai(j)) then
          map(j) = log_mb(i_maps+i-1) 
        end if
        end do
      end do

      return
      end
c
      subroutine qmmm_cons_set(atype,aregion)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
#include "qmmm_cons_data.fh"
      character*(*) atype
      character*(*) aregion
c
      character*32 pname
      logical fixw
      logical fixs
      logical freeqmlink
      logical freeqm
      logical fixqm
      logical fixqmlink
c
      integer i
      logical flag

      if(atype.eq.'fix') then
         flag = .false.
      else if(atype.eq.'free') then
         flag = .true.
      else
         call errquit('unknown constraint type',
     +   0,0)
      end if
c
      if(aregion.eq.'solvent') then
        do i=1,nactw
          log_mb(i_mapw+i-1) = flag
        end do
      else if(aregion.eq.'solute') then
        do i=1,nacts
          log_mb(i_maps+i-1) = flag
        end do
      else if(aregion.eq.'all') then
        do i=1,nacts
          log_mb(i_maps+i-1) = flag
        end do
        do i=1,nactw
          log_mb(i_mapw+i-1) = flag
        end do
      else if(aregion.eq.'mm') then
        do i=1,nactw
          log_mb(i_mapw+i-1) = flag
        end do
        do i=nacts_qm+nacts_link+1,nacts
          log_mb(i_maps+i-1) = flag
        end do
      else if(aregion.eq.'mm_solute') then
        do i=nacts_qm+nacts_link+1,nacts
          log_mb(i_maps+i-1) = flag
        end do
      else if(aregion.eq.'qm') then
        do i=1,nacts_qm
          log_mb(i_maps+i-1) = flag
        end do
      else if(aregion.eq.'qmlink') then
        do i=1,nacts_qm+nacts_link
          log_mb(i_maps+i-1) = flag
        end do
      else if(aregion.eq.'link') then
        do i=nacts_qm+1,nacts_qm+nacts_link
          log_mb(i_maps+i-1) = flag
        end do
       end if

      return
      end
c
      subroutine qmmm_cons_free_exclusive(aregion)
      implicit none
c
#include "qmmm_cons_data.fh"
      character*(*) aregion
c
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",aregion)
c
      return
      end
c
      subroutine qmmm_cons_fixed()
      implicit none
c
#include "mafdecls.fh"
#include "qmmm.fh"
#include "qmmm_cons_data.fh"
      character*32 pname
 
      pname = "qmmm_cons_fixed"
c
      if(qmmm_print_debug()) 
     +  write(*,*) "in",pname
c
      call mm_zero_solute_force(nacts,
     >                          int_mb(i_acts),
     >                          log_mb(i_maps))
c
      call mm_zero_solvent_force(nactw,
     >                          int_mb(i_actw),
     >                          log_mb(i_mapw))
c
      if(qmmm_print_debug()) 
     +  write(*,*) "out",pname
c
      return
      end
c
      subroutine qmmm_cons_set0(irtdb)
      implicit none
c
#include "mm_utils.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
#include "qmmm_cons_data.fh"
      integer irtdb
c
      character*32 pname
      logical fixw
      logical fixs
      logical freeqmlink
      logical freeqm
      logical fixqm
      logical fixqmlink
c
      integer i
      character*26 ma_date
      integer ma_type
      integer nr
      integer h_ires,i_ires

      if (.not. rtdb_get(irtdb,'qmmm:fixsolvent',mt_log,1,fixw))
     $     fixw=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixsolute',mt_log,1,fixs))
     $     fixs=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqmlink',mt_log,1,freeqmlink))
     $     freeqmlink=.false.

      if (.not. rtdb_get(irtdb,'qmmm:freeqm',mt_log,1,freeqm))
     $     freeqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqm',mt_log,1,fixqm))
     $     fixqm=.false.

      if (.not. rtdb_get(irtdb,'qmmm:fixqmlink',mt_log,1,fixqmlink))
     $     fixqmlink=.false.

      if(.not.rtdb_get_info(irtdb, 'qmmm:fixres', ma_type,
     $        nr, ma_date))   
     $   nr = 0
       
      if(nr.ne.0) then
        if(.not.ma_push_get(mt_int,nr,'ires',h_ires,i_ires))
     +   call errquit( pname//'Failed to allocate memory for ires',
     +   nr, MA_ERR)
        if (.not. rtdb_get(irtdb,'qmmm:fixres',
     +            mt_int,nr,int_mb(i_ires)))
     +   call errquit( pname//'Failed ires',0,0)
      end if

      if(fixw) then
        do i=1,nactw
          log_mb(i_mapw+i-1) = .false.
        end do
      end if
c
      if(fixs) then
        do i=nacts_qm+nacts_link+1,nacts
          log_mb(i_maps+i-1) = .false.
        end do
      end if
c
      if(fixqm) then
        do i=1,nacts_qm
          log_mb(i_maps+i-1) = .false.
        end do
      end if
c
      if(fixqmlink) then
        do i=nacts_qm+1,nacts_qm+nacts_link
          log_mb(i_maps+i-1) = .false.
        end do
      end if

      if(nr.ne.0) then
        if(.not.ma_pop_stack(h_ires))
     &   call errquit(pname//'stack deallocation h_ires',nr,
     &         MA_ERR)
      end if
      return
      end
c
      subroutine qmmm_cons_reset()
      implicit none
c
#include "mafdecls.fh"
#include "qmmm_cons_data.fh"
c
      integer i

      do i=1,nacts
        log_mb(i_maps+i-1) = .true.
      end do

      do i=1,nactw
        log_mb(i_mapw+i-1) = .true.
      end do

      return
      end
c
      subroutine qmmm_cons_reaction(irtdb,do_grad,etot,g)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
      logical do_grad
      double precision g(3,*)
      double precision etot
c
      character*32 pname
      integer igeom 
      integer ncent
      integer i_c,h_c
      integer nfg
      integer i_g,i_ifg
      integer ns
      integer i_is,h_is
      integer i_as,h_as
      integer i_rs,h_rs
      double precision k0,r0
      double precision econs,deds,s
c
      if (.not. rtdb_get(irtdb,'qmmm:n_spring',mt_int,1,ns))
     $     ns = 0

      if(ns.eq.0) return

     
      pname = "qmmm_cons_reaction"

      if(qmmm_print_debug()) 
     +  write(*,*) "in",pname

      nfg = qmmm_geom_get_nfg()

      if(.not.ma_alloc_get(mt_int,2*ns,'qmmm rcons ind',h_is,i_is))
     + call errquit(pname//'Failed to allocate memory for is',
     + 2*ns, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,ns,'qmmm rcons r',h_rs,i_rs))
     + call errquit(pname//'Failed to allocate memory for rs',
     + ns, MA_ERR)
        call dfill(ns,0.0d0,dbl_mb(i_rs),1)

      if(.not.ma_alloc_get(mt_dbl,ns,'qmmm rcons coeff',h_as,i_as))
     + call errquit(pname//'Failed to allocate memory for as',
     + ns, MA_ERR)


      if (.not. rtdb_get(irtdb,'qmmm:i_spring',mt_int,
     +                    2*ns,int_mb(i_is)))
     + call errquit(pname//'Failed to allocate memory for is',
     + 3*ns, MA_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:a_spring',mt_dbl,
     +                    ns,dbl_mb(i_as)))
     +  call dfill(ns,1.0d0,dbl_mb(i_as),1)

      if (.not. rtdb_get(irtdb,'qmmm:k_spring',mt_dbl,
     +                    1,k0))
     + call errquit(pname//'Failed to get k0',
     + 1, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:r_spring',mt_dbl,
     +                    1,r0))
     + call errquit(pname//'Failed to get r0',
     + 1, RTDB_ERR)


c    --------------------------------------
c    get coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(nfg.ne.ncent) 
     + call errquit(pname//': geometry mismatch',0, GEOM_ERR)

      if(.not.ma_alloc_get(mt_dbl,3*ncent,'ctmp',h_c,i_c))
     + call errquit( pname//'Failed to allocate memory for c',
     + 3*ncent, MA_ERR)

      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c)))
     + call errquit(pname//': Failed to get geometry',0, GEOM_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)

      call qmmm_cons_spring_egrad(ns,nfg,
     >                       do_grad,
     >                       k0,
     >                       r0,
     >                       econs,
     >                       int_mb(i_is),
     >                       dbl_mb(i_as),
     >                       dbl_mb(i_rs),
     >                       dbl_mb(i_c),
     >                       deds,
     >                       s,
     >                       g)
      etot  = etot + econs
c
      if (.not. rtdb_put(irtdb,'qmmm:spring_e',mt_dbl,
     +                    1,econs))
     + call errquit(pname//'Failed to set spring_e',
     + 1, RTDB_ERR)

      if (.not. rtdb_put(irtdb,'qmmm:spring_s',mt_dbl,
     +                    1,s))
     + call errquit(pname//'Failed to set spring_s',
     + 1, RTDB_ERR)


      if (.not. rtdb_put(irtdb,'qmmm:spring_deds',mt_dbl,
     +                    1,deds))
     + call errquit(pname//'Failed to set spring_deds',
     + 1, RTDB_ERR)

      if(qmmm_master()) then
        write(*,*) "spring e,deds,s",econs,deds,s
      end if
 
      if(.not.ma_free_heap(h_c))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_as))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_rs))
     + call errquit(pname//'Failed to deallocate memory for rs',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_is))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      return
      end
c
      subroutine qmmm_cons_spring_egrad(ns,
     >                       nat,
     >                       do_grad,
     >                       k0,
     >                       r0,
     >                       econs,
     >                       ind,
     >                       a,
     >                       r,
     >                       c,
     >                       deds,
     >                       s,
     >                       g)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "mm_utils.fh"

      integer ns,nat
      logical do_grad
      double precision k0,r0
      integer ind(2,ns)
      double precision econs
      double precision deds
      double precision s
      double precision  a(ns)
      double precision  r(ns)
      double precision g(3,nat)
      double precision c(3,nat)

      integer is,k,i1,i2
      double precision r1(3),r2(3)
      double precision rr,e,g0,ehalf
      character*32 pname

      pname = "qmmm_cons_spring_egrad"
      do is=1,ns
        i1 = ind(1,is)
        i2 = ind(2,is)
        do k=1,3
          r1(k) = c(k,i1)
          r2(k) = c(k,i2)
        end do
        rr = 0.0d0
        do k=1,3
          rr = rr+(r1(k)-r2(k))**2
        end do
        r(is) = sqrt(rr)
      end do
      
      s = 0.0d0
      do is=1,ns
        s = s+a(is)*r(is)
      end do

      ehalf = s-r0
      econs = ehalf*ehalf*k0
      deds = -2.0*k0*ehalf


      if(.not.do_grad) return

      do is=1,ns
        i1 = ind(1,is)
        i2 = ind(2,is)
        do k=1,3
          r1(k) = c(k,i1)
          r2(k) = c(k,i2)
        end do 
        do k=1,3
          g0=2*k0*ehalf*a(is)*(r1(k)-r2(k))/r(is)
          g(k,i1) = g(k,i1)+g0
          g(k,i2) = g(k,i2)-g0
        end do 
      end do

      return 
      end
c
      subroutine qmmm_cons_get_nacts_qm(n)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"

      integer i,n

      n = 0
      do i=1,nacts_qm
        if( log_mb(i_maps+i-1)) n = n + 1
      end do

      end 
c
      subroutine qmmm_cons_get_nacts(n)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"

      integer i,n

      n = 0
      do i=1,nacts
        if( log_mb(i_maps+i-1)) n = n + 1
      end do

      end 
c
      subroutine qmmm_cons_get_nactw(n)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"

      integer i,n

      n = 0
      do i=1,nactw
        if( log_mb(i_mapw+i-1)) n = n + 1
      end do

      end 
c
      subroutine qmmm_cons_get_i_acts(ai)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"

      integer ai(*)

      integer i,n

      n = 0
      do i=1,nacts
        if( log_mb(i_maps+i-1)) then
          n = n + 1
          ai(n) = int_mb(i_acts+i-1)
        end if
      end do

      end 
c
      subroutine qmmm_cons_print_region(un,aform)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      integer un
      character*(*) aform

      integer i,ns,nw,n
      integer i_id,h_id
      character*30 pname

      pname = "qmmm_cons_print_region"

      call qmmm_cons_get_nacts(ns)
      call qmmm_cons_get_nactw(nw)
      n = ns+nw

      if(.not.ma_push_get(mt_int,n,'qmmm_cons i_id',h_id,i_id))
     + call errquit( pname//'Failed to allocate memory for id',
     + n, MA_ERR)
       call ifill(n,0,int_mb(i_id),1)

      call qmmm_cons_get_i_acts(int_mb(i_id))

      call qmmm_region_print(un,ns,nw,int_mb(i_id),aform)
       
      if(.not.ma_pop_stack(h_id))
     + call errquit( pname//'Failed to pop stack for id',
     + n, MA_ERR)
       

      end 

      subroutine qmmm_cons_print()
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"
      integer i

      write(*,*) "solute cons"
      do i=1,nacts
        write(*,*) int_mb(i_acts+i-1),log_mb(i_maps+i-1)
      end do

      write(*,*) "solvent cons"
      do i=1,nactw
        write(*,*)  int_mb(i_actw+i-1),log_mb(i_mapw+i-1)
      end do


      end 

      subroutine qmmm_cons_get_actmaps(n,ai,amap)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"

      integer n,ai(n)
      logical amap(n)

      integer i,j


      do j=1,n
       amap(j) = .false. 
      end do

      do i=1,nacts
        if( log_mb(i_maps+i-1)) then
        do j=1,n
          if(int_mb(i_acts+i-1).eq.ai(j)) then
            amap(j) = .true.
            go to 11
          end if
        end do
        end if
11      continue

      end do

      end 

      subroutine qmmm_cons_get_actmapw(n,ai,amap)
      implicit none
#include "qmmm_cons_data.fh"
#include "mafdecls.fh"

      integer n,ai(n)
      logical amap(n)

      integer i,j,m
      integer mwa

      call mm_get_mwa(mwa)
      do j=1,n
       amap(j) = .false. 
      end do

      do i=1,nactw
        if( log_mb(i_mapw+i-1)) then
        do j=1,n/mwa
          if(int_mb(i_actw+i-1).eq.ai(j)) then
            do m=1,mwa
            amap((j-1)*mwa+m) = .true.
            end do
            go to 11
          end if
        end do
        end if
11      continue

      end do

      end 

      subroutine qmmm_cons_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_cons_data.fh"

      character*32 pname
c
      pname = "qmmm_cons_end"
c
      if(.not.ma_free_heap(h_mapw))
     & call errquit(pname//'Failed to deallocate heap h_mapw',0,
     &       MA_ERR)

      if(.not.ma_free_heap(h_maps))
     & call errquit(pname//'Failed to deallocate heap h_maps',0,
     &       MA_ERR)

      if(.not.ma_free_heap(h_actw))
     & call errquit(pname//'Failed to deallocate heap h_actw',0,
     &       MA_ERR)


      if(.not.ma_free_heap(h_acts))
     & call errquit(pname//'Failed to deallocate heap h_acts',0,
     &       MA_ERR)


       return
      end
