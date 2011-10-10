      subroutine qmmm_check_forces1(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
#include "inp.fh"
c     
      integer rtdb
c
      integer ntot
      integer i_c,h_c
      integer i_c0,h_c0
      integer i_ai,h_ai
      integer i_g,h_g
      double precision c0,gref
      integer ia,i,k
      double precision e0, ep, em, grad, hess, step
      double precision step0
      integer i1,k1
c
      character*32 pname
      character*255 tag
      character*30 region(3)
      integer nregion
c
      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking forces'
         write(6,*)
         call util_flush(6)
      end if
c
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion))
     >      call errquit(pname//tag,0,RTDB_ERR)
      if(nregion.gt.3)
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(rtdb,tag,nregion,region))
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c

c
c     get total number of atoms
c     -------------------------
c      call qmmm_cons_get_nacts(ntot)
      ntot = 1
      if (ga_nodeid() .eq. 0) then
         write(6,*)
         write(6,*) '@ number of atoms',ntot
         write(6,*)
         call util_flush(6)
      end if
c
      if(.not.ma_alloc_get(MT_INT, ntot, 'qmmm ref index array',
     &      h_ai,i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      ntot, MA_ERR)
      call ifill(ntot,-1,int_mb(i_ai),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm ref grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm ref coord array',
     &      h_c,i_c) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_c),1)

c      call qmmm_cons_get_i_acts(int_mb(i_ai))
      int_mb(i_ai) = 10
      call mm_get_solute_coord_raw(ntot,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))

      call md_sp_qmmm()
      call prp_print()
        if (.not. rtdb_get(rtdb,'md:energy', mt_dbl, 1, e0))
     $     call errquit('driver: could not get energy',0, RTDB_ERR)
        call mm_get_solute_force_raw(ntot,
     &                    int_mb(i_ai),
     &                    dbl_mb(i_g))

      if(.not.rtdb_get(rtdb,'qmmm:step',mt_dbl,
     >                 1,step0))
     >   step0=0.01d0
      
       call mm_print_system()

c
c     loop over all atoms
c     -------------------
      if (ga_nodeid() .eq. 0) 
     >   write(6,8) "comp","anal-g","num-g","error","step","ep","em" 
 8       format(1x,'checkgrad ',A5,4A16,2x,2A12)

      i=0
      do ia=1,ntot
        do k=1,3
          i=i+1
          c0=dbl_mb(i_c+i-1)
          gref=dbl_mb(i_g+i-1)
          
          step = step0
10        continue
          dbl_mb(i_c+i-1)=c0+step
          call mm_set_solute_coord_raw(ntot,
     >                     int_mb(i_ai),
     >                     dbl_mb(i_c))

          call md_sp_qmmm()
          if (.not. rtdb_get(rtdb,'md:energy', mt_dbl, 1, ep))
     $     call errquit('driver: could not get energy',0, RTDB_ERR)
          call mm_print_system()
          call prp_print()
          write(*,*) "c,c0,ep",dbl_mb(i_c+i-1),c0,ep

c          if (abs(ep-e0) .lt. 1e-6) then
c             write(6,*) ' Increasing the step ', ep,e0,ep-e0, step
c             step = step*10.0d0
c             goto 10
c          else if (abs(ep-e0) .gt. 1e-2) then
c             write(6,*) ' Decreasing the step ', ep,e0,ep-e0, step
c             step = step/3.0
c             if(step.le.0.00000001d0) goto 12
c             goto 10
c          end if
12        continue

          dbl_mb(i_c+i-1)=c0-step
          call mm_set_solute_coord_raw(ntot,
     >                     int_mb(i_ai),
     >                     dbl_mb(i_c))

          call md_sp_qmmm()
          if (.not. rtdb_get(rtdb,'md:energy', mt_dbl, 1, em))
     $     call errquit('driver: could not get energy',0, RTDB_ERR)

          write(*,*) "c,c0,em",dbl_mb(i_c+i-1),c0,em

          grad = (ep - em) / (2.0d0*step)
c
          if (ga_nodeid() .eq. 0) 
     >      write(6,7) int_mb(i_ai+ia-1), 
     >         gref, grad,abs(gref-grad), step, ep,em
 7        format(1x,'checkgrad ',i5,2f16.8,e16.8,f16.8,2x,2f16.8)
          dbl_mb(i_c+i-1)=c0
        end do
      end do


      end

      subroutine qmmm_check_forces(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
#include "inp.fh"
c     
      integer rtdb
c
      integer ntot
      integer i_c,h_c
      integer i_c0,h_c0
      integer i_ai,h_ai
      integer i_g,h_g
      double precision c0,gref
      integer ia,i,k
      double precision e0, ep, em, grad, hess, step
      double precision step0
      integer i1,k1
c
      logical qmmm_energy_gradient
      logical task_qmmm_gradient
      external qmmm_energy_gradient
      external task_qmmm_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname
      character*255 tag
      character*30 region(3)
      integer nregion
c
      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking forces'
         write(6,*)
         call util_flush(6)
      end if
c
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion))
     >      call errquit(pname//tag,0,RTDB_ERR)
      if(nregion.gt.3)
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(rtdb,tag,nregion,region))
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c

      call mm_print_system()
c
c     get total number of atoms
c     -------------------------
      call qmmm_cons_get_nacts(ntot)
      if (ga_nodeid() .eq. 0) then
         write(6,*)
         write(6,*) '@ number of atoms',ntot
         write(6,*)
         call util_flush(6)
      end if
c
      if(.not.ma_alloc_get(MT_INT, ntot, 'qmmm ref index array',
     &      h_ai,i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      ntot, MA_ERR)
      call ifill(ntot,-1,int_mb(i_ai),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm ref grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm ref coord array',
     &      h_c,i_c) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_c),1)

      call qmmm_cons_get_i_acts(int_mb(i_ai))
      call mm_get_solute_coord(ntot,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))

      call md_sp_qmmm()
      if (.not. qmmm_energy_gradient(rtdb,.true.))
     $    call errquit(pname//'qmmm_energy_gradient failed'
     $                  ,0, GEOM_ERR)
        call qmmm_energy_rtdb_push(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy', mt_dbl, 1, e0))
     $     call errquit('driver: could not get energy',0, RTDB_ERR)
        call qmmm_print_energy(rtdb)
        call mm_get_solute_force(ntot,
     &                    int_mb(i_ai),
     &                    dbl_mb(i_g))

      if(.not.rtdb_get(rtdb,'qmmm:step',mt_dbl,
     >                 1,step0))
     >   step0=0.01d0

c
c     loop over all atoms
c     -------------------
      if (ga_nodeid() .eq. 0) 
     >   write(6,8) "comp","anal-g","num-g","error","step","ep","em" 
 8       format(1x,'checkgrad ',A5,4A16,2x,2A12)

      i=0
      do ia=1,ntot
        do k=1,3
          i=i+1
          c0=dbl_mb(i_c+i-1)
          gref=dbl_mb(i_g+i-1)
          
          step = step0
10        continue
          dbl_mb(i_c+i-1)=c0+step
          call mm_set_solute_coord(ntot,
     >                     int_mb(i_ai),
     >                     dbl_mb(i_c))

          call md_sp_qmmm()
          if (.not. qmmm_energy_gradient(rtdb,.false.))
     $    call errquit(pname//'qmmm_energy_gradient failed'
     $                  ,0, GEOM_ERR)
          call qmmm_energy_rtdb_push(rtdb)
          if (.not. rtdb_get(rtdb,'qmmm:energy', mt_dbl, 1, ep))
     $     call errquit('driver: could not get energy',0, RTDB_ERR)
          call qmmm_print_energy(rtdb)

c          if (abs(ep-e0) .lt. 1e-6) then
c             write(6,*) ' Increasing the step ', ep,e0,ep-e0, step
c             step = step*10.0d0
c             goto 10
c          else if (abs(ep-e0) .gt. 1e-2) then
c             write(6,*) ' Decreasing the step ', ep,e0,ep-e0, step
c             step = step/3.0
c             if(step.le.0.00000001d0) goto 12
c             goto 10
c          end if
12        continue

          dbl_mb(i_c+i-1)=c0-step
          call mm_set_solute_coord(ntot,
     >                     int_mb(i_ai),
     >                     dbl_mb(i_c))

          call md_sp_qmmm()
          if (.not. qmmm_energy_gradient(rtdb,.false.))
     $    call errquit(pname//'qmmm_energy_gradient failed'
     $                  ,0, GEOM_ERR)
          call qmmm_energy_rtdb_push(rtdb)
          if (.not. rtdb_get(rtdb,'qmmm:energy', mt_dbl, 1, em))
     $     call errquit('driver: could not get energy',0, RTDB_ERR)
          call qmmm_print_energy(rtdb)


          grad = (ep - em) / (2.0d0*step)
c
          if (ga_nodeid() .eq. 0) 
     >      write(6,7) int_mb(i_ai+ia-1), 
     >         gref, grad,abs(gref-grad), step, ep,em
 7        format(1x,'checkgrad ',i5,2f16.8,e16.8,f16.8,2x,2f16.8)
          dbl_mb(i_c+i-1)=c0
        end do
      end do


      end

      subroutine qmmm_check_forces0(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "mm_utils.fh"
c     
      integer rtdb
c
      integer ntot
      integer i_c,h_c
      integer i_c0,h_c0
      integer i_ai,h_ai
      integer i_g,h_g
      double precision c0,gref
      integer ia,i,k
      double precision e0, ep, em, grad, hess, step
      double precision step0
      integer i1,k1
c
      logical qmmm_energy_gradient
      external qmmm_energy_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname

      pname = "qmmm_check_forces"
c
       call md_sp()
       if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,e0))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          e0 = e0/cau2kj
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking forces'
         write(6,*)
         call util_flush(6)
      end if
c
c     get total number of atoms
c     -------------------------
      call mm_get_solute_tot_na_gen(ntot,mm_link) 
           if (ga_nodeid() .eq. 0) then
         write(6,*)
         write(6,*) '@ number of atoms',ntot
         write(6,*)
         call util_flush(6)
      end if

c
      if(.not.ma_alloc_get(MT_INT, ntot, 'qmmm ref index array',
     &      h_ai,i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      ntot, MA_ERR)
      call ifill(ntot,-1,int_mb(i_ai),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm ref grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm ref coord array',
     &      h_c,i_c) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_c),1)

      call mm_get_solute_ind_gen(ntot,
     >                         mm_link,
     >                         int_mb(i_ai))

      call mm_get_solute_coord_gen(ntot,
     >                         mm_link,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_c))

      call mm_get_solute_force_gen(ntot,
     >                         mm_link,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_g))
 
      if(.not.rtdb_get(rtdb,'qmmm:step',mt_dbl,
     >                 1,step0))
     >   step0=0.01d0

c      do i=1,3*ntot
c        dbl_mb(i_g+i-1)=dbl_mb(i_g+i-1)*cau2kj/cau2nm
c      end do
c
c     loop over all atoms
c     -------------------
      if (ga_nodeid() .eq. 0) 
     >   write(6,8) "comp","anal-g","num-g","error","step","ep","em" 
 8       format(1x,'checkgrad ',A5,4A16,2x,2A12)

      i=0
      do ia=1,ntot
        do k=1,3
          i=i+1
          c0=dbl_mb(i_c+i-1)
          gref=dbl_mb(i_g+i-1)
          
          step = step0
10        continue
          dbl_mb(i_c+i-1)=c0+step
          call mm_set_solute_coord_gen(ntot,
     >                       mm_link,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))
          call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,ep))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          ep = ep/cau2kj

          if (abs(ep-e0) .lt. 1e-6) then
             write(6,*) ' Increasing the step ', ep,e0,ep-e0, step
             step = step*10.0d0
             goto 10
          else if (abs(ep-e0) .gt. 1e-2) then
             write(6,*) ' Decreasing the step ', ep,e0,ep-e0, step
             step = step/3.0
             if(step.le.0.00000001d0) goto 12
             goto 10
          end if
12        continue

          dbl_mb(i_c+i-1)=c0-step
          call mm_set_solute_coord_gen(ntot,
     >                       mm_link,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))
c
            call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,em))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          em = em/cau2kj
c
          grad = (ep - em) / (2.0d0*step)
c
          if (ga_nodeid() .eq. 0) 
     >      write(6,7) int_mb(i_ai+ia-1), 
     >         gref, grad,abs(gref-grad), step, ep,em
 7        format(1x,'checkgrad ',i5,2f16.8,e16.8,f16.8,2x,2f16.8)
          dbl_mb(i_c+i-1)=c0
        end do
      end do


      end

c $Id$
