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
c     
      integer rtdb
c
      integer nqm
      integer nbq
      integer i_c
      integer i_cbq
      double precision c0,gref
      integer i,ia,k
      double precision e0, ep, em, grad, hess, step
      integer i_ai
      integer i_act
      integer i_g,h_g
c
      logical task_qmmm_gradient
      external task_qmmm_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname

      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking gradient for quantum atoms'
         write(6,*)
         call util_flush(6)
      end if
c
      if(.not.task_qmmm_gradient(rtdb)) 
     >  call errquit('qmmm_interface: failed ',0,0)
      if (.not. rtdb_get(rtdb,'qmmm:energy',mt_dbl,1,e0))
     $   call errquit('qmmm: failed put energy', 0, RTDB_ERR)
      e0=e0*cau2kj


      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking gradient for quantum atoms'
         write(6,*)
         call util_flush(6)
      end if
c
c     get coordinates of quantum atoms
c     -------------------------------
      nqm = qmmm_get_nqm()
      nbq = qmmm_get_nbq()
      i_ai = qmmm_get_i_ai()
      i_act = qmmm_get_i_act()
      i_c = qmmm_get_i_c()
c
      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm ref grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*nqm, MA_ERR)
      call dfill(3*nqm,0,dbl_mb(i_g),1)
c
      call mm_get_solute_force_quant(nqm,
     &                     .false.,
     &                     int_mb(i_ai),
     &                     dbl_mb(i_g))
c
      do i=1,3*nqm
        dbl_mb(i_g+i-1)=dbl_mb(i_g+i-1)*cau2kj/cau2nm
      end do
c
      write(*,*) "energy and gradients"
      do i=1,nqm
        write(*,*) int_mb(i_ai+i-1),(dbl_mb(i_g+(i-1)*3+k-1),k=1,3)
      end do
cc
c     loop over all active quantum atoms
c     ---------------------------------
      if (ga_nodeid() .eq. 0) 
     >   write(6,8) "comp","anal-g","num-g","error","hess","step","de" 
 8       format(1x,'checkgrad ',A5,4A16,2x,2A12)
      i=0
      do ia=1,nqm
        do k=1,3
          i=i+1
          c0=dbl_mb(i_c+i-1)
          gref=dbl_mb(i_g+i-1)
          step = 0.01d0
10        continue
          dbl_mb(i_c+i-1)=c0+step 
          call mm_set_solute_coord_quant(nqm,
     &                   .false.,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))

c          call qmmm_push_qm_coord_mm(rtdb)
          call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,ep))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
c          ep=ep/cau2kj
          if (abs(ep-e0) .lt. 1e-4) then
             write(6,*) ' Increasing the step ', ep-e0, step
             step = step*10.0d0
             goto 10
          else if (abs(ep-e0) .gt. 1e-2) then
             write(6,*) ' Decreasing the step ', ep-e0, step
             step = step/3.0
             goto 10
          end if

          dbl_mb(i_c+i-1)=c0-step 
          call mm_set_solute_coord_quant(nqm,
     &                   .false.,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))

c          call qmmm_push_qm_coord_mm(rtdb)
          call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,em))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
c          em=em/cau2kj
          grad = (ep - em) / (2.d0*step*cau2nm)

          if (ga_nodeid() .eq. 0) 
     >      write(6,7) i, gref, grad,abs(gref-grad), step, ep,em
 7        format(1x,'checkgrad ',i5,2f16.8,e16.8,f16.8,2x,3f12.6)
          dbl_mb(i_c+i-1)=c0
        end do
      end do


      end

      subroutine qmmm_check_forces_all(rtdb)
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
      integer i1,k1
c
      logical task_qmmm_gradient
      external task_qmmm_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname

      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking gradient for quantum atoms'
         write(6,*)
         call util_flush(6)
      end if
c
       call md_sp()
       if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,e0))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
c
      write(*,*) pname,"e0=",e0
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
      call mm_get_solute_tot_na_gen(ntot,mm_identity) 
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
     >                         mm_identity,
     >                         int_mb(i_ai))

      call mm_get_solute_coord_gen(ntot,
     >                         mm_identity,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_c))

      call mm_get_solute_force_gen(ntot,
     >                         mm_identity,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_g))
 
      do i=1,3*ntot
        dbl_mb(i_g+i-1)=dbl_mb(i_g+i-1)*cau2kj/cau2nm
      end do
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
          
          step = 0.000001d0
10        continue
          dbl_mb(i_c+i-1)=c0+step*cnm2au
          call mm_set_solute_coord_gen(ntot,
     >                       mm_identity,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))
          call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,ep))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          ep = ep-e0

          if (abs(ep) .lt. 1e-6) then
             write(6,*) ' Increasing the step ', ep,e0,ep-e0, step
             step = step*10.0d0
             goto 10
          else if (abs(ep) .gt. 1e-2) then
             write(6,*) ' Decreasing the step ', ep,e0,ep-e0, step
             step = step/3.0
             if(step.le.0.00000001d0) goto 12
             goto 10
          end if
12        continue

          dbl_mb(i_c+i-1)=c0-step*cnm2au
          call mm_set_solute_coord_gen(ntot,
     >                       mm_identity,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))
c
            call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,em))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          em=em-e0
c
          grad = (ep - em) / (2.0d0*step)
c
          if (ga_nodeid() .eq. 0) 
     >      write(6,7) int_mb(i_ai+ia-1), 
     >         gref, grad,abs(gref-grad), step*cau2nm, ep,em
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
      logical task_qmmm_gradient
      external task_qmmm_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname

      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking gradient for quantum atoms'
         write(6,*)
         call util_flush(6)
      end if
c
       call md_sp()
       if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,e0))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          e0 = e0/cau2kj
c
      write(*,*) pname,"e0=",e0
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
      call mm_get_solute_tot_na_gen(ntot,mm_quant) 
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
     >                         mm_quant,
     >                         int_mb(i_ai))

      call mm_get_solute_coord_gen(ntot,
     >                         mm_quant,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_c))

      call mm_get_solute_force_gen(ntot,
     >                         mm_quant,
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
     >                       mm_quant,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))
          write(*,*) "debug coord link0",(dbl_mb(i_c+i1-1),i1=1,3*ntot)
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
     >                       mm_quant,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))
          write(*,*) "debug coord link0",(dbl_mb(i_c+i1-1),i1=1,3*ntot)
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

      subroutine qmmm_check_forces_link(rtdb)
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
      logical task_qmmm_gradient
      external task_qmmm_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname

      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking gradient for quantum atoms'
         write(6,*)
         call util_flush(6)
      end if
c
       call md_sp()
       if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,e0))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          e0 = e0/cau2kj
c
      write(*,*) pname,"e0=",e0
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
          write(*,*) "debug coord link0",(dbl_mb(i_c+i1-1),i1=1,3*ntot)
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
          write(*,*) "debug coord link0",(dbl_mb(i_c+i1-1),i1=1,3*ntot)
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

      subroutine qmmm_check_forces3(rtdb)
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
      integer i_ai,h_ai
      integer i_g,h_g
      double precision c0,gref
      integer ia,i,k
      double precision e0, ep, em, grad, hess, step
c
      logical task_qmmm_gradient
      external task_qmmm_gradient
      logical task_qmmm_energy
      external task_qmmm_energy
      character*32 pname

      pname = "qmmm_check_forces"
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) '@ Checking gradient for quantum atoms'
         write(6,*)
         call util_flush(6)
      end if
c
      if(.not.task_qmmm_gradient(rtdb)) 
     >  call errquit('qmmm_interface: failed ',0,0)
c
      if (.not. rtdb_get(rtdb,'qmmm:energy',mt_dbl,1,e0))
     $   call errquit('qmmm: failed put energy', 0, RTDB_ERR)
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
      call mm_get_solute_tot_na_gen(ntot,mm_identity) 
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
     >                         mm_identity,
     >                         int_mb(i_ai))

      call mm_get_solute_coord_gen(ntot,
     >                         mm_identity,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_c))

      call mm_get_solute_force_gen(ntot,
     >                         mm_identity,
     >                         int_mb(i_ai),
     >                         dbl_mb(i_g))
 

c     loop over all atoms
c     -------------------
      if (ga_nodeid() .eq. 0) 
     >   write(6,8) "comp","anal-g","num-g","error","hess","step","de" 
 8       format(1x,'checkgrad ',A5,4A16,2x,2A12)

      i=0
      do ia=1,ntot
        do k=1,3
          i=i+1
          c0=dbl_mb(i_c+i-1)
          gref=dbl_mb(i_g+i-1)
          step = 0.01d0
10        continue
          dbl_mb(i_c+i-1)=c0+step 
          call mm_set_solute_coord_gen(ntot,
     >                       mm_identity,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))

          call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,ep))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          ep=ep/cau2kj


          if (abs(ep-e0) .lt. 1e-4) then
             write(6,*) ' Increasing the step ',ep,e0,ep-e0,step
             step = step*10.0d0
             goto 10
          else if (abs(ep-e0) .gt. 1e-2) then
             write(6,*) ' Decreasing the step ',ep,e0,ep-e0,step
             step = step/3.0
             goto 10
          end if

          dbl_mb(i_c+i-1)=c0-step 
          call mm_set_solute_coord_gen(ntot,
     >                       mm_identity,
     >                       int_mb(i_ai),
     >                       dbl_mb(i_c))

           call md_sp()
          if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,em))
     $      call errquit('qmmm: failed put energy', 0, RTDB_ERR)
          em=em/cau2kj

          grad = (ep - em) / (2.d0*step)

          if (ga_nodeid() .eq. 0) 
     >      write(6,7) int_mb(i_ai+ia-1), 
     >         gref, grad,abs(gref-grad), step, ep,em
 7        format(1x,'checkgrad ',i5,2f16.8,e16.8,f16.8,2x,3f12.6)
          dbl_mb(i_c+i-1)=c0
        end do
      end do


      end

