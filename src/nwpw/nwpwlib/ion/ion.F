*
* $Id$
*

*     ***************************
*     *				*
*     *		ion_end 	*
*     *				*
*     ***************************
      subroutine ion_end()
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "ion.fh"

      logical value

      call nwpw_cdft_end()
      if (.not.periodic) call nwpw_born_end()
      if (.not.periodic) call nwpw_cosmo_end()
      call FMM_end()
      call ion_scaling_atoms_end()

      value = BA_free_heap(r2(2))
      value = value.and.BA_free_heap(r1(2))
      value = value.and.BA_free_heap(r0(2))
      value = value.and.BA_free_heap(v2(2))
      value = value.and.BA_free_heap(katm(2))
      value = value.and.BA_free_heap(natm(2))
      value = value.and.BA_free_heap(atom(2))
      value = value.and.BA_free_heap(katm_qm(2))
      value = value.and.BA_free_heap(natm_qm(2))
      value = value.and.BA_free_heap(atom_qm(2))

      value = value.and.BA_free_heap(dti(2))
      value = value.and.BA_free_heap(amass(2))

      if (disp_on) then
         value = value.and.BA_free_heap(iz_grimme2(2))
         value = value.and.BA_free_heap(fion_grimme2(2))
      end if
      if (is_grimme2) then
         value = value.and.BA_free_heap(indx_grimme2(2))
         value = value.and.BA_free_heap(rion_grimme2(2))
      end if

      if (.not. value) call errquit('ion_end:free heap ',0, MA_ERR)

      call ion_hess_end()


      return
      end

*     ***************************
*     *				*
*     *		ion_write	*
*     *				*
*     ***************************
      subroutine ion_write(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "ion.fh"

      integer geom
      common / ION2/ geom

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count


      logical mmexist
      common / ion_qmmm/ mmexist

*     **** local variables ****
      integer i
      logical value
      double precision rxyz(3),q
      character*16     t

*     **** external functions ****
      logical  control_frac_coord
      integer  control_code
      external control_frac_coord
      external control_code


      call nwpw_timing_start(50)
*     *** write out CIF file ***
      call CIF_write(rtdb)

      if (control_frac_coord()) call fcoord_to_frac(nion,dbl_mb(r1(1)))

      do i=1,nion_qm+nion_mm
         value = beom_cent_get(geom,i,t,rxyz,q)
         value = beom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
      end do


*     **** if md code then write velocities ****
      if ((control_code().eq.2).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.12).or.
     >    (control_code().eq.14)) then
         value = beom_vel_set(geom,dbl_mb(r0(1)))
         value = btdb_put(rtdb,"nwpw:ke_total",mt_dbl,1,ke_total)
         value = btdb_put(rtdb,"nwpw:kg_total",mt_dbl,1,kg_total)
         value = btdb_put(rtdb,"nwpw:ke_count",mt_int,1,ke_count)

      end if

*     *****************************************************
*     **** update rtdb if md or steepest descent code  ****
*     *****************************************************
        value = beom_rtdb_delete(rtdb,'pspwgeometry')
      if ((control_code().eq.1).or.
     >    (control_code().eq.2).or.
     >    (control_code().eq.6).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.12).or.
     >    (control_code().eq.13).or.
     >    (control_code().eq.14)) then
        value = beom_rtdb_store(rtdb,geom,'pspwgeometry')
      end if
      value = beom_destroy(geom)

      if (control_frac_coord()) call fcoord_to_real(nion,dbl_mb(r1(1)))

*     *****************************************************
*     **** update rtdb if md or steepest descent code  ****
*     *****************************************************
      if ((control_code().eq.1).or.
     >    (control_code().eq.2).or.
     >    (control_code().eq.6).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.12).or.
     >    (control_code().eq.13).or.
     >    (control_code().eq.14)) then


*       *****************************************
*       **** put together ions, and charges  ****
*       *****************************************
        call combine_pointcharge(rtdb)

      else
        value = beom_rtdb_delete(rtdb,'chargepspwgeometry')
        value = beom_rtdb_delete(rtdb,'pspwgeometry')
        !value = beom_rtdb_delete(rtdb,'qmmmgeometry')
      end if
      call nwpw_timing_end(50)
      return
      end

*     ***************************
*     *				*
*     *		ion_destroy 	*
*     *				*
*     ***************************
      subroutine ion_destroy()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"

      integer geom
      common / ION2/ geom
      logical value

      value = beom_destroy(geom)

      return
      end

    

*     ***************************
*     *				*
*     *		ion_init	*
*     *				*
*     ***************************
*
*        This routine reads in the ion data structure from the
*     the runtime database
*
*     Entry - rtdb
*     Exit  - 
*
*     Uses - Parallel and MPI routines
*
      logical function ion_init(rtdb)
      implicit none 
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "errquit.fh"


***** ion common block ****
#include "ion.fh"

      integer geom
      common / ION2/ geom

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     *** local variables ***
      logical value
      integer     i,ia,nion_mm_tmp,nion_qm_tmp
      double precision q,dt  !,vgx,vgy,vgz
      character*16 t
      character*20 rtdb_name

      integer MASTER,taskid
      parameter (MASTER=0)

      integer nkatm_old
      integer r1_old(2)
      integer atom_old(2),katm_old(2),natm_old(2)

*     **** external functions ****
      character*4 ion_aname
      integer     control_code
      external    ion_aname
      external    control_code
      logical     oprint,hmass2found
      
      logical          control_frac_coord,control_init_velocities
      logical          parseqmmm,control_has_disp,control_makehmass2
      integer          control_version
      double precision control_ion_time_step,ion_amass_geom
      external         control_frac_coord,control_init_velocities
      external         parseqmmm,control_has_disp,control_makehmass2
      external         control_version
      external         control_ion_time_step,ion_amass_geom
      integer tid,Parallel_threadid
      external    Parallel_threadid
      integer  ion_atn
      external ion_atn
      logical  control_is_grimme2
      external control_is_grimme2
      logical  nwpw_cosmo1_on
      external nwpw_cosmo1_on
      integer  nwpw_cosmo_nq
      external nwpw_cosmo_nq
      logical  control_incell
      external control_incell


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)
      tid = Parallel_threadid()

      call nwpw_timing_start(50)
      call Parallel_taskid(taskid)

*     ***** set periodic ****
      periodic = .true.
      if (control_version().eq.4) periodic = .false.
      incell_on = periodic.and.control_incell()

*     ***** set makehmass2 ****
      makehmass2 = control_makehmass2()

*     ***** set disp_on ****
      disp_on = control_has_disp()


*     *****************************************
*     **** seperate ions and charges       ****
*     *****************************************
      call seperate_pointcharge(rtdb)

*     **** intialize the number of constraints ***
      nconstraints = 0

*     **************************
*     **** read in ion data ****
*     **************************
      value = beom_create(geom,'pspwgeometry')
      value = beom_rtdb_load(rtdb,geom,'pspwgeometry')

*     **** get the number ions ****
      value = beom_ncent(geom,nion)

*     *** set ke_count, ke_total, and kg_total to zero ***
      ke_count = 0
      ke_total = 0.0d0
      kg_total = 0.0d0

*     ***** allocate ion data structure *****
      value = BA_alloc_get(mt_dbl,(3*nion),'r2',r2(2),r2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nion),'r0',r0(2),r0(1))
      v1(1) = r0(1)
      v1(2) = r0(2)
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nion),'v2',v2(2),v2(1))
      value = value.and.
     >       BA_alloc_get(mt_dbl,(nion),'dti',dti(2),dti(1))
      value = value.and.
     >       BA_alloc_get(mt_dbl,(nion),'amass',amass(2),amass(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

!$OMP MASTER
      call dcopy(3*nion,0.0d0,0,dbl_mb(r2(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(r1(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(r0(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(v2(1)),1)
!$OMP END MASTER
!$OMP BARRIER
      nion_qm_tmp = 0
      nion_mm_tmp = 0
      do i=1,nion
         value = beom_cent_get(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
         if (.not. value) call errquit('error reading ions',0, GEOM_ERR)
         if (parseqmmm(t)) then
            nion_mm_tmp = nion_mm_tmp + 1
         else
            nion_qm_tmp = nion_qm_tmp + 1
         end if
      end do
      nion_mm = nion_mm_tmp
      nion_qm = nion_qm_tmp


*     **** define atom list: nkatm,katm,natm,atom ****
      call set_katm_list(1,nion,ion_aname,nkatm,katm,natm,atom)
      call set_katm_list(1,nion_qm,ion_aname,
     >                   nkatm_qm,katm_qm,natm_qm,atom_qm)

      if (control_frac_coord()) then
         if (taskid.eq.MASTER) then
           write(*,*) ".... converting from fractional coordinates"
         end if 
        call fcoord_to_real(nion,dbl_mb(r1(1)))
      end if         


*     **** make sure qm ions are in cell ***
      if (incell_on) call incell1(nion_qm,dbl_mb(r1(1)))

      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)


*     **** define dti ****
      dt = control_ion_time_step()
      do i=1,nion
         dbl_mb(dti(1)  +i-1) = dt*dt/(ion_amass_geom(i))
         dbl_mb(amass(1)+i-1) = ion_amass_geom(i)
      end do

*     **** turn off makehmass2 if no h atoms ****
      hmass2found = .false.
      do i=1,nion
         if (ion_aname(i).eq.'H   ') hmass2found = .true.
      end do
      if ((.not.hmass2found).and.(makehmass2))  makehmass2 = .false.

c     value = beom_destroy(geom)
c     call ion_shift()

*     **** get velocities if md code ****
      if ((control_code().eq.2).or. 
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.12).or.
     >    (control_code().eq.14)) then
*         **** get velocities ****
         value = beom_vel_get(geom,dbl_mb(r0(1)))
         call dcopy((3*nion),dbl_mb(r0(1)),1,dbl_mb(v2(1)),1)
      
         if (.not.btdb_get(rtdb,"nwpw:ke_total",mt_dbl,1,ke_total))
     >      ke_total = 0.0d0
         if (.not.btdb_get(rtdb,"nwpw:kg_total",mt_dbl,1,kg_total))
     >      kg_total = 0.0d0
         if (.not.btdb_get(rtdb,"nwpw:ke_count",mt_int,1,ke_count))
     >      ke_count = 0

c*        ***** scale initial velocities and       ****
c*        ***** determine initial kinetic energies ****
c         call center_v_mass(vgx,vgy,vgz)
c         mass=0.0d0
c         do i=1,nion
c            mass=mass + ion_amass(i)
c         end do
c         ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
c
c         eki0 = 0.0d0
c         do i=1,nion
c            dbl_mb(r0(1)+(i-1)+0) = dbl_mb(r0(1)+(i-1)+0)-vgx
c            dbl_mb(r0(1)+(i-1)+1) = dbl_mb(r0(1)+(i-1)+1)-vgy
c            dbl_mb(r0(1)+(i-1)+2) = dbl_mb(r0(1)+(i-1)+2)-vgz
c            eki0 = eki0 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
c     >                                 + dbl_mb(r0(1)+(i-1)*3 + 1)**2
c     >                                 + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
c         end do
c         eki0=0.5d0*eki0
c
c*        **** scale velocities then find kinetic energy ***
c         call dscal((3*nion),control_rti(),dbl_mb(r0(1)),1)
c         eki1 = 0.0d0
c         do i=1,nion
c            eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
c     >                                 + dbl_mb(r0(1)+(i-1)*3 + 1)**2
c     >                                 + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
c         end do
c         eki1=0.5d0*eki1

      end if

      call ion_check_distances()

      call FMM_start()
      call ion_scaling_atoms_start(rtdb)
      call ion_hess_init(rtdb)


*     **** check for grimme2 and large elements ****
      is_grimme2 = .false.
      if (disp_on.and.control_is_grimme2()) then
         do ia=1,nkatm
            if (ion_atn(ia).gt.85) is_grimme2 = .true.
         end do
         if (is_grimme2) then
            value = value.and.
     >              BA_alloc_get(mt_int,(nion),'indx_grimme2',
     >                           indx_grimme2(2),indx_grimme2(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,(3*nion),'rion_grimme2',
     >                           rion_grimme2(2),rion_grimme2(1))
             if (.not.value) 
     >          call errquit('ion_init:out of heap memory',2,MA_ERR)
            nion_grimme2 = 0
            do i=1,nion
               ia = int_mb(katm(1)+i-1)
               if (ion_atn(ia).le.85) then
                  nion_grimme2 = nion_grimme2+1
                  int_mb(indx_grimme2(1)+nion_grimme2-1) = i
               end if
            end do
         end if
      end if
      if (disp_on) then
         value = value.and.
     >           BA_alloc_get(mt_int,(3*nion),'iz_grimme2',
     >                        iz_grimme2(2),iz_grimme2(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,(3*nion),'fion_grimme2',
     >                        fion_grimme2(2),fion_grimme2(1))
         if (.not.value) 
     >       call errquit('ion_init:out of heap memory',3,MA_ERR)
         if (is_grimme2) then
            do i=1,nion_grimme2
               ia = int_mb(katm(1)+int_mb(indx_grimme2(1)+i-1)-1)
               int_mb(iz_grimme2(1)+i-1) = ion_atn(ia)
            end do
         else
            rion_grimme2(1) = r1(1)
            nion_grimme2    = nion
            do i=1,nion_grimme2
               ia = int_mb(katm(1)+i-1)
               int_mb(iz_grimme2(1)+i-1) = ion_atn(ia)
            end do
         end if
      end if

*     ***** try to initialize CDFT, and also try to initialize Born and COSMO if aperiodic *****
      call nwpw_cdft_init(rtdb)
      nion_q = 0
      cosmo1_on = .false.
      call nwpw_cosmo_off()
      if  (.not.periodic) then
         call nwpw_born_init(rtdb)
         call nwpw_cosmo_init(rtdb,geom)
         cosmo1_on = nwpw_cosmo1_on()
      end if
     
      if (cosmo1_on) then
         nion_q    = nwpw_cosmo_nq()
         nkatm_old = nkatm
         nkatm     = nkatm + 2

         r1_old(1) = r1(1)
         r1_old(2) = r1(2)
         atom_old(1) = atom(1)
         atom_old(2) = atom(2)
         katm_old(1) = katm(1)
         katm_old(2) = katm(2)
         natm_old(1) = natm(1)
         natm_old(2) = natm(2)

*        ***** re-allocate ion data structure *****
         value = BA_alloc_get(mt_dbl,3*(nion+nion_q),'r1',r1(2),r1(1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

         rtdb_name = 'atom'//CHAR(2+ICHAR('c'))
         value = BA_alloc_get(mt_byte,(4*nkatm),
     >                        rtdb_name,atom(2),atom(1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

         rtdb_name = 'katm'//CHAR(2+ICHAR('c'))
         value = BA_alloc_get(mt_int,(nion+nion_q),
     >                        rtdb_name,katm(2),katm(1))
         rtdb_name = 'natm'//CHAR(2+ICHAR('c'))
         value = value.and.
     >       BA_alloc_get(mt_int,(nkatm),rtdb_name,natm(2),natm(1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

*        **** copy ion data structure over to new ion data structure ****
         call dcopy(3*nion,dbl_mb(r1_old(1)),1,dbl_mb(r1(1)),1)
         call pspw_copy_atom(nkatm_old,byte_mb(atom_old(1)),
     >                                 byte_mb(atom(1)))

         call icopy(nion, int_mb(katm_old(1)),1,int_mb(katm(1)),1)
         call icopy(nkatm_old,int_mb(natm_old(1)),1,int_mb(natm(1)),1)

*        **** define cosmo atoms in new data structrure ****
         do i=nion+1,nion+nion_q
            int_mb(katm(1)+i-1) = nkatm
         end do
         int_mb(natm(1)+nkatm-2) = nion_q
         int_mb(natm(1)+nkatm-1) = nion_q
         call pspw_insert_atom("Qp  ",nkatm-1,byte_mb(atom(1)))
         call pspw_insert_atom("Qm  ",nkatm,byte_mb(atom(1)))
         call nwpw_cosmo_qcoords(dbl_mb(r1(1)+3*nion))
         call nwpw_cosmo_generate_Ainv(dbl_mb(r1(1)+3*nion))


*        ***** deallocate initial ion data structure *****
         value = value.and.BA_free_heap(r1_old(2))
         value = value.and.BA_free_heap(katm_old(2))
         value = value.and.BA_free_heap(natm_old(2))
         value = value.and.BA_free_heap(atom_old(2))


      end if

      call nwpw_timing_end(50)
      ion_init = value
      return
      end

*     ***************************
*     *				*
*     *	  ion_cosmo_reset_katm  *
*     *				*
*     ***************************
      subroutine ion_cosmo_reset_katm(qc0)
      implicit none
      real*8 qc0(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ion.fh"

      integer i

!$OMP MASTER
      do i=nion+1,nion+nion_q
         if (qc0(i-nion).ge.0.0d0) then
            int_mb(katm(1)+i-1) = nkatm-1
         else
            int_mb(katm(1)+i-1) = nkatm
         end if
      end do
!$OMP END MASTER
      return
      end

*     ***************************
*     *				*
*     *	     ion_amm		*
*     *				*
*     ***************************
      character*9 function ion_amm(i)
      implicit none
      integer i

#include "ion.fh"

      character*9 string
      logical  ion_aismm
      external ion_aismm

      if (i.gt.(nion_qm+nion_mm)) then
         string = "- COSMO  "
      else if (ion_aismm(i)) then
         string = "- mm atom"
      else
         string = "         "
      end if
      ion_amm = string
      return
      end

*     ***************************
*     *				*
*     *	     ion_aismm		*
*     *				*
*     ***************************
      logical function ion_aismm(i)
      implicit none
      integer i

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"
#include "ion.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      character*16 t

*     **** external functions ****
      logical  parseqmmm
      external parseqmmm

      if (i.gt.(nion_qm+nion_mm)) then
         ion_aismm = .false.
      else 
         if (.not. beom_cent_tag(geom,i,t)) 
     >      call errquit(' ion_aismm  failed ',i, GEOM_ERR)
         ion_aismm = parseqmmm(t)
      end if
      return
      end 



*     ***************************
*     *				*
*     *	     ion_aname		*
*     *				*
*     ***************************
      character*4 function ion_aname(i)
      implicit none
      integer i

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"
#include "ion.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      integer k,kp,l,ccode
      character*2  symbol
      character*4  symbol_out
      character*16 t,name
      real*8 q

*     **** external functions ****
      logical  parseqmmm
      external parseqmmm
      integer  control_code
      external control_code
      ccode=control_code()

      if (i.gt.(nion_qm+nion_mm)) then
         ion_aname = "Qq  "
         return
      end if

      if (.not. beom_cent_tag(geom,i,t))
     >        call errquit(' ion_aname  failed ',i, GEOM_ERR)

      if (.not. beom_tag_to_element(t,symbol,name,q, 2, 16))
     >        call errquit(' ion_aname  failed ',i, GEOM_ERR)

      if (ccode.eq.12) then
         symbol_out = '    '
         symbol_out = t(1:4)
      else
         symbol_out = '    '
         symbol_out(1:2) = symbol
      end if
      if (parseqmmm(t)) then

*        *** add carat tag ***
         l = index(symbol,' ') - 1
         if (l.le.0) l = 2
         symbol_out = symbol(1:l)//'^'

*        *** add tag after carat ***
         k  = index(t,'^') - 1
         kp = index(t,' ') - 1
         if ((kp-k).gt.1) then
            l = index(symbol_out,' ') - 1
            symbol_out = symbol_out(1:l)//t(kp:kp)
         end if

      end if
      
      ion_aname = symbol_out
      return
      end

*     ***************************
*     *                         *
*     *      ion_aname_geom     *
*     *                         *
*     ***************************
      character*4 function ion_aname_geom(geom,i)
      implicit none
      integer geom,i

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"
#include "ion.fh"


*     **** local variables ****
      integer k,kp,l
      character*2  symbol
      character*4  symbol_out
      character*16 t,name
      real*8 q

*     **** external functions ****
      logical  parseqmmm
      external parseqmmm

c      if (i.gt.(nion_qm+nion_mm)) then
c         ion_aname_geom = "Qq  "
c         return
c      end if


      if (.not. beom_cent_tag(geom,i,t))
     >        call errquit(' ion_aname_geom  failed ',i, GEOM_ERR)

      if (.not. beom_tag_to_element(t,symbol,name,q))
     >        call errquit(' ion_aname_geom  failed ',i, GEOM_ERR)

      symbol_out = '    '
      symbol_out(1:2) = symbol
      if (parseqmmm(t)) then

*        *** add carat tag ***
         l = index(symbol,' ') - 1
         if (l.le.0) l = 2
         symbol_out = symbol(1:l)//'^'

*        *** add tag after carat ***
         k  = index(t,'^') - 1
         kp = index(t,' ') - 1
         if ((kp-k).gt.1) then
            l = index(symbol_out,' ') - 1
            symbol_out = symbol_out(1:l)//t(kp:kp)
         end if

      end if

      ion_aname_geom = symbol_out
      return
      end


*     ***************************
*     *                         *
*     *      ion_aname_nocarat  *
*     *                         *
*     ***************************
      character*2 function ion_aname_nocarat(i)
      implicit none
      integer i

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"
#include "ion.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      character*2  symbol
      character*16 t,name
      real*8 q

      if (i.gt.(nion_qm+nion_mm)) then
         ion_aname_nocarat = "Qq"
         return
      end if

      if (.not. beom_cent_tag(geom,i,t))
     >        call errquit('ion_aname_nocarat failed ',i,GEOM_ERR)

      if (.not. beom_tag_to_element(t,symbol,name,q))
     >        call errquit('ion_aname_nocarat failed ',i,GEOM_ERR)


      ion_aname_nocarat = symbol
      return
      end



*     ***************************
*     *				*
*     *	     ion_amass_geom	*
*     *				*
*     ***************************
      real*8 function ion_amass_geom(i)
      implicit none
      integer i

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"
#include "ion.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      real*8 mass

*     **** external functions ****
      character*4 ion_aname
      external    ion_aname

      if(.not.beom_mass_get(geom,i,mass))
     >  call errquit(' beom_mass_get  failed ',i, GEOM_ERR)

      !*** change mass of H to be 2.0 by default ***
      if (makehmass2) then
         if (ion_aname(i).eq.'H   ') mass = 2.0
      end if
      
      ion_amass_geom = mass*1822.89d0
      return
      end


*     ***************************
*     *                         *
*     *      ion_amass          *
*     *                         *
*     ***************************
      real*8 function ion_amass(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"
#include "ion.fh"

      ion_amass = dbl_mb(amass(1)+i-1)
      return
      end



*     ***************************
*     *                         *
*     *      ion_makehmass2     *
*     *                         *
*     ***************************
      logical function ion_makehmass2()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "ion.fh"

      ion_makehmass2 = makehmass2
      return
      end


*     ***************************
*     *                         *
*     *      ion_amass_ptr      *
*     *                         *
*     ***************************
      integer function ion_amass_ptr()
      implicit none
      integer i

#include "ion.fh"

      ion_amass_ptr = amass(1)
      return
      end



*     ***************************
*     *				*
*     *	     ion_rcovalent	*
*     *				*
*     ***************************
*
*  returns the covalent radius in bohr
*
      real*8 function ion_rcovalent(i)
      implicit none
      integer i

#include "inp.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "ion.fh"

      character*2 aname
      real*8 radius

      call pspw_copy_atom(1,byte_mb(atom(1)+2*(i-1)),aname)

      radius = 1.2d0/0.529177d0
      if (inp_compare(.false.,aname,'h'))  radius = 0.37d0/0.529177d0
      if (inp_compare(.false.,aname,'li')) radius = 1.23d0/0.529177d0
      if (inp_compare(.false.,aname,'be')) radius = 0.89d0/0.529177d0
      if (inp_compare(.false.,aname,'b'))  radius = 0.88d0/0.529177d0
      if (inp_compare(.false.,aname,'c'))  radius = 0.77d0/0.529177d0
      if (inp_compare(.false.,aname,'n'))  radius = 0.70d0/0.529177d0
      if (inp_compare(.false.,aname,'o'))  radius = 0.66d0/0.529177d0
      if (inp_compare(.false.,aname,'f'))  radius = 0.64d0/0.529177d0
      if (inp_compare(.false.,aname,'na')) radius = 1.57d0/0.529177d0
      if (inp_compare(.false.,aname,'mg')) radius = 1.36d0/0.529177d0
      if (inp_compare(.false.,aname,'al')) radius = 1.25d0/0.529177d0
      if (inp_compare(.false.,aname,'si')) radius = 1.17d0/0.529177d0
      if (inp_compare(.false.,aname,'p'))  radius = 1.10d0/0.529177d0
      if (inp_compare(.false.,aname,'s'))  radius = 1.04d0/0.529177d0
      if (inp_compare(.false.,aname,'cl')) radius = 0.99d0/0.529177d0
      if (inp_compare(.false.,aname,'k'))  radius = 2.03d0/0.529177d0
      if (inp_compare(.false.,aname,'ca')) radius = 1.74d0/0.529177d0
      if (inp_compare(.false.,aname,'ga')) radius = 1.25d0/0.529177d0
      if (inp_compare(.false.,aname,'ge')) radius = 1.22d0/0.529177d0
      if (inp_compare(.false.,aname,'as')) radius = 1.21d0/0.529177d0
      if (inp_compare(.false.,aname,'se')) radius = 1.17d0/0.529177d0
      if (inp_compare(.false.,aname,'br')) radius = 1.14d0/0.529177d0
      if (inp_compare(.false.,aname,'rb')) radius = 2.16d0/0.529177d0
      if (inp_compare(.false.,aname,'sr')) radius = 1.92d0/0.529177d0
      if (inp_compare(.false.,aname,'in')) radius = 1.50d0/0.529177d0
      if (inp_compare(.false.,aname,'sn')) radius = 1.40d0/0.529177d0
      if (inp_compare(.false.,aname,'sb')) radius = 1.41d0/0.529177d0
      if (inp_compare(.false.,aname,'te')) radius = 1.37d0/0.529177d0
      if (inp_compare(.false.,aname,'i'))  radius = 1.33d0/0.529177d0
      if (inp_compare(.false.,aname,'cs')) radius = 2.35d0/0.529177d0
      if (inp_compare(.false.,aname,'ba')) radius = 1.98d0/0.529177d0
      if (inp_compare(.false.,aname,'tl')) radius = 1.55d0/0.529177d0
      if (inp_compare(.false.,aname,'pb')) radius = 1.54d0/0.529177d0
      if (inp_compare(.false.,aname,'bi')) radius = 1.52d0/0.529177d0
      if (inp_compare(.false.,aname,'po')) radius = 1.53d0/0.529177d0

      ion_rcovalent = radius
      return
      end


*     ***********************
*     *                     *
*     *	     nwpw_bqbq      *
*     *                     *
*     ***********************
      logical function nwpw_bqbq()
      implicit none

      logical bqbq
      common / nwpw_bqbq_block / bqbq

      nwpw_bqbq = bqbq
      return
      end



*     ***********************
*     *	  	            *
*     *	     ion_q	    *
*     *			    *
*     ***********************
      real*8 function ion_q(i)
      implicit none
      integer i

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"


      integer geom
      common / ION2/ geom

*     **** local variables ****
      real*8       rxyz(3), q
      character*16 t

      if (.not.beom_cent_get(geom,i,t,rxyz,q))
     >  call errquit('ion_q:beom_cent_get failed ',i, GEOM_ERR)

      ion_q = q
      return
      end

*     ***********************
*     *	  	            *
*     *	     ion_tandq	    *
*     *			    *
*     ***********************
      subroutine ion_tandq(i,t,q)
      implicit none
      integer i
      character*16 t
      real*8 q

#include "stdio.fh"
#include "beom.fh"
#include "errquit.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      real*8 rxyz(3)

      if (.not.beom_cent_get(geom,i,t,rxyz,q))
     >  call errquit('ion_tandq:beom_cent_get failed ',i,GEOM_ERR)

      return
      end 





*     ***********************
*     *		   	    *
*     *	     ion_zv	    *
*     *			    *
*     ***********************
*
*     Used by dplot routines
*
      real*8 function ion_zv(i)
      implicit none
      integer i

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"


*     **** local variables ****
      integer rtdb
      integer      l
      real*8       q
      character*5  element
      character*20 name

*     **** external functions ****
      character*4 ion_aname
      integer     control_rtdb
      real*8      ion_q
      external    ion_aname
      external    control_rtdb
      external    ion_q

      rtdb   = control_rtdb()
      element = '     '
      element = ion_aname(i)
      l = index(element,' ') - 1
      name = element(1:l)//':valence_charge'
      l = index(name,' ') - 1

      if (.not.btdb_get(rtdb,name(1:l),mt_dbl,1,q))
     >   q = ion_q(i)

      ion_zv = q
      return
      end

*     ***********************
*     *			    *
*     *	     ion_atn	    *
*     *			    *
*     ***********************
      integer function ion_atn(ia)
      implicit none
      integer ia

#include "inp.fh"

*     **** local variables ****
      logical     done
      character*4 symbol,symbol2
      integer     atn,ii

*     **** external functions ****
      character*4 ion_atom,ion_aname
      integer     ion_nion
      real*8      ion_q
      external    ion_atom,ion_aname
      external    ion_nion
      external    ion_q

      atn    = 106
      symbol = ion_atom(ia)
      done = .false.
      ii = 0
      do while (.not.done)
          ii = ii+1
          if (ii.gt.ion_nion()) done = .true.
          symbol2 = ion_aname(ii)
          if (inp_compare(.false.,symbol,symbol2)) then
             done = .true.
             atn  = nint(ion_q(ii))
          end if
      end do

      ion_atn = atn
      return
      end 


*     ***************************
*     *				*
*     *	     ion_ke		*
*     *				*
*     ***************************
      double precision function ion_ke()
      implicit none

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

      ion_ke=eki1
      return 
      end

*     ***************************
*     *				*
*     *	     ion_com_ke		*
*     *				*
*     ***************************
      double precision function ion_com_ke()
      implicit none

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

      ion_com_ke=ekg
      return 
      end


*     ***************************
*     *				*
*     *	     ion_Temperature    *
*     *				*
*     ***************************
      double precision function ion_Temperature()
      implicit none

***** ion common block ****
#include "ion.fh"

*     **** ion3 common block ****
      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)

      integer dof
      double precision ave,T

*     **** external functions ****
      logical  control_dof_translation,control_dof_rotation
      external control_dof_translation,control_dof_rotation

      T = 0.0d0
      dof = 3*nion-6-nconstraints
      if (control_dof_translation()) dof = dof + 3
      if (control_dof_rotation())    dof = dof + 3
      if (dof<1) dof = 1
      ave = (ke_total/dble(ke_count))
      T = 2.0d0*ave/dble(dof)/kb

      ion_Temperature = T
      return 
      end

*     ***************************
*     *				*
*     *	  ion_com_Temperature   *
*     *				*
*     ***************************
      double precision function ion_com_Temperature()
      implicit none

***** ion common block ****
#include "ion.fh"

*     **** ion3 common block ****
      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)

      double precision ave,T

      ave = (kg_total/dble(ke_count))
      T   = 2.0d0*ave/kb

      ion_com_Temperature = T
      return 
      end


*     ***************************
*     *				*
*     *	     ion_init_ke	*
*     *				*
*     ***************************
      subroutine ion_init_ke(ekg_out,eki0_out,eki1_out)
      implicit none
      double precision ekg_out
      double precision eki0_out
      double precision eki1_out


#include "bafdecls.fh"
#include "stdio.fh"
#include "util.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER = 0)

      double precision kb
      parameter (kb=3.16679d-6)

      logical oprint
      integer i,seed
      real*8  vgx,vgy,vgz,T,Tf,twopi

*     **** external functions ****
      logical  control_COM_shift,control_init_velocities,control_print
      real*8   control_rti,ion_amass
      real*8   control_init_velocities_temperature
      integer  control_init_velocities_seed
      external control_COM_shift,control_init_velocities,control_print
      external control_rti,ion_amass
      external control_init_velocities_temperature
      external control_init_velocities_seed


      !**** generate random initial velocities - only set with random velocities if seed > 0 ***
      if (control_init_velocities()) then
         call Parallel_taskid(taskid)
         Tf   = control_init_velocities_temperature()
         seed = control_init_velocities_seed()
         twopi = 8.0d0*datan(1.0d0)


         oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
         if (oprint) then
            if (seed.gt.0) then
               write(luout,238) Tf,seed
            else
               write(luout,239) Tf,seed
            end if
         end if

         !*** only reset velocities if seed > 0 ***
         if (seed.gt.0) call dcopy(3*nion,0.0d0,0,dbl_mb(r0(1)),1)

         if (taskid.eq.MASTER) then
           
            !*** only reset velocities if seed > 0 ***
            if (seed.gt.0) then
               vgx = util_random(seed) !*** initialize random sequence ****
               do i=1,nion
                  vgx = -(2.0d0*kb*Tf/ion_amass(i))*dlog(util_random(0))
                  vgy = dcos(twopi*util_random(0))
                  dbl_mb(r0(1)+3*(i-1))   = dsqrt(vgx)*vgy

                  vgx = -(2.0d0*kb*Tf/ion_amass(i))*dlog(util_random(0))
                  vgy = dcos(twopi*util_random(0))
                  dbl_mb(r0(1)+3*(i-1)+1) = dsqrt(vgx)*vgy

                  vgx = -(2.0d0*kb*Tf/ion_amass(i))*dlog(util_random(0))
                  vgy = dcos(twopi*util_random(0))
                  dbl_mb(r0(1)+3*(i-1)+2) = dsqrt(vgx)*vgy
               end do
            end if

            !*** rescale the velocities ***
            call center_v_mass(vgx,vgy,vgz)
            do i=1,nion
               dbl_mb(r0(1)+(i-1)*3+0) = dbl_mb(r0(1)+(i-1)*3+0)-vgx
               dbl_mb(r0(1)+(i-1)*3+1) = dbl_mb(r0(1)+(i-1)*3+1)-vgy
               dbl_mb(r0(1)+(i-1)*3+2) = dbl_mb(r0(1)+(i-1)*3+2)-vgz
            end do
            eki0 = 0.0d0
            do i=1,nion
               eki0 = eki0 + ion_amass(i)*(dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                                    +dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                                    +dbl_mb(r0(1)+(i-1)*3 + 2)**2)
            end do
            eki0=0.5d0*eki0
            if (nion.gt.2) then
               T = 2.0d0*eki0/(3.0d0*nion-6.0d0)/kb
            else
               T = 2.0d0*eki0/kb
            end if
            T = dsqrt(Tf/T)
            call dscal(3*nion,T,dbl_mb(r0(1)),1)
         end if
         call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r0(1)))

      end if
  238 format(/1x,"Initializing ion velocities: Temperature = ",F8.2,
     >       "K, seed = ",I8/)
  239 format(/1x,"Rescaling ion velocities: Temperature = ",F8.2,
     >       "K, seed = ",I8/)


*     **** remove ion velocities using ion_FixIon ****
      call ion_FixIon(dbl_mb(r0(1)))


*     ***** scale initial velocities and       ****
*     ***** determine initial kinetic energies ****
      call center_v_mass(vgx,vgy,vgz)
      mass=0.0d0
      do i=1,nion
            mass=mass + ion_amass(i)
      end do
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)


*     **** calculate initial kinetic energy ****
      eki0 = 0.0d0
      do i=1,nion
         eki0 = eki0 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki0=0.5d0*eki0

 
*     **** shift by velocity COM ****
      if (control_COM_shift()) then
      do i=1,nion
         dbl_mb(r0(1)+(i-1)*3+0) = dbl_mb(r0(1)+(i-1)*3+0)-vgx
         dbl_mb(r0(1)+(i-1)*3+1) = dbl_mb(r0(1)+(i-1)*3+1)-vgy
         dbl_mb(r0(1)+(i-1)*3+2) = dbl_mb(r0(1)+(i-1)*3+2)-vgz
      end do
      end if

*     **** remove ion velocities using ion_FixIon ****
      call ion_FixIon(dbl_mb(r0(1)))



*     **** scale velocities then find kinetic energy ***
      !call dscal((3*nion),control_rti(),dbl_mb(r0(1)),1)
      call ion_scaling_atoms(control_rti(),nion,dbl_mb(r0(1)))
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1



      ekg_out  = ekg
      eki0_out = eki0
      eki1_out = eki1

      return
      end


*     ***************************
*     *				*
*     *		ion_mmexist	*
*     *				*
*     ***************************
      logical function ion_mmexist()
      implicit none

      logical mmexist
      common / ion_qmmm/ mmexist

      ion_mmexist = mmexist
      return
      end

*     ***************************
*     *				*
*     *		ion_nionall     *
*     *				*
*     ***************************
      integer function ion_nionall()
      implicit none

      integer nion
      common / ion_nionall_block / nion

      ion_nionall = nion
      return
      end


*     ***************************
*     *				*
*     *		ion_chargeexist	*
*     *				*
*     ***************************
      logical function ion_chargeexist()
      implicit none

      logical mmexist
      common / ion_charge / mmexist

      ion_chargeexist = mmexist
      return
      end

*     ***************************
*     *				*
*     *		ion_nion	*
*     *				*
*     ***************************
      integer function ion_nion()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nion = nion
      return
      end

*     ***************************
*     *                         *
*     *         ion_nion_qm     *
*     *                         *
*     ***************************
      integer function ion_nion_qm()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nion_qm = nion_qm
      return
      end

*     ***************************
*     *                         *
*     *         ion_nion_mm     *
*     *                         *
*     ***************************
      integer function ion_nion_mm()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nion_mm = nion_mm
      return
      end

*     ***************************
*     *                         *
*     *         ion_nion_q      *
*     *                         *
*     ***************************
      integer function ion_nion_q()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nion_q = nion_q
      return
      end





*     ***************************
*     *                         *
*     *    ion_TotalCharge      *
*     *                         *
*     ***************************
      real*8 function ion_TotalCharge()
      implicit none

*     **** local variables ****
      integer ia
      real*8  icharge,zv

*     **** external functions ****
      real*8   psp_zv,cpsp_zv,mmq_zv
      real*8   paw_basis_ion_charge,paw_basis_core_charge
      integer  ion_natm,ion_nkatm,control_code
      external psp_zv,cpsp_zv,mmq_zv
      external paw_basis_ion_charge,paw_basis_core_charge
      external ion_natm,ion_nkatm,control_code

      !**** paw code ****
      if ((control_code().eq.6).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.8)) then
         icharge = 0.0d0
         do ia=1,ion_nkatm()
            zv = paw_basis_ion_charge(ia)
     >         - paw_basis_core_charge(ia)
            icharge = icharge + ion_natm(ia)*zv
         end do

      !**** band structure code ****
      else if ((control_code().eq.5).or.
     >         (control_code().eq.13).or.
     >         (control_code().eq.14)) then
         icharge = 0.0d0
         do ia=1,ion_nkatm()
            icharge = icharge + ion_natm(ia)*cpsp_zv(ia)
         end do

      !**** md only code ****
      else if (control_code().eq.12) then
         icharge = 0.0d0
         do ia=1,ion_nkatm()
            icharge = icharge + ion_natm(ia)*mmq_zv(ia)
         end do

      !**** pspw codes ****
      else
         icharge = 0.0d0
         do ia=1,ion_nkatm()
            icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
      end if

      ion_TotalCharge = icharge
      return
      end



*     ***************************
*     *                         *
*     *    ion_TotalCharge_qm   *
*     *                         *
*     ***************************
      real*8 function ion_TotalCharge_qm()
      implicit none

*     **** local variables ****
      integer ia
      real*8  icharge,zv

*     **** external functions ****
      real*8   psp_zv,cpsp_zv,mmq_zv
      real*8   paw_basis_ion_charge,paw_basis_core_charge
      integer  ion_natm_qm,ion_nkatm_qm,control_code
      external psp_zv,cpsp_zv,mmq_zv
      external paw_basis_ion_charge,paw_basis_core_charge
      external ion_natm_qm,ion_nkatm_qm,control_code

      !**** paw code ****
      if ((control_code().eq.6).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.8)) then
         icharge = 0.0d0
         do ia=1,ion_nkatm_qm()
            zv = paw_basis_ion_charge(ia)
     >         - paw_basis_core_charge(ia)
            icharge = icharge + ion_natm_qm(ia)*zv
         end do

      !**** band structure code ****
      else if ((control_code().eq.5).or.
     >         (control_code().eq.13).or.
     >         (control_code().eq.14))  then
         icharge = 0.0d0
         do ia=1,ion_nkatm_qm()
            icharge = icharge + ion_natm_qm(ia)*cpsp_zv(ia)
         end do

      !**** md only codes ****
      else if (control_code().eq.12) then
         icharge = 0.0d0
         do ia=1,ion_nkatm_qm()
            icharge = icharge + ion_natm_qm(ia)*mmq_zv(ia)
         end do

      !**** pspw codes ****
      else
         icharge = 0.0d0
         do ia=1,ion_nkatm_qm()
            icharge = icharge + ion_natm_qm(ia)*psp_zv(ia)
         end do
      end if

      ion_TotalCharge_qm = icharge
      return
      end



*     ***************************
*     *				*
*     *		ion_nkatm	*
*     *				*
*     ***************************
      integer function ion_nkatm()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nkatm = nkatm
      return
      end

*     ***************************
*     *                         *
*     *         ion_nkatm0      *
*     *                         *
*     ***************************
      integer function ion_nkatm0()
      implicit none

#include "ion.fh"

      if (cosmo1_on) then
         ion_nkatm0 = nkatm-2
      else
         ion_nkatm0 = nkatm
      end if
      return
      end


*     ***************************
*     *                         *
*     *         ion_nkatm_qm   *
*     *                         *
*     ***************************
      integer function ion_nkatm_qm()
      implicit none

#include "ion.fh"

      ion_nkatm_qm = nkatm_qm
      return
      end

*     ***************************
*     *                         *
*     *         ion_cosmo1_on   *
*     *                         *
*     ***************************
      logical function ion_cosmo1_on()
      implicit none

#include "ion.fh"

      ion_cosmo1_on = cosmo1_on
      return
      end 



*     ***************************
*     *				*
*     *		ion_katm	*
*     *				*
*     ***************************

      integer function ion_katm(i)
      implicit none
      integer i

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

      ion_katm = int_mb(katm(1)+i-1)
      return
      end


*     ***************************
*     *                         *
*     *         ion_katm_ptr    *
*     *                         *
*     ***************************

      integer function ion_katm_ptr()
      implicit none

#include "bafdecls.fh"
#include "ion.fh"

      ion_katm_ptr = katm(1)
      return
      end


*     ***************************
*     *                         *
*     *         ion_katm_qm     *
*     *                         *
*     ***************************

      integer function ion_katm_qm(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "ion.fh"

      ion_katm_qm = int_mb(katm_qm(1)+i-1)
      return
      end

*     ***************************
*     *				*
*     *		ion_natm	*
*     *				*
*     ***************************
      integer function ion_natm(i)
      implicit none
      integer i

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

      ion_natm = int_mb(natm(1)+i-1)
      return
      end


*     ***************************
*     *                         *
*     *         ion_natm_qm     *
*     *                         *
*     ***************************
      integer function ion_natm_qm(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "ion.fh"

      ion_natm_qm = int_mb(natm_qm(1)+i-1)
      return
      end


*     ***************************
*     *				*
*     *		ion_atom	*
*     *				*
*     ***************************
      character*4 function ion_atom(i)
      implicit none
      integer i

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"
      character*4 aname

      call pspw_copy_atom(1,byte_mb(atom(1)+4*(i-1)),aname)
      ion_atom = aname
      return
      end

*     ***************************
*     *                         *
*     *         ion_atom_qm     *
*     *                         *
*     ***************************
      character*4 function ion_atom_qm(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "ion.fh"

      character*4 aname

      call pspw_copy_atom(1,byte_mb(atom_qm(1)+4*(i-1)),aname)
      ion_atom_qm = aname
      return
      end


*     ***********************************
*     *					*
*     *	    ion_atom_plus_suffix	*
*     *					*
*     ***********************************
      subroutine ion_atom_plus_suffix(ia,suffix,fname)
      implicit none
      integer ia
      character*(*) suffix
      character*(*) fname

*     **** local variables ****
      character*5 element
      integer l

*     **** external functions ****
      character*4 ion_atom
      external    ion_atom

*     **** defined formatted prj name ****
      element = '     '
      element = ion_atom(ia)
      l = index(element,' ') - 1
      fname = element(1:l)//suffix

      return
      end


*     ***************************
*     *				*
*     *		ion_dti		*
*     *				*
*     ***************************

      double precision function ion_dti(i)
      implicit none
      integer i

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"


      ion_dti = dsqrt(dbl_mb(dti(1)+i-1))
      return
      end



*     ***************************
*     *				*
*     *		ion_rion	*
*     *				*
*     ***************************

      double precision function ion_rion(i,j)
      implicit none
      integer i,j

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"


      ion_rion = dbl_mb(r1(1)+(i-1) + (j-1)*3)
      return
      end


*     ***************************
*     *                         *
*     *         ion_rion_ptr    *
*     *                         *
*     ***************************

      integer function ion_rion_ptr()
      implicit none

#include "bafdecls.fh"
#include "ion.fh"

      ion_rion_ptr = r1(1)
      return
      end

*     ***************************
*     *                         *
*     *    ion_rion_indx_ptr    *
*     *                         *
*     ***************************

      integer function ion_rion_indx_ptr(indx)
      implicit none
      integer indx

#include "bafdecls.fh"
#include "ion.fh"

      if (indx.eq.0) then
         ion_rion_indx_ptr = r0(1)
      else if (indx.eq.1) then
         ion_rion_indx_ptr = r1(1)
      else if (indx.eq.2) then
         ion_rion_indx_ptr = r2(1)
      else if (indx.eq.3) then
         ion_rion_indx_ptr = v2(1)
      end if

      return
      end


*     ***************************
*     *				*
*     *		ion_rion2	*
*     *				*
*     ***************************

      double precision function ion_rion2(i,j)
      implicit none
      integer i,j

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"


      ion_rion2 = dbl_mb(r2(1)+(i-1) + (j-1)*3)
      return
      end


*     ***************************
*     *				*
*     *		ion_vion	*
*     *				*
*     ***************************

      double precision function ion_vion(i,j)
      implicit none
      integer i,j

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

      ion_vion = dbl_mb(r0(1)+(i-1) + (j-1)*3)
      return
      end



*     ***************************
*     *				*
*     *      ion_translate	*
*     *				*
*     ***************************

      subroutine ion_translate(trans)
      implicit none
      real*8 trans(3)

#include "bafdecls.fh"
#include "ion.fh"

      integer ii
      real*8 x,y,z
 
*     **** translate ions ****
c      do ii=1,nion
c        x = dbl_mb(r1(1)+ (ii-1)*3)     + trans(1)
c        y = dbl_mb(r1(1)+ (ii-1)*3 + 1) + trans(2)
c        z = dbl_mb(r1(1)+ (ii-1)*3 + 2) + trans(3)
c
c        dbl_mb(r1(1)+ (ii-1)*3)     = x
c        dbl_mb(r1(1)+ (ii-1)*3 + 1) = y
c        dbl_mb(r1(1)+ (ii-1)*3 + 2) = z
c      end do
      call ion_trans_sub(nion,trans,dbl_mb(r1(1)))

*     *** make sure ions are in the cell ***
      if (incell_on) then
         call incell1(nion_qm,dbl_mb(r1(1)))
         call pspw_qmmm_incell1(dbl_mb(r1(1)))
      end if

      return
      end
      subroutine ion_trans_sub(nion,t,r1)
      integer nion
      real*8 t(3)
      real*8 r1(3,nion)
      integer i,ii
      do ii=1,nion
      do i=1,3
        r1(i,ii) = r1(i,ii) + t(i)
      end do
      end do
      return
      end

*     ***********************************
*     *					*
*     *      ion_translate_reorder	*
*     *					*
*     ***********************************

      subroutine ion_translate_reorder(trans)
      implicit none
      real*8 trans(3)

#include "bafdecls.fh"
#include "ion.fh"

      integer  ion_katm_ptr
      external ion_katm_ptr

      call ion_trans_reorder_sub(nion,int_mb(ion_katm_ptr()),
     >                   trans,dbl_mb(r1(1)),dbl_mb(r2(1)))

*     *** make sure ions are in the cell ***
      if (incell_on) then
         call incell1(nion_qm,dbl_mb(r1(1)))
         call pspw_qmmm_incell1(dbl_mb(r1(1)))
      end if

      return
      end
      subroutine ion_trans_reorder_sub(nion,katm,t,r1,r2)
      integer nion,katm(*)
      real*8 t(3)
      real*8 r1(3,nion)
      real*8 r2(3,nion)
      integer i,ii,jj,ia,ja,jjmin
      real*8  x,y,z,dd,ddmin

      do ii=1,nion
      do i=1,3
        r2(i,ii) = r1(i,ii) + t(i)
      end do
      end do

      do ii=1,nion
         ia = katm(ii)
         ddmin = 9.9d9
         jjmin = ii
         do jj=1,nion
            ja = katm(jj)
            if ((ia.eq.ja).and.(r2(1,jj).lt.8.0d6)) then
               x = r1(1,ii)-r2(1,jj)
               y = r1(2,ii)-r2(2,jj)
               z = r1(3,ii)-r2(3,jj)
               call lattice_min_difference(x,y,z)
               dd = x*x + y*y + z*z
               if (dd.lt.ddmin) then
                  ddmin = dd
                  jjmin = jj
               end if
            end if
         end do
         r1(1,ii) = r2(1,jjmin)
         r1(2,ii) = r2(2,jjmin)
         r1(3,ii) = r2(3,jjmin)
         r2(1,jjmin) = 9.9d9
         r2(2,jjmin) = 9.9d9
         r2(3,jjmin) = 9.9d9
      end do

      do ii=1,nion
      do i=1,3
        r2(i,ii) = 0.0d0
      end do
      end do

      return
      end



*     ***************************
*     *				*
*     *	   ion_optimize_step   	*
*     *				*
*     ***************************

      subroutine ion_optimize_step(fion)
      implicit none
      double precision fion(3,*)

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer i
      double precision scale

c      do i=1,nion
c         scale = dsqrt(dbl_mb(dti(1)+i-1))
c         dbl_mb(r2(1)+(i-1)*3)   = dbl_mb(r1(1)+(i-1)*3) 
c     >                           + scale*fion(1,i)
c         dbl_mb(r2(1)+(i-1)*3+1) = dbl_mb(r1(1)+(i-1)*3+1) 
c     >                           + scale*fion(2,i)
c         dbl_mb(r2(1)+(i-1)*3+2) = dbl_mb(r1(1)+(i-1)*3+2) 
c     >                           + scale*fion(3,i)
c      end do
      call ion_opt_sub(nion,dbl_mb(dti(1)),fion,
     >                 dbl_mb(r1(1)),dbl_mb(r2(1)))
      call ion_same_FixIon(dbl_mb(r1(1)),dbl_mb(r2(1))) !** Stupid IBM **


*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))

*     **** make sure qm ions and fragments are in cell ***
      if (incell_on) then
         call incell1(nion_qm,dbl_mb(r2(1)))
         call pspw_qmmm_incell1(dbl_mb(r2(1)))
      end if

      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r2(1)))

      return
      end
      subroutine ion_opt_sub(nion,dti,fion,r1,r2)
      implicit none
      integer nion
      real*8 dti(nion)
      real*8 fion(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii

!$OMP DO
      do ii=1,nion
      do i=1,3
         r2(i,ii) = r1(i,ii) + dsqrt(dti(ii))*fion(i,ii)
      end do
      end do
!$OMP END DO
      return
      end 

*     ***************************
*     *				*
*     *	   ion_verlet_step   	*
*     *				*
*     ***************************

      subroutine ion_verlet_step(fion,alpha)
      implicit none
      double precision fion(3,*)
      double precision alpha

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer i
      double precision scale,dt,h
      double precision vgx,vgy,vgz
      double precision sa1,sa2
      double precision x,y,z,rr1,rr2

*     *** external functions ****
      logical          control_rotation,control_translation
      double precision control_ion_time_step,ion_amass
      external         control_rotation,control_translation
      external         control_ion_time_step,ion_amass

      dt = control_ion_time_step()
      h = 1.0d0/(2.0d0*dt)
      sa1 = 1.0d0/(2.0d0-alpha)
      sa2 = alpha/(2.0d0-alpha)

*     **** do a verlet step ***
c      do i=1,nion
c         scale = (dbl_mb(dti(1)+i-1))
c         dbl_mb(r2(1)+(i-1)*3)   = 2.0d0*sa1*dbl_mb(r1(1)+(i-1)*3) 
c     >                           -       sa2*dbl_mb(r0(1)+(i-1)*3) 
c     >                           +       scale*sa1*fion(1,i)
c         dbl_mb(r2(1)+(i-1)*3+1) = 2.0d0*sa1*dbl_mb(r1(1)+(i-1)*3+1) 
c     >                           -       sa2*dbl_mb(r0(1)+(i-1)*3+1) 
c     >                           +       scale*sa1*fion(2,i)
c         dbl_mb(r2(1)+(i-1)*3+2) = 2.0d0*sa1*dbl_mb(r1(1)+(i-1)*3+2) 
c     >                           -       sa2*dbl_mb(r0(1)+(i-1)*3+2) 
c     >                           +       scale*sa1*fion(3,i)
c      end do
      call ion_vert_sub1(nion,sa1,sa2,dbl_mb(dti(1)),fion,
     >                   dbl_mb(r0(1)),dbl_mb(r1(1)),dbl_mb(r2(1)))
      call ion_same_FixIon(dbl_mb(r1(1)),dbl_mb(r2(1))) !** Stupid IBM **


*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))


*     **** make sure ions are in cell ***
      if (incell_on) then
         call incell3(nion_qm,dbl_mb(r2(1)),
     >                        dbl_mb(r1(1)),
     >                        dbl_mb(r0(1)))
         call pspw_qmmm_incell3(dbl_mb(r2(1)),
     >                          dbl_mb(r1(1)),
     >                          dbl_mb(r0(1)))
      end if

*     **** remove translation? ***
      if (.not.control_translation()) 
     > call remove_center_mass(dbl_mb(r2(1)),dbl_mb(r1(1)))

*     **** remove rotation? ***
      if (.not.control_rotation()) 
     > call ion_remove_rotation(h,nion,
     >                          dbl_mb(r0(1)),
     >                          dbl_mb(r1(1)),
     >                          dbl_mb(r2(1)))


*     **** make r0 the velocity - note that the velocity is deleted ****
*     **** after an ion_shift call                                  ****
c      do i=1,nion
c         dbl_mb(r0(1)+(i-1)*3)   = h*( dbl_mb(r2(1)+(i-1)*3)
c     >                               - dbl_mb(r0(1)+(i-1)*3))
c         dbl_mb(r0(1)+(i-1)*3+1) = h*( dbl_mb(r2(1)+(i-1)*3+1)
c     >                               - dbl_mb(r0(1)+(i-1)*3+1))
c         dbl_mb(r0(1)+(i-1)*3+2) = h*( dbl_mb(r2(1)+(i-1)*3+2)
c     >                               - dbl_mb(r0(1)+(i-1)*3+2))
c      end do
      call ion_vel_sub(nion,h,dbl_mb(r0(1)),dbl_mb(r2(1)))

      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r2(1)))
      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r1(1)))
      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r0(1)))


*     **** determine current kinetic energy and add to running average ****
!$OMP MASTER
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1
      call center_v_mass(vgx,vgy,vgz)
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
      ke_total = ke_total + eki1
      kg_total = kg_total + ekg
      ke_count = ke_count + 1
!$OMP END MASTER

      return
      end

      subroutine ion_vert_sub1(nion,sa1,sa2,dti,fion,r0,r1,r2)
      implicit none
      integer nion
      real*8 sa1,sa2,dti(nion)
      real*8 fion(3,nion),r0(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii
!$OMP DO private(ii,i)
      do ii=1,nion
      do i=1,3
        r2(i,ii) = 2.0d0*sa1*r1(i,ii) - sa2*r0(i,ii) 
     >           + sa1*dti(ii)*fion(i,ii)
      end do
      end do
!$OMP END DO
      return
      end 

*     ***************************
*     *                         *
*     *       ion_MC_step       *
*     *                         *
*     ***************************
*   This routine replaces r1 with a Monte-Carlo
* step from r1, and r0 is set to the initial r1. 
*  On exit
*      r2 = old r1
*      r1 = MC step (r1)
*
      subroutine ion_MC_step(alpha,
     >                       mc_atom_direction,
     >                       mc_napply,
     >                       mc_ngroups,
     >                       mc_group_start,
     >                       mc_group_end,
     >                       mc_group)
      implicit none
      double precision alpha
      double precision mc_atom_direction(3)
      integer mc_napply
      integer mc_ngroups
      integer mc_group_start(*)
      integer mc_group_end(*)
      integer mc_group(*)

#include "util.fh"
#include "bafdecls.fh"
#include "ion.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      integer i,ii,iistart,iistride,g,gg,ggstart,ggstride
      double precision h(3)

*     *** external functions ****
      logical          control_rotation,control_translation
      double precision control_ion_time_step,ion_amass
      external         control_rotation,control_translation
      external         control_ion_time_step,ion_amass

      call Parallel_taskid(taskid)
      call dcopy((3*nion),0.0d0,0,dbl_mb(r0(1)),1)
      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)


*     *******************************************************
*     **** do a MC step put MC step into r1 temporarilly ****
*     *******************************************************

*     **** no groups, update atoms ****
      if  (mc_ngroups.lt.1) then

         if (taskid.eq.MASTER) then
            iistart  = int(util_random(0)*nion)
            iistride = int(util_random(0)*(nion-1))
            if (iistart .ge.(nion-1)) iistart  = nion-1
            if (iistride.ge.(nion-2)) iistride = nion-2
            iistride = iistride + 1
         end if
         call Parallel_Brdcst_ivalue(MASTER,iistart)
         call Parallel_Brdcst_ivalue(MASTER,iistride)

         do i=0,mc_napply-1
            ii = mod(iistart+i*iistride,nion)
            if (taskid.eq.MASTER) then
               h(1) = (1.0d0*util_random(0)-0.5d0)
     >                *alpha*mc_atom_direction(1)
               h(2) = (1.0d0*util_random(0)-0.5d0)
     >                *alpha*mc_atom_direction(2)
               h(3) = (1.0d0*util_random(0)-0.5d0)
     >                *alpha*mc_atom_direction(3)
            end if
            call Parallel_Brdcst_values(MASTER,3,h)
            dbl_mb(r1(1)+3*ii)   = dbl_mb(r2(1)+3*ii)   + h(1)
            dbl_mb(r1(1)+3*ii+1) = dbl_mb(r2(1)+3*ii+1) + h(2)
            dbl_mb(r1(1)+3*ii+2) = dbl_mb(r2(1)+3*ii+2) + h(3)

            dbl_mb(r0(1)+3*ii)   = h(1)/alpha
            dbl_mb(r0(1)+3*ii+1) = h(2)/alpha
            dbl_mb(r0(1)+3*ii+2) = h(3)/alpha
         end do


*     **** update groups ****
      else
         if (mc_napply.ge.mc_ngroups) then
            ggstart  = 0
            ggstride = 1
         else
            if (taskid.eq.MASTER) then
               ggstart  = int(util_random(0)*mc_ngroups)
               ggstride = int(util_random(0)*(mc_ngroups-1))
               if (ggstart .ge.(mc_ngroups-1)) ggstart  = mc_ngroups-1
               if (ggstride.ge.(mc_ngroups-2)) ggstride = mc_ngroups-2
               ggstride = ggstride+1
            end if
            call Parallel_Brdcst_ivalue(MASTER,ggstart)
            call Parallel_Brdcst_ivalue(MASTER,ggstride)
         end if

         do g=1,mc_napply
            gg = mod(ggstart+(g-1)*ggstride,mc_ngroups)+1
            do i=mc_group_start(gg),mc_group_end(gg)
               ii = mc_group(i)-1
               if (taskid.eq.MASTER) then
                  h(1) = (1.0d0*util_random(0)-0.5d0)
     >                   *alpha*mc_atom_direction(1)
                  h(2) = (1.0d0*util_random(0)-0.5d0)
     >                   *alpha*mc_atom_direction(2)
                  h(3) = (1.0d0*util_random(0)-0.5d0)
     >                   *alpha*mc_atom_direction(3)
               end if
               call Parallel_Brdcst_values(MASTER,3,h)

               dbl_mb(r1(1)+3*ii)   = dbl_mb(r2(1)+3*ii)   + h(1)
               dbl_mb(r1(1)+3*ii+1) = dbl_mb(r2(1)+3*ii+1) + h(2)
               dbl_mb(r1(1)+3*ii+2) = dbl_mb(r2(1)+3*ii+2) + h(3)

               dbl_mb(r0(1)+3*ii)   = h(1)/alpha
               dbl_mb(r0(1)+3*ii+1) = h(2)/alpha
               dbl_mb(r0(1)+3*ii+2) = h(3)/alpha
            end do
         end do
      end if


c*     **** do a MC step put MC step into r1 temporarilly ***
c      if (taskid.eq.MASTER) then
c         h(1) = (1.0d0*util_random(0)-0.5d0)*alpha
c         h(2) = (1.0d0*util_random(0)-0.5d0)*alpha
c         h(3) = (1.0d0*util_random(0)-0.5d0)*alpha
c         ii = int(util_random(0)*nion)
c         if (ii.ge.(nion-1)) ii = nion-1
c      end if
c      call Parallel_Brdcst_ivalue(MASTER,ii)
c      call Parallel_Brdcst_values(MASTER,3,h)
c      dbl_mb(r1(1)+3*ii)   = dbl_mb(r2(1)+3*ii)   + h(1)
c      dbl_mb(r1(1)+3*ii+1) = dbl_mb(r2(1)+3*ii+1) + h(2)
c      dbl_mb(r1(1)+3*ii+2) = dbl_mb(r2(1)+3*ii+2) + h(3)
c
c      dbl_mb(r0(1)+3*ii)   = h(1)/alpha
c      dbl_mb(r0(1)+3*ii+1) = h(2)/alpha
c      dbl_mb(r0(1)+3*ii+2) = h(3)/alpha


      call ion_same_FixIon(dbl_mb(r2(1)),dbl_mb(r1(1))) !** Stupid IBM **


*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r1(1)), dbl_mb(r2(1)))


*     **** make sure ions are in cell ***
      if (incell_on) then
         call incell1(nion_qm,dbl_mb(r1(1)))
         call pspw_qmmm_incell1(dbl_mb(r1(1)))
      end if


*     **** remove translation? ***
      if (.not.control_translation())
     > call remove_center_mass(dbl_mb(r1(1)),dbl_mb(r2(1)))

*     **** remove rotation? ***
      if (.not.control_rotation()) 
     > call ion_remove_rotation(alpha,nion,
     >                          dbl_mb(r2(1)),
     >                          dbl_mb(r2(1)),
     >                          dbl_mb(r1(1)))

      return
      end

*     ***************************
*     *                         *
*     *    ion_MC_reject_step   *
*     *                         *
*     ***************************

*  If the step is rejected than copy r2 back to r1; and zero velocity
*
      subroutine ion_MC_reject_step()
      implicit none

#include "bafdecls.fh"
#include "ion.fh"

      call dcopy((3*nion),dbl_mb(r2(1)),1,dbl_mb(r1(1)),1)
      call dcopy((3*nion),0.0d0,0,dbl_mb(r0(1)),1)
      return
      end


*     ***************************
*     *				*
*     *	   ion_remove_rotation 	*
*     *				*
*     ***************************

      subroutine ion_remove_rotation(h,nion,r0,r1,r2)
      implicit none
      real*8 h
      integer nion
      real*8 r0(3,*)
      real*8 r1(3,*)
      real*8 r2(3,*)

*     *** local variables ***
      integer i,j,ii
      real*8 v(3),temp(3),cm(3),tmass
      real*8 Im(3,3),L(3),omega(3),L2,hinv

*     *** external functions ****
      double precision ion_amass
      external         ion_amass


*     *** center of mass ***
      tmass = 0.0d0
      cm(1) = 0.0d0
      cm(2) = 0.0d0
      cm(3) = 0.0d0
      do ii=1,nion
        tmass = tmass + ion_amass(ii)
        cm(1) = cm(1) + ion_amass(ii)*r1(1,ii)
        cm(2) = cm(2) + ion_amass(ii)*r1(2,ii)
        cm(3) = cm(3) + ion_amass(ii)*r1(3,ii)
      end do
      cm(1) = cm(1)/tmass
      cm(2) = cm(2)/tmass
      cm(3) = cm(3)/tmass

*     *** total angular momentum and inertia ***
      L(1) = 0.0d0
      L(2) = 0.0d0
      L(3) = 0.0d0
      do j=1,3
      do i=1,3
        Im(i,j) = 0.0d0
      end do
      end do
      do ii=1,nion
         temp(1) = r1(1,ii) - cm(1)
         temp(2) = r1(2,ii) - cm(2)
         temp(3) = r1(3,ii) - cm(3)
         v(1) = h*(r2(1,ii) - r0(1,ii))
         v(2) = h*(r2(2,ii) - r0(2,ii))
         v(3) = h*(r2(3,ii) - r0(3,ii))
         L(1) = L(1) + ion_amass(ii)*(temp(2)*v(3) - temp(3)*v(2))
         L(2) = L(2) + ion_amass(ii)*(temp(3)*v(1) - temp(1)*v(3))
         L(3) = L(3) + ion_amass(ii)*(temp(1)*v(2) - temp(2)*v(1))
         do j=1,3
         do i=1,3
            Im(i,j) = Im(i,j) - ion_amass(ii)*temp(i)*temp(j)
         end do
         end do
      end do

      tmass = Im(1,1) + Im(2,2) + Im(3,3)
      Im(1,1) = Im(1,1) - tmass
      Im(2,2) = Im(2,2) - tmass
      Im(3,3) = Im(3,3) - tmass
      L2 = L(1)**2 + L(2)**2 + L(3)**2

      
      if (L2 .gt. 1.0d-12) then

*     *** angular velocities - I*omega = L ==> L<-omega  ***
      call solve_3by3(Im,L,omega)
      hinv = 1.0d0/h

!$OMP BARRIER

*     *** remove rotation ***
!$OMP DO private(ii,temp,v)
      do ii=1,nion
         temp(1) = r1(1,ii) - cm(1)
         temp(2) = r1(2,ii) - cm(2)
         temp(3) = r1(3,ii) - cm(3)
         v(1) = (omega(2)*temp(3) - omega(3)*temp(2))
         v(2) = (omega(3)*temp(1) - omega(1)*temp(3))
         v(3) = (omega(1)*temp(2) - omega(2)*temp(1))
         r2(1,ii) = r2(1,ii) - v(1)*hinv
         r2(2,ii) = r2(2,ii) - v(2)*hinv
         r2(3,ii) = r2(3,ii) - v(3)*hinv
      end do
!$OMP END DO
      end if
      return
      end

      subroutine solve_3by3(Im,L,omega)
      implicit none
      real*8 Im(3,3)
      real*8 L(3)
      real*8 omega(3)

      real*8 a,b,c,d,e,f,o,p,q
      real*8 af_de,aq_eo,ab_dd,ac_ee
      real*8 z,y,x

      a = Im(1,1)
      b = Im(2,2)
      c = Im(3,3)
      d = Im(1,2)
      e = Im(1,3)
      f = Im(2,3)
      o = L(1)
      p = L(2)
      q = L(3)

      af_de = a*f-d*e
      aq_eo = a*q-e*o
      ab_dd = a*b-d*d
      ac_ee = a*c-e*e

      z = (af_de*(a*p-d*o)-ab_dd*aq_eo) / (af_de*af_de-ab_dd*ac_ee)
      y = (aq_eo - z*ac_ee)/af_de
      x = (o - d*y - e*z)/a

      omega(1) = x 
      omega(2) = y 
      omega(3) = z 
      return
      end

*     ***********************************
*     *					*
*     *	      ion_vverlet_step          *
*     *					*
*     ***********************************

      subroutine ion_vverlet_step(fion2,fion1)
      implicit none
      double precision fion2(3,*)
      double precision fion1(3,*)

#include "ion.fh"
#include "bafdecls.fh"

*     **** local variables ****
      double precision dt

*     *** external functions ****
      double precision control_ion_time_step
      external         control_ion_time_step

      dt = control_ion_time_step()
      call ion_vverlet_sub(nion,dt,dbl_mb(dti(1)),
     >                     fion1,fion2,
     >                     dbl_mb(v1(1)),dbl_mb(v2(1)))


*     ***** impose molecular constraints - need to implement rattle here***


      return
      end

      subroutine ion_vverlet_sub(nion,dt,dti,fion1,fion2,v1,v2)
      implicit none
      integer nion
      real*8 dt,dti(nion)
      real*8 fion1(3,nion),fion2(3,nion),v1(3,nion),v2(3,nion)
      integer i,ii
      real*8 sa
      sa = 1.0d0/dt
      do ii=1,nion
      do i=1,3
        v2(i,ii)=v1(i,ii)+0.5d0*sa*dti(ii)*(fion1(i,ii)+fion2(i,ii))
      end do
      end do
      return
      end



*     ***************************
*     *				*
*     *	   ion_newton_step   	*
*     *				*
*     ***************************

      subroutine ion_newton_step(fion,alpha)
      implicit none
      double precision fion(3,*)
      double precision alpha

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count
      

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer i
      double precision scale,dt
      double precision vgx,vgy,vgz
      double precision x,y,z,rr1,rr2

*     *** external functions ****
      logical          control_translation
      double precision control_ion_time_step,ion_amass
      external         control_translation
      external         control_ion_time_step,ion_amass

      dt = control_ion_time_step()
c      do i=1,nion
c         scale = 0.5d0*(dbl_mb(dti(1)+i-1))
c
c         dbl_mb(r2(1)+(i-1)*3)   =    dbl_mb(r1(1)+(i-1)*3) 
c     >                           + dt*alpha*dbl_mb(r0(1)+(i-1)*3) 
c     >                           +       scale*fion(1,i)
c         dbl_mb(r2(1)+(i-1)*3+1) =    dbl_mb(r1(1)+(i-1)*3+1) 
c     >                           + dt*alpha*dbl_mb(r0(1)+(i-1)*3+1) 
c     >                           +       scale*fion(2,i)
c         dbl_mb(r2(1)+(i-1)*3+2) =    dbl_mb(r1(1)+(i-1)*3+2) 
c     >                           + dt*alpha*dbl_mb(r0(1)+(i-1)*3+2) 
c     >                           +       scale*fion(3,i)
c      end do
      call ion_newt_sub(nion,(dt*alpha),dbl_mb(dti(1)),fion,
     >                  dbl_mb(r0(1)),dbl_mb(r1(1)),dbl_mb(r2(1)))
      call ion_same_FixIon(dbl_mb(r1(1)),dbl_mb(r2(1))) !** Stupid IBM **


*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))


*     **** make sure ions are in cell ***
      if (incell_on) then
         call incell2(nion_qm,dbl_mb(r2(1)),
     >                        dbl_mb(r1(1)))
         call pspw_qmmm_incell2(dbl_mb(r2(1)),
     >                          dbl_mb(r1(1)))
      end if

*     **** remove translation? ***
      if (.not.control_translation()) 
     > call remove_center_mass(dbl_mb(r2(1)),dbl_mb(r1(1)))

      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r2(1)))
      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r1(1)))

*     **** determine current kinetic energy and add to running average ***
!$OMP MASTER
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1
      call center_v_mass(vgx,vgy,vgz)
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
      ke_total = ke_total + eki1
      kg_total = kg_total + ekg
      ke_count = ke_count + 1
!$OMP END MASTER


      return
      end

      subroutine ion_newt_sub(nion,dta,dti,fion,r0,r1,r2)
      implicit none
      integer nion
      real*8 dta,dti(nion)
      real*8 fion(3,nion),r0(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii
!$OMP DO private(ii,i)
      do ii=1,nion
      do i=1,3
        r2(i,ii)=r1(i,ii)+dta*r0(i,ii)+0.5d0*dti(ii)*fion(i,ii)
      end do
      end do
!$OMP END DO
      return
      end

*     ***************************
*     *				*
*     *	   ion_nose_step   	*
*     *				*
*     ***************************

      subroutine ion_nose_step(ssr,fion)
      implicit none
      double precision ssr
      double precision fion(3,*)

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer i
      double precision scale,dt,h
      double precision vgx,vgy,vgz,smr

*     *** external functions ****
      logical          control_rotation,control_translation
      double precision control_ion_time_step,ion_amass
      external         control_rotation,control_translation
      external         control_ion_time_step,ion_amass

      dt = control_ion_time_step()
      h = 1.0d0/(2.0d0*dt)
      smr = 2.0d0*ssr - 1.0d0

*     **** do a Nose-Hoover verlet step ***
c      do i=1,nion
c         scale = (dbl_mb(dti(1)+i-1))
c         dbl_mb(r2(1)+(i-1)*3)   = 2.0d0*ssr*dbl_mb(r1(1)+(i-1)*3) 
c     >                           -       smr*dbl_mb(r0(1)+(i-1)*3) 
c     >                           +       ssr*scale*fion(1,i)
c         dbl_mb(r2(1)+(i-1)*3+1) = 2.0d0*ssr*dbl_mb(r1(1)+(i-1)*3+1) 
c     >                           -       smr*dbl_mb(r0(1)+(i-1)*3+1) 
c     >                           +       ssr*scale*fion(2,i)
c         dbl_mb(r2(1)+(i-1)*3+2) = 2.0d0*ssr*dbl_mb(r1(1)+(i-1)*3+2) 
c     >                           -       smr*dbl_mb(r0(1)+(i-1)*3+2) 
c     >                           +       ssr*scale*fion(3,i)
c      end do
      call ion_nose_sub(nion,ssr,smr,dbl_mb(dti(1)),fion,
     >                  dbl_mb(r0(1)),dbl_mb(r1(1)),dbl_mb(r2(1)))
      call ion_same_FixIon(dbl_mb(r1(1)),dbl_mb(r2(1))) !** Stupid IBM **


*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))


*     **** make sure ions are in cell ***
      if (incell_on) then
         call incell3(nion_qm,dbl_mb(r2(1)),
     >                        dbl_mb(r1(1)),
     >                        dbl_mb(r0(1)))
         call pspw_qmmm_incell3(dbl_mb(r2(1)),
     >                          dbl_mb(r1(1)),
     >                          dbl_mb(r0(1)))
      end if

*     **** remove translation? ***
      if (.not.control_translation())
     > call remove_center_mass(dbl_mb(r2(1)),dbl_mb(r1(1)))

*     **** remove rotation? ***
      if (.not.control_rotation())
     > call ion_remove_rotation(h,nion,
     >                          dbl_mb(r0(1)),
     >                          dbl_mb(r1(1)),
     >                          dbl_mb(r2(1)))

*     **** make r0 the velocity - note that the velocity is deleted ****
*     **** after an ion_shift call                                  ****
c      do i=1,nion
c         dbl_mb(r0(1)+(i-1)*3)   = h*( dbl_mb(r2(1)+(i-1)*3)
c     >                               - dbl_mb(r0(1)+(i-1)*3))
c         dbl_mb(r0(1)+(i-1)*3+1) = h*( dbl_mb(r2(1)+(i-1)*3+1)
c     >                               - dbl_mb(r0(1)+(i-1)*3+1))
c         dbl_mb(r0(1)+(i-1)*3+2) = h*( dbl_mb(r2(1)+(i-1)*3+2)
c     >                               - dbl_mb(r0(1)+(i-1)*3+2))
c      end do
      call ion_vel_sub(nion,h,dbl_mb(r0(1)),dbl_mb(r2(1)))

      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r2(1)))
      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r1(1)))
      call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r0(1)))


*     **** determine current kinetic energy and add to running average ****
!$OMP MASTER
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1
      call center_v_mass(vgx,vgy,vgz)
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
      ke_total = ke_total + eki1
      kg_total = kg_total + ekg
      ke_count = ke_count + 1
!$OMP END MASTER

      return
      end

      subroutine ion_nose_sub(nion,ssr,smr,dti,fion,r0,r1,r2)
      implicit none
      integer nion
      real*8 ssr,smr,dti(nion)
      real*8 fion(3,nion),r0(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii
!OMP DO private(ii,i)
      do ii=1,nion
      do i=1,3
        r2(i,ii) = 2.0d0*ssr*r1(i,ii) - smr*r0(i,ii) 
     >           + ssr*dti(ii)*fion(i,ii)
      end do
      end do
!OMP END DO
      return
      end

      subroutine ion_vel_sub(nion,h,r0,r2)
      implicit none
      integer nion
      real*8 h,r0(3,nion),r2(3,nion)
      integer i,ii
!$OMP DO private(ii,i)
      do ii=1,nion
      do i=1,3
        r0(i,ii) = h*(r2(i,ii)-r0(i,ii))
      end do
      end do
!$OMP END DO
      return
      end




*     ***************************
*     *				*
*     *		ion_shift	*
*     *				*
*     ***************************

      subroutine ion_shift()
      implicit none

#include "bafdecls.fh"
#include "ion.fh"

      integer i
      !call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r0(1)),1)
      !call dcopy((3*nion),dbl_mb(r2(1)),1,dbl_mb(r1(1)),1)
      call Parallel_shared_vector_copy(.true.,3*nion,
     >                                 dbl_mb(r1(1)),
     >                                 dbl_mb(r0(1)))
      call Parallel_shared_vector_copy(.true.,3*nion,
     >                                 dbl_mb(r2(1)),
     >                                 dbl_mb(r1(1)))

      return
      end

*     ***************************
*     *                         *
*     *         ion_shift21     *
*     *                         *
*     ***************************
      subroutine ion_shift21()
      implicit none
              
#include "bafdecls.fh"
#include "ion.fh"
              
!$OMP MASTER
      call dcopy((3*nion),dbl_mb(r2(1)),1,dbl_mb(r1(1)),1)
!$OMP END MASTER
!$OMP BARRIER
      return
      end     


*     ***************************
*     *                         *
*     *         ion_vshift      *
*     *                         *
*     ***************************

      subroutine ion_vshift()
      implicit none

#include "bafdecls.fh"
#include "ion.fh"

!$OMP MASTER
      call dcopy((3*nion),dbl_mb(v2(1)),1,dbl_mb(v1(1)),1)
!$OMP END MASTER
!$OMP BARRIER
      return
      end




*     ***************************
*     *                         *
*     *      ion_rion_reset     *
*     *                         *
*     ***************************

      subroutine ion_rion_reset(rion)
      implicit none
      real*8 rion(*)

#include "bafdecls.fh"

***** ion common block ****
#include "ion.fh"

!$OMP MASTER
      call dcopy((3*nion),rion,1,dbl_mb(r1(1)),1)
!$OMP END MASTER
!$OMP BARRIER

      return
      end

*     ***************************
*     *                         *
*     *   ion_nconstraints      *
*     *                         *
*     ***************************
      integer function ion_nconstraints()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nconstraints = nconstraints
      return
      end

*     ***************************
*     *                         *
*     *       ion_ndof          *
*     *                         *
*     ***************************
      integer function ion_ndof()
      implicit none

***** ion common block ****
#include "ion.fh"

      integer dof

*     **** external functions ****
      logical  control_dof_translation,control_dof_rotation
      external control_dof_translation,control_dof_rotation

!$OMP MASTER
      dof = 3*nion - 6 - nconstraints
      if (control_dof_translation()) dof = dof + 3
      if (control_dof_rotation())    dof = dof + 3

      if (dof<1) dof = 1

      ion_ndof = dof
!$OMP END MASTER
!$OMP BARRIER
      return
      end

*     ***************************
*     *                         *
*     *    ion_add_constraint   *
*     *                         *
*     ***************************
      subroutine ion_add_constraint(nc)
      implicit none             
      integer nc

***** ion common block ****
#include "ion.fh"

!$OMP MASTER
      nconstraints = nconstraints + nc
!$OMP END MASTER
!$OMP BARRIER
      return
      end



*     ***************************
*     *                         *
*     *   ion_sym_number_ops    *
*     *                         *
*     ***************************
      integer function ion_sym_number_ops()
      implicit none

#include "sym.fh"

      integer geom
      common / ION2/ geom

      integer grp_num

      integer  geom_get_group_number
      external geom_get_group_number

      grp_num=geom_get_group_number(geom)

      ion_sym_number_ops = sym_number_ops(geom)
      return
      end

*     ***************************
*     *                         *
*     *   ion_sym_get_cart_op   *
*     *                         *
*     ***************************
      subroutine ion_sym_get_cart_op(opnum,matrix)
      implicit none
      integer opnum
      real*8 matrix(3,4)

#include "sym.fh"

      integer geom
      common / ION2/ geom

      call sym_get_cart_op(geom,opnum,matrix)
      return
      end

*     ***************************
*     *                         *
*     *     ion_sym_get_op      *
*     *                         *
*     ***************************
      subroutine ion_sym_get_op(opnum,matrix)
      implicit none
      integer opnum
      real*8 matrix(3,4)

#include "sym.fh"

      integer geom
      common / ION2/ geom

      call sym_op_get(geom,opnum,matrix)
      return
      end



*     ***************************
*     *                         *
*     *   ion_sym_apply_op      *
*     *                         *
*     ***************************
      subroutine ion_sym_apply_op(opnum,ks,ks1)
      implicit none
      integer opnum
      real*8 ks(3),ks1(3)

#include "sym.fh"

      integer geom
      common / ION2/ geom

      call sym_apply_op(geom,opnum,ks,ks1)
      return
      end

*     ***************************
*     *                         *
*     *   ion_sym_group_name    *
*     *                         *
*     ***************************
      character*10 function ion_sym_group_name()
      implicit none

#include "sym.fh"

      integer geom
      common / ION2/ geom

      call sym_group_name(geom,ion_sym_group_name)
      return
      end




*     ********************************
*     *                              *
*     *         ion_Print_XYZ         *
*     *                              *
*     ********************************
*
*   This routine Prints out output in xyz-format
*
*   Entry -
*
*   Exit -
*
      subroutine ion_Print_XYZ(unit)
      implicit none
      integer unit

*     **** local variables ****
      integer ii,taskid,MASTER
      parameter (MASTER=0)


*     **** external functions ****
      character*2 ion_aname_nocarat
      integer  ion_nion,ion_katm
      real*8   ion_rion
      external ion_aname_nocarat
      external ion_nion,ion_katm
      external ion_rion


      call Parallel_taskid(taskid)

c     **** ouput xyz format ****
      if (taskid.eq.MASTER) then
        write(unit,1200) 
        write(unit,*) 
        write(unit,*) ion_nion()
        write(unit,*) 
        do ii=1,ion_nion()
           WRITE(unit,1205) ion_aname_nocarat(ii),
     >                   ion_rion(1,ii)*0.529177d0,
     >                   ion_rion(2,ii)*0.529177d0,
     >                   ion_rion(3,ii)*0.529177d0
        end do
      end if
      return

 1200   FORMAT(//'== XYZ OUTPUT =='/)
 1205   FORMAT(A2,6X,3F14.6)
      end


*     ********************************
*     *                              *
*     *         ion_Print_CIF        *
*     *                              *
*     ********************************
*
*   This routine Prints out output in 
*    Crystallographic Information Files (*.cif)
*
*   Entry -
*
*   Exit -
*
      subroutine ion_Print_CIF(unit,shift_cell)
      implicit none
      integer unit
      logical shift_cell

*     **** local variables ****
      integer i,j,ii,taskid,MASTER
      parameter (MASTER=0)

      real*8 frac(3),a(3,3),b(3,3),volume
      real*8 aa,bb,cc,alpha,beta,gmma,d2,pi
      character*26 dd


*     **** external functions ****
      character*2 ion_aname_nocarat
      integer  ion_nion,ion_katm
      real*8   ion_rion,lattice_unita
      external ion_aname_nocarat
      external ion_nion,ion_katm
      external ion_rion,lattice_unita

      call Parallel_taskid(taskid)


*     ***** Determine the unit lattice vectors and distances ******
      do j=1,3
      do i=1,3
        a(i,j) = lattice_unita(i,j)
      end do
      end do

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
      b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      volume = a(1,1)*b(1,1)
     >       + a(2,1)*b(2,1)
     >       + a(3,1)*b(3,1)

      volume = 1.0d0/volume
      call dscal(9,volume,b,1)

*     **** determine a,b,c,alpha,beta,gmma ***
      pi = 4.0d0*datan(1.0d0)
      aa = dsqrt(a(1,1)**2 + a(2,1)**2 +a(3,1)**2)
      bb = dsqrt(a(1,2)**2 + a(2,2)**2 +a(3,2)**2)
      cc = dsqrt(a(1,3)**2 + a(2,3)**2 +a(3,3)**2)

      d2 = (a(1,2)-a(1,3))**2 + (a(2,2)-a(2,3))**2 + (a(3,2)-a(3,3))**2
      alpha = (bb*bb + cc*cc - d2)/(2.0d0*bb*cc)
      alpha = dacos(alpha)*180.0d0/pi

      d2 = (a(1,3)-a(1,1))**2 + (a(2,3)-a(2,1))**2 + (a(3,3)-a(3,1))**2
      beta = (cc*cc + aa*aa - d2)/(2.0d0*cc*aa)
      beta = dacos(beta)*180.0d0/pi

      d2 = (a(1,1)-a(1,2))**2 + (a(2,1)-a(2,2))**2 + (a(3,1)-a(3,2))**2
      gmma = (aa*aa + bb*bb - d2)/(2.0d0*aa*bb)
      gmma = dacos(gmma)*180.0d0/pi


      if (shift_cell) then
      if (taskid.eq.MASTER) then
 
        call util_date(dd)
        write(unit,1200) 
        write(unit,1210) dd(1:24)
        write(unit,1211) 

        write(unit,1220) aa * 0.529177d0
        write(unit,1221) bb * 0.529177d0
        write(unit,1222) cc * 0.529177d0
        write(unit,1223) alpha
        write(unit,1224) beta
        write(unit,1225) gmma

        write(unit,1230)

        write(unit,1240) 
        write(unit,1241) 
c        write(unit,1242)   ! causing problems with mercury??
        write(unit,1243) 
        write(unit,1244) 
        write(unit,1245) 

        do ii=1,ion_nion()
           frac(1) = b(1,1)*ion_rion(1,ii)
     >             + b(2,1)*ion_rion(2,ii)
     >             + b(3,1)*ion_rion(3,ii) + 0.5
           frac(2) = b(1,2)*ion_rion(1,ii)
     >             + b(2,2)*ion_rion(2,ii)
     >             + b(3,2)*ion_rion(3,ii) + 0.5
           frac(3) = b(1,3)*ion_rion(1,ii)
     >             + b(2,3)*ion_rion(2,ii)
     >             + b(3,3)*ion_rion(3,ii) + 0.5
c           WRITE(unit,1250) ion_aname_nocarat(ii),ii,frac
           WRITE(unit,1250) ion_aname_nocarat(ii),frac

        end do
      end if


      else


      if (taskid.eq.MASTER) then

        call util_date(dd)
        write(unit,1200)
        write(unit,1210) dd(1:24)
        write(unit,1211)

        write(unit,1220) aa * 0.529177d0
        write(unit,1221) bb * 0.529177d0
        write(unit,1222) cc * 0.529177d0
        write(unit,1223) alpha
        write(unit,1224) beta
        write(unit,1225) gmma

        write(unit,1230)

        write(unit,1240)
        write(unit,1241)
c        write(unit,1242)   ! causing problems with mercury??
        write(unit,1243)
        write(unit,1244)
        write(unit,1245)

        do ii=1,ion_nion()
           frac(1) = b(1,1)*ion_rion(1,ii)
     >             + b(2,1)*ion_rion(2,ii)
     >             + b(3,1)*ion_rion(3,ii) 
           frac(2) = b(1,2)*ion_rion(1,ii)
     >             + b(2,2)*ion_rion(2,ii)
     >             + b(3,2)*ion_rion(3,ii) 
           frac(3) = b(1,3)*ion_rion(1,ii)
     >             + b(2,3)*ion_rion(2,ii)
     >             + b(3,3)*ion_rion(3,ii) 
           !WRITE(unit,1250) ion_aname_nocarat(ii),ii,frac
           WRITE(unit,1250) ion_aname_nocarat(ii),frac

        end do
      end if

      end if

      return

 1200 FORMAT('data_nwchem_pspw')
 1210 FORMAT(/'_audit_creation_date   ',A)
 1211 FORMAT(
     > '_audit_creation_method    generated by PSPW module of NWChem')

 1220 FORMAT(//'_cell_length_a   ', F16.4)
 1221 FORMAT(  '_cell_length_b   ', F16.4)
 1222 FORMAT(  '_cell_length_c   ', F16.4)
 1223 FORMAT(  '_cell_angle_alpha', F16.4)
 1224 FORMAT(  '_cell_angle_beta ', F16.4)
 1225 FORMAT(  '_cell_angle_gamma', F16.4)

 1230 FORMAT(/'_symmetry_space_group_name_H-M     P1  ')

 1240 FORMAT(/'loop_')
 1241 FORMAT('_atom_site_type_symbol')
 1242 FORMAT('_atom_site_label')
 1243 FORMAT('_atom_site_fract_x')
 1244 FORMAT('_atom_site_fract_y')
 1245 FORMAT('_atom_site_fract_z')
 
c 1250 FORMAT(A2,6X,I4,3x,3F14.6)
 1250 FORMAT(A2,6X,3F14.6)
      end

*     ********************************
*     *                              *
*     *         ion_ecce             *
*     *                              *
*     ********************************

      subroutine ion_ecce()

#include "bafdecls.fh"
#include "beom.fh"
#include "ion.fh"

      integer geom1a
      common / ION2/ geom1a

      logical value
      integer rt(2),qt(2),tt(2),i,j
      double precision unita(3,3)

      double precision angstrom_to_au 
      data angstrom_to_au /1.88972598858d0/

      double precision lattice_unita
      external lattice_unita

      call nwpw_timing_start(50) 

      call dscal(nion*3, 1.0d0/angstrom_to_au, dbl_mb(r1(1)), 1)
      call ecce_print2('cartesian coordinates', mt_dbl,
     >     dbl_mb(r1(1)), 3, 3, nion)
      call dscal(nion*3, angstrom_to_au, dbl_mb(r1(1)), 1)
c

      value = BA_push_get(mt_dbl, (3*nion), 'rt',rt(2),rt(1))
      value = value.and.
     >        BA_push_get(mt_dbl, (nion),   'qt',qt(2),qt(1))
      value = value.and.
     >        BA_push_get(mt_byte,(16*nion),'tt',tt(2),tt(1))
      value = value.and.
     >        beom_cart_get(geom1a,nion,byte_mb(tt(1)),
     >                                dbl_mb(rt(1)),
     >                                dbl_mb(qt(1)))
      if (.not. value) call errquit('error ion_ecce',0,0)

*     *** print out lattice vectors ***
      do j=1,3
      do i=1,3
         unita(i,j) = lattice_unita(i,j)
      end do
      end do
      call ecce_print2('lattice vectors', mt_dbl, unita, 3, 3, 3)


      call ecce_print1('atomic charges', mt_dbl, dbl_mb(qt(1)), nion)
      call ion_ecce_tmp(nion, byte_mb(tt(1)))
      call ecce_print1_char('group name', "C1", 1)

      value = BA_pop_stack(tt(2))
      value = value.and.BA_pop_stack(qt(2))
      value = value.and.BA_pop_stack(rt(2))
      if (.not. value) call errquit('popping stack',0,0)

      call nwpw_timing_end(50) 
c
      return
      end


      subroutine ion_ecce_tmp(nion,tags)
      integer nion
      character*16 tags
      dimension tags(*)
      call ecce_print1_char('atomic tags', tags(1), nion)
      return
      end




      subroutine ion_Print_neighborlist()
      implicit none

      integer ii,jj
      real*8 dx,dy,dz,r,K,Rcut

      integer  ion_nion
      external ion_nion
      real*8   ion_rion
      external ion_rion
      character*4 ion_aname
      external    ion_aname

      K = 0.0d0
      Rcut = 2.5d0/0.529177d0

      write(*,*)
      write(*,*) "bond spring terms:"
      do ii=1,ion_nion()
         do jj=ii+1,ion_nion()
           dx = ion_rion(1,ii) - ion_rion(1,jj)
           dy = ion_rion(2,ii) - ion_rion(2,jj)
           dz = ion_rion(3,ii) - ion_rion(3,jj)
           call lattice_min_difference(dx,dy,dz)
           r = dsqrt(dx*dx + dy*dy + dz*dz)
           if (r.le.Rcut) then
            write(*,'(A,2I4,2F12.6,2A6)') "bond spring ",ii,jj,K,r,
     >                                    ion_aname(ii),ion_aname(jj)
           end if
         end do
      end do
      write(*,*)
      write(*,*) "angle spring terms:"
      write(*,*)

      return
      end


*     ***************************
*     *                         *
*     *      ion_nearest_index  *
*     *                         *
*     ***************************
*
*  returns the index of the ion nearest to rtest

      integer function ion_nearest_index(xx,yy,zz)
      implicit none
      real*8 xx,yy,zz

#include "bafdecls.fh"
#include "ion.fh"

*     **** local variables ****
      integer ii,iimin,i1,i2,i3
      real*8 d,dmin,x,y,z,dx,dy,dz
      real*8 a(3,3)

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita

      do i2=1,3
         do i1=1,3
            a(i1,i2) = lattice_unita(i1,i2)
         end do
      end do

      iimin = -1
      dmin = 9999999990099.0d0
      do ii=1,nion
         dx = dbl_mb(r1(1)+3*(ii-1))   - xx
         dy = dbl_mb(r1(1)+3*(ii-1)+1) - yy
         dz = dbl_mb(r1(1)+3*(ii-1)+2) - zz
         do i3=-1,1
            do i2=-1,1
               do i1=-1,1
                  x = dx + i1*a(1,1) + i2*a(1,2) + i3*a(1,3)
                  y = dy + i1*a(2,1) + i2*a(2,2) + i3*a(2,3)
                  z = dz + i1*a(3,1) + i2*a(3,2) + i3*a(3,3)
                  d = dsqrt(x*x + y*y + z*z)
                  if (d.lt.dmin) then
                     dmin = d
                     iimin = ii
                  end if
               end do
            end do
         end do
      end do

      ion_nearest_index = iimin
      return
      end

*     ***************************
*     *                         *
*     *      ion_rion_load      *
*     *                         *
*     ***************************
      subroutine ion_rion_load(rtdb,name1)
      implicit none
      integer rtdb
      character*(*) name1

#include "stdio.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint
      integer nion1,geom1,r11

*     **** external functions ****
      integer  ion_rion_ptr
      external ion_rion_ptr

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)
      r11 = ion_rion_ptr()

*     **** load name1 geometry ****
      if (.not.beom_create(geom1,name1))
     >   call errquit('ion_rion_load: beom_create?',1,GEOM_ERR)

      if (.not.beom_rtdb_load(rtdb,geom1,name1)) then
         if (oprint)
     >      write(luout,*) "ion_rion_load: Cannot load geometry:",
     >                 name1
         if (.not.beom_destroy(geom1))
     >      call errquit('ion_rion_load:beom_destroy?',2,GEOM_ERR)
         return
      end if

      if (.not.beom_ncent(geom1,nion1))
     >   call errquit('ion_rion_load: beom_ncent?',3,GEOM_ERR)

      call dcopy(3*nion1,0.0d0,0,dbl_mb(r11),1)

      if (.not. beom_cart_coords_get(geom1,
     >                          dbl_mb(r11)))
     >   call errquit('ion_rion_load: geom cart?',4,GEOM_ERR)

      if (.not.beom_destroy(geom1))
     >   call errquit('ion_rion_load: beom_destroy?',5,GEOM_ERR)

      return
      end 


*     ***************************
*     *                         *
*     *      ion_rion_save      *
*     *                         *
*     ***************************
      subroutine ion_rion_save(rtdb,name1)
      implicit none
      integer rtdb
      character*(*) name1

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "errquit.fh"

*     **** local variables ****
      integer nion1,geom1,r11,ii
      real*8 q
      character*16 t


*     **** external functions ****
      integer  ion_rion_ptr,ion_nion
      external ion_rion_ptr,ion_nion

      r11 = ion_rion_ptr()
      nion1 = ion_nion()

*     **** create name1 geometry ****
      if (.not.beom_create(geom1,name1))
     >   call errquit('ion_rion_save: beom_create?',1,GEOM_ERR)

      if (.not.beom_ncent_set(geom1,nion1))
     >   call errquit('ion_rion_save: beom_ncent_set?',2,GEOM_ERR)

      do ii=1,nion1
         call ion_tandq(ii,t,q)
         if (.not.beom_cent_set(geom1,ii,t,dbl_mb(r11+(ii-1)*3),q))
     >      call errquit('ion_rion_save: beom_cent_set?',3,GEOM_ERR)
      end do

      if (.not.beom_rtdb_delete(rtdb,name1))
     >   call errquit('ion_rion_save: beom_rtdb_delete?',4,GEOM_ERR)

      if (.not.beom_rtdb_store(rtdb,geom1,name1))
     >   call errquit('ion_rion_save: beom_rtdb_store?',5,GEOM_ERR)

      if (.not.beom_destroy(geom1))
     >   call errquit('ion_rion_save: beom_destroy?',6,GEOM_ERR)

      return
      end 



*     ***************************
*     *                         *
*     *   ion_load_extra_geom   *
*     *                         *
*     ***************************
      subroutine ion_load_extra_geom(first,rtdb,name1)
      implicit none
      logical first
      integer rtdb
      character*(*) name1

#include "stdio.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint
      integer nion1,geom

*     **** ion_extra common block ****
      integer nion_extra,rion_extra(2,10)
      common / ion_extra / rion_extra,nion_extra

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)
      if (first) nion_extra = 0

*     **** load name1 geometry ****
      if (.not.beom_create(geom,name1))
     >   call errquit('ion_load_extra_geom: beom_create?',1,GEOM_ERR)

      if (.not.beom_rtdb_load(rtdb,geom,name1)) then
         if (oprint)
     >      write(luout,*) "ion_load_extra_geom: Cannot load geometry:",
     >                 name1
         if (.not.beom_destroy(geom))
     >      call errquit('ion_load_extra_geom:beom_destroy?',1,GEOM_ERR)
         return
      end if

      if (.not.beom_ncent(geom,nion1))
     >   call errquit('ion_load_extra_geom: beom_ncent?',1,GEOM_ERR)

      if (.not.BA_alloc_get(mt_dbl,(3*nion1),'rion_extra',
     >                      rion_extra(2,nion_extra+1),
     >                      rion_extra(1,nion_extra+1)))
     >   call errquit('ion_load_extra_geom:out of memory',1,MA_ERR)

      call dcopy(3*nion1,0.0d0,0,dbl_mb(rion_extra(1,nion_extra+1)),1)

      if (.not. beom_cart_coords_get(geom,
     >                          dbl_mb(rion_extra(1,nion_extra+1))))
     >   call errquit('ion_load_extra_geom: geom cart?',1,GEOM_ERR)

      if (.not.beom_destroy(geom))
     >   call errquit('ion_load_extra_geom: beom_destroy?',1,GEOM_ERR)

      nion_extra = nion_extra + 1
      return
      end




*     ***************************
*     *                         *
*     *  ion_delete_extra_geom  *
*     *                         *
*     ***************************
      subroutine ion_delete_extra_geom()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ion_extra common block ****
      integer nion_extra,rion_extra(2,10)
      common / ion_extra / rion_extra,nion_extra

      logical value
      integer i

      value = .true.
      do i=1,nion_extra
         value = value.and.BA_free_heap(rion_extra(2,i))
      end do
      if (.not.value)
     >   call errquit('ion_delete_extra_geom: freeing heap',0,MA_ERR)
      return
      end

*     ********************************************
*     *                                          *
*     *             ion_rion_extra               *
*     *                                          *
*     ********************************************
      real*8 function ion_rion_extra(e,i,j)
      implicit none
      integer e,i,j

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ion_extra common block ****
      integer nion_extra,rion_extra(2,10)
      common / ion_extra / rion_extra,nion_extra

      ion_rion_extra = dbl_mb(rion_extra(1,e)+(i-1) + (j-1)*3)
      return
      end

*     ********************************************
*     *                                          *
*     *             ion_morph_extra              *
*     *                                          *
*     ********************************************
      subroutine ion_morph_extra(e1,e2)
      implicit none
      integer e1,e2

#include "bafdecls.fh"
#include "ion.fh"

*     **** ion_extra common block ****
      integer nion_extra,rion_extra(2,10)
      common / ion_extra / rion_extra,nion_extra

      integer i,j
      real*8 dr(3)

      do j=1,nion
         do i=1,3
            dr(i) = dbl_mb(rion_extra(1,e2)+(i-1)+(j-1)*3)
     >            - dbl_mb(rion_extra(1,e1)+(i-1)+(j-1)*3)
         end do
         call lattice_min_difference(dr(1),dr(2),dr(3))
         do i=1,3
            dbl_mb(r1(1)+(i-1)+(j-1)*3) 
     >        = dbl_mb(rion_extra(1,e1)+(i-1)+(j-1)*3) + 0.5d0*dr(i)
         end do
      end do
      return
      end

*     ********************************************
*     *                                          *
*     *             ion_t_morph_extra            *
*     *                                          *
*     ********************************************
      subroutine ion_t_morph_extra(t,e1,e2)
      implicit none
      real*8  t
      integer e1,e2

#include "bafdecls.fh"
#include "ion.fh"

*     **** ion_extra common block ****
      integer nion_extra,rion_extra(2,10)
      common / ion_extra / rion_extra,nion_extra

      integer i,j
      real*8 dr(3)

      do j=1,nion
         do i=1,3
            dr(i) = dbl_mb(rion_extra(1,e2)+(i-1)+(j-1)*3)
     >            - dbl_mb(rion_extra(1,e1)+(i-1)+(j-1)*3)
         end do
         call lattice_min_difference(dr(1),dr(2),dr(3))
         do i=1,3
            dbl_mb(r1(1)+(i-1)+(j-1)*3)
     >        = dbl_mb(rion_extra(1,e1)+(i-1)+(j-1)*3) + t*dr(i)
         end do
      end do
      return
      end




*     ***************************
*     *                         *
*     *      ion_disp_on        *
*     *                         * 
*     ***************************
      logical function ion_disp_on()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_disp_on = disp_on
      return
      end

*     ***************************
*     *                         *
*     *      ion_disp_energy    *
*     *                         * 
*     ***************************
      real*8 function ion_disp_energy()
      implicit none

***** ion common block ****
#include "bafdecls.fh"
#include "ion.fh"
#include "errquit.fh"

      integer i,j,ii,l,iz(2)
      real*8 edisp
      real*8 lat(3,3)
      real*8 g(3,1)
      real*8 g_lat(3,3)
      character*80 options,options1

      integer  ion_atn,ion_katm,control_version,inp_strlen
      external ion_atn,ion_katm,control_version,inp_strlen
      real*8   lattice_unita
      external lattice_unita
      character*80 control_options_disp
      external     control_options_disp

c      if (.not.BA_push_get(mt_int,nion,'iz',iz(2),iz(1)))
c     >   call errquit("ion_disp_energy:",0,MA_ERR)

      if (is_grimme2) then
!$OMP DO
         do i=1,nion_grimme2
            ii = int_mb(indx_grimme2(1)+i-1)
            dbl_mb(rion_grimme2(1)+3*(i-1))   = dbl_mb(r1(1)+3*(ii-1))
            dbl_mb(rion_grimme2(1)+3*(i-1)+1) = dbl_mb(r1(1)+3*(ii-1)+1)
            dbl_mb(rion_grimme2(1)+3*(i-1)+2) = dbl_mb(r1(1)+3*(ii-1)+2)
         end do
!$OMP END DO
      end if

      options1 = control_options_disp()
      l = inp_strlen(options1)

      if (control_version().eq.4) then
         options = options1
      else
         do j=1,3
         do i=1,3
            lat(i,j) = lattice_unita(i,j)
         end do
         end do

         options = options1(1:l)//' -pbc'
      end if

!$OMP MASTER
      call nwpwxc_vdw3_dftd3(options,nion_grimme2,int_mb(iz_grimme2(1)),
     >                     dbl_mb(rion_grimme2(1)),lat,edisp,g,g_lat)
!$OMP END MASTER
!$OMP BARRIER

c      if (.not.BA_pop_stack(iz(2)))
c     >   call errquit("ion_disp_energy:",1,MA_ERR)

      ion_disp_energy = edisp
      return
      end

*     ***************************
*     *                         *
*     *      ion_disp_force     *
*     *                         * 
*     ***************************
      subroutine ion_disp_force(fion)
      implicit none
      real*8 fion(3,*)


***** ion common block ****
#include "bafdecls.fh"
#include "ion.fh"
#include "errquit.fh"

      integer i,j,ii,l,gtmp(2),iz(2)
      real*8 edisp
      real*8 lat(3,3)
      real*8 g_lat(3,3)
      character*80 options,options1

      integer  ion_atn,ion_katm,control_version,inp_strlen
      external ion_atn,ion_katm,control_version,inp_strlen
      real*8   lattice_unita
      external lattice_unita
      character*80 control_options_disp
      external     control_options_disp

c      if (.not.BA_push_get(mt_int,nion,'iz',iz(2),iz(1)))
c     >   call errquit("ion_disp_energy: pushstack iz",0,MA_ERR)
c
c      if (.not.BA_push_get(mt_dbl,3*nion,'gtmp',gtmp(2),gtmp(1)))
c     >   call errquit("ion_disp_energy: pushstack gtmp",1,MA_ERR)

      if (is_grimme2) then
!$OMP DO
         do i=1,nion_grimme2
            ii = int_mb(indx_grimme2(1)+i-1)
            dbl_mb(rion_grimme2(1)+3*(i-1))   = dbl_mb(r1(1)+3*(ii-1))
            dbl_mb(rion_grimme2(1)+3*(i-1)+1) = dbl_mb(r1(1)+3*(ii-1)+1)
            dbl_mb(rion_grimme2(1)+3*(i-1)+2) = dbl_mb(r1(1)+3*(ii-1)+2)
         end do
!$OMP END DO
      end if
      
      options1= control_options_disp()
      l = inp_strlen(options1)

      if (control_version().eq.4) then
         options = options1(1:l)//' -grad'
      else
         do j=1,3
         do i=1,3
            lat(i,j) = lattice_unita(i,j)
         end do
         end do

         options = options1(1:l)//' -pbc -grad'
      end if

      !call dcopy(3*nion,0.0d0,0,dbl_mb(gtmp(1)),1)
      call Parallel_shared_vector_zero(.true.,3*nion_grimme2,
     >                                 dbl_mb(fion_grimme2(1)))

!$OMP MASTER
      call nwpwxc_vdw3_dftd3(options,nion_grimme2,int_mb(iz_grimme2(1)),
     >               dbl_mb(rion_grimme2(1)),lat,edisp,
     >               dbl_mb(fion_grimme2(1)),g_lat)
!$OMP END MASTER
!$OMP BARRIER

      if (is_grimme2) then
!$OMP DO
         do i=1,nion_grimme2
            ii = int_mb(indx_grimme2(1)+i-1)
            fion(1,ii) = fion(1,ii) - dbl_mb(fion_grimme2(1)+3*(i-1))
            fion(2,ii) = fion(2,ii) - dbl_mb(fion_grimme2(1)+3*(i-1)+1)
            fion(3,ii) = fion(3,ii) - dbl_mb(fion_grimme2(1)+3*(i-1)+2)
         end do
!$OMP END DO
      else
         call daxpy_omp(3*nion_grimme2,-1.0d0,
     >                  dbl_mb(fion_grimme2(1)),1,fion,1)
      end if

c      if (.not.BA_pop_stack(gtmp(2)))
c     >   call errquit("ion_disp_energy: popstak gtmp",1,MA_ERR)
c      if (.not.BA_pop_stack(iz(2)))
c     >   call errquit("ion_disp_energy: popstack iz",1,MA_ERR)


      return
      end


*     ***************************
*     *                         *
*     *      ion_disp_stress    *
*     *                         * 
*     ***************************
      subroutine ion_disp_stress(stress)
      implicit none
      real*8 stress(3,3)

***** ion common block ****
#include "bafdecls.fh"
#include "ion.fh"
#include "errquit.fh"

      integer i,j,ii,l,gtmp(2),iz(2)
      real*8 edisp
      real*8 lat(3,3)
      character*80 options,options1

      integer  ion_atn,ion_katm,control_version,inp_strlen
      external ion_atn,ion_katm,control_version,inp_strlen
      real*8   lattice_unita
      external lattice_unita
      character*80 control_options_disp
      external     control_options_disp

c      if (.not.BA_push_get(mt_int,nion,'iz',iz(2),iz(1)))
c     >   call errquit("ion_disp_energy: pushstack iz",0,MA_ERR)
c
c      if (.not.BA_push_get(mt_dbl,3*nion,'gtmp',gtmp(2),gtmp(1)))
c     >   call errquit("ion_disp_energy: pushstack gtmp",1,MA_ERR)

c      do ii=1,nion
c         int_mb(iz(1)+ii-1) = ion_atn(ion_katm(ii))
c      end do

      if (is_grimme2) then
!$OMP DO
         do i=1,nion_grimme2
            ii = int_mb(indx_grimme2(1)+i-1)
            dbl_mb(rion_grimme2(1)+3*(i-1))   = dbl_mb(r1(1)+3*(ii-1))
            dbl_mb(rion_grimme2(1)+3*(i-1)+1) = dbl_mb(r1(1)+3*(ii-1)+1)
            dbl_mb(rion_grimme2(1)+3*(i-1)+2) = dbl_mb(r1(1)+3*(ii-1)+2)
         end do
!$OMP END DO
      end if

      options1 = control_options_disp()
      l = inp_strlen(options1)
      
      if (control_version().eq.4) then
         options = options1(1:l)//' -grad'
      else
         do j=1,3
         do i=1,3
            lat(i,j) = lattice_unita(i,j)
         end do
         end do

         options = options1(1:l)//' -pbc -grad'
      end if

      !call dcopy(3*nion_grimme2,0.0d0,0,dbl_mb(fion_grimme2(1)),1)
      !call dcopy(9,0.0d0,0,stress,1)
      call Parallel_shared_vector_zero(.false.,3*nion_grimme2,
     >                                 dbl_mb(fion_grimme2(1)))
      call Parallel_shared_vector_zero(.true.,9,stress)

!$OMP MASTER
      call nwpwxc_vdw3_dftd3(options,nion_grimme2,int_mb(iz_grimme2(1)),
     >               dbl_mb(rion_grimme2(1)),lat,edisp,
     >               dbl_mb(fion_grimme2(1)),stress)
!$OMP END MASTER
!$OMP BARRIER

      !call dscal(9,-1.0d0,stress,1)


c      if (.not.BA_pop_stack(gtmp(2)))
c     >   call errquit("ion_disp_stress: popstak gtmp",1,MA_ERR)
c      if (.not.BA_pop_stack(iz(2)))
c     >   call errquit("ion_disp_stress: popstack iz",1,MA_ERR)

      return
      end




*     ***************************
*     *                         *
*     *   ion_check_distances   *
*     *                         * 
*     ***************************
      subroutine ion_check_distances()
      implicit none

***** ion common block ****
#include "bafdecls.fh"
#include "ion.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER, taskid
      parameter (MASTER=0)
      logical tooshort
      integer ii,jj,iiz,jjz,n1,n2,n3,iz(2)
      real*8  dx,dy,dz,dx0,dy0,dz0,r,thr,taux,tauy,tauz

c covalent radii (Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
c values for metals decreased by 10 %
c      data rcov/
c     .  0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67
c     ., 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54
c     ., 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09
c     ., 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39
c     ., 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26
c     ., 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57
c     ., 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53
c     ., 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32
c     ., 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58
c     ., 1.52, 1.53, 1.54, 1.55 /
c scaled by k2=4./3. and converted to atomic units
c autoang=0.52917726d0
      real*8 rcov(111)
      data rcov/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362, 4.00000000,
     . 4.00000000, 4.00000000, 4.00000000, 4.00000000, 4.00000000,
     . 4.00000000, 4.00000000, 4.00000000, 4.00000000, 4.00000000,
     . 4.00000000, 4.00000000, 4.00000000, 4.00000000, 4.00000000,
     . 4.00000000 /

*     ***** external functions ****
      integer  ion_atn,ion_katm,control_version
      external ion_atn,ion_katm,control_version
      real*8   lattice_unita,ion_rion
      external lattice_unita,ion_rion
      integer tid,Parallel_threadid
      external    Parallel_threadid

      call Parallel_taskid(taskid)
      tid = Parallel_threadid()

      if (.not.BA_push_get(mt_int,nion,'iz',iz(2),iz(1)))
     >   call errquit("ion_disp_energy: pushstack iz",0,MA_ERR)

      do ii=1,nion
         int_mb(iz(1)+ii-1) = ion_atn(ion_katm(ii))
      end do

      tooshort = .false.

!$OMP BARRIER
!$OMP SINGLE
      if (control_version().eq.4) then
         do ii=1,   nion-1
         do jj=ii+1,nion
            iiz = int_mb(iz(1)+ii-1)
            jjz = int_mb(iz(1)+jj-1)
            thr = 0.6d0*(rcov(iiz)+rcov(jjz))

            dx=ion_rion(1,ii)-ion_rion(1,jj)
            dy=ion_rion(2,ii)-ion_rion(2,jj)
            dz=ion_rion(3,ii)-ion_rion(3,jj)
            r=dsqrt(dx*dx+dy*dy+dz*dz)
            if (r.lt.thr) tooshort = .true.
           
         end do
         end do
      else
         do ii=1,   nion-1
         do jj=ii+1,nion
            iiz = int_mb(iz(1)+ii-1)
            jjz = int_mb(iz(1)+jj-1)
            thr = 0.6d0*(rcov(iiz)+rcov(jjz))
            dx0=ion_rion(1,ii)-ion_rion(1,jj)
            dy0=ion_rion(2,ii)-ion_rion(2,jj)
            dz0=ion_rion(3,ii)-ion_rion(3,jj)
            do n1=-1,1
            do n2=-1,1
            do n3=-1,1
               taux = n1*lattice_unita(1,1) 
     >              + n2*lattice_unita(1,2) 
     >              + n3*lattice_unita(1,3)
               tauy = n1*lattice_unita(2,1) 
     >              + n2*lattice_unita(2,2) 
     >              + n3*lattice_unita(2,3)
               tauz = n1*lattice_unita(3,1) 
     >              + n2*lattice_unita(3,2) 
     >              + n3*lattice_unita(3,3)
                dx = dx0 + taux
                dy = dy0 + tauy
                dz = dz0 + tauz

                r=dsqrt(dx*dx+dy*dy+dz*dz)
                thr=0.6d0*(rcov(iiz)+rcov(jjz))
                if (r.lt.thr) tooshort = .true.

            end do
            end do
            end do
         end do 
         end do 
      end if
!$OMP END SINGLE copyprivate(tooshort)

      if (.not.BA_pop_stack(iz(2)))
     >   call errquit("ion_disp_stress: popstack iz",1,MA_ERR)

      if (tooshort.and.(taskid.eq.MASTER).and.(tid.eq.MASTER)) then
          write(luout,*)
          write(luout,*)
     >    '======================================================'
          write(luout,*)
     >    '|| Some distances are very short. Check coordinates ||'
          write(luout,*)
     >    '======================================================'
          write(luout,*)
      end if

      !ion_check_distances = (.not.tooshort)
      return
      end


*     *****************************************
*     *                                       *
*     *         ion_hess_init                 *
*     *                                       *
*     *****************************************
*
*  The ion_hess routines are used to model springs.
*
      subroutine ion_hess_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "ion.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,np,MASTER
      parameter (MASTER=0)

      logical oprint,value,notperiodic
      integer iii,nat3,nhess,nhesst,i,k
      integer vcq(2),w1q(2),w2q(2)
      real*8  x
      character*(nw_max_path_len) filehess 

*     **** ion_hess_model common block ****
      logical found_hess
      common /ion_hess_model/ found_hess

*     **** external functions ****
      integer  control_version
      external control_version
      logical  control_hess_model
      external control_hess_model
      real*8   ion_amass
      external ion_amass
      character*4 ion_aname,ion_atom
      external    ion_aname,ion_atom
      character*9 ion_amm
      external    ion_amm



      call Parallel_taskid(taskid)
      oprint = taskid.eq.MASTER

      call control_hess_filename(filehess)
      if (taskid.eq.MASTER) then
        inquire(file=filehess,exist=found_hess)
        if (found_hess) then
           iii = 1
        else
           iii = 0
        end if
      end if
      call Parallel_Brdcst_ivalue(MASTER,iii)
      found_hess  = (iii.eq.1).and.control_hess_model()
      if (.not.found_hess) return

      notperiodic = (control_version().eq.4)

      if (oprint) then
         if (notperiodic) then
            write(luout,100) filehess(1:inp_strlen(filehess)),nion
         else
            write(luout,101) filehess(1:inp_strlen(filehess)),nion
         end if
      end if

      nat3   =  nion*3          ! 3-N (as in degrees of freedom)
      nhess  =  nat3*nat3       ! dimension of hessian
      nhesst =  nat3*(nat3+1)/2 ! dimension of lower triangular hessian

      !**** allocate ihesst, ihess, and ihessp ****
      value =  BA_alloc_get(mt_dbl,nhesst,'ihess',ihesst(2),ihesst(1))
      value = value.and.
     >     BA_alloc_get(mt_dbl,nhess,'ihess',ihess(2),ihess(1))
      value = value.and.
     >     BA_alloc_get(mt_dbl,nhess,'ihessp',ihessp(2),ihessp(1))
      value = value.and.
     >     BA_alloc_get(mt_dbl,nat3,'dq',dq(2),dq(1))
      value = value.and.
     >     BA_alloc_get(mt_dbl,nat3,'rrq0',rrq0(2),rrq0(1))
      if (.not.value)
     >   call errquit("ion_hess_energy_init: out of heap",1,MA_ERR)


      !**** read model from rtdb ****
      if (btdb_get(rtdb,"nwpw:ion_hess:ihess",mt_dbl,
     >                      nhess,dbl_mb(ihess(1))) .and.
     >    btdb_get(rtdb,"nwpw:ion_hess:rrq0",mt_dbl,
     >                       nat3,dbl_mb(rrq0(1)))) then
         if (oprint) write(luout,102)

      !**** read and set model from filehess and geom ****
      else
         !**** allocate tempory data for Eckart ****
         value = value.and.
     >        BA_push_get(mt_dbl,nat3*6,'vcq',vcq(2),vcq(1))
         value = value.and.
     >     BA_push_get(mt_dbl,nhess,'w1q',w1q(2),w1q(1))
         value = value.and.
     >     BA_push_get(mt_dbl,nhess,'w2q',w2q(2),w2q(1))
         if (.not.value)
     >      call errquit("ion_hess_init:out of stack",1,MA_ERR)

         if (taskid.eq.MASTER) then
            open(unit=69,file=filehess,form='formatted',status='old',
     >      err=99900,access='sequential')
            do iii = 0,(nhesst-1)
               read(69,*,err=99901,end=99902) x
               dbl_mb(ihesst(1)+iii) = x
            end do
            close(unit=69,status='keep')
         end if
         call Parallel_Brdcst_values(MASTER,nhesst,dbl_mb(ihesst(1)))
         call ion_dsquar(dbl_mb(ihesst(1)),dbl_mb(ihess(1)),nat3,nat3)
         call ion_eckart(oprint,nion,nat3,notperiodic,
     >                dbl_mb(ihess(1)),dbl_mb(ihessp(1)),
     >                dbl_mb(ihesst(1)),dbl_mb(r1(1)),
     >                dbl_mb(vcq(1)),dbl_mb(w1q(1)),dbl_mb(w2q(1)))

         call dcopy(nat3,dbl_mb(r1(1)),1,dbl_mb(rrq0(1)),1)

         !**** write hess model to rtdb ****
         if (.not.btdb_put(rtdb,"nwpw:ion_hess:ihess",mt_dbl,
     >                    nhess,dbl_mb(ihess(1))))
     >     call errquit("ion_hess_init:failed writing ihess",1,RTDB_ERR)

         if (.not.btdb_put(rtdb,"nwpw:ion_hess:rrq0",mt_dbl,
     >                     nat3,dbl_mb(rrq0(1))))
     >     call errquit("ion_hess_init:failed writing rrq0",1,RTDB_ERR)

         value =           MA_pop_stack(w2q(2))
         value = value.and.MA_pop_stack(w1q(2))
         value = value.and.MA_pop_stack(vcq(2))
         if (.not.value)
     >      call errquit("ion_hess_init:popping stack",1,MA_ERR)

      end if
      if (oprint) then
         write(luout,180)
         do i=1,nion
           write(luout,190) i,ion_aname(i),
     >            (dbl_mb(rrq0(1)+k-1+(i-1)*3),k=1,3),
     >                   ion_amass(i)/1822.89d0,ion_amm(i)
         end do
      end if


      return
  100 format(/" Hessian - Spring Algorithm",
     >    /"     - hessian filename = ",A,
     >    /"     - number of atoms = ",i8,
     >    /"     - aperiodic boundary conditions "/)
  101 format(/" Hessian - Spring Algorithm",
     >    /"     - hessian filename = ",A,
     >    /"     - number of atoms = ",i8,
     >    /"     - periodic boundary conditions "/)
  102 format("     - loading data from rtdb "/)
  180 FORMAT(/' zero position of ions (Bohrs):')
  190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)

99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('ion_hess_init: error opening file: "hess_file"',911,
     >       DISK_ERR)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('ion_hess_init: error reading file: "hess_file"',911,
     >       DISK_ERR)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit
     > ('ion_hess_init: unexpected EOF when reading file: "hess_file"',
     >  911,DISK_ERR)
      end

*     *****************************************
*     *                                       *
*     *         ion_hess_end                  *
*     *                                       *
*     *****************************************
      subroutine ion_hess_end()
      implicit none

#include "bafdecls.fh"
#include "ion.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

*     **** ion_hess_model common block ****
      logical found_hess
      common /ion_hess_model/ found_hess

      !**** deallocate ihesst, ihess, and ihessp ****
      if (found_hess) then
         value =           BA_free_heap(ihesst(2))
         value = value.and.BA_free_heap(ihess(2))
         value = value.and.BA_free_heap(ihessp(2))
         value = value.and.BA_free_heap(dq(2))
         value = value.and.BA_free_heap(rrq0(2))
         if (.not.value)
     >      call errquit("ion_hess_end:freeing heap",1,MA_ERR)

      end if

      return
      end

*     *****************************************
*     *                                       *
*     *             ion_hess_energy           *
*     *                                       *
*     *****************************************
      real*8 function ion_hess_energy()
      implicit none

***** ion common block ****
#include "bafdecls.fh"
#include "ion.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,np

*     **** ion_hess_model common block ****
      logical found_hess
      common /ion_hess_model/ found_hess

c     **** external functions ****
      real*8   ion_hess_energy_sub
      external ion_hess_energy_sub

      if (found_hess) then
         call Parallel_taskid(taskid)
         call Parallel_np(np)
         call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(dq(1)),1)
         call daxpy(3*nion,-1.0d0,dbl_mb(rrq0(1)),1,dbl_mb(dq(1)),1)

         ion_hess_energy = ion_hess_energy_sub(taskid,np,3*nion,
     >                                         dbl_mb(ihess(1)),
     >                                         dbl_mb(dq(1)))
      else
         ion_hess_energy = 0.0d0
      end if
      return
      end

      real*8 function ion_hess_energy_sub(taskid,np,N,A,dq)
      implicit none
      integer taskid,np,N
      real*8 A(N,N),dq(N)

*     **** local variables ****
      integer i,j,icount
      real*8 energy

      icount = 0
      energy = 0.0d0
      do j=1,N
         if (mod(icount,np).eq.taskid) then
            do i=1,N
               energy = energy + dq(i)*A(i,j)*dq(j)
            end do
         end if
         icount = icount + 1
      end do
      call Parallel_SumAll(energy)
      ion_hess_energy_sub = 0.5d0*energy
      return
      end

*     *****************************************
*     *                                       *
*     *             ion_hess_force            *
*     *                                       *
*     *****************************************
      subroutine ion_hess_force(G1)
      implicit none
      real*8 G1(*)

***** ion common block ****
#include "bafdecls.fh"
#include "ion.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,np

*     **** ion_hess_model common block ****
      logical found_hess
      common /ion_hess_model/ found_hess

      if (found_hess) then
         call Parallel_taskid(taskid)
         call Parallel_np(np)
         call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(dq(1)),1)
         call daxpy(3*nion,-1.0d0,dbl_mb(rrq0(1)),1,dbl_mb(dq(1)),1)

         call ion_hess_force_sub(taskid,np,3*nion,
     >                           dbl_mb(ihess(1)),
     >                           dbl_mb(dq(1)),G1)

      end if
      return
      end

      subroutine ion_hess_force_sub(taskid,np,N,A,dq,G1)
      implicit none
      integer taskid,np,N
      real*8 A(N,N),dq(N)
      real*8 G1(N)

*     **** local variables ****
      integer i,j,icount

      icount = 0
      call dcopy(N,0.0d0,0,G1,1)
      do j=1,N
         if (mod(icount,np).eq.taskid) then
            do i=1,N
               G1(i) = G1(i) - A(i,j)*dq(j)
               G1(j) = G1(j) - dq(i)*A(i,j)
            end do
         end if
         icount = icount + 1
      end do
      call Parallel_Vector_SumAll(N,G1)
      return
      end




*     *********************************************
*     *                                           *
*     *               ion_dsquar                  *
*     *                                           *
*     *********************************************
      subroutine ion_dsquar(A,S,N,MS)
C     
C EXPANSION OF LOWER TRIANGULAR PACKED ARRAY TO FULL SQUARE
C INPUT:
C     A   LOWER TRIANGULAR PACKED ARRAY
C     N   ROW ORDER OF A
C     MS  ROW DIMENSION OF S IN CALLING ROUTINE
C OUTPUT:
C     S  SQUARE PACKED ARRAY CORRESPONDING TO A
C
      implicit none
      integer n,ms
      double precision S(*),A(*),VAL
c       
      integer ij,ni,nj,ijs,jis,i,j
c     
      IJ=0
      NI=0
      DO 11 I=1,N
        NJ=0
        DO 10 J=1,I
          IJ=IJ+1
          VAL=A(IJ)
          IJS=I+NJ
          JIS=J+NI
          S(IJS)=VAL
          S(JIS)=VAL
          NJ=NJ+MS
   10   CONTINUE
        NI=NI+MS
   11 CONTINUE
      RETURN
      END

*     *********************************************
*     *                                           *
*     *               ion_eckart                  *
*     *                                           *
*     *********************************************
      subroutine ion_eckart(oprint,NAT,NAT3,removerotation,
     >                      HESS,HESSP,HESST,COORD,VC,
     >                      w1,w2)
      implicit none 
      logical          oprint
      integer          NAT,NAT3
      logical          removerotation
      double precision HESS(NAT3,NAT3)
      double precision HESSP(NAT3,NAT3)
      double precision HESST(NAT3*(NAT3+1)/2)
      double precision COORD(3,NAT)
      double precision VC(NAT3,6)
      double precision w1(nat3,*),w2(nat3,*)
c
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      double precision zero, one
      PARAMETER( ZERO=0.D0, ONE=1.D0 )

      DOUBLE PRECISION UNIVEC(3), TEST(6,6), VNORM, temp, dotval, rnorm
      double precision test_norm
      integer i, j, k, l, m, n, mu, nu, indx, iatom, iaxis, itemp
      integer nhess, nhesst,IMAX
      double precision ddot
      external ddot
C****
      nhess = nat3*nat3
      nhesst =  nat3*(nat3+1)/2 ! dimension of lower triangular hessian

C****
C**** construct translation unit vectors;  these are stored in the
C**** first three columns of array VC, the rotation vectors will
C**** be stored in the other 3 columns
C****
      IMAX=3
      if (removerotation) IMAX=6

      CALL Dfill (IMAX*NAT3,0.0d00,VC,1)
      VNORM = ONE/SQRT(dble(NAT))
      DO 10 I=1,3 ! unit translation vector for x-, y-, and z-directions
        DO 20 IATOM=1,NAT
           VC(3*(IATOM-1)+I,I) = VNORM
   20   CONTINUE
   10 CONTINUE
C****
C**** construct rotation unit vectors; store in last 3 columns of VC
C****
      if (removerotation) then
      DO 40 IAXIS=1,3         ! loop over rotation axes
        DO 50 J=1,3
          IF (J.EQ.IAXIS) THEN
            UNIVEC(J) = ONE
          ELSE
            UNIVEC(J) = ZERO
          END IF
   50   CONTINUE
        DO 60 IATOM=1,NAT       ! loop over atoms
          call vib_CROSS(COORD(1,IATOM),UNIVEC,
     &      VC(3*(IATOM-1)+1,IAXIS+3))
   60   CONTINUE
   40 CONTINUE
      end if
C****
C****   Schmidt orthogonalize the constraint vectors
C****
      DO 70 I=2,IMAX     ! orthogonalize vector I to each of vectors J
        DO 80 J=1,I-1
           dotval=ddot(nat3,vc(1,J),1,vc(1,i),1)
          DO 100 K=1,NAT3
             VC(K,I) = VC(K,I)-dotval*VC(K,J)
  100     CONTINUE
   80   CONTINUE
        RNORM = ddot(nat3,vc(1,I),1,vc(1,I),1)
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
          call dscal(nat3,rnorm,vc(1,I),1)
        ELSE ! dependent constraint vector (linear molecule; set to zero)
          DO 125 K=1,NAT3
            VC(K,I) = ZERO
  125     CONTINUE
          ITEMP = I-3
          if (oprint) then
          write(luout,*)'Dependent rotation vector no.',ITEMP
          write(luout,*)' found in ECKART; assuming linear geometry'
          end if
        END IF
   70 CONTINUE
C****
C**** remove rotations and translations from the Hessian matrix
C****
      call dcopy(nat3*nat3,hess,1,hessp,1)
      call dgemm('n','t',nat3,nat3,IMAX,1d0,vc,nat3,vc,nat3,
     .     0d0,w1,nat3)
      call dgemm('n','n',nat3,nat3,nat3,1d0,w1,nat3,hess,nat3,
     .     0d0,w2,nat3)
      call dgemm('n','n',nat3,nat3,nat3,1d0,w2,nat3,w1,nat3,
     .     1d0,hessp,nat3)
      call dgemm('n','n',nat3,nat3,nat3,-1d0,w1,nat3,hess,nat3,
     .     1d0,hessp,nat3)
      call dgemm('n','n',nat3,nat3,nat3,-1d0,hess,nat3,w1,nat3,
     .     1d0,hessp,nat3)

C****
C**** update triangular hessian to agree with HESS
C****
      DO 190 I=1,NAT3
        DO 200 J=1,I
          indx = i*(i-1)/2 + j
          HESST(indx) = HESSP(I,J)
  200   CONTINUE
  190 CONTINUE
C****
C**** place the projected Hessian in array HESS, and the difference
C**** between projected and unprojected in HESSP
C****
      DO 210 I=1,NAT3
        DO 220 J=1,I
          TEMP = HESS(I,J)
          HESS(I,J) = HESSP(I,J)
          HESSP(I,J) = HESS(I,J)-TEMP
          IF (DABS(HESSP(I,J)).LT.1.D-9) HESSP(I,J) = ZERO
          IF(J.NE.I) then
             HESSP(J,I)=HESSP(I,J)
             HESS(J,I)=HESS(I,J)
          endif
  220   CONTINUE
  210 CONTINUE
C****
C**** output projected Hessian
C****
      if (oprint) then
      if (util_print('eckart',print_debug)) then
        write(luout,*)
     >      'Hessian after projecting out external modes:'
        write(luout,*)HESS
        write(luout,*)
     >      'difference of Hessian with projected hessian:'
        write(luout,*)HESSP
      endif
      endif
C****
C**** construct Hessian in trans-rot subspace (should be zeroes)
C****
      DO 230 MU=1,IMAX
        DO 240 NU=1,IMAX
          TEST(MU,NU) = ZERO
          DO 250 K=1,NAT3
            DO 260 L=1,NAT3
              TEST(MU,NU) = TEST(MU,NU)+VC(K,MU)*HESS(K,L)*VC(L,NU)
  260       CONTINUE
  250     CONTINUE
  240   CONTINUE
  230 CONTINUE
      if (oprint) then
      if (util_print('eckart',print_medium)) then
        test_norm = ddot((IMAX*IMAX),test,1,test,1)
        write(luout,'(a,1pd10.4)')
     >      ' Projected Nuclear Hessian trans-rot subspace norm:',
     >      test_norm
        write(luout,*)
     >      '                        (should be close to zero!) '
      endif
      if (util_print('eckart',print_debug)) then
        write(luout,*)'Hessian projected into trans-rot subspace ',
     >      '(should be zeros):'
        write(luout,*)TEST
      endif
      end if

      return
      end




      





