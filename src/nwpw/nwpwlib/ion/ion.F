*
* $Id: ion.F,v 1.55 2009-03-26 00:18:59 bylaska Exp $
*

*     ***************************
*     *				*
*     *		ion_end 	*
*     *				*
*     ***************************
      subroutine ion_end()
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "ion.fh"

      logical value

      value = MA_free_heap(r2(2))
      value = value.and.MA_free_heap(r1(2))
      value = value.and.MA_free_heap(r0(2))
      value = value.and.MA_free_heap(katm(2))
      value = value.and.MA_free_heap(natm(2))
      value = value.and.MA_free_heap(atom(2))
      value = value.and.MA_free_heap(katm_qm(2))
      value = value.and.MA_free_heap(natm_qm(2))
      value = value.and.MA_free_heap(atom_qm(2))

      value = value.and.MA_free_heap(dti(2))
      value = value.and.MA_free_heap(amass(2))
      if (.not. value) call errquit('ion_end:free heap ',0, MA_ERR)

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

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "ion.fh"

      integer geom
      common / ION2/ geom

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

*     *** write out CIF file ***
      call CIF_write(rtdb)

      if (control_frac_coord()) call fcoord_to_frac(nion,dbl_mb(r1(1)))

      do i=1,nion
         value = geom_cent_get(geom,i,t,rxyz,q)
         value = geom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
      end do

*     **** if md code then write velocities ****
      if ((control_code().eq.2).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.14)) then
         value = geom_vel_set(geom,dbl_mb(r0(1)))
      end if

*     *****************************************************
*     **** update rtdb if md or steepest descent code  ****
*     *****************************************************
        value = geom_rtdb_delete(rtdb,'pspwgeometry')
      if ((control_code().eq.1).or.
     >    (control_code().eq.2).or.
     >    (control_code().eq.6).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.13).or.
     >    (control_code().eq.14)) then
        value = geom_rtdb_store(rtdb,geom,'pspwgeometry')
      end if
      value = geom_destroy(geom)

      if (control_frac_coord()) call fcoord_to_real(nion,dbl_mb(r1(1)))

*     *****************************************************
*     **** update rtdb if md or steepest descent code  ****
*     *****************************************************
      if ((control_code().eq.1).or.
     >    (control_code().eq.2).or.
     >    (control_code().eq.6).or.
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.13).or.
     >    (control_code().eq.14)) then


*       *****************************************
*       **** put together ions, and charges  ****
*       *****************************************
        call combine_pointcharge(rtdb)

      else
        value = geom_rtdb_delete(rtdb,'chargepspwgeometry')
        value = geom_rtdb_delete(rtdb,'pspwgeometry')
        value = geom_rtdb_delete(rtdb,'qmmmgeometry')
      end if
      return
      end

*     ***************************
*     *				*
*     *		ion_destroy 	*
*     *				*
*     ***************************
      subroutine ion_destroy()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer geom
      common / ION2/ geom
      logical value

      value = geom_destroy(geom)

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

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
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
      integer     i
      double precision q,dt  !,vgx,vgy,vgz
      character*16     t
      logical value

      integer MASTER,taskid
      parameter (MASTER=0)

*     **** external functions ****
      character*4 ion_aname
      integer     control_code
      external    ion_aname
      external    control_code
      
      logical          control_frac_coord,control_init_velocities
      logical          parseqmmm
      integer          control_version
      double precision control_ion_time_step,ion_amass_geom
      external         control_frac_coord,control_init_velocities
      external         parseqmmm
      external         control_version
      external         control_ion_time_step,ion_amass_geom

      call Parallel_taskid(taskid)

*     ***** set periodic ****
      periodic = .true.
      if (control_version().eq.4) periodic = .false.

*     *****************************************
*     **** seperate ions and charges       ****
*     *****************************************
      call seperate_pointcharge(rtdb)


*     **************************
*     **** read in ion data ****
*     **************************
      value = geom_create(geom,'pspwgeometry')
      value = geom_rtdb_load(rtdb,geom,'pspwgeometry')

*     **** get the number ions ****
      value = geom_ncent(geom,nion)

*     *** set ke_count, ke_total, and kg_total to zero ***
      ke_count = 0
      ke_total = 0.0d0
      kg_total = 0.0d0

*     ***** allocate ion data structure *****
      value = MA_alloc_get(mt_dbl,(3*nion),'r2',r2(2),r2(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nion),'r0',r0(2),r0(1))
      value = value.and.
     >       MA_alloc_get(mt_dbl,(nion),'dti',dti(2),dti(1))
      value = value.and.
     >       MA_alloc_get(mt_dbl,(nion),'amass',amass(2),amass(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      call dcopy(3*nion,0.0d0,0,dbl_mb(r2(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(r1(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(r0(1)),1)
      nion_qm = 0
      nion_mm = 0
      do i=1,nion
         value = geom_cent_get(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
         if (.not. value) call errquit('error reading ions',0, GEOM_ERR)
         if (parseqmmm(t)) then
            nion_mm = nion_mm + 1
         else
            nion_qm = nion_qm + 1
         end if
      end do

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
      if (periodic) call incell1(nion_qm,dbl_mb(r1(1)))

      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)


*     **** define dti ****
      dt = control_ion_time_step()
      do i=1,nion
         dbl_mb(dti(1)  +i-1) = dt*dt/(ion_amass_geom(i))
         dbl_mb(amass(1)+i-1) = ion_amass_geom(i)
      end do

c     value = geom_destroy(geom)
c     call ion_shift()

*     **** get velocities if md code ****
      if ((control_code().eq.2).or. 
     >    (control_code().eq.7).or.
     >    (control_code().eq.11).or.
     >    (control_code().eq.14)) then
*         **** get velocities ****
         value = geom_vel_get(geom,dbl_mb(r0(1)))

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


      ion_init = value
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
      character*9 string
      logical  ion_aismm
      external ion_aismm

      if (ion_aismm(i)) then
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
#include "geom.fh"
#include "errquit.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      character*16 t

*     **** external functions ****
      logical  parseqmmm
      external parseqmmm

      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' ion_aismm  failed ',i, GEOM_ERR)

      ion_aismm = parseqmmm(t)
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
#include "geom.fh"
#include "errquit.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      integer k,kp,l
      character*2  symbol
      character*4  symbol_out
      character*16 t,name
      real*8 q

*     **** external functions ****
      logical  parseqmmm
      external parseqmmm

      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' ion_aname  failed ',i, GEOM_ERR)

      if (.not. geom_tag_to_element(t,symbol,name,q))
     >        call errquit(' ion_aname  failed ',i, GEOM_ERR)

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
#include "geom.fh"
#include "errquit.fh"


*     **** local variables ****
      integer k,kp,l
      character*2  symbol
      character*4  symbol_out
      character*16 t,name
      real*8 q

*     **** external functions ****
      logical  parseqmmm
      external parseqmmm

      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' ion_aname_geom  failed ',i, GEOM_ERR)

      if (.not. geom_tag_to_element(t,symbol,name,q))
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
#include "geom.fh"
#include "errquit.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      character*2  symbol
      character*16 t,name
      real*8 q


      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit('ion_aname_nocarat failed ',i,GEOM_ERR)

      if (.not. geom_tag_to_element(t,symbol,name,q))
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
#include "geom.fh"
#include "errquit.fh"


      integer geom
      common / ION2/ geom

*     **** local variables ****
      real*8 mass

      if(.not.geom_mass_get(geom,i,mass))
     >  call errquit(' geom_mass_get  failed ',i, GEOM_ERR)

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

#include "mafdecls.fh"
#include "errquit.fh"
#include "ion.fh"

      ion_amass = dbl_mb(amass(1)+i-1)
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
#include "mafdecls.fh"
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
*     *						*
*     *	     nwpw_bqbq		*
*     *						*
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
#include "geom.fh"
#include "errquit.fh"


      integer geom
      common / ION2/ geom

*     **** local variables ****
      real*8       rxyz(3), q
      character*16 t

      if (.not.geom_cent_get(geom,i,t,rxyz,q))
     >  call errquit('ion_q:geom_cent_get failed ',i, GEOM_ERR)

      ion_q = q
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
#include "mafdecls.fh"
#include "rtdb.fh"


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

      if (.not.rtdb_get(rtdb,name(1:l),mt_dbl,1,q))
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

      double precision ave,T

      T = 0.0d0
      ave = (ke_total/dble(ke_count))
      if (nion.gt.2) then
         T = 2.0d0*ave/(3.0d0*nion-6.0d0)/kb
      else
         T = 2.0d0*ave/kb
      end if

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


#include "mafdecls.fh"

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

      integer i
      real*8  vgx,vgy,vgz,T,Tf

*     **** external functions ****
      logical  control_COM_shift,control_init_velocities
      real*8   control_rti,ion_amass,util_random
      real*8   control_init_velocities_temperature
      external control_COM_shift,control_init_velocities
      external control_rti,ion_amass,util_random
      external control_init_velocities_temperature

      !**** generate random initial velocities *****
      if (control_init_velocities()) then
         call Parallel_taskid(taskid)
         call dcopy(3*nion,0.0d0,0,dbl_mb(r0(1)),1)
         Tf = control_init_velocities_temperature()

         if (taskid.eq.MASTER) then
            vgx = util_random(394) !*** initialize random sequence ****
            do i=1,3*nion
               dbl_mb(r0(1)+i-1) = 0.1d0*(0.5d0-util_random(0))
            end do
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
      call dscal((3*nion),control_rti(),dbl_mb(r0(1)),1)
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
*     *				*
*     *		ion_katm	*
*     *				*
*     ***************************

      integer function ion_katm(i)
      implicit none
      integer i

#include "mafdecls.fh"

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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"

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

#include "mafdecls.fh"
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

#include "mafdecls.fh"

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

#include "mafdecls.fh"
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

#include "mafdecls.fh"

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

#include "mafdecls.fh"

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

#include "mafdecls.fh"
#include "ion.fh"

      ion_rion_ptr = r1(1)
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

#include "mafdecls.fh"

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

#include "mafdecls.fh"

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

#include "mafdecls.fh"
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
      if (periodic) then
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



*     ***************************
*     *				*
*     *	   ion_optimize_step   	*
*     *				*
*     ***************************

      subroutine ion_optimize_step(fion)
      implicit none
      double precision fion(3,*)

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

*     **** local variables ****
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
      if (periodic) then
         call incell1(nion_qm,dbl_mb(r2(1)))
         call pspw_qmmm_incell1(dbl_mb(r2(1)))
      end if

      return
      end
      subroutine ion_opt_sub(nion,dti,fion,r1,r2)
      implicit none
      integer nion
      real*8 dti(nion)
      real*8 fion(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii

      do ii=1,nion
      do i=1,3
         r2(i,ii) = r1(i,ii) + dsqrt(dti(ii))*fion(i,ii)
      end do
      end do
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

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
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
      if (periodic) then
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


*     **** determine current kinetic energy and add to running average ****
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

      return
      end

      subroutine ion_vert_sub1(nion,sa1,sa2,dti,fion,r0,r1,r2)
      implicit none
      integer nion
      real*8 sa1,sa2,dti(nion)
      real*8 fion(3,nion),r0(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii
      do ii=1,nion
      do i=1,3
        r2(i,ii) = 2.0d0*sa1*r1(i,ii) - sa2*r0(i,ii) 
     >           + sa1*dti(ii)*fion(i,ii)
      end do
      end do
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
*      r0 = old r1
*      r1 = MC step (r1)
*
      subroutine ion_MC_step(alpha)
      implicit none
      double precision alpha

#include "util.fh"
#include "mafdecls.fh"
#include "ion.fh"

*     **** local variables ****
      integer i
      double precision dt,h

*     *** external functions ****
      logical          control_rotation,control_translation
      double precision control_ion_time_step,ion_amass
      external         control_rotation,control_translation
      external         control_ion_time_step,ion_amass


      dt = control_ion_time_step()   !** not really needed, but used to make r2 = r1 + v*dt 
      h = 1.0d0/(dt)                 !** for rotation removal.

      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r0(1)),1)

*     **** do a MC step put MC step into r2 temporarilly ***
      i = int(util_random(0)*nion) + 1

c only update z direction for dimer test
c      dbl_mb(r1(1)+(i-1)*3)   = dbl_mb(r0(1)+(i-1)*3) 
c     >                        + (2.0d0*util_random(0)-1.0d0)*alpha
c      dbl_mb(r1(1)+(i-1)*3+1) = dbl_mb(r0(1)+(i-1)*3+1) 
c     >                        + (2.0d0*util_random(0)-1.0d0)*alpha
      dbl_mb(r1(1)+(i-1)*3+2) = dbl_mb(r0(1)+(i-1)*3+2) 
     >                        + (2.0d0*util_random(0)-1.0d0)*alpha

      

      call ion_same_FixIon(dbl_mb(r0(1)),dbl_mb(r1(1))) !** Stupid IBM **


*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r1(1)), dbl_mb(r0(1)))



*     **** make sure ions are in cell ***
      if (periodic) then
         call incell1(nion_qm,dbl_mb(r1(1)))
         call pspw_qmmm_incell1(dbl_mb(r1(1)))
      end if


*     **** remove translation? ***
      if (.not.control_translation())
     > call remove_center_mass(dbl_mb(r1(1)),dbl_mb(r0(1)))

*     **** remove rotation? ***
         !r1 = r0 + v*dt = r0+rand_alpha
         !v = rand_alpha/dt = (r1-r0)/dt = h*(r1-r0)
         !v = h*(r1-r0)
         !r0 = r1 - (v/h) = r1 - (h*(r1-r0)/h) = r1 - r1 + r0 = r0
      if (.not.control_rotation()) 
     > call ion_remove_rotation(h,nion,
     >                          dbl_mb(r0(1)),
     >                          dbl_mb(r0(1)),
     >                          dbl_mb(r1(1)))

      return
      end

*     ***************************
*     *                         *
*     *    ion_MC_reject_step   *
*     *                         *
*     ***************************

*  If the step is rejected than copy r0 back to r1.
*
      subroutine ion_MC_reject_step()
      implicit none

#include "mafdecls.fh"
#include "ion.fh"

      call dcopy((3*nion),dbl_mb(r0(1)),1,dbl_mb(r1(1)),1)
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

*     *** remove rotation ***
      hinv = 1.0d0/h
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


*     ***************************
*     *				*
*     *	   ion_newton_step   	*
*     *				*
*     ***************************

      subroutine ion_newton_step(fion,alpha)
      implicit none
      double precision fion(3,*)
      double precision alpha

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count
      

*     **** local variables ****
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
      if (periodic) then
         call incell2(nion_qm,dbl_mb(r2(1)),
     >                        dbl_mb(r1(1)))
         call pspw_qmmm_incell2(dbl_mb(r2(1)),
     >                          dbl_mb(r1(1)))
      end if

*     **** remove translation? ***
      if (.not.control_translation()) 
     > call remove_center_mass(dbl_mb(r2(1)),dbl_mb(r1(1)))

*     **** determine current kinetic energy and add to running average ***
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

      return
      end
      subroutine ion_newt_sub(nion,dta,dti,fion,r0,r1,r2)
      implicit none
      integer nion
      real*8 dta,dti(nion)
      real*8 fion(3,nion),r0(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii
      do ii=1,nion
      do i=1,3
        r2(i,ii)=r1(i,ii)+dta*r0(i,ii)+0.5d0*dti(ii)*fion(i,ii)
      end do
      end do
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

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
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
      if (periodic) then
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


*     **** determine current kinetic energy and add to running average ****
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

      return
      end

      subroutine ion_nose_sub(nion,ssr,smr,dti,fion,r0,r1,r2)
      implicit none
      integer nion
      real*8 ssr,smr,dti(nion)
      real*8 fion(3,nion),r0(3,nion),r1(3,nion),r2(3,nion)
      integer i,ii
      do ii=1,nion
      do i=1,3
        r2(i,ii) = 2.0d0*ssr*r1(i,ii) - smr*r0(i,ii) 
     >           + ssr*dti(ii)*fion(i,ii)
      end do
      end do
      return
      end

      subroutine ion_vel_sub(nion,h,r0,r2)
      implicit none
      integer nion
      real*8 h,r0(3,nion),r2(3,nion)
      integer i,ii
      do ii=1,nion
      do i=1,3
        r0(i,ii) = h*(r2(i,ii)-r0(i,ii))
      end do
      end do
      return
      end




*     ***************************
*     *				*
*     *		ion_shift	*
*     *				*
*     ***************************

      subroutine ion_shift()
      implicit none

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"


      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r0(1)),1)
      call dcopy((3*nion),dbl_mb(r2(1)),1,dbl_mb(r1(1)),1)

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

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      call dcopy((3*nion),rion,1,dbl_mb(r1(1)),1)

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
        write(unit,1242)   ! causing problems with mercury??
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
           WRITE(unit,1250) ion_aname_nocarat(ii),ii,frac

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
        write(unit,1242)   ! causing problems with mercury??
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
           WRITE(unit,1250) ion_aname_nocarat(ii),ii,frac

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
 
 1250 FORMAT(A2,6X,I4,3x,3F14.6)
      end

*     ********************************
*     *                              *
*     *         ion_ecce             *
*     *                              *
*     ********************************

      subroutine ion_ecce()

#include "mafdecls.fh"
#include "geom.fh"
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
 

      call dscal(nion*3, 1.0d0/angstrom_to_au, dbl_mb(r1(1)), 1)
      call ecce_print2('cartesian coordinates', mt_dbl,
     >     dbl_mb(r1(1)), 3, 3, nion)
      call dscal(nion*3, angstrom_to_au, dbl_mb(r1(1)), 1)
c

      value = MA_push_get(mt_dbl, (3*nion), 'rt',rt(2),rt(1))
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),   'qt',qt(2),qt(1))
      value = value.and.
     >        MA_push_get(mt_byte,(16*nion),'tt',tt(2),tt(1))
      value = value.and.
     >        geom_cart_get(geom1a,nion,byte_mb(tt(1)),
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

      value = MA_pop_stack(tt(2))
      value = value.and.MA_pop_stack(qt(2))
      value = value.and.MA_pop_stack(rt(2))
      if (.not. value) call errquit('popping stack',0,0)

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

