c $Id: md_xs.F 25353 2014-03-17 14:12:25Z jhammond $

c     **********************************************
c     *                                            *
c     *                task_md_xs                  *
c     *                                            *
c     **********************************************
      logical function task_md_xs(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"

*     **** local variables ****
      double precision kb
      parameter (kb=3.16679d-6)
      real*8 autoatm
      parameter (autoatm =290.360032539d6)

      integer taskid_i,taskid_j,npi,npj
      integer taskid,np,MASTER
      parameter (MASTER=0)

      character*(nw_max_path_len) filehess         ! name of file storing lower triangular packed hessian
      character*(nw_max_path_len) motion_filename  ! name of motion filename
      logical oprint,value,DW_modeling,ignore,found_hess
      integer nframes,natoms,stride,natoms1,center,nkf,nr,nat3,nhess
      integer iii,i,k,icount,icpu,pto,pfrom,iw,steps(2),nhesst
      real*8  rmax,rcut,x,y,z,vx,vy,vz,time1,volume,unita(3,3),ft(3)
      real*8  xmin,xmax,ymin,ymax,temperature,cx,cy,cz,gx,gy,gz
      real*8  v0,v1,v2,mc_dr,mc_dV,mc_temperature,mc_beta,mc_pressure
      real*8  c0,c1,c2,dA,ddx,ddv,mc_atom_direction(3),aratio

      integer seed,mc_napply,mc_ngroups,mc_group_size,mc_algorithm
      integer ihesst(2),ihess(2),ihessp(2),dq0(2),dq1(2)
      integer mc_group_start(2),mc_group_end(2),mc_group(2)


*     **** external functions ****
      logical  ion_init,ion_q_FixIon,ion_q_xyzFixIon
      external ion_init,ion_q_FixIon,ion_q_xyzFixIon
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm
      external ion_nion,ion_natm,ion_katm,ion_nkatm
      integer  control_mc_seed,control_mc_algorithm
      external control_mc_seed,control_mc_algorithm
      real*8   control_mc_step_size,control_mc_Temperature
      external control_mc_step_size,control_mc_Temperature
      real*8   control_mc_ddx,control_mc_aratio,control_mc_ddv
      external control_mc_ddx,control_mc_aratio,control_mc_ddv
      real*8   control_mc_pressure,control_mc_volume_step
      external control_mc_pressure,control_mc_volume_step
      logical  control_rotation,control_translation
      external control_rotation,control_translation
      integer  control_bo_steps_in,control_bo_steps_out
      external control_bo_steps_in,control_bo_steps_out
      character*14 ion_q_xyzFixIon_label
      external     ion_q_xyzFixIon_label

      real*8   ion_amass,ion_rion
      external ion_amass,ion_rion
      character*4 ion_aname,ion_atom
      external    ion_aname,ion_atom
      character*9 ion_amm
      external    ion_amm


      task_md_xs = .false.
      value = .true.

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call Parallel2d_np_i(npi)
      call Parallel2d_np_j(npj)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

      DW_modeling = .false.

*     **** initialize Ion ****
      if (.not.ion_init(rtdb)) then
         write(luout,*) 'task_md_xs: error reading geometry'
         task_md_xs = .false.
         return
      end if
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)


      natoms = ion_nion()
      nat3   =  natoms*3        ! 3-N (as in degrees of freedom)
      nhess  =  nat3*nat3       ! dimension of hessian
      nhesst =  nat3*(nat3+1)/2 ! dimension of lower triangular hessian

      call util_file_name('DW_motion',.false.,.false.,motion_filename)
      call util_file_name('hess',  .false., .false.,filehess)
      if (taskid.eq.MASTER) then
        inquire(file=filehess,exist=found_hess)
        if (found_hess) then
           iii = 1
        else
           iii = 0
        end if
      end if
      call Parallel_Brdcst_ivalue(MASTER,iii)
      found_hess  = (iii.eq.1)
      DW_modeling = found_hess

      if (DW_modeling) then

         if (.not.rtdb_get(rtdb,'md_xs:DW_Harmonic_Temperature',
     >                     mt_dbl,1,temperature)) then
            temperature = 298.15d0
         end if
         if (.not.rtdb_get(rtdb,'md_xs:DW_Harmonic_steps',
     >                     mt_int,2,steps)) then
            steps(1) = 500
            steps(2) = 500
         end if


         !**** allocate ihesst, ihess, and ihessp ****
         value =  BA_push_get(mt_dbl,nhesst,'ihess',ihesst(2),ihesst(1))
         value = value.and.
     >        BA_push_get(mt_dbl,nhess,'ihess',ihess(2),ihess(1))
         value = value.and.
     >        BA_push_get(mt_dbl,nhess,'ihessp',ihessp(2),ihessp(1))
         value = value.and.
     >        BA_push_get(mt_dbl,nhess,'nat3',dq0(2),dq0(1))
         value = value.and.
     >        BA_push_get(mt_dbl,nhess,'nat3',dq1(2),dq1(1))

         if (.not.value)
     >      call errquit("task_md_xs: out of stack",1,MA_ERR)

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
         call md_xs_dsquar(dbl_mb(ihesst(1)),dbl_mb(ihess(1)),nat3,nat3)


*        **** Initialize Monte-Carlo Parameters ****
         mc_dr  = control_mc_step_size()
         mc_dV  = control_mc_volume_step()
         mc_algorithm  = control_mc_algorithm()
         call control_mc_atom_direction(mc_atom_direction)
         call control_mc_ngroups(mc_napply,mc_ngroups,mc_group_size)
         if (mc_ngroups.lt.1) then

            value = BA_alloc_get(mt_int,1,
     >                       'mc_group_start',
     >                        mc_group_start(2),
     >                        mc_group_start(1))
            value = value.and.BA_alloc_get(mt_int,1,
     >                       'mc_group_end',
     >                        mc_group_end(2),
     >                        mc_group_end(1))
            value = value.and.BA_alloc_get(mt_int,1,
     >                       'mc_group',
     >                        mc_group(2),
     >                        mc_group(1))
            if (.not.value)
     >        call errquit('cgmontecarlov1:out of heap memory',0,MA_ERR)
         else
            value = BA_alloc_get(mt_int,mc_ngroups,
     >                       'mc_group_start',
     >                        mc_group_start(2),
     >                        mc_group_start(1))
            value = value.and.BA_alloc_get(mt_int,mc_ngroups,
     >                       'mc_group_end',
     >                        mc_group_end(2),
     >                        mc_group_end(1))
            value = value.and.BA_alloc_get(mt_int,mc_group_size,
     >                       'mc_group',
     >                        mc_group(2),
     >                        mc_group(1))
            if (.not.value)
     >        call errquit('cgmontecarlov1:out of heap memory',1,MA_ERR)
            call control_mc_groups(int_mb(mc_group_start(1)),
     >                             int_mb(mc_group_end(1)),
     >                             int_mb(mc_group(1)))
         end if

         seed           = control_mc_seed()
         mc_temperature = control_mc_Temperature()
         mc_pressure    = control_mc_pressure()
         mc_beta = util_random(control_mc_seed()) !**seed the random number generator
         mc_beta = 1.0d0/(kb*mc_temperature)

         aratio      = control_mc_aratio()

         ddx = control_mc_ddx()
         dA  = 1.0d0/(aratio*aratio-aratio)
         c0  = 1.0d0-ddx
         c1  = (aratio*aratio*(2*ddx) - (ddx))*dA
         c2  = (-aratio      *(2*ddx) + (ddx))*dA

         if (mc_algorithm.eq.2) then
            ddv = control_mc_ddv()
            v0  = 1.0d0-ddv
            v1  = (aratio*aratio*(2*ddv) - (ddv))*dA
            v2  = (-aratio      *(2*ddv) + (ddv))*dA
         end if

      end if

      if (oprint) then
         write(luout,*)
         write(luout,100)
         write(luout,110)
         write(luout,120)
         write(luout,110)
         write(luout,130)
         write(luout,140)
         write(luout,110)
         write(luout,100)
         write(luout,*)
         if (DW_modeling) then
            write(luout,*) 'Debye-Waller Harmonic Modeling'
            write(luout,*) '   - hessian filename:    ',
     >                     filehess(1:inp_strlen(filehess))
            write(luout,*) '   - DW motion filename:  ',
     >                  motion_filename(1:inp_strlen(motion_filename))
            write(luout,'(A,F8.3)') '    - temperature:         ',
     >                            temperature
            write(luout,'(A,I8)') '    - equilibration steps: ',steps(1)
            write(luout,'(A,I8)') '    - averaging steps:     ',steps(2)
         else
            write(luout,*) 'No Debye-Waller Factors'
            write(luout,*) '   - need to calculate frequencies'
         end if
         write(luout,160)
         write(luout,170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_xyzFixIon(I)) then
           write(luout,194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           end if
         end do
         write(luout,200) cx,cy,cz
         write(luout,210) gx,gy,gz
         write(luout,1310)

         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         if (mc_algorithm.eq.0) then
            write(luout,1311) 'NVE Microcanonical Ensemble - ',
     >                        'HyperVolume Monte Carlo Algorithm'
         else if (mc_algorithm.eq.1) then
            write(luout,1311) 'NVT Canonical Ensemble - ',
     >                        'Metropolis Algorithm'
         else if (mc_algorithm.eq.2) then
            write(luout,1311) 'NPT Canonical Ensemble - ',
     >                        'Metropolis Algorithm'
         end if

         write(luout,1320) mc_dr,
     >    control_bo_steps_in()*control_bo_steps_out(),
     >    control_bo_steps_in(),control_bo_steps_out()

         write(luout,1322) mc_temperature,mc_beta,seed
         if (mc_algorithm.eq.2)
     >      write(luout,1326) mc_pressure,mc_pressure*autoatm,
     >                        ddv,v0,v1,v2
         write(luout,1323) aratio,ddx,c0,c1,c2
         write(luout,1324) mc_napply,mc_atom_direction
         write(luout,1325) mc_ngroups

      end if

*                |***************************|
******************  simple Monte-Carlo loop  **********************
*                |***************************|
      call dcopy(nat3,0.0d0,0,dbl_mb(dq0(1)),1)
      call dcopy(nat3,0.0d0,0,dbl_mb(dq1(1)),1)




      if (DW_modeling) then
c        **** deallocate stack ****
         value =           BA_pop_stack(dq1(2))
         value = value.and.BA_pop_stack(dq0(2))
         value = value.and.BA_pop_stack(ihessp(2))
         value = value.and.BA_pop_stack(ihess(2))
         value = value.and.BA_pop_stack(ihesst(2))
         if (.not.value)
     >      call errquit("task_md_xs: popping stack",2,MA_ERR)

      end if
      call ion_end()
      call ion_end_FixIon()

      !call md_xs_analysis(rtdb,motion_filename)


      task_md_xs = value
      return

  100 FORMAT(10x,' ************************************************')
  110 FORMAT(10x,' *                                              *')
  120 FORMAT(10x,' *           FEFF6L EXAFS Calculation           *')
  130 FORMAT(10x,' *  This interface to FEFF6L was developed by   *')
  140 FORMAT(10x,' *  by Eric J. Bylaska.                         *')

  160 FORMAT(/' atomic composition:')
  170 FORMAT(7(5X,A2,':',I5))
  180 FORMAT(/' initial position of ions (Bohrs):')
  190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)
  191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed ',A)
  193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
  194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)
  200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
  210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1130 FORMAT(/' supercell:')
 1131 FORMAT(5x,' volume : ',F12.1)
 1141 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1142 FORMAT(5x,'             a2=<',3f8.3,' >')
 1143 FORMAT(5x,'             a3=<',3f8.3,' >')
 1270 FORMAT(/' technical parameters for minimizer:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1310 FORMAT(/' Monte-Carlo parameters:')
 1311 FORMAT(5X, ' MC Algorithm= ',A,A)
 1320 FORMAT(5X, ' MC displacement=',F10.2,5X,' iterations=',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1321 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1,
     >       /5X,' iterations=',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1322 FORMAT(5X, ' temperature='F10.3,5x,' beta=',E10.3,' seed=',I10)
 1323 FORMAT(5X, ' targeted accept.=',F10.3,5X,
     >           ' ddx=',F10.3,' c0=',F10.3,' c1=',F10.3,' c2=',F10.3)
 1324 FORMAT(5X, ' number atom changes = 'I10,
     >           ' change vector = <',3F10.3'>')
 1325 FORMAT(5X, ' number of groups = 'I10)
 1326 FORMAT(5X, ' pressure='E10.3,' (',F10.3,' atm)',
     >           ' ddv=',F10.3,' v0=',F10.3,' v1=',F10.3,' v2=',F10.3)
99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('vib_vib: error opening file: "hess_file"',911,
     &       DISK_ERR)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('vib_vib: error reading file: "hess_file"',911,
     &       DISK_ERR)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit
     & ('vib_vib: unexpected EOF when reading file: "hess_file"',911,
     &       DISK_ERR)

      end


c     **********************************************
c     *                                            *
c     *                md_xs_analysis              *
c     *                                            *
c     **********************************************
      subroutine md_xs_analysis(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      real*8  FEFF_RMIN
      parameter (FEFF_RMIN=1.50d0)

      integer taskid_i,taskid_j,npi,npj
      integer taskid,np,MASTER
      parameter (MASTER=0)

      logical oprint,value,runfeff,nohydrogen
      character*20 codeversion
      character*2 tag2
      character*4 tag4
      integer nframes,natoms,stride,natoms1,center,nkf,nr
      integer iii,i,k,icount,icpu,pto,pfrom,iw,first,last
      real*8  rmax,rcut,x,y,z,vx,vy,vz,time1,volume,unita(3,3),ft(3)
      real*8  xmin,xmax,ymin,ymax,rmin,e0,s0
      complex*16 zz

      integer nkf0
      parameter (nkf0 = 601)
      integer kf(2),chif(2),rion(2),tag(2),r(2),chi_r(2),tmpfft(2)

      character*30 spectroscopy,absorber,edge
      character*80 chi_filename
      character*255 header,full_chi_filename

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      real*8   md_xs_find_rmin
      external md_xs_find_rmin
       

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call Parallel2d_np_i(npi)
      call Parallel2d_np_j(npj)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

#ifdef FEFF
      call feff_codeversion(codeversion)

*     **** get feff options from rtdb ****
      if (.not.btdb_cget(rtdb,'md_xs:spectroscopy',1,spectroscopy)) 
     >   spectroscopy = 'exafs'
      if (.not.btdb_cget(rtdb,'md_xs:absorber',1,absorber))
     >   absorber = ''
      if (.not.btdb_cget(rtdb,'md_xs:edge',1,edge))
     >   edge = 'k'
      if (.not.btdb_get(rtdb,'md_xs:center',mt_int,1,center))
     >   center = 1
      if (.not.btdb_get(rtdb,'md_xs:rmax',mt_dbl,1,rmax))
     >   rmax = 10.0d0
      if (.not.btdb_get(rtdb,'md_xs:e0',mt_dbl,1,e0))
     >   e0 = 0.0d0
      if (.not.btdb_get(rtdb,'md_xs:s0',mt_dbl,1,s0))
     >   s0 = 1.0d0
      if (.not.btdb_get(rtdb,'md_xs:rcut',mt_dbl,1,rcut))
     >   rcut = 5.0d0
      if (rcut.lt.1.0d-3) rcut = 0.0d0
      rcut = rcut/0.529177d0
      call util_file_name_noprefix('md_xs',.true.,.true.,header)

      if (.not.btdb_get(rtdb,'md_xs:stride',mt_int,1,stride))
     >   stride = 1

      if (.not.btdb_get(rtdb,'md_xs:first',mt_int,1,first))
     >   first = 1

      if (.not.btdb_get(rtdb,'md_xs:last',mt_int,1,last))
     >   last = -1

      if (.not.btdb_get(rtdb,'md_xs:nohydrogen',mt_log,1,nohydrogen))
     >   nohydrogen = .true.

      if (.not.btdb_get(rtdb,'md_xs:Fourier',mt_dbl,3,ft)) then
         ft(1) = 2.3d0
         ft(2) = 17.0d0
         ft(3) = 1.0d0
      end if
      if (.not.btdb_get(rtdb,'md_xs:Fourier_iw',mt_int,1,iw))
     >   iw = 2

*     **** define the chi data filename ****
      if (.not.btdb_cget(rtdb,'md_xs:chi_filename',1,chi_filename))
     >    call util_file_prefix('chi',chi_filename)
      call util_file_name_noprefix(chi_filename,.false.,
     >                             .false.,
     >                             full_chi_filename)


      if (oprint) then
         write(luout,*) 
     >   "==== ","MD XS analysis - Distribution: ",codeversion," ===="
         write(luout,*)
         write(luout,*) "Copyright (c) [2002] University of Washington"
         write(luout,*) "see ./nwchem/src/nwpw/nwpwlib/md_xs/FEFF6L",
     >                  " LICENSE file for copying details"
          write(luout,*)
          write(luout,*)
         write(luout,*)
     >"************************************************************",
     >"**********"
          write(luout,*) "Distribution:  FEFF6L"
         write(luout,*)
     >"This software was prepared in part with US Government Funding",
     >" under"
      write(luout,*)
     >"DOE contract DE-FG03-97ER45623."
         write(luout,*)
     >"Redistribution and use of this Distribution in source and binary"
         write(luout,*)
     >"formats, with or without modification is permitted, provided the"
         write(luout,*)
     >"following conditions are met:"
         write(luout,*)
         write(luout,*)
     >"Redistributions must retain the above notices and the following",
     >" list"
         write(luout,*)
     >"of conditions and disclaimer "
         write(luout,*)
         write(luout,*)
     >"Modified versions carry the marking"
         write(luout,*)
     >"     Based on or developed using Distribution: FEFF6L"
         write(luout,*)
     >"     Copyright (c) [2002] University of Washington"
         write(luout,*)
         write(luout,*)
     >"Recipient acknowledges the right of the University of",
     >" Washington to"
         write(luout,*)
     >"prepare uses of this Distribution and its modifications that",
     >" may be"
         write(luout,*)
     >"substantially similar or functionally equivalent to"
         write(luout,*)
     >"Recipient-prepared modifications."
         write(luout,*)
         write(luout,*)
     >"Recipient and anyone obtaining access to the Distribution",
     >" through"
         write(luout,*)
     >"recipient's actions accept all risk associated with possession",
     >" and"
         write(luout,*)
     >"use of the Distribution."
         write(luout,*)
         write(luout,*)
     >"THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED"
         write(luout,*)
     >"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED",
     >" WARRANTIES OF"
         write(luout,*)
     >"MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE",
     >" DISCLAIMED."
         write(luout,*)
     >"IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS",
     >" TO THE"
         write(luout,*)
     >"DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,",
     >" SPECIAL,"
         write(luout,*)
     >"EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
     >" TO,"
         write(luout,*)
     >"PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,",
     >" OR"
         write(luout,*)
     >"REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY",
     >" THEORY OF"
         write(luout,*)
     >"LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT",
     >"  (INCLUDING"
         write(luout,*)
     >"NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE",
     >" OF THIS"
         write(luout,*)
     >"SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
         write(luout,*)
     >"************************************************************",
     >"**********"
         write(luout,*)
         write(luout,*)
         write(luout,*) 'MOTION filename:  ',trim(filename)
      end if

      if (taskid.eq.MASTER) then
         call MOTION_nframes(filename,nframes,natoms,volume,unita)
      end if
      call Parallel_Brdcst_ivalue(0,nframes)
      call Parallel_Brdcst_ivalue(0,natoms)
      call Parallel_Brdcst_value(0,volume)
      call Parallel_Brdcst_values(0,9,unita)
      if (last.eq.(-1)) last = nframes

      if (oprint) then
         write(luout,*) 'number of frames: ',nframes
         write(luout,*) 'first:            ',first
         write(luout,*) 'last:             ',last
         write(luout,*) 'stride:           ',stride
         write(luout,*) 'natoms:           ',natoms
         if (nohydrogen) then
            write(luout,*) 'No Hydrogens included.'
         else
            write(luout,*) 'Hydrogens included.'
         end if
         write(luout,1230) 
         write(luout,1241) unita(1,1),unita(2,1),unita(3,1)
         write(luout,1241) unita(1,2),unita(2,2),unita(3,2)
         write(luout,1241) unita(1,3),unita(2,3),unita(3,3)
         write(luout,1231) volume
         write(luout,*)
         write(luout,*) 'spectroscopy:     ',spectroscopy
         write(luout,*) 'absorber:         ',absorber
         write(luout,*) 'edge:             ',edge
         write(luout,*) 'center:           ',center
         write(luout,*) 'e0:               ',e0
         write(luout,*) 's0:               ',s0
         write(luout,*) 'Rmax (Angstroms): ',rmax
         write(luout,*) 'Rcut (Angstroms): ',rcut*0.529177d0
         write(luout,*)
         write(luout,*) 'Fourier Transform Window'
         write(luout,*) 'Kmin (Angstroms): ',ft(1)
         write(luout,*) 'Kmax (Angstroms): ',ft(2)
         write(luout,*) 'dK   (Angstroms): ',ft(3)
         write(luout,*) 'w               : ',iw
      end if


*     **** allocate memory from stack ****
      value =  BA_push_get(mt_dbl,nkf0,'kf',kf(2),kf(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nkf0,'chif',chif(2),chif(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*natoms,'rion',rion(2),rion(1)) 
      value = value.and.
     >        BA_push_get(mt_byte,2*natoms,'tag',tag(2),tag(1)) 
      if (.not.value) 
     >   call errquit("md_xs_analysis: out of stack",1,MA_ERR)


      if (taskid_i.eq.MASTER)
     >   open(unit=19,file=filename,form='formatted')

      call dcopy(nkf0,0.0d0,0,dbl_mb(kf(1)),1)
      call dcopy(nkf0,0.0d0,0,dbl_mb(chif(1)),1)
      icount = 0
      icpu   = 0
      !do k=1,nframes
      do k=1,last
         if (taskid_i.eq.MASTER) then
            read(19,*) time1,natoms1,volume,unita
            do i=1,natoms1
               tag2 = '  '
               read(19,*) iii,tag2,tag4,x,y,z,vx,vy,vz
               byte_mb(tag(1)+2*(i-1))   = tag2(1:1)
               byte_mb(tag(1)+2*(i-1)+1) = tag2(2:2)
               dbl_mb(rion(1)+3*(i-1))   = x
               dbl_mb(rion(1)+3*(i-1)+1) = y
               dbl_mb(rion(1)+3*(i-1)+2) = z
            end do
            rmin = md_xs_find_rmin(unita,natoms1,dbl_mb(rion(1)))
         end if
         call Parallela_Brdcst_value(1,MASTER,rmin)


         runfeff = ((k-first).ge.0)
         runfeff = runfeff.and.(mod(k-first,stride).eq.0)
         runfeff = runfeff.and.(rmin.gt.FEFF_RMIN)

         !if ((mod(k,stride).eq.0).and.(rmin.gt.FEFF_RMIN)) then
         if (runfeff) then
            if (taskid_i.eq.MASTER) write(*,*) "running feff, k=",k
            if (mod(icount,npj).eq.taskid_j) then
                 if (taskid_i.eq.MASTER) then
                    if (icpu.eq.MASTER) then
                       !*** compute md_xs with natoms1,tag,rion +--> kf,chif ***
                       call md_xs(header,spectroscopy,absorber,edge,
     >                         rmax,e0,s0,center,rcut,unita,
     >                         natoms1,byte_mb(tag(1)),dbl_mb(rion(1)),
     >                         nohydrogen,
     >                         nkf,nkf0,dbl_mb(kf(1)),dbl_mb(chif(1)))
                    else
                       !*** send natoms1,tag,rion to icpu ***
                       pto = Parallel2d_convert_taskid_i(icpu)
                       call Parallel_send_ivalues(pto,   4*k, 1,natoms1)
                       call Parallel_send_values(pto,    4*k+1,9,unita)
                       call Parallel_send_characters(pto,4*k+2,
     >                                        2*natoms1,byte_mb(tag(1)))
                       call Parallel_send_values(pto,    4*k+3,
     >                                        3*natoms1,dbl_mb(rion(1)))
                    end if
                 else
                    if (icpu.eq.taskid_i) then
                       !*** receive natoms1,tag,rion from MASTER ***
                       pfrom = Parallel2d_convert_taskid_i(MASTER)
                       call Parallel_recv_ivalues(pfrom,4*k,  1,natoms1)
                       call Parallel_recv_values(pfrom, 4*k+1,9,unita)
                       call Parallel_recv_characters(pfrom,4*k+2,
     >                                        2*natoms1,byte_mb(tag(1)))
                       call Parallel_recv_values(pfrom,    4*k+3,
     >                                        3*natoms1,dbl_mb(rion(1)))

                       !*** compute md_xs with natoms1,tag,rion +--> kf,chif ***
                       call md_xs(header,spectroscopy,absorber,edge,
     >                          rmax,e0,s0,center,rcut,unita,
     >                          natoms1,byte_mb(tag(1)),dbl_mb(rion(1)),
     >                          nohydrogen,
     >                          nkf,nkf0,dbl_mb(kf(1)),dbl_mb(chif(1)))
                    end if
                 end if
                 icpu = mod(icpu+1,npi)
            end if
            icount = icount + 1
         end if
      end do
      if (taskid_i.eq.MASTER) close(19)

      if (icount.gt.0) then
         x = 1.0d0/dble(icount)
      else
         x = 1.0d0
      end if
      call Parallel_Vector_SumAll(nkf0,dbl_mb(chif(1)))
      call dscal(nkf0,x,dbl_mb(chif(1)),1)

c     **** deallocate stack ****
      value =           BA_pop_stack(tag(2))
      value = value.and.BA_pop_stack(rion(2))
      if (.not.value)
     >   call errquit("md_xs_analysis: popping stack",2,MA_ERR)

      
*     ***************************
*     **** ascii plot of chi ****
*     ***************************
      if (taskid.eq.MASTER) then
         write(luout,*)
         write(luout,*)
         xmin =  999.0e12
         ymin =  999.0e12
         xmax = -999.0e12
         ymax = -999.0e12
         do i=1,nkf
            y = dbl_mb(chif(1)+i-1)
            if (y.lt.ymin) ymin = y
            if (y.gt.ymax) ymax = y
         end do
         iii = -9
         xmin = dbl_mb(kf(1))
         xmax = dbl_mb(kf(1)+nkf-1)
         do i=nkf,1,-1
            x = dbl_mb(kf(1)+i-1)
            y = dbl_mb(chif(1)+i-1)
            if ((dabs(y/ymax).gt.0.01d0).and.(iii.lt.0)) then
               iii = i
               xmax = x
            end if
         end do
         call util_ascii_setwindow(xmin,xmax,ymin,ymax)
         !write(luout,'(A,E10.3,A)') "chi_max=",ymax
         call util_ascii_plotter("",luout,iii,
     >                           dbl_mb(kf(1)),
     >                           dbl_mb(chif(1)),'*',
     >                           "Chi(K) - EXAFS Spectra",
     >                           "K (Ang-1)","Chi(K)")
         !write(luout,'(A,E10.3)') "chi_min=",ymin
         !write(luout,*)
         !write(luout,900) xmin,0.5d0*(xmax+xmin),xmax
         write(luout,*)
         write(luout,*)
      end if
c  900 format('kmin=',f6.1,'Ang-1',18x,f6.1,'Ang-1',
c     >                           17x,'kmax=',f6.1,'Ang-1')



*     *********************************
*     **** write out chi data file ****
*     *********************************

*     **** define the chi data filename ****
      if (.not.btdb_cget(rtdb,'md_xs:chi_filename',1,chi_filename))
     >    call util_file_prefix('chi',chi_filename)
      call util_file_name_noprefix(chi_filename,.false.,
     >                             .false.,
     >                             full_chi_filename)
      if (oprint) then
         write(luout,*) "The following file can be used for plotting:"
         write(luout,*) 'Chi_k data filename:',trim(full_chi_filename)
      end if

      if (taskid.eq.MASTER) then
         open(unit=19,file=full_chi_filename,form='formatted')
         write(19,'(A,A)') '#Chi_k data filename:',
     >                     trim(full_chi_filename)
         write(19,'(A,I4)') "#nkf =",nkf
         write(19,'(A,A)') "#spectroscopy:     ",trim(spectroscopy)
         write(19,'(A,A)') "#absorber:         ",trim(absorber)
         write(19,'(A,A)') "#edge:             ",trim(edge)
         write(19,'(A,I6)') "#center =",center
         write(19,'(A,F16.8)') "#rmax =",rmax
         write(19,'(A,F16.8)') "#e0 =",e0
         write(19,'(A,F16.8)') "#s0 =",s0
         write(19,'(A,A)') '# MOTION filename:  ',trim(filename)
         write(19,'(A,I8)') '# number of frames =',nframes
         write(19,'(A,I8)') '# number of frames used, icount =',icount
         write(19,'(A,I8)') '# stride =',stride
         write(19,'(A,I8)') '# natoms =',natoms
         write(19,'(A,F16.8)') '#volume =',volume
         write(19,'(A,F16.8)') '#rcut =',rcut*0.529177d0
         write(19,'(A,1x,A14,4A20)') '#','k(Ang**-1)','chi','k*chi',
     >                               'k**2*chi','k**3*chi'
         do i=1,nkf
            x = dbl_mb(kf(1)+i-1)
            y = dbl_mb(chif(1)+i-1)
            write(19,'(F16.4,4E20.9)') x,y,x*y,x*x*y,x*x*x*y
         end do
         close(19)

      end if


c     **** allocate stack - fft  ****
      !nr = 2*nkf
      nr = 2048
      value = BA_push_get(mt_dbl,(8*nr+30),'tmpfft',
     >                     tmpfft(2),tmpfft(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nr),'chi_r',chi_r(2),chi_r(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nr),'rr_md_xs',r(2),r(1))
      if (.not.value) 
     >   call errquit("md_xs_analysis: out of stack",3,MA_ERR)
      
      call dcffti(nr,dbl_mb(tmpfft(1)))

*     ***************************************
*     **** compute the fourier transform ****
*     ***************************************
      call md_xs_fourier(iw,ft(1),ft(2),ft(3),
     >                    nkf,dbl_mb(kf(1)),dbl_mb(chif(1)),
     >                    nr, dbl_mb(r(1)),dcpl_mb(chi_r(1)),
     >                    dbl_mb(tmpfft(1)))

*     *****************************
*     **** ascii plot of chi_r ****
*     *****************************
      if (taskid.eq.MASTER) then
         write(luout,*)
         write(luout,*)
         ymin = 0.0d0
         ymax = -999.0e12
         do i=1,nr/2
            zz = dcpl_mb(chi_r(1)+i-1)
            y  = dsqrt(dble(zz)**2 + dimag(zz)**2)
            dbl_mb(tmpfft(1)+i-1) = y
            if (y.gt.ymax) ymax = y
         end do
         iii = -9
         xmin = dbl_mb(r(1))
         xmax = dbl_mb(r(1)+nr/2-1)
         do i=nr/2,1,-1
            x  = dbl_mb(r(1)+i-1)
            zz = dcpl_mb(chi_r(1)+i-1)
            y  = dsqrt(dble(zz)**2 + dimag(zz)**2)
            if ((dabs(y/ymax).gt.0.01d0).and.(iii.lt.0)) then
               iii = i
               xmax = x
            end if
         end do
         call util_ascii_setwindow(xmin,xmax,ymin,ymax)
         !write(luout,'(A,E10.3,A)') "chi_r_max=",ymax
         call util_ascii_plotter("",luout,iii,
     >                           dbl_mb(r(1)),
     >                           dbl_mb(tmpfft(1)),'*',
     >                           "|FFT(K^w*Chi(K))| - EXAFS Spectra",
     >                           "R (Ang)", "|Chi(R)|")
         !write(luout,'(A,E10.3)') "chi_r_min=",ymin
         !write(luout,*)
         !write(luout,901) xmin,0.5d0*(xmax+xmin),xmax
         write(luout,*)
         write(luout,*)
      end if
c  901 format('rmin=',f6.1,'Ang',18x,f6.1,'Ang',
c     >                           17x,'rmax=',f6.1,'Ang')

*     ************************************
*     **** write out chi_r data file ****
*     ************************************

*     **** define the chi_r data filename ****
      if (.not.btdb_cget(rtdb,'md_xs:chi_r_filename',1,chi_filename))
     >    call util_file_prefix('chi_r',chi_filename)
      call util_file_name_noprefix(chi_filename,.false.,
     >                             .false.,
     >                             full_chi_filename)
      if (oprint) then
         write(luout,*) 'The following file can be used for plotting:'
         write(luout,*) 'Chi_r data filename:',trim(full_chi_filename)
      end if

      if (taskid.eq.MASTER) then
         open(unit=19,file=full_chi_filename,form='formatted')
         write(19,'(A,A)') '#Chi_r data filename:',
     >                     trim(full_chi_filename)
         write(19,'(A,I4)') "#nr  =",nr
         write(19,'(A,I4)') "#iw  =",iw
         write(19,'(A,F16.8)') "#Kmin =",ft(1)
         write(19,'(A,F16.8)') "#Kmax =",ft(2)
         write(19,'(A,F16.8)') "#dK =",ft(3)
         write(19,'(A,A)') "#spectroscopy:     ",trim(spectroscopy)
         write(19,'(A,A)') "#absorber:         ",trim(absorber)
         write(19,'(A,A)') "#edge:             ",trim(edge)
         write(19,'(A,I6)') "#center =",center
         write(19,'(A,F16.8)') "#rmax =",rmax
         write(19,'(A,F16.8)') "#e0 =",e0
         write(19,'(A,F16.8)') "#s0 =",s0
         write(19,'(A,A)') '# MOTION filename:  ',trim(filename)
         write(19,'(A,I8)') '# number of frames =',nframes
         write(19,'(A,I8)') '# number of frames used, icount =',icount
         write(19,'(A,I8)') '# stride =',stride
         write(19,'(A,I8)') '# natoms =',natoms
         write(19,'(A,F16.8)') '#volume =',volume
         write(19,'(A,F16.8)') '#rcut =',rcut*0.529177d0
         write(19,'(A,1x,A14,3A20)') '#','r(Ang)','|chi_r|',
     >                               'real(chi_r)','cmplx(chi_r)'
         do i=1,nr/2
            x  = dbl_mb(r(1)+i-1)
            zz = dcpl_mb(chi_r(1)+i-1)
            y  = dsqrt(dble(zz)**2 + dimag(zz)**2)
            write(19,'(F16.4,4E20.9)') x,y,zz
         end do
         close(19)
      end if



c     **** deallocate stack ****
      value =           BA_pop_stack(r(2))
      value = value.and.BA_pop_stack(chi_r(2))
      value = value.and.BA_pop_stack(tmpfft(2))
      value = value.and.BA_pop_stack(chif(2))
      value = value.and.BA_pop_stack(kf(2))
      if (.not.value)
     >   call errquit("md_xs_analysis: popping stack",4,MA_ERR)

#else
      if (oprint) then
         write(luout,*) "==== ",
     >     "FEFF library not inlucded in compile: ",
     >     "MD XS analysis not performed"," ===="
      end if
#endif
      return
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F12.1)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
      end

c     **********************************************
c     *                                            *
c     *                md_xs_fourier               *
c     *                                            *
c     **********************************************
      subroutine md_xs_fourier(iw,wmin,wmax,dw,
     >                         nk,k,chi_k,
     >                         nr,r,chi_r,tmpfft)
      integer iw
      real*8 wmin,wmax,dw
      integer nk
      real*8 k(*),chi_k(*)
      integer nr
      real*8 r(*)
      complex*16 chi_r(*)
      real*8     tmpfft(*)
      !complex*16 tmpfft(*)

*     **** local variables ****
      integer i
      real*8 x,y,h,pi
      real*8 dr,dk,scal

      
      pi = 4.0d0*datan(1.0d0)
      dk = k(2)-k(1)
      dr = pi/(dk*dble(nr))
      scal = dk/dsqrt(pi*nr)

      call dcopy(2*nr,0.0d0,0,chi_r,1)

      do i=1,nk
        x = k(i)
        y = chi_k(i)

        !**** Hanning window ****
        if(x.le.(wmin-0.5d0*dw)) then
           h = 0.0d0
        else if ((x.ge.(wmin-0.5d0*dw)).and.(x.lt.(wmin+0.5d0*dw))) then
           h = dsin(0.5d0*pi*(x-wmin+0.5*dw)/dw)**2
        else if ((x.ge.(wmin+0.5d0*dw)).and.(x.le.(wmax-0.5d0*dw))) then
           h = 1.0d0
        else if ((x.gt.(wmax-0.5d0*dw)).and.(x.le.(wmax+0.5d0*dw))) then
           h = dcos(0.5d0*pi*(x-wmax+0.5*dw)/dw)**2
        else if (x.gt.(wmax+0.5d0*dw)) then
           h = 0.0d0
        else
           h = 0.0d0
        end if

        chi_r(i) = dcmplx(h*y*(x**iw),0.0d0)
      end do

      call dcfftb(nr,chi_r,tmpfft)

      do i=1,nr
         r(i)     = (i-1)*dr
         chi_r(i) = scal*chi_r(i)
      end do

      return
      end 

      
        


#ifdef FEFF
c     **********************************************
c     *                                            *
c     *                md_xs                       *
c     *                                            *
c     **********************************************
      subroutine md_xs(header,spectroscopy,absorber,edge,
     >                 rmax,e0,s0,center,rcut,unita,nion,symb,rion,
     >                 nohydrogen,
     >                 nkf,nkf0,kf,chi)
      implicit none
      character*(*) header,spectroscopy,absorber,edge
      real*8  rmax,e0,s0
      integer center
      real*8 rcut
      real*8 unita(3,3)
      integer nion
      character*2 symb(*)
      real*8      rion(3,*)
      logical nohydrogen
      integer nkf,nkf0
      real*8 kf(*)
      real*8 chi(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid,MASTER,np
      parameter (MASTER=0)

      logical value,found
      integer nkf4,j,ii,ia,zi
      integer zion(2),katm(2),zkatm(2),nkatm
      integer zion4(2),katm4(2),zkatm4(2),tchi(2)
      integer rion2(2),indx2(2),nion2,center2
      character*2 item

      character*2 symbols(112)
      data symbols/
     $     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     $     'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     $     'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     $     'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     $     'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     $     'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     $     'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     $     'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     $     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     $     'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/

*     **** external functions ****
      integer  md_xs_size_cluster
      external md_xs_size_cluster


      call Parallel_taskid(taskid)
      call Parallel_np(np)
      
      nion2 = md_xs_size_cluster(center,rcut,unita,nion,rion)

*     ***** allocate stack ****
      value = BA_push_get(mt_int,nion2,'md_xs:zion',zion(2),zion(1))
      value = value.and.
     >        BA_push_get(mt_int,nion2,'md_xs:zkatm',zkatm(2),zkatm(1))
      value = value.and.
     >        BA_push_get(mt_int,nion2,'md_xs:katm',katm(2),katm(1))
      value = value.and.
     >        BA_push_get(mt_int,nion2,'md_xs:zion4',zion4(2),zion4(1))
      value = value.and.
     >      BA_push_get(mt_int,nion2,'md_xs:zkatm4',zkatm4(2),zkatm4(1))
      value = value.and.
     >        BA_push_get(mt_int,nion2,'md_xs:katm4',katm4(2),katm4(1))
      value = value.and.
     >        BA_push_get(mt_int,nion2,'md_xs:indx2',indx2(2),indx2(1))
      value = value.and.
     >       BA_push_get(mt_dbl,3*nion2,'md_xs:rion2',rion2(2),rion2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nkf0,'md_xs:tchi',tchi(2),tchi(1))
      if (.not.value)
     >   call errquit("md_xs: pushing stack",0,MA_ERR)
      call dcopy(nkf0,0.0d0,0,dbl_mb(tchi(1)),1)

      call md_xs_gen_cluster(center,rcut,unita,nion,rion,
     >                       nion2,dbl_mb(rion2(1)),
     >                       int_mb(indx2(1)),center2)

      


      do ii=1,nion2
         item = '  '
         item = symb(int_mb(indx2(1)+ii-1))
         zi = -1
         do j=1,112
            if (item.eq.symbols(j)) zi = j
         end do
         int_mb(zion(1)+ii-1) = zi
      end do

      nkatm = 0
      do ii=1,nion2
         found = .false.
         do j=1,nkatm
            if (int_mb(zion(1)+ii-1).eq.int_mb(zkatm(1)+j-1)) then
               found = .true.
               ia = j
            end if
         end do
         if (found) then
            int_mb(katm(1)+ii-1) = ia
         else
            nkatm = nkatm + 1
            int_mb(zkatm(1)+nkatm-1) = int_mb(zion(1)+ii-1)
            int_mb(katm(1)+ii-1) = nkatm
         end if
      end do



#ifdef FEFF8
      call icopy(nion2,int_mb(zion(1)),1,int_mb(zion4(1)),1)
      call icopy(nion2,int_mb(katm(1)),1,int_mb(katm4(1)),1)
      call icopy(nkatm,int_mb(zkatm(1)),1,int_mb(zkatm4(1)),1)
#else
      call ftupid_icopy8to4(nion2,int_mb(zion(1)),int_mb(zion4(1)))
      call ftupid_icopy8to4(nion2,int_mb(katm(1)),int_mb(katm4(1)))
      call ftupid_icopy8to4(nkatm,int_mb(zkatm(1)),int_mb(zkatm4(1)))
#endif
      call feff_fortran(header,spectroscopy,absorber,edge,
     >                  center2,rmax,e0,s0,
     >                  nkatm,int_mb(katm4(1)),int_mb(zkatm4(1)),
     >                  nion2,int_mb(zion4(1)),dbl_mb(rion2(1)),
     >                  nohydrogen,
     >                  nkf4,kf,dbl_mb(tchi(1)))
#ifdef FEFF8
      nkf = nkf4
#else
      call ftupid_icopy4to8(1,nkf4,nkf)
#endif

*     **** chi += tchi  ****
      call daxpy(nkf0,1.0d0,dbl_mb(tchi(1)),1,chi,1)


*     **** popping stack ****
      value =           BA_pop_stack(tchi(2))
      value = value.and.BA_pop_stack(rion2(2))
      value = value.and.BA_pop_stack(indx2(2))
      value = value.and.BA_pop_stack(katm4(2))
      value = value.and.BA_pop_stack(zkatm4(2))
      value = value.and.BA_pop_stack(zion4(2))
      value = value.and.BA_pop_stack(katm(2))
      value = value.and.BA_pop_stack(zkatm(2))
      value = value.and.BA_pop_stack(zion(2))
      if (.not.value)
     >   call errquit("md_xs: popping stack",1,MA_ERR)
     
      return
      end


cc     **********************************************
cc     *                                            *
cc     *                md_xs_create_json           *
cc     *                                            *
cc     **********************************************
c      subroutine md_xs_create_json(rtdb,header,nion,symbols,rion,json)
c      implicit none
c      integer rtdb
c      character*(*) header
c      integer nion
c      character*2 symbols(*)
c      real*8 rion(3,*)
c      character*(*) json
c#include "bafdecls.fh"
c#include "btdb.fh"
c#include "errquit.fh"
c#include "stdio.fh"
c
c      integer ii,center
c      real*8 rmax
c      character*30 tstr1,tstr2,tstr3
c
c*     **** open json ****
c      json = "{"
c
c*     **** "spectroscopy": ****
c      if (btdb_cget(rtdb,'md_xs:spectroscopy',1,tstr1)) then
c         json = trim(json)//"""spectroscopy"":"
c         json = trim(json)//" "//""""//trim(tstr1)//""","
c      end if
c
c*     **** "absorber": ****
c      if (btdb_cget(rtdb,'md_xs:absorber',1,tstr1)) then
c         json = trim(json)//" ""absorber"": ["
c         json = trim(json)//""""//trim(tstr1)//"""],"
c      end if
c
c*     **** "edge": ****
c      if (btdb_cget(rtdb,'md_xs:edge',1,tstr1)) then
c         json = trim(json)//" ""edge"":"
c         json = trim(json)//" "//""""//trim(tstr1)//""","
c      end if
c
c*     **** "center": ****
c      if (btdb_get(rtdb,'md_xs:center',mt_int,1,center)) then
c         write(tstr1,'(I6)') center
c         json = trim(json)//" ""center"": ["
c         json = trim(json)//trim(adjustl(tstr1))//"],"
c      end if
c
c*     **** "rmax": ****
c      if (btdb_get(rtdb,'md_xs:rmax',mt_dbl,1,rmax)) then
c         write(tstr1,'(F20.6)') rmax
c         json = trim(json)//" ""rmax"":"
c         json = trim(json)//" "//trim(adjustl(tstr1))//","
c      end if
c
c*     **** "geometry": ****
c      json = trim(json)//" ""geometry"": ["
c      do ii=1,nion
c         write(tstr1,'(E22.15)') rion(1,ii)
c         write(tstr2,'(E22.15)') rion(2,ii)
c         write(tstr3,'(E22.15)') rion(3,ii)
c         if (ii.eq.1) then
c            json = trim(json)//"["""//trim(symbols(ii))//""","
c         else
c            json = trim(json)//" ["""//trim(symbols(ii))//""","
c         end if
c         json = trim(json)//" "//trim(adjustl(tstr1))//","
c         json = trim(json)//" "//trim(adjustl(tstr2))//","
c         json = trim(json)//" "//trim(adjustl(tstr3))//"]"
c         if (ii.lt.nion) json = trim(json)//","
c      end do
c      json = trim(json)//"], "
c
c*     **** "scratch_dir": ****
c      json = trim(json)//" ""scratch_dir"":"
c      json = trim(json)//" "//""""//trim(header)//""""
c
c*     **** close json ****
c      json = trim(json)//"}"
c
c      return
c      end



*     **********************************************   
*     *                                            *
*     *              md_xs_gen_cluster             *
*     *                                            *
*     **********************************************   
*
*     The routine generates a cluster spanning rcut from the center, rion(*,center).
*
*     Entry - center: location of center
*             rcut: cluster radius
*             unita: lattice vectors
*             nion: number of ions in geometry
*             rion: ions in geometry
*     Exit  - nion2: number of ions in cluster
*             rion: ions in cluster
*
      subroutine md_xs_gen_cluster(center,rcut,unita,nion,rion,
     >                             nion2,rion2,indx2,center2)
      implicit none
      integer center
      real*8  rcut
      real*8  unita(3,3)
      integer nion
      real*8  rion(3,*)
      integer nion2
      real*8  rion2(3,*)
      integer indx2(*)
      integer center2

*     **** local variables ****
      integer i,j,k,ii,n,n1,n2
      real*8 d,dmin,x,y,z,xc,yc,zc,vol


*     **** cluster center ****
      xc = rion(1,center)
      yc = rion(2,center)
      zc = rion(3,center)


*     **** if rcut<=0 or unit cell volume is small then do a straight copy ****
      vol = unita(1,1)*(unita(2,2)*unita(3,3)-unita(3,2)*unita(2,3))
     >    + unita(2,1)*(unita(3,2)*unita(1,3)-unita(1,2)*unita(3,3))          
     >    + unita(3,1)*(unita(1,2)*unita(2,3)-unita(2,2)*unita(1,3))
      vol = dabs(vol)

      if ((rcut.le.0.0d0).or.(vol.lt.2.0d0)) then
         nion2 = nion
         center2 = 1
         rion2(1,1) = rion(1,center)-xc
         rion2(2,1) = rion(2,center)-yc
         rion2(3,1) = rion(3,center)-zc
         indx2(1)   = center
         do ii=1,nion
            if (ii.ne.center) then
               indx2(ii)   = ii
               rion2(1,ii) = rion(1,ii)-xc
               rion2(2,ii) = rion(2,ii)-yc
               rion2(3,ii) = rion(3,ii)-zc
            end if
         end do
         return
      end if

      dmin = 999999999.999d0
      do i=1,3
         d = dsqrt(unita(1,i)**2 + unita(2,i)**2 + unita(3,i)**2)
         if (d.lt.dmin) dmin = d
      end do
      n = nint((rcut/dmin)+0.5d0)
      if (n.lt.1) n = 1 
      n1 = -n
      n2 =  n


*     **** make center first atom postion in rion2 ****
      nion2 = 0
      nion2 = nion2 + 1
      rion2(1,nion2) = rion(1,center)-xc
      rion2(2,nion2) = rion(2,center)-yc
      rion2(3,nion2) = rion(3,center)-zc
      indx2(nion2)   = center
      center2 = 1

      do k=n1,n2
      do j=n1,n2
      do i=n1,n2
         do ii=1,nion
            x = rion(1,ii) + i*unita(1,1) + j*unita(1,2) + k*unita(1,3)
            y = rion(2,ii) + i*unita(2,1) + j*unita(2,2) + k*unita(2,3)
            z = rion(3,ii) + i*unita(3,1) + j*unita(3,2) + k*unita(3,3)
            d = dsqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)

            !*** d < rcut and ignore center ***
            if ((d.le.rcut).and.(d.gt.1.0d-3)) then
               nion2 = nion2 + 1
               rion2(1,nion2) = x-xc
               rion2(2,nion2) = y-yc
               rion2(3,nion2) = z-zc
               indx2(nion2)   = ii
            end if
         end do
      end do
      end do
      end do

      return
      end


     
*     **********************************************
*     *                                            *
*     *              md_xs_size_cluster            *
*     *                                            *
*     **********************************************
*
*     The routine generates a cluster spanning rcut from the center, rion(*,center).
*
*     Entry - center: location of center
*             rcut: cluster radius
*             unita: lattice vectors
*             nion: number of ions in geometry
*             rion: ions in geometry
*     Exit  - nion2: number of ions in cluster
*
      integer function md_xs_size_cluster(center,rcut,unita,nion,rion)
      implicit none
      integer center
      real*8  rcut
      real*8  unita(3,3)
      integer nion
      real*8  rion(3,*)

*     **** local variables ****
      integer nion2
      integer i,j,k,ii,n,n1,n2
      real*8 d,dmin,x,y,z,xc,yc,zc,vol


      vol = unita(1,1)*(unita(2,2)*unita(3,3)-unita(3,2)*unita(2,3))
     >    + unita(2,1)*(unita(3,2)*unita(1,3)-unita(1,2)*unita(3,3))          
     >    + unita(3,1)*(unita(1,2)*unita(2,3)-unita(2,2)*unita(1,3))
      vol = dabs(vol)

      if ((rcut.le.0.0d0).or.(vol.lt.2.0d0)) then
         md_xs_size_cluster = nion
         return
      end if

      dmin = 999999999.999d0
      do i=1,3
         d = dsqrt(unita(1,i)**2 + unita(2,i)**2 + unita(3,i)**2)
         if (d.lt.dmin) dmin = d
      end do
      n = nint((rcut/dmin)+0.5d0)
      if (n.lt.1) n = 1
      n1 = -n
      n2 =  n

      nion2 = 0
      xc = rion(1,center)
      yc = rion(2,center)
      zc = rion(3,center)
      do k=n1,n2
      do j=n1,n2
      do i=n1,n2
         do ii=1,nion
            x = rion(1,ii) + i*unita(1,1) + j*unita(1,2) + k*unita(1,3)
            y = rion(2,ii) + i*unita(2,1) + j*unita(2,2) + k*unita(2,3)
            z = rion(3,ii) + i*unita(3,1) + j*unita(3,2) + k*unita(3,3)
            d = dsqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
            if (d.le.rcut) nion2 = nion2 + 1
         end do
      end do
      end do
      end do

      md_xs_size_cluster = nion2
      return
      end


      subroutine ftupid_icopy8to4(n,array8,array4)
      implicit none
      integer n
      integer*8 array8(*)
      integer*4 array4(*)
      integer i
      do i=1,n
          array4(i) = array8(i)
      end do
      return
      end

      subroutine ftupid_icopy4to8(n,array4,array8)
      implicit none
      integer n
      integer*4 array4(*)
      integer*8 array8(*)
      integer i
      do i=1,n
          array8(i) = array4(i)
      end do
      return
      end

#endif

      subroutine feff_copyright(luout,codeversion)
      implicit none
      integer luout
      character*(*) codeversion

      write(luout,*)
     >"==== ","MD XS analysis - Distribution: ",codeversion," ===="
      write(luout,*)
      write(luout,*) "Copyright (c) [2002] University of Washington"
      write(luout,*) "see ./nwchem/src/nwpw/nwpwlib/md_xs/FEFF6L",
     >                  " LICENSE file for copying details"
      write(luout,*)
      write(luout,*)
      write(luout,*)
     >"************************************************************",
     >"**********"
      write(luout,*) "Distribution:  FEFF6L"
      write(luout,*)
     >"This software was prepared in part with US Government Funding",
     >" under"
      write(luout,*)
     >"DOE contract DE-FG03-97ER45623."
      write(luout,*)
     >"Redistribution and use of this Distribution in source and binary"
      write(luout,*)
     >"formats, with or without modification is permitted, provided the"
      write(luout,*)
     >"following conditions are met:"
      write(luout,*)
      write(luout,*)
     >"Redistributions must retain the above notices and the following",
     >" list"
      write(luout,*)
     >"of conditions and disclaimer "
      write(luout,*)
      write(luout,*)
     >"Modified versions carry the marking"
      write(luout,*)
     >"     Based on or developed using Distribution: FEFF6L"
      write(luout,*)
     >"     Copyright (c) [2002] University of Washington"
      write(luout,*)
      write(luout,*)
     >"Recipient acknowledges the right of the University of",
     >" Washington to"
      write(luout,*)
     >"prepare uses of this Distribution and its modifications that",
     >" may be"
      write(luout,*)
     >"substantially similar or functionally equivalent to"
      write(luout,*)
     >"Recipient-prepared modifications."
      write(luout,*)
      write(luout,*)
     >"Recipient and anyone obtaining access to the Distribution",
     >" through"
      write(luout,*)
     >"recipient's actions accept all risk associated with possession",
     >" and"
      write(luout,*)
     >"use of the Distribution."
      write(luout,*)
      write(luout,*)
     >"THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED"
      write(luout,*)
     >"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED",
     >" WARRANTIES OF"
      write(luout,*)
     >"MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE",
     >" DISCLAIMED."
      write(luout,*)
     >"IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS",
     >" TO THE"
      write(luout,*)
     >"DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,",
     >" SPECIAL,"
      write(luout,*)
     >"EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
     >" TO,"
      write(luout,*)
     >"PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,",
     >" OR"
      write(luout,*)
     >"REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY",
     >" THEORY OF"
      write(luout,*)
     >"LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT",
     >"  (INCLUDING"
      write(luout,*)
     >"NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE",
     >" OF THIS"
      write(luout,*)
     >"SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
      write(luout,*)
     >"************************************************************",
     >"**********"
      write(luout,*)
      write(luout,*)

      return
      end


*     *****************************************
*     *                                       *
*     *             md_xd_hess_energy         *
*     *                                       *
*     *****************************************

      real*8 function md_xs_hess_energy(N,A,dq)
      implicit none
      integer N
      real*8 A(N,N),dq(N)

*     **** local variables ****
      integer i,j
      real*8 energy

      energy = 0.0d0
      do j=1,N
         do i=1,N
            energy = energy + dq(i)*A(i,j)*dq(j)
         end do
      end do
      md_xs_hess_energy = 0.5d0*energy
      return
      end


      subroutine md_xs_dsquar(A,S,N,MS)
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


*     *****************************************
*     *                                       *
*     *             md_xs_find_rmin           *
*     *                                       *
*     *****************************************
      real*8 function md_xs_find_rmin(unita,nion,rion)
      implicit none
      real*8  unita(3,3)
      integer nion
      real*8  rion(3,nion)

*     **** local variables ****
      integer ii,jj,i,j,k,n1,n2
      real*8  dx,dy,dz,x,y,z,d
      real*8  rmin

      rmin = 999.999d9
      n1 = -1
      n2 =  1
      if ((dabs(unita(1,1)-1.0d0).lt.1.0d-6) .and.
     >    (dabs(unita(2,1))      .lt.1.0d-6) .and.
     >    (dabs(unita(3,1))      .lt.1.0d-6) .and.
     >    (dabs(unita(1,2))      .lt.1.0d-6) .and.
     >    (dabs(unita(2,2)-1.0d0).lt.1.0d-6) .and.
     >    (dabs(unita(3,2))      .lt.1.0d-6) .and.
     >    (dabs(unita(1,3))      .lt.1.0d-6) .and.
     >    (dabs(unita(2,3))      .lt.1.0d-6) .and.
     >    (dabs(unita(3,3)-1.0d0).lt.1.0d-6) ) then
         n1 = 0
         n2 = 0
      end if

      do jj=1,nion-1
         do ii=jj+1,nion
            dx = rion(1,ii)-rion(1,jj)
            dy = rion(2,ii)-rion(2,jj)
            dz = rion(3,ii)-rion(3,jj)
            do k=n1,n2
            do j=n1,n2
            do i=n1,n2
               x = dx + i*unita(1,1) + j*unita(1,2) + k*unita(1,3)
               y = dy + i*unita(2,1) + j*unita(2,2) + k*unita(2,3)
               z = dz + i*unita(3,1) + j*unita(3,2) + k*unita(3,3)
               d = dsqrt(x**2 + y**2 + z**2)
               if (d.lt.rmin) rmin = d
            end do
            end do
            end do
         end do
      end do

      md_xs_find_rmin = rmin
      return
      end


