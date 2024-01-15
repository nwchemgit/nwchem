      subroutine qmd_ext_init(nat,dt_nucl,kbt,ext_forces)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util_params.fh"
#include "stdio.fh"
#include "qmd_ext.fh"
      integer nat              ![in]
      double precision dt_nucl ![in]
      double precision kbt     ![in]
      character*(*) ext_forces
      logical input_exists
c     pass plumed initialization bits
#ifdef USE_PLUMED
      integer ga_comm
      integer*4 ga_comm4
      integer*4 nat4
      if(ext_forces(1:inp_strlen(ext_forces)).ne.'plumed') call
     E     errquit(' qmd_ext works only with plumed ',0,0)
      call ga_mpi_comm_pgroup_default(ga_comm)
      if(ga_nodeid().eq.0) then
      call plumed_f_gcreate()
      call plumed_f_gcmd("setRealPrecision"//char(0),
     c     MA_sizeof(MT_DBL,1,MT_BYTE))
c we have no barostat/virial in qmd      
      call plumed_f_gcmd("setNoVirial"//char(0),0)
c conversion factor between the energy unit used in nwchem and kJ mol-1      
      call plumed_f_gcmd("setMDEnergyUnits"//char(0),cau2kj)
c conversion factor between the length unit used in nwchem and nm      
      call plumed_f_gcmd("setMDLengthUnits"//char(0),cau2nm)
c conversion factor between the time unit used in nwchem and ps      
      call plumed_f_gcmd("setMDTimeUnits"//char(0),cau2fs/1000)
c     kBT
      call plumed_f_gcmd("setKbT"//char(0),kbt)

c     get input file name
      call util_file_name('plumed.dat',.false.,.false.,extfrc_input)
c
c    check presence of plumed input file
c
      inquire(file=extfrc_input(1:inp_strlen(extfrc_input)),
     E     exist=input_exists)
      if(.not.input_exists) then
         write(luout,*) ' Plumed input file ',
     E        extfrc_input(1:inp_strlen(extfrc_input)),
     E ' NOT found '
         call util_file_name_noprefix('plumed.dat',.false.,.false.,
     E        extfrc_input)
         inquire(file=extfrc_input(1:inp_strlen(extfrc_input)),
     E        exist=input_exists)
         if(.not.input_exists) then
            write(luout,*) ' Plumed input file ',
     E           extfrc_input(1:inp_strlen(extfrc_input)),
     E           ' NOT found '
            call errquit(' qmd_ext: missing plumed input',0,0)
         endif
      endif
c     get output file name
      call util_file_name('plumed.out',.false.,.false.,extfrc_output)

      call plumed_f_gcmd("setPlumedDat"//char(0),
     E     extfrc_input(1:inp_strlen(extfrc_input))//char(0))
      call plumed_f_gcmd("setLogFile"//char(0),
     E     extfrc_output(1:inp_strlen(extfrc_output))//char(0))
      nat4=nat
      call plumed_f_gcmd("setNatoms"//char(0),nat4)
      call plumed_f_gcmd("setMDEngine"//char(0),"nwchem")
#ifdef PLUMED_HASMPI 
      call plumed_f_gcmd("setMPIFComm"//char(0),ga_comm4)
#endif
      call plumed_f_gcmd("setTimestep"//char(0),dt_nucl)
      write(luout,*) ' PPPPPPPPPPPPPPPPPPPPPPP '
      write(luout,*) ' NWChem PLUMED Interface '
      write(luout,*) ' PLUMED: input file ',
     I     extfrc_input(1:inp_strlen(extfrc_input))
      write(luout,*) ' PLUMED: output file ',
     I     extfrc_output(1:inp_strlen(extfrc_output))
      call plumed_f_gcmd("init"//char(0),0)
      write(luout,*) ' PLUMED: init done'
      write(luout,*) ' PPPPPPPPPPPPPPPPPPPPPPP '
      endif
#endif      
      return
      end

      subroutine qmd_ext_gradient(nat,istep_nucl,energy,m,g,r,
     E     ext_forces)
      implicit none
#include "global.fh"      
#include "mafdecls.fh"
#include "stdio.fh"
#include "inp.fh"
      integer nat                ! [in]
      integer istep_nucl         ! [in]
      double precision energy    ! [in]
      double precision g(3,nat)  ! gradients [in/out]
      double precision r(3,nat)  ! coordinates [in]
      double precision m(nat)    ! masses [in]
      character*(*) ext_forces
#ifdef USE_PLUMED
      integer i
c     call plumed: r[in], g[out]
c     pointer to current timestep
      if(ext_forces(1:inp_strlen(ext_forces)).ne.'plumed') call
     E     errquit(' qmd_ext works only with plumed ',0,0)
      if(ga_nodeid().eq.0) then
      call plumed_f_gcmd("setStep"//char(0),istep_nucl)
      call plumed_f_gcmd("setMasses"//char(0),m)
c plumed uses forces not gradientes. double flip
c first flip
        call dscal(3*nat, -1d0, g, 1)
       if(ga_nodeid().eq.0) then
          write(luout,*) '### forces before plumed ##'
          do i=1,nat
             write(luout,1) i,g(1,i),g(2,i),g(3,i)
          enddo
       endif
      call plumed_f_gcmd("setForces"//char(0),g)
      call plumed_f_gcmd("setPositions"//char(0),r)

c      call plumed_f_gcmd("setBox"//char(0),at_plumed)
c      call plumed_f_gcmd("setVirial"//char(0),virial)
c      write(6,*) ' E ext', energy
       call plumed_f_gcmd("setEnergy"//char(0),energy)
c     forces before plumed call
       call plumed_f_gcmd("calc"//char(0),0)
c second sign flip
       call dscal(3*nat, -1d0, g, 1)
      endif
c     broadcast from 0
      call ga_brdcst(2020,g,MA_sizeof(MT_DBL,3*nat,MT_BYTE),0)
       if(ga_nodeid().eq.0) then
          write(luout,*) '### forces after plumed ##'
          do i=1,nat
             write(luout,1) i,g(1,i),g(2,i),g(3,i)
          enddo
       endif
 1     format(i5,3f12.6)
c     check sign agreement with plumed ...
#endif

      return
      end
      subroutine qmd_ext_final(ext_forces)
      implicit none
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "qmd_ext.fh"
      character*(*) ext_forces
#ifdef USE_PLUMED
      integer myunit,io_status
      character*1 mychar
      if(ext_forces(1:inp_strlen(ext_forces)).ne.'plumed') call
     E     errquit(' qmd_ext works only with plumed ',0,0)
      if(ga_nodeid().eq.0) then
c     call plumed finalization bits
         call plumed_f_gfinalize()
c
c     copy plumed output to stdout
c
         myunit=2020
c         write(6,*) 'reading ',
c     E        extfrc_output(1:inp_strlen(extfrc_output))
         open(myunit, file=extfrc_output(1:inp_strlen(extfrc_output)),
     F        status='old',iostat=io_status,
     .        access='stream',action='read')
         if(io_status.gt.0) call
     s        errquit('qmd_ext: error reading output ',io_status,
     d        DISK_ERR)
         do
            read(myunit, iostat=io_status) mychar
            if(is_iostat_end(io_status)) exit
            if(io_status.gt.0) call
     s           errquit('qmd_ext: error copying output ',io_status,
     d        DISK_ERR)
            write(luout, '(a)', advance='no') mychar
         end do
      endif
#endif      
      return
      end

      
