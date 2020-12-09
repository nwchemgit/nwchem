      subroutine qmd_ext_init(nat,dt_nucl,kbt,ext_forces)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util_params.fh"
      integer nat             ![in]
      integer dt_nucl         ![in]
      double precision kbt    ![in]
      character*(*) ext_forces
c     pass plumed initialization bits
#ifdef USE_PLUMED
      integer ga_comm
      integer*4 ga_comm4
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
      write(6,*) ' kbt ',kbt
      call plumed_f_gcmd("setKbT"//char(0),kbt)
      call plumed_f_gcmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
      call plumed_f_gcmd("setLogFile"//char(0),"plumed.out"//char(0))
      call plumed_f_gcmd("setNatoms"//char(0),nat)
      call plumed_f_gcmd("setMDEngine"//char(0),"nwchem")
#ifdef PLUMED_HASMPI 
      call plumed_f_gcmd("setMPIFComm"//char(0),ga_comm4)
#endif
      call plumed_f_gcmd("setTimestep"//char(0),dt_nucl)
      write(6,*) ' plumed init '
      call plumed_f_gcmd("init"//char(0),0)
      endif
#endif      
      return
      end

      subroutine qmd_ext_gradient(nat,istep_nucl,energy,m,g,r,
     E     ext_forces)
      implicit none
#include "global.fh"      
#include "mafdecls.fh"
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
      write(6,*) ' qmd_ext_gr istep ',istep_nucl
      call plumed_f_gcmd("setStep"//char(0),istep_nucl)
      call plumed_f_gcmd("setMasses"//char(0),m)
c plumed uses forces not gradientes. double flip
c first flip
        call dscal(3*nat, -1d0, g, 1)
       if(ga_nodeid().eq.0) then
          write(6,*) '### before plumed ##'
          do i=1,nat
             write(6,*) 'frc at ',i,g(1,i),g(2,i),g(3,i)
          enddo
       endif
      call plumed_f_gcmd("setForces"//char(0),g)
      call plumed_f_gcmd("setPositions"//char(0),r)

c      call plumed_f_gcmd("setBox"//char(0),at_plumed)
c      call plumed_f_gcmd("setVirial"//char(0),virial)
      write(6,*) ' E ext', energy
       call plumed_f_gcmd("setEnergy"//char(0),energy)
c     forces before plumed call
       call plumed_f_gcmd("calc"//char(0),0)
      endif
c second sign flip
      call dscal(3*nat, -1d0, g, 1)
c     broadcast from 0
      call ga_brdcst(2020,g,MA_sizeof(MT_DBL,1,MT_BYTE),0)
       if(ga_nodeid().eq.0) then
          write(6,*) '### after plumed ##'
          do i=1,nat
             write(6,*) 'frc at ',i,g(1,i),g(2,i),g(3,i)
          enddo
       endif
c     check sign agreement with plumed ...
#endif

      return
      end
      subroutine qmd_ext_final(ext_forces)
      implicit none
#include "global.fh"
#include "inp.fh"
      character*(*) ext_forces
#ifdef USE_PLUMED
      if(ext_forces(1:inp_strlen(ext_forces)).ne.'plumed') call
     E     errquit(' qmd_ext works only with plumed ',0,0)
      if(ga_nodeid().eq.0) then
c     call plumed finalization bits
         call plumed_f_gfinalize()
      endif
#endif      
      return
      end

      
