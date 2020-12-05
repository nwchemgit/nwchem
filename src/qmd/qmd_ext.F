      subroutine qmd_ext_init(nat,dt_nucl)
      implicit none
#include "mafdecls.fh"
#include "util_params.fh"
      integer nat     ![in]
      integer dt_nucl ![in]
c     pass plumed initialization bits
#ifdef USE_PLUMED
      integer ga_comm
      integer*4 ga_comm4
      double precision kbt
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
#if 0      
c     kBT
c      kbt=1000d0
c     call plumed_f_gcmd("setKbT"//char(0),kbt)
#endif      
      call plumed_f_gcmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
      call plumed_f_gcmd("setLogFile"//char(0),"plumed.out"//char(0))
      call plumed_f_gcmd("setNatoms"//char(0),nat)
      call plumed_f_gcmd("setMDEngine"//char(0),"nwchem")
      call ga_mpi_comm_pgroup_default(ga_comm)
      call plumed_f_gcmd("setMPIFComm"//char(0),ga_comm4)
      call plumed_f_gcmd("setTimestep"//char(0),dt_nucl)
      call plumed_f_gcmd("init"//char(0),0)
#endif      
      return
      end

      subroutine qmd_ext_gradient(nat,istep_nucl,energy,m,g,r)
      implicit none
#include "global.fh"      
      integer nat                ! [in]
      integer istep_nucl         ! [in]
      double precision energy    ! [in]
      double precision g(3,nat)  ! gradients [in/out]
      double precision r(3,nat)  ! coordinates [in]
      double precision m(nat)    ! masses [in]
#ifdef USE_PLUMED
      integer i
c     call plumed: r[in], g[out]
c     pointer to current timestep
      call plumed_f_gcmd("setStep"//char(0),istep_nucl)
      call plumed_f_gcmd("setMasses"//char(0),m)
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
       call plumed_f_gcmd("setEnergy"//char(0),energy)
c     forces before plumed call
      call plumed_f_gcmd("calc"//char(0),0)
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
      subroutine qmd_ext_final()
      implicit none
#ifdef USE_PLUMED      
c     call plumed finalization bits
      call plumed_f_gfinalize()
#endif      
      return
      end

      
