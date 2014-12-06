c     $Id$                       

*     ***************************
*     *				*
*     *	    c_cgsd_excited      *
*     *				*
*     ***************************

      subroutine c_cgsd_excited()
      implicit none

c#include "global.fh"
c#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"

      logical newpsi
      integer MASTER
      parameter (MASTER=0)

      logical stalled,value
      integer taskid
      integer minimizer
      integer i,j,ms,neall,NN,nbrillioun,nb
      integer nexcited(2),ispin
      integer it_in,it_out,icount,bfgscount,vers
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz,f1,f2,f3,f4,f5,f6
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(30)
      integer rtdb

*     **** external functions ****
      logical  control_DOS,control_Mulliken
      logical  ecpsi_initialize,ecpsi_finalize,epsi_filefind
      integer  cpsi_ne_excited,control_excited_ne,control_rtdb
      integer  epsi_get_version
      real*8   cpsi_virtual_brdcst,brillioun_ks_brdcst
      real*8   brillioun_k_brdcst
      external control_DOS,control_Mulliken
      external ecpsi_initialize,ecpsi_finalize,epsi_filefind
      external cpsi_ne_excited,control_excited_ne,control_rtdb
      external epsi_get_version
      external cpsi_virtual_brdcst,brillioun_ks_brdcst
      external brillioun_k_brdcst
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun
      logical  pspw_reformat_c_ewvfnc
      external pspw_reformat_c_ewvfnc
      

      call Parallel_taskid(taskid)
      neall = control_excited_ne(1)+control_excited_ne(2)
      nbrillioun = brillioun_nbrillioun()

      if (neall.gt.0) then

      if (taskid.eq.MASTER) write(luout,1301)

*     **** read in excited wavefunctions and initialize epsi ****
      if (.not.epsi_filefind()) then
        call ecpsi_new()
        newpsi = .true.
      else
         newpsi = .false.

*        **** convert from pspw format to band format ****
         vers = epsi_get_version()
         if ((vers.eq.3).or.(vers.eq.4)) then
           newpsi = .true.
           value = btdb_parallel(.false.)
           if (taskid.eq.MASTER) then
             value= pspw_reformat_c_ewvfnc(1)
           end if 
           value = btdb_parallel(.true.)
         end if
      end if

      call psi_get_ne_excited(ispin,nexcited)
c      call Cnex_init(ispin,nexcited,control_mapping1d())

      value = ecpsi_initialize()
      call cpsi_gen_density_potentials(1)
      call cpsi_minimize_virtual()

      NN=cpsi_ne_excited(1)-cpsi_ne_excited(2)
      EV=27.2116d0

      do nb=1,nbrillioun
        f1 = brillioun_ks_brdcst(1,nb)
        f2 = brillioun_ks_brdcst(2,nb)
        f3 = brillioun_ks_brdcst(3,nb)
        f4 = brillioun_k_brdcst(1,nb)
        f5 = brillioun_k_brdcst(2,nb)
        f6 = brillioun_k_brdcst(3,nb)
        if (taskid.eq.MASTER) then
           write(luout,1508) nb,f1,f2,f3,f4,f5,f6
           write(luout,1500)
        end if
        do i=cpsi_ne_excited(1),(cpsi_ne_excited(2)+1),-1
           f1 = cpsi_virtual_brdcst(nb,1,i)
           if (taskid.eq.MASTER) write(luout,1510) f1,f1*EV
        end do
        do i=cpsi_ne_excited(2),1,-1
           f1 = cpsi_virtual_brdcst(nb,1,i)
           f2 = cpsi_virtual_brdcst(nb,2,i)
           if (taskid.eq.MASTER)  write(luout,1510) f1,f1*EV,f2,f2*EV
        end do
        if (taskid.eq.MASTER) then
           write(luout,*)
           write(luout,*)
        end if
      end do


*     **** calculate Density of States ***
      rtdb = control_rtdb()
      if (control_DOS()) call ecpsi_DOS(rtdb)

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) call ecpsi_Mulliken(rtdb)


*     **** write out excited wavefunctions ***
      value = ecpsi_finalize(.true.)

c      call Cnex_end()
      end if

      return
 1301 FORMAT(//'== Virtual Orbital Calculation =='/)
 1500 FORMAT(/' virtual orbital energies:')
 1508 FORMAT(/' Brillouin zone point: ',i6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')

 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1512 FORMAT(2(E18.7,' (',F8.3,'eV)',A4))
      end


