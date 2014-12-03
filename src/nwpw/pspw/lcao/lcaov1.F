*
* $Id$
*

      logical function lcaov1(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "btdb.fh"
#include "errquit.fh"

      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d

*     **** electronic variables ****
      real*8 icharge,tcharge
      integer ispin
      integer ne(2)


*     **** error variables ****
      integer ierr

*     **** local variables ****
      integer ms,mapping,mapping1d
      real*8  deltae,deltac,deltar
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2
      integer i,j,k,ia,n,nn,nextra_orbs
      integer ii,jj,indx,indx1
      integer icount,it_in,it_out
      character*255 full_filename

      logical value,lcao_print,success
      integer version
      real*8 unita(3,3),t
      integer ind,lcao_iterations,x,xx,mult
      character*50 filename


*     **** external functions ****
      character*50 control_input_psi
      external     control_input_psi
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg

      logical  control_read,ion_init
      external control_read,ion_init
 
      logical      psp_semicore
      external     psp_semicore

      integer      control_np_orbital,control_mapping,control_mapping1d
      external     control_np_orbital,control_mapping,control_mapping1d

      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_ispin
      external control_ngrid,control_ispin
      integer  control_multiplicity
      external control_multiplicity
      real*8   control_totalcharge
      external control_totalcharge
      real*8   ion_TotalCharge_qm
      external ion_TotalCharge_qm
      logical  lcao_wvfnc,control_print,psi_filefind
      external lcao_wvfnc,control_print,psi_filefind
      integer  control_rtdb,control_version
      external control_rtdb,control_version
      integer  control_frac_occ_extra_orbitals
      external control_frac_occ_extra_orbitals


*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.

      call nwpw_timing_init()

*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      

      if (taskid.eq.MASTER) call current_second(cpu1)

      value = control_read(1,rtdb)
      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      mapping = control_mapping()

*     **** initialize psi_data ****
      call psi_data_init(20)

*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_Init()

      call D3dB_pfft_init()
      call ga_sync()

*     **** read ions ****
      value = ion_init(rtdb)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      call phafac()
      if (control_version().eq.3) then
         call ewald_init()
         call ewald_phafac()
      end if

*     **** call the lcao code ****
      filename = control_input_psi()
      version = control_version()
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      unita(1,1) = control_unita(1,1)
      unita(2,1) = control_unita(2,1)
      unita(3,1) = control_unita(3,1)
      unita(1,2) = control_unita(1,2)
      unita(2,2) = control_unita(2,2)
      unita(3,2) = control_unita(3,2)
      unita(1,3) = control_unita(1,3)
      unita(2,3) = control_unita(2,3)
      unita(3,3) = control_unita(3,3)

      ispin   = control_ispin()
      mult    = control_multiplicity()
      tcharge = control_TotalCharge()
      icharge = ion_TotalCharge_qm()
      t = icharge - tcharge       !** total number of electrons **
      x = NINT(t)
      xx = x

*     **** reassign spin to agree with total number of electrons ****
      if ((mod(x,2).ne.0).and.(ispin.eq.1)) then !** odd number of electrons **
         ispin = 2
         call control_ispin_set(ispin)
      end if

*     **** reassign mult to agree with total number of electrons ****
      if ((mod(x,2).ne.0) .and.(mod(mult,2).ne.0)) then !** odd number of electrons and mult odd **
         mult = mult - 1
         do while (mult.gt.(x+1))
            mult = mult - 2
         end do
         if (mult.lt.1) mult = 2
         call control_multiplicity_set(mult)
      end if
      if ((mod(x,2).eq.0) .and.(mod(mult,2).eq.0)) then !** even number of electrons and mult even **
         mult = mult - 1
         do while (mult.gt.(x+1))
            mult = mult - 2
         end do
         if (mult.lt.1) mult = 1
         call control_multiplicity_set(mult)
      end if

*     **** assign number of up and down electrons ****
      if (ispin.eq.1) then
         ne(1) = xx/2
         ne(2) = 0
      else
         x = mult-1
         ne(1) = (xx+x)/2
         ne(2) = (xx-x)/2
      end if

      nextra_orbs = control_frac_occ_extra_orbitals()
      if (ispin.eq.1) then
         ne(1) = ne(1) + nextra_orbs
      else
         ne(1) = ne(1) + nextra_orbs
         ne(2) = ne(2) + nextra_orbs
      end if

      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,
     >         'nwpw:lcao_iterations',mt_int,1,lcao_iterations))
     >   lcao_iterations = 10

      lcao_print = .true.


      success =  lcao_wvfnc(lcao_print,lcao_iterations,
     >                        filename,
     >                        version,
     >                        ngrid,
     >                        unita,
     >                        ispin,
     >                        ne)



*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
      call psp_end()


      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call psi_data_end()


      call Parallel2d_Finalize()
      call Parallel_Finalize()
      lcaov1 = value
      return
      end
