*
* $Id: pspw_translate.F 26433 2014-12-03 23:34:20Z mjacquelin $
*

c     **********************************************
c     *                                            *
c     *                pspw_translate              *
c     *                                            *
c     **********************************************

c Translates the ions and wavefunctions

      logical function pspw_translate(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "errquit.fh"


*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)


      logical value,oprint,reorder
      integer mapping,mapping1d
      integer i,j,k,ms,n,ne(2),ispin,ngrid(3),nemax
      integer neq(2),nemaxq,shift
      integer nfft3d,n2ft3d,npack1
      integer np,nion

      integer psi(2),occ(2)
      real*8  rshift(3)
      real*8  cx,cy,cz,gx,gy,gz

      character*50 geom_name

*     **** psi smearing block ****
      logical fractional
      integer smearoccupation

*     **** external functions ****
      integer ion_katm,ion_nion,ion_nkatm,ion_natm
      real*8  ion_rion,ion_amass
      logical      control_read
      integer      control_ngrid,control_ispin,control_mapping
      integer      control_mapping1d
      integer      pack_nwave,pack_nwave_all
      integer      control_np_orbital
      real*8       lattice_unita,lattice_unitg,control_wannier_timestep
      real*8       lattice_ecut,lattice_wcut,lattice_omega
      external     control_read
      external     control_ngrid,control_ispin,control_mapping
      external     control_mapping1d
      external     pack_nwave,pack_nwave_all
      external     control_np_orbital
      external     lattice_unita,lattice_unitg,control_wannier_timestep
      external     lattice_ecut,lattice_wcut,lattice_omega
      external     ion_katm,ion_nion,ion_nkatm,ion_natm
      external     ion_rion,ion_amass
      character*4 ion_aname,ion_atom
      external    ion_aname,ion_atom
      character*9 ion_amm
      external    ion_amm

      character*50 control_input_psi,control_output_psi
      external     control_input_psi,control_output_psi

      value = .true.

*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      oprint = (taskid.eq.MASTER)

*     ***** print out header ****
      if (oprint) then
         write(luout,1000)
         write(luout,1010)
         write(luout,1020)
         write(luout,1010)
         write(luout,1030)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1050)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
      end if
      
      value = control_read(15,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      mapping = control_mapping()

      call Parallel2d_Init(control_np_orbital())

*     **** initialize psi_data ****
      !call psi_data_init(100)

*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)      
      call D3dB_pfft_init()
      call ga_sync()

      call ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call strfac_init()
      !call psp_init()
      !call psp_readall()


*     ***** allocate psi wavefunctions ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
         fractional = .true.
      else
         fractional = .false.
      end if
      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)
      nemaxq = neq(1)+neq(2)
      value = BA_alloc_get(mt_dcpl,npack1*nemaxq,
     >                     'psi',psi(2),psi(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ',occ(2),occ(1))
      if (.not. value) 
     >   call errquit('pspw_translate:out of heap memory',0,0)

*     *****  read psi wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi(1)),
     >              smearoccupation,dbl_mb(occ(1)))

*     ****  get the translate vector and output geom_name ****
      call control_translate_vector(rshift)
      call control_translate_geom_name(geom_name)
      call control_translate_reorder(reorder)



      if (oprint) then
         write(luout,109) control_input_psi()
         write(luout,110) control_output_psi()
         write(luout,115) geom_name
         if (reorder) then
            write(luout,120) rshift
         else
            write(luout,121) rshift
         end if


         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz

         write(luout,1220) ne(1),ne(ispin),' ( Fourier space)'
         write(luout,1230)
         write(luout,1241) lattice_unita(1,1),
     >                 lattice_unita(2,1),
     >                 lattice_unita(3,1)
         write(luout,1242) lattice_unita(1,2),
     >                 lattice_unita(2,2),
     >                 lattice_unita(3,2)
         write(luout,1243) lattice_unita(1,3),
     >                 lattice_unita(2,3),
     >                 lattice_unita(3,3)
         write(luout,1244) lattice_unitg(1,1),
     >                 lattice_unitg(2,1),
     >                 lattice_unitg(3,1)
         write(luout,1245) lattice_unitg(1,2),
     >                 lattice_unitg(2,2),
     >                 lattice_unitg(3,2)
         write(luout,1246) lattice_unitg(1,3),
     >                 lattice_unitg(2,3),
     >                 lattice_unitg(3,3)
         write(luout,1231) lattice_omega()
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(luout,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)

      end if

*     ****  translate the ions and wavefunctions ****
      if (reorder) then
         call ion_translate_reorder(rshift)
      else
         call ion_translate(rshift)
      end if
      call psi_translate(rshift,
     >                   npack1,(ne(1)+ne(2)),
     >                   dcpl_mb(psi(1)))



      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

      if (oprint) then
         write(luout,1300)
         write(luout,1420)
         do I=1,ion_nion()
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
      end if


*     ***** write psi wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi(1)),
     >               smearoccupation,dbl_mb(occ(1)))

*     ***** write the ions positions *****
      call ion_rion_save(rtdb,geom_name)

*     **** deallocate heap memory ****
      call ion_write(rtdb)
      call ion_end()
      call strfac_end()
      !call psp_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call D3dB_pfft_end()
      call D3dB_end(1)
      call Dne_end()
      !call psi_data_end()


      value = value.and.BA_free_heap(psi(2))
      value = value.and.BA_free_heap(occ(2))
      if (.not. value) 
     >   call errquit('pspw_translate:deallocating heap memory',2,0)


 9000 if (oprint) call nwpw_message(4)
      call Parallel2d_Finalize()
      call Parallel_Finalize()
      pspw_translate = value
      return 


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
  109 format(/' input filename : ',A)
  110 format(' output filename: ',A)
  115 format(/' geometry name: ',A)
  120 format(/,' translate vector=<',3f8.3,' >',
     >         ' - reordering atom indexes')
  121 format(/,' translate vector=<',3f8.3,' >')

 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*       PSPW Tranlate Ions and Wavefunctions       *')
 1030 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #1.00    4/17/17              *')
 1050 FORMAT(10X,'*           Developed by Eric J. Bylaska           *')
 1100 FORMAT(//)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A2,':',I5))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I4,'  spin down=',I4,A)
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1250 FORMAT(5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1300 FORMAT(//)
 1420 FORMAT( ' final position of ions:')
      end
