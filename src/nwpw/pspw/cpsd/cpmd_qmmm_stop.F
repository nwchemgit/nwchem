
      logical function cpmd_qmmm_stop(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
ccc#include "frac_occ.fh"
#include "errquit.fh"
      
      logical value

      real*8 kb
      parameter (kb=3.16679d-6)

      
*     **** parallel variables ****
      integer  taskid,np
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,n2ft3d

*     **** electronic variables ****
      real*8 icharge
      integer n1(2),n2(2),neall
      real*8  en(2)
      real*8 dipole(3)


*     ***** energy variables ****
      real*8  E(30),cv

      integer eig(2)

*     **** cpmd_qmmm_block ****
      integer nfft3d
      integer npack1
      integer ispin
      integer ne(2),nemax
      integer psi0(2),psi1(2),psi2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
      integer hml(2),lmd(2)
      integer icount,icount_shift
      real*8  emotion_time_shift,dt,eave,evar,eke
      logical mulliken
      common /cpmd_qmmm_block/ psi0,psi1,psi2,dn,hml,lmd,Hpsi,psir,
     >                         ispin,ne,npack1,nfft3d,nemax,
     >                         icount,icount_shift,
     >                         emotion_time_shift,dt,eave,evar,eke,
     >                         mulliken




*     **** error variables ****
      integer ierr

*     **** local variables ****
      logical verlet,SA,found
      integer ms
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz
      real*8  ekg,eki0,eki1,sum
      real*8  eke0,eke1
      real*8  EV,pi
      integer i,j,k,ia,n,nn
      integer ii,jj,index,indx
      integer it_in,it_out
      real*8 w,sumall
      real*8 Te_init,Tr_init,Te_new,Tr_new,sa_decay(2),sa_alpha(2)
      integer nfft3,mapping
      parameter (nfft3=32)
      character*50 filename
      character*255 full_filename
  


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      character*4 ion_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    ion_atom,ion_aname


      real*8   control_rti,control_rte,ion_rion
      real*8   ion_vion,ion_com_ke,ion_ke
      real*8   ion_Temperature,ion_com_Temperature
      real*8   pspw_qmmm_Temperature
      external control_rti,control_rte,ion_rion
      external ion_vion,ion_com_ke,ion_ke
      external ion_Temperature,ion_com_Temperature
      external pspw_qmmm_Temperature
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  pspw_qmmm_nion
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pspw_qmmm_nion

      character*12 control_boundry
      external     control_boundry

      logical      psp_semicore
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  control_Nose,control_Mulliken,pspw_qmmm_found
      external control_Nose,control_Mulliken,pspw_qmmm_found

      real*8   control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0
      external control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0

      real*8   control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      external control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      logical      v_psi_filefind
      external     v_psi_filefind
      real*8   nwpw_timing
      external nwpw_timing

      logical  control_out_of_time,control_new_vpsi
      external control_out_of_time,control_new_vpsi

      logical  control_SA,control_Fei,pspw_SIC,pspw_HFX
      real*8   control_SA_decay
      external control_SA,control_Fei,pspw_SIC,pspw_HFX
      external control_SA_decay
      integer  control_mapping
      external control_mapping

      logical  control_translation,control_rotation
      external control_translation,control_rotation
     





*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
        close(unit=32)
        close(unit=34)
      end if

*     *** close fei io ****
      call fei_end()

*     **** close OMOTION file ****
      if (mulliken) call Orb_End()

*     **** ecce print ****
      call ecce_print_module_exit('driver', 'ok')
      call ecce_print_module_exit('task optimize', 'ok')


      if (taskid.eq.MASTER) CALL nwpw_message(3)
      if (taskid.eq.MASTER) call current_second(cpu3)


*         |****************************************|
*********** produce CHECK file and diagonalize hml *****************
*         |****************************************|

*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name('CHECK',.true.,
     >                               .false.,
     >                       full_filename)
         open(unit=17,file=full_filename,form='formatted')
      end if

*     **** check total number of electrons ****
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),sumall)
         en(ms) = sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do
      if (taskid.eq.MASTER) then
         write(17,1320) (en(ms),ms=1,ispin)
      end if

*     **** comparison between hamiltonian an lambda matrix ****
      n = ne(1)
      nn = n*n
      if (taskid.eq.MASTER) then
         write(17,1330)
         do ms=1,ispin
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               do j=n1(ms),n2(ms)
                  jj = j-n1(ms)
                  index = (ii+1) + jj*n +(ms-1)*nn
                  write(17,1340) ms,ii+1,jj+1,
     >                           dbl_mb(hml(1)+index-1),
     >                           dbl_mb(lmd(1)+index-1),
     >             dbl_mb(hml(1)+index-1)-dbl_mb(lmd(1)+index-1)
               end do
            end do
         end do
      end if

*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
      end if

      do ms=1,ispin
         do i=n1(ms),n2(ms)
            ii = i-n1(ms)+1
            do j=i,n2(ms)
               jj = j-n1(ms)+1
               call Pack_cc_dot(1,dcpl_mb(psi1(1)+(i-1)*npack1),
     >                          dcpl_mb(psi1(1)+(j-1)*npack1),
     >                          w)
               if (taskid.eq.MASTER) then
                  write(17,1360) ms,ii,jj,w
               end if
            end do
         end do
      end do

*     ***** do not diagonalize the hamiltonian matrix *****
      if (pspw_SIC()) then
        call dcopy(2*npack1*nemax,
     >             dcpl_mb(psi1(1)),1,
     >             dcpl_mb(psi2(1)),1)
      

*     ***** diagonalize the hamiltonian matrix *****
      else
      n = ne(1)
      nn = n*n
      call dcopy(2*nemax,0.0d0,0,dbl_mb(eig(1)),1)
      do ms=1,ispin
c        call eigen(n,ne(ms),
c    >              dbl_mb(hml(1)+(ms-1)*nn),
c    >              dbl_mb(eig(1)+(ms-1)*n),
c    >              dbl_mb(lmd(1)),ierr)

c        if (fractional) then
c         do jj=1,ne(ms)
c         do ii=1,ne(ms)
c           dbl_mb(hml(1)+(ii-1)+(jj-1)*n+(ms-1)*nn)
c     >     =dbl_mb(hml(1)+(ii-1)+(jj-1)*n+(ms-1)*nn)
c     >     *dbl_mb(fweight(1)+(jj-1)+(ms-1)*n)
c         end do
c         end do
c        end if

        call DSYEV('V','U',ne(ms),
     >              dbl_mb(hml(1)+(ms-1)*nn),n,
     >              dbl_mb(eig(1)+(ms-1)*n),
     >              dbl_mb(lmd(1)),(2*nemax*nemax),
     >              ierr)
        call eigsrt(dbl_mb(eig(1)+(ms-1)*n),
     >              dbl_mb(hml(1)+(ms-1)*nn),
     >              ne(ms),n)
c        if (fractional) then
c        do ii=1,ne(ms)
c           dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >       =dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >       /dbl_mb(fweight(1)+(ii-1)+(ms-1)*n)
c        end do
c        end if

      end do

*     *** rotate current psi ***
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi2(1)),1)
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            jj = j-n1(ms)
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               index = (ii+1) + jj*n + (ms-1)*nn
             
               call Pack_cc_daxpy(1,dbl_mb(hml(1)+index-1),
     >                            dcpl_mb(psi1(1)+(i-1)*npack1),
     >                            dcpl_mb(psi2(1)+(j-1)*npack1))
            end do
         end do
      end do

*     *** rotate current v_psi ***
      call dcopy(2*npack1*nemax,dcpl_mb(psi0(1)),1,
     >                          dcpl_mb(psi1(1)),1)
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi0(1)),1)
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            jj = j-n1(ms)
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               index = (ii+1) + jj*n + (ms-1)*nn
             
               call Pack_cc_daxpy(1,dbl_mb(hml(1)+index-1),
     >                            dcpl_mb(psi1(1)+(i-1)*npack1),
     >                            dcpl_mb(psi0(1)+(j-1)*npack1))
            end do
         end do
      end do
      end if




*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if





*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (taskid.eq.MASTER) then
         call print_elapsed_time(icount*it_in*dt)
         write(6,1300)
         write(6,1410)
         write(6,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(6,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call pspw_qmmm_Print(6)

         write(6,1421)
         write(6,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) vcx,vcy,vcz
         write(6,1210) vgx,vgy,vgz

         call pspw_qmmm_Print2(6)

         write(6,*)
         write(6,1320) en(1),en(ispin),' (real space)'
         write(6,1430) E(2),E(2)/ion_nion()
         write(6,1440) E(5),E(5)/n2(ispin)
         write(6,1450) E(6),E(6)/n2(ispin)
         write(6,1460) E(7),E(7)/n2(ispin)
         if (pspw_SIC()) then
           write(6,1455) E(16),E(16)/n2(ispin)
           write(6,1456) E(17),E(17)/n2(ispin)
         end if
         if (pspw_HFX()) then
           write(6,1457) E(20),E(20)/n2(ispin)
         end if
         write(6,1470) E(8),E(8)/ion_nion()
         write(6,1471) E(3),E(3)/n2(ispin)
         write(6,1472) ion_ke(),ion_ke()/ion_nion()


         if (pspw_qmmm_found()) then
            write(6,1700)
            write(6,1701)
            write(6,1702) E(11)
            write(6,1703) E(12)
            write(6,1704) E(13)
            write(6,1705) E(14)
            write(6,1706) E(15)
         end if


         if (control_Nose()) then
           write(6,1473) E(9),E(9)/n2(ispin)
           write(6,1474) E(10),E(10)/ion_nion()
         end if
         write(6,1226) E(3),ion_ke(),ion_com_ke()
         eke = eke/dble(it_out)
         eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))/pack_nwave_all(1)
         !eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))

*       **** write out Temperatures ****
         write(6,1491) eke
         eki0 = ion_Temperature()
         if (pspw_qmmm_found()) then
            eki1 =pspw_qmmm_Temperature()
            sum = ion_nion() + pspw_qmmm_nion() - 2.0d0
            eki0 = eki0*((ion_nion()-2.0d0)/sum)
     >           + eki1*((pspw_qmmm_nion()-2.0d0)/sum)
         end if
         write(6,1480) eki0
         write(6,1490) ion_com_Temperature()

         eave = E(25)/dble(icount+icount_shift)
         evar = E(26)/dble(icount+icount_shift)
         evar = evar - eave*eave
         cv = (evar)/(kb*ion_Temperature()**2)
         cv = cv/dble(ion_nion())
         write(6,1492) eave
         write(6,1493) evar
         write(6,1494) cv

*        **** write out diagonal <psi|H|psi> matrix ****
         if (pspw_SIC()) then
           n = ne(1)
           nn = n*n
           write(6,1331)
           do ms=1,ispin
           do i=n1(ms),n2(ms)
             ii = i-n1(ms)
             indx = (ii+1) + ii*n +(ms-1)*nn
             write(6,1341) ms,ii+1,ii+1,
     >                       dbl_mb(hml(1)+indx-1)
           end do
           end do

*        **** write out KS eigenvalues ****
         else
         write(6,1500)
         NN=NE(1)-NE(2)
         EV=27.2116d0

c         if (fractional) then
c           do i=1,NN
c             write(6,1511) dbl_mb(EIG(1)+i-1),
c     >                     dbl_mb(EIG(1)+i-1)*EV,
c     >                     dbl_mb(fweight(1)+i-1)
c           end do
c           do i=1,ne(2)
c             write(6,1511) dbl_mb(EIG(1)+i-1+NN),
c     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
c     >                     dbl_mb(fweight(1)+i-1+NN),
c     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
c     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV,
c     >                     dbl_mb(fweight(1)+i-1+n1(2)-1)
c           end do
c         else
           do i=1,NN
             write(6,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
           end do
           do i=1,ne(2)
             write(6,1510) dbl_mb(EIG(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
           end do
c         end if

         end if

*        *** Extra energy output added for QA test ****
         write(6,1600) E(2)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

*     **** calculate spin contamination ****
      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
     >                         dipole)

*     **** calculate the Dipole ***
      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)

      
*     ***** write wavefunctions and v_wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)))
      call v_psi_write(ispin,ne,dcpl_mb(psi0(1)))

*     **** write geometry to rtdb ****
      call pspw_qmmm_write(rtdb)
      call ion_write(rtdb)

*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (controL_version().eq.3) call coulomb_end()
      if (controL_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_qmmm_end()
      call dplot_iteration_end()
ccc      call frac_occ_end()
      if (control_Nose()) call Nose_end()

      value = BA_free_heap(psir(2))
      value = BA_free_heap(Hpsi(2))
      value = BA_free_heap(dn(2))
      value = BA_free_heap(eig(2))
      value = BA_free_heap(hml(2))
      value = BA_free_heap(lmd(2))
      value = BA_free_heap(psi0(2))
      value = BA_free_heap(psi1(2))
      value = BA_free_heap(psi2(2))
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)

*     **** do anaylysis on MOTION files ****
      call cpmd_properties(rtdb)


*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (taskid.eq.MASTER) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(icount*it_in)
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'prologue    : ',T1
         write(6,*) 'main loop   : ',T2
         write(6,*) 'epilogue    : ',T3
         write(6,*) 'total       : ',T4
         write(6,*) 'cputime/step: ',AV
         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Time spent doing:'
         write(6,*) '  FFTs                       : ', nwpw_timing(1)
         write(6,*) '  dot products               : ', nwpw_timing(2)
         write(6,*) '  orthonormalization         : ', nwpw_timing(3)
         write(6,*) '  exchange correlation       : ', nwpw_timing(4)
         write(6,*) '  local pseudopotentials     : ', nwpw_timing(5)
         write(6,*) '  non-local pseudopotentials : ', nwpw_timing(6)
         write(6,*) '  hartree potentials         : ', nwpw_timing(7)
         write(6,*) '  structure factors          : ', nwpw_timing(8)
         write(6,*) '  masking and packing        : ', nwpw_timing(9)
         write(6,*) '  total energy evalulation   : ', nwpw_timing(10)
         write(6,*) '  density                    : ', nwpw_timing(11)
         write(6,*) '  allocate and deallocate    : ', nwpw_timing(12)
         write(6,*)
         CALL nwpw_message(4)
      end if 


      call Parallel_Finalize()
      cpmd_qmmm_stop = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*      [   extended Lagrangian molecular   ]       *')
 1031 FORMAT(10X,'*      [        dynamics simulation        ]       *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #5.00   06/01/99              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1112 FORMAT( ' parallel mapping         : slab')
 1113 FORMAT( ' parallel mapping         : hilbert')
 1115 FORMAT(/' options:')
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i2)
 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge=',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A4  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A4  ,' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
 1192 FORMAT(5X, I4, A4  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1219 FORMAT(/' number of electrons: spin up=',F6.2,'  down=',F6.2,A)
 1220 FORMAT(/' number of electrons: spin up=',I6,  '  down=',I6,A)
 1221 FORMAT( ' number of orbitals : spin up=',I6,  '  down=',I6,A)
 1222 format(5x,' initial kinetic energy: ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)',/50x,
     >                                      e12.5,' (c.o.m.)')
 1223 format(5x,' after scaling:          ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)')
 1224 format(5x,' increased energy:       ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)')
 1226 format(/' final kinetic energy:  ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)',/44x,
     >                                      e12.5,' (c.o.m.)')
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
  
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1290 FORMAT(5X, ' cooling/heatting rates: ',e12.5,' (psi)',2x,
     >                                       e12.5,' (ion)')
 
 1294 format(/' Constant Energy Simulation                     ')
 1295 format(/' Nose-Hoover Simulation - thermostat parameters:')
 1296 format(5x, 'SA decay rates  =',e10.3,' (elc)',e10.3,' (ion)')
 1297 format(5x, 
     > 'Te =',f8.2,' Qe =',e10.3,' 2*pi/we=',e10.3,' Ee0=',e10.3)
 1298 format(5x, 
     > 'Tr =',f8.2,' Qr =',e10.3,' 2*pi/wr=',e10.3,' Er0=',e10.3)
 1299 format(//' Mulliken analysis output generated            ')
 1300 FORMAT(//)
 1305 FORMAT(10X,'============ Car-Parrinello iteration ==============')
 1309 FORMAT(I8,2E19.10,2E14.5,2F9.1)
 1310 FORMAT(I8,2E19.10,2E14.5,F14.2)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Diagonal Elements of Hamiltonian matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I3,2I3,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT(/' Kinetic energy (elc)    :',E19.10,' (',E15.5,'/elc)')
 1472 FORMAT( ' Kinetic energy (ion)    :',E19.10,' (',E15.5,'/ion)')
 1473 FORMAT( ' thermostat energy (elc) :',E19.10,' (',E15.5,'/elc)')
 1474 FORMAT( ' thermostat energy (ion) :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(' Temperature :    ',F10.1,' K (ion)')
 1490 FORMAT('             :    ',F10.1,' K (c.o.m.)')
 1491 FORMAT(' Temperature :    ',F10.1,' K (elc)')
 1492 FORMAT(/' Eaverage           :    ',E19.10)
 1493 FORMAT( ' Evariance          :    ',E19.10)
 1494 FORMAT( ' Cv - f*kb/(2*nion) :    ',E19.10)
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' QM/MM energy           :',E19.10)
 1703 FORMAT( ' MM/MM energy           :',E19.10)
 1704 FORMAT( ' MM Polarization energy :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      cpmd_qmmm_stop = value
      return
      END



c $Id$
