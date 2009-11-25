*
* $Id: band_dplot.F,v 1.12 2009-02-07 03:49:57 bylaska Exp $
*
      logical function band_dplot(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"


      logical value,spin_orbit
      integer taskid,np
      integer MASTER
      parameter (MASTER=0)

*     **** local variables ***

      integer count,i,ia,nbrillioun
      integer ngrid(3)
      integer nfft3d,n2ft3d,npack1,mapping
      integer ispin,ne(2),nemax
      integer psi2(2),dn(2),psir(2)
      real*8  cpu1,cpu2
      real*8  position_tolerance,origin(3)
      character*8 tag

*     **** external functions ****
      logical      control_read,ion_init,control_spin_orbit
      integer      Cram_nwave_all,control_mapping,brillioun_nbrillioun
      integer      control_ngrid,Cram_nwave,control_version
      real*8       lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8       lattice_unitg
      character*12 control_boundry

      external     control_spin_orbit
      external     control_read,ion_init
      external     Cram_nwave_all,control_mapping,brillioun_nbrillioun
      external     control_ngrid,Cram_nwave,control_version
      external     lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external     lattice_unitg
      external     control_boundry

*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.

*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) call current_second(cpu1)

*     ***** print out header ****
      if (taskid.eq.MASTER) then
         write(6,1000)
         write(6,1010)
         write(6,1020)
         write(6,1010)
         write(6,1030)
         write(6,1010)
         write(6,1035)
         write(6,1010)
         write(6,1040)
         write(6,1010)
         write(6,1000)
         write(6,*)
         call nwpw_message(1)
         write(6,1110)
         call util_flush(6)
      end if

*     **** get position_tolerance ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0


*     **** get origin ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:origin',
     >               mt_dbl,3,origin)) then
        origin(1) = 0.0d0 
        origin(2) = 0.0d0 
        origin(3) = 0.0d0 
      end if

*     **** read control file ****
      value = control_read(10,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      mapping   = control_mapping()


*     **** initialize D3dB data structure ****
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call C3dB_nfft3d(1,nfft3d)
      call C3dB_n2ft3d(1,n2ft3d)


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_init()
      call Cram_init()
      call C3dB_pfft_init()
      call Cram_max_npack(npack1)

*     **** initialize ion data structure ****
      value = ion_init(rtdb)


*     ***** allocate psi2 wavefunctions ****
      call psi_get_ne(ispin,ne)
      nemax = ne(1)+ne(2)
      nbrillioun = brillioun_nbrillioun()
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*nemax*nbrillioun,
     >                     'psi2',psi2(2),psi2(1))
      if (.not. value) call errquit('out of heap memory',1, MA_ERR)


*     *****  read psi2 wavefunctions ****
      spin_orbit=control_spin_orbit()
      call cpsi_read(spin_orbit,ispin,ne,nbrillioun,dcpl_mb(psi2(1)))

*     **** allocate other variables *****
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*nemax*nbrillioun,
     >                     'psir',psir(2),psir(1))

      if (.not. value) call errquit('band_dplot:out of heap memory',2,
     &       MA_ERR)

      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call cphafac()
      call cphafac_k()


*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
         write(6,1111) np
         write(6,1115)
         write(6,1121) control_boundry(),control_version()
         write(6,1220) ne(1),ne(ispin),' ( fourier space)'
         write(6,1225) position_tolerance
         write(6,1226) origin

         write(6,1230)
         write(6,1241) lattice_unita(1,1),
     >                 lattice_unita(2,1),
     >                 lattice_unita(3,1)
         write(6,1242) lattice_unita(1,2),
     >                 lattice_unita(2,2),
     >                 lattice_unita(3,2)
         write(6,1243) lattice_unita(1,3),
     >                 lattice_unita(2,3),
     >                 lattice_unita(3,3)
         write(6,1244) lattice_unitg(1,1),
     >                 lattice_unitg(2,1),
     >                 lattice_unitg(3,1)
         write(6,1245) lattice_unitg(1,2),
     >                 lattice_unitg(2,2),
     >                 lattice_unitg(3,2)
         write(6,1246) lattice_unitg(1,3),
     >                 lattice_unitg(2,3),
     >                 lattice_unitg(3,3)
         write(6,1231) lattice_omega()
         ia = brillioun_nbrillioun()
         write(6,1255)
         write(6,1256) ia

         write(6,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 Cram_nwave_all(0),Cram_nwave(0)
         do i=1,ia
          write(6,1251) i,lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 Cram_nwave_all(i),Cram_nwave(i)
         end do

        write(6,*)
        write(6,*)
        call util_flush(6)
      end if



      !**** translate system if origin is not zero ****
      if ((origin(1).ne. 0.0d0).or.
     >    (origin(2).ne. 0.0d0).or .
     >    (origin(3).ne. 0.0d0) )     then

        if (taskid.eq.MASTER) then
          write(*,*) "...translating origin..."
          write(*,*) 
          write(*,*) 
        end if

        origin(1) = -origin(1) !* translate system by -origin *
        origin(2) = -origin(2)
        origin(3) = -origin(3)

        call ion_translate(origin)
        call psi_translate(origin,
     >                     npack1,(ne(1)+ne(2)),
     >                     dcpl_mb(psi2(1)))

        call cphafac()  !** recompute phase factors **

      end if



      call dplot_gen_cpsi_dn(ispin,ne,
     >                npack1,nfft3d,nemax,nbrillioun,
     >                dcpl_mb(psi2(1)),
     >                dbl_mb(dn(1)),
     >                dcpl_mb(psir(1)))



      call band_dplot_loop(rtdb,
     >                ispin,ne,
     >                npack1,nfft3d,nemax,nbrillioun,
     >                dcpl_mb(psi2(1)),
     >                dbl_mb(dn(1)),
     >                dcpl_mb(psir(1)),
     >                .false.,tag)

      value = value.and.rtdb_get(rtdb,'band_dplot:count',
     >                            mt_int,1,count)


     
*     **** deallocate heap memory ****
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call ion_write(rtdb)  !*** can also use call ion_destroy()????
      call ion_end()
      call brillioun_end()
      call Cram_end()
      call c_G_end()
      call Parallel_Finalize()



      value = MA_free_heap(psir(2))
      value = value.and.MA_free_heap(dn(2))
      value = value.and.MA_free_heap(psi2(2))
      call C3dB_pfft_end()
      call C3dB_end(1)


      if (.not. value) call errquit('band_dplot:freeing heap memory',0,
     &       MA_ERR)


      if (taskid.eq.MASTER) then
         call current_second(cpu2)
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'total       : ',(cpu2-cpu1)
         write(6,*) 
         call nwpw_message(4)
      end if

      band_dplot = value

      return 

*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*                   band DPLOT                     *')
 1030 FORMAT(10x,'*    [ Generates density and orbital grids  ]      *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #1.00   08/01/05              *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'============ BAND DPLOT input data =================')
 1111 FORMAT(/' number of processors used:',I3)
 1115 FORMAT(/' options:')
 1121 FORMAT(5X,' boundary conditions   = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
 1225 FORMAT(/' position tolerance = ',E12.6)
 1226 FORMAT(/5x,'      origin=<',3f8.3,' >')

 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice: a1=<',3f8.3,' >')
 1242 FORMAT(5x,'          a2=<',3f8.3,' >')
 1243 FORMAT(5x,'          a3=<',3f8.3,' >')
 1244 FORMAT(5x,'          b1=<',3f8.3,' >')
 1245 FORMAT(5x,'          b2=<',3f8.3,' >')
 1246 FORMAT(5x,'          b3=<',3f8.3,' >')


 1249 FORMAT(/' computational grids:')
 1250 FORMAT(5X,' density     cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc ',I3,' cutoff=',F7.3,
     &        '  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')

 1255 FORMAT(/' brillouin zone:')
 1256 FORMAT(5x,'number of zone points:',I3)
      end





      subroutine band_dplot_loop(rtdb,
     >                      ispin,ne,
     >                      npack1,nfft3d,nemax,nbrill,
     >                      psi,
     >                      dn,
     >                      psi_r,
     >                      add_tag,tag)
      implicit none
      integer    rtdb
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax,nbrill
      complex*16 psi(npack1,nemax,nbrill)
      real*8     dn(nfft3d,2)
      complex*16 psi_r(nfft3d,nemax,nbrill)
      logical add_tag
      character*8 tag

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,grid3d
      integer taskid
      integer MASTER
      parameter (MASTER=0)

      integer n2ft3d
      integer i,count,number,ia,nb,onumber
      integer n1(2),n2(2),nn(2)
      integer nx,ny,nz
      integer rho(2),orb(2)
      real*8  scal1,scal2,kv(3),w
      complex*16 phafac

      character*72 cube_comment
      character*30 name1,name2,name3
      character*30 filename
      integer name1_len,name2_len,name3_len,ind,ind2


*     **** external functions ****
      integer  control_version
      real*8   lattice_omega, brillioun_k
      external control_version
      external lattice_omega, brillioun_k

      
      grid3d = .false.
      if (rtdb_get(rtdb,'band_dplot:3d_grid:nx',mt_int,1,i))
     >  grid3d = .true.

      call Parallel_taskid(taskid)
      value =           MA_push_get(mt_dbl,nfft3d,'rho',rho(2),rho(1))
      value = value.and.MA_push_get(mt_dcpl,nfft3d,'orb',orb(2),orb(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()



*     ********************************************
*     **** loop over orbital and density list ****
*     ********************************************

      value = value.and.rtdb_get(rtdb,'band_dplot:count',
     >                            mt_int,1,count)
      do i=1,count

*       **** define name  - not very elegent and could break if ****
*       ****                      count becomes very large      ****
        ia = ICHAR('a')
        name1 = 'band_dplot:filename'//CHAR(i-1+ia)
        name2 = 'band_dplot:number'//CHAR(i-1+ia)
        name1_len = index(name1,' ') - 1
        name2_len = index(name2,' ') - 1

        value = rtdb_cget(rtdb,name1(1:name1_len),1,filename)
        ind = index(filename,' ') - 1
        if (add_tag) then
          ind2 = index(tag,' ') - 1
          filename = filename(1:ind)//tag(1:ind2)//'.cube'
          ind = index(filename,' ') - 1
        end if


        value = value.and.rtdb_get(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
        if (.not.value)
     >     call errquit(
     >     'band_dplot: rtdb_get failed for orbital', 0, RTDB_ERR)


*       **** outputing density ****
        if (number.lt.0) then
           number = -number

          goto ( 710, 720, 730, 740, 750, 760, 770, 780 ) number
          call errquit(
     >    'band_dplot_loop: unimplemented directive',number,INPUT_ERR)
       
*          *************
*          *** total ***
*          *************
 710       call C3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing total density',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF Total Density"
           goto 790

*          ******************
*          *** difference ***
*          ******************
 720       call C3dB_rr_Sub(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing difference density',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF Spin Density"
           goto 790

*          *************
*          *** alpha ***
*          *************
 730       call C3dB_r_Copy(1,dn(1,1),dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing alpha density',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF Alpha Density"
           goto 790

*          ************
*          *** beta ***
*          ************
 740       call C3dB_r_Copy(1,dn(1,ispin),dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing beta density',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF Beta Density"
           goto 790

*          *****************
*          *** laplacian ***
*          *****************
 750       write(*,*) "laplacian not implemented"
           call C3dB_r_Zero(1,dbl_mb(rho(1)))
           cube_comment = "SCF Laplacian Density"
           goto 790

*          *******************************
*          *** Electrostatic Potential ***
*          *******************************
 760       write(*,*) "Electrostatic potential not implemented"
           call C3dB_r_Zero(1,dbl_mb(rho(1)))
           cube_comment = "SCF Electrostatic Potential"
           goto 790

*          ********************
*          ***      ELF     ***
*          ********************
 770       write(*,*) "ELF not implemented"
           call C3dB_r_Zero(1,dbl_mb(rho(1)))
           cube_comment = "SCF ELF"
           goto 790

 780       write(*,*) "ELF not implemented"
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing down spin ELF',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF ELF"
           goto 790


 790       continue
           number = -number

           if (grid3d) then
             call band_dplot_gcube_write3d(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
           else
             call band_dplot_gcube_write(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
           endif


*       **** outputing wavefunction ****
        else
           name3 = 'band_dplot:nb'//CHAR(i-1+ia)
           name3_len = index(name3,' ') - 1
           if (.not.rtdb_get(rtdb,name3(1:name3_len),mt_int,1,nb)) then
             nb = 1
           end if
           name3 = 'band_dplot:onumber'//CHAR(i-1+ia)
           name3_len = index(name3,' ') - 1
           if (.not.rtdb_get(rtdb,name3(1:name3_len),
     >                       mt_int,1,onumber)) then
             onumber = 1
           end if
           if ((number.gt.nemax).or.(number.lt.1)) then
              if (taskid.eq.MASTER) then
                write(*,*)  '   Bad orbital number ', number,
     >                     ', changing to orbital number 1.'
              end if
              number = 1
           end if
           kv(1) = brillioun_k(1,nb)
           kv(2) = brillioun_k(2,nb)
           kv(3) = brillioun_k(3,nb)
           cube_comment = 'SCF Molecular Orbitals'
       
           if (onumber.eq.1) then
             call C3dB_r_Zero(1,dbl_mb(rho(1)))
             call C3dB_cr_aSqrpy(1, 1.0d0,
     >                         psi_r(1,number,nb),
     >                         dbl_mb(rho(1)))
c             call C3dB_r_SMul(1,scal2,dbl_mb(rho(1)),dbl_mb(rho(1)))
             call C3dB_r_SMul1(1,scal2,dbl_mb(rho(1)))

             if (.not.add_tag) then
             if (taskid.eq.MASTER) then
             write(*,801) number,nb,filename(1:ind)
 801         format('    writing density of orbital ',I4,
     >              ' of brillioun zone ',I4,
     >              ' to filename: ',A)
             end if
             end if
             if (grid3d) then
               call band_dplot_gcube_write3d(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
             else
               call band_dplot_gcube_write(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
             endif

           else if (onumber.eq.2) then
             call C3dB_r_Zero(1,dbl_mb(rho(1)))
             if (.not.add_tag) then
             if (taskid.eq.MASTER) then
             write(*,802) number,nb,filename(1:ind)
 802         format('    writing real part of orbital ',I4,
     >              ' of brillioun zone ',I4,
     >              ' to filename: ',A)
             end if
             end if
             if (grid3d) then
               call band_dplot_gcube_write3d_orb(rtdb,filename,
     >                         number,cube_comment,
     >                         kv,psi_r(1,number,nb),.true.)
             else
               call C3dB_c_Copy(1,psi_r(1,number,nb),dcpl_mb(orb(1)))
               call band_kvtimesorb(kv,dcpl_mb(orb(1)))
               call C3dB_cr_real(1,
     >                           dcpl_mb(orb(1)),
     >                           dbl_mb(rho(1)))
               call band_dplot_gcube_write(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
             endif

           else if (onumber.eq.3) then
             call C3dB_r_Zero(1,dbl_mb(rho(1)))
             if (.not.add_tag) then
             if (taskid.eq.MASTER) then
             write(*,803) number,nb,filename(1:ind)
 803         format('    writing imaginary part of orbital ',I4,
     >              ' of brillioun zone ',I4,
     >              ' to filename: ',A)
             end if
             end if
             if (grid3d) then
               call band_dplot_gcube_write3d_orb(rtdb,filename,
     >                         number,cube_comment,
     >                         kv,psi_r(1,number,nb),.false.)
             else
               call C3dB_c_Copy(1,psi_r(1,number,nb),dcpl_mb(orb(1)))
               call band_kvtimesorb(kv,dcpl_mb(orb(1)))
               call C3dB_cr_imag(1,
     >                           dcpl_mb(orb(1)),
     >                           dbl_mb(rho(1)))
               call band_dplot_gcube_write(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
             endif

           end if

        end if


      end do


*     **** dealocate MA local variables ****
      value =           MA_pop_stack(orb(2))
      value = value.and.MA_pop_stack(rho(2))
      if (.not. value) call errquit('popping of stack memory',0, MA_ERR)


    
      return
      end



 


*     ********************************
*     *                              *
*     *        cpsi_translate        *
*     *                              *
*     ********************************

      subroutine cpsi_translate(trans,npack1,nemax,nbrillioun,psi)
      implicit none
      real*8 trans(3)
      integer npack1,nemax,nbrillioun
      complex*16 psi(npack1,nemax,nbrillioun)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer n,nion,nb
      integer exi(2),rion(2)

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unita
      external ion_nion
      external lattice_unita

      nion  = ion_nion()
      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*nion,'rion', rion(2), rion(1))
      if (.not. value) call errquit('cpsi_translate:push stack',0,0)
      call dscal(3*nion,0.0d0,dbl_mb(rion(1)),1)

      dbl_mb(rion(1)  ) = trans(1)
     >                  - (   lattice_unita(1,1)
     >                      + lattice_unita(1,2)
     >                      + lattice_unita(1,3) )/2.0d0
      dbl_mb(rion(1)+1) = trans(2)
     >                  - (   lattice_unita(2,1)
     >                      + lattice_unita(2,2)
     >                      + lattice_unita(2,3) )/2.0d0
      dbl_mb(rion(1)+2) = trans(3)
     >                  - (   lattice_unita(3,1)
     >                      + lattice_unita(3,2)
     >                      + lattice_unita(3,3) )/2.0d0
      call cphafac_rion(dbl_mb(rion(1)))
 
      do nb=1,nbrillioun
        call cstrfac_pack(nb,1,dcpl_mb(exi(1)))
        do n=1,nemax
c          call Cram_cc_Mul(nb,dcpl_mb(exi(1)),psi(1,n,nb),psi(1,n,nb))
          call Cram_cc_Mul2(nb,dcpl_mb(exi(1)),psi(1,n,nb))
        end do
      end do


      value =           MA_pop_stack(rion(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not. value) call errquit('cpsi_translate:pop stack',1,0)
      return
      end



      subroutine dplot_gen_cpsi_dn(ispin,ne,
     >                      npack1,nfft3d,nemax,nbrill,
     >                      psi,
     >                      dn,
     >                      psi_r)
      implicit none
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax,nbrill
      complex*16 psi(npack1,nemax,nbrill)
      real*8     dn(nfft3d,2)
      complex*16 psi_r(nfft3d,nemax,nbrill)


*     **** local variables ****
      integer taskid
      integer MASTER
      parameter (MASTER=0)

      integer i,ms,n,nbrillioun,neall
      real*8  scal2


*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      
      call Parallel_taskid(taskid)

      neall = ne(1)+ne(2)
      scal2 = 1.0d0/lattice_omega()

*     *******************
*     **** get psi_r ****
*     *******************
      call BGrsm_fh_fftb(psi,psi_r)


*     *********************
*     **** generate dn ****
*     *********************
      call BGrsm_hr_SumSqr(psi_r,dn)
c      call C3dB_r_SMul(1,scal2,dn,dn)
c      if (ispin.gt.1) call C3dB_r_SMul(1,scal2,dn(1,2),dn(1,2))
      call C3dB_r_SMul1(1,scal2,dn)
      if (ispin.gt.1) call C3dB_r_SMul1(1,scal2,dn(1,2))

      return
      end



      subroutine band_dplot_gcube_write(rtdb,filename,
     >                             number,cube_comment,rho)
      implicit none
#include "errquit.fh"
      integer rtdb
      character*30 filename
      integer      number
      character*72 cube_comment
      real*8 rho(*)
     

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer unit_id
      parameter (unit_id = 72)

      integer taskid
      integer MASTER
      parameter (MASTER=0)


*     **** geometry variables ****
c     integer geom1
c     character*16 t
      integer nion,nion2
      real*8 q,rxyz(3)
      real*8 position_tolerance

*     **** lattice variables ****
      integer np1,np2,np3
      real*8 ua(3,3),r0(3)

      integer l,orb_flag,i
      character*255 full_filename
      integer nfft3d,tmp1(2),tmp2(2)

*     **** external functions ****
      integer  ion_nion      
      real*8   lattice_unita,ion_rion,ion_q
      external ion_nion
      external lattice_unita,ion_rion,ion_q

      call parallel_taskid(taskid)

      
*     **** OPEN cube FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         OPEN(unit_id,file=full_filename(1:l),form='formatted')
         WRITE(unit_id,*) 'molecule'
         WRITE(unit_id,*)  cube_comment
      end if

*     **** open geom ****
      nion = ion_nion()
c     value = geom_create(geom1,'geometry')
c     value = value.and.geom_rtdb_load(rtdb,geom1,'geometry')
c     value = value.and.geom_ncent(geom1,nion)
c     if (.not. value) call errquit('opening geometry',0)

*     **** write lattice ****
      call C3dB_nx(1,np1)
      call C3dB_ny(1,np2)
      call C3dB_nz(1,np3)
      do i=1,3
         ua(i,1) = lattice_unita(i,1)/np1
         ua(i,2) = lattice_unita(i,2)/np2
         ua(i,3) = lattice_unita(i,3)/np3

        r0(i) = -( lattice_unita(i,1) 
     >           + lattice_unita(i,2) 
     >           + lattice_unita(i,3) )/2.0d0
      end do

*     **** get position_tolerance and nion2 = nion+special_positions ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0

      nion2 = nion
      do i=1,nion
c        value = geom_cent_get(geom1,i,t,rxyz,q)
c        if (.not. value) call errquit('reading geometry',0)
         rxyz(1) = ion_rion(1,i)
         rxyz(2) = ion_rion(2,i)
         rxyz(3) = ion_rion(3,i)
         call special_position_count(position_tolerance,
     >                                   rxyz,nion2)
      end do


      if (taskid.eq.MASTER) then
         orb_flag = 1
         if (number.gt.0) orb_flag = -1

         write(unit_id,'(I5,3F12.6)') orb_flag*nion2,r0
         write(unit_id,'(I5,3F12.6)') np1,ua(1,1),ua(2,1),ua(3,1)
         write(unit_id,'(I5,3F12.6)') np2,ua(1,2),ua(2,2),ua(3,2)
         write(unit_id,'(I5,3F12.6)') np3,ua(1,3),ua(2,3),ua(3,3)
      end if




*     **** write geometry ****
      do i=1,nion
c        value = geom_cent_get(geom1,i,t,rxyz,q)
c        if (.not. value) call errquit('reading geometry',0)
         rxyz(1) = ion_rion(1,i)
         rxyz(2) = ion_rion(2,i)
         rxyz(3) = ion_rion(3,i)
         q       = ion_q(i)
         if (taskid.eq.MASTER) then
           WRITE(unit_id,'(I5,4F12.6)') int(q),q,rxyz
           call special_position_tolerance(position_tolerance,
     >                                   unit_id,q,rxyz)
         end if
      end do

*     **** write orbital header ****
      if (number.gt.0) then
         if (taskid.eq.MASTER) write(unit_id,*) 1,number 
      end if

*     **** allocate space ****
      call C3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** transpose grid ****
      call C3dB_r_transpose_jk(1,rho,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))

*     **** write grid ****
      call C3dB_r_FormatWrite_reverse(1,unit_id,rho,dcpl_mb(tmp1(1)))

*     **** deallocate space ****
      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)

*     **** close geom ****
c     value = geom_destroy(geom1)
c     if (.not. value) call errquit('closing geometry',0)

*     **** CLOSE cube FILE ****
      if (taskid.eq.MASTER) then
         CLOSE(unit_id)     
      end if

      return
      end 





      subroutine band_dplot_gcube_write3d(rtdb,filename,
     >                             number,cube_comment,rho)
      implicit none
      integer rtdb
      character*30 filename
      integer      number
      character*72 cube_comment
      real*8 rho(*)
     

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer unit_id
      parameter (unit_id = 72)

      integer taskid
      integer MASTER
      parameter (MASTER=0)


*     **** geometry variables ****
c     integer geom1
c     character*16 t
      integer nion,nion2
      real*8 q,rxyz(3),p(3)
      real*8 position_tolerance

*     **** lattice variables ****
      integer np1,np2,np3
      real*8 ua(3,3),r0(3)

      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8 sx,sy,sz,xy,xz,yz,scal1
      real*8 dx,dy,dz,w
      real*8  sizex(2),sizey(2),sizez(2)
      real*8 xaxis(3),yaxis(3),zaxis(3),qq(3)
    

      integer l,orb_flag,i,j,k,ii,jj,kk,indx
      character*255 full_filename
      integer nfft3d,tmp1(2),tmp2(2),tmp3(2)
      integer Gx,Gy,Gz

*     **** external functions ****
      integer  ion_nion      
      real*8   lattice_unita,ion_rion,ion_q
      external ion_nion
      external lattice_unita,ion_rion,ion_q

      integer  c_G_indx
      external c_G_indx


      call parallel_taskid(taskid)

      
*     **** OPEN cube FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         OPEN(unit_id,file=full_filename(1:l),form='formatted')
         WRITE(unit_id,*) 'molecule'
         WRITE(unit_id,*)  cube_comment
      end if

*     **** open geom ****
      nion = ion_nion()

*     **** read the origin ***
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:o',
     >               mt_dbl,3,r0)) then
        r0(1) = 0.0
        r0(2) = 0.0
        r0(3) = 0.0
      end if

*     **** read the x-axis,y-axis,z-axis ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:x',
     >               mt_dbl,3,xaxis)) then
        xaxis(1) = 1.0
        xaxis(2) = 0.0
        xaxis(3) = 0.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:y',
     >               mt_dbl,3,yaxis)) then
        yaxis(1) = 0.0
        yaxis(2) = 1.0
        yaxis(3) = 0.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:z',
     >               mt_dbl,3,zaxis)) then
        zaxis(1) = 0.0
        zaxis(2) = 0.0
        zaxis(3) = 1.0
      end if

*     ** normalization of x **
      sx = dsqrt(xaxis(1)**2+xaxis(2)**2+xaxis(3)**2)
      do i=1,3
        xaxis(i) = xaxis(i) / sx
      end do

*     ** orthogonalization of y wrt x **
      xy = xaxis(1)*yaxis(1)+xaxis(2)*yaxis(2)+xaxis(3)*yaxis(3)
      do i=1,3
        yaxis(i) = yaxis(i) - xy*xaxis(i)
      end do

*     ** normalization of y **
      sy = dsqrt(yaxis(1)**2+yaxis(2)**2+yaxis(3)**2)
      do i=1,3
        yaxis(i) = yaxis(i) / sy
      end do

*     ** orthogonalization of z wrt x **
      xz = xaxis(1)*zaxis(1)+xaxis(2)*zaxis(2)+xaxis(3)*zaxis(3)
      do i=1,3
        zaxis(i) = zaxis(i) - xz*xaxis(i)
      end do

*     ** orthogonalization of z wrt y **
      yz = yaxis(1)*zaxis(1)+yaxis(2)*zaxis(2)+yaxis(3)*zaxis(3)
      do i=1,3
        zaxis(i) = zaxis(i) - yz*yaxis(i)
      end do

*     ** normalization of z **
      sz = dsqrt(zaxis(1)**2+zaxis(2)**2+zaxis(3)**2)
      do i=1,3
        zaxis(i) = zaxis(i) / sz
      end do



*     **** read nx,ny,and nz ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:nx',
     >               mt_int,1,nx))
     >  nx = 32
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:ny',
     >               mt_int,1,ny))
     >  ny = 32
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:nz',
     >               mt_int,1,nz))
     >  nz = 32
       nxh = nx/2
       nyh = ny/2
       nzh = nz/2

*     **** read sizex,sizey,sizez ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:sizex',
     >               mt_dbl,2,sizex)) then
        sizex(1) = 0.0
        sizex(2) = 10.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:sizey',
     >               mt_dbl,2,sizey)) then
        sizey(1) = 0.0
        sizey(2) = 10.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:sizez',
     >               mt_dbl,2,sizez)) then
        sizez(1) = 0.0
        sizez(2) = 10.0
      end if

      if (taskid.eq.MASTER) then
         write(*,*)
         write(*,1240) 
         write(*,1241) r0
         write(*,1242) xaxis,sizex,nx
         write(*,1243) yaxis,sizey,ny
         write(*,1244) zaxis,sizez,nz
         write(*,*)
         write(*,*)
      end if
 1240 FORMAT(5x,'    3d-Grid Generation')
 1241 FORMAT(5x,'    origin=<',3f8.3,' >')
 1242 FORMAT(5x,'    x-axis=<',3f8.3,' >   xmin,xmax=',2F8.3,' nx=',I3)
 1243 FORMAT(5x,'    y-axis=<',3f8.3,' >   ymin,ymax=',2F8.3,' ny=',I3)
 1244 FORMAT(5x,'    z-axis=<',3f8.3,' >   zmin,zmax=',2F8.3,' nz=',I3)


      call C3dB_nx(1,np1)
      call C3dB_ny(1,np2)
      call C3dB_nz(1,np3)

      dx = (sizex(2)-sizex(1))/dble(nx)
      dy = (sizey(2)-sizey(1))/dble(ny)
      dz = (sizez(2)-sizez(1))/dble(nz)
      do i=1,3
         ua(i,1) = dx*xaxis(i)
         ua(i,2) = dy*yaxis(i)
         ua(i,3) = dz*zaxis(i)
      end do

*     ** shift origin **
      qq(1) = r0(1) + 0.5d0*(sizex(1)+sizex(2))
     >       - 0.5d0*(lattice_unita(1,1)
     >               +lattice_unita(1,2)
     >               +lattice_unita(1,3))
      qq(2) = r0(2) + 0.5d0*(sizey(1)+sizey(2))
     >       - 0.5d0*(lattice_unita(2,1)
     >               +lattice_unita(2,2)
     >               +lattice_unita(2,3))
      qq(3) = r0(3) + 0.5d0*(sizez(1)+sizez(2))
     >       - 0.5d0*(lattice_unita(3,1)
     >               +lattice_unita(3,2)
     >               +lattice_unita(3,3))

      r0(1) = r0(1) + sizex(1)
      r0(2) = r0(2) + sizey(1)
      r0(3) = r0(3) + sizez(1)




*     **** get position_tolerance and nion2 = nion+special_positions ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0

      nion2 = nion
      do i=1,nion
c        value = geom_cent_get(geom1,i,t,rxyz,q)
c        if (.not. value) call errquit('reading geometry',0)
         rxyz(1) = ion_rion(1,i)
         rxyz(2) = ion_rion(2,i)
         rxyz(3) = ion_rion(3,i)
         call special_position_count(position_tolerance,
     >                                   rxyz,nion2)
      end do


      if (taskid.eq.MASTER) then
         orb_flag = 1
         if (number.gt.0) orb_flag = -1

         write(unit_id,'(I5,3F12.6)') orb_flag*nion2,r0
         write(unit_id,'(I5,3F12.6)') nx,ua(1,1),ua(2,1),ua(3,1)
         write(unit_id,'(I5,3F12.6)') ny,ua(1,2),ua(2,2),ua(3,2)
         write(unit_id,'(I5,3F12.6)') nz,ua(1,3),ua(2,3),ua(3,3)
      end if


*     **** write geometry ****
      do i=1,nion
c        value = geom_cent_get(geom1,i,t,rxyz,q)
c        if (.not. value) call errquit('reading geometry',0)
         rxyz(1) = ion_rion(1,i)
         rxyz(2) = ion_rion(2,i)
         rxyz(3) = ion_rion(3,i)
         q       = ion_q(i)
         if (taskid.eq.MASTER) then
           WRITE(unit_id,'(I5,4F12.6)') int(q),q,rxyz
           call special_position_tolerance(position_tolerance,
     >                                   unit_id,q,rxyz)
         end if
      end do

*     **** write orbital header ****
      if (number.gt.0) then
         if (taskid.eq.MASTER) write(unit_id,*) 1,number 
      end if

*     **** allocate space ****
      call C3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(nx*ny*nz),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*     **** fft density  ****
      scal1=1.0d0/dble(np1*np2*np3)
      call C3dB_rc_SMul(1,scal1,rho,dcpl_mb(tmp2(1)))
      call C3dB_rc_fft3f(1,dcpl_mb(tmp2(1)))
      call Cram_c_pack(0,dcpl_mb(tmp2(1)))


*     **** density at mesh points ****
      Gx= c_G_indx(1)
      Gy= c_G_indx(2)
      Gz= c_G_indx(3)
      do k=1,nz
      do j=1,ny
      do i=1,nx
        ii = i-1-nxh
        jj = j-1-nyh
        kk = k-1-nzh
        p(1) = ii*ua(1,1) + jj*ua(1,2) + kk*ua(1,3) + qq(1)
        p(2) = ii*ua(2,1) + jj*ua(2,2) + kk*ua(2,3) + qq(2)
        p(3) = ii*ua(3,1) + jj*ua(3,2) + kk*ua(3,3) + qq(3)
        do l=1,nfft3d
          w = dbl_mb(Gx+l-1)*p(1) 
     >      + dbl_mb(Gy+l-1)*p(2) 
     >      + dbl_mb(Gz+l-1)*p(3)
          dcpl_mb(tmp1(1)+l-1)=dcmplx(dcos(w),-dsin(w))
        end do
        call Cram_c_pack(0,dcpl_mb(tmp1(1)))
        call Cram_cc_dot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)),w)
        indx = (i-1) + (j-1)*nx + (k-1)*nx*ny
        dbl_mb(tmp3(1)+indx) = w
      end do
      end do
      end do

*     **** write grid ****
      if (taskid.eq.MASTER) then
      do i=1,nx
      do j=1,ny
        write(unit_id,'(6e13.5)') 
     >    (dbl_mb(tmp3(1)+(i-1) 
     >                   +(j-1)*nx 
     >                   +(k-1)*nx*ny),k=1,nz)
      end do
      end do
      end if

*     **** deallocate space ****
      value = MA_pop_stack(tmp3(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)

   

*     **** close geom ****
c     value = geom_destroy(geom1)
c     if (.not. value) call errquit('closing geometry',0)



*     **** CLOSE cube FILE ****
      if (taskid.eq.MASTER) then
         CLOSE(unit_id)     
      end if

      return
      end 




      subroutine band_dplot_gcube_write3d_orb(rtdb,filename,
     >                                        number,cube_comment,
     >                                        orb,kv,isreal)
      implicit none
      integer rtdb
      character*30 filename
      integer      number
      character*72 cube_comment
      complex*16   orb(*)
      real*8       kv(3)
      logical      isreal
     

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer unit_id
      parameter (unit_id = 72)

      integer taskid
      integer MASTER
      parameter (MASTER=0)


*     **** geometry variables ****
c     integer geom1
c     character*16 t
      integer nion,nion2
      real*8 q,rxyz(3),p(3)
      real*8 position_tolerance

*     **** lattice variables ****
      integer np1,np2,np3
      real*8 ua(3,3),r0(3)

      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8 sx,sy,sz,xy,xz,yz,scal1
      real*8 dx,dy,dz,w
      real*8  sizex(2),sizey(2),sizez(2)
      real*8 xaxis(3),yaxis(3),zaxis(3),qq(3)
    

      integer l,orb_flag,i,j,k,ii,jj,kk,indx
      character*255 full_filename
      integer nfft3d,tmp1(2),tmp2(2),tmp3(2)
      integer Gx,Gy,Gz
      complex*16 phafac,wc

*     **** external functions ****
      integer  ion_nion      
      real*8   lattice_unita,ion_rion,ion_q
      external ion_nion
      external lattice_unita,ion_rion,ion_q

      integer  c_G_indx
      external c_G_indx


      call parallel_taskid(taskid)

      
*     **** OPEN cube FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         OPEN(unit_id,file=full_filename(1:l),form='formatted')
         WRITE(unit_id,*) 'molecule'
         WRITE(unit_id,*)  cube_comment
      end if

*     **** open geom ****
      nion = ion_nion()

*     **** read the origin ***
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:o',
     >               mt_dbl,3,r0)) then
        r0(1) = 0.0
        r0(2) = 0.0
        r0(3) = 0.0
      end if

*     **** read the x-axis,y-axis,z-axis ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:x',
     >               mt_dbl,3,xaxis)) then
        xaxis(1) = 1.0
        xaxis(2) = 0.0
        xaxis(3) = 0.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:y',
     >               mt_dbl,3,yaxis)) then
        yaxis(1) = 0.0
        yaxis(2) = 1.0
        yaxis(3) = 0.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:z',
     >               mt_dbl,3,zaxis)) then
        zaxis(1) = 0.0
        zaxis(2) = 0.0
        zaxis(3) = 1.0
      end if

*     ** normalization of x **
      sx = dsqrt(xaxis(1)**2+xaxis(2)**2+xaxis(3)**2)
      do i=1,3
        xaxis(i) = xaxis(i) / sx
      end do

*     ** orthogonalization of y wrt x **
      xy = xaxis(1)*yaxis(1)+xaxis(2)*yaxis(2)+xaxis(3)*yaxis(3)
      do i=1,3
        yaxis(i) = yaxis(i) - xy*xaxis(i)
      end do

*     ** normalization of y **
      sy = dsqrt(yaxis(1)**2+yaxis(2)**2+yaxis(3)**2)
      do i=1,3
        yaxis(i) = yaxis(i) / sy
      end do

*     ** orthogonalization of z wrt x **
      xz = xaxis(1)*zaxis(1)+xaxis(2)*zaxis(2)+xaxis(3)*zaxis(3)
      do i=1,3
        zaxis(i) = zaxis(i) - xz*xaxis(i)
      end do

*     ** orthogonalization of z wrt y **
      yz = yaxis(1)*zaxis(1)+yaxis(2)*zaxis(2)+yaxis(3)*zaxis(3)
      do i=1,3
        zaxis(i) = zaxis(i) - yz*yaxis(i)
      end do

*     ** normalization of z **
      sz = dsqrt(zaxis(1)**2+zaxis(2)**2+zaxis(3)**2)
      do i=1,3
        zaxis(i) = zaxis(i) / sz
      end do



*     **** read nx,ny,and nz ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:nx',
     >               mt_int,1,nx))
     >  nx = 32
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:ny',
     >               mt_int,1,ny))
     >  ny = 32
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:nz',
     >               mt_int,1,nz))
     >  nz = 32
       nxh = nx/2
       nyh = ny/2
       nzh = nz/2

*     **** read sizex,sizey,sizez ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:sizex',
     >               mt_dbl,2,sizex)) then
        sizex(1) = 0.0
        sizex(2) = 10.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:sizey',
     >               mt_dbl,2,sizey)) then
        sizey(1) = 0.0
        sizey(2) = 10.0
      end if
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:3d_grid:sizez',
     >               mt_dbl,2,sizez)) then
        sizez(1) = 0.0
        sizez(2) = 10.0
      end if

      if (taskid.eq.MASTER) then
         write(*,*)
         write(*,1240) 
         write(*,1241) r0
         write(*,1242) xaxis,sizex,nx
         write(*,1243) yaxis,sizey,ny
         write(*,1244) zaxis,sizez,nz
         write(*,*)
         write(*,*)
      end if
 1240 FORMAT(5x,'    3d-Grid Generation')
 1241 FORMAT(5x,'    origin=<',3f8.3,' >')
 1242 FORMAT(5x,'    x-axis=<',3f8.3,' >   xmin,xmax=',2F8.3,' nx=',I3)
 1243 FORMAT(5x,'    y-axis=<',3f8.3,' >   ymin,ymax=',2F8.3,' ny=',I3)
 1244 FORMAT(5x,'    z-axis=<',3f8.3,' >   zmin,zmax=',2F8.3,' nz=',I3)


      call C3dB_nx(1,np1)
      call C3dB_ny(1,np2)
      call C3dB_nz(1,np3)

      dx = (sizex(2)-sizex(1))/dble(nx)
      dy = (sizey(2)-sizey(1))/dble(ny)
      dz = (sizez(2)-sizez(1))/dble(nz)
      do i=1,3
         ua(i,1) = dx*xaxis(i)
         ua(i,2) = dy*yaxis(i)
         ua(i,3) = dz*zaxis(i)
      end do

*     ** shift origin **
      qq(1) = r0(1) + 0.5d0*(sizex(1)+sizex(2))
     >       - 0.5d0*(lattice_unita(1,1)
     >               +lattice_unita(1,2)
     >               +lattice_unita(1,3))
      qq(2) = r0(2) + 0.5d0*(sizey(1)+sizey(2))
     >       - 0.5d0*(lattice_unita(2,1)
     >               +lattice_unita(2,2)
     >               +lattice_unita(2,3))
      qq(3) = r0(3) + 0.5d0*(sizez(1)+sizez(2))
     >       - 0.5d0*(lattice_unita(3,1)
     >               +lattice_unita(3,2)
     >               +lattice_unita(3,3))

      r0(1) = r0(1) + sizex(1)
      r0(2) = r0(2) + sizey(1)
      r0(3) = r0(3) + sizez(1)




*     **** get position_tolerance and nion2 = nion+special_positions ****
      if (.not.
     >      rtdb_get(rtdb,'band_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0

      nion2 = nion
      do i=1,nion
c        value = geom_cent_get(geom1,i,t,rxyz,q)
c        if (.not. value) call errquit('reading geometry',0)
         rxyz(1) = ion_rion(1,i)
         rxyz(2) = ion_rion(2,i)
         rxyz(3) = ion_rion(3,i)
         call special_position_count(position_tolerance,
     >                                   rxyz,nion2)
      end do


      if (taskid.eq.MASTER) then
         orb_flag = 1
         if (number.gt.0) orb_flag = -1

         write(unit_id,'(I5,3F12.6)') orb_flag*nion2,r0
         write(unit_id,'(I5,3F12.6)') nx,ua(1,1),ua(2,1),ua(3,1)
         write(unit_id,'(I5,3F12.6)') ny,ua(1,2),ua(2,2),ua(3,2)
         write(unit_id,'(I5,3F12.6)') nz,ua(1,3),ua(2,3),ua(3,3)
      end if


*     **** write geometry ****
      do i=1,nion
c        value = geom_cent_get(geom1,i,t,rxyz,q)
c        if (.not. value) call errquit('reading geometry',0)
         rxyz(1) = ion_rion(1,i)
         rxyz(2) = ion_rion(2,i)
         rxyz(3) = ion_rion(3,i)
         q       = ion_q(i)
         if (taskid.eq.MASTER) then
           WRITE(unit_id,'(I5,4F12.6)') int(q),q,rxyz
           call special_position_tolerance(position_tolerance,
     >                                   unit_id,q,rxyz)
         end if
      end do

*     **** write orbital header ****
      if (number.gt.0) then
         if (taskid.eq.MASTER) write(unit_id,*) 1,number 
      end if

*     **** allocate space ****
      call C3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(nx*ny*nz),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*     **** fft density  ****
      scal1=1.0d0/dble(np1*np2*np3)
c      call C3dB_c_SMul(1,scal1,orb,dcpl_mb(tmp2(1)))
      call dcopy(2*nfft3d,orb,1,dcpl_mb(tmp2(1)),1)
      call dscal(2*nfft3d,scal1,dcpl_mb(tmp2(1)),1)
      call C3dB_rc_fft3f(1,dcpl_mb(tmp2(1)))
      call Cram_c_pack(0,dcpl_mb(tmp2(1)))


*     **** density at mesh points ****
      Gx= c_G_indx(1)
      Gy= c_G_indx(2)
      Gz= c_G_indx(3)
      do k=1,nz
      do j=1,ny
      do i=1,nx
        ii = i-1-nxh
        jj = j-1-nyh
        kk = k-1-nzh
        p(1) = ii*ua(1,1) + jj*ua(1,2) + kk*ua(1,3) + qq(1)
        p(2) = ii*ua(2,1) + jj*ua(2,2) + kk*ua(2,3) + qq(2)
        p(3) = ii*ua(3,1) + jj*ua(3,2) + kk*ua(3,3) + qq(3)
        do l=1,nfft3d
          w = dbl_mb(Gx+l-1)*p(1) 
     >      + dbl_mb(Gy+l-1)*p(2) 
     >      + dbl_mb(Gz+l-1)*p(3)
          dcpl_mb(tmp1(1)+l-1)=dcmplx(dcos(w),-dsin(w))
        end do
        call Cram_c_pack(0,dcpl_mb(tmp1(1)))
        call Cram_cc_zdot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)),wc)
        w = kv(1)*p(1) + kv(2)*p(2) + kv(3)*p(3)
        phafac = dcmplx(dcos(w),dsin(w))
        wc = phafac*wc
        indx = (i-1) + (j-1)*nx + (k-1)*nx*ny
        if (isreal) then
           dbl_mb(tmp3(1)+indx) = dreal(wc)
        else
           dbl_mb(tmp3(1)+indx) = dimag(wc)
        end if
      end do
      end do
      end do

*     **** write grid ****
      if (taskid.eq.MASTER) then
      do i=1,nx
      do j=1,ny
        write(unit_id,'(6e13.5)') 
     >    (dbl_mb(tmp3(1)+(i-1) 
     >                   +(j-1)*nx 
     >                   +(k-1)*nx*ny),k=1,nz)
      end do
      end do
      end if

*     **** deallocate space ****
      value = MA_pop_stack(tmp3(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)

   

*     **** close geom ****
c     value = geom_destroy(geom1)
c     if (.not. value) call errquit('closing geometry',0)



*     **** CLOSE cube FILE ****
      if (taskid.eq.MASTER) then
         CLOSE(unit_id)     
      end if

      return
      end 


      subroutine band_kvtimesorb(kv,orb)
      implicit none
      real*8     kv(3)
      complex*16 orb(*)

*     *** local variables ***
      integer    taskid,index,p
      integer    i,j,k,k1,k2,k3
      integer    nx,ny,nz,nhx,nhy,nhz
      real*8     a(3,3),w
      complex*16 phafac

*     *** external functions ***
      real*8   lattice_unita
      external lattice_unita

      call Parallel_taskid(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)

      nhx = nx/2
      nhy = ny/2
      nhz = nz/2

      do i=1,3
         a(i,1) = lattice_unita(i,1)/nx
         a(i,2) = lattice_unita(i,2)/ny
         a(i,3) = lattice_unita(i,3)/nz
      end do

      do k3=-nhz,nhz-1
      do k2=-nhy,nhy-1
      do k1=-nhx,nhx-1
        i = k1 + nhx
        j = k2 + nhy
        k = k3 + nhz
        call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
        if (p.eq.taskid) then
           w = kv(1)*(a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3)
     >       + kv(2)*(a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3)
     >       + kv(3)*(a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3)
           phafac = dcmplx(dcos(w),dsin(w))
           orb(index) = orb(index)*phafac
        end if
      end do
      end do
      end do

      return
      end
