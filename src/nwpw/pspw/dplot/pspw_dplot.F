*
* $Id$
*

*     ***********************************************
*     *                                             *
*     *          pspw_dplot                         *
*     *                                             *
*     ***********************************************

      logical function pspw_dplot(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"


      logical value,fractional
      integer taskid,np
      integer MASTER
      parameter (MASTER=0)

*     **** local variables ***

      integer count,smearoccupation
      integer ngrid(3)
      integer nfft3d,n2ft3d,npack1,mapping,mapping1d
      integer ispin,ne(2),nemax
      integer psi2(2),dn(2),psir(2),occ2(2),ncell(3)
      real*8  cpu1,cpu2
      real*8  position_tolerance,origin(3)
      character*8 tag
      integer psir2(2)

*     **** external functions ****
      logical      control_read,ion_init,control_balance
      integer      pack_nwave_all,control_mapping,control_np_orbital
      integer      control_ngrid,pack_nwave,control_version
      integer      control_mapping1d
      real*8       lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8       lattice_unitg
      character*12 control_boundry

      external     control_read,ion_init,control_balance
      external     pack_nwave_all,control_mapping,control_np_orbital
      external     control_ngrid,pack_nwave,control_version
      external     control_mapping1d
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
     >      btdb_get(rtdb,'pspw_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0

*     **** get ncell ****
      if (.not.btdb_get(rtdb,'pspw_dplot:ncell',mt_int,3,ncell)) then
        ncell(1) = 0
        ncell(2) = 0
        ncell(3) = 0
      end if

*     **** get origin ****
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:origin',
     >               mt_dbl,3,origin)) then
        origin(1) = 0.0d0 
        origin(2) = 0.0d0 
        origin(3) = 0.0d0 
      end if

*     **** read control file ****
      value = control_read(4,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      mapping   = control_mapping()


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
      call Pack_init()
      call Pack_npack(1,npack1)

*     **** initialize ion data structure ****
      value = ion_init(rtdb)

*     ***** allocate psi2 wavefunctions ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
         fractional = .true.
      else
         fractional = .false.
      end if
      nemax = ne(1)+ne(2)
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      if (fractional) then
         value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ2',occ2(2),occ2(1))
      end if

      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      mapping1d   = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)

*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)),
     >              smearoccupation,dbl_mb(occ2(1)))

*     **** allocate other variables *****
      value = value.and.
     >        BA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir',psir(2),psir(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir2',psir2(2),psir2(1))

      if (.not. value) call errquit('pspw_dplot:out of heap memory',0,
     &       MA_ERR)

      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      call phafac()

*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
         write(6,1111) np
         if (mapping.eq.1) write(6,1112)
         if (mapping.eq.2) write(6,1113)
         if (mapping.eq.3) write(6,1118)
         if (control_balance()) then
           write(6,1114)
         else
           write(6,1116)
         end if
         write(6,1115)
         write(6,1121) control_boundry(),control_version()
         write(6,1220) ne(1),ne(ispin),' ( Fourier space)'
         write(6,1224) ncell
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
         write(6,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(6,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)
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

        call phafac()  !** recompute phase factors **

      end if


      call dplot_gen_psi_dn(ispin,ne,
     >                npack1,nfft3d,nemax,
     >                dcpl_mb(psi2(1)),
     >                dbl_mb(dn(1)),
     >                dcpl_mb(psir(1)),
     >                fractional,dbl_mb(occ2(1)))

c  used to localize psi....input needs to be generated for this capabilitiy
c      call psi_dmatrix_localize(ispin,ne,ne,n2ft3d,
c     >                          dcpl_mb(psir(1)),dcpl_mb(psir2(1)))

      call dplot_loop(rtdb,
     >                ispin,ne,
     >                npack1,nfft3d,nemax,
     >                dcpl_mb(psi2(1)),
     >                dbl_mb(dn(1)),
     >                dcpl_mb(psir(1)),
     >                .false.,tag)

      value = value.and.btdb_get(rtdb,'pspw_dplot:count',
     >                            mt_int,1,count)


*     **** deallocate heap memory ****
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call ion_write(rtdb)  !*** can also use call ion_destroy()????
      !call ion_destroy(rtdb)  
      call ion_end()
      call mask_end()
      call Pack_end()
      call G_end()

      value = BA_free_heap(psir2(2))
      value = value.and.BA_free_heap(psir(2))
      value = value.and.BA_free_heap(dn(2))
      value = value.and.BA_free_heap(psi2(2))
      if (fractional) then
         value = BA_free_heap(occ2(2))
      end if

      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call Dne_end()

      if (.not. value) call errquit('pspw_dplot:freeing heap memory',0,
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
      call Parallel_Finalize()

      pspw_dplot = value

      return 

*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*                   pspw DPLOT                     *')
 1030 FORMAT(10x,'*    [ Generates density and orbital grids  ]      *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #1.00   08/22/01              *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'============ PSPW DPLOT input data =================')
 1111 FORMAT(/' number of processors used:',I3)
 1112 FORMAT( ' parallel mapping         :   1d slab')
 1113 FORMAT( ' parallel mapping         :2d hilbert')
 1114 FORMAT( ' parallel mapping         :  balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1118 FORMAT( ' parallel mapping         : 2d hcurve')
 1121 FORMAT(5X,' boundary conditions   = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
 1224 FORMAT(/' ncell              = ',3I2)
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

 1250 FORMAT(5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')

      end

*     ***********************************************
*     *                                             *
*     *          dplot_gen_psi_dn                   *
*     *                                             *
*     ***********************************************
 
      subroutine dplot_gen_psi_dn(ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi,
     >                      dn,
     >                      psi_r, 
     >                      fractional,occ)
      implicit none
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      real*8     psi_r(2*nfft3d,nemax)
      logical fractional
      real*8  occ(*)


*     **** local variables ****
      integer taskid
      integer MASTER
      parameter (MASTER=0)

      integer n2ft3d
      integer i,ms,n
      integer n1(2),n2(2)
      integer nx,ny,nz
      real*8  scal1,scal2


*     **** external functions ****
      integer  control_version
      real*8   lattice_omega
      external control_version
      external lattice_omega

      
      call Parallel_taskid(taskid)
      n2ft3d = 2*nfft3d

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()

*     *******************
*     **** get psi_r ****
*     *******************
      do n=n1(1),n2(ispin)
         call Pack_c_Copy(1,psi(1,n),psi_r(1,n))
         call Pack_c_unpack(1,psi_r(1,n))
         call D3dB_cr_fft3b(1,psi_r(1,n))
         call D3dB_r_Zero_Ends(1,psi_r(1,n))
c         call D3dB_r_SMul(1,dsqrt(scal2),psi_r(1,n),psi_r(1,n))
         call D3dB_r_SMul1(1,dsqrt(scal2),psi_r(1,n))
      end do


*     *********************
*     **** generate dn ****
*     *********************
      call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
      if (fractional) then
         do ms=1,ispin
            do n=n1(ms),n2(ms)
               do i=1,n2ft3d
c                  dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
                  dn(i,ms) = dn(i,ms) + (psi_r(i,n)**2)*occ(n)
               end do
            end do
            call D3dB_r_Zero_Ends(1,dn(1,ms))
         end do
      else
         do ms=1,ispin
            do n=n1(ms),n2(ms)
               do i=1,n2ft3d
c                  dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
                  dn(i,ms) = dn(i,ms) + (psi_r(i,n)**2)
               end do
            end do
            call D3dB_r_Zero_Ends(1,dn(1,ms))
         end do
      end if

      return
      end

*     ***********************************************
*     *                                             *
*     *          dplot_loop                         *
*     *                                             *
*     ***********************************************

      subroutine dplot_loop(rtdb,
     >                      ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi,
     >                      dn,
     >                      psi_r,
     >                      add_tag,tag)
      implicit none
      integer    rtdb
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      real*8     psi_r(2*nfft3d,nemax)
      logical add_tag
      character*8 tag

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,grid3d,grid1d
      integer taskid
      integer MASTER
      parameter (MASTER=0)

      integer n2ft3d
      integer i,i1,i2,i3,count,number,ia,number1,number2,ms
      integer n1(2),n2(2),nn(2)
      integer nx,ny,nz,indx,indxr
      integer rho(2),rhog(2)
      real*8  scal1,scal2,x(3),r,c

      character*72 cube_comment
      character*50 name1,name2
      character*50 filename
      integer name1_len,name2_len,ind,ind2

      integer n_truncate,idx_truncate(2),rgrid(2),trunc(2)

*     **** external functions ****
      integer  control_version
      real*8   lattice_omega
      external control_version
      external lattice_omega

      
      grid3d = .false.
      if (btdb_get(rtdb,'pspw_dplot:3d_grid:nx',mt_int,1,i))
     >  grid3d = .true.

      grid1d = .false.
      if (btdb_get(rtdb,'pspw_dplot:1d_grid:nx',mt_int,1,i))
     >  grid1d = .true.

      call Parallel_taskid(taskid)
      n2ft3d = 2*nfft3d
      value = BA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()

*     **** setup atom_truncate ****
      n_truncate = 0
      if (btdb_get(rtdb,'pspw_dplot:atom_truncate_size',
     >             mt_int,1,n_truncate)) then
         if (n_truncate.gt.0) then
         if (BA_push_get(mt_dbl,n2ft3d,'trunc',trunc(2),trunc(1))) then
           if (BA_push_get(mt_int,n_truncate,'idx_truncate',
     >                     idx_truncate(2),idx_truncate(1))) then
               if (.not.btdb_get(rtdb,'pspw_dplot:atom_truncate',
     >                  mt_int,n_truncate,int_mb(idx_truncate(1)))) then
                  n_truncate = 0
               else
                  if (.not.BA_push_get(mt_dbl,3*n2ft3d,
     >                                 'rgrd',rgrid(2),rgrid(1)))
     >               call errquit('pspw_dplot_loop:push stack',0,MA_ERR)
                  call lattice_r_grid(dbl_mb(rgrid(1)))
                  call Truncating_Function_init(rtdb)
                  call Truncating_Function_index(n_truncate,
     >                  int_mb(idx_truncate(1)),dbl_mb(rgrid(1)),
     >                  dbl_mb(trunc(1)))
                  call Truncating_Function_end()
                  if (.not.BA_pop_stack(rgrid(2)))
     >               call errquit('pspw_dplot_loop:pop stack',0,MA_ERR)
               end if
               if (.not.BA_pop_stack(idx_truncate(2)))
     >             call errquit('error popping idx_truncate',0,MA_ERR)
           end if
         end if
         end if
      end if


*     ********************************************
*     **** loop over orbital and density list ****
*     ********************************************

      value = value.and.btdb_get(rtdb,'pspw_dplot:count',
     >                            mt_int,1,count)
      do i=1,count

*       **** define name  - not very elegent and could break if ****
*       ****                      count becomes very large      ****
        ia = ICHAR('a')
        name1 = 'pspw_dplot:filename'//CHAR(i-1+ia)
        name2 = 'pspw_dplot:number'//CHAR(i-1+ia)
        name1_len = index(name1,' ') - 1
        name2_len = index(name2,' ') - 1

        value = btdb_cget(rtdb,name1(1:name1_len),1,filename)
        ind = index(filename,' ') - 1
        if (add_tag) then
          ind2 = index(tag,' ') - 1
          filename = filename(1:ind)//tag(1:ind2)//'.cube'
          ind = index(filename,' ') - 1
        end if

        if(.not.btdb_get(rtdb,name2(1:name2_len),mt_int,1,number)) then
           number = 99999999
           name1 = 'pspw_dplot:number1'//CHAR(i-1+ia)
           name2 = 'pspw_dplot:number2'//CHAR(i-1+ia)
           name1_len = index(name1,' ') - 1
           name2_len = index(name2,' ') - 1
           value=value.and.btdb_get(rtdb,name1(1:name1_len),mt_int,1,
     >                              number1)
           value=value.and.btdb_get(rtdb,name2(1:name2_len),mt_int,1,
     >                              number2)
        end if
        if (.not.value)
     >     call errquit(
     >     'pspw_dplot: btdb_get failed for orbital', 0, RTDB_ERR)

*       **** outputing density ****
        if (number.lt.0) then
           number = -number

          goto (710,720,730,740,750,760,770,780,785,786 ) number
          call errquit(
     >      'dplot_loop: unimplemented directive', number, INPUT_ERR)
       
*          *************
*          *** total ***
*          *************
 710       call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
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
 720       call D3dB_rr_Sub(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
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
 730       call D3dB_r_Copy(1,dn(1,1),dbl_mb(rho(1)))
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
 740       call D3dB_r_Copy(1,dn(1,ispin),dbl_mb(rho(1)))
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
 750       call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
c           call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dbl_mb(rho(1)))
           call D3dB_r_SMul1(1,scal1,dbl_mb(rho(1)))
           call D3dB_rc_fft3f(1,dbl_mb(rho(1)))
           call Pack_c_pack(1,dbl_mb(rho(1)))
           nn(1) = 1
           nn(2) = 0
           call ke(1,nn,dbl_mb(rho(1)),dbl_mb(rho(1)))
           call Pack_c_unpack(1,dbl_mb(rho(1)))
           call D3dB_cr_fft3b(1,dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing laplacian density',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF Laplacian Density"
           goto 790

*          *******************************
*          *** Electrostatic Potential ***
*          *******************************
 760       call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
           if (control_version().eq.3) then
             call D3dB_rc_fft3f(1,dbl_mb(rho(1)))
             call Pack_c_pack(0,dbl_mb(rho(1)))
c             call Pack_c_SMul(0,scal1,dbl_mb(rho(1)),dbl_mb(rho(1)))
             call Pack_c_SMul1(0,scal1,dbl_mb(rho(1)))
             call pspw_add_core_dng(1.0d0,dbl_mb(rho(1)))
c             call Pack_c_SMul(0,scal2,dbl_mb(rho(1)),dbl_mb(rho(1)))
             call Pack_c_SMul1(0,scal2,dbl_mb(rho(1)))
             call coulomb_v(dbl_mb(rho(1)),dbl_mb(rho(1)))
c             call Pack_c_SMul(0,(1.0d0/scal2),
c     >                        dbl_mb(rho(1)),dbl_mb(rho(1)))
             call Pack_c_unpack(0,dbl_mb(rho(1)))
             call D3dB_cr_fft3b(1,dbl_mb(rho(1)))
           else
             call coulomb2_v(dbl_mb(rho(1)),dbl_mb(rho(1)))
             call pspw_add_core_pot(1.0d0,dbl_mb(rho(1)))
           end if
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing electrostatic potential',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF Electrostatic Potential"
           goto 790

*          ********************
*          ***      ELF     ***
*          ********************
 770       call generate_ELF(npack1,ne(1),
     >                       psi,dn,dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing restricted/up spin ELF',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF ELF"
           goto 790

 780       call generate_ELF(npack1,ne(2),
     >                       psi(1,ne(1)+1),dn(1,2),dbl_mb(rho(1)))
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing down spin ELF',
     >                  ' to filename: ',filename(1:ind)
           end if
           end if
           cube_comment = "SCF ELF"
           goto 790

*          ***************************
*          ***      density matrix ***
*          ***************************
 785       name1 = 'pspw_dplot:dmatrix_ms'//CHAR(i-1+ia)
           name2 = 'pspw_dplot:dmatrix_xyz'//CHAR(i-1+ia)
           name1_len = index(name1,' ') - 1
           name2_len = index(name2,' ') - 1
           if (.not.btdb_get(rtdb,name1(1:name1_len),mt_int,1,ms))
     >        ms = 1
           if (.not.btdb_get(rtdb,name2(1:name2_len),mt_dbl,3,x)) then
              x(1) = 0.0d0
              x(2) = 0.0d0
              x(3) = 0.0d0
           end if
           if (.not.add_tag) then
           if (taskid.eq.MASTER) then
             write(*,*) '   writing density_matrix(ms,:,(x,y,z)',
     >                  ' to filename: ',filename(1:ind)
             write(*,*) "ms   =",ms
             write(*,'(A,3F10.3)') "x,y,z=",x
           end if
           end if
           call  generate_dmatrix_column(ms,x,
     >                                   ispin,ne,n2ft3d,
     >                                   psi_r,dbl_mb(rho(1)))
           cube_comment = 'SCF density matrix'
           goto 790

*          ***************************
*          ***   structure factor  ***
*          ***************************
 786       if(.not.BA_push_get(mt_dcpl,(nfft3d),'rhog',rhog(2),rhog(1)))
     >        call errquit('out of stack memory',0, MA_ERR)

           call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dcpl_mb(rhog(1)))
           call D3dB_rc_fft3f(1,dcpl_mb(rhog(1)))
           do i3=1,nz
           do i2=1,ny
           do i1=1,(nx/2+1)
              indx = (i1-1) + (nx/2+1)*(i2-1) + (nx/2+1)*ny*(i3-1)
              r = dble(dcpl_mb(rhog(1)+indx))
              c = dimag(dcpl_mb(rhog(1)+indx))
              r = dsqrt(r*r + c*c)
              indxr = (i1-1) + (nx+2)*(i2-1) + (nx+2)*ny*(i3-1)
              dbl_mb(rho(1)+indxr) = r
              if (i>1) then
                 indxr = (nx-i1) + (nx+2)*(i2-1) + (nx+2)*ny*(i3-1)
                 dbl_mb(rho(1)+indxr) = r
              end if
           end do
           end do
           end do
           cube_comment = 'Structure Factor'

           if(.not.BA_pop_stack(rhog(2)))
     >        call errquit('popping stack memory',0, MA_ERR)

           goto 790


 790       continue
           number = -number

*       **** outputing wavefunction or wavefunction2 ****
        else
           !*** used to get spin down homo or lumo ****
           if (number.eq.123456789) number = ne(1)+1

           !*** overlap functions ****
           if (number.eq.99999999) then
              number = -99
              if ((number1.gt.nemax).or.(number1.lt.1)) then
              if (taskid.eq.MASTER) then
                write(*,*)  '   Bad orbital number1 ', number1,
     >                     ', changing to orbital number1 1.'
              end if
              number1 = 1
              end if
              if ((number2.gt.nemax).or.(number2.lt.1)) then
              if (taskid.eq.MASTER) then
                write(*,*)  '   Bad orbital number2 ', number2,
     >                     ', changing to orbital number2 1.'
              end if
              number2 = 1
              end if
              call D3dB_rr_Mul(1,psi_r(1,number1),psi_r(1,number2),
     >                         dbl_mb(rho(1)))
              if (.not.add_tag) then
              if (taskid.eq.MASTER) then
                write(*,*) '   writing orbital2 ',number1, number2,
     >                     ' to filename: ',filename(1:ind)
              end if
              end if
              cube_comment = 'SCF Molecular Orbitals squared'

           !*** regular molecular orbital ****
           else
              if ((number.gt.nemax).or.(number.lt.1)) then
              if (taskid.eq.MASTER) then
                write(*,*)  '   Bad orbital number ', number,
     >                     ', changing to orbital number 1.'
              end if
              number = 1
              end if
          
              call D3dB_r_Copy(1,psi_r(1,number),dbl_mb(rho(1)))
              if (.not.add_tag) then
              if (taskid.eq.MASTER) then
                write(*,*) '   writing orbital ',number, 
     >                     ' to filename: ',filename(1:ind)
              end if
              end if
              cube_comment = 'SCF Molecular Orbitals'
           end if
        end if

*       *** atom_truncate ***
        if (n_truncate.gt.0) then
           call D3dB_rr_Mul2(1,dbl_mb(trunc(1)),dbl_mb(rho(1)))
        end if

        if (grid3d) then
          call dplot_gcube_write3d(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
        else if (grid1d) then
          call dplot_gcube_write1d(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
        else
          call dplot_gcube_write(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))
        endif

      end do

*     **** dealocate MA local variables ****
      if (n_truncate.gt.0) then
         if (.not.BA_pop_stack(trunc(2)))
     >      call errquit('error popping trunc',0,MA_ERR)
      end if

      value = BA_pop_stack(rho(2))
      if (.not. value) call errquit('popping of stack memory',0, MA_ERR)

    
      return
      end



*     ***********************************************
*     *                                             *
*     *          dplot_gcube_write                  *
*     *                                             *
*     ***********************************************

      subroutine dplot_gcube_write(rtdb,filename,
     >                             number,cube_comment,rho)
      implicit none
      integer rtdb
      character*50 filename
      integer      number
      character*72 cube_comment
      real*8 rho(*)
     

#include "bafdecls.fh"
#include "btdb.fh"
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
      integer nion,nion2,ncell(3),nl(3),nh(3),n1,n2,n3
      real*8 q,rxyz(3)
      real*8 position_tolerance

*     **** lattice variables ****
      integer np1,np2,np3
      real*8 ua(3,3),r0(3)

      integer l,orb_flag,i
      character*255 full_filename
      integer nfft3d,tmp1(2),tmp2(2)

*     **** external functions ****
      integer  ion_nion,control_mapping
      real*8   lattice_unita,ion_rion,ion_q
      external ion_nion,control_mapping
      external lattice_unita,ion_rion,ion_q

      call Parallel_taskid(taskid)

      
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
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
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
     >      btdb_get(rtdb,'pspw_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0

*     **** get ncell ****
      if (.not.btdb_get(rtdb,'pspw_dplot:ncell',mt_int,3,ncell)) then
        ncell(1) = 0
        ncell(2) = 0
        ncell(3) = 0
      end if
      if (ncell(1).eq.0) then
         nl(1) = 0
         nh(1) = 0
      else if (ncell(1).gt.0) then
         nl(1) = -ncell(1)
         nh(1) =  ncell(1)
      else
         nl(1) =  0
         nh(1) = -ncell(1)
      end if
      if (ncell(2).eq.0) then
         nl(2) = 0
         nh(2) = 0
      else if (ncell(2).gt.0) then
         nl(2) = -ncell(2)
         nh(2) =  ncell(2)
      else
         nl(2) =  0
         nh(2) = -ncell(2)
      end if
      if (ncell(3).eq.0) then
         nl(3) = 0
         nh(3) = 0
      else if (ncell(3).gt.0) then
         nl(3) = -ncell(3)
         nh(3) =  ncell(3)
      else
         nl(3) =  0
         nh(3) = -ncell(3)
      end if
      nion2 = 0
      do n3=nl(3),nh(3)
      do n2=nl(2),nh(2)
      do n1=nl(1),nh(1)
         nion2= nion2+nion
      end do
      end do
      end do


c      nion2 = nion
c      do i=1,nion
cc        value = geom_cent_get(geom1,i,t,rxyz,q)
cc        if (.not. value) call errquit('reading geometry',0)
c         rxyz(1) = ion_rion(1,i)
c         rxyz(2) = ion_rion(2,i)
c         rxyz(3) = ion_rion(3,i)
c         call special_position_count(position_tolerance,
c     >                                   rxyz,nion2)
c      end do

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
         q = ion_q(i)
         do n3=nl(3),nh(3)
         do n2=nl(2),nh(2)
         do n1=nl(1),nh(1)
            rxyz(1) = ion_rion(1,i) + n1*lattice_unita(1,1) 
     >                              + n2*lattice_unita(1,2) 
     >                              + n3*lattice_unita(1,3)
            rxyz(2) = ion_rion(2,i) + n1*lattice_unita(2,1) 
     >                              + n2*lattice_unita(2,2) 
     >                              + n3*lattice_unita(2,3)
            rxyz(3) = ion_rion(3,i) + n1*lattice_unita(3,1) 
     >                              + n2*lattice_unita(3,2) 
     >                              + n3*lattice_unita(3,3)
            if (taskid.eq.MASTER) then
               WRITE(unit_id,'(I5,4F12.6)') int(q),q,rxyz
            end if
         end do
         end do
         end do

c         if (taskid.eq.MASTER) then
c           WRITE(unit_id,'(I5,4F12.6)') int(q),q,rxyz
c           call special_position_tolerance(position_tolerance,
c     >                                   unit_id,q,rxyz)
c         end if
      end do

*     **** write orbital header ****
      if (number.gt.0) then
         if (taskid.eq.MASTER) write(unit_id,*) 1,number 
      end if

*     **** allocate space ****
      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** transpose grid ****
      if (control_mapping().eq.1) 
     > call D3dB_c_transpose_jk(1,rho,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))

*     **** write grid ****
      call D3dB_r_FormatWrite_reverse(1,unit_id,rho,dcpl_mb(tmp1(1)))

*     **** deallocate space ****
      value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
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

*     ***********************************************
*     *                                             *
*     *          special_position_tolerance         *
*     *                                             *
*     ***********************************************

      subroutine special_position_tolerance(position_tolerance,
     >                                      unit,q,r2)
      implicit none
      real*8  position_tolerance
      integer unit
      real*8  q
      real*8  r2(3)

*     **** Local variables defined ****
      real*8  fa1,fa2,fa3
      real*8  a(3,3),b(3,3),volume
      integer i,j
      real*8 rxyz(3)

*      **** external functions ****
       real*8   lattice_unita
       external lattice_unita

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

*      *** Break the Ion positions into the a1, a2, and a3 components ***
       fa1 =  b(1,1) * r2(1)
     >     +  b(2,1) * r2(2)
     >     +  b(3,1) * r2(3)

       fa2 =  b(1,2) * r2(1)
     >     +  b(2,2) * r2(2)
     >     +  b(3,2) * r2(3)

       fa3 =  b(1,3) * r2(1)
     >     +  b(2,3) * r2(2)
     >     +  b(3,3) * r2(3)

 
       if ((fa1+position_tolerance) .GT. (0.5d0)) THEN
         rxyz(1) = r2(1) - lattice_unita(1,1)
         rxyz(2) = r2(2) - lattice_unita(2,1)
         rxyz(3) = r2(3) - lattice_unita(3,1)
         WRITE(unit,'(I5,4F12.6)') int(q),q,rxyz
       end if

       if ((fa1-position_tolerance) .LT. (-0.5d0)) THEN
         rxyz(1) = r2(1) + lattice_unita(1,1)
         rxyz(2) = r2(2) + lattice_unita(2,1)
         rxyz(3) = r2(3) + lattice_unita(3,1)
         WRITE(unit,'(I5,4F12.6)') int(q),q,rxyz
       end if

       if ((fa2+position_tolerance) .GT. (0.5d0)) THEN
         rxyz(1) = r2(1) - lattice_unita(1,2)
         rxyz(2) = r2(2) - lattice_unita(2,2)
         rxyz(3) = r2(3) - lattice_unita(3,2)
         WRITE(unit,'(I5,4F12.6)') int(q),q,rxyz
       end if

       if ((fa2-position_tolerance) .LT. (-0.5d0)) THEN
         rxyz(1) = r2(1) + lattice_unita(1,2)
         rxyz(2) = r2(2) + lattice_unita(2,2)
         rxyz(3) = r2(3) + lattice_unita(3,2)
         WRITE(unit,'(I5,4F12.6)') int(q),q,rxyz
       end if


       if ((fa3+position_tolerance) .GT. (0.5d0)) THEN
         rxyz(1) = r2(1) - lattice_unita(1,3)
         rxyz(2) = r2(2) - lattice_unita(2,3)
         rxyz(3) = r2(3) - lattice_unita(3,3)
         WRITE(unit,'(I5,4F12.6)') int(q),q,rxyz
       end if

       if ((fa3-position_tolerance) .LT. (-0.5d0)) THEN
         rxyz(1) = r2(1) + lattice_unita(1,3)
         rxyz(2) = r2(2) + lattice_unita(2,3)
         rxyz(3) = r2(3) + lattice_unita(3,3)
         WRITE(unit,'(I5,4F12.6)') int(q),q,rxyz
       end if

      return
      end

*     ***********************************************
*     *                                             *
*     *          special_position_count             *
*     *                                             *
*     ***********************************************

      subroutine special_position_count(position_tolerance,
     >                                  r2,count)
      implicit none
      real*8  position_tolerance
      real*8  r2(3)
      integer count

*     **** Local variables defined ****
      real*8  fa1,fa2,fa3
      real*8  a(3,3),b(3,3),volume
      integer i,j

*      **** external functions ****
       real*8   lattice_unita
       external lattice_unita

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

*      *** Break the Ion positions into the a1, a2, and a3 components ***
       fa1 =  b(1,1) * r2(1)
     >     +  b(2,1) * r2(2)
     >     +  b(3,1) * r2(3)

       fa2 =  b(1,2) * r2(1)
     >     +  b(2,2) * r2(2)
     >     +  b(3,2) * r2(3)

       fa3 =  b(1,3) * r2(1)
     >     +  b(2,3) * r2(2)
     >     +  b(3,3) * r2(3)

 
       if ((fa1+position_tolerance) .GT. (0.5d0)) THEN
         count = count+1
       end if

       if ((fa1-position_tolerance) .LT. (-0.5d0)) THEN
         count = count + 1
       end if

       if ((fa2+position_tolerance) .GT. (0.5d0)) THEN
         count = count + 1
       end if

       if ((fa2-position_tolerance) .LT. (-0.5d0)) THEN
         count = count + 1
       end if


       if ((fa3+position_tolerance) .GT. (0.5d0)) THEN
         count = count + 1
       end if

       if ((fa3-position_tolerance) .LT. (-0.5d0)) THEN
         count = count + 1
       end if

      return
      end


*     ********************************
*     *                				 *
*     *        pspw_add_core_dng     *
*     *                 			 *
*     ********************************
      subroutine pspw_add_core_dng(rcut,dng)
      implicit none
#include "errquit.fh"
      real*8     rcut
      complex*16 dng(*)

#include "bafdecls.fh"

*     *** local variables ***
      logical value
      integer nfft3d
      integer i,ii
      integer exi(2),vg(2),G(3)
      real*8  w,scal,gg

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      real*8   ion_zv,lattice_omega
      external G_indx,ion_nion,ion_katm
      external ion_zv,lattice_omega


      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'vg',vg(2),vg(1))
      if (.not. value) call errquit('pspw_add_core_dng:push stack',0,
     &       MA_ERR)

      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)
      w    = 0.25d0*rcut*rcut


      do ii=1,ion_nion()
         scal = -ion_zv(ii)/lattice_omega()

*        *** fourier transform of a gaussian ***
         do i=1,nfft3d
            gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >            + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >            + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

            dbl_mb(vg(1)+i-1) = (scal)*exp(-w*gg)
         end do

*        **** add to dng ***
         call strfac(ii,dcpl_mb(exi(1)))
         call Pack_c_pack(0,dcpl_mb(exi(1)))
         call Pack_t_pack(0,dbl_mb(vg(1)))
c         call Pack_tc_Mul(0,dbl_mb(vg(1)),
c     >                   dcpl_mb(exi(1)),
c     >                   dcpl_mb(exi(1)))
c         call Pack_cc_Sum(0,dng,dcpl_mb(exi(1)),dng)
         call Pack_tc_Mul2(0,dbl_mb(vg(1)),dcpl_mb(exi(1)))
         call Pack_cc_Sum2(0,dcpl_mb(exi(1)),dng)

      end do

      value = BA_pop_stack(vg(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) call errquit('pspw_add_core_dng:pop stack',1,
     &       MA_ERR)
      return
      end
      
*     ********************************
*     *                              *
*     *        pspw_add_core_pot     *
*     *                              *
*     ********************************
      subroutine pspw_add_core_pot(rcut,vh)
      implicit none
#include "errquit.fh"
      real*8  rcut
      real*8 vh(*)

#include "bafdecls.fh"

*     *** local variables ***
      logical value
      integer n2ft3d
      integer i,ii
      integer rgrid(2)
      real*8  c,q,r,xerf,yerf,sqrt_pi
      real*8  xii,yii,zii
      real*8  xi,yi,zi

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   ion_zv,ion_rion,util_erf
      external ion_nion,ion_katm
      external ion_zv,ion_rion,util_erf


      call D3dB_n2ft3d(1,n2ft3d)
      value = BA_push_get(mt_dbl,3*n2ft3d,'rgrid',rgrid(2),rgrid(1))
      if (.not. value) call errquit('pspw_add_core_pot:push stack',0,
     &       MA_ERR)

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      c       = 1.0d0/rcut

      call lattice_r_grid(dbl_mb(rgrid(1)))

      do ii=1,ion_nion()
         q = -ion_zv(ii)
         xii = ion_rion(1,ii)
         yii = ion_rion(2,ii)
         zii = ion_rion(3,ii)

         do i=1,n2ft3d
            xi = dbl_mb(rgrid(1) + 3*(i-1))
            yi = dbl_mb(rgrid(1) + 3*(i-1)+1)
            zi = dbl_mb(rgrid(1) + 3*(i-1)+2)
            r = dsqrt((xii-xi)**2 + (yii-yi)**2 + (zii-zi)**2)

            if (r .gt. 1.0d-15) then
              xerf  = r*c
              yerf  = util_erf(xerf)
              vh(i) = vh(i) + (q/r)*yerf
            else
              vh(i) = vh(i) + 2.0d0*q*c/sqrt_pi
            end if
         end do

      end do

      value = BA_pop_stack(rgrid(2))
      if (.not. value) call errquit('pspw_add_core_pot:pop stack',1,
     &       MA_ERR)
      return
      end
      

*     ********************************
*     *                              *
*     *        psi_translate         *
*     *                              *
*     ********************************
      subroutine psi_translate(trans,npack1,nemax,psi)
      implicit none
      real*8 trans(3)
      integer npack1,nemax
      complex*16 psi(npack1,nemax)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer n,nion
      integer exi(2),rion(2)

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unita
      external ion_nion
      external lattice_unita

      nion  = ion_nion()
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nion,'rion', rion(2), rion(1))
      if (.not. value) call errquit('psi_translate:push stack',0,0)
      call dcopy(3*nion,0.0d0,0,dbl_mb(rion(1)),1)

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
      call phafac_rion(dbl_mb(rion(1)))
      call strfac_pack(1,1,dcpl_mb(exi(1)))
      
 
      do n=1,nemax
c        call Pack_cc_Mul(1,dcpl_mb(exi(1)),psi(1,n),psi(1,n))
        call Pack_cc_Mul2(1,dcpl_mb(exi(1)),psi(1,n))
        !call Pack_cc_conjgMul(1,dcpl_mb(exi(1)),psi(1,n),psi(1,n))
      end do


      value =           BA_pop_stack(rion(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) call errquit('psi_translate:pop stack',1,0)
      return
      end



*
*     ***********************************************
*     *                                             *
*     *          dplot_gcube_write3d                *
*     *                                             *
*     ***********************************************

      subroutine dplot_gcube_write3d(rtdb,filename,
     >                             number,cube_comment,rho)
      implicit none
      integer rtdb
      character*50 filename
      integer      number
      character*72 cube_comment
      real*8 rho(*)
     

#include "bafdecls.fh"
#include "btdb.fh"
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

      integer  G_indx
      external G_indx


      call Parallel_taskid(taskid)

      
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

*     **** read the origin ***
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:o',
     >               mt_dbl,3,r0)) then
        r0(1) = 0.0
        r0(2) = 0.0
        r0(3) = 0.0
      end if

*     **** read the x-axis,y-axis,z-axis ****
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:x',
     >               mt_dbl,3,xaxis)) then
        xaxis(1) = 1.0
        xaxis(2) = 0.0
        xaxis(3) = 0.0
      end if
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:y',
     >               mt_dbl,3,yaxis)) then
        yaxis(1) = 0.0
        yaxis(2) = 1.0
        yaxis(3) = 0.0
      end if
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:z',
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
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:nx',
     >               mt_int,1,nx))
     >  nx = 32
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:ny',
     >               mt_int,1,ny))
     >  ny = 32
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:nz',
     >               mt_int,1,nz))
     >  nz = 32
       nxh = nx/2
       nyh = ny/2
       nzh = nz/2

*     **** read sizex,sizey,sizez ****
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:sizex',
     >               mt_dbl,2,sizex)) then
        sizex(1) = 0.0
        sizex(2) = 10.0
      end if
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:sizey',
     >               mt_dbl,2,sizey)) then
        sizey(1) = 0.0
        sizey(2) = 10.0
      end if
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:3d_grid:sizez',
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


      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

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
     >      btdb_get(rtdb,'pspw_dplot:position_tolerance',
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
      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nx*ny*nz),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*     **** fft density  ****
      scal1=1.0d0/dble(np1*np2*np3)
      call D3dB_r_SMul(1,scal1,rho,dcpl_mb(tmp2(1)))
      call D3dB_r_Zero_Ends(1,dcpl_mb(tmp2(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(tmp2(1)))
      call Pack_c_pack(0,dcpl_mb(tmp2(1)))


*     **** density at mesh points ****
      Gx= G_indx(1)
      Gy= G_indx(2)
      Gz= G_indx(3)
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
        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call Pack_cc_dot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)),w)
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
      value = BA_pop_stack(tmp3(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
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





*
*     ***********************************************
*     *                                             *
*     *          dplot_gcube_write1d                *
*     *                                             *
*     ***********************************************

      subroutine dplot_gcube_write1d(rtdb,filename,
     >                             number,cube_comment,rho)
      implicit none
      integer rtdb
      character*50 filename
      integer      number
      character*72 cube_comment
      real*8 rho(*)
     

#include "bafdecls.fh"
#include "btdb.fh"
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

      integer  G_indx
      external G_indx


      call Parallel_taskid(taskid)

      
*     **** OPEN cube FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         OPEN(unit_id,file=full_filename(1:l),form='formatted')
c         WRITE(unit_id,*) 'molecule'
c         WRITE(unit_id,*)  cube_comment
      end if

*     **** open geom ****
      nion = ion_nion()
c     value = geom_create(geom1,'geometry')
c     value = value.and.geom_rtdb_load(rtdb,geom1,'geometry')
c     value = value.and.geom_ncent(geom1,nion)
c     if (.not. value) call errquit('opening geometry',0)

*     **** read the origin ***
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:1d_grid:o',
     >               mt_dbl,3,r0)) then
        r0(1) = 0.0
        r0(2) = 0.0
        r0(3) = 0.0
      end if

*     **** read the x point ****
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:1d_grid:x',
     >               mt_dbl,3,xaxis)) then
        xaxis(1) = 1.0
        xaxis(2) = 0.0
        xaxis(3) = 0.0
      end if

*     **** read nx ***
      if (.not.
     >      btdb_get(rtdb,'pspw_dplot:1d_grid:nx',
     >               mt_int,1,nx))
     >  nx = 32
       nxh = nx/2

      if (taskid.eq.MASTER) then
         write(*,*)
         write(*,1240) 
         write(*,1241) r0
         write(*,1242) xaxis,nx
         write(*,*)
         write(*,*)
      end if
 1240 FORMAT(5x,'    1d-Grid Generation')
 1241 FORMAT(5x,'    origin=<',3f8.3,' >')
 1242 FORMAT(5x,'    x-point=<',3f8.3,' >  nx=',I3)

*     *** define line ****
      do i=1,3
        xaxis(i) = xaxis(i)-r0(i)
      end do

*     ** normalization of x **
      sx = dsqrt(xaxis(1)**2+xaxis(2)**2+xaxis(3)**2)
      do i=1,3
        xaxis(i) = xaxis(i) / sx
      end do


      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      dx = sx/dble(nx)
      do i=1,3
         ua(i,1) = dx*xaxis(i)
         ua(i,2) = 0.0d0
         ua(i,3) = 0.0d0
      end do

*     ** shift origin **
      qq(1) = r0(1)  
     >       - 0.5d0*(lattice_unita(1,1)
     >               +lattice_unita(1,2)
     >               +lattice_unita(1,3))
      qq(2) = r0(2)  
     >       - 0.5d0*(lattice_unita(2,1)
     >               +lattice_unita(2,2)
     >               +lattice_unita(2,3))
      qq(3) = r0(3)  
     >       - 0.5d0*(lattice_unita(3,1)
     >               +lattice_unita(3,2)
     >               +lattice_unita(3,3))


*     **** allocate space ****
      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*     **** fft density  ****
      scal1=1.0d0/dble(np1*np2*np3)
      call D3dB_r_SMul(1,scal1,rho,dcpl_mb(tmp2(1)))
      call D3dB_r_Zero_Ends(1,dcpl_mb(tmp2(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(tmp2(1)))
      call Pack_c_pack(0,dcpl_mb(tmp2(1)))


*     **** density at mesh points ****
      Gx= G_indx(1)
      Gy= G_indx(2)
      Gz= G_indx(3)
      do i=1,nx
        !ii = i-1-nxh
        ii = i-1
        p(1) = ii*ua(1,1) + qq(1)
        p(2) = ii*ua(2,1) + qq(2)
        p(3) = ii*ua(3,1) + qq(3)
        do l=1,nfft3d
          w = dbl_mb(Gx+l-1)*p(1) 
     >      + dbl_mb(Gy+l-1)*p(2) 
     >      + dbl_mb(Gz+l-1)*p(3)
          dcpl_mb(tmp1(1)+l-1)=dcmplx(dcos(w),-dsin(w))
        end do
        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call Pack_cc_dot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)),w)

        if (taskid.eq.MASTER)
     >     write(unit_id,'(2e13.5)') (i-1)*dx,w
      end do


*     **** deallocate space ****
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)



*     **** CLOSE FILE ****
      if (taskid.eq.MASTER) then
         CLOSE(unit_id)     
      end if

      return
      end 


