      logical function pspw_dplot(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"


      logical value
      integer taskid,np
      integer MASTER
      parameter (MASTER=0)

*     **** local variables ***

      integer count
      integer ngrid(3)
      integer nfft3d,n2ft3d,npack1
      integer ispin,ne(2),nemax
      integer psi2(2),dn(2),psir(2)
      real*8  cpu1,cpu2
      real*8  position_tolerance

*     **** external functions ****
      logical      control_read
      integer      pack_nwave_all
      integer      control_ngrid,pack_nwave,control_version
      real*8       lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8       lattice_unitg
      character*12 control_boundry

      external     control_read
      external     pack_nwave_all
      external     control_ngrid,pack_nwave,control_version
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
         call message(1)
         write(6,1110)
         call util_flush(6)
      end if

*     **** get position_tolerance ****
      if (.not.
     >      rtdb_get(rtdb,'pspw_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0


*     **** read control file ****
      value = control_read(4,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)


      if ((np.gt.(control_ngrid(3)))) then
        if (taskid.eq.MASTER)
     >     write(6,*) 'Error: the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   control_ngrid(3),')'
         call util_flush(6)
         call errquit('PSPW:Bad number of processors:',np)
      end if

*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3))
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)

*     ***** allocate psi2 wavefunctions ****
      call psi_get_ne(ispin,ne)
      nemax = ne(1)+ne(2)
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      if (.not. value) call errquit('out of heap memory',0)


*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)))

*     **** allocate other variables *****
      value = value.and.
     >        MA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir',psir(2),psir(1))

      if (.not. value) call errquit('pspw_dplot:out of heap memory',0)


*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
         write(6,1111) np
         write(6,1115)
         write(6,1121) control_boundry(),control_version()
         write(6,1220) ne(1),ne(ispin),' ( fourier space)'
         write(6,1225) position_tolerance
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

      call dplot_loop(rtdb,
     >                ispin,ne,
     >                npack1,nfft3d,nemax,
     >                dcpl_mb(psi2(1)),
     >                dbl_mb(dn(1)),
     >                dcpl_mb(psir(1)))

      value = value.and.rtdb_get(rtdb,'pspw_dplot:count',
     >                            mt_int,1,count)


*     **** deallocate heap memory ****
      call mask_end()
      call Pack_end()
      call G_end()

      value = MA_free_heap(psir(2))
      value = value.and.MA_free_heap(dn(2))
      value = value.and.MA_free_heap(psi2(2))
      call D3dB_end(1)
      if (.not. value) call errquit('pspw_dplot:freeing heap memory',0)


      if (taskid.eq.MASTER) then
         call current_second(cpu2)
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'total       : ',(cpu2-cpu1)
         write(6,*) 
         call message(4)
      end if

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
 1115 FORMAT(/' options:')
 1121 FORMAT(5X,' boundary conditions   = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
 1225 FORMAT(/' position tolerance = ',E12.6)
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
 

  
      subroutine dplot_loop(rtdb,
     >                      ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi,
     >                      dn,
     >                      psi_r)
      implicit none
      integer    rtdb
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      real*8     psi_r(2*nfft3d,nemax)

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** local variables ****
      logical value
      integer taskid
      integer MASTER
      parameter (MASTER=0)

      integer n2ft3d
      integer i,ms,n,count,number,ia
      integer n1(2),n2(2)
      integer nx,ny,nz
      integer rho(2)
      real*8  scal1,scal2

      character*72 cube_comment
      character*30 name1,name2
      character*30 filename
      integer name1_len,name2_len,ind


*     **** external functions ****
      real*8 lattice_omega
      external lattice_omega

      
      call Parallel_taskid(taskid)
      n2ft3d = 2*nfft3d
      value = MA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      if (.not. value) call errquit('out of stack memory',0)

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
         call D3dB_r_SMul(1,dsqrt(scal2),psi_r(1,n),psi_r(1,n))
      end do


*     *********************
*     **** generate dn ****
*     *********************
      call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            do i=1,n2ft3d
c              dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
               dn(i,ms) = dn(i,ms) + (psi_r(i,n)**2)
            end do
         end do
         call D3dB_r_Zero_Ends(1,dn(1,ms))
      end do

*     ********************************************
*     **** loop over orbital and density list ****
*     ********************************************
      value = value.and.rtdb_get(rtdb,'pspw_dplot:count',
     >                            mt_int,1,count)
      do i=1,count

*       **** define name  - not very elegent and could break if ****
*       ****                      count becomes very large      ****
        ia = ICHAR('a')
        name1 = 'pspw_dplot:filename'//CHAR(i-1+ia)
        name2 = 'pspw_dplot:number'//CHAR(i-1+ia)
        name1_len = index(name1,' ') - 1
        name2_len = index(name2,' ') - 1

        value = rtdb_cget(rtdb,name1(1:name1_len),1,filename)
        ind = index(filename,' ') - 1
        value = value.and.rtdb_get(rtdb,name2(1:name2_len),
     >                           mt_int,1,number)
        if (.not.value)
     >     call errquit(
     >     'pspw_dplot: rtdb_get failed for orbital', 0)


*       **** outputing density ****
        if (number.lt.0) then
           number = -number

          goto ( 710, 720, 730, 740 ) number
          call errquit(
     >      'dplot_loop: unimplemented directive', number)
       
*          *************
*          *** total ***
*          *************
 710       call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
           if (taskid.eq.MASTER) then
             write(*,*) '   writing total density',
     >                  ' to filename: ',filename(1:ind)
           end if
           cube_comment = "SCF Total Density"
           goto 750

*          ******************
*          *** difference ***
*          ******************
 720       call D3dB_rr_Sub(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
           if (taskid.eq.MASTER) then
             write(*,*) '   writing difference density',
     >                  ' to filename: ',filename(1:ind)
           end if
           cube_comment = "SCF Spin Density"
           goto 750

*          *************
*          *** alpha ***
*          *************
 730       call D3dB_r_Copy(1,dn(1,1),dbl_mb(rho(1)))
           if (taskid.eq.MASTER) then
             write(*,*) '   writing alpha density',
     >                  ' to filename: ',filename(1:ind)
           end if
           cube_comment = "SCF Alpha Density"
           goto 750

*          ************
*          *** beta ***
*          ************
 740       call D3dB_r_Copy(1,dn(1,ispin),dbl_mb(rho(1)))
           if (taskid.eq.MASTER) then
             write(*,*) '   writing beta density',
     >                  ' to filename: ',filename(1:ind)
           end if
           cube_comment = "SCF Beta Density"
           goto 750
 750       continue
           number = -number

*       **** outputing wavefunction ****
        else
           if ((number.gt.nemax).or.(number.lt.1)) then
              if (taskid.eq.MASTER) then
                write(*,*)  '   Bad orbital number ', number,
     >                     ', changing to orbital number 1.'
              end if
              number = 1
           end if
          
           call D3dB_r_Copy(1,psi_r(1,number),dbl_mb(rho(1)))
           if (taskid.eq.MASTER) then
             write(*,*) '   writing orbital ',number, 
     >                  ' to filename: ',filename(1:ind)
           end if
          cube_comment = 'SCF Molecular Orbitals'
        end if

        call dplot_gcube_write(rtdb,filename,
     >                         number,cube_comment,dbl_mb(rho(1)))

      end do

*     **** dealocate MA local variables ****
      value = MA_pop_stack(rho(2))
      if (.not. value) call errquit('popping of stack memory',0)

    
      return
      end

      subroutine dplot_gcube_write(rtdb,filename,
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

*     **** local variables ****
      logical value
      integer unit_id
      parameter (unit_id = 72)

      integer taskid
      integer MASTER
      parameter (MASTER=0)


*     **** geometry variables ****
      integer geom1
      character*16 t
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
      real*8   lattice_unita
      external lattice_unita

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
      value = geom_create(geom1,'geometry')
      value = value.and.geom_rtdb_load(rtdb,geom1,'geometry')
      value = value.and.geom_ncent(geom1,nion)
      if (.not. value) call errquit('opening geometry',0)

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
     >      rtdb_get(rtdb,'pspw_dplot:position_tolerance',
     >               mt_dbl,1,position_tolerance))
     >      position_tolerance = 0.0d0

      nion2 = nion
      do i=1,nion
         value = geom_cent_get(geom1,i,t,rxyz,q)
         if (.not. value) call errquit('reading geometry',0)
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
         value = geom_cent_get(geom1,i,t,rxyz,q)
         if (.not. value) call errquit('reading geometry',0)
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
      value = MA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** transpose grid ****
      call D3dB_c_transpose_jk(1,rho,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))

*     **** write grid ****
      call D3dB_r_FormatWrite_reverse(1,unit_id,rho,dcpl_mb(tmp1(1)))

*     **** deallocate space ****
      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0)

   

*     **** close geom ****
      value = geom_destroy(geom1)
      if (.not. value) call errquit('closing geometry',0)



*     **** CLOSE cube FILE ****
      if (taskid.eq.MASTER) then
         CLOSE(unit_id)
      end if

      return
      end 


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
