*
* $Id: pspw_wannier.F,v 1.16 2007-12-17 21:46:01 d3p708 Exp $
*

      logical function pspw_wannier(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"


*     **** local variables ****
      integer MASTER,taskid,tmp_len
      parameter (MASTER=0,tmp_len=140)

      logical value,oprint
      integer mapping,mapping1d
      integer i,j,ms,n,ne(2),ispin,ngrid(3),nemax
      integer neq(2),nemaxq,shift,rank
      integer nfft3d,n2ft3d,npack1
      integer info
      integer np

      integer psi(2),psi2(2),psi_r(2)
      integer X(2,6)
      integer A(2)
      real*8 bv(3,6),wrk(6,6),wts(6),bmat(3,3)
      real*8 b(3),ixmat(3,6)
      real*8 xx,yy,zz,tmp1(tmp_len)
      real*8 dx,dy,dz
      complex*16 arg,wx
      character*2 labels(6)

c      double precision unita(3,3),b(3),alpha

c      complex*16 wx,wy,wz

c      character*30 new_wavefunction_filename
c      character*30 old_wavefunction_filename
c      character*255 full_filename,full_filename2

*     **** external functions ****
      logical      control_read
      integer      control_ngrid,control_ispin,control_mapping
      integer      control_mapping1d
      integer      pack_nwave,pack_nwave_all
      integer      control_np_orbital
      real*8       lattice_unita,lattice_unitg
      real*8       lattice_ecut,lattice_wcut,lattice_omega
      character*30 control_input_psi,control_output_psi
      external     control_read
      external     control_ngrid,control_ispin,control_mapping
      external     control_mapping1d
      external     pack_nwave,pack_nwave_all
      external     control_np_orbital
      external     lattice_unita,lattice_unitg
      external     lattice_ecut,lattice_wcut,lattice_omega
      external     control_input_psi,control_output_psi

      value = .true.

*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      oprint = (taskid.eq.MASTER)

*     ***** print out header ****
      if (oprint) then
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
         write(6,1041)
         write(6,1042)
         write(6,1043)
         write(6,1010)
         write(6,1000)
         call nwpw_message(1)
      end if
      
      value = control_read(9,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      mapping = control_mapping()

      call Parallel2d_Init(control_np_orbital())



      !**** Make sure ngrid is consistent with mapping ***
      if (mapping.eq.1) then
        if ((np.gt.ngrid(3)).or.(ngrid(2).ne.ngrid(3))) then
          if (taskid.eq.MASTER) then
            write(6,*) 'Error: the number of processors must ',
     >                 ' be in the range ( 1 ...ngrid(3)=',
     >                   ngrid(3),')'
           write(6,*) ' and ngrid(2) == ngrid(3), ',
     >                ' ngrid(2)=',ngrid(2),
     >                ' ngrid(3)=',ngrid(3)
          end if
          value = .false.
          go to 9000
        end if
        if (mod(ngrid(1),2).ne.0) then
          if (taskid.eq.MASTER) then
           write(6,*)
     >      'Error: ngrid(1) must be even (ngrid(1) mod 2 == 0)'
           write(6,*) 'Error: ngrid(1)=',ngrid(1)
          end if
          value = .false.
          go to 9000
        end if
      end if
      if (mapping.eq.2) then
        if (np.gt.(ngrid(2)*ngrid(3))) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' (ngrid(1)/2+1)*ngrid(2),',
     >                ' (ngrid(1)/2+1)*ngrid(3))'
           write(6,*) 'Error: np > ngrid(2)*ngrid(3)'
           write(6,*) 'Error: the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   ngrid(2)*ngrid(3),')'
          end if
          value = .false.
          go to 9000
        end if
        if (np.gt.((ngrid(1)/2+1)*ngrid(2))) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' (ngrid(1)/2+1)*ngrid(2),',
     >                ' (ngrid(1)/2+1)*ngrid(3))'
           write(6,*) 'Error: np > (ngrid(1)/2+1)*ngrid(2)'
           write(6,*) 'Error: the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   (ngrid(1)/2+1)*ngrid(2),')'
          end if
          value = .false.
          go to 9000
        end if
        if (np.gt.((ngrid(1)/2+1)*ngrid(3))) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' (ngrid(1)/2+1)*ngrid(2),',
     >                ' (ngrid(1)/2+1)*ngrid(3))'
           write(6,*) 'Error: np > (ngrid(1)/2+1)*ngrid(3)'
           write(6,*) 'Error: the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   (ngrid(1)/2+1)*ngrid(3),')'
          end if
          value = .false.
          go to 9000
        end if
        if (mod(ngrid(1),2).ne.0) then
          if (taskid.eq.MASTER) then
           write(6,*)
     >      'Error: ngrid(1) must be even (ngrid(1) mod 2 == 0)'
           write(6,*) 'Error: ngrid(1)=',ngrid(1)
          end if
          value = .false.
          go to 9000
        end if

      end if


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


*     ***** allocate psi and psi_r wavefunctions ****
      call psi_get_ne(ispin,ne)
      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)
      nemaxq = neq(1)+neq(2)

      
      value = MA_alloc_get(mt_dcpl,npack1*nemaxq,
     >                     'psi',psi(2),psi(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*nemaxq,
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d*nemaxq,
     >                     'psi_r',psi_r(2),psi_r(1))
      if (.not. value) 
     >   call errquit('pspw_wannier:out of heap memory',0,0)


*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi(1)))


      if (oprint) then
  
         write(6,109) control_input_psi()
         write(6,110) control_output_psi()

         write(6,1220) ne(1),ne(ispin),' ( fourier space)'
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

      end if


*     **** psi_r <-- FFT(psi) - complex to real FFT ****
      nemax = ne(1)+ne(2)
      do n=1,nemax
         call Pack_c_Copy(1,dcpl_mb(psi(1)  +(n-1)*npack1),
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      end do
      call Grsm_gh_fftb(nfft3d,nemax,dbl_mb(psi_r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemax,dbl_mb(psi_r(1)))


*     ***** allocate X,Y,Z  ****
      labels(1)='X1'
      labels(2)='X2'
      labels(3)='X3'
      labels(4)='X4'
      labels(5)='X5'
      labels(6)='X6'
      n=ne(1)
      if (n.lt.ne(2)) n=ne(2)
      value=.true.
      do j=1,6
      value = value.and.MA_alloc_get(mt_dcpl,
     >                     n*n,
     >                     labels(j),X(2,j),X(1,j))      
      end do
      value = value.and.
     >        MA_alloc_get(mt_dcpl,n*n,
     >                     'A',A(2),A(1))
      if (.not. value) 
     >   call errquit('pspw_wannier:out of heap memory',1,0)

      ixmat(1,1)=1.0d0
      ixmat(2,1)=0.0d0
      ixmat(3,1)=0.0d0

      ixmat(1,2)=1.0d0
      ixmat(2,2)=1.0d0
      ixmat(3,2)=0.0d0

      ixmat(1,3)=1.0d0
      ixmat(2,3)=0.0d0
      ixmat(3,3)=1.0d0

      ixmat(1,4)=0.0d0
      ixmat(2,4)=1.0d0
      ixmat(3,4)=0.0d0

      ixmat(1,5)=0.0d0
      ixmat(2,5)=1.0d0
      ixmat(3,5)=1.0d0

      ixmat(1,6)=0.0d0
      ixmat(2,6)=0.0d0
      ixmat(3,6)=1.0d0

      do i=1,3
         bmat(i,1)=lattice_unitg(1,i)
         bmat(i,2)=lattice_unitg(2,i)
         bmat(i,3)=lattice_unitg(3,i)
      end do

      do i=1,6
         xx=0.0d0
         yy=0.0d0
         zz=0.0d0
         do j=1,3
           xx=xx+bmat(j,1)*ixmat(j,i)
           yy=yy+bmat(j,2)*ixmat(j,i)
           zz=zz+bmat(j,3)*ixmat(j,i)
         end do
         bv(1,i)=xx
         bv(2,i)=yy
         bv(3,i)=zz
      end do

      do i=1,6
         wrk(1,i)=bv(1,i)*bv(1,i)
         wrk(2,i)=bv(1,i)*bv(2,i)
         wrk(3,i)=bv(1,i)*bv(3,i)
         wrk(4,i)=bv(2,i)*bv(2,i)
         wrk(5,i)=bv(2,i)*bv(3,i)
         wrk(6,i)=bv(3,i)*bv(3,i)
         wts(i)=0.0d0
      end do
      wts(1)=1.0d0
      wts(4)=1.0d0
      wts(6)=1.0d0
      call dgels('N',6,6,1,wrk,6,wts,6,tmp1,tmp_len,info)
      if (info.ne.0) then
        write(*,*)"Illegal argument in call to dgels"
        call flush(6)
      end if
      rank=0
      do i=1,6
         if (dabs(wts(i)).gt.1.e-6) then
           rank=rank+1
           wrk(1,rank)=bv(1,i)
           wrk(2,rank)=bv(2,i)
           wrk(3,rank)=bv(3,i)
           wrk(4,rank)=wts(i)
         end if
      end do
      do i=1,rank
         bv(1,i)=wrk(1,i)
         bv(2,i)=wrk(2,i)
         bv(3,i)=wrk(3,i)
         wts(i)=wrk(4,i)
      end do

      dx=0.0d0
      dy=0.0d0
      dz=0.0d0
      
      do ms=1,ispin

        if ((oprint).and.(ms.eq.1)) write(6,1300)
        if ((oprint).and.(ms.eq.2)) write(6,1301)


        do i=1,rank
          b(1) = bv(1,i)
          b(2) = bv(2,i)
          b(3) = bv(3,i)
          call silvestrelli_overlap(
     >                  b,
     >                  ne(ms),
     >                  dbl_mb(psi_r(1)+(ms-1)*n2ft3d*ne(1)),
     >                  dcpl_mb(X(1,i)))
        end do

        if (oprint) then
          write(6,1380)
          do i=1,ne(ms)
            shift=(i-1)*ne(ms)+(i-1)
            xx=0.0d0
            yy=0.0d0
            zz=0.0d0
            do j=1,rank
               arg=dcpl_mb(X(1,j)+shift)
               arg= -wts(j)*datan2(dimag(arg),dble(arg))
               xx=xx+bv(1,j)*arg
               yy=yy+bv(2,j)*arg
               zz=zz+bv(3,j)*arg
            end do

            if (ms.eq.1)
     >        write(6,1390) i,xx,yy,zz


            if (ms.eq.2)
     >        write(6,1391) i+ne(1),xx,yy,zz
          end do
        end if

        call silvestrelli_minimize(ne(ms),rank,wts,
     >                             dcpl_mb(X(1,1)),dcpl_mb(X(1,2)),
     >                             dcpl_mb(X(1,3)),dcpl_mb(X(1,4)),
     >                             dcpl_mb(X(1,5)),dcpl_mb(X(1,6)),
     >                             dcpl_mb(A(1)))

        if (oprint) then
          write(6,1381)
          do i=1,ne(ms)
            shift=i-1
            shift=shift*ne(ms)+shift
            xx=0.0d0
            yy=0.0d0
            zz=0.0d0
            do j=1,rank
               arg=dcpl_mb(X(1,j)+shift)
               arg= -wts(j)*datan2(dimag(arg),dble(arg))
               xx=xx+bv(1,j)*arg
               yy=yy+bv(2,j)*arg
               zz=zz+bv(3,j)*arg
            end do
            dx=dx+xx
            dy=dy+yy
            dz=dz+zz
            if (ms.eq.1)
     >        write(6,1390) i,xx,yy,zz


            if (ms.eq.2)
     >        write(6,1391) i+ne(1),xx,yy,zz
          end do
        end if




*       **** rotate psi to produce Wannier orbitals ****
        do i=1,ne(ms)
          call Pack_c_Zero(1,
     >                     dcpl_mb(psi2(1)+(i-1) *npack1
     >                                    +(ms-1)*npack1*ne(1)))
          do j=1,ne(ms)
              wx    = dcpl_mb(A(1)+(i-1)+(j-1)*ne(ms) 
     >                            +(ms-1)*ne(1)*ne(1))
c              alpha = dble(wx)

              call Pack_cc_zaxpy(1,wx,
     >                           dcpl_mb(psi(1) +(j-1)*npack1
     >                                          +(ms-1)*npack1*ne(1)),
     >                           dcpl_mb(psi2(1)+(i-1)*npack1
     >                                          +(ms-1)*npack1*ne(1)))
          end do
        end do

      end do !* ms *
cccccccccccccccccccccccccccccccccccccccccccccc
c  Molecular dipoles from Resta's theory!
ccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,1776)
      WRITE(6,1777)dx,dy,dz

*     ***** write psi2 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)))


*     **** deallocate heap memory ****
      call mask_end()
      call Pack_end()
      call G_end()
      call D3dB_pfft_end()
      call D3dB_end(1)
      call Dne_end()

      value =  MA_free_heap(A(2))
      do j=1,6
        value=value.and.MA_free_heap(X(2,j))
      end do
      value = value.and.MA_free_heap(psi_r(2))
      value = value.and.MA_free_heap(psi2(2))
      value = value.and.MA_free_heap(psi(2))
      if (.not. value) 
     >   call errquit('pspw_wannier:deallocating heap memory',2,0)


 9000 if (oprint) call nwpw_message(4)
      call Parallel2d_Finalize()
      call Parallel_Finalize()
      pspw_wannier = value
      return 


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
  109 format(/' input filename : ',A)
  110 format(' output filename: ',A)

 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*              PSPW Wannier Generation             *')
 1030 FORMAT(10X,'*     [        Silvestrelli Algorithm      ]       *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #5.00   10/01/02              *')
 1041 FORMAT(10X,'*    Developed by Eric C. Brown, Kiril Tsemekhman  *')
 1042 FORMAT(10X,'*    (Jonsson group, U. Washington), Eric J.       *')
 1043 FORMAT(10X,'*    Bylaska and Patrick Nichols                   *')
 1044 FORMAT(10X,
 1100 FORMAT(//)
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
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

 1300 FORMAT(/' Wannier orbital generation for spin up orbitals')
 1301 FORMAT(/' Wannier orbital generation for spin down orbitals')
 1380 FORMAT(/' initial position of Wannier centers:')
 1381 FORMAT(/' final position of Wannier centers:')
 1390 FORMAT(5X,'psi',I5,' (',3F11.5,' ) - spin up orbital')
 1391 FORMAT(5X,'psi',I5,' (',3F11.5,' ) - spin down orbital')
 1776 FORMAT("MOLECULAR DIPOLES")
 1777 FORMAT("DX =",F11.5," DY= ",F11.5," DZ= ",F11.5)
      end
