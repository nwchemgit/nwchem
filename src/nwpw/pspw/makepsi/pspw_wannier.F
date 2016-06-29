*
* $Id$
*

c     **********************************************
c     *                                            *
c     *                pspw_wannier                *
c     *                                            *
c     **********************************************

c Computes wannier orbitals using the Silvestrelli algorithm.

      logical function pspw_wannier(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid,tmp_len
      parameter (MASTER=0,tmp_len=140)

      real*8 autoDebye
      parameter (autoDebye=2.5416d0)

      logical value,oprint
      integer mapping,mapping1d,maxiter
      integer i,j,ms,n,ne(2),ispin,ngrid(3),nemax
      integer neq(2),nemaxq,shift,rank
      integer nfft3d,n2ft3d,npack1
      integer info
      integer np,nion,ii

      integer psi(2),psi2(2),psi_r(2),psi_r2(2),Rwav(2)
      integer X(2,6),Xeig(2,6)
      integer A(2),Ar(2)
      real*8 bv(3,6),wrk(6,6),wts(6),bmat(3,3)
      real*8 b(3),ixmat(3,6)
      real*8 xx,yy,zz,tmp1(tmp_len),maxtime,pcharge
      real*8 dx(2),dy(2),dz(2),nx,ny,nz,tx,ty,tz,alpha,scal
      real*8 dix(2),diy(2),diz(2),tix,tiy,tiz
      real*8 drx(2),dry(2),drz(2),trx,try,trz
      complex*16 arg,wx,arg0
      character*2 labels(6)

c      double precision unita(3,3),b(3),alpha

c      complex*16 wx,wy,wz

c      character*50 new_wavefunction_filename
c      character*50 old_wavefunction_filename
c      character*255 full_filename,full_filename2

*     **** external functions ****
      integer ion_katm,ion_nion
      real*8 ion_rion,psp_zv
      logical      control_read
      integer      control_ngrid,control_ispin,control_mapping
      integer      control_mapping1d
      integer      pack_nwave,pack_nwave_all
      integer      control_np_orbital,control_wannier_maxiter
      real*8       lattice_unita,lattice_unitg,control_wannier_timestep
      real*8       lattice_ecut,lattice_wcut,lattice_omega
      character*50 control_input_psi,control_output_psi
      external     control_read
      external     control_ngrid,control_ispin,control_mapping
      external     control_mapping1d
      external     pack_nwave,pack_nwave_all
      external     control_np_orbital,control_wannier_maxiter
      external     lattice_unita,lattice_unitg,control_wannier_timestep
      external     lattice_ecut,lattice_wcut,lattice_omega
      external     control_input_psi,control_output_psi
      external     ion_katm,ion_nion
      external     ion_rion,psp_zv
      logical  Dneall_m_allocate,Dneall_w_allocate,Dneall_m_free
      external Dneall_m_allocate,Dneall_w_allocate,Dneall_m_free
      integer     ion_nearest_index
      external    ion_nearest_index
      character*2 ion_aname_nocarat
      external    ion_aname_nocarat

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
         write(luout,1035)
         write(luout,1010)
         write(luout,1037)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
         write(luout,1042)
         write(luout,1043)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
      end if
      
      value = control_read(9,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      mapping = control_mapping()
      maxiter = control_wannier_maxiter()
      maxtime = control_wannier_timestep()

      call Parallel2d_Init(control_np_orbital())

*     **** initialize psi_data ****
      call psi_data_init(100)

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

      call ion_init(rtdb)
      call psp_init()
      call psp_readall()



*     ***** allocate psi and psi_r wavefunctions ****
      call psi_get_ne(ispin,ne)
      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)
      nemaxq = neq(1)+neq(2)

      
      value = BA_alloc_get(mt_dcpl,npack1*nemaxq,
     >                     'psi',psi(2),psi(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*nemaxq,
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d*nemaxq,
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d*nemaxq,
     >                     'psi_r2',psi_r2(2),psi_r2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,3*(ne(1)+ne(2)),
     >                     'Rwav',Rwav(2),Rwav(1))
      if (.not. value) 
     >   call errquit('pspw_wannier:out of heap memory',0,0)


*     *****  read psi2 wavefunctions ****
      call psi_read_noocc(ispin,ne,dcpl_mb(psi(1)))


      if (oprint) then
  
         write(luout,109) control_input_psi()
         write(luout,110) control_output_psi()

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

         write(luout,1280) maxtime
         write(luout,1281) maxiter
      end if


*     **** psi_r <-- FFT(psi) - complex to real FFT ****
      !nemax = ne(1)+ne(2)
      do n=1,nemaxq
         call Pack_c_Copy(1,dcpl_mb(psi(1)  +(n-1)*npack1),
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      end do
      call Grsm_gh_fftb(nfft3d,nemaxq,dbl_mb(psi_r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemaxq,dbl_mb(psi_r(1)))


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
c      do j=1,6
c      value = value.and.BA_alloc_get(mt_dcpl,
c     >                     n*n,
c     >                     labels(j),X(2,j),X(1,j))      
c      end do
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,n*n,
c     >                     'A',A(2),A(1))

      do j=1,6
         value = value.and.Dneall_w_allocate(1,X(1,j))
         value = value.and.
     >           BA_alloc_get(mt_dcpl,n,'Xeig',Xeig(2,j),Xeig(1,j))
      end do
      value = value.and.Dneall_w_allocate(1,A)
      value = value.and.Dneall_m_allocate(1,Ar)
      if (.not. value) 
     >   call errquit('pspw_wannier:out of heap memory',1,0)

c     *** Silvestrelli G1 ***
      ixmat(1,1)=1.0d0
      ixmat(2,1)=0.0d0
      ixmat(3,1)=0.0d0

c     *** Silvestrelli G4 ***
      ixmat(1,2)=1.0d0
      ixmat(2,2)=1.0d0
      ixmat(3,2)=0.0d0

c     *** Silvestrelli G5 ***
      ixmat(1,3)=1.0d0
      ixmat(2,3)=0.0d0
      ixmat(3,3)=1.0d0

c     *** Silvestrelli G2 ***
      ixmat(1,4)=0.0d0
      ixmat(2,4)=1.0d0
      ixmat(3,4)=0.0d0

c     *** Silvestrelli G6 ***
      ixmat(1,5)=0.0d0
      ixmat(2,5)=1.0d0
      ixmat(3,5)=1.0d0

c     *** Silvestrelli G3 ***
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

*     *** scal=(2*pi/L)**2 ***
      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0*datan(1.0d0)/scal
      scal = scal*scal
      wts(1)=scal
      wts(4)=scal
      wts(6)=scal
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

      nx=0.0d0
      ny=0.0d0
      nz=0.0d0
      pcharge = 0.0d0
      nion=ion_nion()
      do i=1,nion
        j=ion_katm(i)
        nx=nx+psp_zv(j)*ion_rion(1,i)
        ny=ny+psp_zv(j)*ion_rion(2,i)
        nz=nz+psp_zv(j)*ion_rion(3,i)
        pcharge = pcharge + psp_zv(j)
      end do   

      dx(1)=0.0d0
      dx(2)=0.0d0
      dy(1)=0.0d0
      dy(2)=0.0d0
      dz(1)=0.0d0
      dz(2)=0.0d0

      dix(1)=0.0d0
      dix(2)=0.0d0
      diy(1)=0.0d0
      diy(2)=0.0d0
      diz(1)=0.0d0
      diz(2)=0.0d0

      drx(1)=0.0d0
      drx(2)=0.0d0
      dry(1)=0.0d0
      dry(2)=0.0d0
      drz(1)=0.0d0
      drz(2)=0.0d0

      if (oprint) then
         write(*,1260) 
         write(*,1261) rank
         do i=1,rank
            write(*,1262) i,bv(1,i),bv(2,i),bv(3,i),wts(i)
         end do
      end if
      
      do ms=1,ispin

        if ((oprint).and.(ms.eq.1)) write(luout,1300)
        if ((oprint).and.(ms.eq.2)) write(luout,1301)


        do i=1,rank
          b(1) = bv(1,i)
          b(2) = bv(2,i)
          b(3) = bv(3,i)
          call silvestrelli_overlap(
     >                  b,ms,ne,neq,
     >                  dbl_mb(psi_r(1)),
     >                  dbl_mb(psi_r2(1)),
     >                  dcpl_mb(X(1,i)))
          call Dneall_w_eigenvalues(ms,dcpl_mb(X(1,i)),
     >                              dcpl_mb(Xeig(1,i)))
        end do

        if (oprint) then
          write(luout,1380)
          do i=1,ne(ms)
            shift=(i-1)*ne(ms)+(i-1)
            xx=0.0d0
            yy=0.0d0
            zz=0.0d0
            do j=1,rank

               !*** really just want complex eigenvalues of X here ***
               arg0=dcpl_mb(Xeig(1,j)+(i-1))
               arg0= -wts(j)*datan2(dimag(arg0),dble(arg0))
               drx(ms)=drx(ms)+bv(1,j)*arg0/scal
               dry(ms)=dry(ms)+bv(2,j)*arg0/scal
               drz(ms)=drz(ms)+bv(3,j)*arg0/scal

               arg=dcpl_mb(X(1,j)+shift)
               arg= -wts(j)*datan2(dimag(arg),dble(arg))
               xx=xx+bv(1,j)*arg/scal
               yy=yy+bv(2,j)*arg/scal
               zz=zz+bv(3,j)*arg/scal
            end do
            ii = ion_nearest_index(xx,yy,zz)
            dix(ms)=dix(ms)+xx
            diy(ms)=diy(ms)+yy
            diz(ms)=diz(ms)+zz

            if (ms.eq.1)
     >        write(luout,1390) i,xx,yy,zz,ion_aname_nocarat(ii),ii


            if (ms.eq.2)
     >        write(luout,1391) i+ne(1),xx,yy,zz,
     >                          ion_aname_nocarat(ii),ii
          end do
          write(luout,1800)
        end if

        call silvestrelli_minimize(maxiter,ms,ne(ms),rank,wts,
     >                             dcpl_mb(X(1,1)),dcpl_mb(X(1,2)),
     >                             dcpl_mb(X(1,3)),dcpl_mb(X(1,4)),
     >                             dcpl_mb(X(1,5)),dcpl_mb(X(1,6)),
     >                             dcpl_mb(A(1)))

        if (oprint) write(luout,1381)
        do i=1,ne(ms)
           !shift=i-1
           !shift=shift*ne(ms)+shift
           shift=(i-1)*ne(ms)+(i-1)
           xx=0.0d0
           yy=0.0d0
           zz=0.0d0
           do j=1,rank
              arg=dcpl_mb(X(1,j)+shift)
              arg= -wts(j)*datan2(dimag(arg),dble(arg))
              xx=xx+bv(1,j)*arg/scal
              yy=yy+bv(2,j)*arg/scal
              zz=zz+bv(3,j)*arg/scal
           end do
           dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+i-1))   = xx
           dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+i-1)+1) = yy
           dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+i-1)+2) = zz
           ii = ion_nearest_index(xx,yy,zz)
           dx(ms)=dx(ms)+xx
           dy(ms)=dy(ms)+yy
           dz(ms)=dz(ms)+zz
           if (oprint) then
              if (ms.eq.1) then
                write(luout,1390) i,xx,yy,zz,ion_aname_nocarat(ii),ii
                !write(15,*) "H ",xx*0.529177d0,yy*0.529177d0,zz*0.529177d0
              end if

              if (ms.eq.2) then
                write(luout,1391) i+ne(1),xx,yy,zz,
     >                            ion_aname_nocarat(ii),ii
                !write(15,*) "H ",xx*0.529177d0,yy*0.529177d0,zz*0.529177d0
              end if
           end if
        end do

*       **** rotate psi to produce Wannier orbitals ****
        call Dneall_wtom_Real(ms,dcpl_mb(A(1)),dbl_mb(Ar(1)))
        call Dneall_mm_transpose(ms,dbl_mb(Ar(1)),dcpl_mb(A(1)))
        call Dneall_fmf_Multiply(ms,
     >                  dcpl_mb(psi(1)),npack1,
     >                  dcpl_mb(A(1)),1.0d0,
     >                  dcpl_mb(psi2(1)),
     >                  0.0d0)
      end do !* ms *

cccccccccccccccccccccccccccccccccccccccccccccc
c  Molecular dipoles from Resta's theory!
ccccccccccccccccccccccccccccccccccccccccccccc
      tx=nx-dx(1)-dx(ispin)
      ty=ny-dy(1)-dy(ispin)
      tz=nz-dz(1)-dz(ispin)
      xx = dsqrt(tx*tx + ty*ty + tz*tz)

      tix=nx-dix(1)-dix(ispin)
      tiy=ny-diy(1)-diy(ispin)
      tiz=nz-diz(1)-diz(ispin)
      yy = dsqrt(tix*tix + tiy*tiy + tiz*tiz)

      trx=nx-drx(1)-drx(ispin)
      try=ny-dry(1)-dry(ispin)
      trz=nz-drz(1)-drz(ispin)
      zz = dsqrt(trx*trx + try*try + trz*trz)

      if (oprint) then
         write(luout,1771)
         write(luout,1772) 'spin up   ',
     >                     dx(1)/dble(ne(1)),
     >                     dy(1)/dble(ne(1)),
     >                     dz(1)/dble(ne(1))
         if (ne(ispin).ne.0)
     >      write(luout,1772) 'spin down ',
     >                        dx(ispin)/dble(ne(ispin)),
     >                        dy(ispin)/dble(ne(ispin)),
     >                        dz(ispin)/dble(ne(ispin))
         write(luout,1772) 'electronic',
     >                      (dx(1)+dx(ispin))/dble(ne(1)+ne(ispin)),
     >                      (dy(1)+dy(ispin))/dble(ne(1)+ne(ispin)),
     >                      (dz(1)+dz(ispin))/dble(ne(1)+ne(ispin))
         write(luout,1772) 'ionic     ',
     >                      nx/pcharge,
     >                      ny/pcharge,
     >                      nz/pcharge
         write(luout,1773)
         write(luout,1774) tx,ty,tz
         write(luout,1775) xx,xx*autoDebye

         write(luout,1779)
         write(luout,1774) tix,tiy,tiz
         write(luout,1775) yy,yy*autoDebye

         write(luout,1778)
         write(luout,1774) trx,try,trz
         write(luout,1775) zz,zz*autoDebye
      end if


*     ***** write psi2 wavefunctions ****
      call psi_write_noocc(ispin,ne,dcpl_mb(psi2(1)))

c*     ***** write basis functions *****
c      call psi_makebasisfunctions(1,ispin,ne,neq,
c     >                            dbl_mb(Rwav(1)),
c     >                            dcpl_mb(psi2(1)),
c     >                            dcpl_mb(psi(1)))
c

*     **** deallocate heap memory ****
      call ion_write(rtdb)
      call ion_end()
      call psp_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call D3dB_pfft_end()
      call D3dB_end(1)
      call Dne_end()
      call psi_data_end()

c      value =  BA_free_heap(A(2))
c      do j=1,6
c        value=value.and.BA_free_heap(X(2,j))
c      end do
      value = Dneall_m_free(Ar)
      value = value.and.Dneall_m_free(A)
      do j=1,6
      value = value.and.Dneall_m_free(X(1,j))
      value = value.and.BA_free_heap(Xeig(2,j))
      end do
      value = value.and.BA_free_heap(Rwav(2))
      value = value.and.BA_free_heap(psi_r2(2))
      value = value.and.BA_free_heap(psi_r(2))
      value = value.and.BA_free_heap(psi2(2))
      value = value.and.BA_free_heap(psi(2))
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
 1037 FORMAT(10x,'*     Maximally Localized General Symmetry         *')
 1040 FORMAT(10X,'*            version #5.00   10/01/02              *')
 1041 FORMAT(10X,'*    Developed by Eric C. Brown, Kiril Tsemekhman  *')
 1042 FORMAT(10X,'*    (Jonsson group, U. Washington), Eric J.       *')
 1043 FORMAT(10X,'*    Bylaska and Patrick Nichols                   *')
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

 1260 FORMAT(/5x,' Silvestrelli reciprocal lattice vectors:')
 1261 FORMAT(5x,' rank = ',I1)
 1262 FORMAT(5x,' bs',I1,' =<',3f8.3,' > weight=',f8.3)
 1280 FORMAT(5X, ' maximum time step=',F10.4)
 1281 FORMAT(5X, ' maximum iterations =',I8)
 1300 FORMAT(/' Wannier orbital generation for spin up orbitals')
 1301 FORMAT(/' Wannier orbital generation for spin down orbitals')
 1380 FORMAT(/' initial position of Wannier centers:')
 1381 FORMAT(/' final position of Wannier centers:')
 1390 FORMAT(5X,'psi',I5,' (',3F11.5,
     >       ' ) - spin up orbital, nearest ion=',A3,I5)
 1391 FORMAT(5X,'psi',I5,' (',3F11.5,
     >       ' ) - spin down orbital, nearest ion=',A3,I5)
 1771 FORMAT(//'== Center of Charge =='/)
 1772 FORMAT(A10,'  (',F10.4,',',F10.4,',',F10.4,' )')
 1773 FORMAT(//'== Wannier Crystal Dipole =='/)
 1774 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1775 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')
 1776 FORMAT(/"ELECTRONIC DIPOLES")
 1777 FORMAT("DX =",F11.5," DY= ",F11.5," DZ= ",F11.5)
 1778 FORMAT(//'== Resta (not localized) Crystal Dipole =='/)
 1779 FORMAT(//'== Initial Wannier Crystal Dipole =='/)
 1780 FORMAT("NUCLEAR DIPOLES")
 1785 FORMAT("TOTAL DIPOLES")
 1800 FORMAT(/2x,'iter.',8x,'Omega ',6x,' maxgrad',9x,'domega ',5x,
     >       'time_step'/
     >     2x,'------------------------------------------------------',
     >        '---------')

 
      end
