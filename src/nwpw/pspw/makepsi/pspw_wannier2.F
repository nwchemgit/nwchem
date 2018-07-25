*
* $Id: pspw_wannier2.F 26433 2014-12-03 23:34:20Z mjacquelin $
*

c     **********************************************
c     *                                            *
c     *                pspw_wannier2               *
c     *                                            *
c     **********************************************

c Computes wannier orbitals using the Damle-Lin-Ying algorithm.

      logical function pspw_wannier2(rtdb)
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

      integer psi(2),psi_r(2),psi_r2(2),Rwav(2)
      integer X(2,6)
      integer A(2),Ar(2),rotmat1(2),rotmat2(2)
      real*8 bv(3,6),wrk(6,6),wts(6),bmat(3,3)
      real*8 b(3),ixmat(3,6)
      real*8 xx,yy,zz,tmp1(tmp_len),maxtime,pcharge
      real*8 dx(2),dy(2),dz(2),nx,ny,nz,tx,ty,tz,alpha,scal,scal1,scal2
      complex*16 arg,wx
      character*2 labels(6)
      integer xs(6),xsgrid(2),tmp2(2)
      integer q,n1(2),n2(2)
      real*8  aa(3,3)
      real*8  cdx1,cdx2,cdx3,cdy1,cdy2,cdy3,t1,t2,t3,twopi


c      double precision unita(3,3),b(3),alpha

c      complex*16 wx,wy,wz

c      character*50 new_wavefunction_filename
c      character*50 old_wavefunction_filename
c      character*255 full_filename,full_filename2

*     **** external functions ****
      integer ion_katm,ion_nion
      real*8  ion_rion,psp_zv
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
      real*8      Dneall_m_value
      external    Dneall_m_value

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
      call ga_sync()

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
     >        BA_alloc_get(mt_dbl,n2ft3d*nemaxq,
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d*nemaxq,
     >                     'psi_r2',psi_r2(2),psi_r2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,3*(ne(1)+ne(2)),
     >                     'Rwav',Rwav(2),Rwav(1))
      if (.not. value) 
     >   call errquit('pspw_wannier2:out of heap memory',0,0)


*     *****  read psi wavefunctions ****
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
      call D3dB_nx(1,i)
      call D3dB_ny(1,j)
      call D3dB_nz(1,ii)
      scal1 = 1.0d0/dble(i*j*ii)
      call psi_dmatrix_localize(ispin,neq,n2ft3d,
     >                          dbl_mb(psi_r(1)),
     >                          dbl_mb(psi_r2(1)))
      call dcopy(n2ft3d*nemaxq,dbl_mb(psi_r2(1)),1,dbl_mb(psi_r(1)),1)
      call Grsm_hg_fftf0(nfft3d,nemaxq,dbl_mb(psi_r(1)))
      do n=1,nemaxq
         call Pack_c_Copy(1,dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                      dcpl_mb(psi(1) +(n-1)*npack1))
         call Pack_c_SMul1(1,scal1, dcpl_mb(psi(1)+(n-1)*npack1))
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
      

*     ***** compute the centers of the localized orbitals *****
      value = BA_push_get(mt_dbl,n2ft3d,"tmp2",tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,6*n2ft3d,"xsgrid",xsgrid(2),xsgrid(1))
      if (.not.value)
     >    call errquit('pspw_wannier2:push stack',1,MA_ERR)

      twopi = 8.0d0*datan(1.0d0)

      aa(1,1) = lattice_unita(1,1)
      aa(2,1) = lattice_unita(2,1)
      aa(3,1) = lattice_unita(3,1)
      aa(1,2) = lattice_unita(1,2)
      aa(2,2) = lattice_unita(2,2)
      aa(3,2) = lattice_unita(3,2)
      aa(1,3) = lattice_unita(1,3)
      aa(2,3) = lattice_unita(2,3)
      aa(3,3) = lattice_unita(3,3)

      call lattice_p_grid(dbl_mb(xsgrid(1)))

      xs(1) = xsgrid(1)
      xs(2) = xsgrid(1)+n2ft3d
      xs(3) = xsgrid(1)+2*n2ft3d
      xs(4) = xsgrid(1)+3*n2ft3d
      xs(5) = xsgrid(1)+4*n2ft3d
      xs(6) = xsgrid(1)+5*n2ft3d

      n1(1) = 1
      n1(2) = neq(1)+1
      n2(1) = neq(1)
      n2(2) = neq(1)+neq(2)
      call dcopy(3*(ne(1)+ne(2)),0.0d0,0,dbl_mb(Rwav(1)),1)
      do ms=1,ispin
         do q=n1(ms),n2(ms)
            call Dneall_qton(q,n)
            call D3dB_rr_Sqr(1,dbl_mb(psi_r2(1)+(q-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)))
            call D3dB_r_Zero_Ends(1,dbl_mb(tmp2(1)))

            call D3dB_rr_dot(1,dbl_mb(xs(1)),dbl_mb(tmp2(1)),cdx1)
            call D3dB_rr_dot(1,dbl_mb(xs(2)),dbl_mb(tmp2(1)),cdy1)

            call D3dB_rr_dot(1,dbl_mb(xs(3)),dbl_mb(tmp2(1)),cdx2)
            call D3dB_rr_dot(1,dbl_mb(xs(4)),dbl_mb(tmp2(1)),cdy2)

            call D3dB_rr_dot(1,dbl_mb(xs(5)),dbl_mb(tmp2(1)),cdx3)
            call D3dB_rr_dot(1,dbl_mb(xs(6)),dbl_mb(tmp2(1)),cdy3)
            t1 = datan2(cdy1*scal1,cdx1*scal1)/twopi
            t2 = datan2(cdy2*scal1,cdx2*scal1)/twopi
            t3 = datan2(cdy3*scal1,cdx3*scal1)/twopi
            xx = aa(1,1)*t1 + aa(1,2)*t2 + aa(1,3)*t3
            yy = aa(2,1)*t1 + aa(2,2)*t2 + aa(2,3)*t3
            zz = aa(3,1)*t1 + aa(3,2)*t2 + aa(3,3)*t3
            !dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+n-1))   = xx
            !dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+n-1)+1) = yy
            !dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+n-1)+2) = zz
            dbl_mb(Rwav(1)+3*(n-1))   = xx
            dbl_mb(Rwav(1)+3*(n-1)+1) = yy
            dbl_mb(Rwav(1)+3*(n-1)+2) = zz
         end do
      end do
      call D1dB_Vector_SumAll(3*(ne(1)+ne(2)),dbl_mb(Rwav(1)))

      value =           BA_pop_stack(xsgrid(2))
      value = value.and.BA_pop_stack(tmp2(2))
      if (.not.value)
     >    call errquit('pspw_wannier2:pop stack',2,MA_ERR)


      dx(1) = 0.0d0
      dx(2) = 0.0d0
      dy(1) = 0.0d0
      dy(2) = 0.0d0
      dz(1) = 0.0d0
      dz(2) = 0.0d0
      if (oprint) write(luout,1382)
      do ms=1,ispin
         do i=1,ne(ms)
          xx = dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+i-1))  
          yy = dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+i-1)+1)
          zz = dbl_mb(Rwav(1)+3*((ms-1)*ne(1)+i-1)+2)
          ii = ion_nearest_index(xx,yy,zz)
          dx(ms)=dx(ms)+xx
          dy(ms)=dy(ms)+yy
          dz(ms)=dz(ms)+zz
          if (oprint) then
             if (ms.eq.1) then
               write(luout,1390) i,xx,yy,zz,ion_aname_nocarat(ii),ii
               !write(15,*)"H ",xx*0.529177d0,yy*0.529177d0,zz*0.529177d0
             end if
             if (ms.eq.2) then
               write(luout,1391) i+ne(1),xx,yy,zz,
     >                           ion_aname_nocarat(ii),ii
               !write(15,*) "H ",xx*0.529177d0,yy*0.529177d0,zz*0.529177d0
             end if
          end if
         end do
      end do
      if (oprint) then
          write(luout,1778)
          write(luout,1772) 'spin up   ',
     >                      dx(1)/dble(ne(1)),
     >                      dy(1)/dble(ne(1)),
     >                      dz(1)/dble(ne(1))
          if (ne(ispin).ne.0)
     >       write(luout,1772) 'spin down ',
     >                         dx(ispin)/dble(ne(ispin)),
     >                         dy(ispin)/dble(ne(ispin)),
     >                         dz(ispin)/dble(ne(ispin))
          write(luout,1772) 'electronic',
     >                       (dx(1)+dx(ispin))/dble(ne(1)+ne(ispin)),
     >                       (dy(1)+dy(ispin))/dble(ne(1)+ne(ispin)),
     >                       (dz(1)+dz(ispin))/dble(ne(1)+ne(ispin))
          write(luout,1772) 'ionic     ',
     >                       nx/pcharge,
     >                       ny/pcharge,
     >                       nz/pcharge
      end if
      tx=nx-dx(1)-dx(ispin)
      ty=ny-dy(1)-dy(ispin)
      tz=nz-dz(1)-dz(ispin)
      xx = dsqrt(tx*tx + ty*ty + tz*tz)

      if (oprint) then
         write(luout,1779)
         write(luout,1774) tx,ty,tz
         write(luout,1775) xx,xx*autoDebye
         write(luout,*)
         write(luout,*)
      end if


*     ***** write psi2 wavefunctions ****
      call psi_write_noocc(ispin,ne,dcpl_mb(psi(1)))

*     **** deallocate heap memory ****
      !call ion_write(rtdb)
      call ion_destroy()
      call ion_end()
      call psp_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call D3dB_pfft_end()
      call D3dB_end(1)
      call Dne_end()
      call psi_data_end()


      value = value.and.BA_free_heap(Rwav(2))
      value = value.and.BA_free_heap(psi_r2(2))
      value = value.and.BA_free_heap(psi_r(2))
      value = value.and.BA_free_heap(psi(2))
      if (.not. value) 
     >   call errquit('pspw_wannier2:deallocating heap memory',2,0)


 9000 if (oprint) call nwpw_message(4)
      call Parallel2d_Finalize()
      call Parallel_Finalize()
      pspw_wannier2 = value
      return 


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
  109 format(/' input filename : ',A)
  110 format(' output filename: ',A)

 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*              PSPW Wannier Generation             *')
 1030 FORMAT(10X,'*     [        Damle-Lin-Ying Algorithm      ]     *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1037 FORMAT(10x,'*     Maximally Localized General Symmetry         *')
 1040 FORMAT(10X,'*            version #5.00    6/20/16              *')
 1041 FORMAT(10X,'*    Developed by Eric J. Bylaska and Patrick      *')
 1042 FORMAT(10X,'*    Nichols(PNNL), Anil Damle, Lexing Ying        *')
 1043 FORMAT(10X,'*    (Stanford), and Lin Lin (UC Berkeley)         *')
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
 1380 FORMAT(/' position of Wannier centers:')
 1381 FORMAT(/' position of Wannier centers:')
 1382 FORMAT(/' position of localized orbital centers:')
 1390 FORMAT(5X,'psi',I5,' (',3F11.5,
     >       ' ) - spin up orbital, nearest ion=',A3,I5)
 1391 FORMAT(5X,'psi',I5,' (',3F11.5,
     >       ' ) - spin down orbital, nearest ion=',A3,I5)
 
 1771 FORMAT(//'== Center of Charge =='/)
 1772 FORMAT(A10,'  (',F10.4,',',F10.4,',',F10.4,' )')
 1773 FORMAT(//'== Crystal Dipole =='/)
 1774 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1775 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')


 1776 FORMAT(/"ELECTRONIC DIPOLES")
 1777 FORMAT("DX =",F11.5," DY= ",F11.5," DZ= ",F11.5)
 1778 FORMAT(//'== Center of Charge from localized orbitals =='/)
 1779 FORMAT(//'== Crystal Dipole from localized orbitals =='/)
 1780 FORMAT("NUCLEAR DIPOLES")
 1785 FORMAT("TOTAL DIPOLES")
 1800 FORMAT(/2x,'iter.',8x,'Omega ',6x,' maxgrad',9x,'domega ',5x,
     >       'time_step'/
     >     2x,'------------------------------------------------------',
     >        '---------')

 
      end
