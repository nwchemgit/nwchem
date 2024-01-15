*
* $Id$
*

*     **************************************
*     *                                    *
*     *           pawppv1                   *
*     *                                    *
*     **************************************

      logical function pawppv1(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,locp,lmax,rlocal)
      implicit none
      logical          oprint_in
      integer          version
      character*50     psp_filename,formatted_filename
      integer          ngrid(3)
      double precision unita(3,3)
      integer locp,lmax
      real*8  rlocal

#include "errquit.fh"
#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"

*     ***** local variables ****
      integer taskid,MASTER,msglen
      parameter (MASTER=0)


*     **** 1d pseudopotential data ****
      integer psp_type
      integer n1dgrid,nbasis,icut,nmax
      real*8  zv,r1,rmax,core_kin_energy,core_ion_energy
      real*8  zion,sigma,rcore,amass,rc(25)
      real*8  log_amesh,amesh,ecorez,ecorecore
      character*2 atom
      character*80 comment

      integer bprj(2),mprj(2),lprj(2),nprj(2),prj_ps0(2),prj_ps(2)
      integer v_ps(2),core_ps(2),core_ae(2)
      integer core_ps_prime(2),core_ae_prime(2)
      integer dphi_ps(2),phi_ps(2),dphi_ae(2),phi_ae(2)
      integer eig(2),lps(2),nps(2),nae(2)
      integer rgrid(2)
      integer nmaxl(2)

*     **** matrix data ****
      integer Gijl(2)
      integer hartree_matrix(2),comp_charge_matrix(2),comp_pot_matrix(2)

*     ***** ngrid data *****
      integer nproj,nsize,nfft1,nfft2,nfft3
      integer vl(2),vlpaw(2),vnl(2),G_indx,G_hndl
      integer f1(2),f2(2),f3(2),f4(2),cs(2),sn(2)

*     **** ray data ****
      logical filter
      integer nray,G_ray_hndl,tmp_ray_hndl
      integer vnl_ray_hndl,vl_ray_hndl,vlpaw_ray_hndl
      integer G_ray_indx,tmp_ray_indx
      integer vnl_ray_indx,vl_ray_indx,vlpaw_ray_indx
      real*8 filter_scale


*     **** other variables ****
      logical hprint,mprint,oprint,value
      integer idum,l,ii,i,j,ierr
      character*255 full_filename
      real*8 zcore,fourpi,unitg(3,3)

*     **** external functions ****
      logical  control_print,control_kbpp_filter
      external control_print,control_kbpp_filter
      real*8   log_integrate_def,log_coulomb0_energy,log_coulomb_energy
      external log_integrate_def,log_coulomb0_energy,log_coulomb_energy
      integer  kbpp_calc_nray
      external kbpp_calc_nray
      real*8   control_kbpp_filter_scale
      external control_kbpp_filter_scale


      call Parallel_taskid(taskid)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)
      oprint = (oprint_in.or.hprint)


      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')
      read(11,*) psp_type
      read(11,'(A)') atom
      read(11,*) zv
      read(11,*) r1
      read(11,*) rmax
      read(11,*) n1dgrid
      read(11,*) nbasis
      read(11,*) (rc(i),i=1,nbasis)
      read(11,*) icut
      read(11,'(A)') comment
      read(11,*) core_kin_energy
      end if

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,psp_type)
      call Parallel_Brdcst_values(MASTER,msglen,zv)
      call Parallel_Brdcst_values(MASTER,msglen,r1)
      call Parallel_Brdcst_values(MASTER,msglen,rmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,n1dgrid)
      call Parallel_Brdcst_ivalues(MASTER,msglen,nbasis)
      call Parallel_Brdcst_ivalues(MASTER,msglen,icut)
      call Parallel_Brdcst_values(MASTER,msglen,core_kin_energy)
      msglen = nbasis
      call Parallel_Brdcst_values(MASTER,msglen,rc)


*     **** define rgrid ****
      log_amesh = dlog(rmax/r1)/dble(n1dgrid-1)
      amesh     = dexp(log_amesh)
      value = BA_alloc_get(mt_dbl,n1dgrid,'rgrid',rgrid(2),rgrid(1))
      if (.not.value)  call errquit('pawppv1:out of heap',0,MA_ERR)
      dbl_mb(rgrid(1)) = r1
      do i=1,n1dgrid-1
        dbl_mb(rgrid(1)+i) = dbl_mb(rgrid(1)+i-1)*amesh
      end do

*     **** allocate rest of grid data ****
      value =           BA_alloc_get(mt_int,nbasis,'nae',nae(2),nae(1))
      value = value.and.BA_alloc_get(mt_int,nbasis,'nps',nps(2),nps(1))
      value = value.and.BA_alloc_get(mt_int,nbasis,'lps',lps(2),lps(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis,'eig',eig(2),eig(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis*n1dgrid,'phi_ae',
     >                              phi_ae(2),phi_ae(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis*n1dgrid,'dphi_ae',
     >                              dphi_ae(2),dphi_ae(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis*n1dgrid,'phi_ps',
     >                              phi_ps(2),phi_ps(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis*n1dgrid,'dphi_ps',
     >                              dphi_ps(2),dphi_ps(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'core_ae',
     >                              core_ae(2),core_ae(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'core_ps',
     >                              core_ps(2),core_ps(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'core_ae_prime',
     >                              core_ae_prime(2),core_ae_prime(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'core_ps_prime',
     >                              core_ps_prime(2),core_ps_prime(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'v_ps',
     >                              v_ps(2),v_ps(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis*n1dgrid,'prj_ps',
     >                              prj_ps(2),prj_ps(1))
      value = value.and.BA_alloc_get(mt_dbl,nbasis*n1dgrid,'prj_ps0',
     >                              prj_ps0(2),prj_ps0(1))
      if (.not.value)  call errquit('pawppv1:out of heap',1,MA_ERR)


      if (taskid.eq.MASTER) then
         read(11,*) (int_mb(nae(1)+j),
     >                                 dbl_mb(eig(1)+j),
     >                                 int_mb(nps(1)+j),
     >                                 int_mb(lps(1)+j),j=0,nbasis-1)
         read(11,*) ((dbl_mb(phi_ae(1)+i+j*n1dgrid),
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*) ((dbl_mb(dphi_ae(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*) ((dbl_mb(phi_ps(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*) ((dbl_mb(dphi_ps(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*) ((dbl_mb(prj_ps(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*) (dbl_mb(core_ae(1)+i),
     >                                 i=0,n1dgrid-1)
         read(11,*) (dbl_mb(core_ps(1)+i),
     >                                 i=0,n1dgrid-1)
         read(11,*) (dbl_mb(v_ps(1)+i),
     >                                 i=0,n1dgrid-1)
         read(11,*) sigma
         read(11,*) zion
         read(11,*) ((dbl_mb(prj_ps0(1)+i+j*n1dgrid),
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         close(11)
      end if


      msglen = nbasis
      call Parallel_Brdcst_ivalues(MASTER,msglen,int_mb(nae(1)))
      call Parallel_Brdcst_ivalues(MASTER,msglen,int_mb(nps(1)))
      call Parallel_Brdcst_ivalues(MASTER,msglen,int_mb(lps(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(eig(1)))

      msglen = nbasis*n1dgrid
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(phi_ae(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(dphi_ae(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(phi_ps(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(dphi_ps(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(prj_ps(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(prj_ps0(1)))

      msglen = n1dgrid
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(core_ae(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(core_ps(1)))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(v_ps(1)))

      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,sigma)
      call Parallel_Brdcst_values(MASTER,msglen,zion)
      

*     **** define nproj and lmax ****
      locp = -1
      lmax = -1
      nproj = 0
      do ii=1,nbasis
         l    = int_mb(lps(1)+ii-1)
         nproj = nproj + 2*l+1
         if (l.gt.lmax) lmax = l
      end do
      
*     **** define nmax ****
      if(.not. BA_push_get(mt_int,lmax+1,'nmaxl',nmaxl(2),nmaxl(1)))
     > call errquit('pawppv1:out of stack',2,MA_ERR)
      call icopy(lmax+1,0,0,int_mb(nmaxl(1)),1)
      do ii=1,nbasis
         l = int_mb(lps(1)+ii-1)
         int_mb(nmaxl(1)+l) = int_mb(nmaxl(1)+l) + 1
      end do
      nmax = 0
      do l=0,lmax
         if (int_mb(nmaxl(1)+l).gt.nmax) nmax = int_mb(nmaxl(1)+l)
      end do
      if(.not.BA_pop_stack(nmaxl(2)))
     > call errquit('pawppv1:error popping stack',0,MA_ERR)


*     **** allocate Gijl,Sijl,Tijl,vcore,Vpseuo ****
      l = nmax*nmax*(lmax+1)
      value=          BA_alloc_get(mt_dbl,5*l,'Gijl',Gijl(2),Gijl(1))
      !value=value.and.BA_alloc_get(mt_dbl,l,'Tijl',Tijl(2),Tijl(1))
      l = nbasis*nbasis*nbasis*nbasis*(2*lmax+1)
      value=value.and.
     >      BA_alloc_get(mt_dbl,l,'hartree_matrix',
     >                   hartree_matrix(2),hartree_matrix(1))
      l = nbasis*nbasis*(2*lmax+1)
      value=value.and.
     >      BA_alloc_get(mt_dbl,l,'comp_charge_matrix',
     >                   comp_charge_matrix(2),comp_charge_matrix(1))
      l = nbasis*nbasis*(2*lmax+1)
      value=value.and.
     >      BA_alloc_get(mt_dbl,l,'comp_pot_matrix',
     >                   comp_pot_matrix(2),comp_pot_matrix(1))

*     **** allocate nprj,lprj,mprj,bprj ****
      value=value.and.BA_alloc_get(mt_int,nproj,'nprj',nprj(2),nprj(1))
      value=value.and.BA_alloc_get(mt_int,nproj,'lprj',lprj(2),lprj(1))
      value=value.and.BA_alloc_get(mt_int,nproj,'mprj',mprj(2),mprj(1))
      value=value.and.BA_alloc_get(mt_int,nproj,'bprj',bprj(2),bprj(1))
      if (.not.value)  call errquit('pawppv1:out of heap',2,MA_ERR)


*    **** more temporary space ****
      value =           BA_alloc_get(mt_dbl,n1dgrid,'f1',f1(2),f1(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'f2',f2(2),f2(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'f3',f3(2),f3(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'f4',f4(2),f4(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'cs',cs(2),cs(1))
      value = value.and.BA_alloc_get(mt_dbl,n1dgrid,'sn',sn(2),sn(1))
      if (.not.value)call errquit('pawppv1:out of heap',0,MA_ERR)

*     **** allocate vl,vnl,vnlnrm G ****
      nsize = (ngrid(1)/2+1)*ngrid(2)*ngrid(3)
      value = BA_alloc_get(mt_dbl,nsize,'vl',vl(2),vl(1))
      value = value.and.BA_alloc_get(mt_dbl,nsize,
     >                               'vlpaw',vlpaw(2),vlpaw(1))
      value = value.and.BA_alloc_get(mt_dbl,nsize*(nproj),
     >                        'vnl',vnl(2), vnl(1))
      value = value.and.BA_alloc_get(mt_dbl,3*nsize,'G',G_hndl,G_indx)
      if (.not.value)call errquit('pawppv1:out of heap',0,MA_ERR)

*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp(nfft1,nfft2,nfft3,unita,unitg,dbl_mb(G_indx))
      filter = control_kbpp_filter()
      filter_scale = control_kbpp_filter_scale()

      !**** allocate memory for rays ****
      nray = kbpp_calc_nray(nfft1,nfft2,nfft3,unita)

      value =           BA_alloc_get(mt_dbl,nray,
     >                  'G_ray',G_ray_hndl,G_ray_indx)
      value = value.and.BA_alloc_get(mt_dbl,2*nray,
     >                  'vl_ray',vl_ray_hndl,vl_ray_indx)
      value = value.and.BA_alloc_get(mt_dbl,2*nray,
     >                  'vlpaw_ray',vlpaw_ray_hndl,vlpaw_ray_indx)
      value = value.and.BA_alloc_get(mt_dbl,2*nray*(nbasis),
     >                  'vnl_ray',vnl_ray_hndl,vnl_ray_indx)
      value = value.and.BA_alloc_get(mt_dbl,nray,
     >                  'tmp_ray',tmp_ray_hndl,tmp_ray_indx)
      if (.not.value)
     >   call errquit('pawppv1:out of heap memory',0,MA_ERR)

        call kbpp_generate_G_ray(nfft1,nfft2,nfft3,
     >                         unita,dbl_mb(G_ray_indx))

        zcore = zion - zv
        call integrate_pawppv1(version,rlocal,
     >                     n1dgrid,log_amesh,nbasis,nmax,lmax,icut,
     >                     zv,zcore,sigma,
     >                     int_mb(nps(1)),int_mb(lps(1)),
     >                                dbl_mb(v_ps(1)),
     >                                dbl_mb(prj_ps(1)),
     >                                dbl_mb(phi_ae(1)),
     >                                dbl_mb(dphi_ae(1)),
     >                                dbl_mb(phi_ps(1)),
     >                                dbl_mb(dphi_ps(1)),
     >                                dbl_mb(core_ae(1)),
     >                                dbl_mb(core_ps(1)),
     >                                dbl_mb(rgrid(1)),
     >                                dbl_mb(f1(1)),
     >                                dbl_mb(f2(1)),
     >                                dbl_mb(f3(1)),
     >                                dbl_mb(f4(1)),
     >                                dbl_mb(cs(1)),
     >                                dbl_mb(sn(1)),
     >                      nfft1,nfft2,nfft3,nproj,
     >                                dbl_mb(G_indx),
     >                                dbl_mb(vl(1)),
     >                                dbl_mb(vlpaw(1)),
     >                                dbl_mb(vnl(1)),
     >                                int_mb(nprj(1)),
     >                                int_mb(lprj(1)),
     >                                int_mb(mprj(1)),
     >                                int_mb(bprj(1)),
     >                                dbl_mb(Gijl(1)),
     >                                dbl_mb(hartree_matrix(1)),
     >                                dbl_mb(comp_charge_matrix(1)),
     >                                dbl_mb(comp_pot_matrix(1)),
     >                      nray,
     >                                dbl_mb(G_ray_indx),
     >                                dbl_mb(vl_ray_indx),
     >                                dbl_mb(vlpaw_ray_indx),
     >                                dbl_mb(vnl_ray_indx),
     >                                dbl_mb(tmp_ray_indx),
     >                                filter,filter_scale,
     >                      ierr)
        value = BA_free_heap(tmp_ray_hndl)
        value = value.and.BA_free_heap(vl_ray_hndl)
        value = value.and.BA_free_heap(vlpaw_ray_hndl)
        value = value.and.BA_free_heap(vnl_ray_hndl)
        value = value.and.BA_free_heap(G_ray_hndl)
        if (.not.value)
     >   call errquit('pawppv1:Error freeing memory',0,MA_ERR)


*      **** calculate radial derivatives of core densities ****
         call pawppv1_derivative_ngrid(
     >            n1dgrid,
     >            log_amesh,
     >            dbl_mb(rgrid(1)),
     >            dbl_mb(core_ae(1)),
     >            dbl_mb(core_ae_prime(1)))
         call pawppv1_derivative_ngrid(
     >            n1dgrid,
     >            log_amesh,
     >            dbl_mb(rgrid(1)),
     >            dbl_mb(core_ps(1)),
     >            dbl_mb(core_ps_prime(1)))

*     *** integrate core density ***
c      fourpi = 16.0d0*datan(1.0d0)
c      zcore = fourpi*log_integrate_def(0,
c     >                         dbl_mb(core_ae(1)),
c     >                         2,dbl_mb(rgrid(1)),log_amesh,n1dgrid)
c      write(luout,*) "Zcore=",zcore,rmax,dbl_mb(rgrid(1)+n1dgrid-1)
c      zcore = fourpi*log_integrate_def(0,
c     >                         dbl_mb(core_ps(1)),
c     >                         2,dbl_mb(rgrid(1)),log_amesh,n1dgrid)
c      write(luout,*) "Zcoreps=",zcore

c*     **** compute the core_ion_energy = ecorez + ecorecore ****
c      ecorez = -zion*fourpi
c     >              *log_integrate_def(0,dbl_mb(core_ae(1)),
c     >                                 1,dbl_mb(rgrid(1)),
c     >                                 log_amesh,n1dgrid)
      if (dabs(zcore).gt.1.0d-9) then
         ecorez = log_coulomb0_energy(dbl_mb(core_ae(1)),zion-zv,
     >                               dbl_mb(rgrid(1)),n1dgrid,log_amesh,
     >                               zion)
         ecorecore = log_coulomb_energy(dbl_mb(core_ae(1)),zion-zv,
     >                               dbl_mb(rgrid(1)),n1dgrid,log_amesh)
         core_ion_energy = ecorez + ecorecore
      else
         core_ion_energy = 0.0d0
      end if
c      write(*,*) "zv,zion,zcore=",zv,zion,zcore
c      write(*,*) "core_ion_energy=",core_ion_energy,ecorez,ecorecore



      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(formatted_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      if (mprint) then
      write(luout,*)
      write(luout,*) "Generated formatted_filename: ",full_filename(1:l)
      if (filter) write(luout,*) "- filtering pseudopotential -"
      end if
      call openfile(2,full_filename,l,'w',l)

         call cwrite(2,comment,80)
         call iwrite(2,psp_type,1)
         call iwrite(2,version,1)
         call iwrite(2,ngrid,3)
         call dwrite(2,unita,9)
         call cwrite(2,atom,2)
         amass = 0.0d0
         call dwrite(2,amass,1)
         call dwrite(2,zv,1)
         call iwrite(2,lmax,1)
         !call iwrite(2,locp,1)
         call iwrite(2,nbasis,1)

         call iwrite(2,nmax,1)
         call dwrite(2,rc,lmax+1)

         call iwrite(2,nproj,1)
         if (nproj.gt.0) then
         call iwrite(2,int_mb(nprj(1)),nproj)
         call iwrite(2,int_mb(lprj(1)),nproj)
         call iwrite(2,int_mb(mprj(1)),nproj)
         call iwrite(2,int_mb(bprj(1)),nproj)
         call dwrite(2,dbl_mb(Gijl(1)),5*nmax*nmax*(lmax+1))
         end if

         if (version.eq.4) call dwrite(2,rlocal,1)
         rcore = 0.0d0
         call dwrite(2,rcore,1)


*        **** other PAW matrices ****
         l = nbasis*nbasis*nbasis*nbasis*(2*lmax+1)
         call dwrite(2,dbl_mb(hartree_matrix(1)),l)
         l = nbasis*nbasis*(2*lmax+1)
         call dwrite(2,dbl_mb(comp_charge_matrix(1)),l)
         call dwrite(2,dbl_mb(comp_pot_matrix(1)),l)

*        **** miscelaneous PAW energies ****
         call dwrite(2,core_kin_energy,1)
         call dwrite(2,core_ion_energy,1)

*        **** write 1d-wavefunctions ****
         call iwrite(2,n1dgrid,1)
         call iwrite(2,icut,1)
         call dwrite(2,log_amesh,1)
         call dwrite(2,r1,1)
         call dwrite(2,rmax,1)
         call dwrite(2,sigma,1)
         call dwrite(2,zion,1)
         call dwrite(2,dbl_mb(eig(1)),nbasis)
         call iwrite(2,int_mb(nae(1)),nbasis)
         call iwrite(2,int_mb(nps(1)),nbasis)
         call iwrite(2,int_mb(lps(1)),nbasis)

         call dwrite(2,dbl_mb(rgrid(1)),n1dgrid)
         call dwrite(2,dbl_mb(phi_ae(1)),n1dgrid*nbasis)
         call dwrite(2,dbl_mb(dphi_ae(1)),n1dgrid*nbasis)
         call dwrite(2,dbl_mb(phi_ps(1)),n1dgrid*nbasis)
         call dwrite(2,dbl_mb(dphi_ps(1)),n1dgrid*nbasis)
         call dwrite(2,dbl_mb(core_ae(1)),n1dgrid)
         call dwrite(2,dbl_mb(core_ps(1)),n1dgrid)
         call dwrite(2,dbl_mb(core_ae_prime(1)),n1dgrid)
         call dwrite(2,dbl_mb(core_ps_prime(1)),n1dgrid)

         call dwrite(2,dbl_mb(vl(1)),nsize)
         call dwrite(2,dbl_mb(vlpaw(1)),nsize)
         if (nproj.gt.0) then
            call dwrite(2,dbl_mb(vnl(1)),nsize*nproj)
         end if

      call closefile(2)
      end if


      value = .true.
      value = value.and.BA_free_heap(G_hndl)
      value = value.and.BA_free_heap(vnl(2))
      value = value.and.BA_free_heap(vlpaw(2))
      value = value.and.BA_free_heap(vl(2))
      value = value.and.BA_free_heap(sn(2))
      value = value.and.BA_free_heap(cs(2))
      value = value.and.BA_free_heap(f1(2))
      value = value.and.BA_free_heap(f2(2))
      value = value.and.BA_free_heap(f3(2))
      value = value.and.BA_free_heap(f4(2))

      value = value.and.BA_free_heap(Gijl(2))
      value = value.and.BA_free_heap(hartree_matrix(2))
      value = value.and.BA_free_heap(comp_charge_matrix(2))
      value = value.and.BA_free_heap(comp_pot_matrix(2))

      value = value.and.BA_free_heap(bprj(2))
      value = value.and.BA_free_heap(mprj(2))
      value = value.and.BA_free_heap(lprj(2))
      value = value.and.BA_free_heap(nprj(2))

      value = value.and.BA_free_heap(prj_ps0(2))
      value = value.and.BA_free_heap(prj_ps(2))
      value = value.and.BA_free_heap(v_ps(2))
      value = value.and.BA_free_heap(core_ps_prime(2))
      value = value.and.BA_free_heap(core_ae_prime(2))
      value = value.and.BA_free_heap(core_ps(2))
      value = value.and.BA_free_heap(core_ae(2))
      value = value.and.BA_free_heap(dphi_ps(2))
      value = value.and.BA_free_heap(phi_ps(2))
      value = value.and.BA_free_heap(dphi_ae(2))
      value = value.and.BA_free_heap(phi_ae(2))

      value = value.and.BA_free_heap(eig(2))
      value = value.and.BA_free_heap(lps(2))
      value = value.and.BA_free_heap(nps(2))
      value = value.and.BA_free_heap(nae(2))

      value = value.and.BA_free_heap(rgrid(2))
      if (.not.value)  call errquit('pawppv1:freeing heap',5,MA_ERR)

      pawppv1 = value
      return

 9999 call errquit('pawppv1:Error reading psp_filename',0,DISK_ERR)
      pawppv1 = value
      return

      END




c     *************************************************
c     *                                               *
c     *        pawppv1_derivative_ngrid               *
c     *                                               *
c     *************************************************
c
c  This routine computes the seven point derivative of f.
c  where f and df are stored on a logarithmic grid. The
c  dimensions of f and df are, f(1:ng), and df(1:ng)

      subroutine pawppv1_derivative_ngrid(ng,log_amesh,r,f,df)
      implicit none
      integer           ng
      double precision  log_amesh
      double precision  r(ng)
      double precision  f(ng)
      double precision df(ng)

      double precision one_over_60
      parameter (one_over_60 = 1.0d0/60.0d0)

      integer i,n1,n2,m1,m2
      double precision aa

      aa = one_over_60/log_amesh
      n1 = 1
      n2 = ng
      m1 = n1
      m2 = n2


      if (n1.le.3) then
        if ((n1.eq.1).and.(n1.ge.m1).and.(n1.le.m2)) then
          df(1) = aa*(-147.0d0*f(1)
     >               + 360.0d0*f(2)
     >               - 450.0d0*f(3)
     >               + 400.0d0*f(4)
     >               - 225.0d0*f(5)
     >               +  72.0d0*f(6)
     >               -  10.0d0*f(7))/r(1)
          n1 = n1+1
        end if
        if ((n1.eq.2).and.(n1.ge.m1).and.(n1.le.m2)) then
          df(2) = aa*( -10.0d0*f(1)
     >               -  77.0d0*f(2)
     >               + 150.0d0*f(3)
     >               - 100.0d0*f(4)
     >               +  50.0d0*f(5)
     >               -  15.0d0*f(6)
     >               +   2.0d0*f(7))/r(2)
          n1 = n1+1
        end if
        if ((n1.eq.3.and.(n1.ge.m1).and.(n1.le.m2))) then
          df(3) = aa*(  +2.0d0*f(1)
     >               -  24.0d0*f(2)
     >               -  35.0d0*f(3)
     >               +  80.0d0*f(4)
     >               -  30.0d0*f(5)
     >               +   8.0d0*f(6)
     >               -   1.0d0*f(7))/r(3)
          n1 = n1+1
        end if
      end if

      if (n2.ge.(ng-2)) then
        if ((n2.eq.ng).and.(n2.ge.m1).and.(n2.le.m2)) then
          df(ng) = aa*( +147.0d0*f(ng)
     >                - 360.0d0*f(ng-1)
     >                + 450.0d0*f(ng-2)
     >                - 400.0d0*f(ng-3)
     >                + 225.0d0*f(ng-4)
     >                -  72.0d0*f(ng-5)
     >                +  10.0d0*f(ng-6))/r(ng)
          n2 = n2-1
        end if
        if ((n2.eq.(ng-1).and.(n2.ge.m1).and.(n2.le.m2))) then
          df(ng-1) = aa*( +10.0d0*f(ng)
     >                  +  77.0d0*f(ng-1)
     >                  - 150.0d0*f(ng-2)
     >                  + 100.0d0*f(ng-3)
     >                  -  50.0d0*f(ng-4)
     >                  +  15.0d0*f(ng-5)
     >                  -   2.0d0*f(ng-6))/r(ng-1)
          n2 = n2-1
        end if
        if ((n2.eq.(ng-2).and.(n2.ge.m1).and.(n2.le.m2))) then
          df(ng-2) = aa*(  -2.0d0*f(ng)
     >                  +  24.0d0*f(ng-1)
     >                  +  35.0d0*f(ng-2)
     >                  -  80.0d0*f(ng-3)
     >                  +  30.0d0*f(ng-4)
     >                  -   8.0d0*f(ng-5)
     >                  +   1.0d0*f(ng-6))/r(ng-2)
          n2 = n2-1
        end if
      end if

      do i=n1,n2
        df(i) = aa*(  -1.0d0*f(i-3)
     >             +   9.0d0*f(i-2)
     >             -  45.0d0*f(i-1)
     >             +  45.0d0*f(i+1)
     >             -   9.0d0*f(i+2)
     >             +   1.0d0*f(i+3))/r(i)
      end do

      return
      end


      


