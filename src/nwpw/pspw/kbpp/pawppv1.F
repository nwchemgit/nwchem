*
* $Id: pawppv1.F,v 1.4 2007-12-18 02:26:23 bylaska Exp $
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
      character*25     psp_filename,formatted_filename
      integer          ngrid(3)
      double precision unita(3,3)
      integer locp,lmax
      real*8  rlocal

#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"

*     ***** local variables ****
      integer taskid,MASTER,msglen
      parameter (MASTER=0)


*     **** 1d pseudopotential data ****
      integer psp_type
      integer n1dgrid,nbasis,icut,nmax
      real*8  zv,r1,rmax,core_kin_energy,zion,sigma,rcore,amass,rc(25)
      real*8  log_amesh,amesh
      character*2 atom
      character*80 comment

      integer mprj(2),lprj(2),nprj(2),prj_ps0(2),prj_ps(2),v_ps(2)
      integer core_ps(2),core_ae(2)
      integer dphi_ps(2),phi_ps(2),dphi_ae(2),phi_ae(2)
      integer eig(2),lps(2),nps(2),nae(2)
      integer rgrid(2)
      integer nmaxl(2)

*     **** matrix data ****
      integer Gijl(2)

*     ***** ngrid data *****
      integer nproj,nsize,nfft1,nfft2,nfft3
      integer vl(2),vnl(2),G_indx,G_hndl
      integer f(2),cs(2),sn(2)

*     **** ray data ****
      logical filter
      integer nray,G_ray_hndl,tmp_ray_hndl
      integer vnl_ray_hndl,vl_ray_hndl
      integer G_ray_indx,tmp_ray_indx
      integer vnl_ray_indx,vl_ray_indx


*     **** other variables ****
      logical hprint,mprint,oprint,value
      integer idum,l,ii,i,j,ierr
      character*255 full_filename
      real*8 zcore,fourpi,unitg(3,3)

*     **** external functions ****
      logical  control_print,control_kbpp_filter
      external control_print,control_kbpp_filter
      real*8   log_integrate_def
      external log_integrate_def
      integer  kbpp_calc_nray
      external kbpp_calc_nray


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
      read(11,*,ERR=9999,END=9999) psp_type
      read(11,'(A)',ERR=9999,END=9999) atom
      read(11,*,ERR=9999,END=9999) zv
      read(11,*,ERR=9999,END=9999) r1
      read(11,*,ERR=9999,END=9999) rmax
      read(11,*,ERR=9999,END=9999) n1dgrid
      read(11,*,ERR=9999,END=9999) nbasis
      read(11,*,ERR=9999,END=9999) (rc(i),i=1,nbasis)
      read(11,*,ERR=9999,END=9999) icut
      read(11,'(A)',ERR=9999,END=9999) comment
      read(11,*,ERR=9999,END=9999) core_kin_energy
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
      value = MA_alloc_get(mt_dbl,n1dgrid,'rgrid',rgrid(2),rgrid(1))
      if (.not.value)  call errquit('pawppv1:out of heap',0,MA_ERR)
      dbl_mb(rgrid(1)) = r1
      do i=1,n1dgrid-1
        dbl_mb(rgrid(1)+i) = dbl_mb(rgrid(1)+i-1)*amesh
      end do

*     **** allocate rest of grid data ****
      value =           MA_alloc_get(mt_int,nbasis,'nae',nae(2),nae(1))
      value = value.and.MA_alloc_get(mt_int,nbasis,'nps',nps(2),nps(1))
      value = value.and.MA_alloc_get(mt_int,nbasis,'lps',lps(2),lps(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis,'eig',eig(2),eig(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis*n1dgrid,'phi_ae',
     >                              phi_ae(2),phi_ae(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis*n1dgrid,'dphi_ae',
     >                              dphi_ae(2),dphi_ae(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis*n1dgrid,'phi_ps',
     >                              phi_ps(2),phi_ps(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis*n1dgrid,'dphi_ps',
     >                              dphi_ps(2),dphi_ps(1))
      value = value.and.MA_alloc_get(mt_dbl,n1dgrid,'core_ae',
     >                              core_ae(2),core_ae(1))
      value = value.and.MA_alloc_get(mt_dbl,n1dgrid,'core_ps',
     >                              core_ps(2),core_ps(1))
      value = value.and.MA_alloc_get(mt_dbl,n1dgrid,'v_ps',
     >                              v_ps(2),v_ps(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis*n1dgrid,'prj_ps',
     >                              prj_ps(2),prj_ps(1))
      value = value.and.MA_alloc_get(mt_dbl,nbasis*n1dgrid,'prj_ps0',
     >                              prj_ps0(2),prj_ps0(1))
      if (.not.value)  call errquit('pawppv1:out of heap',1,MA_ERR)


      if (taskid.eq.MASTER) then
         read(11,*,ERR=9999,END=9999) (int_mb(nae(1)+j),
     >                                 dbl_mb(eig(1)+j),
     >                                 int_mb(nps(1)+j),
     >                                 int_mb(lps(1)+j),j=0,nbasis-1)
         read(11,*,ERR=9999,END=9999) ((dbl_mb(phi_ae(1)+i+j*n1dgrid),
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*,ERR=9999,END=9999) ((dbl_mb(dphi_ae(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*,ERR=9999,END=9999) ((dbl_mb(phi_ps(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*,ERR=9999,END=9999) ((dbl_mb(dphi_ps(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*,ERR=9999,END=9999) ((dbl_mb(prj_ps(1)+i+j*n1dgrid), 
     >                                  i=0,n1dgrid-1),
     >                                  j=0,nbasis-1)
         read(11,*,ERR=9999,END=9999) (dbl_mb(core_ae(1)+i),
     >                                 i=0,n1dgrid-1)
         read(11,*,ERR=9999,END=9999) (dbl_mb(core_ps(1)+i),
     >                                 i=0,n1dgrid-1)
         read(11,*,ERR=9999,END=9999) (dbl_mb(v_ps(1)+i),
     >                                 i=0,n1dgrid-1)
         read(11,*,ERR=9999,END=9999) sigma
         read(11,*,ERR=9999,END=9999) zion
         read(11,*,ERR=9999,END=9999) ((dbl_mb(prj_ps0(1)+i+j*n1dgrid),
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
      if(.not. MA_push_get(mt_int,lmax+1,'nmaxl',nmaxl(2),nmaxl(1)))
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
      if(.not.MA_pop_stack(nmaxl(2)))
     > call errquit('pawppv1:error popping stack',0,MA_ERR)


*     **** allocate Gijl,Sijl,Tijl,VSRCZ,Vpseuo ****
      l = nmax*nmax*(lmax+1)
      value=          MA_alloc_get(mt_dbl,5*l,'Gijl',Gijl(2),Gijl(1))
      !value=value.and.MA_alloc_get(mt_dbl,l,'Tijl',Tijl(2),Tijl(1))

*     **** allocate nprj,lprj,mprj ****
      value=value.and.MA_alloc_get(mt_int,nproj,'nprj',nprj(2),nprj(1))
      value=value.and.MA_alloc_get(mt_int,nproj,'lprj',lprj(2),lprj(1))
      value=value.and.MA_alloc_get(mt_int,nproj,'mprj',mprj(2),mprj(1))
      if (.not.value)  call errquit('pawppv1:out of heap',2,MA_ERR)


*    **** more temporary space ****
      value =           MA_alloc_get(mt_dbl,n1dgrid,'f',f(2),f(1))
      value = value.and.MA_alloc_get(mt_dbl,n1dgrid,'cs',cs(2),cs(1))
      value = value.and.MA_alloc_get(mt_dbl,n1dgrid,'sn',sn(2),sn(1))
      if (.not.value)call errquit('pawppv1:out of heap',0,MA_ERR)

*     **** allocate vl,vnl,vnlnrm G ****
      nsize = (ngrid(1)/2+1)*ngrid(2)*ngrid(3)
      value = MA_alloc_get(mt_dbl,nsize,'vl',vl(2),vl(1))
      value = value.and.MA_alloc_get(mt_dbl,nsize*(nproj),
     >                        'vnl',vnl(2), vnl(1))
      value = value.and.MA_alloc_get(mt_dbl,3*nsize,'G',G_hndl,G_indx)
      if (.not.value)call errquit('pawppv1:out of heap',0,MA_ERR)

*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp(nfft1,nfft2,nfft3,unita,unitg,dbl_mb(G_indx))
      filter = control_kbpp_filter()

      !**** allocate memory for rays ****
      nray = kbpp_calc_nray(nfft1,nfft2,nfft3,unita)

      value =           MA_alloc_get(mt_dbl,nray,
     >                  'G_ray',G_ray_hndl,G_ray_indx)
      value = value.and.MA_alloc_get(mt_dbl,2*nray,
     >                  'vl_ray',vl_ray_hndl,vl_ray_indx)
      value = value.and.MA_alloc_get(mt_dbl,2*nray*(nbasis),
     >                  'vnl_ray',vnl_ray_hndl,vnl_ray_indx)
      value = value.and.MA_alloc_get(mt_dbl,nray,
     >                  'tmp_ray',tmp_ray_hndl,tmp_ray_indx)
      if (.not.value)
     >   call errquit('pawppv1:out of heap memory',0,MA_ERR)

        call kbpp_generate_G_ray(nfft1,nfft2,nfft3,
     >                         unita,dbl_mb(G_ray_indx))

        zcore = zion - zv
        call integrate_pawppv1(version,rlocal,
     >                     n1dgrid,log_amesh,nbasis,nmax,lmax,zv,zcore,
     >                     int_mb(nps(1)),int_mb(lps(1)),
     >                                dbl_mb(v_ps(1)),
     >                                dbl_mb(prj_ps(1)),
     >                                dbl_mb(phi_ae(1)),
     >                                dbl_mb(dphi_ae(1)),
     >                                dbl_mb(phi_ps(1)),
     >                                dbl_mb(dphi_ps(1)),
     >                                dbl_mb(rgrid(1)),
     >                                dbl_mb(f(1)),
     >                                dbl_mb(cs(1)),
     >                                dbl_mb(sn(1)),
     >                      nfft1,nfft2,nfft3,nproj,
     >                                dbl_mb(G_indx),
     >                                dbl_mb(vl(1)),
     >                                dbl_mb(vnl(1)),
     >                                int_mb(nprj(1)),
     >                                int_mb(lprj(1)),
     >                                int_mb(mprj(1)),
     >                                dbl_mb(Gijl(1)),
     >                      nray,
     >                                dbl_mb(G_ray_indx),
     >                                dbl_mb(vl_ray_indx),
     >                                dbl_mb(vnl_ray_indx),
     >                                dbl_mb(tmp_ray_indx),
     >                                filter,
     >                      ierr)
        value = MA_free_heap(tmp_ray_hndl)
        value = value.and.MA_free_heap(vl_ray_hndl)
        value = value.and.MA_free_heap(vnl_ray_hndl)
        value = value.and.MA_free_heap(G_ray_hndl)
        if (.not.value)
     >   call errquit('pawppv1:Error freeing memory',0,MA_ERR)




c*     *** integrate core density ***
c      fourpi = 16.0d0*datan(1.0d0)
c      zcore = fourpi*log_integrate_def(0,
c     >                         dbl_mb(core_ae(1)),
c     >                         2,dbl_mb(rgrid(1)),log_amesh,n1dgrid)
c      write(*,*) "Zcore=",zcore,rmax,dbl_mb(rgrid(1)+n1dgrid-1)
c      zcore = fourpi*log_integrate_def(0,
c     >                         dbl_mb(core_ps(1)),
c     >                         2,dbl_mb(rgrid(1)),log_amesh,n1dgrid)
c      write(*,*) "Zcoreps=",zcore


      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(formatted_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      if (mprint) then
      write(*,*)
      write(*,*) "Generated formatted_filename: ",full_filename(1:l)
c      if (filter) write(*,*) "- filtering pseudopotential -"
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
         call dwrite(2,dbl_mb(Gijl(1)),5*nmax*nmax*(lmax+1))
         end if

         if (version.eq.4) call dwrite(2,rlocal,1)
         rcore = 0.0d0
         call dwrite(2,rcore,1)

         call dwrite(2,dbl_mb(vl(1)),nsize)
         if (nproj.gt.0) then
            call dwrite(2,dbl_mb(vnl(1)),nsize*nproj)
         end if

*        **** miscelaneous PAW energies ****
         call dwrite(2,core_kin_energy,1)

      call closefile(2)
      end if




      value = value.and.MA_free_heap(G_hndl)
      value = value.and.MA_free_heap(vnl(2))
      value = value.and.MA_free_heap(vl(2))
      value = value.and.MA_free_heap(sn(2))
      value = value.and.MA_free_heap(cs(2))
      value = value.and.MA_free_heap(f(2))

      value = value.and.MA_free_heap(Gijl(2))

      value = value.and.MA_free_heap(mprj(2))
      value = value.and.MA_free_heap(lprj(2))
      value = value.and.MA_free_heap(nprj(2))

      value = value.and.MA_free_heap(prj_ps0(2))
      value = value.and.MA_free_heap(prj_ps(2))
      value = value.and.MA_free_heap(v_ps(2))
      value = value.and.MA_free_heap(core_ps(2))
      value = value.and.MA_free_heap(core_ae(2))
      value = value.and.MA_free_heap(dphi_ps(2))
      value = value.and.MA_free_heap(phi_ps(2))
      value = value.and.MA_free_heap(dphi_ae(2))
      value = value.and.MA_free_heap(phi_ae(2))

      value = value.and.MA_free_heap(eig(2))
      value = value.and.MA_free_heap(lps(2))
      value = value.and.MA_free_heap(nps(2))
      value = value.and.MA_free_heap(nae(2))

      value = value.and.MA_free_heap(rgrid(2))
      if (.not.value)  call errquit('pawppv1:freeing heap',5,MA_ERR)

      pawppv1 = value
      return

 9999 call errquit('pawppv1:Error reading psp_filename',0,DISK_ERR)
      pawppv1 = value
      return

      END

      


