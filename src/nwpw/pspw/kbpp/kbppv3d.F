*
* $Id: kbppv3d.F,v 1.1 2008-03-17 17:53:58 bylaska Exp $
*

*     **************************************
*     *                                    *
*     *              kbppv3d               *
*     *                                    *
*     **************************************

      logical function kbppv3d(oprint_in,version,
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

      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmmax,lmax1,locp1
      double precision rc(0:9),rlocal1
      integer nrho
      double precision drho
      integer rho_indx,vp_indx,wp_indx,sc_r_indx,sc_k_indx
      integer rho_hndl,vp_hndl,wp_hndl,sc_r_hndl,sc_k_hndl

      integer          isemicore
      logical          semicore
      double precision rcore,core_charge

      integer f_indx,cs_indx,sn_indx
      integer n_prj_indx,l_prj_indx,m_prj_indx
      integer f_hndl,cs_hndl,sn_hndl
      integer n_prj_hndl,l_prj_hndl,m_prj_hndl

*     ***** ngrid data *****
      integer vl_indx,vnl_indx,vnlnrm_indx,tmp(2)
      integer vl_hndl,vnl_hndl,vnlnrm_hndl

*     **** ray data ****
      integer nray,G_ray_hndl,tmp_ray_hndl
      integer rho_sc_k_ray_hndl,vnl_ray_hndl,vl_ray_hndl
      integer G_ray_indx,tmp_ray_indx
      integer rho_sc_k_ray_indx,vnl_ray_indx,vl_ray_indx


*     **** other variables ****
      logical value,mprint,hprint,oprint,filter
      integer nsize,i,l,ierr,psp_type
      integer nfft1,nfft2,nfft3
      integer nmax,nprj
      character*255 full_filename

*     **** external functions ****
      logical  control_print,control_kbpp_ray,control_kbpp_filter
      external control_print,control_kbpp_ray,control_kbpp_filter
      double precision simp
      external         simp
      integer  kbpp_calc_nray
      external kbpp_calc_nray
      integer  G_indx
      external G_indx


c      call Parallel_Init()
      call Parallel_taskid(taskid)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)
      oprint = (oprint_in.or.hprint)

      value = .false.
  
*     ***** read in pseudopotential data ****
      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')
      read(11,'(A2)',ERR=9999,END=9999) atom
      read(11,*,ERR=9999,END=9999) zv,amass,lmax0,lmax1,locp1,rlocal1
      read(11,*,ERR=9999,END=9999) (rc(i),i=0,lmax0)
      read(11,*,ERR=9999,END=9999) nrho,drho
      read(11,'(A)',ERR=9999,END=9999) comment
      end if

      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,zv)
      call Parallel_Brdcst_values(MASTER,msglen,amass)
      call Parallel_Brdcst_values(MASTER,msglen,rlocal1)
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax0)
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax1)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp1)

      msglen = lmax0+1
      call Parallel_Brdcst_values(MASTER,msglen,rc)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,nrho)
      call Parallel_Brdcst_values(MASTER,msglen,drho)


*     **** set the maximum angular momentum ****
      if (lmax.eq.-1)    lmax = lmax1
      if (lmax.gt.lmax0) lmax = lmax0
      if (lmax.lt.0)     lmax = lmax0

*     **** set the local potential ****
      if (locp.eq.-1)   locp = locp1
      if (locp.gt.lmax) locp = lmax
      if (locp.lt.0)    locp = lmax

*     **** set the local potential ****
      if (rlocal.eq.1.0d0)   rlocal = rlocal1
      

*     **** allocate rho, vp, and wp ****
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'rho',rho_hndl,rho_indx)
      value = value.and.MA_alloc_get(mt_dbl,nrho*(lmax+1),
     >                        'vp',vp_hndl, vp_indx)
      value = value.and.MA_alloc_get(mt_dbl,nrho*(lmax+1),
     >                        'wp', wp_hndl, wp_indx)
      value = value.and.MA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)
      if (.not.value) 
     > call errquit('kbppv3d:out of heap memory',0,MA_ERR)

      if (taskid.eq.MASTER) then
      call read_vpwp(11,nrho,lmax,dbl_mb(rho_indx),
     >                         dbl_mb(vp_indx),
     >                         dbl_mb(wp_indx))
      call read_semicore(11,isemicore,rcore,nrho,dbl_mb(sc_r_indx))
      close(11)
      end if

      msglen = nrho
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(rho_indx))
      msglen = nrho*(lmax+1)
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(vp_indx))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(wp_indx))

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,isemicore)
      semicore = (isemicore.eq.1)
      if (semicore) then
      msglen = 2*nrho
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(sc_r_indx))
      else
         rcore = 0.0d0
      end if


*    **** more temporary space ****
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'f',f_hndl,f_indx)
      value = value.and.MA_alloc_get(mt_dbl,nrho,
     >                        'cs',cs_hndl,cs_indx)
      value = value.and.MA_alloc_get(mt_dbl,nrho,
     >                        'sn',sn_hndl,sn_indx)
      if (.not.value) 
     > call errquit('kbppv3d:out of heap memory',0,MA_ERR)

*     **** allocate vl,vnl,vnlnrm G ****
      call D3dB_nfft3d(1,nsize)
      lmmax = (lmax+1)**2 - (2*locp+1)
      nprj  = lmmax
      value = MA_alloc_get(mt_dbl,nsize,
     >                        'vl',vl_hndl,vl_indx)
      value = value.and.MA_alloc_get(mt_dbl,nsize*(lmmax),
     >                        'vnl',vnl_hndl, vnl_indx)
      value = value.and.MA_alloc_get(mt_dbl,(lmax+1),
     >                        'vnlnrm', vnlnrm_hndl, vnlnrm_indx)
      value = value.and.MA_alloc_get(mt_dbl,4*nsize,
     >                        'sc_k',sc_k_hndl,sc_k_indx)
      value = value.and.MA_alloc_get(mt_int,nprj,
     >                        'n_prj', n_prj_hndl, n_prj_indx)
      value = value.and.MA_alloc_get(mt_int,nprj,
     >                        'l_prj', l_prj_hndl, l_prj_indx)
      value = value.and.MA_alloc_get(mt_int,nprj,
     >                        'm_prj', m_prj_hndl, m_prj_indx)
      value = value.and.MA_alloc_get(mt_dcpl,nsize,'tmp',tmp(2),tmp(1))
      if (.not.value) 
     > call errquit('kbppv3d:out of heap memory',0,MA_ERR)

      filter = control_kbpp_filter()
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)

      if (control_kbpp_ray()) then
        !**** allocate memory for rays ****
        nray = kbpp_calc_nray(nfft1,nfft2,nfft3,unita)

        value =           MA_alloc_get(mt_dbl,nray,
     >                  'G_ray',G_ray_hndl,G_ray_indx)
        value = value.and.MA_alloc_get(mt_dbl,2*nray,
     >                  'vl_ray',vl_ray_hndl,vl_ray_indx)
        value = value.and.MA_alloc_get(mt_dbl,2*nray*(lmax+1),
     >                  'vnl_ray',vnl_ray_hndl,vnl_ray_indx)
        value = value.and.MA_alloc_get(mt_dbl,2*nray*2,
     >              'rho_sc_k_ray',rho_sc_k_ray_hndl,rho_sc_k_ray_indx)
        value = value.and.MA_alloc_get(mt_dbl,nray,
     >                  'tmp_ray',tmp_ray_hndl,tmp_ray_indx)
        if (.not.value)
     >   call errquit('kbppv3d:out of heap memory',0,MA_ERR)

        call kbpp_generate_G_ray(nfft1,nfft2,nfft3,
     >                         unita,
     >                         dbl_mb(G_ray_indx))
        call integrate_kbppv3d_new(version,rlocal,
     >                      nrho,drho,lmax,locp,zv,
     >                                dbl_mb(vp_indx),
     >                                dbl_mb(wp_indx),
     >                                dbl_mb(rho_indx),
     >                                dbl_mb(f_indx),
     >                                dbl_mb(cs_indx),
     >                                dbl_mb(sn_indx),
     >                      nsize,lmmax,
     >                                dbl_mb(G_indx(1)),
     >                                dbl_mb(vl_indx),
     >                                dbl_mb(vnl_indx),
     >                                int_mb(n_prj_indx),
     >                                int_mb(l_prj_indx),
     >                                int_mb(m_prj_indx),
     >                                dbl_mb(vnlnrm_indx),
     >                                semicore,
     >                                dbl_mb(sc_r_indx),
     >                                dbl_mb(sc_k_indx),
     >                      nray,
     >                                dbl_mb(G_ray_indx),
     >                                dbl_mb(vl_ray_indx),
     >                                dbl_mb(vnl_ray_indx),
     >                                dbl_mb(rho_sc_k_ray_indx),
     >                                dbl_mb(tmp_ray_indx),
     >                                filter,
     >                      ierr)
        value = MA_free_heap(tmp_ray_hndl)
        value = value.and.MA_free_heap(rho_sc_k_ray_hndl)
        value = value.and.MA_free_heap(vnl_ray_hndl)
        value = value.and.MA_free_heap(vl_ray_hndl)
        value = value.and.MA_free_heap(G_ray_hndl)
        if (.not.value)
     >   call errquit('kbppv3d:Error freeing memory',0,MA_ERR)

      else
        call integrate_kbppv3d(version,rlocal,
     >                      nrho,drho,lmax,locp,zv,
     >                                dbl_mb(vp_indx),
     >                                dbl_mb(wp_indx),
     >                                dbl_mb(rho_indx),
     >                                dbl_mb(f_indx),
     >                                dbl_mb(cs_indx),
     >                                dbl_mb(sn_indx),
     >                      nsize,lmmax,
     >                                dbl_mb(G_indx(1)),
     >                                dbl_mb(vl_indx),
     >                                dbl_mb(vnl_indx),
     >                                int_mb(n_prj_indx),
     >                                int_mb(l_prj_indx),
     >                                int_mb(m_prj_indx),
     >                                dbl_mb(vnlnrm_indx),
     >                                semicore,
     >                                dbl_mb(sc_r_indx),
     >                                dbl_mb(sc_k_indx),
     >                      ierr)
      end if



      if ((taskid.eq.MASTER).and.(oprint)) then
      write(*,*) "     ********************************************"
      write(*,*) "     *                                          *"
      write(*,*) "     *     KBPPV3d - Pseudopotential Formatter  *"
      write(*,*) "     *                                          *"
      write(*,*) "     *      version last updated 4/15/99        *"
      write(*,*) "     *                                          *"
      write(*,*) "     * This code was written by Eric J. Bylaska *"
      write(*,*) "     * and was based upon algorithms and code   *"
      write(*,*) "     * developed by the group of                *"
      write(*,*) "     * Prof. John H. Weare                      *"
      write(*,*) "     ********************************************"
      call nwpw_message(1)
      write(*,*)
      write(*,*) "Pseudpotential Data"
      write(*,*) "-------------------"
      write(*,*) "  atom     :",atom
      write(*,*) "  charge   :",zv
      write(*,*) "  mass no. :",amass
      write(*,*) "  highest angular component      :",lmax0
      write(*,*) "  highest angular component used :",lmax
      write(*,*) "  local potential used           :",locp
      if (version.eq.4) 
     >   write(*,*) "  aperiodic cutoff radius        :",rlocal
      write(*,111) "  cutoffs: ",(rc(i), i=0,lmax)
      if (semicore) then
        write(*,*)
        write(*,115) "  semi-core charge included, rcore:",rcore
        do i=1,nrho
           dbl_mb(f_indx+i-1) = dbl_mb(sc_r_indx+i-1)
     >                        * dbl_mb(rho_indx+i-1)**2
        end do
        core_charge=16.0d0*datan(1.0d0)*SIMP(nrho,dbl_mb(f_indx),drho)
        write(*,115) "  semi-core charge                :",core_charge,
     >                dbl_mb(sc_k_indx)
         do i=1,nrho
           dbl_mb(f_indx+i-1) = dbl_mb(sc_r_indx+i-1+nrho)
     >                        * dbl_mb(rho_indx+i-1)**2
         end do
         core_charge=16.0d0*datan(1.0d0)*SIMP(nrho,dbl_mb(f_indx),drho)
         write(*,115) "  Semi-core charge gradient       :",
     >                core_charge

      end if
      write(*,*)
      write(*,*) "Simulation Cell"
      write(*,*) "---------------"
      if (version.eq.3) write(*,112) "  boundry: periodic"
      if (version.eq.4) write(*,112) "  boundry: aperiodic"
      write(*,113) "  ngrid  :",ngrid
      write(*,114) "  unita  :",unita(1,1),unita(2,1),unita(3,1)
      write(*,114) "          ",unita(1,2),unita(2,2),unita(3,2)
      write(*,114) "          ",unita(1,3),unita(2,3),unita(3,3)
      write(*,*)
  111 format(a,10f10.3)
  112 format(a)
  113 format(a,3I4)
  114 format(a,3F10.3)
  115 format(a,2E14.6)
      end if

      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(formatted_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      if (mprint) then
      write(*,*)
      write(*,*) "Generated formatted_filename: ",full_filename(1:l)
      if (filter) write(*,*) "- filtering pseudopotential -"
      !write(*,*)
      end if
      call openfile(2,full_filename,l,'w',l)     

         call cwrite(2,comment,80)
         psp_type = 0
         call iwrite(2,psp_type,1)
         call iwrite(2,version,1)
         call iwrite(2,ngrid,3)
         call dwrite(2,unita,9)
         call cwrite(2,atom,2)
         call dwrite(2,amass,1)
         call dwrite(2,zv,1)
         call iwrite(2,lmax,1)
         call iwrite(2,locp,1)
     
         nmax = 1
         call iwrite(2,nmax,1)
         call dwrite(2,rc,lmax+1)

         call iwrite(2,nprj,1)
         if (nprj.gt.0) then
         call iwrite(2,int_mb(n_prj_indx),nprj)
         call iwrite(2,int_mb(l_prj_indx),nprj)
         call iwrite(2,int_mb(m_prj_indx),nprj)
         call dwrite(2,dbl_mb(vnlnrm_indx),(lmax+1))
         end if

         if (version.eq.4) call dwrite(2,rlocal,1)
         call dwrite(2,rcore,1)
      endif

      call D3dB_t_write(1,2,dbl_mb(vl_indx),dcpl_mb(tmp(1)),0)
      do i=1,nprj
         call D3dB_t_write(1,2,
     >                     dbl_mb(vnl_indx+(i-1)*nsize),
     >                     dcpl_mb(tmp(1)),0)

      end do

      if (semicore) then
      call D3dB_t_write(1,2,dbl_mb(sc_k_indx),        dcpl_mb(tmp(1)),0)
      call D3dB_t_write(1,2,dbl_mb(sc_k_indx+nsize),  dcpl_mb(tmp(1)),0)
      call D3dB_t_write(1,2,dbl_mb(sc_k_indx+2*nsize),dcpl_mb(tmp(1)),0)
      call D3dB_t_write(1,2,dbl_mb(sc_k_indx+3*nsize),dcpl_mb(tmp(1)),0)
      end if
      

      if (taskid.eq.MASTER) then
      call closefile(2)
      end if
     

*     **** free heap space ****
      value = MA_free_heap(rho_hndl)
      value = value.and.MA_free_heap(vp_hndl)
      value = value.and.MA_free_heap(wp_hndl)
      value = value.and.MA_free_heap(sc_r_hndl)
      value = value.and.MA_free_heap(sc_k_hndl)
      value = value.and.MA_free_heap(f_hndl)
      value = value.and.MA_free_heap(cs_hndl)
      value = value.and.MA_free_heap(sn_hndl)

      value = value.and.MA_free_heap(vl_hndl)
      value = value.and.MA_free_heap(vnl_hndl)
      value = value.and.MA_free_heap(vnlnrm_hndl)
      value = value.and.MA_free_heap(n_prj_hndl)
      value = value.and.MA_free_heap(l_prj_hndl)
      value = value.and.MA_free_heap(m_prj_hndl)
      value = value.and.MA_free_heap(tmp(2))
      if (.not.value) 
     > call errquit('kbppv3:Error freeing memory',0,MA_ERR)

      
      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)

      kbppv3d = value
      return

 9999 call errquit('kbppv3d:Error reading psp_filename',0, DISK_ERR)
      kbppv3d = value
      return

      END




