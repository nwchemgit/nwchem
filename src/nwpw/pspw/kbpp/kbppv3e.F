*
* $Id$
*

*     **************************************
*     *                                    *
*     *              kbppv3e               *
*     *                                    *
*     **************************************

      logical function kbppv3e(oprint_in,version,
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

      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmax1,locp1
      double precision rc(0:9),rlocal1
      integer nrho
      double precision drho
      integer rho_indx,vp_indx,wp_indx,up_indx,sc_r_indx,sc_k_indx
      integer rho_hndl,vp_hndl,wp_hndl,up_hndl,sc_r_hndl,sc_k_hndl

      integer          isemicore,ihasae,n_extra,n_expansion(0:9)
      logical          semicore,small_cell
      double precision rcore,core_charge

      integer f_indx,cs_indx,sn_indx
      integer n_prj_indx,l_prj_indx,m_prj_indx,b_prj_indx
      integer f_hndl,cs_hndl,sn_hndl
      integer n_prj_hndl,l_prj_hndl,m_prj_hndl,b_prj_hndl

      integer r3_matrix(2)

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
      logical  control_psp_semicore_small
      external control_psp_semicore_small
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
      ihasae = 0
  
*     ***** read in pseudopotential data ****
      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')

      read(11,'(I2)') ihasae
      read(11,'(A2)') atom
      read(11,*) zv,amass,lmax0,lmax1,locp1,rlocal1
      read(11,*) (rc(i),i=0,lmax0)
     
      read(11,*) n_extra
      read(11,*) (n_expansion(i),i=0,lmax0)
    
      read(11,*) nrho,drho
      read(11,'(A)') comment
      end if


      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,zv)
      call Parallel_Brdcst_values(MASTER,msglen,amass)
      call Parallel_Brdcst_values(MASTER,msglen,rlocal1)
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax0)
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax1)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp1)
      call Parallel_Brdcst_ivalues(MASTER,msglen,n_extra)

      msglen = lmax0+1
      call Parallel_Brdcst_values(MASTER,msglen,rc)
      call Parallel_Brdcst_ivalues(MASTER,msglen,n_expansion)
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
      value = BA_alloc_get(mt_dbl,nrho,
     >                        'rho',rho_hndl,rho_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(lmax+1),
     >                        'vp',vp_hndl, vp_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(lmax+1+n_extra),
     >                        'wp', wp_hndl, wp_indx)

      value = value.and.BA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)
      if (.not.value) 
     > call errquit('kbppv3e:out of heap memory',0,MA_ERR)

      if (taskid.eq.MASTER) then
         call read_vpwp_extra(11,nrho,lmax,lmax0,n_extra,
     >                        dbl_mb(rho_indx),
     >                        dbl_mb(vp_indx),
     >                        dbl_mb(wp_indx))
         call read_semicore(11,isemicore,rcore,nrho,dbl_mb(sc_r_indx))
         close(11)
      end if

      msglen = nrho
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(rho_indx))
      msglen = nrho*(lmax+1)
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(vp_indx))

      msglen = nrho*(lmax+1+n_extra)
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
      value = BA_alloc_get(mt_dbl,nrho,
     >                        'f',f_hndl,f_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho,
     >                        'cs',cs_hndl,cs_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho,
     >                        'sn',sn_hndl,sn_indx)
      if (.not.value) 
     > call errquit('kbppv3e:out of heap memory',0,MA_ERR)

*     **** allocate vl,vnl,vnlnrm G ****
      call D3dB_nfft3d(1,nsize)
      nprj = 0
      nmax = 1
      do l=0,lmax
         if (l.ne.locp) then
            nprj = nprj + n_expansion(l)*(2*l+1)
            if (n_expansion(l).gt.nmax) nmax = n_expansion(l)
         end if
      end do
      !lmmax = (lmax+1)**2 - (2*locp+1)
      value = BA_alloc_get(mt_dbl,nsize,
     >                        'vl',vl_hndl,vl_indx)
      value = value.and.BA_alloc_get(mt_dbl,nsize*(nprj),
     >                        'vnl',vnl_hndl, vnl_indx)
      value = value.and.BA_alloc_get(mt_dbl,nmax*nmax*(lmax+1),
     >                        'vnlnrm', vnlnrm_hndl, vnlnrm_indx)
      value = value.and.BA_alloc_get(mt_dbl,4*nsize,
     >                        'sc_k',sc_k_hndl,sc_k_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'n_prj', n_prj_hndl, n_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'l_prj', l_prj_hndl, l_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'm_prj', m_prj_hndl, m_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'b_prj', b_prj_hndl, b_prj_indx)
      value = value.and.BA_alloc_get(mt_dcpl,nsize,'tmp',tmp(2),tmp(1))
      if (.not.value) 
     > call errquit('kbppv3e:out of heap memory',0,MA_ERR)

      filter = control_kbpp_filter()
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)

      if (control_kbpp_ray()) then
        !**** allocate memory for rays ****
        nray = kbpp_calc_nray(nfft1,nfft2,nfft3,unita)

        value =           BA_alloc_get(mt_dbl,nray,
     >                  'G_ray',G_ray_hndl,G_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,2*nray,
     >                  'vl_ray',vl_ray_hndl,vl_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,2*nray*(lmax+1+n_extra),
     >                  'vnl_ray',vnl_ray_hndl,vnl_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,2*nray*2,
     >              'rho_sc_k_ray',rho_sc_k_ray_hndl,rho_sc_k_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,nray,
     >                  'tmp_ray',tmp_ray_hndl,tmp_ray_indx)
        if (.not.value)
     >   call errquit('kbppv3e:out of heap memory',0,MA_ERR)

        call kbpp_generate_G_ray(nfft1,nfft2,nfft3,
     >                         unita,
     >                         dbl_mb(G_ray_indx))
        call integrate_kbppv3e_new(version,rlocal,
     >                      nrho,drho,lmax,locp,nmax,
     >                      n_extra,n_expansion,zv,
     >                                dbl_mb(vp_indx),
     >                                dbl_mb(wp_indx),
     >                                dbl_mb(rho_indx),
     >                                dbl_mb(f_indx),
     >                                dbl_mb(cs_indx),
     >                                dbl_mb(sn_indx),
     >                      nsize,nprj,
     >                                dbl_mb(G_indx(1)),
     >                                dbl_mb(vl_indx),
     >                                dbl_mb(vnl_indx),
     >                                int_mb(n_prj_indx),
     >                                int_mb(l_prj_indx),
     >                                int_mb(m_prj_indx),
     >                                int_mb(b_prj_indx),
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
        value = BA_free_heap(tmp_ray_hndl)
        value = value.and.BA_free_heap(rho_sc_k_ray_hndl)
        value = value.and.BA_free_heap(vnl_ray_hndl)
        value = value.and.BA_free_heap(vl_ray_hndl)
        value = value.and.BA_free_heap(G_ray_hndl)
        if (.not.value)
     >   call errquit('kbppv3e:Error freeing memory',0,MA_ERR)

      else
        call integrate_kbppv3e(version,rlocal,
     >                      nrho,drho,lmax,locp,nmax,
     >                      n_extra,n_expansion,zv,
     >                                dbl_mb(vp_indx),
     >                                dbl_mb(wp_indx),
     >                                dbl_mb(rho_indx),
     >                                dbl_mb(f_indx),
     >                                dbl_mb(cs_indx),
     >                                dbl_mb(sn_indx),
     >                      nsize,nprj,
     >                                dbl_mb(G_indx(1)),
     >                                dbl_mb(vl_indx),
     >                                dbl_mb(vnl_indx),
     >                                int_mb(n_prj_indx),
     >                                int_mb(l_prj_indx),
     >                                int_mb(m_prj_indx),
     >                                int_mb(b_prj_indx),
     >                                dbl_mb(vnlnrm_indx),
     >                                semicore,
     >                                dbl_mb(sc_r_indx),
     >                                dbl_mb(sc_k_indx),
     >                      ierr)
      end if

      !**** small unitcell semicore - convert rho(G) to chi(G) ****
      small_cell = control_psp_semicore_small().and.semicore
      if (small_cell) then
         call semicore_rho_to_chi(dbl_mb(sc_k_indx))
      end if

*     **** compute the <1/r3> matrix elements ****


      if ((taskid.eq.MASTER).and.(oprint)) then
      write(luout,*) "     ********************************************"
      write(luout,*) "     *                                          *"
      write(luout,*) "     *     KBPPV3d - Pseudopotential Formatter  *"
      write(luout,*) "     *                                          *"
      write(luout,*) "     *      version last updated 4/15/99        *"
      write(luout,*) "     *                                          *"
      write(luout,*) "     * This code was written by Eric J. Bylaska *"
      write(luout,*) "     * and was based upon algorithms and code   *"
      write(luout,*) "     * developed by the group of                *"
      write(luout,*) "     * Prof. John H. Weare                      *"
      write(luout,*) "     ********************************************"
      call nwpw_message(1)
      write(luout,*)
      write(luout,*) "Pseudpotential Data"
      write(luout,*) "-------------------"
      write(luout,*) "  atom     :",atom
      write(luout,*) "  charge   :",zv
      write(luout,*) "  mass no. :",amass
      write(luout,*) "  highest angular component      :",lmax0
      write(luout,*) "  highest angular component used :",lmax
      write(luout,*) "  local potential used           :",locp
      write(luout,*) "  number of projectors           :",nprj
      if (version.eq.4) 
     >   write(luout,*) "  aperiodic cutoff radius        :",rlocal
      write(luout,111) "  cutoffs: ",(rc(i), i=0,lmax)
      if (semicore) then
        write(luout,*)
        write(luout,115) "  semi-core charge included, rcore:",rcore
        do i=1,nrho
           dbl_mb(f_indx+i-1) = dbl_mb(sc_r_indx+i-1)
     >                        * dbl_mb(rho_indx+i-1)**2
        end do
        core_charge=16.0d0*datan(1.0d0)*SIMP(nrho,dbl_mb(f_indx),drho)
        write(luout,115) "  semi-core charge                :",
     >                   core_charge,dbl_mb(sc_k_indx)
         do i=1,nrho
           dbl_mb(f_indx+i-1) = dbl_mb(sc_r_indx+i-1+nrho)
     >                        * dbl_mb(rho_indx+i-1)**2
         end do
         core_charge=16.0d0*datan(1.0d0)*SIMP(nrho,dbl_mb(f_indx),drho)
         write(luout,115) "  Semi-core charge gradient       :",
     >                core_charge

      end if
      write(luout,*)
      write(luout,*) "Simulation Cell"
      write(luout,*) "---------------"
      if (version.eq.3) write(luout,112) "  boundry: periodic"
      if (version.eq.4) write(luout,112) "  boundry: aperiodic"
      write(luout,113) "  ngrid  :",ngrid
      write(luout,114) "  unita  :",unita(1,1),unita(2,1),unita(3,1)
      write(luout,114) "          ",unita(1,2),unita(2,2),unita(3,2)
      write(luout,114) "          ",unita(1,3),unita(2,3),unita(3,3)
      write(luout,*)
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
      write(luout,*)
      write(luout,*) "Generated formatted_filename: ",full_filename(1:l)
      if (filter) write(luout,*) "- filtering pseudopotential -"
      if (small_cell) 
     >  write(luout,*) "- small unit cell semicore correction applied -"
      !write(luout,*)
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
     
         call iwrite(2,nmax,1)
         call dwrite(2,rc,lmax+1)

         call iwrite(2,nprj,1)
         if (nprj.gt.0) then
         call iwrite(2,int_mb(n_prj_indx),nprj)
         call iwrite(2,int_mb(l_prj_indx),nprj)
         call iwrite(2,int_mb(m_prj_indx),nprj)
         call iwrite(2,int_mb(b_prj_indx),nprj)
         call dwrite(2,dbl_mb(vnlnrm_indx),nmax*nmax*(lmax+1))
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
      value = BA_free_heap(rho_hndl)
      value = value.and.BA_free_heap(vp_hndl)
      value = value.and.BA_free_heap(wp_hndl)
      value = value.and.BA_free_heap(sc_r_hndl)
      value = value.and.BA_free_heap(sc_k_hndl)
      value = value.and.BA_free_heap(f_hndl)
      value = value.and.BA_free_heap(cs_hndl)
      value = value.and.BA_free_heap(sn_hndl)

      value = value.and.BA_free_heap(vl_hndl)
      value = value.and.BA_free_heap(vnl_hndl)
      value = value.and.BA_free_heap(vnlnrm_hndl)
      value = value.and.BA_free_heap(n_prj_hndl)
      value = value.and.BA_free_heap(l_prj_hndl)
      value = value.and.BA_free_heap(m_prj_hndl)
      value = value.and.BA_free_heap(b_prj_hndl)
      value = value.and.BA_free_heap(tmp(2))
      if (.not.value) 
     > call errquit('kbppv3:Error freeing memory',0,MA_ERR)

      
      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)

      kbppv3e = value
      return

 9999 call errquit('kbppv3e:Error reading psp_filename',0, DISK_ERR)
      kbppv3e = value
      return

      END


