*
* $Id$
*


*     **************************************
*     *                                    *
*     *           kbppv3e_stress           *
*     *                                    *
*     **************************************

      logical function kbppv3e_stress(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,locp,lmax,rlocal)
      implicit none
#include "errquit.fh"
#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "util.fh"

      logical          oprint_in
      integer          version
      character*50     psp_filename,formatted_filename
      integer          ngrid(3)
      double precision unita(3,3)
      integer locp,lmax
      real*8  rlocal
      character*255 full_filename

      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmax1,locp1,nmax,nprj
      double precision rc(0:9),rlocal1
      integer nrho
      double precision drho
      integer rho_indx,vp_indx,wp_indx,sc_r_indx,sc_k_indx
      integer rho_hndl,vp_hndl,wp_hndl,sc_r_hndl,sc_k_hndl

      integer          isemicore,ihasae,n_extra,n_expansion(0:9)
      logical          semicore
      double precision rcore

      integer f_indx,cs_indx,sn_indx
      integer f_hndl,cs_hndl,sn_hndl

*     ***** ngrid data *****
      integer vl_indx,vnl_indx,tmp(2)
      integer vl_hndl,vnl_hndl

*     **** ray data ****
      integer nray,G_ray_hndl,tmp_ray_hndl
      integer rho_sc_k_ray_hndl,dvnl_ray_hndl,dvl_ray_hndl
      integer G_ray_indx,tmp_ray_indx
      integer rho_sc_k_ray_indx,dvnl_ray_indx,dvl_ray_indx

*     **** other variables ****
      logical value,mprint,hprint,oprint,filter
      double precision unitg(3,3)
      integer nsize,i,l,ierr
      integer nfft1,nfft2,nfft3

*     **** external functions ****
      logical  control_print,control_kbpp_ray,control_kbpp_filter
      external control_print,control_kbpp_ray,control_kbpp_filter
      integer  kbpp_calc_nray
      external kbpp_calc_nray
      integer  G_indx
      external G_indx


c      call Parallel_init()
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
     >  call errquit('kbppv3e_stress: out of heap memory',0,MA_ERR)

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
     >  call errquit('kbppv3e_stress: out of heap memory',0,MA_ERR)

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
      value = value.and.BA_alloc_get(mt_dbl,nsize*(nprj)*3,
     >                        'vnl',vnl_hndl, vnl_indx)
      value = value.and.BA_alloc_get(mt_dbl,4*nsize,
     >                        'sc_k',sc_k_hndl,sc_k_indx)
      value = value.and.BA_alloc_get(mt_dcpl,nsize,'tmp',tmp(2),tmp(1))
      if (.not.value)
     >  call errquit('kbppv3e_stress: out of heap memory',0,MA_ERR)


*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      filter = control_kbpp_filter()

      if (control_kbpp_ray()) then
        !**** allocate memory for rays ****
        nray = kbpp_calc_nray(nfft1,nfft2,nfft3,unita)

        value =           BA_alloc_get(mt_dbl,nray,
     >                  'G_ray',G_ray_hndl,G_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,2*nray,
     >                  'dvl_ray',dvl_ray_hndl,dvl_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,4*nray*(lmax+1+n_extra),
     >                  'dvnl_ray',dvnl_ray_hndl,dvnl_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,2*nray,
     >              'rho_sc_k_ray',rho_sc_k_ray_hndl,rho_sc_k_ray_indx)
        value = value.and.BA_alloc_get(mt_dbl,nray,
     >                  'tmp_ray',tmp_ray_hndl,tmp_ray_indx)
        if (.not.value)
     >   call errquit('kbppv3e_stress:out of heap memory',0,MA_ERR)

        call kbpp_generate_G_ray(nfft1,nfft2,nfft3,
     >                         unita,
     >                         dbl_mb(G_ray_indx))

        call integrate_e_stress_new(version,rlocal,
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
     >                                semicore,
     >                                dbl_mb(sc_r_indx),
     >                                dbl_mb(sc_k_indx),
     >                      nray,
     >                                dbl_mb(G_ray_indx),
     >                                dbl_mb(dvl_ray_indx),
     >                                dbl_mb(dvnl_ray_indx),
     >                                dbl_mb(rho_sc_k_ray_indx),
     >                                dbl_mb(tmp_ray_indx),
     >                                filter,
     >                      ierr)
        value = BA_free_heap(tmp_ray_hndl)
        value = value.and.BA_free_heap(rho_sc_k_ray_hndl)
        value = value.and.BA_free_heap(dvnl_ray_hndl)
        value = value.and.BA_free_heap(dvl_ray_hndl)
        value = value.and.BA_free_heap(G_ray_hndl)
        if (.not.value)
     >   call errquit('kbppv3_stress:Error freeing memory',0,MA_ERR)

      else

      call integrate_e_stress(version,rlocal,
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
     >                                semicore,
     >                                dbl_mb(sc_r_indx),
     >                                dbl_mb(sc_k_indx),
     >                      ierr)
      end if


      if ((taskid.eq.MASTER).and.(oprint)) then
      write(*,*) "     ***********************************************"
      write(*,*) "     *                                             *"
      write(*,*) "     *  KBPPV3e-Stress - Pseudopotential Formatter  *"
      write(*,*) "     *                                             *"
      write(*,*) "     *      version last updated 4/15/99           *"
      write(*,*) "     *                                             *"
      write(*,*) "     ***********************************************"
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
      write(*,*) "  number of projectors           :",nprj

      if (semicore) then
        write(*,*)
        write(*,115) "  semi-core stress derivative added"
      end if
      write(*,*)
      write(*,*) "Simulation Cell"
      write(*,*) "---------------"
      if (version.eq.3) write(*,112) "  boundry: periodic"
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
      write(*,*) "Outputing formatted_stress_filename: ",
     >            full_filename(1:l)
      if (filter) write(*,*) "- filtering pseudopotential -"
      !write(*,*)
      end if
      
         call openfile(2,full_filename,l,'w',l)     
         call iwrite(2,version,1)
         call iwrite(2,ngrid,3)
         call dwrite(2,unita,9)
      endif

      call D3dB_t_write(1,2,dbl_mb(vl_indx),dcpl_mb(tmp(1)),0)
      do l=1,nprj
      do i=1,3
         call D3dB_t_write(1,2,
     >                     dbl_mb(vnl_indx+((i-1)+(l-1)*3)*nsize),
     >                     dcpl_mb(tmp(1)),0)
      end do
      end do
      if (semicore) then
         call D3dB_t_write(1,2,dbl_mb(sc_k_indx),dcpl_mb(tmp(1)),0)
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
      value = value.and.BA_free_heap(tmp(2))
      if (.not.value)
     >  call errquit('kbppv3e_stress: error freeing heap',0,MA_ERR)

      
      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)

      kbppv3e_stress = value
      return

 9999 kbppv3e_stress = .false.
      call errquit('Error reading psp_filename',0, DISK_ERR)
      END




