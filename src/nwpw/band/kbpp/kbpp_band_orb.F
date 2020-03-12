*
* $Id$
*

*     **************************************
*     *                                    *
*     *           kbpp_band_orb            *
*     *                                    *
*     **************************************

      logical function kbpp_band_orb(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,locp,lmax,
     >                  nbrillioun,kvectors,
     >                  rcut,lmbda)
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
      integer nbrillioun
      real*8  kvectors(3,*)
      real*8 rcut,lmbda

*     **** local variables ****
      character*255 full_filename

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
      integer rho_indx,vp_indx,wp_indx,sc_r_indx
      integer rho_hndl,vp_hndl,wp_hndl,sc_r_hndl

      integer          isemicore
      logical          semicore
      double precision rcore,core_charge

      integer f_indx,cs_indx,sn_indx
      integer f_hndl,cs_hndl,sn_hndl

*     ***** ngrid data *****
      integer borbs_indx,G_indx
      integer borbs_hndl,G_hndl


*     **** other variables ****
      logical value,mprint,hprint,oprint,filter,kray
      double precision unitg(3,3)
      integer nsize,i,l,ierr,nb,psp_type
      integer nfft1,nfft2,nfft3
      integer nmax

*     **** external functions ****
      logical  control_print
      external control_print

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
      read(11,*) atom
      read(11,*) zv,amass,lmax0,lmax1,locp1,rlocal1
      read(11,*) (rc(i),i=0,lmax0)
      read(11,*) nrho,drho
      read(11,'(A)') comment
      end if

      psp_type=0
      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,zv)
      call Parallel_Brdcst_values(MASTER,msglen,amass)
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax0)
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax1)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp1)

      msglen = lmax0+1
      call Parallel_Brdcst_values(MASTER,msglen,rc)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,nrho)
      call Parallel_Brdcst_values(MASTER,msglen,drho)

*     **** set the maximum angular momentum ****
      lmax = lmax0

*     **** set the local potential ****
      locp = lmax

*     **** allocate rho, vp, and wp ****
      value = BA_alloc_get(mt_dbl,nrho,
     >                        'rho',rho_hndl,rho_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(lmax+1),
     >                        'vp',vp_hndl, vp_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(lmax+1),
     >                        'wp', wp_hndl, wp_indx)
      value = value.and.BA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)
      if (.not.value)
     > call errquit("kbpp_band_orb: out of stack memory",0,MA_ERR)

      if (taskid.eq.MASTER) then
      call read_vpwp_band(11,nrho,lmax,dbl_mb(rho_indx),
     >                         dbl_mb(vp_indx),
     >                         dbl_mb(wp_indx))
c      call read_semicore_band(11,isemicore,rcore,nrho,dbl_mb(sc_r_indx))
      close(11)
      end if

      msglen = nrho
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(rho_indx))
      msglen = nrho*(lmax+1)
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(vp_indx))
      call Parallel_Brdcst_values(MASTER,msglen,dbl_mb(wp_indx))

      msglen = 1
      isemicore = 0
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
      value = BA_alloc_get(mt_dbl,nrho,
     >                        'cs',cs_hndl,cs_indx)
      value = BA_alloc_get(mt_dbl,nrho,
     >                        'sn',sn_hndl,sn_indx)
      if (.not.value)
     > call errquit("kbpp_band_orb: out of stack memory",0,MA_ERR)

*     **** allocate vl,vnl,vnlnrm G ****
      lmmax = (lmax+1)**2 - (2*locp+1)
      nsize = ngrid(1)*ngrid(2)*ngrid(3)

      value = value.and.BA_alloc_get(mt_dbl,nsize*(lmmax),
     >                        'borbs',borbs_hndl, borbs_indx)
      value = value.and.BA_alloc_get(mt_dbl,nsize*(3),
     >                        'G',G_hndl, G_indx)

      if (.not.value)
     > call errquit("kbpp_band_orb: out of stack memory",0,MA_ERR)

      call damp_kbpp_band_orb(nrho,lmax,rcut,lmbda,
     >                   dbl_mb(rho_indx),dbl_mb(wp_indx))




*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp_band(nfft1,nfft2,nfft3,unita,unitg,dbl_mb(G_indx))
     

      if ((taskid.eq.MASTER).and.(oprint)) then
      write(*,*) "     ********************************************"
      write(*,*) "     *                                          *"
      write(*,*) "     * KBPP_BAND_ORB-ATOMIC orbital Formatter   *"
      write(*,*) "     *                                          *"
      write(*,*) "     *        version last updated 01/01/02     *"
      write(*,*) "     *                                          *"
      write(*,*) "     *        developed by Eric J. Bylaska      *"
      write(*,*) "     *                                          *"
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
      write(*,111) "  cutoffs: ",(rc(i), i=0,lmax)

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
         write(*,*) "Generated formatted_filename: ",full_filename(1:l)
         end if
         call openfile(2,full_filename,l,'w',l)     
         call iwrite(2,version,1)
         call iwrite(2,ngrid,3)
         call dwrite(2,unita,9)
         call cwrite(2,atom,2)
         call iwrite(2,lmax,1)
         call iwrite(2,locp,1)
         call dwrite(2,rcut,1)
         call dwrite(2,lmbda,1)
         call iwrite(2,nbrillioun,1)
         call dwrite(2,kvectors,3*nbrillioun)
      end if

      do nb=1,nbrillioun

         if ((oprint).and.(taskid.eq.MASTER))
     >      write(*,*) "generating brillioun #",nb

         call integrate_kbpp_band_orb(version,
     >                      kvectors(1,nb),
     >                      nrho,drho,lmax,locp,
     >                                dbl_mb(wp_indx),
     >                                dbl_mb(rho_indx),
     >                                dbl_mb(f_indx),
     >                                dbl_mb(cs_indx),
     >                                dbl_mb(sn_indx),
     >                      nfft1,nfft2,nfft3,lmmax,
     >                                dbl_mb(G_indx),
     >                                dbl_mb(borbs_indx),
     >                      ierr)
         if (taskid.eq.MASTER) 
     >     call dwrite(2,dbl_mb(borbs_indx),nsize*lmmax)
      end do
      if (taskid.eq.MASTER) then
         call closefile(2)
      end if
     
*     **** free heap space ****
      value = BA_free_heap(rho_hndl)
      value = value.and.BA_free_heap(vp_hndl)
      value = value.and.BA_free_heap(wp_hndl)
      value = value.and.BA_free_heap(sc_r_hndl)

      value = value.and.BA_free_heap(f_hndl)
      value = value.and.BA_free_heap(cs_hndl)
      value = value.and.BA_free_heap(sn_hndl)

      value = value.and.BA_free_heap(borbs_hndl)
      value = value.and.BA_free_heap(G_hndl)

      if (.not.value)
     > call errquit("kbpp_band_orb: error freeing memory",0,MA_ERR)

      
      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)
      kbpp_band_orb = value
      return

 9999 call errquit('kbpp_band_orb:Error reading pspfilename',0,DISK_ERR)
      kbpp_band_orb = value
      END


         
      subroutine damp_kbpp_band_orb(nrho,lmax,rcut,lmbda,rho,w)
      implicit none
      integer nrho,lmax 
      real*8 rcut,lmbda,rho(nrho),w(nrho,0:lmax)

      integer i,l
    
      if (lmbda.gt.0.0d0) then
         do l=0,lmax
            do i=1,nrho 
              w(i,l) = w(i,l)*dexp(-(rho(i)/rcut)**lmbda)
            end do 
         end do
      end if
      return 
      end

