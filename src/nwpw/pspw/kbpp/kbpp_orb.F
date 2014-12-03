*
* $Id$
*
*     **************************************
*     *                                    *
*     *           kbpp_orb                 *
*     *                                    *
*     **************************************

      logical function kbpp_orb(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,locp,lmax,rcut,lmbda)
      implicit none
      logical          oprint_in
      integer          version
      character*50     psp_filename,formatted_filename
      integer          ngrid(3)
      double precision unita(3,3)
      integer locp,lmax
      real*8 rcut,lmbda

#include "errquit.fh"
#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"

c#include "tcgmsg.fh"
c#include "msgtypesf.h"

      character*255 full_filename

      logical value

      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmmax,k,lmax1,locp1
      double precision rc(0:9),rlocal1
      integer nrho
      double precision drho
      integer rho_indx,vp_indx,wp_indx,sc_r_indx
      integer rho_hndl,vp_hndl,wp_hndl,sc_r_hndl

      integer          isemicore
      logical          semicore
      double precision rcore

      integer f_indx,cs_indx,sn_indx
      integer f_hndl,cs_hndl,sn_hndl

*     ***** ngrid data *****
      integer aorbs_indx,G_indx
      integer aorbs_hndl,G_hndl

*     **** other variables ****
      double precision unitg(3,3),omega
      integer nsize,i,l,ierr
      integer nfft1,nfft2,nfft3
      logical mprint,dprint,oprint

*     **** external functions ****
      logical  control_print
      external control_print
      double precision simp
      external         simp


c      call Parallel_Init()
      call Parallel_taskid(taskid)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)
      dprint = (taskid.eq.MASTER).and.control_print(print_debug)
      oprint = (oprint_in.or.dprint)

      value = .false.
  
*     ***** read in pseudopotential data ****
      if (taskid.eq.MASTER) then
      call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')

      read(11,'(A2)') atom
      read(11,*) zv,amass,lmax0,lmax1,locp1,rlocal1
      read(11,*) (rc(i),i=0,lmax0)
      read(11,*) nrho,drho
      read(11,'(A)') comment
      end if

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
*      if (lmax.eq.-1)    lmax = lmax1
      !if (lmax.eq.-1)    lmax = lmax0  ! make lmax as large as possible
      !if (lmax.gt.lmax0) lmax = lmax0
      !if (lmax.lt.0)     lmax = lmax0

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
      value = value.and.BA_alloc_get(mt_dbl,nrho,
     >                        'f',f_hndl,f_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho,
     >                        'cs',cs_hndl,cs_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho,
     >                        'sn',sn_hndl,sn_indx)

*     **** allocate aorbs and  G ****
      lmmax = (lmax+1)**2 - (2*locp+1)
      nsize = (ngrid(1)/2+1)*ngrid(2)*ngrid(3)

      value = value.and.BA_alloc_get(mt_dbl,nsize*(lmmax),
     >                        'aorbs',aorbs_hndl,aorbs_indx)
      value = value.and.BA_alloc_get(mt_dbl,nsize*(3),
     >                        'G',G_hndl, G_indx)

      call damp_kbpp_orb(nrho,lmax,rcut,lmbda,
     >                   dbl_mb(rho_indx),dbl_mb(wp_indx))

*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp(nfft1,nfft2,nfft3,unita,unitg,dbl_mb(G_indx))

      call integrate_kbpp_orb(version,
     >                      nrho,drho,lmax,locp,
     >                                dbl_mb(wp_indx),
     >                                dbl_mb(rho_indx),
     >                                dbl_mb(f_indx),
     >                                dbl_mb(cs_indx),
     >                                dbl_mb(sn_indx),
     >                      nfft1,nfft2,nfft3,lmmax,
     >                                dbl_mb(G_indx),
     >                                dbl_mb(aorbs_indx),
     >                      ierr)



      if ((taskid.eq.MASTER).and.(oprint)) then
      write(luout,*) "     ********************************************"
      write(luout,*) "     *                                          *"
      write(luout,*) "     *     KBPP_ORB - Atomic Orbital Formatter  *"
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
      write(luout,*) "  Mulliken rcut   used           :",rcut
      write(luout,*) "  Mulliken lmbda  used           :",lmbda

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

*     **** normalize atomic orbitals ****
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      omega = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      omega = dabs(omega)




      if (taskid.eq.MASTER) then

      do k=1,nsize*lmmax
       dbl_mb(aorbs_indx+k-1) = dbl_mb(aorbs_indx+k-1)/dsqrt(omega)
      end do

      call util_file_name_noprefix(formatted_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      if (mprint) then
      write(luout,*)
      write(luout,*) "Generated formatted atomic orbitals, filename: ",
     >            full_filename(1:l)
      !write(*,*)
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
         call dwrite(2,dbl_mb(aorbs_indx),nsize*lmmax)
      call closefile(2)
      end if
     

*     **** free heap space ****
      value = value.and.BA_free_heap(rho_hndl)
      value = value.and.BA_free_heap(vp_hndl)
      value = value.and.BA_free_heap(wp_hndl)
      value = value.and.BA_free_heap(sc_r_hndl)
      value = value.and.BA_free_heap(f_hndl)
      value = value.and.BA_free_heap(cs_hndl)
      value = value.and.BA_free_heap(sn_hndl)

      value = value.and.BA_free_heap(aorbs_hndl)
      value = value.and.BA_free_heap(G_hndl)

      
      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)
      kbpp_orb = value
      return

 9999 call errquit('Error reading psp_filename',0, DISK_ERR)

      kbpp_orb = .FALSE.
      END


      subroutine damp_kbpp_orb(nrho,lmax,rcut,lmbda,rho,w)
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
  
