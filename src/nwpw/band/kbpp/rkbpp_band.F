*
* $Id$
*



*     **************************************
*     *                                    *
*     *           rkbpp_band               *
*     *   this is the full relativistic    *
*     *   KP pseudopotental                *
*     *                                    *
*     **************************************

*
*
      logical function rkbpp_band(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,locp_in,lmax_in,nbrillioun,
     >                  kvectors)
      implicit none
      logical          oprint_in
      integer          version
      character*50     psp_filename,formatted_filename
      integer          ngrid(3),lmax_in,locp_in
      double precision unita(3,3)
      integer nbrillioun
      real*8  kvectors(3,nbrillioun)

#include "errquit.fh"
#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "util.fh"


*     **** local variables ****
      character*255 full_filename
      logical value,mprint,hprint,oprint
      integer locp,nprj_out
      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmax1,locp1,nprj,ierr,nsph
      integer lmax_out
      double precision rc(0:12),rlocal1
      integer nrho
      double precision drho
      integer rho_indx,vp_indx,wp_indx,sc_r_indx,sc_k_indx
      integer rho_hndl,vp_hndl,wp_hndl,sc_r_hndl,sc_k_hndl

      integer          isemicore
      logical          semicore
      double precision rcore,core_charge

      integer f_indx,cs_indx,sn_indx
      integer n_prj_indx,l_prj_indx,m_alf_indx,m_bet_indx
      integer f_hndl,cs_hndl,sn_hndl
      integer n_prj_hndl,l_prj_hndl,m_alf_hndl,m_bet_hndl
      integer k_prj_hndl,k_prj_indx
      integer c_alf_hndl,c_alf_indx
      integer c_bet_hndl,c_bet_indx,pspx,kp

*     ***** ngrid data *****
      integer vl_indx,vnl_indx,vnrm_indx,G_indx
      integer vl_hndl,vnl_hndl,vnrm_hndl,G_hndl
      integer vpso_indx,vpso_hndl
*     **** other variables ****
      double precision unitg(3,3)
      integer nsize,i,l,nb,psp_type
      integer nfft1,nfft2,nfft3
      integer nmax

*     **** external functions ****
      logical  control_print
      external control_print
      double precision simp
      external         simp

      rkbpp_band=.true.
      call Parallel_taskid(taskid)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)
      oprint = (oprint_in.or.hprint)
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
*     ***** read in pseudopotential data ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
         l = index(full_filename,' ') - 1
         open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')
         read(11,*,ERR=9999,END=9999) pspx,atom
         read(11,*,ERR=9999,END=9999) zv,amass,lmax0,lmax1,
     >    locp1,rlocal1
         read(11,*,ERR=9999,END=9999) (rc(i),i=0,lmax0)
         read(11,*,ERR=9999,END=9999) nrho,drho
         read(11,'(A)',ERR=9999,END=9999) comment
      end if
       
      msglen = 1
      call BRDCST(9+MSGINT,pspx,mitob(msglen),MASTER)
      call BRDCST(9+MSGDBL,zv,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,amass,mdtob(msglen),MASTER)
      call BRDCST(9+MSGINT,lmax0,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,lmax1,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,locp1,mitob(msglen),MASTER)
      msglen = lmax0+1
      call BRDCST(9+MSGDBL,rc,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,nrho,mitob(msglen),MASTER)
      call BRDCST(9+MSGDBL,drho,mdtob(msglen),MASTER)


*     **** set the maximum angular momentum ****
      lmax_out=lmax1
      if (lmax_in.lt.0) lmax_out=lmax1
      if (lmax_out.gt.lmax0) lmax_out=lmax0       
        
*     **** set the local potential ****
      locp=locp1
      if (locp.gt.lmax1) locp = lmax1
      if (locp.lt.0)    locp = lmax1

      nprj=0
      do kp=0,lmax0-1
         nprj=nprj+4*kp+2
      end do
      nprj=nprj+2*lmax0
      if (nprj.lt.2) nprj=2
      nprj_out=0
      do kp=0,lmax_out-1
         nprj_out=nprj_out+4*kp+2
      end do
      nprj_out=nprj_out+2*lmax_out
      if (nprj_out.lt.2) nprj_out=2

*     **** allocate rho, vp, and wp ****
      value = BA_alloc_get(mt_dbl,nrho,
     >                     'rho',rho_hndl,rho_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(lmax0+1),
     >                        'vp',vp_hndl,vp_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(lmax0+1),
     >                        'vpso',vpso_hndl,vpso_indx)
      value = value.and.BA_alloc_get(mt_dbl,nrho*(2*lmax0+2),
     >                        'wp', wp_hndl,wp_indx)
      value = value.and.BA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)
      if (.not.value)
     > call errquit("rkbpp_band: out of heap memory",0,MA_ERR)



      if (taskid.eq.MASTER) then
       call read_rvpwp_band(11,nrho,lmax0,dbl_mb(rho_indx),
     >                         dbl_mb(vp_indx),
     >                         dbl_mb(wp_indx))
       call read_rkbpp_so(11,nrho,lmax0,dbl_mb(vpso_indx))
       call read_semicore_band(11,isemicore,rcore,nrho,
     > dbl_mb(sc_r_indx))
       close(11)
      end if

      msglen = nrho
      call BRDCST(9+MSGDBL,dbl_mb(rho_indx),mdtob(msglen),MASTER)
      msglen = nrho*(lmax0+1)
      call BRDCST(9+MSGDBL,dbl_mb(vp_indx),mdtob(msglen),MASTER)
      msglen = nrho*(2*lmax0+2)
      call BRDCST(9+MSGDBL,dbl_mb(wp_indx),mdtob(msglen),MASTER)
      msglen = nrho*(lmax0+1)
      call BRDCST(9+MSGDBL,dbl_mb(vpso_indx),mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,isemicore,mitob(msglen),MASTER)
      semicore = (isemicore.eq.1)
      if (semicore) then
       msglen = 2*nrho
       call BRDCST(9+MSGDBL,dbl_mb(sc_r_indx),mdtob(msglen),MASTER)
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

*     **** allocate vl,vnl,vnrm G ****
      nsize = ngrid(1)*ngrid(2)*ngrid(3)
      value = value.and.BA_alloc_get(mt_dbl,nsize,
     >                        'vl',vl_hndl,vl_indx)
      value = value.and.BA_alloc_get(mt_dcpl,(nsize*2*nprj_out),
     >                        'vnl',vnl_hndl, vnl_indx)
      value = value.and.BA_alloc_get(mt_dbl,nprj_out,
     >                        'vnrm', vnrm_hndl, vnrm_indx)
      value = value.and.BA_alloc_get(mt_dbl,(nsize*3),
     >                        'G',G_hndl, G_indx)
      value = value.and.BA_alloc_get(mt_dbl,(4*nsize),
     >                        'sc_k',sc_k_hndl,sc_k_indx)
      value = value.and.BA_alloc_get(mt_int,nprj_out,
     >                        'n_prj', n_prj_hndl, n_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj_out,
     >                        'l_prj', l_prj_hndl, l_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj_out,
     >                        'm_alf_prj', m_alf_hndl, 
     >                         m_alf_indx)
      value = value.and.BA_alloc_get(mt_int,nprj_out,
     >                        'm_bet_prj', m_bet_hndl, 
     >                         m_bet_indx)
      value = value.and.BA_alloc_get(mt_int,nprj_out,
     >			      'k_prj',k_prj_hndl, k_prj_indx)
      value = value.and.BA_alloc_get(mt_dbl,nprj_out,
     >                        'c_alf_prj', c_alf_hndl, 
     >                         c_alf_indx)
      value = value.and.BA_alloc_get(mt_dbl,nprj_out,
     >			      'c_bet_prj',c_bet_hndl, 
     >                         c_bet_indx)
      if (.not.value)
     > call errquit("rkbpp_band: out of heap memory",1,MA_ERR)

cccccccccc assign the relativistic states ccccccccccccccccccc
      call rkbpp_states(lmax_out,nprj_out,
     > int_mb(n_prj_indx),
     > int_mb(l_prj_indx),int_mb(k_prj_indx),
     > int_mb(m_alf_indx),int_mb(m_bet_indx),
     > dbl_mb(c_alf_indx),dbl_mb(c_bet_indx))

*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp_band(nfft1,nfft2,nfft3,unita,unitg,
     > dbl_mb(G_indx))
      call rkbpp_vlocal(lmax_out,locp,nsize,nrho,nprj_out,
     >       dbl_mb(rho_indx),
     >       dbl_mb(vp_indx),dbl_mb(vpso_indx),
     >       dbl_mb(vnrm_indx),dbl_mb(vl_indx),
     >       dbl_mb(wp_indx),
     >       int_mb(l_prj_indx),
     >       int_mb(k_prj_indx),
     >       dbl_mb(f_indx),dbl_mb(sn_indx),dbl_mb(cs_indx),
     >       dbl_mb(sc_k_indx),dbl_mb(sc_r_indx),
     >       dbl_mb(G_indx),semicore,zv,drho,ierr)

      if ((taskid.eq.MASTER).and.(oprint)) then
      write(*,*) "     ********************************************"
      write(*,*) "     *                                          *"
      write(*,*) "     *   RKBPP_BAND - Pseudopotential Formatter *"
      write(*,*) "     *                                          *"
      write(*,*) "     *        version last updated 03/17/07     *"
      write(*,*) "     *                                          *"
      write(*,*) "     *        developed by Patrick Nichols      *"
      write(*,*) "     *  based on KPBB_BAND by Eric Bylaska      *"
      write(*,*) "     ********************************************"
      call nwpw_message(1)
      write(*,*)
      write(*,*) "Pseudpotential Data"
      write(*,*) "-------------------"
      write(*,*) "  atom     :",atom
      write(*,*) "  charge   :",zv
      write(*,*) "  mass no. :",amass
      write(*,*) "  highest angular component      :",lmax0
      write(*,*) "  highest angular component used :",lmax_out
      write(*,*) "  local potential used           :",locp
      write(*,111) "  cutoffs: ",(rc(i), i=0,lmax_out)
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
      if (version.eq.3) write(*,112) "  boundary: periodic"
      write(*,113) "  ngrid  :",ngrid
      write(*,114) "  unita  :",unita(1,1),unita(2,1),unita(3,1)
      write(*,114) "          ",unita(1,2),unita(2,2),unita(3,2)
      write(*,114) "          ",unita(1,3),unita(2,3),unita(3,3)
      write(*,*)
      call flush(6)
  111 format(a,10f10.3)
  112 format(a)
  113 format(a,3I4)
  114 format(a,3F10.3)
  115 format(a,2E14.6)
      end if

        
      nb=1
      call rkbpp_vnonlocal(lmax_out,nrho,locp,
     >  nfft1,nfft2,nfft3,nprj_out,
     >  dbl_mb(G_indx),kvectors(1,1),
     >  dbl_mb(vp_indx),dbl_mb(vpso_indx),
     >  dbl_mb(wp_indx),dcpl_mb(vnl_indx),
     >  dbl_mb(vnrm_indx),
     >  dbl_mb(cs_indx),dbl_mb(sn_indx),dbl_mb(f_indx),
     >  dbl_mb(rho_indx),int_mb(l_prj_indx),
     >  int_mb(k_prj_indx),
     >  int_mb(m_alf_indx),int_mb(m_bet_indx),
     >  dbl_mb(c_alf_indx),dbl_mb(c_bet_indx),
     >  drho,ierr)

      if ((oprint).and.(taskid.eq.MASTER))
     >      write(*,*) "generating brillioun #",nb

      if (taskid.eq.MASTER) then 
        call util_file_name_noprefix(formatted_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        l = index(full_filename,' ') - 1
        write(*,*)
        write(*,*) "Generated formatted_filename: ",
     >            full_filename(1:l)
        call openfile(2,full_filename,l,'w',l)     
        call cwrite(2,comment,80)
        psp_type = 7
        call iwrite(2,psp_type,1)
        call iwrite(2,version,1)
        call iwrite(2,ngrid,3)
        call dwrite(2,unita,9)
        call cwrite(2,atom,2)
        call dwrite(2,amass,1)
        call dwrite(2,zv,1)
        call iwrite(2,lmax_out,1)
        call iwrite(2,locp,1)
c  nmax = 1 since we have 1 projector per (kappa) states
        nmax = 1
        call iwrite(2,nmax,1)
        call dwrite(2,rc,(lmax_out+1))
        call iwrite(2,nprj_out,1)
        if (nprj_out.gt.0) then
          call iwrite(2,int_mb(n_prj_indx),nprj_out)
          call iwrite(2,int_mb(l_prj_indx),nprj_out)
          call iwrite(2,int_mb(k_prj_indx),nprj_out)
          call dwrite(2,dbl_mb(vnrm_indx),nprj_out)
        end if  
        call dwrite(2,rcore,1)
        call iwrite(2,nbrillioun,1)
        call dwrite(2,kvectors,3*nbrillioun)
        call dwrite(2,dbl_mb(vl_indx),nsize)
        call dwrite(2,dcpl_mb(vnl_indx),nsize*nprj_out*4)
      end if

      do nb=2,nbrillioun

        if ((oprint).and.(taskid.eq.MASTER))
     >      write(*,*) "generating brillioun #",nb

        call rkbpp_vnonlocal(lmax_out,nrho,locp,
     >  nfft1,nfft2,nfft3,nprj_out,
     >  dbl_mb(G_indx),kvectors(1,nb),
     >  dbl_mb(vp_indx),dbl_mb(vpso_indx),
     >  dbl_mb(wp_indx),dcpl_mb(vnl_indx),
     >  dbl_mb(vnrm_indx),
     >  dbl_mb(cs_indx),dbl_mb(sn_indx),dbl_mb(f_indx),
     >  dbl_mb(rho_indx),
     >  int_mb(l_prj_indx),
     >  int_mb(k_prj_indx),
     >  int_mb(m_alf_indx),int_mb(m_bet_indx),
     >  dbl_mb(c_alf_indx),dbl_mb(c_bet_indx),
     >  drho,ierr)

        if (taskid.eq.MASTER) then
          call dwrite(2,dcpl_mb(vnl_indx),nsize*nprj_out*4)
        end if

      end do


      if (taskid.eq.MASTER) then
         if (semicore) then
           call dwrite(2,dbl_mb(sc_k_indx),4*nsize)
         end if
      call closefile(2)
      end if
     
*     **** free heap space ****
      value = BA_free_heap(rho_hndl)
      value = value.and.BA_free_heap(vp_hndl)
      value = value.and.BA_free_heap(vpso_hndl)
      value = value.and.BA_free_heap(wp_hndl)
      value = value.and.BA_free_heap(sc_r_hndl)
      value = value.and.BA_free_heap(sc_k_hndl)
      value = value.and.BA_free_heap(f_hndl)
      value = value.and.BA_free_heap(cs_hndl)
      value = value.and.BA_free_heap(sn_hndl)

      value = value.and.BA_free_heap(vl_hndl)
      value = value.and.BA_free_heap(vnl_hndl)
      value = value.and.BA_free_heap(vnrm_hndl)
      value = value.and.BA_free_heap(G_hndl)
      value = value.and.BA_free_heap(n_prj_hndl)
      value = value.and.BA_free_heap(l_prj_hndl)
      value = value.and.BA_free_heap(k_prj_hndl)
      value = value.and.BA_free_heap(m_alf_hndl)
      value = value.and.BA_free_heap(m_bet_hndl)
      value = value.and.BA_free_heap(c_alf_hndl)
      value = value.and.BA_free_heap(c_bet_hndl)
      if (.not.value)
     > call errquit("rkbpp_band: out of heap memory",2,MA_ERR)

      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)
      rkbpp_band = value

      return

 9999 call errquit('Error reading rel psp_filename',0, DISK_ERR)

      END



*     **********************************************
*     *                                            *
*     *            read_rvpwp_band                 *
*     *                                            *
*     **********************************************

c  Read in pseudowavefunction and pseudopotential

      subroutine read_rvpwp_band(unit,nrho,lmax,rho,vp,wp)
      implicit none
#include "errquit.fh"
      integer unit
      integer nrho,lmax,nprj
      double precision rho(nrho),rx
      double precision vp(nrho,0:lmax)
      double precision wp(nrho,2*lmax+2)
      integer i,j
      nprj= lmax*2+2
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rho(i),(vp(i,j),j=0,lmax)
      end do
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rx,(wp(i,j),j=1,nprj)
      end do
      return
 8000 FORMAT(10(1X,D15.7))
 9999 call errquit('Error reading Rel psp_filename',0, DISK_ERR)
      end


*     *******************************************
*     *                                         *
*     *            read_rkbpp_so                *
*     *                                         *
*     *******************************************

c     Read in spin orbit pseudopot.

      subroutine read_rkbpp_so(unit,nrho,lmax,vso)
      implicit none
      integer unit,nrho,lmax
      double precision vso(nrho,0:lmax)

#include "errquit.fh"

      integer i,j
      double precision rhox

      call dcopy((lmax+1)*nrho,0.0d0,0,vso,1)
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rhox,(vso(i,j),j=1,lmax)
      end do
      return
 9999 call errquit('Error reading Rel psp_so_filename',0, DISK_ERR)
      end


*     *******************************************
*     *                                         *
*     *           rkbpp_vlocal                  *
*     *                                         *
*     *******************************************

      subroutine rkbpp_vlocal(lmax,locl,ngsize,nrho,nprj,rho,
     >                        vp,vpso,vnlnrm,vl,wp,
     >                        l_prj,k_prj,
     >                        f,sinqr,cosqr,rho_sc_k,rho_sc_r,
     >                        gvec,semicore,zv,drho,ierr)
      implicit none
      logical semicore
      integer lmax
      integer nrho
      integer ngsize
      integer locl,ierr,nprj
      integer l_prj(*),k_prj(*)
      real*8 rho(nrho)
      real*8 f(nrho)
      real*8 sinqr(nrho)
      real*8 cosqr(nrho)
      real*8 vp(nrho,0:lmax)
      real*8 vpso(nrho,0:lmax)
      real*8 vl(ngsize)
      real*8 gvec(ngsize,3)
      real*8 wp(nrho,nprj)
      real*8 vnlnrm(nprj)
      real*8 rho_sc_k(ngsize,4)
      real*8 rho_sc_r(nrho,2)
      real*8 zv
      real*8 drho

c     **** local variables ****
      integer nproc,taskid,lp,ip,k,i,task_count,kappa,nw
      real*8 fourpi,vv,a,gx,gy,gz,q,sint,qr,dls
      logical small_cell

c     **** externals functions ****
      real*8   simp
      external simp
      logical  control_psp_semicore_small
      external control_psp_semicore_small

      call Parallel_np(nproc)
      call Parallel_taskid(taskid)

      small_cell = control_psp_semicore_small()

      fourpi=16.0*datan(1.0d0)
      do 10 lp=0,lmax
         if (lp.eq.locl) goto 10
         do k=1,nrho
             vp(k,lp)=vp(k,lp)-vp(k,locl)
         end do
10    continue

c     ***** assign norms to each projector *****
      do 20 ip=1,nprj
         lp=l_prj(ip)
         kappa=k_prj(ip)
         if (kappa.gt.0) then
            nw= 2*lp
         else
            nw= 2*lp+1
         end if
         dls=-0.5d0*(dble(kappa)+1.0d0)
         if (lp.eq.locl) then
            do k=1,nrho
               vv=dls*vpso(k,lp)  
               f(k)= vv*wp(k,nw)*wp(k,nw) 
            end do
            a=simp(nrho,f,drho)
            if (a.lt.1.d-12) then
              vnlnrm(ip)=0.0d0
            else
              vnlnrm(ip)=1.0d0/a
            end if
         else
            do k=1,nrho
              vv= vp(k,lp) + dls*vpso(k,lp)
              f(k)=vv*wp(k,nw)*wp(k,nw)
            end do
            a=simp(nrho,f,drho)
            vnlnrm(ip)=1.0d0/a
         end if
20    continue  

c     ***** calculate local potential in q space *****
      call dcopy(ngsize,0.0d0,0,vl,1)
      call dcopy(4*ngsize,0.0d0,0,rho_sc_k,1)  
      task_count= 0
      do 700 i=2,ngsize
         task_count=task_count+1
         if (mod(task_count,nproc).ne.taskid) then
            goto 700       
         end if 
         gx=gvec(i,1)
         gy=gvec(i,2)
         gz=gvec(i,3)
         q=dsqrt(gx**2+gy**2+gz**2)
         gx=gx/q
         gy=gy/q
         gz=gz/q
         do k=1,nrho
            qr=q*rho(k)  
            cosqr(k)=dcos(qr)
            sinqr(k)=dsin(qr)
            f(k)=rho(k)*vp(k,locl)*sinqr(k)
         end do
         sint=SIMP(nrho,f,drho)
         vl(i)=(sint*fourpi/Q)-(ZV*fourpi*cosqr(nrho)/(Q*Q))
         if (semicore) then
           if (small_cell) then
              do k=1,nrho
                 f(k) = rho(k)*(rho_sc_r(k,1))*sinqr(k)
              end do
           else
              do k=1,nrho
                 f(k) = rho(k)*dsqrt(rho_sc_r(k,1))*sinqr(k)
              end do
           end if
           rho_sc_k(i,1) = simp(nrho,f,drho)*fourpi/q
           do k=1,nrho
             f(k)=(sinqr(k)/(q*rho(k))-
     >            cosqr(k))*rho_sc_r(k,2)*rho(k)
           end do
           a = simp(nrho,f,drho)*fourpi/q
           rho_sc_k(i,2)=a*gx
           rho_sc_k(i,3)=a*gy
           rho_sc_k(i,4)=a*gz
         end if
  700 continue
      call Parallel_Vector_SumAll(4*ngsize,rho_sc_k)
      call Parallel_Vector_SumAll(ngsize,VL)
cccccccccc g=0
      do k=1,nrho
         f(k)=vp(k,locl)*rho(k)**2
      end do
      vl(1)=fourpi*(simp(nrho,f,drho)+0.5d0*zv*rho(nrho)**2)
      if (semicore) then
         if (small_cell) then
            do k=1,nrho
              f(k)=(rho_sc_r(k,1))*(rho(k)**2)
            end do
         else
            do k=1,nrho
              f(k)=dsqrt(rho_sc_r(k,1))*(rho(k)**2)
            end do
         end if
         rho_sc_k(1,1)=simp(nrho,f,drho)*fourpi
         rho_sc_k(1,2)=0.0d0
         rho_sc_k(1,3)=0.0d0
         rho_sc_k(1,4)=0.0d0
      end if
      ierr=0
      return 
      end


*     ****************************************
*     *                                      *
*     *          rkbpp_vnonlocal             *
*     *                                      *
*     ****************************************

      subroutine rkbpp_vnonlocal(lmax,nrho,locl,ngx,ngy,ngz,nprj,
     >                           gvec,kvec,
     >                           vp,vpso,wp,vnl,vnrm,
     >                           cosqr,sinqr,f,rho,
     >                           l_prj,k_prj,malf_prj,mbet_prj,calf_prj,
     >                           cbet_prj,drho,ierr)
      implicit none
      integer lmax,locl,nrho,ierr
      integer ngx,ngy,ngz,nprj
      integer malf_prj(nprj),mbet_prj(nprj),l_prj(nprj),k_prj(nprj)
      real*8 cosqr(nrho),sinqr(nrho),f(nrho),rho(nrho)
      real*8 vp(nrho,0:lmax),vpso(nrho,0:lmax)
      real*8 wp(nrho,2*lmax+2),vnrm(nprj)  
      complex*16 vnl(ngx,ngy,ngz,nprj,2)
      real*8 calf_prj(nprj)
      real*8 cbet_prj(nprj)
      real*8 drho,gvec(ngx,ngy,ngz,3),kvec(3)

#include "bafdecls.fh"

c     **** local variables ****
      complex*16 zzero
      parameter (zzero=(0.0d0,0.0d0))
      integer taskid,nproc,tcount,k
      integer igx,igy,igz,nw,lp,icnt
      integer malfa,mbeta,kappa,ippot,key
      real*8 pi4,q,gx,gy,gz,a,d,vv,qr,pi
      real*8 calfa,cbeta,dls

c     **** external functions ****
      real*8      simp,sjx
      complex*16  sphHarm
      external    sphHarm
      external    simp,sjx  

      call Parallel_taskid(taskid)
      call Parallel_np(nproc)
      pi4=16.0d0*datan(1.0d0)

      call dcopy(nprj*ngx*ngy*ngz*4,0.0d0,0,vnl,1)

      tcount=-1
      do 220 igz=1,ngz
      do 210 igy=1,ngy
      do 200 igx=1,ngx
         tcount=tcount+1
         icnt=tcount-(nproc*(tcount/nproc))
         if (icnt.ne.taskid) goto 200
         gx=gvec(igx,igy,igz,1)+kvec(1)
         gy=gvec(igx,igy,igz,2)+kvec(2)
         gz=gvec(igx,igy,igz,3)+kvec(3)
         q=dsqrt(gx*gx+gy*gy+gz*gz)
         if (dabs(q).lt.1.e-9) then 
           if (locl.ne.0) then
             do k=1,nrho
                f(k)=rho(k)*wp(k,1)*vp(k,0)
             end do
             d=simp(nrho,f,drho)*dsqrt(pi4)
             vnl(igx,igy,igz,1,2)=dcmplx(d,0.0d0)
             vnl(igx,igy,igz,2,1)=dcmplx(d,0.0d0)
           end if    
           goto 200
         end if
         gx=gx/q
         gy=gy/q
         gz=gz/q
         do k=1,nrho
            cosqr(k)=dcos(q*rho(k))
            sinqr(k)=dsin(q*rho(k))
         end do
         do 100 ippot=1,nprj
            lp=l_prj(ippot)
            kappa=k_prj(ippot)
            malfa=malf_prj(ippot)
            mbeta=mbet_prj(ippot)
            calfa=calf_prj(ippot)
            cbeta=cbet_prj(ippot)
            if (kappa.gt.0) then
               nw= 2*lp
            else
               nw= 2*lp+1
            end if
            dls= -0.5d0*(dble(kappa)+1.0d0)
            f(1)=0.0d0
            if (lp.eq.locl) then
               do k=2,nrho
                 vv=dls*vpso(k,lp)
                 vv=vv*wp(k,nw)
                 qr=q*rho(k)
                 a=sjx(qr,lp,cosqr(k),sinqr(k))
                 f(k)=a*vv*rho(k)
               end do        
            else
               if (lp.eq.0) then
                 vv=vp(1,0)
                 vv=vv*wp(1,1)
                 f(1)=vv*sinqr(1)/q
               end if
               do k=2,nrho
                 vv=vp(k,lp)+dls*vpso(k,lp)
      	         vv=vv*wp(k,nw)
                 qr=q*rho(k)
                 a=sjx(qr,lp,cosqr(k),sinqr(k))
                 f(k)=a*vv*rho(k)
               end do
            end if
            d=simp(nrho,f,drho)*pi4
            vnl(igx,igy,igz,ippot,1)=dcmplx(d*calfa,0.0d0)*
     >                 sphharm(lp,malfa,gx,gy,gz)       
            vnl(igx,igy,igz,ippot,2)=dcmplx(d*cbeta,0.0d0)*
     >                 sphharm(lp,mbeta,gx,gy,gz)
100   continue
200   continue 	
210   continue
220   continue
      call Parallel_Vector_SumAll(ngx*ngy*ngz*4*nprj,vnl)
      ierr=0
      return
      end


*     ****************************************
*     *                                      *
*     *              sjx                     *
*     *                                      *
*     ****************************************

c  spherical bessel function for a given L
c  note: no limit of L
c  but L<50 if x is small ....blows up
c  written by Patrick Nichols 6/11/2007
c  based upon subroutine version in patscat 2000

      real*8 function sjx(x,l,cosx,sinx)
      implicit none
      integer l
      real*8  x,cosx,sinx

c     **** local variables ****
      integer n,lstop,lmax
      real*8 tmp,ltest,tnp3,xinv,rho,xfact,xp,c0,c1,c2,c3,c4,c5,c6
      real*8 j0,j1,a1,a2

cccccccccccccccccccccccccccccccccccc   
      if (x.lt.0.1d0) then
        if (x.lt.1.0d-12) then
          sjx=0.0d0
          if (l.eq.0) sjx=1.0d0
          return
        end if
        xfact=x*x*0.5d0
        if (l.eq.0) then
            c1=-xfact*0.3333333333333333d0
            c2=-c1*xfact*0.1000000000000000d0
            c3=-c2*xfact*0.0238095238095238d0
            c4=-c3*xfact*0.0046296296296296d0
            c5=-c4*xfact*0.0007575757575758d0
            c6=-c5*xfact*0.0001068376068376d0
            sjx=(1.0d0+c1+c2+c3+c4+c5+c6)
            return
        end if  
        xp=x**l
        lstop=l+l+1
        c0=1.0d0
        do n=3,lstop,2
            c0=c0*dble(n)
        end do
        c1=-xfact/dble(2*l+3)
        c2=-c1*xfact*0.5d0/dble(2*l+5)
        c3=-c2*xfact*0.1666666666666667d0/dble(2*l+7)
        c4=-c3*xfact*0.0416666666666667d0/dble(2*l+9)
        c5=-c4*xfact*0.0083333333333333d0/dble(2*l+11)
        c6=-c5*xfact*0.0013888888888889d0/dble(2*l+13)
        sjx=xp*(1.0d0+c1+c2+c3+c4+c5+c6)/c0
        return  
      end if 
      xinv=1.0d0/x       
      j0=sinx*xinv
      j1=(sinx*xinv-cosx)*xinv
      if (l.eq.0) then
         sjx=j0
         return
      end if
      if (l.eq.1) then
         sjx=j1
         return
      end if    
      ltest=dble((l*(l+1)))
      if ((x*x).ge.ltest) then
         tnp3=3.0d0
         do n=2,l
           sjx=tnp3*xinv*j1-j0
           j0=j1
           j1=sjx
           tnp3=tnp3+2.0d0
         end do
         return
      end if
      lmax=l+40
      a2=0.0d0
      a1=1.0d0
      tnp3=2.0d0*dble(lmax)+3.0d0
      do n=lmax,l,-1
         sjx=tnp3*xinv*a1-a2
         a2=a1
         a1=sjx
         tnp3=tnp3-2.0d0
      end do
      tmp=sjx
      do n=l-1,0,-1
         sjx=tnp3*xinv*a1-a2
         a2=a1
         a1=sjx
         tnp3=tnp3-2.0d0
      end do
      rho=j0/sjx
      sjx=rho*tmp
      return
      end

*     ********************************************
*     *                                          *
*     *              rkbpp_states                *
*     *                                          *
*     ********************************************

c Assign the quantum numbers for the 2component relativstic 
c projectors. See Sakurai, Adv Quantum Mechanics for these
c numbers...pages 124,128
c               Patrick Jay Nichols

      subroutine rkbpp_states(lmax,nprj,np,lp,kp,ma,mb,ca,cb)
      implicit none
      integer nprj,lmax,np(*),lp(*),kp(*),ma(*),mb(*)
      real*8  ca(*),cb(*)

c     **** local variables ****
      integer kappa,j2,nz,iz,l,cnt
      real*8 ta,tb,td,dj,djz,dl

ccccccccccccccccccccccccc
c initialize counter
ccccccccccccccccccccccccc
      cnt = 0
ccccc do initial l=0 kappa=-1 state j=1/2
      td = 1.0d0
      kappa = -1
      j2 = 1
      dj = 0.5d0
      dl = 0.0d0
      do iz=1,2
        djz=dble(iz-1)-dj
        cnt = cnt+1   
        np(cnt) =  1
        lp(cnt) =  0
        kp(cnt) = -1
        ma(cnt) = int(nint(djz-0.5d0))
        mb(cnt) = int(nint(djz+0.5d0))
        ta = dl+0.5+djz
        tb = dl+0.5-djz
        if (ta.gt.1.d-12) then
           ca(cnt) = dsqrt(ta*td)
        else
           ca(cnt) = 0.0d0
        end if    
        if (tb.gt.1.d-12) then
           cb(cnt) = dsqrt(tb*td)
        else
           cb(cnt) = 0.0d0
        end if    
      end do
ccccccccccccccccccccccccc
      do l=1,lmax-1
         td = 1.0d0/dble(l+l+1)
         kappa = l
         j2 = l+l-1
         nz = j2+1
         dj = dble(j2)*0.5d0 
         dl = dble(l) 
         do iz=1,nz
            djz = dble(iz-1)-dj
            cnt = cnt+1   
            np(cnt) = 1
            lp(cnt) = l
            kp(cnt) = kappa
            ma(cnt) = int(nint(djz-0.5d0))
            mb(cnt) = int(nint(djz+0.5d0))
            ta = dl+0.5-djz
            tb = dl+0.5+djz
            if (ta.gt.1.d-12) then
                 ca(cnt) = -dsqrt(ta*td)
            else
                 ca(cnt) = 0.0d0
            end if    
            if (tb.gt.1.d-12) then
                 cb(cnt) = dsqrt(tb*td)
            else
                 cb(cnt) = 0.0d0
            end if    
         end do
         kappa= -(l+1)
         j2 = l+l+1
         nz = j2+1
         dj = dble(j2)*0.5d0  
         do iz=1,nz
            djz = dble(iz-1)-dj
            cnt = cnt+1   
            np(cnt) = 1
            lp(cnt) = l
            kp(cnt) = kappa
            ma(cnt) = int(nint(djz-0.5d0))
            mb(cnt) = int(nint(djz+0.5d0))
            ta = dl+0.5+djz
            tb = dl+0.5-djz
            if (ta.gt.1.d-12) then
                 ca(cnt) = dsqrt(ta*td)
            else
                 ca(cnt) = 0.0d0
            end if    
            if (tb.gt.1.d-12) then
                 cb(cnt) = dsqrt(tb*td)
            else
                 cb(cnt) = 0.0d0
            end if    
         end do
      end do
cccccccccc final lmax state
cccccccccc kappa= lmax, j=lmax-1/2
      td = 1.0d0/dble(lmax+lmax+1)
      kappa = lmax
      j2 = lmax+lmax-1
      nz = j2+1
      dj = dble(j2)*0.5d0
      dl = dble(lmax)  
      do iz=1,nz
        djz = dble(iz-1)-dj
        cnt = cnt+1   
        np(cnt) = 1
        lp(cnt) = lmax
        kp(cnt) = kappa
        ma(cnt) = int(nint(djz-0.5d0))
        mb(cnt) = int(nint(djz+0.5d0))
        ta = dl+0.5-djz
        tb = dl+0.5+djz
        if (ta.gt.1.d-12) then
          ca(cnt) = -dsqrt(ta*td)
        else
          ca(cnt) = 0.0d0
        end if    
        if (tb.gt.1.d-12) then
          cb(cnt) = dsqrt(tb*td)
        else
          cb(cnt) = 0.0d0
        end if    
      end do
      if (nprj.ne.cnt) then
        write(*,*)"NUMBER OF PROJECTORS GIVEN IS ",nprj
        write(*,*)"NUMBER OF PROJECTORS       IS ",cnt
        write(*,*)"POSSIBLE ERROR IN CODE!!!!!!!!"
        call flush(6)
      end if
      return
      end
