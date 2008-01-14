* $Id: rkbpp_band.F,v 1.19 2008-01-14 18:19:00 d3p708 Exp $
*
*     **************************************
*     *                                    *
*     *           rkbpp_band               *
*     *   this is the full relativistic    *
*     *   KP pseudopotental                *
*     *                                    *
*     **************************************
      logical function rkbpp_band(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,nbrillioun,
     >                  kvectors)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "util.fh"

      logical          oprint_in
      integer          version
      character*25     psp_filename,formatted_filename
      integer          ngrid(3)
      double precision unita(3,3)
      integer nbrillioun
      real*8  kvectors(3,nbrillioun)


*     **** local variables ****
      character*255 full_filename
      logical value,mprint,hprint,oprint
      integer lmax,locp,nprj_out
      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmax1,locp1,nprj,ierr
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
         write(*,*) pspx,atom
         write(*,*) zv,amass,lmax0,lmax1,locp1,rlocal1
         write(*,*) (rc(i),i=0,lmax0)
         write(*,*) nrho,drho
         write(*,*)'NRHO= ',NRHO,' DRHO= ',DRHO
         write(*,'(A)') comment
         write(*,*)"Reading psp full filename"
         call flush(6)
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
      lmax=lmax0
      lmax_out=lmax1
*     **** set the local potential ****
      if (locp.eq.-1)   locp = locp1
      if (locp.gt.lmax_out) locp = lmax_out
      if (locp.lt.0)    locp = lmax
      nprj=0  
      do kp=0,lmax-1
         nprj=nprj+4*kp+4
      end do
      nprj_out=0
      do kp=0,lmax_out-1
         nprj_out=nprj_out+4*kp+4
      end do
*     **** allocate rho, vp, and wp ****
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'rho',rho_hndl,rho_indx)
      value = MA_alloc_get(mt_dbl,nrho*(lmax0+1),
     >                        'vp',vp_hndl,vp_indx)
      value = MA_alloc_get(mt_dbl,nrho*(lmax0+1),
     >                        'vpso',vpso_hndl,vpso_indx)
      value = MA_alloc_get(mt_dbl,nrho*(2*lmax0+2),
     >                        'wp', wp_hndl,wp_indx)
      value = MA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)


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
      msglen = nrho*(lmax+1)
      call BRDCST(9+MSGDBL,dbl_mb(vp_indx),mdtob(msglen),MASTER)
      msglen = nrho*(2*lmax+2)
      call BRDCST(9+MSGDBL,dbl_mb(wp_indx),mdtob(msglen),MASTER)
      msglen = nrho*(lmax+1)
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
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'f',f_hndl,f_indx)
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'cs',cs_hndl,cs_indx)
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'sn',sn_hndl,sn_indx)

*     **** allocate vl,vnl,vnrm G ****
      nsize = ngrid(1)*ngrid(2)*ngrid(3)
      value = MA_alloc_get(mt_dbl,nsize,
     >                        'vl',vl_hndl,vl_indx)
      value = MA_alloc_get(mt_dcpl,(nsize*2*nprj),
     >                        'vnl',vnl_hndl, vnl_indx)
      value = MA_alloc_get(mt_dbl,nprj,
     >                        'vnrm', vnrm_hndl, vnrm_indx)
      value = MA_alloc_get(mt_dbl,(nsize*3),
     >                        'G',G_hndl, G_indx)
      value = MA_alloc_get(mt_dbl,(4*nsize),
     >                        'sc_k',sc_k_hndl,sc_k_indx)
      value = MA_alloc_get(mt_int,nprj,
     >                        'n_prj', n_prj_hndl, n_prj_indx)
      value = MA_alloc_get(mt_int,nprj,
     >                        'l_prj', l_prj_hndl, l_prj_indx)
      value = MA_alloc_get(mt_int,nprj,
     >                        'm_alf_prj', m_alf_hndl, 
     >                         m_alf_indx)
      value = MA_alloc_get(mt_int,nprj,
     >                        'm_bet_prj', m_bet_hndl, 
     >                         m_bet_indx)
      value = MA_alloc_get(mt_int,nprj,
     >			      'k_prj',k_prj_hndl, k_prj_indx)
      value = MA_alloc_get(mt_dbl,nprj,
     >                        'c_alf_prj', c_alf_hndl, 
     >                         c_alf_indx)
      value = MA_alloc_get(mt_dbl,nprj,
     >			      'c_bet_prj',c_bet_hndl, 
     >                         c_bet_indx)
cccccccccc assign the relativistic states ccccccccccccccccccc
      call rkbpp_states(lmax0,nprj,
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
      nsize=nfft1*nfft2*nfft3
      call rkbpp_vlocal(lmax0,locp,nsize,nrho,nprj,
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
      write(*,*) "  highest angular component used :",lmax
      write(*,*) "  local potential used           :",locp
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
      call flush(6)

      nb=1

      call rkbpp_vnonlocal(lmax0,nrho,locp,
     >  nfft1,nfft2,nfft3,nprj,
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
        call iwrite(2,lmax,1)
        call iwrite(2,locp,1)
c  nmax = 1 since we have 1 projector per (kappa) states
        nmax = 1
        call iwrite(2,nmax,1)
        call dwrite(2,rc,lmax+1)
        call iwrite(2,nprj_out,1)
        if (nprj.gt.0) then
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

        call rkbpp_vnonlocal(lmax,nrho,locp,
     >  nfft1,nfft2,nfft3,nprj,
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
      value = MA_free_heap(rho_hndl)
      value = MA_free_heap(vp_hndl)
      value = MA_free_heap(vpso_hndl)
      value = MA_free_heap(wp_hndl)
      value = MA_free_heap(sc_r_hndl)
      value = MA_free_heap(sc_k_hndl)
      value = MA_free_heap(f_hndl)
      value = MA_free_heap(cs_hndl)
      value = MA_free_heap(sn_hndl)

      value = MA_free_heap(vl_hndl)
      value = MA_free_heap(vnl_hndl)
      value = MA_free_heap(vnrm_hndl)
      value = MA_free_heap(G_hndl)
      value = MA_free_heap(n_prj_hndl)
      value = MA_free_heap(l_prj_hndl)
      value = MA_free_heap(k_prj_hndl)
      value = MA_free_heap(m_alf_hndl)
      value = MA_free_heap(m_bet_hndl)
      value = MA_free_heap(k_prj_hndl)
      value = MA_free_heap(c_alf_hndl)
      value = MA_free_heap(c_bet_hndl)

      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)
      rkbpp_band = value
      return

 9999 call errquit('Error reading rel psp_filename',0, DISK_ERR)

      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Read in pseudowavefunction and pseudopotential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
 9999 call errquit('Error reading Rel psp_filename',0, DISK_ERR)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Read in spin orbit pseudopot.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_rkbpp_so(unit,nrho,lmax,vso)
      implicit none
#include "errquit.fh"
      integer i,j,unit,nrho,lmax
      double precision vso(nrho,0:lmax)
      double precision rhox
ccccccccccccccccccc
      call dcopy((lmax+1)*nrho,0.0d0,0,vso,1)
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rhox,(vso(i,j),j=1,lmax)
      end do
      return
 9999 call errquit('Error reading psp_so_filename',0, DISK_ERR)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rkbpp_vlocal(lmax,locl,ngsize,nrho,nprj,rho,
     >       vp,vpso,vnlnrm,vl,wp,
     >       l_prj,k_prj,
     >       f,sinqr,cosqr,rho_sc_k,rho_sc_r,
     >       gvec,semicore,zv,drho,ierr)
      implicit none
      logical semicore
      integer lmax
      integer nrho
      integer ngsize
      integer locl,ierr,nprj
      integer l_prj(*),k_prj(*)
      double precision rho(nrho)
      double precision f(nrho)
      double precision sinqr(nrho)
      double precision cosqr(nrho)
      double precision vp(nrho,0:lmax)
      double precision vpso(nrho,0:lmax)
      double precision vl(ngsize)
      double precision gvec(ngsize,3)
      double precision wp(nrho,(2*lmax+2))
      double precision vnlnrm(nprj)
      double precision rho_sc_k(ngsize,4)
      double precision rho_sc_r(nrho,2)
      double precision zv
      double precision drho
cccccccc locals ccccccccccccccccccccccccccccccccccccccc       
      integer nproc,taskid,lp,ip,k,i,task_count,kappa,nw
      double precision fourpi,vv,a,gx,gy,gz,q,sint
cccccccc externals
      double precision simp
      external simp
ccccccccccccccccccccccccccccccccccccccccc
      call Parallel_np(nproc)
      call Parallel_taskid(taskid)
      fourpi=16.0*datan(1.0d0)
      do 10 lp=0,lmax
         if (lp.eq.locl) goto 10
         do k=1,nrho
             vp(k,lp)=vp(k,lp)-vp(k,locl)
         end do
10    continue
ccccccccccccc assign norms to each projector
      do 20 ip=1,nprj
         lp=l_prj(ip)
         if (lp.eq.locl) then
            vnlnrm(ip)=0.0d0
            goto 20
         end if
         kappa=k_prj(ip)
         if (kappa.gt.0) then
            nw= 2*lp
         else
            nw= 2*lp+1
         end if
         if ((nw.gt.(2*lmax+2)).or.(nw.lt.1)) then
           write(*,*)"ERROR! DANGER! RUN FOR YOUR LIVES!"
           write(*,*)"nw = ",nw," ip=",ip
           call flush(6)
         end if 
         do k=1,nrho
           vv= vp(k,lp) + dble(kappa)*vpso(k,lp)
           f(k)=vv*wp(k,nw)*wp(k,nw)
         end do
         a=simp(nrho,f,drho)
         vnlnrm(ip)=1.0d0/a
c         write(*,*)"ip =",ip," lp= ",lp," kp= ",kappa," nw=",nw
c         write(*,*)" nprj = ",nprj
c         write(*,*)"    vnrm(",ip,")=",vnlnrm(ip)
c         call flush(6)
20    continue  
cccccccccc calculate local potential in q space
      call dcopy(ngsize,0.0d0,0,vl,1)
      call dcopy(4*ngsize,0.0d0,0,rho_sc_k,1)  
      task_count= -1
      do 700 i=1,ngsize
         task_count=task_count+1
         if ((i.eq.1).or.(mod(task_count,nproc).ne.taskid)) then
            goto 700       
         end if 
         gx=gvec(i,1)
         gy=gvec(i,2)
         gz=gvec(i,3)
         q=dsqrt(gx*gx+gy*gy+gz*gz)
         gx=gx/q
         gy=gy/q
         gz=gz/q
         do k=1,nrho
            cosqr(k)=dcos(q*rho(k))
            sinqr(k)=dsin(q*rho(k))
            f(k)=rho(k)*vp(k,locl)*sinqr(k)
         end do
         sint=simp(nrho,f,drho) 
         vl(i)=sint*fourpi/q-zv*fourpi*cosqr(nrho)/q/q
         if (semicore) then
           do k=1,nrho
              f(k) = rho(k)*dsqrt(rho_sc_r(k,1))*sinqr(k)
           end do
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
         f(k)=vp(k,locl)*rho(k)*rho(k)
      end do
      vl(1)=fourpi*simp(nrho,f,drho)+0.5d0*fourpi*zv*rho(nrho)**2
      if (semicore) then
         do k=1,nrho
           f(k)=dsqrt(rho_sc_r(k,1))*rho(k)*rho(k)
         end do
         rho_sc_k(1,1)=fourpi*simp(nrho,f,drho)
         rho_sc_k(1,2)=0.0d0
         rho_sc_k(1,3)=0.0d0
         rho_sc_k(1,4)=0.0d0
      end if
      ierr=0
      return 
      end
cccccccccccccccccccccccccccc
      subroutine rkbpp_vnonlocal(lmax,nrho,locl,ngx,ngy,ngz,nprj,
     >  gvec,kvec,
     >  vp,vpso,wp,vnl,vnrm,
     >  cosqr,sinqr,f,rho,l_prj,k_prj,malf_prj,mbet_prj,calf_prj,
     >  cbet_prj,drho,ierr)
      implicit none
#include "mafdecls.fh"
      integer lmax,locl,nrho,ierr
      integer ngx,ngy,ngz,nprj
      integer malf_prj(*),mbet_prj(*),l_prj(*),k_prj(*)
      real*8 cosqr(nrho),sinqr(nrho),f(nrho),rho(nrho)
      real*8 vp(nrho,0:lmax),vpso(nrho,0:lmax)
      real*8 wp(nrho,2*lmax+2),vnrm(*)  
      complex*16 vnl(ngx,ngy,ngz,2,nprj)
      complex*16 zzero
      real*8 calf_prj(nprj)
      real*8 cbet_prj(nprj)
      real*8 drho,gvec(ngx,ngy,ngz,3),kvec(*)
cccccccc locals
      integer taskid,nproc,tcount,k
      integer igx,igy,igz,nw
      integer malfa,mbeta,kappa,l,ippot,key
      real*8 pi4,q,gx,gy,gz,a,d,vv,qr,pi
      real*8 calfa,cbeta,scof0
cccccc external
      real*8 simp,sjx
      complex*16 sphHarm
      external sphHarm
      external simp,sjx  
cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
      call Parallel_taskid(taskid)
      call Parallel_np(nproc)
      pi=4.0d0*datan(1.0d0)
      pi4=4.0d0*pi
      zzero=dcmplx(0.0d0,0.0d0)
      scof0=dsqrt(pi4)
      call dcopy(nprj*ngx*ngy*ngz*4,0.0d0,0,vnl,1)
      tcount=-1
      do 200 igz=1,ngz
      do 200 igy=1,ngy
      do 200 igx=1,ngx
         tcount=tcount+1
         key=tcount-nproc*(tcount/nproc)
         if (key.ne.taskid) goto 200
         gx=gvec(igx,igy,igz,1)+kvec(1)
         gy=gvec(igx,igy,igz,2)+kvec(2)
         gz=gvec(igx,igy,igz,3)+kvec(3)
         q=dsqrt(gx*gx+gy*gy+gz*gz)
         if (dabs(q).lt.1.e-9) then 
           if (locl.ne.0) then
             do k=1,nrho
               f(k)=rho(k)*wp(k,1)*vp(k,0)
             end do
             vnl(igx,igy,igz,1,1)=simp(nrho,f,drho)*scof0
             vnl(igx,igy,igz,2,1)=vnl(igx,igy,igz,1,1)*vnrm(1)
             vnl(igx,igy,igz,1,1)=vnl(igx,igy,igz,1,1)*vnrm(1)
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
            l=l_prj(ippot)
            if (l.eq.locl) then
              vnl(igx,igy,igz,1,ippot)=zzero
              vnl(igx,igy,igz,2,ippot)=zzero
              goto 100
            end if
            kappa=k_prj(ippot)
            if (kappa.gt.0) then
               nw= 2*l
            else
               nw= 2*l+1
            end if
            malfa=malf_prj(ippot)
            mbeta=mbet_prj(ippot)
            calfa=calf_prj(ippot)
            cbeta=cbet_prj(ippot)
            f(1)=0.0d0
            do k=1,nrho
               vv=vp(k,l)+dble(kappa)*vpso(k,l)
    	       vv=vv*wp(k,nw)
               qr=q*rho(k)
               a=sjx(qr,l,sinqr(k),cosqr(k))*qr
               f(k)=a*vv
            end do
            d=simp(nrho,f,drho)*pi4/q
            vnl(igx,igy,igz,1,ippot)=d*calfa*
     >                 sphHarm(l,malfa,gx,gy,gz)
            vnl(igx,igy,igz,2,ippot)=d*cbeta*
     >                 sphHarm(l,mbeta,gx,gy,gz)
            vnl(igx,igy,igz,1,ippot)=vnl(igx,igy,igz,1,ippot)*
     >              vnrm(ippot)
            vnl(igx,igy,igz,2,ippot)=vnl(igx,igy,igz,2,ippot)*
     >              vnrm(ippot)
100      continue
200   continue 	
      call C3db_Vector_SumAll(nprj*ngx*ngy*ngz*2,vnl)
      ierr=0
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
c  spherical bessel function for a given L
c  note: no limit of L
c  but L<50 if x is small ....blows up
c  written by Patrick Nichols 6/11/2007
c  based upon subroutine version in patscat 2000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
      double precision function SJX(qr,L,csx,snx)
      implicit none
      integer L
      double precision qr,csx,snx
      double precision fourpi
      parameter(fourpi=12.566370614359172d0)
      integer i,n
      integer nstop,nstrt,lmax
      real*8 rfact1,rfact2,rfact3,rfact4,factr
      real*8 xi,hxsq,series,sj0,sj1,fact,sjl,sjlm1,sjlp1
      REAL*8 sav,scale
      sjx=0.0d0
      if (QR.LT.1.D-12) THEN
         IF (L.EQ.0) SJX=1.0D0
         GOTO 777
      END IF   
      IF ( (QR.LT.0.1) .AND. (L.GE.1)) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC SMALL X EXPANSION 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         RFACT1=1.0D0
         NSTOP=L+L+1
         DO I=1,NSTOP,2
           RFACT1=RFACT1*DBLE(I)
         END DO
         HXSQ=QR*QR*0.5D0
         RFACT2=DBLE(L+L+3)
         RFACT3=DBLE(4*L+10)*RFACT2
         RFACT4=DBLE(2*L+7)*3.0D0*RFACT3
         SERIES=1.0D0-(HXSQ/RFACT2)+(HXSQ*HXSQ/RFACT3)
         SERIES=SERIES-(HXSQ**3)/RFACT4
         SJX=(QR**L)*SERIES/RFACT1  
         GOTO 777
      end if 
      xi=1.0d0/qr
      sj0=snx*xi
      sj1=(snx*xi-csx)*xi
      if (L.EQ.0) THEN
        SJX=SJ0
        GOTO 777
      END IF
      IF (L.EQ.1) THEN
        SJX=SJ1
        GOTO 777
      END IF
      HXSQ=QR*QR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC ASCENDING RECURSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (HXSQ.GE.DBLE(L*(L+1))) THEN
        DO I=2,L
          FACT=DBLE(I+I-1)
          SJX=FACT*SJ1*XI-SJ0
          SJ0=SJ1
          SJ1=SJX
        END DO
        GOTO 777
      END IF 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DOWNWARD RECURSE FOLLOWED BY SCALING
C SEE ABRAHAMOWITZ AND STEIGUN FOR THIS..
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      LMAX=L+25
      SJLP1=0.0D0
      SJL=1.D-34
      NSTRT=LMAX-1
      NSTOP=L+1
      DO N=(LMAX-1),(L+1),-1
        FACTR=DBLE(N+N+1)
        SJLM1=FACTR*SJL*XI-SJLP1
        SJLP1=SJL
        SJL=SJLM1
      END DO
      SAV=SJL
      DO N=L,1,-1
        FACTR=DBLE(N+N+1)
        SJLM1=FACTR*SJL*XI-SJLP1
        SJLP1=SJL
        SJL=SJLM1
      END DO
      SCALE=SJ0/SJL
      SJX=SCALE*SAV
777   RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Assign the quantum numbers for the 2component relativstic 
c projectors. See Sakurai, Adv Quantum Mechanics for these
c numbers...
c               PJN
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rkbpp_states(lmax,npr,n_prj,l_prj,k_prj,malf,mbet,
     > calf,cbet)
      implicit none
cccccccc arguments cccccccc
      integer lmax,npr
      integer n_prj(*),l_prj(*),k_prj(*),malf(*),mbet(*)
      real*8 calf(*),cbet(*)
cccccccc  locals cccccccccccc
      integer ip,lp,kp,is,ns,lpend
      real*8 j,jz,ta,tb,td
cccccccccccc s-states
      ip=0 
      lp=0 
      ns=2
      j=0.5d0
      kp= -1
      td=1.0d0
      do is=1,ns
          ip=ip+1
          jz=dble(is-1)-j
          n_prj(ip)=1
          l_prj(ip)=lp
          k_prj(ip)=kp
          ta=dabs(dble(lp)+0.5+jz)
          tb=dabs(dble(lp)+0.5-jz)
          if (ta.lt.1.e-12) then
             malf(ip)= lp+1
             calf(ip)= 0.0d0
          else
             malf(ip)= nint(jz-0.5)
             calf(ip)= dsqrt(ta/td)
          end if
          if (tb.lt.1.e-12) then
             mbet(ip)=lp+1
             cbet(ip)=0.0d0
          else
             mbet(ip)= nint(jz+0.5)
             cbet(ip)= dsqrt(tb/td)
          end if
      end do
cccccccccc do from l=1 to lmax
      lpend=lmax-1
      do lp=1,lpend
        ns=2*lp
        j=dble(lp)-0.5d0
        kp=lp
        td=dble(2*lp+1)
        do is=1,ns
          ip=ip+1
          jz=dble(is-1)-j
          n_prj(ip)=1
          l_prj(ip)=lp
          k_prj(ip)=kp
          ta=dabs(dble(lp)+0.5-jz)
          tb=dabs(dble(lp)+0.5+jz)
          if (ta.lt.1.e-12) then
             malf(ip)=lp+1
             calf(ip)=0.0d0
          else
             malf(ip)=nint(jz-0.5)
             calf(ip)= -dsqrt(ta/td)
          end if
          if (tb.lt.1.e-12) then
             mbet(ip)=lp+1
             cbet(ip)=0.0d0
          else
             mbet(ip)=nint(jz+0.5)
             cbet(ip)= dsqrt(tb/td)
          end if
        end do
        ns=2*lp+2
        j=dble(lp)+0.5d0
        kp= -lp-1
        do is=1,ns
          ip=ip+1
          jz=dble(is-1)-j
          n_prj(ip)=1
          l_prj(ip)=lp
          k_prj(ip)=kp
          ta=dabs(dble(lp)+0.5+jz)
          tb=dabs(dble(lp)+0.5-jz)
          if (ta.lt.1.e-12) then
             malf(ip)= lp+1
             calf(ip)=0.0d0
          else
             malf(ip)= nint(jz-0.5)
             calf(ip)= dsqrt(ta/td)
          end if
          if (tb.lt.1.e-12) then
             mbet(ip)=lp+1
             cbet(ip)=0.0d0
          else
             mbet(ip)= nint(jz+0.5)
             cbet(ip)= dsqrt(tb/td)
          end if
        end do
      end do
cccccccc last state = lmax
      lp=lmax
      ns=2*lp
      j=dble(lp)-0.5d0
      kp=lp
      td=dble(2*lp+1)
      do is=1,ns
          ip=ip+1
          jz=dble(is-1)-j
          n_prj(ip)=1
          l_prj(ip)=lp
          k_prj(ip)=kp
          ta=dabs(dble(lp)+0.5-jz)
          tb=dabs(dble(lp)+0.5+jz)
          if (ta.lt.1.e-12) then
             malf(ip)=lp+1
             calf(ip)=0.0d0
          else
             malf(ip)=nint(jz-0.5)
             calf(ip)= -dsqrt(ta/td)
          end if
          if (tb.lt.1.e-12) then
             mbet(ip)=lp+1
             cbet(ip)=0.0d0
          else
             mbet(ip)=nint(jz+0.5)
             cbet(ip)= dsqrt(tb/td)
          end if
      end do
      if (ip.ne.npr) then
        write(*,*)"NUMEBER OF PROJECTORS GIVEN IS ",npr
        write(*,*)"actual number is ",ip
        write(*,*)"NUMBER OF PROJECTORS DOES NOT EQUAL REAL NUMBER"
        call flush(6)
      END IF
      return
      end
ccccccccccccccccccccc
          
