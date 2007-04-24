
*** These modified files were found in JOE when it aborted on Thu Apr 19 18:21:33 2007
*** JOE was aborted because the terminal closed

*** File 'rkbpp_band.F'
*
* $Id: rkbpp_band.F,v 1.4 2007-04-24 18:38:47 d3p708 Exp $
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
     >                  ngrid,unita,locp,lmax,
     >                  nbrillioun,kvectors)
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
      integer locp,lmax
      integer nbrillioun
      real*8  kvectors(3,*)


*     **** local variables ****
      character*255 full_filename
      character*4 soend
      character*255 pspso_filename
      logical value,mprint,hprint,oprint

      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmmax,lmax1,locp1,nprj
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
      integer c_bet_hndl,c_bet_indx

*     ***** ngrid data *****
      integer vl_indx,vnl_indx,vnrm_indx,G_indx
      integer vl_hndl,vnl_hndl,vnrm_hndl,G_hndl
      integer vnlso_indx,vnlso_hndl
*     **** other variables ****
      double precision unitg(3,3)
      integer nsize,i,l,ierr,nb,psp_type
      integer nfft1,nfft2,nfft3
      integer nmax,nprj2

*     **** external functions ****
      logical  control_print
      external control_print
      double precision simp
      external         simp

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
        pspso_filename=full_filename//'.so'
        open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')
        read(11,'(A2)',ERR=9999,END=9999) atom
        read(11,*,ERR=9999,END=9999) zv,amass,lmax0,lmax1,
     >    locp1,rlocal1
        read(11,*,ERR=9999,END=9999) (rc(i),i=0,lmax0)
        read(11,*,ERR=9999,END=9999) nrho,drho
        read(11,'(A)',ERR=9999,END=9999) comment
      end if

      msglen = 1
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
      if (lmax.eq.-1)    lmax = lmax1
      if (lmax.gt.lmax0) lmax = lmax0
      if (lmax.lt.0)     lmax = lmax0

*     **** set the local potential ****
      if (locp.eq.-1)   locp = locp1
      if (locp.gt.lmax) locp = lmax
      if (locp.lt.0)    locp = lmax

      

*     **** allocate rho, vp, and wp ****
      nprj=2*lmax+1
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'rho',rho_hndl,rho_indx)
      value = MA_alloc_get(mt_dbl,nrho*(lmax+1),
     >                        'vp',vp_hndl, vp_indx)
      value = MA_alloc_get(mt_dbl,nrho*nprj,
     >                        'vpso',vp_hndl, vp_indx)
      value = MA_alloc_get(mt_dbl,nrho*nprj,
     >                        'wp', wp_hndl, wp_indx)
      value = MA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)

      if (taskid.eq.MASTER) then
      call read_rvpwp_band(11,nrho,lmax,dbl_mb(rho_indx),
     >                         dbl_mb(vp_indx),
     >                         dbl_mb(wp_indx))
      call read_semicore_band(11,isemicore,rcore,nrho,
     > dbl_mb(sc_r_indx))
      close(11)
  
      end if
      if (taskid.eq.MASTER) then
         l = index(pspso_filename,' ') - 1
         open(unit=11,file=pspso_filename(1:l),
     >             status='old',form='formatted')
         call read_rkbpp_so(11,nrho,nprj,dbl_mb(vnlso_indx))
      end if

      msglen = nrho
      call BRDCST(9+MSGDBL,dbl_mb(rho_indx),mdtob(msglen),MASTER)
      msglen = nrho*(lmax+1)
      call BRDCST(9+MSGDBL,dbl_mb(vp_indx),mdtob(msglen),MASTER)
      msglen = nrho*nprj
      call BRDCST(9+MSGDBL,dbl_mb(wp_indx),mdtob(msglen),MASTER)
      msglen = nrho*nprj
      call BRDCST(9+MSGDBL,dbl_mb(vnlso_indx),mdtob(msglen),MASTER)
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
      lmmax = 2*(lmax+1)**2 - 1
      nprj  = 2*lmax+1
      nprj2 = nprj+nprj
      nsize = (ngrid(1))*ngrid(2)*ngrid(3)
      value = MA_alloc_get(mt_dbl,nsize*2,
     >                        'vl',vl_hndl,vl_indx)
      value = MA_alloc_get(mt_dbl,nsize*nprj2,
     >                        'vnl',vnl_hndl, vnl_indx)
      value = MA_alloc_get(mt_dbl,nsize*nprj2,
     >                        'vnlso',vnlso_hndl, vnlso_indx)
      value = MA_alloc_get(mt_dbl,nprj,
     >                        'vnrm', vnrm_hndl, vnrm_indx)
      value = MA_alloc_get(mt_dbl,nsize*(3),
     >                        'G',G_hndl, G_indx)
      value = MA_alloc_get(mt_dbl,4*nsize,
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
     >                        c_bet_indx)


*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp_band(nfft1,nfft2,nfft3,unita,unitg,
     > dbl_mb(G_indx))
      nsize=nfft1*nfft2*nfft3
      call init_rkbpp_vlocal(locp,nprj,nsize,nrho,dbl_mb(rho_indx),
     >  dbl_mb(G_indx),dbl_mb(vp_indx),dbl_mb(vl_indx),lmax,
     >  int_mb(n_prj_indx),int_mb(l_prj_indx),int_mb(k_prj_indx),
     >  int_mb(m_alf_indx),
     >  int_mb(m_bet_indx),
     >  dbl_mb(c_alf_indx),
     >  dbl_mb(c_bet_indx),
     >  dbl_mb(sc_r_indx),
     >  dbl_mb(sc_k_indx),semicore,dbl_mb(f_indx),dbl_mb(sn_indx),
     >  dbl_mb(cs_indx),drho,zv)


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
  111 format(a,10f10.3)
  112 format(a)
  113 format(a,3I4)
  114 format(a,3F10.3)
  115 format(a,2E14.6)
      end if


      call init_rkbpp_nonlocal(locp,nsize,lmax,nprj,
     >  int_mb(l_prj_indx),int_mb(k_prj_indx),
     >  int_mb(m_alf_indx),int_mb(m_bet_indx),
     >  dbl_mb(c_alf_indx),dbl_mb(c_bet_indx),
     >  nrho,dbl_mb(rho_indx),dbl_mb(G_indx),
     >  kvectors(1,nb),dcpl_mb(vnl_indx),
     >  dcpl_mb(vnlso_indx),
     >  dbl_mb(vp_indx),dbl_mb(vnrm_indx),
     >  dbl_mb(f_indx),dbl_mb(sn_indx),dbl_mb(cs_indx),drho)

      nb=1

      if ((oprint).and.(taskid.eq.MASTER))
     >      write(*,*) "generating brillioun #",nb

      if (taskid.eq.MASTER) then 
        call util_file_name_noprefix(formatted_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        l = index(full_filename,' ') - 1
        if (mprint) then
           write(*,*)
           write(*,*) "Generated formatted_filename: ",
     >            full_filename(1:l)
        end if
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
        call iwrite(2,nprj,1)

        nmax = 1
        call iwrite(2,nmax,1)
        call dwrite(2,rc,lmax+1)
        call iwrite(2,nprj,1)
        if (nprj.gt.0) then
          call iwrite(2,int_mb(n_prj_indx),nprj)
          call iwrite(2,int_mb(l_prj_indx),nprj)
          call iwrite(2,int_mb(k_prj_indx),nprj)
          call iwrite(2,int_mb(m_alf_indx),nprj)
          call iwrite(2,int_mb(m_bet_indx),nprj)
          call dwrite(2,int_mb(c_alf_indx),nprj)
          call dwrite(2,int_mb(c_bet_indx),nprj)
          call dwrite(2,dbl_mb(vnrm_indx),nprj)
        end if         
        call dwrite(2,rcore,1)
        call iwrite(2,nbrillioun,1)
        call dwrite(2,kvectors,3*nbrillioun) 
        call dwrite(2,dbl_mb(vl_indx),nsize)
        call dwrite(2,dbl_mb(vnl_indx),nsize*nprj*4)
        call dwrite(2,dbl_mb(vnlso_indx),nsize*nprj*4)
      end if


      do nb=2,nbrillioun

        if ((oprint).and.(taskid.eq.MASTER))
     >      write(*,*) "generating brillioun #",nb

        call init_rkbpp_nonlocal(locp,nsize,lmax,nprj,
     >   int_mb(l_prj_indx),int_mb(k_prj_indx),
     >   int_mb(m_alf_indx),int_mb(m_bet_indx),
     >   dbl_mb(c_alf_indx),dbl_mb(c_bet_indx),
     >   nrho,dbl_mb(rho_indx),dbl_mb(G_indx),
     >   kvectors(1,nb),dcpl_mb(vnl_indx),
     >   dcpl_mb(vnlso_indx),
     >   dbl_mb(vp_indx),dbl_mb(vnrm_indx),
     >   dbl_mb(f_indx),dbl_mb(sn_indx),dbl_mb(cs_indx),drho)

        if (taskid.eq.MASTER) then 
           call dwrite(2,dbl_mb(vnl_indx),nsize*nprj*4)
           call dwrite(2,dbl_mb(vnlso_indx),nsize*nprj*4)
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
      value = MA_free_heap(wp_hndl)
      value = MA_free_heap(sc_r_hndl)
      value = MA_free_heap(sc_k_hndl)
      value = MA_free_heap(f_hndl)
      value = MA_free_heap(cs_hndl)
      value = MA_free_heap(sn_hndl)

      value = MA_free_heap(vl_hndl)
      value = MA_free_heap(vnlso_hndl)
      value = MA_free_heap(vnl_hndl)
      value = MA_free_heap(vnrm_hndl)
      value = MA_free_heap(G_hndl)
      value = MA_free_heap(n_prj_hndl)
      value = MA_free_heap(l_prj_hndl)
      value = MA_free_heap(m_alf_hndl)
      value = MA_free_heap(m_bet_hndl)
      value = MA_free_heap(k_prj_hndl)
      value = MA_free_heap(c_alf_hndl)
      value = MA_free_heap(c_bet_hndl)

      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)
      rkbpp_band = value
      return

 9999 call errquit('Error reading psp_filename',0, DISK_ERR)

      END

      subroutine read_rvpwp_band(unit,nrho,lmax,rho,vp,wp)
      implicit none
#include "errquit.fh"
      integer unit
      integer nrho,lmax,npot,nprj
      double precision rho(nrho),rx
      double precision vp(nrho,0:lmax)
      double precision wp(nrho,2*lmax+1)
      integer i,j
      npot= lmax+1
      nprj= lmax*2+1
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rho(i),(vp(i,j),j=0,npot)
      end do
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rx,(wp(i,j),j=1,nprj)
      end do
      
      return

 9999 call errquit('Error reading psp_filename',0, DISK_ERR)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Read in spin orbit pseudopot.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_rkbpp_so_band(unit,nrho,nprj,vso)
      integer i,j,unit
      double precision vso(nrho,nprj)
      double precision rhox
      call dcopy(nprj*nrho,0.0d0,0,vso,1)
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rhox,(vso(i,j),j=2,nprj)
      end do
 9999 call errquit('Error reading psp_so_filename',0, DISK_ERR)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_rel_band_vnonlocal(locp,ngsize,lmax,nprj,
     >  lprj,kprj,malf,mbet,calf,cbet,nrho,rho,gvecs,kvec,
     >  vnl,wp,vnlso,vp,vnrm,f,sn,cs,drho)
      implicit none
      integer nprj,nrho,ngx,ngy,ngz,ngsize,lmax,locp
      integer lprj(nprj)
      integer kprj(nprj)  
      integer malf(nprj)
      integer mbet(nprj)  
      complex*16 vnl(ngsize,nprj,2)
      complex*16 vnlso(ngsize,nprj,2)
      real*8 wp(nrho,nprj)
      real*8 vnrm(nprj)
      real*8 vp(nrho,0:lmax)
      real*8 vso(nrho,nprj)
      real*8 gvecs(3,ngsize)
      real*8 kvec(3)
      real*8 calf(nprj)
      real*8 cbet(nprj)
      real*8 rho(nrho)
      real*8 F(nrho)
      real*8 cs(nrho)
      real*8 sn(nrho)
      real*8 drho
ccccc  local variables ccccccccccccccccccccccccccccccccccccccc
      integer ip,ix,iy,iz,i
      integer malfa,mbeta,l,kappa,l1,locl
      real*8 Fso(nrho),d,calfa,cbeta
      real*8 gx,gy,gz,gr,q,A,p0,qr
      complex*16 ylma,ylmb,spherical_bandx
      real*8 sjx,simp 
      external spherical_bandx
      external sjx,simp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      locl=lprj(locp)
      p0=1.0d0/dsqrt(16.0*datan(1.0d0))
      do 600 ip=1,nprj
        l=lprj(ip)
        kappa=kprj(ip)
        malfa=malf(ip)
        mbeta=mbet(ip)
        calfa=calf(ip)
        cbeta=cbet(ip)
        do 400 ix=1,ngsize
          gx=gvecs(1,ix)+kvec(1)
          gy=gvecs(2,ix)+kvec(2)
          gz=gvecs(3,ix)+kvec(3)
          gr=dsqrt(gx*gx+gy*gy+gz*gz)
          if (gr.ge.1.d-10) then
            gx=gx/gr
            gy=gy/gr
            gz=gz/gr
            do i=1,nrho
              sn(i)= dsin( gr*rho(i) )
              cs(i)= dsqrt( 1.0d0-sn(i)*sn(i) )
            end do
            if (locl.ne.l) then
              F(1)=0.0d0
              Fso(1)=0.0d0
              do i=1,nrho
                qr = gr*rho(i)
                A=sjx(qr,sn(i),cs(i),l)
                F(i)=A*vp(i,l)*wp(i,ip)
                Fso(i)=A*vso(i,ip)*wp(i,ip)
              end do
              d=simp(nrho,f,drho)/gr
              ylma=spherical_bandx(l,malfa,gx,gy,gz)
              ylmb=spherical_bandx(l,mbeta,gx,gy,gz)
              vnl(ix,ip,1)=ylma*(d*calfa)
              vnl(ix,ip,2)=ylmb*(d*cbeta)
              d=simp(nrho,fso,drho)/gr
              vnlso(ix,ip,1)=ylma*(d*calfa*kappa)
              vnlso(ix,ip,2)=ylmb*(d*cbeta*kappa)
            else
              F(1)=0.0d0
              Fso(1)=0.0d0
              do i=2,nrho
                qr=gr*rho(i)
                A=sjx(qr,sn(i),cs(i),l)
                Fso(i)=A*vso(i,ip)*wp(i,ip)
              end do
              ylma=spherical_bandx(l,malfa,gx,gy,gz)
              ylmb=spherical_bandx(l,mbeta,gx,gy,gz)
              d=simp(nrho,fso,drho)/gr
              vnlso(ix,locp,1)=ylma*(d*calfa*kappa)
              vnlso(ix,locp,2)=ylmb*(d*cbeta*kappa)
            end if                                          
          else
            if (l.eq.0 .and. locl.ne.0) then
                F(1)=0.0d0
                do i=2,nrho
                  F(i)=rho(i)*wp(i,1)*vp(i,1)
                end do
                d=simp(nrho,f,drho) 
                if (malfa.le.l) then
                 vnl(ix,ip,1)=d*calfa*p0
                end if
                if (mbeta.le.l) then
                 vnl(ix,ip,2)=d*cbeta*p0
                end if
            end if
          end if
400     continue
600   continue
      call C3db_Vector_Sumall(2*nprj*ngsize,vnl)
      call C3db_Vector_Sumall(2*nprj*ngsize,vnlso)
      do ip=1,nprj
         f(1)=0.d0 
         l1=lprj(ip)
         kappa=kprj(ip)
         do i=1,nrho 
           A=wp(i,ip)*wp(i,ip)
           f(i)=A*(vp(i,l1)+kappa*vso(i,ip))
         end do
         d=simp(nrho,f,drho)
         vnrm(ip)=1.0d0/d
      end do        
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_rkbpp_vlocal(locp,nprj,ngsize,nrho,rho,gvec,
     > vp,vpso,vl,lmax,n_prj,lprj,kprj,malf,mbet,calf,cbet,
     > rho_sc_r,rho_sc_k,semicore,f,sn,cs,drho,zv)
      implicit none
      logical semicore
      integer nrho,nprj,ngsize,locp,lmax
      integer n_prj(*)
      integer lprj(*)
      integer malf(*)
      integer mbet(*)
      integer kprj(*)
      real*8 rho(*),drho,zv   
      real*8 vp(nrho,0:lmax),vl(nrho),vpso(nrho,nprj)
      real*8 gvec(3,ngsize)
      real*8 calf(*)
      real*8 cbet(*)
      real*8 rho_sc_k(ngsize,4)
      real*8 rho_sc_r(nrho,2) 
      real*8 cs(nrho),sn(nrho),f(nrho)
cccccc local
      integer knt,ip,i,k,locl,kappa,lp
      integer task_count,np,taskid,MASTER
      real*8 gx,gy,gz,gr,q,forpi,twopi
      real*8 simp,d,p
      external simp
cccccccccccccccccccccccccccccccccccccccccccc
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      forpi=16.0d0*datan(1.0d0)
      twopi=forpi*0.5d0
      call make_rel_band_states(lmax,n_prj,lprj,kprj,malf,mbet,
     > calf,cbet)
      knt=1
      locl=lprj(locp)
      do lp=0,lmax
        if (lp.ne.locl) then        
          do k=1,nrho
             vp(k,lp)=vp(k,lp)-vp(k,locl)
          end do
        end if
      end do
      call dcopy(ngsize,0.0,0,vl,1)
      call dcopy(4*ngsize,0.0d0,0,rho_sc_k,1)
      task_count=-1
      do 700 i=1,ngsize
        task_count=task_count+1
        if ( (np*(task_count/np)).ne.taskid) cycle
	gx=gvec(1,i)
	gy=gvec(2,i)
	gz=gvec(3,i)
	gr=dsqrt(gx*gx+gy*gy+gz*gz)
	if (i.eq.1) cycle
	gx=gx/gr
        gy=gy/gr
        gz=gz/gr
        do k=1,nrho
          cs(k)=dcos(gr*rho(k))
          sn(k)=dsin(gr*rho(k))
        end do
600     continue
        do k=1,nrho
          f(k)=rho(k)*vp(k,locl)*sn(k)
        end do
        vl(i)=simp(nrho,f,drho)*forpi/gr-Zv*forpi*cs(nrho)/(gr*gr)
        if (semicore) then
          do k=1,nrho  
	    f(k)=rho(k)*dsqrt(rho_sc_r(k,1))*sn(k)
          end do
          rho_sc_k(i,1)=simp(nrho,f,drho)*forpi/gr
          do k=1,nrho
            q=gr*rho(k)
            f(k)=(sn(k)/q-cs(k))*rho_sc_r(k,2)*rho(k)
          end do
          p=simp(nrho,f,drho)*forpi/gr
	  rho_sc_k(i,2)=p*gx
          rho_sc_k(i,3)=p*gy
          rho_sc_k(i,4)=p*gz
        end if  		           	
700   end do
      call D3db_Vector_SumAll(4*ngsize,rho_sc_k)
      call C3db_Vector_SumAll(ngsize,vl)
ccccccccccccccc  G=0
      do k=1,nrho
        f(k)=vp(k,locl)*rho(k)*rho(k)
      end do
      vl(1)=forpi*simp(nrho,f,drho)+twopi*zv*rho(nrho)*rho(nrho)
cccccccccccc semicore density
      if (semicore) then
        do k=1,nrho
          f(k)=dsqrt(rho_sc_r(k,i))*rho(k)*rho(k)
        end do
        rho_sc_k(1,1)=forpi*simp(nrho,f,drho)
        rho_sc_k(1,2)=0.0d0
        rho_sc_k(1,3)=0.0d0
        rho_sc_k(1,4)=0.0d0
      end if
c      ierr=0
      return
      end                       	 
cccccccccccccccccccccccccccccccccccccccccccc         
      complex*16 function spherical_bandx(l,m,gx,gy,gz)
      implicit none
      integer l,m
      real*8 gx,gy,gz
cccccccc local variables
      integer am
      real*8 twopi,pi,forpi,scof(7,0:3)
      real*8 gxy,sin_theta,phi
      complex*16 tmp,expm
ccccccccc local variables   
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi
      scof(1,0)=dsqrt(1.0d0/forpi)
      scof(1,1)=dsqrt(0.75d0/pi)
      scof(2,1)=dsqrt(0.75/2.0d0/forpi)
      scof(3,1)= -scof(1,1)
      scof(1,2)=dsqrt(15.0d0/32.0d0/pi)
      scof(2,2)=dsqrt(15.0d0/8.0d0/pi)
      scof(3,2)=dsqrt(5.0d0/16.0d0/pi)
      scof(4,2)=-scof(2,2)
      scof(5,2)=scof(2,1)
      scof(1,3)=dsqrt(35.0d0/64.0d0/pi)
      scof(2,3)=dsqrt(105.d0/32.0d0/pi)
      scof(3,3)=dsqrt(21.0d0/64.0d0/pi)
      scof(4,3)=dsqrt(7.0d0/16.0d0/pi)
      scof(5,3)= -scof(3,3)
      scof(6,3)=scof(3,2)
      scof(7,3)= -scof(3,1)
      if (m.gt.l .or. m.lt.(-l)) then
         spherical_bandx=dcmplx(0.0d0,0.0d0)
         return
      end if
      if (l.gt.3) then
         write(*,*)"AKKK  l.gt.3 in spherical_band_x"
         spherical_bandx=dcmplx(0.0d0,0.0d0)
         stop
      end if         
      gxy=dsqrt(gx*gx+gy*gy)
      sin_theta=dsin(datan2(gxy,gz))
      phi=datan2(gy,gx)
      expm=dcmplx(dcos(dble(m)*phi),dsin(dble(m)*phi))
      tmp=scof(m+l+1,l)*expm
      am=iabs(m)
      if (l.lt.3) goto 100
      if (am.eq. 3) then
        spherical_bandx=(tmp*(sin_theta*sin_theta*sin_theta)) 
        return
      end if
      if (am.eq.2) then
        spherical_bandx=(tmp*(gz*sin_theta*sin_theta))
        return
      end if
      if (am.eq.1) then
        spherical_bandx=tmp*(gz*gz*5.0d0-1.0)*sin_theta
        return
      end if
      spherical_bandx=tmp*(5.0d0*gz*gz*gz-3.0d0*gz)
      return
100   if (l.lt.2) goto 200
      if (am.eq.2) then
        spherical_bandx=tmp*(sin_theta*sin_theta)   
        return
      end if
      if (am.eq.1) then
        spherical_bandx=tmp*sin_theta*gz
        return
      end if
      spherical_bandx=tmp*(3.0d0*gz*gz-1.0d0)  
      return
200   if (l.lt.1) goto 300
      if (am.eq.1) then
        spherical_bandx=tmp*(sin_theta)
        return
      end if
      spherical_bandx=tmp*gz
      return
300   spherical_bandx=tmp
      return
      end
cccccccccccccccccccccccc
      real*8 function sjx(qr,sn,cs,l)
      implicit none
      real*8 qr,sn,cs,a
      integer l
      if (l.eq.3) then
         a=sn/qr
         sjx=15.d0*(a-cs)/qr/qr-6.0d0*A+cs
         return
      end if
      if (l.eq.2) then
         a=sn/qr
         sjx=3.0d0*(a-cs)/qr-sn
         return
      end if
      if (l.eq.1) then
         a=sn/qr
         sjx=a-cs 
         return
      end if
      sjx=sn
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine make_rel_band_states(lmax,nprj,n_prj,lprj,kprj,
     > malf,mbet,calf,cbet)
      implicit none
      integer lmax
      integer n_prj(nprj),lprj(nprj),kprj(nprj)
      integer malf(nprj),mbet(nprj)
      real*8 calf(nprj),cbet(nprj)
ccccc locals ccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer la,n3,jmax,knt,nprj,kappa
      real*8 j,jz,td,dla,ta,tb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nprj=lmax*2+1
      knt=0
      do la=0,lmax
          kappa= -(la+1.d0)  
          j=dble(la)+0.5d0
          jmax=2*la+2
          do n3=0,jmax
              jz=dble(n3)-j
              n_prj(knt)=1
              lprj(knt)=la
              kprj(knt)=int(kappa)
              td=dble(la+la+1)
              dla=dble(la)
              ta= (dla+0.5+jz)/td
              tb= (dla+0.5-jz)/td
              if (ta.gt.1.d-15) then
                calf(knt)= dsqrt(ta) 
                malf(knt)= int(nint(jz-0.5))
              else
                calf(knt)=0.0d0
                malf(knt)= la+1
              end if
              if (tb.gt.1.d-15) then
                cbet(knt)= dsqrt(tb) 
                mbet(knt)= int(nint(jz+0.5))
              else
                cbet(knt)=0.0d0
                mbet(knt)= la+1
              end if
              knt=knt+1
          end do                
          if (la.ne.0) then
           kappa=la  
           j=dble(la)-0.5d0
           jmax=2*la
           do n3=0,jmax
              jz=dble(n3)-j
              n_prj(knt)=1
              lprj(knt)=la
              kprj(knt)=int(kappa)
              td=dble(la+la+1)
              dla=dble(la)
              ta= (dla+0.5-jz)/td
              tb= (dla+0.5+jz)/td
              if (ta.gt.1.d-15) then
                calf(knt)= -dsqrt(ta) 
                malf(knt)=int(nint(jz-0.5))
              else
                calf(knt)=0.0d0
                malf(knt)=la+1
              end if
              if (tb.gt.1.d-15) then
                cbet(knt)= dsqrt(tb) 
                mbet(knt)=int(nint(jz+0.5))
              else
                cbet(knt)=0.0d0
                mbet(knt)=la+1
              end if
              knt=knt+1
            end do                
          end if
      end do
      if (nprj.ne.knt) then
         write(6,*)"Number of Projectors is wrong!!!!"
         CALL flush(6)
      end if       
      return
      end
cccccccccccccccccccccccc






  
