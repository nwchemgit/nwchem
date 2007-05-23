* $Id: rkbpp_band.F,v 1.11 2007-05-23 00:06:33 d3p708 Exp $
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
      character*255 pspso_filename
      logical value,mprint,hprint,oprint
      integer lmax,locp
      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*2 atom
      character*80 comment
      double precision zv,amass
      integer lmax0,lmmax,lmax1,locp1,nprj,ierr
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
c      value=MA_set_auto_verify(.true.)
cccccccccccccccccccccccccccccccccccccc
c      write(6,*)"DIE KATZE a"
c      CALL FLUSH(6)

      value = .false.
  
*     ***** read in pseudopotential data ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
         l = index(full_filename,' ') - 1
         write(*,*)"psp FILENAME IS ",psp_filename
         call flush(6)
         write(*,*)"full FILENAME IS ",full_filename
         call flush(6)
         pspso_filename=full_filename(1:l)//'.so' 
         write(*,*)"pspso FILENAME IS ",pspso_filename
         call flush(6)
        open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')
        read(11,12,ERR=9999,END=9999)pspx,atom
        read(11,*,ERR=9999,END=9999) zv,amass,lmax0,lmax1,
     >    locp1,rlocal1
        read(11,*,ERR=9999,END=9999) (rc(i),i=0,lmax0)
        read(11,*,ERR=9999,END=9999) nrho,drho
        read(11,'(A)',ERR=9999,END=9999) comment
c        write(*,12)pspx,atom
c        write(*,*) zv,amass,lmax0,lmax1,
c     >    locp1,rlocal1
c        write(*,*) (rc(i),i=0,lmax0)
c        write(*,*) nrho,drho
c        write(*,'(A)') comment
12      format(i1,A2)
        call flush(6)
      end if

c      write(6,*)"DIE KATZE B"
c      CALL FLUSH(6)
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
      lmax=lmax0

*     **** set the local potential ****
      if (locp.eq.-1)   locp = locp1
      if (locp.gt.lmax) locp = lmax
      if (locp.lt.0)    locp = lmax
      
      
      nprj=0
      do kp=0,lmax
        nprj=nprj+kp+1
      end do  
      nprj=nprj*4
*     **** allocate rho, vp, and wp ****
      value = MA_alloc_get(mt_dbl,nrho,
     >                        'rho',rho_hndl,rho_indx)
      value = MA_alloc_get(mt_dbl,nrho*(lmax+2),
     >                        'vp',vp_hndl,vp_indx)
      value = MA_alloc_get(mt_dbl,nrho*(lmax+2),
     >                        'vpso',vpso_hndl,vpso_indx)
      value = MA_alloc_get(mt_dbl,nrho*(2*lmax+3),
     >                        'wp', wp_hndl,wp_indx)
      value = MA_alloc_get(mt_dbl,2*nrho,
     >                        'sc', sc_r_hndl, sc_r_indx)
     

      if (taskid.eq.MASTER) then
       call read_rvpwp_band(11,nrho,lmax,dbl_mb(rho_indx),
     >                         dbl_mb(vp_indx),
     >                         dbl_mb(wp_indx))
       call read_semicore_band(11,isemicore,rcore,nrho,
     > dbl_mb(sc_r_indx))
       close(11)
c      value=MA_verify_allocator_stuff()
       l = index(pspso_filename,' ') - 1
       open(unit=11,file=pspso_filename(1:l),
     >             status='old',form='formatted')
       call read_rkbpp_so(11,nrho,lmax,dbl_mb(vpso_indx))
       close(11)
      end if

      msglen = nrho
      call BRDCST(9+MSGDBL,dbl_mb(rho_indx),mdtob(msglen),MASTER)
      msglen = nrho*(lmax+1)
      call BRDCST(9+MSGDBL,dbl_mb(vp_indx),mdtob(msglen),MASTER)
      msglen = nrho*(2*lmax+1)
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
      value = MA_alloc_get(mt_dcpl,nsize*2*nprj,
     >                        'vnl',vnl_hndl, vnl_indx)
      value = MA_alloc_get(mt_dbl,nprj,
     >                        'vnrm', vnrm_hndl, vnrm_indx)
      value = MA_alloc_get(mt_dbl,nsize*3,
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
     >                         c_bet_indx)
cccccccccc assign the relativistic states ccccccccccccccccccc
      call assign_rel_states(nprj,lmax,int_mb(n_prj_indx),
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
      call rkbpp_vlocal(lmax,locp,nsize,nrho,
     >       dbl_mb(vp_indx),dbl_mb(vpso_indx),
     >       dbl_mb(vnrm_indx),dbl_mb(vl_indx),
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

      nb=1

      call rkbpp_vnonlocal(lmax,nrho,locp,
     >  nfft1,nfft2,nfft3,nprj,
     >  dbl_mb(G_indx),kvectors(1,nb),
     >  dbl_mb(vp_indx),dbl_mb(vpso_indx),
     >  dbl_mb(wp_indx),dcpl_mb(vnl_indx),
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
c        if (mprint) then
           write(*,*)
           write(*,*) "Generated formatted_filename: ",
     >            full_filename(1:l)
c        end if
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
        nmax = 1
        call iwrite(2,nmax,1)
        call dwrite(2,rc,lmax+1)

        call iwrite(2,nprj,1)
        if (nprj.gt.0) then
         call iwrite(2,int_mb(n_prj_indx),nprj)
         call iwrite(2,int_mb(l_prj_indx),nprj)
         call iwrite(2,int_mb(k_prj_indx),nprj)
         call dwrite(2,dbl_mb(vnrm_indx),nprj)
        end if  
        call dwrite(2,rcore,1)
        call iwrite(2,nbrillioun,1)
        call dwrite(2,kvectors,3*nbrillioun)
        call dwrite(2,dbl_mb(vl_indx),nsize)
        call dwrite(2,dcpl_mb(vnl_indx),nsize*nprj*2)
      end if


      do nb=2,nbrillioun

        if ((oprint).and.(taskid.eq.MASTER))
     >      write(*,*) "generating brillioun #",nb

        call rkbpp_vnonlocal(lmax,nrho,locp,
     >  nfft1,nfft2,nfft3,nprj,
     >  dbl_mb(G_indx),kvectors(1,nb),
     >  dbl_mb(vp_indx),dbl_mb(vpso_indx),
     >  dbl_mb(wp_indx),dcpl_mb(vnl_indx),
     >  dbl_mb(cs_indx),dbl_mb(sn_indx),dbl_mb(f_indx),
     >  dbl_mb(rho_indx),int_mb(l_prj_indx),
     >  int_mb(k_prj_indx),
     >  int_mb(m_alf_indx),int_mb(m_bet_indx),
     >  dbl_mb(c_alf_indx),dbl_mb(c_bet_indx),
     >  drho,ierr)

        if (taskid.eq.MASTER) then
          call dwrite(2,dbl_mb(vnl_indx),nsize*nprj*4)
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
      integer nrho,lmax,npot,nprj
      double precision rho(nrho),rx
      double precision vp(nrho,0:(lmax+1))
      double precision wp(nrho,2*lmax+3)
      integer i,j
      npot= lmax+1
      nprj= lmax*2+3
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
      subroutine read_rkbpp_so(unit,nrho,lmax,vso)
      implicit none
#include "errquit.fh"
      integer i,j,unit,nrho,lmax
      double precision vso(nrho,0:lmax)
      double precision rhox
ccccccccccccccccccc
      call dcopy((lmax+1)*nrho,0.0d0,0,vso,1)
      do i=1,nrho
        read(unit,*,ERR=9999,END=9999) rhox,(vso(i,j),j=1,(lmax+1))
      end do
      return
 9999 call errquit('Error reading psp_so_filename',0, DISK_ERR)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Assing the value of l,m,kappa, c_alf, c_beta, etc.
c  for the Dirac Atomic states.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  note m_alfa/m_beta is set equal to lp+1 for states that are zero!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine assign_rel_states(nprj,lmax,n_prj,l_prj,k_prj,
     >  malfa,mbeta,calfa,cbeta)
      implicit none
      integer nprj,lmax,n_prj(*),l_prj(*),k_prj(*),malfa(*),mbeta(*)
      real*8 calfa(*),cbeta(*)
ccccccccccccc local variables cccccccccccccccccc
      integer ip,lp,ik,n3,ns
      real*8 dl,dk,td,ta,tb,jz,j
cccccccccccccccccccccccccccccccccccccccccccccccc      
      ip=1
      do 100 lp=0,lmax
         dk= -dble(lp)-1.0d0
         ik= -lp-1
         dl= dble(lp)
         j=dl+0.5d0
         ns=int(j+j+1.0d0)
         do n3=0,ns
            jz= dble(n3)-j
            n_prj(ip)=1
            l_prj(ip)=lp
            k_prj(ip)=ik
            td=lp+lp+1.0d0
            ta= (dl+0.5d0+jz)/td
            tb= (dl+0.5d0-jz)/td
            if (ta.gt.1.d-15) then
                 calfa(ip)=dsqrt(ta)
                 malfa(ip)=int(nint(jz-0.5d0))
            else
                 calfa(ip)=0.0d0
                 malfa(ip)=lp+1
            end if
            if (tb.gt.1.d-15) then
                 cbeta(ip)=dsqrt(tb)
                 mbeta(ip)=int(nint(jz+0.5d0))
            else
                 cbeta(ip)=0.0d0
                 mbeta(ip)=lp+1
            end if                
            ip=ip+1
         end do
         if (lp.eq.0) goto 100
         dk=dble(lp)
         ik=lp
         j=dl-0.5d0
         ns=int(j+j+1.0d0)
         do n3=0,ns
            jz=dble(n3)-j
            n_prj(ip)=1
            l_prj(ip)=lp
            k_prj(ip)=ik
            td=lp+lp+1.0d0
            ta= (dl+0.5d0-jz)/td
            tb= (dl+0.5d0+jz)/td
            if (ta.gt.1.d-15) then
                 calfa(ip)= -dsqrt(ta)
                 malfa(ip)= int(nint(jz-0.5d0))
            else
                 calfa(ip)=0.0d0
                 malfa(ip)=lp+1
            end if
            if (tb.gt.1.d-15) then
                 cbeta(ip)=dsqrt(tb)
                 mbeta(ip)=int(nint(jz+0.5d0))
            else
                 cbeta(ip)=0.0d0
                 mbeta(ip)=lp+1
            end if                
            ip=ip+1
         end do
100   continue
      nprj=ip
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rkbpp_vlocal(lmax,locl,ngsize,nrho,
     >       vp,vpso,vnlnrm,vl,
     >       f,sinqr,cosqr,rho_sc_k,rho_sc_r,
     >       gvec,semicore,zv,drho,ierr)
      implicit none
      logical semicore
      integer lmax
      integer nrho
      integer ngsize
      integer locl,ierr
      double precision rho(nrho)
      double precision f(nrho)
      double precision sinqr(nrho)
      double precision cosqr(nrho)
      double precision vp(nrho,0:lmax)
      double precision vpso(nrho,0:lmax)
      double precision vl(ngsize)
      double precision gvec(3,ngsize)
      double precision wp(nrho,2*lmax+1)
      double precision vnlnrm(2*lmax+1)
      double precision rho_sc_k(ngsize,4)
      double precision rho_sc_r(nrho,2)
      double precision zv
      double precision drho
cccccccc locals ccccccccccccccccccccccccccccccccccccccc       
      integer nproc,taskid,lp,ip,k,i,task_count
      double precision fourpi,vv,a,gx,gy,gz,q
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
      ip=1
      do 100 lp=0,lmax
         if (lp.eq.locl) then
           vnlnrm(ip)=0.0d0
           ip=ip+1 
           goto 100
         else   
           do k=1,nrho
             vv=vp(k,lp)-dble(lp+1)*vpso(k,lp)
             f(k)=vv*wp(k,ip)*wp(k,ip)            
           end do
           a=simp(nrho,f,drho) 
           vnlnrm(ip)=1.0/a
           ip=ip+1
           do k=1,nrho
             vv=vp(k,lp)+dble(lp)*vpso(k,lp)
             f(k)=vv*wp(k,ip)*wp(k,ip)            
           end do
           a=simp(nrho,f,drho) 
           vnlnrm(ip)=1.0/a
           ip=ip+1
         end if
100   continue
      call dcopy(ngsize,0.0d0,0,vl,1)
      call dcopy(4*ngsize,0.0d0,0,rho_sc_k,1)  
      task_count= -1
      do 700 i=1,ngsize
         task_count=task_count+1
         if ( (nproc*(task_count/nproc)).ne.taskid) goto 700
         if ( i.eq.1) goto 700
         gx=gvec(1,i)
         gy=gvec(2,i)
         gz=gvec(3,i)
         q=dsqrt(gx*gx+gy*gy+gz*gz)
         gx=gx/q
         gy=gy/q
         gz=gz/q
         do k=1,nrho
            cosqr(k)=dcos(q*rho(k))
            sinqr(k)=dsin(q*rho(k))
         end do
         do k=1,nrho
            f(i)=rho(i)*vp(i,locl)*sinqr(i)
         end do
         vl(i)=simp(nrho,f,drho)*fourpi/q-zv*fourpi*cosqr(nrho)/q/q
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
      call D3dB_Vector_SumAll(4*ngsize,rho_sc_k)
      call D3dB_Vector_SumAll(ngsize,VL)
cccccccccc g=0
      do k=1,nrho
         f(k)=vp(k,locl)*rho(k)*rho(k)
         vp(k,locl)=0.0d0
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
     >  vp,vpso,wp,vnl,
     >  cosqr,sinqr,f,rho,l_prj,k_prj,malf_prj,mbet_prj,calf_prj,
     >  cbet_prj,drho,ierr)
      implicit none
#include "mafdecls.fh"
      integer lmax,locl,nrho,ierr
      integer ngx,ngy,ngz,nprj
      integer malf_prj(*),mbet_prj(*),l_prj(*),k_prj(*)
      real*8 cosqr(nrho),sinqr(nrho),f(nrho),rho(nrho)
      real*8 vp(nrho,0:(lmax+1)),vpso(nrho,0:(lmax+1))
      real*8 wp(nrho,0:(2*lmax+2))  
      real*8 vnl(ngx,ngy,ngz,nprj,2)
      real*8 calf_prj(nprj)
      real*8 cbet_prj(nprj)
      real*8 drho,gvec(ngx,ngy,ngz,3),kvec(3)
cccccccc locals
      integer taskid,nproc,tcount,k
      integer igx,igy,igz
      integer malfa,mbeta,kappa,l,ippot,key
      real*8 pi4,q,gx,gy,gz,a,d,vv,qr,pi
      real*8 calfa,cbeta,scof(0:8)
      real*8 simp,sjx
      logical value
      complex*16 sphHarm
      external sphHarm
      external simp,sjx  
cccccccccccccccccccccccccccccccccccccccccccccccc
      value=MA_set_auto_verify(.true.)
      call Parallel_taskid(taskid)
      call Parallel_np(nproc)
      pi=4.0d0*datan(1.0d0)
      pi4=4.0d0*pi
      scof(3)=dsqrt(pi4*105.0d0)
      scof(2)=dsqrt(pi4*15.0d0)
      scof(1)=dsqrt(pi4*3.0d0)
      scof(0)=dsqrt(pi4)
      call dcopy(nprj*ngx*ngy*ngz*4,0.0d0,0,vnl,1)
      value=MA_verify_allocator_stuff()
      write(*,*)"Anl 1"
      CALL FLUSH(6)
      tcount=-1
      do 200 igz=1,ngz
      do 200 igy=1,ngy
      do 200 igx=1,ngx
C         value=MA_verify_allocator_stuff()
C       write(*,*)"Anl TR1",igx,igy,igz
C       CALL FLUSH(6)
         tcount=tcount+1
         key=tcount-nproc*(tcount/nproc)
         if (key.ne.taskid) goto 200
C       write(*,*)"gvec(1)= ",gvec(igx,igy,igz,1)," kvec(1)= ",kvec(1)
C       call flush(6)
C       value=MA_verify_allocator_stuff()
         gx=gvec(igx,igy,igz,1)+kvec(1)
C       value=MA_verify_allocator_stuff()
C       write(*,*)"Anl TR2x"
C       write(*,*)gvec(igx,igy,igz,2)," kvec(2)= ",kvec(2)
C       call flush(6)
C       value=MA_verify_allocator_stuff()
         gy=gvec(igx,igy,igz,2)+kvec(2)
C       value=MA_verify_allocator_stuff()
C       write(*,*)"Anl TR2y"
C       write(*,*)gvec(igx,igy,igz,3)," kvec(3)= ",kvec(3)
C       call flush(6)
C       value=MA_verify_allocator_stuff()
         gz=gvec(igx,igy,igz,3)+kvec(3)
C       value=MA_verify_allocator_stuff()
C       write(*,*)"Anl TR2z"
C       CALL FLUSH(6)
         q=dsqrt(gx*gx+gy*gy+gz*gz)
         if (dabs(q).lt.1.e-9) then 
C      value=MA_verify_allocator_stuff()
C      write(*,*)"Anl 0XX"
C      CALL FLUSH(6)
           if (locl.ne.0) then
             do k=1,nrho
               f(k)=rho(k)*wp(k,0)*vp(k,0)
             end do
             vnl(igx,igy,igz,1,1)=simp(nrho,f,drho)*scof(0)
             do k=1,nrho
               f(k)=rho(k)*wp(k,1)*vp(k,0)
             end do
             vnl(igx,igy,igz,1,2)=simp(nrho,f,drho)*scof(0)
           end if    
C      value=MA_verify_allocator_stuff()
C      write(*,*)"Anl 0XY"
C      CALL FLUSH(6)
           goto 200
         end if
         gx=gx/q
         gy=gy/q
         gz=gz/q
C      value=MA_verify_allocator_stuff()
C      write(*,*)"Anl TR3"
C      CALL FLUSH(6)
         do k=1,nrho
            cosqr(k)=dcos(q*rho(k))
            sinqr(k)=dsin(q*rho(k))
         end do
C      value=MA_verify_allocator_stuff()
C      write(*,*)"Anl TR4"
C      CALL FLUSH(6)
C      value=MA_verify_allocator_stuff()
C      write(*,*)"Anl 0XZ1"
C      CALL FLUSH(6)
         do ippot=1,nprj
            l=l_prj(ippot)
            if (l.eq.locl) then
               vnl(igx,igy,igz,ippot,1)=0.0d0
               vnl(igx,igy,igz,ippot,2)=0.0d0
               cycle
            end if
            kappa=k_prj(ippot)
            malfa=malf_prj(ippot)
            mbeta=mbet_prj(ippot)
            calfa=calf_prj(ippot)
            cbeta=cbet_prj(ippot)
            f(1)=0.0d0
            do k=1,nrho
               vv=vp(k,l)+dble(kappa)*vpso(k,ippot-1)
    	       vv=vv*wp(k,ippot-1)
               qr=q*rho(k)
               a=sjx(qr,l,sinqr(k),cosqr(k))
               f(k)=a*vv
            end do
            d=simp(nrho,f,drho)*scof(l)/q
            vnl(igx,igy,igz,ippot,1)= (d*calfa)*
     >                 sphHarm(l,malfa,gx,gy,gz)
            vnl(igx,igy,igz,ippot,2)= (d*cbeta)*
     >                 sphHarm(l,mbeta,gx,gy,gz)
         end do   
C      value=MA_verify_allocator_stuff()
C      write(*,*)"Anl 0XZ2"
      CALL FLUSH(6)
200   continue 	
      value=MA_verify_allocator_stuff()
      write(*,*)"Anl 2"
      CALL FLUSH(6)
      call C3db_Vector_SumAll(nprj*ngx*ngy*ngz*2,vnl)
      value=MA_verify_allocator_stuff()
      write(*,*)"Anl 3"
      CALL FLUSH(6)
      ierr=0
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      double precision function sjx(qr,L,sn,cs)
      implicit none
      double precision qr,sn,cs
      integer L
      double precision a
      if (L.eq.3) then
        a=sn/qr
        sjx=15.0d0*(a-cs)/qr/qr-6.0d0*a+cs
        return
      end if
      if (L.eq.2) then
        a=sn/qr
        sjx=3.0d0*(a-cs)-cs/qr-sn
        return
      end if
      if (L.eq.1) then
        a=sn/qr
        sjx=a-cs
        return
      end if
      sjx=sn
      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      complex*16 function sphHarm(L,M,x,y,z)
      implicit none
      integer L,M
      real*8 x,y,z,plm,pi,sinth,phi
      complex*16 expm  
      pi=3.14159265358979323846d0
      if (M.gt.L.or.M.lt.(-L)) then
        sphHarm=dcmplx(0.0d0,0.0d0)
        return
      end if
      sinth=dsqrt(x*x+y*y)
      phi=datan2(y,x)
      expm=dcmplx(dcos(phi*dble(M)),dsin(phi*dble(M)))  
      if (L.eq.3) then
        if (m.eq.-3) then
            plm=dsqrt(35.0d0/64.0d0/pi)*sinth*sinth*sinth
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.-2) then
            plm=dsqrt(105.0d0/32.0d0/pi)*sinth*sinth*z
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.-1) then
            plm=dsqrt(21.0d0/64.0d0/pi)*sinth*(z*z*5.0d0-1.0d0)
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.0) then
            plm=dsqrt(7.0d0/16.0d0/pi)*(5.0d0*z*z*z-3.0d0*z)
            sphHarm=plm*expm
            goto 100
        end if
        if (m.eq.1) then
            plm=-dsqrt(21.0d0/64.0d0/pi)*sinth*(z*z*5.0d0-1.0d0)
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.2) then
            plm=dsqrt(105.0d0/32.0d0/pi)*sinth*sinth*z
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.3) then
            plm=-dsqrt(35.0d0/64.0d0/pi)*sinth*sinth*sinth
            sphHarm=plm*expm
            goto 100 
        end if 
      end if 
      if (L.eq.2) then
        if (m.eq.-2) then
            plm=dsqrt(15.0d0/32.0d0/pi)*sinth*sinth
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.-1) then
            plm=dsqrt(15.0d0/8.0d0/pi)*sinth*z
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.0) then
            plm=dsqrt(15.0d0/16.0d0/pi)*(3.0d0*z*z-1.0d0)
            sphHarm=plm*expm
            goto 100
        end if
        if (m.eq.1) then
            plm= -dsqrt(15.0d0/8.0d0/pi)*sinth*z
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.2) then
            plm=dsqrt(15.0d0/32.0d0/pi)*sinth*sinth
            sphHarm=plm*expm
            goto 100 
        end if 
      end if 
      if (L.eq.1) then
        if (m.eq.-1) then
            plm=dsqrt(3.0d0/8.0d0/pi)*sinth
            sphHarm=plm*expm
            goto 100 
        end if 
        if (m.eq.0) then
            plm=dsqrt(3.0d0/4.0d0/pi)*z
            sphHarm=plm*expm
            goto 100
        end if
        if (m.eq.1) then
            plm= -dsqrt(3.0d0/8.0d0/pi)*sinth
            sphHarm=plm*expm
            goto 100 
        end if 
      end if 
      plm= dsqrt(0.75d0/pi)    
      sphHarm=plm*expm
100   return
      end
ccccccccccccccc         









  
