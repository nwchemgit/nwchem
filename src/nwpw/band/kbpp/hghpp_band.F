*
* $Id$
*

*     **************************************
*     *                                    *
*     *           hghpp_band               *
*     *                                    *
*     **************************************

      logical function hghpp_band(oprint_in,version,
     >                  psp_filename,formatted_filename,
     >                  ngrid,unita,
     >                  nbrillioun,kvectors)
      implicit none

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

      logical          oprint_in,filter
      integer          version
      character*50     psp_filename,formatted_filename
      integer          ngrid(3)
      double precision unita(3,3),ecut,wcut,zvec(3)
      integer nbrillioun
      real*8  kvectors(3,*)

*     **** local variables ****
      character*255 full_filename

      logical value,mprint,hprint,oprint
       
      integer taskid,MASTER,msglen
      parameter (MASTER=0)

*     **** 1d pseudopotential data ****
      character*80 comment

      character*2 atom
      integer Zion
      double precision zv,amass

      double precision rloc,C1,C2,C3,C4
      double precision r(0:3),H(3,0:3),K(3,0:3)
      double precision Gijl(9*4),Gtmp(3,3,0:3)
      double precision Kijl(9*4),Ktmp(3,3,0:3)
      double precision rcore

      integer lmax,locp

      integer n_prj_indx,l_prj_indx,m_prj_indx
      integer n_prj_hndl,l_prj_hndl,m_prj_hndl

*     ***** ngrid data *****
      integer Ylm_indx,Slm_indx,vl_indx,vnl_indx,G_indx,vnlso_indx
      integer Ylm_hndl,Slm_hndl,vl_hndl,vnl_hndl,G_hndl,vnlso_hndl

*     **** other variables ****
      double precision unitg(3,3)
      integer nsize,i,j,n,l,m,psp_type,indx,nb
      integer nfft1,nfft2,nfft3
      integer nmax,nprj

*     **** external functions ****
      logical  control_print,control_kbpp_filter
      external control_print,control_kbpp_filter
      real*8   control_ecut,control_wcut
      external control_ecut,control_wcut


      !call Parallel_init()
      call Parallel_taskid(taskid)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)
      oprint = (oprint_in.or.hprint)


      value = .false.
  
*     ***** read in pseudopotential data ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         open(unit=11,file=full_filename(1:l),
     >                status='old',form='formatted')
         read(11,*) psp_type
         read(11,'(A2)') atom
         read(11,*) Zion
         read(11,*) lmax
         read(11,*) rloc,C1,C2,C3,C4

         if (lmax.ge.0) then
            read(11,*) r(0),H(1,0),H(2,0),H(3,0)
            do l=1,lmax
              read(11,*) r(l),H(1,l),H(2,l),H(3,l)
              read(11,*)      K(1,l),K(2,l),K(3,l)
            end do
         end if
         read(11,'(A)') comment
         close(11)
         !**** determine nmax ****
         nmax = 0
         do l=0,lmax
            do i=1,3 
              if (  (dabs(H(i,l)).gt.1.0d-14) .and. 
     >            (i.gt.nmax)) nmax = i
            end do
         end do
      end if

      psp_type=1
      msglen = 1
      call BRDCST(9+MSGINT,psp_type,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,Zion,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,lmax,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,nmax,mitob(msglen),MASTER)

      call BRDCST(9+MSGDBL,rloc,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,C1,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,C2,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,C3,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,C4,mdtob(msglen),MASTER)

      msglen = 4
      call BRDCST(9+MSGDBL,r,mdtob(msglen),MASTER)

      msglen = 12
      call BRDCST(9+MSGDBL,H,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,K,mdtob(msglen),MASTER)



*     **** set the maximum angular momentum ****

*     **** set the local potential ****
      locp = lmax+1

      
      !**** determine nprj ****
      nprj= 0
      do l=0,lmax
            !write(*,*) "???H :",l,(H(i,l),i=1,3)
      do i=1,3 
        if ( dabs(H(i,l)).gt.1.0d-14) nprj = nprj + (2*l+1)
      end do
      end do
      !write(*,*) "???nprj:", nprj
        


*     **** set the projector coeficients ****
*  HGH does not give coefficients for i .ne. j
*  we must use their formula
*******************************************
      call dcopy(9*4,0.0d0,0,Gtmp,1)
      call dcopy(9*4,0.0d0,0,Gijl,1)
      call dcopy(9*4,0.0d0,0,Ktmp,1)
      call dcopy(9*4,0.0d0,0,Kijl,1)
      Gtmp(1,1,0) = H(1,0)
      Gtmp(2,2,0) = H(2,0)
      Gtmp(3,3,0) = H(3,0)
      Gtmp(1,2,0) = -0.5d0*dsqrt(3.0d0/5.0d0)*Gtmp(2,2,0)
      Gtmp(2,1,0) = Gtmp(1,2,0)
      Gtmp(1,3,0) =  0.5d0*dsqrt(5.0d0/21.0d0)*Gtmp(3,3,0)
      Gtmp(3,1,0) = Gtmp(1,3,0)
      Gtmp(2,3,0) = -0.5d0*dsqrt(100.0d0/63.0d0)*Gtmp(3,3,0)
      Gtmp(3,2,0) = Gtmp(2,3,0)


      Gtmp(1,1,1) = H(1,1)
      Gtmp(2,2,1) = H(2,1)
      Gtmp(3,3,1) = H(3,1)
      Gtmp(1,2,1) = -0.5d0*dsqrt(5.0d0/7.0d0)*Gtmp(2,2,1)
      Gtmp(2,1,1) = Gtmp(1,2,1)
      Gtmp(1,3,1) =  (1.0d0/6.0d0)*dsqrt(35.0d0/11.0d0)*Gtmp(3,3,1)
      Gtmp(3,1,1) = Gtmp(1,3,1)
      Gtmp(2,3,1) = -(14.0d0/6.0d0)*dsqrt(1.0d0/11.0d0)*Gtmp(3,3,1)
      Gtmp(3,2,1) = Gtmp(2,3,1)

      Gtmp(1,1,2) = H(1,2)
      Gtmp(2,2,2) = H(2,2)
      Gtmp(3,3,2) = H(3,2)
      Gtmp(1,2,2) = -0.5d0*dsqrt(7.0d0/9.0d0)*Gtmp(2,2,2)
      Gtmp(2,1,2) = Gtmp(1,2,2)
      Gtmp(1,3,2) =  0.5d0*dsqrt(63.0d0/143.0d0)*Gtmp(3,3,2)
      Gtmp(3,1,2) = Gtmp(1,3,2)
      Gtmp(2,3,2) = -0.5d0*18.0d0*dsqrt(1.0d0/143.0d0)*Gtmp(3,3,2)
      Gtmp(3,2,2) = Gtmp(2,3,2)

      Gtmp(1,1,3) = H(1,3)
      Gtmp(2,2,3) = H(2,3)
      Gtmp(3,3,3) = H(3,3)

      Ktmp(1,1,0) = 0.0d0
      Ktmp(2,2,0) = 0.0d0
      Ktmp(3,3,0) = 0.0d0
      Ktmp(1,2,0) = 0.0d0
      Ktmp(2,1,0) = 0.0d0
      Ktmp(1,3,0) = 0.0d0
      Ktmp(3,1,0) = 0.0d0
      Ktmp(2,3,0) = 0.0d0
      Ktmp(3,2,0) = 0.0d0


      Ktmp(1,1,1) = K(1,1)
      Ktmp(2,2,1) = K(2,1)
      Ktmp(3,3,1) = K(3,1)
      Ktmp(1,2,1) = -0.5d0*dsqrt(5.0d0/7.0d0)*Ktmp(2,2,1)
      Ktmp(2,1,1) = Ktmp(1,2,1)
      Ktmp(1,3,1) =  (1.0d0/6.0d0)*dsqrt(35.0d0/11.0d0)*Ktmp(3,3,1)
      Ktmp(3,1,1) = Ktmp(1,3,1)
      Ktmp(2,3,1) = -(14.0d0/6.0d0)*dsqrt(1.0d0/11.0d0)*Ktmp(3,3,1)
      Ktmp(3,2,1) = Ktmp(2,3,1)

      Ktmp(1,1,2) = K(1,2)
      Ktmp(2,2,2) = K(2,2)
      Ktmp(3,3,2) = K(3,2)
      Ktmp(1,2,2) = -0.5d0*dsqrt(7.0d0/9.0d0)*Ktmp(2,2,2)
      Ktmp(2,1,2) = Ktmp(1,2,2)
      Ktmp(1,3,2) =  0.5d0*dsqrt(63.0d0/143.0d0)*Ktmp(3,3,2)
      Ktmp(3,1,2) = Ktmp(1,3,2)
      Ktmp(2,3,2) = -0.5d0*18.0d0*dsqrt(1.0d0/143.0d0)*Ktmp(3,3,2)
      Ktmp(3,2,2) = Ktmp(2,3,2)

      Ktmp(1,1,3) = K(1,3)
      Ktmp(2,2,3) = K(2,3)
      Ktmp(3,3,3) = K(3,3)

      do l=0,lmax
        do i=1,nmax
        do j=1,nmax
           Gijl(i+(j-1)*nmax+l*nmax*nmax) = Gtmp(i,j,l)
           Kijl(i+(j-1)*nmax+l*nmax*nmax) = Ktmp(i,j,l)
        end do
        end do
      end do


*     **** allocate vl, vnl, G ****
      nsize = ngrid(1)*ngrid(2)*ngrid(3)
      value = BA_alloc_get(mt_dbl,nsize,
     >                        'vl',vl_hndl,vl_indx)
      value = value.and.BA_alloc_get(mt_dbl,nsize*(nprj),
     >                        'vnl',vnl_hndl, vnl_indx)
      value = value.and.BA_alloc_get(mt_dcpl,nsize*(nprj),
     >                        'vnlso',vnlso_hndl, vnlso_indx)
      value = value.and.BA_alloc_get(mt_dbl,nsize,
     >                        'Slm',Slm_hndl, Slm_indx)
      value = value.and.BA_alloc_get(mt_dcpl,nsize,
     >                        'Ylm',Ylm_hndl, Ylm_indx)

      value = value.and.BA_alloc_get(mt_dbl,nsize*(3),
     >                        'G',G_hndl, G_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'n_prj', n_prj_hndl, n_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'l_prj', l_prj_hndl, l_prj_indx)
      value = value.and.BA_alloc_get(mt_int,nprj,
     >                        'm_prj', m_prj_hndl, m_prj_indx)
      if(.not.value)
     >    call errquit('hghppv1: out of heap memory', 0, MA_ERR)


      !**** determine n_prj, l_prj, and m_prj arrays ****
      indx = 0
      do l=0,lmax
      do i=1,3 
        if ((H(i,l).ne.0.0d0)) then
          do m=-l,l
            int_mb(n_prj_indx+indx) = i
            int_mb(l_prj_indx+indx) = l
            int_mb(m_prj_indx+indx) = m
            indx = indx + 1
          end do 
        end if
      end do
      end do

      filter = control_kbpp_filter()
      ecut   = control_ecut()
      wcut   = control_wcut()

*     **** preparation of constants ****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      call setup_kbpp_band(nfft1,nfft2,nfft3,unita,unitg,
     >                     dbl_mb(G_indx))
      zv = dble(Zion)

      call HGH_band_local(version,
     >               zv,rloc,C1,C2,C3,C4,
     >               nfft1,nfft2,nfft3,
     >               dbl_mb(G_indx),
     >               dbl_mb(vl_indx))

      if (filter) then
         zvec(1) = 0.0d0
         zvec(2) = 0.0d0
         zvec(3) = 0.0d0
         call kbpp_band_filter(nfft1,nfft2,nfft3,
     >                    zvec,
     >                    dbl_mb(G_indx),
     >                    ecut,
     >                    dbl_mb(vl_indx))
      end if



      if ((taskid.eq.MASTER).and.(oprint)) then
      write(*,*) "     ********************************************"
      write(*,*) "     *                                          *"
      write(*,*) "     *  HGHPP_BAND - Pseudopotential Formatter  *"
      write(*,*) "     *                                          *"
      write(*,*) "     *      version last updated 12/31/06       *"
      write(*,*) "     *                                          *"
      write(*,*) "     *       developed by Eric J. Bylaska       *"
      write(*,*) "     *                                          *"
      write(*,*) "     ********************************************"
      call nwpw_message(1)
      write(*,*)
      write(*,*) "Pseudpotential Data"
      write(*,*) "-------------------"
      write(*,*) "  atom     :",atom
      write(*,*) "  charge   :",Zion
      write(*,*) "  highest angular component used :",lmax
      write(*,*) "  highest radial  component used :",nmax
      write(*,*) "  number of non-local projectors :",nprj
      write(*,111) "   projector cutoffs: ",(r(i), i=0,lmax)
      write(*,*)
      write(*,111) " local psp cutoff       : ",rloc
      write(*,111) " local psp coefficients : ",C1,C2,C3,C4
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
           write(luout,*)
           write(luout,*) "Generated formatted_filename: ",
     >          full_filename(1:l)
           if (filter) write(luout,*) "- filtering pseudopotential -"
         end if
         call openfile(2,full_filename,l,'w',l)     
         call cwrite(2,comment,80)
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
         call dwrite(2,r,lmax+1)

         
         call iwrite(2,nprj,1)
         if (nprj.gt.0) then
          call iwrite(2,int_mb(n_prj_indx),nprj)
          call iwrite(2,int_mb(l_prj_indx),nprj)
          call iwrite(2,int_mb(m_prj_indx),nprj)
          call dwrite(2,Gijl,(nmax*nmax*(lmax+1)))
          call dwrite(2,Kijl,(nmax*nmax*(lmax+1)))
         end if

         rcore = 0.0d0
         call dwrite(2,rcore,1)

         call iwrite(2,nbrillioun,1)
         call dwrite(2,kvectors,3*nbrillioun)

         call dwrite(2,dbl_mb(vl_indx),nsize)
      end if


*     **** generate and write out nonlocal psp ****
      do nb=1,nbrillioun
        do i=1,nprj
           n=int_mb(n_prj_indx+i-1)
           l=int_mb(l_prj_indx+i-1)
           m=int_mb(m_prj_indx+i-1)
           call Tesseral_band(l,m,kvectors(1,nb),
     >                   nfft1,nfft2,nfft3,
     >                   dbl_mb(G_indx),
     >                   dbl_mb(Slm_indx))
           call HGH_band_nonlocal(
     >                       n,l,kvectors(1,nb),
     >                       r(l),
     >                       nfft1,nfft2,nfft3,
     >                       dbl_mb(G_indx),
     >                       dbl_mb(vnl_indx+(i-1)*nsize))
           if ((taskid.eq.MASTER).and.(oprint))
     >        write(*,*) "creating projector:",n,l,m," nb=",nb

           do j=1,nsize
             dbl_mb(vnl_indx+(i-1)*nsize+j-1)
     >        = dbl_mb(vnl_indx+(i-1)*nsize+j-1)*dbl_mb(Slm_indx+j-1)
           end do
           if (filter)
     >        call kbpp_band_filter(nfft1,nfft2,nfft3,
     >                    kvectors(1,nb),
     >                    dbl_mb(G_indx),
     >                    wcut,
     >                    dbl_mb(vnl_indx+(i-1)*nsize))

        end do

        if ((taskid.eq.MASTER).and.(nprj.gt.0)) then
           call dwrite(2,dbl_mb(vnl_indx),nsize*nprj)
        end if
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*     **** generate and write out spin_orbit psp **************************
*     spin orbit is Hermitian complex. it saves us nothing to 
*     use real harmonics
***************************************************************************

      do nb=1,nbrillioun
        do i=1,nprj
           n=int_mb(n_prj_indx+i-1)
           l=int_mb(l_prj_indx+i-1)
           m=int_mb(m_prj_indx+i-1)
           call spherical_band(l,m,kvectors(1,nb),
     >                   nfft1,nfft2,nfft3,
     >                   dbl_mb(G_indx),
     >                   dcpl_mb(Ylm_indx))
           call HGH_band_nonlocal(
     >                       n,l,kvectors(1,nb),
     >                       r(l),
     >                       nfft1,nfft2,nfft3,
     >                       dbl_mb(G_indx),
     >                       dbl_mb(vnl_indx+(i-1)*nsize))
           if ((taskid.eq.MASTER).and.(oprint))
     >        write(*,*) "creating projector:",n,l,m," nb=",nb

           do j=1,nsize
             dcpl_mb(vnlso_indx+(i-1)*nsize+j-1)
     >        = dcmplx(dbl_mb(vnl_indx+(i-1)*nsize+j-1),0.0d0)*
     >             dcpl_mb(Ylm_indx+j-1)
           end do
        end do

        if ((taskid.eq.MASTER).and.(nprj.gt.0)) then
           call dwrite(2,dcpl_mb(vnlso_indx),2*nsize*nprj)
        end if
      end do

      if (taskid.eq.MASTER) call closefile(2)
     
     
*     **** free heap space ****
      value = BA_free_heap(vl_hndl)
      value = value.and.BA_free_heap(vnl_hndl)
      value = value.and.BA_free_heap(vnlso_hndl)
      value = value.and.BA_free_heap(Ylm_hndl)
      value = value.and.BA_free_heap(Slm_hndl)
      value = value.and.BA_free_heap(G_hndl)
      value = value.and.BA_free_heap(n_prj_hndl)
      value = value.and.BA_free_heap(l_prj_hndl)
      value = value.and.BA_free_heap(m_prj_hndl)
      if(.not.value)
     >  call errquit('hghpp_band: deallocating heap memory', 0, MA_ERR)

      
      if ((taskid.eq.MASTER).and.(oprint)) call nwpw_message(4)

      hghpp_band = value

      return
 9999 call errquit('hghpp_band:Error reading psp_filename',0, DISK_ERR)
      END




