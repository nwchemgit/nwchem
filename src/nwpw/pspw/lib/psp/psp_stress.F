#define TCGMSG

      subroutine psp_stress_init()      
*
* $Id: psp_stress.F,v 1.5 2002-04-24 17:44:49 bylaska Exp $
*
      implicit none

#include "mafdecls.fh"

*     **** common blocks ****
      integer    dvl(2)
      integer    dvnl(2)
      integer    npsp

      integer lmmax_max,lmax_max
      parameter (lmax_max=4,lmmax_max=16)
      common / psp2_block / dvl,dvnl,npsp

      integer npack1,npack0
      logical value

*     **** external functions *****
      integer  ion_nkatm
      external ion_nkatm
 

      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      npsp = ion_nkatm()

      value = MA_alloc_get(mt_dbl,(npsp*npack0),
     >                    'dvl',dvl(2),dvl(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*lmmax_max*npsp*npack1),
     >                    'dvnl',dvnl(2),dvnl(1))

      if (.not. value) call errquit('out of heap memory',0)

      call dcopy(npsp*npack0,          0.0d0,0,dbl_mb(dvl(1)), 1)
      call dcopy(3*lmmax_max*npsp*npack1,0.0d0,0,dbl_mb(dvnl(1)),1)

      return
      end




      subroutine psp_stress_end()      
      implicit none

#include "mafdecls.fh"


*     **** common blocks ****
      integer    dvl(2)
      integer    dvnl(2)
      integer    npsp
      common / psp2_block / dvl,dvnl,npsp

      logical value

      value = MA_free_heap(dvl(2))
      value = MA_free_heap(dvnl(2))

      return
      end


*     ***********************************
*     *									*
*     *	 	   v_local_euv  			*
*     *									*
*     ***********************************

      subroutine v_local_euv(dng,euv)
      implicit none
      complex*16 dng(*)
      real*8     euv(3,3)

#include "mafdecls.fh"



*     *** local variables ***
      integer nfft3d,npack0
      integer ii,ia,u,v,s
      integer exi(2),tmp1(2),tmp2(2)
      integer G(2,3),vl(2)
      logical value
      real*8 elocal,ftmp(3)
      real*8 hm(3,3),Bus(3,3)
      real*8 ss,sum,pi,fourpi

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl

*     **** common blocks ****
      integer    dvl(2)
      integer    dvnl(2)
      integer    npsp
      common / psp2_block / dvl,dvnl,npsp

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      real*8   lattice_omega,lattice_unitg
      external G_indx,ion_nion,ion_katm
      external lattice_omega,lattice_unitg

      call nwpw_timing_start(5)

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      ss   = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do v=1,3
      do u=1,3
         hm(u,v) = ss*lattice_unitg(u,v)
      end do
      end do

*     **** average Kohn-Sham v_local energy ****
      value = MA_push_get(mt_dcpl,npack0,'vl',vl(2),vl(1))
      if (.not. value) call errquit('out of stack memory',0)
      call v_local(dcpl_mb(vl(1)),.false.,dng,ftmp)
      call Pack_cc_dot(0,dng,dcpl_mb(vl(1)),elocal)
      value = MA_pop_stack(vl(2))
      if (.not. value) call errquit('error poping stack memory',0)


      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl, npack0,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dbl, npack0,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'G1',G(2,1),G(1,1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'G2',G(2,2),G(1,2))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('out of stack memory',0)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(G(1,1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(G(1,2)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(G(1,3)))
      call Pack_t_pack(0,dbl_mb(G(1,1)))
      call Pack_t_pack(0,dbl_mb(G(1,2)))
      call Pack_t_pack(0,dbl_mb(G(1,3)))


      call dcopy(9,0.0d0,0,Bus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(0,dcpl_mb(exi(1)))
        
*       **** tmp2(G) = Real[ dconjg(dng(G))*exi(G) ] ****
        call Pack_cct_conjgMul(0,dng,
     >                          dcpl_mb(exi(1)),
     >                          dbl_mb(tmp2(1)))

*       **** tmp2(G) = tmp2(G)*(dvl(G))
        call Pack_tt_Mul(0,dbl_mb(tmp2(1)),
     >                     dbl_mb(dvl(1)),
     >                     dbl_mb(tmp2(1)))

*       **** tmp2(G) = tmp2(G)/G ****
        ss     = 1.0d0/fourpi
        call Pack_t_SMul(0,ss,dbl_mb(vc_indx),dbl_mb(tmp1(1)))
        call Pack_tt_Sqrt(0,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
        call Pack_tt_Mul(0,dbl_mb(tmp1(1)),
     >                     dbl_mb(tmp2(1)),
     >                     dbl_mb(tmp2(1)))
        
*       **** Bus = Bus - Sum(G) tmp2(G)*Gu*Gs ***
        do u=1,3
        do s=u,3
          call Pack_tt_Mul(0,dbl_mb(G(1,u)),
     >                       dbl_mb(G(1,s)),
     >                       dbl_mb(tmp1(1)))
          call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)

          Bus(u,s) = Bus(u,s) - sum
        end do
        end do

      end do
      value =           MA_pop_stack(G(2,3))
      value = value.and.MA_pop_stack(G(2,2))
      value = value.and.MA_pop_stack(G(2,1))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not. value) call errquit('error poping stack memory',0)

      do u=1,3
      do s=u+1,3
         Bus(s,u) = Bus(u,s)
      end do
      end do
      do v=1,3
      do u=1,3
         euv(u,v) = -elocal*hm(u,v)
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do


      call nwpw_timing_end(5)
      return 
      end


*     ***********************************
*     *									*
*     *	 	   v_nonlocal_euv_2			*
*     *									*
*     ***********************************

      subroutine v_nonlocal_euv_2(ispin,ne,psi1,euv)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 euv(3,3)

#include "mafdecls.fh"

*     **** common blocks ****
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


*     **** common blocks ****
      integer    dvl(2)
      integer    dvnl(2)
      integer    npsp2
      common / psp2_block / dvl,dvnl,npsp2


*     *** local variables ***
      integer nfft3d,npack1
      integer i,ii,ia,l,nn
      integer s,u,v
      real*8  omega,Bus(3,3),hm(3,3),c_l
      real*8  pi,scal
      integer exi(2),vtmp(2),tmp1(2),sw1(2),sw2(2)
      integer G(2,3)
      logical value

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega,lattice_unitg
      external is_sORd
      external ion_nion,ion_katm,G_indx
      external lattice_omega,lattice_unitg

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, npack1,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',G(2,1),G(1,1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',G(2,2),G(1,2))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',G(2,3),G(1,3))
      value = value.and.
     >        MA_push_get(mt_dbl,nn,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,9*nn,'sw2',sw2(2),sw2(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(G(1,1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(G(1,2)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(G(1,3)))
      call Pack_t_pack(1,dbl_mb(G(1,1)))
      call Pack_t_pack(1,dbl_mb(G(1,2)))
      call Pack_t_pack(1,dbl_mb(G(1,3)))

      omega = lattice_omega()

*     ***********************
*     **** calculate Bus ****
*     ***********************
      call dcopy(9,0.0d0,0,Bus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(1,dcpl_mb(exi(1)))

        do l=1,int_mb(lmmax(1)+ia-1)

*          ****************************
*          **** calculate F^(lm)_I ****
*          ****************************

*          *** current function is s or d ****
           if (is_sORd(l,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
              call Pack_tc_Mul(1,dbl_mb(vnl(1)+(l-1)*npack1
     >                           +(ia-1)*npack1*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Pack_itc_Mul(1,dbl_mb(vnl(1)+(l-1)*npack1
     >                           +(ia-1)*npack1*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
           end if
           call Pack_cc_ndot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)))

*          **********************************
*          **** calculate dF^(lm)_I/dhus ****
*          **********************************
           do s=1,3
           do u=1,3
              call Pack_tt_Mul(1,dbl_mb(dvnl(1)
     >                                +(s-1)*npack1
     >                                +(l-1)*3*npack1
     >                                +(ia-1)*3*npack1*lmmax_max),
     >                         dbl_mb(G(1,u)),
     >                         dbl_mb(tmp1(1)))

*             *** current function is s or d ****
              if (is_sORd(l,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
                call Pack_tc_Mul(1,dbl_mb(tmp1(1)),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))

*             *** current function is p or f ****
              else
                 call Pack_itc_Mul(1,dbl_mb(tmp1(1)),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
              end if

              call Pack_cc_ndot(1,nn,
     >                        psi1,
     >                        dcpl_mb(vtmp(1)),
     >                        dbl_mb(sw2(1)+(u-1)*nn
     >                                     +(s-1)*nn*3))
          end do
          end do

          c_l =(1.0d0/dbl_mb(vnlnrm(1)+(l-1)+(ia-1)*lmmax_max))
          if (ispin.eq.1) c_l = 2.0d0*c_l
          do i=1,nn
            do s=1,3
            do u=1,3

             Bus(u,s) = Bus(u,s) 
     >         - 2.0d0/(omega)
     >          * dbl_mb(sw1(1)+(i-1))
     >          * c_l
     >          * dbl_mb(sw2(1)+(i-1)
     >                         +(u-1)*nn
     >                         +(s-1)*nn*3)
            end do
            end do
          end do

        end do
      end do

      value = MA_pop_stack(sw2(2))
      value = MA_pop_stack(sw1(2))
      value = MA_pop_stack(G(2,3))
      value = MA_pop_stack(G(2,2))
      value = MA_pop_stack(G(2,1))
      value = MA_pop_stack(tmp1(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))



*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     *** calculate euv = Sum(s) hm(s,v)*Bus(u,s)
      call dcopy(9,0.0d0,0,euv,1)
      do u=1,3
      do v=1,3
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do


      call nwpw_timing_end(6)
      return 
      end







*     ***********************************
*     *									*
*     *	 	   psp_stress_read			*
*     *									*
*     ***********************************

      subroutine psp_stress_read(fname,
     >                       version,
     >                       nfft,unita,
     >                       npack0,dvl,
     >                       npack1,lmmax,dvnl,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*20 fname
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      integer npack0
      real*8 dvl(*)
      integer npack1,lmmax
      real*8 dvnl(npack1,3,lmmax)
      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer i,n,l
      integer msglen
      character*255 full_filename


      call Parallel_taskid(taskid)
  
*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
      end if


*     **** send header data to all processors ****
#ifdef MPI
      msglen = 1
      call MPI_BCAST(version,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 3
      call MPI_BCAST(nfft,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
      msglen = 9
      call MPI_BCAST(unita,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,version,mitob(msglen),MASTER)
      msglen = 3
      call BRDCST(9+MSGINT,nfft,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)
#endif

 
*     *** read in vl 3d block ***
      call D3dB_t_read(1,5,tmp2,tmp)
      call Pack_t_pack(0,tmp2)
      call Pack_t_Copy(0,tmp2,dvl)

*     **** read in vnl 3d blocks ****
      do n=1,lmmax
      do i=1,3
         call D3dB_t_read(1,5,tmp2,tmp)
         call Pack_t_pack(1,tmp2)
         call Pack_t_Copy(1,tmp2,dvnl(1,i,n))
      end do
      end do

*     *** close fname binary file ***
      if (taskid.eq.MASTER) then
c       close(11)
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *									*
*     *	 	  psp_stress_readall  		*
*     *									*
*     ***********************************

      subroutine psp_stress_readall()
      implicit none
      
#include "mafdecls.fh"



*     **** psp common block ****
      integer    dvl(2)
      integer    dvnl(2)
      integer    npsp
   
      integer lmmax_max,lmax_max
      parameter (lmax_max=4,lmmax_max=16)
      common / psp2_block / dvl,dvnl,npsp


*     **** local variables ****
      integer ngp(3)
      real*8  unita(3,3)
      integer version,nfft3d,npack1,npack0
      integer ia,l,lmmax
      character*12 boundry
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box
      character*4  element
      character*20 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind
      integer      control_ngrid
      integer      psp_lmmax
      real*8       control_unita
      character*12 control_boundry
      character*2  ion_atom
      external     nwpw_filefind
      external     control_ngrid
      external     psp_lmmax
      external     control_unita
      external     control_boundry
      external     ion_atom
      

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call Parallel_taskid(taskid)

      value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0)

      value = MA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** read pseudopotentials ****
      do ia=1,npsp

*       **** define formatted psp name ****
         element = '    '
         element = ion_atom(ia)
         l = index(element,' ') - 1
         fname = element(1:l)//'.vpp2'
         lmmax = psp_lmmax(ia)

*        **** not finished ****
         found = .false.
         do while (.not. found)
         if (nwpw_filefind(fname)) then
            call psp_stress_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  npack0,
     >                  dbl_mb(dvl(1) + (ia-1)*npack0),
     >                  npack1,lmmax,
     >                  dbl_mb(dvnl(1)+ (ia-1)*3*npack1*lmmax_max),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)


           if (ierr.gt.0) go to 9000
           boundry = control_boundry()
           l =index(boundry,' ') - 1
           if ((boundry(1:l).ne.'periodic') .or. (version.ne.3)) then
              if (taskid.eq.MASTER) then
                 write(6,*) 
     >             "Error reading formatted pseudoptental stress file:"
                 write(6,*) "Formatted pseudopotential was not",
     >                   " generated with periodic boundary condtions"
              end if
           end if

           correct_box = .true.
           if ( (ngp(1).ne.control_ngrid(1)) .or.
     >       (ngp(2).ne.control_ngrid(2)) .or.
     >       (ngp(3).ne.control_ngrid(3)) .or.
     >       (unita(1,1).ne.control_unita(1,1)) .or.
     >       (unita(2,1).ne.control_unita(2,1)) .or.
     >       (unita(3,1).ne.control_unita(3,1)) .or.
     >       (unita(1,2).ne.control_unita(1,2)) .or.
     >       (unita(2,2).ne.control_unita(2,2)) .or.
     >       (unita(3,2).ne.control_unita(3,2)) .or.
     >       (unita(1,3).ne.control_unita(1,3)) .or.
     >       (unita(2,3).ne.control_unita(2,3)) .or.
     >       (unita(3,3).ne.control_unita(3,3)) .or.
     >       ((boundry(1:l).eq.'periodic').and.(version.ne.3)).or.
     >       ((boundry(1:l).eq.'aperiodic').and.(version.ne.4))) then
              correct_box = .false.
              if (taskid.eq.MASTER) then
              write(6,*) "pseudopotential is not correctly formatted:",
     >                    fname
              end if
           end if
           if (correct_box) found = .true.
         end if

*        **** generate formatted pseudopotential atom.vpp2 *****
         if (.not.found) then
             call psp_stress_formatter_auto(ion_atom(ia))
         end if

         end do !*** do while ****
 
      end do
 9000 value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp(2))

      return
      end
