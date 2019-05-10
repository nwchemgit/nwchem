c
c $Id$
c
  
***********************************************************
*                                                         *
*                   PSPW-Charge module                    *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*    -- developed by Eric J. Bylaska on October 18,2001   *
*                                                         *
***********************************************************
*
*
*  

*     **********************************
*     *	                               *
*     *          pspw_charge_init      *
*     *                                *
*     **********************************
 
      subroutine pspw_charge_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"
#include "pspw_charge.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer i,pzero,zero,G(3),npack0,nfft3d,gk(2),gk2(2),taskid
      real*8  gg,temp,temp2,pi,fourpi,epsilon,scal2
      character*16 t

*     **** external functions ****
      logical     ion_chargeexist,control_qmmm
      integer     G_indx
      real*8      lattice_omega
      external    ion_chargeexist,control_qmmm
      external    G_indx
      external    lattice_omega


*     *************************
*     **** read in charge data ****
*     *************************
      bqext  = geom_extbq_on()
      bqgeom = ion_chargeexist()

     
      if(bqext.and.bqgeom)
     +  call errquit("pspw_charge_init: cannot have both bq structures",
     +               0,0)
      if(bqext) then
        nion = geom_extbq_ncenter()
        charge(1) = geom_extbq_charge()
        r1(1)     = geom_extbq_coord()
      else if (bqgeom) then
         value = geom_create(geom,'chargegeometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'chargegeometry')
         value = value.and.geom_ncent(geom,nion)
         if (.not. value) 
     >     call errquit('pspw_charge_init: opening chargegeometry',0,
     &       GEOM_ERR)

*        ***** Allocate pspw_charge geometries, and charges *****
         value = BA_alloc_get(mt_dbl,(3*nion),
     >                        'charge_r1',r1(2),r1(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,(nion),
     >                        'charge_charge',charge(2),charge(1))
         if (.not. value)
     >      call errquit('pspw_charge_init: pushing heap memory',0,
     &       MA_ERR)

*        **** read in charge geometries, and charges ****
         do i=1,nion
             value = value.and.
     >       geom_cent_get(geom,i,t,
     >                     dbl_mb(r1(1)+(i-1)*3),
     >                     dbl_mb(charge(1)+i-1))
         end do
         if (.not. value) 
     >   call errquit('pspw_charge_init:error reading chargegeometry',0,
     &       GEOM_ERR)

      else 
         nion = 0
      end if


*     **** don't do anything if there are no charges atoms ****
      if (nion.gt.0) then

*       **** set charge_type ****
        if (.not.btdb_get(rtdb,
     >         'nwpw:charge_type',mt_int,1,charge_type)) then
          charge_type = 1
          if (control_qmmm()) charge_type = 3 ! set charge type 3 for qmmm
        end if
 
        if (charge_type.eq.2) then
           if (.not.btdb_get(rtdb,
     >         'nwpw:charge_dielectric_d',mt_dbl,1,s_d))
     >       s_d = 0.5d0/0.529177d0

           if (.not.btdb_get(rtdb,
     >         'nwpw:charge_dielectric_rho',mt_dbl,1,s_rho))
     >       s_rho = 2.5d0/0.529177d0

           if (.not.btdb_get(rtdb,
     >         'nwpw:charge_rcut',mt_dbl,1,s_sigma))
     >       s_sigma = 0.8d0/0.529177d0

        end if
        if (charge_type.eq.3) then
           if (.not.btdb_get(rtdb,'nwpw:charge_rcut',mt_dbl,1,s_sigma))
     >       s_sigma = (0.37d0/0.3d0)/0.529177d0  ! should give covalent radii for spc water
           if (.not.btdb_get(rtdb,'nwpw:charge_ncut',mt_int,1,n_sigma))
     >       n_sigma = 4
        end if

        if (charge_type.eq.4) then
           if (.not.btdb_get(rtdb,
     >         'nwpw:charge_rcut',mt_dbl,1,s_sigma))
     >       s_sigma = 0.8d0
           epsilon = 1.0d0/s_sigma


           call Parallel_taskid(taskid)
           pi = 4.0d0*datan(1.0d0)
           fourpi = 4.0d0*pi
           scal2  = 1.0d0/lattice_omega()
           call D3dB_nfft3d(1,nfft3d)
           call Pack_npack(0,npack0)
           G(1) = G_indx(1)
           G(2) = G_indx(2)
           G(3) = G_indx(3)

*          **** allocate vc_charge memory ****
           value =           BA_alloc_get(mt_dbl,npack0,'point_vc',
     >                                    vc_charge(2),vc_charge(1))
           value = value.and.BA_alloc_get(mt_dbl,npack0,'point_vc2',
     >                                    vc2_charge(2),vc2_charge(1))
           if (.not.value) call errquit('out of heap memory',0,MA_ERR)
           value = BA_push_get(mt_dbl,nfft3d,'gk',gk(2),gk(1))
           value = value.and.
     >             BA_push_get(mt_dbl,nfft3d,'gk2',gk2(2),gk2(1))
           if (.not.value) call errquit('out of stack memory',0,MA_ERR)


*          ***** find the G==0 point in the lattice *****
           call D3dB_ijktoindexp(1,1,1,1,zero,pzero)

*          **** long-range and short-range part of Greens function ****
           call dcopy(nfft3d,0.0d0,0,dbl_mb(gk(1)), 1)
           call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
           do i=1,nfft3d

              gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >              + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >              + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )
               if ((pzero.eq.taskid) .and. (i.eq.zero)) then
                 !temp = pi/epsilon**2
                 temp = 0.0d0
              else
               temp = (fourpi/gg)*dexp(-gg/(4.0d0*epsilon**2))
               temp = (fourpi/gg)
              end if
              temp2 = 2.0d0*dexp(-gg/(epsilon**2))

              dbl_mb(gk(1) +i-1) = temp
              dbl_mb(gk2(1)+i-1) = temp2
           end do
           call Pack_t_pack(0,dbl_mb(gk(1)))
           call Pack_t_Copy(0,dbl_mb(gk(1)),dbl_mb(vc_charge(1)))
           call Pack_t_pack(0,dbl_mb(gk2(1)))
           call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(vc2_charge(1)))
c           call Pack_t_SMul(0,scal2,dbl_mb(vc_charge(1)),
c     >                              dbl_mb(vc_charge(1)))
c           call Pack_t_SMul(0,scal2,dbl_mb(vc2_charge(1)),
c     >                              dbl_mb(vc2_charge(1)))
           call Pack_t_SMul1(0,scal2,dbl_mb(vc_charge(1)))
           call Pack_t_SMul1(0,scal2,dbl_mb(vc2_charge(1)))

           value =           BA_pop_stack(gk2(2))
           value = value.and.BA_pop_stack(gk(2))
           if (.not. value)
     >     call errquit('pspw_charge_init_init: poppingstack memory',0,
     >       MA_ERR)

           call qstrfac_init(nion)
           call qphafac_rion(dbl_mb(r1(1)))
        end if
      

      end if

      return
      end



*     **********************************
*     *	                               *
*     *        pspw_charge_destroy     *
*     *                                *
*     **********************************
      subroutine pspw_charge_destroy()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"
#include "pspw_charge.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

      if (bqgeom) then
        value = geom_destroy(geom)

        if (.not. value) 
     >     call errquit('error destroying chargegeometry', 0, GEOM_ERR)
      end if

      return
      end

*     **********************************
*     *	                               *
*     *        pspw_charge_write       *
*     *                                *
*     **********************************

      subroutine pspw_charge_write(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"
#include "errquit.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      logical value
      integer i
      double precision rxyz(3),q
      character*16 t

      integer  control_code
      external control_code

*     *******************************
*     **** write out charge data ****
*     *******************************
      if (bqgeom) then
      if (nion.gt.0) then
        value = .true.
        do i=1,nion
           value = value.and.geom_cent_get(geom,i,t,rxyz,q)
           value = value.and.
     >             geom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
        end do

        value = value.and.geom_rtdb_delete(rtdb,'chargegeometry')
        if ((control_code().eq.1).or.
     >      (control_code().eq.2)) then
        value = value.and.geom_rtdb_store(rtdb,geom,'chargegeometry')
        end if
        value = value.and.geom_destroy(geom)
        if (.not. value) call errquit('error writing chargegeometry', 0,
     &       GEOM_ERR)
      end if
      end if

      return
      end




*     **********************************
*     *	                               *
*     *         pspw_charge_aname        *
*     *                                *
*     **********************************
      character*2 function pspw_charge_aname(i)
      implicit none
#include "errquit.fh"
      integer i

#include "stdio.fh"
#include "geom.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      character*2  symbol
      character*16 t,name
      real*8 q

      if (bqgeom) then
        if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' pspw_charge_aname: failed ',i, GEOM_ERR)

        if (.not. geom_tag_to_element(t,symbol,name,q)) then
          symbol = 'bq'
        end if
      else
         symbol = 'bq'
      end if

      pspw_charge_aname = symbol
      return
      end


*     ***************************
*     *                         *
*     *   pspw_charge_charge    *
*     *                         *
*     ***************************
      real*8 function pspw_charge_charge(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "pspw_charge.fh"


      pspw_charge_charge = dbl_mb(charge(1)+i-1)
      return
      end


*     **********************************
*     *	                               *
*     *         pspw_charge_found        *
*     *                                *
*     **********************************
      logical function pspw_charge_found()
      implicit none

#include "pspw_charge.fh"
 
      logical value

      value = .false.
      if (nion.gt.0) value = .true.

      pspw_charge_found = value
      return
      end
c
      logical function pspw_bqext()
      implicit none

#include "pspw_charge.fh"
 
      logical value

      pspw_bqext = bqext

      return
      end

*     **********************************
*     *	                               *
*     *         pspw_charge_end	       *
*     *                                *
*     **********************************

      subroutine pspw_charge_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_charge.fh"
#include "errquit.fh"

      logical value

      if (nion.gt.0) then
        value = .true.
        if (bqgeom) then
           value = value.and.BA_free_heap(r1(2))
           value = value.and.BA_free_heap(charge(2))
        end if

        if (charge_type.eq.4) then
          value=value.and.BA_free_heap(vc2_charge(2))
          value=value.and.BA_free_heap(vc_charge(2))
          call qstrfac_end()
        end if
        if (.not.value) 
     >    call errquit('pspw_charge_end: freeing heap',0, MA_ERR)
      end if
      return
      end

*     **********************************
*     *	                               *
*     *        pspw_charge_nion        *
*     *                                *
*     **********************************
      integer function pspw_charge_nion()
      implicit none

#include "pspw_charge.fh"

      pspw_charge_nion = nion
      return
      end


*     **********************************
*     *	                               *
*     *       pspw_charge_Print        *
*     *                                *
*     **********************************

      subroutine pspw_charge_Print(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_charge_aname
      real*8      pspw_charge_charge
      external    pspw_charge_aname
      external    pspw_charge_charge

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
        if (charge_type.eq.1) then
          write(unit,1160)   
        end if
        if (charge_type.eq.2) then
          write(unit,1170)   
          write(unit,1171) s_d,s_rho,s_sigma   
        end if
        if (charge_type.eq.3) then
          write(unit,1175)
          write(unit,1176) n_sigma, s_sigma
        end if
        if (charge_type.eq.4) then
          write(unit,1177)
          write(unit,1178) s_sigma
        end if
        write(unit,1180)   
        do ii=1,nion
          write(unit,1190) ii,pspw_charge_aname(ii),
     >                     (dbl_mb(r1(1)+3*(ii-1)+k-1),k=1,3),
     >                       pspw_charge_charge(ii)

        end do
      end if

      end if
      return

 1160 format(/' Point charges used:')
 1170 format(/' Switching charges used:')
 1171 FORMAT(5X, '(d,rho,sigma) = (',3E11.3,' )  ')
 1175 format(/' qmmm charges used:')
 1176 format(5x,'(n_sigma,s_sigma) = (',I3,E11.3,')  ')
 1177 format(/' periodic charges used:')
 1178 format(5x,' s_sigma=',E11.3)
 1180 FORMAT(/' position of point charges (au): ')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - charge = ',F10.3)
      end




*     **********************************
*     *	                               *
*     *       pspw_charge_PrintXYZ       *
*     *                                *
*     **********************************

      subroutine pspw_charge_PrintXYZ(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_charge_aname
      external    pspw_charge_aname

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then  
       do ii=1,nion
        write(unit,*) pspw_charge_aname(ii),'      ',
     >               (dbl_mb(r1(1)+3*(ii-1)+k-1)*0.529177d0,k=1,3)
       end do
      end if

      end if

      return
      end


*     **********************************
*     *	                               *
*     *        pspw_charge_rion        *
*     *                                *
*     **********************************
      real*8 function pspw_charge_rion(i,ii)
      implicit none
      integer i,ii

#include "bafdecls.fh"
#include "pspw_charge.fh"

      pspw_charge_rion = dbl_mb(r1(1)+3*(ii-1)+i-1)
      return
      end






*     **********************************
*     *	                               *
*     *     pspw_charge_Generate_V     *
*     *                                *
*     **********************************

      subroutine pspw_charge_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)

#include "bafdecls.fh"
#include "pspw_charge.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer ii,ia,k,p1,p2,po2,npack0,exi(2),tmp(2),nx,ny,nz
      real*8  x1,y1,z1,q1,r,r2,rc,epsilon,ttt,sss,scal1,scal2
      real*8  rc1,rc2,rr1,rr1p

*     **** external functions ****
      integer  ion_katm
      external ion_katm
      real*8   util_erf,lattice_omega
      external util_erf,lattice_omega


      if (nion.gt.0) then

*        **** point charges ****
         if (charge_type.eq.1) then
           do ii=1,nion
              x1 = dbl_mb(r1(1)+3*(ii-1))
              y1 = dbl_mb(r1(1)+3*(ii-1)+1)
              z1 = dbl_mb(r1(1)+3*(ii-1)+2)
              q1 = dbl_mb(charge(1)+ii-1)
!$OMP DO
              do k=1,n2ft3d
                r = (rgrid(1,k)-x1)**2
     >            + (rgrid(2,k)-y1)**2
     >            + (rgrid(3,k)-z1)**2
                r = dsqrt(r)
                Vqm(k) = Vqm(k) - q1/r
              end do
!$OMP END DO
           end do   

*        **** switching charge potential **** 
         else if (charge_type.eq.2) then
           do ii=1,nion
              x1 = dbl_mb(r1(1)+3*(ii-1))
              y1 = dbl_mb(r1(1)+3*(ii-1)+1)
              z1 = dbl_mb(r1(1)+3*(ii-1)+2)
              q1 = dbl_mb(charge(1)+ii-1)
!$OMP DO
              do k=1,n2ft3d
                r = (rgrid(1,k)-x1)**2
     >            + (rgrid(2,k)-y1)**2
     >            + (rgrid(3,k)-z1)**2
                r = dsqrt(r)

*               **** define dielectric switching function ****
                if (r.le.s_d) then
                   epsilon = 0.0d0
                else if (r.lt.(s_d+s_rho)) then
                   epsilon = 1.0d0-(1.0d0-(r-s_d)**2/s_rho**2)**2        
                else
                   epsilon = 1.0d0
                end if

                Vqm(k) = Vqm(k) - epsilon*q1*util_erf(r/s_sigma)/r
              end do
!$OMP END DO
           end do   

*        **** short-range repulsion ****
         else if (charge_type.eq.3) then
           p1  = n_sigma
           po2 = n_sigma/2
           p2  = n_sigma+1
           do ii=1,nion
              x1 = dbl_mb(r1(1)+3*(ii-1))
              y1 = dbl_mb(r1(1)+3*(ii-1)+1)
              z1 = dbl_mb(r1(1)+3*(ii-1)+2)
              q1 = dbl_mb(charge(1)+ii-1)
              rc = abs(q1)*s_sigma
              if (rc.gt.1.0d0) rc=1.0d0 !stupid check
              !rc1 = abs(q1)*rc**p1
              rc1 = rc**p1
              rc2 = rc**p2
!$OMP DO
              do k=1,n2ft3d
                 r2 = (rgrid(1,k)-x1)**2
     >              + (rgrid(2,k)-y1)**2
     >              + (rgrid(3,k)-z1)**2
                 r = dsqrt(r2)
                 rr1  = r**p1
                 rr1p = r**p2
                 ttt = (rc1 - rr1)
                 sss = (rc2 - rr1p)

                 Vqm(k) = Vqm(k) - q1*(ttt/sss)
              end do
!$OMP END DO
           end do   
           !call sub_coulomb_screened_init_print(Vqm,1.0d0,q1)

         else if (charge_type.eq.4) then

           call D3dB_nx(1,nx)
           call D3dB_ny(1,ny)
           call D3dB_nz(1,nz)
           scal1  = 1.0d0/dble(nx*ny*nz)
           scal2  = 1.0d0/lattice_omega()
           call Pack_npack(0,npack0)
           value = BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
           value = value.and.
     >             BA_push_get(mt_dcpl,npack0,'tmp',tmp(2),tmp(1))
           if (.not.value) call errquit('out of stack memory',0,MA_ERR)

           !call dcopy(n2ft3d,0.0d0,0,Vqm,1)
           !call dcopy(2*npack0,0.0d0,0,dcpl_mb(tmp(1)),1)
           !call dcopy(2*npack0,0.0d0,0,dcpl_mb(exi(1)),1)
           call Parallel_shared_vector_zero(.false.,n2ft3d,Vqm)
           call Parallel_shared_vector_zero(.false.,2*npack0,
     >                                      dcpl_mb(tmp(1)))
           call Parallel_shared_vector_zero(.true.,2*npack0,
     >                                      dcpl_mb(exi(1)))

           do ii=1,nion
             q1 = -dbl_mb(charge(1)+ii-1)
             call qstrfac_pack(0,ii,dcpl_mb(exi(1)))
             call Pack_cc_daxpy(0,q1,dcpl_mb(exi(1)),Vqm)
             call Pack_cc_daxpy(0,dabs(q1),dcpl_mb(exi(1)),
     >                                     dcpl_mb(tmp(1)))
           end do
c           call Pack_tc_Mul(0,dbl_mb(vc_charge(1)),Vqm,Vqm)
c           call Pack_tc_Mul(0,dbl_mb(vc2_charge(1)),dcpl_mb(tmp(1)),
c     >                                              dcpl_mb(tmp(1)))
           call Pack_tc_Mul2(0,dbl_mb(vc_charge(1)),Vqm)
           call Pack_tc_Mul2(0,dbl_mb(vc2_charge(1)),dcpl_mb(tmp(1)))
           !call Pack_cc_Sum(0,dcpl_mb(tmp(1)),Vqm,Vqm)
           call Pack_c_unpack(0,Vqm)
           call D3dB_cr_fft3b(1,Vqm)
           value =           BA_pop_stack(tmp(2))
           value = value.and.BA_pop_stack(exi(2))
           if (.not. value)
     >     call errquit('pspw_charge_Generate_V:poppingstack memory',0,
     >       MA_ERR)

           !call sub_coulomb_screened_init_print(Vqm,1.0d0,q1)


         end if


      end if

      return
      end


*     **********************************
*     *	                               *
*     *     pspw_charge_rho_Fcharge    *
*     *                                *
*     **********************************

      subroutine pspw_charge_rho_Fcharge(n2ft3d,rgrid,rho,dv,fcharge)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fcharge(3,*)

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer i,j,p1,po2,p2
      real*8 rx,ry,rz,c
      real*8 fx,fy,fz,verf,yerf,v
      real*8 x,y,z,q,r,r2,rc,epsilon,depsilon,sqrt_pi
      real*8 ttt,sss,rr1,rr1m,rr1p,rc1,rc2,rc1p
      double precision tttp,sssp

*     **** external functions ****
      real*8   util_erf
      external util_erf

      if (nion.gt.0) then
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

*     **** point charges ****
      if (charge_type.eq.1) then
      do j=1,nion
         x = dbl_mb(r1(1)+3*(j-1))
         y = dbl_mb(r1(1)+3*(j-1)+1)
         z = dbl_mb(r1(1)+3*(j-1)+2)
         q = -dbl_mb(charge(1)+j-1)
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,n2ft3d
            rx = x - rgrid(1,i)
            ry = y - rgrid(2,i)
            rz = z - rgrid(3,i)
            r  = dsqrt( rx**2 + ry**2 + rz**2)

            if (r .gt. 1.0d-8) then
              v    = -q/r**3
            else
              v = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fcharge(1,j) = fcharge(1,j) - fx*dv
         fcharge(2,j) = fcharge(2,j) - fy*dv
         fcharge(3,j) = fcharge(3,j) - fz*dv
      end do


*     **** switching charge potential **** 
      else if (charge_type.eq.2) then
      do j=1,nion
         x = dbl_mb(r1(1)+3*(j-1))
         y = dbl_mb(r1(1)+3*(j-1)+1)
         z = dbl_mb(r1(1)+3*(j-1)+2)
         q = -dbl_mb(charge(1)+j-1)
         c = 1.0d0/s_sigma
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,n2ft3d
            rx = x - rgrid(1,i)
            ry = y - rgrid(2,i)
            rz = z - rgrid(3,i)
            r  = dsqrt( rx**2 + ry**2 + rz**2)

            if (r .gt. 1.0d-8) then

*             **** define dielectric switching function ****
              if (r.le.s_d) then
                 epsilon  = 0.0d0
                 depsilon = 0.0d0
              else if (r.lt.(s_d+s_rho)) then
                 epsilon = 1.0d0-(1.0d0-(r-s_d)**2/s_rho**2)**2
                 depsilon = 4.0d0*((r-s_d)/s_rho**2)
     >                           *(1.0d0-(r-s_d)**2/s_rho**2)        
              else
                 epsilon  = 1.0d0
                 depsilon = 0.0d0
              end if
              yerf=r*c
              verf = util_erf(yerf)
              v    = epsilon*q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                          - verf)/r**3
     >             + depsilon*q*verf/r**2
            else
              v = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fcharge(1,j) = fcharge(1,j) - fx*dv
         fcharge(2,j) = fcharge(2,j) - fy*dv
         fcharge(3,j) = fcharge(3,j) - fz*dv
      end do

*     **** short-range repulsion ****
      else if (charge_type.eq.3) then
      p1  = n_sigma
      p2  = n_sigma+1
      po2 = n_sigma/2
      do j=1,nion
         x = dbl_mb(r1(1)+3*(j-1))
         y = dbl_mb(r1(1)+3*(j-1)+1)
         z = dbl_mb(r1(1)+3*(j-1)+2)

         q = -dbl_mb(charge(1)+j-1)
         rc = abs(q)*s_sigma
         if (rc.gt.1.0d0) rc=1.0d0 !stupid check
         rc1 = rc**p1
         rc2 = rc**p2

         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,n2ft3d
            rx = x - rgrid(1,i)
            ry = y - rgrid(2,i)
            rz = z - rgrid(3,i)
            r  = dsqrt(rx*rx + ry*ry + rz*rz)

            !rr1m = r**(p1-1)
            !rr1  = r*rr1m
            !rr1p = r2*rr1m
            rr1m = r**(p1-1)
            rr1  = r**p1
            rr1p = r**p2

            ttt = (rc1 -  rr1)
            sss = (rc2 - rr1p)

            tttp = -p1*rr1m
            sssp = -p2*rr1

            if (r .gt. 1.0d-8) then
              v  = q*(tttp - (ttt/sss)*sssp)/(sss*r)
            else
              v  = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fcharge(1,j) = fcharge(1,j) - fx*dv
         fcharge(2,j) = fcharge(2,j) - fy*dv
         fcharge(3,j) = fcharge(3,j) - fz*dv
      end do

      end if


      end if

      return
      end






*     **********************************
*     *	                               *
*     *     pspw_charge_Energy_ion     *
*     *                                *
*     **********************************

      real*8 function pspw_charge_Energy_ion()
      implicit none

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,ii,ia
      real*8  qi, xi, yi, zi
      real*8  qii,xii,yii,zii
      real*8  r,energy

*     **** external functions ****
      integer   ion_nion,ion_katm
      real*8    ion_rion,psp_zv
      external  ion_nion,ion_katm
      external  ion_rion,psp_zv

      energy = 0.0d0
      do ii=1,ion_nion()
         ia = ion_katm(ii)
         qii = psp_zv(ia)
         xii = ion_rion(1,ii)
         yii = ion_rion(2,ii)
         zii = ion_rion(3,ii)
         do i=1,nion
            xi = dbl_mb(r1(1)+3*(i-1))
            yi = dbl_mb(r1(1)+3*(i-1)+1)
            zi = dbl_mb(r1(1)+3*(i-1)+2)
            qi = dbl_mb(charge(1)+i-1)
            r = (xi-xii)**2 + (yi-yii)**2 + (zi-zii)**2
            r = dsqrt(r)
            energy = energy + qi*qii/r
         end do
      end do


      pspw_charge_Energy_ion = energy
      return
      end




*     **********************************
*     *	                               *
*     *     pspw_charge_Fion_Fcharge   *
*     *                                *
*     **********************************

*    This routine calculates the forces between the QM ions and point charges.
*
*     Exit - fion: force on the QM ions
*            fcharge: force on the point charges
*
*     Uses - ion_nion,Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_Fion_Fcharge(fion,fcharge)
      implicit none
      real*8 fion(3,*)
      real*8 fcharge(3,*)

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,j,nion_ion
      real*8  qi,ri(3),qj,rj(3)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion

      if (nion.gt.0) then      
     
        nion_ion = ion_nion()

        do i=1,nion_ion
           qi = psp_zv(ion_katm(i))
           ri(1) = ion_rion(1,i)
           ri(2) = ion_rion(2,i)
           ri(3) = ion_rion(3,i)
           do j=1,nion
              qj = dbl_mb(charge(1)+j-1)            
              rj(1) = dbl_mb(r1(1)+3*(j-1))
              rj(2) = dbl_mb(r1(1)+3*(j-1)+1)
              rj(3) = dbl_mb(r1(1)+3*(j-1)+2)
              call Q_Electrostatic_Force(ri,qi,fion(1,i),
     >                                 rj,qj,fcharge(1,j))
           end do
        end do

      end if
      return
      end


*     **********************************
*     *	                               *
*     *    pspw_charge_charge_Fion     *
*     *                                *
*     **********************************

*    This routine calculates the forces between the QM ions and point charges.
*
*     Exit - fion: force on the QM ions
*            fcharge: force on the point charges
*
*     Uses - ion_nion,Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_charge_Fion(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,j,nion_ion
      real*8  qi,ri(3),qj,rj(3)
      real*8  fcharge(3)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion

      if (nion.gt.0) then      
     
        fcharge(1) = 0.0d0
        fcharge(2) = 0.0d0
        fcharge(3) = 0.0d0
        nion_ion = ion_nion()

        do i=1,nion_ion
           qi = psp_zv(ion_katm(i))
           ri(1) = ion_rion(1,i)
           ri(2) = ion_rion(2,i)
           ri(3) = ion_rion(3,i)
           do j=1,nion
              qj = dbl_mb(charge(1)+j-1)            
              rj(1) = dbl_mb(r1(1)+3*(j-1))
              rj(2) = dbl_mb(r1(1)+3*(j-1)+1)
              rj(3) = dbl_mb(r1(1)+3*(j-1)+2)
              call Q_Electrostatic_Force(ri,qi,fion(1,i),
     >                                 rj,qj,fcharge)
           end do
        end do

      end if
      return
      end


*     **********************************
*     *	                               *
*     *     pspw_charge_Energy_charge  *
*     *                                *
*     **********************************

      real*8 function pspw_charge_Energy_charge()
      implicit none

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,ii
      real*8  qi, xi, yi, zi
      real*8  qii,xii,yii,zii
      real*8  r,energy

*     **** external functions ****
      logical  nwpw_bqbq
      external nwpw_bqbq

      energy = 0.0d0
      if ((nion.gt.0).and.nwpw_bqbq()) then
        do ii=1,nion
           qii = dbl_mb(charge(1)+ii-1)            
           xii = dbl_mb(r1(1)+3*(ii-1))
           yii = dbl_mb(r1(1)+3*(ii-1)+1)
           zii = dbl_mb(r1(1)+3*(ii-1)+2)
           qii = dbl_mb(charge(1)+ii-1)
           do i=1,ii-1
              xi = dbl_mb(r1(1)+3*(i-1))
              yi = dbl_mb(r1(1)+3*(i-1)+1)
              zi = dbl_mb(r1(1)+3*(i-1)+2)
              qi = dbl_mb(charge(1)+i-1)
              r = (xi-xii)**2 + (yi-yii)**2 + (zi-zii)**2
              r = dsqrt(r)
              energy = energy + qi*qii/r
           end do
        end do
      end if

      pspw_charge_Energy_charge = energy
      return
      end


*     **********************************
*     *	                               *
*     *     pspw_charge_Fcharge        *
*     *                                *
*     **********************************

*    This routine calculates the force between the point charges.
*
*     Exit - fcharge: force on the point charges
*
*     Uses - Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_Fcharge(fcharge)
      implicit none
      real*8 fcharge(3,*)

#include "bafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,j
      real*8  qi,ri(3),qj,rj(3)

*     **** external functions ****
      logical  nwpw_bqbq
      external nwpw_bqbq

      if ((nion.gt.0).and.(nwpw_bqbq())) then      
     
        do j=1,nion
           qj = dbl_mb(charge(1)+j-1)            
           rj(1) = dbl_mb(r1(1)+3*(j-1))
           rj(2) = dbl_mb(r1(1)+3*(j-1)+1)
           rj(3) = dbl_mb(r1(1)+3*(j-1)+2)
           do i=1,j-1
              qi = dbl_mb(charge(1)+i-1)            
              ri(1) = dbl_mb(r1(1)+3*(i-1))
              ri(2) = dbl_mb(r1(1)+3*(i-1)+1)
              ri(3) = dbl_mb(r1(1)+3*(i-1)+2)
              call Q_Electrostatic_Force(ri,qi,fcharge(1,i),
     >                                   rj,qj,fcharge(1,j))
           end do
        end do

      end if
      return
      end



*********************************************************************************************************




*     **********************************
*     *	                               *
*     *  pspw_charge_Generate_V_qmmm   *
*     *                                *
*     **********************************

      subroutine pspw_charge_Generate_V_qmmm(r1,charge1,n1,
     >                                       n2ft3d,rgrid,Vqm)
      implicit none
      real*8 r1(3,*),charge1(*)
      integer n1
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)


*     ***** local variables ****
      integer ii,ia,k,p1,p2,po2
      real*8 x1,y1,z1,q1,r,r2,rc,ttt,sss
      real*8 rc1,rc2,rr1,rr1p
      real*8  s_sigma
      integer n_sigma

      if (n1.gt.0) then

        s_sigma = ((0.37d0/0.3d0)/0.529177d0)  ! should give covalent radii for spc water
        n_sigma = 4

*        **** short-range repulsion ****
         p1  = n_sigma
         po2 = n_sigma/2
         p2  = n_sigma+1
         do ii=1,n1
            x1 = r1(1,ii)
            y1 = r1(2,ii)
            z1 = r1(3,ii)
            q1 = charge1(ii)
c           a way to define covalent radius
            rc = abs(q1)*s_sigma
            if (rc.gt.1.0d0) rc=1.0d0 !stupid check
            !rc1 = abs(q1)*rc**p1
            rc1 = rc**p1
            rc2 = rc**p2
            do k=1,n2ft3d
               r2 = (rgrid(1,k)-x1)**2
     >            + (rgrid(2,k)-y1)**2
     >            + (rgrid(3,k)-z1)**2
               r = dsqrt(r2)
               rr1  = r**p1
               rr1p = r**p2
               ttt = (rc1 - rr1)
               sss = (rc2 - rr1p)

               Vqm(k) = Vqm(k) - q1*(ttt/sss)
            end do
         end do   


      end if

      return
      end






*     **********************************
*     *	                               *
*     *  pspw_charge_rho_Fcharge_qmmm  *
*     *                                *
*     **********************************
*
*     Calculates the force on the qmmm charges from the electronic density
*

      subroutine pspw_charge_rho_Fcharge_qmmm(r1,qcharge,nion,
     >                                     n2ft3d,rgrid,rho,dv,fcharge)
      implicit none
      real*8  r1(3,*),qcharge(*)
      integer nion
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fcharge(3,*)


*     ***** local variables ****
      integer i,j,p1,po2,p2
      real*8 rx,ry,rz,c
      real*8 fx,fy,fz,v
      real*8 x,y,z,q,r,r2,rc
      real*8 ttt,sss,rr1,rr1m,rr1p,rc1,rc2,rc1p
      double precision tttp,sssp
      integer n_sigma
      real*8  s_sigma

*     **** external functions ****

      if (nion.gt.0) then

      s_sigma = (0.37d0/0.3d0)/0.529177d0  ! should give covalent radii for spc water
      n_sigma = 4
 
      p1  = n_sigma
      p2  = n_sigma+1
      po2 = n_sigma/2
      do j=1,nion
         x = r1(1,j)
         y = r1(2,j)
         z = r1(3,j)
         q = -qcharge(j)
         rc = abs(q)*s_sigma
         if (rc.gt.1.0d0) rc=1.0d0  ! stupid check
         rc1 = rc**p1
         rc2 = rc**p2
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,n2ft3d
            rx = x - rgrid(1,i)
            ry = y - rgrid(2,i)
            rz = z - rgrid(3,i)
            r  = dsqrt(rx*rx + ry*ry + rz*rz)

            !rr1m = r**(p1-1)
            !rr1  = r*rr1m
            !rr1p = r2*rr1m
            rr1m = r**(p1-1)
            rr1  = r**p1
            rr1p = r**p2

            ttt = (rc1  - rr1)
            sss = (rc2 - rr1p)

            tttp = -p1*rr1m
            sssp = -p2*rr1

            if (r .gt. 1.0d-8) then
              v  = q*(tttp - (ttt/sss)*sssp)/(sss*r)
            else
              v  = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fcharge(1,j) = fcharge(1,j) - fx*dv
         fcharge(2,j) = fcharge(2,j) - fy*dv
         fcharge(3,j) = fcharge(3,j) - fz*dv

c         call D3dB_r_dsum(1,rho,sss)
c         write(*,*) "fq-elc a,",j,q,nion,n2ft3d,dv,sss*dv
c         write(*,*) "fq-elc b,",j,r1(1,j),r1(2,j),r1(3,j)
c         write(*,*) "fq-elc c,",j,fcharge(1,j),fcharge(2,j),fcharge(3,j)
c         write(*,*) "fq-elc d,",j,fx,fy,fz
c         write(*,*) "fq-elc e,",j,fx*dv,fy*dv,fz*dv
      end do


      end if

      return
      end





*     **********************************
*     *	                               *
*     *   pspw_charge_Energy_ion_qmmm  *
*     *                                *
*     **********************************

*    This routine calculates the energy between the QM ions and point charges.

      real*8 function pspw_charge_Energy_ion_qmmm(r1,qcharge,nion)
      implicit none
      real*8 r1(3,*),qcharge(*)
      integer nion


*     **** local variables ****
      integer i,ii,ia
      real*8  qi, xi, yi, zi
      real*8  qii,xii,yii,zii
      real*8  r,energy

*     **** external functions ****
      integer   ion_nion,ion_katm
      real*8    ion_rion,psp_zv
      external  ion_nion,ion_katm
      external  ion_rion,psp_zv

      energy = 0.0d0
      do ii=1,ion_nion()
         ia = ion_katm(ii)
         qii = psp_zv(ia)
         xii = ion_rion(1,ii)
         yii = ion_rion(2,ii)
         zii = ion_rion(3,ii)
         do i=1,nion
            xi = r1(1,i)
            yi = r1(2,i)
            zi = r1(3,i)
            qi = qcharge(i)
            r = (xi-xii)**2 + (yi-yii)**2 + (zi-zii)**2
            r = dsqrt(r)
            energy = energy + qi*qii/r
         end do
      end do

      pspw_charge_Energy_ion_qmmm = energy
      return
      end




*     **********************************
*     *	                               *
*     *  pspw_charge_Fion_Fcharge_qmmm *
*     *                                *
*     **********************************

*    This routine calculates the forces between the QM ions and point charges.
*
*     Exit - fion: force on the QM ions
*            fcharge: force on the point charges
*
*     Uses - ion_nion,Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_Fion_Fcharge_qmmm(r1,qcharge,nion,
     >                                         fion,fcharge)
      implicit none
      real*8 r1(3,*),qcharge(*)
      integer nion
      real*8 fion(3,*)
      real*8 fcharge(3,*)

*     **** local variables ****
      integer i,j,nion_ion
      real*8  qi,ri(3),qj,rj(3)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion

      if (nion.gt.0) then      
     
        nion_ion = ion_nion()

        do i=1,nion_ion
           qi = psp_zv(ion_katm(i))
           ri(1) = ion_rion(1,i)
           ri(2) = ion_rion(2,i)
           ri(3) = ion_rion(3,i)
           do j=1,nion
              qj = qcharge(j)            
              rj(1) = r1(1,j)
              rj(2) = r1(2,j)
              rj(3) = r1(3,j)
              call Q_Electrostatic_Force(ri,qi,fion(1,i),
     >                                 rj,qj,fcharge(1,j))
           end do
        end do

      end if
      return
      end

*     **********************************
*     *	                               *
*     *  pspw_charge_ion_Fcharge_qmmm  *
*     *                                *
*     **********************************

*    This routine calculates the forces between the QM ions and point charges.
*
*     Exit - fion: force on the QM ions
*            fcharge: force on the point charges
*
*     Uses - ion_nion,Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_ion_Fcharge_qmmm(r1,qcharge,nion,fcharge)
      implicit none
      real*8 r1(3,*),qcharge(*)
      integer nion
      real*8 fcharge(3,*)

*     **** local variables ****
      integer i,j,nion_ion
      real*8  qi,ri(3),qj,rj(3)
      real*8 fion(3)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion

      fion(1) = 0.0d0
      fion(2) = 0.0d0
      fion(3) = 0.0d0
      if (nion.gt.0) then      
     
        nion_ion = ion_nion()

        do i=1,nion_ion
           qi = psp_zv(ion_katm(i))
           ri(1) = ion_rion(1,i)
           ri(2) = ion_rion(2,i)
           ri(3) = ion_rion(3,i)
           do j=1,nion
              qj = qcharge(j)            
              rj(1) = r1(1,j)
              rj(2) = r1(2,j)
              rj(3) = r1(3,j)
              call Q_Electrostatic_Force(ri,qi,fion,
     >                                 rj,qj,fcharge(1,j))
           end do
        end do

      end if
      return
      end




