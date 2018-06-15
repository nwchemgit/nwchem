*     *****************************************
*     *                                       *
*     *             vdw_DF_init               *
*     *                                       *
*     *****************************************
*
*     *** Order of precedence for choosing name                     ***
*     *** 1) value of NWCHEM_QMESH_DATA environment variable        ***
*     *** 2) value of NWCHEM_QMESH_DATA set in $HOME/.nwchemrc file ***
*     *** 3) value of the compiled in library name                  ***
*
      subroutine vdw_DF_init()
      implicit none

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "bafdecls.fh"

*     **** local variables ****
      logical mprint,hprint,debug,does_it_exist
      logical from_environment,from_compile,from_nwchemrc
      integer iop,lgth,unitf,print_level,i,j
      character*255 qmesh_library,qmesh_data_name

*     **** qmesh common block ****
      logical has_vdw
      integer Nqs,qmesh(2),ya(2),y2a(2)
      common /vdw_df_common/ qmesh,ya,y2a,Nqs,has_vdw

*     **** external functions ***
      logical  control_has_vdw,constrol_is_vdw2
      external control_has_vdw,constrol_is_vdw2


*     **** return if not has_vdw ****
      has_vdw = control_has_vdw()
      is_vdw2 = control_is_vdw2()
      if (.not.has_vdw) return

      call util_print_get_level(print_level)
      mprint = print_medium.le.print_level
      hprint = print_high  .le.print_level
      debug  = print_debug .le.print_level
      from_environment = .false.
      from_nwchemrc    = .false.
      from_compile     = .false.


*     **** Try to get from NWCHEM_QMESH_DATA environment variable ****
      call util_getenv('NWCHEM_QMESH_DATA',qmesh_data_name)
      if (lgth.gt.0) then
         if (util_find_dir(qmesh_data_name)) then
            from_environment = .true.
            goto 99
         else
            write(luout,*)' warning:::::::::::::: from_environment'
            write(luout,*)' NWCHEM_QMESH_DATA set to: <',
     &       qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
            write(luout,*)' but file does not exist !'
            write(luout,*)' using compiled library'
         end if
      end if


*     **** Try to get from NWCHEM_QMESH_DATA nwchemrc ****
*2:   check for NWCHEM_QMESH_DATA defined in users .nwchemrc file
*     assumed structure in .nwchemrc file is variable [whitespace] value
*     one setting per line
*
      qmesh_library='nwchem_qmesh_data'
      call inp_save_state() ! save state of any inp unit
      if(.not.util_nwchemrc_get(qmesh_library,qmesh_data_name)) then
        if (debug) then
          write(luout,*)'util_nwchemrc_get failed'
        endif
      else
        does_it_exist=util_find_dir(qmesh_data_name)
        if (does_it_exist)then
          from_nwchemrc = .true.
          call inp_restore_state() ! restore state of any inp unit
          goto 99
        else
          write(luout,*)' warning:::::::::::::: from_nwchemrc'
          write(luout,*)' NWCHEM_QMESH_DATA set to: <',
     &     qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using compiled in library'
        endif
      endif
      call inp_restore_state() ! restore state of any inp unit



*     **** Try to get from compile ****
      from_compile = .true.
      call util_nwchem_srcdir(qmesh_data_name)
      qmesh_data_name
     >     =qmesh_data_name(1:inp_strlen(qmesh_data_name))
     >     //"/nwpw/pspw/lib/exchange-correlation/vdw-DF/vdw_qmesh.dat"


 99   continue

      if (from_environment) then
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: environment'
      else if (from_nwchemrc) then
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: .nwchemrc'
      else
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: compiled reference'
      endif
      if (debug) then
         write(luout,*) ' NWCHEM_QMESH_DATA set to: <',
     >    qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
      end if

      if(.not.util_io_unit(80,90,unitf))
     >  call errquit("vdw-DF cannot get io unit",0,DISK_ERR)


      open(unit=unitf,file=qmesh_data_name,status='old',
     >     form='formatted',ERR=30)

      read(unitf,*,ERR=99,END=999) Nqs
      value = BA_alloc_get(mt_dbl,Nqs,'qmesh', qmesh(2), qmesh(1))
      value = value.and.BA_alloc_get(mt_dbl,Nqs,'ya', ya(2), ya(1))
      value = value.and.BA_alloc_get(mt_dbl,Nqs,'y2a', y2a(2), y2a(1))
      if (.not.value)
     >    call errquit('vdw_DF_init: out of heap, Nqs=',Nqs, MA_ERR)

      do i=1,Nqs
         read(unitf,*,ERR=999,END=999) dbl_mb(qmesh(1)+i-1)
      end do
      close(unitf)

*    **** initialize interpoation splines ****
      value = BA_push_get(mt_dbl,3*Nqs,'utmp', utmp(2), utmp(1))
      if (.not.value) call errquit('vdw_DF_init:out of stack',1,MA_ERR)
      call vdw_DF_init_poly(Nqs,
     >                      dbl_mb(qmesh(1)),
     >                      dbl_mb(ya(1)),
     >                      dbl_mb(y2a(1)),
     <>                     dbl_mb(utmp(1)))
      value = BA_pop_stack(utmp(2))
      if (.not.value) call errquit('vdw_DF_init:pop stack',2,MA_ERR)

      return

 999  continue
      call errquit('vdw_DF_init: error reading qmesh, Nqs=',Nqs,DISK_ERR)
      return

      end


*     **********************************************
*     *                                            *
*     *                vdw_DF_end                  *
*     *                                            *
*     **********************************************

      subroutine vdw_DF_end()
      implicit none

*     **** qmesh common block ****
      logical has_vdw
      integer Nqs,qmesh(2),ya(2),y2a(2)
      common /vdw_df_common/ qmesh,ya,y2a,Nqs,has_vdw

*     **** local variables ****
      logical value

      if (has_vdw) then
         value = BA_free_heap(qmesh(2))
         value = value.and.BA_free_heap(ya(2))
         value = value.and.BA_free_heap(y2a(2))
         if (.not.value) 
     >      call errquit('vdw_DF_end:free heap failed',0,MA_ERR)
      end if

      return
      end

*     **********************************************
*     *                                            *
*     *                vdw_DF_init_poly            *
*     *                                            *
*     **********************************************
      subroutine vdw_DF_init_poly(n,x,ya,y2a,utmp)
      implicit none
      integer n
      real*8 x(*),ya(n,*),y2a(n,*),utmp(*)

*     **** local variables ****
      integer i,j
      real*8 y1,yp

      y1 = 0.0d0
      yp = 0.0d0

      do i=1,n
         do j=1,n
            if (i.eq.j) then
               ya(i,j) = 1.0d0
            else
               ya(i,j) = 0.0d0
            end if
         end do
      end do
      do j=1,n
         call nwpw_spline(x,ya(1,j),n,yp1,ypn,y2a(1,j),utmp)
      end do

      return
      end


      subroutine langreth_kernel_generation()
      implicit none

c     **** Langreth kernel common block ****
      integer Na
      parameter (Na=256)
      !real*8 amin,amax,a(Na),aweights(Na),a2(Na),cos_a(Na),sin_a(Na)
      !real*8 nu(Na),nu1(Na)
      !real*8 Wab(Na,Na)
      real*8 amin,amax
      integer a(2),aweights(2),a2(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2),Wab(2)
      common /lang_phi_common/ Wab,a,aweights,a2,cos_a,sin_a,nu,nu1,
     >                         amin,amax

      return
      end

      subroutine vdw_DF_init_poly(Nq,qmesh)
      implicit none
      integer Nq
      real*8 qmesh(Nq)

      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grx',grx(2),grx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'gry',gry(2),gry(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grz',grz(2),grz(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'agr',agr(2),agr(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'fn',fn(2),fn(1))
        rhog(1) = fn(1)
        value = value.and.
     >        BA_push_get(mt_dbl, 2*n2ft3d,'fdn',fdn(2),fdn(1))
        tmp(1) = fdn(1)
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      return
      end 

      subroutine vdw_DF_init0(Na,amin,amax,a,a2,aweights,
     >                        cos_a,sin_a,nu,nu1,Wab)
      implicit none
      integer Na
      real*8 amin,amax,a(Na),aweights(Na),a2(Na),cos_a(Na),sin_a(Na)
      real*8 nu(Na),nu1(Na),Wab(Na,Na)

      integer, parameter :: Nqs = 20

       real(dp), dimension(Nqs):: q_mesh = (/  1.00D-5, 0.0449420825586261D0, 0.0975593700991365D0, &
       0.159162633466142D0, 0.231286496836006D0, 0.315727667369529D0, 0.414589693721418D0, &
       0.530335368404141D0, 0.665848079422965D0, 0.824503639537924D0, 1.010254382520950D0, &
       1.227727621364570D0, 1.482340921174910D0, 1.780437058359530D0, 2.129442028133640D0, &
       2.538050036534580D0, 3.016440085356680D0, 3.576529545442460D0, 4.232271035198720D0, &
       5.0D0 /)


      integer i,j
      real*8  D,delta,d1,d2,fourpi
      real*8  phi0,phi5,phi9

*     **** external functions ****
      real*8   phi_value
      external phi_value

      fourpi = 16.0d0*datan(1.0d0)
      amin = 0.0d0
      amax = 64.0d0

      call GaussLegendre(datan(amin),datan(amax),Na,a,aweights)
      do i=1,Na
         a(i) = dtan(a(i))
         a2(i) = a(i)**2
         aweights(i) = aweights(i)*(1.0d0+a2(i))
         cos_a(i) = dcos(a(i))
         sin_a(i) = dsin(a(i))
      end do
      do i = 1, Na
      do j = 1, Na
         Wab(i,j) = 2.0d0 * aweights(i)*aweights(j) *  (
     >       (3.0d0-a2(i))*a(j)*cos_a(j)*sin_a(i)  +
     >       (3.0d0-a2(j))*a(i)*cos_a(i)*sin_a(j)  +
     >       (a2(i)+a2(j)-3.0d0)*sin_a(i)*sin_a(j) -
     >       3.0d0*a(i)*a(j)*cos_a(i)*cos_a(j) )   / (a(i)*a(j))
      end do
      end do


      ! The total number of phi_alpha_beta functions that have to be calculated
      Ntotal = (Nqs**2 + Nqs)/2

      dr = (r_max)/(Nr_points)
      do idx = 1, my_Nqs
      do i = 1, Nr_points

        d1 = q_mesh(indices(idx+my_start_q-1, 1)) * (dr * i)
        d2 = q_mesh(indices(idx+my_start_q-1, 2)) * (dr * i)

        phi(r_i, idx) = phi_value(d1, d2)

      end do
      end do
     



      return
      end

      subroutine vdw_DF_energy(dn,gdn)
      implicit none
      real*8 dn(*),gdn(*)
