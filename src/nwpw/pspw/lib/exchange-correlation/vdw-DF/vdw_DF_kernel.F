
*     *****************************************
*     *                                       *
*     *      vdw_DF_kernel_gen_data           *
*     *                                       *
*     *****************************************
*
      subroutine vdw_DF_kernel_gen_data(filename)
      implicit none
      character*(*) filename

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical mprint,hprint,debug,does_it_exist
      integer unitf,print_level,i,j
      character*255 qmesh_data_name

      integer Nqs,qmesh(2),g(2),phir(2),phik(2)
      integer nr,nr1,nk,nk1
      real*8 dr,dk


c     **** Langreth kernel data ****
      integer Na
      parameter (Na=256)
      !real*8 a(Na),aweights(Na),a2(Na),cos_a(Na),sin_a(Na)
      !real*8 nu(Na),nu1(Na)
      !real*8 Wab(Na,Na)
      integer a(2),aweights(2),a2(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2)
      integer Wab(2)


      call Parallel_taskid(taskid)

      call util_print_get_level(print_level)
      mprint = (print_medium.le.print_level).and.(taskid.eq.MASTER)
      hprint = (print_high  .le.print_level).and.(taskid.eq.MASTER)
      debug  = (print_debug .le.print_level).and.(taskid.eq.MASTER)

*     **** read and allocate qmesh data ****
      if (taskid.eq.MASTER) THEN
         call vdw_DF_get_qmesh_filename(qmesh_data_name)

         if(.not.util_io_unit(80,90,unitf))
     >     call errquit("vdw-DF cannot get io unit",0,DISK_ERR)

         open(unit=unitf,file=qmesh_data_name,status='old',
     >     form='formatted',ERR=999)

         read(unitf,*,ERR=999,END=999) Nqs
      end if

*     **** allocate qmesh ****
      value = BA_push_get(mt_dbl,Nqs,'qmesh',qmesh(2),qmesh(1))
      if (.not.value)
     >   call errquit('vdw_DF_gen_data: out of stack, Nqs=',Nqs, MA_ERR)

      if (taskid.eq.MASTER) then
         do i=1,Nqs
            read(unitf,*,ERR=999,END=999) dbl_mb(qmesh(1)+i-1)
         end do
         close(unitf)
      end if


*     **** initialize r and k grid ****
      !*** r-grid ***
      nr  = 1024
      nr1 = nr+1
      rmax = 100.0d0
      dr = rmax/dble(nr)

      !*** kgrid - maximum g=64 and gg=4096 ... this logic may need to change to reduce the size of the k-grid ***
      nk   = 256
      nk1  = nk+1
      kmax = 64.0
      dk = kmax/dble(nk)    

      value = BA_alloc_get(mt_dbl,nk1,'g',g(2),g(1))
      value = value.and.BA_alloc_get(mt_dbl,nr1,'phir',phir(2),phir(1))
      value = value.and.BA_alloc_get(mt_dbl,nk1,'phik',phik(2),phik(1))
      value = value.and.BA_alloc_get(mt_dbl,Na,'a',a(2),a(1))
      value = value.and.BA_alloc_get(mt_dbl,Na,'a2',a2(2),a2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'aweights',aweights(2),aweights(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'cos_a',cos_a(2),cos_a(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'sin_a',sin_a(2),sin_a(1))
      value = value.and.BA_alloc_get(mt_dbl,Na,'nu',nu(2),nu(1))
      value = value.and.BA_alloc_get(mt_dbl,Na,'nu1',nu1(2),nu1(1))
      value = value.and.BA_alloc_get(mt_dbl,Na*Na,'Wab',Wab(2),Wab(1))
      if (.not.value)
     >   call errquit('vdw_DF_kernel_gen_data: out of stack',3,MA_ERR)


      !**** generate Wab  and other phi kernel variables ****
      call vdw_DF_kernel_gen_Wab(Na,
     >                           dbl_mb(a(1)),
     >                           dbl_mb(a2(1)),
     >                           dbl_mb(aweights(1)),
     >                           dbl_mb(cos_a(1)),
     >                           dbl_mb(sin_a(1)),
     >                           dbl_mb(Wab(1)))


      !**** generate |g| grid ****
      do i=0,nk
         dbl_mb(g(1)+i) = dk*i
      end do


*     **** create vdw data file ****
      full_filename = filename
      call util_file_name_resolve(full_filename, .false.)
      l = index(full_filename,' ') -1
      if (lprint) write(luout,210) full_filename(1:l)

      if (taskid.eq.MASTER) then
         call openfile(6,filename,l,'w',l)
         call iwrite(6,Nqs,1)
         call dwrite(6,dbl_mb(qmesh(1)),Nqs)
         call iwrite(6,nk,1)
         call dwrite(6,kmax,1)
         call dwrite(6,dbl_mb(g(1)),nk1)
      end if


*     **** generate langreth kernels ****
      do j=1,Nqs
         qb = dbl_mb(qmesh(1)+j-1)
         do i=1,Nqs
            qa = dbl_mb(qmesh(1)+i-1)

            !*** generate phik ***
            --- generate phik here ---

            if (taskid.eq.MASTER) then
               call dwrite(6,dbl_mb(phik(1)),nk1)
            end if
         end do
      end do

*     *** close vdw binary file ***
      call ga_sync()
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if
      return

 210   FORMAT(/' vdw data filename:',A)
 999  continue
      call errquit('vdw_DF_kernel_gen_data: error reading qmesh, Nqs=',
     >              Nqs,DISK_ERR)
      return

      end


*     **********************************************
*     *                                            *
*     *         vdw_DF_kernel_gen_Wab              *
*     *                                            *
*     **********************************************

      subroutine vdw_DF_kernel_gen_Wab(Na,a,a2,aweights,cos_a,sin_a,Wab)
      integer Na
      real*8  a(Na),a2(Na),aweights(Na),cos_a(Na),sin_a(Na),Wab(Na,Na)


      **** local variables ****
      integer i,j
      real*8 amin,amax,fourpi

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

      return
      end


*     **********************************************
*     *                                            *
*     *         vdw_DF_kernel_phi_value            *
*     *                                            *
*     **********************************************

      real*8 function vdw_DF_kernel_phi_value(Na,a,a2,nu,nu1,Wab,d1,d2)
      implicit none
      integer Na
      real*8 a(Na),a2(Na)
      real*8 nu(Na),nu1(Na)
      real*8 Wab(Na,Na)
      real*8 d1,d2


*     **** local variables ****
      real*8 small
      parameter (small=1.0d-12)
      integer i,j
      real*8  w,x,y,z,T,pi,gamma,d1s,d2s

      pi    = 4.0d0*datan(1.0d0)
      gamma = 4.0d0*pi/9.0d0
      d1s = d1**2
      d2s = d2**2

      phi_value = 0.0d0
      if ((d1.gt.small).or.(d2.gt.small)) then
         do i=1,Na
            if ((a(i).le.small).and.(d1.gt.small)) then
               nu(i) = 9.0d0/8.0d0*d1s/pi
            else if ((d1.le.small).or.
     >               (((a2(i)*gamma)/d1s).gt.700.0d0)) then
               nu(i) = a2(i)/2.0d0
            else
               nu(i) = a2(i)
     >                /((1.0d0-dexp(-(a2(i)*gamma)/d1s))*2.0d0)
            end if

            if ((a(i).le.small).and.(d2.gt.small)) then
               nu1(i) = 9.0d0/8.0d0*d2s/pi
            else if ((d2.le.small).or.
     >               (((a2(i)*gamma)/d2s).gt.700.0d0)) then
               nu1(i) = a2(i)/2.0d0
            else
               nu1(i) = a2(i)
     >                /((1.0d0-dexp(-(a2(i)*gamma)/d2s))*2.0d0)
            end if
c            write(*,*) "i,nu,nu1=",i,a(i),nu(i),nu1(i),a2(i),
c     >                 (-(a2(i)*gamma)/d1s),
c     >                 (-(a2(i)*gamma)/d2s)
cc     >                 9.0d0/8.0d0*d1s/pi,9.0d0/8.0d0*d2s/pi
         end do
         do i=1,Na
         do j=1,Na
             w = nu(i)
             x = nu(j)
             y = nu1(i)
             z = nu1(j)
             T = (1.0d0/(w+x) + 1.0d0/(y+z))
     >          *(1.0d0/((w+y)*(x+z)) + 1.0d0/((w+z)*(y+x)))
             phi_value = phi_value + T * Wab(i,j)
         end do
         end do
         phi_value = phi_value/(pi**2)
      end if

      return
      end

