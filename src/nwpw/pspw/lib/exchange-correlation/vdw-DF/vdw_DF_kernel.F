
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
      integer taskid,np,MASTER,pcount,msgtype
      parameter (MASTER=0)

      logical value
      logical lprint,mprint,hprint,debug,does_it_exist
      integer unitf,print_level,i,j,l,k
      character*255 qmesh_data_name,fname

      integer Nqs,qmesh(2),g(2),phir(2),phir0(2),phik(2),phik0(2)
      integer rphi(2),sphi(2),utmp(2),xtmp(2),phik2(2)
      integer nr,nr1,nk,nk1,tcount
      real*8 dr,dk,rmax,kmax,r,qa,qb,fourpi,kk,qave,qdelta,yp1,ypn


c     **** Langreth kernel data ****
      integer Na
      parameter (Na=256)
      !real*8 a(Na),aweights(Na),a2(Na),cos_a(Na),sin_a(Na)
      !real*8 nu(Na),nu1(Na)
      !real*8 Wab(Na,Na)
      integer a(2),a2(2),aweights(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2)
      integer Wab(2)

*     **** external functions ****
      logical  util_io_unit
      external util_io_unit
      character*7 c_index_name
      external    c_index_name

      call Parallel_taskid(taskid)
      call Parallel_np(np)

      fourpi = 16.0d0*datan(1.0d0)

      call util_print_get_level(print_level)
      lprint = (print_low.le.print_level).and.(taskid.eq.MASTER)
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
      call Parallel_Brdcst_ivalue(MASTER,Nqs)

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
      call Parallel_Brdcst_values(MASTER,Nqs,dbl_mb(qmesh(1)))


*     **** initialize r and k grid ****
      !*** r-grid ***
      nr  = 2048
      !nr  = 32768
      nr1 = nr+1
      rmax = 100.0d0
      dr = rmax/dble(nr)

      !*** kgrid - maximum g=64 and gg=4096 ... ***
      nk   = 1024
      !nk   = 16384
      nk1  = nk+1
      kmax = 64.0
      dk = kmax/dble(nk)    

      value= BA_push_get(mt_dbl,nk1,'g',g(2),g(1))
      value= value.and.BA_push_get(mt_dbl,nr1,'phir',phir(2),phir(1))
      value= value.and.BA_push_get(mt_dbl,nr1,'phir0',phir0(2),phir0(1))
      value= value.and.BA_push_get(mt_dbl,nr1,'rphi',rphi(2),rphi(1))
      value= value.and.BA_push_get(mt_dbl,nr1,'sphi',sphi(2),sphi(1))
      value= value.and.BA_push_get(mt_dbl,nr1,'utmp',utmp(2),utmp(1))
      value= value.and.BA_push_get(mt_dbl,nr1,'xtmp',xtmp(2),xtmp(1))

      value= value.and.BA_push_get(mt_dbl,nk1,'phik',phik(2),phik(1))
      value=value.and.BA_push_get(mt_dbl,nk1,'phik0',phik0(2),phik0(1))
      value=value.and.BA_push_get(mt_dbl,nk1,'phik2',phik2(2),phik2(1))
      value= value.and.BA_push_get(mt_dbl,Na,'a',a(2),a(1))
      value= value.and.BA_push_get(mt_dbl,Na,'a2',a2(2),a2(1))
      value= value.and.
     >        BA_push_get(mt_dbl,Na,'aweights',aweights(2),aweights(1))
      value= value.and.
     >        BA_push_get(mt_dbl,Na,'cos_a',cos_a(2),cos_a(1))
      value= value.and.
     >        BA_push_get(mt_dbl,Na,'sin_a',sin_a(2),sin_a(1))
      value= value.and.BA_push_get(mt_dbl,Na,'nu',nu(2),nu(1))
      value= value.and.BA_push_get(mt_dbl,Na,'nu1',nu1(2),nu1(1))
      value= value.and.BA_push_get(mt_dbl,Na*Na,'Wab',Wab(2),Wab(1))
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
      call util_file_name_resolve(filename, .false.)
      l = index(filename,' ') -1
      !if (lprint) write(luout,210) filename(1:l)

      if (taskid.eq.MASTER) then
         call openfile(6,filename,l,'w',l)
         call iwrite(6,Nqs,1)
         call iwrite(6,nk,1)
         call dwrite(6,kmax,1)
         call dwrite(6,dbl_mb(qmesh(1)),Nqs)
         !call dwrite(6,dbl_mb(g(1)),nk1)
      end if

*     **** generate langreth kernels ****
      pcount = 0
      tcount = 0
      do j=1,Nqs
         qb = dbl_mb(qmesh(1)+j-1) 
         do i=1,j
            qa = dbl_mb(qmesh(1)+i-1)
            tcount = tcount + 1
            msgtype = pcount+9
            if (pcount.eq.taskid) then

               !*** generate phik ***
                call vdw_DF_kernel_gen_phir(Na,
     >                  dbl_mb(a(1)),dbl_mb(a2(1)),
     >                  dbl_mb(nu(1)),dbl_mb(nu1(1)),dbl_mb(Wab(1)),
     >                  qa,qb,nr,dr,dbl_mb(phir(1)))

                qave = 0.5d0*(qa+qb)
                qdelta = (qb-qa)/(qa+qb)
                call vdw_DF_kernel_bessel(qave,qdelta,
     >                               nr,dr,dbl_mb(phir(1)),
     >                               nk,dk,dbl_mb(phik(1)),
     >                               dbl_mb(xtmp(1)),dbl_mb(rphi(1)),
     >                               dbl_mb(sphi(1)),dbl_mb(utmp(1)))


               if (pcount.ne.MASTER)
     >            call Parallel_send_values(MASTER,msgtype,
     >                                      nk1,dbl_mb(phik(1)))
            end if

            if (taskid.eq.MASTER) then
               if (pcount.ne.MASTER) then
                  call Parallel_recv_values(pcount,msgtype,
     >                                      nk1,dbl_mb(phik0(1)))
                  yp1 = (dbl_mb(phik0(1)+1)-dbl_mb(phik0(1)))/dk
                  ypn = (dbl_mb(phik0(1)+nk)-dbl_mb(phik0(1)+nk-1))/dk
                  call nwpw_spline(dbl_mb(g(1)),dbl_mb(phik0(1)),nk1,
     >                           yp1,ypn,dbl_mb(phik2(1)),
     >                           dbl_mb(utmp(1)))
                  call dwrite(6,dbl_mb(phik0(1)),nk1)
                  call dwrite(6,dbl_mb(phik2(1)),nk1)

               else
                  yp1 = (dbl_mb(phik(1)+1)-dbl_mb(phik(1)))/dk
                  ypn = (dbl_mb(phik(1)+nk)-dbl_mb(phik(1)+nk-1))/dk
                  call nwpw_spline(dbl_mb(g(1)),dbl_mb(phik(1)),nk1,
     >                           yp1,ypn,dbl_mb(phik2(1)),
     >                           dbl_mb(utmp(1)))
                  call dwrite(6,dbl_mb(phik(1)),nk1)
                  call dwrite(6,dbl_mb(phik2(1)),nk1)
               end if
            end if
            pcount = mod(pcount+1,np)
         end do
      end do


      value = BA_pop_stack(Wab(2))
      value = value.and.BA_pop_stack(nu1(2))
      value = value.and.BA_pop_stack(nu(2))
      value = value.and.BA_pop_stack(sin_a(2))
      value = value.and.BA_pop_stack(cos_a(2))
      value = value.and.BA_pop_stack(aweights(2))
      value = value.and.BA_pop_stack(a2(2))
      value = value.and.BA_pop_stack(a(2))
      value = value.and.BA_pop_stack(phik2(2))
      value = value.and.BA_pop_stack(phik0(2))
      value = value.and.BA_pop_stack(phik(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(utmp(2))
      value = value.and.BA_pop_stack(sphi(2))
      value = value.and.BA_pop_stack(rphi(2))
      value = value.and.BA_pop_stack(phir0(2))
      value = value.and.BA_pop_stack(phir(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(qmesh(2))
      if (.not.value) call errquit('vdw_DF_gen_data:pop stack',0,MA_ERR)

*     *** close vdw binary file ***
      call ga_sync()
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if
      return

 210   FORMAT(' vdw data filename:',A)
 999  continue
      call errquit(
     >   '911 vdw_DF_kernel_gen_data: error reading qmesh, Nqs=',
     >              Nqs,DISK_ERR)
      return

      end

*     **********************************************
*     *                                            *
*     *         vdw_DF_kernel_gen_phir             *
*     *                                            *
*     **********************************************
      subroutine vdw_DF_kernel_gen_phir(Na,a,a2,nu,nu1,Wab,
     >                                  q1,q2,nr,dr,phir)
      implicit none
      integer Na
      real*8 a(Na),a2(Na),nu(Na),nu1(Na),Wab(Na,Na)
      real*8 q1,q2
      integer nr
      real*8 dr
      real*8 phir(0:nr)

*     *** local variables ***
      integer i
      real*8 r,d1,d2,qdelta

*     *** external functions ****
      real*8   vdw_DF_kernel_phi_value
      external vdw_DF_kernel_phi_value

      qdelta = (q1-q2)/(q1+q2)
      do i=1,nr
         r = i*dr
         d1 = r * (1.0d0+qdelta)
         d2 = r * (1.0d0-qdelta)
         phir(i) = vdw_DF_kernel_phi_value(Na,a,a2,nu,nu1,Wab,d1,d2)
      end do
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


c      **** local variables ****
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

      vdw_DF_kernel_phi_value = 0.0d0
      !if ((d1.gt.small).or.(d2.gt.small)) then
      if ((d1.eq.0.0d0).and.(d2.eq.0.0d0)) then
         vdw_DF_kernel_phi_value = 0.0d0
      else
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
             if ((w.gt.small).and.(x.gt.small).and.
     >           (y.gt.small).and.(z.gt.small)) then
                T = (1.0d0/(w+x) + 1.0d0/(y+z))
     >             *(1.0d0/((w+y)*(x+z)) + 1.0d0/((w+z)*(y+x)))
                vdw_DF_kernel_phi_value 
     >                = vdw_DF_kernel_phi_value + T * Wab(i,j)
             end if
         end do
         end do
         vdw_DF_kernel_phi_value = vdw_DF_kernel_phi_value/(pi**2)
      end if

      return
      end

*     **********************************************
*     *                                            *
*     *             GaussLegendre                  *
*     *                                            *
*     **********************************************

      subroutine GaussLegendre(amin,amax,Npoints,a,weights)
      implicit none
      real*8 amin,amax
      integer Npoints
      real*8 a(*),weights(*)

c     **** local varaiables ****    
      logical done
      real*8  poly1,poly2,poly3 !! Temporary storage for Legendre Polynomials
      integer i,j,N                !! Indexing variables 
      real*8  root,dpdx,last_root !! The value of the root of a given Legendre polynomial,
      real*8  midpoint,length     !! The middle of the x-range and the length to that point
      real*8  pi

      pi = 4.0d0*datan(1.0d0)
      N = (Npoints + 1)/2
      midpoint = 0.5d0*(amin+amax)
      length   = 0.5d0*(amax-amin)

      do i=1,N
         root = dcos(dble(pi*(i-0.25d0)/(Npoints+0.5d0)))

         done = .false.
         do while (.not.done)
            poly1 = 1.0d0
            poly2 = 0.0d0
            do j = 1, Npoints
               poly3 = poly2
               poly2 = poly1
               poly1 = ((2.0d0*j-1.0d0)*root*poly2-(j-1.0d0)
     >                  *poly3)/j
            end do
            dpdx=Npoints*(root*poly1-poly2)/(root**2-1.0d0)

            last_root = root
            root = last_root - poly1/dpdx
          
            if (dabs(root-last_root).le.1.0d-14) done = .true.
         end do
       
         a(i)           = midpoint-length*root
         a(Npoints+1-i) = midpoint+length*root

         weights(i) = 2.0d0*length/((1.0d0-root**2)*dpdx**2)
         weights(Npoints+1-i) = weights(i)
      end do
    
      return
      end

*     **********************************************
*     *                                            *
*     *             vdw_DF_kernel_bessel           *
*     *                                            *
*     **********************************************

      subroutine vdw_DF_kernel_bessel(q,dq,nr,dr,phir,nk,dk,phik,
     >                                x,rphi,rphipp,utmp)
      implicit none
      real*8  q,dq
      integer nr
      real*8  dr,phir(0:nr)
      integer nk
      real*8  dk,phik(0:nk)
      real*8  x(0:nr),rphi(0:nr),rphipp(0:nr),utmp(0:nr)

      integer i,k
      real*8 r,kk,fourpi
      real*8 rmax,A,yp1,ypn

      real*8   vdw_DF_Fsin
      external vdw_DF_Fsin

      fourpi = 16.0d0*datan(1.0d0)
      call dcopy(nk+1,0.0d0,0,phik,1)

      rmax  = dr*nr
      A     = phir(nr)*rmax**6

      !************************************************************************
      !*** the k=0 term needs to go to zero for the integrals to make sense ***
      !************************************************************************
c      do i=1,nr
c         r = i*dr
c         phik(0) = phik(0) + phir(i)*r**2
c      end do
c      if (dabs(dq).lt.1.0e-5) then
c         phik(0) = 0.0d0
c      else
c         phik(0) = phik(0)*dr + (A/(3.0d0*rmax**3))
c         phik(0) = (fourpi*phik(0)/q**3)
c      end if
c      !write(*,*) "q,dq,A,residual=",q,dq,A,A/(3.0d0*rmax**3)*fourpi
c      !write(*,*) "divide by q**3=",A/(3.0d0*rmax**3)*fourpi/q**3,phik(0)
      phik(0) = 0.0d0

      !*** use a Filon Integration strategy for the other k values ****
      do i=0,nr
         r = i*dr
         rphi(i) = r*phir(i)
         x(i) = r
      end do
      yp1 = (rphi(1)-rphi(0))/dr
      ypn = (rphi(nr)-rphi(nr-1))/dr

      call nwpw_spline(x,rphi,nr+1,yp1,ypn,rphipp,utmp)


      do k=1,nk
         kk = k*dk
         do i=1,nr
            r = i*dr
            phik(k) = phik(k) 
     >               + vdw_DF_Fsin(r-dr,r,
     >                             rphi(i-1),rphi(i),
     >                             rphipp(i-1),rphipp(i),kk/q)
         end do
         phik(k) = (fourpi/(kk*q**2)) * phik(k)
      end do


      return
      end 

*     **********************************************
*     *                                            *
*     *                vdw_DF_Fsin                 *
*     *                                            *
*     **********************************************

      real*8 function vdw_DF_Fsin(x0,x1,y0,y1,ypp0,ypp1,g)
      implicit none
      real*8 x0,x1,y0,y1,ypp0,ypp1,g

      real*8 Asin0, Bsin0,Csin0,Dsin0

      Asin0 = (g*(x0 - x1)*dcos(g*x0) - dsin(g*x0)
     >        + dsin(g*x1))/(g**2*(x0-x1))

      Bsin0 = (-(g*(x0 - x1)*dcos(g*x1)) + dsin(g*x0)
     >        - dsin(g*x1))/(g**2*(x0 - x1))

      Csin0 = -(6.0d0*g*(x0 - x1)*dcos(g*x0)
     >      + 2.0d0*(-3.0d0 + g**2*(x0 - x1)**2)*dsin(g*x0) +
     >     (6.0d0 + g**2*(x0 - x1)**2)*Sin(g*x1))/(6.0d0*g**4*(x0-x1))


      Dsin0 = (-6.0d0*dsin(g*x0) + 6.0d0*dsin(g*x1) + 
     >    g*(x0 - x1)*(6.0d0*dcos(g*x1) - g*(x0 - x1)*(dsin(g*x0)
     >    + 2.0d0*dsin(g*x1))))/ 
     >    (6.0d0*g**4*(x0 - x1))

      vdw_DF_Fsin = Asin0*y0 + Bsin0*y1 + Csin0*ypp0 + Dsin0*ypp1
      return
      end


