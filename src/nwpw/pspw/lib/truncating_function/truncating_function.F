c
c $Id: truncating_function.F,v 1.5 2002-03-05 22:09:23 bylaska Exp $
c

***********************************************************
*                                                         *
*            Truncating Function module                   *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*   -- developed by Eric J. Bylaska on February 27, 2002  *
*                                                         *
***********************************************************
*
*
*

*     **********************************
*     *                                *
*     *   Trunctating_Function_init    *
*     *                                *
*     **********************************

      subroutine Truncating_Function_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "truncating_function.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
c      integer nrho
c      integer rho(2),f(2)

      logical value,found
      integer npack0,nfft3d,G(3),nfft(3),n2ft3d
      integer i,ii,ia,l,j,k,zero,pzero,qzero,nxh,nyh
      real*8  gg,rr,rcut_in,pp_in
      real*8  unita(3,3),r,f,x,y,z
      integer tmp1(2),tmp2(2),rgrid(2)
      character*4  element
      character*20 fname
      character*255 full_filename
      integer nx,ny,nz
      real*8  scal1




*     **** external functions ****
      logical      nwpw_filefind
      character*2 ion_atom
      integer  G_indx,ion_nkatm,ion_atn,control_ngrid
      real*8   Truncating_Function_transform,control_unita
      real*8   Truncating_Function_transform0
      external nwpw_filefind
      external ion_atom
      external G_indx,ion_nkatm,ion_atn,control_ngrid
      external Truncating_Function_transform,control_unita
      external Truncating_Function_transform0


*     ----- van der waals radii in angstroms -Needs to be updated   -----
*
      real*8 vander(106)
      data (vander(i),i=1,106)
     1 /1.20d0,1.20d0,
     2 1.37d0,1.45d0,1.45d0,1.50d0,1.50d0,1.40d0,1.35d0,1.30d0,
     3 1.57d0,1.36d0,1.24d0,1.17d0,1.80d0,1.75d0,1.70d0,1.80d0,
     4  18*1.80d0,
     5  18*1.80d0,
     6  32*1.80d0,
     7  20*1.80d0/



      call Parallel_taskid(taskid)

      if (.not.rtdb_get(rtdb,'nwpw:truncating_rcut_scaling',
     >                  mt_dbl,1,rcut))
     >   rcut = 1.2d0

      if (.not.rtdb_get(rtdb,'nwpw:truncating_power',mt_dbl,1,pp))
     >   pp = 6.0d0

      if (taskid.eq.MASTER) then
          write(*,1160)
          write(*,1161) rcut
          write(*,1162) pp
      end if
 1160 FORMAT(/' Truncating Function definition:')
 1161 FORMAT(5X, '- Truncating rcut scaling  = (', E11.3,' )  ')
 1162 FORMAT(5X, '- Truncating power         = (', E11.3,' )  ')


*     **** define Truncating Functions ****      
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)



      nt = ion_nkatm()

      value = MA_alloc_get(mt_dcpl,(nt*nfft3d),'Trnc',Trnc(2),Trnc(1))
      if (.not. value) 
     >  call errquit('Truncating_Function_init: out of heap memory',0)

      value = MA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,n2ft3d,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*n2ft3d,'rgrid',rgrid(2),rgrid(1))
      if (.not. value) 
     >  call errquit('Truncating_Function_init: pushing stack',2)


      call lattice_r_grid(dbl_mb(rgrid(1)))



*     ***** Generate Fourier transform of Truncating Functions *****
      do ia=1,nt
        ii = ion_atn(ia)
        rr = (vander(ii)/0.529177d0)*rcut
        if (taskid.eq.MASTER) then
          write(*,1163) ion_atom(ia),rr,(vander(ii)/0.529177d0)
        end if
 1163   FORMAT(5X,A2,' radii = ',E11.3,' (', E11.3,' van der Waal)')

*       **** define formatted psp name ****
        element = '    '
        element = ion_atom(ia)
        l = index(element,' ') - 1
        fname = element(1:l)//'.trnc'

        found = .false.

        if (nwpw_filefind(fname)) then

*         ****************************************************
*         **** read Truncating function from  binary file ****
*         ****************************************************
          if (taskid.eq.MASTER) then
            call util_file_name_noprefix(fname,.false.,
     >                                .false.,
     >                          full_filename)
            l = index(full_filename,' ') - 1
            call openfile(5,full_filename,l,'r',l)
            call iread(5,nfft,3)
            call dread(5,unita,9)
            call dread(5,rcut_in,1)
            call dread(5,pp_in,1)
          end if
          call D3dB_r_read(1,5,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)))

*         *** close fname binary file ***
          if (taskid.eq.MASTER) then
             call closefile(5)
          end if

          found = .true.
          if ((nfft(1).ne.control_ngrid(1)) .or.
     >        (nfft(2).ne.control_ngrid(2)) .or.
     >        (nfft(3).ne.control_ngrid(3)) .or.
     >        (unita(1,1).ne.control_unita(1,1)) .or.
     >        (unita(2,1).ne.control_unita(2,1)) .or.
     >        (unita(3,1).ne.control_unita(3,1)) .or.
     >        (unita(1,2).ne.control_unita(1,2)) .or.
     >        (unita(2,2).ne.control_unita(2,2)) .or.
     >        (unita(3,2).ne.control_unita(3,2)) .or.
     >        (unita(1,3).ne.control_unita(1,3)) .or.
     >        (unita(2,3).ne.control_unita(2,3)) .or.
     >        (unita(3,3).ne.control_unita(3,3)) .or.
     >        (rcut_in.ne.rcut).or.
     >        (pp_in.ne.pp)) found = .false.
           
        end if

        if (.not.found) then

          if (taskid.eq.MASTER) then
            write(*,1164) ion_atom(ia)
          end if
 1164     FORMAT(5X,A2,' - generating trunctation function')

c*         ***** find the G==0 point in the lattice *****
c          call D3dB_nx(1,nx)
c          call D3dB_ny(1,ny)
c          nxh=nx/2
c          nyh=ny/2
c          i=0
c          j=0
c          k=0
c          call D3dB_ktoqp(1,k+1,qzero,pzero)
c          zero = (qzero-1)*(nxh+1)*ny
c     >         + j*(nxh+1)
c     >         + i+1
c
c          do i = 1,nfft3d
c            gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
c     >            + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
c     >            + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )
cc
c            if ((pzero.eq.taskid) .and. (i.eq.zero)) then
c             dbl_mb(tmp1(1)+i-1) = 
c     >         Truncating_Function_transform0(rr,pp)
c               pp_in = (16.0d0*datan(1.0d0)/3.0d0)*rr**3
c               write(*,*) "gg:",gg,dbl_mb(tmp1(1)+i-1),pp_in
c            else
c             dbl_mb(tmp1(1)+i-1) = 
c     >         Truncating_Function_transform(gg,rr,pp)
c            end if
c
c          end do
           do i=1,n2ft3d
              x = dbl_mb(rgrid(1) + 3*(i-1))
              y = dbl_mb(rgrid(1) + 3*(i-1)+1)
              z = dbl_mb(rgrid(1) + 3*(i-1)+2)
              r = dsqrt(x*x + y*y + z*z)
              f = (1.0d0-(1.0d0-dexp(-(r/rr)**(pp+2)))**pp)
              dbl_mb(tmp1(1)+i-1) = f
              if (f.gt.0.1d0) write(*,*) "i,f:",i,f
           end do

*         **************************************************
*         **** save Truncating function in  binary file ****
*         **************************************************
          if (taskid.eq.MASTER) then
            call util_file_name_noprefix(fname,.false.,
     >                                .false.,
     >                          full_filename)
            l = index(full_filename,' ') - 1
            nfft(1) = control_ngrid(1)
            nfft(2) = control_ngrid(2)
            nfft(3) = control_ngrid(3)

            unita(1,1) = control_unita(1,1)
            unita(2,1) = control_unita(2,1)
            unita(3,1) = control_unita(3,1)
            unita(1,2) = control_unita(1,2)
            unita(2,2) = control_unita(2,2)
            unita(3,2) = control_unita(3,2)
            unita(1,3) = control_unita(1,3)
            unita(2,3) = control_unita(2,3)
            unita(3,3) = control_unita(3,3)

            call openfile(6,full_filename,l,'w',l)
            call iwrite(6,nfft,3)
            call dwrite(6,unita,9)
            call dwrite(6,rcut,1)
            call dwrite(6,pp,1)
          end if
          call D3dB_r_write(1,6,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)))

*         *** close fname binary file ***
          if (taskid.eq.MASTER) then
             call closefile(6)
          end if

        end if
        
*       ***** pack Truncation function for use ****
        call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
        call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
c        call Pack_c_pack(0,dbl_mb(tmp1(1)))
c        call Pack_c_SMul(0,scal1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
        call D3dB_c_SMul(1,scal1,dbl_mb(tmp1(1)),
     >                     dcpl_mb(Trnc(1)+(ia-1)*nfft3d))

      end do

      value =           MA_pop_stack(rgrid(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) 
     >  call errquit('Truncating_Function_init: popping stack',3)

      return
      end

*     **********************************
*     *                                *
*     *   Truncating_Function_end      *
*     *                                *
*     **********************************

      subroutine Truncating_Function_end()
      implicit none

#include "mafdecls.fh"
#include "truncating_function.fh"

      if (.not. MA_free_heap(Trnc(2))) 
     >  call errquit('Truncating_Function_end: freeing heap',0)
      return
      end

*     ****************************
*     *                          *
*     *   Truncating_Function    *
*     *                          *
*     ****************************
      subroutine Truncating_Function(T)
      implicit none
      complex*16 T(*)

#include "mafdecls.fh"
#include "truncating_function.fh"

*     ***** local variables ****
      logical value
      integer ii,ia
      integer nion,nfft3d,npack0
      integer exi(2)

*     **** external functions ****
      integer  ion_nion,ion_katm
      external ion_nion,ion_katm

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      if (.not. value) 
     >  call errquit('Truncating_Function: pushing stack',0)

      call Pack_c_Zero(0,T)
      call dcopy(2*nfft3d,0.0d0,0,T,1)

      nion = ion_nion()
      nion = 1  ! debug
      do ii=1,nion
        ia = ion_katm(ii)

*       **** structure factor ****
c        call strfac(ii,dcpl_mb(exi(1)))
c        call Pack_c_pack(0,dcpl_mb(exi(1)))
c
c*       **** add to Total Truncating Function ****
c        call Pack_cc_Mul(0,dcpl_mb(Trnc(1)+nfft3d*(ia-1)),
c     >                   dcpl_mb(exi(1)),
c     >                   dcpl_mb(exi(1)))
c        call Pack_cc_Sum(0,T,dcpl_mb(exi(1)),T)
         call Pack_cc_Sum(0,
     >                    T,
     >                    dcpl_mb(Trnc(1)+nfft3d*(ia-1)),
     >                    T)

      end do

      value = MA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('Truncating_Function: popping stack',1)
      return
      end

*     ********************************
*     *   	                         *
*     *   grad_Truncating_Function   *
*     *                              *
*     ********************************
      subroutine grad_Truncating_Function(Amm,fion)
      implicit none
      complex*16 Amm(*)
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "truncating_function.fh"

*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia,nion
      integer exi(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) 
     >  call errquit('grad_Truncating_Function: pushing stack',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))


      nion = ion_nion()
      nion = 1  ! debug
      do ii=1,nion
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(0,dcpl_mb(exi(1)))

*       **** add to local psp ****
        call Pack_cc_Mul(0,dcpl_mb(Trnc(1)+nfft3d*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(exi(1)))

         do i=1,npack0
           dbl_mb(xtmp(1)+i-1)
     >        = dimag(Amm(i))* dble(dcpl_mb(exi(1)+i-1))
     >         - dble(Amm(i))*dimag(dcpl_mb(exi(1)+i-1))
        end do

        call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
        call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
        call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))

      end do
      value = MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      value = value.and.MA_pop_stack(xtmp(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('grad_Truncating_Function: popping stack',1)

      return
      end 


*     *************************************
*     *                                   *
*     *   Truncating_Function_transform   *
*     *                                   *
*     *************************************

*     This function returns the fourier transform of
*
*            T(r) = (1.0d0-(1.0d0-dexp(-(r/rcut)**(pp+2)))**pp)
*
*      May want to use:
*            T(r) = exp(-(r/rcut)**pp)
*
*     Entry - gg: g squared
*             rcut:
*             pp:
*
*     Exit - returns
*                      /infty
*                     |
*      T(g) = (4*pi)* | r**2 * T(r)* j0(gr) dr
*                     |
*                     / 0

      real*8 function Truncating_Function_transform(gg,rcut,pp)
      implicit none
      real*8 gg,rcut,pp

*     **** local variables ****
      integer nrho
      parameter (nrho=5001)

      integer i,no,ne
      real*8 q,fourpi,r,sum,drho,pp2,f

      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho-1)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

      q  = dsqrt(gg)
c      sum = 0.0d0
c      do i=2,nrho-1
c        r = (i-1)*drho
c        sum = sum 
c     >    + r*sin(q*r)*(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c      end do
c      r = drho*(nrho-1)
c      sum = sum + 0.5d0*r*sin(q*r)
c     >     *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)


*     **** simpson integration ****
c      no = nrho/2
c      ne = no+1
c
c      r = drho*(nrho-1)
c      f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c      sum = -0.5d0*r*r*f*sin(q*r)
c      do i=1,no
c         r = 2*(i-1)*drho
c         f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c         sum = sum + 2.0d0*r*r*f*sin(q*r)
c      end do
c      do i=1,ne
c         r = (2*(i-1)+1)*drho
c         f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c         sum = sum + 4.0d0*r*r*f*sin(q*r)
c      end do
c      sum = sum/3.0d0
c      Truncating_Function_transform = (fourpi/q)*sum*drho

      Truncating_Function_transform = 
     >  (fourpi/q**2)*(sin(q*rcut)/q - rcut*cos(q*rcut))
      return
      end




*     *************************************
*     *                                   *
*     *   Truncating_Function_transform0  *
*     *                                   *
*     *************************************

*     This function returns the fourier transform of
*
*            T(r) = (1.0d0-(1.0d0-dexp(-(r/rcut)**(pp+2)))**pp)
*
*      May want to use:
*            T(r) = exp(-(r/rcut)**pp)
*
*     Entry - gg: g squared
*             rcut:
*             pp:
*
*     Exit - returns
*                        /infty
*                       |
*      T(g=0) = (4*pi)* | r**2 * T(r)  dr
*                       |
*                       / 0

      real*8 function Truncating_Function_transform0(rcut,pp)
      implicit none
      real*8 rcut,pp

*     **** local variables ****
      integer nrho,no,ne
      parameter (nrho=5001)

      integer i
      real*8 fourpi,r,sum,drho,pp2,f


      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho-1)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

c      sum = 0.0d0
c      do i=2,nrho-1
c            r = (i-1)*drho
c            f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c            sum = sum + r*r*f
c      end do

*     **** simpson integration ****
c     no = nrho/2
c     ne = no+1
c
c     r = drho*(nrho-1)
c     f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c     sum = -0.5d0*r*r*f
c     do i=1,no
c        r = 2*(i-1)*drho
c        f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c        sum = sum + 2.0d0*r*r*f
c     end do
c     do i=1,ne
c        r = (2*(i-1)+1)*drho
c        f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c        sum = sum + 4.0d0*r*r*f
c     end do
c     sum = sum/3.0d0
c      Truncating_Function_transform0 = (fourpi)*sum*drho
       Truncating_Function_transform0 = (fourpi/3.0d0)*(rcut**3)

      return
      end




