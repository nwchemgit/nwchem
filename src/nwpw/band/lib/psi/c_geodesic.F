* $Id: c_geodesic.F,v 1.9 2006-12-05 00:14:40 bylaska Exp $
*


*     ***********************************
*     *			                *
*     *		c_geodesic_init	        *
*     *					*
*     ***********************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value
      integer npack1,nemax,nbrill
      

*     **** external functions **** 
      logical  Pneb_w_allocate
      integer  cpsi_neq,cpsi_nbrillq
      external Pneb_w_allocate
      external cpsi_neq,cpsi_nbrillq


      call Cram_max_npack(npack1)
      nemax  = cpsi_neq(1)+cpsi_neq(2)
      nbrill = cpsi_nbrillq()

      value = MA_alloc_get(mt_dcpl,npack1*nemax*nbrill,
     >                     'U',U(2),U(1))
      value = value.and.Pneb_w_allocate(0,0,Vt)
      value = value.and.
     >        MA_alloc_get(mt_dbl,nemax*nbrill,
     >                     'S',S(2),S(1))
      if (.not. value) 
     > call errquit('c_geodesic_init:out of heap memory',0, MA_ERR)

      return
      end

*     ***********************************
*     *					*
*     *		c_geodesic_finalize	*
*     *					*
*     ***********************************
*
*     Uses - c_geodesic common block
*
      subroutine c_geodesic_finalize()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value

*     **** external functions ****
      logical  Pneb_w_free
      external Pneb_w_free

      value =           MA_free_heap(S(2))
      value = value.and.Pneb_w_free(Vt)
      value = value.and.MA_free_heap(U(2))
      if (.not. value) 
     > call errquit(
     >  'c_geodesic_finalize:error freeing of heap memory',0, MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		c_geodesic_start        *
*     *					*
*     ***********************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_start(A,max_sigma,dE)
      implicit none
      complex*16 A(*)
      real*8     max_sigma,dE

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      logical spin_orb
      integer npack1,neall,nbrill,ispin
      integer nb,ms,i,j,n,nn,ne(2),indx,indxt
      integer V(2)
      real*8 tmp

*     **** external functions **** 
      logical  Pneb_w_push_get,Pneb_w_pop_stack,cpsi_spin_orbit
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      real*8   c_electron_eorbit
      external Pneb_w_push_get,Pneb_w_pop_stack,cpsi_spin_orbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external c_electron_eorbit


      call nwpw_timing_start(10)   
      call Cram_max_npack(npack1)
      ispin  = cpsi_ispin()
      ne(1)  = cpsi_ne(1)
      ne(2)  = cpsi_ne(2)
      spin_orb=cpsi_spin_orbit()
      neall  = ne(1) + ne(2)
      nbrill = brillioun_nbrillioun()
      n      = ne(1)
      nn     = n*n
      if (spin_orb) neall=ne(1)

*     **** allocate tmp space ****
      if (.not.Pneb_w_push_get(0,0,V)) 
     > call errquit('c_geodesic_start: out of stack memory',0,MA_ERR)

*     **** HomeGrown SVD ****
      call dcopy(nbrill*neall,0.0d0,0,dbl_mb(S(1)),1)
      call Pneb_SVD(0,0,npack1, 
     >              A, 
     >              dcpl_mb(U(1)),
     >              dbl_mb(S(1)),
     >              dcpl_mb(V(1)) )

*     **** find max_sigma ****
      max_sigma = -1.0d200
      if (spin_orb) then
        indx=S(1)
        do i=1,n*nbrill
          tmp=dabs( dbl_mb(indx) )
          if ( tmp.gt.max_sigma) then
            max_sigma=tmp
          end if
          indx=indx+1
        end do
      else
        indx=S(1)
        do i=1,neall*nbrill
          tmp=dabs( dbl_mb(indx) ) 
          if ( tmp.gt.max_sigma) then
            max_sigma=tmp
          end if
          indx=indx+1
        end do
      end if

c      low=0
c      hi=npack1*ne(1)*2*nbrill
*     *** calculate Vt ***
      if (spin_orb) then
        do nb=1,nbrill
          do i=1,n
           do j=1,n
             indx  = (i-1) + (j-1)*n + (nb-1)*nn
             indxt = (j-1) + (i-1)*n + (nb-1)*nn
             dcpl_mb(Vt(1)+indxt) = dconjg(dcpl_mb(V(1)+indx))
           end do
          end do
        end do
      else
        do nb=1,nbrill
         do ms=1,ispin
          do i=1,ne(ms)
           do j=1,ne(ms)
             indx  = (i-1) + (j-1)*n + (ms-1)*nn + (nb-1)*2*nn
             indxt = (j-1) + (i-1)*n + (ms-1)*nn + (nb-1)*2*nn
             dcpl_mb(Vt(1)+indxt) = dconjg(dcpl_mb(V(1)+indx))
           end do
          end do
         end do
        end do
      end if

*     **** calculate 2*<A|H|psi> ****
      dE = 2.0d0*c_electron_eorbit(A)
       
*     **** deallocate tmp space ****     
      if (.not.Pneb_w_pop_stack(V)) 
     > call errquit("c_geo_start:error popping stack memory",0,MA_ERR)
     
      call nwpw_timing_end(10)
      return
      end



*     *******************************
*     *		            	    *
*     *		c_geodesic_get	    *
*     *			            *
*     *******************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_get(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical    value,spin_orb
      integer    npack1,neall,nbrill,ispin,ne(2)
      integer    nb,ms,j,k,n,nn,indx1,indx2,indx,shift1,shift2
      real*8     scal1,scal2
      complex*16 zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)


*     **** external functions ****
      logical cpsi_spin_orbit
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_spin_orbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun

      value=MA_set_auto_verify(.true.)
      call nwpw_timing_start(10)
      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)
      call Cram_max_npack(npack1)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1) + ne(2)
      n     = ne(1)
      nn    = n*n
      nbrill = brillioun_nbrillioun()
      spin_orb=cpsi_spin_orbit()
      if (spin_orb) ispin=2
*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*n*n),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

      if (spin_orb) then
        n=ne(1)
        do nb =1,nbrill
  
         do j=1,n
          indx = S(1) + j - 1 + (nb-1)*n
          scal1 = dcos(dbl_mb(indx)*t)
          scal2 = dsin(dbl_mb(indx)*t)

          do k=1,n
           indx= j-1+(k-1)*n
           indx2=j-1+(k-1)*n+(nb-1)*n*n
           dcpl_mb(tmp1(1)+indx)
     >          = scal1 * dcpl_mb(Vt(1)+indx2)
           dcpl_mb(tmp3(1)+indx)
     >          = scal2 * dcpl_mb(Vt(1)+indx2)
          end do
         end do
         indx=(nb-1)*n*n
         value=MA_verify_allocator_stuff()
         call ZGEMM('C','N',n,n,n,one,
     >              dcpl_mb(Vt(1)+indx),n,
     >              dcpl_mb(tmp1(1)),n,
     >              zero,
     >              dcpl_mb(tmp2(1)),n)

         do k=1,n
          indx=1+(k-1)*npack1+(nb-1)*n*npack1*2
          call zcopy(npack1,zero,0.0d0,Ynew(indx),1)
          indx=indx+n*npack1
          call zcopy(npack1,zero,0.0d0,Ynew(indx),1)
          do j=1,n
            indx=tmp2(1)+j-1+(k-1)*n
            indx1=((j-1)+(nb-1)*n*2)*npack1+1 
            indx2=((k-1)+(nb-1)*n*2)*npack1+1  
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(indx),
     >                         Yold(indx1),
     >                         Ynew(indx2))
            indx1=indx1+n*npack1
            indx2=indx2+n*npack1
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(indx),
     >                         Yold(indx1),
     >                         Ynew(indx2))
            indx=tmp3(1)+j-1+(k-1)*n
            indx1=U(1)+((j-1)+(nb-1)*n*2)*npack1
            indx2=1+((k-1+(nb-1)*n*2)*npack1)
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(indx),
     >                         dcpl_mb(indx1),
     >                         Ynew(indx2))
            indx1=indx1+n*npack1
            indx2=indx2+n*npack1
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(indx),
     >                         dcpl_mb(indx1),
     >                         Ynew(indx2))

          end do
         end do

       end do

      else

       do nb =1,nbrill
        do ms = 1,ispin
         shift1 = (ms-1)*ne(1)
     >         + (nb-1)*neall
         shift2 = (ms-1)*nn
     >         + (nb-1)*2*nn

         do j=1,ne(ms)
          scal1 = dcos(dbl_mb(S(1)+j-1+shift1)*t)
          scal2 = dsin(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,ne(ms)
           dcpl_mb(tmp1(1)+j-1+(k-1)*n+ (ms-1)*nn) 
     >          = scal1 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
           dcpl_mb(tmp3(1)+j-1+(k-1)*n+ (ms-1)*nn) 
     >          = scal2 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
          end do

        end do

        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >              dcpl_mb(  Vt(1)+shift2),n,
     >              dcpl_mb(tmp1(1)+(ms-1)*nn),n,
     >              zero,
     >              dcpl_mb(tmp2(1)+(ms-1)*nn),n)


        do k=1,ne(ms)
          call zcopy(npack1,zero,0,Ynew(1+(k-1+shift1)*npack1),1)
          do j=1,ne(ms)
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(tmp2(1)+j-1+(k-1)*n+(ms-1)*nn),
     >                         Yold(1+(j-1+shift1)*npack1),
     >                         Ynew(1+(k-1+shift1)*npack1))
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(tmp3(1)+j-1+(k-1)*n+(ms-1)*nn),
     >                         dcpl_mb(U(1)+(j-1+shift1)*npack1),
     >                         Ynew(1+(k-1+shift1)*npack1))

          end do
        end do

       end do
      end do
      end if
*     **** Orthonormality Check ****

      call Pneb_orthocheckmake(.true.,0,0,npack1,Ynew)     

*     **** deallocate tmp space ****     
      value = MA_pop_stack(tmp3(2))
      value = value.and.
     >        MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end

*     ***********************************
*     *					*
*     *		c_geodesic_transport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine c_geodesic_transport(t,Yold,Ynew)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical    value
      integer    npack1,neall,nbrill,ispin,ne(2)
      integer    nb,ms,j,k,n,nn,shift1,shift2,soshft
      complex*16 zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2),ix,kx,indx1,indx2

*     **** external functions ****
      logical cpsi_spin_orbit
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_spin_orbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun

      call nwpw_timing_start(10)
      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)

      call Cram_max_npack(npack1)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1) + ne(2)
      n     = ne(1)
      nn    = n*n
      nbrill = brillioun_nbrillioun()
      value = MA_SET_AUTO_VERIFY(.TRUE.)
*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*n*n),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

      if (cpsi_spin_orbit()) then



      do nb=1,nbrill
        n=ne(1)
        soshft=n*npack1
        shift1 = (nb-1)*soshft*2
        shift2 = (nb-1)*n*n
        do j=1,n
          scal3 = dbl_mb(S(1)+j-1+(nb-1)*n)
          scal1 = scal3*dsin(scal3*t)
          scal3 = scal3*dcos(scal3*t)
          do k=1,n
           kx=j-1+(k-1)*n
           ix=Vt(1)+kx+shift2
           dcpl_mb(tmp1(1)+kx) 
     >       = scal1 * dcpl_mb(ix)
           dcpl_mb(tmp3(1)+kx) 
     >       = scal3 * dcpl_mb(ix)
          end do
        end do
        value=ma_verify_allocator_stuff()
     
        call ZGEMM('C','N',n,n,n,one,
     >              dcpl_mb(Vt(1)+shift2), n,
     >              dcpl_mb(tmp1(1)), n,
     >              zero,
     >              dcpl_mb(tmp2(1)), n)

        do k=1,n
          ix=1+((nb-1)*n*2+(k-1))*npack1
          call zcopy(npack1,zero,0,Ynew(ix),1)
          call zcopy(npack1,zero,0,Ynew(ix+n*npack1),1)
          do j=1,n

            indx1=((j-1)+(nb-1)*n*2)*npack1+1 
            indx2=((k-1)+(nb-1)*n*2)*npack1+1  
             
            call Cram_cc_zaxpy(nb,
     >              (-1.0d0*dcpl_mb(tmp2(1)+j-1+(k-1)*n)),
     >              Yold(indx1),
     >              Ynew(indx2))
            call Cram_cc_zaxpy(nb,
     >              (-1.0d0*dcpl_mb(tmp2(1)+j-1+(k-1)*n)),
     >              Yold(indx1+soshft),
     >              Ynew(indx2+soshft))
            call Cram_cc_zaxpy(nb,
     >              dcpl_mb(tmp3(1)+j-1+(k-1)*n),
     >              dcpl_mb(U(1)+indx1-1),
     >              Ynew(indx2))
            call Cram_cc_zaxpy(nb,
     >              dcpl_mb(tmp3(1)+j-1+(k-1)*n),
     >              dcpl_mb(U(1)+indx1-1+soshft),
     >              Ynew(indx2+soshft))

          end do
        end do

      end do


      else


      do nb=1,nbrill
      do ms = 1,ispin
        shift1 = (ms-1)*ne(1) 
     >         + (nb-1)*neall
        shift2 = (ms-1)*nn
     >         + (nb-1)*nn*2

        do j=1,ne(ms)
          scal1 = dbl_mb(S(1)+j-1+shift1)
     >            *dsin(dbl_mb(S(1)+j-1+shift1)*t)
          scal3 = dbl_mb(S(1)+j-1+shift1)
     >            *dcos(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,ne(ms)
           dcpl_mb(tmp1(1)+j-1+(k-1)*n+(ms-1)*nn) 
     >       = scal1 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
           dcpl_mb(tmp3(1)+j-1+(k-1)*n+(ms-1)*nn) 
     >       = scal3 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
          end do
        end do

        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >              dcpl_mb(  Vt(1)+shift2), n,
     >              dcpl_mb(tmp1(1)+(ms-1)*nn), n,
     >              zero,
     >              dcpl_mb(tmp2(1)+(ms-1)*nn), n)


        do k=1,ne(ms)
          call zcopy(npack1,zero,0,Ynew(1+(k-1+shift1)*npack1),1)
          do j=1,ne(ms)
            call Cram_cc_zaxpy(nb,
     >              (-1.0d0*dcpl_mb(tmp2(1)+j-1+(k-1)*n+(ms-1)*nn)),
     >              Yold(1+(j-1+shift1)*npack1),
     >              Ynew(1+(k-1+shift1)*npack1))
            call Cram_cc_zaxpy(nb,
     >              dcpl_mb(tmp3(1)+j-1+(k-1)*n+(ms-1)*nn),
     >              dcpl_mb(U(1)+(j-1+shift1)*npack1),
     >              Ynew(1+(k-1+shift1)*npack1))

          end do
        end do

      end do
      end do

      end if
*     **** deallocate tmp space ****     
      value = MA_pop_stack(tmp3(2))
      value = value.and.
     >        MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *					*
*     *		c_geodesic_Gtransport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine c_geodesic_Gtransport(t,Yold,tG)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 tG(*)

#include "mafdecls.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      logical    value
      integer    npack1,neall,nbrill,ne(2),ispin
      integer    nb,ms,j,k,n,nn,shift1,shift2,soshift
      complex*16 zero,one
      real*8     scal1,scal3,thingie
      integer    tmp1(2),tmp2(2),tmp3(2),indx1,indx2,indxt

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      logical cpsi_spin_orbit
      external cpsi_spin_orbit

      call nwpw_timing_start(10)

      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)

      call Cram_max_npack(npack1)
      ispin  = cpsi_ispin()
      ne(1)  = cpsi_ne(1)
      ne(2)  = cpsi_ne(2)
      neall  = ne(1) + ne(2)
      n      = ne(1)
      nn     = n*n
      nbrill = brillioun_nbrillioun()


*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*nn*nbrill),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*nn*nbrill),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*nn*nbrill),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

*     ***** calculate U^t * G ****

      call Pneb_ffw_multiply(0,0,
     >                  dcpl_mb(U(1)),
     >                  tG,
     >                  npack1,
     >                  dcpl_mb(tmp1(1)))

      call Pneb_w_copy(0,0,
     >           dcpl_mb(tmp1(1)),
     >           dcpl_mb(tmp3(1)))

      if (cpsi_spin_orbit()) then
ccccccccccccccccccccc
      do nb = 1,nbrill
        n=ne(1)
        shift1 =(nb-1)*n
        shift2 = (nb-1)*n*n*2
        soshift= n*npack1
        do j=1,n
          thingie= t*dbl_mb(S(1)+j-1+shift1)
          scal1 =       dsin(thingie)
          scal3 = 1.0d0 - dcos(thingie)

          do k=1,n
           indx1=j-1+(k-1)*n+(nb-1)*n*n
           dcpl_mb(tmp1(1)+indx1) 
     >        = scal1 * dcpl_mb(tmp1(1)+indx1)
           dcpl_mb(tmp3(1)+indx1) 
     >        = scal3 * dcpl_mb(tmp3(1)+indx1)
          end do
        end do

        indx1=(nb-1)*n*n
        call ZGEMM('C','N',n,n,n,one,
     >              dcpl_mb(  Vt(1)+indx1), n,
     >              dcpl_mb(tmp1(1)+indx1), n,
     >              zero,
     >              dcpl_mb(tmp2(1)+indx1), n)

        do k=1,n
          do j=1,n
            indxt=tmp2(1)+j-1+(k-1)*n+(nb-1)*n*n
            indx1=1+((j-1)+(nb-1)*n*2)*npack1
            indx2=1+((k-1)+(nb-1)*n*2)*npack1
            
            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(indx2)),
     >             Yold(indx1),
     >               tG(indx2))
            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(indx2)),
     >             Yold(indx1+soshift),
     >               tG(indx2+soshift))
            indxt=tmp3(1)+(j-1)+(k-1)*n+(nb-1)*n*n
            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(indxt)),
     >             dcpl_mb(U(1)+indx1-1),
     >             tG(indx2))
            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(indxt)),
     >             dcpl_mb(U(1)+indx1-1+soshift),
     >             tG(indx2+soshift))
          end do
        end do

      end do
cccccccccccccccc

      else

cccccccccccccccccccc

      do nb = 1,nbrill
      do ms = 1,ispin
        shift1 = (ms-1)*ne(1)
     >         + (nb-1)*neall
        shift2 = (ms-1)*nn 
     >         + (nb-1)*2*nn

        do j=1,ne(ms)
          scal1 =       dsin(dbl_mb(S(1)+j-1+shift1)*t)
          scal3 = 1.0d0-dcos(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,ne(ms)
           dcpl_mb(tmp1(1)+j-1+(k-1)*n+shift2) 
     >        = scal1 * dcpl_mb(tmp1(1)+j-1+(k-1)*n + shift2)
           dcpl_mb(tmp3(1)+j-1+(k-1)*n+shift2) 
     >        = scal3 * dcpl_mb(tmp3(1)+j-1+(k-1)*n + shift2)
          end do
        end do

        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >              dcpl_mb(  Vt(1)+shift2), n,
     >              dcpl_mb(tmp1(1)+shift2), n,
     >              zero,
     >              dcpl_mb(tmp2(1)+shift2), n)

        do k=1,ne(ms)
          do j=1,ne(ms)

            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(tmp2(1)+j-1+(k-1)*n+shift2)),
     >             Yold(1+(j-1+shift1)*npack1),
     >               tG(1+(k-1+shift1)*npack1))
            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(tmp3(1)+j-1+(k-1)*n+shift2)),
     >             dcpl_mb(U(1)+(j-1+shift1)*npack1),
     >             tG(1+(k-1+shift1)*npack1))
          end do
        end do

      end do
      end do
ccccccccccccc
      end if

*     **** deallocate tmp space ****     
      value =           MA_pop_stack(tmp3(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end
