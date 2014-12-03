! $Id$

!**************************************************
!
!       Name: alloc_paw_basis_data
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      SUBROUTINE  alloc_paw_basis_data(nt,nb,ng)
      implicit none
#include "errquit.fh"
      integer nt
      integer nb(nt)
      integer ng(nt)

#include "paw_basis_data.fh"
#include "bafdecls.fh"      
#include "paw_ma.fh"


      logical ok
      integer it
      integer offset_nb
      integer offset_ngb
      integer offset_ng

      tot_ntype = nt

*     !*** find total size for the arrays ***    
      do it = 1, tot_ntype
        tot_nbasis     = tot_nbasis     + nb(it)
        tot_ngridbasis = tot_ngridbasis + nb(it)*ng(it)
        tot_ngrid      = tot_ngrid     + ng(it)
      end do

c      ok = BA_set_auto_verify(.TRUE.)

      ok = .TRUE.

      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"nbasis",nbasis)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_cut",i_cut)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"ngrid",ngrid)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"mult_l",mult_l)
      ok = ok .AND. my_alloc(MT_DBL,tot_ntype,"r1",r1)
      ok = ok .AND. my_alloc(MT_DBL,tot_ntype,"rmax",rmax)
      ok = ok .AND. my_alloc(MT_DBL,tot_ntype,"log_amesh",log_amesh)
      ok = ok .AND. my_alloc(MT_DBL,tot_ntype,"sigma",sigma)
      ok = ok .AND. my_alloc(MT_DBL,tot_ntype,"ion_charge",ion_charge) 
      ok = ok .AND. my_alloc(MT_DBL,tot_ntype,"core_charge",core_charge)
      ok = ok .AND. 
     >     my_alloc(MT_DBL,tot_ntype,"core_kin_energy",core_kin_energy)

      ok = ok .AND. my_alloc(MT_INT,tot_nbasis,"pr_n_ps",pr_n_ps)
      ok = ok .AND. my_alloc(MT_INT,tot_nbasis,"pr_n",pr_n)
      ok = ok .AND. my_alloc(MT_INT,tot_nbasis,"orb_l",orb_l)
      ok = ok .AND. my_alloc(MT_DBL,tot_nbasis,"eig",eig)

      ok = ok .AND. my_alloc(MT_DBL,tot_ngridbasis,"phi_ae",phi_ae)
      ok = ok .AND. 
     >     my_alloc(MT_DBL,tot_ngridbasis,"phi_ae_prime",phi_ae_prime)
      ok = ok .AND. my_alloc(MT_DBL,tot_ngridbasis,"phi_ps",phi_ps)
      ok = ok .AND. 
     >     my_alloc(MT_DBL,tot_ngridbasis,"phi_ps_prime",phi_ps_prime)
      ok = ok .AND. my_alloc(MT_DBL,tot_ngridbasis,"prj_ps",prj_ps)
      ok = ok .AND. 
     >     my_alloc(MT_DBL,tot_ngrid,"core_density",core_density)
      ok = ok .AND. 
     >     my_alloc(MT_DBL,tot_ngrid,"ps_core_density",ps_core_density)

      ok = ok .AND.
     >     my_alloc(MT_DBL,tot_ngrid,"core_density_prime",
     >              core_density_prime)
      ok = ok .AND. 
     >     my_alloc(MT_DBL,tot_ngrid,"ps_core_density_prime",
     >              ps_core_density_prime)

      ok = ok .AND. my_alloc(MT_DBL,tot_ngrid,"v_ps",v_ps)
      ok = ok .AND. my_alloc(MT_DBL,tot_ngridbasis,"prj_ps0",prj_ps0)
      ok = ok .AND. my_alloc(MT_DBL,tot_ngrid,"rgrid",rgrid)

      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_nb",i_nb)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_ng",i_ng)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_ngb",i_ngb)

      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_start",i_start)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_end",i_end)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"npoints",npoints)

      if(.not. ok) then
        call errquit("failed to allocate paw_basis data ",0, MA_ERR)
      end if

      do it=1,tot_ntype
        int_mb(nbasis(1)-1 + it) = nb(it)
        int_mb(ngrid(1)-1 + it)  = ng(it)
      end do

      int_mb(i_nb(1) )  = 0
      int_mb(i_ng(1) )  = 0
      int_mb(i_ngb(1))  = 0

      do it=1,tot_ntype-1
        int_mb(i_nb(1) +it)  = int_mb(i_nb(1) + it - 1)  
     >                       + int_mb(nbasis(1)-1 + it)
        int_mb(i_ng(1) +it)  = int_mb(i_ng(1) + it - 1)  
     >                       + int_mb(ngrid(1)-1 + it)
        int_mb(i_ngb(1)+it)  = int_mb(i_ngb(1) + it - 1) 
     >                       + int_mb(nbasis(1)-1 + it)
     >                        *int_mb(ngrid(1)-1 + it)
      end do

      END !SUBROUTINE alloc_paw_basis_data


      SUBROUTINE dealloc_paw_basis_data()

      IMPLICIT NONE
#include "errquit.fh"
#include "paw_ma.fh"
#include "paw_basis_data.fh"

      logical ok

        ok = .true.
        ok = ok .and. my_dealloc(i_start)
        ok = ok .and. my_dealloc(i_end)
        ok = ok .and. my_dealloc(npoints)

        ok = ok .and. my_dealloc(i_ngb)
        ok = ok .and. my_dealloc(i_ng)
        ok = ok .and. my_dealloc(i_nb)
        ok = ok .and. my_dealloc(rgrid)
        ok = ok .and. my_dealloc(prj_ps0)
        ok = ok .and. my_dealloc(v_ps)
        ok = ok .and. my_dealloc(ps_core_density)
        ok = ok .and. my_dealloc(core_density)
        ok = ok .and. my_dealloc(ps_core_density_prime)
        ok = ok .and. my_dealloc(core_density_prime)
        ok = ok .and. my_dealloc(prj_ps)
        ok = ok .and. my_dealloc(phi_ps_prime)
        ok = ok .and. my_dealloc(phi_ps)
        ok = ok .and. my_dealloc(phi_ae_prime)
        ok = ok .and. my_dealloc(phi_ae)
        ok = ok .and. my_dealloc(eig)
        ok = ok .and. my_dealloc(orb_l)
        ok = ok .and. my_dealloc(pr_n)
        ok = ok .and. my_dealloc(pr_n_ps)
        ok = ok .and. my_dealloc(core_kin_energy)
        ok = ok .and. my_dealloc(ion_charge) 
        ok = ok .and. my_dealloc(core_charge) 
        ok = ok .and. my_dealloc(sigma)
        ok = ok .and. my_dealloc(log_amesh)
        ok = ok .and. my_dealloc(rmax)
        ok = ok .and. my_dealloc(r1)
        ok = ok .and. my_dealloc(mult_l)
        ok = ok .and. my_dealloc(ngrid)
        ok = ok .and. my_dealloc(i_cut)
        ok = ok .and. my_dealloc(nbasis)

        if(.not. ok) then
          call errquit("failed to deallocate paw_basis data ",0, MA_ERR)
        end if



      END !SUBROUTINE


c     *********************************
c     *                               *
c     *         paw_radgrid_map       *
c     *                               *
c     *********************************

      subroutine paw_radgrid_map(tot_nr,i_start,i_end)
      implicit none
      integer tot_nr
      integer i_start
      integer i_end

c     *** local variables ***
      integer nr,np,taskid
      real*8  tmp

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      tmp = dble(tot_nr)/dble(np)
      nr = dint(tmp)

      i_start = 1 + taskid*nr

      if (taskid .eq. (np-1)) THEN
        i_end = tot_nr
      else
        i_end = i_start + nr
      end if

      return
      end




!**************************************************
!
!       Name: paw_basis_tot_ntype
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_tot_ntype()
      implicit none

#include "paw_basis_data.fh"

      paw_basis_tot_ntype = tot_ntype
      return
      END 

!**************************************************
!
!       Name: paw_basis_tot_nbasis
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_tot_nbasis()
      implicit none

#include "paw_basis_data.fh"

      paw_basis_tot_nbasis = tot_nbasis
      return
      END 

!**************************************************
!
!       Name: paw_basis_tot_ngrid
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_tot_ngrid()
      implicit none

#include "paw_basis_data.fh"

      paw_basis_tot_ngrid = tot_ngrid
      return
      END
 
!**************************************************
!
!       Name: paw_basis_tot_ngridbasis
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_tot_ngridbasis()
      implicit none

#include "paw_basis_data.fh"
 
      paw_basis_tot_ngridbasis = tot_ngridbasis
      return
      END
 
!**************************************************
!
!       Name: paw_basis_i_nbasis
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_nbasis(it)
      implicit none
      integer it
   
#include "paw_basis_data.fh"


      paw_basis_i_nbasis = nbasis(1) +  it-1
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_nbasis
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
        INTEGER FUNCTION  paw_basis_nbasis(it)
        implicit none
        integer it
     
#include "paw_basis_data.fh"
#include "bafdecls.fh"
  
        paw_basis_nbasis = int_mb(nbasis(1) +  it-1)
        return
        END 
  
!**************************************************
!
!       Name: paw_basis_i_ngrid
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_ngrid(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_ngrid = ngrid(1) +  it-1
      return
      END
 
!**************************************************
!
!       Name: paw_basis_ngrid
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      integer function paw_basis_ngrid(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_ngrid = int_mb(ngrid(1) +  it-1)
      return
      END
 
!**************************************************
!
!       Name:
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_mult_l(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"


      paw_basis_i_mult_l = mult_l(1) +  it-1
      return
      END 

!**************************************************
!
!       Name:
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_mult_l(it)
      implicit none
      integer it
    
#include "bafdecls.fh"      
#include "paw_basis_data.fh"


      paw_basis_mult_l = int_mb(mult_l(1) +  it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_r1
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_r1(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_r1 = r1(1) +  it-1
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_rmax
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_rmax(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"


      paw_basis_i_rmax = rmax(1) +  it-1
      return
      END 
!**************************************************
!
!       Name: paw_basis_i_cut
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_cut(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_i_cut = int_mb(i_cut(1) +  it-1)
      return
      END 


!**************************************************
!
!       Name: paw_basis_i_i_cut
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_i_cut(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_i_cut = i_cut(1) +  it-1
      return
      END 


!**************************************************
!
!       Name: paw_basis_i_i_start
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_i_start(it)
      implicit none
      integer it

#include "paw_basis_data.fh"

      paw_basis_i_i_start = i_start(1) +  it-1
      return
      END

!**************************************************
!
!       Name: paw_basis_i_start
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_start(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_i_start = int_mb(i_start(1) +  it-1)
      return
      END


!**************************************************
!
!       Name: paw_basis_i_i_end
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_i_end(it)
      implicit none
      integer it

#include "paw_basis_data.fh"

      paw_basis_i_i_end = i_end(1) +  it-1
      return
      END


!**************************************************
!
!       Name: paw_basis_i_end
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_end(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_i_end = int_mb(i_end(1) +  it-1)
      return
      END


!**************************************************
!
!       Name: paw_basis_i_npoints
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_npoints(it)
      implicit none
      integer it

#include "paw_basis_data.fh"

      paw_basis_i_npoints = npoints(1) +  it-1
      return
      END


!**************************************************
!
!       Name: paw_basis_npoints
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_npoints(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_npoints = int_mb(npoints(1) +  it-1)
      return
      END



!**************************************************
!
!       Name: paw_basis_i_sigma
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_sigma(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_sigma = sigma(1) +  it-1
      return
      END 

!**************************************************
!
!       Name: paw_basis_log_amesh
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      double precision FUNCTION  paw_basis_log_amesh(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_log_amesh = dbl_mb(log_amesh(1)+it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_log_amesh
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_log_amesh(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_log_amesh = log_amesh(1) +  it-1
      return
      END 

!**************************************************
!
!       Name: paw_basis_core_charge
!
!       Purpose
!
!       Created:        8/06/2002
!**************************************************
      subroutine  calc_paw_basis_core_charge(ia,q)
      implicit none
      integer ia
      double precision q

#include "paw_basis_data.fh"
#include "integrate.fh"
#include "paw_basis.fh"
#include "bafdecls.fh"
      
      !*** local variables ***
      real*8 core,fourpi

c      !*** external functions ***
c      integer  paw_basis_i_core_density,paw_basis_i_rgrid
c      integer  paw_basis_i_log_amesh,paw_basis_i_ngrid
c      external paw_basis_i_core_density,paw_basis_i_rgrid
c      external paw_basis_i_log_amesh,paw_basis_i_ngrid

      fourpi = 16.0d0*datan(1.0d0)
      q = fourpi*def_integr(0,
     >                         dbl_mb(paw_basis_i_core_density(ia)),
     >                         2,
     >                         dbl_mb(paw_basis_i_rgrid(ia)),
     >                         dbl_mb(paw_basis_i_log_amesh(ia)),
     >                         int_mb(paw_basis_i_ngrid(ia)))

      return
      end

!**************************************************
!
!       Name: paw_basis_i_core_charge
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_core_charge(it)
      implicit none
      integer it

#include "paw_basis_data.fh"

      paw_basis_i_core_charge = core_charge(1) +  it-1
      return
      END

!**************************************************
!
!       Name: paw_basis_core_charge
!
!       Purpose
!
!       Created:        7/9/2002
!**************************************************
      DOUBLE PRECISION FUNCTION  paw_basis_core_charge(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_core_charge = dbl_mb(core_charge(1) +  it-1)
      return
      END


!**************************************************
!
!       Name: paw_basis_i_ion_charge
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_ion_charge(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_ion_charge = ion_charge(1) +  it-1
      return
      END 

      
!**************************************************
!
!       Name: paw_basis_ion_charge
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      DOUBLE PRECISION FUNCTION  paw_basis_ion_charge(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_ion_charge = dbl_mb(ion_charge(1) +  it-1)
      return
      END 

      
**************************************************
!
!       Name: paw_basis_sphere_radius
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      real*8 function  paw_basis_sphere_radius(ia)
      implicit none
      integer ia
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      


      !*** external functions ***
      integer  paw_basis_i_rgrid,paw_basis_i_i_cut
      external paw_basis_i_rgrid,paw_basis_i_i_cut

      paw_basis_sphere_radius = dbl_mb(paw_basis_i_rgrid(ia)
     >                                 +int_mb(paw_basis_i_i_cut(ia))-1)
      return
      end 


**************************************************
!
!       Name: paw_basis_sigma
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      real*8 function  paw_basis_sigma(ia)
      implicit none
      integer ia
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      !*** external functions ***
      integer  paw_basis_i_sigma
      external paw_basis_i_sigma 

      paw_basis_sigma = dbl_mb(paw_basis_i_sigma(ia))
      return
      end

**************************************************
!
!       Name: paw_tot_nlm_nbasis
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      integer function  paw_tot_nlm_nbasis()

      implicit none

      integer ia
      integer ii
      integer l
    
#include "paw_geom.fh"

      !*** external functions ***
      integer  paw_nlm_nbasis
      external paw_nlm_nbasis

      paw_tot_nlm_nbasis = 0
      do ia=1,ion_nion()
        paw_tot_nlm_nbasis = paw_tot_nlm_nbasis + paw_nlm_nbasis(ia)
      end do

      return
      end


**************************************************
!
!       Name: paw_nlm_nbasis
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      integer function  paw_nlm_nbasis(ia)

      implicit none

      integer ia
      integer ii
      integer l
    
#include "paw_geom.fh"
#include "bafdecls.fh"      

      !*** external functions ***
      integer  paw_basis_i_orb_l,paw_basis_nbasis
      external paw_basis_i_orb_l,paw_basis_nbasis 

      paw_nlm_nbasis = 0
      do ii=1,paw_basis_nbasis(ia)
         l =  int_mb(paw_basis_i_orb_l(ia)+ii-1)
         paw_nlm_nbasis = paw_nlm_nbasis + 2*l+1
      end do

      return
      end


**************************************************
!
!       Name: paw_basis_n
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      integer function  paw_basis_n(ii,ia)
      implicit none
      integer ii,ia
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      !*** external functions ***
      integer  paw_basis_i_pr_n
      external paw_basis_i_pr_n

      paw_basis_n = int_mb(paw_basis_i_pr_n(ia)+ii-1)
      return
      end
     
**************************************************
!
!       Name: paw_basis_n_ps
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      integer function  paw_basis_n_ps(ii,ia)
      implicit none
      integer ii,ia
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      !*** external functions ***
      integer  paw_basis_i_pr_n_ps
      external paw_basis_i_pr_n_ps 

      paw_basis_n_ps = int_mb(paw_basis_i_pr_n_ps(ia)+ii-1)
      return
      end
     
**************************************************
!
!       Name: paw_basis_orb_l
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      integer function  paw_basis_orb_l(ii,ia)
      implicit none
      integer ii,ia
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      !*** external functions ***
      integer  paw_basis_i_orb_l
      external paw_basis_i_orb_l                                

      paw_basis_orb_l = int_mb(paw_basis_i_orb_l(ia)+ii-1)
      return
      end
     
   
**************************************************
!
!       Name: paw_basis_eig
!
!       Purpose 
!
!       Created:        8/06/2002
!**************************************************
      real*8 function  paw_basis_eig(ii,ia)
      implicit none
      integer ii,ia
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      !*** external functions ***
      integer  paw_basis_i_eig
      external paw_basis_i_eig

      paw_basis_eig = dbl_mb(paw_basis_i_eig(ia)+ii-1)
      return
      end
     

 
!**************************************************
!
!       Name: paw_basis_i_core_kin_energy
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_core_kin_energy(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"

      paw_basis_i_core_kin_energy = core_kin_energy(1) +  it-1
      return
      END 


!**************************************************
!
!       Name: paw_basis_core_kin_energy
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      double precision  FUNCTION  paw_basis_core_kin_energy(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_core_kin_energy = dbl_mb(core_kin_energy(1) +  it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_pr_n
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_pr_n(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_pr_n = pr_n(1) + int_mb(i_nb(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_pr_n_ps
!
!       Purpose 
!
!       Created:        7/9/200
!**************************************************
      INTEGER FUNCTION  paw_basis_i_pr_n_ps(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_pr_n_ps = pr_n_ps(1) + int_mb(i_nb(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_orb_l
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_orb_l(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_orb_l = orb_l(1) + int_mb(i_nb(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_eig
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_eig(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_eig = eig(1) + int_mb(i_nb(1) + it-1)
      return
      END   
  
!**************************************************
!
!       Name: paw_basis_i_core_density
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_core_density(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_core_density = core_density(1) +
     +     int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_ps_core_density
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_ps_core_density(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_ps_core_density = ps_core_density(1) 
     >                      + int_mb(i_ng(1) + it-1)
      return
      END 



!**************************************************
!
!       Name: paw_basis_i_core_density_prime
!
!       Purpose - needed for gga's
!
!       Created:        1/28/2006
!**************************************************
      INTEGER FUNCTION  paw_basis_i_core_density_prime(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_i_core_density_prime = core_density_prime(1)
     >                               + int_mb(i_ng(1)+it-1)
      return
      END

!**************************************************
!
!       Name: paw_basis_i_ps_core_density_prime
!
!       Purpose - needed for gga's
!
!       Created:        1/28/2006
!**************************************************
      INTEGER FUNCTION  paw_basis_i_ps_core_density_prime(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"

      paw_basis_i_ps_core_density_prime = ps_core_density_prime(1)
     >                                  + int_mb(i_ng(1)+ it-1)
      return
      END


!**************************************************
!
!       Name: paw_basis_i_v_ps
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_v_ps(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_v_ps = v_ps(1) + int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_rgrid
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_rgrid(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_rgrid = rgrid(1) + int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_phi_ae
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_phi_ae(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_phi_ae = phi_ae(1) + int_mb(i_ngb(1) + it-1)
      return
      END 


!**************************************************
!
!       Name: paw_basis_i_phi_ps
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_phi_ps(it)
      implicit none
      integer it
    
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_phi_ps = phi_ps(1) + int_mb(i_ngb(1) + it-1)
      return
      END 


!**************************************************
!
!       Name: paw_basis_i_phi_ae_prime
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_phi_ae_prime(it)
      implicit none
      integer it
   
#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_phi_ae_prime = phi_ae_prime(1) + 
     +     int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_phi_ps_prime
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_phi_ps_prime(it)
      implicit none
      integer it

#include "paw_basis_data.fh"
#include "bafdecls.fh"      

      paw_basis_i_phi_ps_prime = phi_ps_prime(1) + 
     +     int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_prj_ps
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_prj_ps(it)
      implicit none
      integer it
 
#include "paw_basis_data.fh"
#include "bafdecls.fh"
      
      paw_basis_i_prj_ps = prj_ps(1) + int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!       Name: paw_basis_i_prj_ps0
!
!       Purpose 
!
!       Created:        7/9/2002
!**************************************************
      INTEGER FUNCTION  paw_basis_i_prj_ps0(it)
      implicit none
      integer it
    
#include "bafdecls.fh"      
#include "paw_basis_data.fh"

      paw_basis_i_prj_ps0 = prj_ps0(1) + int_mb(i_ngb(1) + it-1)
      return
      END 

      subroutine set_max_i_cut(ic)
      implicit none
      integer ic
#include "paw_basis_data.fh"

        max_i_cut = ic

       end



      subroutine set_max_mult_l(l)
      implicit none
      integer l
#include "paw_basis_data.fh"

        max_mult_l = l

       end

      integer function paw_basis_max_i_cut()
      implicit none

#include "paw_basis_data.fh"

        paw_basis_max_i_cut= max_i_cut

      end

      integer function paw_basis_max_mult_l()
      implicit none

#include "paw_basis_data.fh"

        paw_basis_max_mult_l= max_mult_l

      end




c     *************************************************
c     *                                               *
c     *        paw_basis_derivative_ngrid             *
c     *                                               *
c     *************************************************
c
c  This routine computes the seven point derivative of f.
c  where f and df are stored on a logarithmic grid. The
c  dimensions of f and df are, f(1:ng), and df(1:ng)

      subroutine paw_basis_derivative_ngrid(ng,log_amesh,r,f,df)
      implicit none
      integer           ng
      double precision  log_amesh
      double precision  r(ng)
      double precision  f(ng)
      double precision df(ng)

      double precision one_over_60
      parameter (one_over_60 = 1.0d0/60.0d0)

      integer i,n1,n2,m1,m2
      double precision aa

      aa = one_over_60/log_amesh
      n1 = 1
      n2 = ng
      m1 = n1
      m2 = n2
     

      if (n1.le.3) then
        if ((n1.eq.1).and.(n1.ge.m1).and.(n1.le.m2)) then
          df(1) = aa*(-147.0d0*f(1)
     >               + 360.0d0*f(2)
     >               - 450.0d0*f(3)
     >               + 400.0d0*f(4)
     >               - 225.0d0*f(5)
     >               +  72.0d0*f(6)
     >               -  10.0d0*f(7))/r(1)
          n1 = n1+1
        end if 
        if ((n1.eq.2).and.(n1.ge.m1).and.(n1.le.m2)) then
          df(2) = aa*( -10.0d0*f(1)
     >               -  77.0d0*f(2)
     >               + 150.0d0*f(3)
     >               - 100.0d0*f(4)
     >               +  50.0d0*f(5)
     >               -  15.0d0*f(6)
     >               +   2.0d0*f(7))/r(2)
          n1 = n1+1
        end if 
        if ((n1.eq.3.and.(n1.ge.m1).and.(n1.le.m2))) then
          df(3) = aa*(  +2.0d0*f(1)
     >               -  24.0d0*f(2)
     >               -  35.0d0*f(3)
     >               +  80.0d0*f(4)
     >               -  30.0d0*f(5)
     >               +   8.0d0*f(6)
     >               -   1.0d0*f(7))/r(3)
          n1 = n1+1
        end if 
      end if

      if (n2.ge.(ng-2)) then
        if ((n2.eq.ng).and.(n2.ge.m1).and.(n2.le.m2)) then
          df(ng) = aa*( +147.0d0*f(ng)
     >                - 360.0d0*f(ng-1)
     >                + 450.0d0*f(ng-2)
     >                - 400.0d0*f(ng-3)
     >                + 225.0d0*f(ng-4)
     >                -  72.0d0*f(ng-5)
     >                +  10.0d0*f(ng-6))/r(ng)
          n2 = n2-1
        end if
        if ((n2.eq.(ng-1).and.(n2.ge.m1).and.(n2.le.m2))) then
          df(ng-1) = aa*( +10.0d0*f(ng)
     >                  +  77.0d0*f(ng-1)
     >                  - 150.0d0*f(ng-2)
     >                  + 100.0d0*f(ng-3)
     >                  -  50.0d0*f(ng-4)
     >                  +  15.0d0*f(ng-5)
     >                  -   2.0d0*f(ng-6))/r(ng-1)
          n2 = n2-1
        end if
        if ((n2.eq.(ng-2).and.(n2.ge.m1).and.(n2.le.m2))) then
          df(ng-2) = aa*(  -2.0d0*f(ng)
     >                  +  24.0d0*f(ng-1)
     >                  +  35.0d0*f(ng-2)
     >                  -  80.0d0*f(ng-3)
     >                  +  30.0d0*f(ng-4)
     >                  -   8.0d0*f(ng-5)
     >                  +   1.0d0*f(ng-6))/r(ng-2)
          n2 = n2-1
        end if
      end if

      do i=n1,n2
        df(i) = aa*(  -1.0d0*f(i-3)
     >             +   9.0d0*f(i-2)
     >             -  45.0d0*f(i-1)
     >             +  45.0d0*f(i+1)
     >             -   9.0d0*f(i+2)
     >             +   1.0d0*f(i+3))/r(i)
      end do

      return
      end
      


