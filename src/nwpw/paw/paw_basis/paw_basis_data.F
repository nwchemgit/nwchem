!
! $Id: paw_basis_data.F,v 1.2 2002-07-30 21:19:56 marat Exp $
!

!**************************************************
!
!	Name: alloc_paw_basis_data
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      SUBROUTINE  alloc_paw_basis_data(nt,nb,ng)
      implicit none
      integer nt
      integer nb(nt)
      integer ng(nt)

#include "paw_basis_data.fh"
#include "myma.fh"


      logical my_alloc
      external my_alloc

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
        tot_ngrid       = tot_ngrid     + ng(it)
      end do

      ok = MA_set_auto_verify(.TRUE.)

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
      ok = ok .AND. my_alloc(MT_DBL,tot_ngrid,"v_ps",v_ps)
      ok = ok .AND. my_alloc(MT_DBL,tot_ngridbasis,"prj_ps0",prj_ps0)
      ok = ok .AND. my_alloc(MT_DBL,tot_ngrid,"rgrid",rgrid)

      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_nb",i_nb)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_ng",i_ng)
      ok = ok .AND. my_alloc(MT_INT,tot_ntype,"i_ngb",i_ngb)

      if(.not. ok) then
        call errquit("failed to allocate paw_basis data ",0)
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

#include "paw_utils.fh"
#include "paw_basis_data.fh"

      logical ok


        ok = .true.
        ok = ok .and. my_dealloc(i_ngb)
        call chkerr("failed to dealloc_paw_basis_data 1",ok)
        ok = ok .and. my_dealloc(i_ng)
        call chkerr("failed to dealloc_paw_basis_data 2",ok)
        ok = ok .and. my_dealloc(i_nb)
        call chkerr("failed to dealloc_paw_basis_data 3",ok)
        ok = ok .and. my_dealloc(rgrid)
        call chkerr("failed to dealloc_paw_basis_data 4",ok)
        ok = ok .and. my_dealloc(prj_ps0)
        call chkerr("failed to dealloc_paw_basis_data 5",ok)
        ok = ok .and. my_dealloc(v_ps)
        call chkerr("failed to dealloc_paw_basis_data 6",ok)
        ok = ok .and. my_dealloc(ps_core_density)
        call chkerr("failed to dealloc_paw_basis_data 7",ok)
        ok = ok .and. my_dealloc(core_density)
        call chkerr("failed to dealloc_paw_basis_data 8",ok)
        ok = ok .and. my_dealloc(prj_ps)
        call chkerr("failed to dealloc_paw_basis_data 9",ok)

        ok = ok .and. my_dealloc(phi_ps_prime)
        call chkerr("failed to dealloc_paw_basis_data 21",ok)
        ok = ok .and. my_dealloc(phi_ps)
        call chkerr("failed to dealloc_paw_basis_data 22",ok)
        ok = ok .and. my_dealloc(phi_ae_prime)
        call chkerr("failed to dealloc_paw_basis_data 23",ok)
        ok = ok .and. my_dealloc(phi_ae)
        call chkerr("failed to dealloc_paw_basis_data 24",ok)
        ok = ok .and. my_dealloc(eig)
        call chkerr("failed to dealloc_paw_basis_data 25",ok)
        ok = ok .and. my_dealloc(orb_l)
        call chkerr("failed to dealloc_paw_basis_data 26",ok)
        ok = ok .and. my_dealloc(pr_n)
        call chkerr("failed to dealloc_paw_basis_data 27",ok)
        ok = ok .and. my_dealloc(pr_n_ps)
        call chkerr("failed to dealloc_paw_basis_data 28",ok)
        ok = ok .and. my_dealloc(core_kin_energy)
        call chkerr("failed to dealloc_paw_basis_data 29",ok)
        ok = ok .and. my_dealloc(ion_charge) 
        call chkerr("failed to dealloc_paw_basis_data 11",ok)
        ok = ok .and. my_dealloc(sigma)
        call chkerr("failed to dealloc_paw_basis_data 12",ok)
        ok = ok .and. my_dealloc(log_amesh)
        call chkerr("failed to dealloc_paw_basis_data 13",ok)
        ok = ok .and. my_dealloc(rmax)
        call chkerr("failed to dealloc_paw_basis_data 14",ok)
        ok = ok .and. my_dealloc(r1)
        call chkerr("failed to dealloc_paw_basis_data 15",ok)
        ok = ok .and. my_dealloc(mult_l)
        call chkerr("failed to dealloc_paw_basis_data 16",ok)
        ok = ok .and. my_dealloc(ngrid)
        call chkerr("failed to dealloc_paw_basis_data 17",ok)
        ok = ok .and. my_dealloc(i_cut)
        call chkerr("failed to dealloc_paw_basis_data 18",ok)
        ok = ok .and. my_dealloc(nbasis)
        call chkerr("failed to dealloc_paw_basis_data 19",ok)


      END SUBROUTINE

!**************************************************
!
!	Name: get_tot_ntype
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_tot_ntype()
      implicit none

#include "paw_basis_data.fh"

      get_tot_ntype = tot_ntype
      return
      END 

!**************************************************
!
!	Name: get_tot_nbasis
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_tot_nbasis()
      implicit none

#include "paw_basis_data.fh"

      get_tot_nbasis = tot_nbasis
      return
      END 

!**************************************************
!
!	Name: get_tot_ngrid
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_tot_ngrid()
      implicit none

#include "paw_basis_data.fh"

      get_tot_ngrid = tot_ngrid
      return
      END
 
!**************************************************
!
!	Name: get_tot_ngridbasis
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_tot_ngridbasis()
      implicit none

#include "paw_basis_data.fh"
 
      get_tot_ngridbasis = tot_ngridbasis
      return
      END
 
!**************************************************
!
!	Name: get_i_nbasis
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_nbasis(it)
      implicit none
      integer it
   
#include "myma.fh"
#include "paw_basis_data.fh"


      get_i_nbasis = nbasis(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_ngrid
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_ngrid(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_ngrid = ngrid(1) +  it-1
      return
      END
 
!**************************************************
!
!	Name:
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_mult_l(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"


      get_i_mult_l = mult_l(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_r1
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_r1(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_r1 = r1(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_rmax
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_rmax(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"


      get_i_rmax = rmax(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_i_cut
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_i_cut(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_i_cut = i_cut(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_sigma
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_sigma(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_sigma = sigma(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_log_amesh
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_log_amesh(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_log_amesh = log_amesh(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_ion_charge
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_ion_charge(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_ion_charge = ion_charge(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_ion_charge
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      DOUBLE PRECISION FUNCTION  get_ion_charge(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_ion_charge = dbl_mb(ion_charge(1) +  it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_core_kin_energy
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_core_kin_energy(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_core_kin_energy = core_kin_energy(1) +  it-1
      return
      END 

!**************************************************
!
!	Name: get_i_pr_n
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_pr_n(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_pr_n = pr_n(1) + int_mb(i_nb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_pr_n_ps
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_pr_n_ps(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_pr_n_ps = pr_n_ps(1) + int_mb(i_nb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_orb_l
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_orb_l(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_orb_l = orb_l(1) + int_mb(i_nb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_eig
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_eig(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_eig = eig(1) + int_mb(i_nb(1) + it-1)
      return
      END   
  
!**************************************************
!
!	Name: get_i_core_density
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_core_density(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_core_density = core_density(1) + int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_ps_core_density
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_ps_core_density(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_ps_core_density = ps_core_density(1) 
     >                      + int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_v_ps
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_v_ps(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_v_ps = v_ps(1) + int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_rgrid
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_rgrid(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_rgrid = rgrid(1) + int_mb(i_ng(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_phi_ae
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_phi_ae(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_phi_ae = phi_ae(1) + int_mb(i_ngb(1) + it-1)
      return
      END 


!**************************************************
!
!	Name: get_i_phi_ps
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_phi_ps(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_phi_ps = phi_ps(1) + int_mb(i_ngb(1) + it-1)
      return
      END 


!**************************************************
!
!	Name: get_i_phi_ae_prime
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_phi_ae_prime(it)
      implicit none
      integer it
   
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_phi_ae_prime = phi_ae_prime(1) + int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_phi_ps_prime
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_phi_ps_prime(it)
      implicit none
      integer it

#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_phi_ps_prime = phi_ps_prime(1) + int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_prj_ps
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_prj_ps(it)
      implicit none
      integer it
 
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_prj_ps = prj_ps(1) + int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_i_prj_ps0
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      INTEGER FUNCTION  get_i_prj_ps0(it)
      implicit none
      integer it
    
#include "myma.fh"
#include "paw_basis_data.fh"

      get_i_prj_ps0 = prj_ps0(1) + int_mb(i_ngb(1) + it-1)
      return
      END 

!**************************************************
!
!	Name: get_nbasis
!
!	Purpose 
!
!	Created:	7/9/2002
!**************************************************
      SUBROUTINE  get_nbasis(nbasis_out)
      implicit none
      integer nbasis_out(*)

#include "myma.fh"
#include "paw_basis_data.fh"

*     !*** local variables ***
      integer it 

      do it=1,tot_ntype
        nbasis_out(it) = int_mb(nbasis(1) + it - 1)
      end do
    
      return
      END !SUBROUTINE get_nbasis
