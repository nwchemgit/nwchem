*
* $Id: paw_basis.F,v 1.12 2003-02-04 01:00:29 marat Exp $
*

!**************************************************
!
!	Name: init_paw_basis
!
!	Purpose 
!
!	Created:	7/11/2002
!**************************************************
      subroutine init_paw_basis()
      implicit none

      include "paw_basis.fh"
#include "mafdecls.fh"
c      include "myma.fh"
  
*     !*** local variables ***
      integer it, ir, nt
      integer  nb,ng
      integer l_nb,l_ng
      integer i_nb,i_ng
      integer i_lo,i_hi
      integer  max_mult_l
      logical ok
      double precision log_amesh,amesh
      character*20 fname

c      !*** external functions ***
      integer ion_nkatm
      external ion_nkatm

      nt = ion_nkatm()
      ok = MA_push_get(MT_INT,nt, "nb_tmp_array", l_nb,i_nb)
      ok = MA_push_get(MT_INT,nt, "ng_tmp_array", l_ng,i_ng)
  
      do it = 1, nt
       call ion_atom_plus_suffix(it,'_basis',fname)
       call read_paw_basis_dims(int_mb(i_nb-1 + it),
     >                          int_mb(i_ng-1 + it),
     >                          fname)
      end do

      call alloc_paw_basis_data(nt,int_mb(i_nb),int_mb(i_ng))


      do it = 1,nt 
         call ion_atom_plus_suffix(it,'_basis',fname)
         call read_paw_basis  ( dbl_mb(paw_basis_i_r1(it)),
     >                         dbl_mb(paw_basis_i_rmax(it)),
     >                         int_mb(paw_basis_i_ngrid(it)),
     >                         int_mb(paw_basis_i_nbasis(it)),
     >                         int_mb(paw_basis_i_i_cut(it)),
     >                         int_mb(paw_basis_i_pr_n(it)),
     >                         int_mb(paw_basis_i_pr_n_ps(it)),
     >                         int_mb(paw_basis_i_orb_l(it)),
     >                         dbl_mb(paw_basis_i_sigma(it)),
     >                         dbl_mb(paw_basis_i_ion_charge(it)),
     >                         dbl_mb(paw_basis_i_core_kin_energy(it)),
     >                         dbl_mb(paw_basis_i_eig(it)),
     >                         dbl_mb(paw_basis_i_phi_ae(it)),
     >                         dbl_mb(paw_basis_i_phi_ae_prime(it)),
     >                         dbl_mb(paw_basis_i_phi_ps(it)),
     >                         dbl_mb(paw_basis_i_phi_ps_prime(it)),
     >                         dbl_mb(paw_basis_i_prj_ps(it)),
     >                         dbl_mb(paw_basis_i_prj_ps0(it)),
     >                         dbl_mb(paw_basis_i_core_density(it)),
     >                         dbl_mb(paw_basis_i_ps_core_density(it)),
     >                         dbl_mb(paw_basis_i_v_ps(it)),
     >                         fname)

  
         log_amesh = LOG(dbl_mb(paw_basis_i_rmax(it))
     >                  /dbl_mb(paw_basis_i_r1(it)))
     >                 /(int_mb(paw_basis_i_ngrid(it))-1)
         amesh     = EXP(log_amesh)
    
         dbl_mb(paw_basis_i_log_amesh(it)) = log_amesh
         dbl_mb(paw_basis_i_rgrid(it)) = dbl_mb(paw_basis_i_r1(it))
  
         do ir=1,int_mb(paw_basis_i_ngrid(it))-1
          dbl_mb(paw_basis_i_rgrid(it)+ir)  = 
     >           amesh*dbl_mb(paw_basis_i_rgrid(it)+ir-1)
         end do
         
         i_lo =  paw_basis_i_orb_l(it)
	 i_hi =  paw_basis_i_orb_l(it)  
     >        + int_mb(paw_basis_i_nbasis(it)) - 1
	 int_mb(paw_basis_i_mult_l(it)) = 
     >        2*dabsmax(i_hi-i_lo+1,int_mb(i_lo))
           
      end do

      ok = MA_chop_stack(l_nb)

c     !*****determine max orbital moment     
      i_lo =  paw_basis_i_orb_l(1)
      i_hi =  paw_basis_i_orb_l(1) +  paw_basis_tot_nbasis() - 1
      max_mult_l = dabsmax(i_hi-i_lo+1,int_mb(i_lo))


c     !*****enter max orbital moment into paw_basis_data structure
      call set_max_mult_l(max_mult_l)

c     !*** set core charges ***
            
      do it = 1,nt 
        call calc_paw_basis_core_charge(
     (        it,dbl_mb(paw_basis_i_core_charge(it)))
      end do

      return
      end 

!**************************************************
!
!	Name: read_paw_basis_dims
!
!	Purpose 
!
!	Created:	7/11/2002
!**************************************************
      subroutine read_paw_basis_dims(nb,ng,filename)
      implicit none
#include "stdio.fh"
      integer          ng
      integer          nb
      character(*)     filename

*     !*** local variables ***
      integer     fp
      double precision dummy_r1
      double precision dummy_rmax

*     !*** external functions ****
      logical  assert_file
      external assert_file

      if( .not. assert_file(filename)) then
         write(luout,*)  ' cannot find ',filename
        call errquit("paswbasis: cannot find file",1)
      end if

      call open_file(filename,"r",fp)
        READ(fp,*) dummy_r1,dummy_rmax,ng,nb
      close(fp)
      return
      end


!**************************************************
!
!	Name: read_paw_basis
!
!	Purpose 
!
!	Created:	7/11/2002
!************************************************** 
      subroutine read_paw_basis(r1,
     >                          rmax,
     >                          ngrid,
     >                          nbasis,
     >                          i_cut,
     >                          prin_n,
     >                          prin_n_ps,
     >                          orb_l,
     >                          sigma,
     >                          ion_charge,
     >                          core_kin_energy,
     >                          eig,
     >                          phi_ae,
     >                          phi_ae_prime,
     >                          phi_ps,
     >                          phi_ps_prime,
     >                          prj_ps,
     >                          prj_ps0,
     >                          core_density,
     >                          ps_core_density,
     >                          v_ps,
     >                          aname)
      implicit none
      double precision r1
      double precision rmax
      integer          ngrid
      integer          nbasis
      integer          i_cut
      integer          prin_n(nbasis)
      integer          prin_n_ps(nbasis)
      integer          orb_l(nbasis)
      double precision sigma
      double precision ion_charge
      double precision core_kin_energy
      double precision eig(nbasis)
      double precision phi_ae(nbasis*ngrid)
      double precision phi_ae_prime(nbasis*ngrid)
      double precision phi_ps(nbasis*ngrid)
      double precision phi_ps_prime(nbasis*ngrid)
      double precision prj_ps(nbasis*ngrid)
      double precision prj_ps0(nbasis*ngrid)
      double precision core_density(ngrid)
      double precision ps_core_density(ngrid)
      double precision v_ps(ngrid)
      character*(*)    aname


*     !*** local variables
      integer i,j
      integer fp
      character*40 filename

*     !*** external functions ***
      logical  assert_file
      external assert_file

      filename = aname
      if( .not. assert_file(filename)) then
        call errquit("cannot find "//filename,0)
      end if

      call open_file(filename,"r",fp)
        READ(fp,*) r1,rmax,ngrid,nbasis
        READ(fp,*) i_cut
        READ(fp,*) core_kin_energy
        READ(fp,*) (prin_n(i),eig(i),
     >              prin_n_ps(i),orb_l(i),i=1,nbasis)
        READ(fp,*) ((phi_ae((i-1)*ngrid+j), j=1,ngrid),
     >                                      i=1,nbasis)
        READ(fp,*) ((phi_ae_prime((i-1)*ngrid+j), j=1,ngrid),
     >                                            i=1,nbasis)
        READ(fp,*) ((phi_ps((i-1)*ngrid+j), j=1,ngrid),
     >                                      i=1,nbasis)
        READ(fp,*) ((phi_ps_prime((i-1)*ngrid+j), j=1,ngrid),
     >                                            i=1,nbasis)
        READ(fp,*) ((prj_ps((i-1)*ngrid+j), j=1,ngrid), 
     >                                      i=1,nbasis)
        READ(fp,*) (core_density(i),   i=1,ngrid)
        READ(fp,*) (ps_core_density(i),i=1,ngrid)
        READ(fp,*) (v_ps(i),           i=1,ngrid)
        READ(fp,*) sigma
        READ(fp,*) ion_charge
        READ(fp,*) ((prj_ps0((i-1)*ngrid+j), j=1,ngrid), 
     >                                       i=1,nbasis)
      close(fp)
      return
      end 




