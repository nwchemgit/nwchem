*
* $Id: paw_basis.F,v 1.27 2004-10-14 00:28:08 bylaska Exp $
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
#include "errquit.fh"

      
#include "paw_basis.fh"
#include "mafdecls.fh"
  
*     !*** local variables ***
      integer it, ir, nt
      integer l_nb,l_ng
      integer i_nb,i_ng
      integer i_lo,i_hi,ii,mult_l
      integer  max_mult_l
      integer max_ic
      logical ok
      double precision log_amesh,amesh
      character*255 fname

c      !*** external functions ***
      integer ion_nkatm
      external ion_nkatm
      
*     *** check for existence of _basis files ***
      call paw_basis_exist_check()

      nt = ion_nkatm()

*     **** allocating memory from stack ****
      ok = MA_push_get(MT_INT,nt, "nb_tmp_array", l_nb,i_nb)
      ok = ok.and.
     >     MA_push_get(MT_INT,nt, "ng_tmp_array", l_ng,i_ng)
      if (.not.ok) 
     > call errquit('init_paw_basis: out of stack memory',0, MA_ERR)
  
      do it = 1, nt
       call ion_atom_plus_suffix(it,'_basis',fname)
       call util_file_name_resolve(fname, .false.)
       call read_paw_basis_dims(int_mb(i_nb-1 + it),
     >                          int_mb(i_ng-1 + it),
     >                          fname)
      end do

      call alloc_paw_basis_data(nt,int_mb(i_nb),int_mb(i_ng))


      do it = 1,nt 
         call ion_atom_plus_suffix(it,'_basis',fname)
         call util_file_name_resolve(fname, .false.)
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
         
c       !***** determine max orbital moment for this kind of atom ****
         i_lo =  paw_basis_i_orb_l(it)
         i_hi =  paw_basis_i_orb_l(it)  
     >        + int_mb(paw_basis_i_nbasis(it)) - 1
         mult_l = int_mb(i_lo)
         do ii=i_lo+1,i_hi
           if (int_mb(ii).gt.mult_l) mult_l = int_mb(ii)
         end do
         int_mb(paw_basis_i_mult_l(it)) = 2*mult_l
           
      end do


c     !***** determine largest max orbital moment ****
      max_mult_l = 0
      do it=1,nt
        if (int_mb(paw_basis_i_mult_l(it)).gt.max_mult_l) 
     >    max_mult_l = int_mb(paw_basis_i_mult_l(it))
      end do
      call set_max_mult_l(max_mult_l)

c     !***** determine largest sphere cutoff ****
      max_ic = 0
      do it=1,nt
        if (int_mb(paw_basis_i_i_cut(it)).gt.max_ic) 
     >    max_ic = int_mb(paw_basis_i_i_cut(it))
      end do
      call set_max_i_cut(max_ic)


c     !*** set core charges ***
      do it = 1,nt 
        call calc_paw_basis_core_charge(
     (        it,dbl_mb(paw_basis_i_core_charge(it)))
      end do


*     **** deallocating memory from stack ****
      ok =        MA_pop_stack(l_ng)
      ok = ok.and.MA_pop_stack(l_nb)
      if (.not.ok) 
     > call errquit('init_paw_basis:error popping stack',1, MA_ERR)

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
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer          ng
      integer          nb
      character*(*)     filename

*     !*** local variables ***
      integer     fp
      double precision dummy_r1
      double precision dummy_rmax
      integer MASTER,taskid
      parameter (MASTER=0)
      integer msglen

*     !*** external functions ****
      logical  assert_file
      external assert_file

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then

      if( .not. assert_file(filename)) then
         write(luout,*)  ' cannot find ',filename
        call errquit("paswbasis: cannot find file",1, DISK_ERR)
      end if

      call open_file(filename,"r",fp)
        READ(fp,*) dummy_r1,dummy_rmax,ng,nb
      close(fp)
      endif
      msglen = 1
      call BRDCST(9+MSGINT,ng,             mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,nb,           mitob(msglen),MASTER)
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
#include "errquit.fh"
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

#include "stdio.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"


*     !*** local variables
      integer MASTER,taskid
      parameter (MASTER=0)
      integer i,j,msglen
      integer fp
      !character*255 filename

*     !*** external functions ***
      logical  assert_file
      external assert_file

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then

      if( .not. assert_file(aname)) then
        write(luout,*)  ' cannot find ',aname
        call errquit('paswbasis: cannot find file',2, DISK_ERR)
      end if

      call open_file(aname,"r",fp)
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
      end if

*     **** send header data to all processors ****
      msglen = 1
      call BRDCST(9+MSGDBL,r1,             mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,rmax,           mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,core_kin_energy,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,sigma,          mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,ion_charge,     mdtob(msglen),MASTER)

      call BRDCST(9+MSGINT,nbasis,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,ngrid, mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,i_cut, mitob(msglen),MASTER)


      msglen = nbasis
      call BRDCST(9+MSGINT,prin_n,   mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,prin_n_ps,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,orb_l,    mitob(msglen),MASTER)

      call BRDCST(9+MSGDBL,eig,      mdtob(msglen),MASTER)


      msglen = ngrid
      call BRDCST(9+MSGDBL,core_density,   mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,ps_core_density,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,v_ps,           mdtob(msglen),MASTER)


      msglen = ngrid*nbasis
      call BRDCST(9+MSGDBL,phi_ae,      mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,phi_ps,      mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,phi_ae_prime,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,phi_ps_prime,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,prj_ps,      mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,prj_ps0,     mdtob(msglen),MASTER)


      return
      end 




!**************************************************
!
!       Name: paw_basis_exist_check
!
!       Purpose:
!
!       Created:        2/12/2002
!**************************************************
      subroutine paw_basis_exist_check()
      implicit none
#include "rtdb.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical oldmode
      logical found
      integer ia,lmax,locp
      character*2 atom
      character*255 fname,libfile
      real*8        sum,rlocal

*     **** external functions ****
      character*2 ion_atom
      integer     ion_nkatm
      external    ion_atom
      external    ion_nkatm

      call Parallel_taskid(taskid)
      oldmode =  rtdb_parallel(.false.)
      if (taskid.eq.MASTER) then
      do ia = 1, ion_nkatm()
       call ion_atom_plus_suffix(ia,'_basis',fname)
       call util_file_name_resolve(fname, .false.)

       inquire(file=fname,exist=found)
       if (.not.found) then

*         **** generate one-dimensional pseudopotential file ****
          atom = ion_atom(ia)
          call psp_generator_auto(2,atom,lmax,locp,rlocal)

c*         *** find location of library file ***
c          atom = ion_atom(ia)
c          call nwpw_libfile(2,libfile,atom)
c
c*         *** copy library file to fname ***
c          call util_file_copy(libfile,fname)
      
       end if

      end do
      end if
      call ga_sync()
      oldmode =  rtdb_parallel(oldmode)
      !*** sync - hack ***
      !sum = 1.0d0
      !call D3dB_SumAll(sum)

      return
      end


