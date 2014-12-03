      subroutine water_pseudopotential_input(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
c
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 8)
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'type:',
     >            'number_water:',
     >            'start_atom_number:',
     >            'lj_ion_parameters:',
     >            'lj_water_parameters:',
     >            'sr_basis:',
     >            'polarization_parameters:',
     >            'end'/
c
      logical value
      integer   water_type,water_number
      integer	start_atom_number
      integer   i,ljcount,tmp_ei(2),tmp_si(2),tmp_ionnum(2)
      integer   maxion,geom
      real*8    ei,si
      real*8    ew,sw

*     *** SR Basis ***
      integer num_basis
      integer tmp_c(2),tmp_g(2),tmp_iatom(2)

*     *** LR parameters ****
      real*8 DELTA,Rcoul,q1,q2,q3
      real*8 Rp,Rpe,alpha


*     *** get maxion from rtdb ****
      value = geom_create(geom,'geometry')
      value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      value = value.and.geom_ncent(geom,maxion)
      value = value.and.geom_destroy(geom)
      if (.not.value)
     > call errquit('water_pseudopotential_input:error get maxion',0,
     &       GEOM_ERR)


       
*     **** allocate MA local variables ****
      value = BA_push_get(mt_int,(maxion),
     >                     'tmp_ionnum',tmp_ionnum(2),tmp_ionnum(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(maxion),
     >                     'tmp_ei',tmp_ei(2),tmp_ei(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(maxion),
     >                     'tmp_si',tmp_si(2),tmp_si(1))

      value = value.and.
     >        BA_push_get(mt_int,(100),
     >                     'tmp_iatom',tmp_iatom(2),tmp_iatom(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(100),
     >                     'tmp_c',tmp_c(2),tmp_c(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(100),
     >                     'tmp_g',tmp_g(2),tmp_g(1))
      if (.not.value)
     > call errquit('water_pseudopotential_input:error push stack',0,
     &       MA_ERR)

*     ***** initializations ****
      water_type        = 0
      water_number      = 0
      start_atom_number = 0
      ljcount           = 0
      num_basis         = 0
      ew=(1.55394259d-1/(27.2116d0*23.06d0))
      sw=(3.16555789d0/(0.529177d0))

      DELTA=0.22183756d0  !!*** BLCJ Coulomb Parameters ****
      Rcoul=0.9d0
      q1=0.6d0
      q2=0.6d0
      q3=-1.2d0 
 
      Rp=1.6d0   !!** BLCJ Polarization parameters ****
      Rpe = 3.0d0
      alpha=9.7446d0

      

 10   if (.not. inp_read()) 
     >  call errquit(
     >  'water_pseudopotential_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'water_pseudopotential_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >  call errquit(
     >  'water_pseudopotential_input: unknown directive', 0, INPUT_ERR)


      goto ( 100, 200, 300, 400, 500, 600, 700,
     >      9999) ind
      call errquit(
     >      'cell_optimize_input: unimplemented directive', ind,
     &       INPUT_ERR)

c     
c  water_type
c     
 100  if (.not. inp_a(test))
     >  call errquit(
     >  'water_pseudopotential_input: failed to read water_type', 0,
     &       INPUT_ERR)
      if (test.eq.'unknown') water_type = 0
      goto 10

c
c  water_number
c
 200  if (.not. inp_i(water_number))
     >  call errquit(
     > 'water_pseudopotential_input: failed to read water_number:', 0,
     &       INPUT_ERR)
      goto 10

c
c  start_atom_number
c
 300  if (.not. inp_i(start_atom_number))
     >  call errquit(
     > 'water_pseudopotential_input:failed read start_atom_number:',0,
     &       INPUT_ERR)
      goto 10
     
c
c  LJ parameters
c
 400  value = inp_i(i)
      value = value.and.inp_f(ei)
      value = value.and.inp_f(si)
      if (.not. value)
     >  call errquit(
     > 'water_pseudopotential_input:failed reading lj_parameters:',0,
     &       INPUT_ERR)

      int_mb(tmp_ionnum(1)+ljcount) = i
      dbl_mb(tmp_ei(1)    +ljcount) = ei
      dbl_mb(tmp_si(1)    +ljcount) = si
      ljcount = ljcount+1
      
      goto 10

 500  value = inp_f(ew)
      value = value.and.inp_f(sw)
      if (.not. value)
     >  call errquit(
     > 'water_pseudopotential_input:failed reading lj_parameters:',0,
     &       INPUT_ERR)

      goto 10

c
c  SR Basis
c
 600  value = inp_f(ei)
      value = value.and.inp_f(si)
      value = value.and.inp_i(i)

      if (.not. value)
     >  call errquit(
     > 'water_pseudopotential_input:failed reading sr_basis:',0,
     &       INPUT_ERR)
      
      int_mb(tmp_iatom(1)+num_basis) = i
      dbl_mb(tmp_c(1)    +num_basis) = ei
      dbl_mb(tmp_g(1)    +num_basis) = si
      num_basis = num_basis+1
     
      goto 10

c
c  Polarization parameters
c
 700  value = inp_f(alpha)
      value = value.and.inp_f(Rp)
      value = value.and.inp_f(Rpe)

      if (.not. value)
     >  call errquit(
     > 'water_pseudopotential_input:failed polarization_parameters:',0,
     &       INPUT_ERR)

      goto 10

*     ***** water_pseudopotential parameters to rtdb ****
 9999 continue

      value = rtdb_put(rtdb,'water_pseudopotential:type',
     >                 mt_int,1,water_type)
      value = value.and. 
     >        rtdb_put(rtdb,'water_pseudopotential:number',
     >                 mt_int,1,water_number)
      value = value.and. 
     >        rtdb_put(rtdb,'water_pseudopotential:atom_start',
     >                 mt_int,1,start_atom_number)


*     **** write LJ parameters to rtdb ****
      value = value.and. 
     >        rtdb_put(rtdb,'water_pseudopotential:lj_count',
     >                 mt_int,1,ljcount)
      if (ljcount.gt.0) then
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:lj_ion_number',
     >                   mt_int,ljcount,int_mb(tmp_ionnum(1)))
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:lj_ei',
     >                   mt_dbl,ljcount,dbl_mb(tmp_ei(1)))
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:lj_si',
     >                   mt_dbl,ljcount,dbl_mb(tmp_si(1)))

        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:lj_ew',
     >                   mt_dbl,1,ew)
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:lj_sw',
     >                   mt_dbl,1,sw)
      end if

*     **** write SR basis to rtdb ****
      value = value.and. 
     >        rtdb_put(rtdb,'water_pseudopotential:num_basis',
     >                 mt_int,1,num_basis)
      if (num_basis.gt.0) then
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:sr_ion_number',
     >                   mt_int,num_basis,int_mb(tmp_iatom(1)))
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:sr_c',
     >                   mt_dbl,num_basis,dbl_mb(tmp_c(1)))
        value = value.and. 
     >          rtdb_put(rtdb,'water_pseudopotential:sr_g',
     >                   mt_dbl,num_basis,dbl_mb(tmp_g(1)))
      end if


*     **** write LR parameters to rtdb ****
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:DELTA',
     >                 mt_dbl,1,DELTA)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:Rcoul',
     >                 mt_dbl,1,Rcoul)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:q1',
     >                 mt_dbl,1,q1)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:q2',
     >                 mt_dbl,1,q2)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:q3',
     >                 mt_dbl,1,q3)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:Rp',
     >                 mt_dbl,1,Rp)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:Rpe',
     >                 mt_dbl,1,Rpe)
      value = value.and.
     >        rtdb_put(rtdb,'water_pseudopotential:alpha',
     >                 mt_dbl,1,alpha)


      if (.not.value)
     >  call errquit(
     > 'water_pseudopotential_input:error writing to rtdb',0, RTDB_ERR)

      value = BA_pop_stack(tmp_g(2))
      value = value.and.BA_pop_stack(tmp_c(2))
      value = value.and.BA_pop_stack(tmp_iatom(2))
      value = value.and.BA_pop_stack(tmp_si(2))
      value = value.and.BA_pop_stack(tmp_ei(2))
      value = value.and.BA_pop_stack(tmp_ionnum(2))
      if (.not.value)
     > call errquit('water_pseudopotential_input:error pop stack',0,
     &       MA_ERR)



      return
      end
c $Id$
