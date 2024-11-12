      MODULE nwc_sym

      IMPLICIT NONE
      SAVE

#include "nwc_const.fh" !nw_max_shells, nw_max_nbf
#include "geomP.fh"     !max_cent, coords

      INTEGER, PARAMETER :: mxsym = 8 	!max number of symmetry groups
      INTEGER, PARAMETER :: mxshel=nw_max_shells!max number of shells
      INTEGER, PARAMETER :: mxbf=nw_max_nbf	!max number of basis functions
      INTEGER, PARAMETER :: mxqn=7
      INTEGER, PARAMETER :: mxaqn=mxqn*(mxqn+1)/2
      INTEGER, PARAMETER :: mxqnm=21
      INTEGER, PARAMETER :: mxaqnm=mxqnm*(mxqnm+1)/2

      INTEGER :: isymax(3,2)		! irreps of (c1) x,y,z, (c2) Rx,Ry,Rz
      INTEGER :: sh_op_map(8,mxshel), bf_op_map(8,mxbf)
      INTEGER :: sh_n_uq_op(mxshel), bf_n_uq_op(mxbf)
      INTEGER :: sh_uq_op(8,mxshel), bf_uq_op(8,mxbf)
      INTEGER :: sh_n_uq, bf_n_uq
      INTEGER :: sh_uq(mxshel), bf_uq(mxbf)
      INTEGER :: sh_nat(2,mxshel), bf_nat(2,mxbf), sh_uq_bf(4,mxshel)
      INTEGER :: bf_per_ir(0:7)
      INTEGER :: bf_per_ir_cum(0:7)
      INTEGER :: bf_so_ir(0:7,mxbf)
      INTEGER :: so_uq_to_lab(2,0:7,mxbf)
      INTEGER :: so_lab1(3,mxbf)
      INTEGER :: so_lab2(mxbf)
      INTEGER :: isymao(mxqn,mxaqn)
      INTEGER :: stabilizer(max_cent) !stabilizer for each center
      INTEGER :: ngen ! number of generating operators
      INTEGER :: parbit(0:7) = (/1,-1,-1,1,-1,1,1,-1/) !parity of 3-bit integer 000-111
      INTEGER :: itran(mxbf,mxbf)
      INTEGER :: irrepmap(0:mxsym-1)=(/0,1,2,3,4,5,6,7/) !map from Dalton to NWChem irrep order
      INTEGER :: opmap(0:mxsym-1)=(/0,1,2,3,4,5,6,7/) !map from Dalton to NWChem operator order
      INTEGER :: wt(mxbf,2) !How many AOs make up each SO
      INTEGER :: nbpsy(1:mxsym) ! number of basis functions per symmetry
      INTEGER :: nop !number of symmetry operators
      INTEGER :: nir !number of irreps
  
      LOGICAL :: oprint 

      DOUBLE PRECISION :: bf_phase(8,mxbf) 
      DOUBLE PRECISION :: char_tab(1:8,0:7)

      CHARACTER*8 :: grp_name 		! name of symmetry group
      CHARACTER*8 :: zir(mxsym)		! irrep names



      END MODULE nwc_sym
