HEADERS   = ecp_nwc.fh
OBJ_GRAD = \
      ecp_gradient.o \
      ecp_grad_a.o \
      ecp_grad_b.o

OBJ_HESS = \
      ecp_hessian.o

OBJ_ENERGY = \
      ecp_3j_prod.o  \
      ecp_angint.o   \
      ecp_angrad.o   \
      ecp_bess_asy.o \
      ecp_bess_pow.o \
      ecp_bessel.o   \
      ecp_c2slower.o \
      ecp_c2snorm.o  \
      ecp_c2srecur.o \
      ecp_cart_xpd.o \
      ecp_contract.o \
      ecp_cstrans.o  \
      ecp_down_m.o   \
      ecp_gaussq.o   \
      ecp_init_c2s.o \
      ecp_init_con.o \
      ecp_int0.o     \
      ecp_int1.o     \
      ecp_int2.o     \
      ecp_integral.o \
      ecp_matpi.o    \
      ecp_matpr.o    \
      ecp_radint0.o  \
      ecp_radint1.o  \
      ecp_radint2.o  \
      ecp_s2ccoef.o  \
      ecp_s2cnorm.o  \
      ecp_sph_tens.o \
      ecp_t1_Q0asy.o \
      ecp_t1_Q0int.o \
      ecp_t1_Q0pow.o \
      ecp_t1_Q1asy.o \
      ecp_t1_Q1int.o \
      ecp_t1_Q1pow.o \
      ecp_t1_Q2int.o \
      ecp_t2_ghq.o   \
      ecp_t2_init0.o \
      ecp_t2_init1.o \
      ecp_t2_init2.o \
      ecp_t2_p0pow.o \
      ecp_t2_p1pow.o \
      ecp_up_k.o     \
      ecp_local.o    \
      int_ecp.o      

OBJ =  

OBJ_OPTIMIZE = $(OBJ_ENERGY) $(OBJ_GRAD)

HEADERS = ecp_nwc.fh

USES_BLAS = ecp_cart_xpd.f \
      ecp_cstrans.f  \
      ecp_int0.F     \
      ecp_int1.F     \
      ecp_int2.F     \
      ecp_integral.F \
      ecp_local.F    \
      ecp_radint0.F  \
      ecp_radint1.F  \
      ecp_radint2.F  \
      ecp_s2ccoef.f  \
      int_ecp.F      

LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h


