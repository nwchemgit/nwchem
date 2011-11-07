#$Id$

OBJ_OPTIMIZE = beta_1_all.o beta_2_all.o beta_3_all.o beta_4_all.o beta_all_offsets.o \

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = beta_1_all.F beta_2_all.F beta_3_all.F beta_4_all.F \
 beta_1_3_3_1_1.F beta_1_4_1_2_1.F beta_1_3_1_4_1.F beta_3_7_5_1_1.F beta_3_1_3_1.F \
 beta_3_2_1_2.F beta_1_4_2_2.F beta_1_3_4_1.F beta_1_4_1_4.F beta_1_3_1_7.F \
 beta_3_7_2_1.F beta_3_2_9_2.F beta_3_3_1.F beta_4_1_2.F beta_3_2_3.F \
 beta_3_4_11_1.F beta_3_4_10_3.F beta_1_4_5.F beta_1_3_7.F beta_1_4_1_1_2_1.F \
 beta_4_3.F beta_3_5.F beta_1_4_2_1_1.F beta_1_4_6_3_1.F beta_3_7_4_3_1.F \
 beta_3_2_2_1.F beta_1_4_3_1.F beta_3_2_1_3.F beta_3_4_10_1_1.F beta_1_4_1_5.F \
 beta_1_3_1_8.F beta_3_6_4_1.F beta_3_4_7_1.F beta_4_2_1.F beta_3_2_4.F \
 beta_1_4_6.F beta_3_4_10_4.F beta_1_3_8.F beta_4_4.F beta_3_6.F \
 beta_3_7_4_5_1_1.F beta_1_1_1_1_1.F beta_1_4_1_3_1.F beta_1_3_1_5_1.F beta_3_6_1_2_1.F \
 beta_3_6_7_1_1.F beta_3_7_5_2_1.F beta_3_7_4_3_2.F beta_3_3_1_1.F beta_3_2_2_2.F \
 beta_3_1_4_1.F beta_3_2_1_4.F beta_1_3_5_1.F beta_1_4_3_2.F beta_3_7_3_1.F \
 beta_3_6_4_2.F beta_3_4_1.F beta_3_4_12_1.F beta_3_2_5.F beta_3_4_10_5.F \
 beta_1_4_6_1_1_1.F beta_3_7.F beta_1_3_4_1_1.F beta_1_4_2_2_1.F beta_1_4_1_3_2.F \
 beta_3_4_4_2_1.F beta_3_7_6_1_1.F beta_3_7_4_4_1.F beta_3_2_3_1.F beta_3_4_10_2_1.F \
 beta_3_2_2_3.F beta_1_4_4_1.F beta_1_4_3_3.F beta_3_6_5_1.F beta_3_4_8_1.F \
 beta_4_3_1.F beta_3_4_2.F beta_3_2_6.F beta_3_4_10_6.F beta_1_1_1_2_1.F \
 beta_1_4_3_1_1.F beta_3_2_1_3_1.F beta_1_4_1_4_1.F beta_1_3_1_6_1.F beta_1_1_1_1.F \
 beta_3_4_11_1_1.F beta_1_3_6_1.F beta_3_2_2_4.F beta_1_4_3_4.F beta_3_7_4_1.F \
 beta_4_3_2.F beta_3_5_1.F beta_3_4_3.F beta_3_4_13_1.F beta_3_2_7.F \
 beta_3_4_10_7.F beta_1_1_2_1_1.F beta_3_7_6_2_1.F beta_3_7_4_5_1.F beta_1_1_1_2.F \
 beta_3_4_1_1.F beta_3_2_4_1.F beta_1_4_5_1.F beta_3_4_10_3_1.F beta_1_3_6_2.F \
 beta_3_7_4_2.F beta_3_6_6_1.F beta_3_4_9_1.F beta_4_4_1.F beta_4_3_3.F \
 beta_3_4_4.F beta_3_4_13_2.F beta_3_2_8.F beta_3_4_10_8.F beta_1_4_6_2_1_1.F \
 beta_1_3_8_1_1_1.F beta_1_2_1_1_1.F beta_1_4_3_2_1.F beta_1_3_5_1_1.F beta_3_2_1_4_1.F \
 beta_1_4_1_5_1.F beta_1_3_1_7_1.F beta_3_7_7_1_1.F beta_1_1_2_1.F beta_3_4_1_2.F \
 beta_1_3_7_1.F beta_3_7_5_1.F beta_3_7_4_3.F beta_1_1_1.F beta_3_6_1.F \
 beta_4_3_4.F beta_3_4_14_1.F beta_3_4_5.F beta_3_2_9.F beta_1_4_1_3_2_1.F \
 beta_1_1_2_2_1.F beta_1_4_4_1_1.F beta_3_2_2_3_1.F beta_3_7_6_3_1.F beta_1_2_1_1.F \
 beta_1_1_2_2.F beta_3_4_2_1.F beta_3_2_5_1.F beta_3_4_12_1_1.F beta_3_4_1_3.F \
 beta_1_4_6_1.F beta_3_4_10_4_1.F beta_1_3_7_2.F beta_3_7_5_2.F beta_3_6_7_1.F \
 beta_3_7_4_4.F beta_1_1_2.F beta_3_6_2.F beta_4_3_5.F beta_3_4_14_2.F \
 beta_3_4_6.F beta_3_7_4_1_2_1.F beta_1_1_3_1_1.F beta_1_4_3_3_1.F beta_1_3_1_8_1.F \
 beta_3_2_9_1_1.F beta_2_1_1_1.F beta_1_1_3_1.F beta_4_3_2_1.F beta_3_5_1_1.F \
 beta_1_4_6_2.F beta_1_3_8_1.F beta_3_7_6_1.F beta_3_7_4_5.F beta_1_2_1.F \
 beta_1_1_3.F beta_3_7_1.F beta_3_6_3.F beta_3_4_7.F beta_1_1.F \
 beta_1_4_6_3_1_1.F beta_1_2_2_1_1.F beta_1_3_6_1_1.F beta_3_2_2_4_1.F beta_3_4_10_8_1_1.F \
 beta_3_7_8_1_1.F beta_3_7_6_4_1.F beta_2_1_1_2.F beta_1_2_2_1.F beta_1_1_3_2.F \
 beta_4_3_2_2.F beta_3_4_3_1.F beta_3_2_6_1.F beta_3_4_10_5_1.F beta_1_4_6_3.F \
 beta_3_7_6_2.F beta_2_1_1.F beta_1_2_2.F beta_3_7_2.F beta_3_6_4.F \
 beta_3_4_8.F beta_1_2.F beta_1_3_1_1_1.F beta_1_1_3_2_1.F beta_1_4_5_1_1.F \
 beta_1_4_3_4_1.F beta_3_2_9_2_1.F beta_2_1_2_1.F beta_1_3_1_1.F beta_4_3_3_1.F \
 beta_3_4_3_2.F beta_3_4_13_1_1.F beta_3_7_7_1.F beta_3_7_6_3.F beta_1_3_1.F \
 beta_2_1_2.F beta_1_2_3.F beta_3_7_3.F beta_3_6_5.F beta_3_4_9.F \
 beta_2_1.F beta_1_3.F beta_1_4_1_5_1_1.F beta_1_3_6_2_1.F beta_3_6_4_2_1.F \
 beta_2_1_2_2.F beta_1_3_1_2.F beta_1_2_3_1.F beta_3_6_1_1.F beta_3_4_4_1.F \
 beta_4_3_3_2.F beta_3_2_7_1.F beta_3_4_10_6_1.F beta_3_7_6_4.F beta_2_1_3.F \
 beta_1_3_2.F beta_3_7_4.F beta_3_6_6.F beta_1_4.F beta_3_4_10.F \
 beta_1_3_1_2_1.F beta_1_2_3_1_1.F beta_3_4_1_3_1.F beta_1_3_7_1_1.F beta_3_1_1_1.F \
 beta_1_3_2_1.F beta_2_1_3_1.F beta_1_3_1_3.F beta_4_3_4_1.F beta_3_6_1_2.F \
 beta_3_4_13_2_1.F beta_3_4_4_2.F beta_3_7_8_1.F beta_3_1_1.F beta_1_4_1.F \
 beta_1_3_3.F beta_2_1_4.F beta_3_7_5.F beta_3_6_7.F beta_3_1.F \
 beta_3_4_11.F beta_1_3_2_1_1.F beta_1_4_6_1_1.F beta_3_7_4_1_1.F beta_1_4_1_1.F \
 beta_3_1_1_2.F beta_1_3_1_4.F beta_3_6_2_1.F beta_3_4_14_1_1.F beta_3_4_5_1.F \
 beta_3_2_8_1.F beta_3_4_10_7_1.F beta_3_1_2.F beta_1_4_2.F beta_1_3_4.F \
 beta_3_7_6.F beta_3_2.F beta_3_4_12.F beta_3_7_4_3_2_1.F beta_1_4_1_1_1.F \
 beta_1_3_1_3_1.F beta_1_3_7_2_1.F beta_3_7_4_1_2.F beta_3_1_2_1.F beta_1_4_1_2.F \
 beta_1_3_3_1.F beta_2_1_4_1.F beta_1_3_1_5.F beta_3_7_1_1.F beta_4_3_5_1.F \
 beta_3_4_5_2.F beta_3_2_1.F beta_3_1_3.F beta_1_4_3.F beta_3_4_10_1.F \
 beta_1_3_5.F beta_3_7_7.F beta_4_1.F beta_3_3.F beta_3_4_13.F \
 beta_1_3_1_8_1_1.F beta_1_4_1_1_2.F beta_1_4_6_2_1.F beta_1_3_8_1_1.F beta_3_7_4_2_1.F \
 beta_3_2_1_1.F beta_1_4_2_1.F beta_3_1_2_2.F beta_1_4_1_3.F beta_1_3_1_6.F \
 beta_3_6_3_1.F beta_3_4_14_2_1.F beta_3_4_6_1.F beta_3_2_9_1.F beta_3_4_10_8_1.F \
 beta_4_1_1.F beta_3_2_2.F beta_3_1_4.F beta_1_4_4.F beta_3_4_10_2.F \
 beta_1_3_6.F beta_3_7_8.F beta_4_2.F beta_3_4.F beta_3_4_14.F \

LIB_DEFINES = -DDEBUG_PRINT
#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

HEADERS = 

include ../../config/makefile.h
include ../../config/makelib.h

