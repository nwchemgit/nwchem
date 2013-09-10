      CHARACTER*4 AO_CENT_FRAG, AO_TYPE_FRAG
      CHARACTER*72 TITMOL_FRAG
      CHARACTER*24 TITLE_FRAG
      CHARACTER*3 CFRAG
      COMMON/FRAGMOL/AO_CENT_FRAG(MXPORB,0:MXPFRAG), 
     &               AO_TYPE_FRAG(MXPORB,0:MXPFRAG),
     &               TITMOL_FRAG(2,0:MXPFRAG),
     &               TITLE_FRAG(0:MXPFRAG),
     &               NBAS_FRAG(MXPOBS,0:MXPFRAG),
     &               NFRAG_IN, IFRAG_TP(0:MXPFRAG),
     &               IB_BAS_SM_FRAG(MXPOBS,MXPFRAG),
     &               IB_BAS_GS_SM_FRAG(MXPNGAS,MXPOBS,MXPFRAG),
     &               N_GS_SM_BAS_FRAG(0:MXPNGAS,MXPOBS,MXPFRAG),
     &               NFRAG_TP, NFRAG_MOL, IFRAG_MOL(MXPFRAG),
     &               CFRAG(MXPFRAG),  !End of old
     &               NSMOB_FRAG(MXPFRAG), 
     &               NEQVGRP_FRAG,LEQVGRP_FRAG(MXPFRAG),
     &               IEQVGRP_FRAG(MXPFRAG,MXPFRAG),
     &               NSYM_FRAG(MXPFRAG)

     


