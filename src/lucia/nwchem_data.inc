      character*8 luciasymname
      character*255 luciascfmo
      integer luciageom, luciabasis
      integer ipnt_nw(mxpirr)
      integer ioff_nw(mxpirr),icount_nw(mxpirr)
      integer kirs_nw(mxporb),reord_nw(mxporb)
      integer nwchem_id(mxporb,mxpirr)
      double precision nwchem_sc(mxporb,mxpirr)
      integer luciaao,luciamo(2)
      integer igalist(50),ignlstnr
      COMMON/NWCHEM_DATAI/luciageom,luciabasis,
     &                    ipnt_nw,ioff_nw,icount_nw,kirs_nw,
     &                    reord_nw,luciaao,luciamo,nwchem_id,
     &                    igalist,igalstnr
      COMMON/NWCHEM_DATAD/nwchem_sc
      COMMON/NWCHEM_DATAC/luciascfmo,luciasymname
