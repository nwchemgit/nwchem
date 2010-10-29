      integer function argos_prepare_atnum(element)
c
c $Id$
c
      implicit none
      character*2 element
c
      character*2 aname(0:105),bname(0:105)
c
      data aname / '  ',
     + ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     + 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     + 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     + 'Ga','Ge','As','Se','Br',' R','Rb','Sr',' Y','Zr',
     + 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     + 'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     + 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     + 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     + 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     + 'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     + 'Md','No','Lr','Rf','Ha' /
      data bname / '  ',
     + ' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
     + 'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA',
     + 'SC','TI',' V','CR','MN','FE','CO','NI','CU','ZN',
     + 'GA','GE','AS','SE','BR',' R','RB','SR',' Y','ZR',
     + 'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
     + 'SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',
     + 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     + 'LU','HF','TA',' W','RE','OS','IR','PT','AU','HG',
     + 'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     + 'PA',' U','Np','PU','AM','CM','BK','CF','ES','FM',
     + 'MD','No','LR','RF','HA' /
c
      integer i
c
      do 1 i=1,105
      if(element.eq.aname(i).or.element.eq.bname(i)) then
      argos_prepare_atnum=i
      return
      endif
    1 continue
c
      argos_prepare_atnum=0
      if(element.eq.'1H') argos_prepare_atnum=1
      if(element.eq.'2H') argos_prepare_atnum=1
      if(element.eq.'3H') argos_prepare_atnum=1
      if(element.eq.'4H') argos_prepare_atnum=1
c
      return
      end
