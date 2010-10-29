C $Id$
      Character*2 Function drdy_atsym(iat)
C
      Implicit None
      Integer iat
C
      Integer i
      Character*2 asymb(103)
C Array of atomic symbols
      Data (asymb(i),i=1,103)  /
     *   'h ','he','li','be','b ','c ','n ','o ','f ','ne', 
     *   'na','mg','al','si','p ','s ','cl','ar',
     *   'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     *        'zn','ga','ge','as','se','br','kr',
     *   'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     *        'cd','in','sn','sb','te','i ','xe',
     *   'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb',
     *        'dy','ho','er','tm','yb','lu','hf','ta','w ','re',
     *        'os','ir','pt','au','hg','tl','pb','bi','po','at','rn',
     *   'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk',
     *        'cf','es','fm','md','no','lr'/
C
      drdy_atsym = asymb(iat)
      return
      end
