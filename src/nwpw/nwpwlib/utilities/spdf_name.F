*
* $Id$
*

      character function spdf_name(l)
      implicit none
      integer l
      character name
      name = '?'
      if (l.eq.0) name = 's'
      if (l.eq.1) name = 'p'
      if (l.eq.2) name = 'd'
      if (l.eq.3) name = 'f'
      if (l.eq.4) name = 'g'
      if (l.eq.5) name = 'h'
      if (l.eq.6) name = 'i'
      if (l.eq.7) name = 'j'
      if (l.eq.8) name = 'k'
      if (l.eq.9) name = 'l'
      spdf_name = name
      return
      end
