      subroutine argos_diana_task
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      call argos_diana_init()
c
    1 continue
c
      call argos_diana_getcmd()
c
c     read reference coordinates
c     --------------------------
c
      if(cmd.eq.'refer ') then
      call argos_diana_getref()
c
c     read topology file
c     ------------------
c
      elseif(cmd.eq.'topol ') then
      call argos_diana_gettop()
c
c     set io mode
c     -----------
c
      elseif(cmd.eq.'iomode') then
      call argos_diana_setiomode()
c
c     read trajectory file header
c     ---------------------------
c
      elseif(cmd.eq.'file  ') then
      call argos_diana_gethdr()
c
c     select which frames to read
c     ---------------------------
c
      elseif(cmd.eq.'frame ') then
      call argos_diana_gettrj()
c
c     atom selection
c     --------------
c
      elseif(cmd.eq.'select') then
      call argos_diana_getsel()
c
c     save atom selection
c     -------------------
c
      elseif(cmd.eq.'savsel') then
      call argos_diana_savsel()
      elseif(cmd.eq.'savsls') then
      call argos_diana_savssl()
c
c     read atom selection
c     -------------------
c
      elseif(card(1:5).eq.'rdsel') then
      call argos_diana_rdsel()
      elseif(cmd.eq.'rdsels') then
      call argos_diana_rdssl()
c
c     geometric selections
c
      elseif(cmd.eq.'distan') then
      call argos_diana_adddistance()
      elseif(cmd.eq.'angle') then
      call argos_diana_addangle()
      elseif(cmd.eq.'torsio') then
      call argos_diana_addtorsion()
      elseif(cmd.eq.'improp') then
      elseif(cmd.eq.'order') then
      call argos_diana_addorder()
c
c     superimpose all stored frames
c     -----------------------------
c
      elseif(card(1:5).eq.'super') then
      call argos_diana_super()
c
c     calculate average coordinates
c     -----------------------------
c
      elseif(cmd.eq.'xaver ') then
      call argos_diana_average()
c
c     calculate one-dimensional root mean square deviation
c     ----------------------------------------------------
c
      elseif(cmd.eq.'rmsdev') then
      call argos_diana_rmsd()
c
c     build covariance matrix
c     -----------------------
c
      elseif(cmd.eq.'covar ') then
      call argos_diana_covar()
c
c     pca projection
c     --------------
c
      elseif(cmd.eq.'projec') then
      call argos_diana_project()
c
      elseif(cmd.eq.'evalua') then
      call argos_diana_evaluate()
c
      elseif(cmd.eq.'espplt') then
      call argos_diana_espplt()
c
      elseif(cmd.eq.'rama  ') then
      call argos_diana_rama()
      else
      goto 9
      endif
      goto 1
c
    9 continue
      call argos_diana_final()
c
      return
      end
