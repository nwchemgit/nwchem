##############################################################################
# Generic MHPCC .cshrc - initialization file for the csh command interpreter
#
# See the man page for more information on command aliases, environments and 
# shell constructs.
#
# Last revised:  1/11/95  Blaise Barney
##############################################################################
umask 022
if ($?prompt) then
  set history=100
  set savehist=500
  set filec
  set prompt =  "`hostname -s`% "
  alias history 'history -r | more'
  alias h 'history -r | more'
  alias mv 'mv -i'
  alias cp 'cp -i'
  alias rm 'rm -i'
  # uncomment next line for MHPCC training room printer
  # alias print 'lpr -h -Pps3'
endif

#
# Set the search path for commands
#
set path = ( $path /usr/lpp/poe/bin /usr/lpp/LoadL/nfs/bin /usr/local/bin /source/local/bin $HOME/bin . /source/pd/gnu/bin/)

setenv MPATH .:/usr/man:/source/local/man:/u/khoros/man
setenv MANPATH .:/usr/man:/source/local/man:/u/khoros/man
setenv WORKSHOP /source/local/etc/httpd/htdocs/training/workshop
setenv WWW_HOME /source/local/etc/httpd/htdocs/mhpcc.html
setenv NNTPSERVER makapu.mhpcc.edu
setenv SP_NAME cws1.mhpcc.edu
setenv NWCHEM_TOP $HOME/nwchem
setenv NWCHEM_TARGET SP1
setenv TARGET SP1
setenv MPI_LOC /s/pd/msg_pass/mpich
alias make $HOME/bin/gmake
set ME = `whoami`
setenv MP_INFOLEVEL 1
setenv MP_EUILIB us 
setenv MP_RESD  YES
setenv SP_NAME cws1.mhpcc.edu 
setenv MP_EUIDEVICE  css0
setenv MP_CSS_INTERRUPT YES
setenv MP_STDINMODE 0
setenv MP_PULSE 0
