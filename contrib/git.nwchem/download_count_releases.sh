#!/bin/bash 
 curl -s https://api.github.com/repos/nwchemgit/nwchem/releases|egrep dow|sed 's/\"browser\_download\_url\":\ \"https:\/\/github.com\/nwchemgit\/nwchem\/releases\/download\///g' | sed 's/\"//g'
# curl -s https://api.github.com/repos/nwchemgit/nwchem/releases|egrep dow|sed 's/\"browser\_download\_url\":\ \"https:\/\/github.com\/nwchemgit\/nwchem\/releases\/download\//g'
