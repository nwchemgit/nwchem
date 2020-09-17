#!/bin/bash
sed \
   -e  's/<img alt=\"\$\$/\$\$/' \
   -e  's/<img alt=\"\\(/\\(/' \
   -e  's/\" src=\"https:\/\/raw.githubusercontent.com\/wiki\/nwchemgit\/nwchem\/svgs.*//' \
   -e  's/\" src=\"svgs.*//' \
$1 
