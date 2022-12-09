#!/bin/bash
git config user.name "nwchemgit"
git config user.email "nwchemgit@gmail.com"
git config advice.updateSparsePath false
#git add -A
git add .
git commit -m website\ changes\ on\ `date +%a_%b_%d_%H:%M:%S_%Y` 
#git diff --stat --cached origin/master

