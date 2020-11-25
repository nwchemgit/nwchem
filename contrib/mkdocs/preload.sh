#!/bin/bash
  cd nwchemgit.github.io/
  git pull
  patch -p1  < ../preload.patch 
  git pull;git commit -m preload index.html ; git push
