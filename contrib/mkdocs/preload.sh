#!/usr/bin/env bash
  rm -rf nwchemgit.github.io_temp
  git clone --depth 1 git@github.com:nwchemgit/nwchemgit.github.io.git nwchemgit.github.io_temp
cd nwchemgit.github.io_temp
  git pull
  patch -p1  < ../preload.patch 
  git pull;git commit -m preload index.html ; git push
  rm -rf nwchemgit.github.io_temp
