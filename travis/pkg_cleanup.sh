#!/usr/bin/env bash
df -h
if [[ $(uname -s) == "Linux" ]]; then
ubuntu_ver=$(cat /etc/os-release | grep VERSION_ID |cut -d \" -f 2)
dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n
sudo apt-get purge -y azure-cli || true
sudo apt-get purge -y google-cloud-cli microsoft-edge-stable dotnet-sdk-7.0 dotnet-sdk-6.0 google-chrome-stable firefox
sudo apt-get purge -y temurin-17-jdk temurin-11-jdk temurin-8-jdk
if [[ $ubuntu_ver == "20.04" ]]; then
sudo apt-get purge -y llvm-12-dev llvm-11-dev llvm-10-dev
sudo apt-get purge -y hhvm
sudo apt-get purge -y libgl1-mesa-dri
fi
if [[ $ubuntu_ver == "22.04" ]]; then
sudo apt-get purge -y llvm-13-dev llvm-14-dev llvm-15-dev
fi
sudo apt-get -y clean
sudo apt-get autoremove -y 
dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n
df -h
fi
