#!/bin/bash
if [ -d "/shared" ]
then
    cd /shared
else
    cd
fi
wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.0-linux-x86_64.tar.gz
tar xvfa julia-1.5.0-linux-x86_64.tar.gz
rm julia-1.5.0-linux-x86_64.tar.gz
if [ -d "/shared" ]
then
    echo PATH=\$PATH:/shared/julia-1.5.0/bin/ >> ~/.bashrc
else
    echo PATH=\$PATH:~/julia-1.5.0/bin/ >> ~/.bashrc
fi
source ~/.bashrc
julia ~/SCEDC-AWS/src/build_environment/julia/add-packages.jl
cd