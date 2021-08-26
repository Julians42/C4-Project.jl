#!/bin/bash

# install unzip - needed for aws cli installation 
sudo apt install unzip

# install aws cli - needed for aws configure (run after script completes to add credentials)
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install


# clone project and set up julia
git clone https://github.com/Julians42/C4-Project.jl.git

cd C4-Project.jl
git checkout Develop

bash BASIN_Project/src/install_julia.sh

# add packages
julia add-packages.jl