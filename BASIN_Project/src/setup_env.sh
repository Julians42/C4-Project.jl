#!/bin/bash
sudo apt install unzip

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install



git clone https://github.com/Julians42/C4-Project.jl.git

cd C4-Project.jl
git checkout Develop

bash BASIN_Project/src/install_julia.sh