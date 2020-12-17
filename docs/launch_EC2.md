# Launch and EC2 Instance in Julia or Python  

See our [Youtube](https://www.youtube.com/watch?v=0hGoK1SdBm4) video for a walkthrough or follow our checklist below!

## Checklist:
1. From the AMI Console choose `Ubuntu Server 18.04 LTS (HVM), SSD Volume Type 64 bit x86`
2. We recommend choosing `r5` and `c5` instances, with SSD-type storage, for memory and computation-heavy workflows.
3. Add disk space. Again job dependent based on I/O use (we run 200+ GB).
4. Create and download `.pem` security key and open terminal in the same directory as the key. 
5. To connect to your virtual machine right click and select connect from the menu. Copy the two commands which look like this into terminal. Note that you might have to change `root` to `ubuntu`, but terminal will prompt:
    -> `chmod 400 my_aws_key.pem`
    -> `ssh -i "my_aws_key.pem" ubuntu@ec2-XXX-XXX-XXX-XXX.us-west-2.compute.amazonaws.com`
    --> chmod command, then ssh
5. Run the following to clone [Tim Clements'](https://github.com/tclements) installer
    ```git clone https://github.com/tclements/SCEDC-AWS.git```
6. After clone is complete, the next step is to install Julia or python from the command line. You'll want to be in the same directory as the `.sh` file. Both the Julia folder and python are under `src` then `build_environment`. For julia `cd` into the `julia` folder then run `sh install_julia.sh`. For the most recent version of Julia (@v1.5) simply copy paste the following into the command line. 
```
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
```
### Python - refer to Tim's (github)[https://github.com/tclements/SCEDC-AWS/tree/master/src/build_environment/python]
6.1. Run `bash miniconda.sh` to install miniconda
6.2. Run `source ~/.bashrc` to update conda
6.3. Run `bash build_env.sh` to install seismology Python3 environment py3
7. Install needed packages using the package manager type `]` then run `add PACKAGE` to add packages
    7.1: To add the SCEDC package you'll need to run:
        `using Pkg; Pkg.add(PackageSpec(url="https://github.com/tclements/SCEDC.jl", rev="master"))`
    7.2: Remember packages don't update on instance! From the package manager use the `status` command to check current package versions

8. To access S3 cloud storage you'll need to add an IAM role:
    8.1: Right click on instance
    8.2: Click Instance Settings
    8.3: Attach/replace AIM role
    8.4: Create New IAM -> Create Role -> (Common) EC2
    8.5: Search for "AmazonS3FullAccess" -> click through tags
    8.6: Name IAM role (ex: "full_S3_access_from_EC2") and Create
    8.7: Add role

9. Create an IAM role to save the installation of Julia or python and packages to speed up launching your next instance. This can be found in the right-click drop down. 


