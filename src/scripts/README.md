# Documentation for Running the C4 job on Batch

In the following, I document some of the technical aspects of creating the C4 workflow, where to watch out for bugs (some bugs are expensive!), and some of the steps to reproduce. 
## Workflow
There are essentially 4 components. I'll outline each below
1. [The Script](BC_ec2.jl) :grin: The exectable which contains the codes for processing, correlating, and stacking the C4 project. This is a good old `julia root/scripts/BC_ec2.jl` into the command line. A few considerations when running with batch is that we need some way of telling each job to run something different. Eg, if we run identical scripts we'll do the same job *n* times! Luckily AWS batch by default generates a unique identifier for batch array jobs in instance environment - which we leverage by getting an index via `INDEX = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]`. Other than that, fine tuning package dependencies and literally 700 lines of code was the other hard part! At one point in the project I developed another package, called [SeisCore.jl](https://github.com/Julians42/SeisCore.jl), which reduced the length of the script to about 50 lines - see [`bc_script.jl`](bc_script.jl) with all the functions wrapped nicely in that package. However, it made debugging too slow to update the package every time something needed to change so I switched back to just one giant file. 
2. [The docker Image](../../Dockefile) Docker images are containerized applications above your operating system in which you build the environment which your script needs. Basically, its the script's mom and dad. In the case of this project, we set up our docker image (available at [`vtskier/seis_core`](https://hub.docker.com/repository/docker/vtskier/seis_core)) to have the latest version of `Julia` (1.5.3 at the time), hold all the packages (deps we used can be found in the [`Project.toml`](../../Project.toml) file), as well as the script itself! Annoyingly, docker images are by default limited to 10 GB of disk space, which is not good enough for the TB throughput we're looking for - more on that below. 

3. The AMI. If we continue with the reference, the AMI is the grandparent, as its the operating system in which the docker image is constructed. Amazon also has has a special AMI that we have to make only a few modifications to before we can run our job. We used Amazon's ECS-Optimized Amazon Linux 2 AMI ([ECS AWS Linux 2](https://aws.amazon.com/marketplace/pp/B07KMLLN73?ref=cns_srchrow)) which is optimized for running container jobs! In order to increase the amount of space available to our docker image, as discussed above, we need to make our own custom AMI, add a volume for docker to write files to, and finally mount the volume during the batch job creation to give our docker image access. We called the volume `/docker_scratch` and mounted to `/scratch` during job creation and so we write and read all our files from there in the script. 
    1. We make a custom AMI by launching an instance with this AMI. When launching in EC2 at `Step 4: Add Storage` we click `Add New Volume` and select `/dev/sdb`. We then arbitrarily adjust the size to whatever we want (eg 750 GB, 5 TB...) and check "Delete on Termination" for both volumes so that whenever we launch a job with those volumes, they are automatically deleted when the new job finishes. 
    2. After connecting to the instance, we update the image and set of the scratch volume via
    ```
    sudo yum -y update
    sudo mkfs -t ext4 /dev/xvdb
    sudo mkdir /docker_scratch
    sudo echo -e '/dev/xvdb\t/docker_scratch\text4\tdefaults\t0\t0' | sudo tee -a /etc/fstab

    sudo mount –a
    sudo systemctl stop ecs
    sudo rm -rf /var/lib/ecs/data/*
    ```
    3. We then create the AMI from the console by right clicking on the running instance, Image and templates, then Create image.

4. Batch - batch is essentially the orchestra conductor for all these jobs - once jobs are launched, batch manages the starting, stopping, and resource cleanup of your jobs. We use it because it scales incredibly well, and unlike single instances, you don't have to even have your computer turned on for jobs to run. For our purposes, batch is a little overdeveloped because it contains three levels (of confusion): `Compute environments`, `Job queues`, and `Job definitions` which are combined when you launch a job. 
    1. When creating a compute environment we used Spot instances as our provisioning model as they are generally 70% cheaper with a small chance of interruption (it's worth the small risk). We are then able to select maximum and desired CPUs based on the job we want to complete. Note that these limits MUST be larger than the total CPUs you want for all your jobs you want to run simultaneously, otherwise your jobs will end up getting stuck in the `runnable` state (fairly common issue). Under additional settings, we also specify the AMI we just created, by checking enable, and then copying the AMI ID over from the EC2 console that we just created. Everything else should be default or straightforward. 
    2. The job queue is easy. Just select the compute environment.
    3. For the job definition we need to specify the `EC2` platform and number of job attempts (generally, 1 or 2). The execution timeout is nice to specify because it tells batch after how long it should kill a job and saves things from bugging and never terminating. Under container properties, we specify the image we created earlier from docker, in this case `vtskier/seis_core:latest` and enter a bash command we want executed (`julia root/scripts/BC_single_instance`). Here we also select the number of CPUs and how much memory we want reserved (32CPUs/200,000MB is what we used). Finally, we add a mount point and volume so that our docker image has access to the larger ECS storage volume. For the mount point we set `Container Path = /scratch` and `Source volume = docker_scratch`. For the volume, we set `Name = docker_scrach` and `Source path = /docker_scratch`. 
    4. Finally, when we launch the job we must also specify the array size if we want to run multiple jobs - for example if we want to run each month separately for a year, our array size is 12, and we map the ENV variable in the script accordingly. 

### Thank You 
This project was advised by [Marine Denolle](https://github.com/mdenolle) and mentored by [Tim Clements](https://github.com/tclements) who taught me how to Amazon and wrote [SeisNoise.jl](https://github.com/tclements/SeisNoise.jl) one of the packages my code is based on! If you're looking to reproduce these steps, be patient, expect a lot of bugs, and be ready to spend big! Feel free to get in touch!



