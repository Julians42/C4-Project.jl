# The BASIN Project

This folder contains the workflows I developed for computing BASIN cross-correlations on AWS, subsequent plotting, and side investigations. For the plotting routines, we primarily used jupyter notebooks for their code and output integration. The correlating scripts will work out of the box but much of the notebooks, which we ran locally, will need filepath edits to run.

## Setting up AWS instance 
We've made it easy to set up an amazon instance for our computing workflow! Simply run [`setup_env.sh`](src/setup_env.sh) (copy/paste) which installs the AWS CLI, clones this github repo, and installs needed packages. You can then add your credentials (we recommend IAM creds) using `aws configure` (`region: us-west-2` and `output: json`). 

### Scripts and the BASIN Module
[`basin_correlations.jl`](scripts/basin_correlations.jl) is the finalized verson of the AWS script for computing basin cross correlations. Finally run the script, e.g. `julia basin_correlations.jl 2019` to run correlations for 2019. 

[`BASIN.jl`](src/BASIN.jl) is the module in which we store much of the functionality for the basin script. It keeps the actual script nice and clean and allows us to differentiate between the larger routines we use, such as our preprocessing routine, and the smaller functions which deal with lower level details such as file pathing or a single correlation.


## Notebooks
I've included 7 notebooks which I used for processing the stacked h5 files and several for looking at daily stacks (I have a few daily files locally so if needed please ping me - they're not anywhere else!). A couple of notes on these:
1. [`make_csv.ipynb`](notebooks/make_csv.ipynb): The notebook performs the intermediary step of extracting amplitude and SNR information from the stacked XCORR files. These are saved to CSVs and whch are then read by the plotting script to make heatplots. This uses some of the functions from the module on C4, should be easy to either copy those files and change filepathing or just copy the few needed functions from github into the notebook.
2. [`max_amp_heatplots.ipynb`](notebooks/max_amp_heatplots.ipynb): Makes maximum amplitude heatplots and SNR plots from the csvs computed in `make_csv.ipnb`. There is additional functionality not currently used which is able to adjust the node lines by the intersections - this is the `ratio` and `get_ratios` functions. There are a couple of other investigations that didn't prove fruitful later on in the notebook and I've left those codes as is. 
3. [`debug_stack.ipynb`](notebooks/debug_stack.ipynb): Scratch work and various investigatons.
4. [`Stack_Investigation.ipynb`](notebooks/Stack_Investigation.ipynb): More stacking stuff
5. [`Stackplots.ipynb`](notebooks/Stackplots.ipynb): Contains several different versions of codes which make single line moveout plots from the BASIN XCORR h5 files. In-line notation and markdown cells might be helpful for deciding what is useful.
6. [`RMS.ipynb`](notebooks/RMS.ipynb): Investigation of the RMS levels across nodes. 
7. [`Ampltudes_through_time.ipynb`](notebooks/Amplitudes_through_time.ipynb): Here we investigated the effects of seasonality on maximum amplitude by looking at daily stacks. As the daily correlations aren't well converged this whole inquiry was somewhat controversial.

## Thank you for Continuing the Project!
Excited to see where you are able to take the project! Please feel free to reach out with questions or open an issue.
## :cloud: :earth_americas: :mount_fuji: :ski: