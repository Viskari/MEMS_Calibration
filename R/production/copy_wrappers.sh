#!/bin/bash

# This script creates copies of "wrapper_01.R" overwriting any existing copies.

# A common use case is running multiple calibrations simultaneously, e.g. with 
# different datasets or samplers: First, copies of "wrapper_01.R" are created
# with this shell script. Then, the copied wrappers are modified to use the 
# desired datasets/sampler. Finally, the wrappers are run simultaneously with 
# the shell script "run_wrappers.sh".

# Once multiple copies of "wrapper_01.R" already exist, all further development
# to the wrappers should be done to "wrapper_01.R", since the changes can be 
# easily updated to the copied wrappers with this script.

# Path to the folder with this script
path_main="/home/users/pusajann/projects/yasso_github/YASSO-calibration/R/production"

# Create copies of the first wrapper overwriting any existing copies
cp $path_main/wrapper_01.R $path_main/wrapper_02.R
cp $path_main/wrapper_01.R $path_main/wrapper_03.R
cp $path_main/wrapper_01.R $path_main/wrapper_04.R
