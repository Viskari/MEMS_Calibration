#!/bin/bash

# This script runs the wrappers in background allowing closure 
# of ssh connection to the server. The logs of the runs are saved to the folder 
# specified in this file and the results of the calibrations to the folders
# specified in the wrapper files.

# A common use case is running multiple calibrations simultaneously, e.g. with
# different datasets or samplers: First, copies of "wrapper_01.R" are created
# with "copy_wrappers.sh". Then, the copied wrappers are modified to use the 
# desired datasets/sampler. Finally, the wrappers are run simultaneously with 
# this shell script.

# Path to the logs folder
path_main="/home/users/pusajann/projects/yasso_github/YASSO-calibration/log"

# Create the logs folder if it does not exist
mkdir -p $path_main

# Ways to use the script:

# Run multiple wrappers in parallel
# nohup Rscript wrapper_01.R > $path_main/wrp1_log.txt 2>&1 </dev/null &
# nohup Rscript wrapper_02.R > $path_main/wrp2_log.txt 2>&1 </dev/null &
# nohup Rscript wrapper_03.R > $path_main/wrp3_log.txt 2>&1 </dev/null &
# nohup Rscript wrapper_04.R > $path_main/wrp4_log.txt 2>&1 </dev/null &

# Run multiple wrappers in parallel and sequentially. This mixed approach is used
# if server CPU/memory is not sufficient to run everything in parallel.
# (nohup Rscript wrapper_01.R > $path_main/wrp1_log.txt 2>&1 </dev/null; nohup Rscript wrapper_03.R > $path_main/wrp3_log.txt 2>&1 </dev/null) &
# (nohup Rscript wrapper_02.R > $path_main/wrp2_log.txt 2>&1 </dev/null; nohup Rscript wrapper_04.R > $path_main/wrp4_log.txt 2>&1 </dev/null) &
