#!/bin/sh
#$ -S /bin/bash
#$ -q p100
#$ -cwd
#$ -l ngpus=1
#$ -e error
#$ -o outpt

source /home/cudasoft/bin/startcuda.sh 

echo $CUDA_VISIBLE_DEVICES
echo $NSLOTS 

source /home/cudasoft/bin/end_cuda.sh
