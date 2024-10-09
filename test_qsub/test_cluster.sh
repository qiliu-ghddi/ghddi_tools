#!/bin/sh
#$ -S /bin/bash
#$ -q v100
#$ -cwd
#$ -l ngpus=4
#$ -e error
#$ -o outpt

source /home/cudasoft/bin/startcuda.sh 

echo ">>> Current path:"
pwd
echo ">>> All envs:"
conda info --envs
conda activate cryodrgn

echo CUDA VISIABLE DEVICES: $CUDA_VISIBLE_DEVICES
echo NSLOTS $NSLOTS 


echo ">>> Current env:"
conda info --envs
conda deactivate cryodrgn

source /home/cudasoft/bin/end_cuda.sh

