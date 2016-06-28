I#!/bin/bash
#$ -N AJK1 
#$ -V -cwd 
#$ -q short.q 
#$ -l mem_free=1G,h_vmem=1.2G
R CMD BATCH /home/eideakuc/model_flu/sero_model/main_model.R