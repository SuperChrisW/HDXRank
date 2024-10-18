#!/bin/bash
#SBATCH --cpus-per-task=4

srun python /home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/BiLSTM+GAT_v2/pepGraph_generation/pepGraph_preprocessing.py