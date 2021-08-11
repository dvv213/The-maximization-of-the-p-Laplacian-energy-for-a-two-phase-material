#!/bin/bash
#SBATCH -J tarea0
#SBATCH -p general
#SBATCH -n 1
#SBATCH --output=tarea0_%j.out
#SBATCH --error=tarea0_%j.err
#SBATCH --mail-user=donato.rolxxx@gmail.com
#SBATCH --mail-type=ALL
cd /home/dvasquez
python ejemplo0.py
