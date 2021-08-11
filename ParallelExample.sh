#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -J MultiFenics
#SBATCH -p slims
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem=8000
#SBATCH --array=1-NumeroArchivos
#SBATCH -o ParallelExampleOut/MultiFenics_%A_%a.out
#SBATCH -e ParallelExampleOut/MultiFenics_%A_%a.err

cd /home/dvasquez

ml icc/2019.2.187-GCC-8.2.0-2.31.1-nlhpc  impi/2019.2.187-nlhpc
ml ifort/2019.2.187-GCC-8.2.0-2.31.1  impi/2019.2.186
ml DOLFIN/2018.1.0.post1-Python-3.7.2
ml MSHR/2018.1.0-Python-3.7.2

# ----------------Comandos--------------------------
file=$(ls Tareas/Tarea_*[0-9].py | sed -n ${SLURM_ARRAY_TASK_ID}p)
srun python $file


