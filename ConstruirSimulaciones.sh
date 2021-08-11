rm -r ParallelExampleOut
mkdir ParallelExampleOut

ml icc/2019.2.187-GCC-8.2.0-2.31.1-nlhpc  impi/2019.2.187-nlhpc
ml ifort/2019.2.187-GCC-8.2.0-2.31.1  impi/2019.2.186
ml DOLFIN/2018.1.0.post1-Python-3.7.2
ml MSHR/2018.1.0-Python-3.7.2

python CrearTareas.py
sbatch EjecutarTareas.sh

