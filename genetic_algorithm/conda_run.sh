#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J GA_v1_test_50_50_jobs_8
## Liczba alokowanych węzłów
#SBATCH -N 1
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --cpus-per-task=16
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=1GB
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=7:00:00

#SBATCH -A plgporphconj-cpu

cd $SLURM_SUBMIT_DIR

module purge
module load miniconda3/4.9.2
conda activate /net/ascratch/people/plgasiomchen/docking_env

python -u main.py --population_size 50 --generations 50 --n_jobs 12


