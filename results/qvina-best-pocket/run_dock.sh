#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J best_128
## Liczba alokowanych węzłów
#SBATCH -N 1
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --cpus-per-task=16
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=300MB
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=10:00
#SBATCH -A plgporphconj-cpu
## Specyfikacja partycji
#SBATCH -p plgrid-testing

cd $SLURM_SUBMIT_DIR
source $1 $2 $3
