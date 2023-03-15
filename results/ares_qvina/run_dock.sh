#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J e192
## Liczba alokowanych węzłów
#SBATCH -N 1
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --cpus-per-task=32
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=1GB
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=10:00:00

#SBATCH -A plgporphconj-cpu

## Specyfikacja partycji
#SBATCH -p plgrid
#SBATCH --mail-user=anton.siomchen@student.uj.edu.pl
#SBATCH --mail-type=ALL


cd $SLURM_SUBMIT_DIR

source $1

